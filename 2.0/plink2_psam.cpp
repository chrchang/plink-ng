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
#include "plink2_psam.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct psam_info_ll_struct {
  // vardata[] starts with 8-byte phenotypes (we don't want to parse the same
  // same numeric string twice), followed by NON-null-terminated sample_id, and
  // then non-terminated paternal and maternal IDs.
  struct psam_info_ll_struct* next;
  uint32_t sample_id_slen;
  uint32_t sid_slen;
  uint32_t paternal_id_slen;
  uint32_t maternal_id_slen;
  uint32_t sex_code; // 0 = unknown, 1 = male, 2 = female
  unsigned char vardata[];
} psam_info_ll_t;

pglerr_t load_psam(const char* psamname, const range_list_t* pheno_range_list_ptr, fam_col_t fam_cols, uint32_t pheno_ct_max, int32_t missing_pheno, uint32_t affection_01, uintptr_t* max_sample_id_blen_ptr, uintptr_t* max_sid_blen_ptr, uintptr_t* max_paternal_id_blen_ptr, uintptr_t* max_maternal_id_blen_ptr, uintptr_t** sample_include_ptr, char** sample_ids_ptr, char** sids_ptr, char** paternal_ids_ptr, char** maternal_ids_ptr, uintptr_t** founder_info_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, pheno_col_t** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* raw_sample_ct_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
  // outparameter pointers assumed to be initialized to nullptr
  //
  // pheno_ct_max should default to something like 0x7fffffff, not 0xffffffffU
  //
  // max_{sample,sid,paternal,maternal}_id_blen are in/out, to support data
  // management operations which change these values
  //
  // permanent allocations are at stack end, not base, to work better with
  // variant_id_htable_find()

  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  
  gzFile gz_infile = nullptr;
  pheno_col_t* pheno_cols = nullptr;
  uintptr_t line_idx = 0;
  uint32_t pheno_ct = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    reterr = gzopen_read_checked(psamname, &gz_infile);
    if (reterr) {
      goto load_psam_ret_1;
    }
    const uintptr_t initial_bigstack_size = bigstack_left();
    uintptr_t loadbuf_size = initial_bigstack_size / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto load_psam_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    // allocated at bottom now, so short string comparsions against end cannot
    // fail
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_first_token;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_psam_ret_READ_FAIL;
	}
	loadbuf_first_token = loadbuf;
	loadbuf_first_token[0] = '\0';
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto load_psam_ret_LONG_LINE;
	}
	goto load_psam_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token) || ((loadbuf_first_token[0] == '#') && strcmp_se(&(loadbuf_first_token[1]), "FID", 3) && strcmp_se(&(loadbuf_first_token[1]), "IID", 3)));
    const uint32_t pheno_name_subset = pheno_range_list_ptr && pheno_range_list_ptr->names;
    uint32_t* col_skips = nullptr;
    uint32_t* col_types = nullptr;
    uint32_t psam_cols_mask = 0;
    ll_str_t* pheno_names_reverse_ll = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    uint32_t relevant_postfid_col_ct = 0;
    g_bigstack_end -= kMaxIdSlen;
    unsigned char* tmp_bigstack_end = g_bigstack_end;
    unsigned char* bigstack_mark2;
    if (loadbuf_first_token[0] == '#') {
      // parse header
      // [-1] = #FID (if present, must be first column)
      // [0] = IID (could also be first column)
      // [1] = SID
      // [2] = PAT
      // [3] = MAT
      // [4] = SEX
      // [5+] = phenotypes
      relevant_postfid_col_ct = count_tokens(loadbuf_first_token);
      if (relevant_postfid_col_ct > pheno_ct_max + 5) {
	relevant_postfid_col_ct = pheno_ct_max + 5;
      }
      if (bigstack_alloc_ui(relevant_postfid_col_ct, &col_skips) ||
	  bigstack_alloc_ui(relevant_postfid_col_ct, &col_types)) {
	goto load_psam_ret_NOMEM;
      }
      bigstack_mark2 = g_bigstack_base;
      uint32_t rpf_col_idx = 0;
      if (loadbuf_first_token[1] == 'I') {
	col_skips[0] = 0;
	col_types[0] = 0;
	++rpf_col_idx;
	psam_cols_mask = 1;
      }
      uint32_t col_idx = 0;
      uint32_t in_interval = 0;
      char* cmdline_pheno_sorted_ids = nullptr;
      uint32_t* cmdline_pheno_id_map = nullptr;
      uintptr_t max_cmdline_pheno_id_blen = 0;
      uintptr_t cmdline_pheno_name_ct = 0;
      if (pheno_name_subset) {
	max_cmdline_pheno_id_blen = pheno_range_list_ptr->name_max_blen;
	cmdline_pheno_name_ct = pheno_range_list_ptr->name_ct;
	uintptr_t* dummy_bitarr;
	// don't bother freeing these before load_psam() is done
	if (bigstack_alloc_c(cmdline_pheno_name_ct * max_cmdline_pheno_id_blen, &cmdline_pheno_sorted_ids) ||
	    bigstack_alloc_ui(cmdline_pheno_name_ct, &cmdline_pheno_id_map) ||
	    bigstack_alloc_ul(BITCT_TO_WORDCT(cmdline_pheno_name_ct), &dummy_bitarr)) {
	  goto load_psam_ret_NOMEM;
	}
	fill_all_bits(cmdline_pheno_name_ct, dummy_bitarr);
	reterr = copy_sort_strbox_subset_noalloc(dummy_bitarr, pheno_range_list_ptr->names, cmdline_pheno_name_ct, max_cmdline_pheno_id_blen, 0, 0, 0, cmdline_pheno_sorted_ids, cmdline_pheno_id_map);
	if (reterr) {
	  goto load_psam_ret_1;
	}
	bigstack_reset(dummy_bitarr);
      }
      char* token_end = &(loadbuf_first_token[4]);
      unsigned char* ll_alloc_base = g_bigstack_base;
      while (1) {
        char* loadbuf_iter = skip_initial_spaces(token_end);
	if (is_eoln_kns(*loadbuf_iter)) {
	  break;
	}
	++col_idx;
	token_end = token_endnn(loadbuf_iter);
	const uint32_t token_slen = (uintptr_t)(token_end - loadbuf_iter);
	if (token_slen == 3) {
	  uint32_t cur_col_type = 0xffffffffU;
	  if (!memcmp(loadbuf_iter, "IID", 3)) {
	    cur_col_type = 0;
	  } else if (!memcmp(loadbuf_iter, "SID", 3)) {
	    cur_col_type = 1;
	  } else if (!memcmp(loadbuf_iter, "PAT", 3)) {
	    cur_col_type = 2;
	  } else if (!memcmp(loadbuf_iter, "MAT", 3)) {
	    cur_col_type = 3;
	  } else if (!memcmp(loadbuf_iter, "SEX", 3)) {
	    cur_col_type = 4;
	  } else if (!memcmp(loadbuf_iter, "FID", 3)) {
	    sprintf(g_logbuf, "Error: 'FID' column header on line %" PRIuPTR " of %s is not at the beginning.\n", line_idx, psamname);
	    goto load_psam_ret_MALFORMED_INPUT_WW;
	  }
	  if (cur_col_type != 0xffffffffU) {
	    const uint32_t cur_col_type_shifted = 1 << cur_col_type;
	    if (psam_cols_mask & cur_col_type_shifted) {
	      *token_end = '\0';
	      sprintf(g_logbuf, "Error: Duplicate column header '%s' on line %" PRIuPTR " of %s.\n", loadbuf_iter, line_idx, psamname);
	      goto load_psam_ret_MALFORMED_INPUT_WW;
	    }
	    psam_cols_mask |= cur_col_type_shifted;
	    col_skips[rpf_col_idx] = col_idx;
	    col_types[rpf_col_idx++] = cur_col_type;
	    continue;
	  }
	}
	if (pheno_ct < pheno_ct_max) {
	  if (pheno_name_subset) {
	    uint32_t cmdline_pos;
	    if (!sorted_idbox_find(loadbuf_iter, cmdline_pheno_sorted_ids, cmdline_pheno_id_map, token_slen, max_cmdline_pheno_id_blen, cmdline_pheno_name_ct, &cmdline_pos)) {
	      // similar to string_range_list_to_bitarr()
	      if (pheno_range_list_ptr->starts_range[cmdline_pos]) {
		if (in_interval) {
		  logerrprint("Error: Overlapping --pheno-name ranges.\n");
		  goto load_psam_ret_INCONSISTENT_INPUT;
		}
		in_interval = 1;
	      } else if (cmdline_pos && pheno_range_list_ptr->starts_range[cmdline_pos - 1]) {
		if (!in_interval) {
		  sprintf(g_logbuf, "Error: --pheno-name range is inconsistent with %s.\n", psamname);
		  goto load_psam_ret_INCONSISTENT_INPUT_WW;
		}
		in_interval = 0;
	      }
	    } else if (!in_interval) {
	      continue;
	    }
	  }
	  const uint32_t tok_blen = token_slen + 1;
	  ll_str_t* ll_str_new = (ll_str_t*)ll_alloc_base;
	  // just word-aligned, not cacheline-aligned
	  ll_alloc_base += round_up_pow2(tok_blen + sizeof(ll_str_t), kBytesPerWord);
	  if (ll_alloc_base > tmp_bigstack_end) {
	    goto load_psam_ret_NOMEM;
	  }
	  ll_str_new->next = pheno_names_reverse_ll;
	  memcpyx(ll_str_new->ss, loadbuf_iter, token_slen, '\0');
	  if (tok_blen > max_pheno_name_blen) {
	    max_pheno_name_blen = tok_blen;
	  }
	  pheno_names_reverse_ll = ll_str_new;
	  col_skips[rpf_col_idx] = col_idx;
	  col_types[rpf_col_idx++] = pheno_ct + 5;
	  ++pheno_ct;
	}
      }
      if (max_pheno_name_blen > kMaxIdBlen) {
	logerrprint("Error: Phenotype/covariate names are limited to " MAX_ID_SLEN_STR " characters.\n");
	goto load_psam_ret_MALFORMED_INPUT;
      }
      g_bigstack_base = (unsigned char*)round_up_pow2((uintptr_t)ll_alloc_base, kCacheline);
      if (!(psam_cols_mask & 1)) {
	sprintf(g_logbuf, "Error: No IID column on line %" PRIuPTR " of %s.\n", line_idx, psamname);
	goto load_psam_ret_MALFORMED_INPUT_WW;
      }
      if (in_interval) {
	sprintf(g_logbuf, "Error: --pheno-name range is inconsistent with %s.\n", psamname);
	goto load_psam_ret_INCONSISTENT_INPUT_WW;
      }
      relevant_postfid_col_ct = rpf_col_idx;
      for (rpf_col_idx = relevant_postfid_col_ct - 1; rpf_col_idx; --rpf_col_idx) {
	col_skips[rpf_col_idx] -= col_skips[rpf_col_idx - 1];
      }
      loadbuf_first_token[0] = '\0'; // forces line to be skipped by main loop
    } else if (loadbuf_first_token[0]) {
      if (pheno_name_subset) {
	logerrprint("Error: --pheno-name requires a --pheno or .psam file with a header.\n");
	goto load_psam_ret_INCONSISTENT_INPUT;
      }
      
      pheno_ct = (fam_cols & kfFamCol6) && pheno_ct_max;
      relevant_postfid_col_ct = ((fam_cols / kfFamCol1) & 1) + ((fam_cols / (kfFamCol34 / 2)) & 2) + ((fam_cols / kfFamCol5) & 1) + pheno_ct;
      // these small allocations can't fail, since kMaxMediumLine <
      // loadbuf_size <= 1/3 of remaining space
      col_skips = (uint32_t*)bigstack_alloc_raw_rd(relevant_postfid_col_ct * sizeof(int32_t));
      col_types = (uint32_t*)bigstack_alloc_raw_rd(relevant_postfid_col_ct * sizeof(int32_t));
      bigstack_mark2 = g_bigstack_base;
      col_skips[0] = fam_cols & 1; // assumes kfFamCol1 == 1
      col_types[0] = 0;
      // psam_cols_mask = 1; // may need this later
      uint32_t rpf_col_idx = 1;
      if (fam_cols & kfFamCol34) {
	col_skips[rpf_col_idx] = 1;
	col_types[rpf_col_idx++] = 2;
	col_skips[rpf_col_idx] = 1;
	col_types[rpf_col_idx++] = 3;
	psam_cols_mask |= 12;
      }
      if (fam_cols & kfFamCol5) {
	col_skips[rpf_col_idx] = 1;
	col_types[rpf_col_idx++] = 4;
	psam_cols_mask |= 0x10;
      }
      if (pheno_ct) {
	col_skips[rpf_col_idx] = 1;
	// col_types[rpf_col_idx++] = 6;
	col_types[rpf_col_idx] = 5;
	ll_str_t* ll_str_new = (ll_str_t*)bigstack_alloc_raw_rd(7 + sizeof(ll_str_t));
	ll_str_new->next = pheno_names_reverse_ll;
	strcpy(ll_str_new->ss, "PHENO1");
	max_pheno_name_blen = 7;
	pheno_names_reverse_ll = ll_str_new;
      }
    }
    if (pheno_ct) {
      char* pheno_names = (char*)malloc(pheno_ct * max_pheno_name_blen);
      if (!pheno_names) {
	goto load_psam_ret_NOMEM;
      }
      *pheno_names_ptr = pheno_names;
      for (uint32_t pheno_idx = pheno_ct; pheno_idx;) {
	--pheno_idx;
	strcpy(&(pheno_names[pheno_idx * max_pheno_name_blen]), pheno_names_reverse_ll->ss);
	pheno_names_reverse_ll = pheno_names_reverse_ll->next;
      }
      if (pheno_ct > 1) {
	if (pheno_ct > kMaxPhenoCt) {
	  // yeah, yeah, this will never come up
	  logerrprint("Error: " PROG_NAME_STR " does not support more than " MAX_PHENO_CT_STR " phenotypes.\n");
	  goto load_psam_ret_MALFORMED_INPUT;
	}
	// verify there are no duplicates
	uint32_t tmp_htable_size;
	uint32_t* htable_tmp;
	if (htable_good_size_alloc(pheno_ct, bigstack_left(), &htable_tmp, &tmp_htable_size)) {
	  goto load_psam_ret_NOMEM;
	}
	const uint32_t duplicate_idx = populate_strbox_htable(pheno_names, pheno_ct, max_pheno_name_blen, tmp_htable_size, htable_tmp);
	if (duplicate_idx) {
	  const char* duplicate_pheno_name = &(pheno_names[duplicate_idx * max_pheno_name_blen]);
	  sprintf(g_logbuf, "Error: Duplicate phenotype/covariate name '%s' on line %" PRIuPTR " of %s.\n", duplicate_pheno_name, line_idx, psamname);
	  goto load_psam_ret_MALFORMED_INPUT_WW;
	}
      }
      // free pheno_names_reverse_ll
      bigstack_reset(bigstack_mark2);
    }
    
    // make sure to handle sample_ct == 0 case properly
    psam_info_ll_t* psam_info_reverse_ll = nullptr;
    const uint32_t sids_present = (psam_cols_mask / 2) & 1;
    const uint32_t paternal_ids_present = psam_cols_mask & 4;
    const uint32_t maternal_ids_present = psam_cols_mask & 8;
    const uint32_t sex_present = psam_cols_mask & 0x10;
    const uint32_t col_type_end = 5 + pheno_ct;
    const uint32_t pheno_ctl = BITCT_TO_WORDCT(pheno_ct);
    const double missing_phenod = (double)missing_pheno;
    const double pheno_ctrld = (double)((int32_t)(1 - affection_01));
    const double pheno_cased = pheno_ctrld + 1.0;
    uintptr_t max_sample_id_blen = *max_sample_id_blen_ptr;
    uintptr_t max_sid_blen = *max_sid_blen_ptr;
    uintptr_t max_paternal_id_blen = *max_paternal_id_blen_ptr;
    uintptr_t max_maternal_id_blen = *max_maternal_id_blen_ptr;
    uint32_t raw_sample_ct = 0;
    uint32_t categorical_pheno_ct = 0;
    
    char** token_ptrs;
    uint32_t* token_slens;
    uintptr_t* categorical_phenos;
    uintptr_t* quantitative_phenos;
    if (bigstack_alloc_cp(col_type_end, &token_ptrs) ||
        bigstack_alloc_ui(col_type_end, &token_slens) ||
	bigstack_calloc_ul(pheno_ctl, &categorical_phenos) ||
	bigstack_calloc_ul(pheno_ctl, &quantitative_phenos)) {
      goto load_psam_ret_NOMEM;
    }
    char* missing_catname = g_missing_catname;
    const uint32_t missing_catname_blen = strlen(missing_catname) + 1;
    const uint32_t missing_catname_hval = hashceil(missing_catname, missing_catname_blen - 1, kCatHtableSize);
    unsigned char* tmp_bigstack_base = g_bigstack_base;
    catname_ll2_t** catname_htable = nullptr;
    catname_ll2_t** pheno_catname_last = nullptr;
    uintptr_t* total_catname_blens = nullptr;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	if (raw_sample_ct == 0x7ffffffe) {
	  logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
	  goto load_psam_ret_MALFORMED_INPUT;
	}
	char* loadbuf_iter = loadbuf_first_token;
	for (uint32_t rpf_col_idx = 0; rpf_col_idx < relevant_postfid_col_ct; ++rpf_col_idx) {
	  const uint32_t cur_col_type = col_types[rpf_col_idx];
	  loadbuf_iter = next_token_multz(loadbuf_iter, col_skips[rpf_col_idx]);
	  if (!loadbuf_iter) {
	    goto load_psam_ret_MISSING_TOKENS;
	  }
	  token_ptrs[cur_col_type] = loadbuf_iter;
	  char* token_end = token_endnn(loadbuf_iter);
	  token_slens[cur_col_type] = (uintptr_t)(token_end - loadbuf_iter);
	  loadbuf_iter = token_end;
	}
	const uint32_t fid_slen = (uintptr_t)(token_endnn(loadbuf_first_token) - loadbuf_first_token);
	const uint32_t iid_slen = token_slens[0];
	const uint32_t sid_slen = sids_present? token_slens[1] : 0;
	const uint32_t paternal_id_slen = paternal_ids_present? token_slens[2] : 1;
	const uint32_t maternal_id_slen = maternal_ids_present? token_slens[3] : 1;
	// phenotypes
	if (!raw_sample_ct) {
	  for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	    if (is_categorical_phenostr_nocsv(token_ptrs[pheno_idx + 5])) {
	      SET_BIT(pheno_idx, categorical_phenos);
	    }
	  }
	  categorical_pheno_ct = popcount_longs(categorical_phenos, pheno_ctl);
	  if (categorical_pheno_ct) {
	    // initialize hash table
	    /*
	    if (categorical_pheno_ct > kCatHtableSize) {
	      // use a larger hash table if/when we ever care for this case
	      logerrprint("Error: " PROG_NAME_STR " does not support more than 2^19 - 1 categorical phenotypes.\n");
	      goto load_psam_ret_MALFORMED_INPUT;
	    }
	    */
	    const uint32_t cat_ul_byte_ct = categorical_pheno_ct * sizeof(intptr_t);
	    const uint32_t htable_byte_ct = kCatHtableSize * sizeof(uintptr_t);
	    const uintptr_t entry_byte_ct = round_up_pow2(offsetof(catname_ll2_t, ss) + missing_catname_blen, sizeof(intptr_t));	    
	    if ((uintptr_t)(tmp_bigstack_end - tmp_bigstack_base) < htable_byte_ct + categorical_pheno_ct * entry_byte_ct + 2 * cat_ul_byte_ct) {
	      goto load_psam_ret_NOMEM;
	    }
	    pheno_catname_last = (catname_ll2_t**)tmp_bigstack_base;
	    tmp_bigstack_base += cat_ul_byte_ct;	    
	    total_catname_blens = (uintptr_t*)tmp_bigstack_base;
	    tmp_bigstack_base += cat_ul_byte_ct;
	    fill_ulong_zero(categorical_pheno_ct, total_catname_blens);
	    catname_htable = (catname_ll2_t**)tmp_bigstack_base;
	    tmp_bigstack_base += htable_byte_ct;
	    for (uint32_t uii = 0; uii < kCatHtableSize; ++uii) {
	      catname_htable[uii] = nullptr;
	    }
	    uint32_t cur_hval = missing_catname_hval;
	    for (uint32_t cat_pheno_idx = 0; cat_pheno_idx < categorical_pheno_ct; ++cat_pheno_idx) {
	      catname_ll2_t* new_entry = (catname_ll2_t*)tmp_bigstack_base;
	      tmp_bigstack_base += entry_byte_ct;
	      pheno_catname_last[cat_pheno_idx] = new_entry;
	      new_entry->cat_idx = 0;
	      new_entry->htable_next = nullptr;
	      new_entry->pheno_next = nullptr;
	      memcpy(new_entry->ss, missing_catname, missing_catname_blen);
	      catname_htable[cur_hval++] = new_entry;
	      if (cur_hval == kCatHtableSize) {
		cur_hval = 0;
	      }
	    }
	  }
	}
	// 1 extra byte for tab between FID and IID; this gets absorbed into
	// the "+ sizeof(intptr_t)" at the end, since that would normally be
	// "+ (sizeof(intptr_t) - 1)"
	// bugfix: pheno_ct * sizeof(intptr_t) -> pheno_ct * 8
	const uint32_t alloc_byte_ct = sizeof(psam_info_ll_t) + sizeof(intptr_t) + round_down_pow2(fid_slen + iid_slen + sid_slen + paternal_id_slen + maternal_id_slen + pheno_ct * 8, sizeof(intptr_t));
	psam_info_ll_t* new_psam_info = (psam_info_ll_t*)tmp_bigstack_base;
	tmp_bigstack_base += alloc_byte_ct;
	if (tmp_bigstack_base > tmp_bigstack_end) {
	  goto load_psam_ret_NOMEM;
	}
	new_psam_info->next = psam_info_reverse_ll;
	char* sample_id_storage = (char*)(&(new_psam_info->vardata[pheno_ct * 8]));
	char* ss_iter = memcpyax(sample_id_storage, loadbuf_first_token, fid_slen, '\t');
	psam_info_reverse_ll = new_psam_info;
	if ((iid_slen == 1) && (token_ptrs[0][0] == '0')) {
	  sprintf(g_logbuf, "Error: Invalid IID '0' on line %" PRIuPTR " of %s.\n", line_idx, psamname);
	  goto load_psam_ret_MALFORMED_INPUT_WW;
	}
	ss_iter = memcpya(ss_iter, token_ptrs[0], iid_slen);
	const uint32_t sample_id_slen = (uintptr_t)(ss_iter - sample_id_storage);
	if (sample_id_slen >= max_sample_id_blen) {
	  max_sample_id_blen = sample_id_slen + 1;
	}
	new_psam_info->sample_id_slen = sample_id_slen;
	if (sids_present) {
	  ss_iter = memcpya(ss_iter, token_ptrs[1], sid_slen);
	  if (sid_slen >= max_sid_blen) {
	    max_sid_blen = sid_slen + 1;
	  }
	}
	if (paternal_ids_present) {
	  if (paternal_id_slen >= max_paternal_id_blen) {
	    max_paternal_id_blen = paternal_id_slen + 1;
	  }
	  ss_iter = memcpya(ss_iter, token_ptrs[2], paternal_id_slen);
	} else {
	  *ss_iter++ = '0';
	}
	new_psam_info->paternal_id_slen = paternal_id_slen;
	if (maternal_ids_present) {
	  if (maternal_id_slen >= max_maternal_id_blen) {
	    max_maternal_id_blen = maternal_id_slen + 1;
	  }
	  ss_iter = memcpya(ss_iter, token_ptrs[3], maternal_id_slen);
	} else {
	  *ss_iter++ = '0';
	}
	new_psam_info->maternal_id_slen = maternal_id_slen;
	uint32_t cur_sex_code = 0;
	// accept 'M'/'F'/'m'/'f' since that's more readable without being any
	// less efficient
	// don't accept "male"/"female", that's overkill
	if (sex_present && (token_slens[4] == 1)) {
	  const unsigned char sex_ucc = token_ptrs[4][0];
	  const unsigned char sex_ucc_upcase = sex_ucc & 0xdfU;
	  if ((sex_ucc == '1') || (sex_ucc_upcase == 'M')) {
	    cur_sex_code = 1;
	  } else if ((sex_ucc == '2') || (sex_ucc_upcase == 'F')) {
	    cur_sex_code = 2;
	  }
	}
	new_psam_info->sex_code = cur_sex_code;
	// phenotypes
	unsigned char* pheno_data = new_psam_info->vardata;
	uint32_t cat_pheno_idx = 0;
	for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  const uint32_t col_type_idx = pheno_idx + 5;
	  char* cur_phenostr = token_ptrs[col_type_idx];
	  double dxx;
	  if (!scanadv_double(cur_phenostr, &dxx)) {
	    const uint32_t slen = token_slens[col_type_idx];
	    if (is_nan_str(cur_phenostr, slen)) {
	      dxx = missing_phenod;
	    } else {
	      if (!IS_SET(categorical_phenos, pheno_idx)) {
		goto load_psam_ret_INCOMPATIBLE_PHENOSTRS;
	      }
	      if (slen > kMaxIdSlen) {
                logerrprint("Error: Categorical phenotypes are limited to " MAX_ID_SLEN_STR " characters.\n");
		goto load_psam_ret_MALFORMED_INPUT;
	      }
	      uint32_t hashval = hashceil(cur_phenostr, slen, kCatHtableSize) + cat_pheno_idx;
	      if (hashval >= kCatHtableSize) {
		hashval -= kCatHtableSize;
	      }
	      uintptr_t htable_idx = 0;
	      catname_ll2_t** cur_entry_ptr = &(catname_htable[hashval]);
	      while (1) {
		catname_ll2_t* cur_entry = *cur_entry_ptr;
		if (!cur_entry) {
		  const uint32_t entry_byte_ct = round_up_pow2(offsetof(catname_ll2_t, ss) + slen + 1, sizeof(intptr_t));
		  htable_idx = pheno_catname_last[cat_pheno_idx]->cat_idx + 1;
		  catname_ll2_t* new_entry = (catname_ll2_t*)tmp_bigstack_base;
		  tmp_bigstack_base += entry_byte_ct;
		  if (tmp_bigstack_base > tmp_bigstack_end) {
		    goto load_psam_ret_NOMEM;
		  }
		  new_entry->htable_next = nullptr;
		  new_entry->pheno_next = nullptr;
		  pheno_catname_last[cat_pheno_idx]->pheno_next = new_entry;
		  pheno_catname_last[cat_pheno_idx] = new_entry;
		  *cur_entry_ptr = new_entry;
		  new_entry->cat_idx = htable_idx;
		  memcpyx(new_entry->ss, cur_phenostr, slen, '\0');
		  total_catname_blens[cat_pheno_idx] += slen + 1;
		  break;
		}
		// safe since we guarantee kMaxIdSlen spare bytes at the end
		// of bigstack
		if ((!memcmp(cur_entry->ss, cur_phenostr, slen)) && (!cur_entry->ss[slen])) {
		  htable_idx = cur_entry->cat_idx;
		  break;
		}
		cur_entry_ptr = &(cur_entry->htable_next);
	      }
	      // don't bother writing top 4 bytes in 32-bit build
	      memcpy(&(pheno_data[pheno_idx * 8]), &htable_idx, sizeof(intptr_t));
	      ++cat_pheno_idx;
	      continue;
	    }
	  }
	  if (IS_SET(categorical_phenos, pheno_idx)) {
	    goto load_psam_ret_INCOMPATIBLE_PHENOSTRS;
	  }
	  if (!IS_SET(quantitative_phenos, pheno_idx)) {
	    if ((dxx != missing_phenod) && (dxx != pheno_ctrld) && (dxx != pheno_cased) && (dxx != 0.0)) {
	      SET_BIT(pheno_idx, quantitative_phenos);
	    }
	  }
	  memcpy(&(pheno_data[pheno_idx * 8]), &dxx, sizeof(double));
	}
	++raw_sample_ct;
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_psam_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto load_psam_ret_LONG_LINE;
	}
	goto load_psam_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, psamname);
	goto load_psam_ret_MALFORMED_INPUT_WW;
      }
    }
    if ((max_sample_id_blen > 2 * kMaxIdBlen) || (max_paternal_id_blen > kMaxIdBlen) || (max_maternal_id_blen > kMaxIdBlen)) {
      logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto load_psam_ret_MALFORMED_INPUT;
    }
    if (gzclose_null(&gz_infile)) {
      goto load_psam_ret_READ_FAIL;
    }
    const uintptr_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    g_bigstack_base = (unsigned char*)round_up_pow2((uintptr_t)tmp_bigstack_base, kCacheline);
    if (pheno_ct) {
      pheno_cols = (pheno_col_t*)malloc(pheno_ct * sizeof(pheno_col_t));
      if (!pheno_cols) {
	goto load_psam_ret_NOMEM;
      }
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	// ensure cleanup works if initialization fails in the middle
	pheno_cols[pheno_idx].nonmiss = nullptr;
      }
      uint32_t cat_pheno_idx = 0;
      pheno_col_t* pheno_cols_iter = pheno_cols;
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	const uintptr_t nonmiss_vec_ct = BITCT_TO_VECCT(raw_sample_ct);
	const uint32_t is_categorical = IS_SET(categorical_phenos, pheno_idx);
	const uint32_t is_qt = IS_SET(quantitative_phenos, pheno_idx);
	uintptr_t data_vec_ct = 0;
	uintptr_t catname_vec_ct = 0;
	uintptr_t catname_storage_vec_ct = 0;
	uint32_t nonnull_catname_ct = 0;
	if (!is_categorical) {
	  pheno_cols_iter->category_names = nullptr;
	  pheno_cols_iter->type_code = (pheno_dtype_t)is_qt;
	  pheno_cols_iter->nonnull_category_ct = 0;
	  if (is_qt) {
	    data_vec_ct = DBLCT_TO_VECCT(raw_sample_ct);
	  } else {
	    data_vec_ct = nonmiss_vec_ct;
	  }
	} else {
	  nonnull_catname_ct = pheno_catname_last[cat_pheno_idx]->cat_idx;
	  data_vec_ct = INT32CT_TO_VECCT(raw_sample_ct);
	  catname_vec_ct = WORDCT_TO_VECCT(nonnull_catname_ct + 1);
	  catname_storage_vec_ct = DIV_UP(total_catname_blens[cat_pheno_idx], kBytesPerVec);
	  pheno_cols_iter->type_code = kPhenoDtypeCat;
	  pheno_cols_iter->nonnull_category_ct = nonnull_catname_ct;
	}
	// pheno_cols_iter->nonmiss = nullptr;
	uintptr_t* new_pheno_data_iter;
	if (vecaligned_malloc((nonmiss_vec_ct + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &new_pheno_data_iter)) {
	  goto load_psam_ret_NOMEM;
	}
	pheno_cols_iter->nonmiss = new_pheno_data_iter;
	fill_ulong_zero(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
	new_pheno_data_iter = &(new_pheno_data_iter[nonmiss_vec_ct * kWordsPerVec]);
	if (is_categorical) {
	  pheno_cols_iter->data.cat = (uint32_t*)new_pheno_data_iter;
	  new_pheno_data_iter = &(new_pheno_data_iter[data_vec_ct * kWordsPerVec]);
	  char** cur_name_ptrs = (char**)new_pheno_data_iter;
	  pheno_cols_iter->category_names = cur_name_ptrs;
	  *cur_name_ptrs++ = missing_catname;
	  char* name_storage_iter = (char*)(&(new_pheno_data_iter[catname_vec_ct * kWordsPerVec]));
	  uint32_t cur_hval = missing_catname_hval + cat_pheno_idx;
	  if (cur_hval >= kCatHtableSize) {
	    cur_hval -= kCatHtableSize;
	  }
	  // make this point to the "NONE" entry for the current phenotype,
	  // which starts the linked list
	  catname_ll2_t* catname_entry_ptr = catname_htable[cur_hval];
	  
	  for (uint32_t catname_idx = 0; catname_idx < nonnull_catname_ct; ++catname_idx) {
	    catname_entry_ptr = catname_entry_ptr->pheno_next;
	    *cur_name_ptrs++ = name_storage_iter;
	    name_storage_iter = strcpyax(name_storage_iter, catname_entry_ptr->ss, '\0');
	  }
	  ++cat_pheno_idx;
	} else if (!is_qt) {
	  pheno_cols_iter->data.cc = new_pheno_data_iter;
	  fill_ulong_zero(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
	} else {
	  pheno_cols_iter->data.qt = (double*)new_pheno_data_iter;
	}
	++pheno_cols_iter;
      }
    }
    // real allocations start here
    // could make these cacheline-aligned?
    g_bigstack_end = bigstack_end_mark;
    const uint32_t aligned_wct = BITCT_TO_ALIGNED_WORDCT(raw_sample_ct);
    if (bigstack_end_alloc_c(raw_sample_ct * max_sample_id_blen, sample_ids_ptr) ||
	bigstack_end_alloc_c(raw_sample_ct * max_paternal_id_blen, paternal_ids_ptr) ||
	bigstack_end_alloc_c(raw_sample_ct * max_maternal_id_blen, maternal_ids_ptr) ||
	bigstack_end_alloc_ul(raw_sample_ctl, sample_include_ptr) ||
	bigstack_end_calloc_ul(aligned_wct, founder_info_ptr) ||
	bigstack_end_calloc_ul(aligned_wct, sex_nm_ptr) ||
	bigstack_end_calloc_ul(aligned_wct, sex_male_ptr)) {
      goto load_psam_ret_NOMEM;
    }
    if (sids_present) {
      if (bigstack_end_alloc_c(raw_sample_ct * max_sid_blen, sids_ptr)) {
	goto load_psam_ret_NOMEM;
      }
    }
    bigstack_end_mark = g_bigstack_end;
    fill_all_bits(raw_sample_ct, *sample_include_ptr);
    // make fill_interleaved_mask_vec() work by default
    fill_ulong_zero(aligned_wct - raw_sample_ctl, &((*sample_include_ptr)[raw_sample_ctl]));
    *raw_sample_ct_ptr = raw_sample_ct;
    *max_sample_id_blen_ptr = max_sample_id_blen;
    *max_sid_blen_ptr = max_sid_blen;
    *max_paternal_id_blen_ptr = max_paternal_id_blen;
    *max_maternal_id_blen_ptr = max_maternal_id_blen;
    *max_pheno_name_blen_ptr = max_pheno_name_blen;
    char* sample_ids = *sample_ids_ptr;
    char* sids = sids_present? (*sids_ptr) : nullptr;
    char* paternal_ids = *paternal_ids_ptr;
    char* maternal_ids = *maternal_ids_ptr;
    uintptr_t* founder_info = *founder_info_ptr;
    uintptr_t* sex_nm = *sex_nm_ptr;
    uintptr_t* sex_male = *sex_male_ptr;
    uint32_t sample_uidx = raw_sample_ct;
    while (sample_uidx) {
      --sample_uidx;
      unsigned char* cur_vardata = psam_info_reverse_ll->vardata;
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	if (IS_SET(categorical_phenos, pheno_idx)) {
	  uint32_t cur_cat;
	  memcpy(&cur_cat, &(cur_vardata[pheno_idx * 8]), sizeof(int32_t));
	  pheno_cols[pheno_idx].data.cat[sample_uidx] = cur_cat;
	  if (cur_cat) {
	    SET_BIT(sample_uidx, pheno_cols[pheno_idx].nonmiss);
	  }
	} else {
	  double dxx;
	  memcpy(&dxx, &(cur_vardata[pheno_idx * 8]), sizeof(double));
	  if (IS_SET(quantitative_phenos, pheno_idx)) {
	    if (dxx != missing_phenod) {
	      SET_BIT(sample_uidx, pheno_cols[pheno_idx].nonmiss);
	      pheno_cols[pheno_idx].data.qt[sample_uidx] = dxx;
	    }
	  } else {
	    if (dxx == pheno_cased) {
	      SET_BIT(sample_uidx, pheno_cols[pheno_idx].data.cc);
	      SET_BIT(sample_uidx, pheno_cols[pheno_idx].nonmiss);
	    } else if (dxx == pheno_ctrld) {
	      SET_BIT(sample_uidx, pheno_cols[pheno_idx].nonmiss);
	    }
	  }
	}
      }
      const uint32_t sample_id_slen = psam_info_reverse_ll->sample_id_slen;
      const uint32_t paternal_id_slen = psam_info_reverse_ll->paternal_id_slen;
      const uint32_t maternal_id_slen = psam_info_reverse_ll->maternal_id_slen;
      const uint32_t sex_code = psam_info_reverse_ll->sex_code;
      char* cur_sample_id = (char*)(&(cur_vardata[pheno_ct * 8]));
      memcpyx(&(sample_ids[sample_uidx * max_sample_id_blen]), cur_sample_id, sample_id_slen, '\0');
      char* cur_paternal_id = &(cur_sample_id[sample_id_slen]);
      if (sids) {
        const uint32_t sid_slen = psam_info_reverse_ll->sid_slen;
	char* cur_sid = cur_paternal_id;
	memcpyx(&(sids[sample_uidx * max_sid_blen]), cur_sid, sid_slen, '\0');
	cur_paternal_id = &(cur_sid[sid_slen]);
      }
      memcpyx(&(paternal_ids[sample_uidx * max_paternal_id_blen]), cur_paternal_id, paternal_id_slen, '\0');
      char* cur_maternal_id = &(cur_paternal_id[paternal_id_slen]);
      if ((paternal_id_slen == 1) && (maternal_id_slen == 1) && (cur_paternal_id[0] == '0') && (cur_maternal_id[0] == '0')) {
	SET_BIT(sample_uidx, founder_info);
      }
      memcpyx(&(maternal_ids[sample_uidx * max_maternal_id_blen]), cur_maternal_id, maternal_id_slen, '\0');
      if (sex_code) {
	SET_BIT(sample_uidx, sex_nm);
	if (sex_code == 1) {
	  SET_BIT(sample_uidx, sex_male);
	}
      }
      psam_info_reverse_ll = psam_info_reverse_ll->next;
    }
    // special case: if there's exactly one phenotype and it's all-missing,
    // discard it
    if ((pheno_ct == 1) && (!popcount_longs(pheno_cols[0].nonmiss, raw_sample_ctl))) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
      cleanup_pheno_cols(1, pheno_cols);
      *pheno_cols_ptr = nullptr;
      *pheno_ct_ptr = 0;
    } else {
      *pheno_cols_ptr = pheno_cols;
      *pheno_ct_ptr = pheno_ct;
    }
  }
  while (0) {
  load_psam_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_psam_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_psam_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  load_psam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  load_psam_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, psamname);
  load_psam_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  load_psam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  load_psam_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, psamname);
    reterr = kPglRetMalformedInput;
    break;
  load_psam_ret_INCOMPATIBLE_PHENOSTRS:
    LOGERRPRINTFWW("Error: Incompatible phenotype values in %s. (Case/control and quantitative phenotypes must be entirely numeric/\"NA\", and categorical phenotypes must be entirely non-numeric.)\n", psamname);
    reterr = kPglRetMalformedInput;
    break;
  }
 load_psam_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  gzclose_cond(gz_infile);
  if (reterr) {
    if (*pheno_names_ptr) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
    }
    cleanup_pheno_cols(pheno_ct, pheno_cols);
    *pheno_ct_ptr = 0;
    *pheno_cols_ptr = nullptr;
  }
  return reterr;
}


typedef struct pheno_info_ll_struct {
  // for categorical phenotypes, phenodata entry should be reinterpreted as
  // uint32_t
  struct pheno_info_ll_struct* next;
  uint32_t sample_uidx;
  double phenodata[];
} pheno_info_ll_t;

// also for loading covariates.  set affection_01 to 2 to prohibit case/control
// and make unnamed variables start with "COVAR" instead of "PHENO"
pglerr_t load_phenos(const char* pheno_fname, const range_list_t* pheno_range_list_ptr, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, int32_t missing_pheno, uint32_t affection_01, pheno_col_t** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  
  gzFile gz_infile = nullptr;
  char* pheno_names = nullptr;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    reterr = gzopen_read_checked(pheno_fname, &gz_infile);
    if (reterr) {
      goto load_phenos_ret_1;
    }
    const uintptr_t initial_bigstack_size = bigstack_left();
    uintptr_t loadbuf_size = initial_bigstack_size / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto load_phenos_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_first_token;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_phenos_ret_READ_FAIL;
	}
	loadbuf_first_token = loadbuf;
	loadbuf_first_token[0] = '\0';
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto load_phenos_ret_LONG_LINE;
	}
	goto load_phenos_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token) || ((loadbuf_first_token[0] == '#') && (((loadbuf_first_token[1] != 'F') && (loadbuf_first_token[1] != 'I')) || memcmp(&(loadbuf_first_token[2]), "ID", 2) || (((unsigned char)loadbuf_first_token[4] > ' ') && (loadbuf_first_token[4] != ',')))));
    if (loadbuf_first_token[0] == '#') {
      ++loadbuf_first_token;
    }    
    const uint32_t old_pheno_ct = *pheno_ct_ptr;
    const uintptr_t old_max_pheno_name_blen = *max_pheno_name_blen_ptr;
    uintptr_t max_pheno_name_blen = old_max_pheno_name_blen;
    uint32_t comma_delim = 0;
    xid_mode_t xid_mode = (loadbuf_first_token[0] == 'I')? kfXidModeIid : kfXidModeFidiid;
    uint32_t* col_skips = nullptr;
    uint32_t new_pheno_ct;
    uint32_t final_pheno_ct;
    uintptr_t final_pheno_names_byte_ct;
    if (((loadbuf_first_token[0] == 'F') || xid_mode) && (!memcmp(&(loadbuf_first_token[1]), "ID", 2)) && (((unsigned char)loadbuf_first_token[3] <= 32) || (loadbuf_first_token[3] == ','))) {
      // treat this as a header line
      // autodetect CSV vs. space/tab-delimited
      // (note that we don't permit CSVs without header lines)
      comma_delim = (loadbuf_first_token[3] == ',');
      char* loadbuf_iter = skip_initial_spaces(&(loadbuf_first_token[3 + comma_delim]));
      if (is_eoln_kns(*loadbuf_iter)) {
	goto load_phenos_ret_MISSING_TOKENS;
      }
      char* iid_end;
      if (!xid_mode) {
        iid_end = comma_or_space_token_end(loadbuf_iter, comma_delim);
	const uintptr_t token_slen = (uintptr_t)(iid_end - loadbuf_iter);
	if ((token_slen != 3) || memcmp(loadbuf_iter, "IID", 3)) {
	  sprintf(g_logbuf, "Error: Second column header in %s must be 'IID'.\n", pheno_fname);
	  goto load_phenos_ret_MALFORMED_INPUT_WW;
	}
	loadbuf_iter = comma_or_space_next_token(iid_end, comma_delim);
      } else {
	iid_end = &(loadbuf_first_token[3]);
      }
      uint32_t pheno_col_ct = 0;
      char* pheno_start = loadbuf_iter;
      while (loadbuf_iter) {
	char* token_end = comma_or_space_token_end(loadbuf_iter, comma_delim);
	const uintptr_t token_slen = (uintptr_t)(token_end - loadbuf_iter);
	if (max_pheno_name_blen <= token_slen) {
	  max_pheno_name_blen = token_slen + 1;
	}
	++pheno_col_ct;
	loadbuf_iter = comma_or_space_next_token(token_end, comma_delim);
      }
      if (max_pheno_name_blen > kMaxIdBlen) {
	logerrprint("Error: Phenotype/covariate names are limited to " MAX_ID_SLEN_STR " characters.\n");
	goto load_phenos_ret_MALFORMED_INPUT;
      }
      if (pheno_range_list_ptr->names && pheno_col_ct) {
	uintptr_t* bitarr;
	reterr = string_range_list_to_bitarr_alloc(pheno_start, pheno_range_list_ptr, "pheno-name", "--pheno file", pheno_col_ct, 0, comma_delim, &bitarr);
	if (reterr) {
	  goto load_phenos_ret_1;
	}
	new_pheno_ct = popcount_longs(bitarr, BITCT_TO_WORDCT(pheno_col_ct));
	if (bigstack_alloc_ui(new_pheno_ct, &col_skips)) {
	  goto load_phenos_ret_NOMEM;
	}
	uint32_t col_uidx = 0;
	int32_t prev_col_uidx = -1;
	for (uint32_t col_idx = 0; col_idx < new_pheno_ct; ++col_idx, ++col_uidx) {
	  next_set_unsafe_ck(bitarr, &col_uidx);
	  col_skips[col_idx] = ((int32_t)col_uidx) - prev_col_uidx;
	  prev_col_uidx = (int32_t)col_uidx;
	}
      } else {
	// usual case, load all phenotypes
	new_pheno_ct = pheno_col_ct;
	if (bigstack_alloc_ui(new_pheno_ct, &col_skips)) {
	  goto load_phenos_ret_NOMEM;
	}
	for (uint32_t col_idx = 0; col_idx < pheno_col_ct; ++col_idx) {
	  col_skips[col_idx] = 1;
	}
      }
      final_pheno_ct = new_pheno_ct + old_pheno_ct;
      final_pheno_names_byte_ct = final_pheno_ct * max_pheno_name_blen;
      pheno_names = (char*)malloc(final_pheno_names_byte_ct);
      if (!pheno_names) {
	goto load_phenos_ret_NOMEM;
      }
      loadbuf_iter = iid_end;
      char* pheno_names_iter = &(pheno_names[old_pheno_ct * max_pheno_name_blen]);
      for (uint32_t new_pheno_idx = 0; new_pheno_idx < new_pheno_ct; ++new_pheno_idx) {
	loadbuf_iter = comma_or_space_next_token_mult(loadbuf_iter, col_skips[new_pheno_idx], comma_delim);
	char* token_end = comma_or_space_token_end(loadbuf_iter, comma_delim);
	const uint32_t name_slen = (uintptr_t)(token_end - loadbuf_iter);
	if (is_reserved_pheno_name(loadbuf_iter, name_slen)) {
	  *token_end = '\0';
	  sprintf(g_logbuf, "Error: '%s' cannot be used as a phenotype/covariate name.\n", loadbuf_iter);
	  goto load_phenos_ret_MALFORMED_INPUT_2;
	}
	memcpyx(pheno_names_iter, loadbuf_iter, name_slen, '\0');
	pheno_names_iter = &(pheno_names_iter[max_pheno_name_blen]);
	loadbuf_iter = token_end;
      }
      
      // forces line to be skipped by main loop
      loadbuf_first_token[0] = '\0';
    } else {
      // no header line
      xid_mode = kfXidModeFidiid;
      if (pheno_range_list_ptr->names) {
	// possible todo: support e.g. "PHENO2-PHENO10"
	sprintf(g_logbuf, "Error: Header line expected in %s (due to --pheno-name/--covar-name).\n", pheno_fname);
	goto load_phenos_ret_INCONSISTENT_INPUT_WW;
      }
      const uint32_t col_ct = count_tokens(loadbuf_first_token);
      if (col_ct < 3) {
	// todo: tolerate col_ct == 2 with --allow-no-phenos
	goto load_phenos_ret_MISSING_TOKENS;
      }
      new_pheno_ct = col_ct - 2;
      final_pheno_ct = new_pheno_ct + old_pheno_ct;
      const uintptr_t max_new_name_blen = 6 + int_slen(final_pheno_ct - 1);
      if (max_new_name_blen > max_pheno_name_blen) {
	max_pheno_name_blen = max_new_name_blen;
      }
      final_pheno_names_byte_ct = final_pheno_ct * max_pheno_name_blen;
      if (bigstack_alloc_ui(new_pheno_ct, &col_skips)) {
	goto load_phenos_ret_NOMEM;
      }
      pheno_names = (char*)malloc(final_pheno_names_byte_ct);
      if (!pheno_names) {
	goto load_phenos_ret_NOMEM;
      }
      for (uint32_t col_idx = 0; col_idx < new_pheno_ct; ++col_idx) {
	col_skips[col_idx] = 1;
      }
      const char* default_prefix = (affection_01 == 2)? "COVAR" : "PHENO";
      for (uint32_t pheno_idx = old_pheno_ct; pheno_idx < final_pheno_ct;) {
	char* write_iter = memcpya(&(pheno_names[pheno_idx * max_pheno_name_blen]), default_prefix, 5);
	++pheno_idx; // 1-based default names, not 0-based
	write_iter = uint32toa(pheno_idx, write_iter);
	*write_iter = '\0';
      }
    }
    for (uint32_t old_pheno_idx = 0; old_pheno_idx < old_pheno_ct; ++old_pheno_idx) {
      strcpy(&(pheno_names[old_pheno_idx * max_pheno_name_blen]), &((*pheno_names_ptr)[old_pheno_idx * old_max_pheno_name_blen]));
    }

    uint32_t tmp_htable_size;
    uint32_t* htable_tmp;
    if (htable_good_size_alloc(final_pheno_ct, bigstack_left(), &htable_tmp, &tmp_htable_size)) {
      goto load_phenos_ret_NOMEM;
    }
    const uint32_t duplicate_idx = populate_strbox_htable(pheno_names, final_pheno_ct, max_pheno_name_blen, tmp_htable_size, htable_tmp);
    if (duplicate_idx) {
      const char* duplicate_pheno_name = &(pheno_names[duplicate_idx * max_pheno_name_blen]);
      sprintf(g_logbuf, "Error: Duplicate phenotype/covariate ID '%s'.\n", duplicate_pheno_name);
      goto load_phenos_ret_MALFORMED_INPUT_WW;
    }
    bigstack_reset(htable_tmp);

    pheno_col_t* new_pheno_cols = (pheno_col_t*)realloc(*pheno_cols_ptr, final_pheno_ct * sizeof(pheno_col_t));
    if (!new_pheno_cols) {
      goto load_phenos_ret_NOMEM;
    }
    // ensure cleanup works if initialization fails in the middle
    for (uint32_t pheno_idx = old_pheno_ct; pheno_idx < final_pheno_ct; ++pheno_idx) {
      new_pheno_cols[pheno_idx].nonmiss = nullptr;
    }
    *pheno_ct_ptr = final_pheno_ct;
    *pheno_cols_ptr = new_pheno_cols;

    // switch to hash table?
    char* sorted_sample_ids;
    uint32_t* sample_id_map;
    // todo: permit duplicates if SIDs are defined
    reterr = copy_sort_strbox_subset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 0, 0, 0, &sorted_sample_ids, &sample_id_map);
    if (reterr) {
      goto load_phenos_ret_1;
    }
    const uintptr_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    char* id_buf;
    uintptr_t* already_seen;
    if (bigstack_alloc_c(max_sample_id_blen, &id_buf) ||
	bigstack_calloc_ul(raw_sample_ctl, &already_seen)) {
      goto load_phenos_ret_NOMEM;
    }

    pheno_info_ll_t* pheno_info_reverse_ll = nullptr;
    const uint32_t pheno_info_alloc_byte_ct = sizeof(pheno_info_ll_t) + new_pheno_ct * sizeof(double);
    const uint32_t new_pheno_ctl = BITCT_TO_WORDCT(new_pheno_ct);
    const double missing_phenod = (double)missing_pheno;
    const double pheno_ctrld = (double)((int32_t)(1 - affection_01));
    const double pheno_cased = pheno_ctrld + 1.0;
    uint32_t categorical_pheno_ct = 0;
    char** token_ptrs;
    uint32_t* token_slens;
    uintptr_t* categorical_phenos;
    uintptr_t* quantitative_phenos;
    if (bigstack_alloc_cp(new_pheno_ct, &token_ptrs) ||
	bigstack_alloc_ui(new_pheno_ct, &token_slens) ||
	bigstack_calloc_ul(new_pheno_ctl, &categorical_phenos) ||
	bigstack_calloc_ul(new_pheno_ctl, &quantitative_phenos)) {
      goto load_phenos_ret_NOMEM;
    }
    const char* missing_catname = g_missing_catname;
    const uint32_t missing_catname_blen = strlen(missing_catname) + 1;
    const uint32_t missing_catname_hval = hashceil(missing_catname, missing_catname_blen - 1, kCatHtableSize);
    unsigned char* bigstack_base_copy = g_bigstack_base;
    unsigned char* tmp_bigstack_end = g_bigstack_end;
    catname_ll2_t** catname_htable = nullptr;
    catname_ll2_t** pheno_catname_last = nullptr;
    uintptr_t* total_catname_blens = nullptr;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	char* loadbuf_iter = loadbuf_first_token;
	uint32_t sample_uidx;
	if (sorted_xidbox_read_find(sorted_sample_ids, sample_id_map, max_sample_id_blen, sample_ct, comma_delim, xid_mode, &loadbuf_iter, &sample_uidx, id_buf)) {
	  if (!loadbuf_iter) {
	    goto load_phenos_ret_MISSING_TOKENS;
	  }
	} else {
	  if (is_set(already_seen, sample_uidx)) {
	    logerrprint("Error: Duplicate sample ID in --pheno/--covar file.\n");
	    goto load_phenos_ret_MALFORMED_INPUT;
	  }
	  set_bit(sample_uidx, already_seen);
	  /*
	  if (comma_delim) {
	  } else {
	  }
	  */
	  for (uint32_t new_pheno_idx = 0; new_pheno_idx < new_pheno_ct; ++new_pheno_idx) {
	    loadbuf_iter = comma_or_space_next_token_mult(loadbuf_iter, col_skips[new_pheno_idx], comma_delim);
	    if (!loadbuf_iter) {
	      goto load_phenos_ret_MISSING_TOKENS;
	    }
	    token_ptrs[new_pheno_idx] = loadbuf_iter;
	    char* token_end = comma_or_space_token_end(loadbuf_iter, comma_delim);
	    token_slens[new_pheno_idx] = (uintptr_t)(token_end - loadbuf_iter);
	    loadbuf_iter = token_end;
	  }
	  if (!pheno_info_reverse_ll) {
	    // first relevant line, detect categorical phenotypes
	    for (uint32_t new_pheno_idx = 0; new_pheno_idx < new_pheno_ct; ++new_pheno_idx) {
	      if (is_categorical_phenostr(token_ptrs[new_pheno_idx])) {
		SET_BIT(new_pheno_idx, categorical_phenos);
	      } else if (affection_01 == 2) {
		SET_BIT(new_pheno_idx, quantitative_phenos);
	      }
	    }
	    categorical_pheno_ct = popcount_longs(categorical_phenos, new_pheno_ctl);
	    if (categorical_pheno_ct) {
	      // initialize hash table
	      if (categorical_pheno_ct > kCatHtableSize) {
		// use a larger hash table if/when we ever care for this case
		logerrprint("Error: " PROG_NAME_STR " does not support more than 2^19 - 1 categorical phenotypes.\n");
		goto load_phenos_ret_MALFORMED_INPUT;
	      }
	      const uint32_t cat_ul_byte_ct = categorical_pheno_ct * sizeof(intptr_t);
	      const uint32_t htable_byte_ct = kCatHtableSize * sizeof(uintptr_t);
	      const uintptr_t entry_byte_ct = round_up_pow2(offsetof(catname_ll2_t, ss) + missing_catname_blen, sizeof(intptr_t));

	      if ((uintptr_t)(tmp_bigstack_end - bigstack_base_copy) < htable_byte_ct + categorical_pheno_ct * entry_byte_ct + 2 * cat_ul_byte_ct) {
		goto load_phenos_ret_NOMEM;
	      }
	      tmp_bigstack_end -= cat_ul_byte_ct;
	      total_catname_blens = (uintptr_t*)tmp_bigstack_end;
	      tmp_bigstack_end -= cat_ul_byte_ct;
	      pheno_catname_last = (catname_ll2_t**)tmp_bigstack_end;
	      fill_ulong_zero(categorical_pheno_ct, total_catname_blens);
	      tmp_bigstack_end -= htable_byte_ct;
	      catname_htable = (catname_ll2_t**)tmp_bigstack_end;
	      for (uint32_t uii = 0; uii < kCatHtableSize; ++uii) {
		catname_htable[uii] = nullptr;
	      }
	      uint32_t cur_hval = missing_catname_hval;
	      for (uint32_t cat_pheno_idx = 0; cat_pheno_idx < categorical_pheno_ct; ++cat_pheno_idx) {
		tmp_bigstack_end -= entry_byte_ct;
		catname_ll2_t* new_entry = (catname_ll2_t*)tmp_bigstack_end;
		pheno_catname_last[cat_pheno_idx] = new_entry;
		new_entry->cat_idx = 0;
		new_entry->htable_next = nullptr;
		new_entry->pheno_next = nullptr;
		memcpy(new_entry->ss, missing_catname, missing_catname_blen);
		catname_htable[cur_hval++] = new_entry;
		if (cur_hval == kCatHtableSize) {
		  cur_hval = 0;
		}
	      }
	    }
	  }
	  if ((uintptr_t)(tmp_bigstack_end - bigstack_base_copy) < pheno_info_alloc_byte_ct) {
	    goto load_phenos_ret_NOMEM;
	  }
	  tmp_bigstack_end -= pheno_info_alloc_byte_ct;

	  pheno_info_ll_t* new_pheno_info = (pheno_info_ll_t*)tmp_bigstack_end;
	  new_pheno_info->next = pheno_info_reverse_ll;
	  new_pheno_info->sample_uidx = sample_uidx;
	  double* pheno_data = new_pheno_info->phenodata;
	  uint32_t cat_pheno_idx = 0;
	  for (uint32_t new_pheno_idx = 0; new_pheno_idx < new_pheno_ct; ++new_pheno_idx) {
	    char* cur_phenostr = token_ptrs[new_pheno_idx];
	    double dxx;
	    if (!scanadv_double(cur_phenostr, &dxx)) {
	      const uint32_t slen = token_slens[new_pheno_idx];
	      if (is_nan_str(cur_phenostr, slen)) {
		// note that, in CSVs, empty string is interpreted as a
		// missing non-categorical phenotype; explicit "NONE" is needed
		// to denote a missing category
		dxx = missing_phenod;
	      } else {
		if (!IS_SET(categorical_phenos, new_pheno_idx)) {
		  goto load_phenos_ret_INCOMPATIBLE_PHENOSTRS;
		}
		uint32_t hashval;
		hashval = hashceil(cur_phenostr, slen, kCatHtableSize) + cat_pheno_idx;
		if (hashval >= kCatHtableSize) {
		  hashval -= kCatHtableSize;
		}
		uintptr_t htable_idx = 0;
		catname_ll2_t** cur_entry_ptr = &(catname_htable[hashval]);
		while (1) {
		  catname_ll2_t* cur_entry = *cur_entry_ptr;
		  if (!cur_entry) {
		    const uint32_t entry_byte_ct = round_up_pow2(offsetof(catname_ll2_t, ss) + slen + 1, sizeof(intptr_t));
		    if ((uintptr_t)(tmp_bigstack_end - bigstack_base_copy) < entry_byte_ct) {
		      goto load_phenos_ret_NOMEM;
		    }
		    tmp_bigstack_end -= entry_byte_ct;
		    htable_idx = pheno_catname_last[cat_pheno_idx]->cat_idx + 1;
		    catname_ll2_t* new_entry = (catname_ll2_t*)tmp_bigstack_end;
		    new_entry->htable_next = nullptr;
		    new_entry->pheno_next = nullptr;
		    pheno_catname_last[cat_pheno_idx]->pheno_next = new_entry;
		    pheno_catname_last[cat_pheno_idx] = new_entry;
		    *cur_entry_ptr = new_entry;
		    new_entry->cat_idx = htable_idx;
		    memcpyx(new_entry->ss, cur_phenostr, slen, '\0');
		    total_catname_blens[cat_pheno_idx] += slen + 1;
		    break;
		  }
		  // safe since hash table entries are in the middle of
		  // bigstack
		  if ((!memcmp(cur_entry->ss, cur_phenostr, slen)) && (!cur_entry->ss[slen])) {
		    htable_idx = cur_entry->cat_idx;
		    break;
		  }
		  cur_entry_ptr = &(cur_entry->htable_next);
		}
		// don't bother writing top 4 bytes in 32-bit build
		memcpy(&(pheno_data[new_pheno_idx]), &htable_idx, sizeof(intptr_t));
		++cat_pheno_idx;
		continue;
	      }
	    }
	    if (IS_SET(categorical_phenos, new_pheno_idx)) {
	      goto load_phenos_ret_INCOMPATIBLE_PHENOSTRS;
	    }
	    if (!IS_SET(quantitative_phenos, new_pheno_idx)) {
	      if ((dxx != missing_phenod) && (dxx != pheno_ctrld) && (dxx != pheno_cased) && (dxx != 0.0)) {
		SET_BIT(new_pheno_idx, quantitative_phenos);
	      }
	    }
	    pheno_data[new_pheno_idx] = dxx;
	  }
	  pheno_info_reverse_ll = new_pheno_info;
	}
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_phenos_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto load_phenos_ret_LONG_LINE;
	}
	goto load_phenos_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, pheno_fname);
	goto load_phenos_ret_MALFORMED_INPUT_WW;
      }
    }    
    if (gzclose_null(&gz_infile)) {
      goto load_phenos_ret_READ_FAIL;
    }
    if (!pheno_info_reverse_ll) {
      if (line_idx == 1) {
	sprintf(g_logbuf, "Error: %s is empty.\n", pheno_fname);
	goto load_phenos_ret_MALFORMED_INPUT_WW;
      }
      // could make this a warning, and automatically delete phenotypes?
      LOGERRPRINTF("Error: No entries in %s correspond to loaded sample IDs.\n", pheno_fname);
      goto load_phenos_ret_INCONSISTENT_INPUT;
    }
    if (new_pheno_ct) {
      const uintptr_t nonmiss_vec_ct = BITCT_TO_VECCT(raw_sample_ct);
      uint32_t cat_pheno_idx = 0;
      pheno_col_t* pheno_cols_iter = &(new_pheno_cols[old_pheno_ct]);
      for (uint32_t new_pheno_idx = 0; new_pheno_idx < new_pheno_ct; ++new_pheno_idx) {
	const uint32_t is_categorical = IS_SET(categorical_phenos, new_pheno_idx);
	const uint32_t is_qt = IS_SET(quantitative_phenos, new_pheno_idx);
	uintptr_t data_vec_ct = 0;
	uintptr_t catname_vec_ct = 0;
	uintptr_t catname_storage_vec_ct = 0;
	uint32_t nonnull_catname_ct = 0;
	if (!is_categorical) {
	  pheno_cols_iter->category_names = nullptr;
	  pheno_cols_iter->type_code = (pheno_dtype_t)is_qt;
	  pheno_cols_iter->nonnull_category_ct = 0;
	  if (is_qt) {
	    data_vec_ct = DBLCT_TO_VECCT(raw_sample_ct);
	  } else {
	    data_vec_ct = nonmiss_vec_ct;
	  }
	} else {
	  nonnull_catname_ct = pheno_catname_last[cat_pheno_idx]->cat_idx;
	  data_vec_ct = INT32CT_TO_VECCT(raw_sample_ct);
	  catname_vec_ct = WORDCT_TO_VECCT(nonnull_catname_ct + 1);
	  catname_storage_vec_ct = DIV_UP(total_catname_blens[cat_pheno_idx], kBytesPerVec);
	  pheno_cols_iter->type_code = kPhenoDtypeCat;
	  pheno_cols_iter->nonnull_category_ct = nonnull_catname_ct;
	}
	// pheno_cols_iter->nonmiss = nullptr;
	uintptr_t* new_pheno_data_iter;
	if (vecaligned_malloc((nonmiss_vec_ct + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &new_pheno_data_iter)) {
	  goto load_phenos_ret_NOMEM;
	}
	pheno_cols_iter->nonmiss = new_pheno_data_iter;
	fill_ulong_zero(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
	new_pheno_data_iter = &(new_pheno_data_iter[nonmiss_vec_ct * kWordsPerVec]);
	if (is_categorical) {
	  // allow nonmiss[] to be ignored in categorical case
	  fill_ulong_zero(data_vec_ct, new_pheno_data_iter);
	  pheno_cols_iter->data.cat = (uint32_t*)new_pheno_data_iter;
	  new_pheno_data_iter = &(new_pheno_data_iter[data_vec_ct * kWordsPerVec]);
	  char** cur_name_ptrs = (char**)new_pheno_data_iter;
	  pheno_cols_iter->category_names = cur_name_ptrs;
	  *cur_name_ptrs++ = g_missing_catname;
	  char* name_storage_iter = (char*)(&(new_pheno_data_iter[catname_vec_ct * kWordsPerVec]));
	  uint32_t cur_hval = missing_catname_hval + cat_pheno_idx;
	  if (cur_hval >= kCatHtableSize) {
	    cur_hval -= kCatHtableSize;
	  }
	  // make this point to the "NONE" entry for the current phenotype,
	  // which starts the linked list
	  catname_ll2_t* catname_entry_ptr = catname_htable[cur_hval];
	  
	  for (uint32_t catname_idx = 0; catname_idx < nonnull_catname_ct; ++catname_idx) {
	    catname_entry_ptr = catname_entry_ptr->pheno_next;
	    *cur_name_ptrs++ = name_storage_iter;
	    name_storage_iter = strcpyax(name_storage_iter, catname_entry_ptr->ss, '\0');
	  }
	  ++cat_pheno_idx;
	} else if (!is_qt) {
	  pheno_cols_iter->data.cc = new_pheno_data_iter;
	  fill_ulong_zero(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
	} else {
	  pheno_cols_iter->data.qt = (double*)new_pheno_data_iter;
	}
	++pheno_cols_iter;
      }
      while (pheno_info_reverse_ll) {
	const uint32_t sample_uidx = pheno_info_reverse_ll->sample_uidx;
	double* pheno_data = pheno_info_reverse_ll->phenodata;
	pheno_cols_iter = &(new_pheno_cols[old_pheno_ct]);
	for (uint32_t new_pheno_idx = 0; new_pheno_idx < new_pheno_ct; ++new_pheno_idx) {
	  if (IS_SET(categorical_phenos, new_pheno_idx)) {
	    uint32_t cur_cat;
	    memcpy(&cur_cat, &(pheno_data[new_pheno_idx]), sizeof(int32_t));
	    pheno_cols_iter->data.cat[sample_uidx] = cur_cat;
	    if (cur_cat) {
	      SET_BIT(sample_uidx, pheno_cols_iter->nonmiss);
	    }
	  } else {
	    double dxx = pheno_data[new_pheno_idx];
	    // bugfix (6 May 2017): forgot to accept 0 as missing value for
	    // case/control
	    if (IS_SET(quantitative_phenos, new_pheno_idx)) {
	      if (dxx != missing_phenod) {
		SET_BIT(sample_uidx, pheno_cols_iter->nonmiss);
		pheno_cols_iter->data.qt[sample_uidx] = dxx;
	      }
	    } else {
	      if (dxx == pheno_cased) {
	        SET_BIT(sample_uidx, pheno_cols_iter->data.cc);
		SET_BIT(sample_uidx, pheno_cols_iter->nonmiss);
	      } else if (dxx == pheno_ctrld) {
		SET_BIT(sample_uidx, pheno_cols_iter->nonmiss);
	      }
	    }
	  }
	  ++pheno_cols_iter;
	}
	pheno_info_reverse_ll = pheno_info_reverse_ll->next;
      }
    }
    *pheno_names_ptr = pheno_names;
    *max_pheno_name_blen_ptr = max_pheno_name_blen;
  }
  while (0) {
  load_phenos_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_phenos_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_phenos_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, pheno_fname);
  load_phenos_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
  load_phenos_ret_MALFORMED_INPUT_2:
    logerrprintb();
  load_phenos_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  load_phenos_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, pheno_fname);
    reterr = kPglRetMalformedInput;
    break;
  load_phenos_ret_INCOMPATIBLE_PHENOSTRS:
    LOGERRPRINTFWW("Error: Incompatible phenotype values in %s. (Case/control and quantitative phenotypes must be entirely numeric/\"NA\", and categorical phenotypes must be entirely non-numeric.)\n", pheno_fname);
    reterr = kPglRetMalformedInput;
    break;
  load_phenos_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  load_phenos_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 load_phenos_ret_1:
  bigstack_reset(bigstack_mark);
  gzclose_cond(gz_infile);
  if (reterr) {
    free_cond(pheno_names);
    if (*pheno_names_ptr) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
    }
    cleanup_pheno_cols(*pheno_ct_ptr, *pheno_cols_ptr);
    *pheno_ct_ptr = 0;
    *pheno_cols_ptr = nullptr;
  }
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
