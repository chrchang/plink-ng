// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

#include <ctype.h>
#include <time.h>
// no more mmap() dependency
// #include <fcntl.h>

#ifndef _WIN32
// #include <sys/mman.h>
#include <unistd.h>
#endif

// #include <sys/stat.h>
#include <sys/types.h>
#include "plink_family.h"
#include "plink_set.h"

#include "bgzf.h"
#include "pigz.h"

#define PHENO_EPSILON 0.000030517578125

int32_t sort_item_ids_nx(char** sorted_ids_ptr, uint32_t** id_map_ptr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len) {
  // Version of sort_item_ids() with no exclusion.
  // (Currently does NOT put id_map on the bottom.)
  uintptr_t ulii;
  char* sorted_ids;
  char* dup_id;
  char* tptr;
  if (bigstack_alloc_c(item_ct * max_id_len, sorted_ids_ptr) ||
      bigstack_alloc_ui(item_ct, id_map_ptr)) {
    return RET_NOMEM;
  }
  sorted_ids = *sorted_ids_ptr;
  if (!item_ct) {
    return 0;
  }
  for (ulii = 0; ulii < item_ct; ulii++) {
    memcpy(&(sorted_ids[ulii * max_id_len]), &(item_ids[ulii * max_id_len]), max_id_len);
    (*id_map_ptr)[ulii] = ulii;
  }
  if (qsort_ext(sorted_ids, item_ct, max_id_len, strcmp_deref, (char*)(*id_map_ptr), sizeof(int32_t))) {
    return RET_NOMEM;
  }
  dup_id = scan_for_duplicate_ids(sorted_ids, item_ct, max_id_len);
  if (dup_id) {
    tptr = strchr(dup_id, '\t');
    if (tptr) {
      *tptr = ' ';
    }
    LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", dup_id);
    return RET_INVALID_FORMAT;
  }
  return 0;
}

int32_t sample_major_to_snp_major(char* sample_major_fname, char* outname, uintptr_t unfiltered_marker_ct, uintptr_t unfiltered_sample_ct, uint64_t fsize) {
  // previously used mmap(); turns out this is more portable without being
  // noticeably slower.
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_marker_ct4 = (unfiltered_marker_ct + 3) / 4;
  uintptr_t unfiltered_marker_ctl2 = QUATERCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t marker_idx_end = 0;
  uint32_t bed_offset = fsize - unfiltered_sample_ct * ((uint64_t)unfiltered_marker_ct4);
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* lptr;
  unsigned char* writebuf;
  unsigned char* ucptr;
  uintptr_t cur_bigstack_left;
  uintptr_t write_marker_ct;
  uintptr_t marker_idx_base;
  uintptr_t marker_idx_block_end;
  uintptr_t marker_idx;
  uintptr_t sample_idx_base;
  uintptr_t sample_idx_end;
  uintptr_t sample_idx;
  uintptr_t cur_word0;
  uintptr_t cur_word1;
  uintptr_t cur_word2;
  uintptr_t cur_word3;
  if (fopen_checked(outname, FOPEN_WB, &outfile)) {
    goto sample_major_to_snp_major_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto sample_major_to_snp_major_ret_WRITE_FAIL;
  }
  if (unfiltered_marker_ct && unfiltered_sample_ct) {
    // could make this allocation a bit smaller in multipass case, but whatever
    if (bigstack_alloc_ul(unfiltered_marker_ctl2 * 4, &loadbuf)) {
      goto sample_major_to_snp_major_ret_NOMEM;
    }
    cur_bigstack_left = bigstack_left();
    if (cur_bigstack_left < unfiltered_sample_ct4) {
      goto sample_major_to_snp_major_ret_NOMEM;
    }
    writebuf = (unsigned char*)g_bigstack_base;
    write_marker_ct = BITCT2 * (cur_bigstack_left / (unfiltered_sample_ct4 * BITCT2));
    if (fopen_checked(sample_major_fname, FOPEN_RB, &infile)) {
      goto sample_major_to_snp_major_ret_OPEN_FAIL;
    }
    loadbuf[unfiltered_marker_ctl2 - 1] = 0;
    loadbuf[2 * unfiltered_marker_ctl2 - 1] = 0;
    loadbuf[3 * unfiltered_marker_ctl2 - 1] = 0;
    loadbuf[4 * unfiltered_marker_ctl2 - 1] = 0;
    do {
      marker_idx_base = marker_idx_end;
      marker_idx_end += write_marker_ct;
      if (marker_idx_end > unfiltered_marker_ct) {
	marker_idx_end = unfiltered_marker_ct;
      }
      if (fseeko(infile, bed_offset, SEEK_SET)) {
	goto sample_major_to_snp_major_ret_READ_FAIL;
      }
      for (sample_idx_end = 0; sample_idx_end < unfiltered_sample_ct;) {
	sample_idx_base = sample_idx_end;
	sample_idx_end = sample_idx_base + 4;
	if (sample_idx_end > unfiltered_sample_ct) {
	  fill_ulong_zero((4 - (unfiltered_sample_ct % 4)) * unfiltered_marker_ctl2, &(loadbuf[(unfiltered_sample_ct % 4) * unfiltered_marker_ctl2]));
	  sample_idx_end = unfiltered_sample_ct;
	}
	lptr = loadbuf;
	for (sample_idx = sample_idx_base; sample_idx < sample_idx_end; sample_idx++) {
	  if (load_raw(unfiltered_marker_ct4, infile, lptr)) {
	    goto sample_major_to_snp_major_ret_READ_FAIL;
	  }
	  lptr = &(lptr[unfiltered_marker_ctl2]);
	}
	lptr = &(loadbuf[marker_idx_base / BITCT2]);
	for (marker_idx_block_end = marker_idx_base; marker_idx_block_end < marker_idx_end; lptr++) {
	  marker_idx = marker_idx_block_end;
	  cur_word0 = *lptr;
	  cur_word1 = lptr[unfiltered_marker_ctl2];
	  cur_word2 = lptr[2 * unfiltered_marker_ctl2];
	  cur_word3 = lptr[3 * unfiltered_marker_ctl2];
	  marker_idx_block_end = marker_idx + BITCT2;
	  if (marker_idx_block_end > marker_idx_end) {
	    marker_idx_block_end = marker_idx_end;
	  }
	  ucptr = &(writebuf[(marker_idx - marker_idx_base) * unfiltered_sample_ct4 + (sample_idx_base / 4)]);
	  while (1) {
	    *ucptr = (unsigned char)((cur_word0 & 3) | ((cur_word1 & 3) << 2) | ((cur_word2 & 3) << 4) | ((cur_word3 & 3) << 6));
	    if (++marker_idx == marker_idx_block_end) {
	      break;
	    }
	    cur_word0 >>= 2;
	    cur_word1 >>= 2;
	    cur_word2 >>= 2;
	    cur_word3 >>= 2;
	    ucptr = &(ucptr[unfiltered_sample_ct4]);
	  }
	}
      }
      if (fwrite_checked(writebuf, (marker_idx_end - marker_idx_base) * unfiltered_sample_ct4, outfile)) {
	goto sample_major_to_snp_major_ret_WRITE_FAIL;
      }
    } while (marker_idx_end < unfiltered_marker_ct);
  }
  if (fclose_null(&outfile)) {
    goto sample_major_to_snp_major_ret_WRITE_FAIL;
  }

  while (0) {
  sample_major_to_snp_major_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sample_major_to_snp_major_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sample_major_to_snp_major_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  sample_major_to_snp_major_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(infile);
  fclose_cond(outfile);
  return retval;
}

int32_t load_map(FILE** mapfile_ptr, char* mapname, uint32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_blen_ptr, uintptr_t** marker_exclude_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, uint32_t** marker_pos_ptr, uint32_t* map_is_unsorted_ptr, uint32_t allow_extra_chroms, uint32_t allow_no_vars) {
  // currently only used by lgen_to_bed()
  // todo: some cleanup
  uintptr_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uintptr_t max_marker_id_blen = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t line_idx = 0;
  uint32_t last_pos = 0;
  int32_t last_chrom = -1;
  int32_t marker_pos_needed = 0;
  uint32_t chroms_encountered_m1 = 0xffffffffU; // deliberate overflow
  int32_t retval = 0;
  uintptr_t loaded_chrom_mask[CHROM_MASK_WORDS];
  uintptr_t* marker_exclude;
  uintptr_t unfiltered_marker_ctl;
  uintptr_t marker_uidx;
  uintptr_t ulii;
  uint32_t cur_pos;
  int32_t ii;
  {
    fill_ulong_zero(CHROM_MASK_WORDS, loaded_chrom_mask);
    if (fopen_checked(mapname, "r", mapfile_ptr)) {
      goto load_map_ret_OPEN_FAIL;
    }
    // first pass: count columns, determine raw marker count, determine maximum
    // marker ID length if necessary.
    g_textbuf[MAXLINELEN - 6] = ' ';
    while (fgets(g_textbuf, MAXLINELEN - 5, *mapfile_ptr)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 6]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .map file is pathologically long.\n", line_idx);
	goto load_map_ret_INVALID_FORMAT_2;
      }
      char* textbuf_first_token = skip_initial_spaces(g_textbuf);
      if (is_eoln_or_comment_kns(*textbuf_first_token)) {
	continue;
      }
      char* textbuf_iter = next_token(textbuf_first_token);
      if (no_more_tokens_kns(textbuf_iter)) {
	goto load_map_ret_MISSING_TOKENS;
      }
      ulii = strlen_se(textbuf_iter) + 1;
      if (ulii > max_marker_id_blen) {
	max_marker_id_blen = ulii;
      }
      if (!unfiltered_marker_ct) {
	// bugfix (24 Jul 2017): this was inappropriately erroring out on
	//   3-column .map files
	textbuf_iter = next_token(textbuf_iter);
	if (!textbuf_iter) {
	  goto load_map_ret_MISSING_TOKENS;
	}
	textbuf_iter = skip_initial_spaces(token_endnn(textbuf_iter));
	if (*textbuf_iter > ' ') {
	  *map_cols_ptr = 4;
	}
      }
      unfiltered_marker_ct++;
    }
    if (!feof(*mapfile_ptr)) {
      goto load_map_ret_READ_FAIL;
    }
    if ((!unfiltered_marker_ct) && (!allow_no_vars)) {
      logerrprint("Error: No variants in .map file.\n");
      goto load_map_ret_INVALID_FORMAT;
    }
    *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
    *max_marker_id_blen_ptr = max_marker_id_blen;
    rewind(*mapfile_ptr);
    unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);

    // unfiltered_marker_ct can be very large, so use bigstack for all
    // allocations that are a multiple of it

    // permanent bigstack allocation #1: marker_exclude
    if (bigstack_calloc_ul(unfiltered_marker_ctl, marker_exclude_ptr)) {
      goto load_map_ret_NOMEM;
    }
    marker_exclude = *marker_exclude_ptr;
    fill_uint_one(MAX_POSSIBLE_CHROM, chrom_info_ptr->chrom_idx_to_foidx);

    // permanent bigstack allocation #2, if needed: marker_pos
    if (marker_pos_needed) {
      if (bigstack_alloc_ui(unfiltered_marker_ct, marker_pos_ptr)) {
	goto load_map_ret_NOMEM;
      }
    }
    if (bigstack_alloc_c(unfiltered_marker_ct * max_marker_id_blen, marker_ids_ptr)) {
      goto load_map_ret_NOMEM;
    }

    // second pass: actually load stuff
    line_idx = 0;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      char* textbuf_first_token;
      if (get_next_noncomment(*mapfile_ptr, &textbuf_first_token, &line_idx)) {
	goto load_map_ret_READ_FAIL;
      }
      char* textbuf_iter = token_endnn(textbuf_first_token);
      if (!(*textbuf_iter)) {
	goto load_map_ret_MISSING_TOKENS;
      }
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code_destructive(".map file", line_idx, allow_extra_chroms, textbuf_first_token, textbuf_iter, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto load_map_ret_1;
      }
      if (cur_chrom_code != last_chrom) {
	if (cur_chrom_code < last_chrom) {
	  *map_is_unsorted_ptr |= UNSORTED_CHROM;
	}
	last_chrom = cur_chrom_code;
	if (is_set(loaded_chrom_mask, cur_chrom_code)) {
	  *map_is_unsorted_ptr |= UNSORTED_SPLIT_CHROM | UNSORTED_BP;
	} else {
	  set_bit(cur_chrom_code, loaded_chrom_mask);
	  chrom_info_ptr->chrom_file_order[++chroms_encountered_m1] = cur_chrom_code;
	  chrom_info_ptr->chrom_fo_vidx_start[chroms_encountered_m1] = marker_uidx;
	  chrom_info_ptr->chrom_idx_to_foidx[(uint32_t)cur_chrom_code] = chroms_encountered_m1;
	}
	last_pos = 0;
      }

      if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	SET_BIT(marker_uidx, marker_exclude);
	marker_exclude_ct++;
      } else {
	textbuf_iter = skip_initial_spaces(&(textbuf_iter[1]));
	if (is_eoln_kns(*textbuf_iter)) {
	  goto load_map_ret_MISSING_TOKENS;
	}
	char* token_end = token_endnn(textbuf_iter);
	memcpyx(&((*marker_ids_ptr)[marker_uidx * max_marker_id_blen]), textbuf_iter, (uintptr_t)(token_end - textbuf_iter), '\0');
	textbuf_iter = next_token_mult(token_end, *map_cols_ptr - 2);
	if (no_more_tokens_kns(textbuf_iter)) {
	  goto load_map_ret_MISSING_TOKENS;
	}
	if (scan_int_abs_defcap(textbuf_iter, &ii)) {
	  sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of .map file.\n", line_idx);
	  goto load_map_ret_INVALID_FORMAT_2;
	}
	if (ii < 0) {
	  SET_BIT(marker_uidx, marker_exclude);
	  marker_exclude_ct++;
	} else {
	  cur_pos = ii;
	  if (cur_pos < last_pos) {
	    *map_is_unsorted_ptr |= UNSORTED_BP;
	  } else {
	    last_pos = cur_pos;
	  }
	  if (marker_pos_needed && cur_chrom_code) {
	    (*marker_pos_ptr)[marker_uidx] = cur_pos;
	  }
	}
      }
    }
    chrom_info_ptr->chrom_ct = ++chroms_encountered_m1;
    *marker_exclude_ct_ptr = marker_exclude_ct;
    if (unfiltered_marker_ct) {
      if (marker_exclude_ct == unfiltered_marker_ct) {
	logerrprint("Error: All variants excluded from .map file.\n");
	goto load_map_ret_ALL_MARKERS_EXCLUDED;
      }
    }
    chrom_info_ptr->chrom_fo_vidx_start[chroms_encountered_m1] = marker_uidx;
  }
  while (0) {
  load_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_map_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .map file has fewer tokens than expected.\n", line_idx);
  load_map_ret_INVALID_FORMAT_2:
    logerrprintb();
  load_map_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  load_map_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 load_map_ret_1:
  return retval;
}

// missing code set to 1 by load_bim()
static uint8_t acgtm_bool_table[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static inline uint32_t is_acgtm(unsigned char ucc) {
  return (uint32_t)(acgtm_bool_table[ucc]);
}

int32_t load_bim(char* bimname, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_blen_ptr, uintptr_t** marker_exclude_ptr, double** set_allele_freqs_ptr, uint32_t** nchrobs_ptr, char*** marker_allele_pp, uintptr_t* max_marker_allele_blen_ptr, char** marker_ids_ptr, char* missing_mid_template, uint32_t new_id_max_allele_slen, const char* missing_marker_id_match, Chrom_info* chrom_info_ptr, double** marker_cms_ptr, uint32_t** marker_pos_ptr, uint64_t misc_flags, uint64_t filter_flags, int32_t marker_pos_start, int32_t marker_pos_end, int32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, uint32_t* map_is_unsorted_ptr, uint32_t marker_pos_needed, uint32_t marker_cms_needed, uint32_t marker_alleles_needed, const char* split_chrom_cmd, const char* ftype_str, uint32_t* max_bim_linelen_ptr) {
  // supports .map now too, to make e.g. --snps + --dosage work
  FILE* bimfile = nullptr;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uintptr_t max_marker_id_blen = *max_marker_id_blen_ptr;
  uintptr_t max_marker_allele_blen = *max_marker_allele_blen_ptr;
  uintptr_t line_idx = 0;
  int32_t prev_chrom = -1;
  uint32_t last_pos = 0;
  double last_cm = -DBL_MAX;
  const uint32_t is_bim = (ftype_str[1] == 'b'); // .map also supported
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t exclude_snp = (filter_flags / FILTER_EXCLUDE_MARKERNAME_SNP) & 1;
  uint32_t snps_only = (filter_flags / FILTER_SNPS_ONLY) & 1;
  uint32_t snps_only_just_acgt = (misc_flags / MISC_SNPS_ONLY_JUST_ACGT) & 1;
  uint32_t from_slen = markername_from? strlen(markername_from) : 0;
  uint32_t to_slen = markername_to? strlen(markername_to) : 0;
  uint32_t snp_slen = markername_snp? strlen(markername_snp) : 0;
  uint32_t slen_check = from_slen || to_slen || snp_slen;
  uint32_t from_chrom = MAX_POSSIBLE_CHROM;
  uint32_t to_chrom = MAX_POSSIBLE_CHROM;
  uint32_t snp_chrom = MAX_POSSIBLE_CHROM;
  uint32_t chroms_encountered_m1 = 0xffffffffU; // intentional overflow
  uint32_t snp_pos = 0;
  uint32_t mcm2 = 1;
  uint32_t split_chrom = 0;
  uint32_t max_bim_linelen = 0;
  uint32_t missing_marker_id_match_len = 0;
  uint32_t missing_ids_set = 0;
  uint32_t missing_template_base_len = 0;
  uint32_t template_insert_ct = 0;
  uint32_t chrom_header_line_present = 0;
  uint32_t uii = 0;
  uint32_t ujj = 0;
  int32_t exclude_window_start = 0;
  int32_t exclude_window_end = -1;
  int32_t retval = 0;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char* loadbuf2 = nullptr; // on heap, second pass
  char* prev_new_id = nullptr;
  char* bufptr2 = nullptr;
  char* bufptr4 = nullptr;
  char* bufptr5 = nullptr;
  char** marker_allele_ptrs = nullptr;
  uintptr_t loaded_chrom_mask[CHROM_MASK_WORDS];
  uint32_t missing_template_seg_len[5];
  uint32_t missing_template_seg_order[4]; // '@', '#', '$1', '$2'
  uint32_t insert_buf_len[4];
  char poscharbuf[12];
  char* missing_template_seg[5];
  char* insert_buf[4];
  uintptr_t* marker_exclude;
  char* loadbuf; // on stack, first pass
  uintptr_t loadbuf_size;
  uintptr_t unfiltered_marker_ctl;
  uintptr_t marker_uidx;
  char* bufptr;
  char* col2_ptr;
  uintptr_t ulii;
  uint32_t ukk;
  uint32_t umm;
  int32_t jj;
  uint32_t cur_pos;
  double cur_cm;
  char cc;
  {
    fill_ulong_zero(CHROM_MASK_WORDS, loaded_chrom_mask);
    insert_buf[0] = nullptr;
    insert_buf[1] = nullptr;
    insert_buf[2] = nullptr;
    insert_buf[3] = nullptr;
    fill_uint_zero(5, missing_template_seg_len);
    missing_template_seg[0] = nullptr;
    missing_template_seg[1] = nullptr;
    missing_template_seg[2] = nullptr;
    missing_template_seg[3] = nullptr;
    missing_template_seg[4] = nullptr;
    if (missing_mid_template) {
      if (!missing_marker_id_match) {
	missing_marker_id_match = &(g_one_char_strs[92]); // '.'
      }
      missing_marker_id_match_len = strlen(missing_marker_id_match);
      bufptr = missing_mid_template;
      missing_template_seg[template_insert_ct] = bufptr; // current segment start
      cc = *bufptr; // template string previously validated
      do {
	if (cc == '@') {
	  ujj = (uintptr_t)(bufptr - missing_template_seg[template_insert_ct]);
	  ukk = 0;
	  goto load_bim_template_match;
	} else if (cc == '#') {
	  ujj = (uintptr_t)(bufptr - missing_template_seg[template_insert_ct]);
	  ukk = 1;
	  goto load_bim_template_match;
	} else if (cc == '$') {
	  ujj = (uintptr_t)(bufptr - missing_template_seg[template_insert_ct]);
	  cc = *(++bufptr);
	  ukk = ((unsigned char)cc) - 47;
	load_bim_template_match:
	  missing_template_seg_len[template_insert_ct] = ujj;
	  missing_template_base_len += ujj;
	  missing_template_seg_order[template_insert_ct++] = ukk;
	  missing_template_seg[template_insert_ct] = &(bufptr[1]);
	}
	cc = *(++bufptr);
      } while (cc);
      ujj = (uintptr_t)(bufptr - missing_template_seg[template_insert_ct]);
      missing_template_seg_len[template_insert_ct] = ujj;
      missing_template_base_len += ujj;
      insert_buf[1] = poscharbuf;
      if (template_insert_ct == 4) {
	insert_buf[2] = (char*)malloc(new_id_max_allele_slen + 1);
	insert_buf[3] = (char*)malloc(new_id_max_allele_slen + 1);
	if ((!insert_buf[2]) || (!insert_buf[3])) {
	  goto load_bim_ret_NOMEM;
	}
      }
    }
    if (fopen_checked(bimname, "r", &bimfile)) {
      goto load_bim_ret_OPEN_FAIL;
    }
    // first pass: count columns, determine raw marker count, determine maximum
    // marker ID length and/or marker allele length if necessary, save
    // nonstandard chromosome names.

    // ensure strcmp_se comparison doesn't read past end of buffer
    loadbuf_size = bigstack_left() - 16;
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto load_bim_ret_NOMEM;
    }
    loadbuf = (char*)g_bigstack_base;
    loadbuf[loadbuf_size - 1] = ' ';
    while (fgets(loadbuf, loadbuf_size, bimfile)) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, ftype_str);
	  goto load_bim_ret_INVALID_FORMAT_2;
	} else {
	  goto load_bim_ret_NOMEM;
	}
      }
      uii = strlen(loadbuf);
      if (uii >= max_bim_linelen) {
	max_bim_linelen = uii + 1;
      }
      char* loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (is_eoln_or_comment_kns(*loadbuf_first_token)) {
	if (is_bim) {
	  if (!strcmp_se(loadbuf_first_token, "#CHROM", 6)) {
	    if (chrom_header_line_present) {
	      sprintf(g_logbuf, "Error: Multiple #CHROM header lines in %s.\n", ftype_str);
	      goto load_bim_ret_INVALID_FORMAT_2;
	    }
	    // support plink 2.x files with default column order; error out on
	    // VCF and other column orders
	    loadbuf_first_token = skip_initial_spaces(&(loadbuf_first_token[6]));
	    if (strcmp_se(loadbuf_first_token, "ID", 2)) {
	      goto load_bim_ret_UNSUPPORTED_COLUMN_ORDER;
	    }
	    loadbuf_first_token = skip_initial_spaces(&(loadbuf_first_token[2]));
	    if (!strcmp_se(loadbuf_first_token, "CM", 2)) {
	      mcm2 = 2;
	      loadbuf_first_token = skip_initial_spaces(&(loadbuf_first_token[2]));
	    }
	    if (strcmp_se(loadbuf_first_token, "POS", 3)) {
	      goto load_bim_ret_UNSUPPORTED_COLUMN_ORDER;
	    }
	    loadbuf_first_token = skip_initial_spaces(&(loadbuf_first_token[3]));
	    if (strcmp_se(loadbuf_first_token, "ALT", 3)) {
	      goto load_bim_ret_UNSUPPORTED_COLUMN_ORDER;
	    }
	    loadbuf_first_token = skip_initial_spaces(&(loadbuf_first_token[3]));
	    if (strcmp_se(loadbuf_first_token, "REF", 3)) {
	      goto load_bim_ret_UNSUPPORTED_COLUMN_ORDER;
	    }
	    chrom_header_line_present = 1;
	  }
	} else if (*loadbuf_first_token == '#') {
	  sprintf(g_logbuf, "Error: Header lines are not permitted in %ss.\n", ftype_str);
	  goto load_bim_ret_INVALID_FORMAT_2;
	}
	continue;
      }
      char* first_token_end = token_endnn(loadbuf_first_token);
      col2_ptr = skip_initial_spaces(first_token_end);
      if (is_eoln_kns(*col2_ptr)) {
	goto load_bim_ret_MISSING_TOKENS;
      }
      const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - loadbuf_first_token);
      *first_token_end = '\0';
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code(loadbuf_first_token, ftype_str, line_idx, chrom_name_slen, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto load_bim_ret_1;
      }

      char* col2_end = token_endnn(col2_ptr);
      ulii = (uintptr_t)(col2_end - col2_ptr);
      if (!unfiltered_marker_ct) {
	if (is_bim) {
	  // .bim: bufptr2 = col 5 start
	  bufptr2 = next_token_mult(col2_end, 3);
	} else {
	  // .map
	  bufptr2 = skip_initial_spaces(col2_end);
	}
	if (no_more_tokens_kns(bufptr2)) {
	  goto load_bim_ret_MISSING_TOKENS;
	}
	// check if CM col exists
	if ((!chrom_header_line_present) && (!is_eoln_kns(*(skip_initial_spaces(token_endnn(bufptr2)))))) {
	  mcm2 = 2;
	}
      }
      if (marker_alleles_needed || missing_marker_id_match_len) {
	bufptr4 = next_token_mult(col2_end, mcm2 + 1);
	bufptr5 = next_token(bufptr4);
	if (no_more_tokens_kns(bufptr5)) {
	  goto load_bim_ret_MISSING_TOKENS;
	}
	uii = strlen_se(bufptr4);
	ujj = strlen_se(bufptr5);
	if (memchr(bufptr4, ',', ujj + ((uintptr_t)(bufptr5 - bufptr4)))) {
	  // this breaks VCF and plink 2.x
	  // may need to add word wrapping if this message is changed
	  sprintf(g_logbuf, "Error: Comma-containing allele code on line %" PRIuPTR " of %s.\n", line_idx, ftype_str);
	  goto load_bim_ret_INVALID_FORMAT_2;
	}
	if (marker_alleles_needed) {
	  if (uii >= max_marker_allele_blen) {
	    max_marker_allele_blen = uii + 1;
	  }
	  if (ujj >= max_marker_allele_blen) {
	    max_marker_allele_blen = ujj + 1;
	  }
	}
      }
      if ((ulii == missing_marker_id_match_len) && (!memcmp(col2_ptr, missing_marker_id_match, missing_marker_id_match_len))) {
	bufptr2 = next_token_mult(col2_end, mcm2);
	if (no_more_tokens_kns(bufptr2)) {
	  goto load_bim_ret_MISSING_TOKENS;
	}
	insert_buf_len[1] = strlen_se(bufptr2);
	if (insert_buf_len[1] > 11) {
	  // permit negative sign and 10 digit number
	  goto load_bim_ret_INVALID_BP_COORDINATE;
	}
	insert_buf_len[0] = chrom_name_slen;
	ulii = missing_template_base_len + insert_buf_len[1] + insert_buf_len[0];
	if (template_insert_ct == 4) {
	  uii = MINV(uii, new_id_max_allele_slen);
	  ujj = MINV(ujj, new_id_max_allele_slen);
	  ulii += uii + ujj;
	}
	if (ulii >= max_marker_id_blen) {
	  if (ulii > MAX_ID_BLEN) {
	    logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  max_marker_id_blen = ulii + 1;
	}
	ulii = 0;
      } else {
	if (ulii >= max_marker_id_blen) {
	  max_marker_id_blen = ulii + 1;
	}
      }
      if (slen_check) {
	if (!ulii) {
	  // --set-missing-var-ids applies
	  // safe to clobber buffer contents
	  insert_buf[0] = loadbuf_first_token;
	  memcpyx(poscharbuf, bufptr2, insert_buf_len[1], '\0');
	  if (template_insert_ct == 4) {
	    bufptr4[uii] = '\0';
	    bufptr5[ujj] = '\0';
	    // ASCII-sort allele names
	    if (strcmp(bufptr4, bufptr5) <= 0) {
	      memcpy(insert_buf[2], bufptr4, uii);
	      insert_buf_len[2] = uii;
	      memcpy(insert_buf[3], bufptr5, ujj);
	      insert_buf_len[3] = ujj;
	    } else {
	      memcpy(insert_buf[3], bufptr4, uii);
	      insert_buf_len[3] = uii;
	      memcpy(insert_buf[2], bufptr5, ujj);
	      insert_buf_len[2] = ujj;
	    }
	  }
	  bufptr4 = col2_ptr;
	  for (uii = 0; uii < template_insert_ct; uii++) {
	    bufptr4 = memcpya(bufptr4, missing_template_seg[uii], missing_template_seg_len[uii]);
	    ujj = missing_template_seg_order[uii];
	    bufptr4 = memcpya(bufptr4, insert_buf[ujj], insert_buf_len[ujj]);
	  }
	  bufptr4 = memcpya(bufptr4, missing_template_seg[uii], missing_template_seg_len[uii]);
	  ulii = (uintptr_t)(bufptr4 - col2_ptr);
	  bufptr2 = poscharbuf;
	} else {
	  bufptr2 = next_token_mult(col2_ptr, mcm2);
	  if (no_more_tokens_kns(bufptr2)) {
	    goto load_bim_ret_MISSING_TOKENS;
	  }
	  bufptr2[strlen_se(bufptr2)] = '\0';
	}
        if ((ulii == from_slen) && (!memcmp(col2_ptr, markername_from, ulii))) {
          if (from_chrom != MAX_POSSIBLE_CHROM) {
            goto load_bim_ret_DUPLICATE_ID;
          }
          from_chrom = cur_chrom_code;
          if (scan_uint_defcap(bufptr2, (uint32_t*)&marker_pos_start)) {
            goto load_bim_ret_INVALID_BP_COORDINATE;
          }
          if (to_chrom != MAX_POSSIBLE_CHROM) {
            if (from_chrom != to_chrom) {
              goto load_bim_ret_FROM_TO_DIFFERENT_CHROM;
            }
          }
          fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->chrom_mask);
          SET_BIT(from_chrom, chrom_info_ptr->chrom_mask);
        }
        if ((ulii == to_slen) && (!memcmp(col2_ptr, markername_to, ulii))) {
          if (to_chrom != MAX_POSSIBLE_CHROM) {
            goto load_bim_ret_DUPLICATE_ID;
          }
          to_chrom = cur_chrom_code;
          if (scan_uint_defcap(bufptr2, (uint32_t*)&marker_pos_end)) {
            goto load_bim_ret_INVALID_BP_COORDINATE;
          }
          if (from_chrom != MAX_POSSIBLE_CHROM) {
            if (to_chrom != from_chrom) {
              goto load_bim_ret_FROM_TO_DIFFERENT_CHROM;
            }
          }
          fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->chrom_mask);
          SET_BIT(to_chrom, chrom_info_ptr->chrom_mask);
        }
        if ((ulii == snp_slen) && (!memcmp(col2_ptr, markername_snp, ulii))) {
          if (snp_chrom != MAX_POSSIBLE_CHROM) {
            goto load_bim_ret_DUPLICATE_ID;
          }
          snp_chrom = cur_chrom_code;
          if (scan_uint_defcap(bufptr2, &snp_pos)) {
            goto load_bim_ret_INVALID_BP_COORDINATE;
          }
          if (!exclude_snp) {
            fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->chrom_mask);
            SET_BIT(snp_chrom, chrom_info_ptr->chrom_mask);
          }
        }
      }
      unfiltered_marker_ct++;
    }
    if (!feof(bimfile)) {
      goto load_bim_ret_READ_FAIL;
    }
    if ((!unfiltered_marker_ct) && (!allow_no_variants)) {
      sprintf(g_logbuf, "Error: No variants in %s.\n", ftype_str);
      goto load_bim_ret_INVALID_FORMAT_2;
    } else if (unfiltered_marker_ct > 2147483645) {
      // maximum prime < 2^32 is 4294967291; quadratic hashing guarantee breaks
      // down past that divided by 2.
      // PLINK/SEQ now supports a 64-bit count here, and few other tools do, so
      // it's appropriate to explicitly recommend it.
      logerrprint("Error: PLINK does not support more than 2^31 - 3 variants.  We recommend other\nsoftware, such as PLINK/SEQ, for very deep studies of small numbers of genomes.\n");
      goto load_bim_ret_INVALID_FORMAT;
    }
    if (from_slen || to_slen) {
      if (from_slen && (from_chrom == MAX_POSSIBLE_CHROM)) {
	LOGPREPRINTFWW("Error: --from variant '%s' not found.\n", markername_from);
	goto load_bim_ret_INVALID_FORMAT_2;
      }
      if (to_slen && (to_chrom == MAX_POSSIBLE_CHROM)) {
	LOGPREPRINTFWW("Error: --to variant '%s' not found.\n", markername_to);
	goto load_bim_ret_INVALID_FORMAT_2;
      }
      if (marker_pos_start == -1) {
	marker_pos_start = 0;
      }
      if (marker_pos_end == -1) {
	marker_pos_end = 0x7fffffff;
      }
      if (marker_pos_start > marker_pos_end) {
	jj = marker_pos_start;
	marker_pos_start = marker_pos_end;
	marker_pos_end = jj;
      }
    }
    if (snp_slen) {
      if (snp_chrom == MAX_POSSIBLE_CHROM) {
	LOGPREPRINTFWW("Error: --%ssnp variant '%s' not found.\n", exclude_snp? "exclude-" : "", markername_snp);
	goto load_bim_ret_INVALID_FORMAT_2;
      }
      if (!exclude_snp) {
	if (snp_window_size == -1) {
	  // no harm in screening on position before variant ID
	  uii = 0;
	} else {
	  uii = snp_window_size;
	}
	if (uii > snp_pos) {
	  marker_pos_start = 0;
	} else {
	  marker_pos_start = snp_pos - uii;
	}
	if (uii > (0x7fffffff - snp_pos)) {
	  marker_pos_end = 0x7fffffff;
	} else {
	  marker_pos_end = snp_pos + uii;
	}
      } else if (snp_window_size != -1) {
	if ((uint32_t)snp_window_size <= snp_pos) {
	  exclude_window_start = snp_pos - snp_window_size;
	}
	if ((uint32_t)snp_window_size > (0x7fffffff - snp_pos)) {
	  exclude_window_end = 0x7fffffff;
	} else {
	  exclude_window_end = snp_pos + snp_window_size;
	}
      }
    }

    if (max_marker_id_blen > MAX_ID_BLEN) {
      logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto load_bim_ret_INVALID_FORMAT;
    }
    *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
    *max_marker_id_blen_ptr = max_marker_id_blen;
    rewind(bimfile);
    unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);

    // unfiltered_marker_ct can be very large, so use bigstack for all
    // allocations that are a multiple of it

    // permanent bigstack allocation #1: marker_exclude
    // permanent bigstack allocation #2: set_allele_freqs
    if (bigstack_calloc_ul(unfiltered_marker_ctl, marker_exclude_ptr)) {
      goto load_bim_ret_NOMEM;
    }
    marker_exclude = *marker_exclude_ptr;
    if (set_allele_freqs_ptr) {
      if (bigstack_alloc_d(unfiltered_marker_ct, set_allele_freqs_ptr)) {
	goto load_bim_ret_NOMEM;
      }
      // leave set_allele_freqs uninitialized
      if (nchrobs_ptr) {
	if (bigstack_alloc_ui(unfiltered_marker_ct, nchrobs_ptr)) {
	  goto load_bim_ret_NOMEM;
	}
	// on the other hand, this is not autocomputed
	fill_uint_one(unfiltered_marker_ct, *nchrobs_ptr);
      }
    }
    fill_uint_one(MAX_POSSIBLE_CHROM, chrom_info_ptr->chrom_idx_to_foidx);
    // permanent bigstack allocation #3, if needed: marker_pos
    if (marker_pos_needed) {
      if (bigstack_alloc_ui(unfiltered_marker_ct, marker_pos_ptr)) {
	goto load_bim_ret_NOMEM;
      }
    }
    if (marker_alleles_needed) {
      if (snps_only) {
	max_marker_allele_blen = 2;
	acgtm_bool_table[(unsigned char)(*missing_geno_ptr)] = 1;
      }
      if (max_marker_allele_blen > NON_BIGSTACK_MIN - 1) {
	// guard against overflows
	LOGERRPRINTF("Error: Alleles are limited to %u characters.\n", NON_BIGSTACK_MIN - 1);
	goto load_bim_ret_INVALID_FORMAT;
      }
      *max_marker_allele_blen_ptr = max_marker_allele_blen;
      marker_allele_ptrs = (char**)bigstack_alloc(unfiltered_marker_ct * 2 * sizeof(intptr_t));
      if (!marker_allele_ptrs) {
	goto load_bim_ret_NOMEM;
      }
      *marker_allele_pp = marker_allele_ptrs;
      ujj = unfiltered_marker_ct * 2;
      for (uii = 0; uii < ujj; uii++) {
	marker_allele_ptrs[uii] = missing_geno_ptr;
      }
    }
    if (bigstack_alloc_c(unfiltered_marker_ct * max_marker_id_blen, marker_ids_ptr)) {
      goto load_bim_ret_NOMEM;
    }
    // todo: check whether marker_cms can be unloaded before
    // marker_ids/marker_alleles, or vice versa
    if (marker_cms_needed & MARKER_CMS_FORCED) {
      if (bigstack_calloc_d(unfiltered_marker_ct, marker_cms_ptr)) {
	goto load_bim_ret_NOMEM;
      }
    }
    if ((filter_flags & FILTER_ZERO_CMS) || (mcm2 == 1)) {
      marker_cms_needed = 0;
    }

    // second pass: actually load stuff
    loadbuf2 = (char*)malloc(max_bim_linelen);
    if (!loadbuf2) {
      goto load_bim_ret_NOMEM;
    }
    if (missing_mid_template) {
      prev_new_id = (char*)malloc(max_marker_id_blen);
      if (!prev_new_id) {
	goto load_bim_ret_NOMEM;
      }
      *prev_new_id = '\0';
    }
    line_idx = 0;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      char* loadbuf_first_token;
      do {
	line_idx++;
	if (!fgets(loadbuf2, max_bim_linelen, bimfile)) {
	  goto load_bim_ret_READ_FAIL;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf2);
      } while (is_eoln_or_comment_kns(*loadbuf_first_token));
      char* first_token_end = token_endnn(loadbuf_first_token);
      col2_ptr = skip_initial_spaces(first_token_end);
      const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - loadbuf_first_token);
      *first_token_end = '\0';
      int32_t cur_chrom_code = get_chrom_code(loadbuf_first_token, chrom_info_ptr, chrom_name_slen);
      if (cur_chrom_code != prev_chrom) {
	if (!split_chrom) {
	  if (cur_chrom_code < prev_chrom) {
	    *map_is_unsorted_ptr |= UNSORTED_CHROM;
	  }
	  prev_chrom = cur_chrom_code;
	  if (is_set(loaded_chrom_mask, cur_chrom_code)) {
	    if (split_chrom_cmd) {
	      sprintf(g_logbuf, "Error: %s has a split chromosome.  Use --%s by itself to\nremedy this.\n", ftype_str, split_chrom_cmd);
	      goto load_bim_ret_INVALID_FORMAT_2;
	    }
	    split_chrom = 1;
	    *map_is_unsorted_ptr |= UNSORTED_CHROM | UNSORTED_BP | UNSORTED_SPLIT_CHROM;
	  } else {
	    chrom_info_ptr->chrom_file_order[++chroms_encountered_m1] = cur_chrom_code;
	    chrom_info_ptr->chrom_fo_vidx_start[chroms_encountered_m1] = marker_uidx;
	    chrom_info_ptr->chrom_idx_to_foidx[(uint32_t)cur_chrom_code] = chroms_encountered_m1;
	  }
	  last_pos = 0;
	  last_cm = -DBL_MAX;
	}
	set_bit(cur_chrom_code, loaded_chrom_mask);
      }

      if (is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	if (no_more_tokens_kns(col2_ptr)) {
	  goto load_bim_ret_MISSING_TOKENS;
	}
	uii = strlen_se(col2_ptr);
	ujj = (uii == missing_marker_id_match_len) && (!memcmp(col2_ptr, missing_marker_id_match, missing_marker_id_match_len));
	if (!ujj) {
	  memcpyx(&((*marker_ids_ptr)[marker_uidx * max_marker_id_blen]), col2_ptr, uii, '\0');
	}
	if (marker_cms_needed) {
	  bufptr = next_token(col2_ptr);
	  if (no_more_tokens_kns(bufptr)) {
	    goto load_bim_ret_MISSING_TOKENS;
	  }
	  if ((*bufptr != '0') || (bufptr[1] > ' ')) {
	    if (!(*marker_cms_ptr)) {
	      if (bigstack_calloc_d(unfiltered_marker_ct, marker_cms_ptr)) {
		goto load_bim_ret_NOMEM;
	      }
	    }
	    if (scan_double(bufptr, &cur_cm)) {
	      sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, ftype_str);
	      goto load_bim_ret_INVALID_FORMAT_2;
	    }
	    if (cur_cm < last_cm) {
	      *map_is_unsorted_ptr |= UNSORTED_CM;
	    } else {
	      last_cm = cur_cm;
	    }
	    (*marker_cms_ptr)[marker_uidx] = cur_cm;
	  }
	  bufptr = next_token(bufptr);
	} else {
	  bufptr = next_token_mult(col2_ptr, mcm2);
	}
	if (no_more_tokens_kns(bufptr)) {
	  goto load_bim_ret_MISSING_TOKENS;
	}
	if (scan_int_abs_defcap(bufptr, (int32_t*)&cur_pos)) {
	  goto load_bim_ret_INVALID_BP_COORDINATE;
	}
	// negative marker positions now have the same effect in .bim as .map
	if ((int32_t)cur_pos < 0) {
	  goto load_bim_skip_marker;
	}
	if (cur_pos < last_pos) {
	  *map_is_unsorted_ptr |= UNSORTED_BP;
	} else {
	  last_pos = cur_pos;
	}
	if ((marker_pos_start != -1) && ((((int32_t)cur_pos) < marker_pos_start) || (((int32_t)cur_pos) > marker_pos_end))) {
	  goto load_bim_skip_marker;
	}
	if (snp_slen) {
	  if (snp_window_size == -1) {
	    if ((uii == snp_slen) && (!memcmp(col2_ptr, markername_snp, snp_slen))) {
	      if (exclude_snp) {
		goto load_bim_skip_marker;
	      }
	    } else if (!exclude_snp) {
	      goto load_bim_skip_marker;
	    }
	  } else if (exclude_snp && ((((int32_t)cur_pos) <= exclude_window_end) && (((int32_t)cur_pos) >= exclude_window_start) && ((uint32_t)cur_chrom_code == snp_chrom))) {
	    goto load_bim_skip_marker;
	  }
	}
	if (marker_pos_needed) {
	  (*marker_pos_ptr)[marker_uidx] = cur_pos;
	}
	if (marker_alleles_needed || ujj) {
	  bufptr4 = next_token(bufptr);
	  bufptr5 = next_token(bufptr4);
	  if (!bufptr5) {
	    goto load_bim_ret_MISSING_TOKENS;
	  }
	  ukk = strlen_se(bufptr4);
	  umm = strlen_se(bufptr5);
	  if (marker_alleles_needed) {
	    if (snps_only) {
	      if ((ukk != 1) || (umm != 1) || (snps_only_just_acgt && ((!is_acgtm(*bufptr4)) || (!is_acgtm(*bufptr5))))) {
		goto load_bim_skip_marker;
	      }
	    }
	    ulii = marker_uidx * 2;
	    if (allele_set(bufptr4, ukk, &(marker_allele_ptrs[ulii]))) {
	      goto load_bim_ret_NOMEM;
	    }
	    ulii++;
	    if (allele_set(bufptr5, umm, &(marker_allele_ptrs[ulii]))) {
	      goto load_bim_ret_NOMEM;
	    }
	  }
	  if (ujj) {
	    // --set-missing-var-ids
	    // bufptr = position string
	    // loadbuf_first_token = chromosome code
	    // bufptr4 and bufptr5: alleles (ok to null-terminate)
	    // ukk and umm: allele lengths
	    insert_buf[0] = loadbuf_first_token;
	    insert_buf_len[0] = chrom_name_slen;
	    insert_buf[1] = bufptr;
	    insert_buf_len[1] = strlen_se(bufptr);
	    if (template_insert_ct == 4) {
	      ukk = MINV(ukk, new_id_max_allele_slen);
	      umm = MINV(umm, new_id_max_allele_slen);
	      bufptr4[ukk] = '\0';
	      bufptr5[umm] = '\0';
	      // ASCII-sort allele names
	      if (strcmp(bufptr4, bufptr5) <= 0) {
		memcpy(insert_buf[2], bufptr4, ukk);
		insert_buf_len[2] = ukk;
		memcpy(insert_buf[3], bufptr5, umm);
		insert_buf_len[3] = umm;
	      } else {
		memcpy(insert_buf[3], bufptr4, ukk);
		insert_buf_len[3] = ukk;
		memcpy(insert_buf[2], bufptr5, umm);
		insert_buf_len[2] = umm;
	      }
	    }
	    bufptr5 = &((*marker_ids_ptr)[marker_uidx * max_marker_id_blen]);
	    bufptr4 = bufptr5;
	    for (uii = 0; uii < template_insert_ct; uii++) {
	      bufptr4 = memcpya(bufptr4, missing_template_seg[uii], missing_template_seg_len[uii]);
	      ujj = missing_template_seg_order[uii];
	      bufptr4 = memcpya(bufptr4, insert_buf[ujj], insert_buf_len[ujj]);
	    }
	    bufptr4 = memcpyax(bufptr4, missing_template_seg[uii], missing_template_seg_len[uii], '\0');
	    if (!strcmp(prev_new_id, bufptr5)) {
	      LOGERRPRINTFWW("Error: Duplicate ID '%s' generated by --set-missing-var-ids.\n", prev_new_id);
	      goto load_bim_ret_INVALID_CMDLINE;
	    }
	    missing_ids_set++;
	    memcpy(prev_new_id, bufptr5, (uintptr_t)(bufptr4 - bufptr5));
	  }
	}
      } else {
      load_bim_skip_marker:
	SET_BIT(marker_uidx, marker_exclude);
	marker_exclude_ct++;
	if (marker_pos_needed) {
	  // support unfiltered marker_pos search
	  (*marker_pos_ptr)[marker_uidx] = last_pos;
	}
	if (marker_cms_needed & MARKER_CMS_FORCED) {
	  (*marker_cms_ptr)[marker_uidx] = last_cm;
	}
      }
    }
    if ((unfiltered_marker_ct == marker_exclude_ct) && (!allow_no_variants)) {
      logerrprint("Error: All variants excluded.\n");
      goto load_bim_ret_ALL_MARKERS_EXCLUDED;
    }
    if (missing_mid_template && ((*map_is_unsorted_ptr) & UNSORTED_BP)) {
      sprintf(g_logbuf, "Error: --set-missing-var-ids requires a sorted %s.  Retry this command\nafter using --make-bed to sort your data.\n", ftype_str);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    for (uii = 0; uii < CHROM_MASK_WORDS; uii++) {
      chrom_info_ptr->chrom_mask[uii] &= loaded_chrom_mask[uii];
    }
    chrom_info_ptr->chrom_ct = ++chroms_encountered_m1;
    chrom_info_ptr->chrom_fo_vidx_start[chroms_encountered_m1] = marker_uidx;
    *marker_exclude_ct_ptr = marker_exclude_ct;
    if (!marker_exclude_ct) {
      LOGPRINTF("%" PRIuPTR " variant%s loaded from %s.\n", unfiltered_marker_ct, (unfiltered_marker_ct == 1)? "" : "s", ftype_str);
    } else {
      LOGPRINTF("%" PRIuPTR " out of %" PRIuPTR " variant%s loaded from %s.\n", unfiltered_marker_ct - marker_exclude_ct, unfiltered_marker_ct, (unfiltered_marker_ct == 1)? "" : "s", ftype_str);
    }
    if (missing_ids_set) {
      LOGPRINTF("%u missing ID%s set.\n", missing_ids_set, (missing_ids_set == 1)? "" : "s");
    }

    if (max_bim_linelen_ptr) {
      *max_bim_linelen_ptr = max_bim_linelen;
    }
  }
  while (0) {
  load_bim_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_bim_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_bim_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  load_bim_ret_UNSUPPORTED_COLUMN_ORDER:
    LOGERRPRINTF("Error: Unsupported column order specified on line %" PRIuPTR " of %s.\n", line_idx, ftype_str);
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_INVALID_BP_COORDINATE:
    LOGERRPRINTF("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, ftype_str);
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, ftype_str);
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_FROM_TO_DIFFERENT_CHROM:
    logerrprint("Error: --from and --to variants are not on the same chromosome.\n");
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_DUPLICATE_ID:
    uii = strlen_se(col2_ptr);
    col2_ptr[uii] = '\0';
    LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in %s.\n", col2_ptr, ftype_str);
  load_bim_ret_INVALID_FORMAT_2:
    logerrprintb();
  load_bim_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 load_bim_ret_1:
  fclose_cond(bimfile);
  free_cond(loadbuf2);
  free_cond(prev_new_id);
  free_cond(insert_buf[2]);
  free_cond(insert_buf[3]);
  return retval;
}

int32_t load_covars(char* covar_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, double missing_phenod, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t gxe_mcovar, uintptr_t* covar_ctx_ptr, char** covar_names_ptr, uintptr_t* max_covar_name_len_ptr, uintptr_t* pheno_nm, uintptr_t** covar_nm_ptr, double** covar_d_ptr, uintptr_t** gxe_covar_nm_ptr, uintptr_t** gxe_covar_c_ptr) {
  // similar to load_clusters() in plink_cluster.c
  // sex_nm and sex_male should be nullptr unless sex is supposed to be added
  // as an extra covariate
  // covar_range_list_ptr is nullptr iff --gxe was specified
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  unsigned char* bigstack_mark2 = nullptr;
  FILE* covar_file = nullptr;
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t covar_raw_ct = 0;
  uintptr_t loaded_sample_ct = 0;
  uintptr_t missing_cov_ct = 0;
  uint32_t* sample_idx_to_uidx = nullptr;
  char* sorted_covar_name_flag_ids = nullptr;
  uint32_t* covar_name_flag_id_map = nullptr;
  int32_t* covar_name_flag_seen_idxs = nullptr;
  char* covar_names = nullptr;
  uintptr_t* covar_nm = nullptr;
  double* covar_d = nullptr;
  uintptr_t* gxe_covar_nm = nullptr;
  uintptr_t* gxe_covar_c = nullptr;
  double* dptr = nullptr;
  char* bufptr = nullptr;
  uintptr_t max_covar_name_len = sex_nm? 4 : 1;
  double dxx = 0.0;
  uint32_t keep_pheno_on_missing_cov = covar_modifier & COVAR_KEEP_PHENO_ON_MISSING_COV;
  int32_t retval = 0;
  uintptr_t covar_raw_ctl;
  uintptr_t covar_ct;
  uintptr_t covar_ctx;
  uintptr_t* covars_active;
  uintptr_t* already_seen;
  char* sorted_ids;
  uint32_t* id_map;
  char* loadbuf;
  uintptr_t loadbuf_size;
  char* bufptr2;
  uintptr_t covar_uidx;
  uintptr_t covar_idx;
  uintptr_t line_idx;
  uintptr_t ulii;
  uint32_t header_absent;
  uint32_t min_covar_col_ct;
  uint32_t uii;
  uint32_t sample_idx;
  uint32_t sample_uidx;
  uint32_t covar_missing;
  int32_t ii;

  if ((!keep_pheno_on_missing_cov) || gxe_mcovar || sex_nm) {
    if (bigstack_end_alloc_ui(sample_ct, &sample_idx_to_uidx)) {
      goto load_covars_ret_NOMEM;
    }
    fill_idx_to_uidx(sample_exclude, unfiltered_sample_ct, sample_ct, sample_idx_to_uidx);
  }
  if (bigstack_end_alloc_c(sample_ct * max_sample_id_len, &sorted_ids) ||
      bigstack_end_alloc_ui(sample_ct, &id_map) ||
      bigstack_end_calloc_ul(sample_ctl, &already_seen)) {
    goto load_covars_ret_NOMEM;
  }
  if (covar_modifier & COVAR_NAME) {
    ulii = covar_range_list_ptr->name_ct;
    if (bigstack_end_alloc_c(ulii * covar_range_list_ptr->name_max_len, &sorted_covar_name_flag_ids) ||
	bigstack_end_alloc_ui(ulii, &covar_name_flag_id_map) ||
        bigstack_end_alloc_i(ulii, &covar_name_flag_seen_idxs)) {
      goto load_covars_ret_NOMEM;
    }

    // kludge to use sort_item_ids_noalloc()
    fill_ulong_zero(BITCT_TO_WORDCT(ulii), (uintptr_t*)covar_name_flag_seen_idxs);
    retval = sort_item_ids_noalloc(ulii, (const uintptr_t*)covar_name_flag_seen_idxs, ulii, covar_range_list_ptr->names, covar_range_list_ptr->name_max_len, 0, 0, strcmp_deref, sorted_covar_name_flag_ids, covar_name_flag_id_map);
    if (retval) {
      if (retval == RET_INVALID_FORMAT) {
	logprint("(in --covar-name parameter sequence)\n");
	retval = RET_INVALID_CMDLINE;
      }
      goto load_covars_ret_1;
    }
    fill_int_one(ulii, covar_name_flag_seen_idxs);
  }
  retval = sort_item_ids_noalloc(unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, 0, 1, strcmp_deref, sorted_ids, id_map);
  if (retval) {
    goto load_covars_ret_1;
  }

  // To simplify loading sequence, guarantee enough space for covars_active[]
  // bitfield on first pass.  Each covariate corresponds to at least 16 bits in
  // the first nonempty line (a value and a space = 2 bytes), so reserving the
  // last 1/17 (rounded up) always works.  (Minor memory leak fix:
  // covars_active no longer remains allocated on function exit.)
  loadbuf_size = (bigstack_left() / 68) * 64;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_covars_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  // was using open_and_load_to_first_token(), but we now don't want to
  // automatically print an error message on an empty file.
  loadbuf[loadbuf_size - 1] = ' ';
  if (fopen_checked(covar_fname, "r", &covar_file)) {
    goto load_covars_ret_OPEN_FAIL;
  }
  line_idx = 0;
  do {
    if (!fgets(loadbuf, loadbuf_size, covar_file)) {
      if (!feof(covar_file)) {
	goto load_covars_ret_READ_FAIL;
      }
      strcpy(g_textbuf, "Empty --covar file.\n");
      goto load_covars_none;
    }
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	LOGERRPRINTF("Error: Line %" PRIuPTR " of --covar file is pathologically long.\n", line_idx);
	goto load_covars_ret_INVALID_FORMAT;
      } else {
	goto load_covars_ret_NOMEM;
      }
    }
    bufptr = skip_initial_spaces(loadbuf);
  } while (is_eoln_kns(*bufptr));
  covar_raw_ct = count_tokens(bufptr);
  if ((covar_raw_ct < 2) || (covar_raw_ct < 2 + gxe_mcovar)) {
    goto load_covars_ret_MISSING_TOKENS;
  }
  covar_raw_ct -= 2;
  if ((!covar_raw_ct) && (!sex_nm)) {
    strcpy(g_textbuf, "No covariate columns in --covar file.\n");
    goto load_covars_none;
  }
  covar_raw_ctl = BITCT_TO_WORDCT(covar_raw_ct);
  if (bigstack_end_alloc_ul(covar_raw_ctl, &covars_active)) {
    goto load_covars_ret_NOMEM;
  }

  // no header line present?
  bufptr2 = next_token(bufptr);
  header_absent = ((strcmp_se(bufptr, "FID", 3) && (strcmp_se(bufptr, "#FID", 4))) || strcmp_se(bufptr2, "IID", 3));
  bufptr = next_token(bufptr2);

  if ((covar_modifier & (COVAR_NAME | COVAR_NUMBER)) && covar_raw_ct) {
    fill_ulong_zero(covar_raw_ctl, covars_active);
    if (covar_modifier & COVAR_NUMBER) {
      if (numeric_range_list_to_bitarr(covar_range_list_ptr, covar_raw_ct, 1, 0, covars_active)) {
	goto load_covars_ret_MISSING_TOKENS;
      }
    } else if (covar_modifier & COVAR_NAME) {
      if (header_absent) {
	logerrprint("Error: --covar file doesn't have a header line for --covar-name.\n");
	goto load_covars_ret_INVALID_FORMAT;
      }
      retval = string_range_list_to_bitarr(bufptr, covar_raw_ct, 0, covar_range_list_ptr, sorted_covar_name_flag_ids, covar_name_flag_id_map, "covar-name", "--covar file header line", covars_active, covar_name_flag_seen_idxs);
      if (retval) {
	goto load_covars_ret_1;
      }
      // can't deallocate --covar-name support here due to covars_active
      // repositioning
    }
    covar_ct = popcount_longs(covars_active, covar_raw_ctl);
  } else if (covar_range_list_ptr) {
    fill_all_bits(covar_raw_ct, covars_active);
    covar_ct = covar_raw_ct;
  } else {
    // --gxe only
    fill_ulong_zero(covar_raw_ctl, covars_active);
    covar_ct = 0;
  }
  covar_ctx = covar_ct + (sex_nm? 1 : 0);
  if ((!covar_ctx) && (!gxe_mcovar)) {
    strcpy(g_textbuf, "No --covar values loaded.\n");
    goto load_covars_none;
  }
  min_covar_col_ct = covar_ct? (last_set_bit(covars_active, covar_raw_ctl) + 1) : 0;
  if (min_covar_col_ct < gxe_mcovar) {
    min_covar_col_ct = gxe_mcovar;
  }
  if (header_absent) {
    max_covar_name_len = 4 + intlen(min_covar_col_ct);
  } else if (min_covar_col_ct) {
    uii = 0;
    while (1) {
      bufptr2 = token_endnn(bufptr);
      if (IS_SET(covars_active, uii)) {
        if (max_covar_name_len <= (uintptr_t)(bufptr2 - bufptr)) {
	  max_covar_name_len = 1 + (uintptr_t)(bufptr2 - bufptr);
	}
      }
      if (++uii == min_covar_col_ct) {
	break;
      }
      bufptr = skip_initial_spaces(&(bufptr2[1]));
    }
  }

  // * covar_nm does NOT have a separate entry per covariate; instead,
  //   if a single covariate is missing for a person, that person's covar_nm
  //   bit is zero.
  // * covar_d is in sample-major order (i.e. the value of covariate m for
  //   sample n is covar_d[n * covar_ctx + m], where m and n are zero-based).
  //   It does track when some covariates are missing and others aren't
  //   (missing covariates are represented as the --missing-phenotype value).
  if (covar_range_list_ptr) {
    if (max_covar_name_len > MAX_ID_BLEN) {
      logerrprint("Error: Covariate names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto load_covars_ret_INVALID_FORMAT;
    }
    // not only --gxe
    *covar_ctx_ptr = covar_ctx;
    *max_covar_name_len_ptr = max_covar_name_len;
    ulii = covar_ctx * sample_ct;
    if (bigstack_alloc_c(covar_ctx * max_covar_name_len, covar_names_ptr) ||
        bigstack_alloc_ul(sample_ctl, covar_nm_ptr) ||
        bigstack_alloc_d(ulii, covar_d_ptr)) {
      goto load_covars_ret_NOMEM;
    }
    covar_names = *covar_names_ptr;
    covar_nm = *covar_nm_ptr;
    covar_d = *covar_d_ptr;
    fill_ulong_zero(sample_ctl, covar_nm);
    for (covar_idx = 0; covar_idx < ulii; covar_idx++) {
      covar_d[covar_idx] = missing_phenod;
    }
  }
  if (gxe_mcovar) {
    if (bigstack_calloc_ul(sample_ctl, gxe_covar_nm_ptr) ||
        bigstack_calloc_ul(sample_ctl, gxe_covar_c_ptr)) {
      goto load_covars_ret_NOMEM;
    }
    gxe_covar_nm = *gxe_covar_nm_ptr;
    gxe_covar_c = *gxe_covar_c_ptr;
  }
  loadbuf_size = bigstack_left();
  if (loadbuf_size <= MAXLINELEN) {
    goto load_covars_ret_NOMEM;
  }
  bigstack_mark2 = g_bigstack_base;
  loadbuf = (char*)g_bigstack_base;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  loadbuf[loadbuf_size - 1] = ' ';

  rewind(covar_file);
  if (header_absent) {
    if (covar_range_list_ptr) {
      for (covar_uidx = 0, covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	covar_uidx = next_set_ul_unsafe(covars_active, covar_uidx);
	uint32toa_x(++covar_uidx, '\0', memcpyl3a(&(covar_names[covar_idx * max_covar_name_len]), "COV"));
      }
    }
    line_idx = 0;
  } else if (covar_ct) {
    covar_idx = 0;
    retval = load_to_first_token(covar_file, loadbuf_size, '\0', "--covar file", loadbuf, &bufptr, &line_idx);
    if (retval) {
      goto load_covars_ret_1;
    }
    if (covar_range_list_ptr) {
      bufptr2 = token_endnn(next_token(bufptr));
      for (uii = 0; uii < min_covar_col_ct; uii++) {
	bufptr = skip_initial_spaces(bufptr2);
	bufptr2 = token_endnn(bufptr);
	if (IS_SET(covars_active, uii)) {
	  memcpyx(&(covar_names[covar_idx * max_covar_name_len]), bufptr, (uintptr_t)(bufptr2 - bufptr), '\0');
	  covar_idx++;
	}
      }
    }
  }
  if (sex_nm) {
    memcpy(&(covar_names[covar_ct * max_covar_name_len]), "SEX", 4);
  }
  while (fgets(loadbuf, loadbuf_size, covar_file)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --covar file is pathologically long.\n", line_idx);
	goto load_covars_ret_INVALID_FORMAT_2;
      } else {
	goto load_covars_ret_NOMEM;
      }
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(bufptr, sorted_ids, max_sample_id_len, sample_ct, &bufptr2, &ii, g_textbuf)) {
      goto load_covars_ret_MISSING_TOKENS;
    }
    if (ii == -1) {
      continue;
    }
    if (is_set(already_seen, ii)) {
      logerrprint("Error: Duplicate sample ID in --covar file.\n");
      goto load_covars_ret_INVALID_FORMAT;
    }
    set_bit(ii, already_seen);
    sample_idx = id_map[(uint32_t)ii];
    bufptr = bufptr2;
    if (min_covar_col_ct) {
      if (min_covar_col_ct > 1) {
	bufptr = next_token_mult(bufptr, min_covar_col_ct - 1);
      }
      if (no_more_tokens_kns(bufptr)) {
	goto load_covars_ret_MISSING_TOKENS;
      }
    }
    if (covar_range_list_ptr) {
      dptr = &(covar_d[sample_idx * covar_ctx]);
    }
    covar_missing = 0;
    for (uii = 0; uii < min_covar_col_ct; uii++) {
      bufptr = skip_initial_spaces(bufptr2);

      // column count already validated
      bufptr2 = token_endnn(bufptr);

      if (IS_SET(covars_active, uii)) {
        if (scan_double(bufptr, &dxx)) {
	  covar_missing = 1;
	  dxx = missing_phenod;
	} else if (dxx == missing_phenod) {
	  covar_missing = 1;
	}
	*dptr++ = dxx;
      }
      if (uii + 1 == gxe_mcovar) {
        sample_uidx = sample_idx_to_uidx[sample_idx];
	if (IS_SET(pheno_nm, sample_uidx)) {
	  if (scan_double(bufptr, &dxx) || (dxx == missing_phenod)) {
	    // PLINK 1.07 quirk: --gxe treats 0 the same as -9/nonnumeric, but
	    // --write-covar does not, so handle 0 separately for backward
	    // compatibility
	    if (!keep_pheno_on_missing_cov) {
	      CLEAR_BIT(sample_uidx, pheno_nm);
	    }
	  } else if (dxx != 0.0) {
	    SET_BIT(sample_idx, gxe_covar_nm);
	    if (dxx == 2.0) {
	      SET_BIT(sample_idx, gxe_covar_c);
	    }
	  }
	}
      }
    }
    if (covar_range_list_ptr) {
      if (sex_nm) {
	sample_uidx = sample_idx_to_uidx[sample_idx];
	if (is_set(sex_nm, sample_uidx)) {
	  *dptr++ = (double)((int32_t)is_set(sex_male, sample_uidx));
	} else {
	  covar_missing = 1;
	  *dptr++ = missing_phenod;
	}
      }
      if (!covar_missing) {
	SET_BIT(sample_idx, covar_nm);
      } else {
	missing_cov_ct++;
	if (!keep_pheno_on_missing_cov) {
	  sample_uidx = sample_idx_to_uidx[sample_idx];
	  if (IS_SET(pheno_nm, sample_uidx)) {
	    CLEAR_BIT(sample_uidx, pheno_nm);
	  }
	}
      }
    }
    loaded_sample_ct++;
  }
  if (!feof(covar_file)) {
    goto load_covars_ret_READ_FAIL;
  }
  if (covar_range_list_ptr) {
    if ((covar_ct + 1 < covar_raw_ct) || ((covar_ct + 1 == covar_raw_ct) && ((!gxe_mcovar) || is_set(covars_active, gxe_mcovar - 1)))) {
      if (gxe_mcovar && (!is_set(covars_active, gxe_mcovar - 1))) {
        sprintf(g_logbuf, "--covar: 1 C/C cov. loaded for --gxe, %" PRIuPTR "/%" PRIuPTR " for other operations.\n", covar_ct, covar_raw_ct);
      } else {
        sprintf(g_logbuf, "--covar: %" PRIuPTR " out of %" PRIuPTR " covariates loaded.\n", covar_ct, covar_raw_ct);
      }
    } else {
      sprintf(g_logbuf, "--covar: %" PRIuPTR " covariate%s loaded.\n", covar_ct, (covar_ct == 1)? "" : "s");
    }
    logprintb();
  } else {
    logprint("--covar: 1 case/control covariate loaded for --gxe.\n");
  }
  ulii = sample_ct - loaded_sample_ct;
  if (missing_cov_ct) {
    LOGPRINTF("%" PRIuPTR " %s had missing value(s)%s", missing_cov_ct, species_str(missing_cov_ct), ulii? ", and " : ".\n");
  }
  if (ulii) {
    LOGPRINTF("%" PRIuPTR " %s %s not seen in the covariate file.\n", ulii, species_str(ulii), (ulii == 1)? "was" : "were");
  }

  if (covar_modifier & COVAR_NO_CONST) {
    if (gxe_mcovar) {
      uii = popcount_longs(gxe_covar_c, sample_ctl);
      if ((!uii) || (uii == popcount_longs(gxe_covar_nm, sample_ctl))) {
	logerrprint("Error: --gxe covariate is constant and --no-const-covar was specified.\n");
	goto load_covars_ret_INVALID_FORMAT;
      }
    }
    if (covar_range_list_ptr) {
      // redefinition
      covar_raw_ctl = BITCT_TO_WORDCT(covar_ctx);
      if (bigstack_calloc_ul(covar_raw_ctl, &already_seen)) {
	goto load_covars_ret_NOMEM;
      }
      // is covariate nonconstant?
      for (covar_idx = 0; covar_idx < covar_ctx; covar_idx++) {
	dptr = &(covar_d[covar_idx]);
	dxx = missing_phenod;
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	  if (dptr[sample_idx * covar_ctx] != missing_phenod) {
	    dxx = dptr[sample_idx * covar_ctx];
	    break;
	  }
	}
	for (; sample_idx < sample_ct; sample_idx++) {
	  if ((dptr[sample_idx * covar_ctx] != missing_phenod) && (dptr[sample_idx * covar_ctx] != dxx)) {
	    break;
	  }
	}
	if (sample_idx < sample_ct) {
	  SET_BIT(covar_idx, already_seen);
	}
      }
      uii = popcount_longs(already_seen, covar_raw_ctl);
      if (!uii) {
	strcpy(g_textbuf, "All covariates are constant.\n");
	goto load_covars_none;
      } else if (uii < covar_ctx) {
	LOGPRINTF("--no-const-covar: %" PRIuPTR " constant covariate%s excluded.\n", covar_ctx - uii, (covar_ctx - uii == 1)? "" : "s");
	*covar_ctx_ptr = uii;
	dptr = covar_d;
        for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	  uii = 0;
	  for (covar_idx = 0; covar_idx < covar_ctx; covar_idx++) {
	    if (IS_SET(already_seen, covar_idx)) {
	      dxx = covar_d[sample_idx * covar_ctx + covar_idx];
	      if (dxx == missing_phenod) {
		uii = 1;
	      }
	      *dptr++ = dxx;
	    }
	  }
	  if (!uii) {
	    // if this sample had some missing covariate values, but all those
	    // covariates were excluded by --no-const-covar, set covar_nm bit
	    SET_BIT(sample_idx, covar_nm);
	  }
	}
	covar_idx = next_unset_unsafe(already_seen, 0);
	uii = covar_idx;
	for (; covar_idx < covar_ctx; covar_idx++) {
	  if (IS_SET(already_seen, covar_idx)) {
	    strcpy(&(covar_names[uii * max_covar_name_len]), &(covar_names[covar_idx * max_covar_name_len]));
	    uii++;
	  }
	}
	// don't worry about memory overallocation for now
      }
    }
  }

  bigstack_reset(bigstack_mark2);
  while (0) {
  load_covars_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_covars_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_covars_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_covars_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Fewer tokens than expected on line %" PRIuPTR " of --covar file.\n", line_idx);
  load_covars_ret_INVALID_FORMAT_2:
    logerrprintb();
  load_covars_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  load_covars_none:
    if (covar_modifier & COVAR_ALLOW_NONE) {
      *covar_ctx_ptr = 0;
      *covar_names_ptr = nullptr;
      *max_covar_name_len_ptr = 1;
      *covar_nm_ptr = nullptr;
      *covar_d_ptr = nullptr;
      // --gxe not possible
      bigstack_reset(bigstack_mark);
      logerrprint("Warning: ");
    } else {
      retval = RET_INVALID_FORMAT;
      logerrprint("Error: ");
    }
    logerrprint(g_textbuf);
  }
 load_covars_ret_1:
  if (retval) {
    bigstack_reset(bigstack_mark);
  }
  bigstack_end_reset(bigstack_end_mark);
  fclose_cond(covar_file);
  return retval;
}

int32_t write_covars(char* outname, char* outname_end, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* outfile = nullptr;
  uint32_t write_pheno = write_covar_modifier & WRITE_COVAR_PHENO;
  uint32_t exclude_parents = write_covar_modifier & WRITE_COVAR_NO_PARENTS;
  uint32_t exclude_sex = write_covar_modifier & WRITE_COVAR_NO_SEX;
  uint32_t female_2 = write_covar_modifier & WRITE_COVAR_FEMALE_2;
  uintptr_t sample_uidx = 0;
  uint32_t* downcoding_level = nullptr;
  uint32_t* downcoding_values = nullptr;
  char* zbuf = nullptr;
  char* out_missing_buf = nullptr;
  uintptr_t omplen_p1 = strlen(output_missing_pheno) + 1;
  uint32_t do_downcoding = (write_covar_modifier & WRITE_COVAR_DUMMY) && (sample_ct > 2);
  uint32_t downcoding_no_round = (write_covar_modifier & WRITE_COVAR_DUMMY_NO_ROUND);
  uintptr_t downcoding_covar_ct = 0;
  uintptr_t covar_nm_ct = 0;
  int32_t retval = 0;
  uint32_t* downcoding_buf_idxs;
  int64_t* sorted_downcoding_intbuf;
  int64_t* category_idx_sort_buf;
  uint32_t* category_remap;
  uint32_t* uiptr;
  double* sorted_downcoding_buf;
  double* dptr;
  char* downcoding_string_buf;
  char* wptr_start;
  char* wptr;
  char* bufptr;
  uintptr_t covar_idx;
  uintptr_t sample_idx;
  uintptr_t sample_idx2;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t slen;
  uint64_t ullii;
  uint64_t ulljj;
  uintptr_t downcode_category_ct;
  uintptr_t downcode_idx;
  double dxx;
  uint32_t uii;
  uint32_t ujj;
  memcpy(outname_end, ".cov", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_covars_ret_OPEN_FAIL;
  }
  if (fputs_checked("FID IID ", outfile)) {
    goto write_covars_ret_WRITE_FAIL;
  }
  if (write_pheno) {
    if (!exclude_parents) {
      fputs("PAT MAT ", outfile);
    }
    if (!exclude_sex) {
      fputs("SEX PHENOTYPE ", outfile);
    } else {
      fputs("PHENOTYPE ", outfile);
    }
  }
  if (do_downcoding) {
    // could make downcoding_values allocation incremental
    // (bigstack_end_alloc() calls have been arranged to make this a simple
    // change; would just need to wrap the qsort_ext() calls)
    if (bigstack_alloc_ui(covar_ct, &downcoding_level) ||
        bigstack_alloc_ui(covar_ct * sample_ct, &downcoding_values)) {
      goto write_covars_ret_NOMEM;
    }
    if (write_covar_dummy_max_categories > sample_ct) {
      write_covar_dummy_max_categories = sample_ct;
    }
    if (bigstack_end_alloc_c(16 * (write_covar_dummy_max_categories + 1), &downcoding_string_buf) ||
        bigstack_end_alloc_ll(write_covar_dummy_max_categories, &category_idx_sort_buf) ||
	bigstack_end_alloc_ui(write_covar_dummy_max_categories, &category_remap)) {
      goto write_covars_ret_NOMEM;
    }
    uiptr = downcoding_values;
    if (!downcoding_no_round) {
      if (bigstack_end_alloc_ll(sample_ct, &sorted_downcoding_intbuf)) {
        goto write_covars_ret_NOMEM;
      }
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
        dptr = &(covar_d[covar_idx]);
	sample_idx2 = 0;
        for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	  dxx = dptr[sample_idx * covar_ct];
	  if (dxx != missing_phenod) {
            // bits 0-30: original position
            // bits 31-62: 2^31 + (int32_t)phenotype
            // top bit zero to defend against >= 2^31 differences
            sorted_downcoding_intbuf[sample_idx2++] = (int64_t)((((uint64_t)(((uint32_t)((int32_t)dxx)) ^ 0x80000000U)) << 31) | ((uint64_t)sample_idx));
	  }
	}
	downcode_category_ct = 0;
	if (sample_idx2 > 2) {
	  covar_nm_ct = sample_idx2;
	  fill_uint_one(sample_ct, uiptr);
#ifdef __cplusplus
	  std::sort(sorted_downcoding_intbuf, &(sorted_downcoding_intbuf[covar_nm_ct]));
#else
          qsort(sorted_downcoding_intbuf, covar_nm_ct, sizeof(int64_t), llcmp);
#endif

          ullii = ((uint64_t)sorted_downcoding_intbuf[0]);
	  ulii = ((uintptr_t)ullii) & (ONELU * 0x7fffffff);
	  ullii >>= 31;
          uiptr[ulii] = 0;
          for (sample_idx2 = 1; sample_idx2 < covar_nm_ct; sample_idx2++) {
	    ulljj = (uint64_t)sorted_downcoding_intbuf[sample_idx2];
	    uljj = ((uintptr_t)ulljj) & (ONELU * 0x7fffffff);
	    ulljj >>= 31;
	    if (ulljj != ullii) {
	      if (downcode_category_ct == write_covar_dummy_max_categories - 1) {
		downcode_category_ct = 0;
		break;
	      }
	      // save phenotype string
	      int32toa_x((int32_t)(((uint32_t)ullii) ^ 0x80000000U), '\0', &(downcoding_string_buf[16 * downcode_category_ct]));

	      // bits 0-31: initial category assignment
	      // bits 32-63: smallest sample_idx2
	      category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)ulii) << 32) | ((uint64_t)downcode_category_ct));
	      downcode_category_ct++;
	      ulii = uljj;
	      ullii = ulljj;
	    } else {
	      if (uljj < ulii) {
		ulii = uljj;
	      }
	    }
	    uiptr[uljj] = downcode_category_ct;
	  }
	}

	// probably want to make this part its own function since it's
	// practically identical when downcoding_no_round is set
        if (downcode_category_ct > 1) {
	  int32toa_x((int32_t)(((uint32_t)ullii) ^ 0x80000000U), '\0', &(downcoding_string_buf[16 * downcode_category_ct]));
	  category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)ulii) << 32) | ((uint64_t)downcode_category_ct));
	  downcode_category_ct++;
          // now recover PLINK 1.07 category order
#ifdef __cplusplus
	  std::sort(category_idx_sort_buf, &(category_idx_sort_buf[downcode_category_ct]));
#else
          qsort(category_idx_sort_buf, downcode_category_ct, sizeof(int64_t), llcmp);
#endif

	  downcoding_level[covar_idx] = downcode_category_ct;
          wptr_start = strcpyax(g_textbuf, &(covar_names[covar_idx * max_covar_name_len]), '_');
	  for (downcode_idx = 0; downcode_idx < downcode_category_ct; downcode_idx++) {
	    uii = (uint32_t)(category_idx_sort_buf[downcode_idx]);
	    if (downcode_idx) {
	      wptr = strcpyax(wptr_start, &(downcoding_string_buf[16 * uii]), ' ');
	      fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	    }
            category_remap[uii] = downcode_idx;
	  }
          for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    if (uiptr[sample_idx] != 0xffffffffU) {
              uiptr[sample_idx] = category_remap[uiptr[sample_idx]];
	    }
	  }
	  downcoding_covar_ct++;
	  uiptr = &(uiptr[sample_ct]);
	} else {
          downcoding_level[covar_idx] = 0;
          fputs(&(covar_names[covar_idx * max_covar_name_len]), outfile);
	  putc_unlocked(' ', outfile);
	}
      }
    } else {
      if (bigstack_end_alloc_ui(sample_ct, &downcoding_buf_idxs) ||
          bigstack_end_alloc_d(sample_ct, &sorted_downcoding_buf)) {
	goto write_covars_ret_NOMEM;
      }
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	dptr = &(covar_d[covar_idx]);
	sample_idx2 = 0;
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	  dxx = dptr[sample_idx * covar_ct];
	  if (dxx != missing_phenod) {
	    sorted_downcoding_buf[sample_idx2] = dxx;
	    downcoding_buf_idxs[sample_idx2++] = sample_idx;
	  }
	}
	downcode_category_ct = 0;
	if (sample_idx2 > 2) {
	  covar_nm_ct = sample_idx2;
	  fill_uint_one(sample_ct, uiptr);
	  if (qsort_ext((char*)sorted_downcoding_buf, covar_nm_ct, sizeof(double), double_cmp_deref, (char*)downcoding_buf_idxs, sizeof(int32_t))) {
	    goto write_covars_ret_NOMEM;
	  }
	  wptr_start = dtoa_g(sorted_downcoding_buf[0], downcoding_string_buf);
	  slen = (uintptr_t)(wptr_start - downcoding_string_buf);
	  *wptr_start = '\0';
	  wptr_start = downcoding_string_buf;
	  uii = downcoding_buf_idxs[0];
	  uiptr[uii] = 0;
	  bufptr = &(downcoding_string_buf[16]);
	  for (sample_idx2 = 1; sample_idx2 < covar_nm_ct; sample_idx2++) {
	    // a bit inefficient, but this is a safe way to achieve "two
	    // doubles are equal if they yield the same printf %g output"
	    // behavior
	    wptr = dtoa_g(sorted_downcoding_buf[sample_idx2], bufptr);
	    ujj = downcoding_buf_idxs[sample_idx2];
	    if (((uintptr_t)(wptr - bufptr) != slen) || memcmp(wptr_start, bufptr, slen)) {
	      *wptr = '\0';
	      if (downcode_category_ct == write_covar_dummy_max_categories - 1) {
		downcode_category_ct = 0;
		break;
	      }
	      category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)uii) << 32) | ((uint64_t)downcode_category_ct));
	      downcode_category_ct++;
	      slen = (uintptr_t)(wptr - bufptr);
	      *wptr = '\0';
	      wptr_start = bufptr;
              uii = ujj;
	      bufptr = &(bufptr[16]);
	    } else {
	      if (ujj < uii) {
		uii = ujj;
	      }
	    }
	    uiptr[ujj] = downcode_category_ct;
	  }
	}
	if (downcode_category_ct > 1) {
	  *wptr = '\0';
	  category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)uii) << 32) | ((uint64_t)downcode_category_ct));
	  downcode_category_ct++;

#ifdef __cplusplus
	  std::sort(category_idx_sort_buf, &(category_idx_sort_buf[downcode_category_ct]));
#else
          qsort(category_idx_sort_buf, downcode_category_ct, sizeof(int64_t), llcmp);
#endif

	  downcoding_level[covar_idx] = downcode_category_ct;
          wptr_start = strcpyax(g_textbuf, &(covar_names[covar_idx * max_covar_name_len]), '_');
          for (downcode_idx = 0; downcode_idx < downcode_category_ct; downcode_idx++) {
	    uii = (uint32_t)(category_idx_sort_buf[downcode_idx]);
	    if (downcode_idx) {
	      wptr = strcpyax(wptr_start, &(downcoding_string_buf[16 * uii]), ' ');
	      fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	    }
	    category_remap[uii] = downcode_idx;
          }
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    if (uiptr[sample_idx] != 0xffffffffU) {
	      uiptr[sample_idx] = category_remap[uiptr[sample_idx]];
	    }
	  }
	  downcoding_covar_ct++;
	  uiptr = &(uiptr[sample_ct]);
	} else {
	  downcoding_level[covar_idx] = 0;
          fputs(&(covar_names[covar_idx * max_covar_name_len]), outfile);
	  putc_unlocked(' ', outfile);
	}
      }
    }
    bigstack_shrink_top(downcoding_values, downcoding_covar_ct * sample_ct * sizeof(int32_t));
    bigstack_end_reset(bigstack_end_mark);

    // (write_covar_dummy_max_categories - 1) columns, then divide by two
    // rounding up; the -1 and +1 cancel
    ujj = write_covar_dummy_max_categories / 2;
    if (bigstack_alloc_c(ujj * sizeof(int32_t), &zbuf) ||
        bigstack_alloc_c((write_covar_dummy_max_categories - 1) * omplen_p1, &out_missing_buf)) {
      goto write_covars_ret_NOMEM;
    }
    uiptr = (uint32_t*)zbuf;
    for (uii = 0; uii < ujj; uii++) {
      *uiptr++ = 0x20302030; // "0 0 "
    }
    bufptr = memcpyax(out_missing_buf, output_missing_pheno, omplen_p1 - 1, ' ');
    for (uii = 2; uii < write_covar_dummy_max_categories; uii++) {
      bufptr = memcpya(bufptr, out_missing_buf, omplen_p1);
    }
  } else {
    for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
      fputs(&(covar_names[covar_idx * max_covar_name_len]), outfile);
      putc_unlocked(' ', outfile);
    }
  }
  if (putc_checked('\n', outfile)) {
    goto write_covars_ret_WRITE_FAIL;
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    sample_uidx = next_unset_ul_unsafe(sample_exclude, sample_uidx);
    bufptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    wptr = (char*)memchr(bufptr, '\t', max_sample_id_len);
    *wptr = ' ';
    if (fputs_checked(bufptr, outfile)) {
      goto write_covars_ret_WRITE_FAIL;
    }
    *wptr = '\t';
    putc_unlocked(' ', outfile);
    if (write_pheno) {
      if (!exclude_parents) {
	fputs(&(paternal_ids[sample_uidx * max_paternal_id_len]), outfile);
	putc_unlocked(' ', outfile);
	fputs(&(maternal_ids[sample_uidx * max_maternal_id_len]), outfile);
	putc_unlocked(' ', outfile);
      }
      if (!exclude_sex) {
	if (!female_2) {
          putc_unlocked((uint32_t)(48 + IS_SET(sex_male, sample_uidx)), outfile);
	} else {
          if (IS_SET(sex_nm, sample_uidx)) {
	    putc_unlocked((uint32_t)(50 - IS_SET(sex_male, sample_uidx)), outfile);
	  } else {
	    putc_unlocked('0', outfile);
	  }
	}
        putc_unlocked(' ', outfile);
      }
      if (!IS_SET(pheno_nm, sample_uidx)) {
        fputs(output_missing_pheno, outfile);
      } else if (pheno_c) {
        putc_unlocked('1' + IS_SET(pheno_c, sample_uidx), outfile);
      } else {
        wptr = dtoa_g(pheno_d[sample_uidx], g_textbuf);
	fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
      }
      putc_unlocked(' ', outfile);
    }
    if (do_downcoding) {
      dptr = &(covar_d[sample_idx * covar_ct]);
      downcode_idx = 0;
      uiptr = &(downcoding_values[sample_idx]);
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	ujj = downcoding_level[covar_idx];
	dxx = *dptr++;
	if (dxx != missing_phenod) {
	  if (ujj) {
	    uii = uiptr[downcode_idx * sample_ct];
	    if (uii) {
	      if (uii > 1) {
		if (fwrite_checked(zbuf, (2 * ONELU) * (uii - 1), outfile)) {
		  goto write_covars_ret_WRITE_FAIL;
		}
	      }
	      fputs("1 ", outfile);
	      if (uii < ujj - 1) {
		if (fwrite_checked(zbuf, (2 * ONELU) * (ujj - 1 - uii), outfile)) {
		  goto write_covars_ret_WRITE_FAIL;
		}
	      }
	    } else {
	      if (fwrite_checked(zbuf, (2 * ONELU) * (ujj - 1), outfile)) {
		goto write_covars_ret_WRITE_FAIL;
	      }
	    }
	    downcode_idx++;
	  } else {
	    wptr = dtoa_gx(dxx, ' ', g_textbuf);
	    fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	  }
	} else {
	  if (ujj) {
	    if (fwrite_checked(out_missing_buf, omplen_p1 * (ujj - 1), outfile)) {
	      goto write_covars_ret_WRITE_FAIL;
	    }
	  } else {
            fputs(output_missing_pheno, outfile);
            putc_unlocked(' ', outfile);
	  }
	}
      }
    } else {
      dptr = &(covar_d[sample_idx * covar_ct]);
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	wptr = dtoa_gx(dptr[covar_idx], ' ', g_textbuf);
        fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
      }
    }
    if (putc_checked('\n', outfile)) {
      goto write_covars_ret_WRITE_FAIL;
    }
  }
  LOGPRINTFWW("Covariates written to %s .\n", outname);
  while (0) {
  write_covars_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_covars_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_covars_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t zero_cluster_init(char* zerofname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t* sample_sort_map, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t** zcdefs, uintptr_t** cluster_zc_masks_ptr) {
  // (todo: patch plink() sequence so that marker ID hash table only needs to
  // be constructed once?)
  //
  // 1. create marker ID hash table and sorted cluster ID list, regular stack
  //    allocation
  // 2. load .zero file, converting to internal indices.  (lines with
  //    unrecognized IDs are skipped; we don't want a header line to cause this
  //    to error out.)  this is bigstack_end_alloc'd.
  // 3. free marker ID/cluster ID lists, sort loaded .zero contents
  // 4. assemble one block bitfield at a time, use save_set_bitfield() to
  //    compress each
  // 5. allocate and initialize cluster_zc_masks
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* zcfile = nullptr;
  uintptr_t marker_ctp2l = (marker_ct + (BITCT + 1)) / BITCT;
  uintptr_t sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  uintptr_t zc_item_ct = 0;
  uintptr_t line_idx = 0;
  uint32_t range_first = marker_ct;
  uint32_t range_last = 0;
  int32_t retval = 0;
  uint32_t* marker_id_htable;
  uint32_t* marker_uidx_to_idx;
  uint32_t* sample_uidx_to_idx;
  uintptr_t* marker_bitfield_tmp;
  uintptr_t* cluster_zc_mask;
  char* bufptr;
  char* bufptr2;
  int64_t* zc_entries;
  int64_t* zc_entries_end;
  uint64_t ullii;
  uintptr_t max_zc_item_ct;
  uintptr_t marker_uidx;
  uint32_t marker_id_htable_size;
  uint32_t cluster_idx;
  uint32_t cur_cluster;
  uint32_t cluster_size;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t uii;
  int32_t ii;
  if (bigstack_end_alloc_ul(marker_ctp2l, &marker_bitfield_tmp) ||
      (!bigstack_end_alloc(sizeof(int64_t)))) {
    goto zero_cluster_init_ret_NOMEM;
  }
#ifdef __LP64__
  fill_ulong_zero(round_up_pow2(marker_ctp2l, 2), marker_bitfield_tmp);
#else
  fill_ulong_zero(round_up_pow2(marker_ctp2l, 4), marker_bitfield_tmp);
#endif
  zc_entries_end = (int64_t*)marker_bitfield_tmp;
  zc_entries = &(zc_entries_end[-1]);
  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, &marker_id_htable_size, &marker_id_htable);
  if (retval) {
    goto zero_cluster_init_ret_1;
  }
  if (bigstack_alloc_ui(unfiltered_marker_ct, &marker_uidx_to_idx)) {
    goto zero_cluster_init_ret_NOMEM;
  }
  fill_uidx_to_idx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_uidx_to_idx);
  // cluster IDs are already natural-sorted

  if (fopen_checked(zerofname, "r", &zcfile)) {
    goto zero_cluster_init_ret_OPEN_FAIL;
  }
  // simplify cluster_idx loop
  *zc_entries = (int64_t)(((uint64_t)cluster_ct) << 32);
  max_zc_item_ct = (((uintptr_t)zc_entries) - ((uintptr_t)g_bigstack_base)) / sizeof(int64_t);
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, zcfile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --zero-cluster file is pathologically long.\n", line_idx);
      goto zero_cluster_init_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = token_endnn(bufptr);
    marker_uidx = id_htable_find(bufptr, (uintptr_t)(bufptr2 - bufptr), marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx != 0xffffffffU) {
      bufptr = skip_initial_spaces(bufptr2);
      if (is_eoln_kns(*bufptr)) {
	goto zero_cluster_init_ret_MISSING_TOKENS;
      }
      bufptr2 = token_endnn(bufptr);
      uii = (uintptr_t)(bufptr2 - bufptr);
      if (uii < max_cluster_id_len) {
	*bufptr2 = '\0';
	ii = bsearch_str_natural(bufptr, cluster_ids, max_cluster_id_len, cluster_ct);
	if (ii != -1) {
	  if (++zc_item_ct > max_zc_item_ct) {
	    goto zero_cluster_init_ret_NOMEM;
	  }
	  // cluster ID in high bits, marker ID in low bits, so sorting works
	  *(--zc_entries) = (int64_t)((((uint64_t)((uint32_t)ii)) << 32) | ((uint64_t)marker_uidx_to_idx[marker_uidx]));
	}
      }
    }
  }
  if (fclose_null(&zcfile)) {
    goto zero_cluster_init_ret_READ_FAIL;
  }
  bigstack_double_reset(marker_id_htable, marker_bitfield_tmp);
  bigstack_end_alloc(zc_item_ct * sizeof(int64_t));
#ifdef __cplusplus
  std::sort(zc_entries, &(zc_entries[zc_item_ct]));
#else
  qsort(zc_entries, zc_item_ct, sizeof(int64_t), llcmp);
#endif
  ullii = (uint64_t)(*zc_entries++);
  cur_cluster = (uint32_t)(ullii >> 32);
  range_first = marker_ct;
  range_last = 0;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    if (range_first < marker_ct) {
      fill_ulong_zero(marker_ctp2l, marker_bitfield_tmp);
      range_first = marker_ct;
      range_last = 0;
      bigstack_end_set(zc_entries);
    }
    if (cur_cluster == cluster_idx) {
      range_first = (uint32_t)ullii;
      do {
	range_last = (uint32_t)ullii;
        SET_BIT(range_last, marker_bitfield_tmp);
        ullii = (uint64_t)(*zc_entries++);
        cur_cluster = (uint32_t)(ullii >> 32);
      } while (cur_cluster == cluster_idx);
    }
    if (save_set_bitfield(marker_bitfield_tmp, marker_ct, range_first, range_last + 1, 0, &(zcdefs[cluster_idx]))) {
      goto zero_cluster_init_ret_NOMEM;
    }
  }
  bigstack_end_reset(bigstack_end_mark);
  if (bigstack_calloc_ul(sample_ctv2 * cluster_ct, cluster_zc_masks_ptr) ||
      bigstack_alloc_ui(unfiltered_sample_ct, &sample_uidx_to_idx)) {
    goto zero_cluster_init_ret_NOMEM;
  }
  cluster_zc_mask = *cluster_zc_masks_ptr;
  if (!sample_sort_map) {
    fill_uidx_to_idx(sample_exclude, unfiltered_sample_ct, sample_ct, sample_uidx_to_idx);
  } else {
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      sample_uidx_to_idx[sample_sort_map[sample_idx]] = sample_idx;
    }
  }
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    // 01 if in cluster, 00 otherwise
    cluster_size = cluster_starts[cluster_idx + 1] - cluster_starts[cluster_idx];
    for (uii = 0; uii < cluster_size; uii++) {
      sample_uidx = *cluster_map++;
      if (!IS_SET(sample_exclude, sample_uidx)) {
	sample_idx = sample_uidx_to_idx[sample_uidx];
        SET_BIT_DBL(sample_idx, cluster_zc_mask);
      }
    }
    cluster_zc_mask = &(cluster_zc_mask[sample_ctv2]);
  }
  bigstack_reset(sample_uidx_to_idx);
  LOGPRINTF("--zero-cluster: %" PRIuPTR " line%s processed.\n", zc_item_ct, (zc_item_ct == 1)? "" : "s");
  while (0) {
  zero_cluster_init_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  zero_cluster_init_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  zero_cluster_init_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  zero_cluster_init_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --zero-cluster file has fewer tokens than expected.\n", line_idx);
  zero_cluster_init_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 zero_cluster_init_ret_1:
  bigstack_end_reset(bigstack_end_mark);
  fclose_cond(zcfile);
  return retval;
}

int32_t write_fam(char* outname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, char delim, uint32_t* sample_sort_map) {
  FILE* outfile = nullptr;
  uintptr_t sample_uidx = 0;
  uintptr_t sample_uidx2 = 0;
  uintptr_t omplen = strlen(output_missing_pheno);
  int32_t retval = 0;
  char* cptr;
  char* bufptr;
  uintptr_t sample_idx;
  uintptr_t clen;
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_fam_ret_OPEN_FAIL;
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    if (sample_sort_map) {
      do {
	sample_uidx = sample_sort_map[sample_uidx2++];
      } while (IS_SET(sample_exclude, sample_uidx));
    } else {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    }
    cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    clen = strlen_se(cptr);
    bufptr = memcpyax(g_textbuf, cptr, clen, delim);
    bufptr = strcpyax(bufptr, &(cptr[clen + 1]), delim);
    bufptr = strcpya(bufptr, &(paternal_ids[sample_uidx * max_paternal_id_len]));
    *bufptr++ = delim;
    bufptr = strcpya(bufptr, &(maternal_ids[sample_uidx * max_maternal_id_len]));
    *bufptr++ = delim;
    *bufptr++ = sexchar(sex_nm, sex_male, sample_uidx);
    *bufptr++ = delim;
    if (!IS_SET(pheno_nm, sample_uidx)) {
      bufptr = memcpya(bufptr, output_missing_pheno, omplen);
    } else if (pheno_c) {
      *bufptr++ = '1' + IS_SET(pheno_c, sample_uidx);
    } else {
      bufptr = dtoa_g(pheno_d[sample_uidx], bufptr);
    }
    *bufptr++ = '\n';
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
      goto write_fam_ret_WRITE_FAIL;
    }
    if (!sample_sort_map) {
      sample_uidx++;
    }
  }
  if (fclose_null(&outfile)) {
    goto write_fam_ret_WRITE_FAIL;
  }
  while (0) {
  write_fam_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_fam_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  return retval;
}

int32_t write_map_or_bim(char* outname, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, char delim, Chrom_info* chrom_info_ptr) {
  // write a .map if marker_allele_ptrs is nullptr, .bim otherwise
  FILE* outfile = nullptr;
  uintptr_t marker_uidx = 0;
  int32_t retval = 0;
  uint32_t chrom_end = 0;
  uint32_t chrom_fo_idx = 0xffffffffU;
  uint32_t chrom_idx = 0;
  const char* output_missing_geno_ptr = g_output_missing_geno_ptr;
  const char* missing_geno_ptr = (g_missing_geno_ptr == g_output_missing_geno_ptr)? nullptr : g_missing_geno_ptr;
  char* buf_start = nullptr;
  uintptr_t marker_idx;
  char* bufptr;
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_map_or_bim_ret_OPEN_FAIL;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    while (marker_uidx >= chrom_end) {
      chrom_idx = chrom_info_ptr->chrom_file_order[++chrom_fo_idx];
      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
      buf_start = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
      *buf_start++ = delim;
    }
    bufptr = strcpyax(buf_start, &(marker_ids[marker_uidx * max_marker_id_len]), delim);
    if (!marker_cms) {
      *bufptr++ = '0';
    } else {
      bufptr = dtoa_g_wxp8(marker_cms[marker_uidx], 1, bufptr);
    }
    *bufptr++ = delim;
    bufptr = uint32toa(marker_pos[marker_uidx], bufptr);
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
      goto write_map_or_bim_ret_WRITE_FAIL;
    }
    if (marker_allele_ptrs) {
      // quasi-bugfix: --missing-genotype + --output-missing-genotype did not
      // have the expected behavior with regular --make-bed (though it *did*
      // already do the right thing when sorting the .bim...)
      putc_unlocked(delim, outfile);
      if (!missing_geno_ptr) {
	fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
	putc_unlocked(delim, outfile);
	fputs(marker_allele_ptrs[2 * marker_uidx + 1], outfile);
      } else {
	fputs(cond_replace(marker_allele_ptrs[2 * marker_uidx], missing_geno_ptr, output_missing_geno_ptr), outfile);
	putc_unlocked(delim, outfile);
	fputs(cond_replace(marker_allele_ptrs[2 * marker_uidx + 1], missing_geno_ptr, output_missing_geno_ptr), outfile);
      }
    }
    putc_unlocked('\n', outfile);
  }
  if (fclose_null(&outfile)) {
    goto write_map_or_bim_ret_WRITE_FAIL;
  }
  while (0) {
  write_map_or_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_map_or_bim_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t load_bim_split_chrom(char* bimname, uintptr_t* marker_exclude, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, int64_t* ll_buf, uint32_t max_bim_linelen) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  char* loadbuf = g_textbuf;
  int32_t retval = 0;
  {
    if (max_bim_linelen > MAXLINELEN) {
      if (bigstack_alloc_c(max_bim_linelen, &loadbuf)) {
	goto load_bim_split_chrom_ret_NOMEM;
      }
    }
    if (fopen_checked(bimname, "r", &infile)) {
      goto load_bim_split_chrom_ret_OPEN_FAIL;
    }
    uint32_t marker_uidx = 0xffffffffU; // deliberate overflow
    for (uintptr_t marker_idx = 0; marker_idx < marker_ct;) {
      if (!fgets(loadbuf, max_bim_linelen, infile)) {
	goto load_bim_split_chrom_ret_READ_FAIL;
      }
      char* loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (is_eoln_or_comment_kns(*loadbuf_first_token)) {
	continue;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	continue;
      }
      // already validated
      char* first_token_end = token_endnn(loadbuf_first_token);
      const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - loadbuf_first_token);
      *first_token_end = '\0';
      uint64_t chrom_idx = (uint32_t)get_chrom_code(loadbuf_first_token, chrom_info_ptr, chrom_name_slen);
      ll_buf[marker_idx] = (int64_t)((chrom_idx << 32) | ((uint64_t)marker_idx));
      ++marker_idx;
    }
  }
  while (0) {
  load_bim_split_chrom_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_bim_split_chrom_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_bim_split_chrom_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

void fill_ll_buf(uintptr_t* marker_exclude, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, int64_t* ll_buf) {
  uint32_t marker_uidx = 0;
  uint32_t chrom_idx_p1 = 0;
  uint32_t chrom_end = 0;
  uint64_t chrom_idx_shifted = 0;
  uintptr_t marker_idx;
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (marker_uidx >= chrom_end) {
      do {
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[++chrom_idx_p1];
      } while (marker_uidx >= chrom_end);
      chrom_idx_shifted = ((uint64_t)(chrom_info_ptr->chrom_file_order[chrom_idx_p1 - 1])) << 32;
    }
    ll_buf[marker_idx] = (int64_t)(chrom_idx_shifted | ((uint64_t)marker_idx));
  }
}

int32_t update_marker_chroms(Two_col_params* update_chr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr, int64_t* ll_buf) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  char skipchar = update_chr->skipchar;
  uint32_t colid_first = (update_chr->colid < update_chr->colx);
  uint32_t marker_ctl = BITCT_TO_WORDCT(marker_ct);
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t marker_uidx = 0;
  uintptr_t line_idx = update_chr->skip;
  uintptr_t* already_seen;
  char* loadbuf;
  char* colid_ptr;
  char* colx_ptr;
  uint32_t* marker_id_htable;
  uint32_t* marker_uidx_to_idx;
  uintptr_t loadbuf_size;
  uint32_t marker_id_htable_size;
  uint32_t colmin;
  uint32_t coldiff;
  uint32_t slen;
  uint32_t marker_idx;
  int32_t retval;
  char cc;
  {
    retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, &marker_id_htable_size, &marker_id_htable);
    if (retval) {
      goto update_marker_chroms_ret_1;
    }
    if (bigstack_calloc_ul(marker_ctl, &already_seen) ||
	bigstack_alloc_ui(unfiltered_marker_ct, &marker_uidx_to_idx)) {
      goto update_marker_chroms_ret_NOMEM;
    }
    fill_uidx_to_idx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_uidx_to_idx);
    loadbuf = (char*)g_bigstack_base;
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    }
    if (loadbuf_size <= MAXLINELEN) {
      goto update_marker_chroms_ret_NOMEM;
    }
    retval = open_and_skip_first_lines(&infile, update_chr->fname, loadbuf, loadbuf_size, update_chr->skip);
    if (retval) {
      goto update_marker_chroms_ret_1;
    }
    if (colid_first) {
      colmin = update_chr->colid - 1;
      coldiff = update_chr->colx - update_chr->colid;
    } else {
      colmin = update_chr->colx - 1;
      coldiff = update_chr->colid - update_chr->colx;
    }
    while (fgets(loadbuf, loadbuf_size, infile)) {
      line_idx++;
      if (!(loadbuf[loadbuf_size - 1])) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-chr file is pathologically long.\n", line_idx);
	  goto update_marker_chroms_ret_INVALID_FORMAT_2;
	} else {
	  goto update_marker_chroms_ret_NOMEM;
	}
      }
      colid_ptr = skip_initial_spaces(loadbuf);
      cc = *colid_ptr;
      if (is_eoln_kns(cc) || (cc == skipchar)) {
	continue;
      }
      if (colid_first) {
	colid_ptr = next_token_multz(colid_ptr, colmin);
	colx_ptr = next_token_mult(colid_ptr, coldiff);
	if (no_more_tokens_kns(colx_ptr)) {
	  goto update_marker_chroms_ret_MISSING_TOKENS;
	}
      } else {
	colx_ptr = next_token_multz(colid_ptr, colmin);
	colid_ptr = next_token_mult(colx_ptr, coldiff);
	if (no_more_tokens_kns(colid_ptr)) {
	  goto update_marker_chroms_ret_MISSING_TOKENS;
	}
      }
      slen = strlen_se(colid_ptr);
      marker_uidx = id_htable_find(colid_ptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
      if (marker_uidx == 0xffffffffU) {
	miss_ct++;
	continue;
      }
      marker_idx = marker_uidx_to_idx[marker_uidx];
      if (is_set(already_seen, marker_idx)) {
	colid_ptr[slen] = '\0';
	LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --update-chr file.\n", colid_ptr);
	goto update_marker_chroms_ret_INVALID_FORMAT_2;
      }
      set_bit(marker_idx, already_seen);
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code_destructive("--update-chr file", line_idx, allow_extra_chroms, colx_ptr, token_endnn(colx_ptr), chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto update_marker_chroms_ret_1;
      }
      ll_buf[marker_idx] = (int64_t)((((uint64_t)((uint32_t)cur_chrom_code)) << 32) | (((uint64_t)ll_buf[marker_idx]) & 0xffffffffLLU));
      hit_ct++;
    }
    if (!feof(infile)) {
      goto update_marker_chroms_ret_READ_FAIL;
    }
    if (miss_ct) {
      sprintf(g_logbuf, "--update-chr: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      sprintf(g_logbuf, "--update-chr: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logprintb();
  }
  while (0) {
  update_marker_chroms_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_chroms_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_chroms_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-chr file has fewer tokens than expected.\n", line_idx);
  update_marker_chroms_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_chroms_ret_1:
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

void sort_marker_chrom_pos(int64_t* ll_buf, uintptr_t marker_ct, uint32_t* pos_buf, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t* unpack_map, uint32_t* chrom_ct_ptr) {
  // Assumes ll_buf is initially filled with chromosome idxs in high 32 bits,
  // and filtered marker indices in low 32 bits.  pos_buf is expected to have
  // base-pair positions; lookup is by filtered_index iff unpack_map is
  // nullptr.
  // After this is finished, ll_buf has marker positions in high bits and
  // filtered original indices in low bits, while chrom_start[] tracks
  // chromosome boundaries.
  uintptr_t marker_idx;
  uint32_t uii;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  if (!marker_ct) {
    chrom_start[0] = 0;
    *chrom_ct_ptr = 0;
    return;
  }
#ifdef __cplusplus
  std::sort(ll_buf, &(ll_buf[marker_ct]));
#else
  qsort(ll_buf, marker_ct, sizeof(int64_t), llcmp);
#endif
  cur_chrom = ll_buf[0] >> 32;
  chrom_ct = 0;
  chrom_start[0] = 0;
  chrom_id[0] = cur_chrom;
  uii = (uint32_t)ll_buf[0];
  if (unpack_map) {
    ll_buf[0] = ((uint64_t)uii) | (((uint64_t)pos_buf[unpack_map[uii]]) << 32);
    for (marker_idx = 1; marker_idx < marker_ct; marker_idx++) {
      if ((ll_buf[marker_idx] >> 32) != cur_chrom) {
	cur_chrom = ll_buf[marker_idx] >> 32;
	chrom_start[++chrom_ct] = marker_idx;
	chrom_id[chrom_ct] = cur_chrom;
      }
      uii = (uint32_t)ll_buf[marker_idx];
      ll_buf[marker_idx] = ((uint64_t)uii) | (((uint64_t)pos_buf[unpack_map[uii]]) << 32);
    }
  } else {
    ll_buf[0] = ((uint64_t)uii) | (((uint64_t)pos_buf[uii]) << 32);
    for (marker_idx = 1; marker_idx < marker_ct; marker_idx++) {
      if ((ll_buf[marker_idx] >> 32) != cur_chrom) {
	cur_chrom = ll_buf[marker_idx] >> 32;
	chrom_start[++chrom_ct] = marker_idx;
	chrom_id[chrom_ct] = cur_chrom;
      }
      uii = (uint32_t)ll_buf[marker_idx];
      ll_buf[marker_idx] = ((uint64_t)uii) | (((uint64_t)pos_buf[uii]) << 32);
    }
  }
  chrom_start[++chrom_ct] = marker_ct;
  for (uii = 0; uii < chrom_ct; uii++) {
#ifdef __cplusplus
    std::sort(&(ll_buf[chrom_start[uii]]), &(ll_buf[chrom_start[uii + 1]]));
#else
    qsort(&(ll_buf[chrom_start[uii]]), chrom_start[uii + 1] - chrom_start[uii], sizeof(int64_t), llcmp);
#endif
  }
  *chrom_ct_ptr = chrom_ct;
}

int32_t sort_and_write_bim(uint32_t* map_reverse, char* outname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, int64_t* ll_buf, Chrom_info* chrom_info_ptr) {
  // caller is expected to pop stuff off stack
  FILE* outfile = nullptr;
  uint32_t max_code = chrom_info_ptr->max_code;
  uint32_t chrom_code_end = max_code + 1 + chrom_info_ptr->name_ct;
  int32_t retval = 0;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  const char* output_missing_geno_ptr = g_output_missing_geno_ptr;
  uint32_t* chrom_start;
  uint32_t* chrom_id;
  uint32_t* unpack_map;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  char* chrom_name_end;
  char* bufptr;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  uint32_t uii;
  uint32_t ujj;
  // There can be a LOT of markers (1000 Genomes files can have ~40-80
  // million), so speeding up the sorting step over just calling qsort_ext()
  // may not be a complete waste of effort.
  // Strategy:
  // 1. fill ll_buf with chromosome idx in high-order bits, original position
  //    in low-order.
  // 2. std::sort() ll_buf, read off chromosome boundaries
  // 3. then replace high-order bits in ll_buf with marker positions, and
  //    std::sort() each chromosome separately.
  // Would be even faster if this was performed in a single sort, in the
  // super-common case where all three numbers can be squeezed together in 64
  // bits.  But we care most about performance when this can't be done, so I
  // haven't bothered with that optimization.
  if (bigstack_alloc_ui(chrom_code_end + 1, &chrom_start) ||
      bigstack_alloc_ui(chrom_code_end, &chrom_id) ||
      bigstack_alloc_ui(marker_ct, &unpack_map)) {
    goto sort_and_write_bim_ret_NOMEM;
  }
  fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, unpack_map);
  sort_marker_chrom_pos(ll_buf, marker_ct, marker_pos, chrom_start, chrom_id, unpack_map, &chrom_ct);
  if (fopen_checked(outname, "w", &outfile)) {
    goto sort_and_write_bim_ret_OPEN_FAIL;
  }

  marker_idx = 0;
  for (uii = 0; uii < chrom_ct; uii++) {
    cur_chrom = chrom_id[uii];
    ujj = chrom_start[uii + 1];
    chrom_name_end = chrom_name_write(chrom_info_ptr, cur_chrom, g_textbuf);
    *chrom_name_end++ = '\t';
    for (; marker_idx < ujj; marker_idx++) {
      marker_uidx = unpack_map[(uint32_t)ll_buf[marker_idx]];
      bufptr = strcpyax(chrom_name_end, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
      if (!marker_cms) {
	*bufptr++ = '0';
      } else {
        bufptr = dtoa_g_wxp8(marker_cms[marker_uidx], 1, bufptr);
      }
      *bufptr++ = '\t';
      bufptr = uint32toa_x((uint32_t)(ll_buf[marker_idx] >> 32), '\t', bufptr);
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto sort_and_write_bim_ret_WRITE_FAIL;
      }
      fputs(cond_replace(marker_allele_ptrs[2 * marker_uidx], missing_geno_ptr, output_missing_geno_ptr), outfile);
      putc_unlocked('\t', outfile);
      fputs(cond_replace(marker_allele_ptrs[2 * marker_uidx + 1], missing_geno_ptr, output_missing_geno_ptr), outfile);
      if (putc_checked('\n', outfile)) {
	goto sort_and_write_bim_ret_WRITE_FAIL;
      }
      map_reverse[marker_uidx] = marker_idx;
    }
  }
  while (0) {
  sort_and_write_bim_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sort_and_write_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sort_and_write_bim_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t load_sort_and_write_map(uint32_t** map_reverse_ptr, FILE* mapfile, uint32_t map_cols, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t max_marker_id_len, int32_t compact_map_reverse, Chrom_info* chrom_info_ptr) {
  // get_chrom_code() cannot fail
  FILE* map_outfile = nullptr;
  int64_t* ll_buf = nullptr;
  uintptr_t line_idx = 0;
  uint32_t orig_zec = chrom_info_ptr->zero_extra_chroms;
  int32_t retval = 0;
  char zstr[2];
  char* marker_ids;
  double* marker_cms;
  uint32_t* pos_buf;
  uint32_t* unpack_map;
  uint32_t* chrom_start;
  uint32_t* chrom_id;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  char* bufptr0;
  char* bufptr;
  uint32_t uii;
  uint32_t ujj;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  {
    // See sort_and_write_bim() for discussion.  Note that marker_ids and
    // marker_cms use filtered instead of unfiltered coordinates, though.
    if (bigstack_alloc_ui(compact_map_reverse? marker_ct : unfiltered_marker_ct, map_reverse_ptr) ||
	bigstack_alloc_ll(marker_ct, &ll_buf) ||
	bigstack_alloc_c(marker_ct * max_marker_id_len, &marker_ids) ||
	bigstack_alloc_d(marker_ct, &marker_cms) ||
	bigstack_alloc_ui(marker_ct, &pos_buf) ||
	bigstack_alloc_ui(marker_ct, &unpack_map) ||
	bigstack_alloc_ui(MAX_POSSIBLE_CHROM + 2, &chrom_start) ||
	bigstack_alloc_ui(MAX_POSSIBLE_CHROM + 1, &chrom_id)) {
      goto load_sort_and_write_map_ret_NOMEM;
    }
    rewind(mapfile);
    marker_idx = 0;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      char* textbuf_first_token;
      if (get_next_noncomment(mapfile, &textbuf_first_token, &line_idx)) {
	goto load_sort_and_write_map_ret_READ_FAIL;
      }
      if (IS_SET(marker_exclude, marker_uidx)) {
	continue;
      }
      char* first_token_end = token_endnn(textbuf_first_token);
      char* textbuf_iter = skip_initial_spaces(first_token_end);
      const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - textbuf_first_token);
      *first_token_end = '\0';
      ll_buf[marker_idx] = (((uint64_t)((uint32_t)get_chrom_code(textbuf_first_token, chrom_info_ptr, chrom_name_slen))) << 32) + marker_idx;
      uii = strlen_se(textbuf_iter);
      memcpyx(&(marker_ids[marker_idx * max_marker_id_len]), textbuf_iter, uii, 0);
      textbuf_iter = next_token(textbuf_iter);
      if (map_cols == 4) {
	if (scan_double(textbuf_iter, &(marker_cms[marker_idx]))) {
	  marker_cms[marker_idx] = 0.0;
	}
	textbuf_iter = next_token(textbuf_iter);
      } else {
	marker_cms[marker_idx] = 0.0;
      }
      unpack_map[marker_idx] = marker_uidx;
      // previously validated
      scan_uint_defcap(textbuf_iter, &(pos_buf[marker_idx++]));
    }
    sort_marker_chrom_pos(ll_buf, marker_ct, pos_buf, chrom_start, chrom_id, nullptr, &chrom_ct);

    strcpy(outname_end, ".map.tmp");
    if (fopen_checked(outname, "w", &map_outfile)) {
      goto load_sort_and_write_map_ret_OPEN_FAIL;
    }

    marker_idx = 0;
    *zstr = '0';
    zstr[1] = '\0';
    chrom_info_ptr->zero_extra_chroms = 0;
    for (uii = 0; uii < chrom_ct; uii++) {
      cur_chrom = chrom_id[uii];
      ujj = chrom_start[uii + 1];
      bufptr0 = chrom_name_write(chrom_info_ptr, cur_chrom, g_textbuf);
      *bufptr0++ = '\t';
      for (; marker_idx < ujj; marker_idx++) {
	marker_idx2 = (uint32_t)ll_buf[marker_idx];
	marker_uidx = unpack_map[marker_idx2];
	bufptr = strcpyax(bufptr0, &(marker_ids[marker_idx2 * max_marker_id_len]), '\t');
	bufptr = dtoa_g_wxp8x(marker_cms[marker_idx2], 1, '\t', bufptr);
	bufptr = uint32toa_x((uint32_t)(ll_buf[marker_idx] >> 32), '\n', bufptr);
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, map_outfile)) {
	  goto load_sort_and_write_map_ret_WRITE_FAIL;
	}
	(*map_reverse_ptr)[compact_map_reverse? marker_idx2 : marker_uidx] = marker_idx;
      }
    }
    if (fclose_null(&map_outfile)) {
      goto load_sort_and_write_map_ret_WRITE_FAIL;
    }
  }
  while (0) {
  load_sort_and_write_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_sort_and_write_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_sort_and_write_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_sort_and_write_map_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  chrom_info_ptr->zero_extra_chroms = orig_zec;
  if (ll_buf) {
    bigstack_reset(ll_buf);
  }
  return retval;
}

int32_t flip_subset_init(char* flip_fname, char* flip_subset_fname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t* sample_sort_map, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* flip_subset_markers, uintptr_t* flip_subset_vec2) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  uintptr_t miss_ct = 0;
  uintptr_t line_idx = 0;
  uint32_t* sample_uidx_to_idx = nullptr;
  uint32_t flip_marker_ct = 0;
  uint32_t flip_sample_ct = 0;
  int32_t retval = 0;
  const char reverse_complements[] = "T\0G\0\0\0C\0\0\0\0\0\0\0\0\0\0\0\0A";
  uint32_t* marker_id_htable;
  char* sorted_sample_ids;
  uint32_t* sample_id_map;
  char* bufptr;
  char* a1ptr;
  char* a2ptr;
  char* id_buf;
  uint32_t marker_id_htable_size;
  uint32_t slen;
  uint32_t marker_uidx;
  uint32_t sample_idx_write;
  int32_t sorted_idx;
  unsigned char ucc;
  // load --flip file, then --flip-subset
  fill_ulong_zero(unfiltered_marker_ctl, flip_subset_markers);
  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, &marker_id_htable_size, &marker_id_htable);
  if (retval) {
    goto flip_subset_init_ret_1;
  }
  if (fopen_checked(flip_fname, "r", &infile)) {
    goto flip_subset_init_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --flip file is pathologically long.\n", line_idx);
      goto flip_subset_init_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    slen = strlen_se(bufptr);
    marker_uidx = id_htable_find(bufptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx == 0xffffffffU) {
      continue;
    }
    a1ptr = marker_allele_ptrs[2 * marker_uidx];
    a2ptr = marker_allele_ptrs[2 * marker_uidx + 1];
    ucc = ((unsigned char)a1ptr[0]) - 'A';
    if (a1ptr[1] || a2ptr[1] || (ucc > 19) || (reverse_complements[ucc] != a2ptr[0])) {
      sprintf(g_logbuf, "Error: Invalid alleles (not reverse complement single bases) on line\n%" PRIuPTR " of --flip file.\n", line_idx);
      goto flip_subset_init_ret_INVALID_FORMAT_2;
    }
    if (is_set(flip_subset_markers, marker_uidx)) {
      bufptr[slen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate marker ID '%s' in --flip file.\n", bufptr);
      goto flip_subset_init_ret_INVALID_FORMAT_2;
    }
    set_bit(marker_uidx, flip_subset_markers);
    flip_marker_ct++;
  }
  if (fclose_null(&infile)) {
    goto flip_subset_init_ret_READ_FAIL;
  }
  bigstack_reset(bigstack_mark);
  retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, unfiltered_sample_ct - sample_ct, sample_ids, max_sample_id_len, 0, 1, strcmp_deref, &sorted_sample_ids, &sample_id_map);
  if (retval) {
    goto flip_subset_init_ret_1;
  }
  if (bigstack_alloc_c(max_sample_id_len, &id_buf)) {
    goto flip_subset_init_ret_NOMEM;
  }
  if (sample_sort_map) {
    if (bigstack_alloc_ui(unfiltered_sample_ct, &sample_uidx_to_idx)) {
      goto flip_subset_init_ret_NOMEM;
    }
    fill_uidx_to_idx(sample_exclude, unfiltered_sample_ct, sample_ct, sample_uidx_to_idx);
  }
  if (fopen_checked(flip_subset_fname, "r", &infile)) {
    goto flip_subset_init_ret_OPEN_FAIL;
  }
  fill_ulong_zero(sample_ctv2, flip_subset_vec2);
  line_idx = 0;
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --flip-subset file is pathologically long.\n", line_idx);
      goto flip_subset_init_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, nullptr, &sorted_idx, id_buf) || (sorted_idx == -1)) {
      miss_ct++;
      continue;
    }
    if (!sample_sort_map) {
      sample_idx_write = sample_id_map[(uint32_t)sorted_idx];
    } else {
      sample_idx_write = sample_uidx_to_idx[sample_sort_map[sample_id_map[(uint32_t)sorted_idx]]];
    }
    if (IS_SET_DBL(flip_subset_vec2, sample_idx_write)) {
      *strchr(id_buf, '\t') = ' ';
      LOGPREPRINTFWW("Error: Duplicate sample ID '%s' in --flip-subset file.\n", id_buf);
      goto flip_subset_init_ret_INVALID_FORMAT_2;
    }
    SET_BIT_DBL(sample_idx_write, flip_subset_vec2);
    flip_sample_ct++;
  }
  if (fclose_null(&infile)) {
    goto flip_subset_init_ret_READ_FAIL;
  }
  LOGPRINTF("--flip-subset: Flipping %u SNP%s for %u %s.\n", flip_marker_ct, (flip_marker_ct == 1)? "" : "s", flip_sample_ct, species_str(flip_sample_ct));
  if (miss_ct) {
    LOGERRPRINTF("Warning: %" PRIuPTR " --flip-subset line%s skipped.\n", miss_ct, (miss_ct == 1)? "" : "s");
  }
  while (0) {
  flip_subset_init_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  flip_subset_init_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  flip_subset_init_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  flip_subset_init_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 flip_subset_init_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(infile);
  return retval;
}

uint32_t merge_or_split_x(uint32_t mergex, uint32_t splitx_bound1, uint32_t splitx_bound2, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, int64_t* ll_buf) {
  // If --merge-x, splitx_bound1 and splitx_bound2 are 0; turns out those
  // settings let --merge-x and --split-x use the same logic.
  // ll_buf[] has chromosome codes in high 32 bits and filtered indices in low
  // 32 bits.
  uint32_t hit_ct = 0;
  uint32_t marker_uidx = 0;
  uint32_t match_chrom;
  uint64_t new_chrom_shifted;
  uint32_t marker_idx;
  uint32_t cur_pos;
  if (mergex) {
    match_chrom = chrom_info_ptr->xymt_codes[XY_OFFSET];
    new_chrom_shifted = ((uint64_t)((uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET])) << 32;
  } else {
    match_chrom = chrom_info_ptr->xymt_codes[X_OFFSET];
    new_chrom_shifted = ((uint64_t)((uint32_t)chrom_info_ptr->xymt_codes[XY_OFFSET])) << 32;
  }

  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++, marker_uidx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if ((((uint64_t)ll_buf[marker_idx]) >> 32) == match_chrom) {
      cur_pos = marker_pos[marker_uidx];
      if ((cur_pos >= splitx_bound2) || (cur_pos <= splitx_bound1)) {
        ll_buf[marker_idx] = (int64_t)(new_chrom_shifted | ((uint64_t)marker_idx));
        hit_ct++;
      }
    }
  }
  if (!hit_ct) {
    return 1;
  }
  LOGPRINTF("--%s: %u chromosome code%s changed.\n", mergex? "merge-x" : "split-x", hit_ct, (hit_ct == 1)? "" : "s");
  return 0;
}

int32_t make_bed_one_marker(FILE* bedfile, uintptr_t* loadbuf, uint32_t unfiltered_sample_ct, uintptr_t unfiltered_sample_ct4, uintptr_t* sample_exclude, uint32_t sample_ct, uint32_t* sample_sort_map, uintptr_t final_mask, uint32_t is_reverse, uintptr_t* writebuf) {
  // final_mask only relevant if unfiltered_sample_ct == sample_ct
  uintptr_t* writeptr = writebuf;
  uintptr_t cur_word = 0;
  uint32_t sample_uidx = 0;
  uint32_t ii_rem = 0;
  uint32_t sample_idx = 0;
  uint32_t sample_uidx2;
  if (sample_sort_map) {
    if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf)) {
      return RET_READ_FAIL;
    }
    for (; sample_idx < sample_ct; sample_idx++) {
      do {
	sample_uidx2 = sample_sort_map[sample_uidx++];
      } while (IS_SET(sample_exclude, sample_uidx2));
      cur_word |= EXTRACT_2BIT_GENO(loadbuf, sample_uidx2) << (ii_rem * 2);
      if (++ii_rem == BITCT2) {
	*writeptr++ = cur_word;
	cur_word = 0;
	ii_rem = 0;
      }
    }
    if (ii_rem) {
      *writeptr = cur_word;
    }
    if (is_reverse) {
      reverse_loadbuf(sample_ct, (unsigned char*)writebuf);
    }
  } else {
    if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, is_reverse, bedfile, loadbuf, writeptr)) {
      return RET_READ_FAIL;
    }
  }
  return 0;
}

int32_t make_bed_me_missing_one_marker(FILE* bedfile, uintptr_t* loadbuf, uint32_t unfiltered_sample_ct, uintptr_t unfiltered_sample_ct4, uintptr_t* sample_exclude, uint32_t sample_ct, uint32_t* sample_sort_map, uintptr_t final_mask, uint32_t unfiltered_sample_ctl2m1, uint32_t is_reverse, uintptr_t* writebuf, uintptr_t* workbuf, uintptr_t* sex_male, uintptr_t* sample_raw_male_include2, uint32_t* trio_lookup, uint32_t trio_ct, uint32_t set_hh_missing, uint32_t is_x, uint32_t multigen, uint64_t* error_ct_ptr) {
  // requires final_mask to be for unfiltered_sample_ct
  uintptr_t* writeptr = writebuf;
  uintptr_t cur_word = 0;
  uint32_t sample_uidx = 0;
  uint32_t ii_rem = 0;
  uint32_t sample_idx = 0;
  uint32_t sample_uidx2;
  if ((!sample_sort_map) && (unfiltered_sample_ct == sample_ct)) {
    loadbuf = writebuf;
  }
  if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf)) {
    return RET_READ_FAIL;
  }
  if (set_hh_missing && is_x) {
    hh_reset((unsigned char*)loadbuf, sample_raw_male_include2, unfiltered_sample_ct);
  }
  if (is_reverse) {
    reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf);
  }
  *error_ct_ptr += erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, sex_male, trio_lookup, trio_ct, is_x, multigen);
  if (sample_sort_map) {
    for (; sample_idx < sample_ct; sample_idx++) {
      do {
	sample_uidx2 = sample_sort_map[sample_uidx++];
      } while (IS_SET(sample_exclude, sample_uidx2));
      cur_word |= EXTRACT_2BIT_GENO(loadbuf, sample_uidx2) << (ii_rem * 2);
      if (++ii_rem == BITCT2) {
	*writeptr++ = cur_word;
	cur_word = 0;
	ii_rem = 0;
      }
    }
    if (ii_rem) {
      *writeptr = cur_word;
    }
  } else if (unfiltered_sample_ct != sample_ct) {
    copy_quaterarr_nonempty_subset_excl(loadbuf, sample_exclude, unfiltered_sample_ct, sample_ct, writebuf);
  }
  return 0;
}

void zeropatch(uintptr_t sample_ctv2, uintptr_t cluster_ct, uintptr_t* cluster_zc_masks, uint32_t** zcdefs, uintptr_t* patchbuf, uintptr_t marker_idx, uintptr_t* writebuf) {
#ifdef __LP64__
  __m128i* writevec = (__m128i*)writebuf;
  __m128i* patchvec = (__m128i*)patchbuf;
  __m128i* patchvec_end = (__m128i*)(&(patchbuf[sample_ctv2]));
  __m128i vec1;
  __m128i vec2;
#else
  uintptr_t* patchbuf_end = &(patchbuf[sample_ctv2]);
  uintptr_t ulii;
#endif
  uintptr_t cluster_idx;
  uint32_t at_least_one_cluster = 0;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    if (in_setdef(zcdefs[cluster_idx], marker_idx)) {
      if (!at_least_one_cluster) {
	at_least_one_cluster = 1;
	fill_ulong_zero(sample_ctv2, patchbuf);
      }
      bitvec_or(&(cluster_zc_masks[cluster_idx * sample_ctv2]), sample_ctv2, patchbuf);
    }
  }
  if (!at_least_one_cluster) {
    return;
  }
#ifdef __LP64__
  do {
    vec1 = *writevec;
    vec2 = *patchvec++;
    vec1 = _mm_andnot_si128(_mm_slli_epi64(vec2, 1), vec1);
    *writevec = _mm_or_si128(vec1, vec2);
    writevec++;
  } while (patchvec < patchvec_end);
#else
  do {
    ulii = *patchbuf++;
    *writebuf = ((*writebuf) & (~(ulii << 1))) | ulii;
    writebuf++;
  } while (patchbuf < patchbuf_end);
#endif
}

void reverse_subset(uintptr_t* writebuf, uintptr_t* subset_vec2, uintptr_t word_ct) {
  // reverse_loadbuf() variant that requires subset_vec2 bit to be set
#ifdef __LP64__
  __m128i* wvec = (__m128i*)writebuf;
  __m128i* svec = (__m128i*)subset_vec2;
  __m128i* wvec_end = (__m128i*)(&(writebuf[word_ct]));
  __m128i vii;
  __m128i vjj;
  while (wvec < wvec_end) {
    vii = *wvec;
    vjj = _mm_andnot_si128(_mm_xor_si128(vii, _mm_srli_epi64(vii, 1)), *svec++);
    vjj = _mm_or_si128(vjj, _mm_slli_epi64(vjj, 1));
    *wvec++ = _mm_xor_si128(vii, vjj);
  }
#else
  uintptr_t* writebuf_end = &(writebuf[word_ct]);
  uintptr_t ulii;
  uintptr_t uljj;
  while (writebuf < writebuf_end) {
    ulii = *writebuf;
    uljj = (*subset_vec2++) & (~(ulii ^ (ulii >> 1)));
    uljj *= 3;
    *writebuf++ = ulii ^ uljj;
  }
#endif
}

void replace_missing_a2(uintptr_t* writebuf, uintptr_t* subset_vec2, uintptr_t word_ct) {
  // 01 -> 11 for each set bit in subset_vec2
#ifdef __LP64__
  __m128i* wvec = (__m128i*)writebuf;
  __m128i* svec = (__m128i*)subset_vec2;
  __m128i* wvec_end = (__m128i*)(&(writebuf[word_ct]));
  __m128i vii;
  __m128i vjj;
  do {
    vii = *wvec;
    vjj = _mm_andnot_si128(vii, _mm_slli_epi64(_mm_and_si128(vii, *svec++), 1));
    *wvec++ = _mm_or_si128(vii, vjj);
  } while (wvec < wvec_end);
#else
  uintptr_t* writebuf_end = &(writebuf[word_ct]);
  uintptr_t ulii;
  uintptr_t uljj;
  do {
    ulii = *writebuf;
    uljj = (~ulii) & (((*subset_vec2++) & ulii) << 1);
    *writebuf++ = ulii | uljj;
  } while (writebuf < writebuf_end);
#endif
}

int32_t make_bed(FILE* bedfile, uintptr_t bed_offset, char* bimname, char* outname, char* outname_end, uint64_t calculation_type, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint32_t* sample_sort_map, uint64_t misc_flags, uint32_t splitx_bound1, uint32_t splitx_bound2, Two_col_params* update_chr, char* flip_fname, char* flip_subset_fname, char* zerofname, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists, Chrom_info* chrom_info_ptr, uint32_t mendel_modifier, uint32_t max_bim_linelen) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t sample_ct4 = (sample_ct + 3) / 4;
  uintptr_t sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uintptr_t trio_ct = 0;
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  FILE* bedoutfile = nullptr;
  int64_t* ll_buf = nullptr;
  uintptr_t* sample_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uintptr_t* sample_raw_male_include2 = nullptr;
  uintptr_t* workbuf = nullptr;
  uintptr_t* cluster_zc_masks = nullptr;
  uintptr_t* patchbuf = nullptr;
  uintptr_t* flip_subset_markers = nullptr;
  uintptr_t* flip_subset_vec2 = nullptr;
  uint64_t* family_list = nullptr;
  uint64_t* trio_list = nullptr;
  uint32_t* trio_lookup = nullptr;
  uint32_t** zcdefs = nullptr;
  uint64_t mendel_error_ct = 0;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t family_ct = 0;
  uint32_t set_hh_missing = (misc_flags / MISC_SET_HH_MISSING) & 1;

  // todo: clear this when no MT markers loaded
  uint32_t set_mixed_mt_missing = (misc_flags / MISC_SET_MIXED_MT_MISSING) & 1;

  uint32_t set_me_missing = ((misc_flags / MISC_SET_ME_MISSING) & 1) && sample_ct;
  uint32_t fill_missing_a2 = (misc_flags / MISC_FILL_MISSING_A2) & 1;
  uint32_t mendel_include_duos = (mendel_modifier / MENDEL_DUOS) & 1;
  uint32_t mendel_multigen = (mendel_modifier / MENDEL_MULTIGEN) & 1;
  uint32_t mergex = (misc_flags / MISC_MERGEX) & 1;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t resort_map = map_is_unsorted || mergex || splitx_bound2 || update_chr;
  uint32_t chrom_end = 0;
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t pct = 1;
  int32_t retval = 0;
  const char errflags[][24] = {"set-hh-missing", "set-mixed-mt-missing", "set-me-missing", "fill-missing-a2"};
  uintptr_t* loadbuf;
  uintptr_t* writebuf;
  uintptr_t* writebuf_ptr;
  uint32_t* map_reverse;
  const char* errptr;
  uintptr_t cur_bigstack_left;
  uintptr_t pass_size;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t loop_end;
  uint32_t pass_ct;
  uint32_t pass_idx;
  uint32_t pass_start;
  uint32_t pass_end;
  uint32_t seek_needed;
  uint32_t markers_done;
  if (flip_subset_fname) {
    if (bigstack_alloc_ul(unfiltered_marker_ctl, &flip_subset_markers) ||
        bigstack_alloc_ul(sample_ctv2, &flip_subset_vec2)) {
      goto make_bed_ret_NOMEM;
    }
    retval = flip_subset_init(flip_fname, flip_subset_fname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, unfiltered_sample_ct, sample_exclude, sample_ct, sample_sort_map, sample_ids, max_sample_id_len, flip_subset_markers, flip_subset_vec2);
    if (retval) {
      goto make_bed_ret_1;
    }
  }
  if (calculation_type & CALC_MAKE_BED) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf)) {
      goto make_bed_ret_NOMEM;
    }

    if (zerofname && cluster_ct) {
      zcdefs = (uint32_t**)bigstack_alloc(cluster_ct * sizeof(intptr_t));
      if (!zcdefs) {
	goto make_bed_ret_NOMEM;
      }
      retval = zero_cluster_init(zerofname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, unfiltered_sample_ct, sample_exclude, sample_ct, sample_sort_map, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, zcdefs, &cluster_zc_masks);
      if (retval) {
	goto make_bed_ret_1;
      }
      if (bigstack_alloc_ul(sample_ctv2, &patchbuf)) {
	goto make_bed_ret_NOMEM;
      }
    }
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(outname, FOPEN_WB, &bedoutfile)) {
      goto make_bed_ret_OPEN_FAIL;
    }

    if (fwrite_checked("l\x1b\x01", 3, bedoutfile)) {
      goto make_bed_ret_WRITE_FAIL;
    }
    fflush(stdout);
    if (resort_map) {
      if (set_hh_missing || set_mixed_mt_missing || set_me_missing || fill_missing_a2) {
	// could remove this restriction if we added a chromosome check to the
	// main loop, but no real need to bother.
	errptr = errflags[set_hh_missing? 0 : (set_mixed_mt_missing? 1 : (set_me_missing? 2 : 3))];
	if (map_is_unsorted) {
	  // assumes UNSORTED_CM was masked out
	  LOGERRPRINTF("Error: --%s cannot be used on an unsorted .bim file.  Use\n--make-bed without --%s to sort by position first; then run\n--make-bed + --%s on the new fileset.\n", errptr, errptr, errptr);
	} else {
	  LOGERRPRINTF("Error: --%s cannot be used with --merge-x/--split-x/--update-chr.\nFinish updating chromosome codes first.\n", errptr);
	}
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto make_bed_ret_1;
      }
      if (bigstack_alloc_ui(unfiltered_marker_ct, &map_reverse) ||
	  bigstack_alloc_ll(marker_ct, &ll_buf)) {
	goto make_bed_ret_NOMEM;
      }
      if ((map_is_unsorted & UNSORTED_SPLIT_CHROM) || mergex || splitx_bound2 || update_chr) {
	if (map_is_unsorted & UNSORTED_SPLIT_CHROM) {
	  retval = load_bim_split_chrom(bimname, marker_exclude, marker_ct, chrom_info_ptr, ll_buf, max_bim_linelen);
	  if (retval) {
	    goto make_bed_ret_1;
	  }
	} else {
	  fill_ll_buf(marker_exclude, marker_ct, chrom_info_ptr, ll_buf);
	}
	if (update_chr) {
	  retval = update_marker_chroms(update_chr, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, allow_extra_chroms, chrom_info_ptr, ll_buf);
	  if (retval) {
	    goto make_bed_ret_1;
	  }
	} else if (mergex || splitx_bound2) {
	  if (splitx_bound2 && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[XY_OFFSET])) {
	    logerrprint("Error: --split-x cannot be used when the dataset already contains an XY region.\n");
	    goto make_bed_ret_INVALID_CMDLINE;
	  }
	  if (merge_or_split_x(mergex, splitx_bound1, splitx_bound2, unfiltered_marker_ct, marker_exclude, marker_ct, marker_pos, chrom_info_ptr, ll_buf)) {
	    if (!(misc_flags & MISC_SPLIT_MERGE_NOFAIL)) {
	      if (mergex) {
		logerrprint("Error: --merge-x requires XY pseudo-autosomal region data.  (Use 'no-fail' to\nforce --make-bed to proceed anyway.\n");
	      } else {
		if (!is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[X_OFFSET])) {
		  logerrprint("Error: --split-x requires X chromosome data.  (Use 'no-fail' to force\n--make-bed to proceed anyway.\n");
		} else {
		  LOGERRPRINTFWW("Error: No X chromosome loci have bp positions <= %u or >= %u. (Use 'no-fail' to force --make-bed to proceed anyway.)\n", splitx_bound1, splitx_bound2);
		}
	      }
	      goto make_bed_ret_INVALID_CMDLINE;
	    }
	  }
	}
      } else {
	fill_ll_buf(marker_exclude, marker_ct, chrom_info_ptr, ll_buf);
      }
      memcpy(outname_end, ".bim", 5);
      retval = sort_and_write_bim(map_reverse, outname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, ll_buf, chrom_info_ptr);
      if (retval) {
	goto make_bed_ret_1;
      }
      bigstack_reset(ll_buf);

      // oops, forgot to multiply by sizeof(intptr_t)!  fortunately, this
      // segfaulted instead of corrupting any data.
      // anyway, it's now time to implement multipass.
      cur_bigstack_left = bigstack_left();
      if (cur_bigstack_left < sample_ctv2 * sizeof(intptr_t)) {
        goto make_bed_ret_NOMEM;
      }
      writebuf = (uintptr_t*)g_bigstack_base;
      // 32-bit bugfix (6 Sep 2017): necessary to cast numerator to 64-bit
      pass_ct = 1 + ((sample_ctv2 * ((uint64_t)marker_ct) * sizeof(intptr_t) - 1) / cur_bigstack_left);
      pass_size = 1 + ((marker_ct - 1) / pass_ct);
      *outname_end = '\0';
      LOGPRINTFWW5("--make-bed to %s.bed + %s.bim + %s.fam ... ", outname, outname, outname);
      fputs("0%", stdout);
      fflush(stdout);
      if (sample_ct) {
	loop_end = marker_ct / 100;
	markers_done = 0;
	for (pass_idx = 0; pass_idx < pass_ct; pass_idx++) {
	  pass_start = pass_idx * pass_size;
	  pass_end = (pass_idx + 1) * pass_size;
	  if (pass_idx + 1 == pass_ct) {
	    pass_end = marker_ct;
	  }
	  seek_needed = 1;
	  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	    if (IS_SET(marker_exclude, marker_uidx)) {
	      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	      seek_needed = 1;
	    }
	    if ((map_reverse[marker_uidx] < pass_start) || (map_reverse[marker_uidx] >= pass_end)) {
	      seek_needed = 1;
	      continue;
	    }
	    writebuf_ptr = &(writebuf[sample_ctv2 * (map_reverse[marker_uidx] - pass_start)]);
	    if (seek_needed) {
	      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
		goto make_bed_ret_READ_FAIL;
	      }
	      seek_needed = 0;
	    }
	    retval = make_bed_one_marker(bedfile, loadbuf, unfiltered_sample_ct, unfiltered_sample_ct4, sample_exclude, sample_ct, sample_sort_map, final_mask, IS_SET(marker_reverse, marker_uidx), writebuf_ptr);
	    if (retval) {
	      goto make_bed_ret_1;
	    }
	    if (zcdefs) {
	      zeropatch(sample_ctv2, cluster_ct, cluster_zc_masks, zcdefs, patchbuf, marker_idx, writebuf_ptr);
	    }
	    if (flip_subset_markers && is_set(flip_subset_markers, marker_uidx)) {
	      reverse_subset(writebuf_ptr, flip_subset_vec2, sample_ctv2);
	    }
	    if (markers_done >= loop_end) {
	      if (pct > 10) {
		putc_unlocked('\b', stdout);
	      }
	      pct = (markers_done * 100LLU) / marker_ct;
	      printf("\b\b%u%%", pct);
	      fflush(stdout);
	      pct++;
	      loop_end = (pct * ((uint64_t)marker_ct)) / 100;
	    }
	    markers_done++;
	  }
	  writebuf_ptr = writebuf;
	  for (marker_idx = pass_start; marker_idx < pass_end; marker_idx++) {
	    if (fwrite_checked(writebuf_ptr, sample_ct4, bedoutfile)) {
	      goto make_bed_ret_WRITE_FAIL;
	    }
	    writebuf_ptr = &(writebuf_ptr[sample_ctv2]);
	  }
	}
      }
    } else {
      if (!hh_exists) {
	set_hh_missing = 0;
      } else if (!set_hh_missing) {
	hh_exists = 0;
      }
      if (fill_missing_a2 || (hh_exists & NXMHH_EXISTS) || set_mixed_mt_missing) {
        if (bigstack_alloc_ul(sample_ctv2, &sample_include2)) {
	  goto make_bed_ret_NOMEM;
	}
        fill_quatervec_55(sample_ct, sample_include2);
      }
      if (alloc_collapsed_haploid_filters(sample_exclude, sex_male, unfiltered_sample_ct, sample_ct, fill_missing_a2? Y_FIX_NEEDED : hh_exists, 0, &sample_include2, &sample_male_include2)) {
	goto make_bed_ret_NOMEM;
      }
      if (set_me_missing) {
	retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, nullptr, nullptr, nullptr, nullptr, &family_list, &family_ct, &trio_list, &trio_ct, &trio_lookup, mendel_include_duos, mendel_multigen);
	if (retval) {
	  goto make_bed_ret_1;
	}
	if (trio_ct) {
	  if (bigstack_alloc_ul(unfiltered_sample_ctp1l2, &workbuf)) {
	    goto make_bed_ret_NOMEM;
	  }
	  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
	  if (set_hh_missing) {
	    if (bigstack_alloc_ul(unfiltered_sample_ctl2, &sample_raw_male_include2)) {
	      goto make_bed_ret_NOMEM;
	    }
	    init_quaterarr_from_inverted_bitarr(sex_male, unfiltered_sample_ct, sample_raw_male_include2);
	  }
	} else {
	  set_me_missing = 0;
	}
      }

      if (bigstack_alloc_ul(sample_ctv2, &writebuf)) {
	goto make_bed_ret_NOMEM;
      }
      if (fseeko(bedfile, bed_offset, SEEK_SET)) {
	goto make_bed_ret_READ_FAIL;
      }
      *outname_end = '\0';
      LOGPRINTFWW5("--make-bed to %s.bed + %s.bim + %s.fam ... ", outname, outname, outname);
      fputs("0%", stdout);
      if (sample_ct) {
	marker_uidx = 0;
	for (pct = 1; pct <= 100; pct++) {
	  loop_end = (pct * ((uint64_t)marker_ct)) / 100;
	  for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	    if (IS_SET(marker_exclude, marker_uidx)) {
	      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
		goto make_bed_ret_READ_FAIL;
	      }
	    }
	    if (marker_uidx >= chrom_end) {
	      chrom_fo_idx++;
	      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	    }
	    if ((!set_me_missing) || (is_haploid && (!is_x)) || is_mt) {
	      retval = make_bed_one_marker(bedfile, loadbuf, unfiltered_sample_ct, unfiltered_sample_ct4, sample_exclude, sample_ct, sample_sort_map, final_mask, IS_SET(marker_reverse, marker_uidx), writebuf);
	      if (is_haploid && set_hh_missing) {
		haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)writebuf);
	      } else if (is_mt && set_mixed_mt_missing) {
		hh_reset((unsigned char*)writebuf, sample_include2, sample_ct);
	      }
	    } else {
	      retval = make_bed_me_missing_one_marker(bedfile, loadbuf, unfiltered_sample_ct, unfiltered_sample_ct4, sample_exclude, sample_ct, sample_sort_map, final_mask, unfiltered_sample_ctl2m1, IS_SET(marker_reverse, marker_uidx), writebuf, workbuf, sex_male, sample_raw_male_include2, trio_lookup, trio_ct, set_hh_missing, is_x, mendel_multigen, &mendel_error_ct);
	    }
	    if (retval) {
	      goto make_bed_ret_1;
	    }
	    if (zcdefs) {
	      zeropatch(sample_ctv2, cluster_ct, cluster_zc_masks, zcdefs, patchbuf, marker_idx, writebuf);
	    }
	    if (flip_subset_markers && is_set(flip_subset_markers, marker_uidx)) {
	      reverse_subset(writebuf, flip_subset_vec2, sample_ctv2);
	    }
	    if (fill_missing_a2) {
	      replace_missing_a2(writebuf, is_y? sample_male_include2 : sample_include2, sample_ctv2);
	    }

	    if (fwrite_checked(writebuf, sample_ct4, bedoutfile)) {
	      goto make_bed_ret_WRITE_FAIL;
	    }
	  }
	  if (pct < 100) {
	    if (pct > 10) {
	      putc_unlocked('\b', stdout);
	    }
	    printf("\b\b%u%%", pct);
	    fflush(stdout);
	  }
	}
      }
    }
    if (fclose_null(&bedoutfile)) {
      goto make_bed_ret_WRITE_FAIL;
    }
  } else if (resort_map && (calculation_type & CALC_MAKE_BIM)) {
    memcpy(outname_end, ".bim", 5);
    if (calculation_type & CALC_MAKE_BIM) {
      LOGPRINTFWW5("--make-just-bim to %s ... ", outname);
      fflush(stdout);
    }
    if (bigstack_alloc_ui(unfiltered_marker_ct, &map_reverse) ||
	bigstack_alloc_ll(marker_ct, &ll_buf)) {
      goto make_bed_ret_NOMEM;
    }
    if (map_is_unsorted & UNSORTED_SPLIT_CHROM) {
      retval = load_bim_split_chrom(bimname, marker_exclude, marker_ct, chrom_info_ptr, ll_buf, max_bim_linelen);
      if (retval) {
	goto make_bed_ret_1;
      }
    } else {
      fill_ll_buf(marker_exclude, marker_ct, chrom_info_ptr, ll_buf);
    }
    retval = sort_and_write_bim(map_reverse, outname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, ll_buf, chrom_info_ptr);
    if (retval) {
      goto make_bed_ret_1;
    }
    bigstack_reset(map_reverse);
    if (calculation_type & CALC_MAKE_BIM) {
      logprint("done.\n");
    }
  }

  if (calculation_type & (CALC_MAKE_BED | CALC_MAKE_FAM)) {
    memcpy(outname_end, ".fam", 5);
    if (calculation_type & CALC_MAKE_FAM) {
      LOGPRINTFWW5("--make-just-fam to %s ... ", outname);
      fflush(stdout);
    }
    retval = write_fam(outname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno, ' ', sample_sort_map);
    if (retval) {
      goto make_bed_ret_1;
    }
    if (calculation_type & CALC_MAKE_FAM) {
      logprint("done.\n");
    }
  }

  if ((!resort_map) && (calculation_type & (CALC_MAKE_BED | CALC_MAKE_BIM))) {
    memcpy(outname_end, ".bim", 5);
    if (calculation_type & CALC_MAKE_BIM) {
      LOGPRINTFWW5("--make-just-bim to %s ... ", outname);
      fflush(stdout);
    }
    retval = write_map_or_bim(outname, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, '\t', chrom_info_ptr);
    if (retval) {
      goto make_bed_ret_1;
    }
    if (calculation_type & CALC_MAKE_BIM) {
      logprint("done.\n");
    }
  }

  if (calculation_type & CALC_MAKE_BED) {
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (set_me_missing) {
      LOGPRINTF("--set-me-missing: %" PRIu64 " error%s addressed.\n", mendel_error_ct, (mendel_error_ct == 1)? "" : "s");
    }
  }
  while (0) {
  make_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  make_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  make_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  make_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  make_bed_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 make_bed_ret_1:
  fclose_cond(bedoutfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t load_fam(char* famname, uint32_t fam_cols, uint32_t tmp_fam_col_6, int32_t missing_pheno, uint32_t affection_01, uintptr_t* unfiltered_sample_ct_ptr, char** sample_ids_ptr, uintptr_t* max_sample_id_len_ptr, char** paternal_ids_ptr, uintptr_t* max_paternal_id_len_ptr, char** maternal_ids_ptr, uintptr_t* max_maternal_id_len_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, uint32_t* affection_ptr, uintptr_t** pheno_nm_ptr, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, uintptr_t** founder_info_ptr, uintptr_t** sample_exclude_ptr, uint32_t allow_no_samples) {
  unsigned char* bigstack_mark = g_bigstack_base;
  double missing_phenod = (double)missing_pheno;
  uintptr_t* pheno_c = nullptr;
  double* pheno_d = nullptr;
  FILE* famfile = nullptr;
  uintptr_t unfiltered_sample_ct = 0;
  uintptr_t max_sample_id_len = *max_sample_id_len_ptr;
  uintptr_t max_paternal_id_len = *max_paternal_id_len_ptr;
  uintptr_t max_maternal_id_len = *max_maternal_id_len_ptr;
  uintptr_t line_idx = 0;
  uint32_t affection = 1;
  int32_t retval = 0;
  double pheno_ctrld = (double)((int32_t)(1 - affection_01));
  double pheno_cased = pheno_ctrld + 1.0;
  uintptr_t* sex_nm;
  uintptr_t* sex_male;
  uintptr_t* pheno_nm;
  uintptr_t* founder_info;
  char* bufptr0;
  char* bufptr;
  char* loadbuf;
  char* sample_ids;
  char* paternal_ids;
  char* maternal_ids;
  char cc;
  uintptr_t loadbuf_size;
  uintptr_t tmp_len;
  uintptr_t tmp_len2;
  uintptr_t sample_uidx;
  uintptr_t unfiltered_sample_ctl;
  double dxx;

  // we want this to work when the file is actually a .ped
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_fam_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  if (fopen_checked(famname, "r", &famfile)) {
    goto load_fam_ret_OPEN_FAIL;
  }
  // ----- .fam read, first pass -----
  // count number of samples, determine maximum person/father/mother ID
  // lengths, affection status, verify all floating point phenotype values are
  // valid
  while (fgets(loadbuf, loadbuf_size, famfile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .fam file is pathologically long.\n", line_idx);
	goto load_fam_ret_INVALID_FORMAT_2;
      } else {
	goto load_fam_ret_NOMEM;
      }
    }
    bufptr0 = skip_initial_spaces(loadbuf);
    if (!is_eoln_or_comment_kns(*bufptr0)) {
      if (fam_cols & FAM_COL_1) {
	bufptr = next_token(bufptr0);
	if (!bufptr) {
	  goto load_fam_ret_MISSING_TOKENS;
	}
      } else {
	bufptr = bufptr0;
      }
      tmp_len = strlen_se(bufptr);
      if ((tmp_len == 1) && (*bufptr == '0')) {
	sprintf(g_logbuf, "Error: Invalid IID '0' on line %" PRIuPTR " of .fam file.\n", line_idx);
	goto load_fam_ret_INVALID_FORMAT_2;
      }
      tmp_len = strlen_se(bufptr0) + strlen_se(bufptr) + 2;
      if (tmp_len > max_sample_id_len) {
	max_sample_id_len = tmp_len;
      }
      if (fam_cols & FAM_COL_34) {
	bufptr0 = next_token(bufptr);
	bufptr = next_token(bufptr0);
	if (!bufptr) {
	  goto load_fam_ret_MISSING_TOKENS;
	}
	tmp_len = strlen_se(bufptr0) + 1;
	if (tmp_len > max_paternal_id_len) {
	  max_paternal_id_len = tmp_len;
	}
	tmp_len = strlen_se(bufptr) + 1;
	if (tmp_len > max_maternal_id_len) {
	  max_maternal_id_len = tmp_len;
	}
      }
      if (fam_cols & FAM_COL_5) {
	bufptr = next_token(bufptr);
      }
      if (tmp_fam_col_6) {
	bufptr = next_token(bufptr);
	if (no_more_tokens_kns(bufptr)) {
	  goto load_fam_ret_MISSING_TOKENS;
	}
	// --1 forces case/control phenotype in plink 1.07, keep that for
	// backward compatibility
	if (affection && (!affection_01)) {
	  affection = eval_affection(bufptr, missing_phenod);
	}
      }
      unfiltered_sample_ct++;
    }
  }
  if (ferror(famfile)) {
    goto load_fam_ret_READ_FAIL;
  }
  if ((!unfiltered_sample_ct) && (!allow_no_samples)) {
    logerrprint("Error: Nobody in .fam file.\n");
    goto load_fam_ret_INVALID_FORMAT;
  }
  // don't yet need to enforce separate FID and IID limits, but in theory this
  // may change
  if ((max_sample_id_len > 2 * MAX_ID_BLEN) || (max_paternal_id_len > MAX_ID_BLEN) || (max_maternal_id_len > MAX_ID_BLEN)) {
    logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
    goto load_fam_ret_INVALID_FORMAT;
  }
  bigstack_reset(bigstack_mark);
  unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  // could make paternal_ids/maternal_ids conditional, but memory footprint is
  // typically negligible
  if (bigstack_alloc_c(unfiltered_sample_ct * max_sample_id_len, sample_ids_ptr) ||
      bigstack_alloc_c(unfiltered_sample_ct * max_paternal_id_len, paternal_ids_ptr) ||
      bigstack_alloc_c(unfiltered_sample_ct * max_maternal_id_len, maternal_ids_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, sex_nm_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, sex_male_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, founder_info_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, sample_exclude_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, pheno_nm_ptr)) {
    goto load_fam_ret_NOMEM;
  }

  // force either pheno_c or pheno_d to be allocated even if there is no
  // phenotype data
  if ((!tmp_fam_col_6) || affection) {
    if (aligned_malloc(unfiltered_sample_ctl * sizeof(intptr_t), pheno_c_ptr)) {
      goto load_fam_ret_NOMEM;
    }
    pheno_c = *pheno_c_ptr;
    fill_ulong_zero(unfiltered_sample_ctl, pheno_c);
  } else {
    pheno_d = (double*)malloc(unfiltered_sample_ct * sizeof(double));
    if (!pheno_d) {
      goto load_fam_ret_NOMEM;
    }
    fill_double_zero(unfiltered_sample_ct, pheno_d);
    *pheno_d_ptr = pheno_d;
  }
  bigstack_mark = g_bigstack_base;
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_fam_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  sample_ids = *sample_ids_ptr;
  paternal_ids = *paternal_ids_ptr;
  maternal_ids = *maternal_ids_ptr;
  if (!(fam_cols & FAM_COL_34)) {
    for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
      memcpy(&(paternal_ids[sample_uidx * max_paternal_id_len]), "0", 2);
    }
    for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
      memcpy(&(maternal_ids[sample_uidx * max_maternal_id_len]), "0", 2);
    }
  }
  sex_nm = *sex_nm_ptr;
  sex_male = *sex_male_ptr;
  pheno_nm = *pheno_nm_ptr;
  founder_info = *founder_info_ptr;
  if (fam_cols & FAM_COL_34) {
    fill_ulong_zero(unfiltered_sample_ctl, founder_info);
  } else {
    fill_all_bits(unfiltered_sample_ct, founder_info);
  }
  fill_ulong_zero(unfiltered_sample_ctl, sex_nm);
  fill_ulong_zero(unfiltered_sample_ctl, sex_male);
  fill_ulong_zero(unfiltered_sample_ctl, *sample_exclude_ptr);
  fill_ulong_zero(unfiltered_sample_ctl, pheno_nm);

  // ----- .fam read, second pass -----
  rewind(famfile);
  sample_uidx = 0;
  while (fgets(loadbuf, loadbuf_size, famfile)) {
    // no need to track line_idx here for now
    if (!loadbuf[loadbuf_size - 1]) {
      // should have been caught on first pass if loadbuf_size is still maximal
      goto load_fam_ret_NOMEM;
    }
    bufptr0 = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment_kns(*bufptr0)) {
      continue;
    }
    if (fam_cols & FAM_COL_1) {
      bufptr = next_token(bufptr0);
    } else {
      bufptr = bufptr0;
    }
    tmp_len = strlen_se(bufptr0);
    memcpyx(memcpyax(&(sample_ids[sample_uidx * max_sample_id_len]), bufptr0, tmp_len, '\t'), bufptr, strlen_se(bufptr), '\0');
    if (fam_cols & FAM_COL_34) {
      bufptr = next_token(bufptr);
      cc = *bufptr;
      tmp_len = strlen_se(bufptr);
      memcpyx(&(paternal_ids[sample_uidx * max_paternal_id_len]), bufptr, tmp_len, '\0');
      bufptr = next_token(bufptr);
      tmp_len2 = strlen_se(bufptr);
      memcpyx(&(maternal_ids[sample_uidx * max_maternal_id_len]), bufptr, tmp_len2, '\0');
      if ((tmp_len == 1) && (tmp_len2 == 1) && (cc == '0') && (*bufptr == '0')) {
	SET_BIT(sample_uidx, founder_info);
      }
    }
    if (fam_cols & FAM_COL_5) {
      bufptr = next_token(bufptr);
      if (strlen_se(bufptr) == 1) {
	if (*bufptr == '1') {
	  SET_BIT(sample_uidx, sex_nm);
	  SET_BIT(sample_uidx, sex_male);
	} else if (*bufptr == '2') {
	  SET_BIT(sample_uidx, sex_nm);
	}
      }
    }
    if (tmp_fam_col_6) {
      bufptr = next_token(bufptr);
      if (!scan_double(bufptr, &dxx)) {
	if (affection) {
	  if (dxx == pheno_ctrld) {
	    SET_BIT(sample_uidx, pheno_nm);
	  } else if (dxx == pheno_cased) {
	    SET_BIT(sample_uidx, pheno_nm);
	    SET_BIT(sample_uidx, pheno_c);
	  }
	} else if (dxx != missing_phenod) {
	  pheno_d[sample_uidx] = dxx;
	  SET_BIT(sample_uidx, pheno_nm);
	}
      }
    }
    sample_uidx++;
  }

  *unfiltered_sample_ct_ptr = unfiltered_sample_ct;
  *max_sample_id_len_ptr = max_sample_id_len;
  *max_paternal_id_len_ptr = max_paternal_id_len;
  *max_maternal_id_len_ptr = max_maternal_id_len;
  *affection_ptr = affection;
  if (fclose_null(&famfile)) {
    goto load_fam_ret_READ_FAIL;
  }
  while (0) {
  load_fam_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_fam_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_fam_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_fam_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .fam file has fewer tokens than expected.\n", line_idx);
  load_fam_ret_INVALID_FORMAT_2:
    logerrprintb();
  load_fam_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(famfile);
  return retval;
}

#define D_EPSILON 0.001220703125

int32_t oxford_to_bed(char* genname, char* samplename, char* outname, char* outname_end, char* single_chr, char* pheno_name, double hard_call_threshold, char* missing_code, int32_t missing_pheno, uint64_t misc_flags, uint32_t is_bgen, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  FILE* outfile_bim = nullptr;
  uintptr_t mc_ct = 0;
  uintptr_t max_mc_len = 0;
  uintptr_t line_idx = 0;
  double hard_call_floor = 1.0 - hard_call_threshold;
  char* loadbuf = nullptr;
  char* sorted_mc = nullptr;
  char* tbuf2 = &(g_textbuf[MAXLINELEN]); // .fam write

  // 0 = not present, otherwise zero-based index (this is fine since first
  //     column has to be FID)
  uint32_t sex_col = 0;

  uint32_t pheno_name_len = 0;
  // load first phenotype (not covariate), if present
  uint32_t pheno_col = 0;

  uint32_t snpid_chr = (misc_flags / MISC_OXFORD_SNPID_CHR) & 1;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t sample_ct = 0;
  uint32_t col_ct = 3;
  uint32_t is_binary_pheno = 0;
  uint32_t is_randomized = (hard_call_threshold == -1);
  uint32_t bgen_hardthresh = 0;
  uint32_t marker_ct = 0;
  int32_t retval = 0;
  uint32_t uint_arr[5];
  char missing_pheno_str[12];
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* wptr;
  uintptr_t* writebuf;
  uintptr_t* ulptr;
  uint16_t* bgen_probs;
  uint16_t* usptr;
  uintptr_t loadbuf_size;
  uintptr_t slen;
  uintptr_t cur_word;
  uintptr_t ulii;
  uintptr_t uljj;
  double dxx;
  double dyy;
  double dzz;
  double drand;
  uLongf zlib_ulongf;
  uint32_t missing_pheno_len;
  uint32_t raw_marker_ct;
  uint32_t marker_uidx;
  uint32_t sample_ct4;
  uint32_t sample_ctl2;
  uint32_t sample_idx;
  uint32_t col_idx;
  uint32_t shiftval;
  uint32_t bgen_compressed;
  uint32_t bgen_multichar_alleles;
  uint32_t identical_alleles;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t ii;
  uint16_t usii;
  uint16_t usjj;
  uint16_t uskk;
  char cc;
  char cc2;
  {
    if (single_chr && (!allow_extra_chroms)) {
      ii = get_chrom_code_raw(single_chr);
      if (!is_set(chrom_info_ptr->chrom_mask, ii)) {
	logerrprint("Error: --oxford-single-chr chromosome code is excluded by chromosome filter.\n");
	goto oxford_to_bed_ret_INVALID_CMDLINE;
      }
    }
    bufptr = int32toa(missing_pheno, missing_pheno_str);
    missing_pheno_len = (uintptr_t)(bufptr - missing_pheno_str);
    if (!missing_code) {
      mc_ct = 1;
      max_mc_len = 3;
      if (bigstack_alloc_c(3, &sorted_mc)) {
	goto oxford_to_bed_ret_NOMEM;
      }
      memcpy(sorted_mc, "NA", 3);
    } else {
      bufptr = missing_code;
      while (*bufptr) {
	while (*bufptr == ',') {
	  bufptr++;
	}
	if (!(*bufptr)) {
	  break;
	}
	mc_ct++;
	bufptr2 = strchr(bufptr, ',');
	if (!bufptr2) {
	  bufptr2 = strchr(bufptr, '\0');
	}
	ulii = (uintptr_t)(bufptr2 - bufptr);
	if (ulii >= max_mc_len) {
	  max_mc_len = ulii + 1;
	}
	bufptr = bufptr2;
      }
      if (mc_ct) {
	if (bigstack_alloc_c(mc_ct * max_mc_len, &sorted_mc)) {
	  goto oxford_to_bed_ret_NOMEM;
	}
	bufptr = missing_code;
	ulii = 0; // current missing-code index
	do {
	  while (*bufptr == ',') {
	    bufptr++;
	  }
	  if (!(*bufptr)) {
	    break;
	  }
	  bufptr2 = strchr(bufptr, ',');
	  if (!bufptr2) {
	    bufptr2 = strchr(bufptr, '\0');
	  }
	  uljj = (uintptr_t)(bufptr2 - bufptr);
	  memcpyx(&(sorted_mc[ulii * max_mc_len]), bufptr, uljj, '\0');
	  bufptr = bufptr2;
	  ulii++;
	} while (*bufptr);
	qsort(sorted_mc, mc_ct, max_mc_len, strcmp_casted);
      }
    }
    if (fopen_checked(samplename, "r", &infile)) {
      goto oxford_to_bed_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".fam", 5);
    if (fopen_checked(outname, "w", &outfile)) {
      goto oxford_to_bed_ret_OPEN_FAIL;
    }
    g_textbuf[MAXLINELEN - 1] = ' ';
    do {
      line_idx++;
      if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	if (ferror(infile)) {
	  goto oxford_to_bed_ret_READ_FAIL;
	}
	logerrprint("Error: Empty --data/--sample file.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (!g_textbuf[MAXLINELEN - 1]) {
	goto oxford_to_bed_ret_SAMPLE_LONG_LINE;
      }
      bufptr = skip_initial_spaces(g_textbuf);
    } while (is_eoln_kns(*bufptr));
    bufptr2 = token_endnn(bufptr);
    if ((((uintptr_t)(bufptr2 - bufptr)) != 4) || memcmp(bufptr, "ID_1", 4)) {
      goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_1;
    }
    bufptr = skip_initial_spaces(bufptr2);
    slen = strlen_se(bufptr);
    if ((slen != 4) || memcmp(bufptr, "ID_2", 4)) {
      goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_1;
    }
    bufptr = skip_initial_spaces(&(bufptr[4]));
    slen = strlen_se(bufptr);
    if ((slen != 7) || (!match_upper_counted(bufptr, "MISSING", 7))) {
      goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_1;
    }
    bufptr = skip_initial_spaces(&(bufptr[7]));
    if (pheno_name) {
      pheno_name_len = strlen(pheno_name);
    }
    while (!is_eoln_kns(*bufptr)) {
      bufptr2 = token_endnn(bufptr);
      ulii = (uintptr_t)(bufptr2 - bufptr);
      // allow "Sex", "SEX", etc.
      if ((ulii == 3) && (tolower(bufptr[0]) == 's') && (tolower(bufptr[1]) == 'e') && (tolower(bufptr[2]) == 'x')) {
	if (sex_col) {
	  goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_1;
	}
	sex_col = col_ct;
      } else if ((ulii == pheno_name_len) && (!memcmp(bufptr, pheno_name, ulii))) {
	if (pheno_col) {
	  goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_1;
	}
	pheno_col = col_ct;
      }
      col_ct++;
      bufptr = skip_initial_spaces(bufptr2);
    }
    if (pheno_name) {
      if (!pheno_col) {
	logerrprint("Error: --oxford-pheno-name parameter not found in .sample file header.\n");
	goto oxford_to_bed_ret_INVALID_CMDLINE;
      } else if (sex_col > pheno_col) {
	logerrprint("Error: .sample phenotype column(s) should be after sex covariate.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
    }
    do {
      line_idx++;
      if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	if (ferror(infile)) {
	  goto oxford_to_bed_ret_READ_FAIL;
	}
	logerrprint("Error: Only one nonempty line in .sample file.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (!g_textbuf[MAXLINELEN - 1]) {
	goto oxford_to_bed_ret_SAMPLE_LONG_LINE;
      }
      bufptr = skip_initial_spaces(g_textbuf);
    } while (is_eoln_kns(*bufptr));
    bufptr2 = token_endnn(bufptr);
    if ((((uintptr_t)(bufptr2 - bufptr)) != 1) || (*bufptr != '0')) {
      goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_2;
    }
    bufptr = skip_initial_spaces(bufptr2);
    slen = strlen_se(bufptr);
    if ((slen != 1) || (*bufptr != '0')) {
      goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_2;
    }
    bufptr = skip_initial_spaces(&(bufptr[1]));
    slen = strlen_se(bufptr);
    if ((slen != 1) || (*bufptr != '0')) {
      goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_2;
    }
    bufptr++;
    col_idx = 3;
    while (col_idx < col_ct) {
      bufptr = skip_initial_spaces(bufptr);
      if (is_eoln_kns(*bufptr)) {
	logerrprint("Error: Second .sample header line has fewer tokens than the first.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (bufptr[1] > ' ') {
	goto oxford_to_bed_ret_INVALID_SAMPLE_HEADER_2;
      }
      cc = *bufptr;
      if ((col_idx == sex_col) && (cc != 'D')) {
	logerrprint("Error: .sample sex column is not of type 'D'.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (!pheno_col) {
	if ((cc == 'B') || (cc == 'P')) {
	  if (sex_col > col_idx) {
	    logerrprint("Error: .sample phenotype column(s) should be after sex covariate.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  pheno_col = col_idx;
	  is_binary_pheno = (cc == 'B');
	  break;
	}
      } else if (col_idx == pheno_col) {
	is_binary_pheno = (cc == 'B');
	if ((!is_binary_pheno) && (cc != 'P')) {
	  logerrprint("Error: --oxford-pheno-name parameter does not refer to a binary or continuous\nphenotype.\n");
	  goto oxford_to_bed_ret_INVALID_CMDLINE;
	}
	break;
      }
      col_idx++;
      bufptr++;
    }
    if (is_binary_pheno) {
      // check for pathological case
      if ((bsearch_str("0", 1, sorted_mc, max_mc_len, mc_ct) != -1) || (bsearch_str("1", 1, sorted_mc, max_mc_len, mc_ct) != -1)) {
	logerrprint("Error: '0' and '1' are unacceptable missing case/control phenotype codes.\n");
	goto oxford_to_bed_ret_INVALID_CMDLINE;
      }
    }
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
	goto oxford_to_bed_ret_SAMPLE_LONG_LINE;
      }
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr2 = token_endnn(bufptr);
      wptr = memcpyax(tbuf2, bufptr, bufptr2 - bufptr, '\t');
      bufptr = skip_initial_spaces(bufptr2);
      if (is_eoln_kns(*bufptr)) {
	goto oxford_to_bed_ret_MISSING_TOKENS;
      }
      bufptr2 = token_endnn(bufptr);
      wptr = memcpya(wptr, bufptr, bufptr2 - bufptr);
      wptr = memcpya(wptr, "\t0\t0\t", 5);
      col_idx = 2;
      bufptr = bufptr2;
      if (sex_col) {
	while (1) {
	  bufptr = skip_initial_spaces(bufptr);
	  if (is_eoln_kns(*bufptr)) {
	    goto oxford_to_bed_ret_MISSING_TOKENS;
	  }
	  if (col_idx == sex_col) {
	    break;
	  }
	  bufptr = token_endnn(bufptr);
	  col_idx++;
	}
	cc = *bufptr++;
	if ((cc < '0') || (cc > '2') || ((*bufptr) > ' ')) {
	  sprintf(g_logbuf, "Error: Invalid sex code on line %" PRIuPTR " of .sample file.\n", line_idx);
	  goto oxford_to_bed_ret_INVALID_FORMAT_2;
	}
	*wptr++ = cc;
	col_idx++;
      } else {
	*wptr++ = '0';
      }
      *wptr++ = '\t';
      if (pheno_col) {
	while (1) {
	  bufptr = skip_initial_spaces(bufptr);
	  if (is_eoln_kns(*bufptr)) {
	    goto oxford_to_bed_ret_MISSING_TOKENS;
	  }
	  if (col_idx == pheno_col) {
	    break;
	  }
	  bufptr = token_endnn(bufptr);
	  col_idx++;
	}
	slen = (uintptr_t)(token_endnn(bufptr) - bufptr);
	if (is_binary_pheno) {
	  cc = *bufptr;
	  if ((slen != 1) || ((cc != '0') && (cc != '1'))) {
	    goto oxford_to_bed_missing_pheno;
	  } else {
	    *wptr++ = cc + 1;
	  }
	} else {
	  if (bsearch_str(bufptr, slen, sorted_mc, max_mc_len, mc_ct) == -1) {
	    if (!scan_double(bufptr, &dxx)) {
	      wptr = memcpya(wptr, bufptr, slen);
	    } else {
	      goto oxford_to_bed_missing_pheno;
	    }
	  } else {
	    goto oxford_to_bed_missing_pheno;
	  }
	}
      } else {
      oxford_to_bed_missing_pheno:
	wptr = memcpya(wptr, missing_pheno_str, missing_pheno_len);
      }
      *wptr++ = '\n';
      if (fwrite_checked(tbuf2, wptr - tbuf2, outfile)) {
	goto oxford_to_bed_ret_WRITE_FAIL;
      }
      sample_ct++;
    }
    if ((!sample_ct) && (!allow_no_samples)) {
      logerrprint("Error: No samples in .sample file.\n");
      goto oxford_to_bed_ret_INVALID_FORMAT;
    }
    if (fclose_null(&infile)) {
      goto oxford_to_bed_ret_READ_FAIL;
    }
    if (fclose_null(&outfile)) {
      goto oxford_to_bed_ret_WRITE_FAIL;
    }
    sample_ct4 = (sample_ct + 3) / 4;
    sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    if (bigstack_alloc_ul(sample_ctl2, &writebuf)) {
      goto oxford_to_bed_ret_NOMEM;
    }
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(outname, "w", &outfile_bim)) {
      goto oxford_to_bed_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto oxford_to_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto oxford_to_bed_ret_WRITE_FAIL;
    }
    if (!is_bgen) {
      loadbuf_size = bigstack_left();
      if (loadbuf_size > MAXLINEBUFLEN) {
	loadbuf_size = MAXLINEBUFLEN;
      } else if (loadbuf_size <= MAXLINELEN) {
	goto oxford_to_bed_ret_NOMEM;
      }
      loadbuf = (char*)g_bigstack_base;
      retval = gzopen_read_checked(genname, &gz_infile);
      if (retval) {
	goto oxford_to_bed_ret_1;
      }
      loadbuf[loadbuf_size - 1] = ' ';
      line_idx = 0;
      while (1) {
	line_idx++;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_infile)) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  break;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  if (loadbuf_size == MAXLINEBUFLEN) {
	    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file is pathologically long.\n", line_idx);
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  goto oxford_to_bed_ret_NOMEM;
	}
	char* loadbuf_first_token = skip_initial_spaces(loadbuf);
	if (is_eoln_kns(*loadbuf_first_token)) {
	  continue;
	}
        char* first_token_end = token_endnn(loadbuf_first_token);
	if (!single_chr) {
          const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - loadbuf_first_token);
          *first_token_end = '\0';
	  int32_t cur_chrom_code;
	  retval = get_or_add_chrom_code(loadbuf_first_token, ".gen file", line_idx, chrom_name_slen, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
	  if (retval) {
	    if ((chrom_name_slen == 3) && (!memcmp(loadbuf_first_token, "---", 3))) {
	      logerrprint("(Did you forget --oxford-single-chr?)\n");
	    }
	    goto oxford_to_bed_ret_1;
	  }
	  if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	    continue;
	  }
          *first_token_end = ' ';
	}
	fill_ulong_zero(sample_ctl2, writebuf);
	if (single_chr) {
	  fputs(single_chr, outfile_bim);
	  putc_unlocked(' ', outfile_bim);
	  bufptr = next_token(first_token_end);
	  bufptr2 = next_token(bufptr);
	} else {
	  bufptr = loadbuf_first_token;
	  bufptr2 = next_token(skip_initial_spaces(first_token_end));
	}
	if (!bufptr2) {
	  goto oxford_to_bed_ret_MISSING_TOKENS_GEN;
	}
	fwrite(bufptr, 1, bufptr2 - bufptr, outfile_bim);
	putc_unlocked('0', outfile_bim);
	if (putc_checked(' ', outfile_bim)) {
	  goto oxford_to_bed_ret_WRITE_FAIL;
	}
	bufptr = next_token(bufptr2);
	bufptr3 = next_token(bufptr);
	if (no_more_tokens_kns(bufptr3)) {
	  goto oxford_to_bed_ret_MISSING_TOKENS_GEN;
	}
	// bufptr2 = pos
	// bufptr  = allele 1
	// bufptr3 = allele 2
	bufptr4 = token_endnn(bufptr3);
	uii = (uintptr_t)(bufptr4 - bufptr3);
	identical_alleles = (strlen_se(bufptr) == uii) && (!memcmp(bufptr, bufptr3, uii));
	if (identical_alleles) {
	  // we treat identical A1 and A2 as a special case, since naive
	  // handling prevents e.g. later data merge.
	  // maybe add a warning?
	  fwrite(bufptr2, 1, strlen_se(bufptr2), outfile_bim);
	  fputs(" 0 ", outfile_bim);
	  fwrite(bufptr3, 1, bufptr4 - bufptr3, outfile_bim);
	} else {
	  fwrite(bufptr2, 1, bufptr4 - bufptr2, outfile_bim);
	}
	if (putc_checked('\n', outfile_bim)) {
	  goto oxford_to_bed_ret_WRITE_FAIL;
	}
	if (sample_ct) {
	  cur_word = 0;
	  shiftval = 0;
	  ulptr = writebuf;
	  bufptr = skip_initial_spaces(bufptr4);
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    if (is_eoln_kns(*bufptr)) {
	      goto oxford_to_bed_ret_MISSING_TOKENS_GEN;
	    }
	    // fast handling of common cases
	    cc = bufptr[1];
	    if ((cc == ' ') || (cc == '\t')) {
	      cc = bufptr[3];
	      cc2 = bufptr[5];
	      if (((cc == ' ') || (cc == '\t')) && ((cc2 == ' ') || (cc2 == '\t'))) {
		cc = *bufptr;
		if (cc == '0') {
		  bufptr2 = &(bufptr[2]);
		  cc = *bufptr2;
		  cc2 = bufptr2[2];
		  if (cc == '0') {
		    if (cc2 == '1') {
		      ulii = 3;
		    } else if (cc2 == '0') {
		      ulii = 1;
		    } else {
		      // could be a space...
		      goto oxford_to_bed_full_parse_2;
		    }
		  } else if ((cc == '1') && (cc2 == '0')) {
		    ulii = 2;
		  } else {
		    goto oxford_to_bed_full_parse_2;
		  }
		} else if ((cc == '1') && (bufptr[2] == '0') && (bufptr[4] == '0')) {
		  ulii = 0;
		} else {
		  goto oxford_to_bed_full_parse;
		}
		bufptr = &(bufptr[6]);
	      } else {
		goto oxford_to_bed_full_parse;
	      }
	    } else {
	      // okay, gotta do things the slow way
	    oxford_to_bed_full_parse:
	      bufptr2 = token_endnn(bufptr);
	    oxford_to_bed_full_parse_2:
	      bufptr2 = skip_initial_spaces(bufptr2);
	      if (is_eoln_kns(*bufptr2)) {
		goto oxford_to_bed_ret_MISSING_TOKENS_GEN;
	      }
	      bufptr3 = token_endnn(bufptr2);
	      dzz = strtod(bufptr3, &bufptr4);
	      if (!is_randomized) {
		if (dzz >= hard_call_floor) {
		  ulii = 3;
		} else {
		  if (bufptr3 == bufptr4) {
		    goto oxford_to_bed_ret_INVALID_DOSAGE;
		  }
		  dyy = strtod(bufptr2, &bufptr3);
		  if (dyy >= hard_call_floor) {
		    ulii = 2;
		  } else {
		    if (bufptr2 == bufptr3) {
		      goto oxford_to_bed_ret_INVALID_DOSAGE;
		    }
		    dxx = strtod(bufptr, &bufptr2);
		    if (dxx >= hard_call_floor) {
		      ulii = 0;
		    } else {
		      if (bufptr == bufptr2) {
			goto oxford_to_bed_ret_INVALID_DOSAGE;
		      }
		      ulii = 1;
		    }
		  }
		}
	      } else {
		drand = rand_unif();
		if (drand < dzz) {
		  ulii = 3;
		} else {
		  if (bufptr3 == bufptr4) {
		    goto oxford_to_bed_ret_INVALID_DOSAGE;
		  }
		  dyy = strtod(bufptr2, &bufptr3) + dzz;
		  if (drand < dyy) {
		    ulii = 2;
		  } else {
		    if (bufptr2 == bufptr3) {
		      goto oxford_to_bed_ret_INVALID_DOSAGE;
		    }
		    dxx = strtod(bufptr, &bufptr2) + dyy;
		    if (drand < dxx) {
		      ulii = 0;
		    } else if (dxx < 1 - D_EPSILON) {
		      ulii = 1;
		    } else {
		      // fully called genotype probabilities may add up to
		      // less than one due to rounding error.  If this
		      // appears to have happened, do NOT make a missing
		      // call; instead rescale everything to add to one and
		      // reinterpret the random number.  (D_EPSILON is
		      // currently set to make 3 decimal place precision safe
		      // to use.)
		      drand *= dxx;
		      if (drand < dzz) {
			ulii = 3;
		      } else if (drand < dyy) {
			ulii = 2;
		      } else {
			ulii = 0;
		      }
		    }
		  }
		}
	      }
	      bufptr = skip_initial_spaces(bufptr4);
	    }
	    cur_word |= ulii << shiftval;
	    shiftval += 2;
	    if (shiftval == BITCT) {
	      *ulptr++ = cur_word;
	      cur_word = 0;
	      shiftval = 0;
	    }
	  }
	  if (shiftval) {
	    *ulptr++ = cur_word;
	  }
	  if (identical_alleles) {
	    // keep missing calls, but convert hom/het A1 to hom A2.
	    for (ulptr = writebuf; ulptr < (&(writebuf[sample_ctl2])); ulptr++) {
	      ulii = *ulptr;
	      *ulptr = ((~ulii) << 1) | ulii | FIVEMASK;
	    }
	    if (sample_ct % 4) {
	      writebuf[sample_ctl2 - 1] &= (ONELU << (2 * (sample_ct % BITCT2))) - ONELU;
	    }
	  }
	  if (fwrite_checked(writebuf, sample_ct4, outfile)) {
	    goto oxford_to_bed_ret_WRITE_FAIL;
	  }
	}
	marker_ct++;
	if (!(marker_ct % 1000)) {
	  printf("\r--data: %uk variants converted.", marker_ct / 1000);
	  fflush(stdout);
	}
      }
      if ((!marker_ct) && (!allow_no_variants)) {
	logerrprint("Error: Empty .gen file.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
    } else {
      if (fopen_checked(genname, FOPEN_RB, &infile)) {
	goto oxford_to_bed_ret_OPEN_FAIL;
      }
      // supports BGEN v1.0 and v1.1.
      bgen_probs = (uint16_t*)bigstack_alloc(6LU * sample_ct);
      if (!bgen_probs) {
	goto oxford_to_bed_ret_NOMEM;
      }
      loadbuf = (char*)g_bigstack_base;
      loadbuf_size = bigstack_left();
      if (loadbuf_size > MAXLINEBUFLEN) {
	loadbuf_size = MAXLINEBUFLEN;
      } else if (loadbuf_size < 3 * 65536) {
	goto oxford_to_bed_ret_NOMEM;
      }
      if (fread(uint_arr, 1, 20, infile) < 20) {
	goto oxford_to_bed_ret_READ_FAIL;
      }
      if (uint_arr[1] > uint_arr[0]) {
	logerrprint("Error: Invalid .bgen header.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      raw_marker_ct = uint_arr[2];
      if ((!raw_marker_ct) && (!allow_no_variants)) {
	logerrprint("Error: .bgen file contains no variants.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (uint_arr[3] != sample_ct) {
	logerrprint("Error: --bgen and --sample files contain different numbers of samples.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (uint_arr[4] && (uint_arr[4] != 0x6e656762)) {
	logerrprint("Error: Invalid .bgen magic number.\n");
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (fseeko(infile, uint_arr[1], SEEK_SET)) {
	goto oxford_to_bed_ret_READ_FAIL;
      }
      if (fread(&uii, 1, 4, infile) < 4) {
	goto oxford_to_bed_ret_READ_FAIL;
      }
      if (uii & (~5)) {
	const uint32_t layout = (uii >> 2) & 15;
	const uint32_t compression_mode = uii & 3;
	if ((layout == 2) && (compression_mode != 3)) {
	  LOGERRPRINTF("Error: BGEN v1.%c input requires PLINK 2.0 (under development as of this\nwriting).  Use qctool2 or bgenix to downcode to BGEN v1.1 if you want to\nprocess this data with PLINK 1.9.\n", (compression_mode == 2)? '3' : '2');
	} else if (layout > 2) {
	  logerrprint("Error: Unrecognized BGEN version.  Use bgenix or a similar tool to downcode to\nBGEN v1.1 if you want to process this data with PLINK 1.9.\n");
	} else {
	  logerrprint("Error: Unrecognized flags in .bgen header.  (PLINK 1.9 only supports BGEN v1.0\nand v1.1.)\n");
	}
	goto oxford_to_bed_ret_INVALID_FORMAT;
      }
      if (fseeko(infile, 4 + uint_arr[0], SEEK_SET)) {
	goto oxford_to_bed_ret_READ_FAIL;
      }
      bgen_compressed = uii & 1;
      bgen_multichar_alleles = (uii >> 2) & 1;
      if ((!bgen_multichar_alleles) && (!snpid_chr) && (chrom_info_ptr->species != SPECIES_HUMAN)) {
	logerrprint("Error: BGEN v1.0 files can only support nonhuman genomes if the SNP ID field is\nused for chromosome codes.\n");
	goto oxford_to_bed_ret_INVALID_CMDLINE;
      }
      if (!is_randomized) {
	bgen_hardthresh = 32768 - (int32_t)(hard_call_threshold * 32768);
      }
      memcpyl3(g_textbuf, " 0 ");
      for (marker_uidx = 0; marker_uidx < raw_marker_ct; marker_uidx++) {
	if (fread(&uii, 1, 4, infile) < 4) {
	  goto oxford_to_bed_ret_READ_FAIL;
	}
	if (uii != sample_ct) {
	  logerrprint("Error: Unexpected number of samples specified in SNP block header.\n");
	  goto oxford_to_bed_ret_INVALID_FORMAT;
	}
	if (bgen_multichar_alleles) {
	  // v1.1
	  if (fread(&usii, 1, 2, infile) < 2) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (!snpid_chr) {
	    if (fseeko(infile, usii, SEEK_CUR)) {
	      goto oxford_to_bed_ret_READ_FAIL;
	    }
	    bufptr = loadbuf;
	  } else {
	    if (!usii) {
	      logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
	      goto oxford_to_bed_ret_INVALID_FORMAT;
	    }
	    if (fread(loadbuf, 1, usii, infile) < usii) {
	      goto oxford_to_bed_ret_READ_FAIL;
	    }
	    loadbuf[usii] = '\0';
	    bufptr = &(loadbuf[usii + 1]);
	  }
	  if (fread(&usjj, 1, 2, infile) < 2) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (!usjj) {
	    logerrprint("Error: Length-0 rsID in .bgen file.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  if (fread(bufptr, 1, usjj, infile) < usjj) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  bufptr2 = &(bufptr[usjj]);
	  if (fread(&uskk, 1, 2, infile) < 2) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (!snpid_chr) {
	    if (!uskk) {
	      logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
	      goto oxford_to_bed_ret_INVALID_FORMAT;
	    }
	    usii = uskk;
	    if (fread(bufptr2, 1, usii, infile) < usii) {
	      goto oxford_to_bed_ret_READ_FAIL;
	    }
	    if ((usii == 2) && (!memcmp(bufptr2, "NA", 2))) {
	      // convert 'NA' to 0
	      usii = 1;
	      memcpy(bufptr2, "0", 2);
	    } else {
	      bufptr2[usii] = '\0';
	    }
	  } else {
	    if (fseeko(infile, uskk, SEEK_CUR)) {
	      goto oxford_to_bed_ret_READ_FAIL;
	    }
	    bufptr2 = loadbuf;
	  }
	  if (fread(uint_arr, 1, 8, infile) < 8) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (!uint_arr[1]) {
	    logerrprint("Error: Length-0 allele ID in .bgen file.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  // bufptr2 (chromosome code) is already zero-terminated, with known
	  // length usii
	  int32_t cur_chrom_code;
	  retval = get_or_add_chrom_code(bufptr2, ".bgen file", 0, usii, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
	  if (retval) {
	    goto oxford_to_bed_ret_1;
	  }
	  if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	    // skip rest of current SNP
	    if (fseeko(infile, uint_arr[1], SEEK_CUR)) {
	      goto oxford_to_bed_ret_READ_FAIL;
	    }
	    if (fread(&uii, 1, 4, infile) < 4) {
	      goto oxford_to_bed_ret_READ_FAIL;
	    }
	    if (bgen_compressed) {
	      if (fseeko(infile, uii, SEEK_CUR)) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	      if (fread(&uii, 1, 4, infile) < 4) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	      if (fseeko(infile, uii, SEEK_CUR)) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	    } else {
	      if (fseeko(infile, uii + ((uint64_t)sample_ct) * 6, SEEK_CUR)) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	    }
	    continue;
	  }
	  fputs(bufptr2, outfile_bim);
	  if (putc_checked(' ', outfile_bim)) {
	    goto oxford_to_bed_ret_WRITE_FAIL;
	  }
	  fwrite(bufptr, 1, usjj, outfile_bim);
	  bufptr = uint32toa_x(uint_arr[0], ' ', &(g_textbuf[3]));
	  fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile_bim);

	  // halve the limit since there are two alleles
	  // (may want to enforce NON_BIGSTACK_MIN allele length limit?)
	  if (uint_arr[1] >= loadbuf_size / 2) {
	    if (loadbuf_size < MAXLINEBUFLEN) {
	      goto oxford_to_bed_ret_NOMEM;
	    }
	    logerrprint("Error: Excessively long allele in .bgen file.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  if (fread(loadbuf, 1, uint_arr[1], infile) < uint_arr[1]) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  loadbuf[uint_arr[1]] = ' ';
	  if (fread(&uii, 1, 4, infile) < 4) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (uii >= loadbuf_size / 2) {
	    if (loadbuf_size < MAXLINEBUFLEN) {
	      goto oxford_to_bed_ret_NOMEM;
	    }
	    logerrprint("Error: Excessively long allele in .bgen file.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  bufptr = &(loadbuf[uint_arr[1] + 1]);
	  if (fread(bufptr, 1, uii, infile) < uii) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  bufptr[uii] = '\n';
	  identical_alleles = (uii == uint_arr[1]) && (!memcmp(loadbuf, bufptr, uii));
	  if (!identical_alleles) {
	    if (fwrite_checked(loadbuf, uint_arr[1] + uii + 2, outfile_bim)) {
	      goto oxford_to_bed_ret_WRITE_FAIL;
	    }
	  } else {
	    fputs("0 ", outfile_bim);
	    if (fwrite_checked(bufptr, uii + 1, outfile_bim)) {
	      goto oxford_to_bed_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  // v1.0
	  uii = 0;
	  if (fread(&uii, 1, 1, infile) < 1) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (fread(loadbuf, 1, 2 * uii + 9, infile) < (2 * uii + 9)) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  // save marker ID length since we might clobber it
	  ukk = (unsigned char)(loadbuf[uii + 1]);
	  int32_t cur_chrom_code;
	  if (!snpid_chr) {
	    cur_chrom_code = ((unsigned char)(loadbuf[2 * uii + 2]));
	    if (cur_chrom_code > 24) {
	      if (cur_chrom_code == 255) {
		// unknown
		cur_chrom_code = 0;
	      } else if (cur_chrom_code > 252) {
		// XY or MT
		cur_chrom_code = cur_chrom_code - 228;
	      } else {
		logerrprint("Error: Invalid chromosome code in BGEN v1.0 file.\n");
		goto oxford_to_bed_ret_INVALID_FORMAT;
	      }
	    }
	    uint32toa_x((uint32_t)cur_chrom_code, '\0', loadbuf);
	    bufptr = loadbuf;
	  } else {
	    ujj = (unsigned char)loadbuf[0];
	    bufptr = &(loadbuf[1]);
	    if ((ujj == 2) && (!memcmp(bufptr, "NA", 2))) {
	      *bufptr = '0';
	      ujj = 1;
	    }
	    bufptr[ujj] = '\0';
	    retval = get_or_add_chrom_code(bufptr, ".bgen file", 0, ujj, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
	    if (retval) {
	      goto oxford_to_bed_ret_1;
	    }
	  }
	  if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	    if (bgen_compressed) {
	      if (fread(&uii, 1, 4, infile) < 4) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	      if (fseeko(infile, uii, SEEK_CUR)) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	    } else {
	      if (fseeko(infile, ((uint64_t)sample_ct) * 6, SEEK_CUR)) {
		goto oxford_to_bed_ret_READ_FAIL;
	      }
	    }
	    continue;
	  }
	  fputs(bufptr, outfile_bim);
	  if (putc_checked(' ', outfile_bim)) {
	    goto oxford_to_bed_ret_WRITE_FAIL;
	  }
	  fwrite(&(loadbuf[uii + 2]), 1, ukk, outfile_bim);
	  memcpy(&ujj, &(loadbuf[2 * uii + 3]), 4);
	  bufptr = uint32toa_x(ujj, ' ', &(g_textbuf[3]));
	  identical_alleles = (loadbuf[2 * uii + 7] == loadbuf[2 * uii + 8]);
	  if (!identical_alleles) {
	    *bufptr++ = loadbuf[2 * uii + 7];
	  } else {
	    *bufptr++ = '0';
	  }
	  *bufptr++ = ' ';
	  *bufptr++ = loadbuf[2 * uii + 8];
	  *bufptr++ = '\n';
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_bim)) {
	    goto oxford_to_bed_ret_WRITE_FAIL;
	  }
	}
	if (bgen_compressed) {
	  if (fread(&uii, 1, 4, infile) < 4) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  if (uii > loadbuf_size) {
	    if (loadbuf_size < MAXLINEBUFLEN) {
	      goto oxford_to_bed_ret_NOMEM;
	    }
	    logerrprint("Error: Excessively long compressed SNP block in .bgen file.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	  if (fread(loadbuf, 1, uii, infile) < uii) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	  zlib_ulongf = 6 * sample_ct;
	  if (uncompress((Bytef*)bgen_probs, &zlib_ulongf, (Bytef*)loadbuf, uii) != Z_OK) {
	    logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
	    goto oxford_to_bed_ret_INVALID_FORMAT;
	  }
	} else {
	  if (fread(bgen_probs, 1, 6 * sample_ct, infile) < 6 * sample_ct) {
	    goto oxford_to_bed_ret_READ_FAIL;
	  }
	}
	cur_word = 0;
	shiftval = 0;
	ulptr = writebuf;
	usptr = bgen_probs;
	if (!is_randomized) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, usptr = &(usptr[3])) {
	    if (usptr[2] >= bgen_hardthresh) {
	      ulii = 3;
	    } else if (usptr[1] >= bgen_hardthresh) {
	      ulii = 2;
	    } else if (usptr[0] >= bgen_hardthresh) {
	      ulii = 0;
	    } else {
	      ulii = 1;
	    }
	    cur_word |= ulii << shiftval;
	    shiftval += 2;
	    if (shiftval == BITCT) {
	      *ulptr++ = cur_word;
	      cur_word = 0;
	      shiftval = 0;
	    }
	  }
	} else {
	  uii = 0;
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, usptr = &(usptr[3])) {
	    // fast handling of common cases
	    ukk = usptr[2];
	    if (ukk >= 32768) {
	      ulii = 3;
	    } else if (usptr[1] >= 32768) {
	      ulii = 2;
	    } else if (usptr[0] >= 32768) {
	      ulii = 0;
	    } else {
	      while (1) {
		uii >>= 16;
		if (!uii) {
		  uii = sfmt_genrand_uint32(&g_sfmt) | 0x80000000U;
		}
		ujj = uii & 32767;
		if (ujj < ukk) {
		  ulii = 3;
		  break;
		} else {
		  ukk += usptr[1];
		  if (ujj < ukk) {
		    ulii = 2;
		    break;
		  } else {
		    ukk += usptr[0];
		    if (ujj < ukk) {
		      ulii = 0;
		      break;
		    } else if (ukk < 32766) {
		      ulii = 1;
		      break;
		    } else {
		      ukk = usptr[2];
		    }
		  }
		}
	      }
	    }
	    cur_word |= ulii << shiftval;
	    shiftval += 2;
	    if (shiftval == BITCT) {
	      *ulptr++ = cur_word;
	      cur_word = 0;
	      shiftval = 0;
	    }
	  }
	}
	if (shiftval) {
	  *ulptr++ = cur_word;
	}
	if (identical_alleles) {
	  for (ulptr = writebuf; ulptr < (&(writebuf[sample_ctl2])); ulptr++) {
	    ulii = *ulptr;
	    *ulptr = ((~ulii) << 1) | ulii | FIVEMASK;
	  }
	  if (sample_ct % 4) {
	    writebuf[sample_ctl2 - 1] &= (ONELU << (2 * (sample_ct % BITCT2))) - ONELU;
	  }
	}
	if (fwrite_checked(writebuf, sample_ct4, outfile)) {
	  goto oxford_to_bed_ret_WRITE_FAIL;
	}
	marker_ct++;
	if (!(marker_ct % 1000)) {
	  if (marker_ct == marker_uidx + 1) {
	    printf("\r--bgen: %uk variants converted.", marker_ct / 1000);
	  } else {
	    printf("\r--bgen: %uk variants converted (out of %u).", marker_ct / 1000, marker_uidx + 1);
	  }
	  fflush(stdout);
	}
      }
      if (fclose_null(&infile)) {
	goto oxford_to_bed_ret_READ_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto oxford_to_bed_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile_bim)) {
      goto oxford_to_bed_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    LOGPRINTFWW("--%s: %s.bed + %s.bim + %s.fam written.\n", is_bgen? "bgen" : "data", outname, outname, outname);
  }
  while (0) {
  oxford_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  oxford_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  oxford_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  oxford_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  oxford_to_bed_ret_INVALID_DOSAGE:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file has an invalid dosage value.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_MISSING_TOKENS_GEN:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_SAMPLE_LONG_LINE:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file is pathologically long.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_INVALID_SAMPLE_HEADER_2:
    logerrprint("Error: Invalid second header line in .sample file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_INVALID_SAMPLE_HEADER_1:
    logerrprint("Error: Invalid first header line in .sample file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_INVALID_FORMAT_2:
    logerrprintb();
  oxford_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  oxford_to_bed_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 oxford_to_bed_ret_1:
  fclose_cond(infile);
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  fclose_cond(outfile_bim);
  bigstack_reset(bigstack_mark);
  return retval;
}

// side effect: initializes textbuf to first nonempty noncomment line of
// .map/.bim
int32_t check_cm_col(FILE* bimfile, char* textbuf, uint32_t is_binary, uint32_t allow_no_variants, uint32_t bufsize, uint32_t* cm_col_exists_ptr, uintptr_t* line_idx_ptr) {
  uintptr_t line_idx = 0;
  char* bufptr;
  while (fgets(textbuf, bufsize, bimfile)) {
    line_idx++;
    bufptr = skip_initial_spaces(textbuf);
    if (is_eoln_or_comment_kns(*bufptr)) {
      continue;
    }
    bufptr = next_token_mult(bufptr, 2 + 2 * is_binary);
    *line_idx_ptr = line_idx;
    if (no_more_tokens_kns(bufptr)) {
      return -1;
    }
    *cm_col_exists_ptr = !no_more_tokens_kns(next_token(bufptr));
    return 0;
  }
  *line_idx_ptr = 0;
  return allow_no_variants? 0 : -1;
}

int32_t incr_text_allele0(char cc, char* marker_alleles, uint32_t* marker_allele_cts) {
  uint32_t uii;
  for (uii = 0; uii < 4; uii++) {
    if (marker_alleles[uii] == '\0') {
      marker_alleles[uii] = cc;
      marker_allele_cts[uii] = 1;
      return 0;
    } else if (marker_alleles[uii] == cc) {
      marker_allele_cts[uii] += 1;
      return 0;
    }
  }
  return -1;
}

typedef struct ll_str_fixed_struct {
  struct ll_str_struct* next;
#ifdef __LP64__
  char ss[8];
#else
  char ss[12];
#endif
} Ll_str_fixed;

int32_t incr_text_allele_str(char* allele_name, uint32_t an_len, Ll_str* allele_list_start, uint32_t* marker_allele_cts) {
  // Start with preallocated array of 16-byte Ll_strs.
  // Ll_str.ss is a null-terminated sequence of ordered, tab-delimited allele
  // names.  If the starting 8 (or 12 bytes, on 32-bit systems) is adequate,
  // Ll_str.next is nullptr.  Otherwise, Ll_str.ss stores the first few (or
  // possibly 0, if the very first allele name is too long) allele names, and
  // Ll_str.next is a pointer to a linked list entry storing the next 1+ allele
  // names.  Worst case, the linked list is of length 4 (beyond that we error
  // out).  Most of the time, it'll be length 0.
  // Allele names are separated by tabs; \0 denotes the end of the allele name
  // list.
  // This type of function will become more important when it's time to parse
  // .vcf files, etc.
  uint32_t allele_num = 0;
  char* cur_allele_name_start = allele_list_start->ss;
  Ll_str* ll_ptr;
  uint32_t slen;
  uintptr_t chars_left;
  if (!(*cur_allele_name_start)) {
    if (!(allele_list_start->next)) {
      if (an_len >= (16 - sizeof(intptr_t))) {
	if (bigstack_end_alloc_llstr(an_len + 1, &ll_ptr)) {
	  return RET_NOMEM;
	}
        allele_list_start->next = ll_ptr;
	ll_ptr->next = nullptr;
	cur_allele_name_start = ll_ptr->ss;
      }
      memcpyx(cur_allele_name_start, allele_name, an_len, '\0');
      marker_allele_cts[0] = 1;
      return 0;
    }
    allele_list_start = allele_list_start->next;
    cur_allele_name_start = allele_list_start->ss;
  }
  while (allele_num < 4) {
    slen = strlen_se(cur_allele_name_start);
    if ((slen == an_len) && (!memcmp(cur_allele_name_start, allele_name, an_len))) {
      marker_allele_cts[allele_num] += 1;
      return 0;
    }
    allele_num++;
    if (cur_allele_name_start[slen] == '\t') {
      cur_allele_name_start = &(cur_allele_name_start[slen + 1]);
    } else if (allele_list_start->next) {
      allele_list_start = allele_list_start->next;
      cur_allele_name_start = allele_list_start->ss;
    } else {
      chars_left = ((0x7ffffff0 - sizeof(intptr_t)) - ((uintptr_t)(&(cur_allele_name_start[slen + 1]) - allele_list_start->ss))) & 15;
      if (chars_left > an_len) {
        cur_allele_name_start[slen] = '\t';
        memcpyx(&(cur_allele_name_start[slen + 1]), allele_name, an_len, '\0');
      } else {
        if (bigstack_end_alloc_llstr(an_len + 1, &ll_ptr)) {
          return RET_NOMEM;
        }
        allele_list_start->next = ll_ptr;
        ll_ptr->next = nullptr;
        cur_allele_name_start = ll_ptr->ss;
        memcpyx(cur_allele_name_start, allele_name, an_len, '\0');
      }
      marker_allele_cts[allele_num] = 1;
      return 0;
    }
  }
  return RET_INVALID_FORMAT;
}

char* get_llstr(Ll_str* ll_ptr, uint32_t allele_idx) {
  char* cptr = ll_ptr->ss;
  if (*cptr == '\0') {
    ll_ptr = ll_ptr->next;
    if (!ll_ptr) {
      return nullptr;
    }
    cptr = ll_ptr->ss;
  }
  while (allele_idx) {
    cptr = token_endnn(cptr);
    allele_idx--;
    if (*cptr) {
      cptr++;
    } else {
      ll_ptr = ll_ptr->next;
      if (!ll_ptr) {
	return nullptr;
      }
      cptr = ll_ptr->ss;
    }
  }
  return cptr;
}

static inline char* write_token_notab(char* read_ptr, FILE* outfile) {
  // assumes read_ptr is at the beginning of an item to write
  uint32_t slen = strlen_se(read_ptr);
  fwrite(read_ptr, 1, slen, outfile);
  return skip_initial_spaces(&(read_ptr[slen]));
}

static inline char* write_token(char* read_ptr, FILE* outfile) {
  uint32_t slen = strlen_se(read_ptr);
  fwrite(read_ptr, 1, slen, outfile);
  putc_unlocked('\t', outfile);
  return skip_initial_spaces(&(read_ptr[slen]));
}

int32_t ped_to_bed_multichar_allele(FILE** pedfile_ptr, FILE** outfile_ptr, char* outname, char* outname_end, FILE** mapfile_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_alleles_f, uint32_t map_is_unsorted, uint32_t fam_cols, uint32_t ped_col_skip_iid, uint32_t ped_col_skip, uint32_t cm_col_exists, uint32_t* map_reverse, int64_t ped_size, char* missing_pheno_str) {
  // maintain allele counts and linked lists of observed alleles at FAR end of
  // bigstack.
  unsigned char* bigstack_end_mark = g_bigstack_end;
  int32_t retval = 0;
  uint32_t ped_buflen = 0;
  uintptr_t sample_ct = 0;
  uintptr_t line_idx = 0;
  uint32_t pct = 1;
  int64_t ped_next_thresh = ped_size / 100;
  uint32_t last_pass = 0;
  int64_t* line_starts = nullptr;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;
  // do NOT convert missing -> output_missing when autoconverting, since the
  // .bim/.fam files are usually read right back in.
  char** marker_allele_ptrs = nullptr;
  FILE* outfile;
  uint32_t pass_ct;
  uintptr_t sample_ct4;
  char* loadbuf;
  char* col1_ptr;
  char* col2_ptr;
  char* bufptr;
  char* bufptr2;
  uintptr_t loadbuf_size;
  uintptr_t cur_slen;
  uintptr_t cur_slen_rdup;
  Ll_str_fixed* marker_alleles_tmp;
  uint32_t* marker_allele_cts;
  uint32_t markers_per_pass;
  uint32_t marker_start;
  uint32_t marker_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uint32_t loop_end;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t alen1;
  uint32_t alen2;
  char* aptr1;
  char* aptr2;
  unsigned char ucc;
  uint32_t ii_shift;
  unsigned char* writebuf;
  unsigned char* wbufptr;
  bigstack_reset(marker_alleles_f);
  loadbuf = (char*)g_bigstack_base;
  if (bigstack_end_calloc_ui(marker_ct * 4, &marker_allele_cts)) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  marker_alleles_tmp = (Ll_str_fixed*)bigstack_end_alloc(marker_ct * sizeof(Ll_str_fixed));
  if (!marker_alleles_tmp) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  memset(marker_alleles_tmp, 0, marker_ct * sizeof(Ll_str_fixed));

  if (fclose_null(outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(outname, "w", outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
  outfile = *outfile_ptr;
  rewind(*pedfile_ptr);
  fputs("Rescanning .ped file... 0%", stdout);
  fflush(stdout);
  while (1) {
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    }
    loadbuf[loadbuf_size - 1] = ' ';
  ped_to_bed_multichar_allele_loop_1_start:
    line_idx++;
    if (!fgets(loadbuf, loadbuf_size, *pedfile_ptr)) {
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      logprint("\n");
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file is pathologically long.\n", line_idx);
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
      } else {
        goto ped_to_bed_multichar_allele_ret_NOMEM;
      }
    }
    cur_slen = strlen(loadbuf);
    ulii = cur_slen + 1;
    if (ulii > ped_buflen) {
      ped_buflen = ulii;
    }
    col1_ptr = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment_kns(*col1_ptr)) {
      goto ped_to_bed_multichar_allele_loop_1_start;
    }
    // check for top-of-stack allocations colliding with load buffer
    cur_slen_rdup = (cur_slen + CACHELINE) & (~(CACHELINE - ONELU));
    if (fam_cols & FAM_COL_1) {
      col2_ptr = next_token(col1_ptr);
    } else {
      col2_ptr = col1_ptr;
    }
    bufptr = next_token_mult(col2_ptr, ped_col_skip_iid);
    if (no_more_tokens_kns(bufptr)) {
      goto ped_to_bed_multichar_allele_ret_MISSING_TOKENS;
    }
    if (fwrite_checked(col1_ptr, strlen_se(col1_ptr), outfile)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    putc_unlocked('\t', outfile);
    bufptr2 = write_token(col2_ptr, outfile);
    if (fam_cols & FAM_COL_34) {
      bufptr2 = write_token(bufptr2, outfile);
      bufptr2 = write_token(bufptr2, outfile);
    } else {
      fputs("0\t0\t", outfile);
    }
    if (fam_cols & FAM_COL_5) {
      bufptr2 = write_token_notab(bufptr2, outfile);
    } else {
      putc_unlocked('0', outfile);
    }
    putc_unlocked('\t', outfile);
    if (fam_cols & FAM_COL_6) {
      uii = strlen_se(bufptr2);
      fwrite(bufptr2, 1, uii, outfile);
    } else {
      fputs(missing_pheno_str, outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    g_bigstack_base += cur_slen_rdup;
    for (marker_uidx = 0, marker_idx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      alen1 = strlen_se(bufptr);
      aptr1 = bufptr;
      bufptr = skip_initial_spaces(&(bufptr[alen1]));
      alen2 = strlen_se(bufptr);
      if (!alen2) {
	g_bigstack_base -= cur_slen_rdup;
	goto ped_to_bed_multichar_allele_ret_MISSING_TOKENS;
      }
      aptr2 = bufptr;
      bufptr = skip_initial_spaces(&(bufptr[alen2]));
      if (IS_SET(marker_exclude, marker_uidx)) {
	continue;
      }
      if ((*aptr1 == missing_geno) && (alen1 == 1)) {
	if ((alen2 != 1) || (*aptr2 != missing_geno)) {
	  goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4;
	}
	marker_idx++;
	continue;
      } else if ((*aptr2 == missing_geno) && (alen2 == 1)) {
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4;
      }
      uii = map_is_unsorted? map_reverse[marker_idx] : marker_idx;
      retval = incr_text_allele_str(aptr1, alen1, (Ll_str*)(&(marker_alleles_tmp[uii])), &(marker_allele_cts[4 * uii]));
      if (retval) {
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_6;
      }
      retval = incr_text_allele_str(aptr2, alen2, (Ll_str*)(&(marker_alleles_tmp[uii])), &(marker_allele_cts[4 * uii]));
      if (retval) {
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_6;
      }
      marker_idx++;
    }
    g_bigstack_base -= cur_slen_rdup;
    if (!is_eoln_kns(*bufptr)) {
      logprint("\n");
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file has more tokens than expected.\n", line_idx);
      goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
    }
    sample_ct++;
    if (ftello(*pedfile_ptr) >= ped_next_thresh) {
      uii = (ftello(*pedfile_ptr) * 100) / ped_size;
      if (pct >= 10) {
	putc_unlocked('\b', stdout);
      }
      printf("\b\b%u%%", uii);
      fflush(stdout);
      pct = uii;
    }
  }
  if (!feof(*pedfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_READ_FAIL;
  }
  putc_unlocked('\r', stdout);
  logprint(".ped scan complete (for binary autoconversion).\n");
  // sample_ct == 0 impossible
  if (fclose_null(outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  marker_allele_ptrs = (char**)bigstack_alloc(marker_ct * 2 * sizeof(intptr_t));
  if (!marker_allele_ptrs) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(outname, "w", outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
  outfile = *outfile_ptr;
  if (map_is_unsorted) {
    memcpy(outname_end, ".map.tmp", 9);
    if (fopen_checked(outname, "r", mapfile_ptr)) {
      goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
    }
  } else {
    rewind(*mapfile_ptr);
  }
  marker_uidx = 0;
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    if (map_is_unsorted) {
      if (!fgets(g_textbuf, MAXLINELEN, *mapfile_ptr)) {
	goto ped_to_bed_multichar_allele_ret_READ_FAIL;
      }
    } else {
      if (get_next_noncomment_excl(marker_exclude, *mapfile_ptr, &bufptr, &line_idx, &marker_uidx)) {
	goto ped_to_bed_multichar_allele_ret_READ_FAIL;
      }
    }
    if (marker_allele_cts[4 * marker_idx + 2]) {
      uii = marker_allele_cts[4 * marker_idx + 3];
      // bugfix (28 May 2018): do NOT look up map_reverse[] in unsorted case
      sprintf(g_logbuf, "Warning: Variant %" PRIuPTR "%s %sallelic; setting rarest missing.\n", marker_idx + 1, map_is_unsorted? " (post-sort)" : "", (uii? "quad" : "tri"));
      logerrprintb();
      get_top_two_ui(&(marker_allele_cts[4 * marker_idx]), uii? 4 : 3, &ulii, &uljj);
    } else {
      ulii = (marker_allele_cts[4 * marker_idx] < marker_allele_cts[4 * marker_idx + 1])? 1 : 0;
      uljj = ulii ^ 1;
    }

    aptr1 = get_llstr((Ll_str*)(&(marker_alleles_tmp[marker_idx])), uljj);
    aptr2 = get_llstr((Ll_str*)(&(marker_alleles_tmp[marker_idx])), ulii);
    if (!aptr1) {
      aptr1 = missing_geno_ptr;
    } else {
      alen1 = strlen_se(aptr1);
      if (aptr1[alen1] == '\t') {
	aptr1[alen1] = '\0';
      }
    }
    if (!aptr2) {
      aptr2 = missing_geno_ptr;
    } else {
      alen2 = strlen_se(aptr2);
      if (aptr2[alen2] == '\t') {
	aptr2[alen2] = '\0';
      }
    }
    marker_allele_ptrs[2 * marker_idx] = aptr1;
    marker_allele_ptrs[2 * marker_idx + 1] = aptr2;
    if (map_is_unsorted) {
      bufptr = (char*)memchr(g_textbuf, '\n', MAXLINELEN);
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
        goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
      }
    } else {
      uii = strlen_se(bufptr);
      if (fwrite_checked(bufptr, uii, outfile)) {
        goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
      }
      putc_unlocked('\t', outfile);
      bufptr = skip_initial_spaces(&(bufptr[uii + 1]));
      bufptr = write_token(bufptr, outfile);
      if (cm_col_exists) {
        ucc = (unsigned char)(*bufptr);
	// should be good enough at detecting nonnumeric values...
	if (((ucc >= '0') && (ucc <= '9')) || (ucc == '-') || (ucc == '+')) {
	  bufptr = write_token_notab(bufptr, outfile);
	} else {
	  putc_unlocked('0', outfile);
	  bufptr = next_token(bufptr);
	}
      } else {
	putc_unlocked('0', outfile);
      }
      putc_unlocked('\t', outfile);
      uii = strlen_se(bufptr);
      fwrite(bufptr, 1, uii, outfile);
    }
    putc_unlocked('\t', outfile);
    fputs(aptr1, outfile);
    putc_unlocked('\t', outfile);
    fputs(aptr2, outfile);
    if (putc_checked('\n', outfile)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    marker_uidx++;
  }
  sample_ct4 = (sample_ct + 3) / 4;
  fclose_null(mapfile_ptr);
  if (map_is_unsorted) {
    unlink(outname);
  }
  fclose_null(outfile_ptr);
  if (bigstack_alloc_c(ped_buflen, &loadbuf)) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  if (bigstack_left() >= marker_ct * sample_ct4) {
    markers_per_pass = marker_ct;
    sprintf(g_logbuf, "Performing single-pass .bed write (%" PRIuPTR " variant%s, %" PRIuPTR " %s).\n", marker_ct, (marker_ct == 1)? "" : "s", sample_ct, species_str(sample_ct));
    pass_ct = 1;
  } else {
    if (!map_is_unsorted) {
      if (bigstack_alloc_ll(sample_ct, &line_starts)) {
	goto ped_to_bed_multichar_allele_ret_NOMEM;
      }
    }
    markers_per_pass = bigstack_left() / sample_ct4;
    if (!markers_per_pass) {
      goto ped_to_bed_multichar_allele_ret_NOMEM;
    }
    pass_ct = (marker_ct + markers_per_pass - 1) / markers_per_pass;
    sprintf(g_logbuf, "Performing %u-pass .bed write (%u/%" PRIuPTR " variant%s/pass, %" PRIuPTR " %s).\n", pass_ct, markers_per_pass, marker_ct, (markers_per_pass == 1)? "" : "s", sample_ct, species_str(sample_ct));
  }
  logprintb();
  writebuf = g_bigstack_base;
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(outname, FOPEN_WB, outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, *outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  rewind(*pedfile_ptr);
  umm = 0;
  for (uii = 0; uii < pass_ct; uii++) {
    marker_start = uii * markers_per_pass;
    if (uii + 1 == pass_ct) {
      ujj = marker_ct - marker_start;
      last_pass = 1;
    } else {
      ujj = markers_per_pass;
    }
    memset(writebuf, 0, ujj * sample_ct4);
    marker_end = marker_start + ujj;
    fputs("0%", stdout);
    sample_idx = 0;
    // 94 instead of 100 due to big fwrite at the end
    for (pct = 1; pct <= 94; pct++) {
      loop_end = (((uint64_t)pct) * sample_ct) / 94LLU;
      for (; sample_idx < loop_end; sample_idx++) {
	if ((!uii) || map_is_unsorted) {
	  do {
	    if (!last_pass) {
	      ped_next_thresh = ftello(*pedfile_ptr);
	    }
	    if (!fgets(loadbuf, ped_buflen, *pedfile_ptr)) {
	      goto ped_to_bed_multichar_allele_ret_READ_FAIL;
	    }
	    col1_ptr = skip_initial_spaces(loadbuf);
	  } while (is_eoln_or_comment_kns(*col1_ptr));
	  bufptr = next_token_mult(col1_ptr, ped_col_skip);
	} else {
	  ped_next_thresh = line_starts[sample_idx];
	  if (fseeko(*pedfile_ptr, line_starts[sample_idx], SEEK_SET)) {
	    goto ped_to_bed_multichar_allele_ret_READ_FAIL;
	  }
	  if (!fgets(loadbuf, ped_buflen, *pedfile_ptr)) {
	    goto ped_to_bed_multichar_allele_ret_READ_FAIL;
	  }
	  bufptr = loadbuf;
	}
	marker_idx = uii * markers_per_pass;
	ii_shift = (sample_idx % 4) * 2;
	wbufptr = &(writebuf[sample_idx / 4]);
	if (map_is_unsorted) {
	  umm = 0;
	  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	    aptr1 = bufptr;
	    // already validated
	    bufptr = token_endnn(bufptr);
	    alen1 = (uintptr_t)(bufptr - aptr1);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr1[alen1] = '\0';
	    aptr2 = bufptr;
	    bufptr = token_endnn(bufptr);
	    alen2 = (uintptr_t)(bufptr - aptr2);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr2[alen2] = '\0';
	    if (IS_SET(marker_exclude, marker_uidx)) {
	      continue;
	    }
	    ukk = map_reverse[umm++];
	    if ((ukk >= marker_start) && (ukk < marker_end)) {
	      ucc = 1;
	      if ((*aptr1 != missing_geno) || (alen1 != 1)) {
		if (!strcmp(aptr1, marker_allele_ptrs[2 * ukk + 1])) {
		  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		    ucc = 3;
		  } else if (!strcmp(aptr2, marker_allele_ptrs[2 * ukk])) {
		    ucc = 2;
		  }
		} else if (!strcmp(aptr1, marker_allele_ptrs[2 * ukk])) {
		  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		    ucc = 0;
		  } else if (!strcmp(aptr2, marker_allele_ptrs[2 * ukk + 1])) {
		    ucc = 2;
		  }
		}
	      }
	      wbufptr[(ukk - marker_start) * sample_ct4] |= ucc << ii_shift;
	      marker_idx++;
	    }
	  }
	} else {
	  for (marker_uidx = umm; marker_idx < marker_end; marker_uidx++) {
	    aptr1 = bufptr;
	    bufptr = token_endnn(bufptr);
	    alen1 = (uintptr_t)(bufptr - aptr1);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr1[alen1] = '\0';
	    aptr2 = bufptr;
	    bufptr = token_endnn(bufptr);
	    alen2 = (uintptr_t)(bufptr - aptr2);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr2[alen2] = '\0';
	    if (is_set(marker_exclude, marker_uidx)) {
	      continue;
	    }
	    ucc = 1;
	    if ((*aptr1 != missing_geno) || (alen1 != 1)) {
	      if (!strcmp(aptr1, marker_allele_ptrs[2 * marker_idx + 1])) {
		if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		  ucc = 3;
		} else if (!strcmp(aptr2, marker_allele_ptrs[2 * marker_idx])) {
		  ucc = 2;
		}
	      } else if (!strcmp(aptr1, marker_allele_ptrs[2 * marker_idx])) {
		if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		  ucc = 0;
		} else if (!strcmp(aptr2, marker_allele_ptrs[2 * marker_idx + 1])) {
		  ucc = 2;
		}
	      }
	    }
	    *wbufptr |= ucc << ii_shift;
	    wbufptr = &(wbufptr[sample_ct4]);
	    marker_idx++;
	  }
	  if (!last_pass) {
	    line_starts[sample_idx] = ped_next_thresh + (uintptr_t)(bufptr - loadbuf);
	  }
	}
      }
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
    if (fwrite_checked(writebuf, ujj * sample_ct4, *outfile_ptr)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    if (!last_pass) {
      printf("\rPass %u:    \b\b\b", uii + 2);
      fflush(stdout);
      if (map_is_unsorted) {
	rewind(*pedfile_ptr);
      } else {
	umm = marker_uidx;
      }
    }
  }

  while (0) {
  ped_to_bed_multichar_allele_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ped_to_bed_multichar_allele_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ped_to_bed_multichar_allele_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ped_to_bed_multichar_allele_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_6:
    g_bigstack_base -= cur_slen_rdup;
    logprint("\n");
    if (retval != RET_NOMEM) {
      LOGERRPRINTF("Error: More than 4 different alleles at variant %u%s.\n", uii + 1, map_is_unsorted? " (post-sort)" : "");
    }
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4:
    g_bigstack_base -= cur_slen_rdup;
    logprint("\n");
    LOGERRPRINTF("Error: Half-missing call in .ped file at variant %" PRIuPTR ", line %" PRIuPTR ".\n", marker_uidx + 1, line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  ped_to_bed_multichar_allele_ret_MISSING_TOKENS:
    logprint("\n");
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file has fewer tokens than expected.\n", line_idx);
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  bigstack_end_reset(bigstack_end_mark);
  // no marker_allele_ptrs free since all strings were allocated on top of
  // stack
  return retval;
}

int32_t ped_to_bed(char* pedname, char* mapname, char* outname, char* outname_end, uint32_t fam_cols, uint64_t misc_flags, int32_t missing_pheno, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* mapfile = nullptr;
  FILE* pedfile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t* marker_exclude;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t marker_ct = 0;
  uintptr_t sample_ct = 0;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t map_is_unsorted = 0;
  int32_t last_chrom = 0;
  uint32_t last_mpos = 0;
  uint32_t ped_buflen = 1;
  int32_t retval = 0;
  uint32_t ped_col_skip_iid_m1 = ((fam_cols & FAM_COL_34) / (FAM_COL_34 / 2)) + ((fam_cols & FAM_COL_5) / FAM_COL_5) + ((fam_cols & FAM_COL_6) / FAM_COL_6);
  uint32_t ped_col_skip = ped_col_skip_iid_m1 + 1 + ((fam_cols & FAM_COL_1) / FAM_COL_1);
  uint32_t last_pass = 0;
  int64_t* line_starts = nullptr;

  uint32_t is_single_char_alleles = 1;
  char missing_geno = *g_missing_geno_ptr;
  char missing_pheno_str[12];

  uint32_t pass_ct;
  uintptr_t sample_ct4;
  uint32_t pct;
  char* marker_alleles_f;
  char* marker_alleles;
  uint32_t* marker_allele_cts;
  uint32_t* map_reverse;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t line_idx;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uint32_t cm_col_exists;
  uint32_t markers_per_pass;
  uint32_t marker_start;
  uint32_t marker_end;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  int32_t jj;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uintptr_t unfiltered_marker_ct_limit;
  char* col1_ptr;
  char* col2_ptr;
  char* bufptr;
  char* bufptr2;
  char cc;
  char cc2;
  unsigned char ucc;
  uint32_t ii_shift;
  unsigned char* writebuf;
  unsigned char* wbufptr;
  int64_t ped_size;
  int64_t ped_next_thresh;
  {
    int32toa_x(missing_pheno, '\0', missing_pheno_str);
    marker_exclude = (uintptr_t*)g_bigstack_base;
    marker_exclude[0] = 0;
    // don't use fopen_checked() here, since we want to customize the error
    // message.
    mapfile = fopen(mapname, "r");
    if (!mapfile) {
      uii = strlen(mapname);
      if ((uii > 8) && ((!memcmp(&(mapname[uii - 8]), ".ped.map", 8)) || (!memcmp(&(mapname[uii - 8]), ".map.map", 8)))) {
	LOGERRPRINTFWW("Error: Failed to open %s. (--file expects a filename *prefix*; '.ped' and '.map' are automatically appended.)\n", mapname);
      } else {
	LOGERRPRINTFWW(g_errstr_fopen, mapname);
      }
      goto ped_to_bed_ret_OPEN_FAIL;
    }
    g_textbuf[MAXLINELEN - 6] = ' ';
    if (check_cm_col(mapfile, g_textbuf, 0, allow_no_variants, MAXLINELEN - 5, &cm_col_exists, &line_idx)) {
      if (line_idx) {
	goto ped_to_bed_ret_MISSING_TOKENS_MAP;
      } else {
	logerrprint("Error: Empty .map file.\n");
	goto ped_to_bed_ret_INVALID_FORMAT;
      }
    }
    if (!line_idx) {
      // no variants
      goto ped_to_bed_empty_map_with_allow_no_vars;
    }
    line_idx--;
    unfiltered_marker_ct_limit = bigstack_left();
    if (unfiltered_marker_ct_limit > 0xfffffff) {
      unfiltered_marker_ct_limit = 0x80000000U;
    } else {
      unfiltered_marker_ct_limit *= 8;
    }
    do {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 6]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .map file is pathologically long.\n", line_idx);
	goto ped_to_bed_ret_INVALID_FORMAT_2;
      }
      col1_ptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_or_comment_kns(*col1_ptr)) {
	continue;
      }
      char* col1_end = token_endnn(col1_ptr);
      col2_ptr = skip_initial_spaces(col1_end);
      bufptr = next_token_mult(col2_ptr, 1 + cm_col_exists);
      if (no_more_tokens_kns(bufptr)) {
	goto ped_to_bed_ret_MISSING_TOKENS_MAP;
      }
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code_destructive(".map file", line_idx, allow_extra_chroms, col1_ptr, col1_end, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto ped_to_bed_ret_1;
      }
      if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	SET_BIT(unfiltered_marker_ct, marker_exclude);
	marker_exclude_ct++;
      } else {
	if (scan_int_abs_defcap(bufptr, &jj)) {
	  sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of .map file.\n", line_idx);
	  goto ped_to_bed_ret_INVALID_FORMAT_2;
	}
	if (jj >= 0) {
	  if (!map_is_unsorted) {
	    if ((cur_chrom_code < last_chrom) || ((cur_chrom_code == last_chrom) && ((uint32_t)jj < last_mpos))) {
	      map_is_unsorted = 1;
	    }
	    last_chrom = cur_chrom_code;
	    last_mpos = (uint32_t)jj;
	  }
	  uii = strlen_se(col2_ptr) + 1;
	  if (uii > max_marker_id_len) {
	    max_marker_id_len = uii;
	  }
	} else {
	  SET_BIT(unfiltered_marker_ct, marker_exclude);
	  marker_exclude_ct++;
	}
      }
      unfiltered_marker_ct++;
      if (unfiltered_marker_ct > 0x7ffffffd) {
	logprint("Error: Too many variants in .map file (max 2147483645).\n");
	goto ped_to_bed_ret_INVALID_FORMAT;
      }
      if (!(unfiltered_marker_ct & (BITCT - 1))) {
	if (unfiltered_marker_ct == unfiltered_marker_ct_limit) {
	  goto ped_to_bed_ret_NOMEM;
	}
	marker_exclude[unfiltered_marker_ct / BITCT] = 0;
      }
    } while (fgets(g_textbuf, MAXLINELEN - 5, mapfile));
    if (!feof(mapfile)) {
      goto ped_to_bed_ret_READ_FAIL;
    }
    marker_ct = unfiltered_marker_ct - marker_exclude_ct;
    if ((!marker_ct) && (!allow_no_variants)) {
      logprint("Error: No variants in current analysis.\n");
      goto ped_to_bed_ret_ALL_MARKERS_EXCLUDED;
    }
   ped_to_bed_empty_map_with_allow_no_vars:
    bigstack_alloc_ul(BITCT_TO_WORDCT(unfiltered_marker_ct), &marker_exclude);

    if (map_is_unsorted) {
      retval = load_sort_and_write_map(&map_reverse, mapfile, 3 + cm_col_exists, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, 1, chrom_info_ptr);
      if (retval) {
	goto ped_to_bed_ret_1;
      }
      cm_col_exists = 1;
      fclose_null(&mapfile);
    }
    // provisionally assume max_marker_allele_blen == 2
    // bugfix: allocate this after map_reverse
    // quasi-bugfix (14 Nov 2017): need to zero-initialize marker_allele_cts
    //   for consistent allele ordering
    if (bigstack_alloc_c(marker_ct * 2, &marker_alleles_f) ||
	bigstack_calloc_c(marker_ct * 4, &marker_alleles) ||
	bigstack_calloc_ui(marker_ct * 4, &marker_allele_cts)) {
      goto ped_to_bed_ret_NOMEM;
    }

    // first .ped scan: count samples, write .fam, note alleles at each locus
    if (fopen_checked(pedname, FOPEN_RB, &pedfile)) {
      goto ped_to_bed_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".fam", 5);
    if (fopen_checked(outname, "w", &outfile)) {
      goto ped_to_bed_ret_OPEN_FAIL;
    }
    loadbuf = (char*)g_bigstack_base;
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto ped_to_bed_ret_NOMEM;
    }
    if (fseeko(pedfile, 0, SEEK_END)) {
      goto ped_to_bed_ret_READ_FAIL;
    }
    ped_size = ftello(pedfile);
    rewind(pedfile);
    logprint("Scanning .ped file...");
    fputs(" 0%", stdout);
    fflush(stdout);
    ped_next_thresh = ped_size / 100;
    loadbuf[loadbuf_size - 1] = ' ';
    pct = 0;
    line_idx = 0;
    while (fgets(loadbuf, loadbuf_size, pedfile)) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  logprint("\n");
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file is pathologically long.\n", line_idx);
	  goto ped_to_bed_ret_INVALID_FORMAT_2;
	} else {
	  goto ped_to_bed_ret_NOMEM;
	}
      }
      col1_ptr = skip_initial_spaces(loadbuf);
      if (is_eoln_or_comment_kns(*col1_ptr)) {
	ulii = strlen(loadbuf) + 1;
	if (ulii > ped_buflen) {
	  ped_buflen = ulii;
	}
	continue;
      }
      if (fam_cols & FAM_COL_1) {
	col2_ptr = next_token(col1_ptr);
      } else {
	col2_ptr = col1_ptr;
      }
      bufptr = next_token_multz(col2_ptr, ped_col_skip_iid_m1);
      if (no_more_tokens_kns(bufptr)) {
	goto ped_to_bed_ret_MISSING_TOKENS_PED;
      }
      bufptr = token_endnn(bufptr);
      if ((bufptr - col1_ptr) > (MAXLINELEN / 2) - 4) {
	logprint("\n");
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file has a pathologically long token.\n", line_idx);
	goto ped_to_bed_ret_INVALID_FORMAT_2;
      }
      if (fwrite_checked(col1_ptr, strlen_se(col1_ptr), outfile)) {
	goto ped_to_bed_ret_WRITE_FAIL;
      }
      putc_unlocked('\t', outfile);
      bufptr2 = write_token(col2_ptr, outfile);
      if (fam_cols & FAM_COL_34) {
	bufptr2 = write_token(bufptr2, outfile);
	bufptr2 = write_token(bufptr2, outfile);
      } else {
	fwrite("0\t0\t", 1, 4, outfile);
      }
      if (fam_cols & FAM_COL_5) {
	bufptr2 = write_token_notab(bufptr2, outfile);
      } else {
	putc_unlocked('0', outfile);
      }
      putc_unlocked('\t', outfile);
      if (fam_cols & FAM_COL_6) {
	fwrite(bufptr2, 1, strlen_se(bufptr2), outfile);
      } else {
	fputs(missing_pheno_str, outfile);
      }
      if (putc_checked('\n', outfile)) {
	goto ped_to_bed_ret_WRITE_FAIL;
      }
      marker_idx = 0;
      bufptr = skip_initial_spaces(bufptr);
      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	cc = *bufptr++;
        // bugfix (23 Nov 2018): must check for any eoln here, not just \0
	if (is_eoln_char(cc)) {
	  goto ped_to_bed_ret_MISSING_TOKENS_PED;
	}
	bufptr = skip_initial_spaces(bufptr);
	cc2 = *bufptr++;
	if (is_eoln_char(cc2)) {
	  goto ped_to_bed_ret_MISSING_TOKENS_PED;
	}
	bufptr = skip_initial_spaces(bufptr);
	if (IS_SET(marker_exclude, marker_uidx)) {
	  continue;
	}
	if (cc == missing_geno) {
	  if (cc2 != missing_geno) {
	    is_single_char_alleles = 0;
	    break;
	  }
	  marker_idx++;
	  continue;
	} else if (cc2 == missing_geno) {
	  is_single_char_alleles = 0;
	  break;
	}
	uii = 4 * (map_is_unsorted? map_reverse[marker_idx] : marker_idx);
	if (incr_text_allele0(cc, &(marker_alleles[uii]), &(marker_allele_cts[uii])) ||
	    incr_text_allele0(cc2, &(marker_alleles[uii]), &(marker_allele_cts[uii]))) {
	  is_single_char_alleles = 0;
	  break;
	}
	marker_idx++;
      }
      if ((!is_single_char_alleles) || (!is_eoln_kns(*bufptr))) {
	// either multi-character alleles, or invalid format.  Restart scan.
	putc_unlocked('\r', stdout);
	logstr("\n");
	if (!marker_ct) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file has more tokens than expected.\n", line_idx);
	  goto ped_to_bed_ret_INVALID_FORMAT_2;
	}
	logprint("Possibly irregular .ped line.  Restarting scan, assuming multichar alleles.\n");
	is_single_char_alleles = 0;
	break;
      }
      ulii = strlen(bufptr) + (uintptr_t)(bufptr - loadbuf) + 1;
      if (ulii > ped_buflen) {
	ped_buflen = ulii;
      }
      sample_ct++;
      if (ftello(pedfile) >= ped_next_thresh) {
	uii = (ftello(pedfile) * 100) / ped_size;
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", uii);
	fflush(stdout);
	pct = uii;
      }
    }
    if (is_single_char_alleles) {
      if (!feof(pedfile)) {
	goto ped_to_bed_ret_READ_FAIL;
      }
      if ((!sample_ct) && (!allow_no_samples)) {
	logprint("\n");
	sprintf(g_logbuf, "Error: No %s in .ped file.\n", g_species_plural);
	goto ped_to_bed_ret_INVALID_FORMAT_2;
      }
      if (fclose_null(&outfile)) {
	goto ped_to_bed_ret_WRITE_FAIL;
      }
      memcpy(outname_end, ".bim", 5);
      if (fopen_checked(outname, "w", &outfile)) {
	goto ped_to_bed_ret_OPEN_FAIL;
      }
      if (map_is_unsorted) {
	memcpy(outname_end, ".map.tmp", 9);
	if (fopen_checked(outname, "r", &mapfile)) {
	  goto ped_to_bed_ret_OPEN_FAIL;
	}
      } else {
	rewind(mapfile);
      }
      logstr(" done.\n");
      fputs("\r.ped scan complete (for binary autoconversion).\n", stdout);
      marker_uidx = 0;
      line_idx = 0;
      for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	if (map_is_unsorted) {
	  if (!fgets(g_textbuf, MAXLINELEN, mapfile)) {
	    goto ped_to_bed_ret_READ_FAIL;
	  }
	} else {
	  if (get_next_noncomment_excl(marker_exclude, mapfile, &bufptr, &line_idx, &marker_uidx)) {
	    goto ped_to_bed_ret_READ_FAIL;
	  }
	}
	if (marker_alleles[marker_idx * 4 + 2]) {
	  cc = marker_alleles[marker_idx * 4 + 3];
          sprintf(g_logbuf, "Warning: Variant %" PRIuPTR "%s %sallelic; setting rarest alleles missing.\n", marker_idx + 1, map_is_unsorted? " (post-sort)" : "", (cc? "quad" : "tri"));
	  logerrprintb();
	  ujj = (cc? 4 : 3);
	  // insertion sort
	  for (uii = 1; uii < ujj; uii++) {
	    ukk = marker_allele_cts[4 * marker_idx + uii];
	    if (marker_allele_cts[4 * marker_idx + uii - 1] < ukk) {
	      cc = marker_alleles[4 * marker_idx + uii];
	      umm = uii;
	      do {
		umm--;
		marker_alleles[4 * marker_idx + umm + 1] = marker_alleles[4 * marker_idx + umm];
		marker_allele_cts[4 * marker_idx + umm + 1] = marker_allele_cts[4 * marker_idx + umm];
	      } while (umm && (marker_allele_cts[4 * marker_idx + umm - 1] < ukk));
	      marker_alleles[4 * marker_idx + umm] = cc;
	      marker_allele_cts[4 * marker_idx + umm] = ukk;
	    }
	  }
	  cc = marker_alleles[marker_idx * 4 + 1];
	  cc2 = marker_alleles[marker_idx * 4];
	} else {
	  if (marker_allele_cts[marker_idx * 4] >= marker_allele_cts[marker_idx * 4 + 1]) {
	    cc = marker_alleles[marker_idx * 4 + 1];
	    cc2 = marker_alleles[marker_idx * 4];
	  } else {
	    cc = marker_alleles[marker_idx * 4];
	    cc2 = marker_alleles[marker_idx * 4 + 1];
	  }
	}
	marker_alleles_f[marker_idx * 2] = cc;
	marker_alleles_f[marker_idx * 2 + 1] = cc2;
	if (!cc) {
	  cc = '0';
	}
	if (!cc2) {
	  cc2 = '0';
	}
	if (map_is_unsorted) {
	  bufptr = (char*)memchr(g_textbuf, '\n', MAXLINELEN);
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto ped_to_bed_ret_WRITE_FAIL;
	  }
	} else {
	  bufptr = write_token(bufptr, outfile);
	  bufptr = write_token(bufptr, outfile);
	  if (cm_col_exists) {
	    ucc = (unsigned char)(*bufptr);
	    if (((ucc >= '0') && (ucc <= '9')) || (ucc == '-') || (ucc == '+')) {
	      bufptr = write_token_notab(bufptr, outfile);
	    } else {
	      putc_unlocked('0', outfile);
	      bufptr = next_token(bufptr);
	    }
	  } else {
	    putc_unlocked('0', outfile);
	  }
	  putc_unlocked('\t', outfile);
	  fwrite(bufptr, 1, strlen_se(bufptr), outfile);
	}
	putc_unlocked('\t', outfile);
	putc_unlocked(cc, outfile);
	putc_unlocked('\t', outfile);
	putc_unlocked(cc2, outfile);
	if (putc_checked('\n', outfile)) {
	  goto ped_to_bed_ret_WRITE_FAIL;
	}
	marker_uidx++;
      }
      sample_ct4 = (sample_ct + 3) / 4;
      bigstack_reset(marker_alleles);
      fclose_null(&mapfile);
      if (map_is_unsorted) {
	unlink(outname);
      }
      fclose_null(&outfile);
      if (bigstack_alloc_c(ped_buflen, &loadbuf)) {
	goto ped_to_bed_ret_NOMEM;
      }
      if (bigstack_left() >= marker_ct * sample_ct4) {
	markers_per_pass = marker_ct;
	sprintf(g_logbuf, "Performing single-pass .bed write (%" PRIuPTR " variant%s, %" PRIuPTR " %s).\n", marker_ct, (marker_ct == 1)? "" : "s", sample_ct, species_str(sample_ct));
	pass_ct = (marker_ct && sample_ct4)? 1 : 0;
      } else {
	if (!map_is_unsorted) {
	  if (bigstack_alloc_ll(sample_ct, &line_starts)) {
	    goto ped_to_bed_ret_NOMEM;
	  }
	}
	markers_per_pass = bigstack_left() / sample_ct4;
	if (!markers_per_pass) {
	  goto ped_to_bed_ret_NOMEM;
	}
	pass_ct = (marker_ct + markers_per_pass - 1) / markers_per_pass;
	sprintf(g_logbuf, "Performing %u-pass .bed write (%u/%" PRIuPTR " variant%s/pass, %" PRIuPTR " %s).\n", pass_ct, markers_per_pass, marker_ct, (markers_per_pass == 1)? "" : "s", sample_ct, species_str(sample_ct));
      }
      logprintb();
      writebuf = g_bigstack_base;
      memcpy(outname_end, ".bed", 5);
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto ped_to_bed_ret_OPEN_FAIL;
      }
      if (fwrite_checked("l\x1b\x01", 3, outfile)) {
	goto ped_to_bed_ret_WRITE_FAIL;
      }
      rewind(pedfile);
      umm = 0;
      for (uii = 0; uii < pass_ct; uii++) {
	marker_start = uii * markers_per_pass;
	if (uii + 1 == pass_ct) {
	  ujj = marker_ct - marker_start;
	  last_pass = 1;
	} else {
	  ujj = markers_per_pass;
	}
	memset(writebuf, 0, ujj * sample_ct4);
	marker_end = marker_start + ujj;
	fputs("0%", stdout);
	sample_idx = 0;
	// 94 instead of 100 due to big fwrite at the end
	for (pct = 1; pct <= 94; pct++) {
	  loop_end = (((uint64_t)pct) * sample_ct) / 94LLU;
	  for (; sample_idx < loop_end; sample_idx++) {
	    if ((!uii) || map_is_unsorted) {
	      do {
		if (!last_pass) {
		  ped_next_thresh = ftello(pedfile);
		}
		if (!fgets(loadbuf, ped_buflen, pedfile)) {
		  goto ped_to_bed_ret_READ_FAIL;
		}
		col1_ptr = skip_initial_spaces(loadbuf);
	      } while (is_eoln_or_comment_kns(*col1_ptr));
	      bufptr = next_token_mult(col1_ptr, ped_col_skip);
	    } else {
	      ped_next_thresh = line_starts[sample_idx];
	      if (fseeko(pedfile, line_starts[sample_idx], SEEK_SET)) {
		goto ped_to_bed_ret_READ_FAIL;
	      }
	      if (!fgets(loadbuf, ped_buflen, pedfile)) {
		goto ped_to_bed_ret_READ_FAIL;
	      }
	      bufptr = loadbuf;
	    }
	    marker_idx = uii * markers_per_pass;
	    ii_shift = (sample_idx % 4) * 2;
	    wbufptr = &(writebuf[sample_idx / 4]);
	    if (map_is_unsorted) {
	      // multipass optimizations are possible, but we won't bother,
	      // especially since the .map should rarely be unsorted in the
	      // first place...
	      umm = 0;
	      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
		cc = *bufptr++;
		bufptr = skip_initial_spaces(bufptr);
		cc2 = *bufptr++;
		bufptr = skip_initial_spaces(bufptr);
		if (IS_SET(marker_exclude, marker_uidx)) {
		  continue;
		}
		ukk = map_reverse[umm++];
		if ((ukk >= marker_start) && (ukk < marker_end)) {
		  ucc = 1;
		  if (cc == marker_alleles_f[2 * ukk + 1]) {
		    if (cc2 == cc) {
		      ucc = 3;
		    } else if (cc2 == marker_alleles_f[2 * ukk]) {
		      ucc = 2;
		    }
		  } else if (cc == marker_alleles_f[2 * ukk]) {
		    if (cc2 == cc) {
		      ucc = 0;
		    } else if (cc2 == marker_alleles_f[2 * ukk + 1]) {
		      ucc = 2;
		    }
		  }
		  wbufptr[(ukk - marker_start) * sample_ct4] |= ucc << ii_shift;
		  marker_idx++;
		}
	      }
	    } else {
	      for (marker_uidx = umm; marker_idx < marker_end; marker_uidx++) {
		cc = *bufptr++;
		bufptr = skip_initial_spaces(bufptr);
		cc2 = *bufptr++;
		bufptr = skip_initial_spaces(bufptr);
		if (IS_SET(marker_exclude, marker_uidx)) {
		  continue;
		}
		ucc = 1;
		if (cc == marker_alleles_f[2 * marker_idx + 1]) {
		  if (cc2 == cc) {
		    ucc = 3;
		  } else if (cc2 == marker_alleles_f[2 * marker_idx]) {
		    ucc = 2;
		  }
		} else if (cc == marker_alleles_f[2 * marker_idx]) {
		  if (cc2 == cc) {
		    ucc = 0;
		  } else if (cc2 == marker_alleles_f[2 * marker_idx + 1]) {
		    ucc = 2;
		  }
		}
		*wbufptr |= ucc << ii_shift;
		wbufptr = &(wbufptr[sample_ct4]);
		marker_idx++;
	      }
	      if (!last_pass) {
		line_starts[sample_idx] = ped_next_thresh + (uintptr_t)(bufptr - loadbuf);
	      }
	    }
	  }
	  if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
	  printf("\b\b%u%%", pct);
	  fflush(stdout);
	}
	if (fwrite_checked(writebuf, ujj * sample_ct4, outfile)) {
	  goto ped_to_bed_ret_WRITE_FAIL;
	}
	if (!last_pass) {
	  printf("\rPass %u:    \b\b\b", uii + 2);
	  fflush(stdout);
	  if (map_is_unsorted) {
	    rewind(pedfile);
	  } else {
	    umm = marker_uidx;
	  }
	}
      }
    } else {
      retval = ped_to_bed_multichar_allele(&pedfile, &outfile, outname, outname_end, &mapfile, unfiltered_marker_ct, marker_exclude, marker_ct, marker_alleles_f, map_is_unsorted, fam_cols, ped_col_skip_iid_m1 + 1, ped_col_skip, cm_col_exists, map_reverse, ped_size, missing_pheno_str);
      if (retval) {
	goto ped_to_bed_ret_1;
      }
    }

    if (fclose_null(&outfile)) {
      goto ped_to_bed_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    LOGPRINTFWW("--file: %s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);
  }

  while (0) {
  ped_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ped_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ped_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ped_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ped_to_bed_ret_MISSING_TOKENS_MAP:
    logprint("\n");
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .map file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  ped_to_bed_ret_MISSING_TOKENS_PED:
    logprint("\n");
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ped file has fewer tokens than expected.\n", line_idx);
  ped_to_bed_ret_INVALID_FORMAT_2:
    logerrprintb();
  ped_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  ped_to_bed_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 ped_to_bed_ret_1:
  fclose_cond(pedfile);
  fclose_cond(mapfile);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

uint32_t realpath_identical(const char* outname, const char* read_realpath, char* write_realpath_buf) {
#ifdef _WIN32
  const uint32_t fname_slen = GetFullPathName(outname, FNAMESIZE, write_realpath_buf, nullptr);
  return (fname_slen && (fname_slen <= FNAMESIZE) && (!strcmp(read_realpath, write_realpath_buf)));
#else
  return (realpath(outname, write_realpath_buf) && (!strcmp(read_realpath, write_realpath_buf)));
#endif
}

int32_t lgen_to_bed(char* lgenname, char* mapname, char* famname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, uint32_t lgen_modifier, char* lgen_reference_fname, Chrom_info* chrom_info_ptr) {
  // This code has not been carefully optimized, and also does not support
  // multipass writes.
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  uint32_t lgen_allele_count = lgen_modifier & LGEN_ALLELE_COUNT;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_vars = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t affection_01 = (misc_flags / MISC_AFFECTION_01) & 1;
  uint32_t map_cols = 3;
  uintptr_t* marker_exclude = nullptr;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t max_sample_id_len = 4;
  uintptr_t marker_ct = 0;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char** marker_allele_ptrs = nullptr;
  char* marker_ids = nullptr;
  uint32_t* marker_pos = nullptr;
  uintptr_t sample_ct = 0;
  char* sample_ids = nullptr;
  char* paternal_ids = nullptr;
  uintptr_t max_paternal_id_len = 2;
  char* maternal_ids = nullptr;
  uintptr_t max_maternal_id_len = 2;
  uintptr_t* sex_nm = nullptr;
  uintptr_t* sex_male = nullptr;
  uint32_t affection = 0;
  uintptr_t* founder_info = nullptr;
  uintptr_t* sample_exclude = nullptr;
  uint32_t map_is_unsorted = 0;
  uint32_t compound_genotypes = 1; // 0 = no, 1 = unresolved, 2 = yes
  char missing_geno = *missing_geno_ptr;
  uintptr_t* pheno_nm = nullptr;
  uintptr_t* pheno_c = nullptr;
  double* pheno_d = nullptr;
  char* sorted_marker_ids;
  uint32_t* marker_id_map;
  uint32_t* map_reverse;
  char* sorted_sample_ids;
  uint32_t* sample_id_map;
  unsigned char* writebuf;
  uintptr_t sample_ct4;
  uintptr_t marker_idx;
  unsigned char ucc;
  char* loadbuf;
  char* id_buf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* a1ptr;
  char* a2ptr;
  char* sptr;
  char* sptr2;
  int64_t lgen_size;
  int64_t lgen_next_thresh;
  uintptr_t loadbuf_size;
  uintptr_t line_idx;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uint32_t a1len;
  uint32_t a2len;
  uint32_t uii;
  uint32_t ujj;
  uint32_t pct;
  int32_t retval;
  int32_t ii;
  if (lgen_modifier == LGEN_ALLELE_COUNT) {
    logerrprint("Error: --allele-count must be used with --reference.\n");
    goto lgen_to_bed_ret_INVALID_CMDLINE;
  }

  retval = load_map(&infile, mapname, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &marker_exclude, &marker_ids, chrom_info_ptr, &marker_pos, &map_is_unsorted, allow_extra_chroms, allow_no_vars);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  retval = sort_item_ids(unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref, &sorted_marker_ids, &marker_id_map);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  if (map_is_unsorted) {
    // Writes a temporary .map which is read later, and then deleted.
    retval = load_sort_and_write_map(&map_reverse, infile, map_cols, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, 0, chrom_info_ptr);
    if (retval) {
      goto lgen_to_bed_ret_1;
    }
    for (uii = 0; uii < marker_ct; uii++) {
      marker_id_map[uii] = map_reverse[marker_id_map[uii]];
    }
  }
  // collapse
  if (bigstack_alloc_ui(unfiltered_marker_ct, &sample_id_map)) {
    goto lgen_to_bed_ret_NOMEM;
  }
  if (marker_ct) {
    fill_uidx_to_idx(marker_exclude, unfiltered_marker_ct, marker_ct, sample_id_map);
    for (uii = 0; uii < marker_ct; uii++) {
      marker_id_map[uii] = sample_id_map[marker_id_map[uii]];
    }
  }
  fclose_null(&infile);
  memcpy(marker_ids, sorted_marker_ids, marker_ct * max_marker_id_len);
  bigstack_reset(sorted_marker_ids);

  retval = load_fam(famname, FAM_COL_13456, 1, missing_pheno, affection_01, &sample_ct, &sample_ids, &max_sample_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &sample_exclude, allow_no_samples);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  retval = sort_item_ids_nx(&sorted_sample_ids, &sample_id_map, sample_ct, sample_ids, max_sample_id_len);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  if (bigstack_alloc_c(MAXV(max_marker_id_len, max_sample_id_len), &id_buf)) {
    goto lgen_to_bed_ret_NOMEM;
  }
  marker_allele_ptrs = (char**)bigstack_alloc(2 * marker_ct * sizeof(char*));
  if (!marker_allele_ptrs) {
    goto lgen_to_bed_ret_NOMEM;
  }
  uii = 2 * marker_ct;
  for (ujj = 0; ujj < uii; ujj++) {
    marker_allele_ptrs[ujj] = missing_geno_ptr;
  }
  sample_ct4 = (sample_ct + 3) / 4;
  if (bigstack_alloc_uc(((uintptr_t)marker_ct) * sample_ct4, &writebuf)) {
    logerrprint("Error: Multipass .lgen -> .bed autoconversions are not yet supported.  Try\nusing --chr and/or --memory (perhaps with a better machine).\n");
    goto lgen_to_bed_ret_CALC_NOT_YET_SUPPORTED;
  }
  if (sample_ct % 4) {
    ucc = 0x15 >> (6 - 2 * (sample_ct % 4));
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      memset(&(writebuf[marker_idx * sample_ct4]), 0x55, sample_ct4 - 1);
      writebuf[(marker_idx + 1) * sample_ct4 - 1] = ucc;
    }
  } else {
    memset(writebuf, 0x55, marker_ct * sample_ct4);
  }
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto lgen_to_bed_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  if (lgen_modifier & LGEN_REFERENCE) {
    if (fopen_checked(lgen_reference_fname, "r", &infile)) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    line_idx = 0;
    while (fgets(loadbuf, loadbuf_size, infile)) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ref file is pathologically long.\n", line_idx);
	  goto lgen_to_bed_ret_INVALID_FORMAT_2;
	}
	goto lgen_to_bed_ret_NOMEM;
      }
      cptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      cptr2 = token_endnn(cptr);
      a1ptr = skip_initial_spaces(cptr2);
      if (no_more_tokens_kns(a1ptr)) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .ref file has fewer tokens than expected.\n", line_idx);
	goto lgen_to_bed_ret_INVALID_FORMAT_2;
      }
      a1len = strlen_se(cptr);
      ii = bsearch_str(cptr, a1len, marker_ids, max_marker_id_len, marker_ct);
      if (ii != -1) {
	marker_idx = marker_id_map[(uint32_t)ii];
	if (marker_allele_ptrs[2 * marker_idx + 1] != missing_geno_ptr) {
	  cptr[a1len] = '\0';
	  LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in .ref file.\n", cptr);
	  goto lgen_to_bed_ret_INVALID_FORMAT_2;
	}
	sptr = token_endnn(a1ptr);
	a2ptr = skip_initial_spaces(sptr);
	a1len = (uintptr_t)(sptr - a1ptr);
	a1ptr[a1len] = '\0';
	if (allele_set(a1ptr, a1len, &(marker_allele_ptrs[2 * marker_idx + 1]))) {
	  goto lgen_to_bed_ret_NOMEM;
	}
	if (no_more_tokens_kns(a2ptr)) {
	  if (lgen_allele_count) {
	    a1ptr[a1len++] = 'v';
	    a1ptr[a1len] = '\0';
	    if (allele_set(a1ptr, a1len, &(marker_allele_ptrs[2 * marker_idx]))) {
	      goto lgen_to_bed_ret_NOMEM;
	    }
	  }
	} else {
	  a2len = strlen_se(a2ptr);
	  a2ptr[a2len] = '\0';
	  if (allele_set(a2ptr, a2len, &(marker_allele_ptrs[2 * marker_idx]))) {
	    goto lgen_to_bed_ret_NOMEM;
	  }
	}
	memset(&(writebuf[marker_idx * sample_ct4]), 0xff, sample_ct / 4);
	if (sample_ct % 4) {
	  writebuf[(marker_idx + 1) * sample_ct4 - 1] = 0x3f >> (6 - 2 * (sample_ct % 4));
	}
      }
    }
    if (!feof(infile)) {
      goto lgen_to_bed_ret_READ_FAIL;
    }
    fclose_null(&infile);
  }
  // PLINK 1.07 reports an error whenever there are 3+ alleles at one locus, so
  // backwards compatibility does not mandate that we worry about that case.
  // Thus we just use the obvious one-pass load, and save proper handling of
  // triallelic sites, etc. for the future .pgen engine.
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(outname, FOPEN_WB, &outfile)) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  if (fopen_checked(lgenname, "r", &infile)) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  if (fseeko(infile, 0, SEEK_END)) {
    goto lgen_to_bed_ret_READ_FAIL;
  }
  lgen_size = ftello(infile);
  rewind(infile);
  logprint("Processing .lgen file... ");
  fputs("0%", stdout);
  fflush(stdout);
  lgen_next_thresh = lgen_size / 100;
  pct = 0;
  if (!lgen_allele_count) {
    line_idx = 0;
    while (fgets(loadbuf, loadbuf_size, infile)) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  goto lgen_to_bed_ret_LONG_LINE;
	}
	goto lgen_to_bed_ret_NOMEM;
      }
      cptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      if (bsearch_read_fam_indiv(cptr, sorted_sample_ids, max_sample_id_len, sample_ct, &cptr3, &ii, id_buf) || is_eoln_kns(*cptr3)) {
	goto lgen_to_bed_ret_MISSING_TOKENS;
      }
      if (ii == -1) {
	goto lgen_to_bed_ret_MISSING_IID;
      }
      sample_idx = sample_id_map[(uint32_t)ii];
      cptr4 = token_endnn(cptr3);
      a1ptr = skip_initial_spaces(cptr4);
      if (is_eoln_kns(*a1ptr)) {
	goto lgen_to_bed_ret_MISSING_TOKENS;
      }
      sptr = token_endnn(a1ptr);
      a2ptr = skip_initial_spaces(sptr);
      if (compound_genotypes == 1) {
	if (is_eoln_kns(*a2ptr)) {
	  compound_genotypes = 2;
	} else {
	  compound_genotypes = 0;
	}
      }
      if (!compound_genotypes) {
	if (is_eoln_kns(*a2ptr)) {
	  goto lgen_to_bed_ret_MISSING_TOKENS;
	}
        a1len = (uintptr_t)(sptr - a1ptr);
	a2len = strlen_se(a2ptr);
      } else {
	if ((uintptr_t)(sptr - a1ptr) != 2) {
	  sprintf(g_logbuf, "Error: Invalid compound genotype on line %" PRIuPTR " of .lgen file.\n", line_idx);
	  goto lgen_to_bed_ret_INVALID_FORMAT_2;
	}
	a1len = 1;
	a2ptr = &(a1ptr[2]);
	*a2ptr = a1ptr[1];
	a2len = 1;
      }
      ii = bsearch_str(cptr3, (uintptr_t)(cptr4 - cptr3), marker_ids, max_marker_id_len, marker_ct);
      if (ii != -1) {
	marker_idx = marker_id_map[(uint32_t)ii];
	sptr = marker_allele_ptrs[2 * marker_idx + 1]; // existing A2
	a1ptr[a1len] = '\0';
	a2ptr[a2len] = '\0';
	if ((*a1ptr == missing_geno) && (a1len == 1)) {
	  if ((*a2ptr == missing_geno) && (a2len == 1)) {
	    uii = 1;
	  } else {
	    goto lgen_to_bed_ret_HALF_MISSING;
	  }
	} else if ((*a2ptr == missing_geno) && (a2len == 1)) {
	  goto lgen_to_bed_ret_HALF_MISSING;
        } else {
          if (sptr == missing_geno_ptr) {
	    if (allele_set(a1ptr, a1len, &(marker_allele_ptrs[2 * marker_idx + 1]))) {
	      goto lgen_to_bed_ret_NOMEM;
	    }
	    if (!strcmp(a1ptr, a2ptr)) {
	      uii = 2;
	    } else {
	      uii = 1;
	      if (allele_set(a2ptr, a2len, &(marker_allele_ptrs[2 * marker_idx]))) {
		goto lgen_to_bed_ret_NOMEM;
	      }
	    }
	  } else {
	    sptr2 = marker_allele_ptrs[2 * marker_idx];
	    if (sptr2 == missing_geno_ptr) {
	      if (!strcmp(a1ptr, sptr)) {
		if (!strcmp(a2ptr, sptr)) {
		  uii = 2;
		} else {
		  uii = 1;
		  if (allele_set(a2ptr, a2len, &(marker_allele_ptrs[2 * marker_idx]))) {
		    goto lgen_to_bed_ret_NOMEM;
		  }
		}
	      } else {
		if (allele_set(a1ptr, a1len, &(marker_allele_ptrs[2 * marker_idx]))) {
		  goto lgen_to_bed_ret_NOMEM;
		}
		if (!strcmp(a2ptr, sptr)) {
		  uii = 1;
		} else if (!strcmp(a2ptr, a1ptr)) {
		  uii = 0;
		} else {
		  printf("\nfail 1\n");
		  goto lgen_to_bed_ret_NOT_BIALLELIC;
		}
	      }
	    } else {
	      if (!strcmp(a1ptr, sptr)) {
		uii = 1;
	      } else if (!strcmp(a1ptr, sptr2)) {
		uii = 0;
	      } else {
		goto lgen_to_bed_ret_NOT_BIALLELIC;
	      }
	      if (!strcmp(a2ptr, sptr)) {
		uii++;
	      } else if (strcmp(a2ptr, sptr2)) {
		goto lgen_to_bed_ret_NOT_BIALLELIC;
	      }
	    }
	  }
	  if (uii) {
	    uii++;
	  }
	}
	ulii = marker_idx * sample_ct4 + (sample_idx / 4);
	ujj = (sample_idx % 4) * 2;
	writebuf[ulii] = (writebuf[ulii] & (~(3 << ujj))) | (uii << ujj);
      }
      if (ftello(infile) >= lgen_next_thresh) {
	uii = (ftello(infile) * 100) / lgen_size;
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", uii);
	fflush(stdout);
	pct = uii;
	lgen_next_thresh = ((pct + 1) * lgen_size) / 100;
      }
    }
  } else {
    line_idx = 0;
    while (fgets(loadbuf, loadbuf_size, infile)) {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  goto lgen_to_bed_ret_LONG_LINE;
	}
	goto lgen_to_bed_ret_NOMEM;
      }
      cptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      if (bsearch_read_fam_indiv(cptr, sorted_sample_ids, max_sample_id_len, sample_ct, &cptr3, &ii, id_buf) || is_eoln_kns(*cptr3)) {
	goto lgen_to_bed_ret_MISSING_TOKENS;
      }
      if (ii == -1) {
	goto lgen_to_bed_ret_MISSING_IID;
      }
      sample_idx = sample_id_map[(uint32_t)ii];
      cptr4 = token_endnn(cptr3);
      a1ptr = skip_initial_spaces(cptr4);
      ucc = *a1ptr;
      if (is_eoln_kns(ucc)) {
	goto lgen_to_bed_ret_MISSING_TOKENS;
      }
      ii = bsearch_str(cptr3, (uintptr_t)(cptr4 - cptr3), marker_ids, max_marker_id_len, marker_ct);
      if (ii != -1) {
	marker_idx = marker_id_map[(uint32_t)ii];
	a1len = strlen_se(a1ptr);
	uii = ((uint32_t)ucc) - 48;
	if ((a1len != 1) || (uii > 2)) {
	  uii = 1;
	} else if (uii) {
	  uii++;
	}
	ulii = marker_idx * sample_ct4 + (sample_idx / 4);
	ujj = (sample_idx % 4) * 2;
	writebuf[ulii] = (writebuf[ulii] & (~(3 << ujj))) | (uii << ujj);
      }
      if (ftello(infile) >= lgen_next_thresh) {
	uii = (ftello(infile) * 100) / lgen_size;
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", uii);
	fflush(stdout);
	pct = uii;
	lgen_next_thresh = ((pct + 1) * lgen_size) / 100;
      }
    }
  }
  if (!feof(infile)) {
    goto lgen_to_bed_ret_READ_FAIL;
  }
  fclose_null(&infile);
  if (pct < 10) {
    fputs("\b\b", stdout);
  } else if (pct < 100) {
    fputs("\b\b\b", stdout);
  } else {
    fputs("\b\b\b\b", stdout);
  }
  logprint("done.\n");
  for (uii = 0; uii < marker_ct; uii++) {
    if (popcount_chars((uintptr_t*)writebuf, uii * sample_ct4, (uii + 1) * sample_ct4) < sample_ct) {
      reverse_loadbuf(sample_ct, &(writebuf[uii * sample_ct4]));
      cptr = marker_allele_ptrs[uii * 2];
      marker_allele_ptrs[uii * 2] = marker_allele_ptrs[uii * 2 + 1];
      marker_allele_ptrs[uii * 2 + 1] = cptr;
    }
  }
  if (fwrite_checked(writebuf, ((uintptr_t)marker_ct) * sample_ct4, outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  if (map_is_unsorted) {
    memcpy(outname_end, ".map.tmp", 9);
    if (fopen_checked(outname, "r", &infile)) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
  } else {
    if (fopen_checked(mapname, "r", &infile)) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  uii = 0;
  marker_idx = 0;
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    if (is_eoln_or_comment_kns(*(skip_initial_spaces(g_textbuf)))) {
      continue;
    }
    if (IS_SET(marker_exclude, uii)) {
      uii++;
      continue;
    }
    cptr = (char*)memchr(g_textbuf, 0, MAXLINELEN);
    if (cptr[-1] == '\n') {
      cptr--;
      if (cptr[-1] == '\r') {
	cptr--;
      }
    }
    *cptr++ = '\t';
    fwrite(g_textbuf, 1, cptr - g_textbuf, outfile);
    fputs(marker_allele_ptrs[marker_idx * 2], outfile);
    putc_unlocked('\t', outfile);
    fputs(marker_allele_ptrs[marker_idx * 2 + 1], outfile);
    if (putc_checked('\n', outfile)) {
      goto lgen_to_bed_ret_WRITE_FAIL;
    }
    uii++;
    marker_idx++;
  }
  if (!feof(infile)) {
    goto lgen_to_bed_ret_READ_FAIL;
  }
  fclose_null(&infile);
  if (map_is_unsorted) {
    memcpy(outname_end, ".map.tmp", 9);
    unlink(outname);
  }
  if (fclose_null(&outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
#ifdef _WIN32
  uii = GetFullPathName(famname, FNAMESIZE, g_textbuf, nullptr);
  if ((!uii) || (uii > FNAMESIZE))
#else
  if (!realpath(famname, g_textbuf))
#endif
  {
    LOGERRPRINTFWW("Error: Failed to open %s.\n", outname);
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  // bugfix (25 Jul 2017): forgot the not
  if (!realpath_identical(outname, g_textbuf, &(g_textbuf[FNAMESIZE + 64]))) {
    if (fopen_checked(famname, "r", &infile)) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      cptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      ulii = strlen(cptr);
      if (cptr[ulii - 1] != '\n') {
	cptr[ulii++] = '\n';
      }
      if (fwrite_checked(cptr, ulii, outfile)) {
	goto lgen_to_bed_ret_WRITE_FAIL;
      }
    }
    if (!feof(infile)) {
      goto lgen_to_bed_ret_READ_FAIL;
    }
  }
  *outname_end = '\0';
  LOGPRINTFWW("--lfile: %s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);

  while (0) {
  lgen_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  lgen_to_bed_ret_CALC_NOT_YET_SUPPORTED:
    retval = RET_CALC_NOT_YET_SUPPORTED;
    break;
  lgen_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  lgen_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  lgen_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  lgen_to_bed_ret_LONG_LINE:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .lgen file is pathologically long.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_HALF_MISSING:
    LOGERRPRINTF("Error: Half-missing genotype on line %" PRIuPTR " of .lgen file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_MISSING_IID:
    LOGERRPRINTF("Error: Sample ID on line %" PRIuPTR " of .lgen file is missing from .fam file.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .lgen file has fewer tokens than expected.\n", line_idx);
  lgen_to_bed_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_NOT_BIALLELIC:
    *cptr4 = '\0';
    LOGERRPRINTFWW("Error: Variant '%s' in .lgen file has 3+ different alleles.\n", cptr3);
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 lgen_to_bed_ret_1:
  cleanup_allele_storage(2, marker_ct * 2, marker_allele_ptrs);
  bigstack_reset(bigstack_mark);
  aligned_free_cond(pheno_c);
  if (infile) {
    fclose(infile);
  }
  if (outfile) {
    fclose(outfile);
  }
  return retval;
}

static inline uint32_t update_tped_alleles_and_cts(uint32_t* allele_tot_ptr, char** alleles, uint32_t* alens, uint32_t* allele_cts, char* ss, uint32_t slen, uint32_t* allele_idx_ptr) {
  uint32_t allele_idx;
  for (allele_idx = 0; allele_idx < (*allele_tot_ptr); allele_idx++) {
    if ((slen == alens[allele_idx]) && (!memcmp(alleles[allele_idx], ss, slen))) {
      allele_cts[allele_idx] += 1;
      *allele_idx_ptr = allele_idx;
      return 0;
    }
  }
  *allele_idx_ptr = allele_idx;
  if (allele_idx < 4) {
    alens[allele_idx] = slen;
    allele_cts[allele_idx] = 1;
    *allele_tot_ptr = allele_idx + 1;
    if (slen < 2) {
      alleles[allele_idx] = (char*)(&(g_one_char_strs[((unsigned char)(*ss)) * 2]));
    } else {
      alleles[allele_idx] = (char*)malloc(slen + 1);
      if (!alleles[allele_idx]) {
        return RET_NOMEM;
      }
      memcpyx(alleles[allele_idx], ss, slen, '\0');
    }
    return 0;
  } else {
    return RET_INVALID_FORMAT;
  }
}

void transposed_to_bed_print_pct(uint32_t pct) {
  printf("Processing .tped file... %u%%", pct);
  fflush(stdout);
}

int32_t transposed_to_bed(char* tpedname, char* tfamname, char* outname, char* outname_end, uint64_t misc_flags, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* infile = nullptr;
  FILE* bimfile = nullptr;
  FILE* outfile = nullptr;
  char** marker_allele_ptrs = nullptr;
  uintptr_t sample_ct = 0;
  uintptr_t line_idx = 0;
  uint32_t no_extra_cols = 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  int32_t retval = 0;
  uint32_t pct = 0;
  uint32_t map_is_unsorted = 0;
  int64_t last_mapval = 0;
  uint32_t allele_tot = 0;
  uintptr_t marker_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t max_marker_allele_blen = 2; // for .bim.tmp reloading
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;

  // We include handling of nonstandard chromosome names here since our own
  // pipeline requires this.  (To set them to zero, --make-bed must be used
  // in addition to this autoconverter.)
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t orig_zec = chrom_info_ptr->zero_extra_chroms;

  // It's unnecessary for now, but we choose to include nice handling of
  // triallelic and quadallelic sites in this routine.
  char* alleles[4];
  char* salleles[4];
  uint32_t alens[4];
  uint32_t allele_cts[4];
  unsigned char writemap[17];
  uintptr_t max_markers;
  uintptr_t sample_ct4;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  int32_t ii;
  int32_t jj;
  unsigned char* writebuf;
  unsigned char* prewritebuf;
  unsigned char ucc;
  uint32_t loadbuf_size;
  char* allele_buf;
  char* loadbuf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* axptr;
  uint32_t* chrom_start;
  uint32_t* chrom_id;
  uint32_t axlen;
  unsigned char* ucptr;
  unsigned char* ucptr2;
  int64_t tped_size;
  int64_t tped_next_thresh;
  int64_t cur_mapval;
  int64_t* mapvals;
  uintptr_t marker_idx;
  uintptr_t marker_uidx;
  int64_t* ll_buf;
  uint32_t* pos_buf;
  char* marker_ids;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  double* marker_cms;
  {
    if (bigstack_alloc_ui(MAX_POSSIBLE_CHROM + 1, &chrom_start) ||
	bigstack_alloc_ui(MAX_POSSIBLE_CHROM, &chrom_id)) {
      goto transposed_to_bed_ret_NOMEM;
    }

    if (fopen_checked(tfamname, "r", &infile)) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".fam", 5);
    if (fopen_checked(outname, "w", &outfile)) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    g_textbuf[MAXLINELEN - 1] = ' ';
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tfam file is pathologically long.\n", line_idx);
	goto transposed_to_bed_ret_INVALID_FORMAT_2R;
      }
      cptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      ulii = strlen(cptr);
      if (cptr[ulii - 1] != '\n') {
	cptr[ulii++] = '\n';
      }
      if (fwrite_checked(cptr, ulii, outfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      sample_ct++;
    }
    if (!feof(infile)) {
      goto transposed_to_bed_ret_READ_FAIL;
    }
    if ((!sample_ct) && (!allow_no_samples)) {
      sprintf(g_logbuf, "Error: No %s in .tfam file.\n", g_species_plural);
      goto transposed_to_bed_ret_INVALID_FORMAT_2R;
    }
    sample_ct4 = (sample_ct + 3) / 4;
    fclose_null(&infile);
    fclose_null(&outfile);

    memcpy(outname_end, ".bim.tmp", 9);
    if (fopen_checked(outname, "w", &bimfile)) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed.tmp", 9);
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    if (bigstack_alloc_uc(sample_ct4, &writebuf) ||
	bigstack_alloc_uc(sample_ct, &prewritebuf)) {
      goto transposed_to_bed_ret_NOMEM;
    }
    if (bigstack_end_alloc_c(NON_BIGSTACK_MIN, &allele_buf)) {
      goto transposed_to_bed_ret_NOMEM;
    }
    max_markers = bigstack_left() / sizeof(int64_t);
    mapvals = (int64_t*)g_bigstack_base;
    writemap[16] = 1;
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }

    if (fopen_checked(tpedname, "r", &infile)) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    if (fseeko(infile, 0, SEEK_END)) {
      goto transposed_to_bed_ret_READ_FAIL;
    }
    logstr("Processing .tped file.\n");
    transposed_to_bed_print_pct(0);
    fflush(stdout);
    tped_size = ftello(infile);
    rewind(infile);
    tped_next_thresh = tped_size / 100;

    line_idx = 0;
    while (1) {
      line_idx++;
      g_textbuf[MAXLINELEN - 1] = ' ';
      if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	break;
      }
      // assume first four fields are within MAXLINELEN characters, but after
      // that, anything goes.  given e.g. 6MB indels in real datasets, there's
      // legitimate reason for a .tped line to be even longer than 2GB, so we
      // use a custom loading loop.
      char* textbuf_first_token = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*textbuf_first_token)) {
	if (!g_textbuf[MAXLINELEN - 1]) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has excessive whitespace.\n", line_idx);
	  goto transposed_to_bed_ret_INVALID_FORMAT_2R;
	}
	continue;
      }
      char* first_token_end = token_endnn(textbuf_first_token);
      const uint32_t chrom_name_slen = (uintptr_t)(first_token_end - textbuf_first_token);
      cptr2 = skip_initial_spaces(first_token_end);
      cptr3 = next_token_mult(cptr2, 2);
      cptr4 = next_token(cptr3);
      if (no_more_tokens_kns(cptr4)) {
	if (!g_textbuf[MAXLINELEN - 1]) {
	  if (chrom_name_slen > MAX_ID_SLEN) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has an excessively long\nchromosome/contig name.  (The " PROG_NAME_CAPS " limit is " MAX_ID_SLEN_STR " characters.)\n", line_idx);
	  } else if (cptr2 && (strlen_se(cptr2) > MAX_ID_SLEN)) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has an excessively long variant ID.\n(The " PROG_NAME_CAPS " limit is " MAX_ID_SLEN_STR " characters.)\n", line_idx);
	  } else if (next_token(cptr2) && (strlen_se(next_token(cptr2)) > MAX_ID_SLEN)) {
	    // far higher bound than necessary; main point is to ensure that if
	    // we fall through to the "excessive whitespace" error message,
	    // that complaint is justified.
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has an excessively long centimorgan\nposition.\n", line_idx);
	  } else if (cptr3 && (strlen_se(cptr3) > MAX_ID_SLEN)) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has an excessively long bp coordinate.\n", line_idx);
	  } else {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has excessive whitespace.\n", line_idx);
	  }
	  goto transposed_to_bed_ret_INVALID_FORMAT_2R;
	} else {
	  goto transposed_to_bed_ret_MISSING_TOKENS;
	}
      }
      if (ftello(infile) >= tped_next_thresh) {
	uii = (ftello(infile) * 100) / tped_size;
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", uii);
	fflush(stdout);
	pct = uii;
	tped_next_thresh = ((pct + 1) * tped_size) / 100;
      }
      *first_token_end = '\0';
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code(textbuf_first_token, ".tped file", line_idx, chrom_name_slen, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto transposed_to_bed_ret_1;
      }

      if (scan_int_abs_defcap(cptr3, &jj)) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has an invalid bp coordinate.\n", line_idx);
	goto transposed_to_bed_ret_INVALID_FORMAT_2R;
      }
      char* textbuf_iter = textbuf_first_token;
      if ((!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) || (jj < 0)) {
	cptr2 = cptr4;
	goto transposed_to_bed_nextline;
      }
      uii = strlen_se(cptr2);
      if (uii >= max_marker_id_len) {
	max_marker_id_len = uii + 1;
      }
      cur_mapval = (int64_t)((((uint64_t)((uint32_t)cur_chrom_code)) << 32) | ((uint32_t)jj));
      if (marker_ct == max_markers) {
	goto transposed_to_bed_ret_NOMEM;
      }
      mapvals[marker_ct++] = cur_mapval;
      if (last_mapval > cur_mapval) {
	map_is_unsorted = 1;
      } else {
	last_mapval = cur_mapval;
      }
      for (uii = 0; uii < 3; uii++) {
	char* token_end = token_endnn(textbuf_iter);
	*token_end++ = '\t';
	fwrite(textbuf_iter, 1, token_end - textbuf_iter, bimfile);
	textbuf_iter = skip_initial_spaces(token_end);
      }
      cptr2 = token_endnn(textbuf_iter);
      *cptr2++ = '\t';
      if (fwrite_checked(textbuf_iter, cptr2 - textbuf_iter, bimfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      cptr2 = cptr4;
      alleles[0] = nullptr;
      alleles[1] = nullptr;
      alleles[2] = nullptr;
      alleles[3] = nullptr;
      fill_uint_zero(4, allele_cts);
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	cptr2 = skip_initial_spaces(cptr2);
	while (cptr2 == &(g_textbuf[MAXLINELEN - 1])) {
	  if (cptr2[-1] == '\n') {
	    goto transposed_to_bed_ret_MISSING_TOKENS;
	  }
	  if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	    if (ferror(infile)) {
	      goto transposed_to_bed_ret_READ_FAIL;
	    }
	    goto transposed_to_bed_ret_MISSING_TOKENS;
	  }
	  cptr2 = skip_initial_spaces(g_textbuf);
	}
	axptr = cptr2;
	axlen = strlen_se(cptr2);
	if (!axlen) {
	  goto transposed_to_bed_ret_MISSING_TOKENS;
	}
	cptr2 = &(axptr[axlen]);
	// only way for this to happen if it isn't at end of buffer is if we're
	// at EOF, which is an error anyway
	if (!(*cptr2)) {
	  cptr3 = memcpya(allele_buf, axptr, axlen);
	  axptr = allele_buf;
	  do {
	    if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	      if (ferror(infile)) {
		goto transposed_to_bed_ret_READ_FAIL;
	      }
	      goto transposed_to_bed_ret_MISSING_TOKENS;
	    }
	    cptr2 = g_textbuf;
	    if (!is_space_or_eoln(*cptr2)) {
	      cptr2 = token_endnn(cptr2);
	    }
	    if ((((uintptr_t)(cptr3 - allele_buf)) + ((uintptr_t)(cptr2 - g_textbuf))) >= NON_BIGSTACK_MIN) {
	      goto transposed_to_bed_ret_NOMEM;
	    }
	    cptr3 = memcpya(cptr3, g_textbuf, cptr2 - g_textbuf);
	  } while (!(*cptr2));
	  axlen = (uintptr_t)(cptr3 - allele_buf);
	}
	if ((*axptr != missing_geno) || (axlen != 1)) {
	  retval = update_tped_alleles_and_cts(&allele_tot, alleles, alens, allele_cts, axptr, axlen, &uii);
	  if (retval) {
	    if (retval == RET_INVALID_FORMAT) {
	      goto transposed_to_bed_ret_TOO_MANY_ALLELES;
	    }
	    goto transposed_to_bed_ret_NOMEM;
	  }
	} else {
	  uii = 4;
	}
	cptr2 = skip_initial_spaces(cptr2);
	while (cptr2 == &(g_textbuf[MAXLINELEN - 1])) {
	  if (cptr2[-1] == '\n') {
	    goto transposed_to_bed_ret_MISSING_TOKENS;
	  }
	  if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	    if (ferror(infile)) {
	      goto transposed_to_bed_ret_READ_FAIL;
	    }
	    goto transposed_to_bed_ret_MISSING_TOKENS;
	  }
	  cptr2 = skip_initial_spaces(g_textbuf);
	}
	axptr = cptr2;
	axlen = strlen_se(cptr2);
	cptr2 = &(axptr[axlen]);
	if (!(*cptr2)) {
	  if (!axlen) {
	    goto transposed_to_bed_ret_MISSING_TOKENS;
	  }
	  cptr3 = memcpya(allele_buf, axptr, axlen);
	  axptr = allele_buf;
	  do {
	    cptr2 = g_textbuf;
	    if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	      if (ferror(infile)) {
		goto transposed_to_bed_ret_READ_FAIL;
	      } else if (sample_idx != sample_ct - 1) {
		goto transposed_to_bed_ret_MISSING_TOKENS;
	      } else {
		g_textbuf[0] = '\0';
		break;
	      }
	    }
	    if (!is_space_or_eoln(*cptr2)) {
	      cptr2 = token_endnn(cptr2);
	    }
	    if ((((uintptr_t)(cptr3 - allele_buf)) + ((uintptr_t)(cptr2 - g_textbuf))) >= NON_BIGSTACK_MIN) {
	      goto transposed_to_bed_ret_NOMEM;
	    }
	    cptr3 = memcpya(cptr3, g_textbuf, cptr2 - g_textbuf);
	  } while (!(*cptr2));
	  axlen = (uintptr_t)(cptr3 - allele_buf);
	}
	if ((*axptr != missing_geno) || (axlen != 1)) {
	  if (uii == 4) {
	    goto transposed_to_bed_ret_HALF_MISSING;
	  }
	  retval = update_tped_alleles_and_cts(&allele_tot, alleles, alens, allele_cts, axptr, axlen, &ujj);
	  if (retval) {
	    if (retval == RET_INVALID_FORMAT) {
	      goto transposed_to_bed_ret_TOO_MANY_ALLELES;
	    }
	    goto transposed_to_bed_ret_NOMEM;
	  }
	  prewritebuf[sample_idx] = uii * 4 + ujj;
	} else {
	  if (uii != 4) {
	    goto transposed_to_bed_ret_HALF_MISSING;
	  }
	  prewritebuf[sample_idx] = 16;
	}
      }

      memcpy(salleles, alleles, 4 * sizeof(intptr_t));
      for (uii = 1; uii < 4; uii++) {
	ujj = allele_cts[uii];
	if (allele_cts[uii - 1] < ujj) {
	  axptr = salleles[uii];
	  ii = uii;
	  do {
	    ii--;
	    salleles[((uint32_t)ii) + 1] = salleles[(uint32_t)ii];
	    allele_cts[((uint32_t)ii) + 1] = allele_cts[(uint32_t)ii];
	  } while (ii && (allele_cts[((uint32_t)ii) - 1] < ujj));
	  salleles[(uint32_t)ii] = axptr;
	  allele_cts[(uint32_t)ii] = ujj;
	}
      }
      if (allele_cts[2]) {
	putc_unlocked('\r', stdout);
	LOGPRINTF("Note: Variant %" PRIuPTR " is %sallelic.  Setting rarest alleles to missing.\n", marker_ct - 1, allele_cts[3]? "quad" : "tri");
	transposed_to_bed_print_pct(pct);
      }
      for (uii = 0; uii < 4; uii++) {
	axptr = alleles[uii];
	ucptr = &(writemap[4 * uii]);
	if (!axptr) {
	  memset(ucptr, 1, 4);
	} else if (axptr == salleles[0]) {
	  for (ujj = 0; ujj < 4; ujj++) {
	    axptr = alleles[ujj];
	    if (!axptr) {
	      *ucptr++ = 1;
	    } else if (axptr == salleles[0]) {
	      *ucptr++ = 3;
	    } else if (axptr == salleles[1]) {
	      *ucptr++ = 2;
	    } else {
	      *ucptr++ = 1;
	    }
	  }
	} else if (axptr == salleles[1]) {
	  for (ujj = 0; ujj < 4; ujj++) {
	    axptr = alleles[ujj];
	    if (!axptr) {
	      *ucptr++ = 1;
	    } else if (axptr == salleles[0]) {
	      *ucptr++ = 2;
	    } else if (axptr == salleles[1]) {
	      *ucptr++ = 0;
	    } else {
	      *ucptr++ = 1;
	    }
	  }
	} else {
	  memset(ucptr, 1, 4);
	}
      }
      uii = sample_ct & (~3U);
      ucptr = writebuf;
      for (ujj = 0; ujj < uii; ujj += 4) {
	*ucptr++ = writemap[prewritebuf[ujj]] | (writemap[prewritebuf[ujj + 1]] << 2) | (writemap[prewritebuf[ujj + 2]] << 4) | (writemap[prewritebuf[ujj + 3]] << 6);
      }
      ucc = 0;
      ucptr2 = &(prewritebuf[uii]);
      uii = sample_ct % 4;
      if (uii) {
	for (ujj = 0; ujj < uii; ujj++) {
	  ucc |= (writemap[*ucptr2++]) << (ujj * 2);
	}
	*ucptr = ucc;
      }
      fwrite(writebuf, 1, sample_ct4, outfile);
      if (!salleles[1]) {
	putc_unlocked(missing_geno, bimfile);
      } else {
	uii = strlen(salleles[1]);
	if (uii >= max_marker_allele_blen) {
	  max_marker_allele_blen = uii + 1;
	}
	fputs(salleles[1], bimfile);
      }
      putc_unlocked('\t', bimfile);
      if (!salleles[0]) {
	putc_unlocked(missing_geno, bimfile);
      } else {
	uii = strlen(salleles[0]);
	if (uii >= max_marker_allele_blen) {
	  max_marker_allele_blen = uii + 1;
	}
	fputs(salleles[0], bimfile);
      }
      for (uii = 0; uii < allele_tot; uii++) {
	if (alleles[uii][1]) {
	  free(alleles[uii]);
	}
      }
      allele_tot = 0;
      if (putc_checked('\n', bimfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      if (no_extra_cols) {
	cptr2 = skip_initial_spaces(cptr2);
	while (cptr2 == &(g_textbuf[MAXLINELEN - 1])) {
	  if (cptr2[-1] == '\n') {
	    break;
	  }
	  cptr2 = g_textbuf;
	  if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	    if (ferror(infile)) {
	      goto transposed_to_bed_ret_READ_FAIL;
	    }
	    g_textbuf[0] = '\0';
	    break;
	  }
	  cptr2 = skip_initial_spaces(cptr2);
	}
	if (!is_space_or_eoln(*cptr2)) {
	  no_extra_cols = 0;
	  putc_unlocked('\r', stdout);
	  logerrprint("Warning: Extra columns in .tped file.  Ignoring.\n");
	  transposed_to_bed_print_pct(pct);
	  goto transposed_to_bed_nextline;
	}
      } else {
      transposed_to_bed_nextline:
	cptr2 = (char*)memchr(cptr2, 0, MAXLINELEN - ((uintptr_t)(cptr2 - g_textbuf)));
	while (cptr2 == &(g_textbuf[MAXLINELEN - 1])) {
	  if (cptr2[-1] == '\n') {
	    break;
	  }
	  if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	    if (ferror(infile)) {
	      goto transposed_to_bed_ret_READ_FAIL;
	    }
	    break;
	  }
	  cptr2 = (char*)memchr(g_textbuf, 0, MAXLINELEN);
	}
      }
    }
    bigstack_end_reset(bigstack_end_mark);
    if (fclose_null(&infile)) {
      goto transposed_to_bed_ret_READ_FAIL;
    }
    if (fclose_null(&bimfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    if ((!marker_ct) && (!allow_no_variants)) {
      fputs("\b\b\b\b\b     \r", stdout);
      logerrprint("Error: Empty .tped file.\n");
      goto transposed_to_bed_ret_INVALID_FORMAT;
    }

    chrom_info_ptr->zero_extra_chroms = 0;
    if (map_is_unsorted) {
      loadbuf_size = 2 * max_marker_allele_blen + MAXLINELEN;
      bigstack_alloc(marker_ct * sizeof(int64_t)); // mapvals

      if (bigstack_alloc_ll(marker_ct, &ll_buf) ||
	  bigstack_alloc_ui(marker_ct, &pos_buf) ||
	  bigstack_alloc_c(marker_ct * max_marker_id_len, &marker_ids) ||
	  bigstack_alloc_d(marker_ct, &marker_cms) ||
	  bigstack_alloc_c(loadbuf_size, &loadbuf)) {
	goto transposed_to_bed_ret_NOMEM;
      }
      marker_allele_ptrs = (char**)bigstack_alloc(marker_ct * 2 * sizeof(intptr_t));
      if (!marker_allele_ptrs) {
	goto transposed_to_bed_ret_NOMEM;
      }
      // prevent cleanup from failing
      uint32_t allele_idx_end = marker_ct * 2;
      for (uint32_t allele_idx = 0; allele_idx < allele_idx_end; ++allele_idx) {
	marker_allele_ptrs[allele_idx] = (char*)missing_geno_ptr;
      }

      for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	pos_buf[marker_idx] = (uint32_t)((uint64_t)mapvals[marker_idx]);
	ll_buf[marker_idx] = (mapvals[marker_idx] & 0xffffffff00000000LLU) | marker_idx;
      }
      sort_marker_chrom_pos(ll_buf, marker_ct, pos_buf, chrom_start, chrom_id, nullptr, &chrom_ct);

      memcpy(outname_end, ".bim.tmp", 9);
      if (fopen_checked(outname, "r", &infile)) {
	goto transposed_to_bed_ret_OPEN_FAIL;
      }
      outname_end[4] = '\0';
      if (fopen_checked(outname, "w", &outfile)) {
	goto transposed_to_bed_ret_OPEN_FAIL;
      }
      marker_idx = 0;
      line_idx = 0;
      while (fgets(loadbuf, loadbuf_size, infile)) {
	line_idx++;
	// .tmp file, guaranteed to be no spaces in front
	cptr = skip_initial_spaces(token_endnn(loadbuf));
	cptr2 = token_endnn(cptr);
	cptr3 = skip_initial_spaces(cptr2);
	cptr4 = next_token_mult(cptr3, 2);
	uii = cptr2 - cptr;
	memcpyx(&(marker_ids[marker_idx * max_marker_id_len]), cptr, uii, '\0');
	if (scan_double(cptr3, &(marker_cms[marker_idx]))) {
	  sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of .tped file\n", line_idx);
	  goto transposed_to_bed_ret_INVALID_FORMAT_2R;
	}
	uii = strlen_se(cptr4);
	if (allele_set(cptr4, uii, &(marker_allele_ptrs[2 * marker_idx]))) {
	  goto transposed_to_bed_ret_NOMEM;
	}
	cptr4 = skip_initial_spaces(&(cptr4[uii + 1]));
	uii = strlen_se(cptr4);
	if (allele_set(cptr4, uii, &(marker_allele_ptrs[2 * marker_idx + 1]))) {
	  goto transposed_to_bed_ret_NOMEM;
	}
	marker_idx++;
      }
      if (!feof(infile)) {
	goto transposed_to_bed_ret_READ_FAIL;
      }
      fclose_null(&infile);
      marker_idx = 0;
      for (uii = 0; uii < chrom_ct; uii++) {
	cur_chrom = chrom_id[uii];
	ujj = chrom_start[uii + 1];
	cptr2 = chrom_name_write(chrom_info_ptr, cur_chrom, &(g_textbuf[MAXLINELEN]));
	*cptr2++ = '\t';
	for (; marker_idx < ujj; marker_idx++) {
	  marker_uidx = (uint32_t)ll_buf[marker_idx];
	  fwrite(&(g_textbuf[MAXLINELEN]), 1, cptr2 - (&(g_textbuf[MAXLINELEN])), outfile);
	  fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
	  g_textbuf[0] = '\t';
	  cptr = dtoa_gx(marker_cms[marker_uidx], '\t', &(g_textbuf[1]));
	  cptr = uint32toa_x((uint32_t)(ll_buf[marker_idx] >> 32), '\t', cptr);
	  if (fwrite_checked(g_textbuf, (uintptr_t)(cptr - g_textbuf), outfile)) {
	    goto transposed_to_bed_ret_WRITE_FAIL;
	  }
	  fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
	  putc_unlocked('\t', outfile);
	  fputs(marker_allele_ptrs[2 * marker_uidx + 1], outfile);
	  if (putc_checked('\n', outfile)) {
	    goto transposed_to_bed_ret_WRITE_FAIL;
	  }
	}
      }
      if (fclose_null(&outfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }

      outname_end[4] = '.';
      unlink(outname);

      outname_end[2] = 'e';
      outname_end[3] = 'd';
      if (fopen_checked(outname, FOPEN_RB, &infile)) {
	goto transposed_to_bed_ret_OPEN_FAIL;
      }
      outname_end[4] = '\0';
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto transposed_to_bed_ret_OPEN_FAIL;
      }
      if (fwrite_checked("l\x1b\x01", 3, outfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      uii = 0xfffffffeU; // last marker uidx
      for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	marker_uidx = (uint32_t)(ll_buf[marker_idx]);
	if (marker_uidx != uii + 1) {
	  if (fseeko(infile, 3 + ((uint64_t)marker_uidx) * sample_ct4, SEEK_SET)) {
	    goto transposed_to_bed_ret_READ_FAIL;
	  }
	}
	if (load_raw(sample_ct4, infile, (uintptr_t*)writebuf)) {
	  goto transposed_to_bed_ret_READ_FAIL;
	}
	if (fwrite_checked(writebuf, sample_ct4, outfile)) {
	  goto transposed_to_bed_ret_WRITE_FAIL;
	}
	uii = marker_uidx;
      }
      fclose_null(&infile);
      outname_end[4] = '.';
      unlink(outname);
      outname_end[4] = '\0';
    } else {
      uii = (outname_end - outname);
      memcpy(outname_end, ".bim.tmp", 9);
      memcpy(g_textbuf, outname, 9 + uii);
      outname_end[4] = '\0';
      if (rename(g_textbuf, outname)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      g_textbuf[uii + 2] = 'e';
      g_textbuf[uii + 3] = 'd';
      outname_end[2] = 'e';
      outname_end[3] = 'd';
      if (rename(g_textbuf, outname)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
    }
    fputs("\rProcessing .tped file... done.\n", stdout);
    *outname_end = '\0';
    LOGPRINTFWW("%s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);
  }

  while (0) {
  transposed_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  transposed_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  transposed_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  transposed_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  transposed_to_bed_ret_MISSING_TOKENS:
    putc_unlocked('\r', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .tped file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_HALF_MISSING:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .tped file has a half-missing call.\n", line_idx);
  transposed_to_bed_ret_INVALID_FORMAT_2R:
    putc_unlocked('\r', stdout);
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_TOO_MANY_ALLELES:
    putc_unlocked('\r', stdout);
    LOGERRPRINTF("Error: More than four alleles at variant %" PRIuPTR ".\n", marker_ct - 1);
    // retval already set
    break;
  }
 transposed_to_bed_ret_1:
  chrom_info_ptr->zero_extra_chroms = orig_zec;
  for (uii = 0; uii < allele_tot; uii++) {
    if (alleles[uii][1]) {
      free(alleles[uii]);
    }
  }
  cleanup_allele_storage(max_marker_allele_blen - 1, marker_ct * 2, marker_allele_ptrs);
  fclose_cond(infile);
  fclose_cond(bimfile);
  fclose_cond(outfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return retval;
}

int32_t vcf_sample_line(char* outname, char* outname_end, int32_t missing_pheno, char* bufptr, char* const_fid, uint32_t double_id, char id_delim, char vcf_idspace_to, char flag_char, uintptr_t* sample_ct_ptr) {
  FILE* outfile = nullptr;
  uintptr_t const_fid_len = 0;
  uintptr_t sample_ct = 0;
  int32_t retval = 0;
  char fam_trailer[20];
  uintptr_t fam_trailer_len;
  char* bufptr2;
  char* bufptr3;
  char* wptr;
  uintptr_t slen;
  bufptr2 = memcpya(fam_trailer, "\t0\t0\t0\t", 7);
  bufptr2 = int32toa_x(missing_pheno, '\n', bufptr2);
  fam_trailer_len = (uintptr_t)(bufptr2 - fam_trailer);

  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto vcf_sample_line_ret_OPEN_FAIL;
  }
  if (const_fid) {
    const_fid_len = strlen(const_fid);
  } else if ((!double_id) && (!id_delim)) {
    // default: --double-id + --id-delim
    double_id = 1;
    id_delim = '_';
  }
  if (id_delim != ' ') {
    bufptr2 = strchr(bufptr, ' ');
    if (bufptr2) {
      if (!vcf_idspace_to) {
	logerrprint("Error: VCF/BCF2 sample ID contains space(s).  Use --vcf-idspace-to to convert\nthem to another character, or \"--id-delim ' '\" to interpret the spaces as\nFID/IID delimiters.\n");
	goto vcf_sample_line_ret_INVALID_FORMAT;
      }
      do {
	*bufptr2 = vcf_idspace_to;
	bufptr2 = strchr(&(bufptr2[1]), ' ');
      } while (bufptr2);
    }
  }
  while (((unsigned char)bufptr[0]) >= ' ') {
    sample_ct++;
    bufptr2 = strchr(bufptr, '\t');
    if (bufptr2) {
      slen = (uintptr_t)(bufptr2 - bufptr);
    } else {
      slen = strlen_se(bufptr);
      bufptr2 = &(bufptr[slen]);
    }
    if (slen > MAX_ID_SLEN) {
      sprintf(g_logbuf, "Error: --%ccf does not support sample IDs longer than " MAX_ID_SLEN_STR " characters.\n", flag_char);
      goto vcf_sample_line_ret_INVALID_FORMAT_2;
    }
    if ((*bufptr == '0') && (slen == 1)) {
      logerrprint("Error: Sample ID cannot be '0'.\n");
      goto vcf_sample_line_ret_INVALID_FORMAT;
    }
    if (id_delim) {
      if (*bufptr == id_delim) {
	sprintf(g_logbuf, "Error: '%c' at beginning of sample ID.\n", id_delim);
	goto vcf_sample_line_ret_INVALID_FORMAT_2;
      } else if (bufptr[slen - 1] == id_delim) {
	sprintf(g_logbuf, "Error: '%c' at end of sample ID.\n", id_delim);
	goto vcf_sample_line_ret_INVALID_FORMAT_2;
      }
      bufptr3 = (char*)memchr(bufptr, (unsigned char)id_delim, slen);
      if (!bufptr3) {
	if (double_id) {
	  goto vcf_sample_line_double_id;
	} else if (const_fid) {
	  goto vcf_sample_line_const_id;
	} else {
	  sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
	  goto vcf_sample_line_ret_INVALID_FORMAT_2;
	}
      }
      if (memchr(&(bufptr3[1]), (unsigned char)id_delim, (uintptr_t)(bufptr2 - &(bufptr3[1])))) {
        LOGERRPRINTF("Error: Multiple instances of '%c' in sample ID.\n", id_delim);
	if (id_delim == '_') {
	  logerrprint("If you do not want '_' to be treated as a FID/IID delimiter, use --double-id or\n--const-fid to choose a different method of converting VCF sample IDs to PLINK\nIDs, or --id-delim to change the FID/IID delimiter.\n");
	}
        goto vcf_sample_line_ret_INVALID_FORMAT;
      }
      wptr = memcpyax(g_textbuf, bufptr, (uintptr_t)(bufptr3 - bufptr), '\t');
      bufptr3++;
      if ((*bufptr3 == '0') && (bufptr2 == &(bufptr3[1]))) {
        sprintf(g_logbuf, "Error: Sample ID ends with \"%c0\", which induces an invalid IID of '0'.\n", id_delim);
        goto vcf_sample_line_ret_INVALID_FORMAT_2;
      }
      wptr = memcpya(wptr, bufptr3, (uintptr_t)(bufptr2 - bufptr3));
    } else {
      if (double_id) {
      vcf_sample_line_double_id:
	wptr = memcpyax(g_textbuf, bufptr, (uintptr_t)(bufptr2 - bufptr), '\t');
      } else {
      vcf_sample_line_const_id:
        wptr = memcpyax(g_textbuf, const_fid, const_fid_len, '\t');
      }
      wptr = memcpya(wptr, bufptr, (uintptr_t)(bufptr2 - bufptr));
    }
    wptr = memcpya(wptr, fam_trailer, fam_trailer_len);
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto vcf_sample_line_ret_WRITE_FAIL;
    }
    if (*bufptr2 != '\t') {
      break;
    }
    bufptr = &(bufptr2[1]);
  }
  if (fclose_null(&outfile)) {
    goto vcf_sample_line_ret_WRITE_FAIL;
  }
  *sample_ct_ptr = sample_ct;
  while (0) {
  vcf_sample_line_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  vcf_sample_line_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  vcf_sample_line_ret_INVALID_FORMAT_2:
    logerrprintb();
  vcf_sample_line_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

uint32_t vcf_gp_invalid(char* bufptr, char* bufptr2, double vcf_min_gp, uint32_t gp_field_pos, uint32_t entry_idx, uint32_t* is_error_ptr) {
  double gp_val;
  uint32_t uii;
  for (uii = 0; uii < gp_field_pos; uii++) {
    bufptr = (char*)memchr(bufptr, ':', (uintptr_t)(bufptr2 - bufptr));
    if (!bufptr) {
      *is_error_ptr = 0;
      return 0;
    }
    bufptr++;
  }
  // do not treat decimal without leading zero as missing value
  if (((*bufptr == '.') || (*bufptr == '?')) && (bufptr[1] < '.')) {
    *is_error_ptr = 0;
    return 0;
  }
  for (uii = 0; uii < entry_idx; uii++) {
    bufptr = (char*)memchr(bufptr, ',', (uintptr_t)(bufptr2 - bufptr));
    if (!bufptr) {
      *is_error_ptr = 1;
      return 1;
    }
    bufptr++;
  }
  if (scan_double(bufptr, &gp_val)) {
    *is_error_ptr = 1;
    return 1;
  }
  *is_error_ptr = 0;
  return (gp_val < vcf_min_gp);
}

uint32_t vcf_gp_diploid_invalid(char* bufptr, char* bufptr2, double vcf_min_gp, uint32_t gp_field_pos, uint32_t idx1, uint32_t idx2, uint32_t* is_error_ptr) {
  uint32_t entry_idx;
  if (idx1 <= idx2) {
    entry_idx = ((idx2 * (idx2 + 1)) / 2) + idx1;
  } else {
    entry_idx = ((idx1 * (idx1 + 1)) / 2) + idx2;
  }
  return vcf_gp_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, entry_idx, is_error_ptr);
}

// oh, what the hell, may as well be liberal (even BCF2 does not support more
// than this)
#define MAX_VCF_ALT 65534

int32_t vcf_to_bed(char* vcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, char vcf_idspace_to, double vcf_min_qual, char* vcf_filter_exceptions_flattened, double vcf_min_gq, double vcf_min_gp, uint32_t vcf_half_call, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  FILE* bimfile = nullptr;
  FILE* skip3file = nullptr;
  char* sorted_fexcepts = nullptr;
  uintptr_t line_idx = 0;
  uintptr_t fexcept_ct = 0;
  uintptr_t max_fexcept_len = 5;
  uintptr_t sample_ct = 0;
  uintptr_t marker_skip_ct = 0;
  uintptr_t missing_gt_ct = 0;
  uint32_t double_id = (misc_flags / MISC_DOUBLE_ID) & 1;
  uint32_t check_qual = (vcf_min_qual != -1);
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t biallelic_only = (misc_flags / MISC_BIALLELIC_ONLY) & 1;
  uint32_t biallelic_strict = (misc_flags / MISC_BIALLELIC_ONLY_STRICT) & 1;
  uint32_t skip3_list = (misc_flags / MISC_BIALLELIC_ONLY_LIST) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t require_gt = (misc_flags / MISC_VCF_REQUIRE_GT) & 1;
  uint32_t marker_ct = 0;
  uint32_t gq_field_pos = 0;
  uint32_t gp_field_pos = 0;
  uint32_t vcf_half_call_explicit_error = (vcf_half_call == VCF_HALF_CALL_ERROR);
  int32_t retval = 0;
  char missing_geno = *g_missing_geno_ptr;
  uint32_t* vcf_alt_cts;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* ref_allele_ptr;
  char* delimiter_ptr;
  char* gq_scan_ptr;
  char* chrom_ptr;
  char* marker_id;
  char* pos_str;
  char* alt_alleles;
  char* geno_start;
  uintptr_t* base_bitfields;
  uintptr_t* alt_bitfield;
  uintptr_t* ref_ptr;
  uintptr_t* alt_ptr;
  uintptr_t final_mask;
  uintptr_t sample_ctl2;
  uintptr_t sample_ctv2;
  uintptr_t sample_ct4;
  uintptr_t sample_idx;
  uintptr_t loadbuf_size;
  uintptr_t slen;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t alt_allele_idx;
  double dxx;
  uint32_t marker_id_len;
  uint32_t alt_idx;
  uint32_t alt_ct;
  uint32_t ref_allele_len;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  char cc;
  {
    if (vcf_half_call_explicit_error) {
      vcf_half_call = 0;
    }
    // don't use gzopen_read_checked() since we want to customize the error
    // message
    gz_infile = gzopen(vcfname, FOPEN_RB);
    if (!gz_infile) {
      uii = strlen(vcfname);
      if ((uii > 4) && (!memcmp(&(vcfname[uii - 4]), ".vcf", 4))) {
	LOGERRPRINTFWW("Error: Failed to open %s.\n", vcfname);
      } else {
	LOGERRPRINTFWW("Error: Failed to open %s. (--vcf expects a complete filename; did you forget '.vcf' at the end?)\n", vcfname);
      }
      goto vcf_to_bed_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto vcf_to_bed_ret_NOMEM;
    }
    if (misc_flags & MISC_VCF_FILTER) {
      // automatically include "." and "PASS"
      fexcept_ct = 2;
      if (vcf_filter_exceptions_flattened) {
	fexcept_ct += count_and_measure_multistr(vcf_filter_exceptions_flattened, &max_fexcept_len);
      }
      if (bigstack_alloc_c(fexcept_ct * max_fexcept_len, &sorted_fexcepts)) {
	goto vcf_to_bed_ret_NOMEM;
      }
      memcpy(sorted_fexcepts, ".", 2);
      memcpy(&(sorted_fexcepts[max_fexcept_len]), "PASS", 5);
      if (vcf_filter_exceptions_flattened) {
	bufptr = vcf_filter_exceptions_flattened;
	for (ulii = 2; ulii < fexcept_ct; ulii++) {
	  slen = strlen(bufptr) + 1;
	  memcpy(&(sorted_fexcepts[ulii * max_fexcept_len]), bufptr, slen);
	  bufptr = &(bufptr[slen]);
	}
	qsort(sorted_fexcepts, fexcept_ct, max_fexcept_len, strcmp_casted);
	fexcept_ct = collapse_duplicate_ids(sorted_fexcepts, fexcept_ct, max_fexcept_len, nullptr);
	// there can't be many filter exceptions, so don't bother to free
	// unused memory in corner case
      }
    }

    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto vcf_to_bed_ret_NOMEM;
    }

    loadbuf = (char*)g_bigstack_base;
    loadbuf[loadbuf_size - 1] = ' ';
    while (1) {
      line_idx++;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	goto vcf_to_bed_ret_READ_FAIL;
      }
      if ((line_idx == 1) && (!memcmp(loadbuf, "BCF", 3))) {
	// this is more informative than "missing header line"...
	if (loadbuf[3] == 2) {
	  LOGPREPRINTFWW("Error: %s appears to be a BCF2 file. Try --bcf instead of --vcf.\n", vcfname);
	  goto vcf_to_bed_ret_INVALID_FORMAT_2;
	} else if (loadbuf[3] == 4) {
	  LOGPREPRINTFWW("Error: %s appears to be a BCF1 file. Use 'bcftools view' to convert it to a PLINK-readable VCF.\n", vcfname);
	  goto vcf_to_bed_ret_INVALID_FORMAT_2;
	}
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  goto vcf_to_bed_ret_LONG_LINE;
	}
	goto vcf_to_bed_ret_NOMEM;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      if (*bufptr != '#') {
	logerrprint("Error: Missing header line in .vcf file.\n");
	goto vcf_to_bed_ret_INVALID_FORMAT;
      }
      if (bufptr[1] != '#') {
	break;
      }
    }
    if (memcmp(bufptr, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", 38)) {
      logerrprint("Error: Improperly formatted .vcf header line.\n");
      goto vcf_to_bed_ret_INVALID_FORMAT;
    }
    bufptr = &(bufptr[38]);
    if (!memcmp(bufptr, "\tFORMAT\t", 8)) {
      retval = vcf_sample_line(outname, outname_end, missing_pheno, &(bufptr[8]), const_fid, double_id, id_delim, vcf_idspace_to, 'v', &sample_ct);
      if (retval) {
	goto vcf_to_bed_ret_1;
      }
    } else if (allow_no_samples) {
      memcpy(outname_end, ".fam", 5);
      if (fopen_checked(outname, "w", &outfile)) {
	goto vcf_to_bed_ret_OPEN_FAIL;
      }
      if (fclose_null(&outfile)) {
	goto vcf_to_bed_ret_WRITE_FAIL;
      }
    }
    if ((!sample_ct) && (!allow_no_samples)) {
      logerrprint("Error: No samples in .vcf file.\n");
      goto vcf_to_bed_ret_INVALID_FORMAT;
    }
    sample_ct4 = (sample_ct + 3) / 4;
    sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
    final_mask = (~ZEROLU) >> (2 * ((0x7fffffe0 - sample_ct) % BITCT2));
    if (bigstack_alloc_ul(sample_ctv2 * 10, &base_bitfields) ||
	bigstack_alloc_ui(MAX_VCF_ALT, &vcf_alt_cts)) {
      goto vcf_to_bed_ret_NOMEM;
    }
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(outname, "w", &bimfile)) {
      goto vcf_to_bed_ret_OPEN_FAIL;
    }
    memcpyl3(&(outname_end[2]), "ed");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto vcf_to_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto vcf_to_bed_ret_WRITE_FAIL;
    }
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto vcf_to_bed_ret_NOMEM;
    }

    loadbuf = (char*)g_bigstack_base;
    loadbuf[loadbuf_size - 1] = ' ';
    while (1) {
      line_idx++;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto vcf_to_bed_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  goto vcf_to_bed_ret_LONG_LINE;
	}
	goto vcf_to_bed_ret_NOMEM;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      // strchr instead of memchr since we explicitly need to catch premature
      // \0 here
      bufptr2 = strchr(bufptr, '\t');
      if (!bufptr2) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code_destructive(".vcf file", line_idx, allow_extra_chroms, bufptr, bufptr2, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto vcf_to_bed_ret_1;
      }
      if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	marker_skip_ct++;
	continue;
      }
      chrom_ptr = bufptr;
      pos_str = ++bufptr2;
      marker_id = strchr(bufptr2, '\t');
      if (!marker_id) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      if ((((unsigned char)(*pos_str)) - '0') >= 10) {
	sprintf(g_logbuf, "Error: Invalid variant bp coordinate on line %" PRIuPTR " of .vcf file.\n", line_idx);
	goto vcf_to_bed_ret_INVALID_FORMAT_2N;
      }
      ref_allele_ptr = strchr(++marker_id, '\t');
      if (!ref_allele_ptr) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      marker_id_len = (uintptr_t)(ref_allele_ptr - marker_id);
      bufptr = strchr(++ref_allele_ptr, '\t');
      // now ref_allele_ptr finally points to the ref allele
      if (!bufptr) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      ref_allele_len = (uintptr_t)(bufptr - ref_allele_ptr);
      alt_ct = 1;
      alt_alleles = ++bufptr;
      cc = *bufptr;
      // ',' < '.'
      while (1) {
	if ((unsigned char)cc <= ',' && (unsigned char)cc != '*') {
	  sprintf(g_logbuf, "Error: Invalid alternate allele on line %" PRIuPTR  " of .vcf file.\n", line_idx);
	  goto vcf_to_bed_ret_INVALID_FORMAT_2N;
	}
	bufptr2 = bufptr;
	do {
	  cc = *(++bufptr);
	  // allow GATK 3.4 <*:DEL> symbolic allele
	} while (((unsigned char)cc > ',') || (cc == '*'));
	if (((uintptr_t)(bufptr - bufptr2) == ref_allele_len) && (!memcmp(ref_allele_ptr, bufptr2, ref_allele_len))) {
	  if ((alt_ct != 1) || (cc == ',')) {
	    sprintf(g_logbuf, "Error: ALT allele duplicates REF allele on line %" PRIuPTR " of .vcf file.\n", line_idx);
	    goto vcf_to_bed_ret_INVALID_FORMAT_2N;
	  }
	  *alt_alleles = '.'; // tolerate SHAPEIT output
	}
	if (cc != ',') {
	  break;
	}
	cc = *(++bufptr);
	alt_ct++;
      }
      if (cc != '\t') {
	sprintf(g_logbuf, "Error: Malformed ALT field on line %" PRIuPTR " of .vcf file.\n", line_idx);
	goto vcf_to_bed_ret_INVALID_FORMAT_2N;
      }
      if (biallelic_strict && (alt_ct > 1)) {
	goto vcf_to_bed_skip3;
      }
      bufptr++;
      bufptr2 = strchr(bufptr, '\t');
      if (!bufptr2) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      if (check_qual) {
	if (*bufptr == '.') {
	  marker_skip_ct++;
	  continue;
	}
	if (scan_double(bufptr, &dxx)) {
	  sprintf(g_logbuf, "Error: Invalid QUAL value on line %" PRIuPTR " of .vcf file.\n", line_idx);
	  goto vcf_to_bed_ret_INVALID_FORMAT_2N;
	}
	if (dxx < vcf_min_qual) {
	  marker_skip_ct++;
	  continue;
	}
      }
      bufptr = &(bufptr2[1]);
      bufptr2 = strchr(bufptr, '\t');
      if (!bufptr2) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      bufptr2++;
      if (fexcept_ct) {
	// bugfix: recognize semicolon delimiter
	bufptr2[-1] = ';';
      vcf_to_bed_check_filter:
	delimiter_ptr = (char*)memchr(bufptr, ';', (uintptr_t)(bufptr2 - bufptr));
	if (bsearch_str(bufptr, (uintptr_t)(delimiter_ptr - bufptr), sorted_fexcepts, max_fexcept_len, fexcept_ct) == -1) {
	  marker_skip_ct++;
	  // if we replace the vcf_to_bed_check_filter goto with a while loop,
	  // can't use "continue" here
	  continue;
	}
	bufptr = &(delimiter_ptr[1]);
	if (bufptr != bufptr2) {
	  goto vcf_to_bed_check_filter;
	}
	bufptr2[-1] = '\t';
      }
      if (!sample_ct) {
	alt_allele_idx = 1;
	goto vcf_to_bed_skip_genotype_write;
      }
      bufptr = bufptr2;
      bufptr2 = strchr(bufptr, '\t');
      if (!bufptr2) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      bufptr = &(bufptr2[1]);
      bufptr2 = strchr(bufptr, '\t');
      if (!bufptr2) {
	goto vcf_to_bed_ret_MISSING_TOKENS;
      }
      if (memcmp(bufptr, "GT", 2)) {
	// We previously always skipped this case, but that's inconsistent with
	// how we now handle zero-sample VCFs.
	if (require_gt) {
	  marker_skip_ct++;
	  continue;
	}
	fill_quatervec_55(sample_ct, base_bitfields);
	missing_gt_ct++;
	alt_allele_idx = 1;
	goto vcf_to_bed_genotype_write;
      }
      bufptr2++;
      if (vcf_min_gq != -1) {
	gq_field_pos = 0;
	bufptr2[-1] = ':';
	gq_scan_ptr = bufptr;
	do {
	  gq_scan_ptr = (char*)memchr(gq_scan_ptr, ':', (uintptr_t)(bufptr2 - gq_scan_ptr));
	  if (++gq_scan_ptr == bufptr2) {
	    gq_field_pos = 0;
	    break;
	  }
	  gq_field_pos++;
	} while (memcmp(gq_scan_ptr, "GQ:", 3));
	bufptr2[-1] = '\t';
      }
      if (vcf_min_gp != -1) {
	gp_field_pos = 0;
	bufptr2[-1] = ':';
	do {
	  bufptr = (char*)memchr(bufptr, ':', (uintptr_t)(bufptr2 - bufptr));
	  if (++bufptr == bufptr2) {
	    gp_field_pos = 0;
	    break;
	  }
	  gp_field_pos++;
	} while (memcmp(bufptr, "GP:", 3));
	bufptr2[-1] = '\t';
      }
      bufptr = bufptr2;
      // okay, finally done with the line header
      if (alt_ct < 10) {
	// slightly faster parsing for the usual case
	fill_ulong_zero((alt_ct + 1) * sample_ctv2, base_bitfields);
	if ((!biallelic_only) || (alt_ct == 1)) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, bufptr = &(bufptr2[1])) {
	    bufptr2 = strchr(bufptr, '\t');
	    if (!bufptr2) {
	      if (sample_idx != sample_ct - 1) {
		goto vcf_to_bed_ret_MISSING_TOKENS;
	      }
	      bufptr2 = &(bufptr[strlen_se(bufptr)]);
	    }
	    uii = (unsigned char)(*bufptr) - '0';
	    // time to provide proper support for VCF import; that means, among
	    // other things, providing a useful error message instead of
	    // segfaulting on an invalid GT field, to help other tool
	    // developers.

	    if (uii <= 9) {
	      // no GQ field with ./. calls, so this check cannot occur earlier
	      if (gq_field_pos) {
		// to test: does splitting this off in an entirely separate
		// loop noticeably speed up common case parsing?  I hope
		// not--this is a predictable branch--but one can never be too
		// paranoid about this sort of performance leak when hundreds
		// of GB are involved...
		gq_scan_ptr = bufptr;
		for (ujj = 0; ujj < gq_field_pos; ujj++) {
		  gq_scan_ptr = (char*)memchr(gq_scan_ptr, ':', (uintptr_t)(bufptr2 - gq_scan_ptr));
		  if (!gq_scan_ptr) {
		    // non-GT fields are allowed to be missing
		    goto vcf_to_bed_missing_gq_1;
		  }
		  gq_scan_ptr++;
		}
		if ((!scan_double(gq_scan_ptr, &dxx)) && (dxx < vcf_min_gq)) {
		  continue;
		}
	      }

	    vcf_to_bed_missing_gq_1:
	      cc = bufptr[1];
	      if ((cc != '/') && (cc != '|')) {
		// haploid
	      vcf_to_bed_haploid_1:
		if (gp_field_pos) {
		  if (vcf_gp_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, &ukk)) {
		    if (ukk) {
		      goto vcf_to_bed_ret_INVALID_GP;
		    }
		    continue;
		  }
		}
		set_bit_ul(sample_idx * 2 + 1, &(base_bitfields[uii * sample_ctv2]));
	      } else {
		cc = bufptr[3];
		if (((cc != '/') && (cc != '|')) || (bufptr[4] == '.')) {
		  // code triploids, etc. as missing
		  // might want to subject handling of 0/0/. to --vcf-half-call
		  // control
		  ujj = ((unsigned char)bufptr[2]) - '0';
		  if (ujj > 9) {
		    if (ujj != (uint32_t)(((unsigned char)'.') - '0')) {
		      goto vcf_to_bed_ret_INVALID_GT;
		    }
		    if (!vcf_half_call) {
		      goto vcf_to_bed_ret_HALF_CALL_ERROR;
		    } else if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
		      goto vcf_to_bed_haploid_1;
		    } else if (vcf_half_call == VCF_HALF_CALL_REFERENCE) {
		      ujj = 0;
		      goto vcf_to_bed_reference_1;
		    }
		    // fall through on VCF_HALF_CALL_MISSING
		  } else {
		  vcf_to_bed_reference_1:
		    if (gp_field_pos) {
		      if (vcf_gp_diploid_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, ujj, &ukk)) {
			if (ukk) {
			  goto vcf_to_bed_ret_INVALID_GP;
			}
			continue;
		      }
		    }
		    set_bit_ul(sample_idx * 2, &(base_bitfields[uii * sample_ctv2]));
		    base_bitfields[ujj * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
		  }
		}
	      }
	    } else if (uii != (uint32_t)(((unsigned char)'.') - '0')) {
	      goto vcf_to_bed_ret_INVALID_GT;
	    } else if (vcf_half_call != VCF_HALF_CALL_MISSING) {
	      cc = bufptr[2];
	      if ((cc != '.') && ((bufptr[1] == '/') || (bufptr[1] == '|'))) {
		uii = ((unsigned char)cc) - '0';
		if (uii > 9) {
		  goto vcf_to_bed_ret_INVALID_GT;
		}
		if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
		  goto vcf_to_bed_haploid_1;
		} else if (!vcf_half_call) {
		  goto vcf_to_bed_ret_HALF_CALL_ERROR;
		} else {
		  // VCF_HALF_CALL_REFERENCE
		  ujj = 0;
		  goto vcf_to_bed_reference_1;
		}
	      }
	    }
	  }
	  alt_allele_idx = 1;
	  if (alt_ct > 1) {
	    ulii = popcount2_longs(&(base_bitfields[sample_ctv2]), sample_ctl2);
	    for (alt_idx = 2; alt_idx <= alt_ct; alt_idx++) {
	      uljj = popcount2_longs(&(base_bitfields[sample_ctv2 * alt_idx]), sample_ctl2);
	      if (uljj > ulii) {
		ulii = uljj;
		alt_allele_idx = alt_idx;
	      }
	    }
	  }
	} else {
	  // --biallelic-only, expect early termination in this case
	  alt_allele_idx = 0;
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, bufptr = &(bufptr2[1])) {
	    bufptr2 = strchr(bufptr, '\t');
	    if (!bufptr2) {
	      if (sample_idx != sample_ct - 1) {
		goto vcf_to_bed_ret_MISSING_TOKENS;
	      }
	      bufptr2 = &(bufptr[strlen_se(bufptr)]);
	    }
	    uii = (unsigned char)(*bufptr) - '0';
	    if (uii && (uii != alt_allele_idx)) {
	      if (uii == (uint32_t)(((unsigned char)'.') - '0')) {
		if (vcf_half_call == VCF_HALF_CALL_MISSING) {
		  continue;
		}
	        cc = bufptr[2];
		if ((cc == '.') || ((bufptr[1] != '/') && (bufptr[1] != '|'))) {
		  continue;
		}
		uii = ((unsigned char)cc) - '0';
		if (uii > 9) {
		  goto vcf_to_bed_ret_INVALID_GT;
		}
	        if (uii) {
		  if (!alt_allele_idx) {
		    alt_allele_idx = uii;
		  } else if (uii != alt_allele_idx) {
		    goto vcf_to_bed_skip3;
		  }
		}
		if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
		  goto vcf_to_bed_haploid_2;
		} else if (!vcf_half_call) {
		  goto vcf_to_bed_ret_HALF_CALL_ERROR;
		} else {
		  // VCF_HALF_CALL_REFERENCE
		  ujj = 0;
		  goto vcf_to_bed_reference_2;
		}
	      } else if (uii > 9) {
		goto vcf_to_bed_ret_INVALID_GT;
	      } else if (alt_allele_idx) {
		goto vcf_to_bed_skip3;
	      }
	      alt_allele_idx = uii;
	    }
	    if (gq_field_pos) {
	      gq_scan_ptr = bufptr;
	      for (ujj = 0; ujj < gq_field_pos; ujj++) {
		gq_scan_ptr = (char*)memchr(gq_scan_ptr, ':', (uintptr_t)(bufptr2 - gq_scan_ptr));
		if (!gq_scan_ptr) {
		  goto vcf_to_bed_missing_gq_2;
		}
		gq_scan_ptr++;
	      }
	      if ((!scan_double(gq_scan_ptr, &dxx)) && (dxx < vcf_min_gq)) {
		continue;
	      }
	    }
	  vcf_to_bed_missing_gq_2:
	    cc = bufptr[1];
	    if ((cc != '/') && (cc != '|')) {
	    vcf_to_bed_haploid_2:
	      if (gp_field_pos) {
		if (vcf_gp_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, &ukk)) {
		  if (ukk) {
		    goto vcf_to_bed_ret_INVALID_GP;
		  }
		  continue;
		}
	      }
	      set_bit_ul(sample_idx * 2 + 1, &(base_bitfields[uii * sample_ctv2]));
	    } else {
	      cc = bufptr[3];
	      if (((cc != '/') && (cc != '|')) || (bufptr[4] == '.')) {
		ujj = ((unsigned char)bufptr[2]) - '0';
		if (ujj && (ujj != alt_allele_idx)) {
		  if (ujj == (uint32_t)(((unsigned char)'.') - '0')) {
		    if (!vcf_half_call) {
		      goto vcf_to_bed_ret_HALF_CALL_ERROR;
		    } else if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
		      goto vcf_to_bed_haploid_2;
		    } else if (vcf_half_call == VCF_HALF_CALL_REFERENCE) {
		      ujj = 0;
		      goto vcf_to_bed_reference_2;
		    }
		    continue;
		  } else if (ujj > 9) {
		    goto vcf_to_bed_ret_INVALID_GT;
		  } else if (alt_allele_idx) {
		    goto vcf_to_bed_skip3;
		  }
		  alt_allele_idx = ujj;
		}
	      vcf_to_bed_reference_2:
		if (gp_field_pos) {
		  if (vcf_gp_diploid_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, ujj, &ukk)) {
		    if (ukk) {
		      goto vcf_to_bed_ret_INVALID_GP;
		    }
		    continue;
		  }
		}
		set_bit_ul(sample_idx * 2, &(base_bitfields[uii * sample_ctv2]));
		base_bitfields[ujj * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
	      }
	    }
	  }
	  if (!alt_allele_idx) {
	    alt_allele_idx = 1;
	  }
	}
	alt_bitfield = &(base_bitfields[alt_allele_idx * sample_ctv2]);
      } else {
	// bleah, multi-digit genotype codes
	// two-pass read: determine most common alt allele, then actually load
	// it
	fill_ulong_zero(2 * sample_ctv2, base_bitfields);
	alt_bitfield = &(base_bitfields[sample_ctv2]);
	fill_uint_zero(alt_ct, vcf_alt_cts);
	geno_start = bufptr;
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, bufptr = &(bufptr2[1])) {
	  bufptr2 = strchr(bufptr, '\t');
	  if (!bufptr2) {
	    if (sample_idx != sample_ct - 1) {
	      goto vcf_to_bed_ret_MISSING_TOKENS;
	    }
	    bufptr2 = &(bufptr[strlen_se(bufptr)]);
	  }
	  uii = (unsigned char)(*bufptr) - '0';
	  if (uii <= 9) {
	    if (gq_field_pos) {
	      gq_scan_ptr = bufptr;
	      for (ujj = 0; ujj < gq_field_pos; ujj++) {
		gq_scan_ptr = (char*)memchr(gq_scan_ptr, ':', (uintptr_t)(bufptr2 - gq_scan_ptr));
		if (!gq_scan_ptr) {
		  goto vcf_to_bed_missing_gq_3;
		}
		gq_scan_ptr++;
	      }
	      if ((!scan_double(gq_scan_ptr, &dxx)) && (dxx < vcf_min_gq)) {
		continue;
	      }
	    }
	  vcf_to_bed_missing_gq_3:
	    while (1) {
	      ujj = ((unsigned char)(*(++bufptr))) - 48;
	      if (ujj > 9) {
		break;
	      }
	      uii = uii * 10 + ujj;
	    }
	    // '/' = ascii 47, '|' = ascii 124
	    if ((ujj != 0xffffffffU) && (ujj != 76)) {
	      // haploid, count 2x
	    vcf_to_bed_haploid_3:
	      if (gp_field_pos) {
		if (vcf_gp_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, &ukk)) {
		  if (ukk) {
		    goto vcf_to_bed_ret_INVALID_GP;
		  }
		  continue;
		}
	      }
	      if (!uii) {
		set_bit_ul(sample_idx * 2 + 1, base_bitfields);
	      } else {
		vcf_alt_cts[uii - 1] += 2;
	      }
	    } else {
	      ujj = (unsigned char)(*(++bufptr)) - '0';
	      if (ujj > 9) {
		if (ujj == (uint32_t)(((unsigned char)'.') - '0')) {
		  if (!vcf_half_call) {
		    goto vcf_to_bed_ret_HALF_CALL_ERROR;
		  } else if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
		    goto vcf_to_bed_haploid_3;
		  } else if (vcf_half_call == VCF_HALF_CALL_REFERENCE) {
		    ujj = 0;
		    goto vcf_to_bed_reference_3;
		  } else {
		    continue;
		  }
		}
		goto vcf_to_bed_ret_INVALID_GT;
	      }
	      while (1) {
		ukk = ((unsigned char)(*(++bufptr))) - 48;
		if (ukk > 9) {
		  break;
		}
		ujj = ujj * 10 + ukk;
	      }
	      if (((ukk != 0xffffffffU) && (ukk != 76)) || (bufptr[1] == '.')) {
		// diploid; triploid+ skipped
	      vcf_to_bed_reference_3:
		if (gp_field_pos) {
		  if (vcf_gp_diploid_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, ujj, &ukk)) {
		    if (ukk) {
		      goto vcf_to_bed_ret_INVALID_GP;
		    }
		    continue;
		  }
		}
		if (!uii) {
		  set_bit_ul(sample_idx * 2, base_bitfields);
		} else {
		  vcf_alt_cts[uii - 1] += 1;
		}
		if (!ujj) {
		  base_bitfields[sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
		} else {
		  vcf_alt_cts[ujj - 1] += 1;
		}
	      }
	    }
	  } else if (uii != (uint32_t)(((unsigned char)'.') - '0')) {
	    goto vcf_to_bed_ret_INVALID_GT;
	  } else if ((vcf_half_call != VCF_HALF_CALL_MISSING) && (bufptr[2] != '.') && ((bufptr[1] == '/') || (bufptr[1] == '|'))) {
	    bufptr = &(bufptr[2]);
	    uii = ((unsigned char)(*bufptr)) - '0';
	    if (uii > 9) {
	      goto vcf_to_bed_ret_INVALID_GT;
	    }
	    while (1) {
	      ujj = ((unsigned char)(*(++bufptr))) - 48;
	      if (ujj > 9) {
		break;
	      }
	      uii = uii * 10 + ujj;
	    }
	    if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
	      goto vcf_to_bed_haploid_3;
	    } else if (!vcf_half_call) {
	      goto vcf_to_bed_ret_HALF_CALL_ERROR;
	    } else {
	      ujj = 0;
	      goto vcf_to_bed_reference_3;
	    }
	  }
	}
	alt_allele_idx = 0;
	uii = vcf_alt_cts[0];
	for (alt_idx = 1; alt_idx < alt_ct; alt_idx++) {
	  ujj = vcf_alt_cts[alt_idx];
	  if (biallelic_only && ujj && uii) {
	    goto vcf_to_bed_skip3;
	  }
	  if (ujj > uii) {
	    alt_allele_idx = alt_idx;
	    uii = vcf_alt_cts[alt_idx];
	  }
	}
	alt_allele_idx++;
	bufptr = geno_start;
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, bufptr = &(bufptr2[1])) {
	  bufptr2 = strchr(bufptr, '\t');
	  if (!bufptr2) {
	    bufptr2 = &(bufptr[strlen_se(bufptr)]);
	  }
	  if (*bufptr == '.') {
	    // validated on first pass
	    if ((vcf_half_call != VCF_HALF_CALL_MISSING) && (bufptr[2] != '.') && ((bufptr[1] == '/') || (bufptr[1] == '|'))) {
	      bufptr = &(bufptr[2]);
	      uii = ((unsigned char)(*bufptr)) - '0';
	      while (1) {
		ujj = ((unsigned char)(*(++bufptr))) - 48;
		if (ujj > 9) {
		  break;
		}
		uii = uii * 10 + ujj;
	      }
	      if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
		goto vcf_to_bed_haploid_4;
	      } else if (!vcf_half_call) {
		goto vcf_to_bed_ret_HALF_CALL_ERROR;
	      } else {
		ujj = 0;
		goto vcf_to_bed_reference_4;
	      }
	    }
	    continue;
	  }
	  if (gq_field_pos) {
	    gq_scan_ptr = bufptr;
	    for (ujj = 0; ujj < gq_field_pos; ujj++) {
	      gq_scan_ptr = (char*)memchr(gq_scan_ptr, ':', (uintptr_t)(bufptr2 - gq_scan_ptr));
	      gq_scan_ptr++;
	    }
	    if ((!scan_double(gq_scan_ptr, &dxx)) && (dxx < vcf_min_gq)) {
	      continue;
	    }
	  }
	  uii = (unsigned char)(*bufptr) - '0';
	  while (1) {
	    ujj = ((unsigned char)(*(++bufptr))) - 48;
	    if (ujj > 9) {
	      break;
	    }
	    uii = uii * 10 + ujj;
	  }
	  if ((ujj != 0xffffffffU) && (ujj != 76)) {
	    if (uii == alt_allele_idx) {
	    vcf_to_bed_haploid_4:
	      if (vcf_gp_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, &ukk)) {
		// no need for ukk check since already validated
		continue;
	      }
	      set_bit_ul(sample_idx * 2 + 1, alt_bitfield);
	    }
	  } else if (*(++bufptr) == '.') {
	    if (uii == alt_allele_idx) {
	      if (vcf_half_call == VCF_HALF_CALL_HAPLOID) {
	        goto vcf_to_bed_haploid_4;
	      } else if (!vcf_half_call) {
		goto vcf_to_bed_ret_HALF_CALL_ERROR;
	      } else if (vcf_half_call == VCF_HALF_CALL_REFERENCE) {
		ujj = 0;
		goto vcf_to_bed_reference_4;
	      }
	    }
	  } else {
	    ujj = (unsigned char)(*bufptr) - '0';
	    while (1) {
	      ukk = ((unsigned char)(*(++bufptr))) - 48;
	      if (ukk > 9) {
		break;
	      }
	      ujj = ujj * 10 + ukk;
	    }
	    if (((ukk != 0xffffffffU) && (ukk != 76)) || (bufptr[1] == '.')) {
	    vcf_to_bed_reference_4:
	      if (vcf_gp_diploid_invalid(bufptr, bufptr2, vcf_min_gp, gp_field_pos, uii, ujj, &ukk)) {
		continue;
	      }
	      if (uii == alt_allele_idx) {
		set_bit_ul(sample_idx * 2, alt_bitfield);
	      }
	      if (ujj == alt_allele_idx) {
		alt_bitfield[sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
	      }
	    }
	  }
	}
      }
      ref_ptr = base_bitfields;
      alt_ptr = alt_bitfield;
      for (sample_idx = 0; sample_idx < sample_ctl2; sample_idx++) {
	// take ref, then:
	// * if ref + alt is not two, force to 01
	// * otherwise, if ref is nonzero, add 1 to match PLINK binary encoding
	ulii = *ref_ptr;
	uljj = *alt_ptr++;
	ulkk = (ulii + uljj) & AAAAMASK;
	uljj = ulii + ((ulii | (ulii >> 1)) & FIVEMASK);
	ulii = ulkk | (ulkk >> 1); // 11 in nonmissing positions
	*ref_ptr++ = (uljj & ulii) | (((~ulkk) >> 1) & FIVEMASK);
      }
      ref_ptr[-1] &= final_mask;
    vcf_to_bed_genotype_write:
      if (fwrite_checked(base_bitfields, sample_ct4, outfile)) {
	goto vcf_to_bed_ret_WRITE_FAIL;
      }
    vcf_to_bed_skip_genotype_write:
      // chrom_ptr already null-terminated
      fputs(chrom_ptr, bimfile);
      putc_unlocked('\t', bimfile);
      fwrite(marker_id, 1, marker_id_len + 1, bimfile);
      putc_unlocked('0', bimfile);
      putc_unlocked('\t', bimfile);
      fwrite(pos_str, 1, marker_id - pos_str, bimfile);

      if (*alt_alleles == '.') {
	putc_unlocked(missing_geno, bimfile);
      } else {
	bufptr = alt_alleles;
	for (alt_idx = 1; alt_idx < alt_allele_idx; alt_idx++) {
	  bufptr = strchr(bufptr, ',');
	  bufptr++;
	}
	bufptr2 = strchr(bufptr, (alt_allele_idx == alt_ct)? '\t' : ',');
	*bufptr2 = '\0';
	fputs(bufptr, bimfile);
      }
      putc_unlocked('\t', bimfile);
      alt_alleles[-1] = '\n';
      *alt_alleles = '\0';
      if (((((unsigned char)ref_allele_ptr[0]) & 0xdf) == 'N') && (ref_allele_ptr[1] == '\t')) {
	*ref_allele_ptr = missing_geno;
      }
      if (fputs_checked(ref_allele_ptr, bimfile)) {
	goto vcf_to_bed_ret_WRITE_FAIL;
      }
      marker_ct++;
      if (!(marker_ct % 1000)) {
	printf("\r--vcf: %uk variants complete.", marker_ct / 1000);
	fflush(stdout);
      }
      continue;
    vcf_to_bed_skip3:
      if (skip3_list) {
	if (!marker_skip_ct) {
	  memcpy(outname_end, ".skip.3allele", 14);
	  if (fopen_checked(outname, "w", &skip3file)) {
	    goto vcf_to_bed_ret_OPEN_FAIL;
	  }
	  memcpy(outname_end, ".bed", 5);
	}
	marker_id[marker_id_len] = '\0';
	if (fputs_checked(marker_id, skip3file)) {
	  goto vcf_to_bed_ret_WRITE_FAIL;
	}
	putc_unlocked('\n', skip3file);
      }
      marker_skip_ct++;
    }
    if (fclose_null(&bimfile) || fclose_null(&outfile)) {
      goto vcf_to_bed_ret_WRITE_FAIL;
    }
    if (skip3file) {
      if (fclose_null(&skip3file)) {
	goto vcf_to_bed_ret_WRITE_FAIL;
      }
    }
    putc_unlocked('\r', stdout);
    if ((!marker_ct) && (!allow_no_variants)) {
      if (marker_skip_ct) {
	logerrprint("Error: All variants in VCF skipped.\n");
	retval = RET_ALL_MARKERS_EXCLUDED;
	goto vcf_to_bed_ret_1;
      } else {
	logerrprint("Error: No variants in VCF file.\n");
	goto vcf_to_bed_ret_INVALID_FORMAT;
      }
    }
    *outname_end = '\0';
    LOGPRINTFWW("--vcf: %s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);
    if (marker_skip_ct) {
      LOGPRINTF("(%" PRIuPTR " variant%s skipped.)\n", marker_skip_ct, (marker_skip_ct == 1)? "" : "s");
    }
    if (missing_gt_ct) {
      LOGERRPRINTF("Warning: %" PRIuPTR " variant record%s had no GT field.\n", missing_gt_ct, (missing_gt_ct == 1)? "" : "s");
    }
  }
  while (0) {
  vcf_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  vcf_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  vcf_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  vcf_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  vcf_to_bed_ret_HALF_CALL_ERROR:
    logprint("\n");
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .vcf file has a GT half-call.\n", line_idx);
    if (!vcf_half_call_explicit_error) {
      logerrprint("Use --vcf-half-call to specify how these should be processed.\n");
    }
    retval = RET_INVALID_FORMAT;
    break;
  vcf_to_bed_ret_INVALID_GP:
    logprint("\n");
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .vcf file has an improperly formatted GP field.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  vcf_to_bed_ret_INVALID_GT:
    logprint("\n");
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .vcf file has an invalid GT field.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  vcf_to_bed_ret_MISSING_TOKENS:
    logprint("\n");
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .vcf file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  vcf_to_bed_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of .vcf file is pathologically long.\n", line_idx);
  vcf_to_bed_ret_INVALID_FORMAT_2N:
    logprint("\n");
  vcf_to_bed_ret_INVALID_FORMAT_2:
    logerrprintb();
  vcf_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 vcf_to_bed_ret_1:
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  fclose_cond(bimfile);
  fclose_cond(skip3file);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t read_bcf_typed_nonnegative_integer(gzFile gz_infile, uint32_t* int_ptr) {
  // errors out on missing and negative values.
  int32_t retval = 0;
  int32_t ii = gzgetc(gz_infile);
  uint32_t uii;
  if (ii == -1) {
    goto read_bcf_typed_nonnegative_integer_ret_READ_OR_FORMAT_FAIL;
  }
  if (ii == 0x11) {
    ii = gzgetc(gz_infile);
    if (((uint32_t)ii) > 127) {
      if (ii == -1) {
	goto read_bcf_typed_nonnegative_integer_ret_READ_OR_FORMAT_FAIL;
      }
      goto read_bcf_typed_nonnegative_integer_ret_INVALID_FORMAT_GENERIC;
    }
    *int_ptr = (uint32_t)ii;
  } else if (ii == 0x12) {
    uii = gzgetc(gz_infile);
    ii = gzgetc(gz_infile);
    if (((uint32_t)ii) > 127) {
      if (ii == -1) {
	goto read_bcf_typed_nonnegative_integer_ret_READ_OR_FORMAT_FAIL;
      }
      goto read_bcf_typed_nonnegative_integer_ret_INVALID_FORMAT_GENERIC;
    }
    *int_ptr = uii | (((uint32_t)ii) << 8);
  } else if (ii == 0x13) {
    if (gzread(gz_infile, int_ptr, 4) < 4) {
      goto read_bcf_typed_nonnegative_integer_ret_READ_OR_FORMAT_FAIL;
    }
  } else {
    goto read_bcf_typed_nonnegative_integer_ret_INVALID_FORMAT_GENERIC;
  }
  while (0) {
  read_bcf_typed_nonnegative_integer_ret_READ_OR_FORMAT_FAIL:
    if (!gzeof(gz_infile)) {
      retval = RET_READ_FAIL;
      break;
    }
  read_bcf_typed_nonnegative_integer_ret_INVALID_FORMAT_GENERIC:
    logerrprint("Error: Improperly formatted .bcf file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int32_t read_bcf_typed_string(gzFile gz_infile, char* readbuf, uint32_t maxlen, uint32_t* len_ptr) {
  int32_t retval = 0;
  uint32_t slen;
  int32_t ii;
  ii = gzgetc(gz_infile);
  if (ii == -1) {
    goto read_bcf_typed_string_ret_READ_OR_FORMAT_FAIL;
  }
  if (((uint32_t)ii) == 0xf7) {
    retval = read_bcf_typed_nonnegative_integer(gz_infile, &slen);
    if (retval) {
      goto read_bcf_typed_string_ret_1;
    }
    if (slen > maxlen) {
      logerrprint("Error: Excessively long typed string in .bcf file.\n");
      goto read_bcf_typed_string_ret_INVALID_FORMAT;
    }
  } else if ((((uint32_t)ii) & 0x0f) == 0x07) {
    slen = ((uint32_t)ii) >> 4;
  } else if (ii) {
    goto read_bcf_typed_string_ret_INVALID_FORMAT_GENERIC;
  } else {
    slen = 0;
  }
  *len_ptr = slen;
  if (slen) {
    if ((uint32_t)((uint64_t)gzread(gz_infile, readbuf, slen)) != slen) {
      goto read_bcf_typed_string_ret_READ_OR_FORMAT_FAIL;
    }
  }
  while (0) {
  read_bcf_typed_string_ret_READ_OR_FORMAT_FAIL:
    if (!gzeof(gz_infile)) {
      retval = RET_READ_FAIL;
      break;
    }
  read_bcf_typed_string_ret_INVALID_FORMAT_GENERIC:
    logerrprint("Error: Improperly formatted .bcf file.\n");
  read_bcf_typed_string_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 read_bcf_typed_string_ret_1:
  return retval;
}

int32_t bcf_header_line_idx_check(char* line_iter, char* line_end, uint32_t* gt_idx_ptr) {
  // adapted from plink 2.0 BcfHeaderLineIdxCheck()
  while (1) {
    char* tag_start = &(line_iter[1]);
    if (!memcmp(tag_start, "IDX=", 4)) {
      char* val_start = &(tag_start[4]);
      if (memchr(val_start, ',', line_end - val_start)) {
        // force this to have the same limitation as plink 2.0; if it's ever
        // worth lifting this restriction, we do so in both programs, not in
        // 1.9 only
        logerrprint("Error: FORMAT:GT line in BCF text header block has IDX= in the center instead\nof the end of the line; this is not currently supported by " PROG_NAME_STR ".  Contact us\nif you need this to work.\n");
        return 1;
      }
      if (scan_posint_defcap(val_start, gt_idx_ptr)) {
        logerrprint("Error: Invalid FORMAT:GT IDX= value in BCF text header block.\n");
        return 1;
      }
      return 0;
    }
    line_iter = strchr(tag_start, '=');
    ++line_iter;
    if (*line_iter != '"') {
      line_iter = (char*)memchr(line_iter, ',', line_end - line_iter);
      if (line_iter) {
        continue;
      }
      return 0;
    }
    // Need to worry about backslash-escaped characters.
    ++line_iter;
    while (1) {
      char cc = *line_iter;
      if (cc == '"') {
        break;
      }
      if ((cc != '\\') && (cc != '\n')) {
        ++line_iter;
        continue;
      }
      if ((cc == '\n') || (line_iter[1] == '\n')) {
        goto bcf_header_line_idx_check_FAIL;
      }
      line_iter = &(line_iter[2]);
    }
    ++line_iter;
    char cc = *line_iter;
    if (cc == ',') {
      continue;
    }
    // bugfix (19 Feb 2020)
    if (cc == '>') {
      return 0;
    }
    break;
  }
 bcf_header_line_idx_check_FAIL:
  logerrprint("Error: FORMAT:GT line in BCF text header block is malformed.\n");
  return 1;
}

int32_t bcf_to_bed(char* bcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, char vcf_idspace_to, double vcf_min_qual, char* vcf_filter_exceptions_flattened, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  FILE* bimfile = nullptr;
  FILE* skip3file = nullptr;
  char* sorted_fexcepts = nullptr;
  uintptr_t* fexcept_bitfield = nullptr;
  uint32_t* fexcept_idxs = nullptr;
  Ll_str* contig_list = nullptr;
  char* tbuf2 = &(g_textbuf[MAXLINELEN]);
  uintptr_t contig_ct = 0;
  uintptr_t max_contig_len = 0;
  uintptr_t max_fexcept_len = 0;
  uintptr_t fexcept_ct = 0;
  uintptr_t marker_skip_ct = 0;
  uintptr_t missing_gt_ct = 0;
  uint32_t double_id = (misc_flags / MISC_DOUBLE_ID) & 1;
  uint32_t check_qual = (vcf_min_qual != -1);
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t biallelic_only = (misc_flags / MISC_BIALLELIC_ONLY) & 1;
  uint32_t biallelic_strict = (misc_flags / MISC_BIALLELIC_ONLY_STRICT) & 1;
  uint32_t skip3_list = (misc_flags / MISC_BIALLELIC_ONLY_LIST) & 1;
  uint32_t vcf_filter = (misc_flags / MISC_VCF_FILTER) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t require_gt = (misc_flags / MISC_VCF_REQUIRE_GT) & 1;
  uint32_t sample_ct = 0;
  uint32_t stringdict_ct = 1;
  uint32_t gt_idx = 0;
  uint32_t marker_ct = 0;
  uint32_t umm = 0;
  uint32_t vcf_min_qualf_compare_bits = 0;
  int32_t retval = 0;
  char missing_geno = *g_missing_geno_ptr;
  uint32_t bcf_var_header[8];
  Ll_str* ll_ptr;
  uintptr_t* contig_bitfield;
  uintptr_t* base_bitfields;
  uintptr_t* ref_ptr;
  uintptr_t* alt_ptr;
  char* loadbuf;
  char* loadbuf_end;
  char* linebuf;
  char* linebuf_end;
  char* contigdict;
  char* bufptr;
  char* bufptr2;
  char* marker_id;
  char* allele_buf;
  char** allele_ptrs;
  uint32_t* allele_lens;
  uint32_t* vcf_alt_cts;
  uint32_t* uiptr;
  uint16_t* ui16ptr;
  unsigned char* ucptr;
  uintptr_t final_mask;
  uintptr_t max_allele_ct;
  uintptr_t slen;
  uintptr_t sample_ctv2;
  uintptr_t alt_allele_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint64_t lastloc;
  uint64_t ullii;
  uint64_t ulljj;
  float vcf_min_qualf;
  uint32_t sample_ct4;
  uint32_t sample_ctl2;
  uint32_t header_size;
  uint32_t marker_id_len;
  uint32_t n_allele;
  uint32_t sample_idx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t ii;
  {
    if (check_qual) {
      if (vcf_min_qual > FLT_MAXD) {
	logerrprint("Error: --vcf-min-qual parameter too large.\n");
	goto bcf_to_bed_ret_INVALID_CMDLINE;
      }
      vcf_min_qualf = (float)vcf_min_qual;
      memcpy(&vcf_min_qualf_compare_bits, &vcf_min_qualf, 4);
      // +infinity = 0x7f800000; this should pass the comparison
      // quiet nan = 0x7f800001; this (and other nans) should fail
      vcf_min_qualf_compare_bits += 0x807fffffU;
    }
    // todo: check if a specialized bgzf reader can do faster forward seeks
    // when we don't have precomputed virtual offsets
    retval = gzopen_read_checked(bcfname, &gz_infile);
    if (retval) {
      goto bcf_to_bed_ret_1;
    }
    if (gzread(gz_infile, g_textbuf, 5) < 5) {
      goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
    }
    if (memcmp(g_textbuf, "BCF\2", 4)) {
      if (memcmp(g_textbuf, "BCF\4", 4)) {
	LOGPREPRINTFWW("Error: %s is not a BCF2 file.\n", bcfname);
      } else {
	LOGPREPRINTFWW("Error: %s appears to be a BCF1 file; --bcf only supports BCF2. Use 'bcftools view' to convert it to a PLINK-readable VCF.\n", bcfname);
      }
      goto bcf_to_bed_ret_INVALID_FORMAT_2;
    }
    if (((unsigned char)(g_textbuf[4])) > 2) {
      // defend against 0x82-0x87 being given a meaning in 8-bit int vectors,
      // etc.
      LOGPREPRINTFWW("Error: %s appears to be formatted as BCFv2.%u; this PLINK build only supports v2.0-2.2. You may need to obtain an updated version of PLINK.\n", bcfname, ((unsigned char)(g_textbuf[4])));
      goto bcf_to_bed_ret_INVALID_FORMAT_2;
    }
    if (gzread(gz_infile, &header_size, 4) < 4) {
      goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
    }
    // must have at least fileformat, and first eight fields of #CHROM line.
    // GT not required with --allow-no-samples, contig not require with
    // --allow-no-vars.
    if (header_size < 59) {
      goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
    }
    if (vcf_filter_exceptions_flattened) {
      // vcf_filter guaranteed to be true
      fexcept_ct = count_and_measure_multistr(vcf_filter_exceptions_flattened, &max_fexcept_len);
      if (bigstack_end_alloc_c(fexcept_ct * max_fexcept_len, &sorted_fexcepts)) {
	goto bcf_to_bed_ret_NOMEM;
      }
      bufptr = vcf_filter_exceptions_flattened;
      for (ulii = 0; ulii < fexcept_ct; ulii++) {
	slen = strlen(bufptr) + 1;
	memcpy(&(sorted_fexcepts[ulii * max_fexcept_len]), bufptr, slen);
	bufptr = &(bufptr[slen]);
      }
      qsort(sorted_fexcepts, fexcept_ct, max_fexcept_len, strcmp_casted);
      fexcept_ct = collapse_duplicate_ids(sorted_fexcepts, fexcept_ct, max_fexcept_len, nullptr);
      if (bigstack_end_calloc_ui(fexcept_ct, &fexcept_idxs)) {
	goto bcf_to_bed_ret_NOMEM;
      }
    }
    if (bigstack_left() <= header_size) {
      goto bcf_to_bed_ret_NOMEM;
    }
    loadbuf = (char*)bigstack_alloc(header_size + 1);
    if ((uint32_t)((uint64_t)gzread(gz_infile, loadbuf, header_size)) != header_size) {
      goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
    }
    if (!(*loadbuf)) {
      goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
    }
    if (!loadbuf[header_size - 1]) {
      loadbuf_end = &(loadbuf[header_size - 2]);
      while (!(*loadbuf_end)) {
	loadbuf_end--;
      }
      loadbuf_end++;
      header_size = (uintptr_t)(loadbuf_end - loadbuf);
    } else {
      loadbuf_end = &(loadbuf[header_size]);
    }
    *loadbuf_end = '\n';
    header_size++;
    linebuf = loadbuf;
    while (1) {
      linebuf_end = (char*)memchr(linebuf, '\n', header_size);
      if (linebuf[0] != '#') {
	goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
      }
      if (linebuf[1] != '#') {
	if (linebuf[1] == 'C') {
	  break; // end of meta-info
	}
	goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
      }
      if (linebuf[2] == 'F') {
	if (!memcmp(&(linebuf[3]), "ORMAT=<ID=", 10)) {
	  if (!memcmp(&(linebuf[13]), "GT,", 3)) {
	    if (gt_idx) {
	      logerrprint("Error: Duplicate GT format specifier in .bcf file.\n");
	      goto bcf_to_bed_ret_INVALID_FORMAT;
	    }
	    if (memcmp(&(linebuf[16]), "Number=1,Type=String,Description=", 33)) {
	      logerrprint("Error: Unrecognized GT field format in .bcf file.\n");
	      goto bcf_to_bed_ret_INVALID_FORMAT;
	    }
	    gt_idx = stringdict_ct;
            // bugfix (17 Feb 2020): BCFv2.2 headers may explicitly specify the
            // index.  Usually it'll be equal to stringdict_ct anyway, but it
            // can naturally diverge when e.g. a FILTER and INFO key are
            // identical.  Search the header line for an "IDX=<#>" entry.
            // &(linebuf[36]) points to the comma before the beginning of the
            // "Description" tag.
            if (bcf_header_line_idx_check(&(linebuf[36]), linebuf_end, &gt_idx)) {
              goto bcf_to_bed_ret_INVALID_FORMAT;
            }
	  }
	  stringdict_ct++;
	} else if (!memcmp(&(linebuf[3]), "ILTER=<ID=", 10)) {
	  bufptr = &(linebuf[13]);
	  bufptr2 = (char*)memchr(bufptr, ',', linebuf_end - bufptr);
	  if (bufptr2 == linebuf_end) {
	    goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
	  }
	  if (memcmp(bufptr, "PASS,", 5)) {
	    if (fexcept_ct) {
	      ii = bsearch_str(bufptr, (uintptr_t)(bufptr2 - bufptr), sorted_fexcepts, max_fexcept_len, fexcept_ct);
	      if (ii != -1) {
		if (fexcept_idxs[(uint32_t)ii]) {
		  goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
		}
		fexcept_idxs[(uint32_t)ii] = stringdict_ct;
	      }
	    }
	    stringdict_ct++;
	  }
	}
      } else if (!memcmp(&(linebuf[2]), "INFO=<ID=", 9)) {
	stringdict_ct++;
      } else if (!memcmp(&(linebuf[2]), "contig=<ID=", 11)) {
	bufptr = &(linebuf[13]);
	bufptr2 = (char*)memchr(bufptr, ',', linebuf_end - bufptr);
	if (bufptr2 == linebuf_end) {
	  goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
	}
	slen = (uintptr_t)(bufptr2 - bufptr);
	if (slen >= max_contig_len) {
	  max_contig_len = slen + 1;
	}
	if (bigstack_end_alloc_llstr(slen + 1, &ll_ptr)) {
	  goto bcf_to_bed_ret_NOMEM;
	}
	ll_ptr->next = contig_list;
	memcpyx(ll_ptr->ss, bufptr, slen, '\0');
	contig_list = ll_ptr;
	contig_ct++;
      }
      linebuf = &(linebuf_end[1]);
      if (linebuf >= loadbuf_end) {
	goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
      }
    }
    if ((!allow_no_variants) && (!contig_ct)) {
      logerrprint("Error: No contig fields in .bcf header.\n");
      goto bcf_to_bed_ret_INVALID_FORMAT;
    }
    if (memcmp(linebuf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", 38)) {
      goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
    }
    if (!memcmp(&(linebuf[38]), "\tFORMAT\t", 8)) {
      *linebuf_end = '\0';
      retval = vcf_sample_line(outname, outname_end, missing_pheno, &(linebuf[46]), const_fid, double_id, id_delim, vcf_idspace_to, 'b', &ulii);
      if (retval) {
	goto bcf_to_bed_ret_1;
      }
      if (ulii >= 0x1000000) {
	// variant records only have 24 bits allocated for n_sample
	logerrprint("Error: .bcf file contains >= 2^24 sample IDs.\n");
	goto bcf_to_bed_ret_INVALID_FORMAT;
      }
      sample_ct = ulii;
    } else if (allow_no_samples) {
      gt_idx = 0;
      memcpy(outname_end, ".fam", 5);
      if (fopen_checked(outname, "w", &outfile)) {
	goto bcf_to_bed_ret_OPEN_FAIL;
      }
      if (fclose_null(&outfile)) {
	goto bcf_to_bed_ret_WRITE_FAIL;
      }
    }
    if ((!sample_ct) && (!allow_no_samples)) {
      logerrprint("Error: No samples in .bcf file.\n");
      goto bcf_to_bed_ret_INVALID_FORMAT;
    }
    sample_ct4 = (sample_ct + 3) / 4;
    sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
    bigstack_reset(loadbuf);
    ulii = BITCT_TO_WORDCT(contig_ct);
    if (bigstack_calloc_ul(ulii, &contig_bitfield) ||
	bigstack_alloc_c(contig_ct * max_contig_len, &contigdict)) {
      goto bcf_to_bed_ret_NOMEM;
    }
    ulii = contig_ct;
    while (ulii) {
      ulii--;
      const uint32_t chrom_name_slen = strlen(contig_list->ss);
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code(contig_list->ss, ".bcf file", 0, chrom_name_slen, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto bcf_to_bed_ret_1;
      }
      if (is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	set_bit_ul(ulii, contig_bitfield);
	strcpy(&(contigdict[ulii * max_contig_len]), contig_list->ss);
      }
      contig_list = contig_list->next;
    }
    if (vcf_filter) {
      uii = BITCT_TO_WORDCT(stringdict_ct);
      if (bigstack_calloc_ul(uii, &fexcept_bitfield)) {
	goto bcf_to_bed_ret_NOMEM;
      }
      fexcept_bitfield[0] = 1; // 'PASS'
      for (ulii = 0; ulii < fexcept_ct; ulii++) {
	// fexcept_idxs[] not dereferenced if --vcf-filter had no parameters
	SET_BIT(fexcept_idxs[ulii], fexcept_bitfield);
      }
    }
    bigstack_end_reset(bigstack_end_mark);

    final_mask = (~ZEROLU) >> (2 * ((0x7fffffe0 - sample_ct) % BITCT2));
    if (bigstack_alloc_c(sample_ct * 12, &loadbuf) ||
	bigstack_alloc_c(65536, &marker_id) ||
	bigstack_alloc_c(NON_BIGSTACK_MIN, &allele_buf) ||
	bigstack_alloc_ui(65535, &allele_lens) ||
	bigstack_alloc_ui(MAX_VCF_ALT, &vcf_alt_cts)) {
      goto bcf_to_bed_ret_NOMEM;
    }
    allele_ptrs = (char**)bigstack_alloc(65535 * sizeof(intptr_t));
    if (!allele_ptrs) {
      goto bcf_to_bed_ret_NOMEM;
    }
    max_allele_ct = bigstack_left() / (sample_ctv2 * sizeof(intptr_t));
    if (max_allele_ct < 3) {
      goto bcf_to_bed_ret_NOMEM;
    } else if (max_allele_ct > 65535) {
      max_allele_ct = 65535;
    }
    bigstack_alloc_ul(sample_ctv2 * max_allele_ct, &base_bitfields);
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(outname, "w", &bimfile)) {
      goto bcf_to_bed_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto bcf_to_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto bcf_to_bed_ret_WRITE_FAIL;
    }
    if ((!gt_idx) && require_gt) {
      if (!allow_no_variants) {
	logerrprint("Error: .bcf header doesn't define FORMAT:GT.\n");
	retval = RET_ALL_MARKERS_EXCLUDED;
	goto bcf_to_bed_ret_1;
      }
      logerrprint("Warning: Skipping all variants since .bcf header doesn't define FORMAT:GT.\n");
      goto bcf_to_bed_skip_all_variants;
    }
    // possible todo: optimize other no-GT cases.  e.g. if no sample
    // information is needed, don't write the .bed or .fam.

    memcpyl3(tbuf2, "\t0\t");
    while (1) {
      lastloc = gztell(gz_infile) + 8;
      if (gzread(gz_infile, bcf_var_header, 32) < 32) {
	break;
      }
      if ((bcf_var_header[0] <= 24) || (bcf_var_header[2] >= contig_ct)) {
	goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
      }
      if ((!bcf_var_header[1]) || (!is_set(contig_bitfield, bcf_var_header[2]))) {
	goto bcf_to_bed_marker_skip;
      }
      if (check_qual) {
	if (bcf_var_header[5] + 0x807fffffU < vcf_min_qualf_compare_bits) {
	  goto bcf_to_bed_marker_skip;
	}
      }
      retval = read_bcf_typed_string(gz_infile, marker_id, 65535, &marker_id_len);
      if (retval) {
	goto bcf_to_bed_ret_1;
      }
      n_allele = bcf_var_header[6] >> 16;
      if (biallelic_strict && (n_allele > 2)) {
	goto bcf_to_bed_skip3;
      }
      if (n_allele > max_allele_ct) {
	goto bcf_to_bed_ret_NOMEM;
      }
      ujj = NON_BIGSTACK_MIN; // remaining allele name buffer space
      bufptr = allele_buf;
      if (n_allele) {
	for (uii = 0; uii < n_allele; uii++) {
	  retval = read_bcf_typed_string(gz_infile, bufptr, ujj, &ukk);
	  if (retval) {
	    goto bcf_to_bed_ret_1;
	  }
	  if ((!uii) && ((!ukk) || ((ukk == 1) && (*bufptr == 'N')))) {
	    // convert ref 'N' or '.' to missing genotype.  ('.' case was
	    // skipped the past, and 'N' was not converted.)
	    allele_lens[0] = 1;
	    allele_ptrs[0] = bufptr;
	    *bufptr++ = missing_geno;
	  } else {
	    allele_lens[uii] = ukk;
	    allele_ptrs[uii] = bufptr;
	    bufptr = &(bufptr[ukk]);
	  }
	}
      } else {
	// n_allele == 0 case was previously skipped, but it might have a place
	// with --allow-no-samples.
	allele_lens[0] = 1;
	allele_ptrs[0] = bufptr;
	*bufptr = missing_geno;
      }
      if (vcf_filter) {
	ii = gzgetc(gz_infile);
	if (ii == -1) {
	  goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
	} else {
	  ujj = ((uint32_t)ii) >> 4;
	  if (ujj == 15) {
	    retval = read_bcf_typed_nonnegative_integer(gz_infile, &ujj);
	    if (retval) {
	      goto bcf_to_bed_ret_1;
	    }
	  }
	  if (ujj) {
	    uii = ((uint32_t)ii) & 0x0f;
	    if ((uii < 1) || (uii > 3)) {
	      goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
	    }
	    if (uii == 1) {
	      if (ujj > 256) {
		goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
	      }
	      ucptr = (unsigned char*)g_textbuf;
	      if ((uint32_t)((uint64_t)gzread(gz_infile, ucptr, ujj)) < ujj) {
		goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
	      }
	      for (ukk = 0; ukk < ujj; ukk++) {
		if (ucptr[ukk] >= stringdict_ct) {
		  goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
		}
		if (!is_set(fexcept_bitfield, ucptr[ukk])) {
		  goto bcf_to_bed_marker_skip;
		}
	      }
	    } else if (uii == 2) {
	      if (ujj > 65536) {
		goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
	      }
	      ui16ptr = (uint16_t*)g_textbuf;
	      if ((uint32_t)((uint64_t)gzread(gz_infile, ui16ptr, ujj * sizeof(int16_t))) < ujj * sizeof(int16_t)) {
		goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
	      }
	      for (ukk = 0; ukk < ujj; ukk++) {
		if (ui16ptr[ukk] >= stringdict_ct) {
		  goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
		}
		if (!is_set(fexcept_bitfield, ui16ptr[ukk])) {
		  goto bcf_to_bed_marker_skip;
		}
	      }
	    } else {
	      // a bit more care required to avoid buffer overflow, if for some
	      // reason there are more than 32k filters...
	      uiptr = (uint32_t*)g_textbuf;
	      do {
		if (ujj > (MAXLINELEN / sizeof(int32_t))) {
		  ukk = MAXLINELEN / sizeof(int32_t);
		} else {
		  ukk = ujj;
		}
		if ((uint32_t)((uint64_t)gzread(gz_infile, uiptr, ukk * sizeof(int32_t))) < ukk * sizeof(int32_t)) {
		  goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
		}
		for (umm = 0; umm < ukk; umm++) {
		  if (uiptr[umm] >= stringdict_ct) {
		    goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
		  }
		  if (!is_set(fexcept_bitfield, uiptr[umm])) {
		    goto bcf_to_bed_marker_skip;
		  }
		}
		ujj -= ukk;
	      } while (ujj);
	    }
	  }
	}
      }
      alt_allele_idx = 1;
      if ((!gt_idx) || (!bcf_var_header[1])) {
	if (require_gt) {
	  goto bcf_to_bed_marker_skip;
	}
	ulljj = gztell(gz_infile);
	ullii = lastloc + bcf_var_header[0] + bcf_var_header[1];
	if (!sample_ct) {
	  goto bcf_to_bed_skip_genotype_write;
	}
	missing_gt_ct++;
	fill_quatervec_55(sample_ct, base_bitfields);
	goto bcf_to_bed_genotype_write;
      }

      // skip INFO
      ullii = lastloc + bcf_var_header[0];
      if (gzseek(gz_infile, ullii, SEEK_SET) == -1) {
	goto bcf_to_bed_ret_READ_FAIL;
      }
      ullii += bcf_var_header[1];
      while (1) {
	retval = read_bcf_typed_nonnegative_integer(gz_infile, &uii);
	if (retval) {
	  goto bcf_to_bed_ret_1;
	}
	ii = gzgetc(gz_infile);
	if (ii == -1) {
	  goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
	}
	ujj = ((uint32_t)ii) >> 4;
	if (ujj == 15) {
	  retval = read_bcf_typed_nonnegative_integer(gz_infile, &ujj);
	  if (retval) {
	    goto bcf_to_bed_ret_1;
	  }
	}
	if (ujj) {
	  ukk = ((uint32_t)ii) & 0x0f;
	  if ((ukk == 3) || (ukk == 5)) {
	    umm = 4; // int32, float = 4 bytes
	  } else if ((!ukk) || (ukk > 2)) {
	    logerrprint("Error: Unrecognized type in .bcf file.\n");
	    goto bcf_to_bed_ret_INVALID_FORMAT;
	  } else {
	    umm = ukk;
	  }
	}
	ulljj = gztell(gz_infile) + ((uint64_t)ujj) * umm * sample_ct;
	// uii = format key
	// ujj = for GT, max ploidy
	// ukk = integer/float/character type code
	// umm = bytes per entry
	if (ulljj > ullii) {
	  goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
	}
	if (uii == gt_idx) {
	  break;
	}
	// possible todo: --vcf-min-gq and --vcf-min-gp support
	if (ujj) {
	  if (gzseek(gz_infile, ((uint64_t)ujj) * umm * sample_ct, SEEK_CUR) == -1) {
	    goto bcf_to_bed_ret_READ_FAIL;
	  }
	  if (ulljj == ullii) {
	    if (require_gt) {
	      goto bcf_to_bed_marker_skip2;
	    } else {
	      missing_gt_ct++;
	      fill_quatervec_55(sample_ct, base_bitfields);
	      goto bcf_to_bed_genotype_write;
	    }
	  }
	}
      }
      if (!ujj) {
	// ploidy zero previously caused the variant to be skipped
	fill_quatervec_55(sample_ct, base_bitfields);
	goto bcf_to_bed_genotype_write;
      }
      if (ukk == 5) {
	logerrprint("Error: GT field cannot contain floating point values.\n");
	goto bcf_to_bed_ret_INVALID_FORMAT;
      }
      if (ujj * umm > 12) {
	// 12 = 12-ploid, or 6-ploid and >= 127 alleles, or triploid and >=
	// 32767 alleles.  this is pretty darn generous.
	logerrprint("Error: --bcf does not support GT vectors requiring >12 bytes per sample.\n");
	goto bcf_to_bed_ret_INVALID_FORMAT;
      }
      // ujj * umm <= 12 and sample_ct < 2^24, so no uint64_t cast needed there
      if ((uint32_t)((uint64_t)gzread(gz_infile, loadbuf, ujj * umm * sample_ct)) < ujj * umm * sample_ct) {
	goto bcf_to_bed_ret_READ_OR_FORMAT_FAIL;
      }
      if (n_allele < 2) {
	fill_ulong_zero(2 * sample_ctv2, base_bitfields);
      } else {
	fill_ulong_zero(n_allele * sample_ctv2, base_bitfields);
      }
      if (ukk == 1) {
	ucptr = (unsigned char*)loadbuf;
	if (ujj == 2) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, ucptr++) {
	    // discard all phase bits for now
	    // missing = 0x80 or 0x81
	    ulii = (*ucptr++) & 0x7e;
	    if (ulii) {
	      ulii = ((ulii / 2) - 1) * sample_ctv2;
	      uljj = (*ucptr) & 0x7e;
	      if (uljj) {
		set_bit(sample_idx * 2, &(base_bitfields[ulii]));
		base_bitfields[((uljj / 2) - 1) * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
	      } else {
		// could be MT or male X.  don't validate for now
		set_bit(sample_idx * 2 + 1, &(base_bitfields[ulii]));
	      }
	    }
	  }
	} else if (ujj == 1) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ulii = (*ucptr++) & 0x7e;
	    if (ulii) {
	      set_bit(sample_idx * 2 + 1, &(base_bitfields[((ulii / 2) - 1) * sample_ctv2]));
	    }
	  }
	} else {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    if (ucptr[2]) {
	      ucptr = &(ucptr[ujj]);
	    } else {
	      ulii = (*ucptr++) & 0x7e;
	      if (ulii) {
		ulii = ((ulii / 2) - 1) * sample_ctv2;
		uljj = (*ucptr) & 0x7e;
		if (uljj) {
		  set_bit(sample_idx * 2, &(base_bitfields[ulii]));
		  base_bitfields[((uljj / 2) - 1) * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
		} else {
		  set_bit(sample_idx * 2 + 1, &(base_bitfields[ulii]));
		}
	      }
	      ucptr = &(ucptr[ujj - 1]);
	    }
	  }
	}
      } else if (ukk == 2) {
	ui16ptr = (uint16_t*)loadbuf;
	// bleah, this should totally use templates instead of cut-and-paste
	if (ujj == 2) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, ui16ptr++) {
	    ulii = (*ui16ptr++) & 0x7ffe;
	    if (ulii) {
	      ulii = ((ulii / 2) - 1) * sample_ctv2;
	      uljj = (*ui16ptr) & 0x7ffe;
	      if (uljj) {
		set_bit(sample_idx * 2, &(base_bitfields[ulii]));
		base_bitfields[((uljj / 2) - 1) * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
	      } else {
		set_bit(sample_idx * 2 + 1, &(base_bitfields[ulii]));
	      }
	    }
	  }
	} else if (ujj == 1) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ulii = (*ui16ptr++) & 0x7ffe;
	    if (ulii) {
	      set_bit(sample_idx * 2 + 1, &(base_bitfields[((ulii / 2) - 1) * sample_ctv2]));
	    }
	  }
	} else {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    if (ui16ptr[2]) {
	      ui16ptr = &(ui16ptr[ujj]);
	    } else {
	      ulii = (*ui16ptr++) & 0x7ffe;
	      if (ulii) {
		ulii = ((ulii / 2) - 1) * sample_ctv2;
		uljj = (*ui16ptr) & 0x7ffe;
		if (uljj) {
		  set_bit(sample_idx * 2, &(base_bitfields[ulii]));
		  base_bitfields[((uljj / 2) - 1) * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
		} else {
		  set_bit(sample_idx * 2 + 1, &(base_bitfields[ulii]));
		}
	      }
	      ui16ptr = &(ui16ptr[ujj - 1]);
	    }
	  }
	}
      } else {
	uiptr = (uint32_t*)loadbuf;
	if (ujj == 2) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, uiptr++) {
	    ulii = (*uiptr++) & 0x7ffffffe;
	    if (ulii) {
	      ulii = ((ulii / 2) - 1) * sample_ctv2;
	      uljj = (*uiptr) & 0x7ffffffe;
	      if (uljj) {
		set_bit(sample_idx * 2, &(base_bitfields[ulii]));
		base_bitfields[((uljj / 2) - 1) * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
	      } else {
		set_bit(sample_idx * 2 + 1, &(base_bitfields[ulii]));
	      }
	    }
	  }
	} else if (ujj == 1) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ulii = (*uiptr++) & 0x7ffffffe;
	    if (ulii) {
	      set_bit(sample_idx * 2 + 1, &(base_bitfields[((ulii / 2) - 1) * sample_ctv2]));
	    }
	  }
	} else {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    if (uiptr[2]) {
	      uiptr = &(uiptr[ujj]);
	    } else {
	      ulii = (*uiptr++) & 0x7ffffffe;
	      if (ulii) {
		ulii = ((ulii / 2) - 1) * sample_ctv2;
		uljj = (*uiptr) & 0x7ffffffe;
		if (uljj) {
		  set_bit(sample_idx * 2, &(base_bitfields[ulii]));
		  base_bitfields[((uljj / 2) - 1) * sample_ctv2 + sample_idx / BITCT2] += ONELU << (2 * (sample_idx % BITCT2));
		} else {
		  set_bit(sample_idx * 2 + 1, &(base_bitfields[ulii]));
		}
	      }
	      uiptr = &(uiptr[ujj - 1]);
	    }
	  }
	}
      }
      if (n_allele > 2) {
	ulii = popcount2_longs(&(base_bitfields[sample_ctv2]), sample_ctl2);
	for (ulkk = 2; ulkk < n_allele; ulkk++) {
	  uljj = popcount2_longs(&(base_bitfields[sample_ctv2 * ulkk]), sample_ctl2);
	  if (!biallelic_only) {
	    if (uljj <= ulii) {
	      continue;
	    }
	  } else {
	    if (!uljj) {
	      continue;
	    }
	    if (ulii) {
	      goto bcf_to_bed_skip3;
	    }
	  }
	  ulii = uljj;
	  alt_allele_idx = ulkk;
	}
      }
      ref_ptr = base_bitfields;
      alt_ptr = &(base_bitfields[alt_allele_idx * sample_ctv2]);
      for (sample_idx = 0; sample_idx < sample_ctl2; sample_idx++) {
	ulii = *ref_ptr;
	uljj = *alt_ptr++;
	ulkk = (ulii + uljj) & AAAAMASK;
	uljj = ulii + ((ulii | (ulii >> 1)) & FIVEMASK);
	ulii = ulkk | (ulkk >> 1); // nonmissing?
	*ref_ptr++ = (uljj & ulii) | (((~ulkk) >> 1) & FIVEMASK);
      }
      ref_ptr[-1] &= final_mask;
    bcf_to_bed_genotype_write:
      if (fwrite_checked(base_bitfields, sample_ct4, outfile)) {
	goto bcf_to_bed_ret_WRITE_FAIL;
      }
    bcf_to_bed_skip_genotype_write:
      fputs(&(contigdict[bcf_var_header[2] * max_contig_len]), bimfile);
      putc_unlocked('\t', bimfile);
      if (marker_id_len) {
	fwrite(marker_id, 1, marker_id_len, bimfile);
      } else {
	putc_unlocked('.', bimfile);
      }
      // bcf2 coordinates are 0-based while vcf is 1-based... (seriously, whose
      // idea was this?  this is basically a bug in the spec due to how e.g.
      // telomeres are supposed to be encoded, but we have to play along)
      bufptr = uint32toa_x(bcf_var_header[3] + 1, '\t', &(tbuf2[3]));
      if (fwrite_checked(tbuf2, bufptr - tbuf2, bimfile)) {
	goto bcf_to_bed_ret_WRITE_FAIL;
      }
      if (n_allele > 1) {
	fwrite(allele_ptrs[alt_allele_idx], 1, allele_lens[alt_allele_idx], bimfile);
      } else {
	putc_unlocked(missing_geno, bimfile);
      }
      putc_unlocked('\t', bimfile);
      fwrite(allele_ptrs[0], 1, allele_lens[0], bimfile);
      if (putc_checked('\n', bimfile)) {
	goto bcf_to_bed_ret_WRITE_FAIL;
      }
      marker_ct++;
      if (!(marker_ct % 1000)) {
	printf("\r--bcf: %uk variants complete.", marker_ct / 1000);
	fflush(stdout);
      }
      if (ulljj < ullii) {
	if (gzseek(gz_infile, ullii, SEEK_SET) == -1) {
	  goto bcf_to_bed_ret_READ_FAIL;
	}
      }
      continue;
    bcf_to_bed_skip3:
      if (skip3_list) {
	if (!marker_skip_ct) {
	  memcpy(outname_end, ".skip.3allele", 14);
	  if (fopen_checked(outname, "w", &skip3file)) {
	    goto bcf_to_bed_ret_OPEN_FAIL;
	  }
	  memcpy(outname_end, ".bed", 5);
	}
	if (marker_id_len) {
	  fwrite(marker_id, 1, marker_id_len, skip3file);
	} else {
	  // up to the user to figure this out...
	  putc_unlocked('.', skip3file);
	}
	if (putc_checked('\n', skip3file)) {
	  goto bcf_to_bed_ret_OPEN_FAIL;
	}
      }
    bcf_to_bed_marker_skip:
      if (gzseek(gz_infile, (lastloc + bcf_var_header[0]) + bcf_var_header[1], SEEK_SET) == -1) {
	goto bcf_to_bed_ret_READ_FAIL;
      }
    bcf_to_bed_marker_skip2:
      marker_skip_ct++;
    }
    if ((!marker_ct) && (!allow_no_variants)) {
      if (marker_skip_ct) {
	logerrprint("Error: All variants in .bcf file skipped.\n");
	retval = RET_ALL_MARKERS_EXCLUDED;
	goto bcf_to_bed_ret_1;
      } else {
	logerrprint("Error: No variants in .bcf file.\n");
	goto bcf_to_bed_ret_INVALID_FORMAT;
      }
    }
    if (gzclose(gz_infile) != Z_OK) {
      gz_infile = nullptr;
      goto bcf_to_bed_ret_READ_FAIL;
    }
    gz_infile = nullptr;
    if (fclose_null(&bimfile)) {
      goto bcf_to_bed_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile)) {
      goto bcf_to_bed_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
  bcf_to_bed_skip_all_variants:
    *outname_end = '\0';
    LOGPRINTFWW("--bcf: %s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);
    if (marker_skip_ct) {
      LOGPRINTF("(%" PRIuPTR " variant%s skipped.)\n", marker_skip_ct, (marker_skip_ct == 1)? "" : "s");
    }
    if (missing_gt_ct) {
      LOGERRPRINTF("Warning: %" PRIuPTR " variant record%s had no GT field.\n", missing_gt_ct, (missing_gt_ct == 1)? "" : "s");
    }
  }
  while (0) {
  bcf_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  bcf_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  bcf_to_bed_ret_READ_OR_FORMAT_FAIL:
    if (gzeof(gz_infile)) {
      goto bcf_to_bed_ret_INVALID_FORMAT_GENERIC;
    }
  bcf_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  bcf_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  bcf_to_bed_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  bcf_to_bed_ret_INVALID_FORMAT_GENERIC:
    logerrprint("Error: Improperly formatted .bcf file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  bcf_to_bed_ret_INVALID_FORMAT_2:
    logprintb();
  bcf_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 bcf_to_bed_ret_1:
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  fclose_cond(bimfile);
  fclose_cond(skip3file);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return retval;
}

uint32_t write_23_cached_chrom(char* write_cache, uint32_t markers_left, char chrom_second_char, FILE* outfile_bed, FILE* outfile_bim) {
  // N.B. must flush writebuf before calling this
  uint32_t uii;
  do {
    uii = strlen(write_cache);
    if (putc_checked('2', outfile_bim)) {
      return 1;
    }
    putc_unlocked(chrom_second_char, outfile_bim);
    if (fwrite_checked(write_cache, uii, outfile_bim)) {
      return 1;
    }
    write_cache = &(write_cache[uii + 1]);
    if (putc_checked(*write_cache++, outfile_bed)) {
      return 1;
    }
  } while (--markers_left);
  return 0;
}

int32_t bed_from_23(char* infile_name, char* outname, char* outname_end, uint32_t modifier_23, char* fid_23, char* iid_23, double pheno_23, uint64_t misc_flags, char* paternal_id_23, char* maternal_id_23, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile_23 = nullptr;
  FILE* outfile_bed = nullptr;
  FILE* outfile_txt = nullptr;
  uintptr_t line_idx = 0;
  uint32_t is_male = modifier_23 & M23_MALE;
  uint32_t is_female = modifier_23 & M23_FEMALE;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  uint32_t x_present = 0;
  uint32_t haploid_x_present = 0;
  uint32_t y_present = 0;
  uint32_t nonmissing_y_present = 0;
  unsigned char* writebuf = (unsigned char*)(&(g_textbuf[MAXLINELEN]));
  int32_t retval = 0;
  uint32_t cur_chrom = 0;
  uint32_t chrom_mask_23 = (uint32_t)(chrom_info_ptr->chrom_mask[0]);
  uint32_t indel_ct = 0;
  unsigned char* writebuf_cur;
  unsigned char* writebuf_end;
  char* id_start;
  char* chrom_start;
  char* pos_start;
  char* allele_start;
  char* writebuf2;
  char* writebuf2_cur;
  uintptr_t id_len;
  uint32_t allele_calls;
  uint32_t null_chrom;
  uint32_t uii;
  char cc;
  char cc2;
  unsigned char ucc;
  {
    if (bigstack_alloc_c(MAXLINELEN, &writebuf2)) {
      goto bed_from_23_ret_NOMEM;
    }
    if (fopen_checked(infile_name, "r", &infile_23)) {
      goto bed_from_23_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(outname, "w", &outfile_txt)) {
      goto bed_from_23_ret_OPEN_FAIL;
    }
    memcpy(&(outname_end[2]), "ed", 2);
    if (fopen_checked(outname, FOPEN_WB, &outfile_bed)) {
      goto bed_from_23_ret_OPEN_FAIL;
    }
    if (bigstack_left() < MAXLINELEN) {
      goto bed_from_23_ret_NOMEM;
    }
    writebuf_cur = (unsigned char*)memcpyl3a((char*)writebuf, "l\x1b\x01");
    writebuf_end = &(writebuf[MAXLINELEN]);
    g_textbuf[MAXLINELEN - 1] = ' ';
    while (fgets(g_textbuf, MAXLINELEN, infile_23)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, infile_name);
	goto bed_from_23_ret_INVALID_FORMAT_2;
      }
      id_start = skip_initial_spaces(g_textbuf);
      if (is_eoln_or_comment_kns(*id_start)) {
	continue;
      }
      chrom_start = token_endnn(id_start);
      id_len = (uintptr_t)(chrom_start - id_start);
      chrom_start = skip_initial_spaces(chrom_start);
      pos_start = next_token(chrom_start);
      allele_start = next_token(pos_start);
      if (no_more_tokens_kns(allele_start)) {
	goto bed_from_23_ret_MISSING_TOKENS;
      }
      allele_calls = strlen_se(allele_start);
      if (allele_calls > 2) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has more allele calls than expected.\n", line_idx, infile_name);
	goto bed_from_23_ret_INVALID_FORMAT_2;
      }
      int32_t cur_chrom_code = get_chrom_code_raw(chrom_start);
      // X/Y/MT bugfix
      if (cur_chrom_code >= MAX_POSSIBLE_CHROM) {
	cur_chrom_code -= MAX_POSSIBLE_CHROM - 23;
      }
      if ((cur_chrom_code < 0) || (cur_chrom_code > 26)) {
	sprintf(g_logbuf, "Error: Invalid chromosome code on line %" PRIuPTR " of %s.\n", line_idx, infile_name);
	goto bed_from_23_ret_INVALID_FORMAT_2;
      }
      uii = (uint32_t)cur_chrom_code;
      if (!(chrom_mask_23 & (1 << uii))) {
	continue;
      }
      if (!uii) {
	null_chrom = 1;
      } else {
	if (uii < cur_chrom) {
	  LOGPREPRINTFWW("Error: Chromosomes in %s are out of order.\n", infile_name);
	  goto bed_from_23_ret_INVALID_FORMAT_2;
	} else if (uii > cur_chrom) {
	  cur_chrom = uii;
	  if (cur_chrom == 23) {
	    x_present = 1;
	  } else if (cur_chrom == 24) {
	    y_present = 1;
	  }
	}
	null_chrom = 0;
      }
      cc2 = allele_start[0];
      if ((cur_chrom == 24) && (!nonmissing_y_present)) {
	if (cc2 != '-') {
	  nonmissing_y_present = 1;
	}
      }
      if (cc2 == '-') {
	ucc = 1;
	cc = '0';
	cc2 = '0';
      } else if (allele_calls == 2) {
	cc = allele_start[1];
	if (cc == cc2) {
	  cc = '0';
	  ucc = 3;
	} else {
	  ucc = 2;
	}
      } else {
	cc = '0';
	ucc = 3;
      }
      if ((cc2 == 'D') || (cc2 == 'I')) {
	indel_ct++;
	cc = (char)(((unsigned char)cc2) ^ 13); // swaps D and I
      }
      if (!null_chrom) {
	if ((cur_chrom == 25) && (allele_calls != 2)) {
	  goto bed_from_23_ret_MISSING_ALLELE_CALLS;
	}
	if ((allele_calls == 1) && (cur_chrom <= 23)) {
	  if ((cur_chrom == 23) && (!is_female)) {
	    is_male = 1;
	    haploid_x_present = 1;
	  } else {
	    goto bed_from_23_ret_MISSING_ALLELE_CALLS;
	  }
	} else if ((cur_chrom == 24) && (cc2 != '0') && (!is_male)) {
	  if (!is_female) {
	    is_male = 1;
	  } else {
	    LOGPREPRINTFWW("Error: Nonmissing female allele call on line %" PRIuPTR " of %s.\n", line_idx, infile_name);
	    goto bed_from_23_ret_INVALID_FORMAT_2;
	  }
	}
      } else if (allele_calls == 1) {
	goto bed_from_23_ret_MISSING_ALLELE_CALLS;
      }
      if (!null_chrom) {
	writebuf2_cur = uint32toa(cur_chrom, writebuf2);
      } else {
	writebuf2[0] = '0';
	writebuf2_cur = &(writebuf2[1]);
      }
      *writebuf2_cur++ = '\t';
      writebuf2_cur = memcpya(writebuf2_cur, id_start, id_len);
      writebuf2_cur = memcpyl3a(writebuf2_cur, "\t0\t");
      writebuf2_cur = memcpyax(writebuf2_cur, pos_start, strlen_se(pos_start), '\t');
      *writebuf2_cur++ = cc;
      *writebuf2_cur++ = '\t';
      *writebuf2_cur++ = cc2;
      *writebuf2_cur++ = '\n';
      if (fwrite_checked(writebuf2, (uintptr_t)(writebuf2_cur - writebuf2), outfile_txt)) {
	goto bed_from_23_ret_WRITE_FAIL;
      }
      if (writebuf_cur == writebuf_end) {
	if (fwrite_checked(writebuf, (uintptr_t)(writebuf_cur - writebuf), outfile_bed)) {
	  goto bed_from_23_ret_WRITE_FAIL;
	}
	writebuf_cur = writebuf;
      }
      *writebuf_cur++ = (char)ucc;
    }
    if (!feof(infile_23)) {
      goto bed_from_23_ret_READ_FAIL;
    }
    if ((writebuf_cur == &(writebuf[3])) && (writebuf[0] == 'l') && (!allow_no_variants)) {
      if (chrom_mask_23 == 0x7ffffff) {
	logerrprint("Error: No --23file variants.\n");
	goto bed_from_23_ret_INVALID_FORMAT;
      } else {
	logerrprint("Error: No --23file variants pass chromosome filter.\n");
	goto bed_from_23_ret_INVALID_CMDLINE;
      }
    }
    if (fwrite_checked(writebuf, (uintptr_t)(writebuf_cur - writebuf), outfile_bed)) {
      goto bed_from_23_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile_txt)) {
      goto bed_from_23_ret_WRITE_FAIL;
    }
    memcpy(outname_end, ".fam", 5);
    if (fopen_checked(outname, "w", &outfile_txt)) {
      goto bed_from_23_ret_OPEN_FAIL;
    }
    if (fid_23) {
      fputs(fid_23, outfile_txt);
      putc_unlocked(' ', outfile_txt);
    } else {
      fputs("FAM001 ", outfile_txt);
    }
    if (iid_23) {
      fputs(iid_23, outfile_txt);
      putc_unlocked(' ', outfile_txt);
    } else {
      fputs("ID001 ", outfile_txt);
    }
    if (paternal_id_23) {
      fputs(paternal_id_23, outfile_txt);
      putc_unlocked(' ', outfile_txt);
    } else {
      fputs("0 ", outfile_txt);
    }
    if (maternal_id_23) {
      fputs(maternal_id_23, outfile_txt);
    } else {
      putc_unlocked('0', outfile_txt);
    }
    if (modifier_23 & M23_FORCE_MISSING_SEX) {
      cc = '0';
    } else if (is_male) {
      cc = '1';
    } else {
      cc = '2';
    }
    fprintf(outfile_txt, " %c %g\n", cc, pheno_23);
    if (fclose_null(&outfile_txt)) {
      goto bed_from_23_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    LOGPRINTFWW("--23file: %s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);
    if (indel_ct) {
      LOGPRINTF("%u variants with indel calls present.  '--snps-only no-DI' or\n--list-23-indels may be useful here.\n", indel_ct);
    }
    if (!(modifier_23 & M23_SEX)) {
      LOGPRINTF("Inferred sex: %smale.\n", is_male? "" : "fe");
    }
    if (modifier_23 & M23_MALE) {
      if (y_present && (!nonmissing_y_present)) {
	if (x_present) {
	  if (!haploid_x_present) {
	    logerrprint("Warning: No explicit haploid calls on X chromosome, and no nonmissing calls on\nY chromosome.  Double-check whether this is really a male sample.\n");
	  }
	} else {
	  logerrprint("Warning: No nonmissing calls on Y chromosome.  Double-check whether this is\nreally a male sample.\n");
	}
      }
    }
  }
  while (0) {
  bed_from_23_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  bed_from_23_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  bed_from_23_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  bed_from_23_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  bed_from_23_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  bed_from_23_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, infile_name);
    retval = RET_INVALID_FORMAT;
    break;
  bed_from_23_ret_MISSING_ALLELE_CALLS:
    LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer allele calls than expected.\n", line_idx, infile_name);
  bed_from_23_ret_INVALID_FORMAT_2:
    logerrprintb();
  bed_from_23_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile_23);
  fclose_cond(outfile_bed);
  fclose_cond(outfile_txt);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t generate_dummy(char* outname, char* outname_end, uint32_t flags, uintptr_t marker_ct, uintptr_t sample_ct, double geno_mrate, double pheno_mrate, int32_t missing_pheno) {
  FILE* outfile = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t sample_ct4 = (sample_ct + 3) / 4;
  uintptr_t urand = 0;
  double missing_phenod = (double)missing_pheno;
  uint32_t dbl_sample_mod4 = 2 * (sample_ct % 4);
  uint32_t four_alleles = 0;
  uint32_t geno_m_check = (geno_mrate >= RECIP_2_32 * 0.5);
  uint32_t geno_m32 = (uint32_t)(geno_mrate * 4294967296.0 - 0.5);
  uint32_t pheno_m_check = (pheno_mrate >= RECIP_2_32 * 0.5);
  uint32_t pheno_m32 = (uint32_t)(pheno_mrate * 4294967296.0 - 0.5);
  uint32_t saved_rnormal = 0;
  int32_t retval = 0;
  char wbuf[64];
  char missing_pheno_str[12];
  char* wptr = &(wbuf[5]);
  char* wptr2;
  double saved_rnormal_val;
  char alleles[13];
  unsigned char* writebuf;
  unsigned char* ucptr;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t pct;
  uint32_t loop_end;
  uint32_t missing_pheno_len;
  unsigned char ucc;
  unsigned char ucc2;
  uint64_t ullii;
  double dxx;
  wptr2 = int32toa(missing_pheno, missing_pheno_str);
  missing_pheno_len = (uintptr_t)(wptr2 - missing_pheno_str);
  *wptr2 = '\0';
  if (flags & DUMMY_ACGT) {
    memcpy(alleles, "ACAGATCGCTGTA", 13);
    four_alleles = 1;
  } else if (flags & DUMMY_1234) {
    memcpy(alleles, "1213142324341", 13);
    four_alleles = 1;
  } else if (flags & DUMMY_12) {
    memcpyl3(alleles, "121");
  } else {
    memcpyl3(alleles, "ABA");
  }
  if (bigstack_alloc_uc(sample_ct4, &writebuf)) {
    goto generate_dummy_ret_NOMEM;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
  memcpy(wbuf, "1\tsnp", 5);
  if (four_alleles) {
    for (uii = 0; uii < marker_ct; uii++) {
      if (!(uii % 8)) {
	do {
	  urand = sfmt_genrand_uint32(&g_sfmt);
	} while (urand < 425132032LU); // 2^32 - 12^8.  heck, why not
      }
      ukk = urand / 12U;
      ujj = urand - (ukk * 12U);
      urand = ukk;
      wptr2 = uint32toa(uii, wptr);
      wptr2 = memcpyl3a(wptr2, "\t0\t");
      wptr2 = uint32toa_x(uii, '\t', wptr2);
      wptr2[0] = alleles[ujj];
      wptr2[1] = '\t';
      wptr2[2] = alleles[ujj + 1];
      wptr2[3] = '\n';
      if (fwrite_checked(wbuf, 4 + (uintptr_t)(wptr2 - wbuf), outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
  } else {
    for (uii = 0; uii < marker_ct; uii++) {
      if (!(uii % 32)) {
	urand = sfmt_genrand_uint32(&g_sfmt);
      }
      ujj = urand & 1;
      urand >>= 1;
      wptr2 = uint32toa(uii, wptr);
      wptr2 = memcpyl3a(wptr2, "\t0\t");
      wptr2 = uint32toa_x(uii, '\t', wptr2);
      wptr2[0] = alleles[ujj];
      wptr2[1] = '\t';
      wptr2[2] = alleles[ujj + 1];
      wptr2[3] = '\n';
      if (fwrite_checked(wbuf, 4 + (uintptr_t)(wptr2 - wbuf), outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto generate_dummy_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
  wptr = memcpyl3a(wbuf, "per");
  if (flags & DUMMY_SCALAR_PHENO) {
    for (uii = 0; uii < sample_ct; uii++) {
      if (pheno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= pheno_m32)) {
	dxx = missing_phenod;
      } else {
	if (saved_rnormal) {
	  dxx = saved_rnormal_val;
	  saved_rnormal = 0;
	} else {
	  dxx = rand_normal(&saved_rnormal_val);
	  saved_rnormal = 1;
	}
      }
      wptr2 = uint32toa(uii, wptr);
      wptr2 = memcpya(wptr2, " per", 4);
      wptr2 = uint32toa(uii, wptr2);
      wptr2 = memcpya(wptr2, " 0 0 2 ", 7);
      wptr2 = dtoa_gx(dxx, '\n', wptr2);
      if (fwrite_checked(wbuf, wptr2 - wbuf, outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
  } else {
    for (uii = 0; uii < sample_ct; uii++) {
      if (!(uii % 32)) {
	urand = sfmt_genrand_uint32(&g_sfmt);
      }
      wptr2 = uint32toa(uii, wptr);
      wptr2 = memcpya(wptr2, " per", 4);
      wptr2 = uint32toa(uii, wptr2);
      wptr2 = memcpya(wptr2, " 0 0 2 ", 7);
      if (pheno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= pheno_m32)) {
	wptr2 = memcpya(wptr2, missing_pheno_str, missing_pheno_len);
      } else {
	*wptr2++ = (char)((urand & 1) + '1');
      }
      *wptr2++ = '\n';
      if (fwrite_checked(wbuf, wptr2 - wbuf, outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
      urand >>= 1;
    }
  }
  if (fclose_null(&outfile)) {
    goto generate_dummy_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(outname, FOPEN_WB, &outfile)) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto generate_dummy_ret_WRITE_FAIL;
  }
  uii = 0;
  ullii = (3 * ONELU) + ((uint64_t)marker_ct) * sample_ct4;
  if (ullii >= 10485760) {
    printf("Writing dummy .bed (%" PRIu64 " MB)... 0%%", ullii >> 20);
  } else {
    fputs("Writing dummy .bed... 0%", stdout);
  }
  fflush(stdout);
  for (pct = 1; pct <= 100; pct++) {
    loop_end = ((uint64_t)(pct * marker_ct)) / 100U;
    for (; uii < loop_end; uii++) {
      ucptr = writebuf;
      for (ujj = 0; ujj < sample_ct4; ujj++) {
	if (!(ujj % 4)) {
	  urand = sfmt_genrand_uint32(&g_sfmt);
	}
	ucc = 0;
	for (ukk = 0; ukk < 8; ukk += 2) {
	  if (geno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= geno_m32)) {
	    ucc2 = 1;
	  } else {
	    ucc2 = urand & 3;
	    if (ucc2 == 1) {
	      ucc2 = 2;
	    }
	  }
	  ucc |= ucc2 << ukk;
	  urand >>= 2;
	}
	*ucptr++ = ucc;
      }
      if (dbl_sample_mod4) {
	ucc = *(--ucptr);
	*ucptr = ucc >> (8 - dbl_sample_mod4);
      }

      ujj = popcount_chars((uintptr_t*)writebuf, 0, sample_ct4);
      if (ujj < sample_ct) {
	reverse_loadbuf(sample_ct, writebuf);
      }
      if (fwrite_checked(writebuf, sample_ct4, outfile)) {
	putc_unlocked('\n', stdout);
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
    if (pct < 100) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  putc_unlocked('\r', stdout);
  *outname_end = '\0';
  LOGPRINTFWW("Dummy data (%" PRIuPTR " %s, %" PRIuPTR " SNP%s) written to %s.bed + %s.bim + %s.fam .\n", sample_ct, species_str(sample_ct), marker_ct, (marker_ct == 1)? "" : "s", outname, outname, outname);
  while (0) {
  generate_dummy_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  generate_dummy_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  generate_dummy_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
  }
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

void simulate_init_freqs_qt(uint32_t tags_or_haps, double dprime, double qt_var, double qt_dom, double missing_freq, double* freqs, uint64_t* thresholds, double* qt_adj) {
  // Initialize frequency table for current SNP.  Similar to instanceSNP_QT()
  // in PLINK 1.07 simul.cpp.  (The difference is that heterozygote frequencies
  // are not stored in a redundant manner.)
  // thresholds[0]: P(causal 11, marker 11) * 2^63
  // thresholds[1]: P(causal 11, marker 11 or missing) * 2^63
  // thresholds[2]: P(causal 11, marker 11 or missing or 12) * 2^63
  // ...
  // thresholds[14]: [1 - P(causal 22, marker 22)] * 2^63
  // We use 2^63 instead of 2^64 to avoid integer overflow headaches.
  double freq = freqs[0];
  double mfreq = freqs[1];
  double nonmissing_freq = 1 - missing_freq;
  double scaled_miss_freq = missing_freq * TWO_63;
  double scaled_nm_freq = nonmissing_freq * TWO_63;
  double ld = freq * (1 - mfreq);

  double qq = 1 - freq;
  double aa = sqrt(qt_var / (2 * freq * qq * ((1 + qt_dom * (qq - freq)) * (1 + qt_dom * (qq - freq)) + qt_dom * 2 * freq * qq * qt_dom)));
  double dd = qt_dom * aa;
  double dxx = aa * (1 - 2 * freq * (1 + qq * qt_dom));

  // Joint causal variant/marker probs, considering one allele at a time.
  // First digit indicates causal allele, second digit indicates marker.
  double h21 = qq * mfreq;
  double h11;
  double h12;
  double h22;

  // first 2 digits indicate causal variant genotype, last 2 indicate marker.
  // this is DIFFERENT from instanceSNP()'s usage.
  double h_11_11;
  double h_11_12;
  double h_11_22;
  double h_12_11;
  double h_12_12;
  double h_12_22;
  double h_22_11;
  double h_22_12;
  double h_22_22;

  // Constraints:
  // 1. mean = (qq^2) * qt_adj[3] + 2 * freq * qq * qt_adj[2] +
  //           (freq^2) * qt_adj[0]
  //         = 0
  // 2. variance = (qq^2) * (qt_adj[3]^2) + 2 * freq * qq * (qt_adj[2]^2) +
  //               (freq^2) * (qt_adj[0]^2)
  //             = qt_var
  // 3. qt_adj[2] = ((qt_adj[0] + qt_adj[3]) / 2) +
  //                qt_dom * ((qt_adj[0] - qt_adj[3]) / 2)
  // PLINK 1.07 computes the correct values, but incorrectly reverses their
  // positions, and also incorrectly uses sp[s].gAB in some places where .gBB
  // should be used and vice versa.
  qt_adj[0] = dxx + aa;
  qt_adj[2] = dxx + dd;
  qt_adj[3] = dxx - aa;

  if (h21 < ld) {
    ld = h21;
  }
  ld *= dprime;
  h11 = freq * mfreq + ld;
  h12 = freq * (1 - mfreq) - ld;
  h21 -= ld;
  h22 = qq * (1 - mfreq) + ld;
  h_11_11 = h11 * h11;
  h_11_12 = h11 * h12 * 2;
  h_11_22 = h12 * h12;
  h_12_11 = h21 * h11 * 2;
  h_12_12 = (h22 * h11 + h21 * h12) * 2;
  h_12_22 = h22 * h12 * 2;
  h_22_11 = h21 * h21;
  h_22_12 = h22 * h21 * 2;
  h_22_22 = h22 * h22;
  // determination of causal variant missing/nonmissing status is deferred, to
  // simplify phenotype updating
  if (!tags_or_haps) {
    thresholds[0] = (uint64_t)((h_11_11 + h_12_11 + h_22_11) * TWO_63);
    thresholds[1] = thresholds[0] + (uint64_t)((h_11_12 + h_12_12 + h_22_12) * TWO_63);
  } else {
    thresholds[0] = (uint64_t)(h_11_11 * scaled_nm_freq);
    thresholds[1] = thresholds[0] + (uint64_t)((h_11_11 + h_11_12 + h_11_22) * scaled_miss_freq);
    thresholds[2] = thresholds[1] + (uint64_t)(h_11_12 * scaled_nm_freq);
    thresholds[3] = thresholds[2] + (uint64_t)(h_11_22 * scaled_nm_freq);
    thresholds[4] = thresholds[3] + (uint64_t)(h_12_11 * scaled_nm_freq);
    thresholds[5] = thresholds[4] + (uint64_t)((h_12_11 + h_12_12 + h_12_22) * scaled_miss_freq);
    thresholds[6] = thresholds[5] + (uint64_t)(h_12_12 * scaled_nm_freq);
    thresholds[7] = thresholds[6] + (uint64_t)(h_12_22 * scaled_nm_freq);
    thresholds[8] = thresholds[7] + (uint64_t)(h_22_11 * scaled_nm_freq);
    thresholds[9] = thresholds[8] + (uint64_t)((h_22_11 + h_22_12 + h_22_22) * scaled_miss_freq);
    thresholds[10] = thresholds[9] + (uint64_t)(h_22_12 * scaled_nm_freq);
  }
}

void simulate_cc_get_conditional_probs(double prevalence, double g0, double g1, double g2, double het_odds, double hom0_odds, double* f0p, double* f1p, double* f2p) {
  // PLINK 1.07 interprets het_odds and hom0_odds as probability ratios instead
  // of odds ratios.  The two are nearly identical if prevalence is small, but
  // it's still better to avoid that approximation.  This requires solving a
  // cubic equation in X := odds(hom2).
  //
  // prevalence = g0 * (hom0_odds * X) / (1 + hom0_odds * X)
  //            + g1 * (het_odds * X) / (1 + het_odds * X)
  //            + g2 * X / (1 + X)
  //
  // prevalence * (1 + X) * (1 + het_odds * X) * (1 + hom0_odds * X) =
  //   X * (g0 * hom0_odds * (1 + het_odds * X) * (1 + X) +
  //        g1 * het_odds * (1 + hom0_odds * X) * (1 + X) +
  //        g2 * (1 + hom0_odds * X) * (1 + het_odds * X))
  //
  //   X^3 * prevalence * het_odds * hom0_odds
  // + X^2 * prevalence * (het_odds * hom0_odds + het_odds + hom0_odds)
  // + X * prevalence * (1 + het_odds + hom0_odds)
  // + prevalence =
  //   X^3 * het_odds * hom0_odds
  // + X^2 * (g0 * hom0_odds * (1 + het_odds) +
  //          g1 * het_odds * (1 + hom0_odds) +
  //          g2 * (hom0_odds + het_odds))
  // + X * (g0 * hom0_odds + g1 * het_odds + g2)
  //
  // 0 = X^3 * het_odds * hom0_odds * (1 - prevalence)
  //   + X^2 * (g0 * hom0_odds * (1 + het_odds) +
  //            g1 * het_odds * (1 + hom0_odds) +
  //            g2 * (hom0_odds + het_odds) -
  //            prevalence * (het_odds * hom0_odds + het_odds + hom0_odds))
  //   + X * (g0 * hom0_odds + g1 * het_odds + g2 -
  //          prevalence * (1 + het_odds + hom0_odds)
  //   - prevalence
  double solutions[3];
  double coef_recip;
  double cur_f2_odds;
  double cur_f0_odds;
  double cur_f1_odds;
  uint32_t root_ct;
  uint32_t root_idx;
  if ((prevalence == 0) || (prevalence == 1)) {
    *f0p = prevalence;
    *f1p = prevalence;
    *f2p = prevalence;
    return;
  }
  coef_recip = 1.0 / (het_odds * hom0_odds * (1.0 - prevalence));
  // this always has a positive solution since f(0) is negative
  root_ct = cubic_real_roots(coef_recip * (g0 * hom0_odds * (1 + het_odds) + g1 * het_odds + (1 + hom0_odds) + g2 * (hom0_odds + het_odds) - prevalence * (het_odds * hom0_odds + het_odds + hom0_odds)), coef_recip * (g0 * hom0_odds + g1 * het_odds + g2 - prevalence * (1 + het_odds + hom0_odds)), coef_recip * (-prevalence), solutions);
  cur_f2_odds = solutions[0];
  root_idx = 0;
  while ((cur_f2_odds <= 0) && (root_idx + 1 < root_ct)) {
    cur_f2_odds = solutions[++root_idx];
  }
  // odds = p / (1 - p)
  // -> p = odds / (1 + odds)
  cur_f0_odds = cur_f2_odds * hom0_odds;
  cur_f1_odds = cur_f2_odds * het_odds;
  *f0p = cur_f0_odds / (1 + cur_f0_odds);
  *f1p = cur_f1_odds / (1 + cur_f1_odds);
  *f2p = cur_f2_odds / (1 + cur_f2_odds);
}

void simulate_init_freqs_cc(uint32_t do_haps, double dprime, double* freqs, double prevalence, double het_odds, double hom0_odds, double missing_freq, uint64_t* ctrl_thresholds, uint64_t* case_thresholds) {
  // Similar to instanceSNP().
  double freq = freqs[0];
  double mfreq = freqs[1];
  double nonmissing_freq = 1 - missing_freq;
  double scaled_missing_freq = missing_freq * TWO_63;
  double scaled_nonmissing_freq = nonmissing_freq * TWO_63;

  // P(causal variant genotype)
  double g0 = freq * freq;
  double g1 = 2 * freq * (1 - freq);
  double g2 = 1 - g0 - g1;
  // P(marker genotype)
  double mg0 = mfreq * mfreq;
  double mg1 = 2 * mfreq * (1 - mfreq);
  double mg2 = 1 - mg0 - mg1;

  double ld = freq * (1 - mfreq);
  double h21 = (1 - freq) * mfreq;
  double h11;
  double h12;
  double h22;

  double f0;
  double f1;
  double f2;
  double mf0;
  double mf1;
  double mf2;
  // first 2 digits indicate causal variant genotype, last 2 indicate marker.
  // this is DIFFERENT from instanceSNP()'s usage.
  double h_11_11;
  double h_11_12;
  double h_11_22;
  double h_12_11;
  double h_12_12;
  double h_12_22;
  double h_22_11;
  double h_22_12;
  double h_22_22;

  double a0;
  double a1;
  double a2;
  double u0;
  double u1;
  double u2;

  double xh_11_11;
  double xh_11_12;
  double xh_11_22;
  double xh_12_11;
  double xh_12_12;
  double xh_12_22;
  double xh_22_11;
  double xh_22_12;
  double xh_22_22;

  double tot_recip;
  double tot_recip_nmsq;
  double tot_recip_1miss;

  if (h21 < ld) {
    ld = h21;
  }
  ld *= dprime;

  // Joint causal variant/marker probs, considering one allele at a time.
  h11 = freq * mfreq + ld;
  h12 = freq * (1 - mfreq) - ld;
  h21 -= ld;
  h22 = (1 - freq) * (1 - mfreq) + ld;


  // Now considering both alleles simultaneously.
  h_11_11 = h11 * h11;
  h_11_12 = h11 * h12 * 2;
  h_11_22 = h12 * h12;
  h_12_11 = h21 * h11 * 2;
  h_12_12 = (h22 * h11 + h21 * h12) * 2;
  h_12_22 = h22 * h12 * 2;
  h_22_11 = h21 * h21;
  h_22_12 = h22 * h21 * 2;
  h_22_22 = h22 * h22;

  // P(case | causal variant genotype)
  simulate_cc_get_conditional_probs(prevalence, g0, g1, g2, het_odds, hom0_odds, &f0, &f1, &f2);
  // P(case | marker genotype)
  mf0 = (f0 * h_11_11 + f1 * h_12_11 + f2 * h_22_11) / mg0;
  mf1 = (f0 * h_11_12 + f1 * h_12_12 + f2 * h_22_12) / mg1;
  mf2 = (f0 * h_11_22 + f1 * h_12_22 + f2 * h_22_22) / mg2;
  // P(marker genotype | affection status)
  a0 = mg0 * mf0;
  a1 = mg1 * mf1;
  a2 = mg2 * mf2;
  tot_recip = 1.0 / (a0 + a1 + a2);
  a0 *= tot_recip;
  a1 *= tot_recip;
  a2 *= tot_recip;
  u0 = mg0 * (1 - mf0);
  u1 = mg1 * (1 - mf1);
  u2 = mg2 * (1 - mf2);
  tot_recip = 1.0 / (u0 + u1 + u2);
  u0 *= tot_recip;
  u1 *= tot_recip;
  u2 *= tot_recip;

  if (!do_haps) {
    case_thresholds[0] = (uint64_t)(a0 * scaled_nonmissing_freq);
    case_thresholds[1] = case_thresholds[0] + (uint64_t)(scaled_missing_freq);
    case_thresholds[2] = case_thresholds[1] + (uint64_t)(a1 * scaled_nonmissing_freq);
    ctrl_thresholds[0] = (uint64_t)(u0 * scaled_nonmissing_freq);
    ctrl_thresholds[1] = ctrl_thresholds[0] + (uint64_t)(scaled_missing_freq);
    ctrl_thresholds[2] = ctrl_thresholds[1] + (uint64_t)(u1 * scaled_nonmissing_freq);
  } else {
    xh_11_11 = h_11_11 * f0;
    xh_11_12 = h_11_12 * f0;
    xh_11_22 = h_11_22 * f0;
    xh_12_11 = h_12_11 * f1;
    xh_12_12 = h_12_12 * f1;
    xh_12_22 = h_12_22 * f1;
    xh_22_11 = h_22_11 * f2;
    xh_22_12 = h_22_12 * f2;
    xh_22_22 = h_22_22 * f2;
    tot_recip = 1.0 / (xh_11_11 + xh_11_12 + xh_11_22 + xh_12_11 + xh_12_12 + xh_12_22 + xh_22_11 + xh_22_12 + xh_22_22);
    tot_recip_nmsq = tot_recip * nonmissing_freq * scaled_nonmissing_freq;
    tot_recip_1miss = tot_recip * nonmissing_freq * scaled_missing_freq;

    case_thresholds[0] = (uint64_t)(xh_11_11 * tot_recip_nmsq);
    case_thresholds[1] = case_thresholds[0] + (uint64_t)((xh_11_11 + xh_11_12 + xh_11_22) * tot_recip_1miss);
    case_thresholds[2] = case_thresholds[1] + (uint64_t)(xh_11_12 * tot_recip_nmsq);
    case_thresholds[3] = case_thresholds[2] + (uint64_t)(xh_11_22 * tot_recip_nmsq);
    case_thresholds[4] = case_thresholds[3] + (uint64_t)((xh_11_11 + xh_12_11 + xh_22_11) * tot_recip_1miss);
    case_thresholds[5] = case_thresholds[4] + (uint64_t)(missing_freq * scaled_missing_freq);
    case_thresholds[6] = case_thresholds[5] + (uint64_t)((xh_11_12 + xh_12_12 + xh_22_12) * tot_recip_1miss);
    case_thresholds[7] = case_thresholds[6] + (uint64_t)((xh_11_22 + xh_12_22 + xh_22_22) * tot_recip_1miss);
    case_thresholds[8] = case_thresholds[7] + (uint64_t)(xh_12_11 * tot_recip_nmsq);
    case_thresholds[9] = case_thresholds[8] + (uint64_t)((xh_12_11 + xh_12_12 + xh_12_22) * tot_recip_1miss);
    case_thresholds[10] = case_thresholds[9] + (uint64_t)(xh_12_12 * tot_recip_nmsq);
    case_thresholds[11] = case_thresholds[10] + (uint64_t)(xh_12_22 * tot_recip_nmsq);
    case_thresholds[12] = case_thresholds[11] + (uint64_t)(xh_22_11 * tot_recip_nmsq);
    case_thresholds[13] = case_thresholds[12] + (uint64_t)((xh_22_11 + xh_22_12 + xh_22_22) * tot_recip_1miss);
    case_thresholds[14] = case_thresholds[13] + (uint64_t)(xh_22_12 * tot_recip_nmsq);

    xh_11_11 = h_11_11 * (1 - f0);
    xh_11_12 = h_11_12 * (1 - f0);
    xh_11_22 = h_11_22 * (1 - f0);
    xh_12_11 = h_12_11 * (1 - f1);
    xh_12_12 = h_12_12 * (1 - f1);
    xh_12_22 = h_12_22 * (1 - f1);
    xh_22_11 = h_22_11 * (1 - f2);
    xh_22_12 = h_22_12 * (1 - f2);
    xh_22_22 = h_22_22 * (1 - f2);
    tot_recip = 1.0 / (xh_11_11 + xh_11_12 + xh_11_22 + xh_12_11 + xh_12_12 + xh_12_22 + xh_22_11 + xh_22_12 + xh_22_22);
    tot_recip_nmsq = tot_recip * nonmissing_freq * scaled_nonmissing_freq;
    tot_recip_1miss = tot_recip * nonmissing_freq * scaled_missing_freq;

    ctrl_thresholds[0] = (uint64_t)(xh_11_11 * tot_recip_nmsq);
    ctrl_thresholds[1] = ctrl_thresholds[0] + (uint64_t)((xh_11_11 + xh_11_12 + xh_11_22) * tot_recip_1miss);
    ctrl_thresholds[2] = ctrl_thresholds[1] + (uint64_t)(xh_11_12 * tot_recip_nmsq);
    ctrl_thresholds[3] = ctrl_thresholds[2] + (uint64_t)(xh_11_22 * tot_recip_nmsq);
    ctrl_thresholds[4] = ctrl_thresholds[3] + (uint64_t)((xh_11_11 + xh_12_11 + xh_22_11) * tot_recip_1miss);
    ctrl_thresholds[5] = ctrl_thresholds[4] + (uint64_t)(missing_freq * scaled_missing_freq);
    ctrl_thresholds[6] = ctrl_thresholds[5] + (uint64_t)((xh_11_12 + xh_12_12 + xh_22_12) * tot_recip_1miss);
    ctrl_thresholds[7] = ctrl_thresholds[6] + (uint64_t)((xh_11_22 + xh_12_22 + xh_22_22) * tot_recip_1miss);
    ctrl_thresholds[8] = ctrl_thresholds[7] + (uint64_t)(xh_12_11 * tot_recip_nmsq);
    ctrl_thresholds[9] = ctrl_thresholds[8] + (uint64_t)((xh_12_11 + xh_12_12 + xh_12_22) * tot_recip_1miss);
    ctrl_thresholds[10] = ctrl_thresholds[9] + (uint64_t)(xh_12_12 * tot_recip_nmsq);
    ctrl_thresholds[11] = ctrl_thresholds[10] + (uint64_t)(xh_12_22 * tot_recip_nmsq);
    ctrl_thresholds[12] = ctrl_thresholds[11] + (uint64_t)(xh_22_11 * tot_recip_nmsq);
    ctrl_thresholds[13] = ctrl_thresholds[12] + (uint64_t)((xh_22_11 + xh_22_12 + xh_22_22) * tot_recip_1miss);
    ctrl_thresholds[14] = ctrl_thresholds[13] + (uint64_t)(xh_22_12 * tot_recip_nmsq);
  }
}

int32_t simulate_dataset(char* outname, char* outname_end, uint32_t flags, char* simulate_fname, uint32_t case_ct, uint32_t ctrl_ct, double prevalence, uint32_t sample_ct, double missing_freq, char* name_prefix) {
  FILE* infile = nullptr;
  FILE* outfile_txt = nullptr;
  FILE* outfile_simfreq = nullptr;
  FILE* outfile_bed = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  double* qt_vals = nullptr;
  char* cur_snp_label = &(g_textbuf[MAXLINELEN]);
  char* marker_freq_lb_ptr = nullptr;
  char* marker_ld_ptr = nullptr;
  uintptr_t* writebuf2 = nullptr;
  double dxx = 0;
  double dyy = 0;
  double qt_totvar = 0;
  uint32_t cur_marker_ct = 0;
  double marker_freq_lb = 0;
  double marker_freq_ub = 0;
  double dprime = 0;
  double het_odds = 0;
  double hom0_odds = 0;
  double qt_var = 0;
  double qt_dom = 0;
  uintptr_t line_idx = 0;
  uint32_t do_tags = flags & SIMULATE_TAGS;
  uint32_t do_haps = flags & SIMULATE_HAPS;
  uint32_t tags_or_haps = do_tags | do_haps;
  uint32_t is_qt = flags & SIMULATE_QT;
  uint32_t randomize_alleles = flags & (SIMULATE_ACGT | SIMULATE_1234 | SIMULATE_12);
  uint32_t simulate_12 = flags & SIMULATE_12;
  uint32_t marker_pos = 1;
  uint32_t snp_label_len = 0;
  uint32_t missing_thresh = (uint32_t)(missing_freq * 4294967296.0);
  uint32_t marker_idx_offset = 0;
  uint32_t pct = 0;
  uint32_t name_prefix_len = 0;
  uint32_t zero_odds_ratio_warning_given = 0;
  int32_t retval = 0;
  char alleles[13];
  char cur_alleles[4];
  double freqs[6];
  uint64_t thresholds[15];
  uint64_t case_thresholds[15];
  double qt_adj[4];
  uintptr_t* writebuf;
  char* cptr;
  char* snp_label_ptr;
  char* freq_lb_ptr;
  char* penult_ptr; // het disease odds for C/C, additive genetic var for QT
  char* last_ptr; // homset disease odds for C/C, dominance/additive for QT
  char* wptr;
  uintptr_t* wbptr;
  uintptr_t* wbptr2;
  sfmt_t* sfmt64p;
  double freq_lb;
  double freq_delta;
  double dzz;
  uint64_t ullii;
  uintptr_t sample_ct4;
  uintptr_t sample_ctl2;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t marker_ct;
  uint32_t sample_idx;
  uint32_t cur_marker_idx;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  char cc;
  if (!is_qt) {
    sample_ct = case_ct + ctrl_ct;
  } else {
    if (bigstack_calloc_d(sample_ct, &qt_vals)) {
      goto simulate_ret_NOMEM;
    }
  }
  sample_ct4 = (sample_ct + 3) / 4;
  sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  if (randomize_alleles) {
    if (flags & SIMULATE_ACGT) {
      memcpy(alleles, "ACAGATCGCTGTA", 13);
    } else if (flags & SIMULATE_1234) {
      memcpy(alleles, "1213142324341", 13);
    } else {
      memcpyl3(alleles, "121");
    }
  } else {
    if (is_qt) {
      memcpy(alleles, "HLAB", 4);
    } else {
      memcpy(alleles, "DdAB", 4);
    }
  }
  if (bigstack_alloc_ul(sample_ctl2, &writebuf)) {
    goto simulate_ret_NOMEM;
  }
  if (do_haps) {
    if (bigstack_alloc_ul(sample_ctl2, &writebuf2)) {
      goto simulate_ret_NOMEM;
    }
  }
  if (fopen_checked(simulate_fname, "r", &infile)) {
    goto simulate_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(outname, "w", &outfile_txt)) {
    goto simulate_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".simfreq", 9);
  if (fopen_checked(outname, "w", &outfile_simfreq)) {
    goto simulate_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(outname, FOPEN_WB, &outfile_bed)) {
    goto simulate_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile_bed)) {
    goto simulate_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  LOGPRINTFWW5("Writing --simulate%s dataset to %s.bed + %s.bim + %s.fam ... ", is_qt? "-qt" : "", outname, outname, outname);
  fputs("0%", stdout);
  fflush(stdout);
  sfmt64p = (sfmt_t*)bigstack_alloc(sizeof(sfmt_t));
  if (!sfmt64p) {
    goto simulate_ret_NOMEM;
  }
  init_sfmt64_from_sfmt32(&g_sfmt, sfmt64p);
  g_textbuf[MAXLINELEN - 1] = ' ';
  // just determine total marker ct in initial scan, for progress indicator
  ullii = 0;
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --simulate%s file is pathologically long.\n", line_idx, is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2N;
    }
    cptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*cptr)) {
      continue;
    }
    if (scan_uint_icap(cptr, &uii)) {
      sprintf(g_logbuf, "Error: Invalid SNP count on line %" PRIuPTR " of --simulate%s input file.\n", line_idx, is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2N;
    }
    ullii += uii;
  }
  if (!feof(infile)) {
    goto simulate_ret_READ_FAIL;
  }
  if (!ullii) {
    sprintf(g_logbuf, "Error: --simulate%s input file specifies zero SNPs.\n", is_qt? "-qt" : "");
    goto simulate_ret_INVALID_FORMAT_2N;
  } else if (ullii > (do_haps? 0x3ffffffe : 0x7ffffffd)) {
    sprintf(g_logbuf, "Error: --simulate%s input file specifies too many SNPs.\n", is_qt? "-qt" : "");
    goto simulate_ret_INVALID_FORMAT_2N;
  }
  marker_ct = ullii;
  loop_end = (marker_ct + 99) / 100;
  rewind(infile);
  line_idx = 0;
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    // already checked for long lines, don't need to repeat
    cptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*cptr)) {
      continue;
    }
    snp_label_ptr = next_token(cptr);
    freq_lb_ptr = next_token(snp_label_ptr);
    if (tags_or_haps) {
      marker_freq_lb_ptr = next_token_mult(freq_lb_ptr, 2);
      marker_ld_ptr = next_token_mult(marker_freq_lb_ptr, 2);
      penult_ptr = next_token(marker_ld_ptr);
    } else {
      penult_ptr = next_token_mult(freq_lb_ptr, 2);
    }
    last_ptr = next_token(penult_ptr);
    if (no_more_tokens_kns(last_ptr)) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --simulate%s file has fewer tokens than expected.\n", line_idx, is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2N;
    }
    if (!no_more_tokens_kns(next_token(last_ptr))) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --simulate%s file has more tokens than expected.\n", line_idx, is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2N;
    }
    scan_uint_icap(cptr, &uii);
    if (!uii) {
      continue;
    }
    cur_marker_ct = uii;
    snp_label_len = strlen_se(snp_label_ptr);
    memcpy(cur_snp_label, snp_label_ptr, snp_label_len);
    cur_snp_label[snp_label_len++] = '_';
    if (scan_two_doubles(freq_lb_ptr, &freq_lb, &freq_delta) || (freq_lb < 0) || (freq_delta < freq_lb) || (freq_delta > 1)) {
      sprintf(g_logbuf, "Error: Invalid allele frequency bound on line %" PRIuPTR " of --simulate%s\nfile.\n", line_idx, is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2N;
    }
    freq_delta -= freq_lb;
    if (tags_or_haps) {
      if (scan_two_doubles(marker_freq_lb_ptr, &marker_freq_lb, &marker_freq_ub) || (marker_freq_lb < 0) || (marker_freq_ub < marker_freq_lb) || (marker_freq_ub > 1)) {
	sprintf(g_logbuf, "Error: Invalid marker allele frequency bound on line %" PRIuPTR " of\n--simulate%s file.\n", line_idx, is_qt? "-qt" : "");
	goto simulate_ret_INVALID_FORMAT_2N;
      }
      if (scan_double(marker_ld_ptr, &dprime) || (dprime < 0) || (dprime > 1)) {
	sprintf(g_logbuf, "Error: Invalid d-prime on line %" PRIuPTR " of --simulate%s input file.\n", line_idx, is_qt? "-qt" : "");
	goto simulate_ret_INVALID_FORMAT_2N;
      }
    } else {
      dprime = 1;
    }
    if (is_qt) {
      if (scan_double(penult_ptr, &qt_var) || (qt_var < 0) || (qt_var > 1)) {
	sprintf(g_logbuf, "Error: Invalid variance value on line %" PRIuPTR " of --simulate-qt file.\n", line_idx);
	goto simulate_ret_INVALID_FORMAT_2N;
      }
      if ((qt_var > 0) && (((freq_delta == 0) && ((freq_lb == 0) || (freq_lb == 1))) || (tags_or_haps && (marker_freq_lb == marker_freq_ub) && ((marker_freq_lb == 0) || (marker_freq_lb == 1))))) {
	sprintf(g_logbuf, "Error: Nonzero variance with fixed 0/1 allele frequency on line %" PRIuPTR " of\n--simulate-qt file.\n", line_idx);
	goto simulate_ret_INVALID_FORMAT_2N;
      }
      qt_totvar += ((intptr_t)cur_marker_ct) * qt_var;
      if (qt_totvar > 1 + EPSILON) {
	logprint("\n");
	logerrprint("Error: --simulate-qt input file specific QTL variance greater than 1.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
      if (scan_double(last_ptr, &qt_dom)) {
	sprintf(g_logbuf, "Error: Invalid dominance deviation value on line %" PRIuPTR " of --simulate-qt\nfile.\n", line_idx);
	goto simulate_ret_INVALID_FORMAT_2N;
      }
    } else {
      if (scan_double(penult_ptr, &het_odds) || (het_odds < 0)) {
	sprintf(g_logbuf, "Error: Invalid heterozygote disease odds ratio on line %" PRIuPTR " of\n--simulate file.\n", line_idx);
	goto simulate_ret_INVALID_FORMAT_2N;
      }
      if ((strlen_se(last_ptr) == 4) && match_upper_counted(last_ptr, "MULT", 4)) {
	hom0_odds = het_odds * het_odds;
      } else if (scan_double(last_ptr, &hom0_odds) || (hom0_odds < 0)) {
	sprintf(g_logbuf, "Error: Invalid homozygote disease odds ratio on line %" PRIuPTR " of --simulate\nfile.\n", line_idx);
	goto simulate_ret_INVALID_FORMAT_2N;
      }
      if ((!zero_odds_ratio_warning_given) && ((het_odds == 0) || (hom0_odds == 0))) {
        putc_unlocked('\r', stdout);
	logstr("\n");
	logerrprint("Warning: Zero odds ratio present in --simulate input file.  Did you mean\n--simulate-qt instead?\n");
	zero_odds_ratio_warning_given = 1;
        printf("%u%%", pct);
      }
    }
    g_textbuf[0] = '1';
    for (cur_marker_idx = 0; cur_marker_idx < cur_marker_ct; cur_marker_idx++) {
      freqs[0] = freq_lb + rand_unif() * freq_delta;
      if (tags_or_haps) {
	// force marker freq/c.v. freq and (1 - marker freq)/(1 - c.v. freq) to
	// not be smaller than dprime or larger than 1 / dprime, unless that
	// causes marker freq to be out of range.
	if (dprime > 0) {
	  dxx = 1 - (1 - freqs[0]) / dprime;
	  dzz = freqs[0] * dprime;
	  if (dxx < dzz) {
	    dxx = dzz;
	  }
	  if (dxx < marker_freq_lb) {
	    dxx = marker_freq_lb;
	  }
	  dyy = 1 - (1 - freqs[0]) * dprime;
	  dzz = freqs[0] / dprime;
	  if (dyy > dzz) {
	    dyy = dzz;
	  }
	  if (dyy > marker_freq_ub) {
	    dyy = marker_freq_ub;
	  }
	  if (dyy < dxx) {
	    if (dyy < marker_freq_lb) {
	      dyy = marker_freq_lb;
	    }
	    dxx = dyy;
	  }
	} else {
	  dxx = marker_freq_lb;
	  dyy = marker_freq_ub;
	}
	freqs[1] = dxx + rand_unif() * (dyy - dxx);
      } else {
	freqs[1] = freqs[0];
      }
      if (is_qt) {
	simulate_init_freqs_qt(do_haps, dprime, qt_var, qt_dom, missing_freq, freqs, thresholds, qt_adj);
      } else {
	simulate_init_freqs_cc(do_haps, dprime, freqs, prevalence, het_odds, hom0_odds, missing_freq, thresholds, case_thresholds);
      }
      wptr = &(g_textbuf[1]);
      *wptr++ = ' ';
      if (cur_marker_ct > 1) {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len);
	wptr = uint32toa(cur_marker_idx, wptr);
      } else {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len - 1);
      }
      *wptr++ = '\t';
      dxx = freqs[0];
      wptr = dtoa_gx(dxx, ' ', wptr);
      wptr = dtoa_gx(dxx, '\t', wptr);
      if (tags_or_haps) {
	dxx = freqs[1];
	wptr = dtoa_gx(dxx, ' ', wptr);
	wptr = dtoa_gx(dxx, '\t', wptr);
	wptr = dtoa_gx(dprime, '\t', wptr);
      }
      if (is_qt) {
	wptr = dtoa_gx(qt_var, '\t', wptr);
	wptr = dtoa_gx(qt_dom, '\n', wptr);
      } else {
	wptr = dtoa_gx(het_odds, '\t', wptr);
	wptr = dtoa_gx(hom0_odds, '\n', wptr);
      }
      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_simfreq)) {
	goto simulate_ret_WRITE_FAIL;
      }
      if (randomize_alleles) {
	if (!simulate_12) {
	  do {
	    uii = sfmt_genrand_uint32(&g_sfmt);
	  } while (uii >= 4294967184U); // largest multiple of 144 < 2^32
	  uii = uii % 144U;
	  ujj = uii / 12;
	  uii -= ujj * 12;
	} else {
	  uii = sfmt_genrand_uint32(&g_sfmt) & 3;
	  ujj = uii >> 1;
	  uii &= 1;
	}
	memcpy(cur_alleles, &(alleles[uii]), 2);
	memcpy(&(cur_alleles[2]), &(alleles[ujj]), 2);
      } else {
	memcpy(cur_alleles, alleles, 4);
      }
      ulii = 0;
      ukk = 0;
      wbptr = writebuf;
      if (!do_haps) {
	if (is_qt) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    if (ullii > thresholds[1]) {
	      ulkk = 3;
	    } else if (ullii > thresholds[0]) {
	      ulkk = 2;
	    } else {
	      ulkk = 0;
	    }
	    qt_vals[sample_idx] += qt_adj[ulkk];
	    if (sfmt_genrand_uint32(&g_sfmt) < missing_thresh) {
	      ulkk = 1;
	    }
	    ulii |= ulkk << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      ulii = 0;
	      ukk = 0;
	    }
	  }
	} else {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    if (sample_idx < case_ct) {
	      if (ullii > case_thresholds[1]) {
		if (ullii > case_thresholds[2]) {
		  ulkk = 3;
		} else {
		  ulkk = 2;
		}
	      } else if (ullii > case_thresholds[0]) {
		ulkk = 1;
	      } else {
		ulkk = 0;
	      }
	    } else {
	      if (ullii > thresholds[1]) {
		if (ullii > thresholds[2]) {
		  ulkk = 3;
		} else {
		  ulkk = 2;
		}
	      } else if (ullii > thresholds[0]) {
		ulkk = 1;
	      } else {
		ulkk = 0;
	      }
	    }
	    ulii |= ulkk << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      ulii = 0;
	      ukk = 0;
	    }
	  }
	}
	if (ukk) {
	  *wbptr = ulii;
	}
      } else {
	uljj = 0;
	wbptr2 = writebuf2;
	if (is_qt) {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    ulkk = uint64arr_greater_than(thresholds, 11, ullii);
	    ulmm = ulkk & 3;
	    ulkk /= 4;
	    ulkk += (ulkk + 1) >> 1;
	    qt_vals[sample_idx] += qt_adj[ulkk];
	    if (sfmt_genrand_uint32(&g_sfmt) < missing_thresh) {
	      ulkk = 1;
	    }
	    ulii |= ulkk << ukk;
	    uljj |= ulmm << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      *wbptr2++ = uljj;
	      ulii = 0;
	      uljj = 0;
	      ukk = 0;
	    }
	  }
	} else {
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    if (sample_idx < case_ct) {
	      ulkk = uint64arr_greater_than(case_thresholds, 15, ullii);
	    } else {
	      ulkk = uint64arr_greater_than(thresholds, 15, ullii);
	    }
	    ulmm = ulkk & 3;
	    ulkk /= 4;
	    ulii |= ulkk << ukk;
	    uljj |= ulmm << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      *wbptr2++ = uljj;
	      ulii = 0;
	      uljj = 0;
	      ukk = 0;
	    }
	  }
	}
	if (ukk) {
	  *wbptr = ulii;
	  *wbptr2 = uljj;
	}
      }
      if (popcount_longs(writebuf, sample_ctl2) < sample_ct) {
	reverse_loadbuf(sample_ct, (unsigned char*)writebuf);
	cc = cur_alleles[0];
	cur_alleles[0] = cur_alleles[1];
	cur_alleles[1] = cc;
      }
      wptr = &(g_textbuf[1]);
      *wptr++ = '\t';
      if (cur_marker_ct > 1) {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len);
	wptr = uint32toa(cur_marker_idx, wptr);
      } else {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len - 1);
      }
      if (do_tags) {
	wptr = memcpya(wptr, "_M", 2);
      }
      wptr = memcpyl3a(wptr, "\t0\t");
      wptr = uint32toa_x(marker_pos++, '\t', wptr);
      *wptr++ = cur_alleles[0];
      *wptr++ = '\t';
      *wptr++ = cur_alleles[1];
      *wptr++ = '\n';
      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_txt)) {
	goto simulate_ret_WRITE_FAIL;
      }
      if (fwrite_checked((unsigned char*)writebuf, sample_ct4, outfile_bed)) {
	goto simulate_ret_WRITE_FAIL;
      }
      if (do_haps) {
	if (popcount_longs(writebuf2, sample_ctl2) < sample_ct) {
	  reverse_loadbuf(sample_ct, (unsigned char*)writebuf2);
	  cc = cur_alleles[2];
	  cur_alleles[2] = cur_alleles[3];
	  cur_alleles[3] = cc;
	}
	wptr = &(g_textbuf[2 + snp_label_len]);
	wptr = uint32toa(cur_marker_idx, wptr);
	wptr = memcpya(wptr, "_M\t0\t", 5);
	wptr = uint32toa_x(marker_pos++, '\t', wptr);
	*wptr++ = cur_alleles[2];
	*wptr++ = '\t';
	*wptr++ = cur_alleles[3];
	*wptr++ = '\n';
	if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_txt)) {
	  goto simulate_ret_WRITE_FAIL;
	}
	if (fwrite_checked((unsigned char*)writebuf2, sample_ct4, outfile_bed)) {
	  goto simulate_ret_WRITE_FAIL;
	}
      }
      if (cur_marker_idx >= loop_end) {
	if (pct > 9) {
	  putc_unlocked('\b', stdout);
	}
	pct = ((cur_marker_idx + marker_idx_offset) * 100LLU) / marker_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = ((marker_ct * (pct + 1LLU) + 99LLU) / 100U) - marker_idx_offset;
      }
    }
    marker_idx_offset += cur_marker_ct;
    loop_end -= cur_marker_ct;
  }
  if (!feof(infile)) {
    goto simulate_ret_READ_FAIL;
  }
  if (fclose_null(&outfile_txt) || fclose_null(&outfile_bed)) {
    goto simulate_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(outname, "w", &outfile_txt)) {
    goto simulate_ret_OPEN_FAIL;
  }
  wptr = g_textbuf;
  if (name_prefix) {
    name_prefix_len = strlen(name_prefix);
    wptr = memcpyax(wptr, name_prefix, name_prefix_len, '-');
    uii = 4 + name_prefix_len;
  } else {
    uii = 3;
  }
  memcpyl3(wptr, "per");
  if (is_qt) {
    if (qt_totvar < 1 - EPSILON) {
      dyy = sqrt(1 - qt_totvar);
    } else {
      dyy = 0;
    }
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    wptr = uint32toa_x(sample_idx, ' ', &(g_textbuf[uii]));
    if (name_prefix_len) {
      wptr = memcpyax(wptr, name_prefix, name_prefix_len, '-');
    }
    wptr = memcpyl3a(wptr, "per");
    wptr = uint32toa(sample_idx, wptr);
    wptr = memcpya(wptr, " 0 0 2 ", 7);
    if (is_qt) {
      if (sample_idx & 1) {
	dzz = qt_vals[sample_idx] + dyy * dxx;
      } else {
	dzz = qt_vals[sample_idx] + dyy * rand_normal(&dxx);
      }
      wptr = dtoa_g(dzz, wptr);
    } else {
      if (sample_idx < case_ct) {
	*wptr++ = '2';
      } else {
	*wptr++ = '1';
      }
    }
    *wptr++ = '\n';
    if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_txt)) {
      goto simulate_ret_WRITE_FAIL;
    }
  }
  *outname_end = '\0';
  if (pct > 9) {
    putc_unlocked('\b', stdout);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  LOGPRINTFWW("Realized simulation parameters saved to %s.simfreq.\n", outname);
  while (0) {
  simulate_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  simulate_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  simulate_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  simulate_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  simulate_ret_INVALID_FORMAT_2N:
    logprint("\n");
    logerrprintb();
  simulate_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  fclose_cond(outfile_txt);
  fclose_cond(outfile_simfreq);
  fclose_cond(outfile_bed);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t recode_allele_load(char* loadbuf, uintptr_t loadbuf_size, char* recode_allele_name, char*** allele_missing_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* recode_allele_reverse, char* recode_allele_extra) {
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* rafile = nullptr;
  uint32_t missing_allele = 0;
  uint32_t marker_id_htable_size = get_id_htable_size(marker_ct);
  uintptr_t rae_size = 0;
  uintptr_t line_idx = 0;
  uintptr_t cur_bigstack_left;
  uint32_t* marker_id_htable;
  char* bufptr;
  char* bufptr2;
  int32_t retval;
  uint32_t slen;
  uint32_t alen;
  uintptr_t marker_uidx;
  if (fopen_checked(recode_allele_name, "r", &rafile)) {
    goto recode_allele_load_ret_OPEN_FAIL;
  }
  if (bigstack_end_alloc_ui(marker_id_htable_size, &marker_id_htable)) {
    goto recode_allele_load_ret_NOMEM;
  }
  retval = populate_id_htable(unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, marker_id_htable_size, marker_id_htable);
  if (retval) {
    goto recode_allele_load_ret_1;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  cur_bigstack_left = bigstack_left();
  while (fgets(loadbuf, loadbuf_size, rafile)) {
    line_idx++;
    if (!loadbuf[loadbuf_size - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --recode-allele file is pathologically long.\n", line_idx);
      goto recode_allele_load_ret_INVALID_FORMAT_3;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    slen = strlen_se(bufptr);
    bufptr2 = skip_initial_spaces(&(bufptr[slen]));
    if (is_eoln_kns(*bufptr2)) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --recode-allele file has fewer tokens than expected.\n", line_idx);
      goto recode_allele_load_ret_INVALID_FORMAT_3;
    }
    alen = strlen_se(bufptr2);
    marker_uidx = id_htable_find(bufptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx != 0xffffffffU) {
      bufptr2[alen++] = '\0';
      if (!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx])) {
	CLEAR_BIT(marker_uidx, recode_allele_reverse);
      } else if (!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx + 1])) {
	SET_BIT(marker_uidx, recode_allele_reverse);
      } else {
	if (rae_size + alen > cur_bigstack_left) {
	  goto recode_allele_load_ret_NOMEM;
	}
	missing_allele = 1;
	(*allele_missing_ptr)[marker_uidx] = &(recode_allele_extra[rae_size]);
	memcpy(&(recode_allele_extra[rae_size]), bufptr2, alen);
	rae_size += alen;
      }
    }
  }
  if (!feof(rafile)) {
    goto recode_allele_load_ret_READ_FAIL;
  }
  while (0) {
  recode_allele_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  recode_allele_load_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  recode_allele_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  recode_allele_load_ret_INVALID_FORMAT_3:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
  }
 recode_allele_load_ret_1:
  fclose_cond(rafile);
  bigstack_end_reset(bigstack_end_mark);
  if (missing_allele) {
    recode_allele_extra = (char*)bigstack_alloc(rae_size);
  } else {
    bigstack_reset(*allele_missing_ptr);
    *allele_missing_ptr = nullptr;
  }
  return retval;
}

uint32_t recode_load_to(unsigned char* loadbuf, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t marker_idx, uintptr_t marker_idx_end, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t* marker_uidx_ptr, uintptr_t unfiltered_sample_ct) {
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uintptr_t marker_uidx_start;
  uintptr_t marker_uidx_stop;
  uintptr_t ulii;
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    return 1;
  }
  while (marker_idx < marker_idx_end) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	return 1;
      }
    }
    if (unfiltered_marker_ct - marker_uidx > marker_idx_end - marker_idx) {
      ulii = next_set_ul_unsafe(marker_exclude, marker_uidx) - marker_uidx;
    } else {
      ulii = unfiltered_marker_ct - marker_uidx;
    }
    marker_uidx_start = marker_uidx;
    marker_uidx_stop = marker_uidx + ulii;
    marker_idx += ulii;
    ulii *= unfiltered_sample_ct4;
    if (fread(loadbuf, 1, ulii, bedfile) < ulii) {
      return 1;
    }
    while (1) {
      next_set_ul_ck(marker_reverse, marker_uidx_stop, &marker_uidx);
      if (marker_uidx == marker_uidx_stop) {
	break;
      }
      reverse_loadbuf(unfiltered_sample_ct, &(loadbuf[(marker_uidx - marker_uidx_start) * unfiltered_sample_ct4]));
      marker_uidx++;
    }
    loadbuf = &(loadbuf[ulii]);
  }
  *marker_uidx_ptr = marker_uidx;
  return 0;
}

static inline int32_t recode_write_first_cols(FILE* outfile, uintptr_t sample_uidx, char delimiter, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, const char* output_missing_pheno) {
  char wbuf[16];
  char* cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
  uintptr_t ulii = strlen_se(cptr);
  fwrite(cptr, 1, ulii, outfile);
  putc_unlocked(delimiter, outfile);
  fputs(&(cptr[ulii + 1]), outfile);
  putc_unlocked(delimiter, outfile);
  fputs(paternal_ids? (&(paternal_ids[sample_uidx * max_paternal_id_len])) : "0", outfile);
  putc_unlocked(delimiter, outfile);
  fputs(maternal_ids? (&(maternal_ids[sample_uidx * max_maternal_id_len])) : "0", outfile);
  putc_unlocked(delimiter, outfile);
  putc_unlocked(sexchar(sex_nm, sex_male, sample_uidx), outfile);
  putc_unlocked(delimiter, outfile);
  if (!IS_SET(pheno_nm, sample_uidx)) {
    fputs(output_missing_pheno, outfile);
  } else if (pheno_c) {
    putc_unlocked('1' + IS_SET(pheno_c, sample_uidx), outfile);
  } else {
    cptr = dtoa_g(pheno_d[sample_uidx], wbuf);
    fwrite(wbuf, 1, cptr - wbuf, outfile);
  }
  if (putc_checked(delimiter, outfile)) {
    return -1;
  }
  return 0;
}

void init_recode_cmax0(char* allele1, char* allele2, char** cur_mk_allelesx, uint32_t* cmalen, char delimiter, char delim2) {
  uintptr_t alen1 = strlen(allele1);
  uintptr_t alen2 = strlen(allele2);
  char* cmaptr = cur_mk_allelesx[0];
  *cmaptr++ = delimiter;
  cmaptr = memcpyax(cmaptr, allele1, alen1, delim2);
  memcpy(cmaptr, allele1, alen1);
  cmalen[0] = alen1 * 2 + 2;
  cmaptr = cur_mk_allelesx[2];
  *cmaptr++ = delimiter;
  cmaptr = memcpyax(cmaptr, allele1, alen1, delim2);
  memcpy(cmaptr, allele2, alen2);
  cmalen[2] = alen1 + alen2 + 2;
  cmaptr = cur_mk_allelesx[3];
  *cmaptr++ = delimiter;
  cmaptr = memcpyax(cmaptr, allele2, alen2, delim2);
  memcpy(cmaptr, allele2, alen2);
  cmalen[3] = alen2 * 2 + 2;
}

void init_recode_cmax(char* allele1, char* allele2, char** cur_mk_allelesx, uint32_t* cmalen, char delimiter, char delim2) {
  uintptr_t alen1 = strlen(allele1);
  uintptr_t alen2 = strlen(allele2);
  char* cmaptr = memcpyax(cur_mk_allelesx[0], allele1, alen1, delim2);
  memcpyx(cmaptr, allele1, alen1, delimiter);
  cmalen[0] = alen1 * 2 + 2;
  cmaptr = memcpyax(cur_mk_allelesx[2], allele1, alen1, delim2);
  memcpyx(cmaptr, allele2, alen2, delimiter);
  cmalen[2] = alen1 + alen2 + 2;
  cmaptr = memcpyax(cur_mk_allelesx[3], allele2, alen2, delim2);
  memcpyx(cmaptr, allele2, alen2, delimiter);
  cmalen[3] = alen2 * 2 + 2;
}

void init_cur_mk_allelesx(char* mk_alleles, uintptr_t max_marker_allele_len, uint32_t do_reverse, char** cur_mk_allelesx, uint32_t* cmalen, char delimiter, char delim2) {
  uint32_t alen = strlen_se(mk_alleles);
  memcpy(cur_mk_allelesx[0], mk_alleles, alen);
  cmalen[0] = alen;
  if (delimiter) {
    cmalen[0] += 1;
    cur_mk_allelesx[0][alen] = delim2;
    memcpy(cur_mk_allelesx[1], mk_alleles, alen);
    cur_mk_allelesx[1][alen] = delimiter;
  } else {
    cur_mk_allelesx[0][alen] = '\0';
  }
  mk_alleles = &(mk_alleles[max_marker_allele_len]);
  alen = strlen_se(mk_alleles);
  cmalen[1] = alen;
  if (delimiter) {
    cmalen[1] += 1;
    memcpy(cur_mk_allelesx[2], mk_alleles, alen);
    cur_mk_allelesx[2][alen] = delim2;
    memcpy(cur_mk_allelesx[3], mk_alleles, alen);
    cur_mk_allelesx[3][alen] = delimiter;
    if (do_reverse) {
      cur_mk_allelesx[4] = cur_mk_allelesx[2];
      cur_mk_allelesx[5] = cur_mk_allelesx[1];
    } else {
      cur_mk_allelesx[4] = cur_mk_allelesx[0];
      cur_mk_allelesx[5] = cur_mk_allelesx[3];
    }
  } else {
    memcpy(cur_mk_allelesx[1], mk_alleles, alen + 1);
    if (do_reverse) {
      cur_mk_allelesx[4] = cur_mk_allelesx[1];
      cur_mk_allelesx[5] = cur_mk_allelesx[0];
    } else {
      cur_mk_allelesx[4] = cur_mk_allelesx[0];
      cur_mk_allelesx[5] = cur_mk_allelesx[1];
    }
  }
  if (do_reverse) {
    cmalen[2] = cmalen[1];
    cmalen[3] = cmalen[0];
  } else {
    cmalen[2] = cmalen[0];
    cmalen[3] = cmalen[1];
  }
}

int32_t recode_beagle_new_chrom(char* outname, char* outname_end2, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uintptr_t* marker_uidx_ptr, uint32_t* chrom_fo_idx_ptr, uint32_t* chrom_idx_ptr, uint32_t* chrom_end_ptr, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_sample_ct4, FILE** datfile_ptr, FILE** mapfile_ptr, char* dat_header, uintptr_t dat_header_len) {
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uint32_t chrom_fo_idx = *chrom_fo_idx_ptr;
  int32_t retval = 0;
  uint32_t chrom_end;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t is_haploid;
  char* wbufptr;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
  if ((chrom_idx > chrom_info_ptr->autosome_ct) && (chrom_idx <= chrom_info_ptr->max_code)) {
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, chrom_end);
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    } while ((chrom_idx > chrom_info_ptr->autosome_ct) && (chrom_idx <= chrom_info_ptr->max_code));
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
      goto recode_beagle_new_chrom_ret_READ_FAIL;
    }
  }
  *marker_uidx_ptr = marker_uidx;
  *chrom_fo_idx_ptr = chrom_fo_idx;
  *chrom_idx_ptr = chrom_idx;
  *chrom_end_ptr = chrom_end;
  if (!mapfile_ptr) {
    return 0;
  }

  wbufptr = chrom_name_write(chrom_info_ptr, chrom_idx, outname_end2);
  memcpy(wbufptr, ".dat", 5);
  if (fopen_checked(outname, "w", datfile_ptr)) {
    goto recode_beagle_new_chrom_ret_OPEN_FAIL;
  }
  memcpy(wbufptr, ".map", 5);
  if (fopen_checked(outname, "w", mapfile_ptr)) {
    goto recode_beagle_new_chrom_ret_OPEN_FAIL;
  }
  if (fwrite_checked(dat_header, dat_header_len, *datfile_ptr)) {
    goto recode_beagle_new_chrom_ret_WRITE_FAIL;
  }
  *wbufptr = '\0';
  LOGPREPRINTFWW("%s.dat + %s.map created.\n", outname, outname);
  logstr(g_logbuf);
  while (0) {
  recode_beagle_new_chrom_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  recode_beagle_new_chrom_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  recode_beagle_new_chrom_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  return retval;
}

int32_t open_and_write_fastphase_header(FILE** outfile_ptr, char* outname, uintptr_t* marker_exclude, uint32_t* marker_pos, uint32_t marker_uidx, uint32_t chrom_size, uint32_t sample_ct) {
  char wbuf[16];
  char* wptr;
  uint32_t marker_idx;
  if (fopen_checked(outname, "w", outfile_ptr)) {
    return RET_OPEN_FAIL;
  }
  wptr = uint32toa_x(sample_ct, '\n', wbuf);
  if (fwrite_checked(wbuf, wptr - wbuf, *outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  wptr = uint32toa(chrom_size, wbuf);
  fwrite(wbuf, 1, wptr - wbuf, *outfile_ptr);
  fputs("\nP ", *outfile_ptr);
  for (marker_idx = 0; marker_idx < chrom_size; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    wptr = uint32toa_x(marker_pos[marker_uidx], ' ', wbuf);
    fwrite(wbuf, 1, wptr - wbuf, *outfile_ptr);
  }
  if (putc_checked('\n', *outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  LOGPREPRINTFWW("%s created.\n", outname);
  logstr(g_logbuf);
  return 0;
}

uint32_t write_ped_lines(FILE* outfile, unsigned char* loadbuf, uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, uintptr_t unfiltered_sample_ct4, uintptr_t* sample_exclude, uintptr_t* sample_uidx_ptr, uintptr_t sample_idx_start, uintptr_t sample_idx_end, uint32_t recode_compound, char** mk_allele_ptrs, char delimiter, char delim2, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char output_missing_geno, const char* output_missing_pheno, char* writebuf) {
  uintptr_t sample_uidx = *sample_uidx_ptr;
  char missing4[4];
  unsigned char* bufptr;
  char* wbufptr;
  char* aptr;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t sample_idx;
  uint32_t shiftval;
  uint32_t alen;
  unsigned char ucc;
  missing4[0] = output_missing_geno;
  missing4[1] = delim2;
  missing4[2] = output_missing_geno;
  missing4[3] = delimiter;
  for (sample_idx = sample_idx_start; sample_idx < sample_idx_end; sample_uidx++, sample_idx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    if (recode_write_first_cols(outfile, sample_uidx, delimiter, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno)) {
      return 1;
    }
    bufptr = &(loadbuf[sample_uidx / 4]);
    wbufptr = writebuf;
    shiftval = (sample_uidx % 4) * 2;
    marker_uidx = marker_uidx_start;
    if (recode_compound) {
      for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	ucc = ((*bufptr) >> shiftval) & 3;
	if (ucc) {
	  if (ucc == 3) {
	    *wbufptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	    *wbufptr = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  } else if (ucc == 2) {
	    *wbufptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	    *wbufptr = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  } else {
	    *wbufptr++ = output_missing_geno;
	    *wbufptr = output_missing_geno;
	  }
	} else {
	  *wbufptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	  *wbufptr = mk_allele_ptrs[2 * marker_uidx][0];
	}
	bufptr = &(bufptr[unfiltered_sample_ct4]);
	wbufptr = &(wbufptr[2]);
      }
      if (fwrite_checked(writebuf, marker_ct * 3, outfile)) {
	return 1;
      }
    } else {
      for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	ucc = ((*bufptr) >> shiftval) & 3;
	if (ucc) {
	  if (ucc == 3) {
	    aptr = mk_allele_ptrs[2 * marker_uidx + 1];
	    alen = strlen(aptr);
	    memcpy(wbufptr, aptr, alen);
	    wbufptr[alen] = delim2;
	    memcpy(&(wbufptr[alen + 1]), aptr, alen);
	    wbufptr[2 * alen + 1] = delimiter;
	    wbufptr = &(wbufptr[2 * alen + 2]);
	  } else if (ucc == 2) {
	    wbufptr = strcpyax(wbufptr, mk_allele_ptrs[2 * marker_uidx], delim2);
	    wbufptr = strcpyax(wbufptr, mk_allele_ptrs[2 * marker_uidx + 1], delimiter);
	  } else {
	    wbufptr = memcpya(wbufptr, missing4, 4);
	  }
	} else {
	  aptr = mk_allele_ptrs[2 * marker_uidx];
	  alen = strlen(aptr);
	  memcpy(wbufptr, aptr, alen);
	  wbufptr[alen] = delim2;
	  memcpy(&(wbufptr[alen + 1]), aptr, alen);
	  wbufptr[2 * alen + 1] = delimiter;
	  wbufptr = &(wbufptr[2 * alen + 2]);
	}
	bufptr = &(bufptr[unfiltered_sample_ct4]);
      }
      if (marker_ct) {
	wbufptr[-1] = '\n';
	if (fwrite_checked(writebuf, wbufptr - writebuf, outfile)) {
	  return 1;
	}
      } else {
	putc_unlocked('\n', outfile);
      }
    }
  }
  *sample_uidx_ptr = sample_uidx;
  return 0;
}

uint32_t write_haploview_map(FILE* outfile, uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos) {
  char wbuf[16];
  char* wptr;
  uintptr_t marker_idx;
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx_start++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx_start);
    fputs(&(marker_ids[marker_uidx_start * max_marker_id_len]), outfile);
    putc_unlocked('\t', outfile);
    wptr = uint32toa_x(marker_pos[marker_uidx_start], '\n', wbuf);
    fwrite(wbuf, 1, wptr - wbuf, outfile);
  }
  if (ferror(outfile)) {
    return 1;
  }
  return 0;
}

uint32_t valid_vcf_allele_code(const char* allele_code) {
  // returns 1 if probably valid (angle-bracket case is not exhaustively
  // checked), 0 if definitely not
  uint32_t uii = (unsigned char)(*allele_code);
  if ((uii == '<') || ((uii == '*') && (!allele_code[1]))) {
    return 1;
  }
  do {
    uii -= 64;
    // A = 1, C = 3, G = 7, N = 14, T = 20, so (0x10408a >> ucc) & 1 works as a
    // set membership test
#ifdef __LP64__
    if ((uii > 63) || (!((0x10408a0010408aLLU >> uii) & 1))) {
      // if '[', ']', or '.', assume breakend
      return ((uii == 27) || (uii == 29) || (uii == 0xffffffeeU))? 1 : 0;
    }
#else
    if ((uii > 63) || (!((0x10408a >> (uii % 32)) & 1))) {
      return ((uii == 27) || (uii == 29) || (uii == 0xffffffeeU))? 1 : 0;
    }
#endif
    uii = (unsigned char)(*(++allele_code));
  } while (uii);
  return 1;
}

int32_t flexbwrite_checked(const void* buf, size_t len, uint32_t output_bgz, FILE* outfile, BGZF* bgz_outfile) {
  if (!output_bgz) {
    return fwrite_checked(buf, len, outfile);
  } else {
    return (bgzf_write(bgz_outfile, buf, len) < 0);
  }
}

int32_t flexbputs_checked(const char* buf, uint32_t output_bgz, FILE* outfile, BGZF* bgz_outfile) {
  if (!output_bgz) {
    return fputs_checked(buf, outfile);
  } else {
    return (bgzf_write(bgz_outfile, buf, strlen(buf)) < 0);
  }
}

int32_t flexbputc_checked(unsigned char ucc, uint32_t output_bgz, FILE* outfile, BGZF* bgz_outfile) {
  if (!output_bgz) {
    return putc_checked(ucc, outfile);
  } else {
    return (bgzf_write(bgz_outfile, &ucc, 1) < 0);
  }
}

int32_t recode(uint32_t recode_modifier, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* recode_allele_name, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uint32_t* marker_pos, uintptr_t* marker_reverse, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint64_t misc_flags, uint32_t hh_exists, Chrom_info* chrom_info_ptr) {
  FILE* outfile = nullptr;
  FILE* outfile2 = nullptr;
  BGZF* bgz_outfile = nullptr;
  char* pzwritep = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  uintptr_t final_mask = get_final_mask(sample_ct);
  uintptr_t cur_final_mask = 0;
  uintptr_t sample_ct_y = 0;
  uintptr_t cur_sample_ct = 0;
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  char delimiter = (recode_modifier & RECODE_TAB)? '\t' : ' ';
  uintptr_t* recode_allele_reverse = nullptr;
  uintptr_t* sample_exclude_y = nullptr;
  uintptr_t* cur_sample_exclude = nullptr;
  char** mk_allele_ptrs = marker_allele_ptrs;
  char** allele_missing = nullptr;
  char* recode_allele_extra = nullptr;
  unsigned char* overflow_buf = nullptr;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char delim2 = delimiter;
  uintptr_t* sample_include2 = nullptr;
  uintptr_t* sample_include2_y = nullptr;
  uintptr_t* cur_sample_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uintptr_t* sample_male_include2_y = nullptr;
  uintptr_t* cur_sample_male_include2 = nullptr;
  uint32_t lgen_ref = (recode_modifier & RECODE_LGEN_REF);
  uint32_t rlist = (recode_modifier & RECODE_RLIST);
  uint32_t beagle_nomap = (recode_modifier & RECODE_BEAGLE_NOMAP);
  uint32_t vcf_not_fid = (recode_modifier & RECODE_VCF) && (!(recode_modifier & RECODE_FID));
  uint32_t vcf_not_iid = (recode_modifier & RECODE_VCF) && (!(recode_modifier & RECODE_IID));
  uint32_t vcf_two_ids = vcf_not_fid && vcf_not_iid;
  uint32_t output_bgz = (recode_modifier / RECODE_BGZ) & 1;
  uint32_t output_gen_gz = (recode_modifier / RECODE_GEN_GZ) & 1;
  uint32_t recode_012 = recode_modifier & (RECODE_01 | RECODE_12);
  uint32_t set_hh_missing = (misc_flags / MISC_SET_HH_MISSING) & 1;
  uint32_t set_mixed_mt_missing = (misc_flags / MISC_SET_MIXED_MT_MISSING) & 1;
  uint32_t real_ref_alleles = (misc_flags / MISC_REAL_REF_ALLELES) & 1;
  uint32_t xmhh_exists_orig = hh_exists & XMHH_EXISTS;
  uint32_t omit_nonmale_y = recode_modifier & RECODE_OMIT_NONMALE_Y;
  uintptr_t header_len = 0;
  uintptr_t max_chrom_size = 0;
  uintptr_t max_fid_len = 0;
  uintptr_t fid_ct = 0;
  uint32_t last_chrom_fo_idx = 0;
  uint32_t onechar_max = (chrom_info_ptr->max_code > 9)? 9 : chrom_info_ptr->max_code;
  uint32_t invalid_allele_code_seen = 0;
  uint32_t last_pos = 0;
  char missing_geno = *g_missing_geno_ptr;
  char output_missing_geno = *g_output_missing_geno_ptr;
  uintptr_t* loadbuf_collapsed = nullptr;
  uintptr_t* loadbuf_collapsed_end = nullptr;
  char* sample_ids_collapsed = nullptr;
  char* sample_ids_collapsed_y = nullptr;
  char* cur_sample_ids_collapsed = nullptr;
  char* writebuf = nullptr;
  char* writebuf2 = nullptr;
  char* writebuf3 = nullptr;
  uint32_t* fid_map = nullptr;
  uint32_t* missing_cts = nullptr;
  char* cur_mk_allelesx_buf = nullptr;
  int32_t retval = 0;
  char* writebufl[4];
  char* writebuflp[4];
  char* writebuflps[4];
  char* cur_mk_allelesx[6];
  char cur_dosage_chars[4];
  uint32_t cmalen[4];
  Pigz_state ps;
  time_t rawtime;
  struct tm *loctime;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uintptr_t autosomal_marker_ct;
  uintptr_t marker_uidx_start;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  unsigned char* loadbuf;
  uintptr_t* ulptr;
  uintptr_t* ulptr2;
  uintptr_t* ulptr_end;
  unsigned char* bufptr;
  char* wbufptr;
  char* cptr;
  char* aptr;
  char* aptr2;
  double dxx;
  uint32_t alen;
  uint32_t alen2;
  unsigned char ucc;
  unsigned char ucc2;
  char cc;
  uint32_t pct;
  uint32_t uidx_cur;
  uint32_t uidx_stop;
  uintptr_t loop_end;
  uintptr_t cur_word;
  uintptr_t ref_word;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t shiftval;
  uint32_t shiftmax;
  uint32_t cur_fid;
  uint32_t uii;
  int32_t ii;
  pzwrite_init_null(&ps);
  if (!hh_exists) {
    set_hh_missing = 0;
  }
  if (set_hh_missing || set_mixed_mt_missing || (recode_modifier & RECODE_VCF)) {
    if (recode_modifier & (RECODE_23 | RECODE_A_TRANSPOSE | RECODE_BEAGLE | RECODE_BEAGLE_NOMAP | RECODE_BIMBAM | RECODE_BIMBAM_1CHR | RECODE_LGEN | RECODE_LGEN_REF | RECODE_LIST | RECODE_OXFORD | RECODE_RLIST | RECODE_TRANSPOSE | RECODE_VCF)) {
      // SNP-major and no need for sample_uidx in inner loop, so we can use
      // collapsed representation
      if (alloc_collapsed_haploid_filters(sample_exclude, sex_male, unfiltered_sample_ct, sample_ct, hh_exists | ((recode_modifier & RECODE_VCF)? XMHH_EXISTS : 0) | (set_mixed_mt_missing? NXMHH_EXISTS : 0), 0, &sample_include2, &sample_male_include2)) {
	goto recode_ret_NOMEM;
      }
    } else {
      // sample-major output (in which case we load large blocks and use
      // haploid_fix_multiple())
      if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists | (set_mixed_mt_missing? NXMHH_EXISTS : 0), 0, sample_exclude, sex_male, &sample_include2, &sample_male_include2)) {
	goto recode_ret_NOMEM;
      }
    }
  }
  if (recode_012) {
    // may as well prevent user from shooting themselves in the foot here
    if ((recode_modifier & RECODE_01) && (output_missing_geno == '0') && (!(recode_modifier & (RECODE_A | RECODE_A_TRANSPOSE | RECODE_AD | RECODE_BIMBAM | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_STRUCTURE)))) {
      logerrprint("Error: The --recode '01' modifier normally has to be used with a nonzero\n--output-missing-genotype setting.\n");
      goto recode_ret_INVALID_CMDLINE;
    }
    mk_allele_ptrs = (char**)bigstack_alloc(unfiltered_marker_ct * 2 * sizeof(intptr_t));
    if (!mk_allele_ptrs) {
      goto recode_ret_NOMEM;
    }
    uii = 96 + (recode_modifier & RECODE_12);
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      mk_allele_ptrs[2 * marker_uidx] = (char*)(&(g_one_char_strs[uii])); // "0" or "1"
      mk_allele_ptrs[2 * marker_uidx + 1] = (char*)(&(g_one_char_strs[uii + 2])); // "1" or "2"
    }
    max_marker_allele_len = 2;
  }
  if (recode_modifier & (RECODE_A_TRANSPOSE | RECODE_BEAGLE | RECODE_BEAGLE_NOMAP | RECODE_BIMBAM | RECODE_BIMBAM_1CHR | RECODE_LGEN | RECODE_LGEN_REF | RECODE_LIST | RECODE_OXFORD | RECODE_RLIST | RECODE_TRANSPOSE | RECODE_VCF)) {
    if (bigstack_alloc_ul(sample_ctv2, &loadbuf_collapsed)) {
      goto recode_ret_NOMEM;
    }
    loadbuf_collapsed_end = &(loadbuf_collapsed[QUATERCT_TO_WORDCT(sample_ct)]);
    if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF | RECODE_LIST | RECODE_RLIST)) {
      // need to collapse sample_ids to remove need for sample_uidx in inner
      // loop
      sample_ids_collapsed = alloc_and_init_collapsed_arr(sample_ids, max_sample_id_len, unfiltered_sample_ct, sample_exclude, sample_ct, (delimiter == '\t'));
      if (!sample_ids_collapsed) {
        goto recode_ret_NOMEM;
      }
      if (omit_nonmale_y) {
	if (bigstack_alloc_ul(unfiltered_sample_ctl, &sample_exclude_y)) {
	  goto recode_ret_NOMEM;
	}
	memcpy(sample_exclude_y, sample_exclude, unfiltered_sample_ctl * sizeof(intptr_t));
	bitvec_ornot(sex_male, unfiltered_sample_ctl, sample_exclude_y);
	zero_trailing_bits(unfiltered_sample_ct, sample_exclude_y);
	sample_ct_y = unfiltered_sample_ct - popcount_longs(sample_exclude_y, unfiltered_sample_ctl);
        uii = QUATERCT_TO_ALIGNED_WORDCT(sample_ct_y);
	if (bigstack_alloc_ul(uii, &sample_include2_y) ||
            bigstack_alloc_ul(uii, &sample_male_include2_y)) {
	  goto recode_ret_NOMEM;
	}
	fill_quatervec_55(sample_ct_y, sample_include2_y);
	fill_quatervec_55(sample_ct_y, sample_male_include2_y);
	sample_ids_collapsed_y = alloc_and_init_collapsed_arr(sample_ids, max_sample_id_len, unfiltered_sample_ct, sample_exclude_y, sample_ct_y, (delimiter == '\t'));
      }
    }
  }
  if (recode_modifier & RECODE_VCF) {
    if (bigstack_alloc_c(sample_ct * 4, &writebuf)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & RECODE_OXFORD) {
    if (bigstack_alloc_uc(PIGZ_BLOCK_SIZE + sample_ct * 6 + 2 * max_marker_allele_len + MAXLINELEN, &overflow_buf) ||
        bigstack_calloc_ui(sample_ct, &missing_cts)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_BEAGLE | RECODE_BEAGLE_NOMAP)) {
    // common header:
    // "P FID " + ... + "\n"
    // "I IID " + ... + "\n"
    // "? PHE " + ... + "\n"
    // ... space requirement is bounded above by
    //   (2 * sample_ct * max_sample_id_len) + 64 * sample_ct
    //
    // per-marker:
    //   "M " + [marker name] + " " + ... + "\n"
    ulii = strlen(output_missing_pheno);
    if (bigstack_alloc_c(2 * ulii + 2, &writebuf) ||
        bigstack_alloc_c(21 + sample_ct * (2 * max_sample_id_len + 64), &writebuf2)) {
      goto recode_ret_NOMEM;
    }
    wbufptr = memcpyax(writebuf, output_missing_pheno, ulii, ' ');
    wbufptr = memcpyax(wbufptr, output_missing_pheno, ulii, ' ');
    ulii = 2 * ulii + 2;
    wbufptr = memcpya(writebuf2, "P FID ", 6);
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_sample_id_len);
      uljj = (uintptr_t)(aptr - cptr);
      wbufptr = memcpyax(wbufptr, cptr, uljj, ' ');
      wbufptr = memcpyax(wbufptr, cptr, uljj, ' ');
    }
    wbufptr = memcpya(wbufptr, "\nI IID ", 7);
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      cptr = (char*)memchr(cptr, '\t', max_sample_id_len);
      cptr++;
      uljj = strlen(cptr);
      wbufptr = memcpyax(wbufptr, cptr, uljj, ' ');
      wbufptr = memcpyax(wbufptr, cptr, uljj, ' ');
    }
    sample_uidx = 0;
    if (pheno_c) {
      wbufptr = memcpya(wbufptr, "\nA PHE ", 7);
      for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
        if (IS_SET(pheno_nm, sample_uidx)) {
	  cc = (unsigned char)('1' + IS_SET(pheno_c, sample_uidx));
	  *wbufptr++ = cc;
	  *wbufptr++ = ' ';
	  *wbufptr++ = cc;
	  *wbufptr++ = ' ';
	} else {
	  wbufptr = memcpya(wbufptr, writebuf, ulii);
	}
      }
    } else {
      wbufptr = memcpya(wbufptr, "\nT PHE ", 7);
      for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
        if (IS_SET(pheno_nm, sample_uidx)) {
	  cptr = dtoa_gx(pheno_d[sample_uidx], ' ', wbufptr);
          wbufptr = memcpya(cptr, wbufptr, (uintptr_t)(cptr - wbufptr));
	} else {
	  wbufptr = memcpya(wbufptr, writebuf, ulii);
	}
      }
    }
    *wbufptr++ = '\n';
    // free unused space, and save header length
    header_len = (uintptr_t)(wbufptr - writebuf2);
    bigstack_shrink_top(writebuf2, header_len);
    cmalen[1] = 4;
    ulii = 2 * max_marker_allele_len;
    if (bigstack_alloc_c(4 * max_marker_allele_len, &cur_mk_allelesx_buf)) {
      goto recode_ret_NOMEM;
    }
    cur_mk_allelesx[0] = cur_mk_allelesx_buf;
    cur_mk_allelesx[1] = &(cur_mk_allelesx_buf[ulii]);
    cur_mk_allelesx[2] = &(cur_mk_allelesx_buf[ulii * 2]);
    cur_mk_allelesx[3] = &(cur_mk_allelesx_buf[ulii * 3]);
    memcpy(cur_mk_allelesx[0], "    ", 4);
    memcpy(cur_mk_allelesx[1], " 0 0", 4);
    memcpy(cur_mk_allelesx[2], "    ", 4);
    memcpy(cur_mk_allelesx[3], "    ", 4);
  } else if (recode_modifier & (RECODE_BIMBAM | RECODE_BIMBAM_1CHR)) {
    if (max_marker_allele_len != 2) {
      logerrprint("Error: --recode bimbam cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = (char*)memchr(&(sample_ids[sample_uidx * max_sample_id_len]), '\t', max_sample_id_len);
      if (strchr(&(cptr[1]), ',')) {
        logerrprint("Error: Comma present in sample ID during --recode bimbam run.\n");
        goto recode_ret_INVALID_FORMAT;
      }
    }
    marker_uidx = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      if (strchr(&(marker_ids[marker_uidx * max_marker_id_len]), ',')) {
        logerrprint("Error: Comma present in SNP ID during --recode bimbam run.\n");
      }
    }
    // +1 because memcpyl3a() copies an extra character
    if (bigstack_alloc_c(3 * sample_ct + 1, &writebuf) ||
        bigstack_alloc_c(32, &writebuf2)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR)) {
    max_chrom_size = get_max_chrom_size(chrom_info_ptr, marker_exclude, &last_chrom_fo_idx);
    if ((recode_modifier & RECODE_FASTPHASE_1CHR) && (max_chrom_size != marker_ct)) {
      logerrprint("Error: --recode fastphase-1chr requires a single-chromosome dataset.  Did you\nmean '--recode fastphase'?  (Note the lack of a dash in the middle.)\n");
      goto recode_ret_INVALID_CMDLINE;
    }
    if (max_marker_allele_len != 2) {
      logerrprint("Error: --recode fastphase cannot be used with multi-character allele names.\n(You can use the '01' or '12' modifier to work around this.)\n");
      goto recode_ret_INVALID_CMDLINE;
    }
    if (recode_012) {
      if (bigstack_alloc_c(8, &writebuf3)) {
        goto recode_ret_NOMEM;
      }
      if (recode_modifier & RECODE_01) {
	memcpy(writebuf3, "0?010?11", 8);
      } else {
	memcpy(writebuf3, "1?121?22", 8);
      }
    } else {
      if (bigstack_alloc_c(max_chrom_size * 2, &writebuf3)) {
	goto recode_ret_NOMEM;
      }
    }
    if (bigstack_left() < ((uint64_t)unfiltered_sample_ct4) * max_chrom_size + 2 * round_up_pow2(max_chrom_size, CACHELINE)) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (bigstack_alloc_c(max_chrom_size, &writebuf) ||
        bigstack_alloc_c(max_chrom_size, &writebuf2)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF)) {
    ulii = 1 + 2 * max_marker_allele_len + max_marker_id_len + max_sample_id_len;
    if (bigstack_alloc_c(4 * ulii, &writebuf)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_LIST | RECODE_RLIST)) {
    // --list:
    // 3 for chromosome and delim
    // + max_marker_id_len
    // + 3, or (2 * max_marker_allele_len - 1)
    // + sample_ct * max_sample_id_len + 1
    //
    // --rlist:
    // max_marker_id_len
    // + 4 for "HOM"/"HET"/"NIL" and delim
    // + 4, or (2 * max_marker_allele_len)
    // + sample_ct * max_sample_id_len + 1
    ulii = 3 + max_marker_id_len + 2 * max_marker_allele_len + sample_ct * max_sample_id_len;
    if ((recode_modifier & RECODE_LIST) && (max_marker_allele_len != 2)) {
      logerrprint("Error: --recode list cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    if (rlist) {
      ulii += 2;
    }
    if (bigstack_alloc_c(ulii * 4, &writebuf)) {
      goto recode_ret_NOMEM;
    }
    writebufl[0] = writebuf;
    writebufl[1] = &(writebuf[ulii]);
    writebufl[2] = &(writebuf[ulii * 2]);
    writebufl[3] = &(writebuf[ulii * 3]);
  } else if (recode_modifier & RECODE_23) {
    if (sample_ct != 1) {
      logerrprint("Error: --recode 23 can only be used on a file with exactly one sample.\n");
      goto recode_ret_INVALID_FORMAT;
    } else if (max_marker_allele_len != 2) {
      logerrprint("Error: --recode 23 cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    // chromosome code, marker position, single-char alleles
    if (bigstack_alloc_c(32, &writebuf)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & RECODE_STRUCTURE) {
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_sample_id_len);
      ulii = (uintptr_t)(aptr - cptr);
      if (ulii >= max_fid_len) {
        max_fid_len = ulii + 1;
	if (max_fid_len == max_sample_id_len - 2) {
	  break;
	}
      }
    }
    if (bigstack_alloc_c(max_fid_len * sample_ct, &writebuf3)) {
      goto recode_ret_NOMEM;
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_fid_len);
      ulii = (uintptr_t)(aptr - cptr);
      memcpy(&(writebuf3[sample_idx * max_fid_len]), cptr, ulii);
      writebuf3[sample_idx * max_fid_len + ulii] = '\0';
    }
    qsort(writebuf3, sample_ct, max_fid_len, strcmp_casted);
    // simpler collapse_duplicate_ids(), probably want to make this a function
    for (ulii = 1; ulii < sample_ct; ulii++) {
      if (!strcmp(&(writebuf3[(ulii - 1) * max_fid_len]), &(writebuf3[ulii * max_fid_len]))) {
	break;
      }
    }
    fid_ct = MINV(ulii, sample_ct);
    while (++ulii < sample_ct) {
      if (strcmp(&(writebuf3[(fid_ct - 1) * max_fid_len]), &(writebuf3[ulii * max_fid_len]))) {
        strcpy(&(writebuf3[fid_ct * max_fid_len]), &(writebuf3[ulii * max_fid_len]));
	fid_ct++;
      }
    }
    bigstack_shrink_top(writebuf3, fid_ct * max_fid_len);
    if (bigstack_calloc_ui(fid_ct, &fid_map) ||
        bigstack_alloc_c(4 * marker_ct, &writebuf) ||
        bigstack_alloc_c(16, &writebuf2)) {
      goto recode_ret_NOMEM;
    }
    fill_uint_zero(fid_ct, fid_map);
  } else {
    if (recode_modifier & RECODE_A_TRANSPOSE) {
      // format is new to PLINK 1.9, so use tab delimiter unless 'spacex'
      // modifier present
      delimiter = ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_DELIMX)? ' ' : '\t';
      if (bigstack_alloc_c(sample_ct * 3 + 1, &writebuf)) {
        goto recode_ret_NOMEM;
      }
    } else {
      if (recode_modifier & RECODE_AD) {
	if (bigstack_alloc_c(32, &writebuf2)) {
	  goto recode_ret_NOMEM;
	}
	memcpy(writebuf2, "2 0     1 1 0 0 NA NA", 21);
	writebuf2[1] = delimiter;
	writebuf2[3] = delimiter;
	writebuf2[9] = delimiter;
	writebuf2[11] = delimiter;
	writebuf2[13] = delimiter;
	writebuf2[15] = delimiter;
	writebuf2[18] = delimiter;
	writebuf2[21] = delimiter;
      }
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	delim2 = ' ';
      }
      if (!(recode_modifier & RECODE_TRANSPOSE)) {
	max_chrom_size = marker_ct;
	if (recode_modifier & RECODE_AD) {
	  ulii = 6;
	} else if (recode_modifier & (RECODE_A | RECODE_COMPOUND)) {
	  if ((max_marker_allele_len != 2) && (recode_modifier & RECODE_COMPOUND)) {
	    logerrprint("Error: --recode compound-genotypes cannot be used with multi-character allele\nnames.\n");
	    goto recode_ret_INVALID_FORMAT;
	  }
	  ulii = 3;
	} else if (recode_modifier & (RECODE_HV | RECODE_HV_1CHR)) {
	  max_chrom_size = get_max_chrom_size(chrom_info_ptr, marker_exclude, &last_chrom_fo_idx);
	  if ((recode_modifier & RECODE_HV_1CHR) && (max_chrom_size != marker_ct)) {
	    logerrprint("Error: --recode HV-1chr requires a single-chromosome dataset.  Did you mean\n'--recode HV'?  (Note the lack of a dash in the middle.)\n");
	    goto recode_ret_INVALID_CMDLINE;
	  }
	  if (max_marker_allele_len == 2) {
	    ulii = max_chrom_size * 4;
	  } else {
	    ulii = 0;
	    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
	      uljj = 0;
	      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
	      for (marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end); marker_uidx < chrom_end;) {
		alen = strlen(mk_allele_ptrs[marker_uidx * 2]);
		alen2 = strlen(mk_allele_ptrs[marker_uidx * 2 + 1]);
		uljj += MAXV(alen, alen2) + 1;
		marker_uidx++;
		next_unset_ul_ck(marker_exclude, chrom_end, &marker_uidx);
	      }
	      if (uljj > ulii) {
		ulii = uljj;
	      }
	    }
	    ulii *= 2;
	  }
	} else {
	  // all chromosomes at once
	  // calculate maximum length of .ped line
	  if (max_marker_allele_len == 2) {
	    ulii = marker_ct * 4;
	  } else {
	    ulii = marker_ct;
	    for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	      alen = strlen(mk_allele_ptrs[marker_uidx * 2]);
	      alen2 = strlen(mk_allele_ptrs[marker_uidx * 2 + 1]);
	      ulii += MAXV(alen, alen2);
	    }
	    ulii *= 2;
	  }
	}
	if (recode_modifier & (RECODE_A | RECODE_AD | RECODE_COMPOUND)) {
	  if (bigstack_alloc_c(max_chrom_size * ulii, &writebuf)) {
	    goto recode_ret_NOMEM;
	  }
	  if ((recode_modifier & RECODE_COMPOUND) && max_chrom_size) {
	    memset(writebuf, delimiter, max_chrom_size * 3 - 1);
	    writebuf[max_chrom_size * 3 - 1] = '\n';
	  }
	} else {
	  // --recode, --recode HV
	  if (bigstack_alloc_c(ulii, &writebuf)) {
	    goto recode_ret_NOMEM;
	  }
	}
      }
    }
    if (recode_allele_name) {
      if (bigstack_calloc_ul(unfiltered_marker_ctl, &recode_allele_reverse)) {
	goto recode_ret_NOMEM;
      }
      // this indicates when we want to report the A2 allele instead of the
      // A1.  (potential double negatives, bleah)
      allele_missing = (char**)bigstack_alloc(unfiltered_marker_ct * sizeof(char**));
      if (!allele_missing) {
	goto recode_ret_NOMEM;
      }
      recode_allele_extra = (char*)g_bigstack_base;
      fill_ulong_zero(unfiltered_marker_ct, (uintptr_t*)allele_missing);
      ulii = round_up_pow2(max_marker_allele_len + MAXLINELEN, END_ALLOC_CHUNK);
      loadbuf = (unsigned char*)bigstack_end_alloc_presized(ulii);
      if (!loadbuf) {
	goto recode_ret_NOMEM;
      }
      // When '12' and 'A'/'AD' are simultaneously present, most sensible
      // behavior is to match against real allele IDs and just apply '12'
      // to the output header line.  If that's not what the user wants,
      // they can do a two-step recode.
      // (--recode12 simply overrode --recodeA/--recodeAD in PLINK 1.07; no
      // need to replicate that.)
      retval = recode_allele_load((char*)loadbuf, ulii, recode_allele_name, &allele_missing, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, max_marker_allele_len, recode_allele_reverse, recode_allele_extra);
      bigstack_end_reset(bigstack_end_mark);
      if (retval) {
	goto recode_ret_1;
      }
    }
  }

  if (!(recode_modifier & (RECODE_A | RECODE_AD | RECODE_BEAGLE | RECODE_BEAGLE_NOMAP | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_LGEN | RECODE_LGEN_REF | RECODE_OXFORD | RECODE_VCF))) {
    if (bigstack_alloc_c(8 * max_marker_allele_len, &cur_mk_allelesx_buf)) {
      goto recode_ret_NOMEM;
    }
    cur_mk_allelesx[0] = cur_mk_allelesx_buf;
    cur_mk_allelesx[1] = &(cur_mk_allelesx_buf[max_marker_allele_len * 2]);
    cur_mk_allelesx[2] = &(cur_mk_allelesx_buf[max_marker_allele_len * 4]);
    cur_mk_allelesx[3] = &(cur_mk_allelesx_buf[max_marker_allele_len * 6]);
  } else if (recode_modifier & RECODE_VCF) {
    if (bigstack_alloc_c(16, &cur_mk_allelesx_buf)) {
      goto recode_ret_NOMEM;
    }
    memcpy(cur_mk_allelesx_buf, "\t1/1\t./.\t0/1\t0/0", 16);
  } else if (recode_modifier & RECODE_OXFORD) {
    if (bigstack_alloc_c(32, &cur_mk_allelesx_buf)) {
      goto recode_ret_NOMEM;
    }
    memcpy(cur_mk_allelesx_buf, " 1 0 0   0 0 0   0 1 0   0 0 1", 30);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto recode_ret_READ_FAIL;
  }
  marker_uidx = 0;
  marker_idx = 0;
  if (bigstack_left() < unfiltered_sample_ct4) {
    goto recode_ret_NOMEM;
  }
  loadbuf = g_bigstack_base;
  chrom_fo_idx = 0;
  if (unfiltered_marker_ct) {
    refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
  } else {
    chrom_end = 0;
    is_x = 0;
    is_y = 0;
    is_mt = 0;
    is_haploid = 0;
  }
  if (recode_modifier & RECODE_TRANSPOSE) {
    strcpy(outname_end, ".tped");
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    *outname_end = '\0';
    LOGPRINTFWW5("--recode transpose to %s.tped + %s.tfam ... ", outname, outname);
    fputs("0%", stdout);
    fflush(stdout);
    cur_mk_allelesx[1][0] = delimiter;
    cur_mk_allelesx[1][1] = output_missing_geno;
    cur_mk_allelesx[1][2] = delim2;
    cur_mk_allelesx[1][3] = output_missing_geno;
    cmalen[1] = 4;

    cptr = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
    *cptr++ = delimiter;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  cptr = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
	  *cptr++ = delimiter;
	}
	wbufptr = strcpyax(cptr, &(marker_ids[marker_uidx * max_marker_id_len]), delimiter);
	if (!marker_cms) {
	  *wbufptr++ = '0';
	} else {
	  wbufptr = dtoa_g(marker_cms[marker_uidx], wbufptr);
	}
	*wbufptr++ = delimiter;
	wbufptr = uint32toa(marker_pos[marker_uidx], wbufptr);
        if (fwrite_checked(g_textbuf, wbufptr - g_textbuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}

	if (sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, sample_include2, sample_ct);
	  }
	  init_recode_cmax0(mk_allele_ptrs[2 * marker_uidx], mk_allele_ptrs[2 * marker_uidx + 1], cur_mk_allelesx, cmalen, delimiter, delim2);
	  ulptr = loadbuf_collapsed;
	  ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
	  shiftmax = BITCT2;
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		fwrite(cur_mk_allelesx[ulii], 1, cmalen[ulii], outfile);
		cur_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = sample_ct % BITCT2;
	  }
	}
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & RECODE_A_TRANSPOSE) {
    strcpy(outname_end, ".traw");
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    fputs((delimiter == '\t')? "CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT" : "CHR SNP (C)M POS COUNTED ALT", outfile);
    shiftval = 0; // repurposed: underscore seen in ID?
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      if (!shiftval) {
	if (strchr(cptr, '_')) {
	  shiftval = 1;
	  logerrprint("Warning: Underscore(s) present in sample IDs.\n");
	}
      }
      aptr = (char*)memchr(cptr, '\t', max_sample_id_len);
      putc_unlocked(delimiter, outfile);
      fwrite(cptr, 1, (uintptr_t)(aptr - cptr), outfile);
      putc_unlocked('_', outfile);
      fputs(&(aptr[1]), outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    LOGPRINTFWW5("--recode A-transpose to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    cptr = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
    *cptr++ = delimiter;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  cptr = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
	  *cptr++ = delimiter;
	}
	wbufptr = strcpyax(cptr, &(marker_ids[marker_uidx * max_marker_id_len]), delimiter);
	if (!marker_cms) {
	  *wbufptr++ = '0';
	} else {
	  wbufptr = dtoa_g(marker_cms[marker_uidx], wbufptr);
	}
	*wbufptr++ = delimiter;
	wbufptr = uint32toa_x(marker_pos[marker_uidx], delimiter, wbufptr);
        if (fwrite_checked(g_textbuf, wbufptr - g_textbuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	uii = IS_NONNULL_AND_SET(recode_allele_reverse, marker_uidx);
	if (allele_missing && allele_missing[marker_uidx]) {
	  fputs(allele_missing[marker_uidx], outfile);
	  putc_unlocked(delimiter, outfile);
          fputs(mk_allele_ptrs[2 * marker_uidx + uii], outfile);
	  putc_unlocked(',', outfile);
	} else {
	  fputs(mk_allele_ptrs[2 * marker_uidx + uii], outfile);
	  putc_unlocked(delimiter, outfile);
	}
	fputs(mk_allele_ptrs[2 * marker_uidx + 1 - uii], outfile);
	wbufptr = writebuf;
	if (sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, uii ^ IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, sample_include2, sample_ct);
	  }
	  ulptr = loadbuf_collapsed;
	  ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
	  shiftmax = BITCT2;
	  if (allele_missing && allele_missing[marker_uidx]) {
	    // all 0s and NAs
	    memcpy(cur_dosage_chars, "0N00", 4);
	  } else {
	    memcpy(cur_dosage_chars, "2N10", 4);
	  }
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		*wbufptr++ = delimiter;
		*wbufptr++ = cur_dosage_chars[ulii];
		if (ulii == 1) {
		  *wbufptr++ = 'A';
		}
		cur_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = sample_ct % BITCT2;
	  }
	}
	*wbufptr++ = '\n';
	if (fwrite_checked(writebuf, wbufptr - writebuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & RECODE_VCF) {
    if (!output_bgz) {
      memcpy(outname_end, ".vcf", 5);
      if (fopen_checked(outname, "w", &outfile)) {
	goto recode_ret_OPEN_FAIL;
      }
    } else {
      memcpy(outname_end, ".vcf.gz", 8);
      bgz_outfile = bgzf_open(outname, "w");
      if (!bgz_outfile) {
	goto recode_ret_OPEN_FAIL;
      }
#ifndef _WIN32
      if (g_thread_ct > 1) {
	bgzf_mt(bgz_outfile, g_thread_ct, 128);
      }
#endif
    }
    wbufptr = memcpya(g_textbuf, "##fileformat=VCFv4.2\n##fileDate=", 32);
    time(&rawtime);
    loctime = localtime(&rawtime);
    wbufptr += strftime(wbufptr, MAXLINELEN, "%Y%m%d", loctime);
    wbufptr = memcpya(wbufptr, "\n##source=PLINKv1.90\n", 21);
    uii = 0; // '0' written already?
    if (flexbwrite_checked(g_textbuf, wbufptr - g_textbuf, output_bgz, outfile, bgz_outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    memcpy(g_textbuf, "##contig=<ID=", 13);
    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if (!IS_SET(chrom_info_ptr->chrom_mask, chrom_idx)) {
	continue;
      }
      cptr = chrom_name_write(chrom_info_ptr, chrom_idx, &(g_textbuf[13]));
      if ((g_textbuf[13] == '0') && (cptr == &(g_textbuf[14]))) {
	if (uii) {
	  continue;
	}
	uii = 1;
	cptr = memcpya(cptr, ",length=2147483645", 18);
      } else {
	*cptr = '\0';
	if (strchr(&(g_textbuf[13]), ':')) {
	  logerrprint("Error: VCF chromosome codes may not include the ':' character.\n");
	  goto recode_ret_INVALID_FORMAT;
	}
        cptr = memcpya(cptr, ",length=", 8);
	if (!(map_is_unsorted & UNSORTED_BP)) {
	  cptr = uint32toa(marker_pos[chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - 1] + 1, cptr);
	} else {
	  cptr = memcpya(cptr, "2147483645", 10); // unknown
	}
      }
      cptr = memcpya(cptr, ">\n", 2);
      if (flexbwrite_checked(g_textbuf, cptr - g_textbuf, output_bgz, outfile, bgz_outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
    if (!real_ref_alleles) {
      if (flexbputs_checked("##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">\n", output_bgz, outfile, bgz_outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
    // todo: include PEDIGREE in header, and make --vcf be able to read it?
    if (flexbputs_checked(
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", output_bgz, outfile, bgz_outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    chrom_fo_idx = 0;
    if (unfiltered_marker_ct) {
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    }
    shiftval = 0; // repurposed: underscore seen in ID?
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      ulii = strlen_se(cptr);
      if (flexbputc_checked('\t', output_bgz, outfile, bgz_outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
      if (vcf_not_iid) {
	if (flexbwrite_checked(cptr, ulii, output_bgz, outfile, bgz_outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	if (vcf_two_ids) {
	  if (!shiftval) {
	    if (strchr(cptr, '_')) {
	      shiftval = 1;
	      logerrprint("Warning: Underscore(s) present in sample IDs.\n");
	    }
	  }
	  if (flexbputc_checked('_', output_bgz, outfile, bgz_outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
      }
      if (vcf_not_fid) {
	if (flexbputs_checked(&(cptr[ulii + 1]), output_bgz, outfile, bgz_outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
    }
    LOGPRINTFWW5("--recode vcf%s%s to %s ... ", vcf_not_iid? (vcf_not_fid? "" : "-fid") : "-iid", output_bgz? " bgz" : "", outname);
    fputs("0%", stdout);
    fflush(stdout);
    g_textbuf[0] = '\n';
    if ((((!hh_exists) || set_hh_missing) && is_haploid && (!is_x)) || (set_mixed_mt_missing && is_mt)) {
      uii = 2;
    } else {
      uii = 4;
    }
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
          marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  if ((((!hh_exists) || set_hh_missing) && is_haploid && (!is_x)) || (set_mixed_mt_missing && is_mt)) {
	    uii = 2;
	  } else {
	    uii = 4;
	  }
	}
	wbufptr = chrom_name_write(chrom_info_ptr, chrom_idx, &(g_textbuf[1]));
	*wbufptr++ = '\t';
	wbufptr = uint32toa_x(marker_pos[marker_uidx], '\t', wbufptr);
	wbufptr = strcpyax(wbufptr, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
	if (flexbwrite_checked(g_textbuf, wbufptr - g_textbuf, output_bgz, outfile, bgz_outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	cptr = mk_allele_ptrs[2 * marker_uidx + 1];
	if (cptr == missing_geno_ptr) {
	  if (flexbputc_checked('N', output_bgz, outfile, bgz_outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	} else {
          if ((!invalid_allele_code_seen) && (!valid_vcf_allele_code(cptr))) {
            invalid_allele_code_seen = 1;
	  }
	  if (flexbputs_checked(cptr, output_bgz, outfile, bgz_outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
	if (flexbputc_checked('\t', output_bgz, outfile, bgz_outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}

	if (sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, sample_include2, sample_ct);
	  }
	}

	cptr = mk_allele_ptrs[2 * marker_uidx];
	if (cptr != missing_geno_ptr) {
          if ((!invalid_allele_code_seen) && (!valid_vcf_allele_code(cptr))) {
	    invalid_allele_code_seen = 1;
	  }
	  // if ALT allele was not actually present in immediate dataset, VCF
	  // spec actually used to require '.'
	  // this restriction was contrary to practice (even by the 1000G team)
	  // and was removed from plink after b3.41.
	  if (flexbputs_checked(cptr, output_bgz, outfile, bgz_outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	} else {
	  if (flexbputc_checked('.', output_bgz, outfile, bgz_outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
	if (flexbputs_checked(real_ref_alleles? "\t.\t.\t.\tGT" : "\t.\t.\tPR\tGT", output_bgz, outfile, bgz_outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	wbufptr = writebuf;
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
	shiftmax = BITCT2;
	if ((!is_x) || (hh_exists && (!set_hh_missing))) {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		wbufptr = memcpya(wbufptr, &(cur_mk_allelesx_buf[ulii * 4]), uii);
		cur_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = sample_ct % BITCT2;
	  }
	} else {
	  ulptr2 = sample_male_include2;
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      ref_word = (*ulptr2++) << 1;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		uljj = ref_word & 3;
		wbufptr = memcpya(wbufptr, &(cur_mk_allelesx_buf[ulii * 4]), 4 - uljj);
		cur_word >>= 2;
		ref_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = sample_ct % BITCT2;
	  }
	}
	if (flexbwrite_checked(writebuf, wbufptr - writebuf, output_bgz, outfile, bgz_outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (flexbputc_checked('\n', output_bgz, outfile, bgz_outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    if (output_bgz) {
      if (bgzf_close(bgz_outfile)) {
	bgz_outfile = nullptr;
	goto recode_ret_WRITE_FAIL;
      }
      bgz_outfile = nullptr;
    }
  } else if (recode_modifier & RECODE_OXFORD) {
    memcpy(outname_end, output_gen_gz? ".gen.gz" : ".gen", output_gen_gz? 8 : 5);
    if (flex_pzwrite_init(output_gen_gz, outname, overflow_buf, 0, &ps)) {
      goto recode_ret_OPEN_FAIL;
    }
    *outname_end = '\0';
    LOGPRINTFWW5("--recode oxford%s to %s.gen%s + %s.sample ... ", output_gen_gz? " gen-gz" : "", outname, output_gen_gz? ".gz" : "", outname);
    pzwritep = (char*)overflow_buf;
    fputs("0%", stdout);
    fflush(stdout);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
        if (IS_SET(marker_exclude, marker_uidx)) {
          marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
          if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
        if (marker_uidx >= chrom_end) {
          chrom_fo_idx++;
          refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  // not clear from documentation whether anything special should be
	  // done for Y/haploid chromosomes
	}
	pzwritep = chrom_name_write(chrom_info_ptr, chrom_idx, pzwritep);
	*pzwritep++ = ' ';
	pzwritep = strcpyax(pzwritep, &(marker_ids[marker_uidx * max_marker_id_len]), ' ');
	pzwritep = uint32toa_x(marker_pos[marker_uidx], ' ', pzwritep);
	pzwritep = strcpyax(pzwritep, mk_allele_ptrs[2 * marker_uidx], ' ');
	pzwritep = strcpya(pzwritep, mk_allele_ptrs[2 * marker_uidx + 1]);
	if (sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, sample_include2, sample_ct);
	  }
	  ulptr = loadbuf_collapsed;
	  ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
	  sample_idx = 0;
	  sample_uidx = BITCT2; // repurposed as stop value
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; sample_idx < sample_uidx; sample_idx++, cur_word >>= 2) {
		ulii = cur_word & 3;
		if (ulii == 1) {
		  missing_cts[sample_idx] += 1;
		}
		pzwritep = memcpya(pzwritep, &(cur_mk_allelesx_buf[ulii * 8]), 6);
	      }
	      sample_uidx += BITCT2;
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    sample_uidx = sample_ct;
	  }
	}
	append_binary_eoln(&pzwritep);
	if (flex_pzwrite(&ps, &pzwritep)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (flex_pzwrite_close_null(&ps, pzwritep)) {
      goto recode_ret_WRITE_FAIL;
    }
    memcpy(outname_end, ".sample", 8);
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    if (fputs_checked(
"ID_1 ID_2 missing sex phenotype\n"
"0 0 0 D "
, outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    putc_unlocked(pheno_d? 'P' : 'B', outfile);
    putc_unlocked('\n', outfile);
    dxx = 1.0 / ((double)((intptr_t)marker_ct));
    for (sample_idx = 0, sample_uidx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_sample_id_len);
      wbufptr = memcpyax(g_textbuf, cptr, aptr - cptr, ' ');
      wbufptr = strcpyax(wbufptr, &(aptr[1]), ' ');
      wbufptr = dtoa_gx(((double)((int32_t)missing_cts[sample_idx])) * dxx, ' ', wbufptr);
      *wbufptr++ = sexchar(sex_nm, sex_male, sample_uidx);
      *wbufptr++ = ' ';
      if (IS_SET(pheno_nm, sample_uidx)) {
        if (pheno_d) {
          wbufptr = dtoa_g(pheno_d[sample_uidx], wbufptr);
        } else {
          *wbufptr++ = '0' + IS_SET(pheno_c, sample_uidx);
	}
      } else {
	wbufptr = memcpya(wbufptr, "NA", 2);
      }
      *wbufptr++ = '\n';
      if (fwrite_checked(g_textbuf, wbufptr - g_textbuf, outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
  } else if (recode_modifier & RECODE_23) {
    memcpy(outname_end, ".txt", 5);
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    LOGPRINTFWW5("--recode 23 to %s ... ", outname);
    time(&rawtime);
    fprintf(outfile, "# This data file generated by " PROG_NAME_CAPS " at: %s", ctime(&rawtime));
    if (fputs_checked(
"#\n"
"# Below is a text version of your data.  Fields are TAB-separated.\n"
"# Each line corresponds to a single SNP.  For each SNP, we provide its\n"
"# identifier, its location on a reference human genome, and the genotype call.\n"
"# For further information (e.g. which reference build was used), consult the\n"
"# original source of your data.\n"
"#\n"
"# rsid\tchromosome\tposition\tgenotype\n"
, outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    writebuf[0] = '\t';
    cptr = chrom_print_human(chrom_idx, &(writebuf[1]));
    *cptr++ = '\t';
    sample_uidx = next_unset_unsafe(sample_exclude, 0);
    ucc = IS_SET(sex_male, sample_uidx);
    ucc2 = ((chrom_idx == 24) || (chrom_idx == 26) || (ucc && (chrom_idx == 23) && (!xmhh_exists_orig)))? 1 : 0;
    // todo: add a way to force two-allele output for entire male
    // pseudo-autosomal region
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        cptr = chrom_print_human(chrom_idx, &(writebuf[1]));
	*cptr++ = '\t';
	ucc2 = ((chrom_idx == 24) || (chrom_idx == 26) || (ucc && (chrom_idx == 23) && (!xmhh_exists_orig)))? 1 : 0;
      }

      if (fseeko(bedfile, bed_offset + (sample_uidx / 4) + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto recode_ret_READ_FAIL;
      }
      ii = getc_unlocked(bedfile);
      if (ii == EOF) {
	goto recode_ret_READ_FAIL;
      }
      cur_word = (((uint32_t)ii) >> ((sample_uidx % 4) * 2)) & 3;
      if (is_haploid && set_hh_missing) {
	// ok, this is a bit silly
	haploid_fix(hh_exists, sample_include2, sample_male_include2, 1, is_x, is_y, (unsigned char*)(&cur_word));
      } else if (is_mt && set_mixed_mt_missing) {
	hh_reset((unsigned char*)(&cur_word), sample_include2, 1);
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	cur_word = cur_word ^ (((~(cur_word ^ (cur_word >> 1))) & 1) * 3);
      }
      aptr = uint32toa_x(marker_pos[marker_uidx], '\t', cptr);
      if (cur_word) {
	if (cur_word == 3) {
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	} else if (cur_word == 2) {
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	} else {
	  aptr = memseta(aptr, '-', 2);
	}
      } else {
	*aptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	*aptr++ = mk_allele_ptrs[2 * marker_uidx][0];
      }
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      if (ucc2 && ((cur_word == 3) || (!cur_word))) {
	aptr--;
      }
      *aptr++ = '\n';
      if (fwrite_checked(writebuf, aptr - writebuf, outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
  } else if (recode_modifier & (RECODE_BEAGLE | RECODE_BEAGLE_NOMAP)) {
    // for backward compatibility, also exclude XY.  don't exclude custom name
    // chromosomes, though, since chromosome 0 was actually processed
    autosomal_marker_ct = marker_ct - count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
    if (chrom_info_ptr->xymt_codes[XY_OFFSET] != -2) {
      autosomal_marker_ct -= count_chrom_markers(chrom_info_ptr, marker_exclude, chrom_info_ptr->xymt_codes[XY_OFFSET]);
    }
    if (!autosomal_marker_ct) {
      // could allow this?
      logerrprint("Error: No autosomal variants for --recode beagle.\n");
      goto recode_ret_ALL_MARKERS_EXCLUDED;
    }
    if (!beagle_nomap) {
      memcpy(outname_end, ".chr-", 6);
      sprintf(g_logbuf, "--recode beagle to %s*.dat + %s*.map... ", outname, outname);
      wordwrapb(5);
      fputs(g_logbuf, stdout);
    } else {
      memcpy(outname_end, ".beagle.dat", 12);
      LOGPRINTFWW5("--recode beagle to %s ... ", outname);
    }
    fputs("0%", stdout);
    fflush(stdout);
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    retval = recode_beagle_new_chrom(outname, &(outname_end[5]), marker_exclude, chrom_info_ptr, &marker_uidx, &chrom_fo_idx, &chrom_idx, &chrom_end, bedfile, bed_offset, unfiltered_sample_ct4, &outfile, beagle_nomap? nullptr : (&outfile2), writebuf2, header_len);
    if (retval) {
      goto recode_ret_1;
    }
    if (beagle_nomap) {
      if (fopen_checked(outname, "w", &outfile)) {
	goto recode_ret_OPEN_FAIL;
      }
      if (fwrite_checked(writebuf2, header_len, outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
      outfile2 = nullptr;
    }
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * autosomal_marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
          if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
        if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  if (outfile2) {
	    if (fclose_null(&outfile) || fclose_null(&outfile2)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  retval = recode_beagle_new_chrom(outname, &(outname_end[5]), marker_exclude, chrom_info_ptr, &marker_uidx, &chrom_fo_idx, &chrom_idx, &chrom_end, bedfile, bed_offset, unfiltered_sample_ct4, &outfile, beagle_nomap? nullptr : (&outfile2), writebuf2, header_len);
	  if (retval) {
	    goto recode_ret_1;
	  }
	}
	if (sample_ct && load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	  goto recode_ret_READ_FAIL;
	}
	cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
        if (fputs_checked("M ", outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
        fputs(cptr, outfile);
        putc_unlocked(' ', outfile);
	if (outfile2) {
	  if (fputs_checked(cptr, outfile2)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  putc_unlocked('\t', outfile2);
	  wbufptr = uint32toa_x(marker_pos[marker_uidx], '\t', g_textbuf);
	  fwrite(g_textbuf, 1, wbufptr - g_textbuf, outfile2);
	}
	aptr = mk_allele_ptrs[2 * marker_uidx];
	aptr2 = mk_allele_ptrs[2 * marker_uidx + 1];
	alen = strlen(aptr);
	alen2 = strlen(aptr2);
	wbufptr = memcpyax(&(cur_mk_allelesx[0][1]), aptr, alen, ' ');
	memcpy(wbufptr, aptr, alen);
	cmalen[0] = 2 * alen + 2;
	wbufptr = memcpyax(&(cur_mk_allelesx[2][1]), aptr, alen, ' ');
	memcpy(wbufptr, aptr2, alen2);
	if (outfile2) {
	  fputs((aptr != missing_geno_ptr)? aptr : "X", outfile2);
	  putc_unlocked('\t', outfile2);
	  fputs((aptr2 != missing_geno_ptr)? aptr2 : "X", outfile2);
	}
	cmalen[2] = alen + alen2 + 2;
	wbufptr = memcpyax(&(cur_mk_allelesx[3][1]), aptr2, alen2, ' ');
	memcpy(wbufptr, aptr2, alen2);
	cmalen[3] = 2 * alen2 + 2;
	if (outfile2) {
	  if (putc_checked('\n', outfile2)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
        shiftmax = BITCT2;
        while (1) {
	  while (ulptr < ulptr_end) {
            cur_word = *ulptr++;
            for (shiftval = 0; shiftval < shiftmax; shiftval++) {
	      ulii = cur_word & 3;
	      fwrite(cur_mk_allelesx[ulii], 1, cmalen[ulii], outfile);
	      cur_word >>= 2;
	    }
	  }
	  if (ulptr == loadbuf_collapsed_end) {
	    break;
	  }
	  ulptr_end++;
	  shiftmax = sample_ct % BITCT2;
	}
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (outfile2) {
      if (fclose_null(&outfile2)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
  } else if (recode_modifier & (RECODE_BIMBAM | RECODE_BIMBAM_1CHR)) {
    if (recode_modifier & RECODE_BIMBAM_1CHR) {
      if (!marker_ct) {
	logerrprint("Error: No variants for --recode bimbam-1chr.\n");
	goto recode_ret_ALL_MARKERS_EXCLUDED;
      }
      ii = single_chrom_start(chrom_info_ptr, marker_exclude, unfiltered_marker_ct);
      if (ii == -1) {
        logerrprint("Error: --recode bimbam-1chr requires a single-chromosome dataset.  Did you mean\n'--recode bimbam'?  (Note the lack of a dash in the middle.)\n");
        goto recode_ret_INVALID_CMDLINE;
      }
      marker_uidx = ((uint32_t)ii);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto recode_ret_READ_FAIL;
      }
    }
    memcpy(outname_end, ".recode.", 9);
    LOGPRINTFWW5("--recode bimbam%s to %sgeno.txt + %spheno.txt + %spos.txt ... ", (recode_modifier & RECODE_BIMBAM_1CHR)? "-1chr" : "", outname, outname, outname);
    fputs("0%", stdout);
    fflush(stdout);
    memcpy(&(outname_end[8]), "pos.txt", 8);
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    writebuf2[0] = ' ';
    ulii = recode_modifier & RECODE_BIMBAM;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      wbufptr = uint32toa(marker_pos[marker_uidx], &(writebuf2[1]));
      if (ulii) {
	if (marker_uidx >= chrom_end) {
          chrom_idx = get_variant_chrom(chrom_info_ptr, marker_uidx);
          chrom_end = get_chrom_end_vidx(chrom_info_ptr, chrom_idx);
	}
        *wbufptr++ = ' ';
        wbufptr = chrom_name_write(chrom_info_ptr, chrom_idx, wbufptr);
      }
      *wbufptr++ = '\n';
      if (fwrite_checked(writebuf2, wbufptr - writebuf2, outfile)) {
        goto recode_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    memcpy(&(outname_end[8]), "pheno.txt", 10);
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    sample_uidx = 0;
    if (pheno_c) {
      for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	if (IS_SET(pheno_nm, sample_uidx)) {
          putc_unlocked('1' + IS_SET(pheno_c, sample_uidx), outfile);
	} else {
	  fputs(output_missing_pheno, outfile);
	}
        if (putc_checked('\n', outfile)) {
          goto recode_ret_WRITE_FAIL;
	}
      }
    } else {
      for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	if (IS_SET(pheno_nm, sample_uidx)) {
          wbufptr = dtoa_g(pheno_d[sample_uidx], writebuf2);
	  fwrite(writebuf2, 1, (uintptr_t)(wbufptr - writebuf2), outfile);
	} else {
          fputs(output_missing_pheno, outfile);
	}
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    memcpy(&(outname_end[8]), "geno.txt", 9);
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    wbufptr = uint32toa_x(sample_ct, '\n', writebuf2);
    wbufptr = uint32toa(marker_ct, wbufptr);
    wbufptr = memcpya(wbufptr, "\nIND", 4);
    if (fwrite_checked(writebuf2, wbufptr - writebuf2, outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      cptr = (char*)memchr(&(sample_ids[sample_uidx * max_sample_id_len]), '\t', max_sample_id_len);
      putc_unlocked(',', outfile);
      fputs(&(cptr[1]), outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    marker_idx = 0;
    chrom_fo_idx = 0;
    if (unfiltered_marker_ct) {
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    }
    writebuf2[0] = ',';
    memcpy(&(writebuf2[4]), ",??", 4);
    writebuf2[8] = ',';
    writebuf2[12] = ',';
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
          if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
            goto recode_ret_READ_FAIL;
	  }
        }
	if (marker_uidx >= chrom_end) {
          chrom_fo_idx++;
          refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
          chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	}
	if (fputs_checked(&(marker_ids[marker_uidx * max_marker_id_len]), outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	if (sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, sample_include2, sample_ct);
	  }
	  ucc = mk_allele_ptrs[2 * marker_uidx][0];
	  ucc2 = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  writebuf2[1] = ucc;
	  writebuf2[2] = ucc;
	  writebuf2[9] = ucc;
	  writebuf2[10] = ucc2;
	  writebuf2[13] = ucc2;
	  writebuf2[14] = ucc2;
	  wbufptr = writebuf;
	  ulptr = loadbuf_collapsed;
	  ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
	  shiftmax = BITCT2;
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		wbufptr = memcpyl3a(wbufptr, &(writebuf2[4 * ulii]));
		cur_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = sample_ct % BITCT2;
	  }
	  if (fwrite_checked(writebuf, 3 * sample_ct, outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
	putc_unlocked('\n', outfile);
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & (RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR)) {
    if (!marker_ct) {
      // why bother
      logerrprint("Error: No variants for --recode fastphase{-1chr}.\n");
      goto recode_ret_ALL_MARKERS_EXCLUDED;
    } else if (!sample_ct) {
      logerrprint("Error: No samples for --recode fastphase{-1chr}.\n");
      goto recode_ret_ALL_SAMPLES_EXCLUDED;
    }
    if (recode_modifier & RECODE_FASTPHASE) {
      memcpy(outname_end, ".chr-*", 7);
    } else {
      *outname_end = '\0';
    }
    sprintf(g_logbuf, "--recode fastphase%s to %s.recode.phase.inp ... ", (recode_modifier & RECODE_FASTPHASE)? "" : "-1chr", outname);
    wordwrapb(15); // strlen("[chromosome 10]")
    fputs(g_logbuf, stdout);
    chrom_fo_idx = 0xffffffffU; // exploit overflow for initialization
    if (recode_modifier & RECODE_FASTPHASE) {
      fputs("[chromosome   ", stdout);
    }
    writebuf3[1] = '?';
    writebuf3[5] = '?';
    do {
      chrom_fo_idx++;
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      ulii = count_chrom_markers(chrom_info_ptr, marker_exclude, chrom_idx);
      if (recode_modifier & RECODE_FASTPHASE) {
        wbufptr = chrom_name_write(chrom_info_ptr, chrom_idx, &(outname_end[5]));
        if (chrom_idx <= chrom_info_ptr->max_code) {
          printf("\b\b%u] \b", chrom_idx);
	} else {
	  // nonstandard chromosome name
          fputs("\b\b**] \b", stdout);
	}
        fflush(stdout);
      } else {
	wbufptr = outname_end;
      }
      memcpy(wbufptr, ".recode.phase.inp", 18);
      retval = open_and_write_fastphase_header(&outfile, outname, marker_exclude, marker_pos, marker_uidx, ulii, sample_ct);
      if (retval) {
	goto recode_ret_1;
      }
      marker_uidx_start = marker_uidx;
      if (recode_load_to(loadbuf, bedfile, bed_offset, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1], 0, ulii, marker_exclude, marker_reverse, &marker_uidx, unfiltered_sample_ct)) {
	goto recode_ret_READ_FAIL;
      }
      if (set_hh_missing || set_mixed_mt_missing) {
	haploid_fix_multiple(marker_exclude, marker_uidx_start, ulii, chrom_info_ptr, hh_exists, set_hh_missing, set_mixed_mt_missing, sample_include2, sample_male_include2, unfiltered_sample_ct, unfiltered_sample_ct4, loadbuf);
      }
      sample_uidx = 0;
      if (!recode_012) {
	uidx_cur = marker_uidx_start;
	cptr = writebuf3;
	for (marker_idx = 0; marker_idx < ulii; uidx_cur++, marker_idx++) {
          next_unset_unsafe_ck(marker_exclude, &uidx_cur);
	  *cptr++ = marker_allele_ptrs[2 * uidx_cur][0];
	  *cptr++ = marker_allele_ptrs[2 * uidx_cur + 1][0];
	}
      }
      for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
        next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	if (fputs_checked("# ID ", outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
        cptr = (char*)memchr(&(sample_ids[sample_uidx * max_sample_id_len]), '\t', max_sample_id_len);
        fputs(&(cptr[1]), outfile);
	putc_unlocked('\n', outfile);
	bufptr = &(loadbuf[sample_uidx / 4]);
	shiftval = (sample_uidx % 4) * 2;
	aptr = writebuf;
	aptr2 = writebuf2;
	uidx_cur = marker_uidx_start;
	marker_idx = 0;
        do {
          uidx_cur = next_unset_unsafe(marker_exclude, uidx_cur);
	  uidx_stop = next_set(marker_exclude, marker_uidx, chrom_end);
	  if (recode_012) {
	    marker_idx += uidx_stop - uidx_cur;
	    do {
	      ucc = ((*bufptr) >> shiftval) & 3;
	      *aptr++ = writebuf3[ucc];
	      *aptr2++ = writebuf3[ucc + 4];
	      bufptr = &(bufptr[unfiltered_sample_ct4]);
	    } while (++uidx_cur < uidx_stop);
	  } else {
	    uljj = marker_idx + uidx_stop - uidx_cur;
	    uidx_cur = uidx_stop;
	    do {
              ucc = ((*bufptr) >> shiftval) & 3;
	      // I'm sure there's a better way to do this, but at least this is
	      // better than derefencing marker_allele_ptrs in the innermost
	      // loop
	      if (ucc != 1) {
		*aptr++ = writebuf3[marker_idx * 2 + (ucc % 2)];
		*aptr2++ = writebuf3[marker_idx * 2 + (ucc / 2)];
	      } else {
		*aptr++ = '?';
		*aptr2++ = '?';
	      }
              bufptr = &(bufptr[unfiltered_sample_ct4]);
	    } while (++marker_idx < uljj);
	  }
	} while (marker_idx < ulii);
	fwrite(writebuf, 1, ulii, outfile);
	putc_unlocked('\n', outfile);
        fwrite(writebuf2, 1, ulii, outfile);
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (recode_modifier & RECODE_FASTPHASE) {
        if (chrom_idx > onechar_max) {
          putc_unlocked('\b', stdout);
	}
      }
      if (fclose_null(&outfile)) {
        goto recode_ret_WRITE_FAIL;
      }
    } while (chrom_fo_idx < last_chrom_fo_idx);
    if (recode_modifier & RECODE_FASTPHASE) {
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b               \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", stdout);
    }
  } else if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF)) {
    if (lgen_ref) {
      strcpy(outname_end, ".ref");
      if (fopen_checked(outname, "w", &outfile2)) {
	goto recode_ret_OPEN_FAIL;
      }
    }
    strcpy(outname_end, ".lgen");
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    if (delimiter == ' ') {
      sample_delim_convert(unfiltered_sample_ct, sample_exclude, sample_ct, max_sample_id_len, '\t', ' ', sample_ids);
    } else {
      if (!(recode_modifier & RECODE_DELIMX)) {
	delim2 = ' ';
      }
    }
    *outname_end = '\0';
    if (lgen_ref) {
      LOGPRINTFWW5("--recode lgen-ref to %s.lgen + %s.map + %s.fam + %s.ref ... ", outname, outname, outname, outname);
    } else {
      LOGPRINTFWW5("--recode lgen to %s.lgen + %s.map + %s.fam ... ", outname, outname, outname);
    }
    fputs("0%", stdout);
    fflush(stdout);
    ref_word = 4;
    writebuf[0] = delimiter;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	}
	if (sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, sample_include2, sample_male_include2, sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, sample_include2, sample_ct);
	  }
	}
	wbufptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	cptr = strcpya(&(writebuf[1]), wbufptr);
	cptr = memseta(cptr, delimiter, 2);
	ulii = (uintptr_t)(cptr - writebuf);
	alen = 2 * max_marker_allele_len + ulii;
	cur_mk_allelesx[0] = cptr;
	cur_mk_allelesx[1] = &(cptr[alen]);
	cur_mk_allelesx[2] = &(cptr[2 * alen]);
	cur_mk_allelesx[3] = &(cptr[3 * alen]);
	memcpy(&(writebuf[alen]), writebuf, ulii);
	memcpy(&(writebuf[2 * alen]), writebuf, ulii);
	memcpy(&(writebuf[3 * alen]), writebuf, ulii);
	cur_mk_allelesx[1][0] = output_missing_geno;
	cur_mk_allelesx[1][1] = delim2;
	cur_mk_allelesx[1][2] = output_missing_geno;
	cur_mk_allelesx[1][3] = '\n';
	init_recode_cmax(mk_allele_ptrs[2 * marker_uidx], mk_allele_ptrs[2 * marker_uidx + 1], cur_mk_allelesx, cmalen, '\n', delim2);
	cmalen[0] += ulii;
	cmalen[1] = 4 + ulii;
	cmalen[2] += ulii;
	cmalen[3] += ulii;
	if (lgen_ref) {
	  aptr = mk_allele_ptrs[2 * marker_uidx + 1];
	  aptr2 = mk_allele_ptrs[2 * marker_uidx];
	  if ((aptr[0] != missing_geno) || aptr[1] || (aptr2[0] != missing_geno) || aptr2[1]) {
	    fputs(wbufptr, outfile2);
	    if ((aptr[0] != missing_geno) || aptr[1]) {
	      putc_unlocked(delimiter, outfile2);
	      fputs(aptr, outfile2);
	    }
	    if ((aptr2[0] != missing_geno) || aptr2[1]) {
	      putc_unlocked(delimiter, outfile2);
	      fputs(aptr2, outfile2);
	    }
	    if (putc_checked('\n', outfile2)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  ref_word = 3;
	}

	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[sample_ct / BITCT2]);
	shiftmax = BITCT2;
	sample_idx = 0;
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; sample_idx < shiftmax; sample_idx++) {
	      ulii = cur_word & 3;
	      if (ulii != ref_word) {
		fputs(&(sample_ids_collapsed[sample_idx * max_sample_id_len]), outfile);
		fwrite(&(writebuf[ulii * alen]), 1, cmalen[ulii], outfile);
	      }
	      cur_word >>= 2;
	    }
	    shiftmax += BITCT2;
	  }
	  if (sample_idx == sample_ct) {
	    break;
	  }
	  ulptr_end++;
	  shiftmax = sample_ct;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    sample_delim_convert(unfiltered_sample_ct, sample_exclude, sample_ct, max_sample_id_len, ' ', '\t', sample_ids);
  } else if (recode_modifier & (RECODE_A | RECODE_AD)) {
    memcpy(outname_end, ".raw", 5);
    if (bigstack_left() < ((uint64_t)unfiltered_sample_ct4) * marker_ct) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    if (fputs_checked((delimiter == ' ')? "FID IID PAT MAT SEX PHENOTYPE" : "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE", outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      uii = IS_NONNULL_AND_SET(recode_allele_reverse, marker_uidx);
      if (allele_missing && allele_missing[marker_uidx]) {
	aptr = allele_missing[marker_uidx];
      } else {
	aptr = mk_allele_ptrs[2 * marker_uidx + uii];
      }
      putc_unlocked(delimiter, outfile);
      fputs(cptr, outfile);
      putc_unlocked('_', outfile);
      fputs(aptr, outfile);
      if (recode_modifier & RECODE_INCLUDE_ALT) {
	putc_unlocked('(', outfile);
	putc_unlocked('/', outfile);
	if (allele_missing && allele_missing[marker_uidx]) {
	  fputs(mk_allele_ptrs[2 * marker_uidx + uii], outfile);
	  putc_unlocked(',', outfile);
	}
	fputs(mk_allele_ptrs[2 * marker_uidx + 1 - uii], outfile);
	putc_unlocked(')', outfile);
      }
      if (recode_modifier & RECODE_AD) {
	putc_unlocked(delimiter, outfile);
	fputs(cptr, outfile);
	fputs("_HET", outfile);
      }
      if (ferror(outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    LOGPRINTFWW5("--recode A%s to %s ... ", (recode_modifier & RECODE_AD)? "D" : "", outname);
    if (!recode_allele_reverse) {
      recode_allele_reverse = marker_reverse;
    } else {
      bitvec_xor(marker_reverse, unfiltered_marker_ctl, recode_allele_reverse);
    }
    if (recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, recode_allele_reverse, &marker_uidx, unfiltered_sample_ct)) {
      goto recode_ret_READ_FAIL;
    }
    if ((set_hh_missing || set_mixed_mt_missing) && marker_ct) {
      haploid_fix_multiple(marker_exclude, 0, marker_ct, chrom_info_ptr, hh_exists, set_hh_missing, set_mixed_mt_missing, sample_include2, sample_male_include2, unfiltered_sample_ct, unfiltered_sample_ct4, loadbuf);
    }
    fputs("0%", stdout);
    sample_uidx = 0;
    sample_idx = 0;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = ((uint64_t)pct * sample_ct) / 100;
      for (; sample_idx < loop_end; sample_uidx++, sample_idx++) {
	sample_uidx = next_unset_ul_unsafe(sample_exclude, sample_uidx);
	if (recode_write_first_cols(outfile, sample_uidx, delimiter, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno)) {
	  goto recode_ret_WRITE_FAIL;
	}
	if (marker_ct) {
	  bufptr = &(loadbuf[sample_uidx / 4]);
	  wbufptr = writebuf;
	  shiftval = (sample_uidx % 4) * 2;
	  marker_uidx = 0;
	  marker_idx = 0;
	  if (recode_modifier & RECODE_A) {
	    do {
	      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	      ulii = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
	      marker_idx += ulii - marker_uidx;
	      do {
		ucc = ((*bufptr) >> shiftval) & 3;
		if (allele_missing && allele_missing[marker_uidx]) {
		  *wbufptr++ = "0N00"[ucc];
		} else {
		  *wbufptr++ = "2N10"[ucc];
		}
		if (ucc == 1) {
		  *wbufptr++ = 'A';
		}
		*wbufptr++ = delimiter;
		bufptr = &(bufptr[unfiltered_sample_ct4]);
	      } while (++marker_uidx < ulii);
	    } while (marker_idx < marker_ct);
	  } else {
	    do {
	      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	      ulii = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
	      marker_idx += ulii - marker_uidx;
	      do {
		ucc = ((*bufptr) >> shiftval) & 3;
		if (ucc != 1) {
		  wbufptr = memcpya(wbufptr, &(writebuf2[4 * ((allele_missing && allele_missing[marker_uidx])? 3 : ucc)]), 4);
		} else {
		  wbufptr = memcpya(wbufptr, &(writebuf2[16]), 6);
		}
		bufptr = &(bufptr[unfiltered_sample_ct4]);
	      } while (++marker_uidx < ulii);
	    } while (marker_idx < marker_ct);
	  }
	  wbufptr[-1] = '\n';
	  ulii = (uintptr_t)(wbufptr - writebuf);
	  if (fwrite_checked(writebuf, ulii, outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	} else {
	  putc_unlocked('\n', outfile);
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & (RECODE_LIST | RECODE_RLIST)) {
    strcpy(outname_end, rlist? ".rlist" : ".list");
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    if (delimiter != '\t') {
      if (bigstack_calloc_ul(sample_ctv2 / 2, &ulptr)) {
	goto recode_ret_NOMEM;
      }
      sample_delim_convert(sample_ct, ulptr, sample_ct, max_sample_id_len, '\t', ' ', sample_ids_collapsed);
      if (omit_nonmale_y) {
        sample_delim_convert(sample_ct_y, ulptr, sample_ct_y, max_sample_id_len, '\t', ' ', sample_ids_collapsed_y);
      }
    }
    if (rlist) {
      *outname_end = '\0';
      sprintf(g_logbuf, "--recode rlist to %s.rlist + %s.map + %s.fam ... ", outname, outname, outname);
    } else {
      sprintf(g_logbuf, "--recode list to %s ... ", outname);
    }
    wordwrapb(5);
    logprintb();
    fputs("0%", stdout);
    fflush(stdout);
    cur_mk_allelesx[1][0] = missing_geno;
    cur_mk_allelesx[1][1] = delimiter;
    cur_mk_allelesx[1][2] = missing_geno;
    cmalen[1] = 3;
    chrom_fo_idx = 0xffffffffU;
    chrom_end = 0;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  if (omit_nonmale_y && is_y) {
	    cur_final_mask = get_final_mask(sample_ct_y);
	    cur_sample_ct = sample_ct_y;
	    cur_sample_exclude = sample_exclude_y;
	    cur_sample_ids_collapsed = sample_ids_collapsed_y;
	    cur_sample_include2 = sample_include2_y;
	    cur_sample_male_include2 = sample_male_include2_y;
	  } else {
	    cur_final_mask = final_mask;
	    cur_sample_ct = sample_ct;
	    cur_sample_exclude = sample_exclude;
	    cur_sample_ids_collapsed = sample_ids_collapsed;
	    cur_sample_include2 = sample_include2;
	    cur_sample_male_include2 = sample_male_include2;
	  }
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	}
	if (cur_sample_ct) {
	  if (load_and_collapse(unfiltered_sample_ct, cur_sample_ct, cur_sample_exclude, cur_final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, (uintptr_t*)loadbuf, loadbuf_collapsed)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (is_haploid && set_hh_missing) {
	    haploid_fix(hh_exists, cur_sample_include2, cur_sample_male_include2, cur_sample_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	  } else if (is_mt && set_mixed_mt_missing) {
	    hh_reset((unsigned char*)loadbuf_collapsed, cur_sample_include2, cur_sample_ct);
	  }
	}
	init_recode_cmax(mk_allele_ptrs[2 * marker_uidx], mk_allele_ptrs[2 * marker_uidx + 1], cur_mk_allelesx, cmalen, '\0', delimiter);
	cmalen[0] -= 1;
	cmalen[2] -= 1;
	cmalen[3] -= 1;
	aptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	alen = strlen(aptr);
	for (ulii = 0; ulii < 4; ulii++) {
	  wbufptr = writebufl[ulii];
	  if (!rlist) {
	    wbufptr = chrom_name_write(chrom_info_ptr, chrom_idx, wbufptr);
	    *wbufptr++ = delimiter;
	  }
	  wbufptr = memcpyax(wbufptr, aptr, alen, delimiter);
	  if (rlist) {
	    switch (ulii) {
	    case 0:
	    case 3:
	      wbufptr = memcpyl3a(wbufptr, "HOM");
	      break;
	    case 1:
	      wbufptr = memcpyl3a(wbufptr, "NIL");
	      break;
	    case 2:
	      wbufptr = memcpyl3a(wbufptr, "HET");
	      break;
	    }
	    *wbufptr++ = delimiter;
	    wbufptr = memcpya(wbufptr, cur_mk_allelesx[ulii], cmalen[ulii]);
	  } else {
	    *wbufptr++ = cur_mk_allelesx[ulii][0];
	    *wbufptr++ = cur_mk_allelesx[ulii][2];
	  }
	  writebuflp[ulii] = wbufptr;
	}
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[cur_sample_ct / BITCT2]);
	sample_idx = 0;
	sample_uidx = BITCT2; // repurposed as stop value
	if (rlist) {
	  for (ulii = 0; ulii < 4; ulii++) {
	    writebuflps[ulii] = writebuflp[ulii];
	  }
	  while (1) {
            while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; sample_idx < sample_uidx; sample_idx++, cur_word >>= 2) {
                ulii = cur_word & 3;
		if (ulii != 3) {
		  *(writebuflp[ulii]++) = delimiter;
		  writebuflp[ulii] = strcpya(writebuflp[ulii], &(cur_sample_ids_collapsed[sample_idx * max_sample_id_len]));
		}
	      }
	      sample_uidx += BITCT2;
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
            ulptr_end++;
	    sample_uidx = cur_sample_ct;
	  }
	  if (writebuflp[2] != writebuflps[2]) {
	    *(writebuflp[2]++) = '\n';
	    if (fwrite_checked(writebufl[2], writebuflp[2] - writebufl[2], outfile)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  if (writebuflp[0] != writebuflps[0]) {
	    *(writebuflp[0]++) = '\n';
	    if (fwrite_checked(writebufl[0], writebuflp[0] - writebufl[0], outfile)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  if (writebuflp[1] != writebuflps[1]) {
	    *(writebuflp[1]++) = '\n';
	    if (fwrite_checked(writebufl[1], writebuflp[1] - writebufl[1], outfile)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  while (1) {
            while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; sample_idx < sample_uidx; sample_idx++, cur_word >>= 2) {
                ulii = cur_word & 3;
		*(writebuflp[ulii]++) = delimiter;
		writebuflp[ulii] = strcpya(writebuflp[ulii], &(cur_sample_ids_collapsed[sample_idx * max_sample_id_len]));
	      }
	      sample_uidx += BITCT2;
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
            ulptr_end++;
	    sample_uidx = cur_sample_ct;
	  }
	  *(writebuflp[0]++) = '\n';
	  *(writebuflp[1]++) = '\n';
	  *(writebuflp[2]++) = '\n';
	  *(writebuflp[3]++) = '\n';
	  if (fwrite_checked(writebufl[0], writebuflp[0] - writebufl[0], outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  if (fwrite_checked(writebufl[2], writebuflp[2] - writebufl[2], outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  if (fwrite_checked(writebufl[3], writebuflp[3] - writebufl[3], outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  if (fwrite_checked(writebufl[1], writebuflp[1] - writebufl[1], outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & (RECODE_HV | RECODE_HV_1CHR)) {
    if (bigstack_left() < ((uint64_t)unfiltered_sample_ct4) * max_chrom_size) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (!marker_ct) {
      logerrprint("Error: No variants for --recode HV{-1chr}.\n");
      goto recode_ret_ALL_MARKERS_EXCLUDED;
    }
    if (recode_modifier & RECODE_HV) {
      memcpy(outname_end, ".chr-", 5);
      sprintf(g_logbuf, "--recode HV to %s*.ped + .info... ", outname);
      wordwrapb(15); // strlen("[chromosome 10]");
      fputs(g_logbuf, stdout);
      fputs("[chromosome   ", stdout);
    } else {
      *outname_end = '\0';
      LOGPRINTFWW5("--recode HV-1chr to %s.ped + %s.info ... ", outname, outname);
      fflush(stdout);
    }
    chrom_fo_idx = 0xffffffffU;
    do {
      chrom_fo_idx++;
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      ulii = count_chrom_markers(chrom_info_ptr, marker_exclude, chrom_idx);
      if (recode_modifier & RECODE_HV) {
        wbufptr = chrom_name_write(chrom_info_ptr, chrom_idx, &(outname_end[5]));
        if (chrom_idx <= chrom_info_ptr->max_code) {
          printf("\b\b%u] \b", chrom_idx);
        } else {
          // nonstandard chromosome name
          fputs("\b\b**] \b", stdout);
        }
        fflush(stdout);
      } else {
        wbufptr = outname_end;
      }
      memcpy(wbufptr, ".ped", 5);
      if (fopen_checked(outname, "w", &outfile)) {
	goto recode_ret_OPEN_FAIL;
      }
      marker_uidx_start = marker_uidx;
      if (recode_load_to(loadbuf, bedfile, bed_offset, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1], 0, ulii, marker_exclude, marker_reverse, &marker_uidx, unfiltered_sample_ct)) {
	goto recode_ret_READ_FAIL;
      }
      if ((set_hh_missing || set_mixed_mt_missing) && marker_ct) {
        haploid_fix_multiple(marker_exclude, marker_uidx_start, ulii, chrom_info_ptr, hh_exists, set_hh_missing, set_mixed_mt_missing, sample_include2, sample_male_include2, unfiltered_sample_ct, unfiltered_sample_ct4, loadbuf);
      }
      sample_uidx = 0;
      // Haploview requires '0' missing phenotype code
      if (write_ped_lines(outfile, loadbuf, marker_exclude, marker_uidx_start, ulii, unfiltered_sample_ct4, sample_exclude, &sample_uidx, 0, sample_ct, 0, mk_allele_ptrs, delimiter, delim2, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_geno, &(g_one_char_strs[96]), writebuf)) {
	goto recode_ret_WRITE_FAIL;
      }
      if (fclose_null(&outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
      memcpy(wbufptr, ".info", 6);
      if (fopen_checked(outname, "w", &outfile)) {
	goto recode_ret_OPEN_FAIL;
      }
      if (write_haploview_map(outfile, marker_exclude, marker_uidx_start, ulii, marker_ids, max_marker_id_len, marker_pos)) {
        goto recode_ret_WRITE_FAIL;
      }
      if (fclose_null(&outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
      if (recode_modifier & RECODE_HV) {
        if (chrom_idx > onechar_max) {
          putc_unlocked('\b', stdout);
	}
	*wbufptr = '\0';
	LOGPREPRINTFWW("%s.ped + %s.info created.\n", outname, outname);
        logstr(g_logbuf);
      }
    } while (chrom_fo_idx < last_chrom_fo_idx);
    if (recode_modifier & RECODE_HV) {
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b               \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", stdout);
    }
  } else if (recode_modifier & RECODE_STRUCTURE) {
    memcpy(outname_end, ".recode.strct_in", 17);
    if (bigstack_left() < ((uint64_t)unfiltered_sample_ct4) * marker_ct) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    LOGPRINTFWW5("--recode structure to %s ... ", outname);
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      putc_unlocked(' ', outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    chrom_fo_idx = 0xffffffffU; // exploit overflow
    chrom_end = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      if (marker_uidx >= chrom_end) {
	do {
          chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1];
	} while (marker_uidx >= chrom_end);
	fputs("-1 ", outfile);
      } else {
        wbufptr = uint32toa_x(marker_pos[marker_uidx] - last_pos, ' ', writebuf2);
        fwrite(writebuf2, 1, wbufptr - writebuf2, outfile);
      }
      last_pos = marker_pos[marker_uidx];
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    if (recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, marker_reverse, &marker_uidx, unfiltered_sample_ct)) {
      goto recode_ret_READ_FAIL;
    }
    if ((set_hh_missing || set_mixed_mt_missing) && marker_ct) {
      haploid_fix_multiple(marker_exclude, 0, marker_ct, chrom_info_ptr, hh_exists, set_hh_missing, set_mixed_mt_missing, sample_include2, sample_male_include2, unfiltered_sample_ct, unfiltered_sample_ct4, loadbuf);
    }
    sample_uidx = 0;
    sample_idx = 0;
    fputs("0%", stdout);
    last_pos = 0; // now unique FIDs encountered so far
    memcpy(writebuf2, " 1 1 0 0 1 2 2 2", 16);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * sample_ct) / 100;
      for (; sample_idx < loop_end; sample_uidx++, sample_idx++) {
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	bufptr = &(loadbuf[sample_uidx / 4]);
	shiftval = (sample_uidx % 4) * 2;
        cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
        aptr = (char*)memchr(cptr, '\t', max_sample_id_len);
        fputs(&(aptr[1]), outfile);
        ii = bsearch_str(cptr, (uintptr_t)(aptr - cptr), writebuf3, max_fid_len, fid_ct);
	cur_fid = fid_map[(uint32_t)ii];
	if (!cur_fid) {
	  cur_fid = ++last_pos;
          fid_map[(uint32_t)ii] = last_pos;
	}
	g_textbuf[0] = ' ';
        wbufptr = uint32toa(cur_fid, &(g_textbuf[1]));
        fwrite(g_textbuf, 1, wbufptr - g_textbuf, outfile);
	marker_uidx = 0;
	wbufptr = writebuf;
        for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
          ucc = ((*bufptr) >> shiftval) & 3;
          wbufptr = memcpya(wbufptr, &(writebuf2[4 * ucc]), 4);
	  bufptr = &(bufptr[unfiltered_sample_ct4]);
	}
	if (fwrite_checked(writebuf, wbufptr - writebuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	putc_unlocked('\n', outfile);
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else {
    memcpy(outname_end, ".ped", 5);
    if (bigstack_left() < ((uint64_t)unfiltered_sample_ct4) * marker_ct) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto recode_ret_OPEN_FAIL;
    }
    *outname_end = '\0';
    LOGPRINTFWW5("--recode ped to %s.ped + %s.map ... ", outname, outname);
    if (recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, marker_reverse, &marker_uidx, unfiltered_sample_ct)) {
      goto recode_ret_READ_FAIL;
    }
    if ((set_hh_missing || set_mixed_mt_missing) && marker_ct) {
      haploid_fix_multiple(marker_exclude, 0, marker_ct, chrom_info_ptr, hh_exists, set_hh_missing, set_mixed_mt_missing, sample_include2, sample_male_include2, unfiltered_sample_ct, unfiltered_sample_ct4, loadbuf);
    }
    sample_uidx = 0;
    sample_idx = 0;
    fputs("0%", stdout);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * sample_ct) / 100;
      if (write_ped_lines(outfile, loadbuf, marker_exclude, 0, marker_ct, unfiltered_sample_ct4, sample_exclude, &sample_uidx, sample_idx, loop_end, recode_modifier & RECODE_COMPOUND, mk_allele_ptrs, delimiter, delim2, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_geno, output_missing_pheno, writebuf)) {
	goto recode_ret_WRITE_FAIL;
      }
      sample_idx = loop_end;
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  }
  if (outfile) {
    if (fclose_null(&outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
  }

  if (recode_modifier & (RECODE_TRANSPOSE | RECODE_LGEN | RECODE_LGEN_REF | RECODE_RLIST)) {
    if (recode_modifier & RECODE_TRANSPOSE) {
      strcpy(outname_end, ".tfam");
    } else {
      strcpy(outname_end, ".fam");
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	delimiter = ' ';
      }
    }
    retval = write_fam(outname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno, delimiter, nullptr);
    if (retval) {
      goto recode_ret_1;
    }
  }

  // bugfix: stop generating a .map file when Oxford-format requested
  if (recode_modifier & (RECODE_COMPOUND | RECODE_LGEN | RECODE_LGEN_REF | RECODE_PED | RECODE_RLIST)) {
    strcpy(outname_end, ".map");
    retval = write_map_or_bim(outname, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, nullptr, ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_DELIMX)? ' ' : '\t', chrom_info_ptr);
    if (retval) {
      goto recode_ret_1;
    }
  }
  if (!(recode_modifier & (RECODE_23 | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_HV | RECODE_HV_1CHR))) {
    fputs("\b\b\b", stdout);
  }
  if (!(recode_modifier & (RECODE_BEAGLE | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_HV))) {
    logprint("done.\n");
    if (invalid_allele_code_seen) {
      logerrprint("Warning: At least one VCF allele code violates the official specification;\nother tools may not accept the file.  (Valid codes must either start with a\n'<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or\nrepresent a breakend.)\n");
    }
  } else {
    fputs("done.\n", stdout);
    if (recode_modifier & RECODE_BEAGLE) {
      logstr("--recode beagle complete.\n");
    } else if (recode_modifier & RECODE_HV) {
      logstr("--recode HV complete.\n");
    } else {
      sprintf(g_logbuf, "--recode fastphase%s complete.\n", (recode_modifier & RECODE_FASTPHASE_1CHR)? "-1chr" : "");
      logstr(g_logbuf);
    }
  }
  while (0) {
  recode_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  recode_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  recode_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  recode_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  recode_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  recode_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  recode_ret_NO_MULTIPASS_YET:
    // probably want to implement this later
    logerrprint("Error: --recode does not yet support multipass recoding of very large files;\ncontact the " PROG_NAME_CAPS " developers if you need this.\nFor now, you can try using a machine with more memory, and/or split the file\ninto smaller pieces and recode them separately.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    break;
  recode_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  recode_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
 recode_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  fclose_cond(outfile2);
  fclose_cond(outfile);
  if (bgz_outfile) {
    bgzf_close(bgz_outfile);
  }
  flex_pzwrite_close_cond(&ps, pzwritep);
  return retval;
}

int32_t sample_sort_file_map(char* sample_sort_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t** sample_sort_map_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  FILE* infile = nullptr;
  // temporary: sample_id_map[ascii-sorted idx] = uidx in input fileset
  uint32_t* sample_id_map = nullptr;
  uintptr_t line_idx = 0;
  uint32_t cur_seq = 0;
  int32_t retval = 0;
  uintptr_t* already_seen;

  // return map: sample_sort_map[final uidx] = uidx in input fileset
  uint32_t* sample_sort_map;

  char* sorted_sample_ids;
  char* bufptr;
  char* idbuf;
  int32_t ii;
  if (sample_exclude) {
    // called from plink()
    if (bigstack_alloc_ui(sample_ct, &sample_sort_map)) {
      goto sample_sort_file_map_ret_NOMEM;
    }
  } else {
    // called from merge_sample_sortf()
    sample_sort_map = *sample_sort_map_ptr;
    sorted_sample_ids = sample_ids;
  }
  if (bigstack_alloc_c(max_sample_id_len, &idbuf) ||
      bigstack_calloc_ul(sample_ctl, &already_seen)) {
    goto sample_sort_file_map_ret_NOMEM;
  }
  if (sample_exclude) {
    retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, unfiltered_sample_ct - sample_ct, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &sorted_sample_ids, &sample_id_map);
    if (retval) {
      goto sample_sort_file_map_ret_1;
    }
  }
  if (fopen_checked(sample_sort_fname, "r", &infile)) {
    goto sample_sort_file_map_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --indiv-sort file is pathologically long.\n", line_idx);
      goto sample_sort_file_map_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, nullptr, &ii, idbuf)) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --indiv-sort file has fewer tokens than expected.\n", line_idx);
      goto sample_sort_file_map_ret_INVALID_FORMAT_2;
    }
    if (ii != -1) {
      if (is_set(already_seen, ii)) {
        strchr(idbuf, '\t')[0] = ' ';
        LOGPREPRINTFWW("Error: Duplicate ID '%s' in --indiv-sort file.\n", idbuf);
        goto sample_sort_file_map_ret_INVALID_FORMAT_2;
      }
      set_bit(ii, already_seen);
      if (sample_id_map) {
        sample_sort_map[cur_seq] = sample_id_map[(uint32_t)ii];
      } else {
	sample_sort_map[cur_seq] = (uint32_t)ii;
      }
      cur_seq++;
    }
  }
  if (fclose_null(&infile)) {
    goto sample_sort_file_map_ret_READ_FAIL;
  }
  if (cur_seq != sample_ct) {
    logerrprint("Error: --indiv-sort file does not contain all loaded sample IDs.\n");
    goto sample_sort_file_map_ret_INVALID_CMDLINE;
  }
  *sample_sort_map_ptr = sample_sort_map;
  bigstack_mark = (unsigned char*)idbuf;
  while (0) {
  sample_sort_file_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sample_sort_file_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sample_sort_file_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  sample_sort_file_map_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  sample_sort_file_map_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 sample_sort_file_map_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(infile);
  return retval;
}

// .fam
typedef struct ll_entry_struct {
  struct ll_entry_struct* next;
  double pheno;
  uint32_t orig_order;
  char idstr[];
} Ll_fam;

// .bim
typedef struct ll_entry2_struct {
  struct ll_entry2_struct* next;
  int64_t pos;
  double cm;
  char* allele[2];
  char idstr[];
} Ll_bim;

static inline int32_t idmatch(char* idtab, char* id0, uint32_t id0_len_p1, char* id1, uint32_t id1_len_p1) {
  return (!(memcmp(idtab, id0, id0_len_p1) || memcmp(&(idtab[id0_len_p1]), id1, id1_len_p1)));
}

static inline uint32_t hashval(char* id1, uint32_t id1_len, char* id2, uint32_t id2_len) {
  // just interpret as little-endian number and take modulo HASHSIZE_S for now.
  // may want to upgrade this later
  unsigned char* ucptr = (unsigned char*)id1;
  unsigned char* ucp_end = &(ucptr[id1_len]);
  uint32_t vv = *ucptr;
  while (++ucptr != ucp_end) {
    vv = ((vv << 8) + (*ucptr)) % HASHSIZE_S;
  }
  vv = ((vv << 8) + 9) % HASHSIZE_S;
  ucptr = (unsigned char*)id2;
  ucp_end = &(ucptr[id2_len]);
  do {
    vv = ((vv << 8) + (*ucptr++)) % HASHSIZE_S;
  } while (ucptr != ucp_end);
  return vv;
}

static inline int32_t bigstack_end_alloc_llfam(uintptr_t idstr_bytes, Ll_fam** llfamp_ptr) {
  *llfamp_ptr = (Ll_fam*)bigstack_end_alloc(idstr_bytes + sizeof(Ll_fam));
  return !(*llfamp_ptr);
}

static inline int32_t bigstack_end_alloc_llbim(uintptr_t idstr_bytes, Ll_bim** llbimp_ptr) {
  *llbimp_ptr = (Ll_bim*)bigstack_end_alloc(idstr_bytes + sizeof(Ll_bim));
  return !(*llbimp_ptr);
}

int32_t merge_fam_id_scan(char* bedname, char* famname, uint32_t allow_no_samples, uintptr_t* max_sample_id_len_ptr, uint32_t* max_sample_full_len_ptr, uint32_t* is_dichot_pheno_ptr, Ll_fam** htable_fam, uint64_t* tot_sample_ct_ptr, uint32_t* ped_buflen_ptr, uint32_t* cur_sample_ct_ptr, uint32_t* orig_idx_ptr) {
  uint64_t tot_sample_ct = *tot_sample_ct_ptr;
  uintptr_t max_sample_id_len = *max_sample_id_len_ptr;
  uintptr_t line_idx = 0;
  FILE* infile = nullptr;
  uint32_t max_sample_full_len = *max_sample_full_len_ptr;
  uint32_t is_dichot_pheno = *is_dichot_pheno_ptr;
  uint32_t orig_idx = *orig_idx_ptr;
  uint32_t cur_sample_ct = 0;
  uint32_t text_file = 0;
  int32_t retval = 0;
  uint32_t col1_len;
  uint32_t col2_len;
  uint32_t col3_len;
  uint32_t col4_len;
  uint32_t tot_len;
  uintptr_t ulii;
  uint32_t uii;
  Ll_fam** llfam_pptr;
  Ll_fam* llfam_ptr;
  char* col2_start_ptr;
  char* col3_start_ptr;
  char* col4_start_ptr;
  char* col5_start_ptr;
  char* col6_start_ptr;
  char* col1_end_ptr;
  char* col2_end_ptr;
  char* col1_start_ptr;
  char* wptr;
  double pheno;
  char cc;
  if (!famname) {
    famname = bedname;
    text_file = 1;
  }
  infile = fopen(famname, "r");
  if (!infile) {
    uii = strlen(famname);
    if ((!orig_idx) && (uii > 8) && ((!memcmp(&(famname[uii - 8]), ".bed.fam", 8)) || (!memcmp(&(famname[uii - 8]), ".bim.fam", 8)) || (!memcmp(&(famname[uii - 8]), ".fam.fam", 8)))) {
      // technically could be --merge-list with no --bfile, but we won't bother
      // with a specialized error message for that case.
      LOGERRPRINTFWW("Error: Failed to open %s. (--bfile expects a filename *prefix*; '.bed', '.bim', and '.fam' are automatically appended.)\n", famname);
    } else {
      LOGERRPRINTFWW(g_errstr_fopen, famname);
    }
    goto merge_fam_id_scan_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    col1_start_ptr = skip_initial_spaces(g_textbuf);
    cc = *col1_start_ptr;
    if (!is_eoln_or_comment_kns(cc)) {
      col1_end_ptr = token_endnn(col1_start_ptr);
      col1_len = col1_end_ptr - col1_start_ptr;
      col2_start_ptr = skip_initial_spaces(col1_end_ptr);
      col2_len = strlen_se(col2_start_ptr);
      col2_end_ptr = &(col2_start_ptr[col2_len]);
      col3_start_ptr = skip_initial_spaces(col2_end_ptr);
      col3_len = strlen_se(col3_start_ptr);
      col4_start_ptr = skip_initial_spaces(&(col3_start_ptr[col3_len]));
      col4_len = strlen_se(col4_start_ptr);
      col5_start_ptr = skip_initial_spaces(&(col4_start_ptr[col4_len]));
      uii = strlen_se(col5_start_ptr);
      col6_start_ptr = skip_initial_spaces(&(col5_start_ptr[uii]));
      if (no_more_tokens_kns(col6_start_ptr)) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, famname);
	if (text_file) {
	  // .ped/.map swap?  give a more helpful error message then
	  uii = strlen(famname);
	  if (!memcmp(&(famname[uii - 4]), ".map", 4)) {
	    logprintb();
	    logprint("(The .ped file must be named before the .map; did you swap them?)\n");
	    goto merge_fam_id_scan_ret_INVALID_FORMAT;
	  }
	}
	goto merge_fam_id_scan_ret_INVALID_FORMAT_2;
      }
      if (uii != 1) {
	*col5_start_ptr = '0';
      }
      *col1_end_ptr = '\t';
      *col2_end_ptr = '\t';
      uii = col1_len + col2_len + 2;
      if (uii > max_sample_id_len) {
	max_sample_id_len = uii;
      }
      tot_len = uii + col3_len + col4_len + 4;
      uii = hashval(col1_start_ptr, col1_len, col2_start_ptr, col2_len);
      llfam_pptr = &(htable_fam[uii]);
      llfam_ptr = *llfam_pptr;
      uii = 1;
      if (is_dichot_pheno) {
	is_dichot_pheno = eval_affection(col6_start_ptr, -9);
      }
      if (scan_double(col6_start_ptr, &pheno)) {
	pheno = -9;
      }
      while (llfam_ptr) {
	if (idmatch(llfam_ptr->idstr, col1_start_ptr, col1_len + 1, col2_start_ptr, col2_len + 1)) {
	  uii = 0;
	  /*
	  // possibly for future: add parental ID/sex merge (not in PLINK 1.07)
	  if (merge_mode == 1) {
	    if (fabs(pheno - llfam_ptr->pheno) > PHENO_EPSILON) {
	      llfam_ptr->pheno = -9;
	    }
	  } else if (merge_mode == 2) {
	    if (llfam_ptr->pheno == -9) {
	      llfam_ptr->pheno = pheno;
	    }
	  } else if ((merge_mode == 5) || ((merge_mode == 3) && (pheno != -9))) {
	    llfam_ptr->pheno = pheno;
	  }
	  */
	  break;
	}
        llfam_pptr = &(llfam_ptr->next);
	llfam_ptr = *llfam_pptr;
      }
      if (uii) {
	if (tot_len > max_sample_full_len) {
	  max_sample_full_len = tot_len;
	}
	if (bigstack_end_alloc_llfam(tot_len, &llfam_ptr)) {
	  goto merge_fam_id_scan_ret_NOMEM;
	}
	llfam_ptr->next = nullptr;
	llfam_ptr->pheno = pheno;
	llfam_ptr->orig_order = orig_idx++;
	wptr = memcpyax(memcpyax(memcpyax(memcpyax(llfam_ptr->idstr, col1_start_ptr, col1_len, '\t'), col2_start_ptr, col2_len, '\t'), col3_start_ptr, col3_len, '\t'), col4_start_ptr, col4_len, '\t');
	*wptr = *col5_start_ptr;
	wptr[1] = '\0';
	*llfam_pptr = llfam_ptr;
	tot_sample_ct++;
      }
      cur_sample_ct++;
    }
    if (!g_textbuf[MAXLINELEN - 1]) {
      if (!text_file) {
	goto merge_fam_id_scan_ret_LONG_LINE;
      }
      ulii = 0;
      do {
	g_textbuf[MAXLINELEN - 1] = ' ';
	if (g_textbuf[MAXLINELEN - 2] == '\n') {
	  break;
	}
	ulii += MAXLINELEN - 1;
	if (ulii >= MAXLINEBUFLEN) {
	  goto merge_fam_id_scan_ret_LONG_LINE;
	}
        if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	  goto merge_fam_id_scan_ret_READ_FAIL;
	}
      } while (!g_textbuf[MAXLINELEN - 1]);
      ulii += strlen(g_textbuf) + 1;
      if (ulii > (*ped_buflen_ptr)) {
	*ped_buflen_ptr = ulii;
      }
    }
  }
  if (!feof(infile)) {
    goto merge_fam_id_scan_ret_READ_FAIL;
  }
  if ((!cur_sample_ct) && (!allow_no_samples)) {
    LOGPREPRINTFWW("Error: No %s in %s.\n", g_species_plural, famname);
    goto merge_fam_id_scan_ret_INVALID_FORMAT_2;
  }
  *max_sample_id_len_ptr = max_sample_id_len;
  *max_sample_full_len_ptr = max_sample_full_len;
  *is_dichot_pheno_ptr = is_dichot_pheno;
  *tot_sample_ct_ptr = tot_sample_ct;
  *cur_sample_ct_ptr = cur_sample_ct;
  *orig_idx_ptr = orig_idx;
  while (0) {
  merge_fam_id_scan_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_fam_id_scan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_fam_id_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_fam_id_scan_ret_LONG_LINE:
    LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, famname);
  merge_fam_id_scan_ret_INVALID_FORMAT_2:
    logerrprintb();
  merge_fam_id_scan_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
  }
  fclose_cond(infile);
  return retval;
}

int32_t merge_sample_sortf(char* sample_sort_fname, char* sample_fids, uintptr_t tot_sample_ct, uintptr_t max_sample_full_len, char* sample_ids, uintptr_t max_sample_id_len, uint32_t* map_reverse) {
  // sample_fids[] is already sorted
  unsigned char* bigstack_mark = g_bigstack_base;
  int32_t retval = 0;
  uintptr_t sample_uidx;
  for (sample_uidx = 0; sample_uidx < tot_sample_ct; sample_uidx++) {
    strcpy(&(sample_ids[sample_uidx * max_sample_id_len]), &(sample_fids[sample_uidx * max_sample_full_len]));
  }
  retval = sample_sort_file_map(sample_sort_fname, tot_sample_ct, nullptr, tot_sample_ct, sample_ids, max_sample_id_len, &map_reverse);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t merge_bim_scan(char* bimname, uint32_t is_binary, uint32_t allow_no_variants, uintptr_t* max_marker_id_len_ptr, Ll_bim** htable_bim, uint32_t* max_bim_linelen_ptr, uint64_t* tot_marker_ct_ptr, uint32_t* cur_marker_ct_ptr, uint64_t* position_warning_ct_ptr, Ll_str** non_biallelics_ptr, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t max_marker_id_len = *max_marker_id_len_ptr;
  uintptr_t loadbuf_size = MAXLINELEN;
  uint32_t max_bim_linelen = *max_bim_linelen_ptr;
  uint64_t tot_marker_ct = *tot_marker_ct_ptr;
  uint64_t position_warning_ct = *position_warning_ct_ptr;
  uint32_t cur_marker_ct = 0;
  double cm = 0.0;
  FILE* infile = nullptr;
  int32_t retval = 0;
  uint32_t alen1 = 1;
  uint32_t alen2 = 1;
  char* cur_alleles[2];
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* aptr1;
  char* aptr2;
  char* new_aptr;
  Ll_bim** llbim_pptr;
  Ll_bim* llbim_ptr;
  Ll_str* llstr_new_ptr;
  int64_t llxx;
  uintptr_t line_idx;
  uint32_t cm_col_exists;
  uint32_t allele_ct;
  uint32_t name_match;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t jj;
  {
    if (fopen_checked(bimname, "r", &infile)) {
      goto merge_bim_scan_ret_OPEN_FAIL;
    }
    if (is_binary) {
      loadbuf_size = (bigstack_left() / 2) & (~(CACHELINE - ONELU));
      if (bigstack_left() > 0x3fffffc0) {
	loadbuf_size = 0x3fffffc0;
      } else if (loadbuf_size <= MAXLINELEN) {
	goto merge_bim_scan_ret_NOMEM;
      }
    }
    bigstack_alloc_c(loadbuf_size, &loadbuf);
    loadbuf[loadbuf_size - 1] = ' ';
    if (check_cm_col(infile, loadbuf, is_binary, allow_no_variants, loadbuf_size, &cm_col_exists, &line_idx)) {
      goto merge_bim_scan_ret_MISSING_TOKENS;
    }
    if (!line_idx) {
      // no variants
      *cur_marker_ct_ptr = 0;
      goto merge_bim_scan_ret_1;
    }
    line_idx--;
    do {
      line_idx++;
      if (!loadbuf[loadbuf_size - 1]) {
	if ((loadbuf_size == 0x3fffffc0) || ((!is_binary) && (loadbuf_size == MAXLINELEN))) {
	  LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, bimname);
	  goto merge_bim_scan_ret_INVALID_FORMAT_2;
	} else {
	  goto merge_bim_scan_ret_NOMEM;
	}
      }
      uii = strlen(loadbuf);
      if (uii >= max_bim_linelen) {
	max_bim_linelen = uii + 1;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_or_comment_kns(*bufptr)) {
	continue;
      }
      char* chrom_token_end = token_endnn(bufptr);
      if (!(*chrom_token_end)) {
	goto merge_bim_scan_ret_MISSING_TOKENS;
      }
      int32_t cur_chrom_code;
      retval = get_or_add_chrom_code_destructive(bimname, line_idx, allow_extra_chroms, bufptr, chrom_token_end, chrom_info_ptr, &cur_chrom_code);
      if (retval) {
	goto merge_bim_scan_ret_1;
      }
      // do not filter on chrom_mask here, since that happens later
      bufptr = skip_initial_spaces(&(chrom_token_end[1]));
      if (is_eoln_kns(*bufptr)) {
	goto merge_bim_scan_ret_MISSING_TOKENS;
      }
      bufptr2 = token_endnn(bufptr);
      uii = bufptr2 - bufptr;
      bufptr2 = skip_initial_spaces(bufptr2);
      if (is_eoln_kns(*bufptr2)) {
	goto merge_bim_scan_ret_MISSING_TOKENS;
      }
      if (cm_col_exists) {
	if (scan_double(bufptr2, &cm)) {
	  cm = 0;
	}
	bufptr2 = next_token(bufptr2);
	if (no_more_tokens_kns(bufptr2)) {
	  goto merge_bim_scan_ret_MISSING_TOKENS;
	}
      }
      if (scan_int_abs_defcap(bufptr2, &jj)) {
	LOGPREPRINTFWW("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, bimname);
	goto merge_bim_scan_ret_INVALID_FORMAT_2;
      }
      if (jj >= 0) {
	if (is_binary) {
	  aptr1 = next_token(bufptr2);
	  aptr2 = next_token(aptr1);
	  if (no_more_tokens_kns(aptr2)) {
	    goto merge_bim_scan_ret_MISSING_TOKENS;
	  }
	  alen1 = strlen_se(aptr1);
	  alen2 = strlen_se(aptr2);
	  aptr1[alen1] = '\0';
	  aptr2[alen2] = '\0';
	  if ((alen1 == 1) && (*aptr1 == '0')) {
	    aptr1 = nullptr;
	  }
	  if (aptr1 && (alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    LOGPREPRINTFWW("Error: Identical A1 and A2 alleles on line %" PRIuPTR " of %s.\n", line_idx, bimname);
	    goto merge_bim_scan_ret_INVALID_FORMAT_2;
	  }
	  if ((alen2 == 1) && (*aptr2 == '0')) {
	    aptr2 = nullptr;
	  }
	} else {
	  aptr1 = nullptr;
	  aptr2 = nullptr;
	}
	llxx = (((uint64_t)((uint32_t)cur_chrom_code)) << 32) + ((uint32_t)jj);
	ujj = hashval2(bufptr, uii);
	llbim_pptr = &(htable_bim[ujj]);
	llbim_ptr = *llbim_pptr;
	name_match = 0;
	bufptr[uii++] = '\0';
	while (llbim_ptr) {
	  if (!strcmp(llbim_ptr->idstr, bufptr)) {
	    if (is_binary) {
	      bufptr2 = llbim_ptr->allele[0];
	      allele_ct = 0;
	      if (bufptr2) {
		cur_alleles[0] = bufptr2;
		allele_ct = 1;
	      }
	      bufptr3 = llbim_ptr->allele[1];
	      if (bufptr3) {
		cur_alleles[allele_ct++] = bufptr3;
	      }
	      if (aptr2) {
		for (ukk = 0; ukk < allele_ct; ukk++) {
		  if (!strcmp(aptr2, cur_alleles[ukk])) {
		    break;
		  }
		}
		if (ukk == allele_ct) {
		  if (allele_ct == 2) {
		    if (bigstack_end_alloc_llstr(uii, &llstr_new_ptr)) {
		      goto merge_bim_scan_ret_NOMEM;
		    }
		    llstr_new_ptr->next = *non_biallelics_ptr;
		    memcpy(llstr_new_ptr->ss, bufptr, uii);
		    *non_biallelics_ptr = llstr_new_ptr;
		  } else {
		    if (allele_set(aptr2, alen2, &new_aptr)) {
		      goto merge_bim_scan_ret_NOMEM;
		    }
		    if (!llbim_ptr->allele[1]) {
		      llbim_ptr->allele[1] = new_aptr;
		    } else {
		      llbim_ptr->allele[0] = new_aptr;
		    }
		    cur_alleles[allele_ct++] = new_aptr;
		  }
		}
	      }
	      if (aptr1) {
		for (ukk = 0; ukk < allele_ct; ukk++) {
		  if (!strcmp(aptr1, cur_alleles[ukk])) {
		    break;
		  }
		}
		if (ukk == allele_ct) {
		  if (allele_ct == 2) {
		    if (bigstack_end_alloc_llstr(uii, &llstr_new_ptr)) {
		      goto merge_bim_scan_ret_NOMEM;
		    }
		    llstr_new_ptr->next = *non_biallelics_ptr;
		    memcpy(llstr_new_ptr->ss, bufptr, uii);
		    *non_biallelics_ptr = llstr_new_ptr;
		  } else {
		    if (allele_set(aptr1, alen1, &new_aptr)) {
		      goto merge_bim_scan_ret_NOMEM;
		    }
		    if (!llbim_ptr->allele[1]) {
		      llbim_ptr->allele[1] = new_aptr;
		    } else {
		      llbim_ptr->allele[0] = new_aptr;
		    }
		    cur_alleles[allele_ct++] = new_aptr;
		  }
		}
	      }
	    }
	    if (llbim_ptr->pos != llxx) {
	      if ((((uint64_t)llbim_ptr->pos) >> 32) == (((uint64_t)llxx) >> 32)) {
		LOGPREPRINTFWW("Warning: Multiple positions seen for variant '%s'.\n", bufptr);
		if (position_warning_ct < 3) {
		  logerrprintb();
		} else {
		  logstr(g_logbuf);
		}
		position_warning_ct++;
	      } else {
		LOGERRPRINTFWW("Warning: Multiple chromosomes seen for variant '%s'.\n", bufptr);
	      }
	    }
	    name_match = 1;
	    break;
	  }
	  llbim_pptr = &(llbim_ptr->next);
	  llbim_ptr = *llbim_pptr;
	}
	if (!name_match) {
	  if (uii > max_marker_id_len) {
	    max_marker_id_len = uii;
	  }
	  if (bigstack_end_alloc_llbim(uii, &llbim_ptr)) {
	    goto merge_bim_scan_ret_NOMEM;
	  }
	  llbim_ptr->next = nullptr;
	  llbim_ptr->pos = llxx;
	  llbim_ptr->cm = cm;
	  if (aptr1) {
	    if (allele_set(aptr1, alen1, &(llbim_ptr->allele[0]))) {
	      goto merge_bim_scan_ret_NOMEM;
	    }
	  } else {
	    llbim_ptr->allele[0] = nullptr;
	  }
	  if (aptr2) {
	    if (allele_set(aptr2, alen2, &(llbim_ptr->allele[1]))) {
	      goto merge_bim_scan_ret_NOMEM;
	    }
	  } else {
	    llbim_ptr->allele[1] = nullptr;
	  }
	  memcpy(llbim_ptr->idstr, bufptr, uii);
	  *llbim_pptr = llbim_ptr;
	  tot_marker_ct++;
	}
	cur_marker_ct++;
      }
    } while (fgets(loadbuf, loadbuf_size, infile));
    if (!feof(infile)) {
      goto merge_bim_scan_ret_READ_FAIL;
    }
    if (!cur_marker_ct) {
      LOGPREPRINTFWW("Error: No variants in %s.\n", bimname);
      goto merge_bim_scan_ret_INVALID_FORMAT_2;
    }
    *max_marker_id_len_ptr = max_marker_id_len;
    *max_bim_linelen_ptr = max_bim_linelen;
    *tot_marker_ct_ptr = tot_marker_ct;
    *cur_marker_ct_ptr = cur_marker_ct;
    *position_warning_ct_ptr = position_warning_ct;
  }

  while (0) {
  merge_bim_scan_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_bim_scan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_bim_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_bim_scan_ret_MISSING_TOKENS:
    LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, bimname);
  merge_bim_scan_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
  }
 merge_bim_scan_ret_1:
  fclose_cond(infile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t report_non_biallelics(char* outname, char* outname_end, Ll_str* non_biallelics) {
  FILE* outfile = nullptr;
  int32_t retval = 0;
  uintptr_t nbmarker_ct_dup = 0;
  uintptr_t nbmarker_ct = 1;
  uintptr_t max_nbmarker_id_len = 0;
  Ll_str* cur_ptr = non_biallelics;
  uintptr_t cur_slen;
  char* id_arr;
  char* id_arr_ptr;
  char* id_arr_ptr_old;
  char* id_arr_end;
  do {
    cur_slen = strlen(cur_ptr->ss);
    if (cur_slen >= max_nbmarker_id_len) {
      max_nbmarker_id_len = cur_slen + 1;
    }
    nbmarker_ct_dup++;
    cur_ptr = cur_ptr->next;
  } while (cur_ptr);
  if (bigstack_alloc_c(nbmarker_ct_dup * max_nbmarker_id_len, &id_arr)) {
    goto report_non_biallelics_ret_NOMEM;
  }
  cur_ptr = non_biallelics;
  id_arr_ptr = id_arr;
  do {
    strcpy(id_arr_ptr, cur_ptr->ss);
    cur_ptr = cur_ptr->next;
    id_arr_ptr = &(id_arr_ptr[max_nbmarker_id_len]);
  } while (cur_ptr);
  qsort(id_arr, nbmarker_ct_dup, max_nbmarker_id_len, strcmp_casted);
  memcpy(outname_end, ".missnp", 8);
  if (fopen_checked(outname, "w", &outfile)) {
    goto report_non_biallelics_ret_OPEN_FAIL;
  }
  id_arr_ptr = id_arr;
  if (fputs_checked(id_arr_ptr, outfile)) {
    goto report_non_biallelics_ret_WRITE_FAIL;
  }
  putc_unlocked('\n', outfile);
  id_arr_end = &(id_arr[(nbmarker_ct_dup - 1) * max_nbmarker_id_len]);
  while (id_arr_ptr != id_arr_end) {
    id_arr_ptr_old = id_arr_ptr;
    id_arr_ptr = &(id_arr_ptr[max_nbmarker_id_len]);
    if (strcmp(id_arr_ptr, id_arr_ptr_old)) {
      nbmarker_ct++;
      fputs(id_arr_ptr, outfile);
      if (putc_checked('\n', outfile)) {
        goto report_non_biallelics_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto report_non_biallelics_ret_WRITE_FAIL;
  }
  LOGERRPRINTF("Error: %" PRIuPTR " variant%s with 3+ alleles present.\n* If you believe this is due to strand inconsistency, try --flip with\n  %s.\n  (Warning: if this seems to work, strand errors involving SNPs with A/T or C/G\n  alleles probably remain in your data.  If LD between nearby SNPs is high,\n  --flip-scan should detect them.)\n* If you are dealing with genuine multiallelic variants, we recommend exporting\n  that subset of the data to VCF (via e.g. '--recode vcf'), merging with\n  another tool/script, and then importing the result; PLINK is not yet suited\n  to handling them.\nSee https://www.cog-genomics.org/plink/1.9/data#merge3 for more discussion.\n", nbmarker_ct, (nbmarker_ct == 1)? "" : "s", outname);
  while (0) {
  report_non_biallelics_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  report_non_biallelics_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  report_non_biallelics_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

void merge_alleles_update_str(char* marker_allele_ptr, char** allele_ptrs, uint32_t* distinct_allele_ctp) {
  uint32_t distinct_allele_ct = *distinct_allele_ctp;
  uint32_t allele_idx = 0;
  if ((marker_allele_ptr == g_missing_geno_ptr) || (distinct_allele_ct == 3)) {
    return;
  }
  while (allele_idx < distinct_allele_ct) {
    if (!strcmp(marker_allele_ptr, allele_ptrs[allele_idx])) {
      return;
    }
    allele_idx++;
  }
  *distinct_allele_ctp = distinct_allele_ct + 1;
  if (distinct_allele_ct == 2) {
    return;
  }
  allele_ptrs[distinct_allele_ct] = marker_allele_ptr;
}

uint32_t merge_equal_pos_alleles(char** marker_allele_ptrs, uint32_t marker_uidx, uint32_t* distinct_allele_ctp, char** allele_ptrs) {
  // reverse order so --keep-allele-order works
  merge_alleles_update_str(marker_allele_ptrs[2 * marker_uidx + 1], allele_ptrs, distinct_allele_ctp);
  merge_alleles_update_str(marker_allele_ptrs[2 * marker_uidx], allele_ptrs, distinct_allele_ctp);
  if (*distinct_allele_ctp > 2) {
    return 1;
  }
  return 0;
}

uint32_t save_equal_pos_alleles(const int64_t* ll_buf, uint32_t read_pos, uint32_t read_pos_stop, char** equal_pos_allele_ptrs, char** marker_allele_ptrs) {
  const uint32_t canonical_marker_uidx = (uint32_t)ll_buf[read_pos];
  // (do we really need the reversal in merge_equal_pos_alleles()?)
  // bugfix (26 May 2018): forgot about plink 1.9's allele memory management
  // strategy.
  const char* new_allele = equal_pos_allele_ptrs[1]? equal_pos_allele_ptrs[1] : g_missing_geno_ptr;
  if (allele_reset(new_allele, strlen(new_allele), &(marker_allele_ptrs[canonical_marker_uidx * 2]))) {
    return 1;
  }
  new_allele = equal_pos_allele_ptrs[0]? equal_pos_allele_ptrs[0] : g_missing_geno_ptr;
  if (allele_reset(new_allele, strlen(new_allele), &(marker_allele_ptrs[canonical_marker_uidx * 2 + 1]))) {
    return 1;
  }
  ++read_pos;
  for (; read_pos < read_pos_stop; ++read_pos) {
    const uint32_t cur_marker_uidx = (uint32_t)ll_buf[read_pos];
    // DIRTY HACK (25 May 2018): save (nullptr, canonical_marker_uidx) tuple,
    // and have merge_main() look up the canonical marker_uidx in this case.
    // something like this is necessary for mixed .bed/.ped merge +
    // --merge-equal-pos to work properly.
    // bugfix (26 May 2018): need to free entries here by hand.
    if (marker_allele_ptrs[cur_marker_uidx * 2][1]) {
      free(marker_allele_ptrs[cur_marker_uidx * 2]);
    }
    marker_allele_ptrs[cur_marker_uidx * 2] = nullptr;
    if (marker_allele_ptrs[cur_marker_uidx * 2 + 1][1]) {
      free(marker_allele_ptrs[cur_marker_uidx * 2 + 1]);
    }
    marker_allele_ptrs[cur_marker_uidx * 2 + 1] = (char*)((uintptr_t)canonical_marker_uidx);
  }
  return 0;
}

static inline uint32_t merge_post_msort_update_maps(char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_map, double* marker_cms, double* marker_cms_tmp, uint32_t* pos_buf, int64_t* ll_buf, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t chrom_ct, uint32_t* dedup_marker_ct_ptr, uint32_t merge_equal_pos, char** marker_allele_ptrs, Chrom_info* chrom_info_ptr) {
  // Input: ll_buf is a sequence of sorted arrays (one per chromosome) with
  // base-pair positions in high 32 bits, and pre-sort indices in low 32 bits.
  // Chromosome boundaries are stored in chrom_start[].
  // Pre-sort indices refer to marker ID ASCII ordering, which is how lookup
  // will be performed mid-merge.  There may be duplicate positions, and
  // markers that don't pass the chromosome filter.

  // Result: Duplicates have been collapsed, with chrom_start[] updated.
  // pos_buf contains sorted base-pair positions,
  // post-chromosome-filtering-and-duplicate-removal.
  // marker_map[n] is the post-filtering position in all other arrays of marker
  // ID n by ASCII ordering.
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  uint32_t read_pos = 0;
  uint32_t write_pos = 0; // may be lower than read_pos due to dups
  uint32_t position_warning_ct = 0;
  char* equal_pos_allele_ptrs[2];
  equal_pos_allele_ptrs[0] = nullptr;
  equal_pos_allele_ptrs[1] = nullptr;
  uint32_t equal_pos_allele_ct = 0;
  uint32_t chrom_idx;
  uint32_t chrom_read_end_idx;
  int64_t llxx;
  uint32_t unplaced;
  uint32_t prev_bp;
  uint32_t cur_bp;
  uint32_t presort_idx;
  for (chrom_idx = 0; chrom_idx < chrom_ct; chrom_idx++) {
    unplaced = chrom_id[chrom_idx]; // initially chromosome code
    if (!IS_SET(chrom_mask, unplaced)) {
      read_pos = chrom_start[chrom_idx + 1];
      chrom_start[chrom_idx + 1] = write_pos;
      continue;
    }
    unplaced = (unplaced == 0) || (chrom_info_ptr->zero_extra_chroms && (unplaced > chrom_info_ptr->max_code));
    const uint32_t cur_merge_equal_pos = merge_equal_pos && (!unplaced);
    chrom_read_end_idx = chrom_start[chrom_idx + 1];
    // ll_buf has base-pair positions in high 32 bits, and pre-sort indices in
    // low 32 bits.
    llxx = ll_buf[read_pos++];
    marker_cms[write_pos] = marker_cms_tmp[(uint32_t)llxx];
    prev_bp = (uint32_t)(((uint64_t)llxx) >> 32);
    pos_buf[write_pos] = prev_bp;
    marker_map[(uint32_t)llxx] = write_pos++;
    uint32_t equal_pos_bp = 0xffffffffU;  // force initial mismatch
    uint32_t equal_pos_read_start = 0;
    for (; read_pos < chrom_read_end_idx; read_pos++) {
      llxx = ll_buf[read_pos];
      presort_idx = (uint32_t)llxx;
      cur_bp = (uint32_t)(llxx >> 32);
      // don't care about position conflicts on chr 0 (unplaced).
      if (cur_merge_equal_pos) {
        // bugfix (25 May 2018): previous merge_alleles() usage was broken for
        // multiple reasons
        if (prev_bp == cur_bp) {
          if (equal_pos_bp != prev_bp) {
            equal_pos_read_start = read_pos - 1;
            equal_pos_bp = prev_bp;
            equal_pos_allele_ptrs[0] = nullptr;
            equal_pos_allele_ptrs[1] = nullptr;
            equal_pos_allele_ct = 0;
            if (merge_equal_pos_alleles(marker_allele_ptrs, (uint32_t)ll_buf[equal_pos_read_start], &equal_pos_allele_ct, equal_pos_allele_ptrs)) {
              LOGERRPRINTFWW("Error: --merge-equal-pos failure: inconsistent alleles at %s:%u.\n", "?", cur_bp);
              return 1;
            }
          }
          if (merge_equal_pos_alleles(marker_allele_ptrs, presort_idx, &equal_pos_allele_ct, equal_pos_allele_ptrs)) {
            LOGERRPRINTFWW("Error: --merge-equal-pos failure: inconsistent alleles at %s:%u.\n", "?", cur_bp);
            return 1;
          }
	  marker_map[presort_idx] = write_pos - 1;
	  continue;
        } else {
          if (equal_pos_bp == prev_bp) {
            if (save_equal_pos_alleles(ll_buf, equal_pos_read_start, read_pos, equal_pos_allele_ptrs, marker_allele_ptrs)) {
              return 1;
            }
          }
          prev_bp = cur_bp;
        }
      } else if ((prev_bp == cur_bp) && (!unplaced)) {
        // shouldn't print warning in --merge-equal-pos case
        LOGPREPRINTFWW("Warning: Variants '%s' and '%s' have the same position.\n", &(marker_ids[max_marker_id_len * presort_idx]), &(marker_ids[max_marker_id_len * ((uint32_t)ll_buf[read_pos - 1])]));
        if (position_warning_ct < 3) {
          logerrprintb();
        } else {
          logstr(g_logbuf);
        }
        position_warning_ct++;
      } else {
        prev_bp = cur_bp;
      }
      marker_map[presort_idx] = write_pos;
      marker_cms[write_pos] = marker_cms_tmp[presort_idx];
      pos_buf[write_pos++] = cur_bp;
    }
    if (equal_pos_bp == prev_bp) {
      if (save_equal_pos_alleles(ll_buf, equal_pos_read_start, read_pos, equal_pos_allele_ptrs, marker_allele_ptrs)) {
        return 1;
      }
    }
    read_pos = chrom_start[chrom_idx + 1];
    chrom_start[chrom_idx + 1] = write_pos;
  }
  if (position_warning_ct > 3) {
    fprintf(stderr, "%u more same-position warning%s: see log file.\n", position_warning_ct - 3, (position_warning_ct == 4)? "" : "s");
  }
  *dedup_marker_ct_ptr = write_pos;
  return 0;
}

static inline int32_t merge_must_track_write(int32_t mm) {
  // modes 6 and 7 can be sped up with early pruning of nonoverlapping
  // markers, but not worth complicating code for this.
  return (mm == 1) || (mm > 5) || (mm == 4);
}

static inline int32_t merge_first_mode(int32_t mm, uint32_t merge_equal_pos) {
  if (merge_equal_pos) {
    return (mm > 5)? 4 : mm;
  } else {
    return merge_must_track_write(mm)? ((mm == 1)? 1 : 4) : 5;
  }
}

int32_t merge_diff_print(FILE* outfile, char* idbuf, char* marker_id, char* sample_id, unsigned char newval, unsigned char oldval, char** marker_allele_ptrs) {
  char* bufptr = token_endnn(sample_id);
  uint32_t slen = strlen_se(marker_id);
  const char* ma1p[4];
  const char* ma2p[4];
  char wbuf[8];
  uint32_t slen2;
  memcpyx(idbuf, marker_id, slen, 0);
  fprintf(outfile, "%20s", idbuf);
  slen = (bufptr++) - sample_id;
  memcpyx(idbuf, sample_id, slen, 0);
  fprintf(outfile, " %20s %20s", idbuf, bufptr);
  ma1p[0] = marker_allele_ptrs[0];
  ma1p[1] = g_missing_geno_ptr;
  ma1p[2] = marker_allele_ptrs[0];
  ma1p[3] = marker_allele_ptrs[1];
  ma2p[0] = marker_allele_ptrs[0];
  ma2p[1] = ma1p[1];
  ma2p[2] = marker_allele_ptrs[1];
  ma2p[3] = marker_allele_ptrs[1];
  slen = strlen(ma1p[newval]);
  slen2 = strlen(ma2p[newval]);
  if (slen + slen2 > 6) {
    fprintf(outfile, " %s/%s", ma1p[newval], ma2p[newval]);
  } else {
    memcpy(wbuf, ma1p[newval], slen);
    wbuf[slen] = '/';
    memcpy(&(wbuf[slen + 1]), ma2p[newval], slen2 + 1);
    fprintf(outfile, " %8s", wbuf);
  }
  slen = strlen(ma1p[oldval]);
  slen2 = strlen(ma2p[oldval]);
  if (slen + slen2 > 6) {
    fprintf(outfile, " %s/%s \n", ma1p[oldval], ma2p[oldval]);
  } else {
    memcpy(wbuf, ma1p[oldval], slen);
    wbuf[slen] = '/';
    memcpy(&(wbuf[slen + 1]), ma2p[oldval], slen2 + 1);
    fprintf(outfile, " %8s \n", wbuf);
  }
  return ferror(outfile);
}

int32_t merge_main(char* bedname, char* bimname, char* famname, char* bim_loadbuf, uint32_t max_bim_linelen, uint32_t tot_sample_ct, uint32_t tot_marker_ct, uint32_t dedup_marker_ct, uint32_t start_marker_idx, uint32_t marker_window_size, char** marker_allele_ptrs, char* marker_ids, uintptr_t max_marker_id_len, char* sample_ids, uintptr_t max_sample_id_len, uint32_t merge_nsort, uint32_t* sample_nsmap, uint32_t* flex_map, uint32_t* marker_map, char* idbuf, unsigned char* readbuf, unsigned char* writebuf, uint32_t merge_mode, uintptr_t* markbuf, FILE* outfile, uint64_t* diff_total_overlap_ptr, uint64_t* diff_not_both_genotyped_ptr, uint64_t* diff_discordant_ptr, uint32_t ped_buflen) {
  // flex_map maps samples for binary filesets, and markers for text filesets.
  uint32_t is_binary = famname? 1 : 0;
  FILE* bedfile = nullptr;
  FILE* infile2 = nullptr;
  int32_t retval = 0;
  // bugfix: there was a potential integer overflow back when these were
  // uint32_t
  uintptr_t tot_sample_ct4 = (tot_sample_ct + 3) / 4;
  uintptr_t tot_sample_ctl = BITCT_TO_WORDCT(tot_sample_ct);
  uint32_t end_marker_idx = start_marker_idx + marker_window_size;
  uint32_t marker_in_idx = 0xffffffffU; // overflow to zero on first add
  uint32_t last_marker_in_idx = 0xfffffffeU;
  uint32_t cur_sample_ct = 0;
  uintptr_t* mbufptr = nullptr; // merge mode 1, 4, 6, 7
  uintptr_t* readbuf_w = nullptr; // used for main binary load
  const char* missing_geno_ptr = g_missing_geno_ptr;
  uint64_t diff_total_overlap = 0;
  uint64_t diff_not_both_genotyped = 0;
  uint64_t diff_discordant = 0;
  uint32_t cur_sample_ct4 = 0;
  uint32_t cur_sample_ctl2 = 0;
  uint32_t is_ped_compound = 1; // 0 = no, 1 = unresolved, 2 = yes
  uint32_t alen1 = 1;
  uint32_t alen2 = 1;
  uintptr_t uljj = 0;
  uintptr_t* mbufptr2;
  uintptr_t* rbufptr;
  uint32_t cm_col_exists;
  char* aptr1;
  char* aptr2;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  unsigned char* wbufptr;
  unsigned char* wbufptr2;
  uintptr_t sample_idx;
  uintptr_t sample_idx_cur_max;
  uintptr_t line_idx;
  uintptr_t cur_word;
  uintptr_t ulii;
  uint32_t marker_out_idx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  int32_t ii;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char ucc3;
  unsigned char ucc4;
  char cc;
  if (is_binary) {
    if (fopen_checked(famname, "r", &infile2)) {
      goto merge_main_ret_OPEN_FAIL;
    }
    while (fgets(g_textbuf, MAXLINELEN, infile2)) {
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr2 = token_endnn(bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr4 = token_endnn(bufptr3); // safe since file was validated
      uii = (bufptr2 - bufptr);
      ujj = (bufptr4 - bufptr3);
      memcpyx(memcpyax(idbuf, bufptr, uii, '\t'), bufptr3, ujj, 0);
      if (merge_nsort) {
        ii = bsearch_str_natural(idbuf, sample_ids, max_sample_id_len, tot_sample_ct);
      } else {
	ii = bsearch_str(idbuf, uii + ujj + 1, sample_ids, max_sample_id_len, tot_sample_ct);
	if (sample_nsmap && (ii != -1)) {
	  ii = sample_nsmap[(uint32_t)ii];
	}
      }
      if (ii == -1) {
	// previously validated, so give read failure error code instead of
	// invalid format
	goto merge_main_ret_READ_FAIL;
      }
      flex_map[cur_sample_ct++] = ii;
    }
    if (!feof(infile2)) {
      goto merge_main_ret_READ_FAIL;
    }
    fclose_null(&infile2);
    cur_sample_ct4 = (cur_sample_ct + 3) / 4;
    cur_sample_ctl2 = QUATERCT_TO_WORDCT(cur_sample_ct);
  } else {
    bim_loadbuf = g_textbuf;
    max_bim_linelen = MAXLINELEN;
  }
  if (fopen_checked(bimname, "r", &infile2)) {
    goto merge_main_ret_OPEN_FAIL;
  }
  if (check_cm_col(infile2, bim_loadbuf, is_binary, 1, max_bim_linelen, &cm_col_exists, &ulii)) {
    goto merge_main_ret_READ_FAIL;
  }
  if (!ulii) {
    bim_loadbuf[0] = '\0';
  }
  if (fopen_checked(bedname, is_binary? FOPEN_RB : "r", &bedfile)) {
    goto merge_main_ret_OPEN_FAIL;
  }
  if (is_binary) {
    if (!start_marker_idx) {
      if (fread(readbuf, 1, 3, bedfile) < 3) {
	goto merge_main_ret_READ_FAIL;
      }
      if (memcmp(readbuf, "l\x1b\x01", 3)) {
	if (!memcmp(readbuf, "l\x1b", 3)) {
	  LOGPREPRINTFWW("Error: %s is an sample-major binary file. Convert to variant-major (with e.g. --make-bed) and then reattempt the merge.\n", bedname);
	} else {
	  LOGPREPRINTFWW("Error: %s is not a PLINK 1 binary file.\n", bedname);
	}
	goto merge_main_ret_INVALID_FORMAT_2N;
      }
    }
    readbuf_w = (uintptr_t*)readbuf;
    readbuf_w[cur_sample_ctl2 - 1] = 0;
  }
  do {
    bufptr = skip_initial_spaces(bim_loadbuf);
    if (is_eoln_or_comment_kns(*bufptr)) {
      continue;
    }
    ++marker_in_idx;
    bufptr = next_token(bufptr);
    bufptr2 = next_token_mult(bufptr, 1 + cm_col_exists);
    if (!bufptr2) {
      goto merge_main_ret_READ_FAIL;
    }
    if (*bufptr2 == '-') {
      if (!is_binary) {
	flex_map[marker_in_idx] = 0xffffffffU;
      }
      continue;
    }
    bufptr3 = token_endnn(bufptr);
    ii = bsearch_str(bufptr, (uintptr_t)(bufptr3 - bufptr), marker_ids, max_marker_id_len, tot_marker_ct);
    if (ii == -1) {
      goto merge_main_ret_READ_FAIL;
    }
    marker_out_idx = marker_map[(uint32_t)ii];
    if ((marker_out_idx < start_marker_idx) || (marker_out_idx >= end_marker_idx)) {
      if (!is_binary) {
	flex_map[marker_in_idx] = 0xffffffffU;
      }
      continue;
    }
    if (is_binary) {
      if (marker_in_idx != last_marker_in_idx + 1) {
	if (fseeko(bedfile, 3 + ((uint64_t)marker_in_idx) * cur_sample_ct4, SEEK_SET)) {
	  goto merge_main_ret_READ_FAIL;
	}
      }
      bufptr2 = next_token(bufptr2);
      bufptr3 = next_token(bufptr2);
      if (no_more_tokens_kns(bufptr3)) {
	goto merge_main_ret_READ_FAIL;
      }
      alen1 = strlen_se(bufptr2);
      bufptr2[alen1] = '\0';
      alen2 = strlen_se(bufptr3);
      bufptr3[alen2] = '\0';
      uint32_t canonical_ascii_uidx = (uint32_t)ii;
      bufptr4 = marker_allele_ptrs[canonical_ascii_uidx * 2];
      bufptr5 = marker_allele_ptrs[canonical_ascii_uidx * 2 + 1];
      if (!bufptr4) {
        // --merge-equal-pos non-canonical index.  Retrieve the canonical
        // index.
        canonical_ascii_uidx = (uintptr_t)bufptr5;
        bufptr4 = marker_allele_ptrs[canonical_ascii_uidx * 2];
        bufptr5 = marker_allele_ptrs[canonical_ascii_uidx * 2 + 1];
      }

      last_marker_in_idx = marker_in_idx;
      if (!cur_sample_ct) {
	continue;
      }
      if (load_raw(cur_sample_ct4, bedfile, readbuf_w)) {
	goto merge_main_ret_READ_FAIL;
      }
      if ((((*bufptr2 != '0') || (alen1 != 1)) && (!strcmp(bufptr2, bufptr5)))  || (((*bufptr3 != '0') || (alen2 != 1)) && (!strcmp(bufptr3, bufptr4)))) {
	// Ack, how did this bug not get caught for so long!
	// Necessary to use reverse_loadbuf here to handle last byte properly
	// (since cur_sample_ct % 4 is not necessarily the same as
	// tot_sample_ct % 4).  And while I'm at it, may as well switch
	// the main loops to be word-based.
	reverse_loadbuf(cur_sample_ct, (unsigned char*)readbuf_w);
      }
      rbufptr = readbuf_w;
      wbufptr = &(writebuf[(marker_out_idx - start_marker_idx) * tot_sample_ct4]);
      if (merge_must_track_write(merge_mode)) {
	mbufptr = &(markbuf[(marker_out_idx - start_marker_idx) * tot_sample_ctl]);
      }
      switch (merge_mode) {
      case 1: // difference -> missing
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ucc = cur_word & 3;
	    // bugfix: do NOT set flag, etc. on missing call
	    if (ucc != 1) {
	      ujj = flex_map[sample_idx];
	      ukk = ujj / BITCT;
	      ulii = ONELU << (ujj % BITCT);
	      wbufptr2 = &(wbufptr[ujj / 4]);
	      umm = (ujj % 4) * 2;
	      unn = 3U << umm;
	      if (mbufptr[ukk] & ulii) {
		ucc2 = *wbufptr2;
		if ((ucc2 ^ (ucc << umm)) & unn) {
		  *wbufptr2 = (ucc2 & (~unn)) | (1U << umm);
		}
	      } else {
		mbufptr[ukk] |= ulii;
		*wbufptr2 = ((*wbufptr2) & (~unn)) | (ucc << umm);
	      }
	    }
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      case 2: // only overwrite originally missing
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ujj = flex_map[sample_idx];
	    ukk = (ujj % 4) * 2;
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    ucc2 = *wbufptr2;
	    if (((ucc2 >> ukk) & 3) == 1) {
	      ucc = cur_word & 3;
	      *wbufptr2 = (ucc2 & (~(3U << ukk))) | (ucc << ukk);
	    }
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      case 3: // only overwrite if nonmissing in new file
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ucc = cur_word & 3;
	    if (ucc != 1) {
	      ujj = flex_map[sample_idx];
	      ukk = (ujj % 4) * 2;
	      wbufptr2 = &(wbufptr[ujj / 4]);
	      *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | (ucc << ukk);
	    }
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      case 4: // never overwrite
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ujj = flex_map[sample_idx];
	    ukk = ujj / BITCT;
	    ulii = ONELU << (ujj % BITCT);
	    if (!(mbufptr[ukk] & ulii)) {
	      mbufptr[ukk] |= ulii;
	      wbufptr2 = &(wbufptr[ujj / 4]);
	      ukk = (ujj % 4) * 2;
	      ucc = cur_word & 3;
	      *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | (ucc << ukk);
	    }
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      case 5: // always overwrite
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ucc = cur_word & 3;
	    ujj = flex_map[sample_idx];
	    ukk = (ujj % 4) * 2;
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | (ucc << ukk);
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      case 6: // report all mismatches
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ujj = flex_map[sample_idx];
	    if (mbufptr[ujj / BITCT] & (ONELU << (ujj % BITCT))) {
	      // would prefer to do this by multiplying sample overlap with
	      // marker overlap, but the same-position automerge screws
	      // with that
	      diff_total_overlap++;
	      ukk = (ujj % 4) * 2;
	      ucc = cur_word & 3;
	      ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	      umm = ((ucc == 1) || (ucc3 == 1));
	      if (umm) {
		diff_not_both_genotyped++;
	      }
	      if (ucc != ucc3) {
		if (!umm) {
		  diff_discordant++;
		}
		if (merge_diff_print(outfile, idbuf, bufptr, &(sample_ids[ujj * max_sample_id_len]), ucc, ucc3, &(marker_allele_ptrs[canonical_ascii_uidx * 2]))) {
		  goto merge_main_ret_WRITE_FAIL;
		}
	      }
	    }
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      case 7: // report nonmissing mismatches
	sample_idx = 0;
	do {
	  sample_idx_cur_max = sample_idx + BITCT2;
	  if (sample_idx_cur_max > cur_sample_ct) {
	    sample_idx_cur_max = cur_sample_ct;
	  }
	  cur_word = *rbufptr++;
	  for (; sample_idx < sample_idx_cur_max; sample_idx++) {
	    ujj = flex_map[sample_idx];
	    if (mbufptr[ujj / BITCT] & (ONELU << (ujj % BITCT))) {
	      diff_total_overlap++;
	      ukk = (ujj % 4) * 2;
	      ucc = cur_word & 3;
	      ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	      if ((ucc == 1) || (ucc3 == 1)) {
		diff_not_both_genotyped++;
	      } else if (ucc != ucc3) {
		diff_discordant++;
		if (merge_diff_print(outfile, idbuf, bufptr, &(sample_ids[ujj * max_sample_id_len]), ucc, ucc3, &(marker_allele_ptrs[((uint32_t)ii) * 2]))) {
		  goto merge_main_ret_WRITE_FAIL;
		}
	      }
	    }
	    cur_word >>= 2;
	  }
	} while (sample_idx < cur_sample_ct);
	break;
      }
    } else {
      flex_map[marker_in_idx] = ii;
    }
  } while (fgets(bim_loadbuf, max_bim_linelen, infile2));
  if (!feof(infile2)) {
    goto merge_main_ret_READ_FAIL;
  }
  if (!is_binary) {
    last_marker_in_idx = marker_in_idx + 1; // total count
    line_idx = 0;
    while (fgets((char*)readbuf, ped_buflen, bedfile)) {
      line_idx++;
      bufptr = skip_initial_spaces((char*)readbuf);
      cc = *bufptr;
      if (is_eoln_or_comment_kns(cc)) {
	continue;
      }
      // only possible to get here if sample_ct and marker_ct are positive
      bufptr2 = token_endnn(bufptr);
      uii = (bufptr2 - bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr2 = token_endnn(bufptr3);
      ujj = (bufptr2 - bufptr3);
      memcpyx(memcpyax(idbuf, bufptr, uii, '\t'), bufptr3, ujj, 0);
      if (merge_nsort) {
        ii = bsearch_str_natural(idbuf, sample_ids, max_sample_id_len, tot_sample_ct);
      } else {
	ii = bsearch_str(idbuf, uii + ujj + 1, sample_ids, max_sample_id_len, tot_sample_ct);
	if (sample_nsmap && (ii != -1)) {
	  ii = sample_nsmap[(uint32_t)ii];
	}
      }
      if (ii == -1) {
	goto merge_main_ret_READ_FAIL;
      }
      bufptr3 = next_token_mult(skip_initial_spaces(bufptr2), 4);
      if (!bufptr3) {
	goto merge_main_ret_READ_FAIL;
      }
      if (is_ped_compound == 1) {
	bufptr4 = bufptr3;
	for (marker_in_idx = 0; marker_in_idx < 2 * last_marker_in_idx; marker_in_idx++) {
          if (is_eoln_kns(*bufptr4)) {
	    is_ped_compound = 2;
	    break;
	  }
	  bufptr4 = skip_initial_spaces(token_endnn(bufptr4));
	}
	if (is_ped_compound == 1) {
	  is_ped_compound = 0;
	}
      }
      wbufptr = &(writebuf[((uint32_t)ii) / 4]);
      ujj = (((uint32_t)ii) % 4) * 2;
      ucc = ~(3U << ujj);
      ucc4 = (1U << ujj);
      if (merge_must_track_write(merge_mode)) {
	mbufptr = &(markbuf[((uint32_t)ii) / BITCT]);
	uljj = ONELU << (((uint32_t)ii) % BITCT);
      }
      for (marker_in_idx = 0; marker_in_idx < last_marker_in_idx; marker_in_idx++) {
	if (!is_ped_compound) {
	  if (is_eoln_kns(*bufptr3)) {
	    goto merge_main_ret_MISSING_TOKENS;
	  }
	  aptr1 = bufptr3;
	  bufptr3 = token_endnn(bufptr3);
	  alen1 = (uintptr_t)(bufptr3 - aptr1);
	  bufptr3 = skip_initial_spaces(bufptr3);
	  aptr1[alen1] = '\0';
	  if (is_eoln_kns(*bufptr3)) {
	    goto merge_main_ret_MISSING_TOKENS;
	  }
	  aptr2 = bufptr3;
	  bufptr3 = token_endnn(bufptr3);
	  alen2 = (uintptr_t)(bufptr3 - aptr2);
	  bufptr3 = skip_initial_spaces(bufptr3);
	  aptr2[alen2] = '\0';
	} else {
	  cc = *bufptr3;
	  if (is_eoln_kns(cc)) {
	    goto merge_main_ret_MISSING_TOKENS;
	  }
          aptr1 = (char*)(&(g_one_char_strs[((unsigned char)cc) * 2]));
          bufptr3 = skip_initial_spaces(&(bufptr3[1]));
          cc = *bufptr3;
          if (is_eoln_kns(cc)) {
	    goto merge_main_ret_MISSING_TOKENS;
	  }
	  aptr2 = (char*)(&(g_one_char_strs[((unsigned char)cc) * 2]));
	  bufptr3 = skip_initial_spaces(&(bufptr3[1]));
	}

	// lexicographic position (or 0xffffffffU skip indicator)
	uii = flex_map[marker_in_idx];
	if (uii == 0xffffffffU) {
	  continue;
	}
	// actual relative position
	marker_out_idx = marker_map[uii] - start_marker_idx;
        // bugfix (25 May 2018): in --merge-equal-pos case, we need to update a
        // single canonical marker_allele_ptrs[].  We now handle this by
        // storing nullptr in [0] and the canonical index in [1] when the
        // current index is non-canonical.
        char** canonical_marker_allele_ptrs = &(marker_allele_ptrs[uii * 2]);
        if (!canonical_marker_allele_ptrs[0]) {
          const uint32_t canonical_ascii_uidx = (uintptr_t)canonical_marker_allele_ptrs[1];
          canonical_marker_allele_ptrs = &(marker_allele_ptrs[canonical_ascii_uidx * 2]);
        }

	if ((*aptr1 == '0') && (alen1 == 1)) {
          if ((*aptr2 != '0') || (alen2 != 1)) {
	    goto merge_main_ret_HALF_MISSING;
	  }
	  ucc2 = 1; // final PLINK encoding
	} else if ((*aptr2 == '0') && (alen2 == 1)) {
	  goto merge_main_ret_HALF_MISSING;
	} else {
	  ucc2 = 0; // A2 count
	  if (!strcmp(aptr1, canonical_marker_allele_ptrs[1])) {
	    ucc2++;
	  } else if (strcmp(aptr1, canonical_marker_allele_ptrs[0])) {
	    // new allele code
	    if (canonical_marker_allele_ptrs[1] == missing_geno_ptr) {
	      // fill A2 first
	      ucc2++;
	      ukk = 1;
	    } else if (canonical_marker_allele_ptrs[0] == missing_geno_ptr) {
	      ukk = 0;
	    } else {
	      goto merge_main_ret_NOT_BIALLELIC;
	    }
	    if (allele_set(aptr1, alen1, &(canonical_marker_allele_ptrs[ukk]))) {
	      goto merge_main_ret_NOMEM;
	    }
	  }
	  if (!strcmp(aptr2, canonical_marker_allele_ptrs[1])) {
	    ucc2++;
	  } else if (strcmp(aptr2, canonical_marker_allele_ptrs[0])) {
	    // put A2 check second since the only way A2 will be unset when A1
	    // is set at this point is if it was specified that way in an
	    // earlier binary file.
	    if (canonical_marker_allele_ptrs[0] == missing_geno_ptr) {
	      ukk = 0;
	    } else if (canonical_marker_allele_ptrs[1] == missing_geno_ptr) {
              ukk = 1;
              // bugfix (14 Nov 2017): forgot to increment the A2 allele count!
              ++ucc2;
	    } else {
	      goto merge_main_ret_NOT_BIALLELIC;
	    }
	    if (allele_set(aptr2, alen2, &(canonical_marker_allele_ptrs[ukk]))) {
	      goto merge_main_ret_NOMEM;
	    }
	  }
	  // transform to final PLINK encoding
	  if (ucc2) {
	    ucc2++;
	  }
	}
	switch (merge_mode) {
	case 1:
	  mbufptr2 = &(mbufptr[marker_out_idx * tot_sample_ctl]);
	  wbufptr2 = &(wbufptr[marker_out_idx * tot_sample_ct4]);
	  if (ucc2 != 1) {
	    if ((*mbufptr2) & uljj) {
	      ucc3 = *wbufptr2;
	      if (((ucc3 >> ujj) & 3) != ucc2) {
		*wbufptr2 = (ucc3 & ucc) | ucc4;
	      }
	    } else {
	      *mbufptr2 |= uljj;
	      *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	    }
	  }
	  break;
	case 2:
	  wbufptr2 = &(wbufptr[marker_out_idx * tot_sample_ct4]);
	  ucc3 = *wbufptr2;
	  if (((ucc3 >> ujj) & 3) == 1) {
	    *wbufptr2 = (ucc3 & ucc) | (ucc2 << ujj);
	  }
	  break;
	case 3:
	  if (ucc2 == 1) {
	    break;
	  }
	  // fall through
	case 5:
	  wbufptr2 = &(wbufptr[marker_out_idx * tot_sample_ct4]);
	  *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  break;
	case 4:
	  mbufptr2 = &(mbufptr[marker_out_idx * tot_sample_ctl]);
	  if (!((*mbufptr2) & uljj)) {
	    *mbufptr2 |= uljj;
	    wbufptr2 = &(wbufptr[marker_out_idx * tot_sample_ct4]);
	    *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  }
	  break;
	case 6:
	  if (mbufptr[marker_out_idx * tot_sample_ctl] & uljj) {
	    diff_total_overlap++;
	    ucc3 = ((wbufptr[marker_out_idx * tot_sample_ct4] >> ujj) & 3);
	    umm = ((ucc2 == 1) || (ucc3 == 1));
	    if (umm) {
	      diff_not_both_genotyped++;
	    }
	    if (ucc2 != ucc3) {
	      if (!umm) {
		diff_discordant++;
	      }
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(sample_ids[((uint32_t)ii) * max_sample_id_len]), ucc2, ucc3, canonical_marker_allele_ptrs)) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	  break;
	case 7:
	  if (mbufptr[marker_out_idx * tot_sample_ctl] & uljj) {
	    diff_total_overlap++;
	    ucc3 = ((wbufptr[marker_out_idx * tot_sample_ct4] >> ujj) & 3);
	    if ((ucc2 == 1) || (ucc3 == 1)) {
	      diff_not_both_genotyped++;
	    } else if (ucc2 != ucc3) {
	      diff_discordant++;
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(sample_ids[((uint32_t)ii) * max_sample_id_len]), ucc2, ucc3, canonical_marker_allele_ptrs)) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      }
      if (!is_eoln_kns(*bufptr3)) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has more tokens than expected.\n", line_idx, bedname);
	goto merge_main_ret_INVALID_FORMAT_2N;
      }
    }
    if (!feof(bedfile)) {
      goto merge_main_ret_READ_FAIL;
    }
  }
  if (merge_mode > 5) {
    *diff_total_overlap_ptr += diff_total_overlap;
    *diff_not_both_genotyped_ptr += diff_not_both_genotyped;
    *diff_discordant_ptr += diff_discordant;
  }
  while (0) {
  merge_main_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_main_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_main_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  merge_main_ret_NOT_BIALLELIC:
    logprint("\n");
    LOGERRPRINTFWW("Error: Variant '%s' is not biallelic. To obtain a full list of merge failures, convert your data to binary format and retry the merge.\n", &(marker_ids[uii * max_marker_id_len]));
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_HALF_MISSING:
    logprint("\n");
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has a half-missing call.\n", line_idx, bedname);
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_MISSING_TOKENS:
    LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, bedname);
  merge_main_ret_INVALID_FORMAT_2N:
    logprint("\n");
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(bedfile);
  fclose_cond(infile2);
  return retval;
}

int32_t merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, char* sample_sort_fname, uint64_t calculation_type, uint32_t merge_type, uint32_t sample_sort, uint64_t misc_flags, Chrom_info* chrom_info_ptr) {
  FILE* mergelistfile = nullptr;
  FILE* outfile = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t max_sample_id_len = 0;
  uintptr_t max_marker_id_len = 0;
  uint32_t max_sample_full_len = 0;
  uint32_t keep_allele_order = (misc_flags / MISC_KEEP_ALLELE_ORDER) & 1;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t is_dichot_pheno = 1;
  uint32_t merge_list = merge_type & MERGE_LIST;
  uint32_t merge_mode = merge_type & MERGE_MODE_MASK;
  uint32_t merge_nsort = ((!sample_sort) || (sample_sort == SAMPLE_SORT_NATURAL))? 1 : 0;
  uint32_t merge_equal_pos = (merge_type / MERGE_EQUAL_POS) & 1;
  uint32_t allow_no_samples = (misc_flags / MISC_ALLOW_NO_SAMPLES) & 1;
  uint32_t allow_no_variants = (misc_flags / MISC_ALLOW_NO_VARS) & 1;
  Ll_str* non_biallelics = nullptr;
  uint32_t ped_buflen = MAXLINELEN;
  uint32_t max_bim_linelen = 0;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char* pheno_c_char = nullptr;
  double* pheno_d = nullptr;
  uint32_t* sample_nsmap = nullptr;
  uint32_t max_cur_sample_ct = 0;
  uint32_t max_cur_marker_text_ct = 0;
  uintptr_t* markbuf = nullptr; // needed for merge modes 1, 4, 6, 7
  uint64_t diff_total_overlap = 0;
  uint64_t diff_not_both_genotyped = 0;
  uint64_t diff_discordant = 0;
  uint64_t position_warning_ct = 0;
  uint32_t orig_idx = 0;
  uint32_t cur_marker_ct = 0;
  uint32_t tot_marker_ct = 0;
  int32_t retval = 0;
  uint32_t* map_reverse = nullptr;
  uintptr_t* reversed = nullptr;
  char* bim_loadbuf = nullptr;
  // N.B. marker_allele_ptrs are ordered by marker_id instead of position
  char** marker_allele_ptrs = nullptr;
  Ll_fam** htable_fam;
  Ll_bim** htable_bim;
  uintptr_t* pcptr;
  uintptr_t markers_per_pass;
  uint32_t pass_ct;
  char* sample_ids;
  char* sample_fids;
  char* marker_ids;
  uint32_t* marker_map;
  uint32_t* flex_map;
  double* marker_cms;
  double* marker_cms_tmp;
  uint32_t* pos_buf;
  int64_t* ll_buf;
  uintptr_t mlpos;
  uintptr_t merge_ct;
  char* idbuf;
  char* mergelist_buf;
  char** mergelist_bed;
  char** mergelist_bim;
  char** mergelist_fam;
  uint32_t cur_sample_ct;
  uint32_t tot_sample_ct;
  uint32_t tot_sample_ct4;
  uint32_t dedup_marker_ct;
  uint64_t ullxx;
  int64_t llxx;
  uintptr_t line_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  Ll_fam* llfam_ptr;
  Ll_bim* llbim_ptr;
  uint32_t* chrom_start;
  uint32_t* chrom_id;
  uint32_t chrom_ct;
  unsigned char* readbuf;
  unsigned char* writebuf;
  unsigned char* ubufptr;
  char cc;
  unsigned char ucc;
  if (bigstack_alloc_ui(MAX_POSSIBLE_CHROM + 1, &chrom_start) ||
      bigstack_alloc_ui(MAX_POSSIBLE_CHROM, &chrom_id)) {
    goto merge_datasets_ret_NOMEM;
  }

  if (!merge_mode) {
    merge_mode = 1;
  }
  if (merge_list) {
    if (fopen_checked(mergename1, "r", &mergelistfile)) {
      goto merge_datasets_ret_READ_FAIL;
    }
    merge_ct = (famname[0] != '\0');
    ullxx = 0;
    // first pass: determine merge_ct, mergelist_buf size, verify no lines have
    // > 3 entries
    g_textbuf[MAXLINELEN - 1] = ' ';
    line_idx = 0;
    while (fgets(g_textbuf, MAXLINELEN, mergelistfile)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --merge-list file is pathologically long.\n", line_idx);
	goto merge_datasets_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(g_textbuf);
      if (no_more_tokens_kns(bufptr)) {
	continue;
      }
      bufptr2 = next_token_mult(bufptr, 3);
      if (!no_more_tokens_kns(bufptr2)) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --merge-list file has more tokens than expected.\n", line_idx);
        goto merge_datasets_ret_INVALID_FORMAT_2;
      }
      if (no_more_tokens_kns(next_token(bufptr))) {
	bufptr2 = token_endnn(bufptr);
	ulii = bufptr2 - bufptr;
	if (ulii > FNAMESIZE - 5) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --merge-list file has an excessively long fileset\nprefix.\n", line_idx);
	  goto merge_datasets_ret_INVALID_FORMAT_2;
	}
	ullxx += 3 * ulii + 15;
      } else {
	do {
	  bufptr2 = token_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  if (ulii > FNAMESIZE - 1) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --merge-list file has an excessively long filename.\n", line_idx);
	    goto merge_datasets_ret_INVALID_FORMAT_2;
	  }
	  ullxx += ulii + 1;
	  bufptr = skip_initial_spaces(bufptr2);
	} while (!no_more_tokens_kns(bufptr));
      }
      merge_ct++;
    }
    if (!feof(mergelistfile)) {
      goto merge_datasets_ret_READ_FAIL;
    }
    if (!merge_ct) {
      logerrprint("Error: --merge-list file is empty, and no other input fileset was specified.\n");
      goto merge_datasets_ret_INVALID_FORMAT;
    } else if (merge_ct == 1) {
      if (famname[0] == '\0') {
        logerrprint("Warning: --merge-list file contains only one entry.\n");
      } else {
        logerrprint("Warning: --merge-list file is empty.\n");
      }
    }
#ifndef __LP64__
    if (ullxx > 0x7fffffff) {
      goto merge_datasets_ret_NOMEM;
    }
#endif
    mergelist_bed = (char**)bigstack_alloc(merge_ct * sizeof(intptr_t));
    mergelist_bim = (char**)bigstack_alloc(merge_ct * sizeof(intptr_t));
    mergelist_fam = (char**)bigstack_alloc(merge_ct * sizeof(intptr_t));
    if (bigstack_alloc_c((uintptr_t)ullxx, &mergelist_buf)) {
      goto merge_datasets_ret_NOMEM;
    }
    rewind(mergelistfile);
    bufptr4 = mergelist_buf;
    mlpos = (famname[0] != '\0');
    while (fgets(g_textbuf, MAXLINELEN, mergelistfile)) {
      bufptr = skip_initial_spaces(g_textbuf);
      if (no_more_tokens_kns(bufptr)) {
	continue;
      }
      bufptr2 = token_endnn(bufptr);
      ulii = (bufptr2 - bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      if (no_more_tokens_kns(bufptr3)) {
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = memcpya(memcpya(bufptr4, bufptr, ulii), ".bed", 5);
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = memcpya(memcpya(bufptr4, bufptr, ulii), ".bim", 5);
	mergelist_fam[mlpos] = bufptr4;
	bufptr4 = memcpya(memcpya(bufptr4, bufptr, ulii), ".fam", 5);
      } else {
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = memcpyax(bufptr4, bufptr, ulii, 0);
	bufptr2 = token_endnn(bufptr3);
	ulii = bufptr2 - bufptr3;
	bufptr = skip_initial_spaces(bufptr2);
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = memcpyax(bufptr4, bufptr3, ulii, 0);
	if (no_more_tokens_kns(bufptr)) {
	  mergelist_fam[mlpos] = nullptr;
	} else {
	  bufptr2 = token_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  mergelist_fam[mlpos] = bufptr4;
	  bufptr4 = memcpyax(bufptr4, bufptr, ulii, 0);
	}
      }
      if (++mlpos == merge_ct) {
	break;
      }
    }
    if ((mlpos < merge_ct) && (!feof(mergelistfile))) {
      goto merge_datasets_ret_READ_FAIL;
    }
    fclose_null(&mergelistfile);
  } else {
    merge_ct = 2;
    mergelist_bed = (char**)bigstack_alloc(2 * sizeof(intptr_t));
    mergelist_bim = (char**)bigstack_alloc(2 * sizeof(intptr_t));
    mergelist_fam = (char**)bigstack_alloc(2 * sizeof(intptr_t));
    mergelist_bed[1] = mergename1;
    mergelist_bim[1] = mergename2;
    mergelist_fam[1] = (merge_type & MERGE_BINARY)? mergename3 : nullptr;
  }
  if (famname[0]) {
    mergelist_bed[0] = bedname;
    mergelist_bim[0] = bimname;
    mergelist_fam[0] = famname;
  }

  // ID counting/duplicate detection strategy:
  // - We do NOT want to scan through .ped files any more times than absolutely
  // necessary.  So we actually use *gasp* a hash table here.
  // - The hash table is positioned at the FAR end of bigstack, automatically
  // sized to ~4MB (or ~2MB on 32-bit systems).  IDs are then stored
  // backwards from there.  This simplifies copying into a sorted list.
  htable_fam = (Ll_fam**)bigstack_end_alloc(HASHSIZE_S * sizeof(intptr_t));
  if (!htable_fam) {
    goto merge_datasets_ret_NOMEM;
  }
  for (uii = 0; uii < HASHSIZE_S; uii++) {
    htable_fam[uii] = nullptr;
  }

  ullxx = 0;
  mlpos = 0;
  for (mlpos = 0; mlpos < merge_ct; mlpos++) {
    retval = merge_fam_id_scan(mergelist_bed[mlpos], mergelist_fam[mlpos], allow_no_samples, &max_sample_id_len, &max_sample_full_len, &is_dichot_pheno, htable_fam, &ullxx, &ped_buflen, &cur_sample_ct, &orig_idx);
    if (retval) {
      goto merge_datasets_ret_1;
    }
    if ((!merge_list) && mlpos) {
      LOGPRINTFWW("%u %s loaded from %s.\n", max_cur_sample_ct, species_str(max_cur_sample_ct), mergelist_fam[0]);
      LOGPRINTFWW("%u %s to be merged from %s.\n", cur_sample_ct, species_str(cur_sample_ct), (merge_type & MERGE_BINARY)? mergelist_fam[1] : mergelist_bed[1]);
      uii = ullxx - max_cur_sample_ct;
      LOGPRINTF("Of these, %u %s new, while %u %s present in the base dataset.\n", uii, (uii == 1)? "is" : "are", cur_sample_ct - uii, (cur_sample_ct - uii == 1)? "is" : "are");
    }
    if (cur_sample_ct > max_cur_sample_ct) {
      max_cur_sample_ct = cur_sample_ct;
    }
  }
#ifdef __LP64__
  if (ullxx > 0x7fffffff) {
    sprintf(g_logbuf, "Error: Too many %s (max 2147483647).\n", g_species_plural);
    goto merge_datasets_ret_INVALID_FORMAT_2;
  }
#else
  // avoid integer overflow in bigstack_alloc calls
  if (ullxx * max_sample_full_len > 0x7fffffff) {
    sprintf(g_logbuf, "Error: Too many %s for 32-bit " PROG_NAME_CAPS ".\n", g_species_plural);
    goto merge_datasets_ret_INVALID_FORMAT_2;
  }
#endif
  if (max_sample_id_len > 2 * MAX_ID_BLEN) {
    logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  tot_sample_ct = ullxx;
  if (sample_sort & (SAMPLE_SORT_NONE | SAMPLE_SORT_FILE)) {
    if (bigstack_alloc_ui(tot_sample_ct, &sample_nsmap)) {
      goto merge_datasets_ret_NOMEM;
    }
  }
  if (bigstack_alloc_c(max_sample_id_len * tot_sample_ct, &sample_ids) ||
      bigstack_alloc_c(max_sample_full_len * tot_sample_ct, &sample_fids)) {
    goto merge_datasets_ret_NOMEM;
  }
  if (is_dichot_pheno) {
    if (bigstack_alloc_c(tot_sample_ct, &pheno_c_char)) {
      goto merge_datasets_ret_NOMEM;
    }
  } else {
    if (bigstack_alloc_d(tot_sample_ct, &pheno_d)) {
      goto merge_datasets_ret_NOMEM;
    }
  }
  if (sample_sort & (SAMPLE_SORT_NONE | SAMPLE_SORT_FILE)) {
    if (bigstack_alloc_ui(tot_sample_ct, &map_reverse)) {
      goto merge_datasets_ret_NOMEM;
    }
  }
  if (sample_sort == SAMPLE_SORT_NONE) {
    for (uii = 0; uii < HASHSIZE_S; uii++) {
      if (htable_fam[uii]) {
	llfam_ptr = htable_fam[uii];
	do {
	  ujj = llfam_ptr->orig_order;
	  strcpy(&(sample_fids[ujj * max_sample_full_len]), llfam_ptr->idstr);
	  if (is_dichot_pheno) {
	    if (llfam_ptr->pheno == -9) {
	      pheno_c_char[ujj] = -1;
	    } else {
	      pheno_c_char[ujj] = llfam_ptr->pheno - 1;
	    }
	  } else {
	    pheno_d[ujj] = llfam_ptr->pheno;
	  }
	  llfam_ptr = llfam_ptr->next;
	} while (llfam_ptr);
      }
    }
    for (ulii = 0; ulii < tot_sample_ct; ulii++) {
      sample_nsmap[ulii] = ulii;
      bufptr = (char*)memchr(&(sample_fids[ulii * max_sample_full_len]), '\t', max_sample_id_len);
      bufptr = (char*)memchr(&(bufptr[1]), '\t', max_sample_id_len);
      *bufptr = '\0';
    }
    if (qsort_ext(sample_fids, tot_sample_ct, max_sample_full_len, strcmp_deref, (char*)sample_nsmap, sizeof(int32_t))) {
      goto merge_datasets_ret_NOMEM;
    }
  } else {
    ulii = 0;
    bufptr = sample_fids;
    for (uii = 0; uii < HASHSIZE_S; uii++) {
      if (htable_fam[uii]) {
	llfam_ptr = htable_fam[uii];
	do {
	  strcpy(bufptr, llfam_ptr->idstr);
	  bufptr = &(bufptr[max_sample_full_len]);
	  if (is_dichot_pheno) {
	    if (llfam_ptr->pheno == -9) {
	      pheno_c_char[ulii] = -1;
	    } else {
	      pheno_c_char[ulii] = llfam_ptr->pheno - 1;
	    }
	  } else {
	    pheno_d[ulii] = llfam_ptr->pheno;
	  }
	  ulii++;
	  llfam_ptr = llfam_ptr->next;
	} while (llfam_ptr);
      }
    }
    // bugfix: parental IDs and phenotype were being used to break sorting
    // ties here, so e.g. "CEU NA07000 0 0 2" was being natural-sorted after
    // "ceu NA07000 0 0 1", resulting in a crash when "CEU NA07000" was looked
    // up later.
    for (ulii = 0; ulii < tot_sample_ct; ulii++) {
      bufptr = (char*)memchr(&(sample_fids[ulii * max_sample_full_len]), '\t', max_sample_id_len);
      bufptr = (char*)memchr(&(bufptr[1]), '\t', max_sample_id_len);
      *bufptr = '\0';
    }
    if (is_dichot_pheno) {
      if (qsort_ext(sample_fids, tot_sample_ct, max_sample_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, pheno_c_char, 1)) {
        goto merge_datasets_ret_NOMEM;
      }
    } else {
      if (qsort_ext(sample_fids, tot_sample_ct, max_sample_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, (char*)pheno_d, sizeof(double))) {
        goto merge_datasets_ret_NOMEM;
      }
    }
    if (sample_sort == SAMPLE_SORT_FILE) {
      retval = merge_sample_sortf(sample_sort_fname, sample_fids, tot_sample_ct, max_sample_full_len, sample_ids, max_sample_id_len, map_reverse);
      if (retval) {
        goto merge_datasets_ret_1;
      }
    }
  }
  bigstack_end_reset(bigstack_end_mark); // deallocate first hash table
  if (merge_mode < 6) {
    memcpy(outname_end, ".fam", 5);
    if (fopen_checked(outname, "w", &outfile)) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
  }
  if (sample_sort == SAMPLE_SORT_NONE) {
    for (ulii = 0; ulii < tot_sample_ct; ulii++) {
      map_reverse[sample_nsmap[ulii]] = ulii;
    }
    for (ulii = 0; ulii < tot_sample_ct; ulii++) {
      ujj = map_reverse[ulii];
      bufptr = &(sample_fids[ujj * max_sample_full_len]);
      bufptr3 = &(sample_ids[ujj * max_sample_id_len]);
      uii = strlen(bufptr) + 1;
      memcpy(bufptr3, bufptr, uii);
      if (merge_mode < 6) {
        // no longer matters whether this is \t or \0, we're about to free it
	bufptr[uii - 1] = '\t';
	uii += strlen(&(bufptr[uii]));
	if (fwrite_checked(bufptr, uii, outfile)) {
	  goto merge_datasets_ret_WRITE_FAIL;
	}
	if (is_dichot_pheno) {
	  cc = pheno_c_char[ulii];
	  fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1");
	} else {
	  fprintf(outfile, "\t%g\n", pheno_d[ulii]);
	}
      }
    }
  } else if (sample_sort != SAMPLE_SORT_FILE) {
    bufptr = sample_fids;
    bufptr3 = sample_ids;
    for (ulii = 0; ulii < tot_sample_ct; ulii++) {
      uii = strlen(bufptr) + 1;
      memcpy(bufptr3, bufptr, uii);
      bufptr3 = &(bufptr3[max_sample_id_len]);
      if (merge_mode < 6) {
	bufptr[uii - 1] = '\t';
	uii += strlen(&(bufptr[uii]));
	if (fwrite_checked(bufptr, uii, outfile)) {
	  goto merge_datasets_ret_WRITE_FAIL;
	}
	if (is_dichot_pheno) {
	  cc = pheno_c_char[ulii];
	  fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1");
	} else {
	  fprintf(outfile, "\t%g\n", pheno_d[ulii]);
	}
      }
      bufptr = &(bufptr[max_sample_full_len]);
    }
  } else {
    // sample_ids already populated
    bufptr = sample_fids;
    if (merge_mode < 6) {
      for (ulii = 0; ulii < tot_sample_ct; ulii++) {
	ujj = map_reverse[ulii];
	bufptr = &(sample_fids[ujj * max_sample_full_len]);
	uii = strlen(bufptr);
	bufptr[uii] = '\t';
	if (fwrite_checked(bufptr, uii + strlen(&(bufptr[uii])), outfile)) {
	  goto merge_datasets_ret_WRITE_FAIL;
	}
	if (is_dichot_pheno) {
          // bugfix (21 Feb 2018): need to use the map_reverse index here,
          // since the phenotypes correspond to ASCII-sorted sample ID order.
	  cc = pheno_c_char[ujj];
	  fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1");
	} else {
	  fprintf(outfile, "\t%g\n", pheno_d[ujj]);
	}
      }
    }
    for (ulii = 0; ulii < tot_sample_ct; ulii++) {
      sample_nsmap[map_reverse[ulii]] = ulii;
    }
  }
  if (merge_mode < 6) {
    if (ferror(outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
  }
  bigstack_reset(sample_fids);
  htable_bim = (Ll_bim**)bigstack_end_alloc(HASHSIZE * sizeof(intptr_t));
  if (!htable_bim) {
    goto merge_datasets_ret_NOMEM;
  }
  for (uii = 0; uii < HASHSIZE; uii++) {
    htable_bim[uii] = nullptr;
  }

  ullxx = 0;
  for (mlpos = 0; mlpos < merge_ct; ++mlpos) {
    retval = merge_bim_scan(mergelist_bim[mlpos], (mergelist_fam[mlpos])? 1 : 0, allow_no_variants, &max_marker_id_len, htable_bim, &max_bim_linelen, &ullxx, &cur_marker_ct, &position_warning_ct, &non_biallelics, allow_extra_chroms, chrom_info_ptr);
    if (retval) {
      goto merge_datasets_ret_1;
    }
    if ((!mlpos) && (ullxx != cur_marker_ct)) {
      // Update (2 May 2020): PLINK 1.07 errored out if the first input fileset
      // had two variants with the same ID.  However, it did *not* do so if
      // this was true of later filesets, so in cases like
      //   https://github.com/chrchang/plink-ng/issues/140
      // where one but not all filesets had a duplicate ID, it would behave in
      // an asymmetric manner.
      // There are valid reasons for permitting duplicate IDs in the first
      // fileset (e.g. there are redundant loci for quality control purposes),
      // so we don't want to copy PLINK 1.07's error-out behavior.  However,
      // there are also common dangers (e.g. there are a whole bunch of
      // variants with ID=. which should be assigned distinct IDs before
      // merge), so printing a warning where there previously was an error is
      // justified.
      // (Obvious todo for PLINK 2.0: also print this warning if the first
      // fileset doesn't have a duplicate ID, but a later fileset does.)
      logerrprint("Warning: First fileset to be merged contains duplicate variant ID(s).  Variants\nwith matching IDs are all merged together; if this is not what you want (e.g.\nyou have a bunch of novel variants, all with ID \".\"), assign distinct IDs to\nthem (with e.g. --set-missing-var-ids) before rerunning this merge.\n");
    }
    if (!merge_list) {
      if (!mlpos) {
	uii = ullxx;
      } else {
	LOGPRINTFWW("%u marker%s loaded from %s.\n", uii, (uii == 1)? "" : "s", mergelist_bim[0]);
	LOGPRINTFWW("%u marker%s to be merged from %s.\n", cur_marker_ct, (cur_marker_ct == 1)? "" : "s", mergelist_bim[1]);
	// bugfix: don't underflow when a single file has duplicate IDs (e.g.
	// '.').
	uii = ullxx - uii;
	LOGPRINTF("Of these, %u %s new, while %u %s present in the base dataset.\n", uii, (uii == 1)? "is" : "are", cur_marker_ct - uii, (cur_marker_ct - uii == 1)? "is" : "are");
      }
    }
    if (!mergelist_fam[mlpos]) {
      if (cur_marker_ct > max_cur_marker_text_ct) {
        max_cur_marker_text_ct = cur_marker_ct;
      }
    }
  }
  if (position_warning_ct > 3) {
    fprintf(stderr, "%" PRIu64 " more multiple-position warning%s: see log file.\n", position_warning_ct - 3, (position_warning_ct == 4)? "" : "s");
  }
#ifdef __LP64__
  if (ullxx > 0x7ffffffd) {
    logerrprint("Error: Too many variants (max 2147483645).\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#else
  if (ullxx * MAXV(max_marker_id_len, 8) > 0x7fffffff) {
    logerrprint("Error: Too many variants for 32-bit " PROG_NAME_CAPS ".\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#endif
  if (max_marker_id_len > MAX_ID_BLEN) {
    logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  } else if (max_marker_id_len > 80) {
    logerrprint("Warning: Unusually long variant ID(s) present.  PLINK 1.9 does not scale well\nto length-80+ variant IDs; consider using a different naming scheme for long\nindels and the like.\n");
  }
  if (non_biallelics) {
    bigstack_reset(bigstack_mark);
    retval = report_non_biallelics(outname, outname_end, non_biallelics);
    if (retval) {
      goto merge_datasets_ret_1;
    }
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  tot_marker_ct = ullxx;
  marker_allele_ptrs = (char**)bigstack_alloc(tot_marker_ct * 2 * sizeof(intptr_t));
  if (!marker_allele_ptrs) {
    goto merge_datasets_ret_NOMEM;
  }
  // prevent cleanup from failing
  for (uii = 0; uii < tot_marker_ct * 2; uii++) {
    marker_allele_ptrs[uii] = missing_geno_ptr;
  }
  if (max_bim_linelen) {
    max_bim_linelen++;
    if (bigstack_alloc_c(max_bim_linelen, &bim_loadbuf)) {
      goto merge_datasets_ret_NOMEM;
    }
  }
  if (bigstack_alloc_c(max_marker_id_len * tot_marker_ct, &marker_ids) ||
      bigstack_alloc_ui(tot_marker_ct, &marker_map) ||
      bigstack_alloc_d(tot_marker_ct, &marker_cms) ||
      bigstack_alloc_ui(tot_marker_ct, &pos_buf) ||
      bigstack_alloc_d(tot_marker_ct, &marker_cms_tmp) ||
      bigstack_alloc_ll(tot_marker_ct, &ll_buf)) {
    goto merge_datasets_ret_NOMEM;
  }
  for (uii = 0; uii < tot_marker_ct; uii++) {
    pos_buf[uii] = uii;
  }
  ulii = 0;
  for (uii = 0; uii < HASHSIZE; uii++) {
    if (htable_bim[uii]) {
      llbim_ptr = htable_bim[uii];
      do {
	strcpy(&(marker_ids[ulii * max_marker_id_len]), llbim_ptr->idstr);
        ulii++;
	llbim_ptr = llbim_ptr->next;
      } while (llbim_ptr);
    }
  }
  // todo: reimplement this in a manner that never performs a variant ID sort.
  // chrom/pos-based sort is of course still needed, but that involves cheaper
  // int64 comparisons.
  if (qsort_ext(marker_ids, tot_marker_ct, max_marker_id_len, strcmp_deref, (char*)pos_buf, sizeof(int32_t))) {
    goto merge_datasets_ret_NOMEM;
  }
  // pos_buf[n] contains the position of lexicographic marker #n in the hash
  // table.  invert this map, then traverse the hash table.
  for (uii = 0; uii < tot_marker_ct; uii++) {
    marker_map[pos_buf[uii]] = uii;
  }
  bigstack_end_reset(bigstack_end_mark); // deallocate second hash table
  ulii = 0;
  for (uii = 0; uii < HASHSIZE; uii++) {
    if (htable_bim[uii]) {
      llbim_ptr = htable_bim[uii];
      do {
	ujj = marker_map[ulii++];
	llxx = llbim_ptr->pos;
	pos_buf[ujj] = (uint32_t)llxx;
	bufptr = llbim_ptr->allele[0];
	if (bufptr) {
          marker_allele_ptrs[ujj * 2] = bufptr;
	}
	// already initialized to missing_geno_ptr otherwise

	bufptr = llbim_ptr->allele[1];
	if (bufptr) {
	  marker_allele_ptrs[ujj * 2 + 1] = bufptr;
	}
	marker_cms_tmp[ujj] = llbim_ptr->cm;
	ll_buf[ujj] = (((uint64_t)llxx) & 0xffffffff00000000LL) | ujj;
	llbim_ptr = llbim_ptr->next;
      } while (llbim_ptr);
    }
  }
  sort_marker_chrom_pos(ll_buf, tot_marker_ct, pos_buf, chrom_start, chrom_id, nullptr, &chrom_ct);
  // bugfix: when chromosomes are filtered out, flag the corresponding markers
  // in marker_map[]
  fill_uint_one(tot_marker_ct, marker_map);
  if (merge_post_msort_update_maps(marker_ids, max_marker_id_len, marker_map, marker_cms, marker_cms_tmp, pos_buf, ll_buf, chrom_start, chrom_id, chrom_ct, &dedup_marker_ct, merge_equal_pos, marker_allele_ptrs, chrom_info_ptr)) {
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  if ((!dedup_marker_ct) && (!allow_no_variants)) {
    logerrprint("Error: No variants in merged file.\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  bigstack_reset((char*)marker_cms_tmp);

  tot_sample_ct4 = (tot_sample_ct + 3) / 4;

  if (!keep_allele_order) {
    if (bigstack_calloc_ul(BITCT_TO_WORDCT(tot_marker_ct), &reversed)) {
      goto merge_datasets_ret_NOMEM;
    }
  }
  if (bigstack_alloc_ui(MAXV(max_cur_sample_ct, max_cur_marker_text_ct), &flex_map) ||
      bigstack_alloc_c(MAXV(max_marker_id_len, max_sample_id_len), &idbuf)) {
    goto merge_datasets_ret_NOMEM;
  }

  ulii = MAXV(tot_sample_ct4, ped_buflen);
  if (bigstack_alloc_uc(MAXV(ulii, 3), &readbuf)) {
    goto merge_datasets_ret_NOMEM;
  }
  if (merge_must_track_write(merge_mode)) {
    ulii = BITCT_TO_WORDCT(tot_sample_ct);
    if (ulii) {
      markers_per_pass = bigstack_left() / (3 * sizeof(intptr_t) * ulii);
      if (markers_per_pass > dedup_marker_ct) {
	markers_per_pass = dedup_marker_ct;
      }
    } else {
      markers_per_pass = dedup_marker_ct;
    }
    bigstack_alloc_ul(markers_per_pass * ulii, &markbuf);
  } else if (tot_sample_ct4) {
    markers_per_pass = bigstack_left() / tot_sample_ct4;
    if (markers_per_pass > dedup_marker_ct) {
      markers_per_pass = dedup_marker_ct;
    }
  } else {
    markers_per_pass = dedup_marker_ct;
  }
  if (dedup_marker_ct) {
    if (!markers_per_pass) {
      goto merge_datasets_ret_NOMEM;
    }
    pass_ct = 1 + ((dedup_marker_ct - 1) / markers_per_pass);
  } else {
    pass_ct = 0;
  }

  writebuf = g_bigstack_base;
  pcptr = (uintptr_t*)g_bigstack_base;
  if (merge_mode < 6) {
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    if (pass_ct == 1) {
      sprintf(g_logbuf, "Performing single-pass merge (%u %s, %u variant%s).\n", tot_sample_ct, species_str(tot_sample_ct), dedup_marker_ct, (dedup_marker_ct == 1)? "" : "s");
    } else {
      sprintf(g_logbuf, "Performing %u-pass merge (%u %s, %" PRIuPTR "/%u variant%s per pass).\n", pass_ct, tot_sample_ct, species_str(tot_sample_ct), markers_per_pass, dedup_marker_ct, (dedup_marker_ct == 1)? "" : "s");
    }
  } else {
    memcpy(outname_end, ".diff", 6);
    if (fopen_checked(outname, "w", &outfile)) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    if (fputs_checked("                 SNP                  FID                  IID      NEW      OLD \n", outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    LOGPREPRINTFWW("Performing %u-pass diff (mode %u), writing results to %s .\n", pass_ct, merge_mode, outname);
  }
  logprintb();
  for (uii = 0; uii < pass_ct; uii++) {
    if (uii + 1 == pass_ct) {
      ujj = dedup_marker_ct - markers_per_pass * uii;
    } else {
      ujj = markers_per_pass;
    }
    if (tot_sample_ct % 4) {
      umm = tot_sample_ct / 4;
      ubufptr = writebuf;
      ucc = 0x15 >> (6 - 2 * (tot_sample_ct % 4));
      for (ukk = 0; ukk < ujj; ukk++) {
        memset(ubufptr, 0x55, umm);
        ubufptr[umm] = ucc;
	ubufptr = &(ubufptr[tot_sample_ct4]);
      }
    } else {
      memset(writebuf, 0x55, ((uintptr_t)ujj) * tot_sample_ct4);
    }
    if (merge_must_track_write(merge_mode)) {
      fill_ulong_zero(ujj * ulii, markbuf);
    }
    for (mlpos = 0; mlpos < merge_ct; mlpos++) {
      retval = merge_main(mergelist_bed[mlpos], mergelist_bim[mlpos], mergelist_fam[mlpos], bim_loadbuf, max_bim_linelen, tot_sample_ct, tot_marker_ct, dedup_marker_ct, uii * markers_per_pass, ujj, marker_allele_ptrs, marker_ids, max_marker_id_len, sample_ids, max_sample_id_len, merge_nsort, sample_nsmap, flex_map, marker_map, idbuf, readbuf, writebuf, mlpos? merge_mode : merge_first_mode(merge_mode, merge_equal_pos), markbuf, outfile, &diff_total_overlap, &diff_not_both_genotyped, &diff_discordant, ped_buflen);
      if (retval) {
	goto merge_datasets_ret_1;
      }
      if (mlpos != merge_ct - 1) {
        printf("\rPass %u: fileset #%" PRIuPTR " complete.", uii + 1, mlpos + 1);
	fflush(stdout);
      }
    }
    if (merge_mode < 6) {
      if (!keep_allele_order) {
	for (ukk = 0; ukk < ujj; ukk++) {
	  uljj = ((uintptr_t)ukk) * tot_sample_ct4;
	  umm = popcount_chars(pcptr, uljj, uljj + tot_sample_ct4);
	  if (umm < tot_sample_ct) {
	    ulkk = (uii * markers_per_pass) + ukk;
	    reversed[ulkk / BITCT] |= (ONELU << (ulkk % BITCT));
	    reverse_loadbuf(tot_sample_ct, &(writebuf[uljj]));
	  }
	}
      }
      if (fwrite_checked(writebuf, ((uintptr_t)ujj) * tot_sample_ct4, outfile)) {
        goto merge_datasets_ret_WRITE_FAIL;
      }
    }
    fputs("\r                                              \r", stdout);
    if (uii + 1 != pass_ct) {
      LOGPRINTF("Pass %u complete.\n", uii + 1);
    }
  }
  if (fclose_null(&outfile)) {
    goto merge_datasets_ret_WRITE_FAIL;
  }
  bigstack_reset(flex_map);
  if (bigstack_alloc_ui(dedup_marker_ct, &map_reverse)) {
    goto merge_datasets_ret_NOMEM;
  }
  if (merge_mode < 6) {
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(outname, "w", &outfile)) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    uii = tot_marker_ct;
    while (uii--) {
      ujj = marker_map[uii];
      if (ujj != 0xffffffffU) {
        map_reverse[ujj] = uii;
      }
    }
    for (ulii = 0; ulii < chrom_ct; ulii++) {
      uii = chrom_start[ulii + 1];
      ujj = chrom_start[ulii];
      ukk = chrom_id[ulii];
      for (; ujj < uii; ujj++) {
	umm = map_reverse[ujj];
	bufptr = chrom_name_write(chrom_info_ptr, ukk, g_textbuf);
	fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile);
	if (keep_allele_order || (!IS_SET(reversed, ujj))) {
	  bufptr = marker_allele_ptrs[2 * umm];
	  bufptr2 = marker_allele_ptrs[2 * umm + 1];
	} else {
	  bufptr = marker_allele_ptrs[2 * umm + 1];
	  bufptr2 = marker_allele_ptrs[2 * umm];
	}
	fprintf(outfile, "\t%s\t%g\t%u\t%s\t%s\n", &(marker_ids[map_reverse[ujj] * max_marker_id_len]), marker_cms[ujj], pos_buf[ujj], bufptr, bufptr2);
      }
      if (ferror(outfile)) {
	goto merge_datasets_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    LOGPRINTFWW("Merged fileset written to %s.bed + %s.bim + %s.fam .\n", outname, outname, outname);
  } else {
    // undo the "not"
    diff_not_both_genotyped = diff_total_overlap - diff_not_both_genotyped;
    LOGPRINTF("%" PRIu64 " overlapping call%s, %" PRIu64 " nonmissing in both filesets.\n%" PRIu64 " concordant, for a concordance rate of %g.\n", diff_total_overlap, (diff_total_overlap == 1)? "" : "s", diff_not_both_genotyped, diff_not_both_genotyped - diff_discordant, 1.0 - (((double)diff_discordant) / ((double)diff_not_both_genotyped)));
  }

  forget_extra_chrom_names(1, chrom_info_ptr);
  while (0) {
  merge_datasets_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_datasets_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_datasets_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_datasets_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  merge_datasets_ret_INVALID_FORMAT_2:
    logerrprintb();
  merge_datasets_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 merge_datasets_ret_1:
  if (merge_equal_pos) {
    cleanup_allele_storage2(tot_marker_ct * 2, marker_allele_ptrs);
  } else {
    cleanup_allele_storage(2, tot_marker_ct * 2, marker_allele_ptrs);
  }
  fclose_cond(mergelistfile);
  fclose_cond(outfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return retval;
}
