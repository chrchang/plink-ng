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

#include "plink_cluster.h"
#include "plink_data.h"
#include "plink_dosage.h"
#include "plink_filter.h"
#include "plink_glm.h"
#include "plink_matrix.h"
#include "plink_misc.h"
#include "pigz.h"

void dosage_init(Dosage_info* doip) {
  doip->fname = nullptr;
  doip->modifier = 0;
  doip->skip0 = 0;
  doip->skip1 = 0;
  doip->skip2 = 0;
  doip->format = 2;
}

void dosage_cleanup(Dosage_info* doip) {
  free_cond(doip->fname);
}

#define DOSAGE_EPSILON 0.000244140625

int32_t dosage_load_score_files(Score_info* sc_ip, char* outname, char* outname_end, uint32_t double_dosage, uintptr_t* score_marker_ct_ptr, uintptr_t* max_score_marker_id_len_ptr, char** score_marker_ids_ptr, char*** score_allele_codes_ptr, double** score_effect_sizes_ptr, uintptr_t** score_qrange_key_exists_ptr, double** score_qrange_keys_ptr, uintptr_t* qrange_ct_ptr, uintptr_t* max_qrange_name_len_ptr, char** score_qrange_names_ptr, double** score_qrange_bounds_ptr) {
  // We don't necessarily have the whole variant ID list in advance, so it
  // makes sense to deviate a bit from score_report().
  //
  // 1. Process main --score file in three passes.  First pass, just determine
  //    variant count and maximum ID length; second pass, save and sort just
  //    the relevant variant IDs and check for duplicates; third pass, save
  //    allele codes and scores.
  // 2. If --q-score-range was specified, load those files.
  FILE* infile = nullptr;
  uintptr_t score_marker_ct = 0;
  uintptr_t max_score_marker_id_len = 0;
  uintptr_t miss_ct = 0;
  uintptr_t marker_idx = 0;
  uintptr_t qrange_ct = 0;
  uintptr_t max_qrange_name_len = 0;
  uint64_t allele_code_buf_len = 0;
  uintptr_t* score_qrange_key_exists = nullptr;
  double* score_qrange_keys = nullptr;
  uint32_t modifier = sc_ip->modifier;
  int32_t retval = 0;
  char* bufptr_arr[3];
  char* score_marker_ids;
  char* score_qrange_names;
  char** score_allele_codes;
  double* score_effect_sizes;
  double* score_qrange_bounds;
  char* allele_code_buf;
  char* loadbuf;
  char* bufptr;
  uintptr_t loadbuf_size;
  uintptr_t score_marker_ctl;
  uintptr_t line_idx;
  uintptr_t ulii;
  double lbound;
  double ubound;
  double dxx;
  uint32_t first_col_m1;
  uint32_t col_01_delta;
  uint32_t col_12_delta;
  uint32_t varid_idx;
  uint32_t allele_idx;
  uint32_t effect_idx;
  uint32_t rangename_len_limit;
  uint32_t slen;
  int32_t ii;
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto dosage_load_score_files_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  retval = open_and_load_to_first_token(&infile, sc_ip->fname, loadbuf_size, '\0', "--score file", loadbuf, &bufptr, &line_idx);
  if (retval) {
    goto dosage_load_score_files_ret_1;
  }

  if (sc_ip->varid_col < sc_ip->allele_col) {
    first_col_m1 = sc_ip->varid_col;
    varid_idx = 0;
    if (sc_ip->allele_col < sc_ip->effect_col) {
      col_01_delta = sc_ip->allele_col - first_col_m1;
      col_12_delta = sc_ip->effect_col - sc_ip->allele_col;
      allele_idx = 1;
      effect_idx = 2;
    } else {
      allele_idx = 2;
      if (sc_ip->varid_col < sc_ip->effect_col) {
        col_01_delta = sc_ip->effect_col - first_col_m1;
        col_12_delta = sc_ip->allele_col - sc_ip->effect_col;
	effect_idx = 1;
      } else {
        first_col_m1 = sc_ip->effect_col;
        col_01_delta = sc_ip->varid_col - first_col_m1;
        col_12_delta = sc_ip->allele_col - sc_ip->varid_col;
        varid_idx = 1;
        effect_idx = 0;
      }
    }
  } else {
    first_col_m1 = sc_ip->allele_col;
    allele_idx = 0;
    if (sc_ip->varid_col < sc_ip->effect_col) {
      col_01_delta = sc_ip->varid_col - first_col_m1;
      col_12_delta = sc_ip->effect_col - sc_ip->varid_col;
      varid_idx = 1;
      effect_idx = 2;
    } else {
      varid_idx = 2;
      if (sc_ip->allele_col < sc_ip->effect_col) {
        col_01_delta = sc_ip->effect_col - first_col_m1;
        col_12_delta = sc_ip->varid_col - sc_ip->effect_col;
	effect_idx = 1;
      } else {
        first_col_m1 = sc_ip->effect_col;
        col_01_delta = sc_ip->allele_col - first_col_m1;
        col_12_delta = sc_ip->varid_col - sc_ip->allele_col;
        allele_idx = 1;
        effect_idx = 0;
      }
    }
  }
  first_col_m1--;
  if (modifier & SCORE_HEADER) {
    goto dosage_load_score_files_next1;
  }
  // first pass: validate, count
  while (1) {
    bufptr_arr[0] = next_token_multz(bufptr, first_col_m1);
    bufptr_arr[1] = next_token_mult(bufptr_arr[0], col_01_delta);
    bufptr_arr[2] = next_token_mult(bufptr_arr[1], col_12_delta);
    if (!bufptr_arr[2]) {
      goto dosage_load_score_files_ret_MISSING_TOKENS;
    }
    if (scan_double(bufptr_arr[effect_idx], &dxx)) {
      miss_ct++;
    } else {
      ulii = strlen_se(bufptr_arr[varid_idx]) + 1;
      if (ulii > max_score_marker_id_len) {
	max_score_marker_id_len = ulii;
      }
      slen = strlen_se(bufptr_arr[allele_idx]);
      if (slen > 1) {
	allele_code_buf_len += slen + 1;
      }
      score_marker_ct++;
    }
  dosage_load_score_files_next1:
    line_idx++;
    if (!fgets(loadbuf, loadbuf_size, infile)) {
      break;
    }
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --score file is pathologically long.\n", line_idx);
        goto dosage_load_score_files_ret_INVALID_FORMAT_2;
      }
      goto dosage_load_score_files_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      goto dosage_load_score_files_next1;
    }
  }
  if (!score_marker_ct) {
    logerrprint("Error: No valid entries in --score file.\n");
    goto dosage_load_score_files_ret_INVALID_FORMAT;
  }
  if (score_marker_ct >= 0x40000000) {
    logerrprint("Error: --score does not support >= 2^30 variants.\n");
    goto dosage_load_score_files_ret_INVALID_FORMAT;
  }
#ifndef __LP64__
  if (allele_code_buf_len > 0x7fffffff) {
    goto dosage_load_score_files_ret_NOMEM;
  }
#endif
  if (bigstack_alloc_c(score_marker_ct * max_score_marker_id_len, score_marker_ids_ptr) ||
      bigstack_alloc_d(score_marker_ct, score_effect_sizes_ptr) ||
      bigstack_alloc_c((uintptr_t)allele_code_buf_len, &allele_code_buf)) {
    goto dosage_load_score_files_ret_NOMEM;
  }
  score_marker_ids = *score_marker_ids_ptr;
  score_effect_sizes = *score_effect_sizes_ptr;
  score_marker_ctl = BITCT_TO_WORDCT(score_marker_ct);
  if (sc_ip->data_fname) {
    if (bigstack_calloc_ul(score_marker_ctl, score_qrange_key_exists_ptr) ||
        bigstack_alloc_d(score_marker_ct, score_qrange_keys_ptr)) {
      goto dosage_load_score_files_ret_NOMEM;
    }
    score_qrange_key_exists = *score_qrange_key_exists_ptr;
    score_qrange_keys = *score_qrange_keys_ptr;
  }
  score_allele_codes = (char**)bigstack_alloc(score_marker_ct * sizeof(intptr_t));
  if (!score_allele_codes) {
    goto dosage_load_score_files_ret_NOMEM;
  }
  *score_allele_codes_ptr = score_allele_codes;
  rewind(infile);
  loadbuf_size = bigstack_left();
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto dosage_load_score_files_ret_NOMEM;
  }
  loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  // pass 2: load and sort variant IDs
  retval = load_to_first_token(infile, loadbuf_size, '\0', "--score file", loadbuf, &bufptr, &line_idx);
  if (modifier & SCORE_HEADER) {
    goto dosage_load_score_files_next2;
  }
  while (1) {
    bufptr_arr[0] = next_token_multz(bufptr, first_col_m1);
    bufptr_arr[1] = next_token_mult(bufptr_arr[0], col_01_delta);
    bufptr_arr[2] = next_token_mult(bufptr_arr[1], col_12_delta);
    if (!scan_double(bufptr_arr[effect_idx], &dxx)) {
      slen = strlen_se(bufptr_arr[varid_idx]);
      memcpyx(&(score_marker_ids[marker_idx * max_score_marker_id_len]), bufptr_arr[varid_idx], slen, '\0');
      if (++marker_idx == score_marker_ct) {
	break;
      }
    }
  dosage_load_score_files_next2:
    line_idx++;
    if (!fgets(loadbuf, loadbuf_size, infile)) {
      goto dosage_load_score_files_ret_READ_FAIL;
    }
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --score file is pathologically long.\n", line_idx);
        goto dosage_load_score_files_ret_INVALID_FORMAT_2;
      }
      goto dosage_load_score_files_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      goto dosage_load_score_files_next2;
    }
  }
  qsort(score_marker_ids, score_marker_ct, max_score_marker_id_len, strcmp_casted);
  bufptr = scan_for_duplicate_ids(score_marker_ids, score_marker_ct, max_score_marker_id_len);
  if (bufptr) {
    LOGPREPRINTFWW("Error: Duplicate variant '%s' in --score file.\n", bufptr);
    goto dosage_load_score_files_ret_INVALID_FORMAT_2;
  }
  rewind(infile);

  // pass 3: load other stuff
  retval = load_to_first_token(infile, loadbuf_size, '\0', "--score file", loadbuf, &bufptr, &line_idx);
  ulii = 0;
  if (modifier & SCORE_HEADER) {
    goto dosage_load_score_files_next3;
  }
  while (1) {
    bufptr_arr[0] = next_token_multz(bufptr, first_col_m1);
    bufptr_arr[1] = next_token_mult(bufptr_arr[0], col_01_delta);
    bufptr_arr[2] = next_token_mult(bufptr_arr[1], col_12_delta);
    if (!scan_double(bufptr_arr[effect_idx], &dxx)) {
      // guaranteed to succeed unless the user is overwriting the file between
      // load passes, which we won't bother defending against
      marker_idx = (uint32_t)bsearch_str(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), score_marker_ids, max_score_marker_id_len, score_marker_ct);
      if (double_dosage) {
	dxx *= 2;
      }
      score_effect_sizes[marker_idx] = dxx;
      slen = strlen_se(bufptr_arr[allele_idx]);
      if (slen == 1) {
	score_allele_codes[marker_idx] = (char*)(&(g_one_char_strs[2 * bufptr_arr[allele_idx][0]]));
      } else {
	score_allele_codes[marker_idx] = allele_code_buf;
        allele_code_buf = memcpyax(allele_code_buf, bufptr_arr[allele_idx], slen, '\0');
      }
      if (++ulii == score_marker_ct) {
	break;
      }
    }
  dosage_load_score_files_next3:
    if (!fgets(loadbuf, loadbuf_size, infile)) {
      goto dosage_load_score_files_ret_READ_FAIL;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      goto dosage_load_score_files_next3;
    }
  }
  if (fclose_null(&infile)) {
    goto dosage_load_score_files_ret_READ_FAIL;
  }
  LOGPRINTFWW("--score: %" PRIuPTR " entr%s loaded from %s.\n", score_marker_ct, (score_marker_ct == 1)? "y" : "ies", sc_ip->fname);
  if (miss_ct) {
    LOGERRPRINTF("Warning: %" PRIuPTR " line%s skipped.\n", miss_ct, (miss_ct == 1)? "" : "s");
  }
  *score_marker_ct_ptr = score_marker_ct;
  *max_score_marker_id_len_ptr = max_score_marker_id_len;
  if (sc_ip->data_fname) {
    retval = open_and_load_to_first_token(&infile, sc_ip->data_fname, loadbuf_size, '\0', "--q-score-range data file", loadbuf, &bufptr, &line_idx);
    if (retval) {
      goto dosage_load_score_files_ret_1;
    }
    if (sc_ip->data_varid_col < sc_ip->data_col) {
      first_col_m1 = sc_ip->data_varid_col;
      col_01_delta = sc_ip->data_col - first_col_m1;
      varid_idx = 0;
    } else {
      first_col_m1 = sc_ip->data_col;
      col_01_delta = sc_ip->data_varid_col - first_col_m1;
      varid_idx = 1;
    }
    first_col_m1--;
    miss_ct = 0;
    if (modifier & SCORE_DATA_HEADER) {
      goto dosage_load_score_files_next4;
    }
    while (1) {
      bufptr_arr[0] = next_token_multz(bufptr, first_col_m1);
      bufptr_arr[1] = next_token_mult(bufptr_arr[0], col_01_delta);
      if (!bufptr_arr[1]) {
        goto dosage_load_score_files_ret_MISSING_TOKENS_Q;
      }
      ii = bsearch_str(bufptr_arr[varid_idx], strlen_se(bufptr_arr[varid_idx]), score_marker_ids, max_score_marker_id_len, score_marker_ct);
      if (ii != -1) {
        if (scan_double(bufptr_arr[1 - varid_idx], &dxx) || (dxx != dxx)) {
	  miss_ct++;
	} else {
	  if (is_set(score_qrange_key_exists, ii)) {
	    LOGPREPRINTFWW("Error: Duplicate variant '%s' in --q-score-range data file.\n", &(score_marker_ids[((uint32_t)ii) * max_score_marker_id_len]));
	    goto dosage_load_score_files_ret_INVALID_FORMAT_2;
	  }
	  score_qrange_keys[(uint32_t)ii] = dxx;
	  set_bit(ii, score_qrange_key_exists);
	}
      } else {
	miss_ct++;
      }
    dosage_load_score_files_next4:
      line_idx++;
      if (!fgets(loadbuf, loadbuf_size, infile)) {
	break;
      }
      if (!(loadbuf[loadbuf_size - 1])) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --q-score-range data file is pathologically long.\n", line_idx);
	  goto dosage_load_score_files_ret_INVALID_FORMAT_2;
	}
	goto dosage_load_score_files_ret_NOMEM;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*bufptr)) {
	goto dosage_load_score_files_next4;
      }
    }
    if (fclose_null(&infile)) {
      goto dosage_load_score_files_ret_READ_FAIL;
    }
    if (miss_ct) {
      LOGERRPRINTF("Warning: %" PRIuPTR " line%s skipped in --q-score-range data file.\n", miss_ct, (miss_ct == 1)? "" : "s");
    }
    if (fopen_checked(sc_ip->range_fname, "r", &infile)) {
      goto dosage_load_score_files_ret_OPEN_FAIL;
    }
    rangename_len_limit = (FNAMESIZE - 10) - ((uintptr_t)(outname_end - outname));
    g_textbuf[MAXLINELEN - 1] = ' ';
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --q-score-range range file is pathologically long.\n", line_idx);
        goto dosage_load_score_files_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
        continue;
      }
      slen = strlen_se(bufptr);
      bufptr_arr[1] = skip_initial_spaces(&(bufptr[slen]));
      bufptr_arr[2] = next_token(bufptr_arr[1]);
      if ((!bufptr_arr[2]) || scan_double(bufptr_arr[1], &lbound) || scan_double(bufptr_arr[2], &ubound) || (lbound != lbound) || (ubound != ubound) || (lbound > ubound)) {
	// count these skips?
	continue;
      }
      if (slen > rangename_len_limit) {
        sprintf(g_logbuf, "Error: Excessively long range name on line %" PRIuPTR " of --q-score-range range\nfile.\n", line_idx);
        goto dosage_load_score_files_ret_INVALID_FORMAT_2;
      }
      qrange_ct++;
      if (slen >= max_qrange_name_len) {
	max_qrange_name_len = slen + 1;
      }
    }
    if (!feof(infile)) {
      goto dosage_load_score_files_ret_READ_FAIL;
    }
    if (!qrange_ct) {
      logerrprint("Error: No valid entries in --q-score-range range file.\n");
      goto dosage_load_score_files_ret_INVALID_FORMAT;
    }
    if (bigstack_alloc_c(qrange_ct * max_qrange_name_len, score_qrange_names_ptr) ||
        bigstack_alloc_d(qrange_ct * 2, score_qrange_bounds_ptr)) {
      goto dosage_load_score_files_ret_NOMEM;
    }
    score_qrange_names = *score_qrange_names_ptr;
    score_qrange_bounds = *score_qrange_bounds_ptr;
    rewind(infile);
    ulii = 0; // range index
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
        continue;
      }
      slen = strlen_se(bufptr);
      bufptr_arr[1] = skip_initial_spaces(&(bufptr[slen]));
      bufptr_arr[2] = next_token(bufptr_arr[1]);
      if ((!bufptr_arr[2]) || scan_double(bufptr_arr[1], &lbound) || scan_double(bufptr_arr[2], &ubound) || (lbound != lbound) || (ubound != ubound) || (lbound > ubound)) {
	continue;
      }
      memcpyx(&(score_qrange_names[ulii * max_qrange_name_len]), bufptr, slen, '\0');
      score_qrange_bounds[2 * ulii] = lbound;
      score_qrange_bounds[2 * ulii + 1] = ubound;
      ulii++;
    }
    if (ulii != qrange_ct) {
      // catches /dev/stdin redirection
      goto dosage_load_score_files_ret_READ_FAIL;
    }
    if (fclose_null(&infile)) {
      goto dosage_load_score_files_ret_READ_FAIL;
    }
    LOGPRINTF("--q-score-range: %" PRIuPTR " range%s loaded.\n", qrange_ct, (qrange_ct == 1)? "" : "s");
    *qrange_ct_ptr = qrange_ct;
    *max_qrange_name_len_ptr = max_qrange_name_len;
  }
  while (0) {
  dosage_load_score_files_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  dosage_load_score_files_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  dosage_load_score_files_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  dosage_load_score_files_ret_MISSING_TOKENS_Q:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --q-score-range data file has fewer tokens than\nexpected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  dosage_load_score_files_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --score file has fewer tokens than expected.\n", line_idx);
  dosage_load_score_files_ret_INVALID_FORMAT_2:
    logerrprintb();
  dosage_load_score_files_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 dosage_load_score_files_ret_1:
  fclose_cond(infile);
  return retval;
}

int32_t plink1_dosage(Dosage_info* doip, char* famname, char* mapname, char* outname, char* outname_end, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* keepfamname, char* removefamname, char* filtername, char* makepheno_str, char* phenoname_str, char* covar_fname, Two_col_params* qual_filter, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* filtervals_flattened, char* filter_attrib_fname, char* filter_attrib_liststr, char* filter_attrib_sample_fname, char* filter_attrib_sample_liststr, double qual_min_thresh, double qual_max_thresh, double thin_keep_prob, uint32_t thin_keep_ct, uint32_t min_bp_space, uint32_t mfilter_col, uint32_t fam_cols, int32_t missing_pheno, char* output_missing_pheno, uint32_t mpheno_col, uint32_t pheno_modifier, Chrom_info* chrom_info_ptr, double tail_bottom, double tail_top, uint64_t misc_flags, uint64_t filter_flags, uint32_t sex_missing_pheno, uint32_t update_sex_col, Cluster_info* cluster_ptr, int32_t marker_pos_start, int32_t marker_pos_end, int32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* snps_range_list_ptr, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t mwithin_col, uint32_t glm_modifier, double glm_vif_thresh, double output_min_p, int32_t known_procs, Score_info* sc_ip) {
  // sucks to duplicate so much, but this code will be thrown out later so
  // there's no long-term maintenance problem
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* phenofile = nullptr;
  FILE* infile = nullptr;
  FILE* profile_outfile = nullptr;
  gzFile* gz_infiles = nullptr;
  char* marker_ids = nullptr;
  char* sample_ids = nullptr;
  char* paternal_ids = nullptr;
  char* maternal_ids = nullptr;
  char* cluster_ids = nullptr;
  char* covar_names = nullptr;
  char* sorted_sample_ids = nullptr;
  char* sep_fnames = nullptr;
  char* cur_marker_id_buf = nullptr;
  char* score_marker_ids = nullptr;
  char* score_qrange_names = nullptr;
  char** score_allele_codes = nullptr;
  char* a1_ptr = nullptr;
  char* a2_ptr = nullptr;
  char* pzwritep = nullptr;
  uintptr_t* marker_exclude = nullptr;
  uintptr_t* sample_exclude = nullptr;
  uintptr_t* sex_nm = nullptr;
  uintptr_t* sex_male = nullptr;
  uintptr_t* pheno_nm = nullptr;
  uintptr_t* pheno_c = nullptr;
  uintptr_t* pheno_nm_collapsed = nullptr;
  uintptr_t* pheno_c_collapsed = nullptr;
  uintptr_t* founder_info = nullptr;
  uintptr_t* covar_nm = nullptr;
  uintptr_t* perm_vec = nullptr;
  uintptr_t* perm_fails = nullptr; // need to enforce alignment
  uintptr_t* score_qrange_key_exists = nullptr;
  double* pheno_d = nullptr;
  double* covar_d = nullptr;
  double* cur_dosages2 = nullptr;
  double* score_effect_sizes = nullptr;
  double* score_qrange_keys = nullptr;
  double* score_qrange_bounds = nullptr;
  double* cur_scores = nullptr;
  double* score_bases = nullptr;
#ifndef NOLAPACK
  double* pheno_d2 = nullptr;
  double* covars_cov_major_buf = nullptr;
  double* covars_sample_major_buf = nullptr;
  double* pheno_d_collapsed = nullptr;
  double* param_2d_buf = nullptr;
  double* param_2d_buf2 = nullptr;
  double* regression_results = nullptr;
  double* dgels_a = nullptr;
  double* dgels_b = nullptr;
  double* dgels_work = nullptr;
  MATRIX_INVERT_BUF1_TYPE* mi_buf = nullptr;
#endif
  float* covar_f = nullptr;
  float* covars_cov_major_f_buf = nullptr;
  float* coef_f = nullptr;
  float* pp_f = nullptr;
  float* sample_1d_buf_f = nullptr;
  float* pheno_buf_f = nullptr;
  float* param_1d_buf_f = nullptr;
  float* param_1d_buf2_f = nullptr;
  float* param_2d_buf_f = nullptr;
  float* param_2d_buf2_f = nullptr;
  float* regression_results_f = nullptr;
  Ll_ctstr_entry** htable = nullptr;
  uint32_t* marker_pos = nullptr;
  uint32_t* cluster_map = nullptr;
  uint32_t* cluster_starts = nullptr;
  uint32_t* marker_id_htable = nullptr;
  uint32_t* sample_id_map = nullptr;
  uint32_t* batch_sizes = nullptr;
  uint32_t* score_range_obs_cts = nullptr;
  uint32_t* score_miss_cts = nullptr;
  uint32_t* uiptr = nullptr;
  uint32_t* uiptr2 = nullptr;
  uint32_t* uiptr3 = nullptr;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t max_marker_id_blen = 0;
  uintptr_t unfiltered_sample_ct = 0;
  uintptr_t sample_exclude_ct = 0;
  uintptr_t max_sample_id_len = 4;
  uintptr_t max_paternal_id_len = 2;
  uintptr_t max_maternal_id_len = 2;
  uintptr_t cluster_ct = 0;
  uintptr_t max_cluster_id_len = 2;
  uintptr_t covar_ct = 0;
  uintptr_t max_covar_name_len = 0;
  uintptr_t max_fn_len = 0;
  uintptr_t max_sepheader_len = 0;
  uintptr_t file_idx_start = 0;
  uintptr_t distinct_id_ct = 0; // occur mode
  uintptr_t max_occur_id_len = 0;
  uintptr_t marker_idx = 0;
  uintptr_t score_marker_ct = 0;
  uintptr_t max_score_marker_id_len = 0;
  uintptr_t qrange_ct = 0;
  uintptr_t max_qrange_name_len = 0;
  uintptr_t ulii = 0;
  double missing_phenod = (double)missing_pheno;
  uint32_t load_map = (mapname[0] != '\0');
  uint32_t do_glm = (doip->modifier / DOSAGE_GLM) & 1;
  uint32_t do_score = (doip->modifier / DOSAGE_SCORE) & 1;
  uint32_t count_occur = doip->modifier & DOSAGE_OCCUR;
  uint32_t sepheader = (doip->modifier / DOSAGE_SEPHEADER) & 1;
  uint32_t noheader = doip->modifier & DOSAGE_NOHEADER;
  uint32_t output_gz = doip->modifier & DOSAGE_ZOUT;
  uint32_t dose1 = doip->modifier & DOSAGE_DOSE1;
  uint32_t score_report_average = doip->modifier & DOSAGE_SCORE_NOSUM;
  uint32_t dosage_score_cnt = doip->modifier & DOSAGE_SCORE_CNT;
  uint32_t sex_covar = doip->modifier & DOSAGE_SEX;
  uint32_t freq_cc = doip->modifier & DOSAGE_FREQ_CC;
  uint32_t skip0 = doip->skip0;
  uint32_t skip1p1 = doip->skip1 + 1;
  uint32_t skip2 = doip->skip2;
  uint32_t format_val = doip->format;
  uint32_t standard_beta = glm_modifier & GLM_STANDARD_BETA;
  uint32_t score_center = sc_ip->modifier & SCORE_CENTER;
  uint32_t score_mean_impute = !(sc_ip->modifier & SCORE_NO_MEAN_IMPUTATION);
  uint32_t map_is_unsorted = 0;
  uint32_t affection = 0;
  uint32_t infile_ct = 0;
  uint32_t pheno_ctrl_ct = 0;
  uint32_t batch_ct = 1;
  uint32_t max_batch_size = 1;
  uint32_t cur_marker_id_len = 0;
  uint32_t marker_id_htable_size = 0;
  uint32_t score_marker_idx = 0;
  uint32_t score_a2_effect = 0;
  uint32_t plink_maxfid = 0;
  uint32_t plink_maxiid = 0;
  uint32_t missing_pheno_len = 0;
  uint32_t a1_len = 0;
  uint32_t a2_len = 0;
  uint32_t ujj = 0;
  int32_t retval = 0;
#ifndef NOLAPACK
  char dgels_trans = 'N';
  __CLPK_integer dgels_m;
  __CLPK_integer dgels_n;
  __CLPK_integer dgels_nrhs;
  __CLPK_integer dgels_ldb;
  __CLPK_integer dgels_info;
  __CLPK_integer dgels_lwork;
#endif
  char missing_pheno_str[32];
  Pigz_state ps;
  unsigned char* bigstack_mark;
  unsigned char* overflow_buf;
  char* fnames;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  char* bufptr6;
  uintptr_t* line_idx_arr;
  uintptr_t* batch_samples;
  uintptr_t* cur_samples;
  double* cur_dosages;
  double* dptr;
  float* fptr;
  Ll_ctstr_entry** ll_pptr;
  Ll_ctstr_entry* ll_ptr;
  uint32_t* file_icts;
  uint32_t* read_idx_to_sample_idx;
  uint32_t* skip_vals;
  uintptr_t marker_ct;
  uintptr_t unfiltered_sample_ctl;
  uintptr_t sample_ct;
  uintptr_t sample_cta4;
  uintptr_t sample_ctl;
  uintptr_t line_idx;
  uintptr_t file_idx;
  uintptr_t sample_idx;
  uintptr_t cur_batch_size;
  uintptr_t loadbuf_size;
  uintptr_t sample_valid_ct;
  uintptr_t param_ct;
  uintptr_t param_cta4;
  uintptr_t qrange_idx;
  uintptr_t uljj;
  double sample_valid_ct_recip;
  double rsq;
  double pval;
  double beta;
  double se;
  double score_cur_effect_size;
  double score_missing_effect;
  double dxx;
  double dyy;
  double dzz;
  uint32_t gender_unk_ct;
  uint32_t pheno_nm_ct;
  uint32_t batch_idx;
  uint32_t read_idx_start;
  uint32_t read_idx;
  uint32_t slen;
  uint32_t sample_uidx;
  uint32_t prev_sample_uidx;
  uint32_t is_valid;
  uint32_t uii;
  uint32_t ukk;
  int32_t ii;
  pzwrite_init_null(&ps);
  if (load_map) {
    retval = load_bim(mapname, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_blen, &marker_exclude, nullptr, nullptr, nullptr, &ulii, &marker_ids, nullptr, 0, nullptr, chrom_info_ptr, nullptr, &marker_pos, misc_flags, filter_flags, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, &map_is_unsorted, do_glm || min_bp_space || (misc_flags & (MISC_EXTRACT_RANGE | MISC_EXCLUDE_RANGE)), 0, 0, nullptr, ".map file", nullptr);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    if (map_is_unsorted & UNSORTED_SPLIT_CHROM) {
      logerrprint("Error: .map file has a split chromosome.\n");
      goto plink1_dosage_ret_INVALID_FORMAT;
    }
  }
  uii = fam_cols & FAM_COL_6;
  if (uii && phenoname) {
    uii = (pheno_modifier & PHENO_MERGE) && (!makepheno_str);
  }
  if (!uii) {
    pheno_modifier &= ~PHENO_MERGE;
  }
  if (update_ids_fname) {
    ulii = 0;
    retval = scan_max_fam_indiv_strlen(update_ids_fname, 3, &max_sample_id_len);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  } else if (update_parents_fname) {
    retval = scan_max_strlen(update_parents_fname, 3, 4, 0, '\0', &max_paternal_id_len, &max_maternal_id_len);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  }
  retval = load_fam(famname, fam_cols, uii, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, &unfiltered_sample_ct, &sample_ids, &max_sample_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &sample_exclude, 0);
  if (retval) {
    goto plink1_dosage_ret_1;
  }
  unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  if (misc_flags & MISC_MAKE_FOUNDERS_FIRST) {
    if (make_founders(unfiltered_sample_ct, unfiltered_sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, sample_exclude, founder_info)) {
      goto plink1_dosage_ret_NOMEM;
    }
  }
  count_genders(sex_nm, sex_male, sample_exclude, unfiltered_sample_ct, &uii, &ujj, &gender_unk_ct);
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (gender_unk_ct) {
    LOGPRINTF("%" PRIuPTR " %s (%u male%s, %u female%s, %u ambiguous) loaded from .fam.\n", unfiltered_sample_ct, species_str(unfiltered_sample_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s", gender_unk_ct);
    retval = write_nosex(outname, outname_end, unfiltered_sample_ct, sample_exclude, sex_nm, gender_unk_ct, sample_ids, max_sample_id_len);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  } else {
    LOGPRINTF("%" PRIuPTR " %s (%d male%s, %d female%s) loaded from .fam.\n", unfiltered_sample_ct, species_str(unfiltered_sample_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s");
  }
  uii = popcount_longs(pheno_nm, unfiltered_sample_ctl);
  if (uii) {
    LOGPRINTF("%u phenotype value%s loaded from .fam.\n", uii, (uii == 1)? "" : "s");
  }
  if (phenoname && fopen_checked(phenoname, "r", &phenofile)) {
    goto plink1_dosage_ret_OPEN_FAIL;
  }
  if (phenofile || update_ids_fname || update_parents_fname || update_sex_fname || (filter_flags & FILTER_TAIL_PHENO)) {
    bigstack_mark = g_bigstack_base;
    retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, 0, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &sorted_sample_ids, &sample_id_map);
    if (retval) {
      goto plink1_dosage_ret_1;
    }

    if (makepheno_str) {
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_sample_ct, sorted_sample_ids, max_sample_id_len, sample_id_map, pheno_nm, &pheno_c);
      if (retval) {
        goto plink1_dosage_ret_1;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_sample_ct, 0, sorted_sample_ids, max_sample_id_len, sample_id_map, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d, nullptr, 0);
      if (retval) {
	if (retval == LOAD_PHENO_LAST_COL) {
	  logprintb();
	  retval = RET_INVALID_FORMAT;
	  bigstack_reset(bigstack_mark);
	}
        goto plink1_dosage_ret_1;
      }
    }
    if (filter_flags & FILTER_TAIL_PHENO) {
      retval = convert_tail_pheno(unfiltered_sample_ct, pheno_nm, &pheno_c, &pheno_d, tail_bottom, tail_top, missing_phenod);
      if (retval) {
        goto plink1_dosage_ret_1;
      }
    }
    bigstack_reset(bigstack_mark);
  }
  if (load_map) {
    uii = update_map || update_name || filter_attrib_fname || qual_filter || snps_range_list_ptr->names;
    if (uii || extractname || excludename) {
      bigstack_mark = g_bigstack_base;
      retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_blen, !uii, &marker_id_htable_size, &marker_id_htable);
      if (retval) {
	goto plink1_dosage_ret_1;
      }

      if (snps_range_list_ptr->names) {
        retval = snps_flag(marker_ids, marker_id_htable, snps_range_list_ptr, unfiltered_marker_ct, max_marker_id_blen, marker_id_htable_size, (filter_flags / FILTER_EXCLUDE_MARKERNAME_SNP) & 1, marker_exclude, &marker_exclude_ct);
        if (retval) {
          goto plink1_dosage_ret_1;
        }
      }
      if (update_map) {
	retval = update_marker_pos(update_map, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, marker_pos, &map_is_unsorted, chrom_info_ptr);
      } else if (update_name) {
	retval = update_marker_names(update_name, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen, unfiltered_marker_ct);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
	if (extractname || excludename) {
	  bigstack_reset(bigstack_mark);
	  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_blen, 0, &marker_id_htable_size, &marker_id_htable);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	}
      }
      if (extractname) {
	if (!(misc_flags & MISC_EXTRACT_RANGE)) {
	  retval = extract_exclude_flag_norange(extractname, marker_id_htable, marker_id_htable_size, 0, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	} else {
	  if (map_is_unsorted & UNSORTED_BP) {
	    logerrprint("Error: '--extract range' requires a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	    goto plink1_dosage_ret_INVALID_CMDLINE;
	  }
	  retval = extract_exclude_range(extractname, marker_pos, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0, 0, chrom_info_ptr);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	  uljj = unfiltered_marker_ct - marker_exclude_ct;
	  LOGPRINTF("--extract range: %" PRIuPTR " variant%s remaining.\n", uljj, (uljj == 1)? "" : "s");
	}
      }
      if (excludename) {
	if (!(misc_flags & MISC_EXCLUDE_RANGE)) {
	  retval = extract_exclude_flag_norange(excludename, marker_id_htable, marker_id_htable_size, 1, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	} else {
	  if (map_is_unsorted & UNSORTED_BP) {
	    logerrprint("Error: '--exclude range' requires a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	    goto plink1_dosage_ret_INVALID_CMDLINE;
	  }
	  retval = extract_exclude_range(excludename, marker_pos, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1, 0, chrom_info_ptr);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	  uljj = unfiltered_marker_ct - marker_exclude_ct;
	  LOGPRINTF("--exclude range: %" PRIuPTR " variant%s remaining.\n", uljj, (uljj == 1)? "" : "s");
	}
      }
      if (filter_attrib_fname) {
	retval = filter_attrib(filter_attrib_fname, filter_attrib_liststr, marker_id_htable, marker_id_htable_size, 0, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      if (qual_filter) {
	retval = filter_qual_scores(qual_filter, qual_min_thresh, qual_max_thresh, marker_id_htable, marker_id_htable_size, 0, marker_ids, max_marker_id_blen, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      bigstack_reset(bigstack_mark);
    }
    if (thin_keep_prob != 1.0) {
      if (random_thin_markers(thin_keep_prob, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0)) {
	goto plink1_dosage_ret_ALL_MARKERS_EXCLUDED;
      }
    } else if (thin_keep_ct != 0xffffffffU) {
      retval = random_thin_markers_ct(thin_keep_ct, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
  }
  if (update_ids_fname || update_parents_fname || update_sex_fname || keepname || keepfamname || removename || removefamname || filter_attrib_sample_fname || filtername) {
    bigstack_mark = g_bigstack_base;
    retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, sample_exclude_ct, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, &sorted_sample_ids, &sample_id_map);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    ulii = unfiltered_sample_ct - sample_exclude_ct;
    if (update_ids_fname) {
      retval = update_sample_ids(update_ids_fname, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, sample_ids);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    } else {
      if (update_parents_fname) {
	retval = update_sample_parents(update_parents_fname, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      if (update_sex_fname) {
        retval = update_sample_sexes(update_sex_fname, update_sex_col, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, sex_nm, sex_male);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
    }
    if (keepfamname) {
      retval = keep_or_remove(keepfamname, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 2, 0);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (keepname) {
      retval = keep_or_remove(keepname, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 0, 0);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (removefamname) {
      retval = keep_or_remove(removefamname, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 3, 0);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (removename) {
      retval = keep_or_remove(removename, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, 1, 0);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (filter_attrib_sample_fname) {
      retval = filter_attrib_sample(filter_attrib_sample_fname, filter_attrib_sample_liststr, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, unfiltered_sample_ct, 0, sample_exclude, &sample_exclude_ct);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_samples_file(filtername, sorted_sample_ids, ulii, max_sample_id_len, sample_id_map, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, filtervals_flattened, mfilter_col, 0);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    bigstack_reset(bigstack_mark);
  }
  if (gender_unk_ct && (!(sex_missing_pheno & ALLOW_NO_SEX))) {
    uii = popcount_longs_exclude(pheno_nm, sex_nm, unfiltered_sample_ctl);
    if (uii) {
      bitvec_and(sex_nm, unfiltered_sample_ctl, pheno_nm);
      logerrprint("Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those\nphenotypes to be ignored, use the --allow-no-sex flag.\n");
    }
  }
  if (do_glm || (filter_flags & FILTER_PRUNE)) {
    ulii = sample_exclude_ct;
    bitvec_ornot(pheno_nm, unfiltered_sample_ctl, sample_exclude);
    zero_trailing_bits(unfiltered_sample_ct, sample_exclude);
    sample_exclude_ct = popcount_longs(sample_exclude, unfiltered_sample_ctl);
    uii = do_glm && (!(filter_flags & FILTER_PRUNE));
    if (sample_exclude_ct == unfiltered_sample_ct) {
      LOGERRPRINTF("Error: All %s removed by %s--prune.\n", g_species_plural, uii? "automatic " : "");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    if ((filter_flags & FILTER_PRUNE) || (sample_exclude_ct != ulii)) {
      LOGPRINTF("%s--prune: %" PRIuPTR " %s remaining.\n", uii? "Automatic " : "", unfiltered_sample_ct - sample_exclude_ct, species_str(unfiltered_sample_ct == sample_exclude_ct + 1));
    }
  }

  if (filter_flags & (FILTER_BINARY_CASES | FILTER_BINARY_CONTROLS)) {
    if (!pheno_c) {
      logerrprint("Error: --filter-cases/--filter-controls requires a case/control phenotype.\n");
      goto plink1_dosage_ret_INVALID_CMDLINE;
    }
    ii = sample_exclude_ct;
    filter_samples_bitfields(unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, pheno_c, (filter_flags / FILTER_BINARY_CASES) & 1, pheno_nm);
    if (sample_exclude_ct == unfiltered_sample_ct) {
      LOGERRPRINTF("Error: All %s removed due to case/control status (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_CASES)? "cases" : "controls");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = sample_exclude_ct - ii;
    LOGPRINTF("%d %s removed due to case/control status (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_CASES)? "cases" : "controls");
  }
  if (filter_flags & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
    ii = sample_exclude_ct;
    filter_samples_bitfields(unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, sex_male, (filter_flags / FILTER_BINARY_MALES) & 1, sex_nm);
    if (sample_exclude_ct == unfiltered_sample_ct) {
      LOGERRPRINTF("Error: All %s removed due to gender filter (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_MALES)? "males" : "females");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = sample_exclude_ct - ii;
    LOGPRINTF("%d %s removed due to gender filter (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_MALES)? "males" : "females");
  }
  if (filter_flags & (FILTER_BINARY_FOUNDERS | FILTER_BINARY_NONFOUNDERS)) {
    ii = sample_exclude_ct;
    filter_samples_bitfields(unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, founder_info, (filter_flags / FILTER_BINARY_FOUNDERS) & 1, nullptr);
    if (sample_exclude_ct == unfiltered_sample_ct) {
      LOGERRPRINTF("Error: All %s removed due to founder status (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = sample_exclude_ct - ii;
    LOGPRINTF("%d %s removed due to founder status (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
  }
  if (cluster_ptr->fname || (misc_flags & MISC_FAMILY_CLUSTERS)) {
    if (cluster_ptr->keep_fname || cluster_ptr->keep_flattened || cluster_ptr->remove_fname || cluster_ptr->remove_flattened) {
      retval = load_clusters(cluster_ptr->fname, unfiltered_sample_ct, sample_exclude, &sample_exclude_ct, sample_ids, max_sample_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, cluster_ptr->keep_fname, cluster_ptr->keep_flattened, cluster_ptr->remove_fname, cluster_ptr->remove_flattened, 0);
      if (retval) {
        goto plink1_dosage_ret_1;
      }
    } else {
      logerrprint("Warning: Ignoring --within/--family since it has no effect.  (PLINK 1.07's\nundocumented Huber-White standard error computation is currently disabled.)\n");
    }
  }
  sample_ct = unfiltered_sample_ct - sample_exclude_ct;
  if (!sample_ct) {
    LOGERRPRINTF("Error: No %s pass QC.\n", g_species_plural);
    goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
  }
  sample_cta4 = round_up_pow2(sample_ct, 4);
  sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uii = do_glm && pheno_d;
  if (g_thread_ct > 1) {
    if (output_gz) {
      LOGPRINTF("Using up to %u threads (change this with --threads).\n", g_thread_ct);
    } else {
      if (uii) {
	logprint("Using 1 thread.\n");
      } else {
        logprint("Using 1 thread (no multithreaded calculations invoked).\n");
      }
    }
  } else {
    logprint("Using 1 thread.\n");
  }
#ifndef NOLAPACK
  if (uii && ((!known_procs) || (known_procs >= 2 * ((int32_t)g_thread_ct)))) {
    logerrprint("Warning: This run includes BLAS/LAPACK linear algebra operations which\ncurrently disregard the --threads limit.  If this is problematic, you may want\nto recompile against single-threaded BLAS/LAPACK.\n");
  }
#endif
  if ((filter_flags & FILTER_MAKE_FOUNDERS) && (!(misc_flags & MISC_MAKE_FOUNDERS_FIRST))) {
    if (make_founders(unfiltered_sample_ct, sample_ct, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, sample_exclude, founder_info)) {
      goto plink1_dosage_ret_NOMEM;
    }
  }
  if (covar_fname) {
    // update this as more covariate-referencing commands are added
    if (!do_glm) {
      logerrprint("Warning: Ignoring --covar since no commands reference the covariates.\n");
    } else {
      retval = load_covars(covar_fname, unfiltered_sample_ct, sample_exclude, sample_ct, sex_covar? sex_nm : nullptr, sex_covar? sex_male : nullptr, sample_ids, max_sample_id_len, missing_phenod, covar_modifier, covar_range_list_ptr, 0, &covar_ct, &covar_names, &max_covar_name_len, pheno_nm, &covar_nm, &covar_d, nullptr, nullptr);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
  } else if (sex_covar) {
    if (bigstack_alloc_c(4, &covar_names) ||
        bigstack_alloc_ul(sample_ctl, &covar_nm) ||
        bigstack_alloc_d(sample_ct, &covar_d)) {
      goto plink1_dosage_ret_NOMEM;
    }
    covar_ct = 1;
    max_covar_name_len = 4;
    memcpy(covar_names, "SEX", 4);
    fill_all_bits(sample_ct, covar_nm);
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(sample_exclude, &sample_uidx);
      if (is_set(sex_nm, sample_uidx)) {
        covar_d[sample_idx] = (double)((int32_t)is_set(sex_male, sample_uidx));
      } else {
	CLEAR_BIT(sample_idx, covar_nm);
        covar_d[sample_idx] = missing_phenod;
      }
    }
  }
  param_ct = covar_ct + 2;
  param_cta4 = round_up_pow2(param_ct, 4);
  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, pheno_nm);
  if (pheno_c) {
    bitvec_and(pheno_nm, unfiltered_sample_ctl, pheno_c);
  }
  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, founder_info);
  bitvec_andnot(sample_exclude, unfiltered_sample_ctl, sex_nm);
  if (gender_unk_ct) {
    gender_unk_ct = sample_ct - popcount_longs(sex_nm, unfiltered_sample_ctl);
  }
  bitvec_and(sex_nm, unfiltered_sample_ctl, sex_male);

  pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);

  if (load_map) {
    if (unfiltered_marker_ct == marker_exclude_ct) {
      logerrprint("Error: No variants remaining.\n");
      goto plink1_dosage_ret_ALL_MARKERS_EXCLUDED;
    }
    if (min_bp_space) {
      if (map_is_unsorted & UNSORTED_BP) {
	logerrprint("Error: --bp-space requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
	goto plink1_dosage_ret_INVALID_FORMAT;
      }
      enforce_min_bp_space(min_bp_space, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, chrom_info_ptr);
    }
    marker_ct = unfiltered_marker_ct - marker_exclude_ct;
    retval = alloc_and_populate_id_htable(unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_blen, 0, &marker_id_htable_size, &marker_id_htable);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    LOGPRINTF("%" PRIuPTR " variant%s and %" PRIuPTR " %s pass filters and QC.\n", marker_ct, (marker_ct == 1)? "" : "s", sample_ct, species_str(sample_ct));
  } else {
    LOGPRINTF("%" PRIuPTR " %s pass%s filters and QC.\n", sample_ct, species_str(sample_ct), (sample_ct == 1)? "es" : "");
  }
  if (!pheno_nm_ct) {
    if (do_glm) {
      logerrprint("Error: No phenotypes present.\n");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    logprint("Note: No phenotypes present.\n");
  } else if (pheno_c) {
    pheno_ctrl_ct = pheno_nm_ct - popcount_longs(pheno_c, unfiltered_sample_ctl);
    if (pheno_nm_ct != sample_ct) {
      sprintf(g_logbuf, "Among remaining phenotypes, %u %s and %u %s.  (%" PRIuPTR " phenotype%s missing.)\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "is a case" : "are cases", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "is a control" : "are controls", sample_ct - pheno_nm_ct, (sample_ct - pheno_nm_ct == 1)? " is" : "s are");
    } else {
      sprintf(g_logbuf, "Among remaining phenotypes, %u %s and %u %s.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "is a case" : "are cases", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "is a control" : "are controls");
    }
    wordwrapb(0);
    logprintb();
    if (standard_beta) {
      logerrprint("Error: --dosage 'standard-beta' modifier cannot be used with a case/control\nphenotype.\n");
      goto plink1_dosage_ret_INVALID_CMDLINE;
    }
  } else {
    logprint("Phenotype data is quantitative.\n");
#ifdef NOLAPACK
    if (do_glm) {
      logerrprint("Error: --dosage linear regression requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
      goto plink1_dosage_ret_INVALID_CMDLINE;
    }
#endif
  }
  if (do_score) {
    retval = dosage_load_score_files(sc_ip, outname, outname_end, (doip->modifier & DOSAGE_SCORE_DOUBLE), &score_marker_ct, &max_score_marker_id_len, &score_marker_ids, &score_allele_codes, &score_effect_sizes, &score_qrange_key_exists, &score_qrange_keys, &qrange_ct, &max_qrange_name_len, &score_qrange_names, &score_qrange_bounds);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    if (qrange_ct) {
      // bugfix (24 Oct 2019): score_range_obs_cts and score_miss_cts need to
      // be zero-initialized
      if (bigstack_calloc_d(sample_ct * qrange_ct, &cur_scores) ||
          bigstack_calloc_d(qrange_ct, &score_bases) ||
          bigstack_calloc_ui(qrange_ct, &score_range_obs_cts) ||
          bigstack_calloc_ui(sample_ct * qrange_ct, &score_miss_cts)) {
	goto plink1_dosage_ret_NOMEM;
      }
      *outname_end = '.';
    } else {
      if (bigstack_calloc_d(sample_ct, &cur_scores) ||
          bigstack_calloc_d(1, &score_bases) ||
          bigstack_calloc_ui(1, &score_range_obs_cts) ||
          bigstack_calloc_ui(sample_ct, &score_miss_cts)) {
	goto plink1_dosage_ret_NOMEM;
      }
    }
    calc_plink_maxfid(unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, &plink_maxfid, &plink_maxiid);
    missing_pheno_len = strlen(output_missing_pheno);
    if (missing_pheno_len < 6) {
      memset(missing_pheno_str, 32, 6 - missing_pheno_len);
      memcpy(&(missing_pheno_str[6 - missing_pheno_len]), output_missing_pheno, missing_pheno_len);
      missing_pheno_len = 6;
    } else {
      memcpy(missing_pheno_str, output_missing_pheno, missing_pheno_len);
    }
  }

  // now the actual --dosage logic begins
  // 1. either load single file, or
  //    a. determine whether batch numbers are present; if yes, determine set
  //       of batch numbers and frequencies, sort
  //    b. determine max file path length(s) on first pass as well
  //    c. actually load on second pass
  // 2. initialize output writer/compressor, header line if necessary
  // 3. loop through batches
  //    a. partial memory reset, loop through file headers (determine sample
  //       counts in each file, skip counts, and ID mappings, error out on FID
  //       duplication)
  //    b. if association analysis, collapse phenotype/covariate data to
  //       account for missing samples in current batch
  //    c. process one variant at a time, all files in batch must have them in
  //       same order
  //    d. to enable parallel gzip, we can append to a write buffer and flush
  //       it when it's within range of overflowing on the next line.  this
  //       requires stable write buffer sizing, though, so we'll just use
  //       serial compression until the rest of the code is finished and we
  //       know how large all our other memory allocations are.
  // 4. final write loop if necessary
  if (doip->modifier & DOSAGE_LIST) {
    retval = open_and_load_to_first_token(&infile, doip->fname, MAXLINELEN, '\0', "--dosage list file", g_textbuf, &bufptr, &line_idx);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    batch_ct = count_tokens(bufptr) - sepheader - 1; // underflow ok
    if (batch_ct) {
      if (batch_ct != 1) {
        logerrprint("Error: Unexpected number of columns in --dosage list file.\n");
	goto plink1_dosage_ret_INVALID_FORMAT;
      }
      batch_sizes = (uint32_t*)g_bigstack_base;
      uiptr = batch_sizes;
      uiptr2 = (uint32_t*)(&(g_bigstack_base[bigstack_left() / 2]));
    } else {
      if (bigstack_alloc_ui(1, &batch_sizes)) {
	goto plink1_dosage_ret_NOMEM;
      }
    }
    // now batch_ct = 0 indicates no batch column (but one actual batch)
    while (1) {
      if (batch_ct) {
	if (uiptr == uiptr2) {
	  goto plink1_dosage_ret_NOMEM;
	}
	bufptr2 = token_endnn(bufptr);
        if (scan_int32(bufptr, (int32_t*)uiptr)) {
	  sprintf(g_logbuf, "Error: Invalid batch number on line %" PRIuPTR " of --dosage list file.\n", line_idx);
          goto plink1_dosage_ret_INVALID_FORMAT_2;
	}
	uiptr++;
	bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
          sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --dosage list file has fewer tokens than expected.\n", line_idx);
	  goto plink1_dosage_ret_INVALID_FORMAT_2;
	}
      }
      bufptr2 = token_endnn(bufptr);
      ulii = (uintptr_t)(bufptr2 - bufptr);
      if (ulii >= max_fn_len) {
	max_fn_len = ulii + 1;
      }
      if (sepheader) {
	bufptr = skip_initial_spaces(bufptr2);
	if (!bufptr) {
          sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --dosage list file has fewer tokens than expected.\n", line_idx);
	  goto plink1_dosage_ret_INVALID_FORMAT_2;
	}
	bufptr2 = token_endnn(bufptr);
	ulii = (uintptr_t)(bufptr2 - bufptr);
	if (ulii >= max_sepheader_len) {
          max_sepheader_len = ulii + 1;
	}
      }
      bufptr = skip_initial_spaces(bufptr2);
      if (!is_eoln_kns(*bufptr)) {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --dosage list file has more tokens than expected.\n", line_idx);
	goto plink1_dosage_ret_INVALID_FORMAT_2;
      }
      infile_ct++;
    plink1_dosage_next_list_line:
      if (!fgets(g_textbuf, MAXLINELEN, infile)) {
	if (ferror(infile)) {
	  goto plink1_dosage_ret_READ_FAIL;
	}
	break;
      }
      line_idx++;
      if (!g_textbuf[MAXLINELEN - 1]) {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --dosage list file is pathologically long.\n", line_idx);
        goto plink1_dosage_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
	goto plink1_dosage_next_list_line;
      }
    }
    // infile_ct must be at least 1 due to use of open_and_load_to_first_token
    if (batch_ct) {
#ifdef __cplusplus
      std::sort((int32_t*)batch_sizes, (int32_t*)uiptr);
#else
      qsort(batch_sizes, infile_ct, sizeof(int32_t), intcmp2);
#endif
      bigstack_alloc_ui(infile_ct, &batch_sizes);

      // temporary batch size buffer
      bigstack_end_alloc_ui(infile_ct, &uiptr3);

      uii = batch_sizes[0];
      uiptr2 = &(batch_sizes[1]);
      uiptr3[0] = 1;
      batch_ct = 1;
      while (uiptr2 < uiptr) {
	ujj = *uiptr2++;
	if (ujj != uii) {
          uiptr3[batch_ct++] = 1;
	  uii = ujj;
	} else {
	  uiptr3[batch_ct - 1] += 1;
	  while (uiptr2 < uiptr) {
	    ujj = *uiptr2++;
	    if (ujj != uii) {
	      uiptr3[batch_ct] = 1;
	      batch_sizes[batch_ct++] = ujj;
	      uii = ujj;
	    } else {
	      uiptr3[batch_ct - 1] += 1;
	    }
	  }
	  break;
	}
      }

      // batch numbers
      if (bigstack_end_alloc_ui(batch_ct, &uiptr)) {
	goto plink1_dosage_ret_NOMEM;
      }
      memcpy(uiptr, batch_sizes, batch_ct * sizeof(int32_t));
      bigstack_shrink_top(batch_sizes, batch_ct * sizeof(int32_t));
      memcpy(batch_sizes, uiptr3, batch_ct * sizeof(int32_t));
      // convert uiptr3 to write offset array
      uii = uiptr3[0];
      uiptr3[0] = 0;
      for (ujj = 1; ujj < batch_ct; ujj++) {
	ukk = uii;
	uii += uiptr3[ujj];
        uiptr3[ujj] = ukk;
      }
    } else {
      uiptr3 = batch_sizes;
      uiptr3[0] = 0;
    }
    if (bigstack_alloc_c(infile_ct * max_fn_len, &fnames)) {
      goto plink1_dosage_ret_NOMEM;
    }
    if (sepheader) {
      if (bigstack_alloc_c(infile_ct * max_sepheader_len, &sep_fnames)) {
	goto plink1_dosage_ret_NOMEM;
      }
    }
    bigstack_end_reset(bigstack_end_mark);
    rewind(infile);
    while (fgets(g_textbuf, MAXLINELEN, infile)) {
      bufptr = skip_initial_spaces(g_textbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      if (batch_ct) {
	scan_int32(bufptr, &ii);
	uii = int32arr_greater_than((int32_t*)uiptr, batch_ct, ii);
	file_idx = uiptr3[uii];
	uiptr3[uii] += 1;
	bufptr = next_token(bufptr);
      } else {
        file_idx = uiptr3[0];
	uiptr3[0] += 1;
      }
      bufptr2 = token_endnn(bufptr);
      memcpyx(&(fnames[file_idx * max_fn_len]), bufptr, (uintptr_t)(bufptr2 - bufptr), '\0');
      if (sepheader) {
	bufptr = skip_initial_spaces(bufptr2);
	bufptr2 = token_endnn(bufptr);
	memcpyx(&(sep_fnames[file_idx * max_sepheader_len]), bufptr, (uintptr_t)(bufptr2 - bufptr), '\0');
      }
    }
    if (fclose_null(&infile)) {
      goto plink1_dosage_ret_READ_FAIL;
    }
    if (batch_ct) {
      for (uii = 0; uii < batch_ct; uii++) {
	if (batch_sizes[uii] > max_batch_size) {
	  max_batch_size = batch_sizes[uii];
	}
      }
    } else {
      batch_ct = 1;
      max_batch_size = infile_ct;
    }
  } else {
    uii = strlen(doip->fname) + 1;
    if (bigstack_alloc_ui(1, &batch_sizes) ||
        bigstack_alloc_c(uii + 1, &fnames)) {
      goto plink1_dosage_ret_NOMEM;
    }
    batch_sizes[0] = 1;
    memcpy(fnames, doip->fname, uii);
    infile_ct = 1;
  }
  if (bigstack_alloc_uc(2 * PIGZ_BLOCK_SIZE, &overflow_buf) ||
      bigstack_alloc_ui(max_batch_size, &file_icts) ||
      bigstack_alloc_ul(max_batch_size, &line_idx_arr) ||
      bigstack_alloc_ul(sample_ctl, &batch_samples) ||
      bigstack_alloc_ul(sample_ctl, &cur_samples) ||
      bigstack_alloc_ui(sample_ct, &read_idx_to_sample_idx) ||
      bigstack_alloc_ui(sample_ct, &skip_vals) ||
      bigstack_alloc_d(sample_ct, &cur_dosages) ||
      bigstack_alloc_c(MAX_ID_SLEN, &cur_marker_id_buf)) {
    goto plink1_dosage_ret_NOMEM;
  }
  gz_infiles = (gzFile*)bigstack_alloc(infile_ct * sizeof(gzFile));
  if (!gz_infiles) {
    infile_ct = 0;
    goto plink1_dosage_ret_NOMEM;
  }
  for (uii = 0; uii < infile_ct; uii++) {
    gz_infiles[uii] = nullptr;
  }
  if (noheader) {
    if (infile_ct != 1) {
      for (batch_idx = 0; batch_idx < batch_ct; batch_idx++) {
	if (batch_sizes[batch_idx] != 1) {
	  logerrprint("Error: --dosage 'noheader' modifier cannot be used with multifile batches.\n");
	  goto plink1_dosage_ret_INVALID_CMDLINE;
	}
      }
    }
    // sorted_sample_ids = nullptr;
  } else {
    retval = sort_item_ids(unfiltered_sample_ct, sample_exclude, sample_exclude_ct, sample_ids, max_sample_id_len, 0, 1, strcmp_deref, &sorted_sample_ids, &sample_id_map);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  }
  if (freq_cc) {
    if (!do_glm) {
      logerrprint("Error: --dosage 'case-control-freqs' modifier can only be used in association\nanalysis mode.\n");
      goto plink1_dosage_ret_INVALID_CMDLINE;
    } else if (!pheno_c) {
      logerrprint("Warning: '--dosage case-control-freqs' is silly with a quantitative phenotype.\n");
    }
  }
  if (do_glm) {
    if (!sample_exclude_ct) {
      pheno_nm_collapsed = pheno_nm;
#ifndef NOLAPACK
      pheno_d_collapsed = pheno_d;
#endif
      pheno_c_collapsed = pheno_c;
    } else {
      if (bigstack_alloc_ul(sample_ctl, &pheno_nm_collapsed)) {
	goto plink1_dosage_ret_NOMEM;
      }
      copy_bitarr_subset_excl(pheno_nm, sample_exclude, unfiltered_sample_ct, sample_ct, pheno_nm_collapsed);
#ifndef NOLAPACK
      if (pheno_d) {
	pheno_d_collapsed = (double*)alloc_and_init_collapsed_arr((char*)pheno_d, sizeof(double), unfiltered_sample_ct, sample_exclude, sample_ct, 0);
	if (!pheno_d_collapsed) {
	  goto plink1_dosage_ret_NOMEM;
	}
      } else {
#endif
	if (bigstack_alloc_ul(sample_ctl, &pheno_c_collapsed)) {
	  goto plink1_dosage_ret_NOMEM;
	}
	copy_bitarr_subset_excl(pheno_c, sample_exclude, unfiltered_sample_ct, sample_ct, pheno_c_collapsed);
#ifndef NOLAPACK
      }
#endif
    }
#ifndef NOLAPACK
    if (pheno_d) {
      if (bigstack_alloc_d(sample_ct, &pheno_d2) ||
	  bigstack_alloc_d(param_ct * sample_ct, &covars_cov_major_buf) ||
	  bigstack_alloc_d(param_ct * sample_ct, &covars_sample_major_buf) ||
	  bigstack_alloc_d(param_ct * param_ct, &param_2d_buf) ||
	  bigstack_alloc_d(param_ct * param_ct, &param_2d_buf2) ||
	  bigstack_alloc_d(param_ct - 1, &regression_results) ||
	  bigstack_alloc_d(param_ct * sample_ct, &dgels_a) ||
	  bigstack_alloc_d(sample_ct, &dgels_b)) {
	goto plink1_dosage_ret_NOMEM;
      }
      mi_buf = (MATRIX_INVERT_BUF1_TYPE*)bigstack_alloc(param_ct * MATRIX_INVERT_BUF1_ELEM_ALLOC);
      if (!mi_buf) {
	goto plink1_dosage_ret_NOMEM;
      }
      // workspace query
      dgels_m = (int32_t)((uint32_t)sample_ct);
      dgels_n = (int32_t)((uint32_t)param_ct);
      dgels_nrhs = 1;
      dgels_ldb = dgels_m;
      dgels_lwork = -1;
      dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, &dxx, &dgels_lwork, &dgels_info);
      if (dxx > 2147483647.0) {
	logerrprint("Error: Multiple linear regression problem too large for current LAPACK version.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto plink1_dosage_ret_1;
      }
      dgels_lwork = (int32_t)dxx;
      if (bigstack_alloc_d(dgels_lwork, &dgels_work)) {
	goto plink1_dosage_ret_NOMEM;
      }
    } else {
#endif
      if (covar_ct) {
	ulii = covar_ct * sample_ct;
	if (bigstack_alloc_f(ulii, &covar_f)) {
	  goto plink1_dosage_ret_NOMEM;
	}
	fptr = covar_f;
	dptr = covar_d;
	for (uljj = 0; uljj < ulii; uljj++) {
	  *fptr++ = (float)(*dptr++);
	}
      }
      if (bigstack_alloc_ul(BITCT_TO_ALIGNED_WORDCT(sample_ct), &perm_vec) ||
	  bigstack_alloc_f(param_ct * sample_cta4, &covars_cov_major_f_buf) ||
	  bigstack_alloc_f(param_cta4, &coef_f) ||
	  bigstack_alloc_f(sample_cta4, &pp_f) ||
	  bigstack_alloc_f(sample_ct, &sample_1d_buf_f) ||
	  bigstack_alloc_f(sample_ct, &pheno_buf_f) ||
	  bigstack_alloc_f(param_ct, &param_1d_buf_f) ||
	  bigstack_alloc_f(param_ct, &param_1d_buf2_f) ||
	  bigstack_alloc_f(param_ct * param_cta4, &param_2d_buf_f) ||
	  bigstack_alloc_f(param_ct * param_cta4, &param_2d_buf2_f) ||
	  bigstack_alloc_f(param_ct - 1, &regression_results_f)) {
	goto plink1_dosage_ret_NOMEM;
      }
#ifndef NOLAPACK
    }
#endif
    if (bigstack_alloc_ul(1, &perm_fails)) {
      goto plink1_dosage_ret_NOMEM;
    }
    if (load_map) {
      bufptr = memcpya(g_textbuf, " CHR         SNP          BP", 28);
    } else {
      bufptr = memcpya(g_textbuf, "         SNP", 12);
    }
    bufptr = memcpya(bufptr, freq_cc? "  A1  A2   FRQ_A   FRQ_U    INFO    " : "  A1  A2     FRQ    INFO    ", freq_cc? 36 : 28);
    bufptr = memcpya(bufptr, pheno_c? "  OR" : "BETA", 4);
    bufptr = memcpya(bufptr, "      SE       P", 16);
    append_binary_eoln(&bufptr);
    bufptr2 = memcpyb(outname_end, ".assoc.dosage", 14);
  } else if (count_occur) {
    // could just use a uint32_t array if .map provided
    htable = (Ll_ctstr_entry**)bigstack_alloc(HASHSIZE * sizeof(intptr_t));
    if (!htable) {
      goto plink1_dosage_ret_NOMEM;
    }
    for (uii = 0; uii < HASHSIZE; uii++) {
      htable[uii] = nullptr;
    }
    bufptr2 = memcpyb(outname_end, ".occur.dosage", 14);
  } else if (!do_score) {
    if (format_val != 1) {
      if (bigstack_alloc_d(sample_ct, &cur_dosages2)) {
	goto plink1_dosage_ret_NOMEM;
      }
    }
    bufptr2 = memcpyb(outname_end, ".out.dosage", 12);
  }

  if (!do_score) {
    if (output_gz) {
      memcpy(bufptr2, ".gz", 4);
    }
    if (flex_pzwrite_init(output_gz, outname, overflow_buf, 0, &ps)) {
      goto plink1_dosage_ret_OPEN_FAIL;
    }
    pzwritep = (char*)overflow_buf;
    if (do_glm) {
      pzwritep = memcpya(pzwritep, g_textbuf, bufptr - g_textbuf);
    } else if (!count_occur) {
      pzwritep = memcpya(pzwritep, "SNP A1 A2 ", 10);
      for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_unsafe_ck(sample_exclude, &sample_uidx);
	bufptr = &(sample_ids[sample_uidx * max_sample_id_len]);
	bufptr2 = strchr(bufptr, '\t');
	pzwritep = memcpya(pzwritep, bufptr, bufptr2 - bufptr);
	*pzwritep++ = ' ';
	pzwritep = strcpyax(pzwritep, &(bufptr2[1]), ' ');
	if (flex_pzwrite(&ps, &pzwritep)) {
	  goto plink1_dosage_ret_WRITE_FAIL;
	}
      }
      append_binary_eoln(&pzwritep);
    }
  }
  bigstack_mark = g_bigstack_base;
  for (batch_idx = 0; batch_idx < batch_ct; batch_idx++, file_idx_start += cur_batch_size) {
    cur_batch_size = batch_sizes[batch_idx];
    read_idx = 0;
    loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto plink1_dosage_ret_NOMEM;
    }
    loadbuf = (char*)g_bigstack_base;
    loadbuf[loadbuf_size - 1] = ' ';
    fill_ulong_zero(sample_ctl, batch_samples);
    bufptr = memcpya(g_logbuf, "--dosage: Reading from ", 23);
    if (cur_batch_size == 1) {
      bufptr = strcpya(bufptr, &(fnames[file_idx_start * max_fn_len]));
    } else if (cur_batch_size == 2) {
      bufptr = strcpya(bufptr, &(fnames[file_idx_start * max_fn_len]));
      bufptr = memcpya(bufptr, " and ", 5);
      bufptr = strcpya(bufptr, &(fnames[(file_idx_start + 1) * max_fn_len]));
    } else {
      for (file_idx = 0; file_idx < cur_batch_size - 1; file_idx++) {
        bufptr = strcpya(bufptr, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	bufptr = memcpya(bufptr, ", ", 2);
      }
      bufptr = memcpya(bufptr, "and ", 4);
      bufptr = strcpya(bufptr, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
    }
    memcpyl3(bufptr, ".\n");
    wordwrapb(0);
    logprintb();
    for (file_idx = 0; file_idx < cur_batch_size; file_idx++) {
      read_idx_start = read_idx;
      if (sepheader) {
	retval = gzopen_read_checked(&(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]), &(gz_infiles[file_idx]));
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
	line_idx = 0;
	uii = 1 + skip2; // current skip value
	while (gzgets(gz_infiles[file_idx], g_textbuf, MAXLINELEN)) {
	  line_idx++;
	  if (!g_textbuf[MAXLINELEN - 1]) {
	    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
          bufptr = skip_initial_spaces(g_textbuf);
          if (is_eoln_kns(*bufptr)) {
	    continue;
	  }
          if (bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, &bufptr2, &ii, &(g_textbuf[MAXLINELEN]))) {
            sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if (ii == -1) {
	    uii += format_val;
	  } else {
	    ii = sample_id_map[(uint32_t)ii];
	    if (is_set(batch_samples, ii)) {
	      bufptr = &(sorted_sample_ids[((uint32_t)ii) * max_sample_id_len]);
	      *strchr(bufptr, '\t') = ' ';
	      sprintf(g_logbuf, "Error: '%s' appears multiple times.\n", bufptr);
	      goto plink1_dosage_ret_INVALID_FORMAT_WW;
	    }
	    set_bit(ii, batch_samples);
	    read_idx_to_sample_idx[read_idx] = (uint32_t)ii;
            skip_vals[read_idx++] = uii;
	    uii = 1 + (format_val == 3);
	  }
	}
	if (gzclose(gz_infiles[file_idx]) != Z_OK) {
	  gz_infiles[file_idx] = nullptr;
	  goto plink1_dosage_ret_READ_FAIL;
	}
	gz_infiles[file_idx] = nullptr;
	if (read_idx_start == read_idx) {
          sprintf(g_logbuf, "Error: %s is empty.\n", &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]));
          goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
      }
      retval = gzopen_read_checked(&(fnames[(file_idx + file_idx_start) * max_fn_len]), &(gz_infiles[file_idx]));
      if (retval) {
	goto plink1_dosage_ret_OPEN_FAIL;
      }
      line_idx = 0;
      if (noheader) {
	sample_uidx = 0;
	next_unset_unsafe_ck(sample_exclude, &sample_uidx);
	skip_vals[0] = 1 + skip2 + format_val * sample_uidx;
	uii = 1 + (format_val == 3);
	for (read_idx = 1; read_idx < sample_ct; read_idx++) {
	  prev_sample_uidx = sample_uidx++;
	  next_unset_unsafe_ck(sample_exclude, &sample_uidx);
	  skip_vals[read_idx] = uii + format_val * (sample_uidx - prev_sample_uidx - 1);
	}
	// bugfix: this is a read_idx -> sample_idx map, not read_idx ->
	// sample_uidx
	for (read_idx = 0; read_idx < sample_ct; read_idx++) {
	  read_idx_to_sample_idx[read_idx] = read_idx;
	}
	fill_all_bits(sample_ct, batch_samples);
      } else if (!sepheader) {
	do {
	  if (!gzgets(gz_infiles[file_idx], loadbuf, loadbuf_size)) {
            sprintf(g_logbuf, "Error: %s is empty.\n", &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  line_idx++;
	  if (!loadbuf[loadbuf_size - 1]) {
	    goto plink1_dosage_ret_LONG_LINE;
	  }
	  bufptr = skip_initial_spaces(loadbuf);
	} while (is_eoln_kns(*bufptr));
	bufptr = next_token_multz(bufptr, skip0);
	bufptr2 = next_token_mult(bufptr, skip1p1);
	if (no_more_tokens_kns(bufptr2)) {
	  goto plink1_dosage_ret_MISSING_TOKENS;
	}
	if (strcmp_se(bufptr, "SNP", 3)) {
	  sprintf(g_logbuf, "Error: Column %u of %s's header isn't 'SNP'.\n", skip0 + 1, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	} else if (strcmp_se(bufptr2, "A1", 2)) {
	  sprintf(g_logbuf, "Error: Column %u of %s's header isn't 'A1'.\n", skip0 + skip1p1 + 1, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
	bufptr = next_token(bufptr2);
	bufptr2 = next_token_multz(bufptr, skip2);
	if (no_more_tokens_kns(bufptr2)) {
	  goto plink1_dosage_ret_MISSING_TOKENS;
	}
	if (strcmp_se(bufptr, "A2", 2)) {
	  sprintf(g_logbuf, "Error: Column %u of %s's header isn't 'A2'.\n", skip0 + skip1p1 + 2, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
	uii = 1 + skip2;
	bufptr = skip_initial_spaces(token_endnn(bufptr2));
	while (!is_eoln_kns(*bufptr)) {
          if (bsearch_read_fam_indiv(bufptr, sorted_sample_ids, max_sample_id_len, sample_ct, &bufptr2, &ii, g_textbuf)) {
	    sprintf(g_logbuf, "Error: Header of %s has an odd number of tokens in the FID/IID section.\n", &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if (ii == -1) {
	    uii += format_val;
	  } else {
	    ii = sample_id_map[(uint32_t)ii];
	    if (is_set(batch_samples, ii)) {
	      bufptr = &(sorted_sample_ids[((uint32_t)ii) * max_sample_id_len]);
	      *strchr(bufptr, '\t') = ' ';
	      sprintf(g_logbuf, "Error: '%s' appears multiple times.\n", bufptr);
	      goto plink1_dosage_ret_INVALID_FORMAT_WW;
	    }
	    set_bit(ii, batch_samples);
	    read_idx_to_sample_idx[read_idx] = (uint32_t)ii;
            skip_vals[read_idx++] = uii;
	    uii = 1 + (format_val == 3);
	  }
	  bufptr = bufptr2;
	}
	if (read_idx_start == read_idx) {
	  sprintf(g_logbuf, "Error: Header of %s has no tokens in the FID/IID section.\n", &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
      }
      file_icts[file_idx] = read_idx - read_idx_start;
      line_idx_arr[file_idx] = line_idx;
    }
    if ((batch_ct == 1) && (!noheader)) {
      ulii = sample_ct - popcount_longs(batch_samples, sample_ctl);
      if (ulii) {
	LOGERRPRINTFWW("Warning: %" PRIuPTR " sample ID%s present in .fam file but missing from dosage file%s.\n", ulii, (ulii == 1)? "" : "s", (cur_batch_size == 1)? "" : "s");
      }
    }

    while (1) {
      read_idx_start = 0;
      memcpy(cur_samples, batch_samples, sample_ctl * sizeof(intptr_t));
      for (file_idx = 0; file_idx < cur_batch_size; file_idx++) {
	line_idx = line_idx_arr[file_idx];
	do {
	  if (!gzgets(gz_infiles[file_idx], loadbuf, loadbuf_size)) {
	    if (file_idx) {
	      logerrprint("Error: Misaligned dosage data files.\n");
	      goto plink1_dosage_ret_INVALID_FORMAT;
	    }
	    goto plink1_dosage_end_loop;
	  }
	  line_idx++;
	  if (!loadbuf[loadbuf_size - 1]) {
	    goto plink1_dosage_ret_LONG_LINE;
	  }
          bufptr = skip_initial_spaces(loadbuf);
	} while (is_eoln_kns(*bufptr));
	line_idx_arr[file_idx] = line_idx;
	bufptr = next_token_multz(bufptr, skip0);
	bufptr3 = next_token_mult(bufptr, skip1p1);
	bufptr5 = next_token(bufptr3);
	if (no_more_tokens_kns(bufptr5)) {
	  goto plink1_dosage_ret_MISSING_TOKENS;
	}
	bufptr2 = token_endnn(bufptr);
	bufptr4 = token_endnn(bufptr3);
	bufptr6 = token_endnn(bufptr5);
        slen = (uintptr_t)(bufptr2 - bufptr);
	if (slen > MAX_ID_SLEN) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has an excessively long variant ID.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
	if (!file_idx) {
	  memcpyx(cur_marker_id_buf, bufptr, slen, '\0');
	  cur_marker_id_len = slen;
	  a1_len = (uintptr_t)(bufptr4 - bufptr3);
	  if (allele_set(bufptr3, a1_len, &a1_ptr)) {
	    goto plink1_dosage_ret_NOMEM;
	  }
	  a2_len = (uintptr_t)(bufptr6 - bufptr5);
	  if (allele_set(bufptr5, a2_len, &a2_ptr)) {
	    goto plink1_dosage_ret_NOMEM;
	  }
	  if (load_map) {
	    marker_idx = id_htable_find(bufptr, slen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_blen);
	    if (marker_idx == 0xffffffffU) {
#ifdef __LP64__
	      marker_idx = ~ZEROLU;
#endif
	      continue;
	    }
	  }
	  if (do_score) {
            score_marker_idx = (uint32_t)bsearch_str(bufptr, slen, score_marker_ids, max_score_marker_id_len, score_marker_ct);
	    if (score_marker_idx == 0xffffffffU) {
	      continue;
	    }
	    if (score_qrange_key_exists && (!is_set(score_qrange_key_exists, score_marker_idx))) {
	      score_marker_idx = 0xffffffffU;
	      continue;
	    }
	    if (!strcmp(a1_ptr, score_allele_codes[score_marker_idx])) {
	      score_a2_effect = 0;
	    } else if (!strcmp(a2_ptr, score_allele_codes[score_marker_idx])) {
	      score_a2_effect = 1;
	    } else {
	      // PLINK 1.07 just skips the marker instead of erroring out here,
	      // so we replicate that
	      score_marker_idx = 0xffffffffU;
	      continue;
	    }
	  }
	} else {
	  if ((slen != cur_marker_id_len) || memcmp(bufptr, cur_marker_id_buf, slen)) {
	    sprintf(g_logbuf, "Error: Variant ID mismatch between line %" PRIuPTR " of %s and line %" PRIuPTR " of %s.\n", line_idx_arr[0], &(fnames[file_idx_start * max_fn_len]), line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if (((uintptr_t)(bufptr4 - bufptr3) != a1_len) || memcmp(bufptr3, a1_ptr, a1_len) || ((uintptr_t)(bufptr6 - bufptr5) != a2_len) || memcmp(bufptr5, a2_ptr, a2_len)) {
	    sprintf(g_logbuf, "Error: Allele code mismatch between line %" PRIuPTR " of %s and line %" PRIuPTR " of %s.\n", line_idx_arr[0], &(fnames[file_idx_start * max_fn_len]), line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if ((marker_idx == ~ZEROLU) || (score_marker_idx == 0xffffffffU)) {
	    continue;
	  }
	}

	if (count_occur) {
          uii = hashval2(cur_marker_id_buf, slen++);
	  ll_pptr = &(htable[uii]);
	  while (1) {
	    ll_ptr = *ll_pptr;
	    if (!ll_ptr) {
	      distinct_id_ct++;
	      ll_ptr = (Ll_ctstr_entry*)bigstack_end_alloc(slen + sizeof(Ll_ctstr_entry));
	      loadbuf_size = bigstack_left();
	      if (loadbuf_size >= MAXLINEBUFLEN) {
		loadbuf_size = MAXLINEBUFLEN;
	      } else if (loadbuf_size > MAXLINELEN) {
                loadbuf[loadbuf_size - 1] = ' ';
	      } else {
		goto plink1_dosage_ret_NOMEM;
	      }
	      ll_ptr->next = nullptr;
	      memcpy(ll_ptr->ss, cur_marker_id_buf, slen);
	      if (slen > max_occur_id_len) {
		max_occur_id_len = slen;
	      }
	      ll_ptr->ct = 1;
	      *ll_pptr = ll_ptr;
	      break;
	    }
            if (!strcmp(ll_ptr->ss, cur_marker_id_buf)) {
	      ll_ptr->ct += 1;
	      break;
	    }
	    ll_pptr = &(ll_ptr->next);
	  }
	} else {
	  *bufptr6 = ' ';
	  bufptr = bufptr5;
	  read_idx = read_idx_start + file_icts[file_idx];
	  if (format_val == 1) {
	    for (; read_idx_start < read_idx; read_idx_start++) {
	      bufptr = next_token_mult(bufptr, skip_vals[read_idx_start]);
	      if (!bufptr) {
		goto plink1_dosage_ret_MISSING_TOKENS;
	      }
	      if (scan_double(bufptr, &dxx)) {
		clear_bit(read_idx_to_sample_idx[read_idx_start], cur_samples);
		continue;
	      }
	      if (!dose1) {
		dxx *= 0.5;
	      }
	      if ((dxx > 1.0 + DOSAGE_EPSILON) || (dxx < 0.0)) {
		clear_bit(read_idx_to_sample_idx[read_idx_start], cur_samples);
		continue;
	      } else if (dxx > 1.0) {
		dxx = 1.0;
	      }
	      cur_dosages[read_idx_to_sample_idx[read_idx_start]] = dxx;
	    }
	  } else {
	    for (; read_idx_start < read_idx; read_idx_start++) {
	      bufptr2 = next_token_mult(bufptr, skip_vals[read_idx_start]);
	      bufptr = next_token(bufptr2);
	      if (no_more_tokens_kns(bufptr)) {
		goto plink1_dosage_ret_MISSING_TOKENS;
	      }
	      if (scan_double(bufptr2, &dxx) || scan_double(bufptr, &dyy)) {
		clear_bit(read_idx_to_sample_idx[read_idx_start], cur_samples);
		continue;
	      }
	      dzz = dxx + dyy;
	      if ((dyy < 0.0) || (dxx < 0.0) || (dzz > 1.0 + DOSAGE_EPSILON)) {
		clear_bit(read_idx_to_sample_idx[read_idx_start], cur_samples);
		continue;
	      } else if (dzz > 1.0) {
		dzz = 1.0 / dzz;
		dxx *= dzz;
		dyy *= dzz;
	      }
	      uii = read_idx_to_sample_idx[read_idx_start];
	      if (!cur_dosages2) {
		dxx += dyy * 0.5;
		cur_dosages[uii] = dxx;
	      } else {
		cur_dosages[uii] = dxx;
		cur_dosages2[uii] = dyy;
	      }
	    }
	  }
	}
      }

      if ((marker_idx != ~ZEROLU) && (score_marker_idx != 0xffffffffU)) {
	if (do_glm) {
	  if (covar_nm) {
	    // it would be more efficient to make covar_nm act as a mask on
	    // sample_exclude earlier, but this is throwaway code
	    bitvec_and(covar_nm, sample_ctl, cur_samples);
	  }
	  sample_valid_ct = popcount_longs(cur_samples, sample_ctl);
	  dxx = 0.0;
	  dyy = 0.0;
	  for (sample_idx = 0, sample_uidx = 0; sample_idx < sample_valid_ct; sample_idx++, sample_uidx++) {
            next_set_unsafe_ck(cur_samples, &sample_uidx);
	    dzz = cur_dosages[sample_uidx];
            dxx += dzz;
	    dyy += dzz * dzz;
	  }
	  sample_valid_ct_recip = 1.0 / ((double)((intptr_t)sample_valid_ct));
	  dzz = dxx * sample_valid_ct_recip;
	  // printf("sample_valid_ct: %lu  dzz: %g\n", sample_valid_ct, dzz);
	  // dosageSSQ = sum(x^2 - bar{x}^2) = sum(x^2) - sum(x) * avg(x)
	  // theoreticalVariance = avg(x) * (1 - avg(x))
	  // empiricalVariance = 2 * (dosageSSQ / cnt)
	  // rsq = (theoretical_variance > 0)? (empirical / theoretical) : 0

	  // dzz is minor allele frequency, dyy is ssq
	  dyy -= dxx * dzz;

          dxx = dzz * (1.0 - dzz); // now dxx = theoretical var
	  dyy = 2 * dyy * sample_valid_ct_recip; // and dyy = empirical
	  rsq = (dxx > 0.0)? (dyy / dxx) : 0.0;
	  if ((dxx >= 0.0098999999999999) && (rsq >= 0.1) && (rsq <= 2.0)) {
#ifndef NOLAPACK
	    if (pheno_d) {
	      is_valid = glm_linear_dosage(sample_ct, cur_samples, sample_valid_ct, pheno_nm_collapsed, pheno_d_collapsed, perm_fails, covar_ct, covar_d, cur_dosages, pheno_d2, covars_cov_major_buf, covars_sample_major_buf, param_2d_buf, mi_buf, param_2d_buf2, regression_results, dgels_a, dgels_b, dgels_work, dgels_lwork, standard_beta, glm_vif_thresh, &beta, &se, &pval);
	      if (is_valid == 2) {
		// NOMEM special case
		goto plink1_dosage_ret_NOMEM;
	      }
	    } else {
#endif
	      is_valid = glm_logistic_dosage(sample_ct, cur_samples, sample_valid_ct, pheno_nm_collapsed, pheno_c_collapsed, perm_vec, perm_fails, covar_ct, covar_f, cur_dosages, coef_f, pp_f, pheno_buf_f, covars_cov_major_f_buf, param_1d_buf_f, param_1d_buf2_f, param_2d_buf_f, param_2d_buf2_f, regression_results_f, sample_1d_buf_f, &beta, &se, &pval);
#ifndef NOLAPACK
	    }
#endif
	  } else {
	    is_valid = 0;
	    if (rsq > 2.0) {
	      rsq = 2.0;
	    }
	  }
	  if (load_map) {
	    pzwritep = width_force(4, pzwritep, chrom_name_write(chrom_info_ptr, get_variant_chrom(chrom_info_ptr, marker_idx), pzwritep));
	    *pzwritep++ = ' ';
	    pzwritep = fw_strcpyn(11, cur_marker_id_len, cur_marker_id_buf, pzwritep);
            pzwritep = memseta(pzwritep, 32, 2);
            pzwritep = uint32toa_w10(marker_pos[marker_idx], pzwritep);
	  } else {
	    *pzwritep++ = ' ';
	    pzwritep = fw_strcpyn(11, cur_marker_id_len, cur_marker_id_buf, pzwritep);
	  }
	  *pzwritep++ = ' ';
	  if (a1_len < 3) {
	    *pzwritep++ = ' ';
	    if (a1_len == 1) {
	      *pzwritep++ = ' ';
	    }
	  }
	  if (flex_pzputs_allele(&ps, &pzwritep, a1_ptr, a1_len)) {
	    goto plink1_dosage_ret_WRITE_FAIL;
	  }
	  *pzwritep++ = ' ';
	  if (a2_len < 3) {
            *pzwritep++ = ' ';
	    if (a2_len == 1) {
	      *pzwritep++ = ' ';
	    }
	  }
	  if (flex_pzputs_allele(&ps, &pzwritep, a2_ptr, a2_len)) {
	    goto plink1_dosage_ret_WRITE_FAIL;
	  }
	  *pzwritep++ = ' ';
	  if (freq_cc && pheno_c) {
	    // compute case and control A1 frequencies
	    dxx = 0.0; // case sum
	    dyy = 0.0; // control sum
	    uii = 0;
	    for (sample_idx = 0, sample_uidx = 0; sample_idx < sample_valid_ct; sample_idx++, sample_uidx++) {
	      next_set_unsafe_ck(cur_samples, &sample_uidx);
	      dzz = cur_dosages[sample_uidx];
	      if (is_set(pheno_c_collapsed, sample_uidx)) {
		dxx += dzz;
		uii++; // could also use a popcount
	      } else {
		dyy += dzz;
	      }
	    }
	    if (uii) {
	      pzwritep = dtoa_f_w7p4x(dxx / ((double)((int32_t)uii)), ' ', pzwritep);
	    } else {
	      pzwritep = memcpya(pzwritep, "     NA ", 8);
	    }
	    uii = sample_valid_ct - uii;
	    if (uii) {
	      pzwritep = dtoa_f_w7p4(dyy / ((double)((int32_t)uii)), pzwritep);
	    } else {
	      pzwritep = memcpya(pzwritep, "     NA", 7);
	    }
	  } else {
            pzwritep = dtoa_f_w7p4(dzz, pzwritep);
	    // remove this kludge once scripts stop depending on it
	    if (freq_cc) {
	      *pzwritep++ = ' ';
	      pzwritep = dtoa_f_w7p4(dzz, pzwritep);
	    }
	  }
	  *pzwritep++ = ' ';
	  pzwritep = dtoa_f_w7p4x(rsq, ' ', pzwritep);
	  if (is_valid) {
	    pzwritep = dtoa_f_w7p4x(pheno_c? exp(beta * 0.5) : (beta * 0.5), ' ', pzwritep);
	    pzwritep = dtoa_f_w7p4x(se * 0.5, ' ', pzwritep);
	    pzwritep = dtoa_g_wxp4(MAXV(pval, output_min_p), 7, pzwritep);
	  } else {
	    pzwritep = memcpya(pzwritep, "     NA      NA      NA", 23);
	  }
	  append_binary_eoln(&pzwritep);
	  if (flex_pzwrite(&ps, &pzwritep)) {
	    goto plink1_dosage_ret_WRITE_FAIL;
	  }
	} else if (do_score) {
	  sample_valid_ct = popcount_longs(cur_samples, sample_ctl);
	  qrange_idx = 0;
	  score_cur_effect_size = score_effect_sizes[score_marker_idx];
	  if (score_a2_effect) {
	    for (sample_idx = 0, sample_uidx = 0; sample_idx < sample_valid_ct; sample_idx++, sample_uidx++) {
              next_set_unsafe_ck(cur_samples, &sample_uidx);
	      cur_dosages[sample_uidx] = 1.0 - cur_dosages[sample_uidx];
	    }
	  }
	  dxx = 0.0; // dosage sum
	  for (sample_idx = 0, sample_uidx = 0; sample_idx < sample_valid_ct; sample_idx++, sample_uidx++) {
	    next_set_unsafe_ck(cur_samples, &sample_uidx);
	    dxx += cur_dosages[sample_uidx];
	  }
	  if (sample_valid_ct) {
            dxx /= (double)((intptr_t)sample_valid_ct); // now average
	  } else {
	    dxx = 0.0;
	  }
	  score_missing_effect = dxx * score_cur_effect_size;
	  do {
	    if (qrange_ct) {
	      if ((score_qrange_bounds[2 * qrange_idx] > score_qrange_keys[score_marker_idx]) || (score_qrange_bounds[2 * qrange_idx + 1] < score_qrange_keys[score_marker_idx])) {
		continue;
	      }
	    }
	    score_range_obs_cts[qrange_idx] += 1;
	    dptr = &(cur_scores[sample_ct * qrange_idx]);
	    if (score_center) {
	      score_bases[qrange_idx] -= score_missing_effect;
	    }
	    uiptr = &(score_miss_cts[sample_ct * qrange_idx]);
	    for (sample_uidx = 0; sample_uidx < sample_ct; sample_uidx++) {
	      if (!is_set(cur_samples, sample_uidx)) {
		uiptr[sample_uidx] += 1;
		if (score_mean_impute) {
		  dptr[sample_uidx] += score_missing_effect;
		}
	      } else {
		dptr[sample_uidx] += cur_dosages[sample_uidx] * score_cur_effect_size;
	      }
	    }
	  } while (++qrange_idx < qrange_ct);
	} else if (!count_occur) {
	  // --write-dosage
	  pzwritep = strcpyax(pzwritep, cur_marker_id_buf, ' ');
	  if (flex_pzputs_allele(&ps, &pzwritep, a1_ptr, a1_len)) {
	    goto plink1_dosage_ret_WRITE_FAIL;
	  }
	  *pzwritep++ = ' ';
	  if (flex_pzputs_allele(&ps, &pzwritep, a2_ptr, a2_len)) {
	    goto plink1_dosage_ret_WRITE_FAIL;
	  }
	  *pzwritep++ = ' ';
	  if (flex_pzwrite(&ps, &pzwritep)) {
	    goto plink1_dosage_ret_WRITE_FAIL;
	  }
	  ulii = 0;
	  // could make output format independent of input format (other than
	  // outformat > 1 can't coexist with format == 1), but kinda pointless
	  // since the code won't be kept for PLINK 2.0.
	  if (format_val == 1) {
	    // bugfix: don't keep resetting sample_idx to 0
	    sample_idx = 0;
	    do {
	      ulii += MAXLINELEN / 16;
	      if (ulii > sample_ct) {
		ulii = sample_ct;
	      }
	      for (; sample_idx < ulii; sample_idx++) {
		if (!is_set(cur_samples, sample_idx)) {
		  pzwritep = memcpyl3a(pzwritep, "NA ");
		} else {
		  pzwritep = dtoa_gx(2 * cur_dosages[sample_idx], ' ', pzwritep);
		}
	      }
	      if (flex_pzwrite(&ps, &pzwritep)) {
		goto plink1_dosage_ret_WRITE_FAIL;
	      }
	    } while (ulii < sample_ct);
	  } else if (format_val == 2) {
	    sample_idx = 0;
	    do {
	      ulii += MAXLINELEN / 32;
	      if (ulii > sample_ct) {
		ulii = sample_ct;
	      }
	      for (; sample_idx < ulii; sample_idx++) {
		if (!is_set(cur_samples, sample_idx)) {
		  pzwritep = memcpya(pzwritep, "NA NA ", 6);
		} else {
		  pzwritep = dtoa_gx(cur_dosages[sample_idx], ' ', pzwritep);
		  pzwritep = dtoa_gx(cur_dosages2[sample_idx], ' ', pzwritep);
		}
	      }
	      if (flex_pzwrite(&ps, &pzwritep)) {
		goto plink1_dosage_ret_WRITE_FAIL;
	      }
	    } while (ulii < sample_ct);
	  } else {
	    sample_idx = 0;
	    do {
	      ulii += MAXLINELEN / 48;
	      if (ulii > sample_ct) {
		ulii = sample_ct;
	      }
	      for (; sample_idx < ulii; sample_idx++) {
		if (!is_set(cur_samples, sample_idx)) {
		  pzwritep = memcpya(pzwritep, "NA NA NA ", 9);
		} else {
		  dxx = cur_dosages[sample_idx];
		  pzwritep = dtoa_gx(dxx, ' ', pzwritep);
		  dyy = cur_dosages2[sample_idx];
		  pzwritep = dtoa_gx(dyy, ' ', pzwritep);
		  dxx = 1.0 - dxx - dyy;
		  if (fabs(dxx) < SMALL_EPSILON) {
		    dxx = 0.0;
		  }
		  pzwritep = dtoa_gx(dxx, ' ', pzwritep);
		}
	      }
	      if (flex_pzwrite(&ps, &pzwritep)) {
		goto plink1_dosage_ret_WRITE_FAIL;
	      }
	    } while (ulii < sample_ct);
	  }
	  append_binary_eoln(&pzwritep);
	}
      }
      if (a1_ptr) {
	if (a1_ptr[1]) {
	  free(a1_ptr);
	}
	a1_ptr = nullptr;
	if (a2_ptr[1]) {
	  free(a2_ptr);
	}
	a2_ptr = nullptr;
      }
    }
  plink1_dosage_end_loop:
    for (file_idx = 0; file_idx < cur_batch_size; file_idx++) {
      if (gzclose(gz_infiles[file_idx]) != Z_OK) {
        gz_infiles[file_idx] = nullptr;
        goto plink1_dosage_ret_READ_FAIL;
      }
      gz_infiles[file_idx] = nullptr;
    }
    bigstack_reset(bigstack_mark);
  }
  if (do_score) {
    qrange_idx = 0;
    do {
      if (qrange_ct) {
	bufptr = strcpya(&(outname_end[1]), &(score_qrange_names[qrange_idx * max_qrange_name_len]));
	memcpy(bufptr, ".profile", 9);
      } else {
	memcpy(outname_end, ".profile", 9);
      }
      // this is not affected by 'gz' in PLINK 1.07; retain that for backward
      // compatibility.
      if (fopen_checked(outname, "w", &profile_outfile)) {
	goto plink1_dosage_ret_OPEN_FAIL;
      }
      sprintf(g_textbuf, "%%%us %%%us  PHENO%s %s\n", plink_maxfid, plink_maxiid, dosage_score_cnt? "    CNT" : "", score_report_average? "   SCORE" : "SCORESUM");
      fprintf(profile_outfile, g_textbuf, "FID", "IID");
      uii = score_range_obs_cts[qrange_idx];
      uiptr = &(score_miss_cts[sample_ct * qrange_idx]);
      dxx = score_bases[qrange_idx];
      dptr = &(cur_scores[sample_ct * qrange_idx]);
      for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_unsafe_ck(sample_exclude, &sample_uidx);
	bufptr = &(sample_ids[sample_uidx * max_sample_id_len]);
	bufptr2 = strchr(bufptr, '\t');
	bufptr = fw_strcpyn(plink_maxfid, (uintptr_t)(bufptr2 - bufptr), bufptr, g_textbuf);
	*bufptr++ = ' ';
	bufptr = fw_strcpy(plink_maxiid, &(bufptr2[1]), bufptr);
	*bufptr++ = ' ';
	if (IS_SET(pheno_nm, sample_uidx)) {
          if (pheno_c) {
	    bufptr = memseta(bufptr, 32, 5);
	    *bufptr++ = '1' + IS_SET(pheno_c, sample_uidx);
	  } else {
	    bufptr = width_force(6, bufptr, dtoa_g(pheno_d[sample_uidx], bufptr));
	  }
	} else {
	  bufptr = memcpya(bufptr, missing_pheno_str, missing_pheno_len);
	}
        *bufptr++ = ' ';
	ujj = uii - score_miss_cts[sample_idx];
	if (dosage_score_cnt) {
	  bufptr = uint32toa_w6x(uii, ' ', bufptr);
	}
        if (score_mean_impute) {
	  ujj = uii;
	}
	dyy = dxx + dptr[sample_idx];
	if (fabs(dyy) < SMALL_EPSILON) {
	  dyy = 0;
	} else if (score_report_average) {
	  dyy /= ((double)((int32_t)ujj));
	}
	bufptr = width_force(8, bufptr, dtoa_g(dyy, bufptr));
	*bufptr++ = '\n';
        if (fwrite_checked(g_textbuf, bufptr - g_textbuf, profile_outfile)) {
	  goto plink1_dosage_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&profile_outfile)) {
	goto plink1_dosage_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("--score: Results written to %s .\n", outname);
    } while (++qrange_idx < qrange_ct);
  } else {
    if (count_occur) {
      max_occur_id_len += sizeof(int32_t) + 1; // null, uint32_t
      if (bigstack_alloc_c(max_occur_id_len * distinct_id_ct, &bufptr)) {
	goto plink1_dosage_ret_NOMEM;
      }
      ulii = 0; // write idx
      ujj = 0; // number of counts > 1
      for (uii = 0; uii < HASHSIZE; uii++) {
	ll_ptr = htable[uii];
	while (ll_ptr) {
	  slen = strlen(ll_ptr->ss) + 1;
	  bufptr2 = memcpya(&(bufptr[ulii * max_occur_id_len]), ll_ptr->ss, slen);
	  ulii++;
	  memcpy(bufptr2, &(ll_ptr->ct), sizeof(int32_t));
	  if (ll_ptr->ct != 1) {
	    ujj++;
	  }
	  ll_ptr = ll_ptr->next;
	}
      }
      qsort(bufptr, distinct_id_ct, max_occur_id_len, strcmp_natural);
      for (ulii = 0; ulii < distinct_id_ct; ulii++) {
	bufptr2 = &(bufptr[ulii * max_occur_id_len]);
	slen = strlen(bufptr2);
	pzwritep = memcpyax(pzwritep, bufptr2, slen, ' ');
	pzwritep = uint32toa(*((uint32_t*)(&(bufptr2[slen + 1]))), pzwritep);
	append_binary_eoln(&pzwritep);
	if (flex_pzwrite(&ps, &pzwritep)) {
	  goto plink1_dosage_ret_WRITE_FAIL;
	}
      }
    }
    if (flex_pzwrite_close_null(&ps, pzwritep)) {
      goto plink1_dosage_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--%sdosage%s: Results saved to %s .\n", (do_glm || count_occur)? "" : "write-", count_occur? " occur" : "", outname);
    if (count_occur) {
      if (ujj) {
	LOGPRINTF("%u variant%s appeared multiple times.\n", ujj, (ujj == 1)? "" : "s");
      } else {
	logprint("No variant appeared multiple times.\n");
      }
    }
  }

  while (0) {
  plink1_dosage_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  plink1_dosage_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  plink1_dosage_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  plink1_dosage_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  plink1_dosage_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  plink1_dosage_ret_LONG_LINE:
    if (loadbuf_size == MAXLINEBUFLEN) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
      retval = RET_INVALID_FORMAT;
      break;
    }
    retval = RET_NOMEM;
    break;
  plink1_dosage_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
  plink1_dosage_ret_INVALID_FORMAT_WW:
    wordwrapb(0);
  plink1_dosage_ret_INVALID_FORMAT_2:
    logerrprintb();
  plink1_dosage_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  plink1_dosage_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  plink1_dosage_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
 plink1_dosage_ret_1:
  aligned_free_cond(pheno_c);
  free_cond(pheno_d);
  fclose_cond(phenofile);
  fclose_cond(infile);
  fclose_cond(profile_outfile);
  flex_pzwrite_close_cond(&ps, pzwritep);
  if (a1_ptr && a1_ptr[1]) {
    free(a1_ptr);
  }
  if (a2_ptr && a2_ptr[1]) {
    free(a2_ptr);
  }
  if (infile_ct) {
    for (uii = 0; uii < infile_ct; uii++) {
      gzclose_cond(gz_infiles[uii]);
    }
  }
  bigstack_end_reset(bigstack_end_mark);

  return retval;
}
