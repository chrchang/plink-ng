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
#include "plink2_misc.h"
#include "plink2_stats.h"

#ifdef __cplusplus
namespace plink2 {
#endif

pglerr_t update_sample_sexes(const char* update_sex_fname, const uintptr_t* sample_include, char* sample_ids, uint32_t raw_sample_ct, uintptr_t sample_ct, uintptr_t max_sample_id_blen, uint32_t update_sex_colm2, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    reterr = gzopen_read_checked(update_sex_fname, &gz_infile);
    if (reterr) {
      goto update_sample_sexes_ret_1;
    }
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    uintptr_t* already_seen;
    char* idbuf;
    if (bigstack_calloc_ul(raw_sample_ctl, &already_seen) ||
	bigstack_alloc_c(max_sample_id_blen, &idbuf)) {
      goto update_sample_sexes_ret_NOMEM;
    }
    uint32_t* id_map;
    char* sorted_idbox;
    if (copy_sort_strbox_subset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 1, 0, 0, &sorted_idbox, &id_map)) {
      goto update_sample_sexes_ret_NOMEM;
    }
    // permit very long lines since this can be pointed at .ped files
    uintptr_t loadbuf_size = bigstack_left();
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else {
      loadbuf_size &= ~(kCacheline - 1);
      if (loadbuf_size <= kMaxMediumLine) {
        goto update_sample_sexes_ret_NOMEM;
      }
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    uint32_t hit_ct = 0;
    uintptr_t miss_ct = 0;
    uintptr_t duplicate_ct = 0;
    while (gzgets(gz_infile, loadbuf, loadbuf_size)) {
      ++line_idx;
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-sex file is pathologically long.\n", line_idx);
	  goto update_sample_sexes_ret_MALFORMED_INPUT_2;
	}
      }
      char* fid_start = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*fid_start)) {
	continue;
      }
      char* fid_end = token_endnn(fid_start);
      char* iid_start = skip_initial_spaces(fid_end);
      if (is_eoln_kns(*iid_start)) {
	goto update_sample_sexes_ret_MISSING_TOKENS;
      }
      char* iid_end = token_endnn(iid_start);
      const uint32_t fid_slen = (uintptr_t)(fid_end - fid_start);
      const uint32_t iid_slen = (uintptr_t)(iid_end - iid_start);
      const uint32_t id_blen = fid_slen + iid_slen + 2;
      if (id_blen > max_sample_id_blen) {
	continue;
      }
      char* idbuf_iter = memcpyax(idbuf, fid_start, fid_slen, '\t');
      idbuf_iter = memcpya(idbuf_iter, iid_start, iid_slen);
      *idbuf_iter = '\0';
      uint32_t lb_idx = bsearch_str_lb(idbuf, sorted_idbox, id_blen, max_sample_id_blen, sample_ct);
      *idbuf_iter = ' ';
      const uint32_t ub_idx = bsearch_str_lb(idbuf, sorted_idbox, id_blen, max_sample_id_blen, sample_ct);
      if (ub_idx == lb_idx) {
	++miss_ct;
	continue;
      }
      char* sex_start = next_token_mult(iid_end, update_sex_colm2);
      if (!sex_start) {
	goto update_sample_sexes_ret_MISSING_TOKENS;
      }
      uint32_t sexval = (unsigned char)(*sex_start);
      const uint32_t ujj = sexval & 0xdfU;
      sexval -= 48;
      if (sexval > 2) {
	if (ujj == 77) {
	  // 'M'/'m'
	  sexval = 1;
	} else if (ujj == 70) {
	  // 'F'/'f'
	  sexval = 2;
	} else {
	  sprintf(g_logbuf, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file. (Acceptable values: 1/M/m = male, 2/F/f = female, 0 = missing.)\n", line_idx);
	  wordwrapb(0);
	  goto update_sample_sexes_ret_MALFORMED_INPUT_2;
	}
      }
      uint32_t sample_uidx = id_map[lb_idx];
      if (IS_SET(already_seen, sample_uidx)) {
	// permit duplicates iff sex value is identical
	const uint32_t old_sexval = IS_SET(sex_nm, sample_uidx) * (2 - IS_SET(sex_male, sample_uidx));
	if (sexval != old_sexval) {
	  idbuf[fid_slen] = ' ';
	  LOGPREPRINTFWW("Error: Duplicate sample ID '%s' with conflicting sex assignments in --update-sex file.\n", idbuf);
	  goto update_sample_sexes_ret_MALFORMED_INPUT_2;
	}
	++duplicate_ct;
	continue;
      }
      SET_BIT(sample_uidx, already_seen);
      while (1) {
	if (sexval) {
	  SET_BIT(sample_uidx, sex_nm);
	  if (sexval == 1) {
	    SET_BIT(sample_uidx, sex_male);
	  } else {
	    CLEAR_BIT(sample_uidx, sex_male);
	  }
	} else {
	  CLEAR_BIT(sample_uidx, sex_nm);
	  CLEAR_BIT(sample_uidx, sex_male);
	}
	++hit_ct;
	if (++lb_idx == ub_idx) {
	  break;
	}
	sample_uidx = id_map[lb_idx];
      }
    }
    if (!gzeof(gz_infile)) {
      goto update_sample_sexes_ret_READ_FAIL;
    }
    if (gzclose_null(&gz_infile)) {
      goto update_sample_sexes_ret_READ_FAIL;
    }
    if (duplicate_ct) {
      LOGPRINTFWW("Note: %" PRIuPTR " duplicate sample ID%s) in --update-sex file.\n", duplicate_ct, (duplicate_ct == 1)? " (with a consistent sex assignment" : "s (with consistent sex assignments");
    }
    if (miss_ct) {
      sprintf(g_logbuf, "--update-sex: %u sample%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      sprintf(g_logbuf, "--update-sex: %u sample%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logprintb();
  }
  while (0) {
  update_sample_sexes_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  update_sample_sexes_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  update_sample_sexes_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --update-sex file has fewer tokens than expected.\n", line_idx);
  update_sample_sexes_ret_MALFORMED_INPUT_2:
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  }
 update_sample_sexes_ret_1:
  gzclose_cond(gz_infile);
  bigstack_reset(bigstack_mark);
  return reterr;
}


pglerr_t process_boundary_token(char* tok_start, char* tok_end, const char* token_source_str, uint32_t max_boundary_ct, pglerr_t err_type, double* prev_boundary_ptr, uint32_t* boundary_ct_ptr, double** freq_bounds_ptr, uint64_t** dosage_bounds_ptr) {
  double cur_boundary;
  char* scan_end = scanadv_double(tok_start, &cur_boundary);
  if ((!scan_end) || (scan_end != tok_end)) {
    LOGERRPRINTF("Error: Invalid token in %s.\n", token_source_str);
    return err_type;
  }
  if (cur_boundary <= (*prev_boundary_ptr)) {
    logerrprint("Error: Invalid bin boundary sequence (must be strictly increasing, and start\nwith a positive number).\n");
    return err_type;
  }
  uint32_t boundary_ct = *boundary_ct_ptr;
  if (boundary_ct == max_boundary_ct) {
#ifdef __LP64__
    if (max_boundary_ct == 0x40000000) {
      logerrprint("Error: Too many bin boundaries.\n");
      return err_type;
    }
#endif
    return kPglRetNomem;
  }
  if (freq_bounds_ptr) {
    if (cur_boundary > 1.0) {
      logerrprint("Error: --freq bin boundary too large (must be <= 1).\n");
      return err_type;
    }
    // strictly-greater-than comparisons
    (*freq_bounds_ptr)[boundary_ct] = cur_boundary * (1 - kSmallEpsilon);
  } else {
    // max 2^31 - 3 variants
    if (cur_boundary > 4294967290.0) {
      logerrprint("Error: --freq counts bin boundary too large.\n");
      return err_type;
    }
    // due to the use of strictly-greater-than for comparison, we round
    // exact multiples of 1/32768 down
    const int64_t int_part = (int64_t)cur_boundary;
    const double cur_boundary_frac_part = cur_boundary - int_part;
    const int64_t int_part_scaled = int_part * kDosageMax;
    if (cur_boundary_frac_part == 0.0) {
      (*dosage_bounds_ptr)[boundary_ct] = int_part_scaled - 1;
    } else {
      (*dosage_bounds_ptr)[boundary_ct] = int_part_scaled + (int64_t)(cur_boundary_frac_part * (kDosageMax * (1 - kSmallEpsilon)));
    }
  }
  *prev_boundary_ptr = cur_boundary;
  *boundary_ct_ptr = boundary_ct + 1;
  return kPglRetSuccess;
}

pglerr_t init_histogram_from_file_or_commalist(const char* binstr, uint32_t is_fname, double** freq_bounds_ptr, uint64_t** dosage_bounds_ptr, uint32_t* boundary_ct_ptr, uint32_t** histogram_ptr) {
  gz_token_stream_t gts;
  gz_token_stream_preinit(&gts);
  uint32_t max_boundary_ct = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    uintptr_t ulii = bigstack_left() & (~(kCacheline - k1LU));
    if (ulii < 2 * kCacheline) {
      goto init_histogram_from_file_or_commalist_ret_NOMEM;
    }
    // 12 = 8 bytes for boundary value + 4 bytes for histogram entry
    ulii = (ulii - 2 * kCacheline) / 12;
#ifdef __LP64__
    max_boundary_ct = MINV(ulii, 0x40000000);
#else
    max_boundary_ct = ulii;
#endif
    if (freq_bounds_ptr) {
      *freq_bounds_ptr = (double*)g_bigstack_base;
    } else {
      *dosage_bounds_ptr = (uint64_t*)g_bigstack_base;
    }
    uint32_t boundary_ct = 0;
    double prev_boundary = 0.0;
    if (is_fname) {
      // we want to accept >100000 numbers on a single line.  this will reject
      // "000...{a million more zeroes}...1"; pretty sure that's okay.
      reterr = gz_token_stream_init(binstr, &gts, g_textbuf);
      if (reterr) {
	goto init_histogram_from_file_or_commalist_ret_1;
      }
      uint32_t token_slen;
      while (1) {
	char* token_start = gz_token_stream_advance(&gts, &token_slen);
	if (!token_start) {
	  break;
	}
	reterr = process_boundary_token(token_start, &(token_start[token_slen]), binstr, max_boundary_ct, kPglRetMalformedInput, &prev_boundary, &boundary_ct, freq_bounds_ptr, dosage_bounds_ptr);
	if (reterr) {
	  goto init_histogram_from_file_or_commalist_ret_1;
	}
      }
      if (token_slen) {
	if (token_slen == 0xffffffffU) {
	  sprintf(g_logbuf, "Error: Excessively long token in %s.\n", binstr);
	  goto init_histogram_from_file_or_commalist_ret_MALFORMED_INPUT_2;
	} else {
	  goto init_histogram_from_file_or_commalist_ret_READ_FAIL;
	}
      }
    } else {
      // const_cast
      char* binstr_iter = (char*)((uintptr_t)binstr);
      while (1) {
	char* tok_end = strchr(binstr_iter, ',');
	const uint32_t is_last_token = (tok_end == nullptr);
	if (is_last_token) {
	  tok_end = (char*)rawmemchr(binstr_iter, '\0');
	}
	reterr = process_boundary_token(binstr_iter, tok_end, "--freq {ref,alt1}bins= list", max_boundary_ct, kPglRetInvalidCmdline, &prev_boundary, &boundary_ct, freq_bounds_ptr, dosage_bounds_ptr);
	if (reterr) {
	  goto init_histogram_from_file_or_commalist_ret_1;
	}
	if (is_last_token) {
	  break;
	}
	binstr_iter = &(tok_end[1]);
      }
    }
    *boundary_ct_ptr = boundary_ct;
    g_bigstack_base += round_up_pow2(boundary_ct * (8 * k1LU), kCacheline);
    *histogram_ptr = (uint32_t*)bigstack_alloc_raw(round_up_pow2((boundary_ct + 1) * sizeof(int32_t), kCacheline));
    fill_uint_zero(boundary_ct + 1, *histogram_ptr);
  }
  while (0) {
  init_histogram_from_file_or_commalist_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  init_histogram_from_file_or_commalist_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  init_histogram_from_file_or_commalist_ret_MALFORMED_INPUT_2:
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  }
 init_histogram_from_file_or_commalist_ret_1:
  gz_token_stream_close(&gts);
  return reterr;
}

pglerr_t write_allele_freqs(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bp, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* founder_allele_dosages, const char* ref_binstr, const char* alt1_binstr, uint32_t variant_ct, uint32_t max_alt_allele_ct, uint32_t max_allele_slen, allele_freq_t allele_freq_modifier, uint32_t nonfounders, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t counts = (allele_freq_modifier / kfAlleleFreqCounts) & 1;
    if (counts) {
      strcpy(outname_end, ".acount");
    } else {
      strcpy(outname_end, ".afreq");
    }
    if (!(allele_freq_modifier & kfAlleleFreqBinsOnly)) {
      const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
      unsigned char* overflow_buf;
      if (bigstack_alloc_uc(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + max_alt_allele_ct * (24 * k1LU) + 2 * max_allele_slen, &overflow_buf)) {
	goto write_allele_freqs_ret_NOMEM;
      }
      const uint32_t output_zst = allele_freq_modifier & kfAlleleFreqZs;
      if (output_zst) {
	strcpy(&(outname_end[6 + counts]), ".zst");
      }
      if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
	goto write_allele_freqs_ret_OPEN_FAIL;
      }
      cswritep = (char*)overflow_buf;
      *cswritep++ = '#';
      const uint32_t chr_col = allele_freq_modifier & kfAlleleFreqColChrom;

      // includes trailing tab
      char* chr_buf = nullptr;
      if (chr_col) {
	if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
	  goto write_allele_freqs_ret_NOMEM;
	}
	cswritep = strcpya(cswritep, "CHROM\t");
      }
      if (allele_freq_modifier & kfAlleleFreqColPos) {
	cswritep = strcpya(cswritep, "POS\t");
      } else {
	variant_bp = nullptr;
      }
      cswritep = strcpya(cswritep, "ID");
      const uint32_t ref_col = allele_freq_modifier & kfAlleleFreqColRef;
      if (ref_col) {
	cswritep = strcpya(cswritep, "\tREF");
      }
      const uint32_t alt1_col = allele_freq_modifier & kfAlleleFreqColAlt1;
      if (alt1_col) {
	cswritep = strcpya(cswritep, "\tALT1");
      }
      const uint32_t alt_col = allele_freq_modifier & kfAlleleFreqColAlt;
      if (alt_col) {
	cswritep = strcpya(cswritep, "\tALT");
      }
      const uint32_t reffreq_col = allele_freq_modifier & kfAlleleFreqColReffreq;
      if (reffreq_col) {
	cswritep = strcpya(cswritep, "\tREF_");
	if (counts) {
	  cswritep = strcpya(cswritep, "CT");
	} else {
	  cswritep = strcpya(cswritep, "FREQ");
	}
      }
      const uint32_t alt1freq_col = allele_freq_modifier & kfAlleleFreqColAlt1freq;
      if (alt1freq_col) {
	cswritep = strcpya(cswritep, "\tALT1_");
	if (counts) {
	  cswritep = strcpya(cswritep, "CT");
	} else {
	  cswritep = strcpya(cswritep, "FREQ");
	}
      }
      const uint32_t freq_col = allele_freq_modifier & (kfAlleleFreqColFreq | kfAlleleFreqColAltfreq);
      const uint32_t commalist_exclude_ref = (allele_freq_modifier & (kfAlleleFreqColAltfreq | kfAlleleFreqColAlteq | kfAlleleFreqColAlteqz | kfAlleleFreqColAltnumeq))? 1 : 0;
      const uint32_t eq_col = allele_freq_modifier & (kfAlleleFreqColEq | kfAlleleFreqColEqz | kfAlleleFreqColAlteq | kfAlleleFreqColAlteqz | kfAlleleFreqColNumeq | kfAlleleFreqColAltnumeq);
      const uint32_t eq_includez = allele_freq_modifier & (kfAlleleFreqColEqz | kfAlleleFreqColAlteqz);
      const uint32_t eq_num = allele_freq_modifier & (kfAlleleFreqColNumeq | kfAlleleFreqColAltnumeq);
      if (freq_col || eq_col) {
	*cswritep++ = '\t';
	if (commalist_exclude_ref) {
	  cswritep = strcpya(cswritep, "ALT_");
	}
	if (counts) {
	  cswritep = memcpyl3a(cswritep, "CTS");
	} else {
	  cswritep = strcpya(cswritep, "FREQS");
	}
      }
      const uint32_t nobs_col = allele_freq_modifier & kfAlleleFreqColNobs;
      if (nobs_col) {
	cswritep = strcpya(cswritep, "\tOBS_CT");
      }
      append_binary_eoln(&cswritep);

      uint32_t variant_uidx = 0;
      uint32_t chr_fo_idx = 0xffffffffU;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t cur_allele_ct = 2;
      printf("--freq%s%s: 0%%", output_zst? " zs" : "", counts? " counts" : "");
      fflush(stdout);
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	if (chr_col) {
	  if (variant_uidx >= chr_end) {
	    do {
	      ++chr_fo_idx;
	      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	    } while (variant_uidx >= chr_end);
	    char* chr_name_end = chr_name_write(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
	    *chr_name_end = '\t';
	    chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	  }
	  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
	}
	if (variant_bp) {
	  cswritep = uint32toa_x(variant_bp[variant_uidx], '\t', cswritep);
	}
	cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
	uintptr_t variant_allele_idx_base = variant_uidx * 2;
	if (variant_allele_idxs) {
	  variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	  cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
	}
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
	      goto write_allele_freqs_ret_WRITE_FAIL;
	    }
	    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	  }
	  --cswritep;
	}
	const uint64_t* cur_allele_dosages = &(founder_allele_dosages[variant_allele_idx_base]);
	uint64_t tot_allele_dosage = cur_allele_dosages[0];
	for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
	  tot_allele_dosage += cur_allele_dosages[allele_idx];
	}
	double tot_allele_dosage_recip = 0.0;
	if (!counts) {
	  tot_allele_dosage_recip = 1.0 / ((double)((int64_t)tot_allele_dosage));
	}
	if (reffreq_col) {
	  *cswritep++ = '\t';
	  const uint64_t ref_dosage = cur_allele_dosages[0];
	  if (counts) {
	    cswritep = print_dosage(ref_dosage, cswritep);
	  } else {
	    cswritep = dtoa_g(((double)((int64_t)ref_dosage)) * tot_allele_dosage_recip, cswritep);
	  }
	}
	if (alt1freq_col) {
	  *cswritep++ = '\t';
	  const uint64_t alt1_dosage = cur_allele_dosages[1];
	  if (counts) {
	    cswritep = print_dosage(alt1_dosage, cswritep);
	  } else {
	    cswritep = dtoa_g(((double)((int64_t)alt1_dosage)) * tot_allele_dosage_recip, cswritep);
	  }
	}
	if (freq_col) {
	  *cswritep++ = '\t';
	  for (uint32_t allele_idx = commalist_exclude_ref; allele_idx < cur_allele_ct; ++allele_idx) {
	    const uint64_t cur_allele_dosage = cur_allele_dosages[allele_idx];
	    if (counts) {
	      cswritep = print_dosage(cur_allele_dosage, cswritep);
	    } else {
	      cswritep = dtoa_g(((double)((int64_t)cur_allele_dosage)) * tot_allele_dosage_recip, cswritep);
	    }
	    *cswritep++ = ',';
	  }
	  --cswritep;
	} else if (eq_col) {
	  *cswritep++ = '\t';
	  uint32_t at_least_one_entry = 0;
	  for (uint32_t allele_idx = commalist_exclude_ref; allele_idx < cur_allele_ct; ++allele_idx) {
	    const uint64_t cur_allele_dosage = cur_allele_dosages[allele_idx];
	    if (eq_includez || cur_allele_dosage) {
	      if (eq_num) {
		cswritep = uint32toa(allele_idx, cswritep);
	      } else {
		if (cswrite(&css, &cswritep)) {
		  goto write_allele_freqs_ret_WRITE_FAIL;
		}
		// possible todo: warning if an allele code includes a '='
		cswritep = strcpya(cswritep, cur_alleles[allele_idx]);
	      }
	      *cswritep++ = '=';
	      if (counts) {
		cswritep = print_dosage(cur_allele_dosage, cswritep);
	      } else {
		cswritep = dtoa_g(((double)((int64_t)cur_allele_dosage)) * tot_allele_dosage_recip, cswritep);
	      }
	      *cswritep++ = ',';
	      at_least_one_entry = 1;
	    }
	  }
	  if (at_least_one_entry) {
	    --cswritep;
	  } else {
	    *cswritep++ = '.';
	  }
	}
	if (nobs_col) {
	  *cswritep++ = '\t';
	  cswritep = print_dosage(tot_allele_dosage, cswritep);
	}
	append_binary_eoln(&cswritep);
	if (cswrite(&css, &cswritep)) {
	  goto write_allele_freqs_ret_WRITE_FAIL;
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
      }
      if (cswrite_close_null(&css, cswritep)) {
	goto write_allele_freqs_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("--freq%s%s: Allele %s (%s) written to %s .\n", output_zst? " zs" : "", counts? " counts" : "", counts? "counts" : "frequencies", nonfounders? "all samples" : "founders only", outname);
    }

    if (ref_binstr || alt1_binstr) {
      bigstack_reset(bigstack_mark);
      double* ref_freq_bounds = nullptr;
      uint64_t* ref_dosage_bounds = nullptr;
      uint32_t* ref_histogram = nullptr;
      uint32_t ref_boundary_ct = 0;
      if (ref_binstr) {
	reterr = init_histogram_from_file_or_commalist(ref_binstr, (allele_freq_modifier / kfAlleleFreqBinsRefFname) & 1, counts? nullptr : (&ref_freq_bounds), counts? (&ref_dosage_bounds) : nullptr, &ref_boundary_ct, &ref_histogram);
	if (reterr) {
	  goto write_allele_freqs_ret_1;
	}
      }
      double* alt1_freq_bounds = nullptr;
      uint64_t* alt1_dosage_bounds = nullptr;
      uint32_t* alt1_histogram = nullptr;
      uint32_t alt1_boundary_ct = 0;
      if (alt1_binstr) {
	reterr = init_histogram_from_file_or_commalist(alt1_binstr, (allele_freq_modifier / kfAlleleFreqBinsAlt1Fname) & 1, counts? nullptr : (&alt1_freq_bounds), counts? (&alt1_dosage_bounds) : nullptr, &alt1_boundary_ct, &alt1_histogram);
	if (reterr) {
	  goto write_allele_freqs_ret_1;
	}
      }
      uint32_t variant_uidx = 0;
      if (!counts) {
        uint32_t cur_allele_ct = 2;
	for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	  next_set_unsafe_ck(variant_include, &variant_uidx);
	  uintptr_t variant_allele_idx_base = variant_uidx * 2;
	  if (variant_allele_idxs) {
	    variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	    cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
	  }
	  const uint64_t* cur_allele_dosages = &(founder_allele_dosages[variant_allele_idx_base]);
	  const uint64_t ref_allele_dosage = cur_allele_dosages[0];
	  const uint64_t alt1_allele_dosage = cur_allele_dosages[1];
	  uint64_t tot_allele_dosage = ref_allele_dosage + alt1_allele_dosage;
	  for (uint32_t allele_idx = 2; allele_idx < cur_allele_ct; ++allele_idx) {
	    tot_allele_dosage += cur_allele_dosages[allele_idx];
	  }
	  const double tot_allele_dosage_recip = 1.0 / ((double)((int64_t)tot_allele_dosage));
	  if (ref_histogram) {
	    ref_histogram[doublearr_greater_than(ref_freq_bounds, ref_boundary_ct, ref_allele_dosage * tot_allele_dosage_recip)] += 1;
	  }
	  if (alt1_histogram) {
	    alt1_histogram[doublearr_greater_than(alt1_freq_bounds, alt1_boundary_ct, alt1_allele_dosage * tot_allele_dosage_recip)] += 1;
	  }
	}
      } else {
	for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	  next_set_unsafe_ck(variant_include, &variant_uidx);
	  uintptr_t variant_allele_idx_base = variant_uidx * 2;
	  if (variant_allele_idxs) {
	    variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	  }
	  const uint64_t* cur_allele_dosages = &(founder_allele_dosages[variant_allele_idx_base]);
	  if (ref_histogram) {
	    ref_histogram[uint64arr_greater_than(ref_dosage_bounds, ref_boundary_ct, cur_allele_dosages[0])] += 1;
	  }
	  if (alt1_histogram) {
	    alt1_histogram[uint64arr_greater_than(alt1_dosage_bounds, alt1_boundary_ct, cur_allele_dosages[1])] += 1;
	  }
	}
      }
      for (uint32_t is_alt1 = 0; is_alt1 < 2; ++is_alt1) {
	const uint32_t* cur_histogram = is_alt1? alt1_histogram : ref_histogram;
	if (!cur_histogram) {
	  continue;
	}
	char* outname_end2 = &(outname_end[6 + counts]);
	if (!is_alt1) {
	  outname_end2 = strcpya(outname_end2, ".ref");
	} else {
	  outname_end2 = strcpya(outname_end2, ".alt1");
	}
	strcpy(outname_end2, ".bins");
	if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	  goto write_allele_freqs_ret_OPEN_FAIL;
	}
	char* textbuf = g_textbuf;
	char* textbuf_flush = &(textbuf[kMaxMediumLine]);
	char* write_iter = strcpya(textbuf, "#BIN_START\tOBS_CT" EOLN_STR);
	const uint32_t cur_boundary_ct = is_alt1? alt1_boundary_ct : ref_boundary_ct;
	if (!counts) {
	  const double* cur_freq_bounds = is_alt1? alt1_freq_bounds : ref_freq_bounds;
	  for (uint32_t bin_idx = 0; bin_idx <= cur_boundary_ct; ++bin_idx) {
	    if (!bin_idx) {
	      *write_iter++ = '0';
	    } else {
	      write_iter = dtoa_g(cur_freq_bounds[bin_idx - 1] * (1.0 / (1 - kSmallEpsilon)), write_iter);
	    }
	    *write_iter++ = '\t';
	    write_iter = uint32toa(cur_histogram[bin_idx], write_iter);
	    append_binary_eoln(&write_iter);
	    if (write_iter >= textbuf_flush) {
	      if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), outfile)) {
		goto write_allele_freqs_ret_WRITE_FAIL;
	      }
	      write_iter = textbuf;
	    }
	  }
	} else {
	  const uint64_t* cur_dosage_bounds = is_alt1? alt1_dosage_bounds : ref_dosage_bounds;
	  for (uint32_t bin_idx = 0; bin_idx <= cur_boundary_ct; ++bin_idx) {
	    if (!bin_idx) {
	      *write_iter++ = '0';
	    } else {
	      write_iter = print_dosage(1 + cur_dosage_bounds[bin_idx - 1], write_iter);
	    }
	    *write_iter++ = '\t';
	    write_iter = uint32toa(cur_histogram[bin_idx], write_iter);
	    append_binary_eoln(&write_iter);
	    if (write_iter >= textbuf_flush) {
	      if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), outfile)) {
		goto write_allele_freqs_ret_WRITE_FAIL;
	      }
	      write_iter = textbuf;
	    }
	  }
	}
	if (write_iter != textbuf) {
	  if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), outfile)) {
	    goto write_allele_freqs_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto write_allele_freqs_ret_WRITE_FAIL;
	}
	const uint32_t cur_is_file = allele_freq_modifier & (is_alt1? kfAlleleFreqBinsAlt1Fname : kfAlleleFreqBinsRefFname);
        LOGPRINTFWW("--freq%s %sbins%s=: Histogram written to %s .\n", counts? " counts" : "", is_alt1? "alt1" : "ref", cur_is_file? "-file" : "", outname);
      }
    }
  }
  while (0) {
  write_allele_freqs_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_allele_freqs_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_allele_freqs_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 write_allele_freqs_ret_1:
  cswrite_close_cond(&css, cswritep);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t write_geno_counts(__attribute__((unused)) const uintptr_t* sample_include, __attribute__((unused)) const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bp, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint32_t* raw_geno_cts, const uint32_t* x_male_geno_cts, __attribute__((unused)) uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t variant_ct, uint32_t x_start, uint32_t max_allele_slen, geno_counts_t geno_counts_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    unsigned char* overflow_buf;
    char* chr_buf;
    if (bigstack_alloc_uc(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + simple_pgrp->fi.max_alt_allele_ct * (24 * k1LU) + 2 * max_allele_slen, &overflow_buf) ||
	bigstack_alloc_c(max_chr_blen, &chr_buf)) {
      goto write_geno_counts_ret_NOMEM;
    }
    /*
    // need the following once multiallelic variants are supported
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t raw_sample_ctv = BITCT_TO_VECCT(raw_sample_ct);
    uintptr_t* sample_include_interleaved_vec = nullptr;
    uint32_t* sample_include_cumulative_popcounts = nullptr;
    uintptr_t* genovec = nullptr;
    uintptr_t* sex_male_interleaved_vec = nullptr;
    uint32_t* sex_male_cumulative_popcounts = nullptr;
    if (simple_pgrp->fi.max_alt_allele_ct > 1) {
      if (bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &sample_include_interleaved_vec) ||
	  bigstack_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
	  bigstack_alloc_ul(QUATERCT_TO_WORDCT(raw_sample_ct), &genovec) ||
	  bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &sex_male_interleaved_vec) ||
	  bigstack_alloc_ui(raw_sample_ctl, &sex_male_cumulative_popcounts)) {
	goto write_geno_counts_ret_NOMEM;
      }
      fill_interleaved_mask_vec(sample_include, raw_sample_ctv, sample_include_interleaved_vec);
      fill_cumulative_popcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      fill_interleaved_mask_vec(sex_male, raw_sample_ctv, sex_male_interleaved_vec);
      fill_cumulative_popcounts(sex_male, raw_sample_ctl, sex_male_cumulative_popcounts);
      pgr_clear_ld_cache(simple_pgrp);
    }
    */
    const uint32_t output_zst = geno_counts_modifier & kfGenoCountsZs;
    char* outname_end2 = strcpya0(outname_end, ".gcount");
    if (output_zst) {
      strcpy(outname_end2, ".zst");
    }
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto write_geno_counts_ret_OPEN_FAIL;
    }
    cswritep = (char*)overflow_buf;
    *cswritep++ = '#';
    const uint32_t chr_col = geno_counts_modifier & kfGenoCountsColChrom;

    // includes trailing tab
    if (chr_col) {
      cswritep = strcpya(cswritep, "CHROM\t");
    }
    if (geno_counts_modifier & kfGenoCountsColPos) {
      cswritep = strcpya(cswritep, "POS\t");
    } else {
      variant_bp = nullptr;
    }
    cswritep = strcpya(cswritep, "ID");
    const uint32_t ref_col = geno_counts_modifier & kfGenoCountsColRef;
    if (ref_col) {
      cswritep = strcpya(cswritep, "\tREF");
    }
    const uint32_t alt1_col = geno_counts_modifier & kfGenoCountsColAlt1;
    if (alt1_col) {
      cswritep = strcpya(cswritep, "\tALT1");
    }
    const uint32_t alt_col = geno_counts_modifier & kfGenoCountsColAlt;
    if (alt_col) {
      cswritep = strcpya(cswritep, "\tALT");
    }
    const uint32_t homref_col = geno_counts_modifier & kfGenoCountsColHomref;
    if (homref_col) {
      cswritep = strcpya(cswritep, "\tHOM_REF_CT");
    }
    const uint32_t refalt1_col = geno_counts_modifier & kfGenoCountsColRefalt1;
    if (refalt1_col) {
      cswritep = strcpya(cswritep, "\tHET_REF_ALT1_CT");
    }
    const uint32_t refalt_col = geno_counts_modifier & kfGenoCountsColRefalt;
    if (refalt_col) {
      cswritep = strcpya(cswritep, "\tHET_REF_ALT_CTS");
    }
    const uint32_t homalt1_col = geno_counts_modifier & kfGenoCountsColHomalt1;
    if (homalt1_col) {
      cswritep = strcpya(cswritep, "\tHOM_ALT1_CT");
    }
    const uint32_t xy_col = geno_counts_modifier & (kfGenoCountsColAltxy | kfGenoCountsColXy);
    const uint32_t xy_col_altonly = (geno_counts_modifier / kfGenoCountsColAltxy) & 1;
    if (xy_col) {
      *cswritep++ = '\t';
      if (xy_col_altonly) {
	cswritep = strcpya(cswritep, "NONREF_");
      }
      cswritep = strcpya(cswritep, "DIPLOID_GENO_CTS");
    }
    const uint32_t hapref_col = geno_counts_modifier & kfGenoCountsColHapref;
    if (hapref_col) {
      cswritep = strcpya(cswritep, "\tHAP_REF_CT");
    }
    const uint32_t hapalt1_col = geno_counts_modifier & kfGenoCountsColHapalt1;
    if (hapalt1_col) {
      cswritep = strcpya(cswritep, "\tHAP_ALT1_CT");
    }
    const uint32_t hap_col = geno_counts_modifier & (kfGenoCountsColHapalt | kfGenoCountsColHap);
    const uint32_t hap_col_altonly = (geno_counts_modifier / kfGenoCountsColHapalt) & 1;
    if (hap_col) {
      if (hap_col_altonly) {
	cswritep = strcpya(cswritep, "\tHAP_ALT_CTS");
      } else {
	cswritep = strcpya(cswritep, "\tHAP_CTS");
      }
    }
    const uint32_t numeq_col = geno_counts_modifier & kfGenoCountsColNumeq;
    if (numeq_col) {
      cswritep = strcpya(cswritep, "\tGENO_CTS");
    }
    const uint32_t missing_col = geno_counts_modifier & kfGenoCountsColMissing;
    if (missing_col) {
      cswritep = strcpya(cswritep, "\tMISSING_CT");
    }
    const uint32_t nobs_col = geno_counts_modifier & kfGenoCountsColNobs;
    if (nobs_col) {
      cswritep = strcpya(cswritep, "\tOBS_CT");
    }
    append_binary_eoln(&cswritep);

    const int32_t x_code = cip->xymt_codes[kChrOffsetX];
    const int32_t y_code = cip->xymt_codes[kChrOffsetY];
    uint32_t is_autosomal_diploid = 0; // also includes MT for now
    uint32_t is_x = 0;
    uint32_t nobs_base = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t variant_uidx = 0;
    uint32_t homref_ct = 0;
    uint32_t het_ct = 0;
    uint32_t homalt1_ct = 0;
    uint32_t hapref_ct = 0;
    uint32_t hapalt1_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    printf("--geno-counts%s: 0%%", output_zst? " zs" : "");
    fflush(stdout);
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	*chr_name_end = '\t';
	chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	is_autosomal_diploid = !is_set(cip->haploid_mask, chr_idx);
	nobs_base = sample_ct;
	is_x = (chr_idx == x_code);
	/*
        if (!is_autosomal_diploid) {
          pgr_clear_ld_cache(simple_pgrp);
          // move chr_idx == y_code check here
          // update cur_sample_include, etc.
        }
	 */
	if (chr_idx == y_code) {
	  nobs_base = male_ct;
	}
	homref_ct = 0;
	het_ct = 0;
	homalt1_ct = 0;
	hapref_ct = 0;
	hapalt1_ct = 0;
      }
      if (chr_col) {
	cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
      }
      if (variant_bp) {
	cswritep = uint32toa_x(variant_bp[variant_uidx], '\t', cswritep);
      }
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
	    goto write_geno_counts_ret_WRITE_FAIL;
	  }
	  cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	}
	--cswritep;
      }
      uint32_t missing_ct;
      if (cur_allele_ct == 2) {
	const uint32_t* cur_raw_geno_cts = &(raw_geno_cts[(3 * k1LU) * variant_uidx]);
	if (is_autosomal_diploid) {
	  homref_ct = cur_raw_geno_cts[0];
	  het_ct = cur_raw_geno_cts[1];
	  homalt1_ct = cur_raw_geno_cts[2];
	  missing_ct = nobs_base - homref_ct - het_ct - homalt1_ct;
	} else {
	  homref_ct = cur_raw_geno_cts[0];
	  het_ct = cur_raw_geno_cts[1];
	  homalt1_ct = cur_raw_geno_cts[2];
	  if (is_x) {
	    if (x_male_geno_cts) {
	      const uint32_t* cur_male_geno_cts = &(x_male_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	      hapref_ct = cur_male_geno_cts[0];
	      homref_ct -= hapref_ct;
	      het_ct -= cur_male_geno_cts[1];
	      hapalt1_ct = cur_male_geno_cts[2];
	      homalt1_ct -= hapalt1_ct;
	    }
	    missing_ct = nobs_base - homref_ct - het_ct - homalt1_ct - hapref_ct - hapalt1_ct;
	  } else {
	    // chrY or haploid
	    hapref_ct = cur_raw_geno_cts[0];
	    hapalt1_ct = cur_raw_geno_cts[2];
	    missing_ct = nobs_base - hapref_ct - hapalt1_ct;
	  }
	}
	if (homref_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(homref_ct, cswritep);
	}
	if (refalt1_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(het_ct, cswritep);
	}
	if (refalt_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(het_ct, cswritep);
	}
	if (homalt1_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(homalt1_ct, cswritep);
	}
	if (xy_col_altonly) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(homalt1_ct, cswritep);
	} else if (xy_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa_x(homref_ct, ',', cswritep);
	  cswritep = uint32toa_x(het_ct, ',', cswritep);
	  cswritep = uint32toa(homalt1_ct, cswritep);
	}
	if (hapref_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(hapref_ct, cswritep);
	}
	if (hapalt1_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(hapalt1_ct, cswritep);
	}
	if (hap_col) {
	  *cswritep++ = '\t';
	  if (!hap_col_altonly) {
	    cswritep = uint32toa_x(hapref_ct, ',', cswritep);
	  }
	  cswritep = uint32toa(hapalt1_ct, cswritep);
	}
	if (numeq_col) {
	  *cswritep++ = '\t';
	  if (homref_ct) {
	    cswritep = strcpya(cswritep, "0/0=");
	    cswritep = uint32toa_x(homref_ct, ',', cswritep);
	  }
	  if (het_ct) {
	    cswritep = strcpya(cswritep, "0/1=");
	    cswritep = uint32toa_x(het_ct, ',', cswritep);
	  }
	  if (homalt1_ct) {
	    cswritep = strcpya(cswritep, "1/1=");
	    cswritep = uint32toa_x(homalt1_ct, ',', cswritep);
	  }
	  if (hapref_ct) {
	    cswritep = strcpya(cswritep, "0=");
	    cswritep = uint32toa_x(hapref_ct, ',', cswritep);
	  }
	  if (hapalt1_ct) {
	    cswritep = strcpya(cswritep, "1=");
	    cswritep = uint32toa_x(hapalt1_ct, ',', cswritep);
	  }
	  if (missing_ct != nobs_base) {
	    --cswritep;
	  } else {
	    *cswritep++ = '.';
	  }
	}
      } else {
	// todo
	missing_ct = 0;
	assert(0);
      }
      if (missing_col) {
	*cswritep++ = '\t';
	cswritep = uint32toa(missing_ct, cswritep);
      }
      if (nobs_col) {
	*cswritep++ = '\t';
	cswritep = uint32toa(nobs_base - missing_ct, cswritep);
      }
      append_binary_eoln(&cswritep);
      if (cswrite(&css, &cswritep)) {
	goto write_geno_counts_ret_WRITE_FAIL;
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
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto write_geno_counts_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    LOGPRINTFWW("--geno-counts%s: Genotype counts written to %s .\n", output_zst? " zs" : "", outname);
  }
  while (0) {
  write_geno_counts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_geno_counts_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_geno_counts_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  // write_geno_counts_ret_1:
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t write_missingness_reports(const uintptr_t* sample_include, const uintptr_t* sex_male, const char* sample_ids, const char* sids, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* sample_missing_geno_cts, const uint32_t* sample_hethap_cts, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bp, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint32_t* variant_missing_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uintptr_t max_allele_slen, uint32_t first_hap_uidx, missing_rpt_t missing_rpt_modifier, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t output_zst = missing_rpt_modifier & kfMissingRptZs;
    if (!(missing_rpt_modifier & kfMissingRptVariantOnly)) {
      unsigned char* overflow_buf;
      if (bigstack_alloc_uc(kCompressStreamBlock + kMaxMediumLine + pheno_ct * 2, &overflow_buf)) {
	goto write_missingness_reports_ret_NOMEM;
      }
      char* outname_end2 = strcpya0(outname_end, ".smiss");
      if (output_zst) {
	strcpy(outname_end2, ".zst");
      }
      if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
	goto write_missingness_reports_ret_OPEN_FAIL;
      }
      cswritep = strcpya((char*)overflow_buf, "#FID\tIID");
      const uint32_t scol_sid = sid_col_required(sample_include, sids, sample_ct, max_sid_blen, missing_rpt_modifier / kfMissingRptScolMaybesid);
      if (scol_sid) {
	cswritep = strcpya(cswritep, "\tSID");
      }
      const uint32_t scol_empty_pheno = (missing_rpt_modifier & kfMissingRptScolMisspheno1) && (!pheno_ct);
      if (scol_empty_pheno) {
	cswritep = strcpya(cswritep, "\tMISS_PHENO1");
      }
      const uint32_t scol_phenos = (missing_rpt_modifier & (kfMissingRptScolMisspheno1 | kfMissingRptScolMissphenos)) && pheno_ct;
      if (scol_phenos) {
	if (!(missing_rpt_modifier & kfMissingRptScolMissphenos)) {
	  pheno_ct = 1;
	}
	for (uintptr_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  *cswritep++ = '\t';
	  cswritep = strcpya(cswritep, &(pheno_names[pheno_idx * max_pheno_name_blen]));
	  if (cswrite(&css, &cswritep)) {
	    goto write_missingness_reports_ret_WRITE_FAIL;
	  }
	}
      }
      const uint32_t scol_nmiss = (missing_rpt_modifier / kfMissingRptScolNmiss) & 1;
      if (scol_nmiss) {
	cswritep = strcpya(cswritep, "\tMISSING_CT");
      }
      const uint32_t scol_nmiss_hh = (missing_rpt_modifier / kfMissingRptScolNmissHh) & 1;
      if (scol_nmiss_hh) {
	cswritep = strcpya(cswritep, "\tMISSING_AND_HETHAP_CT");
      }
      const uint32_t scol_hethap = (missing_rpt_modifier / kfMissingRptScolHethap) & 1;
      if (scol_hethap) {
	cswritep = strcpya(cswritep, "\tHETHAP_CT");
      }
      const uint32_t scol_nobs = (missing_rpt_modifier / kfMissingRptScolNobs) & 1;
      if (scol_nobs) {
	cswritep = strcpya(cswritep, "\tOBS_CT");
      }
      const uint32_t scol_fmiss = (missing_rpt_modifier / kfMissingRptScolFmiss) & 1;
      if (scol_fmiss) {
	cswritep = strcpya(cswritep, "\tF_MISS");
      }
      const uint32_t scol_fmiss_hh = (missing_rpt_modifier / kfMissingRptScolFmissHh) & 1;
      if (scol_fmiss_hh) {
	cswritep = strcpya(cswritep, "\tF_MISS_AND_HETHAP");
      }
      append_binary_eoln(&cswritep);
      uint32_t variant_ct_y = 0;
      int32_t y_code;
      if (xymt_exists(cip, kChrOffsetY, &y_code)) {
	variant_ct_y = count_chr_variants_unsafe(variant_include, cip, y_code);
      }
      const uint32_t variant_ct_nony = variant_ct - variant_ct_y;
      char nobs_strs[2][16];
      uint32_t nobs_slens[2];
      double variant_ct_recips[2];
      char* write_iter = nobs_strs[0];
      *write_iter++ = '\t';
      write_iter = uint32toa(variant_ct_nony, write_iter);
      nobs_slens[0] = (uintptr_t)(write_iter - nobs_strs[0]);
      variant_ct_recips[0] = 1.0 / ((double)((int32_t)variant_ct_nony));
      write_iter = nobs_strs[1];
      *write_iter++ = '\t';
      write_iter = uint32toa(variant_ct, write_iter);
      nobs_slens[1] = (uintptr_t)(write_iter - nobs_strs[1]);
      variant_ct_recips[1] = 1.0 / ((double)((int32_t)variant_ct));
      uintptr_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
        next_set_ul_unsafe_ck(sample_include, &sample_uidx);
	cswritep = strcpya(cswritep, &(sample_ids[sample_uidx * max_sample_id_blen]));
	if (scol_sid) {
	  *cswritep++ = '\t';
	  if (sids) {
	    cswritep = strcpya(cswritep, &(sids[sample_uidx * max_sid_blen]));
	  } else {
	    *cswritep++ = '0';
	  }
	}
	if (scol_phenos) {
	  for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	    *cswritep++ = '\t';
	    // 'Y' - 'N' == 11
	    *cswritep++ = 'Y' - 11 * IS_SET(pheno_cols[pheno_idx].nonmiss, sample_uidx);
	    if (cswrite(&css, &cswritep)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  if (scol_empty_pheno) {
	    cswritep = strcpya(cswritep, "\tY");
	  }
	  if (cswrite(&css, &cswritep)) {
	    goto write_missingness_reports_ret_WRITE_FAIL;
	  }
	}
	const uint32_t cur_missing_geno_base = sample_missing_geno_cts[sample_uidx];
	if (scol_nmiss) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(cur_missing_geno_base, cswritep);
	}
	if (scol_nmiss_hh) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(cur_missing_geno_base + sample_hethap_cts[sample_uidx], cswritep);
	}
	if (scol_hethap) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(sample_hethap_cts[sample_uidx], cswritep);
	}
	const uint32_t is_male = IS_SET(sex_male, sample_uidx);
	if (scol_nobs) {
	  cswritep = memcpya(cswritep, nobs_strs[is_male], nobs_slens[is_male]);
	}
	const double cur_variant_ct_recip = variant_ct_recips[is_male];
	if (scol_fmiss) {
	  *cswritep++ = '\t';
	  cswritep = dtoa_g(((double)((int32_t)cur_missing_geno_base)) * cur_variant_ct_recip, cswritep);
	}
	if (scol_fmiss_hh) {
	  *cswritep++ = '\t';
	  cswritep = dtoa_g(((double)((int32_t)(cur_missing_geno_base + sample_hethap_cts[sample_uidx]))) * cur_variant_ct_recip, cswritep);
	}
	append_binary_eoln(&cswritep);
      }
      if (cswrite_close_null(&css, cswritep)) {
	goto write_missingness_reports_ret_WRITE_FAIL;
      }
      bigstack_reset(bigstack_mark);
      LOGPRINTFWW("--missing: Sample missing data report written to %s .\n", outname);
    }
    if (!(missing_rpt_modifier & kfMissingRptSampleOnly)) {
      const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
      unsigned char* overflow_buf;
      char* chr_buf; // includes trailing tab
      if (bigstack_alloc_uc(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen, &overflow_buf) ||
	  bigstack_alloc_c(max_chr_blen, &chr_buf)) {
	goto write_missingness_reports_ret_NOMEM;
      }
      char* outname_end2 = strcpya0(outname_end, ".vmiss");
      if (output_zst) {
	strcpy(outname_end2, ".zst");
      }
      if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
	goto write_missingness_reports_ret_OPEN_FAIL;
      }
      cswritep = (char*)overflow_buf;
      *cswritep++ = '#';
      const uint32_t chr_col = missing_rpt_modifier & kfMissingRptVcolChrom;

      if (chr_col) {
	cswritep = strcpya(cswritep, "CHROM\t");
      }
      if (missing_rpt_modifier & kfMissingRptVcolPos) {
	cswritep = strcpya(cswritep, "POS\t");
      } else {
	variant_bp = nullptr;
      }
      cswritep = strcpya(cswritep, "ID");
      const uint32_t ref_col = missing_rpt_modifier & kfMissingRptVcolRef;
      if (ref_col) {
	cswritep = strcpya(cswritep, "\tREF");
      }
      const uint32_t alt1_col = missing_rpt_modifier & kfMissingRptVcolAlt1;
      if (alt1_col) {
	cswritep = strcpya(cswritep, "\tALT1");
      }
      const uint32_t alt_col = missing_rpt_modifier & kfMissingRptVcolAlt;
      if (alt_col) {
	cswritep = strcpya(cswritep, "\tALT");
      }
      const uint32_t nmiss_col = (missing_rpt_modifier / kfMissingRptVcolNmiss) & 1;
      if (nmiss_col) {
	cswritep = strcpya(cswritep, "\tMISSING_CT");
      }
      const uint32_t nmiss_hh_col = (missing_rpt_modifier / kfMissingRptVcolNmissHh) & 1;
      if (nmiss_hh_col) {
	cswritep = strcpya(cswritep, "\tMISSING_AND_HETHAP_CT");
      }
      const uint32_t hethap_col = (missing_rpt_modifier / kfMissingRptVcolHethap) & 1;
      if (hethap_col) {
	cswritep = strcpya(cswritep, "\tHETHAP_CT");
      }
      const uint32_t nobs_col = (missing_rpt_modifier / kfMissingRptVcolNobs) & 1;
      if (nobs_col) {
	cswritep = strcpya(cswritep, "\tOBS_CT");
      }
      const uint32_t fmiss_col = (missing_rpt_modifier / kfMissingRptVcolFmiss) & 1;
      if (fmiss_col) {
	cswritep = strcpya(cswritep, "\tF_MISS");
      }
      const uint32_t fmiss_hh_col = (missing_rpt_modifier / kfMissingRptVcolFmissHh) & 1;
      if (fmiss_hh_col) {
	cswritep = strcpya(cswritep, "\tF_MISS_AND_HETHAP");
      }
      const uint32_t fhethap_col = (missing_rpt_modifier / kfMissingRptVcolFhethap) & 1;
      if (fhethap_col) {
	cswritep = strcpya(cswritep, "\tF_HETHAP");
      }
      append_binary_eoln(&cswritep);
      char nobs_str[16];
      nobs_str[0] = '\t';
      const int32_t y_code = cip->xymt_codes[kChrOffsetY];
      uint32_t nobs_slen = 0;
      uint32_t variant_uidx = 0;
      uint32_t chr_fo_idx = 0xffffffffU;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t is_y = 2;
      double nobs_recip = 0.0;
      fputs("--missing variant report: 0%", stdout);
      fflush(stdout);
      uint32_t cur_allele_ct = 2;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	if (variant_uidx >= chr_end) {
	  int32_t chr_idx;
	  do {
	    ++chr_fo_idx;
	    chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	    chr_idx = cip->chr_file_order[chr_fo_idx];
	  } while (variant_uidx >= chr_end);
	  char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	  *chr_name_end = '\t';
	  chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	  const uint32_t new_is_y = (chr_idx == y_code);
	  if (new_is_y != is_y) {
	    is_y = new_is_y;
	    const uint32_t cur_nobs = is_y? male_ct : sample_ct;
	    nobs_recip = 1.0 / ((double)((int32_t)cur_nobs));
	    char* nobs_str_end = uint32toa(cur_nobs, &(nobs_str[1]));
	    nobs_slen = (uintptr_t)(nobs_str_end - nobs_str);
	  }
	}
	if (chr_col) {
	  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
	}
	if (variant_bp) {
	  cswritep = uint32toa_x(variant_bp[variant_uidx], '\t', cswritep);
	}
	cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
	uintptr_t variant_allele_idx_base = variant_uidx * 2;
	if (variant_allele_idxs) {
	  variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	  cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
	}
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	  }
	  --cswritep;
	}
	const uint32_t cur_missing_ct = variant_missing_cts[variant_uidx];
	uint32_t cur_hethap_ct = 0;
	if (variant_uidx >= first_hap_uidx) {
	  cur_hethap_ct = variant_hethap_cts[variant_uidx - first_hap_uidx];
	}
	if (nmiss_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(cur_missing_ct, cswritep);
	}
	if (nmiss_hh_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(cur_missing_ct + cur_hethap_ct, cswritep);
	}
	if (hethap_col) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa(cur_hethap_ct, cswritep);
	}
	if (nobs_col) {
	  cswritep = memcpya(cswritep, nobs_str, nobs_slen);
	}
	if (fmiss_col) {
	  *cswritep++ = '\t';
	  cswritep = dtoa_g(((double)((int32_t)cur_missing_ct)) * nobs_recip, cswritep);
	}
	if (fmiss_hh_col) {
	  *cswritep++ = '\t';
	  cswritep = dtoa_g(((double)((int32_t)(cur_missing_ct + cur_hethap_ct))) * nobs_recip, cswritep);
	}
	if (fhethap_col) {
	  *cswritep++ = '\t';
	  cswritep = dtoa_g(((double)((int32_t)cur_hethap_ct)) * nobs_recip, cswritep);
	}
	append_binary_eoln(&cswritep);
	if (cswrite(&css, &cswritep)) {
	  goto write_missingness_reports_ret_WRITE_FAIL;
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
      }
      if (cswrite_close_null(&css, cswritep)) {
	goto write_missingness_reports_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("--missing: Variant missing data report written to %s .\n", outname);
    }
  }
  while (0) {
  write_missingness_reports_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_missingness_reports_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_missingness_reports_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  bigstack_reset(bigstack_mark);
  cswrite_close_cond(&css, cswritep);
  return reterr;
}

// multithread globals
static const uintptr_t* g_variant_include = nullptr;
static const uint32_t* g_founder_raw_geno_cts = nullptr;
static const uint32_t* g_founder_x_male_geno_cts = nullptr;
static const uint32_t* g_founder_x_nosex_geno_cts = nullptr;
static uint32_t* g_variant_uidx_starts = nullptr;
static double* g_hwe_x_pvals = nullptr;
static uint32_t g_x_start = 0;
static uint32_t g_hwe_x_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_hwe_midp = 0;

THREAD_FUNC_DECL compute_hwe_x_pvals_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t* variant_include = g_variant_include;
  const uint32_t* founder_raw_geno_cts = g_founder_raw_geno_cts;
  const uint32_t* founder_x_male_geno_cts = g_founder_x_male_geno_cts;
  const uint32_t* founder_x_nosex_geno_cts = g_founder_x_nosex_geno_cts;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t x_start = g_x_start;
  const uint32_t hwe_x_ct = g_hwe_x_ct;
  const uint32_t hwe_midp = g_hwe_midp;

  // this needs to be aligned with compute_uidx_start_partition()
  const uint32_t variant_idx_end = (hwe_x_ct * (((uint64_t)tidx) + 1)) / calc_thread_ct;
  uint32_t variant_idx = (hwe_x_ct * ((uint64_t)tidx)) / calc_thread_ct;
  
  double* hwe_x_pvals_iter = &(g_hwe_x_pvals[variant_idx]);
  uint32_t variant_uidx = g_variant_uidx_starts[tidx];
  uint32_t pct = 0;
  uint32_t next_print_variant_idx = variant_idx_end;
  if (!tidx) {
    next_print_variant_idx = variant_idx_end / 100;
  }
  uint32_t male_ref_ct = 0;
  uint32_t male_alt_ct = 0;
  for (; variant_idx < variant_idx_end; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    const uint32_t* cur_raw_geno_cts = &(founder_raw_geno_cts[(3 * k1LU) * variant_uidx]);
    uint32_t female_homref_ct = cur_raw_geno_cts[0];
    uint32_t female_refalt_ct = cur_raw_geno_cts[1];
    uint32_t female_altalt_ct = cur_raw_geno_cts[2];
    if (founder_x_male_geno_cts) {
      const uint32_t* cur_male_geno_cts = &(founder_x_male_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
      male_ref_ct = cur_male_geno_cts[0];
      female_homref_ct -= male_ref_ct;
      female_refalt_ct -= cur_male_geno_cts[1];
      male_alt_ct = cur_male_geno_cts[2];
      female_altalt_ct -= male_alt_ct;
    }
    if (founder_x_nosex_geno_cts) {
      const uint32_t* cur_nosex_geno_cts = &(founder_x_nosex_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
      female_homref_ct -= cur_nosex_geno_cts[0];
      female_refalt_ct -= cur_nosex_geno_cts[1];
      female_altalt_ct -= cur_nosex_geno_cts[2];
    }
    *hwe_x_pvals_iter++ = SNPHWEX(female_refalt_ct, female_homref_ct, female_altalt_ct, male_ref_ct, male_alt_ct, hwe_midp);
    if (variant_idx >= next_print_variant_idx) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      pct = (variant_idx * 100LLU) / variant_idx_end;
      printf("\b\b%u%%", pct++);
      fflush(stdout);
      next_print_variant_idx = (pct * ((uint64_t)variant_idx_end)) / 100;
    }
  }
  if (pct > 10) {
    putc_unlocked('\b', stdout);
  }
  THREAD_RETURN;
}

pglerr_t compute_hwe_x_pvals(const uintptr_t* variant_include, const uint32_t* founder_raw_geno_cts, const uint32_t* founder_x_male_geno_cts, const uint32_t* founder_x_nosex_geno_cts, uint32_t x_start, uint32_t hwe_x_ct, uint32_t hwe_midp, uint32_t calc_thread_ct, double** hwe_x_pvals_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    assert(hwe_x_ct);
    if (bigstack_alloc_d(hwe_x_ct, hwe_x_pvals_ptr)) {
      goto compute_hwe_x_pvals_ret_NOMEM;
    }
    bigstack_mark = g_bigstack_base;
    g_hwe_x_pvals = *hwe_x_pvals_ptr;
    
    if (calc_thread_ct > hwe_x_ct) {
      calc_thread_ct = hwe_x_ct;
    }
    pthread_t* threads = (pthread_t*)bigstack_alloc(calc_thread_ct * sizeof(intptr_t));
    if (!threads) {
      goto compute_hwe_x_pvals_ret_NOMEM;
    }
    if (bigstack_alloc_ui(calc_thread_ct, &g_variant_uidx_starts)) {
      goto compute_hwe_x_pvals_ret_NOMEM;
    }
    compute_uidx_start_partition(variant_include, hwe_x_ct, calc_thread_ct, x_start, g_variant_uidx_starts);
    g_variant_include = variant_include;
    g_founder_raw_geno_cts = founder_raw_geno_cts;
    g_founder_x_male_geno_cts = founder_x_male_geno_cts;
    g_founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
    g_calc_thread_ct = calc_thread_ct;
    g_x_start = x_start;
    g_hwe_x_ct = hwe_x_ct;
    g_hwe_midp = hwe_midp;
    LOGPRINTF("Computing chrX Hardy-Weinberg %sp-values... ", hwe_midp? "mid" : "");
    fputs("0%", stdout);
    fflush(stdout);
    if (spawn_threads(compute_hwe_x_pvals_thread, calc_thread_ct, threads)) {
      goto compute_hwe_x_pvals_ret_THREAD_CREATE_FAIL;
    }
    compute_hwe_x_pvals_thread((void*)0);
    join_threads(calc_thread_ct, threads);
    fputs("\b\b", stdout);
    logprint("done.\n");
  }
  while (0) {
  compute_hwe_x_pvals_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  compute_hwe_x_pvals_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t hardy_report(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bp, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint32_t* founder_raw_geno_cts, const uint32_t* founder_x_male_geno_cts, const uint32_t* founder_x_nosex_geno_cts, const double* hwe_x_pvals, uint32_t variant_ct, uint32_t hwe_x_ct, uint32_t max_allele_slen, double output_min_p, hardy_flags_t hardy_modifier, uint32_t nonfounders, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    if (cip->haploid_mask[0] & 1) {
      logerrprint("Error: --hardy is pointless on an all-haploid genome.\n");
      goto hardy_report_ret_INCONSISTENT_INPUT;
    }
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    const uint32_t chr_code_endl = BITCT_TO_WORDCT(chr_code_end);
    unsigned char* overflow_buf;
    uintptr_t* chr_skips;
    if (bigstack_alloc_uc(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen, &overflow_buf) ||
	bigstack_alloc_ul(chr_code_endl, &chr_skips)) {
      goto hardy_report_ret_NOMEM;
    }
    // skip chrX, chrY, chrM here
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    memcpy(chr_skips, cip->haploid_mask, chr_code_endl * sizeof(intptr_t));
    if (mt_code >= 0) {
      set_bit(mt_code, chr_skips);
    }
    const uint32_t chr_skip_ct = popcount_longs(chr_skips, chr_code_endl);
    uint32_t variant_skip_ct = 0;
    uint32_t chr_uidx = 0;
    for (uint32_t chr_skip_idx = 0; chr_skip_idx < chr_skip_ct; ++chr_skip_idx, ++chr_uidx) {
      next_set_unsafe_ck(chr_skips, &chr_uidx);
      if (is_set(cip->chr_mask, chr_uidx)) {
	const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_uidx];
	variant_skip_ct += popcount_bit_idx(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1]);
      }
    }
    if (variant_skip_ct - hwe_x_ct) {
      LOGPRINTF("--hardy: Skipping %u haploid/chrM variant%s.\n", variant_skip_ct - hwe_x_ct, (variant_skip_ct - hwe_x_ct == 1)? "" : "s");
    }
    variant_ct -= variant_skip_ct;
    const uint32_t output_zst = hardy_modifier & kfHardyZs;
    const uint32_t midp = (hardy_modifier / kfHardyMidp) & 1;
    const uint32_t chr_col = hardy_modifier & kfHardyColChrom;
    const uint32_t ref_col = hardy_modifier & kfHardyColRef;
    const uint32_t alt1_col = hardy_modifier & kfHardyColAlt1;
    const uint32_t alt_col = hardy_modifier & kfHardyColAlt;
    const uint32_t gcounts = hardy_modifier & (kfHardyColGcounts | kfHardyColGcount1col);
    const uint32_t gcount_1col = hardy_modifier & kfHardyColGcount1col;
    const char gcount_delim = gcount_1col? ',' : '\t';
    const uint32_t hetfreq_cols = hardy_modifier & kfHardyColHetfreq;
    const uint32_t p_col = hardy_modifier & kfHardyColP;
    if (variant_ct) {
      char* outname_end2 = strcpya0(outname_end, ".hardy");
      if (output_zst) {
	strcpy(outname_end2, ".zst");
      }
      if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
	goto hardy_report_ret_OPEN_FAIL;
      }
      cswritep = (char*)overflow_buf;
      *cswritep++ = '#';

      // includes trailing tab
      char* chr_buf = nullptr;
      if (chr_col) {
	if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
	  goto hardy_report_ret_NOMEM;
	}
	cswritep = strcpya(cswritep, "CHROM\t");
      }
      if (hardy_modifier & kfHardyColPos) {
	cswritep = strcpya(cswritep, "POS\t");
      } else {
	variant_bp = nullptr;
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
      if (gcounts) {
	if (gcount_1col) {
	  cswritep = strcpya(cswritep, "\tGCOUNTS");
	} else {
	  cswritep = strcpya(cswritep, "\tHOM_REF_CT\tHET_REF_CT\tNONREF_CT");
	}
      }
      if (hetfreq_cols) {
	cswritep = strcpya(cswritep, "\tO(HET_REF)\tE(HET_REF)");
      }
      if (p_col) {
	*cswritep++ = '\t';
	if (midp) {
	  cswritep = strcpya(cswritep, "MIDP");
	} else {
	  *cswritep++ = 'P';
	}
      }
      append_binary_eoln(&cswritep);
      uint32_t variant_uidx = 0;
      uint32_t chr_fo_idx = 0xffffffffU;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      printf("--hardy%s%s: 0%%", output_zst? " zs" : "", midp? " midp" : "");
      fflush(stdout);
      uint32_t cur_allele_ct = 2;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	if (chr_col) {
	  if (variant_uidx >= chr_end) {
	    int32_t chr_idx;
	    do {
	      ++chr_fo_idx;
	      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	      chr_idx = cip->chr_file_order[chr_fo_idx];
	    } while ((variant_uidx >= chr_end) || is_set(chr_skips, chr_idx));
	    variant_uidx = next_set_unsafe(variant_include, cip->chr_fo_vidx_start[chr_fo_idx]);
	    char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	    *chr_name_end = '\t';
	    chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	  }
	  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
	}
	if (variant_bp) {
	  cswritep = uint32toa_x(variant_bp[variant_uidx], '\t', cswritep);
	}
	cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
	uintptr_t variant_allele_idx_base = variant_uidx * 2;
	if (variant_allele_idxs) {
	  variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	  cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
	}
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
	      goto hardy_report_ret_WRITE_FAIL;
	    }
	    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	  }
	  --cswritep;
	}
	const uint32_t* cur_geno_cts = &(founder_raw_geno_cts[(3 * k1LU) * variant_uidx]);
	if (gcounts) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa_x(cur_geno_cts[0], gcount_delim, cswritep);
	  cswritep = uint32toa_x(cur_geno_cts[1], gcount_delim, cswritep);
	  cswritep = uint32toa(cur_geno_cts[2], cswritep);
	}
	if (hetfreq_cols) {
	  *cswritep++ = '\t';
	  const uint32_t tot_obs = cur_geno_cts[0] + cur_geno_cts[1] + cur_geno_cts[2];
	  const double tot_obs_recip = 1.0 / (double)((int32_t)tot_obs);
	  cswritep = dtoa_g(((int32_t)cur_geno_cts[1]) * tot_obs_recip, cswritep);
	  *cswritep++ = '\t';
	  const double dbl_ref_freq = (cur_geno_cts[0] * 2 + cur_geno_cts[1]) * tot_obs_recip;
	  const double expected_het_freq = dbl_ref_freq * (1.0 - dbl_ref_freq * 0.5);
	  cswritep = dtoa_g(expected_het_freq, cswritep);
	}
	if (p_col) {
	  // possible todo: multithread this
	  *cswritep++ = '\t';
	  const double hwe_p = SNPHWE2(cur_geno_cts[1], cur_geno_cts[0], cur_geno_cts[2], midp);
	  cswritep = dtoa_g(MAXV(hwe_p, output_min_p), cswritep);
	}
	append_binary_eoln(&cswritep);
	if (cswrite(&css, &cswritep)) {
	  goto hardy_report_ret_WRITE_FAIL;
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
      }
      if (cswrite_close_null(&css, cswritep)) {
	goto hardy_report_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("--hardy%s%s: Autosomal Hardy-Weinberg report (%s) written to %s .\n", output_zst? " zs" : "", midp? " midp" : "", nonfounders? "all samples" : "founders only", outname);
    }
    if (hwe_x_ct) {
      bigstack_reset(chr_skips);
      char* outname_end2 = strcpya0(outname_end, ".hardy.x");
      if (output_zst) {
	strcpy(outname_end2, ".zst");
      }
      if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
	goto hardy_report_ret_OPEN_FAIL;
      }
      cswritep = (char*)overflow_buf;
      *cswritep++ = '#';

      // includes trailing tab
      char x_name_buf[8];
      uint32_t x_name_blen = 0;
      const int32_t x_code = cip->xymt_codes[kChrOffsetX];
      if (chr_col) {
	cswritep = strcpya(cswritep, "CHROM\t");
	char* write_iter = chr_name_write(cip, x_code, x_name_buf);
	*write_iter++ = '\t';
	x_name_blen = (uintptr_t)(write_iter - x_name_buf);
      }
      if (hardy_modifier & kfHardyColPos) {
	cswritep = strcpya(cswritep, "POS\t");
      } else {
	variant_bp = nullptr;
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
      if (gcounts) {
	if (gcount_1col) {
	  cswritep = strcpya(cswritep, "\tGCOUNTS");
	} else {
	  cswritep = strcpya(cswritep, "\tFEMALE_HOM_REF_CT\tFEMALE_HET_REF_CT\tFEMALE_NONREF_CT\tMALE_REF_CT\tMALE_ALT_CT");
	}
      }
      if (hetfreq_cols) {
	cswritep = strcpya(cswritep, "\tO(FEMALE_HET_REF)\tE(FEMALE_HET_REF)");
      }
      const uint32_t sexaf_cols = hardy_modifier & kfHardyColSexaf;
      if (sexaf_cols) {
	cswritep = strcpya(cswritep, "\tFEMALE_REF_FREQ\tMALE_REF_FREQ");
      }
      const uint32_t femalep_col = hardy_modifier & kfHardyColFemalep;
      if (femalep_col) {
	cswritep = strcpya(cswritep, "\tFEMALE_ONLY_");
	if (midp) {
	  cswritep = strcpya(cswritep, "MIDP");
	} else {
	  *cswritep++ = 'P';
	}
      }
      if (p_col) {
	*cswritep++ = '\t';
	if (midp) {
	  cswritep = strcpya(cswritep, "MIDP");
	} else {
	  *cswritep++ = 'P';
	}
      }
      append_binary_eoln(&cswritep);
      fputs("--hardy: Writing chrX results...", stdout);
      fflush(stdout);
      const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
      const uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
      uint32_t variant_uidx = x_start;
      uint32_t cur_allele_ct = 2;
      uint32_t male_ref_ct = 0;
      uint32_t male_alt_ct = 0;
      for (uint32_t variant_idx = 0; variant_idx < hwe_x_ct; ++variant_idx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	cswritep = memcpya(cswritep, x_name_buf, x_name_blen);
	if (variant_bp) {
	  cswritep = uint32toa_x(variant_bp[variant_uidx], '\t', cswritep);
	}
	cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
	uintptr_t variant_allele_idx_base = variant_uidx * 2;
	if (variant_allele_idxs) {
	  variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	  cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
	}
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
	      goto hardy_report_ret_WRITE_FAIL;
	    }
	    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
	  }
	  --cswritep;
	}
	const uint32_t* cur_geno_cts = &(founder_raw_geno_cts[(3 * k1LU) * variant_uidx]);
	uint32_t female_homref_ct = cur_geno_cts[0];
	uint32_t female_refalt_ct = cur_geno_cts[1];
	uint32_t female_altalt_ct = cur_geno_cts[2];
	if (founder_x_male_geno_cts) {
	  const uint32_t* cur_male_geno_cts = &(founder_x_male_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	  male_ref_ct = cur_male_geno_cts[0];
	  female_homref_ct -= male_ref_ct;
	  female_refalt_ct -= cur_male_geno_cts[1];
	  male_alt_ct = cur_male_geno_cts[2];
	  female_altalt_ct -= male_alt_ct;
	}
	if (founder_x_nosex_geno_cts) {
	  const uint32_t* cur_nosex_geno_cts = &(founder_x_nosex_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	  female_homref_ct -= cur_nosex_geno_cts[0];
	  female_refalt_ct -= cur_nosex_geno_cts[1];
	  female_altalt_ct -= cur_nosex_geno_cts[2];
	}
	if (gcounts) {
	  *cswritep++ = '\t';
	  cswritep = uint32toa_x(female_homref_ct, gcount_delim, cswritep);
	  cswritep = uint32toa_x(female_refalt_ct, gcount_delim, cswritep);
	  cswritep = uint32toa_x(female_altalt_ct, gcount_delim, cswritep);
	  cswritep = uint32toa_x(male_ref_ct, gcount_delim, cswritep);
	  cswritep = uint32toa(male_alt_ct, cswritep);
	}
	if (hetfreq_cols || sexaf_cols) {
	  const uint32_t tot_female_obs = female_homref_ct + female_refalt_ct + female_altalt_ct;
	  const double tot_female_obs_recip = 1.0 / (double)((int32_t)tot_female_obs);
	  const double dbl_ref_freq = (female_homref_ct * 2 + female_refalt_ct) * tot_female_obs_recip;
	  const double ref_freq = dbl_ref_freq * 0.5;
	  if (hetfreq_cols) {
	    *cswritep++ = '\t';
	    cswritep = dtoa_g(((int32_t)female_refalt_ct) * tot_female_obs_recip, cswritep);
	    *cswritep++ = '\t';
	    const double expected_het_freq = dbl_ref_freq * (1.0 - ref_freq);
	    cswritep = dtoa_g(expected_het_freq, cswritep);
	  }
	  if (sexaf_cols) {
	    *cswritep++ = '\t';
	    cswritep = dtoa_g(ref_freq, cswritep);
	    *cswritep++ = '\t';
	    const double male_ref_freq = ((double)((int32_t)male_ref_ct)) / ((double)((int32_t)(male_ref_ct + male_alt_ct)));
	    cswritep = dtoa_g(male_ref_freq, cswritep);
	  }
	}
	if (femalep_col) {
	  *cswritep++ = '\t';
	  const double female_hwe_p = SNPHWE2(female_refalt_ct, female_homref_ct, female_altalt_ct, midp);
	  cswritep = dtoa_g(MAXV(female_hwe_p, output_min_p), cswritep);
	}
	if (p_col) {
	  *cswritep++ = '\t';
	  cswritep = dtoa_g(MAXV(hwe_x_pvals[variant_idx], output_min_p), cswritep);
	}
	append_binary_eoln(&cswritep);
	if (cswrite(&css, &cswritep)) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
      }
      if (cswrite_close_null(&css, cswritep)) {
	goto hardy_report_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("--hardy%s%s: chrX Hardy-Weinberg report (%s) written to %s .\n", output_zst? " zs" : "", midp? " midp" : "", nonfounders? "all samples" : "founders only", outname);
    }
  }
  while (0) {
  hardy_report_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  hardy_report_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  hardy_report_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  hardy_report_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t write_snplist(const uintptr_t* variant_include, char** variant_ids, uint32_t variant_ct, uint32_t output_zst, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    unsigned char* overflow_buf;
    if (bigstack_alloc_uc(kCompressStreamBlock + kMaxIdSlen + 2, &overflow_buf)) {
      goto write_snplist_ret_NOMEM;
    }
    char* outname_end2 = strcpy(outname_end, ".snplist");
    if (output_zst) {
      strcpy(outname_end2, ".zst");
    }
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto write_snplist_ret_OPEN_FAIL;
    }
    cswritep = (char*)overflow_buf;
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      append_binary_eoln(&cswritep);
      if (cswrite(&css, &cswritep)) {
	goto write_snplist_ret_WRITE_FAIL;
      }
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto write_snplist_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--write-snplist%s: Variant IDs written to %s .\n", output_zst? " zs" : "", outname);
  }
  while (0) {
  write_snplist_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_snplist_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_snplist_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
