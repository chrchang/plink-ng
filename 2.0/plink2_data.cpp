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
#include "plink2_data.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"
#include "plink2_random.h"

#include "bgzf.h"
#include "zstd/lib/zstd.h"

#include <time.h>

#ifdef __cplusplus
namespace plink2 {
#endif

void init_plink1_dosage(plink1_dosage_info_t* plink1_dosage_info_ptr) {
  plink1_dosage_info_ptr->flags = kfPlink1Dosage0;
  fill_uint_zero(3, plink1_dosage_info_ptr->skips);
  plink1_dosage_info_ptr->chr_col_idx = 0xffffffffU;
  plink1_dosage_info_ptr->pos_col_idx = 0xffffffffU;
}

void init_gendummy(gendummy_info_t* gendummy_info_ptr) {
  gendummy_info_ptr->flags = kfGenDummy0;
  gendummy_info_ptr->pheno_ct = 1;
  gendummy_info_ptr->geno_mfreq = 0.0;
  gendummy_info_ptr->pheno_mfreq = 0.0;
  gendummy_info_ptr->dosage_freq = 0.0;
}

pglerr_t write_map_or_bim(const char* outname, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst) {
  // set max_allele_slen to zero for .map
  // allele_dosages must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    // includes trailing tab
    char* chr_buf;

    unsigned char* overflow_buf;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
	bigstack_alloc_uc(kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen, &overflow_buf)) {
      goto write_map_or_bim_ret_NOMEM;
    }
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto write_map_or_bim_ret_OPEN_FAIL;
    }
    cswritep = (char*)overflow_buf;

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t variant_uidx = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	char* chr_name_end = chr_name_write(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
	*chr_name_end = delim;
	chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
      }
      cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], delim);
      if (!variant_cms) {
	*cswritep++ = '0';
      } else {
	cswritep = dtoa_g_p8(variant_cms[variant_uidx], cswritep);
      }
      *cswritep++ = delim;
      cswritep = uint32toa(variant_bps[variant_uidx], cswritep);
      if (max_allele_slen) {
	*cswritep++ = delim;
	const uintptr_t variant_allele_idx_base = variant_allele_idxs? variant_allele_idxs[variant_uidx] : (variant_uidx * 2);
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
	// note that VCF ref allele corresponds to A2, not A1
	if (!refalt1_select) {
	  // needs to be revised in multiallelic case
	  if ((!allele_dosages) || allele_dosages[1 + variant_allele_idx_base]) {
	    cswritep = strcpya(cswritep, cur_alleles[1]);
	  } else {
	    *cswritep++ = output_missing_geno_char;
	  }
	  *cswritep++ = delim;
	  cswritep = strcpya(cswritep, cur_alleles[0]);
	} else {
	  const alt_allele_ct_t* cur_refalt1_select = &(refalt1_select[variant_uidx * 2]);
	  if ((!allele_dosages) || allele_dosages[cur_refalt1_select[1] + variant_allele_idx_base]) {
	    cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[1]]);
	  } else {
	    *cswritep++ = output_missing_geno_char;
	  }
	  *cswritep++ = delim;
	  cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[0]]);
	}
      }
      append_binary_eoln(&cswritep);
      if (cswrite(&css, &cswritep)) {
	goto write_map_or_bim_ret_WRITE_FAIL;
      }
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto write_map_or_bim_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_map_or_bim_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_map_or_bim_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_map_or_bim_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t pvar_info_reload_header(const char* pvar_info_reload, gzFile* gz_pvar_reload_ptr, char** loadbuf_ptr, uintptr_t* loadbuf_size_ptr, uint32_t* info_col_idx_ptr) {
  pglerr_t reterr = gzopen_read_checked(pvar_info_reload, gz_pvar_reload_ptr);
  if (reterr) {
    return reterr;
  }
  uintptr_t loadbuf_size = bigstack_left();
  if (loadbuf_size > kMaxLongLine) {
    loadbuf_size = kMaxLongLine;
  } else if (loadbuf_size <= kMaxMediumLine) {
    return kPglRetNomem;
  }
  char* loadbuf = (char*)g_bigstack_base;
  loadbuf[loadbuf_size - 1] = ' ';
  char* loadbuf_iter;
  do {
    // this is a reload, so no need to validate
    if (!gzgets(*gz_pvar_reload_ptr, loadbuf, loadbuf_size)) {
      return kPglRetReadFail;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == kMaxLongLine) {
	return kPglRetReadFail;
      }
      return kPglRetNomem;
    }
    loadbuf_iter = skip_initial_spaces(loadbuf);
  } while (memcmp(loadbuf_iter, "#CHROM", 6));
  uint32_t info_col_idx = 0;
  do {
    loadbuf_iter = next_token(loadbuf_iter);
    ++info_col_idx;
  } while (memcmp(loadbuf_iter, "INFO", 4) || (((unsigned char)loadbuf_iter[4]) > 32));
  *loadbuf_ptr = loadbuf;
  *loadbuf_size_ptr = loadbuf_size;
  *info_col_idx_ptr = info_col_idx;
  return kPglRetSuccess;
}

void pvar_info_write(char* info_token, uint32_t xheader_info_pr, uint32_t is_pr, char** write_iter_ptr) {
  char* info_token_end = token_endnn(info_token);
  uint32_t info_token_slen = (uintptr_t)(info_token_end - info_token);
  char* info_token_pr = nullptr;
  if (xheader_info_pr) {
    info_token_pr = pr_in_info_token(info_token_slen, info_token);
  }
  char* write_iter = *write_iter_ptr;
  if (is_pr || (!info_token_pr))  {
    write_iter = memcpya(write_iter, info_token, info_token_slen);
    if (is_pr && (!info_token_pr)) {
      if ((info_token_slen == 1) && (info_token[0] == '.')) {
	write_iter[-1] = 'P';
	*write_iter++ = 'R';
      } else {
	write_iter = memcpyl3a(write_iter, ";PR");
      }
    }
  } else {
    // currently only possible with --real-ref-alleles
    if (info_token_pr == info_token) {
      if (info_token_slen == 2) {
	*write_iter++ = '.';
      } else {
	write_iter = memcpya(write_iter, &(info_token[3]), info_token_slen - 3);
      }
    } else {
      write_iter = memcpya(write_iter, info_token, ((uintptr_t)(info_token_pr - info_token)) - 1);
      char* pr_end = &(info_token_pr[2]);
      write_iter = memcpya(write_iter, pr_end, (uintptr_t)(info_token_end - pr_end));
    }
  }
  *write_iter_ptr = write_iter;
}

pglerr_t pvar_info_reload_and_write(uintptr_t loadbuf_size, uint32_t xheader_info_pr, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, gzFile gz_pvar_reload, char** write_iter_ptr, uint32_t* gz_variant_uidx_ptr, char* loadbuf) {
  uint32_t gz_variant_uidx = *gz_variant_uidx_ptr;
  char* loadbuf_first_token;
  do {
    do {
      if (!gzgets(gz_pvar_reload, loadbuf, loadbuf_size)) {
	return kPglRetReadFail;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  return kPglRetReadFail;
	}
	return kPglRetNomem;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));
    ++gz_variant_uidx;
  } while (gz_variant_uidx <= variant_uidx);
  char* info_token = next_token_mult(loadbuf_first_token, info_col_idx);
  *gz_variant_uidx_ptr = gz_variant_uidx;
  pvar_info_write(info_token, xheader_info_pr, is_pr, write_iter_ptr);
  return kPglRetSuccess;
}

void append_chrset_line(const chr_info_t* cip, char** write_iter_ptr) {
  char* write_iter = strcpya(*write_iter_ptr, "##chrSet=<");
  if (!(cip->haploid_mask[0] & 1)) {
    write_iter = strcpya(write_iter, "autosomePairCt=");
    write_iter = uint32toa(cip->autosome_ct, write_iter);
    if (cip->xymt_codes[kChrOffsetX] >= 0) {
      write_iter = strcpya(write_iter, ",X");
    }
    if (cip->xymt_codes[kChrOffsetY] >= 0) {
      write_iter = strcpya(write_iter, ",Y");
    }
    if (cip->xymt_codes[kChrOffsetXY] >= 0) {
      write_iter = strcpya(write_iter, ",XY");
    }
    if (cip->xymt_codes[kChrOffsetMT] >= 0) {
      write_iter = strcpya(write_iter, ",M");
    }
    if (cip->xymt_codes[kChrOffsetPAR1] >= 0) {
      write_iter = strcpya(write_iter, ",PAR1");
    }
    if (cip->xymt_codes[kChrOffsetPAR2] >= 0) {
      write_iter = strcpya(write_iter, ",PAR2");
    }
  } else {
    write_iter = strcpya(write_iter, "haploidAutosomeCt=");
    write_iter = uint32toa(cip->autosome_ct, write_iter);
  }
  *write_iter++ = '>';
  *write_iter_ptr = write_iter;
  append_binary_eoln(write_iter_ptr);
}

pglerr_t write_pvar(const char* outname, const char* xheader, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, char** filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, pvar_psam_t pvar_psam_modifier) {
  // allele_dosages must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  gzFile gz_pvar_reload = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    // includes trailing tab
    char* chr_buf;

    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    unsigned char* overflow_buf;
    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
	bigstack_alloc_uc(overflow_buf_size, &overflow_buf) ||
	bigstack_alloc_ul(BITCT_TO_WORDCT(kPglMaxAltAlleleCt), &allele_include)) {
      goto write_pvar_ret_NOMEM;
    }
    const uint32_t output_zst = (pvar_psam_modifier / kfPvarZs) & 1;
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto write_pvar_ret_OPEN_FAIL;
    }
    cswritep = (char*)overflow_buf;
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_modifier & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      uint32_t widx;
      for (widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & nonref_flags[widx]) {
	  break;
	}
      }
      if (widx == raw_variant_ctl) {
	write_info_pr = 0;
      }
    }
    write_info_pr = write_info_pr && write_info;

    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t info_col_idx = 0; // could save this during first load instead
    if (pvar_psam_modifier & kfPvarColXheader) {
      if (csputs_std(xheader, xheader_blen, &css, &cswritep)) {
	goto write_pvar_ret_WRITE_FAIL;
      }
      if (write_info_pr && (!xheader_info_pr)) {
	cswritep = strcpya(cswritep, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
      }
      if (pvar_info_reload) {
	reterr = pvar_info_reload_header(pvar_info_reload, &gz_pvar_reload, &loadbuf, &loadbuf_size, &info_col_idx);
	if (reterr) {
	  goto write_pvar_ret_1;
	}
      }
    }
    if (cip->chrset_source) {
      append_chrset_line(cip, &cswritep);
    }
    cswritep = strcpya(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_modifier & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_modifier & kfPvarColMaybequal) && qual_present) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & qual_present[widx]) {
	  write_qual = 1;
	  break;
	}
      }
    }
    if (write_qual) {
      cswritep = strcpya(cswritep, "\tQUAL");
    }
    
    uint32_t write_filter = 0;
    if (pvar_psam_modifier & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_modifier & kfPvarColMaybefilter) && filter_present) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & filter_present[widx]) {
	  write_filter = 1;
	  break;
	}
      }
    }
    if (write_filter) {
      cswritep = strcpya(cswritep, "\tFILTER");
    }

    if (write_info) {
      cswritep = strcpya(cswritep, "\tINFO");
    }
    
    uint32_t write_cm = 0;
    if (pvar_psam_modifier & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_modifier & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
	// nonzero_cm_present check was performed
	write_cm = 1;
      } else {
	uint32_t variant_uidx = 0;
	for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	  next_set_unsafe_ck(variant_include, &variant_uidx);
	  if (variant_cms[variant_uidx] != 0.0) {
	    write_cm = 1;
	    break;
	  }
	}
      }
    }
    if (write_cm) {
      cswritep = memcpyl3a(cswritep, "\tCM");
    }
    append_binary_eoln(&cswritep);

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t gz_variant_uidx = 0;
    uint32_t variant_uidx = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
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
      cswritep = uint32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
      uintptr_t variant_allele_idx_base;
      if (!variant_allele_idxs) {
	variant_allele_idx_base = variant_uidx * 2;
      } else {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
	alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
      }
      cswritep = strcpyax(cswritep, cur_alleles[ref_allele_idx], '\t');
      if ((!allele_dosages) || allele_dosages[variant_allele_idx_base + alt1_allele_idx]) {
        cswritep = strcpya(cswritep, cur_alleles[alt1_allele_idx]);
      } else {
	*cswritep++ = output_missing_geno_char;
      }
      if (cswrite(&css, &cswritep)) {
	goto write_pvar_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
	fill_all_bits(cur_allele_ct, allele_include);
	CLEAR_BIT(ref_allele_idx, allele_include);
	CLEAR_BIT(alt1_allele_idx, allele_include);
        uint32_t cur_allele_uidx = 0;
	uint32_t alt_allele_idx = 2;
	do {
	  *cswritep++ = ',';
	  next_set_unsafe_ck(allele_include, &cur_allele_uidx);
	  cswritep = strcpya(cswritep, cur_alleles[cur_allele_uidx++]);
	  if (cswrite(&css, &cswritep)) {
	    goto write_pvar_ret_WRITE_FAIL;
	  }
	} while (++alt_allele_idx < cur_allele_ct);
      }

      if (write_qual) {
	*cswritep++ = '\t';
	if (!IS_SET(qual_present, variant_uidx)) {
	  *cswritep++ = '.';
	} else {
	  cswritep = ftoa_g(quals[variant_uidx], cswritep);
	}
      }

      if (write_filter) {
	*cswritep++ = '\t';
	if (!IS_SET(filter_present, variant_uidx)) {
	  *cswritep++ = '.';
	} else if (!IS_SET(filter_npass, variant_uidx)) {
	  cswritep = strcpya(cswritep, "PASS");
	} else {
	  cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
	}
      }

      if (write_info) {
	*cswritep++ = '\t';
	const uint32_t is_pr = all_nonref || (nonref_flags && IS_SET(nonref_flags, variant_uidx));
	if (gz_pvar_reload) {
	  reterr = pvar_info_reload_and_write(loadbuf_size, xheader_info_pr, info_col_idx, variant_uidx, is_pr, gz_pvar_reload, &cswritep, &gz_variant_uidx, loadbuf);
	  if (reterr) {
	    goto write_pvar_ret_1;
	  }
	} else {
	  if (is_pr) {
	    cswritep = strcpya(cswritep, "PR");
	  } else {
	    *cswritep++ = '.';
	  }
	}
      }
      
      if (write_cm) {
        *cswritep++ = '\t';
	if (!variant_cms) {
	  *cswritep++ = '0';
	} else {
	  cswritep = dtoa_g_p8(variant_cms[variant_uidx], cswritep);
	}
      }
      append_binary_eoln(&cswritep);
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto write_pvar_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_pvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_pvar_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_pvar_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 write_pvar_ret_1:
  cswrite_close_cond(&css, cswritep);
  gzclose_cond(gz_pvar_reload);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t write_fam(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, char delim) {
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto write_fam_ret_OPEN_FAIL;
    }
    uintptr_t* pheno_nm = nullptr;
    uintptr_t* pheno_cc = nullptr;
    double* pheno_qt = nullptr;
    // .fam files don't support categorical phenotypes
    const uint32_t pheno_idx = first_cc_or_qt_pheno_idx(pheno_cols, pheno_ct);
    if (pheno_idx != 0xffffffffU) {
      const pheno_dtype_t type_code = pheno_cols[pheno_idx].type_code;
      pheno_nm = pheno_cols[pheno_idx].nonmiss;
      if (type_code == kPhenoDtypeCc) {
	pheno_cc = pheno_cols[pheno_idx].data.cc;
      } else {
	pheno_qt = pheno_cols[pheno_idx].data.qt;
      }
    }
    const char* legacy_output_missing_pheno = g_legacy_output_missing_pheno;
    const uint32_t lomp_slen = strlen(legacy_output_missing_pheno);

    uintptr_t sample_uidx = 0;
    uint32_t sample_uidx2 = 0;
    char* textbuf = g_textbuf;
    char* write_iter = textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      if (!new_sample_idx_to_old) {
	next_set_ul_unsafe_ck(sample_include, &sample_uidx);
      } else {
	do {
	  sample_uidx = new_sample_idx_to_old[sample_uidx2++];
	} while (!IS_SET(sample_include, sample_uidx));
      }
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      if (delim == '\t') {
	write_iter = strcpya(write_iter, cur_sample_id);
      } else {
	const char* fid_end = (const char*)rawmemchr(cur_sample_id, '\t');
	write_iter = memcpyax(write_iter, cur_sample_id, (uintptr_t)(fid_end - cur_sample_id), delim);
	write_iter = strcpya(write_iter, &(fid_end[1]));
      }
      *write_iter++ = delim;
      write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), delim);
      write_iter = strcpyax(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]), delim);
      *write_iter++ = sexchar(sex_nm, sex_male, sample_uidx);
      *write_iter++ = delim;
      if ((!pheno_nm) || (!IS_SET(pheno_nm, sample_uidx))) {
	write_iter = memcpya(write_iter, legacy_output_missing_pheno, lomp_slen);
      } else if (pheno_cc) {
	// do we want to allow user to force 0/1 output?
	*write_iter++ = '1' + IS_SET(pheno_cc, sample_uidx);
      } else {
	write_iter = dtoa_g(pheno_qt[sample_uidx], write_iter);
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= textbuf_flush) {
	if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	  goto write_fam_ret_WRITE_FAIL;
	}
	write_iter = textbuf;
      }
    }
    if (write_iter != textbuf) {
      if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	goto write_fam_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto write_fam_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_fam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_fam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

uint32_t is_parental_info_present(const uintptr_t* sample_include, const char* paternal_ids, const char* maternal_ids, uint32_t sample_ct, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen) {
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include, &sample_uidx);
    if (memcmp(&(paternal_ids[sample_uidx * max_paternal_id_blen]), "0", 2) || memcmp(&(maternal_ids[sample_uidx * max_maternal_id_blen]), "0", 2)) {
      return 1;
    }
  }
  return 0;
}

char* append_pheno_str(const pheno_col_t* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char* write_iter) {
  const pheno_dtype_t type_code = pheno_col->type_code;
  if (type_code <= kPhenoDtypeQt) {
    if (!is_set(pheno_col->nonmiss, sample_uidx)) {
      write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
    } else if (type_code == kPhenoDtypeCc) {
      *write_iter++ = '1' + is_set(pheno_col->data.cc, sample_uidx);
    } else {
      write_iter = dtoa_g(pheno_col->data.qt[sample_uidx], write_iter);
    }
  } else {
    write_iter = strcpya(write_iter, pheno_col->category_names[pheno_col->data.cat[sample_uidx]]);
  }
  return write_iter;
}

pglerr_t write_psam(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, pvar_psam_t pvar_psam_modifier) {
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto write_psam_ret_OPEN_FAIL;
    }
    const char* output_missing_pheno = g_output_missing_pheno;
    const uint32_t omp_slen = strlen(output_missing_pheno);
    
    char* textbuf = g_textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);

    const uint32_t write_sid = sid_col_required(sample_include, sids, sample_ct, max_sid_blen, pvar_psam_modifier / kfPsamColMaybesid);
    uint32_t write_parents = 0;
    if (pvar_psam_modifier & kfPsamColParents) {
      write_parents = 1;
    } else if (pvar_psam_modifier & kfPsamColMaybeparents) {
      write_parents = is_parental_info_present(sample_include, paternal_ids, maternal_ids, sample_ct, max_paternal_id_blen, max_maternal_id_blen);
    }
    const uint32_t write_sex = (pvar_psam_modifier / kfPsamColSex) & 1;
    const uint32_t write_empty_pheno = (pvar_psam_modifier & kfPsamColPheno1) && (!pheno_ct);
    const uint32_t write_phenos = (pvar_psam_modifier & (kfPsamColPheno1 | kfPsamColPhenos)) && pheno_ct;
    if (write_phenos && (!(pvar_psam_modifier & kfPsamColPhenos))) {
      pheno_ct = 1;
    }
    char* write_iter = strcpya(textbuf, "#FID\tIID");
    if (write_sid) {
      write_iter = strcpya(write_iter, "\tSID");
    }
    if (write_parents) {
      write_iter = strcpya(write_iter, "\tPAT\tMAT");
    }
    if (write_sex) {
      write_iter = strcpya(write_iter, "\tSEX");
    }
    if (write_phenos) {
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	*write_iter++ = '\t';
	const char* cur_pheno_name = &(pheno_names[pheno_idx * max_pheno_name_blen]);
	const uint32_t cur_pheno_name_slen = strlen(cur_pheno_name);
	if ((cur_pheno_name_slen == 3) && (!memcmp(cur_pheno_name, "SEX", 3))) {
	  if (write_sex) {
	    logerrprint("Error: .psam file cannot have both a regular SEX column and a phenotype named\n'SEX'.  Exclude or rename one of these columns.\n");
	    goto write_psam_ret_INCONSISTENT_INPUT;
	  }
	  // does this phenotype column conform to the SEX column format?
	  // case/control is always ok, but quantitative or categorical needs
	  // to be checked
	  const pheno_col_t* sex_col = &(pheno_cols[pheno_idx]);
	  if (sex_col->type_code != kPhenoDtypeCc) {
	    // could bitwise-and sample_include and pheno_nm before the loop
	    const uintptr_t* pheno_nm = sex_col->nonmiss;
	    uint32_t sample_uidx = 0;
	    if (sex_col->type_code == kPhenoDtypeQt) {
	      const double* pheno_vals = sex_col->data.qt;
	      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
		next_set_unsafe_ck(sample_include, &sample_uidx);
		if (is_set(pheno_nm, sample_uidx)) {
		  const double dxx = pheno_vals[sample_uidx];
		  // tolerate '-9' and '0' as missing values, and anything in
		  // [1, 2] (could be reasonable to represent XXY, etc. with
		  // decimals).
		  if (((dxx < 1.0) && (dxx != -9.0) && (dxx != 0.0)) || (dxx > 2.0)) {
		    logerrprint("Error: .psam SEX values are expected to be in {-9, 0, 1, 2}.\n");
		    goto write_psam_ret_INCONSISTENT_INPUT;
		  }
		}
	      }
	    } else {
	      assert(sex_col->type_code == kPhenoDtypeCat);
	      const uint32_t nonnull_cat_ct = sex_col->nonnull_category_ct;
	      if (nonnull_cat_ct) {
		char** cur_category_names = sex_col->category_names;
		// tolerate 'M' and 'm' being present simultaneously, etc.
		uint32_t male_cat_idx1 = 0;
		uint32_t male_cat_idx2 = 0;
		uint32_t female_cat_idx1 = 0;
		uint32_t female_cat_idx2 = 0;
		for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
		  const char* cur_cat_name = cur_category_names[cat_idx];
		  if (!cur_cat_name[1]) {
		    uint32_t first_char_code = (unsigned char)cur_cat_name[0];
		    first_char_code &= 0xdf;
		    if (first_char_code == 70) {
		      if (!female_cat_idx1) {
			female_cat_idx1 = cat_idx;
		      } else {
			female_cat_idx2 = cat_idx;
		      }
		    } else if (first_char_code == 77) {
		      if (!male_cat_idx1) {
			male_cat_idx1 = cat_idx;
		      } else {
			male_cat_idx2 = cat_idx;
		      }
		    }
		  }
		}
		if ((uint32_t)((male_cat_idx1 != 0) + (male_cat_idx2 != 0) + (female_cat_idx1 != 0) + (female_cat_idx2 != 0)) < nonnull_cat_ct) {
		  const uint32_t* pheno_vals = sex_col->data.cat;
		  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
		    next_set_unsafe_ck(sample_include, &sample_uidx);
		    if (is_set(pheno_nm, sample_uidx)) {
		      const uint32_t cur_cat_idx = pheno_vals[sample_uidx];
		      if ((cur_cat_idx != male_cat_idx1) && (cur_cat_idx != female_cat_idx1) && (cur_cat_idx != male_cat_idx2) && (cur_cat_idx != female_cat_idx2)) {
			logerrprint("Error: .psam SEX values are expected to be in {'F', 'f', 'M', 'm'}.\n");
			goto write_psam_ret_INCONSISTENT_INPUT;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	write_iter = memcpya(write_iter, cur_pheno_name, cur_pheno_name_slen);
	if (write_iter >= textbuf_flush) {
	  if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	    goto write_psam_ret_WRITE_FAIL;
	  }
	  write_iter = textbuf;
	}
      }
    } else if (write_empty_pheno) {
      write_iter = strcpya(write_iter, "\tPHENO1");
    }
    append_binary_eoln(&write_iter);

    uintptr_t sample_uidx = 0;
    uint32_t sample_uidx2 = 0;
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      if (!new_sample_idx_to_old) {
	next_set_ul_unsafe_ck(sample_include, &sample_uidx);
      } else {
	do {
	  sample_uidx = new_sample_idx_to_old[sample_uidx2++];
	} while (!IS_SET(sample_include, sample_uidx));
      }
      write_iter = strcpya(write_iter, &(sample_ids[max_sample_id_blen * sample_uidx]));
      if (write_sid) {
	*write_iter++ = '\t';
	if (sids) {
	  write_iter = strcpya(write_iter, &(sids[max_sid_blen * sample_uidx]));
	} else {
	  *write_iter++ = '0';
	}
      }
      if (write_parents) {
	*write_iter++ = '\t';
	write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), '\t');
	write_iter = strcpya(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]));
      }
      if (write_sex) {
	*write_iter++ = '\t';
	if (IS_SET(sex_nm, sample_uidx)) {
	  *write_iter++ = '2' - IS_SET(sex_male, sample_uidx);
	} else {
	  // this is better than '0' since it allows the raw column to be used
	  // as --covar input
	  // (can't do this for .fam export, though: not worth the
	  // compatibility issues)
	  write_iter = strcpya(write_iter, "NA");
	}
      }
      if (write_phenos) {
	for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  *write_iter++ = '\t';
	  write_iter = append_pheno_str(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
	  if (write_iter >= textbuf_flush) {
	    if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	      goto write_psam_ret_WRITE_FAIL;
	    }
	    write_iter = textbuf;
	  }
	}
      } else {
	if (write_empty_pheno) {
	  *write_iter++ = '\t';
	  write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
	}
	if (write_iter >= textbuf_flush) {
	  if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	    goto write_psam_ret_WRITE_FAIL;
	  }
	  write_iter = textbuf;
	}	
      }
      append_binary_eoln(&write_iter);
    }
    if (write_iter != textbuf) {
      if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	goto write_psam_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto write_psam_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_psam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_psam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  write_psam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}


pglerr_t vcf_sample_line(const char* preexisting_psamname, const char* const_fid, uint32_t double_id, fam_col_t fam_cols, char id_delim, char idspace_to, char flag_char, char* sample_line_first_id, char* outname, char* outname_end, uintptr_t* sample_ct_ptr) {
  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    uintptr_t const_fid_len = 0;
    if (const_fid) {
      const_fid_len = strlen(const_fid);
    } else if ((!double_id) && (!id_delim)) {
      // default: --double-id + --id-delim
      double_id = 1;
      id_delim = '_';
    }
    const uint32_t double_or_const_fid = double_id || const_fid;
    if (id_delim != ' ') {
      char* sample_line_iter = strchr(sample_line_first_id, ' ');
      if (sample_line_iter) {
	if (!idspace_to) {
	  logerrprint("Error: VCF/BCF2 sample ID contains space(s).  Use --idspace-to to convert them\nto another character, or \"--id-delim ' '\" to interpret the spaces as FID/IID\ndelimiters.\n");
	  goto vcf_sample_line_ret_INCONSISTENT_INPUT;
	}
	do {
	  *sample_line_iter = idspace_to;
	  sample_line_iter = strchr(&(sample_line_iter[1]), ' ');
	} while (sample_line_iter);
      }
    }
    char* sample_line_iter = sample_line_first_id;
    char* textbuf = g_textbuf;
    uintptr_t sample_ct = 0;
    if (!preexisting_psamname) {
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto vcf_sample_line_ret_OPEN_FAIL;
      }
      char* write_iter = strcpya(textbuf, "#FID\tIID");
      uint32_t sid_present = 0;
      if (id_delim) {
	while (((unsigned char)sample_line_iter[0]) >= ' ') {
	  char* token_end = strchr(sample_line_iter, '\t');
	  if (!token_end) {
	    token_end = next_prespace(sample_line_iter);
	  }
	  char* first_delim = (char*)memchr(sample_line_iter, (unsigned char)id_delim, (uintptr_t)(token_end - sample_line_iter));
	  if (first_delim) {
	    sample_line_iter = &(first_delim[1]);
	    if (memchr(sample_line_iter, (unsigned char)id_delim, (uintptr_t)(token_end - sample_line_iter)) != nullptr) {
	      sid_present = 1;
	      write_iter = strcpya(write_iter, "\tSID");
	      break;
	    }
	  }
	  if (*token_end != '\t') {
	    break;
	  }
	  sample_line_iter = &(token_end[1]);
	}
	sample_line_iter = sample_line_first_id;
      }
      write_iter = strcpya(write_iter, "\tSEX");
      append_binary_eoln(&write_iter);
      char* textbuf_flush = &(textbuf[kMaxMediumLine]);
      while (((unsigned char)sample_line_iter[0]) >= ' ') {
	++sample_ct;
	char* token_end = strchr(sample_line_iter, '\t');
	if (!token_end) {
	  token_end = next_prespace(sample_line_iter);
	}
	const uint32_t token_slen = (uintptr_t)(token_end - sample_line_iter);
	if ((*sample_line_iter == '0') && (token_slen == 1)) {
	  logerrprint("Error: Sample ID cannot be '0'.\n");
	  goto vcf_sample_line_ret_MALFORMED_INPUT;
	}
	if (id_delim) {
	  if (*sample_line_iter == id_delim) {
	    sprintf(g_logbuf, "Error: '%c' at beginning of sample ID.\n", id_delim);
	    goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
	  }
	  if (sample_line_iter[token_slen - 1] == id_delim) {
	    sprintf(g_logbuf, "Error: '%c' at end of sample ID.\n", id_delim);
	    goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
	  }
	  char* first_delim = (char*)memchr(sample_line_iter, (unsigned char)id_delim, token_slen);
	  if (!first_delim) {
	    if (double_or_const_fid) {
	      goto vcf_sample_line_nopsam_one_id;
	    }
	    sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
	    goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
	  }
	  char* iid_start = &(first_delim[1]);
	  char* iid_end = (char*)memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(token_end - iid_start));
	  const char* sid_start = &(g_one_char_strs[96]);
	  uint32_t sid_slen = 1;
	  if (iid_end) {
	    if (iid_start == iid_end) {
	      sprintf(g_logbuf, "Error: Consecutive instances of '%c' in sample ID.\n", id_delim);
	      goto vcf_sample_line_ret_INCONSISTENT_INPUT_DELIM;
	    }
	    sid_start = &(iid_end[1]);
	    sid_slen = (uintptr_t)(token_end - sid_start);
	    if (memchr(sid_start, (unsigned char)id_delim, sid_slen)) {
	      sprintf(g_logbuf, "Error: More than two instances of '%c' in sample ID.\n", id_delim);
	      goto vcf_sample_line_ret_INCONSISTENT_INPUT_DELIM;
	    }
	    if (sid_slen > kMaxIdSlen) {
	      logerrprint("Error: SIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
	      goto vcf_sample_line_ret_MALFORMED_INPUT;
	    }
	  } else {
	    iid_end = token_end;
	  }
	  const uint32_t fid_slen = (uintptr_t)(first_delim - sample_line_iter);
	  if (fid_slen > kMaxIdSlen) {
	    // strictly speaking, you could have e.g. a 20k char ID which
	    // splits into a valid FID/IID pair with the right delimiter, but
	    // you're not supposed to have sample IDs anywhere near that length
	    // so I'll classify this as MalformedInput.
	    goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  write_iter = memcpyax(write_iter, sample_line_iter, fid_slen, '\t');
	  const uint32_t iid_slen = (uintptr_t)(iid_end - iid_start);
	  if ((*iid_start == '0') && (iid_slen == 1)) {
	    logerrprint("Error: Sample ID induces an invalid IID of '0'.\n");
	    goto vcf_sample_line_ret_INCONSISTENT_INPUT;
	  }
	  if (iid_slen > kMaxIdSlen) {
	    goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  write_iter = memcpya(write_iter, iid_start, iid_slen);
	  if (sid_present) {
	    *write_iter++ = '\t';
	    write_iter = memcpya(write_iter, sid_start, sid_slen);
	  }
	} else {
	vcf_sample_line_nopsam_one_id:
	  if (token_slen > kMaxIdSlen) {
	    goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  if (double_id) {
	    write_iter = memcpya(write_iter, sample_line_iter, token_slen);
	  } else {
	    write_iter = memcpya(write_iter, const_fid, const_fid_len);
	  }
	  *write_iter++ = '\t';
	  write_iter = memcpya(write_iter, sample_line_iter, token_slen);
	  if (sid_present) {
	    write_iter = strcpya(write_iter, "\t0");
	  }
	}
	// PAT/MAT/PHENO1 not required in .psam file
	// SEX now included, so that --vcf + --out has the same effect as --vcf
	// + --make-pgen + --out
	write_iter = memcpyl3a(write_iter, "\tNA");
	append_binary_eoln(&write_iter);
	if (write_iter >= textbuf_flush) {
	  if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), outfile)) {
	    goto vcf_sample_line_ret_WRITE_FAIL;
	  }
	  write_iter = textbuf;
	}
	if (*token_end != '\t') {
	  break;
	}
	sample_line_iter = &(token_end[1]);
      }
      if (write_iter != textbuf) {
	if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), outfile)) {
	  goto vcf_sample_line_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto vcf_sample_line_ret_WRITE_FAIL;
      }
    } else {
      // check consistency of IIDs between VCF and .psam file.
      reterr = gzopen_read_checked(preexisting_psamname, &gz_infile);
      if (reterr) {
	goto vcf_sample_line_ret_1;
      }
      uintptr_t loadbuf_size = bigstack_left();
      if (loadbuf_size > kMaxLongLine) {
	loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size <= kMaxMediumLine) {
	goto vcf_sample_line_ret_NOMEM;
      }
      char* loadbuf = (char*)g_bigstack_base;
      // not formally allocated for now
      loadbuf[loadbuf_size - 1] = ' ';
      char* loadbuf_first_token;
      do {
	++line_idx;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_infile)) {
	    goto vcf_sample_line_ret_READ_FAIL;
	  }
	  loadbuf_first_token = loadbuf;
	  loadbuf_first_token[0] = '\0';
	  break;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  if (loadbuf_size == kMaxLongLine) {
	    goto vcf_sample_line_ret_LONG_LINE;
	  }
	  goto vcf_sample_line_ret_NOMEM;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token) || ((loadbuf_first_token[0] == '#') && strcmp_se(&(loadbuf_first_token[1]), "FID", 3) && strcmp_se(&(loadbuf_first_token[1]), "IID", 3)));
      uint32_t fid_present;
      if (loadbuf_first_token[0] == '#') {
	// only care about position of IID column
	fid_present = (loadbuf_first_token[1] == 'F');
      } else {
	fid_present = fam_cols & kfFamCol1;
      }
      while (1) {
	if (!is_eoln_kns(*loadbuf_first_token)) {
	  char* psam_iid_start = loadbuf_first_token;
	  if (fid_present) {
	    psam_iid_start = skip_initial_spaces(token_endnn(psam_iid_start));
	    if (is_eoln_kns(*psam_iid_start)) {
	      goto vcf_sample_line_ret_MISSING_TOKENS;
	    }
	  }
	  if (((unsigned char)sample_line_iter[0]) < ' ') {
	    sprintf(g_logbuf, "Error: --%ccf file contains fewer sample IDs than %s.\n", flag_char, preexisting_psamname);
	    goto vcf_sample_line_ret_INCONSISTENT_INPUT_WW;
	  }
	  ++sample_ct;
	  char* sample_line_token_end = strchr(sample_line_iter, '\t');
	  if (!sample_line_token_end) {
	    sample_line_token_end = next_prespace(sample_line_iter);
	  }
	  uint32_t sample_line_token_slen = (uintptr_t)(sample_line_token_end - sample_line_iter);
	  if ((*sample_line_iter == '0') && (sample_line_token_slen == 1)) {
	    logerrprint("Error: Sample ID cannot be '0'.\n");
	    goto vcf_sample_line_ret_MALFORMED_INPUT;
	  }
	  char* sample_line_iid_start = sample_line_iter;
	  if (id_delim) {
	    char* first_delim = (char*)memchr(sample_line_iter, (unsigned char)id_delim, sample_line_token_slen);
	    if (!first_delim) {
	      if (!double_or_const_fid) {
		sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
		goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
	      }
	    } else {
	      sample_line_iid_start = &(first_delim[1]);
	      sample_line_token_slen = (uintptr_t)(sample_line_token_end - sample_line_iid_start);
	      char* sample_line_iid_end = (char*)memchr(sample_line_iid_start, (unsigned char)id_delim, sample_line_token_slen);
	      if (sample_line_iid_end) {
		// don't bother erroring out on >2 instances of delimiter for
		// now
		sample_line_token_slen = (uintptr_t)(sample_line_iid_end - sample_line_iid_start);
	      }
	      if ((*sample_line_iid_start == '0') && (sample_line_token_slen == 1)) {
		logerrprint("Error: Sample ID induces an invalid IID of '0'.\n");
		goto vcf_sample_line_ret_INCONSISTENT_INPUT;
	      }
	    }
	  }
	  if (sample_line_token_slen > kMaxIdSlen) {
	    goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  if (memcmp(sample_line_iid_start, psam_iid_start, sample_line_token_slen) || (((unsigned char)psam_iid_start[sample_line_token_slen]) > 32)) {
	    sprintf(g_logbuf, "Error: Mismatched IDs between --%ccf file and %s.\n", flag_char, preexisting_psamname);
	    goto vcf_sample_line_ret_INCONSISTENT_INPUT_WW;
	  }
	  sample_line_iter = &(sample_line_token_end[1]);
	}
	++line_idx;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_infile)) {
	    goto vcf_sample_line_ret_READ_FAIL;
	  }
	  break;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  if (loadbuf_size == kMaxLongLine) {
	    goto vcf_sample_line_ret_LONG_LINE;
	  }
	  goto vcf_sample_line_ret_NOMEM;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
	if (loadbuf_first_token[0] == '#') {
	  sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, preexisting_psamname);
	  goto vcf_sample_line_ret_MALFORMED_INPUT_WW;
	}
      }
      if (gzclose_null(&gz_infile)) {
	goto vcf_sample_line_ret_READ_FAIL;
      }
      if (((unsigned char)sample_line_iter[0]) >= ' ') {
	sprintf(g_logbuf, "Error: --%ccf file contains more sample IDs than %s.\n", flag_char, preexisting_psamname);
	goto vcf_sample_line_ret_INCONSISTENT_INPUT_WW;
      }
    }
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  vcf_sample_line_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  vcf_sample_line_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  vcf_sample_line_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  vcf_sample_line_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  vcf_sample_line_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, preexisting_psamname);
  vcf_sample_line_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  vcf_sample_line_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  vcf_sample_line_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, preexisting_psamname);
    reterr = kPglRetMalformedInput;
    break;
  vcf_sample_line_ret_INCONSISTENT_INPUT_DELIM:
    logerrprintb();
    if (id_delim == '_') {
      logerrprint("If you do not want '_' to be treated as a FID/IID delimiter, use --double-id or\n--const-fid to choose a different method of converting VCF sample IDs to PLINK\nIDs, or --id-delim to change the FID/IID delimiter.\n");
    }
    reterr = kPglRetInconsistentInput;
    break;
  vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID:
    logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
    reterr = kPglRetMalformedInput;
    break;
  vcf_sample_line_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
  vcf_sample_line_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
  vcf_sample_line_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 vcf_sample_line_ret_1:
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  return reterr;
}

uint32_t vcf_is_het_short(const char* first_gchar_ptr, vcf_half_call_t vcf_half_call) {
  // '.' == ascii 46, '0' == ascii 48
  // if kVcfHalfCallReference, ./0 is not phased, but ./1 is
  const uint32_t first_gchar = (unsigned char)first_gchar_ptr[0];
  const uint32_t second_gchar = (unsigned char)first_gchar_ptr[2];
  return (first_gchar != second_gchar) && (((first_gchar != 46) && (second_gchar != 46)) || ((vcf_half_call == kVcfHalfCallReference) && ((first_gchar > 48) || (second_gchar > 48))));
}

uint32_t get_vcf_format_position(const char* __restrict needle, uint32_t needle_slen, char* format_start, char* format_end) {
  *format_end = '\0';
  // assumes first field is GT
  char* token_start = &(format_start[3]);
  uint32_t field_idx = 0;
  while (1) {
    char* token_end = strchr(token_start, ':');
    ++field_idx;
    if (!token_end) {
      if (((uintptr_t)(format_end - token_start) != needle_slen) || memcmp(token_start, needle, needle_slen)) {
	field_idx = 0;
      }
      break;
    }
    if (((uintptr_t)(token_end - token_start) == needle_slen) && (!memcmp(token_start, needle, needle_slen))) {
      break;
    }
    token_start = &(token_end[1]);
  }
  *format_end = '\t';
  return field_idx;
}

uint32_t vcf_qual_scan_init(int32_t vcf_min_gq, int32_t vcf_min_dp, char* format_start, char* format_end, uint32_t* qual_field_skips, int32_t* qual_thresholds) {
  uint32_t gq_field_idx = 0;
  if (vcf_min_gq >= 0) {
    gq_field_idx = get_vcf_format_position("GQ", 2, format_start, format_end);
  }
  uint32_t qual_field_ct = 0;
  uint32_t dp_field_idx = 0;
  if (vcf_min_dp >= 0) {
    dp_field_idx = get_vcf_format_position("DP", 2, format_start, format_end);
    if (dp_field_idx && ((!gq_field_idx) || (dp_field_idx < gq_field_idx))) {
      qual_field_skips[0] = dp_field_idx;
      qual_thresholds[0] = vcf_min_dp;
      qual_field_ct = 1;
      dp_field_idx = 0;
    }
  }
  if (gq_field_idx) {
    qual_field_skips[qual_field_ct] = gq_field_idx;
    qual_thresholds[qual_field_ct++] = vcf_min_gq;
    if (dp_field_idx) {
      qual_field_skips[qual_field_ct] = dp_field_idx;
      qual_thresholds[qual_field_ct++] = vcf_min_dp;
    }
  }
  if (qual_field_ct == 2) {
    qual_field_skips[1] -= qual_field_skips[0];
  }
  return qual_field_ct;
}

char* vcf_field_advance(char* gtext_iter, char* gtext_end, uint32_t advance_ct) {
  // assumes advance_ct is nonzero
  do {
    char* field_end = (char*)memchr(gtext_iter, ':', gtext_end - gtext_iter);
    if (!field_end) {
      return nullptr;
    }
    gtext_iter = &(field_end[1]);
  } while (--advance_ct);
  return gtext_iter;
}

// returns 1 if a quality check failed
// assumes either 1 or 2 qual fields, otherwise change this to a loop
uint32_t vcf_check_quals(const uint32_t* qual_field_skips, const int32_t* qual_thresholds, char* gtext_iter, char* gtext_end, uint32_t qual_field_ct) {
  gtext_iter = vcf_field_advance(gtext_iter, gtext_end, qual_field_skips[0]);
  if (!gtext_iter) {
    return 0;
  }
  int32_t ii;
  if ((!scan_int32(gtext_iter, &ii)) && (ii < qual_thresholds[0])) {
    return 1;
  }
  if (qual_field_ct == 1) {
    return 0;
  }
  gtext_iter = vcf_field_advance(gtext_iter, gtext_end, qual_field_skips[1]);
  if (!gtext_iter) {
    return 0;
  }
  return (!scan_int32(gtext_iter, &ii)) && (ii < qual_thresholds[1]);
}

boolerr_t parse_vcf_gp(char* gp_iter, uint32_t is_haploid, double import_dosage_certainty, uint32_t* is_missing_ptr, double* alt_dosage_ptr) {
  // P(0/0), P(0/1), P(1/1), etc.
  // assumes is_missing initialized to 0
  double prob_0alt;
  gp_iter = scanadv_double(gp_iter, &prob_0alt);
  if ((!gp_iter) || (prob_0alt < 0.0) || (prob_0alt > 1.0) || (*gp_iter != ',')) {
    return 1;
  }
  double prob_1alt;
  gp_iter = scanadv_double(&(gp_iter[1]), &prob_1alt);
  if ((!gp_iter) || (prob_1alt < 0.0) || (prob_1alt > 1.0)) {
    return 1;
  }
  if (is_haploid) {
    const double denom = prob_0alt + prob_1alt;
    if (denom <= 2 * import_dosage_certainty) {
      if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty)) {
	*is_missing_ptr = 1;
	return 1;
      }
    }
    *alt_dosage_ptr = 2 * prob_1alt / denom;
    return 0;
  }
  double prob_2alt;
  if ((*gp_iter != ',') || (!scanadv_double(&(gp_iter[1]), &prob_2alt)) || (prob_2alt < 0.0) || (prob_2alt > 1.0)) {
    return 1;
  }
  const double denom = prob_0alt + prob_1alt + prob_2alt;
  if (denom <= 3 * import_dosage_certainty) {
    if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty) && (prob_2alt <= import_dosage_certainty)) {
      // treat as missing
      // ok to use <= since we multiplied by (1 - epsilon)
      // during command-line parsing.  this lets us avoid
      // special-casing denom=0.
      *is_missing_ptr = 1;
      return 1; // not really an error
    }
  }
  *alt_dosage_ptr = (prob_1alt + 2 * prob_2alt) / denom;
  return 0;
}

boolerr_t parse_vcf_dosage(char* gtext_iter, char* gtext_end, uint32_t dosage_field_idx, uint32_t is_haploid, uint32_t dosage_is_gp, double import_dosage_certainty, uint32_t* is_missing_ptr, uint32_t* dosage_int_ptr) {
  // assumes is_missing initialized to 0
  // assumes dosage_field_idx != 0
  // returns 1 if missing OR parsing error.
  uint32_t field_idx = 0;
  do {
    char* field_end = (char*)memchr(gtext_iter, ':', gtext_end - gtext_iter);
    if (!field_end) {
      *is_missing_ptr = 1;
      return 1;
    }
    gtext_iter = &(field_end[1]);
  } while (++field_idx < dosage_field_idx);
  if (((gtext_iter[0] == '.') || (gtext_iter[0] == '?')) && (((uint32_t)((unsigned char)gtext_iter[1])) - 48 >= 10)) {
    // missing field (dot/'?' followed by non-digit)
    // could enforce gtext_iter[1] == colon, comma, etc.?
    *is_missing_ptr = 1;
    return 1;
  }
  double alt_dosage;
  if (dosage_is_gp) {
    if (parse_vcf_gp(gtext_iter, is_haploid, import_dosage_certainty, is_missing_ptr, &alt_dosage)) {
      return 1;
    }
  } else {
    if ((!scanadv_double(gtext_iter, &alt_dosage)) || (alt_dosage < 0.0)) {
      return 1;
    }
    if (is_haploid) {
      // possible todo: allow this to be suppressed (maybe upstream of this
      // function); 1000 Genomes phase 1 haploid dosages are still on 0..2
      // scale
      alt_dosage *= 2;
    }
    if (alt_dosage > 2.0) {
      return 1;
    }
  }
  *dosage_int_ptr = (int32_t)(alt_dosage * kDosageMid + 0.5);
  return 0;
}

// dosage_int = 0..2 value in 16384ths
// returns distance from 0.5 or 1.5 in 16384ths, whichever is closer
static inline uint32_t biallelic_dosage_halfdist(uint32_t dosage_int) {
  const uint32_t dosage_int_rem = dosage_int & (kDosageMid - 1);
  return abs_int32(((int32_t)dosage_int_rem) - kDosage4th);
}

static_assert(!kVcfHalfCallReference, "vcf_to_pgen() assumes kVcfHalfCallReference == 0.");
static_assert(kVcfHalfCallHaploid == 1, "vcf_to_pgen() assumes kVcfHalfCallHaploid == 1.");
pglerr_t vcf_to_pgen(const char* vcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, misc_flags_t misc_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, vcf_half_call_t vcf_half_call, fam_col_t fam_cols, char* outname, char* outname_end, chr_info_t* cip) {
  // Now performs a 2-pass load.  Yes, this can be slower than plink 1.9, but
  // it's necessary to use the Pgen_writer classes for now (since we need to
  // know upfront how many variants there are, and whether phase/dosage is
  // present).
  // preexisting_psamname should be nullptr if no such file was specified.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* pvarfile = nullptr;
  uintptr_t line_idx = 0;
  const uint32_t vcf_half_call_explicit_error = (vcf_half_call == kVcfHalfCallError);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  gzFile gz_infile;
  spgw_preinit(&spgw);
  {
    // don't use gzopen_read_checked() since we want to customize the error
    // message
    gz_infile = gzopen(vcfname, FOPEN_RB);
    if (!gz_infile) {
      const uint32_t slen = strlen(vcfname);
      if (((slen > 4) && (!memcmp(&(vcfname[slen - 4]), ".vcf", 4))) || ((slen > 7) && (!memcmp(&(vcfname[slen - 7]), ".vcf.gz", 7)))) {
	LOGERRPRINTFWW(g_errstr_fopen, vcfname);
      } else {
	LOGERRPRINTFWW("Error: Failed to open %s. (--vcf expects a complete filename; did you forget '.vcf' at the end?)\n", vcfname);
      }
      goto vcf_to_pgen_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    uintptr_t loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto vcf_to_pgen_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kEndAllocAlign);
    }
    char* loadbuf = (char*)bigstack_end_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_iter;
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t dosage_import_field_slen = 0;
    if (dosage_import_field) {
      dosage_import_field_slen = strlen(dosage_import_field);
    }
    const uint32_t dosage_is_gp = (dosage_import_field_slen == 2) && (!memcmp(dosage_import_field, "GP", 2));
    uint32_t format_gt_present = 0;
    uint32_t format_gq_relevant = 0;
    uint32_t format_dp_relevant = 0;
    uint32_t format_dosage_relevant = 0;
    uint32_t info_pr_present = 0;
    uint32_t info_nonpr_present = 0;
    uint32_t chrset_present = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto vcf_to_pgen_ret_READ_FAIL;
	}
	logerrprint("Error: No #CHROM header line or variant records in --vcf file.\n");
	goto vcf_to_pgen_ret_MALFORMED_INPUT;
      }
      if ((line_idx == 1) && (!memcmp(loadbuf, "BCF", 3))) {
	// this is more informative than "missing header line"...
	if (loadbuf[3] == 2) {
	  sprintf(g_logbuf, "Error: %s appears to be a BCF2 file. Try --bcf instead of --vcf.\n", vcfname);
	  goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	if (loadbuf[3] == 4) {
	  sprintf(g_logbuf, "Error: %s appears to be a BCF1 file. Use 'bcftools view' to convert it to a PLINK-readable VCF.\n", vcfname);
	  goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
	}
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto vcf_to_pgen_ret_LONG_LINE;
	}
	goto vcf_to_pgen_ret_NOMEM;
      }
      // don't tolerate leading spaces
      loadbuf_iter = loadbuf;
      if (*loadbuf_iter != '#') {
	logerrprint("Error: No #CHROM header line in --vcf file.\n");
	goto vcf_to_pgen_ret_MALFORMED_INPUT;
      }
      if (loadbuf_iter[1] != '#') {
	break;
      }
      // Recognized header lines:
      // ##fileformat: discard (regenerate; todo: conditionally error out)
      // ##fileDate: discard (regenerate)
      // ##source: discard (regenerate)
      // ##contig: conditionally keep
      // ##INFO: note presence of INFO:PR, note presence of at least one non-PR
      //         field, keep data (though, if INFO:PR is the *only* field,
      //         omit it from the .pvar for consistency with --make-pgen
      //         default)
      // ##FORMAT: note presence of FORMAT:GT and FORMAT:GP, discard
      //           (regenerate)
      // ##chrSet: if recognized, perform consistency check and/or update
      //           chr_info
      //
      // Everything else (##FILTER, ##reference, etc.) is passed through
      // unchanged.  FILTER values in the VCF body do not have to be mentioned
      // in the header (since only BCF, not VCF, spec requires that).
      //
      // Because of how ##contig is handled (we only keep the lines which
      // correspond to chromosomes/contigs actually present in the VCF, and not
      // filtered out), we wait until second pass to write the .pvar.
      if (!memcmp(&(loadbuf_iter[2]), "chrSet=<", 8)) {
	if (chrset_present) {
	  logerrprint("Error: Multiple ##chrSet header lines in --vcf file.\n");
	  goto vcf_to_pgen_ret_MALFORMED_INPUT;
	}
	chrset_present = 1;
	// .pvar loader will print a warning if necessary
	reterr = read_chrset_header_line(&(loadbuf_iter[10]), "--vcf file", misc_flags, line_idx, cip);
	if (reterr) {
	  goto vcf_to_pgen_ret_1;
	}
      } else if (!memcmp(&(loadbuf_iter[2]), "FORMAT=<ID=GT,Number=", 21)) {
	if (format_gt_present) {
	  logerrprint("Error: Duplicate FORMAT:GT header line in --vcf file.\n");
	  goto vcf_to_pgen_ret_MALFORMED_INPUT;
	}
	if (memcmp(&(loadbuf_iter[23]), "1,Type=String,Description=", 26)) {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected FORMAT:GT format.\n", line_idx);
	  goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	format_gt_present = 1;
      } else if ((vcf_min_gq != -1) && (!memcmp(&(loadbuf_iter[2]), "FORMAT=<ID=GQ,Number=1,Type=", 21))) {
	if (format_gq_relevant) {
	  logerrprint("Error: Duplicate FORMAT:GQ header line in --vcf file.\n");
	  goto vcf_to_pgen_ret_MALFORMED_INPUT;
	}
	format_gq_relevant = 1;
      } else if ((vcf_min_dp != -1) && (!memcmp(&(loadbuf_iter[2]), "FORMAT=<ID=DP,Number=1,Type=", 21))) {
	if (format_dp_relevant) {
	  logerrprint("Error: Duplicate FORMAT:DP header line in --vcf file.\n");
	  goto vcf_to_pgen_ret_MALFORMED_INPUT;
	}
	format_dp_relevant = 1;
      } else if (dosage_import_field && (!memcmp(&(loadbuf_iter[2]), "FORMAT=<ID=", 11)) && (!memcmp(&(loadbuf_iter[13]), dosage_import_field, dosage_import_field_slen)) && (loadbuf_iter[13 + dosage_import_field_slen] == ',')) {
	if (format_dosage_relevant) {
	  LOGERRPRINTFWW("Error: Duplicate FORMAT:%s header line in --vcf file.\n", dosage_import_field);
	  goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	format_dosage_relevant = 1;
      } else if (!memcmp(&(loadbuf_iter[2]), "INFO=<ID=", 9)) {
	if (!memcmp(&(loadbuf_iter[11]), "PR,Number=", 10)) {
	  if (info_pr_present) {
	    logerrprint("Error: Duplicate INFO:PR header line in --vcf file.\n");
	    goto vcf_to_pgen_ret_MALFORMED_INPUT;
	  }
	  if (memcmp(&(loadbuf_iter[21]), "0,Type=Flag,Description=", 24)) {
	    sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected INFO:PR format.\n", line_idx);
	    goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
	  }
	  info_pr_present = 1;
	} else {
	  info_nonpr_present = 1;
	}
      }
    }
    const uint32_t require_gt = (misc_flags / kfMiscVcfRequireGt) & 1;
    if ((!format_gt_present) && require_gt) {
      // todo: allow_no_variants exception
      logerrprint("Error: No GT field in --vcf file header, when --vcf-require-gt was specified.\n");
      goto vcf_to_pgen_ret_INCONSISTENT_INPUT;
    }
    if ((!format_gq_relevant) && (vcf_min_gq != -1)) {
      logerrprint("Warning: No GQ field in --vcf file header.  --vcf-min-gq ignored.\n");
      vcf_min_gq = -1;
    }
    if ((!format_dp_relevant) && (vcf_min_dp != -1)) {
      logerrprint("Warning: No DP field in --vcf file header.  --vcf-min-dp ignored.\n");
      vcf_min_dp = -1;
    }
    const uint32_t format_gq_or_dp_relevant = format_gq_relevant || format_dp_relevant;
    if ((!format_dosage_relevant) && dosage_import_field) {
      LOGERRPRINTFWW("Warning: No %s field in --vcf file header. Dosages will not be imported.\n", dosage_import_field);
    }
    finalize_chrset(misc_flags, cip);
    // don't call finalize_chr_info here, since this may be followed by
    // --pmerge, etc.
    
    if (memcmp(loadbuf_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", 38)) {
      sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected field sequence after #CHROM.\n", line_idx);
      goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
    }
    loadbuf_iter = &(loadbuf_iter[38]);
    const uint32_t double_id = (misc_flags / kfMiscDoubleId) & 1;
    uintptr_t sample_ct = 0;
    if (!memcmp(loadbuf_iter, "\tFORMAT\t", 8)) {
      reterr = vcf_sample_line(preexisting_psamname, const_fid, double_id, fam_cols, id_delim, idspace_to, 'v', &(loadbuf_iter[8]), outname, outname_end, &sample_ct);
      if (reterr) {
	goto vcf_to_pgen_ret_1;
      }
      /*
    } else if (allow_no_samples) {
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, "w", &outfile)) {
	goto vcf_to_pgen_ret_OPEN_FAIL;
      }
      fputs("#FID\tIID\n", outfile);
      if (fclose_null(&outfile)) {
	goto vcf_to_pgen_ret_WRITE_FAIL;
      }
      */
    }
    // todo: allow_no_samples exception
    if (!sample_ct) {
      logerrprint("Error: No samples in --vcf file.\n");
      goto vcf_to_pgen_ret_INCONSISTENT_INPUT;
    }

    uint32_t variant_ct = 0;
    uint32_t max_alt_ct = 1;
    uintptr_t* variant_allele_idxs = (uintptr_t*)g_bigstack_base;
    uintptr_t max_variant_ct = (uintptr_t)(((uintptr_t*)g_bigstack_end) - variant_allele_idxs);
    max_variant_ct -= BITCT_TO_ALIGNED_WORDCT(max_variant_ct) * kWordsPerVec;
    if (format_dosage_relevant) {
      max_variant_ct -= BITCT_TO_ALIGNED_WORDCT(max_variant_ct) * kWordsPerVec;
    }
    if (info_pr_present) {
      max_variant_ct -= BITCT_TO_ALIGNED_WORDCT(max_variant_ct) * kWordsPerVec;
    }
#ifdef __LP64__
    if (max_variant_ct > 0x7ffffffd) {
      max_variant_ct = 0x7ffffffd;
    }
#endif
    uintptr_t base_chr_present[kChrExcludeWords];
    fill_ulong_zero(kChrExcludeWords, base_chr_present);
    
    const uintptr_t header_line_ct = line_idx;
    const uint32_t max_variant_ctaw = BITCT_TO_ALIGNED_WORDCT(max_variant_ct);
    uintptr_t* phasing_flags = (uintptr_t*)bigstack_end_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t));
    uintptr_t* phasing_flags_iter = phasing_flags;
    uintptr_t* dosage_flags = nullptr;
    if (format_dosage_relevant) {
      dosage_flags = (uintptr_t*)bigstack_end_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t));
    }
    uintptr_t* dosage_flags_iter = dosage_flags;
    uintptr_t* nonref_flags = nullptr;
    if (info_pr_present) {
      nonref_flags = (uintptr_t*)bigstack_end_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t));
    }
    uintptr_t* nonref_flags_iter = nonref_flags;
    if (vcf_half_call == kVcfHalfCallDefault) {
      vcf_half_call = kVcfHalfCallError;
    }
    uintptr_t variant_skip_ct = 0;
    uintptr_t phasing_word = 0;
    uintptr_t dosage_word = 0;
    uintptr_t nonref_word = 0;
    uintptr_t allele_idx_end = 0;
    uint32_t max_allele_slen = 1;
    uint32_t max_qualfilterinfo_slen = 6;
    uint32_t qual_field_ct = 0;

    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

    // temporary kludge
    uintptr_t multiallelic_skip_ct = 0;

    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto vcf_to_pgen_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto vcf_to_pgen_ret_LONG_LINE_N;
	}
	goto vcf_to_pgen_ret_NOMEM;
      }
      // do tolerate trailing newlines
      if ((unsigned char)(*loadbuf) <= 32) {
	if (*loadbuf == ' ') {
	  sprintf(g_logbuf, "Error: Leading space on line %" PRIuPTR " of --vcf file.\n", line_idx);
	  goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
	}
	continue;
      }
      loadbuf_iter = loadbuf;
      char* chr_code_end = strchr(loadbuf, '\t');
      if (!chr_code_end) {
	goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      // QUAL/FILTER enforcement is now postponed till .pvar loading.  only
      // other things we do during the scanning pass are (i) count alt alleles,
      // and (ii) check whether any phased genotype calls are present.

      char* pos_end = strchr(&(chr_code_end[1]), '\t');
      if (!pos_end) {
	goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      
      // may as well check ID length here
      // postpone POS validation till second pass so we only have to parse it
      // once
      char* id_end = strchr(&(pos_end[1]), '\t');
      if (!id_end) {
	goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      if ((uintptr_t)(id_end - pos_end) > kMaxIdBlen) {
	sprintf(g_logbuf, "Error: Invalid ID on line %" PRIuPTR " of --vcf file (max " MAX_ID_SLEN_STR " chars).\n", line_idx);
	goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
      }
      
      // note REF length
      char* ref_allele_start = &(id_end[1]);
      loadbuf_iter = strchr(ref_allele_start, '\t');
      if (!loadbuf_iter) {
	goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      uint32_t cur_max_allele_slen = (uintptr_t)(loadbuf_iter - ref_allele_start);

      uint32_t alt_ct = 1;
      unsigned char ucc;
      // treat ALT=. as if it were an actual allele for now
      while (1) {
	char* cur_allele_start = ++loadbuf_iter;
	ucc = (unsigned char)(*loadbuf_iter);
	if ((ucc <= ',') && (ucc != '*')) {
	  sprintf(g_logbuf, "Error: Invalid alternate allele on line %" PRIuPTR " of --vcf file.\n", line_idx);
	  goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
	}
	do {
	  ucc = (unsigned char)(*(++loadbuf_iter));
	  // allow GATK 3.4 <*:DEL> symbolic allele
	} while ((ucc > ',') || (ucc == '*'));
	const uint32_t cur_allele_slen = (uintptr_t)(loadbuf_iter - cur_allele_start);
	if (cur_allele_slen > cur_max_allele_slen) {
	  cur_max_allele_slen = cur_allele_slen;
	}
	if (ucc != ',') {
	  break;
	}
	++alt_ct;
      }

      // temporary kludge
      if (alt_ct > 1) {
	++multiallelic_skip_ct;
	continue;
      }
      
      if (ucc != '\t') {
	sprintf(g_logbuf, "Error: Malformed ALT field on line %" PRIuPTR " of --vcf file.\n", line_idx);
	goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
      }
      if (alt_ct > max_alt_ct) {
	max_alt_ct = alt_ct;
      }

      // skip QUAL, FILTER
      char* qual_start_m1 = loadbuf_iter;
      for (uint32_t uii = 0; uii < 2; ++uii) {
	loadbuf_iter = strchr(&(loadbuf_iter[1]), '\t');
	if (!loadbuf_iter) {
	  goto vcf_to_pgen_ret_MISSING_TOKENS;
	}
      }
      
      // possibly check for FORMAT:GT before proceeding
      char* info_start = &(loadbuf_iter[1]);
      char* info_end = strchr(info_start, '\t');
      if (!info_end) {
	goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      loadbuf_iter = &(info_end[1]);
      const uint32_t gt_missing = memcmp(loadbuf_iter, "GT", 2) || ((loadbuf_iter[2] != ':') && (loadbuf_iter[2] != '\t'));
      if (require_gt && gt_missing) {
	++variant_skip_ct;
	continue;
      }
      const uint32_t cur_qualfilterinfo_slen = (uintptr_t)(info_end - qual_start_m1);

      // all converters *do* respect chromosome filters
      // wait till this point to apply it, since we don't want to
      // add a contig name to the hash table unless at least one variant on
      // that contig wasn't filtered out for other reasons.
      int32_t cur_chr_code;
      reterr = get_or_add_chr_code_destructive("--vcf file", line_idx, allow_extra_chrs, loadbuf, chr_code_end, cip, &cur_chr_code);
      if (reterr) {
	goto vcf_to_pgen_ret_1;
      }
      if (!is_set(cip->chr_mask, cur_chr_code)) {
	++variant_skip_ct;
	continue;
      }
      if (cur_max_allele_slen > max_allele_slen) {
	max_allele_slen = cur_max_allele_slen;
      }
      if (cur_qualfilterinfo_slen > max_qualfilterinfo_slen) {
	max_qualfilterinfo_slen = cur_qualfilterinfo_slen;
      }
      if ((uint32_t)cur_chr_code <= cip->max_code) {
	set_bit(cur_chr_code, base_chr_present);
      }

      variant_allele_idxs[variant_ct] = allele_idx_end;
      allele_idx_end += alt_ct + 1;
      const uint32_t variant_idx_lowbits = variant_ct % kBitsPerWord;
      if (info_pr_present) {
	if (pr_in_info_token((uintptr_t)(info_end - info_start), info_start)) {
	  nonref_word |= k1LU << variant_idx_lowbits;
	}
	if (variant_idx_lowbits == (kBitsPerWord - 1)) {
	  *nonref_flags_iter++ = nonref_word;
	  nonref_word = 0;
	}
      }
      if (!gt_missing) {
	// todo: sample_ct == 0 case
	// possible todo: import dosages when GT missing
	
	// loadbuf_iter currently points to beginning of FORMAT field
	char* format_end = strchr(loadbuf_iter, '\t');
	if (!format_end) {
	  goto vcf_to_pgen_ret_MISSING_TOKENS;
	}

	uint32_t qual_field_skips[2];
	int32_t qual_thresholds[2];
	if (format_gq_or_dp_relevant) {
	  qual_field_ct = vcf_qual_scan_init(vcf_min_gq, vcf_min_dp, loadbuf_iter, format_end, qual_field_skips, qual_thresholds);
	}
        // if nonzero, 0-based index of dosage field
	uint32_t dosage_field_idx = 0;
	if (format_dosage_relevant) {
	  dosage_field_idx = get_vcf_format_position(dosage_import_field, dosage_import_field_slen, loadbuf_iter, format_end);
	}

	// check if there's at least one phased het call, and/or at least one
	// relevant dosage
	if (alt_ct < 10) {
	  // always check for a phased het
	  char* phasescan_iter = format_end;
	  while (1) {
	    // this should quickly fail if there are no phased calls at all.
	    phasescan_iter = strchr(&(phasescan_iter[1]), '|');
	    if (!phasescan_iter) {
	      break;
	    }
	    if (phasescan_iter[-2] != '\t') {
	      // at least one other gdata field uses the '|' character.
	      // switch to iterating over tabs.
	      while (1) {
		phasescan_iter = strchr(&(phasescan_iter[1]), '\t');
		if (!phasescan_iter) {
		  break;
		}
		if (phasescan_iter[2] == '|') {
		  if (vcf_is_het_short(&(phasescan_iter[1]), vcf_half_call)) {
		    if (qual_field_ct) {
		      char* cur_gtext_end = strchr(phasescan_iter, '\t');
		      if (!cur_gtext_end) {
			cur_gtext_end = next_prespace(phasescan_iter);
		      }
		      if (vcf_check_quals(qual_field_skips, qual_thresholds, phasescan_iter, cur_gtext_end, qual_field_ct)) {
			break;
		      }
		    }
		    phasing_word |= k1LU << variant_idx_lowbits;
		    break;
		  }
		}
	      }
	      break;
	    }
	    if (vcf_is_het_short(&(phasescan_iter[-1]), vcf_half_call)) {
	      if (qual_field_ct) {
		char* cur_gtext_end = strchr(phasescan_iter, '\t');
		if (!cur_gtext_end) {
		  cur_gtext_end = next_prespace(phasescan_iter);
		}
		if (vcf_check_quals(qual_field_skips, qual_thresholds, phasescan_iter, cur_gtext_end, qual_field_ct)) {
		  break;
		}
	      }
	      phasing_word |= k1LU << variant_idx_lowbits;
	      break;
	    }
	    phasescan_iter = strchr(&(phasescan_iter[1]), '\t');
	    if (!phasescan_iter) {
	      break;
	    }
	  }
	  if (dosage_field_idx) {
	    char* dosagescan_iter = format_end;
	    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	      char* cur_gtext_start = ++dosagescan_iter;
	      char* cur_gtext_end = strchr(dosagescan_iter, '\t');
	      if (!cur_gtext_end) {
		if (sample_idx + 1 != sample_ct) {
		  goto vcf_to_pgen_ret_MISSING_TOKENS;
		}
		cur_gtext_end = next_prespace(dosagescan_iter);
	      }
	      dosagescan_iter = cur_gtext_end;
	      if (qual_field_ct) {
		if (vcf_check_quals(qual_field_skips, qual_thresholds, cur_gtext_start, cur_gtext_end, qual_field_ct)) {
		  break;
		}
	      }
	      const uint32_t is_haploid = (cur_gtext_start[1] != '/') && (cur_gtext_start[1] != '|');
	      uint32_t is_missing = 0;
	      uint32_t dosage_int;
	      if (parse_vcf_dosage(cur_gtext_start, cur_gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &is_missing, &dosage_int)) {
		if (is_missing) {
		  continue;
		}
		goto vcf_to_pgen_ret_INVALID_DOSAGE;
	      }
	      const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
	      if (cur_halfdist < dosage_erase_halfdist) {
		goto vcf_to_pgen_dosagescan_hit;
	      }
	    }
	  }
	  if (0) {
	  vcf_to_pgen_dosagescan_hit:
	    dosage_word |= k1LU << variant_idx_lowbits;
	  }
	} else {
	  // alt_ct >= 10
	  // todo
	}
      }
      if (variant_idx_lowbits == (kBitsPerWord - 1)) {
	*phasing_flags_iter++ = phasing_word;
	phasing_word = 0;
	if (dosage_flags_iter) {
	  *dosage_flags_iter++ = dosage_word;
	  dosage_word = 0;
	}
      }
      if (variant_ct++ == max_variant_ct) {
#ifdef __LP64__
	if (variant_ct == 0x7ffffffd) {
	  logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend other\nsoftware, such as PLINK/SEQ, for very deep studies of small numbers of genomes.\n");
	  goto vcf_to_pgen_ret_MALFORMED_INPUT;
	}
#endif
	goto vcf_to_pgen_ret_NOMEM;
      }
      if (!(variant_ct % 1000)) {
	printf("\r--vcf: %uk variants scanned.", variant_ct / 1000);
	fflush(stdout);
      }
    }
    if (variant_ct % kBitsPerWord) {
      *phasing_flags_iter = phasing_word;
      if (dosage_flags_iter) {
	*dosage_flags_iter = dosage_word;
      }
      if (nonref_flags_iter) {
	*nonref_flags_iter = nonref_word;
      }
    } else if (!variant_ct) {
      // todo: allow_no_variants exception
      logerrprint("Error: No variants in --vcf file.\n");
      goto vcf_to_pgen_ret_INCONSISTENT_INPUT;
    }
    
    putc_unlocked('\r', stdout);
    if (!variant_skip_ct) {
      LOGPRINTF("--vcf: %u variant%s scanned.\n", variant_ct, (variant_ct == 1)? "" : "s");
    } else {
      LOGPRINTF("--vcf: %u variant%s scanned (%" PRIuPTR " skipped).\n", variant_ct, (variant_ct == 1)? "" : "s", variant_skip_ct);
    }

    // temporary kludge
    if (multiallelic_skip_ct) {
      LOGERRPRINTFWW("Warning: %" PRIuPTR " multiallelic variant%s %sskipped (not yet supported).\n", multiallelic_skip_ct, (multiallelic_skip_ct == 1)? "" : "s", variant_skip_ct? "also " : "");
    }
    
    if (gzrewind(gz_infile)) {
      goto vcf_to_pgen_ret_READ_FAIL;
    }
    const uintptr_t line_ct = line_idx - 1;

    if (allele_idx_end > 2 * variant_ct) {
      variant_allele_idxs[variant_ct] = allele_idx_end;
      bigstack_finalize_ul(variant_allele_idxs, variant_ct + 1);
    } else {
      variant_allele_idxs = nullptr;
    }
    
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &pvarfile)) {
      goto vcf_to_pgen_ret_OPEN_FAIL;
    }
    line_idx = 0;
    while (1) {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	goto vcf_to_pgen_ret_READ_FAIL;
      }
      if (++line_idx == header_line_ct) {
	break;
      }
      // don't use textbuf here, since header line length could theoretically
      // exceed kMaxMediumLine bytes
      if ((!memcmp(loadbuf, "##fileformat=", 13)) || (!memcmp(loadbuf, "##fileDate=", 11)) || (!memcmp(loadbuf, "##source=", 9)) || (!memcmp(loadbuf, "##FORMAT=", 9)) || (!memcmp(loadbuf, "##chrSet=", 9))) {
	continue;
      }
      const uint32_t line_slen = strlen(loadbuf);
      if (!memcmp(loadbuf, "##contig=<ID=", 13)) {
	char* contig_name_start = &(loadbuf[13]);
	char* contig_name_end = strchr(contig_name_start, ',');
	if (!contig_name_end) {
	  // could search backwards from end of line in this case, but if
	  // contig names are long enough for that to matter we have other
	  // problems...
	  contig_name_end = strchr(contig_name_start, '>');
	  if (!contig_name_end) {
	    sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected ##contig format.\n", line_idx);
	    goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
	  }
	}
        const int32_t cur_chr_code = get_chr_code_counted(cip, (uintptr_t)(contig_name_end - contig_name_start), contig_name_start);
	if (cur_chr_code < 0) {
	  continue;
	}
	if ((uint32_t)cur_chr_code <= cip->max_code) {
	  if (!is_set(base_chr_present, cur_chr_code)) {
	    continue;
	  }
	} else {
	  if (!is_set(cip->chr_mask, cur_chr_code)) {
	    continue;
	  }
	}
	// Note that, when --output-chr is specified, we don't update the
	// ##contig header line chromosome code in the .pvar file, since
	// ##contig is not an explicit part of the .pvar specification, it's
	// just another blob of text as far as the main body of plink2 is
	// concerned.  However, the codes are brought in sync during VCF/BCF
	// export.
      }
      // force OS-appropriate eoln
      // don't need to check for \n since this is pass #2, we already validated
      // file contents, and we can't be on the last line of a file lacking a
      // final eoln since #CHROM is still to come
      char* line_end = &(loadbuf[line_slen - 1]);
      if (line_end[-1] == '\r') {
	--line_end;
      }
      append_binary_eoln(&line_end);
      if (fwrite_checked(loadbuf, (uintptr_t)(line_end - loadbuf), pvarfile)) {
	goto vcf_to_pgen_ret_WRITE_FAIL;
      }
    }
    char* write_iter = g_textbuf;
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    write_iter = strcpya(write_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    if (info_nonpr_present) {
      write_iter = strcpya(write_iter, "\tINFO");
    }
    append_binary_eoln(&write_iter);
    if (fwrite_checked(g_textbuf, write_iter - g_textbuf, pvarfile)) {
      goto vcf_to_pgen_ret_WRITE_FAIL;
    }
    
    const uint32_t variant_ctl = BITCT_TO_WORDCT(variant_ct);
    pgen_global_flags_t phase_dosage_gflags = are_all_words_zero(phasing_flags, variant_ctl)? kfPgenGlobal0 : kfPgenGlobalHardcallPhasePresent;
    if (format_dosage_relevant) {
      if (are_all_words_zero(dosage_flags, variant_ctl)) {
	format_dosage_relevant = 0;
      } else {
        phase_dosage_gflags |= kfPgenGlobalDosagePresent;
      }
    }
    uint32_t nonref_flags_storage = 1;
    if (nonref_flags) {
      const uint32_t variant_ctl_m1 = variant_ctl - 1;
      const uintptr_t last_nonref_flags_word = nonref_flags[variant_ctl_m1];
      if (!last_nonref_flags_word) {
	for (uint32_t widx = 0; widx < variant_ctl_m1; ++widx) {
	  if (nonref_flags[widx]) {
	    nonref_flags_storage = 3;
	    break;
	  }
	}
	// todo: replace kBitsPerWord - MOD_NZ...
      } else if (!((~last_nonref_flags_word) << (kBitsPerWord - MOD_NZ(variant_ct, kBitsPerWord)))) {
	nonref_flags_storage = 2;
	for (uint32_t widx = 0; widx < variant_ctl_m1; ++widx) {
	  if (~nonref_flags[widx]) {
	    nonref_flags_storage = 3;
	    break;
	  }
	}
      } else {
	nonref_flags_storage = 3;
      }
      if (nonref_flags_storage != 3) {
	nonref_flags = nullptr;
	bigstack_end_reset(phasing_flags);
      }
    }
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, variant_allele_idxs, nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto vcf_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (bigstack_alloc_ul(sample_ctl2, &genovec)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    // nothing should go wrong if trailing word is garbage, but keep an eye on
    // this
    // fill_ulong_zero(sample_ctaw2 - sample_ctl2, &(genovec[sample_ctl2]));
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* phasepresent = nullptr;
    uintptr_t* phaseinfo = nullptr;
    if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
      if (bigstack_alloc_ul(sample_ctl, &phasepresent) ||
	  bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
	goto vcf_to_pgen_ret_NOMEM;
      }
    }
    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (phase_dosage_gflags & kfPgenGlobalDosagePresent) {
      if (bigstack_alloc_ul(sample_ctl, &dosage_present) ||
	  bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
	goto vcf_to_pgen_ret_NOMEM;
      }
      dosage_present[sample_ctl - 1] = 0;
    }

    char* writebuf;
    if (bigstack_alloc_c(2 * max_allele_slen + max_qualfilterinfo_slen + kCompressStreamBlock + kMaxIdSlen + 32, &writebuf)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kCompressStreamBlock]);

    if (hard_call_thresh == 0xffffffffU) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;

    uint32_t vidx = 0;
    for (line_idx = header_line_ct + 1; line_idx <= line_ct; ++line_idx) {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	goto vcf_to_pgen_ret_READ_FAIL;
      }
      if ((unsigned char)(*loadbuf) < 32) {
	continue;
      }
      // 1. check if we skip this variant.  chromosome filter, require_gt, and
      //    (temporarily) multiple alt alleles can cause this.
      char* chr_code_end = (char*)rawmemchr(loadbuf, '\t');
      int32_t chr_code_base = get_chr_code_raw(loadbuf);
      if (chr_code_base == -1) {
	// skip hash table lookup if we know we aren't skipping the variant
	if (variant_skip_ct) {
	  *chr_code_end = '\0';
	  const uint32_t chr_code = id_htable_find(loadbuf, cip->nonstd_names, cip->nonstd_id_htable, (uintptr_t)(chr_code_end - loadbuf), kChrHtableSize);
	  if ((chr_code == 0xffffffffU) || (!IS_SET(cip->chr_mask, chr_code))) {
	    continue;
	  }
	  *chr_code_end = '\t';
	}
      } else {
	if (chr_code_base >= ((int32_t)kMaxContigs)) {
	  chr_code_base = cip->xymt_codes[chr_code_base - kMaxContigs];
	}
	if ((chr_code_base < 0) || (!is_set(base_chr_present, chr_code_base))) {
	  assert(variant_skip_ct);
	  continue;
	}
      }
      // chr_code_base is now a proper numeric chromosome index for
      // non-contigs, and -1 if it's a contig name
      char* pos_str = &(chr_code_end[1]);
      char* pos_str_end = (char*)rawmemchr(pos_str, '\t');
      loadbuf_iter = pos_str_end;
      // copy ID, REF verbatim
      for (uint32_t uii = 0; uii < 2; ++uii) {
	loadbuf_iter = (char*)rawmemchr(&(loadbuf_iter[1]), '\t');
      }

      // ALT, QUAL, FILTER, INFO
      char* filter_end = loadbuf_iter;
      for (uint32_t uii = 0; uii < 3; ++uii) {
	filter_end = (char*)rawmemchr(&(filter_end[1]), '\t');
      }
      char* info_end = (char*)rawmemchr(&(filter_end[1]), '\t');
      char* format_start = &(info_end[1]);
      const uint32_t gt_missing = memcmp(format_start, "GT", 2) || ((format_start[2] != ':') && (format_start[2] != '\t'));
      if (require_gt && gt_missing) {
	continue;
      }
      
      // make sure POS starts with an integer, apply --output-chr setting
      uint32_t cur_bp;
      if (scan_uint_defcap(pos_str, &cur_bp)) {
	sprintf(g_logbuf, "Error: Invalid POS on line %" PRIuPTR " of --vcf file.\n", line_idx);
	goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
      }

      // temporary kludge
      char* write_line_start = write_iter;

      if (chr_code_base == -1) {
	write_iter = memcpya(write_iter, loadbuf, (uintptr_t)(chr_code_end - loadbuf));
      } else {
	write_iter = chr_name_write(cip, chr_code_base, write_iter);
      }
      *write_iter++ = '\t';
      write_iter = uint32toa(cur_bp, write_iter);

      uint32_t alt_ct = 1;
      char* copy_start = pos_str_end;
      while (1) {
	++loadbuf_iter;
	unsigned char ucc;
	do {
	  ucc = (unsigned char)(*(++loadbuf_iter));
	  // allow GATK 3.4 <*:DEL> symbolic allele
	} while ((ucc > ',') || (ucc == '*'));

	// temporary kludge	
	if (ucc == ',') {
	  alt_ct = 2;
	  break;
	}
	
	write_iter = memcpya(write_iter, copy_start, (uintptr_t)(loadbuf_iter - copy_start));
	// unsafe to flush for now due to multiallelic kludge
	/*
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	    goto vcf_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	*/
	if (ucc != ',') {
	  break;
	}
	copy_start = loadbuf_iter;
	++alt_ct;
      }

      // temporary kludge
      if (alt_ct > 1) {
	write_iter = write_line_start;
	continue;
      }
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	  goto vcf_to_pgen_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }

      write_iter = memcpya(write_iter, loadbuf_iter, (uintptr_t)((info_nonpr_present ? info_end : filter_end) - loadbuf_iter));
      append_binary_eoln(&write_iter);

      if (gt_missing) {
	fill_all_bits(2 * sample_ct, genovec);
	if (spgw_append_biallelic_genovec(genovec, &spgw)) {
	  goto vcf_to_pgen_ret_WRITE_FAIL;
	}
      } else {
	loadbuf_iter = (char*)rawmemchr(format_start, '\t');
	uint32_t qual_field_skips[2];
	int32_t qual_thresholds[2];
	if (format_gq_or_dp_relevant) {
	  qual_field_ct = vcf_qual_scan_init(vcf_min_gq, vcf_min_dp, format_start, loadbuf_iter, qual_field_skips, qual_thresholds);
	}

	uint32_t dosage_field_idx = 0;
	dosage_t* dosage_vals_iter = dosage_vals;
	if (dosage_flags && IS_SET(dosage_flags, vidx)) {
	  dosage_field_idx = get_vcf_format_position(dosage_import_field, dosage_import_field_slen, format_start, loadbuf_iter);
	}

	// todo: multiallelic variants
	++loadbuf_iter;
	const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
	uint32_t inner_loop_last = kBitsPerWordD2 - 1;
	uint32_t widx = 0;
	if (!IS_SET(phasing_flags, vidx)) {
	  while (1) {
	    if (widx >= sample_ctl2_m1) {
	      if (widx > sample_ctl2_m1) {
		break;
	      }
	      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	    }
	    uintptr_t genovec_word = 0;
	    uint32_t dosage_present_hw = 0;
	    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	      char* cur_gtext_end = strchr(loadbuf_iter, '\t');
	      if (!cur_gtext_end) {
		if ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)) {
		  goto vcf_to_pgen_ret_MISSING_TOKENS;
		}
		cur_gtext_end = next_prespace(loadbuf_iter);
	      }
	      uintptr_t cur_geno;
	      if (qual_field_ct) {
		if (vcf_check_quals(qual_field_skips, qual_thresholds, loadbuf_iter, cur_gtext_end, qual_field_ct)) {
		  goto vcf_to_pgen_force_missing_1;
		}
	      }
	      {
		const uint32_t is_haploid = (loadbuf_iter[1] != '/') && (loadbuf_iter[1] != '|');
		cur_geno = (unsigned char)(*loadbuf_iter) - 48;
		if (cur_geno <= 1) {
		  if (is_haploid) {
		    cur_geno *= 2;
		  } else {
		    const char cc = loadbuf_iter[3];
		    if (((cc != '/') && (cc != '|')) || (loadbuf_iter[4] == '.')) {
		      // code triploids, etc. as missing
		      // might want to subject handling of 0/0/. to
		      // --vcf-half-call control
		      const uintptr_t second_allele_idx = ((unsigned char)loadbuf_iter[2]) - 48;
		      if (second_allele_idx <= 1) {
			cur_geno += second_allele_idx;
		      } else if (second_allele_idx != (uintptr_t)((intptr_t)(-2))) {
			// not '.'
			goto vcf_to_pgen_ret_INVALID_GT;
		      } else if (vcf_half_call == kVcfHalfCallMissing) {
			cur_geno = 3;
		      } else if (vcf_half_call == kVcfHalfCallError) {
			goto vcf_to_pgen_ret_HALF_CALL_ERROR;
		      } else {
			// kVcfHalfCallHaploid, kVcfHalfCallReference
			cur_geno <<= vcf_half_call;
		      }
		    }
		  }
		} else if (cur_geno != (uintptr_t)((intptr_t)(-2))) {
		  // not '.'
		  goto vcf_to_pgen_ret_INVALID_GT;
		} else if (vcf_half_call != kVcfHalfCallMissing) {
		  const char second_allele_char = loadbuf_iter[2];
		  if ((second_allele_char != '.') && ((loadbuf_iter[1] == '/') || (loadbuf_iter[1] == '|'))) {
		    cur_geno = ((unsigned char)second_allele_char) - 48;
		    if (cur_geno > 1) {
		      goto vcf_to_pgen_ret_INVALID_GT;
		    }
		    if (vcf_half_call == kVcfHalfCallError) {
		      goto vcf_to_pgen_ret_HALF_CALL_ERROR;
		    }
		    // kVcfHalfCallHaploid, kVcfHalfCallReference
		    cur_geno <<= vcf_half_call;
		  } else {
		    cur_geno = 3;
		  }
		} else {
		  cur_geno = 3;
		}
		if (dosage_field_idx) {
		  uint32_t is_missing = 0;
		  uint32_t dosage_int;
		  if (!parse_vcf_dosage(loadbuf_iter, cur_gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &is_missing, &dosage_int)) {
		    const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
		    if ((cur_geno == 3) && (cur_halfdist >= hard_call_halfdist)) {
		      // only overwrite GT if (i) it was missing, and (ii) the
		      // dosage's distance from the nearest hardcall doesn't
		      // exceed the --hard-call-threshold value.
		      // (possible todo: warn or error out if dosage and GT are
		      // inconsistent)
		      cur_geno = (dosage_int + (kDosage4th * k1LU)) / kDosageMid;
		    }
		    if (cur_halfdist < dosage_erase_halfdist) {
		      dosage_present_hw |= 1U << sample_idx_lowbits;
		      *dosage_vals_iter++ = dosage_int;
		    }
		  } else if (!is_missing) {
		    goto vcf_to_pgen_ret_INVALID_DOSAGE;
		  }
		}
	      }
	      while (0) {
	      vcf_to_pgen_force_missing_1:
		cur_geno = 3;
	      }
	      genovec_word |= cur_geno << (2 * sample_idx_lowbits);
	      loadbuf_iter = &(cur_gtext_end[1]);
	    }
	    genovec[widx] = genovec_word;
	    if (dosage_field_idx) {
	      ((halfword_t*)dosage_present)[widx] = dosage_present_hw;
	    }
	    ++widx;
	  }
	  if (!dosage_field_idx) {
	    if (spgw_append_biallelic_genovec(genovec, &spgw)) {
	      goto vcf_to_pgen_ret_WRITE_FAIL;
	    }
	  } else {
	    assert(dosage_vals_iter != dosage_vals);
	    if (spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_vals_iter - dosage_vals, &spgw)) {
	      goto vcf_to_pgen_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  while (1) {
	    if (widx >= sample_ctl2_m1) {
	      if (widx > sample_ctl2_m1) {
		break;
	      }
	      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	    }
	    uintptr_t genovec_word = 0;
	    uint32_t phasepresent_halfword = 0;
	    uint32_t phaseinfo_halfword = 0;
	    uint32_t dosage_present_hw = 0;
	    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	      char* cur_gtext_end = strchr(loadbuf_iter, '\t');
	      if (!cur_gtext_end) {
		if ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)) {
		  goto vcf_to_pgen_ret_MISSING_TOKENS;
		}
		cur_gtext_end = next_prespace(loadbuf_iter);
	      }
	      uintptr_t cur_geno;
	      if (qual_field_ct) {
		if (vcf_check_quals(qual_field_skips, qual_thresholds, loadbuf_iter, cur_gtext_end, qual_field_ct)) {
		  goto vcf_to_pgen_force_missing_2;
		}
	      }
	      {
		const uint32_t is_phased = (loadbuf_iter[1] == '|');
		const uint32_t is_haploid = (!is_phased) && (loadbuf_iter[1] != '/');
		cur_geno = (unsigned char)(*loadbuf_iter) - 48;
		if (cur_geno <= 1) {
		  if (is_haploid) {
		    cur_geno *= 2;
		  } else {
		    const char cc = loadbuf_iter[3];
		    if (((cc != '/') && (cc != '|')) || (loadbuf_iter[4] == '.')) {
		      // code triploids, etc. as missing
		      // might want to subject handling of 0/0/. to
		      // --vcf-half-call control
		      const uintptr_t second_allele_idx = ((unsigned char)loadbuf_iter[2]) - 48;
		      if (second_allele_idx <= 1) {
			cur_geno += second_allele_idx;
			// todo: check if this should be less branchy
			if (is_phased && (cur_geno == 1)) {
			  const uint32_t shifted_bit = 1U << sample_idx_lowbits;
			  phasepresent_halfword |= shifted_bit;
			  if (!second_allele_idx) {
			    // 1|0
			    phaseinfo_halfword |= shifted_bit;
			  }
			}
		      } else if (second_allele_idx != (uintptr_t)((intptr_t)(-2))) {
			// not '.'
			goto vcf_to_pgen_ret_INVALID_GT;
		      } else if (vcf_half_call == kVcfHalfCallMissing) {
			cur_geno = 3;
		      } else if (vcf_half_call == kVcfHalfCallError) {
			goto vcf_to_pgen_ret_HALF_CALL_ERROR;
		      } else {
			// kVcfHalfCallHaploid, kVcfHalfCallReference
			cur_geno <<= vcf_half_call;
		      }
		    }
		  }
		} else if (cur_geno != (uintptr_t)((intptr_t)(-2))) {
		  // not '.'
		  goto vcf_to_pgen_ret_INVALID_GT;
		} else if (vcf_half_call != kVcfHalfCallMissing) {
		  const char second_allele_char = loadbuf_iter[2];
		  if ((second_allele_char != '.') && ((loadbuf_iter[1] == '/') || (loadbuf_iter[1] == '|'))) {
		    cur_geno = ((unsigned char)second_allele_char) - 48;
		    if (cur_geno > 1) {
		      goto vcf_to_pgen_ret_INVALID_GT;
		    }
		    if (vcf_half_call == kVcfHalfCallError) {
		      goto vcf_to_pgen_ret_HALF_CALL_ERROR;
		    }
		    // kVcfHalfCallHaploid, kVcfHalfCallReference
		    cur_geno <<= vcf_half_call;
		  } else {
		    cur_geno = 3;
		  }
		} else {
		  cur_geno = 3;
		}
		if (dosage_field_idx) {
		  uint32_t is_missing = 0;
		  uint32_t dosage_int;
		  if (!parse_vcf_dosage(loadbuf_iter, cur_gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &is_missing, &dosage_int)) {
		    const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
		    if ((cur_geno == 3) && (cur_halfdist >= hard_call_halfdist)) {
		      cur_geno = (dosage_int + (kDosage4th * k1LU)) / kDosageMid;
		    }
		    if (cur_halfdist < dosage_erase_halfdist) {
		      dosage_present_hw |= 1U << sample_idx_lowbits;
		      *dosage_vals_iter++ = dosage_int;
		    }
		  } else if (!is_missing) {
		    goto vcf_to_pgen_ret_INVALID_DOSAGE;
		  }
		}
	      }
	      while (0) {
	      vcf_to_pgen_force_missing_2:
		cur_geno = 3;
	      }
	      genovec_word |= cur_geno << (2 * sample_idx_lowbits);
	      loadbuf_iter = &(cur_gtext_end[1]);
	    }
	    genovec[widx] = genovec_word;
	    ((halfword_t*)phasepresent)[widx] = (halfword_t)phasepresent_halfword;
	    ((halfword_t*)phaseinfo)[widx] = (halfword_t)phaseinfo_halfword;
	    if (dosage_field_idx) {
	      ((halfword_t*)dosage_present)[widx] = dosage_present_hw;
	    }
	    ++widx;
	  }
	  if (!dosage_field_idx) {
	    if (spgw_append_biallelic_genovec_hphase(genovec, phasepresent, phaseinfo, &spgw)) {
	      goto vcf_to_pgen_ret_WRITE_FAIL;
	    }
	  } else {
	    if (spgw_append_biallelic_genovec_hphase_dosage16(genovec, phasepresent, phaseinfo, dosage_present, dosage_vals, dosage_vals_iter - dosage_vals, &spgw)) {
	      goto vcf_to_pgen_ret_WRITE_FAIL;
	    }
	  }
	}
      }
      if (!(++vidx % 1000)) {
	printf("\r--vcf: %uk variants converted.", vidx / 1000);
	fflush(stdout);
      }
    }
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	goto vcf_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&pvarfile)) {
      goto vcf_to_pgen_ret_WRITE_FAIL;
    }
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--vcf: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar");
    if (!preexisting_psamname) {
      write_iter = memcpyl3a(write_iter, " + ");
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya(write_iter, ".psam");
    }
    strcpy(write_iter, " written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  vcf_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  vcf_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  vcf_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  vcf_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  vcf_to_pgen_ret_HALF_CALL_ERROR:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --vcf file has a GT half-call.\n", line_idx);
    if (!vcf_half_call_explicit_error) {
      logerrprint("Use --vcf-half-call to specify how these should be processed.\n");
    }
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_INVALID_GT:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --vcf file has an invalid GT field.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --vcf file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_LONG_LINE_N:
    putc_unlocked('\n', stdout);
  vcf_to_pgen_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --vcf file is pathologically long.\n", line_idx);
  vcf_to_pgen_ret_MALFORMED_INPUT_2N:
    logprint("\n");
    logerrprintb();
  vcf_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_INVALID_DOSAGE:
    putc_unlocked('\n', stdout);
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of --vcf file has an invalid %s field.\n", line_idx, dosage_import_field);
  vcf_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 vcf_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  gzclose_cond(gz_infile);
  fclose_cond(pvarfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

pglerr_t ox_sample_to_psam(const char* samplename, const char* ox_missing_code, misc_flags_t misc_flags, char* outname, char* outname_end, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  FILE* psamfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  {
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (misc_flags & kfMiscKeepAutoconv) {
      // must use --output-missing-phenotype parameter, which we've validated
      // to be consistent with --input-missing-phenotype
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      // use "NA" since that's always safe
      memcpy(output_missing_pheno, "NA", 2);
    }
    const char* missing_catname = g_missing_catname;
    uint32_t missing_catname_slen = strlen(missing_catname);
    
    gz_infile = gzopen(samplename, FOPEN_RB);
    if (!gz_infile) {
      const uint32_t slen = strlen(samplename);
      if (((slen > 7) && (!memcmp(&(samplename[slen - 7]), ".sample", 7))) || ((slen > 10) && (!memcmp(&(samplename[slen - 10]), ".sample.gz", 10)))) {
	LOGERRPRINTFWW(g_errstr_fopen, samplename);
      } else {
	LOGERRPRINTFWW("Error: Failed to open %s. (--sample expects a complete filename; did you forget '.sample' at the end?)\n", samplename);
      }
      goto ox_sample_to_psam_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto ox_sample_to_psam_ret_NOMEM;
    }
    uint32_t mc_ct = 0;
    uintptr_t max_mc_blen = 1;
    char* sorted_mc = nullptr;
    if (!ox_missing_code) {
      if (bigstack_alloc_c(3, &sorted_mc)) {
	goto ox_sample_to_psam_ret_NOMEM;
      }
      memcpy(sorted_mc, "NA", 3);
      mc_ct = 1;
      max_mc_blen = 3;
    } else {
      // er, this should use something like
      // count_and_measure_multistr_reverse_alloc()...
      const char* missing_code_iter = ox_missing_code;
      while (*missing_code_iter) {
	while (*missing_code_iter == ',') {
	  ++missing_code_iter;
	}
	if (!(*missing_code_iter)) {
	  break;
	}
	++mc_ct;
	const char* token_end = strchr(missing_code_iter, ',');
	if (!token_end) {
	  token_end = (const char*)rawmemchr(missing_code_iter, '\0');
	}
	uintptr_t token_slen = (uintptr_t)(token_end - missing_code_iter);
	if (token_slen >= max_mc_blen) {
	  max_mc_blen = token_slen + 1;
	}
	missing_code_iter = token_end;
      }
      if (mc_ct) {
	if (bigstack_alloc_c(mc_ct * max_mc_blen, &sorted_mc)) {
	  goto ox_sample_to_psam_ret_NOMEM;
	}
	missing_code_iter = ox_missing_code;
	for (uintptr_t mc_idx = 0; mc_idx < mc_ct; ++mc_idx) {
	  while (*missing_code_iter == ',') {
	    ++missing_code_iter;
	  }
	  const char* token_end = strchr(missing_code_iter, ',');
	  if (!token_end) {
	    token_end = (const char*)rawmemchr(missing_code_iter, '\0');
	  }
	  uintptr_t token_slen = (uintptr_t)(token_end - missing_code_iter);
	  memcpyx(&(sorted_mc[mc_idx * max_mc_blen]), missing_code_iter, token_slen, '\0');
	  missing_code_iter = token_end;
	}
	qsort(sorted_mc, mc_ct, max_mc_blen, strcmp_casted);
      }
    }
    loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto ox_sample_to_psam_ret_NOMEM;
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
	  goto ox_sample_to_psam_ret_READ_FAIL;
	}
	logerrprint("Error: Empty .sample file.\n");
	goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));
    char* token_end = token_endnn(loadbuf_first_token);
    if ((((uintptr_t)(token_end - loadbuf_first_token)) != 4) || memcmp(loadbuf_first_token, "ID_1", 4)) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1;
    }
    // currently accepts tab as delimiter, though .sample spec technically
    // prohibits that
    char* loadbuf_iter = skip_initial_spaces(token_end);
    uint32_t token_slen = strlen_se(loadbuf_iter);
    if ((token_slen != 4) || memcmp(loadbuf_iter, "ID_2", 4)) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1;      
    }
    loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[token_slen]));
    token_slen = strlen_se(loadbuf_iter);
    if ((token_slen != 7) || (!match_upper_counted(loadbuf_iter, "MISSING", 7))) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1;
    }
    loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[token_slen]));

    strcpy(outname_end, ".psam");
    if (fopen_checked(outname, FOPEN_WB, &psamfile)) {
      goto ox_sample_to_psam_ret_OPEN_FAIL;
    }
    // categorical phenotypes are lengthened by 1 character ('C' added in
    // front), so this needs to be 50% larger than loadbuf to handle worst case
    char* writebuf = (char*)bigstack_alloc_raw_rd(loadbuf_size + (loadbuf_size / 2));
    char* write_iter = strcpya(writebuf, "#FID\tIID\tSEX");
    
    // 0 = not present, otherwise zero-based index (this is fine since first
    //     column has to be FID)
    uint32_t sex_col = 0;

    uint32_t col_ct = 3;
    
    while (!is_eoln_kns(*loadbuf_iter)) {      
      token_end = token_endnn(loadbuf_iter);
      token_slen = (uintptr_t)(token_end - loadbuf_iter);
      if ((token_slen == 3) && match_upper_counted(loadbuf_iter, "SEX", 3)) {
	if (sex_col) {
	  logerrprint("Error: Multiple sex columns in .sample file.\n");
	  goto ox_sample_to_psam_ret_MALFORMED_INPUT;
	}
	sex_col = col_ct;
      }
      ++col_ct;
      loadbuf_iter = skip_initial_spaces(token_end);
    }

    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto ox_sample_to_psam_ret_READ_FAIL;
	}
	logerrprint("Error: Only one nonempty line in .sample file.\n");
	goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));

    token_end = token_endnn(loadbuf_first_token);
    if ((((uintptr_t)(token_end - loadbuf_first_token)) != 1) || (*loadbuf_first_token != '0')) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
    }
    loadbuf_iter = skip_initial_spaces(token_end);
    token_slen = strlen_se(loadbuf_iter);
    if ((token_slen != 1) || (*loadbuf_iter != '0')) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
    }
    loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[1]));
    token_slen = strlen_se(loadbuf_iter);
    if ((token_slen != 1) || (*loadbuf_iter != '0')) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
    }
    loadbuf_iter = &(loadbuf_iter[1]);

    const uint32_t col_ctl = BITCT_TO_WORDCT(col_ct);
    uintptr_t* col_is_categorical = (uintptr_t*)bigstack_alloc_raw_rd(col_ctl * sizeof(intptr_t));
    uintptr_t* col_is_qt = (uintptr_t*)bigstack_alloc_raw_rd(col_ctl * sizeof(intptr_t));
    fill_ulong_zero(col_ctl, col_is_categorical);
    fill_ulong_zero(col_ctl, col_is_qt);
    uint32_t at_least_one_binary_pheno = 0;
    for (uint32_t col_idx = 3; col_idx < col_ct; ++col_idx) {
      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
      unsigned char col_type_char = *loadbuf_iter;
      if (is_eoln_kns(col_type_char)) {
	logerrprint("Error: Second .sample header line has fewer tokens than the first.\n");
	goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }
      if (loadbuf_iter[1] > ' ') {
	goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
      }
      if (col_idx == sex_col) {
	if (col_type_char != 'D') {
	  logerrprint("Error: .sample sex column is not of type 'D'.\n");
	  goto ox_sample_to_psam_ret_MALFORMED_INPUT;
	}
      } else {
	if ((col_type_char == 'C') || (col_type_char == 'P')) {
	  set_bit(col_idx, col_is_qt);
	} else if (col_type_char == 'D') {
	  set_bit(col_idx, col_is_categorical);
	} else {
	  at_least_one_binary_pheno = 1;
	  if (col_type_char != 'B') {
	    sprintf(g_logbuf, "Error: Unrecognized .sample variable type '%c'.\n", col_type_char);
	    goto ox_sample_to_psam_ret_MALFORMED_INPUT_2;
	  }
	}
      }
      ++loadbuf_iter;
    }
    if (at_least_one_binary_pheno) {
      // check for pathological case
      if ((bsearch_str("0", sorted_mc, 1, max_mc_blen, mc_ct) != -1) || (bsearch_str("1", sorted_mc, 1, max_mc_blen, mc_ct) != -1)) {
	logerrprint("Error: '0' and '1' are unacceptable missing case/control phenotype codes.\n");
	goto ox_sample_to_psam_ret_INCONSISTENT_INPUT;
      }
    }
    // to make --data and --data --make-pgen consistent, we do a two-pass load,
    // checking for empty phenotypes in the first pass.
    uintptr_t* col_keep = (uintptr_t*)bigstack_alloc_raw(round_up_pow2(col_ct, kCacheline * CHAR_BIT) / CHAR_BIT);
    col_keep[0] = 7;
    fill_ulong_zero(col_ctl - 1, &(col_keep[1]));
    uint32_t uncertain_col_ct = col_ct - 3;
    if (sex_col) {
      // we don't care if sex column is all-NA
      set_bit(sex_col, col_keep);
      --uncertain_col_ct;
    }
    while (uncertain_col_ct) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_iter = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_iter)) {
	continue;
      }

      const uint32_t old_uncertain_col_ct = uncertain_col_ct;
      uint32_t old_col_uidx = 0;
      uint32_t col_uidx = 0;
      for (uint32_t uncertain_col_idx = 0; uncertain_col_idx < old_uncertain_col_ct; ++uncertain_col_idx, ++col_uidx) {
	next_unset_unsafe_ck(col_keep, &col_uidx);
	loadbuf_iter = next_token_mult(loadbuf_iter, col_uidx - old_col_uidx);
	if (!loadbuf_iter) {
	  goto ox_sample_to_psam_ret_MISSING_TOKENS;
	}
	token_end = token_endnn(loadbuf_iter);
	token_slen = (uintptr_t)(token_end - loadbuf_iter);
        if (bsearch_str(loadbuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) == -1) {
	  set_bit(col_uidx, col_keep);
	  --uncertain_col_ct;
	}
	loadbuf_iter = token_end;
	old_col_uidx = col_uidx;
      }
    }

    if (gzrewind(gz_infile)) {
      goto ox_sample_to_psam_ret_READ_FAIL;
    }
    line_idx = 0;
    
    uint32_t sample_ct_p2 = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_first_token)) {
	continue;
      }
      ++sample_ct_p2;
      if (sample_ct_p2 < 3) {
	// header lines
	if (sample_ct_p2 == 1) {
	  loadbuf_iter = next_token_mult(loadbuf_first_token, 3);
	  for (uint32_t col_idx = 3; col_idx < col_ct; ++col_idx) {
	    token_end = token_endnn(loadbuf_iter);
	    if (is_set(col_keep, col_idx) && (col_idx != sex_col)) {
	      *write_iter++ = '\t';
	      write_iter = memcpya(write_iter, loadbuf_iter, (uintptr_t)(token_end - loadbuf_iter));
	    }
	    loadbuf_iter = skip_initial_spaces(token_end);
	  }
	  append_binary_eoln(&write_iter);
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), psamfile)) {
	    goto ox_sample_to_psam_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	continue;
      }
      if (sample_ct_p2 == 0x80000001U) {
	logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
	goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }

      // FID
      token_end = token_endnn(loadbuf_first_token);
      write_iter = memcpyax(writebuf, loadbuf_first_token, (uintptr_t)(token_end - loadbuf_first_token), '\t');

      // IID
      loadbuf_iter = skip_initial_spaces(token_end);
      if (is_eoln_kns(*loadbuf_iter)) {
	goto ox_sample_to_psam_ret_MISSING_TOKENS;
      }
      token_end = token_endnn(loadbuf_iter);
      write_iter = memcpya(write_iter, loadbuf_iter, (uintptr_t)(token_end - loadbuf_iter));
      
      // MISSING
      loadbuf_iter = skip_initial_spaces(token_end);
      if (is_eoln_kns(*loadbuf_iter)) {
	goto ox_sample_to_psam_ret_MISSING_TOKENS;
      }
      token_end = token_endnn(loadbuf_iter);

      // flush now since backfilled sex is variable-length ("NA" vs. "1"/"2")
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), psamfile)) {
	goto ox_sample_to_psam_ret_WRITE_FAIL;
      }
      char* cur_writebuf_start = writebuf;
      write_iter = memcpyl3a(writebuf, "\tNA");
      for (uint32_t col_idx = 3; col_idx < col_ct; ++col_idx) {
	loadbuf_iter = skip_initial_spaces(token_end);
	if (is_eoln_kns(*loadbuf_iter)) {
	  goto ox_sample_to_psam_ret_MISSING_TOKENS;
	}
	token_end = token_endnn(loadbuf_iter);
	if (!is_set(col_keep, col_idx)) {
	  continue;
	}
	token_slen = (uintptr_t)(token_end - loadbuf_iter);
	const uint32_t is_missing = (bsearch_str(loadbuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) != -1);
	if (col_idx == sex_col) {
	  if (!is_missing) {
	    const unsigned char sex_ucc = *loadbuf_iter;
	    if ((token_slen == 1) && ((((uint32_t)sex_ucc) - 49) < 2)) {
	      ++cur_writebuf_start;
	      cur_writebuf_start[0] = '\t';
	      cur_writebuf_start[1] = sex_ucc;
	    } else if ((token_slen != 1) || (sex_ucc != '0')) {
	      // tolerate '0' as a sex-only missing code even when not
	      // explicitly specified
	      *token_end = '\0';
	      sprintf(g_logbuf, "Error: Invalid sex code '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', '2', or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
	      goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
	    }
	  }
	} else {
	  *write_iter++ = '\t';
	  if (is_set(col_is_categorical, col_idx)) {
	    if (!is_missing) {
	      *write_iter++ = 'C';
	      // .sample files are relatively small, so let's go ahead and
	      // (i) validate we have a positive integer < 2^31
	      // (ii) convert e.g. 9000000, 9000000., 9.0e6 all to 9000000
	      double dxx = 0.0;
	      char* num_end = scanadv_double(loadbuf_iter, &dxx);
	      int32_t ii = (int32_t)dxx;
	      if ((num_end != token_end) || (ii <= 0) || (((double)ii) != dxx)) {
		*token_end = '\0';
		sprintf(g_logbuf, "Error: Invalid categorical phenotype '%s' on line %" PRIuPTR ", column %u of .sample file (positive integer < 2^31 or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
		goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
	      }
	      write_iter = uint32toa(ii, write_iter);
	    } else {
	      write_iter = memcpya(write_iter, missing_catname, missing_catname_slen);
	    }
	  } else if (!is_missing) {
	    if (is_set(col_is_qt, col_idx)) {
	      double dxx = 0.0;
	      char* num_end = scanadv_double(loadbuf_iter, &dxx);
	      if (num_end != token_end) {
		*token_end = '\0';
		sprintf(g_logbuf, "Error: Invalid quantitative phenotype '%s' on line %" PRIuPTR ", column %u of .sample file (non-infinite number or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
		goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
	      }
	      // used over memcpy to make --data and --data --make-pgen the
	      // same (could make that conditional on keep_autoconv?)
	      write_iter = dtoa_g(dxx, write_iter);
	    } else {
	      const uint32_t cc_char_m48 = (uint32_t)((unsigned char)(*loadbuf_iter)) - 48;
	      if ((token_slen == 1) && (cc_char_m48 < 2)) {
		*write_iter++ = cc_char_m48 + '1';
	      } else {
		*token_end = '\0';
		sprintf(g_logbuf, "Error: Invalid binary phenotype value '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
		goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
	      }
	    }
	  } else {
	    write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
	  }
	}
      }
      append_binary_eoln(&write_iter);
      if (fwrite_checked(cur_writebuf_start, (uintptr_t)(write_iter - cur_writebuf_start), psamfile)) {
	goto ox_sample_to_psam_ret_WRITE_FAIL;
      }
    }
    if (!gzeof(gz_infile)) {
      goto ox_sample_to_psam_ret_READ_FAIL;
    }

    // no final writebuf flush since we didn't use usual manual-streaming
    // strategy
    if (fclose_null(&psamfile)) {
      goto ox_sample_to_psam_ret_WRITE_FAIL;
    }
    const uint32_t sample_ct = sample_ct_p2 - 2;
    if ((!sample_ct) && (!(misc_flags & kfMiscAllowNoSamples))) {
      logerrprint("Error: No samples in .sample file.\n");
      goto ox_sample_to_psam_ret_INCONSISTENT_INPUT;
    }
    LOGPRINTFWW("%u sample%s imported from .sample file to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  ox_sample_to_psam_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file is pathologically long.\n", line_idx);
      reterr = kPglRetMalformedInput;
      break;
    }
  ox_sample_to_psam_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_sample_to_psam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_sample_to_psam_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_sample_to_psam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ox_sample_to_psam_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_MALFORMED_INPUT_2:
    logerrprintb();
  ox_sample_to_psam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1:
    logerrprint("Error: Invalid first header line in .sample file.\n");
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2:
    logerrprint("Error: Invalid second header line in .sample file.\n");
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  ox_sample_to_psam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  gzclose_cond(gz_infile);
  fclose_cond(psamfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

uint32_t bgen11_dosage_import_check(uint32_t dosage_int_sum_thresh, uint32_t import_dosage_certainty_int, uint32_t dosage_erase_halfdist, uint32_t dosage_int0, uint32_t dosage_int1, uint32_t dosage_int2) {
  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
  if ((dosage_int_sum <= dosage_int_sum_thresh) && (dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
    return 0;
  }
  // ties realistically happen, use banker's rounding
  // 1/65536 -> 0/32768
  // 3/65536 -> 2/32768
  // 5/65536 -> 2/32768
  const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
  uint32_t write_dosage_int;
  if (dosage_int_sum == kDosageMax) {
    // optimize common case
    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
  } else {
    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
  }
  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
  return (halfdist < dosage_erase_halfdist);
}

void bgen11_dosage_import_update(uint32_t dosage_int_sum_thresh, uint32_t import_dosage_certainty_int, uint32_t hard_call_halfdist, uint32_t dosage_erase_halfdist, uint32_t sample_idx_lowbits, uint32_t dosage_int0, uint32_t dosage_int1, uint32_t dosage_int2, uintptr_t* genovec_word_ptr, uint32_t* dosage_present_hw_ptr, dosage_t** dosage_vals_iterp) {
  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
  if (dosage_int_sum <= dosage_int_sum_thresh) {
    if ((dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
      *genovec_word_ptr |= (3 * k1LU) << (2 * sample_idx_lowbits);
      return;
    }
  }
  const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
  uint32_t write_dosage_int;
  if (dosage_int_sum == kDosageMax) {
    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
  } else {
    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
  }
  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
  if (halfdist < hard_call_halfdist) {
    *genovec_word_ptr |= (3 * k1LU) << (2 * sample_idx_lowbits);
  } else {
    *genovec_word_ptr |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
    if (halfdist >= dosage_erase_halfdist) {
      return;
    }
  }
  *dosage_present_hw_ptr |= 1U << sample_idx_lowbits;
  **dosage_vals_iterp = write_dosage_int;
  *dosage_vals_iterp += 1;
}

static_assert(sizeof(dosage_t) == 2, "ox_gen_to_pgen() needs to be updated.");
pglerr_t ox_gen_to_pgen(const char* genname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  FILE* pvarfile = nullptr;
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  spgw_preinit(&spgw);
  {
    uint32_t sample_ct;
    reterr = ox_sample_to_psam(samplename, ox_missing_code, misc_flags, outname, outname_end, &sample_ct);
    if (reterr) {
      goto ox_gen_to_pgen_ret_1;
    }
    if (sample_ct > (kMaxLongLine / 6)) {
      // impossible for a valid .gen line to fit in maximum-length load buffer
      logerrprint("Error: Too many samples for .gen file converter.\n");
      reterr = kPglRetNotYetSupported;
      goto ox_gen_to_pgen_ret_1;
    }
    // Two passes:
    // 1. Count # of (non-chromosome-filtered) variants, write .pvar file, and
    //    check if at least one non-hardcall needs to be saved.
    // 2. Write the .pgen.
    gz_infile = gzopen(genname, FOPEN_RB);
    if (!gz_infile) {
      const uint32_t slen = strlen(genname);
      if (((slen > 4) && (!memcmp(&(genname[slen - 4]), ".gen", 4))) || ((slen > 7) && (!memcmp(&(genname[slen - 7]), ".gen.gz", 7)))) {
	LOGERRPRINTFWW(g_errstr_fopen, genname);
      } else {
	LOGERRPRINTFWW("Error: Failed to open %s. (--gen expects a complete filename; did you forget '.gen' at the end?)\n", genname);
      }
      goto ox_gen_to_pgen_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto ox_gen_to_pgen_ret_NOMEM;
    }    
    loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto ox_gen_to_pgen_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    finalize_chrset(misc_flags, cip);

    char* writebuf = (char*)bigstack_alloc_raw(kCompressStreamBlock + loadbuf_size);
    char* writebuf_flush = &(writebuf[kCompressStreamBlock]);

    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    if (ox_single_chr_str) {
      int32_t chr_code_raw = get_chr_code_raw(ox_single_chr_str);
      if (chr_code_raw == -1) {
	// command-line parser guarantees that allow_extra_chrs is true here
	single_chr_str = ox_single_chr_str;
        single_chr_slen = strlen(ox_single_chr_str);
      } else {
	uint32_t chr_code = chr_code_raw;
	if (chr_code > cip->max_code) {
	  if (chr_code < kMaxContigs) {
	    logerrprint("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
	    goto ox_gen_to_pgen_ret_INVALID_CMDLINE;
	  }
	  chr_code = cip->xymt_codes[chr_code - kMaxContigs];
	  if (((int32_t)chr_code) < 0) {
	    logerrprint("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
	    goto ox_gen_to_pgen_ret_INVALID_CMDLINE;
	  }
	}
	if (!is_set(cip->chr_mask, chr_code)) {
	  // could permit this in --allow-no-vars case, but it's silly
	  logerrprint("Error: --oxford-single-chr chromosome code is excluded by chromosome filter.\n");
	  goto ox_gen_to_pgen_ret_INVALID_CMDLINE;
	}
	char* chr_buf = (char*)bigstack_alloc_raw(kCacheline);
	char* chr_name_end = chr_name_write(cip, chr_code, chr_buf);
	single_chr_str = chr_buf;
        single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
      }
    }

    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &pvarfile)) {
      goto ox_gen_to_pgen_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    write_iter = strcpya(write_iter, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    // Explicit 32768 instead of kDosageMax since this is driven by the BGEN
    // 1.1 format, not plink2's dosage representation.
    // Note that command-line parser multiplies import_dosage_certainty by
    // (1 - kSmallEpsilon), and we want import_dosage_certainty_int to be 1
    // when import_dosage_certainty is zero.
    uint32_t import_dosage_certainty_int = 1 + (int32_t)(import_dosage_certainty * 32768);
    const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
    
    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t dosage_is_present = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto ox_gen_to_pgen_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_gen_to_pgen_ret_LONG_LINE_N;
      }
      char* chr_code_str = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*chr_code_str)) {
	continue;
      }
      char* chr_code_end = token_endnn(chr_code_str);
      char* variant_id_str = skip_initial_spaces(chr_code_end);
      if (is_eoln_kns(*variant_id_str)) {
	goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      
      if (!single_chr_str) {
	int32_t cur_chr_code;
        reterr = get_or_add_chr_code_destructive(".gen file", line_idx, allow_extra_chrs, chr_code_str, chr_code_end, cip, &cur_chr_code);
	if (reterr) {
	  if ((((uintptr_t)(chr_code_end - chr_code_str)) == 3) && (!memcmp(chr_code_str, "---", 3))) {
	    logerrprint("(Did you forget --oxford-single-chr?)\n");
	  }
	  goto ox_gen_to_pgen_ret_1;
	}
	if (!is_set(cip->chr_mask, cur_chr_code)) {
	  ++variant_skip_ct;
	  continue;
	}
	write_iter = chr_name_write(cip, cur_chr_code, write_iter);
      } else {
	write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
      }
      *write_iter++ = '\t';
      ++variant_ct;

      char* variant_id_end = token_endnn(variant_id_str);
      char* pos_str = skip_initial_spaces(variant_id_end);
      if (is_eoln_kns(*pos_str)) {
	goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      char* pos_end = token_endnn(pos_str);
      uint32_t cur_bp;
      if (scan_uint_defcap(pos_str, &cur_bp)) {
	putc_unlocked('\n', stdout);
	sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, genname);
	goto ox_gen_to_pgen_ret_MALFORMED_INPUT_WW;
      }
      
      write_iter = uint32toa_x(cur_bp, '\t', write_iter);
      const uint32_t variant_id_slen = (uintptr_t)(variant_id_end - variant_id_str);
      if (variant_id_slen > kMaxIdSlen) {
	putc_unlocked('\n', stdout);
	logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
	goto ox_gen_to_pgen_ret_MALFORMED_INPUT;
      }
      write_iter = memcpyax(write_iter, variant_id_str, variant_id_slen, '\t');

      // .gen specification does not define which column should be expected to
      // be the reference allele, and which the alternate.  plink 1.9 assumed
      // alt was usually first, but the reverse seems to be more common now.
      // So:
      //   If 'ref-first' or 'ref-second' was specified, we know what to do.
      //   If not, we treat the second allele as the provisional reference, for
      //     backward compatibility.
      char* first_allele_str = skip_initial_spaces(pos_end);
      if (is_eoln_kns(*first_allele_str)) {
	goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      char* first_allele_end = token_endnn(first_allele_str);
      char* second_allele_str = skip_initial_spaces(first_allele_end);
      if (is_eoln_kns(*second_allele_str)) {
	goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      char* loadbuf_iter = token_endnn(second_allele_str);
      if (!prov_ref_allele_second) {
	write_iter = memcpyax(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str), '\t');
	write_iter = memcpya(write_iter, second_allele_str, (uintptr_t)(loadbuf_iter - second_allele_str));
      } else {
	write_iter = memcpyax(write_iter, second_allele_str, (uintptr_t)(loadbuf_iter - second_allele_str), '\t');
	write_iter = memcpya(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str));
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	  goto ox_gen_to_pgen_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
      
      if (!dosage_is_present) {
	for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	  loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	  const char cc = *loadbuf_iter;
	  if (is_eoln_kns(cc)) {
	    goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	  }
	  char cc2 = loadbuf_iter[1];
	  if ((cc2 == ' ') || (cc2 == '\t')) {
	    // fast handling of "1 0 0 ", "0 1 0 ", "0 0 1 ", "0 0 0 " cases
	    cc2 = loadbuf_iter[3];
	    if (((cc2 == ' ') || (cc2 == '\t')) && ((unsigned char)(loadbuf_iter[5]) <= 32)) {
	      const uint32_t uii = ((uint32_t)((unsigned char)cc)) - 48;
	      const uint32_t ujj = ((uint32_t)((unsigned char)loadbuf_iter[2])) - 48;
	      const uint32_t ukk = ((uint32_t)((unsigned char)loadbuf_iter[4])) - 48;
	      if (((uii | ujj | ukk) < 2) && (uii + ujj + ukk < 2)) {
		loadbuf_iter = &(loadbuf_iter[5]);
		continue;
	      }
	    }
	  }
	  double prob_0alt;
	  char* first_dosage_str_end = scanadv_double(loadbuf_iter, &prob_0alt);
	  if (!first_dosage_str_end) {
	    // triple-NA, etc. ok; treat as missing value
	    loadbuf_iter = next_token_mult(loadbuf_iter, 2);
	    if (!loadbuf_iter) {
	      goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	    }
	    loadbuf_iter = token_endnn(loadbuf_iter);
	    continue;
	  }
	  if ((unsigned char)(*first_dosage_str_end) > ' ') {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  loadbuf_iter = skip_initial_spaces(first_dosage_str_end);
	  if (is_eoln_kns(*loadbuf_iter)) {
	    goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	  }
	  double prob_1alt;
	  loadbuf_iter = scanadv_double(loadbuf_iter, &prob_1alt);
	  if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	  if (is_eoln_kns(*loadbuf_iter)) {
	    goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	  }
	  double prob_2alt;
	  loadbuf_iter = scanadv_double(loadbuf_iter, &prob_2alt);
	  if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  // bugfix: forgot the "multiply by 32768" part of "multiply by 32768
	  // and round" .gen -> .bgen conversion.
	  prob_0alt *= 32768;
	  prob_1alt *= 32768;
	  prob_2alt *= 32768;
	  
	  // now treat this identically to bgen-1.1
	  // Compare with 65535.4999999999 instead of 65535.5 since 0.5 +
	  // [first floating point number below 65535.5] may evaluate to 65536.
	  if ((prob_0alt < 0.0) || (prob_0alt >= 65535.4999999999) || (prob_1alt < 0.0) || (prob_1alt >= 65535.4999999999) || (prob_2alt < 0.0) || (prob_2alt >= 65535.4999999999)) {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  const uint32_t dosage_int0 = (int32_t)(prob_0alt + 0.5);
	  const uint32_t dosage_int1 = (int32_t)(prob_1alt + 0.5);
	  const uint32_t dosage_int2 = (int32_t)(prob_2alt + 0.5);
	  dosage_is_present = bgen11_dosage_import_check(dosage_int_sum_thresh, import_dosage_certainty_int, dosage_erase_halfdist, dosage_int0, dosage_int1, dosage_int2);
	  if (dosage_is_present) {
	    break;
	  }
	}
      }
      if (!(variant_ct % 1000)) {
	printf("\r--data/--gen: %uk variants scanned.", variant_ct / 1000);
	fflush(stdout);
      }
    }
    putc_unlocked('\r', stdout);
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	goto ox_gen_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&pvarfile)) {
      goto ox_gen_to_pgen_ret_WRITE_FAIL;
    }
    if (!variant_ct) {
      if (!variant_skip_ct) {
	// permit this in --allow-no-vars case?
	logerrprint("Error: Empty .gen file.\n");
	goto ox_gen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      LOGERRPRINTFWW("Error: All %" PRIuPTR " variant%s in .gen file skipped due to chromosome filter.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s");
      goto ox_gen_to_pgen_ret_INCONSISTENT_INPUT;
    }
    LOGPRINTF("--data/--gen: %u variant%s scanned%s.\n", variant_ct, (variant_ct == 1)? "" : "s", dosage_is_present? "" : " (all hardcalls)");

    // second pass
    bigstack_reset(writebuf);
    if (gzrewind(gz_infile)) {
      goto ox_gen_to_pgen_ret_READ_FAIL;
    }
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto ox_gen_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto ox_gen_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* dosage_present;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
	bigstack_alloc_ul(sample_ctl, &dosage_present)) {
      goto ox_gen_to_pgen_ret_NOMEM;
    }
    dosage_t* dosage_vals = nullptr;
    if (dosage_is_present) {
      if (bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
	goto ox_gen_to_pgen_ret_NOMEM;
      }
    }
    if (hard_call_thresh == 0xffffffffU) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const uintptr_t line_ct = line_idx - 1;
    uint32_t vidx = 0;
    for (line_idx = 1; line_idx <= line_ct; ++line_idx) {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	goto ox_gen_to_pgen_ret_READ_FAIL;
      }
      char* chr_code_str = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*chr_code_str)) {
	continue;
      }
      char* chr_code_end = token_endnn(chr_code_str);
      char* loadbuf_iter = skip_initial_spaces(chr_code_end);
      if (variant_skip_ct) {
	*chr_code_end = '\0';
	const uint32_t chr_code = get_chr_code(chr_code_str, cip, (uintptr_t)(chr_code_end - chr_code_str));
	if (!is_set(cip->chr_mask, chr_code)) {
	  continue;
	}
      }
      loadbuf_iter = next_token_mult(loadbuf_iter, 4);
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      dosage_t* dosage_vals_iter = dosage_vals;
      while (1) {
	if (widx >= sample_ctl2_m1) {
	  if (widx > sample_ctl2_m1) {
	    break;
	  }
	  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	}
	uintptr_t genovec_word = 0;
	uint32_t dosage_present_hw = 0;
	for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	  loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	  const char cc = *loadbuf_iter;
	  if (is_eoln_kns(cc)) {
	    goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	  }
	  char cc2 = loadbuf_iter[1];
	  if ((cc2 == ' ') || (cc2 == '\t')) {
	    cc2 = loadbuf_iter[3];
	    if (((cc2 == ' ') || (cc2 == '\t')) && ((unsigned char)(loadbuf_iter[5]) <= 32)) {
	      const uint32_t uii = ((uint32_t)((unsigned char)cc)) - 48;
	      const uint32_t ujj = ((uint32_t)((unsigned char)loadbuf_iter[2])) - 48;
	      const uint32_t ukk = ((uint32_t)((unsigned char)loadbuf_iter[4])) - 48;
	      const uint32_t all_or = uii | ujj | ukk;
	      if (all_or < 2) {
		if (!all_or) {
	          genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		  loadbuf_iter = &(loadbuf_iter[5]);
		  continue;
		} else if (uii + ujj + ukk == 1) {
		  uintptr_t cur_geno = ukk * 2 + ujj;
		  genovec_word |= cur_geno << (2 * sample_idx_lowbits);
		  loadbuf_iter = &(loadbuf_iter[5]);
		  continue;
		}
	      }
	    }
	  }
	  double prob_0alt;
	  char* first_dosage_str_end = scanadv_double(loadbuf_iter, &prob_0alt);
	  if (!first_dosage_str_end) {
	    // ignore next two tokens if first token in triplet is not numeric,
	    // since we treat this as missing regardless
	    loadbuf_iter = next_token_mult(loadbuf_iter, 2);
	    if (!loadbuf_iter) {
	      goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	    }
	    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	    loadbuf_iter = token_endnn(loadbuf_iter);
	    continue;
	  }
	  if ((unsigned char)(*first_dosage_str_end) > ' ') {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  loadbuf_iter = skip_initial_spaces(first_dosage_str_end);
	  if (is_eoln_kns(*loadbuf_iter)) {
	    goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	  }
	  double prob_1alt;
	  loadbuf_iter = scanadv_double(loadbuf_iter, &prob_1alt);
	  if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	  if (is_eoln_kns(*loadbuf_iter)) {
	    goto ox_gen_to_pgen_ret_MISSING_TOKENS;
	  }
	  double prob_2alt;
	  loadbuf_iter = scanadv_double(loadbuf_iter, &prob_2alt);
	  if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  // bugfix
	  prob_0alt *= 32768;
	  prob_1alt *= 32768;
	  prob_2alt *= 32768;
	  
	  if ((prob_0alt < 0.0) || (prob_0alt >= 65535.4999999999) || (prob_1alt < 0.0) || (prob_1alt >= 65535.4999999999) || (prob_2alt < 0.0) || (prob_2alt >= 65535.4999999999)) {
	    goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
	  }
	  const uint32_t dosage_int0 = (int32_t)(prob_0alt + 0.5);
	  const uint32_t dosage_int1 = (int32_t)(prob_1alt + 0.5);
	  const uint32_t dosage_int2 = (int32_t)(prob_2alt + 0.5);
	  bgen11_dosage_import_update(dosage_int_sum_thresh, import_dosage_certainty_int, hard_call_halfdist, dosage_erase_halfdist, sample_idx_lowbits, dosage_int0, dosage_int1, dosage_int2, &genovec_word, &dosage_present_hw, &dosage_vals_iter);
	}
	genovec[widx] = genovec_word;
	((halfword_t*)dosage_present)[widx] = (halfword_t)dosage_present_hw;
	++widx;
      }
      if (prov_ref_allele_second) {
	genovec_invert_unsafe(sample_ct, genovec);
	zero_trailing_quaters(sample_ct, genovec);
      }
      if (dosage_vals_iter != dosage_vals) {
	const uint32_t dosage_ct = (uintptr_t)(dosage_vals_iter - dosage_vals);
	if (prov_ref_allele_second) {
	  biallelic_dosage16_invert(dosage_ct, dosage_vals);
	}
	if (spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, &spgw)) {
	  goto ox_gen_to_pgen_ret_WRITE_FAIL;
	}
      } else {
	if (spgw_append_biallelic_genovec(genovec, &spgw)) {
	  goto ox_gen_to_pgen_ret_WRITE_FAIL;
	}
      }
      ++vidx;
      if (!(vidx % 1000)) {
	printf("\r--data/--gen: %uk variants converted.", vidx / 1000);
	fflush(stdout);
      }
    }
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--data/--gen: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar");
    strcpy(write_iter, " written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  ox_gen_to_pgen_ret_LONG_LINE_N:
    putc_unlocked('\n', stdout);
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file is pathologically long.\n", line_idx);
      reterr = kPglRetMalformedInput;
      break;
    }
  ox_gen_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_gen_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_gen_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_gen_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ox_gen_to_pgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ox_gen_to_pgen_ret_INVALID_DOSAGE:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file has an invalid dosage value.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ox_gen_to_pgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ox_gen_to_pgen_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  ox_gen_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_gen_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 ox_gen_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  gzclose_cond(gz_infile);
  fclose_cond(pvarfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

// a few multithread globals
static uint16_t** g_bgen_geno_bufs = nullptr; // per-thread


// per-variant (could make compressed_geno_starts per-thread)
static unsigned char** g_compressed_geno_starts[2] = {nullptr, nullptr};
static uintptr_t* g_write_genovecs[2] = {nullptr, nullptr};
static uint32_t* g_write_dosage_cts[2] = {nullptr, nullptr};
static uintptr_t* g_write_dosage_presents[2] = {nullptr, nullptr};
static dosage_t* g_write_dosage_val_bufs[2] = {nullptr, nullptr};

static uint32_t g_sample_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_cur_block_write_ct = 0;
static uint32_t g_hard_call_halfdist = 0;
static uint32_t g_dosage_erase_halfdist = 0;
static uint32_t g_import_dosage_certainty_int = 0;
static uint32_t g_compression_mode = 0;
static uint32_t g_dosage_is_present = 0;
static uint32_t g_prov_ref_allele_second = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

// static uint32_t* g_error_vidxs = nullptr; // per-thread

THREAD_FUNC_DECL bgen11_dosage_scan_thread(void* arg) {
  // this bails as soon as a single non-hardcall is detected.  still
  // multithreaded due to low speed of uncompress() call, practical value of
  // handling the all-hardcall case efficiently, and reduced code complexity
  // (locally more complex, but globally cleaner due to overlap with
  // bgen11_geno_to_pgen_thread()).
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  uint16_t* bgen_geno_buf = g_bgen_geno_bufs[tidx];
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  // hard_call_halfdist irrelevant here
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t import_dosage_certainty_int = g_import_dosage_certainty_int;
  const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
  const uint32_t compression_mode = g_compression_mode;
  // uint32_t vidx_base = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    unsigned char* compressed_geno_iter = g_compressed_geno_starts[parity][vidx];
    uint16_t* bgen_probs = bgen_geno_buf;
    for (; vidx < vidx_end; ++vidx) {
      if (compression_mode) {
	uint32_t compressed_block_byte_ct;
	memcpy(&compressed_block_byte_ct, compressed_geno_iter, 4);
	compressed_geno_iter = &(compressed_geno_iter[4]);
	uLongf zlib_ulongf = 6 * sample_ct;
	if (uncompress((Bytef*)bgen_probs, &zlib_ulongf, (Bytef*)compressed_geno_iter, compressed_block_byte_ct) != Z_OK) {
	  break;
	}
	compressed_geno_iter = &(compressed_geno_iter[compressed_block_byte_ct]);
      } else {
	bgen_probs = (uint16_t*)compressed_geno_iter;
	compressed_geno_iter = &(compressed_geno_iter[6 * sample_ct]);
      }
      const uint16_t* bgen_probs_iter = bgen_probs;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	const uint32_t dosage_int0 = (uint32_t)(*bgen_probs_iter++);
	const uint32_t dosage_int1 = (uint32_t)(*bgen_probs_iter++);
	const uint32_t dosage_int2 = (uint32_t)(*bgen_probs_iter++);
	const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
	if ((dosage_int_sum > dosage_int_sum_thresh) || (dosage_int0 >= import_dosage_certainty_int) || (dosage_int1 >= import_dosage_certainty_int) || (dosage_int2 >= import_dosage_certainty_int)) {
	  // ties realistically happen, use banker's rounding
	  // 1/65536 -> 0/32768
	  // 3/65536 -> 2/32768
	  // 5/65536 -> 2/32768
	  const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
	  uint32_t write_dosage_int;
	  if (dosage_int_sum == kDosageMax) {
	    // optimize common case
	    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
	  } else {
	    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
	    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
	  }
	  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
	  if (halfdist < dosage_erase_halfdist) {
	    goto bgen11_dosage_scan_thread_dosage_found;
	  }
	}
      }
    }
    if (vidx != vidx_end) {
      // g_error_vidxs[tidx] = vidx + vidx_base;
      g_error_ret = kPglRetMalformedInput;
    }
    while (0) {
    bgen11_dosage_scan_thread_dosage_found:
      g_dosage_is_present = 1;
    }
    // vidx_base += cur_block_write_ct;
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static_assert(sizeof(dosage_t) == 2, "bgen11_geno_to_pgen_thread() needs to be updated.");
THREAD_FUNC_DECL bgen11_geno_to_pgen_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t sample_ct = g_sample_ct;
  uint16_t* bgen_geno_buf = g_bgen_geno_bufs[tidx];
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t import_dosage_certainty_int = g_import_dosage_certainty_int;
  const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
  const uint32_t compression_mode = g_compression_mode;
  const uint32_t prov_ref_allele_second = g_prov_ref_allele_second;
  const uintptr_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  // uint32_t vidx_base = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* compressed_geno_iter = g_compressed_geno_starts[parity][vidx];
    uintptr_t* write_genovec_iter = &(g_write_genovecs[parity][vidx * sample_ctaw2]);
    uint32_t* write_dosage_ct_iter = &(g_write_dosage_cts[parity][vidx]);
    uintptr_t* write_dosage_present_iter = &(g_write_dosage_presents[parity][vidx * sample_ctaw]);
    dosage_t* write_dosage_vals_iter = &(g_write_dosage_val_bufs[parity][vidx * sample_ct]);
    uint16_t* bgen_probs = bgen_geno_buf;
    for (; vidx < vidx_end; ++vidx) {
      if (compression_mode) {
	uint32_t compressed_block_byte_ct;
	memcpy(&compressed_block_byte_ct, compressed_geno_iter, 4);
	compressed_geno_iter = &(compressed_geno_iter[4]);
	uLongf zlib_ulongf = 6 * sample_ct;
	if (uncompress((Bytef*)bgen_probs, &zlib_ulongf, (Bytef*)compressed_geno_iter, compressed_block_byte_ct) != Z_OK) {
	  break;
	}
	compressed_geno_iter = &(compressed_geno_iter[compressed_block_byte_ct]);
      } else {
	bgen_probs = (uint16_t*)compressed_geno_iter;
	compressed_geno_iter = &(compressed_geno_iter[6 * sample_ct]);
      }
      const uint16_t* bgen_probs_iter = bgen_probs;
      dosage_t* cur_dosage_vals_iter = write_dosage_vals_iter;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      while (1) {
	if (widx >= sample_ctl2_m1) {
	  if (widx > sample_ctl2_m1) {
	    break;
	  }
	  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	}
	uintptr_t genovec_word = 0;
	uint32_t dosage_present_hw = 0;
	for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	  const uint32_t dosage_int0 = (uint32_t)(*bgen_probs_iter++);
	  const uint32_t dosage_int1 = (uint32_t)(*bgen_probs_iter++);
	  const uint32_t dosage_int2 = (uint32_t)(*bgen_probs_iter++);
	  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
	  if (dosage_int_sum <= dosage_int_sum_thresh) {
	    if ((dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	      continue;
	    }
	  }
	  const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
	  uint32_t write_dosage_int;
	  if (dosage_int_sum == kDosageMax) {
	    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
	  } else {
	    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
	    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
	  }
	  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
	  if (halfdist < hard_call_halfdist) {
	    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	  } else {
	    genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
	    if (halfdist >= dosage_erase_halfdist) {
	      continue;
	    }
	  }
	  dosage_present_hw |= 1U << sample_idx_lowbits;
	  *cur_dosage_vals_iter++ = write_dosage_int;
	}
	write_genovec_iter[widx] = genovec_word;
	((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
	++widx;
      }
      const uint32_t dosage_ct = (uintptr_t)(cur_dosage_vals_iter - write_dosage_vals_iter);
      if (prov_ref_allele_second) {
	genovec_invert_unsafe(sample_ct, write_genovec_iter);
	zero_trailing_quaters(sample_ct, write_genovec_iter);
	if (dosage_ct) {
	  biallelic_dosage16_invert(dosage_ct, write_dosage_vals_iter);
	}
      }
      *write_dosage_ct_iter++ = dosage_ct;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
      write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
    }
    if (vidx != vidx_end) {
      // g_error_vidxs[tidx] = vidx + vidx_base;
      g_error_ret = kPglRetMalformedInput;
    }
    // vidx_base += cur_block_write_ct;
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static unsigned char** g_thread_wkspaces = nullptr;
static uint32_t* g_thread_bidxs[2] = {nullptr, nullptr};
static uint16_t* g_bgen_allele_cts[2] = {nullptr, nullptr};
static uint32_t* g_uncompressed_genodata_byte_cts[2] = {nullptr, nullptr};

// for each bit precision level, how large must
//   max(numerators, 2^{bit_precision} - 1 - [sum of numerators])
// be to avoid throwing out the genotype?
static uint32_t* g_bgen_import_dosage_certainty_thresholds = nullptr;

// Reliably fast division by constants of the form 2^n - 1; see
//   http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html
// The general case also requires a preshift parameter, but it's always zero
// for the odd .bgen denominators.
typedef struct bgen_magic_num_struct {
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
} bgen_magic_num_t;

static const bgen_magic_num_t kBgenMagicNums[25] = {
  {0, 0, 0},
  {1, 0, 0},
  {2863311531U, 33, 0},
  {1227133513U, 33, 1},
  {2290649225U, 35, 0},
  {1108378657U, 35, 1},
  {1090785345U, 36, 1},
  {270549121, 35, 1},
  {2155905153U, 39, 0},
  {134480385, 36, 1},
  {1074791425U, 40, 1},
  {4196353, 33, 1},
  {16781313, 36, 1},
  {67117057, 39, 1},
  {268451841, 42, 1},
  {1073774593U, 45, 1},
  {2147516417U, 47, 0}
  // todo: check whether something similar works for 17-32 bit cases
  /*
  ,{131073, 34, 1},
  {262145, 36, 1},
  {524289, 38, 1},
  {1048577, 40, 1},
  {2097153, 42, 1},
  {4194305, 44, 1},
  {8388609, 46, 1},
  {16777217, 48, 1},
  {33554433, 50, 1},
  {67108865, 52, 1},
  {134217729, 54, 1},
  {268435457, 56, 1},
  {536870913, 58, 1},
  {1073741825U, 60, 1},
  {2147483649U, 62, 1},
  {2147483649U, 63, 0}
  */
};

static_assert(sizeof(dosage_t) == 2, "bgen13_dosage_or_phase_scan_thread() needs to be updated.");
THREAD_FUNC_DECL bgen13_dosage_or_phase_scan_thread(void* arg) {
  // This bails as soon as a single phased or dosage call is detected.  We
  // provisionally assume e.g. phased calls are also present when dosages are,
  // and clean up the relevant header bytes when the assumption is untrue.
  // (Well, that's how it'll work after phased dosages are implemented,
  // anyway.)
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t* bgen_import_dosage_certainty_thresholds = g_bgen_import_dosage_certainty_thresholds;
  const uint32_t compression_mode = g_compression_mode;
  const unsigned char* cur_uncompressed_geno = nullptr;
  if (compression_mode) {
    cur_uncompressed_geno = g_thread_wkspaces[tidx];
  }
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;

    // this is just used as an no-error flag
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    if (cur_block_write_ct) {
      const uint32_t bidx_end = g_thread_bidxs[parity][tidx + 1];
      uint32_t bidx = g_thread_bidxs[parity][tidx];
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[parity];
      const uint16_t* bgen_allele_cts = g_bgen_allele_cts[parity];
      const uint32_t* uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[parity];
      for (; bidx < bidx_end; ++bidx) {
	const unsigned char* compressed_geno_start = compressed_geno_starts[bidx];
	const unsigned char* compressed_geno_end = compressed_geno_starts[bidx + 1];
	uint32_t compressed_byte_ct = (uintptr_t)(compressed_geno_end - compressed_geno_start);
	uint32_t uncompressed_byte_ct;
	if (compression_mode) {
	  uncompressed_byte_ct = uncompressed_genodata_byte_cts[bidx];
	  if (compression_mode == 1) {
	    uLongf zlib_ulongf = uncompressed_byte_ct;
	    // const_cast
	    if (uncompress((Bytef*)((uintptr_t)cur_uncompressed_geno), &zlib_ulongf, (const Bytef*)compressed_geno_start, compressed_byte_ct) != Z_OK) {
	      // possible todo: report variant index
	      goto bgen13_dosage_or_phase_scan_thread_malformed;
	    }
	  } else {
	    // const_cast
	    const uintptr_t extracted_byte_ct = ZSTD_decompress((void*)((uintptr_t)cur_uncompressed_geno), uncompressed_byte_ct, compressed_geno_start, compressed_byte_ct);
	    if (extracted_byte_ct != uncompressed_byte_ct) {
	      // possible todo: inspect error code
	      goto bgen13_dosage_or_phase_scan_thread_malformed;
	    }
	  }
	} else {
	  cur_uncompressed_geno = compressed_geno_start;
	  uncompressed_byte_ct = compressed_byte_ct;
	}
	// 4 bytes: sample_ct
	// 2 bytes: # of alleles, must match bgen_allele_cts[bidx]
	// 1 byte: min ploidy
	// 1 byte: max ploidy
	// sample_ct bytes: low 6 bits = ploidy, top bit = missingness
	// 1 byte: 1 if phased, 0 if not
	// 1 byte: # of bits of probability precision (we just support 8 and 16
	//         for now, add others later)
	if ((uncompressed_byte_ct < 10 + sample_ct) || memcmp(&sample_ct, cur_uncompressed_geno, 4)) {
	  goto bgen13_dosage_or_phase_scan_thread_malformed;
	}
	const uint32_t cur_allele_ct = bgen_allele_cts[bidx];
	if (*((const uint16_t*)(&(cur_uncompressed_geno[4]))) != cur_allele_ct) {
	  goto bgen13_dosage_or_phase_scan_thread_malformed;
	}
	const uint32_t min_ploidy = cur_uncompressed_geno[6];
	const uint32_t max_ploidy = cur_uncompressed_geno[7];
	if ((min_ploidy > max_ploidy) || (max_ploidy > 63)) {
	  goto bgen13_dosage_or_phase_scan_thread_malformed;
	}
	if (max_ploidy > 2) {
	  goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
	}
	const unsigned char* missing_and_ploidy_info = &(cur_uncompressed_geno[8]);
	const unsigned char* uncompressed_geno_iter = &(cur_uncompressed_geno[8 + sample_ct]);
	const uint32_t is_phased = *uncompressed_geno_iter++;
	if (is_phased > 1) {
	  goto bgen13_dosage_or_phase_scan_thread_malformed;
	}
	const uint32_t bit_precision = *uncompressed_geno_iter++;
	if ((!bit_precision) || (bit_precision > 32)) {
	  goto bgen13_dosage_or_phase_scan_thread_malformed;
	}
	if (bit_precision > 16) {
	  goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
	}
	const uint64_t totq_magic = kBgenMagicNums[bit_precision].totq_magic;
	const uint32_t totq_postshift = kBgenMagicNums[bit_precision].totq_postshift;
	uint32_t totq_incr = kBgenMagicNums[bit_precision].totq_incr;
	const uint32_t bytes_per_prob = DIV_UP(bit_precision, CHAR_BIT);

	// also equal to denominator
	const uintptr_t numer_mask = (1U << bit_precision) - 1;

        // diploid (haploid is identical except b is always zero):
        //   round((32768a + 16384b)/(2^{bit precision} - 1))
	//   floor((32768a + 16384b)/(2^{bit_precision} - 1) + 0.5)
	// = floor((32768a + 16384b + 2^{bit_precision - 1})
        //     / (2^{bit_precision} - 1))
	// = (totq_magic * (32768a + 16384b + 2^{bits-1} + totq_incr))
        //     >> totq_postshift
	//
	// This works fine for bit_precision <= 16, anyway.  There are two
	// issues which come up with higher precision:
	// 1. The ridiculous_fish magic numbers assume a 32-bit dividend.  Our
	//    dividend is guaranteed to be divisible by 2^14, but it can be as
	//    large as
	//      (2^{bits} - 1) * 2^15 + 2^{bits-1}.
	//    I would not be surprised if a similar approach still works with
	//    bits > 16, but I'm pretty sure the magic-number-generating
	//    function would need to be different.
	// 2. Relatedly, the current sequence of operations multiplies
	//    totq_magic by (dividend + totq_incr) (where totq_incr is zero or
	//    one); this intermediate result must not overflow a uint64_t.
        //
        // Meanwhile, idempotence is not possible for --import-dosage-certainty
	// anyway, so we apply that check to the pre-conversion numerators.
	totq_incr += 1U << (bit_precision - 1);
	uint32_t numer_certainty_min = 0;
	if (bgen_import_dosage_certainty_thresholds) {
	  numer_certainty_min = bgen_import_dosage_certainty_thresholds[bit_precision];
	}

	if (is_phased) {
	  // todo
	  goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
	} else {
	  if (cur_allele_ct == 2) {
	    if (min_ploidy == max_ploidy) {
	      // faster handling of common cases (no need to keep checking if
	      // we've read past the end)
	      if (uncompressed_byte_ct != (1 + bytes_per_prob * (max_ploidy * k1LU)) * sample_ct + 10) {
		goto bgen13_dosage_or_phase_scan_thread_malformed;
	      }
	      if (max_ploidy == 2) {
		for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
		  const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
		  // treat anything else as missing
		  if (missing_and_ploidy == 2) {
		    const unsigned char* sample_probs_start = &(uncompressed_geno_iter[sample_idx * 2 * bytes_per_prob]);
#ifdef __arm__
  #error "Unaligned accesses in bgen13_dosage_or_phase_scan_thread()."
#endif
		    // this can read slightly past the end of the buffer
		    const uintptr_t numer_aa = (*((const uint32_t*)sample_probs_start)) & numer_mask;
		    const uintptr_t numer_ab = (*((const uint32_t*)(&(sample_probs_start[bytes_per_prob])))) & numer_mask;
		    if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
		      // treat as missing
		      continue;
		    }
		    const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
		    const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		    if (halfdist < dosage_erase_halfdist) {
		      goto bgen13_dosage_or_phase_scan_thread_found;
		    }
		  }
		}
	      } else if (max_ploidy == 1) {
		for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
		  const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
		  if (missing_and_ploidy == 1) {
		    const unsigned char* sample_probs_start = &(uncompressed_geno_iter[sample_idx * bytes_per_prob]);
		    const uintptr_t numer_a = (*((const uint32_t*)sample_probs_start)) & numer_mask;
		    if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
		      continue;
		    }
		    const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
		    const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		    if (halfdist < dosage_erase_halfdist) {
		      goto bgen13_dosage_or_phase_scan_thread_found;
		    }
		  }
		}
	      }
	      // don't need to do anything in all-ploidy-0 case
	    } else {
	      const unsigned char* uncompressed_geno_end = &(cur_uncompressed_geno[uncompressed_byte_ct]);
	      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
		if (uncompressed_geno_iter > uncompressed_geno_end) {
		  goto bgen13_dosage_or_phase_scan_thread_malformed;
		}
		uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
		if (missing_and_ploidy == 2) {
		  const uintptr_t numer_aa = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
		  const uintptr_t numer_ab = (*((const uint32_t*)(&(uncompressed_geno_iter[bytes_per_prob])))) & numer_mask;
		  uncompressed_geno_iter = &(uncompressed_geno_iter[2 * bytes_per_prob]);
		  if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
		    // treat as missing
		    continue;
		  }
		  const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
		  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		  if (halfdist < dosage_erase_halfdist) {
		    goto bgen13_dosage_or_phase_scan_thread_found;
		  }
		} else if (missing_and_ploidy == 1) {
		  const uintptr_t numer_a = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
		  uncompressed_geno_iter = &(uncompressed_geno_iter[bytes_per_prob]);
		  if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
		    continue;
		  }
		  const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
		  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		  if (halfdist < dosage_erase_halfdist) {
		    goto bgen13_dosage_or_phase_scan_thread_found;
		  }
		} else {
		  // treat as missing
		  missing_and_ploidy &= 127;
		  if (missing_and_ploidy > 2) {
		    goto bgen13_dosage_or_phase_scan_thread_malformed;
		  }
		  uncompressed_geno_iter = &(uncompressed_geno_iter[missing_and_ploidy * bytes_per_prob]);
		}
	      }
	    }
	  } else {
	    // todo: unphased multiallelic variants
	    // (shouldn't currently be possible to reach here, I/O thread skips
	    // multiallelics for now)
	    assert(0);
	    goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
	  }
	}
      }
    }
    while (0) {
    bgen13_dosage_or_phase_scan_thread_malformed:
      g_error_ret = kPglRetMalformedInput;
      break;
    bgen13_dosage_or_phase_scan_thread_not_yet_supported:
      g_error_ret = kPglRetNotYetSupported;
      break;
    bgen13_dosage_or_phase_scan_thread_found:
      g_dosage_is_present = 1;
      break;
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static uintptr_t* g_write_phasepresents[2] = {nullptr, nullptr};
static uintptr_t* g_write_phaseinfos[2] = {nullptr, nullptr};
static uintptr_t* g_write_dphase_presents[2] = {nullptr, nullptr};
static uint32_t* g_write_dphase_cts[2] = {nullptr, nullptr};

static_assert(sizeof(dosage_t) == 2, "bgen13_geno_to_pgen_thread() needs to be updated.");
THREAD_FUNC_DECL bgen13_geno_to_pgen_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t sample_ct = g_sample_ct;
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t* bgen_import_dosage_certainty_thresholds = g_bgen_import_dosage_certainty_thresholds;
  const uint32_t compression_mode = g_compression_mode;
  const uint32_t prov_ref_allele_second = g_prov_ref_allele_second;
  const uintptr_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  const unsigned char* cur_uncompressed_geno = nullptr;
  if (compression_mode) {
    cur_uncompressed_geno = g_thread_wkspaces[tidx];
  }
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;

    // this is just used as an no-error flag
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    if (cur_block_write_ct) {
      const uint32_t bidx_end = g_thread_bidxs[parity][tidx + 1];
      uint32_t bidx = g_thread_bidxs[parity][tidx];
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[parity];
      uintptr_t* write_genovec_iter = &(g_write_genovecs[parity][bidx * sample_ctaw2]);
      uint32_t* write_dosage_ct_iter = &(g_write_dosage_cts[parity][bidx]);
      uint32_t* write_dphase_ct_iter = &(g_write_dphase_cts[parity][bidx]);
      uintptr_t* write_dosage_present_iter = &(g_write_dosage_presents[parity][bidx * sample_ctaw]);
      uintptr_t* write_dphase_present_iter = &(g_write_dphase_presents[parity][bidx * sample_ctaw]);
      dosage_t* write_dosage_vals_iter = &(g_write_dosage_val_bufs[parity][bidx * sample_ct * 2]);
      const uint16_t* bgen_allele_ct_iter = &(g_bgen_allele_cts[parity][bidx]);
      const uint32_t* uncompressed_genodata_byte_ct_iter = &(g_uncompressed_genodata_byte_cts[parity][bidx]);
      for (; bidx < bidx_end; ++bidx) {
	const unsigned char* compressed_geno_start = compressed_geno_starts[bidx];
	const unsigned char* compressed_geno_end = compressed_geno_starts[bidx + 1];
	uint32_t compressed_byte_ct = (uintptr_t)(compressed_geno_end - compressed_geno_start);
	uint32_t uncompressed_byte_ct;
	if (compression_mode) {
	  uncompressed_byte_ct = *uncompressed_genodata_byte_ct_iter++;
	  if (compression_mode == 1) {
	    uLongf zlib_ulongf = uncompressed_byte_ct;
	    // const_cast
	    if (uncompress((Bytef*)((uintptr_t)cur_uncompressed_geno), &zlib_ulongf, (const Bytef*)compressed_geno_start, compressed_byte_ct) != Z_OK) {
	      // possible todo: report variant index
	      goto bgen13_geno_to_pgen_thread_malformed;
	    }
	  } else {
            // const_cast
	    const uintptr_t extracted_byte_ct = ZSTD_decompress((void*)((uintptr_t)cur_uncompressed_geno), uncompressed_byte_ct, compressed_geno_start, compressed_byte_ct);
	    if (extracted_byte_ct != uncompressed_byte_ct) {
	      // possible todo: inspect error code
	      goto bgen13_geno_to_pgen_thread_malformed;
	    }
	  }
	} else {
	  cur_uncompressed_geno = compressed_geno_start;
	  uncompressed_byte_ct = compressed_byte_ct;
	}
	// 4 bytes: sample_ct
	// 2 bytes: # of alleles, must match bgen_allele_cts[bidx]
	// 1 byte: min ploidy
	// 1 byte: max ploidy
	// sample_ct bytes: low 6 bits = ploidy, top bit = missingness
	// 1 byte: 1 if phased, 0 if not
	// 1 byte: # of bits of probability precision (we just support 8 and 16
	//         for now, add others later)
	if ((uncompressed_byte_ct < 10 + sample_ct) || memcmp(&sample_ct, cur_uncompressed_geno, 4)) {
	  goto bgen13_geno_to_pgen_thread_malformed;
	}
	const uint32_t cur_allele_ct = *bgen_allele_ct_iter++;
	if (*((const uint16_t*)(&(cur_uncompressed_geno[4]))) != cur_allele_ct) {
	  goto bgen13_geno_to_pgen_thread_malformed;
	}
	const uint32_t min_ploidy = cur_uncompressed_geno[6];
	const uint32_t max_ploidy = cur_uncompressed_geno[7];
	if ((min_ploidy > max_ploidy) || (max_ploidy > 63)) {
	  goto bgen13_geno_to_pgen_thread_malformed;
	}
	if (max_ploidy > 2) {
	  goto bgen13_geno_to_pgen_thread_not_yet_supported;
	}
	const unsigned char* missing_and_ploidy_iter = &(cur_uncompressed_geno[8]);
	const unsigned char* uncompressed_geno_iter = &(cur_uncompressed_geno[8 + sample_ct]);
	const uint32_t is_phased = *uncompressed_geno_iter++;
	if (is_phased > 1) {
	  goto bgen13_geno_to_pgen_thread_malformed;
	}
	const uint32_t bit_precision = *uncompressed_geno_iter++;
	if ((!bit_precision) || (bit_precision > 32)) {
	  goto bgen13_geno_to_pgen_thread_malformed;
	}
	if (bit_precision > 16) {
	  goto bgen13_geno_to_pgen_thread_not_yet_supported;
	}
	const uint64_t totq_magic = kBgenMagicNums[bit_precision].totq_magic;
	const uint32_t totq_postshift = kBgenMagicNums[bit_precision].totq_postshift;
	uint32_t totq_incr = kBgenMagicNums[bit_precision].totq_incr;
	const uint32_t bytes_per_prob = DIV_UP(bit_precision, CHAR_BIT);

	// also equal to denominator
	const uintptr_t numer_mask = (1U << bit_precision) - 1;

	totq_incr += 1U << (bit_precision - 1);
	uint32_t numer_certainty_min = 0;
	if (bgen_import_dosage_certainty_thresholds) {
	  numer_certainty_min = bgen_import_dosage_certainty_thresholds[bit_precision];
	}

	dosage_t* cur_dosage_vals_iter = write_dosage_vals_iter;
	uint32_t inner_loop_last = kBitsPerWordD2 - 1;
	uint32_t widx = 0;
	if (is_phased) {
	  // todo
	  goto bgen13_geno_to_pgen_thread_not_yet_supported;
	} else {
	  // fill_ulong_zero(sample_ctaw, write_dphase_present_iter);
	  if (cur_allele_ct == 2) {
	    if (min_ploidy == max_ploidy) {
	      // faster handling of common cases (no need to keep checking if
	      // we've read past the end)
	      if (uncompressed_byte_ct != (bytes_per_prob * (max_ploidy * k1LU) + 1) * sample_ct + 10) {
		goto bgen13_geno_to_pgen_thread_malformed;
	      }
	      if (max_ploidy == 2) {
		while (1) {
		  if (widx >= sample_ctl2_m1) {
		    if (widx > sample_ctl2_m1) {
		      break;
		    }
		    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
		  }
		  uintptr_t genovec_word = 0;
		  uint32_t dosage_present_hw = 0;
		  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits, uncompressed_geno_iter = &(uncompressed_geno_iter[2 * bytes_per_prob])) {
		    const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
		    if (missing_and_ploidy == 2) {
#ifdef __arm__
  #error "Unaligned accesses in bgen13_geno_to_pgen_thread()."
#endif
		      const uintptr_t numer_aa = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
		      const uintptr_t numer_ab = (*((const uint32_t*)(&(uncompressed_geno_iter[bytes_per_prob])))) & numer_mask;
		      if (numer_aa + numer_ab > numer_mask) {
			goto bgen13_geno_to_pgen_thread_malformed;
		      }
		      if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
			// missing due to --import-dosage-certainty
			goto bgen13_geno_to_pgen_thread_diploid_missing;
		      }
		      const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
		      const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		      if (halfdist < hard_call_halfdist) {
			genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		      } else {
			genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
			if (halfdist >= dosage_erase_halfdist) {
			  continue;
			}
		      }
		      dosage_present_hw |= 1U << sample_idx_lowbits;
		      *cur_dosage_vals_iter++ = write_dosage_int;
		    } else {
		      // (could also validate that missing_and_ploidy == 130)
		    bgen13_geno_to_pgen_thread_diploid_missing:
		      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		    }
		  }
		  write_genovec_iter[widx] = genovec_word;
		  ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
		  ++widx;
		}
	      } else if (max_ploidy == 1) {
		while (1) {
		  if (widx >= sample_ctl2_m1) {
		    if (widx > sample_ctl2_m1) {
		      break;
		    }
		    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
		  }
		  uintptr_t genovec_word = 0;
		  uint32_t dosage_present_hw = 0;
		  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits, uncompressed_geno_iter = &(uncompressed_geno_iter[bytes_per_prob])) {
		    const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
		    if (missing_and_ploidy == 1) {
		      const uintptr_t numer_a = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
		      if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
			goto bgen13_geno_to_pgen_thread_haploid_missing;
		      }
		      const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
		      const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		      if (halfdist < hard_call_halfdist) {
			genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		      } else {
			genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
			if (halfdist >= dosage_erase_halfdist) {
			  continue;
			}
		      }
		      dosage_present_hw |= 1U << sample_idx_lowbits;
		      *cur_dosage_vals_iter++ = write_dosage_int;
		    } else {
		    bgen13_geno_to_pgen_thread_haploid_missing:
		      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		    }
		  }
		  write_genovec_iter[widx] = genovec_word;
		  ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
		  ++widx;
		}
	      }
	      // don't need to do anything in all-ploidy-0 case
	    } else {
	      const unsigned char* uncompressed_geno_end = &(cur_uncompressed_geno[uncompressed_byte_ct]);
	      while (1) {
		if (widx >= sample_ctl2_m1) {
		  if (widx > sample_ctl2_m1) {
		    break;
		  }
		  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
		}
		uintptr_t genovec_word = 0;
		uint32_t dosage_present_hw = 0;
		for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		  if (uncompressed_geno_iter > uncompressed_geno_end) {
		    goto bgen13_geno_to_pgen_thread_malformed;
		  }
		  uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
		  uint32_t write_dosage_int;
		  if (missing_and_ploidy == 2) {
		    const uintptr_t numer_aa = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
		    const uintptr_t numer_ab = (*((const uint32_t*)(&(uncompressed_geno_iter[bytes_per_prob])))) & numer_mask;
		    uncompressed_geno_iter = &(uncompressed_geno_iter[2 * bytes_per_prob]);
		    if (numer_aa + numer_ab > numer_mask) {
		      goto bgen13_geno_to_pgen_thread_malformed;
		    }
		    if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
		      // missing due to --import-dosage-certainty
		      goto bgen13_geno_to_pgen_thread_generic_missing;
		    }
		    write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
		  } else if (missing_and_ploidy == 1) {
		    const uintptr_t numer_a = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
		    uncompressed_geno_iter = &(uncompressed_geno_iter[bytes_per_prob]);
		    if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
		      goto bgen13_geno_to_pgen_thread_generic_missing;
		    }
		    write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
		  } else {
		    missing_and_ploidy &= 127;
		    if (missing_and_ploidy > 2) {
		      goto bgen13_geno_to_pgen_thread_malformed;
		    }
		    uncompressed_geno_iter = &(uncompressed_geno_iter[missing_and_ploidy * bytes_per_prob]);
		  bgen13_geno_to_pgen_thread_generic_missing:
		    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		    continue;
		  }
		  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
		  if (halfdist < hard_call_halfdist) {
		    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		  } else {
		    genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
		    if (halfdist >= dosage_erase_halfdist) {
		      continue;
		    }
		  }
		  dosage_present_hw |= 1U << sample_idx_lowbits;
		  *cur_dosage_vals_iter++ = write_dosage_int;
		}
		write_genovec_iter[widx] = genovec_word;
		((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
		++widx;
	      }
	    }
	    const uint32_t dosage_ct = (uintptr_t)(cur_dosage_vals_iter - write_dosage_vals_iter);
	    // note that this is inverted from bgen-1.1
	    if (!prov_ref_allele_second) {
	      genovec_invert_unsafe(sample_ct, write_genovec_iter);
	      zero_trailing_quaters(sample_ct, write_genovec_iter);
	      if (dosage_ct) {
		biallelic_dosage16_invert(dosage_ct, write_dosage_vals_iter);
	      }
	    }
	    *write_dosage_ct_iter++ = dosage_ct;
	    *write_dphase_ct_iter++ = 0;
	    write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
	    write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
	    write_dphase_present_iter = &(write_dphase_present_iter[sample_ctaw]);
	    write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct * 2]);
	  } else {
	    // todo: unphased multiallelic variants
	    assert(0);
	    goto bgen13_geno_to_pgen_thread_not_yet_supported;
	  }
	}
      }
    }
    while (0) {
    bgen13_geno_to_pgen_thread_malformed:
      g_error_ret = kPglRetMalformedInput;
      break;
    bgen13_geno_to_pgen_thread_not_yet_supported:
      g_error_ret = kPglRetNotYetSupported;
      break;
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static_assert(sizeof(dosage_t) == 2, "ox_bgen_to_pgen() needs to be updated.");
pglerr_t ox_bgen_to_pgen(const char* bgenname, const char* samplename, const char* const_fid, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* bgenfile = nullptr;

  // only if no sample file specified, and .bgen has sample IDs.  (possible
  // todo: consistency check when both sources of sample IDs are present?)
  FILE* psamfile = nullptr;
  
  FILE* pvarfile = nullptr;
  threads_state_t ts;
  init_threads3z(&ts);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    // Pass 1: Determine whether there's at least one non-hardcall needs to be
    //         saved, and if a chromosome filter was specified, count the
    //         number of variants which pass the filter.
    //         For bgen-1.2/1.3, the .pvar is also written in this pass.
    //         For bgen-1.1, we can usually early-bail when no chromosome
    //         filter is involved, so .pvar writing is postponed till the
    //         second pass.
    // Pass 2: Write .pgen file.

    if (fopen_checked(bgenname, FOPEN_RB, &bgenfile)) {
      goto ox_bgen_to_pgen_ret_OPEN_FAIL;
    }
    uint32_t initial_uints[5];
    if (!fread(initial_uints, 20, 1, bgenfile)) {
      // this could be malformed input as well; could distinguish later?
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    if (initial_uints[1] > initial_uints[0]) {
      logerrprint("Error: Invalid .bgen header.\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    const uint32_t raw_variant_ct = initial_uints[2];
    if (!raw_variant_ct) {
      // permit this in --allow-no-vars case?
      logerrprint("Error: Empty .bgen file.\n");
      goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
    }
    const uint32_t sample_ct = initial_uints[3];
    if (initial_uints[4] && (initial_uints[4] != 0x6e656762)) {
      logerrprint("Error: Invalid .bgen magic number.\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }

    if (fseeko(bgenfile, initial_uints[1], SEEK_SET)) {
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    uint32_t header_flags;
    if (!fread(&header_flags, 4, 1, bgenfile)) {
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    const uint32_t compression_mode = header_flags & 3;
    const uint32_t layout = (header_flags >> 2) & 15;
    if (!layout) {
      logerrprint("Error: BGEN v1.0 files are not supported by " PROG_NAME_STR ".\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    if ((compression_mode == 3) || (layout > 2)) {
      logerrprint("Error: Unrecognized BGEN version.  Use gen-convert or a similar tool to\ndowncode to BGEN v1.3 if you want to process this data with " PROG_NAME_STR ".\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    if ((compression_mode == 2) && (layout == 1)) {
      logerrprint("Error: Invalid .bgen header.\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    LOGPRINTF("--bgen: %u variant%s detected, format v1.%c.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s", (layout == 1)? '1' : ((compression_mode == 2)? '3' : '2'));
    if (samplename[0]) {
      uint32_t sfile_sample_ct;
      reterr = ox_sample_to_psam(samplename, ox_missing_code, misc_flags, outname, outname_end, &sfile_sample_ct);
      if (reterr) {
	goto ox_bgen_to_pgen_ret_1;
      }
      if (sfile_sample_ct != sample_ct) {
	LOGERRPRINTF("Error: .sample file has %u sample%s, while .bgen file has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", sample_ct);
	goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      if (header_flags >> 31) {
	uint32_t sample_id_block_byte_ct;
	uint32_t sample_id_block_entry_ct;
	if ((!fread(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
	    (!fread(&sample_id_block_entry_ct, 4, 1, bgenfile))) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if ((((uint64_t)sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
	    (sample_id_block_entry_ct != sample_ct)) {
	  logerrprint("Error: Invalid .bgen header.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
      }
    } else {
      if (!(header_flags >> 31)) {
	logerrprint("Error: .bgen file does not contain sample IDs, and no .sample file was\nspecified.\n");
	goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      // possible todo: optionally error out if sample IDs aren't consistent
      // between .bgen and .sample

      // see vcf_sample_line()
      // probable todo: wrap much of this in its own function
      uint32_t double_id = (misc_flags / kfMiscDoubleId) & 1;
      uintptr_t const_fid_len = 0;
      if (const_fid) {
	const_fid_len = strlen(const_fid);
      } else if ((!double_id) && (!id_delim)) {
	// default: --double-id + --id-delim
	double_id = 1;
	id_delim = '_';
      }
      const uint32_t double_or_const_fid = double_id || const_fid;
      uint32_t sample_id_block_byte_ct;
      uint32_t sample_id_block_entry_ct;
      if ((!fread(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
	  (!fread(&sample_id_block_entry_ct, 4, 1, bgenfile))) {
	goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      if ((sample_id_block_byte_ct < 8) ||
	  (((uint64_t)sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
	  (sample_id_block_entry_ct != sample_ct)) {
	logerrprint("Error: Invalid .bgen header.\n");
	goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
      }
      sample_id_block_byte_ct -= 8;
      unsigned char* sample_id_block_main = bigstack_alloc(sample_id_block_byte_ct);
      if (!sample_id_block_main) {
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      unsigned char* sample_id_block_end = &(sample_id_block_main[sample_id_block_byte_ct]);
      if (fread_checked(sample_id_block_main, sample_id_block_byte_ct, bgenfile)) {
	goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      
      // high 16 bits always zero
      // we don't just use a uint16_t since we add 2
      uint32_t input_id_slen = 0;
      
      // always check if any tab/eoln characters are present, and error out if
      // so
      // if id_delim != ' ', also check if spaces are present; if so, replace
      // with --idspace-to character or error out
      unsigned char* sample_id_block_iter = sample_id_block_main;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	memcpy(&input_id_slen, sample_id_block_iter, 2);

	// need to check this to avoid read-past-the-end indeterminate
	// behavior
	if ((uintptr_t)(sample_id_block_end - sample_id_block_iter) < input_id_slen + 2) {
	  logerrprint("Error: Invalid .bgen header.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	unsigned char* sample_id_iter = &(sample_id_block_iter[2]);
	unsigned char* sample_id_end = &(sample_id_iter[input_id_slen]);
        uint32_t char_code_min = 32 + (id_delim != ' ');
	for (; sample_id_iter != sample_id_end; ++sample_id_iter) {
	  const uint32_t char_code = *sample_id_iter;
	  if (char_code < char_code_min) {
	    if (char_code < 32) {
	      logerrprint("Error: .bgen sample ID contains tabs, newlines, and/or nonprinting characters.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (!idspace_to) {
	      logerrprint("Error: .bgen sample ID contains space(s).  Use --idspace-to to convert them to\nanother character, or \"--id-delim ' '\" to interpret the spaces as FID/IID\ndelimiters.\n");
	      goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	    }
	    *sample_id_iter = idspace_to;
	  }
	}
	sample_id_block_iter = sample_id_end;
      }
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, FOPEN_WB, &psamfile)) {
	goto ox_bgen_to_pgen_ret_OPEN_FAIL;
      }
      char* textbuf = g_textbuf;
      char* write_iter = strcpya(textbuf, "#FID\tIID");
      uint32_t sid_present = 0;
      if (id_delim) {
	// check if three-part IDs are present
	sample_id_block_iter = sample_id_block_main;
	for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	  memcpy(&input_id_slen, sample_id_block_iter, 2);
	  if ((uintptr_t)(sample_id_block_end - sample_id_block_iter) < input_id_slen + 2) {
	    logerrprint("Error: Invalid .bgen header.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	  unsigned char* sample_id_start = &(sample_id_block_iter[2]);
	  unsigned char* sample_id_end = &(sample_id_start[input_id_slen]);
	  unsigned char* first_delim = (unsigned char*)memchr(sample_id_start, (unsigned char)id_delim, (uintptr_t)(sample_id_end - sample_id_start));
	  if (first_delim) {
	    unsigned char* iid_start = &(first_delim[1]);
	    if (memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(sample_id_end - iid_start)) != nullptr) {
	      sid_present = 1;
	      write_iter = strcpya(write_iter, "\tSID");
	      break;
	    }
	  }
	  sample_id_block_iter = sample_id_end;
	}
      }
      write_iter = strcpya(write_iter, "\tSEX");
      append_binary_eoln(&write_iter);
      char* textbuf_flush = &(textbuf[kMaxMediumLine]);
      sample_id_block_iter = sample_id_block_main;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	memcpy(&input_id_slen, sample_id_block_iter, 2);
	if ((uintptr_t)(sample_id_block_end - sample_id_block_iter) < input_id_slen + 2) {
	  logerrprint("Error: Invalid .bgen header.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	unsigned char* sample_id_start = &(sample_id_block_iter[2]);
	if (input_id_slen <= 1) {
	  if (!input_id_slen) {
	    logerrprint("Error: Empty sample ID in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	  if (*sample_id_start == '0') {
	    logerrprint("Error: Sample ID cannot be '0'.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	}
	unsigned char* sample_id_end = &(sample_id_start[input_id_slen]);
	if (id_delim) {
	  if (*sample_id_start == id_delim) {
	    sprintf(g_logbuf, "Error: '%c' at beginning of sample ID.\n", id_delim);
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_2;
	  }
	  unsigned char* first_delim = (unsigned char*)memchr(sample_id_start, (unsigned char)id_delim, input_id_slen);
	  if (!first_delim) {
	    if (double_or_const_fid) {
	      goto ox_bgen_to_pgen_one_sample_id;
	    }
	    sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_2;
	  }
	  unsigned char* iid_start = &(first_delim[1]);
	  unsigned char* iid_end = (unsigned char*)memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(sample_id_end - iid_start));
	  const unsigned char* sid_start = (const unsigned char*)(&(g_one_char_strs[96]));
	  uint32_t sid_slen = 1;
	  if (iid_end) {
	    if (iid_start == iid_end) {
	      sprintf(g_logbuf, "Error: Consecutive instances of '%c' in sample ID.\n", id_delim);
	      goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_DELIM;
	    }
	    sid_start = &(iid_end[1]);
	    sid_slen = (uintptr_t)(sample_id_end - sid_start);
	    if (memchr(sid_start, (unsigned char)id_delim, sid_slen)) {
	      sprintf(g_logbuf, "Error: More than two instances of '%c' in sample ID.\n", id_delim);
	      goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_DELIM;
	    }
	  } else {
	    iid_end = sample_id_end;
	  }
	  const uint32_t fid_slen = (uintptr_t)(first_delim - sample_id_start);
	  if (fid_slen > kMaxIdSlen) {
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  write_iter = memcpyax(write_iter, sample_id_start, fid_slen, '\t');
	  const uint32_t iid_slen = (uintptr_t)(iid_end - iid_start);
	  if ((*iid_start == '0') && (iid_slen == 1)) {
	    logerrprint("Error: Sample ID induces an invalid IID of '0'.\n");
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	  }
	  if (iid_slen > kMaxIdSlen) {
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  write_iter = memcpya(write_iter, iid_start, iid_slen);
	  if (sid_present) {
	    *write_iter++ = '\t';
	    write_iter = memcpya(write_iter, sid_start, sid_slen);
	  }
	} else {
	ox_bgen_to_pgen_one_sample_id:
	  if (input_id_slen > kMaxIdSlen) {
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID;
	  }
	  if (double_id) {
	    write_iter = memcpya(write_iter, sample_id_start, input_id_slen);
	  } else {
	    write_iter = memcpya(write_iter, const_fid, const_fid_len);
	  }
	  *write_iter++ = '\t';
	  write_iter = memcpya(write_iter, sample_id_start, input_id_slen);
	  if (sid_present) {
	    write_iter = strcpya(write_iter, "\t0");
	  }
	}
	// SEX
	write_iter = memcpyl3a(write_iter, "\tNA");
	append_binary_eoln(&write_iter);
	if (write_iter >= textbuf_flush) {
	  if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), psamfile)) {
	    goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = textbuf;
	}
	sample_id_block_iter = sample_id_end;
      }
      if (sample_id_block_iter != &(sample_id_block_main[sample_id_block_byte_ct])) {
	logerrprint("Error: Invalid .bgen header.\n");
	goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
      }
      if (write_iter != textbuf) {
	if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), psamfile)) {
	  goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	}
      }
      bigstack_reset(sample_id_block_main);
      if (fclose_null(&psamfile)) {
	goto ox_bgen_to_pgen_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("--bgen: %u sample ID%s written to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    }
    if (fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET)) {
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    finalize_chrset(misc_flags, cip);
    const uint32_t autosome_ct_p1 = cip->autosome_ct + 1;
    uint32_t chr_filter_present = (popcount_bit_idx(cip->chr_mask, 0, autosome_ct_p1) != autosome_ct_p1) || (allow_extra_chrs && (cip->is_include_stack || cip->incl_excl_name_stack));
    if (!chr_filter_present) {
      for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
	if (cip->xymt_codes[xymt_idx] >= 0) {
	  if (!is_set(cip->chr_mask, autosome_ct_p1 + xymt_idx)) {
	    chr_filter_present = 1;
	    break;
	  }
	}
      }
    }
    
    char* writebuf = (char*)bigstack_alloc_raw(kMaxMediumLine + kCompressStreamBlock + kCacheline);
    char* writebuf_flush = &(writebuf[kCompressStreamBlock]);
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &pvarfile)) {
      goto ox_bgen_to_pgen_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    write_iter = strcpya(write_iter, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    const uint32_t snpid_chr = (oxford_import_flags & kfOxfordImportBgenSnpIdChr);

    // true for both provisional-reference and real-reference second
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);

    if (hard_call_thresh == 0xffffffffU) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
    const uint32_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
    uint32_t dosage_is_present = 0;
    g_sample_ct = sample_ct;
    g_hard_call_halfdist = kDosage4th - hard_call_thresh;
    g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    g_compression_mode = compression_mode;
    g_prov_ref_allele_second = prov_ref_allele_second;
    g_error_ret = kPglRetSuccess;
    g_dosage_is_present = 0;
    if (layout == 1) {
      // v1.1
      uintptr_t loadbuf_size = round_down_pow2(bigstack_left() / 4, kCacheline);
#ifdef __LP64__
      if (loadbuf_size > kMaxLongLine) {
	loadbuf_size = kMaxLongLine;
      }
#endif
      // must have enough space for chromosome and variant IDs
      if (loadbuf_size < 2 * 65536) {
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      unsigned char* loadbuf = (unsigned char*)bigstack_alloc_raw(loadbuf_size);
      g_import_dosage_certainty_int = 1 + (int32_t)(import_dosage_certainty * 32768);
      uintptr_t bgen_geno_max_byte_ct = 6LU * sample_ct;
      if (compression_mode) {
        bgen_geno_max_byte_ct = compressBound(bgen_geno_max_byte_ct);
      }
      if (bgen_geno_max_byte_ct > 0xffffffffU) {
	logerrprint("Error: Too many samples for .bgen format.\n");
	goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
      }
      bgen_geno_max_byte_ct += compression_mode * 4;
      // thread-count-independent:
      //   (everything after "2 *" rounded up to cacheline)
      //   compressed_geno_bufs: 2 * bgen_geno_max_byte_ct * main_block_size
      //   g_compressed_geno_starts: 2 * sizeof(intptr_t) * main_block_size
      //   g_write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t) *
      //                     main_block_size
      //   g_write_dosage_cts: 2 * sizeof(int32_t) * main_block_size
      //   g_write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t) *
      //                            main_block_size
      //   g_write_dosage_val_bufs (main bottleneck): 2 * sample_ct *
      //                                              sizeof(dosage_t)
      // additional requirement per thread:
      //   g_bgen_geno_bufs: sample_ct * 3 * sizeof(int16_t)

      uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      if ((!compression_mode) && (calc_thread_ct > 2)) {
	// computation doesn't seem to saturate when decompression is involved
	calc_thread_ct = 2;
      }
      if (bigstack_alloc_thread(calc_thread_ct, &ts.threads) ||
	  bigstack_alloc_usip(calc_thread_ct, &g_bgen_geno_bufs)) {
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      const uint32_t sample_ct_x3 = sample_ct * 3;
      for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	if (bigstack_alloc_usi(sample_ct_x3, &(g_bgen_geno_bufs[tidx]))) {
	  goto ox_bgen_to_pgen_ret_NOMEM;
	}
      }
      uintptr_t cachelines_avail_m12 = bigstack_left() / kCacheline;
      // reserve 1/8 of remaining memory for writer
      cachelines_avail_m12 -= cachelines_avail_m12 / 8;
      if (cachelines_avail_m12 < 12) {
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      // we're making 12 allocations; be pessimistic re: rounding
      cachelines_avail_m12 -= 12;
      const uintptr_t bytes_req_per_in_block_variant = 2 * (bgen_geno_max_byte_ct + sizeof(intptr_t) + sample_ctaw2 * sizeof(intptr_t) + sizeof(int32_t) + sample_ctaw * sizeof(intptr_t) + sample_ct * sizeof(dosage_t));
      uintptr_t main_block_size = (cachelines_avail_m12 * kCacheline) / bytes_req_per_in_block_variant;
      if (main_block_size > 65536) {
	main_block_size = 65536;
      } else if (main_block_size < 8) {
	// this threshold is arbitrary
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (calc_thread_ct > main_block_size / 8) {
	calc_thread_ct = main_block_size / 8;
      }
      ts.calc_thread_ct = calc_thread_ct;
      g_calc_thread_ct = calc_thread_ct;
      unsigned char* compressed_geno_bufs[2];
      if (bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[0])) ||
	  bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[1])) ||
	  bigstack_alloc_ucp(main_block_size, &(g_compressed_geno_starts[0])) ||
	  bigstack_alloc_ucp(main_block_size, &(g_compressed_geno_starts[1])) ||
	  bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[0])) ||
	  bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[1])) ||
	  bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[0])) ||
	  bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[1])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[0])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[1])) ||
	  bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[0])) ||
	  bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[1]))) {
	// this should be impossible
	assert(0);
	goto ox_bgen_to_pgen_ret_NOMEM;
      }

      // likely cases are (i) non-hardcall near top of the file, and (ii) no
      // non-hardcalls at all.  to handle the first case efficiently, we want
      // the first blocks to be small so we bail quickly; to handle the second
      // case efficiently, we want large blocks on average.  so we start with
      // a minimal block size and then repeatedly double.
      uint32_t variant_ct = 0;
      uint32_t block_vidx = 0;
      uint32_t cur_block_size = calc_thread_ct;
      uint32_t parity = 0;
      uintptr_t compressed_block_byte_ct = 6LU * sample_ct;
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      for (uint32_t variant_uidx = 0; variant_uidx < raw_variant_ct; ) {
	uint32_t uii;
	if (!fread(&uii, 4, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (uii != sample_ct) {
	  logprint("\n");
	  logerrprint("Error: Unexpected number of samples specified in SNP block header.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	uint16_t snpid_slen;
	if (!fread(&snpid_slen, 2, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (!snpid_chr) {
	  if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	} else {
	  if (!snpid_slen) {
	    logprint("\n");
	    logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	  }
	  if (!fread(loadbuf, snpid_slen, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  loadbuf[snpid_slen] = '\0';
	}
	uint16_t rsid_slen;
	if (!fread(&rsid_slen, 2, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (fseeko(bgenfile, rsid_slen, SEEK_CUR)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	uint16_t chr_name_slen;
	if (!fread(&chr_name_slen, 2, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (!snpid_chr) {
	  if (!chr_name_slen) {
	    logprint("\n");
	    logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	  }
	  if (!fread(loadbuf, chr_name_slen, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  if ((chr_name_slen == 2) && (!memcmp(loadbuf, "NA", 2))) {
	    strcpy((char*)loadbuf, "0");
	    chr_name_slen = 1;
	  } else {
	    loadbuf[chr_name_slen] = '\0';
	  }
	} else {
	  if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  chr_name_slen = snpid_slen;
	}
	int32_t cur_chr_code;
	reterr = get_or_add_chr_code_destructive("--bgen file", 0, allow_extra_chrs, (char*)loadbuf, (char*)(&(loadbuf[chr_name_slen])), cip, &cur_chr_code);
	if (reterr) {
	  goto ox_bgen_to_pgen_ret_1;
	}
	const uint32_t skip = !is_set(cip->chr_mask, cur_chr_code);

	uint32_t cur_bp; // ignore in this pass
	if (!fread(&cur_bp, 4, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}

	// allele count always 2 and not stored when layout=1
	for (uint32_t allele_idx = 0; allele_idx < 2; ++allele_idx) {
	  uint32_t allele_slen;
	  if (!fread(&allele_slen, 4, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  if (fseeko(bgenfile, allele_slen, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	}

	if (compression_mode) {
#ifdef __LP64__
	  compressed_block_byte_ct = 0;
#endif	  
	  if (!fread(&compressed_block_byte_ct, 4, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	}
	++variant_uidx;
	if (!(variant_uidx % 1000)) {
	  printf("\r--bgen: %uk variants scanned.", variant_uidx / 1000);
	  fflush(stdout);
	}
	if (dosage_is_present || skip) {
	  if (fseeko(bgenfile, compressed_block_byte_ct, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  // bugfix (25 Jun 2017): block_vidx should be left unchanged here
	  variant_ct += 1 - skip;
	  continue;
	}
	compressed_geno_starts[block_vidx] = bgen_geno_iter;
	if (compression_mode) {
	  memcpy(bgen_geno_iter, &compressed_block_byte_ct, 4);
	  bgen_geno_iter = &(bgen_geno_iter[4]);
	}
	if (fread_checked(bgen_geno_iter, compressed_block_byte_ct, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	bgen_geno_iter = &(bgen_geno_iter[compressed_block_byte_ct]);
	++block_vidx;
	if (block_vidx == cur_block_size) {
	  parity = 1 - parity;
	  if (ts.thread_func_ptr) {
	    // process *previous* block results
	    join_threads3z(&ts);
	    reterr = g_error_ret;
	    if (reterr) {
	      logprint("\n");
	      logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    dosage_is_present = g_dosage_is_present;
	    if (dosage_is_present) {
	      // don't need to scan for any more dosages
	      stop_threads3z(&ts, &g_cur_block_write_ct);
	      if (!chr_filter_present) {
		break;
	      }
	      continue;
	    }
	  }
	  g_cur_block_write_ct = cur_block_size;
	  ts.thread_func_ptr = bgen11_dosage_scan_thread;
	  if (spawn_threads3z(variant_ct, &ts)) {
	    goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	  }
	  compressed_geno_starts = g_compressed_geno_starts[parity];
	  bgen_geno_iter = compressed_geno_bufs[parity];
	  block_vidx = 0;
	  variant_ct += cur_block_size;
	  if (cur_block_size < main_block_size) {
	    cur_block_size *= 2;
	    if (cur_block_size > main_block_size) {
	      cur_block_size = main_block_size;
	    }
	  }
	}
      }

      if (!chr_filter_present) {
	variant_ct = raw_variant_ct;
      } else {
	variant_ct += block_vidx;
        if (!variant_ct) {
	  logprint("\n");
	  LOGERRPRINTFWW("Error: All %u variant%s in .bgen file skipped due to chromosome filter.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
	  goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	}
      }
      if (ts.thread_func_ptr) {
	join_threads3z(&ts);
	reterr = g_error_ret;
	if (reterr) {
	  logprint("\n");
	  logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	if (block_vidx && (!g_dosage_is_present)) {
	  g_cur_block_write_ct = block_vidx;
	} else {
	  g_cur_block_write_ct = 0;
	}
	ts.is_last_block = 1;
	if (spawn_threads3z(1, &ts)) {
	  goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	}
	join_threads3z(&ts);
	dosage_is_present = g_dosage_is_present;
      }

      if (fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET)) {
	goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      strcpy(outname_end, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (reterr) {
	goto ox_bgen_to_pgen_ret_1;
      }
      unsigned char* spgw_alloc;
      if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

      // Main workflow:
      // 1. Set n=0, load genotype data for first main_block_size variants
      //    while writing .pvar
      //
      // 2. Spawn threads processing batch n genotype data
      // 3. If n>0, write results for block (n-1)
      // 4. Increment n by 1
      // 5. Load/write-.pvar for batch (n+1) unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8. Write results for last block
      //
      // (May be better to change this to use one output buffer instead of 2.)
      uint32_t vidx_start = 0;
      uint32_t prev_block_write_ct = 0;
      parity = 0;
      reinit_threads3z(&ts);
      while (1) {
	uint32_t cur_block_write_ct = 0;
	if (!ts.is_last_block) {
	  cur_block_write_ct = MINV(variant_ct - vidx_start, main_block_size);
	  compressed_geno_starts = g_compressed_geno_starts[parity];
          bgen_geno_iter = compressed_geno_bufs[parity];
	  for (block_vidx = 0; block_vidx < cur_block_write_ct;) {
	    uint32_t uii;
	    if (!fread(&uii, 4, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (uii != sample_ct) {
	      logprint("\n");
	      logerrprint("Error: Unexpected number of samples specified in SNP block header.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    uint16_t snpid_slen;
	    if (!fread(&snpid_slen, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    char* rsid_start = (char*)loadbuf;
	    if (!snpid_chr) {
	      if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	    } else {
	      if (!snpid_slen) {
		logprint("\n");
		logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
		goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	      }
	      if (!fread(loadbuf, snpid_slen, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      loadbuf[snpid_slen] = '\0';
	      rsid_start = (char*)(&(loadbuf[snpid_slen + 1]));
	    }
	    uint16_t rsid_slen;
	    if (!fread(&rsid_slen, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (!rsid_slen) {
	      logprint("\n");
	      logerrprint("Error: Length-0 rsID in .bgen file.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (!fread(rsid_start, rsid_slen, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    char* loadbuf_iter = &(rsid_start[rsid_slen]);
	    char* chr_name_start = loadbuf_iter;
	    uint16_t chr_name_slen;
	    if (!fread(&chr_name_slen, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (!snpid_chr) {
	      if (!chr_name_slen) {
		logprint("\n");
		logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
		goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	      }
	      if (!fread(chr_name_start, chr_name_slen, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      if ((chr_name_slen == 2) && (!memcmp(chr_name_start, "NA", 2))) {
		strcpy(chr_name_start, "0");
		chr_name_slen = 1;
	      } else {
		chr_name_start[chr_name_slen] = '\0';
	      }
	    } else {
	      if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      chr_name_start = (char*)loadbuf;
	      chr_name_slen = snpid_slen;
	    }
	    int32_t cur_chr_code;
	    reterr = get_or_add_chr_code_destructive("--bgen file", 0, allow_extra_chrs, (char*)chr_name_start, &(chr_name_start[chr_name_slen]), cip, &cur_chr_code);
	    if (reterr) {
	      goto ox_bgen_to_pgen_ret_1;
	    }
	    const uint32_t skip = !is_set(cip->chr_mask, cur_chr_code);

	    uint32_t cur_bp;
	    uint32_t a1_slen;
	    if ((!fread(&cur_bp, 4, 1, bgenfile)) ||
		(!fread(&a1_slen, 4, 1, bgenfile))) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (skip) {
	      uint32_t a2_slen;
	      if (fseeko(bgenfile, a1_slen, SEEK_CUR) ||
		  (!fread(&a2_slen, 4, 1, bgenfile)) ||
		  fseeko(bgenfile, a2_slen, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      if (compression_mode) {
#ifdef __LP64__
		compressed_block_byte_ct = 0;
#endif
		if (!fread(&compressed_block_byte_ct, 4, 1, bgenfile)) {
		  goto ox_bgen_to_pgen_ret_READ_FAIL;
		}
	      }
	      if (fseeko(bgenfile, compressed_block_byte_ct, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      continue;
	    }
	    char* a1_ptr = loadbuf_iter;
	    if (!a1_slen) {
	      logprint("\n");
	      logerrprint("Error: Empty allele code in .bgen file.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (a1_slen > 1000000000) {
	      logprint("\n");
	      logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (a1_slen + (uintptr_t)(a1_ptr - ((char*)loadbuf)) > loadbuf_size) {
	      goto ox_bgen_to_pgen_ret_NOMEM;
	    }
	    if (!fread(a1_ptr, a1_slen, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    char* a2_ptr = &(a1_ptr[a1_slen]);
	    uint32_t a2_slen;
	    if (!fread(&a2_slen, 4, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (!a2_slen) {
	      logprint("\n");
	      logerrprint("Error: Empty allele code in .bgen file.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (a2_slen > 1000000000) {
	      logprint("\n");
	      logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (a2_slen + (uintptr_t)(a2_ptr - ((char*)loadbuf)) > loadbuf_size) {
	      goto ox_bgen_to_pgen_ret_NOMEM;
	    }
	    if (!fread(a2_ptr, a2_slen, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (compression_mode) {
#ifdef __LP64__
	      compressed_block_byte_ct = 0;
#endif
	      if (!fread(&compressed_block_byte_ct, 4, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	    }
	    write_iter = chr_name_write(cip, cur_chr_code, write_iter);
	    *write_iter++ = '\t';
	    if (cur_bp > 0x7ffffffe) {
	      logprint("\n");
	      logerrprint("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
	      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	    }
	    write_iter = uint32toa_x(cur_bp, '\t', write_iter);
	    write_iter = memcpyax(write_iter, rsid_start, rsid_slen, '\t');
	    if (prov_ref_allele_second) {
	      uint32_t swap_slen = a1_slen;
	      a1_slen = a2_slen;
	      a2_slen = swap_slen;
	      char* swap_ptr = a1_ptr;
	      a1_ptr = a2_ptr;
	      a2_ptr = swap_ptr;
	    }
	    if ((write_iter >= writebuf_flush) || (a1_slen >= kMaxMediumLine)) {
	      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
		goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	      }
	      write_iter = writebuf;
	    }
	    if (a1_slen < kMaxMediumLine) {
	      write_iter = memcpya(write_iter, a1_ptr, a1_slen);
	    } else {
	      if (fwrite_checked(a1_ptr, a1_slen, pvarfile)) {
		goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	      }
	    }
	    *write_iter++ = '\t';
	    if ((write_iter >= writebuf_flush) || (a2_slen >= kMaxMediumLine)) {
	      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
		goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	      }
	      write_iter = writebuf;
	    }
	    if (a2_slen < kMaxMediumLine) {
	      write_iter = memcpya(write_iter, a2_ptr, a2_slen);
	    } else {
	      if (fwrite_checked(a2_ptr, a2_slen, pvarfile)) {
		goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	      }
	    }
	    append_binary_eoln(&write_iter);

	    compressed_geno_starts[block_vidx] = bgen_geno_iter;
	    if (compression_mode) {
	      memcpy(bgen_geno_iter, &compressed_block_byte_ct, 4);
	      bgen_geno_iter = &(bgen_geno_iter[4]);
	    }
	    if (fread_checked(bgen_geno_iter, compressed_block_byte_ct, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    bgen_geno_iter = &(bgen_geno_iter[compressed_block_byte_ct]);
	    ++block_vidx;
	  }
	}
	if (vidx_start) {
	  join_threads3z(&ts);
	  reterr = g_error_ret;
	  if (reterr) {
	    logprint("\n");
	    logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	}
	if (!ts.is_last_block) {
	  g_cur_block_write_ct = cur_block_write_ct;
	  ts.is_last_block = (vidx_start + cur_block_write_ct == variant_ct);
	  ts.thread_func_ptr = bgen11_geno_to_pgen_thread;
	  if (spawn_threads3z(vidx_start, &ts)) {
	    goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	  }
	}
	parity = 1 - parity;
	if (vidx_start) {
	  // write *previous* block results
	  uintptr_t* write_genovec_iter = g_write_genovecs[parity];
	  uint32_t* write_dosage_ct_iter = g_write_dosage_cts[parity];
	  uintptr_t* write_dosage_present_iter = g_write_dosage_presents[parity];
	  dosage_t* write_dosage_vals_iter = g_write_dosage_val_bufs[parity];
	  for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx < vidx_start; ++vidx) {
	    const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
	    if (!cur_dosage_ct) {
	      if (spgw_append_biallelic_genovec(write_genovec_iter, &spgw)) {
		goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	      }
	    } else {
	      if (spgw_append_biallelic_genovec_dosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_vals_iter, cur_dosage_ct, &spgw)) {
		goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	      }
	    }
            write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
	    write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
	    write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
	  }
	}
	if (vidx_start == variant_ct) {
	  break;
	}
	if (vidx_start) {
	  printf("\r--bgen: %uk variants converted.", vidx_start / 1000);
	  if (vidx_start <= main_block_size) {
	    fputs("    \b\b\b\b", stdout);
	  }
	  fflush(stdout);
	}
	vidx_start += cur_block_write_ct;
	prev_block_write_ct = cur_block_write_ct;
      }
    } else {
      // v1.2-1.3

      uintptr_t* allele_idx_offsets;
      if (bigstack_end_alloc_ul(raw_variant_ct + 1, &allele_idx_offsets)) {
	logerrprint("error path 1\n");
	goto ox_bgen_to_pgen_ret_NOMEM;
      }

      g_bgen_import_dosage_certainty_thresholds = nullptr;
      if (import_dosage_certainty > (1.0 - kSmallEpsilon) / 3.0) {
	g_bgen_import_dosage_certainty_thresholds = (uint32_t*)bigstack_alloc_raw_rd(25 * sizeof(int32_t));
	for (uint32_t bit_precision = 1; bit_precision <= 16; ++bit_precision) {
	  const uint32_t denom = (1U << bit_precision) - 1;
	  g_bgen_import_dosage_certainty_thresholds[bit_precision] = 1 + (int32_t)(import_dosage_certainty * ((int32_t)denom));
	}
      }
      // bugfix (2 Jul 2017): if max_thread_ct == 1 but there's >12GB memory,
      //   limit to 1 thread rather than (max_thread_ct - 1)...
      const uint32_t calc_thread_ct_limit = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      if (bigstack_alloc_thread(calc_thread_ct_limit, &ts.threads)) {
	logerrprint("error path 2\n");
	goto ox_bgen_to_pgen_ret_NOMEM;
      }

      g_thread_wkspaces = (unsigned char**)bigstack_alloc_raw_rd(calc_thread_ct_limit * sizeof(intptr_t));
      g_thread_bidxs[0] = (uint32_t*)bigstack_alloc_raw_rd((calc_thread_ct_limit + 1) * sizeof(int32_t));
      g_thread_bidxs[1] = (uint32_t*)bigstack_alloc_raw_rd((calc_thread_ct_limit + 1) * sizeof(int32_t));
      // ***** all allocations from this point on are reset before pass 2 *****
      uintptr_t main_block_size = 65536;
      if (bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[0])) ||
	  bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[1])) ||
	  bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[0])) ||
	  bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[1]))) {
	logerrprint("error path 3\n");
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (compression_mode) {
	if (bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[0])) ||
	    bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[1]))) {
	  logerrprint("error path 4\n");
	  goto ox_bgen_to_pgen_ret_NOMEM;
	}
      } else {
	// defensive
	g_uncompressed_genodata_byte_cts[0] = nullptr;
	g_uncompressed_genodata_byte_cts[1] = nullptr;
      }

      // ploidy >2 is not supported by PLINK 2.  (A future build may have code
      // to treat those calls as missing instead of erroring out, as is done
      // with VCF ploidy >2.  But I'll wait until this case actually comes up
      // in the wild...)
      // But even without that, the diploid worst case of 65535 alleles,
      // unphased 32-bit probabilities blows past the 4GB uncompressed record
      // size limit with just 1 sample!  Consequences:
      // * A simple way to avoid unnecessary NOMEM errors is to give each
      //   thread 4GB of decompression workspace on the first pass.  This may
      //   greatly reduce the number of decompression worker threads we can
      //   deploy, but for the first pass that's acceptable: the worker threads
      //   will usually all exit almost immediately (since we just need to
      //   determine whether *any* phase/dosage info needs to be saved).
      // * Even 1 thread x 4GB won't always be available, especially since we
      //   have a double-buffering workflow which requires additional
      //   allocations summing to more than twice the decompression workspace.
      //   So we need to be able to fall back to a smaller decompression
      //   workspace size, and throw NOMEM when it's insufficient.
      // * Of course, records will almost always be far smaller than 4GB.
      //   During the first pass, we'll see every uncompressed record size
      //   (even if the decompression worker threads terminate early), so we
      //   can usually increase the number of worker threads before the second
      //   pass.
      // Overall memory allocation for first pass:
      //   loadbuf_size (~1/7, up to 2GB) : Chromosome code/variant ID/allele
      //                                    code load buffer.
      //   mainbuf_size (~2/7) : Compressed genotype data buffer 0, up to 4GB
      //                         per decompression thread
      //   mainbuf_size        : Compressed genotype data buffer 1
      //   mainbuf_size        : Decompression thread workspace(s)
      // Second pass:
      //   mainbuf_size (~1/6) : Decompression thread workspaces.
      //   16K                 : .bgen chromosome code load buffer.
      //   remainder (~5/6)    : Compressed genotype data buffers, writer, and
      //                         write buffers.
      uintptr_t loadbuf_size = round_down_pow2(bigstack_left() / 7, kCacheline);
      if (loadbuf_size > kMaxLongLine) {
	loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size < 2 * 65536) {
	// don't want to worry about chromosome/variant ID buffer space checks
	// in inner loop
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      unsigned char* loadbuf = bigstack_alloc_raw(loadbuf_size);

      uintptr_t mainbuf_size = round_down_pow2(bigstack_left() / 3, kCacheline);
      uint32_t calc_thread_ct = 1;
      uintptr_t thread_wkspace_size;
#ifdef __LP64__
      // hard compressed and uncompressed record length limits of 2^31 - 1
      // bytes, since these are represented as uint32s in the file.
      if (mainbuf_size > 0x100000000LLU) {
	thread_wkspace_size = 0x100000000LLU;
	mainbuf_size &= 0xffffffff00000000LLU;
	calc_thread_ct = mainbuf_size >> 32;
	if (calc_thread_ct > calc_thread_ct_limit) {
	  calc_thread_ct = calc_thread_ct_limit;
	  mainbuf_size = ((uintptr_t)calc_thread_ct_limit) << 32;
	}
      } else {
	thread_wkspace_size = mainbuf_size;
      }
#else
      thread_wkspace_size = mainbuf_size;
#endif
      // note that thread_wkspace_size is the size limit for a compressed
      // variant record *and* the uncompressed form

      if (main_block_size > raw_variant_ct + calc_thread_ct - 1) {
	main_block_size = raw_variant_ct + calc_thread_ct - 1;
      }
      uint32_t per_thread_block_limit = main_block_size / calc_thread_ct;
      // may as well guarantee divisibility
      main_block_size = per_thread_block_limit * calc_thread_ct;
      ts.calc_thread_ct = calc_thread_ct;
      g_calc_thread_ct = calc_thread_ct;
      unsigned char* compressed_geno_bufs[2];
      compressed_geno_bufs[0] = bigstack_alloc_raw(mainbuf_size);
      compressed_geno_bufs[1] = bigstack_alloc_raw(mainbuf_size);
      for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	g_thread_wkspaces[tidx] = bigstack_alloc_raw(thread_wkspace_size);
      }

      uint32_t variant_ct = 0;

      uint32_t block_vidx = 0;

      // bgen-1.2 and -1.3 records can vary wildly in size, so we're a bit more
      // careful with load balancing here.
      uint32_t cur_per_thread_block_limit = 1;
      uint32_t cur_thread_block_vidx_limit = 1;
      uint32_t cur_thread_fill_idx = 0;

      uint32_t parity = 0;
      uint32_t* thread_bidxs = g_thread_bidxs[0];
      uint16_t* bgen_allele_cts = g_bgen_allele_cts[0];
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[0];
      uint32_t* uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      unsigned char* cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
      thread_bidxs[0] = 0;
      compressed_geno_starts[0] = bgen_geno_iter;
      uintptr_t* allele_idx_offsets_iter = allele_idx_offsets;
      uintptr_t tot_allele_ct = 0;
      uint32_t max_geno_blen = 0;
      uint32_t uncompressed_genodata_byte_ct = 0;

      // temporary kludge
      uint32_t multiallelic_skip_ct = 0;

      g_cur_block_write_ct = 1; // just used as a flag

      for (uint32_t variant_uidx = 0; variant_uidx < raw_variant_ct; ) {
	// format is mostly identical to bgen 1.1; but there's no sample count,
	// and there is an allele count
	// logic is more similar to the second bgen 1.1 pass since we write the
	// .pvar here.
	uint16_t snpid_slen;
	if (!fread(&snpid_slen, 2, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	char* rsid_start = (char*)loadbuf;
	if (!snpid_chr) {
	  if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	} else {
	  if (!snpid_slen) {
	    logprint("\n");
	    logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	  }
	  if (!fread(loadbuf, snpid_slen, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  loadbuf[snpid_slen] = '\0';
	  rsid_start = (char*)(&(loadbuf[snpid_slen + 1]));
	}
	uint16_t rsid_slen;
	if (!fread(&rsid_slen, 2, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (!rsid_slen) {
	  logprint("\n");
	  logerrprint("Error: Length-0 rsID in .bgen file.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	if (!fread(rsid_start, rsid_slen, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	char* loadbuf_iter = &(rsid_start[rsid_slen]);
	char* chr_name_start = loadbuf_iter;
	uint16_t chr_name_slen;
	if (!fread(&chr_name_slen, 2, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (!snpid_chr) {
	  if (!chr_name_slen) {
	    logprint("\n");
	    logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
	  }
	  if (!fread(chr_name_start, chr_name_slen, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  if ((chr_name_slen == 2) && (!memcmp(chr_name_start, "NA", 2))) {
	    strcpy(chr_name_start, "0");
	    chr_name_slen = 1;
	  } else {
	    chr_name_start[chr_name_slen] = '\0';
	  }
	} else {
	  if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  chr_name_start = (char*)loadbuf;
	  chr_name_slen = snpid_slen;
	}
	// chromosome ID length restriction enforced here, so we don't check
	// earlier
	int32_t cur_chr_code;
	reterr = get_or_add_chr_code_destructive("--bgen file", 0, allow_extra_chrs, (char*)chr_name_start, &(chr_name_start[chr_name_slen]), cip, &cur_chr_code);
	if (reterr) {
	  goto ox_bgen_to_pgen_ret_1;
	}
	uint32_t skip = !is_set(cip->chr_mask, cur_chr_code);

	uint32_t cur_bp;
	uint32_t cur_allele_ct = 0;
	if ((!fread(&cur_bp, 4, 1, bgenfile)) ||
	    (!fread(&cur_allele_ct, 2, 1, bgenfile))) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (cur_allele_ct < 2) {
	  // this is undefined in the 1.3 standard; prohibit for now
	  logprint("\n");
	  logerrprint("Error: .bgen variant has fewer than two alleles.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	++variant_uidx;
	if (!(variant_uidx % 1000)) {
	  printf("\r--bgen: %uk variants scanned.", variant_uidx / 1000);
	  fflush(stdout);
	}

	// the "cur_allele_ct > 2" part is a temporary kludge
	if (skip || (cur_allele_ct > 2)) {
	  if (!skip) {
	    ++multiallelic_skip_ct;
	  }
	  for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
	    uint32_t cur_allele_slen;
	    if ((!fread(&cur_allele_slen, 4, 1, bgenfile)) ||
		fseeko(bgenfile, cur_allele_slen, SEEK_CUR)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	  }
	  uint32_t genodata_byte_ct;
	  if ((!fread(&genodata_byte_ct, 4, 1, bgenfile)) ||
	      fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  continue;
	}
	if (rsid_slen > kMaxIdSlen) {
	  // enforce this iff we aren't skipping
	  logprint("\n");
	  logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	// special handling of first two alleles since either may be
	// reference, so we may need to swap order
	char* a1_ptr = loadbuf_iter;
	uint32_t a1_slen;
	if (!fread(&a1_slen, 4, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (!a1_slen) {
	  logprint("\n");
	  logerrprint("Error: Empty allele code in .bgen file.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	if (a1_slen > 1000000000) {
	  logprint("\n");
	  logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	if (a1_slen + (uintptr_t)(a1_ptr - ((char*)loadbuf)) > loadbuf_size) {
	  goto ox_bgen_to_pgen_ret_NOMEM;
	}
	if (!fread(a1_ptr, a1_slen, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	char* a2_ptr = &(a1_ptr[a1_slen]);
	uint32_t a2_slen;
	if (!fread(&a2_slen, 4, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (!a2_slen) {
	  logprint("\n");
	  logerrprint("Error: Empty allele code in .bgen file.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	if (a2_slen > 1000000000) {
	  logprint("\n");
	  logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	if (a2_slen + (uintptr_t)(a2_ptr - ((char*)loadbuf)) > loadbuf_size) {
	  goto ox_bgen_to_pgen_ret_NOMEM;
	}
	if (!fread(a2_ptr, a2_slen, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	write_iter = chr_name_write(cip, cur_chr_code, write_iter);
	*write_iter++ = '\t';
	if (cur_bp > 0x7ffffffe) {
	  logprint("\n");
	  logerrprint("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
	  goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	}
	write_iter = uint32toa_x(cur_bp, '\t', write_iter);
	write_iter = memcpyax(write_iter, rsid_start, rsid_slen, '\t');
	if (prov_ref_allele_second) {
	  const uint32_t swap_slen = a1_slen;
	  a1_slen = a2_slen;
	  a2_slen = swap_slen;
	  char* swap_ptr = a1_ptr;
	  a1_ptr = a2_ptr;
	  a2_ptr = swap_ptr;
	}
	// allele codes may be too large for write buffer, so we special-case
	// this instead of using fwrite_ck()
	if ((write_iter >= writebuf_flush) || (a1_slen >= kMaxMediumLine)) {
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	    goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	if (a1_slen < kMaxMediumLine) {
	  write_iter = memcpya(write_iter, a1_ptr, a1_slen);
	} else {
	  if (fwrite_checked(a1_ptr, a1_slen, pvarfile)) {
	    goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	  }
	}
	*write_iter++ = '\t';
	if ((write_iter >= writebuf_flush) || (a2_slen >= kMaxMediumLine)) {
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	    goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	if (a2_slen < kMaxMediumLine) {
	  write_iter = memcpya(write_iter, a2_ptr, a2_slen);
	} else {
	  if (fwrite_checked(a2_ptr, a2_slen, pvarfile)) {
	    goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	  }
	}
	for (uint32_t allele_idx = 2; allele_idx < cur_allele_ct; ++allele_idx) {
	  // (can't actually reach here yet since we're skipping multiallelics
	  // for now)
	  // safe to use entire loadbuf for this
	  assert(0);
	  uint32_t cur_allele_slen;
	  if (!fread(&cur_allele_slen, 4, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  if (!cur_allele_slen) {
	    logprint("\n");
	    logerrprint("Error: Empty allele code in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	  if (cur_allele_slen > 1000000000) {
	    logprint("\n");
	    logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	  if (cur_allele_slen > loadbuf_size) {
	    goto ox_bgen_to_pgen_ret_NOMEM;
	  }
	  if (!fread(loadbuf, cur_allele_slen, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  *write_iter++ = ',';
	  if ((write_iter >= writebuf_flush) || (cur_allele_slen >= kMaxMediumLine)) {
	    if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	      goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	    }
	    write_iter = writebuf;
	  }
	  if (cur_allele_slen < kMaxMediumLine) {
	    write_iter = memcpya(write_iter, loadbuf, cur_allele_slen);
	  } else {
	    if (fwrite_checked(loadbuf, cur_allele_slen, pvarfile)) {
	      goto ox_bgen_to_pgen_ret_WRITE_FAIL;
	    }
	  }
	}

	append_binary_eoln(&write_iter);
	*allele_idx_offsets_iter++ = tot_allele_ct;
	tot_allele_ct += cur_allele_ct;
	uint32_t genodata_byte_ct;
	if (!fread(&genodata_byte_ct, 4, 1, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	if (genodata_byte_ct > max_geno_blen) {
	  max_geno_blen = genodata_byte_ct;
	}
	if (uncompressed_genodata_byte_cts) {
	  if (genodata_byte_ct < 4) {
	    logerrprint("Error: Invalid compressed block length in .bgen file.\n");
	    goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
	  }
	  if (!fread(&uncompressed_genodata_byte_ct, 4, 1, bgenfile)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  if (uncompressed_genodata_byte_ct > max_geno_blen) {
	    max_geno_blen = uncompressed_genodata_byte_ct;
	  }
	  genodata_byte_ct -= 4;
	}
	if (dosage_is_present) {
	  if (fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
	    goto ox_bgen_to_pgen_ret_READ_FAIL;
	  }
	  ++variant_ct;
	  continue;
	}

	if ((block_vidx == cur_thread_block_vidx_limit) || ((uintptr_t)(cur_geno_buf_end - bgen_geno_iter) < genodata_byte_ct)) {
	  if (!block_vidx) {
	    goto ox_bgen_to_pgen_ret_NOMEM;
	  }
	  thread_bidxs[++cur_thread_fill_idx] = block_vidx;
	  if (cur_thread_fill_idx == calc_thread_ct) {
	    parity = 1 - parity;
	    if (ts.thread_func_ptr) {
	      // process *previous* block results
	      join_threads3z(&ts);
	      reterr = g_error_ret;
	      if (reterr) {
		goto ox_bgen_to_pgen_ret_bgen13_thread_fail;
	      }
	      dosage_is_present = g_dosage_is_present;
	      if (dosage_is_present) {
		// don't need to scan for any more dosages
		stop_threads3z(&ts, &g_cur_block_write_ct);

		// however, unlike bgen-1.1 case, we can never do full
		// early-exit since we have to scan for multiallelic variants:
		// writer must be initialized with (i) an accurate variant
		// count, which is affected by skipped multiallelic variants,
		// and (ii) when we no longer skip them, the PgenWriter
		// constructor still needs a maximum allele count so it can
		// allocate properly-sized buffers.
		if (fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
		  goto ox_bgen_to_pgen_ret_READ_FAIL;
		}
		++variant_ct;
		continue;
	      }
	    }
	    ts.thread_func_ptr = bgen13_dosage_or_phase_scan_thread;
	    if (spawn_threads3z(variant_ct, &ts)) {
	      goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	    }
	    compressed_geno_starts = g_compressed_geno_starts[parity];
	    uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[parity];
	    thread_bidxs = g_thread_bidxs[parity];
	    bgen_allele_cts = g_bgen_allele_cts[parity];
	    bgen_geno_iter = compressed_geno_bufs[parity];
	    thread_bidxs[0] = 0;
	    compressed_geno_starts[0] = bgen_geno_iter;
	    variant_ct += block_vidx;
	    block_vidx = 0;
	    if (cur_per_thread_block_limit < per_thread_block_limit) {
	      cur_per_thread_block_limit *= 2;
	      if (cur_per_thread_block_limit > per_thread_block_limit) {
		cur_per_thread_block_limit = per_thread_block_limit;
	      }
	    }
	    cur_thread_block_vidx_limit = 0;
	    cur_thread_fill_idx = 0;
	  }
	  cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
	  cur_thread_block_vidx_limit += cur_per_thread_block_limit;
	}
	bgen_allele_cts[block_vidx] = cur_allele_ct;
	if (uncompressed_genodata_byte_cts) {
	  uncompressed_genodata_byte_cts[block_vidx] = uncompressed_genodata_byte_ct;
	}
	if (fread_checked(bgen_geno_iter, genodata_byte_ct, bgenfile)) {
	  goto ox_bgen_to_pgen_ret_READ_FAIL;
	}
	bgen_geno_iter = &(bgen_geno_iter[genodata_byte_ct]);
	compressed_geno_starts[++block_vidx] = bgen_geno_iter;
      }
      variant_ct += block_vidx;
      if (multiallelic_skip_ct) {
	logprint("\n");
        LOGERRPRINTFWW("Warning: %u multiallelic variant%s skipped (not yet supported).\n", multiallelic_skip_ct, (multiallelic_skip_ct == 1)? "" : "s");
      }
      if (!variant_ct) {
	logprint("\n");
	LOGERRPRINTF("Error: All %u variant%s in .bgen file skipped.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
	goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      if (variant_ct == block_vidx) {
	// with multiple threads, there's no guarantee that even the first
	// decompression job has launched (e.g. there's only 1 variant on the
	// relevant chromosome in the entire .bgen, and calc_thread_ct == 2).
	// (this is not an issue with the bgen-1.1 converter because we error
	// out on variant_ct == 0, and the first block size is 1.)
	thread_bidxs[cur_thread_fill_idx + 1] = block_vidx;
	ts.thread_func_ptr = bgen13_dosage_or_phase_scan_thread;
	if (spawn_threads3z(variant_ct, &ts)) {
	  goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	}
	block_vidx = 0;
      }
      chr_filter_present = (variant_ct + multiallelic_skip_ct != raw_variant_ct);
      if (ts.thread_func_ptr) {
	join_threads3z(&ts);
	reterr = g_error_ret;
	if (reterr) {
	  goto ox_bgen_to_pgen_ret_bgen13_thread_fail;
	}
	if ((!block_vidx) || g_dosage_is_present) {
	  // ignore thread_bidxs[] in this case
	  g_cur_block_write_ct = 0;
	} else {
	  for (; cur_thread_fill_idx < calc_thread_ct; ) {
	    // save endpoint for current thread, and tell any leftover threads
	    // to do nothing
	    thread_bidxs[++cur_thread_fill_idx] = block_vidx;
	  }
	}
	ts.is_last_block = 1;
	if (spawn_threads3z(1, &ts)) {
	  goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	}
	join_threads3z(&ts);
	dosage_is_present = g_dosage_is_present;
      }

      if (tot_allele_ct == variant_ct * 2) {
	allele_idx_offsets = nullptr;
	bigstack_end_reset(bigstack_end_mark);
      } else {
	// not yet possible
	assert(0);
	*allele_idx_offsets_iter = tot_allele_ct;
      }
      if (fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET)) {
	goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      strcpy(outname_end, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = spgw_init_phase1(outname, allele_idx_offsets, nullptr, variant_ct, sample_ct, dosage_is_present? (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent) : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (reterr) {
	goto ox_bgen_to_pgen_ret_1;
      }

      bigstack_reset(g_bgen_allele_cts[0]);

      // only needs to fit chromosome codes in second pass
      loadbuf = bigstack_alloc_raw_rd(kMaxIdBlen);
      unsigned char* spgw_alloc;
      if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
	logerrprint("error path 5\n");
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      // Now that we know max_geno_blen, try to increase calc_thread_ct, and
      // resize g_thread_wkspaces[tidx] (and also resize compressed_geno_bufs[]
      // in next step).
      // Additional *6 in denominator since we want to limit these allocations
      // to 1/6 of remaining workspace.
      thread_wkspace_size = round_up_pow2(max_geno_blen, kCacheline);
      // bugfix (16 Jul 2017): was computing cachelines_avail, not bytes_avail
      uintptr_t bytes_avail = round_down_pow2(bigstack_left() / 6, kCacheline);
      if (calc_thread_ct_limit * thread_wkspace_size <= bytes_avail) {
	calc_thread_ct = calc_thread_ct_limit;
      } else {
	calc_thread_ct = bytes_avail / thread_wkspace_size;
	if (!calc_thread_ct) {
	  goto ox_bgen_to_pgen_ret_NOMEM;
	}
      }
      ts.calc_thread_ct = calc_thread_ct;
      for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	g_thread_wkspaces[tidx] = bigstack_alloc_raw(thread_wkspace_size);
      }
      bytes_avail -= thread_wkspace_size * calc_thread_ct;
      // Per-write-buffer-variant allocations:
      //   g_bgen_allele_cts: 2 * sizeof(int16_t)
      //   g_compressed_geno_starts: 2 * sizeof(intptr_t)
      //   g_uncompressed_genodata_byte_cts: 2 * sizeof(int32_t)
      //     (unless compression_mode == 0)
      //   g_write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t)
      //   g_write_phasepresents: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_phaseinfos: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_dphase_presents: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_dosage_cts: 2 * sizeof(int32_t)
      //   g_write_dphase_cts: 2 * sizeof(int32_t)
      //   g_write_dosage_val_bufs (the big one): 2 * sample_ct * 2 *
      //                                          sizeof(dosage_t)
      //     additional factor of 2 here is due to phased-dosage support.  (not
      //     actually implemented yet, but will be soon.)
      //   g_compressed_geno_bufs (the other big one): up to 2 * max_geno_blen
      //     The "up to" here is due to the possibility that a few variants
      //     require much more space than the rest; unlikely now, but will be
      //     important when multiallelic support is added.  To defend against
      //     that possibility, we limit g_compressed_geno_bufs[] to 50% of the
      //     total allocation here.
      uintptr_t cachelines_avail_m24 = bigstack_left() / kCacheline;
      if (cachelines_avail_m24 < 24) {
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      // we're making up to 24 allocations; be pessimistic re: rounding
      // (g_compressed_geno_starts has +1, but we have enough room for error)
      cachelines_avail_m24 -= 24;
      uintptr_t bytes_req_per_in_block_variant = 2 * (sizeof(int16_t) + sizeof(intptr_t) + sample_ctaw2 * sizeof(intptr_t) + sample_ctaw * 4 * sizeof(intptr_t) + 2 * sizeof(int32_t) + sample_ct * 2 * sizeof(dosage_t));
      if (compression_mode) {
	bytes_req_per_in_block_variant += 2 * sizeof(int32_t);
      }
      // 50% cap
      mainbuf_size = MINV(max_geno_blen, bytes_req_per_in_block_variant);
      // bugfix (16 Jul 2017): forgot to include this term
      // (17 Jul 2017): forgot to multiply by 2
      bytes_req_per_in_block_variant += 2 * mainbuf_size;
      main_block_size = (cachelines_avail_m24 * kCacheline) / bytes_req_per_in_block_variant;
      if (main_block_size > 65536) {
	main_block_size = 65536;
      }
      if (main_block_size > raw_variant_ct + calc_thread_ct - 1) {
	main_block_size = raw_variant_ct + calc_thread_ct - 1;
      }
      per_thread_block_limit = main_block_size / calc_thread_ct;
      // may as well guarantee divisibility
      main_block_size = per_thread_block_limit * calc_thread_ct;
      mainbuf_size *= main_block_size;
      if (mainbuf_size < max_geno_blen) {
	// bugfix (2 Jul 2017): don't error out here if the entire .bgen has
	// e.g. only one variant
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[0])) ||
	  bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[1])) ||
	  bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[0])) ||
	  bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[1])) ||
	  bigstack_alloc_uc(mainbuf_size, &(compressed_geno_bufs[0])) ||
	  bigstack_alloc_uc(mainbuf_size, &(compressed_geno_bufs[1])) ||
	  bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[0])) ||
	  bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[1])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phasepresents[0])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phasepresents[1])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phaseinfos[0])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phaseinfos[1])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[0])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[1])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dphase_presents[0])) ||
	  bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dphase_presents[1])) ||
	  bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[0])) ||
	  bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[1])) ||
	  bigstack_alloc_ui(main_block_size, &(g_write_dphase_cts[0])) ||
	  bigstack_alloc_ui(main_block_size, &(g_write_dphase_cts[1])) ||
	  bigstack_alloc_dosage(sample_ct * 2 * main_block_size, &(g_write_dosage_val_bufs[0])) ||
	  bigstack_alloc_dosage(sample_ct * 2 * main_block_size, &(g_write_dosage_val_bufs[1]))) {
	// this should be impossible
	logerrprint("error path 6\n");
	LOGERRPRINTF("main_block_size: %" PRIuPTR "\n", main_block_size);
	LOGERRPRINTF("mainbuf_size: %" PRIuPTR "\n", mainbuf_size);
	LOGERRPRINTF("sample_ctaw: %u\n", sample_ctaw);
	LOGERRPRINTF("cachelines_avail_m24: %" PRIuPTR "\n", cachelines_avail_m24);
	LOGERRPRINTF("bytes_req_per_in_block_variant: %" PRIuPTR "\n", bytes_req_per_in_block_variant);
	assert(0);
	goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (compression_mode) {
	if (bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[0])) ||
	    bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[1]))) {
	  logerrprint("error path 7\n");
	  assert(0);
	  goto ox_bgen_to_pgen_ret_NOMEM;
	}
      }
      spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

      // Main workflow:
      // 1. Set n=0, load genotype data for first main_block_size variants
      //    while writing .pvar
      //
      // 2. Spawn threads processing batch n genotype data
      // 3. If n>0, write results for block (n-1)
      // 4. Increment n by 1
      // 5. Load/write-.pvar for batch (n+1) unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8. Write results for last block
      //
      // (May be better to change this to use one output buffer instead of 2,
      // due to high memory requirement.)
      uint32_t vidx_start = 0;
      uint32_t prev_block_write_ct = 0;
      uint32_t prev_genodata_byte_ct = 0;
      uint32_t prev_allele_ct = 0;
      uint32_t skip = 0;
      parity = 0;
      reinit_threads3z(&ts);
      g_cur_block_write_ct = 1;
      while (1) {
	uint32_t cur_block_write_ct = 0;
	if (!ts.is_last_block) {
	  const uint32_t block_vidx_limit = variant_ct - vidx_start;
	  cur_thread_block_vidx_limit = MINV(block_vidx_limit, per_thread_block_limit);
	  cur_thread_fill_idx = 0;
	  thread_bidxs = g_thread_bidxs[parity];
	  bgen_allele_cts = g_bgen_allele_cts[parity];
	  compressed_geno_starts = g_compressed_geno_starts[parity];
	  uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[parity];
          bgen_geno_iter = compressed_geno_bufs[parity];
	  cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
	  thread_bidxs[0] = 0;
	  compressed_geno_starts[0] = bgen_geno_iter;
	  block_vidx = 0;
	  // strictly speaking, prev_genodata_byte_ct and genodata_byte_ct
	  // can be collapsed into one variable, as well as
	  // {block_vidx, cur_block_write_ct}, but not a big deal if the
	  // compiler fails to see this
	  uint32_t genodata_byte_ct = prev_genodata_byte_ct;
	  uint32_t cur_allele_ct = prev_allele_ct;
	  if (!genodata_byte_ct) {
	    goto ox_bgen_to_pgen_load13_start;
	  }
	  // we may stop before main_block_size due to insufficient space in
	  // compressed_geno_buf.  if so, the file pointer is right before the
	  // genotype data, rather than at the beginning of a variant record.
	  while (1) {
	    bgen_allele_cts[block_vidx] = cur_allele_ct;
	    if (uncompressed_genodata_byte_cts) {
	      uncompressed_genodata_byte_cts[block_vidx] = uncompressed_genodata_byte_ct;
	    }
	    if (fread_checked(bgen_geno_iter, genodata_byte_ct, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    bgen_geno_iter = &(bgen_geno_iter[genodata_byte_ct]);
	    compressed_geno_starts[++block_vidx] = bgen_geno_iter;

	    uint16_t snpid_slen;
	    // true iff this is the last variant we're keeping in the entire
	    // file
	    if (block_vidx == block_vidx_limit) {
	      for (; cur_thread_fill_idx < calc_thread_ct; ) {
		// save endpoint for current thread, and tell any leftover
		// threads to do nothing
		thread_bidxs[++cur_thread_fill_idx] = block_vidx;
	      }
	      break;
	    }
	  ox_bgen_to_pgen_load13_start:
	    if (!fread(&snpid_slen, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }

	    if (!snpid_chr) {
	      if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	    } else {
	      if (!fread(loadbuf, snpid_slen, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      loadbuf[snpid_slen] = '\0';
	    }
	    uint16_t rsid_slen;
	    if (!fread(&rsid_slen, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (fseeko(bgenfile, rsid_slen, SEEK_CUR)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    uint16_t chr_name_slen;
	    if (!fread(&chr_name_slen, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    if (!snpid_chr) {
	      if (!fread(loadbuf, chr_name_slen, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      if ((chr_name_slen == 2) && (!memcmp(loadbuf, "NA", 2))) {
		strcpy((char*)loadbuf, "0");
		chr_name_slen = 1;
	      } else {
		loadbuf[chr_name_slen] = '\0';
	      }
	    } else {
	      if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      chr_name_slen = snpid_slen;
	    }
	    if (chr_filter_present) {
	      const int32_t cur_chr_code = get_chr_code((char*)loadbuf, cip, chr_name_slen);
	      assert(cur_chr_code >= 0); // we scanned all the variants
	      skip = !is_set(cip->chr_mask, cur_chr_code);
	    }

	    uint32_t cur_bp; // ignore in this pass
	    if (!fread(&cur_bp, 4, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }

	    cur_allele_ct = 0;
	    if (!fread(&cur_allele_ct, 2, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }
	    for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
	      uint32_t allele_slen;
	      if (!fread(&allele_slen, 4, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      if (fseeko(bgenfile, allele_slen, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	    }
	    if (!fread(&genodata_byte_ct, 4, 1, bgenfile)) {
	      goto ox_bgen_to_pgen_ret_READ_FAIL;
	    }

	    // "cur_allele_ct > 2" is temporary kludge
	    if (skip || (cur_allele_ct > 2)) {
	      if (fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      goto ox_bgen_to_pgen_load13_start;
	    }
	    if (uncompressed_genodata_byte_cts) {
	      if (!fread(&uncompressed_genodata_byte_ct, 4, 1, bgenfile)) {
		goto ox_bgen_to_pgen_ret_READ_FAIL;
	      }
	      genodata_byte_ct -= 4;
	    }

	    if ((block_vidx == cur_thread_block_vidx_limit) || ((uintptr_t)(cur_geno_buf_end - bgen_geno_iter) < genodata_byte_ct)) {
	      thread_bidxs[++cur_thread_fill_idx] = block_vidx;
	      if (cur_thread_fill_idx == calc_thread_ct) {
		prev_allele_ct = cur_allele_ct;
		prev_genodata_byte_ct = genodata_byte_ct;
		break;
	      }
	      cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
	      cur_thread_block_vidx_limit = MINV(cur_thread_block_vidx_limit + per_thread_block_limit, block_vidx_limit);
	    }
	  }
	  cur_block_write_ct = block_vidx;
	}
	if (vidx_start) {
	  join_threads3z(&ts);
	  reterr = g_error_ret;
	  if (reterr) {
	    goto ox_bgen_to_pgen_ret_bgen13_thread_fail;
	  }
	}
	if (!ts.is_last_block) {
	  ts.is_last_block = (vidx_start + cur_block_write_ct == variant_ct);
	  ts.thread_func_ptr = bgen13_geno_to_pgen_thread;
	  if (spawn_threads3z(vidx_start, &ts)) {
	    goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
	  }
	}
	parity = 1 - parity;
	if (vidx_start) {
	  // write *previous* block results
	  const uintptr_t* write_genovec_iter = g_write_genovecs[parity];
	  // const uintptr_t* write_phasepresents = g_write_phasepresents[parity];
	  // const uintptr_t* write_phaseinfos = g_write_phaseinfos[parity];
	  const uintptr_t* write_dosage_presents = g_write_dosage_presents[parity];
	  // const uintptr_t* write_dphase_presents = g_write_dphase_presents[parity];
	  const uint32_t* write_dosage_cts = g_write_dosage_cts[parity];
	  const uint32_t* write_dphase_cts = g_write_dphase_cts[parity];
	  const dosage_t* write_dosage_val_bufs = g_write_dosage_val_bufs[parity];
	  for (uintptr_t block_vidx = 0; block_vidx < prev_block_write_ct; ++block_vidx) {
	    const uint32_t cur_dosage_ct = write_dosage_cts[block_vidx];
	    const uint32_t cur_dphase_ct = write_dphase_cts[block_vidx];
	    if (!cur_dphase_ct) {
	      if (!cur_dosage_ct) {
		if (spgw_append_biallelic_genovec(write_genovec_iter, &spgw)) {
		  goto ox_bgen_to_pgen_ret_WRITE_FAIL;
		}
	      } else {
		if (spgw_append_biallelic_genovec_dosage16(write_genovec_iter, &(write_dosage_presents[block_vidx * sample_ctaw]), &(write_dosage_val_bufs[block_vidx * 2 * sample_ct]), cur_dosage_ct, &spgw)) {
		  goto ox_bgen_to_pgen_ret_WRITE_FAIL;
		}
	      }
	    } else {
	      // todo
	    }
	    write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
	  }
	}
	if (vidx_start == variant_ct) {
	  break;
	}
	if (vidx_start) {
	  printf("\r--bgen: %uk variants converted.", vidx_start / 1000);
	  if (vidx_start <= main_block_size) {
	    fputs("    \b\b\b\b", stdout);
	  }
	  fflush(stdout);
	}
	vidx_start += cur_block_write_ct;
	prev_block_write_ct = cur_block_write_ct;
      }
    }
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
	goto ox_bgen_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&pvarfile)) {
      goto ox_bgen_to_pgen_ret_WRITE_FAIL;
    }
    
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--bgen: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar");
    write_iter = strcpya(write_iter, " written");
    if (!dosage_is_present) {
      write_iter = strcpya(write_iter, " (only hardcalls)");
    }
    strcpy(write_iter, ".\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  ox_bgen_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_bgen_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_bgen_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_bgen_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ox_bgen_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_DELIM:
    logerrprintb();
    if (id_delim == '_') {
      logerrprint("If you do not want '_' to be treated as a FID/IID delimiter, use --double-id or\n--const-fid to choose a different method of converting .bgen sample IDs to\nPLINK IDs, or --id-delim to change the FID/IID delimiter.\n");
    }
    reterr = kPglRetInconsistentInput;
    break;
  ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID:
    logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
    reterr = kPglRetMalformedInput;
    break;
  ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
  ox_bgen_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  ox_bgen_to_pgen_ret_bgen13_thread_fail:
    if (reterr == kPglRetMalformedInput) {
      logprint("\n");
      logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
    } else if (reterr == kPglRetNotYetSupported) {
      logprint("\n");
      logerrprint("Error: BGEN import doesn't currently support phased variants, >16-bit\nprobability precision, or ploidy > 2.\n");
    }
  }
 ox_bgen_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  fclose_cond(bgenfile);
  fclose_cond(psamfile);
  fclose_cond(pvarfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

boolerr_t import_legend_cols(const char* fname, uintptr_t line_idx, uint32_t prov_ref_allele_second, char** loadbuf_iter_ptr, char** write_iter_ptr, uint32_t* variant_ct_ptr) {
  {
    if (*variant_ct_ptr == 0x7ffffffd) {
      logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend other\nsoftware, such as PLINK/SEQ, for very deep studies of small numbers of genomes.\n");
      return 1;
    }
    *variant_ct_ptr += 1;
    char* write_iter = *write_iter_ptr;
    *write_iter++ = '\t';
    char* id_start = *loadbuf_iter_ptr;
    char* id_end = token_endnn(id_start);
    const uint32_t id_slen = (uintptr_t)(id_end - id_start);
    if (id_slen > kMaxIdSlen) {
      logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      return 1;
    }
    char* pos_str = skip_initial_spaces(id_end);
    if (!pos_str) {
      goto import_legend_cols_ret_MISSING_TOKENS;
    }
    char* pos_end = token_endnn(pos_str);
    uint32_t cur_bp;
    if (scan_uint_defcap(pos_str, &cur_bp)) {
      LOGPREPRINTFWW("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, fname);
      return 1;
    }
    write_iter = uint32toa_x(cur_bp, '\t', write_iter);
    write_iter = memcpyax(write_iter, id_start, id_slen, '\t');
    char* first_allele_str = skip_initial_spaces(pos_end);
    if (is_eoln_kns(*first_allele_str)) {
      goto import_legend_cols_ret_MISSING_TOKENS;
    }
    char* first_allele_end = token_endnn(first_allele_str);
    char* second_allele_str = skip_initial_spaces(first_allele_end);
    if (is_eoln_kns(*second_allele_str)) {
      goto import_legend_cols_ret_MISSING_TOKENS;
    }
    char* second_allele_end = token_endnn(second_allele_str);
    if (!prov_ref_allele_second) {
      write_iter = memcpyax(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str), '\t');
      write_iter = memcpya(write_iter, second_allele_str, (uintptr_t)(second_allele_end - second_allele_str));
    } else {
      write_iter = memcpyax(write_iter, second_allele_str, (uintptr_t)(second_allele_end - second_allele_str), '\t');
      write_iter = memcpya(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str));
    }
    *write_iter_ptr = write_iter;
    append_binary_eoln(write_iter_ptr);
    *loadbuf_iter_ptr = second_allele_end;
    return 0;
  }
  {
  import_legend_cols_ret_MISSING_TOKENS:
    LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, fname);
    return 1;
  }
}

pglerr_t scan_haps_for_het(char* loadbuf_iter, const char* hapsname, uint32_t sample_ct, uint32_t is_haploid_or_mt, uintptr_t line_idx_haps, uint32_t* at_least_one_het_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t first_hap_char_code = (uint32_t)((unsigned char)(*loadbuf_iter));
      const uint32_t first_hap_int = first_hap_char_code - 48;
      // will .haps files ever support triallelic variants?  don't worry about
      // that for now
      char* post_first_hap = &(loadbuf_iter[1]);
      if ((first_hap_int >= 2) || ((unsigned char)(*post_first_hap) > 32)) {
	if (first_hap_char_code <= 32) {
	  goto scan_haps_for_het_ret_MISSING_TOKENS;
	}
	goto scan_haps_for_het_ret_INVALID_TOKEN;
      }
      char* second_hap = skip_initial_spaces(post_first_hap);
      char* post_second_hap = &(second_hap[1]);
      const uint32_t second_hap_char_code = (uint32_t)((unsigned char)(*second_hap));
      const uint32_t second_hap_int = second_hap_char_code - 48;
      const uint32_t post_second_hap_char_code = (uint32_t)((unsigned char)(*post_second_hap));
      if ((second_hap_int >= 2) || (post_second_hap_char_code > 32)) {
	// if haploid or MT, permit '-' in second column
	if ((!is_haploid_or_mt) || (second_hap_char_code != 45)) {
	  if (second_hap_char_code <= 32) {
	    goto scan_haps_for_het_ret_MISSING_TOKENS;
	  }
	  if ((second_hap_char_code == 45) && (post_second_hap_char_code <= 32)) {
	    goto scan_haps_for_het_ret_HAPLOID_TOKEN;
	  }
	  goto scan_haps_for_het_ret_INVALID_TOKEN;
	}
      } else if (first_hap_int != second_hap_int) {
	*at_least_one_het_ptr = 1;
	break;
      }
      loadbuf_iter = skip_initial_spaces(post_second_hap);
    }
  }
  while (0) {
  scan_haps_for_het_ret_HAPLOID_TOKEN:
    sprintf(g_logbuf, "Error: Haploid/MT-only token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    reterr = kPglRetInconsistentInput;
    break;
  scan_haps_for_het_ret_INVALID_TOKEN:
    sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    reterr = kPglRetMalformedInput;
    break;
  scan_haps_for_het_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_haps, hapsname);
    reterr = kPglRetMalformedInput;
    break;
  }
  return reterr;
}

#ifdef __arm__
  #error "Unaligned accesses in ox_hapslegend_to_pgen()."
#endif
pglerr_t ox_hapslegend_to_pgen(const char* hapsname, const char* legendname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_hapsfile = nullptr;
  gzFile gz_legendfile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t line_idx_haps = 0;
  uintptr_t line_idx_legend = 0;
  pglerr_t reterr = kPglRetSuccess;
  st_pgen_writer_t spgw;
  uintptr_t loadbuf_size;
  spgw_preinit(&spgw);
  {
    uint32_t sfile_sample_ct = 0;
    if (samplename[0]) {
      reterr = ox_sample_to_psam(samplename, ox_missing_code, misc_flags, outname, outname_end, &sfile_sample_ct);
      if (reterr) {
	goto ox_hapslegend_to_pgen_ret_1;
      }
      if (sfile_sample_ct > (kMaxLongLine / 4)) {
	logerrprint("Error: Too many samples for .haps file converter.\n");
	reterr = kPglRetNotYetSupported;
	goto ox_hapslegend_to_pgen_ret_1;
      }
    }
    
    reterr = gzopen_read_checked(hapsname, &gz_hapsfile);
    if (reterr) {
      goto ox_hapslegend_to_pgen_ret_1;
    }
    uintptr_t writebuf_size = bigstack_left() / 2;
    if (writebuf_size < kMaxMediumLine + kCompressStreamBlock + kCacheline) {
      return kPglRetNomem;
#ifdef __LP64__
      // in 32-bit case, kMaxLongLine + kCompressStreamBlock overflows
    } else if (writebuf_size > kMaxLongLine + kCompressStreamBlock) {
      writebuf_size = kMaxLongLine + kCompressStreamBlock;
#endif
    } else {
      writebuf_size &= ~(kCacheline - 1);
    }
    loadbuf_size = writebuf_size - kCompressStreamBlock;
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* writebuf = (char*)bigstack_alloc_raw(writebuf_size);
    char* writebuf_flush = &(writebuf[kCompressStreamBlock]);
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto ox_hapslegend_to_pgen_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya(writebuf, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);
    char* loadbuf_first_token;
    do {
      ++line_idx_haps;
      if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_hapsfile)) {
	  goto ox_hapslegend_to_pgen_ret_READ_FAIL;
	}
	sprintf(g_logbuf, "Error: %s is empty.\n", hapsname);
	goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));
    const uint32_t token_ct = count_tokens(loadbuf_first_token);
    // pass 1: count variants, write .pvar file, may as well also verify
    // there's at least one heterozygous call
    finalize_chrset(misc_flags, cip);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t at_least_one_het = 0;
    uintptr_t variant_skip_ct = 0;
    uint32_t variant_ct = 0;
    uint32_t is_haploid_or_mt = 0;
    uint32_t sample_ct;
    // support both .haps + .legend (.haps expected to contain no header
    // columns), and pure .haps
    if (legendname[0]) {
      assert(ox_single_chr_str);
      if (token_ct % 2) {
	sprintf(g_logbuf, "Error: %s has an odd number of tokens in the first line. (With --haps + --legend, the .haps file is expected to have no header columns.)\n", hapsname);
	goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = token_ct / 2;
      if (sfile_sample_ct && (sfile_sample_ct != sample_ct)) {
	sprintf(g_logbuf, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
	goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (gzrewind(gz_hapsfile)) {
	goto ox_hapslegend_to_pgen_ret_READ_FAIL;
      }
      line_idx_haps = 0;
      int32_t chr_code_raw = get_chr_code_raw(ox_single_chr_str);
      const char* single_chr_str = nullptr;
      uint32_t single_chr_slen;
      char chr_buf[8]; // nothing longer than e.g. "chrMT" for now
      if (chr_code_raw == -1) {
	// command-line parser guarantees that allow_extra_chrs is true here
	single_chr_str = ox_single_chr_str;
	single_chr_slen = strlen(ox_single_chr_str);
      } else {
	uint32_t chr_code = chr_code_raw;
	if (chr_code > cip->max_code) {
	  if (chr_code < kMaxContigs) {
	    logerrprint("Error: --legend chromosome code is not in the chromosome set.\n");
	    goto ox_hapslegend_to_pgen_ret_INVALID_CMDLINE;
	  }
	  chr_code = cip->xymt_codes[chr_code - kMaxContigs];
	  if (((int32_t)chr_code) < 0) {
	    logerrprint("Error: --legend chromosome code is not in the chromosome set.\n");
	    goto ox_hapslegend_to_pgen_ret_INVALID_CMDLINE;
	  }
	}
	if (!is_set(cip->chr_mask, chr_code)) {
	  logerrprint("Error: --legend chromosome code is excluded by chromosome filter.\n");
	  goto ox_hapslegend_to_pgen_ret_INVALID_CMDLINE;
	}
	is_haploid_or_mt = is_set(cip->haploid_mask, chr_code) || (((int32_t)chr_code) == mt_code);
	char* chr_name_end = chr_name_write(cip, chr_code, chr_buf);
	single_chr_str = chr_buf;
	single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
      }
      reterr = gzopen_read_checked(legendname, &gz_legendfile);
      if (reterr) {
	goto ox_hapslegend_to_pgen_ret_1;
      }
      do {
	++line_idx_legend;
	if (!gzgets(gz_legendfile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_legendfile)) {
	    goto ox_hapslegend_to_pgen_ret_READ_FAIL;
	  }
	  sprintf(g_logbuf, "Error: %s is empty.\n", legendname);
	  goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  goto ox_hapslegend_to_pgen_ret_LONG_LINE_LEGEND;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token));
      // require at least 4 columns, in ID/pos/A1/A2 order; header text is
      // permitted to vary.  tolerate and ignore extra columns.
      if (!next_token_mult(loadbuf_first_token, 3)) {
	goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_LEGEND;
      }      
      while (1) {
	++line_idx_legend;
	if (!gzgets(gz_legendfile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_legendfile)) {
	    goto ox_hapslegend_to_pgen_ret_READ_FAIL;
	  }
	  break;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  goto ox_hapslegend_to_pgen_ret_LONG_LINE_LEGEND;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
	if (is_eoln_kns(*loadbuf_first_token)) {
	  continue;
	}
	char* loadbuf_iter = loadbuf_first_token;
	write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
	if (import_legend_cols(legendname, line_idx_legend, prov_ref_allele_second, &loadbuf_iter, &write_iter, &variant_ct)) {
	  putc_unlocked('\n', stdout);
	  goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT;
	}
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	    goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	if (!at_least_one_het) {
	  do {
	    ++line_idx_haps;
	    if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
	      if (!gzeof(gz_hapsfile)) {
		goto ox_hapslegend_to_pgen_ret_READ_FAIL;
	      }
	      sprintf(g_logbuf, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
	      goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
	    }
	    if (!loadbuf[loadbuf_size - 1]) {
	      goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
	    }
	    loadbuf_first_token = skip_initial_spaces(loadbuf);
	  } while (is_eoln_kns(*loadbuf_first_token));
	  reterr = scan_haps_for_het(loadbuf_first_token, hapsname, sample_ct, is_haploid_or_mt, line_idx_haps, &at_least_one_het);
	  if (reterr) {
            putc_unlocked('\n', stdout);
	    wordwrapb(0);
	    logerrprintb();
	    goto ox_hapslegend_to_pgen_ret_1;
	  }
	}
      }
      if (gzrewind(gz_legendfile)) {
	goto ox_hapslegend_to_pgen_ret_READ_FAIL;
      }
    } else {
      assert(!legendname[0]);
      if ((token_ct < 7) || (!(token_ct % 2))) {
	sprintf(g_logbuf, "Error: Unexpected token count in line %" PRIuPTR " of %s (should be odd, >5).\n", line_idx_haps, hapsname);
	goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = (token_ct - 5) / 2;
      if (sfile_sample_ct && (sfile_sample_ct != sample_ct)) {
	sprintf(g_logbuf, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
	goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      while (1) {
	if (!is_eoln_kns(*loadbuf_first_token)) {
	  char* chr_code_end = token_endnn(loadbuf_first_token);
	  char* loadbuf_iter = skip_initial_spaces(chr_code_end);
	  if (is_eoln_kns(*loadbuf_iter)) {
	    goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
	  }
	  int32_t cur_chr_code;
	  reterr = get_or_add_chr_code_destructive("--haps file", line_idx_haps, allow_extra_chrs, loadbuf_first_token, chr_code_end, cip, &cur_chr_code);
	  if (reterr) {
	    goto ox_hapslegend_to_pgen_ret_1;
	  }
	  if (!is_set(cip->chr_mask, cur_chr_code)) {
	    ++variant_skip_ct;
	  } else {
	    is_haploid_or_mt = is_set(cip->haploid_mask, cur_chr_code) || (cur_chr_code == mt_code);
	    write_iter = chr_name_write(cip, cur_chr_code, write_iter);
	    if (import_legend_cols(hapsname, line_idx_haps, prov_ref_allele_second, &loadbuf_iter, &write_iter, &variant_ct)) {
	      putc_unlocked('\n', stdout);
	      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT;
	    }
	    if (write_iter >= writebuf_flush) {
	      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
		goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
	      }
	      write_iter = writebuf;
	    }
	    if (!at_least_one_het) {
	      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
	      if (scan_haps_for_het(loadbuf_iter, hapsname, sample_ct, is_haploid_or_mt, line_idx_haps, &at_least_one_het)) {
		// override InconsistentInput return code since chromosome info
		// was also gathered from .haps file
		goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
	      }
	    }
	  }
	}
	++line_idx_haps;
	if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_hapsfile)) {
	    goto ox_hapslegend_to_pgen_ret_READ_FAIL;
	  }
	  break;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
      }
      if (!variant_ct) {
	sprintf(g_logbuf, "Error: All %" PRIuPTR " variant%s in %s skipped due to chromosome filter.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s", hapsname);
	goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
    }
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
    }
    if (!sfile_sample_ct) {
      // create a dummy .psam file with "per0/per0", "per1/per1", etc. IDs,
      // matching --dummy
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto ox_hapslegend_to_pgen_ret_OPEN_FAIL;
      }
      write_iter = strcpya(writebuf, "#FID\tIID\tSEX" EOLN_STR);
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	write_iter = memcpyl3a(write_iter, "per");
	write_iter = uint32toa(sample_idx, write_iter);
	write_iter = strcpya(write_iter, "\tper");
	write_iter = uint32toa(sample_idx, write_iter);
	write_iter = strcpya(write_iter, "\tNA" EOLN_STR);
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	    goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
      }
      if (write_iter != writebuf) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (gzrewind(gz_hapsfile)) {
      goto ox_hapslegend_to_pgen_ret_READ_FAIL;
    }
    line_idx_haps = 0;
    bigstack_reset(writebuf);
    putc_unlocked('\r', stdout);
    LOGPRINTF("--haps%s: %u variant%s scanned.\n", legendname[0]? " + --legend" : "", variant_ct, (variant_ct == 1)? "" : "s");
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, at_least_one_het? kfPgenGlobalHardcallPhasePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto ox_hapslegend_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto ox_hapslegend_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    const uint32_t phaseinfo_match_4char = prov_ref_allele_second? 0x20312030 : 0x20302031;
    const uint32_t phaseinfo_match = 1 + prov_ref_allele_second;
    uintptr_t* genovec;
    uintptr_t* phaseinfo;
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
	bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
      goto ox_hapslegend_to_pgen_ret_NOMEM;
    }

    for (uint32_t vidx = 0; vidx < variant_ct;) {
      ++line_idx_haps;
      if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
	if ((!gzeof(gz_hapsfile)) || (!legendname[0])) {
	  goto ox_hapslegend_to_pgen_ret_READ_FAIL;
	}
	sprintf(g_logbuf, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
	goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_first_token)) {
	continue;
      }
      char* loadbuf_iter = loadbuf_first_token;
      if (!legendname[0]) {
	char* chr_code_end = token_endnn(loadbuf_first_token);
	const int32_t cur_chr_code = get_chr_code_counted(cip, (uintptr_t)(chr_code_end - loadbuf_first_token), loadbuf_first_token);
	if (!is_set(cip->chr_mask, cur_chr_code)) {
	  continue;
	}
	is_haploid_or_mt = is_set(cip->haploid_mask, cur_chr_code) || (cur_chr_code == mt_code);
	loadbuf_iter = next_token_mult(skip_initial_spaces(chr_code_end), 4);
	if (!loadbuf_iter) {
	  goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
	}
      }
      uintptr_t genovec_word_or = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      // optimize common case: autosomal diploid, always exactly one space
      // this loop is time-critical; all my attemps to merge in the haploid/MT
      // case have caused >10% slowdowns
      if ((!is_haploid_or_mt) && ((unsigned char)loadbuf_iter[sample_ct * 4 - 1] < 32)) {
	loadbuf_iter[sample_ct * 4 - 1] = ' ';
	const uint32_t* loadbuf_alias32_iter = (const uint32_t*)loadbuf_iter;
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	  }
	  uintptr_t genovec_word = 0;
	  uint32_t phaseinfo_halfword = 0;
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    // assumes little-endian
	    uint32_t cur_hap_4char = *loadbuf_alias32_iter++;
	    if ((cur_hap_4char & 0xfffefffeU) != 0x20302030) {
	      if ((cur_hap_4char & 0xfffffffeU) == 0x202d2030) {
		// "0 - ", "1 - "
		goto ox_hapslegend_to_pgen_ret_HAPLOID_TOKEN;
	      }
	      // any character < 32?
	      if ((((cur_hap_4char & 0xe0e0e0e0U) * 7) & 0x80808080U) != 0x80808080U) {
		goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
	      }
	      sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
	      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
	    }
	    const uintptr_t new_geno = (cur_hap_4char + (cur_hap_4char >> 16)) & 3;
	    genovec_word |= new_geno << (2 * sample_idx_lowbits);
	    if (cur_hap_4char == phaseinfo_match_4char) {
	      phaseinfo_halfword |= 1U << sample_idx_lowbits;
	    }
	  }
	  genovec[widx] = genovec_word;
	  genovec_word_or |= genovec_word;
	  ((halfword_t*)phaseinfo)[widx] = phaseinfo_halfword;
	  ++widx;
	}
      } else {
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	  }
	  uintptr_t genovec_word = 0;
	  uint32_t phaseinfo_halfword = 0;
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    const uint32_t first_hap_char_code = (uint32_t)((unsigned char)(*loadbuf_iter));
	    const uint32_t first_hap_int = first_hap_char_code - 48;
	    char* post_first_hap = &(loadbuf_iter[1]);
	    if ((first_hap_int >= 2) || ((unsigned char)(*post_first_hap) > 32)) {
	      if (first_hap_char_code <= 32) {
		goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
	      }
	      sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
	      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
	    }
	    char* second_hap = skip_initial_spaces(post_first_hap);
	    char* post_second_hap = &(second_hap[1]);
	    const uint32_t second_hap_char_code = (uint32_t)((unsigned char)(*second_hap));
	    const uint32_t post_second_hap_char_code = (uint32_t)((unsigned char)(*post_second_hap));
	    uint32_t second_hap_int = second_hap_char_code - 48;
	    if ((second_hap_int >= 2) || (post_second_hap_char_code > 32)) {
	      if (is_haploid_or_mt && (second_hap_char_code == 45)) {
		// could require --sample, and require this sample to be male
		// in this case?
		second_hap_int = first_hap_int;
	      } else {
		if (second_hap_char_code <= 32) {
		  goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
		}
		if ((second_hap_char_code == 45) && (post_second_hap_char_code <= 32)) {
		  goto ox_hapslegend_to_pgen_ret_HAPLOID_TOKEN;
		}
		sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
		goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
	      }
	    }
	    genovec_word |= ((uintptr_t)(first_hap_int + second_hap_int)) << (2 * sample_idx_lowbits);
	    if (first_hap_int + 2 * second_hap_int == phaseinfo_match) {
	      phaseinfo_halfword |= 1U << sample_idx_lowbits;
	    }
	    loadbuf_iter = skip_initial_spaces(post_second_hap);
	  }
	  genovec[widx] = genovec_word;
	  genovec_word_or |= genovec_word;
	  ((halfword_t*)phaseinfo)[widx] = phaseinfo_halfword;
	  ++widx;
	}
      }
      if (prov_ref_allele_second) {
	genovec_invert_unsafe(sample_ct, genovec);
	zero_trailing_quaters(sample_ct, genovec);
      }
      if (genovec_word_or & kMask5555) {
	if (spgw_append_biallelic_genovec_hphase(genovec, nullptr, phaseinfo, &spgw)) {
	  goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
	}
      } else {
	if (spgw_append_biallelic_genovec(genovec, &spgw)) {
	  goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
	}
      }
      if (!(++vidx % 1000)) {
	printf("\r--haps%s: %uk variants converted.", legendname[0]? " + --legend" : "", vidx / 1000);
	fflush(stdout);
      }
    }
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--haps");
    if (legendname[0]) {
      write_iter = strcpya(write_iter, " + --legend");
    }
    write_iter = strcpya(write_iter, ": ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    if (!sfile_sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya(write_iter, ".psam + ");
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    strcpy(write_iter, ".pvar written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  ox_hapslegend_to_pgen_ret_LONG_LINE_HAP:
    if (loadbuf_size == kMaxLongLine) {
      sprintf(g_logbuf, "Line %" PRIuPTR " of %s is pathologically long.\n", line_idx_haps, hapsname);
      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
    }
  ox_hapslegend_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_hapslegend_to_pgen_ret_LONG_LINE_LEGEND:
    if (loadbuf_size == kMaxLongLine) {
      sprintf(g_logbuf, "Line %" PRIuPTR " of %s is pathologically long.\n", line_idx_legend, legendname);
      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
    }
    reterr = kPglRetNomem;
    break;
  ox_hapslegend_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_hapslegend_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_hapslegend_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_hapslegend_to_pgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_haps, hapsname);
  ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW:
    putc_unlocked('\n', stdout);
    wordwrapb(0);
    logerrprintb();
  ox_hapslegend_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_hapslegend_to_pgen_ret_MISSING_TOKENS_LEGEND:
    putc_unlocked('\n', stdout);
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_legend, legendname);
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  ox_hapslegend_to_pgen_ret_HAPLOID_TOKEN:
    putc_unlocked('\n', stdout);
    sprintf(g_logbuf, "Error: Haploid/MT-only token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    wordwrapb(0);
    logerrprintb();
    reterr = legendname[0]? kPglRetInconsistentInput : kPglRetMalformedInput;
    break;
  ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW:
    putc_unlocked('\n', stdout);
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 ox_hapslegend_to_pgen_ret_1:
  spgw_cleanup(&spgw);
  fclose_cond(outfile);
  gzclose_cond(gz_legendfile);
  gzclose_cond(gz_hapsfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}


// could add an option to load_pvar() to not require allele columns, but .map
// is easy enough to write a separate loader for...
CONSTU31(kLoadMapBlockSize, 65536);

// assumes finalize_chrset() has already been called.
static_assert(kMaxContigs <= 65536, "load_map() needs to be updated.");
pglerr_t load_map(const char* mapname, misc_flags_t misc_flags, chr_info_t* cip, uint32_t* max_variant_id_slen_ptr, uint16_t** variant_chr_codes_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, double** variant_cms_ptr, uint32_t* variant_ct_ptr) {
  // caller should call forget_extra_chr_names(1, cip) after finishing import.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    reterr = gzopen_read_checked(mapname, &gz_infile);
    if (reterr) {
      goto load_map_ret_1;
    }
    // Workspace used as follows:
    // |--loadbuf--|--temp-->----|----<- variant IDs --|
    //            1/4                                 end
    // loadbuf is overwritten with the main return arrays at the end.
    loadbuf_size = round_down_pow2(bigstack_left() / 4, kCacheline);
    char* loadbuf = (char*)g_bigstack_base;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto load_map_ret_NOMEM;
    }
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_first_token;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_map_ret_READ_FAIL;
	}
	logerrprint("Error: Empty .map file.\n");
	goto load_map_ret_INCONSISTENT_INPUT;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto load_map_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token) || (*loadbuf_first_token == '#'));
    char* loadbuf_iter = next_token_mult(loadbuf_first_token, 2);
    if (!loadbuf_iter) {
      goto load_map_ret_MISSING_TOKENS;
    }
    uint32_t map_cols = 3;
    loadbuf_iter = next_token(loadbuf_iter);
    if (loadbuf_iter) {
      loadbuf_iter = next_token(loadbuf_iter);
      if (!loadbuf_iter) {
	map_cols = 4;
      } else {
	loadbuf_iter = next_token(loadbuf_iter);
	if (loadbuf_iter) {
	  if (next_token(loadbuf_iter)) {
	    // do NOT permit >6 columns, .bim is ok but .pvar is not
	    // (pointless to support .pvar for legacy formats)
	    sprintf(g_logbuf, "Error: %s is not a .map/.bim file (too many columns).\n", mapname);
	    goto load_map_ret_MALFORMED_INPUT_WW;
	  }
	  map_cols = 4;
	}
      }
    }

    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    unsigned char* tmp_alloc_base = (unsigned char*)(&(loadbuf[loadbuf_size]));
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    uint16_t* cur_chr_codes = nullptr;
    uint32_t* cur_bps = nullptr;
    char** cur_ids = nullptr;
    double* cur_cms = nullptr;
    double cur_cm = 0.0;
    uint32_t at_least_one_nzero_cm = 0;
    uint32_t variant_ct = 0;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	// chrom, id, (cm?), pos
	char* loadbuf_iter = token_endnn(loadbuf_first_token);
	if (!(*loadbuf_iter)) {
	  goto load_map_ret_MISSING_TOKENS;
	}
	int32_t cur_chr_code;
	reterr = get_or_add_chr_code_destructive(".map file", line_idx, allow_extra_chrs, loadbuf_first_token, loadbuf_iter, cip, &cur_chr_code);
	if (reterr) {
	  goto load_map_ret_1;
	}
	if (!is_set(cip->chr_mask, cur_chr_code)) {
	  goto load_map_skip_variant;
	}
        loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[1]));
        if (is_eoln_kns(*loadbuf_iter)) {
	  goto load_map_ret_MISSING_TOKENS;
	}
	char* token_end = token_endnn(loadbuf_iter);
	uint32_t id_slen = (uintptr_t)(token_end - loadbuf_iter);
	if (id_slen > max_variant_id_slen) {
	  max_variant_id_slen = id_slen;
	}
	tmp_alloc_end -= id_slen + 1;
	if (tmp_alloc_end < tmp_alloc_base) {
	  goto load_map_ret_NOMEM;
	}
	memcpyx(tmp_alloc_end, loadbuf_iter, id_slen, '\0');
	loadbuf_iter = skip_initial_spaces(token_end);
	if (is_eoln_kns(*loadbuf_iter)) {
	  goto load_map_ret_MISSING_TOKENS;
	}

	if (map_cols == 4) {
	  char* cm_end = scanadv_double(loadbuf_iter, &cur_cm);
	  if (!cm_end) {
	    sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, mapname);
	    goto load_map_ret_MALFORMED_INPUT_WW;
	  }
	  at_least_one_nzero_cm = (cur_cm != 0.0);
	  loadbuf_iter = next_token(cm_end);
	  if (!loadbuf_iter) {
	    goto load_map_ret_MISSING_TOKENS;
	  }
	}
	int32_t cur_bp;
	if (scan_int_abs_defcap(loadbuf_iter, &cur_bp)) {
	  sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, mapname);
	  goto load_map_ret_MALFORMED_INPUT_WW;
	}
	if (cur_bp < 0) {
	  goto load_map_skip_variant;
	}

	const uint32_t variant_idx_lowbits = variant_ct % kLoadMapBlockSize;
	if (!variant_idx_lowbits) {
	  if ((uintptr_t)(tmp_alloc_end - tmp_alloc_base) <= kLoadMapBlockSize * (sizeof(int16_t) + sizeof(int32_t) + sizeof(intptr_t) + sizeof(double))) {
	    goto load_map_ret_NOMEM;
	  }
	  cur_chr_codes = (uint16_t*)tmp_alloc_base;
	  tmp_alloc_base = (unsigned char*)(&(cur_chr_codes[kLoadMapBlockSize]));
	  cur_bps = (uint32_t*)tmp_alloc_base;
	  tmp_alloc_base = (unsigned char*)(&(cur_bps[kLoadMapBlockSize]));
	  cur_ids = (char**)tmp_alloc_base;
	  tmp_alloc_base = (unsigned char*)(&(cur_ids[kLoadMapBlockSize]));
	  cur_cms = (double*)tmp_alloc_base;
	  tmp_alloc_base = (unsigned char*)(&(cur_cms[kLoadMapBlockSize]));
	}
	cur_chr_codes[variant_idx_lowbits] = (uint32_t)cur_chr_code;
	cur_ids[variant_idx_lowbits] = (char*)tmp_alloc_end;
	cur_cms[variant_idx_lowbits] = cur_cm;
	cur_bps[variant_idx_lowbits] = (uint32_t)cur_bp;
	++variant_ct;
      }
    load_map_skip_variant:
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_map_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto load_map_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line.)\n", line_idx, mapname);
	goto load_map_ret_MALFORMED_INPUT_WW;
      }
    }
    if (max_variant_id_slen > kMaxIdSlen) {
      logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto load_map_ret_MALFORMED_INPUT;
    }

    if (!variant_ct) {
      logerrprint("Error: All variants in .map/.bim file skipped due to chromosome filter.\n");
      goto load_map_ret_INCONSISTENT_INPUT;
    }
    tmp_alloc_base = (unsigned char*)(&(loadbuf[loadbuf_size]));
    // true requirement is weaker, but whatever
    g_bigstack_end = tmp_alloc_base;

    if (bigstack_alloc_usi(variant_ct, variant_chr_codes_ptr) ||
	bigstack_alloc_ui(variant_ct, variant_bps_ptr) ||
	bigstack_alloc_cp(variant_ct, variant_ids_ptr)) {
      goto load_map_ret_NOMEM;
    }
    uint16_t* variant_chr_codes = *variant_chr_codes_ptr;
    uint32_t* variant_bps = *variant_bps_ptr;
    char** variant_ids = *variant_ids_ptr;
    double* variant_cms = nullptr;
    if (at_least_one_nzero_cm) {
      if (bigstack_alloc_d(variant_ct, variant_cms_ptr)) {
        goto load_map_ret_NOMEM;
      }
      variant_cms = *variant_cms_ptr;
    } else {
      *variant_cms_ptr = nullptr;
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
    *variant_ct_ptr = variant_ct;
    const uint32_t full_block_ct = variant_ct / kLoadMapBlockSize;
    bigstack_mark = g_bigstack_base;
    bigstack_end_set(tmp_alloc_end);
    bigstack_end_mark = g_bigstack_end;

    unsigned char* read_iter = tmp_alloc_base;
    for (uint32_t block_idx = 0; block_idx < full_block_ct; ++block_idx) {
      memcpy(&(variant_chr_codes[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(int16_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int16_t)]);
      memcpy(&(variant_bps[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(int32_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int32_t)]);
      memcpy(&(variant_ids[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(intptr_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(intptr_t)]);
      if (at_least_one_nzero_cm) {
	memcpy(&(variant_cms[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(double));
      }
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(double)]);
    }
    const uint32_t variant_ct_lowbits = variant_ct % kLoadMapBlockSize;
    memcpy(&(variant_chr_codes[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(int16_t));
    read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int16_t)]);
    memcpy(&(variant_bps[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(int32_t));
    read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int32_t)]);
    memcpy(&(variant_ids[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(intptr_t));
    if (at_least_one_nzero_cm) {
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(intptr_t)]);
      memcpy(&(variant_cms[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(double));
    }
  }
  while (0) {
  load_map_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, mapname);
      reterr = kPglRetMalformedInput;
      break;
    }
  load_map_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_map_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_map_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, mapname);
  load_map_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  load_map_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  load_map_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 load_map_ret_1:
  // forget_extra_chr_names(1, cip);
  gzclose_cond(gz_infile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

static_assert(sizeof(dosage_t) == 2, "plink1_dosage_to_pgen() needs to be updated.");
pglerr_t plink1_dosage_to_pgen(const char* dosagename, const char* famname, const char* mapname, const char* import_single_chr_str, const plink1_dosage_info_t* pdip, misc_flags_t misc_flags, fam_col_t fam_cols, int32_t missing_pheno, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  // these are not allocated on bigstack, and must be explicitly freed
  pheno_col_t* pheno_cols = nullptr;
  char* pheno_names = nullptr;
  uint32_t pheno_ct = 0;

  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    // 1. Read .fam file.  (May as well support most .psam files too, since
    //    it's the same driver function.  However, unless 'noheader' modifier
    //    is present, SID field cannot be used for disambiguation.)
    uintptr_t max_sample_id_blen = 4;
    uintptr_t max_sid_blen = 0;
    uintptr_t max_paternal_id_blen = 2;
    uintptr_t max_maternal_id_blen = 2;
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    char* sample_ids = nullptr;
    char* sids = nullptr;
    char* paternal_ids = nullptr;
    char* maternal_ids = nullptr;
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    reterr = load_psam(famname, nullptr, fam_cols, 0x7fffffff, missing_pheno, (misc_flags / kfMiscAffection01) & 1, &max_sample_id_blen, &max_sid_blen, &max_paternal_id_blen, &max_maternal_id_blen, &sample_include, &sample_ids, &sids, &paternal_ids, &maternal_ids, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
    if (reterr) {
      goto plink1_dosage_to_pgen_ret_1;
    }

    // 2. Read dosage-file header line if it exists, then write new .psam.
    reterr = gzopen_read_checked(dosagename, &gz_infile);
    if (reterr) {
      goto plink1_dosage_to_pgen_ret_1;
    }

    const uint32_t first_data_col_idx = pdip->skips[0] + pdip->skips[1] + pdip->skips[2] + 3;
    uint32_t sample_ct = 0;
    uint32_t* dosage_sample_idx_to_fam_uidx;
    if (bigstack_end_alloc_ui(raw_sample_ct, &dosage_sample_idx_to_fam_uidx)) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    const plink1_dosage_flags_t flags = pdip->flags;
    if (flags & kfPlink1DosageNoheader) {
      sample_ct = raw_sample_ct;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	dosage_sample_idx_to_fam_uidx[sample_idx] = sample_idx;
      }
    } else {
      fill_ulong_zero(BITCT_TO_WORDCT(raw_sample_ct), sample_include);
      const uint32_t tmp_htable_size = get_htable_fast_size(raw_sample_ct);
      uint32_t* htable_tmp;
      char* idbuf;
      if (bigstack_end_alloc_ui(tmp_htable_size, &htable_tmp) ||
	  bigstack_end_alloc_c(max_sample_id_blen, &idbuf)) {
	goto plink1_dosage_to_pgen_ret_NOMEM;
      }
      const uint32_t duplicate_idx = populate_strbox_htable(sample_ids, raw_sample_ct, max_sample_id_blen, tmp_htable_size, htable_tmp);
      if (duplicate_idx) {
	char* duplicate_sample_id = &(sample_ids[duplicate_idx * max_sample_id_blen]);
	char* duplicate_fid_end = (char*)rawmemchr(duplicate_sample_id, '\t');
	*duplicate_fid_end = ' ';
	sprintf(g_logbuf, "Error: Duplicate sample ID '%s' in .fam file.\n", duplicate_sample_id);
	goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
      }

      loadbuf_size = bigstack_left();
      if (loadbuf_size > kMaxLongLine) {
	loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size <= kMaxMediumLine) {
	goto plink1_dosage_to_pgen_ret_NOMEM;
      }
      // not formally allocated
      char* loadbuf = (char*)g_bigstack_base;
      loadbuf[loadbuf_size - 1] = ' ';
      char* loadbuf_first_token;
      do {
	++line_idx;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  if (!gzeof(gz_infile)) {
	    goto plink1_dosage_to_pgen_ret_READ_FAIL;
	  }
	  sprintf(g_logbuf, "Error: %s is empty.\n", dosagename);
	  goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_WW;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  goto plink1_dosage_to_pgen_ret_LONG_LINE;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token));
      char* loadbuf_iter = next_token_mult(loadbuf_first_token, first_data_col_idx);
      if (!loadbuf_iter) {
	goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
      }
      do {
	char* fid_end = token_endnn(loadbuf_iter);
	char* iid_start = skip_initial_spaces(fid_end);
	if (is_eoln_kns(*iid_start)) {
	  goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	}
	char* iid_end = token_endnn(iid_start);
	const uint32_t fid_slen = (uintptr_t)(fid_end - loadbuf_iter);
	const uint32_t iid_slen = (uintptr_t)(iid_end - iid_start);
	const uint32_t cur_id_slen = fid_slen + iid_slen + 1;
	if (cur_id_slen >= max_sample_id_blen) {
	  logerrprint("Error: .fam file does not contain all sample IDs in dosage file.\n");
	  goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
	}
	char* idbuf_iid = memcpyax(idbuf, loadbuf_iter, fid_slen, '\t');
	memcpyx(idbuf_iid, iid_start, iid_slen, '\0');
	uint32_t sample_uidx = strbox_htable_find(idbuf, sample_ids, htable_tmp, max_sample_id_blen, cur_id_slen, tmp_htable_size);
	if (sample_uidx == 0xffffffffU) {
	  logerrprint("Error: .fam file does not contain all sample IDs in dosage file.\n");
	  goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
	}
	if (is_set(sample_include, sample_uidx)) {
	  idbuf_iid[-1] = ' ';
	  sprintf(g_logbuf, "Error: Duplicate sample ID '%s' in dosage file.\n", idbuf);
	  goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	set_bit(sample_uidx, sample_include);
	dosage_sample_idx_to_fam_uidx[sample_ct++] = sample_uidx;
	loadbuf_iter = skip_initial_spaces(iid_end);
      } while (!is_eoln_kns(*loadbuf_iter));
    }

    strcpy(outname_end, ".psam");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto plink1_dosage_to_pgen_ret_OPEN_FAIL;
    }
    char* writebuf = g_textbuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya(writebuf, "#FID\tIID");
    const uint32_t write_sid = sid_col_required(sample_include, sids, sample_ct, max_sid_blen, 1);
    if (write_sid) {
      write_iter = strcpya(write_iter, "\tSID");
    }
    const uint32_t write_parents = is_parental_info_present(sample_include, paternal_ids, maternal_ids, sample_ct, max_paternal_id_blen, max_maternal_id_blen);
    if (write_parents) {
      write_iter = strcpya(write_iter, "\tPAT\tMAT");
    }
    write_iter = strcpya(write_iter, "\tSEX");
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
    }
    append_binary_eoln(&write_iter);
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (misc_flags & kfMiscKeepAutoconv) {
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      memcpy(output_missing_pheno, "NA", 2);
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t sample_uidx = dosage_sample_idx_to_fam_uidx[sample_idx];
      write_iter = strcpya(write_iter, &(sample_ids[sample_uidx * max_sample_id_blen]));
      if (write_sid) {
	*write_iter++ = '\t';
	write_iter = strcpya(write_iter, &(sids[sample_uidx * max_sid_blen]));
      }
      if (write_parents) {
	*write_iter++ = '\t';
	write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), '\t');
	write_iter = strcpya(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]));
      }
      *write_iter++ = '\t';
      if (is_set(sex_nm, sample_uidx)) {
	*write_iter++ = '2' - is_set(sex_male, sample_uidx);
      } else {
	write_iter = strcpya(write_iter, "NA");
      }
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	    goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	*write_iter++ = '\t';
	write_iter = append_pheno_str(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
    }
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
    }
    // Don't need sample info any more.
    bigstack_end_reset(bigstack_end_mark);

    // 3. Read .map file if it exists.
    uint32_t max_variant_id_slen = 1;
    uint16_t* variant_chr_codes = nullptr;
    uint32_t* variant_bps = nullptr;
    char** variant_ids = nullptr;
    double* variant_cms = nullptr;
    uint32_t* variant_id_htable = nullptr;
    uintptr_t* variant_already_seen = nullptr;
    uint32_t variant_id_htable_size = 0;
    uint32_t map_variant_ct = 0;
    finalize_chrset(misc_flags, cip);
    if (mapname) {
      reterr = load_map(mapname, misc_flags, cip, &max_variant_id_slen, &variant_chr_codes, &variant_bps, &variant_ids, &variant_cms, &map_variant_ct);
      if (reterr) {
	goto plink1_dosage_to_pgen_ret_1;
      }
      const uint32_t map_variant_ctl = BITCT_TO_WORDCT(map_variant_ct);
      if (bigstack_alloc_ul(map_variant_ctl, &variant_already_seen)) {
	goto plink1_dosage_to_pgen_ret_NOMEM;
      }
      fill_all_bits(map_variant_ct, variant_already_seen);
      unsigned char* bigstack_end_mark2 = g_bigstack_end;
      g_bigstack_end = &(g_bigstack_base[round_down_pow2(bigstack_left() / 2, kEndAllocAlign)]); // allow hash table to only use half of available memory
      reterr = alloc_and_populate_id_htable_mt(variant_already_seen, variant_ids, map_variant_ct, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size);
      if (reterr) {
	goto plink1_dosage_to_pgen_ret_1;
      }
      g_bigstack_end = bigstack_end_mark2;
      fill_ulong_zero(map_variant_ctl, variant_already_seen);
    }

    // 4. Dosage file pass 1: count variants, check whether any decimal dosages
    //    need to be saved, write .pvar.
    //
    // Lots of overlap with ox_gen_to_pgen().
    loadbuf_size = bigstack_left() / 2;
    if (loadbuf_size <= kMaxMediumLine) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    loadbuf_size -= kMaxMediumLine;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    writebuf = (char*)bigstack_alloc_raw(kMaxMediumLine + loadbuf_size);
    writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    const uint32_t chr_col_idx = pdip->chr_col_idx;
    const uint32_t check_chr_col = (chr_col_idx != 0xffffffffU);
    if (!check_chr_col) {
      if (import_single_chr_str) {
	int32_t chr_code_raw = get_chr_code_raw(import_single_chr_str);
	if (chr_code_raw == -1) {
	  // command-line parser guarantees that allow_extra_chrs is true here
	  single_chr_str = import_single_chr_str;
	  single_chr_slen = strlen(import_single_chr_str);
	} else {
	  uint32_t chr_code = chr_code_raw;
	  if (chr_code > cip->max_code) {
	    if (chr_code < kMaxContigs) {
	      logerrprint("Error: --import-dosage single-chr= code is not in the chromosome set.\n");
	      goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
	    }
	    chr_code = cip->xymt_codes[chr_code - kMaxContigs];
	    if (((int32_t)chr_code) < 0) {
	      logerrprint("Error: --import-dosage single-chr= code is not in the chromosome set.\n");
	      goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
	    }
	  }
	  if (!is_set(cip->chr_mask, chr_code)) {
	    // could permit this in --allow-no-vars case, but it's silly
	    logerrprint("Error: --import-dosage single-chr= code is excluded by chromosome filter.\n");
	    goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
	  }
	  char* chr_buf = (char*)bigstack_alloc_raw(kCacheline);
	  char* chr_name_end = chr_name_write(cip, chr_code, chr_buf);
	  single_chr_str = chr_buf;
	  single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
	}
      } else {
	// default to "chr0"
	if (!is_set(cip->chr_mask, 0)) {
	  logerrprint("Error: No --import-dosage chromosome information specified, and chr0 excluded.\n");
	  goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
	}
	char* chr_buf = (char*)bigstack_alloc_raw(kCacheline);
	char* chr_name_end = chr_name_write(cip, 0, chr_buf);
	single_chr_str = chr_buf;
	single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
      }
    }
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto plink1_dosage_to_pgen_ret_OPEN_FAIL;
    }
    write_iter = strcpya(writebuf, "#CHROM\tPOS\tID\tREF\tALT");
    if (variant_cms) {
      write_iter = memcpyl3a(write_iter, "\tCM");
    }
    append_binary_eoln(&write_iter);
    // types:
    // 0 = #CHROM
    // 1 = POS
    // 2 = ID
    // 3 = REF
    // 4 = ALT
    // 5 = first data column
    // (command-line parser verifies that CHROM/POS don't collide with
    // anything else)
    uint64_t parse_table[6];
    // high bits = col index, low bits = col type
    const uint32_t id_col_idx = pdip->skips[0];
    const uint32_t prov_ref_allele_second = !(flags & kfPlink1DosageRefFirst);
    uint32_t ref_col_idx = id_col_idx + pdip->skips[1] + 1;
    const uint32_t alt_col_idx = ref_col_idx + (!prov_ref_allele_second);
    ref_col_idx += prov_ref_allele_second;
    parse_table[0] = (((uint64_t)id_col_idx) << 32) + 2;
    parse_table[1] = (((uint64_t)ref_col_idx) << 32) + 3;
    parse_table[2] = (((uint64_t)alt_col_idx) << 32) + 4;
    uint32_t relevant_initial_col_ct = 3;
    if (check_chr_col) {
      parse_table[relevant_initial_col_ct++] = ((uint64_t)chr_col_idx) << 32;
    }
    const uint32_t check_pos_col = (pdip->pos_col_idx != 0xffffffffU);
    if (check_pos_col) {
      parse_table[relevant_initial_col_ct++] = (((uint64_t)(pdip->pos_col_idx)) << 32) + 1;
    }
    qsort(parse_table, relevant_initial_col_ct, sizeof(int64_t), uint64cmp);
    uint32_t col_skips[6];
    uint32_t col_types[6];
    for (uint32_t uii = 0; uii < relevant_initial_col_ct; ++uii) {
      const uint64_t parse_table_entry = parse_table[uii];
      col_skips[uii] = parse_table_entry >> 32;
      col_types[uii] = (uint32_t)parse_table_entry;
    }
    col_skips[relevant_initial_col_ct] = first_data_col_idx;
    col_types[relevant_initial_col_ct++] = 5;
    for (uint32_t uii = relevant_initial_col_ct - 1; uii; --uii) {
      col_skips[uii] -= col_skips[uii - 1];
    }

    double dosage_multiplier = kDosageMid;
    double dosage_ceil = 32767.5 / 16384.0;
    if (flags & kfPlink1DosageFormatSingle01) {
      dosage_multiplier = kDosageMax;
      dosage_ceil = 32767.5 / 32768.0;
    }
    const uint32_t format_triple = (flags / kfPlink1DosageFormatTriple) & 1;
    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    uint32_t dosage_is_present = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    uint32_t variant_uidx = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto plink1_dosage_to_pgen_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	goto plink1_dosage_to_pgen_ret_LONG_LINE_N;
      }
      char* loadbuf_iter = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_iter)) {
	continue;
      }
      char* token_ptrs[6];
      uint32_t token_slens[6];
      for (uint32_t ric_col_idx = 0; ric_col_idx < relevant_initial_col_ct; ++ric_col_idx) {
	const uint32_t cur_col_type = col_types[ric_col_idx];
	loadbuf_iter = next_token_multz(loadbuf_iter, col_skips[ric_col_idx]);
	if (!loadbuf_iter) {
	  goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	}
	token_ptrs[cur_col_type] = loadbuf_iter;
	char* token_end = token_endnn(loadbuf_iter);
	token_slens[cur_col_type] = (uintptr_t)(token_end - loadbuf_iter);
	loadbuf_iter = token_end;
      }
      // ID
      const char* variant_id = token_ptrs[2];
      const uint32_t variant_id_slen = token_slens[2];
      if (map_variant_ct) {
	variant_uidx = variant_id_dupflag_htable_find(variant_id, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
	if (variant_uidx >> 31) {
	  if (variant_uidx == 0xffffffffU) {
	    ++variant_skip_ct;
	    continue;
	  }
	  sprintf(g_logbuf, "Error: Variant ID '%s' appears multiple times in .map file.\n", variant_ids[variant_uidx & 0x7fffffff]);
	  goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	if (is_set(variant_already_seen, variant_uidx)) {
	  sprintf(g_logbuf, "Error: Variant ID '%s' appears multiple times in --import-dosage file.\n", variant_ids[variant_uidx]);
	  goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
	}
	// already performed chromosome filtering
	write_iter = chr_name_write(cip, (uint32_t)variant_chr_codes[variant_uidx], write_iter);
	*write_iter++ = '\t';
	write_iter = uint32toa_x(variant_bps[variant_uidx], '\t', write_iter);
	write_iter = memcpya(write_iter, variant_id, variant_id_slen);
      } else {
	if (variant_id_slen > kMaxIdSlen) {
	  putc_unlocked('\n', stdout);
	  logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
	  goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT;
	}
	// #CHROM
	if (check_chr_col) {
	  char* chr_code_str = token_ptrs[0];
	  char* chr_code_end = &(chr_code_str[token_slens[0]]);
	  int32_t cur_chr_code;
	  reterr = get_or_add_chr_code_destructive("--import-dosage file", line_idx, allow_extra_chrs, chr_code_str, chr_code_end, cip, &cur_chr_code);
	  if (reterr) {
	    goto plink1_dosage_to_pgen_ret_1;
	  }
	  if (!is_set(cip->chr_mask, cur_chr_code)) {
	    ++variant_skip_ct;
	    continue;
	  }
	  write_iter = chr_name_write(cip, cur_chr_code, write_iter);
	} else {
	  write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
	}
	*write_iter++ = '\t';
	// POS
	if (check_pos_col) {
	  char* pos_str = token_ptrs[1];
	  // no need to support negative values here
	  uint32_t cur_bp;
	  if (scan_uint_defcap(pos_str, &cur_bp)) {
	    sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, dosagename);
	    goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
	  }
	  write_iter = uint32toa(cur_bp, write_iter);
	} else {
	  *write_iter++ = '0';
	}
	*write_iter++ = '\t';
        write_iter = memcpya(write_iter, variant_id, variant_id_slen);
      }
      ++variant_ct;
      *write_iter++ = '\t';
      // REF, ALT
      write_iter = memcpyax(write_iter, token_ptrs[3], token_slens[3], '\t');
      write_iter = memcpya(write_iter, token_ptrs[4], token_slens[4]);
      if (variant_cms) {
	*write_iter++ = '\t';
	write_iter = dtoa_g(variant_cms[variant_uidx], write_iter);
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
      if (!dosage_is_present) {
	loadbuf_iter = token_ptrs[5];
	if (flags & kfPlink1DosageFormatSingle) {
	  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	    if (!loadbuf_iter) {
	      goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	    }
	    double a1_dosage;
	    char* str_end = scanadv_double(loadbuf_iter, &a1_dosage);
	    if ((!loadbuf_iter) || (a1_dosage < (0.5 / 32768.0)) || (a1_dosage >= dosage_ceil)) {
	      loadbuf_iter = next_token(loadbuf_iter);
	      continue;
	    }
	    a1_dosage *= dosage_multiplier;
	    const uint32_t dosage_int = (uint32_t)(a1_dosage + 0.5);
	    const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
	    if (halfdist < dosage_erase_halfdist) {
	      dosage_is_present = 1;
	      break;
	    }
	    loadbuf_iter = next_token(str_end);
	  }
	} else {
	  // for compatibility with plink 1.x, do not actually parse third
	  // value of each triplet if format=3
	  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	    if (!loadbuf_iter) {
	      goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	    }
	    double prob_2a1;
	    char* str_end = scanadv_double(loadbuf_iter, &prob_2a1);
	    if (!str_end) {
	      loadbuf_iter = next_token_mult(loadbuf_iter, 2 + format_triple);
	      continue;
	    }
	    loadbuf_iter = next_token(str_end);
	    if (!loadbuf_iter) {
	      goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	    }
	    double prob_1a1;
	    str_end = scanadv_double(loadbuf_iter, &prob_1a1);
	    if (!str_end) {
	      loadbuf_iter = next_token_mult(loadbuf_iter, 1 + format_triple);
	      continue;
	    }
	    loadbuf_iter = next_token_mult(str_end, 1 + format_triple);
	    double prob_one_or_two_a1 = prob_2a1 + prob_1a1;
	    if ((prob_2a1 < 0.0) || (prob_1a1 < 0.0) || (prob_one_or_two_a1 > 1.01 * (1 + kSmallEpsilon))) {
	      continue;
	    }
	    if (prob_one_or_two_a1 > 1.0) {
	      const double rescale = 1.0 / prob_one_or_two_a1;
	      prob_2a1 *= rescale;
	      prob_1a1 *= rescale;
	      prob_one_or_two_a1 = 1.0;
	    }
	    const uint32_t dosage_int = (uint32_t)(prob_2a1 * 32768 + prob_1a1 * 16384 + 0.5);
	    const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
	    if ((halfdist < dosage_erase_halfdist) && ((prob_2a1 >= import_dosage_certainty) || (prob_1a1 >= import_dosage_certainty) || (prob_one_or_two_a1 <= 1.0 - import_dosage_certainty))) {
	      dosage_is_present = 1;
	      break;
	    }
	  }
	}
      }
      if (!(variant_ct % 1000)) {
	printf("\r--import-dosage: %uk variants scanned.", variant_ct / 1000);
	fflush(stdout);
      }
    }
    putc_unlocked('\r', stdout);
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
    }
    if (!variant_ct) {
      if (!variant_skip_ct) {
	logerrprint("Error: Empty --import-dosage file.\n");
	goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
      }
      LOGERRPRINTFWW("Error: All %" PRIuPTR " variant%s in --import-dosage file skipped.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s");
      goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
    }
    LOGPRINTF("--import-dosage: %u variant%s scanned%s.\n", variant_ct, (variant_ct == 1)? "" : "s", dosage_is_present? "" : " (all hardcalls)");

    // 5. Dosage file pass 2: write .pgen.
    bigstack_reset(writebuf);
    if (gzrewind(gz_infile)) {
      goto plink1_dosage_to_pgen_ret_READ_FAIL;
    }
    const uintptr_t line_ct = line_idx - 1;
    line_idx = 0;
    if (!(flags & kfPlink1DosageNoheader)) {
      // skip header line again
      char* loadbuf_first_token;
      do {
	++line_idx;
	if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	  goto plink1_dosage_to_pgen_ret_READ_FAIL;
	}
	loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token));
    }
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent: kfPgenGlobal0, (flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto plink1_dosage_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* dosage_present;
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
	bigstack_alloc_ul(sample_ctl, &dosage_present)) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    dosage_t* dosage_vals = nullptr;
    if (dosage_is_present) {
      if (bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
	goto plink1_dosage_to_pgen_ret_NOMEM;
      }
    }
    if (hard_call_thresh == 0xffffffffU) {
      hard_call_thresh = kDosageMid / 10;
    }
    dosage_ceil = 2.02 * (1 + kSmallEpsilon);
    if (flags & kfPlink1DosageFormatSingle01) {
      dosage_ceil = 1.01 * (1 + kSmallEpsilon);
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    uint32_t vidx = 0;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	goto plink1_dosage_to_pgen_ret_READ_FAIL;
      }
      char* loadbuf_iter = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_iter)) {
	continue;
      }
      if (variant_skip_ct) {
	if (map_variant_ct) {
	  char* variant_id = next_token_multz(loadbuf_iter, id_col_idx);
	  const uint32_t variant_id_slen = strlen_se(variant_id);
	  if (variant_id_dupflag_htable_find(variant_id, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen) == 0xffffffffU) {
	    continue;
	  }
	  loadbuf_iter = next_token_mult(variant_id, first_data_col_idx - id_col_idx);
	} else {
	  char* chr_code_str = next_token_multz(loadbuf_iter, chr_col_idx);
	  char* chr_code_end = token_endnn(chr_code_str);
	  loadbuf_iter = next_token_mult(chr_code_end, first_data_col_idx - chr_col_idx);
	  *chr_code_end = '\0';
	  const uint32_t chr_code = get_chr_code(chr_code_str, cip, (uintptr_t)(chr_code_end - chr_code_str));
	  if (!is_set(cip->chr_mask, chr_code)) {
	    continue;
	  }
	}
      } else {
	loadbuf_iter = next_token_mult(loadbuf_iter, first_data_col_idx);
      }
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      dosage_t* dosage_vals_iter = dosage_vals;
      while (1) {
	if (widx >= sample_ctl2_m1) {
	  if (widx > sample_ctl2_m1) {
	    break;
	  }
	  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	}
	uintptr_t genovec_word = 0;
	uint32_t dosage_present_hw = 0;
	if (flags & kfPlink1DosageFormatSingle) {
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    if (!loadbuf_iter) {
	      goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	    }
	    double a1_dosage;
	    char* str_end = scanadv_double(loadbuf_iter, &a1_dosage);
	    if ((!loadbuf_iter) || (a1_dosage < 0.0) || (a1_dosage > dosage_ceil)) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	      loadbuf_iter = next_token(loadbuf_iter);
	      continue;
	    }
	    loadbuf_iter = next_token(str_end);
	    uint32_t dosage_int = (uint32_t)(a1_dosage * dosage_multiplier + 0.5);
	    if (dosage_int > kDosageMax) {
	      dosage_int = kDosageMax;
	    }
	    const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
	    if (cur_halfdist < hard_call_halfdist) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	    } else {
	      genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
	      if (cur_halfdist >= dosage_erase_halfdist) {
		continue;
	      }
	    }
	    dosage_present_hw |= 1U << sample_idx_lowbits;
	    *dosage_vals_iter++ = dosage_int;
	  }
	} else {
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    if (!loadbuf_iter) {
	      goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	    }
	    double prob_2a1;
	    char* str_end = scanadv_double(loadbuf_iter, &prob_2a1);
	    if (!str_end) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	      loadbuf_iter = next_token_mult(loadbuf_iter, 2 + format_triple);
	      continue;
	    }
	    loadbuf_iter = next_token(str_end);
	    if (!loadbuf_iter) {
	      goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
	    }
	    double prob_1a1;
	    str_end = scanadv_double(loadbuf_iter, &prob_1a1);
	    if (!str_end) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	      loadbuf_iter = next_token_mult(loadbuf_iter, 1 + format_triple);
	      continue;
	    }
	    loadbuf_iter = next_token_mult(str_end, 1 + format_triple);
	    double prob_one_or_two_a1 = prob_2a1 + prob_1a1;
	    if ((prob_2a1 < 0.0) || (prob_1a1 < 0.0) || (prob_one_or_two_a1 > 1.01 * (1 + kSmallEpsilon))) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	      continue;
	    }
	    if (prob_one_or_two_a1 > 1.0) {
	      const double rescale = 1.0 / prob_one_or_two_a1;
	      prob_2a1 *= rescale;
	      prob_1a1 *= rescale;
	      prob_one_or_two_a1 = 1.0;
	    }
	    if ((prob_2a1 < import_dosage_certainty) && (prob_1a1 < import_dosage_certainty) && (prob_one_or_two_a1 > 1.0 - import_dosage_certainty)) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	    }
	    const uint32_t dosage_int = (uint32_t)(prob_2a1 * 32768 + prob_1a1 * 16384 + 0.5);
	    const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
	    if (cur_halfdist < hard_call_halfdist) {
	      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
	    } else {
	      genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
	      if (cur_halfdist >= dosage_erase_halfdist) {
		continue;
	      }
	    }
	    dosage_present_hw |= 1U << sample_idx_lowbits;
	    *dosage_vals_iter++ = dosage_int;
	  }
	}
	genovec[widx] = genovec_word;
	((halfword_t*)dosage_present)[widx] = (halfword_t)dosage_present_hw;
	++widx;
      }
      if (!prov_ref_allele_second) {
	genovec_invert_unsafe(sample_ct, genovec);
	zero_trailing_quaters(sample_ct, genovec);
      }
      if (dosage_vals_iter != dosage_vals) {
	const uint32_t dosage_ct = (uintptr_t)(dosage_vals_iter - dosage_vals);
	if (!prov_ref_allele_second) {
	  biallelic_dosage16_invert(dosage_ct, dosage_vals);
	}
	if (spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, &spgw)) {
	  goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
	}
      } else {
	if (spgw_append_biallelic_genovec(genovec, &spgw)) {
	  goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
	}
      }
      ++vidx;
      if (!(vidx % 1000)) {
	printf("\r--import-dosage: %uk variants converted.", vidx / 1000);
	fflush(stdout);
      }
    } while (line_idx < line_ct);
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--import-dosage: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".psam written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  plink1_dosage_to_pgen_ret_LONG_LINE_N:
    putc_unlocked('\n', stdout);
  plink1_dosage_to_pgen_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, dosagename);
      reterr = kPglRetMalformedInput;
      break;
    }
  plink1_dosage_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  plink1_dosage_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  plink1_dosage_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  plink1_dosage_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  plink1_dosage_to_pgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  plink1_dosage_to_pgen_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, dosagename);
  plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    putc_unlocked('\n', stdout);
    logerrprintb();
  plink1_dosage_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 plink1_dosage_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  forget_extra_chr_names(1, cip);
  fclose_cond(outfile);
  gzclose_cond(gz_infile);
  free_cond(pheno_names);
  cleanup_pheno_cols(pheno_ct, pheno_cols);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


// binary search over cdf is faster than (int)(log(drand)/log(q)) for truncated
// geometric distribution
static uint64_t g_geno_missing_geomdist[kBitsPerWordD2];
static uint64_t g_dosage_geomdist[kBitsPerWordD2];
static uint32_t g_geno_missing_invert = 0;
static uint32_t g_dosage_geomdist_max = 0;

static_assert(sizeof(dosage_t) == 2, "generate_dummy_thread() needs to be updated.");
THREAD_FUNC_DECL generate_dummy_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint64_t* geno_missing_geomdist = g_geno_missing_geomdist;
  const uint64_t* dosage_geomdist = g_dosage_geomdist;
  const uint32_t geno_missing_invert = g_geno_missing_invert;
  const uint32_t geno_missing_check = geno_missing_invert || (geno_missing_geomdist[kBitsPerWordD2 - 1] != 0);
  const uint32_t dosage_is_present = (dosage_geomdist[kBitsPerWordD2 - 1] != 0);
  const uint32_t dosage_geomdist_max = g_dosage_geomdist_max;
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uintptr_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint64_t ullrand = sfmt_genrand_uint64(sfmtp);
  uint32_t rand16_left = 4;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    uintptr_t* write_genovec_iter = &(g_write_genovecs[parity][vidx * sample_ctaw2]);
    uint32_t* write_dosage_ct_iter = &(g_write_dosage_cts[parity][vidx]);
    uintptr_t* write_dosage_present_iter = &(g_write_dosage_presents[parity][vidx * sample_ctaw]);
    dosage_t* write_dosage_vals_iter = &(g_write_dosage_val_bufs[parity][vidx * sample_ctaw]);
    for (; vidx < vidx_end; ++vidx) {
      dosage_t* cur_dosage_vals_iter = write_dosage_vals_iter;
      uint32_t loop_len = kBitsPerWordD2;
      uint32_t widx = 0;
      while (1) {
	if (widx >= sample_ctl2_m1) {
	  if (widx > sample_ctl2_m1) {
	    break;
	  }
	  loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
	}
	// sfmt_genrand_uint64 calls can't be mixed with sfmt_genrand_uint32
	// calls, so use it here even in 32-bit build
	uintptr_t genovec_word = sfmt_genrand_uint64(sfmtp);
	genovec_word = genovec_word - ((genovec_word >> 1) & kMask5555);
	if (geno_missing_check) {
	  uintptr_t missing_mask = 0;
	  uint32_t sample_idx_lowbits = 0;
	  while (1) {
	    sample_idx_lowbits += uint64arr_geq(geno_missing_geomdist, kBitsPerWordD2, sfmt_genrand_uint64(sfmtp));
	    if (sample_idx_lowbits >= loop_len) {
	      break;
	    }
	    missing_mask |= (3 * k1LU) << (2 * sample_idx_lowbits);
	    ++sample_idx_lowbits;
	  }
	  if (geno_missing_invert) {
	    missing_mask = ~missing_mask;
	  }
	  genovec_word |= missing_mask;
	}
	uint32_t dosage_present_hw = 0;
	if (dosage_is_present) {
	  // deliberate overflow
	  uint32_t sample_idx_lowbits = 0xffffffffU;
	  while (1) {
	    ++sample_idx_lowbits;
	    if (dosage_geomdist_max) {
	      sample_idx_lowbits += uint64arr_geq(dosage_geomdist, dosage_geomdist_max, sfmt_genrand_uint64(sfmtp));
	    }
	    if (sample_idx_lowbits >= loop_len) {
	      break;
	    }
	    if (((genovec_word >> (2 * sample_idx_lowbits)) & 3) == 3) {
	      continue;
	    }
	    if (!rand16_left) {
	      ullrand = sfmt_genrand_uint64(sfmtp);
	      rand16_left = 4;
	    }
	    const uint32_t dosage_int = ((ullrand & 65535) + 1) / 2;
	    ullrand >>= 16;
	    --rand16_left;
	    const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
	    if (halfdist < dosage_erase_halfdist) {
	      *cur_dosage_vals_iter++ = dosage_int;
	      dosage_present_hw |= 1U << sample_idx_lowbits;
	      if (halfdist < hard_call_halfdist) {
		genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
		continue;
	      }
	    }
	    genovec_word &= ~((3 * k1LU) << (2 * sample_idx_lowbits));
	    genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
	  }
	}
	write_genovec_iter[widx] = genovec_word;
	((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
	++widx;
      }
      zero_trailing_quaters(sample_ct, write_genovec_iter);
      const uint32_t dosage_ct = (uintptr_t)(cur_dosage_vals_iter - write_dosage_vals_iter);
      *write_dosage_ct_iter++ = dosage_ct;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
      write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static_assert(sizeof(dosage_t) == 2, "generate_dummy() needs to be updated.");
pglerr_t generate_dummy(const gendummy_info_t* gendummy_info_ptr, misc_flags_t misc_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  threads_state_t ts;
  init_threads3z(&ts);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    finalize_chrset(misc_flags, cip);
    if (!is_set(cip->chr_mask, 1)) {
      logerrprint("Error: --dummy cannot be used when chromosome 1 is excluded.\n");
      goto generate_dummy_ret_INVALID_CMDLINE;
    }
    if (is_set(cip->haploid_mask, 1)) {
      logerrprint("Error: --dummy cannot be used to generate haploid data.\n");
      goto generate_dummy_ret_INVALID_CMDLINE;
    }
    char chr1_name_buf[5];
    char* chr1_name_end = chr_name_write(cip, 1, chr1_name_buf);
    *chr1_name_end = '\t';
    const uint32_t chr1_name_blen = 1 + (uintptr_t)(chr1_name_end - chr1_name_buf);
    const uint32_t sample_ct = gendummy_info_ptr->sample_ct;
    const uint32_t variant_ct = gendummy_info_ptr->variant_ct;
    // missing pheno string is always "NA"
    const gendummy_flags_t flags = gendummy_info_ptr->flags;
    uint16_t alleles[13];
    uint32_t four_alleles = 0;
    if (flags & kfGenDummyAcgt) {
      memcpy(alleles, "\tA\tC\tA\tG\tA\tT\tC\tG\tC\tT\tG\tT\tA", 26);
      four_alleles = 1;
    } else if (flags & kfGenDummy1234) {
      memcpy(alleles, "\t1\t2\t1\t3\t1\t4\t2\t3\t2\t4\t3\t4\t1", 26);
      four_alleles = 1;
    } else if (flags & kfGenDummy12) {
      memcpy(alleles, "\t1\t2\t1", 6);
    } else {
      memcpy(alleles, "\tA\tB\tA", 6);
    }
    char* textbuf = g_textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto generate_dummy_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya(textbuf, "#CHROM\tPOS\tID\tREF\tALT");
    append_binary_eoln(&write_iter);
    if (four_alleles) {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
	if (!(variant_idx % 8)) {
	  if (write_iter >= textbuf_flush) {
	    if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	      goto generate_dummy_ret_WRITE_FAIL;
	    }
	    write_iter = textbuf;
	  }
	  do {
	    urand = sfmt_genrand_uint32(&g_sfmt);
	  } while (urand < 425132032U); // 2^32 - 12^8
	}
	const uint32_t quotient = urand / 12;
	const uint32_t remainder = urand - (quotient * 12U);
	urand = quotient;
	write_iter = memcpya(write_iter, chr1_name_buf, chr1_name_blen);
	write_iter = uint32toa(variant_idx, write_iter);
	write_iter = strcpya(write_iter, "\tsnp");
	write_iter = uint32toa(variant_idx, write_iter);
	write_iter = memcpya(write_iter, &(alleles[remainder]), 4);
	append_binary_eoln(&write_iter);
      }
    } else {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
	if (!(variant_idx % 32)) {
	  if (write_iter >= textbuf_flush) {
	    if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	      goto generate_dummy_ret_WRITE_FAIL;
	    }
	    write_iter = textbuf;
	  }
	  urand = sfmt_genrand_uint32(&g_sfmt);
	}
	const uint32_t remainder = urand & 1;
	urand >>= 1;
	write_iter = memcpya(write_iter, chr1_name_buf, chr1_name_blen);
	write_iter = uint32toa(variant_idx, write_iter);
	write_iter = strcpya(write_iter, "\tsnp");
	write_iter = uint32toa(variant_idx, write_iter);
	write_iter = memcpya(write_iter, &(alleles[remainder]), 4);
	append_binary_eoln(&write_iter);
      }
    }
    if (write_iter != textbuf) {
      if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto generate_dummy_ret_WRITE_FAIL;
    }

    strcpy(outname_end, ".psam");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto generate_dummy_ret_OPEN_FAIL;
    }
    const uint32_t pheno_ct = gendummy_info_ptr->pheno_ct;
    char* writebuf;
    if (bigstack_alloc_c(kMaxMediumLine + 48 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf)) {
      goto generate_dummy_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (misc_flags & kfMiscKeepAutoconv) {
      // must use --output-missing-phenotype parameter, which we've validated
      // to be consistent with --input-missing-phenotype
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      // use "NA" since that's always safe
      memcpy(output_missing_pheno, "NA", 2);
    }
    write_iter = strcpya(writebuf, "#FID\tIID\tSEX");
    for (uint32_t pheno_idx_p1 = 1; pheno_idx_p1 <= pheno_ct; ++pheno_idx_p1) {
      write_iter = strcpya(write_iter, "\tPHENO");
      write_iter = uint32toa(pheno_idx_p1, write_iter);
    }
    append_binary_eoln(&write_iter);
    const uint32_t pheno_m_check = (gendummy_info_ptr->pheno_mfreq >= kRecip2m32 * 0.5);
    const uint32_t pheno_m32 = (uint32_t)(gendummy_info_ptr->pheno_mfreq * 4294967296.0 - 0.5);
    if ((flags & kfGenDummyScalarPheno) && pheno_ct) {
      uint32_t saved_rnormal = 0;
      double saved_rnormal_val = 0.0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	    goto generate_dummy_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	write_iter = memcpyl3a(write_iter, "per");
	write_iter = uint32toa(sample_idx, write_iter);
	write_iter = strcpya(write_iter, "\tper");
	write_iter = uint32toa(sample_idx, write_iter);
	// could add option to add some males/unknown gender
	write_iter = strcpya(write_iter, "\t2");
	for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  *write_iter++ = '\t';
	  if (pheno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= pheno_m32)) {
	    write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
	  } else {
	    double dxx;
	    if (saved_rnormal) {
	      dxx = saved_rnormal_val;
	    } else {
	      dxx = rand_normal(&g_sfmt, &saved_rnormal_val);
	    }
	    saved_rnormal_val = 1 - saved_rnormal_val;
	    write_iter = dtoa_g(dxx, write_iter);
	  }
	}
	append_binary_eoln(&write_iter);
      }
    } else {
      uint32_t urand = sfmt_genrand_uint32(&g_sfmt);
      uint32_t urand_bits_left = 32;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	    goto generate_dummy_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
	write_iter = memcpyl3a(write_iter, "per");
	write_iter = uint32toa(sample_idx, write_iter);
	write_iter = strcpya(write_iter, "\tper");
	write_iter = uint32toa(sample_idx, write_iter);
	write_iter = strcpya(write_iter, "\t2");
	for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	  *write_iter++ = '\t';
	  if (pheno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= pheno_m32)) {
	    write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
	  } else {
	    if (!urand_bits_left) {
	      urand = sfmt_genrand_uint32(&g_sfmt);
	      urand_bits_left = 32;
	    }
	    *write_iter++ = (char)((urand & 1) + '1');
	    urand >>= 1;
	    --urand_bits_left;
	  }
	}
	append_binary_eoln(&write_iter);
      }
    }
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto generate_dummy_ret_WRITE_FAIL;
    }

    bigstack_reset(writebuf);
    strcpy(outname_end, ".pgen");
    const double geno_mfreq = gendummy_info_ptr->geno_mfreq;
    if (geno_mfreq < kRecip2m53) {
      // beyond this point, 1-x may just be 1
      g_geno_missing_geomdist[kBitsPerWordD2 - 1] = 0;
    } else {
      double remaining_prob = 1.0;
      g_geno_missing_invert = (geno_mfreq > 0.5);
      if (g_geno_missing_invert) {
	for (uint32_t uii = 0; uii < kBitsPerWordD2; ++uii) {
	  remaining_prob *= geno_mfreq;
	  g_geno_missing_geomdist[uii] = -((uint64_t)(remaining_prob * k2m64));
	}
      } else {
        const double geno_nmfreq = 1.0 - geno_mfreq;
	for (uint32_t uii = 0; uii < kBitsPerWordD2; ++uii) {
	  remaining_prob *= geno_nmfreq;
	  g_geno_missing_geomdist[uii] = -((uint64_t)(remaining_prob * k2m64));
	}
      }
    }
    const double dosage_nfreq = 1.0 - gendummy_info_ptr->dosage_freq;
    if (dosage_nfreq >= 1.0) {
      g_dosage_geomdist[kBitsPerWordD2 - 1] = 0;
    } else {
      double remaining_prob = 1.0;
      for (uint32_t uii = 0; uii < kBitsPerWordD2; ++uii) {
	remaining_prob *= dosage_nfreq;
	g_dosage_geomdist[uii] = -((uint64_t)(remaining_prob * k2m64));
      }
      uint32_t dosage_geomdist_max = kBitsPerWordD2;
      for (; dosage_geomdist_max; --dosage_geomdist_max) {
	if (g_dosage_geomdist[dosage_geomdist_max - 1] != 0) {
	  break;
	}
      }
      g_dosage_geomdist_max = dosage_geomdist_max;
    }
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, (dosage_nfreq >= 1.0)? kfPgenGlobal0 : kfPgenGlobalDosagePresent, 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto generate_dummy_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto generate_dummy_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    // thread-count-independent:
    //   (everything after "2 *" rounded up to cacheline)
    //   g_write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t) *
    //                     main_block_size
    //   g_write_dosage_cts: 2 * sizeof(int32_t) * main_block_size
    //   g_write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t) *
    //                            main_block_size
    //   g_write_dosage_val_bufs: 2 * sample_ct * sizeof(dosage_t)
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // saturates around 4 compute threads, both with and without dosage
    // (todo: test this on something other than a MacBook Pro, could just be a
    // hyperthreading artifact)
    if (calc_thread_ct > 4) {
      calc_thread_ct = 4;
    }
    if (bigstack_init_sfmtp(calc_thread_ct, 0)) {
      goto generate_dummy_ret_NOMEM;
    }
    if (bigstack_alloc_thread(calc_thread_ct, &ts.threads)) {
      goto generate_dummy_ret_NOMEM;
    }
    const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
    const uint32_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
    uintptr_t cachelines_avail_m8 = bigstack_left() / kCacheline;
    if (cachelines_avail_m8 < 8) {
      goto generate_dummy_ret_NOMEM;
    }
    // we're making 8 allocations; be pessimistic re: rounding
    cachelines_avail_m8 -= 8;
    const uintptr_t bytes_req_per_in_block_variant = 2 * (sample_ctaw2 * sizeof(intptr_t) + sizeof(int32_t) + sample_ctaw * sizeof(intptr_t) + sample_ct * sizeof(dosage_t));
    uintptr_t main_block_size = (cachelines_avail_m8 * kCacheline) / bytes_req_per_in_block_variant;
    if (main_block_size > 65536) {
      main_block_size = 65536;
    } else if (main_block_size < 8) {
      // this threshold is arbitrary
      goto generate_dummy_ret_NOMEM;
    }
    if (calc_thread_ct > main_block_size / 8) {
      calc_thread_ct = main_block_size / 8;
    }
    ts.calc_thread_ct = calc_thread_ct;
    g_calc_thread_ct = calc_thread_ct;
    g_sample_ct = sample_ct;
    if (bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[0])) ||
	bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[1])) ||
	bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[0])) ||
	bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[1])) ||
	bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[0])) ||
	bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[1])) ||
	bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[0])) ||
	bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[1]))) {
      // this should be impossible
      assert(0);
      goto generate_dummy_ret_NOMEM;
    }
    g_hard_call_halfdist = kDosage4th - hard_call_thresh;
    g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

    // Main workflow:
    // 1. Set n=0
    //
    // 2. Spawn threads generating batch n genotype data
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Join threads
    // 6. Goto step 2 unless eof
    //
    // 7. Write results for last block
    uint32_t vidx_start = 0;
    uint32_t prev_block_write_ct = 0;
    uint32_t parity = 0;
    while (1) {
      uint32_t cur_block_write_ct = 0;
      if (!ts.is_last_block) {
	cur_block_write_ct = MINV(variant_ct - vidx_start, main_block_size);
      }
      if (vidx_start) {
	join_threads3z(&ts);
      }
      if (!ts.is_last_block) {
	g_cur_block_write_ct = cur_block_write_ct;
	ts.is_last_block = (vidx_start + cur_block_write_ct == variant_ct);
	ts.thread_func_ptr = generate_dummy_thread;
	if (spawn_threads3z(vidx_start, &ts)) {
	  goto generate_dummy_ret_THREAD_CREATE_FAIL;
	}
      }
      parity = 1 - parity;
      if (vidx_start) {
	// write *previous* block results
	uintptr_t* write_genovec_iter = g_write_genovecs[parity];
	uint32_t* write_dosage_ct_iter = g_write_dosage_cts[parity];
	uintptr_t* write_dosage_present_iter = g_write_dosage_presents[parity];
	dosage_t* write_dosage_vals_iter = g_write_dosage_val_bufs[parity];
	for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx < vidx_start; ++vidx) {
	  const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
	  if (!cur_dosage_ct) {
	    if (spgw_append_biallelic_genovec(write_genovec_iter, &spgw)) {
	      goto generate_dummy_ret_WRITE_FAIL;
	    }
	  } else {
	    if (spgw_append_biallelic_genovec_dosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_vals_iter, cur_dosage_ct, &spgw)) {
	      goto generate_dummy_ret_WRITE_FAIL;
	    }
	  }
	  write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
	  write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
	  write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
	}
      }
      if (vidx_start == variant_ct) {
	break;
      }
      if (vidx_start) {
	printf("\r--dummy: %uk variants written.", vidx_start / 1000);
	fflush(stdout);
      }
      vidx_start += cur_block_write_ct;
      prev_block_write_ct = cur_block_write_ct;
    }
    spgw_finish(&spgw);

    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    LOGPRINTFWW("Dummy data (%u sample%s, %u SNP%s) written to %s.pgen + %s.pvar + %s.psam .\n", sample_ct, (sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s", outname, outname, outname);
  }
  while (0) {
  generate_dummy_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  generate_dummy_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  generate_dummy_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  generate_dummy_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  generate_dummy_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 generate_dummy_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

/*
#ifdef __arm__
  #error "Unaligned accesses in bitvec_resort()."
#endif
void bitvec_resort(const uintptr_t* bitvec, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, unsigned char* writebuf) {
  const uint32_t sample_ctl_m1 = BITCT_TO_WORDCT(sample_ct) - 1;
  uint32_t widx = 0;
  uint32_t cur_word_entry_ct = kBitsPerWord;
  const uint32_t* new_sample_idx_to_old_base = new_sample_idx_to_old;
  uintptr_t* writebuf_walias = (uintptr_t*)writebuf;
  while (1) {
    if (widx == sample_ctl_m1) {
      cur_word_entry_ct = 1 + ((sample_ct - 1) % kBitsPerWord);
    }
    uintptr_t cur_word = 0;
    for (uint32_t uii = 0; uii < cur_word_entry_ct; ++uii) {
      cur_word |= IS_SET(bitvec, new_sample_idx_to_old_base[uii]) << uii;
    }
    if (widx == sample_ctl_m1) {
      memcpy(&(writebuf_walias[widx]), &cur_word, (cur_word_entry_ct + (CHAR_BIT - 1)) / CHAR_BIT);
      return;
    }
    writebuf_walias[widx++] = cur_word;
    new_sample_idx_to_old_base = &(new_sample_idx_to_old_base[kBitsPerWord]);
  }
}
*/

#ifdef __arm__
  #error "Unaligned accesses in genovec_resort()."
#endif
void genovec_resort(const uintptr_t* genovec, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, unsigned char* writebuf) {
  // writebuf need not be word-aligned
  const uint32_t sample_ctl2_m1 = QUATERCT_TO_WORDCT(sample_ct) - 1;
  uint32_t word_idx = 0;
  uint32_t cur_word_entry_ct = kBitsPerWordD2;
  const uint32_t* new_sample_idx_to_old_iter = new_sample_idx_to_old;
  uintptr_t* writebuf_walias = (uintptr_t*)writebuf;
  while (1) {
    if (word_idx == sample_ctl2_m1) {
      cur_word_entry_ct = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t cur_word = 0;
    for (uint32_t uii = 0; uii < cur_word_entry_ct; ++uii) {
      cur_word |= (GET_QUATERARR_ENTRY(genovec, new_sample_idx_to_old_iter[uii])) << (2 * uii);
    }
    if (word_idx == sample_ctl2_m1) {
      memcpy(&(writebuf_walias[word_idx]), &cur_word, QUATERCT_TO_BYTECT(cur_word_entry_ct));
      return;
    }
    writebuf_walias[word_idx++] = cur_word;
    new_sample_idx_to_old_iter = &(new_sample_idx_to_old_iter[kBitsPerWordD2]);
  }
}

void unpack_hphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, uint32_t raw_sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uintptr_t* phaseraw_iter = phaseraw;
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  uintptr_t phaseraw_word = *phaseraw_iter++;
  uint32_t read_idx_lowbits = 1;
  if (!(phaseraw_word & 1)) {
    // phase always present
    phaseraw_word >>= 1;
    *phasepresent_ptr = nullptr;
    for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
      uintptr_t new_phasepresent_word = all_hets[widx];
      uintptr_t new_phaseinfo_word = 0;
      while (new_phasepresent_word) {
	// this copies over bottom bit of new_phasepresent_word, retaining its
	// position
	if (read_idx_lowbits == kBitsPerWord) {
	  phaseraw_word = *phaseraw_iter++;
	  read_idx_lowbits = 0;
	}
	const uintptr_t new_phasepresent_word_mask = new_phasepresent_word - k1LU;
	const uintptr_t reduced_phasepresent_word = new_phasepresent_word & new_phasepresent_word_mask;
	new_phaseinfo_word += (new_phasepresent_word - reduced_phasepresent_word) * (phaseraw_word & 1);
	++read_idx_lowbits;
	phaseraw_word >>= 1;
	new_phasepresent_word = reduced_phasepresent_word;
      }
      phaseinfo[widx] = new_phaseinfo_word;
    }
  } else {
    phaseraw_word >>= 1;
    uintptr_t* phasepresent = *phasepresent_ptr;
    const uintptr_t* phaseinfo_read_iter = &(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]);
    uintptr_t phaseinfo_read_word = *phaseinfo_read_iter++;
    uint32_t phaseinfo_read_idx_lowbits = 0;
    for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
      uintptr_t cur_hets = all_hets[widx];
      uintptr_t new_phasepresent = 0;
      uintptr_t new_phaseinfo = 0;
      while (cur_hets) {
	if (read_idx_lowbits == kBitsPerWord) {
	  phaseraw_word = *phaseraw_iter++;
	  read_idx_lowbits = 0;
	}
	const uintptr_t cur_hets_mask = cur_hets - k1LU;
	const uintptr_t reduced_hets = cur_hets & cur_hets_mask;
	if (phaseraw_word & 1) {
	  if (phaseinfo_read_idx_lowbits == kBitsPerWord) {
	    phaseinfo_read_word = *phaseinfo_read_iter++;
	    phaseinfo_read_idx_lowbits = 0;
	  }
	  const uintptr_t cur_bit = cur_hets - reduced_hets;
	  new_phasepresent += cur_bit;
	  new_phaseinfo += cur_bit * (phaseinfo_read_word & 1);
	  ++phaseinfo_read_idx_lowbits;
	  phaseinfo_read_word >>= 1;
	}
	++read_idx_lowbits;
	phaseraw_word >>= 1;
	cur_hets = reduced_hets;
      }
      phasepresent[widx] = new_phasepresent;
      phaseinfo[widx] = new_phaseinfo;
    }
  }
}

void unpack_hphase_subset(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uintptr_t* phaseraw_iter = phaseraw;
  uintptr_t* phaseinfo_write_iter = phaseinfo;
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  uintptr_t phaseraw_word = *phaseraw_iter++;
  uintptr_t phaseinfo_write_word = 0;
  uint32_t read_idx_lowbits = 1;
  uint32_t write_idx_lowbits = 0;
  if (!(phaseraw_word & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
      const uintptr_t cur_sample_include = sample_include[widx];
      const uintptr_t geno_hets = all_hets[widx];
      uintptr_t tmp_phaseinfo_write_word = 0;
      if (geno_hets) {
	const uint32_t read_idx_lowbits_end = read_idx_lowbits + popcount_long(geno_hets);
	uintptr_t tmp_phaseinfo_input_word = phaseraw_word >> read_idx_lowbits;
	if (read_idx_lowbits_end >= kBitsPerWord) {
	  // always safe to read an extra word off the end, when
	  // read_idx_lowbits_end == kBitsPerWord and we're at the last word
	  phaseraw_word = *phaseraw_iter++;
	  if (read_idx_lowbits) {
	    tmp_phaseinfo_input_word |= phaseraw_word << (kBitsPerWord - read_idx_lowbits);
	  }
	}
	tmp_phaseinfo_input_word &= (~k0LU) >> (kBitsPerWord + read_idx_lowbits - read_idx_lowbits_end);
	read_idx_lowbits = read_idx_lowbits_end % kBitsPerWord;
	if (tmp_phaseinfo_input_word) {
	  uintptr_t cur_masked_hets = cur_sample_include & geno_hets;
	  while (cur_masked_hets) {
	    const uintptr_t cur_masked_hets_and_arg = cur_masked_hets - k1LU;
	    const uintptr_t lowmask = (cur_masked_hets ^ cur_masked_hets_and_arg) >> 1;
	    const uint32_t read_idx_offset = popcount_long(geno_hets & lowmask);
	    uintptr_t shifted_phaseinfo_input_word = tmp_phaseinfo_input_word >> read_idx_offset;
	    if (shifted_phaseinfo_input_word & 1) {
	      tmp_phaseinfo_write_word |= (k1LU << popcount_long(cur_sample_include & lowmask));
	      if (shifted_phaseinfo_input_word == 1) {
		break;
	      }
	    }
	    cur_masked_hets &= cur_masked_hets_and_arg;
	  }
	}
        phaseinfo_write_word |= tmp_phaseinfo_write_word << write_idx_lowbits;
      }
      const uint32_t write_idx_lowbits_end = write_idx_lowbits + popcount_long(cur_sample_include);
      if (write_idx_lowbits_end >= kBitsPerWord) {
	*phaseinfo_write_iter++ = phaseinfo_write_word;
	if (write_idx_lowbits) {
	  phaseinfo_write_word = tmp_phaseinfo_write_word >> (kBitsPerWord - write_idx_lowbits);
	} else {
	  phaseinfo_write_word = 0;
	}
      }
      write_idx_lowbits = write_idx_lowbits_end % kBitsPerWord;
    }
    if (write_idx_lowbits) {
      *phaseinfo_write_iter = phaseinfo_write_word;
    }
    return;
  }
  const uintptr_t* phaseinfo_read_iter = &(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]);
  uintptr_t* phasepresent_write_iter = *phasepresent_ptr;
  uintptr_t phaseinfo_read_word = *phaseinfo_read_iter++;
  uintptr_t phasepresent_write_word = 0;
  uint32_t phaseinfo_read_idx_lowbits = 0;
  for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
    const uintptr_t cur_sample_include = sample_include[widx];
    const uintptr_t geno_hets = all_hets[widx];
    uintptr_t tmp_phasepresent_write_word = 0;
    uintptr_t tmp_phaseinfo_write_word = 0;
    if (geno_hets) {
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + popcount_long(geno_hets);
      uintptr_t tmp_phasepresent_input_word = phaseraw_word >> read_idx_lowbits;
      if (read_idx_lowbits_end >= kBitsPerWord) {
	// always safe to read an extra word off the end, when
	// read_idx_lowbits_end == kBitsPerWord and we're at the last word
	phaseraw_word = *phaseraw_iter++;
	if (read_idx_lowbits) {
	  tmp_phasepresent_input_word |= phaseraw_word << (kBitsPerWord - read_idx_lowbits);
	}
      }
      tmp_phasepresent_input_word &= (~k0LU) >> (kBitsPerWord + read_idx_lowbits - read_idx_lowbits_end);
      read_idx_lowbits = read_idx_lowbits_end % kBitsPerWord;
      if (tmp_phasepresent_input_word) {
	const uint32_t read_phasepresent_ct = popcount_long(tmp_phasepresent_input_word);
	uintptr_t tmp_phaseinfo_input_word;
	// avoid reading off end of phaseinfo here
	if (phaseinfo_read_idx_lowbits != kBitsPerWord) {
	  const uint32_t phaseinfo_read_idx_lowbits_end = phaseinfo_read_idx_lowbits + read_phasepresent_ct;
	  tmp_phaseinfo_input_word = phaseinfo_read_word >> phaseinfo_read_idx_lowbits;
	  if (phaseinfo_read_idx_lowbits_end < kBitsPerWord) {
	    phaseinfo_read_idx_lowbits = phaseinfo_read_idx_lowbits_end;
	  } else {
	    phaseinfo_read_word = *phaseinfo_read_iter++;
	    tmp_phaseinfo_input_word |= phaseinfo_read_word << (kBitsPerWord - phaseinfo_read_idx_lowbits);
	    phaseinfo_read_idx_lowbits = phaseinfo_read_idx_lowbits_end - kBitsPerWord;
	  }
	} else {
	  // special case, can't right-shift 64
	  phaseinfo_read_word = *phaseinfo_read_iter++;
	  phaseinfo_read_idx_lowbits = read_phasepresent_ct;
	  tmp_phaseinfo_input_word = phaseinfo_read_word;
	}
	tmp_phaseinfo_input_word &= (~k0LU) >> (kBitsPerWord - read_phasepresent_ct);
	
	uintptr_t cur_masked_hets = cur_sample_include & geno_hets;
	while (cur_masked_hets) {
	  const uintptr_t cur_masked_hets_and_arg = cur_masked_hets - k1LU;
	  const uintptr_t lowmask = (cur_masked_hets ^ cur_masked_hets_and_arg) >> 1;
	  const uint32_t read_idx_offset = popcount_long(geno_hets & lowmask);
	  uintptr_t shifted_phasepresent_input_word = tmp_phasepresent_input_word >> read_idx_offset;
	  if (shifted_phasepresent_input_word & 1) {
	    const uintptr_t cur_bit = popcount_long(cur_sample_include & lowmask);
	    tmp_phasepresent_write_word |= cur_bit;
	    tmp_phaseinfo_write_word += cur_bit * ((tmp_phaseinfo_input_word >> (read_phasepresent_ct - popcount_long(shifted_phasepresent_input_word))) & 1);
	    if (shifted_phasepresent_input_word == 1) {
	      break;
	    }
	  }
	  cur_masked_hets &= cur_masked_hets_and_arg;
	}
      }
      phasepresent_write_word |= tmp_phasepresent_write_word << write_idx_lowbits;
      phaseinfo_write_word |= tmp_phaseinfo_write_word << write_idx_lowbits;
    }
    const uint32_t write_idx_lowbits_end = write_idx_lowbits + popcount_long(cur_sample_include);
    if (write_idx_lowbits_end >= kBitsPerWord) {
      *phasepresent_write_iter++ = phasepresent_write_word;
      *phaseinfo_write_iter++ = phaseinfo_write_word;
      if (write_idx_lowbits) {
	const uint32_t rshift = kBitsPerWord - write_idx_lowbits;
	phasepresent_write_word = tmp_phasepresent_write_word >> rshift;
	phaseinfo_write_word = tmp_phaseinfo_write_word >> rshift;
      } else {
	phasepresent_write_word = 0;
	phaseinfo_write_word = 0;
      }
    }
    write_idx_lowbits = write_idx_lowbits_end % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    *phasepresent_write_iter = phasepresent_write_word;
    *phaseinfo_write_iter = phaseinfo_write_word;
  }
}

void unpack_and_resort_hphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* sample_include, const uint32_t* old_sample_idx_to_new, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uintptr_t* phaseraw_iter = phaseraw;
  const uint32_t* old_sample_idx_to_new_iter = old_sample_idx_to_new;
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t phaseraw_word = *phaseraw_iter++;
  uint32_t read_idx_lowbits = 1;
  fill_ulong_zero(sample_ctl, phaseinfo);
  if (!(phaseraw_word & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
      uintptr_t new_phasepresent_word = all_hets[widx];
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + popcount_long(new_phasepresent_word);
      uintptr_t tmp_phaseinfo_input_word = phaseraw_word >> read_idx_lowbits;
      if (read_idx_lowbits_end >= kBitsPerWord) {
	// always safe to read an extra word off the end
	phaseraw_word = *phaseraw_iter++;
	if (read_idx_lowbits) {
	  tmp_phaseinfo_input_word |= phaseraw_word << (kBitsPerWord - read_idx_lowbits);
	}
      }
      // no need to mask off top bits of tmp_phaseinfo_input_word
      read_idx_lowbits = read_idx_lowbits_end % kBitsPerWord;
      if (!sample_include) {
	while (new_phasepresent_word) {
	  const uint32_t sample_uidx_lowbits = CTZLU(new_phasepresent_word);
	  if (tmp_phaseinfo_input_word & 1) {
	    SET_BIT(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
	  }
	  tmp_phaseinfo_input_word >>= 1;
	  new_phasepresent_word &= new_phasepresent_word - k1LU;
	}
      } else {
	uintptr_t masked_phasepresent_word = new_phasepresent_word & sample_include[widx];
	while (masked_phasepresent_word) {
	  const uint32_t sample_uidx_lowbits = CTZLU(masked_phasepresent_word);
	  const uintptr_t lowmask = (k1LU << sample_uidx_lowbits) - k1LU;
	  if ((tmp_phaseinfo_input_word >> popcount_long(new_phasepresent_word & lowmask)) & 1) {
	    SET_BIT(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
	  }
	  masked_phasepresent_word &= masked_phasepresent_word - k1LU;
	}
      }
      old_sample_idx_to_new_iter = &(old_sample_idx_to_new_iter[kBitsPerWord]);
    }
    return;
  }
  uintptr_t* phasepresent = *phasepresent_ptr;
  const uintptr_t* phaseinfo_read_iter = &(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]);
  uintptr_t phaseinfo_read_word = *phaseinfo_read_iter++;
  uint32_t phaseinfo_read_idx_lowbits = 0;
  fill_ulong_zero(sample_ctl, phasepresent);
  for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
    uintptr_t geno_hets = all_hets[widx];
    if (geno_hets) {
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + popcount_long(geno_hets);
      uintptr_t tmp_phasepresent_input_word = phaseraw_word >> read_idx_lowbits;
      if (read_idx_lowbits_end >= kBitsPerWord) {
	// always safe to read an extra word off the end, when
	// read_idx_lowbits_end == kBitsPerWord and we're at the last word
	phaseraw_word = *phaseraw_iter++;
	if (read_idx_lowbits) {
	  tmp_phasepresent_input_word |= phaseraw_word << (kBitsPerWord - read_idx_lowbits);
	}
      }
      tmp_phasepresent_input_word &= (~k0LU) >> (kBitsPerWord + read_idx_lowbits - read_idx_lowbits_end);
      read_idx_lowbits = read_idx_lowbits_end % kBitsPerWord;
      if (tmp_phasepresent_input_word) {
	const uint32_t read_phasepresent_ct = popcount_long(tmp_phasepresent_input_word);
	uintptr_t tmp_phaseinfo_input_word;
	// avoid reading off end of phaseinfo here
	if (phaseinfo_read_idx_lowbits != kBitsPerWord) {
	  const uint32_t phaseinfo_read_idx_lowbits_end = phaseinfo_read_idx_lowbits + read_phasepresent_ct;
	  tmp_phaseinfo_input_word = phaseinfo_read_word >> phaseinfo_read_idx_lowbits;
	  if (phaseinfo_read_idx_lowbits_end < kBitsPerWord) {
	    phaseinfo_read_idx_lowbits = phaseinfo_read_idx_lowbits_end;
	  } else {
	    phaseinfo_read_word = *phaseinfo_read_iter++;
	    tmp_phaseinfo_input_word |= phaseinfo_read_word << (kBitsPerWord - phaseinfo_read_idx_lowbits);
	    phaseinfo_read_idx_lowbits = phaseinfo_read_idx_lowbits_end - kBitsPerWord;
	  }
	} else {
	  // special case, can't right-shift 64
	  phaseinfo_read_word = *phaseinfo_read_iter++;
	  phaseinfo_read_idx_lowbits = read_phasepresent_ct;
	  tmp_phaseinfo_input_word = phaseinfo_read_word;
	}
	// no need to mask off top bits of tmp_phaseinfo_input_word
	if (!sample_include) {
	  while (1) {
	    if (tmp_phasepresent_input_word & 1) {
	      const uint32_t new_sample_idx = old_sample_idx_to_new_iter[CTZLU(geno_hets)];
	      SET_BIT(new_sample_idx, phasepresent);
	      if (tmp_phaseinfo_input_word & 1) {
		SET_BIT(new_sample_idx, phaseinfo);
	      }
	      if (tmp_phasepresent_input_word == 1) {
		break;
	      }
	      tmp_phaseinfo_input_word >>= 1;
	    }
	    tmp_phasepresent_input_word >>= 1;
	    geno_hets &= geno_hets - k1LU;
	  }
	} else {
	  const uintptr_t sample_include_word = sample_include[widx];
	  while (1) {
	    if (tmp_phasepresent_input_word & 1) {
	      const uint32_t sample_uidx_lowbits = CTZLU(geno_hets);
	      if ((sample_include_word >> sample_uidx_lowbits) & 1) {
		const uint32_t new_sample_idx = old_sample_idx_to_new_iter[sample_uidx_lowbits];
		SET_BIT(new_sample_idx, phasepresent);
		if (tmp_phaseinfo_input_word & 1) {
		  SET_BIT(new_sample_idx, phaseinfo);
		}
	      }
	      if (tmp_phasepresent_input_word == 1) {
		break;
	      }
	      tmp_phaseinfo_input_word >>= 1;
	    }
	    tmp_phasepresent_input_word >>= 1;
	    geno_hets &= geno_hets - k1LU;
	  }
	}
      }
    }
    old_sample_idx_to_new_iter = &(old_sample_idx_to_new_iter[kBitsPerWord]);
  }
}

void copy_dosage(const uintptr_t* __restrict dosageraw, uint32_t raw_sample_ct, uintptr_t* __restrict write_dosagepresent, dosage_t* write_dosagevals, uint32_t* write_dosage_ct_ptr) {
  const uint32_t raw_sample_ctaw = BITCT_TO_ALIGNED_WORDCT(raw_sample_ct);
  const uintptr_t* read_dosagepresent = dosageraw;
  const dosage_t* read_dosagevals = (const dosage_t*)(&(dosageraw[raw_sample_ctaw]));
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  const uint32_t dosage_ct = popcount_longs(read_dosagepresent, raw_sample_ctl);
  *write_dosage_ct_ptr = dosage_ct;
  memcpy(write_dosagepresent, read_dosagepresent, raw_sample_ctl * sizeof(intptr_t));
  memcpy(write_dosagevals, read_dosagevals, dosage_ct * sizeof(dosage_t));
}

void copy_dosage_subset(const uintptr_t* __restrict dosageraw, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* __restrict write_dosagepresent, dosage_t* write_dosagevals, uint32_t* __restrict write_dosage_ct_ptr) {
  const uint32_t raw_sample_ctaw = BITCT_TO_ALIGNED_WORDCT(raw_sample_ct);
  const uintptr_t* read_dosagepresent = dosageraw;
  const dosage_t* read_dosagevals = (const dosage_t*)(&(dosageraw[raw_sample_ctaw]));
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  const uint32_t read_dosage_ct = popcount_longs(read_dosagepresent, raw_sample_ctl);
  copy_bitarr_subset(read_dosagepresent, sample_include, sample_ct, write_dosagepresent);
  uint32_t sample_uidx = 0;
  dosage_t* write_dosagevals_iter = write_dosagevals;
  for (uint32_t read_dosage_idx = 0; read_dosage_idx < read_dosage_ct; ++read_dosage_idx, ++sample_uidx) {
    next_set_unsafe_ck(read_dosagepresent, &sample_uidx);
    if (is_set(sample_include, sample_uidx)) {
      *write_dosagevals_iter++ = read_dosagevals[read_dosage_idx];
    }
  }
  *write_dosage_ct_ptr = (uintptr_t)(write_dosagevals_iter - write_dosagevals);
}

void copy_and_resort_dosage(const uintptr_t* __restrict dosageraw, const uint32_t* new_sample_idx_to_old, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* __restrict write_dosagepresent, dosage_t* write_dosagevals, uint32_t* write_dosage_ct_ptr, uint32_t* cumulative_popcount_buf) {
  const uint32_t raw_sample_ctaw = BITCT_TO_ALIGNED_WORDCT(raw_sample_ct);
  const uintptr_t* read_dosagepresent = dosageraw;
  const dosage_t* read_dosagevals = (const dosage_t*)(&(dosageraw[raw_sample_ctaw]));
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  fill_cumulative_popcounts(read_dosagepresent, raw_sample_ctl, cumulative_popcount_buf);
  const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  fill_ulong_zero(sample_ctl, write_dosagepresent);
  dosage_t* write_dosagevals_iter = write_dosagevals;
  for (uint32_t new_sample_idx = 0; new_sample_idx < sample_ct; ++new_sample_idx) {
    const uint32_t old_sample_idx = new_sample_idx_to_old[new_sample_idx];
    if (is_set(read_dosagepresent, old_sample_idx)) {
      set_bit(new_sample_idx, write_dosagepresent);
      const uint32_t old_dosagevals_idx = raw_to_subsetted_pos(read_dosagepresent, cumulative_popcount_buf, old_sample_idx);
      *write_dosagevals_iter++ = read_dosagevals[old_dosagevals_idx];
    }
  }
  *write_dosage_ct_ptr = (uintptr_t)(write_dosagevals_iter - write_dosagevals);
}


// more multithread globals
static pgen_reader_t** g_pgr_ptrs = nullptr;
static uintptr_t** g_genovecs = nullptr;
static uintptr_t** g_dosage_presents = nullptr;
static dosage_t** g_dosage_val_bufs = nullptr;
static uint32_t* g_read_variant_uidx_starts = nullptr; // size calc_thread_ct

static uint64_t* g_allele_dosages = nullptr;
static uint32_t* g_raw_geno_cts = nullptr;
static uint32_t* g_variant_missing_hc_cts = nullptr;
static uint32_t* g_variant_missing_dosage_cts = nullptr;
static uint32_t* g_variant_hethap_cts = nullptr;
static uint64_t* g_founder_allele_dosages = nullptr;
static uint32_t* g_founder_raw_geno_cts = nullptr;
static uint32_t* g_x_male_geno_cts = nullptr;
static uint32_t* g_founder_x_male_geno_cts = nullptr;
static uint32_t* g_x_nosex_geno_cts = nullptr;
static uint32_t* g_founder_x_nosex_geno_cts = nullptr;
static double* g_mach_r2_vals = nullptr;

static unsigned char* g_writebufs[2] = {nullptr, nullptr};

static const uintptr_t* g_variant_include = nullptr;
static const chr_info_t* g_cip = nullptr;
static const uintptr_t* g_sample_include = nullptr;
static uintptr_t* g_sample_include_interleaved_vec = nullptr;
static uint32_t* g_sample_include_cumulative_popcounts = nullptr;
static uintptr_t* g_sex_male = nullptr;
static uintptr_t* g_sex_male_interleaved_vec = nullptr;
static uintptr_t* g_sex_male_collapsed_interleaved = nullptr;
static uintptr_t* g_sex_female_collapsed_interleaved = nullptr;
static uint32_t* g_sex_male_cumulative_popcounts = nullptr;
static uintptr_t* g_nosex_interleaved_vec = nullptr;
static const uintptr_t* g_founder_info = nullptr;
static uintptr_t* g_founder_info_interleaved_vec = nullptr;
static uint32_t* g_founder_info_cumulative_popcounts = nullptr;
static uintptr_t* g_founder_male = nullptr;
static uintptr_t* g_founder_male_interleaved_vec = nullptr;
static uint32_t* g_founder_male_cumulative_popcounts = nullptr;
static uintptr_t* g_founder_nosex_interleaved_vec = nullptr;
static const uintptr_t* g_variant_allele_idxs = nullptr;
static const alt_allele_ct_t* g_refalt1_select = nullptr;
static const uint32_t* g_collapsed_sort_map = nullptr;
static const uint32_t* g_new_sample_idx_to_old = nullptr;
static uint32_t* g_old_sample_idx_to_new = nullptr;
static uint32_t g_raw_sample_ct = 0;
// g_sample_ct, g_calc_thread_ct, g_cur_block_write_ct, g_hard_call_halfdist,
// g_dosage_erase_halfdist, g_error_ret declared earlier
static uint32_t g_founder_ct = 0;
static uint32_t g_male_ct = 0;
static uint32_t g_nosex_ct = 0;
static uint32_t g_founder_male_ct = 0;
static uint32_t g_founder_nosex_ct = 0;
static uint32_t g_first_hap_uidx = 0;
static pgen_global_flags_t g_read_phase_dosage_gflags = kfPgenGlobal0;

// just store the beginning of each vblock for now
// (may want to store record lengths later)
static uintptr_t** g_loadbuf_thread_starts[2] = {nullptr, nullptr};

// phase, dosage
static unsigned char* g_loaded_vrtypes[2] = {nullptr, nullptr};

static vul_t** g_thread_vecaligned_bufs = nullptr;
static uintptr_t** g_thread_write_genovecs = nullptr;
static uintptr_t** g_thread_write_phasepresents = nullptr;
static uintptr_t** g_thread_write_phaseinfos = nullptr;
static uintptr_t** g_thread_all_hets = nullptr;
static uintptr_t** g_thread_write_dosagepresents = nullptr;
static dosage_t** g_thread_write_dosagevals = nullptr;
static uint32_t** g_thread_cumulative_popcount_bufs = nullptr;
static pgen_writer_common_t** g_pwcs = nullptr;

static uintptr_t* g_plink1_smaj_loadbuf_iter = nullptr;
static uint32_t g_stride = 0;

THREAD_FUNC_DECL plink1_smaj_transpose_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t write_batch_ct_m1 = (sample_ct - 1) / kPglQuaterTransposeBatch;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  pgen_writer_common_t* pwcp = g_pwcs[tidx];
  vul_t* vecaligned_buf = g_thread_vecaligned_bufs[tidx];
  uintptr_t* write_genovec = g_thread_write_genovecs[tidx];
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    const uint32_t loadbuf_ul_stride = g_stride;
    uint32_t write_idx = tidx * kPglVblockSize;
    uintptr_t* read_iter = &(g_plink1_smaj_loadbuf_iter[write_idx / kBitsPerWordD2]);
    const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
    while (write_idx < write_idx_end) {
      const uintptr_t* read_iter2 = read_iter;
      // uintptr_t* write_iter = write_genovec;
      const uint32_t vblock_size = MINV(kPglQuaterTransposeBatch, write_idx_end - write_idx);
      uint32_t write_batch_idx = 0;
      uint32_t read_batch_size = kPglQuaterTransposeBatch;
      while (1) {
	if (write_batch_idx >= write_batch_ct_m1) {
	  if (write_batch_idx > write_batch_ct_m1) {
	    break;
	  }
	  read_batch_size = MOD_NZ(sample_ct, kPglQuaterTransposeBatch);
	}
	transpose_quaterblock(read_iter2, loadbuf_ul_stride, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[write_batch_idx * kPglQuaterTransposeWords]), vecaligned_buf);
	read_iter2 = &(read_iter2[kPglQuaterTransposeBatch * loadbuf_ul_stride]);
	++write_batch_idx;
      }
      for (uint32_t uii = 0; uii < vblock_size; ++uii) {
	uintptr_t* cur_write_genovec = &(write_genovec[uii * sample_ctaw2]);
	pgr_plink1_to_plink2_inplace_unsafe(sample_ct, cur_write_genovec);
	zero_trailing_quaters(sample_ct, cur_write_genovec);
	pwc_append_biallelic_genovec(cur_write_genovec, pwcp);
      }
      write_idx += vblock_size;
      read_iter = &(read_iter[kPglQuaterTransposeWords]);
    }
    if ((tidx == calc_thread_ct - 1) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

pglerr_t plink1_sample_major_to_pgen(const char* pgenname, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile) {
  unsigned char* bigstack_mark = g_bigstack_base;
  mt_pgen_writer_t* mpgwp = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    // file size already validated by pgfi_init_phase1()
    LOGPRINTFWW("Sample-major .bed file detected.  Transposing to %s .\n", pgenname);
    fputs("0%", stdout);
    fflush(stdout);
    if ((!variant_ct) || (!sample_ct)) {
      // todo: hardcoded 12-byte write
      logprint("\n");
      logerrprint("Error: Zero-variant/zero-sample .pgen writing is not currently supported.\n");
      reterr = kPglRetNotYetSupported;
      goto plink1_sample_major_to_pgen_ret_1;
    }
    const uint32_t variant_ct4 = QUATERCT_TO_BYTECT(variant_ct);
    unsigned char* raw_loadbuf = nullptr;
    uint32_t raw_load_batch_size = 1;
    if (variant_ct4 < 5120) {
      // assuming 4K block size, fseek won't let us avoid reading many
      // unnecessary disk blocks
      raw_load_batch_size += 131071 / variant_ct4;
      if (bigstack_alloc_uc(raw_load_batch_size * variant_ct4, &raw_loadbuf)) {
	goto plink1_sample_major_to_pgen_ret_NOMEM;
      }
    }
    const uint32_t raw_load_batch_ct_m1 = (sample_ct - 1) / raw_load_batch_size;
    if (!raw_load_batch_ct_m1) {
      raw_load_batch_size = sample_ct;
    }
    const uint32_t raw_load_batch_ct = raw_load_batch_ct_m1 + 1;
    uintptr_t alloc_base_cacheline_ct;
    uint64_t mpgw_per_thread_cacheline_ct;
    uint32_t vrec_len_byte_ct;
    uint64_t vblock_cacheline_ct;
    mpgw_init_phase1(nullptr, variant_ct, sample_ct, kfPgenGlobal0, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);
#ifndef __LP64__
    if ((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (vblock_cacheline_ct > (0x7fffffff / kCacheline))) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
#endif
    
    uint32_t calc_thread_ct = DIV_UP(variant_ct, kPglVblockSize);
    if (calc_thread_ct >= max_thread_ct) {
      calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    }
    mpgwp = (mt_pgen_writer_t*)bigstack_alloc((calc_thread_ct + DIV_UP(sizeof(mt_pgen_writer_t), kBytesPerWord)) * sizeof(intptr_t));
    if (!mpgwp) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
    mpgwp->pgen_outfile = nullptr;
    pthread_t* threads;
    if (bigstack_alloc_thread(calc_thread_ct, &threads) ||
	bigstack_alloc_vp(calc_thread_ct, &g_thread_vecaligned_bufs) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_genovecs)) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
    g_pwcs = &(mpgwp->pwcs[0]);
    uintptr_t cachelines_avail = bigstack_left() / kCacheline;
    // inner loop transposes kPglQuaterTransposeBatch variants at a time
    const uintptr_t transpose_thread_cacheline_ct = kPglQuaterTransposeBufbytes / kCacheline + QUATERCT_TO_VECCT(sample_ct) * (kPglQuaterTransposeBatch / kVecsPerCacheline);
    if (cachelines_avail < calc_thread_ct * ((uint64_t)transpose_thread_cacheline_ct)) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_thread_vecaligned_bufs[tidx] = (vul_t*)bigstack_alloc_raw(kPglQuaterTransposeBufbytes);
      g_thread_write_genovecs[tidx] = (uintptr_t*)bigstack_alloc_raw(QUATERCT_TO_VECCT(sample_ct) * kBytesPerVec * kPglQuaterTransposeBatch);
    }
    cachelines_avail = bigstack_left() / kCacheline;
    // Main workflow:
    // 1. Load next calc_thread_ct * load_multiplier * kPglVblockSize
    //    variants.
    //    calc_thread_ct is reduced as necessary to ensure the compression
    //    write buffers use <= 1/8 of total workspace.
    //    with calc_thread_ct determined, load_multiplier is then chosen to use
    //    as much of the remaining workspace as possible.
    // 2. Repeat load_multiplier times:
    //    a. Spawn threads processing calc_thread_ct vblocks
    //    b. Join threads
    //    c. Flush results
    // 3. Goto step 1 unless eof.  (load_multiplier may be smaller on last
    //    iteration.)
    // No double-buffering here since main bottleneck is how many variants we
    // can load at once.
    if ((cachelines_avail / 8) < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) {
      if ((cachelines_avail / 8) < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct) {
	// possible todo: simple single-threaded fallback
	goto plink1_sample_major_to_pgen_ret_NOMEM;
      }
      calc_thread_ct = ((cachelines_avail / 8) - alloc_base_cacheline_ct) / mpgw_per_thread_cacheline_ct;
    }
    // todo: determine appropriate calc_thread_ct limit.  (should not be less
    // than 7-8.)
    unsigned char* mpgw_alloc = bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline);
    reterr = mpgw_init_phase2(pgenname, nullptr, nullptr, variant_ct, sample_ct, kfPgenGlobal0, 2 - real_ref_alleles, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);
    
    cachelines_avail = bigstack_left() / kCacheline;
    const uint64_t full_load_vecs_req = sample_ct * ((uint64_t)QUATERCT_TO_ALIGNED_WORDCT(variant_ct));
    uintptr_t* plink1_smaj_loadbuf;
    uint32_t load_multiplier;
    uint32_t cur_vidx_ct;
    if (full_load_vecs_req > cachelines_avail * kVecsPerCacheline) {
      // each iteration requires ((kPglVblockSize / 4) * calc_thread_ct *
      //   sample_ct) bytes to be loaded
      load_multiplier = cachelines_avail / ((kPglVblockSize / (4 * kCacheline)) * calc_thread_ct * ((uintptr_t)sample_ct));
      assert(load_multiplier);
      cur_vidx_ct = load_multiplier * calc_thread_ct * kPglVblockSize;
      plink1_smaj_loadbuf = (uintptr_t*)bigstack_alloc_raw_rd((cur_vidx_ct / 4) * ((uintptr_t)sample_ct));
    } else {
      load_multiplier = 1 + ((variant_ct - 1) / (calc_thread_ct * kPglVblockSize));
      cur_vidx_ct = variant_ct;
      plink1_smaj_loadbuf = (uintptr_t*)bigstack_alloc_raw_rd(full_load_vecs_req * kBytesPerVec);
    }
    uint32_t cur_vidx_base = 0;
    uint32_t cur_vidx_ct4 = QUATERCT_TO_BYTECT(cur_vidx_ct);
    uint32_t cur_vidx_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(cur_vidx_ct);
    uint32_t pass_idx = 0;
    const uint32_t pass_ct = 1 + (variant_ct - 1) / cur_vidx_ct;
    g_sample_ct = sample_ct;
    g_stride = QUATERCT_TO_VECCT(cur_vidx_ct) * kWordsPerVec;
    g_calc_thread_ct = calc_thread_ct;
    while (1) {
      uint32_t raw_load_batch_idx = 0;
      uint32_t cur_raw_load_batch_size = raw_load_batch_size;
      uintptr_t* smaj_loadbuf_iter = plink1_smaj_loadbuf;
      ++pass_idx;
      putc_unlocked('\r', stdout);
      printf("Pass %u/%u: loading... 0%%", pass_idx, pass_ct);
      fflush(stdout);
      uint32_t pct = 0;
      uint32_t next_print_idx = raw_load_batch_ct / 100;
      const uint64_t seek_addl_offset = 3 + cur_vidx_base / 4;
      while (1) {
	if (raw_load_batch_size == 1) {
	  if (fseeko(infile, seek_addl_offset + raw_load_batch_idx * ((uint64_t)variant_ct4), SEEK_SET)) {
	    goto plink1_sample_major_to_pgen_ret_READ_FAIL;
	  }
	  if (!fread(smaj_loadbuf_iter, cur_vidx_ct4, 1, infile)) {
	    goto plink1_sample_major_to_pgen_ret_READ_FAIL;
	  }
	  smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
	} else {
	  if (!fread(raw_loadbuf, cur_raw_load_batch_size * variant_ct4, 1, infile)) {
	    goto plink1_sample_major_to_pgen_ret_READ_FAIL;
	  }
	  unsigned char* raw_loadbuf_iter = &(raw_loadbuf[cur_vidx_base / 4]);
	  for (uint32_t uii = 0; uii < cur_raw_load_batch_size; ++uii) {
	    memcpy(smaj_loadbuf_iter, raw_loadbuf_iter, cur_vidx_ct4);
	    raw_loadbuf_iter = &(raw_loadbuf_iter[variant_ct4]);
	    smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
	  }
	}
	++raw_load_batch_idx;
	if (raw_load_batch_idx >= raw_load_batch_ct_m1) {
	  if (raw_load_batch_idx > raw_load_batch_ct_m1) {
	    break;
	  }
	  cur_raw_load_batch_size = sample_ct - raw_load_batch_idx * raw_load_batch_size;
	}
	if (raw_load_batch_idx >= next_print_idx) {
	  if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
	  pct = (raw_load_batch_idx * 100LLU) / raw_load_batch_ct;
	  printf("\b\b%u%%", pct++);
	  fflush(stdout);
	  next_print_idx = (pct * ((uint64_t)raw_load_batch_ct)) / 100;
	}
      }
      const uintptr_t last_tidx = calc_thread_ct - 1;
      uint32_t load_idx = 0;
      g_cur_block_write_ct = calc_thread_ct * kPglVblockSize;
      uint32_t is_last_block;
      putc_unlocked('\r', stdout);
      printf("Pass %u/%u: transposing and compressing... 0%%", pass_idx, pass_ct);
      pct = 0;
      next_print_idx = load_idx / 100;
      do {
	if (load_idx >= next_print_idx) {
	  if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
	  pct = (load_idx * 100LLU) / load_multiplier;
	  printf("\b\b%u%%", pct++);
	  fflush(stdout);
	  next_print_idx = (pct * ((uint64_t)load_multiplier)) / 100;
	}
	g_plink1_smaj_loadbuf_iter = &(plink1_smaj_loadbuf[load_idx * calc_thread_ct * (kPglVblockSize / kBitsPerWordD2)]);
	is_last_block = (++load_idx == load_multiplier);
	if (is_last_block) {
	  g_cur_block_write_ct = cur_vidx_ct - (load_idx - 1) * calc_thread_ct * kPglVblockSize;
	}
	if (last_tidx) {
	  if (spawn_threads2z(plink1_smaj_transpose_thread, last_tidx, is_last_block, threads)) {
	    goto plink1_sample_major_to_pgen_ret_THREAD_CREATE_FAIL;
	  }
	}
	plink1_smaj_transpose_thread((uintptr_t*)last_tidx);
	if (last_tidx) {
	  join_threads2z(last_tidx, is_last_block, threads);
	}
	reterr = mpgw_flush(mpgwp);
	if (reterr) {
	  if (!is_last_block) {
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(plink1_smaj_transpose_thread, last_tidx, threads);
	  }
	  goto plink1_sample_major_to_pgen_ret_WRITE_FAIL;
	}
      } while (!is_last_block);
      cur_vidx_base += cur_vidx_ct;
      if (cur_vidx_base == variant_ct) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	break;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                     ", stdout);
      // assumes pgfi_init_phase1() leaves file pointer at byte 3; otherwise,
      // necessary to put this at top of main loop
      if (fseeko(infile, 3, SEEK_SET)) {
	goto plink1_sample_major_to_pgen_ret_READ_FAIL;
      }
      if (variant_ct - cur_vidx_base <= cur_vidx_ct) {
	cur_vidx_ct = variant_ct - cur_vidx_base;
	cur_vidx_ct4 = QUATERCT_TO_BYTECT(cur_vidx_ct);
        cur_vidx_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(cur_vidx_ct);
        g_stride = QUATERCT_TO_VECCT(cur_vidx_ct) * kWordsPerVec;
        load_multiplier = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
      }
    }
    mpgwp = nullptr;
    fputs("\b\bdone.\n", stdout);
    LOGPRINTF("Transpose complete.\n");
  }
  while (0) {
  plink1_sample_major_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  plink1_sample_major_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  plink1_sample_major_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  plink1_sample_major_to_pgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 plink1_sample_major_to_pgen_ret_1:
  if (mpgw_cleanup(mpgwp) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

THREAD_FUNC_DECL load_allele_and_geno_counts_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const chr_info_t* cip = g_cip;
  const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t subset_ct = (g_founder_info != nullptr) + 1;
  const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  const uint32_t first_hap_uidx = g_first_hap_uidx;
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  dosage_t* dosage_vals = nullptr;
  if (g_dosage_presents) {
    dosage_present = g_dosage_presents[tidx];
    dosage_vals = g_dosage_val_bufs[tidx];
  }
  uint32_t is_y = 0;
  uint32_t is_nonxy_haploid = 0;
  uint32_t x_start = 0;
  int32_t x_code;
  if (xymt_exists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
    x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
  }
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    // no overflow danger since cur_block_write_ct <= 2^16, tidx < (2^16 - 1)
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    const uintptr_t* sample_include = g_sample_include;
    const uintptr_t* sample_include_interleaved_vec = g_sample_include_interleaved_vec;
    const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
    const uintptr_t* sex_male = g_sex_male;
    const uintptr_t* sex_male_interleaved_vec = g_sex_male_interleaved_vec;
    const uint32_t* sex_male_cumulative_popcounts = g_sex_male_cumulative_popcounts;
    const uintptr_t* nosex_interleaved_vec = g_nosex_interleaved_vec;
    uint32_t sample_ct = g_sample_ct;
    uint32_t male_ct = g_male_ct;
    uint32_t nosex_ct = g_nosex_ct;
    uint64_t* allele_dosages = g_allele_dosages;
    uint32_t* raw_geno_cts = g_raw_geno_cts;
    uint32_t* variant_missing_hc_cts = g_variant_missing_hc_cts;
    uint32_t* variant_missing_dosage_cts = g_variant_missing_dosage_cts;
    uint32_t* variant_hethap_cts = g_variant_hethap_cts;
    uint32_t* x_male_geno_cts = g_x_male_geno_cts;
    uint32_t* x_nosex_geno_cts = g_x_nosex_geno_cts;
    double* mach_r2_vals = g_mach_r2_vals;
    uint32_t subset_idx = 0;
    uint32_t dosage_ct = 0;
    while (1) {
      uint32_t cur_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
      uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
      uint32_t chr_end = 0;
      uint32_t is_x_or_y = 0;
      pglerr_t reterr = kPglRetSuccess;
      
      // different approach will be needed for multiallelic case...
      uint32_t genocounts[4];
      uint32_t sex_specific_genocounts[4];
      for (; cur_idx < cur_idx_end; ++cur_idx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	if (variant_uidx >= chr_end) {
	  const uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx);
	  const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	  is_y = 0;
	  is_nonxy_haploid = 0;
	  if (chr_idx == x_code) {
	    is_x_or_y = 1;
	    pgr_clear_ld_cache(pgrp);
	  } else if (chr_idx == y_code) {
	    is_x_or_y = 1;
	    is_y = 1;
	    pgr_clear_ld_cache(pgrp);
	  } else {
	    if (is_x_or_y) {
	      pgr_clear_ld_cache(pgrp);
	    }
	    is_x_or_y = 0;
	    // no way for this to happen now unless everything is haploid?
	    is_nonxy_haploid = is_set(cip->haploid_mask, chr_idx);
	  }
	}
	const uintptr_t cur_variant_allele_idx = variant_allele_idxs? variant_allele_idxs[variant_uidx] : (2 * variant_uidx);
	// insert cur_allele_ct == 2 check here
	uint64_t cur_dosages[2];
	uint32_t hethap_ct;
	if (!is_x_or_y) {
	  // call pgr_get_refalt1_genotype_counts() instead when dosages not
	  // needed?
	  reterr = pgr_get_ref_nonref_genotype_counts_and_dosage16s(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, mach_r2_vals? (&(mach_r2_vals[variant_uidx])) : nullptr, genocounts, cur_dosages);
	  if (reterr) {
	    g_error_ret = reterr;
	    break;
	  }
	  if (!is_nonxy_haploid) {
	    // in multiallelic case, check ref vs. non-ref...	    
	    hethap_ct = 0;
	    if (allele_dosages) {
	      // ...but save all allele counts here.
	      // workhorse multiallelic count function should return both
	      // individual allele counts and (at least if appropriate
	      // parameter is not nullptr) hom-altx total.
	      allele_dosages[cur_variant_allele_idx] = cur_dosages[0] * 2;
	      allele_dosages[cur_variant_allele_idx + 1] = cur_dosages[1] * 2;
	    }
	  } else {
	    hethap_ct = genocounts[1];
	    if (allele_dosages) {
	      allele_dosages[cur_variant_allele_idx] = cur_dosages[0];
	      allele_dosages[cur_variant_allele_idx + 1] = cur_dosages[1];
	    }
	  }
	} else if (is_y) {
	  reterr = pgr_get_ref_nonref_genotype_counts_and_dosage16s(sex_male, sex_male_interleaved_vec, sex_male_cumulative_popcounts, male_ct, variant_uidx, pgrp, nullptr, genocounts, cur_dosages);
	  if (reterr) {
	    g_error_ret = reterr;
	    break;
	  }
	  hethap_ct = genocounts[1];
	  if (allele_dosages) {
	    allele_dosages[cur_variant_allele_idx] = cur_dosages[0];
	    allele_dosages[cur_variant_allele_idx + 1] = cur_dosages[1];
	  }
	} else {
	  // chrX
	  uint32_t is_explicit_alt1;
	  reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(nullptr, nullptr, raw_sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
	  if (reterr) {
	    g_error_ret = reterr;
	    break;
	  }
	  // assert(!is_explicit_alt1);
	  if (sample_ct == raw_sample_ct) {
	    zero_trailing_quaters(raw_sample_ct, genovec);
	    genovec_count_freqs_unsafe(genovec, sample_ct, genocounts);
	  } else {
	    genovec_count_subset_freqs(genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
	  }
	  genovec_count_subset_freqs(genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, sex_specific_genocounts);
	  hethap_ct = sex_specific_genocounts[1];
	  if (allele_dosages) {
	    uint32_t sample_uidx = 0;
	    uintptr_t replaced_ct = 0; // nonmales count twice
	    uintptr_t alt1_ct = 4 * genocounts[2] + 2 * genocounts[1] - 2 * sex_specific_genocounts[2] - hethap_ct; // nonmales count twice
	    uint64_t alt1_dosage = 0; // in 32768ths, nonmales count twice
	    uint32_t included_dosage_ct = 0; // nonmales count twice
	    if (sample_ct == raw_sample_ct) {
	      for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
		next_set_unsafe_ck(dosage_present, &sample_uidx);
		const uintptr_t cur_dosage_val = dosage_vals[dosage_idx];
		const uintptr_t sex_multiplier = 2 - IS_SET(sex_male, sample_uidx);
		alt1_dosage += cur_dosage_val * sex_multiplier;

		// could call genoarr_count_subset_intersect_freqs() twice
		// instead, but since we've already manually extracted the sex
		// bit it probably doesn't help?
		const uintptr_t hardcall_code = GET_QUATERARR_ENTRY(genovec, sample_uidx);
		if (hardcall_code != 3) {
		  alt1_ct -= hardcall_code * sex_multiplier;
		  replaced_ct += sex_multiplier;
		}
	      }
	      included_dosage_ct = 2 * dosage_ct;
	      if (dosage_ct) {
		included_dosage_ct -= popcount_longs_intersect(dosage_present, sex_male, raw_sample_ctl);
	      }
	    } else {
	      for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
		next_set_unsafe_ck(dosage_present, &sample_uidx);
		if (IS_SET(sample_include, sample_uidx)) {
		  const uintptr_t cur_dosage_val = dosage_vals[dosage_idx];
		  const uintptr_t sex_multiplier = 2 - IS_SET(sex_male, sample_uidx);
		  alt1_dosage += cur_dosage_val * sex_multiplier;
		  included_dosage_ct += sex_multiplier;
		  const uintptr_t hardcall_code = GET_QUATERARR_ENTRY(genovec, sample_uidx);
		  if (hardcall_code != 3) {
		    alt1_ct -= hardcall_code * sex_multiplier;
		    replaced_ct += sex_multiplier;
		  }
		}
	      }
	    }
	    const uintptr_t ref_ct = (2 * (sample_ct - genocounts[3]) - male_ct + sex_specific_genocounts[3]) * (2 * k1LU) - alt1_ct;
	    allele_dosages[cur_variant_allele_idx] = included_dosage_ct * ((uint64_t)kDosageMax) - alt1_dosage + (ref_ct * ((uint64_t)kDosageMid));
	    allele_dosages[cur_variant_allele_idx + 1] = alt1_dosage + (alt1_ct * ((uint64_t)kDosageMid));
	  }
	  if (x_male_geno_cts) {
	    uint32_t* cur_x_male_geno_cts = &(x_male_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	    cur_x_male_geno_cts[0] = sex_specific_genocounts[0];
	    cur_x_male_geno_cts[1] = sex_specific_genocounts[1];
	    cur_x_male_geno_cts[2] = sex_specific_genocounts[2];
	    if (x_nosex_geno_cts) {
	      genovec_count_subset_freqs(genovec, nosex_interleaved_vec, raw_sample_ct, nosex_ct, sex_specific_genocounts);
	      uint32_t* cur_nosex_geno_cts = &(x_nosex_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
	      cur_nosex_geno_cts[0] = sex_specific_genocounts[0];
	      cur_nosex_geno_cts[1] = sex_specific_genocounts[1];
	      cur_nosex_geno_cts[2] = sex_specific_genocounts[2];
	    }
	  }
	}
	if (raw_geno_cts) {
	  uint32_t* cur_raw_geno_cts = &(raw_geno_cts[(3 * k1LU) * variant_uidx]);
	  cur_raw_geno_cts[0] = genocounts[0];
	  cur_raw_geno_cts[1] = genocounts[1];
	  cur_raw_geno_cts[2] = genocounts[2];
	}
	if (variant_missing_hc_cts) {
	  variant_missing_hc_cts[variant_uidx] = genocounts[3];
	  if (variant_hethap_cts && (variant_uidx >= first_hap_uidx)) {
	    variant_hethap_cts[variant_uidx - first_hap_uidx] = hethap_ct;
	  }
	}
	if (variant_missing_dosage_cts) {
	  uint32_t missing_dosage_ct;
	  if (!is_x_or_y) {
	    missing_dosage_ct = sample_ct - ((cur_dosages[0] + cur_dosages[1]) / kDosageMax);
	  } else if (is_y) {
	    missing_dosage_ct = male_ct - ((cur_dosages[0] + cur_dosages[1]) / kDosageMax);
	  } else {
	    if (dosage_ct) {
	      zero_trailing_quaters(raw_sample_ct, genovec);
	      missing_dosage_ct = genoarr_count_missing_notsubset_unsafe(genovec, dosage_present, raw_sample_ct);
	    } else {
	      missing_dosage_ct = genocounts[3];
	    }
	  }
	  variant_missing_dosage_cts[variant_uidx] = missing_dosage_ct;
	}
      }
      if ((++subset_idx == subset_ct) || reterr) {
	break;
      }
      sample_include = g_founder_info;
      sample_include_interleaved_vec = g_founder_info_interleaved_vec;
      sample_include_cumulative_popcounts = g_founder_info_cumulative_popcounts;
      sex_male = g_founder_male;
      sex_male_interleaved_vec = g_founder_male_interleaved_vec;
      sex_male_cumulative_popcounts = g_founder_male_cumulative_popcounts;

      nosex_interleaved_vec = g_founder_nosex_interleaved_vec;
      
      sample_ct = g_founder_ct;
      male_ct = g_founder_male_ct;
      nosex_ct = g_founder_nosex_ct;
      allele_dosages = g_founder_allele_dosages;
      variant_missing_hc_cts = nullptr;
      variant_missing_dosage_cts = nullptr;
      raw_geno_cts = g_founder_raw_geno_cts;
      x_male_geno_cts = g_founder_x_male_geno_cts;
      x_nosex_geno_cts = g_founder_x_nosex_geno_cts;
      mach_r2_vals = nullptr;
      pgr_clear_ld_cache(pgrp);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

pglerr_t load_allele_and_geno_counts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uintptr_t* variant_allele_idxs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint64_t* allele_dosages, uint64_t* founder_allele_dosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, uint32_t* raw_geno_cts, uint32_t* founder_raw_geno_cts, uint32_t* x_male_geno_cts, uint32_t* founder_x_male_geno_cts, uint32_t* x_nosex_geno_cts, uint32_t* founder_x_nosex_geno_cts, double* mach_r2_vals) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!variant_ct) {
      goto load_allele_and_geno_counts_ret_1;
    }
    if (variant_allele_idxs) {
      logerrprint("Error: load_allele_and_geno_counts() doesn't support multiallelic variants yet.\n");
      reterr = kPglRetNotYetSupported;
      goto load_allele_and_geno_counts_ret_1;
    }

    // four cases:
    // 1. allele_dosages, raw_geno_cts, and/or variant_missing_{hc,dosage}_cts
    //    required, and that's it
    // 2. founder_allele_dosages and/or founder_raw_geno_cts required, and
    //    that's it
    // 3. both required, and founder_ct != sample_ct.
    // 4. both required, and founder_ct == sample_ct.  caller is expected to
    //    make founder_allele_dosages and allele_dosages point to the same
    //    memory, ditto for founder_raw_geno_cts/raw_geno_cts.
    const uint32_t only_founder_cts_required = (!allele_dosages) && (!raw_geno_cts) && (!variant_missing_hc_cts) && (!variant_missing_dosage_cts);
    const uint32_t two_subsets_required = (founder_ct != sample_ct) && (!only_founder_cts_required) && (founder_allele_dosages || founder_raw_geno_cts);
    g_cip = cip;
    g_sample_include = only_founder_cts_required? founder_info : sample_include;
    g_raw_sample_ct = raw_sample_ct;
    g_sample_ct = only_founder_cts_required? founder_ct : sample_ct;
    g_male_ct = male_ct;
    g_allele_dosages = only_founder_cts_required? founder_allele_dosages : allele_dosages;
    g_raw_geno_cts = only_founder_cts_required? founder_raw_geno_cts : raw_geno_cts;
    g_x_male_geno_cts = only_founder_cts_required? founder_x_male_geno_cts : x_male_geno_cts;
    g_x_nosex_geno_cts = only_founder_cts_required? founder_x_nosex_geno_cts : x_nosex_geno_cts;
    g_mach_r2_vals = mach_r2_vals;
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t raw_sample_ctv = BITCT_TO_VECCT(raw_sample_ct);
    if (bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_sample_include_interleaved_vec) ||
	bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_cumulative_popcounts) ||
	bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_sex_male_interleaved_vec) ||
	bigstack_alloc_ui(raw_sample_ctl, &g_sex_male_cumulative_popcounts)) {
      goto load_allele_and_geno_counts_ret_NOMEM;
    }
    fill_interleaved_mask_vec(g_sample_include, raw_sample_ctv, g_sample_include_interleaved_vec);
    fill_cumulative_popcounts(g_sample_include, raw_sample_ctl, g_sample_include_cumulative_popcounts);
    if ((founder_ct == sample_ct) || (!only_founder_cts_required)) {
      // const_cast
      g_sex_male = (uintptr_t*)((uintptr_t)sex_male);
    } else {
      // no nonfounder counts required
      if (bigstack_alloc_ul(raw_sample_ctl, &g_sex_male)) {
	goto load_allele_and_geno_counts_ret_NOMEM;
      }
      for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
	g_sex_male[widx] = sex_male[widx] & founder_info[widx];
      }
    }
    fill_interleaved_mask_vec(g_sex_male, raw_sample_ctv, g_sex_male_interleaved_vec);
    fill_cumulative_popcounts(g_sex_male, raw_sample_ctl, g_sex_male_cumulative_popcounts);
    if (!(x_nosex_geno_cts || founder_x_nosex_geno_cts)) {
      nosex_ct = 0;
    }
    g_nosex_ct = nosex_ct;
    g_nosex_interleaved_vec = nullptr;
    uintptr_t* nosex_buf = nullptr;
    if (nosex_ct) {
      if (bigstack_end_alloc_ul(raw_sample_ctl, &nosex_buf) ||
          bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_nosex_interleaved_vec)) {
	goto load_allele_and_geno_counts_ret_NOMEM;
      }
      for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
	nosex_buf[widx] = (~sex_nm[widx]) & g_sample_include[widx];
      }
      fill_interleaved_mask_vec(nosex_buf, raw_sample_ctv, g_nosex_interleaved_vec);
    }

    g_variant_missing_hc_cts = variant_missing_hc_cts;
    g_variant_missing_dosage_cts = variant_missing_dosage_cts;
    g_variant_hethap_cts = variant_hethap_cts;
    g_first_hap_uidx = first_hap_uidx;
    
    g_founder_info = nullptr;
    g_founder_info_interleaved_vec = nullptr;
    g_founder_info_cumulative_popcounts = nullptr;
    g_founder_male = nullptr;
    g_founder_male_interleaved_vec = nullptr;
    g_founder_male_cumulative_popcounts = nullptr;
    g_founder_nosex_interleaved_vec = nullptr;
    g_founder_ct = 0;
    g_founder_male_ct = 0;
    g_founder_nosex_ct = 0;
    g_founder_allele_dosages = nullptr;
    g_founder_raw_geno_cts = nullptr;
    g_founder_x_male_geno_cts = nullptr;
    g_founder_x_nosex_geno_cts = nullptr;
    if (two_subsets_required) {
      if (founder_ct) {
	g_founder_info = founder_info;
	if (bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_founder_info_interleaved_vec) ||
	    bigstack_alloc_ui(raw_sample_ctl, &g_founder_info_cumulative_popcounts) ||
	    bigstack_alloc_ul(raw_sample_ctl, &g_founder_male) ||
	    bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_founder_male_interleaved_vec) ||
	    bigstack_alloc_ui(raw_sample_ctl, &g_founder_male_cumulative_popcounts)) {
	  goto load_allele_and_geno_counts_ret_NOMEM;
	}
	fill_interleaved_mask_vec(founder_info, raw_sample_ctv, g_founder_info_interleaved_vec);
	fill_cumulative_popcounts(founder_info, raw_sample_ctl, g_founder_info_cumulative_popcounts);
	for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
	  g_founder_male[widx] = sex_male[widx] & founder_info[widx];
	}
	fill_interleaved_mask_vec(g_founder_male, raw_sample_ctv, g_founder_male_interleaved_vec);
	fill_cumulative_popcounts(g_founder_male, raw_sample_ctl, g_founder_male_cumulative_popcounts);
	g_founder_ct = founder_ct;
	g_founder_male_ct = g_founder_male_cumulative_popcounts[raw_sample_ctl - 1] + popcount_long(g_founder_male[raw_sample_ctl - 1]);
	g_founder_allele_dosages = founder_allele_dosages;
	g_founder_raw_geno_cts = founder_raw_geno_cts;
	g_founder_x_male_geno_cts = founder_x_male_geno_cts;
	if (nosex_ct) {
	  // caller currently responsible for ensuring that when
	  // founder_nosex_ct is zero, founder_x_nosex_geno_cts ==
	  // nullptr
	  if (bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_founder_nosex_interleaved_vec)) {
	    goto load_allele_and_geno_counts_ret_NOMEM;
	  }
	  for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
	    nosex_buf[widx] &= founder_info[widx];
	  }
	  g_founder_nosex_ct = popcount_longs(nosex_buf, raw_sample_ctl);
	  assert(g_founder_nosex_ct);
          fill_interleaved_mask_vec(nosex_buf, raw_sample_ctv, g_founder_nosex_interleaved_vec);
	  g_founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
	}
      } else {
	if (founder_allele_dosages) {
	  fill_ull_zero(variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct), founder_allele_dosages);
	}
	if (founder_raw_geno_cts) {
	  fill_uint_zero((3 * k1LU) * raw_variant_ct, founder_raw_geno_cts);
	}
      }
    } else if (founder_ct == sample_ct) {
      // bugfix: some founder and some nonfounder counts required
      if ((!g_allele_dosages) && founder_allele_dosages) {
	g_allele_dosages = founder_allele_dosages;
      }
      if ((!g_raw_geno_cts) && founder_raw_geno_cts) {
	g_raw_geno_cts = founder_raw_geno_cts;
      }
      if ((!g_x_male_geno_cts) && founder_x_male_geno_cts) {
	g_x_male_geno_cts = founder_x_male_geno_cts;
      }
      if ((!g_x_nosex_geno_cts) && founder_x_nosex_geno_cts) {
	g_x_nosex_geno_cts = founder_x_nosex_geno_cts;
      }
    } else if (only_founder_cts_required) {
      g_male_ct = g_sex_male_cumulative_popcounts[raw_sample_ctl - 1] + popcount_long(g_sex_male[raw_sample_ctl - 1]);
      if (nosex_ct) {
        g_nosex_ct = popcount_longs(nosex_buf, raw_sample_ctl);
      }
    }
    if (!g_sample_ct) {
      if (g_allele_dosages) {
	fill_ull_zero(variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct), g_allele_dosages);
      }
      if (g_raw_geno_cts) {
	fill_uint_zero((3 * k1LU) * raw_variant_ct, g_raw_geno_cts);
      }
      // early exit
      goto load_allele_and_geno_counts_ret_1;
    }
    bigstack_end_reset(bigstack_end_mark); // free nosex_buf

    int32_t ii;
    const uint32_t x_dosages_needed = (allele_dosages || founder_allele_dosages || variant_missing_dosage_cts) && xymt_exists(cip, kChrOffsetX, &ii) && (pgfip->gflags & kfPgenGlobalDosagePresent);
    if (!x_dosages_needed) {
      // defensive
      g_dosage_presents = nullptr;
      g_dosage_val_bufs = nullptr;
    }
    
    // todo: check when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    // todo: check if raw_sample_ct should be replaced with sample_ct here
    if (multithread_load_init(variant_include, raw_sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, x_dosages_needed? (&g_dosage_presents) : nullptr, x_dosages_needed? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto load_allele_and_geno_counts_ret_NOMEM;
    }

    g_variant_include = variant_include;
    g_variant_allele_idxs = variant_allele_idxs;
    g_calc_thread_ct = calc_thread_ct;
    g_error_ret = kPglRetSuccess;

    logprint("Calculating allele frequencies... ");
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
      uintptr_t cur_block_write_ct = 0;
      if (!is_last_block) {
	while (read_block_idx < read_block_ct_m1) {
	  // this uses multithread_load_init's guarantee that read_block_size
	  // is either raw_variant_ct or a multiple of kBitsPerVec
	  cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	  if (cur_block_write_ct) {
	    break;
	  }
	  ++read_block_idx;
	}
	if (read_block_idx == read_block_ct_m1) {
	  cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	  cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	}
	if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
	  if (variant_idx) {
	    join_threads2z(calc_thread_ct, 0, threads);
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(load_allele_and_geno_counts_thread, calc_thread_ct, threads);
	  }
	  goto load_allele_and_geno_counts_ret_READ_FAIL;
	}
      }
      if (variant_idx) {
	join_threads2z(calc_thread_ct, is_last_block, threads);
	reterr = g_error_ret;
	if (reterr) {
	  if (!is_last_block) {
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(load_allele_and_geno_counts_thread, calc_thread_ct, threads);
	  }
	  if (reterr == kPglRetMalformedInput) {
	    logprint("\n");
	    logerrprint("Error: Malformed .pgen file.\n");
	  }
	  goto load_allele_and_geno_counts_ret_1;
	}
      }
      if (!is_last_block) {
	g_cur_block_write_ct = cur_block_write_ct;
	compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	  g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	}
	is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
	if (spawn_threads2z(load_allele_and_geno_counts_thread, calc_thread_ct, is_last_block, threads)) {
	  goto load_allele_and_geno_counts_ret_THREAD_CREATE_FAIL;
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
      variant_idx += cur_block_write_ct;
      // crucially, this is independent of the pgen_reader_t block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  load_allele_and_geno_counts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_allele_and_geno_counts_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_allele_and_geno_counts_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 load_allele_and_geno_counts_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

FLAGSET_DEF_START()
  kfPlink2Write0,
  kfPlink2WriteSetHhMissing = (1 << 0),
  kfPlink2WriteSetMixedMtMissing = (1 << 1),
  kfPlink2WriteMeMissing = (1 << 2),
  kfPlink2WriteZeroCluster = (1 << 3),
  kfPlink2WriteFillRef = (1 << 4),
  kfPlink2WriteLateDosageErase = (1 << 5),
  // no need for sample_sort, determined by g_collapsed_sort_map != nullptr?
  kfPlink2WritePlink1 = (1 << 6)
FLAGSET_DEF_END(plink2_write_flags_t);
// todo: add .pgen-specific stuff

static plink2_write_flags_t g_plink2_write_flags = kfPlink2Write0;

THREAD_FUNC_DECL make_bedlike_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const chr_info_t* cip = g_cip;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed_interleaved = g_sex_female_collapsed_interleaved;
  const uint32_t* collapsed_sort_map = g_collapsed_sort_map;
  const uint32_t set_hh_missing = g_plink2_write_flags & kfPlink2WriteSetHhMissing;
  const uint32_t set_mixed_mt_missing = g_plink2_write_flags & kfPlink2WriteSetMixedMtMissing;
  const uint32_t write_plink1 = g_plink2_write_flags & kfPlink2WritePlink1;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  const uint32_t sample_ctv2 = QUATERCT_TO_VECCT(sample_ct);
  const uint32_t sample_ct4 = QUATERCT_TO_BYTECT(sample_ct);
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(g_writebufs[parity][write_idx * sample_ct4]);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    uint32_t chr_end = 0;
    uint32_t is_x_or_y = 0;
    uint32_t is_y = 0;
    uint32_t is_haploid = 0;
    uint32_t is_mt = 0;
    for (; write_idx < write_idx_end; ++write_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	const uint32_t chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	is_y = (chr_idx == y_code);
	is_x_or_y = is_y || (chr_idx == x_code);
	is_haploid = is_set(cip->haploid_mask, chr_idx);
	is_mt = (chr_idx == mt_code);
      }
      const pglerr_t reterr = pgr_read_refalt1_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec);
      if (reterr) {
	// printf("fail vidx: %u\n", variant_uidx);
	g_error_ret = reterr;
	break;
      }
      // this doesn't work in multiallelic case
      // todo: pgenlib_internal function which takes two allele indexes
      if (refalt1_select && (refalt1_select[2 * variant_uidx] == 1)) {
	genovec_invert_unsafe(sample_ct, genovec);
      }
      if (set_hh_missing && is_haploid) {
	if (is_x_or_y) {
	  // male hets to missing
	  set_male_het_missing(sex_male_collapsed_interleaved, sample_ctv2, genovec);
	  if (is_y) {
	    // all female calls to missing; unknown-sex calls now left alone
	    interleaved_set_missing(sex_female_collapsed_interleaved, sample_ctv2, genovec);
	  }
	} else {
	  // all hets to missing
	  set_het_missing(sample_ctl2, genovec);
	}
      } else if (set_mixed_mt_missing && is_mt) {
	// all hets to missing
	set_het_missing(sample_ctl2, genovec);
      }
      // todo: --set-me-missing, --zero-cluster, --fill-missing-with-ref
      // (--set-me-missing should happen after --set-hh-missing)
      if (write_plink1) {
	pgr_plink2_to_plink1_inplace_unsafe(sample_ct, genovec);
      }
      // trailing bytes don't matter, but trailing bits of last byte may
      zero_trailing_quaters(sample_ct, genovec);
      if (!collapsed_sort_map) {
	writebuf_iter = (unsigned char*)memcpya(writebuf_iter, genovec, sample_ct4);
      } else {
	genovec_resort(genovec, collapsed_sort_map, sample_ct, writebuf_iter);
	writebuf_iter = &(writebuf_iter[sample_ct4]);
      }
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

// combine existing chr_mask/xymt_codes/haploid_mask/chr_idx_to_foidx with new
// collapsed chromosome boundary table
static uint32_t* g_write_chr_fo_vidx_start = nullptr;

static st_pgen_writer_t* g_spgwp = nullptr;

THREAD_FUNC_DECL make_pgen_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t* new_sample_idx_to_old = g_new_sample_idx_to_old;
  const uint32_t* old_sample_idx_to_new = g_old_sample_idx_to_new;
  const chr_info_t* cip = g_cip;
  const uint32_t* write_chr_fo_vidx_start = g_write_chr_fo_vidx_start;
  const alt_allele_ct_t* refalt1_select_iter = g_refalt1_select;
  const uintptr_t* sample_include = g_sample_include;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed_interleaved = g_sex_female_collapsed_interleaved;
  const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  const uint32_t sample_ctv2 = QUATERCT_TO_VECCT(sample_ct);
  const uint32_t raw_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(raw_sample_ct);
  const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];

  const uint32_t set_hh_missing = g_plink2_write_flags & kfPlink2WriteSetHhMissing;
  const uint32_t set_mixed_mt_missing = g_plink2_write_flags & kfPlink2WriteSetMixedMtMissing;
  const uint32_t late_dosage_erase = g_plink2_write_flags & kfPlink2WriteLateDosageErase;
  
  const uint32_t hphase_present = (g_read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
  const uint32_t dosage_present = (g_read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uintptr_t phaseraw_word_ct = kWordsPerVec + round_down_pow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
  // todo: double dosage_vals allocation in phased-dosage case
  const uintptr_t dosageraw_word_ct = kWordsPerVec * (BITCT_TO_VECCT(raw_sample_ct) + DIV_UP(raw_sample_ct, (kBytesPerVec / sizeof(dosage_t))));
  
  st_pgen_writer_t* spgwp = g_spgwp;
  pgen_writer_common_t* pwcp;
  if (spgwp) {
    pwcp = &(spgwp->pwc);
  } else {
    pwcp = g_pwcs[tidx];
  }
  uintptr_t* write_genovec = nullptr;
  // assumes g_sample_include == nullptr if sample_ct == raw_sample_ct
  if (new_sample_idx_to_old || sample_include) {
    write_genovec = g_thread_write_genovecs[tidx];
    write_genovec[sample_ctl2 - 1] = 0;
  }
  uintptr_t* write_phasepresent = nullptr;
  uintptr_t* write_phaseinfo = nullptr;
  uintptr_t* all_hets = nullptr;
  if (hphase_present) {
    all_hets = g_thread_all_hets[tidx];
    write_phasepresent = g_thread_write_phasepresents[tidx];
    write_phaseinfo = g_thread_write_phaseinfos[tidx];
  }
  uintptr_t* write_dosagepresent = nullptr;
  dosage_t* write_dosagevals = nullptr;
  uint32_t* cumulative_popcount_buf = nullptr;
  if (dosage_present) {
    write_dosagepresent = g_thread_write_dosagepresents[tidx];
    write_dosagevals = g_thread_write_dosagevals[tidx];
    if (new_sample_idx_to_old) {
      cumulative_popcount_buf = g_thread_cumulative_popcount_bufs[tidx];
    }
  }
  uint32_t variant_idx_offset = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t write_idx = tidx * kPglVblockSize;
    const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
    uintptr_t* loadbuf_iter = g_loadbuf_thread_starts[parity][tidx];
    unsigned char* loaded_vrtypes = g_loaded_vrtypes[parity];
    uint32_t loaded_vrtype = 0;
    uint32_t chr_end_bidx = 0;
    uint32_t is_x_or_y = 0;
    uint32_t is_y = 0;
    uint32_t is_haploid = 0;
    uint32_t is_mt = 0;
    for (; write_idx < write_idx_end; ++write_idx) {
      if (loaded_vrtypes) {
	loaded_vrtype = loaded_vrtypes[write_idx];
      }
      if (write_idx >= chr_end_bidx) {
	const uint32_t chr_fo_idx = uint32arr_greater_than(&(write_chr_fo_vidx_start[1]), cip->chr_ct, write_idx + variant_idx_offset + 1);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	chr_end_bidx = write_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
	is_y = (chr_idx == y_code);
	is_x_or_y = is_y || (chr_idx == x_code);
	is_haploid = is_set(cip->haploid_mask, chr_idx);
	is_mt = (chr_idx == mt_code);
      }
      uint32_t is_hphase = loaded_vrtype & 0x10;
      const uint32_t is_dosage = loaded_vrtype & 0x60;
      uintptr_t* cur_write_phasepresent = write_phasepresent;
      if (1) {
	// biallelic, no phased-dosage
	uintptr_t* cur_genovec_end = &(loadbuf_iter[raw_sample_ctaw2]);
	uintptr_t* cur_phaseraw = nullptr;
	uintptr_t* cur_dosageraw = nullptr;
	if (is_hphase) {
	  pgr_detect_genovec_hets(loadbuf_iter, raw_sample_ct, all_hets);
	  cur_phaseraw = cur_genovec_end;
	  cur_genovec_end = &(cur_genovec_end[phaseraw_word_ct]);
	}
	if (is_dosage) {
	  cur_dosageraw = cur_genovec_end;
	  cur_genovec_end = &(cur_genovec_end[dosageraw_word_ct]);
	}
	uint32_t write_dosage_ct = 0;
	if (new_sample_idx_to_old) {
	  genovec_resort(loadbuf_iter, new_sample_idx_to_old, sample_ct, (unsigned char*)write_genovec);
	  if (is_hphase) {
	    unpack_and_resort_hphase(all_hets, cur_phaseraw, sample_include, old_sample_idx_to_new, raw_sample_ct, sample_ct, &cur_write_phasepresent, write_phaseinfo);
	  }
	  if (is_dosage) {
	    copy_and_resort_dosage(cur_dosageraw, new_sample_idx_to_old, raw_sample_ct, sample_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct, cumulative_popcount_buf);
	  }
	} else if (sample_include) {
	  copy_quaterarr_nonempty_subset(loadbuf_iter, sample_include, raw_sample_ct, sample_ct, write_genovec);
	  if (is_hphase) {
	    unpack_hphase_subset(all_hets, cur_phaseraw, sample_include, raw_sample_ct, &cur_write_phasepresent, write_phaseinfo);
	  }
	  if (is_dosage) {
	    copy_dosage_subset(cur_dosageraw, sample_include, raw_sample_ct, sample_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct);
	  }
	} else {
	  write_genovec = loadbuf_iter;
	  if (is_hphase) {
	    unpack_hphase(all_hets, cur_phaseraw, sample_ct, &cur_write_phasepresent, write_phaseinfo);
	  }
	  if (is_dosage) {
	    copy_dosage(cur_dosageraw, sample_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct);
	  }
	}
	if (refalt1_select_iter && (refalt1_select_iter[2 * write_idx] == 1)) {
	  genovec_invert_unsafe(sample_ct, write_genovec);
	  if (is_hphase) {
	    bitarr_invert(sample_ctl, write_phaseinfo);
	  }
	  if (write_dosage_ct) {
	    biallelic_dosage16_invert(write_dosage_ct, write_dosagevals);
	  }
	}
	if (write_dosage_ct) {
	  if (hard_call_halfdist) {
	    if (is_hphase && (!cur_write_phasepresent)) {
	      cur_write_phasepresent = write_phasepresent;
	      // unsafe to just copy all_hets, because we may have resorted
	      pgr_detect_genovec_hets(write_genovec, sample_ct, cur_write_phasepresent);
	    }
	    uint32_t sample_uidx = 0;
	    for (uint32_t dosage_idx = 0; dosage_idx < write_dosage_ct; ++dosage_idx, ++sample_uidx) {
	      next_set_unsafe_ck(write_dosagepresent, &sample_uidx);
	      const uint32_t dosage_int = write_dosagevals[dosage_idx];
	      const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
	      const uint32_t widx = sample_uidx / kBitsPerWordD2;
	      uintptr_t prev_geno_word = write_genovec[widx];
	      const uint32_t shift = (sample_uidx % kBitsPerWordD2) * 2;
	      uintptr_t new_geno;
	      if (halfdist < hard_call_halfdist) {
		new_geno = 3;
	      } else {
		new_geno = (dosage_int + kDosage4th) / kDosageMid;
	      }
	      const uintptr_t prev_geno = (prev_geno_word >> shift) & 3;
	      const uintptr_t geno_xor = new_geno ^ prev_geno;
	      if (geno_xor) {
	        if (is_hphase) {
		  // must erase phase here
		  CLEAR_BIT(sample_uidx, cur_write_phasepresent);
		}
		write_genovec[widx] = prev_geno_word ^ (geno_xor << shift);
	      }
	    }
	  }
	  if (dosage_erase_halfdist < kDosage4th) {
	    uint32_t dosage_read_idx = 0;
	    uint32_t sample_uidx = 0;
	    for (; dosage_read_idx < write_dosage_ct; ++dosage_read_idx, ++sample_uidx) {
	      next_set_unsafe_ck(write_dosagepresent, &sample_uidx);
	      const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
	      const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
	      if (halfdist >= dosage_erase_halfdist) {
		clear_bit(sample_uidx, write_dosagepresent);
		++sample_uidx;
		break;
	      }
	    }
	    uint32_t dosage_write_idx = dosage_read_idx;
	    while (++dosage_read_idx < write_dosage_ct) {
	      next_set_unsafe_ck(write_dosagepresent, &sample_uidx);
	      const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
	      const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
	      if (halfdist < dosage_erase_halfdist) {
		write_dosagevals[dosage_write_idx++] = dosage_int;
	      } else {
		clear_bit(sample_uidx, write_dosagepresent);
	      }
	      ++sample_uidx;
	    }
	    write_dosage_ct = dosage_write_idx;
	  } else if (late_dosage_erase) {
	    write_dosage_ct = 0;
	  }
	}
	// moved after --hard-call-threshold, since it makes sense to
	// immediately erase fresh het haploid calls
	if (set_hh_missing && is_haploid) {
	  if (is_x_or_y) {
	    // male hets to missing
	    set_male_het_missing(sex_male_collapsed_interleaved, sample_ctv2, write_genovec);
	    if (is_y) {
	      // all female calls to missing; unknown-sex calls now left alone
	      interleaved_set_missing(sex_female_collapsed_interleaved, sample_ctv2, write_genovec);
	    }
	    if (is_hphase && cur_write_phasepresent) {
	      mask_genovec_hets_unsafe(write_genovec, sample_ctl2, cur_write_phasepresent);
	    }
	  } else {
	    // all hets to missing
	    // may want to move is_hphase zeroing in front
	    set_het_missing(sample_ctl2, write_genovec);
	    is_hphase = 0;
	  }
	} else if (set_mixed_mt_missing && is_mt) {
	  // all hets to missing
	  set_het_missing(sample_ctl2, write_genovec);
	  is_hphase = 0;
	}
	zero_trailing_quaters(sample_ct, write_genovec);
	// todo: --set-me-missing, --zero-cluster, --fill-missing-with-ref
	if (spgwp) {
	  if (pwcp->fwrite_bufp >= &(pwcp->fwrite_buf[kPglFwriteBlockSize])) {
	    const uintptr_t cur_byte_ct = (uintptr_t)(pwcp->fwrite_bufp - pwcp->fwrite_buf);
	    if (fwrite_checked(pwcp->fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
	      g_error_ret = kPglRetWriteFail;
	      break;
	    }
	    // printf("vblock_fpos_offset: %llu\n", pwcp->vblock_fpos_offset);
	    pwcp->vblock_fpos_offset += cur_byte_ct;
	    // printf("%u %llu\n", write_idx + variant_idx_offset, pwcp->vblock_fpos_offset);
	    pwcp->fwrite_bufp = pwcp->fwrite_buf;
	  }
	}
	if (!is_hphase) {
	  pwc_append_biallelic_genovec_dosage16(write_genovec, write_dosagepresent, write_dosagevals, write_dosage_ct, pwcp);
	} else {
	  pwc_append_biallelic_genovec_hphase_dosage16(write_genovec, cur_write_phasepresent, write_phaseinfo, write_dosagepresent, write_dosagevals, write_dosage_ct, pwcp);
	  cur_write_phasepresent = write_phasepresent;
	}
	loadbuf_iter = cur_genovec_end;
      } else {
        // todo: multiallelic write
	// (some trim-alts logic here)
      }
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
    variant_idx_offset += cur_block_write_ct;
    if (refalt1_select_iter) {
      refalt1_select_iter = &(refalt1_select_iter[2 * cur_block_write_ct]);
    }
  }
}

pgen_global_flags_t gflags_vfilter(const uintptr_t* variant_include, const unsigned char* vrtypes, uint32_t raw_variant_ct, pgen_global_flags_t input_gflags) {
  pgen_global_flags_t read_phase_dosage_gflags = kfPgenGlobal0;
  const uintptr_t* vrtypes_alias_iter = (const uintptr_t*)vrtypes;
  const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
  uint32_t mask_multiply = ((input_gflags & kfPgenGlobalHardcallPhasePresent)? 0x10 : 0) + ((input_gflags & kfPgenGlobalDosagePresent)? 0x60 : 0) + ((input_gflags & kfPgenGlobalDosagePhasePresent)? 0x80 : 0);
  uintptr_t vrtypes_or = 0;
  for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
    uintptr_t cur_variant_include_word = variant_include[widx];
    if (cur_variant_include_word) {
#ifdef __LP64__
      for (uint32_t vi_byte_idx = 0; vi_byte_idx < 8; ++vi_byte_idx) {
	// this operation maps binary hgfedcba to h0000000g0000000f...
	//                                        ^       ^       ^
	//                                        |       |       |
	//                                       56      48      40
	// 1. (cur_variant_include_word & 0xfe) gives us hgfedcb0;
	//    necessary to avoid carryover.
	// 2. multiply by the number with bits 7, 14, 21, ..., 49 set, to
	//    get hgfedcbhgfedcbhgf...
	//        ^       ^       ^
	//        |       |       |
	//       56      48      40
	// 3. mask out all but bits 8, 16, 24, ..., 56
	// todo: test if this actually beats the per-character loop...
	const uintptr_t cur_mask = (((cur_variant_include_word & 0xfe) * 0x2040810204080LLU) & kMask0101) | (cur_variant_include_word & 1);
	vrtypes_or |= (*vrtypes_alias_iter++) & (cur_mask * mask_multiply);
	cur_variant_include_word >>= 8;
      }
#else
      for (uint32_t vi_hexa_idx = 0; vi_hexa_idx < 8; ++vi_hexa_idx) {
	// dcba -> d0000000c0000000b0000000a
	const uintptr_t cur_mask = ((cur_variant_include_word & 0xf) * 0x204081) & kMask0101;
	vrtypes_or |= (*vrtypes_alias_iter++) & (cur_mask * mask_multiply);
	cur_variant_include_word >>= 4;
      }
#endif
      if (vrtypes_or) {
	if (vrtypes_or & 0x10) {
	  read_phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent;
	  mask_multiply -= 0x10;
	}
	if (vrtypes_or & 0x60) {
	  read_phase_dosage_gflags |= kfPgenGlobalDosagePresent;
	  mask_multiply -= 0x60;
	}
	if (vrtypes_or & 0x80) {
	  read_phase_dosage_gflags |= kfPgenGlobalDosagePhasePresent;
	  mask_multiply -= 0x80;
	}
	if (!mask_multiply) {
	  return read_phase_dosage_gflags;
	}
      }
    }
  }
  return read_phase_dosage_gflags;
}

// Single-output-thread implementation.  Allows variants to be unsorted.
// (Note that make_plink2_no_vsort() requires enough memory for 64k * 2
// variants per output thread, due to LD compression.  This is faster in the
// common case, but once you have 150k+ samples with dosage data...)
pglerr_t make_pgen_robust(const uintptr_t* sample_include, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const alt_allele_ct_t* refalt1_select, const uint32_t* new_variant_idx_to_old, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  // variant_uidx_new_to_old[] can be nullptr

  // caller responsible for initializing g_cip (may need to be different from
  // initial cip struct)
  unsigned char* bigstack_mark = g_bigstack_base;
  threads_state_t ts;
  init_threads3z(&ts);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    // g_plink2_write_flags assumed to include --set-hh-missing and
    //   --set-mixed-mt-missing
    // g_sex_{fe}male_collapsed_interleaved assumed to be initialized if
    //   necessary

    if (bigstack_alloc_thread(1, &ts.threads)) {
      goto make_pgen_robust_ret_NOMEM;
    }
    ts.calc_thread_ct = 1;
    g_spgwp = &spgw;
    const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
    g_sample_include = subsetting_required? sample_include : nullptr;
    g_new_sample_idx_to_old = new_sample_idx_to_old;
    g_raw_sample_ct = raw_sample_ct;
    g_sample_ct = sample_ct;
    g_error_ret = kPglRetSuccess;
    const uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
    if ((make_plink2_modifier & kfMakeBed) || ((make_plink2_modifier & (kfMakePgen | (kfMakePgenFormatBase * 3))) == (kfMakePgen | kfMakePgenFormatBase))) {
      // g_calc_thread_ct = 1;
      logerrprint("Error: Fixed-width .bed/.pgen output doesn't support sorting yet.  Generate a\nregular sorted .pgen first, and then reformat it.\n");
      reterr = kPglRetNotYetSupported;
      goto make_pgen_robust_ret_1;
    } else {
      const uint32_t input_biallelic = (!variant_allele_idxs);
      const uintptr_t* write_allele_idx_offsets = nullptr;
      if (!input_biallelic) {
	if ((variant_ct < raw_variant_ct) || new_variant_idx_to_old_iter) {
	  uintptr_t* new_allele_idx_offsets;
	  if (bigstack_alloc_ul(variant_ct + 1, &new_allele_idx_offsets)) {
	    goto make_pgen_robust_ret_NOMEM;
	  }
	  uintptr_t cur_offset = 0;
	  // todo: separate trim-alts case
	  if (!new_variant_idx_to_old_iter) {
	    uint32_t variant_uidx = 0;
	    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	      next_set_unsafe_ck(variant_include, &variant_uidx);
	      new_allele_idx_offsets[variant_idx] = cur_offset;
	      cur_offset += variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx];
	    }
	  } else {
	    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
	      const uint32_t variant_uidx = *new_variant_idx_to_old_iter++;
	      new_allele_idx_offsets[variant_idx] = cur_offset;
	      cur_offset += variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx];
	    }
	    new_variant_idx_to_old_iter = new_variant_idx_to_old;
	  }
	  if (cur_offset != 2 * variant_ct) {
	    new_allele_idx_offsets[variant_ct] = cur_offset;
	    write_allele_idx_offsets = new_allele_idx_offsets;
	    logprint("Error: Multiallelic .pgen write is not yet supported.\n");
	    reterr = kPglRetNotYetSupported;
	    goto make_pgen_robust_ret_1;
	  } else {
	    bigstack_reset(new_allele_idx_offsets);
	  }
	} else {
	  write_allele_idx_offsets = variant_allele_idxs;
	}
      }
      if ((variant_ct == raw_variant_ct) || new_variant_idx_to_old_iter) {
	g_write_chr_fo_vidx_start = g_cip->chr_fo_vidx_start;
      } else {
	if (alloc_and_fill_subset_chr_fo_vidx_start(variant_include, g_cip, &g_write_chr_fo_vidx_start)) {
	  goto make_pgen_robust_ret_NOMEM;
	}
      }
      pgen_global_flags_t read_phase_dosage_gflags = simple_pgrp->fi.gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      if (make_plink2_modifier & kfMakePgenErasePhase) {
	read_phase_dosage_gflags &= ~(kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
      }
      if (make_plink2_modifier & kfMakePgenEraseDosage) {
	if (hard_call_thresh == 0xffffffffU) {
	  read_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
	  g_plink2_write_flags |= kfPlink2WriteLateDosageErase;
	}
      }
      if (read_phase_dosage_gflags && (variant_ct < raw_variant_ct)) {
	read_phase_dosage_gflags = gflags_vfilter(variant_include, simple_pgrp->fi.vrtypes, raw_variant_ct, simple_pgrp->fi.gflags);
      }
      g_read_phase_dosage_gflags = read_phase_dosage_gflags;
      g_hard_call_halfdist = (hard_call_thresh == 0xffffffffU)? 0 : (kDosage4th - hard_call_thresh);
      g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_hphase_present = (read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
      const uint32_t read_dosage_present = (read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
      pgen_global_flags_t write_phase_dosage_gflags = read_phase_dosage_gflags;
      if (g_plink2_write_flags & kfPlink2WriteLateDosageErase) {
	write_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      }

      uint32_t nonref_flags_storage = 3;
      if (!simple_pgrp->fi.nonref_flags) {
	nonref_flags_storage = (simple_pgrp->fi.gflags & kfPgenGlobalAllNonref)? 2 : 1;
      } else if (variant_ct < raw_variant_ct) {
	// todo: check if now constant
      }
      strcpy(outname_end, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = spgw_init_phase1(outname, write_allele_idx_offsets, simple_pgrp->fi.nonref_flags, variant_ct, sample_ct, write_phase_dosage_gflags, nonref_flags_storage, g_spgwp, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (reterr) {
	goto make_pgen_robust_ret_1;
      }
      unsigned char* spgw_alloc;
      if (bigstack_alloc_ulp(1, &(g_loadbuf_thread_starts[0])) ||
	  bigstack_alloc_ulp(1, &(g_loadbuf_thread_starts[1])) ||
	  bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
	goto make_pgen_robust_ret_NOMEM;
      }
      spgw_init_phase2(max_vrec_len, g_spgwp, spgw_alloc);

      const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
      if (new_sample_idx_to_old || subsetting_required) {
	if (bigstack_alloc_ulp(1, &g_thread_write_genovecs)) {
	  goto make_pgen_robust_ret_NOMEM;
	}
	if (read_hphase_present && new_sample_idx_to_old) {
	  if (bigstack_alloc_ui(raw_sample_ct, &g_old_sample_idx_to_new)) {
	    goto make_pgen_robust_ret_NOMEM;
	  }
	  for (uint32_t new_sample_idx = 0; new_sample_idx < sample_ct; ++new_sample_idx) {
	    g_old_sample_idx_to_new[new_sample_idx_to_old[new_sample_idx]] = new_sample_idx;
	  }
	}
	if (input_biallelic) {
	  if (bigstack_alloc_ul(sample_ctl2, &(g_thread_write_genovecs[0]))) {
	    goto make_pgen_robust_ret_NOMEM;
	  }
	} else {
	  if (bigstack_alloc_ul(DIV_UP(2 * sample_ct * sizeof(alt_allele_ct_t), kBytesPerWord), &(g_thread_write_genovecs[0]))) {
	    goto make_pgen_robust_ret_NOMEM;
	  }
	}
      } else {
        g_thread_write_genovecs = nullptr;
      }
      const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
      const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
      if (read_hphase_present) {
	if (bigstack_alloc_ulp(1, &g_thread_write_phasepresents) ||
	    bigstack_alloc_ulp(1, &g_thread_write_phaseinfos) ||
	    bigstack_alloc_ulp(1, &g_thread_all_hets) ||
	    bigstack_alloc_ul(sample_ctl, &(g_thread_write_phasepresents[0])) ||
	    bigstack_alloc_ul(sample_ctl, &(g_thread_write_phaseinfos[0])) ||
	    bigstack_alloc_ul(raw_sample_ctl, &(g_thread_all_hets[0]))) {
	  goto make_pgen_robust_ret_NOMEM;
	}
      }
      if (read_dosage_present) {
	if (bigstack_alloc_dosagep(1, &g_thread_write_dosagevals) ||
	    bigstack_alloc_ulp(1, &g_thread_write_dosagepresents) ||
	    bigstack_alloc_ul(sample_ctl, &(g_thread_write_dosagepresents[0])) ||
	    bigstack_alloc_dosage(sample_ct, &(g_thread_write_dosagevals[0]))) {
	  goto make_pgen_robust_ret_NOMEM;
	}
	if (new_sample_idx_to_old) {
	  if (bigstack_alloc_uip(1, &g_thread_cumulative_popcount_bufs) ||
	      bigstack_alloc_ui(raw_sample_ctl, &(g_thread_cumulative_popcount_bufs[0]))) {
	    goto make_pgen_robust_ret_NOMEM;
	  }
	}
      }
      g_refalt1_select = refalt1_select;
      if (refalt1_select) {
	if (variant_ct < raw_variant_ct) {
	  // might want inner loop to map variant uidx -> idx instead
	  g_refalt1_select = (alt_allele_ct_t*)bigstack_alloc(variant_ct * 2 * sizeof(alt_allele_ct_t));
	  if (!g_refalt1_select) {
	    goto make_pgen_robust_ret_NOMEM;
	  }
	  uint32_t variant_uidx = 0;
	  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	    next_set_unsafe_ck(variant_include, &variant_uidx);
	    // const_cast
	    memcpy((alt_allele_ct_t*)((uintptr_t)(&(g_refalt1_select[2 * variant_idx]))), &(refalt1_select[2 * variant_uidx]), 2 * sizeof(alt_allele_ct_t));
	  }
	} else {
	  assert(!new_variant_idx_to_old_iter);
	}
      }

      const uint32_t raw_sample_ctv2 = QUATERCT_TO_VECCT(raw_sample_ct);
      uintptr_t load_variant_vec_ct = raw_sample_ctv2;
      uint32_t loaded_vrtypes_needed = 0;
      if (read_hphase_present || read_dosage_present) {
        loaded_vrtypes_needed = 1;
	if (read_hphase_present) {
	  // phaseraw has two parts:
	  // 1. vec-aligned bitarray of up to (raw_sample_ct + 1) bits.  first
	  //    bit is set iff phasepresent is explicitly stored at all (if not,
	  //    all hets are assumed to be phased), if yes the remaining bits
	  //    store packed phasepresent values for all hets, if no the
	  //    remaining bits store packed phaseinfo values for all hets.
	  // 2. word-aligned bitarray of up to raw_sample_ct bits, storing
	  //    phaseinfo values.  (end of this array is vec-aligned.)
	  const uintptr_t phaseraw_word_ct = kWordsPerVec + round_down_pow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
	  load_variant_vec_ct += WORDCT_TO_VECCT(phaseraw_word_ct);
	}
	if (read_dosage_present) {	
	  // todo: phased dosage

	  // (unphased, biallelic) dosageraw has two parts:
	  // 1. vec-aligned bitarray of up to raw_sample_ct bits, storing which
	  //    samples have dosages.
	  // 2. word-aligned array of uint16s with 0..32768 fixed-point dosages.
	  assert(!(read_phase_dosage_gflags & kfPgenGlobalDosagePhasePresent));
	  const uintptr_t dosageraw_word_ct = kWordsPerVec * (BITCT_TO_VECCT(raw_sample_ct) + DIV_UP(raw_sample_ct, (kBytesPerVec / sizeof(dosage_t))));
	  load_variant_vec_ct += WORDCT_TO_VECCT(dosageraw_word_ct);
	}
      }
      // todo: multiallelic variants

      uintptr_t bytes_left = bigstack_left();
      if (bytes_left < 7 * kCacheline) {
	goto make_pgen_robust_ret_NOMEM;
      }
      bytes_left -= 7 * kCacheline; // defend against adverse rounding
      uintptr_t ulii = bytes_left / (2 * (kBytesPerVec * load_variant_vec_ct + loaded_vrtypes_needed));
      if (!ulii) {
	goto make_pgen_robust_ret_NOMEM;
      }
      if (ulii > MINV(kPglVblockSize, variant_ct)) {
	ulii = MINV(kPglVblockSize, variant_ct);
      }
      const uint32_t write_block_size = ulii;
      uintptr_t* main_loadbufs[2];
      main_loadbufs[0] = (uintptr_t*)bigstack_alloc_raw_rd(load_variant_vec_ct * kBytesPerVec * write_block_size);
      main_loadbufs[1] = (uintptr_t*)bigstack_alloc_raw_rd(load_variant_vec_ct * kBytesPerVec * write_block_size);
      if (loaded_vrtypes_needed) {
	g_loaded_vrtypes[0] = bigstack_alloc_raw_rd(write_block_size);
	g_loaded_vrtypes[1] = bigstack_alloc_raw_rd(write_block_size);
      } else {
	g_loaded_vrtypes[0] = nullptr;
	g_loaded_vrtypes[1] = nullptr;
      }
      
      LOGPRINTFWW5("Writing %s ... ", outname);
      fputs("0%", stdout);
      fflush(stdout);

      // Main workflow:
      // 1. Set n=0, load first write_block_size post-filtering variants
      //
      // 2. Spawn single thread processing batch n
      // 3. Load batch (n+1) unless eof
      // 4. Join thread
      // 5. Increment n by 1
      // 6. Goto step 2 unless eof
      const uint32_t batch_ct_m1 = (variant_ct - 1) / write_block_size;
      uint32_t pct = 0;
      uint32_t parity = 0;
      uint32_t read_batch_idx = 0;
      uint32_t cur_batch_size = write_block_size;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t read_variant_uidx = 0xffffffffU; // deliberate overflow
      pgr_clear_ld_cache(simple_pgrp);
      while (1) {
	if (!ts.is_last_block) {
	  if (read_batch_idx == batch_ct_m1) {
	    cur_batch_size = variant_ct - (read_batch_idx * write_block_size);
	  }
	  uintptr_t* cur_loadbuf = main_loadbufs[parity];
	  uintptr_t* loadbuf_iter = cur_loadbuf;
	  unsigned char* cur_loaded_vrtypes = g_loaded_vrtypes[parity];
	  g_loadbuf_thread_starts[parity][0] = loadbuf_iter;
	  for (uint32_t uii = 0; uii < cur_batch_size; ++uii) {
	    if (!new_variant_idx_to_old_iter) {
	      ++read_variant_uidx;
	      next_set_unsafe_ck(variant_include, &read_variant_uidx);
	    } else {
	      read_variant_uidx = *new_variant_idx_to_old_iter++;
	    }
	    reterr = pgr_read_raw(read_variant_uidx, read_phase_dosage_gflags, simple_pgrp, &loadbuf_iter, cur_loaded_vrtypes? (&(cur_loaded_vrtypes[uii])) : nullptr);
	    if (reterr) {
	      if (reterr == kPglRetMalformedInput) {
		logprint("\n");
		logerrprint("Error: Malformed .pgen file.\n");
	      }
	      goto make_pgen_robust_ret_1;
	    }
	  }
	}
	if (read_batch_idx) {
	  join_threads3z(&ts);
	  reterr = g_error_ret;
	  if (reterr) {
	    goto make_pgen_robust_ret_WRITE_FAIL;
	  }
	}
	if (!ts.is_last_block) {
	  g_cur_block_write_ct = cur_batch_size;
	  ts.is_last_block = (read_batch_idx == batch_ct_m1);
	  ts.thread_func_ptr = make_pgen_thread;
	  if (spawn_threads3z(read_batch_idx, &ts)) {
	    goto make_pgen_robust_ret_THREAD_CREATE_FAIL;
	  }
	}
	parity = 1 - parity;
	if (read_batch_idx) {
	  if (read_batch_idx > batch_ct_m1) {
	    break;
	  }
	  const uint32_t write_idx_end = read_batch_idx * write_block_size;
	  if (write_idx_end >= next_print_variant_idx) {
	    if (pct > 10) {
	      putc_unlocked('\b', stdout);
	    }
	    pct = (write_idx_end * 100LLU) / variant_ct;
	    printf("\b\b%u%%", pct++);
	    fflush(stdout);
	    next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
	  }
	}
	++read_batch_idx;
      }
      spgw_finish(g_spgwp);
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      LOGPRINTF("done.\n");
    }
  }
  while (0) {
  make_pgen_robust_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  make_pgen_robust_ret_WRITE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  make_pgen_robust_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 make_pgen_robust_ret_1:
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t make_plink2_no_vsort(const char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, pvar_psam_t pvar_psam_modifier, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  threads_state_t ts;
  init_threads3z(&ts);
  mt_pgen_writer_t* mpgwp = nullptr;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (make_plink2_modifier & kfMakePlink2MMask) {
      logerrprint("Error: --make-bed/--make-{b}pgen multiallelics= is currently under development.\n");
      reterr = kPglRetNotYetSupported;
      goto make_plink2_no_vsort_ret_1;
    }    
    g_plink2_write_flags = kfPlink2Write0;
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    if (make_plink2_modifier & kfMakePlink2SetHhMissing) {
      const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
      uintptr_t* sex_collapsed_tmp;
      uintptr_t* sex_female;
      if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
	  bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed_interleaved) ||
	  bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_collapsed_tmp) ||
	  bigstack_alloc_ul(raw_sample_ctl, &sex_female)) {
	goto make_plink2_no_vsort_ret_NOMEM;
      }
      copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_collapsed_tmp);
      fill_interleaved_mask_vec(sex_collapsed_tmp, sample_ctv, g_sex_male_collapsed_interleaved);

      bitvec_andnot_copy(sex_nm, sex_male, raw_sample_ctl, sex_female);
      copy_bitarr_subset(sex_female, sample_include, sample_ct, sex_collapsed_tmp);
      fill_interleaved_mask_vec(sex_collapsed_tmp, sample_ctv, g_sex_female_collapsed_interleaved);
      
      bigstack_reset(sex_collapsed_tmp);
      g_plink2_write_flags |= kfPlink2WriteSetHhMissing;
    }
    if (make_plink2_modifier & kfMakePlink2SetMixedMtMissing) {
      g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
    }
    g_cip = cip;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    const uint32_t make_pgen = make_plink2_modifier & kfMakePgen;
    // todo: prohibit .pgen + .bim write when data is multiallelic without
    //   either multiallelic split or erase-alt2+ specified
    //   (--make-bed = automatic erase-alt2+)
    if ((make_plink2_modifier & kfMakeBed) || ((make_plink2_modifier & (kfMakePgen | (kfMakePgenFormatBase * 3))) == (kfMakePgen | kfMakePgenFormatBase))) {
      // fixed-width
      if (make_pgen) {
        strcpy(outname_end, ".pgen");
      } else {
        strcpy(outname_end, ".bed");
      }
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto make_plink2_no_vsort_ret_OPEN_FAIL;
      }
      if (make_pgen) {
	fwrite("l\x1b\x02", 3, 1, outfile);
        fwrite(&variant_ct, 4, 1, outfile);
	fwrite(&sample_ct, 4, 1, outfile);
	if (!pgfip->nonref_flags) {
	  const pgen_global_flags_t gflags = pgfip->gflags;
	  uint32_t uii = 64;
	  if (gflags & kfPgenGlobalAllNonref) {
	    uii = 128;
	  }
	  putc_unlocked(uii, outfile);
	} else {
	  putc_unlocked(192, outfile);
	  fwrite(pgfip->nonref_flags, DIV_UP(variant_ct, CHAR_BIT), 1, outfile);
	}
	if (ferror(outfile)) {
	  goto make_plink2_no_vsort_ret_WRITE_FAIL;
	}
      } else {
	if (fwrite_checked("l\x1b\x01", 3, outfile)) {
	  goto make_plink2_no_vsort_ret_WRITE_FAIL;
	}
      }
      LOGPRINTFWW5("Writing %s ... ", outname);
      fputs("0%", stdout);
      fflush(stdout);
      uint32_t pct = 0;
      if (variant_ct && sample_ct) {
	const uintptr_t sample_ct4 = QUATERCT_TO_BYTECT(sample_ct);
	if (bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_cumulative_popcounts) ||
	    bigstack_alloc_uc(sample_ct4 * kPglVblockSize, &(g_writebufs[0])) ||
	    bigstack_alloc_uc(sample_ct4 * kPglVblockSize, &(g_writebufs[1]))) {
	  // todo: low-memory single-threaded fallback mode
	  goto make_plink2_no_vsort_ret_NOMEM;
	}
	fill_cumulative_popcounts(sample_include, raw_sample_ctl, g_sample_include_cumulative_popcounts);
	// tried more threads, pointless since this is too I/O-bound
	// (exception: reordering samples)
	uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
	g_collapsed_sort_map = new_sample_idx_to_old;
	if (!new_sample_idx_to_old) {
	  // subsetting is most expensive with sample_ct near 2/3 of
	  // raw_sample_ct; up to ~7 compute threads are useful in that case.
	  // (see copy_quaterarr_nonempty_subset().)
	  uint64_t numer;
	  if (sample_ct * (3 * k1LU) <= raw_sample_ct * (2 * k1LU)) {
	    numer = sample_ct * (9 * k1LU);
	  } else {
	    numer = (raw_sample_ct - sample_ct) * (18 * k1LU);
	  }
	  const uint32_t calc_thread_max = 1 + (numer / raw_sample_ct);
	  if (calc_thread_max < calc_thread_ct) {
	    calc_thread_ct = calc_thread_max;
	  }
	} else if (sample_ct < raw_sample_ct) {
	  // const_cast
	  if (bigstack_alloc_ui(sample_ct, (uint32_t**)((uintptr_t)(&g_collapsed_sort_map)))) {
	    goto make_plink2_no_vsort_ret_NOMEM;
	  }
	  uidxs_to_idxs(sample_include, g_sample_include_cumulative_popcounts, sample_ct, (uint32_t*)((uintptr_t)g_collapsed_sort_map));
	}

	if (make_plink2_modifier & kfMakeBed) {
	  g_plink2_write_flags |= kfPlink2WritePlink1;
	}
	
	unsigned char* main_loadbufs[2];
	uint32_t read_block_size;
	if (multithread_load_init(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, nullptr, nullptr, &read_block_size, main_loadbufs, &ts.threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
	  goto make_plink2_no_vsort_ret_NOMEM;
	}

	g_variant_include = variant_include;
	g_refalt1_select = refalt1_select;
	g_sample_include = sample_include;
	g_sample_ct = sample_ct;
	ts.calc_thread_ct = calc_thread_ct;
	g_calc_thread_ct = calc_thread_ct;
	g_error_ret = kPglRetSuccess;

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
	// 8. Write results for last block

	const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
	const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
	uint32_t parity = 0;
	uint32_t read_block_idx = 0;
	uint32_t prev_variant_idx = 0;
	uint32_t variant_idx = 0;
	uint32_t cur_read_block_size = read_block_size;
	uint32_t next_print_variant_idx = variant_ct / 100;
	while (1) {
	  uintptr_t cur_block_write_ct = 0;
	  if (!ts.is_last_block) {
	    while (read_block_idx < read_block_ct_m1) {
	      cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	      if (cur_block_write_ct) {
		break;
	      }
	      ++read_block_idx;
	    }
	    if (read_block_idx == read_block_ct_m1) {
	      cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	      cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	    }
	    if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
	      goto make_plink2_no_vsort_ret_READ_FAIL;
	    }
	  }
	  if (variant_idx) {
	    join_threads3z(&ts);
	    if (reterr) {
	      if (reterr == kPglRetMalformedInput) {
		logprint("\n");
		logerrprint("Error: Malformed .pgen file.\n");
	      }
	      goto make_plink2_no_vsort_ret_1;
	    }
	  }
	  if (!ts.is_last_block) {
	    g_cur_block_write_ct = cur_block_write_ct;
	    compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
	    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	      g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	      g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	    }
	    ts.is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
	    ts.thread_func_ptr = make_bedlike_thread;
	    if (spawn_threads3z(variant_idx, &ts)) {
	      goto make_plink2_no_vsort_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  parity = 1 - parity;
	  if (variant_idx) {
	    // write *previous* block results
	    if (fwrite_checked(g_writebufs[parity], (variant_idx - prev_variant_idx) * sample_ct4, outfile)) {
	      goto make_plink2_no_vsort_ret_WRITE_FAIL;
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
	    prev_variant_idx = variant_idx;
	  }
	  ++read_block_idx;
	  variant_idx += cur_block_write_ct;
	  // crucially, this is independent of the pgen_reader_t block_base
	  // pointers
	  pgfip->block_base = main_loadbufs[parity];
	}
      }
      if (fclose_null(&outfile)) {
	goto make_plink2_no_vsort_ret_WRITE_FAIL;
      }
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      LOGPRINTF("done.\n");
      bigstack_reset(bigstack_mark);
    } else if (make_pgen) {
      // should be straightforward to make this sort variants...
      if ((!variant_ct) || (!sample_ct)) {
	logerrprint("Error: Zero-variant/zero-sample .pgen writing is not currently supported.\n");
	reterr = kPglRetNotYetSupported;
	goto make_plink2_no_vsort_ret_1;
      }
      const uint32_t input_biallelic = (!variant_allele_idxs);
      const uintptr_t* write_allele_idx_offsets = nullptr;
      if (!input_biallelic) {
        if (variant_ct < raw_variant_ct) {
	  uintptr_t* new_allele_idx_offsets;
	  if (bigstack_alloc_ul(variant_ct + 1, &new_allele_idx_offsets)) {
	    goto make_plink2_no_vsort_fallback;
	  }
	  uintptr_t cur_offset = 0;
	  uint32_t variant_uidx = 0;
	  // todo: separate trim-alts case
	  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	    next_set_unsafe_ck(variant_include, &variant_uidx);
	    new_allele_idx_offsets[variant_idx] = cur_offset;
	    cur_offset += variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx];
	  }
	  if (cur_offset != 2 * variant_ct) {
	    new_allele_idx_offsets[variant_ct] = cur_offset;
	    write_allele_idx_offsets = new_allele_idx_offsets;
	    logprint("Error: Multiallelic .pgen write is not yet supported.\n");
	    reterr = kPglRetNotYetSupported;
	    goto make_plink2_no_vsort_ret_1;
	  } else {
	    bigstack_reset(new_allele_idx_offsets);
	  }
	} else {
	  write_allele_idx_offsets = variant_allele_idxs;
	}
      }
      if (variant_ct == raw_variant_ct) {
	g_write_chr_fo_vidx_start = cip->chr_fo_vidx_start;
      } else {
	if (alloc_and_fill_subset_chr_fo_vidx_start(variant_include, cip, &g_write_chr_fo_vidx_start)) {
	  goto make_plink2_no_vsort_fallback;
	}
      }
      pgen_global_flags_t read_phase_dosage_gflags = pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      if (make_plink2_modifier & kfMakePgenErasePhase) {
	read_phase_dosage_gflags &= ~(kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
      }
      if (make_plink2_modifier & kfMakePgenEraseDosage) {
	if (hard_call_thresh == 0xffffffffU) {
	  read_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
	} else {
	  // erase-dosage + --hard-call-threshold currently requires dosages to
	  // be read, and only thrown away at the last minute
	  // (alternatively, we could build --hard-call-threshold directly into
	  // pgr_read_raw?)
	  g_plink2_write_flags |= kfPlink2WriteLateDosageErase;
	}
      }
      if (read_phase_dosage_gflags && (variant_ct < raw_variant_ct)) {
	// did we e.g. filter out all the phased variants?
	read_phase_dosage_gflags = gflags_vfilter(variant_include, pgfip->vrtypes, raw_variant_ct, pgfip->gflags);
      }
      // could check if all the phased samples were also filtered out, but
      // that's already caught by running --make-pgen twice, so not a big deal
      g_read_phase_dosage_gflags = read_phase_dosage_gflags;
      g_hard_call_halfdist = (hard_call_thresh == 0xffffffffU)? 0 : (kDosage4th - hard_call_thresh);
      g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_hphase_present = (read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
      const uint32_t read_dosage_present = (read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
      pgen_global_flags_t write_phase_dosage_gflags = read_phase_dosage_gflags;
      if (g_plink2_write_flags & kfPlink2WriteLateDosageErase) {
	write_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      }
      uintptr_t alloc_base_cacheline_ct;
      uint64_t mpgw_per_thread_cacheline_ct;
      uint32_t vrec_len_byte_ct;
      uint64_t vblock_cacheline_ct;
      // may want to have a load_sample_ct which is raw_sample_ct when e.g.
      // sample_ct > 0.1 * raw_sample_ct, and sample_ct otherwise.
      mpgw_init_phase1(write_allele_idx_offsets, variant_ct, sample_ct, write_phase_dosage_gflags, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);

      // bugfix: each variant currently needs to be vector-aligned
      // bugfix?: need to use raw_sample_ct here, not sample_ct
      const uint32_t raw_sample_ctv2 = QUATERCT_TO_VECCT(raw_sample_ct);
      const uint32_t max_vblock_size = MINV(kPglVblockSize, variant_ct);
      uint64_t load_vblock_cacheline_ct = VECCT_TO_CLCT(((uint64_t)raw_sample_ctv2) * max_vblock_size);

      if (read_hphase_present) {
	// could make this bound tighter when lots of unphased variants are
	// mixed in among the phased variants, but this isn't nearly as
	// important as the analogous multiallelic optimization

	// phaseraw has two parts:
	// 1. vec-aligned bitarray of up to (raw_sample_ct + 1) bits.  first
	//    bit is set iff phasepresent is explicitly stored at all (if not,
	//    all hets are assumed to be phased), if yes the remaining bits
	//    store packed phasepresent values for all hets, if no the
	//    remaining bits store packed phaseinfo values for all hets.
	// 2. word-aligned bitarray of up to raw_sample_ct bits, storing
	//    phaseinfo values.  (end of this array is vec-aligned.)
	const uintptr_t phaseraw_word_ct = kWordsPerVec + round_down_pow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
	load_vblock_cacheline_ct += WORDCT_TO_CLCT(((uint64_t)phaseraw_word_ct) * max_vblock_size);
      }
      if (read_dosage_present) {	
	// todo: phased dosage

	// (unphased, biallelic) dosageraw has two parts:
	// 1. vec-aligned bitarray of up to raw_sample_ct bits, storing which
	//    samples have dosages.
	// 2. word-aligned array of uint16s with 0..32768 fixed-point dosages.
	assert(!(read_phase_dosage_gflags & kfPgenGlobalDosagePhasePresent));
	const uintptr_t dosageraw_word_ct = kWordsPerVec * (BITCT_TO_VECCT(raw_sample_ct) + DIV_UP(raw_sample_ct, (kBytesPerVec / sizeof(dosage_t))));
	load_vblock_cacheline_ct += WORDCT_TO_CLCT(dosageraw_word_ct * ((uint64_t)max_vblock_size));
      }
      // todo: multiallelic variants
      
#ifndef __LP64__
      if ((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (load_vblock_cacheline_ct > (0x7fffffff / kCacheline))) {
	goto make_plink2_no_vsort_fallback;
      }
#endif
      uint32_t calc_thread_ct = DIV_UP(variant_ct, kPglVblockSize);
      if (calc_thread_ct >= max_thread_ct) {
	calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      }
      const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
      if (!new_sample_idx_to_old) {
	// hphase doesn't seem to affect read:write ratio much
	const uint32_t max_calc_thread_ct = 2 + subsetting_required;
	if (calc_thread_ct > max_calc_thread_ct) {
	  calc_thread_ct = max_calc_thread_ct;
	}
      }
      // this is frequently I/O-bound even when resorting, but I'll postpone
      // tuning thread count there
      g_refalt1_select = refalt1_select;
      if (refalt1_select && (variant_ct < raw_variant_ct)) {
	// might want inner loop to map variant uidx -> idx instead
        g_refalt1_select = (alt_allele_ct_t*)bigstack_alloc(variant_ct * 2 * sizeof(alt_allele_ct_t));
	if (!g_refalt1_select) {
	  goto make_plink2_no_vsort_fallback;
	}
	uint32_t variant_uidx = 0;
	for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	  next_set_unsafe_ck(variant_include, &variant_uidx);
	  // const_cast
	  memcpy((alt_allele_ct_t*)((uintptr_t)(&(g_refalt1_select[2 * variant_idx]))), &(refalt1_select[2 * variant_uidx]), 2 * sizeof(alt_allele_ct_t));
	}
      }
      mpgwp = (mt_pgen_writer_t*)bigstack_alloc((calc_thread_ct + DIV_UP(sizeof(mt_pgen_writer_t), kBytesPerWord)) * sizeof(intptr_t));
      if (!mpgwp) {
	goto make_plink2_no_vsort_fallback;
      }
      mpgwp->pgen_outfile = nullptr;
      if (bigstack_alloc_thread(calc_thread_ct, &ts.threads) ||
	  bigstack_alloc_ulp(calc_thread_ct, &(g_loadbuf_thread_starts[0])) ||
	  bigstack_alloc_ulp(calc_thread_ct, &(g_loadbuf_thread_starts[1]))) {
	goto make_plink2_no_vsort_fallback;
      }
      uint32_t nonref_flags_storage = 3;
      if (!pgfip->nonref_flags) {
	nonref_flags_storage = (simple_pgrp->fi.gflags & kfPgenGlobalAllNonref)? 2 : 1;
      } else if (variant_ct < raw_variant_ct) {
	// todo: check if now constant
      }
      g_pwcs = &(mpgwp->pwcs[0]);
      g_new_sample_idx_to_old = new_sample_idx_to_old;
      g_thread_write_genovecs = nullptr;
      uintptr_t other_per_thread_cacheline_ct = 2 * load_vblock_cacheline_ct;
      if (new_sample_idx_to_old || subsetting_required) {
	if (bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_genovecs)) {
	  goto make_plink2_no_vsort_fallback;
	}
	if (read_hphase_present && new_sample_idx_to_old) {
	  if (bigstack_alloc_ui(raw_sample_ct, &g_old_sample_idx_to_new)) {
	    goto make_plink2_no_vsort_fallback;
	  }
	  for (uint32_t new_sample_idx = 0; new_sample_idx < sample_ct; ++new_sample_idx) {
	    g_old_sample_idx_to_new[new_sample_idx_to_old[new_sample_idx]] = new_sample_idx;
	  }
	}
	if (read_dosage_present && new_sample_idx_to_old) {
	  // g_thread_cumulative_popcount_bufs
	  other_per_thread_cacheline_ct += INT32CT_TO_CLCT(raw_sample_ctl);
	}
	// per-thread output buffers required
	if (input_biallelic) {
	  other_per_thread_cacheline_ct += QUATERCT_TO_CLCT(sample_ct);
	} else {
	  other_per_thread_cacheline_ct += DIV_UP(2 * sample_ct * sizeof(alt_allele_ct_t), kCacheline);
	}
      }
      if (read_hphase_present || read_dosage_present) {
	if (read_hphase_present) {
	  if (bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_phasepresents) ||
	      bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_phaseinfos) ||
	      bigstack_alloc_ulp(calc_thread_ct, &g_thread_all_hets)) {
	    goto make_plink2_no_vsort_fallback;
	  }
	  // phasepresent, phaseinfo
	  other_per_thread_cacheline_ct += 2 * BITCT_TO_CLCT(sample_ct);

	  // all_hets
	  other_per_thread_cacheline_ct += BITCT_TO_CLCT(raw_sample_ct);
	}
	if (read_dosage_present) {
	  if (bigstack_alloc_dosagep(calc_thread_ct, &g_thread_write_dosagevals) ||
	      bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_dosagepresents)) {
	    goto make_plink2_no_vsort_fallback;
	  }
	  if (new_sample_idx_to_old) {
	    if (bigstack_alloc_uip(calc_thread_ct, &g_thread_cumulative_popcount_bufs)) {
	      goto make_plink2_no_vsort_fallback;
	    }
	  }
	  // dosage_present
	  other_per_thread_cacheline_ct += BITCT_TO_CLCT(sample_ct);

	  // dosage_vals
	  other_per_thread_cacheline_ct += DIV_UP(sample_ct, (kCacheline / sizeof(dosage_t)));
	}
	// g_loaded_vrtypes
	other_per_thread_cacheline_ct += 2 * (kPglVblockSize / kCacheline);
      }
      const uintptr_t cachelines_avail = bigstack_left() / kCacheline;
      if (cachelines_avail < alloc_base_cacheline_ct + (mpgw_per_thread_cacheline_ct + other_per_thread_cacheline_ct) * calc_thread_ct) {
	if (cachelines_avail < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct + other_per_thread_cacheline_ct) {
	  goto make_plink2_no_vsort_fallback;
	}
	calc_thread_ct = (cachelines_avail - alloc_base_cacheline_ct) / (mpgw_per_thread_cacheline_ct + other_per_thread_cacheline_ct);
      }
      uintptr_t* main_loadbufs[2];
      main_loadbufs[0] = (uintptr_t*)bigstack_alloc_raw(load_vblock_cacheline_ct * calc_thread_ct * kCacheline);
      main_loadbufs[1] = (uintptr_t*)bigstack_alloc_raw(load_vblock_cacheline_ct * calc_thread_ct * kCacheline);
      if (read_hphase_present || read_dosage_present) {
	g_loaded_vrtypes[0] = bigstack_alloc_raw(kPglVblockSize * calc_thread_ct);
	g_loaded_vrtypes[1] = bigstack_alloc_raw(kPglVblockSize * calc_thread_ct);
	const uint32_t bitvec_writebuf_byte_ct = BITCT_TO_CLCT(sample_ct) * kCacheline;
	const uintptr_t dosagevals_writebuf_byte_ct = DIV_UP(sample_ct, (kCacheline / 2)) * kCacheline;
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  if (read_hphase_present) {
	    g_thread_write_phasepresents[tidx] = (uintptr_t*)bigstack_alloc_raw(bitvec_writebuf_byte_ct);
	    g_thread_write_phaseinfos[tidx] = (uintptr_t*)bigstack_alloc_raw(bitvec_writebuf_byte_ct);

	    g_thread_all_hets[tidx] = (uintptr_t*)bigstack_alloc_raw(BITCT_TO_CLCT(raw_sample_ct) * kCacheline);
	  }
	  if (read_dosage_present) {
	    g_thread_write_dosagepresents[tidx] = (uintptr_t*)bigstack_alloc_raw(bitvec_writebuf_byte_ct);
	    g_thread_write_dosagevals[tidx] = (dosage_t*)bigstack_alloc_raw(dosagevals_writebuf_byte_ct);
	    if (new_sample_idx_to_old) {
	      g_thread_cumulative_popcount_bufs[tidx] = (uint32_t*)bigstack_alloc_raw(INT32CT_TO_CLCT(raw_sample_ctl) * kCacheline);
	    }
	  }
	}
      } else {
	g_loaded_vrtypes[0] = nullptr;
	g_loaded_vrtypes[1] = nullptr;
      }
      if (new_sample_idx_to_old || subsetting_required) {
	uintptr_t writebuf_byte_ct = input_biallelic? QUATERCT_TO_BYTECT(sample_ct) : (2 * sample_ct * sizeof(alt_allele_ct_t));
	writebuf_byte_ct = round_up_pow2(writebuf_byte_ct, kCacheline);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_thread_write_genovecs[tidx] = (uintptr_t*)bigstack_alloc_raw(writebuf_byte_ct);
	}
      }
      strcpy(outname_end, ".pgen");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fputs("0%", stdout);
      fflush(stdout);
      unsigned char* mpgw_alloc = bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline);
      assert(g_bigstack_base <= g_bigstack_end);
      reterr = mpgw_init_phase2(outname, write_allele_idx_offsets, pgfip->nonref_flags, variant_ct, sample_ct, write_phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);
      if (reterr) {
	goto make_plink2_no_vsort_ret_1;
      }
      g_sample_include = subsetting_required? sample_include : nullptr;
      g_raw_sample_ct = raw_sample_ct;
      g_sample_ct = sample_ct;
      ts.calc_thread_ct = calc_thread_ct;
      // g_calc_thread_ct = calc_thread_ct;
      g_spgwp = nullptr;
      // g_error_ret = kPglRetSuccess;

      // Main workflow:
      // 1. Set n=0, load first calc_thread_ct * kPglVblockSize
      //    *post-filtering* variants.
      //    This doesn't play well with blockload when any variants are
      //    filtered out, so we don't use it.  (todo: look into special-casing
      //    variant_ct == raw_variant_ct.)
      //
      // 2. Spawn threads processing batch n
      // 3. Load batch (n+1) unless eof
      // 4. Join threads
      // 5. Flush results for batch n (must happen here since we aren't using
      //    two output buffers.  this may be a mistake, revisit this choice...)
      // 6. Increment n by 1
      // 7. Goto step 2 unless eof
      const uint32_t batch_ct_m1 = (variant_ct - 1) / (kPglVblockSize * calc_thread_ct);
      uint32_t pct = 0;
      uint32_t parity = 0;
      uint32_t read_batch_idx = 0;
      uint32_t write_idx_end = 0;
      uint32_t cur_batch_size = kPglVblockSize * calc_thread_ct;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t read_variant_uidx = next_set_unsafe(variant_include, 0);
      pgr_clear_ld_cache(simple_pgrp);
      while (1) {
	if (read_batch_idx) {
	  g_cur_block_write_ct = cur_batch_size;
	  ts.is_last_block = (write_idx_end == variant_ct);
	  ts.thread_func_ptr = make_pgen_thread;
	  if (spawn_threads3z(read_batch_idx - 1, &ts)) {
	    goto make_plink2_no_vsort_ret_THREAD_CREATE_FAIL;
	  }
	}
	if (!ts.is_last_block) {
	  if (read_batch_idx == batch_ct_m1) {
	    cur_batch_size = variant_ct - (read_batch_idx * kPglVblockSize * calc_thread_ct);
	  }
	  uintptr_t* cur_loadbuf = main_loadbufs[parity];
	  uintptr_t* loadbuf_iter = cur_loadbuf;
	  unsigned char* cur_loaded_vrtypes = g_loaded_vrtypes[parity];
	  for (uint32_t uii = 0; uii < cur_batch_size; ++uii, ++read_variant_uidx) {
	    if (!(uii % kPglVblockSize)) {
	      g_loadbuf_thread_starts[parity][uii / kPglVblockSize] = loadbuf_iter;
	    }
	    next_set_unsafe_ck(variant_include, &read_variant_uidx);
	    reterr = pgr_read_raw(read_variant_uidx, read_phase_dosage_gflags, simple_pgrp, &loadbuf_iter, cur_loaded_vrtypes? (&(cur_loaded_vrtypes[uii])) : nullptr);
	    if (reterr) {
	      if (reterr == kPglRetMalformedInput) {
		logprint("\n");
		logerrprint("Error: Malformed .pgen file.\n");
	      }
	      goto make_plink2_no_vsort_ret_1;
	    }
	  }
	}
	if (read_batch_idx) {
	  join_threads3z(&ts);
	}
	parity = 1 - parity;
	if (write_idx_end) {
	  reterr = mpgw_flush(mpgwp);
	  if (reterr) {
	    goto make_plink2_no_vsort_ret_WRITE_FAIL;
	  }
	  if (write_idx_end == variant_ct) {
	    mpgwp = nullptr;
	    break;
	  }
	  if (write_idx_end >= next_print_variant_idx) {
	    if (pct > 10) {
	      putc_unlocked('\b', stdout);
	    }
	    pct = (write_idx_end * 100LLU) / variant_ct;
	    printf("\b\b%u%%", pct++);
	    fflush(stdout);
	    next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
	  }
	}
	++read_batch_idx;
	write_idx_end += cur_batch_size;
      }
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      LOGPRINTF("done.\n");
      bigstack_reset(bigstack_mark);
    } else if (0) {
    make_plink2_no_vsort_fallback:
      mpgwp = nullptr;
      bigstack_reset(bigstack_mark2);
      reterr = make_pgen_robust(sample_include, new_sample_idx_to_old, variant_include, variant_allele_idxs, refalt1_select, nullptr, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, hard_call_thresh, dosage_erase_thresh, make_plink2_modifier, simple_pgrp, outname, outname_end);
      if (reterr) {
	goto make_plink2_no_vsort_ret_1;
      }
    }
    const uint32_t trim_alts = (uint32_t)(make_plink2_modifier & kfMakePlink2TrimAlts);
    // don't bother with trim-alts/set-hh-missing interaction for now
    if (make_plink2_modifier & kfMakeBim) {
      char* bimname_end = strcpya0(outname_end, ".bim");
      const uint32_t bim_zst = (make_plink2_modifier / kfMakeBimZs) & 1;
      if (bim_zst) {
	strcpy(bimname_end, ".zst");
      }
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_map_or_bim(outname, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, '\t', bim_zst);
      if (reterr) {
	goto make_plink2_no_vsort_ret_1;
      }
      logprint("done.\n");
    }
    if (make_plink2_modifier & kfMakePvar) {
      char* pvarname_end = strcpya0(outname_end, ".pvar");
      if (pvar_psam_modifier & kfPvarZs) {
	strcpy(pvarname_end, ".zst");
      }
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t nonref_flags_storage = 3;
      if (!pgfip->nonref_flags) {
	nonref_flags_storage = (pgfip->gflags & kfPgenGlobalAllNonref)? 2 : 1;
      }
      reterr = write_pvar(outname, xheader, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, xheader_info_pr, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_modifier);
      if (reterr) {
	goto make_plink2_no_vsort_ret_1;
      }
      logprint("done.\n");
    }
    if (make_plink2_modifier & kfMakeFam) {
      strcpy(outname_end, ".fam");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_fam(outname, sample_include, sample_ids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, new_sample_idx_to_old, sample_ct, max_sample_id_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, '\t');
      if (reterr) {
	goto make_plink2_no_vsort_ret_1;
      }
      logprint("done.\n");
    }
    if (make_plink2_modifier & kfMakePsam) {
      strcpy(outname_end, ".psam");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_psam(outname, sample_include, sample_ids, sids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, sample_ct, max_sample_id_blen, max_sid_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, max_pheno_name_blen, pvar_psam_modifier);
      if (reterr) {
	goto make_plink2_no_vsort_ret_1;
      }
      logprint("done.\n");
    }
  }
  while (0) {
  make_plink2_no_vsort_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  make_plink2_no_vsort_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  make_plink2_no_vsort_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  make_plink2_no_vsort_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  make_plink2_no_vsort_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 make_plink2_no_vsort_ret_1:
  if (mpgw_cleanup(mpgwp) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  bigstack_reset(bigstack_mark);
  return reterr;
}

/*
pglerr_t write_pvar_resorted(const char* outname, const char* xheader, const uintptr_t* variant_include, const chr_info_t* write_cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, char** filter_storage, const uintptr_t* nonref_flags, char** pvar_info_strs, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, pvar_psam_t pvar_psam_modifier) {
  // allele_dosages must be nullptr unless we're trimming alt alleles

  // The annoying part of this is handling a sequence of INFO strings that
  // don't fit in memory; use a multipass approach for that.
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t max_chr_blen = get_max_chr_slen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    unsigned char* overflow_buf;
    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
	bigstack_alloc_uc(overflow_buf_size, &overflow_buf) ||
	bigstack_alloc_ul(BITCT_TO_WORDCT(kPglMaxAltAlleleCt), &allele_include)) {
      goto write_pvar_resorted_ret_NOMEM;
    }
    const uint32_t output_zst = (pvar_psam_modifier / kfPvarZs) & 1;
    if (cswrite_init(outname, 0, output_zst, overflow_buf, &css)) {
      goto write_pvar_resorted_ret_OPEN_FAIL;
    }
    cswritep = (char*)overflow_buf;
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_modifier & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      uint32_t widx;
      for (widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & nonref_flags[widx]) {
	  break;
	}
      }
      if (widx == raw_variant_ctl) {
	write_info_pr = 0;
      }
    }
    write_info_pr = write_info_pr && write_info;

    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t info_col_idx = 0; // could save this during first load instead
    if (pvar_psam_modifier & kfPvarColXheader) {
      if (csputs_std(xheader, xheader_blen, &css, &cswritep)) {
	goto write_pvar_resorted_ret_WRITE_FAIL;
      }
      if (write_info_pr && (!xheader_info_pr)) {
	cswritep = strcpya(cswritep, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
      }
    }
    if (write_cip->chrset_source) {
      append_chrset_line(write_cip, &cswritep);
    }
    cswritep = strcpya(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_modifier & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_modifier & kfPvarColMaybequal) && qual_present) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & qual_present[widx]) {
	  write_qual = 1;
	  break;
	}
      }
    }
    if (write_qual) {
      cswritep = strcpya(cswritep, "\tQUAL");
    }
    
    uint32_t write_filter = 0;
    if (pvar_psam_modifier & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_modifier & kfPvarColMaybefilter) && filter_present) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & filter_present[widx]) {
	  write_filter = 1;
	  break;
	}
      }
    }
    if (write_filter) {
      cswritep = strcpya(cswritep, "\tFILTER");
    }

    if (write_info) {
      cswritep = strcpya(cswritep, "\tINFO");
    }
    
    uint32_t write_cm = 0;
    if (pvar_psam_modifier & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_modifier & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
	// nonzero_cm_present check was performed
	write_cm = 1;
      } else {
	uint32_t variant_uidx = 0;
	for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	  next_set_unsafe_ck(variant_include, &variant_uidx);
	  if (variant_cms[variant_uidx] != 0.0) {
	    write_cm = 1;
	    break;
	  }
	}
      }
    }
    if (write_cm) {
      cswritep = memcpyl3a(cswritep, "\tCM");
    }
    append_binary_eoln(&cswritep);

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    const uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = *new_variant_idx_to_old_iter++;
      if (variant_idx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = write_cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_idx >= chr_end);
	char* chr_name_end = chr_name_write(write_cip, write_cip->chr_file_order[chr_fo_idx], chr_buf);
	*chr_name_end = '\t';
	chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
      }
      cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
      cswritep = uint32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
      uintptr_t variant_allele_idx_base;
      if (!variant_allele_idxs) {
	variant_allele_idx_base = variant_uidx * 2;
      } else {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
	alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
      }
      cswritep = strcpyax(cswritep, cur_alleles[ref_allele_idx], '\t');
      if ((!allele_dosages) || allele_dosages[variant_allele_idx_base + alt1_allele_idx]) {
        cswritep = strcpya(cswritep, cur_alleles[alt1_allele_idx]);
      } else {
	*cswritep++ = output_missing_geno_char;
      }
      if (cswrite(&css, &cswritep)) {
	goto write_pvar_resorted_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
	fill_all_bits(cur_allele_ct, allele_include);
	CLEAR_BIT(ref_allele_idx, allele_include);
	CLEAR_BIT(alt1_allele_idx, allele_include);
        uint32_t cur_allele_uidx = 0;
	uint32_t alt_allele_idx = 2;
	do {
	  *cswritep++ = ',';
	  next_set_unsafe_ck(allele_include, &cur_allele_uidx);
	  cswritep = strcpya(cswritep, cur_alleles[cur_allele_uidx++]);
	  if (cswrite(&css, &cswritep)) {
	    goto write_pvar_resorted_ret_WRITE_FAIL;
	  }
	} while (++alt_allele_idx < cur_allele_ct);
      }

      if (write_qual) {
	*cswritep++ = '\t';
	if (!IS_SET(qual_present, variant_uidx)) {
	  *cswritep++ = '.';
	} else {
	  cswritep = ftoa_g(quals[variant_uidx], cswritep);
	}
      }

      if (write_filter) {
	*cswritep++ = '\t';
	if (!IS_SET(filter_present, variant_uidx)) {
	  *cswritep++ = '.';
	} else if (!IS_SET(filter_npass, variant_uidx)) {
	  cswritep = strcpya(cswritep, "PASS");
	} else {
	  cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
	}
      }

      if (write_info) {
	*cswritep++ = '\t';
	const uint32_t is_pr = all_nonref || (nonref_flags && IS_SET(nonref_flags, variant_uidx));
	if (pvar_info_strs && pvar_info_strs[variant_uidx]) {
	  pvar_info_write(pvar_info_strs[variant_uidx], xheader_info_pr, is_pr, cswritep);
	} else {
	  if (is_pr) {
	    cswritep = strcpya(cswritep, "PR");
	  } else {
	    *cswritep++ = '.';
	  }
	}
      }
      
      if (write_cm) {
        *cswritep++ = '\t';
	if (!variant_cms) {
	  *cswritep++ = '0';
	} else {
	  cswritep = dtoa_g_p8(variant_cms[variant_uidx], cswritep);
	}
      }
      append_binary_eoln(&cswritep);
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto write_pvar_resorted_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_pvar_resorted_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_pvar_resorted_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_pvar_resorted_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 write_pvar_resorted_ret_1:
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}
*/

/*
pglerr_t make_plink2_vsort(const char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, pvar_psam_t pvar_psam_modifier, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (make_plink2_modifier & kfMakeBim) {
      char* bimname_end = strcpya0(outname_end, ".bim");
      const uint32_t bim_zst = (make_plink2_modifier / kfMakeBimZs) & 1;
      if (bim_zst) {
	strcpy(bimname_end, ".zst");
      }
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      // reterr = write_map_or_bim(outname, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, '\t', bim_zst);
      if (reterr) {
	goto make_plink2_vsort_ret_1;
      }
      logprint("done.\n");
    }
    if (make_plink2_modifier & kfMakePvar) {
      char* pvarname_end = strcpya0(outname_end, ".pvar");
      if (pvar_psam_modifier & kfPvarZs) {
	strcpy(pvarname_end, ".zst");
      }
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t nonref_flags_storage = 3;
      if (!pgfip->nonref_flags) {
	nonref_flags_storage = (pgfip->gflags & kfPgenGlobalAllNonref)? 2 : 1;
      }
      reterr = write_pvar_resorted(outname, xheader, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, xheader_info_pr, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_modifier);
      if (reterr) {
	goto make_plink2_vsort_ret_1;
      }
      logprint("done.\n");
    }
    if (make_plink2_modifier & kfMakeFam) {
      strcpy(outname_end, ".fam");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_fam(outname, sample_include, sample_ids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, new_sample_idx_to_old, sample_ct, max_sample_id_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, '\t');
      if (reterr) {
	goto make_plink2_vsort_ret_1;
      }
      logprint("done.\n");
    }
    if (make_plink2_modifier & kfMakePsam) {
      strcpy(outname_end, ".psam");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_psam(outname, sample_include, sample_ids, sids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, sample_ct, max_sample_id_blen, max_sid_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, max_pheno_name_blen, pvar_psam_modifier);
      if (reterr) {
	goto make_plink2_vsort_ret_1;
      }
      logprint("done.\n");
    }
  }
  while (0) {
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}
*/

pglerr_t sample_sort_file_map(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* sample_sort_fname, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t sid_col_present, uint32_t** new_sample_idx_to_old_ptr) {
  // assumes sample_ct >= 2 (enforced by caller)
  // return strbox is not collapsed
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    char* idbuf;
    uintptr_t* already_seen;
    if (bigstack_alloc_ui(raw_sample_ct, new_sample_idx_to_old_ptr) ||
	bigstack_alloc_c(max_sample_id_blen, &idbuf) ||
	bigstack_calloc_ul(BITCT_TO_WORDCT(raw_sample_ct), &already_seen)) {
      goto sample_sort_file_map_ret_NOMEM;
    }
    
    uintptr_t loadbuf_size = bigstack_left();
    loadbuf_size -= loadbuf_size / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto sample_sort_file_map_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    char* loadbuf_first_token;
    xid_mode_t xid_mode;
    reterr = open_and_load_xid_header(sample_sort_fname, "indiv-sort", sid_col_present? kSidDetectModeForce : (sids? kSidDetectModeLoaded : kSidDetectModeNotLoaded), loadbuf_size, loadbuf, nullptr, &line_idx, &loadbuf_first_token, &gz_infile, &xid_mode);
    if (reterr) {
      if (reterr == kPglRetEmptyFile) {
	logerrprint("Error: --indiv-sort file is empty.\n");
	goto sample_sort_file_map_ret_MALFORMED_INPUT;
      }
      if (reterr == kPglRetLongLine) {
	if (loadbuf_size == kMaxLongLine) {
	  goto sample_sort_file_map_ret_LONG_LINE;
	}
	goto sample_sort_file_map_ret_NOMEM;
      }
      goto sample_sort_file_map_ret_1;
    }
    uint32_t* xid_map;
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = sorted_xidbox_init_alloc(sample_include, sample_ids, sids, sample_ct, max_sample_id_blen, max_sid_blen, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (reterr) {
      goto sample_sort_file_map_ret_1;
    }
    uint32_t* new_sample_idx_to_old_iter = *new_sample_idx_to_old_ptr;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
	char* loadbuf_iter = loadbuf_first_token;
	uint32_t sample_uidx;
	if (!sorted_xidbox_read_find(sorted_xidbox, xid_map, max_xid_blen, sample_ct, 0, xid_mode, &loadbuf_iter, &sample_uidx, idbuf)) {
	  if (IS_SET(already_seen, sample_uidx)) {
	    char* tptr = (char*)rawmemchr(idbuf, '\t');
	    *tptr = ' ';
	    if (xid_mode & kfXidModeFlagSid) {
	      *((char*)rawmemchr(&(tptr[1]), '\t')) = ' ';
	    }
	    sprintf(g_logbuf, "Error: Duplicate sample ID '%s' in --indiv-sort file.\n", idbuf);
	    goto sample_sort_file_map_ret_MALFORMED_INPUT_WW;
	  }
	  SET_BIT(sample_uidx, already_seen);
	  *new_sample_idx_to_old_iter++ = sample_uidx;
	} else if (!loadbuf_iter) {
	  goto sample_sort_file_map_ret_MISSING_TOKENS;
	}
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto sample_sort_file_map_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto sample_sort_file_map_ret_LONG_LINE;
	}
	goto sample_sort_file_map_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --indiv-sort file starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx);
	goto sample_sort_file_map_ret_MALFORMED_INPUT_WW;
      }
    }
    
    if (gzclose_null(&gz_infile)) {
      goto sample_sort_file_map_ret_READ_FAIL;
    }
    if ((uintptr_t)(new_sample_idx_to_old_iter - (*new_sample_idx_to_old_ptr)) != sample_ct) {
      logerrprint("Error: --indiv-sort file does not contain all loaded sample IDs.\n");
      goto sample_sort_file_map_ret_INCONSISTENT_INPUT;
    }
    bigstack_mark = (unsigned char*)idbuf;
  }
  while (0) {
  sample_sort_file_map_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  sample_sort_file_map_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  sample_sort_file_map_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  sample_sort_file_map_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  sample_sort_file_map_ret_LONG_LINE:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --indiv-sort file is pathologically long.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  sample_sort_file_map_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --indiv-sort file has fewer tokens than expected.\n", line_idx);
  sample_sort_file_map_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 sample_sort_file_map_ret_1:
  gzclose_cond(gz_infile);
  bigstack_reset(bigstack_mark);
  return reterr;
}


// assumes rawval is in [0, 163839]
static_assert(kDosageMid == 16384, "print_small_dosage() needs to be updated.");
char* print_small_dosage(uint32_t rawval, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/16384, (n + 0.5)/16384). 
  // E.g. 3277/16384 is 0.20001 when printed with 5-digit precision, but we'd
  // print that as 0.2 since that's still in (3276.5/16384, 3277.5/16384).
  *start++ = '0' + (rawval / 16384);
  rawval = rawval % 16384;
  if (!rawval) {
    return start;
  }
  *start++ = '.';
  // (rawval * 2) is in 32768ths
  // 32768 * 625 = 20480k

  const uint32_t range_top_20480k = (rawval * 2 + 1) * 625;
  // this is technically checking a half-open rather than a fully-open
  // interval, but that's fine since we never hit the boundary points
  if ((range_top_20480k % 2048) < 1250) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_20480k / 2048;
    return uitoa_trunc4(four_decimal_places, start);
  }
  
  // we wish to print (100000 * remainder + 8192) / 16384, left-0-padded.  and
  // may as well banker's round too.
  //
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4:
  //   1/64 = .015625 -> print 0.01562
  //   3/64 = .046875 -> print 0.04688
  //   5/64 = .078125 -> print 0.07812
  const uint32_t five_decimal_places = ((3125 * rawval + 256) / 512) - ((rawval % 1024) == 256);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return uitoa_trunc4(last_four_digits, start);
  }
  return start;
}

#ifdef __arm__
  #error "Unaligned accesses in export_012_vmaj()."
#endif
pglerr_t export_012_vmaj(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const char* sample_ids, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const double* variant_cms, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t variant_ct, uint32_t max_allele_slen, pgen_reader_t* simple_pgrp) {
  // todo: --recode-allele?
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    const uint32_t max_chr_blen = 1 + get_max_chr_slen(cip);
    char* chr_buf; // includes trailing tab
    char* writebuf;
    uintptr_t* genovec;
    // dosages are limited to 7 characters (x.yyyyy)
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
	bigstack_alloc_c(kMaxMediumLine + max_chr_blen + 2 * kMaxIdSlen + 48 + 2 * max_allele_slen + (8 * k1LU) * sample_ct, &writebuf) ||
        bigstack_alloc_ul(sample_ctl2, &genovec)) {
      goto export_012_vmaj_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t dosage_is_present = simple_pgrp->fi.gflags & kfPgenGlobalDosagePresent;
    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (dosage_is_present) {
      if (bigstack_alloc_ul(sample_ctl, &dosage_present) ||
	  bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
	goto export_012_vmaj_ret_NOMEM;
      }
    }
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_012_vmaj_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya(writebuf, "CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT");
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      *write_iter++ = '\t';
      const char* fid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* fid_end = (const char*)rawmemchr(fid_start, '\t');
      write_iter = memcpyax(write_iter, fid_start, fid_end - fid_start, '_');
      write_iter = strcpya(write_iter, &(fid_end[1]));
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto export_012_vmaj_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
    }
    append_binary_eoln(&write_iter);
    LOGPRINTFWW5("--export A-transpose to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_blen = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t cur_allele_ct = 2;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;

    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	*chr_name_end++ = '\t';
	chr_blen = (uintptr_t)(chr_name_end - chr_buf);
      }
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], '\t');
      if (variant_cms) {
	write_iter = dtoa_g(variant_cms[variant_uidx], write_iter);
      } else {
	*write_iter++ = '0';
      }
      *write_iter++ = '\t';
      write_iter = uint32toa_x(variant_bps[variant_uidx], '\t', write_iter);
      // todo: multiallelic case
      uint32_t dosage_ct;
      uint32_t is_explicit_alt1;
      reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
	if (reterr != kPglRetReadFail) {
	  logprint("\n");
	  logerrprint("Error: Malformed .pgen file.\n");
	}
	goto export_012_vmaj_ret_1;
      }
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[2 * variant_uidx];
      }
      if (!ref_allele_idx) {
	// we *usually* invert, since COUNTED = REF.
	genovec_invert_unsafe(sample_ct, genovec);
	biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], '\t');
      const uint32_t first_alt_idx = (ref_allele_idx == 0);
      write_iter = strcpya(write_iter, cur_alleles[first_alt_idx]);
      if (cur_allele_ct > 2) {
	for (uint32_t allele_idx = first_alt_idx + 1; allele_idx < cur_allele_ct; ++allele_idx) {
	  if (allele_idx == ref_allele_idx) {
	    continue;
	  }
	  if (write_iter >= writebuf_flush) {
	    if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	      goto export_012_vmaj_ret_WRITE_FAIL;
	    }
	    write_iter = writebuf;
	  }
	  *write_iter++ = ',';
	  write_iter = strcpya(write_iter, cur_alleles[allele_idx]);
	}
      }
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto export_012_vmaj_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
      uint32_t widx = 0;
      uint32_t loop_len = kBitsPerWordD2;
      if (!dosage_ct) {
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
	  }
	  uintptr_t geno_word = genovec[widx];
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits < loop_len; ++sample_idx_lowbits) {
	    *write_iter++ = '\t';
	    uintptr_t cur_geno = geno_word & 3;
	    if (cur_geno != 3) {
	      *write_iter++ = '0' + cur_geno;
	    } else {
	      write_iter = strcpya(write_iter, "NA");
	    }
	    geno_word >>= 2;
	  }
	  ++widx;
	}
      } else {
	dosage_t* dosage_vals_iter = dosage_vals;
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
	  }
	  uintptr_t geno_word = genovec[widx];
	  uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits < loop_len; ++sample_idx_lowbits) {
	    *write_iter++ = '\t';
	    if (dosage_present_hw & 1) {
	      write_iter = print_small_dosage(*dosage_vals_iter++, write_iter);
	    } else {
	      uintptr_t cur_geno = geno_word & 3;
	      if (cur_geno != 3) {
		*write_iter++ = '0' + cur_geno;
	      } else {
		write_iter = strcpya(write_iter, "NA");
	      }
	    }
	    geno_word >>= 2;
	    dosage_present_hw >>= 1;
	  }
	  ++widx;
	}
      }
      append_binary_eoln(&write_iter);
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
    if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
      goto export_012_vmaj_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile)) {
      goto export_012_vmaj_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  export_012_vmaj_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_012_vmaj_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_012_vmaj_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 export_012_vmaj_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

static uintptr_t* g_vmaj_readbuf = nullptr;

THREAD_FUNC_DECL transpose_to_smaj_read_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  // const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uint32_t read_sample_ct = g_sample_ct;
  const uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
  uintptr_t prev_copy_ct = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_copy_ct = g_cur_block_write_ct;
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_copy_ct) / calc_thread_ct;
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    uint32_t cur_idx = (tidx * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t* vmaj_readbuf_iter = &(g_vmaj_readbuf[(prev_copy_ct + cur_idx) * read_sample_ctaw2]);
    for (; cur_idx < cur_idx_end; ++cur_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      // todo: multiallelic case
      const pglerr_t reterr = pgr_read_refalt1_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, read_sample_ct, variant_uidx, pgrp, vmaj_readbuf_iter);
      if (reterr) {
	g_error_ret = reterr;
	break;
      }
      if (refalt1_select && (refalt1_select[2 * variant_uidx] == 1)) {
	genovec_invert_unsafe(read_sample_ct, vmaj_readbuf_iter);
	// don't need zero_trailing_quaters()
      }
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[read_sample_ctaw2]);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    prev_copy_ct += cur_block_copy_ct;
    THREAD_BLOCK_FINISH(tidx);
  }
}

static uintptr_t* g_smaj_writebufs[2] = {nullptr, nullptr};
static uint32_t g_variant_ct = 0;
static uint32_t g_sample_batch_size = 0;
static uint32_t g_output_calc_thread_ct = 0;

THREAD_FUNC_DECL transpose_to_plink1_smaj_write_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t variant_ct = g_variant_ct;
  const uintptr_t variant_batch_ct = DIV_UP(variant_ct, kPglQuaterTransposeBatch);
  const uintptr_t variant_batch_word_ct = variant_batch_ct * kPglQuaterTransposeWords;
  const uint32_t calc_thread_ct = g_output_calc_thread_ct;
  const uint32_t variant_batch_idx_start = (((uint64_t)tidx) * variant_batch_ct) / calc_thread_ct;
  vul_t* vecaligned_buf = g_thread_vecaligned_bufs[tidx];
  uintptr_t variant_batch_idx_full_end = ((((uint64_t)tidx) + 1) * variant_batch_ct) / calc_thread_ct;
  uint32_t variant_idx_end;
  if (tidx + 1 < calc_thread_ct) {
    variant_idx_end = variant_batch_idx_full_end * kPglQuaterTransposeBatch;
  } else {
    variant_idx_end = variant_ct;
    if (variant_ct % kPglQuaterTransposeBatch) {
      --variant_batch_idx_full_end;
    }
  }
  const uint32_t thread_variant_ct = variant_idx_end - variant_batch_idx_start * kPglQuaterTransposeBatch;
  const uint32_t read_sample_ct = g_sample_ct;
  const uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
  const uintptr_t* vmaj_readbuf = g_vmaj_readbuf;
  uint32_t sample_widx = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    uintptr_t variant_batch_idx = variant_batch_idx_start;
    uint32_t variant_batch_size = kPglQuaterTransposeBatch;
    const uintptr_t* vmaj_readbuf_iter = &(vmaj_readbuf[variant_batch_idx * kPglQuaterTransposeBatch * read_sample_ctaw2 + sample_widx]);
    const uint32_t sample_batch_size = g_sample_batch_size;
    uintptr_t* smaj_writebuf_start = &(g_smaj_writebufs[parity][variant_batch_idx * kPglQuaterTransposeWords]);
    uintptr_t* smaj_writebuf_iter = smaj_writebuf_start;
    while (1) {
      if (variant_batch_idx >= variant_batch_idx_full_end) {
	if (variant_batch_idx * kPglQuaterTransposeBatch >= variant_idx_end) {
	  break;
	}
	variant_batch_size = variant_idx_end - variant_batch_idx * kPglQuaterTransposeBatch;
      }
      transpose_quaterblock(vmaj_readbuf_iter, read_sample_ctaw2, variant_batch_word_ct, variant_batch_size, sample_batch_size, smaj_writebuf_iter, vecaligned_buf);
      smaj_writebuf_iter = &(smaj_writebuf_iter[kPglQuaterTransposeWords]);
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[variant_batch_size * read_sample_ctaw2]);
      ++variant_batch_idx;
    }
    smaj_writebuf_iter = smaj_writebuf_start;
    for (uint32_t sample_idx = 0; sample_idx < sample_batch_size; ++sample_idx) {
      // could fold this into transpose_quaterblock(), but I won't bother,
      // we're already saturating at ~3 threads
      pgr_plink2_to_plink1_inplace_unsafe(thread_variant_ct, smaj_writebuf_iter);
      zero_trailing_quaters(thread_variant_ct, smaj_writebuf_iter);
      smaj_writebuf_iter = &(smaj_writebuf_iter[variant_batch_word_ct]);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
    sample_widx += sample_batch_size / kBitsPerWordD2;
  }
}

pglerr_t export_ind_major_bed(const uintptr_t* orig_sample_include, const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const alt_allele_ct_t* refalt1_select, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    // Possible special case: if the input file is a variant-major .bed, we do
    // not have enough memory to just load the whole file at once, and there
    // are more than ~20k samples, there can be a performance advantage to not
    // loading an entire variant at a time; we can use smaller fread calls and
    // reduce the number of (typically 4096 byte) disk blocks which need to be
    // read on each pass.  But let's get .pgen -> sample-major humming first.
    strcpy(outname_end, ".bed");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ind_major_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\0", 3, outfile)) {
      goto export_ind_major_bed_ret_WRITE_FAIL;
    }
    if (variant_ct && sample_ct) {
      const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
      uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      // todo: if only 1 pass is needed, and no subsetting is happening, this
      // saturates at ~4 threads?
      unsigned char* bigstack_end_mark = g_bigstack_end;
      // restrict multithread_load_init() to half of available workspace
      g_bigstack_end = &(g_bigstack_base[round_up_pow2(bigstack_left() / 2, kCacheline)]);
      unsigned char* main_loadbufs[2];
      pthread_t* threads;
      uint32_t read_block_size;
      if (multithread_load_init(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, nullptr, nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
	goto export_ind_major_bed_ret_NOMEM;
      }
      g_bigstack_end = bigstack_end_mark;
      g_variant_include = variant_include;
      g_variant_allele_idxs = variant_allele_idxs;
      g_refalt1_select = refalt1_select;
      g_calc_thread_ct = calc_thread_ct;
      
      const uintptr_t variant_cacheline_ct = QUATERCT_TO_CLCT(variant_ct);
      uint32_t output_calc_thread_ct = MINV(calc_thread_ct, variant_cacheline_ct);
      if (output_calc_thread_ct > 4) {
	output_calc_thread_ct = 4;
      }
      uintptr_t* sample_include;
      uint32_t* sample_include_cumulative_popcounts;
      if (bigstack_alloc_ul(raw_sample_ctl, &sample_include) ||
	  bigstack_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
	  bigstack_alloc_vp(output_calc_thread_ct, &g_thread_vecaligned_bufs)) {
	goto export_ind_major_bed_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx < output_calc_thread_ct; ++tidx) {
	g_thread_vecaligned_bufs[tidx] = (vul_t*)bigstack_alloc_raw(kPglQuaterTransposeBufbytes);
      }
      // each of the two write buffers should use <= 1/8 of the remaining
      // workspace
      const uintptr_t writebuf_cachelines_avail = bigstack_left() / (kCacheline * 8);
      uint32_t sample_batch_size = kPglQuaterTransposeBatch;
      if (variant_cacheline_ct * kPglQuaterTransposeBatch > writebuf_cachelines_avail) {
	sample_batch_size = round_down_pow2(writebuf_cachelines_avail / variant_cacheline_ct, kBitsPerWordD2);
	if (!sample_batch_size) {
	  goto export_ind_major_bed_ret_NOMEM;
	}
      }
      g_smaj_writebufs[0] = (uintptr_t*)bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size);
      g_smaj_writebufs[1] = (uintptr_t*)bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size);
      const uintptr_t readbuf_vecs_avail = (bigstack_left() / kCacheline) * kVecsPerCacheline;
      if (readbuf_vecs_avail < variant_ct) {
	goto export_ind_major_bed_ret_NOMEM;
      }
      uintptr_t read_sample_ctv2 = readbuf_vecs_avail / variant_ct;
      uint32_t read_sample_ct;
      if (read_sample_ctv2 >= QUATERCT_TO_VECCT(sample_ct)) {
	read_sample_ct = sample_ct;
      } else {
	read_sample_ct = read_sample_ctv2 * kQuatersPerVec;
      }
      uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
      uintptr_t* vmaj_readbuf = (uintptr_t*)bigstack_alloc_raw_rd(variant_ct * read_sample_ctaw2 * kBytesPerWord);
      g_variant_ct = variant_ct;
      g_output_calc_thread_ct = output_calc_thread_ct;
      g_error_ret = kPglRetSuccess;
      uint32_t sample_uidx_start = next_set_unsafe(orig_sample_include, 0);
      const uintptr_t variant_ct4 = QUATERCT_TO_BYTECT(variant_ct);
      const uintptr_t variant_ctaclw2 = variant_cacheline_ct * kWordsPerCacheline;
      const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
      const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
      const uint32_t pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
      for (uint32_t pass_idx = 0; pass_idx < pass_ct; ++pass_idx) {
	memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
	if (sample_uidx_start) {
	  clear_bits_nz(0, sample_uidx_start, sample_include);
	}
	uint32_t sample_uidx_end;
	if (pass_idx + 1 == pass_ct) {
	  read_sample_ct = sample_ct - pass_idx * read_sample_ct;
	  read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
	  sample_uidx_end = raw_sample_ct;
	} else {
	  sample_uidx_end = jump_forward_set_unsafe(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
	  clear_bits_nz(sample_uidx_end, raw_sample_ct, sample_include);
	}
        fill_cumulative_popcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
	g_sample_include = sample_include;
	g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
	g_vmaj_readbuf = vmaj_readbuf;
	g_sample_ct = read_sample_ct;
	if (pass_idx) {
	  pgfip->block_base = main_loadbufs[0];
	  for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	    pgr_clear_ld_cache(g_pgr_ptrs[tidx]);
	    g_pgr_ptrs[tidx]->fi.block_base = main_loadbufs[0];
	    g_pgr_ptrs[tidx]->fi.block_offset = 0;
	  }
	}
	uint32_t parity = 0;
	uint32_t read_block_idx = 0;
        uint32_t variant_idx = 0;
	uint32_t is_last_block = 0;
	uint32_t cur_read_block_size = read_block_size;
	uint32_t pct = 0;
	uint32_t next_print_idx = variant_ct / 100;
	putc_unlocked('\r', stdout);
	printf("--export ind-major-bed pass %u/%u: loading... 0%%", pass_idx + 1, pass_ct);
	fflush(stdout);
	while (1) {
	  uintptr_t cur_block_write_ct = 0;
	  if (!is_last_block) {
	    while (read_block_idx < read_block_ct_m1) {
	      cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	      if (cur_block_write_ct) {
		break;
	      }
	      ++read_block_idx;
	    }
	    if (read_block_idx == read_block_ct_m1) {
	      cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	      cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	    }
	    if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
	      if (variant_idx) {
		join_threads2z(calc_thread_ct, 0, threads);
		g_cur_block_write_ct = 0;
		error_cleanup_threads2z(transpose_to_smaj_read_thread, calc_thread_ct, threads);
	      }
	      goto export_ind_major_bed_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  if (variant_idx) {
	    join_threads2z(calc_thread_ct, is_last_block, threads);
	    reterr = g_error_ret;
	    if (reterr) {
	      if (!is_last_block) {
		g_cur_block_write_ct = 0;
		error_cleanup_threads2z(transpose_to_smaj_read_thread, calc_thread_ct, threads);
	      }
	      if (reterr == kPglRetMalformedInput) {
		logprint("\n");
		logerrprint("Error: Malformed .pgen file.\n");
	      }
	      goto export_ind_major_bed_ret_1;
	    }
	  }
	  if (!is_last_block) {
	    g_cur_block_write_ct = cur_block_write_ct;
	    compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
	    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	      g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	      g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	    }
	    is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
	    if (spawn_threads2z(transpose_to_smaj_read_thread, calc_thread_ct, is_last_block, threads)) {
	      goto export_ind_major_bed_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  parity = 1 - parity;
	  if (variant_idx == variant_ct) {
	    break;
	  }
	  if (variant_idx >= next_print_idx) {
	    if (pct > 10) {
	      putc_unlocked('\b', stdout);
	    }
	    pct = (variant_idx * 100LLU) / variant_ct;
	    printf("\b\b%u%%", pct++);
	    fflush(stdout);
	    next_print_idx = (pct * ((uint64_t)variant_ct)) / 100;
	  }

	  ++read_block_idx;
	  variant_idx += cur_block_write_ct;
	  pgfip->block_base = main_loadbufs[parity];
	}
	// 2. Transpose and write.  (Could parallelize some of the transposing
	//    with the read loop, but since we can't write a single row until
	//    the read loop is done, and both write speed and write buffer
	//    space are bottlenecks, that can't be expected to help much.)
	g_sample_batch_size = sample_batch_size;
	parity = 0;
	is_last_block = 0;
	if (pct > 10) {
	  fputs("\b \b", stdout);
	}
	fputs("\b\b\b\b\b\b\b\b\b\b\b\b\bwriting... 0%", stdout);
	fflush(stdout);
	pct = 0;
	uint32_t flush_sample_idx = 0;
	uint32_t flush_sample_idx_end = 0;
	next_print_idx = read_sample_ct / 100;
	while (1) {
	  if (!is_last_block) {
	    is_last_block = (flush_sample_idx_end + sample_batch_size >= read_sample_ct);
	    if (is_last_block) {
	      g_sample_batch_size = read_sample_ct - flush_sample_idx_end;
	    }
	    if (spawn_threads2z(transpose_to_plink1_smaj_write_thread, output_calc_thread_ct, is_last_block, threads)) {
	      goto export_ind_major_bed_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  if (flush_sample_idx_end) {
	    uintptr_t* smaj_writebuf_iter = g_smaj_writebufs[1 - parity];
	    for (; flush_sample_idx < flush_sample_idx_end; ++flush_sample_idx) {
	      fwrite(smaj_writebuf_iter, variant_ct4, 1, outfile);
	      smaj_writebuf_iter = &(smaj_writebuf_iter[variant_ctaclw2]);
	    }
	    if (flush_sample_idx_end == read_sample_ct) {
	      break;
	    }
	    if (flush_sample_idx_end >= next_print_idx) {
	      if (pct > 10) {
		putc_unlocked('\b', stdout);
	      }
	      pct = (flush_sample_idx_end * 100LLU) / read_sample_ct;
	      printf("\b\b%u%%", pct++);
	      fflush(stdout);
	      next_print_idx = (pct * ((uint64_t)read_sample_ct)) / 100;
	    }
	  }
	  join_threads2z(output_calc_thread_ct, is_last_block, threads);
	  if (ferror(outfile)) {
	    // may as well put this check when there are no threads to clean up
	    goto export_ind_major_bed_ret_WRITE_FAIL;
	  }
	  parity = 1 - parity;
	  flush_sample_idx_end += sample_batch_size;
	  if (flush_sample_idx_end > read_sample_ct) {
	    flush_sample_idx_end = read_sample_ct;
	  }
	}
	if (pct > 10) {
	  fputs("\b \b", stdout);
	}
	sample_uidx_start = sample_uidx_end;
      }
      fputs("\b\bdone.\n", stdout);
    }
    if (fclose_null(&outfile)) {
      goto export_ind_major_bed_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--export ind-major-bed: %s written.\n", outname);
  }
  while (0) {
  export_ind_major_bed_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ind_major_bed_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ind_major_bed_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_ind_major_bed_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_ind_major_bed_ret_1:
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  bigstack_reset(bigstack_mark);
  return reterr;
}

static_assert(kDosageMid == 16384, "print_gen_dosage() needs to be updated.");
char* print_gen_dosage(uint32_t rawval, char* start) {
  // Similar to small_dosage_print(), but it's complicated by .gen import's
  // quantization step (instead of rounding the numbers directly, they're first
  // converted to bgen-1.1 equivalents).  We check
  //   ((n - 0.75)/16384, (n + 0.75)/16384) for even n
  //   ((n - 0.25)/16384, (n + 0.25)/16384) for odd n
  // due to banker's rounding.
  *start++ = '0' + (rawval / 16384);
  rawval = rawval % 16384;
  if (!rawval) {
    return start;
  }
  *start++ = '.';
  const uint32_t halfwidth_65536ths = 3 - (2 * (rawval % 2));
  // (rawval * 4) is in 65536ths
  // 65536 * 625 = 40960k

  const uint32_t range_top_40960k = (rawval * 4 + halfwidth_65536ths) * 625;
  // this is technically checking a half-open rather than a fully-open
  // interval, but that's fine since we never hit the boundary points
  if ((range_top_40960k % 4096) < 1250 * halfwidth_65536ths) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_40960k / 4096;
    return uitoa_trunc4(four_decimal_places, start);
  }
  
  // we wish to print (100000 * remainder + 8192) / 16384, left-0-padded.  and
  // may as well banker's round too.
  //  
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4:
  //   1/64 = .015625 -> print 0.01562
  //   3/64 = .046875 -> print 0.04688
  //   5/64 = .078125 -> print 0.07812
  const uint32_t five_decimal_places = ((3125 * rawval + 256) / 512) - ((rawval % 1024) == 256);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return uitoa_trunc4(last_four_digits, start);
  }
  return start;
}

// note that this is also in plink2_filter; may belong in plink2_common
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

pglerr_t export_ox_gen(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t variant_ct, uint32_t max_allele_slen, exportf_flags_t exportf_modifier, pgen_reader_t* simple_pgrp, uint32_t* sample_missing_geno_cts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* sex_male_collapsed_interleaved;
    uintptr_t* sex_male_collapsed_tmp;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
	
	bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_interleaved) ||
	bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_tmp)) {
      goto export_ox_gen_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed_tmp);
    fill_interleaved_mask_vec(sex_male_collapsed_tmp, sample_ctv, sex_male_collapsed_interleaved);
    bigstack_reset(sex_male_collapsed_tmp);
    
    // See load_sample_missing_cts in plink2_filter.cpp.
    // Yes, this is overkill, but the obvious alternative of incrementing
    // sample_missing_geno_cts[] when writing a missing call requires a bit of
    // custom chrY logic anyway.
    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
    uintptr_t* missing_acc2;
    if (bigstack_calloc_ul(acc2_vec_ct * kWordsPerVec * 23, &missing_acc2)) {
      goto export_ox_gen_ret_NOMEM;
    }
    const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
    const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
    uintptr_t* missing_acc4 = &(missing_acc2[acc2_vec_ct * kWordsPerVec]);
    uintptr_t* missing_acc8 = &(missing_acc4[acc4_vec_ct * kWordsPerVec]);
    uintptr_t* missing_acc32 = &(missing_acc8[acc8_vec_ct * kWordsPerVec]);

    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (simple_pgrp->fi.gflags & kfPgenGlobalDosagePresent) {
      const uint32_t multiallelic_present = (variant_allele_idxs != nullptr);
      if (bigstack_alloc_dosage(sample_ct * (1 + multiallelic_present), &dosage_vals) ||
	  bigstack_alloc_ul(sample_ctl, &dosage_present)) {
	goto export_ox_gen_ret_NOMEM;
      }
    }
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    // if no dosages, all genotypes are 6 bytes (missing = " 0 0 0")
    // with dosages, we print up to 5 digits past the decimal point, so 7 bytes
    //   + space for each number, 24 bytes max
    const uintptr_t max_geno_slen = 6 + (dosage_present != nullptr) * 18;
    char* chr_buf; // includes trailing space
    char* writebuf;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
	bigstack_alloc_c(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 16 + 2 * max_allele_slen + max_geno_slen * sample_ct, &writebuf)) {
      goto export_ox_gen_ret_NOMEM;
    }
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ox_gen_ret_OPEN_FAIL;
    }
    char* writebuf_flush = &(writebuf[kCompressStreamBlock]);
    char* write_iter = writebuf;
    uint32_t variant_uidx = 0;
    uint32_t chr_blen = 0;

    // although we don't support --set-hh-missing, etc. here, we do still want
    // to be aware of chrY so we can exclude nonmales from the
    // sample_missing_geno_cts update there.
    uint32_t is_y = 0;

    uint32_t chr_fo_idx = 0xffffffffU;
    const int32_t y_code = cip->xymt_codes[kChrOffsetY];
    uint32_t chr_end = 0;
    uint32_t vidx_rem3 = 3;
    uint32_t vidx_rem15d3 = 5;
    uint32_t vidx_rem255d15 = 17;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const char hardcall_strs[] = " 1 0 0   0 1 0   0 0 1   0 0 0";
    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	// Oxford spec doesn't seem to require spaces for .gen (only .sample),
	// but in practice spaces always seem to be used, and plink 1.9 doesn't
	// let you toggle this, so let's not worry about supporting tabs here
	*chr_name_end++ = ' ';
	chr_blen = (uintptr_t)(chr_name_end - chr_buf);
	is_y = (chr_idx == y_code);
      }
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
      write_iter = uint32toa_x(variant_bps[variant_uidx], ' ', write_iter);
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      uint32_t is_explicit_alt1;
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
	alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
      }
      // todo: multiallelic case
      uint32_t dosage_ct;
      reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
	if (reterr != kPglRetReadFail) {
	  logprint("\n");
	  logerrprint("Error: Malformed .pgen file.\n");
	}
	goto export_ox_gen_ret_1;
      }
      if (ref_allele_idx + ref_allele_second == 1) {
	assert((!dosage_ct) || (!is_explicit_alt1));
	genovec_invert_unsafe(sample_ct, genovec);
	biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      if (ref_allele_second) {
	write_iter = strcpyax(write_iter, cur_alleles[alt1_allele_idx], ' ');
	write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
      } else {
	write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
	write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
      }
      uint32_t widx = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      if (!dosage_ct) {
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	  }
	  uintptr_t geno_word = genovec[widx];
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    write_iter = memcpya(write_iter, &(hardcall_strs[(geno_word & 3) * 8]), 6);
	    geno_word >>= 2;
	  }
	  ++widx;
	}
      } else {
	const halfword_t* dosage_present_alias = (halfword_t*)dosage_present;
	const dosage_t* dosage_vals_iter = dosage_vals;
	if (!is_explicit_alt1) {
	  while (1) {
	    if (widx >= sample_ctl2_m1) {
	      if (widx > sample_ctl2_m1) {
		break;
	      }
	      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	    }
	    uintptr_t geno_word = genovec[widx];
	    uint32_t dosage_present_hw = dosage_present_alias[widx];
	    if (!dosage_present_hw) {
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		write_iter = memcpya(write_iter, &(hardcall_strs[(geno_word & 3) * 8]), 6);
		geno_word >>= 2;
	      }
	    } else {
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		if (dosage_present_hw & 1) {
		  const uint32_t dosage_int = *dosage_vals_iter++;
		  if (dosage_int <= kDosageMid) {
		    *write_iter++ = ' ';
		    write_iter = print_gen_dosage(kDosageMid - dosage_int, write_iter);
		    *write_iter++ = ' ';
		    write_iter = print_gen_dosage(dosage_int, write_iter);
		    write_iter = strcpya(write_iter, " 0");
		  } else {
		    assert(dosage_int <= kDosageMax);
		    write_iter = memcpyl3a(write_iter, " 0 ");
		    write_iter = print_gen_dosage(kDosageMax - dosage_int, write_iter);
		    *write_iter++ = ' ';
		    write_iter = print_gen_dosage(dosage_int - kDosageMid, write_iter);
		  }
		} else {
		  write_iter = memcpya(write_iter, &(hardcall_strs[(geno_word & 3) * 8]), 6);
		}
		geno_word >>= 2;
		dosage_present_hw >>= 1;
	      }
	    }
	    ++widx;
	  }
	} else {
	  // todo
	  // In multiallelic case, if ref/alt1 dosages sum to less than 2 (but
	  // more than 0), we first internally rescale them to sum to 2, to
	  // make .gen and bgen-1.1 export isomorphic, and bgen-1.2 export as
	  // similar as possible.
	  assert(0);
	}
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto export_ox_gen_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
      if (is_y) {
	interleaved_mask_zero(sex_male_collapsed_interleaved, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, missing_acc2);
      if (!(--vidx_rem3)) {
	unroll_zero_incr_2_4(acc2_vec_ct, missing_acc2, missing_acc4);
	vidx_rem3 = 3;
	if (!(--vidx_rem15d3)) {
	  unroll_zero_incr_4_8(acc4_vec_ct, missing_acc4, missing_acc8);
	  vidx_rem15d3 = 5;
	  if (!(--vidx_rem255d15)) {
	    unroll_zero_incr_8_32(acc8_vec_ct, missing_acc8, missing_acc32);
	    vidx_rem255d15 = 17;
	  }
	}
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
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	goto export_ox_gen_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto export_ox_gen_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    unroll_incr_2_4(missing_acc2, acc2_vec_ct, missing_acc4);
    unroll_incr_4_8(missing_acc4, acc4_vec_ct, missing_acc8);
    unroll_incr_8_32(missing_acc8, acc8_vec_ct, missing_acc32);
    uint32_t* scrambled_missing_cts = (uint32_t*)missing_acc32;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  export_ox_gen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ox_gen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ox_gen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 export_ox_gen_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __arm__
  #error "Unaligned accesses in export_ox_hapslegend()."
#endif
pglerr_t export_ox_hapslegend(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, exportf_flags_t exportf_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  assert(sample_ct);
  assert(variant_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    const uint32_t just_haps = (exportf_modifier / kfExportfHaps) & 1;
    const uint32_t male_ct = popcount_longs(sex_male_collapsed, sample_ctl);
    if (xymt_is_nonempty(variant_include, cip, kChrOffsetY) && (male_ct != sample_ct)) {
      LOGERRPRINTF("Error: '--export haps%s' must exclude chrY unless the dataset is all-male.\n", just_haps? "" : "legend");
      goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
    }
    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    const int32_t x_code = cip->xymt_codes[kChrOffsetX];
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    char* chr_buf = nullptr;
    uint32_t is_x = 0;
    uint32_t is_haploid_or_mt = 0;
    uint32_t variant_uidx = next_set_unsafe(variant_include, 0);
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uintptr_t writebuf_alloc = 0;
    if (!just_haps) {
      // .legend doesn't have a chromosome column, so verify we only need to
      // export a single chromosome
      const uint32_t variant_uidx_start = variant_uidx;
      chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx_start);
      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      if ((chr_end != raw_variant_ct) && (popcount_bit_idx(variant_include, variant_uidx_start, chr_end) != variant_ct)) {
	logerrprint("Error: '--export hapslegend' does not support multiple chromosomes.\n");
	goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
      }
      const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      is_x = (chr_idx == x_code);
      is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == mt_code);
      strcpy(outname_end, ".legend");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto export_ox_hapslegend_ret_OPEN_FAIL;
      }
      char* writebuf;
      if (bigstack_alloc_c(kCompressStreamBlock + kMaxIdSlen + 32 + 2 * max_allele_slen, &writebuf)) {
	goto export_ox_hapslegend_ret_NOMEM;
      }
      char* writebuf_flush = &(writebuf[kCompressStreamBlock]);
      char* write_iter = strcpya(writebuf, "id position a0 a1" EOLN_STR);
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
	next_set_unsafe_ck(variant_include, &variant_uidx);
	write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
	write_iter = uint32toa_x(variant_bps[variant_uidx], ' ', write_iter);
	if (refalt1_select) {
	  ref_allele_idx = refalt1_select[variant_uidx * 2];
	  alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
	}
	uintptr_t variant_allele_idx_base = variant_uidx * 2;
	if (variant_allele_idxs) {
	  variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	}
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
        if (ref_allele_second) {
	  write_iter = strcpyax(write_iter, cur_alleles[alt1_allele_idx], ' ');
	  write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
	} else {
	  write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
	  write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
	}
	append_binary_eoln(&write_iter);
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	    goto export_ox_hapslegend_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
      }
      if (write_iter != writebuf) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto export_ox_hapslegend_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto export_ox_hapslegend_ret_WRITE_FAIL;
      }
      logprint("done.\n");
      variant_uidx = variant_uidx_start;
      bigstack_reset(writebuf);
    } else {
      const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
      if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
	goto export_ox_hapslegend_ret_NOMEM;
      }
      writebuf_alloc = max_chr_blen + kMaxIdSlen + 32 + 2 * max_allele_slen;
    }
    writebuf_alloc += kCompressStreamBlock + (4 * k1LU) * sample_ct + kCacheline;
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    char* writebuf;
    uintptr_t* sex_male_collapsed_interleaved;
    uintptr_t* genovec;
    uintptr_t* phasepresent;
    uintptr_t* phaseinfo;
    if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_interleaved) ||
	bigstack_alloc_c(writebuf_alloc, &writebuf) ||
	bigstack_alloc_ul(sample_ctl2, &genovec) ||
        bigstack_alloc_ul(sample_ctl, &phasepresent) ||
	bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
      goto export_ox_hapslegend_ret_NOMEM;
    }
    fill_interleaved_mask_vec(sex_male_collapsed, sample_ctv, sex_male_collapsed_interleaved);
    // assumes little-endian
    // 3 = 1|0, not missing
    // 4..7 = male chrX
    // user's responsibility to split off PARs
    uint32_t genotext[7];
    genotext[0] = 0x20302030;
    genotext[2] = 0x20312031;
    genotext[4] = 0x202d2030;
    genotext[6] = 0x202d2031;
    if (ref_allele_second) {
      genotext[1] = 0x20302031;
      genotext[3] = 0x20312030;
    } else {
      genotext[1] = 0x20312030;
      genotext[3] = 0x20302031;
    }
#ifndef NDEBUG
    genotext[5] = 0x21475542; // "BUG!"
#endif
    uint32_t* cur_genotext = genotext;
    if (is_haploid_or_mt && (!is_x)) {
      cur_genotext = &(genotext[4]);
    }
    char* writebuf_flush = &(writebuf[kCompressStreamBlock]);
    char* write_iter = writebuf;
    strcpy(outname_end, ".haps");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ox_hapslegend_ret_OPEN_FAIL;
    }
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t chr_blen = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	*chr_name_end++ = ' ';
	chr_blen = (uintptr_t)(chr_name_end - chr_buf);
	is_x = (chr_idx == x_code);
	is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == mt_code);
	if ((!is_haploid_or_mt) || is_x) {
	  cur_genotext = genotext;
	} else {
	  cur_genotext = &(genotext[4]);
	}
      }
      uint32_t phasepresent_ct;
      reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct);
      if (reterr) {
	goto export_ox_hapslegend_ret_PGR_FAIL;
      }
      zero_trailing_quaters(sample_ct, genovec);
      if (!phasepresent_ct) {
	// phaseinfo is NOT cleared in this case
	fill_ulong_zero(sample_ctl, phaseinfo);
      }
      uint32_t genocounts[4];
      genovec_count_freqs_unsafe(genovec, sample_ct, genocounts);
      if (phasepresent_ct != genocounts[1]) {
	logprint("\n");
	LOGERRPRINTF("Error: '--export haps%s' must be used with a fully phased dataset.\n", just_haps? "" : "legend");
	goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
      } else if (genocounts[3]) {
	logprint("\n");
	LOGERRPRINTF("Error: '--export haps%s' cannot be used with missing genotype calls.\n", just_haps? "" : "legend");
	goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
      }
      if (is_haploid_or_mt) {
	// verify that there are no het haploids (treating MT as haploid here)
	if (is_x) {
	  genovec_count_subset_freqs(genovec, sex_male_collapsed_interleaved, sample_ct, male_ct, genocounts);
	}
	if (genocounts[1]) {
	  logprint("\n");
	  LOGERRPRINTFWW("Error: '--export haps%s' cannot be used when heterozygous haploid/MT calls are present.%s\n", just_haps? "" : "legend", (is_x && (variant_bps[variant_uidx] <= 2781479))? " (Did you forget --split-par?)" : "");
	  goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
	}
      }
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
	alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
      }
      // this logic only works in the biallelic case
      if (ref_allele_second + ref_allele_idx == 1) {
	genovec_invert_unsafe(sample_ct, genovec);
	zero_trailing_quaters(sample_ct, genovec);
      }
      if (just_haps) {
	write_iter = memcpya(write_iter, chr_buf, chr_blen);
	write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
	write_iter = uint32toa_x(variant_bps[variant_uidx], ' ', write_iter);
	char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
        if (ref_allele_second) {
	  write_iter = strcpyax(write_iter, cur_alleles[alt1_allele_idx], ' ');
	  write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
	} else {
	  write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
	  write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
	}
	*write_iter++ = ' ';
      }
      uint32_t* write_iter_ui_alias = (uint32_t*)write_iter;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      if (!is_x) {
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	  }
	  uintptr_t genovec_word = genovec[widx];
	  const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    const uintptr_t cur_geno = genovec_word & 3;
	    *write_iter_ui_alias++ = cur_genotext[cur_geno + 2 * ((phaseinfo_halfword >> sample_idx_lowbits) & 1)];
	    genovec_word >>= 2;
	  }
	  ++widx;
	}
      } else {
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	  }
	  uintptr_t genovec_word = genovec[widx];
	  const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
	  const uint32_t male_halfword = ((const halfword_t*)sex_male_collapsed)[widx];
	  
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    const uintptr_t cur_geno = genovec_word & 3;
	    if (cur_geno == 2) {
	      assert(!((phaseinfo_halfword >> sample_idx_lowbits) & 1));
	    }
	    *write_iter_ui_alias++ = cur_genotext[cur_geno + 2 * ((phaseinfo_halfword >> sample_idx_lowbits) & 1) + 4 * ((male_halfword >> sample_idx_lowbits) & 1)];
	    genovec_word >>= 2;
	  }
	  ++widx;
	}
      }
      write_iter = (char*)write_iter_ui_alias;
      decr_append_binary_eoln(&write_iter);
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto export_ox_hapslegend_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
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
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	goto export_ox_hapslegend_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto export_ox_hapslegend_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  export_ox_hapslegend_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ox_hapslegend_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ox_hapslegend_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_ox_hapslegend_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  export_ox_hapslegend_ret_PGR_FAIL:
    if (reterr != kPglRetReadFail) {
      logprint("\n");
      logerrprint("Error: Malformed .pgen file.\n");
    }
  }
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

static uintptr_t** g_missing_acc2 = nullptr;
static uint32_t* g_variant_bytects[2] = {nullptr, nullptr};
static uint32_t g_ref_allele_second = 0;
static uint32_t g_bgen_compressed_buf_max = 0;
static uint32_t g_y_start = 0;
static uint32_t g_y_end = 0;

static const uint16_t bgen11_hardcall_usis[] = {32768, 0, 0, 0,
						0, 32768, 0, 0,
						0, 0, 32768, 0,
						0, 0, 0, 0};

THREAD_FUNC_DECL export_bgen11_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
  const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
  const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
  uintptr_t* missing_acc2 = g_missing_acc2[tidx];
  uintptr_t* missing_acc4 = &(missing_acc2[acc2_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc8 = &(missing_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc32 = &(missing_acc8[acc8_vec_ct * kWordsPerVec]);
  uintptr_t* dosage_present = g_dosage_presents? g_dosage_presents[tidx] : nullptr;
  dosage_t* dosage_vals = dosage_present? g_dosage_val_bufs[tidx] : nullptr;
  uint16_t* bgen_geno_buf = g_bgen_geno_bufs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t sample_ctl2_m1 = QUATERCT_TO_WORDCT(sample_ct) - 1;
  const uint32_t bgen_geno_buf_blen = 6 * sample_ct;
  const uint32_t bgen_compressed_buf_max = g_bgen_compressed_buf_max;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  uint32_t is_y = 0;
  uint32_t y_thresh = g_y_start;
  const uint32_t y_end = g_y_end;
  const uint32_t ref_allele_second = g_ref_allele_second;
  uint32_t vidx_rem3 = 3;
  uint32_t vidx_rem15d3 = 5;
  uint32_t vidx_rem255d15 = 17;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  fill_ulong_zero(acc2_vec_ct * kWordsPerVec * 23, missing_acc2);
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(g_writebufs[parity][write_idx * ((uintptr_t)bgen_compressed_buf_max)]);
    uint32_t* variant_bytect_iter = &(g_variant_bytects[parity][write_idx]);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    for (; write_idx < write_idx_end; ++write_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= y_thresh) {
	if (variant_uidx < y_end) {
	  y_thresh = y_end;
	  is_y = 1;
	} else {
	  y_thresh = 0xffffffffU;
	  is_y = 0;
	}
      }
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
      }
      uint32_t dosage_ct;
      uint32_t is_explicit_alt1;
      pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
	g_error_ret = reterr;
	break;
      }
      if (ref_allele_idx + ref_allele_second == 1) {
	genovec_invert_unsafe(sample_ct, genovec);
	biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      uint32_t widx = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint16_t* bgen_geno_buf_iter = bgen_geno_buf;
      if (!dosage_ct) {
	while (1) {
	  if (widx >= sample_ctl2_m1) {
	    if (widx > sample_ctl2_m1) {
	      break;
	    }
	    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	  }
	  uintptr_t geno_word = genovec[widx];
	  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	    memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
	    bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
	    geno_word >>= 2;
	  }
	  ++widx;
	}
      } else {
	const halfword_t* dosage_present_alias = (halfword_t*)dosage_present;
	const dosage_t* dosage_vals_iter = dosage_vals;
	if (!is_explicit_alt1) {
	  while (1) {
	    if (widx >= sample_ctl2_m1) {
	      if (widx > sample_ctl2_m1) {
		break;
	      }
	      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	    }
	    uintptr_t geno_word = genovec[widx];
	    uint32_t dosage_present_hw = dosage_present_alias[widx];
	    if (!dosage_present_hw) {
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
		bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
		geno_word >>= 2;
	      }
	    } else {
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		if (dosage_present_hw & 1) {
		  uint32_t dosage_int = *dosage_vals_iter++;
		  dosage_int *= 2;
		  if (dosage_int <= kDosageMax) {
		    *bgen_geno_buf_iter++ = kDosageMax - dosage_int;
		    *bgen_geno_buf_iter++ = dosage_int;
		    *bgen_geno_buf_iter++ = 0;
		  } else {
		    dosage_int -= kDosageMax;
		    *bgen_geno_buf_iter++ = 0;
		    *bgen_geno_buf_iter++ = kDosageMax - dosage_int;
		    *bgen_geno_buf_iter++ = dosage_int;
		  }
		} else {
		  memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
		  bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
		}
		geno_word >>= 2;
		dosage_present_hw >>= 1;
	      }
	    }
	    ++widx;
	  }
	} else {
	  // todo
	  assert(0);
	}
      }
      uLongf compressed_blen = bgen_compressed_buf_max;
      if (compress(writebuf_iter, &compressed_blen, (const unsigned char*)bgen_geno_buf, bgen_geno_buf_blen)) {
	// is this actually possible?
	g_error_ret = kPglRetNomem;
	break;
      }
      *variant_bytect_iter++ = compressed_blen;
      writebuf_iter = &(writebuf_iter[bgen_compressed_buf_max]);
      if (is_y) {
	interleaved_mask_zero(sex_male_collapsed_interleaved, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, missing_acc2);
      if (!(--vidx_rem3)) {
	unroll_zero_incr_2_4(acc2_vec_ct, missing_acc2, missing_acc4);
	vidx_rem3 = 3;
	if (!(--vidx_rem15d3)) {
	  unroll_zero_incr_4_8(acc4_vec_ct, missing_acc4, missing_acc8);
	  vidx_rem15d3 = 5;
	  if (!(--vidx_rem255d15)) {
	    unroll_zero_incr_8_32(acc8_vec_ct, missing_acc8, missing_acc32);
	    vidx_rem255d15 = 17;
	  }
	}
      }
    }
    if (is_last_block) {
      unroll_incr_2_4(missing_acc2, acc2_vec_ct, missing_acc4);
      unroll_incr_4_8(missing_acc4, acc4_vec_ct, missing_acc8);
      unroll_incr_8_32(missing_acc8, acc8_vec_ct, missing_acc32);
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

pglerr_t export_bgen11(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, exportf_flags_t exportf_modifier, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_geno_cts) {
  // isomorphic to export_ox_gen().
  assert(sample_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  // use gzip instead of zstd here.
  ZWRAP_useZSTDcompression(0);
  {
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    const uint32_t max_chr_slen = get_max_chr_slen(cip);
    const uintptr_t bgen_compressed_buf_max = compressBound(6LU * sample_ct);
#ifdef __LP64__
    if (bgen_compressed_buf_max > 0xffffffffU) {
      logerrprint("Error: Too many samples for .bgen format.\n");
      goto export_bgen11_ret_INCONSISTENT_INPUT;
    }
#endif
    g_bgen_compressed_buf_max = bgen_compressed_buf_max;
    const uintptr_t writebuf_len = bgen_compressed_buf_max + 2 * max_allele_slen + 2 * kMaxIdSlen + 32;
    char* chr_buf;
    unsigned char* writebuf;
    uintptr_t* sex_male_collapsed_tmp;
    if (bigstack_alloc_c(max_chr_slen, &chr_buf) ||
        bigstack_alloc_uc(writebuf_len, &writebuf) ||
	bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
	bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_tmp)) {
      goto export_bgen11_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed_tmp);
    fill_interleaved_mask_vec(sex_male_collapsed_tmp, sample_ctv, g_sex_male_collapsed_interleaved);
    bigstack_reset(sex_male_collapsed_tmp);

    const uintptr_t max_write_block_byte_ct = bigstack_left() / 4;
    uint32_t max_write_block_size = kPglVblockSize;
    while (1) {
      // limit each write buffer to 1/4 of remaining workspace
      if (((uint64_t)(bgen_compressed_buf_max + sizeof(int32_t))) * max_write_block_size <= max_write_block_byte_ct) {
	break;
      }
      if (max_write_block_size <= kBitsPerVec) {
	goto export_bgen11_ret_NOMEM;
      }
      max_write_block_size /= 2;
    }
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // seems to saturate around this point
    if (calc_thread_ct > 15) {
      calc_thread_ct = 15;
    }
    if (bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[0])) ||
	bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[1])) ||
	bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[0])) ||
	bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[1])) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_missing_acc2) ||
	bigstack_alloc_usip(calc_thread_ct, &g_bgen_geno_bufs)) {
      goto export_bgen11_ret_NOMEM;
    }
    
    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    const uintptr_t track_missing_cacheline_ct = VECCT_TO_CLCT(acc2_vec_ct * 23);
    const uintptr_t bgen_geno_cacheline_ct = DIV_UP(6 * sample_ct, (kCacheline * k1LU));
    const uintptr_t thread_xalloc_cacheline_ct = track_missing_cacheline_ct + bgen_geno_cacheline_ct;
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, pgfip, &calc_thread_ct, &g_genovecs, dosage_is_present? (&g_dosage_presents) : nullptr, dosage_is_present? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto export_bgen11_ret_NOMEM;
    }
    if (read_block_size > max_write_block_size) {
      read_block_size = max_write_block_size;
    }
    
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_bgen11_ret_OPEN_FAIL;
    }
    // bgen 1.1 header
    // note that \xxx character constants are interpreted in octal, so \24 is
    // decimal 20, etc.
    memcpy(writebuf, "\24\0\0\0\24\0\0\0", 8);
    memcpy(&(writebuf[8]), &variant_ct, 4);
    memcpy(&(writebuf[12]), &sample_ct, 4);
    memcpy(&(writebuf[16]), "bgen\5\0\0\0", 8);
    if (fwrite_checked(writebuf, 24, outfile)) {
      goto export_bgen11_ret_WRITE_FAIL;
    }
    
    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_missing_acc2[tidx] = (uintptr_t*)bigstack_alloc_raw(track_missing_cacheline_ct * kCacheline);
      g_bgen_geno_bufs[tidx] = (uint16_t*)bigstack_alloc_raw(bgen_geno_cacheline_ct * kCacheline);
    }
    g_sample_ct = sample_ct;
    g_variant_include = variant_include;
    g_sample_include = sample_include;
    g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    g_calc_thread_ct = calc_thread_ct;
    g_refalt1_select = refalt1_select;
    get_xymt_start_and_end(cip, kChrOffsetY, &g_y_start, &g_y_end);
    g_ref_allele_second = ref_allele_second;
    g_cip = cip;
    
    // 6 bytes present at start of every bgen-1.1 variant record
    memcpy(writebuf, &sample_ct, 4);
    memcpy(&(writebuf[4]), "\0", 2);

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
    // 8. Write results for last block
    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t write_variant_uidx = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_slen = 0;
    
    uint32_t prev_block_write_ct = 0;
    uint32_t variant_idx = 0;
    uint32_t is_last_block = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    while (1) {
      uintptr_t cur_block_write_ct = 0;
      if (!is_last_block) {
	while (read_block_idx < read_block_ct_m1) {
	  cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	  if (cur_block_write_ct) {
	    break;
	  }
	  ++read_block_idx;
	}
	if (read_block_idx == read_block_ct_m1) {
	  cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	  cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	}
	if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
	  if (variant_idx) {
	    join_threads2z(calc_thread_ct, 0, threads);
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(export_bgen11_thread, calc_thread_ct, threads);
	  }
	  goto export_bgen11_ret_READ_FAIL;
	}
      }
      if (variant_idx) {
	join_threads2z(calc_thread_ct, is_last_block, threads);
	reterr = g_error_ret;
	if (reterr) {
	  if (!is_last_block) {
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(export_bgen11_thread, calc_thread_ct, threads);
	  }
	  if (reterr == kPglRetMalformedInput) {
	    logprint("\n");
	    logerrprint("Error: Malformed .pgen file.\n");
	  }
	  goto export_bgen11_ret_1;
	}
      }
      if (!is_last_block) {
	g_cur_block_write_ct = cur_block_write_ct;
	compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	  g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	}
	is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
	if (spawn_threads2z(export_bgen11_thread, calc_thread_ct, is_last_block, threads)) {
	  goto export_bgen11_ret_THREAD_CREATE_FAIL;
	}
      }
      parity = 1 - parity;
      if (variant_idx) {
	// write *previous* block results
	const unsigned char* compressed_data_iter = g_writebufs[parity];
	const uint32_t* variant_bytect_iter = g_variant_bytects[parity];
	for (uint32_t variant_bidx = 0; variant_bidx < prev_block_write_ct; ++variant_bidx, ++write_variant_uidx) {
	  next_set_unsafe_ck(variant_include, &write_variant_uidx);
	  if (write_variant_uidx >= chr_end) {
	    do {
	      ++chr_fo_idx;
	      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	    } while (write_variant_uidx >= chr_end);
	    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	    char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	    chr_slen = (uintptr_t)(chr_name_end - chr_buf);
	  }
	  const char* cur_variant_id = variant_ids[write_variant_uidx];
	  const uint32_t id_slen = strlen(cur_variant_id);
	  memcpy(&(writebuf[6]), &id_slen, 4);
	  // deliberately clobber top two bytes
	  unsigned char* writebuf_iter = (unsigned char*)memcpya(&(writebuf[8]), cur_variant_id, id_slen);
	  memcpy(writebuf_iter, &chr_slen, 4);
	  writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[2]), chr_buf, chr_slen);
	  memcpy(writebuf_iter, &(variant_bps[write_variant_uidx]), 4);
	  writebuf_iter = &(writebuf_iter[4]);
	  uintptr_t variant_allele_idx_base = write_variant_uidx * 2;
	  if (variant_allele_idxs) {
	    variant_allele_idx_base = variant_allele_idxs[write_variant_uidx];
	  }
	  char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
	  if (refalt1_select) {
	    ref_allele_idx = refalt1_select[write_variant_uidx * 2];
	    alt1_allele_idx = refalt1_select[write_variant_uidx * 2 + 1];
	  }
	  const char* first_allele;
	  const char* second_allele;
	  if (ref_allele_second) {
	    first_allele = cur_alleles[alt1_allele_idx];
	    second_allele = cur_alleles[ref_allele_idx];
	  } else {
	    first_allele = cur_alleles[ref_allele_idx];
	    second_allele = cur_alleles[alt1_allele_idx];
	  }
	  uint32_t allele_slen = strlen(first_allele);
	  memcpy(writebuf_iter, &allele_slen, 4);
	  writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), first_allele, allele_slen);
	  allele_slen = strlen(second_allele);
	  memcpy(writebuf_iter, &allele_slen, 4);
	  writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), second_allele, allele_slen);
	  const uint32_t cur_variant_bytect = *variant_bytect_iter++;
	  memcpy(writebuf_iter, &cur_variant_bytect, 4);
	  writebuf_iter = &(writebuf_iter[4]);
	  memcpy(writebuf_iter, compressed_data_iter, cur_variant_bytect);
	  writebuf_iter = &(writebuf_iter[cur_variant_bytect]);
	  compressed_data_iter = &(compressed_data_iter[bgen_compressed_buf_max]);
	  if (fwrite_checked(writebuf, writebuf_iter - writebuf, outfile)) {
	    if (variant_idx < variant_ct) {
	      join_threads2z(calc_thread_ct, is_last_block, threads);
	      if (!is_last_block) {
		g_cur_block_write_ct = 0;
		error_cleanup_threads2z(export_bgen11_thread, calc_thread_ct, threads);
	      }
	    }
	    goto export_bgen11_ret_WRITE_FAIL;
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
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (fclose_null(&outfile)) {
      goto export_bgen11_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    const uint32_t sample_ctav2 = acc2_vec_ct * kQuatersPerVec;
    const uintptr_t acc32_offset = acc2_vec_ct * (7 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx < calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctav2; ++uii) {
	scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  export_bgen11_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_bgen11_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_bgen11_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  export_bgen11_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
#ifdef __LP64__
  export_bgen11_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
#endif
  export_bgen11_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_bgen11_ret_1:
  ZWRAP_useZSTDcompression(1);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

/*
static uintptr_t** g_phasepresents = nullptr;
static uintptr_t** g_phaseinfos = nullptr;
static uintptr_t** g_dphase_presents = nullptr;

static uint32_t g_bgen_bit_precision = 0;
static uint32_t g_bgen_uncompressed_buf_max = 0;

// memcpy(target,
//   &(g_bgen_hardcall_write[cur_geno * biallelic_diploid_byte_ct]),
//   biallelic_diploid_byte_ct) should do the right thing
static unsigned char* g_bgen_hardcall_write = nullptr;

THREAD_FUNC_DECL export_bgen13_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
  const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
  const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
  uintptr_t* missing_acc2 = g_missing_acc2[tidx];
  uintptr_t* missing_acc4 = &(missing_acc2[acc2_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc8 = &(missing_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc32 = &(missing_acc8[acc8_vec_ct * kWordsPerVec]);
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  if (g_phasepresents) {
    phasepresent = g_phasepresents[tidx];
    phaseinfo = g_phaseinfos[tidx];
  }
  uintptr_t* dosage_present = g_dosage_presents? g_dosage_presents[tidx] : nullptr;
  uintptr_t* dphase_present = g_dphase_presents? g_dphase_presents[tidx] : nullptr;
  dosage_t* dosage_vals = dosage_present? g_dosage_val_bufs[tidx] : nullptr;
  unsigned char* uncompressed_bgen_geno_buf = g_thread_wkspaces[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const chr_info_t* cip = g_cip;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const unsigned char* bgen_diploid_hardcall_write = g_bgen_hardcall_write;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t sample_ctl2_m1 = QUATERCT_TO_WORDCT(sample_ct) - 1;
  const uint32_t bit_precision = g_bgen_bit_precision;
  const uint32_t bytes_per_prob = DIV_UP(bit_precision, CHAR_BIT);

  // note that this applies to both unphased and phased output, for different
  // reasons
  const uint32_t biallelic_diploid_byte_ct = 2 * bytes_per_prob;
  const unsigned char* bgen_haploid_hardcall_write = &(bgen_diploid_hardcall_write[4 * biallelic_diploid_byte_ct]);
  ;;;

  const uint32_t bgen_uncompressed_buf_max = g_bgen_uncompressed_buf_max;
  const uint32_t bgen_compressed_buf_max = g_bgen_compressed_buf_max;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  uint32_t chr_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t chr_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;

  uint32_t is_haploid_or_mt = 0; // includes chrX and chrY
  // for bgen-1.2/1.3 and VCF/BCF export, MT ploidy is 1 unless the call is
  //   heterozygous (i.e. it's treated the same way as an ordinary haploid
  //   chromosome); similarly for chrX male ploidy
  // for bgen-1.2/1.3, chrY female (but not unknown-sex) ploidy is 0 when
  //   genotype is missing

  const uint32_t ref_allele_second = g_ref_allele_second;
  uint32_t vidx_rem3 = 3;
  uint32_t vidx_rem15d3 = 5;
  uint32_t vidx_rem255d15 = 17;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  fill_ulong_zero(acc2_vec_ct * kWordsPerVec * 23, missing_acc2);
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(g_writebufs[parity][write_idx * ((uintptr_t)bgen_compressed_buf_max)]);
    uint32_t* variant_bytect_iter = &(g_variant_bytects[parity][write_idx]);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    for (; write_idx < write_idx_end; ++write_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	is_x = (chr_idx == x_code);
	is_y = (chr_idx == y_code);
	is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == mt_code);
      }
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
      }
      // todo: export phase info
      // todo: multiallelic cases
      uint32_t dosage_ct;
      uint32_t is_explicit_alt1;
      pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
	g_error_ret = reterr;
	break;
      }
      if (ref_allele_idx + ref_allele_second == 1) {
	genovec_invert_unsafe(sample_ct, genovec);
	biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      unsigned char* bgen_geno_buf_iter = uncompressed_bgen_geno_buf;
      // 4 bytes: # of samples
      // 2 bytes: # of alleles
      // 1 byte: minimum ploidy
      // 1 byte: maximum ploidy
      // sample_ct bytes: high bit = missing, low bits = ploidy
      // 1 byte: is_phased
      // 1 byte: bit_precision
      bgen_geno_buf_iter = memcpya(bgen_geno_buf_iter, &sample_ct, 4);
      uint32_t widx = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      if (!dosage_ct) {
	if (!is_haploid_or_mt) {
	  // 2 alleles, min ploidy == max ploidy == 2
	  *((uint32_t*)bgen_geno_buf_iter)++ = 0x2020002;
	  unsigned char* sample_ploidy_and_missingness = bgen_geno_buf_iter;
	  bgen_geno_buf_iter = memseta(bgen_geno_buf_iter, 2, sample_ct);
	  *bgen_geno_buf_iter++ = 0; // not phased
	  *bgen_geno_buf_iter++ = bit_precision;
	  while (1) {
	    if (widx >= sample_ctl2_m1) {
	      if (widx > sample_ctl2_m1) {
		break;
	      }
	      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	    }
	    uintptr_t geno_word = genovec[widx];
	    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
	      const uintptr_t cur_geno = geno_word & 3;
	      bgen_geno_buf_iter = memcpya(bgen_geno_buf_iter, &(bgen_hardcall_write[cur_geno * biallelic_diploid_byte_ct]), biallelic_diploid_byte_ct);
	      if (cur_geno == 3) {
		// maybe handle this in a different loop?
		sample_ploidy_and_missingness_iter[sample_idx_lowbits] = 130;
	      }
	      geno_word >>= 2;
	    }
	    ++widx;
	    sample_ploidy_and_missingness_iter = &(sample_ploidy_and_missingness_iter[kBitsPerWordD2]);
	  }
	} else if (is_x) {
	} else {
	  // ...
	}
      } else {
	const halfword_t* dosage_present_alias = (halfword_t*)dosage_present;
	const dosage_t* dosage_vals_iter = dosage_vals;
	if (!is_explicit_alt1) {
	  while (1) {
	    if (widx >= sample_ctl2_m1) {
	      if (widx > sample_ctl2_m1) {
		break;
	      }
	      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	    }
	    uintptr_t geno_word = genovec[widx];
	    uint32_t dosage_present_hw = dosage_present_alias[widx];
	    if (!dosage_present_hw) {
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
		bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
		geno_word >>= 2;
	      }
	    } else {
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		if (dosage_present_hw & 1) {
		  uint32_t dosage_int = *dosage_vals_iter++;
		  dosage_int *= 2;
		  if (dosage_int <= kDosageMax) {
		    *bgen_geno_buf_iter++ = kDosageMax - dosage_int;
		    *bgen_geno_buf_iter++ = dosage_int;
		    *bgen_geno_buf_iter++ = 0;
		  } else {
		    dosage_int -= kDosageMax;
		    *bgen_geno_buf_iter++ = 0;
		    *bgen_geno_buf_iter++ = kDosageMax - dosage_int;
		    *bgen_geno_buf_iter++ = dosage_int;
		  }
		} else {
		  memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
		  bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
		}
		geno_word >>= 2;
		dosage_present_hw >>= 1;
	      }
	    }
	    ++widx;
	  }
	} else {
	  // todo
	  assert(0);
	}
      }
      uLongf compressed_blen = bgen_compressed_buf_max;
      if (compress(writebuf_iter, &compressed_blen, (const unsigned char*)bgen_geno_buf, bgen_geno_buf_blen)) {
	// is this actually possible?
	g_error_ret = kPglRetNomem;
	break;
      }
      *variant_bytect_iter++ = compressed_blen;
      writebuf_iter = &(writebuf_iter[bgen_compressed_buf_max]);
      if (is_y) {
	interleaved_mask_zero(sex_male_collapsed_interleaved, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, missing_acc2);
      if (!(--vidx_rem3)) {
	unroll_zero_incr_2_4(acc2_vec_ct, missing_acc2, missing_acc4);
	vidx_rem3 = 3;
	if (!(--vidx_rem15d3)) {
	  unroll_zero_incr_4_8(acc4_vec_ct, missing_acc4, missing_acc8);
	  vidx_rem15d3 = 5;
	  if (!(--vidx_rem255d15)) {
	    unroll_zero_incr_8_32(acc8_vec_ct, missing_acc8, missing_acc32);
	    vidx_rem255d15 = 17;
	  }
	}
      }
    }
    if (is_last_block) {
      unroll_incr_2_4(missing_acc2, acc2_vec_ct, missing_acc4);
      unroll_incr_4_8(missing_acc4, acc4_vec_ct, missing_acc8);
      unroll_incr_8_32(missing_acc8, acc8_vec_ct, missing_acc32);
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

pglerr_t export_bgen13(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, exportf_flags_t exportf_modifier, uint32_t exportf_bits, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_geno_cts) {
  assert(sample_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  // use gzip iff v1.2
  if (exportf_modifier & kfExportfBgen12) {
    ZWRAP_useZSTDcompression(0);
  }
  {
    if (exportf_bits > 16) {
      logerrprint("Error: bits= parameter is currently limited to 16.  (This is sufficient to\ncapture all information in a .pgen file.)\n");
      reterr = kPglRetNotYetSupported;
      goto export_bgen13_ret_1;
    }
    if (pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent)) {
      logerrprint("Error: Export of phase information to .bgen files is currently under\ndevelopment.\n");
      reterr = kPglRetNotYetSupported;
      goto export_bgen13_ret_1;
    }
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    const uint32_t max_chr_slen = get_max_chr_slen(cip);
    const uintptr_t bgen_compressed_buf_max = compressBound(6LU * sample_ct);
#ifdef __LP64__
    if (bgen_compressed_buf_max > 0xffffffffU) {
      logerrprint("Error: Too many samples for .bgen format.\n");
      goto export_bgen13_ret_INCONSISTENT_INPUT;
    }
#endif
    g_bgen_compressed_buf_max = bgen_compressed_buf_max;
    const uintptr_t writebuf_len = bgen_compressed_buf_max + 2 * max_allele_slen + 2 * kMaxIdSlen + 32;
    char* chr_buf;
    unsigned char* writebuf;
    uintptr_t* sex_male_collapsed_tmp;
    if (bigstack_alloc_c(max_chr_slen, &chr_buf) ||
        bigstack_alloc_uc(writebuf_len, &writebuf) ||
	bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
	bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_tmp)) {
      goto export_bgen13_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed_tmp);
    fill_interleaved_mask_vec(sex_male_collapsed_tmp, sample_ctv, g_sex_male_collapsed_interleaved);
    bigstack_reset(sex_male_collapsed_tmp);

    const uintptr_t max_write_block_byte_ct = bigstack_left() / 4;
    uint32_t max_write_block_size = kPglVblockSize;
    while (1) {
      // limit each write buffer to 1/4 of remaining workspace
      if (((uint64_t)(bgen_compressed_buf_max + sizeof(int32_t))) * max_write_block_size <= max_write_block_byte_ct) {
	break;
      }
      if (max_write_block_size <= kBitsPerVec) {
	goto export_bgen13_ret_NOMEM;
      }
      max_write_block_size /= 2;
    }
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // seems to saturate around this point
    if (calc_thread_ct > 15) {
      calc_thread_ct = 15;
    }
    if (bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[0])) ||
	bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[1])) ||
	bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[0])) ||
	bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[1])) ||
	bigstack_alloc_ulp(calc_thread_ct, &g_missing_acc2) ||
	bigstack_alloc_usip(calc_thread_ct, &g_bgen_geno_bufs)) {
      goto export_bgen13_ret_NOMEM;
    }
    
    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    const uintptr_t track_missing_cacheline_ct = VECCT_TO_CLCT(acc2_vec_ct * 23);
    const uintptr_t bgen_geno_cacheline_ct = DIV_UP(6 * sample_ct, (kCacheline * k1LU));
    const uintptr_t thread_xalloc_cacheline_ct = track_missing_cacheline_ct + bgen_geno_cacheline_ct;
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, sample_ct, raw_variant_ct, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, pgfip, &calc_thread_ct, &g_genovecs, dosage_is_present? (&g_dosage_presents) : nullptr, dosage_is_present? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto export_bgen13_ret_NOMEM;
    }
    if (read_block_size > max_write_block_size) {
      read_block_size = max_write_block_size;
    }
    
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_bgen13_ret_OPEN_FAIL;
    }
    // bgen 1.1 header
    // note that \xxx character constants are interpreted in octal, so \24 is
    // decimal 20, etc.
    memcpy(writebuf, "\24\0\0\0\24\0\0\0", 8);
    memcpy(&(writebuf[8]), &variant_ct, 4);
    memcpy(&(writebuf[12]), &sample_ct, 4);
    memcpy(&(writebuf[16]), "bgen\5\0\0\0", 8);
    if (fwrite_checked(writebuf, 24, outfile)) {
      goto export_bgen13_ret_WRITE_FAIL;
    }
    
    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_missing_acc2[tidx] = (uintptr_t*)bigstack_alloc_raw(track_missing_cacheline_ct * kCacheline);
      g_bgen_geno_bufs[tidx] = (uint16_t*)bigstack_alloc_raw(bgen_geno_cacheline_ct * kCacheline);
    }
    g_sample_ct = sample_ct;
    g_variant_include = variant_include;
    g_sample_include = sample_include;
    g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    g_calc_thread_ct = calc_thread_ct;
    g_refalt1_select = refalt1_select;
    get_xymt_start_and_end(cip, kChrOffsetY, &g_y_start, &g_y_end);
    g_ref_allele_second = ref_allele_second;
    g_cip = cip;
    
    // 6 bytes present at start of every bgen-1.1 variant record
    memcpy(writebuf, &sample_ct, 4);
    memcpy(&(writebuf[4]), "\0", 2);

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
    // 8. Write results for last block
    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t write_variant_uidx = 0;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_slen = 0;
    
    uint32_t prev_block_write_ct = 0;
    uint32_t variant_idx = 0;
    uint32_t is_last_block = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    while (1) {
      uintptr_t cur_block_write_ct = 0;
      if (!is_last_block) {
	while (read_block_idx < read_block_ct_m1) {
	  cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
	  if (cur_block_write_ct) {
	    break;
	  }
	  ++read_block_idx;
	}
	if (read_block_idx == read_block_ct_m1) {
	  cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
	  cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
	}
	if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
	  if (variant_idx) {
	    join_threads2z(calc_thread_ct, 0, threads);
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(export_bgen13_thread, calc_thread_ct, threads);
	  }
	  goto export_bgen13_ret_READ_FAIL;
	}
      }
      if (variant_idx) {
	join_threads2z(calc_thread_ct, is_last_block, threads);
	reterr = g_error_ret;
	if (reterr) {
	  if (!is_last_block) {
	    g_cur_block_write_ct = 0;
	    error_cleanup_threads2z(export_bgen13_thread, calc_thread_ct, threads);
	  }
	  if (reterr == kPglRetMalformedInput) {
	    logprint("\n");
	    logerrprint("Error: Malformed .pgen file.\n");
	  }
	  goto export_bgen13_ret_1;
	}
      }
      if (!is_last_block) {
	g_cur_block_write_ct = cur_block_write_ct;
	compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
	for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
	  g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
	  g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
	}
	is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
	if (spawn_threads2z(export_bgen13_thread, calc_thread_ct, is_last_block, threads)) {
	  goto export_bgen13_ret_THREAD_CREATE_FAIL;
	}
      }
      parity = 1 - parity;
      if (variant_idx) {
	// write *previous* block results
	const unsigned char* compressed_data_iter = g_writebufs[parity];
	const uint32_t* variant_bytect_iter = g_variant_bytects[parity];
	for (uint32_t variant_bidx = 0; variant_bidx < prev_block_write_ct; ++variant_bidx, ++write_variant_uidx) {
	  next_set_unsafe_ck(variant_include, &write_variant_uidx);
	  if (write_variant_uidx >= chr_end) {
	    do {
	      ++chr_fo_idx;
	      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	    } while (write_variant_uidx >= chr_end);
	    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	    char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	    chr_slen = (uintptr_t)(chr_name_end - chr_buf);
	  }
	  const char* cur_variant_id = variant_ids[write_variant_uidx];
	  const uint32_t id_slen = strlen(cur_variant_id);
	  memcpy(&(writebuf[6]), &id_slen, 4);
	  // deliberately clobber top two bytes
	  unsigned char* writebuf_iter = (unsigned char*)memcpya(&(writebuf[8]), cur_variant_id, id_slen);
	  memcpy(writebuf_iter, &chr_slen, 4);
	  writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[2]), chr_buf, chr_slen);
	  memcpy(writebuf_iter, &(variant_bps[write_variant_uidx]), 4);
	  writebuf_iter = &(writebuf_iter[4]);
	  uintptr_t variant_allele_idx_base = write_variant_uidx * 2;
	  if (variant_allele_idxs) {
	    variant_allele_idx_base = variant_allele_idxs[write_variant_uidx];
	  }
	  char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
	  if (refalt1_select) {
	    ref_allele_idx = refalt1_select[write_variant_uidx * 2];
	    alt1_allele_idx = refalt1_select[write_variant_uidx * 2 + 1];
	  }
	  const char* first_allele;
	  const char* second_allele;
	  if (ref_allele_second) {
	    first_allele = cur_alleles[alt1_allele_idx];
	    second_allele = cur_alleles[ref_allele_idx];
	  } else {
	    first_allele = cur_alleles[ref_allele_idx];
	    second_allele = cur_alleles[alt1_allele_idx];
	  }
	  uint32_t allele_slen = strlen(first_allele);
	  memcpy(writebuf_iter, &allele_slen, 4);
	  writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), first_allele, allele_slen);
	  allele_slen = strlen(second_allele);
	  memcpy(writebuf_iter, &allele_slen, 4);
	  writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), second_allele, allele_slen);
	  const uint32_t cur_variant_bytect = *variant_bytect_iter++;
	  memcpy(writebuf_iter, &cur_variant_bytect, 4);
	  writebuf_iter = &(writebuf_iter[4]);
	  memcpy(writebuf_iter, compressed_data_iter, cur_variant_bytect);
	  writebuf_iter = &(writebuf_iter[cur_variant_bytect]);
	  compressed_data_iter = &(compressed_data_iter[bgen_compressed_buf_max]);
	  if (fwrite_checked(writebuf, writebuf_iter - writebuf, outfile)) {
	    if (variant_idx < variant_ct) {
	      join_threads2z(calc_thread_ct, is_last_block, threads);
	      if (!is_last_block) {
		g_cur_block_write_ct = 0;
		error_cleanup_threads2z(export_bgen13_thread, calc_thread_ct, threads);
	      }
	    }
	    goto export_bgen13_ret_WRITE_FAIL;
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
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (fclose_null(&outfile)) {
      goto export_bgen13_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    const uint32_t sample_ctav2 = acc2_vec_ct * kQuatersPerVec;
    const uintptr_t acc32_offset = acc2_vec_ct * (7 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx < calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctav2; ++uii) {
	scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  export_bgen13_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_bgen13_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_bgen13_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  export_bgen13_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
#ifdef __LP64__
  export_bgen13_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
#endif
  export_bgen13_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_bgen13_ret_1:
  ZWRAP_useZSTDcompression(1);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}
*/

pglerr_t export_ox_sample(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const uint32_t* sample_missing_geno_cts, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uint32_t y_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t pheno_ctl = BITCT_TO_WORDCT(pheno_ct);
    char* writebuf;
    uintptr_t* is_basic_categorical;
    // if phenotype is categorical, and all (non-null) category names are of
    // the form P[positive integer], then it's best to emit the positive
    // integer in the name string instead of the internal index.
    if (bigstack_calloc_ul(pheno_ctl, &is_basic_categorical) ||
	bigstack_alloc_c(kMaxMediumLine + max_sample_id_blen + 32 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf)) {
      goto export_ox_sample_ret_NOMEM;
    }
    
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ox_sample_ret_OPEN_FAIL;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya(writebuf, "ID_1 ID_2 missing sex");
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      const pheno_col_t* cur_pheno_col = &(pheno_cols[pheno_idx]);
      if (cur_pheno_col->type_code == kPhenoDtypeCat) {
	const uint32_t nn_cat_ct = cur_pheno_col->nonnull_category_ct;
	char** cur_cat_names = cur_pheno_col->category_names;
	uint32_t cat_idx;
	for (cat_idx = 1; cat_idx <= nn_cat_ct; ++cat_idx) {
	  const char* cat_name_iter = cur_cat_names[cat_idx];
	  if (*cat_name_iter == 'C') {
	    uint32_t char_code = *(++cat_name_iter);
	    if ((char_code - 49) < 9) {
	      uint32_t uii;
	      if (!scan_posint_capped(cat_name_iter, 0x7fffffff, &uii)) {
		continue;
	      }
	    }
	  }
	  break;
	}
	if (cat_idx == nn_cat_ct + 1) {
	  set_bit(pheno_idx, is_basic_categorical);
	}
      }
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto export_ox_sample_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
    }
    append_binary_eoln(&write_iter);

    write_iter = strcpya(write_iter, "0 0 0 D");
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      const pheno_dtype_t cur_type_code = pheno_cols[pheno_idx].type_code;
      if (cur_type_code == kPhenoDtypeCc) {
	*write_iter++ = 'B';
      } else if (cur_type_code == kPhenoDtypeQt) {
	// .psam file does not distinguish between "continuous covariate" and
	// "continuous phenotype", that's lost on round-trip
	*write_iter++ = 'P';
      } else {
	*write_iter++ = 'D';
      }
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto export_ox_sample_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
    }
    append_binary_eoln(&write_iter);

    const double nonmale_geno_ct_recip = 1.0 / ((double)((int32_t)(variant_ct - y_ct)));
    const double male_geno_ct_recip = 1.0 / ((double)((int32_t)variant_ct));
    uintptr_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_ul_unsafe_ck(sample_include, &sample_uidx);
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      const char* fid_end = (const char*)rawmemchr(cur_sample_id, '\t');
      write_iter = memcpyax(write_iter, cur_sample_id, (uintptr_t)(fid_end - cur_sample_id), ' ');
      write_iter = strcpya(write_iter, &(fid_end[1]));
      *write_iter++ = ' ';
      const int32_t cur_missing_geno_ct = sample_missing_geno_cts[sample_idx];
      if (is_set(sex_male, sample_uidx)) {
        write_iter = dtoa_g(cur_missing_geno_ct * male_geno_ct_recip, write_iter);
	write_iter = strcpya(write_iter, " 1");
      } else {
	write_iter = dtoa_g(cur_missing_geno_ct * nonmale_geno_ct_recip, write_iter);
	*write_iter++ = ' ';
	if (is_set(sex_nm, sample_uidx)) {
	  *write_iter++ = '2';
	} else {
	  write_iter = strcpya(write_iter, "NA");
	}
      }
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
	*write_iter++ = ' ';
        const pheno_col_t* cur_pheno_col = &(pheno_cols[pheno_idx]);
	if (!is_set(cur_pheno_col->nonmiss, sample_uidx)) {
	  write_iter = strcpya(write_iter, "NA");
	} else {
	  const pheno_dtype_t cur_type_code = cur_pheno_col->type_code;
	  if (cur_type_code == kPhenoDtypeCc) {
	    *write_iter++ = '0' + is_set(cur_pheno_col->data.cc, sample_uidx);
	  } else if (cur_type_code == kPhenoDtypeQt) {
	    write_iter = dtoa_g(cur_pheno_col->data.qt[sample_uidx], write_iter);
	  } else {
	    const uint32_t cur_cat_idx = cur_pheno_col->data.cat[sample_uidx];
	    if (is_set(is_basic_categorical, pheno_idx)) {
	      write_iter = strcpya(write_iter, &(cur_pheno_col->category_names[cur_cat_idx][1]));
	    } else {
	      write_iter = uint32toa(cur_cat_idx, write_iter);
	    }
	  }
	}
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= writebuf_flush) {
	if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), outfile)) {
	  goto export_ox_sample_ret_WRITE_FAIL;
	}
	write_iter = writebuf;
      }
    }
    if (write_iter != writebuf) {
      if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	goto export_ox_sample_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto export_ox_sample_ret_WRITE_FAIL;
    }
  }
  while (0) {
  export_ox_sample_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ox_sample_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ox_sample_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
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

char* diploid_vcf_dosage_print(uint32_t dosage_int, uint32_t write_ds, char* write_iter) {
  if (write_ds) {
    return print_small_dosage(dosage_int, write_iter);
  }
  if (dosage_int <= kDosageMid) {
    write_iter = print_small_dosage(kDosageMid - dosage_int, write_iter);
    *write_iter++ = ',';
    write_iter = print_small_dosage(dosage_int, write_iter);
    return strcpya(write_iter, ",0");
  }
  write_iter = strcpya(write_iter, "0,");
  write_iter = print_small_dosage(kDosageMax - dosage_int, write_iter);
  *write_iter++ = ',';
  return print_small_dosage(dosage_int - kDosageMid, write_iter);
}

// assumes rawval is in [0, 327679]
static_assert(kDosageMax == 32768, "haploid_dosage_print() needs to be updated.");
char* haploid_dosage_print(uint32_t rawval, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/32768, (n + 0.5)/32768).
  *start++ = '0' + (rawval / 32768);
  rawval = rawval % 32768;
  if (!rawval) {
    // this shouldn't come up for now?
    return start;
  }
  *start++ = '.';

  // (rawval * 2) is in 65536ths
  // 65536 * 625 = 40960k
  const uint32_t range_top_40960k = rawval * 1250 + 625;
  // ok to check half-open interval since we never hit boundary
  if ((range_top_40960k % 4096) < 1250) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_40960k / 4096;
    return uitoa_trunc4(four_decimal_places, start);
  }
  
  // we wish to print (100000 * remainder + 16384) / 32768, left-0-padded.  and
  // may as well banker's round too.
  //
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4.  32768/64 = 512.
  const uint32_t five_decimal_places = ((3125 * rawval + 512) / 1024) - ((rawval % 2048) == 512);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return uitoa_trunc4(last_four_digits, start);
  }
  return start;
}

interr_t flexbwrite_flush(char* buf, size_t len, FILE* outfile, BGZF* bgz_outfile) {
  if (outfile) {
    return fwrite_checked(buf, len, outfile);
  }
  return (bgzf_write(bgz_outfile, buf, len) < 0);
}


// these assume buf_flush - buf = kMaxMediumLine
// outfile should be nullptr iff we're doing bgzf compression
interr_t flexbwrite_flush2(char* buf_flush, FILE* outfile, BGZF* bgz_outfile, char** write_iter_ptr) {
  char* buf = &(buf_flush[-((int32_t)kMaxMediumLine)]);
  char* buf_end = *write_iter_ptr;
  *write_iter_ptr = buf;
  return flexbwrite_flush(buf, (uintptr_t)(buf_end - buf), outfile, bgz_outfile);
}

static inline interr_t flexbwrite_ck(char* buf_flush, FILE* outfile, BGZF* bgz_outfile, char** write_iter_ptr) {
  if ((*write_iter_ptr) < buf_flush) {
    return 0;
  }
  return flexbwrite_flush2(buf_flush, outfile, bgz_outfile, write_iter_ptr);
}


#ifdef __arm__
  #error "Unaligned accesses in export_vcf()."
#endif
pglerr_t export_vcf(char* xheader, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const char* sample_ids, const char* sids, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, __maybe_unused uint32_t max_thread_ct, exportf_flags_t exportf_modifier, idpaste_t exportf_id_paste, char exportf_id_delim, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  gzFile gz_pvar_reload = nullptr;
  BGZF* bgz_outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!(exportf_modifier & kfExportfBgz)) {
      strcpy(outname_end, ".vcf");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto export_vcf_ret_OPEN_FAIL;
      }
    } else {
      strcpy(outname_end, ".vcf.gz");
      bgz_outfile = bgzf_open(outname, "w");
      if (!bgz_outfile) {
	goto export_vcf_ret_OPEN_FAIL;
      }
#ifndef _WIN32
      if (max_thread_ct > 1) {
	// 128 doesn't seem any worse than 256 (and is clearly better than 64)
	// also tried reducing thread count by 1, that seems worse
	if (bgzf_mt2(g_bigstack_end, MINV(128, max_thread_ct), 128, &g_bigstack_base, bgz_outfile)) {
	  goto export_vcf_ret_NOMEM;
	}
      }
#endif
    }
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    // CHROM, POS, ID, REF, one ALT, eoln
    uintptr_t writebuf_blen = kMaxIdSlen + 32 + max_chr_blen + 2 * max_allele_slen;
    // QUAL, FILTER, INFO, FORMAT, genotypes, eoln
    // needs to be larger for >9 alt alleles
    uint32_t write_ds = (exportf_modifier / kfExportfVcfDosageDs) & 1;
    uint32_t write_gp_or_ds = write_ds || (exportf_modifier & kfExportfVcfDosageGp);
    if (write_gp_or_ds && (!(pgfip->gflags & kfPgenGlobalDosagePresent))) {
      write_gp_or_ds = 0;
      LOGERRPRINTF("Warning: No dosage data present.  %s field will not be exported.\n", write_ds? "DS" : "GP");
      write_ds = 0;
    }
    if (writebuf_blen < ((4 * k1LU) + write_gp_or_ds * 24 - write_ds * 16) * sample_ct + 32 + max_filter_slen + info_reload_slen) {
      writebuf_blen = ((4 * k1LU) + write_gp_or_ds * 24 - write_ds * 16) * sample_ct + 32 + max_filter_slen + info_reload_slen;
    }
    writebuf_blen += kCompressStreamBlock;
    char* writebuf;
    if (bigstack_alloc_c(writebuf_blen, &writebuf)) {
      goto export_vcf_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya(writebuf, "##fileformat=VCFv4.3" EOLN_STR "##fileDate=");
    time_t rawtime;
    time(&rawtime);
    struct tm* loctime;
    loctime = localtime(&rawtime);
    write_iter += strftime(write_iter, kMaxMediumLine, "%Y%m%d", loctime);
    write_iter = strcpya(write_iter, EOLN_STR "##source=PLINKv2.00" EOLN_STR);
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    if (flexbwrite_flush(writebuf, write_iter - writebuf, outfile, bgz_outfile)) {
      goto export_vcf_ret_WRITE_FAIL;
    }
    const uint32_t chr_ctl = BITCT_TO_WORDCT(cip->chr_ct);
    uintptr_t* written_contig_header_lines;
    if (bigstack_calloc_ul(chr_ctl, &written_contig_header_lines)) {
      goto export_vcf_ret_NOMEM;
    }
    if (xheader) {
      memcpy(writebuf, "##contig=<ID=", 13);
      char* xheader_iter = xheader;
      char* xheader_end = &(xheader[xheader_blen]);
      char* line_end = xheader;
      while (line_end != xheader_end) {
	xheader_iter = line_end;
	line_end = (char*)rawmemchr(xheader_iter, '\n');
	++line_end;
	const uint32_t slen = (uintptr_t)(line_end - xheader_iter);
	if ((slen > 14) && (!memcmp(xheader_iter, "##contig=<ID=", 13))) {
	  char* contig_name_start = &(xheader_iter[13]);
	  char* contig_name_end = (char*)memchr(contig_name_start, ',', slen - 14);
	  if (!contig_name_end) {
	    // if this line is technically well-formed (ends in '>'), it's
	    // useless anyway, throw it out
	    continue;
	  }
	  const int32_t chr_idx = get_chr_code_counted(cip, contig_name_end - contig_name_start, contig_name_start);
	  if (chr_idx < 0) {
	    continue;
	  }
	  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)chr_idx];
	  if (IS_SET(written_contig_header_lines, chr_fo_idx)) {
	    logerrprint("Error: Duplicate ##contig line in .pvar file.\n");
	    goto export_vcf_ret_MALFORMED_INPUT;
	  }
	  SET_BIT(chr_fo_idx, written_contig_header_lines);
	  // if --output-chr was used at some point, we need to sync the
	  // ##contig chromosome code with the code in the VCF body.
	  write_iter = chr_name_write(cip, chr_idx, &(writebuf[13]));
	  if (flexbwrite_flush(writebuf, write_iter - writebuf, outfile, bgz_outfile)) {
	    goto export_vcf_ret_WRITE_FAIL;
	  }
	  if (flexbwrite_flush(contig_name_end, (uintptr_t)(line_end - contig_name_end), outfile, bgz_outfile)) {
	    goto export_vcf_ret_WRITE_FAIL;
	  }
	} else {
	  if (flexbwrite_flush(xheader_iter, slen, outfile, bgz_outfile)) {
	    goto export_vcf_ret_WRITE_FAIL;
	  }
	}
      }
    }
    write_iter = writebuf;
    // fill in the missing ##contig lines
    uint32_t contig_zero_written = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx < cip->chr_ct; ++chr_fo_idx) {
      if (IS_SET(written_contig_header_lines, chr_fo_idx)) {
	continue;
      }
      const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if ((!IS_SET(cip->chr_mask, chr_idx)) || are_all_bits_zero(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1])) {
	continue;
      }
      char* chr_name_write_start = strcpya(write_iter, "##contig=<ID=");
      char* chr_name_write_end = chr_name_write(cip, chr_idx, chr_name_write_start);
      if ((*chr_name_write_start == '0') && (chr_name_write_end == &(chr_name_write_start[1]))) {
	// --allow-extra-chr 0 special case
	if (contig_zero_written) {
	  continue;
	}
	contig_zero_written = 1;
	write_iter = strcpya(chr_name_write_end, ",length=2147483645");
      } else {
	if (memchr(chr_name_write_start, ':', chr_name_write_end - chr_name_write_start)) {
	  logerrprint("Error: VCF chromosome codes may not include the ':' character.\n");
	  goto export_vcf_ret_MALFORMED_INPUT;
	}
	write_iter = strcpya(chr_name_write_end, ",length=");
	if (1) {
	  write_iter = uint32toa(variant_bps[cip->chr_fo_vidx_start[chr_fo_idx + 1] - 1] + 1, write_iter);
	} else {
	  // todo: unsorted map case
	}
      }
      *write_iter++ = '>';
      append_binary_eoln(&write_iter);
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
	goto export_vcf_ret_WRITE_FAIL;
      }
    }
    bigstack_reset(written_contig_header_lines);
    const uint32_t all_nonref = pgfip->gflags & kfPgenGlobalAllNonref;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uint32_t write_pr = all_nonref;
    if (nonref_flags) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
	if (variant_include[widx] & nonref_flags[widx]) {
	  write_pr = 1;
	  break;
	}
      }
    }
    if (write_pr && (!xheader_info_pr)) {
      write_iter = strcpya(write_iter, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
    }
    if (write_ds) {
      write_iter = strcpya(write_iter, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">" EOLN_STR);
    } else if (write_gp_or_ds) {
      write_iter = strcpya(write_iter, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Phred-scaled Genotype Likelihoods\">" EOLN_STR);
    }
    // possible todo: optionally export .psam information as
    // PEDIGREE/META/SAMPLE lines in header, and make --vcf be able to read it
    write_iter = strcpya(write_iter, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" EOLN_STR "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    uint32_t write_sid = 0;
    // possible for both MAYBESID and SID to be set
    if (exportf_id_paste & kfIdpasteSid) {
      write_sid = 1;
      if (!sids) {
	max_sid_blen = 2;
      }
    } else if ((exportf_id_paste & kfIdpasteMaybesid) && sids) {
      // no nonzero check in load_psam(), so we have to do it here
      uint32_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
	next_set_unsafe_ck(sample_include, &sample_uidx);
	if (memcmp(&(sids[sample_uidx * max_sid_blen]), "0", 2)) {
	  write_sid = 1;
	  break;
	}
      }
    }
    uint32_t sample_uidx = 0;
    uint32_t id_delim_warning = 0;
    char id_delim = exportf_id_delim? exportf_id_delim : '_';
    const uintptr_t max_exported_sample_id_blen = max_sample_id_blen + write_sid * max_sid_blen;
    char* exported_sample_ids;
    const uint32_t exported_id_htable_size = get_htable_min_size(sample_ct);
    uint32_t* exported_id_htable;
    // check for duplicates
    if (bigstack_alloc_c(sample_ct * max_exported_sample_id_blen, &exported_sample_ids) ||
	bigstack_alloc_ui(exported_id_htable_size, &exported_id_htable)) {
      goto export_vcf_ret_NOMEM;
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      const char* orig_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* orig_fid_end = (const char*)rawmemchr(orig_sample_id, '\t');
      char* exported_sample_ids_iter = &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]);
      if (exportf_id_paste & kfIdpasteFid) {
	const uint32_t fid_slen = (uintptr_t)(orig_fid_end - orig_sample_id);
	if ((!id_delim_warning) && memchr(orig_sample_id, id_delim, fid_slen)) {
	  id_delim_warning = 1;
	}
	exported_sample_ids_iter = memcpyax(exported_sample_ids_iter, orig_sample_id, fid_slen, id_delim);
      }
      if (exportf_id_paste & kfIdpasteIid) {
	const char* orig_iid = &(orig_fid_end[1]);
        const uint32_t iid_slen = strlen(orig_iid);
	if ((!id_delim_warning) && memchr(orig_iid, id_delim, iid_slen)) {
	  id_delim_warning = 1;
	}
	exported_sample_ids_iter = memcpyax(exported_sample_ids_iter, orig_iid, iid_slen, id_delim);
      }
      if (write_sid) {
	if (sids) {
	  const char* orig_sid = &(sids[sample_uidx * max_sid_blen]);
	  const uint32_t sid_slen = strlen(orig_sid);
	  if ((!id_delim_warning) && memchr(orig_sid, id_delim, sid_slen)) {
	    id_delim_warning = 1;
	  }
	  exported_sample_ids_iter = memcpya(exported_sample_ids_iter, orig_sid, sid_slen);
	} else {
	  *exported_sample_ids_iter++ = '0';
	}
	++exported_sample_ids_iter;
      }
      exported_sample_ids_iter[-1] = '\0';
    }
    if (id_delim_warning) {
      if (exportf_id_delim) {
	LOGERRPRINTF("Warning: '%c' present in original sample IDs; --vcf will not be able to\nreconstruct them.  Consider rerunning with a different --export id-delim=\nvalue.\n", exportf_id_delim);
      } else {
	logerrprint("Warning: '_' present in original sample IDs; --vcf will not be able to\nreconstruct them.  Consider rerunning with a suitable --export id-delim= value.\n");
      }
    }
    if (populate_strbox_htable(exported_sample_ids, sample_ct, max_exported_sample_id_blen, exported_id_htable_size, exported_id_htable)) {
      logerrprint("Warning: Duplicate sample ID(s) present in exported VCF file.\n");
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]));
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
	goto export_vcf_ret_WRITE_FAIL;
      }
    }
    append_binary_eoln(&write_iter);
    bigstack_reset(exported_sample_ids);

    LOGPRINTFWW5("--export vcf%s to %s ... ", bgz_outfile? " bgz" : "", outname);
    fputs("0%", stdout);
    fflush(stdout);

    // includes trailing tab
    char* chr_buf;

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
	bigstack_alloc_ul(sample_ctl2, &genovec) ||
	bigstack_alloc_ul(BITCT_TO_WORDCT(kPglMaxAltAlleleCt), &allele_include)) {
      goto export_vcf_ret_NOMEM;
    }
    // For now, if phased data is present, each homozygous call is represented
    // as phased iff the previous heterozygous call was phased.  (If no
    // previous heterozygous call exists, it's treated as phased.)  This does
    // the right thing when the entire genome is phased, and it induces about
    // as good a phase set approximation as you can get without explicitly
    // saving that info.  But that approximation is still pretty inaccurate; as
    // soon as we have any use for them, explicit phase set support should be
    // added to pgenlib.
    const uint32_t some_phased = (pgfip->gflags / kfPgenGlobalHardcallPhasePresent) & 1;
    uintptr_t* prev_phased = nullptr;
    uintptr_t* phasepresent = nullptr;
    uintptr_t* phaseinfo = nullptr;
    if (some_phased) {
      if (bigstack_alloc_ul(sample_ctl, &prev_phased) ||
	  bigstack_alloc_ul(sample_ctl, &phasepresent) ||
	  bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
	goto export_vcf_ret_NOMEM;
      }
      fill_all_bits(sample_ct, prev_phased);
    }

    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (write_gp_or_ds) {
      if (bigstack_alloc_ul(sample_ctl, &dosage_present) ||
	  bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
	goto export_vcf_ret_NOMEM;
      }
    }

    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t info_col_idx = 0;
    if (pvar_info_reload) {
      reterr = pvar_info_reload_header(pvar_info_reload, &gz_pvar_reload, &loadbuf, &loadbuf_size, &info_col_idx);
      if (reterr) {
	goto export_vcf_ret_1;
      }
    }

    // assumes little-endian
    uint32_t basic_genotext[4];
    basic_genotext[0] = 0x302f3009; // \t0/0
    basic_genotext[1] = 0x312f3009; // \t0/1
    basic_genotext[2] = 0x312f3109; // \t1/1
    basic_genotext[3] = 0x2e2f2e09; // \t./.
    char haploid_genotext[4][4];
    uint32_t haploid_genotext_blen[8]; // 4..7 = male chrX
    memcpy(haploid_genotext[0], "\t0/0", 4);
    memcpy(haploid_genotext[1], "\t0/1", 4);
    memcpy(haploid_genotext[2], "\t1/1", 4);
    memcpy(haploid_genotext[3], "\t./.", 4);
    haploid_genotext_blen[1] = 4;
    haploid_genotext_blen[4] = 2;
    haploid_genotext_blen[5] = 4;
    haploid_genotext_blen[6] = 2;
    haploid_genotext_blen[7] = 2;
    // don't bother exporting GP/DS for hardcalls
    const char* dot_ptr = &(g_one_char_strs[92]);
    const char* input_missing_geno_ptr = g_input_missing_geno_ptr;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    uint32_t chr_fo_idx = 0xffffffffU;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t variant_uidx = 0;
    uint32_t is_x = 0;
    uint32_t is_haploid_or_mt = 0; // includes chrX and chrY
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t gz_variant_uidx = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      // a lot of this is redundant with write_pvar(), may want to factor the
      // commonalities out
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
	do {
	  ++chr_fo_idx;
	  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
	} while (variant_uidx >= chr_end);
	int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
	is_x = (chr_idx == cip->xymt_codes[kChrOffsetX]);
	is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == cip->xymt_codes[kChrOffsetMT]);
	// forced --merge-par, with diploid male output (is_x NOT set, but
	// chromosome code is X/chrX)
	if ((chr_idx == cip->xymt_codes[kChrOffsetPAR1]) || (chr_idx == cip->xymt_codes[kChrOffsetPAR2])) {
	  chr_idx = cip->xymt_codes[kChrOffsetX];
	}
	char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
	*chr_name_end = '\t';
	chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
	if (is_haploid_or_mt) {
	  if (is_x) {
	    haploid_genotext_blen[0] = 4;
	    haploid_genotext_blen[2] = 4;
	    haploid_genotext_blen[3] = 4;
	  } else {
	    haploid_genotext_blen[0] = 2;
	    haploid_genotext_blen[2] = 2;
	    haploid_genotext_blen[3] = 2;
	  }
	}
      }
      // #CHROM
      write_iter = memcpya(write_iter, chr_buf, chr_buf_blen);

      // POS
      write_iter = uint32toa_x(variant_bps[variant_uidx], '\t', write_iter);

      // ID
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], '\t');

      // REF, ALT
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
	variant_allele_idx_base = variant_allele_idxs[variant_uidx];
	cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      if (refalt1_select) {
	ref_allele_idx = refalt1_select[variant_uidx * 2];
	alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
	// this logic only works in the biallelic case
	assert(cur_allele_ct == 2);
	if (!is_haploid_or_mt) {
	  if (alt1_allele_idx) {
	    basic_genotext[0] = 0x302f3009;
	    basic_genotext[2] = 0x312f3109;
	  } else {
	    basic_genotext[0] = 0x312f3109;
	    basic_genotext[2] = 0x302f3009;
	  }
	} else {
	  if (alt1_allele_idx) {
	    memcpy(haploid_genotext[0], "\t0/0", 4);
	    memcpy(haploid_genotext[2], "\t1/1", 4);
	  } else {
	    memcpy(haploid_genotext[0], "\t1/1", 4);
	    memcpy(haploid_genotext[2], "\t0/0", 4);
	  }
	}
      }
      if ((cur_alleles[ref_allele_idx] != dot_ptr) && (cur_alleles[ref_allele_idx] != input_missing_geno_ptr)) {
        write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
      } else {
	*write_iter++ = 'N';
      }
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
	goto export_vcf_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
	fill_all_bits(cur_allele_ct, allele_include);
	CLEAR_BIT(ref_allele_idx, allele_include);
	CLEAR_BIT(alt1_allele_idx, allele_include);
        uint32_t cur_allele_uidx = 0;
	uint32_t alt_allele_idx = 2;
	do {
	  *write_iter++ = ',';
	  next_set_unsafe_ck(allele_include, &cur_allele_uidx);
	  write_iter = strcpya(write_iter, cur_alleles[cur_allele_uidx++]);
	  if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
	    goto export_vcf_ret_WRITE_FAIL;
	  }
	} while (++alt_allele_idx < cur_allele_ct);
      }

      // QUAL
      *write_iter++ = '\t';
      if ((!pvar_qual_present) || (!IS_SET(pvar_qual_present, variant_uidx))) {
	*write_iter++ = '.';
      } else {
	write_iter = ftoa_g(pvar_quals[variant_uidx], write_iter);
      }

      // FILTER
      *write_iter++ = '\t';
      if ((!pvar_filter_present) || (!IS_SET(pvar_filter_present, variant_uidx))) {
	*write_iter++ = '.';
      } else if (!IS_SET(pvar_filter_npass, variant_uidx)) {
	write_iter = strcpya(write_iter, "PASS");
      } else {
	write_iter = strcpya(write_iter, pvar_filter_storage[variant_uidx]);
      }

      // INFO
      *write_iter++ = '\t';
      const uint32_t is_pr = all_nonref || (nonref_flags && IS_SET(nonref_flags, variant_uidx));
      if (gz_pvar_reload) {
	reterr = pvar_info_reload_and_write(loadbuf_size, xheader_info_pr, info_col_idx, variant_uidx, is_pr, gz_pvar_reload, &write_iter, &gz_variant_uidx, loadbuf);
	if (reterr) {
	  goto export_vcf_ret_1;
	}
      } else {
	if (is_pr) {
	  write_iter = strcpya(write_iter, "PR");
	} else {
	  *write_iter++ = '.';
	}
      }

      // FORMAT
      write_iter = memcpyl3a(write_iter, "\tGT");
      
      uint32_t dosage_ct = 0;
      uint32_t is_explicit_alt1 = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      if (!some_phased) {
	// biallelic, nothing phased in entire file
	if (!write_gp_or_ds) {
	  reterr = pgr_read_refalt1_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec);
	} else {
	  reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
	}
	if (reterr) {
	  goto export_vcf_ret_PGR_FAIL;
	}
	if (!dosage_ct) {
	  if (!is_haploid_or_mt) {
	    // always 4 bytes wide, exploit that
	    uint32_t* write_iter_ui_alias = (uint32_t*)write_iter;
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		*write_iter_ui_alias++ = basic_genotext[genovec_word & 3];
		genovec_word >>= 2;
	      }
	      ++widx;
	    }
	    write_iter = (char*)write_iter_ui_alias;
	  } else {
	    // chrX: male homozygous/missing calls use only one character + tab
	    // other haploid/MT: this is true for nonmales too
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t sex_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uint32_t cur_geno = genovec_word & 3;
		const uint32_t cur_is_male = sex_male_hw & 1;
		write_iter = memcpya(write_iter, haploid_genotext[cur_geno], haploid_genotext_blen[cur_geno + cur_is_male * 4]);
		genovec_word >>= 2;
		sex_male_hw >>= 1;
	      }
	      ++widx;
	    }
	  }
	} else {
	  // some dosages present
	  if (write_ds) {
	    write_iter = memcpyl3a(write_iter, ":DS");
	  } else {
	    write_iter = memcpyl3a(write_iter, ":GP");
	  }
	  if (!alt1_allele_idx) {
	    biallelic_dosage16_invert(dosage_ct, dosage_vals);
	  }
	  dosage_t* dosage_vals_iter = dosage_vals;
          if (!is_haploid_or_mt) {
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uint32_t cur_geno = genovec_word & 3;
		write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
		if (dosage_present_hw & 1) {
		  *write_iter++ = ':';
		  const uint32_t dosage_int = *dosage_vals_iter++;
		  write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
		}
		genovec_word >>= 2;
		dosage_present_hw >>= 1;
	      }
	      ++widx;
	    }
	  } else {
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t sex_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
	      uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uint32_t cur_geno = genovec_word & 3;
		const uint32_t cur_is_male = sex_male_hw & 1;
		const uint32_t cur_genotext_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
		write_iter = memcpya(write_iter, haploid_genotext[cur_geno], cur_genotext_blen);
		if (dosage_present_hw & 1) {
		  *write_iter++ = ':';
		  uint32_t dosage_int = *dosage_vals_iter++;
		  if (cur_genotext_blen == 2) {
		    if (write_ds) {
		      write_iter = haploid_dosage_print(dosage_int, write_iter);
		    } else {
		      write_iter = haploid_dosage_print(kDosageMax - dosage_int, write_iter);
		      *write_iter++ = ',';
		      write_iter = haploid_dosage_print(dosage_int, write_iter);
		    }
		  } else {
		    // het haploid, or female X
		    write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
		  }
		}
		genovec_word >>= 2;
		sex_male_hw >>= 1;
		dosage_present_hw >>= 1;
	      }
	      ++widx;
	    }
	  }
	}
      } else {
	// biallelic, phased
	uint32_t at_least_one_phase_present;
	if (!write_gp_or_ds) {
	  reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, phasepresent, phaseinfo, &at_least_one_phase_present);
	} else {
	  reterr = pgr_read_refalt1_genovec_hphase_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, phasepresent, phaseinfo, &at_least_one_phase_present, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
	}
	if (reterr) {
	  goto export_vcf_ret_PGR_FAIL;
	}
	at_least_one_phase_present = (at_least_one_phase_present != 0);
	if (!dosage_ct) {
	  if (!is_haploid_or_mt) {
	    uint32_t* write_iter_ui_alias = (uint32_t*)write_iter;
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

	      // zero this out if phasepresent_ct == 0
	      const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

	      const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uintptr_t cur_geno = genovec_word & 3;

		// usually "\t0/0", etc.
		uint32_t cur_basic_genotext = basic_genotext[cur_geno];
		if (cur_geno == 1) {
		  const uint32_t cur_shift = (1U << sample_idx_lowbits);
		  if (phasepresent_halfword & cur_shift) {
		    prev_phased_halfword |= cur_shift;
		    if (phaseinfo_halfword & cur_shift) {
		      cur_basic_genotext ^= 0x1000100; // 0|1 -> 1|0
		    }
		  } else {
		    prev_phased_halfword &= ~cur_shift;
		  }
		}
		// '/' = ascii 47, '|' = ascii 124
		*write_iter_ui_alias++ = cur_basic_genotext + 0x4d0000 * ((prev_phased_halfword >> sample_idx_lowbits) & 1);
		genovec_word >>= 2;
	      }
	      ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
	      ++widx;
	    }
	    write_iter = (char*)write_iter_ui_alias;
	  } else {
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t is_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
	      uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

	      // zero this out if phasepresent_ct == 0
	      const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

	      const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uint32_t cur_geno = genovec_word & 3;
		const uint32_t cur_is_male = is_male_hw & 1;
		const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
		write_iter = memcpya(write_iter, haploid_genotext[cur_geno], cur_blen);
		if (cur_blen == 4) {
		  if (cur_geno == 1) {
		    // a bit redundant with how is_male_hw is handled, but
		    // updating this on every loop iteration doesn't seem better
		    const uint32_t cur_shift = (1U << sample_idx_lowbits);
		    if (phasepresent_halfword & cur_shift) {
		      prev_phased_halfword |= cur_shift;
		      if (phaseinfo_halfword & cur_shift) {
			memcpy(&(write_iter[-4]), "\t1|0", 4);
		      } else {
			write_iter[-2] = '|';
		      }
		    } else {
		      prev_phased_halfword &= ~cur_shift;
		    }
		  } else if ((prev_phased_halfword >> sample_idx_lowbits) & 1) {
		    write_iter[-2] = '|';
		  }
		}
		genovec_word >>= 2;
		is_male_hw >>= 1;
	      }
	      ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
	      ++widx;
	    }
	  }
	} else {
	  // both dosage and phase present
	  if (write_ds) {
	    write_iter = memcpyl3a(write_iter, ":DS");
	  } else {
	    write_iter = memcpyl3a(write_iter, ":GP");
	  }
	  if (!alt1_allele_idx) {
	    biallelic_dosage16_invert(dosage_ct, dosage_vals);
	  }
	  dosage_t* dosage_vals_iter = dosage_vals;
	  if (!is_haploid_or_mt) {
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

	      // zero this out if phasepresent_ct == 0
	      const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

	      const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
	      const uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
	      uint32_t cur_shift = 1;
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uint32_t cur_geno = genovec_word & 3;
		write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
		if (cur_geno == 1) {
		  if (phasepresent_halfword & cur_shift) {
		    prev_phased_halfword |= cur_shift;
		    if (phaseinfo_halfword & cur_shift) {
		      memcpy(&(write_iter[-4]), "\t1|0", 4);
		    }
		  } else {
		    prev_phased_halfword &= ~cur_shift;
		  }
		}
		if (prev_phased_halfword & cur_shift) {
		  write_iter[-2] = '|';
		}
		if (dosage_present_hw & cur_shift) {
		  *write_iter++ = ':';
		  const uint32_t dosage_int = *dosage_vals_iter++;
		  write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
		}
		genovec_word >>= 2;
		cur_shift <<= 1;
	      }
	      ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
	      ++widx;
	    }
	  } else {
	    while (1) {
	      if (widx >= sample_ctl2_m1) {
		if (widx > sample_ctl2_m1) {
		  break;
		}
		inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
	      }
	      uintptr_t genovec_word = genovec[widx];
	      uint32_t is_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
	      uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

	      // zero this out if phasepresent_ct == 0
	      const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

	      const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
	      const uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
	      uint32_t cur_shift = 1;
	      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
		const uint32_t cur_geno = genovec_word & 3;
		const uint32_t cur_is_male = is_male_hw & 1;
		const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
		write_iter = memcpya(write_iter, haploid_genotext[cur_geno], cur_blen);
		if (cur_blen == 4) {
		  if (cur_geno == 1) {
		    if (phasepresent_halfword & cur_shift) {
		      prev_phased_halfword |= cur_shift;
		      if (phaseinfo_halfword & cur_shift) {
			memcpy(&(write_iter[-4]), "\t1|0", 4);
		      }
		    } else {
		      prev_phased_halfword &= ~cur_shift;
		    }
		  }
		  if (prev_phased_halfword & cur_shift) {
		    write_iter[-2] = '|';
		  }
		  if (dosage_present_hw & cur_shift) {
		    *write_iter++ = ':';
		    const uint32_t dosage_int = *dosage_vals_iter++;
		    write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
		  }
		} else {
		  if (dosage_present_hw & cur_shift) {
		    *write_iter++ = ':';
		    const uint32_t dosage_int = *dosage_vals_iter++;
		    if (write_ds) {
		      write_iter = haploid_dosage_print(dosage_int, write_iter);
		    } else {
		      write_iter = haploid_dosage_print(kDosageMax - dosage_int, write_iter);
		      *write_iter++ = ',';
		      write_iter = haploid_dosage_print(dosage_int, write_iter);
		    }
		  }
		}
		genovec_word >>= 2;
		is_male_hw >>= 1;
		cur_shift <<= 1;
	      }
	      ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
	      ++widx;
	    }
	  }
	}
      }
      // todo: multiallelic cases (separate out cur_allele_ct <= 10)
      append_binary_eoln(&write_iter);
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
	goto export_vcf_ret_WRITE_FAIL;
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
    if (write_iter != writebuf) {
      if (flexbwrite_flush(writebuf, write_iter - writebuf, outfile, bgz_outfile)) {
	goto export_vcf_ret_WRITE_FAIL;
      }
    }
    if (bgz_outfile) {
      if (bgzf_close(bgz_outfile)) {
	bgz_outfile = nullptr;
	goto export_vcf_ret_WRITE_FAIL;
      }
      bgz_outfile = nullptr;
    } else {
      if (fclose_null(&outfile)) {
	goto export_vcf_ret_WRITE_FAIL;
      }
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  export_vcf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_vcf_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_vcf_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_vcf_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  export_vcf_ret_PGR_FAIL:
    if (reterr != kPglRetReadFail) {
      logprint("\n");
      logerrprint("Error: Malformed .pgen file.\n");
    }
  }
 export_vcf_ret_1:
  fclose_cond(outfile);
  gzclose_cond(gz_pvar_reload);
  if (bgz_outfile) {
    bgzf_close(bgz_outfile);
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t exportf(char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, make_plink2_t make_plink2_modifier, exportf_flags_t exportf_modifier, idpaste_t exportf_id_paste, char exportf_id_delim, __maybe_unused uint32_t exportf_bits, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* sex_male_collapsed;
    if (bigstack_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
	bigstack_alloc_ul(sample_ctaw, &sex_male_collapsed)) {
      goto exportf_ret_NOMEM;
    }
    fill_cumulative_popcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed);
    fill_ulong_zero(sample_ctaw - sample_ctl, sex_male_collapsed);
    uint32_t* sample_missing_geno_cts = nullptr;
    if (exportf_modifier & (kfExportfOxGen | kfExportfHaps | kfExportfHapsLegend | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13)) {
      if (bigstack_alloc_ui(sample_ct, &sample_missing_geno_cts)) {
	goto exportf_ret_NOMEM;
      }
    }
    if (exportf_modifier & (kfExportf01 | kfExportf12)) {
      // todo
    }
    if (exportf_modifier & (kfExportfTypemask - kfExportfIndMajorBed - kfExportfVcf - kfExportfOxGen - kfExportfBgen11 - kfExportfHaps - kfExportfHapsLegend - kfExportfATranspose)) {
      logerrprint("Error: Only VCF, oxford, bgen-1.1, haps, hapslegend, A-transpose, and\nind-major-bed output have been implemented so far.\n");
      reterr = kPglRetNotYetSupported;
      goto exportf_ret_1;
    }
    const char exportf_delim = (exportf_modifier & kfExportfSpaces)? ' ' : '\t';
    if (exportf_modifier & kfExportfATranspose) {
      strcpy(outname_end, ".traw");
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_012_vmaj(outname, sample_include, sample_include_cumulative_popcounts, sample_ids, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, variant_cms, sample_ct, max_sample_id_blen, variant_ct, max_allele_slen, simple_pgrp);
      if (reterr) {
	goto exportf_ret_1;
      }
    }
    if (exportf_modifier & kfExportfIndMajorBed) {
      reterr = export_ind_major_bed(sample_include, variant_include, variant_allele_idxs, refalt1_select, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, outname, outname_end);
      if (reterr) {
	goto exportf_ret_1;
      }
    }
    if (exportf_modifier & kfExportfOxGen) {
      strcpy(outname_end, ".gen");
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_ox_gen(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, variant_ct, max_allele_slen, exportf_modifier, simple_pgrp, sample_missing_geno_cts);
      if (reterr) {
	goto exportf_ret_1;
      }
    }
    if (exportf_modifier & (kfExportfHaps | kfExportfHapsLegend)) {
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_ox_hapslegend(sample_include, sample_include_cumulative_popcounts, sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, exportf_modifier, simple_pgrp, outname, outname_end);
      if (reterr) {
	goto exportf_ret_1;
      }
      fill_uint_zero(sample_ct, sample_missing_geno_cts);
    }
    if (exportf_modifier & kfExportfBgen11) {
      assert(popcount_longs(sample_include, raw_sample_ctl) == sample_ct);
      strcpy(outname_end, ".bgen");
      reterr = export_bgen11(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, exportf_modifier, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (reterr) {
	goto exportf_ret_1;
      }
      /*
    } else if (exportf_modifier & (kfExportfBgen12 | kfExportfBgen13)) {
      strcpy(outname_end, ".bgen");
      reterr = export_bgen13(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, exportf_modifier, exportf_bits, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (reterr) {
	goto exportf_ret_1;
      }
      */
    }
    if (exportf_modifier & (kfExportfOxGen | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13 | kfExportfHaps | kfExportfHapsLegend)) {
      strcpy(outname_end, ".sample");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t y_ct = 0;
      int32_t y_code = cip->xymt_codes[kChrOffsetY];
      if ((y_code >= 0) && is_set(cip->chr_mask, y_code)) {
	y_ct = count_chr_variants_unsafe(variant_include, cip, y_code);
      }
      assert(popcount_longs(sample_include, raw_sample_ctl) == sample_ct);
      reterr = export_ox_sample(outname, sample_include, sample_ids, sample_missing_geno_cts, sex_nm, sex_male, pheno_cols, pheno_names, sample_ct, max_sample_id_blen, pheno_ct, max_pheno_name_blen, variant_ct, y_ct);
      if (reterr) {
	goto exportf_ret_1;
      }
      logprint("done.\n");
    }
    if (exportf_modifier & kfExportfVcf) {
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_vcf(xheader, sample_include, sample_include_cumulative_popcounts, sample_ids, sids, sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pvar_info_reload, xheader_blen, xheader_info_pr, sample_ct, max_sample_id_blen, max_sid_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, max_thread_ct, exportf_modifier, exportf_id_paste, exportf_id_delim, pgfip, simple_pgrp, outname, outname_end);
      if (reterr) {
	goto exportf_ret_1;
      }
    }
    // todo: everything else
    // sample-major output should share a (probably multithreaded) transpose
    // routine

    if ((!(make_plink2_modifier & kfMakeFam)) && (exportf_modifier & kfExportfIndMajorBed)) {
      strcpy(outname_end, ".fam");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_fam(outname, sample_include, sample_ids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, nullptr, sample_ct, max_sample_id_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, exportf_delim);
      if (reterr) {
	goto exportf_ret_1;
      }
      logprint("done.\n");
    }
    if ((!(make_plink2_modifier & kfMakeBim)) && (exportf_modifier & kfExportfIndMajorBed)) {
      strcpy(outname_end, ".bim");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_map_or_bim(outname, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, exportf_delim, 0);
      if (reterr) {
	goto exportf_ret_1;
      }
      logprint("done.\n");
    }
  }
  while (0) {
  exportf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 exportf_ret_1:
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
