// This file is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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
#include "plink2_pvar.h"

#ifdef __cplusplus
namespace plink2 {
#endif

pglerr_t write_map_or_bim(const char* outname, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst, uint32_t thread_ct) {
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
    if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
      goto write_map_or_bim_ret_NOMEM;
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = cswrite_init2(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto write_map_or_bim_ret_1;
    }

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t variant_uidx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
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
  write_map_or_bim_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 write_map_or_bim_ret_1:
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t pvar_info_reload_header(uintptr_t loadbuf_size, gzFile gz_pvar_reload, char* loadbuf, uint32_t* info_col_idx_ptr) {
  char* loadbuf_iter;
  do {
    // this is a reload, so no need to validate
    if (!gzgets(gz_pvar_reload, loadbuf, loadbuf_size)) {
      return kPglRetReadFail;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == kMaxLongLine) {
        return kPglRetReadFail;
      }
      return kPglRetNomem;
    }
    loadbuf_iter = skip_initial_spaces(loadbuf);
  } while (!str_startswith2(loadbuf_iter, "#CHROM"));
  uint32_t info_col_idx = 0;
  do {
    loadbuf_iter = next_token(loadbuf_iter);
    ++info_col_idx;
  } while (!tokequal_k(loadbuf_iter, "INFO"));
  *info_col_idx_ptr = info_col_idx;
  return kPglRetSuccess;
}

pglerr_t pvar_info_open_and_reload_header(const char* pvar_info_reload, gzFile* gz_pvar_reload_ptr, char** loadbuf_ptr, uintptr_t* loadbuf_size_ptr, uint32_t* info_col_idx_ptr) {
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
  *loadbuf_ptr = loadbuf;
  *loadbuf_size_ptr = loadbuf_size;
  return pvar_info_reload_header(loadbuf_size, *gz_pvar_reload_ptr, loadbuf, info_col_idx_ptr);
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

pglerr_t write_pvar(const char* outname, const char* xheader, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, char** filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, pvar_psam_t pvar_psam_modifier, uint32_t thread_ct) {
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

    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
        bigstack_alloc_ul(BITCT_TO_WORDCT(kPglMaxAltAlleleCt), &allele_include)) {
      goto write_pvar_ret_NOMEM;
    }
    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_modifier / kfPvarZs) & 1;
    reterr = cswrite_init2(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto write_pvar_ret_1;
    }
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_modifier & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
        if (variant_include[widx] & nonref_flags[widx]) {
          write_info_pr = 1;
          break;
        }
      }
    }
    write_info_pr = write_info_pr && write_info;
    if (write_info_pr && xheader_info_pr_nonflag) {
      logprint("\n");
      logerrprint("Error: Conflicting INFO:PR fields.  Either fix all REF alleles so that the\n'provisional reference' field is no longer needed, or remove/rename the other\nINFO:PR field.\n");
      goto write_pvar_ret_INCONSISTENT_INPUT;
    }

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
    }
    // bugfix (30 Jul 2017): may be necessary to reload INFO when no ## lines
    // are in the header
    if (pvar_info_reload) {
      reterr = pvar_info_open_and_reload_header(pvar_info_reload, &gz_pvar_reload, &loadbuf, &loadbuf_size, &info_col_idx);
      if (reterr) {
        goto write_pvar_ret_1;
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
    uint32_t chr_fo_idx = UINT32_MAX;
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
  write_pvar_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  write_pvar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
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
    if (pheno_idx != UINT32_MAX) {
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
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
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
      if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
        goto write_fam_ret_WRITE_FAIL;
      }
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
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
    if ((!strequal_k2(&(paternal_ids[sample_uidx * max_paternal_id_blen]), "0")) || (!strequal_k2(&(maternal_ids[sample_uidx * max_maternal_id_blen]), "0"))) {
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

    char* textbuf_flush = &(g_textbuf[kMaxMediumLine]);

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
    char* write_iter = strcpya(g_textbuf, "#FID\tIID");
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
        if (strequal_k(cur_pheno_name, "SEX", cur_pheno_name_slen)) {
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
        if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
          goto write_psam_ret_WRITE_FAIL;
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
          if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
            goto write_psam_ret_WRITE_FAIL;
          }
        }
      } else {
        if (write_empty_pheno) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
        }
        if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
          goto write_psam_ret_WRITE_FAIL;
        }
      }
      append_binary_eoln(&write_iter);
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
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

// a few multithread globals
static uint32_t g_sample_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_cur_block_write_ct = 0;
static uint32_t g_hard_call_halfdist = 0;
static uint32_t g_dosage_erase_halfdist = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

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
      partial_word_store(cur_word, QUATERCT_TO_BYTECT(cur_word_entry_ct), &(writebuf_walias[word_idx]));
      return;
    }
    writebuf_walias[word_idx++] = cur_word;
    new_sample_idx_to_old_iter = &(new_sample_idx_to_old_iter[kBitsPerWordD2]);
  }
}

void unpack_hphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, uint32_t raw_sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  const uint32_t het_ct = popcount_longs(all_hets, raw_sample_ctl);
  if (!(phaseraw[0] & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    expand_bytearr((const unsigned char*)phaseraw, all_hets, raw_sample_ctl, het_ct, 1, phaseinfo);
  } else {
    const unsigned char* phaseinfo_read = (const unsigned char*)(&(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]));
    expand_bytearr_nested(phaseinfo_read, phaseraw, all_hets, raw_sample_ctl, het_ct, 1, *phasepresent_ptr, phaseinfo);
  }
}

void unpack_hphase_subset(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
  const uint32_t het_ct = popcount_longs(all_hets, raw_sample_ctl);
  if (!(phaseraw[0] & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    expand_then_subset_bytearr((const unsigned char*)phaseraw, all_hets, sample_include, het_ct, sample_ct, 1, phaseinfo);
  } else {
    const uint32_t first_half_word_ct = 1 + (het_ct / kBitsPerWord);
    const uint32_t raw_phasepresent_ct = popcount_longs(phaseraw, first_half_word_ct) - 1;
    // see "if (explicit_phasepresent) {}" block in pgr_read_raw().  Could
    // change this convention.
    const unsigned char* phaseinfo_read = (const unsigned char*)(&(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]));
    expand_then_subset_bytearr_nested(phaseinfo_read, phaseraw, all_hets, sample_include, sample_ct, raw_phasepresent_ct, 1, *phasepresent_ptr, phaseinfo);
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
          new_phasepresent_word &= new_phasepresent_word - 1;
        }
      } else {
        uintptr_t masked_phasepresent_word = new_phasepresent_word & sample_include[widx];
        while (masked_phasepresent_word) {
          const uint32_t sample_uidx_lowbits = CTZLU(masked_phasepresent_word);
          const uintptr_t lowmask = (k1LU << sample_uidx_lowbits) - k1LU;
          if ((tmp_phaseinfo_input_word >> popcount_long(new_phasepresent_word & lowmask)) & 1) {
            SET_BIT(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
          }
          masked_phasepresent_word &= masked_phasepresent_word - 1;
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
      tmp_phasepresent_input_word = bzhi_max(tmp_phasepresent_input_word, read_idx_lowbits_end - read_idx_lowbits);
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
            geno_hets &= geno_hets - 1;
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
            geno_hets &= geno_hets - 1;
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
static uintptr_t* g_sex_female_collapsed = nullptr;
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

static uintptr_t** g_thread_write_genovecs = nullptr;
static uintptr_t** g_thread_write_phasepresents = nullptr;
static uintptr_t** g_thread_write_phaseinfos = nullptr;
static uintptr_t** g_thread_all_hets = nullptr;
static uintptr_t** g_thread_write_dosagepresents = nullptr;
static dosage_t** g_thread_write_dosagevals = nullptr;
static uint32_t** g_thread_cumulative_popcount_bufs = nullptr;
static pgen_writer_common_t** g_pwcs = nullptr;

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
      zero_trailing_words(raw_sample_ctl, g_sex_male);
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
      zero_trailing_words(raw_sample_ctl, nosex_buf);
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
        zero_trailing_words(raw_sample_ctl, g_founder_male);
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
          zero_trailing_words(raw_sample_ctl, nosex_buf);
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

void apply_hard_call_thresh(const uintptr_t* dosage_present, const dosage_t* dosage_vals, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec) {
  uint32_t sample_uidx = 0;
  for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
    next_set_unsafe_ck(dosage_present, &sample_uidx);
    const uint32_t dosage_int = dosage_vals[dosage_idx];
    const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
    uintptr_t new_geno;
    if (halfdist < hard_call_halfdist) {
      new_geno = 3;
    } else {
      new_geno = (dosage_int + kDosage4th) / kDosageMid;
    }
    ASSIGN_QUATERARR_ENTRY(sample_uidx, new_geno, genovec);
  }
}

FLAGSET_DEF_START()
  kfPlink2Write0,
  kfPlink2WriteSetHhMissing = (1 << 0),
  kfPlink2WriteSetHhMissingKeepDosage = (1 << 1),
  kfPlink2WriteSetMixedMtMissing = (1 << 2),
  kfPlink2WriteSetMixedMtMissingKeepDosage = (1 << 3),
  kfPlink2WriteMeMissing = (1 << 4),
  kfPlink2WriteZeroCluster = (1 << 5),
  kfPlink2WriteFillRef = (1 << 6),
  kfPlink2WriteLateDosageErase = (1 << 7),
  // no need for sample_sort, determined by g_collapsed_sort_map != nullptr?
  kfPlink2WritePlink1 = (1 << 8)
FLAGSET_DEF_END(plink2_write_flags_t);
// todo: add .pgen-specific stuff

// more multithread globals
static plink2_write_flags_t g_plink2_write_flags = kfPlink2Write0;

THREAD_FUNC_DECL make_bedlike_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  dosage_t* dosage_vals = nullptr;
  uint32_t hard_call_halfdist = 0;
  if (g_dosage_presents) {
    dosage_present = g_dosage_presents[tidx];
    dosage_vals = g_dosage_val_bufs[tidx];
    hard_call_halfdist = g_hard_call_halfdist;
  }
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
      if (!hard_call_halfdist) {
        pglerr_t reterr = pgr_read_refalt1_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec);
        if (reterr) {
          g_error_ret = reterr;
          break;
        }
      } else {
        // quasi-bugfix (4 Dec 2017): it's user-hostile to make
        // --hard-call-threshold not apply here.
        uint32_t dosage_ct;
        uint32_t is_explicit_alt1;
        pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
        if (reterr) {
          g_error_ret = reterr;
          break;
        }
        apply_hard_call_thresh(dosage_present, dosage_vals, dosage_ct, hard_call_halfdist, genovec);
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

// more multithread globals

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

  // er, this global should be named g_sex_male_collapsed...
  const uintptr_t* sex_male_collapsed = g_sex_male;

  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed = g_sex_female_collapsed;
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
  const uint32_t set_hh_missing_keep_dosage = g_plink2_write_flags & kfPlink2WriteSetHhMissingKeepDosage;
  const uint32_t set_mixed_mt_missing = g_plink2_write_flags & kfPlink2WriteSetMixedMtMissing;
  const uint32_t set_mixed_mt_missing_keep_dosage = g_plink2_write_flags & kfPlink2WriteSetMixedMtMissingKeepDosage;
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
  if (dosage_present || set_hh_missing_keep_dosage || set_mixed_mt_missing_keep_dosage) {
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
            unpack_hphase_subset(all_hets, cur_phaseraw, sample_include, raw_sample_ct, sample_ct, &cur_write_phasepresent, write_phaseinfo);
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
            if (!set_hh_missing_keep_dosage) {
              // need to erase dosages associated with the hardcalls we're
              // about to clear

              // male hets to missing
              set_male_het_missing_cleardosage(sex_male_collapsed, sex_male_collapsed_interleaved, sample_ctv2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            } else {
              // need to generate a new unphased dosage for each cleared
              // hardcall lacking a dosage entry
              set_male_het_missing_keepdosage(sex_male_collapsed, sex_male_collapsed_interleaved, sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            }
            if (is_y) {
              // all female calls to missing; unknown-sex calls now left
              // alone
              interleaved_set_missing_cleardosage(sex_female_collapsed, sex_female_collapsed_interleaved, sample_ctv2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
              is_hphase = 0;
            } else if (is_hphase && cur_write_phasepresent) {
              mask_genovec_hets_unsafe(write_genovec, sample_ctl2, cur_write_phasepresent);
            }
          } else {
            // all hets to missing
            // may want to move is_hphase zeroing in front
            if (!set_hh_missing_keep_dosage) {
              set_het_missing_cleardosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            } else {
              set_het_missing_keepdosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            }
            is_hphase = 0;
          }
        } else if (set_mixed_mt_missing && is_mt) {
          if (!set_mixed_mt_missing_keep_dosage) {
            // all hets to missing
            set_het_missing_cleardosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          } else {
            // todo: keep dosage-phase information
            set_het_missing_keepdosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          }
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
        // bugfix (8 Oct 2017): forgot to multiply by kMask0101
        if (vrtypes_or & (0x10 * kMask0101)) {
          read_phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent;
          mask_multiply -= 0x10;
        }
        if (vrtypes_or & (0x60 * kMask0101)) {
          read_phase_dosage_gflags |= kfPgenGlobalDosagePresent;
          mask_multiply -= 0x60;
        }
        if (vrtypes_or & (0x80 * kMask0101)) {
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
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
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
        if (hard_call_thresh == UINT32_MAX) {
          read_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
          g_plink2_write_flags |= kfPlink2WriteLateDosageErase;
        }
      }
      if (read_phase_dosage_gflags && (variant_ct < raw_variant_ct)) {
        read_phase_dosage_gflags = gflags_vfilter(variant_include, simple_pgrp->fi.vrtypes, raw_variant_ct, simple_pgrp->fi.gflags);
      }
      g_read_phase_dosage_gflags = read_phase_dosage_gflags;
      g_hard_call_halfdist = (hard_call_thresh == UINT32_MAX)? 0 : (kDosage4th - hard_call_thresh);
      g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_hphase_present = (read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
      const uint32_t read_dosage_present = (read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
      pgen_global_flags_t write_phase_dosage_gflags = read_phase_dosage_gflags;
      uint32_t read_or_write_dosage_present = read_dosage_present;
      if (g_plink2_write_flags & kfPlink2WriteLateDosageErase) {
        write_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      } else if (make_plink2_modifier & (kfPlink2WriteSetHhMissingKeepDosage | kfPlink2WriteSetMixedMtMissingKeepDosage)) {
        // command-line parser guarantees erase-dosage and
        // --set-hh-missing/--set-mixed-mt-missing keep-dosage aren't used
        // together
        read_or_write_dosage_present = 1;

        // could verify at least one het haploid is present before setting this
        // flag...
        write_phase_dosage_gflags |= kfPgenGlobalDosagePresent;
        if (make_plink2_modifier & kfPlink2WriteSetMixedMtMissingKeepDosage) {
          // todo: keep MT phasing when phased dosage supported
          // write_phase_dosage_gflags |= kfPgenGlobalDosagePhasePresent;
        }
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
      if (read_or_write_dosage_present) {
        if (bigstack_alloc_dosagep(1, &g_thread_write_dosagevals) ||
            bigstack_alloc_ulp(1, &g_thread_write_dosagepresents) ||
            bigstack_alloc_ul(sample_ctl, &(g_thread_write_dosagepresents[0])) ||
            bigstack_alloc_dosage(sample_ct, &(g_thread_write_dosagevals[0]))) {
          goto make_pgen_robust_ret_NOMEM;
        }
        if (read_dosage_present && new_sample_idx_to_old) {
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
          //    bit is set iff phasepresent is explicitly stored at all (if
          //    not, all hets are assumed to be phased), if yes the remaining
          //    bits store packed phasepresent values for all hets, if no the
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
      uint32_t read_variant_uidx = UINT32_MAX; // deliberate overflow
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

pglerr_t make_plink2_no_vsort(const char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, pvar_psam_t pvar_psam_modifier, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
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
      uintptr_t* sex_female;
      if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male) ||
          bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
          bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed) ||
          bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed_interleaved) ||
          bigstack_alloc_ul(raw_sample_ctl, &sex_female)) {
        goto make_plink2_no_vsort_ret_NOMEM;
      }
      copy_bitarr_subset(sex_male, sample_include, sample_ct, g_sex_male);
      zero_trailing_words(BITCT_TO_WORDCT(sample_ct), g_sex_male);
      fill_interleaved_mask_vec(g_sex_male, sample_ctv, g_sex_male_collapsed_interleaved);

      bitvec_andnot_copy(sex_nm, sex_male, raw_sample_ctl, sex_female);
      copy_bitarr_subset(sex_female, sample_include, sample_ct, g_sex_female_collapsed);
      fill_interleaved_mask_vec(g_sex_female_collapsed, sample_ctv, g_sex_female_collapsed_interleaved);

      bigstack_reset(sex_female);
      g_plink2_write_flags |= kfPlink2WriteSetHhMissing;
      if (make_plink2_modifier & kfMakePlink2SetHhMissingKeepDosage) {
        g_plink2_write_flags |= kfPlink2WriteSetHhMissingKeepDosage;
      }
    }
    if (make_plink2_modifier & kfMakePlink2SetMixedMtMissing) {
      g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
      if (make_plink2_modifier & kfMakePlink2SetMixedMtMissingKeepDosage) {
        g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissingKeepDosage;
      }
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
        fwrite_unlocked("l\x1b\x02", 3, 1, outfile);
        fwrite_unlocked(&variant_ct, 4, 1, outfile);
        fwrite_unlocked(&sample_ct, 4, 1, outfile);
        if (!pgfip->nonref_flags) {
          const pgen_global_flags_t gflags = pgfip->gflags;
          uint32_t uii = 64;
          if (gflags & kfPgenGlobalAllNonref) {
            uii = 128;
          }
          putc_unlocked(uii, outfile);
        } else {
          putc_unlocked(192, outfile);
          fwrite_unlocked(pgfip->nonref_flags, DIV_UP(variant_ct, CHAR_BIT), 1, outfile);
        }
        if (ferror_unlocked(outfile)) {
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
          // Without BMI2 instructions, subsetting is most expensive with
          // sample_ct near 2/3 of raw_sample_ct; up to ~7 compute threads are
          // useful in that case.  (See copy_quaterarr_nonempty_subset().)
          // With them, 1-2 compute threads appear to suffice.
#ifdef USE_AVX2
          const uint32_t calc_thread_max = 2;
#else
          uint64_t numer;
          if (sample_ct * (3 * k1LU) <= raw_sample_ct * (2 * k1LU)) {
            numer = sample_ct * (9 * k1LU);
          } else {
            numer = (raw_sample_ct - sample_ct) * (18 * k1LU);
          }
          const uint32_t calc_thread_max = 1 + (numer / raw_sample_ct);
#endif
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

        g_hard_call_halfdist = 0;
        if ((hard_call_thresh != UINT32_MAX) && (pgfip->gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))) {
          g_hard_call_halfdist = kDosage4th - hard_call_thresh;
        }
        unsigned char* main_loadbufs[2];
        uint32_t read_block_size;
        if (multithread_load_init(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, g_hard_call_halfdist? (&g_dosage_presents) : nullptr, g_hard_call_halfdist? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &ts.threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
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
        if (hard_call_thresh == UINT32_MAX) {
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
      g_hard_call_halfdist = (hard_call_thresh == UINT32_MAX)? 0 : (kDosage4th - hard_call_thresh);
      g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_hphase_present = (read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
      const uint32_t read_dosage_present = (read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
      pgen_global_flags_t write_phase_dosage_gflags = read_phase_dosage_gflags;
      uint32_t read_or_write_dosage_present = read_dosage_present;
      if (g_plink2_write_flags & kfPlink2WriteLateDosageErase) {
        write_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      } else if (g_plink2_write_flags & (kfPlink2WriteSetHhMissingKeepDosage | kfPlink2WriteSetMixedMtMissingKeepDosage)) {
        read_or_write_dosage_present = 1;
        write_phase_dosage_gflags |= kfPgenGlobalDosagePresent;
        if (g_plink2_write_flags & kfPlink2WriteSetMixedMtMissingKeepDosage) {
          // todo: keep MT phasing when phased dosage supported
          // write_phase_dosage_gflags |= kfPgenGlobalDosagePhasePresent;
        }
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
#ifdef USE_AVX2
        const uint32_t max_calc_thread_ct = 2;
#else
        const uint32_t max_calc_thread_ct = 2 + subsetting_required;
#endif
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
      if (read_hphase_present || read_or_write_dosage_present) {
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
        if (read_or_write_dosage_present) {
          if (bigstack_alloc_dosagep(calc_thread_ct, &g_thread_write_dosagevals) ||
              bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_dosagepresents)) {
            goto make_plink2_no_vsort_fallback;
          }
          if (read_dosage_present && new_sample_idx_to_old) {
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
      if (read_hphase_present || read_or_write_dosage_present) {
        if (read_hphase_present || read_dosage_present) {
          g_loaded_vrtypes[0] = bigstack_alloc_raw(kPglVblockSize * calc_thread_ct);
          g_loaded_vrtypes[1] = bigstack_alloc_raw(kPglVblockSize * calc_thread_ct);
        }
        const uint32_t bitvec_writebuf_byte_ct = BITCT_TO_CLCT(sample_ct) * kCacheline;
        const uintptr_t dosagevals_writebuf_byte_ct = DIV_UP(sample_ct, (kCacheline / 2)) * kCacheline;
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          if (read_hphase_present) {
            g_thread_write_phasepresents[tidx] = (uintptr_t*)bigstack_alloc_raw(bitvec_writebuf_byte_ct);
            g_thread_write_phaseinfos[tidx] = (uintptr_t*)bigstack_alloc_raw(bitvec_writebuf_byte_ct);

            g_thread_all_hets[tidx] = (uintptr_t*)bigstack_alloc_raw(BITCT_TO_CLCT(raw_sample_ct) * kCacheline);
          }
          if (read_or_write_dosage_present) {
            g_thread_write_dosagepresents[tidx] = (uintptr_t*)bigstack_alloc_raw(bitvec_writebuf_byte_ct);
            g_thread_write_dosagevals[tidx] = (dosage_t*)bigstack_alloc_raw(dosagevals_writebuf_byte_ct);
            if (read_dosage_present && new_sample_idx_to_old) {
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
      reterr = write_map_or_bim(outname, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, '\t', bim_zst, max_thread_ct);
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
      reterr = write_pvar(outname, xheader, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, xheader_info_pr, xheader_info_pr_nonflag, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_modifier, max_thread_ct);
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


boolerr_t sort_chr(const chr_info_t* cip, const uint32_t* chr_idx_to_size, uint32_t use_nsort, chr_info_t* write_cip) {
  // Finishes initialization of write_cip.  Assumes chr_fo_vidx_start is
  // allocated and initialized to all-bits-one, chr_file_order/chr_idx_to_foidx
  // are unallocated, and chr_ct is uninitialized.
  const uint32_t max_code = cip->max_code;
  const uint32_t chr_code_end = max_code + 1 + cip->name_ct;
  uint32_t new_chr_ct = 0;
  for (uint32_t chr_idx = 0; chr_idx < chr_code_end; ++chr_idx) {
    const uint32_t cur_chr_size = chr_idx_to_size[chr_idx];
    if (cur_chr_size) {
      ++new_chr_ct;
    }
  }
  if (bigstack_alloc_ui(new_chr_ct, &(write_cip->chr_file_order)) ||
      bigstack_alloc_ui(new_chr_ct, &(write_cip->chr_fo_vidx_start))) {
    return 1;
  }
  write_cip->chr_ct = new_chr_ct;
  // now for the actual sorting.
  // autosomes and PAR1/X/PAR2/Y/XY/MT come first, then contig names.
  const uint32_t autosome_ct = cip->autosome_ct;
  const uint32_t xymt_ct = max_code - autosome_ct;
  const uint32_t autosome_ct_p1 = autosome_ct + 1;

  const int32_t* xymt_codes = cip->xymt_codes;
  uintptr_t xymt_idx_to_chr_sort_offset[kChrOffsetCt] = {1, 3, 4, 5, 0, 2};

  // chr_sort_idx in high bits, original chr_idx in low
  uint64_t* std_sortbuf;
  uint64_t* std_sortbuf_iter;
  if (bigstack_alloc_ull(max_code + 1, &std_sortbuf)) {
    return 1;
  }
  std_sortbuf_iter = std_sortbuf;
  for (uintptr_t chr_idx = 0; chr_idx <= autosome_ct; ++chr_idx) {
    if (chr_idx_to_size[chr_idx]) {
      *std_sortbuf_iter++ = chr_idx * 0x100000001LLU;
    }
  }
  for (uint32_t xymt_idx = 0; xymt_idx < xymt_ct; ++xymt_idx) {
    const int32_t xymt_code = xymt_codes[xymt_idx];
    if (xymt_code >= 0) {
      if (chr_idx_to_size[xymt_idx + autosome_ct_p1]) {
        *std_sortbuf_iter++ = (((uint64_t)(xymt_idx_to_chr_sort_offset[xymt_idx] + autosome_ct_p1)) << 32) | (xymt_idx + autosome_ct_p1);
      }
    }
  }
  const uint32_t std_sortbuf_len = (uintptr_t)(std_sortbuf_iter - std_sortbuf);
#ifdef __cplusplus
  std::sort(std_sortbuf, &(std_sortbuf[std_sortbuf_len]));
#else
  qsort(std_sortbuf, std_sortbuf_len, sizeof(int64_t), uint64cmp);
#endif
  uint32_t write_vidx = 0;
  write_cip->chr_fo_vidx_start[0] = 0;
  for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx < std_sortbuf_len; ++new_chr_fo_idx) {
    const uint64_t cur_entry = std_sortbuf[new_chr_fo_idx];
    const uintptr_t chr_idx = (uint32_t)cur_entry;
    const uint32_t chr_size = chr_idx_to_size[chr_idx];
    write_cip->chr_file_order[new_chr_fo_idx] = chr_idx;
    write_vidx += chr_size;
    write_cip->chr_fo_vidx_start[new_chr_fo_idx + 1] = write_vidx;
    write_cip->chr_idx_to_foidx[chr_idx] = new_chr_fo_idx;
  }

  const uint32_t new_nonstd_ct = new_chr_ct - std_sortbuf_len;
  if (new_nonstd_ct) {
    str_sort_indexed_deref_t* nonstd_sort_buf = (str_sort_indexed_deref_t*)bigstack_alloc_raw_rd(new_nonstd_ct * sizeof(str_sort_indexed_deref_t));
    if (!nonstd_sort_buf) {
      return 1;
    }
    char** nonstd_names = cip->nonstd_names;
    uint32_t str_idx = 0;
    for (uint32_t chr_idx = max_code + 1; chr_idx < chr_code_end; ++chr_idx) {
      if (chr_idx_to_size[chr_idx]) {
        nonstd_sort_buf[str_idx].strptr = nonstd_names[chr_idx];
        nonstd_sort_buf[str_idx].orig_idx = chr_idx;
        ++str_idx;
      }
    }
    assert(str_idx == new_nonstd_ct);
    strptr_arr_sort_main(new_nonstd_ct, use_nsort, nonstd_sort_buf);
    uint32_t new_chr_fo_idx = std_sortbuf_len;
    for (uint32_t str_idx = 0; str_idx < new_nonstd_ct; ++str_idx, ++new_chr_fo_idx) {
      const uint32_t chr_idx = nonstd_sort_buf[str_idx].orig_idx;
      const uint32_t chr_size = chr_idx_to_size[chr_idx];
      write_cip->chr_file_order[new_chr_fo_idx] = chr_idx;
      write_vidx += chr_size;
      write_cip->chr_fo_vidx_start[new_chr_fo_idx + 1] = write_vidx;
      write_cip->chr_idx_to_foidx[chr_idx] = new_chr_fo_idx;
    }
  }
  bigstack_reset(std_sortbuf);
  return 0;
}

// hybrid of write_map_or_bim() and write_pvar_resorted()
pglerr_t write_bim_resorted(const char* outname, const chr_info_t* write_cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t output_zst, uint32_t thread_ct) {
  // allele_dosages must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t max_chr_blen = get_max_chr_slen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
      goto write_bim_resorted_ret_NOMEM;
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = cswrite_init2(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto write_bim_resorted_ret_1;
    }

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = new_variant_idx_to_old[variant_idx];
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
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
      if (!variant_cms) {
        *cswritep++ = '0';
      } else {
        cswritep = dtoa_g_p8(variant_cms[variant_uidx], cswritep);
      }
      *cswritep++ = '\t';
      cswritep = uint32toa(variant_bps[variant_uidx], cswritep);
      *cswritep++ = '\t';
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
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[0]);
      } else {
        const alt_allele_ct_t* cur_refalt1_select = &(refalt1_select[variant_uidx * 2]);
        if ((!allele_dosages) || allele_dosages[cur_refalt1_select[1] + variant_allele_idx_base]) {
          cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[1]]);
        } else {
          *cswritep++ = output_missing_geno_char;
        }
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[0]]);
      }
      append_binary_eoln(&cswritep);
      if (cswrite(&css, &cswritep)) {
        goto write_bim_resorted_ret_WRITE_FAIL;
      }
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto write_bim_resorted_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_bim_resorted_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_bim_resorted_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 write_bim_resorted_ret_1:
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t pvar_info_reload_interval(const uint32_t* old_variant_uidx_to_new, uintptr_t loadbuf_size, uint32_t variant_idx_start, uint32_t variant_idx_end, uint32_t raw_variant_ct, gzFile gz_pvar_reload, char** pvar_info_strs, char* loadbuf) {
  // We assume the batch size was chosen such that there's no risk of
  // scribbling past g_bigstack_end (barring pathological cases like another
  // process modifying the .pvar file after initial load).
  // We also assume no more dynamic allocations are needed after this;
  // otherwise, str_store_iter should be returned.
  if (gzrewind(gz_pvar_reload)) {
    return kPglRetReadFail;
  }
  const uint32_t cur_batch_size = variant_idx_end - variant_idx_start;
  char* str_store_iter = (char*)g_bigstack_base;
  uint32_t info_col_idx;
  pglerr_t reterr = pvar_info_reload_header(loadbuf_size, gz_pvar_reload, loadbuf, &info_col_idx);
  if (reterr) {
    return reterr;
  }
  for (uint32_t variant_uidx = 0; variant_uidx < raw_variant_ct; ++variant_uidx) {
    char* loadbuf_iter;
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
      loadbuf_iter = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_iter));
    const uint32_t new_variant_idx_offset = old_variant_uidx_to_new[variant_uidx] - variant_idx_start;
    // exploit wraparound, UINT32_MAX null value
    if (new_variant_idx_offset >= cur_batch_size) {
      continue;
    }
    loadbuf_iter = next_token_mult(loadbuf_iter, info_col_idx);
    if (!loadbuf_iter) {
      return kPglRetReadFail;
    }
    const uint32_t info_slen = strlen_se(loadbuf_iter);
    pvar_info_strs[new_variant_idx_offset] = str_store_iter;
    str_store_iter = memcpyax(str_store_iter, loadbuf_iter, info_slen, '\0');
  }
  assert(str_store_iter <= (char*)g_bigstack_end);
  return kPglRetSuccess;
}

// could be boolerr_t
pglerr_t write_pvar_resorted_interval(const chr_info_t* write_cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, char** filter_storage, const uintptr_t* nonref_flags, const double* variant_cms, const uint32_t* new_variant_idx_to_old, char** pvar_info_strs, uint32_t variant_idx_start, uint32_t variant_idx_end, uint32_t xheader_info_pr, uint32_t write_qual, uint32_t write_filter, uint32_t write_info, uint32_t all_nonref, uint32_t write_cm, compress_stream_state_t* cssp, char** cswritepp, uint32_t* chr_fo_idxp, uint32_t* chr_endp, uint32_t* chr_buf_blenp, char* chr_buf, uintptr_t* allele_include) {
  char* cswritep = *cswritepp;
  uint32_t chr_fo_idx = *chr_fo_idxp;
  uint32_t chr_end = *chr_endp;
  uint32_t chr_buf_blen = *chr_buf_blenp;
  pglerr_t reterr = kPglRetSuccess;
  {
    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = variant_idx_start; variant_idx < variant_idx_end; ++variant_idx) {
      const uint32_t variant_uidx = new_variant_idx_to_old[variant_idx];
      if (variant_idx == chr_end) {
        ++chr_fo_idx;
        chr_end = write_cip->chr_fo_vidx_start[chr_fo_idx + 1];
        assert(variant_idx < chr_end);
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
      if (cswrite(cssp, &cswritep)) {
        goto write_pvar_resorted_interval_ret_WRITE_FAIL;
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
          if (cswrite(cssp, &cswritep)) {
            goto write_pvar_resorted_interval_ret_WRITE_FAIL;
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
        if (pvar_info_strs) {
          pvar_info_write(pvar_info_strs[variant_idx - variant_idx_start], xheader_info_pr, is_pr, &cswritep);
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

  }
  while (0) {
  write_pvar_resorted_interval_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  *cswritepp = cswritep;
  *chr_fo_idxp = chr_fo_idx;
  *chr_endp = chr_end;
  *chr_buf_blenp = chr_buf_blen;
  return reterr;
}

// allele_dosages must be nullptr unless we're trimming alt alleles.
//
// The annoying part of this is handling a sequence of INFO strings that don't
// fit in memory; we use a multipass approach for that.  File creation,
// allocation of buffers, and generating the header line occurs directly in
// this function, while loading the next pvar_info_strs batch and writing the
// next .pvar line batch are one level down.
pglerr_t write_pvar_resorted(const char* outname, const char* xheader, const uintptr_t* variant_include, const chr_info_t* write_cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, char** filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, pvar_psam_t pvar_psam_modifier, uint32_t thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  gzFile gz_pvar_reload = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    const uint32_t max_chr_blen = get_max_chr_slen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
        bigstack_alloc_ul(BITCT_TO_WORDCT(kPglMaxAltAlleleCt), &allele_include)) {
      goto write_pvar_resorted_ret_NOMEM;
    }
    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_modifier / kfPvarZs) & 1;
    reterr = cswrite_init2(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto write_pvar_resorted_ret_1;
    }
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_modifier & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
        if (variant_include[widx] & nonref_flags[widx]) {
          write_info_pr = 1;
          break;
        }
      }
    }
    write_info_pr = write_info_pr && write_info;
    if (write_info_pr && xheader_info_pr_nonflag) {
      logprint("\n");
      logerrprint("Error: Conflicting INFO:PR fields.  Either fix all REF alleles so that the\n'provisional reference' field is no longer needed, or remove/rename the other\nINFO:PR field.\n");
      goto write_pvar_resorted_ret_INCONSISTENT_INPUT;
    }

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

    uint32_t* old_variant_uidx_to_new = nullptr;
    char** pvar_info_strs = nullptr;
    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t batch_size = variant_ct;
    uint32_t batch_ct = 1;
    if (pvar_info_reload) {
      if (bigstack_alloc_ui(raw_variant_ct, &old_variant_uidx_to_new)) {
        goto write_pvar_resorted_ret_NOMEM;
      }
      fill_uint_one(raw_variant_ct, old_variant_uidx_to_new);
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
        const uint32_t old_variant_uidx = new_variant_idx_to_old[variant_idx];
        old_variant_uidx_to_new[old_variant_uidx] = variant_idx;
      }

      reterr = gzopen_read_checked(pvar_info_reload, &gz_pvar_reload);
      if (reterr) {
        return reterr;
      }
      loadbuf_size = round_down_pow2(bigstack_left() / 4, kCacheline);
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size <= kMaxMediumLine) {
        return kPglRetNomem;
      }
      loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
      loadbuf[loadbuf_size - 1] = ' ';

      // subtract kCacheline to allow for rounding
      uintptr_t bytes_left = bigstack_left() - kCacheline;
      uint32_t single_variant_byte_ct = info_reload_slen + 1 + sizeof(intptr_t);
      if (variant_ct * single_variant_byte_ct > bytes_left) {
        batch_size = bytes_left / single_variant_byte_ct;
        batch_ct = 1 + (variant_ct - 1) / batch_size;
      }
      pvar_info_strs = (char**)bigstack_alloc_raw_rd(batch_size * sizeof(intptr_t));
    }

    uint32_t variant_idx_start = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t batch_idx = 0; batch_idx < batch_ct; ++batch_idx) {
      uint32_t variant_idx_end = MINV(variant_idx_start + batch_size, variant_ct);
      if (gz_pvar_reload) {
        reterr = pvar_info_reload_interval(old_variant_uidx_to_new, loadbuf_size, variant_idx_start, variant_idx_end, raw_variant_ct, gz_pvar_reload, pvar_info_strs, loadbuf);
        if (reterr) {
          goto write_pvar_resorted_ret_1;
        }
      }
      reterr = write_pvar_resorted_interval(write_cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, allele_dosages, refalt1_select, qual_present, quals, filter_present, filter_npass, filter_storage, nonref_flags, variant_cms, new_variant_idx_to_old, pvar_info_strs, variant_idx_start, variant_idx_end, xheader_info_pr, write_qual, write_filter, write_info, all_nonref, write_cm, &css, &cswritep, &chr_fo_idx, &chr_end, &chr_buf_blen, chr_buf, allele_include);
      if (reterr) {
        goto write_pvar_resorted_ret_1;
      }
      variant_idx_start = variant_idx_end;
    }

    if (cswrite_close_null(&css, cswritep)) {
      goto write_pvar_resorted_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_pvar_resorted_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  write_pvar_resorted_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  write_pvar_resorted_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 write_pvar_resorted_ret_1:
  cswrite_close_cond(&css, cswritep);
  gzclose_cond(gz_pvar_reload);
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t make_plink2_vsort(const char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const uint64_t* allele_dosages, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const chr_idx_t* chr_idxs, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, make_plink2_t make_plink2_modifier, uint32_t use_nsort, pvar_psam_t pvar_psam_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  pglerr_t reterr = kPglRetSuccess;
  {
    // Resort the variants.
    // 1. (todo) Apply --update-chr if necessary.
    // 2. Count number of remaining variants in each chromosome, then sort the
    //    chromosomes.
    // 3. Within each chromosome, sort by position, effectively adding 0.5 for
    //    non-SNPs.  (could multithread this by chromosome, and/or use C++17
    //    multithreaded sort)
    // 4. Scan for position ties, sort on ID (according to --sort-vars setting,
    //    defaults to natural-sort but can be ASCII).
    // 5. Fill new_variant_idx_to_old, free sort buffers.

    // possible todo: put this in a "copy constructor" function
    chr_info_t write_chr_info;

    // const_casts
    write_chr_info.haploid_mask = (uintptr_t*)((uintptr_t)cip->haploid_mask);
    write_chr_info.nonstd_names = (char**)((uintptr_t)cip->nonstd_names);
    write_chr_info.nonstd_id_htable = (uint32_t*)((uintptr_t)cip->nonstd_id_htable);

    write_chr_info.chrset_source = cip->chrset_source;
    memcpy(write_chr_info.chr_exclude, cip->chr_exclude, kChrExcludeWords * sizeof(intptr_t));
    memcpy(write_chr_info.xymt_codes, cip->xymt_codes, kChrOffsetCt * sizeof(int32_t));
    write_chr_info.max_numeric_code = cip->max_numeric_code;
    write_chr_info.max_code = cip->max_code;
    write_chr_info.autosome_ct = cip->autosome_ct;
    write_chr_info.zero_extra_chrs = cip->zero_extra_chrs;
    write_chr_info.name_ct = cip->name_ct;
    // const_cast
    write_chr_info.incl_excl_name_stack = (ll_str_t*)((uintptr_t)cip->incl_excl_name_stack);
    write_chr_info.is_include_stack = cip->is_include_stack;
    write_chr_info.output_encoding = cip->output_encoding;

    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    uint32_t* chr_idx_to_size;
    if (bigstack_calloc_ul(kChrMaskWords, &write_chr_info.chr_mask) ||
        bigstack_alloc_ui(chr_code_end, &write_chr_info.chr_idx_to_foidx) ||
        bigstack_end_calloc_ui(chr_code_end, &chr_idx_to_size)) {
      goto make_plink2_vsort_ret_NOMEM;
    }
    fill_uint_one(chr_code_end, write_chr_info.chr_idx_to_foidx);
    if (chr_idxs) {
      uint32_t variant_uidx = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        next_set_unsafe_ck(variant_include, &variant_uidx);
        chr_idx_to_size[chr_idxs[variant_uidx]] += 1;
      }
      for (uint32_t chr_idx = 0; chr_idx < chr_code_end; ++chr_idx) {
        if (chr_idx_to_size[chr_idx]) {
          set_bit(chr_idx, write_chr_info.chr_mask);
        }
      }
      // bugfix: chr_file_order is invalid
    } else {
      const uint32_t* chr_fo_vidx_start = cip->chr_fo_vidx_start;
      const uint32_t orig_chr_ct = cip->chr_ct;
      uint32_t vidx_start = 0;
      for (uint32_t chr_fo_idx = 0; chr_fo_idx < orig_chr_ct; ++chr_fo_idx) {
        const uint32_t vidx_end = chr_fo_vidx_start[chr_fo_idx + 1];
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_idx_to_size[chr_idx] = popcount_bit_idx(variant_include, vidx_start, vidx_end);
        if (chr_idx_to_size[chr_idx]) {
          set_bit(chr_idx, write_chr_info.chr_mask);
        }
        vidx_start = vidx_end;
      }
    }
    if (sort_chr(cip, chr_idx_to_size, use_nsort, &write_chr_info)) {
      goto make_plink2_vsort_ret_NOMEM;
    }

    uint32_t* new_variant_idx_to_old;

    // pos_vidx_sort_buf has variant_bp in high bits, variant_uidx in low
    uint64_t* pos_vidx_sort_buf;
    if (bigstack_alloc_ui(variant_ct, &new_variant_idx_to_old) ||
        bigstack_alloc_ull(variant_ct + 1, &pos_vidx_sort_buf)) {
      goto make_plink2_vsort_ret_NOMEM;
    }
    pos_vidx_sort_buf[variant_ct] = ~0LLU;
    const uint32_t new_chr_ct = write_chr_info.chr_ct;
    if (chr_idxs) {
      uint32_t* next_write_vidxs;
      if (bigstack_alloc_ui(chr_code_end, &next_write_vidxs)) {
        goto make_plink2_vsort_ret_NOMEM;
      }
      for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx < new_chr_ct; ++new_chr_fo_idx) {
        const uint32_t chr_idx = write_chr_info.chr_file_order[new_chr_fo_idx];
        next_write_vidxs[chr_idx] = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx];
      }
      uint32_t variant_uidx = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        next_set_unsafe_ck(variant_include, &variant_uidx);
        const uint32_t chr_idx = chr_idxs[variant_uidx];
        const uint32_t write_vidx = next_write_vidxs[chr_idx];
        pos_vidx_sort_buf[write_vidx] = (((uint64_t)variant_bps[variant_uidx]) << 32) | variant_uidx;
        next_write_vidxs[chr_idx] += 1;
      }
      bigstack_reset(next_write_vidxs);
    } else {
      uint32_t old_chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t variant_uidx = 0;
      uint32_t chr_idx = 0;
      uint32_t write_vidx = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx, ++write_vidx) {
        next_set_unsafe_ck(variant_include, &variant_uidx);
        if (variant_uidx >= chr_end) {
          do {
            ++old_chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[old_chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          chr_idx = cip->chr_file_order[old_chr_fo_idx];
          write_vidx = write_chr_info.chr_idx_to_foidx[chr_idx];
        }
        pos_vidx_sort_buf[write_vidx] = (((uint64_t)variant_bps[variant_uidx]) << 32) | variant_uidx;
      }
    }

    str_sort_indexed_deref_t* same_pos_sort_buf = (str_sort_indexed_deref_t*)g_bigstack_base;
    const uintptr_t same_pos_sort_buf_size = bigstack_left() / sizeof(str_sort_indexed_deref_t);

    uint32_t vidx_start = 0;
    uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
    for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx < new_chr_ct; ++new_chr_fo_idx) {
      const uint32_t vidx_end = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx + 1];
      const uint32_t chr_size = vidx_end - vidx_start;
      const uint64_t post_entry = pos_vidx_sort_buf[vidx_end];
      pos_vidx_sort_buf[vidx_end] = ~0LLU; // simplify end-of-chromosome logic
      uint64_t* pos_vidx_sort_chr = &(pos_vidx_sort_buf[vidx_start]);

#ifdef __cplusplus
      std::sort(pos_vidx_sort_chr, &(pos_vidx_sort_chr[chr_size]));
#else
      qsort(pos_vidx_sort_chr, chr_size, sizeof(int64_t), uint64cmp);
#endif
      uint32_t prev_pos = pos_vidx_sort_chr[0] >> 32;
      uint32_t prev_variant_uidx = (uint32_t)pos_vidx_sort_chr[0];
      uint32_t prev_cidx = 0;
      uint32_t cidx = 1;
      for (; cidx < chr_size; ++cidx) {
        uint64_t cur_entry = pos_vidx_sort_chr[cidx];
        uint32_t cur_pos = cur_entry >> 32;
        if (cur_pos == prev_pos) {
          same_pos_sort_buf[0].strptr = variant_ids[prev_variant_uidx];
          same_pos_sort_buf[0].orig_idx = prev_variant_uidx;
          uint32_t equal_pos_ct = 1;
          const uint64_t* pos_vidx_sort_chr2 = &(pos_vidx_sort_chr[prev_cidx]);
          do {
            if (equal_pos_ct >= same_pos_sort_buf_size) {
              goto make_plink2_vsort_ret_NOMEM;
            }
            const uint32_t variant_uidx = (uint32_t)cur_entry;
            same_pos_sort_buf[equal_pos_ct].strptr = variant_ids[variant_uidx];
            same_pos_sort_buf[equal_pos_ct].orig_idx = variant_uidx;
            cur_entry = pos_vidx_sort_chr2[++equal_pos_ct];
            cur_pos = cur_entry >> 32;
          } while (cur_pos == prev_pos);
          strptr_arr_sort_main(equal_pos_ct, use_nsort, same_pos_sort_buf);
          for (uint32_t equal_pos_idx = 0; equal_pos_idx < equal_pos_ct; ++equal_pos_idx) {
            *new_variant_idx_to_old_iter++ = same_pos_sort_buf[equal_pos_idx].orig_idx;
          }
          cidx += equal_pos_ct - 1;
        } else {
          *new_variant_idx_to_old_iter++ = prev_variant_uidx;
        }
        prev_pos = cur_pos;
        prev_cidx = cidx;
        prev_variant_uidx = (uint32_t)cur_entry;
      }
      if (cidx == chr_size) {
        // if [cidx - 1] is part of an identical-bp batch, cidx will actually
        // be chr_size + 1 after loop exit.  It's equal to chr_size iff we
        // haven't written the last entry to new_variant_idx_to_old[].
        *new_variant_idx_to_old_iter++ = prev_variant_uidx;
      }
      vidx_start = vidx_end;
      pos_vidx_sort_buf[vidx_end] = post_entry;
    }
    bigstack_reset(pos_vidx_sort_buf);

    const uint32_t trim_alts = (uint32_t)(make_plink2_modifier & kfMakePlink2TrimAlts);
    if (make_plink2_modifier & kfMakeBim) {
      char* bimname_end = strcpya0(outname_end, ".bim");
      const uint32_t bim_zst = (make_plink2_modifier / kfMakeBimZs) & 1;
      if (bim_zst) {
        strcpy(bimname_end, ".zst");
      }
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);

      reterr = write_bim_resorted(outname, &write_chr_info, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, variant_cms, new_variant_idx_to_old, variant_ct, max_allele_slen, bim_zst, max_thread_ct);
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
      if (!simple_pgrp->fi.nonref_flags) {
        nonref_flags_storage = (simple_pgrp->fi.gflags & kfPgenGlobalAllNonref)? 2 : 1;
      }
      reterr = write_pvar_resorted(outname, xheader, variant_include, &write_chr_info, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, simple_pgrp->fi.nonref_flags, pvar_info_reload, variant_cms, new_variant_idx_to_old, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, xheader_info_pr, xheader_info_pr_nonflag, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_modifier, max_thread_ct);
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
    if (make_plink2_modifier & (kfMakeBed | kfMakePgen)) {
      // boilerplate from start of make_plink2_no_vsort()
      if (make_plink2_modifier & kfMakePlink2MMask) {
        logerrprint("Error: --make-bed/--make-{b}pgen multiallelics= is currently under development.\n");
        reterr = kPglRetNotYetSupported;
        goto make_plink2_vsort_ret_1;
      }
      g_plink2_write_flags = kfPlink2Write0;
      const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
      if (make_plink2_modifier & kfMakePlink2SetHhMissing) {
        const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
        uintptr_t* sex_female;
        if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male) ||
            bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
            bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed) ||
            bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed_interleaved) ||
            bigstack_alloc_ul(raw_sample_ctl, &sex_female)) {
          goto make_plink2_vsort_ret_NOMEM;
        }
        copy_bitarr_subset(sex_male, sample_include, sample_ct, g_sex_male);
        zero_trailing_words(BITCT_TO_WORDCT(sample_ct), g_sex_male);
        fill_interleaved_mask_vec(g_sex_male, sample_ctv, g_sex_male_collapsed_interleaved);

        bitvec_andnot_copy(sex_nm, sex_male, raw_sample_ctl, sex_female);
        copy_bitarr_subset(sex_female, sample_include, sample_ct, g_sex_female_collapsed);
        fill_interleaved_mask_vec(g_sex_female_collapsed, sample_ctv, g_sex_female_collapsed_interleaved);

        bigstack_reset(g_sex_female_collapsed);
        g_plink2_write_flags |= kfPlink2WriteSetHhMissing;
      }
      if (make_plink2_modifier & kfMakePlink2SetMixedMtMissing) {
        g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
      }
      g_cip = &write_chr_info;
      reterr = make_pgen_robust(sample_include, new_sample_idx_to_old, variant_include, variant_allele_idxs, refalt1_select, new_variant_idx_to_old, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, hard_call_thresh, dosage_erase_thresh, make_plink2_modifier, simple_pgrp, outname, outname_end);
      if (reterr) {
        goto make_plink2_vsort_ret_1;
      }
    }
  }
  while (0) {
  make_plink2_vsort_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 make_plink2_vsort_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

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
    reterr = sorted_xidbox_init_alloc(sample_include, sample_ids, sids, sample_ct, max_sample_id_blen, max_sid_blen, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (reterr) {
      goto sample_sort_file_map_ret_1;
    }
    uint32_t* new_sample_idx_to_old_iter = *new_sample_idx_to_old_ptr;
    if (*loadbuf_first_token == '#') {
      *loadbuf_first_token = '\0';
    }
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

#ifdef __cplusplus
} // namespace plink2
#endif
