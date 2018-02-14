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

PglErr WriteMapOrBim(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst, uint32_t thread_ct) {
  // set max_allele_slen to zero for .map
  // allele_dosages must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
      goto WriteMapOrBim_ret_NOMEM;
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto WriteMapOrBim_ret_1;
    }

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t variant_uidx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      FindFirst1BitFromU32(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        char* chr_name_end = chrtoa(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end = delim;
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
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
        const char* const* cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
          const AltAlleleCt* cur_refalt1_select = &(refalt1_select[variant_uidx * 2]);
          if ((!allele_dosages) || allele_dosages[cur_refalt1_select[1] + variant_allele_idx_base]) {
            cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[1]]);
          } else {
            *cswritep++ = output_missing_geno_char;
          }
          *cswritep++ = delim;
          cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[0]]);
        }
      }
      AppendBinaryEoln(&cswritep);
      if (Cswrite(&css, &cswritep)) {
        goto WriteMapOrBim_ret_WRITE_FAIL;
      }
    }
    if (CswriteCloseNull(&css, cswritep)) {
      goto WriteMapOrBim_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WriteMapOrBim_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteMapOrBim_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 WriteMapOrBim_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr PvarInfoReloadHeader(uintptr_t loadbuf_size, gzFile gz_pvar_reload, char* loadbuf, uint32_t* info_col_idx_ptr) {
  const char* loadbuf_iter;
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
    loadbuf_iter = FirstNonTspace(loadbuf);
  } while (!StrStartsWithUnsafe(loadbuf_iter, "#CHROM"));
  uint32_t info_col_idx = 0;
  do {
    loadbuf_iter = NextToken(loadbuf_iter);
    ++info_col_idx;
  } while (!tokequal_k(loadbuf_iter, "INFO"));
  *info_col_idx_ptr = info_col_idx;
  return kPglRetSuccess;
}

PglErr PvarInfoOpenAndReloadHeader(const char* pvar_info_reload, gzFile* gz_pvar_reload_ptr, char** loadbuf_ptr, uintptr_t* loadbuf_size_ptr, uint32_t* info_col_idx_ptr) {
  PglErr reterr = gzopen_read_checked(pvar_info_reload, gz_pvar_reload_ptr);
  if (reterr) {
    return reterr;
  }
  uintptr_t loadbuf_size = bigstack_left();
  if (loadbuf_size > kMaxLongLine) {
    loadbuf_size = kMaxLongLine;
  } else if (loadbuf_size <= kMaxMediumLine) {
    return kPglRetNomem;
  }
  char* loadbuf = R_CAST(char*, g_bigstack_base);
  loadbuf[loadbuf_size - 1] = ' ';
  *loadbuf_ptr = loadbuf;
  *loadbuf_size_ptr = loadbuf_size;
  return PvarInfoReloadHeader(loadbuf_size, *gz_pvar_reload_ptr, loadbuf, info_col_idx_ptr);
}

void PvarInfoWrite(uint32_t xheader_info_pr, uint32_t is_pr, char* info_token, char** write_iter_ptr) {
  char* info_token_end = CurTokenEnd(info_token);
  uint32_t info_token_slen = info_token_end - info_token;
  char* info_token_pr = nullptr;
  if (xheader_info_pr) {
    info_token_pr = PrInInfoToken(info_token_slen, info_token);
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
      write_iter = memcpya(write_iter, info_token, S_CAST(uintptr_t, info_token_pr - info_token) - 1);
      const char* pr_end = &(info_token_pr[2]);
      write_iter = memcpya(write_iter, pr_end, info_token_end - pr_end);
    }
  }
  *write_iter_ptr = write_iter;
}

PglErr PvarInfoReloadAndWrite(uintptr_t loadbuf_size, uint32_t xheader_info_pr, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, gzFile gz_pvar_reload, char** write_iter_ptr, uint32_t* gz_variant_uidx_ptr, char* loadbuf) {
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
      loadbuf_first_token = FirstNonTspace(loadbuf);
    } while (IsEolnKns(*loadbuf_first_token));
    ++gz_variant_uidx;
  } while (gz_variant_uidx <= variant_uidx);
  char* info_token = NextTokenMult(loadbuf_first_token, info_col_idx);
  *gz_variant_uidx_ptr = gz_variant_uidx;
  PvarInfoWrite(xheader_info_pr, is_pr, info_token, write_iter_ptr);
  return kPglRetSuccess;
}

void AppendChrsetLine(const ChrInfo* cip, char** write_iter_ptr) {
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
  AppendBinaryEoln(write_iter_ptr);
}

PglErr WritePvar(const char* outname, const char* xheader, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, PvarPsamFlags pvar_psam_flags, uint32_t thread_ct) {
  // allele_dosages must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  gzFile gz_pvar_reload = nullptr;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;

    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
        bigstack_alloc_ul(BitCtToWordCt(kPglMaxAltAlleleCt), &allele_include)) {
      goto WritePvar_ret_NOMEM;
    }
    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_flags / kfPvarZs) & 1;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto WritePvar_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_flags & kfPvarColInfo) || pvar_info_reload;
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
      logputs("\n");
      logerrputs("Error: Conflicting INFO:PR fields.  Either fix all REF alleles so that the\n'provisional reference' field is no longer needed, or remove/rename the other\nINFO:PR field.\n");
      goto WritePvar_ret_INCONSISTENT_INPUT;
    }

    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t info_col_idx = 0;  // could save this during first load instead
    if (pvar_psam_flags & kfPvarColXheader) {
      if (CsputsStd(xheader, xheader_blen, &css, &cswritep)) {
        goto WritePvar_ret_WRITE_FAIL;
      }
      if (write_info_pr && (!xheader_info_pr)) {
        cswritep = strcpya(cswritep, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
      }
    }
    // bugfix (30 Jul 2017): may be necessary to reload INFO when no ## lines
    // are in the header
    if (pvar_info_reload) {
      reterr = PvarInfoOpenAndReloadHeader(pvar_info_reload, &gz_pvar_reload, &loadbuf, &loadbuf_size, &info_col_idx);
      if (reterr) {
        goto WritePvar_ret_1;
      }
    }
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &cswritep);
    }
    cswritep = strcpya(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_flags & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybequal) && qual_present) {
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
    if (pvar_psam_flags & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybefilter) && filter_present) {
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
    if (pvar_psam_flags & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
        // nonzero_cm_present check was performed
        write_cm = 1;
      } else {
        uint32_t variant_uidx = 0;
        for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
          FindFirst1BitFromU32(variant_include, &variant_uidx);
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
    AppendBinaryEoln(&cswritep);

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
      FindFirst1BitFromU32(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        char* chr_name_end = chrtoa(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
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
      const char* const* cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
      if (Cswrite(&css, &cswritep)) {
        goto WritePvar_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
        SetAllBits(cur_allele_ct, allele_include);
        ClearBit(ref_allele_idx, allele_include);
        ClearBit(alt1_allele_idx, allele_include);
        uint32_t cur_allele_uidx = 0;
        uint32_t alt_allele_idx = 2;
        do {
          *cswritep++ = ',';
          FindFirst1BitFromU32(allele_include, &cur_allele_uidx);
          cswritep = strcpya(cswritep, cur_alleles[cur_allele_uidx++]);
          if (Cswrite(&css, &cswritep)) {
            goto WritePvar_ret_WRITE_FAIL;
          }
        } while (++alt_allele_idx < cur_allele_ct);
      }

      if (write_qual) {
        *cswritep++ = '\t';
        if ((!qual_present) || (!IsSet(qual_present, variant_uidx))) {
          *cswritep++ = '.';
        } else {
          cswritep = ftoa_g(quals[variant_uidx], cswritep);
        }
      }

      if (write_filter) {
        *cswritep++ = '\t';
        if ((!filter_present) || (!IsSet(filter_present, variant_uidx))) {
          *cswritep++ = '.';
        } else if (!IsSet(filter_npass, variant_uidx)) {
          cswritep = strcpya(cswritep, "PASS");
        } else {
          cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
        }
      }

      if (write_info) {
        *cswritep++ = '\t';
        const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
        if (gz_pvar_reload) {
          reterr = PvarInfoReloadAndWrite(loadbuf_size, xheader_info_pr, info_col_idx, variant_uidx, is_pr, gz_pvar_reload, &cswritep, &gz_variant_uidx, loadbuf);
          if (reterr) {
            goto WritePvar_ret_1;
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
      AppendBinaryEoln(&cswritep);
    }
    if (CswriteCloseNull(&css, cswritep)) {
      goto WritePvar_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WritePvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WritePvar_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WritePvar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WritePvar_ret_1:
  CswriteCloseCond(&css, cswritep);
  gzclose_cond(gz_pvar_reload);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteFam(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, char delim) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto WriteFam_ret_OPEN_FAIL;
    }
    uintptr_t* pheno_nm = nullptr;
    uintptr_t* pheno_cc = nullptr;
    double* pheno_qt = nullptr;
    // .fam files don't support categorical phenotypes
    const uint32_t pheno_idx = FirstCcOrQtPhenoIdx(pheno_cols, pheno_ct);
    if (pheno_idx != UINT32_MAX) {
      const PhenoDtype type_code = pheno_cols[pheno_idx].type_code;
      pheno_nm = pheno_cols[pheno_idx].nonmiss;
      if (type_code == kPhenoDtypeCc) {
        pheno_cc = pheno_cols[pheno_idx].data.cc;
      } else {
        pheno_qt = pheno_cols[pheno_idx].data.qt;
      }
    }
    const char* legacy_output_missing_pheno = g_legacy_output_missing_pheno;
    const uint32_t lomp_slen = strlen(legacy_output_missing_pheno);

    // possible todo: warning if two sample IDs only differ in SID?  (check for
    // this if any file is being exported that can't have a SID column)
    const char* sample_ids = piip->sii.sample_ids;
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    uintptr_t sample_uidx = 0;
    uint32_t sample_uidx2 = 0;
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      if (!new_sample_idx_to_old) {
        FindFirst1BitFromL(sample_include, &sample_uidx);
      } else {
        do {
          sample_uidx = new_sample_idx_to_old[sample_uidx2++];
        } while (!IsSet(sample_include, sample_uidx));
      }
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      if (delim == '\t') {
        write_iter = strcpya(write_iter, cur_sample_id);
      } else {
        const char* fid_end = S_CAST(const char*, rawmemchr(cur_sample_id, '\t'));
        write_iter = memcpyax(write_iter, cur_sample_id, fid_end - cur_sample_id, delim);
        write_iter = strcpya(write_iter, &(fid_end[1]));
      }
      *write_iter++ = delim;
      write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), delim);
      write_iter = strcpyax(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]), delim);
      *write_iter++ = Sexchar(sex_nm, sex_male, sample_uidx);
      *write_iter++ = delim;
      if ((!pheno_nm) || (!IsSet(pheno_nm, sample_uidx))) {
        write_iter = memcpya(write_iter, legacy_output_missing_pheno, lomp_slen);
      } else if (pheno_cc) {
        // do we want to allow user to force 0/1 output?
        *write_iter++ = '1' + IsSet(pheno_cc, sample_uidx);
      } else {
        write_iter = dtoa_g(pheno_qt[sample_uidx], write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
        goto WriteFam_ret_WRITE_FAIL;
      }
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
      goto WriteFam_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WriteFam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  WriteFam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

uint32_t IsDataFidColRequired(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t maybe_modifier) {
  if (maybe_modifier & 2) {
    return 1;
  }
  if ((!(maybe_modifier & 1)) || (!(siip->flags & kfSampleIdFidPresent))) {
    return 0;
  }
  const char* sample_ids = siip->sample_ids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    FindFirst1BitFromU32(sample_include, &sample_uidx);
    if (memcmp(&(sample_ids[sample_uidx * max_sample_id_blen]), "0\t", 2)) {
      return 1;
    }
  }
  return 0;
}

uint32_t IsDataSidColRequired(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier) {
  // note that MAYBESID and SID can both be set
  if (maybe_modifier & 2) {
    return 1;
  }
  if (sids && (maybe_modifier & 1)) {
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      FindFirst1BitFromU32(sample_include, &sample_uidx);
      if (memcmp(&(sids[sample_uidx * max_sid_blen]), "0", 2)) {
        return 1;
      }
    }
  }
  return 0;
}

uint32_t IsParentalInfoPresent(const uintptr_t* sample_include, const ParentalIdInfo* parental_id_infop, uint32_t sample_ct) {
  const char* paternal_ids = parental_id_infop->paternal_ids;
  const char* maternal_ids = parental_id_infop->maternal_ids;
  const uintptr_t max_paternal_id_blen = parental_id_infop->max_paternal_id_blen;
  const uintptr_t max_maternal_id_blen = parental_id_infop->max_maternal_id_blen;
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    FindFirst1BitFromU32(sample_include, &sample_uidx);
    if ((!strequal_k_unsafe(&(paternal_ids[sample_uidx * max_paternal_id_blen]), "0")) || (!strequal_k_unsafe(&(maternal_ids[sample_uidx * max_maternal_id_blen]), "0"))) {
      return 1;
    }
  }
  return 0;
}

char* AppendPhenoStr(const PhenoCol* pheno_col, const char* output_missing_pheno, uint32_t omp_slen, uint32_t sample_uidx, char* write_iter) {
  const PhenoDtype type_code = pheno_col->type_code;
  if (type_code <= kPhenoDtypeQt) {
    if (!IsSet(pheno_col->nonmiss, sample_uidx)) {
      write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
    } else if (type_code == kPhenoDtypeCc) {
      *write_iter++ = '1' + IsSet(pheno_col->data.cc, sample_uidx);
    } else {
      write_iter = dtoa_g(pheno_col->data.qt[sample_uidx], write_iter);
    }
  } else {
    write_iter = strcpya(write_iter, pheno_col->category_names[pheno_col->data.cat[sample_uidx]]);
  }
  return write_iter;
}

PglErr WritePsam(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, PvarPsamFlags pvar_psam_flags) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto WritePsam_ret_OPEN_FAIL;
    }
    const char* output_missing_pheno = g_output_missing_pheno;
    const uint32_t omp_slen = strlen(output_missing_pheno);

    char* textbuf_flush = &(g_textbuf[kMaxMediumLine]);

    const char* sample_ids = piip->sii.sample_ids;
    const char* sids = piip->sii.sids;
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uintptr_t max_sid_blen = piip->sii.max_sid_blen;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    const uint32_t write_fid = IsDataFidColRequired(sample_include, &(piip->sii), sample_ct, pvar_psam_flags / kfPsamColMaybefid);
    const uint32_t write_sid = IsDataSidColRequired(sample_include, sids, sample_ct, max_sid_blen, pvar_psam_flags / kfPsamColMaybesid);
    uint32_t write_parents = 0;
    if (pvar_psam_flags & kfPsamColParents) {
      write_parents = 1;
    } else if (pvar_psam_flags & kfPsamColMaybeparents) {
      write_parents = IsParentalInfoPresent(sample_include, &(piip->parental_id_info), sample_ct);
    }
    const uint32_t write_sex = (pvar_psam_flags / kfPsamColSex) & 1;
    const uint32_t write_empty_pheno = (pvar_psam_flags & kfPsamColPheno1) && (!pheno_ct);
    const uint32_t write_phenos = (pvar_psam_flags & (kfPsamColPheno1 | kfPsamColPhenos)) && pheno_ct;
    if (write_phenos && (!(pvar_psam_flags & kfPsamColPhenos))) {
      pheno_ct = 1;
    }
    char* write_iter = g_textbuf;
    *write_iter++ = '#';
    if (write_fid) {
      write_iter = strcpya(write_iter, "FID\t");
    }
    write_iter = memcpyl3a(write_iter, "IID");
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
            logerrputs("Error: .psam file cannot have both a regular SEX column and a phenotype named\n'SEX'.  Exclude or rename one of these columns.\n");
            goto WritePsam_ret_INCONSISTENT_INPUT;
          }
          // does this phenotype column conform to the SEX column format?
          // case/control is always ok, but quantitative or categorical needs
          // to be checked
          const PhenoCol* sex_col = &(pheno_cols[pheno_idx]);
          if (sex_col->type_code != kPhenoDtypeCc) {
            // could bitwise-and sample_include and pheno_nm before the loop
            const uintptr_t* pheno_nm = sex_col->nonmiss;
            uint32_t sample_uidx = 0;
            if (sex_col->type_code == kPhenoDtypeQt) {
              const double* pheno_vals = sex_col->data.qt;
              for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
                FindFirst1BitFromU32(sample_include, &sample_uidx);
                if (IsSet(pheno_nm, sample_uidx)) {
                  const double dxx = pheno_vals[sample_uidx];
                  // tolerate '-9' and '0' as missing values, and anything in
                  // [1, 2] (could be reasonable to represent XXY, etc. with
                  // decimals).
                  if (((dxx < 1.0) && (dxx != -9.0) && (dxx != 0.0)) || (dxx > 2.0)) {
                    logerrputs("Error: .psam SEX values are expected to be in {-9, 0, 1, 2}.\n");
                    goto WritePsam_ret_INCONSISTENT_INPUT;
                  }
                }
              }
            } else {
              assert(sex_col->type_code == kPhenoDtypeCat);
              const uint32_t nonnull_cat_ct = sex_col->nonnull_category_ct;
              if (nonnull_cat_ct) {
                const char* const* cur_category_names = sex_col->category_names;
                // tolerate 'M' and 'm' being present simultaneously, etc.
                uint32_t male_cat_idx1 = 0;
                uint32_t male_cat_idx2 = 0;
                uint32_t female_cat_idx1 = 0;
                uint32_t female_cat_idx2 = 0;
                for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
                  const char* cur_cat_name = cur_category_names[cat_idx];
                  if (!cur_cat_name[1]) {
                    uint32_t first_char_code = ctou32(cur_cat_name[0]);
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
                if (S_CAST(uint32_t, (male_cat_idx1 != 0) + (male_cat_idx2 != 0) + (female_cat_idx1 != 0) + (female_cat_idx2 != 0)) < nonnull_cat_ct) {
                  const uint32_t* pheno_vals = sex_col->data.cat;
                  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
                    FindFirst1BitFromU32(sample_include, &sample_uidx);
                    if (IsSet(pheno_nm, sample_uidx)) {
                      const uint32_t cur_cat_idx = pheno_vals[sample_uidx];
                      if ((cur_cat_idx != male_cat_idx1) && (cur_cat_idx != female_cat_idx1) && (cur_cat_idx != male_cat_idx2) && (cur_cat_idx != female_cat_idx2)) {
                        logerrputs("Error: .psam SEX values are expected to be in {'F', 'f', 'M', 'm'}.\n");
                        goto WritePsam_ret_INCONSISTENT_INPUT;
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
          goto WritePsam_ret_WRITE_FAIL;
        }
      }
    } else if (write_empty_pheno) {
      write_iter = strcpya(write_iter, "\tPHENO1");
    }
    AppendBinaryEoln(&write_iter);

    uintptr_t sample_uidx = 0;
    uint32_t sample_uidx2 = 0;
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      if (!new_sample_idx_to_old) {
        FindFirst1BitFromL(sample_include, &sample_uidx);
      } else {
        do {
          sample_uidx = new_sample_idx_to_old[sample_uidx2++];
        } while (!IsSet(sample_include, sample_uidx));
      }
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      if (!write_fid) {
        cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
      }
      write_iter = strcpya(write_iter, cur_sample_id);
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
        if (IsSet(sex_nm, sample_uidx)) {
          *write_iter++ = '2' - IsSet(sex_male, sample_uidx);
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
          write_iter = AppendPhenoStr(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
          if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
            goto WritePsam_ret_WRITE_FAIL;
          }
        }
      } else {
        if (write_empty_pheno) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
        }
        if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
          goto WritePsam_ret_WRITE_FAIL;
        }
      }
      AppendBinaryEoln(&write_iter);
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
      goto WritePsam_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WritePsam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  WritePsam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WritePsam_ret_INCONSISTENT_INPUT:
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
static PglErr g_error_ret = kPglRetSuccess;

/*
#ifdef __arm__
#  error "Unaligned accesses in bitvec_resort()."
#endif
void bitvec_resort(const uintptr_t* bitvec, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, unsigned char* writebuf) {
  const uint32_t sample_ctl_m1 = BitCtToWordCt(sample_ct) - 1;
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
      cur_word |= IsSet(bitvec, new_sample_idx_to_old_base[uii]) << uii;
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
#  error "Unaligned accesses in GenovecResort()."
#endif
void GenovecResort(const uintptr_t* genovec, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, void* writebuf) {
  // writebuf need not be word-aligned
  const uint32_t sample_ctl2_m1 = QuaterCtToWordCt(sample_ct) - 1;
  uint32_t word_idx = 0;
  uint32_t cur_word_entry_ct = kBitsPerWordD2;
  const uint32_t* new_sample_idx_to_old_iter = new_sample_idx_to_old;
  uintptr_t* writebuf_walias = S_CAST(uintptr_t*, writebuf);
  while (1) {
    if (word_idx == sample_ctl2_m1) {
      cur_word_entry_ct = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t cur_word = 0;
    for (uint32_t uii = 0; uii < cur_word_entry_ct; ++uii) {
      cur_word |= (GetQuaterarrEntry(genovec, new_sample_idx_to_old_iter[uii])) << (2 * uii);
    }
    if (word_idx == sample_ctl2_m1) {
      SubwordStore(cur_word, QuaterCtToByteCt(cur_word_entry_ct), &(writebuf_walias[word_idx]));
      return;
    }
    writebuf_walias[word_idx++] = cur_word;
    new_sample_idx_to_old_iter = &(new_sample_idx_to_old_iter[kBitsPerWordD2]);
  }
}

void UnpackHphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, uint32_t raw_sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
  if (!(phaseraw[0] & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    ExpandBytearr(phaseraw, all_hets, raw_sample_ctl, het_ct, 1, phaseinfo);
  } else {
    ExpandBytearrNested(&(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]), phaseraw, all_hets, raw_sample_ctl, het_ct, 1, *phasepresent_ptr, phaseinfo);
  }
}

void UnpackHphaseSubset(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
  if (!(phaseraw[0] & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    ExpandThenSubsetBytearr(phaseraw, all_hets, sample_include, het_ct, sample_ct, 1, phaseinfo);
  } else {
    const uint32_t first_half_word_ct = 1 + (het_ct / kBitsPerWord);
    const uint32_t raw_phasepresent_ct = PopcountWords(phaseraw, first_half_word_ct) - 1;
    // see "if (explicit_phasepresent) {}" block in PgrReadRaw().  Could
    // change this convention.
    ExpandThenSubsetBytearrNested(&(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]), phaseraw, all_hets, sample_include, sample_ct, raw_phasepresent_ct, 1, *phasepresent_ptr, phaseinfo);
  }
}

void UnpackAndResortHphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* sample_include, const uint32_t* old_sample_idx_to_new, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uintptr_t* phaseraw_iter = phaseraw;
  const uint32_t* old_sample_idx_to_new_iter = old_sample_idx_to_new;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uintptr_t phaseraw_word = *phaseraw_iter++;
  uint32_t read_idx_lowbits = 1;
  ZeroUlArr(sample_ctl, phaseinfo);
  if (!(phaseraw_word & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
      uintptr_t new_phasepresent_word = all_hets[widx];
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + PopcountWord(new_phasepresent_word);
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
            SetBit(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
          }
          tmp_phaseinfo_input_word >>= 1;
          new_phasepresent_word &= new_phasepresent_word - 1;
        }
      } else {
        uintptr_t masked_phasepresent_word = new_phasepresent_word & sample_include[widx];
        while (masked_phasepresent_word) {
          const uint32_t sample_uidx_lowbits = CTZLU(masked_phasepresent_word);
          const uintptr_t lowmask = (k1LU << sample_uidx_lowbits) - k1LU;
          if ((tmp_phaseinfo_input_word >> PopcountWord(new_phasepresent_word & lowmask)) & 1) {
            SetBit(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
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
  ZeroUlArr(sample_ctl, phasepresent);
  for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
    uintptr_t geno_hets = all_hets[widx];
    if (geno_hets) {
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + PopcountWord(geno_hets);
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
        const uint32_t read_phasepresent_ct = PopcountWord(tmp_phasepresent_input_word);
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
              SetBit(new_sample_idx, phasepresent);
              if (tmp_phaseinfo_input_word & 1) {
                SetBit(new_sample_idx, phaseinfo);
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
                SetBit(new_sample_idx, phasepresent);
                if (tmp_phaseinfo_input_word & 1) {
                  SetBit(new_sample_idx, phaseinfo);
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

void CopyDosage(const uintptr_t* __restrict dosageraw, uint32_t raw_sample_ct, uintptr_t* __restrict write_dosagepresent, Dosage* write_dosagevals, uint32_t* write_dosage_ct_ptr) {
  const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
  const uintptr_t* read_dosagepresent = dosageraw;
  const Dosage* read_dosagevals = R_CAST(const Dosage*, &(dosageraw[raw_sample_ctaw]));
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t dosage_ct = PopcountWords(read_dosagepresent, raw_sample_ctl);
  *write_dosage_ct_ptr = dosage_ct;
  memcpy(write_dosagepresent, read_dosagepresent, raw_sample_ctl * sizeof(intptr_t));
  memcpy(write_dosagevals, read_dosagevals, dosage_ct * sizeof(Dosage));
}

void CopyDosageSubset(const uintptr_t* __restrict dosageraw, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* __restrict write_dosagepresent, Dosage* write_dosagevals, uint32_t* __restrict write_dosage_ct_ptr) {
  const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
  const uintptr_t* read_dosagepresent = dosageraw;
  const Dosage* read_dosagevals = R_CAST(const Dosage*, &(dosageraw[raw_sample_ctaw]));
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t read_dosage_ct = PopcountWords(read_dosagepresent, raw_sample_ctl);
  CopyBitarrSubset(read_dosagepresent, sample_include, sample_ct, write_dosagepresent);
  uint32_t sample_uidx = 0;
  Dosage* write_dosagevals_iter = write_dosagevals;
  for (uint32_t read_dosage_idx = 0; read_dosage_idx < read_dosage_ct; ++read_dosage_idx, ++sample_uidx) {
    FindFirst1BitFromU32(read_dosagepresent, &sample_uidx);
    if (IsSet(sample_include, sample_uidx)) {
      *write_dosagevals_iter++ = read_dosagevals[read_dosage_idx];
    }
  }
  *write_dosage_ct_ptr = write_dosagevals_iter - write_dosagevals;
}

void CopyAndResortDosage(const uintptr_t* __restrict dosageraw, const uint32_t* new_sample_idx_to_old, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* __restrict write_dosagepresent, Dosage* write_dosagevals, uint32_t* write_dosage_ct_ptr, uint32_t* cumulative_popcount_buf) {
  const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
  const uintptr_t* read_dosagepresent = dosageraw;
  const Dosage* read_dosagevals = R_CAST(const Dosage*, &(dosageraw[raw_sample_ctaw]));
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  FillCumulativePopcounts(read_dosagepresent, raw_sample_ctl, cumulative_popcount_buf);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  ZeroUlArr(sample_ctl, write_dosagepresent);
  Dosage* write_dosagevals_iter = write_dosagevals;
  for (uint32_t new_sample_idx = 0; new_sample_idx < sample_ct; ++new_sample_idx) {
    const uint32_t old_sample_idx = new_sample_idx_to_old[new_sample_idx];
    if (IsSet(read_dosagepresent, old_sample_idx)) {
      SetBit(new_sample_idx, write_dosagepresent);
      const uint32_t old_dosagevals_idx = RawToSubsettedPos(read_dosagepresent, cumulative_popcount_buf, old_sample_idx);
      *write_dosagevals_iter++ = read_dosagevals[old_dosagevals_idx];
    }
  }
  *write_dosage_ct_ptr = write_dosagevals_iter - write_dosagevals;
}


// more multithread globals
static PgenReader** g_pgr_ptrs = nullptr;
static uintptr_t** g_genovecs = nullptr;
static uintptr_t** g_dosage_presents = nullptr;
static Dosage** g_dosage_val_bufs = nullptr;
static uint32_t* g_read_variant_uidx_starts = nullptr;  // size calc_thread_ct

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
static const ChrInfo* g_cip = nullptr;
static const uintptr_t* g_sample_include = nullptr;
static uintptr_t* g_sample_include_interleaved_vec = nullptr;
static uint32_t* g_sample_include_cumulative_popcounts = nullptr;
static const uintptr_t* g_sex_male = nullptr;
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
static const AltAlleleCt* g_refalt1_select = nullptr;
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
static PgenGlobalFlags g_read_phase_dosage_gflags = kfPgenGlobal0;

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
static Dosage** g_thread_write_dosagevals = nullptr;
static uint32_t** g_thread_cumulative_popcount_bufs = nullptr;
static PgenWriterCommon** g_pwcs = nullptr;

THREAD_FUNC_DECL LoadAlleleAndGenoCountsThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  PgenReader* pgrp = g_pgr_ptrs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const ChrInfo* cip = g_cip;
  const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t subset_ct = (g_founder_info != nullptr) + 1;
  const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t first_hap_uidx = g_first_hap_uidx;
  const int32_t y_code = cip->xymt_codes[kChrOffsetY];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_vals = nullptr;
  if (g_dosage_presents) {
    dosage_present = g_dosage_presents[tidx];
    dosage_vals = g_dosage_val_bufs[tidx];
  }
  uint32_t is_y = 0;
  uint32_t is_nonxy_haploid = 0;
  uint32_t x_start = 0;
  int32_t x_code;
  if (XymtExists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[S_CAST(uint32_t, x_code)];
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
      PglErr reterr = kPglRetSuccess;

      // different approach will be needed for multiallelic case...
      uint32_t genocounts[4];
      uint32_t sex_specific_genocounts[4];
      for (; cur_idx < cur_idx_end; ++cur_idx, ++variant_uidx) {
        FindFirst1BitFromU32(variant_include, &variant_uidx);
        if (variant_uidx >= chr_end) {
          const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
          const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          is_y = 0;
          is_nonxy_haploid = 0;
          if (chr_idx == x_code) {
            is_x_or_y = 1;
            PgrClearLdCache(pgrp);
          } else if (chr_idx == y_code) {
            is_x_or_y = 1;
            is_y = 1;
            PgrClearLdCache(pgrp);
          } else {
            if (is_x_or_y) {
              PgrClearLdCache(pgrp);
            }
            is_x_or_y = 0;
            // no way for this to happen now unless everything is haploid?
            is_nonxy_haploid = IsSetI(cip->haploid_mask, chr_idx);
          }
        }
        const uintptr_t cur_variant_allele_idx = variant_allele_idxs? variant_allele_idxs[variant_uidx] : (2 * variant_uidx);
        // insert cur_allele_ct == 2 check here
        uint64_t cur_dosages[2];
        uint32_t hethap_ct;
        if (!is_x_or_y) {
          // call pgr_get_refalt1_genotype_counts() instead when dosages not
          // needed?
          reterr = PgrGetRefNonrefGenotypeCountsAndDosage16s(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, mach_r2_vals? (&(mach_r2_vals[variant_uidx])) : nullptr, genocounts, cur_dosages);
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
          reterr = PgrGetRefNonrefGenotypeCountsAndDosage16s(sex_male, sex_male_interleaved_vec, sex_male_cumulative_popcounts, male_ct, variant_uidx, pgrp, nullptr, genocounts, cur_dosages);
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
          reterr = PgrReadRefalt1GenovecDosage16SubsetUnsafe(nullptr, nullptr, raw_sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
          if (reterr) {
            g_error_ret = reterr;
            break;
          }
          // assert(!is_explicit_alt1);
          if (sample_ct == raw_sample_ct) {
            ZeroTrailingQuaters(raw_sample_ct, genovec);
            GenovecCountFreqsUnsafe(genovec, sample_ct, genocounts);
          } else {
            GenovecCountSubsetFreqs(genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
          }
          GenovecCountSubsetFreqs(genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, sex_specific_genocounts);
          hethap_ct = sex_specific_genocounts[1];
          if (allele_dosages) {
            uint32_t sample_uidx = 0;
            uintptr_t replaced_ct = 0;  // nonmales count twice
            uintptr_t alt1_ct = 4 * genocounts[2] + 2 * genocounts[1] - 2 * sex_specific_genocounts[2] - hethap_ct;  // nonmales count twice
            uint64_t alt1_dosage = 0;  // in 32768ths, nonmales count twice
            uint32_t included_dosage_ct = 0;  // nonmales count twice
            if (sample_ct == raw_sample_ct) {
              for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
                FindFirst1BitFromU32(dosage_present, &sample_uidx);
                const uintptr_t cur_dosage_val = dosage_vals[dosage_idx];
                const uintptr_t sex_multiplier = 2 - IsSet(sex_male, sample_uidx);
                alt1_dosage += cur_dosage_val * sex_multiplier;

                // could call GenoarrCountSubsetIntersectFreqs() twice
                // instead, but since we've already manually extracted the sex
                // bit it probably doesn't help?
                const uintptr_t hardcall_code = GetQuaterarrEntry(genovec, sample_uidx);
                if (hardcall_code != 3) {
                  alt1_ct -= hardcall_code * sex_multiplier;
                  replaced_ct += sex_multiplier;
                }
              }
              included_dosage_ct = 2 * dosage_ct;
              if (dosage_ct) {
                included_dosage_ct -= PopcountWordsIntersect(dosage_present, sex_male, raw_sample_ctl);
              }
            } else {
              for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
                FindFirst1BitFromU32(dosage_present, &sample_uidx);
                if (IsSet(sample_include, sample_uidx)) {
                  const uintptr_t cur_dosage_val = dosage_vals[dosage_idx];
                  const uintptr_t sex_multiplier = 2 - IsSet(sex_male, sample_uidx);
                  alt1_dosage += cur_dosage_val * sex_multiplier;
                  included_dosage_ct += sex_multiplier;
                  const uintptr_t hardcall_code = GetQuaterarrEntry(genovec, sample_uidx);
                  if (hardcall_code != 3) {
                    alt1_ct -= hardcall_code * sex_multiplier;
                    replaced_ct += sex_multiplier;
                  }
                }
              }
            }
            const uintptr_t ref_ct = (2 * (sample_ct - genocounts[3]) - male_ct + sex_specific_genocounts[3]) * (2 * k1LU) - alt1_ct;
            allele_dosages[cur_variant_allele_idx] = included_dosage_ct * S_CAST(uint64_t, kDosageMax) - alt1_dosage + (ref_ct * S_CAST(uint64_t, kDosageMid));
            allele_dosages[cur_variant_allele_idx + 1] = alt1_dosage + (alt1_ct * S_CAST(uint64_t, kDosageMid));
          }
          if (x_male_geno_cts) {
            uint32_t* cur_x_male_geno_cts = &(x_male_geno_cts[(3 * k1LU) * (variant_uidx - x_start)]);
            cur_x_male_geno_cts[0] = sex_specific_genocounts[0];
            cur_x_male_geno_cts[1] = sex_specific_genocounts[1];
            cur_x_male_geno_cts[2] = sex_specific_genocounts[2];
            if (x_nosex_geno_cts) {
              GenovecCountSubsetFreqs(genovec, nosex_interleaved_vec, raw_sample_ct, nosex_ct, sex_specific_genocounts);
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
              ZeroTrailingQuaters(raw_sample_ct, genovec);
              missing_dosage_ct = GenoarrCountMissingInvsubsetUnsafe(genovec, dosage_present, raw_sample_ct);
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
      PgrClearLdCache(pgrp);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* variant_allele_idxs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint64_t* allele_dosages, uint64_t* founder_allele_dosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, uint32_t* raw_geno_cts, uint32_t* founder_raw_geno_cts, uint32_t* x_male_geno_cts, uint32_t* founder_x_male_geno_cts, uint32_t* x_nosex_geno_cts, uint32_t* founder_x_nosex_geno_cts, double* mach_r2_vals) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    if (!variant_ct) {
      goto LoadAlleleAndGenoCounts_ret_1;
    }
    if (variant_allele_idxs) {
      logerrputs("Error: LoadAlleleAndGenoCounts() doesn't support multiallelic variants yet.\n");
      reterr = kPglRetNotYetSupported;
      goto LoadAlleleAndGenoCounts_ret_1;
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
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t raw_sample_ctv = BitCtToVecCt(raw_sample_ct);
    if (bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_sample_include_interleaved_vec) ||
        bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_cumulative_popcounts) ||
        bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_sex_male_interleaved_vec) ||
        bigstack_alloc_ui(raw_sample_ctl, &g_sex_male_cumulative_popcounts)) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    FillInterleavedMaskVec(g_sample_include, raw_sample_ctv, g_sample_include_interleaved_vec);
    FillCumulativePopcounts(g_sample_include, raw_sample_ctl, g_sample_include_cumulative_popcounts);
    if ((founder_ct == sample_ct) || (!only_founder_cts_required)) {
      g_sex_male = sex_male;
    } else {
      // no nonfounder counts required
      uintptr_t* new_sex_male;
      if (bigstack_alloc_ul(raw_sample_ctl, &new_sex_male)) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
        new_sex_male[widx] = sex_male[widx] & founder_info[widx];
      }
      ZeroTrailingWords(raw_sample_ctl, new_sex_male);
      g_sex_male = new_sex_male;
    }
    FillInterleavedMaskVec(g_sex_male, raw_sample_ctv, g_sex_male_interleaved_vec);
    FillCumulativePopcounts(g_sex_male, raw_sample_ctl, g_sex_male_cumulative_popcounts);
    if (!(x_nosex_geno_cts || founder_x_nosex_geno_cts)) {
      nosex_ct = 0;
    }
    g_nosex_ct = nosex_ct;
    g_nosex_interleaved_vec = nullptr;
    uintptr_t* nosex_buf = nullptr;
    if (nosex_ct) {
      if (bigstack_end_alloc_ul(raw_sample_ctl, &nosex_buf) ||
          bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_nosex_interleaved_vec)) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
        nosex_buf[widx] = (~sex_nm[widx]) & g_sample_include[widx];
      }
      ZeroTrailingWords(raw_sample_ctl, nosex_buf);
      FillInterleavedMaskVec(nosex_buf, raw_sample_ctv, g_nosex_interleaved_vec);
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
          goto LoadAlleleAndGenoCounts_ret_NOMEM;
        }
        FillInterleavedMaskVec(founder_info, raw_sample_ctv, g_founder_info_interleaved_vec);
        FillCumulativePopcounts(founder_info, raw_sample_ctl, g_founder_info_cumulative_popcounts);
        for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
          g_founder_male[widx] = sex_male[widx] & founder_info[widx];
        }
        ZeroTrailingWords(raw_sample_ctl, g_founder_male);
        FillInterleavedMaskVec(g_founder_male, raw_sample_ctv, g_founder_male_interleaved_vec);
        FillCumulativePopcounts(g_founder_male, raw_sample_ctl, g_founder_male_cumulative_popcounts);
        g_founder_ct = founder_ct;
        g_founder_male_ct = g_founder_male_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(g_founder_male[raw_sample_ctl - 1]);
        g_founder_allele_dosages = founder_allele_dosages;
        g_founder_raw_geno_cts = founder_raw_geno_cts;
        g_founder_x_male_geno_cts = founder_x_male_geno_cts;
        if (nosex_ct) {
          // caller currently responsible for ensuring that when
          // founder_nosex_ct is zero, founder_x_nosex_geno_cts ==
          // nullptr
          if (bigstack_alloc_ul(raw_sample_ctv * kWordsPerVec, &g_founder_nosex_interleaved_vec)) {
            goto LoadAlleleAndGenoCounts_ret_NOMEM;
          }
          for (uint32_t widx = 0; widx < raw_sample_ctl; ++widx) {
            nosex_buf[widx] &= founder_info[widx];
          }
          g_founder_nosex_ct = PopcountWords(nosex_buf, raw_sample_ctl);
          assert(g_founder_nosex_ct);
          ZeroTrailingWords(raw_sample_ctl, nosex_buf);
          FillInterleavedMaskVec(nosex_buf, raw_sample_ctv, g_founder_nosex_interleaved_vec);
          g_founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
        }
      } else {
        if (founder_allele_dosages) {
          ZeroU64Arr(variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct), founder_allele_dosages);
        }
        if (founder_raw_geno_cts) {
          ZeroUiArr((3 * k1LU) * raw_variant_ct, founder_raw_geno_cts);
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
      g_male_ct = g_sex_male_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(g_sex_male[raw_sample_ctl - 1]);
      if (nosex_ct) {
        g_nosex_ct = PopcountWords(nosex_buf, raw_sample_ctl);
      }
    }
    if (!g_sample_ct) {
      if (g_allele_dosages) {
        ZeroU64Arr(variant_allele_idxs? variant_allele_idxs[raw_variant_ct] : (2 * raw_variant_ct), g_allele_dosages);
      }
      if (g_raw_geno_cts) {
        ZeroUiArr((3 * k1LU) * raw_variant_ct, g_raw_geno_cts);
      }
      // early exit
      goto LoadAlleleAndGenoCounts_ret_1;
    }
    BigstackEndReset(bigstack_end_mark);  // free nosex_buf

    int32_t ii;
    const uint32_t x_dosages_needed = (allele_dosages || founder_allele_dosages || variant_missing_dosage_cts) && XymtExists(cip, kChrOffsetX, &ii) && (pgfip->gflags & kfPgenGlobalDosagePresent);
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
    if (PgenMtLoadInit(variant_include, raw_sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, x_dosages_needed? (&g_dosage_presents) : nullptr, x_dosages_needed? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }

    g_variant_include = variant_include;
    g_variant_allele_idxs = variant_allele_idxs;
    g_calc_thread_ct = calc_thread_ct;
    g_error_ret = kPglRetSuccess;

    logputs("Calculating allele frequencies... ");
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    const uint32_t read_block_sizel = BitCtToWordCt(read_block_size);
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
          cur_block_write_ct = PopcountWords(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
          if (cur_block_write_ct) {
            break;
          }
          ++read_block_idx;
        }
        if (read_block_idx == read_block_ct_m1) {
          cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
          cur_block_write_ct = PopcountWords(&(variant_include[read_block_idx * read_block_sizel]), BitCtToWordCt(cur_read_block_size));
        }
        if (PgfiMultiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
          if (variant_idx) {
            JoinThreads2z(calc_thread_ct, 0, threads);
            g_cur_block_write_ct = 0;
            ErrorCleanupThreads2z(LoadAlleleAndGenoCountsThread, calc_thread_ct, threads);
          }
          goto LoadAlleleAndGenoCounts_ret_READ_FAIL;
        }
      }
      if (variant_idx) {
        JoinThreads2z(calc_thread_ct, is_last_block, threads);
        reterr = g_error_ret;
        if (reterr) {
          if (!is_last_block) {
            g_cur_block_write_ct = 0;
            ErrorCleanupThreads2z(LoadAlleleAndGenoCountsThread, calc_thread_ct, threads);
          }
          if (reterr == kPglRetMalformedInput) {
            logputs("\n");
            logerrputs("Error: Malformed .pgen file.\n");
          }
          goto LoadAlleleAndGenoCounts_ret_1;
        }
      }
      if (!is_last_block) {
        g_cur_block_write_ct = cur_block_write_ct;
        ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
          g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
        }
        is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
        if (SpawnThreads2z(LoadAlleleAndGenoCountsThread, calc_thread_ct, is_last_block, threads)) {
          goto LoadAlleleAndGenoCounts_ret_THREAD_CREATE_FAIL;
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }

      ++read_block_idx;
      variant_idx += cur_block_write_ct;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logprintf("done.\n");
  }
  while (0) {
  LoadAlleleAndGenoCounts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadAlleleAndGenoCounts_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  LoadAlleleAndGenoCounts_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 LoadAlleleAndGenoCounts_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

void ApplyHardCallThresh(const uintptr_t* dosage_present, const Dosage* dosage_vals, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec) {
  uint32_t sample_uidx = 0;
  for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_uidx) {
    FindFirst1BitFromU32(dosage_present, &sample_uidx);
    const uint32_t dosage_int = dosage_vals[dosage_idx];
    const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
    uintptr_t new_geno;
    if (halfdist < hard_call_halfdist) {
      new_geno = 3;
    } else {
      new_geno = (dosage_int + kDosage4th) / kDosageMid;
    }
    AssignQuaterarrEntry(sample_uidx, new_geno, genovec);
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
FLAGSET_DEF_END(Plink2WriteFlags);
// todo: add .pgen-specific stuff

// more multithread globals
static Plink2WriteFlags g_plink2_write_flags = kfPlink2Write0;

THREAD_FUNC_DECL MakeBedlikeThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  PgenReader* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_vals = nullptr;
  uint32_t hard_call_halfdist = 0;
  if (g_dosage_presents) {
    dosage_present = g_dosage_presents[tidx];
    dosage_vals = g_dosage_val_bufs[tidx];
    hard_call_halfdist = g_hard_call_halfdist;
  }
  const uintptr_t* variant_include = g_variant_include;
  const ChrInfo* cip = g_cip;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed_interleaved = g_sex_female_collapsed_interleaved;
  const uint32_t* collapsed_sort_map = g_collapsed_sort_map;
  const uint32_t set_hh_missing = g_plink2_write_flags & kfPlink2WriteSetHhMissing;
  const uint32_t set_mixed_mt_missing = g_plink2_write_flags & kfPlink2WriteSetMixedMtMissing;
  const uint32_t write_plink1 = g_plink2_write_flags & kfPlink2WritePlink1;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  const uint32_t sample_ctv2 = QuaterCtToVecCt(sample_ct);
  const uint32_t sample_ct4 = QuaterCtToByteCt(sample_ct);
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const AltAlleleCt* refalt1_select = g_refalt1_select;
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
    uint32_t is_haploid_nonmt = 0;
    uint32_t is_mt = 0;
    for (; write_idx < write_idx_end; ++write_idx, ++variant_uidx) {
      FindFirst1BitFromU32(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        is_y = (chr_idx == y_code);
        is_x_or_y = is_y || (chr_idx == x_code);
        is_mt = (chr_idx == mt_code);
        is_haploid_nonmt = IsSet(cip->haploid_mask, chr_idx) && (!is_mt);
      }
      if (!hard_call_halfdist) {
        PglErr reterr = PgrReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec);
        if (reterr) {
          g_error_ret = reterr;
          break;
        }
      } else {
        // quasi-bugfix (4 Dec 2017): it's user-hostile to make
        // --hard-call-threshold not apply here.
        uint32_t dosage_ct;
        uint32_t is_explicit_alt1;
        PglErr reterr = PgrReadRefalt1GenovecDosage16SubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
        if (reterr) {
          g_error_ret = reterr;
          break;
        }
        ApplyHardCallThresh(dosage_present, dosage_vals, dosage_ct, hard_call_halfdist, genovec);
      }
      // this doesn't work in multiallelic case
      // todo: pgenlib_internal function which takes two allele indexes
      if (refalt1_select && (refalt1_select[2 * variant_uidx] == 1)) {
        GenovecInvertUnsafe(sample_ct, genovec);
      }
      if (set_hh_missing && is_haploid_nonmt) {
        if (is_x_or_y) {
          // male hets to missing
          SetMaleHetMissing(sex_male_collapsed_interleaved, sample_ctv2, genovec);
          if (is_y) {
            // all female calls to missing; unknown-sex calls now left alone
            InterleavedSetMissing(sex_female_collapsed_interleaved, sample_ctv2, genovec);
          }
        } else {
          // all hets to missing
          SetHetMissing(sample_ctl2, genovec);
        }
      } else if (set_mixed_mt_missing && is_mt) {
        // all hets to missing
        SetHetMissing(sample_ctl2, genovec);
      }
      // todo: --set-me-missing, --zero-cluster, --fill-missing-with-ref
      // (--set-me-missing should happen after --set-hh-missing)
      if (write_plink1) {
        PgrPlink2ToPlink1InplaceUnsafe(sample_ct, genovec);
      }
      // trailing bytes don't matter, but trailing bits of last byte may
      ZeroTrailingQuaters(sample_ct, genovec);
      if (!collapsed_sort_map) {
        writebuf_iter = R_CAST(unsigned char*, memcpya(writebuf_iter, genovec, sample_ct4));
      } else {
        GenovecResort(genovec, collapsed_sort_map, sample_ct, writebuf_iter);
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

static STPgenWriter* g_spgwp = nullptr;

THREAD_FUNC_DECL MakePgenThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t* new_sample_idx_to_old = g_new_sample_idx_to_old;
  const uint32_t* old_sample_idx_to_new = g_old_sample_idx_to_new;
  const ChrInfo* cip = g_cip;
  const uint32_t* write_chr_fo_vidx_start = g_write_chr_fo_vidx_start;
  const AltAlleleCt* refalt1_select_iter = g_refalt1_select;
  const uintptr_t* sample_include = g_sample_include;

  // er, this global should be named g_sex_male_collapsed...
  const uintptr_t* sex_male_collapsed = g_sex_male;

  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed = g_sex_female_collapsed;
  const uintptr_t* sex_female_collapsed_interleaved = g_sex_female_collapsed_interleaved;
  const uint32_t raw_sample_ct = g_raw_sample_ct;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  const uint32_t sample_ctv2 = QuaterCtToVecCt(sample_ct);
  const uint32_t raw_sample_ctaw2 = QuaterCtToAlignedWordCt(raw_sample_ct);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
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
  const uintptr_t phaseraw_word_ct = kWordsPerVec + RoundDownPow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
  // todo: double dosage_vals allocation in phased-dosage case
  const uintptr_t dosageraw_word_ct = kWordsPerVec * (BitCtToVecCt(raw_sample_ct) + DivUp(raw_sample_ct, (kBytesPerVec / sizeof(Dosage))));

  STPgenWriter* spgwp = g_spgwp;
  PgenWriterCommon* pwcp;
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
  Dosage* write_dosagevals = nullptr;
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
    uint32_t is_haploid_nonmt = 0;
    uint32_t is_mt = 0;
    for (; write_idx < write_idx_end; ++write_idx) {
      if (loaded_vrtypes) {
        loaded_vrtype = loaded_vrtypes[write_idx];
      }
      if (write_idx >= chr_end_bidx) {
        const uint32_t chr_fo_idx = CountSortedSmallerUi(&(write_chr_fo_vidx_start[1]), cip->chr_ct, write_idx + variant_idx_offset + 1);
        const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end_bidx = write_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
        is_y = (chr_idx == y_code);
        is_x_or_y = is_y || (chr_idx == x_code);
        is_mt = (chr_idx == mt_code);
        is_haploid_nonmt = IsSetI(cip->haploid_mask, chr_idx) && (!is_mt);
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
          PgrDetectGenovecHets(loadbuf_iter, raw_sample_ct, all_hets);
          cur_phaseraw = cur_genovec_end;
          cur_genovec_end = &(cur_genovec_end[phaseraw_word_ct]);
        }
        if (is_dosage) {
          cur_dosageraw = cur_genovec_end;
          cur_genovec_end = &(cur_genovec_end[dosageraw_word_ct]);
        }
        uint32_t write_dosage_ct = 0;
        if (new_sample_idx_to_old) {
          GenovecResort(loadbuf_iter, new_sample_idx_to_old, sample_ct, write_genovec);
          if (is_hphase) {
            UnpackAndResortHphase(all_hets, cur_phaseraw, sample_include, old_sample_idx_to_new, raw_sample_ct, sample_ct, &cur_write_phasepresent, write_phaseinfo);
          }
          if (is_dosage) {
            CopyAndResortDosage(cur_dosageraw, new_sample_idx_to_old, raw_sample_ct, sample_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct, cumulative_popcount_buf);
          }
        } else if (sample_include) {
          CopyQuaterarrNonemptySubset(loadbuf_iter, sample_include, raw_sample_ct, sample_ct, write_genovec);
          if (is_hphase) {
            UnpackHphaseSubset(all_hets, cur_phaseraw, sample_include, raw_sample_ct, sample_ct, &cur_write_phasepresent, write_phaseinfo);
          }
          if (is_dosage) {
            CopyDosageSubset(cur_dosageraw, sample_include, raw_sample_ct, sample_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct);
          }
        } else {
          write_genovec = loadbuf_iter;
          if (is_hphase) {
            UnpackHphase(all_hets, cur_phaseraw, sample_ct, &cur_write_phasepresent, write_phaseinfo);
          }
          if (is_dosage) {
            CopyDosage(cur_dosageraw, sample_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct);
          }
        }
        if (refalt1_select_iter && (refalt1_select_iter[2 * write_idx] == 1)) {
          GenovecInvertUnsafe(sample_ct, write_genovec);
          if (is_hphase) {
            BitarrInvert(sample_ctl, write_phaseinfo);
          }
          if (write_dosage_ct) {
            BiallelicDosage16Invert(write_dosage_ct, write_dosagevals);
          }
        }
        if (write_dosage_ct) {
          if (hard_call_halfdist) {
            if (is_hphase && (!cur_write_phasepresent)) {
              cur_write_phasepresent = write_phasepresent;
              // unsafe to just copy all_hets, because we may have resorted
              PgrDetectGenovecHets(write_genovec, sample_ct, cur_write_phasepresent);
            }
            uint32_t sample_uidx = 0;
            for (uint32_t dosage_idx = 0; dosage_idx < write_dosage_ct; ++dosage_idx, ++sample_uidx) {
              FindFirst1BitFromU32(write_dosagepresent, &sample_uidx);
              const uint32_t dosage_int = write_dosagevals[dosage_idx];
              const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
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
                  ClearBit(sample_uidx, cur_write_phasepresent);
                }
                write_genovec[widx] = prev_geno_word ^ (geno_xor << shift);
              }
            }
          }
          if (dosage_erase_halfdist < kDosage4th) {
            uint32_t dosage_read_idx = 0;
            uint32_t sample_uidx = 0;
            for (; dosage_read_idx < write_dosage_ct; ++dosage_read_idx, ++sample_uidx) {
              FindFirst1BitFromU32(write_dosagepresent, &sample_uidx);
              const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
              const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
              if (halfdist >= dosage_erase_halfdist) {
                ClearBit(sample_uidx, write_dosagepresent);
                ++sample_uidx;
                break;
              }
            }
            uint32_t dosage_write_idx = dosage_read_idx;
            while (++dosage_read_idx < write_dosage_ct) {
              FindFirst1BitFromU32(write_dosagepresent, &sample_uidx);
              const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
              const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
              if (halfdist < dosage_erase_halfdist) {
                write_dosagevals[dosage_write_idx++] = dosage_int;
              } else {
                ClearBit(sample_uidx, write_dosagepresent);
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
        if (set_hh_missing && is_haploid_nonmt) {
          if (is_x_or_y) {
            if (!set_hh_missing_keep_dosage) {
              // need to erase dosages associated with the hardcalls we're
              // about to clear

              // male hets to missing
              SetMaleHetMissingCleardosage(sex_male_collapsed, sex_male_collapsed_interleaved, sample_ctv2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            } else {
              // need to generate a new unphased dosage for each cleared
              // hardcall lacking a dosage entry
              SetMaleHetMissingKeepdosage(sex_male_collapsed, sex_male_collapsed_interleaved, sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            }
            if (is_y) {
              // all female calls to missing; unknown-sex calls now left
              // alone
              InterleavedSetMissingCleardosage(sex_female_collapsed, sex_female_collapsed_interleaved, sample_ctv2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
              is_hphase = 0;
            } else if (is_hphase && cur_write_phasepresent) {
              MaskGenovecHetsUnsafe(write_genovec, sample_ctl2, cur_write_phasepresent);
            }
          } else {
            // all hets to missing
            // may want to move is_hphase zeroing in front
            if (!set_hh_missing_keep_dosage) {
              SetHetMissingCleardosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            } else {
              SetHetMissingKeepdosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            }
            is_hphase = 0;
          }
        } else if (set_mixed_mt_missing && is_mt) {
          if (!set_mixed_mt_missing_keep_dosage) {
            // all hets to missing
            SetHetMissingCleardosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          } else {
            // todo: keep dosage-phase information
            SetHetMissingKeepdosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          }
          is_hphase = 0;
        }
        ZeroTrailingQuaters(sample_ct, write_genovec);
        // todo: --set-me-missing, --zero-cluster, --fill-missing-with-ref
        if (spgwp) {
          if (pwcp->fwrite_bufp >= &(pwcp->fwrite_buf[kPglFwriteBlockSize])) {
            const uintptr_t cur_byte_ct = pwcp->fwrite_bufp - pwcp->fwrite_buf;
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
          PwcAppendBiallelicGenovecDosage16(write_genovec, write_dosagepresent, write_dosagevals, write_dosage_ct, pwcp);
        } else {
          PwcAppendBiallelicGenovecHphaseDosage16(write_genovec, cur_write_phasepresent, write_phaseinfo, write_dosagepresent, write_dosagevals, write_dosage_ct, pwcp);
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

PgenGlobalFlags GflagsVfilter(const uintptr_t* variant_include, const unsigned char* vrtypes, uint32_t raw_variant_ct, PgenGlobalFlags input_gflags) {
  PgenGlobalFlags read_phase_dosage_gflags = kfPgenGlobal0;
  const uintptr_t* vrtypes_alias_iter = R_CAST(const uintptr_t*, vrtypes);
  const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
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
// (Note that MakePlink2NoVsort() requires enough memory for 64k * 2
// variants per output thread, due to LD compression.  This is faster in the
// common case, but once you have 150k+ samples with dosage data...)
PglErr MakePgenRobust(const uintptr_t* sample_include, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const AltAlleleCt* refalt1_select, const uint32_t* new_variant_idx_to_old, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  // variant_uidx_new_to_old[] can be nullptr

  // caller responsible for initializing g_cip (may need to be different from
  // initial cip struct)
  unsigned char* bigstack_mark = g_bigstack_base;
  ThreadsState ts;
  InitThreads3z(&ts);
  STPgenWriter spgw;
  PglErr reterr = kPglRetSuccess;
  PreinitSpgw(&spgw);
  {
    // g_plink2_write_flags assumed to include --set-hh-missing and
    //   --set-mixed-mt-missing
    // g_sex_{fe}male_collapsed_interleaved assumed to be initialized if
    //   necessary

    if (bigstack_alloc_thread(1, &ts.threads)) {
      goto MakePgenRobust_ret_NOMEM;
    }
    ts.calc_thread_ct = 1;
    g_spgwp = &spgw;
    const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    g_sample_include = subsetting_required? sample_include : nullptr;
    g_new_sample_idx_to_old = new_sample_idx_to_old;
    g_raw_sample_ct = raw_sample_ct;
    g_sample_ct = sample_ct;
    g_error_ret = kPglRetSuccess;
    const uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
    if ((make_plink2_flags & kfMakeBed) || ((make_plink2_flags & (kfMakePgen | (kfMakePgenFormatBase * 3))) == (kfMakePgen | kfMakePgenFormatBase))) {
      // g_calc_thread_ct = 1;
      logerrputs("Error: Fixed-width .bed/.pgen output doesn't support sorting yet.  Generate a\nregular sorted .pgen first, and then reformat it.\n");
      reterr = kPglRetNotYetSupported;
      goto MakePgenRobust_ret_1;
    } else {
      const uint32_t input_biallelic = (!variant_allele_idxs);
      const uintptr_t* write_allele_idx_offsets = nullptr;
      if (!input_biallelic) {
        if ((variant_ct < raw_variant_ct) || new_variant_idx_to_old_iter) {
          uintptr_t* new_allele_idx_offsets;
          if (bigstack_alloc_ul(variant_ct + 1, &new_allele_idx_offsets)) {
            goto MakePgenRobust_ret_NOMEM;
          }
          uintptr_t cur_offset = 0;
          // todo: separate trim-alts case
          if (!new_variant_idx_to_old_iter) {
            uint32_t variant_uidx = 0;
            for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
              FindFirst1BitFromU32(variant_include, &variant_uidx);
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
            logputs("Error: Multiallelic .pgen write is not yet supported.\n");
            reterr = kPglRetNotYetSupported;
            goto MakePgenRobust_ret_1;
          } else {
            BigstackReset(new_allele_idx_offsets);
          }
        } else {
          write_allele_idx_offsets = variant_allele_idxs;
        }
      }
      if ((variant_ct == raw_variant_ct) || new_variant_idx_to_old_iter) {
        g_write_chr_fo_vidx_start = g_cip->chr_fo_vidx_start;
      } else {
        if (AllocAndFillSubsetChrFoVidxStart(variant_include, g_cip, &g_write_chr_fo_vidx_start)) {
          goto MakePgenRobust_ret_NOMEM;
        }
      }
      PgenGlobalFlags read_phase_dosage_gflags = simple_pgrp->fi.gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      if (make_plink2_flags & kfMakePgenErasePhase) {
        read_phase_dosage_gflags &= ~(kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
      }
      if (make_plink2_flags & kfMakePgenEraseDosage) {
        if (hard_call_thresh == UINT32_MAX) {
          read_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
          g_plink2_write_flags |= kfPlink2WriteLateDosageErase;
        }
      }
      if (read_phase_dosage_gflags && (variant_ct < raw_variant_ct)) {
        read_phase_dosage_gflags = GflagsVfilter(variant_include, simple_pgrp->fi.vrtypes, raw_variant_ct, simple_pgrp->fi.gflags);
      }
      g_read_phase_dosage_gflags = read_phase_dosage_gflags;
      g_hard_call_halfdist = (hard_call_thresh == UINT32_MAX)? 0 : (kDosage4th - hard_call_thresh);
      g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_hphase_present = (read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
      const uint32_t read_dosage_present = (read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
      PgenGlobalFlags write_phase_dosage_gflags = read_phase_dosage_gflags;
      uint32_t read_or_write_dosage_present = read_dosage_present;
      if (g_plink2_write_flags & kfPlink2WriteLateDosageErase) {
        write_phase_dosage_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      } else if (make_plink2_flags & (kfPlink2WriteSetHhMissingKeepDosage | kfPlink2WriteSetMixedMtMissingKeepDosage)) {
        // command-line parser guarantees erase-dosage and
        // --set-hh-missing/--set-mixed-mt-missing keep-dosage aren't used
        // together
        read_or_write_dosage_present = 1;

        // could verify at least one het haploid is present before setting this
        // flag...
        write_phase_dosage_gflags |= kfPgenGlobalDosagePresent;
        if (make_plink2_flags & kfPlink2WriteSetMixedMtMissingKeepDosage) {
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
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, write_allele_idx_offsets, simple_pgrp->fi.nonref_flags, variant_ct, sample_ct, write_phase_dosage_gflags, nonref_flags_storage, g_spgwp, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (reterr) {
        goto MakePgenRobust_ret_1;
      }
      unsigned char* spgw_alloc;
      if (bigstack_alloc_ulp(1, &(g_loadbuf_thread_starts[0])) ||
          bigstack_alloc_ulp(1, &(g_loadbuf_thread_starts[1])) ||
          bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
        goto MakePgenRobust_ret_NOMEM;
      }
      SpgwInitPhase2(max_vrec_len, g_spgwp, spgw_alloc);

      const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
      if (new_sample_idx_to_old || subsetting_required) {
        if (bigstack_alloc_ulp(1, &g_thread_write_genovecs)) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if (read_hphase_present && new_sample_idx_to_old) {
          if (bigstack_alloc_ui(raw_sample_ct, &g_old_sample_idx_to_new)) {
            goto MakePgenRobust_ret_NOMEM;
          }
          for (uint32_t new_sample_idx = 0; new_sample_idx < sample_ct; ++new_sample_idx) {
            g_old_sample_idx_to_new[new_sample_idx_to_old[new_sample_idx]] = new_sample_idx;
          }
        }
        if (input_biallelic) {
          if (bigstack_alloc_ul(sample_ctl2, &(g_thread_write_genovecs[0]))) {
            goto MakePgenRobust_ret_NOMEM;
          }
        } else {
          if (bigstack_alloc_ul(DivUp(2 * sample_ct * sizeof(AltAlleleCt), kBytesPerWord), &(g_thread_write_genovecs[0]))) {
            goto MakePgenRobust_ret_NOMEM;
          }
        }
      } else {
        g_thread_write_genovecs = nullptr;
      }
      const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
      if (read_hphase_present) {
        if (bigstack_alloc_ulp(1, &g_thread_write_phasepresents) ||
            bigstack_alloc_ulp(1, &g_thread_write_phaseinfos) ||
            bigstack_alloc_ulp(1, &g_thread_all_hets) ||
            bigstack_alloc_ul(sample_ctl, &(g_thread_write_phasepresents[0])) ||
            bigstack_alloc_ul(sample_ctl, &(g_thread_write_phaseinfos[0])) ||
            bigstack_alloc_ul(raw_sample_ctl, &(g_thread_all_hets[0]))) {
          goto MakePgenRobust_ret_NOMEM;
        }
      }
      if (read_or_write_dosage_present) {
        if (bigstack_alloc_dosagep(1, &g_thread_write_dosagevals) ||
            bigstack_alloc_ulp(1, &g_thread_write_dosagepresents) ||
            bigstack_alloc_ul(sample_ctl, &(g_thread_write_dosagepresents[0])) ||
            bigstack_alloc_dosage(sample_ct, &(g_thread_write_dosagevals[0]))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if (read_dosage_present && new_sample_idx_to_old) {
          if (bigstack_alloc_uip(1, &g_thread_cumulative_popcount_bufs) ||
              bigstack_alloc_ui(raw_sample_ctl, &(g_thread_cumulative_popcount_bufs[0]))) {
            goto MakePgenRobust_ret_NOMEM;
          }
        }
      }
      g_refalt1_select = refalt1_select;
      if (refalt1_select) {
        if (variant_ct < raw_variant_ct) {
          // might want inner loop to map variant uidx -> idx instead
          g_refalt1_select = S_CAST(AltAlleleCt*, bigstack_alloc(variant_ct * 2 * sizeof(AltAlleleCt)));
          if (!g_refalt1_select) {
            goto MakePgenRobust_ret_NOMEM;
          }
          uint32_t variant_uidx = 0;
          for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
            FindFirst1BitFromU32(variant_include, &variant_uidx);
            memcpy(K_CAST(AltAlleleCt*, &(g_refalt1_select[2 * variant_idx])), &(refalt1_select[2 * variant_uidx]), 2 * sizeof(AltAlleleCt));
          }
        } else {
          assert(!new_variant_idx_to_old_iter);
        }
      }

      const uint32_t raw_sample_ctv2 = QuaterCtToVecCt(raw_sample_ct);
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
          const uintptr_t phaseraw_word_ct = kWordsPerVec + RoundDownPow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
          load_variant_vec_ct += WordCtToVecCt(phaseraw_word_ct);
        }
        if (read_dosage_present) {
          // todo: phased dosage

          // (unphased, biallelic) dosageraw has two parts:
          // 1. vec-aligned bitarray of up to raw_sample_ct bits, storing which
          //    samples have dosages.
          // 2. word-aligned array of uint16s with 0..32768 fixed-point dosages.
          assert(!(read_phase_dosage_gflags & kfPgenGlobalDosagePhasePresent));
          const uintptr_t dosageraw_word_ct = kWordsPerVec * (BitCtToVecCt(raw_sample_ct) + DivUp(raw_sample_ct, (kBytesPerVec / sizeof(Dosage))));
          load_variant_vec_ct += WordCtToVecCt(dosageraw_word_ct);
        }
      }
      // todo: multiallelic variants

      uintptr_t bytes_left = bigstack_left();
      if (bytes_left < 7 * kCacheline) {
        goto MakePgenRobust_ret_NOMEM;
      }
      bytes_left -= 7 * kCacheline;  // defend against adverse rounding
      uintptr_t ulii = bytes_left / (2 * (kBytesPerVec * load_variant_vec_ct + loaded_vrtypes_needed));
      if (!ulii) {
        goto MakePgenRobust_ret_NOMEM;
      }
      if (ulii > MINV(kPglVblockSize, variant_ct)) {
        ulii = MINV(kPglVblockSize, variant_ct);
      }
      const uint32_t write_block_size = ulii;
      uintptr_t* main_loadbufs[2];
      main_loadbufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(load_variant_vec_ct * kBytesPerVec * write_block_size));
      main_loadbufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(load_variant_vec_ct * kBytesPerVec * write_block_size));
      if (loaded_vrtypes_needed) {
        g_loaded_vrtypes[0] = S_CAST(unsigned char*, bigstack_alloc_raw_rd(write_block_size));
        g_loaded_vrtypes[1] = S_CAST(unsigned char*, bigstack_alloc_raw_rd(write_block_size));
      } else {
        g_loaded_vrtypes[0] = nullptr;
        g_loaded_vrtypes[1] = nullptr;
      }

      logprintfww5("Writing %s ... ", outname);
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
      uint32_t read_variant_uidx = UINT32_MAX;  // deliberate overflow
      PgrClearLdCache(simple_pgrp);
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
              FindFirst1BitFromU32(variant_include, &read_variant_uidx);
            } else {
              read_variant_uidx = *new_variant_idx_to_old_iter++;
            }
            reterr = PgrReadRaw(read_variant_uidx, read_phase_dosage_gflags, simple_pgrp, &loadbuf_iter, cur_loaded_vrtypes? (&(cur_loaded_vrtypes[uii])) : nullptr);
            if (reterr) {
              if (reterr == kPglRetMalformedInput) {
                logputs("\n");
                logerrputs("Error: Malformed .pgen file.\n");
              }
              goto MakePgenRobust_ret_1;
            }
          }
        }
        if (read_batch_idx) {
          JoinThreads3z(&ts);
          reterr = g_error_ret;
          if (reterr) {
            goto MakePgenRobust_ret_WRITE_FAIL;
          }
        }
        if (!ts.is_last_block) {
          g_cur_block_write_ct = cur_batch_size;
          ts.is_last_block = (read_batch_idx == batch_ct_m1);
          ts.thread_func_ptr = MakePgenThread;
          if (SpawnThreads3z(read_batch_idx, &ts)) {
            goto MakePgenRobust_ret_THREAD_CREATE_FAIL;
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
            next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
          }
        }
        ++read_batch_idx;
      }
      SpgwFinish(g_spgwp);
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logprintf("done.\n");
    }
  }
  while (0) {
  MakePgenRobust_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakePgenRobust_ret_WRITE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  MakePgenRobust_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 MakePgenRobust_ret_1:
  CleanupThreads3z(&ts, &g_cur_block_write_ct);
  if (SpgwCleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr MakePlink2NoVsort(const char* xheader, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  ThreadsState ts;
  InitThreads3z(&ts);
  MTPgenWriter* mpgwp = nullptr;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (make_plink2_flags & kfMakePlink2MMask) {
      logerrputs("Error: --make-bed/--make-{b}pgen multiallelics= is currently under development.\n");
      reterr = kPglRetNotYetSupported;
      goto MakePlink2NoVsort_ret_1;
    }
    g_plink2_write_flags = kfPlink2Write0;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (make_plink2_flags & kfMakePlink2SetHhMissing) {
      const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
      uintptr_t* new_sex_male;
      uintptr_t* sex_female;
      if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &new_sex_male) ||
          bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
          bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed) ||
          bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed_interleaved) ||
          bigstack_alloc_ul(raw_sample_ctl, &sex_female)) {
        goto MakePlink2NoVsort_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, sample_include, sample_ct, new_sex_male);
      ZeroTrailingWords(BitCtToWordCt(sample_ct), new_sex_male);
      g_sex_male = new_sex_male;
      FillInterleavedMaskVec(g_sex_male, sample_ctv, g_sex_male_collapsed_interleaved);

      BitvecAndNotCopy(sex_nm, sex_male, raw_sample_ctl, sex_female);
      CopyBitarrSubset(sex_female, sample_include, sample_ct, g_sex_female_collapsed);
      FillInterleavedMaskVec(g_sex_female_collapsed, sample_ctv, g_sex_female_collapsed_interleaved);

      BigstackReset(sex_female);
      g_plink2_write_flags |= kfPlink2WriteSetHhMissing;
      if (make_plink2_flags & kfMakePlink2SetHhMissingKeepDosage) {
        g_plink2_write_flags |= kfPlink2WriteSetHhMissingKeepDosage;
      }
    }
    if (make_plink2_flags & kfMakePlink2SetMixedMtMissing) {
      g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
      if (make_plink2_flags & kfMakePlink2SetMixedMtMissingKeepDosage) {
        g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissingKeepDosage;
      }
    }
    g_cip = cip;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    const uint32_t make_pgen = make_plink2_flags & kfMakePgen;
    // todo: prohibit .pgen + .bim write when data is multiallelic without
    //   either multiallelic split or erase-alt2+ specified
    //   (--make-bed = automatic erase-alt2+)
    if ((make_plink2_flags & kfMakeBed) || ((make_plink2_flags & (kfMakePgen | (kfMakePgenFormatBase * 3))) == (kfMakePgen | kfMakePgenFormatBase))) {
      // fixed-width
      if (make_pgen) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      } else {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".bed");
      }
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
        goto MakePlink2NoVsort_ret_OPEN_FAIL;
      }
      if (make_pgen) {
        fwrite_unlocked("l\x1b\x02", 3, 1, outfile);
        fwrite_unlocked(&variant_ct, 4, 1, outfile);
        fwrite_unlocked(&sample_ct, 4, 1, outfile);
        if (!pgfip->nonref_flags) {
          const PgenGlobalFlags gflags = pgfip->gflags;
          uint32_t uii = 64;
          if (gflags & kfPgenGlobalAllNonref) {
            uii = 128;
          }
          putc_unlocked(uii, outfile);
        } else {
          putc_unlocked(192, outfile);
          fwrite_unlocked(pgfip->nonref_flags, DivUp(variant_ct, CHAR_BIT), 1, outfile);
        }
        if (ferror_unlocked(outfile)) {
          goto MakePlink2NoVsort_ret_WRITE_FAIL;
        }
      } else {
        if (fwrite_checked("l\x1b\x01", 3, outfile)) {
          goto MakePlink2NoVsort_ret_WRITE_FAIL;
        }
      }
      logprintfww5("Writing %s ... ", outname);
      fputs("0%", stdout);
      fflush(stdout);
      uint32_t pct = 0;
      if (variant_ct && sample_ct) {
        const uintptr_t sample_ct4 = QuaterCtToByteCt(sample_ct);
        if (bigstack_alloc_ui(raw_sample_ctl, &g_sample_include_cumulative_popcounts) ||
            bigstack_alloc_uc(sample_ct4 * kPglVblockSize, &(g_writebufs[0])) ||
            bigstack_alloc_uc(sample_ct4 * kPglVblockSize, &(g_writebufs[1]))) {
          // todo: low-memory single-threaded fallback mode
          goto MakePlink2NoVsort_ret_NOMEM;
        }
        FillCumulativePopcounts(sample_include, raw_sample_ctl, g_sample_include_cumulative_popcounts);
        // tried more threads, pointless since this is too I/O-bound
        // (exception: reordering samples)
        uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
        g_collapsed_sort_map = new_sample_idx_to_old;
        if (!new_sample_idx_to_old) {
          // Without BMI2 instructions, subsetting is most expensive with
          // sample_ct near 2/3 of raw_sample_ct; up to ~7 compute threads are
          // useful in that case.  (See CopyQuaterarrNonemptySubset().)
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
          uint32_t* new_collapsed_sort_map;
          if (bigstack_alloc_ui(sample_ct, &new_collapsed_sort_map)) {
            goto MakePlink2NoVsort_ret_NOMEM;
          }
          UidxsToIdxs(sample_include, g_sample_include_cumulative_popcounts, sample_ct, new_collapsed_sort_map);
          g_collapsed_sort_map = new_collapsed_sort_map;
        }

        if (make_plink2_flags & kfMakeBed) {
          g_plink2_write_flags |= kfPlink2WritePlink1;
        }

        g_hard_call_halfdist = 0;
        if ((hard_call_thresh != UINT32_MAX) && (pgfip->gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))) {
          g_hard_call_halfdist = kDosage4th - hard_call_thresh;
        }
        unsigned char* main_loadbufs[2];
        uint32_t read_block_size;
        if (PgenMtLoadInit(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, g_hard_call_halfdist? (&g_dosage_presents) : nullptr, g_hard_call_halfdist? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &ts.threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
          goto MakePlink2NoVsort_ret_NOMEM;
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

        const uint32_t read_block_sizel = BitCtToWordCt(read_block_size);
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
              cur_block_write_ct = PopcountWords(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
              if (cur_block_write_ct) {
                break;
              }
              ++read_block_idx;
            }
            if (read_block_idx == read_block_ct_m1) {
              cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
              cur_block_write_ct = PopcountWords(&(variant_include[read_block_idx * read_block_sizel]), BitCtToWordCt(cur_read_block_size));
            }
            if (PgfiMultiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
              goto MakePlink2NoVsort_ret_READ_FAIL;
            }
          }
          if (variant_idx) {
            JoinThreads3z(&ts);
            if (reterr) {
              if (reterr == kPglRetMalformedInput) {
                logputs("\n");
                logerrputs("Error: Malformed .pgen file.\n");
              }
              goto MakePlink2NoVsort_ret_1;
            }
          }
          if (!ts.is_last_block) {
            g_cur_block_write_ct = cur_block_write_ct;
            ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
            for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
              g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
              g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
            }
            ts.is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
            ts.thread_func_ptr = MakeBedlikeThread;
            if (SpawnThreads3z(variant_idx, &ts)) {
              goto MakePlink2NoVsort_ret_THREAD_CREATE_FAIL;
            }
          }
          parity = 1 - parity;
          if (variant_idx) {
            // write *previous* block results
            if (fwrite_checked(g_writebufs[parity], (variant_idx - prev_variant_idx) * sample_ct4, outfile)) {
              goto MakePlink2NoVsort_ret_WRITE_FAIL;
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
              next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
            }
            prev_variant_idx = variant_idx;
          }
          ++read_block_idx;
          variant_idx += cur_block_write_ct;
          // crucially, this is independent of the PgenReader block_base
          // pointers
          pgfip->block_base = main_loadbufs[parity];
        }
      }
      if (fclose_null(&outfile)) {
        goto MakePlink2NoVsort_ret_WRITE_FAIL;
      }
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logprintf("done.\n");
      BigstackReset(bigstack_mark);
    } else if (make_pgen) {
      if ((!variant_ct) || (!sample_ct)) {
        logerrputs("Error: Zero-variant/zero-sample .pgen writing is not currently supported.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePlink2NoVsort_ret_1;
      }
      const uint32_t input_biallelic = (!variant_allele_idxs);
      const uintptr_t* write_allele_idx_offsets = nullptr;
      if (!input_biallelic) {
        if (variant_ct < raw_variant_ct) {
          uintptr_t* new_allele_idx_offsets;
          if (bigstack_alloc_ul(variant_ct + 1, &new_allele_idx_offsets)) {
            goto MakePlink2NoVsort_fallback;
          }
          uintptr_t cur_offset = 0;
          uint32_t variant_uidx = 0;
          // todo: separate trim-alts case
          for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
            FindFirst1BitFromU32(variant_include, &variant_uidx);
            new_allele_idx_offsets[variant_idx] = cur_offset;
            cur_offset += variant_allele_idxs[variant_uidx + 1] - variant_allele_idxs[variant_uidx];
          }
          if (cur_offset != 2 * variant_ct) {
            new_allele_idx_offsets[variant_ct] = cur_offset;
            write_allele_idx_offsets = new_allele_idx_offsets;
            logputs("Error: Multiallelic .pgen write is not yet supported.\n");
            reterr = kPglRetNotYetSupported;
            goto MakePlink2NoVsort_ret_1;
          } else {
            BigstackReset(new_allele_idx_offsets);
          }
        } else {
          write_allele_idx_offsets = variant_allele_idxs;
        }
      }
      if (variant_ct == raw_variant_ct) {
        g_write_chr_fo_vidx_start = cip->chr_fo_vidx_start;
      } else {
        if (AllocAndFillSubsetChrFoVidxStart(variant_include, cip, &g_write_chr_fo_vidx_start)) {
          goto MakePlink2NoVsort_fallback;
        }
      }
      PgenGlobalFlags read_phase_dosage_gflags = pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      if (make_plink2_flags & kfMakePgenErasePhase) {
        read_phase_dosage_gflags &= ~(kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
      }
      if (make_plink2_flags & kfMakePgenEraseDosage) {
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
        read_phase_dosage_gflags = GflagsVfilter(variant_include, pgfip->vrtypes, raw_variant_ct, pgfip->gflags);
      }
      // could check if all the phased samples were also filtered out, but
      // that's already caught by running --make-pgen twice, so not a big deal
      g_read_phase_dosage_gflags = read_phase_dosage_gflags;
      g_hard_call_halfdist = (hard_call_thresh == UINT32_MAX)? 0 : (kDosage4th - hard_call_thresh);
      g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_hphase_present = (read_phase_dosage_gflags / kfPgenGlobalHardcallPhasePresent) & 1;
      const uint32_t read_dosage_present = (read_phase_dosage_gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))? 1 : 0;
      PgenGlobalFlags write_phase_dosage_gflags = read_phase_dosage_gflags;
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
      MpgwInitPhase1(write_allele_idx_offsets, variant_ct, sample_ct, write_phase_dosage_gflags, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);

      // bugfix: each variant currently needs to be vector-aligned
      // bugfix?: need to use raw_sample_ct here, not sample_ct
      const uint32_t raw_sample_ctv2 = QuaterCtToVecCt(raw_sample_ct);
      const uint32_t max_vblock_size = MINV(kPglVblockSize, variant_ct);
      uint64_t load_vblock_cacheline_ct = VecCtToCachelineCtU64(S_CAST(uint64_t, raw_sample_ctv2) * max_vblock_size);

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
        const uintptr_t phaseraw_word_ct = kWordsPerVec + RoundDownPow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
        load_vblock_cacheline_ct += WordCtToCachelineCtU64(S_CAST(uint64_t, phaseraw_word_ct) * max_vblock_size);
      }
      if (read_dosage_present) {
        // todo: phased dosage

        // (unphased, biallelic) dosageraw has two parts:
        // 1. vec-aligned bitarray of up to raw_sample_ct bits, storing which
        //    samples have dosages.
        // 2. word-aligned array of uint16s with 0..32768 fixed-point dosages.
        assert(!(read_phase_dosage_gflags & kfPgenGlobalDosagePhasePresent));
        const uintptr_t dosageraw_word_ct = kWordsPerVec * (BitCtToVecCt(raw_sample_ct) + DivUp(raw_sample_ct, kBytesPerVec / sizeof(Dosage)));
        load_vblock_cacheline_ct += WordCtToCachelineCtU64(dosageraw_word_ct * S_CAST(uint64_t, max_vblock_size));
      }
      // todo: multiallelic variants

#ifndef __LP64__
      if ((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (load_vblock_cacheline_ct > (0x7fffffff / kCacheline))) {
        goto MakePlink2NoVsort_fallback;
      }
#endif
      uint32_t calc_thread_ct = DivUp(variant_ct, kPglVblockSize);
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
        g_refalt1_select = S_CAST(AltAlleleCt*, bigstack_alloc(variant_ct * 2 * sizeof(AltAlleleCt)));
        if (!g_refalt1_select) {
          goto MakePlink2NoVsort_fallback;
        }
        uint32_t variant_uidx = 0;
        for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
          FindFirst1BitFromU32(variant_include, &variant_uidx);
          memcpy(K_CAST(AltAlleleCt*, &(g_refalt1_select[2 * variant_idx])), &(refalt1_select[2 * variant_uidx]), 2 * sizeof(AltAlleleCt));
        }
      }
      mpgwp = S_CAST(MTPgenWriter*, bigstack_alloc((calc_thread_ct + DivUp(sizeof(MTPgenWriter), kBytesPerWord)) * sizeof(intptr_t)));
      if (!mpgwp) {
        goto MakePlink2NoVsort_fallback;
      }
      mpgwp->pgen_outfile = nullptr;
      if (bigstack_alloc_thread(calc_thread_ct, &ts.threads) ||
          bigstack_alloc_ulp(calc_thread_ct, &(g_loadbuf_thread_starts[0])) ||
          bigstack_alloc_ulp(calc_thread_ct, &(g_loadbuf_thread_starts[1]))) {
        goto MakePlink2NoVsort_fallback;
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
          goto MakePlink2NoVsort_fallback;
        }
        if (read_hphase_present && new_sample_idx_to_old) {
          if (bigstack_alloc_ui(raw_sample_ct, &g_old_sample_idx_to_new)) {
            goto MakePlink2NoVsort_fallback;
          }
          for (uint32_t new_sample_idx = 0; new_sample_idx < sample_ct; ++new_sample_idx) {
            g_old_sample_idx_to_new[new_sample_idx_to_old[new_sample_idx]] = new_sample_idx;
          }
        }
        if (read_dosage_present && new_sample_idx_to_old) {
          // g_thread_cumulative_popcount_bufs
          other_per_thread_cacheline_ct += Int32CtToCachelineCt(raw_sample_ctl);
        }
        // per-thread output buffers required
        if (input_biallelic) {
          other_per_thread_cacheline_ct += QuaterCtToCachelineCt(sample_ct);
        } else {
          other_per_thread_cacheline_ct += DivUp(2 * sample_ct * sizeof(AltAlleleCt), kCacheline);
        }
      }
      if (read_hphase_present || read_or_write_dosage_present) {
        if (read_hphase_present) {
          if (bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_phasepresents) ||
              bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_phaseinfos) ||
              bigstack_alloc_ulp(calc_thread_ct, &g_thread_all_hets)) {
            goto MakePlink2NoVsort_fallback;
          }
          // phasepresent, phaseinfo
          other_per_thread_cacheline_ct += 2 * BitCtToCachelineCt(sample_ct);

          // all_hets
          other_per_thread_cacheline_ct += BitCtToCachelineCt(raw_sample_ct);
        }
        if (read_or_write_dosage_present) {
          if (bigstack_alloc_dosagep(calc_thread_ct, &g_thread_write_dosagevals) ||
              bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_dosagepresents)) {
            goto MakePlink2NoVsort_fallback;
          }
          if (read_dosage_present && new_sample_idx_to_old) {
            if (bigstack_alloc_uip(calc_thread_ct, &g_thread_cumulative_popcount_bufs)) {
              goto MakePlink2NoVsort_fallback;
            }
          }
          // dosage_present
          other_per_thread_cacheline_ct += BitCtToCachelineCt(sample_ct);

          // dosage_vals
          other_per_thread_cacheline_ct += DivUp(sample_ct, (kCacheline / sizeof(Dosage)));
        }
        // g_loaded_vrtypes
        other_per_thread_cacheline_ct += 2 * (kPglVblockSize / kCacheline);
      }
      const uintptr_t cachelines_avail = bigstack_left() / kCacheline;
      if (cachelines_avail < alloc_base_cacheline_ct + (mpgw_per_thread_cacheline_ct + other_per_thread_cacheline_ct) * calc_thread_ct) {
        if (cachelines_avail < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct + other_per_thread_cacheline_ct) {
          goto MakePlink2NoVsort_fallback;
        }
        calc_thread_ct = (cachelines_avail - alloc_base_cacheline_ct) / (mpgw_per_thread_cacheline_ct + other_per_thread_cacheline_ct);
      }
      uintptr_t* main_loadbufs[2];
      main_loadbufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw(load_vblock_cacheline_ct * calc_thread_ct * kCacheline));
      main_loadbufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw(load_vblock_cacheline_ct * calc_thread_ct * kCacheline));
      if (read_hphase_present || read_or_write_dosage_present) {
        if (read_hphase_present || read_dosage_present) {
          g_loaded_vrtypes[0] = S_CAST(unsigned char*, bigstack_alloc_raw(kPglVblockSize * calc_thread_ct));
          g_loaded_vrtypes[1] = S_CAST(unsigned char*, bigstack_alloc_raw(kPglVblockSize * calc_thread_ct));
        }
        const uint32_t bitvec_writebuf_byte_ct = BitCtToCachelineCt(sample_ct) * kCacheline;
        const uintptr_t dosagevals_writebuf_byte_ct = DivUp(sample_ct, (kCacheline / 2)) * kCacheline;
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          if (read_hphase_present) {
            g_thread_write_phasepresents[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));
            g_thread_write_phaseinfos[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));

            g_thread_all_hets[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(BitCtToCachelineCt(raw_sample_ct) * kCacheline));
          }
          if (read_or_write_dosage_present) {
            g_thread_write_dosagepresents[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));
            g_thread_write_dosagevals[tidx] = S_CAST(Dosage*, bigstack_alloc_raw(dosagevals_writebuf_byte_ct));
            if (read_dosage_present && new_sample_idx_to_old) {
              g_thread_cumulative_popcount_bufs[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(Int32CtToCachelineCt(raw_sample_ctl) * kCacheline));
            }
          }
        }
      } else {
        g_loaded_vrtypes[0] = nullptr;
        g_loaded_vrtypes[1] = nullptr;
      }
      if (new_sample_idx_to_old || subsetting_required) {
        uintptr_t writebuf_byte_ct = input_biallelic? QuaterCtToByteCt(sample_ct) : (2 * sample_ct * sizeof(AltAlleleCt));
        writebuf_byte_ct = RoundUpPow2(writebuf_byte_ct, kCacheline);
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          g_thread_write_genovecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(writebuf_byte_ct));
        }
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      logprintfww5("Writing %s ... ", outname);
      fputs("0%", stdout);
      fflush(stdout);
      unsigned char* mpgw_alloc = S_CAST(unsigned char*, bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline));
      assert(g_bigstack_base <= g_bigstack_end);
      reterr = MpgwInitPhase2(outname, write_allele_idx_offsets, pgfip->nonref_flags, variant_ct, sample_ct, write_phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);
      if (reterr) {
        goto MakePlink2NoVsort_ret_1;
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
      uint32_t read_variant_uidx = FindFirst1BitFrom(variant_include, 0);
      PgrClearLdCache(simple_pgrp);
      while (1) {
        if (read_batch_idx) {
          g_cur_block_write_ct = cur_batch_size;
          ts.is_last_block = (write_idx_end == variant_ct);
          ts.thread_func_ptr = MakePgenThread;
          if (SpawnThreads3z(read_batch_idx - 1, &ts)) {
            goto MakePlink2NoVsort_ret_THREAD_CREATE_FAIL;
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
            FindFirst1BitFromU32(variant_include, &read_variant_uidx);
            reterr = PgrReadRaw(read_variant_uidx, read_phase_dosage_gflags, simple_pgrp, &loadbuf_iter, cur_loaded_vrtypes? (&(cur_loaded_vrtypes[uii])) : nullptr);
            if (reterr) {
              if (reterr == kPglRetMalformedInput) {
                logputs("\n");
                logerrputs("Error: Malformed .pgen file.\n");
              }
              goto MakePlink2NoVsort_ret_1;
            }
          }
        }
        if (read_batch_idx) {
          JoinThreads3z(&ts);
        }
        parity = 1 - parity;
        if (write_idx_end) {
          reterr = MpgwFlush(mpgwp);
          if (reterr) {
            goto MakePlink2NoVsort_ret_WRITE_FAIL;
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
            next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
          }
        }
        ++read_batch_idx;
        write_idx_end += cur_batch_size;
      }
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logprintf("done.\n");
      BigstackReset(bigstack_mark);
    } else if (0) {
    MakePlink2NoVsort_fallback:
      mpgwp = nullptr;
      BigstackReset(bigstack_mark2);
      reterr = MakePgenRobust(sample_include, new_sample_idx_to_old, variant_include, variant_allele_idxs, refalt1_select, nullptr, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, hard_call_thresh, dosage_erase_thresh, make_plink2_flags, simple_pgrp, outname, outname_end);
      if (reterr) {
        goto MakePlink2NoVsort_ret_1;
      }
    }
    const uint32_t trim_alts = S_CAST(uint32_t, make_plink2_flags & kfMakePlink2TrimAlts);
    // don't bother with trim-alts/set-hh-missing interaction for now
    if (make_plink2_flags & kfMakeBim) {
      const uint32_t bim_zst = (make_plink2_flags / kfMakeBimZs) & 1;
      OutnameZstSet(".bim", bim_zst, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteMapOrBim(outname, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, '\t', bim_zst, max_thread_ct);
      if (reterr) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePvar) {
      OutnameZstSet(".pvar", pvar_psam_flags & kfPvarZs, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t nonref_flags_storage = 3;
      if (!pgfip->nonref_flags) {
        nonref_flags_storage = (pgfip->gflags & kfPgenGlobalAllNonref)? 2 : 1;
      }
      reterr = WritePvar(outname, xheader, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, xheader_info_pr, xheader_info_pr_nonflag, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_flags, max_thread_ct);
      if (reterr) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakeFam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, new_sample_idx_to_old, sample_ct, pheno_ct, '\t');
      if (reterr) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePsam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WritePsam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, sample_ct, pheno_ct, max_pheno_name_blen, pvar_psam_flags);
      if (reterr) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
  }
  while (0) {
  MakePlink2NoVsort_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakePlink2NoVsort_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  MakePlink2NoVsort_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  MakePlink2NoVsort_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  MakePlink2NoVsort_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 MakePlink2NoVsort_ret_1:
  if (MpgwCleanup(mpgwp) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  CleanupThreads3z(&ts, &g_cur_block_write_ct);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackReset(bigstack_mark);
  return reterr;
}


BoolErr SortChr(const ChrInfo* cip, const uint32_t* chr_idx_to_size, uint32_t use_nsort, ChrInfo* write_cip) {
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
        *std_sortbuf_iter++ = (S_CAST(uint64_t, xymt_idx_to_chr_sort_offset[xymt_idx] + autosome_ct_p1) << 32) | (xymt_idx + autosome_ct_p1);
      }
    }
  }
  const uint32_t std_sortbuf_len = std_sortbuf_iter - std_sortbuf;
#ifdef __cplusplus
  std::sort(std_sortbuf, &(std_sortbuf[std_sortbuf_len]));
#else
  qsort(std_sortbuf, std_sortbuf_len, sizeof(int64_t), uint64cmp);
#endif
  uint32_t write_vidx = 0;
  write_cip->chr_fo_vidx_start[0] = 0;
  for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx < std_sortbuf_len; ++new_chr_fo_idx) {
    const uint64_t cur_entry = std_sortbuf[new_chr_fo_idx];
    const uintptr_t chr_idx = S_CAST(uint32_t, cur_entry);
    const uint32_t chr_size = chr_idx_to_size[chr_idx];
    write_cip->chr_file_order[new_chr_fo_idx] = chr_idx;
    write_vidx += chr_size;
    write_cip->chr_fo_vidx_start[new_chr_fo_idx + 1] = write_vidx;
    write_cip->chr_idx_to_foidx[chr_idx] = new_chr_fo_idx;
  }

  const uint32_t new_nonstd_ct = new_chr_ct - std_sortbuf_len;
  if (new_nonstd_ct) {
    StrSortIndexedDeref* nonstd_sort_buf = S_CAST(StrSortIndexedDeref*, bigstack_alloc_raw_rd(new_nonstd_ct * sizeof(StrSortIndexedDeref)));
    if (!nonstd_sort_buf) {
      return 1;
    }
    const char** nonstd_names = cip->nonstd_names;
    uint32_t str_idx = 0;
    for (uint32_t chr_idx = max_code + 1; chr_idx < chr_code_end; ++chr_idx) {
      if (chr_idx_to_size[chr_idx]) {
        nonstd_sort_buf[str_idx].strptr = nonstd_names[chr_idx];
        nonstd_sort_buf[str_idx].orig_idx = chr_idx;
        ++str_idx;
      }
    }
    assert(str_idx == new_nonstd_ct);
    StrptrArrSortMain(new_nonstd_ct, use_nsort, nonstd_sort_buf);
    uint32_t new_chr_fo_idx = std_sortbuf_len;
    for (str_idx = 0; str_idx < new_nonstd_ct; ++str_idx, ++new_chr_fo_idx) {
      const uint32_t chr_idx = nonstd_sort_buf[str_idx].orig_idx;
      const uint32_t chr_size = chr_idx_to_size[chr_idx];
      write_cip->chr_file_order[new_chr_fo_idx] = chr_idx;
      write_vidx += chr_size;
      write_cip->chr_fo_vidx_start[new_chr_fo_idx + 1] = write_vidx;
      write_cip->chr_idx_to_foidx[chr_idx] = new_chr_fo_idx;
    }
  }
  BigstackReset(std_sortbuf);
  return 0;
}

// hybrid of WriteMapOrBim() and write_pvar_resorted()
PglErr WriteBimResorted(const char* outname, const ChrInfo* write_cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t output_zst, uint32_t thread_ct) {
  // allele_dosages must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
      goto WriteBimResorted_ret_NOMEM;
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto WriteBimResorted_ret_1;
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
        char* chr_name_end = chrtoa(write_cip, write_cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
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
      const char* const* cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
        const AltAlleleCt* cur_refalt1_select = &(refalt1_select[variant_uidx * 2]);
        if ((!allele_dosages) || allele_dosages[cur_refalt1_select[1] + variant_allele_idx_base]) {
          cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[1]]);
        } else {
          *cswritep++ = output_missing_geno_char;
        }
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[0]]);
      }
      AppendBinaryEoln(&cswritep);
      if (Cswrite(&css, &cswritep)) {
        goto WriteBimResorted_ret_WRITE_FAIL;
      }
    }
    if (CswriteCloseNull(&css, cswritep)) {
      goto WriteBimResorted_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WriteBimResorted_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteBimResorted_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 WriteBimResorted_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr PvarInfoReloadInterval(const uint32_t* old_variant_uidx_to_new, uintptr_t loadbuf_size, uint32_t variant_idx_start, uint32_t variant_idx_end, uint32_t raw_variant_ct, gzFile gz_pvar_reload, char** pvar_info_strs, char* loadbuf) {
  // We assume the batch size was chosen such that there's no risk of
  // scribbling past g_bigstack_end (barring pathological cases like another
  // process modifying the .pvar file after initial load).
  // We also assume no more dynamic allocations are needed after this;
  // otherwise, str_store_iter should be returned.
  if (gzrewind(gz_pvar_reload)) {
    return kPglRetReadFail;
  }
  const uint32_t cur_batch_size = variant_idx_end - variant_idx_start;
  char* str_store_iter = R_CAST(char*, g_bigstack_base);
  uint32_t info_col_idx;
  PglErr reterr = PvarInfoReloadHeader(loadbuf_size, gz_pvar_reload, loadbuf, &info_col_idx);
  if (reterr) {
    return reterr;
  }
  for (uint32_t variant_uidx = 0; variant_uidx < raw_variant_ct; ++variant_uidx) {
    const char* loadbuf_iter;
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
      loadbuf_iter = FirstNonTspace(loadbuf);
    } while (IsEolnKns(*loadbuf_iter));
    const uint32_t new_variant_idx_offset = old_variant_uidx_to_new[variant_uidx] - variant_idx_start;
    // exploit wraparound, UINT32_MAX null value
    if (new_variant_idx_offset >= cur_batch_size) {
      continue;
    }
    loadbuf_iter = NextTokenMult(loadbuf_iter, info_col_idx);
    if (!loadbuf_iter) {
      return kPglRetReadFail;
    }
    const uint32_t info_slen = strlen_se(loadbuf_iter);
    pvar_info_strs[new_variant_idx_offset] = str_store_iter;
    str_store_iter = memcpyax(str_store_iter, loadbuf_iter, info_slen, '\0');
  }
  assert(str_store_iter <= R_CAST(char*, g_bigstack_end));
  return kPglRetSuccess;
}

// could be BoolErr
PglErr WritePvarResortedInterval(const ChrInfo* write_cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t variant_idx_start, uint32_t variant_idx_end, uint32_t xheader_info_pr, uint32_t write_qual, uint32_t write_filter, uint32_t write_info, uint32_t all_nonref, uint32_t write_cm, char** pvar_info_strs, CompressStreamState* cssp, char** cswritepp, uint32_t* chr_fo_idxp, uint32_t* chr_endp, uint32_t* chr_buf_blenp, char* chr_buf, uintptr_t* allele_include) {
  char* cswritep = *cswritepp;
  uint32_t chr_fo_idx = *chr_fo_idxp;
  uint32_t chr_end = *chr_endp;
  uint32_t chr_buf_blen = *chr_buf_blenp;
  PglErr reterr = kPglRetSuccess;
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
        char* chr_name_end = chrtoa(write_cip, write_cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
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
      const char* const* cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
      if (Cswrite(cssp, &cswritep)) {
        goto WritePvarResortedInterval_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
        SetAllBits(cur_allele_ct, allele_include);
        ClearBit(ref_allele_idx, allele_include);
        ClearBit(alt1_allele_idx, allele_include);
        uint32_t cur_allele_uidx = 0;
        uint32_t alt_allele_idx = 2;
        do {
          *cswritep++ = ',';
          FindFirst1BitFromU32(allele_include, &cur_allele_uidx);
          cswritep = strcpya(cswritep, cur_alleles[cur_allele_uidx++]);
          if (Cswrite(cssp, &cswritep)) {
            goto WritePvarResortedInterval_ret_WRITE_FAIL;
          }
        } while (++alt_allele_idx < cur_allele_ct);
      }

      if (write_qual) {
        *cswritep++ = '\t';
        if (!IsSet(qual_present, variant_uidx)) {
          *cswritep++ = '.';
        } else {
          cswritep = ftoa_g(quals[variant_uidx], cswritep);
        }
      }

      if (write_filter) {
        *cswritep++ = '\t';
        if (!IsSet(filter_present, variant_uidx)) {
          *cswritep++ = '.';
        } else if (!IsSet(filter_npass, variant_uidx)) {
          cswritep = strcpya(cswritep, "PASS");
        } else {
          cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
        }
      }

      if (write_info) {
        *cswritep++ = '\t';
        const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
        if (pvar_info_strs) {
          PvarInfoWrite(xheader_info_pr, is_pr, pvar_info_strs[variant_idx - variant_idx_start], &cswritep);
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
      AppendBinaryEoln(&cswritep);
    }

  }
  while (0) {
  WritePvarResortedInterval_ret_WRITE_FAIL:
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
PglErr WritePvarResorted(const char* outname, const char* xheader, const uintptr_t* variant_include, const ChrInfo* write_cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, PvarPsamFlags pvar_psam_flags, uint32_t thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  gzFile gz_pvar_reload = nullptr;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
        bigstack_alloc_ul(BitCtToWordCt(kPglMaxAltAlleleCt), &allele_include)) {
      goto WritePvarResorted_ret_NOMEM;
    }
    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_flags / kfPvarZs) & 1;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto WritePvarResorted_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_flags & kfPvarColInfo) || pvar_info_reload;
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
      logputs("\n");
      logerrputs("Error: Conflicting INFO:PR fields.  Either fix all REF alleles so that the\n'provisional reference' field is no longer needed, or remove/rename the other\nINFO:PR field.\n");
      goto WritePvarResorted_ret_INCONSISTENT_INPUT;
    }

    if (pvar_psam_flags & kfPvarColXheader) {
      if (CsputsStd(xheader, xheader_blen, &css, &cswritep)) {
        goto WritePvarResorted_ret_WRITE_FAIL;
      }
      if (write_info_pr && (!xheader_info_pr)) {
        cswritep = strcpya(cswritep, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
      }
    }
    if (write_cip->chrset_source) {
      AppendChrsetLine(write_cip, &cswritep);
    }
    cswritep = strcpya(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_flags & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybequal) && qual_present) {
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
    if (pvar_psam_flags & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybefilter) && filter_present) {
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
    if (pvar_psam_flags & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
        // nonzero_cm_present check was performed
        write_cm = 1;
      } else {
        uint32_t variant_uidx = 0;
        for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
          FindFirst1BitFromU32(variant_include, &variant_uidx);
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
    AppendBinaryEoln(&cswritep);

    uint32_t* old_variant_uidx_to_new = nullptr;
    char** pvar_info_strs = nullptr;
    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t batch_size = variant_ct;
    uint32_t batch_ct = 1;
    if (pvar_info_reload) {
      if (bigstack_alloc_ui(raw_variant_ct, &old_variant_uidx_to_new)) {
        goto WritePvarResorted_ret_NOMEM;
      }
      SetAllUiArr(raw_variant_ct, old_variant_uidx_to_new);
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
        const uint32_t old_variant_uidx = new_variant_idx_to_old[variant_idx];
        old_variant_uidx_to_new[old_variant_uidx] = variant_idx;
      }

      reterr = gzopen_read_checked(pvar_info_reload, &gz_pvar_reload);
      if (reterr) {
        return reterr;
      }
      loadbuf_size = RoundDownPow2(bigstack_left() / 4, kCacheline);
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size <= kMaxMediumLine) {
        return kPglRetNomem;
      }
      loadbuf = S_CAST(char*, bigstack_alloc_raw(loadbuf_size));
      loadbuf[loadbuf_size - 1] = ' ';

      // subtract kCacheline to allow for rounding
      uintptr_t bytes_left = bigstack_left() - kCacheline;
      uint32_t single_variant_byte_ct = info_reload_slen + 1 + sizeof(intptr_t);
      if (variant_ct * single_variant_byte_ct > bytes_left) {
        batch_size = bytes_left / single_variant_byte_ct;
        batch_ct = 1 + (variant_ct - 1) / batch_size;
      }
      pvar_info_strs = S_CAST(char**, bigstack_alloc_raw_rd(batch_size * sizeof(intptr_t)));
    }

    uint32_t variant_idx_start = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t batch_idx = 0; batch_idx < batch_ct; ++batch_idx) {
      uint32_t variant_idx_end = MINV(variant_idx_start + batch_size, variant_ct);
      if (gz_pvar_reload) {
        reterr = PvarInfoReloadInterval(old_variant_uidx_to_new, loadbuf_size, variant_idx_start, variant_idx_end, raw_variant_ct, gz_pvar_reload, pvar_info_strs, loadbuf);
        if (reterr) {
          goto WritePvarResorted_ret_1;
        }
      }
      reterr = WritePvarResortedInterval(write_cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, allele_dosages, refalt1_select, qual_present, quals, filter_present, filter_npass, filter_storage, nonref_flags, variant_cms, new_variant_idx_to_old, variant_idx_start, variant_idx_end, xheader_info_pr, write_qual, write_filter, write_info, all_nonref, write_cm, pvar_info_strs, &css, &cswritep, &chr_fo_idx, &chr_end, &chr_buf_blen, chr_buf, allele_include);
      if (reterr) {
        goto WritePvarResorted_ret_1;
      }
      variant_idx_start = variant_idx_end;
    }

    if (CswriteCloseNull(&css, cswritep)) {
      goto WritePvarResorted_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WritePvarResorted_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WritePvarResorted_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WritePvarResorted_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WritePvarResorted_ret_1:
  CswriteCloseCond(&css, cswritep);
  gzclose_cond(gz_pvar_reload);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr MakePlink2Vsort(const char* xheader, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const uint64_t* allele_dosages, const AltAlleleCt* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ChrIdx* chr_idxs, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, uint32_t use_nsort, PvarPsamFlags pvar_psam_flags, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
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
    ChrInfo write_chr_info;

    write_chr_info.haploid_mask = K_CAST(uintptr_t*, cip->haploid_mask);
    write_chr_info.nonstd_names = K_CAST(const char**, cip->nonstd_names);
    write_chr_info.nonstd_id_htable = K_CAST(uint32_t*, cip->nonstd_id_htable);
    write_chr_info.chrset_source = cip->chrset_source;
    memcpy(write_chr_info.chr_exclude, cip->chr_exclude, kChrExcludeWords * sizeof(intptr_t));
    memcpy(write_chr_info.xymt_codes, cip->xymt_codes, kChrOffsetCt * sizeof(int32_t));
    write_chr_info.max_numeric_code = cip->max_numeric_code;
    write_chr_info.max_code = cip->max_code;
    write_chr_info.autosome_ct = cip->autosome_ct;
    write_chr_info.zero_extra_chrs = cip->zero_extra_chrs;
    write_chr_info.name_ct = cip->name_ct;
    write_chr_info.incl_excl_name_stack = K_CAST(LlStr*, cip->incl_excl_name_stack);
    write_chr_info.is_include_stack = cip->is_include_stack;
    write_chr_info.output_encoding = cip->output_encoding;

    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    uint32_t* chr_idx_to_size;
    if (bigstack_calloc_ul(kChrMaskWords, &write_chr_info.chr_mask) ||
        bigstack_alloc_ui(chr_code_end, &write_chr_info.chr_idx_to_foidx) ||
        bigstack_end_calloc_ui(chr_code_end, &chr_idx_to_size)) {
      goto MakePlink2Vsort_ret_NOMEM;
    }
    SetAllUiArr(chr_code_end, write_chr_info.chr_idx_to_foidx);
    if (chr_idxs) {
      uint32_t variant_uidx = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        FindFirst1BitFromU32(variant_include, &variant_uidx);
        chr_idx_to_size[chr_idxs[variant_uidx]] += 1;
      }
      for (uint32_t chr_idx = 0; chr_idx < chr_code_end; ++chr_idx) {
        if (chr_idx_to_size[chr_idx]) {
          SetBit(chr_idx, write_chr_info.chr_mask);
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
        chr_idx_to_size[chr_idx] = PopcountBitRange(variant_include, vidx_start, vidx_end);
        if (chr_idx_to_size[chr_idx]) {
          SetBit(chr_idx, write_chr_info.chr_mask);
        }
        vidx_start = vidx_end;
      }
    }
    if (SortChr(cip, chr_idx_to_size, use_nsort, &write_chr_info)) {
      goto MakePlink2Vsort_ret_NOMEM;
    }

    uint32_t* new_variant_idx_to_old;

    // pos_vidx_sort_buf has variant_bp in high bits, variant_uidx in low
    uint64_t* pos_vidx_sort_buf;
    if (bigstack_alloc_ui(variant_ct, &new_variant_idx_to_old) ||
        bigstack_alloc_ull(variant_ct + 1, &pos_vidx_sort_buf)) {
      goto MakePlink2Vsort_ret_NOMEM;
    }
    pos_vidx_sort_buf[variant_ct] = ~0LLU;
    const uint32_t new_chr_ct = write_chr_info.chr_ct;
    if (chr_idxs) {
      uint32_t* next_write_vidxs;
      if (bigstack_alloc_ui(chr_code_end, &next_write_vidxs)) {
        goto MakePlink2Vsort_ret_NOMEM;
      }
      for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx < new_chr_ct; ++new_chr_fo_idx) {
        const uint32_t chr_idx = write_chr_info.chr_file_order[new_chr_fo_idx];
        next_write_vidxs[chr_idx] = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx];
      }
      uint32_t variant_uidx = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        FindFirst1BitFromU32(variant_include, &variant_uidx);
        const uint32_t chr_idx = chr_idxs[variant_uidx];
        const uint32_t write_vidx = next_write_vidxs[chr_idx];
        pos_vidx_sort_buf[write_vidx] = (S_CAST(uint64_t, variant_bps[variant_uidx]) << 32) | variant_uidx;
        next_write_vidxs[chr_idx] += 1;
      }
      BigstackReset(next_write_vidxs);
    } else {
      uint32_t old_chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t variant_uidx = 0;
      uint32_t chr_idx = 0;
      uint32_t write_vidx = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx, ++write_vidx) {
        FindFirst1BitFromU32(variant_include, &variant_uidx);
        if (variant_uidx >= chr_end) {
          do {
            ++old_chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[old_chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          chr_idx = cip->chr_file_order[old_chr_fo_idx];
          write_vidx = write_chr_info.chr_idx_to_foidx[chr_idx];
        }
        pos_vidx_sort_buf[write_vidx] = (S_CAST(uint64_t, variant_bps[variant_uidx]) << 32) | variant_uidx;
      }
    }

    StrSortIndexedDeref* same_pos_sort_buf = R_CAST(StrSortIndexedDeref*, g_bigstack_base);
    const uintptr_t same_pos_sort_buf_size = bigstack_left() / sizeof(StrSortIndexedDeref);

    uint32_t vidx_start = 0;
    uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
    for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx < new_chr_ct; ++new_chr_fo_idx) {
      const uint32_t vidx_end = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx + 1];
      const uint32_t chr_size = vidx_end - vidx_start;
      const uint64_t post_entry = pos_vidx_sort_buf[vidx_end];
      pos_vidx_sort_buf[vidx_end] = ~0LLU;  // simplify end-of-chromosome logic
      uint64_t* pos_vidx_sort_chr = &(pos_vidx_sort_buf[vidx_start]);

#ifdef __cplusplus
      std::sort(pos_vidx_sort_chr, &(pos_vidx_sort_chr[chr_size]));
#else
      qsort(pos_vidx_sort_chr, chr_size, sizeof(int64_t), uint64cmp);
#endif
      uint32_t prev_pos = pos_vidx_sort_chr[0] >> 32;
      uint32_t prev_variant_uidx = S_CAST(uint32_t, pos_vidx_sort_chr[0]);
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
              goto MakePlink2Vsort_ret_NOMEM;
            }
            const uint32_t variant_uidx = S_CAST(uint32_t, cur_entry);
            same_pos_sort_buf[equal_pos_ct].strptr = variant_ids[variant_uidx];
            same_pos_sort_buf[equal_pos_ct].orig_idx = variant_uidx;
            cur_entry = pos_vidx_sort_chr2[++equal_pos_ct];
            cur_pos = cur_entry >> 32;
          } while (cur_pos == prev_pos);
          StrptrArrSortMain(equal_pos_ct, use_nsort, same_pos_sort_buf);
          for (uint32_t equal_pos_idx = 0; equal_pos_idx < equal_pos_ct; ++equal_pos_idx) {
            *new_variant_idx_to_old_iter++ = same_pos_sort_buf[equal_pos_idx].orig_idx;
          }
          cidx += equal_pos_ct - 1;
        } else {
          *new_variant_idx_to_old_iter++ = prev_variant_uidx;
        }
        prev_pos = cur_pos;
        prev_cidx = cidx;
        prev_variant_uidx = S_CAST(uint32_t, cur_entry);
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
    BigstackReset(pos_vidx_sort_buf);

    const uint32_t trim_alts = S_CAST(uint32_t, make_plink2_flags & kfMakePlink2TrimAlts);
    if (make_plink2_flags & kfMakeBim) {
      const uint32_t bim_zst = (make_plink2_flags / kfMakeBimZs) & 1;
      OutnameZstSet(".bim", bim_zst, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);

      reterr = WriteBimResorted(outname, &write_chr_info, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, variant_cms, new_variant_idx_to_old, variant_ct, max_allele_slen, bim_zst, max_thread_ct);
      if (reterr) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePvar) {
      OutnameZstSet(".pvar", pvar_psam_flags & kfPvarZs, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t nonref_flags_storage = 3;
      if (!simple_pgrp->fi.nonref_flags) {
        nonref_flags_storage = (simple_pgrp->fi.gflags & kfPgenGlobalAllNonref)? 2 : 1;
      }
      reterr = WritePvarResorted(outname, xheader, variant_include, &write_chr_info, variant_bps, variant_ids, variant_allele_idxs, allele_storage, trim_alts? allele_dosages : nullptr, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, simple_pgrp->fi.nonref_flags, pvar_info_reload, variant_cms, new_variant_idx_to_old, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, xheader_info_pr, xheader_info_pr_nonflag, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_flags, max_thread_ct);
      if (reterr) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakeFam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, new_sample_idx_to_old, sample_ct, pheno_ct, '\t');
      if (reterr) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePsam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WritePsam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, sample_ct, pheno_ct, max_pheno_name_blen, pvar_psam_flags);
      if (reterr) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
      // boilerplate from start of MakePlink2NoVsort()
      if (make_plink2_flags & kfMakePlink2MMask) {
        logerrputs("Error: --make-bed/--make-{b}pgen multiallelics= is currently under development.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePlink2Vsort_ret_1;
      }
      g_plink2_write_flags = kfPlink2Write0;
      const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      if (make_plink2_flags & kfMakePlink2SetHhMissing) {
        const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
        uintptr_t* new_sex_male;
        uintptr_t* sex_female;
        if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &new_sex_male) ||
            bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
            bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed) ||
            bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_female_collapsed_interleaved) ||
            bigstack_alloc_ul(raw_sample_ctl, &sex_female)) {
          goto MakePlink2Vsort_ret_NOMEM;
        }
        CopyBitarrSubset(sex_male, sample_include, sample_ct, new_sex_male);
        ZeroTrailingWords(BitCtToWordCt(sample_ct), new_sex_male);
        g_sex_male = new_sex_male;
        FillInterleavedMaskVec(g_sex_male, sample_ctv, g_sex_male_collapsed_interleaved);

        BitvecAndNotCopy(sex_nm, sex_male, raw_sample_ctl, sex_female);
        CopyBitarrSubset(sex_female, sample_include, sample_ct, g_sex_female_collapsed);
        FillInterleavedMaskVec(g_sex_female_collapsed, sample_ctv, g_sex_female_collapsed_interleaved);

        BigstackReset(g_sex_female_collapsed);
        g_plink2_write_flags |= kfPlink2WriteSetHhMissing;
      }
      if (make_plink2_flags & kfMakePlink2SetMixedMtMissing) {
        g_plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
      }
      g_cip = &write_chr_info;
      reterr = MakePgenRobust(sample_include, new_sample_idx_to_old, variant_include, variant_allele_idxs, refalt1_select, new_variant_idx_to_old, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, hard_call_thresh, dosage_erase_thresh, make_plink2_flags, simple_pgrp, outname, outname_end);
      if (reterr) {
        goto MakePlink2Vsort_ret_1;
      }
    }
  }
  while (0) {
  MakePlink2Vsort_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 MakePlink2Vsort_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

PglErr SampleSortFileMap(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* sample_sort_fname, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t** new_sample_idx_to_old_ptr) {
  // assumes sample_ct >= 2 (enforced by caller)
  // return strbox is not collapsed
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  {
    char* idbuf;
    uintptr_t* already_seen;
    if (bigstack_alloc_ui(raw_sample_ct, new_sample_idx_to_old_ptr) ||
        bigstack_alloc_c(siip->max_sample_id_blen, &idbuf) ||
        bigstack_calloc_ul(BitCtToWordCt(raw_sample_ct), &already_seen)) {
      goto SampleSortFileMap_ret_NOMEM;
    }

    uintptr_t loadbuf_size = bigstack_left();
    loadbuf_size -= loadbuf_size / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto SampleSortFileMap_ret_NOMEM;
    } else {
      loadbuf_size = RoundUpPow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = S_CAST(char*, bigstack_alloc_raw(loadbuf_size));
    char* loadbuf_first_token;
    XidMode xid_mode;
    reterr = OpenAndLoadXidHeader(sample_sort_fname, "indiv-sort", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, loadbuf_size, loadbuf, nullptr, &line_idx, &loadbuf_first_token, &gz_infile, &xid_mode);
    if (reterr) {
      if (reterr == kPglRetEmptyFile) {
        logerrputs("Error: --indiv-sort file is empty.\n");
        goto SampleSortFileMap_ret_MALFORMED_INPUT;
      }
      if (reterr == kPglRetLongLine) {
        if (loadbuf_size == kMaxLongLine) {
          goto SampleSortFileMap_ret_LONG_LINE;
        }
        goto SampleSortFileMap_ret_NOMEM;
      }
      goto SampleSortFileMap_ret_1;
    }
    uint32_t* xid_map;
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (reterr) {
      goto SampleSortFileMap_ret_1;
    }
    uint32_t* new_sample_idx_to_old_iter = *new_sample_idx_to_old_ptr;
    if (*loadbuf_first_token == '#') {
      *loadbuf_first_token = '\0';
    }
    while (1) {
      if (!IsEolnKns(*loadbuf_first_token)) {
        const char* loadbuf_iter = loadbuf_first_token;
        uint32_t sample_uidx;
        if (!SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, sample_ct, 0, xid_mode, &loadbuf_iter, &sample_uidx, idbuf)) {
          if (IsSet(already_seen, sample_uidx)) {
            char* tab_iter = S_CAST(char*, rawmemchr(idbuf, '\t'));
            *tab_iter = ' ';
            if (xid_mode & kfXidModeFlagSid) {
              *S_CAST(char*, rawmemchr(&(tab_iter[1]), '\t')) = ' ';
            }
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID '%s' in --indiv-sort file.\n", idbuf);
            goto SampleSortFileMap_ret_MALFORMED_INPUT_WW;
          }
          SetBit(sample_uidx, already_seen);
          *new_sample_idx_to_old_iter++ = sample_uidx;
        } else if (!loadbuf_iter) {
          goto SampleSortFileMap_ret_MISSING_TOKENS;
        }
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto SampleSortFileMap_ret_READ_FAIL;
        }
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        if (loadbuf_size == kMaxLongLine) {
          goto SampleSortFileMap_ret_LONG_LINE;
        }
        goto SampleSortFileMap_ret_NOMEM;
      }
      loadbuf_first_token = FirstNonTspace(loadbuf);
      if (loadbuf_first_token[0] == '#') {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --indiv-sort file starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx);
        goto SampleSortFileMap_ret_MALFORMED_INPUT_WW;
      }
    }

    if (gzclose_null(&gz_infile)) {
      goto SampleSortFileMap_ret_READ_FAIL;
    }
    if (S_CAST(uintptr_t, new_sample_idx_to_old_iter - (*new_sample_idx_to_old_ptr)) != sample_ct) {
      logerrputs("Error: --indiv-sort file does not contain all loaded sample IDs.\n");
      goto SampleSortFileMap_ret_INCONSISTENT_INPUT;
    }
    bigstack_mark = R_CAST(unsigned char*, idbuf);
  }
  while (0) {
  SampleSortFileMap_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SampleSortFileMap_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  SampleSortFileMap_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  SampleSortFileMap_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  SampleSortFileMap_ret_LONG_LINE:
    logerrprintf("Error: Line %" PRIuPTR " of --indiv-sort file is pathologically long.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  SampleSortFileMap_ret_MISSING_TOKENS:
    logerrprintf("Error: Line %" PRIuPTR " of --indiv-sort file has fewer tokens than expected.\n", line_idx);
  SampleSortFileMap_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 SampleSortFileMap_ret_1:
  gzclose_cond(gz_infile);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
