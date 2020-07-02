// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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


#include "include/pgenlib_write.h"
#include "plink2_compress_stream.h"
#include "plink2_data.h"
#include "plink2_pvar.h"

#include <time.h>

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr WriteMapOrBim(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const double* variant_cms, uint32_t variant_ct, uint32_t max_allele_slen, char delim, uint32_t output_zst, uint32_t thread_ct) {
  // - Normally generates a .bim file.  Set max_allele_slen to zero to generate
  //   a .map.
  // - allele_presents must be nullptr unless we're trimming alt alleles.
  // - Errors out when writing .bim if any remaining variant is multiallelic
  //   and refalt1_select is nullptr.
  // - Multiallelic-split case is handled by WriteBimSplit().
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WriteMapOrBim_ret_NOMEM;
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WriteMapOrBim_ret_1;
    }

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
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
      cswritep = u32toa(variant_bps[variant_uidx], cswritep);
      if (max_allele_slen) {
        *cswritep++ = delim;
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          if (!refalt1_select) {
            const uintptr_t allele_idx_offset_end = allele_idx_offsets[variant_uidx + 1];
            if (allele_idx_offset_end != allele_idx_offset_base + 2) {
              // not actually unlikely at this point, but simplest to stay
              // consistent
              if (unlikely((!allele_presents) || (!AllBitsAreZero(allele_presents, 2 + allele_idx_offset_base, allele_idx_offset_end)))) {
                logputs("\n");
                logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
                goto WriteMapOrBim_ret_INCONSISTENT_INPUT;
              }
            }
          }
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        // note that VCF ref allele corresponds to A2, not A1
        if (!refalt1_select) {
          if ((!allele_presents) || IsSet(allele_presents, 1 + allele_idx_offset_base)) {
            cswritep = strcpya(cswritep, cur_alleles[1]);
          } else {
            *cswritep++ = output_missing_geno_char;
          }
          *cswritep++ = delim;
          cswritep = strcpya(cswritep, cur_alleles[0]);
        } else {
          STD_ARRAY_KREF(AlleleCode, 2) cur_refalt1_select = refalt1_select[variant_uidx];
          if ((!allele_presents) || IsSet(allele_presents, cur_refalt1_select[1] + allele_idx_offset_base)) {
            cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[1]]);
          } else {
            *cswritep++ = output_missing_geno_char;
          }
          *cswritep++ = delim;
          cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[0]]);
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto WriteMapOrBim_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
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
  WriteMapOrBim_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WriteMapOrBim_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr PvarInfoReloadHeader(TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* info_col_idx_ptr) {
  char* line_iter;
  do {
    PglErr reterr = TextNextLineLstrip(pvar_reload_txsp, &line_iter);
    if (unlikely(reterr)) {
      return reterr;
    }
  } while (!StrStartsWithUnsafe(line_iter, "#CHROM"));
  uint32_t info_col_idx = 0;
  do {
    line_iter = NextToken(line_iter);
    ++info_col_idx;
  } while (!tokequal_k(line_iter, "INFO"));
  *line_iterp = line_iter;
  *info_col_idx_ptr = info_col_idx;
  return kPglRetSuccess;
}

// May use all remaining workspace memory.
PglErr PvarInfoOpenAndReloadHeader(const char* pvar_info_reload, uint32_t calc_thread_ct, TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* info_col_idx_ptr) {
  PglErr reterr = SizeAndInitTextStream(pvar_info_reload, bigstack_left(), calc_thread_ct, pvar_reload_txsp);
  if (unlikely(reterr)) {
    return reterr;
  }
  return PvarInfoReloadHeader(pvar_reload_txsp, line_iterp, info_col_idx_ptr);
}

void PvarInfoWrite(uint32_t info_pr_flag_present, uint32_t is_pr, char* info_token, char** write_iter_ptr) {
  char* info_token_end = CurTokenEnd(info_token);
  uint32_t info_token_slen = info_token_end - info_token;
  char* info_token_pr = nullptr;
  if (info_pr_flag_present) {
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
        write_iter = strcpya_k(write_iter, ";PR");
      }
    }
  } else {
    // possible with --real-ref-alleles/--ref-from-fa
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

PglErr PvarInfoReload(uint32_t info_col_idx, uint32_t variant_uidx, TextStream* pvar_reload_txsp, char** line_iterp, uint32_t* trs_variant_uidx_ptr) {
  uint32_t trs_variant_uidx = *trs_variant_uidx_ptr;
  char* line_iter = AdvPastDelim(*line_iterp, '\n');
  if (trs_variant_uidx < variant_uidx) {
    TextSetPos(line_iter, pvar_reload_txsp);
    PglErr reterr = TextSkipNz(variant_uidx - trs_variant_uidx, pvar_reload_txsp);
    if (unlikely(reterr)) {
      return reterr;
    }
    line_iter = TextLineEnd(pvar_reload_txsp);
    trs_variant_uidx = variant_uidx;
  }
  PglErr reterr = TextNextLineLstripUnsafe(pvar_reload_txsp, &line_iter);
  if (unlikely(reterr)) {
    return reterr;
  }
  *line_iterp = NextTokenMultFar(line_iter, info_col_idx);

  // index *after* just-loaded line.
  *trs_variant_uidx_ptr = trs_variant_uidx + 1;
  return kPglRetSuccess;
}

PglErr PvarInfoReloadAndWrite(uint32_t info_pr_flag_present, uint32_t info_col_idx, uint32_t variant_uidx, uint32_t is_pr, TextStream* pvar_reload_txsp, char** line_iterp, char** write_iter_ptr, uint32_t* trs_variant_uidx_ptr) {
  PglErr reterr = PvarInfoReload(info_col_idx, variant_uidx, pvar_reload_txsp, line_iterp, trs_variant_uidx_ptr);
  if (unlikely(reterr)) {
    return reterr;
  }
  PvarInfoWrite(info_pr_flag_present, is_pr, *line_iterp, write_iter_ptr);
  return kPglRetSuccess;
}

void AppendChrsetLine(const ChrInfo* cip, char** write_iter_ptr) {
  char* write_iter = strcpya_k(*write_iter_ptr, "##chrSet=<");
  if (!(cip->haploid_mask[0] & 1)) {
    write_iter = strcpya_k(write_iter, "autosomePairCt=");
    write_iter = u32toa(cip->autosome_ct, write_iter);
    if (!IsI32Neg(cip->xymt_codes[kChrOffsetX])) {
      write_iter = strcpya_k(write_iter, ",X");
    }
    if (!IsI32Neg(cip->xymt_codes[kChrOffsetY])) {
      write_iter = strcpya_k(write_iter, ",Y");
    }
    if (!IsI32Neg(cip->xymt_codes[kChrOffsetXY])) {
      write_iter = strcpya_k(write_iter, ",XY");
    }
    if (!IsI32Neg(cip->xymt_codes[kChrOffsetMT])) {
      write_iter = strcpya_k(write_iter, ",M");
    }
    if (!IsI32Neg(cip->xymt_codes[kChrOffsetPAR1])) {
      write_iter = strcpya_k(write_iter, ",PAR1");
    }
    if (!IsI32Neg(cip->xymt_codes[kChrOffsetPAR2])) {
      write_iter = strcpya_k(write_iter, ",PAR2");
    }
  } else {
    write_iter = strcpya_k(write_iter, "haploidAutosomeCt=");
    write_iter = u32toa(cip->autosome_ct, write_iter);
  }
  *write_iter++ = '>';
  *write_iter_ptr = write_iter;
  AppendBinaryEoln(write_iter_ptr);
}

// fileformat, fileDate, source
void AppendVcfHeaderStart(uint32_t v43, char** cswritepp) {
  char* cswritep = *cswritepp;
  cswritep = strcpya_k(cswritep, "##fileformat=VCFv4.");
  *cswritep++ = v43 + '2';
  cswritep = strcpya_k(cswritep, EOLN_STR "##fileDate=");
  time_t rawtime;
  time(&rawtime);
  const struct tm* loctime = localtime(&rawtime);
  cswritep += strftime(cswritep, kMaxMediumLine, "%Y%m%d", loctime);
  cswritep = strcpya_k(cswritep, EOLN_STR "##source=PLINKv2.00" EOLN_STR);
  *cswritepp = cswritep;
  return;
}

// Note that the order-of-operations page lists this as happening right after
// the filtering performed by LoadPvar().  Which is effectively true, since we
// ignore variant_include (this is safe since LoadPvar() always initializes
// all variant_bps[] and allele_storage[] entries appropriately).
// possible todo: ChrInfo can have a length field, which is initialized by the
// ##contig header line when possible, but when that doesn't exist LoadPvar()
// can conditionally detect INFO:END and take that into account.  (Or a reason
// to keep the entire info_end array in memory may emerge.)
uint32_t ChrLenLbound(const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* new_variant_idx_to_old, uint32_t chr_fo_idx, uint32_t max_allele_slen, UnsortedVar vpos_sortstatus) {
  const uint32_t vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  assert(vidx_start != vidx_end);
  if (!(vpos_sortstatus & kfUnsortedVarBp)) {
    if (!new_variant_idx_to_old) {
      if (max_allele_slen == 1) {
        return variant_bps[vidx_end - 1];
      }
      uint32_t bp_end = 0;
      for (uint32_t vidx = vidx_end; vidx != vidx_start; ) {
        --vidx;
        const uint32_t cur_bp = variant_bps[vidx];
        if (cur_bp + max_allele_slen <= bp_end) {
          break;
        }
        uintptr_t allele_idx_offset_base = vidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[vidx];
        }
        // We only care about reference-allele length.
        const uint32_t cur_bp_end = cur_bp + strlen(allele_storage[allele_idx_offset_base]) - 1;
        if (cur_bp_end > bp_end) {
          bp_end = cur_bp_end;
        }
      }
      return bp_end;
    }
    if (max_allele_slen == 1) {
      return variant_bps[new_variant_idx_to_old[vidx_end - 1]];
    }
    uint32_t bp_end = 0;
    for (uint32_t new_vidx = vidx_end; new_vidx != vidx_start; ) {
      --new_vidx;
      const uint32_t old_vidx = new_variant_idx_to_old[new_vidx];
      const uint32_t cur_bp = variant_bps[old_vidx];
      if (cur_bp + max_allele_slen <= bp_end) {
        break;
      }
      uintptr_t allele_idx_offset_base = old_vidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[old_vidx];
      }
      const uint32_t cur_bp_end = cur_bp + strlen(allele_storage[allele_idx_offset_base]) - 1;
      if (cur_bp_end > bp_end) {
        bp_end = cur_bp_end;
      }
    }
    return bp_end;
  }
  uint32_t bp_end = U32ArrMax(&(variant_bps[vidx_start]), vidx_end - vidx_start);
  if (max_allele_slen == 1) {
    return bp_end;
  }
  uint32_t min_check_bp = 0;
  if (bp_end >= max_allele_slen) {
    min_check_bp = bp_end + 1 - max_allele_slen;
  }
  for (uint32_t vidx = vidx_start; vidx != vidx_end; ++vidx) {
    const uint32_t cur_bp = variant_bps[vidx];
    if (cur_bp < min_check_bp) {
      continue;
    }
    uintptr_t allele_idx_offset_base = vidx * 2;
    if (allele_idx_offsets) {
      allele_idx_offset_base = allele_idx_offsets[vidx];
    }
    const uint32_t cur_bp_end = cur_bp + strlen(allele_storage[allele_idx_offset_base]) - 1;
    if (cur_bp_end > bp_end) {
      bp_end = cur_bp_end;
    }
  }
  return bp_end;
}

PglErr PvarXheaderWrite(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* new_variant_idx_to_old, uintptr_t xheader_blen, uint32_t vcfheader, uint32_t write_filter, uint32_t write_info, uint32_t append_info_pr_header_line, uint32_t max_allele_slen, UnsortedVar vpos_sortstatus, char* xheader, CompressStreamState* css_ptr, char** cswritepp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!vcfheader) {
      if (write_filter && write_info) {
        if (unlikely(CsputsStd(xheader, xheader_blen, css_ptr, cswritepp))) {
          goto PvarXheaderWrite_ret_WRITE_FAIL;
        }
      } else {
        // Filter out FILTER/INFO definitions iff the corresponding column has
        // been removed.
        const char* copy_start = xheader;
        const char* xheader_end = &(xheader[xheader_blen]);
        for (const char* xheader_iter = xheader; xheader_iter != xheader_end; ) {
          const char* next_line_start = AdvPastDelim(xheader_iter, '\n');
          if (((!write_filter) && StrStartsWithUnsafe(xheader_iter, "##FILTER=<ID=")) ||
              ((!write_info) && StrStartsWithUnsafe(xheader_iter, "##INFO=<ID="))) {
            if (copy_start != xheader_iter) {
              if (unlikely(CsputsStd(copy_start, xheader_iter - copy_start, css_ptr, cswritepp))) {
                goto PvarXheaderWrite_ret_WRITE_FAIL;
              }
            }
            copy_start = next_line_start;
          }
          xheader_iter = next_line_start;
        }
        if (copy_start != xheader_end) {
          if (unlikely(CsputsStd(copy_start, xheader_end - copy_start, css_ptr, cswritepp))) {
            goto PvarXheaderWrite_ret_WRITE_FAIL;
          }
        }
      }
    } else {
      // See the start of ExportVcf().
      AppendVcfHeaderStart(1, cswritepp);
      const uint32_t chr_ctl = BitCtToWordCt(cip->chr_ct);
      uintptr_t* written_contig_header_lines;
      if (unlikely(bigstack_calloc_w(chr_ctl, &written_contig_header_lines))) {
        goto PvarXheaderWrite_ret_NOMEM;
      }
      uint32_t contig_zero_written = 0;
      char* cswritep = *cswritepp;
      // ExportVcf() has to perform a customized --merge-par operation, so it
      // has special handling of chrX/PAR1/PAR2 ##contig header lines.  We omit
      // that here.
      char* xheader_end = &(xheader[xheader_blen]);
      for (char* line_end = xheader; line_end != xheader_end; ) {
        char* line_start = line_end;
        line_end = AdvPastDelim(line_start, '\n');
        const uint32_t slen = line_end - line_start;
        if ((slen > 14) && StrStartsWithUnsafe(line_start, "##contig=<ID=")) {
          char* contig_name_start = &(line_start[13]);
          char* contig_name_end = S_CAST(char*, memchr(contig_name_start, ',', slen - 14));
          if (!contig_name_end) {
            // if this line is technically well-formed (ends in '>'), it's
            // useless anyway, throw it out
            continue;
          }
          const uint32_t chr_idx = GetChrCodeCounted(cip, contig_name_end - contig_name_start, contig_name_start);
          if (IsI32Neg(chr_idx) || (!IsSet(cip->chr_mask, chr_idx))) {
            continue;
          }
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
          if (unlikely(IsSet(written_contig_header_lines, chr_fo_idx))) {
            logerrputs("Error: Duplicate ##contig line in .pvar file.\n");
            goto PvarXheaderWrite_ret_MALFORMED_INPUT;
          }
          SetBit(chr_fo_idx, written_contig_header_lines);
          // if --output-chr was used at some point, we need to sync the
          // ##contig chromosome code with the code in the .pvar body.
          char* chr_name_write_start = strcpya_k(cswritep, "##contig=<ID=");
          char* chr_name_write_end = chrtoa(cip, chr_idx, chr_name_write_start);
          if ((*chr_name_write_start == '0') && (chr_name_write_end == &(chr_name_write_start[1]))) {
            // --allow-extra-chr 0 special case
            // note that cswritep has *not* been advanced
            contig_zero_written = 1;  // technically we write this a bit later
            continue;
          }
          cswritep = chr_name_write_end;
          if (unlikely(Cswrite(css_ptr, &cswritep))) {
            goto PvarXheaderWrite_ret_WRITE_FAIL;
          }
          if (unlikely(CsputsStd(contig_name_end, line_end - contig_name_end, css_ptr, &cswritep))) {
            goto PvarXheaderWrite_ret_WRITE_FAIL;
          }
        } else {
          if (!write_filter) {
            if (StrStartsWithUnsafe(line_start, "##FILTER=<ID=")) {
              continue;
            }
          }
          if (!write_info) {
            if (StrStartsWithUnsafe(line_start, "##INFO=<ID=")) {
              continue;
            }
          }
          if (unlikely(CsputsStd(line_start, slen, css_ptr, &cswritep))) {
            goto PvarXheaderWrite_ret_WRITE_FAIL;
          }
        }
      }
      // fill in the missing ##contig lines
      if (contig_zero_written) {
        cswritep = strcpya_k(cswritep, "##contig=<ID=0,length=2147483645>" EOLN_STR);
      }
      for (uint32_t chr_fo_idx = 0; chr_fo_idx != cip->chr_ct; ++chr_fo_idx) {
        if (IsSet(written_contig_header_lines, chr_fo_idx)) {
          continue;
        }
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        // AllBitsAreZero() doesn't do what we want in the --sort-vars case,
        // but fortunately we don't need it there.
        if ((!IsSet(cip->chr_mask, chr_idx)) || (variant_include && AllBitsAreZero(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1]))) {
          continue;
        }
        char* chr_name_write_start = strcpya_k(cswritep, "##contig=<ID=");
        char* chr_name_write_end = chrtoa(cip, chr_idx, chr_name_write_start);
        if ((*chr_name_write_start == '0') && (chr_name_write_end == &(chr_name_write_start[1]))) {
          // --allow-extra-chr 0 special case
          if (contig_zero_written) {
            continue;
          }
          contig_zero_written = 1;
          cswritep = strcpya_k(chr_name_write_end, ",length=2147483645");
        } else {
          cswritep = strcpya_k(chr_name_write_end, ",length=");
          const uint32_t pos_end = ChrLenLbound(cip, variant_bps, allele_idx_offsets, allele_storage, new_variant_idx_to_old, chr_fo_idx, max_allele_slen, vpos_sortstatus);
          cswritep = u32toa(pos_end, cswritep);
        }
        *cswritep++ = '>';
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(css_ptr, &cswritep))) {
          goto PvarXheaderWrite_ret_WRITE_FAIL;
        }
      }
      *cswritepp = cswritep;
    }
    if (append_info_pr_header_line) {
      *cswritepp = strcpya_k(*cswritepp, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
    }
  }
  while (0) {
  PvarXheaderWrite_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PvarXheaderWrite_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  PvarXheaderWrite_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WritePvar(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, UnsortedVar vpos_sortstatus, PvarPsamFlags pvar_psam_flags, uint32_t thread_ct, char* xheader) {
  // allele_presents must be nullptr unless we're trimming alt alleles
  // split/join cases handled by WritePvarSplit() and WritePvarJoin()
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  TextStream pvar_reload_txs;
  PreinitCstream(&css);
  PreinitTextStream(&pvar_reload_txs);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WritePvar_ret_NOMEM;
    }
    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_flags / kfPvarZs) & 1;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WritePvar_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_flags & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      write_info_pr = !IntersectionIsEmpty(variant_include, nonref_flags, raw_variant_ctl);
    }
    write_info_pr = write_info_pr && write_info;
    if (unlikely(write_info_pr && (info_flags & kfInfoPrNonflagPresent))) {
      logputs("\n");
      logerrputs("Error: Conflicting INFO:PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO:PR key.\n");
      goto WritePvar_ret_INCONSISTENT_INPUT;
    }

    uint32_t write_filter = 0;
    if (pvar_psam_flags & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybefilter) && filter_present) {
      write_filter = !IntersectionIsEmpty(variant_include, filter_present, raw_variant_ctl);
    }
    char* pvar_info_line_iter = nullptr;
    uint32_t info_col_idx = 0;  // could save this during first load instead
    const uint32_t info_pr_flag_present = (info_flags / kfInfoPrFlagPresent) & 1;
    if (pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader)) {
      reterr = PvarXheaderWrite(variant_include, cip, variant_bps, allele_idx_offsets, allele_storage, nullptr, xheader_blen, (pvar_psam_flags / kfPvarColVcfheader) & 1, write_filter, write_info, write_info_pr && (!info_pr_flag_present), max_allele_slen, vpos_sortstatus, xheader, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WritePvar_ret_1;
      }
    }
    // bugfix (30 Jul 2017): may be necessary to reload INFO when no ## lines
    // are in the header... er, should we still allow this?
    if (pvar_info_reload) {
      reterr = PvarInfoOpenAndReloadHeader(pvar_info_reload, 1 + (thread_ct > 1), &pvar_reload_txs, &pvar_info_line_iter, &info_col_idx);
      if (unlikely(reterr)) {
        goto WritePvar_ret_TSTREAM_FAIL;
      }
    }
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &cswritep);
    }
    cswritep = strcpya_k(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_flags & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybequal) && qual_present) {
      write_qual = !IntersectionIsEmpty(variant_include, qual_present, raw_variant_ctl);
    }
    if (write_qual) {
      cswritep = strcpya_k(cswritep, "\tQUAL");
    }

    if (write_filter) {
      cswritep = strcpya_k(cswritep, "\tFILTER");
    }

    if (write_info) {
      cswritep = strcpya_k(cswritep, "\tINFO");
    }

    uint32_t write_cm = 0;
    if (pvar_psam_flags & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
        // nonzero_cm_present check was performed
        write_cm = 1;
      } else {
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          if (variant_cms[variant_uidx] != 0.0) {
            write_cm = 1;
            break;
          }
        }
      }
    }
    if (write_cm) {
      cswritep = strcpya_k(cswritep, "\tCM");
    }
    AppendBinaryEoln(&cswritep);

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t trs_variant_uidx = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t cur_allele_ct = 2;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
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
      cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx][0];
        alt1_allele_idx = refalt1_select[variant_uidx][1];
      }
      cswritep = strcpyax(cswritep, cur_alleles[ref_allele_idx], '\t');
      uint32_t alt_allele_written = 0;
      if ((!allele_presents) || IsSet(allele_presents, allele_idx_offset_base + alt1_allele_idx)) {
        cswritep = strcpya(cswritep, cur_alleles[alt1_allele_idx]);
        alt_allele_written = 1;
      }
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto WritePvar_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if ((allele_idx == ref_allele_idx) || (allele_idx == alt1_allele_idx) || (allele_presents && (!IsSet(allele_presents, allele_idx_offset_base + allele_idx)))) {
            continue;
          }
          if (alt_allele_written) {
            *cswritep++ = ',';
          }
          alt_allele_written = 1;
          cswritep = strcpya(cswritep, cur_alleles[allele_idx]);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WritePvar_ret_WRITE_FAIL;
          }
        }
      }
      if (!alt_allele_written) {
        *cswritep++ = output_missing_geno_char;
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
          cswritep = strcpya_k(cswritep, "PASS");
        } else {
          cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
        }
      }

      if (write_info) {
        *cswritep++ = '\t';
        const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
        if (pvar_info_line_iter) {
          reterr = PvarInfoReloadAndWrite(info_pr_flag_present, info_col_idx, variant_uidx, is_pr, &pvar_reload_txs, &pvar_info_line_iter, &cswritep, &trs_variant_uidx);
          if (unlikely(reterr)) {
            goto WritePvar_ret_TSTREAM_FAIL;
          }
        } else {
          if (is_pr) {
            cswritep = strcpya_k(cswritep, "PR");
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
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WritePvar_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  WritePvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WritePvar_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_reload_txs);
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
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteFam(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, char delim) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
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
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    uint32_t sample_uidx2 = 0;
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      uintptr_t sample_uidx;
      if (!new_sample_idx_to_old) {
        sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      } else {
        do {
          sample_uidx = new_sample_idx_to_old[sample_uidx2++];
        } while (!IsSet(sample_include, sample_uidx));
      }
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      if (delim == '\t') {
        write_iter = strcpya(write_iter, cur_sample_id);
      } else {
        const char* fid_end = AdvToDelim(cur_sample_id, '\t');
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
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto WriteFam_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
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

uint32_t DataFidColIsRequired(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t maybe_modifier) {
  if (maybe_modifier & 2) {
    return 1;
  }
  if ((!(maybe_modifier & 1)) || (!(siip->flags & kfSampleIdFidPresent))) {
    return 0;
  }
  const char* sample_ids = siip->sample_ids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    if (!memequal_k(&(sample_ids[sample_uidx * max_sample_id_blen]), "0\t", 2)) {
      return 1;
    }
  }
  return 0;
}

uint32_t DataSidColIsRequired(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier) {
  // note that MAYBESID and SID can both be set
  if (maybe_modifier & 2) {
    return 1;
  }
  if (sids && (maybe_modifier & 1)) {
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      if (!memequal_k(&(sids[sample_uidx * max_sid_blen]), "0", 2)) {
        return 1;
      }
    }
  }
  return 0;
}

uint32_t DataParentalColsAreRequired(const uintptr_t* sample_include, const PedigreeIdInfo* piip, uint32_t sample_ct, uint32_t maybe_modifier) {
  if (maybe_modifier & 2) {
    return 1;
  }
  if ((!(maybe_modifier & 1)) || (!(piip->sii.flags & kfSampleIdParentsPresent))) {
    return 0;
  }
  const char* paternal_ids = piip->parental_id_info.paternal_ids;
  const char* maternal_ids = piip->parental_id_info.maternal_ids;
  const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
  const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
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
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
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
    const uint32_t write_fid = DataFidColIsRequired(sample_include, &(piip->sii), sample_ct, pvar_psam_flags / kfPsamColMaybefid);
    const uint32_t write_sid = DataSidColIsRequired(sample_include, sids, sample_ct, max_sid_blen, pvar_psam_flags / kfPsamColMaybesid);
    const uint32_t write_parents = DataParentalColsAreRequired(sample_include, piip, sample_ct, pvar_psam_flags / kfPsamColMaybeparents);
    const uint32_t write_sex = (pvar_psam_flags / kfPsamColSex) & 1;
    const uint32_t write_empty_pheno = (pvar_psam_flags & kfPsamColPheno1) && (!pheno_ct);
    const uint32_t write_phenos = (pvar_psam_flags & (kfPsamColPheno1 | kfPsamColPhenos)) && pheno_ct;
    if (write_phenos && (!(pvar_psam_flags & kfPsamColPhenos))) {
      pheno_ct = 1;
    }
    char* write_iter = g_textbuf;
    *write_iter++ = '#';
    if (write_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    if (write_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    if (write_parents) {
      write_iter = strcpya_k(write_iter, "\tPAT\tMAT");
    }
    if (write_sex) {
      write_iter = strcpya_k(write_iter, "\tSEX");
    }
    if (write_phenos) {
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        *write_iter++ = '\t';
        const char* cur_pheno_name = &(pheno_names[pheno_idx * max_pheno_name_blen]);
        const uint32_t cur_pheno_name_slen = strlen(cur_pheno_name);
        if (strequal_k(cur_pheno_name, "SEX", cur_pheno_name_slen)) {
          if (unlikely(write_sex)) {
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
            uintptr_t sample_uidx_base = 0;
            uintptr_t cur_bits = sample_include[0];
            if (sex_col->type_code == kPhenoDtypeQt) {
              const double* pheno_vals = sex_col->data.qt;
              for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
                if (IsSet(pheno_nm, sample_uidx)) {
                  const double dxx = pheno_vals[sample_uidx];
                  // tolerate '-9' and '0' as missing values, and anything in
                  // [1, 2] (could be reasonable to represent XXY, etc. with
                  // decimals).
                  if (unlikely(((dxx < 1.0) && (dxx != -9.0) && (dxx != 0.0)) || (dxx > 2.0))) {
                    logerrputs("Error: .psam numeric SEX values are expected to be in {-9, 0, 1, 2}.\n");
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
                  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
                    if (IsSet(pheno_nm, sample_uidx)) {
                      const uint32_t cur_cat_idx = pheno_vals[sample_uidx];
                      if (unlikely((cur_cat_idx != male_cat_idx1) && (cur_cat_idx != female_cat_idx1) && (cur_cat_idx != male_cat_idx2) && (cur_cat_idx != female_cat_idx2))) {
                        logerrputs("Error: .psam alphabetic SEX values are expected to be in {'F', 'f', 'M', 'm'}.\n");
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
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto WritePsam_ret_WRITE_FAIL;
        }
      }
    } else if (write_empty_pheno) {
      write_iter = strcpya_k(write_iter, "\tPHENO1");
    }
    AppendBinaryEoln(&write_iter);

    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    uint32_t sample_uidx2 = 0;
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      uintptr_t sample_uidx;
      if (!new_sample_idx_to_old) {
        sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
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
          write_iter = strcpya_k(write_iter, "NA");
        }
      }
      if (write_phenos) {
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          write_iter = AppendPhenoStr(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
          if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
            goto WritePsam_ret_WRITE_FAIL;
          }
        }
      } else {
        if (write_empty_pheno) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
        }
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto WritePsam_ret_WRITE_FAIL;
        }
      }
      AppendBinaryEoln(&write_iter);
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
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

/*
#ifdef __arm__
#  error "Unaligned accesses in BitvecResort()."
#endif
void BitvecResort(const uintptr_t* bitvec, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, unsigned char* writebuf) {
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
    for (uint32_t uii = 0; uii != cur_word_entry_ct; ++uii) {
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
  const uint32_t sample_ctl2_m1 = NypCtToWordCt(sample_ct) - 1;
  const uint32_t* new_sample_idx_to_old_iter = new_sample_idx_to_old;
  uintptr_t* writebuf_walias = S_CAST(uintptr_t*, writebuf);
  for (uint32_t widx = 0; widx != sample_ctl2_m1; ++widx) {
    uintptr_t cur_word = 0;
    // this is noticeably better than the ascending loop
    for (uint32_t uii = kBitsPerWordD2 - 1; ; --uii) {
      cur_word |= GetNyparrEntry(genovec, new_sample_idx_to_old_iter[uii]);
      if (!uii) {
        break;
      }
      cur_word = cur_word << 2;
    }
    writebuf_walias[widx] = cur_word;
    new_sample_idx_to_old_iter = &(new_sample_idx_to_old_iter[kBitsPerWordD2]);
  }
  const uint32_t cur_word_entry_ct = ModNz(sample_ct, kBitsPerWordD2);
  uintptr_t cur_word = 0;
  for (uint32_t uii = cur_word_entry_ct - 1; ; --uii) {
    cur_word |= GetNyparrEntry(genovec, new_sample_idx_to_old_iter[uii]);
    if (!uii) {
      break;
    }
    cur_word = cur_word << 2;
  }
  SubwordStore(cur_word, NypCtToByteCt(cur_word_entry_ct), &(writebuf_walias[sample_ctl2_m1]));
}

// Revised phaseraw:
//   4 byte het_ct, 4 byte explicit_phasepresent_ct
//   first half, up to (1 + (het_ct / kBitsPerWord)) words
//   second half, rounded up to vector boundary
void UnpackHphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, uint32_t raw_sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t het_ct = S_CAST(uint32_t, phaseraw[0]);
  const uintptr_t* aux2a = &(phaseraw[8 / kBytesPerWord]);
  if (!(aux2a[0] & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    ExpandBytearr(aux2a, all_hets, raw_sample_ctl, het_ct, 1, phaseinfo);
  } else {
    // bugfix (4 Mar 2018): need to pass raw_phasepresent_ct, not het_ct
#ifdef __LP64__
    const uint32_t raw_phasepresent_ct = phaseraw[0] >> 32;
#else
    const uint32_t raw_phasepresent_ct = phaseraw[1];
#endif
    const uintptr_t* aux2b = &(aux2a[1 + (het_ct / kBitsPerWord)]);
    ExpandBytearrNested(aux2b, aux2a, all_hets, raw_sample_ctl, raw_phasepresent_ct, 1, *phasepresent_ptr, phaseinfo);
  }
}

void UnpackHphaseSubset(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* __restrict sample_include, uint32_t sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  // const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  // const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
  const uint32_t het_ct = S_CAST(uint32_t, phaseraw[0]);
  const uintptr_t* aux2a = &(phaseraw[8 / kBytesPerWord]);
  if (!(aux2a[0] & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    ExpandThenSubsetBytearr(aux2a, all_hets, sample_include, het_ct, sample_ct, 1, phaseinfo);
  } else {
    const uint32_t first_half_word_ct = 1 + (het_ct / kBitsPerWord);
    // const uint32_t raw_phasepresent_ct = PopcountWords(phaseraw, first_half_word_ct) - 1;
#ifdef __LP64__
    const uint32_t raw_phasepresent_ct = phaseraw[0] >> 32;
#else
    const uint32_t raw_phasepresent_ct = phaseraw[1];
#endif
    const uintptr_t* aux2b = &(aux2a[first_half_word_ct]);

    // see "if (explicit_phasepresent) {}" block in PgrGetRaw().  Could
    // change this convention.
    ExpandThenSubsetBytearrNested(aux2b, aux2a, all_hets, sample_include, sample_ct, raw_phasepresent_ct, 1, *phasepresent_ptr, phaseinfo);
  }
}

void UnpackAndResortHphase(const uintptr_t* __restrict all_hets, const uintptr_t* __restrict phaseraw, const uintptr_t* sample_include, const uint32_t* old_sample_idx_to_new, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t** phasepresent_ptr, uintptr_t* __restrict phaseinfo) {
  const uintptr_t* aux2a_iter = &(phaseraw[8 / kBytesPerWord]);
  const uint32_t* old_sample_idx_to_new_iter = old_sample_idx_to_new;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uintptr_t aux2a_word = *aux2a_iter++;
  uint32_t read_idx_lowbits = 1;
  ZeroWArr(sample_ctl, phaseinfo);
  if (!(aux2a_word & 1)) {
    // phase always present
    *phasepresent_ptr = nullptr;
    for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
      uintptr_t new_phasepresent_word = all_hets[widx];
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + PopcountWord(new_phasepresent_word);
      uintptr_t tmp_phaseinfo_input_word = aux2a_word >> read_idx_lowbits;
      if (read_idx_lowbits_end >= kBitsPerWord) {
        // always safe to read an extra word off the end
        aux2a_word = *aux2a_iter++;
        if (read_idx_lowbits) {
          tmp_phaseinfo_input_word |= aux2a_word << (kBitsPerWord - read_idx_lowbits);
        }
      }
      // no need to mask off top bits of tmp_phaseinfo_input_word
      read_idx_lowbits = read_idx_lowbits_end % kBitsPerWord;
      if (!sample_include) {
#ifdef USE_AVX2
        uintptr_t phaseinfo_bits_to_set = _pdep_u64(tmp_phaseinfo_input_word, new_phasepresent_word);
        while (phaseinfo_bits_to_set) {
          const uint32_t sample_uidx_lowbits = ctzw(phaseinfo_bits_to_set);
          SetBit(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
          phaseinfo_bits_to_set &= phaseinfo_bits_to_set - 1;
        }
#else
        while (new_phasepresent_word) {
          const uint32_t sample_uidx_lowbits = ctzw(new_phasepresent_word);
          if (tmp_phaseinfo_input_word & 1) {
            SetBit(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
          }
          tmp_phaseinfo_input_word >>= 1;
          new_phasepresent_word &= new_phasepresent_word - 1;
        }
#endif
      } else {
#ifdef USE_AVX2
        uintptr_t phaseinfo_bits_to_set = _pdep_u64(tmp_phaseinfo_input_word, new_phasepresent_word) & sample_include[widx];
        while (phaseinfo_bits_to_set) {
          const uint32_t sample_uidx_lowbits = ctzw(phaseinfo_bits_to_set);
          SetBit(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
          phaseinfo_bits_to_set &= phaseinfo_bits_to_set - 1;
        }
#else
        uintptr_t masked_phasepresent_word = new_phasepresent_word & sample_include[widx];
        while (masked_phasepresent_word) {
          const uint32_t sample_uidx_lowbits = ctzw(masked_phasepresent_word);
          const uintptr_t lowmask = (k1LU << sample_uidx_lowbits) - k1LU;
          if ((tmp_phaseinfo_input_word >> PopcountWord(new_phasepresent_word & lowmask)) & 1) {
            SetBit(old_sample_idx_to_new_iter[sample_uidx_lowbits], phaseinfo);
          }
          masked_phasepresent_word &= masked_phasepresent_word - 1;
        }
#endif
      }
      old_sample_idx_to_new_iter = &(old_sample_idx_to_new_iter[kBitsPerWord]);
    }
    return;
  }
  uintptr_t* phasepresent = *phasepresent_ptr;
  const uint32_t het_ct = S_CAST(uint32_t, phaseraw[0]);
  const uintptr_t* phaseinfo_read_iter = &(phaseraw[(8 / kBytesPerWord) + 1 + (het_ct / kBitsPerWord)]);
  uintptr_t phaseinfo_read_word = *phaseinfo_read_iter++;
  uint32_t phaseinfo_read_idx_lowbits = 0;
  ZeroWArr(sample_ctl, phasepresent);
  for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
    uintptr_t geno_hets = all_hets[widx];
    if (geno_hets) {
      const uint32_t read_idx_lowbits_end = read_idx_lowbits + PopcountWord(geno_hets);
      uintptr_t tmp_phasepresent_input_word = aux2a_word >> read_idx_lowbits;
      if (read_idx_lowbits_end >= kBitsPerWord) {
        // always safe to read an extra word off the end, when
        // read_idx_lowbits_end == kBitsPerWord and we're at the last word
        aux2a_word = *aux2a_iter++;
        if (read_idx_lowbits) {
          tmp_phasepresent_input_word |= aux2a_word << (kBitsPerWord - read_idx_lowbits);
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
#ifdef USE_AVX2
          for (uintptr_t phasepresent_bits_to_set = _pdep_u64(tmp_phasepresent_input_word, geno_hets); ; ) {
            const uint32_t new_sample_idx = old_sample_idx_to_new_iter[ctzw(phasepresent_bits_to_set)];
            const uint32_t new_sample_widx = new_sample_idx / kBitsPerWord;
            const uint32_t new_sample_lowbits = new_sample_idx % kBitsPerWord;
            const uintptr_t shifted_bit = k1LU << new_sample_lowbits;
            phasepresent[new_sample_widx] |= shifted_bit;
            if (tmp_phaseinfo_input_word & 1) {
              phaseinfo[new_sample_widx] |= shifted_bit;
            }
            // branchless version doesn't seem to be any better here; probably
            // due to additional random memory access.
            // phaseinfo[new_sample_widx] |= (tmp_phaseinfo_input_word & 1) << new_sample_lowbits;

            phasepresent_bits_to_set &= phasepresent_bits_to_set - 1;
            if (!phasepresent_bits_to_set) {
              break;
            }
            tmp_phaseinfo_input_word >>= 1;
          }
#else
          for (; ; tmp_phasepresent_input_word >>= 1) {
            if (tmp_phasepresent_input_word & 1) {
              const uint32_t new_sample_idx = old_sample_idx_to_new_iter[ctzw(geno_hets)];
              const uint32_t new_sample_widx = new_sample_idx / kBitsPerWord;
              const uint32_t new_sample_lowbits = new_sample_idx % kBitsPerWord;
              const uintptr_t shifted_bit = k1LU << new_sample_lowbits;
              phasepresent[new_sample_widx] |= shifted_bit;
              if (tmp_phaseinfo_input_word & 1) {
                phaseinfo[new_sample_widx] |= shifted_bit;
              }
              if (tmp_phasepresent_input_word == 1) {
                break;
              }
              tmp_phaseinfo_input_word >>= 1;
            }
            geno_hets &= geno_hets - 1;
          }
#endif
        } else {
          const uintptr_t sample_include_word = sample_include[widx];
#ifdef USE_AVX2
          const uintptr_t phasepresent_word_expanded = _pdep_u64(tmp_phasepresent_input_word, geno_hets);
          uintptr_t phasepresent_bits_to_set = phasepresent_word_expanded & sample_include_word;
          if (phasepresent_bits_to_set) {
            // tmp_phaseinfo_input_word gives us the phasing state of the
            // positions in phasepresent_word_expanded.
            // However, we're only iterating over the positions in
            // (phasepresent_word_expanded & sample_include_word).
            // (can replace sample_include_word with phasepresent_bits_to_set
            // in this expression)
            uintptr_t collapsed_phaseinfo_input_word = _pext_u64(tmp_phaseinfo_input_word, _pext_u64(sample_include_word, phasepresent_word_expanded));
            while (1) {
              const uint32_t new_sample_idx = old_sample_idx_to_new_iter[ctzw(phasepresent_bits_to_set)];
              const uint32_t new_sample_widx = new_sample_idx / kBitsPerWord;
              const uint32_t new_sample_lowbits = new_sample_idx % kBitsPerWord;
              const uintptr_t shifted_bit = k1LU << new_sample_lowbits;
              phasepresent[new_sample_widx] |= shifted_bit;
              if (collapsed_phaseinfo_input_word & 1) {
                phaseinfo[new_sample_widx] |= shifted_bit;
              }

              phasepresent_bits_to_set &= phasepresent_bits_to_set - 1;
              if (!phasepresent_bits_to_set) {
                break;
              }
              collapsed_phaseinfo_input_word >>= 1;
            }
          }
#else
          for (; ; tmp_phasepresent_input_word >>= 1) {
            if (tmp_phasepresent_input_word & 1) {
              const uintptr_t geno_hets_lowbit = geno_hets & (-geno_hets);
              if (sample_include_word & geno_hets_lowbit) {
                const uint32_t sample_uidx_lowbits = ctzw(geno_hets_lowbit);
                const uint32_t new_sample_idx = old_sample_idx_to_new_iter[sample_uidx_lowbits];
                const uint32_t new_sample_widx = new_sample_idx / kBitsPerWord;
                const uint32_t new_sample_lowbits = new_sample_idx % kBitsPerWord;
                const uintptr_t shifted_bit = k1LU << new_sample_lowbits;
                phasepresent[new_sample_widx] |= shifted_bit;
                if (tmp_phaseinfo_input_word & 1) {
                  phaseinfo[new_sample_widx] |= shifted_bit;
                }
              }
              if (tmp_phasepresent_input_word == 1) {
                break;
              }
              tmp_phaseinfo_input_word >>= 1;
            }
            geno_hets &= geno_hets - 1;
          }
#endif
        }
      }
    }
    old_sample_idx_to_new_iter = &(old_sample_idx_to_new_iter[kBitsPerWord]);
  }
}


// these also work on dphaseraw
void CopyDosage(const uintptr_t* __restrict read_dosagepresent, const Dosage* read_dosagevals, uint32_t raw_sample_ct, uint32_t dosage_ct, uintptr_t* __restrict write_dosagepresent, Dosage* write_dosagevals, uint32_t* write_dosage_ct_ptr) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  *write_dosage_ct_ptr = dosage_ct;
  memcpy(write_dosagepresent, read_dosagepresent, raw_sample_ctl * sizeof(intptr_t));
  memcpy(write_dosagevals, read_dosagevals, dosage_ct * sizeof(Dosage));
}

uint32_t CopyAndResort8bit(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uint32_t* __restrict new_sample_idx_to_old, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals, uint32_t* __restrict cumulative_popcount_buf) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  FillCumulativePopcounts(src_subset, raw_sample_ctl, cumulative_popcount_buf);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl, dst_subset);
  const unsigned char* src_vals_uc = S_CAST(const unsigned char*, src_vals);
  unsigned char* dst_vals_uc = S_CAST(unsigned char*, dst_vals);
  unsigned char* dst_vals_iter = dst_vals_uc;
  // Tried word-based loop, was significantly worse
  for (uint32_t new_sample_idx = 0; new_sample_idx != sample_ct; ++new_sample_idx) {
    const uint32_t old_sample_idx = new_sample_idx_to_old[new_sample_idx];
    if (IsSet(src_subset, old_sample_idx)) {
      SetBit(new_sample_idx, dst_subset);
      const uint32_t old_dosagevals_idx = RawToSubsettedPos(src_subset, cumulative_popcount_buf, old_sample_idx);
      *dst_vals_iter++ = src_vals_uc[old_dosagevals_idx];
    }
  }
  return dst_vals_iter - dst_vals_uc;
}

uint32_t CopyAndResort16bit(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uint32_t* __restrict new_sample_idx_to_old, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals, uint32_t* __restrict cumulative_popcount_buf) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  FillCumulativePopcounts(src_subset, raw_sample_ctl, cumulative_popcount_buf);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl, dst_subset);
  const uint16_t* src_vals_u16 = S_CAST(const uint16_t*, src_vals);
  uint16_t* dst_vals_u16 = S_CAST(uint16_t*, dst_vals);
  uint16_t* dst_vals_iter = dst_vals_u16;
  // Tried word-based loop, was significantly worse
  for (uint32_t new_sample_idx = 0; new_sample_idx != sample_ct; ++new_sample_idx) {
    const uint32_t old_sample_idx = new_sample_idx_to_old[new_sample_idx];
    if (IsSet(src_subset, old_sample_idx)) {
      SetBit(new_sample_idx, dst_subset);
      const uint32_t old_dosagevals_idx = RawToSubsettedPos(src_subset, cumulative_popcount_buf, old_sample_idx);
      *dst_vals_iter++ = src_vals_u16[old_dosagevals_idx];
    }
  }
  return dst_vals_iter - dst_vals_u16;
}

// Requires trailing bits of genovec to be zeroed out.
// "Flat" = don't separate one_cts and two_cts.
void GetMFlatCounts64(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const PgenVariant* pgvp, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t allele_ct, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* all_dosages) {
  if (sample_ct == raw_sample_ct) {
    GenoarrCountFreqsUnsafe(pgvp->genovec, sample_ct, genocounts);
  } else {
    GenoarrCountSubsetFreqs(pgvp->genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  all_dosages[0] = 2 * genocounts[0] + genocounts[1];
  all_dosages[1] = 2 * genocounts[2] + genocounts[1];
  ZeroU64Arr(allele_ct - 2, &(all_dosages[2]));
  const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
  const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (sample_ct == raw_sample_ct) {
    all_dosages[1] -= patch_01_ct + 2 * patch_10_ct;
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      all_dosages[patch_01_vals[uii]] += 1;
    }
    const uint32_t patch_10_ct_x2 = patch_10_ct * 2;
    for (uint32_t uii = 0; uii != patch_10_ct_x2; ++uii) {
      all_dosages[patch_10_vals[uii]] += 1;
    }
  } else {
    if (patch_01_ct) {
      const uintptr_t* patch_01_set = pgvp->patch_01_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_01_bits = patch_01_set[0];
      uint32_t subsetted_patch_01_ct = 0;
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_01_set, &sample_widx, &patch_01_bits);
        if (sample_include[sample_widx] & lowbit) {
          all_dosages[patch_01_vals[uii]] += 1;
          ++subsetted_patch_01_ct;
        }
      }
      all_dosages[1] -= subsetted_patch_01_ct;
    }
    if (patch_10_ct) {
      const uintptr_t* patch_10_set = pgvp->patch_10_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_10_bits = patch_10_set[0];
      uint32_t subsetted_patch_10_ct = 0;
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
        if (sample_include[sample_widx] & lowbit) {
          all_dosages[patch_10_vals[2 * uii]] += 1;
          all_dosages[patch_10_vals[2 * uii + 1]] += 1;
          ++subsetted_patch_10_ct;
        }
      }
      all_dosages[1] -= 2 * subsetted_patch_10_ct;
    }
  }
}

void GetMCounts64(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const PgenVariant* pgvp, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t allele_ct, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict one_cts, uint64_t* __restrict two_cts) {
  // This mirrors GetMultiallelicCountsAndDosage16s().
  if (sample_ct == raw_sample_ct) {
    GenoarrCountFreqsUnsafe(pgvp->genovec, sample_ct, genocounts);
  } else {
    GenoarrCountSubsetFreqs(pgvp->genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  one_cts[0] = genocounts[1];
  one_cts[1] = genocounts[1];
  ZeroU64Arr(allele_ct - 2, &(one_cts[2]));
  two_cts[0] = genocounts[0];
  two_cts[1] = genocounts[2];
  ZeroU64Arr(allele_ct - 2, &(two_cts[2]));
  const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
  const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (sample_ct == raw_sample_ct) {
    one_cts[1] -= patch_01_ct;
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      one_cts[patch_01_vals[uii]] += 1;
    }
    two_cts[1] -= patch_10_ct;
    const AlleleCode* patch_10_vals_iter = patch_10_vals;
    for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
      const AlleleCode code_lo = *patch_10_vals_iter++;
      const AlleleCode code_hi = *patch_10_vals_iter++;
      if (code_lo == code_hi) {
        two_cts[code_lo] += 1;
      } else {
        one_cts[code_lo] += 1;
        one_cts[code_hi] += 1;
      }
    }
  } else {
    if (patch_01_ct) {
      const uintptr_t* patch_01_set = pgvp->patch_01_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_01_bits = patch_01_set[0];
      uint32_t subsetted_patch_01_ct = 0;
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_01_set, &sample_widx, &patch_01_bits);
        if (sample_include[sample_widx] & lowbit) {
          one_cts[patch_01_vals[uii]] += 1;
          ++subsetted_patch_01_ct;
        }
      }
      one_cts[1] -= subsetted_patch_01_ct;
    }
    if (patch_10_ct) {
      const uintptr_t* patch_10_set = pgvp->patch_10_set;
      uintptr_t sample_widx = 0;
      uintptr_t patch_10_bits = patch_10_set[0];
      uint32_t subsetted_patch_10_ct = 0;
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
        if (sample_include[sample_widx] & lowbit) {
          ++subsetted_patch_10_ct;
          const AlleleCode code_lo = patch_10_vals[2 * uii];
          const AlleleCode code_hi = patch_10_vals[2 * uii + 1];
          if (code_lo == code_hi) {
            two_cts[code_lo] += 1;
          } else {
            one_cts[code_lo] += 1;
            one_cts[code_hi] += 1;
          }
        }
      }
      two_cts[1] -= subsetted_patch_10_ct;
    }
  }
}

typedef struct LoadAlleleAndGenoCountsCtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  uintptr_t* sample_include_interleaved_vec;
  uint32_t* sample_include_cumulative_popcounts;
  const uintptr_t* sex_male;
  uintptr_t* sex_male_interleaved_vec;
  uint32_t* sex_male_cumulative_popcounts;
  uintptr_t* nosex_interleaved_vec;
  const uintptr_t* founder_info;
  uintptr_t* founder_info_interleaved_vec;
  uint32_t* founder_info_cumulative_popcounts;
  uintptr_t* founder_male;
  uintptr_t* founder_male_interleaved_vec;
  uint32_t* founder_male_cumulative_popcounts;
  uintptr_t* founder_nosex_interleaved_vec;
  uint32_t raw_sample_ct;
  uint32_t sample_ct;
  uint32_t founder_ct;
  uint32_t male_ct;
  uint32_t nosex_ct;
  uint32_t founder_male_ct;
  uint32_t founder_nosex_ct;
  uint32_t first_hap_uidx;
  uint32_t is_minimac3_r2;

  PgenReader** pgr_ptrs;

  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uint64_t** all_dosages;
  uint32_t* read_variant_uidx_starts;

  // shouldn't need array, or errno storage, since kPglRetMalformedInput is the
  // only possible error for now
  PglErr reterr;

  uint32_t cur_block_size;

  unsigned char* allele_presents_bytearr;
  uint64_t* allele_ddosages;
  STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts);
  uint32_t* variant_missing_hc_cts;
  uint32_t* variant_missing_dosage_cts;
  uint32_t* variant_hethap_cts;
  uint64_t* founder_allele_ddosages;
  STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts);
  STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts);
  double* imp_r2_vals;
} LoadAlleleAndGenoCountsCtx;

THREAD_FUNC_DECL LoadAlleleAndGenoCountsThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  LoadAlleleAndGenoCountsCtx* ctx = S_CAST(LoadAlleleAndGenoCountsCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = ctx->cip;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t subset_ct = (ctx->founder_info != nullptr) + 1;
  const uint32_t raw_sample_ct = ctx->raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t first_hap_uidx = ctx->first_hap_uidx;
  const uint32_t is_minimac3_r2 = ctx->is_minimac3_r2;
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = ctx->genovecs[tidx];
  SetPgvThreadMhcNull(raw_sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (ctx->dosage_presents) {
    pgv.dosage_present = ctx->dosage_presents[tidx];
    pgv.dosage_main = ctx->dosage_mains[tidx];
  }
  uint64_t* all_dosages = nullptr;
  if (ctx->all_dosages) {
    all_dosages = ctx->all_dosages[tidx];
  }
  uint32_t is_y = 0;
  uint32_t is_nonxy_haploid = 0;
  uint32_t x_start = 0;
  uint32_t x_code;
  if (XymtExists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
    x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
  }
  uint32_t allele_ct = 2;
  do {
    const uintptr_t cur_block_size = ctx->cur_block_size;
    // no overflow danger since cur_block_size <= 2^16, tidx < (2^16 - 1)
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_size) / thread_ct;
    const uintptr_t* sample_include = ctx->sample_include;
    const uintptr_t* sample_include_interleaved_vec = ctx->sample_include_interleaved_vec;
    const uint32_t* sample_include_cumulative_popcounts = ctx->sample_include_cumulative_popcounts;
    const uintptr_t* sex_male = ctx->sex_male;
    const uintptr_t* sex_male_interleaved_vec = ctx->sex_male_interleaved_vec;
    const uint32_t* sex_male_cumulative_popcounts = ctx->sex_male_cumulative_popcounts;
    const uintptr_t* nosex_interleaved_vec = ctx->nosex_interleaved_vec;
    uint32_t sample_ct = ctx->sample_ct;
    uint32_t male_ct = ctx->male_ct;
    uint32_t nosex_ct = ctx->nosex_ct;
    unsigned char* allele_presents_bytearr = ctx->allele_presents_bytearr;
    uint64_t* allele_ddosages = ctx->allele_ddosages;
    STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts) = ctx->raw_geno_cts;
    uint32_t* variant_missing_hc_cts = ctx->variant_missing_hc_cts;
    uint32_t* variant_missing_dosage_cts = ctx->variant_missing_dosage_cts;
    uint32_t* variant_hethap_cts = ctx->variant_hethap_cts;
    STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts) = ctx->x_male_geno_cts;
    STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts) = ctx->x_nosex_geno_cts;
    double* imp_r2_vals = ctx->imp_r2_vals;
    pgv.dosage_ct = 0;
    for (uint32_t subset_idx = 0; ; ) {
      // bugfix (29 Dec 2019): this boolean can change with subset_idx
      const uint32_t no_multiallelic_branch = (!variant_hethap_cts) && (!allele_presents_bytearr) && (!allele_ddosages) && (!imp_r2_vals);
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, pgrp, &pssi);
      uint32_t cur_idx = (tidx * cur_block_size) / thread_ct;
      uintptr_t variant_uidx_base;
      uintptr_t variant_include_bits;
      BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
      uint32_t chr_end = 0;
      uint32_t is_x_or_y = 0;
      PglErr reterr = kPglRetSuccess;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      STD_ARRAY_DECL(uint32_t, 4, sex_specific_genocounts);
      for (; cur_idx != cur_idx_end; ++cur_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (variant_uidx >= chr_end) {
          const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          is_y = 0;
          is_nonxy_haploid = 0;
          if (chr_idx == x_code) {
            is_x_or_y = 1;
            PgrClearSampleSubsetIndex(pgrp, &pssi);
          } else if (chr_idx == y_code) {
            is_x_or_y = 1;
            is_y = 1;
            // ugh
            if ((!allele_presents_bytearr) || (sample_ct == male_ct)) {
              PgrSetSampleSubsetIndex(sex_male_cumulative_popcounts, pgrp, &pssi);
            } else {
              PgrClearSampleSubsetIndex(pgrp, &pssi);
            }
          } else {
            if (is_x_or_y) {
              PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, pgrp, &pssi);
            }
            is_x_or_y = 0;
            // true for MT
            is_nonxy_haploid = IsSet(cip->haploid_mask, chr_idx);
          }
        }
        uintptr_t cur_allele_idx_offset;
        if (!allele_idx_offsets) {
          cur_allele_idx_offset = 2 * variant_uidx;
        } else {
          cur_allele_idx_offset = allele_idx_offsets[variant_uidx];
          allele_ct = allele_idx_offsets[variant_uidx + 1] - cur_allele_idx_offset;
        }
        uint32_t hethap_ct;
        if ((allele_ct == 2) || no_multiallelic_branch) {
          uint64_t cur_dosages[2];
          if (!is_x_or_y) {
            reterr = PgrGetDCounts(sample_include, sample_include_interleaved_vec, pssi, sample_ct, variant_uidx, is_minimac3_r2, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, genocounts, cur_dosages);
            if (unlikely(reterr)) {
              ctx->reterr = reterr;
              break;
            }
            if (allele_presents_bytearr) {
              if (cur_dosages[0]) {
                allele_presents_bytearr[cur_allele_idx_offset] = 128;
              }
              if (cur_dosages[1]) {
                allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
              }
            }
            if (!is_nonxy_haploid) {
              hethap_ct = 0;
              if (allele_ddosages) {
                // ...but save all allele counts here.
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0] * 2;
                allele_ddosages[cur_allele_idx_offset + 1] = cur_dosages[1] * 2;
              }
            } else {
              // this hethap_ct can be inaccurate in multiallelic case
              hethap_ct = genocounts[1];
              if (imp_r2_vals && (!is_minimac3_r2)) {
                // Assuming the input data isn't malformed "phased haploid",
                // minimac3-r2 is independent of haploid/diploid state; only
                // mach-r2 requires a haploid correction.
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_ddosages) {
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0];
                allele_ddosages[cur_allele_idx_offset + 1] = cur_dosages[1];
              }
            }
          } else if (is_y) {
            if ((!allele_presents_bytearr) || (sample_ct == male_ct)) {
              reterr = PgrGetDCounts(sex_male, sex_male_interleaved_vec, pssi, male_ct, variant_uidx, 0, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, genocounts, cur_dosages);
              if (unlikely(reterr)) {
                ctx->reterr = reterr;
                break;
              }
              hethap_ct = genocounts[1];
              if (imp_r2_vals && (!is_minimac3_r2)) {
                // note that female/unknown-sex are not counted here
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_presents_bytearr) {
                if (cur_dosages[0]) {
                  allele_presents_bytearr[cur_allele_idx_offset] = 128;
                }
                if (cur_dosages[1]) {
                  allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                }
              }
              if (allele_ddosages) {
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0];
                allele_ddosages[cur_allele_idx_offset + 1] = cur_dosages[1];
              }
            } else {
              // ugh, need to count female/unknown-sex for allele_presents and
              // ignore elsewhere
              reterr = PgrGetD(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
              if (unlikely(reterr)) {
                ctx->reterr = reterr;
                break;
              }
              const uint32_t dosage_is_relevant = pgv.dosage_ct && ((sample_ct == raw_sample_ct) || (!IntersectionIsEmpty(sample_include, pgv.dosage_present, raw_sample_ctl)));
              if (dosage_is_relevant) {
                // at least one dosage value is present, that's all we need to
                // know
                allele_presents_bytearr[cur_allele_idx_offset] = 128;
                allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
              } else {
                // only hardcalls matter
                // bugfix (31 Jul 2018): forgot to initialize genocounts here
                // possible todo: use a specialized function which just checks
                // which alleles exist
                if (sample_ct == raw_sample_ct) {
                  ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
                  GenoarrCountFreqsUnsafe(pgv.genovec, sample_ct, genocounts);
                } else {
                  GenoarrCountSubsetFreqs(pgv.genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
                }
                if (genocounts[0] || genocounts[1]) {
                  allele_presents_bytearr[cur_allele_idx_offset] = 128;
                }
                if (genocounts[1] || genocounts[2]) {
                  allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                }
              }
              GenoarrCountSubsetFreqs(pgv.genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, genocounts);
              hethap_ct = genocounts[1];
              // x2, x4 since this is haploid
              uintptr_t alt1_ct_x2 = genocounts[2] * 2 + hethap_ct;
              uintptr_t alt1_sq_sum_x4 = genocounts[2] * (4 * k1LU) + hethap_ct;
              uint64_t alt1_ddosage = 0;  // in 32768ths
              uint64_t alt1_ddosage_sq_sum = 0;
              uint32_t additional_dosage_ct = 0;
              if (dosage_is_relevant) {
                uintptr_t sample_widx = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                uint32_t sample_uidx = 0;
                for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                  const uintptr_t lowbit = BitIter1y(pgv.dosage_present, &sample_widx, &dosage_present_bits);
                  if (sample_include[sample_widx] & lowbit) {
                    const uintptr_t cur_dosage_val = pgv.dosage_main[dosage_idx];
                    alt1_ddosage += cur_dosage_val;
                    alt1_ddosage_sq_sum += cur_dosage_val * cur_dosage_val;
                    const uintptr_t hardcall_code = GetNyparrEntry(pgv.genovec, sample_uidx);
                    if (hardcall_code != 3) {
                      alt1_ct_x2 -= hardcall_code;
                      alt1_sq_sum_x4 -= hardcall_code * hardcall_code;
                    } else {
                      ++additional_dosage_ct;
                    }
                  }
                }
              }
              const uintptr_t obs_ct = male_ct + additional_dosage_ct - genocounts[3];
              alt1_ddosage += alt1_ct_x2 * S_CAST(uint64_t, kDosageMid);
              alt1_ddosage_sq_sum += alt1_sq_sum_x4 * 0x10000000LLU;
              cur_dosages[0] = obs_ct * S_CAST(uint64_t, kDosageMax) - alt1_ddosage;
              cur_dosages[1] = alt1_ddosage;
              if (imp_r2_vals) {
                // minimac3-r2 and mach-r2 are identical in haploid case
                const double dosage_sumd = u63tod(alt1_ddosage);
                const double dosage_avg = dosage_sumd / u31tod(obs_ct);
                const double dosage_variance = u63tod(alt1_ddosage_sq_sum) - dosage_sumd * dosage_avg;
                imp_r2_vals[variant_uidx] = dosage_variance / (dosage_sumd * (32768 - dosage_avg));
              }
              if (allele_ddosages) {
                allele_ddosages[cur_allele_idx_offset] = cur_dosages[0];
                allele_ddosages[cur_allele_idx_offset + 1] = alt1_ddosage;
              }
            }
          } else {
            // chrX
            reterr = PgrGetD(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
            if (unlikely(reterr)) {
              ctx->reterr = reterr;
              break;
            }
            if (sample_ct == raw_sample_ct) {
              ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
              GenoarrCountFreqsUnsafe(pgv.genovec, sample_ct, genocounts);
            } else {
              GenoarrCountSubsetFreqs(pgv.genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
            }
            GenoarrCountSubsetFreqs(pgv.genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, sex_specific_genocounts);
            hethap_ct = sex_specific_genocounts[1];
            // Could compute imputation r2 iff there are no unknown-sex
            // samples, but probably not worth it since larger datasets could
            // have a small number of Klinefelter syndrome cases, etc. coded as
            // unknown-sex, and we don't want to discourage their inclusion;
            // let's delegate that chrX filter to other software for now.

            if (allele_presents_bytearr) {
              if (pgv.dosage_ct && ((sample_ct == raw_sample_ct) || (!IntersectionIsEmpty(sample_include, pgv.dosage_present, raw_sample_ctl)))) {
                // at least one dosage value is present, that's all we need to
                // know
                allele_presents_bytearr[cur_allele_idx_offset] = 128;
                allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
              } else {
                // only hardcalls matter
                if (genocounts[0] || genocounts[1]) {
                  allele_presents_bytearr[cur_allele_idx_offset] = 128;
                }
                if (genocounts[1] || genocounts[2]) {
                  allele_presents_bytearr[cur_allele_idx_offset + 1] = 128;
                }
              }
            }
            if (allele_ddosages) {
              uintptr_t alt1_ct = 4 * genocounts[2] + 2 * genocounts[1] - 2 * sex_specific_genocounts[2] - hethap_ct;  // nonmales count twice
              uint64_t alt1_ddosage = 0;  // in 32768ths, nonmales count twice
              uint32_t additional_dosage_ct = 0;  // missing hardcalls only; nonmales count twice
              // bugfix (12 Jul 2018): dosage_present may be null if dosage_ct
              // == 0
              if (pgv.dosage_ct) {
                uintptr_t sample_uidx_base = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                if (sample_ct == raw_sample_ct) {
                  for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                    const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
                    const uintptr_t cur_dosage_val = pgv.dosage_main[dosage_idx];
                    const uintptr_t sex_multiplier = 2 - IsSet(sex_male, sample_uidx);
                    alt1_ddosage += cur_dosage_val * sex_multiplier;

                    // could call GenoarrCountSubsetIntersectFreqs() twice
                    // instead, but since we've already manually extracted the
                    // sex bit it probably doesn't help?
                    const uintptr_t hardcall_code = GetNyparrEntry(pgv.genovec, sample_uidx);
                    if (hardcall_code != 3) {
                      alt1_ct -= hardcall_code * sex_multiplier;
                    } else {
                      additional_dosage_ct += sex_multiplier;
                    }
                  }
                } else {
                  for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                    const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
                    if (IsSet(sample_include, sample_uidx)) {
                      const uintptr_t cur_dosage_val = pgv.dosage_main[dosage_idx];
                      const uintptr_t sex_multiplier = 2 - IsSet(sex_male, sample_uidx);
                      alt1_ddosage += cur_dosage_val * sex_multiplier;
                      const uintptr_t hardcall_code = GetNyparrEntry(pgv.genovec, sample_uidx);
                      if (hardcall_code != 3) {
                        alt1_ct -= hardcall_code * sex_multiplier;
                      } else {
                        additional_dosage_ct += sex_multiplier;
                      }
                    }
                  }
                }
              }
              alt1_ddosage += alt1_ct * S_CAST(uint64_t, kDosageMid);

              // bugfix (14 May 2018): this didn't correctly distinguish
              // between missing vs. 'replaced' hardcalls
              const uintptr_t weighted_obs_ct = (2 * (sample_ct - genocounts[3]) - male_ct + sex_specific_genocounts[3] + additional_dosage_ct) * (2 * k1LU);

              allele_ddosages[cur_allele_idx_offset] = weighted_obs_ct * S_CAST(uint64_t, kDosageMid) - alt1_ddosage;
              allele_ddosages[cur_allele_idx_offset + 1] = alt1_ddosage;
            }
            if (x_male_geno_cts) {
              STD_ARRAY_REF(uint32_t, 3) cur_x_male_geno_cts = x_male_geno_cts[variant_uidx - x_start];
              cur_x_male_geno_cts[0] = sex_specific_genocounts[0];
              cur_x_male_geno_cts[1] = sex_specific_genocounts[1];
              cur_x_male_geno_cts[2] = sex_specific_genocounts[2];
              if (x_nosex_geno_cts) {
                GenoarrCountSubsetFreqs(pgv.genovec, nosex_interleaved_vec, raw_sample_ct, nosex_ct, sex_specific_genocounts);
                STD_ARRAY_REF(uint32_t, 3) cur_nosex_geno_cts = x_nosex_geno_cts[variant_uidx - x_start];
                cur_nosex_geno_cts[0] = sex_specific_genocounts[0];
                cur_nosex_geno_cts[1] = sex_specific_genocounts[1];
                cur_nosex_geno_cts[2] = sex_specific_genocounts[2];
              }
            }
          }
          if (variant_missing_dosage_cts) {
            uint32_t missing_dosage_ct;
            if (!is_x_or_y) {
              missing_dosage_ct = sample_ct - ((cur_dosages[0] + cur_dosages[1]) / kDosageMax);
            } else if (is_y) {
              missing_dosage_ct = male_ct - ((cur_dosages[0] + cur_dosages[1]) / kDosageMax);
            } else {
              if (pgv.dosage_ct) {
                ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
                missing_dosage_ct = GenoarrCountMissingInvsubsetUnsafe(pgv.genovec, pgv.dosage_present, raw_sample_ct);
              } else {
                missing_dosage_ct = genocounts[3];
              }
            }
            variant_missing_dosage_cts[variant_uidx] = missing_dosage_ct;
          }
        } else {
          // multiallelic cases
          if (!is_x_or_y) {
            reterr = PgrGetMDCounts(sample_include, sample_include_interleaved_vec, pssi, sample_ct, variant_uidx, is_minimac3_r2, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, &hethap_ct, genocounts, all_dosages);
            if (unlikely(reterr)) {
              ctx->reterr = reterr;
              break;
            }
            if (allele_presents_bytearr) {
              for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                if (all_dosages[aidx]) {
                  allele_presents_bytearr[cur_allele_idx_offset + aidx] = 128;
                }
              }
            }
            if (!is_nonxy_haploid) {
              hethap_ct = 0;
              if (allele_ddosages) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  allele_ddosages[cur_allele_idx_offset + aidx] = all_dosages[aidx] * 2;
                }
              }
            } else {
              if (imp_r2_vals && (!is_minimac3_r2)) {
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_ddosages) {
                memcpy(&(allele_ddosages[cur_allele_idx_offset]), all_dosages, allele_ct * sizeof(int64_t));
              }
            }
          } else if (is_y) {
            if ((!allele_presents_bytearr) || (sample_ct == male_ct)) {
              reterr = PgrGetMDCounts(sex_male, sex_male_interleaved_vec, pssi, male_ct, variant_uidx, 0, pgrp, imp_r2_vals? (&(imp_r2_vals[variant_uidx])) : nullptr, &hethap_ct, genocounts, all_dosages);
              if (unlikely(reterr)) {
                ctx->reterr = reterr;
                break;
              }
              if (imp_r2_vals && (!is_minimac3_r2)) {
                imp_r2_vals[variant_uidx] *= 0.5;
              }
              if (allele_presents_bytearr) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  if (all_dosages[aidx]) {
                    allele_presents_bytearr[cur_allele_idx_offset + aidx] = 128;
                  }
                }
              }
              if (allele_ddosages) {
                memcpy(&(allele_ddosages[cur_allele_idx_offset]), all_dosages, allele_ct * sizeof(int64_t));
              }
            } else {
              // need to count female/unknown-sex for allele_presents and
              // ignore elsewhere
              reterr = PgrGetM(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, &pgv);
              if (unlikely(reterr)) {
                ctx->reterr = reterr;
                break;
              }
              // possible todo: use a specialized function which just checks
              // which alleles exist
              ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
              GetMFlatCounts64(sample_include, sample_include_interleaved_vec, &pgv, raw_sample_ct, sample_ct, allele_ct, genocounts, all_dosages);
              for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                if (all_dosages[aidx]) {
                  allele_presents_bytearr[cur_allele_idx_offset + aidx] = 128;
                }
              }

              uint64_t* two_cts = &(all_dosages[allele_ct]);
              GetMCounts64(sex_male, sex_male_interleaved_vec, &pgv, raw_sample_ct, male_ct, allele_ct, genocounts, all_dosages, two_cts);
              uintptr_t hethap_x2 = 0;
              for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
                hethap_x2 += all_dosages[aidx];
              }
              hethap_ct = hethap_x2 / 2;
              if (allele_ddosages) {
                for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                  allele_ddosages[cur_allele_idx_offset + aidx] = all_dosages[aidx] * kDosageMid + two_cts[aidx] * kDosageMax;
                }
              }
              if (imp_r2_vals) {
                for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
                  const uint64_t one_ct = allele_ddosages[aidx];
                  const uint64_t two_ct = two_cts[aidx];
                  // now sums
                  allele_ddosages[aidx] = one_ct * kDosageMid + two_ct * kDosageMax;
                  // now ssqs
                  two_cts[aidx] = one_ct * kDosageMid * kDosageMid + two_ct * kDosageMax * kDosageMax;
                }
                imp_r2_vals[variant_uidx] = 0.5 * MultiallelicDiploidMachR2(all_dosages, two_cts, male_ct - genocounts[3], allele_ct);
              }
            }
          } else {
            // chrX
            // multiallelic dosages not supported yet
            reterr = PgrGetM(nullptr, pssi, raw_sample_ct, variant_uidx, pgrp, &pgv);
            if (unlikely(reterr)) {
              ctx->reterr = reterr;
              break;
            }
            ZeroTrailingNyps(raw_sample_ct, pgv.genovec);
            // We don't attempt to compute imp_r2 on chrX, so flat counts are
            // fine.
            GetMFlatCounts64(sample_include, sample_include_interleaved_vec, &pgv, raw_sample_ct, sample_ct, allele_ct, genocounts, all_dosages);

            // Double all counts, then subtract male counts.
            for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
              all_dosages[aidx] *= 2;
            }
            GenoarrCountSubsetFreqs(pgv.genovec, sex_male_interleaved_vec, raw_sample_ct, male_ct, sex_specific_genocounts);
            hethap_ct = sex_specific_genocounts[1];
            if (male_ct) {
              all_dosages[0] -= 2 * sex_specific_genocounts[0] + hethap_ct;

              // may underflow
              all_dosages[1] -= 2 * sex_specific_genocounts[2] + hethap_ct;

              if (pgv.patch_01_ct) {
                uintptr_t sample_widx = 0;
                uintptr_t patch_01_bits = pgv.patch_01_set[0];
                uint32_t male_patch_01_ct = 0;
                for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                  const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &patch_01_bits);
                  if (sex_male[sample_widx] & lowbit) {
                    ++male_patch_01_ct;
                    all_dosages[pgv.patch_01_vals[uii]] -= 1;
                  }
                }
                all_dosages[1] += male_patch_01_ct;
              }
              if (pgv.patch_10_ct) {
                uintptr_t sample_widx = 0;
                uintptr_t patch_10_bits = pgv.patch_10_set[0];
                uint32_t male_patch_10_ct = 0;
                for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                  const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &patch_10_bits);
                  if (sex_male[sample_widx] & lowbit) {
                    ++male_patch_10_ct;
                    const AlleleCode code_lo = pgv.patch_10_vals[2 * uii];
                    const AlleleCode code_hi = pgv.patch_10_vals[2 * uii + 1];
                    all_dosages[code_lo] -= 1;
                    all_dosages[code_hi] -= 1;
                    hethap_ct += (code_lo != code_hi);
                  }
                }
                all_dosages[1] += male_patch_10_ct * 2;
              }
            }
            if (allele_presents_bytearr) {
              for (uintptr_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if (all_dosages[allele_idx]) {
                  allele_presents_bytearr[cur_allele_idx_offset + allele_idx] = 128;
                }
              }
            }
            if (allele_ddosages) {
              for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
                allele_ddosages[cur_allele_idx_offset + aidx] = all_dosages[aidx] * kDosageMid;
              }
            }
            if (x_male_geno_cts) {
              STD_ARRAY_REF(uint32_t, 3) cur_x_male_geno_cts = x_male_geno_cts[variant_uidx - x_start];
              cur_x_male_geno_cts[0] = sex_specific_genocounts[0];
              cur_x_male_geno_cts[1] = sex_specific_genocounts[1];
              cur_x_male_geno_cts[2] = sex_specific_genocounts[2];
              if (x_nosex_geno_cts) {
                GenoarrCountSubsetFreqs(pgv.genovec, nosex_interleaved_vec, raw_sample_ct, nosex_ct, sex_specific_genocounts);
                STD_ARRAY_REF(uint32_t, 3) cur_nosex_geno_cts = x_nosex_geno_cts[variant_uidx - x_start];
                cur_nosex_geno_cts[0] = sex_specific_genocounts[0];
                cur_nosex_geno_cts[1] = sex_specific_genocounts[1];
                cur_nosex_geno_cts[2] = sex_specific_genocounts[2];
              }
            }
          }
          if (variant_missing_dosage_cts) {
            // multiallelic dosage not supported yet
            variant_missing_dosage_cts[variant_uidx] = genocounts[3];
          }
        }
        if (raw_geno_cts) {
          STD_ARRAY_REF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
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
      }
      if ((++subset_idx == subset_ct) || reterr) {
        break;
      }
      sample_include = ctx->founder_info;
      sample_include_interleaved_vec = ctx->founder_info_interleaved_vec;
      sample_include_cumulative_popcounts = ctx->founder_info_cumulative_popcounts;
      sex_male = ctx->founder_male;
      sex_male_interleaved_vec = ctx->founder_male_interleaved_vec;
      sex_male_cumulative_popcounts = ctx->founder_male_cumulative_popcounts;

      nosex_interleaved_vec = ctx->founder_nosex_interleaved_vec;

      sample_ct = ctx->founder_ct;
      male_ct = ctx->founder_male_ct;
      nosex_ct = ctx->founder_nosex_ct;
      allele_presents_bytearr = nullptr;
      allele_ddosages = ctx->founder_allele_ddosages;
      variant_missing_hc_cts = nullptr;
      variant_missing_dosage_cts = nullptr;
      raw_geno_cts = ctx->founder_raw_geno_cts;
      x_male_geno_cts = ctx->founder_x_male_geno_cts;
      x_nosex_geno_cts = ctx->founder_x_nosex_geno_cts;
      imp_r2_vals = nullptr;
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr LoadAlleleAndGenoCounts(const uintptr_t* sample_include, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t first_hap_uidx, uint32_t is_minimac3_r2, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* allele_presents, uint64_t* allele_ddosages, uint64_t* founder_allele_ddosages, uint32_t* variant_missing_hc_cts, uint32_t* variant_missing_dosage_cts, uint32_t* variant_hethap_cts, STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, x_nosex_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), double* imp_r2_vals) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  LoadAlleleAndGenoCountsCtx ctx;
  {
    if (!variant_ct) {
      goto LoadAlleleAndGenoCounts_ret_1;
    }

    // four cases:
    // 1. allele_ddosages, raw_geno_cts, and/or variant_missing_{hc,dosage}_cts
    //    required, and that's it
    // 2. founder_allele_ddosages and/or founder_raw_geno_cts required, and
    //    that's it
    // 3. both required, and founder_ct != sample_ct.
    // 4. both required, and founder_ct == sample_ct.  caller is expected to
    //    make founder_allele_ddosages and allele_ddosages point to the same
    //    memory, ditto for founder_raw_geno_cts/raw_geno_cts.
    const uint32_t only_founder_cts_required = (!allele_presents) && (!allele_ddosages) && (!raw_geno_cts) && (!variant_missing_hc_cts) && (!variant_missing_dosage_cts);
    const uint32_t two_subsets_required = (founder_ct != sample_ct) && (!only_founder_cts_required) && (founder_allele_ddosages || founder_raw_geno_cts);
    ctx.cip = cip;
    ctx.sample_include = only_founder_cts_required? founder_info : sample_include;
    ctx.raw_sample_ct = raw_sample_ct;
    ctx.sample_ct = only_founder_cts_required? founder_ct : sample_ct;
    ctx.male_ct = male_ct;
    ctx.allele_ddosages = only_founder_cts_required? founder_allele_ddosages : allele_ddosages;
    ctx.raw_geno_cts = only_founder_cts_required? founder_raw_geno_cts : raw_geno_cts;
    ctx.x_male_geno_cts = only_founder_cts_required? founder_x_male_geno_cts : x_male_geno_cts;
    ctx.x_nosex_geno_cts = only_founder_cts_required? founder_x_nosex_geno_cts : x_nosex_geno_cts;
    ctx.imp_r2_vals = imp_r2_vals;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t raw_sample_ctv = BitCtToVecCt(raw_sample_ct);
    if (unlikely(
            bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.sample_include_interleaved_vec) ||
            bigstack_alloc_u32(raw_sample_ctl, &ctx.sample_include_cumulative_popcounts) ||
            bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.sex_male_interleaved_vec) ||
            bigstack_alloc_u32(raw_sample_ctl, &ctx.sex_male_cumulative_popcounts))) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    FillInterleavedMaskVec(ctx.sample_include, raw_sample_ctv, ctx.sample_include_interleaved_vec);
    FillCumulativePopcounts(ctx.sample_include, raw_sample_ctl, ctx.sample_include_cumulative_popcounts);
    if ((founder_ct == sample_ct) || (!only_founder_cts_required)) {
      ctx.sex_male = sex_male;
    } else {
      // no nonfounder counts required
      uintptr_t* new_sex_male;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &new_sex_male))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      BitvecAndCopy(sex_male, founder_info, raw_sample_ctl, new_sex_male);
      ZeroTrailingWords(raw_sample_ctl, new_sex_male);
      ctx.sex_male = new_sex_male;
    }
    FillInterleavedMaskVec(ctx.sex_male, raw_sample_ctv, ctx.sex_male_interleaved_vec);
    FillCumulativePopcounts(ctx.sex_male, raw_sample_ctl, ctx.sex_male_cumulative_popcounts);
    if (!(x_nosex_geno_cts || founder_x_nosex_geno_cts)) {
      nosex_ct = 0;
    }
    ctx.nosex_ct = nosex_ct;
    ctx.nosex_interleaved_vec = nullptr;
    uintptr_t* nosex_buf = nullptr;
    if (nosex_ct) {
      if (unlikely(
              bigstack_end_alloc_w(raw_sample_ctl, &nosex_buf) ||
              bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.nosex_interleaved_vec))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      BitvecInvmaskCopy(ctx.sample_include, sex_nm, raw_sample_ctl, nosex_buf);
      ZeroTrailingWords(raw_sample_ctl, nosex_buf);
      FillInterleavedMaskVec(nosex_buf, raw_sample_ctv, ctx.nosex_interleaved_vec);
    }

    ctx.variant_missing_hc_cts = variant_missing_hc_cts;
    ctx.variant_missing_dosage_cts = variant_missing_dosage_cts;
    ctx.variant_hethap_cts = variant_hethap_cts;
    ctx.first_hap_uidx = first_hap_uidx;
    ctx.is_minimac3_r2 = is_minimac3_r2;

    ctx.founder_info = nullptr;
    ctx.founder_info_interleaved_vec = nullptr;
    ctx.founder_info_cumulative_popcounts = nullptr;
    ctx.founder_male = nullptr;
    ctx.founder_male_interleaved_vec = nullptr;
    ctx.founder_male_cumulative_popcounts = nullptr;
    ctx.founder_nosex_interleaved_vec = nullptr;
    ctx.founder_ct = 0;
    ctx.founder_male_ct = 0;
    ctx.founder_nosex_ct = 0;
    ctx.founder_allele_ddosages = nullptr;
    ctx.founder_raw_geno_cts = nullptr;
    ctx.founder_x_male_geno_cts = nullptr;
    ctx.founder_x_nosex_geno_cts = nullptr;
    if (two_subsets_required) {
      if (founder_ct) {
        ctx.founder_info = founder_info;
        if (unlikely(
                bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_info_interleaved_vec) ||
                bigstack_alloc_u32(raw_sample_ctl, &ctx.founder_info_cumulative_popcounts) ||
                bigstack_alloc_w(raw_sample_ctl, &ctx.founder_male) ||
                bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_male_interleaved_vec) ||
                bigstack_alloc_u32(raw_sample_ctl, &ctx.founder_male_cumulative_popcounts))) {
          goto LoadAlleleAndGenoCounts_ret_NOMEM;
        }
        FillInterleavedMaskVec(founder_info, raw_sample_ctv, ctx.founder_info_interleaved_vec);
        FillCumulativePopcounts(founder_info, raw_sample_ctl, ctx.founder_info_cumulative_popcounts);
        BitvecAndCopy(sex_male, founder_info, raw_sample_ctl, ctx.founder_male);
        ZeroTrailingWords(raw_sample_ctl, ctx.founder_male);
        FillInterleavedMaskVec(ctx.founder_male, raw_sample_ctv, ctx.founder_male_interleaved_vec);
        FillCumulativePopcounts(ctx.founder_male, raw_sample_ctl, ctx.founder_male_cumulative_popcounts);
        ctx.founder_ct = founder_ct;
        ctx.founder_male_ct = ctx.founder_male_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(ctx.founder_male[raw_sample_ctl - 1]);
        ctx.founder_allele_ddosages = founder_allele_ddosages;
        ctx.founder_raw_geno_cts = founder_raw_geno_cts;
        ctx.founder_x_male_geno_cts = founder_x_male_geno_cts;
        if (nosex_ct) {
          // caller currently responsible for ensuring that when
          // founder_nosex_ct is zero, founder_x_nosex_geno_cts ==
          // nullptr
          if (unlikely(bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &ctx.founder_nosex_interleaved_vec))) {
            goto LoadAlleleAndGenoCounts_ret_NOMEM;
          }
          BitvecAnd(founder_info, raw_sample_ctl, nosex_buf);
          ctx.founder_nosex_ct = PopcountWords(nosex_buf, raw_sample_ctl);
          assert(ctx.founder_nosex_ct);
          ZeroTrailingWords(raw_sample_ctl, nosex_buf);
          FillInterleavedMaskVec(nosex_buf, raw_sample_ctv, ctx.founder_nosex_interleaved_vec);
          ctx.founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
        }
      } else {
        if (founder_allele_ddosages) {
          ZeroU64Arr(allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct), founder_allele_ddosages);
        }
        if (founder_raw_geno_cts) {
          memset(founder_raw_geno_cts, 0, raw_variant_ct * (3 * sizeof(int32_t)));
        }
      }
    } else if (founder_ct == sample_ct) {
      // bugfix: some founder and some nonfounder counts required
      if ((!ctx.allele_ddosages) && founder_allele_ddosages) {
        ctx.allele_ddosages = founder_allele_ddosages;
      }
      if ((!ctx.raw_geno_cts) && founder_raw_geno_cts) {
        ctx.raw_geno_cts = founder_raw_geno_cts;
      }
      if ((!ctx.x_male_geno_cts) && founder_x_male_geno_cts) {
        ctx.x_male_geno_cts = founder_x_male_geno_cts;
      }
      if ((!ctx.x_nosex_geno_cts) && founder_x_nosex_geno_cts) {
        ctx.x_nosex_geno_cts = founder_x_nosex_geno_cts;
      }
    } else if (only_founder_cts_required) {
      ctx.male_ct = ctx.sex_male_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(ctx.sex_male[raw_sample_ctl - 1]);
      if (nosex_ct) {
        ctx.nosex_ct = PopcountWords(nosex_buf, raw_sample_ctl);
      }
    }
    const uintptr_t raw_allele_ct = allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct);
    if (!ctx.sample_ct) {
      if (allele_presents) {
        ZeroWArr(BitCtToWordCt(raw_allele_ct), allele_presents);
      }
      if (ctx.allele_ddosages) {
        ZeroU64Arr(raw_allele_ct, ctx.allele_ddosages);
      }
      if (ctx.raw_geno_cts) {
        memset(ctx.raw_geno_cts, 0, raw_variant_ct * (3 * sizeof(int32_t)));
      }
      // early exit
      goto LoadAlleleAndGenoCounts_ret_1;
    }
    BigstackEndReset(bigstack_end_mark);  // free nosex_buf
    if (allele_presents) {
      const uintptr_t raw_allele_ct_a64 = RoundUpPow2(raw_allele_ct, kCacheline);
      if (unlikely(bigstack_left() < raw_allele_ct_a64)) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      // fill byte-array instead of bitarray so multithreading works
      ctx.allele_presents_bytearr = S_CAST(unsigned char*, bigstack_alloc_raw(raw_allele_ct_a64));
      memset(ctx.allele_presents_bytearr, 0, raw_allele_ct_a64);
    } else {
      ctx.allele_presents_bytearr = nullptr;
    }

    uint32_t unused_chr_code;
    uint32_t unused_chr_code2;
    uint32_t xy_complications_present = ((allele_presents || allele_ddosages || founder_allele_ddosages || variant_missing_dosage_cts) && XymtExists(cip, kChrOffsetX, &unused_chr_code)) || (allele_presents && (sample_ct != male_ct) && XymtExists(cip, kChrOffsetY, &unused_chr_code2));
    const uint32_t xy_dosages_needed = (pgfip->gflags & kfPgenGlobalDosagePresent) && xy_complications_present;

    // todo: check when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    const uint32_t max_allele_ct = pgfip->max_allele_ct;
    uint32_t mhc_needed = 0;
    ctx.thread_read_mhc = nullptr;
    if ((max_allele_ct > 2) && (variant_hethap_cts || allele_presents || allele_ddosages || founder_allele_ddosages || imp_r2_vals)) {
      if (unlikely(
              bigstack_alloc_u64p(calc_thread_ct, &ctx.all_dosages))) {
        goto LoadAlleleAndGenoCounts_ret_NOMEM;
      }
      mhc_needed = (xy_complications_present || ((variant_hethap_cts || imp_r2_vals) && XymtExists(cip, kChrOffsetX, &unused_chr_code)));
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        // double allocation size, to leave room for chrY ssqs
        if (unlikely(
                bigstack_alloc_u64(max_allele_ct * 2, &(ctx.all_dosages[tidx])))) {
          goto LoadAlleleAndGenoCounts_ret_NOMEM;
        }
      }
    } else {
      ctx.all_dosages = nullptr;
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    // defensive
    ctx.dosage_presents = nullptr;
    ctx.dosage_mains = nullptr;
    uint32_t read_block_size;
    // todo: check if raw_sample_ct should be replaced with sample_ct here
    if (unlikely(PgenMtLoadInit(variant_include, raw_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, mhc_needed? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, xy_dosages_needed? (&ctx.dosage_presents) : nullptr, xy_dosages_needed? (&ctx.dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto LoadAlleleAndGenoCounts_ret_NOMEM;
    }
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.reterr = kPglRetSuccess;
    SetThreadFuncAndData(LoadAlleleAndGenoCountsThread, &ctx, &tg);

    logputs("Calculating allele frequencies... ");
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto LoadAlleleAndGenoCounts_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = ctx.reterr;
        if (unlikely(reterr)) {
          goto LoadAlleleAndGenoCounts_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
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
      variant_idx += cur_block_size;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (allele_presents) {
      const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
      allele_presents[raw_allele_ctl - 1] = 0;
#ifdef __LP64__
      const uintptr_t vec_ct = DivUp(raw_allele_ct, kBytesPerVec);
      VecUc* bytearr_alias = R_CAST(VecUc*, ctx.allele_presents_bytearr);
      Vec8thUint* allele_presents_alias = R_CAST(Vec8thUint*, allele_presents);
      for (uintptr_t vec_idx = 0; vec_idx != vec_ct; ++vec_idx) {
        allele_presents_alias[vec_idx] = vecuc_movemask(bytearr_alias[vec_idx]);
      }
#else
      const uintptr_t twovec_ct = DivUp(raw_allele_ct, 8);
      uintptr_t* bytearr_iter = R_CAST(uintptr_t*, ctx.allele_presents_bytearr);
      unsigned char* allele_presents_iter = R_CAST(unsigned char*, allele_presents);
      unsigned char* allele_presents_stop = &(allele_presents_iter[twovec_ct]);
      for (; allele_presents_iter != allele_presents_stop; ++allele_presents_iter) {
        // 31,23,15,7 -> 3,2,1,0: multiply by number with bits 0,7,14,21 set,
        // then right-shift
        uintptr_t cur_word = ((*bytearr_iter++) * 0x204081) >> 28;
        cur_word |= ((*bytearr_iter++) * 0x204081) >> 24;
        *allele_presents_iter = cur_word;
      }
#endif
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  LoadAlleleAndGenoCounts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadAlleleAndGenoCounts_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  LoadAlleleAndGenoCounts_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 LoadAlleleAndGenoCounts_ret_1:
  CleanupThreads(&tg);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

void ApplyHardCallThresh(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec) {
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    const uint32_t dosage_int = dosage_main[dosage_idx];
    const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
    const uintptr_t widx = sample_uidx / kBitsPerWordD2;
    uintptr_t prev_geno_word = genovec[widx];
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
      genovec[widx] = prev_geno_word ^ (geno_xor << shift);
    }
  }
}

void FillMissingFromBiallelicDosage(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uintptr_t* genovec) {
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    const uintptr_t widx = sample_uidx / kBitsPerWordD2;
    const uint32_t shift = 2 * (sample_uidx % kBitsPerWordD2);
    uintptr_t geno_word = genovec[widx];
    if (((geno_word >> shift) & 3) != 3) {
      // not missing, don't recalculate
      continue;
    }
    const uint32_t dosage_int = dosage_main[dosage_idx];
    // Ties are rounded in favor of the lower-index allele, which reduces to
    // rounding down in the biallelic case.
    // uintptr_t new_geno = (dosage_int + kDosage4th - 1) / kDosageMid;
    const uintptr_t three_minus_new_geno = (7 * kDosage4th - dosage_int) / kDosageMid;
    geno_word ^= three_minus_new_geno << shift;
    genovec[widx] = geno_word;
  }
}

uint32_t ApplyHardCallThreshPhased(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uint32_t hard_call_halfdist, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uintptr_t* dphase_present, SDosage* dphase_delta, SDosage* tmp_dphase_delta) {
  // Generate new hphase values when we're converting a hardcall from
  // missing/hom to het, and abs(dphase_delta) > 0.5.  Erase explicit dphase in
  // that case if dphase_delta is maximal.
  //
  // Erase hphase value when we're converting a hardcall from het to
  // missing/hom.  If hardcall was previously phased and no explicit dphase
  // value existed, add it.
  //
  // Since both insertions and deletions are possible, we write the updated
  // dphase_delta to a buffer and copy it back, instead of editing in place.
  //
  // Returns final dphase_ct.
  //
  // Some extraneous phaseinfo bits may be set on return.
  const SDosage* dphase_read_iter = dphase_delta;
  SDosage* dphase_write_iter = tmp_dphase_delta;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    const uint32_t dosage_int = dosage_main[dosage_idx];
    const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
    const uintptr_t widx = sample_uidx / kBitsPerWordD2;
    uintptr_t prev_geno_word = genovec[widx];
    const uint32_t shift = (sample_uidx % kBitsPerWordD2) * 2;
    uintptr_t new_geno;
    if (halfdist < hard_call_halfdist) {
      new_geno = 3;
    } else {
      new_geno = (dosage_int + kDosage4th) / kDosageMid;
    }
    const uintptr_t prev_geno = (prev_geno_word >> shift) & 3;
    const uintptr_t geno_xor = new_geno ^ prev_geno;
    const uint32_t cur_hphase_present = IsSet(phasepresent, sample_uidx);
    if (IsSet(dphase_present, sample_uidx)) {
      int32_t dphase_delta_val = *dphase_read_iter++;
      *dphase_write_iter++ = dphase_delta_val;
      if (geno_xor) {
        if (new_geno == 1) {
          const uint32_t neg_sign_bit = -(S_CAST(uint32_t, dphase_delta_val) >> 31);
          const uint32_t abs_dphase_delta_val = (S_CAST(uint32_t, dphase_delta_val) ^ neg_sign_bit) - neg_sign_bit;
          if (abs_dphase_delta_val > kDosage4th) {
            SetBit(sample_uidx, phasepresent);
            AssignBit(sample_uidx, neg_sign_bit + 1, phaseinfo);
            // is dphase_delta maximal?
            if ((abs_dphase_delta_val == dosage_int) || (abs_dphase_delta_val + dosage_int == kDosageMax)) {
              ClearBit(sample_uidx, dphase_present);
              --dphase_write_iter;
            }
          }
        } else {
          ClearBit(sample_uidx, phasepresent);
        }
        genovec[widx] = prev_geno_word ^ (geno_xor << shift);
      }
    } else {
      if (geno_xor) {
        if (cur_hphase_present) {
          assert(new_geno != 1);
          ClearBit(sample_uidx, phasepresent);
          SetBit(sample_uidx, dphase_present);
          int32_t new_dphase_delta_val = DosageHomdist(dosage_int);
          if (!IsSet(phaseinfo, sample_uidx)) {
            new_dphase_delta_val = -new_dphase_delta_val;
          }
          *dphase_write_iter++ = new_dphase_delta_val;
        }
        genovec[widx] = prev_geno_word ^ (geno_xor << shift);
      }
    }
  }
  const uint32_t dphase_ct = dphase_write_iter - tmp_dphase_delta;
  memcpy(dphase_delta, tmp_dphase_delta, dphase_ct * sizeof(Dosage));
  return dphase_ct;
}

uint32_t FillMissingFromBiallelicDosagePhased(const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t dosage_ct, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uintptr_t* dphase_present, SDosage* dphase_delta) {
  // Generate new hphase values when we're creating a new het hardcall, and
  // abs(dphase_delta) > 0.5.  Erase explicit dphase in that case if
  // dphase_delta is maximal.
  //
  // Only deletions to dphase_delta[] are possible, so we edit in place.
  //
  // Returns final dphase_ct.
  uint32_t dphase_read_idx = 0;
  uint32_t dphase_write_idx = 0;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  uint32_t dosage_idx = 0;
  for (; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    const uintptr_t widx = sample_uidx / kBitsPerWordD2;
    const uint32_t shift = 2 * (sample_uidx % kBitsPerWordD2);
    uintptr_t geno_word = genovec[widx];
    if (((geno_word >> shift) & 3) != 3) {
      if (IsSet(dphase_present, sample_uidx)) {
        if (dphase_read_idx != dphase_write_idx) {
          dphase_delta[dphase_write_idx] = dphase_delta[dphase_read_idx];
        }
        ++dphase_read_idx;
        ++dphase_write_idx;
      }
      continue;
    }
    const uint32_t dosage_int = dosage_main[dosage_idx];
    const uintptr_t three_minus_new_geno = (7 * kDosage4th - dosage_int) / kDosageMid;
    geno_word ^= three_minus_new_geno << shift;
    genovec[widx] = geno_word;
    if (!IsSet(dphase_present, sample_uidx)) {
      continue;
    }
    const int32_t dphase_delta_val = dphase_delta[dphase_read_idx];
    if (dphase_write_idx != dphase_read_idx) {
      dphase_delta[dphase_write_idx] = dphase_delta_val;
    }
    ++dphase_read_idx;
    ++dphase_write_idx;
    if (three_minus_new_geno != 2) {
      continue;
    }
    const uint32_t neg_sign_bit = -(S_CAST(uint32_t, dphase_delta_val) >> 31);
    const uint32_t abs_dphase_delta_val = (S_CAST(uint32_t, dphase_delta_val) ^ neg_sign_bit) - neg_sign_bit;
    if (abs_dphase_delta_val > kDosage4th) {
      SetBit(sample_uidx, phasepresent);
      AssignBit(sample_uidx, neg_sign_bit + 1, phaseinfo);
      // is dphase_delta maximal?
      if ((abs_dphase_delta_val == dosage_int) || (abs_dphase_delta_val + dosage_int == kDosageMax)) {
        ClearBit(sample_uidx, dphase_present);
        --dphase_write_idx;
      }
    }
  }
  return dphase_write_idx;
}

uintptr_t InitWriteAlleleIdxOffsets(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uint32_t* new_variant_idx_to_old, uint32_t variant_ct, uintptr_t* new_allele_idx_offsets) {
  uintptr_t cur_offset = 0;
  if (allele_presents) {
    uint32_t ref_allele_idx = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = 0;
    if (!new_variant_idx_to_old) {
      cur_bits = variant_include[0];
    }
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      uint32_t variant_uidx;
      if (new_variant_idx_to_old) {
        variant_uidx = new_variant_idx_to_old[variant_idx];
      } else {
        variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      }
      new_allele_idx_offsets[variant_idx] = cur_offset;
      const uintptr_t old_offset_start = allele_idx_offsets[variant_uidx];
      const uintptr_t old_offset_end = allele_idx_offsets[variant_uidx + 1];
      uint32_t cur_allele_ct = old_offset_end - old_offset_start;
      if (cur_allele_ct > 2) {
        cur_allele_ct = PopcountBitRange(allele_presents, old_offset_start, old_offset_end);
        if (refalt1_select) {
          ref_allele_idx = refalt1_select[variant_uidx][0];
        }
        if (!IsSet(allele_presents, old_offset_start + ref_allele_idx)) {
          ++cur_allele_ct;
        }
        if (cur_allele_ct < 2) {
          cur_allele_ct = 2;
        }
      }
      cur_offset += cur_allele_ct;
    }
  } else if (!new_variant_idx_to_old) {
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      new_allele_idx_offsets[variant_idx] = cur_offset;
      cur_offset += allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
    }
  } else {
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = new_variant_idx_to_old[variant_idx];
      new_allele_idx_offsets[variant_idx] = cur_offset;
      cur_offset += allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
    }
  }
  return cur_offset;
}

// Join behavior:
// - Require sorted input .pvar for now, though it won't be difficult to lift
//   this restriction later.  Ok for input to contain multiallelic variants.
// - Don't need to do anything different when chr:pos only appears once in a
//   biallelic variant, or in a multiallelic variant in "+any" mode.  For
//   multiallelic variants in +both/+snps mode, error out if the variant is
//   mixed SNP/non-SNP.  In +both case, also error out if the variant is mixed
//   symbolic/non-symbolic, or it's symbolic and does not satisfy the REF or
//   INFO:END constraints.
// - When multiple variants have the same chr:pos:
//   - In +snps mode, also don't need to do anything different for non-SNPs.
//   - Otherwise, create up to three linked lists of input variant records
//     (need to distinguish SNPs from non-SNPs in +both mode, and within the
//     non-SNP category, symbolic alleles are separate from non-symbolic).
//     The variant with symbolic alleles has additional constraints: a warning
//     is printed if INFO:END isn't defined, and an error occurs if either REF
//     is multi-character, or there's an INFO:END mismatch.
//   - Error out if REF alleles aren't all consistent, or any ALT allele is
//     duplicated (note that --pmerge must support the latter).
//   - For joined not-entirely-SNP non-symbolic variants, the final REF is the
//     longest of the original REFs; ALT alleles have bases added to the end if
//     necessary.  (Yes, this causes SNPs to stop being visible to a strlen ==
//     1 check in +any mode, which is why + is interpreted as +both instead.)
//   - Final ALT allele order is based on allele frequency (highest first),
//     with ties broken by natural-sort.
//   - ID: 1. If --set-all-var-ids specified, apply template.
//         2. Otherwise, keep original variant ID if all sources identical and
//            nonmissing.
//         3. Otherwise, if vid-join specified, check if all original IDs are
//            nonmissing and contain exactly one ';' per extra ALT allele.  If
//            so, join on ';' in final ALT allele order.
//         4. Otherwise, set to --set-missing-var-ids template or missing code.
//   - QUAL is minimum of inputs ('.' treated as positive infinity).  FILTER is
//     natural-sorted union of non-PASS values; if all inputs were '.'/PASS,
//     output is PASS unless all inputs were missing.
//   - INFO join is based on the Number field in the key's header line.  (If
//     there's no header line corresponding to a key, we error out.)  For
//     Number=A, we join in the obvious manner.  For Number=R (or Number=G in
//     the haploid case), we replace the reference allele entry with . iff
//     there's any string mismatch (e.g. '13' and '13.0' will be treated as
//     unequal), and print a warning (with more than 3 warnings, later warnings
//     are only written to log file).  For diploid Number=G, we do the same for
//     the hom-ref entry.  For Number=0, we error out if the header line Type
//     isn't Flag, and the final flag is set iff any of the original variants
//     have the flag set.  For the other cases (Number=<fixed constant> or
//     '.'), we replace the value with '.' if there's any string mismatch.
//     The key won't appear at all iff it doesn't appear in any of the original
//     variants (not even a '.').
//   - We error out if the joined variant has total ALT dosage > ~2.02, or
//     either side of the total ALT phased dosage > ~1.01.  (We scale the
//     components down if there's a <1% overflow.)
// - Missing ALTs... ugh.  Have to permit this in biallelic case, but don't
//   allow it elsewhere.  This forces a SNP (if REF is single-char) and/or
//   non-SNP (if REF is multichar) entry to be written when no corresponding
//   regular ALT is present, genotype writing/merging will error out if a
//   missing allele has any dosage, and QUAL/FILTER/INFO is merged as usual
//   when there are multiple missing-ALT same-type variants at the same
//   position.  However, when a regular ALT is present, the missing-ALT
//   variants are completely ignored.

// The simplest design involves precomputing the entire (new variant idx, new
// allele idx) -> (old variant uidx, old allele idx) mapping here, and
// referring to that in both the .pvar and .pgen writers.  Unfortunately, that
// has a rather high memory requirement of 5 bytes per allele (assuming
// sizeof(AlleleCode) == 1).  While that's smaller than the 8 bytes/allele we
// pay for allele_storage[], and also practically always smaller than the
// 21+[variant ID len] per variant we pay for variant_bps + variant_ids +
// allele_idx_offsets, it's still worth some effort to avoid; in particular, we
// want an 8 GiB workspace to be sufficient for most operations on the full
// 1000 Genomes phase 3 variant set (~84.8 million), and that's barely true
// right now, so a bit of additional complexity to avoid losing ~850 MB is
// justified.
//
// Thus, we only save the number of alleles in each new variant here, and force
// the which-allele-comes-from-where computation (as well as SNP vs. non-SNP
// vs. symbolic handling) to be repeated in the .pvar and .pgen writers.  This
// sucks, but being forced to perform an ordinary analysis on a remote machine
// rather than locally sucks a bit more.
//
// Incidentally, another place to look, if it's important to further reduce
// memory requirements, is internal representation of the FILTER field.  The
// current design gains a bit of extra speed by simply storing the non-./PASS
// strings without parsing them further; but we could put them into a temporary
// storage location and then convert to bitarray + string table at the end of
// LoadPvar().  (Note that we already use only bitarrays when all FILTER values
// are ./PASS, though.)

ENUM_U31_DEF_START()
  kJoinVtypeError,
  kJoinVtypeSnp,
  kJoinVtypeNonsnp,
  kJoinVtypeMixedSnpNonsnp,
  kJoinVtypeSymbolic,
  kJoinVtypeEnd
ENUM_U31_DEF_END(JoinVtype);

typedef struct JoinCountsStruct {
  uintptr_t snp_ct;
  uintptr_t nonsnp_ct;
  uintptr_t symbolic_ct;
  uint32_t missalt_snp_ct;
  uint32_t missalt_nonsnp_ct;
} JoinCounts;

JoinVtype JoinCount(const char* const* cur_alleles, uintptr_t allele_ct, JoinCounts* jcp) {
  jcp->snp_ct = 0;
  jcp->symbolic_ct = 0;
  jcp->missalt_snp_ct = 0;
  jcp->missalt_nonsnp_ct = 0;
  if (cur_alleles[0][1] == '\0') {
    jcp->nonsnp_ct = 0;
    for (uintptr_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
      const char* cur_allele = cur_alleles[allele_idx];
      if (cur_allele[0] == '<') {
        jcp->symbolic_ct += 1;
      } else if (cur_allele[1] == '\0') {
        if (cur_allele[0] == '.') {
          if (allele_ct == 2) {
            jcp->missalt_snp_ct = 1;
            return kJoinVtypeSnp;
          }
          return kJoinVtypeError;
        }
        jcp->snp_ct += 1;
      } else {
        jcp->nonsnp_ct += 1;
      }
    }
    if (jcp->symbolic_ct) {
      return (jcp->symbolic_ct == allele_ct - 1)? kJoinVtypeSymbolic : kJoinVtypeError;
    }
    if (jcp->nonsnp_ct) {
      return jcp->snp_ct? kJoinVtypeMixedSnpNonsnp : kJoinVtypeNonsnp;
    }
    return kJoinVtypeSnp;
  }
  for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
    const char* cur_allele = cur_alleles[allele_idx];
    if (cur_allele[0] == '<') {
      return kJoinVtypeError;
    }
    if (memequal_k(cur_allele, ".", 2)) {
      if (allele_ct == 2) {
        jcp->nonsnp_ct = 0;
        jcp->missalt_nonsnp_ct = 1;
        return kJoinVtypeNonsnp;
      }
      return kJoinVtypeError;
    }
  }
  jcp->nonsnp_ct = allele_ct - 1;
  return kJoinVtypeNonsnp;
}

void PlanJoinOne(uint32_t cur_alt_allele_ct, uintptr_t** write_allele_idx_offsets_iterp, uintptr_t* cur_offsetp, uint32_t* max_write_allele_ctp) {
  const uint32_t cur_write_allele_ct = 1 + MAXV(1, cur_alt_allele_ct);
  if (cur_write_allele_ct > (*max_write_allele_ctp)) {
    *max_write_allele_ctp = cur_write_allele_ct;
  }
  *cur_offsetp += cur_write_allele_ct;
  uintptr_t* write_allele_idx_offsets_iter = *write_allele_idx_offsets_iterp;
  *write_allele_idx_offsets_iter++ = *cur_offsetp;
  *write_allele_idx_offsets_iterp = write_allele_idx_offsets_iter;
}

void PlanJoinFlushPos(const JoinCounts* jcp, MakePlink2Flags join_mode, uintptr_t** write_allele_idx_offsets_iterp, uintptr_t* cur_offsetp, uint32_t* max_write_allele_ctp, uint32_t* max_missalt_ctp) {
  if (join_mode == kfMakePlink2MJoinSnps) {
    if (!(jcp->snp_ct || jcp->missalt_snp_ct)) {
      // all non-SNPs at this position, which were already accounted for
      return;
    }
    PlanJoinOne(jcp->snp_ct, write_allele_idx_offsets_iterp, cur_offsetp, max_write_allele_ctp);
    if ((!jcp->snp_ct) && (jcp->missalt_snp_ct > (*max_missalt_ctp))) {
      *max_missalt_ctp = jcp->missalt_snp_ct;
    }
    return;
  }
  if (join_mode == kfMakePlink2MJoinBoth) {
    if (jcp->snp_ct || jcp->missalt_snp_ct) {
      PlanJoinOne(jcp->snp_ct, write_allele_idx_offsets_iterp, cur_offsetp, max_write_allele_ctp);
      if ((!jcp->snp_ct) && (jcp->missalt_snp_ct > (*max_missalt_ctp))) {
        *max_missalt_ctp = jcp->missalt_snp_ct;
      }
    }
    if (jcp->nonsnp_ct || jcp->missalt_nonsnp_ct) {
      PlanJoinOne(jcp->nonsnp_ct, write_allele_idx_offsets_iterp, cur_offsetp, max_write_allele_ctp);
      if ((!jcp->nonsnp_ct) && (jcp->missalt_nonsnp_ct > (*max_missalt_ctp))) {
        *max_missalt_ctp = jcp->missalt_nonsnp_ct;
      }
    }
  } else {
    if (jcp->snp_ct || jcp->nonsnp_ct || jcp->missalt_snp_ct || jcp->missalt_nonsnp_ct) {
      const uint32_t alt_allele_ct = jcp->snp_ct + jcp->nonsnp_ct;
      PlanJoinOne(alt_allele_ct, write_allele_idx_offsets_iterp, cur_offsetp, max_write_allele_ctp);
      const uint32_t missalt_ct = jcp->missalt_snp_ct + jcp->missalt_nonsnp_ct;
      if ((missalt_ct > (*max_missalt_ctp)) && (!alt_allele_ct)) {
        *max_missalt_ctp = missalt_ct;
      }
    }
  }
  if (jcp->symbolic_ct) {
    PlanJoinOne(jcp->symbolic_ct, write_allele_idx_offsets_iterp, cur_offsetp, max_write_allele_ctp);
  }
}


// *write_allele_idx_offsetsp is assumed to be initialized to nullptr.
// *max_missalt_ctp is assumed to be initialized to 0.
PglErr PlanMultiallelicJoin(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, MakePlink2Flags flags, uint32_t* write_variant_ctp, const uintptr_t** write_allele_idx_offsetsp, uint32_t* max_write_allele_ctp, uint32_t* max_missalt_ctp) {
  uint32_t variant_uidx = 0;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t variant_ct = *write_variant_ctp;
    uintptr_t* write_allele_idx_offsets = R_CAST(uintptr_t*, g_bigstack_base);
    uintptr_t* write_allele_idx_offsets_stop = R_CAST(uintptr_t*, RoundDownPow2(R_CAST(uintptr_t, g_bigstack_end), kCacheline));
    if (write_allele_idx_offsets == write_allele_idx_offsets_stop) {
      goto PlanMultiallelicJoin_ret_NOMEM;
    }
    write_allele_idx_offsets_stop = &(write_allele_idx_offsets_stop[-4]);
    const MakePlink2Flags join_mode = flags & kfMakePlink2MMask;
    uintptr_t* write_allele_idx_offsets_iter = write_allele_idx_offsets;
    *write_allele_idx_offsets_iter++ = 0;
    uintptr_t cur_offset = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t prev_bp = 0;
    uint32_t allele_ct = 2;
    uint32_t max_write_allele_ct = 2;
    JoinCounts jc;
    // possible todo: track max_write_allele_ct for each subcase, instead of
    // having a single value
    jc.snp_ct = 0;
    jc.nonsnp_ct = 0;
    jc.symbolic_ct = 0;
    jc.missalt_snp_ct = 0;
    jc.missalt_nonsnp_ct = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        prev_bp = UINT32_MAX;
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      JoinCounts cur_jc;
      JoinVtype jvt = JoinCount(cur_alleles, allele_ct, &cur_jc);
      if (unlikely(jvt == kJoinVtypeError)) {
        goto PlanMultiallelicJoin_ret_MIXED_SYMBOLIC;
      }
      if (unlikely((join_mode != kfMakePlink2MJoinAny) && (jvt == kJoinVtypeMixedSnpNonsnp))) {
        logerrprintfww("Error: Variant '%s' is mixed SNP/non-SNP; multiallelics=+both and +snps don't permit this.\n", variant_ids[variant_uidx]);
        goto PlanMultiallelicJoin_ret_INCONSISTENT_INPUT;
      }
      const uint32_t cur_bp = variant_bps[variant_uidx];
      if (cur_bp != prev_bp) {
        PlanJoinFlushPos(&jc, join_mode, &write_allele_idx_offsets_iter, &cur_offset, &max_write_allele_ct, max_missalt_ctp);
        if (join_mode == kfMakePlink2MJoinSnps) {
          if (cur_jc.nonsnp_ct || cur_jc.symbolic_ct || cur_jc.missalt_nonsnp_ct) {
            // Flush non-SNP immediately.
            const uint32_t cur_write_allele_ct = 1 + cur_jc.nonsnp_ct + cur_jc.symbolic_ct;
            cur_offset += cur_write_allele_ct;
            if (cur_write_allele_ct > max_write_allele_ct) {
              max_write_allele_ct = cur_write_allele_ct;
            }
            *write_allele_idx_offsets_iter++ = cur_offset;
            // Also need to reinitialize.
            jc.snp_ct = 0;
            jc.missalt_snp_ct = 0;
          } else {
            jc.snp_ct = cur_jc.snp_ct;
            jc.missalt_snp_ct = cur_jc.missalt_snp_ct;
          }
        } else {
          jc = cur_jc;
        }
        prev_bp = cur_bp;
      } else if ((join_mode == kfMakePlink2MJoinSnps) && (cur_jc.nonsnp_ct || cur_jc.symbolic_ct)) {
        // Flush non-SNP immediately.
        const uint32_t cur_write_allele_ct = 1 + cur_jc.nonsnp_ct + cur_jc.symbolic_ct;
        cur_offset += cur_write_allele_ct;
        if (cur_write_allele_ct > max_write_allele_ct) {
          max_write_allele_ct = cur_write_allele_ct;
        }
        *write_allele_idx_offsets_iter++ = cur_offset;
      } else {
        jc.snp_ct += cur_jc.snp_ct;
        jc.nonsnp_ct += cur_jc.nonsnp_ct;
        jc.symbolic_ct += cur_jc.symbolic_ct;
        jc.missalt_snp_ct += cur_jc.missalt_snp_ct;
        jc.missalt_nonsnp_ct += cur_jc.missalt_nonsnp_ct;
        continue;
      }
      if (write_allele_idx_offsets_iter > write_allele_idx_offsets_stop) {
        goto PlanMultiallelicJoin_ret_NOMEM;
      }
    }
    // Flush last position.
    PlanJoinFlushPos(&jc, join_mode, &write_allele_idx_offsets_iter, &cur_offset, &max_write_allele_ct, max_missalt_ctp);
    if (max_write_allele_ct > kPglMaxAltAlleleCt + 1) {
      goto PlanMultiallelicJoin_ret_TOO_MANY_ALTS;
    }
    *write_variant_ctp = S_CAST(uintptr_t, write_allele_idx_offsets_iter - write_allele_idx_offsets) - 1;
    *max_write_allele_ctp = max_write_allele_ct;
    if (max_write_allele_ct > 2) {
      BigstackBaseSet(write_allele_idx_offsets_iter);
      *write_allele_idx_offsetsp = write_allele_idx_offsets;
    }
  }
  while (0) {
  PlanMultiallelicJoin_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PlanMultiallelicJoin_ret_MIXED_SYMBOLIC:
    logerrprintfww("Error: Variant '%s' mixes symbolic and non-symbolic alleles in an unsupported manner.\n", variant_ids[variant_uidx]);
  PlanMultiallelicJoin_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  PlanMultiallelicJoin_ret_TOO_MANY_ALTS:
    logerrprintf("Error: Variant-join would create a variant with too many ALT alleles for this\nplink2 build.\n");
    reterr = kPglRetNotYetSupported;
    break;
  }
  return reterr;
}

PglErr PlanMultiallelicSplit(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t max_allele_ct, MakePlink2Flags flags, uint32_t* write_variant_ctp, const uintptr_t** write_allele_idx_offsetsp) {
  uint32_t variant_uidx = 0;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t variant_ct = *write_variant_ctp;
    const uint32_t only_split_snps = ((flags & kfMakePlink2MMask) == kfMakePlink2MSplitSnps);
    uintptr_t* write_allele_idx_offsets = nullptr;
    uintptr_t* write_allele_idx_offsets_stop = nullptr;
    uintptr_t* write_allele_idx_offsets_iter = nullptr;
    if (only_split_snps) {
      write_allele_idx_offsets = R_CAST(uintptr_t*, g_bigstack_base);
      write_allele_idx_offsets_stop = R_CAST(uintptr_t*, RoundDownPow2(R_CAST(uintptr_t, g_bigstack_end), kCacheline));
      if (S_CAST(uintptr_t, write_allele_idx_offsets_stop - write_allele_idx_offsets) <= max_allele_ct) {
        goto PlanMultiallelicSplit_ret_NOMEM;
      }
      write_allele_idx_offsets_stop -= max_allele_ct;
      write_allele_idx_offsets_iter = write_allele_idx_offsets;
      *write_allele_idx_offsets_iter++ = 0;
    }
    uintptr_t cur_offset = 0;
    uintptr_t write_variant_ct = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      const uintptr_t allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      const uint32_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      if (allele_ct == 2) {
        if (only_split_snps) {
          cur_offset += 2;
          *write_allele_idx_offsets_iter++ = cur_offset;
        }
        ++write_variant_ct;
      } else {
        if (only_split_snps) {
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t do_split = 1;
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            if (cur_alleles[allele_idx][1] != '\0') {
              do_split = 0;
              break;
            }
          }
          if (do_split) {
            for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
              cur_offset += 2;
              *write_allele_idx_offsets_iter++ = cur_offset;
            }
            write_variant_ct += allele_ct - 1;
          } else {
            cur_offset += allele_ct;
            *write_allele_idx_offsets_iter++ = cur_offset;
            ++write_variant_ct;
          }
          if (write_allele_idx_offsets_iter > write_allele_idx_offsets_stop) {
            goto PlanMultiallelicSplit_ret_NOMEM;
          }
        } else {
          write_variant_ct += allele_ct - 1;
        }
      }
    }
    if (write_variant_ct > 0x7ffffffd) {
      logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
      goto PlanMultiallelicSplit_ret_INCONSISTENT_INPUT;
    }
    *write_variant_ctp = write_variant_ct;
    if (only_split_snps && (cur_offset != 2 * write_variant_ct)) {
      assert(cur_offset > 2 * write_variant_ct);
      BigstackBaseSet(write_allele_idx_offsets_iter);
      *write_allele_idx_offsetsp = write_allele_idx_offsets;
    }
  }
  while (0) {
  PlanMultiallelicSplit_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PlanMultiallelicSplit_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  return reterr;
}

// Returns 1 iff there are exactly (allele_ct - 2) semicolons in
// orig_variant_id, and no two are adjacent (or leading/trailing).
uint32_t VaridSplitOk(const char* orig_variant_id, uint32_t allele_ct) {
  const char* id_iter = orig_variant_id;
  for (uint32_t aidx = 2; aidx != allele_ct; ++aidx) {
    const char* tok_end = strchr(id_iter, ';');
    if ((!tok_end) || (tok_end == id_iter)) {
      return 0;
    }
    id_iter = &(tok_end[1]);
  }
  return (*id_iter != '\0') && (!strchr(id_iter, ';'));
}

// Similar to WriteMapOrBim(), but there are enough small differences to
// justify making this a separate function instead of clogging the original
// with more conditionals.
PglErr WriteBimSplit(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* variant_cms, const char* varid_template_str, const char* missing_varid_match, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t new_variant_id_max_allele_slen, uint32_t varid_split, uint32_t varid_dup, MiscFlags misc_flags, uint32_t output_zst, uint32_t thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WriteBimSplit_ret_NOMEM;
    }
    const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
    const uint32_t varid_dup_nosplit = varid_dup && (!varid_split);
    VaridTemplate* varid_templatep = nullptr;
    uint32_t missing_varid_slen = 0;
    uint32_t missing_varid_match_blen = 0; // nonzero iff --set-missing-var-ids
    if (varid_template_str) {
      if (!missing_varid_match) {
        missing_varid_match = &(g_one_char_strs[92]); // '.'
      }
      missing_varid_slen = strlen(missing_varid_match);
      if (misc_flags & kfMiscSetMissingVarIds) {
        missing_varid_match_blen = missing_varid_slen + 1;
      }
      if (unlikely(BIGSTACK_ALLOC_X(VaridTemplate, 1, &varid_templatep))) {
        goto WriteBimSplit_ret_NOMEM;
      }
      const uint32_t overflow_substitute_blen = new_variant_id_overflow_missing? (missing_varid_slen + 1) : 0;
      VaridTemplateInit(varid_template_str, missing_varid_match, chr_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, varid_templatep);
      if (varid_dup) {
        for (uint32_t uii = 0; uii != varid_templatep->insert_ct; ++uii) {
          const uint32_t insert_type = varid_templatep->insert_types[uii];
          if ((insert_type == 3) || ((insert_type == 2) && (varid_templatep->alleles_needed & 4))) {
            // Could define what takes precedence here, but simpler to prohibit
            // this combination.
            logerrputs("Error: 'vid-[split-]dup' cannot be used with a --set-all-var-ids or\n--set-missing-var-ids template string containing a non-REF allele.\n");
            goto WriteBimSplit_ret_INVALID_CMDLINE;
          }
        }
      }
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WriteBimSplit_ret_1;
    }

    const VaridTemplate* cur_varid_templatep = nullptr;
    const char* varid_token_start = nullptr; // for vid-split
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        char* chr_name_end = chrtoa(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end = '\t';
        const uint32_t chr_slen = chr_name_end - chr_buf;
        chr_buf_blen = 1 + chr_slen;
        if (varid_templatep) {
          const int32_t chr_slen_delta = chr_slen - varid_templatep->chr_slen;
          varid_templatep->chr_slen = chr_slen;
          varid_templatep->base_len += chr_slen_delta;
        }
      }
      const uintptr_t allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      const uint32_t orig_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      const char* orig_variant_id = variant_ids[variant_uidx];
      const char* ref_allele = cur_alleles[0];
      const uint32_t ref_allele_slen = strlen(ref_allele);
      uint32_t keep_orig_id = 1;
      if ((orig_allele_ct > 2) && (!varid_dup_nosplit)) {
        keep_orig_id = 0;
        if (varid_templatep && (!missing_varid_match_blen)) {
          cur_varid_templatep = varid_templatep;
        } else {
          cur_varid_templatep = nullptr;
          if (varid_split) {
            if (VaridSplitOk(orig_variant_id, orig_allele_ct)) {
              varid_token_start = orig_variant_id;
            } else if (varid_dup) {
              keep_orig_id = 1;
            } else {
              varid_token_start = nullptr;
            }
          }
          if ((!varid_token_start) && varid_templatep) {
            // --set-missing-var-ids usually applies here when it's specified;
            // the exceptions are when vid-split was also specified and the
            // split succeeded, or vid-split-dup was specified.
            // (In the latter case, this value is ignored anyway.)
            cur_varid_templatep = varid_templatep;
          }
        }
      }
      const uint32_t cur_bp = variant_bps[variant_uidx];
      // We already verified that no variants to be written have >2 alleles, so
      // we don't need to distinguish between '-' and '-snps'.
      for (uint32_t alt_allele_idx = 1; alt_allele_idx != orig_allele_ct; ++alt_allele_idx) {
        cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
        const char* cur_alt_allele = cur_alleles[alt_allele_idx];
        const uint32_t cur_alt_allele_slen = strlen(cur_alt_allele);
        if (keep_orig_id) {
          cswritep = strcpyax(cswritep, orig_variant_id, '\t');
        } else {
          if (cur_varid_templatep) {
            // Always true in --set-all-var-ids case.  True in
            // --set-missing-var-ids case when vid-split unspecified, or split
            // failed.
            cswritep = VaridTemplateWrite(cur_varid_templatep, ref_allele, cur_alt_allele, cur_bp, ref_allele_slen, 0, cur_alt_allele_slen, cswritep);
            *cswritep++ = '\t';
          } else if (varid_token_start) {
            const char* varid_token_end = strchrnul(varid_token_start, ';');
            // If substring matches missing code and --set-missing-var-ids is
            // specified, we replace it.
            if (varid_templatep && (S_CAST(uintptr_t, varid_token_end - varid_token_start) == missing_varid_slen) && memequal(varid_token_start, missing_varid_match, missing_varid_slen)) {
              cswritep = VaridTemplateWrite(varid_templatep, ref_allele, cur_alt_allele, cur_bp, ref_allele_slen, 0, cur_alt_allele_slen, cswritep);
            } else {
              cswritep = memcpya(cswritep, varid_token_start, varid_token_end - varid_token_start);
            }
            *cswritep++ = '\t';
            varid_token_start = &(varid_token_end[1]);
          } else {
            cswritep = memcpyax(cswritep, missing_varid_match, missing_varid_slen, '\t');
          }
        }
        if (!variant_cms) {
          *cswritep++ = '0';
        } else {
          cswritep = dtoa_g_p8(variant_cms[variant_uidx], cswritep);
        }
        *cswritep++ = '\t';
        cswritep = u32toa(cur_bp, cswritep);
        *cswritep++ = '\t';
        // note that VCF ref allele corresponds to A2, not A1
        cswritep = memcpya(cswritep, cur_alt_allele, cur_alt_allele_slen);
        *cswritep++ = '\t';
        cswritep = memcpya(cswritep, ref_allele, ref_allele_slen);
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto WriteBimSplit_ret_WRITE_FAIL;
        }
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WriteBimSplit_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WriteBimSplit_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteBimSplit_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WriteBimSplit_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 WriteBimSplit_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

// We only need to distinguish between the following INFO-value-type cases:
// Number=0 (flag), Number=<positive integer>, Number=., Number=A, Number=R,
// and Number=G.  We use negative numbers to represent the last 4 cases in
// InfoVtype.
CONSTI32(kInfoVtypeUnknown, -1);
CONSTI32(kInfoVtypeA, -2);
CONSTI32(kInfoVtypeR, -3);
CONSTI32(kInfoVtypeG, -4);

// Main fixed data structure when splitting/joining INFO is a hashmap of keys.
// Behavior when splitting:
// - Field order in the original variant is retained.
// - Number >= 0 and Number=. don't require any special handling, just copy the
//   entire key=value pair (or lone key, in the Flag case).
// - Number=A and Number=R require splitting the value on ',' and verifying the
//   comma count is correct, but is otherwise straightforward since alleles
//   can't be permuted.
// - Number=G requires a bit more work but isn't fundamentally different from
//   A/R.
// When joining:
// - Field order is determined by header line order.
// - Number=. and Number>0 just require a buffer of size ~info_reload_slen, and
//   a boolean indicating whether no mismatch has been found.
// - Number=0 (Flag) requires a single boolean, we perform an or operation.
// - Number=A/R/G are the messy ones: we need to have enough space for
//   max_write_allele_ct (or that minus 1) comma-separated values in the =A and
//   =R cases, and max_write_allele_ct * (max_write_allele_ct + 1) / 2 in the
//   diploid =G case.
//   Since we permit already-multiallelic variants to be part of a join, the =G
//   case may require a lot of working memory to handle.  We reserve up to 1/16
//   of remaining workspace memory for this when we cannot prove that we can
//   get by with less.

typedef struct InfoVtypeStruct {
  NONCOPYABLE(InfoVtypeStruct);
  int32_t num;
  char key[];
} InfoVtype;

// info_keys[] entries point to the (variable-size) key[] member of InfoVtype
// structs.  We use [const_]container_of(x)->num to look up the associated
// Number= value.
PglErr ParseInfoHeader(const char* xheader, uintptr_t xheader_blen, const char* const** info_keys_ptr, uint32_t* info_key_ctp, uint32_t** info_keys_htablep, uint32_t* info_keys_htable_sizep) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    // Parsing loop is similar to that in ExportVcf().
    const char* xheader_iter = xheader;
    const char* xheader_end = &(xheader[xheader_blen]);
    const char* line_end = xheader;
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    const char** info_keys = R_CAST(const char**, bigstack_mark);
    const char** info_keys_iter = info_keys;
    while (line_end != xheader_end) {
      xheader_iter = line_end;
      line_end = AdvPastDelim(xheader_iter, '\n');
      const uint32_t slen = line_end - xheader_iter;
      if ((slen <= 12) || (!StrStartsWithUnsafe(xheader_iter, "##INFO=<ID="))) {
        continue;
      }
      const char* key_start = &(xheader_iter[11]);
      const char* key_end = S_CAST(const char*, memchr(key_start, ',', slen - 12));
      if (unlikely((!key_end) || (!StrStartsWithUnsafe(key_end, ",Number=")))) {
        goto ParseInfoHeader_ret_MALFORMED_INFO_HEADER_LINE;
      }
      const uint32_t key_slen = key_end - key_start;
      if (key_slen > kMaxInfoKeySlen) {
        logerrputs("Error: " PROG_NAME_STR " does not support INFO keys longer than " MAX_INFO_KEY_SLEN_STR " characters.\n");
        // VCF spec doesn't specify a limit, so this isn't "malformed input".
        // We enforce a limit so we can safely print INFO keys in error
        // messages, etc.; it's trivial to increase the limit if it's ever
        // necessary.
        reterr = kPglRetNotYetSupported;
        goto ParseInfoHeader_ret_1;
      }
      const uintptr_t entry_byte_ct = RoundUpPow2(offsetof(InfoVtype, key) + 1 + key_slen, sizeof(intptr_t));
      if (S_CAST(uintptr_t, tmp_alloc_end - R_CAST(unsigned char*, info_keys_iter)) < entry_byte_ct + 8) {
        goto ParseInfoHeader_ret_NOMEM;
      }
      tmp_alloc_end -= entry_byte_ct;
      InfoVtype* new_entry = R_CAST(InfoVtype*, tmp_alloc_end);
      memcpyx(new_entry->key, key_start, key_slen, '\0');
      *info_keys_iter++ = new_entry->key;

      const char* num_iter = &(key_end[8]);
      const unsigned char first_num_char = num_iter[0];
      if (first_num_char < '1') {
        if (first_num_char == '0') {
          // don't see a reason to tolerate Number=01, etc.
          if (unlikely(!StrStartsWithUnsafe(num_iter, "0,Type=Flag,"))) {
            goto ParseInfoHeader_ret_MALFORMED_INFO_HEADER_LINE;
          }
          new_entry->num = 0;
        } else if (likely(first_num_char == '.')) {
          new_entry->num = kInfoVtypeUnknown;
        } else {
          goto ParseInfoHeader_ret_MALFORMED_INFO_HEADER_LINE;
        }
      } else if (first_num_char > '9') {
        if (first_num_char == 'A') {
          new_entry->num = kInfoVtypeA;
        } else if (first_num_char == 'R') {
          new_entry->num = kInfoVtypeR;
        } else if (likely(first_num_char == 'G')) {
          new_entry->num = kInfoVtypeG;
        } else {
          goto ParseInfoHeader_ret_MALFORMED_INFO_HEADER_LINE;
        }
      } else {
        uint32_t val;
        if (unlikely(ScanmovPosintCapped(UINT32_MAX, &num_iter, &val) || (num_iter[0] != ','))) {
          goto ParseInfoHeader_ret_MALFORMED_INFO_HEADER_LINE;
        }
        new_entry->num = val;
      }
    }
    const uintptr_t info_key_ct = info_keys_iter - info_keys;
#ifdef __LP64__
    if (unlikely(info_key_ct > 0x7ffffffdU)) {
      logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 INFO keys.\n");
      reterr = kPglRetMalformedInput;
      goto ParseInfoHeader_ret_1;
    }
#endif
    assert(info_key_ct);
    *info_key_ctp = info_key_ct;
    BigstackBaseSet(info_keys_iter);
    BigstackEndSet(tmp_alloc_end);
    bigstack_end_mark = g_bigstack_end;
    const uintptr_t info_key_ctl = BitCtToWordCt(info_key_ct);
    uintptr_t* dummy_include;
    if (unlikely(
            (g_bigstack_base > g_bigstack_end) ||
            bigstack_end_alloc_w(info_key_ctl, &dummy_include))) {
      goto ParseInfoHeader_ret_NOMEM;
    }
    SetAllBits(info_key_ct, dummy_include);
    reterr = AllocAndPopulateIdHtableMt(dummy_include, info_keys, info_key_ct, bigstack_left() / 32, 1, info_keys_htablep, nullptr, info_keys_htable_sizep, nullptr);
    if (unlikely(reterr)) {
      goto ParseInfoHeader_ret_1;
    }
    *info_keys_ptr = info_keys;
    bigstack_mark = g_bigstack_base;
  }
  while (0) {
  ParseInfoHeader_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ParseInfoHeader_ret_MALFORMED_INFO_HEADER_LINE:
    logputs("\n");
    logerrputs("Error: Malformed or unrecognized INFO header line.\n");
    reterr = kPglRetMalformedInput;
    break;
  }
 ParseInfoHeader_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

PglErr WritePvarSplit(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, const char* varid_template_str, const char* missing_varid_match, const char* const* info_keys, const uint32_t* info_keys_htable, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t new_variant_id_max_allele_slen, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, UnsortedVar vpos_sortstatus, uint32_t info_key_ct, uint32_t info_keys_htable_size, MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uint32_t thread_ct, char* xheader) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  TextStream pvar_reload_txs;
  PreinitCstream(&css);
  PreinitTextStream(&pvar_reload_txs);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WritePvarSplit_ret_NOMEM;
    }
    const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
    const uint32_t varid_dup = (make_plink2_flags / kfMakePlink2VaridDup) & 1;
    VaridTemplate* varid_templatep = nullptr;
    if (!missing_varid_match) {
      missing_varid_match = &(g_one_char_strs[92]); // '.'
    }
    uint32_t missing_varid_slen = strlen(missing_varid_match);
    uint32_t missing_varid_match_blen = 0; // nonzero iff --set-missing-var-ids
    if (varid_template_str) {
      if (misc_flags & kfMiscSetMissingVarIds) {
        missing_varid_match_blen = missing_varid_slen + 1;
      }
      if (unlikely(BIGSTACK_ALLOC_X(VaridTemplate, 1, &varid_templatep))) {
        goto WritePvarSplit_ret_NOMEM;
      }
      const uint32_t overflow_substitute_blen = new_variant_id_overflow_missing? (missing_varid_slen + 1) : 0;
      VaridTemplateInit(varid_template_str, missing_varid_match, chr_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, varid_templatep);
      if (varid_dup) {
        for (uint32_t uii = 0; uii != varid_templatep->insert_ct; ++uii) {
          const uint32_t insert_type = varid_templatep->insert_types[uii];
          if ((insert_type == 3) || ((insert_type == 2) && (varid_templatep->alleles_needed & 4))) {
            // Could define what takes precedence here, but simpler to prohibit
            // this combination.
            logerrputs("Error: 'vid-[split-]dup' cannot be used with a --set-all-var-ids or\n--set-missing-var-ids template string containing a non-REF allele.\n");
            goto WritePvarSplit_ret_INVALID_CMDLINE;
          }
        }
      }
    }

    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_flags / kfPvarZs) & 1;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WritePvarSplit_ret_1;
    }

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_flags & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      write_info_pr = !IntersectionIsEmpty(variant_include, nonref_flags, raw_variant_ctl);
    }
    write_info_pr = write_info_pr && write_info;
    if (unlikely(write_info_pr && (info_flags & kfInfoPrNonflagPresent))) {
      logputs("\n");
      logerrputs("Error: Conflicting INFO:PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO:PR key.\n");
      goto WritePvarSplit_ret_INCONSISTENT_INPUT;
    }

    char* pvar_info_line_iter = nullptr;
    uint32_t write_filter = 0;
    if (pvar_psam_flags & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybefilter) && filter_present) {
      write_filter = !IntersectionIsEmpty(variant_include, filter_present, raw_variant_ctl);
    }
    uint32_t info_col_idx = 0;  // could save this during first load instead
    const uint32_t info_pr_flag_present = (info_flags / kfInfoPrFlagPresent) & 1;
    if (pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader)) {
      reterr = PvarXheaderWrite(variant_include, cip, variant_bps, allele_idx_offsets, allele_storage, nullptr, xheader_blen, (pvar_psam_flags / kfPvarColVcfheader) & 1, write_filter, write_info, write_info_pr && (!info_pr_flag_present), max_allele_slen, vpos_sortstatus, xheader, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WritePvarSplit_ret_1;
      }
    }
    // could also make this an array-of-structs
    uint32_t* info_key_order = nullptr;
    const char** info_starts = nullptr;
    const char** info_ends = nullptr;
    const char** info_curs = nullptr;
    uint32_t* info_ref_blens = nullptr;
    if (pvar_info_reload) {
      if (unlikely(
              bigstack_alloc_u32(info_key_ct, &info_key_order) ||
              bigstack_alloc_kcp(info_key_ct, &info_starts) ||
              bigstack_alloc_kcp(info_key_ct, &info_ends) ||
              bigstack_alloc_kcp(info_key_ct, &info_curs) ||
              bigstack_alloc_u32(info_key_ct, &info_ref_blens))) {
        goto WritePvarSplit_ret_NOMEM;
      }
      reterr = PvarInfoOpenAndReloadHeader(pvar_info_reload, 1 + (thread_ct > 1), &pvar_reload_txs, &pvar_info_line_iter, &info_col_idx);
      if (unlikely(reterr)) {
        goto WritePvarSplit_ret_TSTREAM_FAIL;
      }
    }
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &cswritep);
    }
    cswritep = strcpya_k(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_flags & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybequal) && qual_present) {
      write_qual = !IntersectionIsEmpty(variant_include, qual_present, raw_variant_ctl);
    }
    if (write_qual) {
      cswritep = strcpya_k(cswritep, "\tQUAL");
    }
    if (write_filter) {
      cswritep = strcpya_k(cswritep, "\tFILTER");
    }
    if (write_info) {
      cswritep = strcpya_k(cswritep, "\tINFO");
    }

    uint32_t write_cm = 0;
    if (pvar_psam_flags & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
        // nonzero_cm_present check was performed
        write_cm = 1;
      } else {
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          if (variant_cms[variant_uidx] != 0.0) {
            write_cm = 1;
            break;
          }
        }
      }
    }
    if (write_cm) {
      cswritep = strcpya_k(cswritep, "\tCM");
    }
    AppendBinaryEoln(&cswritep);

    const VaridTemplate* cur_varid_templatep = nullptr;
    const char* varid_token_start = nullptr; // for vid-split
    const uint32_t varid_split = (make_plink2_flags / kfMakePlink2VaridSemicolon) & 1;
    const uint32_t varid_dup_nosplit = varid_dup && (!varid_split);
    const uint32_t split_just_snps = ((make_plink2_flags & (kfMakePlink2MSplitBase * 3)) == kfMakePlink2MSplitSnps);
    uint32_t trs_variant_uidx = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t orig_allele_ct = 2;
    uint32_t cur_info_key_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        char* chr_name_end = chrtoa(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end = '\t';
        const uint32_t chr_slen = chr_name_end - chr_buf;
        chr_buf_blen = 1 + chr_slen;
        if (varid_templatep) {
          const int32_t chr_slen_delta = chr_slen - varid_templatep->chr_slen;
          varid_templatep->chr_slen = chr_slen;
          varid_templatep->base_len += chr_slen_delta;
        }
      }
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        orig_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      const char* orig_variant_id = variant_ids[variant_uidx];
      const char* ref_allele = cur_alleles[0];
      const uint32_t ref_allele_slen = strlen(ref_allele);
      const uint32_t cur_bp = variant_bps[variant_uidx];
      uint32_t split_ct_p1 = orig_allele_ct;
      uint32_t keep_orig_id = 1;
      if (orig_allele_ct > 2) {
        if (!varid_dup_nosplit) {
          keep_orig_id = 0;
          if (varid_templatep && (!missing_varid_match_blen)) {
            cur_varid_templatep = varid_templatep;
          } else {
            cur_varid_templatep = nullptr;
            if (varid_split) {
              if (VaridSplitOk(orig_variant_id, orig_allele_ct)) {
                varid_token_start = orig_variant_id;
              } else if (varid_dup) {
                keep_orig_id = 1;
              } else {
                varid_token_start = nullptr;
              }
            }
            if ((!varid_token_start) && varid_templatep) {
              // Note that --set-missing-var-ids almost always applies here
              // when it's specified; only exception is when vid-split was also
              // specified and the split succeeded.
              cur_varid_templatep = varid_templatep;
            }
          }
        }
        // Necessary to distinguish between '-' and '-snps' here.
        if (split_just_snps) {
          for (uint32_t uii = 0; uii != orig_allele_ct; ++uii) {
            if (cur_alleles[uii][1]) {
              split_ct_p1 = 2;
              break;
            }
          }
        }
        if ((split_ct_p1 != 2) && pvar_info_line_iter) {
          reterr = PvarInfoReload(info_col_idx, variant_uidx, &pvar_reload_txs, &pvar_info_line_iter, &trs_variant_uidx);
          if (unlikely(reterr)) {
            goto WritePvarSplit_ret_TSTREAM_FAIL;
          }
          char* info_subtoken_iter = pvar_info_line_iter;
          pvar_info_line_iter = CurTokenEnd(pvar_info_line_iter);
          cur_info_key_ct = 0;
          // special case: if entire info field is '.', treat as zero keys
          if ((info_subtoken_iter[0] != '.') || (pvar_info_line_iter != &(info_subtoken_iter[1]))) {
            while (1) {
              if (unlikely(cur_info_key_ct == info_key_ct)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Too many INFO keys for variant ID '%s'.\n", orig_variant_id);
                goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
              }
              char* info_subtoken_end = AdvToDelimOrEnd(info_subtoken_iter, pvar_info_line_iter, ';');
              char* key_end = AdvToDelimOrEnd(info_subtoken_iter, info_subtoken_end, '=');
              const uint32_t key_slen = key_end - info_subtoken_iter;
              const uint32_t kidx = IdHtableFindNnt(info_subtoken_iter, info_keys, info_keys_htable, key_slen, info_keys_htable_size);
              if (unlikely(kidx == UINT32_MAX)) {
                snprintf(g_logbuf, kLogbufSize, "Error: INFO key for variant ID '%s' missing from header.\n", orig_variant_id);
                goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
              }
              info_key_order[cur_info_key_ct] = kidx;
              const int32_t knum = const_container_of(info_keys[kidx], InfoVtype, key)->num;
              if (key_end == info_subtoken_end) {
                if (unlikely(knum)) {
                  snprintf(g_logbuf, kLogbufSize, "Error: INFO key '%s' for variant ID '%s' does not have an accompanying value.\n", info_keys[kidx], orig_variant_id);
                  goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
                }
              } else {
                if (unlikely(!knum)) {
                  snprintf(g_logbuf, kLogbufSize, "Error: INFO key '%s' for variant ID '%s' has an accompanying value, despite being of type Flag.\n", info_keys[kidx], orig_variant_id);
                  goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
                }
                info_subtoken_iter = &(key_end[1]);

                // don't actually need this for Number=A case
                info_starts[cur_info_key_ct] = info_subtoken_iter;

                info_ends[cur_info_key_ct] = info_subtoken_end;
                if (knum <= kInfoVtypeA) {
                  // (Don't need to do anything else for kInfoVtypeUnknown or
                  // positive; we unconditionally copy all the text in those
                  // cases.)
                  if (knum == kInfoVtypeA) {
                    info_curs[cur_info_key_ct] = info_subtoken_iter;
                  } else {
                    char* ref_value_end = S_CAST(char*, memchr(info_subtoken_iter, ',', info_subtoken_end - info_subtoken_iter));
                    if (unlikely(!ref_value_end)) {
                      snprintf(g_logbuf, kLogbufSize, "Error: Too few values for INFO key '%s', variant ID '%s'.\n", info_keys[kidx], orig_variant_id);
                      goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
                    }
                    ++ref_value_end;
                    info_ref_blens[cur_info_key_ct] = ref_value_end - info_subtoken_iter;
                    info_curs[cur_info_key_ct] = ref_value_end;
                  }
                }
              }
              ++cur_info_key_ct;
              if (info_subtoken_end == pvar_info_line_iter) {
                break;
              }
              info_subtoken_iter = &(info_subtoken_end[1]);
            }
          }
        }
      }
      for (uint32_t alt_allele_idx = 1; alt_allele_idx != split_ct_p1; ++alt_allele_idx) {
        cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
        cswritep = u32toa_x(cur_bp, '\t', cswritep);
        const char* cur_alt_allele = cur_alleles[alt_allele_idx];
        const uint32_t cur_alt_allele_slen = strlen(cur_alt_allele);
        if ((split_ct_p1 == 2) || keep_orig_id) {
          cswritep = strcpyax(cswritep, orig_variant_id, '\t');
          cswritep = memcpyax(cswritep, ref_allele, ref_allele_slen, '\t');
          cswritep = memcpya(cswritep, cur_alt_allele, cur_alt_allele_slen);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WritePvarSplit_ret_WRITE_FAIL;
          }
          if ((orig_allele_ct > 2) && (split_ct_p1 == 2)) {
            // -snps non-split case
            for (uint32_t allele_idx = 2; allele_idx != orig_allele_ct; ++allele_idx) {
              *cswritep++ = ',';
              cswritep = strcpya(cswritep, cur_alleles[allele_idx]);
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto WritePvarSplit_ret_WRITE_FAIL;
              }
            }
          }
        } else {
          if (cur_varid_templatep) {
            // Always true in --set-all-var-ids case.  True in
            // --set-missing-var-ids case when vid-split unspecified, or split
            // failed.
            cswritep = VaridTemplateWrite(cur_varid_templatep, ref_allele, cur_alt_allele, cur_bp, ref_allele_slen, 0, cur_alt_allele_slen, cswritep);
            *cswritep++ = '\t';
          } else if (varid_token_start) {
            const char* varid_token_end = strchrnul(varid_token_start, ';');
            // If substring matches missing code and --set-missing-var-ids is
            // specified, we replace it.
            if (varid_templatep && (S_CAST(uintptr_t, varid_token_end - varid_token_start) == missing_varid_slen) && memequal(varid_token_start, missing_varid_match, missing_varid_slen)) {
              cswritep = VaridTemplateWrite(varid_templatep, ref_allele, cur_alt_allele, cur_bp, ref_allele_slen, 0, cur_alt_allele_slen, cswritep);
            } else {
              cswritep = memcpya(cswritep, varid_token_start, varid_token_end - varid_token_start);
            }
            *cswritep++ = '\t';
            varid_token_start = &(varid_token_end[1]);
          } else {
            cswritep = memcpyax(cswritep, missing_varid_match, missing_varid_slen, '\t');
          }
          cswritep = memcpyax(cswritep, ref_allele, ref_allele_slen, '\t');
          cswritep = memcpya(cswritep, cur_alt_allele, cur_alt_allele_slen);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WritePvarSplit_ret_WRITE_FAIL;
          }
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
            cswritep = strcpya_k(cswritep, "PASS");
          } else {
            cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
          }
        }

        if (write_info) {
          *cswritep++ = '\t';
          const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
          if (pvar_info_line_iter) {
            if (split_ct_p1 == 2) {
              reterr = PvarInfoReloadAndWrite(info_pr_flag_present, info_col_idx, variant_uidx, is_pr, &pvar_reload_txs, &pvar_info_line_iter, &cswritep, &trs_variant_uidx);
              if (unlikely(reterr)) {
                goto WritePvarSplit_ret_TSTREAM_FAIL;
              }
            } else {
              if (!cur_info_key_ct) {
                *cswritep++ = '.';
              } else {
                const uint32_t is_last_allele = (alt_allele_idx + 1 == split_ct_p1);
                for (uint32_t kpos = 0; kpos != cur_info_key_ct; ++kpos) {
                  const uint32_t kidx = info_key_order[kpos];
                  const char* cur_key_str = info_keys[kidx];
                  cswritep = strcpya(cswritep, cur_key_str);
                  const int32_t knum = const_container_of(info_keys[kidx], InfoVtype, key)->num;
                  if (knum) {
                    *cswritep++ = '=';
                    const char* cur_info_start = info_starts[kpos];
                    const char* cur_info_end = info_ends[kpos];
                    if (knum >= kInfoVtypeUnknown) {
                      cswritep = memcpya(cswritep, cur_info_start, cur_info_end - cur_info_start);
                    } else {
                      if (knum != kInfoVtypeA) {
                        cswritep = memcpya(cswritep, cur_info_start, info_ref_blens[kpos]);
                      }
                      // okay, this needs a better name
                      const char* cur_info_cur = info_curs[kpos];

                      const char* subtoken_end = AdvToDelimOrEnd(cur_info_cur, cur_info_end, ',');
                      if (knum == kInfoVtypeG) {
                        if (unlikely(subtoken_end == cur_info_end)) {
                          snprintf(g_logbuf, kLogbufSize, "Error: Too few values for INFO key '%s', variant ID '%s'.\n", cur_key_str, orig_variant_id);
                          goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
                        }
                        cswritep = memcpya(cswritep, cur_info_cur, 1 + S_CAST(uintptr_t, subtoken_end - cur_info_cur));
                        cur_info_cur = subtoken_end;
                        const uint32_t skip_ct = alt_allele_idx - 1;
                        if (skip_ct) {
                          cur_info_cur = AdvToNthDelimChecked(&(cur_info_cur[1]), cur_info_end, skip_ct, ',');
                          if (unlikely(!subtoken_end)) {
                            snprintf(g_logbuf, kLogbufSize, "Error: Too few values for INFO key '%s', variant ID '%s'.\n", cur_key_str, orig_variant_id);
                            goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
                          }
                        }
                        ++cur_info_cur;
                        subtoken_end = AdvToDelimOrEnd(cur_info_cur, cur_info_end, ',');
                      }
                      if (unlikely((subtoken_end == cur_info_end) != is_last_allele)) {
                        snprintf(g_logbuf, kLogbufSize, "Error: Wrong number of values for INFO key '%s', variant ID '%s'.\n", cur_key_str, orig_variant_id);
                        goto WritePvarSplit_ret_MALFORMED_INPUT_WW;
                      }
                      cswritep = memcpya(cswritep, cur_info_cur, subtoken_end - cur_info_cur);
                      info_curs[kpos] = &(subtoken_end[1]);
                    }
                  }
                  *cswritep++ = ';';
                }
                --cswritep;
              }
            }
          } else {
            if (is_pr) {
              cswritep = strcpya_k(cswritep, "PR");
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
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WritePvarSplit_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  WritePvarSplit_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WritePvarSplit_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_reload_txs);
    break;
  WritePvarSplit_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WritePvarSplit_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  WritePvarSplit_ret_MALFORMED_INPUT_WW:
    logputs("\n");
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  WritePvarSplit_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WritePvarSplit_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// Final filter_keys is natural-sorted.
// Return values are allocated on bottom of bigstack.
// Caller must initialize all return values to correspond to the null table.
PglErr MakeFilterHtable(const uintptr_t* variant_include, const uintptr_t* filter_npass, const char* const* filter_storage, uint32_t variant_ct, const char*** filter_keys_ptr, uint32_t** filter_keys_htable_ptr, uint32_t* filter_key_ct_ptr, uint32_t* filter_keys_htable_size_ptr) {
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    // Start with empty size-128 table, which will practically always be enough
    // while still being small relative to L1 cache.  Double table size
    // whenever load factor reaches 0.25; there shouldn't be *that* many
    // distinct filters.
    // possible todo: multithread this scan, merge results at the end; can also
    // separate this stage from the rest of the function.
    uint32_t table_size = 128;
    uint32_t hash_shift = 25; // 32 - log2(table_size)
    uint32_t filter_key_ct = 0;
    char** filter_tokens;
    if (unlikely(
            bigstack_end_calloc_cp(table_size, &filter_tokens))) {
      goto MakeFilterHtable_ret_NOMEM;
    }

    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    uintptr_t variant_widx = 0;
    uintptr_t cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t lowbit = BitIter1y(variant_include, &variant_widx, &cur_bits);
      if (lowbit & filter_npass[variant_widx]) {
        const char* filter_iter = filter_storage[variant_widx * kBitsPerWord + ctzw(lowbit)];
        while (1) {
          const char* token_end = strchrnul(filter_iter, ';');
          const uint32_t cur_id_slen = token_end - filter_iter;
          for (uint32_t hashval = Hash32(filter_iter, cur_id_slen) >> hash_shift; ; ) {
            char* cur_token_ptr = filter_tokens[hashval];
            if (!cur_token_ptr) {
              char* storage_loc;
              if (StoreStringAtBase(tmp_alloc_end, filter_iter, cur_id_slen, &tmp_alloc_base, &storage_loc)) {
                goto MakeFilterHtable_ret_NOMEM;
              }
              ++filter_key_ct;
              if (filter_key_ct * 4 < table_size) {
                filter_tokens[hashval] = storage_loc;
                break;
              }
#ifdef __LP64__
              if (unlikely(hash_shift == 1)) {
                // this is technically "not yet supported", but I fail to see a
                // valid use case for >536 million distinct FILTER keys...
                logerrprintf("Error: Too many distinct FILTER keys (max 2^29 - 1).\n");
                goto MakeFilterHtable_ret_MALFORMED_INPUT;
              }
#endif
              // It's fine for the new table to overlap the old table, since we
              // can iterate through all the strings by walking forward from
              // g_bigstack_base.
              const uintptr_t extra_byte_ct = table_size * sizeof(intptr_t);
              if (unlikely(S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) < extra_byte_ct)) {
                goto MakeFilterHtable_ret_NOMEM;
              }
              tmp_alloc_end -= extra_byte_ct;
              filter_tokens = R_CAST(char**, tmp_alloc_end);
              memset(filter_tokens, 0, 2 * extra_byte_ct);
              table_size *= 2;
              --hash_shift;
              char* rehash_iter = R_CAST(char*, g_bigstack_base);
              for (uint32_t uii = 0; uii != filter_key_ct; ++uii) {
                char* rehash_token_end = strnul(rehash_iter);
                const uint32_t rehash_id_slen = rehash_token_end - rehash_iter;
                for (uint32_t rehashval = Hash32(rehash_iter, rehash_id_slen) >> hash_shift; ; ) {
                  if (!filter_tokens[rehashval]) {
                    filter_tokens[rehashval] = rehash_iter;
                    break;
                  }
                  if (++rehashval == table_size) {
                    rehashval = 0;
                  }
                }
                rehash_iter = &(rehash_token_end[1]);
              }
              break;
            }
            if ((!memcmp(filter_iter, cur_token_ptr, cur_id_slen)) && (!cur_token_ptr[cur_id_slen])) {
              break;
            }
            if (++hashval == table_size) {
              hashval = 0;
            }
          }
          if (!(*token_end)) {
            break;
          }
          filter_iter = &(token_end[1]);
        }
      }
    }
    if (!filter_key_ct) {
      // All nonpassing variants were already filtered out.
      // Caller already initialized null table.
      goto MakeFilterHtable_ret_1;
    }
    char* token_iter = R_CAST(char*, g_bigstack_base);
    BigstackBaseSet(tmp_alloc_base);
    const uint32_t filter_keys_htable_size = GetHtableFastSize(filter_key_ct);
    if (unlikely(
            bigstack_alloc_kcp(filter_key_ct, filter_keys_ptr) ||
            bigstack_alloc_u32(filter_keys_htable_size, filter_keys_htable_ptr))) {
      goto MakeFilterHtable_ret_NOMEM;
    }
    const char** filter_keys = *filter_keys_ptr;
    for (uint32_t uii = 0; uii != filter_key_ct; ++uii) {
      filter_keys[uii] = token_iter;
      char* token_end = strnul(token_iter);
      token_iter = &(token_end[1]);
    }
    StrptrArrNsort(filter_key_ct, filter_keys);
    *filter_key_ct_ptr = filter_key_ct;
    *filter_keys_htable_size_ptr = filter_keys_htable_size;
    uint32_t* filter_keys_htable = *filter_keys_htable_ptr;
    SetAllU32Arr(filter_keys_htable_size, filter_keys_htable);
    for (uint32_t uii = 0; uii != filter_key_ct; ++uii) {
      for (uint32_t hashval = Hashceil(filter_keys[uii], strlen(filter_keys[uii]), filter_keys_htable_size); ; ) {
        if (filter_keys_htable[hashval] == UINT32_MAX) {
          filter_keys_htable[hashval] = uii;
        }
        if (++hashval == filter_keys_htable_size) {
          hashval = 0;
        }
      }
    }
  }
  while (0) {
  MakeFilterHtable_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
#ifdef __LP64__
  MakeFilterHtable_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
#endif
  }
 MakeFilterHtable_ret_1:
  BigstackEndReset(bigstack_end_mark);
  return reterr;
}

/*
PglErr WritePvarJoin(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, const char* varid_template_str, const char* missing_varid_match, const char* const* info_keys, const uint32_t* info_keys_htable, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t new_variant_id_max_allele_slen, uint32_t max_write_allele_ct, uint32_t max_missalt_ct, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, UnsortedVar vpos_sortstatus, uint32_t info_key_ct, uint32_t info_keys_htable_size, MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uint32_t thread_ct, char* xheader) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  TextStream pvar_reload_txs;
  PreinitCstream(&css);
  PreinitTextStream(&pvar_reload_txs);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WritePvarJoin_ret_NOMEM;
    }
    const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
    const uint32_t varid_dup = (make_plink2_flags / kfMakePlink2VaridDup) & 1;
    VaridTemplate* varid_templatep = nullptr;
    if (!missing_varid_match) {
      missing_varid_match = &(g_one_char_strs[92]); // '.'
    }
    uint32_t missing_varid_slen = strlen(missing_varid_match);
    uint32_t missing_varid_match_blen = 0; // nonzero iff --set-missing-var-ids
    if (varid_template_str) {
      if (misc_flags & kfMiscSetMissingVarIds) {
        missing_varid_match_blen = missing_varid_slen + 1;
      }
      if (unlikely(BIGSTACK_ALLOC_X(VaridTemplate, 1, &varid_templatep))) {
        goto WritePvarJoin_ret_NOMEM;
      }
      const uint32_t overflow_substitute_blen = new_variant_id_overflow_missing? (missing_varid_slen + 1) : 0;
      VaridTemplateInit(varid_template_str, missing_varid_match, chr_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, varid_templatep);
      if (varid_dup) {
        for (uint32_t uii = 0; uii != varid_templatep->insert_ct; ++uii) {
          const uint32_t insert_type = varid_templatep->insert_types[uii];
          if ((insert_type == 3) || ((insert_type == 2) && (varid_templatep->alleles_needed & 4))) {
            // Could define what takes precedence here, but simpler to prohibit
            // this combination.
            logerrputs("Error: 'vid-[split-]dup' cannot be used with a --set-all-var-ids or\n--set-missing-var-ids template string containing a non-REF allele.\n");
            goto WritePvarJoin_ret_INVALID_CMDLINE;
          }
        }
      }
    }

    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + S_CAST(uintptr_t, info_reload_slen) * (max_write_allele_ct - 1);
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_flags / kfPvarZs) & 1;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WritePvarJoin_ret_1;
    }

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_flags & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      write_info_pr = !IntersectionIsEmpty(variant_include, nonref_flags, raw_variant_ctl);
    }
    write_info_pr = write_info_pr && write_info;
    if (unlikely(write_info_pr && (info_flags & kfInfoPrNonflagPresent))) {
      logputs("\n");
      logerrputs("Error: Conflicting INFO:PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO:PR key.\n");
      goto WritePvarJoin_ret_INCONSISTENT_INPUT;
    }

    char* pvar_info_line_iter = nullptr;
    uint32_t write_filter = 0;
    if (pvar_psam_flags & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybefilter) && filter_present) {
      write_filter = !IntersectionIsEmpty(variant_include, filter_present, raw_variant_ctl);
    }
    uint32_t info_col_idx = 0;  // could save this during first load instead
    const uint32_t info_pr_flag_present = (info_flags / kfInfoPrFlagPresent) & 1;
    if (pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader)) {
      reterr = PvarXheaderWrite(variant_include, cip, variant_bps, allele_idx_offsets, allele_storage, nullptr, xheader_blen, (pvar_psam_flags / kfPvarColVcfheader) & 1, write_filter, write_info, write_info_pr && (!info_pr_flag_present), max_allele_slen, vpos_sortstatus, xheader, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WritePvarJoin_ret_1;
      }
    }
    const uint32_t join_mode = (make_plink2_flags & (kfMakePlink2MSplitBase * 7));
    uintptr_t info_cache_size = max_missalt_ct + max_write_allele_ct - 1;
    if (join_mode != kfMakePlink2MJoinSnps) {
      info_cache_size *= 3;
    }
#ifndef __LP64__
    if (S_CAST(uint64_t, info_cache_size) * info_key_ct * sizeof(intptr_t) > 0x7fffffff) {
      goto WritePvarJoin_ret_NOMEM;
    }
#endif

    if (cip->chrset_source) {
      AppendChrsetLine(cip, &cswritep);
    }
    cswritep = strcpya_k(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_flags & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybequal) && qual_present) {
      write_qual = !IntersectionIsEmpty(variant_include, qual_present, raw_variant_ctl);
    }
    if (write_qual) {
      cswritep = strcpya_k(cswritep, "\tQUAL");
    }
    const char** filter_keys = nullptr;
    uint32_t* filter_keys_htable = nullptr;
    uintptr_t* cur_filter_keys = nullptr;
    uint32_t filter_keys_htable_size = 0;
    uint32_t filter_key_ct = 0;
    uint32_t filter_key_ctl = 0;
    if (write_filter) {
      // The VCF spec doesn't require ##FILTER= header lines, and unlike the
      // case with INFO Number=A/R/G, we can join correctly without header
      // information.  It's slightly computationally more expensive, but INFO
      // and genotype joining costs are more significant.
      if (filter_npass) {
        reterr = MakeFilterHtable(variant_include, filter_npass, filter_storage, variant_ct, &filter_keys, &filter_keys_htable, &filter_key_ct, &filter_keys_htable_size);
        if (unlikely(reterr)) {
          goto WritePvarJoin_ret_1;
        }
        if (filter_key_ct) {
          filter_key_ctl = BitCtToWordCt(filter_key_ct);
          if (unlikely(bigstack_alloc_w(filter_key_ctl, &cur_filter_keys))) {
            goto WritePvarJoin_ret_1;
          }
        }
      }
      cswritep = strcpya_k(cswritep, "\tFILTER");
    }

    char** info_bufs = nullptr;
    const char** info_starts = nullptr;
    const char** info_ends = nullptr;  // ugh, this is not related to INFO:END
    const char** info_curs = nullptr;
    uint32_t info_end_key_idx = UINT32_MAX;
    if (pvar_info_reload) {
      if (unlikely(
              bigstack_alloc_cp(info_cache_size, &info_bufs) ||
              bigstack_alloc_kcp(info_key_ct * info_cache_size, &info_starts) ||
              bigstack_alloc_kcp(info_key_ct * info_cache_size, &info_ends) ||
              bigstack_alloc_kcp(info_key_ct * info_cache_size, &info_curs))) {
        goto WritePvarJoin_ret_NOMEM;
      }
      reterr = PvarInfoOpenAndReloadHeader(pvar_info_reload, 1 + (thread_ct > 1), &pvar_reload_txs, &pvar_info_line_iter, &info_col_idx);
      if (unlikely(reterr)) {
        goto WritePvarJoin_ret_TSTREAM_FAIL;
      }
      info_end_key_idx = IdHtableFind("END", info_keys, info_keys_htable, strlen("END"), info_keys_htable_size);
      if (info_end_key_idx != UINT32_MAX) {
        const int32_t knum = const_container_of(info_keys[info_end_key_idx], InfoVtype, key)->num;
        if ((knum != 1) && (knum != kInfoVtypeUnknown)) {
          // TODO: verify type instead.
          // but if number is not . or 1, this is not the INFO:END we're
          // looking for.
          info_end_key_idx = UINT32_MAX;
        }
      }
    }
    if (write_info) {
      cswritep = strcpya_k(cswritep, "\tINFO");
    }

    uint32_t write_cm = 0;
    if (pvar_psam_flags & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
        // nonzero_cm_present check was performed
        write_cm = 1;
      } else {
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          if (variant_cms[variant_uidx] != 0.0) {
            write_cm = 1;
            break;
          }
        }
      }
    }
    if (write_cm) {
      cswritep = strcpya_k(cswritep, "\tCM");
    }
    AppendBinaryEoln(&cswritep);

    const VaridTemplate* cur_varid_templatep = nullptr;
    const char* varid_token_start = nullptr; // for vid-split
    const uint32_t varid_split = (make_plink2_flags / kfMakePlink2VaridSemicolon) & 1;
    const uint32_t varid_dup_nosplit = varid_dup && (!varid_split);
    uint32_t next_variant_idx = 0;
    uint32_t trs_variant_uidx = 0;
    uint32_t next_variant_uidx = 0;
    uintptr_t next_variant_uidx_base = 0;
    uintptr_t next_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t prev_bp = 0;
    uint32_t cur_bp = 0;
    uint32_t bp_start_variant_idx = 0;
    uint32_t bp_start_variant_uidx = 0;
    uintptr_t bp_start_variant_uidx_base = 0;
    uintptr_t bp_start_bits = variant_include[0];
    uint32_t allele_ct = 2;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    JoinCounts jc;
    jc.snp_ct = 0;
    jc.nonsnp_ct = 0;
    jc.symbolic_ct = 0;
    jc.missalt_snp_ct = 0;
    jc.missalt_nonsnp_ct = 0;
    JoinCounts next_jc = jc;
    fputs("0%", stdout);
    fflush(stdout);
    while (1) {
      for (; next_variant_idx != variant_ct; ++next_variant_idx) {
        next_variant_uidx = BitIter1(variant_include, &next_variant_uidx_base, &next_bits);
        if (next_variant_uidx >= chr_end) {
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          } while (next_variant_uidx >= chr_end);
          char* chr_name_end = chrtoa(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
          *chr_name_end = '\t';
          const uint32_t chr_slen = chr_name_end - chr_buf;
          chr_buf_blen = 1 + chr_slen;
          if (varid_templatep) {
            const int32_t chr_slen_delta = chr_slen - varid_templatep->chr_slen;
            varid_templatep->chr_slen = chr_slen;
            varid_templatep->base_len += chr_slen_delta;
          }
          prev_bp = UINT32_MAX;
        }
        cur_bp = variant_bps[next_variant_uidx];
        if (cur_bp != prev_bp) {
          break;
        }
        uintptr_t allele_idx_offset_base;
        if (!allele_idx_offsets) {
          allele_idx_offset_base = next_variant_uidx * 2;
        } else {
          allele_idx_offset_base = allele_idx_offsets[next_variant_uidx];
          allele_ct = allele_idx_offsets[next_variant_uidx + 1] - allele_idx_offset_base;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        JoinVtype jvt = JoinCount(cur_alleles, allele_ct, &next_jc);
        // previously validated
        // if ((join_mode == kfMakePlink2MJoinSnps) && ()) {
        // }

        // TODO
        jc.snp_ct += next_jc.snp_ct;
        jc.nonsnp_ct += next_jc.nonsnp_ct;
        jc.symbolic_ct += next_jc.symbolic_ct;
        jc.missalt_snp_ct += next_jc.missalt_snp_ct;
        jc.missalt_nonsnp_ct += next_jc.missalt_nonsnp_ct;
      }
      if (next_variant_idx == bp_start_variant_idx + 1) {
        // No join needed.  This is usually the common case, so we duplicate a
        // bunch of code for the sake of avoiding slowdown here.
        cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
        cswritep = u32toa_x(variant_bps[bp_start_variant_uidx], '\t', cswritep);
        cswritep = strcpyax(cswritep, variant_ids[bp_start_variant_uidx], '\t');
        uintptr_t allele_idx_offset_base;
        if (!allele_idx_offsets) {
          allele_idx_offset_base = bp_start_variant_uidx * 2;
        } else {
          allele_idx_offset_base = allele_idx_offsets[bp_start_variant_uidx];
          allele_ct = allele_idx_offsets[bp_start_variant_uidx + 1] - allele_idx_offset_base;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        cswritep = strcpyax(cswritep, cur_alleles[0], '\t');
        cswritep = strcpya(cswritep, cur_alleles[1]);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto WritePvarJoin_ret_WRITE_FAIL;
        }
        for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
          *cswritep++ = ',';
          cswritep = strcpya(cswritep, cur_alleles[allele_idx]);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WritePvarJoin_ret_WRITE_FAIL;
          }
        }

        if (write_qual) {
          *cswritep++ = '\t';
          if ((!qual_present) || (!IsSet(qual_present, bp_start_variant_uidx))) {
            *cswritep++ = '.';
          } else {
            cswritep = ftoa_g(quals[bp_start_variant_uidx], cswritep);
          }
        }

        if (write_filter) {
          *cswritep++ = '\t';
          if ((!filter_present) || (!IsSet(filter_present, bp_start_variant_uidx))) {
            *cswritep++ = '.';
          } else if (!IsSet(filter_npass, bp_start_variant_uidx)) {
            cswritep = strcpya_k(cswritep, "PASS");
          } else {
            cswritep = strcpya(cswritep, filter_storage[bp_start_variant_uidx]);
          }
        }

        if (write_info) {
          *cswritep++ = '\t';
          const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, bp_start_variant_uidx));
          if (pvar_info_line_iter) {
            reterr = PvarInfoReloadAndWrite(info_pr_flag_present, info_col_idx, bp_start_variant_uidx, is_pr, &pvar_reload_txs, &pvar_info_line_iter, &cswritep, &trs_variant_uidx);
            if (unlikely(reterr)) {
              goto WritePvarJoin_ret_TSTREAM_FAIL;
            }
          } else {
            if (is_pr) {
              cswritep = strcpya_k(cswritep, "PR");
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
            cswritep = dtoa_g_p8(variant_cms[bp_start_variant_uidx], cswritep);
          }
        }
        AppendBinaryEoln(&cswritep);
        // next_jc guaranteed to be zero-initialized
      } else if (next_variant_idx) {
        // TODO
        ;;;;
      const char* orig_variant_id = variant_ids[variant_uidx];
      const char* ref_allele = cur_alleles[0];
      const uint32_t ref_allele_slen = strlen(ref_allele);
      uint32_t split_ct_p1 = allele_ct;
      if (allele_ct > 2) {
        if (!varid_dup) {
          if (varid_templatep && (!missing_varid_match_blen)) {
            cur_varid_templatep = varid_templatep;
          } else {
            cur_varid_templatep = nullptr;
            if (varid_split) {
              if (VaridSplitOk(orig_variant_id, allele_ct)) {
                varid_token_start = orig_variant_id;
              } else {
                varid_token_start = nullptr;
              }
            }
            if ((!varid_token_start) && varid_templatep) {
              // Note that --set-missing-var-ids almost always applies here
              // when it's specified; only exception is when vid-split was also
              // specified and the split succeeded.
              cur_varid_templatep = varid_templatep;
            }
          }
        }
      }
       ;;;;
        next_jc.snp_ct = 0;
        next_jc.nonsnp_ct = 0;
        next_jc.symbolic_ct = 0;
        next_jc.missalt_snp_ct = 0;
        next_jc.missalt_nonsnp_ct = 0;
      }
      if (next_variant_idx == variant_ct) {
        break;
      }
      // this_pos_write_variant_ct = 0;
      jc = next_jc;
      prev_bp = cur_bp;
      bp_start_variant_idx = next_variant_idx;
      bp_start_variant_uidx = next_variant_uidx;
      if (next_variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (next_variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WritePvarJoin_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  WritePvarJoin_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WritePvarJoin_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_reload_txs);
    break;
  WritePvarJoin_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WritePvarJoin_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  WritePvarJoin_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WritePvarJoin_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}
*/

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
  kfPlink2WriteFillMissingFromDosage = (1 << 8),
  // no need for sample_sort, determined by collapsed_sort_map != nullptr?
  kfPlink2WritePlink1 = (1 << 9)
FLAGSET_DEF_END(Plink2WriteFlags);
// todo: add .pgen-specific stuff

typedef struct MakeCommonStruct {
  const ChrInfo* cip;
  const uintptr_t* sample_include;
  uintptr_t* sex_male_collapsed_interleaved;
  uintptr_t* sex_female_collapsed_interleaved;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select);
  uint32_t raw_sample_ct;
  uint32_t sample_ct;
  Plink2WriteFlags plink2_write_flags;
  uint32_t hard_call_halfdist;
} MakeCommon;

typedef struct MakeBedlikeCtxStruct {
  const MakeCommon* mcp;

  const uintptr_t* variant_include;
  uint32_t* sample_include_cumulative_popcounts;
  const uint32_t* collapsed_sort_map;

  PgenReader** pgr_ptrs;

  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_write_ct;

  uintptr_t** genovecs;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;

  unsigned char* writebufs[2];
  PglErr reterr;  // can only be kPglRetMalformedInput for now
} MakeBedlikeCtx;


THREAD_FUNC_DECL MakeBedlikeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  MakeBedlikeCtx* ctx = S_CAST(MakeBedlikeCtx*, arg->sharedp->context);

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  uintptr_t* genovec = ctx->genovecs[tidx];
  const MakeCommon* mcp = ctx->mcp;
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  uint32_t hard_call_halfdist = 0;
  if (ctx->dosage_presents) {
    dosage_present = ctx->dosage_presents[tidx];
    dosage_main = ctx->dosage_mains[tidx];
    hard_call_halfdist = mcp->hard_call_halfdist;
  }
  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = mcp->cip;
  const uintptr_t* sample_include = mcp->sample_include;
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  const uintptr_t* sex_male_collapsed_interleaved = mcp->sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed_interleaved = mcp->sex_female_collapsed_interleaved;
  const uint32_t* collapsed_sort_map = ctx->collapsed_sort_map;
  const Plink2WriteFlags plink2_write_flags = mcp->plink2_write_flags;
  const uint32_t set_hh_missing = plink2_write_flags & kfPlink2WriteSetHhMissing;
  const uint32_t set_mixed_mt_missing = plink2_write_flags & kfPlink2WriteSetMixedMtMissing;
  const uint32_t write_plink1 = plink2_write_flags & kfPlink2WritePlink1;
  const uint32_t sample_ct = mcp->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t sample_ctv2 = NypCtToVecCt(sample_ct);
  const uint32_t sample_ct4 = NypCtToByteCt(sample_ct);
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select) = mcp->refalt1_select;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
  uint32_t parity = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(ctx->writebufs[parity][write_idx * sample_ct4]);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    uint32_t chr_end = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t is_haploid_nonmt = 0;
    uint32_t is_mt = 0;
    for (; write_idx != write_idx_end; ++write_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
        is_mt = (chr_idx == mt_code);
        is_haploid_nonmt = IsSet(cip->haploid_mask, chr_idx) && (!is_mt);
      }
      // todo: Multiallelic -> two-specific-alleles downcode.
      // This is pretty straightforward if we're just saving hardcalls:
      // with 1 copy of one allele and zero copies of the other allele, we
      // default to saving a missing call (in the diploid case).
      // If dosages are involved, things are a bit less obvious: what if the
      // unincluded alleles have a total dosage of 0.1?  0.5?  It'll be
      // necessary to define a new flag allowing this threshold to be
      // configured.
      // I'm currently inclined to set unincluded dosage >= 0.5 to missing, and
      // otherwise the dosages are scaled up to sum to 2.
      // (Note that the multiallelic split operation won't work this way; it
      // has to use the convention that REF = anything other than the current
      // ALT allele.  Probably also want to support that here.)
      if (!hard_call_halfdist) {
        // if multiallelic:
        //   if split: call PgrGet1()
        //   otherwise, if erase-alt2+: call PgrGet2()
        //   otherwise, error out
        PglErr reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, pgrp, genovec);
        if (unlikely(reterr)) {
          ctx->reterr = reterr;
          break;
        }
      } else {
        // this isn't fully implemented yet.

        // quasi-bugfix (4 Dec 2017): it's user-hostile to make
        // --hard-call-threshold not apply here.
        uint32_t dosage_ct;
        // if multiallelic:
        //    if split: call PgrGet1D()
        //    otherwise, if refalt1_select + erase-alt2+: call PgrGetMD(),
        //      rescale
        //    otherwise, error out
        PglErr reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
        if (unlikely(reterr)) {
          ctx->reterr = reterr;
          break;
        }
        ApplyHardCallThresh(dosage_present, dosage_main, dosage_ct, hard_call_halfdist, genovec);
      }
      // remove this when proper multiallelic logic implemented
      if (refalt1_select && (refalt1_select[variant_uidx][0] == 1)) {
        GenovecInvertUnsafe(sample_ct, genovec);
      }
      if (set_hh_missing && is_haploid_nonmt) {
        if (is_x) {
          SetMaleHetMissing(sex_male_collapsed_interleaved, sample_ctv2, genovec);
        } else {
          // all hets to missing
          SetHetMissing(sample_ctl2, genovec);
          if (is_y) {
            InterleavedSetMissing(sex_female_collapsed_interleaved, sample_ctv2, genovec);
          }
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
      ZeroTrailingNyps(sample_ct, genovec);
      if (!collapsed_sort_map) {
        writebuf_iter = memcpyua(writebuf_iter, genovec, sample_ct4);
      } else {
        GenovecResort(genovec, collapsed_sort_map, sample_ct, writebuf_iter);
        writebuf_iter = &(writebuf_iter[sample_ct4]);
      }
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// initialized mcp fields: cip, sex_male_collapsed_interleaved,
// sex_female_collapsed_interleaved, raw_sample_ct, sample_ct,
// plink2_write_flags
PglErr MakeBedlikeMain(const uintptr_t* sample_include, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_thread_ct, uint32_t hard_call_thresh, MakePlink2Flags make_plink2_flags, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, MakeCommon* mcp, char* outname, char* outname_end) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  MakeBedlikeCtx ctx;
  {
    assert(variant_ct);
    const uint32_t sample_ct = mcp->sample_ct;
    assert(sample_ct);
    if (make_plink2_flags & kfMakePlink2MMask) {
      logerrputs("Error: Multiallelic-split fixed-width output is not implemented yet.\n");
      reterr = kPglRetNotYetSupported;
      goto MakeBedlikeMain_ret_1;
    }
    // fixed-width
    const uint32_t make_pgen = make_plink2_flags & kfMakePgen;
    if (make_pgen) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    } else {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".bed");
    }
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto MakeBedlikeMain_ret_OPEN_FAIL;
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
      if (unlikely(ferror_unlocked(outfile))) {
        goto MakeBedlikeMain_ret_WRITE_FAIL;
      }
    } else {
      if (unlikely(fwrite_checked("l\x1b\x01", 3, outfile))) {
        goto MakeBedlikeMain_ret_WRITE_FAIL;
      }
    }
    logprintfww5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;
    const uint32_t raw_sample_ct = mcp->raw_sample_ct;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uintptr_t sample_ct4 = NypCtToByteCt(sample_ct);
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &ctx.sample_include_cumulative_popcounts))) {
      goto MakeBedlikeMain_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, ctx.sample_include_cumulative_popcounts);
    // tried more threads, pointless since this is too I/O-bound
    // (exception: reordering samples)
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    ctx.collapsed_sort_map = new_sample_idx_to_old;
    if (!new_sample_idx_to_old) {
      // Without BMI2 instructions, subsetting is most expensive with
      // sample_ct near 2/3 of raw_sample_ct; up to ~7 compute threads are
      // useful in that case.  (See CopyNyparrNonemptySubset().)
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
      if (unlikely(bigstack_alloc_u32(sample_ct, &new_collapsed_sort_map))) {
        goto MakeBedlikeMain_ret_NOMEM;
      }
      // bugfix (26 Mar 2018): forgot to initialize this
      memcpy(new_collapsed_sort_map, new_sample_idx_to_old, sample_ct * sizeof(int32_t));
      UidxsToIdxs(sample_include, ctx.sample_include_cumulative_popcounts, sample_ct, new_collapsed_sort_map);
      ctx.collapsed_sort_map = new_collapsed_sort_map;
    }

    if (make_plink2_flags & kfMakeBed) {
      mcp->plink2_write_flags |= kfPlink2WritePlink1;
    }

    mcp->hard_call_halfdist = 0;
    if ((hard_call_thresh != UINT32_MAX) && (pgfip->gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))) {
      mcp->hard_call_halfdist = kDosage4th - hard_call_thresh;
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    ctx.dosage_presents = nullptr;
    ctx.dosage_mains = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, 0, 2 * (sample_ct4 + 1), 0, pgfip, &calc_thread_ct, &ctx.genovecs, nullptr, nullptr, nullptr, mcp->hard_call_halfdist? (&ctx.dosage_presents) : nullptr, mcp->hard_call_halfdist? (&ctx.dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto MakeBedlikeMain_ret_NOMEM;
    }
    if (unlikely(
            bigstack_alloc_uc(sample_ct4 * read_block_size, &(ctx.writebufs[0])) ||
            bigstack_alloc_uc(sample_ct4 * read_block_size, &(ctx.writebufs[1])))) {
      // shouldn't be possible for this to fail
      goto MakeBedlikeMain_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto MakeBedlikeMain_ret_NOMEM;
    }

    ctx.variant_include = variant_include;
    mcp->refalt1_select = refalt1_select;
    mcp->sample_include = sample_include;
    mcp->sample_ct = sample_ct;
    ctx.mcp = mcp;
    ctx.reterr = kPglRetSuccess;
    SetThreadFuncAndData(MakeBedlikeThread, &ctx, &tg);

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
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t prev_variant_idx = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto MakeBedlikeMain_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = ctx.reterr;
        if (unlikely(reterr)) {
          // this should only be possible in MakePgenRobust()
          assert(reterr != kPglRetWriteFail);
          goto MakeBedlikeMain_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_write_ct = cur_block_write_ct;
        ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_write_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto MakeBedlikeMain_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        if (unlikely(fwrite_checked(ctx.writebufs[parity], (variant_idx - prev_variant_idx) * sample_ct4, outfile))) {
          goto MakeBedlikeMain_ret_WRITE_FAIL;
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
      // crucially, this is independent of the PgenReader block_base pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(fclose_null(&outfile))) {
      goto MakeBedlikeMain_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    // BigstackReset(bigstack_mark);
  }
  while (0) {
  MakeBedlikeMain_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakeBedlikeMain_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  MakeBedlikeMain_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  MakeBedlikeMain_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  MakeBedlikeMain_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 MakeBedlikeMain_ret_1:
  CleanupThreads(&tg);
  fclose_cond(outfile);
  // parent will free memory
  return reterr;
}

typedef struct MakePgenCtxStruct {
  MakeCommon* mcp;

  const uint32_t* new_sample_idx_to_old;
  uint32_t* old_sample_idx_to_new;
  // combine existing chr_mask/xymt_codes/haploid_mask/chr_idx_to_foidx with
  // new collapsed chromosome boundary table
  uint32_t* write_chr_fo_vidx_start;
  const uintptr_t* write_allele_idx_offsets;
  const uintptr_t* sex_male_collapsed;
  uintptr_t* sex_female_collapsed;
  uint32_t dosage_erase_halfdist;

  uintptr_t** loadbuf_thread_starts[2];
  // phase, dosage
  unsigned char* loaded_vrtypes[2];

  uint32_t cur_block_write_ct;

  STPgenWriter* spgwp;
  PgenWriterCommon** pwcs;
  uintptr_t** thread_write_genovecs;
  uintptr_t** thread_write_mhc;
  // AlleleCode** thread_ac_rotate;
  uintptr_t** thread_write_phasepresents;
  uintptr_t** thread_write_phaseinfos;
  uintptr_t** thread_all_hets;
  uintptr_t** thread_write_dosagepresents;
  Dosage** thread_write_dosagevals;
  uintptr_t** thread_write_dphasepresents;
  SDosage** thread_write_dphasedeltas;
  uint32_t** thread_cumulative_popcount_bufs;
  PglErr write_reterr;
  int32_t write_errno;
} MakePgenCtx;

// One-thread-per-vblock is sensible for possibly-phased biallelic data, where
// subsetting and LD-compression are a substantial fraction of processing time,
// and memory requirements tend to be low enough that it's actually reasonable
// for each thread job to comprise 64k variants.
// Beyond that... the VCF/.pgen division of labor looks nice, but far too much
// of the work is usually being done in the initial PgrGetRaw() call, so just
// fall back on single-threaded invocation of the same function; only
// difference is that the worker thread owns the writer object.
THREAD_FUNC_DECL MakePgenThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  MakePgenCtx* ctx = S_CAST(MakePgenCtx*, arg->sharedp->context);

  const uint32_t* new_sample_idx_to_old = ctx->new_sample_idx_to_old;
  const uint32_t* old_sample_idx_to_new = ctx->old_sample_idx_to_new;
  const MakeCommon* mcp = ctx->mcp;
  const ChrInfo* cip = mcp->cip;
  const uint32_t* write_chr_fo_vidx_start = ctx->write_chr_fo_vidx_start;
  const uintptr_t* write_allele_idx_offsets = ctx->write_allele_idx_offsets;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select_iter) = mcp->refalt1_select;
  const uintptr_t* sample_include = mcp->sample_include;

  const uintptr_t* sex_male_collapsed = ctx->sex_male_collapsed;

  const uintptr_t* sex_male_collapsed_interleaved = mcp->sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed = ctx->sex_female_collapsed;
  const uintptr_t* sex_female_collapsed_interleaved = mcp->sex_female_collapsed_interleaved;
  const uint32_t raw_sample_ct = mcp->raw_sample_ct;
  const uint32_t sample_ct = mcp->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t sample_ctv2 = NypCtToVecCt(sample_ct);
  const uint32_t raw_sample_ctaw2 = NypCtToAlignedWordCt(raw_sample_ct);
  const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];

  const Plink2WriteFlags plink2_write_flags = mcp->plink2_write_flags;
  const uint32_t set_hh_missing = plink2_write_flags & kfPlink2WriteSetHhMissing;
  const uint32_t set_hh_missing_keep_dosage = plink2_write_flags & kfPlink2WriteSetHhMissingKeepDosage;
  const uint32_t set_mixed_mt_missing = plink2_write_flags & kfPlink2WriteSetMixedMtMissing;
  const uint32_t set_mixed_mt_missing_keep_dosage = plink2_write_flags & kfPlink2WriteSetMixedMtMissingKeepDosage;
  const uint32_t late_dosage_erase = plink2_write_flags & kfPlink2WriteLateDosageErase;
  const uint32_t fill_missing_from_dosage = plink2_write_flags & kfPlink2WriteFillMissingFromDosage;

  const uint32_t hard_call_halfdist = mcp->hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = ctx->dosage_erase_halfdist;
  const uintptr_t dosageraw_word_ct = kWordsPerVec * (BitCtToVecCt(raw_sample_ct) + DivUp(raw_sample_ct, (kBytesPerVec / sizeof(Dosage))));

  STPgenWriter* spgwp = ctx->spgwp;
  PgenWriterCommon* pwcp;
  if (spgwp) {
    // make this function stand out as an intrusive one
    pwcp = &GET_PRIVATE(*spgwp, pwc);
  } else {
    pwcp = ctx->pwcs[tidx];
  }
  uintptr_t* write_genovec = nullptr;
  // assumes sample_include == nullptr if sample_ct == raw_sample_ct
  if (new_sample_idx_to_old || sample_include) {
    write_genovec = ctx->thread_write_genovecs[tidx];
    write_genovec[sample_ctl2 - 1] = 0;
  }
  uintptr_t* write_patch_01_set = nullptr;
  AlleleCode* write_patch_01_vals = nullptr;
  uintptr_t* write_patch_10_set = nullptr;
  AlleleCode* write_patch_10_vals = nullptr;
  if (ctx->thread_write_mhc) {
    ExpandMhc(sample_ct, ctx->thread_write_mhc[tidx], &write_patch_01_set, &write_patch_01_vals, &write_patch_10_set, &write_patch_10_vals);
  }
  uintptr_t* write_phasepresent = nullptr;
  uintptr_t* write_phaseinfo = nullptr;
  uintptr_t* all_hets = nullptr;
  if (ctx->thread_write_phasepresents) {
    write_phasepresent = ctx->thread_write_phasepresents[tidx];
    write_phaseinfo = ctx->thread_write_phaseinfos[tidx];
    if (ctx->thread_all_hets) {
      all_hets = ctx->thread_all_hets[tidx];
    }
  }
  uintptr_t* write_dosagepresent = nullptr;
  Dosage* write_dosagevals = nullptr;
  uintptr_t* write_dphasepresent = nullptr;
  SDosage* write_dphasedeltas = nullptr;
  SDosage* tmp_dphasedeltas = nullptr;
  uint32_t* cumulative_popcount_buf = nullptr;
  if (ctx->thread_write_dosagepresents) {
    write_dosagepresent = ctx->thread_write_dosagepresents[tidx];
    write_dosagevals = ctx->thread_write_dosagevals[tidx];
    if (ctx->thread_write_dphasepresents) {
      write_dphasepresent = ctx->thread_write_dphasepresents[tidx];
      write_dphasedeltas = ctx->thread_write_dphasedeltas[tidx];
      tmp_dphasedeltas = &(write_dphasedeltas[RoundUpPow2(sample_ct, kCacheline / 2)]);
    }
  }
  if ((ctx->thread_write_mhc || ctx->thread_write_dosagepresents) && new_sample_idx_to_old) {
    cumulative_popcount_buf = ctx->thread_cumulative_popcount_bufs[tidx];
  }
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t parity = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t write_idx = tidx * kPglVblockSize;
    const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
    uintptr_t* loadbuf_iter = ctx->loadbuf_thread_starts[parity][tidx];
    unsigned char* loaded_vrtypes = ctx->loaded_vrtypes[parity];
    uint32_t loaded_vrtype = 0;
    uint32_t chr_end_bidx = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t is_haploid_nonmt = 0;
    uint32_t is_mt = 0;
    // write_idx may start larger than write_idx_end
    for (; write_idx < write_idx_end; ++write_idx) {
      if (loaded_vrtypes) {
        loaded_vrtype = loaded_vrtypes[write_idx];
      }
      if (write_idx >= chr_end_bidx) {
        const uint32_t chr_fo_idx = CountSortedSmallerU32(&(write_chr_fo_vidx_start[1]), cip->chr_ct, write_idx + variant_idx_offset + 1);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end_bidx = write_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
        is_mt = (chr_idx == mt_code);
        is_haploid_nonmt = IsSet(cip->haploid_mask, chr_idx) && (!is_mt);
      }
      uintptr_t* cur_genovec_end = &(loadbuf_iter[raw_sample_ctaw2]);
      if (write_allele_idx_offsets) {
        allele_ct = write_allele_idx_offsets[write_idx + variant_idx_offset + 1] - write_allele_idx_offsets[write_idx + variant_idx_offset];
      }
      const uint32_t is_mhc = loaded_vrtype & 8;
      uint32_t read_rare01_ct = 0;
      uint32_t read_rare10_ct = 0;
      uintptr_t* read_patch_01_set = nullptr;
      AlleleCode* read_patch_01_vals = nullptr;
      uintptr_t* read_patch_10_set = nullptr;
      AlleleCode* read_patch_10_vals = nullptr;
      if (is_mhc) {
        assert(allele_ct > 2);
        read_rare01_ct = cur_genovec_end[0];
        read_rare10_ct = cur_genovec_end[1];
        cur_genovec_end = &(cur_genovec_end[RoundUpPow2(2, kWordsPerVec)]);
        if (read_rare01_ct) {
          read_patch_01_set = cur_genovec_end;
          cur_genovec_end = &(cur_genovec_end[raw_sample_ctl]);
          read_patch_01_vals = R_CAST(AlleleCode*, cur_genovec_end);
          cur_genovec_end = &(cur_genovec_end[DivUp(read_rare01_ct, kBytesPerWord / sizeof(AlleleCode))]);
          VecAlignUp64(&cur_genovec_end);
        }
        if (read_rare10_ct) {
          read_patch_10_set = cur_genovec_end;
          cur_genovec_end = &(cur_genovec_end[raw_sample_ctl]);
          read_patch_10_vals = R_CAST(AlleleCode*, cur_genovec_end);
          cur_genovec_end = &(cur_genovec_end[DivUp(read_rare10_ct, kBytesPerWord / (2 * sizeof(AlleleCode)))]);
          VecAlignUp64(&cur_genovec_end);
        }
      }
      uint32_t is_hphase = loaded_vrtype & 0x10;
      uintptr_t* cur_phaseraw = nullptr;
      if (is_hphase) {
        // tried skipping this and using ExpandThenSubsetBytearr in simplest
        // case, not worthwhile
        if (!read_rare10_ct) {
          PgrDetectGenoarrHets(loadbuf_iter, raw_sample_ct, all_hets);
        } else {
          PgrDetectGenoarrHetsMultiallelic(loadbuf_iter, read_patch_10_set, read_patch_10_vals, raw_sample_ct, all_hets);
        }
        cur_phaseraw = cur_genovec_end;
        const uint32_t het_ct = S_CAST(uint32_t, cur_phaseraw[0]);
#ifdef __LP64__
        const uint32_t explicit_phasepresent_ct = cur_phaseraw[0] >> 32;
#else
        const uint32_t explicit_phasepresent_ct = cur_phaseraw[1];
#endif
        const uint32_t phaseraw_word_ct = (8 / kBytesPerWord) + 1 + (het_ct / kBitsPerWord) + DivUp(explicit_phasepresent_ct, kBitsPerWord);
        cur_genovec_end = &(cur_genovec_end[RoundUpPow2(phaseraw_word_ct, kWordsPerVec)]);
      }
      const uint32_t is_dosage = loaded_vrtype & 0x60;
      const uint32_t is_dphase = loaded_vrtype & 0x80;
      uintptr_t* cur_write_phasepresent = write_phasepresent;
      uintptr_t* cur_dosagepresent = nullptr;
      Dosage* cur_dosagevals = nullptr;
      uintptr_t* cur_dphasepresent = nullptr;
      SDosage* cur_dphasedelta = nullptr;
      uint32_t read_dosage_ct = 0;
      uint32_t read_dphase_ct = 0;
      if (is_dosage) {
        // multiallelic dosage not implemented yet
        assert(allele_ct == 2);

        // this should have length dependent on dosage_ct
        cur_dosagepresent = cur_genovec_end;
        cur_dosagevals = R_CAST(Dosage*, &(cur_dosagepresent[raw_sample_ctaw]));
        read_dosage_ct = PopcountWords(cur_dosagepresent, raw_sample_ctl);

        // temporary
        cur_genovec_end = &(cur_genovec_end[dosageraw_word_ct]);

        if (is_dphase) {
          cur_dphasepresent = cur_genovec_end;
          cur_dphasedelta = R_CAST(SDosage*, &(cur_dphasepresent[raw_sample_ctaw]));
          read_dphase_ct = PopcountWords(cur_dphasepresent, raw_sample_ctl);

          // temporary
          cur_genovec_end = &(cur_genovec_end[dosageraw_word_ct]);
        }
      }
      uint32_t write_rare01_ct = 0;
      uint32_t write_rare10_ct = 0;
      uint32_t write_dosage_ct = 0;
      uint32_t write_dphase_ct = 0;
      if (new_sample_idx_to_old) {
        GenovecResort(loadbuf_iter, new_sample_idx_to_old, sample_ct, write_genovec);
        if (read_rare01_ct) {
          write_rare01_ct = CopyAndResort8bit(read_patch_01_set, read_patch_01_vals, new_sample_idx_to_old, raw_sample_ct, sample_ct, write_patch_01_set, write_patch_01_vals, cumulative_popcount_buf);
        }
        if (read_rare10_ct) {
          write_rare10_ct = CopyAndResort16bit(read_patch_10_set, read_patch_10_vals, new_sample_idx_to_old, raw_sample_ct, sample_ct, write_patch_10_set, write_patch_10_vals, cumulative_popcount_buf);
        }
        if (is_hphase) {
          UnpackAndResortHphase(all_hets, cur_phaseraw, sample_include, old_sample_idx_to_new, raw_sample_ct, sample_ct, &cur_write_phasepresent, write_phaseinfo);
        }
        if (is_dosage) {
          write_dosage_ct = CopyAndResort16bit(cur_dosagepresent, cur_dosagevals, new_sample_idx_to_old, raw_sample_ct, sample_ct, write_dosagepresent, write_dosagevals, cumulative_popcount_buf);
          if (is_dphase) {
            write_dphase_ct = CopyAndResort16bit(cur_dphasepresent, cur_dphasedelta, new_sample_idx_to_old, raw_sample_ct, sample_ct, write_dphasepresent, write_dphasedeltas, cumulative_popcount_buf);
          }
        }
      } else if (sample_include) {
        CopyNyparrNonemptySubset(loadbuf_iter, sample_include, raw_sample_ct, sample_ct, write_genovec);
        if (is_mhc) {
          write_rare01_ct = Copy1bit8Subset(read_patch_01_set, read_patch_01_vals, sample_include, read_rare01_ct, sample_ct, write_patch_01_set, write_patch_01_vals);
          write_rare10_ct = Copy1bit16Subset(read_patch_10_set, read_patch_10_vals, sample_include, read_rare10_ct, sample_ct, write_patch_10_set, write_patch_10_vals);
        }
        if (is_hphase) {
          UnpackHphaseSubset(all_hets, cur_phaseraw, sample_include, sample_ct, &cur_write_phasepresent, write_phaseinfo);
        }
        if (is_dosage) {
          write_dosage_ct = Copy1bit16Subset(cur_dosagepresent, cur_dosagevals, sample_include, read_dosage_ct, sample_ct, write_dosagepresent, write_dosagevals);
          if (is_dphase) {
            write_dphase_ct = Copy1bit16Subset(cur_dphasepresent, cur_dphasedelta, sample_include, read_dphase_ct, sample_ct, write_dphasepresent, write_dphasedeltas);
          }
        }
      } else {
        write_genovec = loadbuf_iter;
        if (is_mhc) {
          // this doesn't work in refalt1_select case
          write_patch_01_set = read_patch_01_set;
          write_patch_01_vals = read_patch_01_vals;
          write_patch_10_set = read_patch_10_set;
          write_patch_10_vals = read_patch_10_vals;
          write_rare01_ct = read_rare01_ct;
          write_rare10_ct = read_rare10_ct;
        }
        if (is_hphase) {
          UnpackHphase(all_hets, cur_phaseraw, sample_ct, &cur_write_phasepresent, write_phaseinfo);
        }
        if (is_dosage) {
          CopyDosage(cur_dosagepresent, cur_dosagevals, sample_ct, read_dosage_ct, write_dosagepresent, write_dosagevals, &write_dosage_ct);
          if (is_dphase) {
            CopyDosage(cur_dphasepresent, R_CAST(Dosage*, cur_dphasedelta), sample_ct, read_dphase_ct, write_dphasepresent, R_CAST(Dosage*, write_dphasedeltas), &write_dphase_ct);
          }
        }
      }
      // multiallelic -> biallelic split:
      //   main thread will probably compute split mapping in advance (bitarray
      //   with filtered-and-split variant indices, set bit = unsplit variant
      //   or last variant in a split group)?  same pre-split variant can be
      //   loaded multiple times.
      // biallelic -> multiallelic merge:
      //   could require no multiallelic variants in remainder of dataset?
      //   if handled with pgenlib, PgfiMultiread,
      //   PgfiMultireadGetCachelineReq, and PgrGetRaw would need to be
      //   extended to take a merge-info parameter.  +both will be tricky...
      //   probably not worth it.
      //   compute merge pattern in MakePgenRobust() before main loop instead.
      //   main thread also performs the actual merge.
      // both should require sorted .pvar.
      // neither should require any handling in this function.
      if (refalt1_select_iter && (refalt1_select_iter[write_idx][0] || (refalt1_select_iter[write_idx][1] != 1))) {
        if (allele_ct == 2) {
          GenovecInvertUnsafe(sample_ct, write_genovec);
          if (is_hphase) {
            // trailing bits don't matter
            BitvecInvert(sample_ctl, write_phaseinfo);
          }
          if (write_dosage_ct) {
            BiallelicDosage16Invert(write_dosage_ct, write_dosagevals);
            if (write_dphase_ct) {
              BiallelicDphase16Invert(write_dphase_ct, write_dphasedeltas);
            }
          }
        } else {
          exit(S_CAST(int32_t, kPglRetNotYetSupported));
          // this is the fun case
          // 1. fill length-(2 * sample_ct) AlleleCode[] buffer with codes
          // 2. fill lookup table describing remapping
          // 3. replace elements of table
          // 4. normalize order of each code pair, inverting a phaseinfo bit on
          //    each swap
          // 5. call PglMultiallelicDenseToSparse to write back
        }
      }
      if (write_dosage_ct) {
        assert((!write_rare01_ct) && (!write_rare10_ct));
        if (hard_call_halfdist || fill_missing_from_dosage || (dosage_erase_halfdist < kDosage4th)) {
          if (is_hphase && (!cur_write_phasepresent)) {
            // explicit phasepresent required for these
            cur_write_phasepresent = write_phasepresent;
            // unsafe to just copy all_hets, because we may have resorted
            // todo: multiallelic dosage
            PgrDetectGenoarrHets(write_genovec, sample_ct, write_phasepresent);
          }
          if (write_dphasepresent && is_hphase && (!write_dphase_ct)) {
            // bugfix (29 Apr 2019): write_dphasepresent not guaranteed to be
            // non-null.
            ZeroWArr(sample_ctl, write_dphasepresent);
          }
          if (hard_call_halfdist) {
            if ((!is_hphase) && (!write_dphase_ct)) {
              ApplyHardCallThresh(write_dosagepresent, write_dosagevals, write_dosage_ct, hard_call_halfdist, write_genovec);
            } else {
              if (!is_hphase) {
                ZeroWArr(sample_ctl, write_phasepresent);
              }
              write_dphase_ct = ApplyHardCallThreshPhased(write_dosagepresent, write_dosagevals, write_dosage_ct, hard_call_halfdist, write_genovec, write_phasepresent, write_phaseinfo, write_dphasepresent, write_dphasedeltas, tmp_dphasedeltas);
              is_hphase = !AllWordsAreZero(write_phasepresent, sample_ctl);
            }
          }
          if (fill_missing_from_dosage) {
            if ((!is_hphase) && (!write_dphase_ct)) {
              FillMissingFromBiallelicDosage(write_dosagepresent, write_dosagevals, write_dosage_ct, write_genovec);
            } else {
              if (!is_hphase) {
                ZeroWArr(sample_ctl, write_phasepresent);
              }
              write_dphase_ct = FillMissingFromBiallelicDosagePhased(write_dosagepresent, write_dosagevals, write_dosage_ct, write_genovec, write_phasepresent, write_phaseinfo, write_dphasepresent, write_dphasedeltas);
              is_hphase = !AllWordsAreZero(write_phasepresent, sample_ctl);
            }
          }
          if (dosage_erase_halfdist < kDosage4th) {
            if (!is_hphase) {
              ZeroWArr(sample_ctl, write_phasepresent);
            }
            uint32_t dosage_read_idx = 0;
            uintptr_t sample_widx = 0;
            uintptr_t cur_bits = write_dosagepresent[0];
            uint32_t dosage_write_idx;
            if (!write_dphase_ct) {
              // If hardcall-phase and dosage present, threshold/2 applies
              // thanks to implicit dosage-phase value
              // const uint32_t dosage_erase_halfdist2 = (dosage_erase_halfdist + kDosage4th + 1) / 2;
              const uint32_t halfdist_extra = (kDosage4th + 1 - dosage_erase_halfdist) / 2;
              for (; dosage_read_idx != write_dosage_ct; ++dosage_read_idx) {
                const uint32_t sample_uidx_lowbits = BitIter1x(write_dosagepresent, &sample_widx, &cur_bits);
                const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
                const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
                if (halfdist >= dosage_erase_halfdist + ((write_phasepresent[sample_widx] >> sample_uidx_lowbits) & 1) * halfdist_extra) {
                  write_dosagepresent[sample_widx] ^= k1LU << sample_uidx_lowbits;
                  break;
                }
              }
              dosage_write_idx = dosage_read_idx;
              while (++dosage_read_idx < write_dosage_ct) {
                const uint32_t sample_uidx_lowbits = BitIter1x(write_dosagepresent, &sample_widx, &cur_bits);
                const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
                const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
                if (halfdist < dosage_erase_halfdist + ((write_phasepresent[sample_widx] >> sample_uidx_lowbits) & 1) * halfdist_extra) {
                  write_dosagevals[dosage_write_idx++] = dosage_int;
                } else {
                  write_dosagepresent[sample_widx] ^= k1LU << sample_uidx_lowbits;
                }
              }
            } else {
              // Only erase dosage if both sides are less than threshold/2
              // away from an integer.
              const uint32_t halfdist_extra = (kDosage4th + 1 - dosage_erase_halfdist) / 2;
              const uint32_t dosage_erase_halfdist2 = dosage_erase_halfdist + halfdist_extra;
              uint32_t dphase_read_idx = 0;
              uintptr_t lowbit = 0;
              for (; dosage_read_idx != write_dosage_ct; ++dosage_read_idx) {
                lowbit = BitIter1y(write_dosagepresent, &sample_widx, &cur_bits);
                const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
                if (!(write_dphasepresent[sample_widx] & lowbit)) {
                  // necessary for this to be separate to handle odd
                  // dosage_int, missing phase case correctly
                  const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
                  if (halfdist >= dosage_erase_halfdist + ((write_phasepresent[sample_widx] & lowbit) != 0) * halfdist_extra) {
                    break;
                  }
                } else {
                  const int32_t dphase_delta = write_dphasedeltas[dphase_read_idx++];
                  const uint32_t halfdist1 = HaploidDosageHalfdist((dosage_int + dphase_delta) >> 1);
                  const uint32_t halfdist2 = HaploidDosageHalfdist((dosage_int - dphase_delta) >> 1);
                  if ((halfdist1 >= dosage_erase_halfdist2) && (halfdist2 >= dosage_erase_halfdist2)) {
                    break;
                  }
                }
              }
              dosage_write_idx = dosage_read_idx;
              if (dosage_read_idx < write_dosage_ct) {
                uint32_t dphase_write_idx = dphase_read_idx;
                if (write_dphasepresent[sample_widx] & lowbit) {
                  --dphase_write_idx;
                  write_dphasepresent[sample_widx] ^= lowbit;
                }
                write_dosagepresent[sample_widx] ^= lowbit;
                while (++dosage_read_idx < write_dosage_ct) {
                  lowbit = BitIter1y(write_dosagepresent, &sample_widx, &cur_bits);
                  const uint32_t dosage_int = write_dosagevals[dosage_read_idx];
                  if (!(write_dphasepresent[sample_widx] & lowbit)) {
                    const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
                    if (halfdist < dosage_erase_halfdist + ((write_phasepresent[sample_widx] & lowbit) != 0) * halfdist_extra) {
                      write_dosagevals[dosage_write_idx++] = dosage_int;
                    } else {
                      write_dosagepresent[sample_widx] ^= lowbit;
                    }
                  } else {
                    const int32_t dphase_delta = write_dphasedeltas[dphase_read_idx++];
                    const uint32_t halfdist1 = HaploidDosageHalfdist((dosage_int + dphase_delta) >> 1);
                    const uint32_t halfdist2 = HaploidDosageHalfdist((dosage_int - dphase_delta) >> 1);
                    if ((halfdist1 < dosage_erase_halfdist2) || (halfdist2 < dosage_erase_halfdist2)) {
                      write_dosagevals[dosage_write_idx++] = dosage_int;
                      write_dphasedeltas[dphase_write_idx++] = dphase_delta;
                    } else {
                      write_dosagepresent[sample_widx] ^= lowbit;
                      write_dphasepresent[sample_widx] ^= lowbit;
                    }
                  }
                }
                write_dphase_ct = dphase_write_idx;
              }
            }
            write_dosage_ct = dosage_write_idx;
          }
        }
        if (late_dosage_erase) {
          write_dosage_ct = 0;
          write_dphase_ct = 0;
        }
      }
      // moved after --hard-call-threshold, since it makes sense to
      // immediately erase fresh het haploid calls
      if (set_hh_missing && is_haploid_nonmt) {
        if (is_x) {
          EraseMaleDphases(sex_male_collapsed, &write_dphase_ct, write_dphasepresent, write_dphasedeltas);
          if (!set_hh_missing_keep_dosage) {
            // need to erase dosages associated with the hardcalls we're
            // about to clear

            // male 0/x hets to missing
            SetMaleHetMissingCleardosage(sex_male_collapsed, sex_male_collapsed_interleaved, sample_ctv2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
            // male x/y hets to missing
            if (write_rare10_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t patch_10_bits = write_patch_10_set[0];
              uint32_t read_patch_10_idx = 0;
              for (; read_patch_10_idx != write_rare10_ct; ++read_patch_10_idx) {
                uintptr_t lowbit = BitIter1y(write_patch_10_set, &sample_widx, &patch_10_bits);
                AlleleCode lo_code = write_patch_10_vals[read_patch_10_idx * 2];
                AlleleCode hi_code = write_patch_10_vals[read_patch_10_idx * 2 + 1];
                if ((sex_male_collapsed[sample_widx] & lowbit) && (lo_code != hi_code)) {
                  write_patch_10_set[sample_widx] ^= lowbit;
                  uint32_t write_patch_10_idx = read_patch_10_idx;
                  ++read_patch_10_idx;
                  for (; read_patch_10_idx != write_rare10_ct; ++read_patch_10_idx) {
                    lowbit = BitIter1y(write_patch_10_set, &sample_widx, &patch_10_bits);
                    lo_code = write_patch_10_vals[read_patch_10_idx * 2];
                    hi_code = write_patch_10_vals[read_patch_10_idx * 2 + 1];
                    if ((sex_male_collapsed[sample_widx] & lowbit) && (lo_code != hi_code)) {
                      write_patch_10_set[sample_widx] ^= lowbit;
                    } else {
                      write_patch_10_vals[write_patch_10_idx * 2] = lo_code;
                      write_patch_10_vals[write_patch_10_idx * 2 + 1] = hi_code;
                      ++write_patch_10_idx;
                    }
                  }
                  write_rare10_ct = write_patch_10_idx;
                  break;
                }
              }
            }
          } else {
            assert(!write_rare01_ct);
            assert(!write_rare10_ct);
            // need to generate a new unphased dosage for each cleared
            // hardcall lacking a dosage entry
            SetMaleHetMissingKeepdosage(sex_male_collapsed, sex_male_collapsed_interleaved, sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          }
          if (is_hphase && cur_write_phasepresent) {
            // bugfix (28 Jul 2018): I was on crack when I moved this code
            // before SetMaleHetMissing{Clear,Keep}dosage() on 31 Mar
            if (!write_rare10_ct) {
              MaskGenoarrHetsUnsafe(write_genovec, sample_ctl2, cur_write_phasepresent);
            } else {
              MaskGenoarrHetsMultiallelicUnsafe(write_genovec, write_patch_10_set, write_patch_10_vals, sample_ctl2, cur_write_phasepresent);
            }
            is_hphase = !AllWordsAreZero(write_phasepresent, sample_ctl);
          }
          if (write_rare01_ct) {
            ClearGenoarrMissing1bit8Unsafe(write_genovec, &write_rare01_ct, write_patch_01_set, write_patch_01_vals);
          }
          if (write_rare10_ct) {
            ClearGenoarrMissing1bit16Unsafe(write_genovec, &write_rare10_ct, write_patch_10_set, write_patch_10_vals);
          }
        } else {
          // all hets to missing
          // may want to move is_hphase zeroing in front
          if (!set_hh_missing_keep_dosage) {
            SetHetMissingCleardosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          } else {
            SetHetMissingKeepdosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          }
          if (is_y) {
            InterleavedSetMissingCleardosage(sex_female_collapsed, sex_female_collapsed_interleaved, sample_ctv2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
          }
          is_hphase = 0;
          write_rare01_ct = 0;
          if (write_rare10_ct) {
            ClearGenoarrMissing1bit16Unsafe(write_genovec, &write_rare10_ct, write_patch_10_set, write_patch_10_vals);
          }
          write_dphase_ct = 0;
        }
      } else if (set_mixed_mt_missing && is_mt) {
        if (!set_mixed_mt_missing_keep_dosage) {
          // all hets to missing
          SetHetMissingCleardosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
        } else {
          SetHetMissingKeepdosage(sample_ctl2, write_genovec, &write_dosage_ct, write_dosagepresent, write_dosagevals);
        }
        is_hphase = 0;
        write_rare01_ct = 0;
        if (write_rare10_ct) {
          ClearGenoarrMissing1bit16Unsafe(write_genovec, &write_rare10_ct, write_patch_10_set, write_patch_10_vals);
        }
        write_dphase_ct = 0;
      }
      ZeroTrailingNyps(sample_ct, write_genovec);
      // todo: --set-me-missing, --zero-cluster, --fill-missing-with-ref
      if (spgwp) {
        if (pwcp->fwrite_bufp >= &(pwcp->fwrite_buf[kPglFwriteBlockSize])) {
          const uintptr_t cur_byte_ct = pwcp->fwrite_bufp - pwcp->fwrite_buf;
          if (unlikely(fwrite_checked(pwcp->fwrite_buf, cur_byte_ct, GET_PRIVATE(*spgwp, pgen_outfile)))) {
            ctx->write_reterr = kPglRetWriteFail;
            ctx->write_errno = errno;
            break;
          }
          // printf("vblock_fpos_offset: %llu\n", pwcp->vblock_fpos_offset);
          pwcp->vblock_fpos_offset += cur_byte_ct;
          // printf("%u %llu\n", write_idx + variant_idx_offset, pwcp->vblock_fpos_offset);
          pwcp->fwrite_bufp = pwcp->fwrite_buf;
        }
      }
      if ((!write_rare01_ct) && (!write_rare10_ct)) {
        if ((!is_hphase) && (!write_dphase_ct)) {
          if (unlikely(PwcAppendBiallelicGenovecDosage16(write_genovec, write_dosagepresent, write_dosagevals, write_dosage_ct, pwcp))) {
            ctx->write_reterr = kPglRetVarRecordTooLarge;
            break;
          }
        } else {
          if (!is_hphase) {
            ZeroWArr(sample_ctl, write_phasepresent);
          }
          // extraneous phaseinfo bits may be set
          if (unlikely(PwcAppendBiallelicGenovecDphase16(write_genovec, cur_write_phasepresent, write_phaseinfo, write_dosagepresent, write_dphasepresent, write_dosagevals, write_dphasedeltas, write_dosage_ct, write_dphase_ct, pwcp))) {
            ctx->write_reterr = kPglRetVarRecordTooLarge;
            break;
          }
        }
      } else {
        // multiallelic dosage not supported
        if (!is_hphase) {
          if (unlikely(PwcAppendMultiallelicSparse(write_genovec, write_patch_01_set, write_patch_01_vals, write_patch_10_set, write_patch_10_vals, write_rare01_ct, write_rare10_ct, pwcp))) {
            ctx->write_reterr = kPglRetVarRecordTooLarge;
            break;
          }
        } else {
          if (unlikely(PwcAppendMultiallelicGenovecHphase(write_genovec, write_patch_01_set, write_patch_01_vals, write_patch_10_set, write_patch_10_vals, cur_write_phasepresent, write_phaseinfo, write_rare01_ct, write_rare10_ct, pwcp))) {
            ctx->write_reterr = kPglRetVarRecordTooLarge;
            break;
          }
        }
      }
      loadbuf_iter = cur_genovec_end;
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_write_ct;
    if (refalt1_select_iter) {
      refalt1_select_iter = &(refalt1_select_iter[cur_block_write_ct]);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PgenGlobalFlags GflagsVfilter(const uintptr_t* variant_include, const unsigned char* vrtypes, uint32_t raw_variant_ct, PgenGlobalFlags input_gflags) {
  PgenGlobalFlags read_gflags = kfPgenGlobal0;
  const uintptr_t* vrtypes_alias = R_CAST(const uintptr_t*, vrtypes);
  const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
  uint32_t mask_multiply = ((input_gflags & kfPgenGlobalHardcallPhasePresent)? 0x10 : 0) + ((input_gflags & kfPgenGlobalDosagePresent)? 0x60 : 0) + ((input_gflags & kfPgenGlobalDosagePhasePresent)? 0x80 : 0);
  uintptr_t vrtypes_or = 0;
  // todo: try changing loop to be vec-based, use movemask to extract
  // information from vrtypes in 64-bit cases
  for (uint32_t widx = 0; widx != raw_variant_ctl; ++widx) {
    uintptr_t cur_variant_include_word = variant_include[widx];
    if (cur_variant_include_word) {
      // bugfix (20 Aug 2018): this needs to advance on every variant_include
      // word, not just the nonzero ones
      const uintptr_t* cur_vrtypes = &(vrtypes_alias[8 * widx]);
#ifdef __LP64__
      for (uint32_t vi_byte_idx = 0; vi_byte_idx != 8; ++vi_byte_idx) {
#  ifdef USE_AVX2
        // this doesn't seem to be much faster than non-AVX2 code on my Mac...
        // inverse-movemask shouldn't be better than regular movemask here
        const uintptr_t cur_mask = _pdep_u64(cur_variant_include_word, kMask0101);
#  else
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
#  endif
        vrtypes_or |= cur_vrtypes[vi_byte_idx] & (cur_mask * mask_multiply);
        cur_variant_include_word >>= 8;
      }
#else
      for (uint32_t vi_hexa_idx = 0; vi_hexa_idx != 8; ++vi_hexa_idx) {
        // dcba -> d0000000c0000000b0000000a
        const uintptr_t cur_mask = ((cur_variant_include_word & 0xf) * 0x204081) & kMask0101;
        vrtypes_or |= cur_vrtypes[vi_hexa_idx] & (cur_mask * mask_multiply);
        cur_variant_include_word >>= 4;
      }
#endif
      if (vrtypes_or) {
        // bugfix (8 Oct 2017): forgot to multiply by kMask0101
        if (vrtypes_or & (0x10 * kMask0101)) {
          read_gflags |= kfPgenGlobalHardcallPhasePresent;
          mask_multiply -= 0x10;
        }
        if (vrtypes_or & (0x60 * kMask0101)) {
          read_gflags |= kfPgenGlobalDosagePresent;
          mask_multiply -= 0x60;
        }
        if (vrtypes_or & (0x80 * kMask0101)) {
          read_gflags |= kfPgenGlobalDosagePhasePresent;
          mask_multiply -= 0x80;
        }
        if (!mask_multiply) {
          return read_gflags;
        }
      }
    }
  }
  return read_gflags;
}

void SplitNonrefFlags() {
  logerrputs("Provisional-reference flag split is not implemented yet.\n");
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
}

void JoinNonrefFlags() {
  logerrputs("Provisional-reference flag join is not implemented yet.\n");
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
}

// Single-output-thread implementation.  Allows variants to be unsorted.
// (Note that MakePlink2NoVsort() currently requires enough memory for 64k * 2
// variants per output thread, due to LD compression.  This is faster in the
// common case, but once you have 150k+ samples with dosage data...)
//
// initialized mcp fields: cip, sex_male_collapsed_interleaved,
// sex_female_collapsed_interleaved, raw_sample_ct, sample_ct,
// plink2_write_flags
PglErr MakePgenRobust(const uintptr_t* sample_include, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, __maybe_unused const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* write_allele_idx_offsets, const uint32_t* new_variant_idx_to_old, const uintptr_t* sex_male_collapsed, uintptr_t* sex_female_collapsed, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t write_variant_ct, uint32_t max_read_allele_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, MakeCommon* mcp, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  // variant_uidx_new_to_old[] can be nullptr

  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  STPgenWriter spgw;
  PreinitSpgw(&spgw);
  MakePgenCtx ctx;
  {
    // plink2_write_flags assumed to include --set-hh-missing and
    //   --set-mixed-mt-missing
    // sex_{fe}male_collapsed_interleaved assumed to be initialized if
    //   necessary

    if (unlikely(SetThreadCt(1, &tg))) {
      goto MakePgenRobust_ret_NOMEM;
    }
    ctx.spgwp = &spgw;
    const uint32_t raw_sample_ct = mcp->raw_sample_ct;
    const uint32_t sample_ct = mcp->sample_ct;
    const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    mcp->sample_include = subsetting_required? sample_include : nullptr;
    ctx.new_sample_idx_to_old = new_sample_idx_to_old;
    ctx.sex_male_collapsed = sex_male_collapsed;
    ctx.sex_female_collapsed = sex_female_collapsed;
    ctx.write_reterr = kPglRetSuccess;
    if ((make_plink2_flags & kfMakeBed) || ((make_plink2_flags & (kfMakePgen | (kfMakePgenFormatBase * 3))) == (kfMakePgen | kfMakePgenFormatBase))) {
      logerrputs("Error: Fixed-width .bed/.pgen output doesn't support sorting yet.  Generate a\nregular sorted .pgen first, and then reformat it.\n");
      reterr = kPglRetNotYetSupported;
      goto MakePgenRobust_ret_1;
    } else {
      const uint32_t input_biallelic = (!allele_idx_offsets);
      // output_biallelic: test write_allele_idx_offsets equality to null
      ctx.write_allele_idx_offsets = write_allele_idx_offsets;
      if ((variant_ct == raw_variant_ct) || new_variant_idx_to_old) {
        ctx.write_chr_fo_vidx_start = mcp->cip->chr_fo_vidx_start;
      } else {
        if (unlikely(AllocAndFillSubsetChrFoVidxStart(variant_include, mcp->cip, &ctx.write_chr_fo_vidx_start))) {
          goto MakePgenRobust_ret_NOMEM;
        }
      }
      PgenGlobalFlags read_gflags = PgrGetGflags(simple_pgrp) & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      if (make_plink2_flags & (kfMakePlink2MJoin | kfMakePlink2EraseAlt2Plus)) {
        logerrputs("Error: multiallelic-join and 'erase-alt2+' modifiers are under development.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePgenRobust_ret_1;
      }
      if (make_plink2_flags & kfMakePgenErasePhase) {
        read_gflags &= ~(kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
      }
      const uint32_t fill_missing_from_dosage = (make_plink2_flags & kfMakePgenFillMissingFromDosage) && (read_gflags & kfPgenGlobalDosagePresent);
      if (fill_missing_from_dosage) {
        mcp->plink2_write_flags |= kfPlink2WriteFillMissingFromDosage;
      }
      if (make_plink2_flags & kfMakePgenEraseDosage) {
        if ((hard_call_thresh == UINT32_MAX) && (!fill_missing_from_dosage)) {
          read_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
        } else {
          // bugfix (11 Apr 2018): this was in the wrong branch
          mcp->plink2_write_flags |= kfPlink2WriteLateDosageErase;
        }
      }
      if (read_gflags && (variant_ct < raw_variant_ct)) {
        read_gflags &= GflagsVfilter(variant_include, PgrGetVrtypes(simple_pgrp), raw_variant_ct, PgrGetGflags(simple_pgrp));
      }
      if (!input_biallelic) {
        // todo: conditional erase-alt2+ exception
        read_gflags |= kfPgenGlobalMultiallelicHardcallFound;
      }
      const uint32_t read_dosage_present = (read_gflags / kfPgenGlobalDosagePresent) & 1;
      // bugfix (25 Jul 2018): left expression needs ||, not &&
      mcp->hard_call_halfdist = ((hard_call_thresh == UINT32_MAX) || (!read_dosage_present))? 0 : (kDosage4th - hard_call_thresh);
      ctx.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      // bugfix/simplification (10 Mar 2020): it is possible for dosage-phase
      // to be present in the input without hardcall-phase.  Don't try to treat
      // that differently than the usual scenario where hardcall-phase is
      // present.
      const uint32_t read_phase_present = !!(read_gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent));
      const uint32_t read_dphase_present = (read_gflags / kfPgenGlobalDosagePhasePresent) & 1;
      PgenGlobalFlags write_gflags = read_gflags;
      // When --hard-call-threshold is specified, if either hphase or dphase
      // values exist, the other can be generated.
      uint32_t read_or_write_phase_present = read_phase_present;
      uint32_t read_or_write_dphase_present = read_dphase_present;
      if ((mcp->hard_call_halfdist || fill_missing_from_dosage) && (read_phase_present || read_or_write_dphase_present)) {
        read_or_write_phase_present = 1;
        read_or_write_dphase_present = 1;
        write_gflags |= kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent;
      } else if (dosage_erase_thresh && read_dosage_present) {
        // need write_phasepresent, pretty harmless to allocate write_phaseinfo
        read_or_write_phase_present = 1;
      }
      uint32_t read_or_write_dosage_present = read_dosage_present;
      if (mcp->plink2_write_flags & kfPlink2WriteLateDosageErase) {
        write_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      } else if (mcp->plink2_write_flags & (kfPlink2WriteSetHhMissingKeepDosage | kfPlink2WriteSetMixedMtMissingKeepDosage)) {
        // bugfix (25 Jul 2018): this needs to check plink2_write_flags, not
        // make_plink2_flags

        // command-line parser guarantees erase-dosage and
        // --set-hh-missing/--set-mixed-mt-missing keep-dosage aren't used
        // together
        read_or_write_dosage_present = 1;

        // could verify at least one het haploid is present before setting this
        // flag...
        write_gflags |= kfPgenGlobalDosagePresent;
      }
      if ((write_gflags & (kfPgenGlobalMultiallelicHardcallFound | kfPgenGlobalDosagePresent)) == (kfPgenGlobalMultiallelicHardcallFound | kfPgenGlobalDosagePresent)) {
        logerrputs("Error: Multiallelic dosages aren't supported yet.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePgenRobust_ret_1;
      }

      uint32_t nonref_flags_storage = 3;
      uintptr_t* nonref_flags_write = PgrGetNonrefFlags(simple_pgrp);
      if (!nonref_flags_write) {
        nonref_flags_storage = (PgrGetGflags(simple_pgrp) & kfPgenGlobalAllNonref)? 2 : 1;
      } else if (variant_ct < raw_variant_ct) {
        const uint32_t write_variant_ctl = BitCtToWordCt(write_variant_ct);
        uintptr_t* old_nonref_flags = nonref_flags_write;
        if (bigstack_alloc_w(write_variant_ctl, &nonref_flags_write)) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if ((variant_ct == write_variant_ct) && (!new_variant_idx_to_old)) {
          CopyBitarrSubset(old_nonref_flags, variant_include, variant_ct, nonref_flags_write);
        } else {
          ZeroWArr(write_variant_ctl, nonref_flags_write);
          if (variant_ct == write_variant_ct) {
            for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
              const uintptr_t variant_uidx = new_variant_idx_to_old[variant_idx];
              if (IsSet(old_nonref_flags, variant_uidx)) {
                SetBit(variant_idx, nonref_flags_write);
              }
            }
          } else if (!write_allele_idx_offsets) {
            SplitNonrefFlags();
          } else {
            JoinNonrefFlags();
          }
        }
        if (nonref_flags_write[0] & 1) {
          if (AllBitsAreOne(nonref_flags_write, write_variant_ct)) {
            BigstackReset(nonref_flags_write);
            nonref_flags_write = nullptr;
            nonref_flags_storage = 2;
          }
        } else if (AllWordsAreZero(nonref_flags_write, write_variant_ctl)) {
          BigstackReset(nonref_flags_write);
          nonref_flags_write = nullptr;
          nonref_flags_storage = 1;
        }
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, write_allele_idx_offsets, nonref_flags_write, write_variant_ct, sample_ct, write_gflags, nonref_flags_storage, ctx.spgwp, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto MakePgenRobust_ret_1;
      }
      unsigned char* spgw_alloc;
      if (unlikely(
              bigstack_alloc_wp(1, &(ctx.loadbuf_thread_starts[0])) ||
              bigstack_alloc_wp(1, &(ctx.loadbuf_thread_starts[1])) ||
              bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
        goto MakePgenRobust_ret_NOMEM;
      }
      SpgwInitPhase2(max_vrec_len, ctx.spgwp, spgw_alloc);

      const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
      const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
      ctx.thread_write_genovecs = nullptr;
      uint32_t write_mhc_needed = 0;
      if (new_sample_idx_to_old || subsetting_required) {
        if (unlikely(bigstack_alloc_wp(1, &ctx.thread_write_genovecs))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if (read_phase_present && new_sample_idx_to_old) {
          if (unlikely(bigstack_alloc_u32(raw_sample_ct, &ctx.old_sample_idx_to_new))) {
            goto MakePgenRobust_ret_NOMEM;
          }
          for (uint32_t new_sample_idx = 0; new_sample_idx != sample_ct; ++new_sample_idx) {
            ctx.old_sample_idx_to_new[new_sample_idx_to_old[new_sample_idx]] = new_sample_idx;
          }
        }
        if (unlikely(bigstack_alloc_w(sample_ctl2, &(ctx.thread_write_genovecs[0])))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        write_mhc_needed = 1;
      }
      ctx.thread_write_mhc = nullptr;
      if (write_mhc_needed) {
        if (unlikely(bigstack_alloc_wp(1, &ctx.thread_write_mhc))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        // todo: refalt1_select
        const uintptr_t mhcwrite_word_ct = GetMhcWordCt(sample_ct);
        if (unlikely(bigstack_alloc_w(mhcwrite_word_ct, &(ctx.thread_write_mhc[0])))) {
          goto MakePgenRobust_ret_NOMEM;
        }
      }
      ctx.thread_write_phasepresents = nullptr;
      ctx.thread_all_hets = nullptr;
      if (read_or_write_phase_present) {
        if (unlikely(
                bigstack_alloc_wp(1, &ctx.thread_write_phasepresents) ||
                bigstack_alloc_wp(1, &ctx.thread_write_phaseinfos) ||
                bigstack_alloc_w(sample_ctl, &(ctx.thread_write_phasepresents[0])) ||
                bigstack_alloc_w(sample_ctl, &(ctx.thread_write_phaseinfos[0])))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if (read_phase_present) {
          if (unlikely(
                  bigstack_alloc_wp(1, &ctx.thread_all_hets) ||
                  bigstack_alloc_w(raw_sample_ctl, &(ctx.thread_all_hets[0])))) {
            goto MakePgenRobust_ret_NOMEM;
          }
        }
      }
      ctx.thread_write_dosagepresents = nullptr;
      ctx.thread_write_dphasepresents = nullptr;
      if (read_or_write_dosage_present) {
        if (unlikely(
                bigstack_alloc_wp(1, &ctx.thread_write_dosagepresents) ||
                bigstack_alloc_dosagep(1, &ctx.thread_write_dosagevals) ||
                bigstack_alloc_w(sample_ctl, &(ctx.thread_write_dosagepresents[0])) ||
                bigstack_alloc_dosage(sample_ct, &(ctx.thread_write_dosagevals[0])))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if (read_or_write_dphase_present) {
          if (unlikely(
                  bigstack_alloc_wp(1, &ctx.thread_write_dphasepresents) ||
                  bigstack_alloc_dphasep(1, &ctx.thread_write_dphasedeltas) ||
                  bigstack_alloc_w(sample_ctl, &(ctx.thread_write_dphasepresents[0])) ||
                  bigstack_alloc_dphase(sample_ct + RoundUpPow2(sample_ct, kCacheline / 2), &(ctx.thread_write_dphasedeltas[0])))) {
            goto MakePgenRobust_ret_NOMEM;
          }
        }
      }
      if ((write_mhc_needed || read_dosage_present) && new_sample_idx_to_old) {
        if (unlikely(
                bigstack_alloc_u32p(1, &ctx.thread_cumulative_popcount_bufs) ||
                bigstack_alloc_u32(raw_sample_ctl, &(ctx.thread_cumulative_popcount_bufs[0])))) {
          goto MakePgenRobust_ret_NOMEM;
        }
      }
      mcp->refalt1_select = refalt1_select;
      if (refalt1_select) {
        if (write_allele_idx_offsets) {
          // this will require write_mhc and an additional AlleleCode buffer
          logerrputs("Error: Multiallelic allele rotation is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto MakePgenRobust_ret_1;
        }
        if (new_variant_idx_to_old || (variant_ct < raw_variant_ct)) {
          // might want inner loop to map variant uidx -> idx instead
          STD_ARRAY_PTR_DECL(AlleleCode, 2, tmp_refalt1_select);
          if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(AlleleCode, 2, variant_ct, &tmp_refalt1_select))) {
            goto MakePgenRobust_ret_NOMEM;
          }
          if (new_variant_idx_to_old) {
            for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
              const uintptr_t variant_uidx = new_variant_idx_to_old[variant_idx];
              STD_ARRAY_COPY(refalt1_select[variant_uidx], 2, tmp_refalt1_select[variant_idx]);
            }
          } else {
            uintptr_t variant_uidx_base = 0;
            uintptr_t cur_bits = variant_include[0];
            for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
              const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
              STD_ARRAY_COPY(refalt1_select[variant_uidx], 2, tmp_refalt1_select[variant_idx]);
            }
          }
          mcp->refalt1_select = tmp_refalt1_select;
        }
      }
      ctx.mcp = mcp;
      const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
      PgenVariant pgv;
      PreinitPgv(&pgv);
      uint32_t* alt_regular_one_cts = nullptr;
      uint32_t* alt_invphase_one_cts = nullptr;
      uint32_t* alt_two_cts = nullptr;
      uint32_t* alt_sample_idx_buf = nullptr;
      uint32_t** alt_regular_one_sample_idx_starts = nullptr;
      uint32_t** alt_invphase_one_sample_idx_starts = nullptr;
      uint32_t** alt_two_sample_idx_starts = nullptr;
      if (make_plink2_flags & (kfMakePlink2MSplitBase * 7)) {
        // split or join
        // this is currently for split with no dosages
        if (unlikely(
                bigstack_alloc_w(raw_sample_ctl2, &pgv.genovec) ||
                bigstack_alloc_w(raw_sample_ctl, &pgv.patch_01_set) ||
                bigstack_alloc_ac(raw_sample_ct, &pgv.patch_01_vals) ||
                bigstack_alloc_w(raw_sample_ctl, &pgv.patch_10_set) ||
                bigstack_alloc_ac(2 * raw_sample_ct, &pgv.patch_10_vals) ||
                bigstack_alloc_u32(max_read_allele_ct, &alt_regular_one_cts) ||
                bigstack_alloc_u32(max_read_allele_ct, &alt_two_cts) ||
                bigstack_alloc_u32(2 * raw_sample_ct + 1, &alt_sample_idx_buf) ||
                bigstack_alloc_u32p(max_read_allele_ct + 1, &alt_regular_one_sample_idx_starts) ||
                bigstack_alloc_u32p(max_read_allele_ct + 1, &alt_two_sample_idx_starts))) {
          goto MakePgenRobust_ret_NOMEM;
        }
        if (read_phase_present) {
          if (unlikely(
                  bigstack_alloc_w(raw_sample_ctl, &pgv.phasepresent) ||
                  bigstack_alloc_w(raw_sample_ctl, &pgv.phaseinfo) ||
                  bigstack_alloc_u32(max_read_allele_ct, &alt_invphase_one_cts) ||
                  bigstack_alloc_u32p(max_read_allele_ct + 1, &alt_invphase_one_sample_idx_starts))) {
            goto MakePgenRobust_ret_NOMEM;
          }
        }
        if (read_dosage_present) {
          logerrputs("Error: Multiallelic dosages aren't supported yet.\n");
          reterr = kPglRetNotYetSupported;
          goto MakePgenRobust_ret_1;
        }
      }

      const uint32_t raw_sample_ctv2 = NypCtToVecCt(raw_sample_ct);
      uintptr_t load_variant_vec_ct = raw_sample_ctv2;
      uint32_t loaded_vrtypes_needed = (read_gflags & kfPgenGlobalMultiallelicHardcallFound)? 1 : 0;
      if (read_phase_present || read_dosage_present) {
        loaded_vrtypes_needed = 1;
        if (read_phase_present) {
          // phaseraw has three parts:
          // 1. het_ct as uint32_t, and explicit_phasepresent_ct as uint32_t.
          // 2. vec-aligned bitarray of up to (raw_sample_ct + 1) bits.  first
          //    bit is set iff phasepresent is explicitly stored at all (if
          //    not, all hets are assumed to be phased), if yes the remaining
          //    bits store packed phasepresent values for all hets, if no the
          //    remaining bits store packed phaseinfo values for all hets.
          // 3. word-aligned bitarray of up to raw_sample_ct bits, storing
          //    phaseinfo values.  (end of this array is vec-aligned.)
          const uintptr_t phaseraw_word_ct = (8 / kBytesPerWord) + kWordsPerVec + RoundDownPow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
          load_variant_vec_ct += WordCtToVecCt(phaseraw_word_ct);
        }
        if (read_dosage_present) {
          // biallelic dosageraw has two parts:
          // 1. vec-aligned bitarray of up to raw_sample_ct bits, storing which
          //    samples have dosages.
          // 2. word-aligned array of uint16s with 0..32768 fixed-point
          //    dosages.
          // dphaseraw has the same structure, with the uint16s replaced with
          // an int16 array of (left - right) values.
          const uintptr_t dosageraw_word_ct = kWordsPerVec * (BitCtToVecCt(raw_sample_ct) + DivUp(raw_sample_ct, (kBytesPerVec / sizeof(Dosage))));
          load_variant_vec_ct += WordCtToVecCt(dosageraw_word_ct) * (1 + read_dphase_present);
        }
      }

      uintptr_t bytes_left = bigstack_left();
      if (unlikely(bytes_left < 7 * kCacheline)) {
        goto MakePgenRobust_ret_NOMEM;
      }
      bytes_left -= 7 * kCacheline;  // defend against adverse rounding
      uintptr_t ulii = bytes_left / (2 * (kBytesPerVec * load_variant_vec_ct + loaded_vrtypes_needed));
      if (unlikely(!ulii)) {
        goto MakePgenRobust_ret_NOMEM;
      }
      if (ulii > MINV(kPglVblockSize, write_variant_ct)) {
        ulii = MINV(kPglVblockSize, write_variant_ct);
      }
      const uint32_t write_block_size = ulii;
      uintptr_t* main_loadbufs[2];
      main_loadbufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(load_variant_vec_ct * kBytesPerVec * write_block_size));
      main_loadbufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(load_variant_vec_ct * kBytesPerVec * write_block_size));

      // todo: multiallelic trim-alts support

      if (loaded_vrtypes_needed) {
        ctx.loaded_vrtypes[0] = S_CAST(unsigned char*, bigstack_alloc_raw_rd(write_block_size));
        ctx.loaded_vrtypes[1] = S_CAST(unsigned char*, bigstack_alloc_raw_rd(write_block_size));
      } else {
        ctx.loaded_vrtypes[0] = nullptr;
        ctx.loaded_vrtypes[1] = nullptr;
      }
      SetThreadFuncAndData(MakePgenThread, &ctx, &tg);

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
      const uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
      const uintptr_t* cur_write_allele_idx_offsets = nullptr;
      const uint32_t batch_ct_m1 = (write_variant_ct - 1) / write_block_size;
      uint32_t pct = 0;
      uint32_t parity = 0;
      uint32_t cur_batch_size = write_block_size;
      uint32_t next_print_write_variant_idx = write_variant_ct / 100;
      uint32_t cur_read_allele_ct = 2;
      uint32_t cur_write_allele_ct = 2;
      uint32_t cur_het_ct = 0;
      uintptr_t read_variant_uidx_base = 0;

      // now need to retain these across loop iterations in case a split is
      // interrupted by batch-end
      uint32_t read_variant_uidx = 0;
      uint32_t write_aidx = 1;

      uintptr_t cur_bits = variant_include[0];
      PgrSampleSubsetIndex null_pssi;
      PgrClearSampleSubsetIndex(simple_pgrp, &null_pssi);
      for (uint32_t read_batch_idx = 0; ; ++read_batch_idx) {
        if (!IsLastBlock(&tg)) {
          if (read_batch_idx == batch_ct_m1) {
            cur_batch_size = write_variant_ct - (read_batch_idx * write_block_size);
          }
          uintptr_t* cur_loadbuf = main_loadbufs[parity];
          uintptr_t* loadbuf_iter = cur_loadbuf;
          unsigned char* cur_loaded_vrtypes = ctx.loaded_vrtypes[parity];
          ctx.loadbuf_thread_starts[parity][0] = loadbuf_iter;
          if (write_allele_idx_offsets) {
            cur_write_allele_idx_offsets = &(write_allele_idx_offsets[read_batch_idx * write_block_size]);
          }
          for (uint32_t block_widx = 0; block_widx != cur_batch_size; ) {
            if (write_aidx == 1) {
              if (!new_variant_idx_to_old_iter) {
                read_variant_uidx = BitIter1(variant_include, &read_variant_uidx_base, &cur_bits);
              } else {
                read_variant_uidx = *new_variant_idx_to_old_iter++;
              }
              // todo: multiallelic trim-alts
              // todo: multiallelic merge
              // split: load to buffer instead of loadbuf_iter, have function
              //        for writing to loadbuf_iter given buffer contents, this
              //        should work if split is 'interrupted' by batch boundary
              //        in middle
              // merge: track loadbuf_iter location at beginning of
              //        same-position block... (finish writing this later)
              if (allele_idx_offsets) {
                cur_read_allele_ct = allele_idx_offsets[read_variant_uidx + 1] - allele_idx_offsets[read_variant_uidx];
              }
            }
            if (cur_write_allele_idx_offsets) {
              cur_write_allele_ct = cur_write_allele_idx_offsets[block_widx + 1] - cur_write_allele_idx_offsets[block_widx];
            }
            if (cur_read_allele_ct == cur_write_allele_ct) {
              reterr = PgrGetRaw(read_variant_uidx, read_gflags, simple_pgrp, &loadbuf_iter, cur_loaded_vrtypes? (&(cur_loaded_vrtypes[block_widx])) : nullptr);
              if (unlikely(reterr)) {
                goto MakePgenRobust_ret_PGR_FAIL;
              }
              ++block_widx;
              continue;
            } else if (cur_write_allele_ct == 2) {
              if (write_aidx == 1) {
                // 1. read into normal, not raw representation
                if (read_phase_present) {
                  reterr = PgrGetMDp(nullptr, null_pssi, raw_sample_ct, read_variant_uidx, simple_pgrp, &pgv);
                } else {
                  reterr = PgrGetMD(nullptr, null_pssi, raw_sample_ct, read_variant_uidx, simple_pgrp, &pgv);
                }
                if (unlikely(reterr)) {
                  goto MakePgenRobust_ret_PGR_FAIL;
                }

                // 2a. count # of each alt
                // 2b. create het and hom lists for each alt
                uintptr_t* genovec = pgv.genovec;
                ZeroTrailingNyps(raw_sample_ct, genovec);
                uint32_t raw_01_ct;
                uint32_t raw_10_ct;
                GenovecCount12Unsafe(genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
                ZeroU32Arr(cur_read_allele_ct, alt_regular_one_cts);
                alt_regular_one_cts[1] = raw_01_ct - pgv.patch_01_ct;
                for (uint32_t rarealt_idx = 0; rarealt_idx != pgv.patch_01_ct; ++rarealt_idx) {
                  alt_regular_one_cts[pgv.patch_01_vals[rarealt_idx]] += 1;
                }
                ZeroU32Arr(cur_read_allele_ct, alt_two_cts);
                if (!pgv.phasepresent_ct) {
                  for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                    const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                    const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                    if (ac0 == ac1) {
                      alt_two_cts[ac0] += 1;
                    } else {
                      alt_regular_one_cts[ac0] += 1;
                      alt_regular_one_cts[ac1] += 1;
                    }
                  }
                } else {
                  ZeroU32Arr(cur_read_allele_ct, alt_invphase_one_cts);
                  for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                    const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                    const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                    if (ac0 == ac1) {
                      alt_two_cts[ac0] += 1;
                    } else {
                      alt_invphase_one_cts[ac0] += 1;
                      alt_regular_one_cts[ac1] += 1;
                    }
                  }
                }

                alt_two_cts[1] = raw_10_ct - pgv.patch_10_ct;
                cur_het_ct = raw_01_ct + pgv.patch_10_ct;
                for (uint32_t aidx = 2; aidx != cur_read_allele_ct; ++aidx) {
                  cur_het_ct -= alt_two_cts[aidx];
                }

                uint32_t* sample_idx_buf_iter = alt_sample_idx_buf;
                alt_regular_one_sample_idx_starts[0] = alt_sample_idx_buf;
                for (uint32_t aidx = 1; aidx != cur_read_allele_ct; ++aidx) {
                  alt_regular_one_sample_idx_starts[aidx] = sample_idx_buf_iter;
                  sample_idx_buf_iter = &(sample_idx_buf_iter[alt_regular_one_cts[aidx]]);
                }
                alt_regular_one_sample_idx_starts[cur_read_allele_ct] = sample_idx_buf_iter;
                if (pgv.phasepresent_ct) {
                  alt_invphase_one_sample_idx_starts[0] = sample_idx_buf_iter;
                  for (uint32_t aidx = 1; aidx != cur_read_allele_ct - 1; ++aidx) {
                    alt_invphase_one_sample_idx_starts[aidx] = sample_idx_buf_iter;
                    sample_idx_buf_iter = &(sample_idx_buf_iter[alt_invphase_one_cts[aidx]]);
                  }
                  alt_invphase_one_sample_idx_starts[cur_read_allele_ct - 1] = sample_idx_buf_iter;
                  alt_invphase_one_sample_idx_starts[cur_read_allele_ct] = sample_idx_buf_iter;
                }
                alt_two_sample_idx_starts[0] = sample_idx_buf_iter;
                for (uint32_t aidx = 1; aidx != cur_read_allele_ct; ++aidx) {
                  alt_two_sample_idx_starts[aidx] = sample_idx_buf_iter;
                  sample_idx_buf_iter = &(sample_idx_buf_iter[alt_two_cts[aidx]]);
                }
                alt_two_sample_idx_starts[cur_read_allele_ct] = sample_idx_buf_iter;

                Halfword* patch_01_set_alias = R_CAST(Halfword*, pgv.patch_01_set);
                Halfword* patch_10_set_alias = R_CAST(Halfword*, pgv.patch_10_set);
                uint32_t idx_01 = 0;
                uint32_t idx_10 = 0;
                for (uint32_t widx = 0; widx != raw_sample_ctl2; ++widx) {
                  const uintptr_t geno_word = genovec[widx];
                  const uint32_t sample_idx_offset = widx * kBitsPerWordD2;
                  uintptr_t geno_01 = Word01(geno_word);
                  if (geno_01) {
                    if (!pgv.patch_01_ct) {
                      // patch_01_set not initialized in this case
                      do {
                        const uint32_t sample_idx = sample_idx_offset + ctzw(geno_01) / 2;
                        alt_regular_one_sample_idx_starts[1][0] = sample_idx;
                        alt_regular_one_sample_idx_starts[1] += 1;
                        geno_01 &= geno_01 - 1;
                      } while (geno_01);
                    } else {
                      uint32_t geno_01_hw = PackWordToHalfword(geno_01);
                      const uint32_t patch_01_hw = patch_01_set_alias[widx];
                      do {
                        const uint32_t lowbit = geno_01_hw & (-geno_01_hw);
                        const uint32_t sample_idx = sample_idx_offset + ctzu32(lowbit);
                        if (lowbit & patch_01_hw) {
                          AlleleCode ac = pgv.patch_01_vals[idx_01];
                          alt_regular_one_sample_idx_starts[ac][0] = sample_idx;
                          alt_regular_one_sample_idx_starts[ac] += 1;
                          ++idx_01;
                        } else {
                          alt_regular_one_sample_idx_starts[1][0] = sample_idx;
                          alt_regular_one_sample_idx_starts[1] += 1;
                        }
                        geno_01_hw ^= lowbit;
                      } while (geno_01_hw);
                    }
                  }
                  uintptr_t geno_10 = Word10(geno_word);
                  if (geno_10) {
                    if (!pgv.patch_10_ct) {
                      // patch_10_set not initialized in this case
                      do {
                        const uint32_t sample_idx = sample_idx_offset + ctzw(geno_10) / 2;
                        alt_two_sample_idx_starts[1][0] = sample_idx;
                        alt_two_sample_idx_starts[1] += 1;
                        geno_10 &= geno_10 - 1;
                      } while (geno_10);
                    } else {
                      uint32_t geno_10_hw = PackWordToHalfword(geno_10);
                      const uint32_t patch_10_hw = patch_10_set_alias[widx];
                      if (!pgv.phasepresent_ct) {
                        do {
                          const uint32_t lowbit = geno_10_hw & (-geno_10_hw);
                          const uint32_t sample_idx = sample_idx_offset + ctzu32(lowbit);
                          if (lowbit & patch_10_hw) {
                            AlleleCode ac0 = pgv.patch_10_vals[2 * idx_10];
                            AlleleCode ac1 = pgv.patch_10_vals[2 * idx_10 + 1];
                            if (ac0 == ac1) {
                              alt_two_sample_idx_starts[ac0][0] = sample_idx;
                              alt_two_sample_idx_starts[ac0] += 1;
                            } else {
                              alt_regular_one_sample_idx_starts[ac0][0] = sample_idx;
                              alt_regular_one_sample_idx_starts[ac0] += 1;
                              alt_regular_one_sample_idx_starts[ac1][0] = sample_idx;
                              alt_regular_one_sample_idx_starts[ac1] += 1;
                            }
                            ++idx_10;
                          } else {
                            alt_two_sample_idx_starts[1][0] = sample_idx;
                            alt_two_sample_idx_starts[1] += 1;
                          }
                          geno_10_hw ^= lowbit;
                        } while (geno_10_hw);
                      } else {
                        do {
                          const uint32_t lowbit = geno_10_hw & (-geno_10_hw);
                          const uint32_t sample_idx = sample_idx_offset + ctzu32(lowbit);
                          if (lowbit & patch_10_hw) {
                            AlleleCode ac0 = pgv.patch_10_vals[2 * idx_10];
                            AlleleCode ac1 = pgv.patch_10_vals[2 * idx_10 + 1];
                            if (ac0 == ac1) {
                              alt_two_sample_idx_starts[ac0][0] = sample_idx;
                              alt_two_sample_idx_starts[ac0] += 1;
                            } else {
                              alt_invphase_one_sample_idx_starts[ac0][0] = sample_idx;
                              alt_invphase_one_sample_idx_starts[ac0] += 1;
                              alt_regular_one_sample_idx_starts[ac1][0] = sample_idx;
                              alt_regular_one_sample_idx_starts[ac1] += 1;
                            }
                            ++idx_10;
                          } else {
                            alt_two_sample_idx_starts[1][0] = sample_idx;
                            alt_two_sample_idx_starts[1] += 1;
                          }
                          geno_10_hw ^= lowbit;
                        } while (geno_10_hw);
                      }
                    }
                  }
                }
                for (uint32_t aidx = cur_read_allele_ct - 1; aidx; --aidx) {
                  alt_regular_one_sample_idx_starts[aidx] = alt_regular_one_sample_idx_starts[aidx - 1];
                  alt_two_sample_idx_starts[aidx] = alt_two_sample_idx_starts[aidx - 1];
                }
                if (pgv.phasepresent_ct) {
                  for (uint32_t aidx = cur_read_allele_ct - 1; aidx; --aidx) {
                    alt_invphase_one_sample_idx_starts[aidx] = alt_invphase_one_sample_idx_starts[aidx - 1];
                  }
                }
                // todo: multiallelic dosage

                for (uint32_t widx = 0; widx != raw_sample_ctl2; ++widx) {
                  // keep 3s, set 1s and 2s to 0
                  genovec[widx] = Word11(genovec[widx]) * 3;
                }
              }
              const uint32_t split_stop = MINV(cur_batch_size + 1 - block_widx, cur_read_allele_ct);
              for (; write_aidx != split_stop; ++write_aidx, ++block_widx) {
                // 3. synthesize raw
                //   (save to loaded_vrtypes if necessary)
                //   genovec, vector-aligned
                //   if hphase present and relevant:
                //     (compute het_ct; het_ctdl := het_ct / kBitsPerWord)
                //     (first_half_byte_ct := 1 + (het_ct / CHAR_BIT))
                //     <uint32 het_ct>
                //     <uint32 raw_phasepresent_ct if explicit>
                //     <first_half_byte_ct phasepresent or phaseinfo bytes>
                //     <0-pad up to word boundary, to make popcount safe>
                //     [if explicit phasepresent, i.e. lowest bit set:
                //       (second_half_byte_ct := DivUp(raw_phasepresent_ct, 8))
                //       <second_half_byte_ct phaseinfo contents>
                //     ]
                //     align up to vector boundary
                uintptr_t* new_genovec = loadbuf_iter;
                memcpy(new_genovec, pgv.genovec, raw_sample_ctl2 * sizeof(intptr_t));
                loadbuf_iter = &(loadbuf_iter[raw_sample_ctv2 * kWordsPerVec]);
                uint32_t new_phasepresent_ct = 0;
                uint32_t new_het_ct = 0;
                uint32_t* regular_stop = alt_regular_one_sample_idx_starts[write_aidx + 1];
                if (pgv.phasepresent_ct) {
                  uint32_t* regular_iter = alt_regular_one_sample_idx_starts[write_aidx];
                  uint32_t* invphase_iter = alt_invphase_one_sample_idx_starts[write_aidx];
                  uint32_t* invphase_stop = alt_invphase_one_sample_idx_starts[write_aidx + 1];
                  new_het_ct = (regular_stop - regular_iter) + (invphase_stop - invphase_iter);
                  if (pgv.phasepresent_ct == cur_het_ct) {
                    new_phasepresent_ct = new_het_ct;
                  } else {
                    uintptr_t* phasepresent = pgv.phasepresent;
                    for (; regular_iter != regular_stop; ++regular_iter) {
                      new_phasepresent_ct += IsSet(phasepresent, *regular_iter);
                    }
                    for (; invphase_iter != invphase_stop; ++invphase_iter) {
                      new_phasepresent_ct += IsSet(phasepresent, *invphase_iter);
                    }
                  }
                }
                uint32_t* two_stop = alt_two_sample_idx_starts[write_aidx + 1];
                for (uint32_t* two_iter = alt_two_sample_idx_starts[write_aidx]; two_iter != two_stop; ++two_iter) {
                  const uint32_t sample_uidx = *two_iter;
                  SetBit(sample_uidx * 2 + 1, new_genovec);
                }
                uint32_t* regular_iter = alt_regular_one_sample_idx_starts[write_aidx];
                if (!new_phasepresent_ct) {
                  for (; regular_iter != regular_stop; ++regular_iter) {
                    const uint32_t sample_uidx = *regular_iter;
                    SetBit(sample_uidx * 2, new_genovec);
                  }
                  if (pgv.phasepresent_ct) {
                    uint32_t* invphase_stop = alt_invphase_one_sample_idx_starts[write_aidx + 1];
                    for (uint32_t* invphase_iter = alt_invphase_one_sample_idx_starts[write_aidx]; invphase_iter != invphase_stop; ++invphase_iter) {
                      const uint32_t sample_uidx = *invphase_iter;
                      SetBit(sample_uidx * 2, new_genovec);
                    }
                  }
                  if (cur_loaded_vrtypes) {
                    cur_loaded_vrtypes[block_widx] = 0;
                  }
                } else {
                  // need to write raw hphase
                  const uint32_t het_ctdl = new_het_ct / kBitsPerWord;
                  uintptr_t* shifted_part1 = &(loadbuf_iter[8 / kBytesPerWord]);
                  uintptr_t* part1_end = &(shifted_part1[1 + het_ctdl]);
                  uint32_t* invphase_iter = alt_invphase_one_sample_idx_starts[write_aidx];
                  uint32_t* invphase_stop = alt_invphase_one_sample_idx_starts[write_aidx + 1];
                  const uint32_t orig_regular_end = *regular_stop;
                  const uint32_t orig_invphase_end = *invphase_stop;
                  // sentinel value to simplify the next loop.
                  *invphase_stop = UINT32_MAX;
                  // must grab this before setting *regular_stop, in case
                  // they overlap; and after setting *invphase_stop, in case
                  // this list is empty
                  uint32_t invphase_idx = *invphase_iter++;

                  *regular_stop = UINT32_MAX;
                  uint32_t regular_idx = *regular_iter++;
                  uint32_t shifted_het_idx = 1;
                  if (new_phasepresent_ct == new_het_ct) {
                    loadbuf_iter[0] = new_het_ct;
#ifndef __LP64__
                    loadbuf_iter[1] = 0;
#endif
                    // shifted_part1 is phaseinfo
                    shifted_part1[0] = 0;
                    shifted_part1[het_ctdl] = 0;
                    while (regular_idx != invphase_idx) {
                      uintptr_t is_inverted = (invphase_idx < regular_idx);
                      uint32_t sample_uidx;
                      if (is_inverted) {
                        sample_uidx = invphase_idx;
                        invphase_idx = *invphase_iter++;
                      } else {
                        sample_uidx = regular_idx;
                        regular_idx = *regular_iter++;
                      }
                      SetBit(sample_uidx * 2, new_genovec);
                      AssignBit(shifted_het_idx, is_inverted ^ IsSet(pgv.phaseinfo, sample_uidx), shifted_part1);
                      ++shifted_het_idx;
                    }
                    assert(shifted_het_idx == new_het_ct + 1);
                    loadbuf_iter = part1_end;
                  } else {
#ifdef __LP64__
                    loadbuf_iter[0] = new_het_ct | (S_CAST(uint64_t, new_phasepresent_ct) << 32);
#else
                    loadbuf_iter[0] = new_het_ct;
                    loadbuf_iter[1] = new_phasepresent_ct;
#endif
                    shifted_part1[0] = 1;
                    memset(shifted_part1, 0, (1 + het_ctdl) * sizeof(intptr_t));
                    // shifted_part1 is phasepresent
                    // part1_end is start of phaseinfo
                    const uint32_t new_phasepresent_ctl = BitCtToWordCt(new_phasepresent_ct);
                    part1_end[new_phasepresent_ctl - 1] = 0;
                    uint32_t phasepresent_idx = 0;
                    while (regular_idx != invphase_idx) {
                      uintptr_t is_inverted = (invphase_idx < regular_idx);
                      uint32_t sample_uidx;
                      if (is_inverted) {
                        sample_uidx = invphase_idx;
                        invphase_idx = *invphase_iter++;
                      } else {
                        sample_uidx = regular_idx;
                        regular_idx = *regular_iter++;
                      }
                      SetBit(sample_uidx * 2, new_genovec);
                      if (IsSet(pgv.phasepresent, sample_uidx)) {
                        SetBit(shifted_het_idx, shifted_part1);
                        AssignBit(phasepresent_idx, is_inverted ^ IsSet(pgv.phaseinfo, sample_uidx), part1_end);
                        ++phasepresent_idx;
                      }
                      ++shifted_het_idx;
                    }
                    assert(phasepresent_idx == new_phasepresent_ct);
                  }
                  assert(regular_idx == UINT32_MAX);
                  *regular_stop = orig_regular_end;
                  *invphase_stop = orig_invphase_end;
                  VecAlignUp(&loadbuf_iter);
                  if (cur_loaded_vrtypes) {
                    cur_loaded_vrtypes[block_widx] = 0x10;
                  }
                }
              }
              if (split_stop != cur_read_allele_ct) {
                break;
              }
              write_aidx = 1;
            } else {
              // merge; todo
            }
          }
        }
        if (read_batch_idx) {
          JoinThreads(&tg);
          reterr = ctx.write_reterr;
          if (unlikely(reterr)) {
            if (reterr == kPglRetWriteFail) {
              errno = ctx.write_errno;
            }
            goto MakePgenRobust_ret_1;
          }
        }
        if (!IsLastBlock(&tg)) {
          ctx.cur_block_write_ct = cur_batch_size;
          if (read_batch_idx == batch_ct_m1) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto MakePgenRobust_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (read_batch_idx) {
          if (read_batch_idx > batch_ct_m1) {
            break;
          }
          const uint32_t write_idx_end = read_batch_idx * write_block_size;
          if (write_idx_end >= next_print_write_variant_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (write_idx_end * 100LLU) / write_variant_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_write_variant_idx = (pct * S_CAST(uint64_t, write_variant_ct)) / 100;
          }
        }
      }
      SpgwFinish(ctx.spgwp);
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logputs("done.\n");
    }
  }
  while (0) {
  MakePgenRobust_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakePgenRobust_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  MakePgenRobust_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 MakePgenRobust_ret_1:
  CleanupThreads(&tg);
  CleanupSpgw(&spgw, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// allele_presents should be nullptr iff trim_alts not true
PglErr MakePlink2NoVsort(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const char* varid_template_str, __maybe_unused const char* varid_multi_template_str, __maybe_unused const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, UnsortedVar vpos_sortstatus, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t new_variant_id_max_allele_slen, MiscFlags misc_flags, MakePlink2Flags make_plink2_flags, PvarPsamFlags pvar_psam_flags, uintptr_t pgr_alloc_cacheline_ct, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  MTPgenWriter* mpgwp = nullptr;
  MakePgenCtx ctx;
  {
    if (make_plink2_flags & kfMakeFam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, new_sample_idx_to_old, sample_ct, pheno_ct, '\t');
      if (unlikely(reterr)) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePsam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WritePsam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, sample_ct, pheno_ct, max_pheno_name_blen, pvar_psam_flags);
      if (unlikely(reterr)) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
    const uint32_t input_biallelic = (!allele_idx_offsets);
    // output_biallelic: test write_allele_idx_offsets equality to null
    PgenGlobalFlags read_gflags = pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
    if (!input_biallelic) {
      // Can only skip this when there are actually zero copies of alt2+.
      // Otherwise, even with erase-alt2+, we still need to distinguish alt1
      // from alt2 so we can set calls involving the latter to missing.
      read_gflags |= kfPgenGlobalMultiallelicHardcallFound;
    }
    const uintptr_t* write_allele_idx_offsets = nullptr;
    uint32_t write_variant_ct = variant_ct;
    uint32_t max_write_allele_ct = max_allele_ct;
    uint32_t max_missalt_ct = 0;
    if (make_plink2_flags & kfMakePlink2MMask) {
      // TODO: enforce on command line
      assert((!refalt1_select) && (!allele_presents));
      if (make_plink2_flags & kfMakePlink2MJoin) {
        reterr = PlanMultiallelicJoin(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, make_plink2_flags, &write_variant_ct, &write_allele_idx_offsets, &max_write_allele_ct, &max_missalt_ct);
      } else if (!allele_idx_offsets) {
        // no splitting to do
        logputs("Note: All variants are biallelic; nothing to split.\n");
      } else {
        reterr = PlanMultiallelicSplit(variant_include, allele_idx_offsets, allele_storage, max_allele_ct, make_plink2_flags, &write_variant_ct, &write_allele_idx_offsets);
      }
      if (unlikely(reterr)) {
        goto MakePlink2NoVsort_ret_1;
      }
    } else if (allele_idx_offsets) {
      if (allele_presents) {
        fputs("multiallelic variants + trim-alts not yet supported\n", stderr);
        exit(S_CAST(int32_t, kPglRetNotYetSupported));
      }
      if (variant_ct < raw_variant_ct) {
        uintptr_t* new_allele_idx_offsets;
        if (bigstack_alloc_w(variant_ct + 1, &new_allele_idx_offsets)) {
          goto MakePlink2NoVsort_ret_NOMEM;
        }
        const uintptr_t final_offset = InitWriteAlleleIdxOffsets(variant_include, allele_idx_offsets, nullptr, refalt1_select, nullptr, variant_ct, new_allele_idx_offsets);
        if (final_offset != 2 * variant_ct) {
          new_allele_idx_offsets[variant_ct] = final_offset;
          write_allele_idx_offsets = new_allele_idx_offsets;
        } else {
          BigstackReset(new_allele_idx_offsets);
        }
      } else {
        write_allele_idx_offsets = allele_idx_offsets;
      }
    }
    if (make_plink2_flags & kfMakeBim) {
      const uint32_t bim_zst = (make_plink2_flags / kfMakeBimZs) & 1;
      OutnameZstSet(".bim", bim_zst, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      if (unlikely(write_allele_idx_offsets)) {
        logputs("\n");
        logerrprintf("Error: %s cannot contain multiallelic variants.\n", outname);
        goto MakePlink2NoVsort_ret_INCONSISTENT_INPUT;
      }
      if (write_variant_ct == variant_ct) {
        reterr = WriteMapOrBim(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_presents, refalt1_select, variant_cms, variant_ct, max_allele_slen, '\t', bim_zst, max_thread_ct);
      } else {
        assert(write_variant_ct > variant_ct);
        reterr = WriteBimSplit(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, variant_cms, varid_template_str, missing_varid_match, variant_ct, max_allele_slen, new_variant_id_max_allele_slen, (make_plink2_flags / kfMakePlink2VaridSemicolon) & 1, (make_plink2_flags / kfMakePlink2VaridDup) & 1, misc_flags, bim_zst, max_thread_ct);
      }
      if (unlikely(reterr)) {
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
      if (write_variant_ct == variant_ct) {
        reterr = WritePvar(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_presents, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, info_flags, nonref_flags_storage, max_filter_slen, info_reload_slen, vpos_sortstatus, pvar_psam_flags, max_thread_ct, xheader);
      } else {
        const char* const* info_keys = nullptr;
        uint32_t info_key_ct = 0;
        uint32_t* info_keys_htable = nullptr;
        uint32_t info_keys_htable_size = 0;
        if (pvar_info_reload) {
          reterr = ParseInfoHeader(xheader, xheader_blen, &info_keys, &info_key_ct, &info_keys_htable, &info_keys_htable_size);
          if (reterr) {
            goto MakePlink2NoVsort_ret_1;
          }
        }
        if (write_variant_ct > variant_ct) {
          reterr = WritePvarSplit(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, varid_template_str, missing_varid_match, info_keys, info_keys_htable, raw_variant_ct, variant_ct, max_allele_slen, new_variant_id_max_allele_slen, xheader_blen, info_flags, nonref_flags_storage, max_filter_slen, info_reload_slen, vpos_sortstatus, info_key_ct, info_keys_htable_size, misc_flags, make_plink2_flags, pvar_psam_flags, max_thread_ct, xheader);
        } else {
          logerrputs("Error: Multiallelic join is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto MakePlink2NoVsort_ret_1;
          // reterr = WritePvarJoin(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pgfip->nonref_flags, pvar_info_reload, variant_cms, varid_template_str, missing_varid_match, info_keys, info_keys_htable, raw_variant_ct, variant_ct, max_allele_slen, new_variant_id_max_allele_slen, max_write_allele_ct, max_missalt_ct, xheader_blen, info_flags, nonref_flags_storage, max_filter_slen, info_reload_slen, vpos_sortstatus, info_key_ct, info_keys_htable_size, misc_flags, make_plink2_flags, pvar_psam_flags, max_thread_ct, xheader);
        }
      }
      if (unlikely(reterr)) {
        goto MakePlink2NoVsort_ret_1;
      }
      logputs("done.\n");
    }
    MakeCommon mc;
    mc.plink2_write_flags = kfPlink2Write0;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    ctx.sex_male_collapsed = nullptr;  // defensive
    if (make_plink2_flags & kfMakePlink2SetHhMissing) {
      const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
      uintptr_t* new_sex_male;
      if (unlikely(
              bigstack_alloc_w(sample_ctv * kWordsPerVec, &new_sex_male) ||
              bigstack_alloc_w(sample_ctv * kWordsPerVec, &mc.sex_male_collapsed_interleaved) ||
              bigstack_alloc_w(sample_ctv * kWordsPerVec, &ctx.sex_female_collapsed) ||
              bigstack_alloc_w(sample_ctv * kWordsPerVec, &mc.sex_female_collapsed_interleaved))) {
        goto MakePlink2NoVsort_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, sample_include, sample_ct, new_sex_male);
      ZeroTrailingWords(sample_ctl, new_sex_male);
      ctx.sex_male_collapsed = new_sex_male;
      FillInterleavedMaskVec(ctx.sex_male_collapsed, sample_ctv, mc.sex_male_collapsed_interleaved);

      CopyBitarrSubset(sex_nm, sample_include, sample_ct, ctx.sex_female_collapsed);
      BitvecInvmask(new_sex_male, sample_ctl, ctx.sex_female_collapsed);
      ZeroTrailingWords(sample_ctl, ctx.sex_female_collapsed);
      FillInterleavedMaskVec(ctx.sex_female_collapsed, sample_ctv, mc.sex_female_collapsed_interleaved);

      mc.plink2_write_flags |= kfPlink2WriteSetHhMissing;
      if (make_plink2_flags & kfMakePlink2SetHhMissingKeepDosage) {
        mc.plink2_write_flags |= kfPlink2WriteSetHhMissingKeepDosage;
      }
    } else {
      // defensive
      mc.sex_male_collapsed_interleaved = nullptr;
      mc.sex_female_collapsed_interleaved = nullptr;
    }
    if (make_plink2_flags & kfMakePlink2SetMixedMtMissing) {
      mc.plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
      if (make_plink2_flags & kfMakePlink2SetMixedMtMissingKeepDosage) {
        mc.plink2_write_flags |= kfPlink2WriteSetMixedMtMissingKeepDosage;
      }
    }
    mc.cip = cip;
    mc.raw_sample_ct = raw_sample_ct;
    mc.sample_ct = sample_ct;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    const uint32_t make_pgen = make_plink2_flags & kfMakePgen;
    // todo: prohibit .pgen + .bim write when data is multiallelic without
    //   either multiallelic split or erase-alt2+ specified
    //   (--make-bed = automatic erase-alt2+?)
    if ((make_plink2_flags & kfMakeBed) || ((make_plink2_flags & (kfMakePgen | (kfMakePgenFormatBase * 3))) == (kfMakePgen | kfMakePgenFormatBase))) {
      reterr = MakeBedlikeMain(sample_include, new_sample_idx_to_old, variant_include, refalt1_select, raw_variant_ct, variant_ct, max_thread_ct, hard_call_thresh, make_plink2_flags, pgr_alloc_cacheline_ct, pgfip, &mc, outname, outname_end);
    } else if (make_pgen) {
      assert(variant_ct);
      assert(sample_ct);
      if (make_plink2_flags & (kfMakePlink2MSplitBase * 7)) {
        // don't duplicate complicated multiallelic split/merge/trim-alts logic
        // here for now.
        // (also auto-punt multiallelic dosage?)
        goto MakePlink2NoVsort_fallback;
      }
      ctx.write_allele_idx_offsets = write_allele_idx_offsets;
      if (variant_ct == raw_variant_ct) {
        ctx.write_chr_fo_vidx_start = cip->chr_fo_vidx_start;
      } else {
        if (AllocAndFillSubsetChrFoVidxStart(variant_include, cip, &ctx.write_chr_fo_vidx_start)) {
          goto MakePlink2NoVsort_fallback;
        }
      }
      if (make_plink2_flags & kfMakePgenErasePhase) {
        read_gflags &= ~(kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
      }
      const uint32_t fill_missing_from_dosage = (make_plink2_flags & kfMakePgenFillMissingFromDosage) && (read_gflags & kfPgenGlobalDosagePresent);
      if (fill_missing_from_dosage) {
        mc.plink2_write_flags |= kfPlink2WriteFillMissingFromDosage;
      }
      if (make_plink2_flags & kfMakePgenEraseDosage) {
        if ((hard_call_thresh == UINT32_MAX) && (!fill_missing_from_dosage)) {
          read_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
        } else {
          // erase-dosage + --hard-call-threshold currently requires dosages to
          // be read, and only thrown away at the last minute
          // (alternatively, we could build --hard-call-threshold directly into
          // pgr_read_raw?)
          mc.plink2_write_flags |= kfPlink2WriteLateDosageErase;
        }
      }
      if (read_gflags && (variant_ct < raw_variant_ct)) {
        // did we e.g. filter out all the phased variants?
        // do not check for multiallelic-hc here for now
        // (write_allele_idx_offsets check above serves the same purpose)
        read_gflags &= kfPgenGlobalMultiallelicHardcallFound | GflagsVfilter(variant_include, pgfip->vrtypes, raw_variant_ct, pgfip->gflags);
      }
      // could check if all the phased samples were also filtered out, but
      // that's already caught by running --make-pgen twice, so not a big deal

      const uint32_t read_dosage_present = (read_gflags / kfPgenGlobalDosagePresent) & 1;
      mc.hard_call_halfdist = ((hard_call_thresh == UINT32_MAX) || (!read_dosage_present))? 0 : (kDosage4th - hard_call_thresh);
      ctx.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
      const uint32_t read_phase_present = !!(read_gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent));
      const uint32_t read_dphase_present = (read_gflags / kfPgenGlobalDosagePhasePresent) & 1;
      PgenGlobalFlags write_gflags = read_gflags;
      uint32_t read_or_write_phase_present = read_phase_present;
      uint32_t read_or_write_dphase_present = read_dphase_present;
      if ((mc.hard_call_halfdist || fill_missing_from_dosage) && (read_phase_present || read_or_write_dphase_present)) {
        read_or_write_phase_present = 1;
        read_or_write_dphase_present = 1;
        write_gflags |= kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent;
      } else if (dosage_erase_thresh && read_dosage_present) {
        read_or_write_phase_present = 1;
      }
      uint32_t read_or_write_dosage_present = read_dosage_present;
      if (mc.plink2_write_flags & kfPlink2WriteLateDosageErase) {
        write_gflags &= ~(kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent);
      } else if (mc.plink2_write_flags & (kfPlink2WriteSetHhMissingKeepDosage | kfPlink2WriteSetMixedMtMissingKeepDosage)) {
        read_or_write_dosage_present = 1;
        write_gflags |= kfPgenGlobalDosagePresent;
      }
      if ((write_gflags & (kfPgenGlobalMultiallelicHardcallFound | kfPgenGlobalDosagePresent)) == (kfPgenGlobalMultiallelicHardcallFound | kfPgenGlobalDosagePresent)) {
        logerrputs("Error: Multiallelic dosages aren't supported yet.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePlink2NoVsort_ret_1;
      }
      write_gflags &= ~kfPgenGlobalMultiallelicHardcallFound;
      uintptr_t alloc_base_cacheline_ct;
      uint64_t mpgw_per_thread_cacheline_ct;
      uint32_t vrec_len_byte_ct;
      uint64_t vblock_cacheline_ct;
      // may want to have a load_sample_ct which is raw_sample_ct when e.g.
      // sample_ct > 0.1 * raw_sample_ct, and sample_ct otherwise.
      MpgwInitPhase1(write_allele_idx_offsets, variant_ct, sample_ct, write_gflags, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);

      // bugfix: each variant currently needs to be vector-aligned
      // bugfix?: need to use raw_sample_ct here, not sample_ct
      const uint32_t raw_sample_ctv2 = NypCtToVecCt(raw_sample_ct);
      const uint32_t max_vblock_size = MINV(kPglVblockSize, variant_ct);
      uint64_t load_vblock_cacheline_ct = VecCtToCachelineCtU64(S_CAST(uint64_t, raw_sample_ctv2) * max_vblock_size);

      if (make_plink2_flags & ((kfMakePlink2MSplitBase * 7) | kfMakePlink2EraseAlt2Plus)) {
        logerrputs("Error: 'multiallelics=' and 'erase-alt2+' modifiers are under development.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePlink2NoVsort_ret_NOMEM;
      }
      if (read_gflags & kfPgenGlobalMultiallelicHardcallFound) {
        // raw multiallelic hardcall track has three parts:
        // 1. two words with rare01_ct and rare10_ct.
        // 2. (vector-aligned) patch_01_set and patch_01_vals.
        // 3. (vector-aligned) patch_10_set and patch_10_vals.
        const uintptr_t mhcraw_word_ct = RoundUpPow2(2, kWordsPerVec) + GetMhcWordCt(raw_sample_ct);
        load_vblock_cacheline_ct += WordCtToCachelineCtU64(S_CAST(uint64_t, mhcraw_word_ct) * max_vblock_size);
      }
      if (read_phase_present) {
        // could make this bound tighter when lots of unphased variants are
        // mixed in among the phased variants, but this isn't nearly as
        // important as the analogous multiallelic optimization

        // phaseraw has three parts:
        // 1. het_ct as uint32_t, and explicit_phasepresent_ct as uint32_t.
        // 2. vec-aligned bitarray of up to (raw_sample_ct + 1) bits.  first
        //    bit is set iff phasepresent is explicitly stored at all (if not,
        //    all hets are assumed to be phased), if yes the remaining bits
        //    store packed phasepresent values for all hets, if no the
        //    remaining bits store packed phaseinfo values for all hets.
        // 3. word-aligned bitarray of up to raw_sample_ct bits, storing
        //    phaseinfo values.  (end of this array is vec-aligned.)
        const uintptr_t phaseraw_word_ct = (8 / kBytesPerWord) + kWordsPerVec + RoundDownPow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec);
        load_vblock_cacheline_ct += WordCtToCachelineCtU64(S_CAST(uint64_t, phaseraw_word_ct) * max_vblock_size);
      }
      if (read_dosage_present) {
        // biallelic dosageraw has two parts:
        // 1. vec-aligned bitarray of up to raw_sample_ct bits, storing which
        //    samples have dosages.
        // 2. word-aligned array of uint16s with 0..32768 fixed-point dosages.
        // dphaseraw has the same structure, with the uint16s replaced with an
        // int16 array of (left - right) values.
        const uintptr_t dosageraw_word_ct = kWordsPerVec * (BitCtToVecCt(raw_sample_ct) + DivUp(raw_sample_ct, kBytesPerVec / sizeof(Dosage)));
        load_vblock_cacheline_ct += WordCtToCachelineCtU64(dosageraw_word_ct * S_CAST(uint64_t, max_vblock_size)) * (1 + read_dphase_present);
      }

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
      mc.refalt1_select = refalt1_select;
      if (refalt1_select) {
        if (write_allele_idx_offsets) {
          // this will require write_mhc and an additional AlleleCode buffer
          logerrputs("Error: Multiallelic allele rotation is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto MakePlink2NoVsort_ret_1;
        }
        if (variant_ct < raw_variant_ct) {
          // might want inner loop to map variant uidx -> idx instead
          STD_ARRAY_PTR_DECL(AlleleCode, 2, tmp_refalt1_select);
          if (BIGSTACK_ALLOC_STD_ARRAY(AlleleCode, 2, variant_ct, &tmp_refalt1_select)) {
            goto MakePlink2NoVsort_fallback;
          }
          uintptr_t variant_uidx_base = 0;
          uintptr_t cur_bits = variant_include[0];
          for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
            const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
            STD_ARRAY_COPY(refalt1_select[variant_uidx], 2, tmp_refalt1_select[variant_idx]);
          }
          mc.refalt1_select = tmp_refalt1_select;
        }
      }
      mpgwp = S_CAST(MTPgenWriter*, bigstack_alloc((calc_thread_ct + DivUp(sizeof(MTPgenWriter), kBytesPerWord)) * sizeof(intptr_t)));
      if (!mpgwp) {
        goto MakePlink2NoVsort_fallback;
      }
      mpgwp->pgen_outfile = nullptr;
      if (bigstack_alloc_wp(calc_thread_ct, &(ctx.loadbuf_thread_starts[0])) ||
          bigstack_alloc_wp(calc_thread_ct, &(ctx.loadbuf_thread_starts[1]))) {
        goto MakePlink2NoVsort_fallback;
      }
      uint32_t nonref_flags_storage = 3;
      uintptr_t* nonref_flags_write = pgfip->nonref_flags;
      if (!nonref_flags_write) {
        nonref_flags_storage = (pgfip->gflags & kfPgenGlobalAllNonref)? 2 : 1;
      } else if (variant_ct < raw_variant_ct) {
        const uint32_t write_variant_ctl = BitCtToWordCt(write_variant_ct);
        uintptr_t* old_nonref_flags = nonref_flags_write;
        if (bigstack_alloc_w(write_variant_ctl, &nonref_flags_write)) {
          goto MakePlink2NoVsort_fallback;
        }
        if (variant_ct == write_variant_ct) {
          CopyBitarrSubset(old_nonref_flags, variant_include, variant_ct, nonref_flags_write);
        } else {
          ZeroWArr(write_variant_ctl, nonref_flags_write);
          if (!write_allele_idx_offsets) {
            SplitNonrefFlags();
          } else {
            JoinNonrefFlags();
          }
        }
        if (nonref_flags_write[0] & 1) {
          if (AllBitsAreOne(nonref_flags_write, write_variant_ct)) {
            BigstackReset(nonref_flags_write);
            nonref_flags_write = nullptr;
            nonref_flags_storage = 2;
          }
        } else if (AllWordsAreZero(nonref_flags_write, write_variant_ctl)) {
          BigstackReset(nonref_flags_write);
          nonref_flags_write = nullptr;
          nonref_flags_storage = 1;
        }
      }
      ctx.pwcs = &(mpgwp->pwcs[0]);
      ctx.new_sample_idx_to_old = new_sample_idx_to_old;
      ctx.thread_write_genovecs = nullptr;
      ctx.thread_write_mhc = nullptr;

      // Each worker thread handles with 64k loaded variants at a time, while
      // the I/O thread loads the next (64k * thread_ct).
      uintptr_t other_per_thread_cacheline_ct = 2 * load_vblock_cacheline_ct;

      uint32_t write_mhc_needed = 0;
      if (new_sample_idx_to_old || subsetting_required) {
        if (bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_genovecs)) {
          goto MakePlink2NoVsort_fallback;
        }
        if (read_phase_present && new_sample_idx_to_old) {
          if (bigstack_alloc_u32(raw_sample_ct, &ctx.old_sample_idx_to_new)) {
            goto MakePlink2NoVsort_fallback;
          }
          for (uint32_t new_sample_idx = 0; new_sample_idx != sample_ct; ++new_sample_idx) {
            ctx.old_sample_idx_to_new[new_sample_idx_to_old[new_sample_idx]] = new_sample_idx;
          }
        }
        // ctx.thread_write_genovecs
        other_per_thread_cacheline_ct += NypCtToCachelineCt(sample_ct);
        write_mhc_needed = 1;
      }
      uintptr_t write_mhcraw_cacheline_ct = 0;
      if (write_mhc_needed) {
        if (bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_mhc)) {
          goto MakePlink2NoVsort_fallback;
        }
        // todo: refalt1_select
        const uintptr_t mhcwrite_word_ct = GetMhcWordCt(sample_ct);
        write_mhcraw_cacheline_ct = DivUp(mhcwrite_word_ct, kWordsPerCacheline);
        other_per_thread_cacheline_ct += write_mhcraw_cacheline_ct;
      }
      if ((write_mhc_needed || read_dosage_present) && new_sample_idx_to_old) {
        // ctx.thread_cumulative_popcount_bufs
        other_per_thread_cacheline_ct += Int32CtToCachelineCt(raw_sample_ctl);
        if (bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_cumulative_popcount_bufs)) {
          goto MakePlink2NoVsort_fallback;
        }
      }
      ctx.thread_write_phasepresents = nullptr;
      ctx.thread_all_hets = nullptr;
      ctx.thread_write_dosagepresents = nullptr;
      ctx.thread_write_dphasepresents = nullptr;
      if (read_or_write_phase_present || read_or_write_dosage_present) {
        if (read_or_write_phase_present) {
          if (bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_phasepresents) ||
              bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_phaseinfos)) {
            goto MakePlink2NoVsort_fallback;
          }
          if (read_phase_present) {
            if (bigstack_alloc_wp(calc_thread_ct, &ctx.thread_all_hets)) {
              goto MakePlink2NoVsort_fallback;
            }
            other_per_thread_cacheline_ct += BitCtToCachelineCt(raw_sample_ct);
          }
          // phasepresent, phaseinfo
          other_per_thread_cacheline_ct += 2 * BitCtToCachelineCt(sample_ct);
        }
        if (read_or_write_dosage_present) {
          if (bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_dosagepresents) ||
              bigstack_alloc_dosagep(calc_thread_ct, &ctx.thread_write_dosagevals)) {
            goto MakePlink2NoVsort_fallback;
          }
          if (read_or_write_dphase_present) {
            if (bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_dphasepresents) ||
                bigstack_alloc_dphasep(calc_thread_ct, &ctx.thread_write_dphasedeltas)) {
              goto MakePlink2NoVsort_fallback;
            }
          }
          // dosage_present, dphase_present
          other_per_thread_cacheline_ct += BitCtToCachelineCt(sample_ct) * (1 + read_or_write_dphase_present);

          // dosage_main, dphase_delta
          other_per_thread_cacheline_ct += DivUp(sample_ct, (kCacheline / sizeof(Dosage))) * (1 + 2 * read_or_write_dphase_present);

          // todo: multiallelic dosage
        }
      }
      if (read_or_write_phase_present || read_dosage_present || (read_gflags & kfPgenGlobalMultiallelicHardcallFound)) {
        // ctx.loaded_vrtypes
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
      ctx.loaded_vrtypes[0] = nullptr;
      ctx.loaded_vrtypes[1] = nullptr;
      if (read_or_write_phase_present || read_dosage_present || (read_gflags & kfPgenGlobalMultiallelicHardcallFound)) {
        ctx.loaded_vrtypes[0] = S_CAST(unsigned char*, bigstack_alloc_raw(kPglVblockSize * calc_thread_ct));
        ctx.loaded_vrtypes[1] = S_CAST(unsigned char*, bigstack_alloc_raw(kPglVblockSize * calc_thread_ct));
      }
      if (read_or_write_phase_present || read_or_write_dosage_present) {
        const uint32_t bitvec_writebuf_byte_ct = BitCtToCachelineCt(sample_ct) * kCacheline;
        const uintptr_t dosagevals_writebuf_byte_ct = DivUp(sample_ct, (kCacheline / 2)) * kCacheline;
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          if (read_or_write_phase_present) {
            ctx.thread_write_phasepresents[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));
            ctx.thread_write_phaseinfos[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));

            if (read_phase_present) {
              ctx.thread_all_hets[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(BitCtToCachelineCt(raw_sample_ct) * kCacheline));
            }
          }
          if (read_or_write_dosage_present) {
            ctx.thread_write_dosagepresents[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));
            ctx.thread_write_dosagevals[tidx] = S_CAST(Dosage*, bigstack_alloc_raw(dosagevals_writebuf_byte_ct));
            if (read_or_write_dphase_present) {
              ctx.thread_write_dphasepresents[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitvec_writebuf_byte_ct));
              ctx.thread_write_dphasedeltas[tidx] = S_CAST(SDosage*, bigstack_alloc_raw(2 * dosagevals_writebuf_byte_ct));
            }
          }
        }
      }
      if (new_sample_idx_to_old || subsetting_required) {
        uintptr_t writebuf_byte_ct = input_biallelic? NypCtToByteCt(sample_ct) : (2 * sample_ct * sizeof(AlleleCode));
        writebuf_byte_ct = RoundUpPow2(writebuf_byte_ct, kCacheline);
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          ctx.thread_write_genovecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(writebuf_byte_ct));
          if (write_mhc_needed) {
            ctx.thread_write_mhc[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(write_mhcraw_cacheline_ct * kCacheline));
          }
          if ((write_mhc_needed || read_dosage_present) && new_sample_idx_to_old) {
            ctx.thread_cumulative_popcount_bufs[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(Int32CtToCachelineCt(raw_sample_ctl) * kCacheline));
          }
        }
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      logprintfww5("Writing %s ... ", outname);
      fputs("0%", stdout);
      fflush(stdout);
      unsigned char* mpgw_alloc = S_CAST(unsigned char*, bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline));
      assert(g_bigstack_base <= g_bigstack_end);
      reterr = MpgwInitPhase2(outname, write_allele_idx_offsets, nonref_flags_write, variant_ct, sample_ct, write_gflags, nonref_flags_storage, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logputs("\n");
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto MakePlink2NoVsort_ret_1;
      }
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto MakePlink2NoVsort_ret_NOMEM;
      }
      mc.sample_include = subsetting_required? sample_include : nullptr;
      ctx.mcp = &mc;
      ctx.spgwp = nullptr;
      ctx.write_reterr = kPglRetSuccess;
      SetThreadFuncAndData(MakePgenThread, &ctx, &tg);

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
      uint32_t cur_batch_size = kPglVblockSize * calc_thread_ct;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uintptr_t read_variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      PgrClearLdCache(simple_pgrp);
      for (uint32_t write_idx_end = 0; ; ++read_batch_idx, write_idx_end += cur_batch_size) {
        if (read_batch_idx) {
          ctx.cur_block_write_ct = cur_batch_size;
          if (write_idx_end == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto MakePlink2NoVsort_ret_THREAD_CREATE_FAIL;
          }
        }
        if (!IsLastBlock(&tg)) {
          if (read_batch_idx == batch_ct_m1) {
            cur_batch_size = variant_ct - (read_batch_idx * kPglVblockSize * calc_thread_ct);
          }
          uintptr_t* cur_loadbuf = main_loadbufs[parity];
          uintptr_t* loadbuf_iter = cur_loadbuf;
          unsigned char* cur_loaded_vrtypes = ctx.loaded_vrtypes[parity];
          for (uint32_t uii = 0; uii != cur_batch_size; ++uii) {
            if (!(uii % kPglVblockSize)) {
              ctx.loadbuf_thread_starts[parity][uii / kPglVblockSize] = loadbuf_iter;
            }
            const uintptr_t read_variant_uidx = BitIter1(variant_include, &read_variant_uidx_base, &cur_bits);
            reterr = PgrGetRaw(read_variant_uidx, read_gflags, simple_pgrp, &loadbuf_iter, cur_loaded_vrtypes? (&(cur_loaded_vrtypes[uii])) : nullptr);
            if (unlikely(reterr)) {
              goto MakePlink2NoVsort_ret_PGR_FAIL;
            }
          }
        }
        if (read_batch_idx) {
          JoinThreads(&tg);
          reterr = ctx.write_reterr;
          if (unlikely(reterr)) {
            // only possible error is kPglRetVarRecordTooLarge?
            goto MakePlink2NoVsort_ret_1;
          }
        }
        parity = 1 - parity;
        if (write_idx_end) {
          reterr = MpgwFlush(mpgwp);
          if (unlikely(reterr)) {
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
      }
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logputs("done.\n");
      // BigstackReset(bigstack_mark);
    } else if (0) {
    MakePlink2NoVsort_fallback:
      g_failed_alloc_attempt_size = 0;
      mpgwp = nullptr;
      BigstackReset(bigstack_mark2);
      reterr = MakePgenRobust(sample_include, new_sample_idx_to_old, variant_include, allele_idx_offsets, allele_presents, refalt1_select, write_allele_idx_offsets, nullptr, ctx.sex_male_collapsed, ctx.sex_female_collapsed, raw_variant_ct, variant_ct, write_variant_ct, max_allele_ct, hard_call_thresh, dosage_erase_thresh, make_plink2_flags, &mc, simple_pgrp, outname, outname_end);
      if (unlikely(reterr)) {
        goto MakePlink2NoVsort_ret_1;
      }
      if (variant_ct != write_variant_ct) {
        logprintfww("Multiallelic %s: %u variant%s written.\n", (variant_ct < write_variant_ct)? "split" : "join", write_variant_ct, (write_variant_ct == 1)? "" : "s");
      }
    }
  }
  while (0) {
  MakePlink2NoVsort_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakePlink2NoVsort_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  MakePlink2NoVsort_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  MakePlink2NoVsort_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  MakePlink2NoVsort_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 MakePlink2NoVsort_ret_1:
  CleanupMpgw(mpgwp, &reterr);
  CleanupThreads(&tg);
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
  for (uint32_t chr_idx = 0; chr_idx != chr_code_end; ++chr_idx) {
    const uint32_t cur_chr_size = chr_idx_to_size[chr_idx];
    if (cur_chr_size) {
      ++new_chr_ct;
    }
  }
  // bugfix (25 Nov 2019): must add 1 for chr_fo_vidx_start
  if (bigstack_alloc_u32(new_chr_ct, &(write_cip->chr_file_order)) ||
      bigstack_alloc_u32(new_chr_ct + 1, &(write_cip->chr_fo_vidx_start))) {
    return 1;
  }
  write_cip->chr_ct = new_chr_ct;
  // now for the actual sorting.
  // autosomes and PAR1/X/PAR2/Y/XY/MT come first, then contig names.
  const uint32_t autosome_ct = cip->autosome_ct;
  const uint32_t xymt_ct = max_code - autosome_ct;
  const uint32_t autosome_ct_p1 = autosome_ct + 1;

  STD_ARRAY_KREF(uint32_t, kChrOffsetCt) xymt_codes = cip->xymt_codes;
  const uintptr_t xymt_idx_to_chr_sort_offset[kChrOffsetCt] = {1, 3, 4, 5, 0, 2};

  // chr_sort_idx in high bits, original chr_idx in low
  uint64_t* std_sortbuf;
  uint64_t* std_sortbuf_iter;
  if (bigstack_alloc_u64(max_code + 1, &std_sortbuf)) {
    return 1;
  }
  std_sortbuf_iter = std_sortbuf;
  for (uintptr_t chr_idx = 0; chr_idx <= autosome_ct; ++chr_idx) {
    if (chr_idx_to_size[chr_idx]) {
      *std_sortbuf_iter++ = chr_idx * 0x100000001LLU;
    }
  }
  for (uint32_t xymt_idx = 0; xymt_idx != xymt_ct; ++xymt_idx) {
    const uint32_t xymt_code = xymt_codes[xymt_idx];
    if (!IsI32Neg(xymt_code)) {
      if (chr_idx_to_size[xymt_idx + autosome_ct_p1]) {
        *std_sortbuf_iter++ = (S_CAST(uint64_t, xymt_idx_to_chr_sort_offset[xymt_idx] + autosome_ct_p1) << 32) | (xymt_idx + autosome_ct_p1);
      }
    }
  }
  const uint32_t std_sortbuf_len = std_sortbuf_iter - std_sortbuf;
  STD_SORT(std_sortbuf_len, u64cmp, std_sortbuf);
  uint32_t write_vidx = 0;
  write_cip->chr_fo_vidx_start[0] = 0;
  for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx != std_sortbuf_len; ++new_chr_fo_idx) {
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
    for (uint32_t chr_idx = max_code + 1; chr_idx != chr_code_end; ++chr_idx) {
      if (chr_idx_to_size[chr_idx]) {
        nonstd_sort_buf[str_idx].strptr = nonstd_names[chr_idx];
        nonstd_sort_buf[str_idx].orig_idx = chr_idx;
        ++str_idx;
      }
    }
    assert(str_idx == new_nonstd_ct);
    // nonstd_names are not allocated in main workspace, so can't overread.
    StrptrArrSortMain(new_nonstd_ct, 0, use_nsort, nonstd_sort_buf);
    uint32_t new_chr_fo_idx = std_sortbuf_len;
    for (str_idx = 0; str_idx != new_nonstd_ct; ++str_idx, ++new_chr_fo_idx) {
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
PglErr WriteBimResorted(const char* outname, const ChrInfo* write_cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t output_zst, uint32_t thread_ct) {
  // allele_presents must be nullptr unless we're trimming alt alleles
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WriteBimResorted_ret_NOMEM;
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WriteBimResorted_ret_1;
    }

    const char output_missing_geno_char = *g_output_missing_geno_ptr;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
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
      cswritep = u32toa(variant_bps[variant_uidx], cswritep);
      *cswritep++ = '\t';
      const uintptr_t allele_idx_offset_base = allele_idx_offsets? allele_idx_offsets[variant_uidx] : (variant_uidx * 2);
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      // note that VCF ref allele corresponds to A2, not A1
      if (!refalt1_select) {
        if ((!allele_presents) || IsSet(allele_presents, 1 + allele_idx_offset_base)) {
          cswritep = strcpya(cswritep, cur_alleles[1]);
        } else {
          *cswritep++ = output_missing_geno_char;
        }
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[0]);
      } else {
        STD_ARRAY_KREF(AlleleCode, 2) cur_refalt1_select = refalt1_select[variant_uidx];
        if ((!allele_presents) || IsSet(allele_presents, cur_refalt1_select[1] + allele_idx_offset_base)) {
          cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[1]]);
        } else {
          *cswritep++ = output_missing_geno_char;
        }
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[cur_refalt1_select[0]]);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto WriteBimResorted_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
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

PglErr PvarInfoReloadInterval(const uint32_t* old_variant_uidx_to_new, uint32_t variant_idx_start, uint32_t variant_idx_end, TextStream* pvar_reload_txsp, char** pvar_info_strs) {
  // We assume the batch size was chosen such that there's no risk of
  // scribbling past g_bigstack_end (barring pathological cases like another
  // process modifying the .pvar file after initial load).
  // We also assume no more dynamic allocations are needed after this;
  // otherwise, str_store_iter should be returned.
  char* line_iter;
  // probable todo: avoid rewind when one batch is entirely after the previous
  // batch (this is likely when input was already almost-sorted, and just a few
  // coordinates changed due to e.g. --normalize)
  PglErr reterr = TextRewind(pvar_reload_txsp);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uint32_t cur_batch_size = variant_idx_end - variant_idx_start;
  char* str_store_iter = R_CAST(char*, g_bigstack_base);
  uint32_t info_col_idx;
  reterr = PvarInfoReloadHeader(pvar_reload_txsp, &line_iter, &info_col_idx);
  if (unlikely(reterr)) {
    return reterr;
  }
  uint32_t variant_idx = 0;
  for (uint32_t variant_uidx = 0; ; ++variant_uidx) {
    reterr = TextNextLineLstrip(pvar_reload_txsp, &line_iter);
    if (unlikely(reterr)) {
      return reterr;
    }
    const uint32_t new_variant_idx_offset = old_variant_uidx_to_new[variant_uidx] - variant_idx_start;
    // exploit wraparound, UINT32_MAX null value
    if (new_variant_idx_offset >= cur_batch_size) {
      continue;
    }
    line_iter = NextTokenMultFar(line_iter, info_col_idx);
    if (!line_iter) {
      return kPglRetRewindFail;
    }
    char* info_end = CurTokenEnd(line_iter);
    const uint32_t info_slen = info_end - line_iter;
    pvar_info_strs[new_variant_idx_offset] = str_store_iter;
    str_store_iter = memcpyax(str_store_iter, line_iter, info_slen, '\0');
    line_iter = info_end;
    if (++variant_idx == cur_batch_size) {
      break;
    }
  }
  assert(str_store_iter <= R_CAST(char*, g_bigstack_end));
  return kPglRetSuccess;
}

// could be BoolErr
PglErr WritePvarResortedInterval(const ChrInfo* write_cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t variant_idx_start, uint32_t variant_idx_end, uint32_t info_pr_flag_present, uint32_t write_qual, uint32_t write_filter, uint32_t write_info, uint32_t all_nonref, uint32_t write_cm, char** pvar_info_strs, CompressStreamState* cssp, char** cswritepp, uint32_t* chr_fo_idxp, uint32_t* chr_endp, uint32_t* chr_buf_blenp, char* chr_buf) {
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
    for (uint32_t variant_idx = variant_idx_start; variant_idx != variant_idx_end; ++variant_idx) {
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
      cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx][0];
        alt1_allele_idx = refalt1_select[variant_uidx][1];
      }
      cswritep = strcpyax(cswritep, cur_alleles[ref_allele_idx], '\t');
      uint32_t alt_allele_written = 0;
      if ((!allele_presents) || IsSet(allele_presents, allele_idx_offset_base + alt1_allele_idx)) {
        cswritep = strcpya(cswritep, cur_alleles[alt1_allele_idx]);
        alt_allele_written = 1;
      }
      if (unlikely(Cswrite(cssp, &cswritep))) {
        goto WritePvarResortedInterval_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if ((allele_idx == ref_allele_idx) || (allele_idx == alt1_allele_idx) || (allele_presents && (!IsSet(allele_presents, allele_idx_offset_base + allele_idx)))) {
            continue;
          }
          if (alt_allele_written) {
            *cswritep++ = ',';
          }
          alt_allele_written = 1;
          cswritep = strcpya(cswritep, cur_alleles[allele_idx]);
          if (unlikely(Cswrite(cssp, &cswritep))) {
            goto WritePvarResortedInterval_ret_WRITE_FAIL;
          }
        }
      }
      if (!alt_allele_written) {
        *cswritep++ = output_missing_geno_char;
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
          cswritep = strcpya_k(cswritep, "PASS");
        } else {
          cswritep = strcpya(cswritep, filter_storage[variant_uidx]);
        }
      }

      if (write_info) {
        *cswritep++ = '\t';
        const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
        if (pvar_info_strs) {
          PvarInfoWrite(info_pr_flag_present, is_pr, pvar_info_strs[variant_idx - variant_idx_start], &cswritep);
        } else {
          if (is_pr) {
            cswritep = strcpya_k(cswritep, "PR");
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

// allele_presents must be nullptr unless we're trimming alt alleles.
//
// The annoying part of this is handling a sequence of INFO strings that don't
// fit in memory; we use a multipass approach for that.  File creation,
// allocation of buffers, and generating the header line occurs directly in
// this function, while loading the next pvar_info_strs batch and writing the
// next .pvar line batch are one level down.
PglErr WritePvarResorted(const char* outname, const uintptr_t* variant_include, const ChrInfo* write_cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* qual_present, const float* quals, const uintptr_t* filter_present, const uintptr_t* filter_npass, const char* const* filter_storage, const uintptr_t* nonref_flags, const char* pvar_info_reload, const double* variant_cms, const uint32_t* new_variant_idx_to_old, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t nonref_flags_storage, uint32_t max_filter_slen, uint32_t info_reload_slen, PvarPsamFlags pvar_psam_flags, uint32_t thread_ct, char* xheader) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  TextStream pvar_reload_txs;
  PreinitCstream(&css);
  PreinitTextStream(&pvar_reload_txs);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(write_cip) + 1;
    // includes trailing tab
    char* chr_buf;

    if (unlikely(
                 bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WritePvarResorted_ret_NOMEM;
    }
    uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_filter_slen + info_reload_slen;
    if (overflow_buf_size < 2 * kCompressStreamBlock) {
      overflow_buf_size = 2 * kCompressStreamBlock;
    }
    const uint32_t output_zst = (pvar_psam_flags / kfPvarZs) & 1;
    reterr = InitCstreamAlloc(outname, 0, output_zst, thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WritePvarResorted_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t all_nonref = (nonref_flags_storage == 2);
    uint32_t write_info_pr = all_nonref;
    uint32_t write_info = (pvar_psam_flags & kfPvarColInfo) || pvar_info_reload;
    if (write_info && nonref_flags) {
      write_info_pr = !IntersectionIsEmpty(variant_include, nonref_flags, raw_variant_ctl);
    }
    write_info_pr = write_info_pr && write_info;
    if (unlikely(write_info_pr && (info_flags & kfInfoPrNonflagPresent))) {
      logputs("\n");
      logerrputs("Error: Conflicting INFO:PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO:PR key.\n");
      goto WritePvarResorted_ret_INCONSISTENT_INPUT;
    }

    uint32_t write_filter = 0;
    if (pvar_psam_flags & kfPvarColFilter) {
      write_filter = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybefilter) && filter_present) {
      write_filter = !IntersectionIsEmpty(variant_include, filter_present, raw_variant_ctl);
    }
    const uint32_t info_pr_flag_present = (info_flags / kfInfoPrFlagPresent) & 1;
    if (pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader)) {
      reterr = PvarXheaderWrite(nullptr, write_cip, variant_bps, allele_idx_offsets, allele_storage, new_variant_idx_to_old, xheader_blen, (pvar_psam_flags / kfPvarColVcfheader) & 1, write_filter, write_info, write_info_pr && (!info_pr_flag_present), max_allele_slen, kfUnsortedVar0, xheader, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WritePvarResorted_ret_1;
      }
    }
    if (write_cip->chrset_source) {
      AppendChrsetLine(write_cip, &cswritep);
    }
    cswritep = strcpya_k(cswritep, "#CHROM\tPOS\tID\tREF\tALT");

    uint32_t write_qual = 0;
    if (pvar_psam_flags & kfPvarColQual) {
      write_qual = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybequal) && qual_present) {
      write_qual = !IntersectionIsEmpty(variant_include, qual_present, raw_variant_ctl);
    }
    if (write_qual) {
      cswritep = strcpya_k(cswritep, "\tQUAL");
    }

    if (write_filter) {
      cswritep = strcpya_k(cswritep, "\tFILTER");
    }

    if (write_info) {
      cswritep = strcpya_k(cswritep, "\tINFO");
    }

    uint32_t write_cm = 0;
    if (pvar_psam_flags & kfPvarColCm) {
      write_cm = 1;
    } else if ((pvar_psam_flags & kfPvarColMaybecm) && variant_cms) {
      if (raw_variant_ct == variant_ct) {
        // nonzero_cm_present check was performed
        write_cm = 1;
      } else {
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          if (variant_cms[variant_uidx] != 0.0) {
            write_cm = 1;
            break;
          }
        }
      }
    }
    if (write_cm) {
      cswritep = strcpya_k(cswritep, "\tCM");
    }
    AppendBinaryEoln(&cswritep);

    uint32_t* old_variant_uidx_to_new = nullptr;
    char** pvar_info_strs = nullptr;
    uint32_t batch_size = variant_ct;
    uint32_t batch_ct = 1;
    if (pvar_info_reload) {
      if (unlikely(bigstack_alloc_u32(raw_variant_ct, &old_variant_uidx_to_new))) {
        goto WritePvarResorted_ret_NOMEM;
      }
      SetAllU32Arr(raw_variant_ct, old_variant_uidx_to_new);
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t old_variant_uidx = new_variant_idx_to_old[variant_idx];
        old_variant_uidx_to_new[old_variant_uidx] = variant_idx;
      }

      uint32_t decompress_thread_ct = 1;
      if (!output_zst) {
        decompress_thread_ct = thread_ct - 1;
        if (!decompress_thread_ct) {
          decompress_thread_ct = 1;
        }
      }
      reterr = SizeAndInitTextStream(pvar_info_reload, bigstack_left() / 4, decompress_thread_ct, &pvar_reload_txs);
      if (unlikely(reterr)) {
        goto WritePvarResorted_ret_TSTREAM_FAIL;
      }

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
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t batch_idx = 0; batch_idx != batch_ct; ++batch_idx) {
      if (variant_idx_start >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx_start * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
      uint32_t variant_idx_end = MINV(variant_idx_start + batch_size, variant_ct);
      if (pvar_info_reload) {
        reterr = PvarInfoReloadInterval(old_variant_uidx_to_new, variant_idx_start, variant_idx_end, &pvar_reload_txs, pvar_info_strs);
        if (unlikely(reterr)) {
          goto WritePvarResorted_ret_TSTREAM_FAIL;
        }
      }
      reterr = WritePvarResortedInterval(write_cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_presents, refalt1_select, qual_present, quals, filter_present, filter_npass, filter_storage, nonref_flags, variant_cms, new_variant_idx_to_old, variant_idx_start, variant_idx_end, info_pr_flag_present, write_qual, write_filter, write_info, all_nonref, write_cm, pvar_info_strs, &css, &cswritep, &chr_fo_idx, &chr_end, &chr_buf_blen, chr_buf);
      if (unlikely(reterr)) {
        goto WritePvarResorted_ret_1;
      }
      variant_idx_start = variant_idx_end;
    }

    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WritePvarResorted_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  WritePvarResorted_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WritePvarResorted_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_reload_txs);
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
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr MakePlink2Vsort(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* new_sample_idx_to_old, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* allele_presents, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ChrIdx* chr_idxs, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, MakePlink2Flags make_plink2_flags, uint32_t use_nsort, PvarPsamFlags pvar_psam_flags, char* xheader, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    // Resort the variants.
    // 1. (todo) Apply --update-chr if necessary.
    // 2. Count number of remaining variants in each chromosome, then sort the
    //    chromosomes.
    // 3. Within each chromosome, sort by position.  Could add 0.5 for
    //    non-SNPs (not currently implemented)?  Could multithread this by
    //    chromosome, and/or use C++17 multithreaded sort, but INFO-reload is a
    //    much bigger bottleneck in practice.
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
    STD_ARRAY_COPY(cip->xymt_codes, kChrOffsetCt, write_chr_info.xymt_codes);
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
    if (unlikely(
            bigstack_calloc_w(kChrMaskWords, &write_chr_info.chr_mask) ||
            bigstack_alloc_u32(chr_code_end, &write_chr_info.chr_idx_to_foidx) ||
            bigstack_end_calloc_u32(chr_code_end, &chr_idx_to_size))) {
      goto MakePlink2Vsort_ret_NOMEM;
    }
    SetAllU32Arr(chr_code_end, write_chr_info.chr_idx_to_foidx);
    if (chr_idxs) {
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_base = variant_include[0];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_base);
        chr_idx_to_size[chr_idxs[variant_uidx]] += 1;
      }
      for (uint32_t chr_idx = 0; chr_idx != chr_code_end; ++chr_idx) {
        if (chr_idx_to_size[chr_idx]) {
          SetBit(chr_idx, write_chr_info.chr_mask);
        }
      }
      // bugfix: chr_file_order is invalid
    } else {
      const uint32_t* chr_fo_vidx_start = cip->chr_fo_vidx_start;
      const uint32_t orig_chr_ct = cip->chr_ct;
      uint32_t vidx_start = 0;
      for (uint32_t chr_fo_idx = 0; chr_fo_idx != orig_chr_ct; ++chr_fo_idx) {
        const uint32_t vidx_end = chr_fo_vidx_start[chr_fo_idx + 1];
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_idx_to_size[chr_idx] = PopcountBitRange(variant_include, vidx_start, vidx_end);
        if (chr_idx_to_size[chr_idx]) {
          SetBit(chr_idx, write_chr_info.chr_mask);
        }
        vidx_start = vidx_end;
      }
    }
    if (unlikely(SortChr(cip, chr_idx_to_size, use_nsort, &write_chr_info))) {
      goto MakePlink2Vsort_ret_NOMEM;
    }

    uint32_t* new_variant_idx_to_old;

    // pos_vidx_sort_buf has variant_bp in high bits, variant_uidx in low
    uint64_t* pos_vidx_sort_buf;
    if (unlikely(
            bigstack_alloc_u32(variant_ct, &new_variant_idx_to_old) ||
            bigstack_alloc_u64(variant_ct + 1, &pos_vidx_sort_buf))) {
      goto MakePlink2Vsort_ret_NOMEM;
    }
    pos_vidx_sort_buf[variant_ct] = ~0LLU;
    const uint32_t new_chr_ct = write_chr_info.chr_ct;
    if (chr_idxs) {
      uint32_t* next_write_vidxs;
      if (unlikely(bigstack_alloc_u32(chr_code_end, &next_write_vidxs))) {
        goto MakePlink2Vsort_ret_NOMEM;
      }
      for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx != new_chr_ct; ++new_chr_fo_idx) {
        const uint32_t chr_idx = write_chr_info.chr_file_order[new_chr_fo_idx];
        next_write_vidxs[chr_idx] = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx];
      }
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        const uint32_t chr_idx = chr_idxs[variant_uidx];
        const uint32_t write_vidx = next_write_vidxs[chr_idx];
        pos_vidx_sort_buf[write_vidx] = (S_CAST(uint64_t, variant_bps[variant_uidx]) << 32) | variant_uidx;
        next_write_vidxs[chr_idx] += 1;
      }
      BigstackReset(next_write_vidxs);
    } else {
      uint32_t old_chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t chr_idx = 0;
      uint32_t write_vidx = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx, ++write_vidx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        if (variant_uidx >= chr_end) {
          do {
            ++old_chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[old_chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          chr_idx = cip->chr_file_order[old_chr_fo_idx];
          // bugfix (8 Sep 2018): write_vidx was set to the wrong value here
          const uint32_t new_chr_fo_idx = write_chr_info.chr_idx_to_foidx[chr_idx];
          write_vidx = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx];
        }
        pos_vidx_sort_buf[write_vidx] = (S_CAST(uint64_t, variant_bps[variant_uidx]) << 32) | variant_uidx;
      }
    }

    StrSortIndexedDeref* same_pos_sort_buf = R_CAST(StrSortIndexedDeref*, g_bigstack_base);
    const uintptr_t same_pos_sort_buf_size = bigstack_left() / sizeof(StrSortIndexedDeref);

    uint32_t vidx_start = 0;
    uint32_t* new_variant_idx_to_old_iter = new_variant_idx_to_old;
    for (uint32_t new_chr_fo_idx = 0; new_chr_fo_idx != new_chr_ct; ++new_chr_fo_idx) {
      const uint32_t vidx_end = write_chr_info.chr_fo_vidx_start[new_chr_fo_idx + 1];
      const uint32_t chr_size = vidx_end - vidx_start;
      const uint64_t post_entry = pos_vidx_sort_buf[vidx_end];
      pos_vidx_sort_buf[vidx_end] = ~0LLU;  // simplify end-of-chromosome logic
      uint64_t* pos_vidx_sort_chr = &(pos_vidx_sort_buf[vidx_start]);
      STD_SORT_PAR_UNSEQ(chr_size, u64cmp, pos_vidx_sort_chr);
      uint32_t prev_pos = pos_vidx_sort_chr[0] >> 32;
      uint32_t prev_variant_uidx = S_CAST(uint32_t, pos_vidx_sort_chr[0]);
      uint32_t prev_cidx = 0;
      uint32_t cidx = 1;
      // is chr_size == 0 possible here?  document if this code is revisited.
      for (; cidx < chr_size; ++cidx) {
        uint64_t cur_entry = pos_vidx_sort_chr[cidx];
        uint32_t cur_pos = cur_entry >> 32;
        if (cur_pos == prev_pos) {
          same_pos_sort_buf[0].strptr = variant_ids[prev_variant_uidx];
          same_pos_sort_buf[0].orig_idx = prev_variant_uidx;
          uint32_t equal_pos_ct = 1;
          const uint64_t* pos_vidx_sort_chr2 = &(pos_vidx_sort_chr[prev_cidx]);
          do {
            if (unlikely(equal_pos_ct >= same_pos_sort_buf_size)) {
              goto MakePlink2Vsort_ret_NOMEM;
            }
            const uint32_t variant_uidx = S_CAST(uint32_t, cur_entry);
            same_pos_sort_buf[equal_pos_ct].strptr = variant_ids[variant_uidx];
            same_pos_sort_buf[equal_pos_ct].orig_idx = variant_uidx;
            cur_entry = pos_vidx_sort_chr2[++equal_pos_ct];
            cur_pos = cur_entry >> 32;
          } while (cur_pos == prev_pos);
          StrptrArrSortMain(equal_pos_ct, 1, use_nsort, same_pos_sort_buf);
          for (uint32_t equal_pos_idx = 0; equal_pos_idx != equal_pos_ct; ++equal_pos_idx) {
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

    if (make_plink2_flags & kfMakeBim) {
      const uint32_t bim_zst = (make_plink2_flags / kfMakeBimZs) & 1;
      OutnameZstSet(".bim", bim_zst, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);

      reterr = WriteBimResorted(outname, &write_chr_info, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_presents, refalt1_select, variant_cms, new_variant_idx_to_old, variant_ct, max_allele_slen, bim_zst, max_thread_ct);
      if (unlikely(reterr)) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePvar) {
      OutnameZstSet(".pvar", pvar_psam_flags & kfPvarZs, outname_end);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t nonref_flags_storage = 3;
      if (!PgrGetNonrefFlags(simple_pgrp)) {
        nonref_flags_storage = (PgrGetGflags(simple_pgrp) & kfPgenGlobalAllNonref)? 2 : 1;
      }
      reterr = WritePvarResorted(outname, variant_include, &write_chr_info, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_presents, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, PgrGetNonrefFlags(simple_pgrp), pvar_info_reload, variant_cms, new_variant_idx_to_old, raw_variant_ct, variant_ct, max_allele_slen, xheader_blen, info_flags, nonref_flags_storage, max_filter_slen, info_reload_slen, pvar_psam_flags, max_thread_ct, xheader);
      if (unlikely(reterr)) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakeFam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, new_sample_idx_to_old, sample_ct, pheno_ct, '\t');
      if (unlikely(reterr)) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & kfMakePsam) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WritePsam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, pheno_names, new_sample_idx_to_old, sample_ct, pheno_ct, max_pheno_name_blen, pvar_psam_flags);
      if (unlikely(reterr)) {
        goto MakePlink2Vsort_ret_1;
      }
      logputs("done.\n");
    }
    if (make_plink2_flags & (kfMakeBed | kfMakePgen)) {
      // boilerplate from start of MakePlink2NoVsort()
      if (make_plink2_flags & kfMakePlink2MMask) {
        logerrputs("Error: --make-bed/--make-[b]pgen multiallelics= is currently under development.\n");
        reterr = kPglRetNotYetSupported;
        goto MakePlink2Vsort_ret_1;
      }
      MakeCommon mc;
      mc.plink2_write_flags = kfPlink2Write0;
      mc.raw_sample_ct = raw_sample_ct;
      mc.sample_ct = sample_ct;
      uintptr_t* sex_male_collapsed = nullptr;
      uintptr_t* sex_female_collapsed = nullptr;
      if (make_plink2_flags & kfMakePlink2SetHhMissing) {
        const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
        const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
        uintptr_t* new_sex_male;
        if (unlikely(
                bigstack_alloc_w(sample_ctv * kWordsPerVec, &new_sex_male) ||
                bigstack_alloc_w(sample_ctv * kWordsPerVec, &mc.sex_male_collapsed_interleaved) ||
                bigstack_alloc_w(sample_ctv * kWordsPerVec, &sex_female_collapsed) ||
                bigstack_alloc_w(sample_ctv * kWordsPerVec, &mc.sex_female_collapsed_interleaved))) {
          goto MakePlink2Vsort_ret_NOMEM;
        }
        CopyBitarrSubset(sex_male, sample_include, sample_ct, new_sex_male);
        ZeroTrailingWords(sample_ctl, new_sex_male);
        sex_male_collapsed = new_sex_male;
        FillInterleavedMaskVec(sex_male_collapsed, sample_ctv, mc.sex_male_collapsed_interleaved);

        CopyBitarrSubset(sex_nm, sample_include, sample_ct, sex_female_collapsed);
        BitvecInvmask(new_sex_male, sample_ctl, sex_female_collapsed);
        ZeroTrailingWords(sample_ctl, sex_female_collapsed);
        FillInterleavedMaskVec(sex_female_collapsed, sample_ctv, mc.sex_female_collapsed_interleaved);

        mc.plink2_write_flags |= kfPlink2WriteSetHhMissing;
      } else {
        // defensive
        mc.sex_male_collapsed_interleaved = nullptr;
        mc.sex_female_collapsed_interleaved = nullptr;
      }
      if (make_plink2_flags & kfMakePlink2SetMixedMtMissing) {
        mc.plink2_write_flags |= kfPlink2WriteSetMixedMtMissing;
      }
      mc.cip = &write_chr_info;
      const uintptr_t* write_allele_idx_offsets = nullptr;
      if (allele_idx_offsets && (!(make_plink2_flags & kfMakePlink2EraseAlt2Plus))) {
        if ((variant_ct < raw_variant_ct) || new_variant_idx_to_old) {
          uintptr_t* new_allele_idx_offsets;
          if (unlikely(bigstack_alloc_w(variant_ct + 1, &new_allele_idx_offsets))) {
            goto MakePlink2Vsort_ret_NOMEM;
          }
          const uintptr_t final_offset = InitWriteAlleleIdxOffsets(variant_include, allele_idx_offsets, allele_presents, refalt1_select, new_variant_idx_to_old, variant_ct, new_allele_idx_offsets);
          if (final_offset != 2 * variant_ct) {
            new_allele_idx_offsets[variant_ct] = final_offset;
            write_allele_idx_offsets = new_allele_idx_offsets;
          } else {
            BigstackReset(new_allele_idx_offsets);
          }
        } else {
          write_allele_idx_offsets = allele_idx_offsets;
        }
      }
      reterr = MakePgenRobust(sample_include, new_sample_idx_to_old, variant_include, allele_idx_offsets, allele_presents, refalt1_select, write_allele_idx_offsets, new_variant_idx_to_old, sex_male_collapsed, sex_female_collapsed, raw_variant_ct, variant_ct, variant_ct, max_allele_ct, hard_call_thresh, dosage_erase_thresh, make_plink2_flags, &mc, simple_pgrp, outname, outname_end);
      if (unlikely(reterr)) {
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
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    char* idbuf;
    uintptr_t* already_seen;
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ct, new_sample_idx_to_old_ptr) ||
            bigstack_alloc_c(siip->max_sample_id_blen, &idbuf) ||
            bigstack_calloc_w(BitCtToWordCt(raw_sample_ct), &already_seen))) {
      goto SampleSortFileMap_ret_NOMEM;
    }

    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() - (bigstack_left() / 4), &max_line_blen))) {
      goto SampleSortFileMap_ret_NOMEM;
    }
    char* line_start;
    XidMode xid_mode;
    reterr = OpenAndLoadXidHeader(sample_sort_fname, "indiv-sort", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, max_line_blen, &txs, &xid_mode, &line_idx, &line_start, nullptr);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: --indiv-sort file is empty.\n");
        goto SampleSortFileMap_ret_MALFORMED_INPUT;
      }
      goto SampleSortFileMap_ret_TSTREAM_XID_FAIL;
    }
    uint32_t* xid_map;
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto SampleSortFileMap_ret_1;
    }
    uint32_t* new_sample_idx_to_old_iter = *new_sample_idx_to_old_ptr;
    if (*line_start == '#') {
      ++line_idx;
      line_start = TextGet(&txs);
    }
    for (; line_start; ++line_idx, line_start = TextGet(&txs)) {
      if (unlikely(line_start[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --indiv-sort file starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx);
        goto SampleSortFileMap_ret_MALFORMED_INPUT_WW;
      }
      const char* linebuf_iter = line_start;
      uint32_t sample_uidx;
      if (!SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx, idbuf)) {
        if (unlikely(IsSet(already_seen, sample_uidx))) {
          char* tab_iter = AdvToDelim(idbuf, '\t');
          *tab_iter = ' ';
          if (xid_mode & kfXidModeFlagSid) {
            *AdvToDelim(&(tab_iter[1]), '\t') = ' ';
          }
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID '%s' in --indiv-sort file.\n", idbuf);
          goto SampleSortFileMap_ret_MALFORMED_INPUT_WW;
        }
        SetBit(sample_uidx, already_seen);
        *new_sample_idx_to_old_iter++ = sample_uidx;
      } else if (unlikely(!linebuf_iter)) {
        goto SampleSortFileMap_ret_MISSING_TOKENS;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto SampleSortFileMap_ret_TSTREAM_FAIL;
    }
    if (unlikely(S_CAST(uintptr_t, new_sample_idx_to_old_iter - (*new_sample_idx_to_old_ptr)) != sample_ct)) {
      logerrputs("Error: --indiv-sort file does not contain all loaded sample IDs.\n");
      goto SampleSortFileMap_ret_INCONSISTENT_INPUT;
    }
    bigstack_mark = R_CAST(unsigned char*, idbuf);
  }
  while (0) {
  SampleSortFileMap_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SampleSortFileMap_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  SampleSortFileMap_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  SampleSortFileMap_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  SampleSortFileMap_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--indiv-sort file", &txs);
    break;
  SampleSortFileMap_ret_MISSING_TOKENS:
    logerrprintf("Error: Line %" PRIuPTR " of --indiv-sort file has fewer tokens than expected.\n", line_idx);
  SampleSortFileMap_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 SampleSortFileMap_ret_1:
  CleanupTextStream2("--indiv-sort file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
