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

#include "plink2_adjust.h"
#include "plink2_compress_stream.h"
#include "plink2_stats.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitAdjust(AdjustInfo* adjust_info_ptr, AdjustFileInfo* adjust_file_info_ptr) {
  adjust_info_ptr->flags = kfAdjust0;
  adjust_info_ptr->lambda = 0.0;
  adjust_file_info_ptr->base.flags = kfAdjust0;
  adjust_file_info_ptr->base.lambda = 0.0;
  adjust_file_info_ptr->fname = nullptr;
  adjust_file_info_ptr->test_name = nullptr;
  adjust_file_info_ptr->chr_field = nullptr;
  adjust_file_info_ptr->pos_field = nullptr;
  adjust_file_info_ptr->id_field = nullptr;
  adjust_file_info_ptr->ref_field = nullptr;
  adjust_file_info_ptr->alt_field = nullptr;
  adjust_file_info_ptr->test_field = nullptr;
  adjust_file_info_ptr->p_field = nullptr;
}

void CleanupAdjust(AdjustFileInfo* adjust_file_info_ptr) {
  free_cond(adjust_file_info_ptr->alt_field);
  free_cond(adjust_file_info_ptr->chr_field);
  if (adjust_file_info_ptr->fname) {
    free(adjust_file_info_ptr->fname);
    free_cond(adjust_file_info_ptr->pos_field);
    free_cond(adjust_file_info_ptr->id_field);
    free_cond(adjust_file_info_ptr->ref_field);
    free_cond(adjust_file_info_ptr->test_field);
    free_cond(adjust_file_info_ptr->p_field);
  }
}

typedef struct AdjAssocResultStruct {
  double ln_pval;
  double chisq;  // do we really need this?...
  uint32_t variant_uidx;
#ifdef __cplusplus
  bool operator<(const struct AdjAssocResultStruct& rhs) const {
    return ln_pval < rhs.ln_pval;
  }
#endif
} AdjAssocResult;

static inline void adjust_print_ln(const char* output_min_p_str, double ln_pval, double output_min_ln, uint32_t output_min_p_slen, uint32_t is_neglog10, char** bufpp) {
  **bufpp = '\t';
  *bufpp += 1;
  if (ln_pval <= output_min_ln) {
    *bufpp = memcpya(*bufpp, output_min_p_str, output_min_p_slen);
  } else {
    if (!is_neglog10) {
      *bufpp = lntoa_g(ln_pval, *bufpp);
    } else {
      *bufpp = dtoa_g(ln_pval * (-1.0 / kLn10), *bufpp);
    }
  }
}

// Now based around ln_pvals, to allow useful comparisons < 5e-324.
PglErr Multcomp(const uintptr_t* variant_include, const ChrInfo* cip, const char* const* chr_ids, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const AdjustInfo* adjust_info_ptr, const double* ln_pvals, const double* chisqs, uint32_t orig_variant_ct, uint32_t max_allele_slen, double ln_pfilter, double output_min_ln, uint32_t skip_gc, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    AdjAssocResult* sortbuf = S_CAST(AdjAssocResult*, bigstack_alloc(orig_variant_ct * sizeof(AdjAssocResult)));
    if (!sortbuf) {
      goto Multcomp_ret_NOMEM;
    }
    uint32_t valid_variant_ct = 0;
    if (chisqs) {
      uint32_t variant_uidx = 0;
      if (ln_pvals) {
        for (uint32_t vidx = 0; vidx < orig_variant_ct; ++vidx, ++variant_uidx) {
          MovU32To1Bit(variant_include, &variant_uidx);
          const double cur_chisq = chisqs[vidx];
          if (cur_chisq >= 0.0) {
            sortbuf[valid_variant_ct].chisq = cur_chisq;
            sortbuf[valid_variant_ct].ln_pval = ln_pvals[vidx];
            sortbuf[valid_variant_ct].variant_uidx = variant_uidx;
            ++valid_variant_ct;
          }
        }
      } else {
        for (uint32_t vidx = 0; vidx < orig_variant_ct; ++vidx, ++variant_uidx) {
          MovU32To1Bit(variant_include, &variant_uidx);
          const double cur_chisq = chisqs[vidx];
          if (cur_chisq >= 0.0) {
            sortbuf[valid_variant_ct].chisq = cur_chisq;
            sortbuf[valid_variant_ct].ln_pval = ChisqToLnP(cur_chisq, 1);
            sortbuf[valid_variant_ct].variant_uidx = variant_uidx;
            ++valid_variant_ct;
          }
        }
      }
    } else {
      uint32_t variant_uidx = 0;
      for (uint32_t vidx = 0; vidx < orig_variant_ct; ++vidx, ++variant_uidx) {
        MovU32To1Bit(variant_include, &variant_uidx);
        const double cur_ln_pval = ln_pvals[vidx];
        if (cur_ln_pval <= 0.0) {
          sortbuf[valid_variant_ct].chisq = LnPToChisq(cur_ln_pval);
          sortbuf[valid_variant_ct].ln_pval = cur_ln_pval;
          sortbuf[valid_variant_ct].variant_uidx = variant_uidx;
          ++valid_variant_ct;
        }
      }
    }
    if (!valid_variant_ct) {
      logputs("Zero valid tests; --adjust skipped.\n");
      goto Multcomp_ret_1;
    }
    BigstackShrinkTop(sortbuf, valid_variant_ct * sizeof(AdjAssocResult));

    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * kMaxIdSlen + 512 + 2 * max_allele_slen;
    const AdjustFlags flags = adjust_info_ptr->flags;
    const uint32_t output_zst = flags & kfAdjustZs;
    OutnameZstSet(".adjusted", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto Multcomp_ret_1;
    }
    *cswritep++ = '#';
    const uint32_t chr_col = flags & kfAdjustColChrom;
    if (chr_col) {
      cswritep = strcpya(cswritep, "CHROM\t");
    }
    if (flags & kfAdjustColPos) {
      cswritep = strcpya(cswritep, "POS\t");
    } else {
      variant_bps = nullptr;
    }
    cswritep = memcpyl3a(cswritep, "ID\t");
    const uint32_t ref_col = flags & kfAdjustColRef;
    if (ref_col) {
      cswritep = strcpya(cswritep, "REF\t");
    }
    const uint32_t alt1_col = flags & kfAdjustColAlt1;
    if (alt1_col) {
      cswritep = strcpya(cswritep, "ALT1\t");
    }
    const uint32_t alt_col = flags & kfAdjustColAlt;
    if (alt_col) {
      cswritep = strcpya(cswritep, "ALT\t");
    }
    const uint32_t is_neglog10 = flags & kfAdjustLog10;
    const uint32_t unadj_col = flags & kfAdjustColUnadj;
    if (unadj_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "UNADJ\t");
    }
    const uint32_t gc_col = (flags & kfAdjustColGc) && (!skip_gc);
    if (gc_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = memcpyl3a(cswritep, "GC\t");
    }
    const uint32_t qq_col = flags & kfAdjustColQq;
    if (qq_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = memcpyl3a(cswritep, "QQ\t");
    }
    const uint32_t bonf_col = flags & kfAdjustColBonf;
    if (bonf_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "BONF\t");
    }
    const uint32_t holm_col = flags & kfAdjustColHolm;
    if (holm_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "HOLM\t");
    }
    const uint32_t sidakss_col = flags & kfAdjustColSidakss;
    if (sidakss_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "SIDAK_SS\t");
    }
    const uint32_t sidaksd_col = flags & kfAdjustColSidaksd;
    if (sidaksd_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "SIDAK_SD\t");
    }
    const uint32_t fdrbh_col = flags & kfAdjustColFdrbh;
    if (fdrbh_col) {
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "FDR_BH\t");
    }
    double* ln_pv_by = nullptr;
    if (flags & kfAdjustColFdrby) {
      if (bigstack_alloc_d(valid_variant_ct, &ln_pv_by)) {
        goto Multcomp_ret_NOMEM;
      }
      if (is_neglog10) {
        cswritep = strcpya(cswritep, "LOG10_");
      }
      cswritep = strcpya(cswritep, "FDR_BY\t");
    }
    DecrAppendBinaryEoln(&cswritep);

    // reverse-order calculations
    double* ln_pv_bh;
    double* ln_pv_gc;
    double* unadj_sorted_ln_pvals;
    if (bigstack_alloc_d(valid_variant_ct, &ln_pv_bh) ||
        bigstack_alloc_d(valid_variant_ct, &ln_pv_gc) ||
        bigstack_alloc_d(valid_variant_ct, &unadj_sorted_ln_pvals)) {
      goto Multcomp_ret_NOMEM;
    }

#ifdef __cplusplus
    std::sort(sortbuf, &(sortbuf[valid_variant_ct]));
#else
    qsort(sortbuf, valid_variant_ct, sizeof(AdjAssocResult), double_cmp);
#endif

    double lambda_recip = 1.0;
    if (!skip_gc) {
      if (adjust_info_ptr->lambda != 0.0) {
        lambda_recip = 1.0 / adjust_info_ptr->lambda;
      } else {
        const uint32_t valid_variant_ct_d2 = valid_variant_ct / 2;
        double lambda = sortbuf[valid_variant_ct_d2].chisq;
        if (!(valid_variant_ct % 2)) {
          lambda = (lambda + sortbuf[valid_variant_ct_d2 - 1].chisq) * 0.5;
        }
        lambda = lambda / 0.456;
        if (lambda < 1.0) {
          lambda = 1.0;
        }
        logprintf("--adjust: Genomic inflation est. lambda (based on median chisq) = %g.\n", lambda);
        lambda_recip = 1.0 / lambda;
      }
    }
    double* sorted_ln_pvals = unadj_sorted_ln_pvals;
    for (uint32_t vidx = 0; vidx < valid_variant_ct; ++vidx) {
      ln_pv_gc[vidx] = ChisqToLnP(sortbuf[vidx].chisq * lambda_recip, 1);
      unadj_sorted_ln_pvals[vidx] = sortbuf[vidx].ln_pval;
    }
    if ((flags & kfAdjustGc) && (!skip_gc)) {
      sorted_ln_pvals = ln_pv_gc;
    }

    const uint32_t valid_variant_ct_m1 = valid_variant_ct - 1;
    const double valid_variant_ctd = u31tod(valid_variant_ct);
    const double ln_valid_variant_ct = log(valid_variant_ctd);
    double bh_ln_pval_min = sorted_ln_pvals[valid_variant_ct_m1];
    ln_pv_bh[valid_variant_ct_m1] = bh_ln_pval_min;
    double harmonic_sum = 1.0;
    for (uint32_t vidx = valid_variant_ct_m1; vidx; --vidx) {
      const double harmonic_term = valid_variant_ctd / u31tod(vidx);
      harmonic_sum += harmonic_term;
      const double bh_ln_pval = sorted_ln_pvals[vidx - 1] + log(harmonic_term);
      if (bh_ln_pval_min > bh_ln_pval) {
        bh_ln_pval_min = bh_ln_pval;
      }
      ln_pv_bh[vidx - 1] = bh_ln_pval_min;
    }

    if (ln_pv_by) {
      const double ln_harmonic_sum = log(harmonic_sum);
      double by_ln_pval_min = sorted_ln_pvals[valid_variant_ct_m1] - ln_valid_variant_ct + ln_harmonic_sum;
      if (by_ln_pval_min > 0.0) {
        by_ln_pval_min = 0.0;
      }
      ln_pv_by[valid_variant_ct_m1] = by_ln_pval_min;
      for (uint32_t vidx = valid_variant_ct_m1; vidx; --vidx) {
        const double by_ln_pval = sorted_ln_pvals[vidx - 1] - log(u31tod(vidx)) + ln_harmonic_sum;
        if (by_ln_pval_min > by_ln_pval) {
          by_ln_pval_min = by_ln_pval;
        }
        ln_pv_by[vidx - 1] = by_ln_pval_min;
      }
    }

    char output_min_p_buf[24];
    uint32_t output_min_p_slen;
    if (!is_neglog10) {
      char* str_end = lntoa_g(output_min_ln, output_min_p_buf);
      output_min_p_slen = str_end - output_min_p_buf;
    } else {
      // -log10(p) output ignores --output-min-p now.
      // instead, set to maximum int32 to distinguish from plink 1.x's
      // much-less-extreme 'inf'.
      strcpy(output_min_p_buf, "2147483647");
      output_min_p_slen = 10;
    }
    const double valid_variant_ct_recip = 1.0 / valid_variant_ctd;
    double ln_pv_sidak_sd = -DBL_MAX;
    double ln_pv_holm = -DBL_MAX;
    uint32_t cur_allele_ct = 2;
    uint32_t vidx = 0;
    for (; vidx < valid_variant_ct; ++vidx) {
      double ln_pval = sorted_ln_pvals[vidx];
      if (ln_pval > ln_pfilter) {
        break;
      }
      const uint32_t variant_uidx = sortbuf[vidx].variant_uidx;
      if (chr_col) {
        if (cip) {
          cswritep = chrtoa(cip, GetVariantChr(cip, variant_uidx), cswritep);
        } else {
          cswritep = strcpya(cswritep, chr_ids[variant_uidx]);
        }
        *cswritep++ = '\t';
      }
      if (variant_bps) {
        cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      }
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
        variant_allele_idx_base = variant_allele_idxs[variant_uidx];
        cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      }
      const char* const* cur_alleles = &(allele_storage[variant_allele_idx_base]);
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
          if (Cswrite(&css, &cswritep)) {
            goto Multcomp_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        --cswritep;
      }
      if (unadj_col) {
        adjust_print_ln(output_min_p_buf, unadj_sorted_ln_pvals[vidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (gc_col) {
        adjust_print_ln(output_min_p_buf, ln_pv_gc[vidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (qq_col) {
        *cswritep++ = '\t';
        double qq_val = (u31tod(vidx) + 0.5) * valid_variant_ct_recip;
        if (is_neglog10) {
          qq_val = -log10(qq_val);
        }
        cswritep = dtoa_g(qq_val, cswritep);
      }
      if (bonf_col) {
        double bonf_ln_pval = ln_pval + ln_valid_variant_ct;
        if (bonf_ln_pval > 0.0) {
          bonf_ln_pval = 0.0;
        }
        adjust_print_ln(output_min_p_buf, bonf_ln_pval, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (holm_col) {
        if (ln_pv_holm < 0.0) {
          const double ln_pv_holm_new = ln_pval + log(u31tod(valid_variant_ct - vidx));
          if (ln_pv_holm_new > 0.0) {
            ln_pv_holm = 0.0;
          } else if (ln_pv_holm < ln_pv_holm_new) {
            ln_pv_holm = ln_pv_holm_new;
          }
        }
        adjust_print_ln(output_min_p_buf, ln_pv_holm, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (sidakss_col) {
        // avoid catastrophic cancellation for small p-values
        // 1 - (1-p)^c = 1 - e^{c log(1-p)}
        // 2^{-7} threshold is arbitrary
        // 2^{-90} corresponds to cp + (cp)^2/2! == cp in double-precision
        // arithmetic, with several bits to spare
        double ln_pv_sidak_ss;
        if (ln_pval > -90 * kLn2) {
          const double pval = exp(ln_pval);
          double pv_sidak_ss;
          if (ln_pval >= -7 * kLn2) {
            pv_sidak_ss = 1 - pow(1 - pval, valid_variant_ctd);
          } else {
            pv_sidak_ss = 1 - exp(valid_variant_ctd * log1p(-pval));
          }
          ln_pv_sidak_ss = log(pv_sidak_ss);
        } else {
          // log(1-x) = -x - x^2/2 - x^3/3 + ...
          // 1 - exp(x) = -x - x^2/2! - x^3/3! - ...
          // if p <= 2^{-90},
          //   log(1-p) is -p
          //   1 - e^{-cp} is cp
          ln_pv_sidak_ss = ln_pval + ln_valid_variant_ct;
        }
        adjust_print_ln(output_min_p_buf, ln_pv_sidak_ss, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (sidaksd_col) {
        double ln_pv_sidak_sd_new;
        if (ln_pval > -90 * kLn2) {
          const double pval = exp(ln_pval);
          double pv_sidak_sd_new;
          if (ln_pval >= -7 * kLn2) {
            pv_sidak_sd_new = 1 - pow(1 - pval, valid_variant_ctd - u31tod(vidx));
          } else {
            const double cur_exp = valid_variant_ctd - u31tod(vidx);
            pv_sidak_sd_new = 1 - exp(cur_exp * log1p(-pval));
          }
          ln_pv_sidak_sd_new = log(pv_sidak_sd_new);
        } else {
          ln_pv_sidak_sd_new = ln_pval + log(valid_variant_ctd - u31tod(vidx));
        }
        if (ln_pv_sidak_sd < ln_pv_sidak_sd_new) {
          ln_pv_sidak_sd = ln_pv_sidak_sd_new;
        }
        adjust_print_ln(output_min_p_buf, ln_pv_sidak_sd, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (fdrbh_col) {
        adjust_print_ln(output_min_p_buf, ln_pv_bh[vidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (ln_pv_by) {
        adjust_print_ln(output_min_p_buf, ln_pv_by[vidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      AppendBinaryEoln(&cswritep);
      if (Cswrite(&css, &cswritep)) {
        goto Multcomp_ret_WRITE_FAIL;
      }
    }
    if (CswriteCloseNull(&css, cswritep)) {
      goto Multcomp_ret_WRITE_FAIL;
    }
    // don't use valid_variant_ct due to --pfilter
    logprintfww("--adjust%s values (%u variant%s) written to %s .\n", cip? "" : "-file", vidx, (vidx == 1)? "" : "s", outname);
  }
  while (0) {
  Multcomp_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Multcomp_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 Multcomp_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr AdjustFile(const AdjustFileInfo* afip, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* in_fname = afip->fname;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  ReadLineStream adjust_rls;
  PreinitRLstream(&adjust_rls);
  {
    // Two-pass load.
    // 1. Parse header line, count # of variants.
    // intermission. Allocate top-level arrays.
    // 2. Rewind and fill arrays.
    // (some overlap with LoadPvar(), though that's one-pass.)
    const char* line_iter;
    reterr = SizeAndInitRLstreamRawK(in_fname, bigstack_left() / 4, &adjust_rls, &line_iter);
    if (reterr) {
      goto AdjustFile_ret_1;
    }

    do {
      ++line_idx;
      reterr = RlsNextLstripK(&adjust_rls, &line_iter);
      if (reterr) {
        if (reterr == kPglRetEof) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", in_fname);
          goto AdjustFile_ret_MALFORMED_INPUT_WW;
        }
        goto AdjustFile_ret_READ_RLSTREAM;
      }
    } while (strequal_k_unsafe(line_iter, "##"));
    const char* linebuf_first_token = line_iter;
    if (*linebuf_first_token == '#') {
      ++linebuf_first_token;
    }

    const AdjustFlags flags = afip->base.flags;
    // [0] = CHROM
    // [1] = POS
    // [2] = ID (required)
    // [3] = REF
    // [4] = ALT
    // [5] = TEST (always scan)
    // [6] = P (required)
    const char* col_search_order[7];
    const uint32_t need_chr = (flags & kfAdjustColChrom);
    const uint32_t need_pos = (flags & kfAdjustColPos);
    const uint32_t need_ref = (flags & kfAdjustColRef);
    const uint32_t need_alt = (flags & (kfAdjustColAlt1 | kfAdjustColAlt));
    const uint32_t alt_comma_truncate = (need_alt == kfAdjustColAlt1);
    if (need_alt == (kfAdjustColAlt1 | kfAdjustColAlt)) {
      // Could theoretically support this later (allocate variant_allele_idxs
      // and count # of multiallelic variants in first pass, etc.), but
      // unlikely to be relevant.
      // (For now, we abuse allele_storage by storing all comma-separated alt
      // alleles in a single string, since Multcomp() is ok with that.)
      logerrputs("Error: --adjust-file does not currently support simultaneous alt1 and alt\ncolumn output.\n");
      goto AdjustFile_ret_INVALID_CMDLINE;
    }
    col_search_order[0] = need_chr? (afip->chr_field? afip->chr_field : "CHROM\0CHR\0") : "";
    col_search_order[1] = need_pos? (afip->pos_field? afip->pos_field : "POS\0BP\0") : "";
    col_search_order[2] = afip->id_field? afip->id_field : "ID\0SNP\0";
    col_search_order[3] = need_ref? (afip->ref_field? afip->ref_field : "REF\0A2\0") : "";
    col_search_order[4] = need_alt? (afip->alt_field? afip->alt_field : "ALT\0ALT1\0A1\0") : "";
    col_search_order[5] = afip->test_field? afip->test_field : "TEST\0";
    const uint32_t input_log10 = (flags & kfAdjustInputLog10);
    col_search_order[6] = afip->p_field? afip->p_field : (input_log10? "LOG10_P\0LOG10_UNADJ\0P\0UNADJ\0" : "P\0UNADJ\0");

    uint32_t col_skips[7];
    uint32_t col_types[7];
    uint32_t relevant_col_ct;
    uint32_t found_type_bitset;
    reterr = SearchHeaderLine(linebuf_first_token, col_search_order, "adjust-file", 7, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
    if (reterr) {
      goto AdjustFile_ret_1;
    }
    if ((found_type_bitset & 0x44) != 0x44) {
      logerrputs("Error: --adjust-file requires ID and P columns.\n");
      goto AdjustFile_ret_INCONSISTENT_INPUT;
    }
    const char* test_name = afip->test_name;
    uint32_t test_name_slen = 0;
    uint32_t test_col_idx = 0;
    if (test_name) {
      test_name_slen = strlen(test_name);
      // this duplicates a bit of work done in SearchHeaderLine(), but not a
      // big deal
      uint32_t relevant_col_idx = 0;
      while (1) {
        test_col_idx += col_skips[relevant_col_idx];
        if (col_types[relevant_col_idx] == 5) {
          break;
        }
        ++relevant_col_idx;
      }
    } else if (found_type_bitset & 0x20) {
      snprintf(g_logbuf, kLogbufSize, "Error: TEST column present in %s, but no test= parameter was provided to --adjust-file.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (need_chr && (!(found_type_bitset & 0x1))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No chromosome column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (need_pos && (!(found_type_bitset & 0x2))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No bp coordinate column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (need_ref && (!(found_type_bitset & 0x8))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No REF column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (need_alt && (!(found_type_bitset & 0x10))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No ALT column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }

    uintptr_t variant_ct = 0;
    while (1) {
      reterr = RlsNextNonemptyLstripK(&adjust_rls, &line_idx, &line_iter);
      if (reterr) {
        if (reterr == kPglRetEof) {
          // reterr = kPglRetSuccess;
          break;
        }
        goto AdjustFile_ret_READ_RLSTREAM;
      }
      if (test_name) {
        // Don't count different-test variants.
        line_iter = NextTokenMult0(line_iter, test_col_idx);
        if (!line_iter) {
          goto AdjustFile_ret_MISSING_TOKENS;
        }
        const char* test_name_start = line_iter;
        const uint32_t cur_test_slen = strlen_se(line_iter);
        line_iter = &(line_iter[cur_test_slen]);
        if ((cur_test_slen != test_name_slen) || memcmp(test_name_start, test_name, test_name_slen)) {
          continue;
        }
      }
      ++variant_ct;
    }
#ifdef __LP64__
    if (variant_ct > 0xffffffffU) {
      logerrputs("Error: Too many variants for --adjust-file.\n");
      goto AdjustFile_ret_INCONSISTENT_INPUT;
    }
#endif

    RewindRLstreamRawK(&adjust_rls, &line_iter);
    const uintptr_t line_ct = line_idx;
    line_idx = 0;
    do {
      ++line_idx;
      reterr = RlsNextLstripK(&adjust_rls, &line_iter);
      if (reterr) {
        goto AdjustFile_ret_READ_FAIL;
      }
    } while (strequal_k_unsafe(line_iter, "##"));
    linebuf_first_token = line_iter;

    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    uintptr_t* variant_include_dummy;
    if (bigstack_alloc_w(variant_ctl, &variant_include_dummy)) {
      goto AdjustFile_ret_NOMEM;
    }
    SetAllBits(variant_ct, variant_include_dummy);
    char** chr_ids;
    if (need_chr) {
      if (bigstack_alloc_cp(variant_ct, &chr_ids)) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      chr_ids = nullptr;
    }
    uint32_t* variant_bps;
    if (need_pos) {
      if (bigstack_alloc_u32(variant_ct, &variant_bps)) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      variant_bps = nullptr;
    }
    char** variant_ids;
    double* ln_pvals;
    if (bigstack_alloc_cp(variant_ct, &variant_ids) ||
        bigstack_alloc_d(variant_ct, &ln_pvals)) {
      goto AdjustFile_ret_NOMEM;
    }
    char** allele_storage;
    if (need_ref || need_alt) {
      if (bigstack_alloc_cp(variant_ct * 2, &allele_storage)) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      allele_storage = nullptr;
    }
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    uint32_t max_allele_slen = 1;
    uintptr_t variant_idx = 0;
    while (line_idx < line_ct) {
      ++line_idx;
      if (RlsNextLstripK(&adjust_rls, &line_iter)) {
        goto AdjustFile_ret_READ_FAIL;
      }
      if (IsEolnKns(*line_iter)) {
        continue;
      }
      const char* token_ptrs[7];
      uint32_t token_slens[7];
      line_iter = TokenLexK0(line_iter, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens);
      if (!line_iter) {
        goto AdjustFile_ret_MISSING_TOKENS;
      }
      if (test_name) {
        if ((token_slens[5] != test_name_slen) || memcmp(token_ptrs[5], test_name, test_name_slen)) {
          continue;
        }
      }
      if (chr_ids) {
        const uint32_t cur_slen = token_slens[0];
        if (cur_slen >= S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base)) {
          goto AdjustFile_ret_NOMEM;
        }
        chr_ids[variant_idx] = R_CAST(char*, tmp_alloc_base);
        memcpyx(tmp_alloc_base, token_ptrs[0], cur_slen, '\0');
        tmp_alloc_base = &(tmp_alloc_base[cur_slen + 1]);
      }
      if (variant_bps) {
        if (ScanUintDefcap(token_ptrs[1], &(variant_bps[variant_idx]))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, in_fname);
          goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
        }
      }
      const uint32_t id_slen = token_slens[2];
      if (id_slen >= S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base)) {
        goto AdjustFile_ret_NOMEM;
      }
      variant_ids[variant_idx] = R_CAST(char*, tmp_alloc_base);
      memcpyx(tmp_alloc_base, token_ptrs[2], id_slen, '\0');
      tmp_alloc_base = &(tmp_alloc_base[id_slen + 1]);
      if (need_ref) {
        const uint32_t cur_slen = token_slens[3];
        if (cur_slen >= S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base)) {
          goto AdjustFile_ret_NOMEM;
        }
        allele_storage[2 * variant_idx] = R_CAST(char*, tmp_alloc_base);
        memcpyx(tmp_alloc_base, token_ptrs[3], cur_slen, '\0');
        tmp_alloc_base = &(tmp_alloc_base[cur_slen + 1]);
      }
      if (need_alt) {
        const char* alt_str = token_ptrs[4];
        uint32_t cur_slen = token_slens[4];
        if (alt_comma_truncate) {
          const char* alt_comma = S_CAST(const char*, memchr(alt_str, ',', cur_slen));
          if (alt_comma) {
            cur_slen = alt_comma - alt_str;
          }
        }
        if (cur_slen >= S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base)) {
          goto AdjustFile_ret_NOMEM;
        }
        allele_storage[2 * variant_idx + 1] = R_CAST(char*, tmp_alloc_base);
        memcpyx(tmp_alloc_base, alt_str, cur_slen, '\0');
        tmp_alloc_base = &(tmp_alloc_base[cur_slen + 1]);
      }
      const char* pval_str = token_ptrs[6];
      double ln_pval;
      if (!input_log10) {
        if (!ScanadvLn(pval_str, &ln_pval)) {
          uint32_t cur_slen;
        AdjustFile_alphabetic_pval:
          cur_slen = token_slens[6];
          if (IsNanStr(pval_str, cur_slen)) {
            ln_pval = kLnPvalError;
          } else if (strequal_k(pval_str, "INF", cur_slen) ||
                     (input_log10 && strequal_k(pval_str, "inf", cur_slen))) {
            // From plink 1.x, could be anything smaller than log(5e-324).
            // Just fill with log(5e-324) for now.
            ln_pval = kLnDenormalMin;
          } else {
            goto AdjustFile_ret_INVALID_PVAL;
          }
        }
      } else {
        double neglog10_pval;
        if (!ScanadvDouble(pval_str, &neglog10_pval)) {
          goto AdjustFile_alphabetic_pval;
        }
        ln_pval = neglog10_pval * (-kLn10);
      }
      if (ln_pval > 0.0) {
        goto AdjustFile_ret_INVALID_PVAL;
      }
      ln_pvals[variant_idx] = ln_pval;
      ++variant_idx;
    }
    BigstackEndReset(bigstack_end_mark);
    BigstackBaseSet(tmp_alloc_base);
    reterr = Multcomp(variant_include_dummy, nullptr, TO_CONSTCPCONSTP(chr_ids), variant_bps, TO_CONSTCPCONSTP(variant_ids), nullptr, TO_CONSTCPCONSTP(allele_storage), &(afip->base), ln_pvals, nullptr, variant_ct, max_allele_slen, ln_pfilter, output_min_ln, 0, max_thread_ct, outname, outname_end);
    if (reterr) {
      goto AdjustFile_ret_1;
    }
  }
  while (0) {
  AdjustFile_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  AdjustFile_ret_READ_RLSTREAM:
    RLstreamErrPrint(in_fname, &adjust_rls, &reterr);
    break;
  AdjustFile_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  AdjustFile_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  AdjustFile_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  AdjustFile_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, in_fname);
  AdjustFile_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  AdjustFile_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  AdjustFile_ret_INVALID_PVAL:
    logerrprintfww("Error: Invalid p-value on line %" PRIuPTR " of %s.\n", line_idx, in_fname);
    reterr = kPglRetInconsistentInput;
    break;
  }
 AdjustFile_ret_1:
  CleanupRLstream(&adjust_rls);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
