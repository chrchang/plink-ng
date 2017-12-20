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

#include "plink2_adjust.h"
#include "plink2_compress_stream.h"
#include "plink2_stats.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void init_adjust(adjust_info_t* adjust_info_ptr) {
  adjust_info_ptr->flags = kfAdjust0;
  adjust_info_ptr->lambda = 0.0;
}

typedef struct adj_assoc_result_struct {
  double chisq;
  double pval;
  uint32_t variant_uidx;
#ifdef __cplusplus
  bool operator<(const struct adj_assoc_result_struct& rhs) const {
    // avoids p-value underflow issue, for what it's worth
    return chisq > rhs.chisq;
  }
#endif
} adj_assoc_result_t;

static inline void adjust_print(const char* output_min_p_str, double pval, double output_min_p, uint32_t output_min_p_slen, uint32_t is_log10, char** bufpp) {
  **bufpp = '\t';
  *bufpp += 1;
  if (pval <= output_min_p) {
    *bufpp = memcpya(*bufpp, output_min_p_str, output_min_p_slen);
  } else {
    if (is_log10) {
      pval = -log10(pval);
    }
    *bufpp = dtoa_g(pval, *bufpp);
  }
}

pglerr_t multcomp(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const adjust_info_t* adjust_info_ptr, double* pvals, double* chisqs, uint32_t orig_variant_ct, uint32_t max_allele_slen, double pfilter, double output_min_p, uint32_t skip_gc, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  compress_stream_state_t css;
  pglerr_t reterr = kPglRetSuccess;
  cswrite_init_null(&css);
  {
    adj_assoc_result_t* sortbuf = (adj_assoc_result_t*)bigstack_alloc(orig_variant_ct * sizeof(adj_assoc_result_t));
    if (!sortbuf) {
      goto multcomp_ret_NOMEM;
    }
    uint32_t valid_variant_ct = 0;
    if (chisqs) {
      uint32_t variant_uidx = 0;
      if (pvals) {
        for (uint32_t vidx = 0; vidx < orig_variant_ct; ++vidx, ++variant_uidx) {
          next_set_unsafe_ck(variant_include, &variant_uidx);
          const double cur_chisq = chisqs[vidx];
          if (cur_chisq >= 0.0) {
            sortbuf[valid_variant_ct].chisq = cur_chisq;
            sortbuf[valid_variant_ct].pval = pvals[vidx];
            sortbuf[valid_variant_ct].variant_uidx = variant_uidx;
            ++valid_variant_ct;
          }
        }
      } else {
        for (uint32_t vidx = 0; vidx < orig_variant_ct; ++vidx, ++variant_uidx) {
          next_set_unsafe_ck(variant_include, &variant_uidx);
          const double cur_chisq = chisqs[vidx];
          if (cur_chisq >= 0.0) {
            sortbuf[valid_variant_ct].chisq = cur_chisq;
            sortbuf[valid_variant_ct].pval = chiprob_p(cur_chisq, 1);
            sortbuf[valid_variant_ct].variant_uidx = variant_uidx;
            ++valid_variant_ct;
          }
        }
      }
    } else {
      uint32_t variant_uidx = 0;
      for (uint32_t vidx = 0; vidx < orig_variant_ct; ++vidx, ++variant_uidx) {
        next_set_unsafe_ck(variant_include, &variant_uidx);
        const double cur_pval = pvals[vidx];
        if (cur_pval >= 0.0) {
          sortbuf[valid_variant_ct].chisq = (cur_pval == 0.0)? kMaxInverseChiprob1df : inverse_chiprob(cur_pval, 1);
          sortbuf[valid_variant_ct].pval = cur_pval;
          sortbuf[valid_variant_ct].variant_uidx = variant_uidx;
          ++valid_variant_ct;
        }
      }
    }
    if (!valid_variant_ct) {
      logprint("Zero valid tests; --adjust skipped.\n");
      goto multcomp_ret_1;
    }
    bigstack_shrink_top(sortbuf, valid_variant_ct * sizeof(adj_assoc_result_t));

    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * kMaxIdSlen + 256 + 2 * max_allele_slen;
    const adjust_flags_t flags = adjust_info_ptr->flags;
    const uint32_t output_zst = flags & kfAdjustZs;
    strcpy(outname_end, output_zst? ".adjusted.zst" : ".adjusted");
    reterr = cswrite_init2(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto multcomp_ret_1;
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
    cswritep = strcpya(cswritep, "ID");
    const uint32_t ref_col = flags & kfAdjustColRef;
    if (ref_col) {
      cswritep = strcpya(cswritep, "\tREF");
    }
    const uint32_t alt1_col = flags & kfAdjustColAlt1;
    if (alt1_col) {
      cswritep = strcpya(cswritep, "\tALT1");
    }
    const uint32_t alt_col = flags & kfAdjustColAlt;
    if (alt_col) {
      cswritep = strcpya(cswritep, "\tALT");
    }
    const uint32_t unadj_col = flags & kfAdjustColUnadj;
    if (unadj_col) {
      cswritep = strcpya(cswritep, "\tUNADJ");
    }
    const uint32_t gc_col = (flags & kfAdjustColGc) && (!skip_gc);
    if (gc_col) {
      cswritep = strcpya(cswritep, "\tGC");
    }
    const uint32_t qq_col = flags & kfAdjustColQq;
    if (qq_col) {
      cswritep = strcpya(cswritep, "\tQQ");
    }
    const uint32_t bonf_col = flags & kfAdjustColBonf;
    if (bonf_col) {
      cswritep = strcpya(cswritep, "\tBONF");
    }
    const uint32_t holm_col = flags & kfAdjustColHolm;
    if (holm_col) {
      cswritep = strcpya(cswritep, "\tHOLM");
    }
    const uint32_t sidakss_col = flags & kfAdjustColSidakss;
    if (sidakss_col) {
      cswritep = strcpya(cswritep, "\tSIDAK_SS");
    }
    const uint32_t sidaksd_col = flags & kfAdjustColSidaksd;
    if (sidaksd_col) {
      cswritep = strcpya(cswritep, "\tSIDAK_SD");
    }
    const uint32_t fdrbh_col = flags & kfAdjustColFdrbh;
    if (fdrbh_col) {
      cswritep = strcpya(cswritep, "\tFDR_BH");
    }
    double* pv_by = nullptr;
    if (flags & kfAdjustColFdrby) {
      if (bigstack_alloc_d(valid_variant_ct, &pv_by)) {
        goto multcomp_ret_NOMEM;
      }
      cswritep = strcpya(cswritep, "\tFDR_BY");
    }
    append_binary_eoln(&cswritep);

    // reverse-order calculations
    double* pv_bh;
    double* pv_gc;
    double* unadj_sorted_pvals;
    if (bigstack_alloc_d(valid_variant_ct, &pv_bh) ||
        bigstack_alloc_d(valid_variant_ct, &pv_gc) ||
        bigstack_alloc_d(valid_variant_ct, &unadj_sorted_pvals)) {
      goto multcomp_ret_NOMEM;
    }

#ifdef __cplusplus
    std::sort(sortbuf, &(sortbuf[valid_variant_ct]));
#else
    qsort(sortbuf, valid_variant_ct, sizeof(adj_assoc_result_t), double_cmp_decr);
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
        LOGPRINTF("--adjust: Genomic inflation est. lambda (based on median chisq) = %g.\n", lambda);
        lambda_recip = 1.0 / lambda;
      }
    }
    double* sorted_pvals = unadj_sorted_pvals;
    for (uint32_t vidx = 0; vidx < valid_variant_ct; ++vidx) {
      pv_gc[vidx] = chiprob_p(sortbuf[vidx].chisq * lambda_recip, 1);
      unadj_sorted_pvals[vidx] = sortbuf[vidx].pval;
    }
    if ((flags & kfAdjustGc) && (!skip_gc)) {
      sorted_pvals = pv_gc;
    }

    const uint32_t valid_variant_ct_m1 = valid_variant_ct - 1;
    const double valid_variant_ctd = (double)((int32_t)valid_variant_ct);
    double bh_pval_min = sorted_pvals[valid_variant_ct_m1];
    pv_bh[valid_variant_ct_m1] = bh_pval_min;
    double harmonic_sum = 1.0;
    for (uint32_t vidx = valid_variant_ct_m1; vidx; --vidx) {
      const double harmonic_term = valid_variant_ctd / ((double)((int32_t)vidx));
      harmonic_sum += harmonic_term;
      const double bh_pval = harmonic_term * sorted_pvals[vidx - 1];
      if (bh_pval_min > bh_pval) {
        bh_pval_min = bh_pval;
      }
      pv_bh[vidx - 1] = bh_pval_min;
    }

    const double valid_variant_ct_recip = 1.0 / valid_variant_ctd;
    if (pv_by) {
      double by_pval_min = harmonic_sum * valid_variant_ct_recip * sorted_pvals[valid_variant_ct_m1];
      if (by_pval_min > 1.0) {
        by_pval_min = 1.0;
      }
      pv_by[valid_variant_ct_m1] = by_pval_min;
      for (uint32_t vidx = valid_variant_ct_m1; vidx; --vidx) {
        double by_pval = (harmonic_sum / ((double)((int32_t)vidx))) * sorted_pvals[vidx - 1];
        if (by_pval_min > by_pval) {
          by_pval_min = by_pval;
        }
        pv_by[vidx - 1] = by_pval_min;
      }
    }

    const uint32_t is_log10 = flags & kfAdjustLog10;
    char output_min_p_buf[16];
    uint32_t output_min_p_slen;
    if (!is_log10) {
      char* str_end = dtoa_g(output_min_p, output_min_p_buf);
      output_min_p_slen = (uintptr_t)(str_end - output_min_p_buf);
    } else if (output_min_p > 0.0) {
      char* str_end = dtoa_g(-log10(output_min_p), output_min_p_buf);
      output_min_p_slen = (uintptr_t)(str_end - output_min_p_buf);
    } else {
      memcpyl3(output_min_p_buf, "inf");
      output_min_p_slen = 3;
    }
    double pv_sidak_sd = 0.0;
    double pv_holm = 0.0;
    uint32_t cur_allele_ct = 2;
    uint32_t vidx = 0;
    for (; vidx < valid_variant_ct; ++vidx) {
      double pval = sorted_pvals[vidx];
      if (pval > pfilter) {
        break;
      }
      const uint32_t variant_uidx = sortbuf[vidx].variant_uidx;
      if (chr_col) {
        cswritep = chr_name_write(cip, get_variant_chr(cip, variant_uidx), cswritep);
        *cswritep++ = '\t';
      }
      if (variant_bps) {
        cswritep = uint32toa_x(variant_bps[variant_uidx], '\t', cswritep);
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
            goto multcomp_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        --cswritep;
      }
      if (unadj_col) {
        adjust_print(output_min_p_buf, unadj_sorted_pvals[vidx], output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (gc_col) {
        adjust_print(output_min_p_buf, pv_gc[vidx], output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (qq_col) {
        *cswritep++ = '\t';
        cswritep = dtoa_g((((double)((int32_t)vidx)) + 0.5) * valid_variant_ct_recip, cswritep);
      }
      if (bonf_col) {
        const double bonf_pval = MINV(pval * valid_variant_ctd, 1.0);
        adjust_print(output_min_p_buf, bonf_pval, output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (holm_col) {
        if (pv_holm < 1.0) {
          const double pv_holm_new = (double)((int32_t)(valid_variant_ct - vidx)) * pval;
          if (pv_holm_new > 1.0) {
            pv_holm = 1.0;
          } else if (pv_holm < pv_holm_new) {
            pv_holm = pv_holm_new;
          }
        }
        adjust_print(output_min_p_buf, pv_holm, output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (sidakss_col) {
        // avoid catastrophic cancellation for small p-values
        // 1 - (1-p)^c = 1 - (1 - cp + (c(c-1) / 2)p^2 - (c(c-1)(c-2) / 6)p^3 +
        //               ...)
        //             = cp - (c(c-1) / 2)p^2 + (c(c-1)(c-2) / 6)p^3 -
        //               [stuff smaller than (c^4p^4)/24]
        // current threshold is chosen to ensure at least 6 digits of precision
        // in (1 - pval) if pow() is used, since 6 significant digits are
        // printed in the .adjusted file.  but in theory we should take
        // valid_variant_ct into account too: small valid_variant_ct lets us
        // use a higher threshold.
        double pv_sidak_ss;
        if (pval >= kRecip2m53 * 2097152) {
          pv_sidak_ss = 1 - pow(1 - pval, valid_variant_ctd);
        } else {
          pv_sidak_ss = pval * valid_variant_ctd * (1 + pval * (valid_variant_ctd - 1) * (-0.5 + pval * (valid_variant_ctd - 2) * (1.0 / 6.0)));
        }
        adjust_print(output_min_p_buf, pv_sidak_ss, output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (sidaksd_col) {
        double pv_sidak_sd_new;
        if (pval >= kRecip2m53 * 2097152) {
          pv_sidak_sd_new = 1 - pow(1 - pval, valid_variant_ctd - ((double)((int32_t)vidx)));
        } else {
          // for very large valid_variant_ct, we might want to include the p^3
          // term of the binomial expansion as well
          const double cur_exp = valid_variant_ctd - (double)((int32_t)vidx);
          pv_sidak_sd_new = pval * cur_exp * (1 + pval * (cur_exp - 1) * (-0.5 + pval * (cur_exp - 2) * (1.0 / 6.0)));
        }
        if (pv_sidak_sd < pv_sidak_sd_new) {
          pv_sidak_sd = pv_sidak_sd_new;
        }
        adjust_print(output_min_p_buf, pv_sidak_sd, output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (fdrbh_col) {
        adjust_print(output_min_p_buf, pv_bh[vidx], output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      if (pv_by) {
        adjust_print(output_min_p_buf, pv_by[vidx], output_min_p, output_min_p_slen, is_log10, &cswritep);
      }
      append_binary_eoln(&cswritep);
      if (cswrite(&css, &cswritep)) {
        goto multcomp_ret_WRITE_FAIL;
      }
    }
    if (cswrite_close_null(&css, cswritep)) {
      goto multcomp_ret_WRITE_FAIL;
    }
    // don't use valid_variant_ct due to --pfilter
    LOGPRINTFWW("--adjust values (%u variant%s) written to %s .\n", vidx, (vidx == 1)? "" : "s", outname);
  }
  while (0) {
  multcomp_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  multcomp_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 multcomp_ret_1:
  cswrite_close_cond(&css, cswritep);
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
