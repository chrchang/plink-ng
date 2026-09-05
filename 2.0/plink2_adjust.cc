// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "include/plink2_float.h"
#include "include/plink2_simd.h"
#include "include/plink2_stats.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"

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
  adjust_file_info_ptr->provref_field = nullptr;
  adjust_file_info_ptr->a1_field = nullptr;
  adjust_file_info_ptr->test_field = nullptr;
  adjust_file_info_ptr->p_field = nullptr;
}

void CleanupAdjust(AdjustFileInfo* adjust_file_info_ptr) {
  free_cond(adjust_file_info_ptr->a1_field);
  free_cond(adjust_file_info_ptr->alt_field);
  free_cond(adjust_file_info_ptr->chr_field);
  if (adjust_file_info_ptr->fname) {
    free(adjust_file_info_ptr->fname);
    free_cond(adjust_file_info_ptr->pos_field);
    free_cond(adjust_file_info_ptr->id_field);
    free_cond(adjust_file_info_ptr->ref_field);
    free_cond(adjust_file_info_ptr->provref_field);
    free_cond(adjust_file_info_ptr->test_field);
    free_cond(adjust_file_info_ptr->p_field);
  }
}

typedef struct AdjAssocResultStruct {
  double ln_pval;
  double chisq;  // do we really need this?...
  uint32_t variant_uidx;
  uint32_t allele_idx;
#ifdef __cplusplus
  bool operator<(const struct AdjAssocResultStruct& rhs) const {
    if (ln_pval != rhs.ln_pval) {
      return ln_pval < rhs.ln_pval;
    }
    // update (11 Mar 2026): make pval tie-handling deterministic
    if (variant_uidx != rhs.variant_uidx) {
      return variant_uidx < rhs.variant_uidx;
    }
    return allele_idx < rhs.allele_idx;
  }
#endif
} AdjAssocResult;

#ifndef __cplusplus
int32_t AdjAssocCmp(const void* aa, const void* bb) {
  const AdjAssocResult* aar1 = S_CAST(const AdjAssocResult*, aa);
  const AdjAssocResult* aar2 = S_CAST(const AdjAssocResult*, bb);
  const double ln_pval1 = aar1->ln_pval;
  const double ln_pval2 = aar2->ln_pval;
  if (ln_pval1 != ln_pval2) {
    return (ln_pval1 < ln_pval2)? -1 : 1;
  }
  const uint32_t variant_uidx1 = aar1->variant_uidx;
  const uint32_t variant_uidx2 = aar2->variant_uidx;
  if (variant_uidx1 != variant_uidx2) {
    return (variant_uidx1 < variant_uidx2)? -1 : 1;
  }
  return S_CAST(int32_t, aar1->allele_idx) - S_CAST(int32_t, aar2->allele_idx);
}
#endif

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

// Now based around ln_pvals, to allow useful comparisons < 2.23e-308.
PglErr Multcomp(const uintptr_t* variant_include, const ChrInfo* cip, const char* const* chr_ids, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_include, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const char* const* loaded_a1, const AdjustInfo* adjust_info_ptr, const double* ln_pvals, const double* chisqs, uint32_t raw_variant_ct, uintptr_t orig_allele_ct, uint32_t max_allele_slen, PgenGlobalFlags gflags, double ln_pfilter, double output_min_ln, uint32_t skip_gc, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    AdjAssocResult* sortbuf;
    if (unlikely(BIGSTACK_ALLOC_X(AdjAssocResult, orig_allele_ct, &sortbuf))) {
      goto Multcomp_ret_NOMEM;
    }
    uintptr_t valid_allele_ct = 0;
    uintptr_t allele_uidx_base = 0;
    uintptr_t allele_include_bits = allele_include[0];
    if (!allele_idx_offsets) {
      if (chisqs) {
        if (ln_pvals) {
          for (uintptr_t aidx = 0; aidx != orig_allele_ct; ++aidx) {
            const uintptr_t allele_uidx = BitIter1(allele_include, &allele_uidx_base, &allele_include_bits);
            const double cur_chisq = chisqs[aidx];
            if (cur_chisq >= 0.0) {
              sortbuf[valid_allele_ct].chisq = cur_chisq;
              sortbuf[valid_allele_ct].ln_pval = ln_pvals[aidx];
              sortbuf[valid_allele_ct].variant_uidx = allele_uidx / 2;
              sortbuf[valid_allele_ct].allele_idx = allele_uidx % 2;
              ++valid_allele_ct;
            }
          }
        } else {
          for (uintptr_t aidx = 0; aidx != orig_allele_ct; ++aidx) {
            const uintptr_t allele_uidx = BitIter1(allele_include, &allele_uidx_base, &allele_include_bits);
            const double cur_chisq = chisqs[aidx];
            if (cur_chisq >= 0.0) {
              sortbuf[valid_allele_ct].chisq = cur_chisq;
              sortbuf[valid_allele_ct].ln_pval = ChisqToLnP(cur_chisq, 1);
              sortbuf[valid_allele_ct].variant_uidx = allele_uidx / 2;
              sortbuf[valid_allele_ct].allele_idx = allele_uidx % 2;
              ++valid_allele_ct;
            }
          }
        }
      } else {
        for (uintptr_t aidx = 0; aidx != orig_allele_ct; ++aidx) {
          const uintptr_t allele_uidx = BitIter1(allele_include, &allele_uidx_base, &allele_include_bits);
          const double cur_ln_pval = ln_pvals[aidx];
          // In --adjust-file case, possible for cur_ln_pval == kLnPvalError
          // (which is intentionally positive) when allele_include bit set.
          // (Don't think there are any other cases yet where valid_allele_ct
          // can be less than orig_allele_ct?)
          if (cur_ln_pval <= 0.0) {
            sortbuf[valid_allele_ct].chisq = LnPToChisq(cur_ln_pval);
            sortbuf[valid_allele_ct].ln_pval = cur_ln_pval;
            sortbuf[valid_allele_ct].variant_uidx = allele_uidx / 2;
            sortbuf[valid_allele_ct].allele_idx = allele_uidx % 2;
            ++valid_allele_ct;
          }
        }
      }
    } else {
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      uintptr_t allele_idx_offset_start = allele_idx_offsets[variant_uidx];
      uintptr_t allele_idx_offset_end = allele_idx_offsets[variant_uidx + 1];
      if (chisqs) {
        if (ln_pvals) {
          for (uintptr_t aidx = 0; aidx != orig_allele_ct; ++aidx) {
            const uintptr_t allele_uidx = BitIter1(allele_include, &allele_uidx_base, &allele_include_bits);
            if (allele_uidx >= allele_idx_offset_end) {
              variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
              allele_idx_offset_start = allele_idx_offsets[variant_uidx];
              allele_idx_offset_end = allele_idx_offsets[variant_uidx + 1];
            }
            const double cur_chisq = chisqs[aidx];
            if (cur_chisq >= 0.0) {
              sortbuf[valid_allele_ct].chisq = cur_chisq;
              sortbuf[valid_allele_ct].ln_pval = ln_pvals[aidx];
              sortbuf[valid_allele_ct].variant_uidx = variant_uidx;
              sortbuf[valid_allele_ct].allele_idx = allele_uidx - allele_idx_offset_start;
              ++valid_allele_ct;
            }
          }
        } else {
          for (uintptr_t aidx = 0; aidx != orig_allele_ct; ++aidx) {
            const uintptr_t allele_uidx = BitIter1(allele_include, &allele_uidx_base, &allele_include_bits);
            if (allele_uidx >= allele_idx_offset_end) {
              variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
              allele_idx_offset_start = allele_idx_offsets[variant_uidx];
              allele_idx_offset_end = allele_idx_offsets[variant_uidx + 1];
            }
            const double cur_chisq = chisqs[aidx];
            if (cur_chisq >= 0.0) {
              sortbuf[valid_allele_ct].chisq = cur_chisq;
              sortbuf[valid_allele_ct].ln_pval = ChisqToLnP(cur_chisq, 1);
              sortbuf[valid_allele_ct].variant_uidx = variant_uidx;
              sortbuf[valid_allele_ct].allele_idx = allele_uidx - allele_idx_offset_start;
              ++valid_allele_ct;
            }
          }
        }
      } else {
        for (uintptr_t aidx = 0; aidx != orig_allele_ct; ++aidx) {
          const uintptr_t allele_uidx = BitIter1(allele_include, &allele_uidx_base, &allele_include_bits);
          if (allele_uidx >= allele_idx_offset_end) {
            variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
            allele_idx_offset_start = allele_idx_offsets[variant_uidx];
            allele_idx_offset_end = allele_idx_offsets[variant_uidx + 1];
          }
          const double cur_ln_pval = ln_pvals[aidx];
          if (cur_ln_pval <= 0.0) {
            sortbuf[valid_allele_ct].chisq = LnPToChisq(cur_ln_pval);
            sortbuf[valid_allele_ct].ln_pval = cur_ln_pval;
            sortbuf[valid_allele_ct].variant_uidx = variant_uidx;
            sortbuf[valid_allele_ct].allele_idx = allele_uidx - allele_idx_offset_start;
            ++valid_allele_ct;
          }
        }
      }
    }
    if (!valid_allele_ct) {
      logputs("Zero valid tests; --adjust skipped.\n");
      goto Multcomp_ret_1;
    }
    BigstackShrinkTop(sortbuf, valid_allele_ct * sizeof(AdjAssocResult));

    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * kMaxIdSlen + 512 + 2 * max_allele_slen;
    const AdjustFlags flags = adjust_info_ptr->flags;
    const uint32_t output_zst = flags & kfAdjustZs;
    OutnameZstSet(".adjusted", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto Multcomp_ret_1;
    }
    *cswritep++ = '#';
    const uint32_t chr_col = flags & kfAdjustColChrom;
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (flags & kfAdjustColPos) {
      cswritep = strcpya_k(cswritep, "POS\t");
    } else {
      variant_bps = nullptr;
    }
    cswritep = strcpya_k(cswritep, "ID\t");
    const uint32_t ref_col = flags & kfAdjustColRef;
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "REF\t");
    }
    const uint32_t alt1_col = flags & kfAdjustColAlt1;
    if (alt1_col) {
      cswritep = strcpya_k(cswritep, "ALT1\t");
    }
    const uint32_t alt_col = flags & kfAdjustColAlt;
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "ALT\t");
    }
    const uint32_t all_nonref = (gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, flags / kfAdjustColMaybeprovref, raw_variant_ct, all_nonref);
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "PROVISIONAL_REF?\t");
    }
    const uint32_t a1_col = (flags & kfAdjustColA1) && (loaded_a1 || cip);
    if (a1_col) {
      cswritep = strcpya_k(cswritep, "A1\t");
    }
    const uint32_t is_neglog10 = flags & kfAdjustLog10;
    const uint32_t unadj_col = flags & kfAdjustColUnadj;
    if (unadj_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "UNADJ\t");
    }
    const uint32_t gc_col = (flags & kfAdjustColGc) && (!skip_gc);
    if (gc_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "GC\t");
    }
    const uint32_t qq_col = flags & kfAdjustColQq;
    if (qq_col) {
      cswritep = strcpya_k(cswritep, "QQ\t");
    }
    const uint32_t bonf_col = flags & kfAdjustColBonf;
    if (bonf_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "BONF\t");
    }
    const uint32_t holm_col = flags & kfAdjustColHolm;
    if (holm_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "HOLM\t");
    }
    const uint32_t sidakss_col = flags & kfAdjustColSidakss;
    if (sidakss_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "SIDAK_SS\t");
    }
    const uint32_t sidaksd_col = flags & kfAdjustColSidaksd;
    if (sidaksd_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "SIDAK_SD\t");
    }
    const uint32_t fdrbh_col = flags & kfAdjustColFdrbh;
    if (fdrbh_col) {
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "FDR_BH\t");
    }
    double* ln_pv_by = nullptr;
    if (flags & kfAdjustColFdrby) {
      if (unlikely(bigstack_alloc_d(valid_allele_ct, &ln_pv_by))) {
        goto Multcomp_ret_NOMEM;
      }
      if (is_neglog10) {
        cswritep = strcpya_k(cswritep, "NEG_LOG10_");
      }
      cswritep = strcpya_k(cswritep, "FDR_BY\t");
    }
    DecrAppendBinaryEoln(&cswritep);

    // reverse-order calculations
    double* ln_pv_bh;
    double* ln_pv_gc;
    double* unadj_sorted_ln_pvals;
    if (unlikely(bigstack_alloc_d(valid_allele_ct, &ln_pv_bh) ||
                 bigstack_alloc_d(valid_allele_ct, &ln_pv_gc) ||
                 bigstack_alloc_d(valid_allele_ct, &unadj_sorted_ln_pvals))) {
      goto Multcomp_ret_NOMEM;
    }

    STD_SORT_PAR_UNSEQ(valid_allele_ct, AdjAssocCmp, sortbuf);

    double lambda_recip = 1.0;
    if (!skip_gc) {
      if (adjust_info_ptr->lambda != 0.0) {
        lambda_recip = 1.0 / adjust_info_ptr->lambda;
      } else {
        const uintptr_t valid_allele_ct_d2 = valid_allele_ct / 2;
        double lambda = sortbuf[valid_allele_ct_d2].chisq;
        if (!(valid_allele_ct % 2)) {
          lambda = (lambda + sortbuf[valid_allele_ct_d2 - 1].chisq) * 0.5;
        }
        lambda = lambda / 0.456;
        logprintf("--adjust: Genomic inflation est. lambda (based on median chisq) = %g.\n", lambda);
        if (lambda < 1.0) {
          logprintf("(Treating lambda as 1 in GC-corrected p-value calculation.)\n");
          lambda = 1.0;
        }
        lambda_recip = 1.0 / lambda;
      }
    }
    double* sorted_ln_pvals = unadj_sorted_ln_pvals;
    for (uintptr_t aidx = 0; aidx != valid_allele_ct; ++aidx) {
      ln_pv_gc[aidx] = ChisqToLnP(sortbuf[aidx].chisq * lambda_recip, 1);
      unadj_sorted_ln_pvals[aidx] = sortbuf[aidx].ln_pval;
    }
    if ((flags & kfAdjustGc) && (!skip_gc)) {
      sorted_ln_pvals = ln_pv_gc;
    }

    const uint32_t valid_allele_ct_m1 = valid_allele_ct - 1;
    const double valid_allele_ctd = swtod(valid_allele_ct);
    const double ln_valid_allele_ct = log(valid_allele_ctd);
    double bh_ln_pval_min = sorted_ln_pvals[valid_allele_ct_m1];
    ln_pv_bh[valid_allele_ct_m1] = bh_ln_pval_min;
    double harmonic_sum = 1.0;
    for (uint32_t aidx = valid_allele_ct_m1; aidx; --aidx) {
      const double harmonic_term = valid_allele_ctd / u31tod(aidx);
      harmonic_sum += harmonic_term;
      const double bh_ln_pval = sorted_ln_pvals[aidx - 1] + log(harmonic_term);
      if (bh_ln_pval_min > bh_ln_pval) {
        bh_ln_pval_min = bh_ln_pval;
      }
      ln_pv_bh[aidx - 1] = bh_ln_pval_min;
    }

    if (ln_pv_by) {
      const double ln_harmonic_sum = log(harmonic_sum);
      double by_ln_pval_min = sorted_ln_pvals[valid_allele_ct_m1] - ln_valid_allele_ct + ln_harmonic_sum;
      if (by_ln_pval_min > 0.0) {
        by_ln_pval_min = 0.0;
      }
      ln_pv_by[valid_allele_ct_m1] = by_ln_pval_min;
      for (uint32_t aidx = valid_allele_ct_m1; aidx; --aidx) {
        const double by_ln_pval = sorted_ln_pvals[aidx - 1] - log(u31tod(aidx)) + ln_harmonic_sum;
        if (by_ln_pval_min > by_ln_pval) {
          by_ln_pval_min = by_ln_pval;
        }
        ln_pv_by[aidx - 1] = by_ln_pval_min;
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
      strcpy_k(output_min_p_buf, "2147483647");
      output_min_p_slen = 10;
    }
    const double valid_allele_ct_recip = 1.0 / valid_allele_ctd;
    double ln_pv_sidak_sd = -DBL_MAX;
    double ln_pv_holm = -DBL_MAX;
    uint32_t cur_allele_ct = 2;
    uint32_t aidx = 0;
    for (; aidx < valid_allele_ct; ++aidx) {
      double ln_pval = sorted_ln_pvals[aidx];
      if (ln_pval > ln_pfilter) {
        break;
      }
      const uint32_t variant_uidx = sortbuf[aidx].variant_uidx;
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
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
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
        for (uint32_t allele_idx = 1; allele_idx != cur_allele_ct; ++allele_idx) {
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto Multcomp_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        --cswritep;
      }
      if (provref_col) {
        *cswritep++ = '\t';
        *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
      }
      if (a1_col) {
        *cswritep++ = '\t';
        const char* cur_allele;
        if (loaded_a1) {
          // --adjust-file hack
          cur_allele = loaded_a1[variant_uidx];
        } else {
          cur_allele = cur_alleles[sortbuf[aidx].allele_idx];
        }
        cswritep = strcpya(cswritep, cur_allele);
      }
      if (unadj_col) {
        adjust_print_ln(output_min_p_buf, unadj_sorted_ln_pvals[aidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (gc_col) {
        adjust_print_ln(output_min_p_buf, ln_pv_gc[aidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      const double aidx_d = swtod(aidx);
      if (qq_col) {
        *cswritep++ = '\t';
        double qq_val = (aidx_d + 0.5) * valid_allele_ct_recip;
        cswritep = dtoa_g(qq_val, cswritep);
      }
      if (bonf_col) {
        double bonf_ln_pval = ln_pval + ln_valid_allele_ct;
        if (bonf_ln_pval > 0.0) {
          bonf_ln_pval = 0.0;
        }
        adjust_print_ln(output_min_p_buf, bonf_ln_pval, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (holm_col) {
        if (ln_pv_holm < 0.0) {
          const double ln_pv_holm_new = ln_pval + log(u31tod(valid_allele_ct - aidx));
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
            pv_sidak_ss = 1 - pow(1 - pval, valid_allele_ctd);
          } else {
            pv_sidak_ss = 1 - exp(valid_allele_ctd * log1p(-pval));
          }
          ln_pv_sidak_ss = log(pv_sidak_ss);
        } else {
          // log(1-x) = -x - x^2/2 - x^3/3 + ...
          // 1 - exp(x) = -x - x^2/2! - x^3/3! - ...
          // if p <= 2^{-90},
          //   log(1-p) is -p
          //   1 - e^{-cp} is cp
          ln_pv_sidak_ss = ln_pval + ln_valid_allele_ct;
        }
        adjust_print_ln(output_min_p_buf, ln_pv_sidak_ss, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (sidaksd_col) {
        double ln_pv_sidak_sd_new;
        if (ln_pval > -90 * kLn2) {
          const double pval = exp(ln_pval);
          double pv_sidak_sd_new;
          if (ln_pval >= -7 * kLn2) {
            pv_sidak_sd_new = 1 - pow(1 - pval, valid_allele_ctd - aidx_d);
          } else {
            const double cur_exp = valid_allele_ctd - aidx_d;
            pv_sidak_sd_new = 1 - exp(cur_exp * log1p(-pval));
          }
          ln_pv_sidak_sd_new = log(pv_sidak_sd_new);
        } else {
          ln_pv_sidak_sd_new = ln_pval + log(valid_allele_ctd - aidx_d);
        }
        if (ln_pv_sidak_sd < ln_pv_sidak_sd_new) {
          ln_pv_sidak_sd = ln_pv_sidak_sd_new;
        }
        adjust_print_ln(output_min_p_buf, ln_pv_sidak_sd, output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (fdrbh_col) {
        adjust_print_ln(output_min_p_buf, ln_pv_bh[aidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      if (ln_pv_by) {
        adjust_print_ln(output_min_p_buf, ln_pv_by[aidx], output_min_ln, output_min_p_slen, is_neglog10, &cswritep);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto Multcomp_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto Multcomp_ret_WRITE_FAIL;
    }
    // don't use valid_allele_ct due to --pfilter
    logprintfww("--adjust%s values (%" PRIuPTR " test%s) written to %s .\n", cip? "" : "-file", aidx, (aidx == 1)? "" : "s", outname);
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

void InitMeta(MetaInfo* mip) {
  mip->fnames = nullptr;
  mip->chr_field = nullptr;
  mip->snp_field = nullptr;
  mip->bp_field = nullptr;
  mip->a1_field = nullptr;
  mip->a2_field = nullptr;
  mip->p_field = nullptr;
  mip->se_field = nullptr;
  mip->ess_field = nullptr;
  mip->flags = kfMeta0;
}

void CleanupMeta(MetaInfo* mip) {
  free_cond(mip->fnames);
  free_cond(mip->chr_field);
  free_cond(mip->snp_field);
  free_cond(mip->bp_field);
  free_cond(mip->a1_field);
  free_cond(mip->a2_field);
  free_cond(mip->p_field);
  free_cond(mip->se_field);
  free_cond(mip->ess_field);
}

// --meta-analysis: inverse-variance meta-analysis over several variant-based
// association reports.  Rows are matched by variant ID; when an A2 column is
// present an allele flip is corrected rather than discarded, as PLINK 1.x
// does.
//
// Each file is read twice, once to size the record array and string arena
// exactly and once to fill them.  Growing either at an unknown rate is what
// the bigstack allocator has no good answer for.
typedef struct MetaRecordStruct {
  double beta;
  double se;
  double ln_pval;   // only filled for 'weighted-z'
  double ess;       // effective sample size, ditto
  const char* id;
  const char* a1;  // nullptr when alleles are unused
  const char* a2;  // nullptr when the file has no A2 column
  const char* chr;  // nullptr when 'no-map'
  uint32_t bp;
  uint32_t file_idx;
#ifdef __cplusplus
  bool operator<(const struct MetaRecordStruct& rhs) const {
    const int32_t id_cmp = strcmp(id, rhs.id);
    if (id_cmp) {
      return (id_cmp < 0);
    }
    return (file_idx < rhs.file_idx);
  }
#endif
} MetaRecord;

// ID-major, so same-variant rows are adjacent; file index breaks ties so a
// group's order is deterministic.
// Chromosomes are compared numerically when both sides are numeric, so the
// usual 1..22 ordering falls out without needing a finalized ChrInfo, which a
// standalone batch job does not have.
static int32_t MetaChrCmp(const char* chr1, const char* chr2) {
  if ((!chr1) || (!chr2)) {
    return 0;
  }
  uint32_t val1;
  uint32_t val2;
  const uint32_t numeric1 = !ScanUintDefcapx(chr1, &val1);
  const uint32_t numeric2 = !ScanUintDefcapx(chr2, &val2);
  if (numeric1 && numeric2) {
    if (val1 != val2) {
      return (val1 < val2)? -1 : 1;
    }
    return 0;
  }
  if (numeric1 != numeric2) {
    return numeric1? -1 : 1;
  }
  return strcmp(chr1, chr2);
}

typedef struct MetaGroupStruct {
  const char* chr;
  uint32_t bp;
  uint32_t first_rec;
  uint32_t rec_ct;
#ifdef __cplusplus
  bool operator<(const struct MetaGroupStruct& rhs) const {
    const int32_t chr_cmp = MetaChrCmp(chr, rhs.chr);
    if (chr_cmp) {
      return (chr_cmp < 0);
    }
    if (bp != rhs.bp) {
      return (bp < rhs.bp);
    }
    return (first_rec < rhs.first_rec);
  }
#endif
} MetaGroup;

// Fills col_skips/col_types for one input file; shared by both passes so the
// two cannot disagree about which columns they are reading.
static PglErr MetaSetupCols(const MetaInfo* mip, const char* header_start, const char* fname, uint32_t no_map, uint32_t no_allele, uint32_t input_is_beta, uint32_t weighted_z, uint32_t* col_skips, uint32_t* col_types, uint32_t* relevant_col_ct_ptr, uint32_t* has_a2_ptr) {
  // [0] CHR, [1] SNP, [2] BP, [3] A1, [4] A2, [5] effect, [6] SE, [7] P,
  // [8] effective sample size
  const char* col_search_order[9];
  col_search_order[0] = no_map? "" : (mip->chr_field? mip->chr_field : "CHR\0CHROM\0");
  col_search_order[1] = mip->snp_field? mip->snp_field : "SNP\0ID\0";
  col_search_order[2] = no_map? "" : (mip->bp_field? mip->bp_field : "BP\0POS\0");
  col_search_order[3] = no_allele? "" : (mip->a1_field? mip->a1_field : "A1\0");
  col_search_order[4] = no_allele? "" : (mip->a2_field? mip->a2_field : "A2\0");
  col_search_order[5] = input_is_beta? "BETA\0" : "OR\0";
  col_search_order[6] = mip->se_field? mip->se_field : "SE\0";
  col_search_order[7] = weighted_z? (mip->p_field? mip->p_field : "P\0") : "";
  col_search_order[8] = weighted_z? (mip->ess_field? mip->ess_field : "NMISS\0OBS_CT\0") : "";
  uint32_t found_type_bitset;
  PglErr reterr = SearchHeaderLine(header_start, col_search_order, "--meta-analysis", 9, relevant_col_ct_ptr, &found_type_bitset, col_skips, col_types);
  if (unlikely(reterr)) {
    return reterr;
  }
  if (unlikely(!(found_type_bitset & 0x2))) {
    snprintf(g_logbuf, kLogbufSize, "Error: No variant ID field found in %s.\n", fname);
    return kPglRetInconsistentInput;
  }
  if (unlikely(!(found_type_bitset & 0x20))) {
    snprintf(g_logbuf, kLogbufSize, "Error: No %s field found in %s.\n", input_is_beta? "BETA" : "OR", fname);
    return kPglRetInconsistentInput;
  }
  if (unlikely(!(found_type_bitset & 0x40))) {
    snprintf(g_logbuf, kLogbufSize, "Error: No SE field found in %s.\n", fname);
    return kPglRetInconsistentInput;
  }
  if (unlikely((!no_map) && ((found_type_bitset & 0x5) != 0x5))) {
    snprintf(g_logbuf, kLogbufSize, "Error: No chromosome and/or position field found in %s.  ('no-map' ignores them.)\n", fname);
    return kPglRetInconsistentInput;
  }
  if (unlikely((!no_allele) && (!(found_type_bitset & 0x8)))) {
    snprintf(g_logbuf, kLogbufSize, "Error: No A1 field found in %s.  ('no-allele' ignores it.)\n", fname);
    return kPglRetInconsistentInput;
  }
  if (unlikely(weighted_z && ((found_type_bitset & 0x180) != 0x180))) {
    snprintf(g_logbuf, kLogbufSize, "Error: 'weighted-z' requires p-value and effective sample size fields, and %s is\nmissing at least one.\n", fname);
    return kPglRetInconsistentInput;
  }
  *has_a2_ptr = (found_type_bitset >> 4) & 1;
  return kPglRetSuccess;
}

PglErr MetaAnalysis(const MetaInfo* mip, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PreinitCstream(&css);
  TextStream txs;
  PreinitTextStream(&txs);
  uintptr_t line_idx = 0;
  const char* cur_fname = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const MetaFlags flags = mip->flags;
    const uint32_t qt_mode = (flags / kfMetaQt) & 1;
    const uint32_t no_map = (flags / kfMetaNoMap) & 1;
    const uint32_t no_allele = no_map || ((flags / kfMetaNoAllele) & 1);
    const uint32_t study_col = (flags / kfMetaStudy) & 1;
    const uint32_t report_all = (flags / kfMetaReportAll) & 1;
    const uint32_t output_zst = (flags / kfMetaZs) & 1;
    const uint32_t weighted_z = (flags / kfMetaWeightedZ) & 1;
    // Effect sizes are betas internally; an OR column is logged on the way in.
    const uint32_t input_is_beta = ((flags / kfMetaLogscale) & 1) || qt_mode;

    uint32_t file_ct = 0;
    {
      const char* fname_iter = mip->fnames;
      while (*fname_iter) {
        ++file_ct;
        fname_iter = strnul(fname_iter);
        ++fname_iter;
      }
    }
    if (unlikely(file_ct < 2)) {
      logerrputs("Error: --meta-analysis requires at least two files.\n");
      goto MetaAnalysis_ret_INCONSISTENT_INPUT;
    }

    uintptr_t record_ct = 0;
    uintptr_t str_byte_ct = 0;
    MetaRecord* records = nullptr;
    char* str_arena = nullptr;
    char* str_iter = nullptr;
    for (uint32_t pass_idx = 0; pass_idx != 2; ++pass_idx) {
      if (pass_idx == 1) {
        if (unlikely(!record_ct)) {
          logerrputs("Error: --meta-analysis: No usable rows.\n");
          goto MetaAnalysis_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(BIGSTACK_ALLOC_X(MetaRecord, record_ct, &records) ||
                     bigstack_alloc_c(str_byte_ct, &str_arena))) {
          goto MetaAnalysis_ret_NOMEM;
        }
        str_iter = str_arena;
        record_ct = 0;
      }
      const char* fname_iter = mip->fnames;
      for (uint32_t file_idx = 0; file_idx != file_ct; ++file_idx) {
        cur_fname = fname_iter;
        line_idx = 0;
        reterr = SizeAndInitTextStream(cur_fname, bigstack_left() / 8, max_thread_ct, &txs);
        if (unlikely(reterr)) {
          goto MetaAnalysis_ret_TSTREAM_FAIL;
        }
        const char* header_start;
        do {
          ++line_idx;
          header_start = TextGet(&txs);
          if (unlikely(!header_start)) {
            reterr = TextStreamRawErrcode(&txs);
            if (reterr == kPglRetEof) {
              snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", cur_fname);
              goto MetaAnalysis_ret_MALFORMED_INPUT_WW;
            }
            goto MetaAnalysis_ret_TSTREAM_FAIL;
          }
        } while (strequal_k_unsafe(header_start, "##"));
        if (*header_start == '#') {
          ++header_start;
        }
        uint32_t col_skips[9];
        uint32_t col_types[9];
        uint32_t relevant_col_ct;
        uint32_t has_a2;
        reterr = MetaSetupCols(mip, header_start, cur_fname, no_map, no_allele, input_is_beta, weighted_z, col_skips, col_types, &relevant_col_ct, &has_a2);
        if (unlikely(reterr)) {
          if (reterr == kPglRetInconsistentInput) {
            goto MetaAnalysis_ret_INCONSISTENT_INPUT_WW;
          }
          goto MetaAnalysis_ret_1;
        }
        while (1) {
          ++line_idx;
          const char* line_start = TextGet(&txs);
          if (!line_start) {
            break;
          }
          const char* token_ptrs[9];
          uint32_t token_slens[9];
          if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
            goto MetaAnalysis_ret_MISSING_TOKENS;
          }
          double effect;
          if (!ScantokDouble(token_ptrs[5], &effect)) {
            continue;  // NA rows are skipped, as in PLINK 1.x
          }
          double se;
          if ((!ScantokDouble(token_ptrs[6], &se)) || (se <= 0.0)) {
            continue;
          }
          double cur_ln_pval = 0.0;
          double cur_ess = 0.0;
          if (weighted_z) {
            double cur_pval;
            if ((!ScantokDouble(token_ptrs[7], &cur_pval)) || (cur_pval <= 0.0) || (cur_pval > 1.0)) {
              continue;
            }
            cur_ln_pval = log(cur_pval);
            if ((!ScantokDouble(token_ptrs[8], &cur_ess)) || (cur_ess <= 0.0)) {
              continue;
            }
          }
          if (!input_is_beta) {
            if (effect <= 0.0) {
              continue;
            }
            effect = log(effect);
          }
          uint32_t bp = 0;
          if (!no_map) {
            if (ScanUintDefcapx(token_ptrs[2], &bp)) {
              continue;
            }
          }
          const uint32_t id_slen = token_slens[1];
          uintptr_t cur_str_bytes = id_slen + 1;
          if (!no_map) {
            cur_str_bytes += token_slens[0] + 1;
          }
          if (!no_allele) {
            cur_str_bytes += token_slens[3] + 1;
            if (has_a2) {
              cur_str_bytes += token_slens[4] + 1;
            }
          }
          if (!pass_idx) {
            str_byte_ct += cur_str_bytes;
          } else {
            MetaRecord* cur_rec = &(records[record_ct]);
            cur_rec->beta = effect;
            cur_rec->se = se;
            cur_rec->ln_pval = cur_ln_pval;
            cur_rec->ess = cur_ess;
            cur_rec->bp = bp;
            cur_rec->file_idx = file_idx;
            cur_rec->id = str_iter;
            str_iter = memcpyax(str_iter, token_ptrs[1], id_slen, '\0');
            if (no_map) {
              cur_rec->chr = nullptr;
            } else {
              cur_rec->chr = str_iter;
              str_iter = memcpyax(str_iter, token_ptrs[0], token_slens[0], '\0');
            }
            if (no_allele) {
              cur_rec->a1 = nullptr;
              cur_rec->a2 = nullptr;
            } else {
              cur_rec->a1 = str_iter;
              str_iter = memcpyax(str_iter, token_ptrs[3], token_slens[3], '\0');
              if (has_a2) {
                cur_rec->a2 = str_iter;
                str_iter = memcpyax(str_iter, token_ptrs[4], token_slens[4], '\0');
              } else {
                cur_rec->a2 = nullptr;
              }
            }
          }
          ++record_ct;
        }
        if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
          goto MetaAnalysis_ret_TSTREAM_FAIL;
        }
        if (unlikely(CleanupTextStream2(cur_fname, &txs, &reterr))) {
          goto MetaAnalysis_ret_1;
        }
        fname_iter = strnul(fname_iter);
        ++fname_iter;
      }
    }
    cur_fname = nullptr;

    STD_SORT(record_ct, MetaRecordCmp, records);

    // Group by ID, then order the groups genomically.
    uintptr_t group_ct = 0;
    for (uintptr_t rec_idx = 0; rec_idx != record_ct; ) {
      uintptr_t rec_idx_end = rec_idx + 1;
      while ((rec_idx_end != record_ct) && (!strcmp(records[rec_idx_end].id, records[rec_idx].id))) {
        ++rec_idx_end;
      }
      ++group_ct;
      rec_idx = rec_idx_end;
    }
    MetaGroup* groups;
    if (unlikely(BIGSTACK_ALLOC_X(MetaGroup, group_ct, &groups))) {
      goto MetaAnalysis_ret_NOMEM;
    }
    {
      uintptr_t group_idx = 0;
      for (uintptr_t rec_idx = 0; rec_idx != record_ct; ) {
        uintptr_t rec_idx_end = rec_idx + 1;
        while ((rec_idx_end != record_ct) && (!strcmp(records[rec_idx_end].id, records[rec_idx].id))) {
          ++rec_idx_end;
        }
        groups[group_idx].chr = records[rec_idx].chr;
        groups[group_idx].bp = records[rec_idx].bp;
        groups[group_idx].first_rec = rec_idx;
        groups[group_idx].rec_ct = rec_idx_end - rec_idx;
        ++group_idx;
        rec_idx = rec_idx_end;
      }
    }
    STD_SORT(group_ct, MetaGroupCmp, groups);

    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 256 + file_ct * 32;
    OutnameZstSet(".meta", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto MetaAnalysis_ret_1;
    }
    *cswritep++ = '#';
    if (!no_map) {
      cswritep = strcpya_k(cswritep, "CHROM\tPOS\t");
    }
    cswritep = strcpya_k(cswritep, "ID");
    if (!no_allele) {
      cswritep = strcpya_k(cswritep, "\tA1\tA2");
    }
    cswritep = strcpya_k(cswritep, "\tN\tP\tP_R\t");
    cswritep = strcpya(cswritep, qt_mode? "BETA\tBETA_R" : "OR\tOR_R");
    cswritep = strcpya_k(cswritep, "\tQ\tI2");
    if (weighted_z) {
      cswritep = strcpya_k(cswritep, "\tWEIGHTED_Z\tP_WZ");
    }
    if (study_col) {
      for (uint32_t file_idx = 0; file_idx != file_ct; ++file_idx) {
        cswritep = strcpya_k(cswritep, "\tF");
        cswritep = u32toa(file_idx, cswritep);
      }
    }
    AppendBinaryEoln(&cswritep);

    double* study_betas = nullptr;
    if (study_col) {
      if (unlikely(bigstack_alloc_d(file_ct, &study_betas))) {
        goto MetaAnalysis_ret_NOMEM;
      }
    }
    uintptr_t written_ct = 0;
    uintptr_t flip_discard_ct = 0;
    for (uintptr_t group_idx = 0; group_idx != group_ct; ++group_idx) {
      const MetaGroup* cur_group = &(groups[group_idx]);
      const MetaRecord* group_recs = &(records[cur_group->first_rec]);
      const uint32_t group_rec_ct = cur_group->rec_ct;
      // The first record fixes the reference orientation; later ones either
      // agree, flip, or are dropped.
      const char* ref_a1 = group_recs[0].a1;
      const char* ref_a2 = group_recs[0].a2;
      double sum_w = 0.0;
      double sum_wbeta = 0.0;
      double sum_w2 = 0.0;
      // METAL's weighted Z: each study's signed Z, weighted by the square root
      // of its effective sample size.
      double sum_wz = 0.0;
      double sum_wz_w2 = 0.0;
      uint32_t used_ct = 0;
      if (study_col) {
        for (uint32_t file_idx = 0; file_idx != file_ct; ++file_idx) {
          study_betas[file_idx] = 0.0 / 0.0;
        }
      }
      for (uint32_t rec_idx = 0; rec_idx != group_rec_ct; ++rec_idx) {
        const MetaRecord* cur_rec = &(group_recs[rec_idx]);
        double cur_beta = cur_rec->beta;
        if (!no_allele) {
          if (strcmp(cur_rec->a1, ref_a1)) {
            // An A1/A2 swap is recoverable when both files name A2.
            if (ref_a2 && cur_rec->a2 && (!strcmp(cur_rec->a1, ref_a2)) && (!strcmp(cur_rec->a2, ref_a1))) {
              cur_beta = -cur_beta;
            } else {
              ++flip_discard_ct;
              continue;
            }
          }
        }
        const double cur_w = 1.0 / (cur_rec->se * cur_rec->se);
        sum_w += cur_w;
        sum_w2 += cur_w * cur_w;
        sum_wbeta += cur_w * cur_beta;
        if (weighted_z) {
          const double abs_z = sqrt(LnPToChisq(cur_rec->ln_pval));
          const double cur_z = (cur_beta < 0.0)? (-abs_z) : abs_z;
          const double cur_wz_w = sqrt(cur_rec->ess);
          sum_wz += cur_wz_w * cur_z;
          sum_wz_w2 += cur_rec->ess;
        }
        if (study_col) {
          study_betas[cur_rec->file_idx] = cur_beta;
        }
        ++used_ct;
      }
      if (!used_ct) {
        continue;
      }
      if ((used_ct == 1) && (!report_all)) {
        continue;
      }
      const double beta_fe = sum_wbeta / sum_w;
      const double se_fe = sqrt(1.0 / sum_w);
      const double z_fe = beta_fe / se_fe;
      const double ln_p_fe = ZscoreToLnP(z_fe);
      // Cochran's Q, and the DerSimonian-Laird random-effects estimate built
      // from it.
      double qq = 0.0;
      for (uint32_t rec_idx = 0; rec_idx != group_rec_ct; ++rec_idx) {
        const MetaRecord* cur_rec = &(group_recs[rec_idx]);
        double cur_beta = cur_rec->beta;
        if (!no_allele) {
          if (strcmp(cur_rec->a1, ref_a1)) {
            if (ref_a2 && cur_rec->a2 && (!strcmp(cur_rec->a1, ref_a2)) && (!strcmp(cur_rec->a2, ref_a1))) {
              cur_beta = -cur_beta;
            } else {
              continue;
            }
          }
        }
        const double cur_w = 1.0 / (cur_rec->se * cur_rec->se);
        const double dxx = cur_beta - beta_fe;
        qq += cur_w * dxx * dxx;
      }
      double beta_re = beta_fe;
      double se_re = se_fe;
      double i2 = 0.0;
      if (used_ct > 1) {
        const double df = u31tod(used_ct - 1);
        if (qq > df) {
          i2 = 100.0 * (qq - df) / qq;
          const double denom = sum_w - (sum_w2 / sum_w);
          if (denom > 0.0) {
            const double tau2 = (qq - df) / denom;
            double sum_wr = 0.0;
            double sum_wrbeta = 0.0;
            for (uint32_t rec_idx = 0; rec_idx != group_rec_ct; ++rec_idx) {
              const MetaRecord* cur_rec = &(group_recs[rec_idx]);
              double cur_beta = cur_rec->beta;
              if (!no_allele) {
                if (strcmp(cur_rec->a1, ref_a1)) {
                  if (ref_a2 && cur_rec->a2 && (!strcmp(cur_rec->a1, ref_a2)) && (!strcmp(cur_rec->a2, ref_a1))) {
                    cur_beta = -cur_beta;
                  } else {
                    continue;
                  }
                }
              }
              const double cur_wr = 1.0 / (cur_rec->se * cur_rec->se + tau2);
              sum_wr += cur_wr;
              sum_wrbeta += cur_wr * cur_beta;
            }
            beta_re = sum_wrbeta / sum_wr;
            se_re = sqrt(1.0 / sum_wr);
          }
        }
      }
      const double ln_p_re = ZscoreToLnP(beta_re / se_re);

      if (!no_map) {
        cswritep = strcpyax(cswritep, cur_group->chr, '\t');
        cswritep = u32toa_x(cur_group->bp, '\t', cswritep);
      }
      cswritep = strcpya(cswritep, group_recs[0].id);
      if (!no_allele) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, ref_a1);
        *cswritep++ = '\t';
        if (ref_a2) {
          cswritep = strcpya(cswritep, ref_a2);
        } else {
          *cswritep++ = '.';
        }
      }
      *cswritep++ = '\t';
      cswritep = u32toa_x(used_ct, '\t', cswritep);
      cswritep = lntoa_g(ln_p_fe, cswritep);
      *cswritep++ = '\t';
      cswritep = lntoa_g(ln_p_re, cswritep);
      *cswritep++ = '\t';
      cswritep = dtoa_g(qt_mode? beta_fe : exp(beta_fe), cswritep);
      *cswritep++ = '\t';
      cswritep = dtoa_g(qt_mode? beta_re : exp(beta_re), cswritep);
      *cswritep++ = '\t';
      if (used_ct > 1) {
        // PLINK 1.x's Q column is the heterogeneity test's p-value, not the
        // statistic itself.
        cswritep = lntoa_g(ChisqToLnP(qq, used_ct - 1), cswritep);
        *cswritep++ = '\t';
        cswritep = dtoa_g(i2, cswritep);
      } else {
        cswritep = strcpya_k(cswritep, "NA\tNA");
      }
      if (weighted_z) {
        const double wz = sum_wz / sqrt(sum_wz_w2);
        *cswritep++ = '\t';
        cswritep = dtoa_g(wz, cswritep);
        *cswritep++ = '\t';
        cswritep = lntoa_g(ZscoreToLnP(wz), cswritep);
      }
      if (study_col) {
        for (uint32_t file_idx = 0; file_idx != file_ct; ++file_idx) {
          *cswritep++ = '\t';
          const double cur_beta = study_betas[file_idx];
          if (cur_beta != cur_beta) {
            cswritep = strcpya_k(cswritep, "NA");
          } else {
            cswritep = dtoa_g(qt_mode? cur_beta : exp(cur_beta), cswritep);
          }
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto MetaAnalysis_ret_WRITE_FAIL;
      }
      ++written_ct;
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto MetaAnalysis_ret_WRITE_FAIL;
    }
    logprintfww("--meta-analysis: %" PRIuPTR " variant%s written to %s .\n", written_ct, (written_ct == 1)? "" : "s", outname);
    if (flip_discard_ct) {
      logerrprintf("Warning: %" PRIuPTR " row%s discarded due to an unresolvable allele mismatch.\n", flip_discard_ct, (flip_discard_ct == 1)? "" : "s");
    }
  }
  while (0) {
  MetaAnalysis_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MetaAnalysis_ret_TSTREAM_FAIL:
    TextStreamErrPrint(cur_fname, &txs);
    break;
  MetaAnalysis_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  MetaAnalysis_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, cur_fname);
  MetaAnalysis_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  MetaAnalysis_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  MetaAnalysis_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 MetaAnalysis_ret_1:
  CswriteCloseCond(&css, cswritep);
  if (cur_fname) {
    CleanupTextStream2(cur_fname, &txs, &reterr);
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr AdjustFile(const AdjustFileInfo* afip, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* in_fname = afip->fname;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream adjust_txs;
  PreinitTextStream(&adjust_txs);
  {
    // Two-pass load.
    // 1. Parse header line, count # of variants.
    // intermission. Allocate top-level arrays.
    // 2. Rewind and fill arrays.
    // (some overlap with LoadPvar(), though that's one-pass.)
    reterr = SizeAndInitTextStream(in_fname, bigstack_left() / 4, max_thread_ct, &adjust_txs);
    if (unlikely(reterr)) {
      goto AdjustFile_ret_TSTREAM_FAIL;
    }

    const char* header_start;
    do {
      ++line_idx;
      header_start = TextGet(&adjust_txs);
      if (unlikely(!header_start)) {
        reterr = TextStreamRawErrcode(&adjust_txs);
        if (reterr == kPglRetEof) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", in_fname);
          goto AdjustFile_ret_MALFORMED_INPUT_WW;
        }
        goto AdjustFile_ret_TSTREAM_FAIL;
      }
    } while (strequal_k_unsafe(header_start, "##"));
    if (*header_start == '#') {
      ++header_start;
    }

    const AdjustFlags flags = afip->base.flags;
    // [0] = CHROM
    // [1] = POS
    // [2] = ID (required)
    // [3] = REF
    // [4] = ALT
    // [5] = PROVISIONAL_REF?
    // [6] = A1
    // [7] = TEST (always scan)
    // [8] = P (required)
    const char* col_search_order[9];
    const uint32_t need_chr = (flags & kfAdjustColChrom);
    const uint32_t need_pos = (flags & kfAdjustColPos);
    const uint32_t need_ref = (flags & kfAdjustColRef);
    const uint32_t need_alt = (flags & (kfAdjustColAlt1 | kfAdjustColAlt));
    const uint32_t need_provref = (flags & (kfAdjustColMaybeprovref | kfAdjustColProvref));
    uint32_t check_a1 = (flags & kfAdjustColA1);
    const uint32_t alt_comma_truncate = (need_alt == kfAdjustColAlt1);
    if (unlikely(need_alt == (kfAdjustColAlt1 | kfAdjustColAlt))) {
      // Could theoretically support this later (allocate allele_idx_offsets
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
    col_search_order[4] = need_alt? (afip->alt_field? afip->alt_field : "ALT\0ALT1\0") : "";
    col_search_order[5] = need_provref? (afip->provref_field? afip->provref_field : "PROVISIONAL_REF?\0") : "";
    col_search_order[6] = check_a1? (afip->a1_field? afip->a1_field : "A1\0") : "";
    col_search_order[7] = afip->test_field? afip->test_field : "TEST\0";
    const uint32_t input_log10 = (flags & kfAdjustInputLog10);
    col_search_order[8] = afip->p_field? afip->p_field : (input_log10? "LOG10_P\0NEG_LOG10_P\0LOG10_UNADJ\0NEG_LOG10_UNADJ\0P\0UNADJ\0" : "P\0UNADJ\0");

    uint32_t col_skips[9];
    uint32_t col_types[9];
    uint32_t relevant_col_ct;
    uint32_t found_type_bitset;
    reterr = SearchHeaderLine(header_start, col_search_order, "--adjust-file", 9, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
    if (unlikely(reterr)) {
      goto AdjustFile_ret_1;
    }
    if (unlikely((found_type_bitset & 0x104) != 0x104)) {
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
      for (uint32_t relevant_col_idx = 0; ; ++relevant_col_idx) {
        test_col_idx += col_skips[relevant_col_idx];
        if (col_types[relevant_col_idx] == 7) {
          break;
        }
      }
    } else if (unlikely(found_type_bitset & 0x80)) {
      snprintf(g_logbuf, kLogbufSize, "Error: TEST column present in %s, but no test= parameter was provided to --adjust-file.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(need_chr && (!(found_type_bitset & 0x1)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No chromosome column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(need_pos && (!(found_type_bitset & 0x2)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No bp coordinate column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(need_ref && (!(found_type_bitset & 0x8)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No REF column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(need_alt && (!(found_type_bitset & 0x10)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No ALT column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(need_provref && (!(found_type_bitset & 0x20)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: No PROVISIONAL_REF? column in %s.\n", in_fname);
      goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
    }
    if (check_a1 && (!(found_type_bitset & 0x40))) {
      snprintf(g_logbuf, kLogbufSize, "Warning: No A1 column in %s. Omitting from output.\n", in_fname);
      check_a1 = 0;
    }

    uintptr_t entry_ct = 0;
    while (1) {
      ++line_idx;
      const char* line_start = TextGet(&adjust_txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&adjust_txs, &reterr))) {
          break;
        }
        goto AdjustFile_ret_TSTREAM_FAIL;
      }
      if (test_name) {
        // Don't count different-test entries.
        const char* test_name_start = NextTokenMult0(line_start, test_col_idx);
        if (unlikely(!test_name_start)) {
          goto AdjustFile_ret_MISSING_TOKENS;
        }
        const uint32_t cur_test_slen = strlen_se(test_name_start);
        if ((cur_test_slen != test_name_slen) || (!memequal(test_name_start, test_name, test_name_slen))) {
          continue;
        }
      }
      ++entry_ct;
    }
#ifdef __LP64__
    // probably want to permit this soon
    if (entry_ct > 0xffffffffU) {
      logerrputs("Error: Too many entries for --adjust-file.\n");
      reterr = kPglRetNotYetSupported;
      goto AdjustFile_ret_1;
    }
#endif

    reterr = TextRewind(&adjust_txs);
    if (unlikely(reterr)) {
      goto AdjustFile_ret_TSTREAM_FAIL;
    }
    const uintptr_t line_ct = line_idx - 1;
    line_idx = 0;
    do {
      ++line_idx;
      reterr = TextNextLineLstripK(&adjust_txs, &header_start);
      if (unlikely(reterr)) {
        goto AdjustFile_ret_TSTREAM_REWIND_FAIL;
      }
    } while (strequal_k_unsafe(header_start, "##"));

    const uintptr_t entry_ctl = BitCtToWordCt(entry_ct);
    const uintptr_t entry_ctl2 = NypCtToWordCt(entry_ct);
    uintptr_t* variant_include_dummy;
    uintptr_t* allele_include_dummy;
    if (unlikely(bigstack_alloc_w(entry_ctl, &variant_include_dummy) ||
                 bigstack_alloc_w(entry_ctl2, &allele_include_dummy))) {
      goto AdjustFile_ret_NOMEM;
    }
    SetAllBits(entry_ct, variant_include_dummy);
    for (uintptr_t ulii = 0; ulii != entry_ct / kBitsPerWordD2; ++ulii) {
      allele_include_dummy[ulii] = kMaskAAAA;
    }
    const uint32_t remainder = entry_ct % kBitsPerWordD2;
    if (remainder) {
      allele_include_dummy[entry_ct / kBitsPerWordD2] = kMaskAAAA >> (2 * (kBitsPerWordD2 - remainder));
    }
    char** chr_ids;
    if (need_chr) {
      if (unlikely(bigstack_alloc_cp(entry_ct, &chr_ids))) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      chr_ids = nullptr;
    }
    uint32_t* variant_bps;
    if (need_pos) {
      if (unlikely(bigstack_alloc_u32(entry_ct, &variant_bps))) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      variant_bps = nullptr;
    }
    char** variant_ids;
    double* ln_pvals;
    if (unlikely(bigstack_alloc_cp(entry_ct, &variant_ids) ||
                 bigstack_alloc_d(entry_ct, &ln_pvals))) {
      goto AdjustFile_ret_NOMEM;
    }
    char** allele_storage;
    if (need_ref || need_alt) {
      if (unlikely(bigstack_alloc_cp(entry_ct * 2, &allele_storage))) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      allele_storage = nullptr;
    }
    uintptr_t* nonref_flags = nullptr;
    if (need_provref) {
      if (unlikely(bigstack_calloc_w(entry_ctl, &nonref_flags))) {
        goto AdjustFile_ret_NOMEM;
      }
    }
    char** a1_storage;
    if (check_a1) {
      if (unlikely(bigstack_alloc_cp(entry_ct, &a1_storage))) {
        goto AdjustFile_ret_NOMEM;
      }
    } else {
      a1_storage = nullptr;
    }
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = BigstackEndRoundedDown();
    uint32_t max_allele_slen = 1;
    uintptr_t variant_idx = 0;
    while (line_idx < line_ct) {
      ++line_idx;
      const char* line_start = nullptr;  // gcc 14 warning
      reterr = TextNextLineLstripK(&adjust_txs, &line_start);
      if (unlikely(reterr)) {
        goto AdjustFile_ret_TSTREAM_REWIND_FAIL;
      }
      const char* token_ptrs[9];
      uint32_t token_slens[9];
      if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
        goto AdjustFile_ret_MISSING_TOKENS;
      }
      if (test_name) {
        if ((token_slens[7] != test_name_slen) || (!memequal(token_ptrs[7], test_name, test_name_slen))) {
          continue;
        }
      }
      if (chr_ids) {
        const uint32_t cur_slen = token_slens[0];
        if (StoreStringAtBase(tmp_alloc_end, token_ptrs[0], cur_slen, &tmp_alloc_base, &(chr_ids[variant_idx]))) {
          goto AdjustFile_ret_NOMEM;
        }
      }
      if (variant_bps) {
        if (unlikely(ScanUintDefcap(token_ptrs[1], &(variant_bps[variant_idx])))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, in_fname);
          goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
        }
      }
      const uint32_t id_slen = token_slens[2];
      if (StoreStringAtBase(tmp_alloc_end, token_ptrs[2], id_slen, &tmp_alloc_base, &(variant_ids[variant_idx]))) {
        goto AdjustFile_ret_NOMEM;
      }
      if (need_ref) {
        const uint32_t cur_slen = token_slens[3];
        if (StoreStringAtBase(tmp_alloc_end, token_ptrs[3], cur_slen, &tmp_alloc_base, &(allele_storage[2 * variant_idx]))) {
          goto AdjustFile_ret_NOMEM;
        }
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
        if (StoreStringAtBase(tmp_alloc_end, alt_str, cur_slen, &tmp_alloc_base, &(allele_storage[2 * variant_idx + 1]))) {
          goto AdjustFile_ret_NOMEM;
        }
      }
      if (nonref_flags) {
        const char provref_char = token_ptrs[5][0];
        const uint32_t cur_slen = token_slens[5];
        if ((provref_char == 'Y') && (cur_slen == 1)) {
          SetBit(variant_idx, nonref_flags);
        } else if (unlikely((provref_char != 'N') || (cur_slen != 1))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid provref entry on line %" PRIuPTR " of %s.\n", line_idx, in_fname);
          goto AdjustFile_ret_INCONSISTENT_INPUT_WW;
        }
      }
      if (check_a1) {
        const uint32_t cur_slen = token_slens[6];
        if (StoreStringAtBase(tmp_alloc_end, token_ptrs[6], cur_slen, &tmp_alloc_base, &(a1_storage[variant_idx]))) {
          goto AdjustFile_ret_NOMEM;
        }
      }
      const char* pval_str = token_ptrs[8];
      double ln_pval;
      if (!input_log10) {
        if (!ScantokLn(pval_str, &ln_pval)) {
          uint32_t cur_slen;
        AdjustFile_alphabetic_pval:
          cur_slen = token_slens[8];
          if (IsNanStr(pval_str, cur_slen)) {
            ln_pval = kLnPvalError;
          } else if (likely(strequal_k(pval_str, "INF", cur_slen) ||
                            (input_log10 && strequal_k(pval_str, "inf", cur_slen)))) {
            // From plink 1.x, could be anything smaller than log(5e-324).
            // Just fill with log(2.23e-308) for now.
            ln_pval = kLnNormalMin;
          } else {
            goto AdjustFile_ret_INVALID_PVAL;
          }
        }
      } else {
        double neglog10_pval;
        if (!ScantokDouble(pval_str, &neglog10_pval)) {
          goto AdjustFile_alphabetic_pval;
        }
        ln_pval = neglog10_pval * (-kLn10);
        if (unlikely(ln_pval > 0.0)) {
          goto AdjustFile_ret_INVALID_PVAL;
        }
      }
      ln_pvals[variant_idx] = ln_pval;
      ++variant_idx;
    }
    BigstackEndReset(bigstack_end_mark);
    BigstackBaseSet(tmp_alloc_base);
    reterr = Multcomp(variant_include_dummy, nullptr, TO_CONSTCPCONSTP(chr_ids), variant_bps, TO_CONSTCPCONSTP(variant_ids), allele_include_dummy, nullptr, TO_CONSTCPCONSTP(allele_storage), nonref_flags, TO_CONSTCPCONSTP(a1_storage), &(afip->base), ln_pvals, nullptr, entry_ct, entry_ct, max_allele_slen, kfPgenGlobal0, ln_pfilter, output_min_ln, 0, max_thread_ct, outname, outname_end);
    if (unlikely(reterr)) {
      goto AdjustFile_ret_1;
    }
  }
  while (0) {
  AdjustFile_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  AdjustFile_ret_TSTREAM_FAIL:
    TextStreamErrPrint(in_fname, &adjust_txs);
    break;
  AdjustFile_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(in_fname, &adjust_txs, &reterr);
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
  CleanupTextStream2(in_fname, &adjust_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
