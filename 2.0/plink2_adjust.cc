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
#include "include/plink2_htable.h"
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
    reterr = SearchHeaderLine(header_start, col_search_order, "adjust-file", 9, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
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

void InitAcat(AcatInfo* acat_info_ptr) {
  acat_info_ptr->flags = kfAcat0;
  acat_info_ptr->fname = nullptr;
  acat_info_ptr->test_name = nullptr;
  acat_info_ptr->id_field = nullptr;
  acat_info_ptr->test_field = nullptr;
  acat_info_ptr->p_field = nullptr;
  acat_info_ptr->freq_field = nullptr;
  acat_info_ptr->beta_a1 = 1.0;
  acat_info_ptr->beta_a2 = 25.0;
}

void CleanupAcat(AcatInfo* acat_info_ptr) {
  free_cond(acat_info_ptr->fname);
  free_cond(acat_info_ptr->test_name);
  free_cond(acat_info_ptr->id_field);
  free_cond(acat_info_ptr->test_field);
  free_cond(acat_info_ptr->p_field);
  free_cond(acat_info_ptr->freq_field);
}

// Beta(maf; a1, a2) density, the standard rare-variant weight.  Squared and
// multiplied by maf * (1 - maf) it becomes the ACAT-V variant weight, which is
// what makes a rarer variant count for more.
static double BetaDensity(double xx, double a1, double a2) {
  if ((xx <= 0.0) || (xx >= 1.0)) {
    return 0.0;
  }
  const double ln_beta = lgamma(a1) + lgamma(a2) - lgamma(a1 + a2);
  return exp((a1 - 1.0) * log(xx) + (a2 - 1.0) * log1p(-xx) - ln_beta);
}

PglErr AcatSets(const AcatInfo* acip, const char* set_fname, double output_min_ln, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* in_fname = acip->fname;
  uintptr_t line_idx = 0;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  TextStream txs;
  PreinitCstream(&css);
  PreinitTextStream(&txs);
  {
    if (unlikely(!set_fname)) {
      logerrputs("Error: --acat-file requires --set-list.\n");
      goto AcatSets_ret_INCONSISTENT_INPUT;
    }
    // Pass 1 counts the rows we will keep and measures the longest variant ID;
    // pass 2 fills the arrays.  Same shape as AdjustFile().
    reterr = SizeAndInitTextStream(in_fname, bigstack_left() / 4, max_thread_ct, &txs);
    if (unlikely(reterr)) {
      goto AcatSets_ret_TSTREAM_FAIL;
    }
    const char* header_start;
    do {
      ++line_idx;
      header_start = TextGet(&txs);
      if (unlikely(!header_start)) {
        reterr = TextStreamRawErrcode(&txs);
        if (reterr == kPglRetEof) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", in_fname);
          goto AcatSets_ret_MALFORMED_INPUT_WW;
        }
        goto AcatSets_ret_TSTREAM_FAIL;
      }
    } while (strequal_k_unsafe(header_start, "##"));
    if (*header_start == '#') {
      ++header_start;
    }
    const uint32_t input_log10 = (acip->flags / kfAcatInputLog10) & 1;
    // [0] = ID, [1] = TEST, [2] = P, [3] = A1_FREQ
    const char* col_search_order[4];
    col_search_order[0] = acip->id_field? acip->id_field : "ID\0SNP\0";
    col_search_order[1] = acip->test_field? acip->test_field : "TEST\0";
    col_search_order[2] = acip->p_field? acip->p_field : (input_log10? "LOG10_P\0NEG_LOG10_P\0P\0" : "P\0UNADJ\0");
    col_search_order[3] = acip->freq_field? acip->freq_field : "A1_FREQ\0MAF\0FREQ\0";
    uint32_t col_skips[4];
    uint32_t col_types[4];
    uint32_t relevant_col_ct;
    uint32_t found_type_bitset;
    reterr = SearchHeaderLine(header_start, col_search_order, "--acat-file", 4, &relevant_col_ct, &found_type_bitset, col_skips, col_types);
    if (unlikely(reterr)) {
      goto AcatSets_ret_1;
    }
    if (unlikely((found_type_bitset & 5) != 5)) {
      logerrputs("Error: --acat-file requires ID and P columns.\n");
      goto AcatSets_ret_INCONSISTENT_INPUT;
    }
    const uint32_t have_freq = (found_type_bitset >> 3) & 1;
    const char* test_name = acip->test_name;
    const uint32_t test_name_slen = test_name? strlen(test_name) : 0;
    if (unlikely(test_name && (!((found_type_bitset >> 1) & 1)))) {
      logerrputs("Error: --acat-file test= was specified, but the file has no TEST column.\n");
      goto AcatSets_ret_INCONSISTENT_INPUT;
    }

    uintptr_t variant_ct = 0;
    uintptr_t max_id_blen = 2;
    while (1) {
      ++line_idx;
      const char* line_start = TextGet(&txs);
      if (!line_start) {
        break;
      }
      const char* token_ptrs[4];
      uint32_t token_slens[4];
      if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
        goto AcatSets_ret_MISSING_TOKENS;
      }
      if (test_name) {
        if ((token_slens[1] != test_name_slen) || (!memequal(token_ptrs[1], test_name, test_name_slen))) {
          continue;
        }
      }
      if (token_slens[0] >= max_id_blen) {
        max_id_blen = token_slens[0] + 1;
      }
      ++variant_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto AcatSets_ret_TSTREAM_FAIL;
    }
    if (unlikely(!variant_ct)) {
      logerrputs("Error: --acat-file: no association results to combine.\n");
      goto AcatSets_ret_INCONSISTENT_INPUT;
    }
#ifdef __LP64__
    if (unlikely(variant_ct > 0xffffffffU)) {
      logerrputs("Error: Too many variants for --acat-file.\n");
      goto AcatSets_ret_MALFORMED_INPUT;
    }
#endif
    char* variant_ids;
    double* ln_pvals;
    double* weights;
    if (unlikely(bigstack_alloc_c(variant_ct * max_id_blen, &variant_ids) ||
                 bigstack_alloc_d(variant_ct, &ln_pvals) ||
                 bigstack_alloc_d(variant_ct, &weights))) {
      goto AcatSets_ret_NOMEM;
    }

    reterr = TextRewind(&txs);
    if (unlikely(reterr)) {
      goto AcatSets_ret_TSTREAM_FAIL;
    }
    line_idx = 0;
    do {
      ++line_idx;
      reterr = TextNextLineLstripK(&txs, &header_start);
      if (unlikely(reterr)) {
        goto AcatSets_ret_TSTREAM_REWIND_FAIL;
      }
    } while (strequal_k_unsafe(header_start, "##"));

    const double beta_a1 = acip->beta_a1;
    const double beta_a2 = acip->beta_a2;
    const double ln_ten = kLn10;
    uintptr_t variant_idx = 0;
    while (variant_idx < variant_ct) {
      ++line_idx;
      const char* line_start = TextGet(&txs);
      if (unlikely(!line_start)) {
        break;
      }
      const char* token_ptrs[4];
      uint32_t token_slens[4];
      if (unlikely(!TokenLexK0(line_start, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens))) {
        goto AcatSets_ret_MISSING_TOKENS;
      }
      if (test_name) {
        if ((token_slens[1] != test_name_slen) || (!memequal(token_ptrs[1], test_name, test_name_slen))) {
          continue;
        }
      }
      memcpyx(&(variant_ids[variant_idx * max_id_blen]), token_ptrs[0], token_slens[0], '\0');

      const char* pval_str = token_ptrs[2];
      double cur_ln_pval;
      if (IsNanStr(pval_str, token_slens[2])) {
        // A variant with no p-value contributes nothing rather than poisoning
        // its set.
        cur_ln_pval = 1.0;  // sentinel; positive is impossible for a log p
      } else {
        double dxx;
        if (unlikely(!ScantokDouble(pval_str, &dxx))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid p-value on line %" PRIuPTR " of %s.\n", line_idx, in_fname);
          goto AcatSets_ret_MALFORMED_INPUT_WW;
        }
        if (input_log10) {
          cur_ln_pval = (-dxx) * ln_ten;
        } else if (dxx > 0.0) {
          cur_ln_pval = log(dxx);
        } else {
          // An exact zero in the input file is a p-value that underflowed on
          // the way out of whatever wrote it.  Treat it as the smallest
          // representable rather than -inf, which would make the combination
          // degenerate.
          cur_ln_pval = -745.0;
        }
        if (cur_ln_pval > 0.0) {
          cur_ln_pval = 0.0;
        }
      }
      ln_pvals[variant_idx] = cur_ln_pval;

      double cur_weight = 1.0;
      if (have_freq) {
        double freq;
        if (ScantokDouble(token_ptrs[3], &freq) && (freq > 0.0) && (freq < 1.0)) {
          const double maf = (freq > 0.5)? (1.0 - freq) : freq;
          const double beta_wt = BetaDensity(maf, beta_a1, beta_a2);
          cur_weight = beta_wt * beta_wt * maf * (1.0 - maf);
          if (!(cur_weight > 0.0)) {
            cur_weight = 0.0;
          }
        } else {
          cur_weight = 0.0;
        }
      }
      weights[variant_idx] = cur_weight;
      ++variant_idx;
    }
    // No TextStreamErrcode2() here on purpose: this loop stops on the variant
    // count rather than on end-of-file, so the stream is normally still mid-
    // file, and TextStreamErrcode2() reports anything that is not EOF as an
    // error.
    const uint32_t final_variant_ct = variant_idx;
    CleanupTextStream2(in_fname, &txs, &reterr);
    if (unlikely(reterr)) {
      goto AcatSets_ret_1;
    }

    // Variant-ID lookup.  Duplicate IDs are allowed here: a --glm output can
    // legitimately carry one row per ALT allele, and every such row is a
    // p-value for the set.  PopulateStrboxHtable() keeps the first, so
    // duplicates are resolved by a linear rescan below.
    const uint32_t id_htable_size = GetHtableFastSize(final_variant_ct);
    uint32_t* id_htable;
    if (unlikely(bigstack_alloc_u32(id_htable_size, &id_htable))) {
      goto AcatSets_ret_NOMEM;
    }
    PopulateStrboxHtable(variant_ids, final_variant_ct, max_id_blen, id_htable_size, id_htable);

    // A set can name at most every variant in the results file, so one
    // allocation up front removes any need to grow these later.  That matters:
    // the parsed set line points into the text stream's own buffer, so
    // reallocating underneath it would leave those pointers dangling.
    double* set_ln_pvals;
    double* set_weights;
    if (unlikely(bigstack_alloc_d(final_variant_ct, &set_ln_pvals) ||
                 bigstack_alloc_d(final_variant_ct, &set_weights))) {
      goto AcatSets_ret_NOMEM;
    }

    // Second stream: the set definitions.
    reterr = SizeAndInitTextStream(set_fname, bigstack_left() / 4, MAXV(max_thread_ct, 1), &txs);
    if (unlikely(reterr)) {
      goto AcatSets_ret_TSTREAM_SET_FAIL;
    }
    OutnameZstSet(".acat", acip->flags & kfAcatZs, outname_end);
    reterr = InitCstreamAlloc(outname, 0, acip->flags & kfAcatZs, MAXV(max_thread_ct, 1), kCompressStreamBlock + kMaxMediumLine, &css, &cswritep);
    if (unlikely(reterr)) {
      goto AcatSets_ret_1;
    }
    cswritep = strcpya_k(cswritep, "#SET\tCHROM\tPOS\tNVAR\tNVAR_TESTED\tP" EOLN_STR);

    uint32_t set_ct = 0;
    uint64_t skipped_set_ct = 0;
    line_idx = 0;
    while (1) {
      ++line_idx;
      const char* line_start = TextGet(&txs);
      if (!line_start) {
        break;
      }
      if ((*line_start == '#') || (*line_start == '\0')) {
        continue;
      }
      // <set name> <chromosome> <position> <comma-separated variant IDs>
      const char* set_name = line_start;
      const char* set_name_end = CurTokenEnd(set_name);
      const char* chr_str = FirstNonTspace(set_name_end);
      if (unlikely(IsEolnKns(*chr_str))) {
        goto AcatSets_ret_SET_MISSING_TOKENS;
      }
      const char* chr_end = CurTokenEnd(chr_str);
      const char* pos_str = FirstNonTspace(chr_end);
      if (unlikely(IsEolnKns(*pos_str))) {
        goto AcatSets_ret_SET_MISSING_TOKENS;
      }
      const char* pos_end = CurTokenEnd(pos_str);
      const char* id_list = FirstNonTspace(pos_end);
      if (unlikely(IsEolnKns(*id_list))) {
        goto AcatSets_ret_SET_MISSING_TOKENS;
      }
      const char* id_list_end = CurTokenEnd(id_list);

      // Count members, growing the per-set buffers as needed.
      uint32_t member_ct = 1;
      for (const char* scan = id_list; scan != id_list_end; ++scan) {
        if (*scan == ',') {
          ++member_ct;
        }
      }
      uint32_t tested_ct = 0;
      const char* id_iter = id_list;
      while (1) {
        const char* id_end = id_iter;
        while ((id_end != id_list_end) && (*id_end != ',')) {
          ++id_end;
        }
        const uint32_t id_slen = id_end - id_iter;
        if (id_slen) {
          // Nnt: the ID runs up to a comma or the end of the line, so it is not
          // null-terminated the way StrboxHtableFind() requires.
          uint32_t variant_uidx = StrboxHtableFindNnt(id_iter, variant_ids, id_htable, max_id_blen, id_slen, id_htable_size);
          while (variant_uidx != UINT32_MAX) {
            const double cur_ln_p = ln_pvals[variant_uidx];
            if (cur_ln_p <= 0.0) {
              set_ln_pvals[tested_ct] = cur_ln_p;
              set_weights[tested_ct] = weights[variant_uidx];
              ++tested_ct;
            }
            variant_uidx = UINT32_MAX;
          }
        }
        if (id_end == id_list_end) {
          break;
        }
        id_iter = &(id_end[1]);
      }

      if (!tested_ct) {
        ++skipped_set_ct;
        continue;
      }
      cswritep = memcpyax(cswritep, set_name, set_name_end - set_name, '\t');
      cswritep = memcpyax(cswritep, chr_str, chr_end - chr_str, '\t');
      cswritep = memcpyax(cswritep, pos_str, pos_end - pos_str, '\t');
      cswritep = u32toa_x(member_ct, '\t', cswritep);
      cswritep = u32toa_x(tested_ct, '\t', cswritep);
      double set_ln_p = AcatCombineLnP(set_ln_pvals, set_weights, tested_ct);
      if (set_ln_p > output_min_ln) {
        cswritep = lntoa_g(set_ln_p, cswritep);
      } else {
        cswritep = lntoa_g(output_min_ln, cswritep);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto AcatSets_ret_WRITE_FAIL;
      }
      ++set_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto AcatSets_ret_TSTREAM_SET_FAIL;
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto AcatSets_ret_WRITE_FAIL;
    }
    logprintfww("--acat-file: %u set%s written to %s .\n", set_ct, (set_ct == 1)? "" : "s", outname);
    if (skipped_set_ct) {
      logerrprintfww("Warning: %" PRIu64 " set%s skipped, with no member p-value in %s.\n", skipped_set_ct, (skipped_set_ct == 1)? "" : "s", in_fname);
    }
  }
  while (0) {
  AcatSets_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  AcatSets_ret_TSTREAM_FAIL:
    TextStreamErrPrint(in_fname, &txs);
    break;
  AcatSets_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(in_fname, &txs, &reterr);
    break;
  AcatSets_ret_TSTREAM_SET_FAIL:
    TextStreamErrPrint(set_fname, &txs);
    break;
  AcatSets_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, in_fname);
  AcatSets_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  AcatSets_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  AcatSets_ret_SET_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, set_fname);
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  AcatSets_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  AcatSets_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 AcatSets_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2(in_fname, &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
