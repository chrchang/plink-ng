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

#include "plink_assoc.h"
#include "plink_cluster.h"
#include "plink_ld.h"
#include "plink_matrix.h"
#include "plink_perm.h"
#include "plink_stats.h"

void aperm_init(Aperm_info* apip) {
  apip->min = 6;
  apip->max = 1000000;
  apip->alpha = 0;
  apip->beta = 0.0001;
  apip->init_interval = 1;
  apip->interval_slope = 0.001;
}

void single_marker_cc_freqs(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
  // Counts the number of A2 alleles and missing calls for both cases and
  // controls, for an autosomal marker.  (The caller is expected to calculate
  // the A1 allele count.)
  // See single_marker_freqs_and_hwe() for discussion.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & (~B)
  // missing: popcount(C)
  // A2: [popcount(A) + popcount(B)] - popcount(C)
  uint32_t tot_ctrl_ab = 0;
  uint32_t tot_ctrl_c = 0;
  uint32_t tot_case_ab = 0;
  uint32_t tot_case_c = 0;
  uintptr_t* lptr_end = &(lptr[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
#ifdef __LP64__
  uintptr_t cur_decr = 60;
  uintptr_t* lptr_6x_end;
  sample_ctl2 -= sample_ctl2 % 6;
  while (sample_ctl2 >= 60) {
  single_marker_cc_freqs_loop:
    lptr_6x_end = &(lptr[cur_decr]);
    count_2freq_dbl_960b((__m128i*)lptr, (__m128i*)lptr_6x_end, (__m128i*)ctrl_include2, (__m128i*)case_include2, &tot_ctrl_ab, &tot_ctrl_c, &tot_case_ab, &tot_case_c);
    lptr = lptr_6x_end;
    ctrl_include2 = &(ctrl_include2[cur_decr]);
    case_include2 = &(case_include2[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto single_marker_cc_freqs_loop;
  }
#else
  uintptr_t* lptr_six_end = &(lptr[sample_ctl2 - (sample_ctl2 % 6)]);
  while (lptr < lptr_six_end) {
    count_2freq_dbl_24b(lptr, ctrl_include2, case_include2, &tot_ctrl_ab, &tot_ctrl_c, &tot_case_ab, &tot_case_c);
    lptr = &(lptr[6]);
    ctrl_include2 = &(ctrl_include2[6]);
    case_include2 = &(case_include2[6]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *ctrl_include2++;
    loader3 = loader >> 1;
    loader4 = loader2 & loader;
    tot_ctrl_ab += popcount2_long(loader4 + (loader3 & loader2));
    tot_ctrl_c += popcount2_long(loader4 & (~loader3));
    loader2 = *case_include2++;
    loader4 = loader2 & loader;
    tot_case_ab += popcount2_long(loader4 + (loader3 & loader2));
    tot_case_c += popcount2_long(loader4 & (~loader3));
  }
  *ctrl_missingp = tot_ctrl_c;
  *ctrl_setp = tot_ctrl_ab - tot_ctrl_c;
  *case_missingp = tot_case_c;
  *case_setp = tot_case_ab - tot_case_c;
}

void haploid_single_marker_cc_freqs(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
  // Counts the number of A1 and A2 alleles for both cases and controls, for a
  // haploid marker.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := B ^ A
  // missing: popcount(C)
  // A2: popcount(A & B)
  uint32_t tot_ctrl_ab = 0;
  uint32_t tot_ctrl_c = 0;
  uint32_t tot_case_ab = 0;
  uint32_t tot_case_c = 0;
  uintptr_t* lptr_end = &(lptr[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = loader >> 1;
    loader3 = loader2 ^ loader;
    loader &= loader2;
    loader2 = *ctrl_include2++;
    tot_ctrl_ab += popcount2_long(loader & loader2);
    tot_ctrl_c += popcount2_long(loader3 & loader2);
    loader2 = *case_include2++;
    tot_case_ab += popcount2_long(loader & loader2);
    tot_case_c += popcount2_long(loader3 & loader2);
  }
  *ctrl_setp = tot_ctrl_ab;
  *ctrl_missingp = tot_ctrl_c;
  *case_setp = tot_case_ab;
  *case_missingp = tot_case_c;
}

void single_marker_cc_3freqs(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_hom2p, uint32_t* ctrl_hetp, uint32_t* ctrl_missingp, uint32_t* case_hom2p, uint32_t* case_hetp, uint32_t* case_missingp) {
  // Counts the number of heterozygotes, A2 homozygotes, and missing calls for
  // both cases and controls.  Assumes marker is diploid.  The caller is
  // expected to calculate the A1 allele count.
  // See single_marker_freqs_and_hwe() for discussion.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  // hom2: popcount(C)
  // het: popcount(B) - popcount(C)
  // missing: popcount(A) - popcount(C)
  uint32_t tot_ctrl_a = 0;
  uint32_t tot_ctrl_b = 0;
  uint32_t tot_ctrl_c = 0;
  uint32_t tot_case_a = 0;
  uint32_t tot_case_b = 0;
  uint32_t tot_case_c = 0;
  uintptr_t* lptr_end = &(lptr[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  sample_ctl2 -= sample_ctl2 % 12;
  while (sample_ctl2 >= 120) {
  single_marker_cc_3freqs_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_1920b((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)ctrl_include2, &tot_ctrl_a, &tot_ctrl_b, &tot_ctrl_c);
    count_3freq_1920b((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)case_include2, &tot_case_a, &tot_case_b, &tot_case_c);
    lptr = lptr_12x_end;
    ctrl_include2 = &(ctrl_include2[cur_decr]);
    case_include2 = &(case_include2[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto single_marker_cc_3freqs_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[sample_ctl2 - (sample_ctl2 % 12)]);
  while (lptr < lptr_twelve_end) {
    count_3freq_48b(lptr, ctrl_include2, &tot_ctrl_a, &tot_ctrl_b, &tot_ctrl_c);
    count_3freq_48b(lptr, case_include2, &tot_case_a, &tot_case_b, &tot_case_c);
    lptr = &(lptr[12]);
    ctrl_include2 = &(ctrl_include2[12]);
    case_include2 = &(case_include2[12]);
  }
#endif
  while (lptr < lptr_end) {
    //   A := genotype & 0x5555...
    //   B := (genotype >> 1) & 0x5555...
    //   C := A & B
    //   popcount(C) = homozyg major ct
    //   popcount(B) = het ct + homozyg major ct
    //   popcount(A) = missing_ct + homozyg major ct
    loader = *lptr++;
    loader2 = *ctrl_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_ctrl_a += popcount2_long(loader2);
    tot_ctrl_b += popcount2_long(loader3);
    tot_ctrl_c += popcount2_long(loader2 & loader3);
    loader2 = *case_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_case_a += popcount2_long(loader2);
    tot_case_b += popcount2_long(loader3);
    tot_case_c += popcount2_long(loader2 & loader3);
  }
  *ctrl_hom2p = tot_ctrl_c;
  *ctrl_hetp = tot_ctrl_b - tot_ctrl_c;
  *ctrl_missingp = tot_ctrl_a - tot_ctrl_c;
  *case_hom2p = tot_case_c;
  *case_hetp = tot_case_b - tot_case_c;
  *case_missingp = tot_case_a - tot_case_c;
}

static inline void adjust_print(double pval, double output_min_p, const char* output_min_p_str, uint32_t output_min_p_strlen, char** bufpp) {
  if (pval < 0) {
    *bufpp = memcpya(*bufpp, "        NA ", 11);
  } else if (pval <= output_min_p) {
    *bufpp = memcpya(*bufpp, output_min_p_str, output_min_p_strlen);
  } else {
    *bufpp = dtoa_g_wxp4x(pval, 10, ' ', *bufpp);
  }
}

static inline void adjust_print_log10(double pval, double output_min_p, const char* output_min_logp_str, uint32_t output_min_logp_strlen, char** bufpp) {
  if (pval < 0) {
    *bufpp = memcpya(*bufpp, "        NA ", 11);
  } else if (pval <= output_min_p) {
    *bufpp = memcpya(*bufpp, output_min_logp_str, output_min_logp_strlen);
  } else if (pval < 1) {
    *bufpp = dtoa_g_wxp4x(-log10(pval), 10, ' ', *bufpp);
  } else {
    *bufpp = memcpya(*bufpp, "         0 ", 11);
  }
}

int32_t multcomp(char* outname, char* outname_end, uint32_t* marker_uidxs, uintptr_t chi_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, double* chi, double pfilter, double output_min_p, uint32_t mtest_adjust, uint32_t skip_gc, double adjust_lambda, uint32_t* tcnt, double* pvals) {
  // Association statistics can be provided in three ways:
  // 1. Just p-values (pvals[]).
  // 2. T statistics (in chi[]) and dfs (in tcnt[]).
  // 3. 1df chi-square stats (in chi[]).
  unsigned char* bigstack_mark = g_bigstack_base;
  uint32_t is_log10 = mtest_adjust & ADJUST_LOG10;
  uint32_t qq_plot = mtest_adjust & ADJUST_QQ;
  FILE* outfile = nullptr;
  double pv_holm = 0.0;
  double pv_sidak_sd = 0;
  int32_t retval = 0;
  uint32_t is_set_test = !plink_maxsnp;
  uint32_t adjust_gc = (mtest_adjust & ADJUST_GC) && (!skip_gc);
  uint32_t output_min_p_strlen = 11;
  uint32_t uii = 0;
  uint32_t* new_tcnt = nullptr;
  double* unadj = nullptr;
  char output_min_p_str[16];
  uint32_t pct;
  double* sp;
  double* schi;
  double* pv_bh;
  double* pv_by;
  uint32_t* new_order;
  uint32_t cur_idx;
  uintptr_t marker_uidx;
  double dxx;
  double dyy;
  double dzz;
  double harmonic_sum;
  double dct;
  double pval;
  double unadj_pval;
  double* pv_gc;
  double lambda_recip;
  double bonf;
  double pv_sidak_ss;
  char* bufptr;
  uint32_t ujj;
  uint32_t loop_end;

  if (bigstack_alloc_d(chi_ct, &sp) ||
      bigstack_alloc_d(chi_ct, &schi) ||
      bigstack_alloc_ui(chi_ct, &new_order)) {
    goto multcomp_ret_NOMEM;
  }
  if (pvals) {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      dxx = pvals[cur_idx];
      dyy = inverse_chiprob(dxx, 1);
      if (dyy >= 0) {
	sp[uii] = dxx;
	new_order[uii] = marker_uidxs[cur_idx];
	schi[uii] = dyy;
	uii++;
      }
    }
  } else if (tcnt) {
    if (bigstack_alloc_ui(chi_ct, &new_tcnt)) {
      goto multcomp_ret_NOMEM;
    }
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      ujj = tcnt[cur_idx];
      if (ujj) {
	dxx = chi[cur_idx]; // not actually squared
	dyy = calc_tprob(dxx, ujj);
	if (dyy > -1) {
	  sp[uii] = dyy;
	  new_order[uii] = marker_uidxs[cur_idx];
	  schi[uii] = dxx * dxx;
	  new_tcnt[uii] = ujj;
	  uii++;
	}
      }
    }
  } else {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      dxx = chi[cur_idx];
      if (dxx >= 0) {
	dyy = chiprob_p(dxx, 1);
	if (dyy > -1) {
	  sp[uii] = dyy;
	  new_order[uii] = marker_uidxs[cur_idx];
	  schi[uii] = dxx;
	  uii++;
	}
      }
    }
  }
  chi_ct = uii;
  if (!chi_ct) {
    logprint("Zero valid tests; --adjust skipped.\n");
    goto multcomp_ret_1;
  }
  if (qsort_ext((char*)sp, chi_ct, sizeof(double), double_cmp_deref, (char*)new_order, sizeof(int32_t))) {
    goto multcomp_ret_NOMEM;
  }
  if (tcnt) {
    if (qsort_ext((char*)schi, chi_ct, sizeof(double), double_cmp_deref, (char*)new_tcnt, sizeof(int32_t))) {
      goto multcomp_ret_NOMEM;
    }
  } else {
#ifdef __cplusplus
    std::sort(schi, &(schi[chi_ct]));
#else
    qsort(schi, chi_ct, sizeof(double), double_cmp);
#endif
  }
  dct = chi_ct;

  // get lambda...
  if (skip_gc) {
    lambda_recip = 1.0;
  } else if (mtest_adjust & ADJUST_LAMBDA) {
    lambda_recip = adjust_lambda;
  } else {
    if (chi_ct & 1) {
      lambda_recip = schi[(chi_ct - 1) / 2];
    } else {
      lambda_recip = (schi[chi_ct / 2 - 1] + schi[chi_ct / 2]) * 0.5;
    }
    lambda_recip = lambda_recip / 0.456;
    if (lambda_recip < 1.0) {
      lambda_recip = 1.0;
    }
    LOGPRINTF("--adjust: Genomic inflation est. lambda (based on median chisq) = %g.\n", lambda_recip);
  }
  // ...now take the reciprocal (bugfix: forgot to do this with --lambda)
  if (lambda_recip > 1.0) {
    lambda_recip = 1.0 / lambda_recip;
  }

  // handle reverse-order calculations
  if (bigstack_alloc_d(chi_ct, &pv_bh) ||
      bigstack_alloc_d(chi_ct, &pv_by) ||
      bigstack_alloc_d(chi_ct, &pv_gc)) {
    goto multcomp_ret_NOMEM;
  }
  if (adjust_gc) {
    unadj = sp;
    sp = pv_gc;
  }
  uii = chi_ct;
  if (tcnt) {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      uii--;
      pv_gc[cur_idx] = calc_tprob(sqrt(schi[uii] * lambda_recip), new_tcnt[uii]);
    }
  } else {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      pv_gc[cur_idx] = chiprob_p(schi[--uii] * lambda_recip, 1);
    }
  }

  dyy = sp[chi_ct - 1];
  pv_bh[chi_ct - 1] = dyy;
  harmonic_sum = 1.0;
  for (cur_idx = chi_ct - 1; cur_idx > 0; cur_idx--) {
    dzz = dct / ((double)((int32_t)cur_idx));
    harmonic_sum += dzz;
    dxx = dzz * sp[cur_idx - 1];
    if (dyy > dxx) {
      dyy = dxx;
    }
    pv_bh[cur_idx - 1] = dyy;
  }

  dzz = 1.0 / dct;
  harmonic_sum *= dzz;

  dyy = harmonic_sum * sp[chi_ct - 1];
  if (dyy >= 1) {
    dyy = 1;
  }
  pv_by[chi_ct - 1] = dyy;
  harmonic_sum *= dct;
  for (cur_idx = chi_ct - 1; cur_idx > 0; cur_idx--) {
    dxx = (harmonic_sum / ((double)((int32_t)cur_idx))) * sp[cur_idx - 1];
    if (dyy > dxx) {
      dyy = dxx;
    }
    pv_by[cur_idx - 1] = dyy;
  }

  uii = strlen(outname_end);
  memcpy(&(outname_end[uii]), ".adjusted", 10);
  if (fopen_checked(outname, "w", &outfile)) {
    goto multcomp_ret_OPEN_FAIL;
  }
  if (!is_set_test) {
    sprintf(g_textbuf, " CHR %%%us      UNADJ %s", plink_maxsnp, skip_gc? "" : "        GC ");
    fprintf(outfile, g_textbuf, "SNP");
  } else {
    plink_maxsnp = max_marker_id_len - 1;
    if (plink_maxsnp < 3) {
      plink_maxsnp = 3;
    }
    sprintf(g_textbuf, " %%%us      UNADJ ", plink_maxsnp);
    fprintf(outfile, g_textbuf, "SET");
  }
  if (qq_plot) {
    fputs("        QQ ", outfile);
  }
  if (fputs_checked("      BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY\n", outfile)) {
    goto multcomp_ret_WRITE_FAIL;
  }
  fputs("0%", stdout);
  fflush(stdout);
  cur_idx = 0;
  if (!is_log10) {
    if (output_min_p == 0.0) {
      memcpy(output_min_p_str, "       INF ", 11);
    } else {
      bufptr = dtoa_g_wxp4x(output_min_p, 10, ' ', output_min_p_str);
      output_min_p_strlen = (uintptr_t)(bufptr - output_min_p_str);
    }
  } else {
    if (output_min_p == 0.0) {
      memcpy(output_min_p_str, "       INF ", 11);
    } else {
      bufptr = dtoa_g_wxp4x(-log10(output_min_p), 10, ' ', output_min_p_str);
      output_min_p_strlen = (uintptr_t)(bufptr - output_min_p_str);
    }
  }
  for (pct = 1; pct <= 100; pct++) {
    loop_end = (((uint64_t)pct) * chi_ct) / 100LLU;
    for (; cur_idx < loop_end; cur_idx++) {
      pval = sp[cur_idx];
      // if --pfilter specified, filter out both nan and negative pvals, since
      // both are currently used by upstream functions
      if ((pfilter != 2.0) && ((!(pval >= 0.0)) || (pval > pfilter))) {
	continue;
      }
      if (adjust_gc) {
        unadj_pval = unadj[cur_idx];
      } else {
	unadj_pval = pval;
      }
      marker_uidx = new_order[cur_idx];
      if (!is_set_test) {
        bufptr = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, get_variant_chrom(chrom_info_ptr, marker_uidx), g_textbuf));
      } else {
        bufptr = g_textbuf;
      }
      *bufptr++ = ' ';
      bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
      *bufptr++ = ' ';
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto multcomp_ret_WRITE_FAIL;
      }
      bonf = pval * dct;
      if (bonf > 1) {
	bonf = 1;
      }
      if (pv_holm < 1) {
	dyy = (chi_ct - cur_idx) * pval;
	if (dyy > 1) {
	  pv_holm = 1;
	} else if (pv_holm < dyy) {
	  pv_holm = dyy;
	}
      }
      // avoid catastrophic cancellation for small p-values
      // 1 - (1-p)^c = 1 - e^{c log(1-p)}
      // 2^{-7} threshold is arbitrary
      if (pval >= 0.0078125) {
	pv_sidak_ss = 1 - pow(1 - pval, dct);
	dyy = 1 - pow(1 - pval, dct - ((double)((int32_t)cur_idx)));
      } else {
	pv_sidak_ss = 1 - exp(dct * log1p(-pval));
	dyy = dct - (double)((int32_t)cur_idx);
	dyy = 1 - exp(dyy * log1p(-pval));
      }
      if (pv_sidak_sd < dyy) {
	pv_sidak_sd = dyy;
      }

      bufptr = g_textbuf;
      if (!is_log10) {
	adjust_print(unadj_pval, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	if (!skip_gc) {
	  adjust_print(pv_gc[cur_idx], output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	}
	if (qq_plot) {
	  bufptr = dtoa_g_wxp4x((((double)((int32_t)cur_idx)) + 0.5) * dzz, 10, ' ', bufptr);
	}
	adjust_print(bonf, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print(pv_holm, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print(pv_sidak_ss, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print(pv_sidak_sd, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print(pv_bh[cur_idx], output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print(pv_by[cur_idx], output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
      } else {
	adjust_print_log10(pval, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	if (!is_set_test) {
	  adjust_print_log10(pv_gc[cur_idx], output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	}
	if (qq_plot) {
          // quasi-bugfix (23 Mar 2018): this should be logscale, both for
          // consistency with plink 1.07 and because it makes more sense
	  bufptr = dtoa_g_wxp4x(-log10((((double)((int32_t)cur_idx)) + 0.5) * dzz), 10, ' ', bufptr);
	}
	adjust_print_log10(bonf, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print_log10(pv_holm, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print_log10(pv_sidak_ss, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print_log10(pv_sidak_sd, output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print_log10(pv_bh[cur_idx], output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
	adjust_print_log10(pv_by[cur_idx], output_min_p, output_min_p_str, output_min_p_strlen, &bufptr);
      }
      *bufptr++ = '\n';
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto multcomp_ret_WRITE_FAIL;
      }
    }
    if (pct < 100) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  fputs("\b\b\b", stdout);
  LOGPRINTFWW("--adjust values (%" PRIuPTR " %s%s) written to %s .\n", chi_ct, is_set_test? "nonempty set" : "variant", (chi_ct == 1)? "" : "s", outname);

  while (0) {
  multcomp_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  multcomp_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  multcomp_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 multcomp_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

char* model_assoc_tna(uint32_t model_fisher, char* wptr) {
  // write terminal NAs to buffer
  if (model_fisher) {
    return memcpya(wptr, "          NA\n", 13);
  } else {
    return memcpya(wptr, "          NA   NA           NA\n", 31);
  }
}

void calc_git(uint32_t pheno_nm_ct, uint32_t perm_vec_ct, uintptr_t* __restrict__ loadbuf, uint32_t* perm_vecst, uint32_t* results_bufs, uint32_t* thread_wkspace) {
  // Brian Browning's genotype indexing algorithm for low-MAF (and low missing
  // frequency) markers.
  // We accelerate it by using a special interleaved permutation representation
  // which supports vector addition without occupying extra space: see
  // generate_cc_perm_vec().  Counting the number of e.g. case heterozygote
  // genotypes across all permutations then proceeds as follows:
  // 1. For the first 15 heterozygote samples, just use 4-bit accumulators.
  //    This allows the inner loop to increment 32 counters simultaneously.
  // 2. Right before they'd otherwise be at risk of overflowing, we unfold the
  //    4-bit accumulators into a larger buffer of 8-bit accumulators.  Then we
  //    zero out the 4-bit accumulators, and restart the inner loop.
  // 3. This can happen up to 17 times before the 8-bit accumulators risk
  //    overflow.  Then, they are unfolded into the final output array of
  //    32-bit ints, zeroed out, and the second loop restarts.
  // Note that results_bufs[] is assumed to be zeroed out before this function
  // is called.
  uint32_t pheno_nm_ctl2x = QUATERCT_TO_WORDCT(pheno_nm_ct);
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct128x4 = perm_ct128 * 4;
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct16x4 = 4 * perm_ct16;
  __m128i* permsv = (__m128i*)perm_vecst;
  __m128i* gitv[9];
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct32x4 = perm_ct32 * 4;
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t perm_ct16x16 = 16 * perm_ct16;
  uintptr_t* permsv = (uintptr_t*)perm_vecst;
  uintptr_t* gitv[9];
#endif
  uint32_t cur_cts[3];
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t sample_type;
#ifdef __LP64__
  // 4- and 8-bit partial counts
  gitv[0] = (__m128i*)thread_wkspace;
  gitv[1] = &(((__m128i*)thread_wkspace)[perm_ct128x4]);
  gitv[2] = &(((__m128i*)thread_wkspace)[2 * perm_ct128x4]);
  gitv[3] = &(((__m128i*)thread_wkspace)[3 * perm_ct128x4]);
  gitv[4] = &(((__m128i*)thread_wkspace)[3 * perm_ct128x4 + 2 * perm_ct32]);
  gitv[5] = &(((__m128i*)thread_wkspace)[3 * perm_ct128x4 + 4 * perm_ct32]);
  gitv[6] = &(((__m128i*)results_bufs)[2 * perm_ct16x4]);
  gitv[7] = &(((__m128i*)results_bufs)[perm_ct16x4]);
  gitv[8] = (__m128i*)results_bufs;
#else
  gitv[0] = (uintptr_t*)thread_wkspace;
  gitv[1] = (uintptr_t*)(&(thread_wkspace[perm_ct32x4]));
  gitv[2] = (uintptr_t*)(&(thread_wkspace[2 * perm_ct32x4]));
  gitv[3] = (uintptr_t*)(&(thread_wkspace[3 * perm_ct32x4]));
  gitv[4] = (uintptr_t*)(&(thread_wkspace[3 * perm_ct32x4 + 2 * perm_ct8]));
  gitv[5] = (uintptr_t*)(&(thread_wkspace[3 * perm_ct32x4 + 4 * perm_ct8]));
  gitv[6] = (uintptr_t*)(&(results_bufs[2 * perm_ct16x16]));
  gitv[7] = (uintptr_t*)(&(results_bufs[perm_ct16x16]));
  gitv[8] = (uintptr_t*)results_bufs;
#endif
  cur_cts[0] = 0;
  cur_cts[1] = 0;
  cur_cts[2] = 0;
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulii = ~(*loadbuf++);
    if (uii + 1 == pheno_nm_ctl2x) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulii &= (ONELU << (ujj * 2)) - ONELU;
      }
    }
    while (ulii) {
      ujj = CTZLU(ulii) & (BITCT - 2); // get pos of next non-[hom A2] sample
      sample_type = ((ulii >> ujj) & 3) - 1;
      ukk = cur_cts[sample_type] + 1;
      cur_cts[sample_type] = ukk;
#ifdef __LP64__
      unroll_incr_1_4(&(permsv[(ujj / 2) * perm_ct128]), gitv[sample_type], perm_ct128);
      if (!(ukk % 15)) {
	unroll_zero_incr_4_8(gitv[sample_type], gitv[sample_type + 3], perm_ct32);
	if (!(ukk % 255)) {
	  unroll_zero_incr_8_32(gitv[sample_type + 3], gitv[sample_type + 6], perm_ct16);
	}
      }
#else
      unroll_incr_1_4(&(permsv[(ujj / 2) * perm_ct32]), gitv[sample_type], perm_ct32);
      if (!(ukk % 15)) {
	unroll_zero_incr_4_8(gitv[sample_type], gitv[sample_type + 3], perm_ct8);
	if (!(ukk % 255)) {
	  unroll_zero_incr_8_32(gitv[sample_type + 3], gitv[sample_type + 6], perm_ct4);
	}
      }
#endif
      ulii &= ~((3 * ONELU) << ujj);
    }
#ifdef __LP64__
    permsv = &(permsv[BITCT2 * perm_ct128]);
#else
    permsv = &(permsv[BITCT2 * perm_ct32]);
#endif
  }
  for (sample_type = 0; sample_type < 3; sample_type++) {
    uii = cur_cts[sample_type];
#ifdef __LP64__
    if (uii % 15) {
      unroll_incr_4_8(gitv[sample_type], gitv[sample_type + 3], perm_ct32);
    }
    if (uii % 255) {
      unroll_incr_8_32(gitv[sample_type + 3], gitv[sample_type + 6], perm_ct16);
    }
#else
    if (uii % 15) {
      unroll_incr_4_8(gitv[sample_type], gitv[sample_type + 3], perm_ct8);
    }
    if (uii % 255) {
      unroll_incr_8_32(gitv[sample_type + 3], gitv[sample_type + 6], perm_ct4);
    }
#endif
  }
}

void calc_qgit(uint32_t pheno_nm_ct, uintptr_t perm_vec_ctcl8m, uint32_t num_perms_now, uintptr_t* __restrict__ loadbuf, double* perm_vecstd, double* thread_bufs) {
  uint32_t pheno_nm_ctl2x = QUATERCT_TO_WORDCT(pheno_nm_ct);
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the sample idx
  uint32_t row_mult = perm_vec_ctcl8m / 4;

  uint32_t loop_len = (num_perms_now + 1) / 2;
  __m128d* permsv = (__m128d*)perm_vecstd;
  __m128d* __restrict__ perm_readv;
  __m128d* __restrict__ git_writev;
  __m128d* __restrict__ git_write2v;
  __m128d vxx;
#else
  uint32_t row_mult = perm_vec_ctcl8m / 2;
  double* __restrict__ perm_read;
  double* __restrict__ git_write;
  double* __restrict__ git_write2;
  double dxx;
#endif
  uintptr_t ulii;
  uint32_t sample_type;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulii = ~(*loadbuf++);
    if (uii + 1 == pheno_nm_ctl2x) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulii &= (ONELU << (ujj * 2)) - ONELU;
      }
    }
    while (ulii) {
      ujj = CTZLU(ulii) & (BITCT - 2);
      sample_type = (ulii >> ujj) & 3;
#ifdef __LP64__
      // note that the gain from using SSE2 for double-precision arithmetic is
      // typically minimal because modern cores tend to have two FPUs, so we
      // should only use it opportunistically.  it's painless here, though.
      perm_readv = &(permsv[ujj * row_mult]);
      if (sample_type == 1) {
	git_writev = (__m128d*)thread_bufs;
	for (ukk = 0; ukk < loop_len; ukk++) {
	  *git_writev = _mm_add_pd(*git_writev, *perm_readv++);
	  git_writev++;
	}
      } else if (sample_type == 3) {
	// hom rare
	git_writev = (__m128d*)thread_bufs;
	for (ukk = 0; ukk < loop_len; ukk++) {
	  vxx = *perm_readv++;
	  *git_writev = _mm_add_pd(*git_writev, _mm_add_pd(vxx, vxx));
	  git_writev++;
	}
      } else {
	// missing
	git_writev = (__m128d*)(&(thread_bufs[perm_vec_ctcl8m]));
	git_write2v = (__m128d*)(&(thread_bufs[2 * perm_vec_ctcl8m]));
	for (ukk = 0; ukk < loop_len; ukk++) {
	  vxx = *perm_readv++;
	  *git_writev = _mm_add_pd(*git_writev, vxx);
	  git_writev++;
	  *git_write2v = _mm_add_pd(*git_write2v, _mm_mul_pd(vxx, vxx));
	  git_write2v++;
	}
      }
#else
      perm_read = &(perm_vecstd[ujj * row_mult]);
      if (sample_type == 1) {
	git_write = thread_bufs;
	for (ukk = 0; ukk < num_perms_now; ukk++) {
	  *git_write += *perm_read++;
	  git_write++;
	}
      } else if (sample_type == 3) {
	git_write = thread_bufs;
	for (ukk = 0; ukk < num_perms_now; ukk++) {
	  dxx = *perm_read++;
	  *git_write += dxx * 2;
	  git_write++;
	}
      } else {
	git_write = &(thread_bufs[perm_vec_ctcl8m]);
	git_write2 = &(thread_bufs[2 * perm_vec_ctcl8m]);
	for (ukk = 0; ukk < num_perms_now; ukk++) {
	  dxx = *perm_read++;
	  *git_write += dxx;
	  git_write++;
	  *git_write2 += dxx * dxx;
	  git_write2++;
	}
      }
#endif
      ulii &= ~((3 * ONELU) << ujj);
    }
#ifdef __LP64__
    permsv = &(permsv[(BITCT2 / 2) * perm_vec_ctcl8m]);
#else
    perm_vecstd = &(perm_vecstd[BITCT2 * perm_vec_ctcl8m]);
#endif
  }
}

void calc_qgit_lin(uint32_t pheno_nm_ct, uintptr_t perm_vec_ctcl8m, uint32_t num_perms_now, uintptr_t* __restrict__ loadbuf, double* perm_vecstd, double* thread_bufs) {
  uint32_t pheno_nm_ctl2x = QUATERCT_TO_WORDCT(pheno_nm_ct);
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the sample idx
  uint32_t row_mult = perm_vec_ctcl8m / 4;

  uint32_t loop_len = (num_perms_now + 1) / 2;
  __m128d* permsv = (__m128d*)perm_vecstd;
  __m128d* __restrict__ perm_readv;
  __m128d* __restrict__ git_writev;
  __m128d* __restrict__ git_write2v;
  __m128d vxx;
#else
  uint32_t row_mult = perm_vec_ctcl8m / 2;
  double* __restrict__ perm_read;
  double* __restrict__ git_write;
  double* __restrict__ git_write2;
  double dxx;
#endif
  uintptr_t ulii;
  uint32_t sample_type;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulii = ~(*loadbuf++);
    if (uii + 1 == pheno_nm_ctl2x) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulii &= (ONELU << (ujj * 2)) - ONELU;
      }
    }
    while (ulii) {
      ujj = CTZLU(ulii) & (BITCT - 2);
      sample_type = (ulii >> ujj) & 3;
#ifdef __LP64__
      perm_readv = &(permsv[ujj * row_mult]);
      if (sample_type == 1) {
	git_writev = (__m128d*)thread_bufs;
	git_write2v = (__m128d*)(&(thread_bufs[perm_vec_ctcl8m]));
      } else if (sample_type == 3) {
	// hom rare
	git_writev = (__m128d*)(&(thread_bufs[2 * perm_vec_ctcl8m]));
	git_write2v = (__m128d*)(&(thread_bufs[3 * perm_vec_ctcl8m]));
      } else {
	// missing
	git_writev = (__m128d*)(&(thread_bufs[4 * perm_vec_ctcl8m]));
	git_write2v = (__m128d*)(&(thread_bufs[5 * perm_vec_ctcl8m]));
      }
      for (ukk = 0; ukk < loop_len; ukk++) {
	vxx = *perm_readv++;
	*git_writev = _mm_add_pd(*git_writev, vxx);
	git_writev++;
	*git_write2v = _mm_add_pd(*git_write2v, _mm_mul_pd(vxx, vxx));
	git_write2v++;
      }
#else
      perm_read = &(perm_vecstd[ujj * row_mult]);
      if (sample_type == 1) {
	git_write = thread_bufs;
	git_write2 = &(thread_bufs[perm_vec_ctcl8m]);
      } else if (sample_type == 3) {
	git_write = &(thread_bufs[2 * perm_vec_ctcl8m]);
	git_write2 = &(thread_bufs[3 * perm_vec_ctcl8m]);
      } else {
	git_write = &(thread_bufs[4 * perm_vec_ctcl8m]);
	git_write2 = &(thread_bufs[5 * perm_vec_ctcl8m]);
      }
      for (ukk = 0; ukk < num_perms_now; ukk++) {
	dxx = *perm_read++;
	*git_write += dxx;
	git_write++;
	*git_write2 += dxx * dxx;
	git_write2++;
      }
#endif
      ulii &= ~((3 * ONELU) << ujj);
    }
#ifdef __LP64__
    permsv = &(permsv[(BITCT2 / 2) * perm_vec_ctcl8m]);
#else
    perm_vecstd = &(perm_vecstd[BITCT2 * perm_vec_ctcl8m]);
#endif
  }
}

#ifdef __LP64__
uintptr_t rem_cost_60v(__m128i* vec1, __m128i* vend, __m128i* vec2) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i xor_vec;
  __m128i detect_homcom;
  __m128i result_a;
  __m128i acc_a;
  __m128i acc_b;
  __univec acc;
  acc.vi = _mm_setzero_si128();
  do {
    loader = *vec1++;
    loader2 = *vec2++;
    xor_vec = _mm_xor_si128(loader, loader2);
    detect_homcom = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(loader, 1), loader), _mm_and_si128(_mm_srli_epi64(loader2, 1), loader2));
    acc_a = _mm_and_si128(_mm_or_si128(xor_vec, _mm_srli_epi64(xor_vec, 1)), m1);
    acc_b = _mm_andnot_si128(detect_homcom, acc_a);

    loader = *vec1++;
    loader2 = *vec2++;
    xor_vec = _mm_xor_si128(loader, loader2);
    detect_homcom = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(loader, 1), loader), _mm_and_si128(_mm_srli_epi64(loader2, 1), loader2));
    result_a = _mm_and_si128(_mm_or_si128(xor_vec, _mm_srli_epi64(xor_vec, 1)), m1);
    acc_a = _mm_add_epi64(acc_a, result_a);
    acc_b = _mm_add_epi64(acc_b, _mm_andnot_si128(detect_homcom, result_a));

    loader = *vec1++;
    loader2 = *vec2++;
    xor_vec = _mm_xor_si128(loader, loader2);
    detect_homcom = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(loader, 1), loader), _mm_and_si128(_mm_srli_epi64(loader2, 1), loader2));
    result_a = _mm_and_si128(_mm_or_si128(xor_vec, _mm_srli_epi64(xor_vec, 1)), m1);
    acc_a = _mm_add_epi64(acc_a, result_a);
    acc_b = _mm_add_epi64(acc_b, _mm_andnot_si128(detect_homcom, result_a));
    acc_a = _mm_add_epi64(_mm_and_si128(acc_a, m2), _mm_and_si128(_mm_srli_epi64(acc_a, 2), m2));
    acc_a = _mm_add_epi64(acc_a, _mm_add_epi64(_mm_and_si128(acc_b, m2), _mm_and_si128(_mm_srli_epi64(acc_b, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(acc_a, m4), _mm_and_si128(_mm_srli_epi64(acc_a, 4), m4)));
  } while (vec1 < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  return ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
}

uintptr_t qrem_cost2_40v(__m128i* vec1, __m128i* vend, __m128i* vec2) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i xor_vec;
  __m128i detect_missing;
  __m128i result_a;
  __m128i result_b;
  __m128i result_c;
  __m128i inner_acc;
  __univec acc;
  acc.vi = _mm_setzero_si128();
  do {
    loader = *vec1++;
    loader2 = *vec2++;
    xor_vec = _mm_xor_si128(loader, loader2);
    detect_missing = _mm_or_si128(_mm_andnot_si128(_mm_srli_epi64(loader, 1), loader), _mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2));
    result_a = _mm_and_si128(_mm_or_si128(xor_vec, _mm_srli_epi64(xor_vec, 1)), m1);
    result_b = _mm_and_si128(result_a, detect_missing);
    inner_acc = _mm_and_si128(result_b, xor_vec);
    inner_acc = _mm_add_epi64(_mm_add_epi64(result_a, result_b), inner_acc);
    inner_acc = _mm_add_epi64(_mm_and_si128(inner_acc, m2), _mm_and_si128(_mm_srli_epi64(inner_acc, 2), m2));
    loader = *vec1++;
    loader2 = *vec2++;
    xor_vec = _mm_xor_si128(loader, loader2);
    detect_missing = _mm_or_si128(_mm_andnot_si128(_mm_srli_epi64(loader, 1), loader), _mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2));
    result_a = _mm_and_si128(_mm_or_si128(xor_vec, _mm_srli_epi64(xor_vec, 1)), m1);
    result_b = _mm_and_si128(result_a, detect_missing);
    result_c = _mm_and_si128(result_b, xor_vec);
    result_c = _mm_add_epi64(_mm_add_epi64(result_a, result_b), result_c);
    inner_acc = _mm_add_epi64(inner_acc, _mm_add_epi64(_mm_and_si128(result_c, m2), _mm_and_si128(_mm_srli_epi64(result_c, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(inner_acc, m4), _mm_and_si128(_mm_srli_epi64(inner_acc, 4), m4)));
  } while (vec1 < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  return ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
uintptr_t rem_cost_6(uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
  uintptr_t loader = *loadbuf1++;
  uintptr_t loader2 = *loadbuf2++;
  uintptr_t xor_word = loader ^ loader2;
  uintptr_t detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
  uintptr_t acc_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  uintptr_t acc_b = acc_a & (~detect_homcom);
  uintptr_t result_a;
  uintptr_t acc;

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  acc_a += result_a;
  acc_b += result_a & (~detect_homcom);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  acc_a += result_a;
  acc_b += result_a & (~detect_homcom);
  acc_a = (acc_a & 0x33333333) + ((acc_a >> 2) & 0x33333333);
  acc_a += (acc_b & 0x33333333) + ((acc_b >> 2) & 0x33333333);
  acc = (acc_a & 0x0f0f0f0f) + ((acc_a >> 4) & 0x0f0f0f0f);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
  acc_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  acc_b = acc_a & (~detect_homcom);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  acc_a += result_a;
  acc_b += result_a & (~detect_homcom);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  acc_a += result_a;
  acc_b += result_a & (~detect_homcom);
  acc_a = (acc_a & 0x33333333) + ((acc_a >> 2) & 0x33333333);
  acc_a += (acc_b & 0x33333333) + ((acc_b >> 2) & 0x33333333);
  acc += (acc_a & 0x0f0f0f0f) + ((acc_a >> 4) & 0x0f0f0f0f);
  return (acc * 0x01010101) >> 24;
}

uintptr_t qrem_cost2_4(uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
  uintptr_t loader = *loadbuf1++;
  uintptr_t loader2 = *loadbuf2++;
  uintptr_t xor_word = loader ^ loader2;
  uintptr_t detect_missing = (loader & (~(loader >> 1))) | (loader2 & (~(loader2 >> 1)));
  uintptr_t result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  uintptr_t result_b = result_a & detect_missing;
  uintptr_t inner_acc = result_b & xor_word;
  uintptr_t result_c;
  uintptr_t acc;
  inner_acc += result_a + result_b;
  inner_acc = (inner_acc & 0x33333333) + ((inner_acc >> 2) & 0x33333333);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_missing = (loader & (~(loader >> 1))) | (loader2 & (~(loader2 >> 1)));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  result_b = result_a & detect_missing;
  result_c = result_b & xor_word;
  result_c += result_a + result_b;
  inner_acc += (result_c & 0x33333333) + ((result_c >> 2) & 0x33333333);
  acc = (inner_acc & 0x0f0f0f0f) + ((inner_acc >> 4) & 0x0f0f0f0f);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_missing = (loader & (~(loader >> 1))) | (loader2 & (~(loader2 >> 1)));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  result_b = result_a & detect_missing;
  inner_acc = result_b & xor_word;
  inner_acc += result_a + result_b;
  inner_acc = (inner_acc & 0x33333333) + ((inner_acc >> 2) & 0x33333333);

  loader = *loadbuf1++;
  loader2 = *loadbuf2++;
  xor_word = loader & loader2;
  detect_missing = (loader & (~(loader >> 1))) | (loader2 & (~(loader2 >> 1)));
  result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
  result_b = result_a & detect_missing;
  result_c = result_b & xor_word;
  result_c += result_a + result_b;
  inner_acc += (result_c & 0x33333333) + ((result_c >> 2) & 0x33333333);
  acc += (inner_acc & 0x0f0f0f0f) + ((inner_acc >> 4) & 0x0f0f0f0f);
  return (acc * 0x01010101) >> 24;
}
#endif

uintptr_t rem_cost(uintptr_t sample_ctv2, uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
  // Cost: 2 * (<-> neither side homcom) + (<-> homcom)
  //
  // We can efficiently calculate this as follows:
  //   xor = vec1 ^ vec2
  //   detect_homcom = (vec1 & (vec1 >> 1)) | (vec2 & (vec2 >> 1))
  //   A := (xor | (xor >> 1)) & 0x5555...
  //   B := A & (~detect_homcom)
  //   cost += popcount2(A + B)
  uintptr_t* lptr_end = &(loadbuf1[sample_ctv2]);
  uintptr_t cost = 0;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t xor_word;
  uintptr_t detect_homcom;
  uintptr_t result_a;
  uintptr_t result_b;
#ifdef __LP64__
  uintptr_t cur_decr = 60;
  uintptr_t* lptr_6x_end;
  sample_ctv2 -= sample_ctv2 % 6;
  while (sample_ctv2 >= 60) {
  rem_cost_loop:
    lptr_6x_end = &(loadbuf1[cur_decr]);
    cost += rem_cost_60v((__m128i*)loadbuf1, (__m128i*)lptr_6x_end, (__m128i*)loadbuf2);
    loadbuf1 = lptr_6x_end;
    loadbuf2 = &(loadbuf2[cur_decr]);
    sample_ctv2 -= cur_decr;
  }
  if (sample_ctv2) {
    cur_decr = sample_ctv2;
    goto rem_cost_loop;
  }
#else
  uintptr_t* lptr_six_end = &(loadbuf1[sample_ctv2 - (sample_ctv2 % 6)]);
  while (loadbuf1 < lptr_six_end) {
    cost += rem_cost_6(loadbuf1, loadbuf2);
    loadbuf1 = &(loadbuf1[6]);
    loadbuf2 = &(loadbuf2[6]);
  }
#endif
  while (loadbuf1 < lptr_end) {
    loader = *loadbuf1++;
    loader2 = *loadbuf2++;
    xor_word = loader ^ loader2;
    detect_homcom = (loader & (loader >> 1)) | (loader2 & (loader2 >> 1));
    result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
    result_b = result_a & (~detect_homcom);
    cost += popcount2_long(result_a + result_b);
  }
  return cost;
}

uintptr_t qrem_cost2(uintptr_t sample_ctl2, uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
  // Cost: 3 + 3 * (missing <-> homrar/het) + 2 * (missing <-> homcom) +
  //       (homrar <-> het/homcom) + (het <-> homcom)
  //
  // xor 01: 3 if 00-01, 1 of 10-11
  // xor 10: 2 if 01-11, 1 if 00-10
  // xor 11: 3 if 01-10, 1 if 00-11
  //
  // We can efficiently calculate this as follows:
  //   xor = vec1 ^ vec2
  //   detect_missing = (vec1 & (~(vec1 >> 1))) | (vec2 & (~(vec2 >> 1)))
  //   A := (xor | (xor >> 1)) & 0x5555...
  //   B := A & detect_missing
  //   C := B & xor
  //   cost += popcount2(A + B + C)
  // (I would not be surprised if a few operations could be shaved from this.)
  uintptr_t* lptr_end = &(loadbuf1[sample_ctl2]);
  uintptr_t cost = 3;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t xor_word;
  uintptr_t detect_missing;
  uintptr_t result_a;
  uintptr_t result_b;
  uintptr_t result_c;
#ifdef __LP64__
  uintptr_t cur_decr = 40;
  uintptr_t* lptr_4x_end;
  sample_ctl2 &= ~3LLU;
  while (sample_ctl2 >= 40) {
  qrem_cost2_loop:
    lptr_4x_end = &(loadbuf1[cur_decr]);
    cost += qrem_cost2_40v((__m128i*)loadbuf1, (__m128i*)lptr_4x_end, (__m128i*)loadbuf2);
    loadbuf1 = lptr_4x_end;
    loadbuf2 = &(loadbuf2[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto qrem_cost2_loop;
  }
#else
  uintptr_t* lptr_four_end = &(loadbuf1[sample_ctl2 & (~3)]);
  while (loadbuf1 < lptr_four_end) {
    cost += qrem_cost2_4(loadbuf1, loadbuf2);
    loadbuf1 = &(loadbuf1[4]);
    loadbuf2 = &(loadbuf2[4]);
  }
#endif
  while (loadbuf1 < lptr_end) {
    loader = *loadbuf1++;
    loader2 = *loadbuf2++;
    xor_word = loader ^ loader2;
    detect_missing = (loader & (~(loader >> 1))) | (loader2 & (~(loader2 >> 1)));
    result_a = (xor_word | (xor_word >> 1)) & FIVEMASK;
    result_b = result_a & detect_missing;
    result_c = result_b & xor_word;
    cost += popcount2_long(result_a + result_b + result_c);
  }
  return cost;
}

#ifdef __LP64__
static inline void calc_rem_merge4_two(uint32_t perm_ct128, __m128i* __restrict__ perm_ptr, __m128i* __restrict__ rem_merge4a, __m128i* __restrict__ rem_merge4b) {
  const __m128i m1x4 = {0x1111111111111111LLU, 0x1111111111111111LLU};
  __m128i loader;
  __m128i loader2;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct128; pbidx++) {
    loader = *perm_ptr++;
    loader2 = _mm_and_si128(loader, m1x4);
    rem_merge4a[0] = _mm_add_epi64(rem_merge4a[0], loader2);
    rem_merge4b[0] = _mm_add_epi64(rem_merge4b[0], loader2);
    loader = _mm_srli_epi64(loader, 1);
    loader2 = _mm_and_si128(loader, m1x4);
    rem_merge4a[1] = _mm_add_epi64(rem_merge4a[1], loader2);
    rem_merge4b[1] = _mm_add_epi64(rem_merge4b[1], loader2);
    loader = _mm_srli_epi64(loader, 1);
    loader2 = _mm_and_si128(loader, m1x4);
    rem_merge4a[2] = _mm_add_epi64(rem_merge4a[2], loader2);
    rem_merge4b[2] = _mm_add_epi64(rem_merge4b[2], loader2);
    loader = _mm_srli_epi64(loader, 1);
    loader2 = _mm_and_si128(loader, m1x4);
    rem_merge4a[3] = _mm_add_epi64(rem_merge4a[3], loader2);
    rem_merge4b[3] = _mm_add_epi64(rem_merge4b[3], loader2);
    rem_merge4a = &(rem_merge4a[4]);
    rem_merge4b = &(rem_merge4b[4]);
  }
}

static inline void calc_rem_merge32_minus(uint32_t perm_ct16, __m128i* __restrict__ rem_merge8, __m128i* rem_write) {
  // temporary integer underflow is possible here, but by the end of the
  // calculation it should be reversed
  const __m128i m8x32 = {0x000000ff000000ffLLU, 0x000000ff000000ffLLU};
  __m128i loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct16; pbidx++) {
    loader = *rem_merge8;
    rem_write[0] = _mm_sub_epi64(rem_write[0], _mm_and_si128(loader, m8x32));
    loader = _mm_srli_epi64(loader, 8);
    rem_write[1] = _mm_sub_epi64(rem_write[1], _mm_and_si128(loader, m8x32));
    loader = _mm_srli_epi64(loader, 8);
    rem_write[2] = _mm_sub_epi64(rem_write[2], _mm_and_si128(loader, m8x32));
    loader = _mm_srli_epi64(loader, 8);
    rem_write[3] = _mm_sub_epi64(rem_write[3], _mm_and_si128(loader, m8x32));
    rem_write = &(rem_write[4]);
    *rem_merge8++ = _mm_setzero_si128();
  }
}
#else
static inline void calc_rem_merge4_two(uint32_t perm_ct32, uintptr_t* __restrict__ perm_ptr, uintptr_t* __restrict__ rem_merge4a, uintptr_t* __restrict__ rem_merge4b) {
  uintptr_t loader;
  uintptr_t loader2;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct32; pbidx++) {
    loader = *perm_ptr++;
    loader2 = loader & 0x11111111;
    rem_merge4a[0] += loader2;
    rem_merge4b[0] += loader2;
    loader2 = (loader >> 1) & 0x11111111;
    rem_merge4a[1] += loader2;
    rem_merge4b[1] += loader2;
    loader2 = (loader >> 2) & 0x11111111;
    rem_merge4a[2] += loader2;
    rem_merge4b[2] += loader2;
    loader2 = (loader >> 3) & 0x11111111;
    rem_merge4a[3] += loader2;
    rem_merge4b[3] += loader2;
    rem_merge4a = &(rem_merge4a[4]);
    rem_merge4b = &(rem_merge4b[4]);
  }
}

static inline void calc_rem_merge32_minus(uint32_t perm_ct4, uintptr_t* __restrict__ rem_merge8, uintptr_t* __restrict__ rem_write) {
  uintptr_t loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct4; pbidx++) {
    loader = *rem_merge8;
    rem_write[0] -= (uint8_t)loader;
    loader >>= 8;
    rem_write[1] -= (uint8_t)loader;
    loader >>= 8;
    rem_write[2] -= (uint8_t)loader;
    loader >>= 8;
    rem_write[3] -= loader;
    rem_write = &(rem_write[4]);
    *rem_merge8++ = 0;
  }
}
#endif

void calc_rem(uint32_t pheno_nm_ct, uintptr_t perm_vec_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_ref, uint32_t* perm_vecst, uint32_t* results_bufs, uint32_t* thread_wkspace) {
  uint32_t pheno_nm_ctl2x = QUATERCT_TO_WORDCT(pheno_nm_ct);
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  // [cur_xor - 1][cur_raw]
  // low 8 bits give index of first remv[] array to increment; next 8 bits give
  // second index if nonzero, or indicate its absence
  const uint32_t idx_table[3][4] = {{0x300, 0x102, 4, 5}, {0x500, 2, 0x104, 3}, {0, 0x502, 0x304, 1}};
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct128x4 = perm_ct128 * 4;
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct16x4 = 4 * perm_ct16;
  __m128i* permsv = (__m128i*)perm_vecst;
  // 0, 2, 4: homrar, missing, het ct increment
  // 1, 3, 5: homrar, missing, het ct decrement
  __m128i* remv[15];
  __m128i* __restrict__ perm_ptr;
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct32x4 = perm_ct32 * 4;
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t perm_ct16x16 = 16 * perm_ct16;
  uintptr_t* permsv = (uintptr_t*)perm_vecst;
  uintptr_t* remv[15];
  uintptr_t* perm_ptr;
#endif

  uint32_t cur_cts[6];
  uintptr_t ulraw1;
  uintptr_t ulxor;
  uint32_t cur_xor;
  uint32_t cur_raw;
  uint32_t idx1;
  uint32_t idx2;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
#ifdef __LP64__
  for (uii = 0; uii < 6; uii++) {
    remv[uii] = &(((__m128i*)thread_wkspace)[uii * perm_ct128x4]);
  }
  for (uii = 0; uii < 6; uii++) {
    remv[uii + 6] = &(((__m128i*)thread_wkspace)[6 * perm_ct128x4 + 2 * uii * perm_ct32]);
  }
  remv[12] = (__m128i*)results_bufs;
  remv[13] = &(((__m128i*)results_bufs)[perm_ct16x4]);
  remv[14] = &(((__m128i*)results_bufs)[2 * perm_ct16x4]);
#else
  for (uii = 0; uii < 6; uii++) {
    remv[uii] = (uintptr_t*)(&(thread_wkspace[uii * perm_ct32x4]));
  }
  for (uii = 0; uii < 6; uii++) {
    remv[uii + 6] = (uintptr_t*)(&(thread_wkspace[6 * perm_ct32x4 + 2 * uii * perm_ct8]));
  }
  remv[12] = (uintptr_t*)results_bufs;
  remv[13] = (uintptr_t*)(&(results_bufs[perm_ct16x16]));
  remv[14] = (uintptr_t*)(&(results_bufs[2 * perm_ct16x16]));
#endif

  for (uii = 0; uii < 6; uii++) {
    cur_cts[uii] = 0;
  }
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulraw1 = *loadbuf++;
    ulxor = ulraw1 ^ (*loadbuf_ref++);
    if (uii + 1 == pheno_nm_ctl2x) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulxor &= (ONELU << (ujj * 2)) - ONELU;
      }
    }
    while (ulxor) {
      ujj = CTZLU(ulxor) & (BITCT - 2);
      cur_xor = (ulxor >> ujj) & 3;
      cur_raw = (ulraw1 >> ujj) & 3;
      idx1 = idx_table[cur_xor - 1][cur_raw];
      idx2 = idx1 >> 8;
      idx1 &= 255;
#ifdef __LP64__
      perm_ptr = &(permsv[(ujj / 2) * perm_ct128]);
      if (!idx2) {
	unroll_incr_1_4(perm_ptr, remv[idx1], perm_ct128);
      } else {
	calc_rem_merge4_two(perm_ct128, perm_ptr, remv[idx1], remv[idx2]);
	ukk = cur_cts[idx2] + 1;
	cur_cts[idx2] = ukk;
	if (!(ukk % 15)) {
	  unroll_zero_incr_4_8(remv[idx2], remv[idx2 + 6], perm_ct32);
	  if (!(ukk % 255)) {
	    calc_rem_merge32_minus(perm_ct16, remv[idx2 + 6], remv[(idx2 / 2) + 12]);
	  }
	}
      }
      ukk = cur_cts[idx1] + 1;
      cur_cts[idx1] = ukk;
      if (!(ukk % 15)) {
	unroll_zero_incr_4_8(remv[idx1], remv[idx1 + 6], perm_ct32);
	if (!(ukk % 255)) {
	  if (!(idx1 & 1)) {
	    unroll_zero_incr_8_32(remv[idx1 + 6], remv[(idx1 / 2) + 12], perm_ct16);
	  } else {
	    calc_rem_merge32_minus(perm_ct16, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
	  }
	}
      }
#else
      perm_ptr = &(permsv[(ujj / 2) * perm_ct32]);
      if (!idx2) {
	unroll_incr_1_4(perm_ptr, remv[idx1], perm_ct32);
      } else {
	calc_rem_merge4_two(perm_ct32, perm_ptr, remv[idx1], remv[idx2]);
	ukk = cur_cts[idx2] + 1;
	cur_cts[idx2] = ukk;
	if (!(ukk % 15)) {
	  unroll_zero_incr_4_8(remv[idx2], remv[idx2 + 6], perm_ct8);
	  if (!(ukk % 255)) {
	    calc_rem_merge32_minus(perm_ct4, remv[idx2 + 6], remv[(idx2 / 2) + 12]);
	  }
	}
      }
      ukk = cur_cts[idx1] + 1;
      cur_cts[idx1] = ukk;
      if (!(ukk % 15)) {
	unroll_zero_incr_4_8(remv[idx1], remv[idx1 + 6], perm_ct8);
	if (!(ukk % 255)) {
	  if (!(idx1 & 1)) {
	    unroll_zero_incr_8_32(remv[idx1 + 6], remv[(idx1 / 2) + 12], perm_ct4);
	  } else {
	    calc_rem_merge32_minus(perm_ct4, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
	  }
	}
      }
#endif
      ulxor &= ~((3 * ONELU) << ujj);
    }
#ifdef __LP64__
    permsv = &(permsv[BITCT2 * perm_ct128]);
#else
    permsv = &(permsv[BITCT2 * perm_ct32]);
#endif
  }
  for (idx1 = 0; idx1 < 6; idx1++) {
    uii = cur_cts[idx1];
#ifdef __LP64__
    if (uii % 15) {
      // todo: check if zeroing needed
      unroll_zero_incr_4_8(remv[idx1], remv[idx1 + 6], perm_ct32);
    }
    if (uii % 255) {
      if (!(idx1 & 1)) {
	unroll_zero_incr_8_32(remv[idx1 + 6], remv[(idx1 / 2) + 12], perm_ct16);
      } else {
	calc_rem_merge32_minus(perm_ct16, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
      }
    }
#else
    if (uii % 15) {
      unroll_zero_incr_4_8(remv[idx1], remv[idx1 + 6], perm_ct8);
    }
    if (uii % 255) {
      if (!(idx1 & 1)) {
	unroll_zero_incr_8_32(remv[idx1 + 6], remv[(idx1 / 2) + 12], perm_ct4);
      } else {
	calc_rem_merge32_minus(perm_ct4, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
      }
    }
#endif
  }
}

void calc_qrem(uint32_t pheno_nm_ct, uintptr_t perm_vec_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_ref, double* perm_vecstd, double* outbufs) {
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  uint32_t pheno_nm_ctl2x = QUATERCT_TO_WORDCT(pheno_nm_ct);
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the sample idx
  uint32_t row_mult = perm_vec_ctcl8m / 4;

  uint32_t loop_len = (perm_vec_ct + 1) / 2;
  __m128d* permsv = (__m128d*)perm_vecstd;
  __m128d* __restrict__ perm_readv;
  __m128d* __restrict__ rem_writev;
  __m128d* __restrict__ rem_write2v;
  __m128d* __restrict__ rem_write3v;
  __m128d vxx;
#else
  uint32_t row_mult = perm_vec_ctcl8m / 2;
  double* __restrict__ perm_read;
  double* __restrict__ rem_write;
  double* __restrict__ rem_write2;
  double* __restrict__ rem_write3;
  double dxx;
#endif
  uintptr_t ulraw1;
  uintptr_t ulxor;
  uint32_t cur_xor;
  uint32_t cur_raw;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulraw1 = *loadbuf++;
    ulxor = ulraw1 ^ (*loadbuf_ref++);
    if (uii + 1 == pheno_nm_ctl2x) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulxor &= (ONELU << (ujj * 2)) - ONELU;
      }
    }
    while (ulxor) {
      ujj = CTZLU(ulxor) & (BITCT - 2);
      cur_xor = (ulxor >> ujj) & 3;
      cur_raw = (ulraw1 >> ujj) & 3;
#ifdef __LP64__
      perm_readv = &(permsv[ujj * row_mult]);
      rem_writev = (__m128d*)outbufs;
      rem_write2v = (__m128d*)(&(outbufs[perm_vec_ctcl8m]));
      rem_write3v = (__m128d*)(&(outbufs[2 * perm_vec_ctcl8m]));
      if (cur_raw == 3) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, vxx);
	    rem_writev++;
	  }
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, _mm_add_pd(vxx, vxx));
	    rem_writev++;
	  }
        } else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write2v = _mm_sub_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write3v = _mm_sub_pd(*rem_write3v, _mm_mul_pd(vxx, vxx));
	    rem_write3v++;
	  }
        }
      } else if (cur_raw == 2) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, vxx);
	    rem_writev++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, vxx);
	    rem_writev++;
	  }
	} else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write2v = _mm_sub_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write3v = _mm_sub_pd(*rem_write3v, _mm_mul_pd(vxx, vxx));
	    rem_write3v++;
	  }
	}
      } else if (!cur_raw) {
	if (cur_xor == 3) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, _mm_add_pd(vxx, vxx));
	    rem_writev++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, vxx);
	    rem_writev++;
	  }
	} else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, _mm_add_pd(vxx, vxx));
	    rem_writev++;
	    *rem_write2v = _mm_sub_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write3v = _mm_sub_pd(*rem_write3v, _mm_mul_pd(vxx, vxx));
	    rem_write3v++;
	  }
	}
      } else {
	if (cur_xor == 2) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write2v = _mm_add_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write3v = _mm_add_pd(*rem_write3v, _mm_mul_pd(vxx, vxx));
	    rem_write3v++;
	  }
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write2v = _mm_add_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write3v = _mm_add_pd(*rem_write3v, _mm_mul_pd(vxx, vxx));
	    rem_write3v++;
	  }
	} else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, _mm_add_pd(vxx, vxx));
	    rem_writev++;
	    *rem_write2v = _mm_add_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write3v = _mm_add_pd(*rem_write3v, _mm_mul_pd(vxx, vxx));
	    rem_write3v++;
	  }
	}
      }
#else
      perm_read = &(perm_vecstd[ujj * row_mult]);
      rem_write = outbufs;
      rem_write2 = &(outbufs[perm_vec_ctcl8m]);
      rem_write3 = &(outbufs[2 * perm_vec_ctcl8m]);
      if (cur_raw == 3) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= dxx;
	    rem_write++;
	  }
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= 2 * dxx;
	    rem_write++;
	  }
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write2 -= dxx;
	    rem_write2++;
	    *rem_write3 -= dxx * dxx;
	    rem_write3++;
	  }
	}
      } else if (cur_raw == 2) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += dxx;
	    rem_write++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= dxx;
	    rem_write++;
	  }
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += dxx;
	    rem_write++;
	    *rem_write2 -= dxx;
	    rem_write2++;
	    *rem_write3 -= dxx * dxx;
	    rem_write3++;
	  }
	}
      } else if (!cur_raw) {
	if (cur_xor == 3) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += 2 * dxx;
	    rem_write++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += dxx;
	    rem_write++;
	  }
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += 2 * dxx;
	    rem_write++;
	    *rem_write2 -= dxx;
	    rem_write2++;
	    *rem_write3 -= dxx * dxx;
	    rem_write3++;
	  }
	}
      } else {
	if (cur_xor == 2) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write2 += dxx;
	    rem_write2++;
	    *rem_write3 += dxx * dxx;
	    rem_write3++;
	  }
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= dxx;
	    rem_write++;
	    *rem_write2 += dxx;
	    rem_write2++;
	    *rem_write3 += dxx * dxx;
	    rem_write3++;
	  }
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= 2 * dxx;
	    rem_write++;
	    *rem_write2 += dxx;
	    rem_write2++;
	    *rem_write3 += dxx * dxx;
	    rem_write3++;
	  }
	}
      }
#endif
      ulxor &= ~((3 * ONELU) << ujj);
    }
#ifdef __LP64__
    permsv = &(permsv[(BITCT2 / 2) * perm_vec_ctcl8m]);
#else
    perm_vecstd = &(perm_vecstd[BITCT2 * perm_vec_ctcl8m]);
#endif
  }
}

void calc_qrem_lin(uint32_t pheno_nm_ct, uintptr_t perm_vec_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_ref, double* perm_vecstd, double* outbufs) {
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  uint32_t pheno_nm_ctl2x = QUATERCT_TO_WORDCT(pheno_nm_ct);
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the sample idx
  uint32_t row_mult = perm_vec_ctcl8m / 4;

  uint32_t loop_len = (perm_vec_ct + 1) / 2;
  __m128d* permsv = (__m128d*)perm_vecstd;
  __m128d* __restrict__ perm_readv;
  __m128d* __restrict__ rem_writev = (__m128d*)outbufs;
  __m128d* __restrict__ rem_write2v = (__m128d*)(&(outbufs[perm_vec_ctcl8m]));
  __m128d* __restrict__ rem_write3v = (__m128d*)(&(outbufs[2 * perm_vec_ctcl8m]));
  __m128d* __restrict__ rem_write4v = (__m128d*)(&(outbufs[3 * perm_vec_ctcl8m]));
  __m128d* __restrict__ rem_write5v = (__m128d*)(&(outbufs[4 * perm_vec_ctcl8m]));
  __m128d* __restrict__ rem_write6v = (__m128d*)(&(outbufs[5 * perm_vec_ctcl8m]));
  __m128d vxx;
#else
  uint32_t row_mult = perm_vec_ctcl8m / 2;
  double* __restrict__ perm_read;
  double* __restrict__ rem_write = outbufs;
  double* __restrict__ rem_write2 = &(outbufs[perm_vec_ctcl8m]);
  double* __restrict__ rem_write3 = &(outbufs[2 * perm_vec_ctcl8m]);
  double* __restrict__ rem_write4 = &(outbufs[3 * perm_vec_ctcl8m]);
  double* __restrict__ rem_write5 = &(outbufs[4 * perm_vec_ctcl8m]);
  double* __restrict__ rem_write6 = &(outbufs[5 * perm_vec_ctcl8m]);
  double dxx;
#endif
  uintptr_t ulraw1;
  uintptr_t ulxor;
  uint32_t cur_xor;
  uint32_t cur_raw;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulraw1 = *loadbuf++;
    ulxor = ulraw1 ^ (*loadbuf_ref++);
    if (uii + 1 == pheno_nm_ctl2x) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulxor &= (ONELU << (ujj * 2)) - ONELU;
      }
    }
    while (ulxor) {
      ujj = CTZLU(ulxor) & (BITCT - 2);
      cur_xor = (ulxor >> ujj) & 3;
      cur_raw = (ulraw1 >> ujj) & 3;
#ifdef __LP64__
      perm_readv = &(permsv[ujj * row_mult]);
      if (cur_raw == 3) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write2v = _mm_sub_pd(*rem_write2v, _mm_mul_pd(vxx, vxx));
	    rem_write2v++;
	  }
	  rem_writev = (__m128d*)outbufs;
	  rem_write2v = (__m128d*)(&(outbufs[perm_vec_ctcl8m]));
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write3v = _mm_sub_pd(*rem_write3v, vxx);
	    rem_write3v++;
	    *rem_write4v = _mm_sub_pd(*rem_write4v, _mm_mul_pd(vxx, vxx));
	    rem_write4v++;
	  }
	  rem_write3v = (__m128d*)(&(outbufs[2 * perm_vec_ctcl8m]));
	  rem_write4v = (__m128d*)(&(outbufs[3 * perm_vec_ctcl8m]));
        } else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write5v = _mm_sub_pd(*rem_write5v, vxx);
	    rem_write5v++;
	    *rem_write6v = _mm_sub_pd(*rem_write6v, _mm_mul_pd(vxx, vxx));
	    rem_write6v++;
	  }
	  rem_write5v = (__m128d*)(&(outbufs[4 * perm_vec_ctcl8m]));
	  rem_write6v = (__m128d*)(&(outbufs[5 * perm_vec_ctcl8m]));
        }
      } else if (cur_raw == 2) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write2v = _mm_add_pd(*rem_write2v, _mm_mul_pd(vxx, vxx));
	    rem_write2v++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write3v = _mm_sub_pd(*rem_write3v, vxx);
	    rem_write3v++;
	    vxx = _mm_mul_pd(vxx, vxx);
	    *rem_write2v = _mm_add_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write4v = _mm_sub_pd(*rem_write4v, vxx);
	    rem_write4v++;
	  }
	  rem_write3v = (__m128d*)(&(outbufs[2 * perm_vec_ctcl8m]));
	  rem_write4v = (__m128d*)(&(outbufs[3 * perm_vec_ctcl8m]));
	} else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_add_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write5v = _mm_sub_pd(*rem_write5v, vxx);
	    rem_write5v++;
	    vxx = _mm_mul_pd(vxx, vxx);
	    *rem_write2v = _mm_add_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write6v = _mm_sub_pd(*rem_write6v, vxx);
	    rem_write6v++;
	  }
	  rem_write5v = (__m128d*)(&(outbufs[4 * perm_vec_ctcl8m]));
	  rem_write6v = (__m128d*)(&(outbufs[5 * perm_vec_ctcl8m]));
	}
	rem_writev = (__m128d*)outbufs;
	rem_write2v = (__m128d*)(&(outbufs[perm_vec_ctcl8m]));
      } else if (!cur_raw) {
	if (cur_xor == 3) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write3v = _mm_add_pd(*rem_write3v, vxx);
	    rem_write3v++;
	    *rem_write4v = _mm_add_pd(*rem_write4v, _mm_mul_pd(vxx, vxx));
	    rem_write4v++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write3v = _mm_add_pd(*rem_write3v, vxx);
	    rem_write3v++;
	    vxx = _mm_mul_pd(vxx, vxx);
	    *rem_write2v = _mm_sub_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write4v = _mm_add_pd(*rem_write4v, vxx);
	    rem_write4v++;
	  }
	  rem_writev = (__m128d*)outbufs;
	  rem_write2v = (__m128d*)(&(outbufs[perm_vec_ctcl8m]));
	} else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write3v = _mm_add_pd(*rem_write3v, vxx);
	    rem_write3v++;
	    *rem_write5v = _mm_sub_pd(*rem_write5v, vxx);
	    rem_write5v++;
	    vxx = _mm_mul_pd(vxx, vxx);
	    *rem_write4v = _mm_add_pd(*rem_write4v, vxx);
	    rem_write4v++;
	    *rem_write6v = _mm_sub_pd(*rem_write6v, vxx);
	    rem_write6v++;
	  }
	  rem_write5v = (__m128d*)(&(outbufs[4 * perm_vec_ctcl8m]));
	  rem_write6v = (__m128d*)(&(outbufs[5 * perm_vec_ctcl8m]));
	}
	rem_write3v = (__m128d*)(&(outbufs[2 * perm_vec_ctcl8m]));
	rem_write4v = (__m128d*)(&(outbufs[3 * perm_vec_ctcl8m]));
      } else {
	if (cur_xor == 2) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write5v = _mm_add_pd(*rem_write5v, vxx);
	    rem_write5v++;
	    *rem_write6v = _mm_add_pd(*rem_write6v, _mm_mul_pd(vxx, vxx));
	    rem_write6v++;
	  }
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_writev = _mm_sub_pd(*rem_writev, vxx);
	    rem_writev++;
	    *rem_write5v = _mm_add_pd(*rem_write5v, vxx);
	    rem_write5v++;
	    vxx = _mm_mul_pd(vxx, vxx);
	    *rem_write2v = _mm_sub_pd(*rem_write2v, vxx);
	    rem_write2v++;
	    *rem_write6v = _mm_add_pd(*rem_write6v, vxx);
	    rem_write6v++;
	  }
	  rem_writev = (__m128d*)outbufs;
	  rem_write2v = (__m128d*)(&(outbufs[perm_vec_ctcl8m]));
	} else {
	  for (ukk = 0; ukk < loop_len; ukk++) {
	    vxx = *perm_readv++;
	    *rem_write3v = _mm_sub_pd(*rem_write3v, vxx);
	    rem_write3v++;
	    *rem_write5v = _mm_add_pd(*rem_write5v, vxx);
	    rem_write5v++;
	    vxx = _mm_mul_pd(vxx, vxx);
	    *rem_write4v = _mm_sub_pd(*rem_write4v, vxx);
	    rem_write4v++;
	    *rem_write6v = _mm_add_pd(*rem_write6v, vxx);
	    rem_write6v++;
	  }
	  rem_write3v = (__m128d*)(&(outbufs[2 * perm_vec_ctcl8m]));
	  rem_write4v = (__m128d*)(&(outbufs[3 * perm_vec_ctcl8m]));
	}
	rem_write5v = (__m128d*)(&(outbufs[4 * perm_vec_ctcl8m]));
	rem_write6v = (__m128d*)(&(outbufs[5 * perm_vec_ctcl8m]));
      }
#else
      perm_read = &(perm_vecstd[ujj * row_mult]);
      if (cur_raw == 3) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= dxx;
	    rem_write++;
	    *rem_write2 -= dxx * dxx;
	    rem_write2++;
	  }
          rem_write = outbufs;
          rem_write2 = &(outbufs[perm_vec_ctcl8m]);
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write3 -= dxx;
	    rem_write3++;
	    *rem_write4 -= dxx * dxx;
	    rem_write4++;
	  }
          rem_write3 = &(outbufs[2 * perm_vec_ctcl8m]);
          rem_write4 = &(outbufs[3 * perm_vec_ctcl8m]);
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write5 -= dxx;
	    rem_write5++;
	    *rem_write6 -= dxx * dxx;
	    rem_write6++;
	  }
          rem_write5 = &(outbufs[4 * perm_vec_ctcl8m]);
          rem_write6 = &(outbufs[5 * perm_vec_ctcl8m]);
	}
      } else if (cur_raw == 2) {
	if (cur_xor == 1) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += dxx;
	    rem_write++;
	    *rem_write2 += dxx * dxx;
	    rem_write2++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += dxx;
	    rem_write++;
	    *rem_write3 -= dxx;
	    rem_write3++;
	    dxx *= dxx;
	    *rem_write2 += dxx;
	    rem_write2++;
	    *rem_write4 -= dxx;
	    rem_write4++;
	  }
          rem_write3 = &(outbufs[2 * perm_vec_ctcl8m]);
          rem_write4 = &(outbufs[3 * perm_vec_ctcl8m]);
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write += dxx;
	    rem_write++;
	    *rem_write5 -= dxx;
	    rem_write5++;
	    dxx *= dxx;
	    *rem_write2 += dxx;
	    rem_write2++;
	    *rem_write6 -= dxx;
	    rem_write6++;
	  }
          rem_write5 = &(outbufs[4 * perm_vec_ctcl8m]);
          rem_write6 = &(outbufs[5 * perm_vec_ctcl8m]);
	}
	rem_write = outbufs;
	rem_write2 = &(outbufs[perm_vec_ctcl8m]);
      } else if (!cur_raw) {
	if (cur_xor == 3) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write3 += dxx;
	    rem_write3++;
	    *rem_write4 += dxx * dxx;
	    rem_write4++;
	  }
	} else if (cur_xor == 2) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= dxx;
	    rem_write++;
	    *rem_write3 += dxx;
	    rem_write3++;
	    dxx *= dxx;
	    *rem_write2 -= dxx;
	    rem_write2++;
	    *rem_write4 += dxx;
	    rem_write4++;
	  }
	  rem_write = outbufs;
	  rem_write2 = &(outbufs[perm_vec_ctcl8m]);
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write3 += dxx;
	    rem_write3++;
	    *rem_write5 -= dxx;
	    rem_write5++;
	    dxx *= dxx;
	    *rem_write4 += dxx;
	    rem_write4++;
	    *rem_write6 -= dxx;
	    rem_write6++;
	  }
          rem_write5 = &(outbufs[4 * perm_vec_ctcl8m]);
          rem_write6 = &(outbufs[5 * perm_vec_ctcl8m]);
	}
	rem_write3 = &(outbufs[2 * perm_vec_ctcl8m]);
	rem_write4 = &(outbufs[3 * perm_vec_ctcl8m]);
      } else {
	if (cur_xor == 2) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write5 += dxx;
	    rem_write5++;
	    *rem_write6 += dxx * dxx;
	    rem_write6++;
	  }
	} else if (cur_xor == 3) {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write -= dxx;
	    rem_write++;
	    *rem_write5 += dxx;
	    rem_write5++;
	    dxx *= dxx;
	    *rem_write2 -= dxx;
	    rem_write2++;
	    *rem_write6 += dxx;
	    rem_write6++;
	  }
	  rem_write = outbufs;
	  rem_write2 = &(outbufs[perm_vec_ctcl8m]);
	} else {
	  for (ukk = 0; ukk < perm_vec_ct; ukk++) {
	    dxx = *perm_read++;
	    *rem_write3 -= dxx;
	    rem_write3++;
	    *rem_write5 += dxx;
	    rem_write5++;
	    dxx *= dxx;
	    *rem_write4 -= dxx;
	    rem_write4++;
	    *rem_write6 += dxx;
	    rem_write6++;
	  }
	  rem_write3 = &(outbufs[2 * perm_vec_ctcl8m]);
	  rem_write4 = &(outbufs[3 * perm_vec_ctcl8m]);
	}
	rem_write5 = &(outbufs[4 * perm_vec_ctcl8m]);
	rem_write6 = &(outbufs[5 * perm_vec_ctcl8m]);
      }
#endif
      ulxor &= ~((3 * ONELU) << ujj);
    }
#ifdef __LP64__
    permsv = &(permsv[(BITCT2 / 2) * perm_vec_ctcl8m]);
#else
    perm_vecstd = &(perm_vecstd[BITCT2 * perm_vec_ctcl8m]);
#endif
  }
}

void check_for_better_rem_cost(uintptr_t best_cost, uint32_t maxt_block_base, uint32_t maxt_block_base2, uint32_t maxt_block_base3, uintptr_t marker_idx, uint32_t* __restrict__ missing_cts, uint32_t* __restrict__ homcom_cts, uint32_t* __restrict__ het_cts, uint16_t* ldrefs, uint32_t pheno_nm_ct, int32_t missing_ct, int32_t het_ct, int32_t homcom_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_cur, uint32_t* ldrefp) {
  // Check if PERMORY-style LD exploitation is better than genotype indexing
  // algorithm.
  //
  // Effective inner loop iterations required for LD exploitation:
  //   2 * (<-> neither side homcom) + (<-> homcom) + constant
  // Simple lower bound:
  //   max(delta(homcom), delta(non-homcom)) + constant
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t marker_idx_tmp = maxt_block_base;
  uint32_t loop_ceil = maxt_block_base2;
  int32_t homrar_ct = pheno_nm_ct - missing_ct - het_ct - homcom_ct;
  int32_t missing_ct_tmp;
  int32_t het_ct_tmp;
  int32_t homcom_ct_tmp;
  int32_t homrar_ct_tmp;
  uint32_t marker_bidx2;
  uintptr_t homcom_delta;
  uintptr_t cur_cost;
  do {
    if (marker_idx_tmp == maxt_block_base2) {
      marker_idx_tmp = maxt_block_base3;
      loop_ceil = marker_idx;
    }
    for (; marker_idx_tmp < loop_ceil; marker_idx_tmp++) {
      if (ldrefs[marker_idx_tmp] != 65535) {
	missing_ct_tmp = missing_cts[marker_idx_tmp];
	homcom_ct_tmp = homcom_cts[marker_idx_tmp];
	het_ct_tmp = het_cts[marker_idx_tmp];
	homrar_ct_tmp = pheno_nm_ct - missing_ct_tmp - het_ct_tmp - homcom_ct_tmp;
	homcom_delta = labs(((int32_t)homcom_ct) - homcom_ct_tmp);
	cur_cost = labs(((int32_t)missing_ct) - missing_ct_tmp) + labs(((int32_t)homrar_ct) - homrar_ct_tmp) + labs(((int32_t)het_ct) - het_ct_tmp);
	cur_cost = MAXV(homcom_delta, cur_cost);
	if (cur_cost < best_cost) {
	  marker_bidx2 = marker_idx_tmp - maxt_block_base;
	  cur_cost = rem_cost(pheno_nm_ctv2, &(loadbuf[marker_bidx2 * pheno_nm_ctv2]), loadbuf_cur);
	  if (cur_cost < best_cost) {
	    *ldrefp = marker_bidx2;
	    best_cost = cur_cost;
	  }
	}
      }
    }
  } while (marker_idx_tmp < marker_idx);
}

// multithread globals
static double* g_orig_pvals;
static double* g_orig_chisq;
static double* g_mperm_save_all;

// A separated-low-and-high-bit format was tried, and found to not really be
// any better than the usual PLINK 2-bit format.
static uintptr_t* g_loadbuf;

static uint32_t* g_perm_vecst; // genotype indexing support
static uint32_t* g_thread_git_wkspace;
static uint32_t* g_resultbuf;

// always use genotype indexing for QT --assoc
static double* g_thread_git_qbufs;
static double* g_qresultbuf;
static double g_pheno_sum;
static double g_pheno_ssq;
static uint16_t* g_ldrefs;
static double* g_orig_linsq; // square of Lin t-statistic

// maximum number of precomputed table entries per marker
static uint32_t g_precomp_width;
// precomputed table contains entries for missing_cts ranging from
//   g_precomp_start[marker_bidx] to
//   (g_precomp_start[marker_bidx] + g_precomp_width - 1).
static uint32_t g_precomp_start[MODEL_BLOCKSIZE];

// Space for precomputed tables to accelerate permutation p-value computations.
// The sizing and usage of this space varies depending on the permutation
// analysis requested.  (The main objective here is to bring Fisher 2x2 exact
// p-values to the masses.  There's a very minor chi-square speedup as well;
// it's really only present because it allowed for simpler debugging of parts
// of the Fisher logic.)
//
// In what follows,
//   n := (g_precomp_width * marker_bidx) + missing_ct -
//        g_precomp_start[marker_bidx].
//
// For --assoc perm/--model {dom|rec|trend} perm:
//   g_precomp_ui[4n] and [4n + 1] define the interval with less extreme
//     p-values than the original.  [4n + 2] and [4n + 3] define the
//     interval with less or equally extreme p-values.
//
// For --assoc mperm fisher/--model {dom|rec} fisher:
//   g_precomp_ui[6n]...[6n + 3] is as in --assoc perm.
//   g_precomp_ui[6n + 4] and [6n + 5] are the floor and offset for the
//     range of case_set_cts where Fisher p-value calculation is unnecessary.
//   g_precomp_d[2n] and [2n + 1] are tot_prob and right_prob for
//     fisher22_tail_pval().  (This is almost irrelevant.)
//
// For --assoc mperm/--model {dom|rec|trend} mperm:
//   g_precomp_ui is as in --assoc mperm fisher.
//   g_precomp_d[2n] and [2n + 1] are expm11 and recip_sum from
//     chi22_get_coeffs()/ca_trend_get_coeffs().
//
// For --model perm-gen:
//   No precomputation at all.
//
// For regular --model perm:
//   g_precomp_ui[12n] to [12n + 3] cover the allelic test, [12n + 4] to
//     [12n + 7] cover the dom test, and [12n + 8] to [12n + 11] cover rec.
//     [12n + 4] is 0xffffffff if the dom and rec tests should be skipped.
//
// For regular --model mperm fisher:
//   g_precomp_ui[18n] to [18n + 5] cover the allelic test, etc.
//   g_precomp_d[6n] to [6n + 1] are fisher22_tail_pval() coefficients for the
//     allelic test, etc.
//
// For regular --model mperm:
//   g_precomp_ui as in --model mperm fisher.
//   g_precomp_d[6n] and [6n + 1] are expm11 and recip_sum for the allelic
//     test, etc.
//
static uint32_t* g_precomp_ui;
static double* g_precomp_d;

// X-chromosome: number of missing allele observations per marker relative to
//   *all female* case (so all males automatically contribute at least 1)
// elsewhere: number of missing samples for each marker
static uint32_t* g_missing_cts;

static uint32_t* g_set_cts;
static uint32_t* g_het_cts;
static uint32_t* g_homcom_cts;

// This is *twice* the number of successes, because PLINK 1.07 counts tie as
// 0.5.  (Actually, it randomizes instead of deterministically adding 0.5; this
// randomization just adds noise so we don't replicate it.)
static uint32_t* g_perm_2success_ct;
static uint32_t* g_perm_attempt_ct;
static double* g_maxt_extreme_stat;
static double* g_maxt_thread_results;

// to avoid pathological multithreading issues, this is not a bitset
static unsigned char* g_perm_adapt_stop;

static uint32_t g_adapt_m_table[MODEL_BLOCKSIZE];
static uintptr_t* g_sample_nonmale_include2;
static uintptr_t* g_sample_male_include2;
static uintptr_t* g_is_invalid_bitfield;
static uint32_t g_model_fisher;
static uint32_t g_fisher_midp;
static uint32_t g_assoc_thread_ct;
static uint32_t g_maxt_block_base;
static uint32_t g_block_start;
static uint32_t g_qblock_start;
static uint32_t g_block_diff;
static uint32_t g_perms_done;
static uint32_t g_first_adapt_check;
static uint32_t g_male_ct;
static double g_adaptive_intercept;
static double g_adaptive_slope;
static double g_aperm_alpha;
static double g_adaptive_ci_zt;
static uint32_t g_is_x;
static uint32_t g_is_y;

// X, Y, MT.  note that X, and now MT as well, have max ploidy 2
static uint32_t g_min_ploidy_1;

static int32_t g_is_model_prec;

static uint32_t* g_male_case_cts;

THREAD_RET_TYPE assoc_adapt_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_perm_case_ct;
  uintptr_t* __restrict__ male_vec = g_sample_male_include2;
  uintptr_t* __restrict__ nonmale_vec = g_sample_nonmale_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t* __restrict__ loadbuf;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ set_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t min_ploidy_1;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  uint32_t min_ploidy;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t missing_start;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  while (1) {
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto assoc_adapt_thread_skip_all;
      }
      marker_bidx = g_block_start + tidx;
      marker_bceil = marker_bidx + 1;
    } else {
      marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    min_ploidy_1 = g_min_ploidy_1;
    loadbuf = g_loadbuf;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    missing_cts = g_missing_cts;
    set_cts = g_set_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    is_x = g_is_x;
    is_y = g_is_y;
    if (min_ploidy_1) {
      min_ploidy = 1;
    } else {
      min_ploidy = 2;
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      // guaranteed during loading that g_perm_adapt_stop[] is not set yet
      marker_idx = g_adapt_m_table[marker_bidx];
      next_adapt_check = first_adapt_check;
      col1_sum = set_cts[marker_idx];
      if (is_x) {
	row1x_sum = 2 * case_ct;
	tot_obs = 2 * pheno_nm_ct - missing_cts[marker_idx];
      } else {
	row1x_sum = min_ploidy * case_ct;
	tot_obs = min_ploidy * (pheno_nm_ct - missing_cts[marker_idx]);
      }
      col2_sum = tot_obs - col1_sum;
      missing_start = precomp_start[marker_bidx];
      gpui = &(precomp_ui[4 * precomp_width * marker_bidx]);
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      if (orig_pvals[marker_idx] == -9) {
        perm_adapt_stop[marker_idx] = 1;
        perm_attempt_ct[marker_idx] = next_adapt_check;
        perm_2success_ct[marker_idx] = next_adapt_check;
        continue;
      }
      if (model_fisher) {
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      for (pidx = 0; pidx < perm_vec_ct;) {
	if (!min_ploidy_1) {
	  genovec_set_freq(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	} else if (is_x) {
	  genovec_set_freq_x(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), male_vec, pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	} else if (!is_y) {
	  genovec_3freq(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), pheno_nm_ctv2, &case_missing_ct, &uii, &case_set_ct);
	  case_missing_ct += uii;
	} else {
	  genovec_set_freq_y(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), nonmale_vec, pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	}
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_set_ct < gpui[4 * uii]) {
	    if (case_set_ct < gpui[4 * uii + 2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (case_set_ct >= gpui[4 * uii + 1]) {
	      if (case_set_ct >= gpui[4 * uii + 3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	} else {
	  uii = row1x_sum - case_missing_ct * min_ploidy; // row1_sum
	  if (model_fisher) {
	    dxx = fisher22(case_set_ct, uii - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - uii, fisher_midp);
	    if (dxx < stat_low) {
	      success_2incr += 2;
	    } else if (dxx <= stat_high) {
	      success_2incr++;
	    }
	  } else {
	    dxx = chi22_eval(case_set_ct, uii, col1_sum, tot_obs);
	    if (dxx > stat_high) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  }
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    dyy = pval - dxx; // lower bound
	    dzz = pval + dxx; // upper bound
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  assoc_adapt_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE assoc_maxt_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;

  // currently safe for this to be uint32_t since perm_vec_ct < 2^30
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_perm_case_ct;
  uintptr_t* __restrict__ male_vec = g_sample_male_include2;
  uintptr_t* __restrict__ nonmale_vec = g_sample_nonmale_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ mperm_save_all = g_mperm_save_all;
  double* msa_ptr = nullptr;
  uintptr_t* __restrict__ loadbuf;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ set_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  double* __restrict__ precomp_d;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint16_t* ldrefs;
  uintptr_t* loadbuf_cur;
  uint32_t* resultbuf;
  uint32_t* gpui;
  double* gpd;
  uintptr_t pidx;
  uintptr_t marker_idx;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t block_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t is_x;
  uint32_t is_x_or_y;
  uint32_t min_ploidy_1;
  uint32_t min_ploidy;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double stat_high;
  double stat_low;
  double sval;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t ldref;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto assoc_maxt_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + block_start;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    is_x = g_is_x;
    is_x_or_y = is_x || g_is_y;
    min_ploidy_1 = g_min_ploidy_1;
    memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
    if (min_ploidy_1) {
      min_ploidy = 1;
    } else {
      min_ploidy = 2;
    }
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    set_cts = g_set_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    precomp_d = g_precomp_d;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    resultbuf = g_resultbuf;
    ldrefs = g_ldrefs;

    if (mperm_save_all) {
      msa_ptr = &(mperm_save_all[marker_idx * perm_vec_ct]);
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      if (orig_pvals[marker_idx] == -9) {
        if (msa_ptr) {
          for (pidx = 0; pidx < perm_vec_ct; pidx++) {
            *msa_ptr++ = -9;
          }
        }
        perm_2success_ct[marker_idx++] += perm_vec_ct;
        continue;
      }
      if (model_fisher) {
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      col1_sum = set_cts[marker_idx];
      missing_ct = missing_cts[marker_idx];
      if (is_x) {
	row1x_sum = 2 * case_ct;
	tot_obs = 2 * pheno_nm_ct - missing_ct;
      } else {
	row1x_sum = min_ploidy * case_ct;
	tot_obs = min_ploidy * (pheno_nm_ct - missing_ct);
      }
      col2_sum = tot_obs - col1_sum;
      gpui = &(precomp_ui[6 * precomp_width * marker_bidx]);
      missing_start = precomp_start[marker_bidx];
      success_2incr = 0;
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      if (!is_x_or_y) {
	ldref = ldrefs[marker_idx];
	if (!min_ploidy_1) {
	  het_ct = het_cts[marker_idx];
	  homcom_ct = (col1_sum - het_ct) / 2;
	} else {
	  het_ct = 0;
	  homcom_ct = col1_sum;
	}
	git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
	git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
	git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
	if (ldref == 65535) {
	  ldref = marker_bidx;
	  if (pheno_nm_ct - homcom_ct > 50) {
	    check_for_better_rem_cost(pheno_nm_ct - homcom_ct - 50, maxt_block_base, maxt_block_base2, maxt_block_base3, marker_idx, missing_cts, homcom_cts, het_cts, ldrefs, pheno_nm_ct, missing_ct, het_ct, homcom_ct, loadbuf, loadbuf_cur, &ldref);
	  }
	  ldrefs[marker_idx] = ldref;
	}
	if (ldref == marker_bidx) {
	  fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
	  calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
	  fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
	} else {
	  memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
	  calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
	}
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	if (!is_x_or_y) {
	  if (!min_ploidy_1) {
	    case_missing_ct = git_missing_cts[pidx];
	    case_set_ct = row1x_sum - (git_het_cts[pidx] + 2 * (case_missing_ct + git_homrar_cts[pidx]));
	  } else {
	    case_missing_ct = git_missing_cts[pidx] + git_het_cts[pidx];
	    case_set_ct = row1x_sum - case_missing_ct - git_homrar_cts[pidx];
	  }
	} else {
	  if (is_x) {
	    genovec_set_freq_x(loadbuf_cur, &(perm_vecs[pidx * pheno_nm_ctv2]), male_vec, pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	  } else {
	    genovec_set_freq_y(loadbuf_cur, &(perm_vecs[pidx * pheno_nm_ctv2]), nonmale_vec, pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	  }
	}
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_set_ct < gpui[6 * uii]) {
	    if (case_set_ct < gpui[6 * uii + 2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (case_set_ct >= gpui[6 * uii + 1]) {
	      if (case_set_ct >= gpui[6 * uii + 3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	  ukk = gpui[6 * uii + 4];
	  ujj = (uint32_t)(case_set_ct - ukk); // deliberate underflow
	  if (ujj >= gpui[6 * uii + 5]) {
	    if (model_fisher) {
	      ujj = row1x_sum - case_missing_ct * min_ploidy;
	      // sval = fisher22(case_set_ct, ujj - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - ujj);
	      sval = fisher22_tail_pval(ukk, ujj - ukk, col1_sum - ukk, col2_sum + ukk - ujj, gpui[6 * uii + 5] - 1, gpd[2 * uii], gpd[2 * uii + 1], fisher_midp, case_set_ct);
	      if (results[pidx] > sval) {
		results[pidx] = sval;
	      }
	    } else {
	      sval = ((double)((intptr_t)case_set_ct)) - gpd[2 * uii];
	      sval = sval * sval * gpd[2 * uii + 1];
	      if (results[pidx] < sval) {
		results[pidx] = sval;
	      }
	    }
	  }
	} else {
	  uii = row1x_sum - case_missing_ct * min_ploidy;
	  if (model_fisher) {
	    sval = fisher22(case_set_ct, uii - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - uii, fisher_midp);
	    if (sval < stat_low) {
	      success_2incr += 2;
	    } else if (sval <= stat_high) {
	      success_2incr++;
	    }
	    if (results[pidx] > sval) {
	      results[pidx] = sval;
	    }
	  } else {
	    sval = chi22_eval(case_set_ct, uii, col1_sum, tot_obs);
	    if (sval > stat_high) {
	      success_2incr += 2;
	    } else if (sval > stat_low) {
	      success_2incr++;
	    }
	    if (results[pidx] < sval) {
	      results[pidx] = sval;
	    }
	  }
	  if (msa_ptr) {
	    *msa_ptr++ = sval;
	  }
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  assoc_maxt_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE assoc_set_thread(void* arg) {
  // Basically a simplified version of what assoc_maxt_thread() does; we save
  // chi-square stats for the given number of permutations for all still-active
  // variants.  Adaptive pruning, if applicable, happens outside this loop.
  //
  // LD-exploitation should be added if this sees significant usage.
  // (possible todo: permit Fisher test, converting p-values into equivalent
  // chi-square stats?)
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t case_ct = g_perm_case_ct;
  uintptr_t* __restrict__ male_vec = g_sample_male_include2;
  uintptr_t* __restrict__ nonmale_vec = g_sample_nonmale_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  double* msa_ptr = nullptr;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_cur;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ set_allele_cts;
  uintptr_t pidx;
  uintptr_t marker_idx;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t tot_obs;
  uint32_t block_start;
  uint32_t marker_bidx_start;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t is_x;
  uint32_t is_x_or_y;
  uint32_t min_ploidy_1;
  uint32_t min_ploidy;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t missing_ct;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
	goto assoc_set_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    marker_bidx = marker_bidx_start;
    is_x = g_is_x;
    is_x_or_y = is_x || g_is_y;
    min_ploidy_1 = g_min_ploidy_1;
    min_ploidy = 2;
    if (min_ploidy_1) {
      min_ploidy = 1;
    }
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    set_allele_cts = g_set_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      msa_ptr = &(g_mperm_save_all[marker_bidx * perm_vec_ct]);
      col1_sum = set_allele_cts[marker_idx];
      missing_ct = missing_cts[marker_idx];
      if (is_x) {
	row1x_sum = 2 * case_ct;
	tot_obs = 2 * pheno_nm_ct - missing_ct;
      } else {
	row1x_sum = min_ploidy * case_ct;
	tot_obs = min_ploidy * (pheno_nm_ct - missing_ct);
      }
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      if (!is_x_or_y) {
	git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
	git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
	git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
	fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
	calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
	fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	if (!is_x_or_y) {
	  if (!min_ploidy_1) {
	    case_missing_ct = git_missing_cts[pidx];
	    case_set_ct = row1x_sum - (git_het_cts[pidx] + 2 * (case_missing_ct + git_homrar_cts[pidx]));
	  } else {
	    case_missing_ct = git_missing_cts[pidx] + git_het_cts[pidx];
	    case_set_ct = row1x_sum - case_missing_ct - git_homrar_cts[pidx];
	  }
	} else {
	  if (is_x) {
	    genovec_set_freq_x(loadbuf_cur, &(perm_vecs[pidx * pheno_nm_ctv2]), male_vec, pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	  } else {
	    genovec_set_freq_y(loadbuf_cur, &(perm_vecs[pidx * pheno_nm_ctv2]), nonmale_vec, pheno_nm_ctv2, &case_set_ct, &case_missing_ct);
	  }
	}
	// Fisher's exact test not supported since we are adding raw chi-square
	// stats, so little to gain from precomputation
	*msa_ptr++ = chi22_eval(case_set_ct, row1x_sum - case_missing_ct * min_ploidy, col1_sum, tot_obs);
      }
    }
  assoc_set_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE qassoc_adapt_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t max_thread_ct = g_assoc_thread_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* git_qt_g_prod = &(g_thread_git_qbufs[perm_vec_ctcl8m * tidx * 3]);
  double* git_qt_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 3 + 1)]);
  double* git_qt_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 3 + 2)]);
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  unsigned char* perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t* __restrict__ loadbuf;
  double* __restrict__ orig_chiabs;
  uintptr_t next_cqg;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uintptr_t ulii;
  intptr_t geno_sum;
  intptr_t geno_ssq;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  uint32_t nanal;
  uint32_t cur_thread_ct;
  uint32_t next_adapt_check;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t uii;
  double nanal_recip;
  double nanal_m1_recip;
  double geno_mean;
  double geno_var;
  double qt_sum;
  double qt_ssq;
  double qt_g_prod;
  double qt_mean;
  double qt_var;
  double qt_g_covar;
  double beta;
  double betasq;
  double dxx;
  double dyy;
  double dzz;
  double stat_high;
  double stat_low;
  double sval;
  while (1) {
    cur_thread_ct = g_block_diff / CACHELINE_DBL;
    if (cur_thread_ct > max_thread_ct) {
      cur_thread_ct = max_thread_ct;
    } else if (!cur_thread_ct) {
      cur_thread_ct = 1;
    }
    if (cur_thread_ct <= tidx) {
      goto qassoc_adapt_thread_skip_all;
    }
    marker_bidx = g_qblock_start + (((uint64_t)tidx) * g_block_diff) / cur_thread_ct;
    marker_bceil = g_qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / cur_thread_ct;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    orig_chiabs = g_orig_chisq;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      next_adapt_check = first_adapt_check;
      missing_ct = missing_cts[marker_idx];
      nanal = pheno_nm_ct - missing_ct;
      homcom_ct = homcom_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      homrar_ct = nanal - het_ct - homcom_ct;
      if ((nanal < 3) || (homcom_ct == nanal) || (het_ct == nanal) || (homrar_ct == nanal)) {
	// the current code might otherwise report a spurious association if
	// geno_var is zero, so we explicitly check for it here.
	perm_adapt_stop[marker_idx] = 1;
	perm_attempt_ct[marker_idx] = 0;
	continue;
      }
      sval = orig_chiabs[marker_idx];
      // tstat = beta / vbeta_sqrt
      // tstat^2 = beta * beta / vbeta;
      //         = beta^2 * (nanal - 2) / ((qt_var / geno_var) - beta^2)
      // [stop here for max(T) since nanal varies across markers]
      // tstat^2 / (nanal - 2) = beta^2 / ((qt_var / geno_var) - beta^2)
      //                       = beta^2 * geno_var / (qt_var - beta^2 * geno_var)
      // Larger values of this last statistic monotonically result in smaller
      // P-values, so this is what we use for comparison (this saves a few
      // floating point operations at the end).
      sval = sval * sval / ((double)(((int32_t)nanal) - 2));
      stat_high = sval + EPSILON;
      stat_low = sval - EPSILON;
      geno_sum = 2 * homrar_ct + het_ct;
      geno_ssq = 4 * homrar_ct + het_ct;
      nanal_recip = 1.0 / ((double)((int32_t)nanal));
      nanal_m1_recip = 1.0 / ((double)(((int32_t)nanal) - 1));
      geno_mean = ((double)geno_sum) * nanal_recip;
      geno_var = (((double)geno_ssq) - geno_sum * geno_mean) * nanal_m1_recip;
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      next_cqg = 0;
      for (pidx = 0; pidx < perm_vec_ct;) {
	if (pidx == next_cqg) {
	  next_cqg = next_adapt_check;
	  ulii = pidx + pidx_offset;
	  if (next_cqg < ulii + (ulii >> 2)) {
	    // increase ~25% at a time
	    next_cqg = ulii + (ulii >> 2);
	  }
	  next_cqg -= pidx_offset;
	  next_cqg = round_up_pow2(next_cqg, CACHELINE_DBL);
	  if (next_cqg > perm_vec_ct) {
	    next_cqg = perm_vec_ct;
	  }
	  calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, next_cqg - pidx, &(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecstd[pidx]), &(git_qt_g_prod[pidx]));
	}
	qt_sum = pheno_sum - git_qt_sum[pidx];
	qt_ssq = pheno_ssq - git_qt_ssq[pidx];
	qt_g_prod = git_qt_g_prod[pidx];
	qt_mean = qt_sum * nanal_recip;
	qt_var = (qt_ssq - qt_sum * qt_mean) * nanal_m1_recip;
	qt_g_covar = (qt_g_prod - qt_sum * geno_mean) * nanal_m1_recip;
	dxx = 1.0 / geno_var;
	beta = qt_g_covar * dxx;
	betasq = beta * beta;
	sval = betasq / (qt_var * dxx - betasq);
	if (sval > stat_high) {
	  success_2incr += 2;
	} else if (sval > stat_low) {
	  success_2incr++;
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    sval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(sval * (1 - sval) / ((int32_t)next_adapt_check));
	    dyy = sval - dxx; // lower bound
	    dzz = sval + dxx; // upper bound
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      fill_double_zero(next_cqg, git_qt_g_prod);
	      fill_double_zero(next_cqg, git_qt_sum);
	      fill_double_zero(next_cqg, git_qt_ssq);
	      goto qassoc_adapt_thread_lesszero;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      fill_double_zero(perm_vec_ctcl8m * 3, git_qt_g_prod);
    qassoc_adapt_thread_lesszero:
      perm_2success_ct[marker_idx] += success_2incr;
    }
  qassoc_adapt_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE qassoc_adapt_lin_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t max_thread_ct = g_assoc_thread_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* git_qt_het_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * tidx * 6]);
  double* git_qt_het_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 1)]);
  double* git_qt_homrar_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 2)]);
  double* git_qt_homrar_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 3)]);
  double* git_qt_missing_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 4)]);
  double* git_qt_missing_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 5)]);
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  unsigned char* perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t* __restrict__ loadbuf;
  double* __restrict__ orig_linsq;
  uintptr_t next_cqg;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uintptr_t ulii;
  intptr_t geno_sum;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  uint32_t nanal;
  uint32_t cur_thread_ct;
  double het_ctd;
  double homrar_ctd;
  double nanal_recip;
  double geno_mean;
  double geno_mean_sq;
  double geno_mean_coeff2;
  double geno_mean_coeff3;
  double qt_sum;
  double qt_ssq;
  double qt_het_sum;
  double qt_het_ssq;
  double qt_homrar_sum;
  double qt_homrar_ssq;
  double qt_g_prod;
  double qt_mean;
  double qt_g_prod_centered;
  double dxx;
  double dyy;
  double dzz;
  uint32_t next_adapt_check;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double sval;
  while (1) {
    cur_thread_ct = g_block_diff / CACHELINE_DBL;
    if (cur_thread_ct > max_thread_ct) {
      cur_thread_ct = max_thread_ct;
    } else if (!cur_thread_ct) {
      cur_thread_ct = 1;
    }
    if (cur_thread_ct <= tidx) {
      goto qassoc_adapt_lin_thread_skip_all;
    }
    marker_bidx = g_qblock_start + (((uint64_t)tidx) * g_block_diff) / cur_thread_ct;
    marker_bceil = g_qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / cur_thread_ct;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    orig_linsq = g_orig_linsq;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      next_adapt_check = first_adapt_check;
      missing_ct = missing_cts[marker_idx];
      nanal = pheno_nm_ct - missing_ct;
      homcom_ct = homcom_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      if ((nanal < 3) || (homcom_ct == nanal) || (het_ct == nanal)) {
	perm_adapt_stop[marker_idx] = 1;
	perm_attempt_ct[marker_idx] = 0;
	continue;
      }
      homrar_ct = nanal - het_ct - homcom_ct;
      sval = orig_linsq[marker_idx];
      stat_high = sval + EPSILON;
      stat_low = sval - EPSILON;
      geno_sum = 2 * homrar_ct + het_ct;
      nanal_recip = 1.0 / ((double)((int32_t)nanal));
      het_ctd = het_ct;
      homrar_ctd = homrar_ct;
      geno_mean = ((double)geno_sum) * nanal_recip;
      geno_mean_sq = geno_mean * geno_mean;
      geno_mean_coeff2 = 1 - 2 * geno_mean;
      geno_mean_coeff3 = 4 - 4 * geno_mean;
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      next_cqg = 0;
      for (pidx = 0; pidx < perm_vec_ct;) {
	if (pidx == next_cqg) {
	  next_cqg = next_adapt_check;
	  ulii = pidx + pidx_offset;
	  if (next_cqg < ulii + (ulii >> 2)) {
	    // increase ~25% at a time
	    next_cqg = ulii + (ulii >> 2);
	  }
	  next_cqg -= pidx_offset;
	  next_cqg = round_up_pow2(next_cqg, CACHELINE_DBL);
	  if (next_cqg > perm_vec_ct) {
	    next_cqg = perm_vec_ct;
	  }
	  calc_qgit_lin(pheno_nm_ct, perm_vec_ctcl8m, next_cqg - pidx, &(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecstd[pidx]), &(git_qt_het_sum[pidx]));
	}
	qt_sum = pheno_sum - git_qt_missing_sum[pidx];
	qt_ssq = pheno_ssq - git_qt_missing_ssq[pidx];
	qt_het_sum = git_qt_het_sum[pidx];
	qt_het_ssq = git_qt_het_ssq[pidx];
	qt_homrar_sum = git_qt_homrar_sum[pidx];
	qt_homrar_ssq = git_qt_homrar_ssq[pidx];
	qt_g_prod = qt_het_sum + 2 * qt_homrar_sum;
	qt_mean = qt_sum * nanal_recip;
	qt_g_prod_centered = qt_g_prod - qt_sum * geno_mean;
	sval = qt_g_prod_centered * qt_g_prod_centered / (geno_mean_sq * (qt_ssq + (qt_mean - 2) * qt_sum) + geno_mean_coeff2 * (qt_het_ssq + qt_mean * (qt_mean * het_ctd - 2 * qt_het_sum)) + geno_mean_coeff3 * (qt_homrar_ssq + qt_mean * (qt_mean * homrar_ctd - 2 * qt_homrar_sum)));
	if (sval > stat_high) {
	  success_2incr += 2;
	} else if (sval > stat_low) {
	  success_2incr++;
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    sval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(sval * (1 - sval) / ((int32_t)next_adapt_check));
	    dyy = sval - dxx;
	    dzz = sval + dxx;
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      fill_double_zero(next_cqg, git_qt_het_sum);
	      fill_double_zero(next_cqg, git_qt_het_ssq);
	      fill_double_zero(next_cqg, git_qt_homrar_sum);
	      fill_double_zero(next_cqg, git_qt_homrar_ssq);
	      fill_double_zero(next_cqg, git_qt_missing_sum);
	      fill_double_zero(next_cqg, git_qt_missing_ssq);
	      goto qassoc_adapt_lin_thread_lesszero;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      fill_double_zero(perm_vec_ctcl8m * 6, git_qt_het_sum);
    qassoc_adapt_lin_thread_lesszero:
      perm_2success_ct[marker_idx] += success_2incr;
    }
  qassoc_adapt_lin_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE qassoc_maxt_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  uint32_t max_thread_ct = g_assoc_thread_ct;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* msa_ptr = nullptr;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  double* git_qt_g_prod;
  double* git_qt_sum;
  double* git_qt_ssq;
  double* qresultbuf;
  double* __restrict__ orig_chiabs;
  uintptr_t* loadbuf;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint16_t* ldrefs;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  uintptr_t marker_idx;
  uint32_t qblock_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t marker_bidx2;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  intptr_t geno_sum;
  intptr_t geno_ssq;
  uint32_t nanal;
  double nanal_recip;
  double nanal_m1_recip;
  double geno_mean;
  double geno_var;
  double qt_sum;
  double qt_ssq;
  double qt_g_prod;
  double qt_mean;
  double qt_var;
  double qt_g_covar;
  double nanal_m2d;
  double beta;
  double betasq;
  double dxx;
  uint32_t success_2incr;
  double stat_high;
  double stat_low;
  double sval;
  uintptr_t best_cost;
  uint32_t cur_thread_ct;
  uint32_t marker_idx_tmp;
  int32_t missing_ct_tmp;
  int32_t het_ct_tmp;
  int32_t homcom_ct_tmp;
  int32_t homrar_ct_tmp;
  uint32_t loop_ceil;
  uintptr_t cur_cost;
  uint32_t ldref;
  while (1) {
    cur_thread_ct = g_block_diff / CACHELINE_DBL;
    if (cur_thread_ct > max_thread_ct) {
      cur_thread_ct = max_thread_ct;
    } else if (!cur_thread_ct) {
      cur_thread_ct = 1;
    }
    if (cur_thread_ct <= tidx) {
      goto qassoc_maxt_thread_skip_all;
    }
    qblock_start = g_qblock_start;
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + qblock_start;
    marker_bidx_start = qblock_start + (((uint64_t)tidx) * g_block_diff) / cur_thread_ct;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    marker_bceil = qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / cur_thread_ct;
    memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
    if (g_mperm_save_all) {
      msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
    }
    loadbuf = g_loadbuf;
    qresultbuf = g_qresultbuf;
    orig_chiabs = g_orig_chisq;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    ldrefs = g_ldrefs;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      missing_ct = missing_cts[marker_idx];
      nanal = pheno_nm_ct - missing_ct;
      homcom_ct = homcom_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      if ((nanal < 3) || (homcom_ct == nanal) || (het_ct == nanal)) {
	perm_2success_ct[marker_idx++] += perm_vec_ct;
	if (msa_ptr) {
	  for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	    *msa_ptr++ = -9;
	  }
	}
	continue;
      }
      homrar_ct = nanal - het_ct - homcom_ct;
      sval = orig_chiabs[marker_idx];
      sval = sval * sval;
      stat_high = sval + EPSILON;
      stat_low = sval - EPSILON;
      geno_sum = 2 * homrar_ct + het_ct;
      geno_ssq = 4 * homrar_ct + het_ct;
      nanal_recip = 1.0 / ((double)((int32_t)nanal));
      nanal_m1_recip = 1.0 / ((double)(((int32_t)nanal) - 1));
      nanal_m2d = nanal - 2;
      geno_mean = ((double)geno_sum) * nanal_recip;
      geno_var = (((double)geno_ssq) - geno_sum * geno_mean) * nanal_m1_recip;
      success_2incr = 0;
      git_qt_g_prod = &(qresultbuf[3 * marker_bidx * perm_vec_ctcl8m]);
      git_qt_sum = &(git_qt_g_prod[perm_vec_ctcl8m]);
      git_qt_ssq = &(git_qt_g_prod[2 * perm_vec_ctcl8m]);
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      ldref = ldrefs[marker_idx];
      if (ldref == 65535) {
	// Addition loops required for genotype indexing:
	//   het_ct + homrar_ct + 2 * missing_ct
	//
	// Addition/initial copy loops required for LD exploitation:
	//   3 + 3 * (missing <-> homrar/het) + 2 * (missing <-> homcom) +
	//   (homrar <-> het/homcom) + (het <-> homcom)
	// Simple lower bound (may allow us to skip full LD cost calculation):
	//   (delta(homrar) + 2*delta(missing) + delta(het) + delta(homcom)) / 2
	best_cost = het_ct + homrar_ct + 2 * missing_ct;
	ldref = marker_bidx;
	marker_idx_tmp = maxt_block_base;
	loop_ceil = maxt_block_base2;
	do {
	  if (marker_idx_tmp == maxt_block_base2) {
	    marker_idx_tmp = maxt_block_base3;
	    loop_ceil = marker_idx;
	  }
	  for (; marker_idx_tmp < loop_ceil; marker_idx_tmp++) {
	    if (ldrefs[marker_idx_tmp] != 65535) {
	      missing_ct_tmp = missing_cts[marker_idx_tmp];
	      homcom_ct_tmp = homcom_cts[marker_idx_tmp];
	      het_ct_tmp = het_cts[marker_idx_tmp];
	      homrar_ct_tmp = pheno_nm_ct - missing_ct_tmp - het_ct_tmp - homcom_ct_tmp;
	      cur_cost = labs(((int32_t)missing_ct) - missing_ct_tmp) + (labs(((int32_t)homrar_ct) - homrar_ct_tmp) + labs(((int32_t)het_ct) - het_ct_tmp) + labs(((int32_t)homcom_ct) - homcom_ct_tmp) + 7) / 2;
	      if (cur_cost < best_cost) {
		marker_bidx2 = marker_idx_tmp - maxt_block_base;
		cur_cost = qrem_cost2(pheno_nm_ctv2, &(loadbuf[marker_bidx2 * pheno_nm_ctv2]), loadbuf_cur);
		if (cur_cost < best_cost) {
		  ldref = marker_bidx2;
		  best_cost = cur_cost;
		}
	      }
	    }
	  }
	} while (marker_idx_tmp < marker_idx);
	ldrefs[marker_idx] = ldref;
      }
      if (ldref == marker_bidx) {
	fill_double_zero(perm_vec_ctcl8m * 3, git_qt_g_prod);
	calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, perm_vec_ct, loadbuf_cur, perm_vecstd, git_qt_g_prod);
      } else {
	memcpy(git_qt_g_prod, &(qresultbuf[3 * ldref * perm_vec_ctcl8m]), 3 * perm_vec_ctcl8m * sizeof(double));
	calc_qrem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecstd, git_qt_g_prod);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	qt_sum = pheno_sum - git_qt_sum[pidx];
	qt_ssq = pheno_ssq - git_qt_ssq[pidx];
	qt_g_prod = git_qt_g_prod[pidx];
	qt_mean = qt_sum * nanal_recip;
	qt_var = (qt_ssq - qt_sum * qt_mean) * nanal_m1_recip;
	qt_g_covar = (qt_g_prod - qt_sum * geno_mean) * nanal_m1_recip;
	dxx = 1.0 / geno_var;
	beta = qt_g_covar * dxx;
	betasq = beta * beta;
	sval = betasq * nanal_m2d / (qt_var * dxx - betasq);
	if (sval > stat_high) {
	  success_2incr += 2;
	} else if (sval > stat_low) {
	  success_2incr++;
	}
	if (results[pidx] < sval) {
	  results[pidx] = sval;
	}
	if (msa_ptr) {
	  *msa_ptr++ = sval;
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  qassoc_maxt_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE qassoc_maxt_lin_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  uint32_t max_thread_ct = g_assoc_thread_ct;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* msa_ptr = nullptr;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  double* git_qt_het_sum;
  double* git_qt_het_ssq;
  double* git_qt_homrar_sum;
  double* git_qt_homrar_ssq;
  double* git_qt_missing_sum;
  double* git_qt_missing_ssq;
  uintptr_t* loadbuf;
  double* qresultbuf;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint16_t* ldrefs;
  double* __restrict__ orig_linsq;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  uintptr_t marker_idx;
  uint32_t qblock_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  intptr_t geno_sum;
  uint32_t nanal;
  uint32_t success_2incr;
  uint32_t cur_thread_ct;
  double het_ctd;
  double homrar_ctd;
  double nanal_recip;
  double geno_mean;
  double geno_mean_sq;
  double geno_mean_coeff2;
  double geno_mean_coeff3;
  double qt_sum;
  double qt_ssq;
  double qt_het_sum;
  double qt_het_ssq;
  double qt_homrar_sum;
  double qt_homrar_ssq;
  double qt_g_prod;
  double qt_mean;
  double qt_g_prod_centered;
  double stat_high;
  double stat_low;
  double sval;
  uint32_t ldref;
  while (1) {
    cur_thread_ct = g_block_diff / CACHELINE_DBL;
    if (cur_thread_ct > max_thread_ct) {
      cur_thread_ct = max_thread_ct;
    } else if (!cur_thread_ct) {
      cur_thread_ct = 1;
    }
    if (cur_thread_ct <= tidx) {
      goto qassoc_maxt_lin_thread_skip_all;
    }
    qblock_start = g_qblock_start;
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + qblock_start;
    marker_bidx_start = qblock_start + (((uint64_t)tidx) * g_block_diff) / cur_thread_ct;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    marker_bceil = qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / cur_thread_ct;
    memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
    if (g_mperm_save_all) {
      msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
    }
    loadbuf = g_loadbuf;
    qresultbuf = g_qresultbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    ldrefs = g_ldrefs;
    orig_linsq = g_orig_linsq;

    for (; marker_bidx < marker_bceil; marker_bidx++) {
      missing_ct = missing_cts[marker_idx];
      nanal = pheno_nm_ct - missing_ct;
      homcom_ct = homcom_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      if ((nanal < 3) || (homcom_ct == nanal) || (het_ct == nanal)) {
	perm_2success_ct[marker_idx++] += perm_vec_ct;
	if (msa_ptr) {
	  for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	    *msa_ptr++ = -9;
	  }
	}
	continue;
      }
      homrar_ct = nanal - het_ct - homcom_ct;
      sval = orig_linsq[marker_idx];
      stat_high = sval + EPSILON;
      stat_low = sval - EPSILON;
      geno_sum = 2 * homrar_ct + het_ct;
      nanal_recip = 1.0 / ((double)((int32_t)nanal));
      het_ctd = het_ct;
      homrar_ctd = homrar_ct;
      geno_mean = ((double)geno_sum) * nanal_recip;
      geno_mean_sq = geno_mean * geno_mean;
      geno_mean_coeff2 = 1 - 2 * geno_mean;
      geno_mean_coeff3 = 4 - 4 * geno_mean;
      success_2incr = 0;
      git_qt_het_sum = &(qresultbuf[6 * marker_bidx * perm_vec_ctcl8m]);
      git_qt_het_ssq = &(git_qt_het_sum[perm_vec_ctcl8m]);
      git_qt_homrar_sum = &(git_qt_het_sum[2 * perm_vec_ctcl8m]);
      git_qt_homrar_ssq = &(git_qt_het_sum[3 * perm_vec_ctcl8m]);
      git_qt_missing_sum = &(git_qt_het_sum[4 * perm_vec_ctcl8m]);
      git_qt_missing_ssq = &(git_qt_het_sum[5 * perm_vec_ctcl8m]);
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      ldref = ldrefs[marker_idx];
      if (ldref == 65535) {
	// 2x addition loops required for genotype indexing:
	//   het_ct + homrar_ct + missing_ct
	//
	// 2x addition/initial copy loops required for LD exploitation:
	//   3 + 2 * (<-> neither side homcom) + (<-> homcom)
	// Simple lower bound (may allow us to skip full LD cost calculation):
	//   3 + delta(homcom) if delta(homcom) >= sum of other deltas
	//   3 + delta(non-homcom) otherwise
	ldref = marker_bidx;
	if (pheno_nm_ct - homcom_ct > 3) {
	  check_for_better_rem_cost(pheno_nm_ct - homcom_ct - 3, maxt_block_base, maxt_block_base2, maxt_block_base3, marker_idx, missing_cts, homcom_cts, het_cts, ldrefs, pheno_nm_ct, missing_ct, het_ct, homcom_ct, loadbuf, loadbuf_cur, &ldref);
	}
	ldrefs[marker_idx] = ldref;
      }
      if (ldref == marker_bidx) {
	fill_double_zero(perm_vec_ctcl8m * 6, git_qt_het_sum);
	calc_qgit_lin(pheno_nm_ct, perm_vec_ctcl8m, perm_vec_ct, loadbuf_cur, perm_vecstd, git_qt_het_sum);
      } else {
	memcpy(git_qt_het_sum, &(qresultbuf[6 * ldref * perm_vec_ctcl8m]), 6 * perm_vec_ctcl8m * sizeof(double));
	calc_qrem_lin(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecstd, git_qt_het_sum);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	qt_sum = pheno_sum - git_qt_missing_sum[pidx];
	qt_ssq = pheno_ssq - git_qt_missing_ssq[pidx];
	qt_het_sum = git_qt_het_sum[pidx];
	qt_het_ssq = git_qt_het_ssq[pidx];
	qt_homrar_sum = git_qt_homrar_sum[pidx];
	qt_homrar_ssq = git_qt_homrar_ssq[pidx];
	qt_g_prod = qt_het_sum + 2 * qt_homrar_sum;
	qt_mean = qt_sum * nanal_recip;
	qt_g_prod_centered = qt_g_prod - qt_sum * geno_mean;
	sval = qt_g_prod_centered * qt_g_prod_centered / (geno_mean_sq * (qt_ssq + (qt_mean - 2) * qt_sum) + geno_mean_coeff2 * (qt_het_ssq + qt_mean * (qt_mean * het_ctd - 2 * qt_het_sum)) + geno_mean_coeff3 * (qt_homrar_ssq + qt_mean * (qt_mean * homrar_ctd - 2 * qt_homrar_sum)));
	if (sval > stat_high) {
	  success_2incr += 2;
	} else if (sval > stat_low) {
	  success_2incr++;
	}
	if (results[pidx] < sval) {
	  results[pidx] = sval;
	}
	if (msa_ptr) {
	  *msa_ptr++ = sval;
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  qassoc_maxt_lin_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE qassoc_set_thread(void* arg) {
  // Simplified version of qassoc_adapt/maxt_thread(), except we need to save
  // actual t-statistics.
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t max_thread_ct = g_assoc_thread_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* git_qt_g_prod = &(g_thread_git_qbufs[perm_vec_ctcl8m * tidx * 3]);
  double* git_qt_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 3 + 1)]);
  double* git_qt_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 3 + 2)]);
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t* __restrict__ loadbuf;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t geno_sum;
  intptr_t geno_ssq;
  double* msa_ptr;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  uint32_t nanal;
  uint32_t cur_thread_ct;
  double nanal_recip;
  double nanal_m1_recip;
  double nanal_m2_recip;
  double geno_mean;
  double geno_var_recip;
  double qt_sum;
  double qt_ssq;
  double qt_g_prod;
  double qt_mean;
  double qt_var;
  double qt_g_covar;
  double beta;
  double vbeta_sqrt;
  while (1) {
    cur_thread_ct = g_block_diff / CACHELINE_DBL;
    if (cur_thread_ct > max_thread_ct) {
      cur_thread_ct = max_thread_ct;
    } else if (!cur_thread_ct) {
      cur_thread_ct = 1;
    }
    if (cur_thread_ct <= tidx) {
      goto qassoc_set_thread_skip_all;
    }
    marker_bidx = (((uint64_t)tidx) * g_block_diff) / cur_thread_ct;
    marker_bceil = (((uint64_t)tidx + 1) * g_block_diff) / cur_thread_ct;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      msa_ptr = &(g_mperm_save_all[marker_bidx * perm_vec_ct]);
      missing_ct = missing_cts[marker_idx];
      nanal = pheno_nm_ct - missing_ct;
      homcom_ct = homcom_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      homrar_ct = nanal - het_ct - homcom_ct;
      geno_sum = 2 * homrar_ct + het_ct;
      geno_ssq = 4 * homrar_ct + het_ct;
      nanal_recip = 1.0 / ((double)((int32_t)nanal));
      nanal_m1_recip = 1.0 / ((double)(((int32_t)nanal) - 1));
      nanal_m2_recip = 1.0 / ((double)(((int32_t)nanal) - 2));
      geno_mean = ((double)geno_sum) * nanal_recip;
      geno_var_recip = 1.0 / ((((double)geno_ssq) - geno_sum * geno_mean) * nanal_m1_recip);
      calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, perm_vec_ct, &(loadbuf[marker_bidx * pheno_nm_ctv2]), perm_vecstd, git_qt_g_prod);
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	qt_sum = pheno_sum - git_qt_sum[pidx];
	qt_ssq = pheno_ssq - git_qt_ssq[pidx];
	qt_g_prod = git_qt_g_prod[pidx];
	qt_mean = qt_sum * nanal_recip;
	qt_var = (qt_ssq - qt_sum * qt_mean) * nanal_m1_recip;
	qt_g_covar = (qt_g_prod - qt_sum * geno_mean) * nanal_m1_recip;
	beta = qt_g_covar * geno_var_recip;
	vbeta_sqrt = sqrt((qt_var * geno_var_recip - beta * beta) * nanal_m2_recip);
	*msa_ptr++ = fabs(beta / vbeta_sqrt);
      }
      fill_double_zero(perm_vec_ctcl8m * 3, git_qt_g_prod);
    }
  qassoc_set_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_adapt_domrec_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_perm_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t* __restrict__ loadbuf;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t missing_start;
  uint32_t case_homx_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  while (1) {
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_adapt_domrec_thread_skip_all;
      }
      marker_bidx = g_block_start + tidx;
      marker_bceil = marker_bidx + 1;
    } else {
      marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    loadbuf = g_loadbuf;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      if (model_fisher) {
        if (orig_pvals[marker_idx] == -9) {
          perm_adapt_stop[marker_idx] = 1;
          perm_attempt_ct[marker_idx] = 0;
          continue;
        }
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
        if (orig_chisq[marker_idx] == -9) {
          perm_adapt_stop[marker_idx] = 1;
          perm_attempt_ct[marker_idx] = 0;
          continue;
        }
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      next_adapt_check = first_adapt_check;
      tot_obs = pheno_nm_ct - missing_cts[marker_idx];
      if (is_model_prec) {
	col2_sum = homcom_cts[marker_idx] + het_cts[marker_idx];
	col1_sum = tot_obs - col2_sum;
      } else {
	col1_sum = homcom_cts[marker_idx];
	col2_sum = tot_obs - col1_sum;
      }
      missing_start = precomp_start[marker_bidx];
      gpui = &(precomp_ui[4 * precomp_width * marker_bidx]);
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      for (pidx = 0; pidx < perm_vec_ct;) {
	genovec_3freq(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), pheno_nm_ctv2, &case_missing_ct, &uii, &case_homx_ct);
	if (is_model_prec) {
	  case_homx_ct = case_ct - case_homx_ct - case_missing_ct - uii;
	}
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_homx_ct < gpui[4 * uii]) {
	    if (case_homx_ct < gpui[4 * uii + 2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (case_homx_ct >= gpui[4 * uii + 1]) {
	      if (case_homx_ct >= gpui[4 * uii + 3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	} else {
	  uii = case_ct - case_missing_ct;
	  if (model_fisher) {
	    dxx = fisher22(case_homx_ct, uii - case_homx_ct, col1_sum - case_homx_ct, col2_sum + case_homx_ct - uii, fisher_midp);
	    if (dxx < stat_low) {
	      success_2incr += 2;
	    } else if (dxx <= stat_high) {
	      success_2incr++;
	    }
	  } else {
	    dxx = chi22_eval(case_homx_ct, uii, col1_sum, tot_obs);
	    if (dxx > stat_high) {
	      success_2incr += 2;
	    } else if (dxx > stat_low) {
	      success_2incr++;
	    }
	  }
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    dyy = pval - dxx; // lower bound
	    dzz = pval + dxx; // upper bound
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  model_adapt_domrec_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_maxt_domrec_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_perm_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ mperm_save_all = g_mperm_save_all;
  double* msa_ptr = nullptr;
  uintptr_t* __restrict__ loadbuf;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  double* __restrict__ precomp_d;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint16_t* ldrefs;
  uintptr_t* loadbuf_cur;
  uint32_t* resultbuf;
  uint32_t* gpui;
  double* gpd;
  uintptr_t pidx;
  uintptr_t marker_idx;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t block_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_homx_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double stat_high;
  double stat_low;
  double sval;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t ldref;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_maxt_domrec_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + block_start;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    precomp_d = g_precomp_d;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    resultbuf = g_resultbuf;
    ldrefs = g_ldrefs;
    memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
    if (mperm_save_all) {
      msa_ptr = &(mperm_save_all[marker_idx * perm_vec_ct]);
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      if (model_fisher) {
	if (orig_pvals[marker_idx] == -9) {
	model_maxt_domrec_thread_skip_marker:
	  marker_idx++;
	  if (msa_ptr) {
	    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	      *msa_ptr++ = -9;
	    }
	  }
	  continue;
	}
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
	if (orig_chisq[marker_idx] == -9) {
	  goto model_maxt_domrec_thread_skip_marker;
	}
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      missing_ct = missing_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      if (is_model_prec) {
	col2_sum = homcom_ct + het_ct;
	col1_sum = tot_obs - col2_sum;
      } else {
	col1_sum = homcom_ct;
	col2_sum = tot_obs - col1_sum;
      }
      missing_start = precomp_start[marker_bidx];
      gpui = &(precomp_ui[6 * precomp_width * marker_bidx]);
      success_2incr = 0;
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      ldref = ldrefs[marker_idx];
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      if (ldref == 65535) {
	ldref = marker_bidx;
	if (pheno_nm_ct - homcom_ct > 50) {
	  check_for_better_rem_cost(pheno_nm_ct - homcom_ct - 50, maxt_block_base, maxt_block_base2, maxt_block_base3, marker_idx, missing_cts, homcom_cts, het_cts, ldrefs, pheno_nm_ct, missing_ct, het_ct, homcom_ct, loadbuf, loadbuf_cur, &ldref);
	}
	ldrefs[marker_idx] = ldref;
      }
      if (ldref == marker_bidx) {
	fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
	calc_git(pheno_nm_ct, perm_vec_ct, &(loadbuf[marker_bidx * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
	fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      } else {
	memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
	calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
	if (is_model_prec) {
	  case_homx_ct = git_homrar_cts[pidx];
	} else {
	  case_homx_ct = case_ct - case_missing_ct - git_homrar_cts[pidx] - git_het_cts[pidx];
	}
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_homx_ct < gpui[6 * uii]) {
	    if (case_homx_ct < gpui[6 * uii + 2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (case_homx_ct >= gpui[6 * uii + 1]) {
	      if (case_homx_ct >= gpui[6 * uii + 3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	  ukk = gpui[6 * uii + 4];
	  ujj = (uint32_t)(case_homx_ct - ukk); // deliberate underflow
	  if (ujj >= gpui[6 * uii + 5]) {
	    if (model_fisher) {
	      ujj = case_ct - case_missing_ct;
	      sval = fisher22_tail_pval(ukk, ujj - ukk, col1_sum - ukk, col2_sum + ukk - ujj, gpui[6 * uii + 5] - 1, gpd[2 * uii], gpd[2 * uii + 1], fisher_midp, case_homx_ct);
	      if (results[pidx] > sval) {
		results[pidx] = sval;
	      }
	    } else {
	      sval = ((double)((intptr_t)case_homx_ct)) - gpd[2 * uii];
	      sval = sval * sval * gpd[2 * uii + 1];
	      if (results[pidx] < sval) {
		results[pidx] = sval;
	      }
	    }
	  }
	} else {
	  uii = case_ct - case_missing_ct;
	  if (model_fisher) {
	    sval = fisher22(case_homx_ct, uii - case_homx_ct, col1_sum - case_homx_ct, col2_sum + case_homx_ct - uii, fisher_midp);
	    if (sval < stat_low) {
	      success_2incr += 2;
	    } else if (sval <= stat_high) {
	      success_2incr++;
	    }
	    if (results[pidx] > sval) {
	      results[pidx] = sval;
	    }
	  } else {
	    sval = chi22_eval(case_homx_ct, uii, col1_sum, tot_obs);
	    if (sval > stat_high) {
	      success_2incr += 2;
	    } else if (sval > stat_low) {
	      success_2incr++;
	    }
	    if (results[pidx] < sval) {
	      results[pidx] = sval;
	    }
	  }
	  if (msa_ptr) {
	    *msa_ptr++ = sval;
	  }
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  model_maxt_domrec_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_set_domrec_thread(void* arg) {
  // Similar to assoc_set_thread().
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t case_ct = g_perm_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  double* msa_ptr = nullptr;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_cur;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t pidx;
  uintptr_t marker_idx;
  intptr_t col1_sum;
  intptr_t tot_obs;
  uint32_t block_start;
  uint32_t marker_bidx_start;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t case_homx_ct;
  uint32_t case_missing_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
	goto model_set_domrec_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    marker_bidx = marker_bidx_start;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      msa_ptr = &(g_mperm_save_all[marker_bidx * perm_vec_ct]);
      missing_ct = missing_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      if (is_model_prec) {
	col1_sum = tot_obs - homcom_ct - het_ct;
      } else {
        col1_sum = homcom_ct;
      }
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
      calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
      fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
        if (is_model_prec) {
	  case_homx_ct = git_homrar_cts[pidx];
	} else {
	  case_homx_ct = case_ct - case_missing_ct - git_homrar_cts[pidx] - git_het_cts[pidx];
	}
	*msa_ptr++ = chi22_eval(case_homx_ct, case_ct - case_missing_ct, col1_sum, tot_obs);
      }
    }
  model_set_domrec_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_adapt_trend_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_perm_case_ct;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t* __restrict__ loadbuf;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  intptr_t tot_obs;
  uint32_t missing_start;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t case_com_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  double chisq_high;
  double chisq_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  while (1) {
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_adapt_trend_thread_skip_all;
      }
      marker_bidx = g_block_start + tidx;
      marker_bceil = marker_bidx + 1;
    } else {
      marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    loadbuf = g_loadbuf;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      next_adapt_check = first_adapt_check;
      if (orig_pvals[marker_idx] == -9) {
	perm_adapt_stop[marker_idx] = 1;
	perm_attempt_ct[marker_idx] = next_adapt_check;
	perm_2success_ct[marker_idx] = next_adapt_check;
	continue;
      }
      tot_obs = pheno_nm_ct - missing_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      missing_start = precomp_start[marker_bidx];
      gpui = &(precomp_ui[4 * precomp_width * marker_bidx]);
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      chisq_high = orig_chisq[marker_idx] + EPSILON;
      chisq_low = orig_chisq[marker_idx] - EPSILON;
      for (pidx = 0; pidx < perm_vec_ct;) {
	genovec_set_freq(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), pheno_nm_ctv2, &case_com_ct, &case_missing_ct);
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_com_ct < gpui[4 * uii]) {
	    if (case_com_ct < gpui[4 * uii + 2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (case_com_ct >= gpui[4 * uii + 1]) {
	      if (case_com_ct >= gpui[4 * uii + 3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	} else {
	  uii = case_ct - case_missing_ct;
	  dxx = ca_trend_eval(case_com_ct, uii, het_ct, homcom_ct, tot_obs);
	  if (dxx > chisq_high) {
	    success_2incr += 2;
	  } else if (dxx > chisq_low) {
	    success_2incr++;
	  }
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    dyy = pval - dxx; // lower bound
	    dzz = pval + dxx; // upper bound
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  model_adapt_trend_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_maxt_trend_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_perm_case_ct;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ mperm_save_all = g_mperm_save_all;
  double* msa_ptr = nullptr;
  uintptr_t* __restrict__ loadbuf;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  double* __restrict__ precomp_d;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint16_t* ldrefs;
  uintptr_t* loadbuf_cur;
  uint32_t* resultbuf;
  uint32_t* gpui;
  double* gpd;
  uint32_t block_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uintptr_t marker_idx;
  uint32_t marker_bceil;
  uintptr_t pidx;
  intptr_t tot_obs;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t ldref;
  uint32_t case_com_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double chisq_high;
  double chisq_low;
  double chisq;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_maxt_trend_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + block_start;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    precomp_d = g_precomp_d;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    resultbuf = g_resultbuf;
    ldrefs = g_ldrefs;
    memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
    if (mperm_save_all) {
      msa_ptr = &(mperm_save_all[marker_idx * perm_vec_ct]);
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      if (orig_pvals[marker_idx] == -9) {
	perm_2success_ct[marker_idx++] += perm_vec_ct;
	continue;
      }
      missing_ct = missing_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      missing_start = precomp_start[marker_bidx];
      gpui = &(precomp_ui[6 * precomp_width * marker_bidx]);
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      chisq_high = orig_chisq[marker_idx] + EPSILON;
      chisq_low = orig_chisq[marker_idx] - EPSILON;
      success_2incr = 0;
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      ldref = ldrefs[marker_idx];
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      if (ldref == 65535) {
	ldref = marker_bidx;
	if (pheno_nm_ct - homcom_ct > 50) {
	  check_for_better_rem_cost(pheno_nm_ct - homcom_ct - 50, maxt_block_base, maxt_block_base2, maxt_block_base3, marker_idx, missing_cts, homcom_cts, het_cts, ldrefs, pheno_nm_ct, missing_ct, het_ct, homcom_ct, loadbuf, loadbuf_cur, &ldref);
	}
	ldrefs[marker_idx] = ldref;
      }
      if (ldref == marker_bidx) {
	fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
	calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
	fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      } else {
	memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
	calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
	case_com_ct = 2 * (case_ct - case_missing_ct - git_homrar_cts[pidx]) - git_het_cts[pidx];
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_com_ct < gpui[6 * uii]) {
	    if (case_com_ct < gpui[6 * uii + 2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (case_com_ct >= gpui[6 * uii + 1]) {
	      if (case_com_ct >= gpui[6 * uii + 3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	  ukk = gpui[6 * uii + 4];
	  ujj = (uint32_t)(case_com_ct - ukk); // deliberate underflow
	  if (ujj >= gpui[6 * uii + 5]) {
	    chisq = ((double)((intptr_t)case_com_ct)) - gpd[2 * uii];
	    chisq = chisq * chisq * gpd[2 * uii + 1];
	    if (results[pidx] < chisq) {
	      results[pidx] = chisq;
	    }
	  }
	} else {
	  chisq = ca_trend_eval(case_com_ct, case_ct - case_missing_ct, het_ct, homcom_ct, tot_obs);
	  if (chisq > chisq_high) {
	    success_2incr += 2;
	  } else if (chisq > chisq_low) {
	    success_2incr++;
	  }
	  if (results[pidx] < chisq) {
	    results[pidx] = chisq;
	  }
	  if (msa_ptr) {
	    *msa_ptr++ = chisq;
	  }
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  model_maxt_trend_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_set_trend_thread(void* arg) {
  // Similar to model_set_domrec_thread().  (In fact, it's so similar that it
  // may be appropriate to merge the functions.)
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t case_ct = g_perm_case_ct;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  double* msa_ptr = nullptr;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_cur;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t pidx;
  uintptr_t marker_idx;
  intptr_t tot_obs;
  uint32_t block_start;
  uint32_t marker_bidx_start;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t case_com_ct;
  uint32_t case_missing_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
	goto model_set_trend_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    marker_bidx = marker_bidx_start;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      msa_ptr = &(g_mperm_save_all[marker_bidx * perm_vec_ct]);
      missing_ct = missing_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
      calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
      fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
	case_com_ct = 2 * (case_ct - case_missing_ct - git_homrar_cts[pidx]) - git_het_cts[pidx];
	*msa_ptr++ = ca_trend_eval(case_com_ct, case_ct - case_missing_ct, het_ct, homcom_ct, tot_obs);
      }
    }
  model_set_trend_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_adapt_gen_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_perm_case_ct;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t* __restrict__ loadbuf;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  uint32_t missing_col;
  intptr_t tot_obs;
  intptr_t homcom_ct;
  intptr_t homrar_ct;
  intptr_t het_ct;
  uint32_t case_missing_ct;
  uint32_t case_het_ct;
  uint32_t case_homcom_ct;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  while (1) {
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_adapt_gen_thread_skip_all;
      }
      marker_bidx = g_block_start + tidx;
      marker_bceil = marker_bidx + 1;
    } else {
      marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    loadbuf = g_loadbuf;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      if (model_fisher) {
        if (orig_pvals[marker_idx] == -9) {
          perm_adapt_stop[marker_idx] = 1;
          perm_attempt_ct[marker_idx] = 0;
          continue;
        }
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
        if (orig_chisq[marker_idx] == -9) {
          perm_adapt_stop[marker_idx] = 1;
          perm_attempt_ct[marker_idx] = 0;
          continue;
        }
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      next_adapt_check = first_adapt_check;
      het_ct = het_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      homrar_ct = tot_obs - het_ct - homcom_ct;
      if (!homcom_ct) {
	missing_col = 3;
      } else if ((het_ct + homcom_ct == tot_obs) || (!het_ct)) {
	missing_col = 2; // either no hom A1s or no hets (no need to distinguish)
      } else {
	missing_col = 0;
      }
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      for (pidx = 0; pidx < perm_vec_ct;) {
	genovec_3freq(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), pheno_nm_ctv2, &case_missing_ct, &case_het_ct, &case_homcom_ct);
	if (model_fisher) {
	  uii = case_ct - case_het_ct - case_homcom_ct - case_missing_ct;
	  // this is very slow.  a precomputed 2-dimensional table could
	  // improve matters, but I doubt it's worth the effort for now.
	  dxx = fisher23(case_homcom_ct, case_het_ct, uii, homcom_ct - case_homcom_ct, het_ct - case_het_ct, homrar_ct - uii, fisher_midp);
	  if (dxx < stat_low) {
	    success_2incr += 2;
	  } else if (dxx <= stat_high) {
	    success_2incr++;
	  }
	} else {
	  if (!missing_col) {
	    dxx = chi23_eval(case_homcom_ct, case_het_ct, case_ct - case_missing_ct, homcom_ct, het_ct, tot_obs);
	  } else if (missing_col == 3) {
	    dxx = chi22_eval(case_het_ct, case_ct - case_missing_ct, het_ct, tot_obs);
	  } else {
	    dxx = chi22_eval(case_homcom_ct, case_ct - case_missing_ct, homcom_ct, tot_obs);
	  }
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  }
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    dyy = pval - dxx; // lower bound
	    dzz = pval + dxx; // upper bound
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  model_adapt_gen_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_maxt_gen_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t case_ct = g_perm_case_ct;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ mperm_save_all = g_mperm_save_all;
  double* msa_ptr = nullptr;
  uintptr_t* __restrict__ loadbuf;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint16_t* ldrefs;
  uintptr_t* loadbuf_cur;
  uint32_t* resultbuf;
  uintptr_t pidx;
  uint32_t missing_col;
  intptr_t tot_obs;
  uintptr_t marker_idx;
  uint32_t block_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  int32_t missing_ct;
  intptr_t homcom_ct;
  intptr_t homrar_ct;
  intptr_t het_ct;
  uint32_t ldref;
  uint32_t success_2incr;
  uint32_t case_missing_ct;
  uint32_t case_het_ct;
  uint32_t case_homcom_ct;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double sval;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_maxt_gen_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + block_start;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    resultbuf = g_resultbuf;
    ldrefs = g_ldrefs;

    memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
    if (mperm_save_all) {
      msa_ptr = &(mperm_save_all[marker_idx * perm_vec_ct]);
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      if (model_fisher) {
	if (orig_pvals[marker_idx] == -9) {
	model_maxt_gen_thread_skip_marker:
	  marker_idx++;
	  if (msa_ptr) {
	    for (pidx = 0; pidx < perm_vec_ct; ++pidx) {
	      *msa_ptr++ = -9;
	    }
	  }
	  continue;
	}
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
	if (orig_chisq[marker_idx] == -9) {
	  goto model_maxt_gen_thread_skip_marker;
	}
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      missing_ct = missing_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      homcom_ct = homcom_cts[marker_idx];
      homrar_ct = tot_obs - het_ct - homcom_ct;
      if (!homcom_ct) {
	missing_col = 3;
      } else if ((het_ct + homcom_ct == tot_obs) || (!het_ct)) {
	missing_col = 2;
      } else {
	missing_col = 0;
      }
      success_2incr = 0;
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      ldref = ldrefs[marker_idx];
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      if (ldref == 65535) {
	ldref = marker_bidx;
	if (pheno_nm_ct - homcom_ct > 50) {
	  check_for_better_rem_cost(pheno_nm_ct - homcom_ct - 50, maxt_block_base, maxt_block_base2, maxt_block_base3, marker_idx, missing_cts, homcom_cts, het_cts, ldrefs, pheno_nm_ct, missing_ct, het_ct, homcom_ct, loadbuf, loadbuf_cur, &ldref);
	}
	ldrefs[marker_idx] = ldref;
      }
      if (ldref == marker_bidx) {
	fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
	calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
	fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      } else {
	memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
	calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
	case_het_ct = git_het_cts[pidx];
	case_homcom_ct = case_ct - case_missing_ct - case_het_ct - git_homrar_cts[pidx];
	if (model_fisher) {
	  uii = case_ct - case_het_ct - case_homcom_ct - case_missing_ct;
	  sval = fisher23(case_homcom_ct, case_het_ct, uii, homcom_ct - case_homcom_ct, het_ct - case_het_ct, homrar_ct - uii, fisher_midp);
	  if (sval < stat_low) {
	    success_2incr += 2;
	  } else if (sval <= stat_high) {
	    success_2incr++;
	  }
	  if (results[pidx] > sval) {
	    results[pidx] = sval;
	  }
	} else {
	  if (!missing_col) {
	    sval = chi23_eval(case_homcom_ct, case_het_ct, case_ct - case_missing_ct, homcom_ct, het_ct, tot_obs);
	  } else if (missing_col == 3) {
	    sval = chi22_eval(case_het_ct, case_ct - case_missing_ct, het_ct, tot_obs);
	  } else {
	    sval = chi22_eval(case_homcom_ct, case_ct - case_missing_ct, homcom_ct, tot_obs);
	  }
	  if (sval > stat_high) {
	    success_2incr += 2;
	  } else if (sval > stat_low) {
	    success_2incr++;
	  }
	  if (results[pidx] < sval) {
	    results[pidx] = sval;
	  }
	}
	if (msa_ptr) {
	  *msa_ptr++ = sval;
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  model_maxt_gen_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_adapt_best_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_perm_case_ct;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t* __restrict__ loadbuf;
  uintptr_t* is_invalid;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  intptr_t tot_obs;
  intptr_t com_ct;
  intptr_t het_ct;
  intptr_t homrar_ct;
  intptr_t homcom_ct;
  uint32_t missing_start;
  uint32_t case_homrar_ct;
  uint32_t case_homcom_ct;
  uint32_t case_het_ct;
  uint32_t case_missing_ct;
  uint32_t case_com_ct;
  uint32_t skip_domrec;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  while (1) {
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_adapt_best_thread_skip_all;
      }
      marker_bidx = g_block_start + tidx;
      marker_bceil = marker_bidx + 1;
    } else {
      marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    loadbuf = g_loadbuf;
    is_invalid = g_is_invalid_bitfield;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;

    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      if (model_fisher) {
        if (orig_pvals[marker_idx] == -9) {
	  perm_adapt_stop[marker_idx] = 1;
	  perm_attempt_ct[marker_idx] = 0;
	  continue;
        }
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      } else {
	if (orig_chisq[marker_idx] == -9) {
	  perm_adapt_stop[marker_idx] = 1;
	  perm_attempt_ct[marker_idx] = 0;
	  continue;
	}
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
      }
      next_adapt_check = first_adapt_check;
      tot_obs = pheno_nm_ct - missing_cts[marker_idx];
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      com_ct = homcom_ct * 2 + het_ct;
      homrar_ct = tot_obs - het_ct - homcom_ct;
      missing_start = precomp_start[marker_bidx];
      skip_domrec = IS_SET(is_invalid, marker_idx);
      gpui = &(precomp_ui[12 * precomp_width * marker_bidx]);
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      for (pidx = 0; pidx < perm_vec_ct;) {
	genovec_3freq(&(loadbuf[marker_bidx * pheno_nm_ctv2]), &(perm_vecs[pidx * pheno_nm_ctv2]), pheno_nm_ctv2, &case_missing_ct, &case_het_ct, &case_homcom_ct);
	case_homrar_ct = case_ct - case_missing_ct - case_het_ct - case_homcom_ct;
	case_com_ct = case_het_ct + 2 * case_homcom_ct;
	ujj = 0; // best increment so far
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  if (case_com_ct < gpui[12 * uii]) {
	    if (case_com_ct < gpui[12 * uii + 2]) {
	      goto model_adapt_best_thread_betterstat;
	    } else {
	      ujj = 1;
	    }
	  } else {
	    if (case_com_ct >= gpui[12 * uii + 1]) {
	      if (case_com_ct >= gpui[12 * uii + 3]) {
		goto model_adapt_best_thread_betterstat;
	      } else {
		ujj = 1;
	      }
	    }
	  }
	  if (!skip_domrec) {
	    if (case_homcom_ct < gpui[12 * uii + 4]) {
	      if (case_homcom_ct < gpui[12 * uii + 6]) {
		goto model_adapt_best_thread_betterstat;
	      } else {
		ujj = 1;
	      }
	    } else {
	      if (case_homcom_ct >= gpui[12 * uii + 5]) {
		if (case_homcom_ct >= gpui[12 * uii + 7]) {
		  goto model_adapt_best_thread_betterstat;
		} else {
		  ujj = 1;
		}
	      }
	    }
	    if (case_homrar_ct < gpui[12 * uii + 8]) {
	      if (case_homrar_ct < gpui[12 * uii + 10]) {
		goto model_adapt_best_thread_betterstat;
	      } else {
		ujj = 1;
	      }
	    } else {
	      if (case_homrar_ct >= gpui[12 * uii + 9]) {
		if (case_homrar_ct >= gpui[12 * uii + 11]) {
		  goto model_adapt_best_thread_betterstat;
		} else {
		  ujj = 1;
		}
	      }
	    }
	  }
	} else if (1) {
	  uii = case_ct - case_missing_ct; // nonmissing cases
	  if (model_fisher) {
	    ukk = tot_obs - uii; // nonmissing controls
	    dxx = fisher22(case_com_ct, 2 * uii - case_com_ct, com_ct - case_com_ct, 2 * ukk + case_com_ct - com_ct, fisher_midp);
	    if (dxx < stat_low) {
	      goto model_adapt_best_thread_betterstat;
	    } else if (dxx <= stat_high) {
	      ujj = 1;
	    }
	    if (!skip_domrec) {
	      dxx = fisher22(case_homcom_ct, uii - case_homcom_ct, homcom_ct - case_homcom_ct, ukk + case_homcom_ct - homcom_ct, fisher_midp);
	      if (dxx < stat_low) {
		goto model_adapt_best_thread_betterstat;
	      } else if (dxx <= stat_high) {
		ujj = 1;
	      }
	      dxx = fisher22(case_homrar_ct, uii - case_homrar_ct, homrar_ct - case_homrar_ct, ukk + case_homrar_ct - homrar_ct, fisher_midp);
	      if (dxx < stat_low) {
		goto model_adapt_best_thread_betterstat;
	      } else if (dxx <= stat_high) {
		ujj = 1;
	      }
	    }
	  } else {
	    dxx = chi22_eval(case_com_ct, 2 * uii, com_ct, 2 * tot_obs);
	    if (dxx > stat_high) {
	      goto model_adapt_best_thread_betterstat;
	    } else if (dxx > stat_low) {
	      ujj = 1;
	    }
	    if (!skip_domrec) {
	      dxx = chi22_eval(case_homcom_ct, uii, homcom_ct, tot_obs);
	      if (dxx > stat_high) {
		goto model_adapt_best_thread_betterstat;
	      } else if (dxx > stat_low) {
		ujj = 1;
	      }
	      dxx = chi22_eval(case_homrar_ct, uii, homrar_ct, tot_obs);
	      if (dxx > stat_high) {
		goto model_adapt_best_thread_betterstat;
	      } else if (dxx > stat_low) {
		ujj = 1;
	      }
	    }
	  }
	} else {
	model_adapt_best_thread_betterstat:
	  ujj = 2;
	}
	success_2incr += ujj;
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    dyy = pval - dxx; // lower bound
	    dzz = pval + dxx; // upper bound
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  model_adapt_best_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_maxt_best_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t fisher_midp = g_fisher_midp;
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_perm_case_ct;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ mperm_save_all = g_mperm_save_all;
  double* msa_ptr = nullptr;
  uintptr_t* __restrict__ loadbuf;
  uintptr_t* is_invalid;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uint32_t* __restrict__ precomp_start;
  uint32_t* __restrict__ precomp_ui;
  double* __restrict__ precomp_d;
  double* __restrict__ orig_pvals;
  double* __restrict__ orig_chisq;
  uint16_t* ldrefs;
  uintptr_t* loadbuf_cur;
  uint32_t* resultbuf;
  uint32_t* gpui;
  double* gpd;
  uintptr_t pidx;
  uintptr_t marker_idx;
  int32_t missing_ct;
  intptr_t tot_obs;
  intptr_t com_ct;
  intptr_t rar_ct;
  intptr_t het_ct;
  intptr_t homrar_ct;
  intptr_t homcom_ct;
  uint32_t block_start;
  uint32_t maxt_block_base;
  uint32_t maxt_block_base2;
  uint32_t marker_bidx_start;
  uint32_t maxt_block_base3;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t ldref;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_homrar_ct;
  uint32_t case_homcom_ct;
  uint32_t case_het_ct;
  uint32_t case_missing_ct;
  uint32_t case_com_ct;
  uint32_t skip_domrec;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t cur_add;
  double stat_high;
  double stat_low;
  double sval;
  double best_stat;
  double default_best_stat;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
        goto model_maxt_best_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    maxt_block_base = g_maxt_block_base;
    maxt_block_base2 = maxt_block_base + block_start;
    maxt_block_base3 = maxt_block_base + marker_bidx_start;
    marker_bidx = marker_bidx_start;
    marker_idx = maxt_block_base3;
    loadbuf = g_loadbuf;
    is_invalid = g_is_invalid_bitfield;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    precomp_start = g_precomp_start;
    precomp_ui = g_precomp_ui;
    precomp_d = g_precomp_d;
    orig_pvals = g_orig_pvals;
    orig_chisq = g_orig_chisq;
    resultbuf = g_resultbuf;
    ldrefs = g_ldrefs;

    memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
    if (mperm_save_all) {
      msa_ptr = &(mperm_save_all[marker_idx * perm_vec_ct]);
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      if (model_fisher) {
	if (orig_pvals[marker_idx] == -9) {
	  marker_idx++;
	  continue;
	}
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
	default_best_stat = 1;
      } else {
	if (orig_chisq[marker_idx] == -9) {
	  marker_idx++;
	  continue;
	}
	stat_high = orig_chisq[marker_idx] + EPSILON;
	stat_low = orig_chisq[marker_idx] - EPSILON;
	default_best_stat = 0;
      }
      gpd = &(precomp_d[6 * precomp_width * marker_bidx]);
      missing_ct = missing_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      com_ct = 2 * homcom_ct + het_ct;
      rar_ct = tot_obs * 2 - com_ct;
      homrar_ct = tot_obs - homcom_ct - het_ct;
      missing_start = precomp_start[marker_bidx];
      skip_domrec = IS_SET(is_invalid, marker_idx);
      gpui = &(precomp_ui[18 * precomp_width * marker_bidx]);
      success_2incr = 0;
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      ldref = ldrefs[marker_idx];
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      if (ldref == 65535) {
	ldref = marker_bidx;
	if (pheno_nm_ct - homcom_ct > 50) {
	  check_for_better_rem_cost(pheno_nm_ct - homcom_ct - 50, maxt_block_base, maxt_block_base2, maxt_block_base3, marker_idx, missing_cts, homcom_cts, het_cts, ldrefs, pheno_nm_ct, missing_ct, het_ct, homcom_ct, loadbuf, loadbuf_cur, &ldref);
	}
	ldrefs[marker_idx] = ldref;
      }
      if (ldref == marker_bidx) {
	fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
	calc_git(pheno_nm_ct, perm_vec_ct, &(loadbuf[marker_bidx * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
	fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      } else {
	memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
	calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctv2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
      }
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
	case_het_ct = git_het_cts[pidx];
	case_homrar_ct = git_homrar_cts[pidx];
	case_homcom_ct = case_ct - case_missing_ct - case_het_ct - case_homrar_ct;
	case_com_ct = case_het_ct + 2 * case_homcom_ct;
	cur_add = 0;
	// deliberate underflow
	uii = (uint32_t)(case_missing_ct - missing_start);
	if (uii < precomp_width) {
	  best_stat = default_best_stat;
	  if (case_com_ct < gpui[18 * uii]) {
	    if (case_com_ct < gpui[18 * uii + 2]) {
	      cur_add = 2;
	    } else {
	      cur_add = 1;
	    }
	  } else {
	    if (case_com_ct >= gpui[18 * uii + 1]) {
	      if (case_com_ct >= gpui[18 * uii + 3]) {
		cur_add = 2;
	      } else {
		cur_add = 1;
	      }
	    }
	  }
	  ukk = gpui[18 * uii + 4];
	  ujj = (uint32_t)(case_com_ct - ukk); // deliberate underflow
	  if (ujj >= gpui[18 * uii + 5]) {
	    if (model_fisher) {
	      ujj = 2 * (case_ct - case_missing_ct);
	      best_stat = fisher22_tail_pval(ukk, ujj - ukk, com_ct - ukk, rar_ct + ukk - ujj, gpui[18 * uii + 5] - 1, gpd[6 * uii], gpd[6 * uii + 1], fisher_midp, case_com_ct);
	    } else {
	      best_stat = ((double)((intptr_t)case_com_ct)) - gpd[6 * uii];
	      best_stat = best_stat * best_stat * gpd[6 * uii + 1];
	    }
	  }
	  if (!skip_domrec) {
	    if (cur_add != 2) {
	      if (case_homcom_ct < gpui[18 * uii + 6]) {
		if (case_homcom_ct < gpui[18 * uii + 8]) {
		  goto model_maxt_best_thread_domrec2;
		} else {
		  cur_add = 1;
		}
	      } else {
		if (case_homcom_ct >= gpui[18 * uii + 7]) {
		  if (case_homcom_ct >= gpui[18 * uii + 9]) {
		    goto model_maxt_best_thread_domrec2;
		  } else {
		    cur_add = 1;
		  }
		}
	      }
	      if (1) {
		if (case_homrar_ct < gpui[18 * uii + 12]) {
		  if (case_homrar_ct < gpui[18 * uii + 14]) {
		    goto model_maxt_best_thread_domrec2;
		  } else {
		    cur_add = 1;
		  }
		} else {
		  if (case_homrar_ct >= gpui[18 * uii + 13]) {
		    if (case_homrar_ct >= gpui[18 * uii + 15]) {
		      goto model_maxt_best_thread_domrec2;
		    } else {
		      cur_add = 1;
		    }
		  }
		}
	      } else {
	      model_maxt_best_thread_domrec2:
		cur_add = 2;
	      }
	    }
	    ukk = gpui[18 * uii + 10];
	    ujj = (uint32_t)(case_homcom_ct - ukk); // deliberate underflow
	    if (ujj >= gpui[18 * uii + 11]) {
	      if (model_fisher) {
		ujj = case_ct - case_missing_ct;
		sval = fisher22_tail_pval(ukk, ujj - ukk, homcom_ct - ukk, homrar_ct + het_ct + ukk - ujj, gpui[18 * uii + 11] - 1, gpd[6 * uii + 2], gpd[6 * uii + 3], fisher_midp, case_homcom_ct);
		if (sval < best_stat) {
		  best_stat = sval;
		}
	      } else {
		sval = ((double)((intptr_t)case_homcom_ct)) - gpd[6 * uii + 2];
		sval = sval * sval * gpd[6 * uii + 3];
		if (sval > best_stat) {
		  best_stat = sval;
		}
	      }
	    }
	    ukk = gpui[18 * uii + 16];
	    ujj = (uint32_t)(case_homrar_ct - ukk); // deliberate underflow
	    if (ujj >= gpui[18 * uii + 17]) {
	      if (model_fisher) {
		ujj = case_ct - case_missing_ct;
		sval = fisher22_tail_pval(ukk, ujj - ukk, homrar_ct - ukk, homcom_ct + het_ct + ukk - ujj, gpui[18 * uii + 17] - 1, gpd[6 * uii + 4], gpd[6 * uii + 5], fisher_midp, case_homrar_ct);
		if (sval < best_stat) {
		  best_stat = sval;
		}
	      } else {
		sval = ((double)((intptr_t)case_homrar_ct)) - gpd[6 * uii + 4];
		sval = sval * sval * gpd[6 * uii + 5];
		if (sval > best_stat) {
		  best_stat = sval;
		}
	      }
	    }
	  }
	} else {
	  uii = case_ct - case_missing_ct;
	  if (model_fisher) {
	    ukk = tot_obs - uii;
	    best_stat = fisher22(case_com_ct, 2 * uii - case_com_ct, com_ct - case_com_ct, 2 * ukk + case_com_ct - com_ct, fisher_midp);
	    if (!skip_domrec) {
	      sval = fisher22(case_homcom_ct, uii - case_homcom_ct, homcom_ct - case_homcom_ct, ukk + case_homcom_ct - homcom_ct, fisher_midp);
	      if (sval < best_stat) {
		best_stat = sval;
	      }
	      sval = fisher22(case_homrar_ct, uii - case_homrar_ct, homrar_ct - case_homrar_ct, ukk + case_homrar_ct - homrar_ct, fisher_midp);
	      if (sval < best_stat) {
		best_stat = sval;
	      }
	    }
	    if (best_stat < stat_low) {
	      cur_add = 2;
	    } else if (best_stat <= stat_high) {
	      cur_add = 1;
	    }
	  } else {
	    best_stat = chi22_eval(case_com_ct, 2 * uii, com_ct, 2 * tot_obs);
	    if (!skip_domrec) {
	      sval = chi22_eval(case_homcom_ct, uii, homcom_ct, tot_obs);
	      if (sval > best_stat) {
		best_stat = sval;
	      }
	      sval = chi22_eval(case_homrar_ct, uii, homrar_ct, tot_obs);
	      if (sval > best_stat) {
		best_stat = sval;
	      }
	    }
	    if (best_stat > stat_high) {
	      cur_add = 2;
	    } else if (best_stat > stat_low) {
	      cur_add = 1;
	    }
	  }
	  if (msa_ptr) {
	    *msa_ptr++ = best_stat;
	  }
	}
	success_2incr += cur_add;
	if (model_fisher) {
	  if (results[pidx] > best_stat) {
	    results[pidx] = best_stat;
	  }
	} else {
	  if (results[pidx] < best_stat) {
	    results[pidx] = best_stat;
	  }
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  model_maxt_best_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE model_set_best_thread(void* arg) {
  // Similar to model_set_domrec_thread().
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uint32_t assoc_thread_ct = g_assoc_thread_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uint32_t perm_ctvc = BITCT_TO_VECCT(perm_vec_ct);
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ctvc * 144 * BYTECT4]);
  uint32_t* git_homrar_cts = nullptr;
  uint32_t* git_missing_cts = nullptr;
  uint32_t* git_het_cts = nullptr;
  uintptr_t perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t case_ct = g_perm_case_ct;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  double* msa_ptr = nullptr;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_cur;
  uintptr_t* is_invalid;
  uint32_t* __restrict__ missing_cts;
  uint32_t* __restrict__ het_cts;
  uint32_t* __restrict__ homcom_cts;
  uintptr_t pidx;
  uintptr_t marker_idx;
  intptr_t tot_obs;
  intptr_t com_ct;
  intptr_t het_ct;
  intptr_t homrar_ct;
  intptr_t homcom_ct;
  double best_stat;
  double sval;
  uint32_t block_start;
  uint32_t marker_bidx_start;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t case_homrar_ct;
  uint32_t case_homcom_ct;
  uint32_t case_het_ct;
  uint32_t case_missing_ct;
  uint32_t case_com_ct;
  uint32_t skip_domrec;
  uint32_t uii;
  int32_t missing_ct;
  while (1) {
    block_start = g_block_start;
    if (g_block_diff <= assoc_thread_ct) {
      if (g_block_diff <= tidx) {
	goto model_set_best_thread_skip_all;
      }
      marker_bidx_start = block_start + tidx;
      marker_bceil = marker_bidx_start + 1;
    } else {
      marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / assoc_thread_ct;
      marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / assoc_thread_ct;
    }
    marker_bidx = marker_bidx_start;
    loadbuf = g_loadbuf;
    is_invalid = g_is_invalid_bitfield;
    missing_cts = g_missing_cts;
    het_cts = g_het_cts;
    homcom_cts = g_homcom_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      msa_ptr = &(g_mperm_save_all[marker_bidx * perm_vec_ct]);
      missing_ct = missing_cts[marker_idx];
      tot_obs = pheno_nm_ct - missing_ct;
      het_ct = het_cts[marker_idx];
      homcom_ct = homcom_cts[marker_idx];
      com_ct = 2 * homcom_ct + het_ct;
      homrar_ct = tot_obs - homcom_ct - het_ct;
      skip_domrec = IS_SET(is_invalid, marker_idx);
      loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctv2]);
      git_homrar_cts = &(resultbuf[3 * marker_bidx * perm_vec_ctcl4m]);
      git_missing_cts = &(git_homrar_cts[perm_vec_ctcl4m]);
      git_het_cts = &(git_homrar_cts[2 * perm_vec_ctcl4m]);
      fill_uint_zero(3 * perm_vec_ctcl4m, git_homrar_cts);
      calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
      fill_uint_zero(perm_ctvc * 72 * BYTECT4, thread_git_wkspace);
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	case_missing_ct = git_missing_cts[pidx];
	case_het_ct = git_het_cts[pidx];
	case_homrar_ct = git_homrar_cts[pidx];
	case_homcom_ct = case_ct - case_missing_ct - case_het_ct - case_homrar_ct;
	case_com_ct = case_het_ct + 2 * case_homcom_ct;
	uii = case_ct - case_missing_ct;
	best_stat = chi22_eval(case_com_ct, 2 * uii, com_ct, 2 * tot_obs);
	if (!skip_domrec) {
          sval = chi22_eval(case_homcom_ct, uii, homcom_ct, tot_obs);
	  if (sval > best_stat) {
	    best_stat = sval;
	  }
	  sval = chi22_eval(case_homrar_ct, uii, homrar_ct, tot_obs);
          if (sval > best_stat) {
            best_stat = sval;
	  }
	}
	*msa_ptr++ = best_stat;
      }
    }
  model_set_best_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

int32_t model_assoc_set_test(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* outname_end2, uint32_t model_modifier, uint32_t model_mperm_val, double pfilter, double output_min_p, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, uintptr_t* marker_exclude_mid, uintptr_t marker_ct_mid, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sex_male, Aperm_info* apip, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* founder_pnm, uint32_t gender_req, uint32_t ld_ignore_x, uint32_t hh_exists, uint32_t perm_batch_size, Set_info* sip, uintptr_t* loadbuf_raw) {
  // Could reuse more of the code in model_assoc() since there's considerable
  // overlap, but there are enough differences between the regular and set
  // permutation tests that separating this out and doing a fair bit of
  // cut-and-paste is justifiable (especially for the first version of this
  // function).

  // There are three levels of marker subsets here.
  // 1. marker_exclude_orig refers to all markers which passed QC filters, etc.
  //    This is needed to interpret the main set data structure.
  // 2. marker_exclude_mid refers to all markers contained in at least one set.
  //    This is a subset of marker_exclude_orig.  (They are identical if
  //    --gene-all was specified.)  It was used during the single-marker
  //    association test phase, and describes which markers orig_chisq[],
  //    g_missing_cts[], etc. elements initially refer to.
  // 3. Finally, the marker_exclude used for set-based permutation testing
  //    refers to all markers contained in at least one *significant* set.
  //    orig_chisq is collapsed before permutation to be congruent to this
  //    marker_exclude.
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t* marker_exclude = marker_exclude_mid;
  uintptr_t* unstopped_markers = nullptr;
  uintptr_t* loadbuf = g_loadbuf;
  uintptr_t* sample_male_include2 = g_sample_male_include2;
  uintptr_t* perm_adapt_set_unstopped = nullptr;
  char* tbuf2 = &(g_textbuf[MAXLINELEN]);
  double* orig_chisq = g_orig_chisq;
  double* sorted_chisq_buf = nullptr;
  uint32_t* marker_idx_to_uidx = nullptr;
  uint32_t* sorted_marker_idx_buf = nullptr;
  uint32_t* proxy_arr = nullptr;
  uint32_t* perm_2success_ct = nullptr;
  uint32_t* perm_attempt_ct = nullptr;
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t marker_ct = marker_ct_mid;
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uintptr_t ulii = 0;
  double adaptive_ci_zt = 0.0;
  uint32_t model_assoc = model_modifier & MODEL_ASSOC;
  uint32_t perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t model_perm_best = !(model_modifier & MODEL_PMASK);
  uint32_t max_thread_ct = g_thread_ct;
  uint32_t perms_done = 0;
  int32_t x_code = chrom_info_ptr->xymt_codes[X_OFFSET];
  int32_t retval = 0;
  uintptr_t* set_incl;
  uintptr_t* loadbuf_ptr;
  double* orig_set_scores;
  double* chisq_pmajor;
  double* chisq_ptr;
  double* read_dptr;
  double* write_dptr;
  unsigned char* bigstack_mark2;
  uint32_t** setdefs;
  uint32_t** ld_map;
  uintptr_t marker_uidx;
  uintptr_t marker_midx;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t set_ct;
  uintptr_t set_idx;
  uintptr_t perm_vec_ct;
  uintptr_t perm_vec_ctcl4m;
  uintptr_t pidx;
  double chisq_threshold;
  double dxx;
  uint32_t perms_total;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t assoc_thread_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t min_ploidy_1;
  uint32_t marker_unstopped_ct;
  uint32_t is_last_block;
  uint32_t first_adapt_check;
  uint32_t max_sigset_size;
  uint32_t marker_bidx;
  uint32_t uii;
  if (sip->set_test_lambda > 1.0) {
    dxx = 1.0 / sip->set_test_lambda;
    chisq_ptr = orig_chisq;
    for (marker_midx = 0; marker_midx < marker_ct; marker_midx++) {
      *chisq_ptr *= dxx;
      chisq_ptr++;
    }
  }
  ulii = (uintptr_t)(outname_end - outname);
  // don't want to overwrite .assoc extension, etc.
  memcpy(tbuf2, outname, ulii);
  retval = set_test_common_init(threads, bedfile, bed_offset, tbuf2, &(tbuf2[ulii]), unfiltered_marker_ct, marker_exclude_orig, marker_ct_orig, marker_ids, max_marker_id_len, marker_reverse, orig_chisq, sip, chrom_info_ptr, unfiltered_sample_ct, sex_male, founder_pnm, ld_ignore_x, hh_exists, "--assoc/--model", &marker_ct, &marker_exclude, &set_incl, &marker_idx_to_uidx, &setdefs, &set_ct, &max_sigset_size, &ld_map, &chisq_threshold, &orig_set_scores, &sorted_chisq_buf, &sorted_marker_idx_buf, &proxy_arr, &perm_adapt_set_unstopped, &perm_2success_ct, &perm_attempt_ct, &unstopped_markers);
  if (retval) {
    goto model_assoc_set_test_ret_1;
  }
  if (!set_ct) {
    goto model_assoc_set_test_write;
  }
  if (marker_ct_mid != marker_ct) {
    // collapse these arrays so the permutation inner loop is faster
    inplace_delta_collapse_arr((char*)g_missing_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
    if (model_assoc) {
      inplace_delta_collapse_arr((char*)g_set_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
    } else {
      inplace_delta_collapse_arr((char*)g_het_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
      inplace_delta_collapse_arr((char*)g_homcom_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
      if (model_perm_best) {
	inplace_delta_collapse_bitfield(g_is_invalid_bitfield, marker_ct, marker_exclude_mid, marker_exclude);
      }
    }
  }

  if (model_modifier & MODEL_PERM) {
    perms_total = apip->max;
    first_adapt_check = (apip->min < apip->init_interval)? ((int32_t)apip->init_interval) : apip->min;
    adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)set_ct)));
  } else {
    perms_total = model_mperm_val;
    first_adapt_check = perms_total + 1;
  }
  for (uii = 0; uii < set_ct; uii++) {
    perm_attempt_ct[uii] = perms_total;
  }
  if (max_thread_ct > perms_total) {
    max_thread_ct = perms_total;
  }
  if (bigstack_init_sfmtp(max_thread_ct)) {
    goto model_assoc_set_test_ret_NOMEM;
  }
  marker_unstopped_ct = marker_ct;
  g_block_start = 0; // will be nonzero sometimes after LD-exploitation added

  // generate a permutation batch, efficiently compute chi-square stats for all
  // variants in at least one tested set, compute set score, compare to base
  // set score.
  bigstack_mark2 = g_bigstack_base;
 model_assoc_set_test_more_perms:
  if (perms_done) {
    uii = apip->init_interval;
    while (first_adapt_check <= perms_done) {
      first_adapt_check += (int32_t)(uii + ((int32_t)first_adapt_check) * apip->interval_slope);
    }
  }
  // perm_vec_ct memory allocation dependencies:
  //   g_perm_vecst: 16 * ((perm_vec_ct + 127) / 128) * pheno_nm_ct
  //   g_thread_git_wkspace: ((perm_vec_ct + 127) / 128) * 1152 * thread_ct
  //   g_resultbuf: MODEL_BLOCKSIZE * (4 * perm_vec_ct, CL-aligned) * 3
  //   g_perm_vecs: pheno_nm_ctv2 * sizeof(intptr_t) * perm_vec_ct
  //   g_mperm_save_all: MODEL_BLOCKSIZE * 8 * perm_vec_ct
  //   chisq_pmajor: marker_ct * 8 * perm_vec_ct
  // If we force perm_vec_ct to be a multiple of 128, then we have
  //   perm_vec_ct * (9 * max_thread_ct + 20 * MODEL_BLOCKSIZE +
  //                    pheno_nm_ct / 8 + sizeof(intptr_t) * pheno_nm_ctv2
  //                    + marker_ct * sizeof(double))
  perm_vec_ct = 128 * (bigstack_left() / (128LL * sizeof(intptr_t) * pheno_nm_ctv2 + 1152LL * max_thread_ct + 2560LL * MODEL_BLOCKSIZE + 16LL * pheno_nm_ct + 128LL * sizeof(double) * marker_ct));
  if (perm_vec_ct > perm_batch_size) {
    perm_vec_ct = perm_batch_size;
  }
  if (perm_vec_ct > perms_total - perms_done) {
    perm_vec_ct = perms_total - perms_done;
  } else if (!perm_vec_ct) {
    goto model_assoc_set_test_ret_NOMEM;
  }
  perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
  perms_done += perm_vec_ct;
  g_perms_done = perms_done;
  g_perm_vec_ct = perm_vec_ct;
  bigstack_alloc_ul(perm_vec_ct * pheno_nm_ctv2, &g_perm_vecs);
  g_perm_generation_thread_ct = MINV(max_thread_ct, perm_vec_ct);
  ulii = 0;
  if (!g_perm_cluster_starts) {
    if (spawn_threads(threads, &generate_cc_perms_thread, g_perm_generation_thread_ct)) {
      goto model_assoc_set_test_ret_THREAD_CREATE_FAIL;
    }
    generate_cc_perms_thread((void*)ulii);
  } else {
    if (spawn_threads(threads, &generate_cc_cluster_perms_thread, g_perm_generation_thread_ct)) {
      goto model_assoc_set_test_ret_THREAD_CREATE_FAIL;
    }
    generate_cc_cluster_perms_thread((void*)ulii);
  }
  join_threads(threads, g_perm_generation_thread_ct);
  g_assoc_thread_ct = max_thread_ct;
  bigstack_alloc_ui(perm_vec_ctcl4m * 3 * MODEL_BLOCKSIZE, &g_resultbuf);
#ifdef __LP64__
  ulii = ((perm_vec_ct + 127) / 128) * 4;
  bigstack_alloc_ui(ulii * pheno_nm_ct, &g_perm_vecst);
#else
  ulii = (perm_vec_ct + 31) / 32;
  bigstack_alloc_ui(ulii * pheno_nm_ct, &g_perm_vecst);
  ulii = ((perm_vec_ct + 63) / 64) * 2;
#endif
  bigstack_calloc_ui(ulii * 72 * max_thread_ct, &g_thread_git_wkspace);
  transpose_perms(g_perm_vecs, perm_vec_ct, pheno_nm_ct, g_perm_vecst);
  bigstack_alloc_d(MODEL_BLOCKSIZE * perm_vec_ct, &g_mperm_save_all);
  bigstack_alloc_d(marker_ct * perm_vec_ct, &chisq_pmajor);
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto model_assoc_set_test_ret_READ_FAIL;
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  do {
    if (marker_uidx >= chrom_end) {
      if (model_assoc) {
	// exploit overflow
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
	min_ploidy_1 |= uii; // treat MT as haploid
	g_min_ploidy_1 = min_ploidy_1;
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	g_is_y = is_y;
      } else {
	// no need to skip MT/haploid here, since we error out on that case
	// earlier
	do {
	  chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	is_x = (uii == (uint32_t)x_code);
      }
      g_is_x = is_x;
    }
    block_size = g_block_start;
    block_end = marker_unstopped_ct - marker_idx;
    if (block_end > MODEL_BLOCKSIZE) {
      block_end = MODEL_BLOCKSIZE;
    }
    do {
      if (!IS_SET(unstopped_markers, marker_idx2)) {
	do {
	  marker_uidx++;
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && (!IS_SET(unstopped_markers, marker_idx2)));
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto model_assoc_set_test_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(loadbuf[block_size * pheno_nm_ctv2]);
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf_ptr)) {
	goto model_assoc_set_test_ret_READ_FAIL;
      }
      g_adapt_m_table[block_size] = marker_idx2++;
      if (is_x && (!model_assoc)) {
	force_missing((unsigned char*)(&(loadbuf[block_size * pheno_nm_ctv2])), sample_male_include2, pheno_nm_ct);
      }
      block_size++;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto model_assoc_set_test_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (!block_size) {
      continue;
    }
    g_block_diff = block_size;
    assoc_thread_ct = g_block_diff;
    if (assoc_thread_ct > max_thread_ct) {
      assoc_thread_ct = max_thread_ct;
    }
    is_last_block = (marker_idx + block_size == marker_unstopped_ct);
    ulii = 0;
    if (model_assoc) {
      if (spawn_threads2(threads, &assoc_set_thread, max_thread_ct, is_last_block)) {
	goto model_assoc_set_test_ret_THREAD_CREATE_FAIL;
      }
      assoc_set_thread((void*)ulii);
    } else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
      if (spawn_threads2(threads, &model_set_domrec_thread, max_thread_ct, is_last_block)) {
	goto model_assoc_set_test_ret_THREAD_CREATE_FAIL;
      }
      model_set_domrec_thread((void*)ulii);
    } else if (model_modifier & MODEL_PTREND) {
      if (spawn_threads2(threads, &model_set_trend_thread, max_thread_ct, is_last_block)) {
	goto model_assoc_set_test_ret_THREAD_CREATE_FAIL;
      }
      model_set_trend_thread((void*)ulii);
    } else {
      if (spawn_threads2(threads, &model_set_best_thread, max_thread_ct, is_last_block)) {
	goto model_assoc_set_test_ret_THREAD_CREATE_FAIL;
      }
      model_set_best_thread((void*)ulii);
    }
    join_threads2(threads, max_thread_ct, is_last_block);
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      // transpose
      read_dptr = &(g_mperm_save_all[pidx]);
      write_dptr = &(chisq_pmajor[pidx * marker_ct]);
      for (marker_bidx = 0; marker_bidx < block_size; marker_bidx++) {
	write_dptr[g_adapt_m_table[marker_bidx]] = read_dptr[marker_bidx * perm_vec_ct];
      }
    }
    marker_idx += block_size;
  } while (marker_idx < marker_unstopped_ct);
  compute_set_scores(marker_ct, perm_vec_ct, set_ct, chisq_pmajor, orig_set_scores, sorted_chisq_buf, sorted_marker_idx_buf, proxy_arr, setdefs, ld_map, apip, chisq_threshold, adaptive_ci_zt, first_adapt_check, perms_done, sip->set_max, perm_adapt_set_unstopped, perm_2success_ct, perm_attempt_ct);
  bigstack_reset(bigstack_mark2);
  if (perms_done < perms_total) {
    if (model_modifier & MODEL_PERM) {
      if (!extract_set_union(setdefs, set_ct, perm_adapt_set_unstopped, unstopped_markers, marker_ct)) {
	perms_done = 0;
	for (set_idx = 0; set_idx < set_ct; set_idx++) {
	  if (perms_done < perm_attempt_ct[set_idx]) {
	    perms_done = perm_attempt_ct[set_idx];
	  }
	}
	goto model_assoc_set_test_perms_done;
      }
      // bugfix (7 Aug 2018): forgot to update marker_unstopped_ct
      marker_unstopped_ct = popcount_longs(unstopped_markers, (marker_ct + BITCT - 1) / BITCT);
    }
    printf("\r%u permutation%s complete.", perms_done, (perms_done != 1)? "s" : "");
    fflush(stdout);
    goto model_assoc_set_test_more_perms;
  }
 model_assoc_set_test_perms_done:
  putc_unlocked('\r', stdout);
  LOGPRINTF("%u permutation%s complete.\n", perms_done, (perms_done != 1)? "s" : "");
 model_assoc_set_test_write:
  if (model_modifier & MODEL_PERM) {
    memcpy(outname_end2, ".set.perm", 10);
  } else {
    memcpy(outname_end2, ".set.mperm", 11);
  }
  retval = write_set_test_results(outname, &(outname_end2[4]), sip, ld_map, setdefs, set_incl, set_ct, marker_ct_orig, marker_ct, marker_idx_to_uidx, marker_ids, max_marker_id_len, perm_2success_ct, perm_attempt_ct, mtest_adjust, perm_count, pfilter, output_min_p, chisq_threshold, orig_chisq, sorted_chisq_buf, sorted_marker_idx_buf, proxy_arr);
  while (0) {
  model_assoc_set_test_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  model_assoc_set_test_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  model_assoc_set_test_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 model_assoc_set_test_ret_1:
  bigstack_reset(bigstack_mark);
  return retval;
}

void get_model_assoc_precomp_bounds(uint32_t missing_ct, uint32_t is_model, uint32_t* minp, uint32_t* ctp) {
  // Estimate which case missing counts are most common.
  // Expected value = (g_perm_case_ct * missing_ct / g_perm_pheno_nm_ct)
  // If X-chromosome and (!is_model):
  //   Lower bound = max(0, missing_ct - 2 * (g_perm_pheno_nm_ct -
  //                 g_perm_case_ct))
  //   Upper bound = min(g_perm_case_ct * 2, missing_ct)
  //   (Could be a bit more precise if we tracked missing male and female
  //    counts separately, but whatever)
  //   Each male automatically contributes 1 to initial missing_ct!
  // Otherwise:
  //   Lower bound = max(0, missing_ct - (g_perm_pheno_nm_ct - g_perm_case_ct))
  //   Upper bound = min(g_perm_case_ct, missing_ct)
  double xval = ((double)(g_perm_case_ct * ((int64_t)missing_ct))) / ((double)((intptr_t)g_perm_pheno_nm_ct));
  intptr_t lbound = (intptr_t)(xval + EPSILON + 1 - ((double)((intptr_t)g_precomp_width)) * 0.5);
  intptr_t ctrl_ct = g_perm_pheno_nm_ct - g_perm_case_ct;
  intptr_t ubound = missing_ct;
  intptr_t lii;
  if (lbound < 0) {
    lbound = 0;
  }
  if (g_is_x && (!is_model)) {
    lii = missing_ct - (2 * ctrl_ct);
    if (((uintptr_t)ubound) > g_perm_case_ct * 2) {
      ubound = g_perm_case_ct * 2;
    }
  } else {
    lii = missing_ct - ctrl_ct;
    if (((uintptr_t)ubound) > g_perm_case_ct) {
      ubound = g_perm_case_ct;
    }
  }
  if (lii > lbound) {
    lbound = lii;
  }
  *minp = lbound;
  if ((intptr_t)(lbound + g_precomp_width) > ubound) {
    *ctp = ubound + 1 - lbound;
  } else {
    *ctp = g_precomp_width;
  }
}

int32_t model_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_male, uint32_t hh_exists, uint32_t ld_ignore_x, uint32_t perm_batch_size, Set_info* sip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  int32_t retval = 0;
  FILE* outfile = nullptr;
  FILE* outfile_msa = nullptr;
  uintptr_t* marker_exclude = marker_exclude_orig;
  uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  uintptr_t marker_ct = marker_ct_orig;
  uintptr_t perm_vec_ct = 0;
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uint32_t model_assoc = model_modifier & MODEL_ASSOC;
  uint32_t model_perms = model_modifier & (MODEL_PERM | MODEL_MPERM);
  uint32_t is_set_test = model_modifier & MODEL_SET_TEST;
  uint32_t model_adapt_nst = (model_modifier & MODEL_PERM) && (!is_set_test);
  uint32_t model_maxt_nst = (model_modifier & MODEL_MPERM) && (!is_set_test);
  uint32_t model_perms_nst = model_perms && (!is_set_test);
  uint32_t model_trendonly = model_modifier & MODEL_TRENDONLY;
  uint32_t model_perm_best = !(model_modifier & MODEL_PMASK);
  uint32_t model_perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t assoc_counts = model_modifier & MODEL_ASSOC_COUNTS;
  uint32_t display_ci = (ci_size > 0);
  uint32_t perms_total = 0;
  uint32_t male_ct = 0;
  uint32_t nonmale_ct = 0;
  uint32_t ctrl_male_ct = 0;
  uint32_t case_male_ct = 0;
  uint32_t ctrl_nonmale_ct = 0;
  uint32_t case_nonmale_ct = 0;
  uint32_t load_ctrl_ct = 0;
  uint32_t load_case_ct = 0;
  uint32_t precomp_width = 0;
  uint32_t is_y = 0;
  int32_t x_code = chrom_info_ptr->xymt_codes[X_OFFSET];
  int32_t y_code = chrom_info_ptr->xymt_codes[Y_OFFSET];
  int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  uintptr_t* sample_nonmale_ctrl_include2 = nullptr;
  uintptr_t* sample_nonmale_case_include2 = nullptr;
  uintptr_t* sample_male_ctrl_include2 = nullptr;
  uintptr_t* sample_male_case_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uintptr_t* cur_ctrl_include2 = nullptr;
  uintptr_t* cur_case_include2 = nullptr;
  uintptr_t* is_invalid_bitfield = nullptr;
  uintptr_t* founder_pnm = nullptr;
  uint32_t* perm_2success_ct = nullptr;
  uint32_t* perm_attempt_ct = nullptr;
  uint32_t* set_cts = nullptr;
  uint32_t* het_cts = nullptr;
  uint32_t* homcom_cts = nullptr;
  uint32_t* precomp_ui = nullptr;
  double* orig_chisq = nullptr;
  double* maxt_extreme_stat = nullptr;
  double* orig_odds = nullptr;
  double* precomp_d = nullptr;
  unsigned char* perm_adapt_stop = nullptr;
  double dxx = 0.0;
  double dww = 0.0;
  double dvv = 0.0;
  double mult_p = 0.0;
  double gen_p = 0.0;
  double dom_p = 0.0;
  double rec_p = 0.0;
  double ca_chisq = 0.0;
  double maxt_cur_extreme_stat = 0;
  uint32_t pct = 0;
  uint32_t max_thread_ct = g_thread_ct;
  uint32_t perm_pass_idx = 0;
  uintptr_t perm_vec_ctcl4m = 0;
  uint32_t model_fisher = model_modifier & MODEL_FISHER;
  uint32_t model_fisherx = model_fisher && (!(model_modifier & MODEL_PTREND));
  uint32_t fisher_midp = model_modifier & MODEL_FISHER_MIDP;
  char* writebuf = g_textbuf;
  char* chrom_name_ptr = nullptr;
  uint32_t chrom_name_len = 0;
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_SLEN];
  uint32_t mu_table[MODEL_BLOCKSIZE];
  uint32_t uibuf[4];
  char wbuf[48];
  char* wptr_start;
  char* wptr;
  char* wptr2;
  char* wptr_mid;
  char* wptr_mid2;
  char* outname_end2;
  uint32_t assoc_thread_ct;
  uint32_t fill_orig_chisq;
  uint32_t marker_unstopped_ct;
  uint32_t gender_req;
  uint32_t case_ct;
  uint32_t ctrl_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t marker_bidx;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t perms_done;
  uintptr_t marker_uidx; // loading
  uintptr_t marker_uidx2; // writing
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uint32_t* marker_idx_to_uidx;
  uint32_t* missp;
  uint32_t* setp;
  uint32_t* hetp;
  uint32_t* missing_cts;
  double* orig_pvals;
  double* orig_pvals_ptr;
  double* ooptr;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_ptr;
  uintptr_t* sample_ctrl_include2;
  uintptr_t* sample_case_include2;
  uint32_t load_sample_ct;
  uintptr_t ulii;
  uint32_t min_ploidy_1;
  uint32_t is_x;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  uint32_t upp;
  uint32_t uqq;
  uint32_t urr;
  uint32_t uss;
  uint32_t is_invalid;
  uint32_t marker_ctl;
  double pval;
  double dyy;
  double dzz;
  double da1;
  double da2;
  double du1;
  double du2;
  double ca_p;
  char* a1ptr;
  char* a2ptr;
  uint32_t loop_end;
  if (pheno_nm_ct < 2) {
    logerrprint("Warning: Skipping --assoc/--model since less than two phenotypes are present.\n");
    goto model_assoc_ret_1;
  }
  if (max_marker_allele_len > MAXLINELEN) {
    if (bigstack_alloc_c(2 * max_marker_allele_len + MAXLINELEN, &writebuf)) {
      goto model_assoc_ret_NOMEM;
    }
  }
  g_model_fisher = model_fisher;
  g_fisher_midp = fisher_midp;
  g_perm_pheno_nm_ct = pheno_nm_ct;
  perms_done = 0;
  g_is_model_prec = model_modifier / MODEL_PREC;
  g_perm_is_1bit = 0;
  g_mperm_save_all = nullptr;
  g_sample_male_include2 = nullptr;
  if (is_set_test) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl, &founder_pnm)) {
      goto model_assoc_ret_NOMEM;
    }
    memcpy(founder_pnm, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
    bitvec_and(founder_info, unfiltered_sample_ctl, founder_pnm);
    if (extract_set_union_unfiltered(sip, nullptr, unfiltered_marker_ct, marker_exclude_orig, &marker_exclude, &marker_ct)) {
      goto model_assoc_ret_NOMEM;
    }
  }
  if (model_maxt_nst) {
    perms_total = model_mperm_val;
    if (bigstack_alloc_d(perms_total, &maxt_extreme_stat)) {
      goto model_assoc_ret_NOMEM;
    }
    g_maxt_extreme_stat = maxt_extreme_stat;
    if (model_fisherx) {
      for (uii = 0; uii < perms_total; uii++) {
	maxt_extreme_stat[uii] = 1;
      }
    } else {
      fill_double_zero(perms_total, maxt_extreme_stat);
    }
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(outname, "w", &outfile_msa)) {
	goto model_assoc_ret_OPEN_FAIL;
      }
      LOGPRINTFWW("Dumping all permutation %svalues to %s .\n", model_fisherx? "p-" : "chi-square ", outname);
    }
  } else {
    mperm_save = 0;
    if (model_adapt_nst) {
      g_aperm_alpha = apip->alpha;
      perms_total = apip->max;
      if (apip->min < apip->init_interval) {
	g_first_adapt_check = (int32_t)(apip->init_interval);
      } else {
	g_first_adapt_check = apip->min;
      }
      g_adaptive_intercept = apip->init_interval;
      g_adaptive_slope = apip->interval_slope;
    }
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw)) {
    goto model_assoc_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  if (model_assoc) {
    if (model_fisher) {
      outname_end2 = memcpyb(outname_end, ".assoc.fisher", 14);
    } else {
      outname_end2 = memcpyb(outname_end, ".assoc", 7);
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(g_logbuf, "Writing C/C --assoc report to %s ... ", outname);
    wordwrapb(25); // strlen("[generating permutations]")
    logprintb();
    fflush(stdout);
    sprintf(g_textbuf, " CHR %%%us         BP   A1 ", plink_maxsnp);
    fprintf(outfile, g_textbuf, "SNP");
    if (assoc_counts) {
      fputs("     C_A      C_U   A2 ", outfile);
    } else {
      fputs("     F_A      F_U   A2 ", outfile);
    }
    if (!model_fisher) {
      fputs("       CHISQ ", outfile);
    }
    if (fputs_checked("           P           OR ", outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (display_ci) {
      uii = (uint32_t)((int32_t)(ci_size * (100 + EPSILON)));
      if (uii >= 10) {
	fprintf(outfile, "          SE          L%u          U%u ", uii, uii);
      } else {
	fprintf(outfile, "          SE           L%u           U%u ", uii, uii);
      }
    }
    if (putc_checked('\n', outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
  } else {
    if (is_set(chrom_info_ptr->haploid_mask, 0)) {
      logerrprint("Error: --model cannot be used on haploid genomes.\n");
      goto model_assoc_ret_INVALID_CMDLINE;
    }
    uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0, 1);
    if (uii) {
      if (is_set_test) {
	// given how the data structures are currently designed, and how easy
	// the command-line fix is, this is not worth the trouble of supporting
	// (this problem illustrates why core data structures should use
	// unfiltered indexes when possible, though)
	logerrprint("Error: --model set-test cannot be used with sets containing MT/haploid\nvariants.  (You can use e.g. '--not-chr y, mt' to exclude them.)\n");
	goto model_assoc_ret_INVALID_CMDLINE;
      }
      LOGPRINTF("Excluding %u MT/haploid variant%s from --model analysis.\n", uii, (uii == 1)? "" : "s");
      marker_ct -= uii;
      if (!marker_ct) {
	logerrprint("Error: No variants remaining for --model analysis.\n");
	goto model_assoc_ret_INVALID_CMDLINE;
      }
    }
    outname_end2 = memcpyb(outname_end, ".model", 7);
    if (fopen_checked(outname, "w", &outfile)) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(g_logbuf, "Writing --model report to %s ... ", outname);
    wordwrapb(25);
    logprintb();
    fflush(stdout);
    if (model_perm_best && model_perms) {
      outname_end2 = memcpyb(outname_end2, ".best", 6);
    } else if ((model_modifier & MODEL_PGEN) && model_perms) {
      outname_end2 = memcpyb(outname_end2, ".gen", 5);
    } else if (model_modifier & MODEL_PDOM) {
      outname_end2 = memcpyb(outname_end2, ".dom", 5);
    } else if (model_modifier & MODEL_PREC) {
      outname_end2 = memcpyb(outname_end2, ".rec", 5);
    } else if (model_modifier & MODEL_PTREND) {
      outname_end2 = memcpyb(outname_end2, ".trend", 7);
    }
    sprintf(g_textbuf, " CHR %%%us   A1   A2     TEST            AFF          UNAFF ", plink_maxsnp);
    fprintf(outfile, g_textbuf, "SNP");
    if (!model_fisher) {
      fputs("       CHISQ   DF ", outfile);
    } else {
      outname_end2 = memcpyb(outname_end2, ".fisher", 8);
    }
    if (fputs_checked("           P\n", outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
  }
  marker_ctl = BITCT_TO_WORDCT(marker_ct);
  g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)marker_ct)));
  if (bigstack_alloc_ul(MODEL_BLOCKSIZE * pheno_nm_ctv2, &loadbuf) ||
      bigstack_alloc_d(marker_ct, &orig_pvals) ||
      bigstack_alloc_ui(marker_ct, &missing_cts)) {
    goto model_assoc_ret_NOMEM;
  }
  g_loadbuf = loadbuf;
  g_orig_pvals = orig_pvals;
  g_missing_cts = missing_cts;
  if (model_assoc) {
    if (bigstack_alloc_d(marker_ct, &orig_odds) ||
        bigstack_alloc_ui(marker_ct, &set_cts)) {
      goto model_assoc_ret_NOMEM;
    }
    g_set_cts = set_cts;
  }
  if ((!model_assoc) || model_maxt_nst) {
    if (bigstack_alloc_ui(marker_ct, &het_cts) ||
        bigstack_alloc_ui(marker_ct, &homcom_cts)) {
      goto model_assoc_ret_NOMEM;
    }
    g_het_cts = het_cts;
    g_homcom_cts = homcom_cts;
  }
  gender_req = ((x_code != -2) && is_set(chrom_info_ptr->chrom_mask, x_code)) || (model_assoc && (((y_code != -2) && is_set(chrom_info_ptr->chrom_mask, y_code))));
  if (gender_req) {
    if (bigstack_alloc_ul(pheno_nm_ctv2, &g_sample_nonmale_include2) ||
	bigstack_alloc_ul(pheno_nm_ctv2, &sample_male_include2)) {
      goto model_assoc_ret_NOMEM;
    }
    g_sample_male_include2 = sample_male_include2;
    quaterarr_collapse_init(sex_male, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, sample_male_include2);
    male_ct = popcount01_longs(sample_male_include2, pheno_nm_ctv2);
    quatervec_01_init_invert(sample_male_include2, pheno_nm_ct, g_sample_nonmale_include2);
    nonmale_ct = pheno_nm_ct - male_ct;
  }
  // Set test does not support Fisher stats, so currently guaranteed to be
  // true there.  Will need to modify this expression if we ever support
  // generation of synthetic chi-square stats from Fisher p-values.
  fill_orig_chisq = (!model_fisherx) || (mtest_adjust && (!model_fisher));
  if (fill_orig_chisq) {
    if (bigstack_calloc_d(marker_ct, &orig_chisq)) {
      goto model_assoc_ret_NOMEM;
    }
  }
  g_orig_chisq = orig_chisq;

  if (model_perms) {
    if (cluster_starts) {
      retval = cluster_include_and_reindex(unfiltered_sample_ct, pheno_nm, 1, pheno_c, pheno_nm_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_perm_cluster_ct, &g_perm_cluster_map, &g_perm_cluster_starts, &g_perm_cluster_case_cts, &g_perm_cluster_cc_preimage);
      if (retval) {
	goto model_assoc_ret_1;
      }
      if (!g_perm_cluster_ct) {
        logerrprint("Error: No size 2+ clusters for permutation test.\n");
	goto model_assoc_ret_INVALID_CMDLINE;
      }
      retval = cluster_alloc_and_populate_magic_nums(g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts, &g_perm_tot_quotients, &g_perm_totq_magics, &g_perm_totq_preshifts, &g_perm_totq_postshifts, &g_perm_totq_incrs);
      if (retval) {
        goto model_assoc_ret_1;
      }
    } else {
      g_perm_cluster_starts = nullptr;
    }
    if (!is_set_test) {
      if (max_thread_ct > perms_total) {
	max_thread_ct = perms_total;
      }
      if (bigstack_init_sfmtp(max_thread_ct)) {
	goto model_assoc_ret_NOMEM;
      }
    }
    if (model_perm_best) {
      if (bigstack_calloc_ul(marker_ctl, &is_invalid_bitfield)) {
	goto model_assoc_ret_NOMEM;
      }
      g_is_invalid_bitfield = is_invalid_bitfield;
    }

    if (!is_set_test) {
      g_ldrefs = (uint16_t*)bigstack_alloc(marker_ct * sizeof(int16_t));
      if (!g_ldrefs) {
	goto model_assoc_ret_NOMEM;
      }
#ifdef __LP64__
      fill_ulong_one((marker_ct + 3) / 4, (uintptr_t*)g_ldrefs);
#else
      fill_ulong_one((marker_ct + 1) / 2, (uintptr_t*)g_ldrefs);
#endif
      if (!(mperm_save & MPERM_DUMP_ALL)) {
	// 5.65686 = roughly 4 * sqrt(2), corresponding to 4 stdevs.  this is
	// a somewhat arbitrary choice.
	// currently just need this to never exceed (2^32 - 1) / (12 * 1024),
	// to avoid uint32_t overflow.
	precomp_width = (1 + (int32_t)(sqrt(pheno_nm_ct) * EXPECTED_MISSING_FREQ * 5.65686));
      } else {
	precomp_width = 0;
      }
      g_precomp_width = precomp_width;
      if (bigstack_calloc_ui(marker_ct, &perm_2success_ct)) {
	goto model_assoc_ret_NOMEM;
      }
      if (model_maxt_nst) {
	if (model_fisherx) {
	  if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC))) {
	    if (bigstack_alloc_ui(precomp_width * 6 * MODEL_BLOCKSIZE, &precomp_ui) ||
		bigstack_alloc_d(precomp_width * 2 * MODEL_BLOCKSIZE, &precomp_d)) {
	      goto model_assoc_ret_NOMEM;
	    }
	  } else if (model_perm_best) {
	    if (bigstack_alloc_ui(precomp_width * 18 * MODEL_BLOCKSIZE, &precomp_ui) ||
		bigstack_alloc_d(precomp_width * 6 * MODEL_BLOCKSIZE, &precomp_d)) {
	      goto model_assoc_ret_NOMEM;
	    }
	  }
	} else if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
	  if (bigstack_alloc_ui(precomp_width * 6 * MODEL_BLOCKSIZE, &precomp_ui) ||
	      bigstack_alloc_d(precomp_width * 2 * MODEL_BLOCKSIZE, &precomp_d)) {
	    goto model_assoc_ret_NOMEM;
	  }
	} else if (model_perm_best) {
	  if (bigstack_alloc_ui(precomp_width * 18 * MODEL_BLOCKSIZE, &precomp_ui) ||
	      bigstack_alloc_d(precomp_width * 6 * MODEL_BLOCKSIZE, &precomp_d)) {
	    goto model_assoc_ret_NOMEM;
	  }
	}
      } else if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
	if (bigstack_alloc_ui(precomp_width * 4 * MODEL_BLOCKSIZE, &precomp_ui)) {
	  goto model_assoc_ret_NOMEM;
	}
      } else if (model_perm_best) {
	if (bigstack_alloc_ui(precomp_width * 12 * MODEL_BLOCKSIZE, &precomp_ui)) {
	  goto model_assoc_ret_NOMEM;
	}
      }
      g_perm_2success_ct = perm_2success_ct;
      if (model_adapt_nst) {
	if (bigstack_alloc_ui(marker_ct, &perm_attempt_ct) ||

	    // we need to zero out trailing bytes of the last word
	    bigstack_calloc_uc(round_up_pow2(marker_ct, BYTECT), &perm_adapt_stop)) {
	  goto model_assoc_ret_NOMEM;
	}
	g_perm_attempt_ct = perm_attempt_ct;
	g_perm_adapt_stop = perm_adapt_stop;
	ujj = apip->max;
	for (uii = 0; uii < marker_ct; uii++) {
	  perm_attempt_ct[uii] = ujj;
	}
      }
    }
    if (!cluster_starts) {
      g_perm_tot_quotient = 0x100000000LLU / pheno_nm_ct;
      magic_num(g_perm_tot_quotient, &g_perm_totq_magic, &g_perm_totq_preshift, &g_perm_totq_postshift, &g_perm_totq_incr);
    }
  }
  g_precomp_ui = precomp_ui;
  g_precomp_d = precomp_d;
  if (bigstack_alloc_ul(pheno_nm_ctv2, &sample_ctrl_include2) ||
      bigstack_alloc_ul(pheno_nm_ctv2, &sample_case_include2)) {
    goto model_assoc_ret_NOMEM;
  }
  quaterarr_collapse_init(pheno_c, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, sample_case_include2);
  case_ct = popcount01_longs(sample_case_include2, pheno_nm_ctv2);
  g_perm_case_ct = case_ct;
  quatervec_01_init_invert(sample_case_include2, pheno_nm_ct, sample_ctrl_include2);
  ctrl_ct = pheno_nm_ct - case_ct;
  if (gender_req) {
    // todo: get rid of these and just use the functions called by the
    // permutation tests
    if (bigstack_alloc_ul(pheno_nm_ctv2, &sample_nonmale_ctrl_include2) ||
	bigstack_alloc_ul(pheno_nm_ctv2, &sample_nonmale_case_include2) ||
	bigstack_alloc_ul(pheno_nm_ctv2, &sample_male_ctrl_include2) ||
	bigstack_alloc_ul(pheno_nm_ctv2, &sample_male_case_include2)) {
      goto model_assoc_ret_NOMEM;
    }
    quaterarr_collapse_init(sex_male, unfiltered_sample_ct, pheno_nm, pheno_nm_ct, sample_male_case_include2);
    bitvec_and(sample_case_include2, pheno_nm_ctv2, sample_male_case_include2);
    case_male_ct = popcount01_longs(sample_male_case_include2, pheno_nm_ctv2);
    bitvec_andnot_copy(sample_male_include2, sample_male_case_include2, pheno_nm_ctv2, sample_male_ctrl_include2);
    bitvec_andnot_copy(sample_case_include2, sample_male_case_include2, pheno_nm_ctv2, sample_nonmale_case_include2);
    bitvec_andnot_copy(sample_ctrl_include2, sample_male_ctrl_include2, pheno_nm_ctv2, sample_nonmale_ctrl_include2);
    ctrl_male_ct = male_ct - case_male_ct;
    case_nonmale_ct = case_ct - case_male_ct;
    ctrl_nonmale_ct = ctrl_ct - ctrl_male_ct;
  }

  for (uii = 1; uii <= MODEL_BLOCKSIZE; uii++) {
    loadbuf[uii * pheno_nm_ctv2 - 2] = 0;
    loadbuf[uii * pheno_nm_ctv2 - 1] = 0;
  }
  if (model_perms) {
    if (bigstack_left() < pheno_nm_ctv2 * sizeof(intptr_t)) {
      goto model_assoc_ret_NOMEM;
    }
  }
  marker_unstopped_ct = marker_ct;

  // ----- begin main loop -----
 model_assoc_more_perms:
  if (model_perms_nst) {
    if (!perm_pass_idx) {
      fputs("[generating permutations]", stdout);
      fflush(stdout);
    }
    if (model_adapt_nst) {
      if (perm_pass_idx) {
	uii = g_first_adapt_check;
	ujj = apip->init_interval;
	while (uii <= perms_done) {
	  // APERM_MAX prevents infinite loop here
	  uii += (int32_t)(ujj + ((int32_t)uii) * apip->interval_slope);
	}
	g_first_adapt_check = uii;
      }
      perm_vec_ct = bigstack_left() / (pheno_nm_ctv2 * sizeof(intptr_t));
    } else {
      // perm_vec_ct memory allocation dependencies:
      //   g_maxt_thread_results: (8 * perm_vec_ct, cacheline-aligned) *
      //     max_thread_ct
      //   g_perm_vecst: 16 * ((perm_vec_ct + 127) / 128) * pheno_nm_ct
      //   g_thread_git_wkspace: ((perm_vec_ct + 127) / 128) * 1152 * thread_ct
      //   g_resultbuf: MODEL_BLOCKSIZE * (4 * perm_vec_ct, CL-aligned) * 3
      //   g_perm_vecs: pheno_nm_ctv2 * sizeof(intptr_t) * perm_vec_ct
      //   g_mperm_save_all (if needed): marker_ct * 8 * perm_vec_ct
      // If we force perm_vec_ct to be a multiple of 128, then we have
      //   perm_vec_ct * (17 * max_thread_ct + 12 * MODEL_BLOCKSIZE +
      //                    pheno_nm_ct / 8 + sizeof(intptr_t) * pheno_nm_ctv2
      //                    [+ marker_ct * sizeof(double) * mperm_save_all])
      //
      // Each max(T) thread has six buffers to support rapid execution of the
      // genotype indexing and LD exploiter algorithms:
      //   six with 4-bit accumulators, each has size perm_vec_ct / 2 bytes
      //   six with 8-bit accumulators, each has size perm_vec_ct bytes
      // The initial 6 multiplier is to allow heterozygote, homozygote minor,
      // and missing genotype increments and decrements to be counted
      // simultaneously.
      // Adding all this up, we have 9 * perm_vec_ct bytes, and multiplying
      // by 128 yields 1152.  The other thread_ct dependence contributes
      // 8 * perm_vec_ct bytes, multiplying by 128 yields 1024, and
      // 1152 + 1024 = 2176.
      if (mperm_save & MPERM_DUMP_ALL) {
        perm_vec_ct = 128 * (bigstack_left() / (128LL * sizeof(intptr_t) * pheno_nm_ctv2 + 2176LL * max_thread_ct + 1536LL * MODEL_BLOCKSIZE + 16LL * pheno_nm_ct + 128LL * sizeof(double) * marker_ct));
      } else {
        perm_vec_ct = 128 * (bigstack_left() / (128LL * sizeof(intptr_t) * pheno_nm_ctv2 + 2176LL * max_thread_ct + 1536LL * MODEL_BLOCKSIZE + 16LL * pheno_nm_ct));
      }
    }
    if (perm_vec_ct > perms_total - perms_done) {
      perm_vec_ct = perms_total - perms_done;
    } else if (!perm_vec_ct) {
      goto model_assoc_ret_NOMEM;
    }
    perm_vec_ctcl4m = round_up_pow2(perm_vec_ct, CACHELINE_INT32);
    perms_done += perm_vec_ct;
    g_perms_done = perms_done;
    g_perm_vec_ct = perm_vec_ct;
    bigstack_alloc_ul(perm_vec_ct * pheno_nm_ctv2, &g_perm_vecs);
    g_perm_generation_thread_ct = MINV(max_thread_ct, perm_vec_ct);
    ulii = 0;
    if (!cluster_starts) {
      if (spawn_threads(threads, &generate_cc_perms_thread, g_perm_generation_thread_ct)) {
	goto model_assoc_ret_THREAD_CREATE_FAIL;
      }
      generate_cc_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &generate_cc_cluster_perms_thread, g_perm_generation_thread_ct)) {
	goto model_assoc_ret_THREAD_CREATE_FAIL;
      }
      generate_cc_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_perm_generation_thread_ct);
    g_assoc_thread_ct = max_thread_ct;
    if (!model_adapt_nst) {
      bigstack_alloc_d(max_thread_ct * round_up_pow2(perm_vec_ct, CACHELINE_DBL), &g_maxt_thread_results);
      bigstack_alloc_ui(perm_vec_ctcl4m * 3 * MODEL_BLOCKSIZE, &g_resultbuf);
#ifdef __LP64__
      ulii = ((perm_vec_ct + 127) / 128) * 4;
      bigstack_alloc_ui(ulii * pheno_nm_ct, &g_perm_vecst);
#else
      ulii = (perm_vec_ct + 31) / 32;
      bigstack_alloc_ui(ulii * pheno_nm_ct, &g_perm_vecst);
      ulii = ((perm_vec_ct + 63) / 64) * 2;
#endif
      bigstack_calloc_ui(ulii * 72 * max_thread_ct, &g_thread_git_wkspace);
      transpose_perms(g_perm_vecs, perm_vec_ct, pheno_nm_ct, g_perm_vecst);
      if (mperm_save & MPERM_DUMP_ALL) {
	bigstack_alloc_d(marker_ct * perm_vec_ct, &g_mperm_save_all);
      }
    }
    if (!perm_pass_idx) {
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                         \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", stdout);
    }
  }
  if (!perm_pass_idx) {
    fputs("0%", stdout);
    fflush(stdout);
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto model_assoc_ret_READ_FAIL;
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  loop_end = marker_ct / 100;
  do {
    if (marker_uidx >= chrom_end) {
      g_block_start = 0;
      if (model_assoc) {
	// exploit overflow
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
	min_ploidy_1 |= uii; // treat MT as haploid
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if (min_ploidy_1 && (!is_x)) {
	  if (is_y) {
	    cur_ctrl_include2 = sample_male_ctrl_include2;
	    cur_case_include2 = sample_male_case_include2;
	    load_sample_ct = male_ct;
	    load_case_ct = case_male_ct;
	  } else {
	    cur_ctrl_include2 = sample_ctrl_include2;
	    cur_case_include2 = sample_case_include2;
	    load_sample_ct = pheno_nm_ct;
	    load_case_ct = case_ct;
	  }
	  load_ctrl_ct = load_sample_ct - load_case_ct;
	}
	g_min_ploidy_1 = min_ploidy_1;
	g_is_y = is_y;
      } else {
	while (1) {
	  do {
	    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
	  } while (marker_uidx >= chrom_end);
	  uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_x = (uii == (uint32_t)x_code);
	  if (((!IS_SET(haploid_mask, uii)) && (uii != (uint32_t)mt_code)) || is_x) {
	    break;
	  }
	  marker_uidx = next_unset_unsafe(marker_exclude, chrom_end);
	}
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
      }
      g_is_x = is_x;
      chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, uii, &chrom_name_len, chrom_name_buf);
    } else if (model_maxt_nst) {
      marker_idx -= MODEL_BLOCKKEEP;
      if (marker_idx) { // max(T) initial block special case, see below
        memcpy(loadbuf, &(loadbuf[(MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * pheno_nm_ctv2]), MODEL_BLOCKKEEP * pheno_nm_ctv2 * sizeof(intptr_t));
        memcpy(g_resultbuf, &(g_resultbuf[3 * (MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * perm_vec_ctcl4m]), MODEL_BLOCKKEEP * perm_vec_ctcl4m * 3 * sizeof(int32_t));
      }
      g_block_start = MODEL_BLOCKKEEP;
    } else {
      g_block_start = 0;
    }
    block_size = g_block_start;
    block_end = marker_unstopped_ct - marker_idx;
    if ((!marker_idx) && (!block_size) && model_maxt_nst) {
      // For max(T) permutation tests, minimize how long we have to work with
      // crappy precomputed values.  Most important when using Fisher exact
      // test p-values.
      if (block_end > MODEL_BLOCKKEEP) {
	block_end = MODEL_BLOCKKEEP;
      }
    } else if (block_end > MODEL_BLOCKSIZE) {
      block_end = MODEL_BLOCKSIZE;
    }
    do {
      if (model_adapt_nst && perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx++;
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(loadbuf[block_size * pheno_nm_ctv2]);
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf_ptr)) {
	goto model_assoc_ret_READ_FAIL;
      }
      if (model_adapt_nst) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      if (is_x && (!model_assoc)) {
	force_missing((unsigned char*)(&(loadbuf[block_size * pheno_nm_ctv2])), sample_male_include2, pheno_nm_ct);
      }
      // no need for usual haploid_fix since the popcount routines here
      // interpret het. haploids as missing anyway
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (block_size == g_block_start) {
      continue;
    }
    if (!perm_pass_idx) {
      // basic --assoc/--model
      orig_pvals_ptr = &(orig_pvals[marker_idx + g_block_start]);
      missp = &(missing_cts[marker_idx + g_block_start]);
      if (model_assoc) {
	setp = &(set_cts[marker_idx + g_block_start]);
	ooptr = &(orig_odds[marker_idx + g_block_start]);
	for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	  marker_uidx2 = mu_table[marker_bidx];
	  if (!min_ploidy_1) {
	    if (model_maxt_nst) {
	      single_marker_cc_3freqs(pheno_nm_ctv2, &(loadbuf[marker_bidx * pheno_nm_ctv2]), sample_ctrl_include2, sample_case_include2, &unn, &uoo, &ujj, &upp, &uqq, &umm);
	      het_cts[marker_idx + marker_bidx] = uoo + uqq;
	      homcom_cts[marker_idx + marker_bidx] = unn + upp;
	      uii = 2 * unn + uoo;
	      ukk = 2 * upp + uqq;
	    } else {
	      single_marker_cc_freqs(pheno_nm_ctv2, &(loadbuf[marker_bidx * pheno_nm_ctv2]), sample_ctrl_include2, sample_case_include2, &uii, &ujj, &ukk, &umm);
	    }
	    *missp = ujj + umm;
	    *setp = uii + ukk;
	    ujj = 2 * (ctrl_ct - ujj) - uii;
	    umm = 2 * (case_ct - umm) - ukk;
	  } else if (is_x) {
	    single_marker_cc_freqs(pheno_nm_ctv2, &(loadbuf[marker_bidx * pheno_nm_ctv2]), sample_nonmale_ctrl_include2, sample_nonmale_case_include2, &uii, &ujj, &ukk, &umm);
	    *missp = 2 * (ujj + umm);
	    *setp = uii + ukk;
	    ujj = 2 * (ctrl_nonmale_ct - ujj) - uii;
	    umm = 2 * (case_nonmale_ct - umm) - ukk;
	    haploid_single_marker_cc_freqs(pheno_nm_ctv2, &(loadbuf[marker_bidx * pheno_nm_ctv2]), sample_male_ctrl_include2, sample_male_case_include2, &unn, &uoo, &upp, &uqq);
	    *missp += uoo + uqq + male_ct;
	    *setp += unn + upp;
	    uoo = ctrl_male_ct - uoo - unn;
	    uqq = case_male_ct - uqq - upp;
	    uii += unn;
	    ujj += uoo;
	    ukk += upp;
	    umm += uqq;
	  } else {
	    haploid_single_marker_cc_freqs(pheno_nm_ctv2, &(loadbuf[marker_bidx * pheno_nm_ctv2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm);
	    *missp = ujj + umm;
	    *setp = uii + ukk;
	    ujj = load_ctrl_ct - ujj - uii;
	    umm = load_case_ct - umm - ukk;
	    if (is_y) {
	      *missp += nonmale_ct;
	    } else if (model_maxt_nst) {
	      het_cts[marker_idx + marker_bidx] = 0;
	      homcom_cts[marker_idx + marker_bidx] = *setp;
	    }
	  }
	  da1 = umm;
	  da2 = ukk;
	  du1 = ujj;
	  du2 = uii;
	  if (model_fisher) {
            // bugfix (12 Jun 2018): If MAF is zero, test should not be
            // considered valid for --adjust or permutation testing purposes.
            // plink 1.07 got this right, but in a wrong way: it considered
            // *all* Fisher's-exact-test p-values of 1 to be invalid tests.  So
            // we don't generally want to match its output (even before
            // considering the problems with its fisher22 routine).
            if ((umm + ujj) && (ukk + uii)) {
              pval = fisher22(uii, ujj, ukk, umm, fisher_midp);
            } else {
              pval = -9;
            }
	    *orig_pvals_ptr = pval;
	  } else {
	    if ((umm + ujj) && (ukk + uii)) {
	      dxx = chi22_eval(ukk, ukk + umm, uii + ukk, uii + ujj + ukk + umm);
	      pval = chiprob_p(dxx, 1);
	      *orig_pvals_ptr = pval;
	      if (fill_orig_chisq) {
		orig_chisq[marker_idx + marker_bidx] = dxx;
	      }
	    } else {
	      *orig_pvals_ptr = -9;
	      pval = -1;
	      dxx = 0;
	      if (fill_orig_chisq) {
		orig_chisq[marker_idx + marker_bidx] = -9;
	      }
            }
	  }
	  *ooptr = (da1 * du2) / (du1 * da2);
	  if ((pfilter == 2.0) || ((pval <= pfilter) && (pval >= 0.0))) {
	    a1ptr = marker_allele_ptrs[2 * marker_uidx2];
	    a2ptr = marker_allele_ptrs[2 * marker_uidx2 + 1];
	    wptr = memcpyax(writebuf, chrom_name_ptr, chrom_name_len, ' ');
	    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	    *wptr++ = ' ';
	    wptr = uint32toa_w10x(marker_pos[marker_uidx2], ' ', wptr);
	    wptr = fw_strcpy(4, a1ptr, wptr);
	    *wptr++ = ' ';
	    if (umm + ukk) {
	      if (assoc_counts) {
		wptr = uint32toa_w8(umm, wptr);
	      } else {
		wptr = dtoa_g_wxp4(da1 / (da1 + da2), 8, wptr);
	      }
	      *wptr++ = ' ';
	    } else {
	      wptr = memcpya(wptr, "      NA ", 9);
	    }
	    if (ujj + uii) {
	      if (assoc_counts) {
		wptr = uint32toa_w8(ujj, wptr);
	      } else {
		wptr = dtoa_g_wxp4(du1 / (du1 + du2), 8, wptr);
	      }
	    } else {
	      wptr = memcpya(wptr, "      NA", 8);
	    }
	    *wptr = ' ';
	    wptr = fw_strcpy(4, a2ptr, &(wptr[1]));
	    *wptr++ = ' ';
	    if (model_fisher) {
              if (pval == -9) {
                wptr = memcpya(wptr, "           1", 12);
              } else {
                wptr = dtoa_g_wxp4(MAXV(pval, output_min_p), 12, wptr);
              }
	    } else {
	      if (pval > -1) {
		wptr = dtoa_g_wxp4x(dxx, 12, ' ', wptr);
		wptr = dtoa_g_wxp4(MAXV(pval, output_min_p), 12, wptr);
	      } else {
		wptr = memcpya(wptr, "          NA           NA", 25);
	      }
	    }
	    *wptr++ = ' ';
	    if (du1 * da2 == 0.0) {
	      wptr = memcpya(wptr, "          NA", 12);
	      if (display_ci) {
		wptr = memcpya(wptr, "           NA           NA           NA", 39);
	      }
	    } else {
	      wptr = dtoa_g_wxp4(*ooptr, 12, wptr);
	      if (display_ci) {
		dxx = log(*ooptr);
		dyy = sqrt(1 / da1 + 1 / da2 + 1 / du1 + 1 / du2);
		dzz = ci_zt * dyy;
		dww = exp(dxx - dzz);
		dvv = exp(dxx + dzz);
		*wptr++ = ' ';
		wptr = dtoa_g_wxp4x(dyy, 12, ' ', wptr);
		wptr = dtoa_g_wxp4x(dww, 12, ' ', wptr);
		wptr = dtoa_g_wxp4(dvv, 12, wptr);
	      }
	    }
	    wptr = memcpya(wptr, " \n", 2);
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	  }
	  missp++;
	  setp++;
	  orig_pvals_ptr++;
	  ooptr++;
	}
      } else {
	// repurpose setp as homcom_cts pointer
	setp = &(homcom_cts[marker_idx + g_block_start]);
	hetp = &(het_cts[marker_idx + g_block_start]);
	for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	  marker_uidx2 = mu_table[marker_bidx];
	  single_marker_cc_3freqs(pheno_nm_ctv2, &(loadbuf[marker_bidx * pheno_nm_ctv2]), sample_ctrl_include2, sample_case_include2, &uii, &ujj, &ukk, &umm, &unn, &uoo);
	  *missp = ukk + uoo;
	  *setp = uii + umm;
	  ukk = pheno_nm_ct - case_ct - uii - ujj - ukk;
	  uoo = case_ct - umm - unn - uoo;
	  *hetp = ujj + unn;
	  is_invalid = (uoo < model_cell_ct) || (unn < model_cell_ct) || (umm < model_cell_ct) || (ukk < model_cell_ct) || (ujj < model_cell_ct) || (uii < model_cell_ct);
	  a1ptr = marker_allele_ptrs[2 * marker_uidx2];
	  a2ptr = marker_allele_ptrs[2 * marker_uidx2 + 1];
	  wptr = memcpya(writebuf, chrom_name_ptr, chrom_name_len);
	  *wptr++ = ' ';
	  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	  *wptr++ = ' ';
	  wptr = fw_strcpy(4, a1ptr, wptr);
	  *wptr++ = ' ';
	  wptr = fw_strcpy(4, a2ptr, wptr);
	  memset(wptr, 32, 2);
	  wptr = &(wptr[2]);
	  wptr_mid = wptr;
	  if (!model_trendonly) {
	    memcpy(wptr, "   GENO ", 8);
	    wptr2 = uint32toa_x(uoo, '/', wbuf);
	    wptr2 = uint32toa_x(unn, '/', wptr2);
	    wptr2 = uint32toa(umm, wptr2);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr[8]));
	    *wptr++ = ' ';
	    wptr2 = uint32toa_x(ukk, '/', wbuf);
	    wptr2 = uint32toa_x(ujj, '/', wptr2);
	    wptr2 = uint32toa(uii, wptr2);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	    *wptr++ = ' ';
	    if (is_invalid) {
	      gen_p = -9;
	      if (fill_orig_chisq && (model_modifier & MODEL_PGEN)) {
		orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      if (model_fisher) {
		gen_p = fisher23(uii, ujj, ukk, umm, unn, uoo, fisher_midp);
	      } else {
		chi23_evalx(uii, ujj, ukk, umm, unn, uoo, &dvv, &upp);
		gen_p = chiprob_px(dvv, upp);
		if (fill_orig_chisq && (model_modifier & MODEL_PGEN)) {
		  if (dvv != -9) {
		    orig_chisq[marker_idx + marker_bidx] = dvv;
		  } else {
		    orig_chisq[marker_idx + marker_bidx] = 0;
		  }
		}
	      }
	    }
	    if (gen_p < -1) {
	      wptr = model_assoc_tna(model_fisher, wptr);
	    } else {
	      if (!model_fisher) {
		wptr = dtoa_g_wxp4(dvv, 12, wptr);
		wptr = memcpya(wptr, "    ", 4);
		*wptr++ = '0' + upp;
		*wptr++ = ' ';
	      }
	      wptr = dtoa_g_wxp4x(MAXV(gen_p, output_min_p), 12, '\n', wptr);
	    }
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	  }
	  memcpy(wptr_mid, "  TREND ", 8);
	  wptr2 = uint32toa_x(uoo * 2 + unn, '/', wbuf);
	  wptr2 = uint32toa(umm * 2 + unn, wptr2);
	  wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr_mid[8]));
	  *wptr++ = ' ';
	  wptr2 = uint32toa_x(ukk * 2 + ujj, '/', wbuf);
	  wptr2 = uint32toa(uii * 2 + ujj, wptr2);
	  wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	  *wptr++ = ' ';
	  wptr_mid2 = wptr; // save this for next line
	  ca_chisq = ca_trend_evalx(umm * 2 + unn, umm + unn + uoo, ujj + unn, uii + umm, uii + ujj + ukk + umm + unn + uoo);
	  ca_p = chiprob_px(ca_chisq, 1);
	  if (fill_orig_chisq && (model_modifier & MODEL_PTREND)) {
	    if (ca_chisq != -9) {
	      orig_chisq[marker_idx + marker_bidx] = ca_chisq;
	    } else {
	      orig_chisq[marker_idx + marker_bidx] = 0;
	    }
	  }
	  if (ca_p > -1) {
	    if (!model_fisher) {
	      wptr = dtoa_g_wxp4(ca_chisq, 12, wptr);
	      wptr = memcpya(wptr, "    1 ", 6);
	    }
	    wptr = dtoa_g_wxp4x(MAXV(ca_p, output_min_p), 12, '\n', wptr);
	  } else {
	    wptr = model_assoc_tna(model_fisher, wptr);
	  }
	  if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	  if (!model_trendonly) {
	    memcpy(wptr_mid, "ALLELIC", 7);
	    wptr = wptr_mid2;
	    if (model_fisher) {
	      mult_p = fisher22(2 * uoo + unn, 2 * umm + unn, 2 * ukk + ujj, 2 * uii + ujj, fisher_midp);
	    } else {
	      dww = chi22_evalx(2 * uoo + unn, 2 * (uoo + unn + umm), 2 * (uoo + ukk) + unn + ujj, 2 * (uoo + unn + umm + ukk + ujj + uii));
	      mult_p = chiprob_px(dww, 1);
	    }
	    if (mult_p > -1) {
	      if (!model_fisher) {
		wptr = dtoa_g_wxp4(dww, 12, wptr);
		wptr = memcpya(wptr, "    1 ", 6);
	      }
	      wptr = dtoa_g_wxp4x(MAXV(mult_p, output_min_p), 12, '\n', wptr);
	    } else {
	      wptr = model_assoc_tna(model_fisher, wptr);
	    }
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    memcpy(wptr_mid, "    DOM", 7);
	    wptr2 = uint32toa_x(uoo + unn, '/', wbuf);
	    wptr2 = uint32toa(umm, wptr2);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr_mid[8]));
	    *wptr++ = ' ';
	    wptr2 = uint32toa_x(ukk + ujj, '/', wbuf);
	    wptr2 = uint32toa(uii, wptr2);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	    *wptr++ = ' ';
	    if (is_invalid) {
	      dom_p = -9;
	      if (fill_orig_chisq && (model_modifier & MODEL_PDOM)) {
		orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      if (model_fisher) {
		dom_p = fisher22(uoo + unn, umm, ukk + ujj, uii, fisher_midp);
	      } else {
		dww = chi22_evalx(uoo + unn, uoo + unn + umm, uoo + unn + ukk + ujj, uoo + unn + umm + ukk + ujj + uii);
		dom_p = chiprob_px(dww, 1);
		if (fill_orig_chisq && (model_modifier & MODEL_PDOM)) {
		  if (dww != -9) {
		    orig_chisq[marker_idx + marker_bidx] = dww;
		  } else {
		    orig_chisq[marker_idx + marker_bidx] = 0;
		  }
		}
	      }
	    }
	    if (dom_p < -1) {
	      wptr = model_assoc_tna(model_fisher, wptr);
	    } else {
	      if (!model_fisher) {
		wptr = dtoa_g_wxp4(dww, 12, wptr);
		wptr = memcpya(wptr, "    1 ", 6);
	      }
	      wptr = dtoa_g_wxp4x(MAXV(dom_p, output_min_p), 12, '\n', wptr);
	    }
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    memcpy(&(wptr_mid[4]), "REC", 3);
	    wptr2 = uint32toa_x(uoo, '/', wbuf);
	    wptr2 = uint32toa(unn + umm, wptr2);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr_mid[8]));
	    *wptr++ = ' ';
	    wptr2 = uint32toa_x(ukk, '/', wbuf);
	    wptr2 = uint32toa(ujj + uii, wptr2);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	    *wptr++ = ' ';
	    if (is_invalid) {
	      rec_p = -9;
	      if (fill_orig_chisq && (model_modifier & MODEL_PREC)) {
		orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      if (model_fisher) {
		rec_p = fisher22(uoo, unn + umm, ukk, ujj + uii, fisher_midp);
	      } else {
		dww = chi22_evalx(uoo, uoo + unn + umm, uoo + ukk, uoo + unn + umm + ukk + ujj + uii);
		rec_p = chiprob_px(dww, 1);
		if (fill_orig_chisq && (model_modifier & MODEL_PREC)) {
		  if (dww != -9) {
		    orig_chisq[marker_idx + marker_bidx] = dww;
		  } else {
		    orig_chisq[marker_idx + marker_bidx] = 0;
		  }
		}
	      }
	    }
	    if (rec_p < -1) {
	      wptr = model_assoc_tna(model_fisher, wptr);
	    } else {
	      if (!model_fisher) {
		wptr = dtoa_g_wxp4(dww, 12, wptr);
		wptr = memcpya(wptr, "    1 ", 6);
	      }
	      wptr = dtoa_g_wxp4x(MAXV(rec_p, output_min_p), 12, '\n', wptr);
	    }
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	  }
	  if (model_perm_best) {
	    dxx = mult_p;
	    if (!is_invalid) {
	      if ((dom_p < dxx) && (dom_p >= 0)) {
		dxx = dom_p;
	      }
	      if ((rec_p < dxx) && (rec_p >= 0)) {
		dxx = rec_p;
	      }
	    }
	    if (model_perms && is_invalid) {
	      set_bit_ul(marker_idx + marker_bidx, is_invalid_bitfield);
	    }
	    if (fill_orig_chisq) {
	      if (dxx != -9) {
		orig_chisq[marker_idx + marker_bidx] = inverse_chiprob(dxx, 1);
	      } else {
		orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    }
	  } else if (model_modifier & MODEL_PGEN) {
	    dxx = (gen_p >= 0)? gen_p : -9;
	  } else if (model_modifier & MODEL_PDOM) {
	    dxx = (dom_p >= 0)? dom_p : -9;
	  } else if (model_modifier & MODEL_PREC) {
	    dxx = (rec_p >= 0)? rec_p : -9;
	  } else if (model_modifier & MODEL_PTREND) {
	    dxx = (ca_p >= 0)? ca_p : -9;
	  }
	  missp++;
	  setp++;
	  hetp++;
	  *orig_pvals_ptr++ = dxx;
	}
      }
    }
    if (model_perms_nst) {
      g_block_diff = block_size - g_block_start;
      assoc_thread_ct = g_block_diff;
      if (assoc_thread_ct > max_thread_ct) {
	assoc_thread_ct = max_thread_ct;
      }
      if (model_maxt_nst) {
	if (model_fisherx) {
	  maxt_cur_extreme_stat = maxt_extreme_stat[0];
	  for (uii = 1; uii < perm_vec_ct; uii++) {
	    dxx = maxt_extreme_stat[uii];
	    if (dxx > maxt_cur_extreme_stat) {
	      maxt_cur_extreme_stat = dxx;
	    }
	  }
	} else {
	  maxt_cur_extreme_stat = maxt_extreme_stat[0];
	  for (uii = 1; uii < perm_vec_ct; uii++) {
	    dxx = maxt_extreme_stat[uii];
	    if (dxx < maxt_cur_extreme_stat) {
	      maxt_cur_extreme_stat = dxx;
	    }
	  }
	}
      }
      if (model_assoc) {
	if (min_ploidy_1) {
	  uqq = 1;
	} else {
	  uqq = 2;
	}
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt_nst) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 0, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  uoo = set_cts[urr];
	  if (is_x) {
	    unn = 2 * case_ct;
	    upp = 2 * pheno_nm_ct - upp;
	  } else {
	    unn = uqq * case_ct;
	    upp = uqq * (pheno_nm_ct - upp);
	  }
	  ujj *= uqq;
	  ukk += uii * precomp_width;
	  if (model_fisher) {
	    dxx = orig_pvals[urr];
	    if (model_adapt_nst) {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		fisher22_precomp_pval_bounds(dxx, fisher_midp, unn - ujj, uoo, upp, &(precomp_ui[umm * 4]), nullptr);
		ujj += uqq;
	      }
	    } else {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		fisher22_precomp_pval_bounds(dxx, fisher_midp, unn - ujj, uoo, upp, &(precomp_ui[umm * 6]), nullptr);
		fisher22_precomp_pval_bounds(maxt_cur_extreme_stat, fisher_midp, unn - ujj, uoo, upp, uibuf, &(precomp_d[umm * 2]));
		precomp_ui[umm * 6 + 4] = uibuf[2];
		precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
		ujj += uqq;
	      }
	    }
	  } else {
	    dxx = orig_chisq[urr];
	    if (model_adapt_nst) {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - ujj, uoo, upp, &(precomp_ui[umm * 4]), nullptr);
		ujj += uqq;
	      }
	    } else {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - ujj, uoo, upp, &(precomp_ui[umm * 6]), nullptr);
		chi22_precomp_val_bounds(maxt_cur_extreme_stat, unn - ujj, uoo, upp, uibuf, &(precomp_d[umm * 2]));
		precomp_ui[umm * 6 + 4] = uibuf[2];
		precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
		ujj += uqq;
	      }
	    }
	  }
	}
      } else if (model_perm_best) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt_nst) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 1, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  unn = 2 * case_ct;
	  uqq = 2 * (pheno_nm_ct - upp);
	  uoo = 2 * homcom_cts[urr] + het_cts[urr];
	  ukk += uii * precomp_width;
	  uss = 2 * ujj;
	  if (model_fisher) {
	    dxx = orig_pvals[urr];
	    if (model_adapt_nst) {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, fisher_midp, unn - uss, uoo, uqq, &(precomp_ui[umm * 12]), nullptr);
		uss += 2;
	      }
	    } else {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, fisher_midp, unn - uss, uoo, uqq, &(precomp_ui[umm * 18]), nullptr);
	        fisher22_precomp_pval_bounds(maxt_cur_extreme_stat, fisher_midp, 2 * case_ct - uss, uoo, uqq, uibuf, &(precomp_d[umm * 6]));
		precomp_ui[umm * 18 + 4] = uibuf[2];
		precomp_ui[umm * 18 + 5] = uibuf[3] - uibuf[2];
		uss += 2;
	      }
	    }
	    if (!IS_SET(is_invalid_bitfield, urr)) {
	      upp = pheno_nm_ct - upp;
	      uoo = homcom_cts[urr];
	      uqq = upp - uoo - het_cts[urr];
	      ujj = case_ct - ujj;
	      if (model_adapt_nst) {
		for (umm = uii * precomp_width; umm < ukk; umm++) {
		  fisher22_precomp_pval_bounds(dxx, fisher_midp, ujj, uoo, upp, &(precomp_ui[umm * 12 + 4]), nullptr);
		  fisher22_precomp_pval_bounds(dxx, fisher_midp, ujj, uqq, upp, &(precomp_ui[umm * 12 + 8]), nullptr);
		  ujj--;
		}
	      } else {
		for (umm = uii * precomp_width; umm < ukk; umm++) {
		  fisher22_precomp_pval_bounds(dxx, fisher_midp, ujj, uoo, upp, &(precomp_ui[umm * 18 + 6]), nullptr);
		  fisher22_precomp_pval_bounds(maxt_cur_extreme_stat, fisher_midp, ujj, uoo, upp, uibuf, &(precomp_d[umm * 6 + 2]));
		  precomp_ui[umm * 18 + 10] = uibuf[2];
		  precomp_ui[umm * 18 + 11] = uibuf[3] - uibuf[2];
		  fisher22_precomp_pval_bounds(dxx, fisher_midp, ujj, uqq, upp, &(precomp_ui[umm * 18 + 12]), nullptr);
		  fisher22_precomp_pval_bounds(maxt_cur_extreme_stat, fisher_midp, ujj, uqq, upp, uibuf, &(precomp_d[umm * 6 + 4]));
		  precomp_ui[umm * 18 + 16] = uibuf[2];
		  precomp_ui[umm * 18 + 17] = uibuf[3] - uibuf[2];
		  ujj--;
		}
	      }
	    }
	  } else {
	    dxx = orig_chisq[urr];
	    if (model_adapt_nst) {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - uss, uoo, uqq, &(precomp_ui[umm * 12]), nullptr);
		uss += 2;
	      }
	    } else {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - uss, uoo, uqq, &(precomp_ui[umm * 18]), nullptr);
		chi22_precomp_val_bounds(maxt_cur_extreme_stat, unn - uss, uoo, uqq, uibuf, &(precomp_d[umm * 6]));
		precomp_ui[umm * 18 + 4] = uibuf[2];
		precomp_ui[umm * 18 + 5] = uibuf[3] - uibuf[2];
		uss += 2;
	      }
	    }
	    if (!IS_SET(is_invalid_bitfield, urr)) {
	      upp = pheno_nm_ct - upp;
	      uoo = homcom_cts[urr];
	      uqq = upp - uoo - het_cts[urr];
	      ujj = case_ct - ujj;
	      if (model_adapt_nst) {
		for (umm = uii * precomp_width; umm < ukk; umm++) {
		  chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(precomp_ui[umm * 12 + 4]), nullptr);
		  chi22_precomp_val_bounds(dxx, ujj, uqq, upp, &(precomp_ui[umm * 12 + 8]), nullptr);
		  ujj--;
		}
	      } else {
		for (umm = uii * precomp_width; umm < ukk; umm++) {
		  chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(precomp_ui[umm * 18 + 6]), nullptr);
		  chi22_precomp_val_bounds(maxt_cur_extreme_stat, ujj, uoo, upp, uibuf, &(precomp_d[umm * 6 + 2]));
		  precomp_ui[umm * 18 + 10] = uibuf[2];
		  precomp_ui[umm * 18 + 11] = uibuf[3] - uibuf[2];
		  chi22_precomp_val_bounds(dxx, ujj, uqq, upp, &(precomp_ui[umm * 18 + 12]), nullptr);
		  chi22_precomp_val_bounds(maxt_cur_extreme_stat, ujj, uqq, upp, uibuf, &(precomp_d[umm * 6 + 4]));
		  precomp_ui[umm * 18 + 16] = uibuf[2];
		  precomp_ui[umm * 18 + 17] = uibuf[3] - uibuf[2];
		  ujj--;
		}
	      }
	    }
	  }
	}
      } else if (model_modifier & MODEL_PTREND) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt_nst) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 1, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  unn = het_cts[urr];
	  upp = pheno_nm_ct - upp; // tot_obs
	  uoo = homcom_cts[urr];
	  ukk += uii * precomp_width;
	  ujj = case_ct - ujj;
	  dxx = orig_chisq[urr];
	  if (model_adapt_nst) {
	    for (umm = uii * precomp_width; umm < ukk; umm++) {
	      ca_trend_precomp_val_bounds(dxx, ujj--, unn, uoo, upp, &(precomp_ui[umm * 4]), nullptr);
	    }
	  } else {
	    for (umm = uii * precomp_width; umm < ukk; umm++) {
	      ca_trend_precomp_val_bounds(dxx, ujj, unn, uoo, upp, &(precomp_ui[umm * 6]), nullptr);
              ca_trend_precomp_val_bounds(maxt_cur_extreme_stat, ujj--, unn, uoo, upp, uibuf, &(precomp_d[umm * 2]));
	      precomp_ui[umm * 6 + 4] = uibuf[2];
	      precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
	    }
	  }
	}
      } else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt_nst) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 1, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  upp = pheno_nm_ct - upp; // tot_obs
	  if (model_modifier & MODEL_PREC) {
	    uoo = upp - homcom_cts[urr] - het_cts[urr]; // col1_sum
	  } else {
	    uoo = homcom_cts[urr];
	  }
	  ukk += uii * precomp_width;
	  ujj = case_ct - ujj;
	  if (model_fisher) {
	    dxx = orig_pvals[urr];
	    if (model_adapt_nst) {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, fisher_midp, ujj--, uoo, upp, &(precomp_ui[umm * 4]), nullptr);
	      }
	    } else {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, fisher_midp, ujj, uoo, upp, &(precomp_ui[umm * 6]), nullptr);
	        fisher22_precomp_pval_bounds(maxt_cur_extreme_stat, fisher_midp, ujj--, uoo, upp, uibuf, &(precomp_d[umm * 2]));
		precomp_ui[umm * 6 + 4] = uibuf[2];
		precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
	      }
	    }
	  } else {
	    dxx = orig_chisq[urr];
	    if (model_adapt_nst) {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, ujj--, uoo, upp, &(precomp_ui[umm * 4]), nullptr);
	      }
	    } else {
	      for (umm = uii * precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(precomp_ui[umm * 6]), nullptr);
		chi22_precomp_val_bounds(maxt_cur_extreme_stat, ujj--, uoo, upp, uibuf, &(precomp_d[umm * 2]));
		precomp_ui[umm * 6 + 4] = uibuf[2];
		precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
	      }
	    }
	  }
	}
      }
      is_last_block = (marker_idx + block_size == marker_unstopped_ct);
      ulii = 0;
      if (model_adapt_nst) {
	if (model_assoc) {
	  if (spawn_threads2(threads, &assoc_adapt_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  assoc_adapt_thread((void*)ulii);
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads2(threads, &model_adapt_domrec_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_domrec_thread((void*)ulii);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads2(threads, &model_adapt_trend_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_trend_thread((void*)ulii);
	} else if (model_modifier & MODEL_PGEN) {
	  if (spawn_threads2(threads, &model_adapt_gen_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_gen_thread((void*)ulii);
	} else {
	  if (spawn_threads2(threads, &model_adapt_best_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_best_thread((void*)ulii);
	}
	join_threads2(threads, max_thread_ct, is_last_block);
      } else {
	g_maxt_block_base = marker_idx;
	if (model_assoc) {
	  if (spawn_threads2(threads, &assoc_maxt_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  assoc_maxt_thread((void*)ulii);
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads2(threads, &model_maxt_domrec_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_domrec_thread((void*)ulii);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads2(threads, &model_maxt_trend_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_trend_thread((void*)ulii);
	} else if (model_modifier & MODEL_PGEN) {
	  if (spawn_threads2(threads, &model_maxt_gen_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_gen_thread((void*)ulii);
	} else {
	  if (spawn_threads2(threads, &model_maxt_best_thread, max_thread_ct, is_last_block)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_best_thread((void*)ulii);
	}
	join_threads2(threads, max_thread_ct, is_last_block);
	ulii = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
        if (model_fisherx) {
	  for (uii = 0; uii < assoc_thread_ct; uii++) {
	    ooptr = &(g_maxt_thread_results[uii * ulii]);
	    for (ujj = perms_done - perm_vec_ct; ujj < perms_done; ujj++) {
	      dxx = *ooptr++;
	      if (dxx < maxt_extreme_stat[ujj]) {
		maxt_extreme_stat[ujj] = dxx;
	      }
	    }
	  }
	} else {
	  for (uii = 0; uii < assoc_thread_ct; uii++) {
	    ooptr = &(g_maxt_thread_results[uii * ulii]);
	    for (ujj = perms_done - perm_vec_ct; ujj < perms_done; ujj++) {
	      dxx = *ooptr++;
	      if (dxx > maxt_extreme_stat[ujj]) {
		maxt_extreme_stat[ujj] = dxx;
	      }
	    }
	  }
	}
      }
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_unstopped_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = (marker_idx * 100LLU) / marker_unstopped_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = (((uint64_t)pct + 1LLU) * marker_unstopped_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (model_perms_nst) {
      bigstack_reset(g_perm_vecs);
    }
    if (fclose_null(&outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (!is_set_test) {
      if (mtest_adjust) {
        if (bigstack_alloc_ui(marker_ct, &marker_idx_to_uidx)) {
	  goto model_assoc_ret_NOMEM;
        }
        fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
        retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, model_fisher? nullptr : orig_chisq, pfilter, output_min_p, mtest_adjust, (!model_assoc) && (!(model_modifier & MODEL_PTREND)), adjust_lambda, nullptr, model_fisher? orig_pvals : nullptr);
        if (retval) {
	  goto model_assoc_ret_1;
        }
        bigstack_reset(marker_idx_to_uidx);
      }
      if (mperm_save & MPERM_DUMP_ALL) {
	g_textbuf[0] = '0';
	wptr = &(g_textbuf[1]);
	a1ptr = &(g_textbuf[MAXLINELEN]);
	if (model_fisherx) {
	  for (uii = 0; uii < marker_ct; uii++) {
	    *wptr++ = ' ';
	    dxx = orig_pvals[uii];
	    if (dxx >= 0) {
	      wptr = dtoa_g(dxx, wptr);
	    } else {
	      wptr = memcpya(wptr, "NA", 2);
	    }
	    if (wptr >= a1ptr) {
	      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
		goto model_assoc_ret_WRITE_FAIL;
	      }
	      wptr = g_textbuf;
	    }
	  }
	} else {
	  for (uii = 0; uii < marker_ct; uii++) {
	    *wptr++ = ' ';
	    dxx = orig_chisq[uii];
	    if (dxx >= 0) {
	      wptr = dtoa_g(dxx, wptr);
	    } else {
	      wptr = memcpya(wptr, "NA", 2);
	    }
	    if (wptr >= a1ptr) {
	      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
		goto model_assoc_ret_WRITE_FAIL;
	      }
	      wptr = g_textbuf;
	    }
	  }
	}
	*wptr++ = '\n';
	if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      }
    } else {
      retval = model_assoc_set_test(threads, bedfile, bed_offset, outname, outname_end, outname_end2, model_modifier, model_mperm_val, pfilter, output_min_p, mtest_adjust, unfiltered_marker_ct, marker_exclude_orig, marker_ct_orig, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, sex_male, apip, pheno_nm_ct, pheno_nm, founder_pnm, gender_req, ld_ignore_x, hh_exists, perm_batch_size, sip, loadbuf_raw);
      if (retval) {
        goto model_assoc_ret_1;
      }
    }
  }
  if (model_perms_nst) {
    if (mperm_save & MPERM_DUMP_ALL) {
      if (perm_pass_idx) {
	putc_unlocked(' ', stdout);
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = perm_vec_ct;
      ujj = 1 + perms_done - ulii;
      wptr = g_textbuf;
      a1ptr = &(g_textbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32toa(uii + ujj, wptr);
        orig_pvals_ptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = orig_pvals_ptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = dtoa_g(dxx, wptr);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= a1ptr) {
	    if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    wptr = g_textbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	goto model_assoc_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    bigstack_reset(g_perm_vecs);
    if (perms_done < perms_total) {
      if (model_adapt_nst) {
	marker_unstopped_ct = marker_ct - popcount01_longs((uintptr_t*)perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
	if (!marker_unstopped_ct) {
	  goto model_assoc_adapt_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", perms_done, (perms_done != 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto model_assoc_more_perms;
    }
    if (model_adapt_nst) {
    model_assoc_adapt_perm_count:
      perms_done = 0;
      for (uii = 0; uii < marker_ct; uii++) {
	if (perm_attempt_ct[uii] > perms_done) {
	  perms_done = perm_attempt_ct[uii];
	  if (perms_done == perms_total) {
	    break;
	  }
	}
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("%u %s permutation%s complete.\n", perms_done, model_maxt_nst? "max(T)" : "(adaptive)", (perms_done != 1)? "s" : "");
    if (model_fisher && (model_modifier & MODEL_PTREND)) {
      outname_end2 -= 7; // remove ".fisher"
    }
    if (model_adapt_nst) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	if (bigstack_alloc_c(FNAMESIZE, &a1ptr)) {
	  goto model_assoc_ret_NOMEM;
	}
	ulii = outname_end - outname;
	memcpy(a1ptr, outname, ulii);
	memcpy(&(a1ptr[ulii]), ".mperm.dump.best", 17);
	LOGPRINTFWW("Dumping best permutation %svalues to %s .\n", model_fisherx? "p-" : "chi-square ", a1ptr);
	if (fopen_checked(a1ptr, "w", &outfile)) {
	  goto model_assoc_ret_OPEN_FAIL;
	}
	dxx = 0;
	if (model_fisherx) {
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    if ((orig_pvals[marker_idx] != -9) && (orig_pvals[marker_idx] < dxx)) {
	      dxx = orig_pvals[marker_idx];
	    }
	  }
	  dxx = 1 - dxx;
	} else {
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    if (orig_chisq[marker_idx] > dxx) {
	      dxx = orig_chisq[marker_idx];
	    }
	  }
	}
        memcpy(g_textbuf, "0 ", 2);
	wptr = dtoa_gx(dxx, '\n', &(g_textbuf[2]));
	if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile)) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
	for (uii = 0; uii < perms_total; uii++) {
	  wptr = uint32toa_x(uii + 1, ' ', g_textbuf);
	  wptr = dtoa_gx(maxt_extreme_stat[uii], '\n', wptr);
	  if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    if (model_adapt_nst) {
      sprintf(g_textbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(g_textbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(maxt_extreme_stat, &(maxt_extreme_stat[perms_total]));
#else
      qsort(maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
    }
    /*
    if (model_maxt_nst) {
      printf("extreme stats: %g %g\n", maxt_extreme_stat[0], maxt_extreme_stat[perms_total - 1]);
    }
    */
    fprintf(outfile, g_textbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    while (1) {
      while (1) {
	do {
          chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	is_x = (uii == (uint32_t)x_code);
	if (model_assoc || (((!IS_SET(haploid_mask, uii)) && (uii != (uint32_t)mt_code)) || is_x)) {
	  break;
	}
	marker_uidx = next_unset_unsafe(marker_exclude, chrom_end);
      }
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, uii, g_textbuf));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	if (model_adapt_nst) {
	  pval = ((double)(perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(perm_2success_ct[marker_idx] + 2)) * dxx;
	}
        if (pval <= pfilter) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	  wptr = &(wptr_start[1 + plink_maxsnp]);
	  if ((!model_assoc) && ((model_adapt_nst && (!perm_attempt_ct[marker_idx])) || ((!model_adapt_nst) && ((model_fisherx && (orig_pvals[marker_idx] == -9)) || ((!model_fisherx) && (orig_chisq[marker_idx] == -9)))))) {
	    // invalid
            wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    if (!model_perm_count) {
	      wptr = dtoa_g_wxp4x(pval, 12, ' ', wptr);
	    } else {
	      wptr = dtoa_g_wxp4x(((double)perm_2success_ct[marker_idx]) * 0.5, 12, ' ', wptr);
	    }
	    if (model_adapt_nst) {
	      wptr = memseta(wptr, 32, 2);
	      wptr = uint32toa_w10(perm_attempt_ct[marker_idx], wptr);
	    } else {
	      if (model_fisherx) {
		// minimum p-value
		dzz = (int32_t)(doublearr_greater_than(maxt_extreme_stat, perms_total, orig_pvals[marker_idx] * (1.0 + EPSILON)) + 1);
	      } else {
		// maximum chisq
		dzz = (int32_t)(perms_total - doublearr_greater_than(maxt_extreme_stat, perms_total, orig_chisq[marker_idx] - EPSILON) + 1);
	      }
	      if (!model_perm_count) {
		wptr = dtoa_g_wxp4(dzz * dyy, 12, wptr);
	      } else {
		wptr = dtoa_g_wxp4(dzz - 1, 12, wptr);
	      }
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto model_assoc_loop_end;
	}
	marker_uidx++;
        next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      }
    }
  model_assoc_loop_end:
    if (fclose_null(&outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
  }

  while (0) {
  model_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  model_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  model_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  model_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  model_assoc_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  model_assoc_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 model_assoc_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
}

int32_t qassoc_set_test(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_mperm_val, double pfilter, double output_min_p, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, uintptr_t* marker_exclude_mid, uintptr_t marker_ct_mid, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sex_male, Aperm_info* apip, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* founder_pnm, uintptr_t* sample_include2, uintptr_t* sample_male_include2, uint32_t ld_ignore_x, uint32_t hh_exists, uint32_t hh_or_mt_exists, uint32_t perm_batch_size, Set_info* sip, uint32_t* tcnt, uintptr_t* loadbuf_raw) {
  // Similar to glm_linear_assoc_set_test().
  // Side effect: t-statistics in g_orig_chisq[] are clobbered and replaced
  // with same-p-value 1df chi-square statistics.
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t* marker_exclude = marker_exclude_mid;
  uintptr_t* unstopped_markers = nullptr;
  uintptr_t* loadbuf = g_loadbuf;
  uintptr_t* perm_adapt_set_unstopped = nullptr;
  uintptr_t* regression_skip = nullptr;
  double* orig_stats = g_orig_chisq; // initially contains t-statistics
  double* sorted_chisq_buf = nullptr;
  uint32_t* marker_idx_to_uidx = nullptr;
  uint32_t* sorted_marker_idx_buf = nullptr;
  uint32_t* proxy_arr = nullptr;
  uint32_t* perm_2success_ct = nullptr;
  uint32_t* perm_attempt_ct = nullptr;
  uintptr_t marker_ct = marker_ct_mid;
  uintptr_t set_ct = 0;
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  double adaptive_ci_zt = 0.0;
  uint32_t max_thread_ct = g_thread_ct;
  uint32_t perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t perms_done = 0;
  int32_t retval = 0;
  unsigned char* bigstack_mark2;
  uintptr_t* set_incl;
  uintptr_t* loadbuf_ptr;
  double* orig_set_scores;
  double* chisq_pmajor;
  double* read_dptr;
  double* write_dptr;
  uint32_t** setdefs;
  uint32_t** ld_map;
  uintptr_t marker_ctl;
  uintptr_t marker_midx;
  uintptr_t set_idx;
  uintptr_t perm_vec_ct;
  uintptr_t perm_vec_ctcl8m;
  uintptr_t pidx;
  uintptr_t ulii;
  double chisq_threshold;
  double dxx;
  double dyy;
  uint32_t perms_total;
  uint32_t max_sigset_size;
  uint32_t marker_unstopped_ct;
  uint32_t is_last_block;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t first_adapt_check;
  uint32_t marker_uidx;
  uint32_t marker_idx;
  uint32_t marker_idx2;
  uint32_t marker_bidx;
  uint32_t skip_ct;
  uint32_t uii;
  if (sip->set_test_lambda > 1.0) {
    dxx = 1.0 / sip->set_test_lambda;
  } else {
    dxx = 1.0;
  }
  for (marker_midx = 0; marker_midx < marker_ct; marker_midx++) {
    dyy = calc_tprob(orig_stats[marker_midx], tcnt[marker_midx]);
    if (dyy == 0.0) {
      dyy = MAX_INVERSE_CHIPROB_1DF * dxx;
    } else {
      orig_stats[marker_midx] = inverse_chiprob(dyy, 1) * dxx;
    }
  }
  retval = set_test_common_init(threads, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude_orig, marker_ct_orig, marker_ids, max_marker_id_len, marker_reverse, orig_stats, sip, chrom_info_ptr, unfiltered_sample_ct, sex_male, founder_pnm, ld_ignore_x, hh_exists, "QT --assoc", &marker_ct, &marker_exclude, &set_incl, &marker_idx_to_uidx, &setdefs, &set_ct, &max_sigset_size, &ld_map, &chisq_threshold, &orig_set_scores, &sorted_chisq_buf, &sorted_marker_idx_buf, &proxy_arr, &perm_adapt_set_unstopped, &perm_2success_ct, &perm_attempt_ct, &unstopped_markers);
  if (retval) {
    goto qassoc_set_test_ret_1;
  }
  if (!set_ct) {
    goto qassoc_set_test_write;
  }
  marker_ctl = BITCT_TO_WORDCT(marker_ct);
  if (marker_ct_mid != marker_ct) {
    inplace_delta_collapse_arr((char*)tcnt, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
    inplace_delta_collapse_arr((char*)g_missing_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
    inplace_delta_collapse_arr((char*)g_het_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
    inplace_delta_collapse_arr((char*)g_homcom_cts, sizeof(int32_t), marker_ct_mid, marker_ct, marker_exclude_mid, marker_exclude);
  }
  if (bigstack_calloc_ul(marker_ctl, &regression_skip)) {
    goto qassoc_set_test_ret_NOMEM;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    // nanal
    uii = tcnt[marker_idx] + 2;
    if ((uii == 2) || (g_homcom_cts[marker_idx] == uii) || (g_het_cts[marker_idx] == uii) || (g_het_cts[marker_idx] + g_homcom_cts[marker_idx] == 0)) {
      // 0 df or no genotype variation, regression always fails
      SET_BIT(marker_idx, regression_skip);
    }
  }
  if (model_modifier & MODEL_PERM) {
    perms_total = apip->max;
    first_adapt_check = (apip->min < apip->init_interval)? ((int32_t)apip->init_interval) : apip->min;
    adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)set_ct)));
  } else {
    perms_total = model_mperm_val;
    first_adapt_check = perms_total + 1;
  }
  for (uii = 0; uii < set_ct; uii++) {
    perm_attempt_ct[uii] = perms_total;
  }
  if (max_thread_ct > perms_total) {
    max_thread_ct = perms_total;
  }
  if (bigstack_init_sfmtp(max_thread_ct)) {
    goto qassoc_set_test_ret_NOMEM;
  }

  bigstack_mark2 = g_bigstack_base;
 qassoc_set_test_more_perms:
  bitvec_and(unstopped_markers, marker_ctl, regression_skip);
  bitvec_andnot(regression_skip, marker_ctl, unstopped_markers);
  skip_ct = popcount_longs(regression_skip, marker_ctl);
  marker_unstopped_ct = popcount_longs(unstopped_markers, marker_ctl);

  if (perms_done) {
    uii = apip->init_interval;
    while (first_adapt_check <= perms_done) {
      first_adapt_check += (int32_t)(uii + ((int32_t)first_adapt_check) * apip->interval_slope);
    }
  }
  perm_vec_ct = perm_batch_size;
  // possible todo: split first batch to reduce adaptive overshoot
  if (perm_vec_ct > perms_total - perms_done) {
    perm_vec_ct = perms_total - perms_done;
  }
  g_perm_vec_ct = perm_vec_ct;
  if (perm_vec_ct >= CACHELINE_INT32 * max_thread_ct) {
    g_perm_generation_thread_ct = max_thread_ct;
  } else {
    g_perm_generation_thread_ct = MAXV(perm_vec_ct / CACHELINE_INT32, 1);
  }
  perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  if (bigstack_alloc_d(perm_vec_ctcl8m * pheno_nm_ct, &g_perm_vecstd) ||
      bigstack_calloc_d(perm_vec_ctcl8m * 3 * max_thread_ct, &g_thread_git_qbufs)) {
    goto qassoc_set_test_ret_NOMEM;
  }

  ulii = 0;
  if (!g_perm_cluster_ct) {
    if (spawn_threads(threads, &generate_qt_perms_smajor_thread, g_perm_generation_thread_ct)) {
      goto qassoc_set_test_ret_THREAD_CREATE_FAIL;
    }
    generate_qt_perms_smajor_thread((void*)ulii);
  } else {
    if (spawn_threads(threads, &generate_qt_cluster_perms_smajor_thread, g_perm_generation_thread_ct)) {
      goto qassoc_set_test_ret_THREAD_CREATE_FAIL;
    }
    generate_qt_cluster_perms_smajor_thread((void*)ulii);
  }
  join_threads(threads, g_perm_generation_thread_ct);
  if (bigstack_alloc_d(MODEL_BLOCKSIZE * perm_vec_ct, &g_mperm_save_all) ||
      bigstack_alloc_d(marker_ct * perm_vec_ct, &chisq_pmajor)) {
    goto qassoc_set_test_ret_NOMEM;
  }
  for (pidx = 0; pidx < perm_vec_ct; pidx++) {
    write_dptr = &(chisq_pmajor[pidx * marker_ct]);
    for (marker_idx = 0, marker_idx2 = 0; marker_idx < skip_ct; marker_idx++, marker_idx2++) {
      next_set_unsafe_ck(regression_skip, &marker_idx2);
      write_dptr[marker_idx2] = -9;
    }
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto qassoc_set_test_ret_READ_FAIL;
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  do {
    if (marker_uidx >= chrom_end) {
      // exploit overflow
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &uii, &g_min_ploidy_1);
      g_min_ploidy_1 |= uii; // treat MT as haploid
    }
    block_size = 0;
    block_end = marker_unstopped_ct - marker_idx;
    if (block_end > MODEL_BLOCKSIZE) {
      block_end = MODEL_BLOCKSIZE;
    }
    do {
      if (!IS_SET(unstopped_markers, marker_idx2)) {
        do {
	  marker_uidx++;
	  next_unset_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
        } while ((marker_uidx < chrom_end) && (!IS_SET(unstopped_markers, marker_idx2)));
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto qassoc_set_test_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(loadbuf[block_size * pheno_nm_ctv2]);
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf_ptr)) {
	goto qassoc_set_test_ret_READ_FAIL;
      }
      if (g_min_ploidy_1 && hh_or_mt_exists) {
	haploid_fix(hh_or_mt_exists, sample_include2, sample_male_include2, pheno_nm_ct, g_is_x, g_is_y, (unsigned char*)loadbuf_ptr);
      }
      g_adapt_m_table[block_size] = marker_idx2++;
      block_size++;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto qassoc_set_test_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (!block_size) {
      continue;
    }
    is_last_block = (marker_idx + block_size >= marker_unstopped_ct);
    g_block_diff = block_size;
    ulii = 0;
    if (spawn_threads2(threads, &qassoc_set_thread, max_thread_ct, is_last_block)) {
      goto qassoc_set_test_ret_THREAD_CREATE_FAIL;
    }
    qassoc_set_thread((void*)ulii);
    join_threads2(threads, max_thread_ct, is_last_block);

    // convert to equivalent chi-square stats and transpose
    // (conversion has to be done here since dcdflib is not thread-safe)
    read_dptr = g_mperm_save_all;
    for (marker_bidx = 0; marker_bidx < block_size; marker_bidx++) {
      uii = g_adapt_m_table[marker_bidx];
      write_dptr = &(chisq_pmajor[uii]);
      uii = tcnt[uii];
      dyy = inverse_tprob(sip->set_p, uii);
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	dxx = *read_dptr++;
	if (dxx < dyy) {
	  dxx = -9;
	} else {
	  dxx = calc_tprob(dxx, uii);
	  if (dxx == 0.0) {
	    dxx = MAX_INVERSE_CHIPROB_1DF;
	  } else {
	    dxx = inverse_chiprob(dxx, 1);
	  }
	}
	// this is cache-unfriendly, may want to update in-place instead and
	// separate out the transpose
	write_dptr[pidx * marker_ct] = dxx;
      }
    }
    marker_idx += block_size;
  } while (marker_idx < marker_unstopped_ct);
  perms_done += perm_vec_ct;
  compute_set_scores(marker_ct, perm_vec_ct, set_ct, chisq_pmajor, orig_set_scores, sorted_chisq_buf, sorted_marker_idx_buf, proxy_arr, setdefs, ld_map, apip, chisq_threshold, adaptive_ci_zt, first_adapt_check, perms_done, sip->set_max, perm_adapt_set_unstopped, perm_2success_ct, perm_attempt_ct);
  bigstack_reset(bigstack_mark2);
  if (perms_done < perms_total) {
    if (model_modifier & MODEL_PERM) {
      if (!extract_set_union(setdefs, set_ct, perm_adapt_set_unstopped, unstopped_markers, marker_ct)) {
	perms_done = 0;
	for (set_idx = 0; set_idx < set_ct; set_idx++) {
          if (perms_done < perm_attempt_ct[set_idx]) {
	    perms_done = perm_attempt_ct[set_idx];
	  }
	}
	goto qassoc_set_test_perms_done;
      }
    }
    printf("\r%u permutation%s complete.", perms_done, (perms_done != 1)? "s" : "");
    fflush(stdout);
    goto qassoc_set_test_more_perms;
  }
 qassoc_set_test_perms_done:
  putc_unlocked('\r', stdout);
  LOGPRINTF("%u permutation%s complete.\n", perms_done, (perms_done != 1)? "s" : "");
 qassoc_set_test_write:
  if (model_modifier & MODEL_PERM) {
    memcpy(outname_end, ".qassoc.set.perm", 17);
  } else {
    memcpy(outname_end, ".qassoc.set.mperm", 18);
  }
  retval = write_set_test_results(outname, &(outname_end[11]), sip, ld_map, setdefs, set_incl, set_ct, marker_ct_orig, marker_ct, marker_idx_to_uidx, marker_ids, max_marker_id_len, perm_2success_ct, perm_attempt_ct, mtest_adjust, perm_count, pfilter, output_min_p, chisq_threshold, orig_stats, sorted_chisq_buf, sorted_marker_idx_buf, proxy_arr);
  while (0) {
  qassoc_set_test_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  qassoc_set_test_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  qassoc_set_test_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 qassoc_set_test_ret_1:
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t qassoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_mperm_val, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* founder_info, uintptr_t* sex_male, uint32_t hh_exists, uint32_t ld_ignore_x, uint32_t perm_batch_size, Set_info* sip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t marker_ct = marker_ct_orig;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uintptr_t perm_vec_ctcl8m = 0;
  FILE* outfile = nullptr;
  FILE* outfile_qtm = nullptr;
  FILE* outfile_msa = nullptr;
  uint32_t is_set_test = model_modifier & MODEL_SET_TEST;
  uint32_t perm_adapt_nst = (model_modifier & MODEL_PERM) && (!is_set_test);
  uint32_t perm_maxt_nst = (model_modifier & MODEL_MPERM) && (!is_set_test);
  uint32_t do_perms = model_modifier & (MODEL_PERM | MODEL_MPERM);
  uint32_t do_perms_nst = do_perms && (!is_set_test);
  uint32_t qt_means = model_modifier & MODEL_QT_MEANS;
  uint32_t do_lin = model_modifier & MODEL_LIN;
  uint32_t qt_means_or_lin = qt_means || do_lin;
  uint32_t perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t fill_orig_chiabs = do_perms || mtest_adjust;
  uint32_t perms_total = 0;
  uint32_t pct = 0;
  uint32_t max_thread_ct = g_thread_ct;
  uint32_t perm_pass_idx = 0;
  uint32_t mt_exists = (chrom_info_ptr->xymt_codes[MT_OFFSET] != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[MT_OFFSET]);
  uint32_t hh_or_mt_exists = hh_exists | (mt_exists * NXMHH_EXISTS);
  int32_t retval = 0;
  double x11 = 0;
  double x12 = 0;
  double x22 = 0;
  uintptr_t* marker_exclude = marker_exclude_orig;
  uintptr_t* founder_pnm = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uint32_t* tcnt = nullptr;
  char* chrom_name_ptr = nullptr;
  uint32_t chrom_name_len = 0;
  char chrom_name_buf[5];
  uint32_t mu_table[MODEL_BLOCKSIZE];
  char numbuf[16]; // ' -1.23456e-200\0' fits, barely
  char spacebuf[8];
  char* outname_end2;
  char* wptr_start;
  char* wptr;
  char* wptr_restart;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* lbptr2;
  uintptr_t* sample_include2;
  double* ooptr;
  double* dptr;
  double* dptr2;
  double* dptr3;
  uint32_t* marker_idx_to_uidx;
  uint32_t marker_unstopped_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t marker_bidx;
  uintptr_t marker_uidx; // loading
  uintptr_t marker_uidx2; // writing
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t sample_uidx;
  uintptr_t sample_uidx_stop;
  uintptr_t sample_idx;
  uintptr_t ulii;
  intptr_t geno_sum;
  intptr_t nanal;
  intptr_t geno_ssq;
  double nanal_recip;
  double qt_sum;
  double qt_ssq;
  double qt_g_prod;
  double qt_g_prod_centered;
  double qt_mean;
  double geno_mean;
  double qt_var;
  double geno_var;
  double qt_g_covar;
  double beta;
  double vbeta_sqrt;
  double tstat;
  double tp;
  double rsq;
  double qt_het_sum;
  double qt_het_ssq;
  double qt_homrar_sum;
  double qt_homrar_ssq;
  double qt_homcom_sum;
  double dxx;
  double dyy;
  double dzz;
  double pval;
  uint32_t homrar_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t is_last_block;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  char* a1ptr;
  char* a2ptr;
  if (pheno_nm_ct < 2) {
    logerrprint("Warning: Skipping QT --assoc since less than two phenotypes are present.\n");
    goto qassoc_ret_1;
  }
  if (is_set_test) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl, &founder_pnm)) {
      goto qassoc_ret_NOMEM;
    }
    memcpy(founder_pnm, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
    bitvec_and(founder_info, unfiltered_sample_ctl, founder_pnm);
    if (extract_set_union_unfiltered(sip, nullptr, unfiltered_marker_ct, marker_exclude_orig, &marker_exclude, &marker_ct)) {
      goto qassoc_ret_NOMEM;
    }
  }
  memset(spacebuf, 32, 8);
  g_perm_pheno_nm_ct = pheno_nm_ct;
  g_perms_done = 0;
  g_mperm_save_all = nullptr;
  numbuf[0] = ' ';
  if (perm_maxt_nst) {
    perms_total = model_mperm_val;
    // square of t-stat
    if (bigstack_calloc_d(perms_total, &g_maxt_extreme_stat)) {
      goto qassoc_ret_NOMEM;
    }
    g_ldrefs = (uint16_t*)bigstack_alloc(marker_ct * sizeof(int16_t));
    if (!g_ldrefs) {
      goto qassoc_ret_NOMEM;
    }
#ifdef __LP64__
    fill_ulong_one((marker_ct + 3) / 4, (uintptr_t*)g_ldrefs);
#else
    fill_ulong_one((marker_ct + 1) / 2, (uintptr_t*)g_ldrefs);
#endif
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(outname, "w", &outfile_msa)) {
	goto qassoc_ret_OPEN_FAIL;
      }
      if (putc_checked('0', outfile_msa)) {
	goto qassoc_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("Dumping all permutation squared %sstats to %s .\n", do_lin? "Lin " : "Wald t-", outname);
    }
  } else {
    mperm_save = 0;
    if (perm_adapt_nst) {
      g_aperm_alpha = apip->alpha;
      perms_total = apip->max;
      if (bigstack_alloc_ui(marker_ct, &g_perm_attempt_ct) ||
	  bigstack_calloc_uc(round_up_pow2(marker_ct, BYTECT), &g_perm_adapt_stop)) {
	goto qassoc_ret_NOMEM;
      }
      ujj = apip->max;
      for (uii = 0; uii < marker_ct; uii++) {
	g_perm_attempt_ct[uii] = ujj;
      }
      g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)marker_ct)));
      if (apip->min < apip->init_interval) {
	g_first_adapt_check = (int32_t)(apip->init_interval);
      } else {
	g_first_adapt_check = apip->min;
      }
      g_adaptive_intercept = apip->init_interval;
      g_adaptive_slope = apip->interval_slope;
    }
  }
  outname_end2 = memcpyb(outname_end, ".qassoc", 8);
  if (bigstack_alloc_ul(unfiltered_sample_ctv2, &loadbuf_raw)) {
    goto qassoc_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctv2 - 2] = 0;
  loadbuf_raw[unfiltered_sample_ctv2 - 1] = 0;
  if (fill_orig_chiabs) {
    if (bigstack_alloc_d(marker_ct, &g_orig_chisq)) {
      goto qassoc_ret_NOMEM;
    }
    if (mtest_adjust || is_set_test) {
      if (bigstack_alloc_ui(marker_ct, &tcnt)) {
	goto qassoc_ret_NOMEM;
      }
    }
  }
  if (fopen_checked(outname, "w", &outfile)) {
    goto qassoc_ret_OPEN_FAIL;
  }
  if (qt_means) {
    memcpy(outname_end2, ".means", 7);
    if (fopen_checked(outname, "w", &outfile_qtm)) {
      goto qassoc_ret_OPEN_FAIL;
    }
    sprintf(g_textbuf, " CHR %%%us  VALUE      G11      G12      G22\n", plink_maxsnp);
    fprintf(outfile_qtm, g_textbuf, "SNP");
    *outname_end2 = '\0';
  }
  if (haploid_chrom_present(chrom_info_ptr) || mt_exists) {
    logerrprint("Warning: QT --assoc doesn't handle X/Y/MT/haploid variants normally (try\n--linear).\n");
  }
  LOGPRINTFWW5("Writing QT --assoc report to %s ... ", outname);
  fflush(stdout);
  sprintf(g_textbuf, " CHR %%%us         BP    NMISS       BETA         SE         R2        T            P ", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  if (do_lin) {
    fputs("         LIN        LIN_P ", outfile);
  }
  if (putc_checked('\n', outfile)) {
    goto qassoc_ret_WRITE_FAIL;
  }
  if (do_perms) {
    if (model_modifier & MODEL_PERM) {
      if (perm_batch_size > apip->max) {
	perm_batch_size = apip->max;
      }
    } else {
      if (perm_batch_size > model_mperm_val) {
	perm_batch_size = model_mperm_val;
      }
    }
    uii = MINV(perm_batch_size, perms_total) / CACHELINE_DBL;
    if (max_thread_ct > uii) {
      max_thread_ct = MAXV(uii, 1);
    }
    if (cluster_starts) {
      retval = cluster_include_and_reindex(unfiltered_sample_ct, pheno_nm, 1, nullptr, pheno_nm_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_perm_cluster_ct, &g_perm_cluster_map, &g_perm_cluster_starts, nullptr, nullptr);
      if (retval) {
	goto qassoc_ret_1;
      }
      if (!g_perm_cluster_ct) {
        logerrprint("Error: No size 2+ clusters for permutation test.\n");
        goto qassoc_ret_INVALID_CMDLINE;
      }
      if (bigstack_alloc_ui(pheno_nm_ct, &g_perm_sample_to_cluster) ||
          bigstack_alloc_ui(max_thread_ct * round_up_pow2(g_perm_cluster_ct, CACHELINE_INT32), &g_perm_qt_cluster_thread_wkspace)) {
	goto qassoc_ret_NOMEM;
      }
      fill_unfiltered_sample_to_cluster(pheno_nm_ct, g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts, g_perm_sample_to_cluster);
    }
    if (bigstack_alloc_ui(marker_ct, &g_missing_cts) ||
	bigstack_alloc_ui(marker_ct, &g_het_cts) ||
	bigstack_alloc_ui(marker_ct, &g_homcom_cts)) {
      goto qassoc_ret_NOMEM;
    }
    if (!is_set_test) {
      if (bigstack_init_sfmtp(max_thread_ct)) {
	goto qassoc_ret_NOMEM;
      }
      if (bigstack_calloc_ui(marker_ct, &g_perm_2success_ct)) {
	goto qassoc_ret_NOMEM;
      }
    }
  }
  if (do_lin) {
    if (bigstack_alloc_d(marker_ct, &g_orig_linsq)) {
      goto qassoc_ret_NOMEM;
    }
  }
  if (bigstack_alloc_ul(MODEL_BLOCKSIZE * pheno_nm_ctv2, &g_loadbuf) ||
      bigstack_alloc_ui(marker_ct, &marker_idx_to_uidx) ||
      bigstack_alloc_ul(pheno_nm_ctv2, &sample_include2)) {
    goto qassoc_ret_NOMEM;
  }
  fill_quatervec_55(pheno_nm_ct, sample_include2);
  if (alloc_collapsed_haploid_filters(pheno_nm, sex_male, unfiltered_sample_ct, pheno_nm_ct, hh_or_mt_exists, 1, &sample_include2, &sample_male_include2)) {
    goto qassoc_ret_NOMEM;
  }
  marker_unstopped_ct = marker_ct;
  if (bigstack_alloc_d(pheno_nm_ct, &g_perm_pheno_d2)) {
    goto qassoc_ret_NOMEM;
  }
  g_pheno_sum = 0;
  g_pheno_ssq = 0;
  sample_uidx = 0;
  sample_idx = 0;
  dptr = g_perm_pheno_d2;
  do {
    sample_uidx = next_set_ul_unsafe(pheno_nm, sample_uidx);
    sample_uidx_stop = next_unset_ul(pheno_nm, sample_uidx, unfiltered_sample_ct);
    sample_idx += sample_uidx_stop - sample_uidx;
    dptr2 = &(pheno_d[sample_uidx]);
    sample_uidx = sample_uidx_stop;
    dptr3 = &(pheno_d[sample_uidx_stop]);
    do {
      dxx = *dptr2++;
      *dptr++ = dxx;
      g_pheno_sum += dxx;
      g_pheno_ssq += dxx * dxx;
    } while (dptr2 < dptr3);
  } while (sample_idx < pheno_nm_ct);
  fputs("0%", stdout);
  fflush(stdout);

  // ----- begin main loop -----
 qassoc_more_perms:
  if (do_perms_nst) {
    if (perm_adapt_nst && perm_pass_idx) {
      while (g_first_adapt_check <= g_perms_done) {
	// APERM_MAX prevents infinite loop here
	g_first_adapt_check += (int32_t)(apip->init_interval + ((int32_t)g_first_adapt_check) * apip->interval_slope);
      }
    }
    // g_perm_vec_ct memory allocation dependencies:
    //   g_maxt_thread_results: (8 * perm_vec_ct, CL-aligned) * thread_ct
    //   g_perm_vecstd: (8 * perm_vec_ct, CL-aligned) * pheno_nm_ct
    //   g_mperm_save_all (if needed): marker_ct * 8 * perm_vec_ct
    //   adaptive, Wald:
    //     g_thread_git_qbufs: (8 * perm_vec_ct, CL-aligned) * 3 * thread_ct
    //   adaptive, Lin:
    //     g_thread_git_qbufs: (8 * perm_vec_ct, CL-aligned) * 6 * thread_ct
    //   max(T), Wald:
    //     g_qresultbuf: MODEL_BLOCKSIZE * (8 * perm_vec_ct, CL-aligned) * 3
    //   max(T), Lin:
    //     g_qresultbuf: MODEL_BLOCKSIZE * (8 * perm_vec_ct, CL-aligned) * 6
    g_perm_vec_ct = perm_batch_size;
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    }
    perm_vec_ctcl8m = round_up_pow2(g_perm_vec_ct, CACHELINE_DBL);
    if (bigstack_alloc_d(perm_vec_ctcl8m * pheno_nm_ct, &g_perm_vecstd)) {
      goto qassoc_ret_NOMEM;
    }
    ulii = do_lin? 6 : 3;
    if (perm_maxt_nst) {
      if (bigstack_alloc_d(max_thread_ct * perm_vec_ctcl8m, &g_maxt_thread_results) ||
	  bigstack_alloc_d(ulii * MODEL_BLOCKSIZE * perm_vec_ctcl8m, &g_qresultbuf)) {
	goto qassoc_ret_NOMEM;
      }
      if (mperm_save & MPERM_DUMP_ALL) {
	if (bigstack_alloc_d(marker_ct * g_perm_vec_ct, &g_mperm_save_all)) {
	  goto qassoc_ret_NOMEM;
	}
      }
    } else {
      if (bigstack_calloc_d(perm_vec_ctcl8m * ulii * max_thread_ct, &g_thread_git_qbufs)) {
	goto qassoc_ret_NOMEM;
      }
    }
    g_perms_done += g_perm_vec_ct;
    if (g_perm_vec_ct >= CACHELINE_DBL * max_thread_ct) {
      g_perm_generation_thread_ct = max_thread_ct;
    } else {
      g_perm_generation_thread_ct = MAXV(g_perm_vec_ct / CACHELINE_DBL, 1);
    }
    ulii = 0;
    if (!cluster_starts) {
      if (spawn_threads(threads, &generate_qt_perms_smajor_thread, g_perm_generation_thread_ct)) {
	goto qassoc_ret_THREAD_CREATE_FAIL;
      }
      generate_qt_perms_smajor_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &generate_qt_cluster_perms_smajor_thread, g_perm_generation_thread_ct)) {
	goto qassoc_ret_THREAD_CREATE_FAIL;
      }
      generate_qt_cluster_perms_smajor_thread((void*)ulii);
    }
    join_threads(threads, g_perm_generation_thread_ct);
    g_assoc_thread_ct = max_thread_ct;
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto qassoc_ret_READ_FAIL;
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  loop_end = marker_ct / 100;
  do {
    if (marker_uidx >= chrom_end) {
      g_qblock_start = 0;
      // exploit overflow
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &uii, &g_min_ploidy_1);
      g_min_ploidy_1 |= uii; // treat MT as haploid
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, uii, &chrom_name_len, chrom_name_buf);
    } else if (perm_maxt_nst) {
      marker_idx -= MODEL_BLOCKKEEP;
      memcpy(g_loadbuf, &(g_loadbuf[(MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * pheno_nm_ctv2]), MODEL_BLOCKKEEP * pheno_nm_ctv2 * sizeof(intptr_t));
      if (!do_lin) {
	memcpy(g_qresultbuf, &(g_qresultbuf[3 * (MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * perm_vec_ctcl8m]), MODEL_BLOCKKEEP * perm_vec_ctcl8m * 3 * sizeof(double));
      } else {
	memcpy(g_qresultbuf, &(g_qresultbuf[6 * (MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * perm_vec_ctcl8m]), MODEL_BLOCKKEEP * perm_vec_ctcl8m * 6 * sizeof(double));
      }
      g_qblock_start = MODEL_BLOCKKEEP;
    } else {
      g_qblock_start = 0;
    }
    block_size = g_qblock_start;
    block_end = marker_unstopped_ct - marker_idx;
    if (block_end > MODEL_BLOCKSIZE) {
      block_end = MODEL_BLOCKSIZE;
    }
    do {
      if (perm_adapt_nst && g_perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx++;
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto qassoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(g_loadbuf[block_size * pheno_nm_ctv2]);
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf_ptr)) {
	goto qassoc_ret_READ_FAIL;
      }
      if (g_min_ploidy_1 && hh_or_mt_exists) {
	haploid_fix(hh_or_mt_exists, sample_include2, sample_male_include2, pheno_nm_ct, g_is_x, g_is_y, (unsigned char*)loadbuf_ptr);
      }
      if (perm_adapt_nst) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto qassoc_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (block_size == g_qblock_start) {
      continue;
    }
    if (!perm_pass_idx) {
      for (marker_bidx = g_qblock_start; marker_bidx < block_size; marker_bidx++) {
	marker_uidx2 = mu_table[marker_bidx];
        marker_idx_to_uidx[marker_idx + marker_bidx] = marker_uidx2;
	loadbuf_ptr = &(g_loadbuf[marker_bidx * pheno_nm_ctv2]);
	genovec_3freq(loadbuf_ptr, sample_include2, pheno_nm_ctv2, &missing_ct, &het_ct, &homcom_ct);
	nanal = pheno_nm_ct - missing_ct;
	wptr = memcpya(g_textbuf, chrom_name_ptr, chrom_name_len);
	*wptr++ = ' ';
        wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	*wptr++ = ' ';
	wptr = uint32toa_w10x(marker_pos[marker_uidx2], ' ', wptr);
	wptr = uint32toa_w8x(nanal, ' ', wptr);
	homrar_ct = nanal - het_ct - homcom_ct;
	if (do_perms) {
	  g_missing_cts[marker_idx + marker_bidx] = missing_ct;
	  g_homcom_cts[marker_idx + marker_bidx] = homcom_ct;
	  g_het_cts[marker_idx + marker_bidx] = het_ct;
	}
	geno_sum = 2 * homrar_ct + het_ct;
	geno_ssq = 4 * homrar_ct + het_ct;
	qt_sum = g_pheno_sum;
	qt_g_prod = 0;
	qt_ssq = g_pheno_ssq;
	lbptr2 = loadbuf_ptr;
	uii = 0;
	qt_het_sum = 0;
	qt_het_ssq = 0;
	qt_homrar_sum = 0;
	qt_homrar_ssq = 0;
	do {
	  ulii = ~(*lbptr2++);
	  if (uii + BITCT2 > pheno_nm_ct) {
	    ulii &= (ONELU << ((pheno_nm_ct & (BITCT2 - 1)) * 2)) - ONELU;
	  }
	  while (ulii) {
	    ujj = CTZLU(ulii) & (BITCT - 2);
	    ukk = (ulii >> ujj) & 3;
	    sample_idx = uii + (ujj / 2);
	    dxx = g_perm_pheno_d2[sample_idx];
	    if (ukk == 1) {
	      qt_g_prod += dxx;
	      if (qt_means_or_lin) {
		qt_het_sum += dxx;
		qt_het_ssq += dxx * dxx;
	      }
	    } else if (ukk == 3) {
	      qt_g_prod += 2 * dxx;
	      if (qt_means_or_lin) {
		qt_homrar_sum += dxx;
		qt_homrar_ssq += dxx * dxx;
	      }
	    } else {
	      qt_sum -= dxx;
	      qt_ssq -= dxx * dxx;
	    }
	    ulii &= ~((3 * ONELU) << ujj);
	  }
	  uii += BITCT2;
	} while (uii < pheno_nm_ct);
	nanal_recip = 1.0 / ((double)nanal);
	qt_mean = qt_sum * nanal_recip;
	geno_mean = ((double)geno_sum) * nanal_recip;
	dxx = 1.0 / ((double)(nanal - 1));
	qt_var = (qt_ssq - qt_sum * qt_mean) * dxx;
	geno_var = (((double)geno_ssq) - geno_sum * geno_mean) * dxx;
	qt_g_prod_centered = qt_g_prod - qt_sum * geno_mean;
	qt_g_covar = qt_g_prod_centered * dxx;

	dxx = 1.0 / geno_var;
	beta = qt_g_covar * dxx;
	vbeta_sqrt = sqrt((qt_var * dxx - beta * beta) / ((double)(nanal - 2)));
	tstat = beta / vbeta_sqrt;
	if (fill_orig_chiabs) {
	  g_orig_chisq[marker_idx + marker_bidx] = tstat;
	  if (tcnt) {
	    tcnt[marker_idx + marker_bidx] = (nanal > 2)? (nanal - 2) : 0;
	  }
	}
	if (do_lin) {
	  // Square of Lin statistic:
	  //   \frac{(\sum_{i=1}^nU_{ji})^2}{\sum_{i=1}^nU_{ji}^2}
	  // where U_{ji} = (Y_i - \bar{Y_{\dot}})(X_{ji} - \bar{X_{j\dot}}),
	  // Y_{\dot}s are phenotypes, and X_{\dot\dot}s are genotypes.
	  //
	  // We evaluate the denominator by separating the sum into three
	  // components (one for each possible genotype value), each of which
	  // can be computed from the partial sums/sums-of-squares we already
	  // have.
	  g_orig_linsq[marker_idx + marker_bidx] = qt_g_prod_centered * qt_g_prod_centered / (geno_mean * geno_mean * (qt_ssq - 2 * qt_sum + qt_mean * qt_sum) + (1 - 2 * geno_mean) * (qt_het_ssq - 2 * qt_het_sum * qt_mean + qt_mean * qt_mean * ((intptr_t)het_ct)) + (4 - 4 * geno_mean) * (qt_homrar_ssq - 2 * qt_homrar_sum * qt_mean + qt_mean * qt_mean * ((intptr_t)homrar_ct)));
	}
	if (nanal > 1) {
	  tp = calc_tprob(tstat, nanal - 2);
	  rsq = (qt_g_covar * qt_g_covar) / (qt_var * geno_var);
	  if (mperm_save & MPERM_DUMP_ALL) {
	    if (!do_lin) {
	      if (tp >= 0) {
		dtoa_gx(tstat * tstat, '\0', &(numbuf[1]));
		fputs(numbuf, outfile_msa);
	      } else {
		fputs(" NA", outfile_msa);
	      }
	    } else {
	      dxx = g_orig_linsq[marker_idx + marker_bidx];
	      if ((nanal > 2) && realnum(dxx)) {
		dtoa_gx(dxx, '\0', &(numbuf[1]));
		fputs(numbuf, outfile_msa);
	      } else {
		fputs(" NA", outfile_msa);
	      }
	    }
	  }
	  if ((pfilter != 2.0) && ((tp > pfilter) || (tp == -9))) {
	    continue;
	  }
	  if (!realnum(beta)) {
	    wptr = memcpya(wptr, "        NA         NA         NA ", 33);
	  } else {
	    wptr = dtoa_g_wxp4x(beta, 10, ' ', wptr);
	    wptr = dtoa_g_wxp4x(vbeta_sqrt, 10, ' ', wptr);
	    wptr = dtoa_g_wxp4x(rsq, 10, ' ', wptr);
	  }
	  if (tp >= 0) {
	    wptr = dtoa_g_wxp4x(tstat, 8, ' ', wptr);
	    wptr = dtoa_g_wxp4(MAXV(tp, output_min_p), 12, wptr);
	  } else {
	    wptr = memcpya(wptr, "      NA           NA", 21);
	  }
	  if (do_lin && (nanal > 2)) {
	    dxx = g_orig_linsq[marker_idx + marker_bidx];
	    if (realnum(dxx)) {
	      *wptr++ = ' ';
	      dxx = sqrt(dxx);
	      wptr = dtoa_g_wxp4x(dxx, 12, ' ', wptr);
	      dxx = calc_tprob(dxx, nanal - 2);
	      wptr = dtoa_g_wxp4(MAXV(dxx, output_min_p), 12, wptr);
	    } else {
	      wptr = memcpya(wptr, "           NA           NA", 26);
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	} else if (pfilter != 2.0) {
	  continue;
	} else {
	  wptr = memcpya(wptr, "        NA         NA         NA       NA           NA ", 55);
	  if (mperm_save & MPERM_DUMP_ALL) {
	    fputs(" NA", outfile_msa);
	  }
	  if (do_lin) {
	    wptr = memcpya(wptr, "          NA           NA ", 26);
	  }
	  *wptr++ = '\n';
	}
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto qassoc_ret_WRITE_FAIL;
	}
	if (qt_means) {
	  wptr_restart = &(g_textbuf[2 + chrom_name_len + plink_maxsnp]);
	  wptr = memcpya(wptr_restart, "  GENO ", 7);
	  a1ptr = marker_allele_ptrs[2 * marker_uidx2];
	  a2ptr = marker_allele_ptrs[2 * marker_uidx2 + 1];
	  uii = strlen(a1ptr);
	  ujj = strlen(a2ptr);
	  if (uii < 4) {
	    wptr = memseta(wptr, 32, 7 - 2 * uii);
	  }
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  fputs(a1ptr, outfile_qtm);
	  putc_unlocked('/', outfile_qtm);
	  fputs(a1ptr, outfile_qtm);
	  putc_unlocked(' ', outfile_qtm);
	  if (uii + ujj < 7) {
	    fwrite(spacebuf, 1, 7 - uii - ujj, outfile_qtm);
	  }
	  fputs(a1ptr, outfile_qtm);
	  putc_unlocked('/', outfile_qtm);
	  fputs(a2ptr, outfile_qtm);
	  putc_unlocked(' ', outfile_qtm);
	  if (ujj < 4) {
	    fwrite(spacebuf, 1, 7 - 2 * ujj, outfile_qtm);
	  }
          fputs(a2ptr, outfile_qtm);
	  putc_unlocked('/', outfile_qtm);
          fputs(a2ptr, outfile_qtm);
	  putc_unlocked('\n', outfile_qtm);
	  wptr = memcpya(wptr_restart, "COUNTS ", 7);
	  wptr = uint32toa_w8x(homrar_ct, ' ', wptr);
	  wptr = uint32toa_w8x(het_ct, ' ', wptr);
	  wptr = uint32toa_w8x(homcom_ct, '\n', wptr);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "  FREQ ", 7);
	  wptr = dtoa_g_wxp4x(nanal_recip * ((intptr_t)homrar_ct), 8, ' ', wptr);
	  wptr = dtoa_g_wxp4x(nanal_recip * ((intptr_t)het_ct), 8, ' ', wptr);
	  wptr = dtoa_g_wxp4x(nanal_recip * ((intptr_t)homcom_ct), 8, '\n', wptr);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "  MEAN ", 7);
	  qt_homcom_sum = qt_sum - qt_homrar_sum - qt_het_sum;
	  if (homrar_ct) {
	    x11 = qt_homrar_sum / ((double)homrar_ct);
	    wptr = dtoa_g_wxp4(x11, 8, wptr);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (het_ct) {
	    x12 = qt_het_sum / ((double)het_ct);
	    wptr = dtoa_g_wxp4(x12, 8, wptr);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (homcom_ct) {
	    x22 = qt_homcom_sum / ((double)homcom_ct);
	    wptr = dtoa_g_wxp4(x22, 8, wptr);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = '\n';
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "    SD ", 7);
	  if (homrar_ct > 1) {
            dxx = sqrt((qt_homrar_ssq - qt_homrar_sum * x11) / ((double)((intptr_t)homrar_ct - 1)));
	    wptr = dtoa_g_wxp4(dxx, 8, wptr);
	  } else if (homrar_ct == 1) {
	    wptr = memcpya(wptr, "       0", 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (het_ct > 1) {
            dxx = sqrt((qt_het_ssq - qt_het_sum * x12) / ((double)((intptr_t)het_ct - 1)));
	    wptr = dtoa_g_wxp4(dxx, 8, wptr);
	  } else if (het_ct == 1) {
	    wptr = memcpya(wptr, "       0", 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (homcom_ct > 1) {
            dxx = sqrt((qt_ssq - qt_het_ssq - qt_homrar_ssq - qt_homcom_sum * x22) / ((double)((intptr_t)homcom_ct - 1)));
	    wptr = dtoa_g_wxp4(dxx, 8, wptr);
	  } else if (homcom_ct == 1) {
	    wptr = memcpya(wptr, "       0", 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = '\n';
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
      }
    }
    if (do_perms_nst) {
      is_last_block = (marker_idx + block_size >= marker_unstopped_ct);
      g_block_diff = block_size - g_qblock_start;
      ulii = 0;
      if (perm_maxt_nst) {
	g_maxt_block_base = marker_idx;
	// don't actually use maxt_cur_extreme_stat here?...
	if (!do_lin) {
	  if (spawn_threads2(threads, &qassoc_maxt_thread, max_thread_ct, is_last_block)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_maxt_thread((void*)ulii);
	} else {
	  if (spawn_threads2(threads, &qassoc_maxt_lin_thread, max_thread_ct, is_last_block)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_maxt_lin_thread((void*)ulii);
	}
        join_threads2(threads, max_thread_ct, is_last_block);
	ukk = g_block_diff / CACHELINE_DBL;
	if (ukk > max_thread_ct) {
	  ukk = max_thread_ct;
	} else if (!ukk) {
	  ukk = 1;
	}
	ulii = round_up_pow2(g_perm_vec_ct, CACHELINE_DBL);
	for (uii = 0; uii < ukk; uii++) {
	  ooptr = &(g_maxt_thread_results[uii * ulii]);
	  for (ujj = g_perms_done - g_perm_vec_ct; ujj < g_perms_done; ujj++) {
	    dxx = *ooptr++;
	    if (dxx > g_maxt_extreme_stat[ujj]) {
	      g_maxt_extreme_stat[ujj] = dxx;
	    }
	  }
	}
      } else {
	if (!do_lin) {
	  if (spawn_threads2(threads, &qassoc_adapt_thread, max_thread_ct, is_last_block)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_adapt_thread((void*)ulii);
	} else {
	  if (spawn_threads2(threads, &qassoc_adapt_lin_thread, max_thread_ct, is_last_block)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_adapt_lin_thread((void*)ulii);
	}
        join_threads2(threads, max_thread_ct, is_last_block);
      }
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_unstopped_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = (marker_idx * 100LLU) / marker_unstopped_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = (((uint64_t)pct + 1LLU) * marker_unstopped_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (qt_means) {
      LOGPRINTFWW("QT means report saved to %s.means .\n", outname);
      if (fclose_null(&outfile_qtm)) {
	goto qassoc_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto qassoc_ret_WRITE_FAIL;
    }
    if (!is_set_test) {
      if (do_perms_nst) {
	bigstack_reset(g_perm_vecstd);
      }
      if (mtest_adjust) {
	if (do_lin) {
	  for (uii = 0; uii < marker_ct; uii++) {
	    g_orig_chisq[uii] = sqrt(g_orig_linsq[uii]);
	  }
	}
	retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, g_orig_chisq, pfilter, output_min_p, mtest_adjust, 0, adjust_lambda, tcnt, nullptr);
	if (retval) {
	  goto qassoc_ret_1;
	}
      }
      if (mperm_save & MPERM_DUMP_ALL) {
	if (putc_checked('\n', outfile_msa)) {
	  goto qassoc_ret_WRITE_FAIL;
	}
      }
    } else {
      retval = qassoc_set_test(threads, bedfile, bed_offset, outname, outname_end, model_modifier, model_mperm_val, pfilter, output_min_p, mtest_adjust, unfiltered_marker_ct, marker_exclude_orig, marker_ct_orig, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, sex_male, apip, pheno_nm_ct, pheno_nm, founder_pnm, sample_include2, sample_male_include2, ld_ignore_x, hh_exists, hh_or_mt_exists, perm_batch_size, sip, tcnt, loadbuf_raw);
      if (retval) {
	goto qassoc_ret_1;
      }
    }
  }
  if (do_perms_nst) {
    if (mperm_save & MPERM_DUMP_ALL) {
      if (perm_pass_idx) {
	putc_unlocked(' ', stdout);
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = g_perm_vec_ct;
      ujj = 1 + g_perms_done - ulii;
      wptr = g_textbuf;
      a1ptr = &(g_textbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32toa(uii + ujj, wptr);
	ooptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = ooptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = dtoa_g(dxx, wptr);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= a1ptr) {
	    if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	      goto qassoc_ret_WRITE_FAIL;
	    }
	    wptr = g_textbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	goto qassoc_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    bigstack_reset(g_perm_vecstd);
    if (g_perms_done < perms_total) {
      if (perm_adapt_nst) {
	marker_unstopped_ct = marker_ct - popcount01_longs((uintptr_t*)g_perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
	if (!marker_unstopped_ct) {
	  goto qassoc_adapt_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done != 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto qassoc_more_perms;
    }
    if (perm_adapt_nst) {
    qassoc_adapt_perm_count:
      g_perms_done = 0;
      for (uii = 0; uii < marker_ct; uii++) {
	if (g_perm_attempt_ct[uii] > g_perms_done) {
	  g_perms_done = g_perm_attempt_ct[uii];
	  if (g_perms_done == perms_total) {
	    break;
	  }
	}
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("%u %s permutation%s complete.\n", g_perms_done, perm_maxt_nst? "max(T)" : "(adaptive)", (g_perms_done != 1)? "s" : "");

    if (perm_adapt_nst) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	memcpy(outname_end, ".mperm.dump.best", 17);
	LOGPRINTFWW("Dumping best permutation squared %sstats to %s .\n", do_lin? "Lin " : "Wald t-", outname);
	if (fopen_checked(outname, "w", &outfile)) {
	  goto qassoc_ret_OPEN_FAIL;
	}
	dxx = 0;
	if (!do_lin) {
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    if (fabs(g_orig_chisq[marker_idx]) > dxx) {
	      dxx = fabs(g_orig_chisq[marker_idx]);
	    }
	  }
	  dxx = dxx * dxx;
	} else {
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    if (g_orig_linsq[marker_idx] > dxx) {
	      dxx = g_orig_linsq[marker_idx];
	    }
	  }
	}
        memcpy(g_textbuf, "0 ", 2);
	wptr = dtoa_gx(dxx, '\n', &(g_textbuf[2]));
	if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile)) {
	  goto qassoc_ret_WRITE_FAIL;
	}
	for (uii = 0; uii < perms_total; uii++) {
	  wptr = uint32toa_x(uii + 1, ' ', g_textbuf);
	  wptr = dtoa_gx(g_maxt_extreme_stat[uii], '\n', wptr);
	  if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto qassoc_ret_WRITE_FAIL;
	}
	memcpy(outname_end, ".qassoc", 7); // deliberately not null-terminated
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto qassoc_ret_OPEN_FAIL;
    }
    if (perm_adapt_nst) {
      sprintf(g_textbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(g_textbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(g_maxt_extreme_stat, &(g_maxt_extreme_stat[perms_total]));
#else
      qsort(g_maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
    }
    // (debugging)
    // if (perm_maxt) {
    //   printf("extreme stats: %g %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[(perms_total - 1) / 2], g_maxt_extreme_stat[perms_total - 1]);
    // }
    fprintf(outfile, g_textbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, uii, g_textbuf));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	if (perm_adapt_nst) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) * dxx;
	}
        if (pval <= pfilter) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	  wptr = &(wptr_start[1 + plink_maxsnp]);
	  if (perm_adapt_nst && (!g_perm_attempt_ct[marker_idx])) {
	    // invalid
            wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    if (!perm_count) {
	      wptr = dtoa_g_wxp4x(pval, 12, ' ', wptr);
	    } else {
	      wptr = dtoa_g_wxp4x(((double)g_perm_2success_ct[marker_idx]) * 0.5, 12, ' ', wptr);
	    }
	    if (perm_adapt_nst) {
	      wptr = memseta(wptr, 32, 2);
	      wptr = uint32toa_w10(g_perm_attempt_ct[marker_idx], wptr);
	    } else {
	      // maximum chisq
	      // N.B. numbers in maxt_extreme_stat[] have been pre-squared
	      // while orig_chisq[] has not been
	      if (do_lin) {
		dzz = g_orig_linsq[marker_idx];
	      } else {
		dzz = g_orig_chisq[marker_idx] * g_orig_chisq[marker_idx];
	      }
	      dzz = (int32_t)(perms_total - doublearr_greater_than(g_maxt_extreme_stat, perms_total, dzz - EPSILON) + 1);
	      if (!perm_count) {
		wptr = dtoa_g_wxp4(dzz * dyy, 12, wptr);
	      } else {
		wptr = dtoa_g_wxp4(dzz - 1, 12, wptr);
	      }
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto qassoc_loop_end;
	}
	marker_uidx++;
        next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      }
    }
  qassoc_loop_end:
    if (fclose_null(&outfile)) {
      goto qassoc_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
  }

  while (0) {
  qassoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  qassoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  qassoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  qassoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  qassoc_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  qassoc_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 qassoc_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_qtm);
  fclose_cond(outfile_msa);
  return retval;
}

int32_t gxe_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, uintptr_t* sample_exclude, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* gxe_covar_nm, uintptr_t* gxe_covar_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t covar_nm_ct = popcount_longs(gxe_covar_nm, sample_ctl);
  uintptr_t covar_nm_ctl = BITCT_TO_WORDCT(covar_nm_ct);
  // gxe_covar_c has opposite truth value from ->bcovar in PLINK 1.07 gxe.cpp;
  // see lines 50-58 in gxe.cpp
  uintptr_t group2_size = popcount_longs(gxe_covar_c, sample_ctl);
  uintptr_t group1_size = covar_nm_ct - group2_size;
  uintptr_t male_ct = 0;
  uintptr_t male_ctl = 0;
  uintptr_t group1_size_male = 0;
  uintptr_t group2_size_male = 0;
  uintptr_t marker_uidx = 0;
  uintptr_t final_mask = 0;
  uintptr_t* sample_include2 = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uintptr_t* sample_male_all_include2 = nullptr;
  uintptr_t* group1_include2 = nullptr;
  uintptr_t* group2_include2 = nullptr;
  uintptr_t* group1_male_include2 = nullptr;
  uintptr_t* group2_male_include2 = nullptr;
  uintptr_t* covar_nm_raw = nullptr;
  uintptr_t* covar_nm_male_raw = nullptr;
  uintptr_t* cur_sample_i2 = nullptr;
  uintptr_t* cur_sample_male_i2 = nullptr;
  uintptr_t* cur_group1_i2 = nullptr;
  uintptr_t* cur_group2_i2 = nullptr;
  uintptr_t* cur_covar_nm_raw = nullptr;
  double* pheno_d_collapsed = nullptr;
  double* pheno_d_male_collapsed = nullptr;
  double* cur_pheno_d = nullptr;
  char* wptr_start = nullptr;
  uintptr_t cur_sample_ct = 0;
  uintptr_t cur_sample_ctv2 = 0;
  uintptr_t cur_group1_size = 0;
  uintptr_t cur_group2_size = 0;
  uint32_t y_exists = (chrom_info_ptr->xymt_codes[Y_OFFSET] != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[Y_OFFSET]);
  uint32_t mt_exists = (chrom_info_ptr->xymt_codes[MT_OFFSET] != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[MT_OFFSET]);
  uint32_t skip_y = 0;
  double pheno_sum_g1 = 0;
  double pheno_ssq_g1 = 0;
  double pheno_sum_g2 = 0;
  double pheno_ssq_g2 = 0;
  double pheno_sum_male_g1 = 0;
  double pheno_ssq_male_g1 = 0;
  double pheno_sum_male_g2 = 0;
  double pheno_ssq_male_g2 = 0;
  double base_pheno_sum_g1 = 0;
  double base_pheno_ssq_g1 = 0;
  double base_pheno_sum_g2 = 0;
  double base_pheno_ssq_g2 = 0;
  int32_t retval = 0;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_ptr;
  uintptr_t* cgr_ptr;
  char* wptr;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uintptr_t loop_end;
  uintptr_t marker_idx;
  uintptr_t sample_uidx;
  uintptr_t sample_uidx_stop;
  uintptr_t sample_idx;
  uintptr_t sample_idx2;
  uintptr_t sample_idx2_offset;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uintptr_t ulnn;
  double dxx;
  double qt_sum1;
  double qt_ssq1;
  double qt_g_prod1;
  double nanal_recip1;
  double nanal_m1_recip1;
  double geno_mean1;
  double g_var1;
  double qt_var1;
  double qt_g_covar1;
  double beta1;
  double vbeta1;

  double qt_sum2;
  double qt_ssq2;
  double qt_g_prod2;
  double nanal_recip2;
  double nanal_m1_recip2;
  double geno_mean2;
  double g_var2;
  double qt_var2;
  double qt_g_covar2;
  double beta2;
  double vbeta2;

  double zval;

  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t min_ploidy_1;
  uint32_t pct;

  uint32_t missing_ct1;
  uint32_t het_ct1;
  uint32_t homcom_ct1;
  uint32_t homrar_ct1;
  uint32_t nanal1;
  uint32_t geno_sum1;
  uint32_t geno_ssq1;

  uint32_t missing_ct2;
  uint32_t het_ct2;
  uint32_t homcom_ct2;
  uint32_t homrar_ct2;
  uint32_t nanal2;
  uint32_t geno_sum2;
  uint32_t geno_ssq2;

  if (group1_size < 3) {
    logerrprint("Error: First --gxe group has fewer than three members.\n");
    goto gxe_assoc_ret_INVALID_CMDLINE;
  } else if (group2_size < 3) {
    logerrprint("Error: Second --gxe group has fewer than three members.\n");
    goto gxe_assoc_ret_INVALID_CMDLINE;
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl * 2, &loadbuf_raw) ||
      bigstack_alloc_ul(covar_nm_ctl * 2, &loadbuf) ||
      bigstack_calloc_ul(unfiltered_sample_ctl, &covar_nm_raw) ||
      bigstack_alloc_d(covar_nm_ct, &pheno_d_collapsed)) {
    goto gxe_assoc_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl * 2 - 1] = 0;

  sample_uidx = 0;
  sample_idx = 0;
  sample_idx2 = 0;
  do {
    sample_uidx = next_unset_ul_unsafe(sample_exclude, sample_uidx);
    sample_uidx_stop = next_set_ul(sample_exclude, sample_uidx, unfiltered_sample_ct);
    do {
      if (IS_SET(gxe_covar_nm, sample_idx)) {
        SET_BIT(sample_uidx, covar_nm_raw);
        dxx = pheno_d[sample_uidx];
        if (IS_SET(gxe_covar_c, sample_idx)) {
	  pheno_sum_g2 += dxx;
	  pheno_ssq_g2 += dxx * dxx;
	} else {
	  pheno_sum_g1 += dxx;
	  pheno_ssq_g1 += dxx * dxx;
	}
	pheno_d_collapsed[sample_idx2++] = dxx;
      }
      sample_idx++;
    } while (++sample_uidx < sample_uidx_stop);
  } while (sample_idx < sample_ct);

  if (bigstack_alloc_ul(covar_nm_ctl * 2, &group1_include2) ||
      bigstack_calloc_ul(covar_nm_ctl * 2, &group2_include2)) {
    goto gxe_assoc_ret_NOMEM;
  }
  fill_quatervec_55(covar_nm_ct, group1_include2);
  sample_idx = 0;
  sample_idx2 = 0;
  do {
    sample_idx = next_set_ul_unsafe(gxe_covar_nm, sample_idx);
    sample_uidx_stop = next_unset_ul(gxe_covar_nm, sample_idx, sample_ct);
    do {
      if (IS_SET(gxe_covar_c, sample_idx)) {
	SET_BIT_DBL(sample_idx2, group2_include2);
      }
      sample_idx2++;
    } while (++sample_idx < sample_uidx_stop);
  } while (sample_idx2 < covar_nm_ct);
  bitvec_andnot(group2_include2, covar_nm_ctl * 2, group1_include2);

  hh_or_mt_exists |= mt_exists * NXMHH_EXISTS;
  if ((hh_or_mt_exists & NXMHH_EXISTS) || y_exists) {
    if (bigstack_alloc_ul(covar_nm_ctl * 2, &sample_include2)) {
      goto gxe_assoc_ret_NOMEM;
    }
    fill_quatervec_55(covar_nm_ct, sample_include2);
  }
  if ((hh_or_mt_exists & XMHH_EXISTS) || y_exists) {
    if (bigstack_calloc_ul(covar_nm_ctl * 2, &sample_male_include2)) {
      goto gxe_assoc_ret_NOMEM;
    }
    sample_uidx = 0;
    sample_idx = 0;
    sample_idx2 = 0;
    do {
      sample_uidx = next_unset_ul_unsafe(sample_exclude, sample_uidx);
      sample_uidx_stop = next_set_ul(sample_exclude, sample_uidx, unfiltered_sample_ct);
      do {
        if (IS_SET(gxe_covar_nm, sample_idx)) {
          if (IS_SET(sex_male, sample_uidx)) {
	    SET_BIT_DBL(sample_idx2, sample_male_include2);
	    male_ct++;
	  }
	  sample_idx2++;
	}
	sample_idx++;
      } while (++sample_uidx < sample_uidx_stop);
    } while (sample_idx < sample_ct);
    male_ctl = BITCT_TO_WORDCT(male_ct);
    if (y_exists) {
      group1_size_male = popcount_longs_exclude(sample_male_include2, group2_include2, covar_nm_ctl * 2);
      group2_size_male = male_ct - group1_size_male;
      if ((group1_size_male < 3) || (group2_size_male < 3)) {
        logerrprint("Warning: Skipping Y chromosome for --gxe since a group has less than 3 males.\n");
	skip_y = 1;
      }
      // currently still need to initialize covar_nm_male_raw even on skip_y
      if (bigstack_alloc_ul(male_ctl * 2, &sample_male_all_include2) ||
          bigstack_alloc_ul(male_ctl * 2, &group1_male_include2) ||
	  bigstack_calloc_ul(male_ctl * 2, &group2_male_include2) ||
	  bigstack_alloc_d(male_ct, &pheno_d_male_collapsed) ||
	  bigstack_alloc_ul(unfiltered_sample_ctl, &covar_nm_male_raw)) {
	goto gxe_assoc_ret_NOMEM;
      }
      fill_quatervec_55(male_ct, sample_male_all_include2);
      fill_quatervec_55(male_ct, group1_male_include2);
      sample_idx = 0;
      for (sample_idx2 = 0; sample_idx2 < covar_nm_ct; sample_idx2++) {
	if (IS_SET_DBL(sample_male_include2, sample_idx2)) {
	  dxx = pheno_d_collapsed[sample_idx2];
	  if (IS_SET_DBL(group2_include2, sample_idx2)) {
	    SET_BIT_DBL(sample_idx, group2_male_include2);
	    pheno_sum_male_g2 += dxx;
	    pheno_ssq_male_g2 += dxx * dxx;
	  } else {
	    pheno_sum_male_g1 += dxx;
            pheno_ssq_male_g1 += dxx * dxx;
	  }
	  pheno_d_male_collapsed[sample_idx++] = dxx;
	}
      }
      bitvec_andnot(group2_male_include2, male_ctl * 2, group1_male_include2);
      for (ulii = 0; ulii < unfiltered_sample_ctl; ulii++) {
	covar_nm_male_raw[ulii] = covar_nm_raw[ulii] & sex_male[ulii];
      }
    }
  }

  memcpy(outname_end, ".qassoc.gxe", 12);
  if (fopen_checked(outname, "w", &outfile)) {
    goto gxe_assoc_ret_OPEN_FAIL;
  }
  if (haploid_chrom_present(chrom_info_ptr) || mt_exists) {
    logerrprint("Warning: --gxe doesn't currently handle X/Y/MT/haploid variants properly.\n");
  }
  LOGPRINTFWW5("Writing --gxe report to %s ... ", outname);
  fputs("0%", stdout);
  fflush(stdout);
  sprintf(g_textbuf, " CHR %%%us   NMISS1      BETA1        SE1   NMISS2      BETA2        SE2    Z_GXE        P_GXE \n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");

  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto gxe_assoc_ret_READ_FAIL;
  }
  // exploit overflow for initialization
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = 0;
  marker_idx = 0;
  chrom_end = 0;
  for (pct = 1; pct <= 100; pct++) {
    loop_end = (((uint64_t)pct) * marker_ct) / 100;
    for (; marker_idx < loop_end; marker_idx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto gxe_assoc_ret_READ_FAIL;
	}
      }
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &min_ploidy_1);
	min_ploidy_1 |= is_mt;
	if (!is_y) {
	  cur_sample_ct = covar_nm_ct;
	  cur_group1_size = group1_size;
          cur_group2_size = group2_size;
	  base_pheno_sum_g1 = pheno_sum_g1;
	  base_pheno_ssq_g1 = pheno_ssq_g1;
          base_pheno_sum_g2 = pheno_sum_g2;
          base_pheno_ssq_g2 = pheno_ssq_g2;
          cur_sample_i2 = sample_include2;
          cur_sample_male_i2 = sample_male_include2;
	  cur_group1_i2 = group1_include2;
          cur_group2_i2 = group2_include2;
          cur_pheno_d = pheno_d_collapsed;
	  cur_covar_nm_raw = covar_nm_raw;
	} else {
	  cur_sample_ct = male_ct;
	  cur_group1_size = group1_size_male;
          cur_group2_size = group2_size_male;
          base_pheno_sum_g1 = pheno_sum_male_g1;
	  base_pheno_ssq_g1 = pheno_ssq_male_g1;
          base_pheno_sum_g2 = pheno_sum_male_g2;
	  base_pheno_ssq_g2 = pheno_ssq_male_g2;
          cur_sample_i2 = sample_male_all_include2;
          cur_sample_male_i2 = sample_male_all_include2;
          cur_group1_i2 = group1_male_include2;
          cur_group2_i2 = group2_male_include2;
          cur_pheno_d = pheno_d_male_collapsed;
	  cur_covar_nm_raw = covar_nm_male_raw;
	}
	wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], g_textbuf));
	*wptr_start++ = ' ';
	cur_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(cur_sample_ct);
        loadbuf[cur_sample_ctv2 - 1] = 0;
	final_mask = get_final_mask(cur_sample_ct);
      }

      if (load_and_collapse_incl(unfiltered_sample_ct, cur_sample_ct, cur_covar_nm_raw, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf)) {
	goto gxe_assoc_ret_READ_FAIL;
      }
      if (is_y && skip_y) {
	marker_uidx++;
	continue;
      }
      if (min_ploidy_1) {
	haploid_fix(hh_or_mt_exists, cur_sample_i2, cur_sample_male_i2, cur_sample_ct, is_x, is_y, (unsigned char*)loadbuf);
      }

      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
      *wptr++ = ' ';

      // We are interested in the following quantities:
      //   qt_var{1,2}: (qt_ssq - (qt_sum^2 / N)) / (N-1)
      //   g_var{1,2}: (geno_ssq - (geno_sum^2 / N)) / (N-1)
      //   qt_g_covar{1,2}: (qt_g_prod - ((qt_sum * geno_sum) / N)) / (N-1)

      single_marker_cc_3freqs(cur_sample_ctv2, loadbuf, cur_group1_i2, cur_group2_i2, &homcom_ct1, &het_ct1, &missing_ct1, &homcom_ct2, &het_ct2, &missing_ct2);
      nanal1 = ((uint32_t)cur_group1_size) - missing_ct1;
      nanal2 = ((uint32_t)cur_group2_size) - missing_ct2;
      homrar_ct1 = nanal1 - (het_ct1 + homcom_ct1);
      homrar_ct2 = nanal2 - (het_ct2 + homcom_ct2);
      geno_sum1 = 2 * homrar_ct1 + het_ct1;
      geno_sum2 = 2 * homrar_ct2 + het_ct2;
      geno_ssq1 = 4 * homrar_ct1 + het_ct1;
      geno_ssq2 = 4 * homrar_ct2 + het_ct2;

      if ((nanal1 > 2) && (nanal2 > 2)) {
	nanal_recip1 = 1.0 / ((int32_t)nanal1);
	nanal_recip2 = 1.0 / ((int32_t)nanal2);
	nanal_m1_recip1 = 1.0 / ((int32_t)(nanal1 - 1));
	nanal_m1_recip2 = 1.0 / ((int32_t)(nanal2 - 1));
	geno_mean1 = geno_sum1 * nanal_recip1;
	g_var1 = (geno_ssq1 - geno_sum1 * geno_mean1) * nanal_m1_recip1;
        geno_mean2 = geno_sum2 * nanal_recip2;
        g_var2 = (geno_ssq2 - geno_sum2 * geno_mean2) * nanal_m1_recip2;
	if ((g_var1 == 0) || (g_var2 == 0)) {
	  goto gxe_assoc_nan_line;
	}
	qt_sum1 = base_pheno_sum_g1;
	qt_ssq1 = base_pheno_ssq_g1;
	qt_sum2 = base_pheno_sum_g2;
	qt_ssq2 = base_pheno_ssq_g2;
	qt_g_prod1 = 0;
	qt_g_prod2 = 0;
	sample_idx2_offset = 0;
	loadbuf_ptr = loadbuf;
	cgr_ptr = cur_group2_i2;
	do {
	  ulmm = ~(*loadbuf_ptr++);
	  if (sample_idx2_offset + BITCT2 > cur_sample_ct) {
	    ulmm &= (ONELU << ((cur_sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
	  }
	  if (ulmm) {
	    ulnn = (*cgr_ptr) * 3;
            ulii = ulmm & (~ulnn);
            while (ulii) {
	      uljj = CTZLU(ulii) & (BITCT - 2);
	      ulkk = (ulii >> uljj) & 3;
	      sample_idx2 = sample_idx2_offset + (uljj / 2);
	      dxx = cur_pheno_d[sample_idx2];
	      if (ulkk == 1) {
		// het
		qt_g_prod1 += dxx;
	      } else if (ulkk == 3) {
		// hom rare
		qt_g_prod1 += 2 * dxx;
	      } else {
		// missing
		qt_sum1 -= dxx;
		qt_ssq1 -= dxx * dxx;
	      }
	      ulii &= ~((3 * ONELU) << uljj);
	    }
	    ulii = ulmm & ulnn;
	    while (ulii) {
	      uljj = CTZLU(ulii) & (BITCT - 2);
	      ulkk = (ulii >> uljj) & 3;
	      sample_idx2 = sample_idx2_offset + (uljj / 2);
	      dxx = cur_pheno_d[sample_idx2];
	      if (ulkk == 1) {
		qt_g_prod2 += dxx;
	      } else if (ulkk == 3) {
		qt_g_prod2 += 2 * dxx;
	      } else {
		qt_sum2 -= dxx;
		qt_ssq2 -= dxx * dxx;
	      }
	      ulii &= ~((3 * ONELU) << uljj);
	    }
	  }
	  cgr_ptr++;
	  sample_idx2_offset += BITCT2;
	} while (sample_idx2_offset < cur_sample_ct);
        qt_var1 = (qt_ssq1 - (qt_sum1 * qt_sum1 * nanal_recip1)) * nanal_m1_recip1;
        qt_var2 = (qt_ssq2 - (qt_sum2 * qt_sum2 * nanal_recip2)) * nanal_m1_recip2;
	qt_g_covar1 = (qt_g_prod1 - (qt_sum1 * geno_mean1)) * nanal_m1_recip1;
        qt_g_covar2 = (qt_g_prod2 - (qt_sum2 * geno_mean2)) * nanal_m1_recip2;
	beta1 = qt_g_covar1 / g_var1;
        beta2 = qt_g_covar2 / g_var2;
        vbeta1 = (qt_var1 / g_var1 - (qt_g_covar1 * qt_g_covar1) / (g_var1 * g_var1)) / ((double)(((int32_t)nanal1) - 2));

        vbeta2 = (qt_var2 / g_var2 - (qt_g_covar2 * qt_g_covar2) / (g_var2 * g_var2)) / ((double)(((int32_t)nanal2) - 2));
        if (vbeta1 + vbeta2 <= 0) {
	  goto gxe_assoc_nan_line;
	}
        zval = (beta1 - beta2) / sqrt(vbeta1 + vbeta2);
        wptr = uint32toa_w8x(nanal1, ' ', wptr);
        wptr = dtoa_g_wxp4x(beta1, 10, ' ', wptr);
        wptr = dtoa_g_wxp4x(sqrt(vbeta1), 10, ' ', wptr);
        wptr = uint32toa_w8x(nanal2, ' ', wptr);
        wptr = dtoa_g_wxp4x(beta2, 10, ' ', wptr);
        wptr = dtoa_g_wxp4x(sqrt(vbeta2), 10, ' ', wptr);
        wptr = dtoa_g_wxp4x(zval, 8, ' ', wptr);
	dxx = chiprob_p(zval * zval, 1);
        wptr = dtoa_g_wxp4x(MAXV(dxx, output_min_p), 12, '\n', wptr);
      } else {
      gxe_assoc_nan_line:
        wptr = memcpya(wptr, "      NA         NA         NA       NA         NA         NA       NA           NA\n", 84);
      }
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto gxe_assoc_ret_WRITE_FAIL;
      }
      marker_uidx++;
    }
    if (pct < 100) {
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  if (fclose_null(&outfile)) {
    goto gxe_assoc_ret_WRITE_FAIL;
  }
  if (pct >= 10) {
    putc_unlocked('\b', stdout);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");

  while (0) {
  gxe_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  gxe_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  gxe_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  gxe_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  gxe_assoc_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

void calc_git_missing(uint32_t pheno_nm_ct, uint32_t perm_vec_ct, uintptr_t* __restrict__ loadbuf, uint32_t* perm_vecst, uint32_t* thread_wkspace) {
  // Simplified calc_git() for when we only need to distinguish between missing
  // and nonmissing.
  // thread_wkspace[] is assumed to be zeroed out before this function is
  // called.
  uint32_t pheno_nm_ctl = BITCT_TO_WORDCT(pheno_nm_ct);
#ifdef __LP64__
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct128x4 = perm_ct128 * 4;
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  __m128i* permsv = (__m128i*)perm_vecst;
  __m128i* gitv[3];
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct32x4 = perm_ct32 * 4;
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uintptr_t* permsv = (uintptr_t*)perm_vecst;
  uintptr_t* gitv[3];
#endif
  uint32_t cur_ct;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  // 4- and 8-bit partial counts
  gitv[0] = &(((__m128i*)thread_wkspace)[8 * perm_ct128x4]);
  gitv[1] = &(((__m128i*)thread_wkspace)[9 * perm_ct128x4]);
  gitv[2] = (__m128i*)thread_wkspace;
#else
  gitv[0] = (uintptr_t*)(&(thread_wkspace[8 * perm_ct32x4]));
  gitv[1] = (uintptr_t*)(&(thread_wkspace[9 * perm_ct32x4]));
  gitv[2] = (uintptr_t*)thread_wkspace;
#endif
  cur_ct = 0;
  for (uii = 0; uii < pheno_nm_ctl; uii++) {
    ulii = *loadbuf++;
    if (uii + 1 == pheno_nm_ctl) {
      ujj = pheno_nm_ct & (BITCT2 - 1);
      if (ujj) {
	ulii &= (ONELU << ujj) - ONELU;
      }
    }
    while (ulii) {
      ujj = CTZLU(ulii);
      cur_ct++;
#ifdef __LP64__
      unroll_incr_1_4(&(permsv[ujj * perm_ct128]), gitv[0], perm_ct128);
      if (!(cur_ct % 15)) {
	unroll_zero_incr_4_8(gitv[0], gitv[1], perm_ct32);
	if (!(cur_ct % 255)) {
	  unroll_zero_incr_8_32(gitv[1], gitv[2], perm_ct16);
	}
      }
#else
      unroll_incr_1_4(&(permsv[ujj * perm_ct32]), gitv[0], perm_ct32);
      if (!(cur_ct % 15)) {
	unroll_zero_incr_4_8(gitv[0], gitv[1], perm_ct8);
	if (!(cur_ct % 255)) {
	  unroll_zero_incr_8_32(gitv[1], gitv[2], perm_ct4);
	}
      }
#endif
      ulii &= ulii - 1;
    }
#ifdef __LP64__
    permsv = &(permsv[BITCT * perm_ct128]);
#else
    permsv = &(permsv[BITCT * perm_ct32]);
#endif
  }
#ifdef __LP64__
  if (cur_ct % 15) {
    unroll_incr_4_8(gitv[0], gitv[1], perm_ct32);
  }
  if (cur_ct % 255) {
    unroll_incr_8_32(gitv[1], gitv[2], perm_ct16);
  }
#else
  if (cur_ct % 15) {
    unroll_incr_4_8(gitv[0], gitv[1], perm_ct8);
  }
  if (cur_ct % 255) {
    unroll_incr_8_32(gitv[1], gitv[2], perm_ct4);
  }
#endif
}

THREAD_RET_TYPE testmiss_adapt_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t pheno_nm_ctl = BITCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t pheno_nm_ctv = round_up_pow2(pheno_nm_ctl, VEC_WORDS);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t max_thread_ct = g_assoc_thread_ct;
  uint32_t pidx_offset = g_perms_done;
  uint32_t is_midp = g_fisher_midp;
  uint32_t first_adapt_check = g_first_adapt_check;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  // this can be cached since testmiss() computes all raw p-values before
  // starting permutation test
  double* __restrict__ orig_pvals = g_orig_pvals;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double stat_high = 0;
  double stat_low = 0;
  uint32_t missing_sum = 0;
  uint32_t nm_sum = 0;
  uint32_t* male_case_cts = nullptr;
  uintptr_t* __restrict__ loadbuf;
  uintptr_t* loadbuf_ptr;
  uint32_t* __restrict__ precomp_ui;
  uint32_t* __restrict__ missing_cts;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t is_y;
  uint32_t valid_obs_ct;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  uint32_t missing_case_ct;
  uint32_t case_ct;
  uint32_t uii;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  while (1) {
    if (g_block_diff <= max_thread_ct) {
      if (g_block_diff <= tidx) {
	goto testmiss_adapt_thread_skip_all;
      }
      marker_bidx = tidx;
      marker_bceil = tidx + 1;
    } else {
      marker_bidx = (((uint64_t)tidx) * g_block_diff) / max_thread_ct;
      marker_bceil = (((uint64_t)tidx + 1) * g_block_diff) / max_thread_ct;
    }
    is_y = 0;
    if (g_is_y) {
      valid_obs_ct = g_male_ct;
      if (valid_obs_ct != pheno_nm_ct) {
	is_y = 1; // if all male, can pretend as if this isn't Ychr
	male_case_cts = g_male_case_cts;
      }
    } else {
      valid_obs_ct = pheno_nm_ct;
    }
    loadbuf = g_loadbuf;
    precomp_ui = g_precomp_ui;
    missing_cts = g_missing_cts;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      next_adapt_check = first_adapt_check;
      gpui = &(precomp_ui[4 * marker_bidx]);
      if (is_y) {
	missing_sum = missing_cts[marker_idx];
	nm_sum = valid_obs_ct - missing_sum;
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      }
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      loadbuf_ptr = &(loadbuf[marker_bidx * pheno_nm_ctv]);
      for (pidx = 0; pidx < perm_vec_ct;) {
	missing_case_ct = popcount_longs_intersect(loadbuf_ptr, &(perm_vecs[pidx * pheno_nm_ctv]), pheno_nm_ctl);
	if (!is_y) {
	  if (missing_case_ct < gpui[0]) {
	    if (missing_case_ct < gpui[2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (missing_case_ct >= gpui[1]) {
	      if (missing_case_ct >= gpui[3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	} else {
	  case_ct = male_case_cts[pidx];
	  pval = fisher22(missing_case_ct, case_ct - missing_case_ct, missing_sum - missing_case_ct, nm_sum + missing_case_ct - case_ct, is_midp);
	  if (pval < stat_low) {
	    success_2incr += 2;
	  } else if (pval <= stat_high) {
	    success_2incr++;
	  }
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	    dyy = pval - dxx;
	    dzz = pval + dxx;
	    if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	      perm_adapt_stop[marker_idx] = 1;
	      perm_attempt_ct[marker_idx] = next_adapt_check;
	      break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  testmiss_adapt_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE testmiss_maxt_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uint32_t is_midp = g_fisher_midp;
  uint32_t max_thread_ct = g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl = BITCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t pheno_nm_ctv = round_up_pow2(pheno_nm_ctl, VEC_WORDS);
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct128 * 176]);
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  double* __restrict__ orig_pvals = g_orig_pvals;
  double* msa_ptr = nullptr;
  uint32_t* male_case_cts = nullptr;
  uint32_t* gpui = nullptr;
  double* gpd = nullptr;
  double stat_high = 0;
  double stat_low = 0;
  uint32_t case_ct = g_perm_case_ct;
  uint32_t cur_case_ct = case_ct;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_ptr;
  uint32_t* precomp_ui;
  uint32_t* __restrict__ missing_cts;
  double* __restrict__ precomp_d;
  uintptr_t pidx;
  uintptr_t marker_idx;
  double pval;
  uint32_t marker_bidx_start;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t is_y;
  uint32_t valid_obs_ct;
  uint32_t missing_sum;
  uint32_t nm_sum;
  uint32_t success_2incr;
  uint32_t missing_case_ct;
  uint32_t uii;
  uint32_t ujj;
  while (1) {
    if (g_block_diff <= max_thread_ct) {
      if (g_block_diff <= tidx) {
	goto testmiss_maxt_thread_skip_all;
      }
      marker_bidx_start = tidx;
      marker_bceil = tidx + 1;
    } else {
      marker_bidx_start = (((uint64_t)tidx) * g_block_diff) / max_thread_ct;
      marker_bceil = (((uint64_t)tidx + 1) * g_block_diff) / max_thread_ct;
    }
    marker_bidx = marker_bidx_start;
    marker_idx = g_maxt_block_base + marker_bidx_start;
    memcpy(results, &(g_maxt_extreme_stat[g_perms_done]), perm_vec_ct * sizeof(double));
    is_y = 0;
    if (g_is_y) {
      valid_obs_ct = g_male_ct;
      if (valid_obs_ct != pheno_nm_ct) {
	is_y = 1;
	male_case_cts = g_male_case_cts;
	precomp_ui = nullptr;
      }
    } else {
      valid_obs_ct = pheno_nm_ct;
    }
    loadbuf = g_loadbuf;
    missing_cts = g_missing_cts;
    precomp_d = g_precomp_d;
    if (g_mperm_save_all) {
      msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
      precomp_ui = nullptr;
    } else {
      precomp_ui = g_precomp_ui;
    }
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      missing_sum = missing_cts[marker_idx];
      nm_sum = valid_obs_ct - missing_sum;
      if (precomp_ui) {
	gpui = &(precomp_ui[6 * marker_bidx]);
	gpd = &(precomp_d[2 * marker_bidx]);
      } else {
	stat_high = orig_pvals[marker_idx] * (1.0 + EPSILON);
	stat_low = orig_pvals[marker_idx] * (1.0 - EPSILON);
      }
      loadbuf_ptr = &(loadbuf[marker_bidx * pheno_nm_ctv]);
      success_2incr = 0;
      fill_uint_zero(perm_ct128 * 176, thread_git_wkspace);
      calc_git_missing(pheno_nm_ct, perm_vec_ct, loadbuf_ptr, perm_vecst, thread_git_wkspace);
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	missing_case_ct = thread_git_wkspace[pidx];
	if (precomp_ui) {
	  if (missing_case_ct < gpui[0]) {
	    if (missing_case_ct < gpui[2]) {
	      success_2incr += 2;
	    } else {
	      success_2incr++;
	    }
	  } else {
	    if (missing_case_ct >= gpui[1]) {
	      if (missing_case_ct >= gpui[3]) {
		success_2incr += 2;
	      } else {
		success_2incr++;
	      }
	    }
	  }
	  ujj = gpui[4];
	  uii = (uint32_t)(missing_case_ct - ujj); // deliberate underflow
	  if (uii >= gpui[5]) {
	    pval = fisher22_tail_pval(ujj, missing_sum - ujj, case_ct - ujj, nm_sum + ujj - case_ct, gpui[5] - 1, gpd[0], gpd[1], is_midp, missing_case_ct);
	    if (results[pidx] > pval) {
	      results[pidx] = pval;
	    }
	  }
	} else {
	  if (is_y) {
	    cur_case_ct = male_case_cts[pidx];
	  }
	  pval = fisher22(missing_case_ct, missing_sum - missing_case_ct, cur_case_ct - missing_case_ct, nm_sum + missing_case_ct - cur_case_ct, is_midp);
	  if (pval < stat_low) {
	    success_2incr += 2;
	  } else if (pval <= stat_high) {
	    success_2incr++;
	  }
	  if (results[pidx] > pval) {
	    results[pidx] = pval;
	  }
	  if (msa_ptr) {
	    *msa_ptr++ = pval;
	  }
	}
      }
      perm_2success_ct[marker_idx++] += success_2incr;
    }
  testmiss_maxt_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

int32_t testmiss(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t testmiss_mperm_val, uint32_t testmiss_modifier, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_exists) {
  // Simple variant of model_assoc().
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t pheno_nm_ctl = BITCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t cur_sample_ctl = pheno_nm_ctl;
  uintptr_t pheno_nm_ctv = round_up_pow2(pheno_nm_ctl, VEC_WORDS);
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t marker_uidx = next_unset_unsafe(marker_exclude_orig, 0);
  double maxt_cur_extreme_stat = 0;
  FILE* outfile = nullptr;
  FILE* outfile_msa = nullptr;
  uintptr_t* sample_hh_include2 = nullptr;
  uintptr_t* sample_hh_male_include2 = nullptr;
  uintptr_t* pheno_male_nm2 = nullptr;
  uintptr_t* pheno_c_collapsed_male = nullptr;
  uintptr_t* sex_male_collapsed = nullptr;
  char* wptr_start = nullptr;
  char* tbuf2 = &(g_textbuf[MAXLINELEN]);
  uint32_t perm_adapt = testmiss_modifier & TESTMISS_PERM;
  uint32_t perm_maxt = testmiss_modifier & TESTMISS_MPERM;
  uint32_t perm_count = testmiss_modifier & TESTMISS_PERM_COUNT;
  uint32_t midp = testmiss_modifier & TESTMISS_MIDP;
  uint32_t do_perms = perm_adapt | perm_maxt;
  uint32_t perms_total = 0;
  uint32_t chrom_fo_idx = 0xffffffffU;
  uint32_t is_x = 0;
  // don't treat MT heterozygous call as missing
  uint32_t is_haploid = 0;
  uint32_t skip_y = 0;
  uint32_t cur_pheno_nm_ct = pheno_nm_ct;
  uint32_t case_ct = popcount_longs(pheno_c, unfiltered_sample_ctl);
  uint32_t ctrl_ct = pheno_nm_ct - case_ct;
  uint32_t male_ct = popcount_longs_intersect(sex_male, pheno_nm, unfiltered_sample_ctl);
  uint32_t case_ct_y = popcount_longs_intersect(sex_male, pheno_c, unfiltered_sample_ctl);
  uint32_t ctrl_ct_y = male_ct - case_ct_y;
  uint32_t cur_case_ct = case_ct;
  uint32_t cur_ctrl_ct = ctrl_ct;
  uint32_t chrom_end = 0;
  uint32_t mperm_dump_all = 0;
  uint32_t max_thread_ct = g_thread_ct;
  uintptr_t pheno_male_nm_ctl = BITCT_TO_WORDCT(male_ct);
  int32_t y_code = chrom_info_ptr->xymt_codes[Y_OFFSET];
  int32_t retval = 0;
  uint32_t uibuf[4];
  uintptr_t* loadbuf_raw;
  uintptr_t* pheno_nm2;
  uintptr_t* cur_pheno_nm2;
  uintptr_t* pheno_c_collapsed;
  uintptr_t* cur_pheno_c_collapsed;
  uintptr_t* missing_bitfield;
  uintptr_t* marker_exclude;
  uintptr_t* loadbuf_ptr;
  double* dptr;
  uint32_t* marker_idx_to_uidx;
  char* outname_end2;
  char* wptr;
  uintptr_t marker_uidx_end;
  uintptr_t marker_ct;
  uintptr_t marker_unstopped_ct;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t block_size;
  uintptr_t block_end;
  uintptr_t perm_idx;
  uintptr_t ulii;
  double pval;
  double cur_case_ct_recip;
  double cur_ctrl_ct_recip;
  double dxx;
  double dyy;
  double dzz;
  uint32_t missing_ct;
  uint32_t marker_cidx;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  if ((!case_ct) || (!ctrl_ct)) {
    logerrprint("Warning: Skipping --test-missing since at least one case and one control is\nrequired.\n");
    goto testmiss_ret_1;
  }
  cur_case_ct_recip = 1.0 / ((double)((int32_t)case_ct));
  cur_ctrl_ct_recip = 1.0 / ((double)((int32_t)ctrl_ct));
  // Y chromosome requires special handling--only male genotypes should be
  // considered.
  if ((y_code == -2) || (!is_set(chrom_info_ptr->chrom_mask, y_code))) {
    skip_y = 1;
  } else if ((!case_ct_y) || (!ctrl_ct_y)) {
    logerrprint("Warning: --test-missing is skipping Y chromosome since at least one male case\nand one male control are necessary.\n");
    skip_y = 1;
  }
  if (perm_maxt) {
    mperm_dump_all = mperm_save & MPERM_DUMP_ALL;
    perms_total = testmiss_mperm_val;
    if (bigstack_alloc_d(perms_total, &g_maxt_extreme_stat)) {
      goto testmiss_ret_NOMEM;
    }
    for (uii = 0; uii < perms_total; uii++) {
      g_maxt_extreme_stat[uii] = 1;
    }
    if (mperm_dump_all) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(outname, "w", &outfile_msa)) {
        goto testmiss_ret_OPEN_FAIL;
      }
      LOGPRINTFWW("Dumping all permutation p-values to %s .\n", outname);
    }
  } else {
    mperm_save = 0;
    if (perm_adapt) {
      g_aperm_alpha = apip->alpha;
      perms_total = apip->max;
    }
  }
  // Sites with no (or all) missing calls are now excluded from the permutation
  // test.  Since it's likely that many such sites exist, we postpone the
  // associated memory allocations until after the basic .missing report is
  // generated.
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &pheno_nm2) ||
      bigstack_alloc_ul(pheno_nm_ctl, &pheno_c_collapsed) ||
      bigstack_alloc_ul(pheno_nm_ctl, &missing_bitfield) ||
      bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude)) {
    goto testmiss_ret_NOMEM;
  }
  memcpy(marker_exclude, marker_exclude_orig, unfiltered_marker_ctl * sizeof(intptr_t));
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  init_quaterarr_from_bitarr(pheno_nm, unfiltered_sample_ct, pheno_nm2);
  cur_pheno_nm2 = pheno_nm2;
  copy_bitarr_subset(pheno_c, pheno_nm, unfiltered_sample_ct, pheno_nm_ct, pheno_c_collapsed);
  cur_pheno_c_collapsed = pheno_c_collapsed;
  if (!skip_y) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl2, &pheno_male_nm2) ||
        bigstack_alloc_ul(pheno_male_nm_ctl, &pheno_c_collapsed_male)) {
      goto testmiss_ret_NOMEM;
    }
    // temporary non-excluded male bitfield
    memcpy(pheno_male_nm2, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
    bitvec_and(sex_male, unfiltered_sample_ctl, pheno_male_nm2);
    copy_bitarr_subset(pheno_c, pheno_male_nm2, unfiltered_sample_ct, male_ct, pheno_c_collapsed_male);
    memcpy(pheno_male_nm2, pheno_nm2, unfiltered_sample_ctl2 * sizeof(intptr_t));
    apply_bitarr_mask_to_quaterarr_01(sex_male, unfiltered_sample_ct, pheno_male_nm2);
  }
  outname_end2 = memcpyb(outname_end, ".missing", 9);
  if (fopen_checked(outname, "w", &outfile)) {
    goto testmiss_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing --test-missing report to %s ... ", outname);
  fflush(stdout);
  sprintf(g_textbuf, " CHR %%%us     F_MISS_A     F_MISS_U            P \n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  if (ferror(outfile)) {
    goto testmiss_ret_WRITE_FAIL;
  }
  // technically this part could be even faster with some custom code, but not
  // worth the additional maintenance
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists, 1, pheno_nm, sex_male, &sample_hh_include2, &sample_hh_male_include2)) {
    goto testmiss_ret_NOMEM;
  }
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto testmiss_ret_READ_FAIL;
  }
  chrom_end = 0;
  // must be last allocation
  if (bigstack_alloc_d(marker_ct_orig, &g_orig_pvals)) {
    goto testmiss_ret_NOMEM;
  }
  dptr = g_orig_pvals;
  for (marker_idx = 0; marker_idx < marker_ct_orig; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude_orig, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude_orig, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
        goto testmiss_ret_NOMEM;
      }
    }
    if (marker_uidx >= chrom_end) {
      // exploit overflow
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &g_is_y, &uii, &is_haploid);
      if (!skip_y) {
        if (!g_is_y) {
          cur_pheno_nm2 = pheno_nm2;
          cur_pheno_nm_ct = pheno_nm_ct;
          cur_sample_ctl = pheno_nm_ctl;
          cur_case_ct = case_ct;
          cur_ctrl_ct = ctrl_ct;
	  cur_pheno_c_collapsed = pheno_c_collapsed;
	} else {
          cur_pheno_nm2 = pheno_male_nm2;
          cur_pheno_nm_ct = male_ct;
          cur_sample_ctl = pheno_male_nm_ctl;
          cur_case_ct = case_ct_y;
          cur_ctrl_ct = ctrl_ct_y;
	  cur_pheno_c_collapsed = pheno_c_collapsed_male;
	}
        cur_case_ct_recip = 1.0 / ((double)((int32_t)cur_case_ct));
        cur_ctrl_ct_recip = 1.0 / ((double)((int32_t)cur_ctrl_ct));
      } else if (g_is_y) {
        fill_bits(marker_uidx, chrom_end - marker_uidx, marker_exclude);
	marker_idx += chrom_end - marker_uidx - 1 - popcount_bit_idx(marker_exclude_orig, marker_uidx, chrom_end);
	marker_uidx = chrom_end - 1;
	continue;
      }
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, uii, g_textbuf));
      *wptr_start++ = ' ';
    }
    if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf_raw)) {
      goto testmiss_ret_READ_FAIL;
    }
    if (is_haploid && hh_exists) {
      haploid_fix(hh_exists, sample_hh_include2, sample_hh_male_include2, unfiltered_sample_ct, is_x, g_is_y, (unsigned char*)loadbuf_raw);
    }
    extract_collapsed_missing_bitfield(loadbuf_raw, unfiltered_sample_ct, cur_pheno_nm2, cur_pheno_nm_ct, missing_bitfield);
    missing_ct = popcount_longs(missing_bitfield, cur_sample_ctl);
    if ((!missing_ct) || (missing_ct == cur_pheno_nm_ct)) {
      SET_BIT(marker_uidx, marker_exclude);
      continue;
    }
    uii = popcount_longs_intersect(missing_bitfield, cur_pheno_c_collapsed, cur_sample_ctl);
    ujj = missing_ct - uii;
    pval = fisher22(uii, ujj, cur_case_ct - uii, cur_ctrl_ct - ujj, midp);
    *dptr++ = pval;
    if (!(pval <= pfilter)) {
      continue;
    }
    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
    *wptr++ = ' ';
    wptr = dtoa_g_wxp4x(((int32_t)uii) * cur_case_ct_recip, 12, ' ', wptr);
    wptr = dtoa_g_wxp4x(((int32_t)ujj) * cur_ctrl_ct_recip, 12, ' ', wptr);
    wptr = dtoa_g_wxp4x(MAXV(pval, output_min_p), 12, '\n', wptr);
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto testmiss_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto testmiss_ret_WRITE_FAIL;
  }
  logprint("done.\n");
  marker_ct = (uintptr_t)(dptr - g_orig_pvals);
  bigstack_shrink_top(g_orig_pvals, marker_ct * sizeof(double));
  if (mtest_adjust) {
    if (bigstack_alloc_ui(marker_ct, &marker_idx_to_uidx)) {
      goto testmiss_ret_NOMEM;
    }
    fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
    retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, nullptr, pfilter, output_min_p, mtest_adjust, 1, 0.0, nullptr, g_orig_pvals);
    if (retval) {
      goto testmiss_ret_1;
    }
  }
  if (do_perms) {
    if (!marker_ct) {
      logprint("Skipping --test-missing permutation test since all loci are degenerate.\n");
      goto testmiss_ret_1;
    }
    LOGPRINTF("Including %" PRIuPTR " loc%s in --test-missing permutation test.\n", marker_ct, (marker_ct == 1)? "us" : "i");
    if (mperm_dump_all) {
      g_textbuf[0] = '0';
      wptr = &(g_textbuf[1]);
      for (uii = 0; uii < marker_ct; uii++) {
        *wptr++ = ' ';
        dxx = g_orig_pvals[uii];
	if (dxx >= 0) {
	  wptr = dtoa_g(dxx, wptr);
	} else {
	  wptr = memcpya(wptr, "NA", 2);
	}
	if (wptr >= tbuf2) {
	  if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	    goto testmiss_ret_WRITE_FAIL;
	  }
	  wptr = g_textbuf;
	}
      }
      *wptr++ = '\n';
      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	goto testmiss_ret_WRITE_FAIL;
      }
    }

    if (!skip_y) {
      // maybe all Y chromosome markers had no missing calls?
      uii = get_chrom_start_vidx(chrom_info_ptr, (uint32_t)y_code);
      ujj = get_chrom_end_vidx(chrom_info_ptr, (uint32_t)y_code);
      if (popcount_bit_idx(marker_exclude, uii, ujj) == ujj - uii) {
	skip_y = 1;
      } else {
	if (bigstack_alloc_ul(pheno_nm_ctl, &sex_male_collapsed)) {
	  goto testmiss_ret_NOMEM;
	}
	copy_bitarr_subset(sex_male, pheno_nm, unfiltered_sample_ct, pheno_nm_ct, sex_male_collapsed);
      }
    }

    if (cluster_starts) {
      retval = cluster_include_and_reindex(unfiltered_sample_ct, pheno_nm, 1, pheno_c, pheno_nm_ct, 1, cluster_ct, cluster_map, cluster_starts, &g_perm_cluster_ct, &g_perm_cluster_map, &g_perm_cluster_starts, &g_perm_cluster_case_cts, &g_perm_cluster_cc_preimage);
      if (retval) {
	goto testmiss_ret_1;
      }
      if (!g_perm_cluster_ct) {
	logerrprint("Error: No size 2+ clusters for permutation test.\n");
	goto testmiss_ret_INVALID_CMDLINE;
      }
      retval = cluster_alloc_and_populate_magic_nums(g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts, &g_perm_tot_quotients, &g_perm_totq_magics, &g_perm_totq_preshifts, &g_perm_totq_postshifts, &g_perm_totq_incrs);
      if (retval) {
	goto testmiss_ret_1;
      }
    } else {
      g_perm_cluster_starts = nullptr;
    }
    if (max_thread_ct > perms_total) {
      max_thread_ct = perms_total;
    }
    if (bigstack_init_sfmtp(max_thread_ct)) {
      goto testmiss_ret_NOMEM;
    }
    if (bigstack_alloc_ul(MODEL_BLOCKSIZE * pheno_nm_ctv, &g_loadbuf) ||
	bigstack_calloc_ui(marker_ct, &g_perm_2success_ct) ||
	bigstack_alloc_ui(marker_ct, &g_missing_cts)) {
      goto testmiss_ret_NOMEM;
    }
    for (uii = 1; uii <= MODEL_BLOCKSIZE; uii++) {
      g_loadbuf[uii * pheno_nm_ctv - 2] = 0;
      g_loadbuf[uii * pheno_nm_ctv - 1] = 0;
    }
    uii = marker_ct;
    if (perm_maxt) {
      if (!mperm_dump_all) {
	if (bigstack_alloc_ui(6 * MODEL_BLOCKSIZE, &g_precomp_ui) ||
	    bigstack_alloc_d(2 * MODEL_BLOCKSIZE, &g_precomp_d)) {
	  goto testmiss_ret_NOMEM;
	}
      }
    } else {
      if (bigstack_alloc_ui(marker_ct, &g_perm_attempt_ct) ||
	  bigstack_calloc_uc(round_up_pow2(marker_ct, BYTECT), &g_perm_adapt_stop) ||
	  bigstack_alloc_ui(4 * MODEL_BLOCKSIZE, &g_precomp_ui)) {
	goto testmiss_ret_NOMEM;
      }
      for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	g_perm_attempt_ct[marker_idx] = perms_total;
      }
      g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)marker_ct)));
    }
    if (!cluster_starts) {
      g_perm_tot_quotient = 0x100000000LLU / pheno_nm_ct;
      magic_num(g_perm_tot_quotient, &g_perm_totq_magic, &g_perm_totq_preshift, &g_perm_totq_postshift, &g_perm_totq_incr);
    }
    marker_unstopped_ct = marker_ct;
    g_perm_is_1bit = 1;
    g_perms_done = 0;
    g_perm_pheno_nm_ct = pheno_nm_ct;
    g_perm_case_ct = case_ct;
    g_male_ct = male_ct;
    g_fisher_midp = midp;
    g_mperm_save_all = nullptr;
    // ----- begin main loop -----
  testmiss_more_perms:
    if (perm_adapt) {
      if (g_perms_done) {
	while (g_first_adapt_check <= g_perms_done) {
	  g_first_adapt_check += (int32_t)(apip->init_interval + ((int32_t)g_first_adapt_check) * apip->interval_slope);
	}
      } else {
	if (apip->min < apip->init_interval) {
	  g_first_adapt_check = (int32_t)(apip->init_interval);
	} else {
	  g_first_adapt_check = apip->min;
	}
	g_adaptive_intercept = apip->init_interval;
	g_adaptive_slope = apip->interval_slope;
      }
      g_perm_vec_ct = (bigstack_left() - CACHELINE + sizeof(int32_t)) / (pheno_nm_ctv * sizeof(intptr_t) + (1 - skip_y) * sizeof(int32_t));
    } else {
      // g_perm_vec_ct memory allocation dependencies:
      //   g_maxt_thread_results: (8 * g_perm_vec_ct, cacheline-aligned) *
      //     max_thread_ct
      //   g_perm_vecst: 16 * ((g_perm_vec_ct + 127) / 128) * pheno_nm_ct
      //   g_thread_git_wkspace: ((perm_vec_ct + 127) / 128) * 704 *
      //     max_thread_ct
      //   g_perm_vecs: pheno_nm_ctv * sizeof(intptr_t) * g_perm_vec_ct
      //   g_male_case_cts (if needed): sizeof(int32_t) * g_perm_vec_ct
      //   g_mperm_save_all (if needed): marker_ct * 8 * g_perm_vec_ct
      // Forcing g_perm_vec_ct to be a multiple of 128, total is
      //   g_perm_vec_ct * (13.5 * max_thread_ct + pheno_nm_ct / 8 + 4 +
      //                    sizeof(intptr_t) * pheno_nm_ctv
      //                    [+ marker_ct * sizeof(double) * mperm_save_all])
      if (mperm_dump_all) {
	g_perm_vec_ct = 128 * (bigstack_left() / (128 * sizeof(intptr_t) * pheno_nm_ctv + 1728LL * max_thread_ct + 16LL * pheno_nm_ct + 512 * (1 - skip_y) + 128LL * sizeof(double) * marker_ct));
      } else {
	g_perm_vec_ct = 128 * (bigstack_left() / (128 * sizeof(intptr_t) * pheno_nm_ctv + 1728LL * max_thread_ct + 16LL * pheno_nm_ct + 512 * (1 - skip_y)));
      }
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    } else if (!g_perm_vec_ct) {
      goto testmiss_ret_NOMEM;
    }
    bigstack_alloc_ul(g_perm_vec_ct * pheno_nm_ctv, &g_perm_vecs);
    g_perm_generation_thread_ct = MINV(max_thread_ct, g_perm_vec_ct);
    ulii = 0;
    if (!cluster_starts) {
      if (spawn_threads(threads, &generate_cc_perms_thread, g_perm_generation_thread_ct)) {
	goto testmiss_ret_THREAD_CREATE_FAIL;
      }
      generate_cc_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &generate_cc_cluster_perms_thread, g_perm_generation_thread_ct)) {
	goto testmiss_ret_THREAD_CREATE_FAIL;
      }
      generate_cc_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_perm_generation_thread_ct);
    g_assoc_thread_ct = max_thread_ct;
    if (perm_maxt) {
      bigstack_alloc_d(max_thread_ct * round_up_pow2(g_perm_vec_ct, CACHELINE_DBL), &g_maxt_thread_results);
#ifdef __LP64__
      ulii = ((g_perm_vec_ct + 127) / 128) * 4;
      bigstack_alloc_ui(ulii * pheno_nm_ct, &g_perm_vecst);
#else
      ulii = (g_perm_vec_ct + 31) / 32;
      bigstack_alloc_ui(ulii * pheno_nm_ct, &g_perm_vecst);
      ulii = ((g_perm_vec_ct + 127) / 128) * 4; // force 64-byte align
#endif
      bigstack_calloc_ui(ulii * 44 * max_thread_ct, &g_thread_git_wkspace);
      transpose_perm1s(g_perm_vecs, g_perm_vec_ct, pheno_nm_ct, g_perm_vecst);
      if (mperm_dump_all) {
	bigstack_alloc_d(marker_ct * g_perm_vec_ct, &g_mperm_save_all);
      }
    }
    if (!skip_y) {
      bigstack_alloc_ui(g_perm_vec_ct, &g_male_case_cts);
      for (perm_idx = 0; perm_idx < g_perm_vec_ct; perm_idx++) {
	g_male_case_cts[perm_idx] = popcount_longs_intersect(sex_male_collapsed, &(g_perm_vecs[perm_idx * pheno_nm_ctv]), pheno_nm_ctl);
      }
    }
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
      goto testmiss_ret_READ_FAIL;
    }
    marker_idx = 0;
    marker_idx2 = 0;
    chrom_end = 0;
    // only forced to terminate block at Y chromosome boundaries
    if (!skip_y) {
      marker_uidx_end = get_chrom_start_vidx(chrom_info_ptr, (uint32_t)y_code);
      pheno_male_nm_ctl = round_up_pow2(pheno_male_nm_ctl, 2);
    } else {
      marker_uidx_end = unfiltered_marker_ct;
    }
    do {
      block_size = 0;
      block_end = marker_unstopped_ct - marker_idx;
      if ((!marker_idx) && perm_maxt) {
	if (block_end > MODEL_BLOCKKEEP) {
	  block_end = MODEL_BLOCKKEEP;
	}
      } else if (block_end > MODEL_BLOCKSIZE) {
	block_end = MODEL_BLOCKSIZE;
      }
      if (marker_uidx >= marker_uidx_end) {
	marker_uidx_end = get_chrom_end_vidx(chrom_info_ptr, (uint32_t)y_code);
	if (marker_uidx >= marker_uidx_end) {
	  marker_uidx_end = unfiltered_marker_ct;
	}
      }
      do {
	if (perm_adapt && g_perm_adapt_stop[marker_idx2]) {
	  do {
	    marker_uidx++;
	    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	    marker_idx2++;
	  } while ((marker_uidx < marker_uidx_end) && g_perm_adapt_stop[marker_idx2]);
	  if (marker_uidx >= marker_uidx_end) {
	    break;
	  }
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto testmiss_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  // exploit overflow
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &uii, &is_haploid);
	  if (!g_is_y) {
	    g_perm_case_ct = case_ct;
	  } else {
	    g_perm_case_ct = case_ct_y;
	  }
	}
	if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf_raw)) {
	  goto testmiss_ret_READ_FAIL;
	}
	if (is_haploid && hh_exists) {
	  haploid_fix(hh_exists, sample_hh_include2, sample_hh_male_include2, unfiltered_sample_ct, is_x, g_is_y, (unsigned char*)loadbuf_raw);
	}
	loadbuf_ptr = &(g_loadbuf[block_size * pheno_nm_ctv]);
	extract_collapsed_missing_bitfield(loadbuf_raw, unfiltered_sample_ct, pheno_nm2, pheno_nm_ct, loadbuf_ptr);
	if (g_is_y) {
	  bitvec_and(sex_male_collapsed, pheno_nm_ctl, loadbuf_ptr);
	}
	if (!g_perms_done) {
	  missing_ct = popcount_longs(loadbuf_ptr, pheno_nm_ctl);
	  if (perm_adapt) {
	    g_missing_cts[marker_idx2] = missing_ct;
	  } else {
	    g_missing_cts[marker_idx + block_size] = missing_ct;
	  }
	}
	if (perm_adapt) {
	  g_adapt_m_table[block_size] = marker_idx2++;
	}
	block_size++;
	if (marker_idx + block_size == marker_unstopped_ct) {
	  break;
	}
	marker_uidx++;
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto testmiss_ret_READ_FAIL;
	  }
	}
      } while ((block_size < block_end) && (marker_uidx < marker_uidx_end));
      g_block_diff = block_size;
      if ((!mperm_dump_all) && ((!g_is_y) || (male_ct == pheno_nm_ct))) {
	if (perm_maxt) {
	  maxt_cur_extreme_stat = g_maxt_extreme_stat[0];
	  for (uii = 1; uii < g_perm_vec_ct; uii++) {
	    dxx = g_maxt_extreme_stat[uii];
	    if (dxx > maxt_cur_extreme_stat) {
	      maxt_cur_extreme_stat = dxx;
	    }
	  }
	}
	// need raw p-values for --mperm-save-all
	// valid case/control counts differ between permutations on Y
	// chromosome, and I won't bother with g_precomp_width just for that
	for (uii = 0; uii < block_size; uii++) {
	  if (perm_adapt) {
	    marker_cidx = g_adapt_m_table[uii];
	  } else {
	    marker_cidx = marker_idx + uii;
	  }
	  pval = g_orig_pvals[marker_cidx];
	  missing_ct = g_missing_cts[marker_cidx];
	  if (perm_adapt) {
	    fisher22_precomp_pval_bounds(pval, midp, case_ct, missing_ct, pheno_nm_ct, &(g_precomp_ui[uii * 4]), nullptr);
	  } else {
	    fisher22_precomp_pval_bounds(pval, midp, case_ct, missing_ct, pheno_nm_ct, &(g_precomp_ui[uii * 6]), nullptr);
	    fisher22_precomp_pval_bounds(maxt_cur_extreme_stat, midp, case_ct, missing_ct, pheno_nm_ct, uibuf, &(g_precomp_d[uii * 2]));
	    g_precomp_ui[uii * 6 + 4] = uibuf[2];
	    g_precomp_ui[uii * 6 + 5] = uibuf[3] - uibuf[2];
	  }
	}
      }
      ulii = 0;
      is_last_block = (marker_idx + block_size >= marker_unstopped_ct);
      if (perm_adapt) {
	if (spawn_threads2(threads, &testmiss_adapt_thread, max_thread_ct, is_last_block)) {
	  goto testmiss_ret_THREAD_CREATE_FAIL;
	}
	testmiss_adapt_thread((void*)ulii);
	join_threads2(threads, max_thread_ct, is_last_block);
      } else {
	g_maxt_block_base = marker_idx;
	if (spawn_threads2(threads, &testmiss_maxt_thread, max_thread_ct, is_last_block)) {
	  goto testmiss_ret_THREAD_CREATE_FAIL;
	}
	testmiss_maxt_thread((void*)ulii);
	join_threads2(threads, max_thread_ct, is_last_block);
	ulii = round_up_pow2(g_perm_vec_ct, CACHELINE_DBL);
	umm = block_size;
	if (umm > max_thread_ct) {
	  umm = max_thread_ct;
	}
	for (uii = 0; uii < max_thread_ct; uii++) {
	  dptr = &(g_maxt_thread_results[uii * ulii]);
	  ujj = g_perms_done;
	  ukk = ujj + g_perm_vec_ct;
	  for (; ujj < ukk; ujj++) {
	    dxx = *dptr++;
	    if (dxx < g_maxt_extreme_stat[ujj]) {
	      g_maxt_extreme_stat[ujj] = dxx;
	    }
	  }
	}
      }
      marker_idx += block_size;
    } while (marker_idx < marker_unstopped_ct);
    if (mperm_dump_all) {
      if (g_perms_done) {
	putc_unlocked(' ', stdout);
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = g_perm_vec_ct;
      ujj = 1 + g_perms_done;
      wptr = g_textbuf;
      tbuf2 = &(g_textbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32toa(uii + ujj, wptr);
	dptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = dptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = dtoa_g(dxx, wptr);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= tbuf2) {
	    if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	      goto testmiss_ret_WRITE_FAIL;
	    }
	    wptr = g_textbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile_msa)) {
	goto testmiss_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    // really should postpone this for --assoc/--model too
    g_perms_done += g_perm_vec_ct;
    bigstack_reset(g_perm_vecs);
    if (g_perms_done < perms_total) {
      if (perm_adapt) {
	marker_unstopped_ct = marker_ct - popcount01_longs((uintptr_t*)g_perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
	if (!marker_unstopped_ct) {
	  goto testmiss_adapt_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done != 1)? "s" : "");
      fflush(stdout);
      goto testmiss_more_perms;
    }
    if (perm_adapt) {
    testmiss_adapt_perm_count:
      g_perms_done = 0;
      for (ulii = 0; ulii < marker_ct; ulii++) {
	if (g_perm_attempt_ct[ulii] > g_perms_done) {
	  g_perms_done = g_perm_attempt_ct[ulii];
	  if (g_perms_done == perms_total) {
	    break;
	  }
	}
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("%u %s permutation%s complete.\n", g_perms_done, perm_maxt? "max(T)" : "(adaptive)", (g_perms_done != 1)? "s" : "");
    if (perm_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	ulii = outname_end - outname;
	memcpy(outname_end, ".mperm.dump.best", 17);
	LOGPRINTFWW("Dumping best permutation p-values to %s .\n", outname);
	if (fopen_checked(outname, "w", &outfile)) {
	  goto testmiss_ret_OPEN_FAIL;
	}
	dxx = 1.0;
	for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	  if (g_orig_pvals[marker_idx] < dxx) {
	    dxx = g_orig_pvals[marker_idx];
	  }
	}
	memcpy(g_textbuf, "0 ", 2);
	wptr = dtoa_gx(dxx, '\n', &(g_textbuf[2]));
	if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile)) {
	  goto testmiss_ret_WRITE_FAIL;
	}
	for (uii = 0; uii < perms_total; uii++) {
	  wptr = uint32toa_x(uii + 1, ' ', g_textbuf);
	  wptr = dtoa_gx(g_maxt_extreme_stat[uii], '\n', wptr);
	  if (fwrite_checked(g_textbuf, (uintptr_t)(wptr - g_textbuf), outfile)) {
	    goto testmiss_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto testmiss_ret_WRITE_FAIL;
	}
	memcpy(outname_end, ".missing", 8);
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto testmiss_ret_OPEN_FAIL;
    }
    if (perm_adapt) {
      sprintf(g_textbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(g_textbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(g_maxt_extreme_stat, &(g_maxt_extreme_stat[perms_total]));
#else
      qsort(g_maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
    }
    /*
    if (perm_maxt) {
      printf("extreme stats: %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[perms_total - 1]);
    }
    */
    fprintf(outfile, g_textbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, uii, g_textbuf));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	if (perm_adapt) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) * dxx;
	}
	if (pval <= pfilter) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	  wptr = &(wptr_start[1 + plink_maxsnp]);
	  if (!perm_count) {
	    wptr = dtoa_g_wxp4x(pval, 12, ' ', wptr);
	  } else {
	    wptr = dtoa_g_wxp4x(((double)g_perm_2success_ct[marker_idx]) * 0.5, 12, ' ', wptr);
	  }
	  if (perm_adapt) {
	    wptr = memseta(wptr, 32, 2);
	    wptr = uint32toa_w10(g_perm_attempt_ct[marker_idx], wptr);
	  } else {
	    // minimum p-value
	    dzz = (int32_t)(doublearr_greater_than(g_maxt_extreme_stat, perms_total, g_orig_pvals[marker_idx] * (1.0 + EPSILON)) + 1);
	    if (!perm_count) {
	      wptr = dtoa_g_wxp4(dzz * dyy, 12, wptr);
	    } else {
	      wptr = dtoa_g_wxp4(dzz - 1, 12, wptr);
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto testmiss_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto testmiss_loop_end;
	}
	marker_uidx++;
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      }
    }
  testmiss_loop_end:
    if (fclose_null(&outfile)) {
      goto testmiss_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
  }

  while (0) {
  testmiss_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  testmiss_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  testmiss_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  testmiss_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  testmiss_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  testmiss_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 testmiss_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
}

int32_t cluster_assoc_init(const char* flag_name, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uintptr_t* cluster_bitfield, uintptr_t** pheno_nm_11_ptr, uintptr_t** pheno_nm_nonmale_11_ptr, uintptr_t** pheno_nm_male_11_ptr, uint32_t** sample_to_cluster_pheno_ptr, uint32_t** cluster_pheno_gtots_ptr, uint32_t** cur_cluster_pheno_gtots_ptr, uint32_t** cluster_geno_cts_ptr, uintptr_t** loadbuf_raw_ptr, uint32_t* cluster_ct2_ptr) {
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uint32_t cluster_ct2 = 0;
  uint32_t sample_ct = 0;
  uint32_t cluster_end = 0;
  uint32_t case_ct_total = 0;
  uint32_t is_mh2 = (flag_name[4] == '2'); // yeah, this is a hack
  uintptr_t* pheno_nm_nonmale_11 = nullptr;
  uintptr_t* pheno_nm_male_11 = nullptr;
  uintptr_t* pheno_nm_11;
  uint32_t* sample_to_cluster_pheno;
  uint32_t* cluster_pheno_gtots;
  uint32_t cluster_idx;
  uint32_t sample_uidx;
  uint32_t ctrl_ct;
  uint32_t case_ct;
  uint32_t ctrl_male_ct;
  uint32_t case_male_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  if (cluster_ct < 2) {
    LOGERRPRINTF("Error: %s requires at least two valid clusters.\n", flag_name);
    return RET_INVALID_CMDLINE;
  }
  // 1. Identify clusters with at least one case and one control, and create
  //    new cluster data structures describing only these.
  // 2. Main loop efficiently skips homozygous A2s via use of CTZLU, and skips
  //    samples not in a valid cluster via application of the pheno_nm_11
  //    bitmask.  sample_to_cluster_pheno[] maps sample_uidx to (valid) cluster
  //    index (high 31 bits) and case/control status (low bit).
  if (bigstack_calloc_ul(unfiltered_sample_ctl2, pheno_nm_11_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, pheno_nm_nonmale_11_ptr) ||
      bigstack_calloc_ul(unfiltered_sample_ctl2, pheno_nm_male_11_ptr) ||
      bigstack_alloc_ui(unfiltered_sample_ct, sample_to_cluster_pheno_ptr) ||
      bigstack_alloc_ui(cluster_ct * 4, cluster_pheno_gtots_ptr)) {
    return RET_NOMEM;
  }
  pheno_nm_11 = *pheno_nm_11_ptr;
  pheno_nm_nonmale_11 = *pheno_nm_nonmale_11_ptr;
  pheno_nm_male_11 = *pheno_nm_male_11_ptr;
  sample_to_cluster_pheno = *sample_to_cluster_pheno_ptr;
  cluster_pheno_gtots = *cluster_pheno_gtots_ptr;
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    uii = cluster_end;
    cluster_end = cluster_starts[cluster_idx + 1];
    for (; uii < cluster_end; uii++) {
      sample_uidx = cluster_map[uii];
      if (is_set(pheno_nm, sample_uidx)) {
	if (is_mh2) {
	  goto cluster_assoc_init_valid;
	}
	if (is_set(pheno_c, sample_uidx)) {
	  // we have a case, check for a control
	  while (++uii < cluster_end) {
	    sample_uidx = cluster_map[uii];
	    if (is_set(pheno_nm, sample_uidx) && (!is_set(pheno_c, sample_uidx))) {
	      goto cluster_assoc_init_valid;
	    }
	  }
	  continue;
	} else {
	  // we have a control, check for a case
	  while (++uii < cluster_end) {
	    sample_uidx = cluster_map[uii];
	    if (is_set(pheno_c, sample_uidx)) {
	      goto cluster_assoc_init_valid;
	    }
	  }
	  continue;
	}
      }
    }
    continue;
  cluster_assoc_init_valid:
    for (uii = cluster_starts[cluster_idx], ctrl_ct = 0, ctrl_male_ct = 0, case_ct = 0, case_male_ct = 0; uii < cluster_end; uii++) {
      sample_uidx = cluster_map[uii];
      if (is_set(pheno_nm, sample_uidx)) {
        pheno_nm_11[sample_uidx / BITCT2] |= (3 * ONELU) << (2 * (sample_uidx % BITCT2));
	ukk = is_set(sex_male, sample_uidx);
	if (ukk) {
	  pheno_nm_male_11[sample_uidx / BITCT2] |= (3 * ONELU) << (2 * (sample_uidx % BITCT2));
	}
	ujj = is_set(pheno_c, sample_uidx);
	sample_to_cluster_pheno[sample_uidx] = 2 * cluster_ct2 + ujj;
	if (ujj) {
	  case_ct++;
	  case_male_ct += ukk;
	} else {
	  ctrl_ct++;
	  ctrl_male_ct += ukk;
	}
      }
    }
    cluster_pheno_gtots[4 * cluster_ct2] = ctrl_ct;
    cluster_pheno_gtots[4 * cluster_ct2 + 1] = ctrl_male_ct;
    cluster_pheno_gtots[4 * cluster_ct2 + 2] = case_ct;
    cluster_pheno_gtots[4 * cluster_ct2 + 3] = case_male_ct;
    sample_ct += ctrl_ct + case_ct;
    case_ct_total += case_ct;
    if (cluster_bitfield) {
      SET_BIT(cluster_idx, cluster_bitfield);
    }
    cluster_ct2++;
  }
  bitvec_andnot_copy(pheno_nm_11, pheno_nm_male_11, unfiltered_sample_ctl2, pheno_nm_nonmale_11);
  bigstack_shrink_top(cluster_pheno_gtots, cluster_ct2 * 4 * sizeof(int32_t));
  if (bigstack_alloc_ui(cluster_ct2 * 2, cur_cluster_pheno_gtots_ptr) ||
      bigstack_alloc_ui(cluster_ct2 * 4, cluster_geno_cts_ptr) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, loadbuf_raw_ptr)) {
    return RET_NOMEM;
  }
  if (cluster_ct2 < 2) {
    LOGERRPRINTF("Error: %s requires at least two valid clusters.\n", flag_name);
    return RET_INVALID_CMDLINE;
  } else if (sample_ct >= 0x40000000) {
    // silly, but I'll document this
    LOGERRPRINTF("Error: %s does not support >= 2^30 samples.\n", flag_name);
    return RET_INVALID_CMDLINE;
  }
  LOGPRINTF("%s: %u valid clusters, with a total of %u cases and %u controls.\n", flag_name, cluster_ct2, case_ct_total, sample_ct - case_ct_total);
  (*loadbuf_raw_ptr)[unfiltered_sample_ctl2 - 1] = 0;
  *cluster_ct2_ptr = cluster_ct2;
  return 0;
}

int32_t cluster_assoc_load_one(FILE* bedfile, uintptr_t bed_offset, uintptr_t* marker_exclude, uintptr_t unfiltered_sample_ct, uintptr_t* sample_hh_include2, uintptr_t* sample_hh_male_include2, uintptr_t* loadbuf_raw, uintptr_t* pheno_nm_11, uintptr_t* pheno_nm_nonmale_11, uintptr_t* pheno_nm_male_11, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uint32_t hh_or_mt_exists, char* chrom_name_buf, uint32_t cluster_ct2, uint32_t* sample_to_cluster_pheno, uint32_t* cluster_pheno_gtots, uint32_t* cur_cluster_pheno_gtots, uint32_t* cluster_geno_cts, uintptr_t* marker_uidx_ptr, uint32_t* chrom_end_ptr, uint32_t* chrom_fo_idx_ptr, uint32_t* min_ploidy_1_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr, char** chrom_name_pp, uint32_t* chrom_name_len_ptr) {
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uint32_t min_ploidy_1 = *min_ploidy_1_ptr;
  uintptr_t cur_word;
  uintptr_t* ulptr;
  uintptr_t* ulptr2;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t chrom_idx;
  uint32_t cpidx;
  uint32_t sample_uidx_base;
  uint32_t sample_uidx;
  uint32_t uii;
  uint32_t ujj;
  if (IS_SET(marker_exclude, marker_uidx)) {
    marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
    *marker_uidx_ptr = marker_uidx;
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
      return RET_READ_FAIL;
    }
  }
  if (marker_uidx >= (*chrom_end_ptr)) {
    chrom_fo_idx = *chrom_fo_idx_ptr;
    do {
      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
    } while (marker_uidx >= chrom_end);
    *chrom_end_ptr = chrom_end;
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    min_ploidy_1 = is_set(chrom_info_ptr->haploid_mask, chrom_idx) || (chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[MT_OFFSET]);
    *chrom_fo_idx_ptr = chrom_fo_idx;
    *min_ploidy_1_ptr = min_ploidy_1;
    *is_x_ptr = (chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET]);
    *is_y_ptr = (chrom_idx == (uint32_t)chrom_info_ptr->xymt_codes[Y_OFFSET]);
    if (!min_ploidy_1) {
      for (cpidx = 0; cpidx < 2 * cluster_ct2; cpidx++) {
	cur_cluster_pheno_gtots[cpidx] = 2 * cluster_pheno_gtots[cpidx * 2];
      }
    } else if (*is_x_ptr) {
      for (cpidx = 0; cpidx < 2 * cluster_ct2; cpidx++) {
	cur_cluster_pheno_gtots[cpidx] = 2 * cluster_pheno_gtots[cpidx * 2] - cluster_pheno_gtots[cpidx * 2 + 1];
      }
    } else if (*is_y_ptr) {
      for (cpidx = 0; cpidx < 2 * cluster_ct2; cpidx++) {
	cur_cluster_pheno_gtots[cpidx] = cluster_pheno_gtots[cpidx * 2 + 1];
      }
    } else {
      for (cpidx = 0; cpidx < 2 * cluster_ct2; cpidx++) {
	cur_cluster_pheno_gtots[cpidx] = cluster_pheno_gtots[cpidx * 2];
      }
    }
    if (chrom_name_len_ptr) {
      *chrom_name_pp = chrom_name_buf5w4write(chrom_info_ptr, chrom_idx, chrom_name_len_ptr, chrom_name_buf);
    } else {
      // --mh2
      // chrom_name_buf = g_textbuf in this case, and we return wptr_start
      *chrom_name_pp = chrom_name_write(chrom_info_ptr, chrom_idx, chrom_name_buf);
      *(*chrom_name_pp)++ = '\t';
    }
  }
  if (load_raw(unfiltered_sample_ct4, bedfile, loadbuf_raw)) {
    return RET_READ_FAIL;
  }
  if (IS_SET(marker_reverse, marker_uidx)) {
    reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf_raw);
  }
  if (min_ploidy_1 && hh_or_mt_exists) {
    haploid_fix(hh_or_mt_exists, sample_hh_include2, sample_hh_male_include2, unfiltered_sample_ct, *is_x_ptr, *is_y_ptr, (unsigned char*)loadbuf_raw);
  }
  fill_uint_zero(4 * cluster_ct2, cluster_geno_cts);
  ulptr = loadbuf_raw;
  ulptr2 = pheno_nm_11;
  if ((!min_ploidy_1) || (*is_x_ptr)) {
    if (*is_x_ptr) {
      ulptr2 = pheno_nm_nonmale_11;
    }
    for (sample_uidx_base = 0; sample_uidx_base < unfiltered_sample_ct; sample_uidx_base += BITCT2) {
      cur_word = (~(*ulptr++)) & (*ulptr2++);
      while (cur_word) {
	uii = CTZLU(cur_word) & (BITCT - 2);
	ujj = (cur_word >> uii) & 3;
	sample_uidx = sample_uidx_base + (uii / 2);
	cpidx = sample_to_cluster_pheno[sample_uidx];
	// this does the following branchlessly:
	// 1. increment A1 count by one for heterozygous calls (ujj == 1)
	// 2. increment missing count by two when ujj == 2
	// 3. increment A1 count by two when ujj == 3
	cluster_geno_cts[cpidx * 2 + (ujj == 2)] += 2 - (ujj == 1);
	cur_word &= ~((3 * ONELU) << uii);
      }
    }
  }
  if (min_ploidy_1) {
    ulptr = loadbuf_raw;
    if ((*is_x_ptr) || (*is_y_ptr)) {
      ulptr2 = pheno_nm_male_11;
    }
    for (sample_uidx_base = 0; sample_uidx_base < unfiltered_sample_ct; sample_uidx_base += BITCT2) {
      cur_word = (~(*ulptr++)) & (*ulptr2++);
      while (cur_word) {
	uii = CTZLU(cur_word) & (BITCT - 2);
	ujj = (cur_word >> uii) & 3;
	sample_uidx = sample_uidx_base + (uii / 2);
	cpidx = sample_to_cluster_pheno[sample_uidx];
	// increments A1 count by one, or missing count by one
	cluster_geno_cts[cpidx * 2 + 3 - ujj] += 1;
	cur_word &= ~((3 * ONELU) << uii);
      }
    }
  }
  return 0;
}

int32_t cmh_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t cmh_mperm_val, uint32_t cmh_modifier, double ci_size, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists, Set_info* sip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  FILE* outfile_msa = nullptr;
  uintptr_t* sample_hh_include2 = nullptr;
  uintptr_t* sample_hh_male_include2 = nullptr;
  uint32_t* orig_df = nullptr;
  char* chrom_name_ptr = nullptr;
  uint32_t breslow_day = cmh_modifier & CLUSTER_CMH_BD;
  uint32_t perm_bd = cmh_modifier & CLUSTER_CMH_PERM_BD;
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t chrom_end = 0;
  uint32_t chrom_name_len = 0;
  uint32_t pct = 0;
  uint32_t min_ploidy_1 = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  int32_t retval = 0;
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_SLEN];
  uintptr_t* pheno_nm_11;
  uintptr_t* pheno_nm_nonmale_11;
  uintptr_t* pheno_nm_male_11;
  uintptr_t* loadbuf_raw;
  double* orig_chisq;
  double* dptr;
  char* wptr;
  uint32_t* sample_to_cluster_pheno;
  uint32_t* cluster_pheno_gtots;
  uint32_t* cur_cluster_pheno_gtots;
  uint32_t* cluster_geno_cts;
  uint32_t* marker_idx_to_uidx;
  uint32_t* uiptr;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  double ci_zt;
  double allele_ct_recip;
  double allele_ctm1_recip;
  double ctrl_ctd;
  double case_ctd;
  double ctrl_a1_ctd;
  double ctrl_a2_ctd;
  double case_a1_ctd;
  double case_a2_ctd;
  double a1_ctd;
  double a2_ctd;
  double mean_case_a1d;
  double var_case_a1d;
  double cmh_stat;
  double cmh_denom;
  double r2;
  double s2;
  double rtot;
  double stot;
  double v1;
  double v2;
  double v3;
  double odds_ratio;
  double se;
  double log_or;
  double pval;
  double one_minus_odds_ratio;
  double double_1mor_recip;
  double bdx2;
  double amax;
  double bb;
  double discrim;
  double as_plus;
  double as_minus;
  double a_star;
  double b_star;
  double c_star;
  double d_star;
  double dxx;
  double dyy;
  uint32_t cluster_idx;
  uint32_t loop_end;
  uint32_t ctrl_ct;
  uint32_t case_ct;
  uint32_t cluster_ct2;
  uint32_t allele_ct;
  uint32_t uii;
  int32_t cur_df;

  // The best data structures for permutation testing are somewhat different
  // from those for the single-pass computation, so we separate the logic.

  retval = cluster_assoc_init("--mh/--bd", unfiltered_sample_ct, pheno_nm, pheno_c, sex_male, cluster_ct, cluster_map, cluster_starts, nullptr, &pheno_nm_11, &pheno_nm_nonmale_11, &pheno_nm_male_11, &sample_to_cluster_pheno, &cluster_pheno_gtots, &cur_cluster_pheno_gtots, &cluster_geno_cts, &loadbuf_raw, &cluster_ct2);
  if (retval) {
    goto cmh_assoc_ret_1;
  }
  if (breslow_day && (cluster_ct2 > 10) && (!perm_bd)) {
    logerrprint("Warning: Breslow-Day statistics are unreliable with a large number of small\nclusters.  You may want to look at empirical p-values from the 'perm-bd'\nadaptive permutation test.\n");
  }

  memcpy(outname_end, ".cmh", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto cmh_assoc_ret_OPEN_FAIL;
  }
  if (ci_size == 0.0) {
    ci_size = 0.95;
  }
  ci_zt = ltqnorm(1 - (1 - ci_size) / 2);
  LOGPRINTFWW5("Writing report to %s ... ", outname);
  fputs("0%", stdout);
  fflush(stdout);
  sprintf(g_textbuf, " CHR %%%us         BP   A1      MAF   A2      CHISQ          P         OR         SE        ", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  uii = (uint32_t)((int32_t)(ci_size * (100 + EPSILON)));
  if (uii >= 10) {
    fprintf(outfile, "L%u        U%u ", uii, uii);
  } else {
    fprintf(outfile, " L%u         U%u ", uii, uii);
  }
  if (breslow_day) {
    fputs("  CHISQ_BD       P_BD ", outfile);
  }
  if (putc_checked('\n', outfile)) {
    goto cmh_assoc_ret_WRITE_FAIL;
  }
  if ((chrom_info_ptr->xymt_codes[MT_OFFSET] != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[MT_OFFSET])) {
    hh_or_mt_exists |= NXMHH_EXISTS;
  }
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_or_mt_exists, 1, pheno_nm, sex_male, &sample_hh_include2, &sample_hh_male_include2)) {
    goto cmh_assoc_ret_NOMEM;
  }
  if (bigstack_alloc_d(marker_ct, &orig_chisq)) {
    goto cmh_assoc_ret_NOMEM;
  }
  if (perm_bd) {
    if (bigstack_alloc_ui(marker_ct, &orig_df)) {
      goto cmh_assoc_ret_NOMEM;
    }
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto cmh_assoc_ret_READ_FAIL;
  }
  dptr = orig_chisq;
  loop_end = marker_ct / 100;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (cluster_assoc_load_one(bedfile, bed_offset, marker_exclude, unfiltered_sample_ct, sample_hh_include2, sample_hh_male_include2, loadbuf_raw, pheno_nm_11, pheno_nm_nonmale_11, pheno_nm_male_11, marker_reverse, chrom_info_ptr, hh_or_mt_exists, chrom_name_buf, cluster_ct2, sample_to_cluster_pheno, cluster_pheno_gtots, cur_cluster_pheno_gtots, cluster_geno_cts, &marker_uidx, &chrom_end, &chrom_fo_idx, &min_ploidy_1, &is_x, &is_y, &chrom_name_ptr, &chrom_name_len)) {
      goto cmh_assoc_ret_READ_FAIL;
    }
    cmh_stat = 0.0;
    cmh_denom = 0.0;
    rtot = 0.0;
    stot = 0.0;
    v1 = 0.0;
    v2 = 0.0;
    v3 = 0.0;
    for (cluster_idx = 0, uiptr = cluster_geno_cts; cluster_idx < cluster_ct2; cluster_idx++, uiptr = &(uiptr[4])) {
      ctrl_ct = cur_cluster_pheno_gtots[cluster_idx * 2] - uiptr[1];
      case_ct = cur_cluster_pheno_gtots[cluster_idx * 2 + 1] - uiptr[3];
      // skip cluster if all controls missing, or all cases missing
      if (ctrl_ct && case_ct) {
	allele_ct = ctrl_ct + case_ct;
	allele_ct_recip = 1.0 / ((double)((int32_t)allele_ct));
	allele_ctm1_recip = 1.0 / ((double)((int32_t)(allele_ct - 1)));
	ctrl_ctd = (double)((int32_t)ctrl_ct);
	case_ctd = (double)((int32_t)case_ct);
	ctrl_a1_ctd = (double)((int32_t)uiptr[0]);
	ctrl_a2_ctd = ctrl_ctd - ctrl_a1_ctd;
	case_a1_ctd = (double)((int32_t)uiptr[2]);
	case_a2_ctd = case_ctd - case_a1_ctd;
	a1_ctd = ctrl_a1_ctd + case_a1_ctd;
	a2_ctd = ctrl_a2_ctd + case_a2_ctd;
        mean_case_a1d = case_ctd * a1_ctd * allele_ct_recip;
	var_case_a1d = ctrl_ctd * case_ctd * a1_ctd * a2_ctd * allele_ct_recip * allele_ct_recip * allele_ctm1_recip;
	cmh_stat += case_a1_ctd - mean_case_a1d;
        cmh_denom += var_case_a1d;
	r2 = case_a1_ctd * ctrl_a2_ctd * allele_ct_recip;
	s2 = case_a2_ctd * ctrl_a1_ctd * allele_ct_recip;
        rtot += r2;
        stot += s2;
	v1 += allele_ct_recip * r2 * (case_a1_ctd + ctrl_a2_ctd);
	v2 += allele_ct_recip * s2 * (case_a2_ctd + ctrl_a1_ctd);
        v3 += allele_ct_recip * ((case_a1_ctd + ctrl_a2_ctd) * s2 + (case_a2_ctd + ctrl_a1_ctd) * r2);
      }
    }
    cmh_stat *= cmh_stat / cmh_denom;
    odds_ratio = rtot / stot;
    se = sqrt(v1 / (2 * rtot * rtot) + v2 / (2 * stot * stot) + v3 / (2 * rtot * stot));
    log_or = log(odds_ratio);
    pval = chiprob_p(cmh_stat, 1);
    if (cmh_stat >= 0.0) {
      *dptr++ = cmh_stat;
    } else {
      *dptr++ = -9;
    }
    if ((pfilter == 2.0) || ((pval <= pfilter) && (pval != -9))) {
      wptr = memcpyax(g_textbuf, chrom_name_ptr, chrom_name_len, ' ');
      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr);
      *wptr++ = ' ';
      wptr = uint32toa_w10x(marker_pos[marker_uidx], ' ', wptr);
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto cmh_assoc_ret_WRITE_FAIL;
      }
      fputs_w4(marker_allele_ptrs[marker_uidx * 2], outfile);
      g_textbuf[0] = ' ';
      wptr = dtoa_g_wxp4x(1.0 - set_allele_freqs[marker_uidx], 8, ' ', &(g_textbuf[1]));
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto cmh_assoc_ret_WRITE_FAIL;
      }
      fputs_w4(marker_allele_ptrs[marker_uidx * 2 + 1], outfile);
      if (realnum(cmh_stat)) {
	g_textbuf[0] = ' ';
	wptr = dtoa_g_wxp4x(cmh_stat, 10, ' ', &(g_textbuf[1]));
	wptr = dtoa_g_wxp4x(MAXV(pval, output_min_p), 10, ' ', wptr);
      } else {
        wptr = memcpya(g_textbuf, "         NA         NA ", 23);
      }
      if (realnum(odds_ratio)) {
        wptr = dtoa_g_wxp4x(odds_ratio, 10, ' ', wptr);
      } else {
	wptr = memcpya(wptr, "        NA ", 11);
      }
      if (realnum(se)) {
        wptr = dtoa_g_wxp4x(se, 10, ' ', wptr);
	dxx = ci_zt * se;
	dyy = exp(log_or - dxx);
	if (realnum(dyy)) {
          wptr = dtoa_g_wxp4x(dyy, 10, ' ', wptr);
	} else {
	  wptr = memcpya(wptr, "        NA ", 11);
	}
	dyy = exp(log_or + dxx);
        if (realnum(dyy)) {
          wptr = dtoa_g_wxp4x(dyy, 10, ' ', wptr);
	} else {
	  wptr = memcpya(wptr, "        NA ", 11);
	}
      } else {
	wptr = memcpya(wptr, "        NA         NA         NA ", 33);
      }
      if (breslow_day) {
	if (realnum(odds_ratio) && (odds_ratio != 1.0)) {
	  one_minus_odds_ratio = 1.0 - odds_ratio;
          double_1mor_recip = 0.5 / one_minus_odds_ratio;
	  bdx2 = 0.0;
	  cur_df = -1;
	  for (cluster_idx = 0, uiptr = cluster_geno_cts; cluster_idx < cluster_ct2; cluster_idx++, uiptr = &(uiptr[4])) {
	    ctrl_ct = cur_cluster_pheno_gtots[cluster_idx * 2] - uiptr[1];
	    case_ct = cur_cluster_pheno_gtots[cluster_idx * 2 + 1] - uiptr[3];
	    if (ctrl_ct && case_ct) {
	      cur_df++;
	      ctrl_ctd = (double)((int32_t)ctrl_ct);
	      case_ctd = (double)((int32_t)case_ct);
	      ctrl_a1_ctd = (double)((int32_t)uiptr[0]);
	      case_a1_ctd = (double)((int32_t)uiptr[2]);
	      a1_ctd = ctrl_a1_ctd + case_a1_ctd;
	      amax = MINV(case_ctd, a1_ctd);
	      bb = ctrl_ctd + case_ctd * odds_ratio - a1_ctd * one_minus_odds_ratio;
	      discrim = sqrt(bb * bb + 4 * one_minus_odds_ratio * odds_ratio * case_ctd * a1_ctd);
	      as_plus = (-bb + discrim) * double_1mor_recip;
	      as_minus = (-bb - discrim) * double_1mor_recip;
	      a_star = ((as_minus <= amax) && (as_minus >= 0))? as_minus : as_plus;
              b_star = case_ctd - a_star;
              c_star = a1_ctd - a_star;
              d_star = ctrl_ctd - a1_ctd + a_star;

              // concordance fix (25 May 2018): print NA,NA instead of inf,0
              if ((a_star == 0.0) || (b_star == 0.0) || (c_star == 0.0) || (d_star == 0.0)) {
                goto cmh_assoc_bd_fail;
              }

	      // inverse variance
              dxx = 1.0 / a_star + 1.0 / b_star + 1.0 / c_star + 1.0 / d_star;

	      dyy = case_a1_ctd - a_star;
	      bdx2 += dyy * dyy * dxx;
	    }
	  }
	  pval = chiprob_p(bdx2, cur_df);
	  if (pval > -1) {
	    wptr = dtoa_g_wxp4x(bdx2, 10, ' ', wptr);
	    wptr = dtoa_g_wxp4x(MAXV(pval, output_min_p), 10, ' ', wptr);
	  } else {
	    goto cmh_assoc_bd_fail;
	  }
	} else {
	cmh_assoc_bd_fail:
	  wptr = memcpya(wptr, "        NA         NA ", 22);
	}
      }
      *wptr++ = '\n';
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto cmh_assoc_ret_WRITE_FAIL;
      }
    }
    if (marker_idx >= loop_end) {
      if (marker_idx < marker_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = (marker_idx * 100LLU) / marker_ct;
        printf("\b\b%u%%", pct);
        fflush(stdout);
        loop_end = (((uint64_t)pct + 1LLU) * marker_ct) / 100;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto cmh_assoc_ret_WRITE_FAIL;
  }
  if (pct >= 10) {
    putc_unlocked('\b', stdout);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  if (mtest_adjust) {
    if (bigstack_alloc_ui(marker_ct, &marker_idx_to_uidx)) {
      goto cmh_assoc_ret_NOMEM;
    }
    fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
    retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, orig_chisq, pfilter, output_min_p, mtest_adjust, 0, adjust_lambda, nullptr, nullptr);
  }

  if (cmh_modifier & (CLUSTER_CMH_PERM | CLUSTER_CMH_MPERM)) {
    logerrprint("Error: --mh/--bd permutation tests are currently under development.\n");
    goto cmh_assoc_ret_INVALID_CMDLINE;
  }

  // Given the genotypes at a marker, the following quantities are invariant
  // through permutations:
  // * set of possibly valid clusters (2+ nonmissing genotypes)
  // * allele counts in each cluster
  // while the following quantities need to be recomputed:
  // * [case x A1] and [case x A2] counts in each cluster (control counts can
  //   then be determined via subtraction; but note that [case x A2] CANNOT
  //   generally be determined from [case x A1] because the number of cases
  //   with missing genotypes may vary, though we could special-case
  //   no-missing-genotypes if this ever becomes popular enough to justify the
  //   complexity)
  //
  // To handle both large and small clusters efficiently without too much
  // special-casing, we preprocess the raw data so that each cluster's
  // genotypes occupy separate words.  (Exception: on 64-bit systems, clusters
  // of size <= 16 are stuffed into 4 bytes, to improve memory efficiency.)
  // This allows the inner loops to be based on bitwise operations and
  // sequential memory accessses.  We also scan for clusters containing only a
  // single genotype, or less than two nonmissing genotypes, and exclude them
  // from the main loop.

  // ...

  while (0) {
  cmh_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cmh_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cmh_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cmh_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  cmh_assoc_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 cmh_assoc_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
}

int32_t cmh2_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t* sample_hh_include2 = nullptr;
  uintptr_t* sample_hh_male_include2 = nullptr;
  char* wptr_start = nullptr;
  uint32_t chrom_fo_idx = 0xffffffffU;
  uint32_t chrom_end = 0;
  uint32_t pct = 0;
  uint32_t min_ploidy_1 = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t cluster_ct1 = 0;
  uint32_t ctrl_ct = 0;
  uint32_t case_ct = 0;
  int32_t retval = 0;
  uintptr_t* pheno_nm_11;
  uintptr_t* pheno_nm_nonmale_11;
  uintptr_t* pheno_nm_male_11;
  uintptr_t* loadbuf_raw;
  char* wptr;
  MATRIX_INVERT_BUF1_TYPE* mi_buf;
  double* ty_ctrl;
  double* ty_case;
  double* n0;
  double* u0;
  double* v0;
  double* dbl_2d_buf;
  double* dptr;
  double* dptr2;
  uint32_t* sample_to_cluster_pheno;
  uint32_t* cluster_pheno_gtots;
  uint32_t* cur_cluster_pheno_gtots;
  uint32_t* cluster_geno_cts;
  uint32_t* uiptr;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  double ctrl_a1_ctd; // Tx[.][]
  double case_a1_ctd;
  double ctrl_ctd; // T[.]
  double case_ctd;
  double cur_ty_ctrl;
  double cur_ty_case;
  double ctrl_umult;
  double case_umult;
  double ctrl_vmult; // (Tx[] * (T[] - Tx[])) / (T[] * T[] * (T[]-1))
  double case_vmult;
  double cur_ctrl_vmult;
  double cur_case_vmult;
  double chisq;
  double dxx;
  uint32_t cur_ctrl_ct;
  uint32_t cur_case_ct;
  uint32_t cluster_ctrl_ct;
  uint32_t cluster_case_ct;
  uint32_t ctrl_a1_ct;
  uint32_t case_a1_ct;
  uint32_t cur_cluster_ct;
  uint32_t cur_cluster_ctm1;
  uint32_t cluster_idx;
  uint32_t loop_end;
  uint32_t uii;
  // no reason to keep X/Y/MT/haploid restriction
  retval = cluster_assoc_init("--mh2", unfiltered_sample_ct, pheno_nm, pheno_c, sex_male, cluster_ct, cluster_map, cluster_starts, nullptr, &pheno_nm_11, &pheno_nm_nonmale_11, &pheno_nm_male_11, &sample_to_cluster_pheno, &cluster_pheno_gtots, &cur_cluster_pheno_gtots, &cluster_geno_cts, &loadbuf_raw, &cluster_ct1);
  for (cluster_idx = 0; cluster_idx < cluster_ct1; cluster_idx++) {
    ctrl_ct += cluster_pheno_gtots[4 * cluster_idx];
    case_ct += cluster_pheno_gtots[4 * cluster_idx + 2];
  }
  if ((ctrl_ct < 2) || (case_ct < 2)) {
    logerrprint("Error: --mh2 requires at least two cases and two controls.\n");
    goto cmh2_assoc_ret_INVALID_CMDLINE;
  }
#ifdef __LP64__
  if (cluster_ct1 > 46341) {
    // might actually be ok, but play it safe in case LAPACK matrix inversion
    // routine has an integer overflow here
    // (if/when we do permit this, will need to switch a few variables to type
    // uintptr_t)
    logerrprint("Error: --mh2 does not currently support more than 46341 clusters.\n");
    goto cmh2_assoc_ret_INVALID_CMDLINE;
  }
#endif
  if (bigstack_alloc_d(cluster_ct1, &ty_ctrl) ||
      bigstack_alloc_d(cluster_ct1, &ty_case) ||
      bigstack_alloc_d(cluster_ct1 - 1, &n0) ||
      bigstack_alloc_d(cluster_ct1 - 1, &u0) ||
      bigstack_alloc_d((cluster_ct1 - 1) * (cluster_ct1 - 1), &v0) ||
      bigstack_alloc_d((cluster_ct1 - 1) * (cluster_ct1 - 1), &dbl_2d_buf)) {
    goto cmh2_assoc_ret_NOMEM;
  }
  mi_buf = (MATRIX_INVERT_BUF1_TYPE*)bigstack_alloc((cluster_ct1 - 1) * MATRIX_INVERT_BUF1_ELEM_ALLOC);
  if (!mi_buf) {
    goto cmh2_assoc_ret_NOMEM;
  }
  if ((chrom_info_ptr->xymt_codes[MT_OFFSET] != -2) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->xymt_codes[MT_OFFSET])) {
    hh_or_mt_exists |= NXMHH_EXISTS;
  }
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_or_mt_exists, 1, pheno_nm, sex_male, &sample_hh_include2, &sample_hh_male_include2)) {
    goto cmh2_assoc_ret_NOMEM;
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto cmh2_assoc_ret_READ_FAIL;
  }
  memcpy(outname_end, ".cmh2", 6);
  if (fopen_checked(outname, "w", &outfile)) {
    goto cmh2_assoc_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing report to %s ... ", outname);
  fputs("0%", stdout);
  fflush(stdout);
  if (fputs_checked("CHR\tSNP\tCHISQ\tDF\tP\n", outfile)) {
    goto cmh2_assoc_ret_WRITE_FAIL;
  }
  loop_end = marker_ct / 100;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (cluster_assoc_load_one(bedfile, bed_offset, marker_exclude, unfiltered_sample_ct, sample_hh_include2, sample_hh_male_include2, loadbuf_raw, pheno_nm_11, pheno_nm_nonmale_11, pheno_nm_male_11, marker_reverse, chrom_info_ptr, hh_or_mt_exists, g_textbuf, cluster_ct1, sample_to_cluster_pheno, cluster_pheno_gtots, cur_cluster_pheno_gtots, cluster_geno_cts, &marker_uidx, &chrom_end, &chrom_fo_idx, &min_ploidy_1, &is_x, &is_y, &wptr_start, nullptr)) {
      goto cmh2_assoc_ret_READ_FAIL;
    }
    wptr = strcpyax(wptr_start, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
    cur_ctrl_ct = 0;
    cur_case_ct = 0;
    ctrl_a1_ct = 0;
    case_a1_ct = 0;
    cur_cluster_ct = 0;
    for (cluster_idx = 0, uiptr = cluster_geno_cts; cluster_idx < cluster_ct1; cluster_idx++, uiptr = &(uiptr[4])) {
      cluster_ctrl_ct = cur_cluster_pheno_gtots[cluster_idx * 2] - uiptr[1];
      cluster_case_ct = cur_cluster_pheno_gtots[cluster_idx * 2 + 1] - uiptr[3];
      uii = cluster_ctrl_ct + cluster_case_ct;
      if (uii) {
	// don't count toward cur_cluster_ct if all observations are missing
        n0[cur_cluster_ct] = (double)((int32_t)(uiptr[0] + uiptr[2]));
	ctrl_a1_ct += uiptr[0];
        case_a1_ct += uiptr[2];
	cur_ctrl_ct += cluster_ctrl_ct;
        cur_case_ct += cluster_case_ct;
	ty_ctrl[cur_cluster_ct] = (double)((int32_t)cluster_ctrl_ct);
	ty_case[cur_cluster_ct] = (double)((int32_t)cluster_case_ct);
	cur_cluster_ct++;
      }
    }

    // This is always a 2xJx2 test (where J = cluster ct), so we can omit PLINK
    // 1.07 calcMantelHaenszel_IxJxK code which only comes into play for larger
    // I/K values.
    if (((!cur_ctrl_ct) && cur_case_ct) || ((!cur_case_ct) && cur_ctrl_ct) || (cur_cluster_ct == 1)) {
      // may as well distinguish 0df from other problems
      wptr = memcpya(wptr, "0\t0\tNA\n", 7);
      goto cmh2_assoc_fail2;
    } else if ((cur_ctrl_ct < 2) || (cur_case_ct < 2) || (!cur_cluster_ct)) {
      goto cmh2_assoc_fail;
    }
    cur_cluster_ctm1 = cur_cluster_ct - 1;
    ctrl_ctd = (double)((int32_t)cur_ctrl_ct);
    case_ctd = (double)((int32_t)cur_case_ct);
    ctrl_a1_ctd = (double)((int32_t)ctrl_a1_ct);
    case_a1_ctd = (double)((int32_t)case_a1_ct);
    ctrl_umult = ctrl_a1_ctd / ctrl_ctd;
    case_umult = case_a1_ctd / case_ctd;
    ctrl_vmult = ctrl_umult * (ctrl_ctd - ctrl_a1_ctd) / (ctrl_ctd * (ctrl_ctd - 1));
    case_vmult = case_umult * (case_ctd - case_a1_ctd) / (case_ctd * (case_ctd - 1));
    for (cluster_idx = 0; cluster_idx < cur_cluster_ctm1; cluster_idx++) {
      // instead of a two-step process where e.g. U[][] is filled first, and
      // then columnwise sums are saved to U0, we just fill U0 directly.
      cur_ty_ctrl = ty_ctrl[cluster_idx];
      cur_ty_case = ty_case[cluster_idx];
      u0[cluster_idx] = cur_ty_ctrl * ctrl_umult + cur_ty_case * case_umult;
      cur_ctrl_vmult = -cur_ty_ctrl * ctrl_vmult;
      cur_case_vmult = -cur_ty_case * case_vmult;
      dptr = &(v0[cluster_idx * cur_cluster_ct]);
      // should be guaranteed to be nonnegative, no need for fabs()?
      *dptr++ = (cur_ty_ctrl - ctrl_ctd) * cur_ctrl_vmult + (cur_ty_case - case_ctd) * cur_case_vmult;
      for (uii = cluster_idx + 1; uii < cur_cluster_ctm1; uii++) {
	*dptr++ = ty_ctrl[uii] * cur_ctrl_vmult + ty_case[uii] * cur_case_vmult;
      }
    }
    for (cluster_idx = 0; cluster_idx < cur_cluster_ctm1; cluster_idx++) {
      dptr = &(v0[cluster_idx * cur_cluster_ctm1]);
      dptr2 = &(v0[cluster_idx]);
      for (uii = 0; uii < cluster_idx; uii++) {
	*dptr++ = dptr2[uii * cur_cluster_ctm1];
      }
    }

    if (!invert_matrix(cur_cluster_ctm1, v0, mi_buf, dbl_2d_buf)) {
      // Q = G'V{-1}G
      chisq = 0.0;
      for (cluster_idx = 0; cluster_idx < cur_cluster_ctm1; cluster_idx++) {
	dbl_2d_buf[cluster_idx] = n0[cluster_idx] - u0[cluster_idx];
      }
      dptr = v0;
      for (cluster_idx = 0; cluster_idx < cur_cluster_ctm1; cluster_idx++) {
	dxx = 0.0;
	dptr2 = dbl_2d_buf;
	for (uii = 0; uii < cur_cluster_ctm1; uii++) {
	  dxx += (*dptr++) * (*dptr2++);
	}
	chisq += dxx * (dbl_2d_buf[cluster_idx]);
      }
      wptr = dtoa_gx(chisq, '\t', wptr);
      wptr = uint32toa_x(cur_cluster_ctm1, '\t', wptr);
      dxx = chiprob_p(chisq, (int32_t)cur_cluster_ctm1);
      wptr = dtoa_gx(MAXV(dxx, output_min_p), '\n', wptr);
    } else {
    cmh2_assoc_fail:
      wptr = memcpya(wptr, "NA\tNA\tNA\n", 9);
    }
  cmh2_assoc_fail2:
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto cmh2_assoc_ret_WRITE_FAIL;
    }
    if (marker_idx >= loop_end) {
      if (marker_idx < marker_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
        pct = (marker_idx * 100LLU) / marker_ct;
        printf("\b\b%u%%", pct);
        fflush(stdout);
        loop_end = (((uint64_t)pct + 1LLU) * marker_ct) / 100;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto cmh2_assoc_ret_WRITE_FAIL;
  }
  if (pct >= 10) {
    putc_unlocked('\b', stdout);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  while (0) {
  cmh2_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  cmh2_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  cmh2_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  cmh2_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  cmh2_assoc_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t homog_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* outfile = nullptr;
  uintptr_t* sample_hh_include2 = nullptr;
  uintptr_t* sample_hh_male_include2 = nullptr;
  char* writebuf = g_textbuf;
  char* chrom_name_ptr = nullptr;
  uint32_t cluster_ct2 = 0;
  uint32_t chrom_fo_idx = 0xffffffffU;
  uint32_t chrom_end = 0;
  uint32_t chrom_name_len = 0;
  uint32_t pct = 0;
  uint32_t min_ploidy_1 = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  int32_t retval = 0;
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_SLEN];
  uintptr_t* cluster_bitfield;
  uintptr_t* pheno_nm_11;
  uintptr_t* pheno_nm_nonmale_11;
  uintptr_t* pheno_nm_male_11;
  uintptr_t* loadbuf_raw;
  double* cluster_tables;
  double* cluster_chisq;
  double* cluster_or;
  double* dptr;
  char* cluster_ids_collapsed;
  char* wptr_start;
  char* wptr;
  uint32_t* sample_to_cluster_pheno;
  uint32_t* cluster_pheno_gtots;
  uint32_t* cur_cluster_pheno_gtots;
  uint32_t* cluster_geno_cts;
  uint32_t* uiptr;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t ulii;
  double cluster_ct2d;
  double cluster_ct2m1d;
  double case_ctd;
  double ctrl_ctd;
  double case_a1_ctd;
  double case_a2_ctd;
  double ctrl_a1_ctd;
  double ctrl_a2_ctd;
  double case_a2_recip;
  double ctrl_a1_recip;
  double ln_or;
  double se_sq_recip;
  double x_total;
  double x_assoc1;
  double x_assoc2;
  double x_assoc;
  double dxx;
  uint32_t cluster_idx;
  uint32_t loop_end;
  ulii = 2 * max_marker_allele_len + MAX_ID_SLEN + max_marker_id_len + max_cluster_id_len + 256;
  if (ulii > MAXLINELEN) {
    if (bigstack_alloc_c(ulii, &writebuf)) {
      goto homog_assoc_ret_NOMEM;
    }
  }
  if (bigstack_end_calloc_ul(BITCT_TO_WORDCT(cluster_ct), &cluster_bitfield)) {
    goto homog_assoc_ret_NOMEM;
  }
  // Factor out common initialization with cmh_assoc().
  retval = cluster_assoc_init("--homog", unfiltered_sample_ct, pheno_nm, pheno_c, sex_male, cluster_ct, cluster_map, cluster_starts, cluster_bitfield, &pheno_nm_11, &pheno_nm_nonmale_11, &pheno_nm_male_11, &sample_to_cluster_pheno, &cluster_pheno_gtots, &cur_cluster_pheno_gtots, &cluster_geno_cts, &loadbuf_raw, &cluster_ct2);
  if (retval) {
    goto homog_assoc_ret_1;
  }
  if (cluster_ct == cluster_ct2) {
    cluster_ids_collapsed = cluster_ids;
  } else {
    if (bigstack_alloc_c(cluster_ct2 * max_cluster_id_len, &cluster_ids_collapsed)) {
      goto homog_assoc_ret_NOMEM;
    }
    for (ulii = 0, cluster_idx = 0; cluster_idx < cluster_ct2; ulii++, cluster_idx++) {
      next_set_ul_unsafe_ck(cluster_bitfield, &ulii);
      memcpy(&(cluster_ids_collapsed[cluster_idx * max_cluster_id_len]), &(cluster_ids[ulii * max_cluster_id_len]), max_cluster_id_len);
    }
  }
  bigstack_end_reset(bigstack_end_mark);
  cluster_ct2d = (double)((int32_t)cluster_ct2);
  cluster_ct2m1d = (double)((int32_t)cluster_ct2 - 1);
  if (bigstack_alloc_d(cluster_ct2 * 4, &cluster_tables) ||
      bigstack_alloc_d(cluster_ct2, &cluster_or) ||
      bigstack_alloc_d(cluster_ct2, &cluster_chisq)) {
    goto homog_assoc_ret_NOMEM;
  }
  if (cluster_ct2 > 10) {
    logerrprint("Warning: --homog statistics can be unreliable with small clusters.\n");
  }

  memcpy(outname_end, ".homog", 7);
  if (fopen_checked(outname, "w", &outfile)) {
    goto homog_assoc_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing report to %s ... ", outname);
  fputs("0%", stdout);
  fflush(stdout);
  // misaligned for backward compatibility
  sprintf(g_textbuf, " CHR %%%us   A1   A2      F_A      F_U      N_A      N_U     TEST      CHISQ   DF          P         OR\n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  if (chrom_info_ptr->xymt_codes[MT_OFFSET] != -2) {
    hh_or_mt_exists |= NXMHH_EXISTS;
  }
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_or_mt_exists, 1, pheno_nm, sex_male, &sample_hh_include2, &sample_hh_male_include2)) {
    goto homog_assoc_ret_NOMEM;
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto homog_assoc_ret_READ_FAIL;
  }
  loop_end = marker_ct / 100;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (cluster_assoc_load_one(bedfile, bed_offset, marker_exclude, unfiltered_sample_ct, sample_hh_include2, sample_hh_male_include2, loadbuf_raw, pheno_nm_11, pheno_nm_nonmale_11, pheno_nm_male_11, marker_reverse, chrom_info_ptr, hh_or_mt_exists, chrom_name_buf, cluster_ct2, sample_to_cluster_pheno, cluster_pheno_gtots, cur_cluster_pheno_gtots, cluster_geno_cts, &marker_uidx, &chrom_end, &chrom_fo_idx, &min_ploidy_1, &is_x, &is_y, &chrom_name_ptr, &chrom_name_len)) {
      goto homog_assoc_ret_READ_FAIL;
    }
    dptr = cluster_tables;
    x_total = 0.0;
    x_assoc1 = 0.0;
    x_assoc2 = 0.0;
    for (cluster_idx = 0, uiptr = cluster_geno_cts; cluster_idx < cluster_ct2; cluster_idx++, uiptr = &(uiptr[4])) {
      ctrl_ctd = (double)((int32_t)(1 + cur_cluster_pheno_gtots[cluster_idx * 2] - uiptr[1]));
      case_ctd = (double)((int32_t)(1 + cur_cluster_pheno_gtots[cluster_idx * 2 + 1] - uiptr[3]));
      ctrl_a1_ctd = (double)((int32_t)uiptr[0]) + 0.5;
      ctrl_a2_ctd = ctrl_ctd - ctrl_a1_ctd;
      case_a1_ctd = (double)((int32_t)uiptr[2]) + 0.5;
      case_a2_ctd = case_ctd - case_a1_ctd;
      *dptr++ = case_a1_ctd;
      *dptr++ = case_a2_ctd;
      *dptr++ = ctrl_a1_ctd;
      *dptr++ = ctrl_a2_ctd;
      case_a2_recip = 1.0 / case_a2_ctd;
      ctrl_a1_recip = 1.0 / ctrl_a1_ctd;
      dxx = case_a1_ctd * ctrl_a2_ctd * case_a2_recip * ctrl_a1_recip;
      cluster_or[cluster_idx] = dxx;
      ln_or = log(dxx);
      se_sq_recip = 1.0 / ((1.0 / case_a1_ctd) + (1.0 / ctrl_a2_ctd) + case_a2_recip + ctrl_a1_recip);
      x_assoc2 += se_sq_recip;
      dxx = ln_or * se_sq_recip;
      x_assoc1 += dxx;
      dxx *= ln_or;
      cluster_chisq[cluster_idx] = dxx;
      x_total += dxx;
    }
    x_assoc = x_assoc1 * x_assoc1 / x_assoc2;
    wptr_start = memcpyax(writebuf, chrom_name_ptr, chrom_name_len, ' ');
    wptr_start = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
    *wptr_start++ = ' ';
    wptr_start = fw_strcpy(4, marker_allele_ptrs[marker_uidx * 2], wptr_start);
    *wptr_start++ = ' ';
    wptr_start = fw_strcpy(4, marker_allele_ptrs[marker_uidx * 2 + 1], wptr_start);
    *wptr_start++ = ' ';
    wptr_start = memcpya(wptr_start, "      NA       NA       NA       NA ", 36);
    wptr = memcpya(wptr_start, " TOTAL ", 7);
    wptr = dtoa_g_wxp4x(x_total, 10, ' ', wptr);
    wptr = uint32toa_w4x(cluster_ct2, ' ', wptr);
    wptr = dtoa_g_wxp4x(chiprob_p(x_total, cluster_ct2d), 10, ' ', wptr);
    wptr = memcpya(wptr, "        NA\n", 11);
    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
      goto homog_assoc_ret_WRITE_FAIL;
    }
    wptr = memcpya(wptr_start, " ASSOC ", 7);
    wptr = dtoa_g_wxp4(x_assoc, 10, wptr);
    wptr = memcpya(wptr, "    1 ", 6);
    wptr = dtoa_g_wxp4x(chiprob_p(x_assoc, 1), 10, ' ', wptr);
    wptr = memcpya(wptr, "        NA\n", 11);
    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
      goto homog_assoc_ret_WRITE_FAIL;
    }
    dxx = x_total - x_assoc;
    wptr = memcpya(wptr_start, " HOMOG ", 7);
    wptr = dtoa_g_wxp4x(dxx, 10, ' ', wptr);
    wptr = uint32toa_w4x(cluster_ct2 - 1, ' ', wptr);
    wptr = dtoa_g_wxp4x(chiprob_p(dxx, cluster_ct2m1d), 10, ' ', wptr);
    wptr = memcpya(wptr, "        NA\n", 11);
    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
      goto homog_assoc_ret_WRITE_FAIL;
    }
    wptr_start = &(wptr_start[-36]);
    for (cluster_idx = 0, dptr = cluster_tables; cluster_idx < cluster_ct2; cluster_idx++, dptr = &(dptr[4])) {
      case_ctd = dptr[0] + dptr[1];
      ctrl_ctd = dptr[2] + dptr[3];
      if ((case_ctd < 1.5) || (ctrl_ctd < 1.5)) {
	wptr = memcpya(wptr_start, "      NA       NA ", 18);
	wptr = dtoa_g_wxp4x(case_ctd - 1, 8, ' ', wptr);
	wptr = dtoa_g_wxp4x(ctrl_ctd - 1, 8, ' ', wptr);
	wptr = fw_strcpy(6, &(cluster_ids_collapsed[cluster_idx * max_cluster_id_len]), wptr);
        wptr = memcpya(wptr, "         NA   NA         NA         NA\n", 39);
      } else {
        wptr = dtoa_g_wxp4x(dptr[0] / case_ctd, 8, ' ', wptr_start);
        wptr = dtoa_g_wxp4x(dptr[2] / ctrl_ctd, 8, ' ', wptr);
	wptr = dtoa_g_wxp4x(case_ctd - 1, 8, ' ', wptr);
	wptr = dtoa_g_wxp4x(ctrl_ctd - 1, 8, ' ', wptr);
	wptr = fw_strcpy(6, &(cluster_ids_collapsed[cluster_idx * max_cluster_id_len]), wptr);
	*wptr++ = ' ';
	dxx = cluster_chisq[cluster_idx];
	if (dxx < SMALL_EPSILON * SMALL_EPSILON) {
	  // probably rounding error
	  dxx = 0;
	}
        wptr = dtoa_g_wxp4(dxx, 10, wptr);
        wptr = memcpya(wptr, "    1 ", 6);
	wptr = dtoa_g_wxp4x(MAXV(chiprob_p(dxx, 1), output_min_p), 10, ' ', wptr);
	dxx = cluster_or[cluster_idx];
        if (realnum(dxx)) {
          wptr = dtoa_g_wxp4x(dxx, 10, '\n', wptr);
	} else {
	  wptr = memcpya(wptr, "        NA\n", 11);
	}
	if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	  goto homog_assoc_ret_WRITE_FAIL;
	}
      }
    }
    if (marker_idx >= loop_end) {
      if (marker_idx < marker_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
        pct = (marker_idx * 100LLU) / marker_ct;
        printf("\b\b%u%%", pct);
        fflush(stdout);
        loop_end = (((uint64_t)pct + 1LLU) * marker_ct) / 100;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto homog_assoc_ret_WRITE_FAIL;
  }
  if (pct >= 10) {
    putc_unlocked('\b', stdout);
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  while (0) {
  homog_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  homog_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  homog_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  homog_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 homog_assoc_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  fclose_cond(outfile);
  return retval;
}
