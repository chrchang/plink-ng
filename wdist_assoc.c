#include "wdist_cluster.h"
#include "wdist_matrix.h"
#include "wdist_stats.h"

// load markers in blocks to enable multithreading and, for quantitative
// phenotypes, PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 1024
#define MODEL_BLOCKKEEP 64

void single_marker_cc_freqs(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
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
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_6x_end;
  indiv_ctl2 -= indiv_ctl2 % 6;
  while (indiv_ctl2 >= 60) {
    cur_decr = 60;
  single_marker_cc_freqs_loop:
    lptr_6x_end = &(lptr[cur_decr]);
    count_2freq_dbl_60v((__m128i*)lptr, (__m128i*)lptr_6x_end, (__m128i*)ctrl_include2, (__m128i*)case_include2, &tot_ctrl_ab, &tot_ctrl_c, &tot_case_ab, &tot_case_c);
    lptr = lptr_6x_end;
    ctrl_include2 = &(ctrl_include2[cur_decr]);
    case_include2 = &(case_include2[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto single_marker_cc_freqs_loop;
  }
#else
  uintptr_t* lptr_six_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 6)]);
  while (lptr < lptr_six_end) {
    count_2freq_dbl_6(lptr, ctrl_include2, case_include2, &tot_ctrl_ab, &tot_ctrl_c, &tot_case_ab, &tot_case_c);
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

void haploid_single_marker_cc_freqs(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
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
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
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

void single_marker_cc_3freqs(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_hom2p, uint32_t* ctrl_hetp, uint32_t* ctrl_missingp, uint32_t* case_hom2p, uint32_t* case_hetp, uint32_t* case_missingp) {
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
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  indiv_ctl2 -= indiv_ctl2 % 12;
  while (indiv_ctl2 >= 120) {
    cur_decr = 120;
  single_marker_cc_3freqs_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)ctrl_include2, &tot_ctrl_a, &tot_ctrl_b, &tot_ctrl_c);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)case_include2, &tot_case_a, &tot_case_b, &tot_case_c);
    lptr = lptr_12x_end;
    ctrl_include2 = &(ctrl_include2[cur_decr]);
    case_include2 = &(case_include2[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto single_marker_cc_3freqs_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 12)]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, ctrl_include2, &tot_ctrl_a, &tot_ctrl_b, &tot_ctrl_c);
    count_3freq_12(lptr, case_include2, &tot_case_a, &tot_case_b, &tot_case_c);
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

static inline void adjust_print(FILE* outfile, double pval) {
  char wbuf[16];
  char* wpos;
  if (pval == 0) {
    fputs("       INF ", outfile);
  } else if (pval < 0) {
    fputs("        NA ", outfile);
  } else {
    wpos = double_g_writewx4x(wbuf, pval, 10, ' ');
    fwrite(wbuf, 1, wpos - wbuf, outfile);
  }
}

static inline void adjust_print_log10(FILE* outfile, double pval) {
  char wbuf[16];
  char* wpos;
  if (pval == 0) {
    fputs("       INF ", outfile);
  } else if (pval < 0) {
    fputs("        NA ", outfile);
  } else if (pval < 1) {
    wpos = double_g_writewx4x(wbuf, -log10(pval), 10, ' ');
    fwrite(wbuf, 1, wpos - wbuf, outfile);
  } else {
    fputs("         0 ", outfile);
  }
}

int32_t multcomp(char* outname, char* outname_end, uint32_t* marker_uidxs, uintptr_t chi_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, double* chi, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t non_chi, uint32_t* tcnt) {
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t adjust_gc = mtest_adjust & ADJUST_GC;
  uint32_t is_log10 = mtest_adjust & ADJUST_LOG10;
  uint32_t qq_plot = mtest_adjust & ADJUST_QQ;
  FILE* outfile = NULL;
  double pv_holm = 0.0;
  double pv_sidak_sd = 0;
  int32_t retval = 0;
  uint32_t uii = 0;
  uint32_t* new_tcnt = NULL;
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
  double* pv_gc;
  double lambda;
  double bonf;
  double pv_sidak_ss;
  char* bufptr;
  uint32_t ujj;
  uint32_t loop_end;

  if (wkspace_alloc_d_checked(&sp, chi_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&schi, chi_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&new_order, chi_ct * sizeof(uint32_t))) {
    goto multcomp_ret_NOMEM;
  }
  if (tcnt) {
    if (wkspace_alloc_ui_checked(&new_tcnt, chi_ct * sizeof(uint32_t))) {
      goto multcomp_ret_NOMEM;
    }
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      ujj = tcnt[cur_idx];
      if (ujj > 2) {
	dxx = chi[cur_idx]; // not actually squared
	dyy = calc_tprob(dxx, ujj - 2);
	if (dyy > -1) {
	  sp[uii] = dyy;
	  new_order[uii] = marker_uidxs[cur_idx];
	  schi[uii] = dxx * dxx;
	  new_tcnt[uii] = ujj - 2;
	  uii++;
	}
      }
    }
  } else if (non_chi) {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      dxx = 1 - chi[cur_idx];
      if (dxx < 1) {
	dyy = inverse_chiprob(dxx, 1);
	if (dyy >= 0) {
	  sp[uii] = dxx;
	  new_order[uii] = marker_uidxs[cur_idx];
	  schi[uii] = dyy;
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
  if (qsort_ext((char*)sp, chi_ct, sizeof(double), double_cmp_deref, (char*)new_order, sizeof(uint32_t))) {
    goto multcomp_ret_NOMEM;
  }
  if (tcnt) {
    if (qsort_ext((char*)schi, chi_ct, sizeof(double), double_cmp_deref, (char*)new_tcnt, sizeof(uint32_t))) {
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

  if (mtest_adjust & ADJUST_LAMBDA) {
    lambda = adjust_lambda;
  } else {
    if (chi_ct & 1) {
      lambda = schi[(chi_ct - 1) / 2];
    } else {
      lambda = (schi[chi_ct / 2 - 1] + schi[chi_ct / 2]) / 2.0;
    }
    lambda = lambda / 0.456;
    if (lambda < 1) {
      lambda = 1.0;
    } else {
      lambda = 1.0 / lambda;
    }
  }

  // handle reverse-order calculations
  if (wkspace_alloc_d_checked(&pv_bh, chi_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&pv_by, chi_ct * sizeof(double))) {
    goto multcomp_ret_NOMEM;
  }
  if (adjust_gc) {
    pv_gc = sp;
  } else {
    if (wkspace_alloc_d_checked(&pv_gc, chi_ct * sizeof(double))) {
      goto multcomp_ret_NOMEM;
    }
  }
  uii = chi_ct;
  if (tcnt) {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      uii--;
      pv_gc[cur_idx] = calc_tprob(sqrt(schi[uii] * lambda), new_tcnt[uii]);
    }
  } else {
    for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
      pv_gc[cur_idx] = chiprob_p(schi[--uii] * lambda, 1);
    }
  }

  dyy = sp[chi_ct - 1];
  pv_bh[chi_ct - 1] = dyy;
  harmonic_sum = 1.0;
  for (cur_idx = chi_ct - 1; cur_idx > 0; cur_idx--) {
    dzz = dct / ((double)cur_idx);
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
    dxx = (harmonic_sum / ((double)cur_idx)) * sp[cur_idx - 1];
    if (dyy > dxx) {
      dyy = dxx;
    }
    pv_by[cur_idx - 1] = dyy;
  }

  uii = strlen(outname_end);
  memcpy(&(outname_end[uii]), ".adjusted", 10);
  if (fopen_checked(&outfile, outname, "w")) {
    goto multcomp_ret_OPEN_FAIL;
  }
  uii = sprintf(tbuf, " CHR %%%us      UNADJ         GC ", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  if (qq_plot) {
    fputs("        QQ ", outfile);
  }
  if (fputs_checked("      BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY\n", outfile)) {
    goto multcomp_ret_WRITE_FAIL;
  }
  fputs("0%", stdout);
  fflush(stdout);
  cur_idx = 0;
  for (pct = 1; pct <= 100; pct++) {
    loop_end = (((uint64_t)pct) * chi_ct) / 100LLU;
    for (; cur_idx < loop_end; cur_idx++) {
      pval = sp[cur_idx];
      if (pval > pfilter) {
	continue;
      }
      marker_uidx = new_order[cur_idx];
      bufptr = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
      *bufptr++ = ' ';
      bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
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
      pv_sidak_ss = 1 - pow(1 - pval, dct);
      dyy = 1 - pow(1 - pval, dct - ((double)cur_idx));
      if (pv_sidak_sd < dyy) {
	pv_sidak_sd = dyy;
      }

      if (!is_log10) {
	adjust_print(outfile, pval);
	adjust_print(outfile, pv_gc[cur_idx]);
	if (qq_plot) {
	  adjust_print(outfile, (((double)cur_idx) + 0.5) * dzz);
	}
	adjust_print(outfile, bonf);
	adjust_print(outfile, pv_holm);
	adjust_print(outfile, pv_sidak_ss);
	adjust_print(outfile, pv_sidak_sd);
	adjust_print(outfile, pv_bh[cur_idx]);
	adjust_print(outfile, pv_by[cur_idx]);
      } else {
	adjust_print_log10(outfile, pval);
	adjust_print_log10(outfile, pv_gc[cur_idx]);
	if (qq_plot) {
	  adjust_print_log10(outfile, (((double)cur_idx) + 0.5) * dzz);
	}
	adjust_print_log10(outfile, bonf);
	adjust_print_log10(outfile, pv_holm);
	adjust_print_log10(outfile, pv_sidak_ss);
	adjust_print_log10(outfile, pv_sidak_sd);
	adjust_print_log10(outfile, pv_bh[cur_idx]);
	adjust_print_log10(outfile, pv_by[cur_idx]);
      }
      if (putc_checked('\n', outfile)) {
	goto multcomp_ret_WRITE_FAIL;
      }
    }
    if (pct < 100) {
      if (pct > 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  fputs("\b\b\b", stdout);
  sprintf(logbuf, "--adjust values (%" PRIuPTR " marker%s) written to %s.\n", chi_ct, (chi_ct == 1)? "" : "s", outname);
  logprintb();

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
  wkspace_reset(wkspace_mark);
  return retval;
}

void generate_cc_perm_vec(uint32_t tot_ct, uint32_t set_ct, uint32_t tot_quotient, uint64_t totq_magic, uint32_t totq_preshift, uint32_t totq_postshift, uint32_t totq_incr, uintptr_t* perm_vec, sfmt_t* sfmtp) {
  // Assumes tot_quotient is 4294967296 / tot_ct, and
  // totq_magic/totq_preshift/totq_postshift/totq_incr have been precomputed
  // from magic_num().
  uint32_t num_set = 0;
  uint32_t upper_bound = tot_ct * tot_quotient - 1;
  uintptr_t widx;
  uintptr_t wcomp;
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  if (set_ct * 2 < tot_ct) {
    fill_ulong_zero(perm_vec, 2 * ((tot_ct + (BITCT - 1)) / BITCT));
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT2;
	wcomp = ONELU << (2 * (uii % BITCT2));
	pv_val = perm_vec[widx];
      } while (pv_val & wcomp);
      perm_vec[widx] = pv_val | wcomp;
    }
  } else {
    fill_vec_55(perm_vec, tot_ct);
    set_ct = tot_ct - set_ct;
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT2;
	wcomp = ONELU << (2 * (uii % BITCT2));
	pv_val = perm_vec[widx];
      } while (!(pv_val & wcomp));
      perm_vec[widx] = pv_val - wcomp;
    }
  }
}

void generate_cc_cluster_perm_vec(uint32_t tot_ct, uintptr_t* preimage, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uint32_t* tot_quotients, uint64_t* totq_magics, uint32_t* totq_preshifts, uint32_t* totq_postshifts, uint32_t* totq_incrs, uintptr_t* perm_vec, sfmt_t* sfmtp) {
  uint32_t tot_ctl2 = 2 * ((tot_ct + (BITCT - 1)) / BITCT);
  uint32_t cluster_idx;
  uint32_t target_ct;
  uint32_t cluster_end;
  uint32_t cluster_size;
  uint32_t* map_ptr;
  uint32_t num_swapped;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uintptr_t widx;
  uintptr_t wcomp;
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  memcpy(perm_vec, preimage, tot_ctl2 * sizeof(intptr_t));
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    target_ct = cluster_case_cts[cluster_idx];
    cluster_end = cluster_starts[cluster_idx + 1];
    cluster_size = cluster_end - cluster_starts[cluster_idx];
    if (target_ct && (target_ct != cluster_size)) {
      upper_bound = tot_ct * tot_quotients[cluster_idx] - 1;
      totq_magic = totq_magics[cluster_idx];
      totq_preshift = totq_preshifts[cluster_idx];
      totq_postshift = totq_postshifts[cluster_idx];
      totq_incr = totq_incrs[cluster_idx];
      map_ptr = &(cluster_map[cluster_starts[cluster_idx]]);
      if (target_ct * 2 < cluster_size) {
	for (num_swapped = 0; num_swapped < target_ct; num_swapped++) {
	  do {
	    do {
	      urand = sfmt_genrand_uint32(sfmtp);
	    } while (urand > upper_bound);
	    uii = map_ptr[(uint32_t)((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift)];
	    widx = uii / BITCT2;
	    wcomp = ONELU << (2 * (uii % BITCT2));
	    pv_val = perm_vec[widx];
	  } while (pv_val & wcomp);
	  perm_vec[widx] = pv_val | wcomp;
	}
      } else {
	target_ct = cluster_size - target_ct;
	for (num_swapped = 0; num_swapped < target_ct; num_swapped++) {
	  do {
	    do {
	      urand = sfmt_genrand_uint32(sfmtp);
	    } while (urand > upper_bound);
	    uii = map_ptr[(uint32_t)((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift)];
	    widx = uii / BITCT2;
	    wcomp = ONELU << (2 * (uii % BITCT2));
	    pv_val = perm_vec[widx];
	  } while (!(pv_val & wcomp));
	  perm_vec[widx] = pv_val - wcomp;
	}
      }
    }
  }
}

void transpose_perms(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uint32_t* perm_vecst) {
  // Transpose permutations so PRESTO/PERMORY-style genotype indexing can work.
  //
  // We used a 32-ply interleaved format, to allow counts up to the uint32_t
  // limit without giving up highly parallel adds in the calc_git() inner loop.
  // The index order used here is:
  // 64-bit build:
  //   first 16 bytes: 0 32 64 96 16 48 80 112 4 36 68 100 20 52 84 116
  //     8 40 72 104 24 56 88 120 12 44 76 108 28 60 92 124 1...
  //   next 16 bytes: 128 160 192...
  //
  // 32-bit build:
  //   first 4 bytes: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18...
  //   next 4 bytes: 32 40 48...
  uintptr_t indiv_idx = 0;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  uint32_t wbuf[4];
  uint32_t* wbptr;
#else
  uint32_t wval;
#endif
  uint32_t rshift;
  uint32_t wshift;
  uintptr_t* pvptr;
  uintptr_t perm_idx;
  for (; indiv_idx < pheno_nm_ct; indiv_idx++) {
    perm_idx = 0;
    pvptr = &(perm_vecs[indiv_idx / BITCT2]);
    rshift = 2 * (indiv_idx % BITCT2);
    goto transpose_perms_loop_start;
#ifdef __LP64__
    do {
      if (!(perm_idx % 4)) {
	if (perm_idx % 128) {
	  wshift = ((perm_idx & 96) >> 5) | ((perm_idx & 16) >> 2) | ((perm_idx & 12) << 1);
	} else {
	  memcpy(perm_vecst, wbuf, 16);
	  perm_vecst = &(perm_vecst[4]);
	transpose_perms_loop_start:
	  fill_ulong_zero((uintptr_t*)wbuf, 2);
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctl2] >> rshift) & 1) << wshift;
      wbptr++;
    } while (++perm_idx < perm_vec_ct);
    memcpy(perm_vecst, wbuf, 16);
    perm_vecst = &(perm_vecst[4]);
#else
    do {
      if (perm_idx % 32) {
	wshift = ((perm_idx & 24) >> 3) | (perm_idx & 4) | ((perm_idx & 3) << 3);
      } else {
	*perm_vecst++ = wval;
      transpose_perms_loop_start:
	wval = 0;
	wshift = 0;
      }
      wval |= ((pvptr[perm_idx * pheno_nm_ctl2] >> rshift) & 1) << wshift;
    } while (++perm_idx < perm_vec_ct);
    *perm_vecst++ = wval;
#endif
  }
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
  // 1. For the first 15 heterozygote individuals, just use 4-bit accumulators.
  //    This allows the inner loop to increment 32 counters simultaneously.
  // 2. Right before they'd otherwise be at risk of overflowing, we unfold the
  //    4-bit accumulators into a larger buffer of 8-bit accumulators.  Then we
  //    zero out the 4-bit accumulators, and restart the inner loop.
  // 3. This can happen up to 17 times before the 8-bit accumulators risk
  //    overflow.  Then, they are unfolded into the final output array of
  //    32-bit ints, zeroed out, and the second loop restarts.
  // Note that thread_bufs[] is assumed to be zeroed out before this function
  // is called.
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct128x4 = perm_ct128 * 4;
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct16x4 = 4 * perm_ct16;
  const __m128i m1x4 = {0x1111111111111111LLU, 0x1111111111111111LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8x32 = {0x000000ff000000ffLLU, 0x000000ff000000ffLLU};
  __m128i* permsv = (__m128i*)perm_vecst;
  __m128i* gitv[9];
  __m128i* __restrict__ git_merge4; // no conflicts, please
  __m128i* __restrict__ git_merge8;
  __m128i* __restrict__ git_write;
  __m128i* __restrict__ perm_ptr;
  __m128i loader;
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct32x4 = perm_ct32 * 4;
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t perm_ct16x16 = 16 * perm_ct16;
  uint32_t* permsv = perm_vecst;
  uint32_t* gitv[9];
  uint32_t* git_merge4;
  uint32_t* git_merge8;
  uint32_t* git_write;
  uint32_t* perm_ptr;
  uintptr_t loader;
#endif
  uint32_t cur_cts[3];
  uintptr_t ulii;
  uint32_t pbidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t indiv_type;
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
  gitv[0] = thread_wkspace;
  gitv[1] = &(thread_wkspace[perm_ct32x4]);
  gitv[2] = &(thread_wkspace[2 * perm_ct32x4]);
  gitv[3] = &(thread_wkspace[3 * perm_ct32x4]);
  gitv[4] = &(thread_wkspace[3 * perm_ct32x4 + 2 * perm_ct8]);
  gitv[5] = &(thread_wkspace[3 * perm_ct32x4 + 4 * perm_ct8]);
  gitv[6] = &(results_bufs[2 * perm_ct16x16]);
  gitv[7] = &(results_bufs[perm_ct16x16]);
  gitv[8] = results_bufs;
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
      ujj = CTZLU(ulii) & (BITCT - 2); // get pos of next non-[hom A2] indiv
      indiv_type = ((ulii >> ujj) & 3) - 1;
      git_merge4 = gitv[indiv_type];
#ifdef __LP64__
      perm_ptr = &(permsv[(ujj / 2) * perm_ct128]);
      for (pbidx = 0; pbidx < perm_ct128; pbidx++) {
	loader = *perm_ptr++;
	git_merge4[0] = _mm_add_epi64(git_merge4[0], _mm_and_si128(loader, m1x4));
	git_merge4[1] = _mm_add_epi64(git_merge4[1], _mm_and_si128(_mm_srli_epi64(loader, 1), m1x4));
	git_merge4[2] = _mm_add_epi64(git_merge4[2], _mm_and_si128(_mm_srli_epi64(loader, 2), m1x4));
	git_merge4[3] = _mm_add_epi64(git_merge4[3], _mm_and_si128(_mm_srli_epi64(loader, 3), m1x4));
	git_merge4 = &(git_merge4[4]);
      }
      ukk = cur_cts[indiv_type] + 1;
      cur_cts[indiv_type] = ukk;
      if (!(ukk % 15)) {
	git_merge4 = gitv[indiv_type];
	git_merge8 = gitv[indiv_type + 3];
	for (pbidx = 0; pbidx < perm_ct32; pbidx++) {
	  loader = *git_merge4;
	  git_merge8[0] = _mm_add_epi64(git_merge8[0], _mm_and_si128(loader, m4));
	  git_merge8[1] = _mm_add_epi64(git_merge8[1], _mm_and_si128(_mm_srli_epi64(loader, 4), m4));
	  git_merge8 = &(git_merge8[2]);
	  *git_merge4++ = _mm_setzero_si128();
	}
	if (!(ukk % 255)) {
	  git_merge8 = gitv[indiv_type + 3];
	  git_write = gitv[indiv_type + 6];
	  for (pbidx = 0; pbidx < perm_ct16; pbidx++) {
	    loader = *git_merge8;
	    git_write[0] = _mm_add_epi64(git_write[0], _mm_and_si128(loader, m8x32));
	    git_write[1] = _mm_add_epi64(git_write[1], _mm_and_si128(_mm_srli_epi64(loader, 8), m8x32));
	    git_write[2] = _mm_add_epi64(git_write[2], _mm_and_si128(_mm_srli_epi64(loader, 16), m8x32));
	    git_write[3] = _mm_add_epi64(git_write[3], _mm_and_si128(_mm_srli_epi64(loader, 24), m8x32));
	    git_write = &(git_write[4]);
	    *git_merge8++ = _mm_setzero_si128();
	  }
	}
      }
#else
      perm_ptr = &(permsv[(ujj / 2) * perm_ct32]);
      for (pbidx = 0; pbidx < perm_ct32; pbidx++) {
	loader = *perm_ptr++;
	git_merge4[0] += loader & 0x11111111;
	git_merge4[1] += (loader >> 1) & 0x11111111;
	git_merge4[2] += (loader >> 2) & 0x11111111;
	git_merge4[3] += (loader >> 3) & 0x11111111;
	git_merge4 = &(git_merge4[4]);
      }
      ukk = cur_cts[indiv_type] + 1;
      cur_cts[indiv_type] = ukk;
      if (!(ukk % 15)) {
	git_merge4 = gitv[indiv_type];
	git_merge8 = gitv[indiv_type + 3];
	for (pbidx = 0; pbidx < perm_ct8; pbidx++) {
	  loader = *git_merge4;
	  git_merge8[0] += loader & 0x0f0f0f0f;
	  git_merge8[1] += (loader >> 4) & 0x0f0f0f0f;
	  git_merge8 = &(git_merge8[2]);
	  *git_merge4++ = 0;
	}
	if (!(ukk % 255)) {
	  git_merge8 = gitv[indiv_type + 3];
	  git_write = gitv[indiv_type + 6];
	  for (pbidx = 0; pbidx < perm_ct4; pbidx++) {
	    loader = *git_merge8;
	    git_write[0] += loader & 0x000000ff;
	    git_write[1] += (loader >> 8) & 0x000000ff;
	    git_write[2] += (loader >> 16) & 0x000000ff;
	    git_write[3] += loader >> 24;
	    git_write = &(git_write[4]);
	    *git_merge8++ = 0;
	  }
	}
      }
#endif
      ulii &= ~(3 * (ONELU << ujj));
    }
#ifdef __LP64__
    permsv = &(permsv[BITCT2 * perm_ct128]);
#else
    permsv = &(permsv[BITCT2 * perm_ct32]);
#endif
  }
  for (indiv_type = 0; indiv_type < 3; indiv_type++) {
    uii = cur_cts[indiv_type];
#ifdef __LP64__
    if (uii % 15) {
      git_merge4 = gitv[indiv_type];
      git_merge8 = gitv[indiv_type + 3];
      for (pbidx = 0; pbidx < perm_ct32; pbidx++) {
	loader = *git_merge4++;
	git_merge8[0] = _mm_add_epi64(git_merge8[0], _mm_and_si128(loader, m4));
	git_merge8[1] = _mm_add_epi64(git_merge8[1], _mm_and_si128(_mm_srli_epi64(loader, 4), m4));
	git_merge8 = &(git_merge8[2]);
      }
    }
    if (uii % 255) {
      git_merge8 = gitv[indiv_type + 3];
      git_write = gitv[indiv_type + 6];
      for (pbidx = 0; pbidx < perm_ct16; pbidx++) {
	loader = *git_merge8++;
	git_write[0] = _mm_add_epi64(git_write[0], _mm_and_si128(loader, m8x32));
	git_write[1] = _mm_add_epi64(git_write[1], _mm_and_si128(_mm_srli_epi64(loader, 8), m8x32));
	git_write[2] = _mm_add_epi64(git_write[2], _mm_and_si128(_mm_srli_epi64(loader, 16), m8x32));
	git_write[3] = _mm_add_epi64(git_write[3], _mm_and_si128(_mm_srli_epi64(loader, 24), m8x32));
	git_write = &(git_write[4]);
      }
    }
#else
    if (uii % 15) {
      git_merge4 = gitv[indiv_type];
      git_merge8 = gitv[indiv_type + 3];
      for (pbidx = 0; pbidx < perm_ct8; pbidx++) {
	loader = *git_merge4++;
	git_merge8[0] += loader & 0x0f0f0f0f;
	git_merge8[1] += (loader >> 4) & 0x0f0f0f0f;
	git_merge8 = &(git_merge8[2]);
      }
    }
    if (uii % 255) {
      git_merge8 = gitv[indiv_type + 3];
      git_write = gitv[indiv_type + 6];
      for (pbidx = 0; pbidx < perm_ct4; pbidx++) {
	loader = *git_merge8++;
	git_write[0] += loader & 0x000000ff;
	git_write[1] += (loader >> 8) & 0x000000ff;
	git_write[2] += (loader >> 16) & 0x000000ff;
	git_write[3] += loader >> 24;
	git_write = &(git_write[4]);
      }
    }
#endif
  }
}

void calc_qgit(uint32_t pheno_nm_ct, uintptr_t perm_vec_ctcl8m, uint32_t num_perms_now, uintptr_t* __restrict__ loadbuf, double* perm_vecstd, double* thread_bufs) {
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the indiv idx
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
  uint32_t indiv_type;
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
      indiv_type = (ulii >> ujj) & 3;
#ifdef __LP64__
      // note that the gain from using SSE2 for double-precision arithmetic is
      // typically minimal because modern cores tend to have two FPUs, so we
      // should only use it opportunistically.  it's painless here, though.
      perm_readv = &(permsv[ujj * row_mult]);
      if (indiv_type == 1) {
	git_writev = (__m128d*)thread_bufs;
	for (ukk = 0; ukk < loop_len; ukk++) {
	  *git_writev = _mm_add_pd(*git_writev, *perm_readv++);
	  git_writev++;
	}
      } else if (indiv_type == 3) {
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
      if (indiv_type == 1) {
	git_write = thread_bufs;
	for (ukk = 0; ukk < num_perms_now; ukk++) {
	  *git_write += *perm_read++;
	  git_write++;
	}
      } else if (indiv_type == 3) {
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
      ulii &= ~(3 * (ONELU << ujj));
    }
#ifdef __LP64__
    permsv = &(permsv[(BITCT2 / 2) * perm_vec_ctcl8m]);
#else
    perm_vecstd = &(perm_vecstd[BITCT2 * perm_vec_ctcl8m]);
#endif
  }
}

void calc_qgit_lin(uint32_t pheno_nm_ct, uintptr_t perm_vec_ctcl8m, uint32_t num_perms_now, uintptr_t* __restrict__ loadbuf, double* perm_vecstd, double* thread_bufs) {
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the indiv idx
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
  uint32_t indiv_type;
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
      indiv_type = (ulii >> ujj) & 3;
#ifdef __LP64__
      perm_readv = &(permsv[ujj * row_mult]);
      if (indiv_type == 1) {
	git_writev = (__m128d*)thread_bufs;
	git_write2v = (__m128d*)(&(thread_bufs[perm_vec_ctcl8m]));
      } else if (indiv_type == 3) {
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
      if (indiv_type == 1) {
	git_write = thread_bufs;
	git_write2 = &(thread_bufs[perm_vec_ctcl8m]);
      } else if (indiv_type == 3) {
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
      ulii &= ~(3 * (ONELU << ujj));
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
  __uni16 acc;
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
  __uni16 acc;
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

uintptr_t rem_cost(uintptr_t indiv_ctl2, uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
  // Cost: 2 * (<-> neither side homcom) + (<-> homcom)
  //
  // We can efficiently calculate this as follows:
  //   xor = vec1 ^ vec2
  //   detect_homcom = (vec1 & (vec1 >> 1)) | (vec2 & (vec2 >> 1))
  //   A := (xor | (xor >> 1)) & 0x5555...
  //   B := A & (~detect_homcom)
  //   cost += popcount2(A + B)
  uintptr_t* lptr_end = &(loadbuf1[indiv_ctl2]);
  uintptr_t cost = 0;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t xor_word;
  uintptr_t detect_homcom;
  uintptr_t result_a;
  uintptr_t result_b;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_6x_end;
  indiv_ctl2 -= indiv_ctl2 % 6;
  while (indiv_ctl2 >= 60) {
    cur_decr = 60;
  rem_cost_loop:
    lptr_6x_end = &(loadbuf1[cur_decr]);
    cost += rem_cost_60v((__m128i*)loadbuf1, (__m128i*)lptr_6x_end, (__m128i*)loadbuf2);
    loadbuf1 = lptr_6x_end;
    loadbuf2 = &(loadbuf2[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto rem_cost_loop;
  }
#else
  uintptr_t* lptr_six_end = &(loadbuf1[indiv_ctl2 - (indiv_ctl2 % 6)]);
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

uintptr_t qrem_cost2(uintptr_t indiv_ctl2, uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
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
  uintptr_t* lptr_end = &(loadbuf1[indiv_ctl2]);
  uintptr_t cost = 3;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t xor_word;
  uintptr_t detect_missing;
  uintptr_t result_a;
  uintptr_t result_b;
  uintptr_t result_c;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_4x_end;
  indiv_ctl2 &= ~3LLU;
  while (indiv_ctl2 >= 40) {
    cur_decr = 40;
  qrem_cost2_loop:
    lptr_4x_end = &(loadbuf1[cur_decr]);
    cost += qrem_cost2_40v((__m128i*)loadbuf1, (__m128i*)lptr_4x_end, (__m128i*)loadbuf2);
    loadbuf1 = lptr_4x_end;
    loadbuf2 = &(loadbuf2[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto qrem_cost2_loop;
  }
#else
  uintptr_t* lptr_four_end = &(loadbuf1[indiv_ctl2 & (~3U)]);
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
static inline void calc_rem_merge4_one(uint32_t perm_ct128, __m128i* __restrict__ perm_ptr, __m128i* __restrict__ rem_merge4) {
  const __m128i m1x4 = {0x1111111111111111LLU, 0x1111111111111111LLU};
  __m128i loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct128; pbidx++) {
    loader = *perm_ptr++;
    rem_merge4[0] = _mm_add_epi64(rem_merge4[0], _mm_and_si128(loader, m1x4));
    rem_merge4[1] = _mm_add_epi64(rem_merge4[1], _mm_and_si128(_mm_srli_epi64(loader, 1), m1x4));
    rem_merge4[2] = _mm_add_epi64(rem_merge4[2], _mm_and_si128(_mm_srli_epi64(loader, 2), m1x4));
    rem_merge4[3] = _mm_add_epi64(rem_merge4[3], _mm_and_si128(_mm_srli_epi64(loader, 3), m1x4));
    rem_merge4 = &(rem_merge4[4]);
  }
}

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
    loader2 = _mm_and_si128(_mm_srli_epi64(loader, 1), m1x4);
    rem_merge4a[1] = _mm_add_epi64(rem_merge4a[1], loader2);
    rem_merge4b[1] = _mm_add_epi64(rem_merge4b[1], loader2);
    loader2 = _mm_and_si128(_mm_srli_epi64(loader, 2), m1x4);
    rem_merge4a[2] = _mm_add_epi64(rem_merge4a[2], loader2);
    rem_merge4b[2] = _mm_add_epi64(rem_merge4b[2], loader2);
    loader2 = _mm_and_si128(_mm_srli_epi64(loader, 3), m1x4);
    rem_merge4a[3] = _mm_add_epi64(rem_merge4a[3], loader2);
    rem_merge4b[3] = _mm_add_epi64(rem_merge4b[3], loader2);
    rem_merge4a = &(rem_merge4a[4]);
    rem_merge4b = &(rem_merge4b[4]);
  }
}

static inline void calc_rem_merge8(uint32_t perm_ct32, __m128i* __restrict__ rem_merge4, __m128i* __restrict__ rem_merge8) {
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct32; pbidx++) {
    loader = *rem_merge4;
    rem_merge8[0] = _mm_add_epi64(rem_merge8[0], _mm_and_si128(loader, m4));
    rem_merge8[1] = _mm_add_epi64(rem_merge8[1], _mm_and_si128(_mm_srli_epi64(loader, 4), m4));
    rem_merge8 = &(rem_merge8[2]);
    *rem_merge4++ = _mm_setzero_si128();
  }
}

static inline void calc_rem_merge32_plus(uint32_t perm_ct16, __m128i* __restrict__ rem_merge8, __m128i* rem_write) {
  const __m128i m8x32 = {0x000000ff000000ffLLU, 0x000000ff000000ffLLU};
  __m128i loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct16; pbidx++) {
    loader = *rem_merge8;
    rem_write[0] = _mm_add_epi64(rem_write[0], _mm_and_si128(loader, m8x32));
    rem_write[1] = _mm_add_epi64(rem_write[1], _mm_and_si128(_mm_srli_epi64(loader, 8), m8x32));
    rem_write[2] = _mm_add_epi64(rem_write[2], _mm_and_si128(_mm_srli_epi64(loader, 16), m8x32));
    rem_write[3] = _mm_add_epi64(rem_write[3], _mm_and_si128(_mm_srli_epi64(loader, 24), m8x32));
    rem_write = &(rem_write[4]);
    *rem_merge8++ = _mm_setzero_si128();
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
    rem_write[1] = _mm_sub_epi64(rem_write[1], _mm_and_si128(_mm_srli_epi64(loader, 8), m8x32));
    rem_write[2] = _mm_sub_epi64(rem_write[2], _mm_and_si128(_mm_srli_epi64(loader, 16), m8x32));
    rem_write[3] = _mm_sub_epi64(rem_write[3], _mm_and_si128(_mm_srli_epi64(loader, 24), m8x32));
    rem_write = &(rem_write[4]);
    *rem_merge8++ = _mm_setzero_si128();
  }
}
#else
static inline void calc_rem_merge4_one(uint32_t perm_ct32, uintptr_t* __restrict__ perm_ptr, uintptr_t* __restrict__ rem_merge4) {
  uintptr_t loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct32; pbidx++) {
    loader = *perm_ptr++;
    rem_merge4[0] += loader & 0x11111111;
    rem_merge4[1] += (loader >> 1) & 0x11111111;
    rem_merge4[2] += (loader >> 2) & 0x11111111;
    rem_merge4[3] += (loader >> 3) & 0x11111111;
    rem_merge4 = &(rem_merge4[4]);
  }
}

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

static inline void calc_rem_merge8(uint32_t perm_ct8, uintptr_t* __restrict__ rem_merge4, uintptr_t* __restrict__ rem_merge8) {
  uintptr_t loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct8; pbidx++) {
    loader = *rem_merge4;
    rem_merge8[0] += loader & 0x0f0f0f0f;
    rem_merge8[1] += (loader >> 4) & 0x0f0f0f0f;
    rem_merge8 = &(rem_merge8[2]);
    *rem_merge4++ = 0;
  }
}

static inline void calc_rem_merge32_plus(uint32_t perm_ct4, uintptr_t* __restrict__ rem_merge8, uintptr_t* __restrict__ rem_write) {
  uintptr_t loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct4; pbidx++) {
    loader = *rem_merge8;
    rem_write[0] += loader & 0x000000ff;
    rem_write[1] += (loader >> 8) & 0x000000ff;
    rem_write[2] += (loader >> 16) & 0x000000ff;
    rem_write[3] += loader >> 24;
    rem_write = &(rem_write[4]);
    *rem_merge8++ = 0;
  }
}

static inline void calc_rem_merge32_minus(uint32_t perm_ct4, uintptr_t* __restrict__ rem_merge8, uintptr_t* __restrict__ rem_write) {
  uintptr_t loader;
  uint32_t pbidx;
  for (pbidx = 0; pbidx < perm_ct4; pbidx++) {
    loader = *rem_merge8;
    rem_write[0] -= loader & 0x000000ff;
    rem_write[1] -= (loader >> 8) & 0x000000ff;
    rem_write[2] -= (loader >> 16) & 0x000000ff;
    rem_write[3] -= loader >> 24;
    rem_write = &(rem_write[4]);
    *rem_merge8++ = 0;
  }
}
#endif

void calc_rem(uint32_t pheno_nm_ct, uintptr_t perm_vec_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_ref, uint32_t* perm_vecst, uint32_t* results_bufs, uint32_t* thread_wkspace) {
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
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
	calc_rem_merge4_one(perm_ct128, perm_ptr, remv[idx1]);
      } else {
	calc_rem_merge4_two(perm_ct128, perm_ptr, remv[idx1], remv[idx2]);
	ukk = cur_cts[idx2] + 1;
	cur_cts[idx2] = ukk;
	if (!(ukk % 15)) {
	  calc_rem_merge8(perm_ct32, remv[idx2], remv[idx2 + 6]);
	  if (!(ukk % 255)) {
	    calc_rem_merge32_minus(perm_ct16, remv[idx2 + 6], remv[(idx2 / 2) + 12]);
	  }
	}
      }
      ukk = cur_cts[idx1] + 1;
      cur_cts[idx1] = ukk;
      if (!(ukk % 15)) {
	calc_rem_merge8(perm_ct32, remv[idx1], remv[idx1 + 6]);
	if (!(ukk % 255)) {
	  if (!(idx1 & 1)) {
	    calc_rem_merge32_plus(perm_ct16, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
	  } else {
	    calc_rem_merge32_minus(perm_ct16, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
	  }
	}
      }
#else
      perm_ptr = &(permsv[(ujj / 2) * perm_ct32]);
      if (!idx2) {
	calc_rem_merge4_one(perm_ct32, perm_ptr, remv[idx1]);
      } else {
	calc_rem_merge4_two(perm_ct32, perm_ptr, remv[idx1], remv[idx2]);
	ukk = cur_cts[idx2] + 1;
	cur_cts[idx2] = ukk;
	if (!(ukk % 15)) {
	  calc_rem_merge8(perm_ct8, remv[idx2], remv[idx2 + 6]);
	  if (!(ukk % 255)) {
	    calc_rem_merge32_minus(perm_ct4, remv[idx2 + 6], remv[(idx2 / 2) + 12]);
	  }
	}
      }
      ukk = cur_cts[idx1] + 1;
      cur_cts[idx1] = ukk;
      if (!(ukk % 15)) {
	calc_rem_merge8(perm_ct8, remv[idx1], remv[idx1 + 6]);
	if (!(ukk % 255)) {
	  if (!(idx1 & 1)) {
	    calc_rem_merge32_plus(perm_ct4, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
	  } else {
	    calc_rem_merge32_minus(perm_ct4, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
	  }
	}
      }
#endif
      ulxor &= ~(3 * (ONELU << ujj));
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
      calc_rem_merge8(perm_ct32, remv[idx1], remv[idx1 + 6]);
    }
    if (uii % 255) {
      if (!(idx1 & 1)) {
	calc_rem_merge32_plus(perm_ct16, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
      } else {
	calc_rem_merge32_minus(perm_ct16, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
      }
    }
#else
    if (uii % 15) {
      calc_rem_merge8(perm_ct8, remv[idx1], remv[idx1 + 6]);
    }
    if (uii % 255) {
      if (!(idx1 & 1)) {
	calc_rem_merge32_plus(perm_ct4, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
      } else {
	calc_rem_merge32_minus(perm_ct4, remv[idx1 + 6], remv[(idx1 / 2) + 12]);
      }
    }
#endif
  }
}

void calc_qrem(uint32_t pheno_nm_ct, uintptr_t perm_vec_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_ref, double* perm_vecstd, double* outbufs) {
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the indiv idx
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
      ulxor &= ~(3 * (ONELU << ujj));
    }
#ifdef __LP64__
    permsv = &(permsv[(BITCT2 / 2) * perm_vec_ctcl8m]);
#else
    perm_vecstd = &(perm_vecstd[BITCT2 * perm_vec_ctcl8m]);
#endif
  }
}

void calc_qrem_lin(uint32_t pheno_nm_ct, uintptr_t perm_vec_ct, uintptr_t* loadbuf, uintptr_t* loadbuf_ref, double* perm_vecstd, double* outbufs) {
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
#ifdef __LP64__
  // halve for 8 bytes vs. 16, halve again for ujj being double the indiv idx
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
      ulxor &= ~(3 * (ONELU << ujj));
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
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
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
	  cur_cost = rem_cost(pheno_nm_ctl2, &(loadbuf[marker_bidx2 * pheno_nm_ctl2]), loadbuf_cur);
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
static double* g_orig_1mpval;
static double* g_orig_chisq;
static double* g_mperm_save_all;

// A separated-low-and-high-bit format was tried, and found to not really be
// any better than the usual PLINK 2-bit format.
static uintptr_t* g_loadbuf;

static uintptr_t* g_perm_vecs;

static uint32_t* g_perm_vecst; // genotype indexing support
static uint32_t* g_thread_git_wkspace;
static uint32_t* g_resultbuf;

// always use genotype indexing for quantitative traits
static double* g_perm_vecstd;
static double* g_thread_git_qbufs;
static double* g_qresultbuf;
static double* g_pheno_d2;
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
// For --assoc perm/--model [dom/rec/trend] perm:
//   g_precomp_ui[4n] and [4n + 1] define the interval with less extreme
//     p-values than the original.  [4n + 2] and [4n + 3] define the
//     interval with less or equally extreme p-values.
//
// For --assoc mperm fisher/--model [dom/rec] fisher:
//   g_precomp_ui[6n]...[6n + 3] is as in --assoc perm.
//   g_precomp_ui[6n + 4] and [6n + 5] are the floor and offset for the
//     range of case_set_cts where Fisher p-value calculation is unnecessary.
//   g_precomp_d[3n], [3n + 1], and [3n + 2] are tot_prob, right_prob, and
//     tail_sum, respectively, for fisher22_tail_pval().  (This is almost
//     irrelevant.)
//
// For --assoc mperm/--model [dom/rec/trend] mperm:
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
//   g_precomp_d[9n] to [9n + 2] are fisher22_tail_pval() coefficients for the
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
// elsewhere: number of missing individuals for each marker
static uint32_t* g_missing_cts;

static uint32_t* g_set_cts;
static uint32_t* g_het_cts;
static uint32_t* g_homcom_cts;

// This is *twice* the number of successes, because PLINK counts tie as 0.5.
// (Actually, PLINK randomizes instead of deterministically adding 0.5; this
// randomization just adds noise so we don't replicate it.)
static uint32_t* g_perm_2success_ct;
static uint32_t* g_perm_attempt_ct;
static double* g_maxt_extreme_stat;
static double* g_maxt_thread_results;
static double g_maxt_cur_extreme_stat;

// to avoid pathological multithreading issues, this is not a bitset
static unsigned char* g_perm_adapt_stop;

static uint32_t g_adapt_m_table[MODEL_BLOCKSIZE];
static uintptr_t* g_indiv_nonmale_include2;
static uintptr_t* g_indiv_male_include2 = NULL;
static uintptr_t* g_is_invalid;
static uint32_t* g_marker_uidxs;
static uint32_t g_model_fisher;
static uint32_t g_assoc_thread_ct;
static uintptr_t g_perm_vec_ct;
static uint32_t g_thread_block_ctl;
static uint32_t g_maxt_block_base;
static uint32_t g_block_start;
static uint32_t g_qblock_start;
static uint32_t g_block_diff;
static uint32_t g_perms_done;
static uint32_t g_first_adapt_check;
static uint32_t g_pheno_nm_ct;
static uint32_t g_case_ct;
static double g_adaptive_intercept;
static double g_adaptive_slope;
static double g_aperm_alpha;
static double g_adaptive_ci_zt;
static uint32_t g_is_x;
static uint32_t g_is_y;
static uint32_t g_is_haploid;
static int32_t g_is_model_prec;

static uint32_t g_tot_quotient;
static uint64_t g_totq_magic;
static uint32_t g_totq_preshift;
static uint32_t g_totq_postshift;
static uint32_t g_totq_incr;
static sfmt_t** g_sfmtp_arr;

static uint32_t g_cluster_ct;
static uint32_t* g_cluster_map;
static uint32_t* g_cluster_starts;
static uint32_t* g_cluster_case_cts;

// per-cluster magic number sets
static uintptr_t* g_cluster_cc_perm_preimage;
static uint32_t* g_tot_quotients;
static uint64_t* g_totq_magics;
static uint32_t* g_totq_preshifts;
static uint32_t* g_totq_postshifts;
static uint32_t* g_totq_incrs;

static uint32_t* g_indiv_to_cluster;
static uint32_t* g_qassoc_cluster_thread_wkspace;

THREAD_RET_TYPE model_assoc_gen_perms_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uint32_t case_ct = g_case_ct;
  uint32_t tot_quotient = g_tot_quotient;
  uint64_t totq_magic = g_totq_magic;
  uint32_t totq_preshift = g_totq_preshift;
  uint32_t totq_postshift = g_totq_postshift;
  uint32_t totq_incr = g_totq_incr;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  sfmt_t* __restrict__ sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_assoc_thread_ct;
  for (; pidx < pmax; pidx++) {
    generate_cc_perm_vec(pheno_nm_ct, case_ct, tot_quotient, totq_magic, totq_preshift, totq_postshift, totq_incr, &(perm_vecs[pidx * pheno_nm_ctl2]), sfmtp);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE model_assoc_gen_cluster_perms_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  sfmt_t* __restrict__ sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t cluster_ct = g_cluster_ct;
  uint32_t* cluster_map = g_cluster_map;
  uint32_t* cluster_starts = g_cluster_starts;
  uint32_t* cluster_case_cts = g_cluster_case_cts;
  uintptr_t* cluster_cc_perm_preimage = g_cluster_cc_perm_preimage;
  uint32_t* tot_quotients = g_tot_quotients;
  uint64_t* totq_magics = g_totq_magics;
  uint32_t* totq_preshifts = g_totq_preshifts;
  uint32_t* totq_postshifts = g_totq_postshifts;
  uint32_t* totq_incrs = g_totq_incrs;
  for (; pidx < pmax; pidx++) {
    generate_cc_cluster_perm_vec(pheno_nm_ct, cluster_cc_perm_preimage, cluster_ct, cluster_map, cluster_starts, cluster_case_cts, tot_quotients, totq_magics, totq_preshifts, totq_postshifts, totq_incrs, &(perm_vecs[pidx * pheno_nm_ctl2]), sfmtp);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE qassoc_gen_perms_thread(void* arg) {
  // Used by QT --assoc, --linear, and --logistic.
  //
  // Takes an array of phenotype values in g_pheno_d2 of length g_pheno_nm_ct,
  // and populates g_perm_vecstd[] with permutations of those values.  Also
  // requires g_sfmtp_arr[] and g_assoc_thread_ct to be initialized.
  //
  // g_perm_vecstd is individual-major.  The nth permutation is stored across
  //   g_perm_vecstd[n]
  //   g_perm_vecstd[n + perm_vec_ctcl8m]
  //   g_perm_vecstd[n + 2 * perm_vec_ctcl8m]
  //   ...
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ctcl8 = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  uintptr_t perm_vec_ctcl8m = perm_vec_ctcl8 * CACHELINE_DBL;
  double* pheno_d2 = g_pheno_d2;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint32_t pmin = CACHELINE_DBL * ((((uint64_t)tidx) * perm_vec_ctcl8) / g_assoc_thread_ct);
  uint32_t pmax = CACHELINE_DBL * ((((uint64_t)tidx + 1) * perm_vec_ctcl8) / g_assoc_thread_ct);
  double* perm_vecstd = &(g_perm_vecstd[pmin]);
  uint32_t pdiff = pmax - pmin;
  uint32_t poffset = 0;
  uint32_t indiv_idx = 1;
  uint32_t tot_quotient;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uint32_t urand;
  uint32_t uii;
  double* wptr;
  double* wptr2;
  double* wptr3;
  double cur_source;
  if (((uintptr_t)tidx) + 1 == g_assoc_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  cur_source = *pheno_d2++;
  wptr = perm_vecstd;
  for (; poffset < pdiff; poffset++) {
    *wptr++ = cur_source;
  }
  for (; indiv_idx < pheno_nm_ct; indiv_idx++) {
    tot_quotient = 4294967296LLU / (indiv_idx + 1);
    upper_bound = (indiv_idx + 1) * tot_quotient - 1;
    magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
    cur_source = *pheno_d2++;
    wptr = &(perm_vecstd[indiv_idx * perm_vec_ctcl8m]);
    wptr2 = perm_vecstd;
    for (poffset = 0; poffset < pdiff; poffset++) {
      do {
	urand = sfmt_genrand_uint32(sfmtp);
      } while (urand > upper_bound);
      uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
      wptr3 = &(wptr2[uii * perm_vec_ctcl8m]);
      *wptr++ = *wptr3;
      *wptr3 = cur_source;
      wptr2++;
    }
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE qassoc_gen_cluster_perms_thread(void* arg) {
  // Variant of qassoc_gen_perms_thread() which restricts permutations to be
  // within-cluster.
  // On top of the qassoc_gen_perms_thread requirements, this also needs
  // g_cluster_ct, g_cluster_map, g_cluster_starts,
  // g_qassoc_cluster_thread_wkspace, and g_indiv_to_cluster to be initialized.
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ctcl8 = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  uintptr_t perm_vec_ctcl8m = perm_vec_ctcl8 * CACHELINE_DBL;
  double* pheno_d2 = g_pheno_d2;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint32_t pmin = CACHELINE_DBL * ((((uint64_t)tidx) * perm_vec_ctcl8) / g_assoc_thread_ct);
  uint32_t pmax = CACHELINE_DBL * ((((uint64_t)tidx + 1) * perm_vec_ctcl8) / g_assoc_thread_ct);
  double* perm_vecstd = &(g_perm_vecstd[pmin]);
  uint32_t cluster_ct = g_cluster_ct;
  uint32_t cluster_ctcl = (cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  uint32_t* cluster_map = g_cluster_map;
  uint32_t* cluster_starts = g_cluster_starts;
  uint32_t* in_cluster_positions = &(g_qassoc_cluster_thread_wkspace[tidx * cluster_ctcl * CACHELINE_INT32]);
  uint32_t* indiv_to_cluster = g_indiv_to_cluster;
  uint32_t pdiff = pmax - pmin;
  uint32_t poffset = 0;
  uint32_t indiv_idx = 0;
  uint32_t cluster_idx;
  uint32_t cur_in_cluster_pos;
  uint32_t* cur_map_start;
  uint32_t tot_quotient;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uint32_t urand;
  uint32_t uii;
  double* wptr;
  double* wptr2;
  double* wptr3;
  double cur_source;
  if (((uintptr_t)tidx) + 1 == g_assoc_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  fill_uint_zero(in_cluster_positions, cluster_ct);
  for (; indiv_idx < pheno_nm_ct; indiv_idx++) {
    cur_source = *pheno_d2++;
    cluster_idx = indiv_to_cluster[indiv_idx];
    if (cluster_idx == 0xffffffffU) {
      cur_in_cluster_pos = 0;
    } else {
      cur_in_cluster_pos = in_cluster_positions[cluster_idx];
      in_cluster_positions[cluster_idx] += 1;
    }
    wptr = &(perm_vecstd[indiv_idx * perm_vec_ctcl8m]);
    if (!cur_in_cluster_pos) {
      for (poffset = 0; poffset < pdiff; poffset++) {
        *wptr++ = cur_source;
      }
    } else {
      cur_map_start = &(cluster_map[cluster_starts[cluster_idx]]);
      tot_quotient = 4294967296LLU / (cur_in_cluster_pos + 1);
      upper_bound = (cur_in_cluster_pos + 1) * tot_quotient - 1;
      magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
      wptr2 = perm_vecstd;
      for (poffset = 0; poffset < pdiff; poffset++) {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	wptr3 = &(wptr2[cur_map_start[uii] * perm_vec_ctcl8m]);
	*wptr++ = *wptr3;
	*wptr3 = cur_source;
	wptr2++;
      }
    }
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE assoc_adapt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t is_x = g_is_x;
  uint32_t is_y = g_is_y;
  uint32_t is_haploid = g_is_haploid;
  uint32_t min_ploidy = 2;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ nonmale_vec = g_indiv_nonmale_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uint32_t* __restrict__ gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
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
  if (is_haploid) { // includes g_is_x
    min_ploidy = 1;
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    // guaranteed during loading that g_perm_adapt_stop[] is not set yet
    marker_idx = adapt_m_table[marker_bidx];
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
    if (model_fisher) {
      stat_high = orig_1mpval[marker_idx] + EPSILON;
      stat_low = orig_1mpval[marker_idx] - EPSILON;
    } else {
      if (orig_1mpval[marker_idx] == -9) {
	perm_adapt_stop[marker_idx] = 1;
	perm_attempt_ct[marker_idx] = next_adapt_check;
	perm_2success_ct[marker_idx] = next_adapt_check;
	continue;
      }
      stat_high = orig_chisq[marker_idx] + EPSILON;
      stat_low = orig_chisq[marker_idx] - EPSILON;
    }
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!is_haploid) {
	vec_set_freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
      } else if (is_x) {
        vec_set_freq_x(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_set_ct, &case_missing_ct);
      } else if (!is_y) {
	vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &uii, &case_set_ct);
	case_missing_ct += uii;
      } else {
	vec_set_freq_y(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), nonmale_vec, &case_set_ct, &case_missing_ct);
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
	  dxx = 1.0 - fisher22(case_set_ct, uii - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - uii);
	} else {
	  dxx = chi22_eval(case_set_ct, uii, col1_sum, tot_obs);
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
	  pval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
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
  THREAD_RETURN;
}

THREAD_RET_TYPE assoc_maxt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uint32_t is_x = g_is_x;
  uint32_t is_x_or_y = is_x || g_is_y;
  uint32_t is_haploid = g_is_haploid;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t block_start = g_block_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + block_start;
  uint32_t marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct128 * 288]);
#else
  uint32_t perm_ct64 = (perm_vec_ct + 63) / 64;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct64 * 144]);
#endif
  uint32_t* git_homrar_cts = NULL;
  uint32_t* git_missing_cts = NULL;
  uint32_t* git_het_cts = NULL;
  uintptr_t perm_vec_ctcl4m = (perm_vec_ct + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1));
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t min_ploidy = 2;
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  uintptr_t* loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ nonmale_vec = g_indiv_nonmale_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double* msa_ptr = NULL;
  uint16_t* ldrefs = g_ldrefs;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
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
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (is_haploid) { // includes g_is_x
    min_ploidy = 1;
  }
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (model_fisher) {
      gpd = &(precomp_d[3 * precomp_width * marker_bidx]);
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
    } else {
      if (g_orig_1mpval[marker_idx] == -9) {
	if (msa_ptr) {
	  for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	    *msa_ptr++ = -9;
	  }
	}
	perm_2success_ct[marker_idx++] += perm_vec_ct;
	continue;
      }
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      stat_high = orig_chisq[marker_idx] + EPSILON;
      stat_low = orig_chisq[marker_idx] - EPSILON;
    }
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
    if (!is_x_or_y) {
      ldref = ldrefs[marker_idx];
      if (!is_haploid) {
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
#ifdef __LP64__
        fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * (perm_vec_ctcl4m / 2));
#else
        fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * perm_vec_ctcl4m);
#endif
        calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
#ifdef __LP64__
        fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct128 * 72);
#else
        fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct64 * 72);
#endif
      } else {
	memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
	calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
      }
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!is_x_or_y) {
	if (!is_haploid) {
	  case_missing_ct = git_missing_cts[pidx];
	  case_set_ct = row1x_sum - (git_het_cts[pidx] + 2 * (case_missing_ct + git_homrar_cts[pidx]));
	} else {
	  case_missing_ct = git_missing_cts[pidx] + git_het_cts[pidx];
	  case_set_ct = row1x_sum - case_missing_ct - git_homrar_cts[pidx];
	}
      } else {
	if (is_x) {
	  vec_set_freq_x(pheno_nm_ctl2, loadbuf_cur, &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_set_ct, &case_missing_ct);
	} else {
	  vec_set_freq_y(pheno_nm_ctl2, loadbuf_cur, &(perm_vecs[pidx * pheno_nm_ctl2]), nonmale_vec, &case_set_ct, &case_missing_ct);
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
	    sval = fisher22_tail_pval(ukk, ujj - ukk, col1_sum - ukk, col2_sum + ukk - ujj, gpui[6 * uii + 5] - 1, gpd[3 * uii], gpd[3 * uii + 1], gpd[3 * uii + 2], case_set_ct);
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
	  sval = fisher22(case_set_ct, uii - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - uii);
	  if (sval < stat_low) {
	    success_2incr += 2;
	  } else if (sval < stat_high) {
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
  THREAD_RETURN;
}

THREAD_RET_TYPE qassoc_adapt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t marker_bidx = g_qblock_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uint32_t first_adapt_check = g_first_adapt_check;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* git_qt_g_prod = &(g_thread_git_qbufs[perm_vec_ctcl8m * tidx * 3]);
  double* git_qt_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 3 + 1)]);
  double* git_qt_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 3 + 2)]);
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ orig_chiabs = g_orig_chisq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  uintptr_t next_cqg;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uintptr_t ulii;
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
  double beta;
  double betasq;
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
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    next_adapt_check = first_adapt_check;
    missing_ct = missing_cts[marker_idx];
    nanal = pheno_nm_ct - missing_ct;
    homcom_ct = homcom_cts[marker_idx];
    het_ct = het_cts[marker_idx];
    if ((nanal < 3) || (homcom_ct == nanal) || (het_ct == nanal)) {
      // the current code might otherwise report a spurious association if
      // geno_var is zero, so we explicitly check for it here.
      // yes, that could also be caused by all-homrars if the user actively
      // tries to shoot themselves in the foot.  that's not my problem.
      perm_adapt_stop[marker_idx] = 1;
      perm_attempt_ct[marker_idx] = 0;
      continue;
    }
    homrar_ct = nanal - het_ct - homcom_ct;
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
	next_cqg = (next_cqg + CACHELINE_DBL - 1) & (~(CACHELINE_DBL - 1));
	if (next_cqg > perm_vec_ct) {
	  next_cqg = perm_vec_ct;
	}
	calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, next_cqg - pidx, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecstd[pidx]), &(git_qt_g_prod[pidx]));
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
	  sval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	  dxx = adaptive_ci_zt * sqrt(sval * (1 - sval) / ((int32_t)next_adapt_check));
	  dyy = sval - dxx; // lower bound
	  dzz = sval + dxx; // upper bound
	  if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	    perm_adapt_stop[marker_idx] = 1;
	    perm_attempt_ct[marker_idx] = next_adapt_check;
	    fill_double_zero(git_qt_g_prod, next_cqg);
	    fill_double_zero(git_qt_sum, next_cqg);
	    fill_double_zero(git_qt_ssq, next_cqg);
	    goto qassoc_adapt_thread_lesszero;
	  }
	}
	next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
      }
    }
    fill_double_zero(git_qt_g_prod, perm_vec_ctcl8m * 3);
  qassoc_adapt_thread_lesszero:
    perm_2success_ct[marker_idx] += success_2incr;
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE qassoc_adapt_lin_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t marker_bidx = g_qblock_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uint32_t first_adapt_check = g_first_adapt_check;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* git_qt_het_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * tidx * 6]);
  double* git_qt_het_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 1)]);
  double* git_qt_homrar_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 2)]);
  double* git_qt_homrar_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 3)]);
  double* git_qt_missing_sum = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 4)]);
  double* git_qt_missing_ssq = &(g_thread_git_qbufs[perm_vec_ctcl8m * (tidx * 6 + 5)]);
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ orig_linsq = g_orig_linsq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  uintptr_t next_cqg;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uintptr_t ulii;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  intptr_t geno_sum;
  uint32_t nanal;
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
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
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
	next_cqg = (next_cqg + CACHELINE_DBL - 1) & (~(CACHELINE_DBL - 1));
	if (next_cqg > perm_vec_ct) {
	  next_cqg = perm_vec_ct;
	}
	calc_qgit_lin(pheno_nm_ct, perm_vec_ctcl8m, next_cqg - pidx, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecstd[pidx]), &(git_qt_het_sum[pidx]));
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
	  sval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	  dxx = adaptive_ci_zt * sqrt(sval * (1 - sval) / ((int32_t)next_adapt_check));
	  dyy = sval - dxx;
	  dzz = sval + dxx;
	  if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	    perm_adapt_stop[marker_idx] = 1;
	    perm_attempt_ct[marker_idx] = next_adapt_check;
	    fill_double_zero(git_qt_het_sum, next_cqg);
	    fill_double_zero(git_qt_het_ssq, next_cqg);
	    fill_double_zero(git_qt_homrar_sum, next_cqg);
	    fill_double_zero(git_qt_homrar_ssq, next_cqg);
	    fill_double_zero(git_qt_missing_sum, next_cqg);
	    fill_double_zero(git_qt_missing_ssq, next_cqg);
	    goto qassoc_adapt_lin_thread_lesszero;
	  }
	}
	next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
      }
    }
    fill_double_zero(git_qt_het_sum, perm_vec_ctcl8m * 6);
  qassoc_adapt_lin_thread_lesszero:
    perm_2success_ct[marker_idx] += success_2incr;
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE qassoc_maxt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t qblock_start = g_qblock_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + qblock_start;
  uint32_t marker_bidx_start = qblock_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  double* qresultbuf = g_qresultbuf;
  uintptr_t* loadbuf = g_loadbuf;
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  uint16_t* ldrefs = g_ldrefs;
  double* __restrict__ orig_chiabs = g_orig_chisq;
  double* msa_ptr = NULL;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  double* git_qt_g_prod;
  double* git_qt_sum;
  double* git_qt_ssq;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
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
  uint32_t marker_idx_tmp;
  int32_t missing_ct_tmp;
  int32_t het_ct_tmp;
  int32_t homcom_ct_tmp;
  int32_t homrar_ct_tmp;
  uint32_t loop_ceil;
  uintptr_t cur_cost;
  uint32_t ldref;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
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
	      cur_cost = qrem_cost2(pheno_nm_ctl2, &(loadbuf[marker_bidx2 * pheno_nm_ctl2]), loadbuf_cur);
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
      fill_double_zero(git_qt_g_prod, perm_vec_ctcl8m * 3);
      calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, perm_vec_ct, loadbuf_cur, perm_vecstd, git_qt_g_prod);
    } else {
      memcpy(git_qt_g_prod, &(qresultbuf[3 * ldref * perm_vec_ctcl8m]), 3 * perm_vec_ctcl8m * sizeof(double));
      calc_qrem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecstd, git_qt_g_prod);
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
  THREAD_RETURN;
}

THREAD_RET_TYPE qassoc_maxt_lin_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t qblock_start = g_qblock_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + qblock_start;
  uint32_t marker_bidx_start = qblock_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  double* qresultbuf = g_qresultbuf;
  uintptr_t* loadbuf = g_loadbuf;
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* msa_ptr = NULL;
  uint16_t* ldrefs = g_ldrefs;
  double* __restrict__ orig_linsq = g_orig_linsq;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  double* git_qt_het_sum;
  double* git_qt_het_ssq;
  double* git_qt_homrar_sum;
  double* git_qt_homrar_ssq;
  double* git_qt_missing_sum;
  double* git_qt_missing_ssq;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  intptr_t geno_sum;
  uint32_t nanal;
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
  uint32_t success_2incr;
  double stat_high;
  double stat_low;
  double sval;
  uint32_t ldref;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
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
      fill_double_zero(git_qt_het_sum, perm_vec_ctcl8m * 6);
      calc_qgit_lin(pheno_nm_ct, perm_vec_ctcl8m, perm_vec_ct, loadbuf_cur, perm_vecstd, git_qt_het_sum);
    } else {
      memcpy(git_qt_het_sum, &(qresultbuf[6 * ldref * perm_vec_ctcl8m]), 6 * perm_vec_ctcl8m * sizeof(double));
      calc_qrem_lin(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecstd, git_qt_het_sum);
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_adapt_domrec_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uint32_t* __restrict__ gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
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
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    if (model_fisher) {
      if (orig_1mpval[marker_idx] == -9) {
	perm_adapt_stop[marker_idx] = 1;
	perm_attempt_ct[marker_idx] = 0;
	continue;
      }
      stat_high = orig_1mpval[marker_idx] + EPSILON;
      stat_low = orig_1mpval[marker_idx] - EPSILON;
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
      vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &uii, &case_homx_ct);
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
	  dxx = 1.0 - fisher22(case_homx_ct, uii - case_homx_ct, col1_sum - case_homx_ct, col2_sum + case_homx_ct - uii);
	} else {
	  dxx = chi22_eval(case_homx_ct, uii, col1_sum, tot_obs);
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
	  pval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_domrec_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t block_start = g_block_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + block_start;
  uint32_t marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct128 * 288]);
#else
  uint32_t perm_ct64 = (perm_vec_ct + 63) / 64;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct64 * 144]);
#endif
  uint32_t* git_homrar_cts = NULL;
  uint32_t* git_missing_cts = NULL;
  uint32_t* git_het_cts = NULL;
  uintptr_t perm_vec_ctcl4m = (perm_vec_ct + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1));
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double* msa_ptr = NULL;
  uint16_t* ldrefs = g_ldrefs;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
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
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (model_fisher) {
      if (orig_1mpval[marker_idx] == -9) {
      model_maxt_domrec_thread_skip_marker:
	marker_idx++;
	if (msa_ptr) {
	  for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	    *msa_ptr++ = -9;
	  }
	}
	continue;
      }
      gpd = &(precomp_d[3 * precomp_width * marker_bidx]);
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
    } else {
      if (orig_chisq[marker_idx] == -9) {
	goto model_maxt_domrec_thread_skip_marker;
      }
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      stat_high = orig_chisq[marker_idx] + EPSILON;
      stat_low = orig_chisq[marker_idx] - EPSILON;
    }
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
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
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * (perm_vec_ctcl4m / 2));
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * perm_vec_ctcl4m);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct128 * 72);
#else
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct64 * 72);
#endif
    } else {
      memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
      calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
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
	    sval = fisher22_tail_pval(ukk, ujj - ukk, col1_sum - ukk, col2_sum + ukk - ujj, gpui[6 * uii + 5] - 1, gpd[3 * uii], gpd[3 * uii + 1], gpd[3 * uii + 2], case_homx_ct);
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
	  sval = fisher22(case_homx_ct, uii - case_homx_ct, col1_sum - case_homx_ct, col2_sum + case_homx_ct - uii);
	  if (sval < stat_low) {
	    success_2incr += 2;
	  } else if (sval < stat_high) {
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_adapt_trend_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uint32_t* __restrict__ gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
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
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    next_adapt_check = first_adapt_check;
    if (orig_1mpval[marker_idx] == -9) {
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
      vec_set_freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_com_ct, &case_missing_ct);
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
	  pval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_trend_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t block_start = g_block_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + block_start;
  uint32_t marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct128 * 288]);
#else
  uint32_t perm_ct64 = (perm_vec_ct + 63) / 64;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct64 * 144]);
#endif
  uint32_t* git_homrar_cts = NULL;
  uint32_t* git_missing_cts = NULL;
  uint32_t* git_het_cts = NULL;
  uintptr_t perm_vec_ctcl4m = (perm_vec_ct + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1));
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double* msa_ptr = NULL;
  uint16_t* ldrefs = g_ldrefs;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t* loadbuf_cur;
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
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (orig_1mpval[marker_idx] == -9) {
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
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
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * (perm_vec_ctcl4m / 2));
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * perm_vec_ctcl4m);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct128 * 72);
#else
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct64 * 72);
#endif
    } else {
      memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
      calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_adapt_gen_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t marker_idx;
  uintptr_t pidx;
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
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    if (model_fisher) {
      if (orig_1mpval[marker_idx] == -9) {
	perm_adapt_stop[marker_idx] = 1;
	perm_attempt_ct[marker_idx] = 0;
	continue;
      }
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
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
      vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &case_het_ct, &case_homcom_ct);
      if (model_fisher) {
        uii = case_ct - case_het_ct - case_homcom_ct - case_missing_ct;
	// this is very slow.  a precomputed 2-dimensional table could improve
	// matters, but I doubt it's worth the effort for now.
	dxx = fisher23(case_homcom_ct, case_het_ct, uii, homcom_ct - case_homcom_ct, het_ct - case_het_ct, homrar_ct - uii);
	if (dxx < stat_low) {
	  success_2incr += 2;
	} else if (dxx < stat_high) {
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
	  pval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_gen_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t block_start = g_block_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + block_start;
  uint32_t marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct128 * 288]);
#else
  uint32_t perm_ct64 = (perm_vec_ct + 63) / 64;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct64 * 144]);
#endif
  uint32_t* git_homrar_cts = NULL;
  uint32_t* git_missing_cts = NULL;
  uint32_t* git_het_cts = NULL;
  uintptr_t perm_vec_ctcl4m = (perm_vec_ct + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1));
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double* msa_ptr = NULL;
  uint16_t* ldrefs = g_ldrefs;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  uint32_t missing_col;
  intptr_t tot_obs;
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
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (model_fisher) {
      if (orig_1mpval[marker_idx] == -9) {
	marker_idx++;
	continue;
      }
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
    } else {
      if (orig_chisq[marker_idx] == -9) {
	marker_idx++;
	continue;
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
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
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * (perm_vec_ctcl4m / 2));
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * perm_vec_ctcl4m);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, loadbuf_cur, perm_vecst, git_homrar_cts, thread_git_wkspace);
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct128 * 72);
#else
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct64 * 72);
#endif
    } else {
      memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
      calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      case_missing_ct = git_missing_cts[pidx];
      case_het_ct = git_het_cts[pidx];
      case_homcom_ct = case_ct - case_missing_ct - case_het_ct - git_homrar_cts[pidx];
      if (model_fisher) {
        uii = case_ct - case_het_ct - case_homcom_ct - case_missing_ct;
	sval = fisher23(case_homcom_ct, case_het_ct, uii, homcom_ct - case_homcom_ct, het_ct - case_het_ct, homrar_ct - uii);
	if (sval < stat_low) {
	  success_2incr += 2;
	} else if (sval < stat_high) {
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_adapt_best_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* is_invalid = g_is_invalid;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  uint32_t* __restrict__ adapt_m_table = g_adapt_m_table;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uint32_t* __restrict__ gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
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
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    if (model_fisher) {
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 + EPSILON - orig_1mpval[marker_idx];
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
      vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &case_het_ct, &case_homcom_ct);
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
	  dxx = fisher22(case_com_ct, 2 * uii - case_com_ct, com_ct - case_com_ct, 2 * ukk + case_com_ct - com_ct);
	  if (dxx < stat_low) {
	    goto model_adapt_best_thread_betterstat;
	  } else if (dxx < stat_high) {
	    ujj = 1;
	  }
	  if (!skip_domrec) {
	    dxx = fisher22(case_homcom_ct, uii - case_homcom_ct, homcom_ct - case_homcom_ct, ukk + case_homcom_ct - homcom_ct);
	    if (dxx < stat_low) {
	      goto model_adapt_best_thread_betterstat;
	    } else if (dxx < stat_high) {
	      ujj = 1;
	    }
	    dxx = fisher22(case_homrar_ct, uii - case_homrar_ct, homrar_ct - case_homrar_ct, ukk + case_homrar_ct - homrar_ct);
	    if (dxx < stat_low) {
	      goto model_adapt_best_thread_betterstat;
	    } else if (dxx < stat_high) {
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
	  pval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
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
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_best_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t block_start = g_block_start;
  uint32_t maxt_block_base = g_maxt_block_base;
  uint32_t maxt_block_base2 = maxt_block_base + block_start;
  uint32_t marker_bidx_start = block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t maxt_block_base3 = maxt_block_base + marker_bidx_start;
  uint32_t marker_bidx = marker_bidx_start;
  uintptr_t marker_idx = maxt_block_base3;
  uint32_t marker_bceil = block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct128 * 288]);
#else
  uint32_t perm_ct64 = (perm_vec_ct + 63) / 64;
  uint32_t* thread_git_wkspace = &(g_thread_git_wkspace[tidx * perm_ct64 * 144]);
#endif
  uint32_t* git_homrar_cts = NULL;
  uint32_t* git_missing_cts = NULL;
  uint32_t* git_het_cts = NULL;
  uintptr_t perm_vec_ctcl4m = (perm_vec_ct + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1));
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  uint32_t* resultbuf = g_resultbuf;
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* is_invalid = g_is_invalid;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  double* msa_ptr = NULL;
  uint16_t* ldrefs = g_ldrefs;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t* loadbuf_cur;
  uintptr_t pidx;
  int32_t missing_ct;
  intptr_t tot_obs;
  intptr_t com_ct;
  intptr_t rar_ct;
  intptr_t het_ct;
  intptr_t homrar_ct;
  intptr_t homcom_ct;
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
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (g_mperm_save_all) {
    msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (model_fisher) {
      gpd = &(precomp_d[9 * precomp_width * marker_bidx]);
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
      default_best_stat = 1;
    } else {
      if (orig_chisq[marker_idx] == -9) {
	marker_idx++;
	continue;
      }
      gpd = &(precomp_d[6 * precomp_width * marker_bidx]);
      stat_high = orig_chisq[marker_idx] + EPSILON;
      stat_low = orig_chisq[marker_idx] - EPSILON;
      default_best_stat = 0;
    }
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
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
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
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * (perm_vec_ctcl4m / 2));
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 3 * perm_vec_ctcl4m);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct128 * 72);
#else
      fill_ulong_zero((uintptr_t*)thread_git_wkspace, perm_ct64 * 72);
#endif
    } else {
      memcpy(git_homrar_cts, &(resultbuf[3 * ldref * perm_vec_ctcl4m]), 3 * perm_vec_ctcl4m * sizeof(int32_t));
      calc_rem(pheno_nm_ct, perm_vec_ct, loadbuf_cur, &(loadbuf[ldref * pheno_nm_ctl2]), perm_vecst, git_homrar_cts, thread_git_wkspace);
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
	    best_stat = fisher22_tail_pval(ukk, ujj - ukk, com_ct - ukk, rar_ct + ukk - ujj, gpui[18 * uii + 5] - 1, gpd[9 * uii], gpd[9 * uii + 1], gpd[9 * uii + 2], case_com_ct);
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
	      sval = fisher22_tail_pval(ukk, ujj - ukk, homcom_ct - ukk, homrar_ct + het_ct + ukk - ujj, gpui[18 * uii + 11] - 1, gpd[9 * uii + 3], gpd[9 * uii + 4], gpd[9 * uii + 5], case_homcom_ct);
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
	      sval = fisher22_tail_pval(ukk, ujj - ukk, homrar_ct - ukk, homcom_ct + het_ct + ukk - ujj, gpui[18 * uii + 17] - 1, gpd[9 * uii + 6], gpd[9 * uii + 7], gpd[9 * uii + 8], case_homrar_ct);
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
	  best_stat = fisher22(case_com_ct, 2 * uii - case_com_ct, com_ct - case_com_ct, 2 * ukk + case_com_ct - com_ct);
	  if (!skip_domrec) {
	    sval = fisher22(case_homcom_ct, uii - case_homcom_ct, homcom_ct - case_homcom_ct, ukk + case_homcom_ct - homcom_ct);
	    if (sval < best_stat) {
	      best_stat = sval;
	    }
	    sval = fisher22(case_homrar_ct, uii - case_homrar_ct, homrar_ct - case_homrar_ct, ukk + case_homrar_ct - homrar_ct);
	    if (sval < best_stat) {
	      best_stat = sval;
	    }
	  }
	  if (best_stat < stat_low) {
	    cur_add = 2;
	  } else if (best_stat < stat_high) {
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
  THREAD_RETURN;
}

void get_model_assoc_precomp_bounds(uint32_t missing_ct, uint32_t is_model, uint32_t* minp, uint32_t* ctp) {
  // Estimate which case missing counts are most common.
  // Expected value = (g_case_ct * missing_ct / g_pheno_nm_ct)
  // If X-chromosome and (!is_model):
  //   Lower bound = max(0, missing_ct - 2 * (g_pheno_nm_ct - g_case_ct))
  //   Upper bound = min(g_case_ct * 2, missing_ct)
  //   (Could be a bit more precise if we tracked missing male and female
  //    counts separately, but whatever)
  //   Each male automatically contributes 1 to initial missing_ct!
  // Otherwise:
  //   Lower bound = max(0, missing_ct - (g_pheno_nm_ct - g_case_ct))
  //   Upper bound = min(g_case_ct, missing_ct)
  double xval = ((double)(g_case_ct * ((int64_t)missing_ct))) / ((double)((intptr_t)g_pheno_nm_ct));
  intptr_t lbound = (intptr_t)(xval + EPSILON + 1 - ((double)((intptr_t)g_precomp_width)) * 0.5);
  intptr_t ctrl_ct = g_pheno_nm_ct - g_case_ct;
  intptr_t ubound = missing_ct;
  intptr_t lii;
  if (lbound < 0) {
    lbound = 0;
  }
  if (g_is_x && (!is_model)) {
    lii = missing_ct - (2 * ctrl_ct);
    if (((uintptr_t)ubound) > g_case_ct * 2) {
      ubound = g_case_ct * 2;
    }
  } else {
    lii = missing_ct - ctrl_ct;
    if (((uintptr_t)ubound) > g_case_ct) {
      ubound = g_case_ct;
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

int32_t model_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  int32_t retval = 0;
  FILE* outfile = NULL;
  FILE* outfile_msa = NULL;
  uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  uint32_t model_assoc = model_modifier & MODEL_ASSOC;
  uint32_t model_adapt = model_modifier & MODEL_PERM;
  uint32_t model_maxt = model_modifier & MODEL_MPERM;
  uint32_t model_perms = model_adapt | model_maxt;
  uint32_t model_trendonly = model_modifier & MODEL_TRENDONLY;
  uint32_t model_perm_best = !(model_modifier & MODEL_PMASK);
  uint32_t model_perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t assoc_counts = model_modifier & MODEL_ASSOC_COUNTS;
  uint32_t assoc_p2 = model_modifier & MODEL_ASSOC_P2;
  uint32_t display_ci = (ci_size > 0);
  uint32_t perms_total = 0;
  int32_t x_code = -1;
  int32_t y_code = -1;
  uint32_t male_ct = 0;
  uint32_t nonmale_ct = 0;
  uint32_t ctrl_male_ct = 0;
  uint32_t case_male_ct = 0;
  uint32_t ctrl_nonmale_ct = 0;
  uint32_t case_nonmale_ct = 0;
  uint32_t load_ctrl_ct = 0;
  uint32_t load_case_ct = 0;
  uintptr_t* indiv_nonmale_ctrl_include2 = NULL;
  uintptr_t* indiv_nonmale_case_include2 = NULL;
  uintptr_t* indiv_male_ctrl_include2 = NULL;
  uintptr_t* indiv_male_case_include2 = NULL;
  uintptr_t* cur_ctrl_include2 = NULL;
  uintptr_t* cur_case_include2 = NULL;
  double* orig_odds = NULL;
  double dxx = 0.0;
  double dww = 0.0;
  double dvv = 0.0;
  double mult_p = 0.0;
  double gen_p = 0.0;
  double dom_p = 0.0;
  double rec_p = 0.0;
  double ca_chisq = 0.0;
  uint32_t pct = 0;
  uint32_t perm_pass_idx = 0;
  uintptr_t perm_vec_ctcl4m = 0;
  uint32_t model_fisherx = (model_modifier & MODEL_FISHER) && (!(model_modifier & MODEL_PTREND));
  char* chrom_name_ptr = NULL;
  uint32_t chrom_name_len = 0;
  char chrom_name_buf[4];
  uint32_t mu_table[MODEL_BLOCKSIZE];
  uint32_t uibuf[4];
  char wbuf[48];
  char* wptr_start;
  char* wptr;
  char* wptr2;
  char* wptr_mid;
  char* wptr_mid2;
  char* outname_end2;
  uint32_t fill_orig_chisq;
  uint32_t marker_unstopped_ct;
  uint32_t gender_req;
  uint32_t ctrl_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t marker_bidx;
  uintptr_t marker_uidx; // loading
  uintptr_t marker_uidx2; // writing
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t indiv_uidx;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uint32_t* missp;
  uint32_t* setp;
  uint32_t* hetp;
  double* o1mpptr;
  double* ooptr;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* indiv_ctrl_include2;
  uintptr_t* indiv_case_include2;
  uint32_t load_indiv_ct;
  uintptr_t ulii;
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
    logprint("Warning: Skipping --assoc/--model since less than two phenotypes are present.\n");
    goto model_assoc_ret_1;
  }
  g_orig_chisq = NULL;
  g_model_fisher = model_modifier & MODEL_FISHER;
  g_pheno_nm_ct = pheno_nm_ct;
  g_perms_done = 0;
  g_aperm_alpha = aperm_alpha;
  g_is_model_prec = (model_modifier & MODEL_PREC)? 1 : 0;
  g_is_y = 0;
  g_mperm_save_all = NULL;

  if (assoc_p2) {
    logprint("Error: --assoc p2 not yet supported.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
  }
  if (model_maxt) {
    perms_total = model_mperm_val;
    if (wkspace_alloc_d_checked(&g_maxt_extreme_stat, sizeof(double) * perms_total)) {
      goto model_assoc_ret_NOMEM;
    }
    if (model_fisherx) {
      for (uii = 0; uii < perms_total; uii++) {
	g_maxt_extreme_stat[uii] = 1;
      }
    } else {
      fill_double_zero(g_maxt_extreme_stat, perms_total);
    }
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(&outfile_msa, outname, "w")) {
	goto model_assoc_ret_OPEN_FAIL;
      }
      sprintf(logbuf, "Dumping all permutation %svalues to %s.\n", model_fisherx? "p-" : "chi-square ", outname);
      logprintb();
    }
  } else {
    mperm_save = 0;
    if (model_adapt) {
      perms_total = aperm_max;
    }
  }
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto model_assoc_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctl2 - 1] = 0;
  if (model_assoc) {
    if (g_model_fisher) {
      outname_end2 = memcpyb(outname_end, ".assoc.fisher", 14);
    } else {
      outname_end2 = memcpyb(outname_end, ".assoc", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "Writing C/C --assoc report to %s...", outname);
    logprintb();
    fflush(stdout);
    sprintf(tbuf, " CHR %%%us         BP   A1 ", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
    if (assoc_counts) {
      fputs("     C_A      C_U   A2 ", outfile);
    } else {
      fputs("     F_A      F_U   A2 ", outfile);
    }
    if (!g_model_fisher) {
      fputs("       CHISQ ", outfile);
    }
    if (fputs_checked("           P           OR ", outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (display_ci) {
      uii = (uint32_t)((int32_t)(ci_size * 100));
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
    uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0);
    if (uii) {
      sprintf(logbuf, "Excluding %u haploid marker%s from --model analysis.\n", uii, (uii == 1)? "" : "s");
      logprintb();
    }
    marker_ct -= uii;
    if (!marker_ct) {
      logprint("Error: No markers remaining for --model analysis.\n");
      goto model_assoc_ret_INVALID_CMDLINE;
    }
    outname_end2 = memcpyb(outname_end, ".model", 7);
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "Writing --model report to %s...", outname);
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
    sprintf(tbuf, " CHR %%%us   A1   A2     TEST            AFF          UNAFF ", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
    if (!g_model_fisher) {
      fputs("       CHISQ   DF ", outfile);
    } else {
      outname_end2 = memcpyb(outname_end2, ".fisher", 8);
    }
    if (fputs_checked("           P\n", outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
  }
  marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  g_adaptive_ci_zt = ltqnorm(1 - aperm_beta / (2.0 * marker_ct));
  if (wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * pheno_nm_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&g_orig_1mpval, marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&g_missing_cts, marker_ct * sizeof(uint32_t))) {
    goto model_assoc_ret_NOMEM;
  }
  if (model_assoc) {
    if (wkspace_alloc_d_checked(&orig_odds, marker_ct * sizeof(double)) ||
        wkspace_alloc_ui_checked(&g_set_cts, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
  }
  if ((!model_assoc) || model_maxt) {
    if (wkspace_alloc_ui_checked(&g_het_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_homcom_cts, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
  }
  x_code = chrom_info_ptr->x_code;
  y_code = chrom_info_ptr->y_code;
  gender_req = ((x_code != -1) && is_set(chrom_info_ptr->chrom_mask, x_code)) || (model_assoc && (((y_code != -1) && is_set(chrom_info_ptr->chrom_mask, y_code))));
  if (gender_req) {
    if (wkspace_alloc_ul_checked(&g_indiv_nonmale_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&g_indiv_male_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
      goto model_assoc_ret_NOMEM;
    }
    fill_ulong_zero(g_indiv_male_include2, pheno_nm_ctl2);
    indiv_uidx = 0;
    indiv_idx = 0;
    do {
      indiv_uidx = next_set_ul_unsafe(pheno_nm, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
      do {
        if (IS_SET(sex_male, indiv_uidx)) {
	  SET_BIT_DBL(g_indiv_male_include2, indiv_idx);
	  male_ct++;
	}
	indiv_idx++;
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < pheno_nm_ct);
    vec_init_invert(pheno_nm_ct, g_indiv_nonmale_include2, g_indiv_male_include2);
    nonmale_ct = pheno_nm_ct - male_ct;
  }
  fill_orig_chisq = (!model_fisherx) || (mtest_adjust && (!g_model_fisher));
  if (fill_orig_chisq) {
    if (wkspace_alloc_d_checked(&g_orig_chisq, marker_ct * sizeof(double))) {
      goto model_assoc_ret_NOMEM;
    }
    fill_double_zero(g_orig_chisq, marker_ct);
  }

  if (model_perms) {
    if (cluster_starts) {
      retval = cluster_include_and_reindex(unfiltered_indiv_ct, pheno_nm, 1, pheno_c, pheno_nm_ct, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, &g_cluster_case_cts, &g_cluster_cc_perm_preimage);
      if (retval) {
	goto model_assoc_ret_1;
      }
      if (!g_cluster_ct) {
        logprint("Error: No size 2+ clusters for permutation test.\n");
	goto model_assoc_ret_INVALID_CMDLINE;
      }
      retval = cluster_alloc_and_populate_magic_nums(g_cluster_ct, g_cluster_map, g_cluster_starts, &g_tot_quotients, &g_totq_magics, &g_totq_preshifts, &g_totq_postshifts, &g_totq_incrs);
      if (retval) {
        goto model_assoc_ret_1;
      }
    }
    g_ldrefs = (uint16_t*)wkspace_alloc(marker_ct * sizeof(uint16_t));
    if (!g_ldrefs) {
      goto model_assoc_ret_NOMEM;
    }
#ifdef __LP64__
    fill_ulong_one((uintptr_t*)g_ldrefs, (marker_ct + 3) / 4);
#else
    fill_ulong_one((uintptr_t*)g_ldrefs, (marker_ct + 1) / 2);
#endif
    g_sfmtp_arr = (sfmt_t**)wkspace_alloc(g_thread_ct * sizeof(intptr_t));
    if (!g_sfmtp_arr) {
      goto model_assoc_ret_NOMEM;
    }
    g_sfmtp_arr[0] = &sfmt;
    if (g_thread_ct > 1) {
      // make the permutation analysis reproducible with --seed + --threads
      for (uii = 1; uii < g_thread_ct; uii++) {
	g_sfmtp_arr[uii] = (sfmt_t*)wkspace_alloc(sizeof(sfmt_t));
	if (!g_sfmtp_arr[uii]) {
	  goto model_assoc_ret_NOMEM;
	}
	for (ujj = 0; ujj < 4; ujj++) {
	  uibuf[ujj] = sfmt_genrand_uint32(&sfmt);
	}
	sfmt_init_by_array(g_sfmtp_arr[uii], uibuf, 4);
      }
    }

    if (!(mperm_save & MPERM_DUMP_ALL)) {
      g_precomp_width = (1 + (int32_t)(sqrt(pheno_nm_ct) * EXPECTED_MISSING_FREQ * 5.65686));
    } else {
      g_precomp_width = 0;
    }
    if (wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
    if (model_maxt) {
      if (model_fisherx) {
	if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC))) {
	  if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 6 * MODEL_BLOCKSIZE * sizeof(uint32_t)) ||
	      wkspace_alloc_d_checked(&g_precomp_d, g_precomp_width * 3 * MODEL_BLOCKSIZE * sizeof(double))) {
	    goto model_assoc_ret_NOMEM;
	  }
	} else if (model_perm_best) {
	  if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 18 * MODEL_BLOCKSIZE * sizeof(uint32_t)) ||
	      wkspace_alloc_d_checked(&g_precomp_d, g_precomp_width * 9 * MODEL_BLOCKSIZE * sizeof(double))) {
	    goto model_assoc_ret_NOMEM;
	  }
	}
      } else if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
	if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 6 * MODEL_BLOCKSIZE * sizeof(uint32_t)) ||
	    wkspace_alloc_d_checked(&g_precomp_d, g_precomp_width * 2 * MODEL_BLOCKSIZE * sizeof(double))) {
	  goto model_assoc_ret_NOMEM;
	}
      } else if (model_perm_best) {
	if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 18 * MODEL_BLOCKSIZE * sizeof(uint32_t)) ||
	    wkspace_alloc_d_checked(&g_precomp_d, g_precomp_width * 6 * MODEL_BLOCKSIZE * sizeof(double))) {
	  goto model_assoc_ret_NOMEM;
	}
      }
    } else if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
      if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 4 * MODEL_BLOCKSIZE * sizeof(uint32_t))) {
	goto model_assoc_ret_NOMEM;
      }
    } else if (model_perm_best) {
      if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 12 * MODEL_BLOCKSIZE * sizeof(uint32_t))) {
	goto model_assoc_ret_NOMEM;
      }
    }
    if (model_perm_best) {
      if (wkspace_alloc_ul_checked(&g_is_invalid, marker_ctl * sizeof(intptr_t))) {
	goto model_assoc_ret_NOMEM;
      }
      fill_ulong_zero(g_is_invalid, marker_ctl);
    }
    fill_uint_zero(g_perm_2success_ct, marker_ct);
    g_thread_block_ctl = (MODEL_BLOCKSIZE + (g_thread_ct - 1)) / g_thread_ct;
    if (model_adapt) {
      if (wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_ct * sizeof(uint32_t)) ||
          wkspace_alloc_uc_checked(&g_perm_adapt_stop, marker_ct)) {
	goto model_assoc_ret_NOMEM;
      }
      for (uii = 0; uii < marker_ct; uii++) {
	g_perm_attempt_ct[uii] = aperm_max;
      }
      fill_ulong_zero((uintptr_t*)g_perm_adapt_stop, (marker_ct + sizeof(uintptr_t) - 1) / sizeof(uintptr_t));
    }
    if (!cluster_starts) {
      g_tot_quotient = 4294967296LLU / pheno_nm_ct;
      magic_num(g_tot_quotient, &g_totq_magic, &g_totq_preshift, &g_totq_postshift, &g_totq_incr);
    }
  }
  if (wkspace_alloc_ui_checked(&g_marker_uidxs, marker_ct * sizeof(uint32_t))) {
    goto model_assoc_ret_NOMEM;
  }
  g_case_ct = 0;
  if (wkspace_alloc_ul_checked(&indiv_ctrl_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_case_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
    goto model_assoc_ret_NOMEM;
  }
  fill_ulong_zero(indiv_case_include2, pheno_nm_ctl2);
  indiv_uidx = 0;
  indiv_idx = 0;
  do {
    indiv_uidx = next_set_ul_unsafe(pheno_nm, indiv_uidx);
    indiv_uidx_stop = next_unset_ul(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    do {
      if (IS_SET(pheno_c, indiv_uidx)) {
	SET_BIT_DBL(indiv_case_include2, indiv_idx);
	g_case_ct++;
      }
      indiv_idx++;
    } while (++indiv_uidx < indiv_uidx_stop);
  } while (indiv_idx < pheno_nm_ct);
  vec_init_invert(pheno_nm_ct, indiv_ctrl_include2, indiv_case_include2);
  ctrl_ct = pheno_nm_ct - g_case_ct;
  if (gender_req) {
    // todo: get rid of these and just use the functions called by the
    // permutation tests
    if (wkspace_alloc_ul_checked(&indiv_nonmale_ctrl_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&indiv_nonmale_case_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&indiv_male_ctrl_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&indiv_male_case_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
      goto model_assoc_ret_NOMEM;
    }
    fill_ulong_zero(indiv_male_case_include2, pheno_nm_ctl2);
    indiv_uidx = 0;
    indiv_idx = 0;
    do {
      indiv_uidx = next_set_ul_unsafe(pheno_nm, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
      do {
	if (IS_SET(sex_male, indiv_uidx) & IS_SET(pheno_c, indiv_uidx)) {
	  SET_BIT_DBL(indiv_male_case_include2, indiv_idx);
	  case_male_ct++;
	}
	indiv_idx++;
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < pheno_nm_ct);
    vec_init_andnot(pheno_nm_ctl2, indiv_male_ctrl_include2, g_indiv_male_include2, indiv_male_case_include2);
    vec_init_andnot(pheno_nm_ctl2, indiv_nonmale_case_include2, indiv_case_include2, indiv_male_case_include2);
    vec_init_andnot(pheno_nm_ctl2, indiv_nonmale_ctrl_include2, indiv_ctrl_include2, indiv_male_ctrl_include2);
    ctrl_male_ct = male_ct - case_male_ct;
    case_nonmale_ct = g_case_ct - case_male_ct;
    ctrl_nonmale_ct = ctrl_ct - ctrl_male_ct;
  }

  for (uii = 1; uii <= MODEL_BLOCKSIZE; uii++) {
    g_loadbuf[uii * pheno_nm_ctl2 - 2] = 0;
    g_loadbuf[uii * pheno_nm_ctl2 - 1] = 0;
  }
  if (model_perms) {
    if (wkspace_left < pheno_nm_ctl2 * sizeof(intptr_t)) {
      goto model_assoc_ret_NOMEM;
    }
  }
  marker_unstopped_ct = marker_ct;

  // ----- begin main loop -----
 model_assoc_more_perms:
  if (model_perms) {
    if (!perm_pass_idx) {
      fputs(" [generating permutations]", stdout);
      fflush(stdout);
    }
    if (model_adapt) {
      if (perm_pass_idx) {
	while (g_first_adapt_check <= g_perms_done) {
	  // APERM_MAX prevents infinite loop here
	  g_first_adapt_check += (int32_t)(aperm_init_interval + ((int32_t)g_first_adapt_check) * aperm_interval_slope);
	}
      } else {
	if (aperm_min < aperm_init_interval) {
	  g_first_adapt_check = (int32_t)aperm_init_interval;
	} else {
	  g_first_adapt_check = aperm_min;
	}
	g_adaptive_intercept = aperm_init_interval;
	g_adaptive_slope = aperm_interval_slope;
      }
      g_perm_vec_ct = wkspace_left / (pheno_nm_ctl2 * sizeof(intptr_t));
    } else {
      // g_perm_vec_ct memory allocation dependencies:
      //   g_maxt_thread_results: (8 * g_perm_vec_ct, cacheline-aligned) *
      //     g_thread_ct
      //   g_perm_vecst: 16 * ((g_perm_vec_ct + 127) / 128) * pheno_nm_ct
      //   g_thread_git_wkspace: ((perm_vec_ct + 127) / 128) * 1152 * thread_ct
      //   g_resultbuf: MODEL_BLOCKSIZE * (4 * perm_vec_ct, CL-aligned) * 3
      //   g_perm_vecs: pheno_nm_ctl2 * sizeof(intptr_t) * g_perm_vec_ct
      //   g_mperm_save_all (if needed): marker_ct * 8 * g_perm_vec_ct
      // If we force g_perm_vec_ct to be a multiple of 128, then we have
      //   g_perm_vec_ct * (17 * g_thread_ct + 12 * MODEL_BLOCKSIZE +
      //                    pheno_nm_ct + sizeof(intptr_t) * pheno_nm_ctl2
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
        g_perm_vec_ct = 128 * (wkspace_left / (128LL * sizeof(intptr_t) * pheno_nm_ctl2 + 2176LL * g_thread_ct + 1536LL + 16LL * pheno_nm_ct + 128LL * sizeof(double) * marker_ct));
      } else {
        g_perm_vec_ct = 128 * (wkspace_left / (128LL * sizeof(intptr_t) * pheno_nm_ctl2 + 2176LL * g_thread_ct + 1536LL + 16LL * pheno_nm_ct));
      }
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    } else if (!g_perm_vec_ct) {
      goto model_assoc_ret_NOMEM;
    }
    perm_vec_ctcl4m = (g_perm_vec_ct + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1));
    g_perms_done += g_perm_vec_ct;
    g_perm_vecs = (uintptr_t*)wkspace_alloc(g_perm_vec_ct * pheno_nm_ctl2 * sizeof(intptr_t));
    if (g_perm_vec_ct > g_thread_ct) {
      g_assoc_thread_ct = g_thread_ct;
    } else {
      g_assoc_thread_ct = g_perm_vec_ct;
    }
    if (!cluster_starts) {
      if (spawn_threads(threads, &model_assoc_gen_perms_thread, g_assoc_thread_ct)) {
	goto model_assoc_ret_THREAD_CREATE_FAIL;
      }
      ulii = 0;
      model_assoc_gen_perms_thread((void*)ulii);
      join_threads(threads, g_assoc_thread_ct);
    } else {
      if (spawn_threads(threads, &model_assoc_gen_cluster_perms_thread, g_assoc_thread_ct)) {
	goto model_assoc_ret_THREAD_CREATE_FAIL;
      }
      ulii = 0;
      model_assoc_gen_cluster_perms_thread((void*)ulii);
      join_threads(threads, g_assoc_thread_ct);
    }
    if (!model_adapt) {
      ulii = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
      g_maxt_thread_results = (double*)wkspace_alloc(g_thread_ct * ulii * CACHELINE);
      g_resultbuf = (uint32_t*)wkspace_alloc(perm_vec_ctcl4m * 3 * MODEL_BLOCKSIZE * sizeof(int32_t));
#ifdef __LP64__
      ulii = ((g_perm_vec_ct + 127) / 128) * 16;
      g_perm_vecst = (uint32_t*)wkspace_alloc(ulii * pheno_nm_ct);
#else
      ulii = ((g_perm_vec_ct + 31) / 32) * 4;
      g_perm_vecst = (uint32_t*)wkspace_alloc(ulii * pheno_nm_ct);
      ulii = ((g_perm_vec_ct + 63) / 64) * 8;
#endif
      g_thread_git_wkspace = (uint32_t*)wkspace_alloc(ulii * 72 * g_thread_ct);
      transpose_perms(g_perm_vecs, g_perm_vec_ct, pheno_nm_ct, g_perm_vecst);
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)g_thread_git_wkspace, ulii * 9 * g_thread_ct);
#else
      fill_ulong_zero((uintptr_t*)g_thread_git_wkspace, ulii * 18 * g_thread_ct);
#endif
      if (mperm_save & MPERM_DUMP_ALL) {
	g_mperm_save_all = (double*)wkspace_alloc(marker_ct * g_perm_vec_ct * sizeof(double));
      }
    }
    if (!perm_pass_idx) {
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                         \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", stdout);
    }
  }
  if (!perm_pass_idx) {
    fputs(" 0%", stdout);
    fflush(stdout);
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
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
	refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &g_is_haploid);
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if (g_is_haploid && (!g_is_x)) {
	  if (g_is_y) {
	    cur_ctrl_include2 = indiv_male_ctrl_include2;
	    cur_case_include2 = indiv_male_case_include2;
	    load_indiv_ct = male_ct;
	    load_case_ct = case_male_ct;
	  } else {
	    cur_ctrl_include2 = indiv_ctrl_include2;
	    cur_case_include2 = indiv_case_include2;
	    load_indiv_ct = pheno_nm_ct;
	    load_case_ct = g_case_ct;
	  }
	  load_ctrl_ct = load_indiv_ct - load_case_ct;
	}
      } else {
	while (1) {
	  do {
	    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	  } while (marker_uidx >= chrom_end);
	  uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  g_is_x = (uii == (uint32_t)x_code);
	  if ((!IS_SET(haploid_mask, uii)) || g_is_x) {
	    break;
	  }
	  marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
	}
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
      }
      chrom_name_ptr = chrom_name_buf;
      chrom_name_len = 4;
      if (uii <= chrom_info_ptr->max_code) {
	memset(chrom_name_buf, 32, 2);
        intprint2(&(chrom_name_buf[2]), uii);
      } else if (zero_extra_chroms) {
	memcpy(chrom_name_buf, "   0", 4);
      } else {
	ujj = strlen(chrom_info_ptr->nonstd_names[uii]);
	if (ujj < 4) {
	  fw_strcpyn(4, ujj, chrom_info_ptr->nonstd_names[uii], chrom_name_buf);
	} else {
	  chrom_name_ptr = chrom_info_ptr->nonstd_names[uii];
	  chrom_name_len = ujj;
	}
      }
    } else if (model_maxt) {
      marker_idx -= MODEL_BLOCKKEEP;
      if (marker_idx) { // max(T) initial block special case, see below
        memcpy(g_loadbuf, &(g_loadbuf[(MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * pheno_nm_ctl2]), MODEL_BLOCKKEEP * pheno_nm_ctl2 * sizeof(intptr_t));
        memcpy(g_resultbuf, &(g_resultbuf[3 * (MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * perm_vec_ctcl4m]), MODEL_BLOCKKEEP * perm_vec_ctcl4m * 3 * sizeof(int32_t));
      }
      g_block_start = MODEL_BLOCKKEEP;
    } else {
      g_block_start = 0;
    }
    block_size = g_block_start;
    block_end = marker_unstopped_ct - marker_idx;
    if ((!marker_idx) && model_maxt) {
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
      if (model_adapt && g_perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(g_loadbuf[block_size * pheno_nm_ctl2]);
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_ptr, pheno_nm_ct, pheno_nm, IS_SET(marker_reverse, marker_uidx))) {
	goto model_assoc_ret_READ_FAIL;
      }
      if (model_adapt) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      if ((!model_assoc) && g_is_x) {
	force_missing((unsigned char*)(&(g_loadbuf[block_size * pheno_nm_ctl2])), g_indiv_male_include2, pheno_nm_ct);
      }
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (block_size == g_block_start) {
      continue;
    }
    if (!perm_pass_idx) {
      // basic --assoc/--model
      o1mpptr = &(g_orig_1mpval[marker_idx + g_block_start]);
      missp = &(g_missing_cts[marker_idx + g_block_start]);
      if (model_assoc) {
	setp = &(g_set_cts[marker_idx + g_block_start]);
	ooptr = &(orig_odds[marker_idx + g_block_start]);
	for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	  marker_uidx2 = mu_table[marker_bidx];
	  g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	  if (!g_is_haploid) {
	    if (model_maxt) {
	      single_marker_cc_3freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &unn, &uoo, &ujj, &upp, &uqq, &umm);
	      g_het_cts[marker_idx + marker_bidx] = uoo + uqq;
	      g_homcom_cts[marker_idx + marker_bidx] = unn + upp;
	      uii = 2 * unn + uoo;
	      ukk = 2 * upp + uqq;
	    } else {
	      single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &uii, &ujj, &ukk, &umm);
	    }
	    *missp = ujj + umm;
	    *setp = uii + ukk;
	    ujj = 2 * (ctrl_ct - ujj) - uii;
	    umm = 2 * (g_case_ct - umm) - ukk;
	  } else if (g_is_x) {
	    single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_nonmale_ctrl_include2, indiv_nonmale_case_include2, &uii, &ujj, &ukk, &umm);
	    *missp = 2 * (ujj + umm);
	    *setp = uii + ukk;
	    ujj = 2 * (ctrl_nonmale_ct - ujj) - uii;
	    umm = 2 * (case_nonmale_ct - umm) - ukk;
	    haploid_single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_male_ctrl_include2, indiv_male_case_include2, &unn, &uoo, &upp, &uqq);
	    *missp += uoo + uqq + male_ct;
	    *setp += unn + upp;
	    uoo = ctrl_male_ct - uoo - unn;
	    uqq = case_male_ct - uqq - upp;
	    uii += unn;
	    ujj += uoo;
	    ukk += upp;
	    umm += uqq;
	  } else if (g_is_haploid) {
	    haploid_single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm);
	    *missp = ujj + umm;
	    *setp = uii + ukk;
	    ujj = load_ctrl_ct - ujj - uii;
	    umm = load_case_ct - umm - ukk;
	    if (g_is_y) {
	      *missp += nonmale_ct;
	    } else if (model_maxt) {
	      g_het_cts[marker_idx + marker_bidx] = 0;
	      g_homcom_cts[marker_idx + marker_bidx] = *setp;
	    }
	  }
	  da1 = umm;
	  da2 = ukk;
	  du1 = ujj;
	  du2 = uii;
	  if (g_model_fisher) {
	    pval = fisher22(uii, ujj, ukk, umm);
	    *o1mpptr = 1 - pval;
	  } else if (!assoc_p2) {
	    if ((!umm) && (!ujj)) {
	      *o1mpptr = -9;
	      pval = -1;
	      dxx = 0;
	      if (fill_orig_chisq) {
		g_orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      dxx = chi22_eval(ukk, ukk + umm, uii + ukk, uii + ujj + ukk + umm);
	      pval = chiprob_p(dxx, 1);
	      *o1mpptr = 1 - pval;
	      if (fill_orig_chisq) {
		g_orig_chisq[marker_idx + marker_bidx] = dxx;
	      }
	    }
	  } else {
	    // --p2 not written yet
	  }
	  *ooptr = (da1 * du2) / (du1 * da2);
	  if (pval <= pfilter) {
	    a1ptr = &(marker_alleles[(2 * marker_uidx2) * max_marker_allele_len]);
	    a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1) * max_marker_allele_len]);
	    wptr = memcpya(tbuf, chrom_name_ptr, chrom_name_len);
	    *wptr++ = ' ';
	    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	    *wptr++ = ' ';
	    wptr = uint32_writew10(wptr, marker_pos[marker_uidx2]);
	    if (max_marker_allele_len == 1) {
	      memset(wptr, 32, 4);
	      wptr[4] = *a1ptr;
	      wptr = &(wptr[5]);
	    } else {
	      *wptr = ' ';
	      wptr = fw_strcpy(4, a1ptr, &(wptr[1]));
	    }
	    *wptr++ = ' ';
	    if (umm + ukk) {
	      if (assoc_counts) {
		wptr = uint32_writew8(wptr, umm);
	      } else {
		wptr = double_g_writewx4(wptr, da1 / (da1 + da2), 8);
	      }
	      *wptr++ = ' ';
	    } else {
	      wptr = memcpya(wptr, "      NA ", 9);
	    }
	    if (ujj + uii) {
	      if (assoc_counts) {
		wptr = uint32_writew8(wptr, ujj);
	      } else {
		wptr = double_g_writewx4(wptr, du1 / (du1 + du2), 8);
	      }
	    } else {
	      wptr = memcpya(wptr, "      NA", 8);
	    }
	    if (max_marker_allele_len == 1) {
	      memset(wptr, 32, 4);
	      wptr[4] = *a2ptr;
	      wptr = &(wptr[5]);
	    } else {
	      *wptr = ' ';
	      wptr = fw_strcpy(4, a2ptr, &(wptr[1]));
	    }
	    *wptr++ = ' ';
	    if (g_model_fisher) {
	      wptr = double_g_writewx4(wptr, pval, 12);
	    } else {
	      if (pval > -1) {
		wptr = double_g_writewx4(double_g_writewx4x(wptr, dxx, 12, ' '), pval, 12);
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
	      wptr = double_g_writewx4(wptr, *ooptr, 12);
	      if (display_ci) {
		dxx = log(*ooptr);
		dyy = sqrt(1 / da1 + 1 / da2 + 1 / du1 + 1 / du2);
		dzz = ci_zt * dyy;
		dww = exp(dxx - dzz);
		dvv = exp(dxx + dzz);
		*wptr++ = ' ';
		wptr = double_g_writewx4(double_g_writewx4x(double_g_writewx4x(wptr, dyy, 12, ' '), dww, 12, ' '), dvv, 12);
	      }
	    }
	    wptr = memcpya(wptr, " \n", 2);
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	  }
	  missp++;
	  setp++;
	  o1mpptr++;
	  ooptr++;
	}
      } else {
	// repurpose setp as homcom_cts pointer
	setp = &(g_homcom_cts[marker_idx + g_block_start]);
	hetp = &(g_het_cts[marker_idx + g_block_start]);
	for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	  marker_uidx2 = mu_table[marker_bidx];
	  g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	  single_marker_cc_3freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &uii, &ujj, &ukk, &umm, &unn, &uoo);
	  *missp = ukk + uoo;
	  *setp = uii + umm;
	  ukk = pheno_nm_ct - g_case_ct - uii - ujj - ukk;
	  uoo = g_case_ct - umm - unn - uoo;
	  *hetp = ujj + unn;
	  is_invalid = (uoo < model_cell_ct) || (unn < model_cell_ct) || (umm < model_cell_ct) || (ukk < model_cell_ct) || (ujj < model_cell_ct) || (uii < model_cell_ct);
	  a1ptr = &(marker_alleles[(2 * marker_uidx2) * max_marker_allele_len]);
	  a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1) * max_marker_allele_len]);
	  wptr = memcpya(tbuf, chrom_name_ptr, chrom_name_len);
	  *wptr++ = ' ';
	  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	  if (max_marker_allele_len == 1) {
	    memset(wptr, 32, 4);
	    wptr[4] = *a1ptr;
	    memset(&(wptr[5]), 32, 4);
	    wptr[9] = *a2ptr;
	    wptr = &(wptr[10]);
	  } else {
	    *wptr++ = ' ';
	    wptr = fw_strcpy(4, a1ptr, wptr);
	    *wptr++ = ' ';
	    wptr = fw_strcpy(4, a2ptr, wptr);
	  }
	  memset(wptr, 32, 2);
	  wptr = &(wptr[2]);
	  wptr_mid = wptr;
	  if (!model_trendonly) {
	    memcpy(wptr, "   GENO ", 8);
	    wptr2 = uint32_write(uint32_writex(uint32_writex(wbuf, uoo, '/'), unn, '/'), umm);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr[8]));
	    *wptr++ = ' ';
	    wptr2 = uint32_write(uint32_writex(uint32_writex(wbuf, ukk, '/'), ujj, '/'), uii);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	    *wptr++ = ' ';
	    if (is_invalid) {
	      gen_p = -9;
	      if (fill_orig_chisq && (model_modifier & MODEL_PGEN)) {
		g_orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      if (g_model_fisher) {
		gen_p = fisher23(uii, ujj, ukk, umm, unn, uoo);
	      } else {
		chi23_evalx(uii, ujj, ukk, umm, unn, uoo, &dvv, &upp);
		gen_p = chiprob_px(dvv, upp);
		if (fill_orig_chisq && (model_modifier & MODEL_PGEN)) {
		  if (dvv != -9) {
		    g_orig_chisq[marker_idx + marker_bidx] = dvv;
		  } else {
		    g_orig_chisq[marker_idx + marker_bidx] = 0;
		  }
		}
	      }
	    }
	    if (gen_p < -1) {
	      wptr = model_assoc_tna(g_model_fisher, wptr);
	    } else {
	      if (!g_model_fisher) {
		wptr = memcpya(double_g_writewx4(wptr, dvv, 12), "    ", 4);
		*wptr++ = '0' + upp;
		*wptr++ = ' ';
	      }
	      wptr = double_g_writewx4x(wptr, gen_p, 12, '\n');
	    }
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	  }
	  memcpy(wptr_mid, "  TREND ", 8);
	  wptr2 = uint32_write(uint32_writex(wbuf, uoo * 2 + unn, '/'), umm * 2 + unn);
	  wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr_mid[8]));
	  *wptr++ = ' ';
	  wptr2 = uint32_write(uint32_writex(wbuf, ukk * 2 + ujj, '/'), uii * 2 + ujj);
	  wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	  *wptr++ = ' ';
	  wptr_mid2 = wptr; // save this for next line
	  ca_chisq = ca_trend_evalx(umm * 2 + unn, umm + unn + uoo, ujj + unn, uii + umm, uii + ujj + ukk + umm + unn + uoo);
	  ca_p = chiprob_px(ca_chisq, 1);
	  if (fill_orig_chisq && (model_modifier & MODEL_PTREND)) {
	    if (ca_chisq != -9) {
	      g_orig_chisq[marker_idx + marker_bidx] = ca_chisq;
	    } else {
	      g_orig_chisq[marker_idx + marker_bidx] = 0;
	    }
	  }
	  if (ca_p > -1) {
	    if (!g_model_fisher) {
	      wptr = memcpya(double_g_writewx4(wptr, ca_chisq, 12), "    1 ", 6);
	    }
	    wptr = double_g_writewx4x(wptr, ca_p, 12, '\n');
	  } else {
	    wptr = model_assoc_tna(g_model_fisher, wptr);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	  if (!model_trendonly) {
	    memcpy(wptr_mid, "ALLELIC", 7);
	    wptr = wptr_mid2;
	    if (g_model_fisher) {
	      mult_p = fisher22(2 * uoo + unn, 2 * umm + unn, 2 * ukk + ujj, 2 * uii + ujj);
	    } else {
	      dww = chi22_evalx(2 * uoo + unn, 2 * (uoo + unn + umm), 2 * (uoo + ukk) + unn + ujj, 2 * (uoo + unn + umm + ukk + ujj + uii));
	      mult_p = chiprob_px(dww, 1);
	    }
	    if (mult_p > -1) {
	      if (!g_model_fisher) {
		wptr = memcpya(double_g_writewx4(wptr, dww, 12), "    1 ", 6);
	      }
	      wptr = double_g_writewx4x(wptr, mult_p, 12, '\n');
	    } else {
	      wptr = model_assoc_tna(g_model_fisher, wptr);
	    }
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    memcpy(wptr_mid, "    DOM", 7);
	    wptr2 = uint32_write(uint32_writex(wbuf, uoo + unn, '/'), umm);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr_mid[8]));
	    *wptr++ = ' ';
	    wptr2 = uint32_write(uint32_writex(wbuf, ukk + ujj, '/'), uii);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	    *wptr++ = ' ';
	    if (is_invalid) {
	      dom_p = -9;
	      if (fill_orig_chisq && (model_modifier & MODEL_PDOM)) {
		g_orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      if (g_model_fisher) {
		dom_p = fisher22(uoo + unn, umm, ukk + ujj, uii);
	      } else {
		dww = chi22_evalx(uoo + unn, uoo + unn + umm, uoo + unn + ukk + ujj, uoo + unn + umm + ukk + ujj + uii);
		dom_p = chiprob_px(dww, 1);
		if (fill_orig_chisq && (model_modifier & MODEL_PDOM)) {
		  if (dww != -9) {
		    g_orig_chisq[marker_idx + marker_bidx] = dww;
		  } else {
		    g_orig_chisq[marker_idx + marker_bidx] = 0;
		  }
		}
	      }
	    }
	    if (dom_p < -1) {
	      wptr = model_assoc_tna(g_model_fisher, wptr);
	    } else {
	      if (!g_model_fisher) {
		wptr = memcpya(double_g_writewx4(wptr, dww, 12), "    1 ", 6);
	      }
	      wptr = double_g_writewx4x(wptr, dom_p, 12, '\n');
	    }
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    memcpy(&(wptr_mid[4]), "REC", 3);
	    wptr2 = uint32_write(uint32_writex(wbuf, uoo, '/'), unn + umm);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, &(wptr_mid[8]));
	    *wptr++ = ' ';
	    wptr2 = uint32_write(uint32_writex(wbuf, ukk, '/'), ujj + uii);
	    wptr = fw_strcpyn(14, wptr2 - wbuf, wbuf, wptr);
	    *wptr++ = ' ';
	    if (is_invalid) {
	      rec_p = -9;
	      if (fill_orig_chisq && (model_modifier & MODEL_PREC)) {
		g_orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    } else {
	      if (g_model_fisher) {
		rec_p = fisher22(uoo, unn + umm, ukk, ujj + uii);
	      } else {
		dww = chi22_evalx(uoo, uoo + unn + umm, uoo + ukk, uoo + unn + umm + ukk + ujj + uii);
		rec_p = chiprob_px(dww, 1);
		if (fill_orig_chisq && (model_modifier & MODEL_PREC)) {
		  if (dww != -9) {
		    g_orig_chisq[marker_idx + marker_bidx] = dww;
		  } else {
		    g_orig_chisq[marker_idx + marker_bidx] = 0;
		  }
		}
	      }
	    }
	    if (rec_p < -1) {
	      wptr = model_assoc_tna(g_model_fisher, wptr);
	    } else {
	      if (!g_model_fisher) {
		wptr = memcpya(double_g_writewx4(wptr, dww, 12), "    1 ", 6);
	      }
	      wptr = double_g_writewx4x(wptr, rec_p, 12, '\n');
	    }
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
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
	      set_bit_ul(g_is_invalid, marker_idx + marker_bidx);
	    }
	    if (fill_orig_chisq) {
	      if (dxx != -9) {
		g_orig_chisq[marker_idx + marker_bidx] = inverse_chiprob(dxx, 1);
	      } else {
		g_orig_chisq[marker_idx + marker_bidx] = -9;
	      }
	    }
	    if (dxx != -9) {
	      dxx = 1 - dxx;
	    }
	  } else if (model_modifier & MODEL_PGEN) {
	    dxx = (gen_p >= 0)? (1 - gen_p) : -9;
	  } else if (model_modifier & MODEL_PDOM) {
	    dxx = (dom_p >= 0)? (1 - dom_p) : -9;
	  } else if (model_modifier & MODEL_PREC) {
	    dxx = (rec_p >= 0)? (1 - rec_p) : -9;
	  } else if (model_modifier & MODEL_PTREND) {
	    dxx = (ca_p >= 0)? ca_chisq : -9;
	  }
	  missp++;
	  setp++;
	  hetp++;
	  *o1mpptr++ = dxx;
	}
      }
    }
    if (model_perms) {
      g_block_diff = block_size - g_block_start;
      g_assoc_thread_ct = g_block_diff;
      if (g_assoc_thread_ct > g_thread_ct) {
	g_assoc_thread_ct = g_thread_ct;
      }
      if (model_maxt) {
	if (model_fisherx) {
	  g_maxt_cur_extreme_stat = g_maxt_extreme_stat[0];
	  for (uii = 1; uii < g_perm_vec_ct; uii++) {
	    dxx = g_maxt_extreme_stat[uii];
	    if (dxx > g_maxt_cur_extreme_stat) {
	      g_maxt_cur_extreme_stat = dxx;
	    }
	  }
	} else {
	  g_maxt_cur_extreme_stat = g_maxt_extreme_stat[0];
	  for (uii = 1; uii < g_perm_vec_ct; uii++) {
	    dxx = g_maxt_extreme_stat[uii];
	    if (dxx < g_maxt_cur_extreme_stat) {
	      g_maxt_cur_extreme_stat = dxx;
	    }
	  }
	}
      }
      if (model_assoc) {
	if (g_is_haploid) { // includes g_is_x
	  uqq = 1;
	} else {
	  uqq = 2;
	}
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = g_missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 0, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  uoo = g_set_cts[urr];
	  if (g_is_x) {
	    unn = 2 * g_case_ct;
	    upp = 2 * pheno_nm_ct - upp;
	  } else {
	    unn = uqq * g_case_ct;
	    upp = uqq * (pheno_nm_ct - upp);
	  }
	  ujj *= uqq;
	  ukk += uii * g_precomp_width;
	  if (g_model_fisher) {
	    dxx = 1 - g_orig_1mpval[urr];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		fisher22_precomp_pval_bounds(dxx, unn - ujj, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
		ujj += uqq;
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		fisher22_precomp_pval_bounds(dxx, unn - ujj, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
		fisher22_precomp_pval_bounds(g_maxt_cur_extreme_stat, unn - ujj, uoo, upp, uibuf, &(g_precomp_d[umm * 3]));
		g_precomp_ui[umm * 6 + 4] = uibuf[2];
		g_precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
		ujj += uqq;
	      }
	    }
	  } else {
	    dxx = g_orig_chisq[urr];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - ujj, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
		ujj += uqq;
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - ujj, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
		chi22_precomp_val_bounds(g_maxt_cur_extreme_stat, unn - ujj, uoo, upp, uibuf, &(g_precomp_d[umm * 2]));
		g_precomp_ui[umm * 6 + 4] = uibuf[2];
		g_precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
		ujj += uqq;
	      }
	    }
	  }
	}
      } else if (model_perm_best) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = g_missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 1, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  unn = 2 * g_case_ct;
	  uqq = 2 * (pheno_nm_ct - upp);
	  uoo = 2 * g_homcom_cts[urr] + g_het_cts[urr];
	  ukk += uii * g_precomp_width;
	  uss = 2 * ujj;
	  if (g_model_fisher) {
	    dxx = 1 - g_orig_1mpval[urr];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, unn - uss, uoo, uqq, &(g_precomp_ui[umm * 12]), NULL);
		uss += 2;
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, unn - uss, uoo, uqq, &(g_precomp_ui[umm * 18]), NULL);
	        fisher22_precomp_pval_bounds(g_maxt_cur_extreme_stat, 2 * g_case_ct - uss, uoo, uqq, uibuf, &(g_precomp_d[umm * 9]));
		g_precomp_ui[umm * 18 + 4] = uibuf[2];
		g_precomp_ui[umm * 18 + 5] = uibuf[3] - uibuf[2];
		uss += 2;
	      }
	    }
	    if (!IS_SET(g_is_invalid, urr)) {
	      upp = pheno_nm_ct - upp;
	      uoo = g_homcom_cts[urr];
	      uqq = upp - uoo - g_het_cts[urr];
	      ujj = g_case_ct - ujj;
	      if (model_adapt) {
		for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		  fisher22_precomp_pval_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 12 + 4]), NULL);
		  fisher22_precomp_pval_bounds(dxx, ujj, uqq, upp, &(g_precomp_ui[umm * 12 + 8]), NULL);
		  ujj--;
		}
	      } else {
		for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		  fisher22_precomp_pval_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 18 + 6]), NULL);
		  fisher22_precomp_pval_bounds(g_maxt_cur_extreme_stat, ujj, uoo, upp, uibuf, &(g_precomp_d[umm * 9 + 3]));
		  g_precomp_ui[umm * 18 + 10] = uibuf[2];
		  g_precomp_ui[umm * 18 + 11] = uibuf[3] - uibuf[2];
		  fisher22_precomp_pval_bounds(dxx, ujj, uqq, upp, &(g_precomp_ui[umm * 18 + 12]), NULL);
		  fisher22_precomp_pval_bounds(g_maxt_cur_extreme_stat, ujj, uqq, upp, uibuf, &(g_precomp_d[umm * 9 + 6]));
		  g_precomp_ui[umm * 18 + 16] = uibuf[2];
		  g_precomp_ui[umm * 18 + 17] = uibuf[3] - uibuf[2];
		  ujj--;
		}
	      }
	    }
	  } else {
	    dxx = g_orig_chisq[urr];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - uss, uoo, uqq, &(g_precomp_ui[umm * 12]), NULL);
		uss += 2;
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, unn - uss, uoo, uqq, &(g_precomp_ui[umm * 18]), NULL);
		chi22_precomp_val_bounds(g_maxt_cur_extreme_stat, unn - uss, uoo, uqq, uibuf, &(g_precomp_d[umm * 6]));
		g_precomp_ui[umm * 18 + 4] = uibuf[2];
		g_precomp_ui[umm * 18 + 5] = uibuf[3] - uibuf[2];
		uss += 2;
	      }
	    }
	    if (!IS_SET(g_is_invalid, urr)) {
	      upp = pheno_nm_ct - upp;
	      uoo = g_homcom_cts[urr];
	      uqq = upp - uoo - g_het_cts[urr];
	      ujj = g_case_ct - ujj;
	      if (model_adapt) {
		for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		  chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 12 + 4]), NULL);
		  chi22_precomp_val_bounds(dxx, ujj, uqq, upp, &(g_precomp_ui[umm * 12 + 8]), NULL);
		  ujj--;
		}
	      } else {
		for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		  chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 18 + 6]), NULL);
		  chi22_precomp_val_bounds(g_maxt_cur_extreme_stat, ujj, uoo, upp, uibuf, &(g_precomp_d[umm * 6 + 2]));
		  g_precomp_ui[umm * 18 + 10] = uibuf[2];
		  g_precomp_ui[umm * 18 + 11] = uibuf[3] - uibuf[2];
		  chi22_precomp_val_bounds(dxx, ujj, uqq, upp, &(g_precomp_ui[umm * 18 + 12]), NULL);
		  chi22_precomp_val_bounds(g_maxt_cur_extreme_stat, ujj, uqq, upp, uibuf, &(g_precomp_d[umm * 6 + 4]));
		  g_precomp_ui[umm * 18 + 16] = uibuf[2];
		  g_precomp_ui[umm * 18 + 17] = uibuf[3] - uibuf[2];
		  ujj--;
		}
	      }
	    }
	  }
	}
      } else if (model_modifier & MODEL_PTREND) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = g_missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 1, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  unn = g_het_cts[urr];
	  upp = pheno_nm_ct - upp; // tot_obs
	  uoo = g_homcom_cts[urr];
	  ukk += uii * g_precomp_width;
	  ujj = g_case_ct - ujj;
	  dxx = g_orig_chisq[urr];
	  if (model_adapt) {
	    for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	      ca_trend_precomp_val_bounds(dxx, ujj--, unn, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
	    }
	  } else {
	    for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	      ca_trend_precomp_val_bounds(dxx, ujj, unn, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
              ca_trend_precomp_val_bounds(g_maxt_cur_extreme_stat, ujj--, unn, uoo, upp, uibuf, &(g_precomp_d[umm * 2]));
	      g_precomp_ui[umm * 6 + 4] = uibuf[2];
	      g_precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
	    }
	  }
	}
      } else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  if (model_adapt) {
	    urr = g_adapt_m_table[uii];
	  } else {
	    urr = marker_idx + uii;
	  }
	  upp = g_missing_cts[urr];
	  get_model_assoc_precomp_bounds(upp, 1, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  upp = pheno_nm_ct - upp; // tot_obs
	  if (model_modifier & MODEL_PREC) {
	    uoo = upp - g_homcom_cts[urr] - g_het_cts[urr]; // col1_sum
	  } else {
	    uoo = g_homcom_cts[urr];
	  }
	  ukk += uii * g_precomp_width;
	  ujj = g_case_ct - ujj;
	  if (g_model_fisher) {
	    dxx = 1 - g_orig_1mpval[urr];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, ujj--, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
	        fisher22_precomp_pval_bounds(g_maxt_cur_extreme_stat, ujj--, uoo, upp, uibuf, &(g_precomp_d[umm * 3]));
		g_precomp_ui[umm * 6 + 4] = uibuf[2];
		g_precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
	      }
	    }
	  } else {
	    dxx = g_orig_chisq[urr];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, ujj--, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
		chi22_precomp_val_bounds(g_maxt_cur_extreme_stat, ujj--, uoo, upp, uibuf, &(g_precomp_d[umm * 2]));
		g_precomp_ui[umm * 6 + 4] = uibuf[2];
		g_precomp_ui[umm * 6 + 5] = uibuf[3] - uibuf[2];
	      }
	    }
	  }
	}
      }
      if (model_adapt) {
	ulii = 0;
	if (model_assoc) {
	  if (spawn_threads(threads, &assoc_adapt_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  assoc_adapt_thread((void*)ulii);
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads(threads, &model_adapt_domrec_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_domrec_thread((void*)ulii);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads(threads, &model_adapt_trend_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_trend_thread((void*)ulii);
	} else if (model_modifier & MODEL_PGEN) {
	  if (spawn_threads(threads, &model_adapt_gen_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_gen_thread((void*)ulii);
	} else {
	  if (spawn_threads(threads, &model_adapt_best_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_best_thread((void*)ulii);
	}
	join_threads(threads, g_assoc_thread_ct);
      } else {
	g_maxt_block_base = marker_idx;
	ulii = 0;
	if (model_assoc) {
	  if (spawn_threads(threads, &assoc_maxt_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  assoc_maxt_thread((void*)ulii);
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads(threads, &model_maxt_domrec_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_domrec_thread((void*)ulii);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads(threads, &model_maxt_trend_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_trend_thread((void*)ulii);
	} else if (model_modifier & MODEL_PGEN) {
	  if (spawn_threads(threads, &model_maxt_gen_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_gen_thread((void*)ulii);
	} else {
	  if (spawn_threads(threads, &model_maxt_best_thread, g_assoc_thread_ct)) {
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_best_thread((void*)ulii);
	}
	join_threads(threads, g_assoc_thread_ct);
	ulii = CACHELINE_DBL * ((g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL);
	if (model_fisherx) {
	  for (uii = 0; uii < g_assoc_thread_ct; uii++) {
	    ooptr = &(g_maxt_thread_results[uii * ulii]);
	    for (ujj = g_perms_done - g_perm_vec_ct; ujj < g_perms_done; ujj++) {
	      dxx = *ooptr++;
	      if (dxx < g_maxt_extreme_stat[ujj]) {
		g_maxt_extreme_stat[ujj] = dxx;
	      }
	    }
	  }
	} else {
	  for (uii = 0; uii < g_assoc_thread_ct; uii++) {
	    ooptr = &(g_maxt_thread_results[uii * ulii]);
	    for (ujj = g_perms_done - g_perm_vec_ct; ujj < g_perms_done; ujj++) {
	      dxx = *ooptr++;
	      if (dxx > g_maxt_extreme_stat[ujj]) {
		g_maxt_extreme_stat[ujj] = dxx;
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
	  putchar('\b');
	}
	pct = (marker_idx * 100LLU) / marker_unstopped_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = ((uint64_t)(pct + 1LLU) * marker_unstopped_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putchar('\b');
    }
    fputs("\b\b\b", stdout);
    logprint(" done.\n");
    if (model_perms) {
      wkspace_reset((unsigned char*)g_perm_vecs);
    }
    if (fclose_null(&outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (mtest_adjust) {
      retval = multcomp(outname, outname_end, g_marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, g_model_fisher? g_orig_1mpval : g_orig_chisq, pfilter, mtest_adjust, adjust_lambda, g_model_fisher, NULL);
      if (retval) {
	goto model_assoc_ret_1;
      }
    }
    if (mperm_save & MPERM_DUMP_ALL) {
      tbuf[0] = '0';
      wptr = &(tbuf[1]);
      a1ptr = &(tbuf[MAXLINELEN]);
      if (model_fisherx) {
	for (uii = 0; uii < marker_ct; uii++) {
	  *wptr++ = ' ';
	  dxx = g_orig_1mpval[uii];
	  if (dxx >= 0) {
	    wptr = double_g_write(wptr, 1 - dxx);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= a1ptr) {
	    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
      } else {
	for (uii = 0; uii < marker_ct; uii++) {
	  *wptr++ = ' ';
	  dxx = g_orig_chisq[uii];
	  if (dxx >= 0) {
	    wptr = double_g_write(wptr, dxx);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= a1ptr) {
	    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
      }
      *wptr++ = '\n';
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    }
  }
  if (model_perms) {
    if (mperm_save & MPERM_DUMP_ALL) {
      if (perm_pass_idx) {
	putchar(' ');
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = g_perm_vec_ct;
      ujj = 1 + g_perms_done - ulii;
      wptr = tbuf;
      a1ptr = &(tbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32_write(wptr, uii + ujj);
        o1mpptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = o1mpptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = double_g_write(wptr, dxx);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= a1ptr) {
	    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	goto model_assoc_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    wkspace_reset((unsigned char*)g_perm_vecs);
    if (g_perms_done < perms_total) {
      if (model_adapt) {
	marker_unstopped_ct = marker_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, 0, (marker_ct + sizeof(uintptr_t) - 1) / sizeof(uintptr_t));
	if (!marker_unstopped_ct) {
	  goto model_assoc_adapt_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done > 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto model_assoc_more_perms;
    }
    if (model_adapt) {
    model_assoc_adapt_perm_count:
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
    putchar('\r');
    sprintf(logbuf, "%u %s permutation%s complete.\n", g_perms_done, model_maxt? "max(T)" : "(adaptive)", (g_perms_done > 1)? "s" : "");
    logprintb();
    if (g_model_fisher && (model_modifier & MODEL_PTREND)) {
      outname_end2 -= 7; // remove ".fisher"
    }
    if (model_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	if (wkspace_alloc_c_checked(&a1ptr, FNAMESIZE)) {
	  goto model_assoc_ret_NOMEM;
	}
	ulii = outname_end - outname;
	memcpy(a1ptr, outname, ulii);
	memcpy(&(a1ptr[ulii]), ".mperm.dump.best", 17);
	sprintf(logbuf, "Dumping best permutation %svalues to %s.\n", model_fisherx? "p-" : "chi-square ", a1ptr);
	logprintb();
	if (fopen_checked(&outfile, a1ptr, "w")) {
	  goto model_assoc_ret_OPEN_FAIL;
	}
	dxx = 0;
	if (model_fisherx) {
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    if (g_orig_1mpval[marker_idx] > dxx) {
	      dxx = g_orig_1mpval[marker_idx];
	    }
	  }
	  dxx = 1 - dxx;
	} else {
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    if (g_orig_chisq[marker_idx] > dxx) {
	      dxx = g_orig_chisq[marker_idx];
	    }
	  }
	}
        memcpy(tbuf, "0 ", 2);
	wptr = double_g_writex(&(tbuf[2]), dxx, '\n');
	if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
	for (uii = 0; uii < perms_total; uii++) {
	  wptr = uint32_writex(tbuf, uii + 1, ' ');
	  wptr = double_g_writex(wptr, g_maxt_extreme_stat[uii], '\n');
	  if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    if (model_adapt) {
      sprintf(tbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(tbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(g_maxt_extreme_stat, &(g_maxt_extreme_stat[perms_total]));
#else
      qsort(g_maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
    }
    /*
    if (model_maxt) {
      printf("extreme stats: %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[perms_total - 1]);
    }
    */
    fprintf(outfile, tbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_non_set_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    while (1) {
      while (1) {
	do {
          chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	g_is_x = (uii == (uint32_t)x_code);
	if (model_assoc || ((!IS_SET(haploid_mask, uii)) || g_is_x)) {
	  break;
	}
	marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
      }
      wptr_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, uii, zero_extra_chroms));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	if (model_adapt) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) * dxx;
	}
        if (pval <= pfilter) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	  wptr = &(wptr_start[1 + plink_maxsnp]);
	  if ((!model_assoc) && ((model_adapt && (!g_perm_attempt_ct[marker_idx])) || ((!model_adapt) && ((model_fisherx && (g_orig_1mpval[marker_idx] == -9)) || ((!model_fisherx) && (g_orig_chisq[marker_idx] == -9)))))) {
	    // invalid
            wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    if (!model_perm_count) {
	      wptr = double_g_writewx4x(wptr, pval, 12, ' ');
	    } else {
	      wptr = double_g_writewx4x(wptr, ((double)g_perm_2success_ct[marker_idx]) / 2.0, 12, ' ');
	    }
	    if (model_adapt) {
	      if (g_perm_attempt_ct[marker_idx]) {
		wptr = memseta(wptr, 32, 2);
		wptr = uint32_writew10(wptr, g_perm_attempt_ct[marker_idx]);
	      }
	    } else {
	      if (model_fisherx) {
		// minimum p-value
		dzz = (int32_t)(doublearr_greater_than(g_maxt_extreme_stat, perms_total, 1.0 + EPSILON - g_orig_1mpval[marker_idx]) + 1);
	      } else {
		// maximum chisq
		dzz = (int32_t)(perms_total - doublearr_greater_than(g_maxt_extreme_stat, perms_total, g_orig_chisq[marker_idx] - EPSILON) + 1);
	      }
	      if (!model_perm_count) {
		wptr = double_g_writewx4(wptr, dzz * dyy, 12);
	      } else {
		wptr = double_g_writewx4(wptr, dzz, 12);
	      }
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto model_assoc_loop_end;
	}
        marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      }
    }
  model_assoc_loop_end:
    if (fclose_null(&outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "Permutation test report written to %s.\n", outname);
    logprintb();
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
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 model_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
}

int32_t qassoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_mperm_val, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctv2 = 2 * ((unfiltered_indiv_ct + BITCT - 1) / BITCT);
  uintptr_t pheno_nm_ctv2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  FILE* outfile = NULL;
  FILE* outfile_qtm = NULL;
  FILE* outfile_msa = NULL;
  uint32_t perm_adapt = model_modifier & MODEL_PERM;
  uint32_t perm_maxt = model_modifier & MODEL_MPERM;
  uint32_t do_perms = perm_adapt | perm_maxt;
  uint32_t qt_means = model_modifier & MODEL_QT_MEANS;
  uint32_t do_lin = model_modifier & MODEL_LIN;
  uint32_t qt_means_or_lin = qt_means || do_lin;
  uint32_t perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t fill_orig_chiabs = do_perms || mtest_adjust;
  uint32_t perms_total = 0;
  uint32_t pct = 0;
  uint32_t perm_pass_idx = 0;
  uintptr_t perm_vec_ctcl8m = 0;
  int32_t retval = 0;
  double x11 = 0;
  double x12 = 0;
  double x22 = 0;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t* tcnt = NULL;
  char* chrom_name_ptr = NULL;
  uint32_t chrom_name_len = 0;
  char chrom_name_buf[4];
  uint32_t mu_table[MODEL_BLOCKSIZE];
  uint32_t uibuf[4];
  char* outname_end2;
  char* wptr_start;
  char* wptr;
  char* wptr_restart;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* cur_loadbuf;
  uintptr_t* lbptr2;
  uintptr_t* indiv_include2;
  double* ooptr;
  double* dptr;
  double* dptr2;
  double* dptr3;
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
  uintptr_t indiv_uidx;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  intptr_t geno_sum;
  intptr_t nanal;
  double nanal_recip;
  double qt_sum;
  double qt_ssq;
  intptr_t geno_ssq;
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
  uint32_t homrar_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uintptr_t ulii;
  double dxx;
  double dyy;
  double dzz;
  double pval;
  uint32_t loop_end;
  char* a1ptr;
  char* a2ptr;
  if (pheno_nm_ct < 2) {
    logprint("Warning: Skipping QT --assoc since less than two phenotypes are present.\n");
    goto qassoc_ret_1;
  }
  g_pheno_nm_ct = pheno_nm_ct;
  g_aperm_alpha = aperm_alpha;
  g_perms_done = 0;
  g_mperm_save_all = NULL;
  if (perm_maxt) {
    perms_total = model_mperm_val;
    if (wkspace_alloc_d_checked(&g_maxt_extreme_stat, sizeof(double) * perms_total)) {
      goto qassoc_ret_NOMEM;
    }
    fill_double_zero(g_maxt_extreme_stat, perms_total); // square of t-stat
    g_ldrefs = (uint16_t*)wkspace_alloc(marker_ct * sizeof(uint16_t));
    if (!g_ldrefs) {
      goto qassoc_ret_NOMEM;
    }
#ifdef __LP64__
    fill_ulong_one((uintptr_t*)g_ldrefs, (marker_ct + 3) / 4);
#else
    fill_ulong_one((uintptr_t*)g_ldrefs, (marker_ct + 1) / 2);
#endif
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(&outfile_msa, outname, "w")) {
	goto qassoc_ret_OPEN_FAIL;
      }
      if (putc_checked('0', outfile_msa)) {
	goto qassoc_ret_WRITE_FAIL;
      }
      sprintf(logbuf, "Dumping all permutation squared %sstats to %s.\n", do_lin? "Lin " : "Wald t-", outname);
      logprintb();
    }
  } else {
    mperm_save = 0;
    if (perm_adapt) {
      perms_total = aperm_max;
      if (wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_ct * sizeof(uint32_t)) ||
	  wkspace_alloc_uc_checked(&g_perm_adapt_stop, marker_ct)) {
	goto qassoc_ret_NOMEM;
      }
      for (uii = 0; uii < marker_ct; uii++) {
	g_perm_attempt_ct[uii] = aperm_max;
      }
      fill_ulong_zero((uintptr_t*)g_perm_adapt_stop, (marker_ct + sizeof(uintptr_t) - 1) / sizeof(uintptr_t));
    }
  }
  outname_end2 = memcpyb(outname_end, ".qassoc", 8);
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto qassoc_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctv2 - 2] = 0;
  loadbuf_raw[unfiltered_indiv_ctv2 - 1] = 0;
  cur_loadbuf = loadbuf_raw;
  if (fill_orig_chiabs) {
    if (wkspace_alloc_d_checked(&g_orig_chisq, marker_ct * sizeof(double))) {
      goto qassoc_ret_NOMEM;
    }
    if (mtest_adjust) {
      if (wkspace_alloc_ui_checked(&tcnt, marker_ct * sizeof(int32_t))) {
	goto qassoc_ret_NOMEM;
      }
    }
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto qassoc_ret_OPEN_FAIL;
  }
  if (qt_means) {
    memcpy(outname_end2, ".means", 7);
    if (fopen_checked(&outfile_qtm, outname, "w")) {
      goto qassoc_ret_OPEN_FAIL;
    }
    sprintf(tbuf, " CHR %%%us  VALUE      G11      G12      G22\n", plink_maxsnp);
    fprintf(outfile_qtm, tbuf, "SNP");
    *outname_end2 = '\0';
  }
  if (haploid_chrom_present(chrom_info_ptr)) {
    logprint("Warning: QT --assoc doesn't handle X/Y/haploid markers normally (try --linear).\n");
  }
  sprintf(logbuf, "Writing QT --assoc report to %s...", outname);
  logprintb();
  fflush(stdout);
  sprintf(tbuf, " CHR %%%us         BP    NMISS       BETA         SE         R2        T            P ", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  if (do_lin) {
    fputs("         LIN        LIN_P ", outfile);
  }
  if (putc_checked('\n', outfile)) {
    goto qassoc_ret_WRITE_FAIL;
  }
  if (do_perms) {
    if (cluster_starts) {
      retval = cluster_include_and_reindex(unfiltered_indiv_ct, pheno_nm, 1, NULL, pheno_nm_ct, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, NULL, NULL);
      if (retval) {
	goto qassoc_ret_1;
      }
      if (!g_cluster_ct) {
        logprint("Error: No size 2+ clusters for permutation test.\n");
        goto qassoc_ret_INVALID_CMDLINE;
      }
      if (wkspace_alloc_ui_checked(&g_indiv_to_cluster, pheno_nm_ct * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&g_qassoc_cluster_thread_wkspace, g_thread_ct * ((g_cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32) * CACHELINE)) {
	goto qassoc_ret_NOMEM;
      }
      fill_unfiltered_indiv_to_cluster(pheno_nm_ct, g_cluster_ct, g_cluster_map, g_cluster_starts, g_indiv_to_cluster);
    }
    g_sfmtp_arr = (sfmt_t**)wkspace_alloc(g_thread_ct * sizeof(intptr_t));
    if (!g_sfmtp_arr) {
      goto qassoc_ret_NOMEM;
    }
    g_sfmtp_arr[0] = &sfmt;
    if (g_thread_ct > 1) {
      for (uii = 1; uii < g_thread_ct; uii++) {
	g_sfmtp_arr[uii] = (sfmt_t*)wkspace_alloc(sizeof(sfmt_t));
	if (!g_sfmtp_arr[uii]) {
	  goto qassoc_ret_NOMEM;
	}
	for (ujj = 0; ujj < 4; ujj++) {
	  uibuf[ujj] = sfmt_genrand_uint32(&sfmt);
	}
	sfmt_init_by_array(g_sfmtp_arr[uii], uibuf, 4);
      }
    }
    if (wkspace_alloc_ui_checked(&g_missing_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_het_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_homcom_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_ct * sizeof(uint32_t))) {
      goto qassoc_ret_NOMEM;
    }
    fill_uint_zero(g_perm_2success_ct, marker_ct);
  }
  if (do_lin) {
    if (wkspace_alloc_d_checked(&g_orig_linsq, marker_ct * sizeof(double))) {
      goto qassoc_ret_NOMEM;
    }
  }
  g_adaptive_ci_zt = ltqnorm(1 - aperm_beta / (2.0 * marker_ct));
  if (wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * pheno_nm_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&g_orig_1mpval, marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&g_marker_uidxs, marker_ct * sizeof(uint32_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, pheno_nm_ctv2 * sizeof(intptr_t))) {
    goto qassoc_ret_NOMEM;
  }
  fill_vec_55(indiv_include2, pheno_nm_ct);
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, pheno_nm_ct, hh_exists, 1, pheno_nm, sex_male, &indiv_include2, &indiv_male_include2)) {
    goto qassoc_ret_NOMEM;
  }
  marker_unstopped_ct = marker_ct;
  if (wkspace_alloc_d_checked(&g_pheno_d2, pheno_nm_ct * sizeof(double))) {
    goto qassoc_ret_NOMEM;
  }
  g_pheno_sum = 0;
  g_pheno_ssq = 0;
  indiv_uidx = 0;
  indiv_idx = 0;
  dptr = g_pheno_d2;
  do {
    indiv_uidx = next_set_ul_unsafe(pheno_nm, indiv_uidx);
    indiv_uidx_stop = next_unset_ul(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    dptr2 = &(pheno_d[indiv_uidx]);
    indiv_uidx = indiv_uidx_stop;
    dptr3 = &(pheno_d[indiv_uidx_stop]);
    do {
      dxx = *dptr2++;
      *dptr++ = dxx;
      g_pheno_sum += dxx;
      g_pheno_ssq += dxx * dxx;
    } while (dptr2 < dptr3);
  } while (indiv_idx < pheno_nm_ct);
  fputs(" 0%", stdout);
  fflush(stdout);

  // ----- begin main loop -----
 qassoc_more_perms:
  if (do_perms) {
    if (perm_adapt) {
      if (perm_pass_idx) {
	while (g_first_adapt_check <= g_perms_done) {
	  // APERM_MAX prevents infinite loop here
	  g_first_adapt_check += (int32_t)(aperm_init_interval + ((int32_t)g_first_adapt_check) * aperm_interval_slope);
	}
      } else {
	if (aperm_min < aperm_init_interval) {
	  g_first_adapt_check = (int32_t)aperm_init_interval;
	} else {
	  g_first_adapt_check = aperm_min;
	}
	g_adaptive_intercept = aperm_init_interval;
	g_adaptive_slope = aperm_interval_slope;
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
    if (perm_batch_size) {
      g_perm_vec_ct = perm_batch_size;
    } else {
      // this seems to work better than pregenerating as many permutations as
      // possible.  It might be best to auto-tune this, but we need performance
      // data across a wide variety of machines to do this intelligently.
      g_perm_vec_ct = 512;
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    }
    perm_vec_ctcl8m = (g_perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
    if (wkspace_alloc_d_checked(&g_perm_vecstd, perm_vec_ctcl8m * sizeof(double) * pheno_nm_ct)) {
      goto qassoc_ret_NOMEM;
    }
    ulii = do_lin? 6 : 3;
    if (perm_maxt) {
      if (wkspace_alloc_d_checked(&g_maxt_thread_results, g_thread_ct * perm_vec_ctcl8m * sizeof(double)) ||
	  wkspace_alloc_d_checked(&g_qresultbuf, ulii * MODEL_BLOCKSIZE * perm_vec_ctcl8m * sizeof(double))) {
	goto qassoc_ret_NOMEM;
      }
      if (mperm_save & MPERM_DUMP_ALL) {
	if (wkspace_alloc_d_checked(&g_mperm_save_all, marker_ct * sizeof(double) * g_perm_vec_ct)) {
	  goto qassoc_ret_NOMEM;
	}
      }
    } else {
      if (wkspace_alloc_d_checked(&g_thread_git_qbufs, perm_vec_ctcl8m * sizeof(double) * ulii * g_thread_ct)) {
	goto qassoc_ret_NOMEM;
      }
      fill_double_zero(g_thread_git_qbufs, ulii * g_thread_ct * perm_vec_ctcl8m);
    }
    g_perms_done += g_perm_vec_ct;
    if (g_perm_vec_ct >= CACHELINE_DBL * g_thread_ct) {
      g_assoc_thread_ct = g_thread_ct;
    } else {
      g_assoc_thread_ct = g_perm_vec_ct / CACHELINE_DBL;
      if (!g_assoc_thread_ct) {
	g_assoc_thread_ct = 1;
      }
    }
    ulii = 0;
    if (!cluster_starts) {
      if (spawn_threads(threads, &qassoc_gen_perms_thread, g_assoc_thread_ct)) {
	goto qassoc_ret_THREAD_CREATE_FAIL;
      }
      qassoc_gen_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &qassoc_gen_cluster_perms_thread, g_assoc_thread_ct)) {
	goto qassoc_ret_THREAD_CREATE_FAIL;
      }
      qassoc_gen_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_assoc_thread_ct);
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
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
      refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &g_is_haploid);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      chrom_name_ptr = chrom_name_buf;
      chrom_name_len = 4;
      if (uii <= chrom_info_ptr->max_code) {
	memset(chrom_name_buf, 32, 2);
        intprint2(&(chrom_name_buf[2]), uii);
      } else if (zero_extra_chroms) {
	memcpy(chrom_name_buf, "   0", 4);
      } else {
	ujj = strlen(chrom_info_ptr->nonstd_names[uii]);
	if (ujj < 4) {
	  fw_strcpyn(4, ujj, chrom_info_ptr->nonstd_names[uii], chrom_name_buf);
	} else {
	  chrom_name_ptr = chrom_info_ptr->nonstd_names[uii];
	  chrom_name_len = ujj;
	}
      }
    } else if (perm_maxt) {
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
      if (perm_adapt && g_perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx + 1);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto qassoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(g_loadbuf[block_size * pheno_nm_ctv2]);
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_ptr, pheno_nm_ct, pheno_nm, IS_SET(marker_reverse, marker_uidx))) {
	goto qassoc_ret_READ_FAIL;
      }
      if (g_is_haploid && hh_exists) {
	haploid_fix(hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, g_is_x, g_is_y, (unsigned char*)loadbuf_ptr);
      }
      if (perm_adapt) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
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
        g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	loadbuf_ptr = &(g_loadbuf[marker_bidx * pheno_nm_ctv2]);
	vec_3freq(pheno_nm_ctv2, loadbuf_ptr, indiv_include2, &missing_ct, &het_ct, &homcom_ct);
	nanal = pheno_nm_ct - missing_ct;
	wptr = memcpya(tbuf, chrom_name_ptr, chrom_name_len);
	*wptr++ = ' ';
        wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	*wptr++ = ' ';
	wptr = uint32_writew10(wptr, marker_pos[marker_uidx2]);
	*wptr++ = ' ';
	wptr = uint32_writew8(wptr, nanal);
	*wptr++ = ' ';
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
	    indiv_idx = uii + (ujj / 2);
	    dxx = g_pheno_d2[indiv_idx];
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
	    ulii &= ~(3 * (ONELU << ujj));
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
	  if (mtest_adjust) {
	    tcnt[marker_idx + marker_bidx] = nanal;
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
		fprintf(outfile_msa, " %g", tstat * tstat);
	      } else {
		fputs(" NA", outfile_msa);
	      }
	    } else {
	      dxx = g_orig_linsq[marker_idx + marker_bidx];
	      if ((nanal > 2) && realnum(dxx)) {
		fprintf(outfile_msa, " %g", dxx);
	      } else {
		fputs(" NA", outfile_msa);
	      }
	    }
	  }
	  if ((pfilter == 1.0) || ((tp != -9) && (tp <= pfilter))) {
	    if (!realnum(beta)) {
	      wptr = memcpya(wptr, "        NA         NA         NA ", 33);
	    } else {
	      wptr = double_g_writewx4x(wptr, beta, 10, ' ');
	      wptr = double_g_writewx4x(wptr, vbeta_sqrt, 10, ' ');
	      wptr = double_g_writewx4x(wptr, rsq, 10, ' ');
	    }
	    if (tp >= 0) {
	      wptr = double_g_writewx4x(wptr, tstat, 8, ' ');
	      wptr = double_g_writewx4(wptr, tp, 12);
	    } else {
	      wptr = memcpya(wptr, "      NA           NA", 21);
	    }
	    if (do_lin && (nanal > 2)) {
	      dxx = g_orig_linsq[marker_idx + marker_bidx];
	      if (realnum(dxx)) {
		*wptr++ = ' ';
		dxx = sqrt(dxx);
		wptr = double_g_writewx4x(wptr, dxx, 12, ' ');
		wptr = double_g_writewx4(wptr, calc_tprob(dxx, nanal - 2), 12);
	      } else {
		wptr = memcpya(wptr, "           NA           NA", 26);
	      }
	    }
	    wptr = memcpya(wptr, " \n", 2);
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto qassoc_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  wptr = memcpya(wptr, "        NA         NA         NA       NA           NA ", 55);
	  if (mperm_save & MPERM_DUMP_ALL) {
	    fputs(" NA", outfile_msa);
	  }
	  if (do_lin) {
	    wptr = memcpya(wptr, "          NA           NA ", 26);
	  }
	  *wptr++ = '\n';
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
	if (qt_means) {
	  wptr_restart = &(tbuf[2 + chrom_name_len + plink_maxsnp]);
	  wptr = memcpya(wptr_restart, "  GENO ", 7);
	  a1ptr = &(marker_alleles[(2 * marker_uidx2) * max_marker_allele_len]);
	  a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1) * max_marker_allele_len]);
	  if (max_marker_allele_len == 1) {
	    memset(wptr, 32, 5);
	    wptr[5] = *a1ptr;
	    wptr[6] = '/';
	    wptr[7] = *a1ptr;
	    memset(&(wptr[8]), 32, 6);
            wptr[14] = *a1ptr;
	    wptr[15] = '/';
	    wptr[16] = *a2ptr;
	    memset(&(wptr[17]), 32, 6);
	    wptr[23] = *a2ptr;
	    wptr[24] = '/';
	    wptr[25] = *a2ptr;
	    wptr = &(wptr[26]);
	  } else {
	    uii = strlen(a1ptr);
	    ujj = strlen(a2ptr);
	    if (uii < 4) {
	      wptr = memseta(wptr, 32, 7 - 2 * uii);
	    }
	    wptr = memcpyax(wptr, a1ptr, uii, '/');
	    wptr = memcpyax(wptr, a1ptr, uii, ' ');
	    if (uii + ujj < 7) {
	      wptr = memseta(wptr, 32, 7 - uii - ujj);
	    }
	    wptr = memcpyax(wptr, a1ptr, uii, '/');
	    wptr = memcpyax(wptr, a2ptr, ujj, ' ');
	    if (ujj < 4) {
	      wptr = memseta(wptr, 32, 7 - 2 * ujj);
	    }
	    wptr = memcpyax(wptr, a2ptr, ujj, '/');
	    wptr = memcpya(wptr, a2ptr, ujj);
	  }
	  *wptr++ = '\n';
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "COUNTS ", 7);
	  wptr = uint32_writew8(wptr, homrar_ct);
	  *wptr++ = ' ';
	  wptr = uint32_writew8(wptr, het_ct);
	  *wptr++ = ' ';
	  wptr = uint32_writew8(wptr, homcom_ct);
	  *wptr++ = '\n';
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "  FREQ ", 7);
	  wptr = double_g_writewx4x(wptr, nanal_recip * ((intptr_t)homrar_ct), 8, ' ');
	  wptr = double_g_writewx4x(wptr, nanal_recip * ((intptr_t)het_ct), 8, ' ');
	  wptr = double_g_writewx4x(wptr, nanal_recip * ((intptr_t)homcom_ct), 8, '\n');
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "  MEAN ", 7);
	  qt_homcom_sum = qt_sum - qt_homrar_sum - qt_het_sum;
	  if (homrar_ct) {
	    x11 = qt_homrar_sum / ((double)homrar_ct);
	    wptr = double_g_writewx4(wptr, x11, 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (het_ct) {
	    x12 = qt_het_sum / ((double)het_ct);
	    wptr = double_g_writewx4(wptr, x12, 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (homcom_ct) {
	    x22 = qt_homcom_sum / ((double)homcom_ct);
	    wptr = double_g_writewx4(wptr, x22, 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = '\n';
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	  wptr = memcpya(wptr_restart, "    SD ", 7);
	  if (homrar_ct > 1) {
            dxx = sqrt((qt_homrar_ssq - qt_homrar_sum * x11) / ((double)((intptr_t)homrar_ct - 1)));
	    wptr = double_g_writewx4(wptr, dxx, 8);
	  } else if (homrar_ct == 1) {
	    wptr = memcpya(wptr, "       0", 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (het_ct > 1) {
            dxx = sqrt((qt_het_ssq - qt_het_sum * x12) / ((double)((intptr_t)het_ct - 1)));
	    wptr = double_g_writewx4(wptr, dxx, 8);
	  } else if (het_ct == 1) {
	    wptr = memcpya(wptr, "       0", 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = ' ';
	  if (homcom_ct > 1) {
            dxx = sqrt((qt_ssq - qt_het_ssq - qt_homrar_ssq - qt_homcom_sum * x22) / ((double)((intptr_t)homcom_ct - 1)));
	    wptr = double_g_writewx4(wptr, dxx, 8);
	  } else if (homcom_ct == 1) {
	    wptr = memcpya(wptr, "       0", 8);
	  } else {
	    wptr = memcpya(wptr, "      NA", 8);
	  }
	  *wptr++ = '\n';
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile_qtm)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
      }
    }
    if (do_perms) {
      g_block_diff = block_size - g_qblock_start;
      g_assoc_thread_ct = g_block_diff / CACHELINE_DBL;
      if (g_assoc_thread_ct > g_thread_ct) {
	g_assoc_thread_ct = g_thread_ct;
      } else if (!g_assoc_thread_ct) {
	g_assoc_thread_ct = 1;
      }
      ulii = 0;
      if (perm_maxt) {
	g_maxt_block_base = marker_idx;
	g_maxt_cur_extreme_stat = g_maxt_extreme_stat[0];
	for (uii = 1; uii < g_perm_vec_ct; uii++) {
	  dxx = g_maxt_extreme_stat[uii];
	  if (dxx > g_maxt_cur_extreme_stat) {
	    g_maxt_cur_extreme_stat = dxx;
	  }
	}
	if (!do_lin) {
	  if (spawn_threads(threads, &qassoc_maxt_thread, g_assoc_thread_ct)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_maxt_thread((void*)ulii);
	} else {
	  if (spawn_threads(threads, &qassoc_maxt_lin_thread, g_assoc_thread_ct)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_maxt_lin_thread((void*)ulii);
	}
        join_threads(threads, g_assoc_thread_ct);
	ulii = CACHELINE_DBL * ((g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL);
	for (uii = 0; uii < g_assoc_thread_ct; uii++) {
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
	  if (spawn_threads(threads, &qassoc_adapt_thread, g_assoc_thread_ct)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_adapt_thread((void*)ulii);
	} else {
	  if (spawn_threads(threads, &qassoc_adapt_lin_thread, g_assoc_thread_ct)) {
	    goto qassoc_ret_THREAD_CREATE_FAIL;
	  }
	  qassoc_adapt_lin_thread((void*)ulii);
	}
        join_threads(threads, g_assoc_thread_ct);
      }
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_unstopped_ct) {
	if (pct >= 10) {
	  putchar('\b');
	}
	pct = (marker_idx * 100LLU) / marker_unstopped_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = ((uint64_t)(pct + 1LLU) * marker_unstopped_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putchar('\b');
    }
    fputs("\b\b\b", stdout);
    logprint(" done.\n");
    if (qt_means) {
      sprintf(logbuf, "QT means report saved to %s.means.\n", outname);
      logprintb();
      if (fclose_null(&outfile_qtm)) {
	goto qassoc_ret_WRITE_FAIL;
      }
    }
    if (do_perms) {
      wkspace_reset((unsigned char*)g_perm_vecstd);
    }
    if (fclose_null(&outfile)) {
      goto qassoc_ret_WRITE_FAIL;
    }
    if (mtest_adjust) {
      if (do_lin) {
	for (uii = 0; uii < marker_ct; uii++) {
	  g_orig_chisq[uii] = sqrt(g_orig_linsq[uii]);
	}
      }
      retval = multcomp(outname, outname_end, g_marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, g_orig_chisq, pfilter, mtest_adjust, adjust_lambda, 1, tcnt);
      if (retval) {
	goto qassoc_ret_1;
      }
    }
    if (mperm_save & MPERM_DUMP_ALL) {
      if (putc_checked('\n', outfile_msa)) {
	goto qassoc_ret_WRITE_FAIL;
      }
    }
  }
  if (do_perms) {
    if (mperm_save & MPERM_DUMP_ALL) {
      if (perm_pass_idx) {
	putchar(' ');
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = g_perm_vec_ct;
      ujj = 1 + g_perms_done - ulii;
      wptr = tbuf;
      a1ptr = &(tbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32_write(wptr, uii + ujj);
	ooptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = ooptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = double_g_write(wptr, dxx);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= a1ptr) {
	    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	      goto qassoc_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	goto qassoc_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    wkspace_reset((unsigned char*)g_perm_vecstd);
    if (g_perms_done < perms_total) {
      if (perm_adapt) {
	marker_unstopped_ct = marker_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, 0, (marker_ct + sizeof(uintptr_t) - 1) / sizeof(uintptr_t));
	if (!marker_unstopped_ct) {
	  goto qassoc_adapt_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done > 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto qassoc_more_perms;
    }
    if (perm_adapt) {
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
    putchar('\r');
    sprintf(logbuf, "%u %s permutation%s complete.\n", g_perms_done, perm_maxt? "max(T)" : "(adaptive)", (g_perms_done > 1)? "s" : "");
    logprintb();

    if (perm_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	memcpy(outname_end, ".mperm.dump.best", 17);
	sprintf(logbuf, "Dumping best permutation squared %sstats to %s.\n", do_lin? "Lin " : "Wald t-", outname);
	logprintb();
	if (fopen_checked(&outfile, outname, "w")) {
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
        memcpy(tbuf, "0 ", 2);
	wptr = double_g_writex(&(tbuf[2]), dxx, '\n');
	if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	  goto qassoc_ret_WRITE_FAIL;
	}
	for (uii = 0; uii < perms_total; uii++) {
	  wptr = uint32_writex(tbuf, uii + 1, ' ');
	  wptr = double_g_writex(wptr, g_maxt_extreme_stat[uii], '\n');
	  if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto qassoc_ret_WRITE_FAIL;
	}
	memcpy(outname_end, ".qassoc", 7);
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto qassoc_ret_OPEN_FAIL;
    }
    if (perm_adapt) {
      sprintf(tbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(tbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
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
    fprintf(outfile, tbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_non_set_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, uii, zero_extra_chroms));
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
	  if (perm_adapt && (!g_perm_attempt_ct[marker_idx])) {
	    // invalid
            wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    if (!perm_count) {
	      wptr = double_g_writewx4x(wptr, pval, 12, ' ');
	    } else {
	      wptr = double_g_writewx4x(wptr, ((double)g_perm_2success_ct[marker_idx]) / 2.0, 12, ' ');
	    }
	    if (perm_adapt) {
	      if (g_perm_attempt_ct[marker_idx]) {
		wptr = memseta(wptr, 32, 2);
		wptr = uint32_writew10(wptr, g_perm_attempt_ct[marker_idx]);
	      }
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
		wptr = double_g_writewx4(wptr, dzz * dyy, 12);
	      } else {
		wptr = double_g_writewx4(wptr, dzz, 12);
	      }
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto qassoc_loop_end;
	}
        marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      }
    }
  qassoc_loop_end:
    if (fclose_null(&outfile)) {
      goto qassoc_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "Permutation test report written to %s.\n", outname);
    logprintb();

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
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 qassoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_qtm);
  fclose_cond(outfile_msa);
  return retval;
}

int32_t gxe_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* gxe_covar_nm, uintptr_t* gxe_covar_c, uintptr_t* sex_male, uint32_t hh_exists) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t covar_nm_ct = popcount_longs(gxe_covar_nm, 0, indiv_ctl);
  uintptr_t covar_nm_ctl = (covar_nm_ct + (BITCT - 1)) / BITCT;
  // gxe_covar_c has opposite truth value from ->bcovar in PLINK gxe.cpp; see
  // lines 50-58 in gxe.cpp
  uintptr_t group2_size = popcount_longs(gxe_covar_c, 0, indiv_ctl);
  uintptr_t group1_size = covar_nm_ct - group2_size;
  uintptr_t male_ct = 0;
  uintptr_t male_ctl = 0;
  uintptr_t group1_size_male = 0;
  uintptr_t group2_size_male = 0;
  uintptr_t marker_uidx = 0;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* indiv_male_all_include2 = NULL;
  uintptr_t* group1_include2 = NULL;
  uintptr_t* group2_include2 = NULL;
  uintptr_t* group1_male_include2 = NULL;
  uintptr_t* group2_male_include2 = NULL;
  uintptr_t* covar_nm_raw = NULL;
  uintptr_t* covar_nm_male_raw = NULL;
  uintptr_t* cur_indiv_i2 = NULL;
  uintptr_t* cur_indiv_male_i2 = NULL;
  uintptr_t* cur_group1_i2 = NULL;
  uintptr_t* cur_group2_i2 = NULL;
  uintptr_t* cur_covar_nm_raw = NULL;
  double* pheno_d_collapsed = NULL;
  double* pheno_d_male_collapsed = NULL;
  double* cur_pheno_d = NULL;
  char* wptr_start = NULL;
  uintptr_t cur_indiv_ct = 0;
  uintptr_t cur_indiv_ctl2 = 0;
  uintptr_t cur_group1_size = 0;
  uintptr_t cur_group2_size = 0;
  uint32_t y_exists = (chrom_info_ptr->y_code != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->y_code);
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
  uintptr_t indiv_uidx;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t indiv_idx2;
  uintptr_t indiv_idx2_offset;
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
  uint32_t is_haploid;
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
    logprint("Error: First --gxe group has fewer than three members.\n");
    goto gxe_assoc_ret_INVALID_CMDLINE;
  } else if (group2_size < 3) {
    logprint("Error: Second --gxe group has fewer than three members.\n");
    goto gxe_assoc_ret_INVALID_CMDLINE;
  }
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctl * 2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, covar_nm_ctl * 2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&covar_nm_raw, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&pheno_d_collapsed, covar_nm_ct * sizeof(double))) {
    goto gxe_assoc_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctl * 2 - 2] = 0;
  loadbuf_raw[unfiltered_indiv_ctl * 2 - 1] = 0;

  fill_ulong_zero(covar_nm_raw, unfiltered_indiv_ctl);
  indiv_uidx = 0;
  indiv_idx = 0;
  indiv_idx2 = 0;
  do {
    indiv_uidx = next_unset_ul_unsafe(indiv_exclude, indiv_uidx);
    indiv_uidx_stop = next_set_ul(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
    do {
      if (IS_SET(gxe_covar_nm, indiv_idx)) {
        SET_BIT(covar_nm_raw, indiv_uidx);
        dxx = pheno_d[indiv_uidx];
        if (IS_SET(gxe_covar_c, indiv_idx)) {
	  pheno_sum_g2 += dxx;
	  pheno_ssq_g2 += dxx * dxx;
	} else {
	  pheno_sum_g1 += dxx;
	  pheno_ssq_g1 += dxx * dxx;
	}
	pheno_d_collapsed[indiv_idx2++] = dxx;
      }
      indiv_idx++;
    } while (++indiv_uidx < indiv_uidx_stop);
  } while (indiv_idx < indiv_ct);

  if (wkspace_alloc_ul_checked(&group1_include2, covar_nm_ctl * 2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&group2_include2, covar_nm_ctl * 2 * sizeof(intptr_t))) {
    goto gxe_assoc_ret_NOMEM;
  }
  fill_vec_55(group1_include2, covar_nm_ct);
  fill_ulong_zero(group2_include2, covar_nm_ctl * 2);
  indiv_idx = 0;
  indiv_idx2 = 0;
  do {
    indiv_idx = next_set_ul_unsafe(gxe_covar_nm, indiv_idx);
    indiv_uidx_stop = next_unset_ul(gxe_covar_nm, indiv_idx, indiv_ct);
    do {
      if (IS_SET(gxe_covar_c, indiv_idx)) {
	SET_BIT_DBL(group2_include2, indiv_idx2);
      }
      indiv_idx2++;
    } while (++indiv_idx < indiv_uidx_stop);
  } while (indiv_idx2 < covar_nm_ct);
  bitfield_andnot(group1_include2, group2_include2, covar_nm_ctl * 2);

  if ((hh_exists & NXMHH_EXISTS) || y_exists) {
    if (wkspace_alloc_ul_checked(&indiv_include2, covar_nm_ctl * 2 * sizeof(intptr_t))) {
      goto gxe_assoc_ret_NOMEM;
    }
    fill_vec_55(indiv_include2, covar_nm_ct);
  }
  if ((hh_exists & XMHH_EXISTS) || y_exists) {
    if (wkspace_alloc_ul_checked(&indiv_male_include2, covar_nm_ctl * 2 * sizeof(intptr_t))) {
      goto gxe_assoc_ret_NOMEM;
    }
    fill_ulong_zero(indiv_male_include2, covar_nm_ctl * 2);
    indiv_uidx = 0;
    indiv_idx = 0;
    indiv_idx2 = 0;
    do {
      indiv_uidx = next_unset_ul_unsafe(indiv_exclude, indiv_uidx);
      indiv_uidx_stop = next_set_ul(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
      do {
        if (IS_SET(gxe_covar_nm, indiv_idx)) {
          if (IS_SET(sex_male, indiv_uidx)) {
	    SET_BIT_DBL(indiv_male_include2, indiv_idx2);
	    male_ct++;
	  }
	  indiv_idx2++;
	}
	indiv_idx++;
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_ct);
    male_ctl = (male_ct + (BITCT - 1)) / BITCT;
    if (y_exists) {
      group1_size_male = popcount_longs_exclude(indiv_male_include2, group2_include2, covar_nm_ctl * 2);
      group2_size_male = male_ct - group1_size_male;
      if ((group1_size_male < 3) || (group2_size_male < 3)) {
        logprint("Warning: Skipping Y chromosome for --gxe since a group has less than 3 males.\n");
	skip_y = 1;
      }
      // currently still need to initialize covar_nm_male_raw even on skip_y
      if (wkspace_alloc_ul_checked(&indiv_male_all_include2, male_ctl * 2 * sizeof(intptr_t)) ||
          wkspace_alloc_ul_checked(&group1_male_include2, male_ctl * 2 * sizeof(intptr_t)) ||
	  wkspace_alloc_ul_checked(&group2_male_include2, male_ctl * 2 * sizeof(intptr_t)) ||
	  wkspace_alloc_d_checked(&pheno_d_male_collapsed, male_ct * sizeof(double)) ||
	  wkspace_alloc_ul_checked(&covar_nm_male_raw, unfiltered_indiv_ctl * sizeof(intptr_t))) {
	goto gxe_assoc_ret_NOMEM;
      }
      fill_vec_55(indiv_male_all_include2, male_ct);
      fill_vec_55(group1_male_include2, male_ct);
      fill_ulong_zero(group2_male_include2, male_ctl * 2);
      indiv_idx = 0;
      for (indiv_idx2 = 0; indiv_idx2 < covar_nm_ct; indiv_idx2++) {
	if (IS_SET_DBL(indiv_male_include2, indiv_idx2)) {
	  dxx = pheno_d_collapsed[indiv_idx2];
	  if (IS_SET_DBL(group2_include2, indiv_idx2)) {
	    SET_BIT_DBL(group2_male_include2, indiv_idx);
	    pheno_sum_male_g2 += dxx;
	    pheno_ssq_male_g2 += dxx * dxx;
	  } else {
	    pheno_sum_male_g1 += dxx;
            pheno_ssq_male_g1 += dxx * dxx;
	  }
	  pheno_d_male_collapsed[indiv_idx++] = dxx;
	}
      }
      bitfield_andnot(group1_male_include2, group2_male_include2, male_ctl * 2);
      for (ulii = 0; ulii < unfiltered_indiv_ctl; ulii++) {
	covar_nm_male_raw[ulii] = covar_nm_raw[ulii] & sex_male[ulii];
      }
    }
  }

  memcpy(outname_end, ".qassoc.gxe", 12);
  if (fopen_checked(&outfile, outname, "w")) {
    goto gxe_assoc_ret_OPEN_FAIL;
  }
  if (haploid_chrom_present(chrom_info_ptr)) {
    logprint("Warning: --gxe doesn't currently handle X/Y/haploid markers properly.\n");
  }
  sprintf(logbuf, "Writing --gxe report to %s...", outname);
  logprintb();
  fputs(" 0%", stdout);
  fflush(stdout);
  sprintf(tbuf, " CHR %%%us   NMISS1      BETA1        SE1   NMISS2      BETA2        SE2    Z_GXE        P_GXE \n", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");

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
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto gxe_assoc_ret_READ_FAIL;
	}
      }
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	if (!is_y) {
	  cur_indiv_ct = covar_nm_ct;
	  cur_group1_size = group1_size;
          cur_group2_size = group2_size;
	  base_pheno_sum_g1 = pheno_sum_g1;
	  base_pheno_ssq_g1 = pheno_ssq_g1;
          base_pheno_sum_g2 = pheno_sum_g2;
          base_pheno_ssq_g2 = pheno_ssq_g2;
          cur_indiv_i2 = indiv_include2;
          cur_indiv_male_i2 = indiv_male_include2;
	  cur_group1_i2 = group1_include2;
          cur_group2_i2 = group2_include2;
          cur_pheno_d = pheno_d_collapsed;
	  cur_covar_nm_raw = covar_nm_raw;
	} else {
	  cur_indiv_ct = male_ct;
	  cur_group1_size = group1_size_male;
          cur_group2_size = group2_size_male;
          base_pheno_sum_g1 = pheno_sum_male_g1;
	  base_pheno_ssq_g1 = pheno_ssq_male_g1;
          base_pheno_sum_g2 = pheno_sum_male_g2;
	  base_pheno_ssq_g2 = pheno_ssq_male_g2;
          cur_indiv_i2 = indiv_male_all_include2;
          cur_indiv_male_i2 = indiv_male_all_include2;
          cur_group1_i2 = group1_male_include2;
          cur_group2_i2 = group2_male_include2;
          cur_pheno_d = pheno_d_male_collapsed;
	  cur_covar_nm_raw = covar_nm_male_raw;
	}
	wptr_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
	*wptr_start++ = ' ';
	cur_indiv_ctl2 = ((cur_indiv_ct + (BITCT - 1)) / BITCT) * 2;
        loadbuf[cur_indiv_ctl2 - 2] = 0;
        loadbuf[cur_indiv_ctl2 - 1] = 0;
      }

      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf, cur_indiv_ct, cur_covar_nm_raw, IS_SET(marker_reverse, marker_uidx))) {
	goto gxe_assoc_ret_READ_FAIL;
      }
      if (is_y && skip_y) {
	marker_uidx++;
	continue;
      }
      if (is_haploid) {
	haploid_fix(hh_exists, cur_indiv_i2, cur_indiv_male_i2, cur_indiv_ct, is_x, is_y, (unsigned char*)loadbuf);
      }

      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
      *wptr++ = ' ';

      // We are interested in the following quantities:
      //   qt_var{1,2}: (qt_ssq - (qt_sum^2 / N)) / (N-1)
      //   g_var{1,2}: (geno_ssq - (geno_sum^2 / N)) / (N-1)
      //   qt_g_covar{1,2}: (qt_g_prod - ((qt_sum * geno_sum) / N)) / (N-1)

      single_marker_cc_3freqs(cur_indiv_ctl2, loadbuf, cur_group1_i2, cur_group2_i2, &homcom_ct1, &het_ct1, &missing_ct1, &homcom_ct2, &het_ct2, &missing_ct2);
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
	indiv_idx2_offset = 0;
	loadbuf_ptr = loadbuf;
	cgr_ptr = cur_group2_i2;
	do {
	  ulmm = ~(*loadbuf_ptr++);
	  if (indiv_idx2_offset + BITCT2 > cur_indiv_ct) {
	    ulmm &= (ONELU << ((cur_indiv_ct & (BITCT2 - 1)) * 2)) - ONELU;
	  }
	  if (ulmm) {
	    ulnn = (*cgr_ptr) * 3;
            ulii = ulmm & (~ulnn);
            while (ulii) {
	      uljj = CTZLU(ulii) & (BITCT - 2);
	      ulkk = (ulii >> uljj) & 3;
	      indiv_idx2 = indiv_idx2_offset + (uljj / 2);
	      dxx = cur_pheno_d[indiv_idx2];
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
	      ulii &= ~(3 * (ONELU << uljj));
	    }
	    ulii = ulmm & ulnn;
	    while (ulii) {
	      uljj = CTZLU(ulii) & (BITCT - 2);
	      ulkk = (ulii >> uljj) & 3;
	      indiv_idx2 = indiv_idx2_offset + (uljj / 2);
	      dxx = cur_pheno_d[indiv_idx2];
	      if (ulkk == 1) {
		qt_g_prod2 += dxx;
	      } else if (ulkk == 3) {
		qt_g_prod2 += 2 * dxx;
	      } else {
		qt_sum2 -= dxx;
		qt_ssq2 -= dxx * dxx;
	      }
	      ulii &= ~(3 * (ONELU << uljj));
	    }
	  }
	  cgr_ptr++;
	  indiv_idx2_offset += BITCT2;
	} while (indiv_idx2_offset < cur_indiv_ct);
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
        wptr = uint32_writew8x(wptr, nanal1, ' ');
        wptr = double_g_writewx4x(wptr, beta1, 10, ' ');
        wptr = double_g_writewx4x(wptr, sqrt(vbeta1), 10, ' ');
        wptr = uint32_writew8x(wptr, nanal2, ' ');
        wptr = double_g_writewx4x(wptr, beta2, 10, ' ');
        wptr = double_g_writewx4x(wptr, sqrt(vbeta2), 10, ' ');
        wptr = double_g_writewx4x(wptr, zval, 8, ' ');
        wptr = double_g_writewx4x(wptr, chiprob_p(zval * zval, 1), 12, '\n');
      } else {
      gxe_assoc_nan_line:
        wptr = memcpya(wptr, "      NA         NA         NA       NA         NA         NA       NA           NA\n", 84);
      }
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto gxe_assoc_ret_WRITE_FAIL;
      }
      marker_uidx++;
    }
    if (pct < 100) {
      if (pct > 10) {
        putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }

  if (pct >= 10) {
    putchar('\b');
  }
  fputs("\b\b\b", stdout);
  logprint(" done.\n");
  fclose_cond(outfile);

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
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

uint32_t glm_init_load_mask(uintptr_t* indiv_exclude, uintptr_t* pheno_nm, uintptr_t* covar_nm, uint32_t indiv_ct, uintptr_t unfiltered_indiv_ctv2, uintptr_t** load_mask_ptr) {
  uint32_t indiv_uidx = 0;
  uintptr_t* load_mask;
  uint32_t indiv_idx;
  if (wkspace_alloc_ul_checked(load_mask_ptr, unfiltered_indiv_ctv2 * (sizeof(intptr_t) / 2))) {
    return 1;
  }
  load_mask = *load_mask_ptr;
  fill_ulong_zero(load_mask, unfiltered_indiv_ctv2 / 2);
  if (covar_nm) {
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      if (IS_SET(indiv_exclude, indiv_uidx)) {
        indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
      }
      if (IS_SET(pheno_nm, indiv_uidx) & IS_SET(covar_nm, indiv_idx)) {
	SET_BIT(load_mask, indiv_uidx);
      }
    }
  } else {
    memcpy(load_mask, pheno_nm, unfiltered_indiv_ctv2 * (sizeof(intptr_t) / 2));
  }
  return 0;
}

int32_t glm_scan_conditions(char* condition_mname, char* condition_fname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uintptr_t* loadbuf_raw, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_indiv_ct, uintptr_t* sex_male, uintptr_t* load_mask, uintptr_t* indiv_valid_ct_ptr, uintptr_t* condition_ct_ptr, uint32_t** condition_uidxs_ptr, uintptr_t* indiv_raw_include2, uintptr_t* indiv_raw_male_include2) {
  // side effects: load_mask and indiv_valid_ct potentially updated,
  //   condition_ct should be changed, condition_uidxs should be malloc'd
  unsigned char* wkspace_mark = wkspace_base;
  FILE* condition_file = NULL;
  uint32_t* condition_uidxs = NULL;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctv2 = 2 * ((unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t indiv_valid_ct = *indiv_valid_ct_ptr;
  uintptr_t miss_ct = 0;
  uintptr_t condition_ct = 0;
  int32_t retval = 0;
#ifdef __LP64__
  __m128i* loadbuf_vptr;
  __m128i* loadbuf_mask_vptr;
  __m128i* loadbuf_vend;
  __m128i vii;
#else
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t* loadbuf_end;
#endif
  uintptr_t* loadbuf_ptr;
  uintptr_t* loadbuf_mask_ptr;
  char* sorted_ids;
  uint32_t* id_map;
  uintptr_t* already_seen;
  uintptr_t* loadbuf_mask_orig;
  uintptr_t* loadbuf_mask;
  uint32_t* marker_idx_to_uidx;
  uint32_t* condition_uidxs_tmp;
  char* bufptr;
  char* bufptr2;
  uintptr_t condition_ct_max;
  uintptr_t condition_idx;
  uintptr_t indiv_uidx_offset;
  uintptr_t slen;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t marker_uidx;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  int32_t ii;
  char cc;
  if (condition_mname) {
    ii = get_uidx_from_unsorted(condition_mname, marker_exclude, marker_ct, marker_ids, max_marker_id_len);
    if (ii == -1) {
      sprintf(logbuf, "Warning: --condition marker ID '%s' not found.\n", condition_mname);
      logprintb();
      return 0;
    }
    condition_ct = 1;
    condition_uidxs = (uint32_t*)malloc(sizeof(int32_t));
    condition_uidxs[0] = (uint32_t)ii;
  } else {
    if (wkspace_alloc_c_checked(&sorted_ids, marker_ct * max_marker_id_len) ||
        wkspace_alloc_ui_checked(&id_map, marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t)) ||
        wkspace_alloc_ui_checked(&marker_idx_to_uidx, marker_ct * sizeof(int32_t))) {
      goto glm_scan_conditions_ret_NOMEM;
    }
    fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
    fill_ulong_zero(already_seen, marker_ctl);
    retval = sort_item_ids_noalloc(sorted_ids, id_map, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, 1, strcmp_deref);
    if (retval) {
      goto glm_scan_conditions_ret_1;
    }
    condition_uidxs_tmp = (uint32_t*)wkspace_base;
    if (wkspace_left > marker_ct * sizeof(int32_t)) {
      condition_ct_max = marker_ct;
    } else {
      condition_ct_max = wkspace_left / sizeof(int32_t);
    }
    if (fopen_checked(&condition_file, condition_fname, "r")) {
      goto glm_scan_conditions_ret_OPEN_FAIL;
    }
    tbuf[MAXLINELEN - 1] = ' ';
    while (fgets(tbuf, MAXLINELEN, condition_file)) {
      if (!tbuf[MAXLINELEN - 1]) {
        logprint("Error: Pathologically long line in --condition-list file.\n");
        goto glm_scan_conditions_ret_INVALID_FORMAT;
      }
      bufptr = skip_initial_spaces(tbuf);
      while (!is_eoln_kns(*bufptr)) {
        bufptr2 = item_endnn(bufptr);
        slen = (uintptr_t)(bufptr2 - bufptr);
	if (slen < max_marker_id_len) {
	  cc = *bufptr2;
	  *bufptr2 = '\0';
	  ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, marker_ct - 1);
	  if (ii == -1) {
            miss_ct++;
	  }
	  if (is_set(already_seen, ii)) {
	    sprintf(logbuf, "Error: Duplicate marker %s in --condition-list file.\n", bufptr);
            logprintb();
	    goto glm_scan_conditions_ret_INVALID_FORMAT;
	  }
	  if (condition_ct == condition_ct_max) {
	    goto glm_scan_conditions_ret_NOMEM;
	  }
          set_bit(already_seen, ii);
          condition_uidxs_tmp[condition_ct++] = marker_idx_to_uidx[id_map[(uint32_t)ii]];
	  *bufptr2 = cc;
	}
        bufptr = skip_initial_spaces(bufptr2);
      }
    }
    if (!feof(condition_file)) {
      goto glm_scan_conditions_ret_READ_FAIL;
    }
    if (condition_ct) {
      condition_uidxs = (uint32_t*)malloc(condition_ct * sizeof(int32_t));
      memcpy(condition_uidxs, condition_uidxs_tmp, condition_ct * sizeof(int32_t));
    } else if (!miss_ct) {
      logprint("Warning: --condition-list file is empty.\n");
      goto glm_scan_conditions_ret_1;
    }
    if (miss_ct) {
      sprintf(logbuf, "--condition-list: %" PRIuPTR " of %" PRIuPTR " marker ID%s loaded from %s.\n", condition_ct, condition_ct + miss_ct, (condition_ct + miss_ct == 1)? "" : "s", condition_fname);
    } else {
      sprintf(logbuf, "--condition-list: %" PRIuPTR " marker ID%s loaded from %s.\n", condition_ct, (condition_ct == 1)? "" : "s", condition_fname);
    }
    logprintb();
  }
  if (condition_ct) {
    if (wkspace_alloc_ul_checked(&loadbuf_mask_orig, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
        wkspace_alloc_ul_checked(&loadbuf_mask, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
      goto glm_scan_conditions_ret_NOMEM;
    }
    vec_include_init(unfiltered_indiv_ct, loadbuf_mask_orig, load_mask);
    memcpy(loadbuf_mask, loadbuf_mask_orig, unfiltered_indiv_ctv2 * sizeof(intptr_t));
#ifdef __LP64__
    loadbuf_vend = (__m128i*)(&(loadbuf_raw[unfiltered_indiv_ctv2]));
#else
    loadbuf_end = &(loadbuf_raw[unfiltered_indiv_ctl2]);
#endif
    for (condition_idx = 0; condition_idx < condition_ct; condition_idx++) {
      // scan for missing values to finalize indiv_valid_ct
      marker_uidx = condition_uidxs[condition_idx];
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto glm_scan_conditions_ret_READ_FAIL;
      }
      // don't use load_and_collapse since collapse bitmask not finalized
      if (fread(loadbuf_raw, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto glm_scan_conditions_ret_READ_FAIL;
      }
      chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
      if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
	is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
	is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
	haploid_fix(hh_exists, indiv_raw_include2, indiv_raw_male_include2, unfiltered_indiv_ct, is_x, is_y, (unsigned char*)loadbuf_raw);
      }
      // clear loadbuf_mask bits where loadbuf is 01.
#ifdef __LP64__
      loadbuf_vptr = (__m128i*)loadbuf_raw;
      loadbuf_mask_vptr = (__m128i*)loadbuf_mask;
      do {
        vii = *loadbuf_vptr++;
        vii = _mm_andnot_si128(_mm_srli_epi64(vii, 1), vii);
        *loadbuf_mask_vptr = _mm_andnot_si128(vii, *loadbuf_mask_vptr);
        loadbuf_mask_vptr++;
      } while (loadbuf_vptr < loadbuf_vend);
#else
      loadbuf_ptr = loadbuf_raw;
      loadbuf_mask_ptr = loadbuf_mask;
      do {
        ulii = *loadbuf_ptr++;
        ulii = ((~ulii) >> 1) & ulii;
        *loadbuf_mask_ptr &= ~ulii;
	loadbuf_mask_ptr++;
      } while (loadbuf_ptr < loadbuf_end);
#endif
    }
    loadbuf_ptr = loadbuf_mask_orig;
    loadbuf_mask_ptr = loadbuf_mask;
    for (indiv_uidx_offset = 0; indiv_uidx_offset < unfiltered_indiv_ct; indiv_uidx_offset += BITCT2) {
      ulii = (*loadbuf_ptr++) & (~(*loadbuf_mask_ptr++));
      while (ulii) {
        uljj = CTZLU(ulii);
        clear_bit_ul(load_mask, indiv_uidx_offset + (uljj / 2));
        indiv_valid_ct--;
        ulii &= ulii - 1;
      }
    }
    *condition_uidxs_ptr = condition_uidxs;
  }
  *condition_ct_ptr = condition_ct;
  while (0) {
  glm_scan_conditions_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_scan_conditions_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_scan_conditions_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_scan_conditions_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 glm_scan_conditions_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(condition_file);
  return retval;
}

uint32_t glm_loadbuf_to_doubles(uintptr_t* loadbuf_collapsed, uint32_t indiv_valid_ct, double* covar_row, double* geno_map, uintptr_t* cur_missing) {
  // ok for cur_missing to be NULL if there can't possibly be any missing
  uintptr_t* ulptr_end = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uint32_t cur_missing_ct = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_idx_stop = BITCT2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  while (1) {
    while (loadbuf_collapsed < ulptr_end) {
      cur_word = *loadbuf_collapsed++;
      for (; indiv_idx < indiv_idx_stop; indiv_idx++, covar_row++, cur_word >>= 2) {
        cur_genotype = cur_word & 3;
        if (cur_genotype != 1) {
          *covar_row = geno_map[cur_genotype];
	} else {
          SET_BIT(cur_missing, indiv_idx);
	}
      }
      indiv_idx_stop += BITCT2;
    }
    if (indiv_idx == indiv_valid_ct) {
      return cur_missing_ct;
    }
    ulptr_end++;
    indiv_idx_stop = indiv_valid_ct;
  }
}

uint32_t glm_loadbuf_to_doubles_x(uintptr_t* loadbuf_collapsed, uintptr_t* sex_male_collapsed, uint32_t indiv_valid_ct, double* covar_row, double* geno_map, uintptr_t* cur_missing) {
  uintptr_t* ulptr_end = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uint32_t cur_missing_ct = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_idx_stop = BITCT2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  while (1) {
    while (loadbuf_collapsed < ulptr_end) {
      cur_word = *loadbuf_collapsed++;
      for (; indiv_idx < indiv_idx_stop; indiv_idx++, covar_row++, cur_word >>= 2) {
        cur_genotype = cur_word & 3;
        if (cur_genotype != 1) {
          *covar_row = geno_map[cur_genotype + 4 * IS_SET(sex_male_collapsed, indiv_idx)];
	} else {
          SET_BIT(cur_missing, indiv_idx);
	}
      }
      indiv_idx_stop += BITCT2;
    }
    if (indiv_idx == indiv_valid_ct) {
      return cur_missing_ct;
    }
    ulptr_end++;
    indiv_idx_stop = indiv_valid_ct;
  }
}

int32_t glm_check_vif(double vif_thresh, uintptr_t param_ct, uintptr_t indiv_valid_ct, double* covars_collapsed, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* param_2d_buf2) {
  __CLPK_integer dim = ((uint32_t)param_ct) - 1;
  uintptr_t param_ct_m1 = param_ct - 1;
  double indiv_ct_d = (double)((intptr_t)indiv_valid_ct);
  double indiv_ct_recip = 1.0 / indiv_ct_d;
  double indiv_ct_m1_recip = 1.0 / ((double)((intptr_t)(indiv_valid_ct - 1)));
  double* dptr;
  double* dptr2;
  double* dptr3;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  uintptr_t indiv_idx;
  double dxx;
  double dyy;
  int32_t ii;
  for (param_idx = 1; param_idx < param_ct; param_idx++) {
    dyy = 0; // sum
    dptr = &(covars_collapsed[param_idx * indiv_valid_ct]);
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      dxx = *dptr++;
      dyy += dxx;
    }
    param_2d_buf2[param_idx] = dyy * indiv_ct_recip;
  }
  for (param_idx = 1; param_idx < param_ct; param_idx++) {
    dptr = &(param_2d_buf[(param_idx - 1) * param_ct]);
    dyy = param_2d_buf2[param_idx];
    for (param_idx2 = param_idx; param_idx2 < param_ct; param_idx2++) {
      dxx = 0;
      dptr2 = &(covars_collapsed[param_idx * indiv_valid_ct]);
      dptr3 = &(covars_collapsed[param_idx2 * indiv_valid_ct]);
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	dxx += (*dptr2++) * (*dptr3++);
      }
      dxx -= dyy * param_2d_buf2[param_idx2] * indiv_ct_d;
      *dptr++ = dxx * indiv_ct_m1_recip;
    }
  }
  for (param_idx = 0; param_idx < param_ct_m1; param_idx++) {
    dxx = param_2d_buf[param_idx * param_ct]; // diagonal element
    if ((dxx == 0) || (!realnum(dxx))) {
      return 1;
    }
    param_2d_buf2[param_idx] = 1.0 / sqrt(dxx); // now inverse sqrt of diagonal
  }
  for (param_idx = 1; param_idx < param_ct_m1; param_idx++) {
    dyy = param_2d_buf2[param_idx];
    dptr = &(param_2d_buf[param_idx * param_ct_m1]);
    dptr3 = param_2d_buf2;
    for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
      dptr2 = &(param_2d_buf[param_idx2 * param_ct_m1 + param_idx]);
      dxx = (*dptr2) * dyy * (*dptr3++);
      if (dxx > 0.999) {
	return 1;
      }
      *dptr2 = dxx;
      *dptr++ = dxx;
    }
  }
  for (param_idx = 0; param_idx < param_ct_m1; param_idx++) {
    param_2d_buf[param_idx * param_ct] = 1;
  }
  
  ii = invert_matrix(dim, param_2d_buf, mi_buf, param_2d_buf2);
  if (ii) {
    return ii;
  }
  for (param_idx = 0; param_idx < param_ct_m1; param_idx++) {
    if (param_2d_buf[param_idx * param_ct] > vif_thresh) {
      return 1;
    }
  }
  return 0;
}

uint32_t glm_linear_robust_cluster_covar(uintptr_t cur_batch_size, uintptr_t param_ct, uintptr_t indiv_valid_ct, double* covars_cov_major, double* covars_indiv_major, double* pheno_perms, uintptr_t pheno_perms_stride, double* coef, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* param_2d_buf2, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, double* cluster_param_buf, double* cluster_param_buf2, double* indiv_1d_buf, double* linear_results, uint32_t* perm_fail_ct_ptr, uintptr_t* perm_fails) {
  // See the second half of PLINK linear.cpp fitLM(), and validParameters().
  // Diagonals of the final covariance matrices (not including the intercept
  // element) are saved to linear_results[(perm_idx * (param_ct - 1))..
  // ((perm_idx + 1) * (param_ct - 1) - 1)].
  // If not all permutations yield valid results, bits in perm_fails[] are set.
  // A return value of 1 reports that ALL permutations failed.  In this case,
  // perm_fails is not necessarily updated.
  uintptr_t param_ct_p1 = param_ct + 1; // diagonals of param * param matrix
  uintptr_t param_ct_m1 = param_ct - 1;
  uint32_t cluster_ct1_p1 = cluster_ct1 + 1;
  uint32_t perm_fail_ct = 0;
  double* dptr;
  double* dptr2;
  double* dptr3;
  double* pheno_ptr;
  double partial;
  double min_sigma;
  uintptr_t indiv_idx;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  uintptr_t perm_idx;

  // +1 since "no cluster" is represented as -1
  // 32-bit since that -1 is actually 0xffffffffU
  uint32_t cluster_idx_p1;

  double dxx;
  double dyy;
  fill_ulong_zero(perm_fails, ((cur_batch_size + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
  col_major_matrix_multiply((uint32_t)param_ct, (uint32_t)param_ct, (uint32_t)indiv_valid_ct, covars_indiv_major, covars_cov_major, param_2d_buf);
  if (invert_matrix((uint32_t)param_ct, param_2d_buf, mi_buf, param_2d_buf2)) {
    return 1;
  }

  if (!cluster_ct1) {
    // only need to perform S[i][j] / sqrt(S[i][i] * S[j][j]) check once, since
    // it's independent of sigma

    // may as well cache inverse square roots
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      param_2d_buf2[param_idx] = 1.0 / sqrt(param_2d_buf[param_idx * param_ct_p1]);
    }
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      dxx = param_2d_buf2[param_idx]; // 1 / sqrt(S[i][i])
      dptr = &(param_2d_buf[param_idx * param_ct]); // S[i][j]
      dptr2 = param_2d_buf2; // 1 / sqrt(S[j][j])
      for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
	if ((*dptr++) * (*dptr2++) * dxx > 0.99999) {
	  return 1;
	}
      }
    }
    // now determine min(S[i][i]), i >= 1
    min_sigma = param_2d_buf[param_ct_p1];
    for (param_idx = 2; param_idx < param_ct; param_idx++) {
      dxx = param_2d_buf[param_idx * param_ct_p1];
      if (min_sigma > dxx) {
	min_sigma = dxx;
      }
    }
    if (min_sigma <= 0) {
      return 1;
    }
    // S[i][i] * sigma < 1e-20 iff sigma < 1e-20 / S[i][i]
    min_sigma = 1e-20 / min_sigma;

    // now temporarily save sigmas in linear_results[perm_idx * param_ct_m1]
    fill_double_zero(linear_results, cur_batch_size * param_ct_m1);
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      pheno_ptr = &(pheno_perms[indiv_idx * pheno_perms_stride]);
      dptr = &(covars_indiv_major[indiv_idx * param_ct]);
      for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
	partial = 0;
	// not coef[perm_idx * param_ct] due to how dgels() works
	dptr2 = &(coef[perm_idx * indiv_valid_ct]);
	dptr3 = dptr;
	for (param_idx = 0; param_idx < param_ct; param_idx++) {
	  partial += (*dptr2++) * (*dptr3++);
	}
	partial -= *pheno_ptr++;
	linear_results[perm_idx * param_ct_m1] += partial * partial;
      }
    }
    dxx = 1.0 / ((double)((intptr_t)(indiv_valid_ct - param_ct)));
    dptr2 = param_2d_buf2;
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      *dptr2++ = param_2d_buf[param_idx * param_ct_p1]; // S0[i][i]
    }
    dptr = linear_results;
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      dyy = dxx * (*dptr); // sigma = (previous sigma) / nind-np
      if (dyy < min_sigma) {
	dyy = 0;
	perm_fail_ct++;
	dptr = &(dptr[param_ct_m1]);
	SET_BIT(perm_fails, perm_idx);
      }
      dptr2 = param_2d_buf2;
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	*dptr++ = (*dptr2++) * dyy; // S[i][i] = S0[i][i] * sigma
      }
    }
  } else {
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      fill_double_zero(cluster_param_buf, cluster_ct1_p1 * param_ct);
      dptr2 = covars_indiv_major; // X[i][j]
      pheno_ptr = &(pheno_perms[perm_idx]); // Y[i]
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	dptr3 = &(coef[perm_idx * indiv_valid_ct]);
	partial = 0;
	for (param_idx = 0; param_idx < param_ct; param_idx++) {
	  partial += dptr2[param_idx] * (*dptr3++);
	}
	partial -= pheno_ptr[indiv_idx * pheno_perms_stride];
	cluster_idx_p1 = indiv_to_cluster1[indiv_idx] + 1;
	dptr3 = &(cluster_param_buf[cluster_idx_p1 * param_ct]);
        for (param_idx = 0; param_idx < param_ct; param_idx++) {
          // sc[clst[i]][j] += partial * X[i][j]
	  *dptr3 += partial * (*dptr2++);
	  dptr3++;
	}
      }
      dptr = cluster_param_buf2;
      // transpose
      for (param_idx = 0; param_idx < param_ct; param_idx++) {
        dptr2 = &(cluster_param_buf[param_idx]);
        for (cluster_idx_p1 = 0; cluster_idx_p1 < cluster_ct1_p1; cluster_idx_p1++) {
          *dptr++ = dptr2[cluster_idx_p1 * param_ct];
	}
      }
      // can't overwrite param_2d_buf (= S0), everything else is fair game
      col_major_matrix_multiply(param_ct, param_ct, cluster_ct1_p1, cluster_param_buf, cluster_param_buf2, param_2d_buf2);
      col_major_matrix_multiply(param_ct, param_ct, param_ct, param_2d_buf, param_2d_buf2, cluster_param_buf); // multMatrix (S0, meat, tmp1)
      col_major_matrix_multiply(param_ct, param_ct, param_ct, cluster_param_buf, param_2d_buf, param_2d_buf2); // multMatrix (tmp1, S0, S)

      // now do validParameters() check.  validate and cache 1/sqrt(S[i][i])...
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	dxx = param_2d_buf2[param_idx * param_ct_p1];
        if ((dxx < 1e-20) || (!realnum(dxx))) {
	  break;
	}
        cluster_param_buf[param_idx] = 1.0 / sqrt(dxx);
      }
      if (param_idx == param_ct) {
	cluster_param_buf[0] = 1.0 / sqrt(param_2d_buf2[0]);
        for (param_idx = 1; param_idx < param_ct; param_idx++) {
	  dxx = cluster_param_buf[param_idx];
	  dptr = &(param_2d_buf2[param_idx * param_ct]); // S[i][j]
	  dptr2 = cluster_param_buf;
	  for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
	    if ((*dptr++) * (*dptr2++) * dxx > 0.99999) {
	      goto glm_linear_robust_cluster_covar_multicollinear;
	    }
	  }
	}
	dptr = &(linear_results[perm_idx * param_ct_m1]);
        for (param_idx = 1; param_idx < param_ct; param_idx++) {
	  *dptr++ = param_2d_buf2[param_idx * param_ct_p1];
	}
      } else {
      glm_linear_robust_cluster_covar_multicollinear:
	// technically may not need to fill with zeroes
        fill_double_zero(&(linear_results[perm_idx * param_ct_m1]), param_ct_m1);
	SET_BIT(perm_fails, perm_idx);
	perm_fail_ct++;
      }
    }
  }
  *perm_fail_ct_ptr = perm_fail_ct;
  return 0;
}

// make this configurable on command line?
#define LOGISTIC_MAX_ITERS 20

/*
uint32_t glm_logistic_robust_cluster_covar(uintptr_t cur_batch_size, uintptr_t param_ct, uintptr_t indiv_valid_ct, double* covars_cov_major, double* covars_indiv_major, double* pheno_perms, uintptr_t pheno_perms_stride, double* coef, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* param_2d_buf2, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, double* cluster_param_buf, double* cluster_param_buf2, double* indiv_1d_buf, double* linear_results, uint32_t* perm_fail_ct_ptr, uintptr_t* perm_fails) {
*/
uint32_t glm_logistic_robust_cluster_covar(uintptr_t cur_batch_size, uintptr_t param_ct, uintptr_t indiv_valid_ct, double* covars_cov_major, double* pheno_perms, uintptr_t pheno_perms_stride, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, double* logistic_results, uint32_t* perm_fail_ct_ptr, uintptr_t* perm_fails) {
  // See PLINK logistic.cpp fitLM().
  // todo
  uint32_t iters = 0;
  uintptr_t param_idx;
  double delta;
  do {
    delta = 0;
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      // delta += fabs();
    }
    if (delta < 0.000001) {
      break;
    }
  } while (++iters < LOGISTIC_MAX_ITERS);
  return 0;
}

/*
int32_t glm_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctv2 = 2 * ((unfiltered_indiv_ct + BITCT - 1) / BITCT);
  FILE* outfile = NULL;
  FILE* outfile_msa = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t indiv_valid_ct = 0;
  uint32_t perm_adapt = glm_modifier & GLM_PERM;
  uint32_t perm_maxt = glm_modifier & GLM_MPERM;
  uint32_t do_perms = glm_modifier & (GLM_PERM | GLM_MPERM);
  uint32_t perm_count = glm_modifier & GLM_PERM_COUNT;
  uint32_t standard_beta = glm_modifier & GLM_STANDARD_BETA;
  uint32_t report_odds = pheno_c && (!(glm_modifier & GLM_BETA));
  uint32_t display_ci = (ci_size > 0);
  uint32_t fill_orig_chiabs = do_perms || mtest_adjust;
  uint32_t perms_total = 0;
  uint32_t perm_pass_idx = 0;
  uint32_t np_base = 2;
  uint32_t variation_in_sex = 0; // no need to initialize if no-x-sex specified
  uint32_t diploid_present = 0;
  uint32_t diploid_np = 0;
  uint32_t x_present = 0;
  uint32_t sex_always = 0;
  uint32_t sex_sometimes = 0;
  uint32_t sex_np = 0;
  uint32_t male_ct = 0;
  int32_t retval = 0;
  char* outname_end2;
  uintptr_t indiv_idx;
  int32_t ii;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    if (IS_SET(pheno_nm, indiv_uidx) & IS_SET(covar_nm, indiv_idx)) {
      if (is_set(sex_male, indiv_uidx)) {
	male_ct++;
      }
      indiv_valid_ct++;
    }
    indiv_uidx++;
  }
  ii = chrom_info_ptr->x_code;
  if ((ii != -1) && count_chrom_markers(chrom_info_ptr, (uint32_t)ii, marker_exclude)) {
    x_present = 1;
  }
  if (!(glm_modifier & GLM_NO_X_SEX)) {
    variation_in_sex = ((!male_ct) || (male_ct == indiv_valid_ct))? 0 : 1;
    if (variation_in_sex) {
      sex_always = (glm_modifier & GLM_SEX)? 1 : 0;
      sex_sometimes = (sex_always || x_present)? 1 : 0;
      sex_np = ;
    }
  }
  
  if (indiv_valid_ct < np_base + diploid_np + sex_np) {
    logprint("Warning: Skipping --linear/--logistic since number of variables > number of samples.\n");
    if (pheno_nm_ct >= np_base + diploid_np + sex_np) {
      logprint("(Check your covariate file--all samples with at least one missing covariate are\nexcluded from this analysis.)\n");
    }
    goto glm_assoc_ret_1;
  }
  g_mperm_save_all = NULL;
  if (perm_maxt) {
    perms_total = glm_mperm_val;
    if (wkspace_alloc_d_checked(&g_maxt_extreme_stat, perms_total * sizeof(double))) {
      goto glm_assoc_ret_NOMEM;
    }
    fill_double_zero(g_maxt_extreme_stat, perms_total);
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(&outfile_msa, outname, "w")) {
        goto glm_assoc_ret_OPEN_FAIL;
      }
      if (putc_checked('0', outfile_msa)) {
        goto glm_assoc_ret_WRITE_FAIL;
      }
      sprintf(logbuf, "Dumping all permutation squared t-stats to %s.\n", outname);
      logprintb();
    }
  } else {
    mperm_save = 0;
    if (perm_adapt) {
    }
  }
  if (pheno_d) {
    outname_end2 = memcpyb(outname_end, ".assoc.linear", 14);
  } else {
    outname_end2 = memcpyb(outname_end, ".assoc.logistic", 16);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto glm_assoc_ret_OPEN_FAIL;
  }
  sprintf(logbuf, "Writing %s model association results to %s...", pheno_d? "linear" : "logistic", outname);
  logprintb();
  fflush(stdout);
  sprintf(tbuf, " CHR %%%us         BP   A1       TEST    NMISS       %s ", plink_maxsnp, report_odds? "  OR" : "BETA");
  fprintf(outfile, tbuf, "SNP");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "      SE      L%u      U%u ", uii, uii);
    } else {
      fprintf(outfile, "      SE       L%u       U%u ", uii, uii);
    }
  }
  fputs("        STAT            P \n", outfile);
  logprint("\nError: --linear and --logistic are currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  while (0) {
  glm_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 glm_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
}
*/

int32_t glm_assoc_nosnp(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double ci_zt, double pfilter, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctv2 = 2 * ((unfiltered_indiv_ct + BITCT - 1) / BITCT);
  FILE* outfile = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t topsize = 0;
  uintptr_t max_param_name_len = 2;
  uintptr_t perm_vec_ctcl8m = 0;
  uintptr_t param_raw_ct = 1;
  uintptr_t condition_ct = 0;
  uint32_t cur_batch_size = 1;
  uint32_t cluster_ct1 = 0;
  uint32_t do_perms = (glm_modifier & GLM_MPERM)? 1 : 0;
  uint32_t perm_count = glm_modifier & GLM_PERM_COUNT;
  uint32_t report_odds = pheno_c && (!(glm_modifier & GLM_BETA));
  uint32_t display_ci = (ci_size > 0);
  // uint32_t perm_pass_idx = 0;
  uint32_t variation_in_sex = 0; // not initialized without sex modifier
  uint32_t male_ct = 0;
  // uint32_t perm_attempt_ct = 0;
  uint32_t perms_done = 0;
  int32_t retval = 0;
  uint32_t perm_fail_total = 0;
#ifndef NOLAPACK
  uint32_t perm_fail_ct = 0;
  char dgels_trans = 'N';
  int32_t dgels_m = 0;
  int32_t dgels_n = 0;
  int32_t dgels_nrhs = 0;
  double* dgels_a = NULL;
  double* dgels_b = NULL;
  double* dgels_work = NULL;
  int32_t dgels_ldb = 0;
  int32_t dgels_lwork = -1;
  int32_t dgels_info;
  double* linear_results = NULL;
#endif
  double* cluster_param_buf = NULL; // guaranteed to be size >= param_ct^2
  double* cluster_param_buf2 = NULL; // not guaranteed
  double* indiv_1d_buf = NULL; // currently needed only with clusters
  double* mperm_save_stats = NULL;
  uintptr_t* loadbuf_raw = NULL;
  uintptr_t* loadbuf_collapsed = NULL;
  uintptr_t* load_mask = NULL;
  uintptr_t* sex_male_collapsed = NULL;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* active_params = NULL;
  uintptr_t* perm_fails = NULL;
  uint32_t* perm_2success_ct = NULL;
  uint32_t* condition_uidxs = NULL;
  uint32_t* cluster_map1 = NULL;
  uint32_t* cluster_starts1 = NULL;
  uint32_t* indiv_to_cluster1 = NULL;
  double geno_map[12];
  uint32_t uibuf[4];
  double* geno_map_ptr;
  double* param_2d_buf;
  double* param_2d_buf2;
  MATRIX_INVERT_BUF1_TYPE* mi_buf;
  char* param_names;
  double* covars_collapsed; // covariate-major
  double* covars_transposed_buf; // individual-major

  // linear: pval = calc_tprob(orig_stats[i], indiv_valid_ct - param_ct)
  // logistic: pval = chiprob_p(orig_stats[i] * orig_stats[i], 1);
  double* orig_stats;

  char* outname_end2;
  char* wptr;
  char* wptr2;
  double* dptr;
  double* dptr2;
  uintptr_t indiv_valid_ct;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t param_raw_ctl;
  uintptr_t param_ct;
  uintptr_t param_idx;
  uintptr_t ulii;
  // warning-suppressing stopgap
#ifndef NOLAPACK
  double* msa_ptr;
  double* dptr3;
  double se;
  double zval;
#endif
  double pval;
  double dxx;
  double dyy;
  double dzz;
  uint32_t perm_idx;
  uint32_t marker_uidx;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t uii;
  uint32_t ujj;
  uint32_t slen;
  int32_t ii;
  if (glm_init_load_mask(indiv_exclude, pheno_nm, covar_nm, indiv_ct, unfiltered_indiv_ctv2, &load_mask)) {
    goto glm_assoc_nosnp_ret_NOMEM;
  }
  indiv_valid_ct = popcount_longs(load_mask, 0, unfiltered_indiv_ctv2 / 2);
  if (condition_mname || condition_fname) {
    loadbuf_raw = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
    if (!loadbuf_raw) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    loadbuf_raw[unfiltered_indiv_ctv2 - 2] = 0;
    loadbuf_raw[unfiltered_indiv_ctv2 - 1] = 0;
    ulii = topsize;

    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_include2) {
        goto glm_assoc_nosnp_ret_NOMEM;
      }
      fill_vec_55(indiv_include2, unfiltered_indiv_ct); // harmless
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_male_include2) {
	goto glm_assoc_nosnp_ret_NOMEM;
      }
      fill_ulong_zero(indiv_male_include2, unfiltered_indiv_ctv2);
      vec_include_init(unfiltered_indiv_ct, indiv_male_include2, sex_male);
    }
    wkspace_left -= topsize;
    retval = glm_scan_conditions(condition_mname, condition_fname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, chrom_info_ptr, hh_exists, loadbuf_raw, bedfile, bed_offset, unfiltered_indiv_ct, sex_male, load_mask, &indiv_valid_ct, &condition_ct, &condition_uidxs, indiv_include2, indiv_male_include2);
    wkspace_left += topsize;
    if (retval) {
      goto glm_assoc_nosnp_ret_1;
    }

    // now initialize collapsed indiv_include2, indiv_male_include2, sex_male
    topsize = ulii;
    ulii = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);
    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, ulii * sizeof(intptr_t));
      fill_vec_55(indiv_include2, indiv_valid_ct);
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, ulii * sizeof(intptr_t));
      alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_valid_ct, hh_exists, 1, load_mask, sex_male, &indiv_include2, &indiv_male_include2);
    }
    loadbuf_collapsed = (uintptr_t*)top_alloc(&topsize, ulii * sizeof(intptr_t));
    if (!loadbuf_collapsed) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    loadbuf_collapsed[ulii - 1] = 0;
    sex_male_collapsed = (uintptr_t*)top_alloc(&topsize, ulii * (sizeof(intptr_t) / 2));
    if (!sex_male_collapsed) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    collapse_copy_bitarr_incl(unfiltered_indiv_ct, sex_male, load_mask, indiv_valid_ct, sex_male_collapsed);

    param_raw_ct += condition_ct;
  }
  param_raw_ct += covar_ct;
  if (glm_modifier & GLM_SEX) {
    indiv_uidx = 0;
    indiv_idx = 0;
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      do {
	if (IS_SET(sex_male, indiv_uidx)) {
	  male_ct++;
	}
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_valid_ct);
    variation_in_sex = ((!male_ct) || (male_ct == indiv_valid_ct))? 0 : 1;
    if (variation_in_sex) {
      param_raw_ct++;
    }
  }
  param_raw_ctl = (param_raw_ct + BITCT - 1) / BITCT;
  active_params = (uintptr_t*)malloc(param_raw_ctl * sizeof(intptr_t));
  if (!active_params) {
    goto glm_assoc_nosnp_ret_NOMEM;
  }
  if (parameters_range_list_ptr->name_ct) {
    fill_ulong_zero(active_params, param_raw_ctl);
    active_params[0] = 1;
    numeric_range_list_to_bitfield(parameters_range_list_ptr, param_raw_ct, active_params, 0, 1);
    param_ct = popcount_longs(active_params, 0, param_raw_ctl);
  } else {
    active_params[param_raw_ctl - 1] = 0;
    fill_bits(active_params, 0, param_raw_ct);
    param_ct = param_raw_ct;
  }
  if (param_ct == 1) {
    logprint("Warning: Skipping --linear/--logistic since the intercept is the only variable.\n");
    goto glm_assoc_nosnp_ret_1;
  } else if (indiv_valid_ct <= param_ct) {
    logprint("Warning: Skipping --linear/--logistic since # variables >= # samples.\n");
    if (pheno_nm_ct > param_ct) {
      logprint("(Check your covariates--all samples with at least one missing covariate are\nexcluded from this analysis.)\n");
    }
    goto glm_assoc_nosnp_ret_1;
  }
  // parameter sequence:
  // 1. intercept
  // 2. --condition-list (1 to condition_ct in raw list)
  // 3. --covar (condition_ct + 1 to condition_ct + covar_ct in raw list)
  // 4. sex (condition_ct + covar_ct + 1 in raw list)
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + 1)) {
      slen = strlen(&(marker_ids[condition_uidxs[uii] * max_marker_id_len]));
      if (max_param_name_len <= slen) {
        max_param_name_len = slen + 1;
      }
    }
  }
  uii = condition_ct + 1;
  for (ujj = 0; ujj < covar_ct; ujj++) {
    if (is_set(active_params, ujj + uii)) {
      slen = strlen(&(covar_names[ujj * max_covar_name_len]));
      if (max_param_name_len <= slen) {
	max_param_name_len = slen + 1;
      }
    }
  }

  ulii = condition_ct + covar_ct + 1;
  if (variation_in_sex && IS_SET(active_params, ulii)) {
    if (max_param_name_len < 4) {
      max_param_name_len = 4;
    }
  }
  wkspace_left -= topsize;
  if (wkspace_alloc_c_checked(&param_names, param_ct * max_param_name_len) ||
      wkspace_alloc_d_checked(&covars_collapsed, param_ct * indiv_valid_ct * sizeof(double))) {
    goto glm_assoc_nosnp_ret_NOMEM2;
  }
  wkspace_left += topsize;
  for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
    covars_collapsed[indiv_idx] = 1;
  }
  param_idx = 1;
  // 0..3: diploid chromosomes, X chromosome female
  // 4..7: X chromosome male
  // 8..11: haploid
  fill_double_zero(geno_map, 12);
  geno_map[2] = 1;
  geno_map[3] = 1;
  geno_map[7] = 1;
  geno_map[11] = 1;
  if (glm_modifier & GLM_CONDITION_RECESSIVE) {
    geno_map[2] = 0;
  } else if (!(glm_modifier & GLM_CONDITION_DOMINANT)) {
    geno_map[3] = 2;
    if (glm_xchr_model == 2) {
      geno_map[7] = 2;
    }
  }
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + 1)) {
      marker_uidx = condition_uidxs[uii];
      chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto glm_assoc_nosnp_ret_READ_FAIL;
      }
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_collapsed, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
	goto glm_assoc_nosnp_ret_READ_FAIL;
      }
      chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
      geno_map_ptr = geno_map;
      if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
	is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
	is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
	if (hh_exists) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	if (!is_x) {
	  geno_map_ptr = &(geno_map[8]);
	}
      } else {
	is_x = 0;
      }
      if (!is_x) {
        glm_loadbuf_to_doubles(loadbuf_collapsed, indiv_valid_ct, &(covars_collapsed[param_idx * indiv_valid_ct]), geno_map_ptr, NULL);
      } else {
	glm_loadbuf_to_doubles_x(loadbuf_collapsed, sex_male_collapsed, indiv_valid_ct, &(covars_collapsed[param_idx * indiv_valid_ct]), geno_map_ptr, NULL);
      }
      strcpy(&(param_names[param_idx * max_param_name_len]), &(marker_ids[marker_uidx * max_marker_id_len]));
      param_idx++;
    }
  }
  // topsize = 0;
  if (wkspace_alloc_d_checked(&covars_transposed_buf, param_ct * indiv_valid_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&g_pheno_d2, indiv_valid_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&perm_2success_ct, (param_ct - 1) * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&orig_stats, (param_ct - 1) * sizeof(double)) ||
      wkspace_alloc_d_checked(&param_2d_buf, param_ct * param_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&param_2d_buf2, param_ct * param_ct * sizeof(double))) {
    goto glm_assoc_nosnp_ret_NOMEM;
  }
  mi_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(param_ct * sizeof(MATRIX_INVERT_BUF1_TYPE));
  if (!mi_buf) {
    goto glm_assoc_nosnp_ret_NOMEM;
  }
  fill_uint_zero(perm_2success_ct, param_ct - 1);
  g_pheno_sum = 0;
  g_pheno_ssq = 0;
  if (pheno_d) {
#ifdef NOLAPACK
    logprint("Warning: Skipping --logistic on --all-pheno QT since this is a no-LAPACK " PROG_NAME_CAPS "\nbuild.\n");
    goto glm_assoc_nosnp_ret_1;
#else
    indiv_uidx = 0;
    indiv_idx = 0;
    dptr = g_pheno_d2;
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      dptr2 = &(pheno_d[indiv_uidx]);
      indiv_uidx = indiv_uidx_stop;
      dptr3 = &(pheno_d[indiv_uidx_stop]);
      do {
	dxx = *dptr2++;
        *dptr++ = dxx;
	g_pheno_sum += dxx;
	g_pheno_ssq += dxx * dxx;
      } while (dptr2 < dptr3);
    } while (indiv_idx < indiv_valid_ct);
    if (g_pheno_ssq == g_pheno_sum * (g_pheno_sum / ((double)((intptr_t)indiv_valid_ct)))) {
      logprint("Warning: Skipping --linear/--logistic since phenotype is constant.\n");
      goto glm_assoc_nosnp_ret_1;
    }
#endif
  }
  uii = condition_ct + 1;
  for (ujj = 0; ujj < covar_ct; ujj++) {
    if (is_set(active_params, ujj + uii)) {
      // indiv-major to covariate-major
      indiv_uidx = 0;
      dptr = &(covars_collapsed[param_idx * indiv_valid_ct]);
      dptr2 = &(covar_d[ujj]);
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	if (IS_SET(load_mask, indiv_uidx)) {
	  *dptr++ = dptr2[indiv_idx * covar_ct];
	}
	indiv_uidx++;
      }
      strcpy(&(param_names[param_idx * max_param_name_len]), &(covar_names[ujj * max_covar_name_len]));
      param_idx++;
    }
  }
  ulii = condition_ct + covar_ct + 1;;
  if (variation_in_sex && IS_SET(active_params, ulii)) {
    indiv_uidx = 0;
    dptr = &(covars_collapsed[param_idx * indiv_valid_ct]);
    fill_double_zero(dptr, indiv_valid_ct);
    dptr2 = &(dptr[indiv_valid_ct]);
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      do {
        if (IS_SET(sex_male, indiv_uidx)) {
	  *dptr = 1;
	}
	dptr++;
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (dptr < dptr2);
    memcpy(&(param_names[param_idx * max_param_name_len]), "SEX", 4);
  }
  if (glm_modifier & GLM_STANDARD_BETA) {
    // with no SNPs, only need to do this once
    dzz = g_pheno_sum / ((double)((intptr_t)indiv_valid_ct)); // mean
    dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (g_pheno_ssq - g_pheno_sum * dzz)); // 1/stdev
    dptr = g_pheno_d2;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      *dptr = ((*dptr) - dzz) * dyy;
      dptr++;
    }
    dxx = 0; // sum
    dyy = 0; // ssq
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      dptr = &(covars_collapsed[param_idx * indiv_valid_ct]);
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
        dzz = *dptr++;
        dxx += dzz;
        dyy += dzz * dzz;
      }
      dptr = &(covars_collapsed[param_idx * indiv_valid_ct]);
      dzz = dxx / ((double)((intptr_t)indiv_valid_ct));
      dyy = sqrt((dyy - dxx * dzz) / ((double)((intptr_t)(indiv_valid_ct - 1))));
      if (dyy == 0) {
	fill_double_zero(dptr, indiv_valid_ct);
      } else {
	dyy = 1.0 / dyy;
        for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
          *dptr = ((*dptr) - dzz) * dyy;
          dptr++;
        }
      }
    }
  }
  ii = glm_check_vif(glm_vif_thresh, param_ct, indiv_valid_ct, covars_collapsed, param_2d_buf, mi_buf, param_2d_buf2);
  if (ii == -1) {
    goto glm_assoc_nosnp_ret_NOMEM;
  } else if (ii == 1) {
    logprint("Warning: Skipping --linear/--logistic no-snp since VIF check failed.\n");
    goto glm_assoc_nosnp_ret_1;
  }

  // required for multithreaded permutation generation
  g_cluster_ct = 0;
  g_pheno_nm_ct = indiv_valid_ct;

  perms_done = 0;

  if (do_perms) {
    if (!perm_batch_size) {
      // er, maybe this should be initialized in main()
      perm_batch_size = 512;
    }
    // not actually max(T), just fixed permutation count.
    if (perm_batch_size > glm_mperm_val) {
      perm_batch_size = glm_mperm_val;
    }
  } else {
    perm_batch_size = 1;
    mperm_save = 0;
  }
  ulii = (perm_batch_size + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(&perm_fails, ulii * sizeof(intptr_t))) {
    goto glm_assoc_nosnp_ret_NOMEM;
  }

  if (pheno_d) {
#ifndef NOLAPACK
    // may want to put a multiple linear regression wrapper function in
    // wdist_dmatrix, perhaps with the PLINK 1.07 svdcmp/svbksb no-LAPACK
    // fallback

    // multiple linear regression-specific allocations and preprocessing
    dgels_m = (int32_t)((uint32_t)indiv_valid_ct);
    dgels_n = (int32_t)((uint32_t)param_ct);
    dgels_nrhs = perm_batch_size;
    dgels_ldb = dgels_m;
    if (wkspace_alloc_d_checked(&dgels_a, param_ct * indiv_valid_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&dgels_b, perm_batch_size * indiv_valid_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&linear_results, perm_batch_size * (param_ct - 1) * sizeof(double)) ||
        wkspace_alloc_d_checked(&dgels_work, sizeof(double))) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    fill_double_zero(linear_results, param_ct - 1);
    memcpy(dgels_a, covars_collapsed, param_ct * indiv_valid_ct * sizeof(double));
    memcpy(dgels_b, g_pheno_d2, indiv_valid_ct * sizeof(double));
    dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
    if (dgels_work[0] > 2147483647.0) {
      // maybe this can't actually happen, but just in case...
      // (todo: see if there is any way to do cross-platform linking of a
      // 64-bit integer LAPACK in the near future)
      logprint("Error: Multiple linear regression problem too large for current LAPACK version.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto glm_assoc_nosnp_ret_1;
    }
    dgels_lwork = (int32_t)dgels_work[0];
    wkspace_reset((unsigned char*)dgels_work);
    if (wkspace_alloc_d_checked(&dgels_work, dgels_lwork * sizeof(double))) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    dgels_nrhs = 1;
    dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
    if (dgels_info) {
      logprint("Warning: Skipping --linear/--logistic no-snp since regression failed.\n");
      goto glm_assoc_nosnp_ret_1;
    }
#endif
  } else {
    // todo: logistic regression-specific allocations and preprocessing
    logprint("\nError: --linear and --logistic are currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto glm_assoc_nosnp_ret_1;
  }

  if (cluster_starts) {
    // If there are any size-1 clusters, we actually want two cluster indexes:
    // - one for permutation, which excludes the size-1 clusters, and
    // - one for use by the robust cluster variance estimator, which does not.
    retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 0, NULL, indiv_valid_ct, cluster_ct, cluster_map, cluster_starts, &cluster_ct1, &cluster_map1, &cluster_starts1, NULL, NULL);
    if (retval) {
      goto glm_assoc_nosnp_ret_1;
    }
    if (cluster_ct1) {
      ulii = MAXV(cluster_ct1 + 1, param_ct);
      if (wkspace_alloc_d_checked(&cluster_param_buf, ulii * param_ct * sizeof(double)) ||
          wkspace_alloc_d_checked(&cluster_param_buf2, (cluster_ct1 + 1) * param_ct * sizeof(double)) ||
          wkspace_alloc_d_checked(&indiv_1d_buf, indiv_valid_ct * sizeof(double)) ||
          wkspace_alloc_ui_checked(&indiv_to_cluster1, indiv_valid_ct * sizeof(int32_t))) {
	goto glm_assoc_nosnp_ret_NOMEM;
      }
      fill_unfiltered_indiv_to_cluster(indiv_valid_ct, cluster_ct1, cluster_map1, cluster_starts1, indiv_to_cluster1);
      if (do_perms) {
	retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 1, NULL, indiv_valid_ct, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, NULL, NULL);
	if (retval) {
	  goto glm_assoc_nosnp_ret_1;
	}
	if (!g_cluster_ct) {
	  logprint("Error: No size 2+ clusters for permutation test.\n");
	  goto glm_assoc_nosnp_ret_INVALID_CMDLINE;
	}
	if (wkspace_alloc_ui_checked(&g_indiv_to_cluster, indiv_valid_ct * sizeof(int32_t)) ||
	    wkspace_alloc_ui_checked(&g_qassoc_cluster_thread_wkspace, g_thread_ct * ((g_cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32) * CACHELINE)) {
	  goto glm_assoc_nosnp_ret_NOMEM;
	}
	fill_unfiltered_indiv_to_cluster(indiv_valid_ct, g_cluster_ct, g_cluster_map, g_cluster_starts, g_indiv_to_cluster);
      }
    }
  }
  if (do_perms) {
    g_sfmtp_arr = (sfmt_t**)wkspace_alloc(g_thread_ct * sizeof(intptr_t));
    if (!g_sfmtp_arr) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    g_sfmtp_arr[0] = &sfmt;
    if (g_thread_ct > 1) {
      for (uii = 1; uii < g_thread_ct; uii++) {
	g_sfmtp_arr[uii] = (sfmt_t*)wkspace_alloc(sizeof(sfmt_t));
	if (!g_sfmtp_arr[uii]) {
	  goto glm_assoc_nosnp_ret_NOMEM;
	}
	for (ujj = 0; ujj < 4; ujj++) {
	  uibuf[ujj] = sfmt_genrand_uint32(&sfmt);
	}
	sfmt_init_by_array(g_sfmtp_arr[uii], uibuf, 4);
      }
    }
  }

#ifndef NOLAPACK
  if (pheno_d) {
    dptr = covars_transposed_buf;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      dptr2 = &(covars_collapsed[indiv_idx]);
      for (param_idx = 0; param_idx < param_ct; param_idx++) {
        *dptr++ = dptr2[param_idx * indiv_valid_ct];
      }
    }
    if (glm_linear_robust_cluster_covar(1, param_ct, indiv_valid_ct, covars_collapsed, covars_transposed_buf, g_pheno_d2, 1, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, linear_results, &perm_fail_ct, perm_fails) || perm_fail_ct) {
      logprint("Warning: Skipping --linear/--logistic no-snp due to multicollinearity.\n");
      goto glm_assoc_nosnp_ret_1;
    }
  }
#endif

  if (do_perms) {
    // Note that, for now, the main nosnp regression loop is not multithreaded;
    // only the permutation generation process is.
    perm_vec_ctcl8m = (perm_batch_size + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
    if (wkspace_alloc_d_checked(&g_perm_vecstd, perm_vec_ctcl8m * sizeof(double) * indiv_valid_ct)) {
      goto glm_assoc_nosnp_ret_NOMEM;
    }
    if (mperm_save) {
      // --mperm-save prevented during command-line parsing, so must be
      // --mperm-save-all
      if (wkspace_alloc_d_checked(&mperm_save_stats, glm_mperm_val * (param_ct - 1) * sizeof(double))) {
	goto glm_assoc_nosnp_ret_NOMEM;
      }
      *outname_end = '\0';
      sprintf(logbuf, "Dumping all permutation absolute t-stats to %s.[test ID].mperm.dump.all.\n", outname);
      logprintb();
    }
  }
  if (pheno_d) {
    outname_end2 = memcpyb(outname_end, ".assoc.linear", 14);
  } else {
    outname_end2 = memcpyb(outname_end, ".assoc.logistic", 16);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto glm_assoc_nosnp_ret_OPEN_FAIL;
  }
  sprintf(logbuf, "Writing %s model association results to %s...", pheno_d? "linear" : "logistic", outname);
  logprintb();
  fflush(stdout);
  fprintf(outfile, "      TEST    NMISS       %s ", report_odds? "  OR" : "BETA");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "      SE      L%u      U%u ", uii, uii);
    } else {
      fprintf(outfile, "      SE       L%u       U%u ", uii, uii);
    }
  }
  if (fputs_checked("        STAT            P \n", outfile)) {
    goto glm_assoc_nosnp_ret_WRITE_FAIL;
  }

  if (pheno_d) {
#ifndef NOLAPACK
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      dxx = dgels_b[param_idx]; // coef[p]
      se = sqrt(linear_results[param_idx - 1]);
      zval = dxx / se;
      orig_stats[param_idx - 1] = fabs(zval);
      pval = calc_tprob(zval, indiv_valid_ct - param_ct);
      if (pval <= pfilter) {
	wptr = fw_strcpy(10, &(param_names[param_idx * max_param_name_len]), tbuf);
	*wptr++ = ' ';
        wptr = uint32_writew8x(wptr, (uint32_t)indiv_valid_ct, ' ');
        wptr = double_g_writewx4x(wptr, dxx, 10, ' ');
	if (display_ci) {
	  dyy = ci_zt * se;
	  wptr = double_g_writewx4x(wptr, se, 8, ' ');
	  wptr = double_g_writewx4x(wptr, se - dyy, 8, ' ');
	  wptr = double_g_writewx4x(wptr, se + dyy, 8, ' ');
	}
        wptr = double_g_writewx4x(wptr, zval, 12, ' ');
	wptr = double_g_writewx4x(wptr, pval, 12, '\n');
	if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto glm_assoc_nosnp_ret_WRITE_FAIL;
	}
      }      
    }
#endif
  } else {
  }
  if (fclose_null(&outfile)) {
    goto glm_assoc_nosnp_ret_WRITE_FAIL;
  }
  logprint(" done.\n");

#ifndef NOLAPACK
  msa_ptr = mperm_save_stats;
#endif
  while (perms_done < glm_mperm_val) {
    cur_batch_size = perm_batch_size;
    if (cur_batch_size > glm_mperm_val - perms_done) {
      cur_batch_size = glm_mperm_val - perms_done;
    }
    g_perm_vec_ct = cur_batch_size;
    perm_vec_ctcl8m = (cur_batch_size + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
    if (cur_batch_size >= CACHELINE_DBL * g_thread_ct) {
      g_assoc_thread_ct = g_thread_ct;
    } else {
      g_assoc_thread_ct = cur_batch_size / CACHELINE_DBL;
      if (!g_assoc_thread_ct) {
        g_assoc_thread_ct = 1;
      }
    }
    ulii = 0;
    if (!g_cluster_ct) {
      if (spawn_threads(threads, &qassoc_gen_perms_thread, g_assoc_thread_ct)) {
	goto glm_assoc_nosnp_ret_THREAD_CREATE_FAIL;
      }
      qassoc_gen_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &qassoc_gen_cluster_perms_thread, g_assoc_thread_ct)) {
        goto glm_assoc_nosnp_ret_THREAD_CREATE_FAIL;
      }
      qassoc_gen_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_assoc_thread_ct);

    if (pheno_d) {
#ifndef NOLAPACK
      dgels_nrhs = cur_batch_size;
      fill_double_zero(linear_results, (param_ct - 1) * cur_batch_size);
      memcpy(dgels_a, covars_collapsed, param_ct * indiv_valid_ct * sizeof(double));
      // bleah, have to transpose
      dptr = dgels_b;
      for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
        dptr2 = &(g_perm_vecstd[perm_idx]);
        for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	  *dptr++ = dptr2[indiv_idx * perm_vec_ctcl8m];
	}
      }
      dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
      if (glm_linear_robust_cluster_covar(cur_batch_size, param_ct, indiv_valid_ct, covars_collapsed, covars_transposed_buf, g_perm_vecstd, perm_vec_ctcl8m, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, linear_results, &perm_fail_ct, perm_fails)) {
	perm_fail_ct = cur_batch_size;
      }
      perm_fail_total += perm_fail_ct;
      ulii = param_ct - 1;
      for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
	if (IS_SET(perm_fails, perm_idx)) {
	  if (mperm_save) {
	    for (param_idx = 0; param_idx < ulii; param_idx++) {
	      *msa_ptr++ = -9;
	    }
	  }
	  continue;
	}
        // permutation-major regression coefficients
        dptr = &(dgels_b[perm_idx * indiv_valid_ct + 1]);
	// permutation-major variances
	dptr2 = &(linear_results[perm_idx * ulii]);
	dptr3 = orig_stats;
        for (param_idx = 0; param_idx < ulii; param_idx++) {
	  dxx = *dptr++;
	  se = sqrt(*dptr2++);
          zval = fabs(dxx / se);
	  dyy = *dptr3++;
	  if (zval > dyy + EPSILON) {
	    perm_2success_ct[param_idx] += 2;
	  } else if (zval > dyy - EPSILON) {
	    perm_2success_ct[param_idx] += 1;
	  }
	  if (mperm_save) {
	    *msa_ptr++ = zval;
	  }
	}
      }
#endif
    } else {
      // todo: logistic per-permutation-block stuff
    }
    perms_done += cur_batch_size;
    putchar('\r');
    sprintf(logbuf, "%u permutation%s complete.", perms_done, (perms_done > 1)? "s" : "");
    logprintb();
    fflush(stdout);
  }
  if (do_perms) {
    putchar('\n');
    memcpy(outname_end2, ".mperm", 7);
    if (fopen_checked(&outfile, outname, "w")) {
      goto glm_assoc_nosnp_ret_OPEN_FAIL;
    }
    if (fputs_checked("      TEST         EMP1           NP \n", outfile)) {
      goto glm_assoc_nosnp_ret_WRITE_FAIL;
    }
    dxx = 0.5 / ((double)((int32_t)(glm_mperm_val - perm_fail_total) + 1));
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      wptr = fw_strcpy(10, &(param_names[param_idx * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      pval = ((double)(perm_2success_ct[param_idx - 1] + 2)) * dxx;
      if (pval <= pfilter) {
	if (!perm_count) {
	  wptr = double_g_writewx4(wptr, pval, 12);
	} else {
          wptr = double_g_writewx4(wptr, ((double)perm_2success_ct[param_idx - 1]) / 2.0, 12);
	}
        wptr = memseta(wptr, 32, 3);
        wptr = uint32_writew10(wptr, glm_mperm_val - perm_fail_total);
        wptr = memcpya(wptr, " \n", 2);
	if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto glm_assoc_nosnp_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto glm_assoc_nosnp_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "Permutation test report written to %s.\n", outname);
    logprintb();
    if (mperm_save) {
      *outname_end = '.';
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	wptr = strcpya(&(outname_end[1]), &(param_names[param_idx * max_param_name_len]));
	memcpy(wptr, ".mperm.dump.all", 17);
	if (fopen_checked(&outfile, outname, "w")) {
	  goto glm_assoc_nosnp_ret_OPEN_FAIL;
	}
	ulii = param_ct - 1;
	if (pheno_d) {
	  wptr = memcpya(tbuf, "0 ", 2);
	  wptr = double_g_writex(wptr, orig_stats[param_idx - 1], '\n');
	  wptr2 = &(tbuf[MAXLINELEN]);
	  dptr = &(mperm_save_stats[param_idx - 1]);
	  for (perm_idx = 0; perm_idx < glm_mperm_val; perm_idx++) {
	    wptr = uint32_writex(wptr, perm_idx + 1, ' ');
	    dxx = dptr[perm_idx * ulii];
	    if (dxx >= 0) {
	      wptr = double_g_writex(wptr, dxx, '\n');
	    } else {
	      wptr = memcpyl3a(wptr, "NA\n");
	    }
	    if (wptr >= wptr2) {
	      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
		goto glm_assoc_nosnp_ret_WRITE_FAIL;
	      }
	      wptr = tbuf;
	    }
	  }
	  if (fwrite_checkedz(tbuf, wptr - tbuf, outfile)) {
	    goto glm_assoc_nosnp_ret_WRITE_FAIL;
	  }
	} else {
	  // todo
	}
	if (fclose_null(&outfile)) {
	  goto glm_assoc_nosnp_ret_WRITE_FAIL;
	}
      }
    }
  }
  while (0) {
  glm_assoc_nosnp_ret_NOMEM2:
    wkspace_left += topsize;
  glm_assoc_nosnp_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_assoc_nosnp_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_assoc_nosnp_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_assoc_nosnp_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  glm_assoc_nosnp_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  glm_assoc_nosnp_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 glm_assoc_nosnp_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  free_cond(active_params);
  free_cond(condition_uidxs);
  return retval;
}
