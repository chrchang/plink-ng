#include "wdist_common.h"
#include "wdist_stats.h"

// load markers in blocks to enable multithreading and, for quantitative
// phenotypes, PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 256
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

int32_t multcomp(char* outname, char* outname_end, uint32_t* marker_uidxs, uintptr_t chi_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, double* chi, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t non_chi, uint32_t* tcnt) {
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
	dyy = tprob(dxx, ujj - 2);
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
  qsort_ext((char*)sp, chi_ct, sizeof(double), double_cmp_deref, (char*)new_order, sizeof(uint32_t));
  if (tcnt) {
    qsort_ext((char*)schi, chi_ct, sizeof(double), double_cmp_deref, (char*)new_tcnt, sizeof(uint32_t));
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
      pv_gc[cur_idx] = tprob(sqrt(schi[uii] * lambda), new_tcnt[uii]);
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
  if (fprintf(outfile, tbuf, "SNP") < 0) {
    goto multcomp_ret_WRITE_FAIL;
  }
  if (qq_plot) {
    fputs("        QQ ", outfile);
  }
  tbuf[1] = ' ';
  tbuf[uii - 22] = '\0';
  if (fputs("      BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY\n", outfile) == EOF) {
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
      intprint2(&(tbuf[2]), (uint32_t)get_marker_chrom(chrom_info_ptr, marker_uidx));
      if (fprintf(outfile, tbuf, &(marker_ids[marker_uidx * max_marker_id_len])) < 0) {
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
      if (putc('\n', outfile) == EOF) {
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
  // Assumes perm_vec is initially zeroed out, tot_quotient is
  // 4294967296 / tot_ct, and totq_magic/totq_preshift/totq_postshift/totq_incr
  // have been precomputed from magic_num().
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

void calc_git(uint32_t pheno_nm_ct, uint32_t perm_vec_ct, uint32_t do_reverse, uintptr_t* __restrict__ loadbuf, uint32_t* perm_vecst, uint32_t* thread_bufs) {
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
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct128x4 = perm_ct128 * 4;
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
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
  uint32_t perm_ct4x4 = 4 * perm_ct4;
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
  uintptr_t xormask;
  uint32_t pbidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t indiv_type;
#ifdef __LP64__
  // thread_bufs structure:
  // first perm_ct8 16-byte blocks: homa1 cts
  // next perm_ct8 blocks: missing cts
  // next perm_ct8 blocks: het cts
  // afterward: free workspace, which is used for 4- and 8-bit partial counts
  gitv[0] = &(((__m128i*)thread_bufs)[3 * perm_ct16x4]);
  gitv[1] = &(((__m128i*)thread_bufs)[3 * perm_ct16x4 + perm_ct128x4]);
  gitv[2] = &(((__m128i*)thread_bufs)[3 * perm_ct16x4 + 2 * perm_ct128x4]);
  gitv[3] = &(((__m128i*)thread_bufs)[3 * (perm_ct16x4 + perm_ct128x4)]);
  gitv[4] = &(((__m128i*)thread_bufs)[3 * (perm_ct16x4 + perm_ct128x4) + 2 * perm_ct32]);
  gitv[5] = &(((__m128i*)thread_bufs)[3 * (perm_ct16x4 + perm_ct128x4) + 4 * perm_ct32]);
  if (!do_reverse) {
    gitv[6] = &(((__m128i*)thread_bufs)[2 * perm_ct16x4]);
    gitv[7] = &(((__m128i*)thread_bufs)[perm_ct16x4]);
    xormask = ~0LLU;
  } else {
    gitv[6] = &(((__m128i*)thread_bufs)[perm_ct16x4]);
    gitv[7] = &(((__m128i*)thread_bufs)[2 * perm_ct16x4]);
    xormask = 0;
  }
  gitv[8] = (__m128i*)thread_bufs;
#else
  // first perm_ct 4-byte blocks: homa1_cts
  // ...
  gitv[0] = &(thread_bufs[3 * perm_ct4x4]);
  gitv[1] = &(thread_bufs[3 * perm_ct4x4 + perm_ct32x4]);
  gitv[2] = &(thread_bufs[3 * perm_ct4x4 + 2 * perm_ct32x4]);
  gitv[3] = &(thread_bufs[3 * (perm_ct4x4 + perm_ct32x4)]);
  gitv[4] = &(thread_bufs[3 * (perm_ct4x4 + perm_ct32x4) + 2 * perm_ct8]);
  gitv[5] = &(thread_bufs[3 * (perm_ct4x4 + perm_ct32x4) + 4 * perm_ct8]);
  if (!do_reverse) {
    gitv[6] = &(thread_bufs[2 * perm_ct4x4]);
    gitv[7] = &(thread_bufs[perm_ct4x4]);
    xormask = ~0LU;
  } else {
    gitv[6] = &(thread_bufs[perm_ct4x4]);
    gitv[7] = &(thread_bufs[2 * perm_ct4x4]);
    xormask = 0;
  }
  gitv[8] = thread_bufs;
#endif
  cur_cts[0] = 0;
  cur_cts[1] = 0;
  cur_cts[2] = 0;
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulii = *loadbuf++;
    ulii ^= xormask;
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

void calc_qgit(uint32_t pheno_nm_ct, uintptr_t perm_vec_ctcl8m, uint32_t num_perms_now, uint32_t is_reverse, uintptr_t* __restrict__ loadbuf, double* perm_vecstd, double* thread_bufs) {
  uint32_t pheno_nm_ctl2x = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t xormask = ~ZEROLU;
  uint32_t het_code_postxor = 1;
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
  if (is_reverse) {
    xormask = 0;
    het_code_postxor = 2;
  }
  for (uii = 0; uii < pheno_nm_ctl2x; uii++) {
    ulii = *loadbuf++;
    ulii ^= xormask;
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
      if (indiv_type == het_code_postxor) {
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
      if (indiv_type == het_code_postxor) {
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

#ifdef __LP64__
uintptr_t qgit_ld_cost_40v(__m128i* vec1, __m128i* vend, __m128i* vec2) {
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
uintptr_t qgit_ld_cost_4(uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
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

uintptr_t qgit_ld_cost(uintptr_t indiv_ctl2, uintptr_t* loadbuf1, uintptr_t* loadbuf2) {
  // Cost: 3 * (missing <-> homrar/het) + 2 * (missing <-> homcom) +
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
  uintptr_t cost = 0;
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
  qgit_ld_cost_loop:
    lptr_4x_end = &(loadbuf1[cur_decr]);
    cost += qgit_ld_cost_40v((__m128i*)loadbuf1, (__m128i*)lptr_4x_end, (__m128i*)loadbuf2);
    loadbuf1 = lptr_4x_end;
    loadbuf2 = &(loadbuf2[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto qgit_ld_cost_loop;
  }
#else
  uintptr_t* lptr_four_end = &(loadbuf1[indiv_ctl2 & (~3U)]);
  while (loadbuf1 < lptr_four_end) {
    cost += qgit_ld_cost_4(loadbuf1, loadbuf2);
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

// multithread globals
static double* g_orig_1mpval;
static double* g_orig_chisq;

// A separated-low-and-high-bit format was tried, and found to not really be
// any better than the usual PLINK 2-bit format.
static uintptr_t* g_loadbuf;

static uintptr_t* g_perm_vecs;

static uint32_t* g_perm_vecst; // genotype indexing support
static uint32_t* g_thread_git_cts;

// always use genotype indexing for quantitative traits
static double* g_perm_vecstd;
static double* g_thread_git_qbufs;
static double* g_qresultbuf;
static double* g_pheno_d2;
static double g_pheno_sum;
static double g_pheno_ssq;

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

// QT --assoc uses this instead of g_set_cts
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
static uintptr_t* g_reverse;
static uint32_t g_model_fisher;
static uint32_t g_assoc_thread_ct;
static uintptr_t g_perm_vec_ct;
static uint32_t g_thread_block_ctl;
static uint32_t g_maxt_block_base;
static const uint32_t g_block_start = 0;
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

THREAD_RET_TYPE qassoc_gen_perms_thread(void* arg) {
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
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t* git_homclear_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 16 * perm_ct16]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 32 * perm_ct16]);
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t* git_homclear_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 4 * perm_ct4]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 8 * perm_ct4]);
#endif
  uintptr_t perm_vec_ctcl8 = (perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  uint32_t min_ploidy = 2;
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ nonmale_vec = g_indiv_nonmale_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t marker_idx;
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
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  if (is_haploid) { // includes g_is_x
    min_ploidy = 1;
  }
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
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
    gpui = &(precomp_ui[6 * precomp_width * marker_bidx]);
    if (model_fisher) {
      gpd = &(precomp_d[3 * precomp_width * marker_bidx]);
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
    } else {
      if (g_orig_1mpval[marker_idx] == -9) {
	perm_2success_ct[marker_idx++] += perm_vec_ct;
	continue;
      }
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      stat_high = orig_chisq[marker_idx] + EPSILON;
      stat_low = orig_chisq[marker_idx] - EPSILON;
    }
    success_2incr = 0;
    if (!is_x_or_y) {
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homclear_cts, perm_ct128 * 264);
#else
      fill_ulong_zero((uintptr_t*)git_homclear_cts, 132 * perm_ct32);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, 0, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homclear_cts);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!is_x_or_y) {
	if (!is_haploid) {
	  case_missing_ct = git_missing_cts[pidx];
	  case_set_ct = row1x_sum - (git_het_cts[pidx] + 2 * (case_missing_ct + git_homclear_cts[pidx]));
	} else {
	  case_missing_ct = git_missing_cts[pidx] + git_het_cts[pidx];
	  case_set_ct = row1x_sum - case_missing_ct - git_homclear_cts[pidx];
	}
      } else {
	if (is_x) {
	  vec_set_freq_x(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_set_ct, &case_missing_ct);
	} else {
	  vec_set_freq_y(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), nonmale_vec, &case_set_ct, &case_missing_ct);
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
  uintptr_t* reverse = g_reverse;
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
  uint32_t is_reverse;
  intptr_t g_sum;
  intptr_t g_ssq;
  uint32_t nanal;
  double nanal_recip;
  double nanal_m1_recip;
  double g_mean;
  double g_var;
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
      // g_var is zero, so we explicitly check for it here.
      // yes, that could also be caused by all-homrars if the user actively
      // tries to shoot themselves in the foot.  that's not my problem.
      perm_adapt_stop[marker_idx] = 1;
      perm_attempt_ct[marker_idx] = 0;
      continue;
    }
    is_reverse = is_set(reverse, marker_idx);
    homrar_ct = pheno_nm_ct - missing_ct - het_ct - homcom_ct;
    sval = orig_chiabs[marker_idx];
    // tstat = beta / vbeta_sqrt
    // tstat^2 = beta * beta / vbeta;
    //         = beta^2 * (nanal - 2) / ((qt_var / g_var) - beta^2)
    // [stop here for max(T) since nanal varies across markers]
    // tstat^2 / (nanal - 2) = beta^2 / ((qt_var / g_var) - beta^2)
    //                       = beta^2 * g_var / (qt_var - beta^2 * g_var)
    // Larger values of this last statistic monotonically result in smaller
    // P-values, so this is what we use for comparison (this saves a few
    // floating point operations at the end).
    sval = sval * sval / ((double)(((int32_t)nanal) - 2));
    stat_high = sval + EPSILON;
    stat_low = sval - EPSILON;
    g_sum = 2 * homrar_ct + het_ct;
    g_ssq = 4 * homrar_ct + het_ct;
    nanal_recip = 1.0 / ((double)((int32_t)nanal));
    nanal_m1_recip = 1.0 / ((double)(((int32_t)nanal) - 1));
    g_mean = ((double)g_sum) * nanal_recip;
    g_var = (((double)g_ssq) - g_sum * g_mean) * nanal_m1_recip;
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
	calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, next_cqg - pidx, is_reverse, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecstd[pidx]), &(git_qt_g_prod[pidx]));
      }
      qt_sum = pheno_sum - git_qt_sum[pidx];
      qt_ssq = pheno_ssq - git_qt_ssq[pidx];
      qt_g_prod = git_qt_g_prod[pidx];
      qt_mean = qt_sum * nanal_recip;
      qt_var = (qt_ssq - qt_sum * qt_mean) * nanal_m1_recip;
      qt_g_covar = (qt_g_prod - qt_sum * g_mean) * nanal_m1_recip;
      dxx = 1.0 / g_var;
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
	    fill_double_zero(git_qt_g_prod, pidx);
	    fill_double_zero(git_qt_sum, pidx);
	    fill_double_zero(git_qt_ssq, pidx);
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
  uint32_t marker_bceil = qblock_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  double* qresultbuf = g_qresultbuf;
  uintptr_t* loadbuf = g_loadbuf;
  uintptr_t* reverse = g_reverse;
  double* __restrict__ perm_vecstd = g_perm_vecstd;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  uint32_t* __restrict__ homcom_cts = g_homcom_cts;
  double* __restrict__ orig_chiabs = g_orig_chisq;
  double pheno_sum = g_pheno_sum;
  double pheno_ssq = g_pheno_ssq;
  double* git_qt_g_prod;
  double* git_qt_sum;
  double* git_qt_ssq;
  uintptr_t* loadbuf_cur;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t marker_bidx2;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t homrar_ct;
  uint32_t is_reverse;
  intptr_t g_sum;
  intptr_t g_ssq;
  uint32_t nanal;
  double nanal_recip;
  double nanal_m1_recip;
  double g_mean;
  double g_var;
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
  int32_t nanal_tmp;
  int32_t het_ct_tmp;
  int32_t homcom_ct_tmp;
  int32_t homrar_ct_tmp;
  uint32_t loop_ceil;
  uintptr_t cur_cost;
  uint32_t ldref;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  marker_idx = maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    missing_ct = missing_cts[marker_idx];
    nanal = pheno_nm_ct - missing_ct;
    homcom_ct = homcom_cts[marker_idx];
    het_ct = het_cts[marker_idx];
    if ((nanal < 3) || (homcom_ct == nanal) || (het_ct == nanal)) {
      perm_2success_ct[marker_idx] += perm_vec_ct;
      continue;
    }
    is_reverse = is_set(reverse, marker_idx);
    homrar_ct = nanal - het_ct - homcom_ct;
    sval = orig_chiabs[marker_idx];
    sval = sval * sval;
    stat_high = sval + EPSILON;
    stat_low = sval - EPSILON;
    g_sum = 2 * homrar_ct + het_ct;
    g_ssq = 4 * homrar_ct + het_ct;
    nanal_recip = 1.0 / ((double)((int32_t)nanal));
    nanal_m1_recip = 1.0 / ((double)(((int32_t)nanal) - 1));
    nanal_m2d = nanal - 2;
    g_mean = ((double)g_sum) * nanal_recip;
    g_var = (((double)g_ssq) - g_sum * g_mean) * nanal_m1_recip;
    success_2incr = 0;
    git_qt_g_prod = &(qresultbuf[3 * marker_bidx * perm_vec_ctcl8m]);
    git_qt_sum = &(git_qt_g_prod[perm_vec_ctcl8m]);
    git_qt_ssq = &(git_qt_g_prod[perm_vec_ctcl8m]);

    // Check if PERMORY-style LD exploitation is better than genotype indexing
    // algorithm.
    // Addition loops required for genotype indexing:
    //   het_ct + homrar_ct + 2 * missing_ct
    //
    // Addition loops required for LD exploitation:
    //   3 * (missing <-> homrar/het) + 2 * (missing <-> homcom) +
    //   (homrar <-> het/homcom) + (het <-> homcom)
    // Simple lower bound (may allow us to skip full LD cost calculation):
    //   (delta(homrar) + 2 * delta(missing) + delta(het) + delta(homcom)) / 2
    best_cost = het_ct + homrar_ct + 2 * missing_ct;
    ldref = marker_bidx;
    marker_idx_tmp = maxt_block_base;
    loop_ceil = maxt_block_base2;
    loadbuf_cur = &(loadbuf[marker_bidx * pheno_nm_ctl2]);
    do {
      if (marker_idx_tmp == maxt_block_base2) {
	marker_idx_tmp = maxt_block_base3;
	loop_ceil = marker_idx;
      }
      for (; marker_idx_tmp < loop_ceil; marker_idx_tmp++) {
	missing_ct_tmp = missing_cts[marker_idx_tmp];
	nanal_tmp = pheno_nm_ct - missing_ct_tmp;
	homcom_ct_tmp = homcom_cts[marker_idx_tmp];
	het_ct_tmp = het_cts[marker_idx_tmp];
	homrar_ct_tmp = nanal_tmp - het_ct_tmp - homcom_ct_tmp;
	if ((nanal_tmp >= 3) && (homcom_ct_tmp < nanal_tmp) && (het_ct_tmp < nanal_tmp)) {
	  cur_cost = labs(((int32_t)missing_ct) - missing_ct_tmp) + (labs(((int32_t)homrar_ct) - homrar_ct_tmp) + labs(((int32_t)het_ct) - het_ct_tmp) + labs(((int32_t)homcom_ct) - homcom_ct_tmp) + 1) / 2;
	  if (cur_cost < best_cost) {
	    marker_bidx2 = marker_idx_tmp - marker_idx;
	    cur_cost = qgit_ld_cost(pheno_nm_ctl2, &(loadbuf[marker_bidx2 * pheno_nm_ctl2]), loadbuf_cur);
	    if (cur_cost < best_cost) {
	      ldref = marker_bidx2;
	      best_cost = cur_cost;
	    }
	  }
	}
      }
    } while (marker_idx_tmp == maxt_block_base2);
    if (1) {
    // if (ldref == marker_bidx) {
      fill_double_zero(git_qt_g_prod, perm_vec_ctcl8m * 3);
      calc_qgit(pheno_nm_ct, perm_vec_ctcl8m, perm_vec_ct, is_reverse, loadbuf_cur, perm_vecstd, git_qt_g_prod);
    } else {
      /*
      printf("%lu %u %u\n", );
      exit(1);
      */
      // calc_qrem();
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      qt_sum = pheno_sum - git_qt_sum[pidx];
      qt_ssq = pheno_ssq - git_qt_ssq[pidx];
      qt_g_prod = git_qt_g_prod[pidx];
      qt_mean = qt_sum * nanal_recip;
      qt_var = (qt_ssq - qt_sum * qt_mean) * nanal_m1_recip;
      qt_g_covar = (qt_g_prod - qt_sum * g_mean) * nanal_m1_recip;
      dxx = 1.0 / g_var;
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
  uint32_t is_x = g_is_x;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* reverse = g_reverse;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
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
  uint32_t cur_xor;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    next_adapt_check = first_adapt_check;
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    cur_xor = is_model_prec ^ is_set(reverse, marker_idx);
    if (cur_xor) {
      col2_sum = (set_cts[marker_idx] + het_cts[marker_idx]) / 2;
      col1_sum = tot_obs - col2_sum;
    } else {
      col1_sum = (set_cts[marker_idx] - het_cts[marker_idx]) / 2;
      col2_sum = tot_obs - col1_sum;
    }
    missing_start = precomp_start[marker_bidx];
    gpui = &(precomp_ui[4 * precomp_width * marker_bidx]);
    success_2start = perm_2success_ct[marker_idx];
    success_2incr = 0;
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
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!is_x) {
	vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &uii, &case_homx_ct);
      } else {
	vec_3freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_missing_ct, &uii, &case_homx_ct);
      }
      if (cur_xor) {
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
  uint32_t is_x = g_is_x;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 16 * perm_ct16]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 32 * perm_ct16]);
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 4 * perm_ct4]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 8 * perm_ct4]);
#endif
  uintptr_t perm_vec_ctcl8 = (perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  int32_t is_model_prec = g_is_model_prec;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* reverse = g_reverse;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_homx_ct;
  uint32_t case_missing_ct;
  uint32_t cur_reverse;
  uint32_t cur_xor;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double stat_high;
  double stat_low;
  double sval;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    cur_reverse = is_set(reverse, marker_idx);
    cur_xor = is_model_prec ^ cur_reverse;
    if (cur_xor) {
      col2_sum = (set_cts[marker_idx] + het_cts[marker_idx]) / 2;
      col1_sum = tot_obs - col2_sum;
    } else {
      col1_sum = (set_cts[marker_idx] - het_cts[marker_idx]) / 2;
      col2_sum = tot_obs - col1_sum;
    }
    missing_start = precomp_start[marker_bidx];
    gpui = &(precomp_ui[6 * precomp_width * marker_bidx]);
    if (model_fisher) {
      if (orig_1mpval[marker_idx] == -9) {
	marker_idx++;
	continue;
      }
      gpd = &(precomp_d[3 * precomp_width * marker_bidx]);
      stat_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
      stat_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
    } else {
      if (orig_chisq[marker_idx] == -9) {
	marker_idx++;
	continue;
      }
      gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
      stat_high = orig_chisq[marker_idx] + EPSILON;
      stat_low = orig_chisq[marker_idx] - EPSILON;
    }
    success_2incr = 0;
    if (!is_x) {
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, perm_ct128 * 264);
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 132 * perm_ct32);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, cur_reverse, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homrar_cts);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!is_x) {
	case_missing_ct = git_missing_cts[pidx];
	if (is_model_prec) {
	  case_homx_ct = git_homrar_cts[pidx];
	} else {
	  case_homx_ct = case_ct - case_missing_ct - git_homrar_cts[pidx] - git_het_cts[pidx];
	}
      } else {
	vec_3freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_missing_ct, &uii, &case_homx_ct);
	if (cur_xor) {
	  case_homx_ct = case_ct - case_homx_ct - case_missing_ct - uii;
	}
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
  uint32_t is_x = g_is_x;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* __restrict__ reverse = g_reverse;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
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
  uint32_t cur_reverse;
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
    cur_reverse = is_set(reverse, marker_idx);
    if (!cur_reverse) {
      homcom_ct = (set_cts[marker_idx] - het_ct) / 2;
    } else {
      homcom_ct = tot_obs - ((set_cts[marker_idx] + het_ct) / 2);
    }
    missing_start = precomp_start[marker_bidx];
    gpui = &(precomp_ui[4 * precomp_width * marker_bidx]);
    success_2start = perm_2success_ct[marker_idx];
    success_2incr = 0;
    chisq_high = orig_chisq[marker_idx] + EPSILON;
    chisq_low = orig_chisq[marker_idx] - EPSILON;
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!is_x) {
	vec_set_freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_com_ct, &case_missing_ct);
      } else {
        vec_set_freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_com_ct, &case_missing_ct);
      }
      if (cur_reverse) {
	case_com_ct = (case_ct - case_missing_ct) * 2 - case_com_ct;
      }
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
  uint32_t is_x = g_is_x;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 16 * perm_ct16]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 32 * perm_ct16]);
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 4 * perm_ct4]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 8 * perm_ct4]);
#endif
  uintptr_t perm_vec_ctcl8 = (perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* __restrict__ reverse = g_reverse;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t tot_obs;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t cur_reverse;
  uint32_t case_com_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double chisq_high;
  double chisq_low;
  double chisq;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (orig_1mpval[marker_idx] == -9) {
      perm_2success_ct[marker_idx++] += perm_vec_ct;
      continue;
    }
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    het_ct = het_cts[marker_idx];
    cur_reverse = is_set(reverse, marker_idx);
    if (!cur_reverse) {
      homcom_ct = (set_cts[marker_idx] - het_ct) / 2;
    } else {
      homcom_ct = tot_obs - ((set_cts[marker_idx] + het_ct) / 2);
    }
    missing_start = precomp_start[marker_bidx];
    gpui = &(precomp_ui[6 * precomp_width * marker_bidx]);
    gpd = &(precomp_d[2 * precomp_width * marker_bidx]);
    chisq_high = orig_chisq[marker_idx] + EPSILON;
    chisq_low = orig_chisq[marker_idx] - EPSILON;
    success_2incr = 0;
    if (!is_x) {
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, perm_ct128 * 264);
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 132 * perm_ct32);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, cur_reverse, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homrar_cts);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!is_x) {
	case_missing_ct = git_missing_cts[pidx];
	case_com_ct = 2 * (case_ct - case_missing_ct - git_homrar_cts[pidx]) - git_het_cts[pidx];
      } else {
	vec_set_freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_com_ct, &case_missing_ct);
	if (cur_reverse) {
	  case_com_ct = 2 * (case_ct - case_missing_ct) - case_com_ct;
	}
      }
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
  uint32_t is_x = g_is_x;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
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
  intptr_t homset_ct;
  intptr_t homclear_ct;
  intptr_t het_ct;
  uint32_t case_missing_ct;
  uint32_t case_het_ct;
  uint32_t case_homset_ct;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = adapt_m_table[marker_bidx];
    next_adapt_check = first_adapt_check;
    het_ct = het_cts[marker_idx];
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    homset_ct = (set_cts[marker_idx] - het_ct) / 2;
    homclear_ct = tot_obs - het_ct - homset_ct;
    if (!homset_ct) {
      missing_col = 3;
    } else if ((het_ct + homset_ct == tot_obs) || (!het_ct)) {
      missing_col = 2; // either no hom A1s or no hets (no need to distinguish)
    } else {
      missing_col = 0;
    }
    success_2start = perm_2success_ct[marker_idx];
    success_2incr = 0;
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
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!is_x) {
	vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &case_het_ct, &case_homset_ct);
      } else {
	vec_3freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_missing_ct, &case_het_ct, &case_homset_ct);
      }
      if (model_fisher) {
        uii = case_ct - case_het_ct - case_homset_ct - case_missing_ct;
	// this is very slow.  a precomputed 2-dimensional table could improve
	// matters, but I doubt it's worth the effort for now.
	dxx = fisher23(case_homset_ct, case_het_ct, uii, homset_ct - case_homset_ct, het_ct - case_het_ct, homclear_ct - uii);
	if (dxx < stat_low) {
	  success_2incr += 2;
	} else if (dxx < stat_high) {
	  success_2incr++;
	}
      } else {
	if (!missing_col) {
	  dxx = chi23_eval(case_homset_ct, case_het_ct, case_ct - case_missing_ct, homset_ct, het_ct, tot_obs);
	} else if (missing_col == 3) {
	  dxx = chi22_eval(case_het_ct, case_ct - case_missing_ct, het_ct, tot_obs);
	} else {
	  dxx = chi22_eval(case_homset_ct, case_ct - case_missing_ct, homset_ct, tot_obs);
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
  uint32_t is_x = g_is_x;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 16 * perm_ct16]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 32 * perm_ct16]);
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 4 * perm_ct4]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 8 * perm_ct4]);
#endif
  uintptr_t perm_vec_ctcl8 = (perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* __restrict__ reverse = g_reverse;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t missing_col;
  intptr_t tot_obs;
  intptr_t homset_ct;
  intptr_t homclear_ct;
  intptr_t het_ct;
  uint32_t success_2incr;
  uint32_t case_missing_ct;
  uint32_t case_het_ct;
  uint32_t case_homset_ct;
  uint32_t cur_reverse;
  uint32_t uii;
  double stat_high;
  double stat_low;
  double sval;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - perm_vec_ct]), perm_vec_ct * sizeof(double));
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    het_ct = het_cts[marker_idx];
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    homset_ct = (set_cts[marker_idx] - het_ct) / 2;
    homclear_ct = tot_obs - het_ct - homset_ct;
    if (!homset_ct) {
      missing_col = 3;
    } else if ((het_ct + homset_ct == tot_obs) || (!het_ct)) {
      missing_col = 2;
    } else {
      missing_col = 0;
    }
    cur_reverse = is_set(reverse, marker_idx);
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
    success_2incr = 0;
    if (!is_x) {
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, perm_ct128 * 264);
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 132 * perm_ct32);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, cur_reverse, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homrar_cts);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!is_x) {
	case_missing_ct = git_missing_cts[pidx];
	case_het_ct = git_het_cts[pidx];
	if (!cur_reverse) {
	  case_homset_ct = case_ct - case_missing_ct - case_het_ct - git_homrar_cts[pidx];
	} else {
	  case_homset_ct = git_homrar_cts[pidx];
	}
      } else {
	vec_3freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_missing_ct, &case_het_ct, &case_homset_ct);
      }
      if (model_fisher) {
        uii = case_ct - case_het_ct - case_homset_ct - case_missing_ct;
	sval = fisher23(case_homset_ct, case_het_ct, uii, homset_ct - case_homset_ct, het_ct - case_het_ct, homclear_ct - uii);
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
	  sval = chi23_eval(case_homset_ct, case_het_ct, case_ct - case_missing_ct, homset_ct, het_ct, tot_obs);
	} else if (missing_col == 3) {
	  sval = chi22_eval(case_het_ct, case_ct - case_missing_ct, het_ct, tot_obs);
	} else {
	  sval = chi22_eval(case_homset_ct, case_ct - case_missing_ct, homset_ct, tot_obs);
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
  uint32_t is_x = g_is_x;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* __restrict__ reverse = g_reverse;
  uintptr_t* is_invalid = g_is_invalid;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
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
  uint32_t cur_reverse;
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
    next_adapt_check = first_adapt_check;
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    het_ct = het_cts[marker_idx];
    cur_reverse = is_set(reverse, marker_idx);
    if (!cur_reverse) {
      com_ct = set_cts[marker_idx];
    } else {
      com_ct = tot_obs * 2 - set_cts[marker_idx];
    }
    homrar_ct = tot_obs - ((com_ct + het_ct) / 2);
    homcom_ct = (com_ct - het_ct) / 2;
    missing_start = precomp_start[marker_bidx];
    skip_domrec = is_set(is_invalid, marker_idx);
    gpui = &(precomp_ui[12 * precomp_width * marker_bidx]);
    success_2start = perm_2success_ct[marker_idx];
    success_2incr = 0;
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
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!is_x) {
	vec_3freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), &case_missing_ct, &case_het_ct, &case_homcom_ct);
      } else {
	vec_3freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_missing_ct, &case_het_ct, &case_homcom_ct);
      }
      case_homrar_ct = case_ct - case_missing_ct - case_het_ct - case_homcom_ct;
      if (cur_reverse) {
	uii = case_homrar_ct;
	case_homrar_ct = case_homcom_ct;
	case_homcom_ct = uii;
      }
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
  uint32_t is_x = g_is_x;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t model_fisher = g_model_fisher;
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 16 * perm_ct16]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct128 * 528 + 32 * perm_ct16]);
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t* git_homrar_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132]);
  uint32_t* git_missing_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 4 * perm_ct4]);
  uint32_t* git_het_cts = &(g_thread_git_cts[tidx * perm_ct32 * 132 + 8 * perm_ct4]);
#endif
  uintptr_t perm_vec_ctcl8 = (perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  uint32_t precomp_width = g_precomp_width;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ male_vec = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* reverse = g_reverse;
  uintptr_t* is_invalid = g_is_invalid;
  uint32_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  uint32_t* __restrict__ het_cts = g_het_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  double* __restrict__ orig_chisq = g_orig_chisq;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t tot_obs;
  intptr_t com_ct;
  intptr_t rar_ct;
  intptr_t het_ct;
  intptr_t homrar_ct;
  intptr_t homcom_ct;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_homrar_ct;
  uint32_t case_homcom_ct;
  uint32_t case_het_ct;
  uint32_t case_missing_ct;
  uint32_t case_com_ct;
  uint32_t cur_reverse;
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
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    tot_obs = pheno_nm_ct - missing_cts[marker_idx];
    het_ct = het_cts[marker_idx];
    cur_reverse = is_set(reverse, marker_idx);
    if (!cur_reverse) {
      com_ct = set_cts[marker_idx];
      rar_ct = tot_obs * 2 - com_ct;
    } else {
      rar_ct = set_cts[marker_idx];
      com_ct = tot_obs * 2 - rar_ct;
    }
    homrar_ct = tot_obs - ((com_ct + het_ct) / 2);
    homcom_ct = (com_ct - het_ct) / 2;
    missing_start = precomp_start[marker_bidx];
    skip_domrec = is_set(is_invalid, marker_idx);
    gpui = &(precomp_ui[18 * precomp_width * marker_bidx]);
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
    success_2incr = 0;
    if (!is_x) {
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)git_homrar_cts, perm_ct128 * 264);
#else
      fill_ulong_zero((uintptr_t*)git_homrar_cts, 132 * perm_ct32);
#endif
      calc_git(pheno_nm_ct, perm_vec_ct, cur_reverse, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, git_homrar_cts);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!is_x) {
	case_missing_ct = git_missing_cts[pidx];
	case_het_ct = git_het_cts[pidx];
	case_homrar_ct = git_homrar_cts[pidx];
	case_homcom_ct = case_ct - case_missing_ct - case_het_ct - case_homrar_ct;
      } else {
	vec_3freq_xx(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctl2]), male_vec, &case_missing_ct, &case_het_ct, &case_homcom_ct);
	case_homrar_ct = case_ct - case_missing_ct - case_het_ct - case_homcom_ct;
	if (cur_reverse) {
	  uii = case_homrar_ct;
	  case_homrar_ct = case_homcom_ct;
	  case_homcom_ct = uii;
	}
      }
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

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  int32_t retval = 0;
  FILE* outfile = NULL;
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
  uint32_t model_fisherx = (model_modifier & MODEL_FISHER) && (!(model_modifier & MODEL_PTREND));
  uint32_t mu_table[MODEL_BLOCKSIZE];
  uint32_t uibuf[4];
  char wprefix[5];
  char wbuf[48];
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
  uintptr_t indiv_idx;
  uint32_t* missp;
  uint32_t* setp;
  uint32_t* hetp;
  double* o1mpptr;
  double* ooptr;
  unsigned char* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* indiv_ctrl_include2;
  uintptr_t* indiv_case_include2;
  uint32_t load_indiv_ct;
  uint32_t is_reverse;
  uint64_t ullii;
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
  g_orig_chisq = NULL;
  g_model_fisher = model_modifier & MODEL_FISHER;
  g_pheno_nm_ct = pheno_nm_ct;
  g_perms_done = 0;
  g_aperm_alpha = aperm_alpha;
  g_reverse = marker_reverse;
  g_is_model_prec = (model_modifier & MODEL_PREC)? 1 : 0;
  g_is_y = 0;

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
  } else if (model_adapt) {
    perms_total = aperm_max;
  }
  if (wkspace_alloc_uc_checked(&loadbuf_raw, unfiltered_indiv_ct4)) {
    goto model_assoc_ret_NOMEM;
  }
  memset(wprefix, 32, 5);
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
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (assoc_counts) {
      if (fputs("     C_A      C_U   A2 ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    } else {
      if (fputs("     F_A      F_U   A2 ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    }
    if (!g_model_fisher) {
      if (fputs("       CHISQ ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    }
    if (fputs("           P           OR ", outfile) == EOF) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (display_ci) {
      uii = (uint32_t)(ci_size * 100);
      if (uii >= 10) {
	if (fprintf(outfile, "          SE          L%u          U%u ", uii, uii) < 0) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      } else {
	if (fprintf(outfile, "          SE           L%u           U%u ", uii, uii) < 0) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      }
    }
    if (putc('\n', outfile) == EOF) {
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
      retval = RET_INVALID_CMDLINE;
      goto model_assoc_ret_1;
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
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (!g_model_fisher) {
      if (fputs("       CHISQ   DF ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    } else {
      outname_end2 = memcpyb(outname_end2, ".fisher", 8);
    }
    if (fputs("           P\n", outfile) == EOF) {
      goto model_assoc_ret_WRITE_FAIL;
    }
  }
  marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  g_adaptive_ci_zt = ltqnorm(1 - aperm_beta / (2.0 * marker_ct));
  if (wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * pheno_nm_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&g_orig_1mpval, marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&g_missing_cts, marker_ct * sizeof(uint32_t)) ||
      wkspace_alloc_ui_checked(&g_set_cts, marker_ct * sizeof(uint32_t))) {
    goto model_assoc_ret_NOMEM;
  }
  if (model_assoc) {
    if (wkspace_alloc_d_checked(&orig_odds, marker_ct * sizeof(double))) {
      goto model_assoc_ret_NOMEM;
    }
  } else {
    if (wkspace_alloc_ui_checked(&g_het_cts, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
  }
  x_code = species_x_code[chrom_info_ptr->species];
  y_code = species_y_code[chrom_info_ptr->species];
  ullii = chrom_info_ptr->chrom_mask;
  gender_req = ((x_code != -1) && (ullii & (1LLU << x_code))) || (model_assoc && (((y_code != -1) && (ullii & (1LLU << y_code)))));
  if (gender_req) {
    if (wkspace_alloc_ul_checked(&g_indiv_nonmale_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&g_indiv_male_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
      goto model_assoc_ret_NOMEM;
    }
    fill_ulong_zero(g_indiv_male_include2, pheno_nm_ctl2);
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
      indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
      if (is_set(sex_male, indiv_uidx)) {
	set_bit_noct(g_indiv_male_include2, indiv_idx * 2);
	male_ct++;
      }
      indiv_uidx++;
    }
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

    g_precomp_width = (1 + (int32_t)(sqrt(pheno_nm_ct) * EXPECTED_MISSING_FREQ * 5.65686));
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
    g_tot_quotient = 4294967296LLU / pheno_nm_ct;
    magic_num(g_tot_quotient, &g_totq_magic, &g_totq_preshift, &g_totq_postshift, &g_totq_incr);
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
  for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
    indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
    if (is_set(pheno_c, indiv_uidx)) {
      set_bit_noct(indiv_case_include2, indiv_idx * 2);
      g_case_ct++;
    }
    indiv_uidx++;
  }
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
    for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
      indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
      if (is_set(sex_male, indiv_uidx) && is_set(pheno_c, indiv_uidx)) {
	set_bit_noct(indiv_male_case_include2, indiv_idx * 2);
	case_male_ct++;
      }
      indiv_uidx++;
    }
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
      //   g_thread_git_cts: ((g_perm_vec_ct + 127) / 128) * 2112 * thread_ct
      //   g_perm_vecs: pheno_nm_ctl2 * sizeof(intptr_t) * g_perm_vec_ct
      // If we force g_perm_vec_ct to be a multiple of 128, then we have
      //   g_perm_vec_ct * (24.5 * g_thread_ct + pheno_nm_ct +
      //                    sizeof(intptr_t) * pheno_nm_ctl2)
      //
      // "2112?  24.5?  wtf?"
      // Excellent questions.  Each max(T) thread has nine buffers to support
      // rapid execution of the genotype indexing algorithm:
      //   three with 4-bit accumulators, total size perm_vec_ct / 2 bytes
      //   three with 8-bit accumulators, total size perm_vec_ct bytes
      //   three with 32-bit accumulators, total size 4 * perm_vec_ct bytes
      // The initial 3 multiplier is to allow heterozygotes, homozygote minors,
      // and missing genotypes to be counted simultaneously.
      // Adding all this up, we have 16.5 * perm_vec_ct bytes, and multiplying
      // by 128 yields 2112.  The other thread_ct dependence contributes
      // 8 * perm_vec_ct bytes, multiplying by 128 yields 1024, and
      // 2112 + 1024 = 3136.
      g_perm_vec_ct = 128 * (wkspace_left / (128LL * sizeof(intptr_t) * pheno_nm_ctl2 + 3136LL * g_thread_ct + 16LL * pheno_nm_ct));
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    } else if (!g_perm_vec_ct) {
      goto model_assoc_ret_NOMEM;
    }
    g_perms_done += g_perm_vec_ct;
    g_perm_vecs = (uintptr_t*)wkspace_alloc(g_perm_vec_ct * pheno_nm_ctl2 * sizeof(intptr_t));
    if (g_perm_vec_ct > g_thread_ct) {
      g_assoc_thread_ct = g_thread_ct;
    } else {
      g_assoc_thread_ct = g_perm_vec_ct;
    }
    if (spawn_threads(threads, &model_assoc_gen_perms_thread, g_assoc_thread_ct)) {
      logprint(errstr_thread_create);
      goto model_assoc_ret_THREAD_CREATE_FAIL;
    }
    ulii = 0;
    model_assoc_gen_perms_thread((void*)ulii);
    join_threads(threads, g_assoc_thread_ct);
    if (!model_adapt) {
      ulii = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
      g_maxt_thread_results = (double*)wkspace_alloc(g_thread_ct * ulii * CACHELINE);
#ifdef __LP64__
      ulii = ((g_perm_vec_ct + 127) / 128) * 16;
#else
      ulii = ((g_perm_vec_ct + 31) / 32) * 4;
#endif
      g_perm_vecst = (uint32_t*)wkspace_alloc(ulii * pheno_nm_ct);
      g_thread_git_cts = (uint32_t*)wkspace_alloc(ulii * 132 * g_thread_ct);
      transpose_perms(g_perm_vecs, g_perm_vec_ct, pheno_nm_ct, g_perm_vecst);
#ifdef __LP64__
      fill_ulong_zero((uintptr_t*)g_thread_git_cts, (ulii * 33 * g_thread_ct) / 2);
#else
      fill_ulong_zero((uintptr_t*)g_thread_git_cts, ulii * 33 * g_thread_ct);
#endif
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
      // g_block_start = 0;
      if (model_assoc) {
	// exploit overflow
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_haploid);
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	g_is_y = (uii == (uint32_t)y_code);
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
	  if ((!((species_haploid_mask[chrom_info_ptr->species] >> uii) & 1LLU)) || g_is_x) {
	    break;
	  }
	  marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
	}
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto model_assoc_ret_READ_FAIL;
	}
	if (!g_is_x) {
	  cur_ctrl_include2 = indiv_ctrl_include2;
	  cur_case_include2 = indiv_case_include2;
	  load_indiv_ct = pheno_nm_ct;
	  load_case_ct = g_case_ct;
	} else {
	  cur_ctrl_include2 = indiv_nonmale_ctrl_include2;
	  cur_case_include2 = indiv_nonmale_case_include2;
	  load_indiv_ct = nonmale_ct;
	  load_case_ct = case_nonmale_ct;
	}
	load_ctrl_ct = load_indiv_ct - load_case_ct;
      }
      intprint2(&(wprefix[2]), uii);
      /*
    } else if (model_maxt) {
      g_block_start = 0;
      // todo: copy LD table, etc.
      // could conditionally do this for adaptive permutation, but I doubt
      // it's worth the trouble
    } else {
      g_block_start = 0;
      */
    }
    block_size = g_block_start;
    block_end = marker_unstopped_ct + g_block_start - marker_idx;
    if (block_end > MODEL_BLOCKSIZE) {
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
      if (fread(loadbuf_raw, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto model_assoc_ret_READ_FAIL;
      }
      indiv_uidx = 0;
      ukk = pheno_nm_ct / BITCT2;
      loadbuf_ptr = &(g_loadbuf[block_size * pheno_nm_ctl2]);
      for (uii = 0; uii < ukk; uii++) {
	ulii = 0;
	for (ujj = 0; ujj < BITCT; ujj += 2) {
	  indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	  ulii |= ((uintptr_t)(((loadbuf_raw[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3))) << ujj;
	  indiv_uidx++;
	}
	*loadbuf_ptr++ = ulii;
      }
      ujj = 2 * (pheno_nm_ct & (BITCT2 - 1));
      if (ujj) {
	ulii = 0;
	for (uii = 0; uii < ujj; uii += 2) {
	  indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	  ulii |= ((uintptr_t)(((loadbuf_raw[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3))) << uii;
	  indiv_uidx++;
	}
	*loadbuf_ptr = ulii;
      }
      if (model_adapt) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      if (is_set(marker_exclude, ++marker_uidx)) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
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
      setp = &(g_set_cts[marker_idx + g_block_start]);
      if (model_assoc) {
	ooptr = &(orig_odds[marker_idx + g_block_start]);
	for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	  marker_uidx2 = mu_table[marker_bidx];
	  g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	  is_reverse = is_set(marker_reverse, marker_uidx2);
	  if (!g_is_haploid) {
	    if (is_reverse) {
	      single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &ujj, &uii, &umm, &ukk);
	      *missp = uii + ukk;
	      *setp = ujj + umm;
	      uii = 2 * (ctrl_ct - uii) - ujj;
	      ukk = 2 * (g_case_ct - ukk) - umm;
	    } else {
	      single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &uii, &ujj, &ukk, &umm);
	      *missp = ujj + umm;
	      *setp = uii + ukk;
	      ujj = 2 * (ctrl_ct - ujj) - uii;
	      umm = 2 * (g_case_ct - umm) - ukk;
	    }
	  } else if (g_is_x) {
	    if (is_reverse) {
	      single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_nonmale_ctrl_include2, indiv_nonmale_case_include2, &ujj, &uii, &umm, &ukk);
	      *missp = 2 * (uii + ukk);
	      *setp = ujj + umm;
	      uii = 2 * (ctrl_nonmale_ct - uii) - ujj;
	      ukk = 2 * (case_nonmale_ct - ukk) - umm;
	      haploid_single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_male_ctrl_include2, indiv_male_case_include2, &uoo, &unn, &uqq, &upp);
	      *missp += unn + upp + male_ct;
	      *setp += uoo + uqq;
	      unn = ctrl_male_ct - unn - uoo;
	      upp = case_male_ct - upp - uqq;
	    } else {
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
	    }
	    uii += unn;
	    ujj += uoo;
	    ukk += upp;
	    umm += uqq;
	  } else if (g_is_haploid) {
	    if (is_reverse) {
	      haploid_single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &ujj, &uii, &umm, &ukk);
	      *missp = uii + ukk;
	      *setp = ujj + umm;
	      uii = load_ctrl_ct - uii - ujj;
	      ukk = load_case_ct - ukk - umm;
	    } else {
	      haploid_single_marker_cc_freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm);
	      *missp = ujj + umm;
	      *setp = uii + ukk;
	      ujj = load_ctrl_ct - ujj - uii;
	      umm = load_case_ct - umm - ukk;
	    }
	    if (g_is_y) {
	      *missp += nonmale_ct;
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
	  if ((pval <= pfilter) || (pfilter == 1.0)) {
	    a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
	    a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
	    memcpy(tbuf, wprefix, 5);
	    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), &(tbuf[5]));
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
	    memcpy(wptr, " \n", 2);
	    if (fwrite_checked(tbuf, 2 + (wptr - tbuf), outfile)) {
	      goto model_assoc_ret_WRITE_FAIL;
	    }
	  }
	  missp++;
	  setp++;
	  o1mpptr++;
	  ooptr++;
	}
      } else {
	hetp = &(g_het_cts[marker_idx + g_block_start]);
	for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	  marker_uidx2 = mu_table[marker_bidx];
	  g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	  is_reverse = is_set(marker_reverse, marker_uidx2);
	  if (is_reverse) {
	    single_marker_cc_3freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &ukk, &ujj, &uii, &uoo, &unn, &umm);
	    *missp = uii + umm;
	    *setp = ukk + uoo;
	    uii = load_ctrl_ct - uii - ujj - ukk;
	    umm = load_case_ct - umm - unn - uoo;
	  } else {
	    single_marker_cc_3freqs(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm, &unn, &uoo);
	    *missp = ukk + uoo;
	    *setp = uii + umm;
	    ukk = load_ctrl_ct - uii - ujj - ukk;
	    uoo = load_case_ct - umm - unn - uoo;
	  }
	  if (g_is_x) {
	    *missp += male_ct;
	  }
	  *hetp = ujj + unn;
	  *setp = 2 * (*setp) + (*hetp);
	  is_invalid = (uoo < model_cell_ct) || (unn < model_cell_ct) || (umm < model_cell_ct) || (ukk < model_cell_ct) || (ujj < model_cell_ct) || (uii < model_cell_ct);
	  a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
	  a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
	  memcpy(tbuf, wprefix, 5);
	  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), &(tbuf[5]));
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
	      set_bit_noct(g_is_invalid, marker_idx + marker_bidx);
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
	  is_reverse = is_set(marker_reverse, urr);
	  if (!is_reverse) {
	    uoo = g_set_cts[urr];
	  } else {
	    uoo = uqq - g_set_cts[urr];
	  }
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
	    if (!is_set(g_is_invalid, urr)) {
	      upp = pheno_nm_ct - upp;
	      if (!is_reverse) {
	        uoo = (g_set_cts[urr] - g_het_cts[urr]) / 2;
		uqq = upp - ((g_set_cts[urr] + g_het_cts[urr]) / 2);
	      } else {
		uoo = upp - ((g_set_cts[urr] + g_het_cts[urr]) / 2);
		uqq = (g_set_cts[urr] - g_het_cts[urr]) / 2;
	      }
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
	    if (!is_set(g_is_invalid, urr)) {
	      upp = pheno_nm_ct - upp;
	      if (!is_reverse) {
		uoo = (g_set_cts[urr] - g_het_cts[urr]) / 2;
		uqq = upp - ((g_set_cts[urr] + g_het_cts[urr]) / 2);
	      } else {
		uoo = upp - ((g_set_cts[urr] + g_het_cts[urr]) / 2);
		uqq = (g_set_cts[urr] - g_het_cts[urr]) / 2;
	      }
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
	  is_reverse = is_set(marker_reverse, urr);
	  if (!is_reverse) {
	    uoo = (g_set_cts[urr] - unn) / 2; // homcom_ct
	  } else {
	    uoo = upp - ((g_set_cts[urr] + unn) / 2);
	  }
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
	  is_reverse = is_set(marker_reverse, urr);
	  if ((model_modifier & MODEL_PREC) ^ (MODEL_PREC * is_reverse)) {
	    uoo = upp - ((g_set_cts[urr] + g_het_cts[urr]) / 2); // col1_sum
	  } else {
	    uoo = (g_set_cts[urr] - g_het_cts[urr]) / 2;
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
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  assoc_adapt_thread((void*)ulii);
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads(threads, &model_adapt_domrec_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_domrec_thread((void*)ulii);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads(threads, &model_adapt_trend_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	} else if (model_modifier & MODEL_PGEN) {
	  if (spawn_threads(threads, &model_adapt_gen_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_gen_thread((void*)ulii);
	} else {
	  if (spawn_threads(threads, &model_adapt_best_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_adapt_best_thread((void*)ulii);
	}
	join_threads(threads, g_assoc_thread_ct);
      } else {
	g_maxt_block_base = marker_idx - g_block_start;
	ulii = 0;
	if (model_assoc) {
	  if (spawn_threads(threads, &assoc_maxt_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  assoc_maxt_thread((void*)ulii);
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads(threads, &model_maxt_domrec_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_domrec_thread((void*)ulii);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads(threads, &model_maxt_trend_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_trend_thread((void*)ulii);
	} else if (model_modifier & MODEL_PGEN) {
	  if (spawn_threads(threads, &model_maxt_gen_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_gen_thread((void*)ulii);
	} else {
	  if (spawn_threads(threads, &model_maxt_best_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
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
    marker_idx += block_size - g_block_start;
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
      retval = multcomp(outname, outname_end, g_marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, g_model_fisher? g_orig_1mpval : g_orig_chisq, pfilter, mtest_adjust, adjust_lambda, g_model_fisher, NULL);
      if (retval) {
	goto model_assoc_ret_1;
      }
    }
  }
  if (model_perms) {
    wkspace_reset((unsigned char*)g_perm_vecs);
    if (g_perms_done < perms_total) {
      if (model_adapt) {
	marker_unstopped_ct = marker_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, 0, (marker_ct + sizeof(uintptr_t) - 1) / sizeof(uintptr_t));
	if (!marker_unstopped_ct) {
	  goto model_assoc_adapt_perm_count;
	}
      }
      printf("\r%u permutations complete.", g_perms_done);
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
    sprintf(logbuf, "%u %s permutations complete.\n", g_perms_done, model_maxt? "max(T)" : "(adaptive)");
    logprintb();
    if (g_model_fisher && (model_modifier & MODEL_PTREND)) {
      outname_end2 -= 7; // remove ".fisher"
      g_model_fisher = 0;
    }
    if (model_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
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
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_non_set_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    memset(tbuf, 32, 5);
    tbuf[5 + plink_maxsnp] = ' ';
    while (1) {
      while (1) {
	do {
          chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
	uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	g_is_x = (uii == (uint32_t)x_code);
	if (model_assoc || ((!((species_haploid_mask[chrom_info_ptr->species] >> uii) & 1LLU)) || g_is_x)) {
	  break;
	}
	marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
      }
      intprint2(&(tbuf[2]), uii);
      for (; marker_uidx < chrom_end;) {
	if (model_adapt) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) * dxx;
	}
        if ((pval <= pfilter) || (pfilter == 1.0)) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), &(tbuf[5]));
	  wptr = &(tbuf[6 + plink_maxsnp]);
	  if ((!model_assoc) && ((model_adapt && (!g_perm_attempt_ct[marker_idx])) || ((!model_adapt) && ((g_model_fisher && (g_orig_1mpval[marker_idx] == -9)) || ((!g_model_fisher) && (g_orig_chisq[marker_idx] == -9)))))) {
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
	      if (g_model_fisher) {
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
	  memcpy(wptr, " \n", 3);
	  if (fwrite_checked(tbuf, 2 + (uintptr_t)(wptr - tbuf), outfile)) {
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
  model_assoc_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 model_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t qassoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_mperm_val, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t xmhh_exists, uint32_t nxmhh_exists, uint32_t perm_batch_size) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = 2 * ((unfiltered_indiv_ct + BITCT - 1) / BITCT);
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  FILE* outfile = NULL;
  FILE* outfile_qtm = NULL;
  uint32_t perm_adapt = model_modifier & MODEL_PERM;
  uint32_t perm_maxt = model_modifier & MODEL_MPERM;
  uint32_t do_perms = perm_adapt | perm_maxt;
  uint32_t qt_means = model_modifier & MODEL_QT_MEANS;
  uint32_t perm_count = model_modifier & MODEL_PERM_COUNT;
  uint32_t fill_orig_chiabs = do_perms || mtest_adjust;
  char* outname_end2 = memcpyb(outname_end, ".qassoc", 8);
  uint32_t perms_total = 0;
  uint32_t pct = 0;
  uint32_t perm_pass_idx = 0;
  uintptr_t perm_vec_ctcl8m = 0;
  int32_t retval = 0;
  double x11 = 0;
  double x12 = 0;
  double x22 = 0;
  uintptr_t* indiv_raw_male_include2 = NULL;
  uint32_t* tcnt = NULL;
  uint32_t mu_table[MODEL_BLOCKSIZE];
  uint32_t uibuf[4];
  char wprefix[5];
  char* wptr;
  char* wptr_restart;
  unsigned char* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* lbptr2;
  uintptr_t* indiv_raw_include2;
  uintptr_t* indiv_include2;
  double* ooptr;
  uint32_t marker_unstopped_ct;
  uint32_t gender_req;
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
  uintptr_t indiv_idx;
  intptr_t g_sum;
  intptr_t nanal;
  double nanal_recip;
  double qt_sum;
  double qt_ssq;
  intptr_t g_ssq;
  double qt_g_prod;
  double qt_mean;
  double g_mean;
  double qt_var;
  double g_var;
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
  int32_t x_code;
  int32_t y_code;
  uint32_t is_reverse;
  uint32_t homrar_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homcom_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t het_code_postxor;
  uintptr_t xor_homcom;
  uintptr_t ulii;
  uint64_t ullii;
  double dxx;
  double dyy;
  double dzz;
  double pval;
  uint32_t loop_end;
  char* a1ptr;
  char* a2ptr;
  if (pheno_nm_ct < 2) {
    logprint("Skipping QT --assoc since less than two phenotypes are present.\n");
    goto qassoc_ret_1;
  }
  g_pheno_nm_ct = pheno_nm_ct;
  g_aperm_alpha = aperm_alpha;
  g_reverse = marker_reverse;
  g_perms_done = 0;
  if (perm_maxt) {
    perms_total = model_mperm_val;
    if (wkspace_alloc_d_checked(&g_maxt_extreme_stat, sizeof(double) * perms_total)) {
      goto qassoc_ret_NOMEM;
    }
    fill_double_zero(g_maxt_extreme_stat, perms_total); // square of t-stat
  } else if (perm_adapt) {
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
  if (wkspace_alloc_uc_checked(&loadbuf_raw, unfiltered_indiv_ct4)) {
    goto qassoc_ret_NOMEM;
  }
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
    if (fprintf(outfile_qtm, tbuf, "SNP") < 0) {
      goto qassoc_ret_WRITE_FAIL;
    }
    *outname_end2 = '\0';
  }
  if (species_haploid_mask[chrom_info_ptr->species] & chrom_info_ptr->chrom_mask) {
    logprint("Warning: QT --assoc doesn't handle X/Y/haploid markers normally (try --linear).\n");
  }
  sprintf(logbuf, "Writing QT --assoc report to %s...", outname);
  logprintb();
  fflush(stdout);
  sprintf(tbuf, " CHR %%%us         BP    NMISS       BETA         SE         R2        T            P \n", plink_maxsnp);
  if (fprintf(outfile, tbuf, "SNP") < 0) {
    goto qassoc_ret_WRITE_FAIL;
  }
  if (do_perms) {
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
    g_tot_quotient = 4294967296LLU / pheno_nm_ct;
    magic_num(g_tot_quotient, &g_totq_magic, &g_totq_preshift, &g_totq_postshift, &g_totq_incr);
    if (wkspace_alloc_ui_checked(&g_missing_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_het_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_homcom_cts, marker_ct * sizeof(uint32_t)) ||
        wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_ct * sizeof(uint32_t))) {
      goto qassoc_ret_NOMEM;
    }
    fill_uint_zero(g_perm_2success_ct, marker_ct);
  }
  g_adaptive_ci_zt = ltqnorm(1 - aperm_beta / (2.0 * marker_ct));
  if (wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * pheno_nm_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&g_orig_1mpval, marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&g_marker_uidxs, marker_ct * sizeof(uint32_t)) ||
      wkspace_alloc_ul_checked(&indiv_raw_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
    goto qassoc_ret_NOMEM;
  }
  memset(wprefix, 32, 5);
  exclude_to_vec_include(unfiltered_indiv_ct, indiv_raw_include2, pheno_nm);
  fill_vec_55(indiv_include2, pheno_nm_ct);
  x_code = species_x_code[chrom_info_ptr->species];
  y_code = species_y_code[chrom_info_ptr->species];
  ullii = chrom_info_ptr->chrom_mask;
  gender_req = ((x_code != -1) && (ullii & (1LLU << x_code))) || ((y_code != -1) && (ullii & (1LLU << y_code)));
  if (gender_req) {
    if (wkspace_alloc_ul_checked(&indiv_raw_male_include2, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
      goto qassoc_ret_NOMEM;
    }
  }
  marker_unstopped_ct = marker_ct;
  indiv_uidx = 0;
  if (wkspace_alloc_d_checked(&g_pheno_d2, pheno_nm_ct * sizeof(double))) {
    goto qassoc_ret_NOMEM;
  }
  g_pheno_sum = 0;
  g_pheno_ssq = 0;
  for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
    indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
    dxx = pheno_d[indiv_uidx++];
    g_pheno_d2[indiv_idx] = dxx;
    g_pheno_sum += dxx;
    g_pheno_ssq += dxx * dxx;
  }
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
    //   adaptive:
    //     g_thread_git_qbufs: (8 * perm_vec_ct, CL-aligned) * 3 * thread_ct
    //   max(T):
    //     g_qresultbuf: MODEL_BLOCKSIZE * (8 * perm_vec_ct, CL-aligned) * 3
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
    } else if (!g_perm_vec_ct) {
      goto qassoc_ret_NOMEM;
    }
    perm_vec_ctcl8m = (g_perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
    if (wkspace_alloc_d_checked(&g_perm_vecstd, perm_vec_ctcl8m * sizeof(double) * pheno_nm_ct)) {
      goto qassoc_ret_NOMEM;
    }
    if (perm_maxt) {
      if (wkspace_alloc_d_checked(&g_maxt_thread_results, g_thread_ct * perm_vec_ctcl8m * sizeof(double)) ||
	  wkspace_alloc_d_checked(&g_qresultbuf, 3 * MODEL_BLOCKSIZE * perm_vec_ctcl8m * sizeof(double))) {
	goto qassoc_ret_NOMEM;
      }
    } else {
      if (wkspace_alloc_d_checked(&g_thread_git_qbufs, perm_vec_ctcl8m * sizeof(double) * 3 * g_thread_ct)) {
	goto qassoc_ret_NOMEM;
      }
      fill_double_zero(g_thread_git_qbufs, 3 * g_thread_ct * perm_vec_ctcl8m);
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
    if (spawn_threads(threads, &qassoc_gen_perms_thread, g_assoc_thread_ct)) {
      logprint(errstr_thread_create);
      goto qassoc_ret_THREAD_CREATE_FAIL;
    }
    ulii = 0;
    qassoc_gen_perms_thread((void*)ulii);
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
      refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_haploid);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      g_is_y = (uii == (uint32_t)y_code);
      intprint2(&(wprefix[2]), uii);
    } else if (perm_maxt) {
      memcpy(g_loadbuf, &(g_loadbuf[(MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * pheno_nm_ctl2]), MODEL_BLOCKKEEP * pheno_nm_ctl2 * sizeof(intptr_t));
      memcpy(g_qresultbuf, &(g_qresultbuf[3 * (MODEL_BLOCKSIZE - MODEL_BLOCKKEEP) * perm_vec_ctcl8m]), MODEL_BLOCKKEEP * perm_vec_ctcl8m * 3 * sizeof(double));
      g_qblock_start = MODEL_BLOCKKEEP;
    } else {
      g_qblock_start = 0;
    }
    block_size = g_qblock_start;
    block_end = marker_unstopped_ct + g_qblock_start - marker_idx;
    if (block_end > MODEL_BLOCKSIZE) {
      block_end = MODEL_BLOCKSIZE;
    }
    do {
      if (perm_adapt && g_perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto qassoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      if (fread(loadbuf_raw, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto qassoc_ret_READ_FAIL;
      }
      if (g_is_haploid) {
	if (g_is_x) {
	  if (xmhh_exists) {
	    hh_reset(loadbuf_raw, indiv_raw_male_include2, unfiltered_indiv_ct);
	  }
	} else if (g_is_y) {
	  hh_reset_y(loadbuf_raw, indiv_raw_include2, indiv_raw_male_include2, unfiltered_indiv_ct);
	} else if (nxmhh_exists) {
	  hh_reset(loadbuf_raw, indiv_raw_include2, unfiltered_indiv_ct);
	}
      }
      indiv_uidx = 0;
      ukk = pheno_nm_ct / BITCT2;
      loadbuf_ptr = &(g_loadbuf[block_size * pheno_nm_ctl2]);
      for (uii = 0; uii < ukk; uii++) {
	ulii = 0;
	for (ujj = 0; ujj < BITCT; ujj += 2) {
	  indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	  ulii |= ((uintptr_t)(((loadbuf_raw[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3))) << ujj;
	  indiv_uidx++;
	}
	*loadbuf_ptr++ = ulii;
      }
      ujj = 2 * (pheno_nm_ct & (BITCT2 - 1));
      if (ujj) {
	ulii = 0;
	for (uii = 0; uii < ujj; uii += 2) {
	  indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	  ulii |= ((uintptr_t)(((loadbuf_raw[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3))) << uii;
	  indiv_uidx++;
	}
	*loadbuf_ptr = ulii;
      }
      if (perm_adapt) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      if (is_set(marker_exclude, ++marker_uidx)) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
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
        is_reverse = is_set(marker_reverse, marker_uidx2);
	loadbuf_ptr = &(g_loadbuf[marker_bidx * pheno_nm_ctl2]);
	vec_3freq(pheno_nm_ctl2, loadbuf_ptr, indiv_include2, &missing_ct, &het_ct, &homcom_ct);
	nanal = pheno_nm_ct - missing_ct;
	memcpy(tbuf, wprefix, 5);
        wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), &(tbuf[5]));
	*wptr++ = ' ';
	wptr = uint32_writew10(wptr, marker_pos[marker_uidx2]);
	*wptr++ = ' ';
	wptr = uint32_writew8(wptr, nanal);
	*wptr++ = ' ';
	if (!is_reverse) {
	  homrar_ct = nanal - het_ct - homcom_ct;
	  xor_homcom = ~ZEROLU;
	  het_code_postxor = 1;
	} else {
	  homrar_ct = homcom_ct;
	  homcom_ct = nanal - het_ct - homcom_ct;
	  xor_homcom = 0;
	  het_code_postxor = 2;
	}
	if (do_perms) {
	  g_missing_cts[marker_idx + marker_bidx] = missing_ct;
	  g_homcom_cts[marker_idx + marker_bidx] = homcom_ct;
	  g_het_cts[marker_idx + marker_bidx] = het_ct;
	}
	g_sum = 2 * homrar_ct + het_ct;
	g_ssq = 4 * homrar_ct + het_ct;
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
	  ulii = (*lbptr2++) ^ xor_homcom;
	  if (uii + BITCT2 > pheno_nm_ct) {
	    ulii &= (ONELU << ((pheno_nm_ct & (BITCT2 - 1)) * 2)) - ONELU;
	  }
	  while (ulii) {
	    ujj = CTZLU(ulii) & (BITCT - 2);
	    ukk = (ulii >> ujj) & 3;
	    indiv_idx = uii + (ujj / 2);
	    dxx = g_pheno_d2[indiv_idx];
	    if (ukk == het_code_postxor) {
	      qt_g_prod += dxx;
	      if (qt_means) {
		qt_het_sum += dxx;
		qt_het_ssq += dxx * dxx;
	      }
	    } else if (ukk == 3) {
	      qt_g_prod += 2 * dxx;
	      if (qt_means) {
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
	g_mean = ((double)g_sum) * nanal_recip;
	dxx = 1.0 / ((double)(nanal - 1));
	qt_var = (qt_ssq - qt_sum * qt_mean) * dxx;
	g_var = (((double)g_ssq) - g_sum * g_mean) * dxx;
	qt_g_covar = (qt_g_prod - qt_sum * g_mean) * dxx;

	dxx = 1.0 / g_var;
	beta = qt_g_covar * dxx;
	vbeta_sqrt = sqrt((qt_var * dxx - beta * beta) / ((double)(nanal - 2)));
	tstat = beta / vbeta_sqrt;
	if (fill_orig_chiabs) {
	  g_orig_chisq[marker_idx + marker_bidx] = tstat;
	  if (mtest_adjust) {
	    tcnt[marker_idx + marker_bidx] = nanal;
	  }
	}
	if (nanal > 1) {
	  tp = tprob(tstat, nanal - 2);
	  rsq = (qt_g_covar * qt_g_covar) / (qt_var * g_var);
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
	    memcpy(wptr, " \n", 2);
	    if (fwrite_checked(tbuf, 2 + (wptr - tbuf), outfile)) {
	      goto qassoc_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  wptr = memcpya(wptr, "        NA         NA         NA       NA           NA \n", 56);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto qassoc_ret_WRITE_FAIL;
	  }
	}
	if (qt_means) {
	  wptr_restart = &(tbuf[6 + plink_maxsnp]);
	  wptr = memcpya(wptr_restart, "  GENO ", 7);
	  a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
	  a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
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
	g_maxt_block_base = marker_idx - g_qblock_start;
	g_maxt_cur_extreme_stat = g_maxt_extreme_stat[0];
	for (uii = 1; uii < g_perm_vec_ct; uii++) {
	  dxx = g_maxt_extreme_stat[uii];
	  if (dxx > g_maxt_cur_extreme_stat) {
	    g_maxt_cur_extreme_stat = dxx;
	  }
	}
	if (spawn_threads(threads, &qassoc_maxt_thread, g_assoc_thread_ct)) {
	  logprint(errstr_thread_create);
	  goto qassoc_ret_THREAD_CREATE_FAIL;
	}
	qassoc_maxt_thread((void*)ulii);
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
	if (spawn_threads(threads, &qassoc_adapt_thread, g_assoc_thread_ct)) {
	  logprint(errstr_thread_create);
	  goto qassoc_ret_THREAD_CREATE_FAIL;
	}
	qassoc_adapt_thread((void*)ulii);
        join_threads(threads, g_assoc_thread_ct);
      }
    }
    marker_idx += block_size - g_qblock_start;
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
      retval = multcomp(outname, outname_end, g_marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, g_orig_chisq, pfilter, mtest_adjust, adjust_lambda, 1, tcnt);
      if (retval) {
	goto qassoc_ret_1;
      }
    }
  }
  if (do_perms) {
    wkspace_reset((unsigned char*)g_perm_vecstd);
    if (g_perms_done < perms_total) {
      if (perm_adapt) {
	marker_unstopped_ct = marker_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, 0, (marker_ct + sizeof(uintptr_t) - 1) / sizeof(uintptr_t));
	if (!marker_unstopped_ct) {
	  goto qassoc_adapt_perm_count;
	}
      }
      printf("\r%u permutations complete.", g_perms_done);
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
    sprintf(logbuf, "%u %s permutations complete.\n", g_perms_done, perm_maxt? "max(T)" : "(adaptive)");
    logprintb();

    if (perm_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
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
    /*
    if (perm_maxt) {
      printf("extreme stats: %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[perms_total - 1]);
    }
    */
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto qassoc_ret_WRITE_FAIL;
    }
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_non_set_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    memset(tbuf, 32, 5);
    tbuf[5 + plink_maxsnp] = ' ';
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      intprint2(&(tbuf[2]), uii);
      for (; marker_uidx < chrom_end;) {
	if (perm_adapt) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) * dxx;
	}
        if ((pval <= pfilter) || (pfilter == 1.0)) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), &(tbuf[5]));
	  wptr = &(tbuf[6 + plink_maxsnp]);
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
	      dzz = (int32_t)(perms_total - doublearr_greater_than(g_maxt_extreme_stat, perms_total, g_orig_chisq[marker_idx] * g_orig_chisq[marker_idx] - EPSILON) + 1);
	      if (!perm_count) {
		wptr = double_g_writewx4(wptr, dzz * dyy, 12);
	      } else {
		wptr = double_g_writewx4(wptr, dzz, 12);
	      }
	    }
	  }
	  memcpy(wptr, " \n", 3);
	  if (fwrite_checked(tbuf, 2 + (uintptr_t)(wptr - tbuf), outfile)) {
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
  qassoc_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 qassoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_qtm);
  return retval;
}
