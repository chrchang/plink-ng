#include "wdist_common.h"
#include "wdist_stats.h"

// load markers in blocks to enable PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 256
#define MODEL_BLOCKKEEP 64

void single_marker_cc_freqs(uintptr_t indiv_ct, uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
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
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_ctrl_ab += popcount2_long(loader2 + loader3);
    tot_ctrl_c += popcount2_long(loader2 & (~loader3));
    loader2 = *case_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_case_ab += popcount2_long(loader2 + loader3);
    tot_case_c += popcount2_long(loader2 & (~loader3));
  }
  *ctrl_missingp = tot_ctrl_c;
  *ctrl_setp = tot_ctrl_ab - tot_ctrl_c;
  *case_missingp = tot_case_c;
  *case_setp = tot_case_ab - tot_case_c;
}

void haploid_single_marker_cc_freqs(uintptr_t indiv_ct, uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
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

void single_marker_cc_3freqs(uintptr_t indiv_ct, uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_hom2p, uint32_t* ctrl_hetp, uint32_t* ctrl_missingp, uint32_t* case_hom2p, uint32_t* case_hetp, uint32_t* case_missingp) {
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
  if (pval == 0) {
    fputs("       INF ", outfile);
  } else if (pval < 0) {
    fputs("        NA ", outfile);
  } else {
    fprintf(outfile, "%10.4g ", pval);
  }
}

static inline void adjust_print_log10(FILE* outfile, double pval) {
  if (pval == 0) {
    fputs("       INF ", outfile);
  } else if (pval < 0) {
    fputs("        NA ", outfile);
  } else if (pval < 1) {
    fprintf(outfile, "%10.4g ", -log10(pval));
  } else {
    fputs("         0 ", outfile);
  }
}

int32_t multcomp(char* outname, char* outname_end, uint32_t* marker_uidxs, uint32_t chi_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, double* chi, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t non_chi) {
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t adjust_gc = mtest_adjust & ADJUST_GC;
  uint32_t is_log10 = mtest_adjust & ADJUST_LOG10;
  uint32_t qq_plot = mtest_adjust & ADJUST_QQ;
  FILE* outfile = NULL;
  double pv_holm = 0.0;
  double pv_sidak_sd = 0;
  int32_t retval = 0;
  uint32_t pct;
  double* sp;
  double* schi;
  double* pv_bh;
  double* pv_by;
  uint32_t* new_order;
  uint32_t cur_idx;
  uint32_t uii;
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
  uint32_t loop_end;

  if (wkspace_alloc_d_checked(&sp, chi_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&schi, chi_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&new_order, chi_ct * sizeof(uint32_t))) {
    goto multcomp_ret_NOMEM;
  }
  if (non_chi) {
    uii = 0;
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
    uii = 0;
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
    goto multcomp_ret_1;
  }
  qsort_ext((char*)sp, chi_ct, sizeof(double), double_cmp_deref, (char*)new_order, sizeof(uint32_t));
#ifdef __cplusplus
  std::sort(schi, &(schi[chi_ct]));
#else
  qsort(schi, chi_ct, sizeof(double), double_cmp);
#endif
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
  for (cur_idx = 0; cur_idx < chi_ct; cur_idx++) {
    pv_gc[cur_idx] = chiprob_p(schi[--uii] * lambda, 1);
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
  sprintf(logbuf, "--adjust values (%u marker%s) written to %s.\n", chi_ct, (chi_ct == 1)? "" : "s", outname);
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

void generate_cc_perm_vec(uint32_t tot_ct, uint32_t set_ct, uint32_t lower_bound, uintptr_t* perm_vec) {
  // assumes perm_vec is initially zeroed out, and lower_bound is set to
  // 4294967296 % tot_ct
  uint32_t num_set = 0;
  uint32_t urand;
  uint32_t uii;
  if (set_ct * 2 < tot_ct) {
    fill_ulong_zero(perm_vec, 2 * (tot_ct + (BITCT - 1) / BITCT));
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(&sfmt);
	} while (urand < lower_bound);
	uii = 2 * (urand % tot_ct);
      } while (is_set(perm_vec, uii));
      set_bit_noct(perm_vec, uii);
    }
  } else {
    fill_vec_55(perm_vec, tot_ct);
    set_ct = tot_ct - set_ct;
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(&sfmt);
	} while (urand < lower_bound);
	uii = 2 * (urand % tot_ct);
      } while (!is_set(perm_vec, uii));
      clear_bit_noct(perm_vec, uii);
    }
  }
}

// multithread globals
// static uint32_t g_thread_start[MAX_THREADS_P1];
static double* g_orig_stats;
static uintptr_t* g_loadbuf;
static uintptr_t* g_perm_vecs;
static uint32_t* g_perm_success_ct;
static uint32_t* g_perm_attempt_ct;
static uintptr_t* g_perm_adapt_stop;
// static uintptr_t* g_perm_adapt_stop_update;
static uint32_t* g_adapt_m_table;
static uintptr_t* g_indiv_nonmale_include2;
static uintptr_t* g_indiv_male_include2 = NULL;
// static uint32_t g_assoc_thread_ct;
static uint32_t g_perm_vec_ct;
static uint32_t g_thread_block_ctl;

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  // uint32_t rand_lbound = (uint32_t)(4294967296LLU % ((uint64_t)pheno_nm_ct));
  int32_t retval = 0;
  FILE* outfile = NULL;
  uint32_t model_assoc = model_modifier & MODEL_ASSOC;
  uint32_t model_adapt = model_modifier & MODEL_PERM;
  uint32_t model_maxt = model_modifier & MODEL_MPERM;
  uint32_t model_fisher = model_modifier & MODEL_FISHER;
  uint32_t model_trendonly = model_modifier & MODEL_TRENDONLY;
  uint32_t model_perm_best = !(model_modifier & MODEL_PMASK);
  uint32_t assoc_counts = model_modifier & MODEL_ASSOC_COUNTS;
  uint32_t assoc_p2 = model_modifier & MODEL_ASSOC_P2;
  uint32_t display_ci = (ci_size > 0);
  uint32_t perms_done = 0;
  uint32_t perms_total = 0;
  int32_t x_code = -1;
  int32_t y_code = -1;
  uint32_t case_ct = 0;
  uint32_t male_ct = 0;
  uint32_t nonmale_ct = 0;
  uint32_t ctrl_male_ct = 0;
  uint32_t case_male_ct = 0;
  uint32_t is_y = 0;
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
  uint32_t mu_table[MODEL_BLOCKSIZE];
  char wformat[64];
  uintptr_t marker_ctl;
  unsigned char* wkspace_mark2;
  char* outname_end2;
  uint32_t* marker_uidxs;
  uint32_t gender_req;
  uint32_t ctrl_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_start;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t marker_bidx;
  // uint32_t ld_check_idx;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx2;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  double* osptr;
  double* ooptr;
  unsigned char* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  // uintptr_t* perm_exclude;
  uintptr_t* indiv_ctrl_include2;
  uintptr_t* indiv_case_include2;
  uint32_t load_indiv_ct;
  uint32_t is_x;
  uint32_t is_haploid;
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
  double pval;
  double dyy;
  double dzz;
  double da1;
  double da12;
  double da2;
  double du1;
  double du12;
  double du2;
  double obs_t;
  double obs_1;
  double obs_2;
  double obs_11;
  double obs_12;
  double obs_22;
  double obs_a1;
  double obs_a2;
  double obs_u1;
  double obs_u2;
  double exp_a1;
  double exp_a2;
  double exp_a11;
  double exp_a12;
  double exp_a22;
  double exp_u1;
  double exp_u2;
  double exp_u11;
  double exp_u12;
  double exp_u22;
  double ca_p;
  char* a1ptr;
  char* a2ptr;
  uint32_t loop_end;

  if (assoc_p2) {
    logprint("Error: --assoc p2 not yet supported.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
  }
  if (!pheno_c) {
    logprint("Error: --assoc does not yet support quantitative phenotypes.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
  }
  if (model_maxt) {
    logprint("Error: max(T) permutation is not yet implemented.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
    perms_total = model_mperm_val;
  } else if (model_adapt) {
    logprint("Error: adaptive permutation is currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
    perms_total = aperm_max;
  }
  if (wkspace_alloc_uc_checked(&loadbuf_raw, unfiltered_indiv_ct4)) {
    goto model_assoc_ret_NOMEM;
  }
  if (model_assoc) {
    if (model_fisher) {
      memcpy(outname_end, ".assoc.fisher", 14);
    } else {
      memcpy(outname_end, ".assoc", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "Writing --assoc report to %s...", outname);
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
    if (!model_fisher) {
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
    if (max_marker_allele_len == 1) {
      sprintf(wformat, "     %%%us %%10u    %%c ", plink_maxsnp);
    } else {
      sprintf(wformat, "     %%%us %%10u %%4s ", plink_maxsnp);
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
    memcpy(outname_end, ".model", 7);
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "Writing --model report to %s...", outname);
    logprintb();
    fflush(stdout);
    outname_end2 = &(outname_end[6]);
    if (model_perm_best && perms_total) {
      memcpy(outname_end2, ".best", 6);
      outname_end2 = &(outname_end2[5]);
    } else if ((model_modifier & MODEL_PGEN) && perms_total) {
      memcpy(outname_end2, ".gen", 5);
      outname_end2 = &(outname_end2[4]);
    } else if (model_modifier & MODEL_PDOM) {
      memcpy(outname_end2, ".dom", 5);
      outname_end2 = &(outname_end2[4]);
    } else if (model_modifier & MODEL_PREC) {
      memcpy(outname_end2, ".rec", 5);
      outname_end2 = &(outname_end2[4]);
    } else if (model_modifier & MODEL_PTREND) {
      memcpy(outname_end2, ".trend", 7);
      outname_end2 = &(outname_end2[6]);
    }
    sprintf(tbuf, " CHR %%%us   A1   A2     TEST            AFF          UNAFF ", plink_maxsnp);
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (!model_fisher) {
      if (fputs("       CHISQ   DF ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    } else {
      memcpy(outname_end2, ".fisher", 8);
      outname_end2 = &(outname_end2[7]);
    }
    if (fputs("           P\n", outfile) == EOF) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (max_marker_allele_len == 1) {
      sprintf(wformat, "     %%%us    %%c    %%c  ", plink_maxsnp);
    } else {
      sprintf(wformat, "     %%%us %%4s %%4s  ", plink_maxsnp);
    }
  }
  marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * pheno_nm_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&g_orig_stats, marker_ct * sizeof(double))) {
    goto model_assoc_ret_NOMEM;
  }
  if (model_assoc) {
    if (wkspace_alloc_d_checked(&orig_odds, marker_ct * sizeof(double))) {
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
      if (is_set(sex_nm, indiv_uidx) && is_set(sex_male, indiv_uidx)) {
	set_bit_noct(g_indiv_male_include2, indiv_idx * 2);
	male_ct++;
      }
      indiv_uidx++;
    }
    vec_init_invert(pheno_nm_ct, g_indiv_nonmale_include2, g_indiv_male_include2);
    nonmale_ct = pheno_nm_ct - male_ct;
  }
  if (model_adapt || model_maxt) {
    if (wkspace_alloc_ui_checked(&g_perm_success_ct, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
    fill_uint_zero(g_perm_success_ct, marker_ct);
    g_thread_block_ctl = (MODEL_BLOCKSIZE + (g_thread_ct - 1)) / g_thread_ct;
    if (model_adapt) {
      if (wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_ct * sizeof(uint32_t)) ||
          wkspace_alloc_ul_checked(&g_perm_adapt_stop, marker_ctl * sizeof(uintptr_t))) {
	goto model_assoc_ret_NOMEM;
      }
      fill_uint_zero(g_perm_attempt_ct, marker_ct);
      fill_ulong_zero(g_perm_adapt_stop, marker_ctl);
    }
  }
  wkspace_mark2 = wkspace_base;
  if (wkspace_alloc_ui_checked(&marker_uidxs, marker_ct * sizeof(uint32_t))) {
    goto model_assoc_ret_NOMEM;
  }
  fputs(" 0%", stdout);
  if (pheno_c) {
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
	case_ct++;
      }
      indiv_uidx++;
    }
    vec_init_invert(pheno_nm_ct, indiv_ctrl_include2, indiv_case_include2);
    ctrl_ct = pheno_nm_ct - case_ct;
  }
  if (gender_req && pheno_c) {
    // replace this with use of special X/Y/generic haploid functions soon
    // since they need to be written for permutation testing anyway, argh
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
      if (is_set(sex_nm, indiv_uidx) && is_set(sex_male, indiv_uidx) && is_set(pheno_c, indiv_uidx)) {
	set_bit_noct(indiv_male_case_include2, indiv_idx * 2);
	case_male_ct++;
      }
      indiv_uidx++;
    }
    vec_init_andnot(pheno_nm_ctl2, indiv_male_ctrl_include2, g_indiv_male_include2, indiv_male_case_include2);
    vec_init_andnot(pheno_nm_ctl2, indiv_nonmale_case_include2, indiv_case_include2, indiv_male_case_include2);
    vec_init_andnot(pheno_nm_ctl2, indiv_nonmale_ctrl_include2, indiv_ctrl_include2, indiv_male_ctrl_include2);
    ctrl_male_ct = male_ct - case_male_ct;
    case_nonmale_ct = case_ct - case_male_ct;
    ctrl_nonmale_ct = ctrl_ct - ctrl_male_ct;
  }

  for (uii = 0; uii < MODEL_BLOCKSIZE; uii++) {
    g_loadbuf[(uii + 1) * pheno_nm_ctl2 - 2] = 0;
    g_loadbuf[(uii + 1) * pheno_nm_ctl2 - 1] = 0;
  }
  if (model_assoc || model_maxt) {
    if (pheno_c && (wkspace_left < pheno_nm_ctl2 * sizeof(intptr_t))) {
      goto model_assoc_ret_NOMEM;
    }
  }
  marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
    goto model_assoc_ret_READ_FAIL;
  }
  do {
    if (pheno_c && (model_assoc || model_maxt)) {
      g_perm_vec_ct = wkspace_left / (pheno_nm_ctl2 * sizeof(intptr_t));
      if (g_perm_vec_ct > perms_total - perms_done) {
	g_perm_vec_ct = perms_total - perms_done;
      }
      perms_done += g_perm_vec_ct;
      g_perm_vecs = (uintptr_t*)wkspace_alloc(g_perm_vec_ct * pheno_nm_ctl2 * sizeof(intptr_t));
    }
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = 0;
    marker_idx = 0;
    marker_idx2 = 0;
    chrom_end = 0;
    loop_end = marker_idx / 100LLU;
    do {
      if (marker_uidx >= chrom_end) {
	block_start = 0;
	if (model_assoc) {
	  // exploit overflow
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_y = (uii == (uint32_t)y_code);
	  if (is_haploid && (!is_x)) {
	    if (is_y) {
	      cur_ctrl_include2 = indiv_male_ctrl_include2;
	      cur_case_include2 = indiv_male_case_include2;
	      load_indiv_ct = male_ct;
	      load_case_ct = case_male_ct;
	    } else {
	      cur_ctrl_include2 = indiv_ctrl_include2;
	      cur_case_include2 = indiv_case_include2;
	      load_indiv_ct = pheno_nm_ct;
	      load_case_ct = case_ct;
	    }
	    load_ctrl_ct = load_indiv_ct - load_case_ct;
	  }
	} else {
	  while (1) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	    } while (marker_uidx >= chrom_end);
	    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    is_x = (uii == (uint32_t)x_code);
	    if ((!((species_haploid_mask[chrom_info_ptr->species] >> uii) & 1LLU)) || is_x) {
	      break;
	    }
	    marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
	  }
	  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto model_assoc_ret_READ_FAIL;
	  }
	  if (!is_x) {
	    cur_ctrl_include2 = indiv_ctrl_include2;
	    cur_case_include2 = indiv_case_include2;
	    load_indiv_ct = pheno_nm_ct;
	    load_case_ct = case_ct;
	  } else {
	    cur_ctrl_include2 = indiv_nonmale_ctrl_include2;
	    cur_case_include2 = indiv_nonmale_case_include2;
	    load_indiv_ct = nonmale_ct;
	    load_case_ct = case_nonmale_ct;
	  }
	  load_ctrl_ct = load_indiv_ct - load_case_ct;
	}
	intprint2(&(wformat[2]), uii);
      } else if (model_maxt) {
	// todo: copy LD table, etc.
	// could conditionally do this for adaptive permutation, but I doubt
	// it's worth the trouble
      } else {
	block_start = 0;
      }
      block_size = block_start;
      block_end = marker_ct + block_start - marker_idx;
      if (block_end > MODEL_BLOCKSIZE) {
	block_end = MODEL_BLOCKSIZE;
      }
      do {
	if ((!model_adapt) || (!is_set(g_perm_adapt_stop, marker_idx2))) {
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
	  if (is_set(marker_exclude, ++marker_uidx)) {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	    if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	      goto model_assoc_ret_READ_FAIL;
	    }
	  }
	} else {
	  do {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  } while (is_set(g_perm_adapt_stop, ++marker_idx2) && (marker_uidx < chrom_end));
	  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto model_assoc_ret_READ_FAIL;
	  }
	}
      } while ((block_size < block_end) && (marker_uidx < chrom_end));
      if (!perm_pass_idx) {
	// basic --assoc/--model
	osptr = &(g_orig_stats[marker_idx + block_start]);
	if (pheno_c) {
	  if (model_assoc) {
	    ooptr = &(orig_odds[marker_idx + block_start]);
	    for (marker_bidx = block_start; marker_bidx < block_size; marker_bidx++) {
	      marker_uidx2 = mu_table[marker_bidx];
	      marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	      is_reverse = is_set(marker_reverse, marker_uidx2);
	      if (is_x) {
		if (is_reverse) {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_nonmale_ctrl_include2, indiv_nonmale_case_include2, &ujj, &uii, &umm, &ukk);
		  uii = 2 * (ctrl_nonmale_ct - uii) - ujj;
		  ukk = 2 * (case_nonmale_ct - ukk) - umm;
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_male_ctrl_include2, indiv_male_case_include2, &uoo, &unn, &uqq, &upp);
		  unn = ctrl_male_ct - unn - uoo;
		  upp = case_male_ct - upp - uqq;
		} else {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_nonmale_ctrl_include2, indiv_nonmale_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = 2 * (ctrl_nonmale_ct - ujj) - uii;
		  umm = 2 * (case_nonmale_ct - umm) - ukk;
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_male_ctrl_include2, indiv_male_case_include2, &unn, &uoo, &upp, &uqq);
		  uoo = ctrl_male_ct - uoo - unn;
		  uqq = case_male_ct - uqq - upp;
		}
		uii += unn;
		ujj += uoo;
		ukk += upp;
		umm += uqq;
	      } else if (is_haploid) {
		if (is_reverse) {
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &ujj, &uii, &umm, &ukk);
		  uii = load_ctrl_ct - uii - ujj;
		  ukk = load_case_ct - ukk - umm;
		} else {
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = load_ctrl_ct - ujj - uii;
		  umm = load_case_ct - umm - ukk;
		}
	      } else {
		if (is_reverse) {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &ujj, &uii, &umm, &ukk);
		  uii = 2 * (ctrl_ct - uii) - ujj;
		  ukk = 2 * (case_ct - ukk) - umm;
		} else {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = 2 * (ctrl_ct - ujj) - uii;
		  umm = 2 * (case_ct - umm) - ukk;
		}
	      }
	      da1 = umm;
	      da2 = ukk;
	      du1 = ujj;
	      du2 = uii;
	      if (model_fisher) {
		pval = fisher22(uii, ujj, ukk, umm);
		*osptr = 1 - pval;
	      } else if (!assoc_p2) {
		if ((((uint64_t)umm) * (ujj + uii)) == (((uint64_t)ujj) * (umm + ukk))) {
		  if ((!umm) && (!ujj)) {
		    *osptr = -9;
		    pval = -1;
		  } else {
		    *osptr = 0;
		    pval = 1;
		  }
		} else {
		  dxx = 1.0 / ((double)(da1 + da2 + du1 + du2));
		  dzz = (da1 + du1) * dxx;
		  dvv = (da2 + du2) * dxx;
		  dyy = (da1 + da2) * dzz;
		  dzz *= du1 + du2;
		  dww = (da1 + da2) * dvv;
		  dvv *= du1 + du2;
		  *osptr = (da1 - dyy) * (da1 - dyy) / dyy + (du1 - dzz) * (du1 - dzz) / dzz + (da2 - dww) * (da2 - dww) / dww + (du2 - dvv) * (du2 - dvv) / dvv;
		  pval = chiprob_p(*osptr, 1);
		}
	      } else {
		// todo
	      }
	      *ooptr = (da1 * du2) / (du1 * da2);
	      if ((pval <= pfilter) || (pfilter == 1.0)) {
		a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
		a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
		if (max_marker_allele_len == 1) {
		  if (fprintf(outfile, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), marker_pos[marker_uidx2], *a1ptr) < 0) {
		    goto model_assoc_ret_WRITE_FAIL;
		  }
		  if (umm + ukk) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u ", umm);
		    } else {
		      fprintf(outfile, "%8.4g ", da1 / (da1 + da2));
		    }
		  } else {
		    fputs("      NA ", outfile);
		  }
		  if (ujj + uii) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u    ", ujj);
		    } else {
		      fprintf(outfile, "%8.4g    ", du1 / (du1 + du2));
		    }
		  } else {
		    fputs("      NA    ", outfile);
		  }
		  putc(*a2ptr, outfile);
		} else {
		  if (fprintf(outfile, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), marker_pos[marker_uidx2], a1ptr) < 0) {
		    goto model_assoc_ret_WRITE_FAIL;
		  }
		  if (umm + ukk) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u ", umm);
		    } else {
		      fprintf(outfile, "%8.4g ", da1 / (da1 + da2));
		    }
		  } else {
		    fputs("      NA ", outfile);
		  }
		  if (ujj + uii) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u", ujj);
		    } else {
		      fprintf(outfile, "%8.4g", du1 / (du1 + du2));
		    }
		  } else {
		    fputs("      NA", outfile);
		  }
		  fprintf(outfile, " %4s", a2ptr);
		}
		if (model_fisher) {
		  fprintf(outfile, " %12.4g ", pval);
		} else {
		  if (pval > -1) {
		    fprintf(outfile, " %12.4g %12.4g ", *osptr, pval);
		  } else {
		    fputs("           NA           NA ", outfile);
		  }
		}
		if (du1 * da2 == 0.0) {
		  fputs("          NA ", outfile);
		  if (display_ci) {
		    fputs("          NA           NA           NA ", outfile);
		  }
		} else {
		  fprintf(outfile, "%12.4g ", *ooptr);
		  if (display_ci) {
		    dxx = log(*ooptr);
		    dyy = sqrt(1 / da1 + 1 / da2 + 1 / du1 + 1 / du2);
		    dzz = ci_zt * dyy;
		    dww = exp(dxx - dzz);
		    dvv = exp(dxx + dzz);
		    fprintf(outfile, "%12.4g %12.4g %12.4g ", dyy, dww, dvv);
		  }
		}
		if (putc('\n', outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
	      }
	      osptr++;
	      ooptr++;
	    }
	  } else {
	    for (marker_bidx = block_start; marker_bidx < block_size; marker_bidx++) {
	      marker_uidx2 = mu_table[marker_bidx];
	      marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	      is_reverse = is_set(marker_reverse, marker_uidx2);
	      if (is_reverse) {
		single_marker_cc_3freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &ukk, &ujj, &uii, &uoo, &unn, &umm);
		uii = load_ctrl_ct - uii - ujj - ukk;
		umm = load_case_ct - umm - unn - uoo;
	      } else {
		single_marker_cc_3freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm, &unn, &uoo);
		ukk = load_ctrl_ct - uii - ujj - ukk;
		uoo = load_case_ct - umm - unn - uoo;
	      }
	      da1 = uoo;
	      da12 = unn;
	      da2 = umm;
	      du1 = ukk;
	      du12 = ujj;
	      du2 = uii;
	      // invalid?
	      upp = (uoo < model_cell_ct) || (unn < model_cell_ct) || (umm < model_cell_ct) || (ukk < model_cell_ct) || (ujj < model_cell_ct) || (uii < model_cell_ct);

	      dxx = da1 + da12 + da2; // obs_A
	      dyy = du1 + du12 + du2; // obs_U
	      obs_t = dxx + dyy;
	      dzz = 1.0 / obs_t;
	      obs_11 = da1 + du1;
	      obs_12 = da12 + du12;
	      obs_22 = da2 + du2;
	      obs_1 = 2 * obs_11 + obs_12;
	      obs_2 = 2 * obs_22 + obs_12;
	      dvv = dzz * dxx;
	      exp_a11 = dvv * obs_11;
	      exp_a12 = dvv * obs_12;
	      exp_a22 = dvv * obs_22;
	      dvv = dzz * dyy;
	      exp_u11 = dvv * obs_11;
	      exp_u12 = dvv * obs_12;
	      exp_u22 = dvv * obs_22;

	      a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
	      a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
	      if (max_marker_allele_len == 1) {
		uss = sprintf(tbuf, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), *a1ptr, *a2ptr);
	      } else {
		uss = sprintf(tbuf, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), a1ptr, a2ptr);
	      }
	      if (!model_trendonly) {
		memcpy(&(tbuf[uss]), "   GENO          ", 17);
		uqq = uss + 8;
		urr = sprintf(&(tbuf[uqq + 43]), "%u/%u/%u          ", uoo, unn, umm);
		if (urr < 24) {
		  memcpy(&(tbuf[uqq + 24 - urr]), &(tbuf[uqq + 43]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 43]), urr);
		  uqq += urr - 9;
		}
		urr = sprintf(&(tbuf[uqq + 34]), "%u/%u/%u ", ukk, ujj, uii);
		if (urr < 15) {
		  memcpy(&(tbuf[uqq + 15 - urr]), &(tbuf[uqq + 34]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 34]), urr);
		  uqq += urr;
		}
		if (upp) {
		  gen_p = -9;
		} else {
		  if (model_fisher) {
		    gen_p = fisher23(uii, ujj, ukk, umm, unn, uoo);
		  } else {
		    dvv = ((da1 - exp_a11) * (da1 - exp_a11)) / exp_a11 +
		      ((da12 - exp_a12) * (da12 - exp_a12)) / exp_a12 +
		      ((da2 - exp_a22) * (da2 - exp_a22)) / exp_a22 +
		      ((du1 - exp_u11) * (du1 - exp_u11)) / exp_u11 +
		      ((du12 - exp_u12) * (du12 - exp_u12)) / exp_u12 +
		      ((du2 - exp_u22) * (du2 - exp_u22)) / exp_u22;
		    gen_p = chiprob_p(dvv, 2);
		  }
		}
		if (gen_p < -1) {
		  if (model_fisher) {
		    memcpy(&(tbuf[uqq]), "          NA\n", 14);
		  } else {
		    memcpy(&(tbuf[uqq]), "          NA   NA           NA\n", 32);
		  }
		} else {
		  if (model_fisher) {
		    sprintf(&(tbuf[uqq]), "%12.4g\n", gen_p);
		  } else {
		    sprintf(&(tbuf[uqq]), "%12.4g    2 %12.4g\n", dvv, gen_p);
		  }
		}
		if (fputs(tbuf, outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
	      }
	      memcpy(&(tbuf[uss]), "  TREND            ", 19);
	      uqq = uss + 8;
	      urr = sprintf(&(tbuf[uqq + 34]), "%u/%u            ", uoo * 2 + unn, umm * 2 + unn);
	      if (urr < 26) {
		memcpy(&(tbuf[uqq + 26 - urr]), &(tbuf[uqq + 34]), urr);
		uqq += 15;
	      } else {
		memcpy(&(tbuf[uqq]), &(tbuf[uqq + 34]), urr);
		uqq += urr - 11;
	      }
	      urr = sprintf(&(tbuf[uqq + 24]), "%u/%u ", ukk * 2 + ujj, uii * 2 + ujj);
	      if (urr < 15) {
		memcpy(&(tbuf[uqq + 15 - urr]), &(tbuf[uqq + 24]), urr);
		uqq += 15;
	      } else {
		memcpy(&(tbuf[uqq]), &(tbuf[uqq + 24]), urr);
		uqq += urr;
	      }
	      urr = uqq; // save this for next line
	      if ((((uint64_t)(uoo * 2 + unn)) * (uii * 2 + ujj)) == (((uint64_t)(umm * 2 + unn)) * (ukk * 2 + ujj))) {
		if (uoo * 2 + unn) {
		  ca_chisq = 0;
		  ca_p = 1;
		} else {
		  ca_p = -9;
		}
	      } else {
		// CA
		dww = (dyy * dzz * da12 - dxx * dzz * du12) + 2 * (dyy * dzz * da2 - dxx * dzz * du2);
		// varCA
		dvv = dxx * dyy * (obs_t * (obs_12 + 4 * obs_22) - (obs_12 + 2 * obs_22) * (obs_12 + 2 * obs_22)) * dzz * dzz * dzz;
		ca_chisq = (dww * dww) / dvv;
		ca_p = chiprob_p(ca_chisq, 1);
	      }
	      if (ca_p > -1) {
		if (model_fisher) {
		  sprintf(&(tbuf[uqq]), "%12.4g\n", ca_p);
		} else {
		  sprintf(&(tbuf[uqq]), "%12.4g    1 %12.4g\n", ca_chisq, ca_p);
		}
	      } else {
		if (!model_fisher) {
		  memcpy(&(tbuf[uqq]), "          NA   NA ", 18);
		  uqq += 18;
		}
		memcpy(&(tbuf[uqq]), "          NA\n", 14);
	      }
	      if (fputs(tbuf, outfile) == EOF) {
		goto model_assoc_ret_WRITE_FAIL;
	      }
	      if (!model_trendonly) {
		memcpy(&(tbuf[uss]), "ALLELIC", 7);
		uqq = urr;
		obs_a1 = 2 * da1 + da12;
		obs_a2 = 2 * da2 + da12;
		obs_u1 = 2 * du1 + du12;
		obs_u2 = 2 * du2 + du12;
		if (model_fisher) {
		  mult_p = fisher22(2 * uoo + unn, 2 * umm + unn, 2 * ukk + ujj, 2 * uii + ujj);
		} else {
		  if ((((uint64_t)(uoo * 2 + unn)) * (uii * 2 + ujj)) == (((uint64_t)(umm * 2 + unn)) * (ukk * 2 + ujj))) {
		    if (uoo * 2 + unn) {
		      dww = 0;
		      mult_p = 1;
		    } else {
		      mult_p = -9;
		    }
		  } else {
		    exp_a1 = dxx * obs_1 * dzz;
		    exp_a2 = dxx * obs_2 * dzz;
		    exp_u1 = dyy * obs_1 * dzz;
		    exp_u2 = dyy * obs_2 * dzz;
		    dww = ((obs_a1 - exp_a1) * (obs_a1 - exp_a1)) / exp_a1 + ((obs_a2 - exp_a2) * (obs_a2 - exp_a2)) / exp_a2 + ((obs_u1 - exp_u1) * (obs_u1 - exp_u1)) / exp_u1 + ((obs_u2 - exp_u2) * (obs_u2 - exp_u2)) / exp_u2;
		    mult_p = chiprob_p(dww, 1);
		  }
		}
		if (mult_p > -1) {
		  if (model_fisher) {
		    sprintf(&(tbuf[uqq]), "%12.4g\n", mult_p);
		  } else {
		    sprintf(&(tbuf[uqq]), "%12.4g    1 %12.4g\n", dww, mult_p);
		  }
		} else {
		  if (!model_fisher) {
		    memcpy(&(tbuf[uqq]), "          NA   NA ", 18);
		    uqq += 18;
		  }
		  memcpy(&(tbuf[uqq]), "          NA\n", 14);
		}
		if (fputs(tbuf, outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
		memcpy(&(tbuf[uss]), "    DOM            ", 19);
		uqq = uss + 8;
		urr = sprintf(&(tbuf[uqq + 34]), "%u/%u            ", uoo + unn, umm);
		if (urr < 26) {
		  memcpy(&(tbuf[uqq + 26 - urr]), &(tbuf[uqq + 34]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 34]), urr);
		  uqq += urr - 11;
		}
		urr = sprintf(&(tbuf[uqq + 24]), "%u/%u ", ukk + ujj, uii);
		if (urr < 15) {
		  memcpy(&(tbuf[uqq + 15 - urr]), &(tbuf[uqq + 24]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 24]), urr);
		  uqq += urr;
		}
		if (upp) {
		  dom_p = -9;
		} else {
		  if (model_fisher) {
		    dom_p = fisher22(uoo + unn, umm, ukk + ujj, uii);
		  } else if (((uint64_t)umm) * (ukk + ujj) == ((uint64_t)uii) * (uoo + unn)) {
		    dww = 0;
		    dom_p = 1;
		  } else {
		    dww = exp_a11 + exp_a12;
		    dvv = exp_u11 + exp_u12;
		    dww = (((da1 + da12) - dww) * ((da1 + da12) - dww)) / dww + ((da2 - exp_a22) * (da2 - exp_a22)) / exp_a22 + (((du1 + du12) - dvv) * ((du1 + du12) - dvv)) / dvv + ((du2 - exp_u22) * (du2 - exp_u22)) / exp_u22;
		    dom_p = chiprob_p(dww, 1);
		  }
		}
		if (dom_p < -1) {
		  if (model_fisher) {
		    memcpy(&(tbuf[uqq]), "          NA\n", 14);
		  } else {
		    memcpy(&(tbuf[uqq]), "          NA   NA           NA\n", 32);
		  }
		} else {
		  if (model_fisher) {
		    sprintf(&(tbuf[uqq]), "%12.4g\n", dom_p);
		  } else {
		    sprintf(&(tbuf[uqq]), "%12.4g    1 %12.4g\n", dww, dom_p);
		  }
		}
		if (fputs(tbuf, outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
		memcpy(&(tbuf[uss]), "    REC            ", 19);
		uqq = uss + 8;
		urr = sprintf(&(tbuf[uqq + 34]), "%u/%u            ", uoo, unn + umm);
		if (urr < 26) {
		  memcpy(&(tbuf[uqq + 26 - urr]), &(tbuf[uqq + 34]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 34]), urr);
		  uqq += urr - 11;
		}
		urr = sprintf(&(tbuf[uqq + 24]), "%u/%u ", ukk, ujj + uii);
		if (urr < 15) {
		  memcpy(&(tbuf[uqq + 15 - urr]), &(tbuf[uqq + 24]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 24]), urr);
		  uqq += urr;
		}
		if (upp) {
		  rec_p = -9;
		} else {
		  if (model_fisher) {
		    rec_p = fisher22(uoo, unn + umm, ukk, ujj + uii);
		  } else if (((uint64_t)uoo) * (ujj + uii) == ((uint64_t)ukk) * (unn + umm)) {
		    dww = 0;
		    rec_p = 1;
		  } else {
		    dww = exp_a22 + exp_a12;
		    dvv = exp_u22 + exp_u12;
		    dww = (((da2 + da12) - dww) * ((da2 + da12) - dww)) / dww + ((da1 - exp_a11) * (da1 - exp_a11)) / exp_a11 + (((du2 + du12) - dvv) * ((du2 + du12) - dvv)) / dvv + ((du1 - exp_u11) * (du1 - exp_u11)) / exp_u11;
		    rec_p = chiprob_p(dww, 1);
		  }
		}

		if (rec_p < -1) {
		  if (model_fisher) {
		    memcpy(&(tbuf[uqq]), "          NA\n", 14);
		  } else {
		    memcpy(&(tbuf[uqq]), "          NA   NA           NA\n", 32);
		  }
		} else {
		  if (model_fisher) {
		    sprintf(&(tbuf[uqq]), "%12.4g\n", rec_p);
		  } else {
		    sprintf(&(tbuf[uqq]), "%12.4g    1 %12.4g\n", dww, rec_p);
		  }
		}
		if (fputs(tbuf, outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
	      }
	      if (model_perm_best) {
		dxx = mult_p;
		if (!upp) {
		  if ((dom_p < dxx) && (dom_p >= 0)) {
		    dxx = dom_p;
		  }
		  if ((rec_p < dxx) && (rec_p >= 0)) {
		    dxx = rec_p;
		  }
		}
		dxx = 1 - dxx;
	      } else if (model_modifier & MODEL_PGEN) {
		dxx = (gen_p >= 0)? (1 - gen_p) : -9;
	      } else if (model_modifier & MODEL_PDOM) {
		dxx = (dom_p >= 0)? (1 - dom_p) : -9;
	      } else if (model_modifier & MODEL_PREC) {
		dxx = (rec_p >= 0)? (1 - rec_p) : -9;
	      } else if (model_modifier & MODEL_PTREND) {
		dxx = (ca_p >= 0)? ca_chisq : -9;
	      }
	      *osptr++ = dxx;
	    }
	  }
	}
      }
      if (perms_total) {
	if (model_assoc) {
	  // TODO
	} else {
	  // TODO
	}
      }
      marker_idx += block_size - block_start;
      if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
	if (marker_idx < marker_ct) {
	  if (pct >= 10) {
	    putchar('\b');
	  }
	  pct = (marker_idx * 100LLU) / marker_ct;
	  printf("\b\b%u%%", pct);
	  fflush(stdout);
	  loop_end = ((uint64_t)(pct + 1LLU) * marker_ct) / 100;
	}
      }
    } while (marker_idx < marker_ct);
    if (!perm_pass_idx) {
      if (pct >= 10) {
	putchar('\b');
      }
      fputs("\b\b\b", stdout);
      logprint(" done.\n");
      if (model_adapt || model_maxt) {
        wkspace_reset((unsigned char*)g_perm_vecs);
      }
      if (mtest_adjust) {
	retval = multcomp(outname, outname_end, marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, g_orig_stats, pfilter, mtest_adjust, adjust_lambda, model_fisher);
	if (retval) {
	  goto model_assoc_ret_1;
	}
      }
    }
    if (model_adapt || model_maxt) {
      printf("\r%u permutations complete.", perms_done);
      fflush(stdout);
      wkspace_reset(wkspace_mark2);
    }
  } while (perms_done < perms_total);
  if (model_adapt || model_maxt) {
    putchar('\r');
    sprintf(logbuf, "%u permutations complete.\n", perms_done);
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
  }
 model_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
