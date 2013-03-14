#include "wdist_common.h"
#include "wdist_stats.h"

// load markers in blocks to enable multithreading and, perhaps later,
// PERMORY-style LD exploitation
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

void generate_cc_perm_vec(uint32_t tot_ct, uint32_t set_ct, uint32_t tot_quotient, uint64_t totq_magic, uint32_t totq_preshift, uint32_t totq_postshift, uint32_t totq_incr, uintptr_t* perm_vec) {
  // assumes perm_vec is initially zeroed out, tot_quotient is
  // 4294967296 / tot_ct, and totq_magic/totq_preshift/totq_postshift/totq_incr
  // have been precomputed from magic_num();
  uint32_t num_set = 0;
  uint32_t upper_bound = tot_ct * tot_quotient;
  uint32_t urand;
  uint32_t uii;
  if (set_ct * 2 < tot_ct) {
    fill_ulong_zero(perm_vec, 2 * (tot_ct + (BITCT - 1) / BITCT));
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(&sfmt);
	} while (urand >= upper_bound);
	uii = 2 * ((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift);
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
	} while (urand >= upper_bound);
	uii = 2 * ((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift);
      } while (!is_set(perm_vec, uii));
      clear_bit_noct(perm_vec, uii);
    }
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

// multithread globals
static double* g_orig_stats;
static uintptr_t* g_loadbuf;
static uintptr_t* g_perm_vecs;
// For --assoc:
// g_compare_bounds[4n] and [4n + 1] define the range of case set allele counts
// which define a failure (Dijkstra interval semantics: minimum, then
// (1 + maximum)).  [4n + 2] and [4n + 3] define ties.  (I.e. usually A2
// counts, but if is_reverse is true, then A1.)
//
// For --model perm-gen:
// [3n] = homozygote clear ct, [3n + 1] = het total, [3n + 2] = hom set total
//
// For --model perm-dom and --model perm-rec:
// Same as --assoc, except the relevant variable is case homozygote a2 or case
// homozygote a1 count, respectively.
static uint32_t* g_compare_bounds;
// For --model perm-trend: fabs(CA) threshold
static double* g_compare_d;
// This is *twice* the number of successes, because PLINK counts tie as 0.5.
// (Actually, PLINK randomizes instead of deterministically adding 0.5; this
// randomization just adds noise so we don't replicate it.)
static uint32_t* g_perm_2success_ct;
static uint32_t* g_perm_attempt_ct;
static double* g_maxt_minp;
static uint32_t* g_maxt_thread_results;
static uintptr_t* g_perm_adapt_stop;
static uint32_t g_adapt_m_table[MODEL_BLOCKSIZE];
static uintptr_t* g_indiv_nonmale_include2;
static uintptr_t* g_indiv_male_include2 = NULL;
static uint32_t* g_marker_uidxs;
static uintptr_t* g_reverse;
static uint32_t g_assoc_thread_ct;
static uint32_t g_perm_vec_ct;
static uint32_t g_thread_block_ctl;
static uint32_t g_block_start;
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
static uint32_t g_is_haploid;

THREAD_RET_TYPE assoc_adapt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uint32_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx_offset = g_perms_done - g_perm_vec_ct;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t marker_idx;
  uint32_t pidx;
  uint32_t next_adapt_check;
  uint32_t compare_min;
  uint32_t compare_maxp1;
  uint32_t compare_tiel;
  uint32_t compare_tieh;
  uint32_t case_set_ct;
  uint32_t uii;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    // guaranteed during loading that g_perm_adapt_stop[] is not set
    // yet
    marker_idx = g_adapt_m_table[marker_bidx];
    next_adapt_check = g_first_adapt_check;
    compare_min = g_compare_bounds[4 * marker_idx];
    compare_maxp1 = g_compare_bounds[4 * marker_idx + 1];
    compare_tiel = g_compare_bounds[4 * marker_idx + 2];
    compare_tieh = g_compare_bounds[4 * marker_idx + 3];
    success_2start = g_perm_2success_ct[marker_idx];
    success_2incr = 0;
    for (pidx = 0; pidx < g_perm_vec_ct;) {
      if (g_is_x) {
	// TODO
	case_set_ct = 0;
      } else if (g_is_haploid) {
	case_set_ct = vec_set_freq_haploid(g_pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]));
      } else {
	case_set_ct = vec_set_freq(g_pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]));
	// if (marker_idx < 10) {
	//   printf("%u %u %u %u %u\n", case_set_ct, compare_min, compare_tiel, compare_maxp1, compare_tieh);
	// }
      }
      if (case_set_ct < compare_min) {
	if (case_set_ct == compare_tiel) {
	  success_2incr++;
	} else {
	  success_2incr += 2;
	}
      } else if (case_set_ct >= compare_maxp1) {
	if (case_set_ct == compare_tieh) {
	  success_2incr++;
	} else {
	  success_2incr += 2;
	}
      }
      if (++pidx == next_adapt_check - pidx_offset) {
	uii = success_2start + success_2incr;
	if (uii) {
	  pval = ((double)((int64_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	  dxx = g_adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	  dyy = pval - dxx; // lower bound
	  dzz = pval + dxx; // upper bound
	  if ((dyy > g_aperm_alpha) || (dzz < g_aperm_alpha)) {
	    set_bit_noct(g_perm_adapt_stop, marker_idx);
	    g_perm_attempt_ct[marker_idx] = next_adapt_check;
	    g_perm_2success_ct[marker_idx] = uii;
	    break;
	  }
	}
	next_adapt_check += (int32_t)(g_adaptive_intercept + ((int32_t)next_adapt_check) * g_adaptive_slope);
      }
    }
    g_perm_2success_ct[marker_idx] = success_2start + success_2incr;
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE assoc_maxt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uint32_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t success_2incr;
  uint32_t marker_idx;
  uint32_t pidx;
  uint32_t compare_min;
  uint32_t compare_maxp1;
  uint32_t compare_tiel;
  uint32_t compare_tieh;
  uint32_t case_set_ct;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = g_adapt_m_table[marker_bidx];
    compare_min = g_compare_bounds[4 * marker_idx];
    compare_maxp1 = g_compare_bounds[4 * marker_idx + 1];
    compare_tiel = g_compare_bounds[4 * marker_idx + 2];
    compare_tieh = g_compare_bounds[4 * marker_idx + 3];
    success_2incr = 0;
    for (pidx = 0; pidx < g_perm_vec_ct;) {
      if (g_is_x) {
	// TODO
	case_set_ct = 0;
      } else if (g_is_haploid) {
	case_set_ct = vec_set_freq_haploid(g_pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]));
      } else {
	case_set_ct = vec_set_freq(g_pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]));
      }
      if (case_set_ct < compare_min) {
	if (case_set_ct == compare_tiel) {
	  success_2incr++;
	} else {
	  success_2incr += 2;
	}
      } else if (case_set_ct >= compare_maxp1) {
	if (case_set_ct == compare_tieh) {
	  success_2incr++;
	} else {
	  success_2incr += 2;
	}
      }
    }
    g_perm_2success_ct[marker_idx] += success_2incr;
  }
  THREAD_RETURN;
}

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  int32_t retval = 0;
  FILE* outfile = NULL;
  uint32_t model_assoc = model_modifier & MODEL_ASSOC;
  uint32_t model_adapt = model_modifier & MODEL_PERM;
  uint32_t model_maxt = model_modifier & MODEL_MPERM;
  uint32_t model_perms = model_adapt | model_maxt;
  uint32_t model_fisher = model_modifier & MODEL_FISHER;
  uint32_t model_trendonly = model_modifier & MODEL_TRENDONLY;
  uint32_t model_perm_best = !(model_modifier & MODEL_PMASK);
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
  uint32_t tot_quotient = 0;
  uint32_t mu_table[MODEL_BLOCKSIZE];
  char wprefix[5];
  char wbuf[48];
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  char* wptr;
  char* wptr2;
  char* wptr_mid;
  char* wptr_mid2;
  uintptr_t marker_ctl;
  unsigned char* wkspace_mark2;
  char* outname_end2;
  uint32_t marker_unstopped_ct;
  uint32_t gender_req;
  uint32_t ctrl_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
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
  uint32_t is_invalid;
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
  g_pheno_nm_ct = pheno_nm_ct;
  g_perms_done = 0;
  g_aperm_alpha = aperm_alpha;
  g_reverse = marker_reverse;

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
    if (wkspace_alloc_d_checked(&g_maxt_minp, perms_total * sizeof(double)) ||
        wkspace_alloc_ui_checked(&g_maxt_thread_results, g_thread_ct * MODEL_BLOCKSIZE * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
    for (uii = 0; uii < perms_total; uii++) {
      g_maxt_minp[uii] = 1.0;
    }
  } else if (model_adapt) {
    perms_total = aperm_max;
  }
  if (wkspace_alloc_uc_checked(&loadbuf_raw, unfiltered_indiv_ct4)) {
    goto model_assoc_ret_NOMEM;
  }
  memset(wprefix, 32, 5);
  if (model_assoc) {
    if (model_fisher) {
      outname_end2 = memcpyb(outname_end, ".assoc.fisher", 14);
    } else {
      outname_end2 = memcpyb(outname_end, ".assoc", 7);
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
    if (model_perm_best && perms_total) {
      outname_end2 = memcpyb(outname_end2, ".best", 6);
    } else if ((model_modifier & MODEL_PGEN) && perms_total) {
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
    if (!model_fisher) {
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
  g_adaptive_ci_zt = ltqnorm(1 - aperm_beta / (2.0 * marker_ct));
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
  if (model_perms) {
    if (wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
    if ((!model_assoc) && (!model_perm_best)) {
      if (model_modifier & MODEL_PGEN) {
	if (wkspace_alloc_ui_checked(&g_compare_bounds, 3 * marker_ct * sizeof(uint32_t))) {
	  goto model_assoc_ret_NOMEM;
	}
      } else if (model_modifier & MODEL_PTREND) {
	if (wkspace_alloc_d_checked(&g_compare_d, marker_ct * sizeof(double))) {
	  goto model_assoc_ret_NOMEM;
	}
      } else {
	if (wkspace_alloc_ui_checked(&g_compare_bounds, 4 * marker_ct * sizeof(uint32_t))) {
	  goto model_assoc_ret_NOMEM;
	}
      }
    } else if (model_assoc) {
      if (wkspace_alloc_ui_checked(&g_compare_bounds, 4 * marker_ct * sizeof(uint32_t))) {
	goto model_assoc_ret_NOMEM;
      }
    }
    fill_uint_zero(g_perm_2success_ct, marker_ct);
    g_thread_block_ctl = (MODEL_BLOCKSIZE + (g_thread_ct - 1)) / g_thread_ct;
    if (model_adapt) {
      if (wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_ct * sizeof(uint32_t)) ||
          wkspace_alloc_ul_checked(&g_perm_adapt_stop, marker_ctl * sizeof(uintptr_t))) {
	goto model_assoc_ret_NOMEM;
      }
      for (uii = 0; uii < marker_ct; uii++) {
	g_perm_attempt_ct[uii] = aperm_max;
      }
      fill_ulong_zero(g_perm_adapt_stop, marker_ctl);
    }
    tot_quotient = 4294967296LLU / pheno_nm_ct;
    magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
  }
  wkspace_mark2 = wkspace_base;
  if (wkspace_alloc_ui_checked(&g_marker_uidxs, marker_ct * sizeof(uint32_t))) {
    goto model_assoc_ret_NOMEM;
  }
  fputs(" 0%", stdout);
  g_case_ct = 0;
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
	g_case_ct++;
      }
      indiv_uidx++;
    }
    vec_init_invert(pheno_nm_ct, indiv_ctrl_include2, indiv_case_include2);
    ctrl_ct = pheno_nm_ct - g_case_ct;
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
    case_nonmale_ct = g_case_ct - case_male_ct;
    ctrl_nonmale_ct = ctrl_ct - ctrl_male_ct;
  }

  for (uii = 0; uii < MODEL_BLOCKSIZE; uii++) {
    g_loadbuf[(uii + 1) * pheno_nm_ctl2 - 2] = 0;
    g_loadbuf[(uii + 1) * pheno_nm_ctl2 - 1] = 0;
  }
  if (model_perms) {
    if (pheno_c && (wkspace_left < pheno_nm_ctl2 * sizeof(intptr_t))) {
      goto model_assoc_ret_NOMEM;
    }
  }
  marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
    goto model_assoc_ret_READ_FAIL;
  }
  do {
    marker_unstopped_ct = marker_ct;
    if (pheno_c && model_perms) {
      if (model_adapt) {
	if (perm_pass_idx) {
	  marker_unstopped_ct = marker_ct - popcount_longs(g_perm_adapt_stop, 0, marker_ctl);
	  if (!marker_unstopped_ct) {
	    break;
	  }
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
      g_perm_vec_ct = wkspace_left / (pheno_nm_ctl2 * sizeof(intptr_t));
      if (g_perm_vec_ct > perms_total - g_perms_done) {
	g_perm_vec_ct = perms_total - g_perms_done;
      }
      g_perms_done += g_perm_vec_ct;
      g_perm_vecs = (uintptr_t*)wkspace_alloc(g_perm_vec_ct * pheno_nm_ctl2 * sizeof(intptr_t));
      for (uii = 0; uii < g_perm_vec_ct; uii++) {
        generate_cc_perm_vec(pheno_nm_ct, g_case_ct, tot_quotient, totq_magic, totq_preshift, totq_postshift, totq_incr, &(g_perm_vecs[uii * pheno_nm_ctl2]));
      }
    }
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = 0;
    marker_idx = 0;
    marker_idx2 = 0;
    chrom_end = 0;
    loop_end = marker_idx / 100LLU;
    do {
      if (marker_uidx >= chrom_end) {
	g_block_start = 0;
	if (model_assoc) {
	  // exploit overflow
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_haploid);
	  uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_y = (uii == (uint32_t)y_code);
	  if (g_is_haploid && (!g_is_x)) {
	    if (is_y) {
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
      } else if (model_maxt) {
	// todo: copy LD table, etc.
	// could conditionally do this for adaptive permutation, but I doubt
	// it's worth the trouble
      } else {
	g_block_start = 0;
      }
      block_size = g_block_start;
      block_end = marker_unstopped_ct + g_block_start - marker_idx;
      if (block_end > MODEL_BLOCKSIZE) {
	block_end = MODEL_BLOCKSIZE;
      }
      do {
	if (model_adapt && is_set(g_perm_adapt_stop, marker_idx2)) {
	  do {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	    marker_idx2++;
	  } while ((marker_uidx < chrom_end) && is_set(g_perm_adapt_stop, marker_idx2));
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
	osptr = &(g_orig_stats[marker_idx + g_block_start]);
	if (pheno_c) {
	  if (model_assoc) {
	    ooptr = &(orig_odds[marker_idx + g_block_start]);
	    for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	      marker_uidx2 = mu_table[marker_bidx];
	      g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	      is_reverse = is_set(marker_reverse, marker_uidx2);
	      if (g_is_x) {
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
	      } else if (g_is_haploid) {
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
		  ukk = 2 * (g_case_ct - ukk) - umm;
		} else {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = 2 * (ctrl_ct - ujj) - uii;
		  umm = 2 * (g_case_ct - umm) - ukk;
		}
	      }
	      da1 = umm;
	      da2 = ukk;
	      du1 = ujj;
	      du2 = uii;
	      if (model_fisher) {
		pval = fisher22(uii, ujj, ukk, umm);
		if (model_perms) {
		  ulii = marker_idx + marker_bidx;
		  ulii *= 4;
		  if (is_reverse) {
		    fisher22_precomp_thresh(umm, ukk, ujj, uii, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		  } else {
		    fisher22_precomp_thresh(ukk, umm, uii, ujj, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		  }
		}
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
		if (model_perms) {
		  ulii = marker_idx + marker_bidx;
		  ulii *= 4;
		  if (is_reverse) {
		    chi22_precomp_thresh(umm, ukk, ujj, uii, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		  } else {
		    chi22_precomp_thresh(ukk, umm, uii, ujj, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		  }
		}
	      } else {
		// todo
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
		if (model_fisher) {
		  wptr = double_g_writewx4(wptr, pval, 12);
		} else {
		  if (pval > -1) {
		    wptr = double_g_writewx4(double_g_writewx4x(wptr, *osptr, 12, ' '), pval, 12);
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
	      osptr++;
	      ooptr++;
	    }
	  } else {
	    for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	      marker_uidx2 = mu_table[marker_bidx];
	      g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
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
	      is_invalid = (uoo < model_cell_ct) || (unn < model_cell_ct) || (umm < model_cell_ct) || (ukk < model_cell_ct) || (ujj < model_cell_ct) || (uii < model_cell_ct);

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
		  if (model_perms && (model_modifier & MODEL_PGEN)) {
		    ulii = marker_idx + marker_bidx;
		    ulii *= 3;
		    g_compare_bounds[ulii + 1] = ujj + unn;
		    if (is_reverse) {
		      g_compare_bounds[ulii] = uii + umm;
		      g_compare_bounds[ulii + 2] = ukk + uoo;
		    } else {
		      g_compare_bounds[ulii] = ukk + uoo;
		      g_compare_bounds[ulii + 2] = uii + umm;
		    }
		  }
		}
		if (gen_p < -1) {
		  wptr = model_assoc_tna(model_fisher, wptr);
		} else {
		  if (!model_fisher) {
		    wptr = memcpya(double_g_writewx4(wptr, dvv, 12), "    2 ", 6);
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
	      if ((((uint64_t)(uoo * 2 + unn)) * (uii * 2 + ujj)) == (((uint64_t)(umm * 2 + unn)) * (ukk * 2 + ujj))) {
		dww = 0;
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
	      if (model_perms && (model_modifier & MODEL_PTREND)) {
		// varCA is constant, and obs_A/obs_U/obs_T are constant across
		// all markers.  So we just save fabs(CA).
		g_compare_d[marker_idx + marker_bidx] = fabs(dww);
	      }
	      if (ca_p > -1) {
		if (!model_fisher) {
		  wptr = memcpya(double_g_writewx4(wptr, ca_chisq, 12), "    1 ", 6);
		}
		wptr = double_g_writewx4x(wptr, ca_p, 12, '\n');
	      } else {
		wptr = model_assoc_tna(model_fisher, wptr);
	      }
	      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
		goto model_assoc_ret_WRITE_FAIL;
	      }
	      if (!model_trendonly) {
		memcpy(wptr_mid, "ALLELIC", 7);
		wptr = wptr_mid2;
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
		  if (!model_fisher) {
		    wptr = memcpya(double_g_writewx4(wptr, dww, 12), "    1 ", 6);
		  }
		  wptr = double_g_writewx4x(wptr, mult_p, 12, '\n');
		} else {
		  wptr = model_assoc_tna(model_fisher, wptr);
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
		} else {
		  if (model_fisher) {
		    dom_p = fisher22(uoo + unn, umm, ukk + ujj, uii);
		    if (model_perms && (model_modifier & MODEL_PDOM)) {
		      ulii = marker_idx + marker_bidx;
		      ulii *= 4;
		      fisher22_precomp_thresh(umm, unn + uoo, uii, ujj + ukk, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		    }
		  } else {
		    if (((uint64_t)umm) * (ukk + ujj) == ((uint64_t)uii) * (uoo + unn)) {
		      dww = 0;
		      dom_p = 1;
		    } else {
		      dww = exp_a11 + exp_a12;
		      dvv = exp_u11 + exp_u12;
		      dww = (((da1 + da12) - dww) * ((da1 + da12) - dww)) / dww + ((da2 - exp_a22) * (da2 - exp_a22)) / exp_a22 + (((du1 + du12) - dvv) * ((du1 + du12) - dvv)) / dvv + ((du2 - exp_u22) * (du2 - exp_u22)) / exp_u22;
		      dom_p = chiprob_p(dww, 1);
		    }
		    if (model_perms && (model_modifier & MODEL_PDOM)) {
		      ulii = marker_idx + marker_bidx;
		      ulii *= 4;
		      chi22_precomp_thresh(umm, unn + uoo, uii, ujj + ukk, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		    }
		  }
		}
		if (dom_p < -1) {
		  wptr = model_assoc_tna(model_fisher, wptr);
		} else {
		  if (!model_fisher) {
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
		} else {
		  if (model_fisher) {
		    rec_p = fisher22(uoo, unn + umm, ukk, ujj + uii);
		    if (model_perms && (model_modifier & MODEL_PREC)) {
		      ulii = marker_idx + marker_bidx;
		      ulii *= 4;
		      fisher22_precomp_thresh(uoo, unn + umm, ukk, ujj + uii, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		    }
		  } else {
		    if (((uint64_t)uoo) * (ujj + uii) == ((uint64_t)ukk) * (unn + umm)) {
		      // special-case zero?
		      dww = 0;
		      rec_p = 1;
		    } else {
		      dww = exp_a22 + exp_a12;
		      dvv = exp_u22 + exp_u12;
		      dww = (((da2 + da12) - dww) * ((da2 + da12) - dww)) / dww + ((da1 - exp_a11) * (da1 - exp_a11)) / exp_a11 + (((du2 + du12) - dvv) * ((du2 + du12) - dvv)) / dvv + ((du1 - exp_u11) * (du1 - exp_u11)) / exp_u11;
		      rec_p = chiprob_p(dww, 1);
		    }
		    if (model_perms && (model_modifier & MODEL_PREC)) {
		      ulii = marker_idx + marker_bidx;
		      ulii *= 4;
		      chi22_precomp_thresh(uoo, unn + umm, ukk, ujj + uii, &(g_compare_bounds[ulii]), &(g_compare_bounds[ulii + 1]), &(g_compare_bounds[ulii + 2]), &(g_compare_bounds[ulii + 3]));
		    }
		  }
		}
		if (rec_p < -1) {
		  wptr = model_assoc_tna(model_fisher, wptr);
		} else {
		  if (!model_fisher) {
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
	g_block_diff = block_size - g_block_start;
	if (g_block_diff >= g_thread_ct) {
	  g_assoc_thread_ct = g_thread_ct;
	} else {
	  g_assoc_thread_ct = g_block_diff;
	}
	if (model_adapt) {
	  if (model_assoc) {
	    if (spawn_threads(threads, &assoc_adapt_thread, g_assoc_thread_ct)) {
	      logprint(errstr_thread_create);
	      goto model_assoc_ret_THREAD_CREATE_FAIL;
	    }
	    ulii = 0;
	    assoc_adapt_thread((void*)ulii);
	    join_threads(threads, g_assoc_thread_ct);
	  } else {
	    // TODO
	  }
	} else {
	  if (model_assoc) {
	    if (spawn_threads(threads, &assoc_maxt_thread, g_assoc_thread_ct)) {
	      logprint(errstr_thread_create);
	      goto model_assoc_ret_THREAD_CREATE_FAIL;
	    }
	    ulii = 0;
	    assoc_maxt_thread((void*)ulii);
	    join_threads(threads, g_assoc_thread_ct);
	  } else {
	    // TODO
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
	retval = multcomp(outname, outname_end, g_marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, g_orig_stats, pfilter, mtest_adjust, adjust_lambda, model_fisher);
	if (retval) {
	  goto model_assoc_ret_1;
	}
      }
    }
    if (model_perms) {
      printf("\r%u permutations complete.", g_perms_done);
      fflush(stdout);
      wkspace_reset(wkspace_mark2);
    }
  } while (g_perms_done < perms_total);
  if (model_perms) {
    putchar('\n');
    sprintf(logbuf, "%u permutations complete.\n", g_perms_done);
    logstr(logbuf);
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
      std::sort(g_maxt_minp, &(g_maxt_minp[perms_total]));
#else
      qsort(g_maxt_minp, perms_total, sizeof(double), double_cmp);
#endif
    }
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
    do {
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
      for (; marker_uidx < chrom_end; marker_idx++) {
        marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	if (model_adapt) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) * dxx;
	}
        if ((pval <= pfilter) || (pfilter == 1.0)) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), &(tbuf[5]));
	  wptr = &(tbuf[6 + plink_maxsnp]);
	  if (!assoc_counts) {
	    wptr = double_g_writewx4x(wptr, pval, 12, ' ');
	  } else {
	    wptr = double_g_writewx4x(wptr, ((double)g_perm_2success_ct[marker_idx]) / 2.0, 12, ' ');
	  }
	  if (model_adapt) {
	    wptr = memseta(wptr, 32, 2);
	    wptr = uint32_writew10(wptr, g_perm_attempt_ct[marker_idx]);
	  } else {
	    dzz = ((int32_t)(doublearr_greater_than(g_maxt_minp, perms_total, 1.0 + EPSILON - g_orig_stats[marker_idx]) + 1));
	    wptr = double_g_writewx4(wptr, dzz * dyy, 12);
	  }
	  memcpy(wptr, " \n", 3);
	  if (fwrite_checked(tbuf, 2 + (uintptr_t)(wptr - tbuf), outfile)) {
	    goto model_assoc_ret_WRITE_FAIL;
	  }
	}
	marker_uidx++;
      }
    } while (marker_idx < marker_ct);
    if (fclose_null(&outfile)) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "Report written to %s.\n", outname);
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
