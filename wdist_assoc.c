#include "wdist_common.h"
#include "wdist_stats.h"

// load markers in blocks to enable multithreading and, perhaps later,
// PERMORY-style LD exploitation (though it might be pointless since our
// hom/het counting routines are already very fast?)
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

int32_t multcomp(char* outname, char* outname_end, uint32_t* marker_uidxs, uintptr_t chi_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, double* chi, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t non_chi) {
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
  //
  // Permutations are stored in an interleaved format, to efficiently bridge
  // the gap between 1 bit permutation statuses and 2-bit PLINK genotypes.  The
  // 64-bit build stores the statuses in 16-byte blocks, in the following
  // order (low bit first):
  //   first 16 bytes: 0 64 1 65 2 66 ... 63 127
  //   second 16 bytes: 128 192 129 193 130 194 ... 191 255
  //   third 16 bytes: 256 320 257 ...
  // whereas the 32-bit build uses 4-byte blocks:
  //   first 4 bytes: 0 16 1 17 2 18 ... 15 31
  //   second 4 bytes: 32 48 33 49 34 50 ... 47 63
  //   third 4 bytes: 64 80 65 ...
  // This forces an additional mask/shift-and-mask to be used in the main
  // loops, but that's easily worth the halving of memory traffic.
  //
  // An alternative approach is to reshuffle the genotype data instead, to use
  // pairs of 16-byte blocks where the first block of each pair contains 128
  // original low bits and the second block contains all the original high
  // bits.  Actually, I'd expect that to be asymptotically faster than the
  // current approach.  But the advantage (if it exists at all) is very small
  // because the popcount algorithms spend so much of their time in
  // 2-bit-accumulator mode, so I'm not bothering for now.  (Maybe that could
  // be thrown in if there's ever an AVX-exploiting rewrite a few years down
  // the line.)
  //
  // Note that the interleaving is almost irrelevant during permutation
  // generation: we act as if there is no interleaving, and then just reshuffle
  // the last block.
  uint32_t num_set = 0;
  uint32_t upper_bound = tot_ct * tot_quotient;
  uintptr_t widx;
  uintptr_t wcomp;
#ifdef __LP64__
  uintptr_t wtmp2;
#endif
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  if (set_ct * 2 < tot_ct) {
#ifdef __LP64__
    fill_ulong_zero(perm_vec, 2 * ((tot_ct + 127) / 128));
#else
    fill_ulong_zero(perm_vec, (tot_ct + 31) / 32);
#endif
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand >= upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT;
	wcomp = ONELU << (uii % BITCT);
	pv_val = perm_vec[widx];
      } while (pv_val & wcomp);
      perm_vec[widx] = pv_val | wcomp;
    }
  } else {
#ifdef __LP64__
    widx = 2 * (tot_ct / 128);
    fill_ulong_one(perm_vec, widx);
    if (tot_ct % 128) {
      if (tot_ct & 64) {
	perm_vec[widx] = ~0LLU;
	perm_vec[widx + 1] = (1LLU << (tot_ct % 64)) - 1LLU;
      } else {
	perm_vec[widx] = (1LLU << (tot_ct % 64)) - 1LLU;
	perm_vec[widx + 1] = 0;
      }
    }
#else
    widx = tot_ct / 32;
    fill_ulong_one(perm_vec, widx);
    if (tot_ct % 32) {
      perm_vec[widx] = (ONELU << (widx % 32)) - 1LLU;
    }
#endif
    set_ct = tot_ct - set_ct;
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand >= upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT;
	wcomp = ONELU << (uii % BITCT);
	pv_val = perm_vec[widx];
      } while (!(pv_val & wcomp));
      perm_vec[widx] = pv_val & (~wcomp);
    }
  }
#ifdef __LP64__
  if (tot_ct % 128) {
    widx = 2 * (tot_ct / 128);
    pv_val = (uint32_t)perm_vec[widx];
    wcomp = 0;
    wtmp2 = 0;
    while (pv_val) {
      uii = CTZLU(pv_val);
      wcomp |= ONELU << (uii * 2);
      pv_val &= pv_val - 1;
    }
    pv_val = perm_vec[widx] >> 32;
    while (pv_val) {
      uii = CTZLU(pv_val);
      wtmp2 |= ONELU << (uii * 2);
      pv_val &= pv_val - 1;
    }
    pv_val = (uint32_t)perm_vec[widx + 1];
    while (pv_val) {
      uii = CTZLU(pv_val);
      wcomp |= ONELU << (uii * 2 + 1);
      pv_val &= pv_val - 1;
    }
    pv_val = perm_vec[widx + 1] >> 32;
    while (pv_val) {
      uii = CTZLU(pv_val);
      wtmp2 |= ONELU << (uii * 2 + 1);
      pv_val &= pv_val - 1;
    }
    perm_vec[widx] = wcomp;
    perm_vec[widx + 1] = wtmp2;
  }
#else
  if (tot_ct % 32) {
    widx = tot_ct / 32;
    pv_val = perm_vec[widx] & 0xffff;
    wcomp = 0;
    while (pv_val) {
      uii = CTZLU(pv_val);
      wcomp |= ONELU << (uii * 2);
      pv_val &= pv_val - 1;
    }
    pv_val = perm_vec[widx] >> 16;
    while (pv_val) {
      uii = CTZLU(pv_val);
      wcomp |= ONELU << (uii * 2 + 1);
      pv_val &= pv_val - 1;
    }
    perm_vec[widx] = wcomp;
  }
#endif
}

void transpose_perms(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uintptr_t* perm_vecst) {
  // Transpose permutations so PRESTO/PERMORY-style genotype indexing can work.
  //
  // We used a 16-ply interleaved format, to allow counts up to 65535 while
  // limiting memory traffic; see generate_cc_perm_vec() for discussion.  The
  // index order used here is:
  // 64-bit build:
  //   first 16 bytes: 0 32 64 96 16 48 80 112 8 40 72 104 24 56 88 120 1... 
  //   next 16 bytes: 128 160 192...
  // 32-bit build:
  //   first 4 bytes: 0 8 16 24 4 12 20 28 2 10 18 26 6 14 22 30 1 9 17...
  //   next 4 bytes: 32 40 48...
  //
  // Yes, converting between two interleaved formats is a bit ugly, but hey,
  // this function only needs to be written once.
  uintptr_t indiv_idx = 0;
#ifdef __LP64__
  uintptr_t pheno_nm_ctlv = 2 * ((pheno_nm_ct + 127) / 128);
  uint16_t wbuf[8];
  uint32_t uii;
#else
  uintptr_t pheno_nm_ctlv = (pheno_nm_ct + 31) / 32;
  uint16_t wbuf[2];
#endif
  uint32_t rshift;
  uint32_t wshift;
  uintptr_t* pvptr;
  uint16_t* wbptr;
  uintptr_t perm_idx;
  for (; indiv_idx < pheno_nm_ct; indiv_idx++) {
    perm_idx = 0;
#ifdef __LP64__
    uii = (indiv_idx / 32) & 1;
    pvptr = &(perm_vecs[2 * (indiv_idx / 128) + uii]);
    rshift = ((indiv_idx % 32) * 2) | ((indiv_idx / 64) & 1);
    goto transpose_perms_loop_start;
    do {
      if (!(perm_idx % 8)) {
	if (perm_idx % 128) {
	  wshift = ((perm_idx & 96) >> 5) | ((perm_idx & 16) >> 2) | (perm_idx & 8);
	} else {
	  memcpy(perm_vecst, wbuf, 16);
	  perm_vecst = &(perm_vecst[2]);
	transpose_perms_loop_start:
	  fill_ulong_zero((uintptr_t*)wbuf, 2);
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctlv] >> rshift) & 1) << wshift;
      wbptr++;
    } while (++perm_idx < perm_vec_ct);
    memcpy(perm_vecst, wbuf, 16);
    perm_vecst = &(perm_vecst[2]);
#else
    pvptr = &(perm_vecs[indiv_idx / 32]);
    rshift = ((indiv_idx % 16) * 2) | ((indiv_idx / 16) & 1);
    goto transpose_perms_loop_start;
    do {
      if (!(perm_idx % 2)) {
	if (perm_idx % 32) {
	  wshift = ((perm_idx & 24) >> 3) | (perm_idx & 4) | ((perm_idx & 2) << 2);
	} else {
	  memcpy(perm_vecst, wbuf, 4);
	  perm_vecst++;
	transpose_perms_loop_start:
	  wbuf[0] = 0;
	  wbuf[1] = 0;
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctlv] >> rshift) & 1) << wshift;
      wbptr++;
    } while (++perm_idx < perm_vec_ct);
    memcpy(perm_vecst, wbuf, 4);
    perm_vecst++;
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

void calc_git(uint32_t pheno_nm_ct, uint32_t perm_vec_ct, uint32_t do_reverse, uintptr_t* __restrict__ loadbuf, uintptr_t* perm_vecst, uintptr_t* thread_bufs) {
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
  //    16-bit short ints, zeroed out, and the second loop restarts.
  // This function is never called when a count might exceed 65535 (the short
  // int limit).
  uint32_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct16 = (perm_vec_ct + 15) / 16;
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  const __m128i m4b = {0x1111111111111111LLU, 0x1111111111111111LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i* permsv = (__m128i*)perm_vecst;
  __m128i* gitv[9];
  __m128i* __restrict__ git_merge4; // no conflicts, please
  __m128i* __restrict__ git_merge8;
  __m128i* __restrict__ git_write;
  __m128i* __restrict__ perm_ptr;
  __m128i loader;
  uint32_t cur_cts[3];
  uintptr_t xormask;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t iidx;
  uint32_t pbidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t indiv_type;
  // thread_bufs structure:
  // first perm_ct8 16-byte blocks: homa1 cts
  // next perm_ct8 blocks: missing cts
  // next perm_ct8 blocks: het cts
  // afterward: free workspace
  gitv[0] = (__m128i*)(&(thread_bufs[6 * perm_ct8]));
  gitv[1] = (__m128i*)(&(thread_bufs[6 * perm_ct8 + 2 * perm_ct32]));
  gitv[2] = (__m128i*)(&(thread_bufs[6 * perm_ct8 + 4 * perm_ct32]));
  gitv[3] = (__m128i*)(&(thread_bufs[6 * (perm_ct8 + perm_ct32)]));
  gitv[4] = (__m128i*)(&(thread_bufs[6 * (perm_ct8 + perm_ct32) + 2 * perm_ct16]));
  gitv[5] = (__m128i*)(&(thread_bufs[6 * (perm_ct8 + perm_ct32) + 4 * perm_ct16]));
  if (do_reverse) {
    xormask = 0;
    gitv[6] = (__m128i*)(&(thread_bufs[2 * perm_ct8]));
    gitv[7] = (__m128i*)(&(thread_bufs[4 * perm_ct8]));
  } else {
    xormask = ~0LLU;
    gitv[6] = (__m128i*)(&(thread_bufs[4 * perm_ct8]));
    gitv[7] = (__m128i*)(&(thread_bufs[2 * perm_ct8]));
  }
  gitv[8] = (__m128i*)thread_bufs;
  cur_cts[0] = 0;
  cur_cts[1] = 0;
  cur_cts[2] = 0;
  for (indiv_type = 0; indiv_type < 3; indiv_type++) {
    git_merge4 = gitv[indiv_type];
    for (uii = 0; uii < perm_ct32; uii++) {
      *git_merge4++ = _mm_setzero_si128();
    }
    git_merge8 = gitv[indiv_type + 3];
    for (uii = 0; uii < perm_ct16; uii++) {
      *git_merge8++ = _mm_setzero_si128();
    }
    git_write = gitv[indiv_type + 6];
    for (uii = 0; uii < perm_ct8; uii++) {
      *git_write++ = _mm_setzero_si128();
    }
  }
  for (uii = 0; uii < pheno_nm_ctl2; uii++) {
    ulii = loadbuf[uii] ^ xormask;
    while (ulii) {
      ujj = CTZLU(ulii) & (BITCT - 2);
      uljj = ulii >> ujj;
      indiv_type = ((ulii >> ujj) & 3) - 1;
      git_merge4 = gitv[indiv_type];
      iidx = uii * BITCT2 + ujj / 2;
      if (iidx >= pheno_nm_ct) {
	break;
      }
      perm_ptr = &(permsv[iidx * perm_ct128]);
      for (pbidx = 0; pbidx < perm_ct128; pbidx++) {
	loader = *perm_ptr++;
	git_merge4[0] = _mm_add_epi64(git_merge4[0], _mm_and_si128(loader, m4b));
	git_merge4[1] = _mm_add_epi64(git_merge4[1], _mm_and_si128(_mm_srli_epi64(loader, 1), m4b));
	git_merge4[2] = _mm_add_epi64(git_merge4[2], _mm_and_si128(_mm_srli_epi64(loader, 2), m4b));
	git_merge4[3] = _mm_add_epi64(git_merge4[3], _mm_and_si128(_mm_srli_epi64(loader, 3), m4b));
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
	    git_write[0] = _mm_add_epi64(git_write[0], _mm_and_si128(loader, m8));
	    git_write[1] = _mm_add_epi64(git_write[1], _mm_and_si128(_mm_srli_epi64(loader, 8), m8));
	    git_write = &(git_write[2]);
	    *git_merge8++ = _mm_setzero_si128();
	  }
	}
      }
      ulii = (uljj & (~3LLU)) << ujj;
    }
  }
  for (indiv_type = 0; indiv_type < 3; indiv_type++) {
    uii = cur_cts[indiv_type];
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
	git_write[0] = _mm_add_epi64(git_write[0], _mm_and_si128(loader, m8));
	git_write[1] = _mm_add_epi64(git_write[1], _mm_and_si128(_mm_srli_epi64(loader, 8), m8));
	git_write = &(git_write[2]);
      }
    }
  }
#else
  uint32_t perm_ct32 = (perm_vec_ct + 31) / 32;
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  uint32_t perm_ct4 = (perm_vec_ct + 3) / 4;
  uint32_t perm_ct2 = (perm_vec_ct + 1) / 2;
  uintptr_t* permsv = perm_vecst;
  uintptr_t* gitv[9];
  uintptr_t* git_merge4;
  uintptr_t* git_merge8;
  uintptr_t* git_write;
  uintptr_t* perm_ptr;
  uintptr_t loader;
  uint32_t cur_cts[3];
  uintptr_t xormask;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t iidx;
  uint32_t pbidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t indiv_type;
  // first perm_ct2 4-byte blocks: homa1_cts
  // ...
  gitv[0] = &(thread_bufs[3 * perm_ct2]);
  gitv[1] = &(thread_bufs[3 * perm_ct2 + perm_ct8]);
  gitv[2] = &(thread_bufs[3 * perm_ct2 + 2 * perm_ct8]);
  gitv[3] = &(thread_bufs[3 * (perm_ct2 + perm_ct8)]);
  gitv[4] = &(thread_bufs[3 * (perm_ct2 + perm_ct8) + perm_ct4]);
  gitv[5] = &(thread_bufs[3 * (perm_ct2 + perm_ct8) + 2 * perm_ct4]);
  if (do_reverse) {
    xormask = 0;
    gitv[6] = &(thread_bufs[perm_ct2]);
    gitv[7] = &(thread_bufs[2 * perm_ct2]);
  } else {
    xormask = ~0LLU;
    gitv[6] = &(thread_bufs[2 * perm_ct2]);
    gitv[7] = &(thread_bufs[perm_ct2]);
  }
  gitv[8] = thread_bufs;
  cur_cts[0] = 0;
  cur_cts[1] = 0;
  cur_cts[2] = 0;
  fill_ulong_zero(thread_bufs, 3 * (perm_ct2 + perm_ct4 + perm_ct8));
  for (uii = 0; uii < pheno_nm_ctl2; uii++) {
    ulii = loadbuf[uii] ^ xormask;
    while (ulii) {
      ujj = CTZLU(ulii) & (BITCT - 2);
      uljj = ulii >> ujj;
      indiv_type = ((ulii >> ujj) & 3) - 1;
      git_merge4 = gitv[indiv_type];
      iidx = uii * BITCT2 + ujj / 2;
      if (iidx >= pheno_nm_ct) {
	break;
      }
      perm_ptr = &(permsv[iidx * perm_ct32]);
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
	    git_write[0] += loader & 0x00ff00ff;
	    git_write[1] += (loader >> 8) & 0x00ff00ff;
	    git_write = &(git_write[2]);
	    *git_merge8++ = 0;
	  }
	}
      }
      ulii = (uljj & (~3LLU)) << ujj;
    }
  }
  for (indiv_type = 0; indiv_type < 3; indiv_type++) {
    uii = cur_cts[indiv_type];
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
	git_write[0] += loader & 0x00ff00ff;
	git_write[1] += (loader >> 8) & 0x00ff00ff;
	git_write = &(git_write[2]);
      }
    }
  }
#endif
}

// multithread globals
static double* g_orig_1mpval;
static double* g_orig_chisq;
static uintptr_t* g_loadbuf;
static uintptr_t* g_perm_vecs;

static uintptr_t* g_perm_vecst; // genotype indexing support
static uintptr_t* g_thread_git_cts;
static uint32_t g_git_thresh;

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
//   p-values than the original.  [4n + 2] and [4n + 3] define the
//   interval with less or equally extreme p-values.
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
//   chi22_get_coeffs()/ca_trend_get_coeffs().
//
// ...fix later...
// For --model perm-gen:
// [3n] = homozygote clear ct, [3n + 1] = het total, [3n + 2] = hom set total
//
static uint32_t* g_precomp_ui;
static double* g_precomp_d;

// X-chromosome: number of missing allele observations per marker relative to
//   *all female* case (so all males automatically contribute at least 1)
// elsewhere: number of missing individuals for each marker
static uint32_t* g_missing_cts;

static uint32_t* g_set_cts;
static uint32_t* g_het_cts;

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
static uint32_t* g_marker_uidxs;
static uintptr_t* g_reverse;
static uint32_t g_model_fisher;
static uint32_t g_assoc_thread_ct;
static uint32_t g_perm_vec_ct;
static uint32_t g_thread_block_ctl;
static uint32_t g_maxt_block_base;
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
#ifdef __LP64__
  uintptr_t pheno_nm_ctlv = 2 * ((pheno_nm_ct + 127) / 128);
#else
  uintptr_t pheno_nm_ctlv = (pheno_nm_ct + 31) / 32;
#endif
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_assoc_thread_ct;
  for (; pidx < pmax; pidx++) {
    generate_cc_perm_vec(pheno_nm_ct, case_ct, tot_quotient, totq_magic, totq_preshift, totq_postshift, totq_incr, &(perm_vecs[pidx * pheno_nm_ctlv]), sfmtp);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE assoc_adapt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  uintptr_t pheno_nm_ctlv = 2 * ((pheno_nm_ct + 127) / 128);
#else
  uintptr_t pheno_nm_ctlv = (pheno_nm_ct + 31) / 32;
#endif
  uint32_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done - perm_vec_ct;
  uint32_t model_fisher = g_model_fisher;
  uint32_t is_x = g_is_x;
  uint32_t is_haploid = g_is_haploid;
  uint32_t min_ploidy = 2;
  uint32_t precomp_width = g_precomp_width;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ indiv_nonmale_include2 = g_indiv_nonmale_include2;
  uintptr_t* __restrict__ indiv_male_include2 = g_indiv_male_include2;
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
      if (is_x) {
        ivec_set_freq_x(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctlv]), indiv_nonmale_include2, indiv_male_include2, &case_set_ct, &case_missing_ct);
      } else if (is_haploid) {
	ivec_set_freq_haploid(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctlv]), &case_set_ct, &case_missing_ct);
      } else {
	ivec_set_freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctlv]), &case_set_ct, &case_missing_ct);
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

THREAD_RET_TYPE assoc_maxt_fisher_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uint32_t is_x = g_is_x;
  uint32_t is_haploid = g_is_haploid;
  uint32_t perm_vec_ct = g_perm_vec_ct;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t perm_ct128 = (perm_vec_ct + 127) / 128;
#ifdef __LP64__
  uintptr_t pheno_nm_ctlv = 2 * ((pheno_nm_ct + 127) / 128);
  uint32_t perm_ct8 = (perm_vec_ct + 7) / 8;
  uint16_t* git_homa1_cts = (uint16_t*)(&(g_thread_git_cts[tidx * perm_ct128 * 168]));
  uint16_t* git_missing_cts = (uint16_t*)(&(g_thread_git_cts[tidx * perm_ct128 * 168 + 2 * perm_ct8]));
  uint16_t* git_het_cts = (uint16_t*)(&(g_thread_git_cts[tidx * perm_ct128 * 168 + 4 * perm_ct8]));
#else
  uintptr_t pheno_nm_ctlv = (pheno_nm_ct + 31) / 32;
  uint32_t perm_ct2 = (perm_vec_ct + 1) / 2;
  uint16_t* git_homa1_cts = (uint16_t*)(&(g_thread_git_cts[tidx * perm_ct128 * 336]));
  uint16_t* git_missing_cts = (uint16_t*)(&(g_thread_git_cts[tidx * perm_ct128 * 336 + perm_ct2]));
  uint16_t* git_het_cts = (uint16_t*)(&(g_thread_git_cts[tidx * perm_ct128 * 336 + 2 * perm_ct2]));
#endif
  uintptr_t perm_vec_ctcl8 = (perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  uint32_t min_ploidy = 2;
  uint32_t precomp_width = g_precomp_width;
  uint32_t git_thresh = g_git_thresh;
  uint32_t case_ct = g_case_ct;
  uintptr_t* __restrict__ loadbuf = g_loadbuf;
  uintptr_t* __restrict__ indiv_nonmale_include2 = g_indiv_nonmale_include2;
  uintptr_t* __restrict__ indiv_male_include2 = g_indiv_male_include2;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  uintptr_t* __restrict__ perm_vecst = g_perm_vecst;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uint32_t* __restrict__ precomp_ui = g_precomp_ui;
  uint32_t* __restrict__ precomp_start = g_precomp_start;
  uint32_t* __restrict__ missing_cts = g_missing_cts;
  uint32_t* __restrict__ set_cts = g_set_cts;
  double* __restrict__ precomp_d = g_precomp_d;
  double* __restrict__ orig_1mpval = g_orig_1mpval;
  uint32_t* __restrict__ gpui;
  double* __restrict__ gpd;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t do_git;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double pv_high;
  double pv_low;
  double pval;
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
    gpd = &(precomp_d[3 * precomp_width * marker_bidx]);
    pv_high = 1.0 + EPSILON - orig_1mpval[marker_idx];
    pv_low = 1.0 - EPSILON - orig_1mpval[marker_idx];
    success_2incr = 0;
    // overly conservative because it combines hets and homa1s, but this isn't
    // a big deal
    // do_git = (col2_sum <= git_thresh) && (missing_cts[marker_idx] <= git_thresh) && (!is_x);
    do_git = (!is_x);
    if (do_git) {
      calc_git(pheno_nm_ct, perm_vec_ct, 0, &(loadbuf[marker_bidx * pheno_nm_ctl2]), perm_vecst, (uintptr_t*)git_homa1_cts);
    }
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!do_git) {
	if (is_x) {
	  ivec_set_freq_x(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctlv]), indiv_nonmale_include2, indiv_male_include2, &case_set_ct, &case_missing_ct);
	} else if (is_haploid) {
	  ivec_set_freq_haploid(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctlv]), &case_set_ct, &case_missing_ct);
	} else {
	  ivec_set_freq(pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), &(perm_vecs[pidx * pheno_nm_ctlv]), &case_set_ct, &case_missing_ct);
	}
      } else {
	if (!is_haploid) {
	  case_missing_ct = git_missing_cts[pidx];
	  case_set_ct = row1x_sum - ((uintptr_t)git_het_cts[pidx]) - 2 * (case_missing_ct + ((uint32_t)git_homa1_cts[pidx]));
	} else {
	  case_missing_ct = ((uint32_t)git_missing_cts[pidx]) + ((uint32_t)git_het_cts[pidx]);
	  case_set_ct = row1x_sum - case_missing_ct - ((uintptr_t)git_homa1_cts[pidx]);
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
	if (ujj > gpui[6 * uii + 5]) {
	  ujj = row1x_sum - case_missing_ct * min_ploidy;
	  // pval = fisher22(case_set_ct, ujj - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - ujj);
	  pval = fisher22_tail_pval(ukk, ujj - ukk, col1_sum - ukk, col2_sum + ukk - ujj, gpui[6 * uii + 5], gpd[3 * uii], gpd[3 * uii + 1], gpd[3 * uii + 2], case_set_ct);
	  if (results[pidx] > pval) {
	    results[pidx] = pval;
	  }
	}
      } else {
	uii = row1x_sum - case_missing_ct * min_ploidy;
	pval = fisher22(case_set_ct, uii - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - uii);
        if (pval < pv_low) {
	  success_2incr += 2;
	} else if (pval < pv_high) {
	  success_2incr++;
	}
	if (results[pidx] > pval) {
	  results[pidx] = pval;
	}
      }
    }
    perm_2success_ct[marker_idx++] += success_2incr;
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE assoc_maxt_thread(void* arg) {
  /*
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t perm_vec_ct16 = (g_perm_vec_ct + 15) / 16;
  uintptr_t perm_vec_ct16d = perm_vec_ct16 * 16;
  uintptr_t perm_vec_ctcl8 = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  unsigned char* git_het_cts = &(g_thread_git_cts[3 * tidx * perm_vec_ct16d]);
  unsigned char* git_homa1_cts = &(g_thread_git_cts[(3 * tidx + 1) * perm_vec_ct16d]);
  unsigned char* git_missing_cts = &(g_thread_git_cts[(3 * tidx + 2) * perm_vec_ct16d]);
  uint32_t min_ploidy = 2;
  uint32_t* gpui;
  double* gpd;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t do_git;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double chisq_high;
  double chisq_low;
  double chisq;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - g_perm_vec_ct]), g_perm_vec_ct * sizeof(double));
  if (g_is_haploid) { // includes g_is_x
    min_ploidy = 1;
  }
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (g_orig_1mpval[marker_idx] == -9) {
      g_perm_2success_ct[marker_idx++] += g_perm_vec_ct;
      continue;
    }
    col1_sum = g_set_cts[marker_idx];
    if (g_is_x) {
      row1x_sum = 2 * g_case_ct;
      tot_obs = 2 * g_pheno_nm_ct - g_missing_cts[marker_idx];
    } else {
      row1x_sum = min_ploidy * g_case_ct;
      tot_obs = min_ploidy * (g_pheno_nm_ct - g_missing_cts[marker_idx]);
    }
    col2_sum = tot_obs - col1_sum;
    missing_start = g_precomp_start[marker_bidx];
    gpui = &(g_precomp_ui[6 * g_precomp_width * marker_bidx]);
    gpd = &(g_precomp_d[2 * g_precomp_width * marker_bidx]);
    chisq_high = g_orig_chisq[marker_idx] + EPSILON;
    chisq_low = g_orig_chisq[marker_idx] - EPSILON;
    success_2incr = 0;
    do_git = (!g_is_x) && (col2_sum <= g_git_thresh) && (g_missing_cts[marker_idx] <= g_git_thresh);
    if (do_git) {
      calc_git(pheno_nm_ctl2, 0, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), git_het_cts, git_homa1_cts, git_missing_cts);
    }
    for (pidx = 0; pidx < g_perm_vec_ct; pidx++) {
      if (!do_git) {
	if (g_is_x) {
	  vec_set_freq_x(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_nonmale_include2, g_indiv_male_include2, &case_set_ct, &case_missing_ct);
	} else if (g_is_haploid) {
	  vec_set_freq_haploid(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
	} else {
	  vec_set_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
	}
      } else {
	if (!g_is_haploid) {
	  case_missing_ct = git_missing_cts[pidx];
	  case_set_ct = row1x_sum - ((uintptr_t)git_het_cts[pidx]) - 2 * (case_missing_ct + ((uint32_t)git_homa1_cts[pidx]));
	} else {
	  case_missing_ct = ((uint32_t)git_missing_cts[pidx]) + ((uint32_t)git_het_cts[pidx]);
	  case_set_ct = row1x_sum - case_missing_ct - ((uintptr_t)git_homa1_cts[pidx]);
	}
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	if (ujj > gpui[6 * uii + 5]) {
	  chisq = ((double)((intptr_t)case_set_ct)) - gpd[2 * uii];
	  chisq = chisq * chisq * gpd[2 * uii + 1];
	  if (results[pidx] < chisq) {
	    results[pidx] = chisq;
	  }
	}
      } else {
	chisq = chi22_eval(case_set_ct, row1x_sum - case_missing_ct * min_ploidy, col1_sum, tot_obs);
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
    g_perm_2success_ct[marker_idx++] += success_2incr;
  }
  */
  THREAD_RETURN;
}

THREAD_RET_TYPE model_adapt_domrec_thread(void* arg) {
  /*
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx_offset = g_perms_done - g_perm_vec_ct;
  uint32_t* gpui;
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
    marker_idx = g_adapt_m_table[marker_bidx];
    next_adapt_check = g_first_adapt_check;
    row1x_sum = g_case_ct;
    tot_obs = g_pheno_nm_ct - g_missing_cts[marker_idx];
    if (g_is_model_prec) {
      col2_sum = (g_set_cts[marker_idx] + g_het_cts[marker_idx]) / 2;
      col1_sum = tot_obs - col2_sum;
    } else {
      col1_sum = (g_set_cts[marker_idx] - g_het_cts[marker_idx]) / 2;
      col2_sum = tot_obs - col1_sum;
    }
    missing_start = g_precomp_start[marker_bidx];
    gpui = &(g_precomp_ui[4 * g_precomp_width * marker_bidx]);
    success_2start = g_perm_2success_ct[marker_idx];
    success_2incr = 0;
    if (g_model_fisher) {
      if (g_orig_1mpval[marker_idx] == -9) {
	g_perm_adapt_stop[marker_idx] = 1;
	g_perm_attempt_ct[marker_idx] = 0;
	continue;
      }
      stat_high = g_orig_1mpval[marker_idx] + EPSILON;
      stat_low = g_orig_1mpval[marker_idx] - EPSILON;
    } else {
      if (g_orig_chisq[marker_idx] == -9) {
	g_perm_adapt_stop[marker_idx;
	g_perm_attempt_ct[marker_idx] = 0;
	continue;
      }
      stat_high = g_orig_chisq[marker_idx] + EPSILON;
      stat_low = g_orig_chisq[marker_idx] - EPSILON;
    }
    cur_xor = g_is_model_prec ^ is_set(g_reverse, marker_idx);
    for (pidx = 0; pidx < g_perm_vec_ct;) {
      if (g_is_x) {
	if (cur_xor) {
	  vec_homclear_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_homx_ct, &case_missing_ct);
	} else {
	  vec_homset_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_homx_ct, &case_missing_ct);
	}
      } else {
	if (cur_xor) {
	  vec_homclear_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_homx_ct, &case_missing_ct);
	} else {
	  vec_homset_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_homx_ct, &case_missing_ct);
	}
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	uii = g_case_ct - case_missing_ct;
        if (g_model_fisher) {
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
	  dxx = g_adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	  dyy = pval - dxx; // lower bound
	  dzz = pval + dxx; // upper bound
	  if ((dyy > g_aperm_alpha) || (dzz < g_aperm_alpha)) {
	    set_bit_noct(g_perm_adapt_stop, marker_idx);
	    g_perm_attempt_ct[marker_idx] = next_adapt_check;
	    break;
	  }
	}
	next_adapt_check += (int32_t)(g_adaptive_intercept + ((int32_t)next_adapt_check) * g_adaptive_slope);
      }
    }
    g_perm_2success_ct[marker_idx] += success_2incr;
  }
  */
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_domrec_fisher_thread(void* arg) {
  /*
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t perm_vec_ct16 = (g_perm_vec_ct + 15) / 16;
  uintptr_t perm_vec_ct16d = perm_vec_ct16 * 16;
  uintptr_t perm_vec_ctcl8 = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  double* results = &(g_maxt_thread_results[perm_vec_ctcl8 * CACHELINE_DBL * tidx]);
  unsigned char* git_het_cts = &(g_thread_git_cts[3 * tidx * perm_vec_ct16d]);
  unsigned char* git_homa1_cts = &(g_thread_git_cts[(3 * tidx + 1) * perm_vec_ct16d]);
  unsigned char* git_missing_cts = &(g_thread_git_cts[(3 * tidx + 2) * perm_vec_ct16d]);
  uint32_t min_ploidy = 2;
  uint32_t* gpui;
  double* gpd;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t row1x_sum;
  intptr_t col1_sum;
  intptr_t col2_sum;
  intptr_t tot_obs;
  uint32_t do_git;
  uint32_t cur_xor;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_homx_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double pv_high;
  double pv_low;
  double pval;
  memcpy(results, &(g_maxt_extreme_stat[g_perms_done - g_perm_vec_ct]), g_perm_vec_ct * sizeof(double));
  if (g_is_haploid) { // includes g_is_x
    min_ploidy = 1;
  }
  marker_idx = g_maxt_block_base + marker_bidx;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    if (g_orig_1mpval[marker_idx] == -9) {
      marker_idx++;
      continue;
    }
    tot_obs = g_pheno_nm_ct - g_missing_cts[marker_idx];
    if (g_is_model_prec) {
      col2_sum = (g_set_cts[marker_idx] + g_het_cts[marker_idx]) / 2;
      col1_sum = tot_obs - col2_sum;
      do_git = (col1_sum <= g_git_thresh);
    } else {
      col1_sum = (g_set_cts[marker_idx] - g_het_cts[marker_idx]) / 2;
      col2_sum = tot_obs - col1_sum;
      do_git = (col2_sum <= g_git_thresh);
    }
    do_git = do_git && (!g_is_x) && (g_het_cts[marker_idx] <= g_git_thresh) && (g_missing_cts[marker_idx] <= g_git_thresh);
    cur_xor = g_is_model_prec ^ is_set(g_reverse, marker_idx);
    cur_reverse = is_set(g_reverse, marker_idx);
    missing_start = g_precomp_start[marker_bidx];
    gpui = &(g_precomp_ui[6 * g_precomp_width * marker_bidx]);
    gpd = &(g_precomp_d[3 * g_precomp_width * marker_bidx]);
    pv_high = 1.0 + EPSILON - g_orig_1mpval[marker_idx];
    pv_low = 1.0 - EPSILON - g_orig_1mpval[marker_idx];
    success_2incr = 0;
    if (do_git) {
      calc_git(pheno_nm_ctl2, 0, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), git_het_cts, git_homa1_cts, git_missing_cts);
    }
    for (pidx = 0; pidx < g_perm_vec_ct; pidx++) {
      if (!do_git) {
	if (g_is_x) {
	  if (cur_xor) {
	    vec_homclear_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_homx_ct, &case_missing_ct);
	  } else {
	  }
	} else {
	  if (cur_xor) {
	  } else {
	  }
	  vec_set_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
	}
      } else {
	if (!g_is_haploid) {
	  case_missing_ct = git_missing_cts[pidx];
	  case_set_ct = row1x_sum - ((uintptr_t)git_het_cts[pidx]) - 2 * (case_missing_ct + ((uint32_t)git_homa1_cts[pidx]));
	} else {
	  case_missing_ct = ((uint32_t)git_missing_cts[pidx]) + ((uint32_t)git_het_cts[pidx]);
	  case_set_ct = row1x_sum - case_missing_ct - ((uintptr_t)git_homa1_cts[pidx]);
	}
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	if (ujj > gpui[6 * uii + 5]) {
	  ujj = row1x_sum - case_missing_ct * min_ploidy;
	  // pval = fisher22(case_set_ct, ujj - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - ujj);
	  pval = fisher22_tail_pval(ukk, ujj - ukk, col1_sum - ukk, col2_sum + ukk - ujj, gpui[6 * uii + 5], gpd[3 * uii], gpd[3 * uii + 1], gpd[3 * uii + 2], case_set_ct);
	  if (results[pidx] > pval) {
	    results[pidx] = pval;
	  }
	}
      } else {
	uii = row1x_sum - case_missing_ct * min_ploidy;
	pval = fisher22(case_set_ct, uii - case_set_ct, col1_sum - case_set_ct, col2_sum + case_set_ct - uii);
        if (pval < pv_low) {
	  success_2incr += 2;
	} else if (pval < pv_high) {
	  success_2incr++;
	}
	if (results[pidx] > pval) {
	  results[pidx] = pval;
	}
      }
    }
    g_perm_2success_ct[marker_idx++] += success_2incr;
  }
  */
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_domrec_thread(void* arg) {
  /*
  intptr_t tidx = (intptr_t)arg;
  uint32_t* results = &(g_maxt_thread_results[MODEL_BLOCKSIZE * tidx]);
  uint32_t marker_bceil = g_block_start + g_block_diff;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t perm_vec_ct16 = (g_perm_vec_ct + 15) / 16;
  uint32_t pidx_start = 16 * ((((int64_t)tidx) * perm_vec_ct16) / g_assoc_thread_ct);
  uint32_t pidx_end = 16 * ((((int64_t)tidx + 1) * perm_vec_ct16) / g_assoc_thread_ct);
  double min_chisq = g_maxt_extreme_stat[pidx_start];
  uint32_t uibuf[4];
  uint32_t* gpui;
  uint32_t* tpui;
  double* gpd;
  uintptr_t marker_bidx;
  uintptr_t marker_idx;
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
  double chisq_high;
  double chisq_low;
  double chisq;
  if (tidx == g_assoc_thread_ct - 1) {
    pidx_end = g_perm_vec_ct;
  }
  uii = (g_precomp_width * 2 + CACHELINE_INT32 - 1) / CACHELINE_INT32;
  tpui = &(g_thread_precomp_ui[tidx * uii * CACHELINE_INT32]);
  for (pidx = pidx_start + 1; pidx < pidx_end; pidx++) {
    if (g_maxt_extreme_stat[pidx] < min_chisq) {
      min_chisq = g_maxt_extreme_stat[pidx];
    }
  }
  marker_idx = g_maxt_block_base + g_block_start;
  for (marker_bidx = g_block_start; marker_bidx < marker_bceil; marker_bidx++) {
    if (g_orig_chisq[marker_idx] == -9) {
      marker_idx++;
      continue;
    }
    tot_obs = g_pheno_nm_ct - g_missing_cts[marker_idx];
    if (g_is_model_prec) {
      col2_sum = (g_set_cts[marker_idx] + g_het_cts[marker_idx]) / 2;
      col1_sum = tot_obs - col2_sum;
    } else {
      col1_sum = (g_set_cts[marker_idx] - g_het_cts[marker_idx]) / 2;
      col2_sum = tot_obs - col1_sum;
    }
    missing_start = g_precomp_start[marker_bidx];
    chisq_high = g_orig_chisq[marker_idx] + EPSILON;
    chisq_low = g_orig_chisq[marker_idx] - EPSILON;
    success_2incr = 0;
    gpui = &(g_precomp_ui[4 * g_precomp_width * marker_bidx]);
    gpd = &(g_precomp_d[2 * g_precomp_width * marker_bidx]);
    case_missing_ct = g_case_ct - missing_start; // actually nonmissing ct here
    ujj = case_missing_ct + 1;
    if (ujj > g_precomp_width) {
      ujj = g_precomp_width;
    }
    for (uii = 0; uii < ujj; uii++) {
      chi22_precomp_val_bounds(min_chisq, case_missing_ct--, col1_sum, tot_obs, uibuf, NULL);
      tpui[uii * 2] = uibuf[2];
      tpui[uii * 2 + 1] = uibuf[3] - uibuf[2];
    }
    for (pidx = pidx_start; pidx < pidx_end; pidx++) {
      if (g_is_x) {
	if (g_is_model_prec ^ is_set(g_reverse, marker_idx)) {
	  vec_homclear_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_homx_ct, &case_missing_ct);
	} else {
	  vec_homset_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_homx_ct, &case_missing_ct);
	}
      } else {
	if (g_is_model_prec ^ is_set(g_reverse, marker_idx)) {
	  vec_homclear_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_homx_ct, &case_missing_ct);
	} else {
	  vec_homset_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_homx_ct, &case_missing_ct);
	}
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	ujj = (uint32_t)(case_homx_ct - tpui[2 * uii]); // deliberate underflow
	if (ujj >= tpui[2 * uii + 1]) {
	  chisq = ((double)((intptr_t)case_homx_ct)) - gpd[2 * uii];
	  chisq = chisq * chisq * gpd[2 * uii + 1];
	  if (g_maxt_extreme_stat[pidx] < chisq) {
	    g_maxt_extreme_stat[pidx] = chisq;
	  }
	}
      } else {
        chisq = chi22_eval(case_homx_ct, g_case_ct - case_missing_ct, col1_sum, tot_obs);
	if (chisq > chisq_high) {
	  success_2incr += 2;
	} else if (chisq > chisq_low) {
	  success_2incr++;
	}
	if (g_maxt_extreme_stat[pidx] < chisq) {
	  g_maxt_extreme_stat[pidx] = chisq;
	}
      }
    }
    results[marker_bidx] = success_2incr;
    marker_idx++;
  }
  */
  THREAD_RETURN;
}

THREAD_RET_TYPE model_adapt_trend_thread(void* arg) {
  /*
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx_offset = g_perms_done - g_perm_vec_ct;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  intptr_t het_ct;
  intptr_t homa2_ct;
  intptr_t tot_obs;
  uint32_t missing_start;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  double chisq_high;
  double chisq_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = g_adapt_m_table[marker_bidx];
    next_adapt_check = g_first_adapt_check;
    if (g_orig_1mpval[marker_idx] == -9) {
      set_bit_noct(g_perm_adapt_stop, marker_idx);
      g_perm_attempt_ct[marker_idx] = next_adapt_check;
      g_perm_2success_ct[marker_idx] = next_adapt_check;
      continue;
    }
    tot_obs = g_pheno_nm_ct - g_missing_cts[marker_idx];
    het_ct = g_het_cts[marker_idx];
    homa2_ct = (g_set_cts[marker_idx] - het_ct) / 2;
    missing_start = g_precomp_start[marker_bidx];
    gpui = &(g_precomp_ui[4 * g_precomp_width * marker_bidx]);
    success_2start = g_perm_2success_ct[marker_idx];
    success_2incr = 0;
    chisq_high = g_orig_chisq[marker_idx] + EPSILON;
    chisq_low = g_orig_chisq[marker_idx] - EPSILON;
    for (pidx = 0; pidx < g_perm_vec_ct;) {
      if (g_is_x) {
	vec_set_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_set_ct, &case_missing_ct);
      } else {
	vec_set_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
      }
      if (is_set(g_reverse, marker_idx)) {
	case_set_ct = (g_case_ct - case_missing_ct) * 2 - case_set_ct;
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	uii = g_case_ct - case_missing_ct;
	dxx = ca_trend_eval(case_set_ct, uii, het_ct, homa2_ct, tot_obs);
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
	  dxx = g_adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	  dyy = pval - dxx; // lower bound
	  dzz = pval + dxx; // upper bound
	  if ((dyy > g_aperm_alpha) || (dzz < g_aperm_alpha)) {
	    set_bit_noct(g_perm_adapt_stop, marker_idx);
	    g_perm_attempt_ct[marker_idx] = next_adapt_check;
	    break;
	  }
	}
	next_adapt_check += (int32_t)(g_adaptive_intercept + ((int32_t)next_adapt_check) * g_adaptive_slope);
      }
    }
    g_perm_2success_ct[marker_idx] += success_2incr;
  }
  */
  THREAD_RETURN;
}

THREAD_RET_TYPE model_maxt_trend_thread(void* arg) {
  /*
  intptr_t tidx = (intptr_t)arg;
  uint32_t* results = &(g_maxt_thread_results[MODEL_BLOCKSIZE * tidx]);
  uint32_t marker_bceil = g_block_start + g_block_diff;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t perm_vec_ct16 = (g_perm_vec_ct + 15) / 16;
  uint32_t pidx_start = 16 * ((((int64_t)tidx) * perm_vec_ct16) / g_assoc_thread_ct);
  uint32_t pidx_end = 16 * ((((int64_t)tidx + 1) * perm_vec_ct16) / g_assoc_thread_ct);
  double min_chisq = g_maxt_extreme_stat[pidx_start];
  uint32_t uibuf[4];
  uint32_t* gpui;
  uint32_t* tpui;
  double* gpd;
  uintptr_t marker_bidx;
  uintptr_t marker_idx;
  uintptr_t pidx;
  intptr_t het_ct;
  intptr_t homa2_ct;
  intptr_t tot_obs;
  uint32_t success_2incr;
  uint32_t missing_start;
  uint32_t case_set_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  uint32_t ujj;
  double chisq_high;
  double chisq_low;
  double chisq;
  if (tidx == g_assoc_thread_ct - 1) {
    pidx_end = g_perm_vec_ct;
  }
  uii = (g_precomp_width * 2 + CACHELINE_INT32 - 1) / CACHELINE_INT32;
  tpui = &(g_thread_precomp_ui[tidx * uii * CACHELINE_INT32]);
  for (pidx = pidx_start + 1; pidx < pidx_end; pidx++) {
    if (g_maxt_extreme_stat[pidx] < min_chisq) {
      min_chisq = g_maxt_extreme_stat[pidx];
    }
  }
  marker_idx = g_maxt_block_base + g_block_start;
  for (marker_bidx = g_block_start; marker_bidx < marker_bceil; marker_bidx++) {
    if (g_orig_1mpval[marker_idx] == -9) {
      results[marker_bidx] = pidx_end - pidx_start;
      marker_idx++;
      continue;
    }
    tot_obs = g_pheno_nm_ct - g_missing_cts[marker_idx];
    het_ct = g_het_cts[marker_idx];
    homa2_ct = (g_set_cts[marker_idx] - het_ct) / 2;
    missing_start = g_precomp_start[marker_bidx];
    chisq_high = g_orig_chisq[marker_idx] + EPSILON;
    chisq_low = g_orig_chisq[marker_idx] - EPSILON;
    success_2incr = 0;
    gpui = &(g_precomp_ui[4 * g_precomp_width * marker_bidx]);
    gpd = &(g_precomp_d[2 * g_precomp_width * marker_bidx]);
    case_missing_ct = g_case_ct - missing_start; // actually nonmissing ct here
    ujj = case_missing_ct + 1;
    if (ujj > g_precomp_width) {
      ujj = g_precomp_width;
    }
    for (uii = 0; uii < ujj; uii++) {
      ca_trend_precomp_val_bounds(min_chisq, case_missing_ct--, het_ct, homa2_ct, tot_obs, uibuf, NULL);
      tpui[uii * 2] = uibuf[2];
      tpui[uii * 2 + 1] = uibuf[3] - uibuf[2];
    }
    for (pidx = pidx_start; pidx < pidx_end; pidx++) {
      if (g_is_x) {
	vec_set_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_set_ct, &case_missing_ct);
      } else {
	vec_set_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
      }
      if (is_set(g_reverse, marker_idx)) {
	case_set_ct = (g_case_ct - case_missing_ct) * 2 - case_set_ct;
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	ujj = (uint32_t)(case_set_ct - tpui[2 * uii]); // deliberate underflow
	if (ujj >= tpui[2 * uii + 1]) {
	  chisq = ((double)((intptr_t)case_set_ct)) - gpd[2 * uii];
	  chisq = chisq * chisq * gpd[2 * uii + 1];
	  if (g_maxt_extreme_stat[pidx] < chisq) {
	    g_maxt_extreme_stat[pidx] = chisq;
	  }
	}
      } else {
        chisq = ca_trend_eval(case_set_ct, g_case_ct - case_missing_ct, het_ct, homa2_ct, tot_obs);
	if (chisq > chisq_high) {
	  success_2incr += 2;
	} else if (chisq > chisq_low) {
	  success_2incr++;
	}
	if (g_maxt_extreme_stat[pidx] < chisq) {
	  g_maxt_extreme_stat[pidx] = chisq;
	}
      }
    }
    results[marker_bidx] = success_2incr;
    marker_idx++;
  }
  */
  THREAD_RETURN;
}

/*
THREAD_RET_TYPE model_adapt_gen_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t marker_bidx = g_block_start + (((uint64_t)tidx) * g_block_diff) / g_assoc_thread_ct;
  uint32_t marker_bceil = g_block_start + (((uint64_t)tidx + 1) * g_block_diff) / g_assoc_thread_ct;
  uintptr_t pheno_nm_ctl2 = 2 * ((g_pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx_offset = g_perms_done - g_perm_vec_ct;
  uint32_t* gpui;
  uintptr_t marker_idx;
  uintptr_t pidx;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  intptr_t het_ct;
  intptr_t homa2_ct;
  intptr_t tot_obs;
  uint32_t missing_start;
  uint32_t case_hom2_ct;
  uint32_t case_het_ct;
  uint32_t case_missing_ct;
  uint32_t uii;
  double chisq_high;
  double chisq_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = g_adapt_m_table[marker_bidx];
    if (g_orig_1mpval[marker_idx] == -9) {
      set_bit_noct(g_perm_adapt_stop, marker_idx);
      g_perm_attempt_ct[marker_idx] = 0;
      continue;
    }
    next_adapt_check = g_first_adapt_check;
    tot_obs = g_pheno_nm_ct - g_missing_cts[marker_idx];
    het_ct = g_het_cts[marker_idx];
    homa2_ct = (g_set_cts[marker_idx] - het_ct) / 2;
    missing_start = g_precomp_start[marker_bidx];
    gpui = &(g_precomp_ui[4 * g_precomp_width * marker_bidx]);
    success_2start = g_perm_2success_ct[marker_idx];
    success_2incr = 0;
    chisq_high = g_orig_chisq[marker_idx] + EPSILON;
    chisq_low = g_orig_chisq[marker_idx] - EPSILON;
    for (pidx = 0; pidx < g_perm_vec_ct;) {
      if (g_is_x) {
	vec_set_freq_xx(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), g_indiv_male_include2, &case_set_ct, &case_missing_ct);
      } else {
	vec_set_freq(pheno_nm_ctl2, &(g_loadbuf[marker_bidx * pheno_nm_ctl2]), &(g_perm_vecs[pidx * pheno_nm_ctl2]), &case_set_ct, &case_missing_ct);
      }
      if (is_set(g_reverse, marker_idx)) {
	case_set_ct = (g_case_ct - case_missing_ct) * 2 - case_set_ct;
      }
      // deliberate underflow
      uii = (uint32_t)(case_missing_ct - missing_start);
      if (uii < g_precomp_width) {
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
	uii = g_case_ct - case_missing_ct;
	dxx = ca_trend_eval(case_set_ct, uii, het_ct, homa2_ct, tot_obs);
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
	  dxx = g_adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
	  dyy = pval - dxx; // lower bound
	  dzz = pval + dxx; // upper bound
	  if ((dyy > g_aperm_alpha) || (dzz < g_aperm_alpha)) {
	    set_bit_noct(g_perm_adapt_stop, marker_idx);
	    g_perm_attempt_ct[marker_idx] = next_adapt_check;
	    break;
	  }
	}
	next_adapt_check += (int32_t)(g_adaptive_intercept + ((int32_t)next_adapt_check) * g_adaptive_slope);
      }
    }
    g_perm_2success_ct[marker_idx] += success_2incr;
  }
  THREAD_RETURN;
}
*/

void get_model_assoc_precomp_bounds(uint32_t missing_ct, uint32_t* minp, uint32_t* ctp) {
  // Estimate which case missing counts are most common.
  // Expected value = (g_case_ct * missing_ct / g_pheno_nm_ct)
  // If X-chromosome:
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
  if (g_is_x) {
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

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  // number of machine words per permutation
  uintptr_t pheno_nm_ctlv = 2 * ((pheno_nm_ct + 127) / 128);
#else
  uintptr_t pheno_nm_ctlv = (pheno_nm_ct + 31) / 32;
#endif
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
  uint32_t model_fisherx = (model_modifier & MODEL_FISHER) && (!(model_modifier & MODEL_PTREND));
  uint32_t mu_table[MODEL_BLOCKSIZE];
  uint32_t uibuf[4];
  char wprefix[5];
  char wbuf[48];
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
  uint32_t* missp;
  uint32_t* setp;
  uint32_t* hetp;
  double* o1mpptr;
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
  g_git_thresh = pheno_nm_ct / 4;
  if (g_git_thresh > 65535) {
    g_git_thresh = 65535;
  }

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
    if (model_perm_best && perms_total) {
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
  g_adaptive_ci_zt = ltqnorm(1 - aperm_beta / (2.0 * marker_ct));
  marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
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

    g_precomp_width = (1 + (int32_t)(sqrt(pheno_nm_ct) * EXPECTED_MISSING_FREQ * 4.24264));
    if (wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_ct * sizeof(uint32_t))) {
      goto model_assoc_ret_NOMEM;
    }
    if (!model_fisherx) {
      if (wkspace_alloc_d_checked(&g_orig_chisq, marker_ct * sizeof(double))) {
	goto model_assoc_ret_NOMEM;
      }
      fill_double_zero(g_orig_chisq, marker_ct);
    }
    if (model_maxt) {
      if (model_fisherx) {
	if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC))) {
	  if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 6 * MODEL_BLOCKSIZE * sizeof(uint32_t)) ||
	      wkspace_alloc_d_checked(&g_precomp_d, g_precomp_width * 3 * MODEL_BLOCKSIZE * sizeof(double))) {
	    goto model_assoc_ret_NOMEM;
	  }
	} else {
	  // TODO
	}
      } else if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
	if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 6 * MODEL_BLOCKSIZE * sizeof(uint32_t)) ||
	    wkspace_alloc_d_checked(&g_precomp_d, g_precomp_width * 2 * MODEL_BLOCKSIZE * sizeof(double))) {
	  goto model_assoc_ret_NOMEM;
	}
      } else {
	// TODO
      }
    } else if (model_assoc || (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
      if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 4 * MODEL_BLOCKSIZE * sizeof(uint32_t))) {
	goto model_assoc_ret_NOMEM;
      }
    } else {
      /*
      if (model_modifier & MODEL_PGEN) {
	if (wkspace_alloc_ui_checked(&g_precomp_ui, 3 * marker_ct * sizeof(uint32_t))) {
	  goto model_assoc_ret_NOMEM;
	}
      } else {
	if (wkspace_alloc_ui_checked(&g_precomp_ui, g_precomp_width * 3 * MODEL_BLOCKSIZE * sizeof(uint32_t))) {
	  goto model_assoc_ret_NOMEM;
	}
      }
      */
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
  marker_unstopped_ct = marker_ct;
 model_assoc_more_perms:
  if (pheno_c && model_perms) {
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
      g_perm_vec_ct = wkspace_left / (pheno_nm_ctlv * sizeof(intptr_t));
    } else {
      // g_perm_vec_ct memory allocation dependencies:
      //   g_maxt_thread_results: (8 * g_perm_vec_ct, cacheline-aligned) *
      //     g_thread_ct
      //   g_perm_vecst: 16 * ((g_perm_vec_ct + 127) / 128) * pheno_nm_ct
      //   g_thread_git_cts: ((g_perm_vec_ct + 127) / 128) * 1344 * thread_ct
      //   g_perm_vecs: pheno_nm_ctlv * g_perm_vec_ct
      // If we force g_perm_vec_ct to be a multiple of 128, then we have
      //   g_perm_vec_ct * (11 * g_thread_ct + g_pheno_nm_ct + pheno_nm_ctl2)
      //
      // "1344?  2368?  wtf?"
      // Excellent questions.  Each max(T) thread has nine buffers to support
      // rapid execution of the genotype indexing algorithm:
      //   three with 4-bit accumulators, total size perm_vec_ct / 2 bytes
      //   three with 8-bit accumulators, total size perm_vec_ct bytes
      //   three with 16-bit accumulators, total size 2 * perm_vec_ct bytes
      // The initial 3 multiplier is to allow heterozygotes, homozygote minors,
      // and missing genotypes to be counted simultaneously.
      // Adding all this up, we have 10.5 * perm_vec_ct bytes, and multiplying
      // by 128 yields 1344.  The other thread_ct dependence contributes
      // 8 * perm_vec_ct bytes, multiplying by 128 yields 1024, and
      // 1344 + 1024 = 2368.
      g_perm_vec_ct = 128 * (wkspace_left / (128LL * pheno_nm_ctlv + 2368LL * g_thread_ct + 16LL * pheno_nm_ct));
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    } else if (!g_perm_vec_ct) {
      goto model_assoc_ret_NOMEM;
    }
    g_perms_done += g_perm_vec_ct;
    g_perm_vecs = (uintptr_t*)wkspace_alloc(g_perm_vec_ct * pheno_nm_ctlv * sizeof(intptr_t));
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
      ulii = ((g_perm_vec_ct + 127) / 128) * 16;
      g_perm_vecst = (uintptr_t*)wkspace_alloc(ulii * pheno_nm_ct);
      g_thread_git_cts = (uintptr_t*)wkspace_alloc(ulii * 84 * g_thread_ct);
      transpose_perms(g_perm_vecs, g_perm_vec_ct, pheno_nm_ct, g_perm_vecst);
      memset(g_thread_git_cts, 0, ulii * 84 * g_thread_ct);
    }
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
    goto model_assoc_ret_READ_FAIL;
  }
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
      g_block_start = 0;
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
      if (pheno_c) {
	if (model_assoc) {
	  ooptr = &(orig_odds[marker_idx + g_block_start]);
	  for (marker_bidx = g_block_start; marker_bidx < block_size; marker_bidx++) {
	    marker_uidx2 = mu_table[marker_bidx];
	    g_marker_uidxs[marker_idx + marker_bidx] = marker_uidx2;
	    is_reverse = is_set(marker_reverse, marker_uidx2);
	    if (g_is_x) {
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
	    } else {
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
	      } else {
		dxx = chi22_eval(ukk, ukk + umm, uii + ukk, uii + ujj + ukk + umm);
		pval = chiprob_p(dxx, 1);
		*o1mpptr = 1 - pval;
	      }
	      if (model_perms) {
		g_orig_chisq[marker_idx + marker_bidx] = dxx;
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
	    da1 = uoo;
	    da12 = unn;
	    da2 = umm;
	    du1 = ukk;
	    du12 = ujj;
	    du2 = uii;
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
	      } else {
		if (g_model_fisher) {
		  gen_p = fisher23(uii, ujj, ukk, umm, unn, uoo);
		} else {
		  chi23_evalx(uii, ujj, ukk, umm, unn, uoo, &dvv, &upp);
		  gen_p = chiprob_px(dvv, upp);
		  if (model_perms && (model_modifier & MODEL_PGEN)) {
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
	    if (model_perms && (model_modifier & MODEL_PTREND)) {
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
		if (model_perms && (!g_model_fisher) && (model_modifier & MODEL_PDOM)) {
		  g_orig_chisq[marker_idx + marker_bidx] = -9;
		}
	      } else {
		if (g_model_fisher) {
		  dom_p = fisher22(uoo + unn, umm, ukk + ujj, uii);
		} else {
		  dww = chi22_evalx(uoo + unn, uoo + unn + umm, uoo + unn + ukk + ujj, uoo + unn + umm + ukk + ujj + uii);
		  dom_p = chiprob_px(dww, 1);
		  if (model_perms && (model_modifier & MODEL_PDOM)) {
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
		if (model_perms && (!g_model_fisher) && (model_modifier & MODEL_PREC)) {
		  g_orig_chisq[marker_idx + marker_bidx] = -9;
		}
	      } else {
		if (g_model_fisher) {
		  rec_p = fisher22(uoo, unn + umm, ukk, ujj + uii);
		} else {
		  dww = chi22_evalx(uoo, uoo + unn + umm, uoo + ukk, uoo + unn + umm + ukk + ujj + uii);
		  rec_p = chiprob_px(dww, 1);
		  if (model_perms && (model_modifier & MODEL_PREC)) {
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
	    missp++;
	    setp++;
	    hetp++;
	    *o1mpptr++ = dxx;
	  }
	}
      }
    }
    if (perms_total) {
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
	  for (uii = 0; uii < g_perm_vec_ct; uii++) {
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
	  upp = g_missing_cts[marker_idx + uii];
	  get_model_assoc_precomp_bounds(upp, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  uoo = g_set_cts[marker_idx + uii];
	  if (g_is_x) {
	    unn = 2 * g_case_ct;
	    upp = 2 * g_pheno_nm_ct - upp;
	  } else {
	    unn = uqq * g_case_ct;
	    upp = uqq * (g_pheno_nm_ct - upp);
	  }
	  ujj *= uqq;
	  ukk += uii * g_precomp_width;
	  if (g_model_fisher) {
	    dxx = 1 - g_orig_1mpval[marker_idx + uii];
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
	    dxx = g_orig_chisq[marker_idx + uii];
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
      } else if (model_modifier & MODEL_PGEN) {
      } else if (model_modifier & MODEL_PTREND) {
	for (uii = g_block_start; uii < block_size; uii++) {
	  upp = g_missing_cts[marker_idx + uii];
	  get_model_assoc_precomp_bounds(upp, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  unn = g_het_cts[marker_idx + uii];
	  upp = g_pheno_nm_ct - upp; // tot_obs
	  uoo = (g_set_cts[marker_idx + uii] - unn) / 2; // homa2_ct
	  ukk += uii * g_precomp_width;
	  ujj = g_case_ct - ujj;
	  dxx = g_orig_chisq[marker_idx + uii];
	  if (model_adapt) {
	    for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	      ca_trend_precomp_val_bounds(dxx, ujj--, unn, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
	    }
	  } else {
	    for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	      ca_trend_precomp_val_bounds(dxx, ujj, unn, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
              ca_trend_precomp_val_bounds(g_maxt_cur_extreme_stat, ujj--, unn, uoo, upp, uibuf, &(g_precomp_d[umm * 2]));
	      g_precomp_ui[umm * 6 + 4] = uibuf[2];
	      g_precomp_ui[umm * 6 + 5] = uibuf[3];
	    }
	  }
	}
      } else {
	for (uii = g_block_start; uii < block_size; uii++) {
	  upp = g_missing_cts[marker_idx + uii];
	  get_model_assoc_precomp_bounds(upp, &ujj, &ukk);
	  g_precomp_start[uii] = ujj;
	  upp = g_pheno_nm_ct - upp; // tot_obs
	  if (model_modifier & MODEL_PREC) {
	    uoo = upp - ((g_set_cts[marker_idx + uii] + g_het_cts[marker_idx + uii]) / 2); // col1_sum
	  } else {
	    uoo = (g_set_cts[marker_idx + uii] - g_het_cts[marker_idx + uii]) / 2;
	  }
	  ukk += uii * g_precomp_width;
	  ujj = g_case_ct - ujj;
	  if (g_model_fisher) {
	    dxx = 1 - g_orig_1mpval[marker_idx + uii];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, ujj--, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
	        fisher22_precomp_pval_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
	        fisher22_precomp_pval_bounds(g_maxt_cur_extreme_stat, ujj--, uoo, upp, uibuf, &(g_precomp_d[umm * 3]));
		g_precomp_ui[umm * 6 + 4] = uibuf[2];
		g_precomp_ui[umm * 6 + 5] = uibuf[3];
	      }
	    }
	  } else {
	    dxx = g_orig_chisq[marker_idx + uii];
	    if (model_adapt) {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, ujj--, uoo, upp, &(g_precomp_ui[umm * 4]), NULL);
	      }
	    } else {
	      for (umm = uii * g_precomp_width; umm < ukk; umm++) {
		chi22_precomp_val_bounds(dxx, ujj, uoo, upp, &(g_precomp_ui[umm * 6]), NULL);
		chi22_precomp_val_bounds(g_maxt_cur_extreme_stat, ujj--, uoo, upp, uibuf, &(g_precomp_d[umm * 2]));
		g_precomp_ui[umm * 6 + 4] = uibuf[2];
		g_precomp_ui[umm * 6 + 5] = uibuf[3];
	      }
	    }
	  }
	}
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
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (spawn_threads(threads, &model_adapt_domrec_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  ulii = 0;
	  model_adapt_domrec_thread((void*)ulii);
	  join_threads(threads, g_assoc_thread_ct);
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads(threads, &model_adapt_trend_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  ulii = 0;
	  model_adapt_trend_thread((void*)ulii);
	  join_threads(threads, g_assoc_thread_ct);
	} else {
	  // TODO
	}
      } else {
	g_maxt_block_base = marker_idx;
	ulii = 0;
	if (model_assoc) {
	  if (g_model_fisher) {
	    if (spawn_threads(threads, &assoc_maxt_fisher_thread, g_assoc_thread_ct)) {
	      logprint(errstr_thread_create);
	      goto model_assoc_ret_THREAD_CREATE_FAIL;
	    }
	    assoc_maxt_fisher_thread((void*)ulii);
	  } else {
	    if (spawn_threads(threads, &assoc_maxt_thread, g_assoc_thread_ct)) {
	      logprint(errstr_thread_create);
	      goto model_assoc_ret_THREAD_CREATE_FAIL;
	    }
	    assoc_maxt_thread((void*)ulii);
	  }
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC)) {
	  if (g_model_fisher) {
	    if (spawn_threads(threads, &model_maxt_domrec_fisher_thread, g_assoc_thread_ct)) {
	      logprint(errstr_thread_create);
	      goto model_assoc_ret_THREAD_CREATE_FAIL;
	    }
	    model_maxt_domrec_fisher_thread((void*)ulii);
	  } else {
	    if (spawn_threads(threads, &model_maxt_domrec_thread, g_assoc_thread_ct)) {
	      logprint(errstr_thread_create);
	      goto model_assoc_ret_THREAD_CREATE_FAIL;
	    }
	    model_maxt_domrec_thread((void*)ulii);
	  }
	} else if (model_modifier & MODEL_PTREND) {
	  if (spawn_threads(threads, &model_maxt_trend_thread, g_assoc_thread_ct)) {
	    logprint(errstr_thread_create);
	    goto model_assoc_ret_THREAD_CREATE_FAIL;
	  }
	  model_maxt_trend_thread((void*)ulii);
	} else {
	  // TODO
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
      retval = multcomp(outname, outname_end, g_marker_uidxs, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, g_orig_1mpval, pfilter, mtest_adjust, adjust_lambda, g_model_fisher);
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
    // printf("extreme stats: %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[perms_total - 1]);
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
