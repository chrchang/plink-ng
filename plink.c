// PLINK 1.90
// Copyright (C) 2005-2014 Shaun Purcell, Christopher Chang

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Uncomment "#define NOLAPACK" in plink_common.h to build without LAPACK.

#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include "plink_common.h"

#ifdef __APPLE__
  #include <sys/types.h>
  #include <sys/sysctl.h>
#endif

#include "plink_assoc.h"
#include "plink_calc.h"
#include "plink_cnv.h"
#include "plink_data.h"
// #include "plink_dosage.h"
#include "plink_filter.h"
#include "plink_glm.h"
#include "plink_help.h"
#include "plink_homozyg.h"
#include "plink_lasso.h"
#include "plink_ld.h"
#include "plink_set.h"
#include "plink_stats.h"
#include "pigz.h"

#define HWE_MIDP 1
#define HWE_THRESH_MIDP 2
#define HWE_THRESH_ALL 4

// default jackknife iterations
#define ITERS_DEFAULT 100000
#define DEFAULT_PPC_GAP 500000
#define MAX_PCS_DEFAULT 20

#define DEFAULT_IBS_TEST_PERMS 100000

#define LOAD_RARE_GRM 1
#define LOAD_RARE_GRM_BIN 2
#define LOAD_RARE_LGEN 4
#define LOAD_RARE_TRANSPOSE 8
#define LOAD_RARE_TPED 0x10
#define LOAD_RARE_TFAM 0x20
#define LOAD_RARE_TRANSPOSE_MASK (LOAD_RARE_TRANSPOSE | LOAD_RARE_TPED | LOAD_RARE_TFAM)
#define LOAD_RARE_DUMMY 0x40
#define LOAD_RARE_SIMULATE 0x80
#define LOAD_RARE_CNV 0x100
#define LOAD_RARE_GVAR 0x200
#define LOAD_RARE_23 0x400
// er, this won't actually be rare...
#define LOAD_RARE_VCF 0x800
#define LOAD_RARE_BCF 0x1000

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^31 - 1
#define PARALLEL_MAX 32768

const char ver_str[] =
#ifdef STABLE_BUILD
  "PLINK v1.90a"
#else
  "PLINK v1.90b1p"
#endif
#ifdef NOLAPACK
  "NL"
#endif
#ifdef __LP64__
  " 64-bit"
#else
  " 32-bit"
#endif
  // include trailing space if day < 10, so character length stays the same
  " (23 Jan 2014)";
const char ver_str2[] =
#ifdef STABLE_BUILD
  "  "
#endif
#ifndef NOLAPACK
  "  "
#endif
  "      https://www.cog-genomics.org/plink2\n"
  "(C) 2005-2014 Shaun Purcell, Christopher Chang   GNU General Public License v3\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
const char errstr_filter_format[] = "Error: Improperly formatted filter file.\n";
const char errstr_freq_format[] = "Error: Improperly formatted frequency file.\n";
const char null_calc_str[] = "Warning: No output requested.  Exiting.\n";
#ifdef STABLE_BUILD
const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --merge-list, --write-snplist, --freqx,\n--missing, --hardy, --ibc, --impute-sex, --indep, --r2, --distance, --genome,\n--homozyg, --make-rel, --make-grm-gz, --rel-cutoff, --cluster, --neighbour,\n--ibs-test, --regress-distance, --model, --gxe, --logistic, --lasso, and\n--fast-epistasis.\n\n'" PROG_NAME_STR " --help | more' describes all functions (warning: long).\n";
#else
  #ifndef NOLAPACK
const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --merge-list, --write-snplist, --freqx,\n--missing, --test-mishap, --hardy, --ibc, --impute-sex, --indep, --r2, --clump,\n--distance, --genome, --homozyg, --make-rel, --make-grm-gz, --rel-cutoff,\n--cluster, --neighbour, --ibs-test, --regress-distance, --model, --gxe,\n--logistic, --lasso, --test-missing, --unrelated-heritability, and\n--fast-epistasis.\n\n'" PROG_NAME_STR " --help | more' describes all functions (warning: long).\n";
  #else
const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --merge-list, --write-snplist, --freqx,\n--missing, --test-mishap, --hardy, --ibc, --impute-sex, --indep, --r2, --clump,\n--distance, --genome, --homozyg, --make-rel, --make-grm-gz, --rel-cutoff,\n--cluster, --neighbour, --ibs-test, --regress-distance, --model, --gxe,\n--logistic, --lasso, --test-missing, and --fast-epistasis.\n\n'" PROG_NAME_STR " --help | more' describes all functions (warning: long).\n";
  #endif
#endif

intptr_t malloc_size_mb = 0;

unsigned char* wkspace;

void dispmsg(int32_t retval) {
  switch (retval) {
  case RET_NOMEM:
    logprint("\nError: Out of memory.  Try the --memory and/or --parallel flags.\n");
    break;
  case RET_WRITE_FAIL:
    logprint("\nError: File write failure.\n");
    break;
  case RET_READ_FAIL:
    logprint("\nError: File read failure.\n");
    break;
  }
}

// back to our regular program

uint32_t random_thin_markers(double thin_keep_prob, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t removed_ct = 0;
  uint32_t uint32_thresh = (uint32_t)(thin_keep_prob * 4294967296.0 + 0.5);
  uint32_t marker_uidx_stop;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      if (sfmt_genrand_uint32(&sfmt) >= uint32_thresh) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
    } while (++marker_uidx < marker_uidx_stop);
  }
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed by --thin.  Try a higher probability.\n");
    return 1;
  }
  sprintf(logbuf, "--thin: %u variant%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", marker_ct - removed_ct);
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

int32_t write_nosex(char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_nm, uintptr_t gender_unk_ct, char* person_ids, uintptr_t max_person_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_uidx = 0;
  int32_t retval = 0;
  uintptr_t* sex_missing;
  uintptr_t indiv_idx;
  if (wkspace_alloc_ul_checked(&sex_missing, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto write_nosex_ret_NOMEM;
  }
  bitfield_exclude_to_include(indiv_exclude, sex_missing, unfiltered_indiv_ct);
  bitfield_andnot(sex_missing, sex_nm, unfiltered_indiv_ctl);
  memcpy(outname_end, ".nosex", 7);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_nosex_ret_OPEN_FAIL;
  }
  for (indiv_idx = 0; indiv_idx < gender_unk_ct; indiv_idx++, indiv_uidx++) {
    next_set_ul_unsafe_ck(sex_missing, &indiv_uidx);
    fputs(&(person_ids[indiv_uidx * max_person_id_len]), outfile);
    putc('\n', outfile);
  }
  if (fclose_null(&outfile)) {
    goto write_nosex_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Ambiguous sex ID%s written to %s.\n", (gender_unk_ct == 1)? "" : "s", outname);
  logprintb();
  while (0) {
  write_nosex_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_nosex_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_nosex_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

double calc_wt_mean(double exponent, int32_t lhi, int32_t lli, int32_t hhi) {
  double lcount = (double)lli + ((double)lhi * 0.5);
  int64_t tot = lhi + lli + hhi;
  double dtot = (double)tot;
  int64_t subcount = lli; // avoid 32-bit integer overflow
  double weight;
  if ((!lhi) && ((!lli) || (!hhi))) {
    return 0.0;
  }
  weight = pow(2 * lcount * (dtot - lcount) / (dtot * dtot), -exponent);
  subcount = lhi * (subcount + hhi) + 2 * subcount * hhi;
  return (subcount * weight * 2) / (double)(tot * tot);
}

int32_t sexcheck(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t do_impute, double check_sex_fthresh, double check_sex_mthresh, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* gender_unk_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t x_variant_ct = 0;
  double nei_sum = 0.0;
  uint32_t gender_unk_ct = 0;
  uint32_t problem_ct = 0;
  int32_t x_code = chrom_info_ptr->x_code;
  int32_t retval = 0;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  // We wish to compute three quantities for each individual:
  // 1. Observed homozygous polymorphic Xchr sites.
  // 2. Observed nonmissing polymorphic Xchr sites.
  // 3. Nei's unbiased estimator of the expected quantity #1.  (This has an
  //    N/(N-1) term which makes it slightly different from GCTA's Fhat2.
  //    Todo: check whether --ibc Fhat2 calculation should be revised to be
  //    consistent with --het...)
  // edit: Actually, forget about the 2N/(2N-1) multiplier for now since it
  // doesn't play well with --freq, and PLINK 1.07 only succeeds in applying it
  // when N=1 due to use of integer division.  Maybe let it be used with
  // inferred MAFs/--freqx/--freq counts later.
  uintptr_t* lptr;
  uint32_t* het_cts;
  uint32_t* missing_cts;
  double* nei_offsets;
  char* fid_ptr;
  char* iid_ptr;
  char* wptr;
  double dpp;
  double dtot;
  double cur_nei;
  double dff;
  double dee;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx_end;
  uintptr_t marker_idxs_left;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t cur_missing_ct;
  uintptr_t allele_obs_ct;
  uintptr_t cur_word;
  uintptr_t ulii;
  uint32_t orig_sex_code;
  uint32_t imputed_sex_code;
  if ((x_code == -1) || (!is_set(chrom_info_ptr->chrom_mask, (uint32_t)x_code))) {
    goto sexcheck_ret_INVALID_CMDLINE;
  }
  marker_uidx_end = chrom_info_ptr->chrom_end[(uint32_t)x_code];
  marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_start[(uint32_t)x_code], marker_uidx_end);
  if (marker_uidx == marker_uidx_end) {
    goto sexcheck_ret_INVALID_CMDLINE;
  }
  marker_idxs_left = marker_uidx_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, marker_uidx_end);
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, indiv_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&het_cts, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&nei_offsets, indiv_ct * sizeof(double))) {
    goto sexcheck_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctl2 - 1] = 0;
  loadbuf[indiv_ctl2 - 1] = 0;
  fill_uint_zero(het_cts, indiv_ct);
  fill_uint_zero(missing_cts, indiv_ct);
  fill_double_zero(nei_offsets, indiv_ct);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
    goto sexcheck_ret_READ_FAIL;
  }
  for (; marker_idxs_left; marker_idxs_left--, marker_uidx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
        goto sexcheck_ret_READ_FAIL;
      }
    }
    if (load_and_collapse(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf, indiv_ct, indiv_exclude, 0)) {
      goto sexcheck_ret_READ_FAIL;
    }
    cur_missing_ct = count_01(loadbuf, indiv_ctl2);
    allele_obs_ct = 2 * (indiv_ct - cur_missing_ct);
    dpp = set_allele_freqs[marker_uidx];
    // skip monomorphic sites
    if ((!allele_obs_ct) || (dpp < 1e-8) || (dpp > (1 - 1e-8))) {
      continue;
    }
    cur_nei = 1.0 - 2 * dpp * (1 - dpp);
    x_variant_ct++;
    if (cur_missing_ct) {
      // iterate through missing calls
      lptr = loadbuf;
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
        cur_word = *lptr++;
        cur_word = (~(cur_word >> 1)) & cur_word & FIVEMASK;
        while (cur_word) {
          ulii = indiv_idx + CTZLU(cur_word) / 2;
          missing_cts[ulii] += 1;
          nei_offsets[ulii] += cur_nei;
          cur_word &= cur_word - 1;
	}
      }
    }
    lptr = loadbuf;
    // iterate through heterozygous calls
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx += BITCT2) {
      cur_word = *lptr++;
      cur_word = (cur_word >> 1) & (~cur_word) & FIVEMASK;
      while (cur_word) {
        het_cts[indiv_idx + CTZLU(cur_word) / 2] += 1;
        cur_word &= cur_word - 1;
      }
    }
    nei_sum += cur_nei;
  }
  if (!x_variant_ct) {
    goto sexcheck_ret_INVALID_CMDLINE;
  }
  memcpy(outname_end, ".sexcheck", 10);
  if (fopen_checked(&outfile, outname, "w")) {
    goto sexcheck_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us       PEDSEX       SNPSEX       STATUS            F\n", plink_maxfid, plink_maxiid);
  fprintf(outfile, tbuf, "FID", "IID");
  indiv_uidx = 0;
  if (do_impute) {
    bitfield_andnot(sex_nm, indiv_exclude, unfiltered_indiv_ctl);
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++, indiv_uidx++) {
    next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    fid_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    iid_ptr = (char*)memchr(fid_ptr, '\t', max_person_id_len);
    wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(iid_ptr - fid_ptr), fid_ptr, tbuf);
    *wptr++ = ' ';
    wptr = fw_strcpy(plink_maxiid, &(iid_ptr[1]), wptr);
    if (!IS_SET(sex_nm, indiv_uidx)) {
      orig_sex_code = 0;
    } else {
      orig_sex_code = 2 - IS_SET(sex_male, indiv_uidx);
    }
    wptr = memseta(wptr, 32, 12);
    *wptr++ = '0' + orig_sex_code;
    wptr = memseta(wptr, 32, 12);
    ulii = x_variant_ct - missing_cts[indiv_idx];
    if (ulii) {
      dee = nei_sum - nei_offsets[indiv_idx];
      dtot = (double)((intptr_t)ulii) - dee;
      dff = (dtot - ((double)((intptr_t)(het_cts[indiv_idx])))) / dtot;
      if (dff > check_sex_mthresh) {
        imputed_sex_code = 1;
      } else if (dff < check_sex_fthresh) {
        imputed_sex_code = 2;
      } else {
        imputed_sex_code = 0;
      }
      *wptr++ = '0' + imputed_sex_code;
      if (orig_sex_code && (orig_sex_code == imputed_sex_code)) {
        wptr = memcpya(wptr, "           OK ", 14);
      } else {
        wptr = memcpya(wptr, "      PROBLEM ", 14);
	problem_ct++;
      }
      wptr = double_g_writewx4x(wptr, dff, 12, '\n');
    } else {
      imputed_sex_code = 0;
      wptr = memcpya(wptr, "0      PROBLEM          nan\n", 28);
      problem_ct++;
    }
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto sexcheck_ret_WRITE_FAIL;
    }
    if (do_impute) {
      if (imputed_sex_code) {
	SET_BIT(sex_nm, indiv_uidx);
	if (imputed_sex_code == 1) {
	  SET_BIT(sex_male, indiv_uidx);
	} else {
	  CLEAR_BIT(sex_male, indiv_uidx);
	}
      } else {
	CLEAR_BIT(sex_nm, indiv_uidx);
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto sexcheck_ret_WRITE_FAIL;
  }
  if (do_impute) {
    bitfield_and(sex_male, sex_nm, unfiltered_indiv_ctl);
    gender_unk_ct = indiv_ct - popcount_longs(sex_nm, unfiltered_indiv_ctl);
    if (!gender_unk_ct) {
      sprintf(logbuf, "--impute-sex: %" PRIuPTR " Xchr variant%s scanned, all sexes imputed.\n", x_variant_ct, (x_variant_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "--impute-sex: %" PRIuPTR " Xchr variant%s scanned, %" PRIuPTR "/%" PRIuPTR " sex%s imputed.\n", x_variant_ct, (x_variant_ct == 1)? "" : "s", (indiv_ct - gender_unk_ct), indiv_ct, (indiv_ct == 1)? "" : "es");
    }
    *gender_unk_ct_ptr = gender_unk_ct;
  } else {
    if (!problem_ct) {
      sprintf(logbuf, "--check-sex: %" PRIuPTR " Xchr variant%s scanned, no problems detected.\n", x_variant_ct, (x_variant_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "--check-sex: %" PRIuPTR " Xchr variant%s scanned, %u problem%s detected.\n", x_variant_ct, (x_variant_ct == 1)? "" : "s", problem_ct, (problem_ct == 1)? "" : "s");
    }
  }
  logprintb();
  sprintf(logbuf, "Report written to %s.\n", outname);
  logprintb();
  while (0) {
  sexcheck_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sexcheck_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sexcheck_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  sexcheck_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  sexcheck_ret_INVALID_CMDLINE:
    logprint("Error: --check-sex/--impute-sex requires at least one polymorphic X chromosome\nsite.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  // possible todo: if any families have been entirely filtered out, don't
  // construct pedigree for them
  unsigned char* wkspace_mark;
  unsigned char* wkspace_mark2;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t initial_family_blocks;
  uintptr_t ulii;
  uintptr_t indiv_uidx;
  uint64_t ullii;
  char* family_ids;
  char* cur_person_id;
  char* last_family_id = NULL;
  char* cur_family_id;
  char* id_ptr;
  uint32_t* family_sizes;
  uint32_t* uiptr;
  uint32_t* uiptr2 = NULL;
  uint32_t fidx;
  int32_t family_size;
  uint32_t* remaining_indiv_idxs;
  int32_t* remaining_indiv_parent_idxs; // -1 = no parent (or nonshared)
  uint32_t remaining_indiv_ct;
  uint32_t indiv_idx_write;
  uintptr_t max_family_id_len = 0;
  char* indiv_ids;
  uint32_t* indiv_id_lookup;
  uintptr_t max_indiv_id_len = 0;
  uintptr_t max_pm_id_len;
  uint32_t family_id_ct;
  uint32_t* fis_ptr;
  char* stray_parent_ids;
  intptr_t stray_parent_ct;
  uintptr_t* processed_indivs;
  uint32_t founder_ct;
  int32_t max_family_nf = 0;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctlm = unfiltered_indiv_ctl * BITCT;
  uint32_t* complete_indiv_idxs;
  uintptr_t complete_indiv_idx_ct;
  double* rs_ptr;
  double* rel_writer;
  double dxx;
  double* tmp_rel_space = NULL;
  double* tmp_rel_writer = NULL;

  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    ujj = strlen_se(&(person_ids[indiv_uidx * max_person_id_len])) + 1;
    if (ujj > max_family_id_len) {
      max_family_id_len = ujj;
    }
    ujj = strlen(&(person_ids[indiv_uidx * max_person_id_len + ujj]));
    if (ujj >= max_indiv_id_len) {
      max_indiv_id_len = ujj + 1;
    }
  }
  if (max_paternal_id_len > max_maternal_id_len) {
    max_pm_id_len = max_paternal_id_len;
  } else {
    max_pm_id_len = max_maternal_id_len;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_space), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_rel_nf_idxs), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_idxs), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&family_ids, unfiltered_indiv_ct * max_family_id_len) ||
      wkspace_alloc_ui_checked(&family_sizes, unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }

  // copy all the items over in order, then qsort, then eliminate duplicates
  // and count family sizes.
  cur_family_id = family_ids;
  cur_person_id = person_ids;
  uiptr = family_sizes;
  *uiptr = 1;
  jj = strlen_se(cur_person_id);
  memcpyx(cur_family_id, cur_person_id, jj, 0);
  for (indiv_uidx = 1; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    cur_person_id = &(cur_person_id[max_person_id_len]);
    mm = strlen_se(cur_person_id);
    if ((jj != mm) || memcmp(cur_family_id, cur_person_id, mm)) {
      cur_family_id = &(cur_family_id[max_family_id_len]);
      memcpyx(cur_family_id, cur_person_id, mm, 0);
      jj = mm;
      *(++uiptr) = 1;
    } else {
      *uiptr += 1;
    }
  }
  initial_family_blocks = 1 + (uint32_t)(uiptr - family_sizes);
  if (qsort_ext(family_ids, initial_family_blocks, max_family_id_len, strcmp_deref, (char*)family_sizes, sizeof(int32_t))) {
    return RET_NOMEM;
  }

  last_family_id = family_ids;
  cur_family_id = &(family_ids[max_family_id_len]);
  family_id_ct = 1;
  uii = 1; // read idx
  if (initial_family_blocks != 1) {
    uiptr = family_sizes;
    while (strcmp(cur_family_id, last_family_id)) {
      family_id_ct++;
      uiptr++;
      last_family_id = cur_family_id;
      cur_family_id = &(cur_family_id[max_family_id_len]);
      uii++;
      if (uii == initial_family_blocks) {
	break;
      }
    }
    if (uii < initial_family_blocks) {
      uiptr2 = uiptr; // family_sizes read pointer
      *uiptr += *(++uiptr2);
      uii++;
      cur_family_id = &(cur_family_id[max_family_id_len]); // read pointer
      while (uii < initial_family_blocks) {
	while (!strcmp(cur_family_id, last_family_id)) {
	  *uiptr += *(++uiptr2);
	  uii++;
	  if (uii == initial_family_blocks) {
	    break;
	  }
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
	if (uii < initial_family_blocks) {
	  *(++uiptr) = *(++uiptr2);
	  last_family_id = &(last_family_id[max_family_id_len]);
	  strcpy(last_family_id, cur_family_id);
	  family_id_ct++;
	  uii++;
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
      }
    }
  }

  uiptr = family_sizes;
  if (family_id_ct < unfiltered_indiv_ct) {
    wkspace_reset(family_ids);
    family_ids = (char*)wkspace_alloc(family_id_ct * max_family_id_len);
    family_sizes = (uint32_t*)wkspace_alloc(family_id_ct * sizeof(int32_t));
    for (uii = 0; uii < family_id_ct; uii++) {
      family_sizes[uii] = *uiptr++;
    }
  }
  pri_ptr->family_ids = family_ids;
  pri_ptr->family_id_ct = family_id_ct;
  pri_ptr->max_family_id_len = max_family_id_len;
  pri_ptr->family_sizes = family_sizes;

  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_offsets), (family_id_ct + 1) * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&(pri_ptr->family_rel_space_offsets), (family_id_ct + 1) * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_founder_cts), family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_int_zero((int32_t*)(pri_ptr->family_founder_cts), family_id_ct);

  ii = 0; // running family_info offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_info_offsets[fidx] = ii;
    ii += family_size;
  }

  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_uint_zero(uiptr, family_id_ct);

  // Fill family_idxs, family_founder_cts, and founder portion of
  // family_rel_nf_idxs.
  cur_person_id = person_ids;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    kk = bsearch_str(cur_person_id, strlen_se(cur_person_id), family_ids, max_family_id_len, family_id_ct);
    pri_ptr->family_idxs[indiv_uidx] = kk;
    if (IS_SET(founder_info, indiv_uidx)) {
      pri_ptr->family_founder_cts[(uint32_t)kk] += 1;
      pri_ptr->family_rel_nf_idxs[indiv_uidx] = uiptr[(uint32_t)kk];
      uiptr[kk] += 1;
    }
    cur_person_id = &(cur_person_id[max_person_id_len]);
  }
  wkspace_reset(uiptr);
  ulii = 0; // running rel_space offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_rel_space_offsets[fidx] = ulii;
    kk = pri_ptr->family_founder_cts[fidx];
    if (family_size - kk > max_family_nf) {
      max_family_nf = family_size - kk;
    }
    // No need to explicitly store the (kk * (kk - 1)) / 2 founder-founder
    // relationships.
    ulii += (((int64_t)family_size) * (family_size - 1) - ((int64_t)kk) * (kk - 1)) / 2;
  }

  // make it safe to determine size of blocks by subtracting from the next
  // offset, even if we're at the last family
  pri_ptr->family_info_offsets[family_id_ct] = unfiltered_indiv_ct;
  pri_ptr->family_rel_space_offsets[family_id_ct] = ulii;
  if (wkspace_alloc_d_checked(&(pri_ptr->rel_space), ulii * sizeof(double))) {
    return RET_NOMEM;
  }

  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  // populate family_info_space
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    uiptr[fidx] = pri_ptr->family_info_offsets[fidx];
  }
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    fidx = pri_ptr->family_idxs[indiv_uidx];
    pri_ptr->family_info_space[uiptr[fidx]] = indiv_uidx;
    uiptr[fidx] += 1;
  }
  wkspace_reset(wkspace_mark);

  if (wkspace_alloc_ul_checked(&processed_indivs, (unfiltered_indiv_ctl + (max_family_nf + (BITCT2 - 1)) / BITCT2) * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_one(&(processed_indivs[unfiltered_indiv_ctl]), (max_family_nf + (BITCT2 - 1)) / BITCT2);

  wkspace_mark2 = wkspace_base;
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = family_sizes[fidx];
    founder_ct = pri_ptr->family_founder_cts[fidx];
    remaining_indiv_ct = family_size - founder_ct;
    stray_parent_ct = 0;
    if (remaining_indiv_ct) {
      memcpy(processed_indivs, founder_info, unfiltered_indiv_ctl * sizeof(intptr_t));
      if (wkspace_alloc_ui_checked(&complete_indiv_idxs, family_size * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&remaining_indiv_idxs, remaining_indiv_ct * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&indiv_ids, family_size * max_indiv_id_len) ||
          wkspace_alloc_ui_checked(&indiv_id_lookup, family_size * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&remaining_indiv_parent_idxs, remaining_indiv_ct * 2 * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&stray_parent_ids, remaining_indiv_ct * 2 * max_pm_id_len)) {
	return RET_NOMEM;
      }
      ii = pri_ptr->family_info_offsets[fidx];
      fis_ptr = &(pri_ptr->family_info_space[ii]);
      rs_ptr = &(pri_ptr->rel_space[pri_ptr->family_rel_space_offsets[fidx]]);
      rel_writer = rs_ptr;
      cur_person_id = indiv_ids;
      for (ii = 0; ii < family_size; ii++) {
	kk = fis_ptr[(uint32_t)ii];
	jj = strlen_se(&(person_ids[kk * max_person_id_len])) + 1;
	strcpy(cur_person_id, &(person_ids[kk * max_person_id_len + jj]));
	cur_person_id = &(cur_person_id[max_indiv_id_len]);
	indiv_id_lookup[(uint32_t)ii] = ii;
      }

      if (qsort_ext(indiv_ids, family_size, max_indiv_id_len, strcmp_deref, (char*)indiv_id_lookup, sizeof(int32_t))) {
	return RET_NOMEM;
      }
      // Compile list of non-founder family member indices, and identify
      // parents who are referred to by individual ID but are NOT in the
      // dataset.
      ii = 0; // family_info_space index
      complete_indiv_idx_ct = 0;
      cur_person_id = stray_parent_ids;
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	while (IS_SET(founder_info, fis_ptr[ii])) {
	  complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[ii];
	  ii++;
	}
	kk = fis_ptr[ii++];

	// does not track sex for now
	id_ptr = &(paternal_ids[((uint32_t)kk) * max_paternal_id_len]);
	if (memcmp("0", id_ptr, 2)) {
	  ujj = strlen(id_ptr);
	  mm = bsearch_str(id_ptr, ujj, indiv_ids, max_indiv_id_len, family_size);
	  if (mm == -1) {
	    memcpy(cur_person_id, id_ptr, ujj + 1);
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2] = -2;
	  } else {
            remaining_indiv_parent_idxs[uii * 2] = fis_ptr[indiv_id_lookup[(uint32_t)mm]];
	  }
	} else {
          remaining_indiv_parent_idxs[uii * 2] = -1;
	}
	id_ptr = &(maternal_ids[((uint32_t)kk) * max_maternal_id_len]);
	if (memcmp("0", id_ptr, 2)) {
	  ujj = strlen(id_ptr);
          mm = bsearch_str(id_ptr, ujj, indiv_ids, max_indiv_id_len, family_size);
	  if (mm == -1) {
	    memcpy(cur_person_id, id_ptr, ujj + 1);
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2 + 1] = -2;
	  } else {
	    remaining_indiv_parent_idxs[uii * 2 + 1] = fis_ptr[indiv_id_lookup[(uint32_t)mm]];
	  }
	} else {
	  remaining_indiv_parent_idxs[uii * 2 + 1] = -1;
	}
        remaining_indiv_idxs[uii] = kk;
      }
      while (ii < family_size) {
	complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[(uint32_t)ii];
	ii++;
      }
      qsort(stray_parent_ids, stray_parent_ct, max_pm_id_len, strcmp_casted);
      cur_person_id = stray_parent_ids;
      ii = 0; // read idx
      jj = 0; // write idx

      // Now filter out all such parents who aren't referenced at least twice.
      while (ii + 1 < stray_parent_ct) {
        if (strcmp(&(stray_parent_ids[ii * max_pm_id_len]), &(stray_parent_ids[(ii + 1) * max_pm_id_len]))) {
	  ii++;
	  continue;
	}
	ii++;
	strcpy(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len]));
	do {
	  ii++;
        } while (!(strcmp(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len])) || (ii > stray_parent_ct)));
        cur_person_id = &(cur_person_id[max_pm_id_len]);
	jj++;
      }
      stray_parent_ct = jj;

      // Now allocate temporary relatedness table between nonfounders and
      // stray parents with multiple references.
      if (stray_parent_ct) {
        if (wkspace_alloc_d_checked(&tmp_rel_space, (family_size - founder_ct) * stray_parent_ct * sizeof(double))) {
	  return RET_NOMEM;
        }
	tmp_rel_writer = tmp_rel_space;
      }

      // Now fill in remainder of remaining_indiv_parent_idxs.
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	jj = remaining_indiv_idxs[uii];
	if (remaining_indiv_parent_idxs[uii * 2] == -2) {
	  kk = bsearch_str_nl(&(paternal_ids[((uint32_t)jj) * max_paternal_id_len]), stray_parent_ids, max_pm_id_len, stray_parent_ct);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2] = kk;
	}
	if (remaining_indiv_parent_idxs[uii * 2 + 1] == -2) {
	  kk = bsearch_str_nl(&(maternal_ids[((uint32_t)jj) * max_maternal_id_len]), stray_parent_ids, max_pm_id_len, stray_parent_ct);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2 + 1] = kk;
	}
      }
      ullii = ((uint64_t)founder_ct) * (founder_ct - 1);
      while (remaining_indiv_ct) {
	indiv_idx_write = 0;
	for (uii = 0; uii < remaining_indiv_ct; uii++) {
	  kk = remaining_indiv_parent_idxs[uii * 2];
	  mm = remaining_indiv_parent_idxs[uii * 2 + 1];
	  jj = remaining_indiv_idxs[uii];
	  if (((kk == -1) || is_set(processed_indivs, kk)) && ((mm == -1) || is_set(processed_indivs, mm))) {
	    for (ujj = 0; ujj < founder_ct; ujj++) {
	      // relationship between kk and ujjth founder
	      if ((kk >= (int32_t)unfiltered_indiv_ct) || (kk == -1)) {
		dxx = 0.0;
	      } else if (is_set(founder_info, kk)) {
		if (kk == (int32_t)complete_indiv_idxs[ujj]) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[(uint32_t)kk];
                dxx = 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
	      }
	      if (is_set(founder_info, mm)) {
		if (mm == (int32_t)complete_indiv_idxs[ujj]) {
		  dxx += 0.5;
		}
	      } else if ((mm != -1) && (mm < (int32_t)unfiltered_indiv_ct)) {
		ukk = pri_ptr->family_rel_nf_idxs[(uint32_t)mm];
		dxx += 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
	      }
	      *rel_writer++ = dxx;
	    }
	    for (; ujj < complete_indiv_idx_ct; ujj++) {
	      if (kk == -1) {
		dxx = 0.0;
	      } else if (kk >= (int32_t)unfiltered_indiv_ct) {
		dxx = 0.5 * tmp_rel_space[(ujj - founder_ct) * stray_parent_ct + kk - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, kk)) {
                dxx = 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + pri_ptr->family_rel_nf_idxs[kk]];
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[kk];
		if (ukk == ujj) {
		  dxx = 0.5;
		} else if (ukk < ujj) {
		  dxx = 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + ukk];
		} else {
		  dxx = 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
		}
	      }
	      if (mm >= (int32_t)unfiltered_indiv_ct) {
		dxx += 0.5 * tmp_rel_space[(ujj - founder_ct) * stray_parent_ct + mm - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, mm)) {
		dxx += 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + pri_ptr->family_rel_nf_idxs[mm]];
	      } else if (mm != -1) {
		ukk = pri_ptr->family_rel_nf_idxs[mm];
		if (ukk == ujj) {
		  dxx += 0.5;
		} else if (ukk < ujj) {
		  dxx += 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + ukk];
		} else {
		  dxx += 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
		}
	      }
	      *rel_writer++ = dxx;
	    }
	    for (ujj = 0; ujj < (uintptr_t)stray_parent_ct; ujj++) {
	      if (kk >= (int32_t)unfiltered_indiv_ct) {
		if ((uint32_t)kk == ujj + unfiltered_indiv_ctlm) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else if (kk == -1) {
                dxx = 0.0;
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[kk];
		if (ukk < founder_ct) {
		  dxx = 0.0;
		} else {
                  dxx = 0.5 * tmp_rel_space[(ukk - founder_ct) * stray_parent_ct + ujj];
		}
	      }
	      if (mm >= (int32_t)unfiltered_indiv_ct) {
		if ((uint32_t)mm == ujj + unfiltered_indiv_ctlm) {
		  dxx += 0.5;
		}
	      } else if (mm != -1) {
		ukk = pri_ptr->family_rel_nf_idxs[mm];
		if (ukk >= founder_ct) {
		  dxx += 0.5 * tmp_rel_space[(ukk - founder_ct) * stray_parent_ct + ujj];
		}
	      }
	      *tmp_rel_writer++ = dxx;
	    }
	    pri_ptr->family_rel_nf_idxs[jj] = complete_indiv_idx_ct;
	    complete_indiv_idxs[complete_indiv_idx_ct++] = jj;
	    set_bit(processed_indivs, jj);
	  } else {
            remaining_indiv_parent_idxs[indiv_idx_write * 2] = kk;
	    remaining_indiv_parent_idxs[indiv_idx_write * 2 + 1] = mm;
	    remaining_indiv_idxs[indiv_idx_write++] = jj;
	  }
	}
	if (indiv_idx_write == remaining_indiv_ct) {
	  logprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
	  return RET_INVALID_FORMAT;
	}
	remaining_indiv_ct = indiv_idx_write;
      }
      wkspace_reset(wkspace_mark2);
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t make_founders(uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t require_two, uintptr_t* indiv_exclude, uintptr_t* founder_info) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uint32_t new_founder_ct = 0;
  int32_t retval = 0;
  char* sorted_ids;
  char* id_buf;
  char* wptr;
  char* pat_ptr;
  char* mat_ptr;
  uintptr_t* nf_bitarr;
  uintptr_t indiv_uidx;
  uint32_t fam_len_p1;
  uint32_t missing_parent_ct;
  uint32_t cur_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&nf_bitarr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto make_founders_ret_NOMEM;
  }
  bitfield_exclude_to_include(indiv_exclude, nf_bitarr, unfiltered_indiv_ct);
  bitfield_andnot(nf_bitarr, founder_info, unfiltered_indiv_ctl);
  indiv_uidx = next_set(nf_bitarr, 0, unfiltered_indiv_ct);
  if (indiv_uidx == unfiltered_indiv_ct) {
    logprint("Note: Skipping --make-founders since there are no nonfounders.\n");
    goto make_founders_ret_1;
  }
  sorted_ids = alloc_and_init_collapsed_arr(person_ids, max_person_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_ct, 0);
  if (!sorted_ids) {
    goto make_founders_ret_NOMEM;
  }
  qsort(sorted_ids, indiv_ct, max_person_id_len, strcmp_casted);
  do {
    pat_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    fam_len_p1 = strlen_se(pat_ptr) + 1;
    wptr = memcpya(id_buf, pat_ptr, fam_len_p1);
    missing_parent_ct = 0;
    pat_ptr = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
    cur_len = strlen(pat_ptr);
    if (cur_len + fam_len_p1 >= max_person_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, pat_ptr, cur_len);
      if (bsearch_str(id_buf, cur_len + fam_len_p1, sorted_ids, max_person_id_len, indiv_ct) == -1) {
	missing_parent_ct++;
      }
    }
    mat_ptr = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
    cur_len = strlen(mat_ptr);
    if (cur_len + fam_len_p1 >= max_person_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, mat_ptr, cur_len);
      if (bsearch_str(id_buf, cur_len + fam_len_p1, sorted_ids, max_person_id_len, indiv_ct) == -1) {
	missing_parent_ct++;
      }
    }
    if (missing_parent_ct > require_two) {
      SET_BIT(founder_info, indiv_uidx);
      memcpy(pat_ptr, "0", 2);
      memcpy(mat_ptr, "0", 2);
      new_founder_ct++;
    }
    indiv_uidx++;
    next_set_ul_ck(nf_bitarr, &indiv_uidx, unfiltered_indiv_ct);
  } while (indiv_uidx < unfiltered_indiv_ct);
  sprintf(logbuf, "--make-founders: %u individual%s affected.\n", new_founder_ct, (new_founder_ct == 1)? "" : "s");
  logprintb();
  while (0) {
  make_founders_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
 make_founders_ret_1:
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_indiv_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr) {
  uint32_t mp_strlen = strlen(makepheno_str);
  uint32_t makepheno_all = ((mp_strlen == 1) && (makepheno_str[0] == '*'));
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* pheno_c = *pheno_c_ptr;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int32_t ii;
  uint32_t person_idx;
  uint32_t tmp_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  if (!pheno_c) {
    pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
    *pheno_c_ptr = pheno_c;
  }
  if (makepheno_all) {
    fill_all_bits(pheno_nm, unfiltered_indiv_ct);
  }
  // probably want to permit long lines here
  tbuf[MAXLINELEN - 1] = ' '; 
  while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in phenotype file.\n");
      return RET_INVALID_FORMAT;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    if (bsearch_read_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, unfiltered_indiv_ct, bufptr0, &bufptr, &ii)) {
      logprint(errstr_phenotype_format);
      return RET_INVALID_FORMAT;
    }
    if (ii != -1) {
      person_idx = id_map[(uint32_t)ii];
      if (makepheno_all) {
	SET_BIT(pheno_c, person_idx);
      } else {
	SET_BIT(pheno_nm, person_idx);
        tmp_len = strlen_se(bufptr);
	if ((tmp_len == mp_strlen) && (!memcmp(bufptr, makepheno_str, mp_strlen))) {
	  SET_BIT(pheno_c, person_idx);
	}
      }
    }
  }
  if (!feof(phenofile)) {
    return RET_READ_FAIL;
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t convert_tail_pheno(uint32_t unfiltered_indiv_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod) {
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  uint32_t indiv_uidx;
  uint32_t indiv_uidx_stop;
  double dxx;
  if (!(*pheno_d_ptr)) {
    logprint("Error: --tail-pheno requires scalar phenotype data.\n");
    return RET_INVALID_FORMAT;
  }
  indiv_uidx = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (!pheno_c) {
    pheno_c = (uintptr_t*)malloc(indiv_uidx * sizeof(intptr_t));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    *pheno_c_ptr = pheno_c;
  }
  fill_ulong_zero(pheno_c, indiv_uidx);
  indiv_uidx = 0;
  do {
    indiv_uidx = next_set(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    indiv_uidx_stop = next_unset(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    for (; indiv_uidx < indiv_uidx_stop; indiv_uidx++) {
      dxx = pheno_d[indiv_uidx];
      if (dxx > tail_bottom) {
        if (dxx > tail_top) {
          SET_BIT(pheno_c, indiv_uidx);
        } else {
	  CLEAR_BIT(pheno_nm, indiv_uidx);
        }
      }
    }
  } while (indiv_uidx_stop < unfiltered_indiv_ct);
  free(pheno_d);
  *pheno_d_ptr = NULL;
  return 0;
}

const unsigned char acgt_reverse_arr[] = "1B2DEF3HIJKLMNOPQRS4";
const unsigned char acgt_arr[] = "ACGT";
// g_one_char_strs offsets (double)
const unsigned char acgt_reverse_arr1[] = "b\204d\210\212\214f\220\222\224\226\230\232\234\236\240\242\244\246h";
const unsigned char acgt_arr1[] = "\202\206\216\250";
// const unsigned char acgt_reverse_arr1[] = "\"D$HJL&PRTVXZ\\^`bdf(";
// const unsigned char acgt_arr1[] = "BFNh";

static inline unsigned char conditional_convert(unsigned char diff, unsigned char ucc2_max, const unsigned char* convert_arr, unsigned char ucc) {
  unsigned char ucc2 = ucc - diff;
  return (ucc2 < ucc2_max)? convert_arr[ucc2] : ucc;
}

static inline void conditional_convert1(unsigned char diff, unsigned char ucc2_max, const unsigned char* convert_arr, char** allele_ptr) {
  unsigned char ucc2 = ((unsigned char)(**allele_ptr)) - diff;
  if (ucc2 < ucc2_max) {
    *allele_ptr = (char*)(&(g_one_char_strs[convert_arr[ucc2]]));
  }
}

void allelexxxx_recode(uint32_t allelexxxx, char** marker_allele_ptrs, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t recode_multichar = allelexxxx & ALLELE_RECODE_MULTICHAR;
  const unsigned char* convert_arr;
  const unsigned char* convert_arr1;
  char** map_ptr;
  char** map_ptr_stop;
  char* cptr;
  uint32_t marker_uidx_stop;
  unsigned char diff;
  unsigned char ucc2_max;
  unsigned char ucc;
  if (allelexxxx & ALLELE_RECODE_ACGT) {
    diff = 49;
    ucc2_max = 4;
    convert_arr = acgt_arr;
    convert_arr1 = acgt_arr1;
  } else {
    diff = 65;
    ucc2_max = 20;
    convert_arr = acgt_reverse_arr;
    convert_arr1 = acgt_reverse_arr1;
  }
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    map_ptr = &(marker_allele_ptrs[marker_uidx * 2]);
    map_ptr_stop = &(marker_allele_ptrs[marker_uidx_stop * 2]);
    marker_uidx = marker_uidx_stop;
    if (recode_multichar) {
      do {
	cptr = *map_ptr;
	if (!cptr[1]) {
	  conditional_convert1(diff, ucc2_max, convert_arr1, map_ptr);
	} else {
	  ucc = *cptr;
	  do {
	    *cptr = conditional_convert(diff, ucc2_max, convert_arr, ucc);
	    ucc = *(++cptr);
	  } while (ucc);
	}
      } while (++map_ptr < map_ptr_stop);
    } else {
      do {
        if (!(map_ptr[0][1])) {
          conditional_convert1(diff, ucc2_max, convert_arr1, map_ptr);
	}
      } while (++map_ptr < map_ptr_stop);
    }
  }
}

void calc_plink_maxfid(uint32_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uint32_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t* plink_maxfid_ptr, uint32_t* plink_maxiid_ptr) {
  uintptr_t plink_maxfid = 4;
  uintptr_t plink_maxiid = 4;
  uint32_t indiv_uidx = 0;
  uint32_t indivs_done = 0;
  char* cptr;
  char* cptr2;
  char* cptr_end;
  uintptr_t slen;
  uint32_t indiv_uidx_stop;
  // imitate PLINK 1.07 behavior (see Plink::prettyPrintLengths() in
  // helper.cpp), to simplify testing and avoid randomly breaking existing
  // scripts
  do {
    indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
    indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
    indivs_done += indiv_uidx_stop - indiv_uidx;
    cptr = &(person_ids[indiv_uidx * max_person_id_len]);
    cptr_end = &(person_ids[indiv_uidx_stop * max_person_id_len]);
    indiv_uidx = indiv_uidx_stop;
    do {
      cptr2 = (char*)memchr(cptr, '\t', max_person_id_len);
      slen = (uintptr_t)(cptr2 - cptr);
      if (slen > plink_maxfid) {
	plink_maxfid = slen + 2;
      }
      slen = strlen(&(cptr2[1]));
      if (slen > plink_maxiid) {
        plink_maxiid = slen + 2;
      }
      cptr = &(cptr[max_person_id_len]);
    } while (cptr < cptr_end);
  } while (indivs_done < indiv_ct);
  *plink_maxfid_ptr = plink_maxfid;
  *plink_maxiid_ptr = plink_maxiid;
}

uint32_t calc_plink_maxsnp(uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len) {
  uintptr_t plink_maxsnp = 4;
  uintptr_t max_marker_id_len_m1 = max_marker_id_len - 1;
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  char* cptr;
  char* cptr_end;
  uintptr_t slen;
  uint32_t marker_uidx_stop;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
    cptr_end = &(marker_ids[marker_uidx_stop * max_marker_id_len]);
    marker_uidx = marker_uidx_stop;
    do {
      slen = strlen(cptr);
      if (slen > plink_maxsnp) {
	plink_maxsnp = slen + 2;
	if (plink_maxsnp >= max_marker_id_len_m1) {
	  return plink_maxsnp;
	}
      }
      cptr = &(cptr[max_marker_id_len]);
    } while (cptr < cptr_end);
  }
  return plink_maxsnp;
}

// aptr1 = minor, aptr2 = major
int32_t load_one_freq(uint32_t alen1, const char* aptr1, uint32_t alen2, const char* aptr2, double maf, double* set_allele_freq_ptr, char* mastr1, char* mastr2, char missing_geno) {
  uint32_t malen1 = strlen(mastr1);
  uint32_t malen2 = strlen(mastr2);
  uint32_t uii;
  const char* aptr;
  if (maf > 0.5) {
    aptr = aptr2;
    uii = alen2;
    aptr2 = aptr1;
    alen2 = alen1;
    aptr1 = aptr;
    alen1 = uii;
    maf = 1.0 - maf;
  }
  if ((malen1 == alen1) && (!memcmp(mastr1, aptr1, alen1))) {
    if ((malen2 == alen2) && (!memcmp(mastr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0 - maf;
    } else {
      return -1;
    }
  } else if ((malen2 == alen1) && (!memcmp(mastr2, aptr1, alen1))) {
    if ((malen1 == alen2) && (!memcmp(mastr1, aptr2, alen2))) {
      *set_allele_freq_ptr = maf;
    } else {
      return -1;
    }
  } else if ((*aptr1 == missing_geno) && (alen1 == 1) && (maf == 0.0)) {
    if ((malen1 == alen2) && (!memcmp(mastr1, aptr2, alen2))) {
      *set_allele_freq_ptr = 0.0;
    } else if ((malen2 == alen2) && (!memcmp(mastr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

int32_t read_external_freqs(char* freqname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs, uint32_t* marker_allele_cts, double* set_allele_freqs, uint32_t maf_succ, double exponent, uint32_t wt_needed, double* marker_weights) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* freqfile = NULL;
  uint32_t freq_counts = 0;
  uint32_t alen1 = 0;
  uint32_t alen2 = 0;
  char* aptr1 = NULL;
  char* aptr2 = NULL;
  int32_t retval = 0;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;
  char* loadbuf;
  char* sorted_ids;
  uint32_t* id_map;
  uintptr_t loadbuf_size;
  uint32_t chrom_idx;
  uint32_t marker_uidx;
  uint32_t uii;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  double maf;
  int32_t c_hom_a1;
  int32_t c_het;
  int32_t c_hom_a2;
  int32_t c_hap_a1;
  int32_t c_hap_a2;
  int32_t ii;
  if (fopen_checked(&freqfile, freqname, "r")) {
    goto read_external_freqs_ret_OPEN_FAIL;
  }
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto read_external_freqs_ret_1;
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto read_external_freqs_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  if (fgets(loadbuf, loadbuf_size, freqfile) == NULL) {
    logprint("Error: Empty --read-freq file.\n");
    goto read_external_freqs_ret_INVALID_FORMAT_2;
  }
  if (!memcmp(loadbuf, " CHR  ", 6)) {
    uii = strlen(loadbuf);
    if (loadbuf[uii - 2] == '0') { // --counts makes G0 the last column header
      freq_counts = 1;
    } else if (loadbuf[uii - 2] != 'S') { // NCHROBS
      goto read_external_freqs_ret_INVALID_FORMAT;
    }
    while (fgets(loadbuf, loadbuf_size, freqfile) != NULL) {
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
      ii = get_chrom_code(chrom_info_ptr, bufptr);
      if (ii == -1) {
	goto read_external_freqs_ret_INVALID_FORMAT;
      }
      chrom_idx = ii;
      bufptr = next_item(bufptr); // now at beginning of marker name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      ii = bsearch_str(bufptr, strlen_se(bufptr), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
        if ((chrom_idx == get_marker_chrom(chrom_info_ptr, marker_uidx)) || (!chrom_idx) || (!get_marker_chrom(chrom_info_ptr, marker_uidx))) {
	  alen1 = strlen_se(bufptr2);
	  aptr1 = bufptr2;
	  bufptr2 = next_item(bufptr2);
	  if (no_more_items_kns(bufptr2)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  alen2 = strlen_se(bufptr2);
	  aptr2 = bufptr2;
	  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr2);
	  if (no_more_items_kns(bufptr)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  if (freq_counts) {
	    if (no_more_items_kns(next_item(bufptr))) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	    c_hom_a1 = atoi(bufptr);
	    c_hom_a2 = atoi(next_item(bufptr));
	    maf = ((double)c_hom_a1 + maf_succ) / ((double)(c_hom_a1 + c_hom_a2 + 2 * maf_succ));
	  } else {
	    if (scan_double(bufptr, &maf)) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  }
	  if (load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	  }
        }
      }
    }
    if (freq_counts) {
      logprint(".frq.count file loaded.\n");
    } else {
      logprint(".frq file loaded.\n");
    }
  } else if (!memcmp(loadbuf, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)", 71)) {
    // changed from strcmp to avoid eoln problems
    // known --freqx format, v0.15.3 or later
    while (fgets(loadbuf, loadbuf_size, freqfile) != NULL) {
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      ii = get_chrom_code(chrom_info_ptr, loadbuf);
      if (ii == -1) {
	goto read_external_freqs_ret_INVALID_FORMAT;
      }
      chrom_idx = ii;
      bufptr = next_item(loadbuf); // now at beginning of marker name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      ii = bsearch_str(bufptr, strlen_se(bufptr), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
        if ((chrom_idx == get_marker_chrom(chrom_info_ptr, marker_uidx)) || (!chrom_idx) || (!get_marker_chrom(chrom_info_ptr, marker_uidx))) {
	  alen1 = strlen_se(bufptr2);
	  aptr1 = bufptr2;
	  bufptr2 = next_item(bufptr2);
	  if (no_more_items_kns(bufptr2)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  alen2 = strlen_se(bufptr2);
	  aptr2 = bufptr2;
	  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr2);
	  bufptr2 = next_item(bufptr);
	  bufptr3 = next_item(bufptr2);
	  bufptr4 = next_item(bufptr3);
	  bufptr5 = next_item(bufptr4);
	  if (no_more_items_kns(bufptr5)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  c_hom_a1 = atoi(bufptr);
	  c_het = atoi(bufptr2);
	  c_hom_a2 = atoi(bufptr3);
	  c_hap_a1 = atoi(bufptr4);
	  c_hap_a2 = atoi(bufptr5);
	  maf = ((double)(c_hom_a1 * 2 + c_het + c_hap_a1 + maf_succ)) / ((double)(2 * (c_hom_a1 + c_het + c_hom_a2 + maf_succ) + c_hap_a1 + c_hap_a2));
	  if (load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    if (c_hap_a1 || c_hap_a2) {
	      marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	    } else {
	      marker_weights[marker_uidx] = calc_wt_mean(exponent, c_het, c_hom_a1, c_hom_a2);
	    }
	  }
        }
      }
    }
    logprint(".frqx file loaded.\n");
  } else {
    // Also support GCTA-style frequency files:
    // [marker ID]\t[reference allele]\t[frequency of reference allele]\n
    do {
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr = next_item(bufptr);
      if (!bufptr) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      ii = bsearch_str(loadbuf, strlen_se(loadbuf), sorted_ids, max_marker_id_len, unfiltered_marker_ct - marker_exclude_ct);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
	alen1 = strlen_se(bufptr);
	aptr1 = bufptr;
        bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr)) {
          goto read_external_freqs_ret_INVALID_FORMAT;
	}
	if (scan_double(bufptr, &maf)) {
          goto read_external_freqs_ret_INVALID_FORMAT;
        }
	if (load_one_freq(1, missing_geno_ptr, alen1, aptr1, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	  goto read_external_freqs_ret_ALLELE_MISMATCH;
	}
	if (wt_needed) {
	  marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	}
      } else {
	// if there aren't exactly 3 columns, this isn't a GCTA .freq file
	bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr) || (!no_more_items_kns(next_item(bufptr)))) {
	  goto read_external_freqs_ret_INVALID_FORMAT;
	}
      }
    } while (fgets(loadbuf, loadbuf_size, freqfile) != NULL);
    logprint("GCTA-formatted .freq file loaded.\n");
  }
  while (0) {
  read_external_freqs_ret_TOO_LONG_LINE:
    if (loadbuf_size == MAXLINEBUFLEN) {
      logprint("Error: Pathologically long line in --freq{x} file.\n");
      retval = RET_INVALID_FORMAT;
      break;
    }
  read_external_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  read_external_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_external_freqs_ret_INVALID_FORMAT:
    logprint(errstr_freq_format);
  read_external_freqs_ret_INVALID_FORMAT_2:
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_ALLELE_MISMATCH:
    sprintf(logbuf, "Error: Mismatch between .bim/.ped and --read-freq alleles at %s.\n", next_item(skip_initial_spaces(loadbuf)));
    logprintb();
    retval = RET_ALLELE_MISMATCH;
    break;
  }
 read_external_freqs_ret_1:
  fclose_cond(freqfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t write_freqs(char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t indiv_f_ct, uint32_t indiv_f_male_ct, uint64_t misc_flags, uintptr_t* marker_reverse) {
  FILE* outfile = NULL;
  uint32_t reverse = 0;
  uint32_t freq_counts = (misc_flags / MISC_FREQ_COUNTS) & 1;
  uint32_t freqx = (misc_flags / MISC_FREQX) & 1;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  char* tbuf2 = &(tbuf[MAXLINELEN]);
  int32_t retval = 0;
  char* minor_ptr;
  char* major_ptr;
  char* bufptr;
  uint32_t chrom_end;
  uint32_t marker_uidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t missing_ct;
  int32_t chrom_idx;
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_freqs_ret_OPEN_FAIL;
  }
  if (freqx) {
    if (fputs_checked("CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)\n", outfile)) {
      goto write_freqs_ret_WRITE_FAIL;
    }
  } else if (plink_maxsnp < 5) {
    if (freq_counts) {
      if (fputs_checked(" CHR  SNP   A1   A2     C1     C2     G0\n", outfile)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
      strcpy(tbuf, " %4s %4s %4s %6u %6u %6u\n");
    } else {
      if (fputs_checked(" CHR  SNP   A1   A2          MAF  NCHROBS\n", outfile)) {
        goto write_freqs_ret_WRITE_FAIL;
      }
      strcpy(tbuf, " %4s %4s %4s %12.4g %8d\n");
    }
  } else if (freq_counts) {
    sprintf(tbuf, " CHR %%%us   A1   A2     C1     C2     G0\n", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
    sprintf(tbuf, " %%%us %%4s %%4s %%6u %%6u %%6u\n", plink_maxsnp);
  } else {
    sprintf(tbuf, " CHR %%%us   A1   A2          MAF  NCHROBS\n", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
    sprintf(tbuf, " %%%us %%4s %%4s %%12.4g %%8d\n", plink_maxsnp);
  }
  if (ferror(outfile)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->x_code);
    is_y = (chrom_idx == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end);
    while (marker_uidx < chrom_end) {
      reverse = IS_SET(marker_reverse, marker_uidx);
      major_ptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      minor_ptr = marker_allele_ptrs[marker_uidx * 2];
      if (freq_counts || freqx) {
	if (is_x) {
	  missing_ct = indiv_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx] + hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	} else if (is_haploid) {
	  if (is_y) {
	    missing_ct = indiv_f_male_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  } else {
	    missing_ct = indiv_f_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  }
	} else {
	  missing_ct = indiv_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]);
	}
	if (freqx) {
	  bufptr = chrom_name_write(tbuf2, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms);
	  fwrite(tbuf2, 1, bufptr - tbuf2, outfile);
	  fprintf(outfile, "\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\n", &(marker_ids[marker_uidx * max_marker_id_len]), minor_ptr, major_ptr, reverse? hh_cts[marker_uidx] : ll_cts[marker_uidx], lh_cts[marker_uidx], reverse? ll_cts[marker_uidx] : hh_cts[marker_uidx], reverse? haph_cts[marker_uidx] : hapl_cts[marker_uidx], reverse? hapl_cts[marker_uidx] : haph_cts[marker_uidx], missing_ct);
	} else {
	  bufptr = width_force(4, tbuf2, chrom_name_write(tbuf2, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
	  fwrite(tbuf2, 1, bufptr - tbuf2, outfile);
	  fprintf(outfile, tbuf, &(marker_ids[marker_uidx * max_marker_id_len]), minor_ptr, major_ptr, 2 * ll_cts[marker_uidx] + lh_cts[marker_uidx] + hapl_cts[marker_uidx], 2 * hh_cts[marker_uidx] + lh_cts[marker_uidx] + haph_cts[marker_uidx], missing_ct);
	}
      } else {
	bufptr = width_force(4, tbuf2, chrom_name_write(tbuf2, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
	fwrite(tbuf2, 1, bufptr - tbuf2, outfile);
	fprintf(outfile, tbuf, &(marker_ids[marker_uidx * max_marker_id_len]), minor_ptr, major_ptr, (1.0 - set_allele_freqs[marker_uidx]) * (1 + SMALL_EPSILON), 2 * (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]) + hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
      }
      if (ferror(outfile)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
      marker_uidx = next_unset(marker_exclude, marker_uidx + 1, chrom_end);
    }
  }
  if (fclose_null(&outfile)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Allele frequencies written to %s.\n", outname);
  logprintb();
  while (0) {
  write_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t write_stratified_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uint32_t indiv_f_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uint32_t indiv_f_male_ct, uintptr_t* marker_reverse, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uint32_t* cur_cluster_map = cluster_map;
  uint32_t* cur_cluster_starts = cluster_starts;
  uint32_t* cluster_map_nonmale = NULL;
  uint32_t* cluster_starts_nonmale = NULL;
  uint32_t* cluster_map_male = NULL;
  uint32_t* cluster_starts_male = NULL;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  uint32_t cslen = 10;
  int32_t retval = 0;
  uint32_t cur_cts[4];
  uintptr_t* readbuf;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  char* csptr;
  char* col_2_start;
  char* wptr_start;
  char* wptr;
  char* sptr;
  uintptr_t chrom_end;
  uintptr_t marker_uidx;
  int32_t chrom_idx;
  uintptr_t clidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t clmpos;
  uint32_t a1_obs;
  uint32_t tot_obs;
  uint32_t uii;
  if (wkspace_alloc_ul_checked(&readbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  if ((indiv_ct > indiv_f_ct) && (!nonfounders)) {
    if (wkspace_alloc_ui_checked(&cur_cluster_starts, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cur_cluster_map, indiv_f_ct * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cur_cluster_starts[0] = 0;
    uiptr = cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cluster_map[cluster_starts[clidx + 1]]);
      do {
	uii = *uiptr;
	if (IS_SET(founder_info, uii)) {
          cur_cluster_map[clmpos++] = uii;
	}
      } while ((++uiptr) < uiptr2);
      cur_cluster_starts[clidx + 1] = clmpos;
    }
  }
  chrom_idx = chrom_info_ptr->x_code;
  if ((chrom_idx != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
    if (wkspace_alloc_ui_checked(&cluster_starts_nonmale, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_map_nonmale, (indiv_f_ct - indiv_f_male_ct) * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cluster_starts_nonmale[0] = 0;
    uiptr = cur_cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
      while (uiptr < uiptr2) {
	uii = *uiptr++;
	if (!IS_SET(sex_male, uii)) {
          cluster_map_nonmale[clmpos++] = uii;
	}
      }
      cluster_starts_nonmale[clidx + 1] = clmpos;
    }
  }
  chrom_idx = chrom_info_ptr->y_code;
  if (cluster_map_nonmale || ((chrom_idx != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_idx))) {
    if (wkspace_alloc_ui_checked(&cluster_starts_male, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_map_male, indiv_f_male_ct * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cluster_starts_male[0] = 0;
    uiptr = cur_cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
      while (uiptr < uiptr2) {
	uii = *uiptr++;
	if (IS_SET(sex_male, uii)) {
          cluster_map_male[clmpos++] = uii;
	}
      }
      cluster_starts_male[clidx + 1] = clmpos;
    }
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_stratified_freqs_ret_OPEN_FAIL;
  }
  sprintf(tbuf, " CHR %%%ds     CLST   A1   A2      MAF    MAC  NCHROBS\n", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  if (wkspace_alloc_c_checked(&csptr, 2 * max_marker_allele_len + 16)) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  memset(csptr, 32, 10);
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->x_code);
    is_y = (chrom_idx == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
    marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end);
    if (marker_uidx >= chrom_end) {
      continue;
    }
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto write_stratified_freqs_ret_READ_FAIL;
    }
    col_2_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms));
    *col_2_start++ = ' ';
    do {
      sptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      uii = strlen(sptr);
      wptr_start = memseta(col_2_start, 32, plink_maxsnp - uii);
      wptr_start = memcpyax(wptr_start, sptr, uii, ' ');
      sptr = marker_allele_ptrs[marker_uidx * 2];
      wptr = fw_strcpy(4, sptr, &(csptr[1]));
      *wptr++ = ' ';
      sptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      wptr = fw_strcpy(4, sptr, wptr);
      *wptr++ = ' ';
      cslen = (uintptr_t)(wptr - csptr);

      if (fread(readbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto write_stratified_freqs_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf((unsigned char*)readbuf, unfiltered_indiv_ct);
      }
      if (is_x) {
	uiptr = cluster_map_nonmale;
	uiptr2 = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr3 = &(cluster_map_nonmale[cluster_starts_nonmale[clidx + 1]]);
	  while (uiptr < uiptr3) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs = 2 * cur_cts[0] + cur_cts[2];
	  tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  fill_uint_zero(cur_cts, 4);
	  uiptr3 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr2 < uiptr3) {
	    uii = *uiptr2++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs += cur_cts[0];
	  tot_obs += cur_cts[0] + cur_cts[3];
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else if (is_y) {
	uiptr = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr2 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr < uiptr2) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  if (is_haploid) {
	    a1_obs = cur_cts[0];
	    tot_obs = cur_cts[0] + cur_cts[3];
	  } else {
	    a1_obs = 2 * cur_cts[0] + cur_cts[2];
	    tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  }
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else {
        uiptr = cur_cluster_map;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
	  while (uiptr < uiptr2) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  if (is_haploid) {
	    a1_obs = cur_cts[0];
	    tot_obs = cur_cts[0] + cur_cts[3];
	  } else {
	    a1_obs = 2 * cur_cts[0] + cur_cts[2];
	    tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  }
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto write_stratified_freqs_ret_READ_FAIL;
	}
      }
    } while (marker_uidx < chrom_end);
  }
  if (fclose_null(&outfile)) {
    goto write_stratified_freqs_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Cluster-stratified allele frequencies written to %s.\n", outname);
  logprintb();
  while (0) {
  write_stratified_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_stratified_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_stratified_freqs_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_stratified_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

void calc_marker_reverse_bin(uintptr_t* marker_reverse, uintptr_t* marker_exclude, uint32_t unfiltered_marker_ct, uint32_t marker_ct, double* set_allele_freqs) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  double dxx;
  do {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    for (; marker_uidx < marker_uidx_stop; marker_uidx++) {
      dxx = set_allele_freqs[marker_uidx];
      if (dxx < 0.5) {
	SET_BIT(marker_reverse, marker_uidx);
	set_allele_freqs[marker_uidx] = 1.0 - dxx;
      }
    }
  } while (markers_done < marker_ct);
}

int32_t hardy_report_write_line(FILE* outfile, char* prefix_buf, uint32_t prefix_len, uint32_t reverse, uint32_t ll_ct, uint32_t lh_ct, uint32_t hh_ct, char* midbuf_ptr, double pvalue) {
  char wbuf[48];
  char* cptr;
  uint32_t denom;
  double drecip;
  double minor_freq;
  fwrite(prefix_buf, 1, prefix_len, outfile);
  if (reverse) {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, hh_ct, '/'), lh_ct, '/'), ll_ct);
  } else {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, ll_ct, '/'), lh_ct, '/'), hh_ct);
  }
  cptr = fw_strcpyn(20, cptr - wbuf, wbuf, midbuf_ptr);
  *cptr++ = ' ';
  denom = (ll_ct + lh_ct + hh_ct) * 2;
  if (denom) {
    drecip = 1.0 / ((double)denom);
    minor_freq = (2 * ll_ct + lh_ct) * drecip;
    cptr = double_g_writewx4x(double_g_writewx4x(double_g_writewx4x(cptr, (lh_ct * 2) * drecip, 8, ' '), minor_freq * (2 * hh_ct + lh_ct) * drecip * 2, 8, ' '), pvalue, 12, '\n');
  } else {
    cptr = memcpya(cptr, "     nan      nan           NA\n", 31);
  }
  return fwrite_checked(midbuf_ptr, (cptr - midbuf_ptr), outfile);
}

int32_t hardy_report(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, int32_t* hwe_ll_cases, int32_t* hwe_lh_cases, int32_t* hwe_hh_cases, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, uint32_t pheno_nm_ct, uintptr_t* pheno_c, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr) {
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t hwe_midp = hwe_modifier & HWE_MIDP;
  int32_t retval = 0;
  uint32_t skip_chrom = 0;
  uint32_t pct = 0;
  uint32_t prefix_len;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t report_type;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t reverse;
  double* p_values;
  char* writebuf;
  char* cptr0;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5;
  if (pheno_nm_ct) {
    report_type = pheno_c? 0 : 1;
  } else {
    report_type = 2;
  }
  uii = report_type? 1 : 3;
  if (wkspace_alloc_d_checked(&p_values, uii * marker_ct * sizeof(double)) ||
      wkspace_alloc_c_checked(&writebuf, 2 * max_marker_allele_len + MAXLINELEN)) {
    goto hardy_report_ret_NOMEM;
  }

  // todo: multithread?
  if (report_type) {
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      p_values[marker_idx] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], hwe_midp);
    }
  } else {
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      p_values[marker_idx * 3] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], hwe_midp);
      p_values[marker_idx * 3 + 1] = SNPHWE2(hwe_lh_cases[marker_uidx], hwe_ll_cases[marker_uidx], hwe_hh_cases[marker_uidx], hwe_midp);
      p_values[marker_idx * 3 + 2] = SNPHWE2(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_midp);
    }
  }
  marker_uidx = 0;
  marker_idx = 0;

  memcpy(outname_end, ".hwe", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto hardy_report_ret_OPEN_FAIL;
  }
  sprintf(logbuf, "Writing Hardy-Weinberg report to %s...", outname);
  logprintb();
  fputs(" 0%", stdout);
  sprintf(writebuf, " CHR %%%us     TEST   A1   A2                 GENO   O(HET)   E(HET)            P \n", plink_maxsnp);
  fprintf(outfile, writebuf, "SNP");
 
  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
  skip_chrom = (is_haploid && (!is_x)) || is_mt;
  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
  *cptr0++ = ' ';
  cptr = &(cptr0[10 + plink_maxsnp]);
  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
  if (report_type) {
    if (report_type == 1) {
      memcpy(&(cptr0[plink_maxsnp]), "  ALL(QT)           ", 20);
    } else {
      memcpy(&(cptr0[plink_maxsnp]), "  ALL(NP)           ", 20);
    }
    cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  skip_chrom = (is_haploid && (!is_x)) || is_mt;
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
	  *cptr0++ = ' ';
	  cptr = &(cptr0[10 + plink_maxsnp]);
	  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
	  if (report_type == 1) {
	    memcpy(&(cptr0[plink_maxsnp]), "  ALL(QT)           ", 20);
	  } else {
	    memcpy(&(cptr0[plink_maxsnp]), "  ALL(NP)           ", 20);
	  }
	  cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
	}
        if (skip_chrom) {
	  continue;
	}
	fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr0);
	reverse = IS_SET(marker_reverse, marker_uidx);
	cptr3 = marker_allele_ptrs[2 * marker_uidx];
	cptr4 = marker_allele_ptrs[2 * marker_uidx + 1];
	cptr5 = fw_strcpy(4, cptr3, cptr);
	*cptr5 = ' ';
	cptr5 = fw_strcpy(4, cptr4, &(cptr5[1]));
	*cptr5 = ' ';
	prefix_len = 1 + (cptr5 - writebuf);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[marker_idx])) {
	  goto hardy_report_ret_WRITE_FAIL;
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
  } else {
    memset(&(cptr0[plink_maxsnp]), 32, 20);
    cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  skip_chrom = (is_haploid && (!is_x)) || is_mt;
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
          memset(&(cptr0[plink_maxsnp]), 32, 20);
	  cptr = &(cptr0[10 + plink_maxsnp]);
	  cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
	  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
	}
	if (skip_chrom) {
	  continue;
	}
	fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr0);
	memcpy(&(cptr0[4 + plink_maxsnp]), "  ALL", 5);
	reverse = IS_SET(marker_reverse, marker_uidx);
	cptr3 = marker_allele_ptrs[2 * marker_uidx];
	cptr4 = marker_allele_ptrs[2 * marker_uidx + 1];
	cptr5 = fw_strcpy(4, cptr3, cptr);
	*cptr5 = ' ';
	cptr5 = fw_strcpy(4, cptr4, &(cptr5[1]));
	*cptr5 = ' ';
	prefix_len = 1 + (cptr5 - writebuf);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[3 * marker_idx])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[7 + plink_maxsnp]), "FF", 2);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_cases[marker_uidx], hwe_lh_cases[marker_uidx], hwe_hh_cases[marker_uidx], cptr2, p_values[3 * marker_idx + 1])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[4 + plink_maxsnp]), "UN", 2);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_lls[marker_uidx], hwe_lhs[marker_uidx], hwe_hhs[marker_uidx], cptr2, p_values[3 * marker_idx + 2])) {
	  goto hardy_report_ret_WRITE_FAIL;
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
  }
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");

  while (0) {
  hardy_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  hardy_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  hardy_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}


uint32_t enforce_hwe_threshold(double hwe_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, Chrom_info* chrom_info_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t hwe_all = hwe_modifier & HWE_THRESH_ALL;
  uint32_t hwe_thresh_midp = hwe_modifier & HWE_THRESH_MIDP;
  uint32_t min_obs = 0xffffffffU;
  uint32_t max_obs = 0;
  int32_t mt_code = chrom_info_ptr->mt_code;
  uint32_t mt_start = 0;
  uint32_t mt_end = 0;
  uint32_t markers_done;
  uint32_t cur_obs;
  hwe_thresh *= 1 + SMALL_EPSILON;
  if (hwe_all) {
    hwe_lhs = hwe_lh_allfs;
    hwe_lls = hwe_ll_allfs;
    hwe_hhs = hwe_hh_allfs;
  }
  if ((mt_code != -1) && is_set(chrom_info_ptr->chrom_mask, mt_code)) {
    mt_start = chrom_info_ptr->chrom_start[(uint32_t)mt_code];
    mt_end = chrom_info_ptr->chrom_end[(uint32_t)mt_code];
  }
  if (hwe_thresh_midp) {
    for (markers_done = 0; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      if ((marker_uidx < mt_end) && (marker_uidx >= mt_start)) {
        continue;
      }
      if (SNPHWE_midp_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
      cur_obs = hwe_lhs[marker_uidx] + hwe_lls[marker_uidx] + hwe_hhs[marker_uidx];
      if (cur_obs < min_obs) {
	min_obs = cur_obs;
      }
      if (cur_obs > max_obs) {
	max_obs = cur_obs;
      }
    }
  } else {
    for (markers_done = 0; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      if ((marker_uidx < mt_end) && (marker_uidx >= mt_start)) {
        continue;
      }
      if (SNPHWE_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
      cur_obs = hwe_lhs[marker_uidx] + hwe_lls[marker_uidx] + hwe_hhs[marker_uidx];
      if (cur_obs < min_obs) {
	min_obs = cur_obs;
      }
      if (cur_obs > max_obs) {
	max_obs = cur_obs;
      }
    }
  }
  if (((uint64_t)max_obs) * 9 > ((uint64_t)min_obs) * 10) {
    logprint("Warning: --hwe observation counts vary by more than 10%.  Consider using\n--geno, and/or applying different p-value thresholds to distinct subsets of\nyour data.\n");
  }
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed due to Hardy-Weinberg exact test (--hwe).\n");
    return 1;
  }
  sprintf(logbuf, "%u variant%s removed due to Hardy-Weinberg exact test (--hwe).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

uint32_t enforce_maf_threshold(double min_maf, double max_maf, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, double* set_allele_freqs) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  double dxx;
  min_maf *= 1 - SMALL_EPSILON;
  max_maf *= 1 + SMALL_EPSILON;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      dxx = get_maf(set_allele_freqs[marker_uidx]);
      if ((dxx < min_maf) || (dxx > max_maf)) {
        SET_BIT(marker_exclude, marker_uidx);
        removed_ct++;
      }
    } while (++marker_uidx < marker_uidx_stop);
  }
  if (marker_ct == removed_ct) {
    logprint("Error: All variants removed due to MAF threshold(s) (--maf/--max-maf).\n");
    return 1;
  }
  sprintf(logbuf, "%u variant%s removed due to MAF threshold(s) (--maf/--max-maf).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
  return 0;
}

void enforce_min_bp_space(int32_t min_bp_space, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, Chrom_info* chrom_info_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - (uint32_t)(*marker_exclude_ct_ptr);
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t removed_ct = 0;
  uint32_t chrom_end = 0;
  uint32_t marker_uidx = next_unset(marker_exclude, 0, unfiltered_marker_ct);
  uint32_t chrom_fo_idx_p1 = 0;
  uint32_t marker_uidx_stop;
  int32_t last_pos;
  int32_t cur_pos;
  for (chrom_fo_idx_p1 = 1; chrom_fo_idx_p1 <= chrom_ct; chrom_fo_idx_p1++) {
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx_p1];
    if (marker_uidx >= chrom_end) {
      continue;
    }
    last_pos = -2147483647;
    do {
      marker_uidx_stop = next_set(marker_exclude, marker_uidx, chrom_end);
      do {
        cur_pos = marker_pos[marker_uidx];
        if (cur_pos < last_pos + min_bp_space) {
          SET_BIT(marker_exclude, marker_uidx);
	  removed_ct++;
	} else {
	  last_pos = cur_pos;
	}
      } while (++marker_uidx < marker_uidx_stop);
      marker_uidx = next_unset(marker_exclude, marker_uidx, unfiltered_marker_ct);
    } while (marker_uidx < chrom_end);
  }
  sprintf(logbuf, "--bp-space: %u variant%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", marker_ct - removed_ct);
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
}

void calc_marker_weights(double exponent, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, double* marker_weights) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  do {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      if (marker_weights[marker_uidx] < 0.0) {
	marker_weights[marker_uidx] = calc_wt_mean(exponent, lh_cts[marker_uidx], ll_cts[marker_uidx], hh_cts[marker_uidx]);
      }
    } while (++marker_uidx < marker_uidx_stop);
  } while (markers_done < marker_ct);
}

int32_t load_ax_alleles(Two_col_params* axalleles, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uint32_t is_a2) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = axalleles->skipchar;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  uint32_t colid_first = (axalleles->colid < axalleles->colx);
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t slen;
  char* sorted_marker_ids;
  char* colid_ptr;
  char* colx_ptr;
  uint32_t* marker_id_map;
  uint32_t colmin;
  uint32_t coldiff;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  char cc;
  int32_t retval;
  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto load_ax_alleles_ret_1;
  }
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto load_ax_alleles_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_ax_alleles_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  retval = open_and_skip_first_lines(&infile, axalleles->fname, loadbuf, loadbuf_size, axalleles->skip);
  if (colid_first) {
    colmin = axalleles->colid - 1;
    coldiff = axalleles->colx - axalleles->colid;
  } else {
    colmin = axalleles->colx - 1;
    coldiff = axalleles->colid - axalleles->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Pathologically long line in %s.\n", axalleles->fname);
        logprintb();
	goto load_ax_alleles_ret_INVALID_FORMAT;
      } else {
        goto load_ax_alleles_ret_NOMEM;
      }
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      if (colmin) {
        colid_ptr = next_item_mult(colid_ptr, colmin);
      }
      colx_ptr = next_item_mult(colid_ptr, coldiff);
    } else {
      colx_ptr = colid_ptr;
      if (colmin) {
	colx_ptr = next_item_mult(colx_ptr, colmin);
      }
      colid_ptr = next_item_mult(colx_ptr, coldiff);
    }
    slen = strlen_se(colid_ptr);
    sorted_idx = bsearch_str(colid_ptr, slen, sorted_marker_ids, max_marker_id_len, marker_ct);
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      colid_ptr[slen] = '\0';
      sprintf(logbuf, "Error: Duplicate variant %s in --a%c-allele file.\n", colid_ptr, is_a2? '2' : '1');
      logprintb();
      goto load_ax_alleles_ret_INVALID_FORMAT;
    }
    SET_BIT(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    slen = strlen_se(colx_ptr);
    colx_ptr[slen] = '\0';
    if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + is_a2])) {
      if (IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        CLEAR_BIT(marker_reverse, marker_uidx);
      }
    } else if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + 1 - is_a2])) {
      if (!IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        SET_BIT(marker_reverse, marker_uidx);
      }
    } else if (marker_allele_ptrs[marker_uidx * 2 + is_a2] == missing_geno_ptr) {
      if (allele_reset(&(marker_allele_ptrs[marker_uidx * 2 + is_a2]), colx_ptr, slen)) {
	goto load_ax_alleles_ret_NOMEM;
      }
      if (slen >= max_marker_allele_len) {
	max_marker_allele_len = slen + 1;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
        set_allele_freqs[marker_uidx] = 1.0 - set_allele_freqs[marker_uidx];
        CLEAR_BIT(marker_reverse, marker_uidx);
      }
    } else {
      colid_ptr[slen] = '\0';
      sprintf(logbuf, "Warning: Impossible A%c allele assignment for variant %s.\n", is_a2? '2' : '1', colid_ptr);
      logprintb();
    }
  }
  if (!feof(infile)) {
    goto load_ax_alleles_ret_READ_FAIL;
  }
  *max_marker_allele_len_ptr = max_marker_allele_len;
  while (0) {
  load_ax_alleles_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_ax_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_ax_alleles_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_ax_alleles_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

void swap_reversed_marker_alleles(uintptr_t unfiltered_marker_ct, uintptr_t* marker_reverse, char** marker_allele_ptrs) {
  uintptr_t marker_uidx = 0;
  char* swap_ptr;
  while (1) {
    next_set_ul_ck(marker_reverse, &marker_uidx, unfiltered_marker_ct);
    if (marker_uidx == unfiltered_marker_ct) {
      return;
    }
    swap_ptr = marker_allele_ptrs[marker_uidx * 2];
    marker_allele_ptrs[marker_uidx * 2] = marker_allele_ptrs[marker_uidx * 2 + 1];
    marker_allele_ptrs[marker_uidx * 2 + 1] = swap_ptr;
    marker_uidx++;
  }
}

int32_t write_snplist(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t list_23_indels) {
  FILE* outfile;
  uintptr_t marker_uidx = 0;
  uintptr_t markers_done = 0;
  int32_t retval = 0;
  const char* a0ptr = g_missing_geno_ptr;
  const char* adptr = &(g_one_char_strs[136]); // "D"
  const char* aiptr = &(g_one_char_strs[146]); // "I"
  char* a1ptr;
  char* a2ptr;
  char* cptr;
  char* cptr_end;
  uintptr_t marker_uidx_stop;
  if (!list_23_indels) {
    memcpy(outname_end, ".snplist", 9);
  } else {
    memcpy(outname_end, ".indel", 7);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_snplist_ret_OPEN_FAIL;
  }
  if (!list_23_indels) {
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      markers_done += marker_uidx_stop - marker_uidx;
      cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      cptr_end = &(marker_ids[marker_uidx_stop * max_marker_id_len]);
      marker_uidx = marker_uidx_stop;
      do {
	fputs(cptr, outfile);
	if (putc_checked('\n', outfile)) {
          goto write_snplist_ret_WRITE_FAIL;
	}
        cptr = &(cptr[max_marker_id_len]);
      } while (cptr < cptr_end);
    } while (markers_done < marker_ct);
  } else {
    for (; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      a1ptr = marker_allele_ptrs[2 * marker_uidx];
      a2ptr = marker_allele_ptrs[2 * marker_uidx + 1];
      if ((a1ptr != adptr) && (a1ptr != aiptr)) {
        if ((a1ptr != a0ptr) || ((a2ptr != adptr) && (a2ptr != aiptr))) {
	  continue;
	}
      } else if ((a2ptr != adptr) && (a2ptr != aiptr) && (a2ptr != a0ptr)) {
	continue;
      }
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      if (putc_checked('\n', outfile)) {
	goto write_snplist_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto write_snplist_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "List of %svariant IDs written to %s.\n", list_23_indels? "indel " : "" , outname);
  logprintb();
  while (0) {
  write_snplist_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_snplist_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

static inline uint32_t are_marker_pos_needed(uint64_t calculation_type, char* cm_map_fname, char* set_fname, uint32_t min_bp_space, uint32_t genome_skip_write, uint32_t ld_modifier, uint32_t epi_modifier) {
  return (calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_GENOME | CALC_HOMOZYG | CALC_LD_PRUNE | CALC_REGRESS_PCS | CALC_MODEL | CALC_GLM | CALC_CLUMP)) || cm_map_fname || set_fname || min_bp_space || genome_skip_write || ((calculation_type & CALC_LD) && (!(ld_modifier & LD_MATRIX_SHAPEMASK))) || ((calculation_type & CALC_EPI) && (epi_modifier & EPI_FAST_CASE_ONLY));
}

static inline uint32_t are_marker_cms_needed(uint64_t calculation_type, char* cm_map_fname, Two_col_params* update_cm) {
  if (calculation_type & (CALC_MAKE_BED | CALC_RECODE)) {
    if (cm_map_fname || update_cm) {
      return MARKER_CMS_FORCED;
    } else {
      return MARKER_CMS_OPTIONAL;
    }
  }
  return 0;
}

static inline uint32_t are_marker_alleles_needed(uint64_t calculation_type, char* freqname, Homozyg_info* homozyg_ptr, Two_col_params* a1alleles, Two_col_params* a2alleles, uint32_t ld_modifier, uint32_t snp_only, uint32_t clump_modifier) {
  return (freqname || (calculation_type & (CALC_FREQ | CALC_HARDY | CALC_MAKE_BED | CALC_RECODE | CALC_REGRESS_PCS | CALC_MODEL | CALC_GLM | CALC_LASSO | CALC_LIST_23_INDELS | CALC_EPI | CALC_TESTMISHAP)) || ((calculation_type & CALC_HOMOZYG) && (homozyg_ptr->modifier & HOMOZYG_GROUP_VERBOSE)) || ((calculation_type & CALC_LD) && (!(ld_modifier & LD_MATRIX_SHAPEMASK))) || a1alleles || a2alleles || snp_only || (clump_modifier & (CLUMP_VERBOSE | CLUMP_BEST)));
}

inline int32_t relationship_or_ibc_req(uint64_t calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

inline int32_t distance_wt_req(uint64_t calculation_type, char* read_dists_fname, uint32_t dist_calc_type) {
  return (((calculation_type & CALC_DISTANCE) || ((!read_dists_fname) && ((calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))))) && (!(dist_calc_type & DISTANCE_FLAT_MISSING)));
}

int32_t plink(char* outname, char* outname_end, char* pedname, char* mapname, char* famname, char* cm_map_fname, char* cm_map_chrname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* keepfamname, char* removefamname, char* filtername, char* freqname, char* read_dists_fname, char* read_dists_id_fname, char* evecname, char* mergename1, char* mergename2, char* mergename3, char* makepheno_str, char* phenoname_str, Two_col_params* a1alleles, Two_col_params* a2alleles, char* recode_allele_name, char* covar_fname, char* update_alleles_fname, char* read_genome_fname, Two_col_params* update_chr, Two_col_params* update_cm, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* loop_assoc_fname, char* flip_fname, char* flip_subset_fname, char* filtervals_flattened, char* condition_mname, char* condition_fname, char* filter_attrib_fname, char* filter_attrib_liststr, char* filter_attrib_indiv_fname, char* filter_attrib_indiv_liststr, double thin_keep_prob, uint32_t min_bp_space, uint32_t mfilter_col, uint32_t filter_binary, uint32_t fam_cols, int32_t missing_pheno, char* output_missing_pheno, uint32_t mpheno_col, uint32_t pheno_modifier, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, double check_sex_fthresh, double check_sex_mthresh, double exponent, double min_maf, double max_maf, double geno_thresh, double mind_thresh, double hwe_thresh, double rel_cutoff, double tail_bottom, double tail_top, uint64_t misc_flags, uint64_t calculation_type, uint32_t rel_calc_type, uint32_t dist_calc_type, uintptr_t groupdist_iters, uint32_t groupdist_d, uintptr_t regress_iters, uint32_t regress_d, uintptr_t regress_rel_iters, uint32_t regress_rel_d, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_covr, int32_t ibc_type, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t splitx_bound1, uint32_t splitx_bound2, uint32_t ppc_gap, uint32_t sex_missing_pheno, uint32_t hwe_modifier, uint32_t genome_modifier, double genome_min_pi_hat, double genome_max_pi_hat, Homozyg_info* homozyg_ptr, Cluster_info* cluster_ptr, uint32_t neighbor_n1, uint32_t neighbor_n2, Set_info* sip, Ld_info* ldip, Epi_info* epi_ip, Clump_info* clump_ip, uint32_t regress_pcs_modifier, uint32_t max_pcs, uint32_t recode_modifier, uint32_t allelexxxx, uint32_t merge_type, uint32_t indiv_sort, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* snps_range_list_ptr, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uint32_t mwithin_col, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t gxe_mcovar, Aperm_info* apip, uint32_t mperm_save, uint32_t ibs_test_perms, uint32_t perm_batch_size, double lasso_h2, double lasso_minlambda, Range_list* lasso_select_covars_range_list_ptr, uint32_t testmiss_modifier, uint32_t testmiss_mperm_val, Ll_str** file_delete_list_ptr) {
  FILE* bedfile = NULL;
  FILE* famfile = NULL;
  FILE* phenofile = NULL;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t* marker_exclude = NULL;
  uintptr_t max_marker_id_len = 0;
  // set_allele_freqs = frequency of allele corresponding to set bits in .bed
  //   (i.e. A2), or frequency of MAJOR allele in middle of text loading.
  double* set_allele_freqs = NULL;
  uintptr_t topsize = 0;
  uintptr_t unfiltered_indiv_ct = 0;
  uintptr_t* indiv_exclude = NULL;
  uintptr_t indiv_exclude_ct = 0;
  uint32_t* indiv_sort_map = NULL;
  uintptr_t* founder_info = NULL;
  uintptr_t* sex_nm = NULL;
  uintptr_t* sex_male = NULL;
  uint32_t genome_skip_write = (cluster_ptr->ppc != 0.0) && (!(calculation_type & CALC_GENOME)) && (!read_genome_fname);
  uint32_t marker_pos_needed = are_marker_pos_needed(calculation_type, cm_map_fname, sip->fname, min_bp_space, genome_skip_write, ldip->modifier, epi_ip->modifier);
  uint32_t marker_cms_needed = are_marker_cms_needed(calculation_type, cm_map_fname, update_cm);
  uint32_t marker_alleles_needed = are_marker_alleles_needed(calculation_type, freqname, homozyg_ptr, a1alleles, a2alleles, ldip->modifier, (misc_flags / MISC_SNPS_ONLY) & 1, clump_ip->modifier);
  uint32_t zero_extra_chroms = (misc_flags / MISC_ZERO_EXTRA_CHROMS) & 1;
  uint32_t uii = 0;
  int64_t llxx = 0;
  uint32_t nonfounders = (misc_flags / MISC_NONFOUNDERS) & 1;
  uint32_t pheno_all = pheno_modifier & PHENO_ALL;
  char* marker_ids = NULL;
  double* marker_cms = NULL;
  // marker_allele_ptrs[2 * i] is id of A1 (usually minor) allele at marker i
  // marker_allele_ptrs[2 * i + 1] is id of A2 allele
  // Single-character allele names point to g_one_char_strs[]; otherwise
  // string allocation occurs on the heap.
  char** marker_allele_ptrs = NULL;
  uintptr_t max_marker_allele_len = 2; // includes trailing null
  uintptr_t* marker_reverse = NULL;
  int32_t retval = 0;
  uint32_t map_is_unsorted = 0;
  uint32_t map_cols = 3;
  uint32_t affection = 0;
  uintptr_t* pheno_nm = NULL;
  uintptr_t* pheno_nm_datagen = NULL; // --make-bed/--recode/--write-covar only
  uintptr_t* orig_pheno_nm = NULL; // --all-pheno + --pheno-merge
  uintptr_t* pheno_c = NULL;
  uintptr_t* orig_pheno_c = NULL;
  double* pheno_d = NULL;
  double* orig_pheno_d = NULL;
  double* marker_weights = NULL;
  uint32_t marker_weight_sum = 0;
  uint32_t* marker_weights_i = NULL;
  char* person_ids = NULL;
  uintptr_t max_person_id_len = 4;
  char* paternal_ids = NULL;
  uintptr_t max_paternal_id_len = 2;
  char* maternal_ids = NULL;
  uintptr_t max_maternal_id_len = 2;
  unsigned char* wkspace_mark = NULL;
  uintptr_t cluster_ct = 0;
  uint32_t* cluster_map = NULL; // unfiltered indiv IDs
  // index for cluster_map, length (cluster_ct + 1)
  // cluster_starts[n+1] - cluster_starts[n] = length of cluster n (0-based)
  uint32_t* cluster_starts = NULL;
  char* cluster_ids = NULL;
  uintptr_t max_cluster_id_len = 2;
  double* mds_plot_dmatrix_copy = NULL;
  uintptr_t* cluster_merge_prevented = NULL;
  double* cluster_sorted_ibs = NULL;
  char* cptr = NULL;
  uint64_t dists_alloc = 0;
  uintptr_t marker_exclude_ct = 0;
  char* pid_list = NULL;
  char* id_list = NULL;
  double missing_phenod = (double)missing_pheno;
  double ci_zt = 0.0;
  uint32_t missing_pheno_len = intlen(missing_pheno);
  uint32_t wt_needed = distance_wt_req(calculation_type, read_dists_fname, dist_calc_type);
  uintptr_t bed_offset = 3;
  uint32_t* marker_pos = NULL;
  uint32_t hh_exists = 0;
  uint32_t pheno_ctrl_ct = 0;
  uintptr_t covar_ct = 0;
  char* covar_names = NULL;
  uintptr_t max_covar_name_len = 0;
  uintptr_t* covar_nm = NULL;
  double* covar_d = NULL;
  uintptr_t* gxe_covar_nm = NULL;
  uintptr_t* gxe_covar_c = NULL;
  uint32_t plink_maxfid = 0;
  uint32_t plink_maxiid = 0;
  unsigned char* wkspace_mark2 = NULL;
  unsigned char* wkspace_mark_precluster = NULL;
  unsigned char* wkspace_mark_postcluster = NULL;
  pthread_t threads[MAX_THREADS];
  uintptr_t unfiltered_indiv_ct4;
  uintptr_t unfiltered_indiv_ctl;
  uint32_t* marker_allele_cts;
  uint32_t* uiptr;
  double* dptr;
  double* dptr2;
  double* rel_ibc;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t ujj;
  uint32_t ukk;
  double dxx;
  char* outname_end2;
  uintptr_t marker_ct;
  uint32_t plink_maxsnp;
  int32_t ii;
  int64_t llyy;
  int32_t* hwe_lls;
  int32_t* hwe_lhs;
  int32_t* hwe_hhs;
  int32_t* hwe_ll_cases;
  int32_t* hwe_lh_cases;
  int32_t* hwe_hh_cases;
  int32_t* hwe_ll_allfs;
  int32_t* hwe_lh_allfs;
  int32_t* hwe_hh_allfs;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uint32_t indiv_male_ct;
  uint32_t indiv_f_ct;
  uint32_t indiv_f_male_ct;
  uint32_t pheno_nm_ct;
  Pedigree_rel_info pri;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx_stop;
  uintptr_t marker_idx;
  uint32_t gender_unk_ct;

  if ((cm_map_fname || update_cm) && (!marker_cms_needed)) {
    sprintf(logbuf, "Error: --%s results would never be used.  (Did you forget --make-bed?)\n", cm_map_fname? "cm-map" : "update-cm");
    logprintb();
    goto plink_ret_INVALID_CMDLINE;
  } else if (update_map && (!marker_pos_needed)) {
    logprint("Error: --update-map results would never be used.  (Did you forget --make-bed?)\n");
    goto plink_ret_INVALID_CMDLINE;
  }
  if (ci_size != 0.0) {
    ci_zt = ltqnorm(1 - (1 - ci_size) / 2);
  }
  if (rel_calc_type & REL_CALC_COV) {
    ibc_type = -1;
  }

  if (calculation_type & CALC_MAKE_BED) {
#if _WIN32
    uii = GetFullPathName(pedname, FNAMESIZE, tbuf, NULL);
    if ((!uii) || (uii > FNAMESIZE))
#else
    if (!realpath(pedname, tbuf))
#endif
    {
      sprintf(logbuf, "Error: Failed to open %s.\n", pedname);
      logprintb();
      goto plink_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    // if file doesn't exist, realpath returns NULL on Linux instead of what
    // the path would be.
#if _WIN32
    uii = GetFullPathName(outname, FNAMESIZE, &(tbuf[FNAMESIZE + 64]), NULL);
    if (uii && (uii <= FNAMESIZE) && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64]))))
#else
    cptr = realpath(outname, &(tbuf[FNAMESIZE + 64]));
    if (cptr && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64]))))
#endif
    {
      logprint("Note: --make-bed input and output filenames match.  Appending '~' to input\nfilenames.\n");
      uii = strlen(pedname);
      memcpy(tbuf, pedname, uii + 1);
      memcpy(&(pedname[uii]), "~", 2);
      if (rename(tbuf, pedname)) {
	logprint("Error: Failed to append '~' to input .bed filename.\n");
	goto plink_ret_OPEN_FAIL;
      }
      uii = strlen(mapname);
      memcpy(tbuf, mapname, uii + 1);
      memcpy(&(mapname[uii]), "~", 2);
      if (rename(tbuf, mapname)) {
	logprint("Error: Failed to append '~' to input .bim filename.\n");
	goto plink_ret_OPEN_FAIL;
      }
      uii = strlen(famname);
      memcpy(tbuf, famname, uii + 1);
      memcpy(&(famname[uii]), "~", 2);
      if (rename(tbuf, famname)) {
	logprint("Error: Failed to append '~' to input .fam filename.\n");
	goto plink_ret_OPEN_FAIL;
      }
    }
  }

  if (calculation_type & CALC_MERGE) {
    if (!(((fam_cols & FAM_COL_13456) == FAM_COL_13456) && (!(misc_flags & MISC_AFFECTION_01)) && (missing_pheno == -9))) {
      logprint("Error: --merge/--bmerge/--merge-list cannot be used with an irregularly\nformatted reference fileset (--no-fid, --no-parents, --no-sex, --no-pheno,\n--1).  Use --make-bed first.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    // Only append -merge to the filename stem if --make-bed or --recode lgen
    // is specified.
    ulii = bed_suffix_conflict(calculation_type, recode_modifier);
    if (ulii) {
      memcpy(outname_end, "-merge", 7);
    }
    retval = merge_datasets(pedname, mapname, famname, outname, ulii? &(outname_end[6]) : outname_end, mergename1, mergename2, mergename3, calculation_type, merge_type, indiv_sort, misc_flags, chrom_info_ptr);
    if (retval || (!(calculation_type & (~CALC_MERGE)))) {
      goto plink_ret_1;
    }
    uljj = (uintptr_t)(outname_end - outname) + (ulii? 6 : 0);
    memcpy(memcpya(pedname, outname, uljj), ".bed", 5);
    memcpy(memcpya(famname, pedname, uljj), ".fam", 5);
    memcpy(memcpya(mapname, pedname, uljj), ".bim", 5);
    if ((calculation_type & CALC_MAKE_BED) && ulii) {
      if (push_ll_str(file_delete_list_ptr, pedname) || push_ll_str(file_delete_list_ptr, famname) || push_ll_str(file_delete_list_ptr, mapname)) {
	goto plink_ret_NOMEM;
      }
    }
  }

  if (fopen_checked(&bedfile, pedname, "rb")) {
    goto plink_ret_OPEN_FAIL;
  }
  if (fopen_checked(&famfile, famname, "rb")) {
    goto plink_ret_OPEN_FAIL;
  }
  // load .bim, count markers, filter chromosomes
  if (update_name) {
    ulii = 0;
    retval = scan_max_strlen(update_name->fname, update_name->colid, update_name->colx, update_name->skip, update_name->skipchar, &max_marker_id_len, &ulii);
    if (retval) {
      goto plink_ret_1;
    }
    if (ulii > 80) {
      // only warn on long new marker ID, since if there's a long old marker ID
      // and no long new one, it's reasonable to infer that the user is fixing
      // the problem, so we shouldn't spam them.
      logprint("Warning: Unusually long new variant ID(s) in --update-name file.  Double-check\nyour file and command-line parameters, and consider changing your naming\nscheme if you encounter memory problems.\n");
    }
    if (ulii > max_marker_id_len) {
      max_marker_id_len = ulii;
    }
  }
  if (!marker_alleles_needed) {
    allelexxxx = 0;
  }
  retval = load_bim(mapname, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &marker_exclude, &set_allele_freqs, &marker_allele_ptrs, &max_marker_allele_len, &marker_ids, chrom_info_ptr, &marker_cms, &marker_pos, freqname, calculation_type, misc_flags, recode_modifier, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, snps_range_list_ptr, &map_is_unsorted, marker_pos_needed, marker_cms_needed, marker_alleles_needed, "bim", ((calculation_type == CALC_MAKE_BED) && (mind_thresh == 1.0) && (geno_thresh == 1.0) && (!update_map) && (!freqname))? NULL : "make-bed");
  if (retval) {
    goto plink_ret_1;
  }
  if (map_is_unsorted & UNSORTED_SPLIT_CHROM) {
    if ((hwe_thresh > 0.0) || (calculation_type & CALC_HARDY)) {
      logprint("Error: --hardy/--hwe cannot be used on a .bim file with split chromosome(s).\nRetry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_FORMAT;
    }
  }

  // load .fam, count indivs
  uii = fam_cols & FAM_COL_6;
  if (uii && phenoname) {
    uii = (pheno_modifier & PHENO_MERGE) && (!makepheno_str);
  }
  if (!uii) {
    pheno_modifier &= ~PHENO_MERGE; // nothing to merge
  }
  if (update_ids_fname) {
    ulii = 0;
    retval = scan_max_fam_indiv_strlen(update_ids_fname, 3, &max_person_id_len);
    if (retval) {
      goto plink_ret_1;
    }
  } else if (update_parents_fname) {
    retval = scan_max_strlen(update_parents_fname, 3, 4, 0, '\0', &max_paternal_id_len, &max_maternal_id_len);
    if (retval) {
      goto plink_ret_1;
    }
  }

  retval = load_fam(famfile, MAXLINELEN, fam_cols, uii, missing_pheno, missing_pheno_len, (misc_flags / MISC_AFFECTION_01) & 1, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &indiv_exclude);
  if (retval) {
    goto plink_ret_1;
  }

  unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;

  if (misc_flags & MISC_MAKE_FOUNDERS_FIRST) {
    if (make_founders(unfiltered_indiv_ct, unfiltered_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, indiv_exclude, founder_info)) {
      goto plink_ret_NOMEM;
    }
  }

  if (pheno_modifier & PHENO_MERGE) {
    if (pheno_all) {
      orig_pheno_nm = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      if (!orig_pheno_nm) {
	goto plink_ret_NOMEM;
      }
      memcpy(orig_pheno_nm, pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
      if (pheno_c) {
	orig_pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
	if (!orig_pheno_c) {
	  goto plink_ret_NOMEM;
	}
	memcpy(orig_pheno_c, pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t));
      } else if (pheno_d) {
	orig_pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
	if (!orig_pheno_d) {
	  goto plink_ret_NOMEM;
	}
	memcpy(orig_pheno_d, pheno_d, unfiltered_indiv_ct * sizeof(double));
      }
    }
  }
  count_genders(sex_nm, sex_male, unfiltered_indiv_ct, indiv_exclude, &uii, &ujj, &gender_unk_ct);
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (gender_unk_ct) {
    sprintf(logbuf, "%" PRIuPTR " variant%s and %" PRIuPTR " %s (%d male%s, %d female%s, %u ambiguous) loaded.\n", marker_ct, (marker_ct == 1)? "" : "s", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s", gender_unk_ct);
    logprintb();
    retval = write_nosex(outname, outname_end, unfiltered_indiv_ct, indiv_exclude, sex_nm, gender_unk_ct, person_ids, max_person_id_len);
    if (retval) {
      goto plink_ret_1;
    }
  } else {
    sprintf(logbuf, "%" PRIuPTR " variant%s and %" PRIuPTR " %s (%d male%s, %d female%s) loaded.\n", marker_ct, (marker_ct == 1)? "" : "s", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s");
    logprintb();
  }

  if (phenoname && fopen_checked(&phenofile, phenoname, "r")) {
    goto plink_ret_OPEN_FAIL;
  }

  if (phenofile || update_ids_fname || update_parents_fname || update_sex_fname || (misc_flags & MISC_TAIL_PHENO)) {
    wkspace_mark = wkspace_base;
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, 0, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink_ret_1;
    }

    if (makepheno_str) {
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_indiv_ct, cptr, max_person_id_len, uiptr, pheno_nm, &pheno_c);
      if (retval) {
	goto plink_ret_1;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_indiv_ct, 0, cptr, max_person_id_len, uiptr, missing_pheno, missing_pheno_len, (misc_flags / MISC_AFFECTION_01) & 1, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d, NULL, 0);
      if (retval) {
	if (retval == LOAD_PHENO_LAST_COL) {
	  logprint(errstr_phenotype_format);
	  logprint("Fewer tokens than expected in line.\n");
	  retval = RET_INVALID_FORMAT;
	  wkspace_reset(wkspace_mark);
	}
	goto plink_ret_1;
      }
    }
    if (misc_flags & MISC_TAIL_PHENO) {
      retval = convert_tail_pheno(unfiltered_indiv_ct, pheno_nm, &pheno_c, &pheno_d, tail_bottom, tail_top, missing_phenod);
      if (retval) {
	goto plink_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if (pheno_c) {
    /*
    if (calculation_type & (CALC_REGRESS_PCS | CALC_REGRESS_PCS_DISTANCE)) {
      sprintf(logbuf, "Error: --regress-pcs%s requires a scalar phenotype.\n", (calculation_type & CALC_REGRESS_PCS_DISTANCE)? "-distance" : "");
      goto plink_ret_INVALID_CMDLINE_2;
    */
    if (calculation_type & (CALC_REGRESS_DISTANCE | CALC_UNRELATED_HERITABILITY | CALC_GXE)) {
      if (calculation_type & CALC_REGRESS_DISTANCE) {
        logprint("Error: --regress-distance calculation requires a scalar phenotype.\n");
      } else if (calculation_type & CALC_UNRELATED_HERITABILITY) {
        logprint("Error: --unrelated-heritability requires a scalar phenotype.\n");
      } else if (calculation_type & CALC_GXE) {
        logprint("Error: --gxe requires a scalar phenotype.\n");
      }
      goto plink_ret_INVALID_CMDLINE;
    }
  } else {
    if ((calculation_type & CALC_MODEL) && (!(model_modifier & MODEL_ASSOC))) {
      logprint("Error: --model requires a case/control phenotype.\n");
      goto plink_ret_INVALID_CMDLINE;
    } else if (calculation_type & CALC_CLUSTER) {
      if (cluster_ptr->modifier & CLUSTER_CC) {
        logprint("Error: --cc requires a case/control phenotype.\n");
        goto plink_ret_INVALID_CMDLINE;
      } else if ((cluster_ptr->max_cases != 0xffffffffU) || (cluster_ptr->max_ctrls != 0xffffffffU)) {
        logprint("Error: --mcc requires a case/control phenotype.\n");
        goto plink_ret_INVALID_CMDLINE;
      }
    } else if ((calculation_type & CALC_EPI) && (epi_ip->modifier & EPI_FAST)) {
      logprint("Error: --fast-epistasis requires a case/control phenotype.\n");
      goto plink_ret_INVALID_CMDLINE;
    } else if (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_CMH | CALC_HOMOG | CALC_TESTMISS)) {
      if (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST)) {
        logprint("Error: --ibs-test and --groupdist calculations require a case/control\nphenotype.\n");
      } else if (calculation_type & CALC_CMH) {
        logprint("Error: --mh and --mh2 require a case/control phenotype.\n");
      } else if (calculation_type & CALC_HOMOG) {
        logprint("Error: --homog requires a case/control phenotype.\n");
      } else {
        logprint("Error: --test-missing requires a case/control phenotype.\n");
      }
      goto plink_ret_INVALID_CMDLINE;
    }
  }

  if ((calculation_type & CALC_GLM) && (!(pheno_modifier & PHENO_ALL))) {
    if (glm_modifier & GLM_LOGISTIC) {
      if (!pheno_c) {
	logprint("Error: --logistic without --all-pheno requires a case/control phenotype.\n");
	goto plink_ret_INVALID_CMDLINE;
      }
    } else if (!pheno_d) {
      logprint("Error: --linear without --all-pheno requires a scalar phenotype.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
  }

  if (cm_map_fname) {
    // need sorted bps, but not marker IDs
    if (map_is_unsorted & (UNSORTED_BP | UNSORTED_SPLIT_CHROM)) {
      logprint("Error: --cm-map requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_FORMAT;
    }
    retval = apply_cm_map(cm_map_fname, cm_map_chrname, unfiltered_marker_ct, marker_exclude, marker_pos, marker_cms, chrom_info_ptr);
    if (retval) {
      goto plink_ret_1;
    }
  }

  uii = update_cm || update_map || update_name || (marker_alleles_needed && (update_alleles_fname || (flip_fname && (!flip_subset_fname)))) || filter_attrib_fname;
  if (uii || extractname || excludename) {
    wkspace_mark = wkspace_base;
    // only permit duplicate marker IDs for --extract/--exclude
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, !uii, 0, strcmp_deref);
    if (retval) {
      goto plink_ret_1;
    }
    // length of sorted list is NOT necessarily equal to unfiltered_marker_ct -
    // marker_exclude_ct for --exclude, since marker_exclude_ct may first
    // change from --update-map or --extract
    ulii = unfiltered_marker_ct - marker_exclude_ct;

    if (update_cm) {
      retval = update_marker_cms(update_cm, cptr, ulii, max_marker_id_len, uiptr, marker_cms);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (update_map) {
      retval = update_marker_pos(update_map, cptr, ulii, max_marker_id_len, uiptr, marker_exclude, &marker_exclude_ct, marker_pos, &map_is_unsorted, chrom_info_ptr);
    }
    if (update_name) {
      retval = update_marker_names(update_name, cptr, ulii, max_marker_id_len, uiptr, marker_ids);
      if (retval) {
	goto plink_ret_1;
      }
      if (update_alleles_fname || (marker_alleles_needed && flip_fname && (!flip_subset_fname)) || extractname || excludename) {
	wkspace_reset(wkspace_mark);
	retval = sort_item_ids(&cptr, &uiptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
	if (retval) {
	  goto plink_ret_1;
	}
	ulii = unfiltered_marker_ct - marker_exclude_ct;
      }
    }
    if (marker_alleles_needed) {
      if (update_alleles_fname) {
        retval = update_marker_alleles(update_alleles_fname, cptr, ulii, max_marker_id_len, uiptr, marker_allele_ptrs, &max_marker_allele_len, outname, outname_end);
        if (retval) {
	  goto plink_ret_1;
        }
      }
      if (flip_fname && (!flip_subset_fname)) {
        retval = flip_strand(flip_fname, cptr, ulii, max_marker_id_len, uiptr, marker_allele_ptrs);
        if (retval) {
	  goto plink_ret_1;
        }
      }
    }
    if (extractname) {
      retval = include_or_exclude(extractname, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
      if (retval) {
        goto plink_ret_1;
      }
    }
    if (excludename) {
      retval = include_or_exclude(excludename, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1);
      if (retval) {
        goto plink_ret_1;
      }
    }
    if (filter_attrib_fname) {
      retval = filter_attrib(filter_attrib_fname, filter_attrib_liststr, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
      if (retval) {
	goto plink_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if (allelexxxx) {
    allelexxxx_recode(allelexxxx, marker_allele_ptrs, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct);
  }

  if (thin_keep_prob != 1.0) {
    if (random_thin_markers(thin_keep_prob, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct)) {
      goto plink_ret_ALL_MARKERS_EXCLUDED;
    }
  }

  if (fseeko(bedfile, 0, SEEK_END)) {
    goto plink_ret_READ_FAIL;
  }
  llxx = ftello(bedfile);
  llyy = llxx - ((uint64_t)unfiltered_indiv_ct4) * unfiltered_marker_ct;
  rewind(bedfile);
  if (llyy == 3LL) {
    // v1.00 or later
    if (fread(tbuf, 1, 3, bedfile) < 3) {
      goto plink_ret_READ_FAIL;
    }
    if (memcmp(tbuf, "l\x1b\x01", 3)) {
      if (memcmp(tbuf, "l\x1b", 3)) {
	if ((tbuf[0] == '#') || (!memcmp(tbuf, "chr", 3))) {
          logprint("Error: Invalid header bytes in PLINK 1 .bed file.  (Is this a UCSC Genome\nBrowser BED file instead?)\n");
	} else {
	  logprint("Error: Invalid header bytes in PLINK 1 .bed file.\n");
	}
	goto plink_ret_INVALID_FORMAT;
      }
      bed_offset = 2;
    }
  } else if (llyy == 1LL) {
    // v0.99
    if (fread(tbuf, 1, 1, bedfile) != 1) {
      goto plink_ret_READ_FAIL;
    }
    if (*tbuf == '\x01') {
      bed_offset = 1;
    } else if (*tbuf == '\0') {
      bed_offset = 2;
    } else {
      logprint("Error: Invalid header bytes in .bed file.\n");
      goto plink_ret_INVALID_FORMAT;
    }
  } else if (llyy != 0LL) {
    sprintf(logbuf, "Error: Invalid .bed file size (expected %" PRIu64 " bytes).\n", 3LLU + ((uint64_t)unfiltered_indiv_ct4) * unfiltered_marker_ct);
    goto plink_ret_INVALID_FORMAT_2;
  } else {
    // pre-0.99, no magic number, indiv-major
    bed_offset = 2;
  }
  if (bed_offset == 2) {
    strcpy(outname_end, ".bed.tmp"); // not really temporary
    logprint("Individual-major .bed file detected.  Transposing to SNP-major form.\n");
    fclose(bedfile);
    retval = indiv_major_to_snp_major(pedname, outname, unfiltered_marker_ct);
    if (retval) {
      goto plink_ret_1;
    }
    strcpy(pedname, outname);
    if (fopen_checked(&bedfile, pedname, "rb")) {
      goto plink_ret_OPEN_FAIL;
    }
    bed_offset = 3;
  }

  if (update_ids_fname || update_parents_fname || update_sex_fname || keepname || keepfamname || removename || removefamname || filter_attrib_indiv_fname || om_ip->marker_fname || filtername) {
    wkspace_mark = wkspace_base;
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink_ret_1;
    }
    ulii = unfiltered_indiv_ct - indiv_exclude_ct;
    if (update_ids_fname) {
      retval = update_indiv_ids(update_ids_fname, cptr, ulii, max_person_id_len, uiptr, person_ids);
      if (retval) {
	goto plink_ret_1;
      }
    } else {
      if (update_parents_fname) {
	retval = update_indiv_parents(update_parents_fname, cptr, ulii, max_person_id_len, uiptr, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (update_sex_fname) {
        retval = update_indiv_sexes(update_sex_fname, cptr, ulii, max_person_id_len, uiptr, sex_nm, sex_male);
	if (retval) {
	  goto plink_ret_1;
	}
      }
    }
    if (keepfamname) {
      retval = include_or_exclude(keepfamname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 6);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (keepname) {
      retval = include_or_exclude(keepname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 2);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (removefamname) {
      retval = include_or_exclude(removefamname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 7);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (removename) {
      retval = include_or_exclude(removename, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 3);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (filter_attrib_indiv_fname) {
      retval = filter_attrib(filter_attrib_indiv_fname, filter_attrib_indiv_liststr, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 1);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (om_ip->marker_fname) {
      // would rather do this with pre-sorted markers, but that might break
      // order-of-operations assumptions in existing pipelines
      retval = load_oblig_missing(bedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, sex_male, chrom_info_ptr, om_ip);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_indivs_file(filtername, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, filtervals_flattened, mfilter_col);
      if (retval) {
	goto plink_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if (gender_unk_ct && popcount_longs_exclude(pheno_nm, sex_nm, unfiltered_indiv_ctl)) {
    if (!(sex_missing_pheno & ALLOW_NO_SEX)) {
      if ((!sex_missing_pheno) && (calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_WRITE_COVAR))) {
	if (calculation_type & (~(CALC_MAKE_BED | CALC_RECODE | CALC_WRITE_COVAR))) {
	  logprint("Error: When ambiguous-sex samples with phenotype data are present,\n--make-bed/--recode/--write-covar usually cannot be combined with other\ncommands.  Split them across multiple PLINK runs, or use\n--allow-no-sex/--must-have-sex.\n");
	  goto plink_ret_INVALID_CMDLINE;
	}
      } else {
	// either --must-have-sex without --allow-no-sex, or no data generation
	// command
        bitfield_and(pheno_nm, sex_nm, unfiltered_indiv_ctl);
      }
    }
  }
  if (misc_flags & MISC_PRUNE) {
    bitfield_ornot(indiv_exclude, pheno_nm, unfiltered_indiv_ctl);
    zero_trailing_bits(indiv_exclude, unfiltered_indiv_ct);
    indiv_exclude_ct = popcount_longs(indiv_exclude, unfiltered_indiv_ctl);
  }

  if (filter_binary & (FILTER_BINARY_CASES | FILTER_BINARY_CONTROLS)) {
    if (!pheno_c) {
      logprint("Error: --filter-cases/--filter-controls requires a case/control phenotype.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    ii = indiv_exclude_ct;
    // fcc == 1: exclude all zeroes in pheno_c
    // fcc == 2: exclude all ones in pheno_c
    // -> flip on fcc == 1
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, pheno_c, (filter_binary & FILTER_BINARY_CASES)? 1 : 0, pheno_nm);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      sprintf(logbuf, "Error: All %s removed due to case/control status (--filter-%s).\n", g_species_plural, (filter_binary & FILTER_BINARY_CASES)? "cases" : "controls");
      logprintb();
      goto plink_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to case/control status (--filter-%s).\n", ii, species_str(ii), (filter_binary & FILTER_BINARY_CASES)? "cases" : "controls");
    logprintb();
  }
  if (filter_binary & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_male, (filter_binary & FILTER_BINARY_MALES)? 1 : 0, sex_nm);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      sprintf(logbuf, "Error: All %s removed due to gender filter (--filter-%s).\n", g_species_plural, (filter_binary & FILTER_BINARY_MALES)? "males" : "females");
      logprintb();
      goto plink_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to gender filter (--filter-%s).\n", ii, species_str(ii), (filter_binary & FILTER_BINARY_MALES)? "males" : "females");
    logprintb();
  }
  if (filter_binary & (FILTER_BINARY_FOUNDERS | FILTER_BINARY_NONFOUNDERS)) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, founder_info, (filter_binary & FILTER_BINARY_FOUNDERS)? 1 : 0, NULL);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      sprintf(logbuf, "Error: All %s removed due to founder status (--filter-%s).\n", g_species_plural, (filter_binary & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
      logprintb();
      goto plink_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to founder status (--filter-%s).\n", ii, species_str(ii), (filter_binary & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
    logprintb();
  }

  if (mind_thresh < 1.0) {
    retval = mind_filter(bedfile, bed_offset, outname, outname_end, mind_thresh, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, sex_male, chrom_info_ptr, om_ip);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (cluster_ptr->fname) {
    retval = load_clusters(cluster_ptr->fname, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, cluster_ptr->keep_fname, cluster_ptr->keep_flattened, cluster_ptr->remove_fname, cluster_ptr->remove_flattened);
    if (retval) {
      goto plink_ret_1;
    }
  }
  // er, obvious todo: convert this to a local variable
  g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  if (!g_indiv_ct) {
    // defensive; currently shouldn't happen since we're actually checking at
    // every filter
    sprintf(logbuf, "Error: No %s pass QC.\n", g_species_plural);
    logprintb();
    goto plink_ret_ALL_SAMPLES_EXCLUDED;
  }
  if ((g_indiv_ct == 1) && (relationship_or_ibc_req(calculation_type) || distance_req(calculation_type, read_dists_fname) || (calculation_type & (CALC_GENOME | CALC_CLUSTER | CALC_NEIGHBOR)))) {
    sprintf(logbuf, "Error: More than 1 %s required for pairwise analysis.\n", g_species_singular);
    goto plink_ret_INVALID_CMDLINE_2;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  bitfield_andnot(pheno_nm, indiv_exclude, unfiltered_indiv_ctl);
  if (pheno_c) {
    bitfield_and(pheno_c, pheno_nm, unfiltered_indiv_ctl);
  }
  bitfield_andnot(founder_info, indiv_exclude, unfiltered_indiv_ctl);

  if ((parallel_tot > 1) && (parallel_tot > g_indiv_ct / 2)) {
    sprintf(logbuf, "Error: Too many --parallel jobs (maximum %" PRIuPTR "/2 = %" PRIuPTR ").\n", g_indiv_ct, g_indiv_ct / 2);
    goto plink_ret_INVALID_CMDLINE_2;
  }
  if (g_thread_ct > 1) {
    if ((calculation_type & (CALC_RELATIONSHIP | CALC_REL_CUTOFF | CALC_GDISTANCE_MASK | CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE | CALC_GENOME | CALC_REGRESS_REL | CALC_UNRELATED_HERITABILITY | CALC_LD)) || ((calculation_type & CALC_MODEL) && (model_modifier & (MODEL_PERM | MODEL_MPERM))) || ((calculation_type & CALC_GLM) && (glm_modifier & (GLM_PERM | GLM_MPERM))) || ((calculation_type & CALC_TESTMISS) && (testmiss_modifier & (TESTMISS_PERM | TESTMISS_MPERM))) || ((calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) && (!read_genome_fname) && ((cluster_ptr->ppc != 0.0) || (!read_dists_fname))) || ((calculation_type & CALC_EPI) && (epi_ip->modifier & EPI_FAST))) {
      sprintf(logbuf, "Using %u threads (change this with --threads).\n", g_thread_ct);
      logprintb();
    } else {
      logprint("Using 1 thread (no multithreaded calculations invoked).\n");
    }
  }

  if ((calculation_type & (CALC_SEXCHECK | CALC_MISSING_REPORT | CALC_GENOME | CALC_HOMOZYG)) || cluster_ptr->mds_dim_ct) {
    calc_plink_maxfid(unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, &plink_maxfid, &plink_maxiid);
  }
  plink_maxsnp = calc_plink_maxsnp(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len);

  if (indiv_sort & (INDIV_SORT_NATURAL | INDIV_SORT_ASCII)) {
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, (indiv_sort & INDIV_SORT_NATURAL)? strcmp_natural_deref : strcmp_deref);
    if (retval) {
      goto plink_ret_1;
    }
    indiv_sort_map = uiptr;
    wkspace_reset((unsigned char*)cptr);
  }

  if ((misc_flags & (MISC_MAKE_FOUNDERS | MISC_MAKE_FOUNDERS_FIRST)) == MISC_MAKE_FOUNDERS) {
    if (make_founders(unfiltered_indiv_ct, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, indiv_exclude, founder_info)) {
      goto plink_ret_NOMEM;
    }
  }

  if (calculation_type & CALC_WRITE_CLUSTER) {
    retval = write_clusters(outname, outname_end, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, (misc_flags / MISC_WRITE_CLUSTER_OMIT_UNASSIGNED) & 1, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER))))) {
      goto plink_ret_1;
    }
  }

  // this currently has to come last since covar data structures refer to
  // filtered individual indices.
  if (covar_fname) {
    // update this as more covariate-referencing commands are added
    if (!(calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_WRITE_COVAR | CALC_GXE | CALC_GLM | CALC_LASSO))) {
      logprint("Warning: Ignoring --covar since no commands reference the covariates.\n");
    } else {
      // if only --gxe, ignore --covar-name/--covar-number
      uii = (calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_WRITE_COVAR | CALC_GLM | CALC_LASSO))? 1 : 0;
      retval = load_covars(covar_fname, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, missing_phenod, uii? covar_modifier : (covar_modifier & COVAR_KEEP_PHENO_ON_MISSING_COV), uii? covar_range_list_ptr : NULL, gxe_mcovar, &covar_ct, &covar_names, &max_covar_name_len, pheno_nm, &covar_nm, &covar_d, &gxe_covar_nm, &gxe_covar_c);
      if (retval) {
	goto plink_ret_1;
      }
    }
  }

  pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_indiv_ctl);
  if (!pheno_nm_ct) {
    logprint("Note: No phenotypes present.\n");
    hwe_modifier |= HWE_THRESH_ALL;
  } else if (pheno_c) {
    pheno_ctrl_ct = popcount_longs_exclude(pheno_nm, pheno_c, unfiltered_indiv_ctl);
    if (pheno_nm_ct != g_indiv_ct) {
      sprintf(logbuf, "%u case%s, %u control%s, and %" PRIuPTR " missing phenotype%s present.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "" : "s", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "" : "s", g_indiv_ct - pheno_nm_ct, (g_indiv_ct - pheno_nm_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "%u case%s and %u control%s present.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "" : "s", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "" : "s");
    }
    logprintb();
    if (!pheno_ctrl_ct) {
      hwe_modifier |= HWE_THRESH_ALL;
    }
  } else {
    if (pheno_nm_ct != g_indiv_ct) {
      sprintf(logbuf, "%u quantitative phenotype%s present (%" PRIuPTR " missing).\n", pheno_nm_ct, (pheno_nm_ct == 1)? "" : "s", g_indiv_ct - pheno_nm_ct);
    } else {
      sprintf(logbuf, "%u quantitative phenotype%s present.\n", pheno_nm_ct, (pheno_nm_ct == 1)? "" : "s");
    }
    logprintb();
    hwe_modifier |= HWE_THRESH_ALL;
  }

  if (unfiltered_marker_ct == marker_exclude_ct) {
    // defensive
    logprint("Error: No variants remaining.\n");
    goto plink_ret_ALL_MARKERS_EXCLUDED;
  }
  retval = calc_freqs_and_hwe(bedfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, founder_info, nonfounders, (misc_flags / MISC_MAF_SUCC) & 1, set_allele_freqs, &marker_reverse, &marker_allele_cts, bed_offset, (hwe_thresh > 0.0) || (calculation_type & CALC_HARDY) || (geno_thresh < 1.0), hwe_modifier & HWE_THRESH_ALL, (pheno_nm_ct && pheno_c)? (calculation_type & CALC_HARDY) : 0, pheno_nm, pheno_nm_ct? pheno_c : NULL, &hwe_lls, &hwe_lhs, &hwe_hhs, &hwe_ll_cases, &hwe_lh_cases, &hwe_hh_cases, &hwe_ll_allfs, &hwe_lh_allfs, &hwe_hh_allfs, &hwe_hapl_allfs, &hwe_haph_allfs, &indiv_male_ct, &indiv_f_ct, &indiv_f_male_ct, wt_needed, &topsize, &marker_weights, exponent, chrom_info_ptr, sex_nm, sex_male, map_is_unsorted & UNSORTED_SPLIT_CHROM, &hh_exists);
  if (retval) {
    goto plink_ret_1;
  }

  if (freqname) {
    retval = read_external_freqs(freqname, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, chrom_info_ptr, marker_allele_ptrs, marker_allele_cts, set_allele_freqs, (misc_flags / MISC_MAF_SUCC) & 1, exponent, wt_needed, marker_weights);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (!(misc_flags & MISC_KEEP_ALLELE_ORDER)) {
    calc_marker_reverse_bin(marker_reverse, marker_exclude, unfiltered_marker_ct, unfiltered_marker_ct - marker_exclude_ct, set_allele_freqs);
  }
  if (a1alleles || a2alleles) {
    retval = load_ax_alleles(a1alleles? a1alleles : a2alleles, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_allele_ptrs, &max_marker_allele_len, marker_reverse, marker_ids, max_marker_id_len, set_allele_freqs, a2alleles? 1 : 0);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (marker_allele_ptrs) {
    swap_reversed_marker_alleles(unfiltered_marker_ct, marker_reverse, marker_allele_ptrs);
  }

  // contrary to the PLINK 1.07 flowchart, --freq effectively resolves before
  // --geno.
  if (calculation_type & CALC_FREQ) {
    if (cluster_ct && (!(misc_flags & MISC_FREQX))) {
      memcpy(outname_end, ".frq.strat", 11);
      retval = write_stratified_freqs(bedfile, bed_offset, outname, plink_maxsnp, unfiltered_marker_ct, marker_exclude, zero_extra_chroms, chrom_info_ptr, marker_ids, max_marker_id_len, marker_allele_ptrs, max_marker_allele_len, unfiltered_indiv_ct, g_indiv_ct, indiv_f_ct, founder_info, nonfounders, sex_male, indiv_f_male_ct, marker_reverse, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len);
    } else {
      if (misc_flags & MISC_FREQX) {
	memcpy(outname_end, ".frqx", 6);
      } else if (misc_flags & MISC_FREQ_COUNTS) {
	memcpy(outname_end, ".frq.count", 11);
      } else {
	memcpy(outname_end, ".frq", 5);
      }
      retval = write_freqs(outname, plink_maxsnp, unfiltered_marker_ct, marker_exclude, set_allele_freqs, zero_extra_chroms, chrom_info_ptr, marker_ids, max_marker_id_len, marker_allele_ptrs, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, hwe_hapl_allfs, hwe_haph_allfs, indiv_f_ct, indiv_f_male_ct, misc_flags, marker_reverse);
    }
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ))))) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_MISSING_REPORT) {
    retval = write_missingness_reports(bedfile, bed_offset, outname, outname_end, plink_maxfid, plink_maxiid, plink_maxsnp, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, zero_extra_chroms, chrom_info_ptr, om_ip, marker_ids, max_marker_id_len, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, pheno_nm, sex_male, indiv_male_ct, person_ids, max_person_id_len, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, hh_exists);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_MISSING_REPORT))))) {
      goto plink_ret_1;
    }
  }
  if (geno_thresh < 1.0) {
    ulii = binary_geno_filter(geno_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, g_indiv_ct, indiv_male_ct, hwe_lls, hwe_lhs, hwe_hhs, chrom_info_ptr, om_ip);
    if (marker_exclude_ct == unfiltered_marker_ct) {
      logprint("Error: All variants removed due to missing genotype data (--geno).\n");
      goto plink_ret_ALL_MARKERS_EXCLUDED;
    }
    sprintf(logbuf, "%" PRIuPTR " variant%s removed due to missing genotype data (--geno).\n", ulii, (ulii == 1)? "" : "s");
    logprintb();
  }
  oblig_missing_cleanup(om_ip);
  wkspace_reset(marker_allele_cts);
  marker_allele_cts = NULL;
  if (calculation_type & CALC_HARDY) {
    retval = hardy_report(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, hwe_lls, hwe_lhs, hwe_hhs, hwe_modifier, hwe_ll_cases, hwe_lh_cases, hwe_hh_cases, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, pheno_nm_ct, pheno_c, zero_extra_chroms, chrom_info_ptr);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_HARDY))))) {
      goto plink_ret_1;
    }
  }
  if (hwe_thresh > 0.0) {
    if (enforce_hwe_threshold(hwe_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, hwe_lls, hwe_lhs, hwe_hhs, hwe_modifier, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, chrom_info_ptr)) {
      goto plink_ret_ALL_MARKERS_EXCLUDED;
    }
  }
  if ((min_maf != 0.0) || (max_maf != 0.5)) {
    if (enforce_maf_threshold(min_maf, max_maf, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, set_allele_freqs)) {
      goto plink_ret_ALL_MARKERS_EXCLUDED;
    }
  }
  if (min_bp_space) {
    if (map_is_unsorted & (UNSORTED_BP | UNSORTED_SPLIT_CHROM)) {
      logprint("Error: --bp-space requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
      goto plink_ret_INVALID_FORMAT;
    }
    enforce_min_bp_space(min_bp_space, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, chrom_info_ptr);
  }

  if (wt_needed) {
    calc_marker_weights(exponent, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, marker_weights);
  }
  wkspace_reset(hwe_lls);

  if (sip->fname) {
    if (map_is_unsorted & (UNSORTED_BP | UNSORTED_SPLIT_CHROM)) {
      logprint("Error: --set/--make-set requires a sorted .bim file.  Retry this command after\nusing --make-bed to sort your data.\n");
      goto plink_ret_INVALID_FORMAT;
    }
    retval = define_sets(sip, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, marker_ids, max_marker_id_len, chrom_info_ptr);
    if (retval) {
      goto plink_ret_1;
    }
  }

  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (!marker_ct) {
    // defensive
    logprint("Error: All variants fail QC.\n");
    goto plink_ret_ALL_MARKERS_EXCLUDED;
  }
  sprintf(logbuf, "%" PRIuPTR " variant%s and %" PRIuPTR " %s pass filters and QC%s.\n", marker_ct, (marker_ct == 1)? "" : "s", g_indiv_ct, species_str(g_indiv_ct), (calculation_type & CALC_REL_CUTOFF)? " (before --rel-cutoff)": "");
  logprintb();

  if (wt_needed) {
    // normalize included marker weights to add to just under 2^32.  (switch to
    // 2^64 if/when 32-bit performance becomes less important than accuracy on
    // 50+ million marker sets.)
    dxx = 0.0;
    marker_uidx = 0;
    marker_idx = 0;
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      marker_idx += marker_uidx_stop - marker_uidx;
      dptr = &(marker_weights[marker_uidx]);
      dptr2 = &(marker_weights[marker_uidx_stop]);
      marker_uidx = marker_uidx_stop;
      do {
        dxx += *dptr++;
      } while (dptr < dptr2);
    } while (marker_idx < marker_ct);
    // subtract marker_ct to guard against marker_weight_sum overflow from
    // rounding
    dxx = (4294967296.0 - ((double)((intptr_t)marker_ct))) / dxx;
    if (wkspace_alloc_ui_checked(&marker_weights_i, marker_idx * sizeof(int32_t))) {
      goto plink_ret_NOMEM;
    }
    marker_uidx = 0;
    marker_idx = 0;
    uiptr = marker_weights_i;
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      marker_idx += marker_uidx_stop - marker_uidx;
      dptr = &(marker_weights[marker_uidx]);
      dptr2 = &(marker_weights[marker_uidx_stop]);
      marker_uidx = marker_uidx_stop;
      do {
        uii = (uint32_t)((*dptr++) * dxx + 0.5);
        marker_weight_sum += uii;
        *uiptr++ = uii;
      } while (dptr < dptr2);
    } while (marker_idx < marker_ct);
    wkspace_left += topsize;
    topsize = 0;
  }

  if (calculation_type & CALC_SEXCHECK) {
    retval = sexcheck(bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, sex_nm, sex_male, (misc_flags / MISC_IMPUTE_SEX) & 1, check_sex_fthresh, check_sex_mthresh, chrom_info_ptr, set_allele_freqs, &gender_unk_ct);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if ((calculation_type & CALC_GENOME) || genome_skip_write) {
    retval = populate_pedigree_rel_info(&pri, unfiltered_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_WRITE_SET) {
    retval = write_set(sip, outname, outname_end, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, marker_pos, zero_extra_chroms, chrom_info_ptr);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_WRITE_SNPLIST) {
    retval = write_snplist(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, NULL, 0);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_LIST_23_INDELS) {
    retval = write_snplist(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, 1);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE)) {
    if (gender_unk_ct && (sex_missing_pheno & MUST_HAVE_SEX)) {
      pheno_nm_datagen = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      memcpy(pheno_nm_datagen, pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
      bitfield_and(pheno_nm_datagen, sex_nm, unfiltered_indiv_ctl);
    }
    if (covar_fname) {
      retval = write_covars(outname, outname_end, write_covar_modifier, write_covar_dummy_max_categories, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm_datagen? pheno_nm_datagen : pheno_nm, pheno_c, pheno_d, missing_phenod, output_missing_pheno, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d);
      if (retval) {
	goto plink_ret_1;
      }
    }
    if (calculation_type & CALC_MAKE_BED) {
      retval = make_bed(bedfile, bed_offset, mapname, map_cols, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, marker_reverse, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm_datagen? pheno_nm_datagen : pheno_nm, pheno_c, pheno_d, output_missing_pheno, map_is_unsorted, indiv_sort_map, misc_flags, splitx_bound1, splitx_bound2, update_chr, flip_subset_fname, cluster_ptr->zerofname, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, hh_exists, chrom_info_ptr);
      if (retval) {
        goto plink_ret_1;
      }
    }
    if (calculation_type & CALC_RECODE) {
      retval = recode(recode_modifier, bedfile, bed_offset, outname, outname_end, recode_allele_name, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, marker_ids, max_marker_id_len, marker_cms, marker_allele_ptrs, max_marker_allele_len, marker_pos, marker_reverse, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm_datagen? pheno_nm_datagen : pheno_nm, pheno_c, pheno_d, output_missing_pheno, map_is_unsorted, misc_flags, hh_exists, chrom_info_ptr);
      if (retval) {
        goto plink_ret_1;
      }
    }
    if (pheno_nm_datagen) {
      free(pheno_nm_datagen);
      pheno_nm_datagen = NULL;
    }
  }

  if ((calculation_type & CALC_EPI) && epi_ip->twolocus_mkr1) {
    retval = twolocus(epi_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, chrom_info_ptr, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, pheno_nm, pheno_nm_ct, pheno_ctrl_ct, pheno_c, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_HOMOZYG) {
    if (map_is_unsorted & UNSORTED_BP) {
      logprint("Error: Run-of-homozygosity scanning requires a sorted .map/.bim.  Retry this\ncommand after using --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = calc_homozyg(homozyg_ptr, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, zero_extra_chroms, chrom_info_ptr, marker_pos, g_indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, outname, outname_end, pheno_nm, pheno_c, pheno_d, missing_pheno, sex_male);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_LD_PRUNE) {
    if (map_is_unsorted & UNSORTED_BP) {
      logprint("Error: LD-based pruning requires a sorted .map/.bim.  Retry this command after\nusing --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = ld_prune(ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_indiv_ct, founder_info, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if ((calculation_type & CALC_EPI) && epi_ip->ld_mkr1) {
    retval = twolocus(epi_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, chrom_info_ptr, unfiltered_indiv_ct, founder_info, 0, NULL, 0, 0, NULL, sex_male, NULL, NULL, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_LD) {
    if ((!(ldip->modifier & (LD_MATRIX_SHAPEMASK & LD_INTER_CHR))) && (map_is_unsorted & UNSORTED_BP)) {
      logprint("Error: Windowed --r/--r2 runs require a sorted .map/.bim.  Retry this command\nafter using --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = ld_report(threads, ldip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, plink_maxsnp, set_allele_freqs, zero_extra_chroms, chrom_info_ptr, marker_pos, unfiltered_indiv_ct, founder_info, parallel_idx, parallel_tot, sex_male, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_TESTMISHAP) {
    logprint("Error: --test-mishap is currently under development.\n.");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto plink_ret_1;
  }

  /*
  if (calculation_type & CALC_REGRESS_PCS) {
    retval = calc_regress_pcs(evecname, regress_pcs_modifier, max_pcs, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, marker_allele_ptrs, zero_extra_chroms, chrom_info_ptr, marker_pos, g_indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, sex_nm, sex_male, pheno_d, missing_phenod, outname, outname_end, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }
  */

  // sometimes no more need for marker_ids/marker_allele_ptrs, conditional
  // unload to clear space for IBS matrix, etc.?  (probably want to initially
  // load at far end of stack to make this workable...)

  if (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) {
    wkspace_mark_postcluster = wkspace_base;
    ulii = (g_indiv_ct * (g_indiv_ct - 1)) >> 1;
    if (cluster_ptr->mds_dim_ct) {
#ifndef __LP64__
      // catch 32-bit intptr_t overflow
      if (g_indiv_ct > 23169) {
        goto plink_ret_NOMEM;
      }
#endif
      if (((!read_dists_fname) && (!read_genome_fname)) || (cluster_ptr->modifier & CLUSTER_MISSING)) {
	if ((!(cluster_ptr->modifier & CLUSTER_MDS)) || (!cluster_ct)) {
          if (wkspace_alloc_d_checked(&mds_plot_dmatrix_copy, ulii * sizeof(double))) {
            goto plink_ret_NOMEM;
          }
	} else {
	  ulii = cluster_ct + g_indiv_ct - cluster_starts[cluster_ct];
          if (wkspace_alloc_d_checked(&mds_plot_dmatrix_copy, (ulii * (ulii - 1)) * (sizeof(double) / 2))) {
            goto plink_ret_NOMEM;
          }
	}
      }
    }

    if (cluster_ct) {
      ulii = cluster_ct + g_indiv_ct - cluster_starts[cluster_ct];
#ifndef __LP64__
      if (ulii > 23169) {
	goto plink_ret_NOMEM;
      }
#endif
      ulii = (ulii * (ulii - 1)) >> 1;
#ifndef __LP64__
    } else if (g_indiv_ct > 23169) {
      goto plink_ret_NOMEM;
#endif
    }
    if (wkspace_alloc_ul_checked(&cluster_merge_prevented, ((ulii + (BITCT - 1)) / BITCT) * sizeof(intptr_t))) {
      goto plink_ret_NOMEM;
    }
    if (cluster_ct || (calculation_type & CALC_GENOME) || genome_skip_write) {
      if (wkspace_alloc_d_checked(&cluster_sorted_ibs, ulii * sizeof(double))) {
	goto plink_ret_NOMEM;
      }
      if (cluster_ptr->modifier & CLUSTER_GROUP_AVG) {
        fill_double_zero(cluster_sorted_ibs, ulii);
      } else {
	for (uljj = 0; uljj < ulii; uljj++) {
	  cluster_sorted_ibs[uljj] = 1.0;
	}
      }
    }
    wkspace_mark_precluster = wkspace_base;
  }

  wkspace_mark2 = wkspace_base;

  if (relationship_or_ibc_req(calculation_type)) {
    if (rel_calc_type & REL_CALC_SINGLE_PREC) {
      retval = calc_rel_f(threads, parallel_idx, parallel_tot, calculation_type, rel_calc_type, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, ibc_type, (float)rel_cutoff, set_allele_freqs, chrom_info_ptr);
    } else {
      retval = calc_rel(threads, parallel_idx, parallel_tot, calculation_type, rel_calc_type, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, ibc_type, rel_cutoff, set_allele_freqs, &rel_ibc, chrom_info_ptr);
    }
    if (retval) {
      goto plink_ret_1;
    }

    if (calculation_type & CALC_REGRESS_REL) {
      retval = regress_rel_main(unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, regress_rel_iters, regress_rel_d, threads, pheno_d);
      if (retval) {
	goto plink_ret_1;
      }
    }
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      retval = calc_unrelated_herit(calculation_type, ibc_type, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, pheno_d, rel_ibc, (misc_flags / MISC_UNRELATED_HERITABILITY_STRICT) & 1, unrelated_herit_covg, unrelated_herit_covr, unrelated_herit_tol);
      if (retval) {
	goto plink_ret_1;
      }
    }
#endif
    wkspace_reset(g_indiv_missing_unwt);
    g_indiv_missing_unwt = NULL;
    g_missing_dbl_excluded = NULL;
  }

  /*
  if (calculation_type & CALC_REGRESS_PCS_DISTANCE) {
    logprint("Error: --regress-pcs-distance has not yet been written.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto plink_ret_1;
  } else
  */
  if (distance_req(calculation_type, read_dists_fname)) {
    retval = calc_distance(threads, parallel_idx, parallel_tot, bedfile, bed_offset, outname, outname_end, calculation_type, dist_calc_type, marker_exclude, marker_ct, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, chrom_info_ptr, wt_needed, marker_weight_sum, marker_weights_i, exponent);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (read_dists_fname && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))) {
    // use delayed and specialized load for --cluster/--neighbour, since PPC
    // values may be needed, and user may want to process a distance matrix too
    // large to be loaded in memory by doing some pre-clustering
    dists_alloc = (g_indiv_ct * (g_indiv_ct - 1)) * (sizeof(double) / 2);
    if (wkspace_alloc_d_checked(&g_dists, dists_alloc)) {
      goto plink_ret_NOMEM;
    }
    retval = read_dists(read_dists_fname, read_dists_id_fname, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, 0, NULL, NULL, 0, 0, g_dists, 0, NULL, NULL);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & CALC_IBS_TEST) {
    retval = ibs_test_calc(threads, read_dists_fname, unfiltered_indiv_ct, indiv_exclude, ibs_test_perms, pheno_nm_ct, pheno_ctrl_ct, pheno_nm, pheno_c);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_GROUPDIST) {
    retval = groupdist_calc(threads, unfiltered_indiv_ct, indiv_exclude, groupdist_iters, groupdist_d, pheno_nm_ct, pheno_ctrl_ct, pheno_nm, pheno_c);
    if (retval) {
      goto plink_ret_1;
    }
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    retval = regress_distance(calculation_type, g_dists, pheno_d, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, g_thread_ct, regress_iters, regress_d);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (read_dists_fname && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))) {
    wkspace_reset((unsigned char*)g_dists);
    g_dists = NULL;
  }

  if ((calculation_type & CALC_GENOME) || genome_skip_write) {
    wkspace_reset(wkspace_mark2);
    g_dists = NULL;
    retval = calc_genome(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_pos, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info, parallel_idx, parallel_tot, outname, outname_end, nonfounders, calculation_type, genome_modifier, ppc_gap, genome_min_pi_hat, genome_max_pi_hat, pheno_nm, pheno_c, pri, genome_skip_write);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) {
    retval = calc_cluster_neighbor(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, read_dists_fname, read_dists_id_fname, read_genome_fname, outname, outname_end, calculation_type, cluster_ct, cluster_map, cluster_starts, cluster_ptr, missing_pheno, neighbor_n1, neighbor_n2, ppc_gap, pheno_c, mds_plot_dmatrix_copy, cluster_merge_prevented, cluster_sorted_ibs, wkspace_mark_precluster, wkspace_mark_postcluster);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if ((calculation_type & CALC_EPI) && (epi_ip->modifier & (EPI_FAST | EPI_REG))) {
    if ((map_is_unsorted & UNSORTED_BP) && (epi_ip->modifier & EPI_FAST_CASE_ONLY)) {
      logprint("Error: --fast-epistasis case-only requires a sorted .map/.bim.  Retry this\ncommand after using --make-bed to sort your data.\n");
      goto plink_ret_INVALID_CMDLINE;
    }
    retval = epistasis_report(threads, epi_ip, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, marker_pos, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, pheno_nm, pheno_nm_ct, pheno_ctrl_ct, pheno_c, pheno_d, parallel_idx, parallel_tot, outname, outname_end, sip);
    if (retval) {
      goto plink_ret_1;
    }
  }

  if (calculation_type & (CALC_MODEL | CALC_GXE | CALC_GLM | CALC_LASSO | CALC_CMH | CALC_HOMOG | CALC_TESTMISS)) {
    // can't use pheno_ctrl_ct in here since new phenotypes may be loaded, and
    // we don't bother updating it...
    if ((!pheno_all) && (!loop_assoc_fname)) {
      outname_end2 = outname_end;
      goto plink_skip_all_pheno;
    }
    uii = 0; // phenotype/cluster number
    *outname_end = '.';
    if (loop_assoc_fname) {
      retval = load_clusters(loop_assoc_fname, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, NULL, NULL, NULL, NULL);
      if (retval) {
	goto plink_ret_1;
      }
      if (pheno_d) {
	free(pheno_d);
	pheno_d = NULL;
      }
      if (!pheno_c) {
	pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      }
    } else {
      wkspace_mark = wkspace_base;
      retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
      if (retval) {
	goto plink_ret_1;
      }
    }
    do {
      if (loop_assoc_fname) {
	if (uii == cluster_ct) {
	  break;
	}
	outname_end2 = strcpya(&(outname_end[1]), &(cluster_ids[uii * max_cluster_id_len]));
	fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
	ukk = cluster_starts[uii + 1];
	for (ujj = cluster_starts[uii]; ujj < ukk; ujj++) {
	  SET_BIT(pheno_c, cluster_map[ujj]);
	}
	uii++;
      } else {
	// --all-pheno
	if (pheno_modifier & PHENO_MERGE) {
	  memcpy(pheno_nm, orig_pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
	  if (orig_pheno_c) {
	    if (!pheno_c) {
	      free(pheno_d);
	      pheno_d = NULL;
	      pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
	    }
	    memcpy(pheno_c, orig_pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t));
	  } else {
	    memcpy(pheno_d, orig_pheno_d, unfiltered_indiv_ct * sizeof(double));
	  }
	} else {
	  fill_ulong_zero(pheno_nm, unfiltered_indiv_ctl);
	  if (pheno_c) {
	    free(pheno_c);
	    pheno_c = NULL;
	  }
	  if (pheno_d) {
	    free(pheno_d);
	    pheno_d = NULL;
	  }
	}
	uii++;
      plink_skip_empty_pheno:
	rewind(phenofile);
	outname_end[1] = '\0';
	retval = load_pheno(phenofile, unfiltered_indiv_ct, indiv_exclude_ct, cptr, max_person_id_len, uiptr, missing_pheno, missing_pheno_len, (misc_flags / MISC_AFFECTION_01) & 1, uii, NULL, pheno_nm, &pheno_c, &pheno_d, &(outname_end[1]), (uintptr_t)((&(outname[FNAMESIZE - 32])) - outname_end));
	if (retval == LOAD_PHENO_LAST_COL) {
	  wkspace_reset(wkspace_mark);
	  break;
	} else if (retval) {
	  goto plink_ret_1;
	}
	bitfield_andnot(pheno_nm, indiv_exclude, unfiltered_indiv_ctl);
	if (gender_unk_ct && (!(sex_missing_pheno & ALLOW_NO_SEX))) {
	  bitfield_and(pheno_nm, sex_nm, unfiltered_indiv_ctl);
	}
	pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_indiv_ctl);
	if (!pheno_nm_ct) {
	  goto plink_skip_empty_pheno;
	}
	if (!outname_end[1]) {
	  outname_end[1] = 'P';
	  outname_end2 = uint32_write(&(outname_end[2]), uii);
	} else {
          outname_end2 = (char*)memchr(&(outname_end[1]), '\0', FNAMESIZE);
	}
      }
      *outname_end2 = '\0';
    plink_skip_all_pheno:
      if (calculation_type & CALC_MODEL) {
	if (pheno_d) {
	  retval = qassoc(threads, bedfile, bed_offset, outname, outname_end2, model_modifier, model_mperm_val, pfilter, mtest_adjust, adjust_lambda, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_allele_ptrs, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, cluster_ct, cluster_map, cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_d, sex_male, hh_exists, perm_batch_size, sip);
	} else {
	  retval = model_assoc(threads, bedfile, bed_offset, outname, outname_end2, model_modifier, model_cell_ct, model_mperm_val, ci_size, ci_zt, pfilter, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_len, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, cluster_ct, cluster_map, loop_assoc_fname? NULL : cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, sex_male, sip);
	}
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (calculation_type & CALC_GLM) {
	if (!(glm_modifier & GLM_NO_SNP)) {
          retval = glm_assoc(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_len, marker_reverse, zero_extra_chroms, condition_mname, condition_fname, chrom_info_ptr, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, cluster_ct, cluster_map, cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, perm_batch_size, sip);
	} else {
	  retval = glm_assoc_nosnp(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, cluster_ct, cluster_map, cluster_starts, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, perm_batch_size, sip);
	}
	if (retval) {
	  goto plink_ret_1;
	}
      }
      // if case/control phenotype loaded with --all-pheno, skip --gxe
      if ((calculation_type & CALC_GXE) && pheno_d) {
	retval = gxe_assoc(bedfile, bed_offset, outname, outname_end, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, pheno_nm, pheno_d, gxe_covar_nm, gxe_covar_c, sex_male, hh_exists);
	if (retval) {
	  goto plink_ret_1;
	}
      }
      if (calculation_type & CALC_LASSO) {
	retval = lasso(threads, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, pheno_nm_ct, lasso_h2, lasso_minlambda, lasso_select_covars_range_list_ptr, misc_flags, pheno_nm, pheno_c, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_male, hh_exists);
        if (retval) {
	  goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_CMH) && pheno_c) {
        retval = cmh_assoc();
        if (retval) {
          goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_HOMOG) && pheno_c) {
	retval = homog_assoc();
        if (retval) {
          goto plink_ret_1;
	}
      }
      if ((calculation_type & CALC_TESTMISS) && pheno_c) {
        retval = testmiss(threads, bedfile, bed_offset, outname, outname_end, testmiss_mperm_val, testmiss_modifier, pfilter, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, cluster_ct, cluster_map, loop_assoc_fname? NULL : cluster_starts, apip, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, sex_male, hh_exists);
        if (retval) {
	  goto plink_ret_1;
	}
      }
    } while (pheno_all || loop_assoc_fname);
  }
  if (calculation_type & CALC_CLUMP) {
    retval = clump_reports(bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, clump_ip, sex_male, hh_exists);
    if (retval) {
      goto plink_ret_1;
    }
  }
  while (0) {
  plink_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  plink_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  plink_ret_INVALID_FORMAT_2:
    logprintb();
  plink_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  plink_ret_INVALID_CMDLINE_2:
    logprintb();
  plink_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  plink_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  plink_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  plink_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 plink_ret_1:
  if (topsize) {
    wkspace_left += topsize;
  }
  free_cond(pheno_nm_datagen);
  free_cond(orig_pheno_d);
  free_cond(orig_pheno_c);
  free_cond(orig_pheno_nm);
  free_cond(pheno_d);
  free_cond(pheno_c);
  free_cond(id_list);
  free_cond(pid_list);
  fclose_cond(phenofile);
  fclose_cond(famfile);
  fclose_cond(bedfile);
  if (marker_allele_ptrs && (max_marker_allele_len > 2)) {
    ulii = unfiltered_marker_ct * 2;
    for (marker_uidx = 0; marker_uidx < ulii; marker_uidx++) {
      cptr = marker_allele_ptrs[marker_uidx];
      if ((cptr < g_one_char_strs) || (cptr >= &(g_one_char_strs[512]))) {
	free(cptr);
      }
    }
  }
  return retval;
}

// output-missing-phenotype + terminating null
#define MAX_FLAG_LEN 25

static inline int32_t is_flag(char* param) {
  char cc = param[1];
  return ((*param == '-') && ((cc > '9') || (cc < '0'))); 
}

static inline char* is_flag_start(char* param) {
  char cc = param[1];
  if ((*param == '-') && ((cc > '9') || (cc < '0'))) {
    return (cc == '-')? (&(param[2])) : (&(param[1]));
  }
  return NULL;
}

uint32_t param_count(int32_t argc, char** argv, int32_t flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any parameter not beginning with "--" as optional.
  int32_t opt_params = 0;
  int32_t cur_idx = flag_idx + 1;
  while (cur_idx < argc) {
    if (is_flag(argv[cur_idx])) {
      break;
    }
    opt_params++;
    cur_idx++;
  }
  return opt_params;
}

int32_t enforce_param_ct_range(uint32_t param_ct, char* flag_name, uint32_t min_ct, uint32_t max_ct) {
  if (param_ct > max_ct) {
    if (max_ct > min_ct) {
      sprintf(logbuf, "Error: %s accepts at most %d parameter%s.%s", flag_name, max_ct, (max_ct == 1)? "" : "s", errstr_append);
    } else {
      sprintf(logbuf, "Error: %s only accepts %d parameter%s.%s", flag_name, max_ct, (max_ct == 1)? "" : "s", errstr_append);
    }
    return -1;
  } else if (param_ct < min_ct) {
    if (min_ct == 1) {
      sprintf(logbuf, "Error: Missing %s parameter.%s", flag_name, errstr_append);
    } else {
      sprintf(logbuf, "Error: %s requires %s%d parameters.%s", flag_name, (min_ct < max_ct)? "at least " : "", min_ct, errstr_append);
    }
    return -1;
  }
  return 0;
}

int32_t parse_next_range(uint32_t param_ct, char range_delim, char** argv, uint32_t* cur_param_idx_ptr, char** cur_arg_pptr, char** range_start_ptr, uint32_t* rs_len_ptr, char** range_end_ptr, uint32_t* re_len_ptr) {
  // Starts reading from argv[cur_param_idx][cur_pos].  If a valid range is
  // next, range_start + rs_len + range_end + re_len are updated.  If only a
  // single item is next, range_end is set to NULL and range_start + rs_len are
  // updated.  If there are no items left, range_start is set to NULL.  If
  // the input is not well-formed, -1 is returned instead of 0.
  uint32_t cur_param_idx = *cur_param_idx_ptr;
  char* cur_arg_ptr = *cur_arg_pptr;
  char cc;
  if (cur_param_idx > param_ct) {
    *cur_arg_pptr = NULL;
    return 0;
  }
  while (1) {
    cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
	*range_start_ptr = NULL;
	return 0;
      }
      cur_arg_ptr = argv[cur_param_idx];
      cc = *cur_arg_ptr;
    }
    if (cc == range_delim) {
      return -1;
    }
    if (cc != ',') {
      break;
    }
    cur_arg_ptr++;
  }
  *range_start_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if ((!cc) || (cc == ',')) {
      *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
      *cur_arg_pptr = cur_arg_ptr;
      *range_end_ptr = NULL;
      return 0;
    }
  } while (cc != range_delim);
  *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
  cc = *(++cur_arg_ptr);
  if ((!cc) || (cc == ',') || (cc == range_delim)) {
    return -1;
  }
  *range_end_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if (cc == range_delim) {
      return -1;
    }
  } while (cc && (cc != ','));
  *re_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_end_ptr));
  *cur_arg_pptr = cur_arg_ptr;
  return 0;
}

int32_t parse_chrom_ranges(uint32_t param_ct, char range_delim, char** argv, uintptr_t* chrom_mask, Chrom_info* chrom_info_ptr, uint32_t allow_extra_chroms, char* cur_flag_str) {
  uint32_t argct = 0;
  uint32_t cur_param_idx = 1;
  int32_t retval = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  int32_t chrom_code_start;
  int32_t chrom_code_end;
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", cur_flag_str, argv[cur_param_idx], errstr_append);
	goto parse_chrom_ranges_ret_INVALID_CMDLINE;
      }
      if (!range_start) {
	break;
      }
      chrom_code_start = get_chrom_code2(chrom_info_ptr, range_start, rs_len);
      if (chrom_code_start == -1) {
	range_start[rs_len] = '\0';
	if (!allow_extra_chroms) {
	  sprintf(logbuf, "Error: Invalid --%s chromosome code '%s'.%s", cur_flag_str, range_start, errstr_append);
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE;
	} else if (range_end) {
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE_2;
	}
        if (push_ll_str(&(chrom_info_ptr->incl_excl_name_stack), range_start)) {
	  goto parse_chrom_ranges_ret_NOMEM;
	}
      } else if (range_end) {
        chrom_code_end = get_chrom_code2(chrom_info_ptr, range_end, re_len);
	if (chrom_code_end == -1) {
	  if (!allow_extra_chroms) {
	    range_end[re_len] = '\0';
	    sprintf(logbuf, "Error: Invalid --%s chromosome code '%s'.%s", cur_flag_str, range_end, errstr_append);
	    goto parse_chrom_ranges_ret_INVALID_CMDLINE;
	  } else {
	    goto parse_chrom_ranges_ret_INVALID_CMDLINE_2;
	  }
	}
        if (chrom_code_end <= chrom_code_start) {
	  range_start[rs_len] = '\0';
	  range_end[re_len] = '\0';
	  sprintf(logbuf, "Error: --%s chromosome code '%s' is not greater than '%s'.%s", cur_flag_str, range_end, range_start, errstr_append);
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE;
	}
	fill_bits(chrom_mask, chrom_code_start, chrom_code_end + 1 - chrom_code_start);
      } else {
        SET_BIT(chrom_mask, chrom_code_start);
      }
      argct++;
    }
  }
  if (!argct) {
    sprintf(logbuf, "Error: --%s requires at least one value.%s", cur_flag_str, errstr_append);
    logprintb();
    return -1;
  }
  while (0) {
  parse_chrom_ranges_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  parse_chrom_ranges_ret_INVALID_CMDLINE_2:
    logprint("Error: Chromosome ranges cannot include nonstandard names.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  parse_chrom_ranges_ret_INVALID_CMDLINE:
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

int32_t parse_name_ranges(uint32_t param_ct, char range_delim, char** argv, Range_list* range_list_ptr, uint32_t require_posint) {
  uint32_t name_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t name_max_len = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  char* cur_name_str;
  char* dup_check;
  unsigned char* cur_name_starts_range;
  int32_t last_val;
  int32_t cur_val;
  // two passes.  first pass: count parameters, determine name_max_len;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(logbuf, "Error: Invalid %s parameter '%s'.%s", argv[0], argv[cur_param_idx], errstr_append);
        logprintb();
        return RET_INVALID_CMDLINE;
      }
      if (!range_start) {
	break;
      }
      name_ct++;
      if (rs_len > name_max_len) {
	name_max_len = rs_len; // does NOT include trailing null yet
      }
      if (range_end) {
	name_ct++;
	if (re_len > name_max_len) {
	  name_max_len = re_len;
	}
      }
    }
  }
  if (!name_ct) {
    sprintf(logbuf, "Error: %s requires at least one value.%s", argv[0], errstr_append);
    logprintb();
    return RET_INVALID_CMDLINE;
  }
  range_list_ptr->name_max_len = ++name_max_len;
  range_list_ptr->name_ct = name_ct;
  range_list_ptr->names = (char*)malloc(name_ct * ((uintptr_t)name_max_len));
  if (!range_list_ptr->names) {
    return RET_NOMEM;
  }
  range_list_ptr->starts_range = (unsigned char*)malloc(name_ct * sizeof(char));
  if (!range_list_ptr->starts_range) {
    return RET_NOMEM;
  }
  cur_name_str = range_list_ptr->names;
  cur_name_starts_range = range_list_ptr->starts_range;
  cur_param_idx = 1;
  cur_arg_ptr = argv[1];
  while (1) {
    parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      if (require_posint) {
	last_val = 0;
	for (cur_param_idx = 0; cur_param_idx < name_ct; cur_param_idx++) {
	  cur_name_str = &(range_list_ptr->names[cur_param_idx * ((uintptr_t)name_max_len)]);
	  dup_check = cur_name_str; // actually a numeric check
	  do {
	    if (is_not_digit(*dup_check)) {
	      sprintf(logbuf, "Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	      logprintb();
	      return RET_INVALID_CMDLINE;
	    }
	  } while (*(++dup_check));
	  cur_val = atoi(cur_name_str);
	  if (cur_val < 1) {
	    sprintf(logbuf, "Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	    logprintb();
	    return RET_INVALID_CMDLINE;
	  }
	  if (range_list_ptr->starts_range[cur_param_idx]) {
	    last_val = cur_val;
	  } else {
	    if (cur_val <= last_val) {
	      sprintf(logbuf, "Error: Invalid %s range '%s-%s'.\n", argv[0], &(range_list_ptr->names[(cur_param_idx - 1) * name_max_len]), cur_name_str);
	      logprintb();
	      return RET_INVALID_CMDLINE;
	    }
	    last_val = 0;
	  }
	}
      }
      return 0;
    }
    memcpyx(cur_name_str, range_start, rs_len, 0);
    dup_check = range_list_ptr->names;
    while (dup_check < cur_name_str) {
      if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	sprintf(logbuf, "Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	logprintb();
	return RET_INVALID_CMDLINE;
      }
      dup_check = &(dup_check[name_max_len]);
    }
    cur_name_str = &(cur_name_str[name_max_len]);
    if (range_end) {
      *cur_name_starts_range++ = 1;
      memcpyx(cur_name_str, range_end, re_len, 0);
      dup_check = range_list_ptr->names;
      while (dup_check < cur_name_str) {
	if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	  sprintf(logbuf, "Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	  logprintb();
	  return RET_INVALID_CMDLINE;
	}
        dup_check = &(dup_check[name_max_len]);
      }
      cur_name_str = &(cur_name_str[name_max_len]);
      *cur_name_starts_range++ = 0;
    } else {
      *cur_name_starts_range++ = 0;
    }
  }
}

void invalid_arg(char* argv) {
  sprintf(logbuf, "Error: Unrecognized flag ('%s').%s%s", argv, (argv[0] == '-')? "" : "  All flags must be preceded by 1-2 dashes.", errstr_append);
}

void print_ver() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

char extract_char_param(char* ss) {
  // maps c, 'c', and "c" to c, and anything else to the null char.  This is
  // intended to support e.g. always using '#' to designate a # parameter
  // without worrying about differences between shells.
  char cc = ss[0];
  if (((cc == '\'') || (cc == '"')) && (ss[1]) && (ss[2] == cc) && (!ss[3])) {
    return ss[1];
  } else if (cc && (!ss[1])) {
    return cc;
  } else {
    return '\0';
  }
}

int32_t alloc_string(char** sbuf, const char* source) {
  uint32_t slen = strlen(source) + 1;
  *sbuf = (char*)malloc(slen * sizeof(char));
  if (!(*sbuf)) {
    return -1;
  }
  memcpy(*sbuf, source, slen);
  return 0;
}

int32_t alloc_fname(char** fnbuf, char* source, char* argptr, uint32_t extra_size) {
  uint32_t slen = strlen(source) + 1;
  if (slen > (FNAMESIZE - extra_size)) {
    sprintf(logbuf, "Error: --%s filename too long.\n", argptr);
    logprintb();
    return RET_OPEN_FAIL;
  }
  *fnbuf = (char*)malloc((slen + extra_size) * sizeof(char));
  if (!(*fnbuf)) {
    return RET_NOMEM;
  }
  memcpy(*fnbuf, source, slen);
  return 0;
}

int32_t alloc_and_flatten(char** flattened_buf_ptr, char** sources, uint32_t param_ct) {
  uint32_t totlen = 1;
  char* bufptr;
  uint32_t param_idx;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    totlen += 1 + strlen(sources[param_idx]);
  }
  bufptr = (char*)malloc(totlen);
  if (!bufptr) {
    return RET_NOMEM;
  }
  *flattened_buf_ptr = bufptr;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    bufptr = strcpyax(bufptr, sources[param_idx], '\0');
  }
  *bufptr = '\0';
  return 0;
}

int32_t alloc_2col(Two_col_params** tcbuf, char** params_ptr, char* argptr, uint32_t param_ct) {
  uint32_t slen = strlen(*params_ptr) + 1;
  int32_t ii;
  char cc;
  if (slen > FNAMESIZE) {
    sprintf(logbuf, "Error: --%s filename too long.\n", argptr);
    logprintb();
    return RET_OPEN_FAIL;
  }
  *tcbuf = (Two_col_params*)malloc(sizeof(Two_col_params) + slen);
  if (!(*tcbuf)) {
    return RET_NOMEM;
  }
  memcpy((*tcbuf)->fname, params_ptr[0], slen);
  (*tcbuf)->skip = 0;
  (*tcbuf)->skipchar = '\0';
  if (param_ct > 1) {
    ii = atoi(params_ptr[1]);
    if (ii < 1) {
      sprintf(logbuf, "Error: Invalid --%s column number.\n", argptr);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    (*tcbuf)->colx = ii;
    if (param_ct > 2) {
      ii = atoi(params_ptr[2]);
      if (ii < 1) {
	sprintf(logbuf, "Error: Invalid --%s variant ID column number.\n", argptr);
	logprintb();
	return RET_INVALID_FORMAT;
      }
      (*tcbuf)->colid = ii;
      if (param_ct == 4) {
	cc = params_ptr[3][0];
	if ((cc < '0') || (cc > '9')) {
	  cc = extract_char_param(params_ptr[3]);
	  if (!cc) {
            goto alloc_2col_invalid_skip;
	  }
	  (*tcbuf)->skipchar = cc;
	} else {
	  if (atoiz(params_ptr[3], &ii)) {
	  alloc_2col_invalid_skip:
	    sprintf(logbuf, "Error: Invalid --%s skip parameter.  This needs to either be a\nsingle character (usually '#') which, when present at the start of a line,\nindicates it should be skipped; or the number of initial lines to skip.  (Note\nthat in shells such as bash, '#' is a special character that must be\nsurrounded by single- or double-quotes to be parsed correctly.)\n", argptr);
	    logprintb();
	    return RET_INVALID_FORMAT;
	  }
	  (*tcbuf)->skip = ii;
	}
      }
    } else {
      (*tcbuf)->colid = 1;
    }
    if ((*tcbuf)->colx == (*tcbuf)->colid) {
      sprintf(logbuf, "Error: Column numbers for --%s cannot be equal.%s", argptr, errstr_append);
      logprintb();
      return RET_INVALID_FORMAT;
    }
  } else {
    (*tcbuf)->colx = 2;
    (*tcbuf)->colid = 1;
  }
  return 0;
}

int32_t flag_match(const char* to_match, uint32_t* cur_flag_ptr, uint32_t flag_ct, char* flag_buf) {
  int32_t ii;
  while (*cur_flag_ptr < flag_ct) {
    ii = strcmp(to_match, &(flag_buf[(*cur_flag_ptr) * MAX_FLAG_LEN]));
    if (ii < 0) {
      return 0;
    }
    *cur_flag_ptr += 1;
    if (!ii) {
      flag_buf[((*cur_flag_ptr) - 1) * MAX_FLAG_LEN] = '\0';
      return 1;
    }
  }
  return 0;
}

uint32_t species_flag(uint32_t* species_code_ptr, uint32_t new_code) {
  if (*species_code_ptr) {
    logprint("Error: Multiple chromosome set flags.\n");
    return 1;
  }
  *species_code_ptr = new_code;
  return 0;
}

// these need global scope to stay around on all systems
const char species_singular_constants[][7] = {"person", "cow", "dog", "horse", "mouse", "plant", "sheep", "sample"};
const char species_plural_constants[][8] = {"people", "cows", "dogs", "horses", "mice", "plants", "sheep", "samples"};

int32_t init_delim_and_species(uint32_t flag_ct, char* flag_buf, uint32_t* flag_map, int32_t argc, char** argv, char* range_delim_ptr, Chrom_info* chrom_info_ptr) {
  // human: 22, X, Y, XY, MT
  // cow: 29, X, Y
  // dog: 38, X, Y, XY
  // horse: 31, X, Y
  // mouse: 19, X, Y
  // rice: 12
  // sheep: 26, X, Y
  const int32_t species_x_code[] = {23, 30, 39, 32, 20, -1, 27};
  const int32_t species_y_code[] = {24, 31, 40, 33, 21, -1, 28};
  const int32_t species_xy_code[] = {25, -1, 41, -1, -1, -1, -1};
  const int32_t species_mt_code[] = {26, -1, -1, -1, -1, -1, -1};
  const uint32_t species_max_code[] = {26, 31, 41, 33, 21, 12, 28};
  uint32_t species_code = SPECIES_HUMAN;
  uint32_t flag_idx = 0;
  uint32_t retval = 0;
  int32_t cur_arg;
  uint32_t param_ct;
  int32_t ii;
  uint32_t param_idx;
  fill_ulong_zero(chrom_info_ptr->haploid_mask, CHROM_MASK_WORDS);
  fill_ulong_zero(chrom_info_ptr->chrom_mask, CHROM_MASK_WORDS);
  if (flag_match("autosome-num", &flag_idx, flag_ct, flag_buf)) {
    species_code = SPECIES_UNKNOWN;
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    ii = atoi(argv[cur_arg + 1]);
    if ((ii < 1) || (ii > 59)) {
      sprintf(logbuf, "Error: Invalid --autosome-num parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    chrom_info_ptr->x_code = ii + 1;
    chrom_info_ptr->y_code = -1;
    chrom_info_ptr->xy_code = -1;
    chrom_info_ptr->mt_code = -1;
    chrom_info_ptr->max_code = ii + 1;
    chrom_info_ptr->autosome_ct = ii;
    set_bit(chrom_info_ptr->haploid_mask, ii + 1);
  }
  if (flag_match("chr-set", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_UNKNOWN)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 5)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    ii = atoi(argv[cur_arg + 1]);
    if ((!ii) || (ii > 59) || (ii < -59)) {
      sprintf(logbuf, "Error: Invalid --chr-set parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    if (ii < 0) {
      if (param_ct > 1) {
	sprintf(logbuf, "Error: --chr-set does not accept multiple parameters in haploid mode.%s", errstr_append);
	goto init_delim_and_species_ret_INVALID_CMDLINE_2;
      }
      ii = -ii;
      chrom_info_ptr->autosome_ct = ii;
      chrom_info_ptr->x_code = -1;
      chrom_info_ptr->y_code = -1;
      chrom_info_ptr->xy_code = -1;
      chrom_info_ptr->mt_code = -1;
      chrom_info_ptr->max_code = ii;
      fill_all_bits(chrom_info_ptr->haploid_mask, ((uint32_t)ii) + 1);
    } else {
      chrom_info_ptr->autosome_ct = ii;
      chrom_info_ptr->x_code = ii + 1;
      chrom_info_ptr->y_code = ii + 2;
      chrom_info_ptr->xy_code = ii + 3;
      chrom_info_ptr->mt_code = ii + 4;
      set_bit(chrom_info_ptr->haploid_mask, ii + 1);
      set_bit(chrom_info_ptr->haploid_mask, ii + 2);
      for (param_idx = 2; param_idx <= param_ct; param_idx++) {
	if (!strcmp(argv[cur_arg + param_idx], "no-x")) {
	  chrom_info_ptr->x_code = -1;
	  clear_bit(chrom_info_ptr->haploid_mask, ii + 1);
	} else if (!strcmp(argv[cur_arg + param_idx], "no-y")) {
	  chrom_info_ptr->y_code = -1;
	  clear_bit(chrom_info_ptr->haploid_mask, ii + 2);
	} else if (!strcmp(argv[cur_arg + param_idx], "no-xy")) {
	  chrom_info_ptr->xy_code = -1;
	} else if (!strcmp(argv[cur_arg + param_idx], "no-mt")) {
	  chrom_info_ptr->mt_code = -1;
	} else {
	  sprintf(logbuf, "Error: Invalid --chr-set parameter '%s'.%s", argv[cur_arg + param_idx], errstr_append);
	  goto init_delim_and_species_ret_INVALID_CMDLINE_2;
	}
      }
      if (chrom_info_ptr->mt_code != -1) {
	chrom_info_ptr->max_code = ii + 4;
      } else if (chrom_info_ptr->xy_code != -1) {
	chrom_info_ptr->max_code = ii + 3;
      } else if (chrom_info_ptr->y_code != -1) {
	chrom_info_ptr->max_code = ii + 2;
      } else if (chrom_info_ptr->x_code != -1) {
	chrom_info_ptr->max_code = ii + 1;
      } else {
	chrom_info_ptr->max_code = ii;
      }
    }
  }
  if (flag_match("cow", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_COW)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --cow doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("d", &flag_idx, flag_ct, flag_buf)) {
    // moved here to support --covar-name + --d
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    *range_delim_ptr = extract_char_param(argv[cur_arg + 1]);
    if (!(*range_delim_ptr)) {
      sprintf(logbuf, "Error: --d parameter too long (must be a single character).%s", errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    } else if ((*range_delim_ptr == '-') || (*range_delim_ptr == ',')) {
      sprintf(logbuf, "Error: --d parameter cannot be '-' or ','.%s", errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
  }
  if (flag_match("dog", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_DOG)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --dog doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("horse", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_HORSE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --horse doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("mouse", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_MOUSE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --mouse doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("rice", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_RICE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --rice doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("sheep", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_SHEEP)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --sheep doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  chrom_info_ptr->species = species_code;
  chrom_info_ptr->is_include_stack = 0;
  if (species_code != SPECIES_UNKNOWN) {
    chrom_info_ptr->x_code = species_x_code[species_code];
    chrom_info_ptr->y_code = species_y_code[species_code];
    chrom_info_ptr->xy_code = species_xy_code[species_code];
    chrom_info_ptr->mt_code = species_mt_code[species_code];
    chrom_info_ptr->max_code = species_max_code[species_code];
  }
  g_species_singular = species_singular_constants[species_code];
  g_species_plural = species_plural_constants[species_code];
  switch (species_code) {
  case SPECIES_HUMAN:
    chrom_info_ptr->autosome_ct = 22;
    chrom_info_ptr->haploid_mask[0] = 0x1800000;
    break;
  case SPECIES_COW:
    chrom_info_ptr->autosome_ct = 29;
    chrom_info_ptr->haploid_mask[0] = 0xc0000000LU;
    break;
  case SPECIES_DOG:
    chrom_info_ptr->autosome_ct = 38;
#ifdef __LP64__
    chrom_info_ptr->haploid_mask[0] = 0x18000000000LLU;
#else
    chrom_info_ptr->haploid_mask[1] = 0x180;
#endif
    break;
  case SPECIES_HORSE:
    chrom_info_ptr->autosome_ct = 31;
#ifdef __LP64__
    chrom_info_ptr->haploid_mask[0] = 0x300000000LLU;
#else
    chrom_info_ptr->haploid_mask[1] = 3;
#endif
    break;
  case SPECIES_MOUSE:
    chrom_info_ptr->autosome_ct = 19;
    chrom_info_ptr->haploid_mask[0] = 0x300000;
    break;
  case SPECIES_RICE:
    chrom_info_ptr->autosome_ct = 12;
    chrom_info_ptr->haploid_mask[0] = 0x1fff;
    break;
  case SPECIES_SHEEP:
    chrom_info_ptr->autosome_ct = 26;
    chrom_info_ptr->haploid_mask[0] = 0x18000000;
    break;
  }
  while (0) {
  init_delim_and_species_ret_INVALID_CMDLINE_2:
    logprintb();
  init_delim_and_species_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

void fill_chrom_mask(Chrom_info* chrom_info_ptr) {
  if (chrom_info_ptr->species != SPECIES_UNKNOWN) {
    fill_all_bits(chrom_info_ptr->chrom_mask, chrom_info_ptr->max_code + 1);
  } else {
    fill_all_bits(chrom_info_ptr->chrom_mask, chrom_info_ptr->autosome_ct + 1);
    // --chr-set support
    if (chrom_info_ptr->x_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->x_code);
    }
    if (chrom_info_ptr->y_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->y_code);
    }
    if (chrom_info_ptr->xy_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->xy_code);
    }
    if (chrom_info_ptr->mt_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->mt_code);
    }
  }
}

int32_t recode_type_set(uint32_t* recode_modifier_ptr, uint32_t cur_code) {
  if (*recode_modifier_ptr & (RECODE_TYPEMASK - cur_code)) {
    sprintf(logbuf, "Error: Conflicting --recode modifiers.%s", errstr_append);
    return -1;
  }
  *recode_modifier_ptr |= cur_code;
  return 0;
}

int32_t main(int32_t argc, char** argv) {
  char* outname_end = NULL;
  char** subst_argv = NULL;
  char* script_buf = NULL;
  char* rerun_buf = NULL;
  char* flag_buf = NULL;
  uint32_t* flag_map = NULL;
  char* makepheno_str = NULL;
  char* phenoname_str = NULL;
  Two_col_params* a1alleles = NULL;
  Two_col_params* a2alleles = NULL;
  char* filtervals_flattened = NULL;
  char* evecname = NULL;
  char* filtername = NULL;
  char* read_dists_fname = NULL;
  char* read_dists_id_fname = NULL;
  char* freqname = NULL;
  char* extractname = NULL;
  char* excludename = NULL;
  char* keepname = NULL;
  char* removename = NULL;
  char* keepfamname = NULL;
  char* removefamname = NULL;
  char* cm_map_fname = NULL;
  char* cm_map_chrname = NULL;
  char* phenoname = NULL;
  char* recode_allele_name = NULL;
  char* lgen_reference_fname = NULL;
  char* covar_fname = NULL;
  char* update_alleles_fname = NULL;
  Two_col_params* update_chr = NULL;
  Two_col_params* update_cm = NULL;
  Two_col_params* update_map = NULL;
  Two_col_params* update_name = NULL;
  char* update_ids_fname = NULL;
  char* update_parents_fname = NULL;
  char* update_sex_fname = NULL;
  char* loop_assoc_fname = NULL;
  char* flip_fname = NULL;
  char* flip_subset_fname = NULL;
  char* read_genome_fname = NULL;
  char* condition_mname = NULL;
  char* condition_fname = NULL;
  char* filter_attrib_fname = NULL;
  char* filter_attrib_liststr = NULL;
  char* filter_attrib_indiv_fname = NULL;
  char* filter_attrib_indiv_liststr = NULL;
  char* const_fid = NULL;
  char* vcf_filter_exceptions_flattened = NULL;
  double vcf_min_qual = -INFINITY;
  char id_delim = '\0';
  int32_t retval = 0;
  uint32_t load_params = 0; // describes what file parameters have been provided
  uint32_t load_rare = 0;
  uint32_t fam_cols = FAM_COL_13456;
  uint32_t mpheno_col = 0;
  uint32_t mwithin_col = 0;
  uint64_t misc_flags = 0;
  double thin_keep_prob = 1.0;
  uint32_t min_bp_space = 0;
  double check_sex_fthresh = 0.2;
  double check_sex_mthresh = 0.8;
  double exponent = 0.0;
  double min_maf = 0.0;
  double max_maf = 0.5;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  double rel_cutoff = 0.025;
  uint32_t cur_arg = 1;
  uint64_t calculation_type = 0;
  uint32_t rel_calc_type = 0;
  uint32_t dist_calc_type = 0;
  uint32_t mfilter_col = 0;
  uint32_t pheno_modifier = 0;
  int32_t missing_pheno = -9;
  uintptr_t groupdist_iters = ITERS_DEFAULT;
  uint32_t groupdist_d = 0;
  uintptr_t regress_iters = ITERS_DEFAULT;
  uint32_t regress_d = 0;
  uintptr_t regress_rel_iters = ITERS_DEFAULT;
  uint32_t regress_rel_d = 0;
  double unrelated_herit_tol = 0.0000001;
  double unrelated_herit_covg = 0.45;
  double unrelated_herit_covr = 0.55;
  int32_t ibc_type = 0; // -1 for cov
  uint32_t parallel_idx = 0;
  uint32_t parallel_tot = 1;
  uint32_t splitx_bound1 = 0;
  uint32_t splitx_bound2 = 0;
  uint32_t sex_missing_pheno = 0;
  uint32_t hwe_modifier = 0;
  uint32_t write_covar_modifier = 0;
  uint32_t write_covar_dummy_max_categories = 49;
  uint32_t model_modifier = 0;
  int32_t model_cell_ct = -1;
  uint32_t gxe_mcovar = 0;
  uint32_t glm_modifier = 0;
  double glm_vif_thresh = 50.0;
  uint32_t glm_xchr_model = 1;
  uint32_t ppc_gap = DEFAULT_PPC_GAP;
  uint32_t* rseeds = NULL;
  uint32_t rseed_ct = 0;
  uint32_t genome_modifier = 0;
  double genome_min_pi_hat = -1.0;
  double genome_max_pi_hat = 1.0;
  FILE* scriptfile = NULL;
  uint32_t filter_binary = 0;
  uint32_t regress_pcs_modifier = 0;
  uint32_t max_pcs = MAX_PCS_DEFAULT;
  uint32_t recode_modifier = 0;
  uint32_t allelexxxx = 0;
  uint32_t merge_type = 0;
  uint32_t indiv_sort = 0;
  uint32_t cur_flag = 0;
  uint32_t flag_ct = 0;
  uint32_t dummy_marker_ct = 0;
  uint32_t dummy_indiv_ct = 0;
  uint32_t dummy_flags = 0;
  double dummy_missing_geno = 0.0;
  double dummy_missing_pheno = 0.0;
  char* simulate_fname = NULL;
  uint32_t simulate_flags = 0;
  uint32_t simulate_cases = 1000;
  uint32_t simulate_controls = 1000;
  double simulate_prevalence = 0.01;
  char* simulate_label = NULL;
  double simulate_missing = 0.0;
  uint32_t simulate_qt_indivs = 1000;
  char* markername_from = NULL;
  char* markername_to = NULL;
  char* markername_snp = NULL;
  uint32_t snp_window_size = 0;
  int32_t marker_pos_start = -1;
  int32_t marker_pos_end = -1;
  uint32_t lgen_modifier = 0;
  uint32_t covar_modifier = 0;
  uint32_t update_map_modifier = 0;
  uint32_t model_mperm_val = 0;
  uint32_t glm_mperm_val = 0;
  uint32_t mperm_save = 0;
  uint32_t mperm_val = 0;
  double ci_size = 0.0;
  double pfilter = 1.0;
  uint32_t perm_batch_size = 0;
  uint32_t mtest_adjust = 0;
  double adjust_lambda = 0.0;
  uint32_t ibs_test_perms = DEFAULT_IBS_TEST_PERMS;
  uint32_t neighbor_n1 = 0;
  uint32_t neighbor_n2 = 0;
  uint32_t cnv_calc_type = 0;
  uint32_t cnv_indiv_mperms = 0;
  uint32_t cnv_test_mperms = 0;
  uint32_t cnv_test_region_mperms = 0;
  uint32_t cnv_enrichment_test_mperms = 0;
  uint32_t cnv_min_seglen = 0;
  uint32_t cnv_max_seglen = 0xffffffffU;
  double cnv_min_score = -INFINITY;
  double cnv_max_score = INFINITY;
  uint32_t cnv_min_sites = 0;
  uint32_t cnv_max_sites = 0xffffffffU;
  uint32_t cnv_intersect_filter_type = 0;
  char* cnv_intersect_filter_fname = NULL;
  char* cnv_subset_fname = NULL;
  uint32_t cnv_overlap_type = 0;
  double cnv_overlap_val = 0.0;
  uint32_t cnv_freq_type = 0;
  uint32_t cnv_freq_val = 0;
  double cnv_freq_val2 = 0.0;
  uint32_t cnv_test_window = 0;
  uint32_t segment_modifier = 0;
  uint32_t matrix_flag_state = 0; // 1 = present and unclaimed, 2 = claimed
  double hard_call_threshold = 0.1;
  double tail_bottom = 0.0;
  double tail_top = 0.0;
  double lasso_h2 = 0.0;
  double lasso_minlambda = -INFINITY;
  uint32_t testmiss_modifier = 0;
  uint32_t testmiss_mperm_val = 0;
  char* segment_spanning_fname = NULL;
  char* missing_code = NULL;
  char range_delim = '-';
  uint32_t modifier_23 = 0;
  double pheno_23 = INFINITY;
  char* fid_23 = NULL;
  char* iid_23 = NULL;
  char* paternal_id_23 = NULL;
  char* maternal_id_23 = NULL;
  Ll_str* file_delete_list = NULL;
  uint32_t chrom_flag_present = 0;
  uintptr_t chrom_exclude[CHROM_MASK_INITIAL_WORDS];
  // er, except for first four, these should not be preallocated...
  char outname[FNAMESIZE];
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char mergename1[FNAMESIZE];
  char mergename2[FNAMESIZE];
  char mergename3[FNAMESIZE];
  char output_missing_pheno[32];
#ifdef __APPLE__
  int32_t mib[2];
  size_t sztmp;
#endif
  unsigned char* wkspace_ua;
  char** subst_argv2;
  uint32_t param_ct;
  time_t rawtime;
  char* argptr;
  char* sptr;
  char* bubble;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t num_params;
  int32_t in_param;
  Chrom_info chrom_info;
  Oblig_missing_info oblig_missing_info;
  Aperm_info aperm;
  Cluster_info cluster;
  Set_info set_info;
  Homozyg_info homozyg;
  Ld_info ld_info;
  Epi_info epi_info;
  Clump_info clump_info;
  Range_list snps_range_list;
  Range_list covar_range_list;
  Range_list lasso_select_covars_range_list;
  Range_list parameters_range_list;
  Range_list tests_range_list;
  char* argptr2;
  char* flagptr;
  double dxx;
  char cc;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  intptr_t default_alloc_mb;
  int64_t llxx;
  Ll_str* ll_str_ptr;
#if _WIN32
  SYSTEM_INFO sysinfo;
  MEMORYSTATUSEX memstatus;
  DWORD windows_dw; // why the f*** does uint32_t not work?
#endif
  oblig_missing_init(&oblig_missing_info);
  aperm_init(&aperm);
  cluster_init(&cluster);
  set_init(&set_info);
  homozyg_init(&homozyg);
  ld_epi_init(&ld_info, &epi_info, &clump_info);
  range_list_init(&snps_range_list);
  range_list_init(&covar_range_list);
  range_list_init(&lasso_select_covars_range_list);
  range_list_init(&parameters_range_list);
  range_list_init(&tests_range_list);

  chrom_info.name_ct = 0;
  chrom_info.incl_excl_name_stack = NULL;
  for (uii = 1; uii < (uint32_t)argc; uii++) {
    if ((!strcmp("-script", argv[uii])) || (!strcmp("--script", argv[uii]))) {
      ujj = param_count(argc, argv, uii);
      if (enforce_param_ct_range(ujj, argv[uii], 1, 1)) {
	print_ver();
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      for (ujj = uii + 2; ujj < (uint32_t)argc; ujj++) {
	if ((!strcmp("-script", argv[ujj])) || (!strcmp("--script", argv[ujj]))) {
	  print_ver();
	  printf("Error: Multiple --script flags.  Merge the files into one.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      // logging not yet active, so don't use fopen_checked()
      scriptfile = fopen(argv[uii + 1], "rb");
      if (!scriptfile) {
	print_ver();
	printf(errstr_fopen, argv[uii + 1]);
	goto main_ret_OPEN_FAIL;
      }
      if (fseeko(scriptfile, 0, SEEK_END)) {
	print_ver();
	goto main_ret_READ_FAIL;
      }
      llxx = ftello(scriptfile);
      if (llxx == -1) {
	print_ver();
	goto main_ret_READ_FAIL;
      } else if (llxx > 0x7fffffff) {
	// could actually happen if user enters parameters in the wrong order,
	// so may as well catch it and print a somewhat informative error msg
	print_ver();
        fputs("Error: --script file too large.\n", stdout);
        goto main_ret_NOMEM;
      }
      rewind(scriptfile);
      ujj = (uint32_t)((uint64_t)llxx);
      script_buf = (char*)malloc(ujj);
      if (!script_buf) {
	print_ver();
	goto main_ret_NOMEM;
      }
      ukk = fread(script_buf, 1, ujj, scriptfile);
      if (ukk < ujj) {
	print_ver();
	goto main_ret_READ_FAIL;
      }
      fclose_null(&scriptfile);
      num_params = 0;
      in_param = 0;
      for (ukk = 0; ukk < ujj; ukk++) {
	if (is_space_or_eoln(script_buf[ukk])) {
	  in_param = 0;
	} else if (!in_param) {
	  num_params++;
	  in_param = 1;
	}
      }
      subst_argv = (char**)malloc((num_params + argc - 3) * sizeof(char*));
      num_params = 0;
      in_param = 0;
      for (ukk = 1; ukk < uii; ukk++) {
        subst_argv[num_params++] = argv[ukk];
      }
      for (ukk = 0; ukk < ujj; ukk++) {
	if (is_space_or_eoln(script_buf[ukk])) {
	  if (in_param) {
	    script_buf[ukk] = '\0';
	    in_param = 0;
	  }
	} else if (!in_param) {
	  subst_argv[num_params++] = &(script_buf[ukk]);
	  in_param = 1;
	}
      }
      for (ujj = uii + 2; ujj < (uint32_t)argc; ujj++) {
	subst_argv[num_params++] = argv[ujj];
      }
      argc = num_params;
      cur_arg = 0;
      argv = subst_argv;
    }
  }
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    if ((!strcmp("-rerun", argv[uii])) || (!strcmp("--rerun", argv[uii]))) {
      ujj = param_count(argc, argv, uii);
      if (enforce_param_ct_range(ujj, argv[uii], 0, 1)) {
	print_ver();
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      for (ukk = uii + ujj + 1; ukk < (uint32_t)argc; ukk++) {
	if ((!strcmp("-rerun", argv[ukk])) || (!strcmp("--rerun", argv[ukk]))) {
	  print_ver();
	  fputs("Error: Duplicate --rerun flag.\n", stdout);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      if (ujj) {
	scriptfile = fopen(argv[uii + 1], "r");
      } else {
        scriptfile = fopen(PROG_NAME_STR ".log", "r");
      }
      if (!scriptfile) {
	print_ver();
	goto main_ret_OPEN_FAIL;
      }
      if (!fgets(tbuf, MAXLINELEN, scriptfile)) {
	print_ver();
	fputs("Error: Empty log file for --rerun.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      if (!fgets(tbuf, MAXLINELEN, scriptfile)) {
	print_ver();
	fputs("Error: Only one line in --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      fclose_null(&scriptfile);
      kk = atoi(tbuf);
      if ((kk < 1) || (kk > MAXLINELEN)) {
	print_ver();
	fputs("Error: Improperly formatted --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      ukk = strlen(tbuf) + 1;
      if (ukk == MAXLINELEN) {
	print_ver();
	fputs("Error: Second line too long in --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      rerun_buf = (char*)malloc(ukk);
      memcpy(rerun_buf, tbuf, ukk);

      memset(tbuf, 1, (uint32_t)kk);
      sptr = next_item_mult(rerun_buf, 2);
      umm = 0;
      ukk = 0;
      do {
	if (no_more_items_kns(sptr)) {
	  print_ver();
	  fputs("Error: Improperly formatted --rerun log file.\n", stdout);
	  goto main_ret_INVALID_FORMAT;
	}
	argptr = is_flag_start(sptr);
	if (argptr) {
          for (unn = cur_arg; unn < (uint32_t)argc; unn++) {
	    argptr2 = is_flag_start(argv[unn]);
	    if (argptr2) {
	      if (!strcmp(argptr, argptr2)) {
		unn = 0xffffffffU;
		break;
	      }
	    }
	  }
          if (unn == 0xffffffffU) {
	    // matching flag, override --rerun
            do {
	      ukk++;
	      tbuf[umm++] = 0;
	      if (umm == (uint32_t)kk) {
		break;
	      }
	      sptr = next_item(sptr);
	    } while (!is_flag(sptr));
	  } else {
	    umm++;
	    sptr = next_item(sptr);
	  }
	} else {
	  umm++;
          sptr = next_item(sptr);
	}
      } while (umm < (uint32_t)kk);
      subst_argv2 = (char**)malloc((argc + kk - ukk - ujj - 1 - cur_arg) * sizeof(char*));
      if (!subst_argv2) {
	print_ver();
	goto main_ret_NOMEM;
      }
      unn = 0;
      for (umm = cur_arg; umm < uii; umm++) {
	subst_argv2[unn++] = argv[umm];
      }
      sptr = next_item_mult(rerun_buf, 2);
      for (umm = 0; umm < (uint32_t)kk; umm++) {
        if (tbuf[umm]) {
	  ukk = strlen_se(sptr);
	  subst_argv2[unn++] = sptr;
	  sptr[ukk] = '\0';
	  if (umm != ((uint32_t)kk) - 1) {
	    sptr = skip_initial_spaces(&(sptr[ukk + 1]));
	  }
	} else {
	  sptr = next_item(sptr);
	}
      }
      for (umm = uii + ujj + 1; umm < (uint32_t)argc; umm++) {
	subst_argv2[unn++] = argv[umm];
      }
      cur_arg = 0;
      argc = unn;
      if (subst_argv) {
	free(subst_argv);
      }
      subst_argv = subst_argv2;
      argv = subst_argv2;
      subst_argv2 = NULL;
    }
  }
  if ((cur_arg < (uint32_t)argc) && (!is_flag(argv[cur_arg]))) {
    print_ver();
    printf("Error: First parameter must be a flag.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE;
  }
  flag_ct = 0;
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    argptr = is_flag_start(argv[uii]);
    if (argptr) {
      if (!strcmp("help", argptr)) {
	print_ver();
	if ((cur_arg != 1) || (uii != 1) || subst_argv) {
	  fputs("--help present, ignoring other flags.\n", stdout);
	}
	retval = disp_help(argc - uii - 1, &(argv[uii + 1]));
	goto main_ret_1;
      }
      if ((!strcmp("h", argptr)) || (!strcmp("?", argptr))) {
	// these just act like the no-parameter case
	print_ver();
	if ((cur_arg != 1) || (uii != 1) || subst_argv) {
	  printf("-%s present, ignoring other flags.\n", argptr);
	}
	fputs(cmdline_format_str, stdout);
	fputs(notestr_null_calc2, stdout);
        retval = RET_HELP;
	goto main_ret_1;
      }
      if (strlen(argptr) >= MAX_FLAG_LEN) {
	print_ver();
	invalid_arg(argv[uii]);
	fputs(logbuf, stdout);
        goto main_ret_INVALID_CMDLINE;
      }
      flag_ct++;
    }
  }
  if (!flag_ct) {
    print_ver();
    fputs(cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
    goto main_ret_1;
  }
  flag_buf = (char*)malloc(flag_ct * MAX_FLAG_LEN * sizeof(char));
  flag_map = (uint32_t*)malloc(flag_ct * sizeof(int32_t));
  if ((!flag_buf) || (!flag_map)) {
    print_ver();
    goto main_ret_NOMEM;
  }
  flagptr = flag_buf;
  umm = 0; // parameter count increase due to aliases
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    argptr = is_flag_start(argv[uii]);
    if (argptr) {
      ukk = strlen(argptr) + 1;
      // handle aliases now, so sorting will have the desired effects
      switch (*argptr) {
      case 'Z':
	if (!strcmp(argptr, "Z-genome")) {
	  memcpy(flagptr, "genome gz", 10);
	  umm++;
	  break;
	}
	goto main_flag_copy;
      case 'a':
	if ((ukk == 11) && (!memcmp(argptr, "allele", 6))) {
	  if (match_upper(&(argptr[6]), "ACGT")) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	} else if ((ukk == 12) && (!memcmp(argptr, "allele-", 7))) {
          if (!memcmp(&(argptr[7]), "1234", 4)) {
	    memcpy(flagptr, "allele1234", 11);
	    break;
	  } else if (match_upper(&(argptr[7]), "ACGT")) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	}
	goto main_flag_copy;
      case 'b':
        if (!strcmp(argptr, "border")) {
          memcpy(flagptr, "make-set-border", 16);
	  break;
	}
	goto main_flag_copy;
      case 'c':
        if (!strcmp(argptr, "chr-excl")) {
          logprint("Note: --chr-excl flag has been renamed to --not-chr.\n");
	  memcpy(flagptr, "not-chr", 8);
	  break;
	} else if (!strcmp(argptr, "cmh")) {
	  memcpy(flagptr, "mh", 3);
	  break;
	}
	goto main_flag_copy;
      case 'e':
	if (!strcmp(argptr, "extract-snp")) {
	  memcpy(flagptr, "snp", 4);
	  break;
	} else if (!strcmp(argptr, "exponent")) {
	  logprint("Note: --exponent flag has been renamed to --distance-exp.\n");
	  memcpy(flagptr, "distance-exp", 13);
	  break;
	}
	goto main_flag_copy;
      case 'f':
	if (!strcmp(argptr, "frqx")) {
	  memcpy(flagptr, "freqx", 6);
	  break;
	}
	goto main_flag_copy;
      case 'g':
	if (!strcmp(argptr, "grm-cutoff")) {
          memcpy(flagptr, "rel-cutoff", 11);
	  break;
	}
	goto main_flag_copy;
      case 'k':
	if (!memcmp(argptr, "k", 2)) {
	  memcpy(flagptr, "K", 2);
	  break;
	}
	goto main_flag_copy;
      case 'l':
	if (!strcmp(argptr, "list")) {
	  memcpy(flagptr, "recode list", 12);
	  printf("Note: --list flag deprecated.  Use '--recode list' instead.\n");
	  recode_modifier |= RECODE_LIST;
	  misc_flags |= MISC_SET_HH_MISSING;
	  break;
	} else if (!strcmp(argptr, "load-dists")) {
          memcpy(flagptr, "read-dists", 11);
          printf("Note: --load-dists flag has been renamed to --read-dists.\n");
          break;
	}
	goto main_flag_copy;
      case 'm':
	if (!strcmp(argptr, "missing_code")) {
	  memcpy(flagptr, "missing-code", 13);
	  break;
	} else if (!strcmp(argptr, "mh1")) {
	  memcpy(flagptr, "mh", 3);
	  break;
	} else if (!strcmp(argptr, "make-set-collapse-all")) {
	  memcpy(flagptr, "set-collapse-all", 17);
	  break;
	}
	goto main_flag_copy;
      case 'n':
	if (!strcmp(argptr, "neighbor")) {
	  memcpy(flagptr, "neighbour", 10);
	  break;
	} else if (!strcmp(argptr, "num_threads")) {
	  memcpy(flagptr, "threads", 8);
	  break;
	}
	goto main_flag_copy;
      case 'r':
	if ((ukk >= 8) && (!memcmp(argptr, "recode", 6))) {
	  ujj = 0; // alias match?
	  argptr2 = &(argptr[6]);
          switch (ukk) {
	  case 8:
            if (tolower(*argptr2) == 'a') {
	      memcpy(flagptr, "recode A", 9);
	      recode_modifier |= RECODE_A;
	      ujj = 1;
	    }
	    break;
	  case 9:
	    if (!memcmp(argptr2, "12", 2)) {
              memcpy(flagptr, "recode 12", 10);
	      recode_modifier |= RECODE_12;
	      ujj = 1;
            } else if (match_upper(argptr2, "AD")) {
              memcpy(flagptr, "recode AD", 10);
	      recode_modifier |= RECODE_AD;
	      ujj = 1;
	    } else if (match_upper(argptr2, "HV")) {
	      memcpy(flagptr, "recode HV-1chr", 15);
	      recode_modifier |= RECODE_HV_1CHR;
              printf("Note: --recodeHV flag deprecated.  Use '--recode HV' or '--recode HV-1chr'.\n");
	      ujj = 2;
	    }
	    break;
	  case 11:
	    if (!memcmp(argptr2, "-vcf", 4)) {
	      memcpy(flagptr, "recode vcf", 11);
	      recode_modifier |= RECODE_VCF;
	      ujj = 1;
	    }
	    break;
          case 12:
            if (!memcmp(argptr2, "-lgen", 5)) {
              memcpy(flagptr, "recode lgen", 12);
	      recode_modifier |= RECODE_LGEN;
	      // backwards compatibility
	      misc_flags |= MISC_SET_HH_MISSING;
              ujj = 1;
	    }
	    break;
	  case 13:
	    if (!memcmp(argptr2, "-rlist", 6)) {
	      memcpy(flagptr, "recode rlist", 13);
	      recode_modifier |= RECODE_RLIST;
	      misc_flags |= MISC_SET_HH_MISSING;
	      ujj = 1;
	    }
	    break;
	  case 14:
	    if (!memcmp(argptr2, "-beagle", 7)) {
	      memcpy(flagptr, "recode beagle", 14);
	      recode_modifier |= RECODE_BEAGLE;
	      ujj = 1;
	    } else if (!memcmp(argptr2, "-bimbam", 7)) {
	      memcpy(flagptr, "recode bimbam-1chr", 19);
	      recode_modifier |= RECODE_BIMBAM_1CHR;
	      misc_flags |= MISC_SET_HH_MISSING;
	      printf("Note: --recode-bimbam flag deprecated.  Use '--recode bimbam' or\n'--recode bimbam-1chr'.\n");
	      ujj = 2;
	    }
	    break;
	  case 17:
	    if (!memcmp(argptr2, "-fastphase", 10)) {
	      memcpy(flagptr, "recode fastphase-1chr", 22);
	      recode_modifier |= RECODE_FASTPHASE_1CHR;
	      misc_flags |= MISC_SET_HH_MISSING;
	      printf("Note: --recode-fastphase flag deprecated.  Use '--recode fastphase' or\n'--recode fastphase-1chr'.\n");
	      ujj = 2;
	    } else if (!memcmp(argptr2, "-structure", 10)) {
	      memcpy(flagptr, "recode structure", 17);
	      recode_modifier |= RECODE_STRUCTURE;
	      misc_flags |= MISC_SET_HH_MISSING;
	      ujj = 1;
	    }
	    break;
	  }
	  if (ujj) {
	    if (ujj == 1) {
	      printf("Note: --%s flag deprecated.  Use '%s ...'.\n", argptr, flagptr);
	    }
	    umm++;
	    break;
	  }
	} else if (!strcmp(argptr, "reference-allele")) {
	  memcpy(flagptr, "a1-allele", 10);
	  break;
	}
	goto main_flag_copy;
      case 't':
        if (!strcmp(argptr, "thread-num")) {
	  memcpy(flagptr, "threads", 8);
	  break;
	}
	goto main_flag_copy;
      case 'u':
	if (!strcmp(argptr, "update-freq")) {
	  memcpy(flagptr, "read-freq", 10);
	  break;
	} else if (!strcmp(argptr, "update-ref-allele")) {
	  // GCTA alias
	  memcpy(flagptr, "a1-allele", 10);
	  break;
	}
	goto main_flag_copy;
      case 'v':
	if (!strcmp(argptr, "version")) {
	  fputs(ver_str, stdout);
	  putchar('\n');
          goto main_ret_1;
	}
	// fall through
      default:
      main_flag_copy:
	memcpy(flagptr, argptr, ukk);
      }
      flagptr = &(flagptr[MAX_FLAG_LEN]);
      flag_map[cur_flag++] = uii;
    }
  }
  sptr = (char*)malloc(flag_ct * MAX_FLAG_LEN);
  if (!sptr) {
    print_ver();
    goto main_ret_NOMEM;
  }
  qsort_ext2(flag_buf, flag_ct, MAX_FLAG_LEN, strcmp_deref, (char*)flag_map, sizeof(int32_t), sptr, MAX_FLAG_LEN);
  free(sptr);
  ujj = strlen_se(flag_buf);
  for (cur_flag = 1; cur_flag < flag_ct; cur_flag++) {
    ukk = strlen_se(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    if ((ujj == ukk) && (!memcmp(&(flag_buf[(cur_flag - 1) * MAX_FLAG_LEN]), &(flag_buf[cur_flag * MAX_FLAG_LEN]), ukk))) {
      flag_buf[cur_flag * MAX_FLAG_LEN + ukk] = '\0'; // just in case of aliases
      print_ver();
      printf("Error: Duplicate --%s flag.\n", &(flag_buf[cur_flag * MAX_FLAG_LEN]));
      goto main_ret_INVALID_CMDLINE;
    }
    ujj = ukk;
  }

  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    if ((!memcmp("silent", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 7)) || (!memcmp("gplink", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 7))) {
      freopen("/dev/null", "w", stdout);
      break;
    }
  }
  print_ver();
  uii = 5;
  memcpy(outname, PROG_NAME_STR, 6);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    ii = memcmp("out", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 4);
    if (!ii) {
      ujj = flag_map[cur_flag];
      ukk = param_count(argc, argv, ujj);
      if (enforce_param_ct_range(ukk, argv[ujj], 1, 1)) {
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      if (strlen(argv[ujj + 1]) > (FNAMESIZE - MAX_POST_EXT)) {
	fputs("Error: --out parameter too long.\n", stdout);
	goto main_ret_OPEN_FAIL;
      }
      uii = strlen(argv[ujj + 1]);
      memcpy(outname, argv[ujj + 1], uii + 1);
      outname_end = &(outname[uii]);
    }
    if (ii <= 0) {
      break;
    }
  }
  memcpy(&(outname[uii]), ".log", 5);
  logfile = fopen(outname, "w");
  if (!logfile) {
    printf("Error: Failed to open %s.  Try ", outname);
    if (!memcmp(outname, PROG_NAME_STR, 6)) {
      printf("using --out.\n");
    } else {
      printf("changing the --out parameter.\n");
    }
    goto main_ret_OPEN_FAIL;
  }
  printf("Logging to %s.\n", outname);
  outname[uii] = '\0';

  logstr(ver_str);
  sprintf(logbuf, "\n%d argument%s:", argc + umm - cur_arg, (argc + umm - cur_arg == 1)? "" : "s");
  logstr(logbuf);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    logstr(" --");
    logstr(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    ii = flag_map[cur_flag] + 1;
    while ((ii < argc) && (!is_flag(argv[ii]))) {
      logstr(" ");
      logstr(argv[ii++]);
    }
  }
#if _WIN32
  windows_dw = 4 * MAXLINELEN + 256;
  if (GetComputerName(tbuf, &windows_dw))
#else
  if (gethostname(tbuf, 4 * MAXLINELEN + 256) != -1)
#endif
  {
    logstr("\nHostname: ");
    logstr(tbuf);
  }
  logstr("\nWorking directory: ");
  getcwd(tbuf, FNAMESIZE);
  logstr(tbuf);
  logstr("\nStart time: ");
  time(&rawtime);
  logstr(ctime(&rawtime));
  logstr("\n");

#if _WIN32
  GetSystemInfo(&sysinfo);
  g_thread_ct = sysinfo.dwNumberOfProcessors;
#else
  ii = sysconf(_SC_NPROCESSORS_ONLN);
  if (ii == -1) {
    g_thread_ct = 1;
  } else {
    g_thread_ct = ii;
  }
#endif
  if (g_thread_ct > 8) {
    if (g_thread_ct > MAX_THREADS) {
      g_thread_ct = MAX_THREADS;
    } else {
      g_thread_ct--;
    }
  }
  memcpy(mapname, PROG_NAME_STR ".map", 10);
  memcpy(pedname, PROG_NAME_STR ".ped", 10);
  famname[0] = '\0';
  memcpyl3(output_missing_pheno, "-9");
  // stuff that must be processed before regular alphabetical loop
  retval = init_delim_and_species(flag_ct, flag_buf, flag_map, argc, argv, &range_delim, &chrom_info);
  if (retval) {
    goto main_ret_1;
  }
  fill_ulong_zero(chrom_exclude, CHROM_MASK_INITIAL_WORDS);
  cur_flag = 0;
  do {
    argptr = &(flag_buf[cur_flag * MAX_FLAG_LEN]);
    if (!(*argptr)) {
      // preprocessed
      continue;
    }
    argptr2 = &(argptr[1]);
    cur_arg = flag_map[cur_flag];
    param_ct = param_count(argc, argv, cur_arg);
    switch (*argptr) {
    case '1':
      if (*argptr2 == '\0') {
	misc_flags |= MISC_AFFECTION_01;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case '2':
      if (!memcmp(argptr2, "3file", 6)) {
	if (chrom_info.species != SPECIES_HUMAN) {
	  logprint("Error: --23file cannot be used with a nonhuman species flag.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 7)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = strlen(argv[cur_arg + 1]);
	if (ii > FNAMESIZE - 1) {
	  logprint("Error: --23file filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
        memcpy(pedname, argv[cur_arg + 1], ii + 1);
	if (param_ct > 1) {
	  if (strchr(argv[cur_arg + 2], ' ')) {
	    logprint("Error: Space present in --23file family ID.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (alloc_string(&fid_23, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	  if (param_ct > 2) {
	    if (strchr(argv[cur_arg + 3], ' ')) {
	      logprint("Error: Space present in --23file individual ID.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (alloc_string(&iid_23, argv[cur_arg + 3])) {
	      goto main_ret_NOMEM;
	    }
	    if (param_ct > 3) {
	      cc = extract_char_param(argv[cur_arg + 4]);
	      if ((cc == 'M') || (cc == 'm') || (cc == '1')) {
		modifier_23 |= M23_MALE;
	      } else if ((cc == 'F') || (cc == 'f') || (cc == '2')) {
		modifier_23 |= M23_FEMALE;
	      } else if (cc == '0') {
		modifier_23 |= M23_FORCE_MISSING_SEX;
	      } else if ((cc != 'I') && (cc != 'i')) {
		logprint("Error: Invalid --23file sex parameter (M or 1 = male, F or 2 = female,\nI = infer from data, 0 = force missing).\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (param_ct > 4) {
		if (scan_double(argv[cur_arg + 5], &pheno_23)) {
		  sprintf(logbuf, "Error: Invalid --23file phenotype '%s'.%s", argv[cur_arg + 5], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
		if (param_ct > 5) {
		  if (strchr(argv[cur_arg + 6], ' ')) {
		    logprint("Error: Space present in --23file paternal ID.\n");
		    goto main_ret_INVALID_CMDLINE;
		  }
		  if (alloc_string(&paternal_id_23, argv[cur_arg + 6])) {
		    goto main_ret_NOMEM;
		  }
		  if (param_ct > 6) {
		    if (strchr(argv[cur_arg + 7], ' ')) {
		      logprint("Error: Space present in --23file maternal ID.\n");
		      goto main_ret_INVALID_CMDLINE;
		    }
		    if (alloc_string(&maternal_id_23, argv[cur_arg + 7])) {
		      goto main_ret_NOMEM;
		    }
		  }
		}
	      }
	    }
	  }
	}
	load_rare = LOAD_RARE_23;
      } else if ((!memcmp(argptr2, "3file-convert-xy", 17)) || (!memcmp(argptr2, "3file-make-xylist", 18))) {
        sprintf(logbuf, "Error: --%s has been retired due to brain-damaged design.  Use\n--split-x instead.%s", argptr, errstr_append);
        goto main_ret_INVALID_CMDLINE_3;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'D':
      if (*argptr2 == '\0') {
	logprint("Note: --D flag deprecated.  Use e.g. '--r2 dprime'.\n");
	ld_info.modifier |= LD_DPRIME;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'K':
      if (*argptr2 == '\0') {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
          sprintf(logbuf, "Error: Invalid --K cluster count '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.min_ct = ii;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'a':
      if (!memcmp(argptr2, "utosome", 8)) {
	fill_bits(chrom_info.chrom_mask, 1, chrom_info.autosome_ct);
	chrom_info.is_include_stack = 1;
	chrom_flag_present = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "utosome-xy", 11)) {
	if (chrom_flag_present) {
          logprint("Error: --autosome-xy cannot be used with --autosome.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (chrom_info.xy_code == -1) {
	  sprintf(logbuf, "Error: --autosome-xy used with a species lacking an XY region.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	fill_bits(chrom_info.chrom_mask, 1, chrom_info.autosome_ct);
	set_bit(chrom_info.chrom_mask, chrom_info.xy_code);
	chrom_info.is_include_stack = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llow-extra-chr", 15)) {
	if (load_rare == LOAD_RARE_23) {
	  logprint("Error: --allow-extra-chr cannot currently be used with --23file.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  if (memcmp("0", argv[cur_arg + 1], 2)) {
            sprintf(logbuf, "Error: Invalid --allow-extra-chr parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_ZERO_EXTRA_CHROMS;
	}
	misc_flags |= MISC_ALLOW_EXTRA_CHROMS;
      } else if (!memcmp(argptr2, "llow-no-sex", 12)) {
        sex_missing_pheno |= ALLOW_NO_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ll", 3)) {
	logprint("Note: --all flag has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llele1234", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  if (strcmp("multichar", argv[cur_arg + 1])) {
	    sprintf(logbuf, "Error: Invalid --allele1234 parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  allelexxxx = ALLELE_RECODE_MULTICHAR;
	} else {
	  allelexxxx = ALLELE_RECODE;
	}
      } else if (!memcmp(argptr2, "lleleACGT", 9)) {
	if (allelexxxx) {
	  logprint("Error: --allele1234 and --alleleACGT cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  if (strcmp("multichar", argv[cur_arg + 1])) {
	    sprintf(logbuf, "Error: Invalid --alleleACGT parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  allelexxxx = ALLELE_RECODE_ACGT | ALLELE_RECODE_MULTICHAR;
	} else {
	  allelexxxx = ALLELE_RECODE_ACGT;
	}
      } else if (!memcmp(argptr2, "llele-count", 12)) {
	lgen_modifier |= LGEN_ALLELE_COUNT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ll-pheno", 9)) {
	pheno_modifier |= PHENO_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ssoc", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 7)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "counts")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with 'counts'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_ASSOC_COUNTS;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with 'fisher'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher-midp")) {
            if (model_modifier & MODEL_QMASK) {
              sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with 'fisher-midp'.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
            model_modifier |= MODEL_FISHER | MODEL_FISHER_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: --assoc 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with 'genedrop'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "p2")) {
	    model_modifier |= MODEL_ASSOC_P2;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (model_modifier & MODEL_PERM) {
	      sprintf(logbuf, "Error: --assoc 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --assoc 'mperm' modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    kk = atoi(&(argv[cur_arg + uii][6]));
	    if (kk < 1) {
	      sprintf(logbuf, "Error: Invalid --assoc mperm parameter '%s'.%s", &(argv[cur_arg + uii][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_mperm_val = (uint32_t)kk;
	    model_modifier |= MODEL_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "qt-means")) {
	    if (model_modifier & MODEL_DMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' does not make sense with a case/control-specific\nmodifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_QT_MEANS;
	  } else if (!strcmp(argv[cur_arg + uii], "lin")) {
	    if (model_modifier & MODEL_DMASK) {
	      sprintf(logbuf, "Error: --assoc 'lin' does not make sense with a case/control-specific modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_LIN;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logprint("Error: Improper --assoc mperm syntax.  (Use '--assoc mperm=[value]'.)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    model_modifier |= MODEL_SET_TEST;
	  } else {
	    sprintf(logbuf, "Error: Invalid --assoc parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	model_modifier |= MODEL_ASSOC;
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "djust", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mtest_adjust = 1;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "gc")) {
	    mtest_adjust |= ADJUST_GC;
	  } else if (!strcmp(argv[cur_arg + uii], "log10")) {
	    mtest_adjust |= ADJUST_LOG10;
	  } else if (!strcmp(argv[cur_arg + uii], "qq-plot")) {
	    mtest_adjust |= ADJUST_QQ;
	  } else {
	    sprintf(logbuf, "Error: Invalid --adjust parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 6)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || ((param_ct == 1) && (ii >= ((int32_t)aperm.max) - 1))) {
	  sprintf(logbuf, "Error: Invalid --aperm min permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	aperm.min = ii + 1;
	if (param_ct > 1) {
	  ii = atoi(argv[cur_arg + 2]);
	  // may as well disallow equality since there's no reason not to use
	  // max(T) then...
	  if ((ii <= (int32_t)(aperm.min)) || (ii > APERM_MAX)) {
	    sprintf(logbuf, "Error: Invalid --aperm max permutation count '%s'.%s", argv[cur_arg + 2], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  aperm.max = ii;
	}
	if (param_ct > 2) {
	  if (scan_double(argv[cur_arg + 3], &aperm.alpha)) {
	    sprintf(logbuf, "Error: Invalid --aperm alpha threshold '%s'.%s", argv[cur_arg + 3], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct > 3) {
	    if (scan_double(argv[cur_arg + 4], &aperm.beta) || (aperm.beta <= 0)) {
	      sprintf(logbuf, "Error: Invalid --aperm beta '%s'.%s", argv[cur_arg + 4], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (param_ct > 4) {
	      if (scan_double(argv[cur_arg + 5], &aperm.init_interval)) {
		sprintf(logbuf, "Error: Invalid --aperm initial pruning interval '%s'.%s", argv[cur_arg + 5], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if ((aperm.init_interval < 1) || (aperm.init_interval > 1000000)) {
		sprintf(logbuf, "Error: Invalid --aperm initial pruning interval '%s'.%s", argv[cur_arg + 5], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if (param_ct == 6) {
		if (scan_double(argv[cur_arg + 6], &aperm.interval_slope)) {
		  sprintf(logbuf, "Error: Invalid --aperm pruning interval slope '%s'.%s", argv[cur_arg + 6], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
		if ((aperm.interval_slope < 0) || (aperm.interval_slope > 1)) {
		  sprintf(logbuf, "Error: Invalid --aperm pruning interval slope '%s'.%s", argv[cur_arg + 6], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
	      }
	    }
	  }
	}
      } else if (!memcmp(argptr2, "1-allele", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_2col(&a1alleles, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "2-allele", 9)) {
	if (a1alleles) {
	  logprint("Error: --a2-allele cannot be used with --a1-allele.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_2col(&a2alleles, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'b':
      if (!memcmp(argptr2, "file", 5)) {
	load_params |= 8;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 5)) {
	    logprint("Error: --bfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  sptr = (char*)PROG_NAME_STR;
	}
	if (!(load_params & 16)) {
	  memcpy(strcpya(pedname, sptr), ".bed", 5);
	}
	memcpy(strcpya(mapname, sptr), ".bim", 5);
	memcpy(strcpya(famname, sptr), ".fam", 5);
      } else if (!memcmp(argptr2, "ed", 3)) {
	load_params |= 0x10;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --bed parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "im", 3)) {
	load_params |= 0x20;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --bim parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "merge", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 2) {
	  sprintf(logbuf, "Error: --bmerge must have exactly 1 or 3 parameters.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (param_ct == 3) {
	  if (++jj > FNAMESIZE) {
	    logprint("Error: --bmerge .bed filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --bmerge .bim filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	  jj = strlen(argv[cur_arg + 3]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --bmerge .fam filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename3, argv[cur_arg + 3], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logprint("Error: --bmerge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(mergename1, argv[cur_arg + 1], jj), ".bed", 5);
	  memcpy(memcpya(mergename2, argv[cur_arg + 1], jj), ".bim", 5);
	  memcpy(memcpya(mergename3, argv[cur_arg + 1], jj), ".fam", 5);
	}
	calculation_type |= CALC_MERGE;
	merge_type |= MERGE_BINARY;
      } else if (!memcmp(argptr2, "p-space", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --bp-space minimum bp distance '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        min_bp_space = ii;
      } else if (!memcmp(argptr2, "eta", 4)) {
	logprint("Note: --beta flag deprecated.  Use e.g. '--logistic beta'.\n");
	glm_modifier |= GLM_BETA;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "d", 2)) {
	logprint("Note: --bd flag deprecated.  Use '--mh bd'.\n");
	calculation_type |= CALC_CMH;
	misc_flags |= MISC_CMH_BD;
      } else if (!memcmp(argptr2, "iallelic-only", 14)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "strict")) {
	    misc_flags |= MISC_BIALLELIC_ONLY_STRICT;
	  } else if (!strcmp(argv[cur_arg + uii], "list")) {
	    misc_flags |= MISC_BIALLELIC_ONLY_LIST;
	  } else {
	    sprintf(logbuf, "Error: Invalid --biallelic-only modifier '%s'.%s", argv[cur_arg + uii], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        misc_flags |= MISC_BIALLELIC_ONLY;
      } else if (!memcmp(argptr2, "cf", 3)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = strlen(argv[cur_arg + 1]);
	if (uii > FNAMESIZE - 1) {
	  logprint("Error: --bcf filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], uii + 1);
	load_rare = LOAD_RARE_BCF;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'c':
      if (!memcmp(argptr2, "hr", 3)) {
	if (chrom_flag_present) {
	  sprintf(logbuf, "Error: --chr cannot be used with --autosome{-xy}.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = parse_chrom_ranges(param_ct, '-', &(argv[cur_arg]), chrom_info.chrom_mask, &chrom_info, (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1, argptr);
	if (retval) {
	  goto main_ret_1;
	}
	chrom_info.is_include_stack = 1;
	chrom_flag_present = 1;
      } else if (!memcmp(argptr2, "ompound-genotypes", 18)) {
	logprint("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nand .lgen files are optional if all alleles are single-character).\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ompress", 8)) {
	logprint("Error: --compress flag retired.  Use e.g. 'gzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ounts", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  if (model_modifier & MODEL_QMASK) {
	    sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with --counts.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  logprint("Note: --counts flag deprecated.  Use '--assoc counts' instead.\n");
          model_modifier |= MODEL_ASSOC_COUNTS;
	} else {
	  logprint("Note: --counts flag deprecated.  Use '--freq counts' or --freqx instead.\n");
	}
	misc_flags |= MISC_FREQ_COUNTS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ovar", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "keep-pheno-on-missing-cov")) {
	    uii = 2;
	  } else if (strcmp(argv[cur_arg + 2], "keep-pheno-on-missing-cov")) {
	    sprintf(logbuf, "Error: Invalid --covar parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
          covar_modifier |= COVAR_KEEP_PHENO_ON_MISSING_COV;
	}
	retval = alloc_fname(&covar_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ovar-name", 10)) {
	if (!covar_fname) {
	  logprint("Error: --covar-name must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &covar_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
	covar_modifier |= COVAR_NAME;
      } else if (!memcmp(argptr2, "ovar-number", 12)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (covar_modifier & COVAR_NAME) {
	  logprint("Error: --covar-number cannot be used with --covar-name.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (!covar_fname) {
	  logprint("Error: --covar-number must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &covar_range_list, 1);
	if (retval) {
	  goto main_ret_1;
	}
	covar_modifier |= COVAR_NUMBER;
      } else if (!memcmp(argptr2, "ell", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &model_cell_ct)) {
	  sprintf(logbuf, "Error: Invalid --cell parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "i", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ci parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < 0.01) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --ci confidence interval size s must satisfy 0.01 <= s < 1.%s", errstr_append);
	}
	ci_size = dxx;
      } else if (!memcmp(argptr2, "luster", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cc")) {
            cluster.modifier |= CLUSTER_CC;
	  } else if (!strcmp(argv[cur_arg + uii], "group-avg")) {
	    if (cluster.modifier & CLUSTER_OLD_TIEBREAKS) {
              sprintf(logbuf, "Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    cluster.modifier |= CLUSTER_GROUP_AVG;
	  } else if (!strcmp(argv[cur_arg + uii], "missing")) {
	    cluster.modifier |= CLUSTER_MISSING;
	  } else if (!strcmp(argv[cur_arg + uii], "only2")) {
	    cluster.modifier |= CLUSTER_ONLY2;
	  } else if (!strcmp(argv[cur_arg + uii], "old-tiebreaks")) {
	    if (cluster.modifier & CLUSTER_GROUP_AVG) {
              sprintf(logbuf, "Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
	    cluster.modifier |= CLUSTER_OLD_TIEBREAKS;
	  } else {
            sprintf(logbuf, "Error: Invalid --cluster parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        calculation_type |= CALC_CLUSTER;
      } else if (!memcmp(argptr2, "c", 2)) {
        logprint("Note: --cc flag deprecated.  Use '--cluster cc'.\n");
        cluster.modifier |= CLUSTER_CC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "luster-missing", 15)) {
	if (calculation_type & CALC_CLUSTER) {
	  sprintf(logbuf, "Error: --cluster-missing cannot be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --cluster-missing flag deprecated.  Use '--cluster missing'.\n");
        calculation_type |= CALC_CLUSTER;
        cluster.modifier |= CLUSTER_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "file", 5)) {
        UNSTABLE;
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	sptr = argv[cur_arg + 1];
	uii = strlen(sptr);
	if (uii > (FNAMESIZE - 9)) {
	  logprint("Error: --cfile parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(memcpya(pedname, sptr, uii), ".cnv", 5);
	memcpy(memcpya(famname, sptr, uii), ".fam", 5);
	memcpy(memcpya(mapname, sptr, uii), ".cnv.map", 9);
	load_rare = LOAD_RARE_CNV;
      } else if (!memcmp(argptr2, "nv-count", 9)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_COUNT;
      } else if (!memcmp(argptr2, "nv-del", 7)) {
	UNSTABLE;
	cnv_calc_type |= CNV_DEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-disrupt", 11)) {
	UNSTABLE;
	cnv_overlap_type = CNV_DISRUPT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-dup", 7)) {
	UNSTABLE;
	if (cnv_calc_type & CNV_DEL) {
	  sprintf(logbuf, "Error: --cnv-dup cannot be used with --cnv-del.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_calc_type |= CNV_DUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-enrichment-test", 19)) {
	UNSTABLE;
	if (!cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-enrichment-test must be used with --cnv-count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --cnv-enrichment-test permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_enrichment_test_mperms = ii;
	}
	cnv_calc_type |= CNV_ENRICHMENT_TEST;
      } else if (!memcmp(argptr2, "nv-exclude", 11)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-exclude cannot be used with --cnv-count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_EXCLUDE;
      } else if (!memcmp(argptr2, "nv-exclude-off-by-1", 20)) {
	UNSTABLE;
        cnv_calc_type |= CNV_EXCLUDE_OFF_BY_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-freq-exclude-above", 22)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-exclude-above parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_ABOVE;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-exclude-below", 22)) {
	UNSTABLE;
	if (cnv_freq_type) {
	  logprint("Error: --cnv-freq-exclude-below cannot be used with --cnv-freq-exclude-above.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-exclude-below parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_BELOW;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-exclude-exact", 22)) {
	UNSTABLE;
	if (cnv_freq_type) {
	  logprint("Error: --cnv-freq-exclude-exact cannot be used with\n--cnv-freq-exclude-above/-below.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-exclude-exact parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_EXACT;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-include-exact", 22)) {
	UNSTABLE;
	if (cnv_freq_type) {
	  logprint("Error: --cnv-freq-include-exact cannot be used with\n--cnv-freq-exclude-above/-below/-exact.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-include-exact parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_INCLUDE_EXACT;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-method2", 16)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &cnv_freq_val2) || (cnv_freq_val2 < 0) || (cnv_freq_val2 > 1)) {
	    sprintf(logbuf, "Error: Invalid --cnv-freq-method2 parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	cnv_freq_type |= CNV_FREQ_METHOD2;
	if (cnv_freq_val2 == 0) {
	  // allow >= comparison to be used
	  cnv_freq_val2 = SMALLISH_EPSILON;
	}
      } else if (!memcmp(argptr2, "nv-freq-overlap", 16)) {
	UNSTABLE;
	if (!(cnv_freq_type & CNV_FREQ_FILTER)) {
	  logprint("Error: --cnv-freq-overlap must be used with --cnv-freq-include-exact or\n--cnv-freq-exclude-above/-below/-exact.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_freq_type & CNV_FREQ_METHOD2) {
	  logprint("Error: --cnv-freq-overlap cannot be used with --cnv-freq-method2.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &cnv_freq_val2) || (cnv_freq_val2 < 0) || (cnv_freq_val2 > 1)) {
	    sprintf(logbuf, "Error: Invalid --cnv-freq-overlap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (cnv_freq_val2 == 0) {
	  cnv_freq_val2 = SMALLISH_EPSILON;
	}
	cnv_freq_type |= CNV_FREQ_OVERLAP;
      } else if (!memcmp(argptr2, "nv-indiv-perm", 14)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --cnv-indiv-perm permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_indiv_mperms = ii;
	}
	cnv_calc_type |= CNV_INDIV_PERM;
      } else if (!memcmp(argptr2, "nv-intersect", 13)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-intersect cannot be used with --cnv-count/--cnv-exclude.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_INTERSECT;
      } else if (!memcmp(argptr2, "nv-kb", 6)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx > 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --cnv-kb size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_min_seglen = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
      } else if (!memcmp(argptr2, "nv-list", 8)) {
	UNSTABLE;
	if ((load_rare & (~LOAD_RARE_CNV)) || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	strcpya(pedname, argv[cur_arg + 1]);
	load_rare = LOAD_RARE_CNV;
      } else if (!memcmp(argptr2, "nv-make-map", 12)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-make-map cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "short")) {
            sprintf(logbuf, "Error: Invalid --cnv-make-map parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_calc_type |= CNV_MAKE_MAP;
	} else {
	  cnv_calc_type |= CNV_MAKE_MAP | CNV_MAKE_MAP_LONG;
	}
      } else if (!memcmp(argptr2, "nv-max-kb", 10)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-max-kb cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx > 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --cnv-max-kb size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_max_seglen = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	if (cnv_min_seglen > cnv_max_seglen) {
	  logprint("Error: --cnv-max-kb value cannot be smaller than --cnv-kb value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-max-score", 13)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-max-score cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_max_score)) {
	  sprintf(logbuf, "Error: Invalid --cnv-max-score value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "nv-max-sites", 13)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-max-sites cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], (int32_t*)(&cnv_max_sites))) {
	  sprintf(logbuf, "Error: Invalid --cnv-max-sites parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "nv-overlap", 11)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type == CNV_DISRUPT) {
	  sprintf(logbuf, "Error: --cnv-overlap cannot be used with --cnv-disrupt.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val < 0) || (cnv_overlap_val > 1))  {
	  sprintf(logbuf, "Error: Invalid --cnv-overlap value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_overlap_val != 0) {
	  // ignore --cnv-overlap 0
	  cnv_overlap_type = CNV_OVERLAP;
	}
	if ((cnv_freq_type & CNV_FREQ_FILTER) && (!(cnv_freq_type & (CNV_FREQ_OVERLAP | CNV_FREQ_METHOD2)))) {
	  logprint("Note: --cnv-overlap + --cnv-freq-... deprecated.  Use --cnv-freq-overlap.\n");
	  if (cnv_overlap_val != 0) {
	    cnv_freq_type |= CNV_FREQ_OVERLAP;
	    cnv_freq_val2 = cnv_overlap_val;
	  }
	}
      } else if (!memcmp(argptr2, "nv-region-overlap", 18)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-region-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type) {
	  sprintf(logbuf, "Error: --cnv-region-overlap cannot be used with --cnv-overlap/-disrupt.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val <= 0) || (cnv_overlap_val > 1))  {
	  sprintf(logbuf, "Error: Invalid --cnv-region-overlap value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_overlap_type = CNV_OVERLAP_REGION;
      } else if (!memcmp(argptr2, "nv-score", 9)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-score cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_min_score)) {
	  sprintf(logbuf, "Error: Invalid --cnv-score value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_min_score > cnv_max_score) {
	  logprint("Error: --cnv-score value cannot be greater than --cnv-max-score value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-sites", 9)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-sites cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], (int32_t*)(&cnv_min_sites))) {
	  sprintf(logbuf, "Error: Invalid --cnv-sites parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_min_sites > cnv_max_sites) {
	  logprint("Error: --cnv-sites value cannot be greater than --cnv-max-sites value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-subset", 10)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-subset cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-subset must be used with --cnv-intersect/-exclude/-count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "nv-test", 8)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-test cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "1sided")) {
	    uii = 2;
	    cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
	  } else if (!strcmp(argv[cur_arg + 1], "2sided")) {
	    uii = 2;
	    cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
	  } else {
	    uii = 1;
	    if (!strcmp(argv[cur_arg + 2], "1sided")) {
	      cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
	    } else if (!strcmp(argv[cur_arg + 2], "2sided")) {
	      cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
	    } else {
	      sprintf(logbuf, "Error: Invalid --cnv-test parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	} else {
	  uii = 1;
	}
	ii = atoi(argv[cur_arg + uii]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-test permutation count '%s'.%s", argv[cur_arg + uii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_test_mperms = ii;
	cnv_calc_type |= CNV_TEST;
      } else if (!memcmp(argptr2, "nv-test-1sided", 15)) {
	UNSTABLE;
	if (cnv_calc_type & CNV_TEST_FORCE_2SIDED) {
	  logprint("Error: --cnv-test cannot be both 1-sided and 2-sided at the same time.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --cnv-test-1sided flag deprecated.  Use '--cnv-test 1sided'.\n");
	cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
      } else if (!memcmp(argptr2, "nv-test-2sided", 15)) {
	UNSTABLE;
	if (cnv_calc_type & CNV_TEST_FORCE_1SIDED) {
	  logprint("Error: --cnv-test cannot be both 1-sided and 2-sided at the same time.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --cnv-test-2sided flag deprecated.  Use '--cnv-test 2sided'.\n");
	cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
      } else if (!memcmp(argptr2, "nv-test-region", 15)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-test-region cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --cnv-test-region permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_test_region_mperms = ii;
	}
	cnv_calc_type |= CNV_TEST_REGION;
      } else if (!memcmp(argptr2, "nv-test-window", 15)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-test-window cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001)) {
	  sprintf(logbuf, "Error: Invalid --cnv-test-window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dxx *= 1000;
	if (dxx > 2147483647) {
	  cnv_test_window = 0x7fffffff;
	} else {
	  cnv_test_window = (int32_t)(dxx * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "nv-union-overlap", 17)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-union-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type) {
	  sprintf(logbuf, "Error: --cnv-union-overlap cannot be used with --cnv-{region-}overlap/-disrupt.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val <= 0) || (cnv_overlap_val > 1)) {
	  sprintf(logbuf, "Error: Invalid --cnv-union-overlap value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_overlap_type = CNV_OVERLAP_UNION;
      } else if (!memcmp(argptr2, "nv-write", 9)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-write cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "freq")) {
            sprintf(logbuf, "Error: Invalid --cnv-write parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (!(cnv_freq_val & CNV_FREQ_METHOD2)) {
	    sprintf(logbuf, "Error: --cnv-write 'freq' modifier must be used with --cnv-freq-method2.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_calc_type |= CNV_WRITE_FREQ;
	}
	cnv_calc_type |= CNV_WRITE;
      } else if (!memcmp(argptr2, "nv-write-freq", 14)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-write freq cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!(cnv_freq_val & CNV_FREQ_METHOD2)) {
	  sprintf(logbuf, "Error: --cnv-write 'freq' modifier must be used with --cnv-freq-method2.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (!(cnv_calc_type & CNV_WRITE)) {
	  sprintf(logbuf, "Error: --cnv-write-freq must be used with --cnv-write.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --cnv-write-freq flag deprecated.  Use '--cnv-write freq'.\n");
	cnv_calc_type |= CNV_WRITE_FREQ;
      } else if (!memcmp(argptr2, "onsensus-match", 15)) {
        logprint("Note: --consensus-match flag deprecated.  Use '--homozyg consensus-match'.\n");
	homozyg.modifier |= HOMOZYG_CONSENSUS_MATCH;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ondition", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp("dominant", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_DOMINANT;
	  } else if (!strcmp("recessive", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_RECESSIVE;
	  } else {
	    uii = 2;
            if (!strcmp("dominant", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_DOMINANT;
	    } else if (!strcmp("recessive", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_RECESSIVE;
	    } else {
	      sprintf(logbuf, "Error: Invalid --condition parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	if (alloc_string(&condition_mname, argv[cur_arg + uii])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ondition-list", 14)) {
	if (condition_mname) {
	  logprint("Error: --condition-list cannot be used with --condition.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 2) {
	  if (!strcmp("dominant", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_DOMINANT;
	  } else if (!strcmp("recessive", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_RECESSIVE;
	  } else {
	    uii = 2;
            if (!strcmp("dominant", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_DOMINANT;
	    } else if (!strcmp("recessive", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_RECESSIVE;
	    } else {
	      sprintf(logbuf, "Error: Invalid --condition-list parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	retval = alloc_fname(&condition_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "omplement-sets", 15)) {
        set_info.modifier |= SET_COMPLEMENTS | SET_C_PREFIX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "onst-fid", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&const_fid, param_ct? argv[cur_arg + 1] : "0")) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ase-only", 9)) {
	logprint("Note: --case-only flag deprecated.  Use '--fast-epistasis case-only'.\n");
        epi_info.modifier |= EPI_FAST_CASE_ONLY;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "m-map", 6)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  // must contain exactly one '#'
          sptr = strchr(argv[cur_arg + 1], '#');
	  if (!sptr) {
            sprintf(logbuf, "Error: --cm-map requires either a '#' in the filename pattern, or a chromosome\ncode as the second parameter.%s", errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
          if (strchr(&(sptr[1]), '#')) {
	    sprintf(logbuf, "Error: Multiple '#'s in --cm-map filename pattern.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          if (alloc_string(&cm_map_fname, argv[cur_arg + 1])) {
	    goto main_ret_NOMEM;
	  }
	} else {
          retval = alloc_fname(&cm_map_fname, argv[cur_arg + 1], argptr, 0);
          if (retval) {
            goto main_ret_1;
	  }
          if (alloc_string(&cm_map_chrname, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "heck-sex", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
          if (scan_double(argv[cur_arg + 1], &check_sex_fthresh) || (check_sex_fthresh <= 0.0) || ((param_ct == 1) && (check_sex_fthresh > check_sex_mthresh))) {
	    sprintf(logbuf, "Error: Invalid --check-sex female F-statistic ceiling '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    if (scan_double(argv[cur_arg + 2], &check_sex_mthresh) || (check_sex_mthresh < check_sex_fthresh) || (check_sex_mthresh >= 1.0)) {
	      sprintf(logbuf, "Error: Invalid --check-sex male F-statistic floor '%s'.%s", argv[cur_arg + 1], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	} else {
	  logprint("Warning: --check-sex will use default 0.2 and 0.8 thresholds.  This is not\nrecommended.\n");
	}
        calculation_type |= CALC_SEXCHECK;
      } else if (!memcmp(argptr2, "lump", 5)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_and_flatten(&clump_info.fnames_flattened, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
	calculation_type |= CALC_CLUMP;
        clump_info.fname_ct = param_ct;
      } else if (!memcmp(argptr2, "lump-allow-overlap", 19)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-allow-overlap must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        clump_info.modifier |= CLUMP_ALLOW_OVERLAP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-annotate", 14)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-annotate must be used with --clump.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_and_flatten(&clump_info.annotate_flattened, &(argv[cur_arg + 1]), param_ct);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "lump-best", 10)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-best must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        clump_info.modifier |= CLUMP_BEST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-field", 11)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-field must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&clump_info.field_name, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "lump-index-first", 17)) {
        if (clump_info.fname_ct < 2) {
	  sprintf(logbuf, "Error: --clump-index-first requires multiple --clump files.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        clump_info.modifier |= CLUMP_INDEX_FIRST;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-kb", 8)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-kb must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001)) {
	  sprintf(logbuf, "Error: Invalid --clump-kb parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dxx *= 1000;
        if (dxx > 2147483647) {
	  clump_info.bp_radius = 0x7fffffff;
	} else {
	  clump_info.bp_radius = ((int32_t)(dxx * (1 + SMALL_EPSILON))) - 1;
	}
      } else if (!memcmp(argptr2, "lump-p1", 8)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-p1 must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0) || (dxx > 1)) {
	  sprintf(logbuf, "Error: Invalid --clump-p1 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	clump_info.p1 = dxx;
      } else if (!memcmp(argptr2, "lump-p2", 8)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-p2 must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < clump_info.p1) || (dxx > 1)) {
	  sprintf(logbuf, "Error: Invalid --clump-p2 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	clump_info.p2 = dxx;
      } else if (!memcmp(argptr2, "lump-r2", 8)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-r2 must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx >= 1)) {
	  sprintf(logbuf, "Error: Invalid --clump-r2 parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	clump_info.r2 = dxx;
      } else if (!memcmp(argptr2, "lump-range", 11)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-range must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&clump_info.range_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "lump-range-border", 18)) {
	if (!clump_info.range_fname) {
	  logprint("Error: --clump-range-border must be used with --clump-range.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(logbuf, "Error: Invalid --clump-range-border parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (dxx > 2147483.647) {
	  clump_info.range_border = 0x7fffffff;
	} else {
	  clump_info.range_border = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "lump-replicate", 16)) {
        if (clump_info.fname_ct < 2) {
	  sprintf(logbuf, "Error: --clump-replicate requires multiple --clump files.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        clump_info.modifier |= CLUMP_REPLICATE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "lump-verbose", 14)) {
        if (!clump_info.fname_ct) {
	  logprint("Error: --clump-verbose must be used with --clump.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        clump_info.modifier |= CLUMP_VERBOSE;
        goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'd':
      if (!memcmp(argptr2, "ebug", 5)) {
	g_debug_on = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ata", 4)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x80;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 8)) {
	    logprint("Error: --data parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  sptr = (char*)PROG_NAME_STR;
	}
	memcpy(strcpya(pedname, sptr), ".gen", 5);
	// cheating: this is of course more like a .fam file
	memcpy(strcpya(mapname, sptr), ".sample", 8);
      } else if (!memcmp(argptr2, "ecompress", 10)) {
	logprint("Error: --decompress flag retired.  Use e.g. 'gunzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "istance", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (dist_calc_type & DISTANCE_BIN) {
	      sprintf(logbuf, "Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (dist_calc_type & DISTANCE_GZ) {
	      sprintf(logbuf, "Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "ibs")) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logprint("Error: Duplicate --distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "1-ibs")) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logprint("Error: Duplicate --distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "allele-ct")) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logprint("Error: Duplicate --distance 'allele-ct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!strcmp(argv[cur_arg + uii], "flat-missing")) {
	    dist_calc_type |= DISTANCE_FLAT_MISSING;
	  } else {
	    sprintf(logbuf, "Error: Invalid --distance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(dist_calc_type & DISTANCE_TYPEMASK)) {
	  dist_calc_type |= DISTANCE_ALCT;
	}
	calculation_type |= CALC_DISTANCE;
      } else if (!memcmp(argptr2, "istance-exp", 12)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &exponent)) {
	  sprintf(logbuf, "Error: Invalid --distance-exp parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "istance-matrix", 15)) {
	if (exponent != 0.0) {
	  logprint("Error: --distance-matrix cannot be used with --distance-exp.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	  sprintf(logbuf, "Error: --distance-matrix flag cannot be used with '--distance 1-ibs'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_PLINK1_DISTANCE_MATRIX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ummy", 5)) {
	if (load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 6)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy individual count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_indiv_ct = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy variant count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_marker_ct = ii;
        for (uii = 3; uii <= param_ct; uii++) {
	  if (match_upper(argv[cur_arg + uii], "ACGT")) {
	    if (dummy_flags & (DUMMY_1234 | DUMMY_12)) {
	      sprintf(logbuf, "Error: --dummy 'acgt' modifier cannot be used with '1234' or '12'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_ACGT;
	  } else if (!strcmp(argv[cur_arg + uii], "1234")) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_12)) {
	      sprintf(logbuf, "Error: --dummy '1234' modifier cannot be used with 'acgt' or '12'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_1234;
	  } else if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_1234)) {
	      sprintf(logbuf, "Error: --dummy '12' modifier cannot be used with 'acgt' or '1234'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_12;
	  } else if (!strcmp(argv[cur_arg + uii], "scalar-pheno")) {
	    dummy_flags |= DUMMY_SCALAR_PHENO;
	  } else {
	    if ((dummy_flags & DUMMY_MISSING_PHENO) || scan_double(argv[cur_arg + uii], &dxx) || (dxx < 0.0) || (dxx > 1.0)) {
	      sprintf(logbuf, "Error: Invalid --dummy parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (dummy_flags & DUMMY_MISSING_GENO) {
	      dummy_missing_pheno = dxx;
	      dummy_flags |= DUMMY_MISSING_PHENO;
	    } else {
	      dummy_missing_geno = dxx;
	      dummy_flags |= DUMMY_MISSING_GENO;
	    }
	  }
	}
	load_rare = LOAD_RARE_DUMMY;
      } else if (!memcmp(argptr2, "ummy-coding", 12)) {
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --dummy-coding cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "no-round")) {
	    uii = 2;
	    write_covar_modifier |= WRITE_COVAR_DUMMY_NO_ROUND;
	  } else {
	    if (param_ct == 2) {
              if (strcmp(argv[cur_arg + 2], "no-round")) {
		sprintf(logbuf, "Error: Invalid --dummy-coding parameter sequence.%s", errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      write_covar_modifier |= WRITE_COVAR_DUMMY_NO_ROUND;
	    }
	    uii = 1;
	  }
	  if (uii <= param_ct) {
            ii = atoi(argv[cur_arg + uii]);
	    if (ii < 3) {
	      sprintf(logbuf, "Error: Invalid --dummy-coding max categories parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            write_covar_dummy_max_categories = ii;
	  }
	}
	write_covar_modifier |= WRITE_COVAR_DUMMY;
      } else if (!memcmp(argptr2, "ominant", 8)) {
	logprint("Note: --dominant flag deprecated.  Use e.g. '--linear dominant' (and\n'--condition-list [filename] dominant' to change covariate coding).\n");
	glm_modifier |= GLM_DOMINANT | GLM_CONDITION_DOMINANT;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ouble-id", 9)) {
        if (const_fid) {
	  sprintf(logbuf, "Error: --double-id cannot be used with --const-fid.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        misc_flags |= MISC_DOUBLE_ID;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "prime", 6)) {
	logprint("Note: --dprime flag deprecated.  Use e.g. '--r2 dprime'.\n");
	ld_info.modifier |= LD_DPRIME;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'e':
      if (!memcmp(argptr2, "xtract", 7)) {
	if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "--extract cannot be used with a .cnv fileset.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&extractname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "xclude", 7)) {
	if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "--exclude cannot be used with a .cnv fileset.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&excludename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "xclude-snp", 11)) {
	if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --exclude-snp cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb\nand --to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
        misc_flags |= MISC_EXCLUDE_MARKERNAME_SNP;
      } else if (!memcmp(argptr2, "xclude-snps", 12)) {
	if (markername_snp) {
	  sprintf(logbuf, "Error: --exclude-snps cannot be used with --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --exclude-snps cannot be used with a .cnv fileset.  Use\n--from-bp/-kb/-mb and --to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &snps_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
        misc_flags |= MISC_EXCLUDE_MARKERNAME_SNP;
      } else if (!memcmp(argptr2, "pistasis", 9)) {
	if (epi_info.modifier & EPI_FAST_CASE_ONLY) {
	  sprintf(logbuf, "Error: --epistasis cannot be used with --case-only.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
          if (!strcmp(argv[cur_arg + 1], "set-by-set")) {
	    epi_info.modifier |= EPI_SET_BY_SET;
	  } else if (!strcmp(argv[cur_arg + 1], "set-by-all")) {
	    epi_info.modifier |= EPI_SET_BY_ALL;
	  } else {
	    sprintf(logbuf, "Error: Invalid --epistasis modifier '%s'.%s", argv[cur_arg + 1], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        epi_info.modifier |= EPI_REG;
	calculation_type |= CALC_EPI;
      } else if (!memcmp(argptr2, "pistasis-summary-merge", 23)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	// .summary.32767 = 14 chars
	retval = alloc_fname(&epi_info.summary_merge_prefix, argv[cur_arg + 1], argptr, 14);
	if (retval) {
	  goto main_ret_1;
	}
        ii = atoi(argv[cur_arg + 2]);
	if ((ii < 2) || (ii > PARALLEL_MAX)) {
	  sprintf(logbuf, "Error: Invalid --epistasis-summary-merge job count '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        epi_info.summary_merge_ct = ii;
      } else if (!memcmp(argptr2, "pi1", 4)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0)) {
	  sprintf(logbuf, "Error: Invalid --epi1 parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	epi_info.epi1 = dxx;
      } else if (!memcmp(argptr2, "pi2", 4)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0) || (dxx >= 1)) {
	  sprintf(logbuf, "Error: Invalid --epi2 parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	epi_info.epi2 = dxx;
      } else if (!memcmp(argptr2, "xclude-before-extract", 22)) {
        logprint("Note: --exclude-before-extract has no effect.\n");
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'f':
      if (!memcmp(argptr2, "ile", 4)) {
	if (load_params & 0x3f9) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 1;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 5)) {
	    logprint("Error: --file parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  if (!(load_params & 2)) {
	    memcpy(strcpya(pedname, argv[cur_arg + 1]), ".ped", 5);
	  }
	  if (!(load_params & 4)) {
	    memcpy(strcpya(mapname, argv[cur_arg + 1]), ".map", 5);
	  }
	}
      } else if (!memcmp(argptr2, "am", 3)) {
	if (load_params & 0x3c7) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x40;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --fam parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(famname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "ilter", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&filtername, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        retval = alloc_and_flatten(&filtervals_flattened, &(argv[cur_arg + 2]), param_ct - 1);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ilter-cases", 12)) {
	filter_binary |= FILTER_BINARY_CASES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-controls", 15)) {
	if (filter_binary & FILTER_BINARY_CASES) {
	  logprint("Error: --filter-cases and --filter-controls cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_binary |= FILTER_BINARY_CONTROLS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-females", 14)) {
	filter_binary |= FILTER_BINARY_FEMALES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-males", 12)) {
	if (filter_binary & FILTER_BINARY_FEMALES) {
	  logprint("Error: --filter-males and --filter-females cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_binary |= FILTER_BINARY_MALES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-founders", 15)) {
	filter_binary |= FILTER_BINARY_FOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-nonfounders", 18)) {
	if (filter_binary & FILTER_BINARY_FOUNDERS) {
	  logprint("Error: --filter-founders and --filter-nonfounders cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_binary |= FILTER_BINARY_NONFOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "req", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "counts")) {
            sprintf(logbuf, "Error: Invalid --freq parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  misc_flags |= MISC_FREQ_COUNTS;
	}
	calculation_type |= CALC_FREQ;
	if (misc_flags & MISC_FREQ_COUNTS) {
	  // --keep-allele-order also set for backward compatibility
	  misc_flags |= MISC_KEEP_ALLELE_ORDER;
	}
      } else if (!memcmp(argptr2, "reqx", 5)) {
	if (calculation_type & CALC_FREQ) {
	  sprintf(logbuf, "Error: --freqx cannot be used with --freq.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_FREQ;
	misc_flags |= MISC_FREQX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rom", 4)) {
	if (chrom_flag_present) {
	  sprintf(logbuf, "Error: --from cannot be used with --autosome{-xy} or --{not-}chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --from cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb\ninstead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_from, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if ((!memcmp(argptr2, "rom-bp", 7)) || (!memcmp(argptr2, "rom-kb", 7)) || (!memcmp(argptr2, "rom-mb", 7))) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --from-bp/-kb/-mb cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argptr2[4];
	if (cc == 'b') {
	  if (atoiz(argv[cur_arg + 1], &marker_pos_start)) {
	    sprintf(logbuf, "Error: Invalid --from-bp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  if (marker_pos_start != -1) {
	    logprint("Error: Multiple --from-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (scan_double(argv[cur_arg + 1], &dxx)) {
	    sprintf(logbuf, "Error: Invalid --from-kb/-mb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    marker_pos_start = 0;
	  } else if (dxx > 2147483647) {
	    marker_pos_start = 0x7fffffff;
	  } else {
	    marker_pos_start = (int32_t)(dxx * (1 + SMALL_EPSILON));
	  }
	}
      } else if (!memcmp(argptr2, "isher", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --fisher cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --fisher flag deprecated.  Use '--assoc fisher' or '--model fisher'.\n");
	model_modifier |= MODEL_ASSOC | MODEL_FISHER | MODEL_ASSOC_FDEPR;
	calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "id", 3)) {
        logprint("Note: --fid flag deprecated.  Use '--recode vcf-fid'.\n");
	recode_modifier |= RECODE_FID;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lip", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&flip_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "lip-subset", 11)) {
	UNSTABLE;
	if (!flip_fname) {
          sprintf(logbuf, "Error: --flip-subset must be used with --flip.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (allelexxxx) {
	  // fix for this is too messy to be worthwhile
	  logprint("Error: --flip-subset cannot be used with --allele1234/ACGT.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&flip_subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ilter-attrib", 13)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&filter_attrib_fname, argv[cur_arg + 1], argptr, 0);
        if (retval) {
          goto main_ret_1;
	}
	if (param_ct == 2) {
	  // force comma-terminated string to simplify parsing
	  uii = strlen(argv[cur_arg + 2]);
	  filter_attrib_liststr = (char*)malloc(uii + 2);
	  if (!filter_attrib_liststr) {
	    goto main_ret_NOMEM;
	  }
          memcpy(filter_attrib_liststr, argv[cur_arg + 2], uii);
	  memcpy(&(filter_attrib_liststr[uii]), ",", 2);
	}
      } else if (!memcmp(argptr2, "ilter-attrib-indiv", 19)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&filter_attrib_indiv_fname, argv[cur_arg + 1], argptr, 0);
        if (retval) {
          goto main_ret_1;
	}
	if (param_ct == 2) {
	  uii = strlen(argv[cur_arg + 2]);
	  filter_attrib_indiv_liststr = (char*)malloc(uii + 2);
	  if (!filter_attrib_indiv_liststr) {
	    goto main_ret_NOMEM;
	  }
          memcpy(filter_attrib_indiv_liststr, argv[cur_arg + 2], uii);
	  memcpy(&(filter_attrib_indiv_liststr[uii]), ",", 2);
	}
      } else if (!memcmp(argptr2, "ast-epistasis", 14)) {
	if (epi_info.modifier & EPI_REG) {
	  logprint("Error: --fast-epistasis cannot be used with --epistasis.\n");
          goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-ueki")) {
	    if (epi_info.modifier & (EPI_FAST_BOOST | EPI_FAST_JOINT_EFFECTS)) {
	      sprintf(logbuf, "Error: --fast-epistasis 'no-ueki' modifier cannot be used with '%s'.%s", (epi_info.modifier & EPI_FAST_BOOST)? "boost" : "joint-effects", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    epi_info.modifier |= EPI_FAST_NO_UEKI;
	  } else if (!strcmp(argv[cur_arg + uii], "boost")) {
	    if (epi_info.modifier & (EPI_FAST_NO_UEKI | EPI_FAST_JOINT_EFFECTS)) {
	      sprintf(logbuf, "Error: --fast-epistasis 'boost' modifier cannot be used with '%s'.%s", (epi_info.modifier & EPI_FAST_NO_UEKI)? "no-ueki" : "joint-effects", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (epi_info.modifier & EPI_FAST_CASE_ONLY) {
	      logprint("Error: --fast-epistasis boost does not have a case-only mode.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    epi_info.modifier |= EPI_FAST_BOOST;
	  } else if (!strcmp(argv[cur_arg + uii], "joint-effects")) {
	    if (epi_info.modifier & (EPI_FAST_NO_UEKI | EPI_FAST_BOOST)) {
	      sprintf(logbuf, "Error: --fast-epistasis 'joint-effects' modifier cannot be used with '%s'.%s", (epi_info.modifier & EPI_FAST_NO_UEKI)? "no-ueki" : "boost", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    epi_info.modifier |= EPI_FAST_JOINT_EFFECTS;
	  } else if (!strcmp(argv[cur_arg + uii], "case-only")) {
	    if (epi_info.modifier & EPI_FAST_BOOST) {
              logprint("Error: --fast-epistasis boost does not have a case-only mode.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    epi_info.modifier |= EPI_FAST_CASE_ONLY;
	  } else if (!strcmp(argv[cur_arg + uii], "set-by-set")) {
            if (!(epi_info.modifier & EPI_SET_BY_ALL)) {
	      epi_info.modifier |= EPI_SET_BY_SET;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "set-by-all")) {
	    if (epi_info.modifier & EPI_SET_BY_SET) {
	      epi_info.modifier -= EPI_SET_BY_SET;
	    }
            epi_info.modifier |= EPI_SET_BY_ALL;
	  } else if (!strcmp(argv[cur_arg + uii], "nop")) {
	    epi_info.modifier |= EPI_FAST_NO_P_VALUE;
	  } else {
	    sprintf(logbuf, "Error: Invalid --fast-epistasis modifier '%s'.%s", argv[cur_arg + uii], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        epi_info.modifier |= EPI_FAST;
	calculation_type |= CALC_EPI;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'g':
      if (!memcmp(argptr2, "eno", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &geno_thresh)) {
	    sprintf(logbuf, "Error: Invalid --geno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((geno_thresh < 0.0) || (geno_thresh > 1.0)) {
	    sprintf(logbuf, "Error: Invalid --geno parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  geno_thresh = 0.1;
	}
      } else if (!memcmp(argptr2, "en", 3)) {
	if (load_rare || (load_params & 0x7f)) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x100;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --gen parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if ((!memcmp(argptr2, "enome", 6)) || (!memcmp(argptr2, "enome gz", 9))) {
	if (argptr2[5] == ' ') {
	  kk = 1;
	  genome_modifier |= GENOME_OUTPUT_GZ;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5 - kk)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "gz")) {
            genome_modifier |= GENOME_OUTPUT_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "rel-check")) {
            genome_modifier |= GENOME_REL_CHECK;
	  } else if (!strcmp(argv[cur_arg + uii], "full")) {
	    genome_modifier |= GENOME_OUTPUT_FULL;
	  } else if (!strcmp(argv[cur_arg + uii], "unbounded")) {
	    genome_modifier |= GENOME_IBD_UNBOUNDED;
	  } else if (!strcmp(argv[cur_arg + uii], "nudge")) {
            genome_modifier |= GENOME_NUDGE;
	  } else {
	    sprintf(logbuf, "Error: Invalid --genome parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_GENOME;
      } else if (!memcmp(argptr2, "enome-full", 11)) {
	if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --genome-full must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --genome-full flag deprecated.  Use '--genome full'.\n");
	genome_modifier |= GENOME_OUTPUT_FULL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "roupdist", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  groupdist_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((groupdist_iters < 2) || (groupdist_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --groupdist jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    ii = atoi(argv[cur_arg + 2]);
	    if (ii <= 0) {
	      sprintf(logbuf, "Error: Invalid --groupdist jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    groupdist_d = ii;
	  }
	}
	calculation_type |= CALC_GROUPDIST;
      } else if (!memcmp(argptr2, "rm", 3)) {
	logprint("Error: --grm has been retired due to inconsistent meaning across GCTA versions.\nUse --grm-gz or --grm-bin.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "rm-gz", 6)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 8)) {
	    logprint("Error: --grm-gz parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
        load_rare = LOAD_RARE_GRM;
      } else if (!memcmp(argptr2, "rm-bin", 7)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 11)) {
	    logprint("Error: --grm-bin parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
        load_rare = LOAD_RARE_GRM_BIN;
      } else if (!memcmp(argptr2, "xe", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  logprint("Error: --gxe must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --gxe parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  gxe_mcovar = ii;
	} else {
	  gxe_mcovar = 1;
	}
	calculation_type |= CALC_GXE;
      } else if (!memcmp(argptr2, "enedrop", 8)) {
	if (model_modifier & MODEL_QMASK) {
	  sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with --genedrop.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --genedrop flag deprecated.  Use e.g. '--model genedrop'.\n");
	model_modifier |= MODEL_GENEDROP;
	glm_modifier |= GLM_GENEDROP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "c", 2)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --gc must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --gc flag deprecated.  Use '--adjust gc'.\n");
	mtest_adjust |= ADJUST_GC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "file", 5)) {
	UNSTABLE;
	if (load_rare || (load_params & (~0x40))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	sptr = argv[cur_arg + 1];
	uii = strlen(sptr);
	if (uii > (FNAMESIZE - 6)) {
	  logprint("Error: --gfile parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(memcpya(pedname, sptr, uii), ".gvar", 6);
	if (!(load_params & 0x40)) {
	  memcpy(memcpya(famname, sptr, uii), ".fam", 5);
	}
	memcpy(memcpya(mapname, sptr, uii), ".map", 5);
	load_rare = LOAD_RARE_GVAR;
      } else if (!memcmp(argptr2, "enome-lists", 12)) {
	logprint("Error: --genome-lists flag retired.  Use --parallel.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "enome-minimal", 14)) {
	logprint("Error: --genome-minimal flag retired.  Use '--genome gz'.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if ((!memcmp(argptr2, "roup-avg", 9)) || (!memcmp(argptr2, "roup-average", 13))) {
        if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --group-avg must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (cluster.modifier & CLUSTER_OLD_TIEBREAKS) {
	  sprintf(logbuf, "Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        sprintf(logbuf, "Note: --%s flag deprecated.  Use '--cluster group-avg'.\n", argptr);
	logprintb();
	cluster.modifier |= CLUSTER_GROUP_AVG;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "enotypic", 9)) {
	if (glm_modifier & GLM_DOMINANT) {
	  logprint("Error: --genotypic cannot be used with --dominant.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --genotypic flag deprecated.  Use e.g. '--linear genotypic'.\n");
	glm_modifier |= GLM_GENOTYPIC;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ene", 4)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_and_flatten(&set_info.genekeep_flattened, &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ene-all", 8)) {
	if (set_info.genekeep_flattened) {
          sprintf(logbuf, "Error: --gene-all cannot be used with --gene.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        set_info.modifier |= SET_GENE_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ap", 3)) {
	if ((epi_info.modifier & (EPI_FAST | EPI_FAST_CASE_ONLY)) != (EPI_FAST | EPI_FAST_CASE_ONLY)) {
	  sprintf(logbuf, "Error: --gap must be used with '--fast-epistasis case-only'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(logbuf, "Error: Invalid --gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (dxx > 2147483.647) {
	  epi_info.case_only_gap = 2147483647;
	} else {
          epi_info.case_only_gap = (int32_t)(dxx * 1000 * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "plink", 6)) {
        misc_flags |= MISC_GPLINK;
      } else if (!memcmp(argptr2, "ates", 5)) {
        logprint("Error: --gates is not implemented yet.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
        goto main_ret_1;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'h':
      if (!memcmp(argptr2, "we", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ujj = 0;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "midp")) {
	    hwe_modifier |= HWE_THRESH_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "include-nonctrl")) {
	    hwe_modifier |= HWE_THRESH_ALL;
	  } else {
	    if (scan_double(argv[cur_arg + uii], &hwe_thresh) || ujj) {
	      sprintf(logbuf, "Error: Invalid --hwe parameter sequence.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
            ujj = 1;
            if ((hwe_thresh < 0.0) || (hwe_thresh >= 1.0)) {
	      sprintf(logbuf, "Error: Invalid --hwe threshold '%s' (must be between 0 and 1).%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	if (!ujj) {
	  sprintf(logbuf, "Error: --hwe now requires a p-value threshold.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if ((hwe_modifier & HWE_THRESH_MIDP) && (hwe_thresh >= 0.5)) {
	  sprintf(logbuf, "Error: --hwe threshold must be smaller than 0.5 when using mid-p adjustment.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "we-all", 7)) {
	logprint("Note: --hwe-all flag deprecated.  Use '--hwe include-nonctrl'.\n");
	hwe_modifier |= HWE_THRESH_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et", 3)) {
	sprintf(logbuf, "Error: --het provisionally retired.  Contact us if --ibc is unsatisfactory.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "ardy", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "midp")) {
            sprintf(logbuf, "Error: Invalid --hardy parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
          hwe_modifier |= HWE_MIDP;
	}
	calculation_type |= CALC_HARDY;
      } else if (!memcmp(argptr2, "we2", 4)) {
	sprintf(logbuf, "Error: --hwe2 retired.  Use the --hwe exact test.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "ardy2", 6)) {
	sprintf(logbuf, "Error: --hardy2 retired.  Use the exact test-based --hardy report.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "omozyg", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "group")) {
	    if (homozyg.modifier & HOMOZYG_GROUP_VERBOSE) {
	      logprint("Error: --homozyg 'group' and 'group-verbose' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    homozyg.modifier |= HOMOZYG_GROUP;
	  } else if (!strcmp(argv[cur_arg + uii], "group-verbose")) {
	    if (homozyg.modifier & HOMOZYG_GROUP) {
	      logprint("Error: --homozyg 'group' and 'group-verbose' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    homozyg.modifier |= HOMOZYG_GROUP_VERBOSE;
	  } else if (!strcmp(argv[cur_arg + uii], "consensus-match")) {
	    homozyg.modifier |= HOMOZYG_CONSENSUS_MATCH;
	  } else if (!strcmp(argv[cur_arg + uii], "extend")) {
	    homozyg.modifier |= HOMOZYG_EXTEND;
	  } else if (!strcmp(argv[cur_arg + uii], "subtract-1-from-lengths")) {
            homozyg.modifier |= HOMOZYG_OLD_LENGTHS;
	  } else {
	    sprintf(logbuf, "Error: Invalid --homozyg parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-snp", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --homozyg-snp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_HOMOZYG;
	homozyg.min_snp = ii;
      } else if (!memcmp(argptr2, "omozyg-kb", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < SMALL_EPSILON) || (dxx >= (2147483.647 * (1 + SMALL_EPSILON)))) {
	  sprintf(logbuf, "Error: Invalid --homozyg-kb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_HOMOZYG;
	// round up
	homozyg.min_bases = 1 + (uint32_t)((int32_t)(dxx * 1000 * (1 - SMALL_EPSILON)));
      } else if (!memcmp(argptr2, "omozyg-density", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx >= 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-density parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.max_bases_per_snp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
      } else if (!memcmp(argptr2, "omozyg-gap", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx >= 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.max_gap = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
      } else if (!memcmp(argptr2, "omozyg-het", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-het parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii && (homozyg.modifier & HOMOZYG_EXTEND)) {
	  sprintf(logbuf, "Error: --homozyg-het with a nonzero parameter cannot be used with --homozyg\nextend.  For fine-grained control over heterozygote frequency, use\n--homozyg-window-snp and --homozyg-window-het instead.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_HOMOZYG;
        homozyg.max_hets = ii;
      } else if (!memcmp(argptr2, "omozyg-window-snp", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-snp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.window_size = ii;
      } else if (!memcmp(argptr2, "omozyg-window-kb", 17)) {
        logprint("Error: --homozyg-window-kb flag provisionally retired, since it had no effect\nin PLINK 1.07.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "omozyg-window-het", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-het parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.window_max_hets = ii;
      } else if (!memcmp(argptr2, "omozyg-window-missing", 22)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-missing parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.window_max_missing = ii;
      } else if (!memcmp(argptr2, "omozyg-window-threshold", 24)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-threshold parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.hit_threshold = dxx;
      } else if (!memcmp(argptr2, "omozyg-match", 13)) {
	if (!(homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE))) {
          homozyg.modifier |= HOMOZYG_GROUP;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-match parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	homozyg.overlap_min = dxx;
      } else if (!memcmp(argptr2, "omozyg-group", 13)) {
	if (homozyg.modifier & HOMOZYG_GROUP_VERBOSE) {
	  logprint("Note: --homozyg-group deprecated, and superseded by --homozyg group-verbose.\n");
	} else {
	  logprint("Note: --homozyg-group flag deprecated.  Use '--homozyg group'.\n");
	  homozyg.modifier |= HOMOZYG_GROUP;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omozyg-verbose", 15)) {
	if (!(homozyg.modifier & HOMOZYG_GROUP)) {
	  logprint("Error: --homozyg-verbose must be used with --homozyg group.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --homozyg-verbose flag deprecated.  Use '--homozyg group-verbose'.\n");
	homozyg.modifier = (homozyg.modifier & (~HOMOZYG_GROUP)) | HOMOZYG_GROUP_VERBOSE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omozyg-include-missing", 23)) {
        logprint("Error: --homozyg-include-missing flag provisionally retired, since it had no\neffect in PLINK 1.07.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ethom", 6)) {
	if (!(glm_modifier & GLM_GENOTYPIC)) {
	  logprint("Error: --hethom must be used with --genotypic.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --hethom flag deprecated.  Use e.g. '--linear hethom' (and\n'--condition-list [filename] recessive' to change covariate coding).\n");
	glm_modifier |= GLM_HETHOM | GLM_CONDITION_RECESSIVE;
	glm_modifier -= GLM_GENOTYPIC;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ide-covar", 10)) {
	logprint("Note: --hide-covar flag deprecated.  Use e.g. '--linear hide-covar'.\n");
	glm_modifier |= GLM_HIDE_COVAR;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ard-call-threshold", 19)) {
        if (!(load_params & 0x180)) {
	  sprintf(logbuf, "Error: --hard-call-threshold must be used with --data or --gen.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!strcmp(argv[cur_arg + 1], "random")) {
	  hard_call_threshold = -1;
	} else {
	  if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.0) || (dxx > (0.5 - SMALL_EPSILON))) {
	    sprintf(logbuf, "Error: Invalid --hard-call-threshold parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  hard_call_threshold = dxx * (1 + SMALL_EPSILON);
	}
      } else if (!memcmp(argptr2, "omog", 5)) {
        calculation_type |= CALC_HOMOG;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'i':
      if (!memcmp(argptr2, "bc", 3)) {
	calculation_type |= CALC_IBC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ndep-pairwise", 14)) {
	if (calculation_type & CALC_LD_PRUNE) {
	  logprint("Error: --indep-pairwise cannot be used with --indep.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 3, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || ((ii == 1) && (param_ct == 3))) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_info.prune_window_size = ii;
	if (param_ct == 4) {
	  if (!match_upper(argv[cur_arg + 2], "KB")) {
	    sprintf(logbuf, "Error: Invalid --indep-pairwise parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ld_info.prune_window_kb = 1;
	} else {
	  jj = strlen(argv[cur_arg + 1]);
	  if ((jj > 2) && match_upper(&(argv[cur_arg + 1][jj - 2]), "KB")) {
	    ld_info.prune_window_kb = 1;
	  }
	}
	ii = atoi(argv[cur_arg + param_ct - 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid increment '%s' for --indep-pairwise.%s", argv[cur_arg + param_ct - 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_info.prune_window_incr = ii;
	if (scan_double(argv[cur_arg + param_ct], &ld_info.prune_last_param) || (ld_info.prune_last_param < 0.0) || (ld_info.prune_last_param >= 1.0)) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise r^2 threshold '%s'.%s", argv[cur_arg + param_ct], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_LD_PRUNE;
        ld_info.modifier |= LD_PRUNE_PAIRWISE;
      } else if (!memcmp(argptr2, "ndep", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 3, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || ((ii == 1) && (param_ct == 3))) {
	  sprintf(logbuf, "Error: Invalid --indep window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_info.prune_window_size = ii;
	if (param_ct == 4) {
	  if (!match_upper(argv[cur_arg + 2], "KB")) {
	    sprintf(logbuf, "Error: Invalid --indep parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ld_info.prune_window_kb = 1;
	} else {
	  jj = strlen(argv[cur_arg + 1]);
	  if ((jj > 2) && match_upper(&(argv[cur_arg + 1][jj - 2]), "KB")) {
	    ld_info.prune_window_kb = 1;
	  }
	}
	ii = atoi(argv[cur_arg + param_ct - 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid increment '%s' for --indep.%s", argv[cur_arg + param_ct - 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_info.prune_window_incr = ii;
	if (scan_double(argv[cur_arg + param_ct], &ld_info.prune_last_param)) {
	  sprintf(logbuf, "Error: Invalid --indep VIF threshold '%s'.%s", argv[cur_arg + param_ct], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ld_info.prune_last_param < 1.0) {
	  sprintf(logbuf, "Error: --indep VIF threshold '%s' too small (must be >= 1).%s", argv[cur_arg + param_ct], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_LD_PRUNE;
      } else if (!memcmp(argptr2, "ndiv-sort", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = (argv[cur_arg + 1][1] == '\0');
        if ((!strcmp(argv[cur_arg + 1], "none")) || ((argv[cur_arg + 1][0] == '0') && jj)) {
	  indiv_sort = INDIV_SORT_NONE;
	} else if ((!strcmp(argv[cur_arg + 1], "natural")) || ((tolower(argv[cur_arg + 1][0]) == 'n') && jj)) {
	  indiv_sort = INDIV_SORT_NATURAL;
	} else if ((!strcmp(argv[cur_arg + 1], "ascii")) || ((tolower(argv[cur_arg + 1][0]) == 'a') && jj)) {
	  indiv_sort = INDIV_SORT_ASCII;
	} else {
	  sprintf(logbuf, "Error: '%s' is not a valid mode for --indiv-sort.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "bs-test", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --ibs-test permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (ii < MAX_THREADS * 2) {
	    sprintf(logbuf, "Error: --ibs test permutation count '%s' too small (min %u).%s", argv[cur_arg + 1], MAX_THREADS * 2, errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ibs_test_perms = ii;
	}
	calculation_type |= CALC_IBS_TEST;
      } else if (!memcmp(argptr2, "id", 3)) {
        logprint("Note: --iid flag deprecated.  Use '--recode vcf-iid'.\n");
	recode_modifier |= RECODE_IID;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "bm", 3)) {
        if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --ibm must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ibm parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(logbuf, "Error: --ibm threshold must be in (0, 1].%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.min_ibm = dxx;
      } else if (!memcmp(argptr2, "mpossible", 10)) {
	logprint("Error: --impossible flag retired.  Use '--genome nudge', or explicitly validate\nZ0/Z1/Z2/PI_HAT in your script.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "nteraction", 11)) {
	logprint("Note: --interaction flag deprecated.  Use e.g. '--linear interaction'.\n");
	glm_modifier |= GLM_INTERACTION;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "bs-matrix", 10)) {
        calculation_type |= CALC_PLINK1_IBS_MATRIX;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "d-delim", 8)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
          id_delim = extract_char_param(argv[cur_arg + 1]);
	  if (!id_delim) {
	    sprintf(logbuf, "Error: --id-delim parameter must be a single character.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (((unsigned char)id_delim) <= ' ') {
	    logprint("Error: --id-delim parameter must be a nonspace character.\n");
            goto main_ret_INVALID_CMDLINE;
	  }
	} else {
          id_delim = '_';
	}
      } else if (!memcmp(argptr2, "nter-chr", 9)) {
        logprint("Note: --inter-chr flag deprecated.  Use e.g. '--r2 inter-chr'.\n");
	ld_info.modifier |= LD_INTER_CHR;
      } else if (!memcmp(argptr2, "nd-major", 9)) {
	logprint("Error: --ind-major is retired, to discourage creation of .bed files that\nconstantly have to be transposed back.  --recode exports individual-major files\nwhich are good enough for smaller jobs; we suggest transposing small data\nwindows on the fly when tackling large jobs.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "mpute-sex", 10)) {
	if (calculation_type & CALC_SEXCHECK) {
	  logprint("Error: --check-sex is redundant with --impute-sex.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
          if (scan_double(argv[cur_arg + 1], &check_sex_fthresh) || (check_sex_fthresh <= 0.0) || ((param_ct == 1) && (check_sex_fthresh > check_sex_mthresh))) {
	    sprintf(logbuf, "Error: Invalid --impute-sex female F-statistic ceiling '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    if (scan_double(argv[cur_arg + 2], &check_sex_mthresh) || (check_sex_mthresh < check_sex_fthresh) || (check_sex_mthresh >= 1.0)) {
	      sprintf(logbuf, "Error: Invalid --impute-sex male F-statistic floor '%s'.%s", argv[cur_arg + 1], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	} else{
	  logprint("Warning: --impute-sex will use default 0.2 and 0.8 thresholds.  This is not\nrecommended.\n");
	}
        calculation_type |= CALC_SEXCHECK;
        misc_flags |= MISC_IMPUTE_SEX;
	sex_missing_pheno |= ALLOW_NO_SEX;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'k':
      if (!memcmp(argptr2, "eep", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&keepname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-fam", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&keepfamname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-allele-order", 17)) {
	misc_flags |= MISC_KEEP_ALLELE_ORDER;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-before-remove", 18)) {
        logprint("Note: --keep-before-remove has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-autoconv", 13)) {
        misc_flags |= MISC_KEEP_AUTOCONV;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-clusters", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&(cluster.keep_fname), argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-cluster-names", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_and_flatten(&(cluster.keep_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'l':
      if (!memcmp(argptr2, "file", 5)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 6) {
	    logprint("Error: --lfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, argv[cur_arg + 1]);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
	load_rare = LOAD_RARE_LGEN;
      } else if (!memcmp(argptr2, "oop-assoc", 10)) {
	if (pheno_modifier & PHENO_ALL) {
	  sprintf(logbuf, "Error: --loop-assoc cannot be used with --all-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if ((strlen(argv[cur_arg + 1]) == 7) && (!memcmp(argv[cur_arg + 1], "keep-", 5)) && match_upper(&(argv[cur_arg + 1][5]), "NA")) {
	    uii = 2;
	  } else if ((strlen(argv[cur_arg + 2]) != 7) || memcmp(argv[cur_arg + 2], "keep-", 5) || (!match_upper(&(argv[cur_arg + 2][5]), "NA"))) {
            sprintf(logbuf, "Error: Invalid --loop-assoc parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_LOAD_CLUSTER_KEEP_NA;
	}
	retval = alloc_fname(&loop_assoc_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "og10", 5)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --log10 must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --log10 flag deprecated.  Use '--adjust log10'.\n");
	mtest_adjust |= ADJUST_LOG10;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ambda", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --lambda must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &adjust_lambda)) {
	  sprintf(logbuf, "Error: Invalid --lambda parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (adjust_lambda < 1) {
	  logprint("Note: --lambda parameter set to 1.\n");
	  adjust_lambda = 1;
	}
	mtest_adjust |= ADJUST_LAMBDA;
      } else if (!memcmp(argptr2, "ist-23-indels", 14)) {
        calculation_type |= CALC_LIST_23_INDELS;
      } else if ((!memcmp(argptr2, "inear", 6)) || (!memcmp(argptr2, "ogistic", 8))) {
#ifndef NOLAPACK
        if (calculation_type & CALC_GLM) {
	  logprint("Error: --logistic cannot be used with --linear.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
	if (*argptr2 == 'o') {
	  glm_modifier |= GLM_LOGISTIC;
#ifdef NOLAPACK
	} else {
	  logprint("Error: --linear requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	  goto main_ret_INVALID_CMDLINE;
#endif
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 10)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (glm_modifier & GLM_MPERM) {
	      sprintf(logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            glm_modifier |= GLM_PERM;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
            if (glm_modifier & GLM_PERM) {
	      sprintf(logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (glm_modifier & GLM_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --%s 'mperm' modifier.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ii = atoi(&(argv[cur_arg + uii][6]));
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --%s mperm parameter '%s'.%s", argptr, &(argv[cur_arg + uii][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            glm_mperm_val = (uint32_t)ii;
            glm_modifier |= GLM_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    glm_modifier |= GLM_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    glm_modifier |= GLM_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "genotypic")) {
	    if (glm_modifier & (GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_GENOTYPIC;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "hethom")) {
	    if (glm_modifier & (GLM_GENOTYPIC | GLM_DOMINANT | GLM_RECESSIVE)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_HETHOM;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "dominant")) {
	    if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_RECESSIVE)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_DOMINANT;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "recessive")) {
	    if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_RECESSIVE;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "no-snp")) {
	    if (mtest_adjust) {
	      sprintf(logbuf, "Error: --%s no-snp cannot be used with --adjust.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    // defer the rest of the check
	    glm_modifier |= GLM_NO_SNP;
	  } else if (!strcmp(argv[cur_arg + uii], "hide-covar")) {
	    glm_modifier |= GLM_HIDE_COVAR;
	  } else if (!strcmp(argv[cur_arg + uii], "sex")) {
	    if (glm_modifier & GLM_NO_X_SEX) {
	      sprintf(logbuf, "Error: --%s 'sex' and 'no-x-sex' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "no-x-sex")) {
	    if (glm_modifier & GLM_SEX) {
	      sprintf(logbuf, "Error: --%s 'sex' and 'no-x-sex' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_NO_X_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "interaction")) {
	    glm_modifier |= GLM_INTERACTION;
	  } else if (!strcmp(argv[cur_arg + uii], "standard-beta")) {
	    if (glm_modifier & GLM_LOGISTIC) {
	      sprintf(logbuf, "Error: --logistic does not have a 'standard-beta' modifier.  (Did you mean\n--linear or 'beta'?)%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_STANDARD_BETA;
	  } else if (!strcmp(argv[cur_arg + uii], "beta")) {
	    glm_modifier |= GLM_BETA;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    glm_modifier |= GLM_SET_TEST;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    sprintf(logbuf, "Error: Improper --%s mperm syntax.  (Use '--%s mperm=[value]'.)\n", argptr, argptr);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", argptr, argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if ((glm_modifier & GLM_NO_SNP) && (glm_modifier & GLM_NO_SNP_EXCL)) {
	  sprintf(logbuf, "Error: --%s 'no-snp' modifier conflicts with another modifier.%s", argptr, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_GLM;
      } else if (!memcmp(argptr2, "d-xchr", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '3') || (argv[cur_arg + 1][1] != '\0')) {
	  sprintf(logbuf, "Error: Invalid --ld-xchr parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (cc == '2') {
          ld_info.modifier |= LD_IGNORE_X;
	} else if (cc == '3') {
	  ld_info.modifier |= LD_WEIGHTED_X;
	}
      } else if (!memcmp(argptr2, "asso", 5)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &lasso_h2) || (lasso_h2 > 1) || (lasso_h2 <= 0)) {
	  sprintf(logbuf, "Error: Invalid --lasso heritability estimate '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "report-zeroes")) {
	    misc_flags |= MISC_LASSO_REPORT_ZEROES;
#ifndef STABLE_BUILD
	  } else if (!strcmp(argv[cur_arg + uii], "no-geno-std")) {
	    misc_flags |= MISC_LASSO_NO_GENO_STD;
#endif
	  } else if (lasso_minlambda > 0) {
            sprintf(logbuf, "Error: Invalid --lasso parameter sequence.%s", errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  } else if (scan_double(argv[cur_arg + uii], &lasso_minlambda) || (lasso_minlambda <= 0)) {
	    sprintf(logbuf, "Error: Invalid --lasso minimum lambda '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        calculation_type |= CALC_LASSO;
      } else if (!memcmp(argptr2, "asso-select-covars", 19)) {
	if (!(calculation_type & CALC_LASSO)) {
	  logprint("Error: --lasso-select-covars must be used with --lasso.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!covar_fname) {
	  logprint("Error: --lasso-select-covars must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        misc_flags |= MISC_LASSO_SELECT_COVARS;
	if (param_ct) {
	  retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &lasso_select_covars_range_list, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "d-window", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ii = atoi(argv[cur_arg + 1]);
        if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --ld-window window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ld_info.window_size = ii;
      } else if (!memcmp(argptr2, "d-window-kb", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(logbuf, "Error: Invalid --ld-window-kb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (dxx > 2147483.647) {
	  ld_info.window_bp = 2147483647;
	} else {
	  ld_info.window_bp = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "d-window-r2", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0) || (dxx > 1)) {
	  sprintf(logbuf, "Error: Invalid --ld-window-r2 parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ld_info.window_r2 = dxx;
      } else if (!memcmp(argptr2, "d-snp", 6)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&(ld_info.snpstr), argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "d-snps", 7)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &(ld_info.snps_rl), 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "d-snp-list", 11)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&(ld_info.snpstr), argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	ld_info.modifier |= LD_SNP_LIST_FILE;
      } else if (!memcmp(argptr2, "d", 2)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (alloc_string(&(epi_info.ld_mkr1), argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	if (alloc_string(&(epi_info.ld_mkr2), argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
	if (param_ct == 3) {
	  if (strcmp(argv[cur_arg + 3], "hwe-midp")) {
	    sprintf(logbuf, "Error: Invalid --ld parameter '%s'.%s", argv[cur_arg + 3], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
          epi_info.modifier |= EPI_HWE_MIDP;
	}
        calculation_type |= CALC_EPI;
      } else if ((!memcmp(argptr2, "ookup", 6)) ||
                 (!memcmp(argptr2, "ookup-list", 11)) ||
                 (!memcmp(argptr2, "ookup-gene", 11)) ||
                 (!memcmp(argptr2, "ookup-gene-kb", 14)) ||
                 (!memcmp(argptr2, "ookup-gene-list", 16))) {
        logprint("Error: --lookup commands have been retired since the Sullivan Lab web database\nis no longer operational.  Use e.g. PLINK/SEQ's lookup command instead.\n");
        goto main_ret_INVALID_CMDLINE;
      } else {
        goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'm':
      if (!memcmp(argptr2, "ap", 3)) {
	if ((load_params & 0x3fc) || (load_rare & (~(LOAD_RARE_CNV | LOAD_RARE_GVAR)))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 4;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --map parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "issing-genotype", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((argv[cur_arg + 1][1] != '\0') || (((unsigned char)cc) <= ' ') || ((cc > '0') && (cc <= '4')) || (cc == 'A') || (cc == 'C') || (cc == 'G') || (cc == 'T')) {
	  sprintf(logbuf, "Error: Invalid --missing-genotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	g_missing_geno_ptr = &(g_one_char_strs[((unsigned char)cc) * 2]);
	g_output_missing_geno_ptr = g_missing_geno_ptr;
      } else if (!memcmp(argptr2, "issing-phenotype", 17)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	missing_pheno = atoi(argv[cur_arg + 1]);
	if ((missing_pheno == 0) || (missing_pheno == 1)) {
	  sprintf(logbuf, "Error: Invalid --missing-phenotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if ((!memcmp(argptr2, "issing-code", 12))) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  missing_code = argv[cur_arg + 1];
	} else {
	  missing_code = (char*)"";
	}
      } else if (!memcmp(argptr2, "ake-pheno", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (alloc_string(&makepheno_str, argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
	if (((argv[cur_arg + 2][0] == '\'') || (argv[cur_arg + 2][0] == '"')) && (argv[cur_arg + 2][1] == '*') && (argv[cur_arg + 2][2] == argv[cur_arg + 2][0]) && (!argv[cur_arg + 2][3])) {
	  memcpy(makepheno_str, "*", 2);
	}
      } else if (!memcmp(argptr2, "pheno", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mpheno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mpheno_col = ii;
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mfilter parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mfilter_col = ii;
      } else if (!memcmp(argptr2, "emory", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	malloc_size_mb = atoi(argv[cur_arg + 1]);
	if (malloc_size_mb < WKSPACE_MIN_MB) {
	  if (malloc_size_mb > 0) {
	    sprintf(logbuf, "Error: Invalid --memory parameter '%s' (minimum %u).%s", argv[cur_arg + 1], WKSPACE_MIN_MB, errstr_append);
	  } else {
	    sprintf(logbuf, "Error: Invalid --memory parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  }
	  goto main_ret_INVALID_CMDLINE_3;
	}
#ifndef __LP64__
	if (malloc_size_mb > 2944) {
	  logprint("Error: --memory parameter too large for 32-bit version (max 2944).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
      } else if (!memcmp(argptr2, "af", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &min_maf)) {
	    sprintf(logbuf, "Error: Invalid --maf parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (min_maf <= 0.0) {
	    sprintf(logbuf, "Error: --maf parameter '%s' too small (must be > 0).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (min_maf > max_maf) {
	    sprintf(logbuf, "Error: --maf parameter '%s' too large (must be <= %g).%s", argv[cur_arg + 1], max_maf, errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  min_maf = 0.01;
	}
      } else if (!memcmp(argptr2, "ax-maf", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &max_maf)) {
	  sprintf(logbuf, "Error: Invalid --max-maf parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (max_maf < min_maf) {
	  sprintf(logbuf, "Error: --max-maf parameter '%s' too small (must be >= %g).%s", argv[cur_arg + 1], min_maf, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (max_maf >= 0.5) {
	  sprintf(logbuf, "Error: --max-maf parameter '%s' too large (must be < 0.5).%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ind", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &mind_thresh)) {
	    sprintf(logbuf, "Error: Invalid --mind parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((mind_thresh < 0.0) || (mind_thresh > 1.0)) {
	    sprintf(logbuf, "Error: Invalid --mind parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  mind_thresh = 0.1;
	}
      } else if (!memcmp(argptr2, "ake-grm", 8)) {
	logprint("Error: --make-grm has been retired due to inconsistent meaning across GCTA\nversions.  Use --make-grm-gz or --make-grm-bin.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ake-grm-gz", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rel_calc_type |= REL_CALC_GZ | REL_CALC_GRM;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-grm-gz 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm-gz 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + uii], "no-gz")) {
	    rel_calc_type &= ~REL_CALC_GZ;
	  } else if ((!strcmp(argv[cur_arg + uii], "ibc2")) || (!strcmp(argv[cur_arg + uii], "ibc3"))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-grm-gz 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm-gz '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + uii][3] - '0';
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    rel_calc_type |= REL_CALC_SINGLE_PREC;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-grm-gz parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-grm-bin", 12)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  sprintf(logbuf, "Error: --make-grm-bin cannot be used with --make-grm-gz.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rel_calc_type |= REL_CALC_GRM_BIN | REL_CALC_SINGLE_PREC;
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-grm-bin 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if ((!strcmp(argv[cur_arg + 1], "ibc2")) || (!strcmp(argv[cur_arg + 1], "ibc3"))) {
	    ibc_type = argv[cur_arg + 1][3] - '0';
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-grm-bin parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-rel", 8)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  sprintf(logbuf, "Error: --make-rel cannot be used with --make-grm-gz/--make-grm-bin.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (rel_calc_type & REL_CALC_BIN) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (rel_calc_type & REL_CALC_GZ) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_TRI;
	  } else if ((!strcmp(argv[cur_arg + uii], "ibc2")) || (!strcmp(argv[cur_arg + uii], "ibc3"))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + uii][3] - '0';
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    rel_calc_type |= REL_CALC_SINGLE_PREC;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-rel parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(rel_calc_type & REL_CALC_SHAPEMASK)) {
	  rel_calc_type |= (rel_calc_type & REL_CALC_BIN)? REL_CALC_SQ : REL_CALC_TRI;
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "atrix", 6)) {
	logprint("Note: --matrix flag deprecated.  Migrate to '--distance ibs flat-missing',\n'--r2 square', etc.\n");
        matrix_flag_state = 1;
	if (calculation_type & CALC_CLUSTER) {
	  calculation_type |= CALC_PLINK1_IBS_MATRIX;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "af-succ", 8)) {
	misc_flags |= MISC_MAF_SUCC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ap3", 4)) {
	logprint("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-bed", 8)) {
        if (misc_flags & MISC_KEEP_AUTOCONV) {
	  sprintf(logbuf, "Error: --make-bed cannot be used with --keep-autoconv.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  sprintf(logbuf, "Error: --%s doesn't accept parameters.%s%s", argptr, ((param_ct == 1) && (!outname_end))? "  (Did you forget '--out'?)" : "", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_MAKE_BED;
      } else if (!memcmp(argptr2, "erge", 5)) {
	if (calculation_type & CALC_MERGE) {
	  sprintf(logbuf, "Error: --merge cannot be used with --bmerge.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare & LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --merge does not currently support .cnv filesets.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (param_ct == 2) {
	  if (++jj > FNAMESIZE) {
	    logprint("Error: --merge .ped filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --merge .map filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logprint("Error: --merge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(mergename1, argv[cur_arg + 1], jj), ".ped", 5);
	  memcpy(memcpya(mergename2, argv[cur_arg + 1], jj), ".map", 5);
	}
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-list", 10)) {
	if (calculation_type & CALC_MERGE) {
	  logprint("Error: --merge-list cannot be used with --merge or --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]) + 1;
	if (jj > FNAMESIZE) {
	  logprint("Error: --merge-list filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(mergename1, argv[cur_arg + 1], jj);
	merge_type |= MERGE_LIST;
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-mode", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '7') || (argv[cur_arg + 1][1] != '\0')) {
          sprintf(logbuf, "Error: Invalid --merge-mode parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((merge_type & MERGE_LIST) && (cc > '5')) {
	  sprintf(logbuf, "Error: --merge-mode 6-7 cannot be used with --merge-list.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        merge_type |= cc - '0';
      } else if (!memcmp(argptr2, "erge-equal-pos", 15)) {
	merge_type |= MERGE_EQUAL_POS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ust-have-sex", 13)) {
        sex_missing_pheno |= MUST_HAVE_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "covar", 6)) {
        if (!(calculation_type & CALC_GXE)) {
	  logprint("Error: --mcovar must be used with --covar and --gxe.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (gxe_mcovar > 1) {
	  logprint("Error: --mcovar cannot be used with a --gxe parameter.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mcovar parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --mcovar flag deprecated.  Use '--gxe [covariate index]'.\n");
	gxe_mcovar = ii;
      } else if (!memcmp(argptr2, "odel", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_ASSOC_FDEPR) {
	  model_modifier &= ~(MODEL_ASSOC | MODEL_ASSOC_FDEPR);
	} else if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "fisher")) {
	    if (model_modifier & MODEL_TRENDONLY) {
	      sprintf(logbuf, "Error: --model 'fisher' and 'trend-only' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher-midp")) {
	    if (model_modifier & MODEL_TRENDONLY) {
	      sprintf(logbuf, "Error: --model 'fisher-midp' and 'trend-only' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER | MODEL_FISHER_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: --model 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "dom")) {
	    if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PDOM;
	  } else if (!strcmp(argv[cur_arg + uii], "rec")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PREC;
	  } else if (!strcmp(argv[cur_arg + uii], "gen")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (mtest_adjust) {
	      sprintf(logbuf, "Error: --model perm-gen cannot be used with --adjust.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PGEN;
	  } else if (!strcmp(argv[cur_arg + uii], "trend")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PTREND;
	  } else if (!strcmp(argv[cur_arg + uii], "trend-only")) {
	    if (model_modifier & (MODEL_FISHER | MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (model_modifier & MODEL_PERM) {
	      sprintf(logbuf, "Error: --model 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --model 'mperm' modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    kk = atoi(&(argv[cur_arg + uii][6]));
	    if (kk < 1) {
	      sprintf(logbuf, "Error: Invalid --model mperm parameter '%s'.%s", &(argv[cur_arg + uii][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_mperm_val = (uint32_t)kk;
	    model_modifier |= MODEL_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logprint("Error: Improper --model mperm syntax.  (Use '--model mperm=[value]'.)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else if (!strcmp(argv[cur_arg + uii], "set-test")) {
	    model_modifier |= MODEL_SET_TEST;
	  } else {
	    sprintf(logbuf, "Error: Invalid --model parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "odel-dom", 9)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model-dom cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-dom flag deprecated.  Use '--model dom'.\n");
	model_modifier |= MODEL_PDOM;
	calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-gen", 9)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model-gen cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-gen flag deprecated.  Use '--model gen'.\n");
	model_modifier |= MODEL_PGEN;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-rec", 9)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model-rec cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-rec flag deprecated.  Use '--model rec'.\n");
	model_modifier |= MODEL_PREC;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-trend", 11)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model-trend cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PREC)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-trend flag deprecated.  Use '--model trend'.\n");
	model_modifier |= MODEL_PTREND;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (model_modifier & (MODEL_PERM | MODEL_MPERM)) {
	  sprintf(logbuf, "Error: --mperm cannot be used with --%s %sperm.%s", (model_modifier & MODEL_ASSOC)? "assoc" : "model", (model_modifier & MODEL_PERM)? "" : "m", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_PERM | GLM_MPERM)) {
	  sprintf(logbuf, "Error: --mperm cannot be used with --%s %sperm.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_PERM)? "" : "m", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mperm parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mperm_val = (uint32_t)ii;
	if (load_rare & LOAD_RARE_CNV) {
	  if ((cnv_calc_type & CNV_INDIV_PERM) && (!cnv_indiv_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-indiv-perm [perm. count]'.\n");
	    cnv_indiv_mperms = mperm_val;
	  } else if ((cnv_calc_type & CNV_TEST_REGION) && (!cnv_test_region_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-test-region [perm. count]'.\n");
	  } else if ((cnv_calc_type & CNV_ENRICHMENT_TEST) && (!cnv_enrichment_test_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-enrichment-test [perm. count]'.\n");
	  } else {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-test [permutation count]'.\n");
            if (!(cnv_calc_type & (CNV_INDIV_PERM | CNV_ENRICHMENT_TEST | CNV_TEST | CNV_TEST_REGION))) {
	      cnv_calc_type |= CNV_TEST;
	    }
	    cnv_test_mperms = mperm_val;
	  }
	  // if e.g. --cnv-test-region had a valid parameter, don't clobber it
	  if (!cnv_test_region_mperms) {
	    cnv_test_region_mperms = mperm_val;
	  }
	  if (!cnv_enrichment_test_mperms) {
	    cnv_enrichment_test_mperms = mperm_val;
	  }
	} else {
	  logprint("Note: --mperm flag deprecated.  Use e.g. '--model mperm=[value]'.\n");
	  model_mperm_val = mperm_val;
	  model_modifier |= MODEL_MPERM;
	  glm_mperm_val = mperm_val;
	  glm_modifier |= GLM_MPERM;
          testmiss_mperm_val = mperm_val;
          testmiss_modifier |= TESTMISS_MPERM;
	}
      } else if (!memcmp(argptr2, "perm-save", 10)) {
	if (glm_modifier & GLM_NO_SNP) {
          sprintf(logbuf, "Error: --mperm-save cannot be used with --linear/--logistic no-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mperm_save |= MPERM_DUMP_BEST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "perm-save-all", 14)) {
	mperm_save |= MPERM_DUMP_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "c", 2)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mc must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --mc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.max_size = ii;
      } else if (!memcmp(argptr2, "cc", 2)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mcc must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mcc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (((uint32_t)ii) > cluster.max_size) {
          logprint("Error: --mcc parameter exceeds --mc parameter.\n");
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.max_cases = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mcc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (((uint32_t)ii) > cluster.max_size) {
          logprint("Error: --mcc parameter exceeds --mc parameter.\n");
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.max_ctrls = ii;
      } else if (!memcmp(argptr2, "atch", 5)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --match must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cluster.match_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
	  if (alloc_string(&cluster.match_missing_str, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "atch-type", 10)) {
	if (!cluster.match_fname) {
	  sprintf(logbuf, "Error: --match-type must be used with --match.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cluster.match_type_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ds-plot", 8)) {
#ifdef NOLAPACK
	// PLINK 1.07's SVD-based non-LAPACK implementation does not conform to
	// classical MDS, so we do not replicate it.
        logprint("Error: --mds-plot requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mds-plot must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cluster.mds_dim_ct = 0;
        for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "by-cluster")) {
	    cluster.modifier |= CLUSTER_MDS;
	  } else if (!strcmp(argv[cur_arg + uii], "eigvals")) {
	    cluster.modifier |= CLUSTER_MDS_EIGVALS;
	  } else {
	    ii = atoi(argv[cur_arg + uii]);
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --mds-plot parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    } else if (cluster.mds_dim_ct) {
	      sprintf(logbuf, "Error: Invalid --mds-plot parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    cluster.mds_dim_ct = ii;
	  }
	}
#endif
      } else if (!memcmp(argptr2, "ds-cluster", 11)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mds-cluster must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --mds-cluster flag deprecated.  Use '--mds-plot by-cluster'.\n");
        cluster.modifier |= CLUSTER_MDS;
      } else if (!memcmp(argptr2, "within", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mwithin parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        mwithin_col = ii;
      } else if (!memcmp(argptr2, "in", 3)) {
        if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --min must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --min parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < -1.0) || (dxx > 1.0)) {
          sprintf(logbuf, "Error: --min threshold must be between -1 and 1 inclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (dxx > genome_max_pi_hat) {
	  logprint("Error: --min value cannot be greater than --max value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        genome_modifier |= GENOME_FILTER_PI_HAT;
	genome_min_pi_hat = dxx;
      } else if (!memcmp(argptr2, "ax", 3)) {
        if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --max must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --max parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < -1.0) || (dxx > 1.0)) {
          sprintf(logbuf, "Error: --max threshold must be between -1 and 1 inclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	genome_modifier |= GENOME_FILTER_PI_HAT;
	genome_max_pi_hat = dxx;
      } else if (!memcmp(argptr2, "ake-founders", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "require-2-missing")) {
	    misc_flags |= MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING;
	  } else if (!strcmp(argv[cur_arg + uii], "first")) {
	    misc_flags |= MISC_MAKE_FOUNDERS_FIRST;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-founders parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        misc_flags |= MISC_MAKE_FOUNDERS;
      } else if (!memcmp(argptr2, "issing", 7)) {
	calculation_type |= CALC_MISSING_REPORT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "h", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "bd")) {
	    sprintf(logbuf, "Error: Invalid --mh parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  misc_flags |= MISC_CMH_BD;
	}
	calculation_type |= CALC_CMH;
	logprint("Error: --mh is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "h2", 3)) {
	if (calculation_type & CALC_CMH) {
	  logprint("Error: --mh2 cannot be used with --mh.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	calculation_type |= CALC_CMH;
        misc_flags |= MISC_CMH2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-set", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&set_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        set_info.modifier |= SET_MAKE_FROM_RANGES;
      } else if (!memcmp(argptr2, "ake-set-border", 15)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --make-set-border must be used with --make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0)) {
	  sprintf(logbuf, "Error: Invalid --make-set-border parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (dxx > 2147483.647) {
	  set_info.make_set_border = 2147483647;
	} else {
	  set_info.make_set_border = ((int32_t)(dxx * 1000 * (1 + SMALL_EPSILON)));
	}
      } else if (!memcmp(argptr2, "ake-set-collapse-group", 23)) {
        if (!set_info.fname) {
	  sprintf(logbuf, "Error: --make-set-collapse-group must be used with --make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        set_info.modifier |= SET_MAKE_COLLAPSE_GROUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-set-complement-all", 23)) {
	if (set_info.modifier & SET_COMPLEMENTS) {
	  sprintf(logbuf, "Error: --make-set-complement-all cannot be used with --complement-sets.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&set_info.merged_set_name, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
	set_info.modifier |= SET_COMPLEMENTS;
      } else if (!memcmp(argptr2, "ake-set-complement-group", 25)) {
        if (!set_info.fname) {
	  sprintf(logbuf, "Error: --make-set-complement-group must be used with --make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (set_info.modifier & (SET_COMPLEMENTS | SET_MAKE_COLLAPSE_GROUP)) {
	  sprintf(logbuf, "Error: --make-set-complement-group cannot be used with --complement-sets,\n--make-set-collapse-group, or --make-set-complement-all.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        set_info.modifier |= SET_COMPLEMENTS | SET_C_PREFIX | SET_MAKE_COLLAPSE_GROUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "erge-x", 7)) {
	if ((chrom_info.x_code == -1) || (chrom_info.xy_code == -1)) {
	  sprintf(logbuf, "Error: --merge-x must be used with a chromosome set containing X and XY codes.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	misc_flags |= MISC_MERGEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lma", 4)) {
        logprint("Error: --mlma is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "lma-loco", 9)) {
        logprint("Error: --mlma-loco is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "lma-no-adj-covar", 17)) {
        logprint("Error: --mlma-no-adj-covar is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ishap-window", 13)) {
        logprint("Error: --mishap-window is provisionally retired.  Contact the developers if you\nneed this function.\n");
        goto main_ret_INVALID_CMDLINE;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'n':
      if (!memcmp(argptr2, "o-fid", 6)) {
	fam_cols &= ~FAM_COL_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-parents", 10)) {
	fam_cols &= ~FAM_COL_34;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-sex", 6)) {
	if (filter_binary & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
	  logprint("Error: --filter-males/--filter-females cannot be used with --no-sex.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	fam_cols &= ~FAM_COL_5;
	sex_missing_pheno |= ALLOW_NO_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-pheno", 8)) {
	fam_cols &= ~FAM_COL_6;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "onfounders", 11)) {
	misc_flags |= MISC_NONFOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eighbour", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --neighbour parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	neighbor_n1 = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --neighbour parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	neighbor_n2 = ii;
	if (neighbor_n2 < neighbor_n1) {
	  sprintf(logbuf, "Error: Second --neighbour parameter cannot be smaller than first parameter.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_NEIGHBOR;
      } else if (!memcmp(argptr2, "ot-chr", 7)) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --from cannot be used with --autosome{-xy} or --{not-}chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	// allowed:
	//   --allow-extra-chr --chr 5-22 bobs_chrom --not-chr 17
	// allowed:
	//   --allow-extra-chr --not-chr 12-17 bobs_chrom
	// does not make sense, disallowed:
	//   --allow-extra-chr --chr 5-22 --not-chr bobs_chrom

	// --allow-extra-chr present, --chr/--autosome{-xy} not present
	uii = ((misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1) && (!chrom_info.is_include_stack);
	retval = parse_chrom_ranges(param_ct, '-', &(argv[cur_arg]), chrom_exclude, &chrom_info, uii, argptr);
	if (retval) {
	  goto main_ret_1;
	}
	if (chrom_info.is_include_stack) {
	  fill_chrom_mask(&chrom_info);
	}
	for (uii = 0; uii < CHROM_MASK_INITIAL_WORDS; uii++) {
	  chrom_info.chrom_mask[uii] &= ~chrom_exclude[uii];
	}
	if (all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS) && ((!((misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1)) || (chrom_info.is_include_stack && (!chrom_info.incl_excl_name_stack)))) {
	  sprintf(logbuf, "Error: All chromosomes excluded.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	chrom_flag_present = 1;
      } else if (!memcmp(argptr2, "udge", 5)) {
        if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --nudge must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --nudge flag deprecated.  Use '--genome nudge'.\n");
        genome_modifier |= GENOME_NUDGE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "o-snp", 6)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --no-snp must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (mtest_adjust) {
	  sprintf(logbuf, "Error: --no-snp cannot be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (mperm_save & MPERM_DUMP_BEST) {
	  sprintf(logbuf, "Error: --no-snp cannot be used with --mperm-save.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((glm_modifier & (GLM_NO_SNP_EXCL - GLM_HETHOM - GLM_DOMINANT)) || ((glm_modifier & (GLM_HETHOM | GLM_DOMINANT)) && (!(glm_modifier & (GLM_CONDITION_DOMINANT | GLM_CONDITION_RECESSIVE))))) {
	  sprintf(logbuf, "Error: --no-snp conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --no-snp flag deprecated.  Use e.g. '--linear no-snp'.\n");
        glm_modifier |= GLM_NO_SNP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-x-sex", 8)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --no-x-sex must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_NO_SNP | GLM_SEX)) {
	  sprintf(logbuf, "Error: --no-x-sex conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --no-x-sex flag deprecated.  Use e.g. '--linear no-x-sex'.\n");
	glm_modifier |= GLM_NO_X_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "op", 3)) {
	if (!(epi_info.modifier & EPI_FAST)) {
	  sprintf(logbuf, "Error: --nop must be used with --fast-epistasis.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --nop flag deprecated.  Use '--fast-epistasis nop'.\n");
        epi_info.modifier |= EPI_FAST_NO_P_VALUE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "oweb", 5)) {
        logprint("Note: --noweb has no effect since no web check is implemented yet.\n");
	goto main_param_zero;
      } else {
        goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'o':
      if (!memcmp(argptr2, "utput-missing-genotype", 23)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((argv[cur_arg + 1][1] != '\0') || (((unsigned char)cc) <= ' ')) {
	  sprintf(logbuf, "Error: Invalid --output-missing-genotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	g_output_missing_geno_ptr = &(g_one_char_strs[((unsigned char)cc) * 2]);
      } else if (!memcmp(argptr2, "utput-missing-phenotype", 24)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > 31) {
	  logprint("Error: --output-missing-phenotype string too long (max 31 chars).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  logprint("Error: --output-missing-phenotype parameter currently must be numeric.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	memcpy(output_missing_pheno, argv[cur_arg + 1], jj + 1);
      } else if (!memcmp(argptr2, "blig-clusters", 14)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&oblig_missing_info.indiv_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	logprint("Note: --oblig-clusters flag deprecated.  Use just --oblig-missing.\n");
      } else if (!memcmp(argptr2, "blig-missing", 13)) {
	if ((geno_thresh == 1.0) && (mind_thresh == 1.0) && (!(calculation_type & CALC_MISSING_REPORT))) {
	  sprintf(logbuf, "Error: --oblig-missing must be used with --geno, --mind, or --missing.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!oblig_missing_info.indiv_fname) {
          if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else if (param_ct != 1) {
          sprintf(logbuf, "Error: --oblig-missing requires exactly one parameter when --oblig-clusters is\nalso present.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&oblig_missing_info.marker_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (param_ct == 2) {
	  retval = alloc_fname(&oblig_missing_info.indiv_fname, argv[cur_arg + 2], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (memcmp(argptr2, "ut", 3)) {
	// --out is a special case due to logging
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'p':
      if (!memcmp(argptr2, "ed", 3)) {
	if ((load_params & 0x3fa) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 2;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --ped parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "heno", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (makepheno_str) {
	  logprint("Error: --pheno and --make-pheno flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "heno-name", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (mpheno_col != 0) {
	  logprint("Error: --mpheno and --pheno-name flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (alloc_string(&phenoname_str, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "heno-merge", 11)) {
	pheno_modifier |= PHENO_MERGE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rune", 5)) {
	misc_flags |= MISC_PRUNE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "arallel", 8)) {
	if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	  sprintf(logbuf, "Error: --parallel cannot be used with '--distance square'.  Use '--distance\nsquare0' or plain --distance instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((dist_calc_type & DISTANCE_BIN) && (!(dist_calc_type & DISTANCE_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --parallel cannot be used with plain '--distance bin'.  Use '--distance\nbin square0' or '--distance bin triangle' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	  sprintf(logbuf, "Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain '--make-rel' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((rel_calc_type & REL_CALC_BIN) && (!(rel_calc_type & REL_CALC_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --parallel cannot be used with plain '--make-rel bin'.  Use '--make-rel\nbin square0' or '--make-rel bin triangle' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_PLINK1_DISTANCE_MATRIX) {
	  sprintf(logbuf, "Error: --parallel and --distance-matrix cannot be used together.  Use\n--distance instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_GROUPDIST) {
	  sprintf(logbuf, "Error: --parallel and --groupdist cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_CLUSTER) {
	  sprintf(logbuf, "Error: --parallel and --cluster cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_NEIGHBOR) {
	  sprintf(logbuf, "Error: --parallel and --neighbour cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || (ii > PARALLEL_MAX)) {
	  sprintf(logbuf, "Error: Invalid --parallel job index '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	parallel_idx = ii - 1; // internal 0..(n-1) indexing
	ii = atoi(argv[cur_arg + 2]);
	if ((ii < 2) || (ii > PARALLEL_MAX) || (((uint32_t)ii) < parallel_idx)) {
	  sprintf(logbuf, "Error: Invalid --parallel total job count '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        parallel_tot = ii;
      } else if (!memcmp(argptr2, "pc-gap", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ppc-gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dxx *= 1000;
	if (dxx < 0) {
	  ppc_gap = 0;
	} else if (dxx > 2147483647) {
	  ppc_gap = 0x7fffffff;
	} else {
	  ppc_gap = (int32_t)(dxx * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "erm", 4)) {
	if ((model_modifier & MODEL_MPERM) && (calculation_type & CALC_MODEL)) {
	  sprintf(logbuf, "Error: --perm cannot be used with --%s mperm.%s", (model_modifier & MODEL_ASSOC)? "assoc" : "model", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((calculation_type & CALC_GLM) && (glm_modifier & (GLM_MPERM + GLM_NO_SNP))) {
	  sprintf(logbuf, "Error: --perm cannot be used with --%s %s.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_MPERM)? "mperm" : "no-snp", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (model_modifier & MODEL_MPERM) {
	  sprintf(logbuf, "Error: --perm cannot be used with --mperm.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	model_modifier |= MODEL_PERM;
        glm_modifier |= GLM_PERM;
        testmiss_modifier |= TESTMISS_PERM;
	logprint("Note: --perm flag deprecated.  Use e.g. '--model perm'.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "erm-count", 10)) {
	model_modifier |= MODEL_PERM_COUNT;
	glm_modifier |= GLM_PERM_COUNT;
        testmiss_modifier |= TESTMISS_PERM_COUNT;
	logprint("Note: --perm-count flag deprecated.  Use e.g. '--model perm-count'.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "2", 2)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  logprint("Error: --p2 must be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & MODEL_QMASK) {
	  sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with --p2.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --p2 flag deprecated.  Use '--assoc p2 ...'.\n");
	model_modifier |= MODEL_ASSOC_P2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --pfilter parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --pfilter threshold must be between 0 and 1 exclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	pfilter = dxx;
      } else if (!memcmp(argptr2, "erm-batch-size", 1)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	perm_batch_size = atoi(argv[cur_arg + 1]);
	if ((perm_batch_size < 1) || (perm_batch_size > 0x7fffffff)) {
	  sprintf(logbuf, "Error: Invalid --perm-batch-size parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "pc", 3)) {
	if (!(calculation_type & (CALC_NEIGHBOR | CALC_CLUSTER))) {
          sprintf(logbuf, "Error: --ppc must be used with --cluster or --neigbour.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ppc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --ppc threshold must be between 0 and 1 exclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.ppc = dxx;
      } else if (!memcmp(argptr2, "ool-size", 9)) {
	if (!(homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE))) {
          logprint("Error: --pool-size must be used with --homozyg group{-verbose}.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --pool-size parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	homozyg.pool_size_min = ii;
      } else if (!memcmp(argptr2, "arameters", 10)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --parameters must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &parameters_range_list, 1);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "roxy-assoc", 11)) ||
                 (!memcmp(argptr2, "roxy-drop", 10)) ||
                 (!memcmp(argptr2, "roxy-impute", 12)) ||
                 (!memcmp(argptr2, "roxy-impute-threshold", 22)) ||
                 (!memcmp(argptr2, "roxy-genotypic-concordance", 27)) ||
                 (!memcmp(argptr2, "roxy-show-proxies", 18)) ||
                 (!memcmp(argptr2, "roxy-dosage", 12)) ||
                 (!memcmp(argptr2, "roxy-replace", 13)) ||
                 (!memcmp(argptr2, "roxy-verbose", 13))) {
        logprint("Error: PLINK 1 imputation commands have been retired due to poor accuracy.\n(See Nothnagel M et al. (2009) A comprehensive evaluation of SNP genotype\nimputation.)  We suggest using another tool, such as BEAGLE 4 or IMPUTE2, for\nimputation instead.  ('--recode vcf' and --vcf can be used to exchange data\nwith BEAGLE 4, while '--recode oxford' and --data let you work with IMPUTE2.)\n");
        goto main_ret_INVALID_CMDLINE;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'q':
      if (!memcmp(argptr2, "t-means", 8)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  sprintf(logbuf, "Error: --qt-means must be used with --assoc.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_DMASK) {
	  sprintf(logbuf, "Error: --qt-means does not make sense with a case/control-specific --assoc\nmodifier.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --qt-means flag deprecated.  Use '--assoc qt-means ...'.\n");
	model_modifier |= MODEL_QT_MEANS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "q-plot", 7)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --qq-plot must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --qq-plot flag deprecated.  Use '--adjust qq-plot'.\n");
	mtest_adjust |= ADJUST_QQ;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "match", 6)) {
        if (!(calculation_type & CALC_CLUSTER)) {
          sprintf(logbuf, "Error: --qmatch must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&cluster.qmatch_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
	  if (alloc_string(&cluster.qmatch_missing_str, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "t", 2)) {
        if (!cluster.qmatch_fname) {
          sprintf(logbuf, "Error: --qt must be used with --qmatch.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&cluster.qt_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'r':
      if (!memcmp(argptr2, "emove", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&removename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "emove-fam", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&removefamname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "emove-clusters", 15)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&(cluster.remove_fname), argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "emove-cluster-names", 20)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_and_flatten(&(cluster.remove_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "el-cutoff", 10)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel cannot be used with --rel-cutoff.  (Use a combination of\n--make-rel, --keep/--remove, and a filtering script.)%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &rel_cutoff) || (rel_cutoff <= 0.0) || (rel_cutoff >= 1.0)) {
	    sprintf(logbuf, "Error: Invalid --rel-cutoff parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_REL_CUTOFF;
      } else if (!memcmp(argptr2, "egress-distance", 16)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --regress-distance cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  regress_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((regress_iters < 2) || (regress_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --regress-distance jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    ii = atoi(argv[cur_arg + 2]);
	    if (ii <= 0) {
	      sprintf(logbuf, "Error: Invalid --regress-distance jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    regress_d = ii;
	  }
	}
	calculation_type |= CALC_REGRESS_DISTANCE;
      } else if (!memcmp(argptr2, "egress-rel", 11)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --regress-rel flags cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (rel_calc_type & REL_CALC_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --regress-rel cannot currently be used with a single-precision\nrelationship matrix.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  regress_rel_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((regress_rel_iters < 2) || (regress_rel_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --regress-rel jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    ii = atoi(argv[cur_arg + 2]);
	    if (ii <= 0) {
	      sprintf(logbuf, "Error: Invalid --regress-rel jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    regress_rel_d = ii;
	  }
	}
	calculation_type |= CALC_REGRESS_REL;
      } else if ((!memcmp(argptr2, "egress-pcs", 11)) || (!memcmp(argptr2, "egress-pcs-distance", 20))) {
	logprint("Error: --regress-pcs has been temporarily disabled.  Contact the developers if\nyou need a build with the old implementation unlocked.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
	/*
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 5)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&evecname, argv[cur_arg + 1], argptr, 9);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "normalize-pheno")) {
	    regress_pcs_modifier |= REGRESS_PCS_NORMALIZE_PHENO;
	  } else if (!strcmp(argv[cur_arg + uii], "sex-specific")) {
	    regress_pcs_modifier |= REGRESS_PCS_SEX_SPECIFIC;
	  } else if (!strcmp(argv[cur_arg + uii], "clip")) {
	    regress_pcs_modifier |= REGRESS_PCS_CLIP;
	  } else if ((max_pcs != MAX_PCS_DEFAULT) || (argv[cur_arg + uii][0] < '0') || (argv[cur_arg + uii][0] > '9')) {
	    sprintf(logbuf, "Error: Invalid --regress-pcs parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    ii = atoi(argv[cur_arg + uii]);
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --regress-pcs maximum principal component count '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    max_pcs = ii;
	  }
	}
	calculation_type |= CALC_REGRESS_PCS;
      } else if (!memcmp(argptr2, "egress-pcs-distance", 20)) {
	if (calculation_type & CALC_REGRESS_PCS) {
	  sprintf(logbuf, "Error: --regress-pcs-distance cannot be used with --regress-pcs.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_DISTANCE) {
	  sprintf(logbuf, "Error: --regress-pcs-distance cannot be used with --distance.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 11)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&evecname, argv[cur_arg + 1], argptr, 9);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "normalize-pheno")) {
	    regress_pcs_modifier |= REGRESS_PCS_NORMALIZE_PHENO;
	  } else if (!strcmp(argv[cur_arg + uii], "sex-specific")) {
	    regress_pcs_modifier |= REGRESS_PCS_SEX_SPECIFIC;
	  } else if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (parallel_tot > 1) {
	      sprintf(logbuf, "Error: --parallel cannot be used with '--regress-pcs-distance square'.  Use\nthe square0 or triangle shape instead.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square0' and 'triangle' modifiers can't coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square0' and 'triangle' modifiers can't coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (dist_calc_type & DISTANCE_BIN) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (dist_calc_type & DISTANCE_GZ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "ibs")) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logprint("Error: Duplicate --regress-pcs-distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "1-ibs")) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logprint("Error: Duplicate --regress-pcs-distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "allele-ct")) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logprint("Error: Duplicate --regress-pcs-distance 'allele-ct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!strcmp(argv[cur_arg + uii], "flat-missing")) {
	    dist_calc_type |= DISTANCE_FLAT_MISSING;
	  } else if ((max_pcs != MAX_PCS_DEFAULT) || (argv[cur_arg + uii][0] < '0') || (argv[cur_arg + uii][0] > '9')) {
	    sprintf(logbuf, "Error: Invalid --regress-pcs-distance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    ii = atoi(argv[cur_arg + uii]);
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --regress-pcs-distance maximum PC count '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    max_pcs = ii;
	  }
	}
	if (!(dist_calc_type & DISTANCE_TYPEMASK)) {
	  dist_calc_type |= DISTANCE_ALCT;
	}
	calculation_type |= CALC_REGRESS_PCS_DISTANCE;
	*/
      } else if (!memcmp(argptr2, "ead-freq", 9)) {
	if (calculation_type & CALC_FREQ) {
	  sprintf(logbuf, "Error: --freq and --read-freq flags cannot coexist.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&freqname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "ecode", 6)) || (!memcmp(argptr2, "ecode ", 6))) {
	if (argptr2[5] == ' ') {
	  kk = 1;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4 - kk)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (recode_modifier & (RECODE_A | RECODE_AD)) {
	      sprintf(logbuf, "Error: --recode '12' modifier cannot be used with 'A' or 'AD'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_12;
	  } else if (!strcmp(argv[cur_arg + uii], "compound-genotypes")) {
	    if (recode_modifier & RECODE_STRUCTURE) {
              logprint("Error: --recode 'compound-genotypes' modifier cannot be used with 'structure'.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    recode_modifier |= RECODE_COMPOUND;
	  } else if (!strcmp(argv[cur_arg + uii], "23")) {
	    if (recode_type_set(&recode_modifier, RECODE_23)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
          } else if ((!argv[cur_arg + uii][1]) && (tolower(argv[cur_arg + uii][0]) == 'a')) {
	    if (recode_type_set(&recode_modifier, RECODE_A)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if ((!argv[cur_arg + uii][2]) && match_upper(argv[cur_arg + uii], "AD")) {
	    if (recode_type_set(&recode_modifier, RECODE_AD)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (match_upper(argv[cur_arg + uii], "HV")) {
	    if (!argv[cur_arg + uii][2]) {
	      if (recode_type_set(&recode_modifier, RECODE_HV)) {
	        goto main_ret_INVALID_CMDLINE_3;
	      }
	    } else if (!strcmp(&(argv[cur_arg + uii][2]), "-1chr")) {
	      if (recode_type_set(&recode_modifier, RECODE_HV_1CHR)) {
	        goto main_ret_INVALID_CMDLINE_3;
	      }
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "tab")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TAB;
	  } else if (!strcmp(argv[cur_arg + uii], "tabx")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TAB | RECODE_DELIMX;
	  } else if (!strcmp(argv[cur_arg + uii], "spacex")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_DELIMX;
	  } else if (!strcmp(argv[cur_arg + uii], "beagle")) {
	    if (recode_type_set(&recode_modifier, RECODE_BEAGLE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "bimbam")) {
	    if (recode_type_set(&recode_modifier, RECODE_BIMBAM)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "bimbam-1chr")) {
	    if (recode_type_set(&recode_modifier, RECODE_BIMBAM_1CHR)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "fastphase")) {
	    if (recode_type_set(&recode_modifier, RECODE_FASTPHASE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "fastphase-1chr")) {
	    if (recode_type_set(&recode_modifier, RECODE_FASTPHASE_1CHR)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "lgen")) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "lgen-ref")) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN_REF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "list")) {
	    if (recode_type_set(&recode_modifier, RECODE_LIST)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "oxford")) {
            if (recode_type_set(&recode_modifier, RECODE_OXFORD)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "rlist")) {
	    if (recode_type_set(&recode_modifier, RECODE_RLIST)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "structure")) {
	    if (recode_modifier & RECODE_COMPOUND) {
              logprint("Error: --recode 'compound-genotypes' modifier cannot be used with 'structure'.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_STRUCTURE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "transpose")) {
	    if (recode_type_set(&recode_modifier, RECODE_TRANSPOSE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "vcf")) {
	    if (recode_modifier & RECODE_VCF) {
	      sprintf(logbuf, "Error: Conflicting or redundant --recode modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "vcf-fid")) {
	    if (recode_modifier & (RECODE_VCF | RECODE_IID)) {
	      sprintf(logbuf, "Error: Conflicting or redundant --recode modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_FID;
	  } else if (!strcmp(argv[cur_arg + uii], "vcf-iid")) {
	    if (recode_modifier & (RECODE_VCF | RECODE_FID)) {
	      sprintf(logbuf, "Error: Conflicting or redundant --recode modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_IID;
	  } else {
	    sprintf(logbuf, "Error: Invalid --recode parameter '%s'.%s%s", argv[cur_arg + uii], ((uii == param_ct) && (!outname_end))? "  (Did you forget '--out'?)" : "", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RECODE;
      } else if (!memcmp(argptr2, "ecode-whap", 11)) {
        logprint("Error: --recode-whap flag retired since WHAP is no longer supported.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ecode-allele", 13)) {
	if (!(recode_modifier & (RECODE_A | RECODE_AD))) {
	  sprintf(logbuf, "Error: --recode-allele must be used with --recode A or --recode AD.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (recode_modifier & RECODE_12) {
	  sprintf(logbuf, "Error: --recode-allele cannot be used with --recode 12.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&recode_allele_name, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eference", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&lgen_reference_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	lgen_modifier |= LGEN_REFERENCE;
      } else if (!memcmp(argptr2, "ead-genome", 11)) {
	if (calculation_type & CALC_GENOME) {
          sprintf(logbuf, "Error: --read-genome cannot be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (!(calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR))) {
          sprintf(logbuf, "Error: --read-genome cannot be used without --cluster or --neighbour.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	} else if ((cluster.ppc == 0.0) && ((calculation_type & (CALC_DISTANCE | CALC_PLINK1_DISTANCE_MATRIX)) || read_dists_fname)) {
          sprintf(logbuf, "Error: --read-genome is pointless with --distance, --distance-matrix, and\n--read-dists unless --ppc is also present.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&read_genome_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ead-genome-list", 19)) {
	logprint("Error: --read-genome-list flag retired.  Use --parallel + Unix cat instead.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ead-genome-minimal", 19)) {
	logprint("Error: --read-genome-minimal flag retired.  Use --genome gz + --read-genome\ninstead.");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ead-dists", 10)) {
	if (calculation_type & (CALC_DISTANCE | CALC_PLINK1_DISTANCE_MATRIX)) {
	  sprintf(logbuf, "Error: --read-dists cannot be used with a distance matrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (cluster.modifier & CLUSTER_MISSING) {
          sprintf(logbuf, "Error: --read-dists cannot be used with '--cluster missing'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&read_dists_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
          retval = alloc_fname(&read_dists_id_fname, argv[cur_arg + 2], argptr, 0);
          if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "el-check", 9)) {
        if (!(calculation_type & CALC_GENOME)) {
          sprintf(logbuf, "Error: --rel-check must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --rel-check flag deprecated.  Use '--genome rel-check'.\n");
        genome_modifier |= GENOME_REL_CHECK;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ecessive", 9)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --recessive must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT)) {
	  sprintf(logbuf, "Error: --recessive conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --recessive flag deprecated.  Use e.g. '--linear recessive' (and\n'--condition-list [filename] recessive' to change covariate coding).\n");
	glm_modifier |= GLM_RECESSIVE | GLM_CONDITION_RECESSIVE;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if ((*argptr2 == '\0') || (!memcmp(argptr2, "2", 2))) {
	if (calculation_type & CALC_LD) {
          logprint("Error: --r and --r2 cannot be used together.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (*argptr2 == '2') {
          ld_info.modifier |= LD_R2;
	} else if (ld_info.window_r2 != 0.2) {
	  logprint("Error: --ld-window-r2 flag cannot be used with --r.\n");
	}
	if (matrix_flag_state) {
	  matrix_flag_state = 2;
	  ld_info.modifier |= LD_MATRIX_SQ | LD_MATRIX_SPACES;
	}
	for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "square")) {
	    if (ld_info.modifier & LD_MATRIX_SHAPEMASK) {
	      logprint("Error: Multiple --r/--r2 shape modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & (LD_INTER_CHR | LD_DPRIME)) {
	    main_r2_matrix_conflict:
              sprintf(logbuf, "Error: --r/--r2 '%s' cannot be used with matrix output.%s", (ld_info.modifier & LD_INTER_CHR)? "inter-chr" : "dprime", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ld_info.modifier |= LD_MATRIX_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if (ld_info.modifier & LD_MATRIX_SHAPEMASK) {
	      logprint("Error: Multiple --r/--r2 shape modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & (LD_INTER_CHR | LD_DPRIME)) {
	      goto main_r2_matrix_conflict;
	    }
	    ld_info.modifier |= LD_MATRIX_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if (ld_info.modifier & LD_MATRIX_SHAPEMASK) {
	      logprint("Error: Multiple --r/--r2 shape modifiers.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & (LD_INTER_CHR | LD_DPRIME)) {
	      goto main_r2_matrix_conflict;
	    }
	    ld_info.modifier |= LD_MATRIX_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "inter-chr")) {
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
            ld_info.modifier |= LD_INTER_CHR;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (ld_info.modifier & LD_MATRIX_BIN) {
	      logprint("Error: --r/--r2 'gz' and 'bin' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_REPORT_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (ld_info.modifier & (LD_INTER_CHR | LD_DPRIME)) {
	      goto main_r2_matrix_conflict;
	    } else if (ld_info.modifier & LD_REPORT_GZ) {
	      logprint("Error: --r/--r2 'gz' and 'bin' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    } else if (ld_info.modifier & LD_MATRIX_SPACES) {
	      logprint("Error: --r/--r2 'bin' and 'spaces' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_MATRIX_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    // yeah, as a practical matter this should probably be the default
            // since there are no long chains of floating point calculations...
	    ld_info.modifier |= LD_SINGLE_PREC;
	  } else if (!strcmp(argv[cur_arg + uii], "spaces")) {
	    if (ld_info.modifier & (LD_INTER_CHR | LD_DPRIME)) {
	      goto main_r2_matrix_conflict;
	    } else if (ld_info.modifier & LD_MATRIX_BIN) {
	      logprint("Error: --r/--r2 'bin' and 'spaces' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    ld_info.modifier |= LD_MATRIX_SPACES;
	  } else if (!strcmp(argv[cur_arg + uii], "dprime")) {
            if (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_MATRIX_BIN | LD_MATRIX_SPACES)) {
	      goto main_r2_matrix_conflict;
	    }
	    ld_info.modifier |= LD_DPRIME;
	  } else if (!strcmp(argv[cur_arg + uii], "with-freqs")) {
	    ld_info.modifier |= LD_WITH_FREQS;
	  } else if (!strcmp(argv[cur_arg + uii], "yes-really")) {
	    ld_info.modifier |= LD_YES_REALLY;
	  } else {
	    sprintf(logbuf, "Error: Invalid --r/--r2 parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if ((ld_info.modifier & (LD_SINGLE_PREC | LD_MATRIX_BIN)) == LD_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --r/--r2 'single-prec' modifier currently must be used with 'bin'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if ((ld_info.modifier & LD_MATRIX_BIN) && (!(ld_info.modifier & LD_MATRIX_SHAPEMASK))) {
          ld_info.modifier |= LD_MATRIX_SQ;
	}
	if ((ld_info.modifier & LD_MATRIX_SPACES) && (!(ld_info.modifier & LD_MATRIX_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --r/--r2 'spaces' modifier must be used with a shape modifier.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	if ((ld_info.modifier & LD_WEIGHTED_X) && (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_INTER_CHR))) {
	  sprintf(logbuf, "Error: --ld-xchr 3 cannot be used with --r/--r2 non-windowed reports.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	} else if (ld_info.modifier & LD_WEIGHTED_X) {
	  logprint("Error: --r/--r2 + --ld-xchr 3 has not been implemented yet.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	calculation_type |= CALC_LD;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 's':
      if (!memcmp(argptr2, "eed", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rseed_ct = param_ct;
	rseeds = (uint32_t*)malloc(param_ct * sizeof(int32_t));
	for (uii = 1; uii <= param_ct; uii++) {
	  if (strtoui32(argv[cur_arg + uii], &(rseeds[uii - 1]))) {
	    sprintf(logbuf, "Error: Invalid --seed parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
      } else if (!memcmp(argptr2, "ample", 6)) {
	if ((load_params & 0x7f) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	} else if (!(load_params & 0x180)) {
	  sprintf(logbuf, "Error: --sample cannot be used without --data or --gen.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	load_params |= 0x200;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --sample parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "np", 3)) {
        if (markername_from) {
	  sprintf(logbuf, "Error: --snp cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (marker_pos_start != -1) {
	  sprintf(logbuf, "Error: --snp cannot be used with --from-bp/-kb/-mb.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((!all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS)) || chrom_info.incl_excl_name_stack) {
	  sprintf(logbuf, "Error: --snp cannot be used with --autosome{-xy} or --{not-}chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
          sprintf(logbuf, "Error: --snp cannot be used with --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
          sprintf(logbuf, "Error: --snp cannot be used with --exclude-snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --snp cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb and\n--to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "nps", 4)) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --snps cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (marker_pos_start != -1) {
	  sprintf(logbuf, "Error: --snps cannot be used with --from-bp/-kb/-mb.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
	  sprintf(logbuf, "Error: --snps cannot be used with --snp or --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
	  sprintf(logbuf, "Error: --snps cannot be used with --exclude-snps.%s", errstr_append);
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --snps cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb and\n--to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	// mise well allow --snps + --autosome/--autosome-xy/--chr/--not-chr
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &snps_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "et-hh-missing", 14)) {
	misc_flags |= MISC_SET_HH_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et", 3)) {
	if (set_info.fname) {
	  sprintf(logbuf, "Error: --set cannot be used with --make-set.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&set_info.fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "et-collapse-all", 16)) {
	if (!set_info.fname) {
	  logprint("Error: --set-collapse-all must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (set_info.modifier & SET_MAKE_COLLAPSE_GROUP) {
	  logprint("Error: --set-collapse-all cannot be used with --make-set-collapse-group or\n--make-set-complement-group.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (set_info.merged_set_name) {
	  logprint("Error: --set-collapse-all cannot be used with --make-set-complement-all.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&set_info.merged_set_name, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "et-names", 9)) {
	if (!set_info.fname) {
	  logprint("Error: --set-names must be used with --set/--make-set.\n");
          goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_and_flatten(&(set_info.setnames_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ubset", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!set_info.fname) {
	  logprint("Error: --subset must be used with --set/--make-set.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&set_info.subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "imulate", 8)) || (!memcmp(argptr2, "imulate-qt", 11))) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&simulate_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (argptr2[7] == '-') {
	  simulate_flags |= SIMULATE_QT;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "tags")) {
	    if (simulate_flags & SIMULATE_HAPS) {
	      sprintf(logbuf, "Error: --%s 'tags' and 'haps' modifiers cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    simulate_flags |= SIMULATE_TAGS;
	  } else if (!strcmp(argv[cur_arg + uii], "haps")) {
	    if (simulate_flags & SIMULATE_TAGS) {
	      sprintf(logbuf, "Error: --%s 'tags' and 'haps' modifiers cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    simulate_flags |= SIMULATE_HAPS;
	  } else if (match_upper(argv[cur_arg + uii], "ACGT")) {
	    if (simulate_flags & (SIMULATE_1234 | SIMULATE_12)) {
	      sprintf(logbuf, "Error: --%s 'acgt' modifier cannot be used with '1234' or '12'.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            simulate_flags |= SIMULATE_ACGT;
	  } else if (!strcmp(argv[cur_arg + uii], "1234")) {
	    if (simulate_flags & (SIMULATE_ACGT | SIMULATE_12)) {
	      sprintf(logbuf, "Error: --%s '1234' modifier cannot be used with 'acgt' or '12'.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            simulate_flags |= SIMULATE_1234;
	  } else if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (simulate_flags & (SIMULATE_ACGT | SIMULATE_1234)) {
	      sprintf(logbuf, "Error: --%s '12' modifier cannot be used with 'acgt' or '1234'.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            simulate_flags |= SIMULATE_12;
	  } else {
	    sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", argptr, argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	load_rare = LOAD_RARE_SIMULATE;
      } else if (!memcmp(argptr2, "imulate-ncases", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare != LOAD_RARE_SIMULATE) {
	  sprintf(logbuf, "Error: --simulate-ncases must be used with --simulate.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --simulate-ncases parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	simulate_cases = ii;
      } else if (!memcmp(argptr2, "imulate-ncontrols", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare != LOAD_RARE_SIMULATE) {
	  sprintf(logbuf, "Error: --simulate-ncontrols must be used with --simulate.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --simulate-ncontrols parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((!ii) && (!simulate_cases)) {
	  logprint("Error: '--simulate-ncases 0' cannot be used with '--simulate-ncontrols 0'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	simulate_controls = ii;
      } else if (!memcmp(argptr2, "imulate-prevalence", 19)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &simulate_prevalence) || (simulate_prevalence < 0) || (simulate_prevalence > 1)) {
	  sprintf(logbuf, "Error: Invalid --simulate-prevalence parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "imulate-label", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&simulate_label, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "imulate-missing", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &simulate_missing) || (simulate_missing < 0) || (simulate_missing > 1)) {
	  sprintf(logbuf, "Error: Invalid --simulate-missing parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "imulate-n", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare == LOAD_RARE_SIMULATE) {
	  sprintf(logbuf, "Error: --simulate-n must be used with --simulate-qt, not --simulate.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --simulate-n parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	simulate_qt_indivs = ii;
      } else if (!memcmp(argptr2, "imulate-haps", 13)) {
	if (simulate_flags & SIMULATE_TAGS) {
	  sprintf(logbuf, "Error: --simulate-tags cannot be used with --simulate-haps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --simulate-haps flag deprecated.  Use e.g. '--simulate haps'.\n");
	simulate_flags |= SIMULATE_HAPS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "imulate-tags", 13)) {
	if (simulate_flags & SIMULATE_HAPS) {
	  sprintf(logbuf, "Error: --simulate-tags cannot be used with --simulate-haps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --simulate-tags flag deprecated.  Use e.g. '--simulate tags'.\n");
	simulate_flags |= SIMULATE_TAGS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ex", 3)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --sex must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & GLM_NO_X_SEX) {
	  sprintf(logbuf, "Error: --sex conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --sex flag deprecated.  Use e.g. '--linear sex'.\n");
	glm_modifier |= GLM_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "tandard-beta", 13)) {
	if ((!(calculation_type & CALC_GLM)) || (glm_modifier & GLM_LOGISTIC)) {
	  sprintf(logbuf, "Error: --standard-beta must be used wtih --linear.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --standard-beta flag deprecated.  Use '--linear standard-beta'.\n");
	glm_modifier |= GLM_STANDARD_BETA;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et-table", 9)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --set-table must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_WRITE_SET;
	set_info.modifier |= SET_WRITE_TABLE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "et-test", 8)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --set-test must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --set-test flag deprecated.  Use e.g. '--assoc perm set-test'.\n");
	if (calculation_type & CALC_MODEL) {
          model_modifier |= MODEL_SET_TEST;
	}
	if (calculation_type & CALC_GLM) {
	  model_modifier |= GLM_SET_TEST;
	}
	if ((epi_info.modifier & (EPI_FAST | EPI_REG)) && (!(epi_info.modifier & EPI_SET_BY_ALL))) {
	  epi_info.modifier |= EPI_SET_BY_SET;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et-p", 5)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --set-p must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0) || (dxx > 1)) {
	  sprintf(logbuf, "Error: Invalid --set-p parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	set_info.set_p = dxx;
      } else if (!memcmp(argptr2, "et-r2", 5)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --set-r2 must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        uii = 1;
        if (!strcmp(argv[cur_arg + 1], "write")) {
          uii = 2;
          set_info.modifier |= SET_R2_WRITE;
	} else if (param_ct == 2) {
          if (!strcmp(argv[cur_arg + 2], "write")) {
            set_info.modifier |= SET_R2_WRITE;
	  } else {
            sprintf(logbuf, "Error: Invalid --set-r2 parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        if (uii <= param_ct) {
	  if (scan_double(argv[cur_arg + uii], &dxx) || (dxx < 0.0)) {
	    sprintf(logbuf, "Error: Invalid --set-r2 parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (dxx > 0.0) {
	    // greater than 1 = no LD check.  (it still happens with a
	    // parameter of 1.)
	    set_info.set_r2 = dxx;
	  } else {
	    set_info.set_max = 1;
	  }
	}
      } else if (!memcmp(argptr2, "et-max", 5)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --set-max must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --set-max parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        set_info.set_max = ii;
      } else if (!memcmp(argptr2, "et-by-all", 10)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --set-by-all must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!(epi_info.modifier & (EPI_FAST | EPI_REG))) {
	  sprintf(logbuf, "Error: --set-by-all must be used with --{fast-}epistasis.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --set-by-all flag deprecated.  Use e.g. '--fast-epistasis set-by-all'.\n");
        epi_info.modifier |= EPI_SET_BY_ALL;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "nps-only", 9)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if ((strlen(argv[cur_arg + 1]) != 5) || (memcmp(argv[cur_arg + 1], "no-", 3)) || (!match_upper(&(argv[cur_arg + 1][3]), "DI"))) {
	    sprintf(logbuf, "Error: Invalid --snps-only parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_SNPS_ONLY_NO_DI;
	}
        misc_flags |= MISC_SNPS_ONLY;
      } else if (!memcmp(argptr2, "plit-x", 7)) {
	if (misc_flags & MISC_MERGEX) {
          sprintf(logbuf, "Error: --split-x cannot be used with --merge-x.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	} else if ((chrom_info.x_code == -1) || (chrom_info.xy_code == -1)) {
	  sprintf(logbuf, "Error: --split-x must be used with a chromosome set containing X and XY codes.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
          if ((!strcmp(argv[cur_arg + 1], "b36")) || (!strcmp(argv[cur_arg + 1], "hg18"))) {
            splitx_bound1 = 2709521;
            splitx_bound2 = 154584237;
	  } else if ((!strcmp(argv[cur_arg + 1], "b37")) || (!strcmp(argv[cur_arg + 1], "hg19"))) {
            splitx_bound1 = 2699520;
            splitx_bound2 = 154931044;
	  } else if ((!strcmp(argv[cur_arg + 1], "b38")) || (!strcmp(argv[cur_arg + 1], "hg20"))) {
            splitx_bound1 = 2781479;
            splitx_bound1 = 155701383;
	  } else {
            sprintf(logbuf, "Error: Unrecognized --split-x build code '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (chrom_info.species != SPECIES_HUMAN) {
	    logprint("Error: --split-x build codes cannot be used with nonhuman chromosome sets.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	} else {
	  if (atoiz(argv[cur_arg + 1], &ii)) {
	    sprintf(logbuf, "Error: Invalid --split-x parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  splitx_bound1 = (uint32_t)ii;
	  ii = atoi(argv[cur_arg + 2]);
	  if (ii <= ((int32_t)splitx_bound1)) {
	    sprintf(logbuf, "Error: Invalid --split-x parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  splitx_bound2 = (uint32_t)ii;
	}
      } else if (!memcmp(argptr2, "kato", 5)) {
	logprint("Error: --skato is not implemented yet.  Use e.g. PLINK/SEQ to perform this test\nfor now.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
      } else if (memcmp(argptr2, "ilent", 6)) {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 't':
      if (!memcmp(argptr2, "ail-pheno", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (makepheno_str) {
	  sprintf(logbuf, "Error: --tail-pheno cannot be used with --make-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &tail_bottom)) {
	  sprintf(logbuf, "Error: Invalid --tail-pheno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  tail_top = tail_bottom;
	} else {
	  if (scan_double(argv[cur_arg + 2], &tail_top)) {
	    sprintf(logbuf, "Error: Invalid --tail-pheno parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (tail_bottom > tail_top) {
	  sprintf(logbuf, "Error: Ltop cannot be larger than Hbottom for --tail-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	misc_flags |= MISC_TAIL_PHENO;
      } else if (!memcmp(argptr2, "hreads", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --threads parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (ii > MAX_THREADS) {
	  sprintf(logbuf, "Note: Reducing --threads parameter to %u.  (If this is not large enough,\nrecompile with a larger MAX_THREADS setting.)\n", MAX_THREADS);
	  logprintb();
	  ii = MAX_THREADS;
	}
	g_thread_ct = ii;
      } else if (!memcmp(argptr2, "ab", 3)) {
	logprint("Note: --tab flag deprecated.  Use '--recode tab ...'.\n");
	if (recode_modifier & RECODE_DELIMX) {
	  sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	recode_modifier |= RECODE_TAB;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ranspose", 9)) {
	logprint("Note: --transpose flag deprecated.  Use '--recode transpose ...'.\n");
	if (recode_modifier & RECODE_LGEN) {
	  sprintf(logbuf, "Error: --recode 'transpose' and 'lgen' modifiers cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	recode_modifier |= RECODE_TRANSPOSE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "fam", 4)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logprint("Error: --tfam filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(famname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TFAM;
      } else if (!memcmp(argptr2, "file", 5)) {
	if (load_params || (load_rare & (~LOAD_RARE_TFAM))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  jj = strlen(argv[cur_arg + 1]);
	  if (jj > FNAMESIZE - 6) {
	    logprint("Error: --tfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(pedname, argv[cur_arg + 1], jj), ".tped", 6);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(memcpya(famname, argv[cur_arg + 1], jj), ".tfam", 6);
	  }
	} else {
	  memcpy(pedname, PROG_NAME_STR ".tped", 11);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(famname, PROG_NAME_STR ".tfam", 11);
	  }
	}
	load_rare |= LOAD_RARE_TRANSPOSE;
      } else if (!memcmp(argptr2, "ped", 4)) {
	if (load_params || (load_rare & (~(LOAD_RARE_TRANSPOSE | LOAD_RARE_TFAM)))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logprint("Error: --tped filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TPED;
      } else if (!memcmp(argptr2, "o", 2)) {
	if ((!all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS)) || chrom_info.incl_excl_name_stack) {
	  sprintf(logbuf, "Error: --to cannot be used with --autosome{-xy} or --{not-}chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
	  sprintf(logbuf, "Error: --to cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
	  sprintf(logbuf, "Error: --to cannot be used with --snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "Error: --to cannot be used with a .cnv fileset.  Use --to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_to, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if ((!memcmp(argptr2, "o-bp", 5)) || (!memcmp(argptr2, "o-kb", 5)) || (!memcmp(argptr2, "o-mb", 5))) {
	if (markername_snp && (!(misc_flags & MISC_EXCLUDE_MARKERNAME_SNP))) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_to) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --to.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((!markername_from) && (marker_pos_start == -1)) {
	  marker_pos_start = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argptr2[2];
	if (cc == 'b') {
	  if (atoiz(argv[cur_arg + 1], &ii)) {
	    sprintf(logbuf, "Error: Invalid --to-bp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  if (marker_pos_end != -1) {
	    logprint("Error: Multiple --to-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (scan_double(argv[cur_arg + 1], &dxx)) {
	    sprintf(logbuf, "Error: Invalid --to-kb/-mb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    ii = 0;
	  } else if (dxx > 2147483647) {
	    ii = 0x7fffffff;
	  } else {
	    ii = (int32_t)(dxx * (1 + SMALL_EPSILON));
	  }
	}
	if (ii < marker_pos_start) {
	  marker_pos_end = marker_pos_start;
	  marker_pos_start = ii;
	} else {
	  marker_pos_end = ii;
	}
      } else if (!memcmp(argptr2, "rend", 5)) {
	if (model_modifier & MODEL_ASSOC) {
	  sprintf(logbuf, "Error: --trend cannot be used with --assoc.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_FISHER) {
	  sprintf(logbuf, "Error: --trend cannot be used with --model fisher.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	  sprintf(logbuf, "Error: --trend cannot be used with --model dom/rec/gen.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --trend flag deprecated.  Use '--model trend ...'.\n");
	model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
        calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "hin", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &thin_keep_prob)) {
	  sprintf(logbuf, "Error: Invalid --thin variant retention probability '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (thin_keep_prob < (0.5 / 4294967296.0)) {
	  sprintf(logbuf, "Error: --thin variant retention probability too small.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (thin_keep_prob >= (4294967295.5 / 4294967296.0)) {
	  sprintf(logbuf, "Error: --thin variant retention probability too large.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ests", 5)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --tests must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((param_ct == 1) && (!strcmp(argv[cur_arg + 1], "all"))) {
	  glm_modifier |= GLM_TEST_ALL;
	} else {
	  if (glm_modifier & GLM_TEST_ALL) {
	    sprintf(logbuf, "Error: --test-all cannot be used with --tests.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &tests_range_list, 1);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "est-all", 8)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --test-all must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --test-all flag deprecated.  Use '--tests all'.\n");
	glm_modifier |= GLM_TEST_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "wolocus", 8)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&(epi_info.twolocus_mkr1), argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
        if (alloc_string(&(epi_info.twolocus_mkr2), argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
	if (!strcmp(epi_info.twolocus_mkr1, epi_info.twolocus_mkr2)) {
          logprint("Error: --twolocus parameters cannot match.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        calculation_type |= CALC_EPI;
      } else if (!memcmp(argptr2, "est-missing", 12)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "perm")) {
            if (testmiss_modifier & TESTMISS_MPERM) {
              sprintf(logbuf, "Error: --test-missing 'mperm' and 'perm' cannot be used together.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
            testmiss_modifier |= TESTMISS_PERM;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
            if (testmiss_modifier & TESTMISS_PERM) {
              sprintf(logbuf, "Error: --test-missing 'mperm' and 'perm' cannot be used together.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    } else if (testmiss_modifier & TESTMISS_MPERM) {
	      // when --mperm and --test-missing mperm=[] are used together
              sprintf(logbuf, "Error: Duplicate --test-missing 'mperm' modifier.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
            ii = atoi(&(argv[cur_arg + uii][6]));
            if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --test-missing mperm parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
            testmiss_mperm_val = (uint32_t)ii;
            testmiss_modifier |= TESTMISS_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
            testmiss_modifier |= TESTMISS_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "midp")) {
            testmiss_modifier |= TESTMISS_MIDP;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
            logprint("Error: Improper --test-missing mperm syntax.  (Use\n'--test-missing mperm=[value]'.)\n");
            goto main_ret_INVALID_CMDLINE;
	  } else {
            sprintf(logbuf, "Error: Invalid --test-missing parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        calculation_type |= CALC_TESTMISS;
      } else if (!memcmp(argptr2, "est-mishap", 11)) {
        calculation_type |= CALC_TESTMISHAP;
        goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'u':
      if (!memcmp(argptr2, "nrelated-heritability", 22)) {
#ifdef NOLAPACK
        logprint("Error: --unrelated-heritability requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	UNSTABLE;
	if (rel_calc_type & REL_CALC_COV) {
	  sprintf(logbuf, "Error: --unrelated-heritability flag cannot coexist with a covariance\nmatrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --unrelated-heritability cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (rel_calc_type & REL_CALC_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --unrelated-heritability flag cannot be used with a single-precision\nrelationship matrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
	  if (calculation_type & CALC_REL_CUTOFF) {
	    sprintf(logbuf, "Error: --unrelated-heritability + --grm-gz/--grm-bin cannot be used with\n--rel-cutoff.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (!phenoname) {
	    sprintf(logbuf, "Error: --unrelated-heritability + --grm-gz/--grm-bin requires --pheno as well.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "strict")) {
	    misc_flags |= MISC_UNRELATED_HERITABILITY_STRICT;
	    uii = 2;
	  } else {
	    uii = 1;
	  }
	  if (param_ct >= uii) {
	    if (scan_double(argv[cur_arg + uii], &unrelated_herit_tol)) {
	      sprintf(logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (unrelated_herit_tol <= 0.0) {
	      sprintf(logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (param_ct > uii) {
	      if (scan_double(argv[cur_arg + uii + 1], &unrelated_herit_covg)) {
		sprintf(logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.%s", argv[cur_arg + uii + 1], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if ((unrelated_herit_covg <= 0.0) || (unrelated_herit_covg > 1.0)) {
		sprintf(logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.%s", argv[cur_arg + uii + 1], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if (param_ct == uii + 2) {
		if (scan_double(argv[cur_arg + uii + 2], &unrelated_herit_covr)) {
		  sprintf(logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.%s", argv[cur_arg + uii + 2], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
		if ((unrelated_herit_covr <= 0.0) || (unrelated_herit_covr > 1.0)) {
		  sprintf(logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.%s", argv[cur_arg + uii + 2], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
	      } else {
		unrelated_herit_covr = 1.0 - unrelated_herit_covg;
	      }
	    }
	  }
	}
	calculation_type |= CALC_UNRELATED_HERITABILITY;
#endif
      } else if (!memcmp(argptr2, "pdate-alleles", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&update_alleles_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "pdate-chr", 10)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-chr cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((misc_flags & MISC_MERGEX) || splitx_bound2) {
          sprintf(logbuf, "--update-chr cannot be used with --split-x or --merge-x.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!param_ct) {
	  update_map_modifier = 1;
	} else {
	  retval = alloc_2col(&update_chr, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-cm", 9)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "Error: --update-cm cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cm_map_fname) {
	  sprintf(logbuf, "Error: --update-cm cannot be used with --cm-map.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_map_modifier) {
	  logprint("Error: --update-map 'cm' modifier cannot be used with 'chr'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!param_ct) {
	  update_map_modifier = 2;
	} else {
	  retval = alloc_2col(&update_cm, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-ids", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&update_ids_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "pdate-map", 10)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-map cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_map_modifier) {
	  if (param_ct != 1) {
	    sprintf(logbuf, "Error: Multi-parameter --update-map cannot be used with deprecated\nparameter-free --update-%s.\n", (update_map_modifier == 1)? "chr" : "cm");
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  retval = alloc_2col((update_map_modifier == 1)? (&update_chr) : (&update_cm), &(argv[cur_arg + 1]), argptr, 1);
	  if (retval) {
	    goto main_ret_1;
	  }
	  sprintf(logbuf, "Note: --update-map [filename] + parameter-free --update-%s deprecated.  Use\n--update-%s [filename] instead.\n", (update_map_modifier == 1)? "chr" : "cm", (update_map_modifier == 1)? "chr" : "cm");
	  logprintb();
	  update_map_modifier = 0;
	} else {
	  retval = alloc_2col(&update_map, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-name", 11)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-name cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_chr) {
	  logprint("Error: --update-name cannot be used with --update-chr.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (update_cm) {
	  logprint("Error: --update-name cannot be used with --update-cm.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!param_ct) {
	  if (!update_map) {
	    sprintf(logbuf, "Error: Deprecated parameter-free --update-name cannot be used without\n--update-map.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((update_map->colx != 2) || (update_map->colid != 1) || (update_map->skip) || (update_map->skipchar)) {
	    logprint("Error: Multi-parameter --update-map cannot be used with deprecated\nparameter-free --update-name.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  logprint("Note: --update-map [filename] + parameter-free --update-name deprecated.  Use\n--update-name [filename] instead.\n");
	  update_name = update_map;
	  update_map = NULL;
	} else {
	  if (update_map) {
	    // no point in explaining the deprecated exception to this in the
	    // error message
	    logprint("Error: --update-name cannot be used with --update-map.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  retval = alloc_2col(&update_name, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-parents", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_ids_fname) {
	  logprint("Error: --update-parents cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&update_parents_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "pdate-sex", 10)) {
	if (update_ids_fname) {
	  logprint("Error: --update-sex cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&update_sex_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "nbounded", 9)) {
	if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --unbounded must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --unbounded flag deprecated.  Use '--genome unbounded'.\n");
	genome_modifier |= GENOME_IBD_UNBOUNDED;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'v':
      if (!memcmp(argptr2, "if", 3)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --vif must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &glm_vif_thresh)) {
	  sprintf(logbuf, "Error: Invalid --linear/--logistic VIF threshold '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (glm_vif_thresh < 1.0) {
	  sprintf(logbuf, "Error: --linear/--logistic VIF threshold '%s' too small (must be >= 1).%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "egas", 5)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --vegas must be used with --set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Error: --vegas is currently under development.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto main_ret_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "cf", 3)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = strlen(argv[cur_arg + 1]);
	if (uii > FNAMESIZE - 1) {
	  logprint("Error: --vcf filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], uii + 1);
	load_rare = LOAD_RARE_VCF;
      } else if (!memcmp(argptr2, "cf-min-qual", 12)) {
	if (!(load_rare & (LOAD_RARE_VCF | LOAD_RARE_BCF))) {
	  logprint("Error: --vcf-min-qual must be used with --vcf/--bcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &vcf_min_qual)) {
	  sprintf(logbuf, "Error: Invalid --vcf-min-qual parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	vcf_min_qual *= 1 - SMALL_EPSILON;
      } else if (!memcmp(argptr2, "cf-filter", 10)) {
	if (!(load_rare & (LOAD_RARE_VCF | LOAD_RARE_BCF))) {
	  logprint("Error: --vcf-filter must be used with --vcf/--bcf.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (param_ct) {
	  retval = alloc_and_flatten(&vcf_filter_exceptions_flattened, &(argv[cur_arg + 1]), param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
	misc_flags |= MISC_VCF_FILTER;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'w':
      if (!memcmp(argptr2, "rite-snplist", 13)) {
	calculation_type |= CALC_WRITE_SNPLIST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "indow", 6)) {
        if (!markername_snp) {
	  sprintf(logbuf, "Error: --window must be used with --snp or --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --window parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        dxx *= 500;
	if (dxx < 1) {
	  sprintf(logbuf, "Error: --window parameter '%s' too small.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (dxx > 2147483647) {
	  snp_window_size = 0x7fffffff;
	} else {
	  snp_window_size = (int32_t)(dxx * (1 + SMALL_EPSILON));
	}
      } else if (!memcmp(argptr2, "ithin", 6)) {
        if (loop_assoc_fname) {
	  sprintf(logbuf, "Error: --within cannot be used with --loop-assoc.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((calculation_type & CALC_FREQ) && (misc_flags & MISC_FREQ_COUNTS)) {
	  sprintf(logbuf, "Error: --within cannot be used with '--freq counts'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if ((strlen(argv[cur_arg + 1]) == 7) && (!memcmp(argv[cur_arg + 1], "keep-", 5)) && match_upper(&(argv[cur_arg + 1][5]), "NA")) {
	    uii = 2;
	  } else if ((strlen(argv[cur_arg + 2]) != 7) || memcmp(argv[cur_arg + 2], "keep-", 5) || (!match_upper(&(argv[cur_arg + 2][5]), "NA"))) {
            sprintf(logbuf, "Error: Invalid --within parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_LOAD_CLUSTER_KEEP_NA;
	}
	retval = alloc_fname(&cluster.fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ith-phenotype", 14)) {
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --with-phenotype cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-parents")) {
	    write_covar_modifier |= WRITE_COVAR_NO_PARENTS;
	  } else if (!strcmp(argv[cur_arg + uii], "no-sex")) {
	    if (write_covar_modifier & WRITE_COVAR_FEMALE_2) {
	      sprintf(logbuf, "Error: --with-phenotype 'female-2' modifier cannot be used with 'no-sex'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    write_covar_modifier |= WRITE_COVAR_NO_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "female-2")) {
	    if (write_covar_modifier & WRITE_COVAR_NO_SEX) {
	      sprintf(logbuf, "Error: --with-phenotype 'female-2' modifier cannot be used with 'no-sex'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    write_covar_modifier |= WRITE_COVAR_FEMALE_2;
	  } else {
	    sprintf(logbuf, "Error: Invalid --with-phenotype parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	write_covar_modifier |= WRITE_COVAR_PHENO;
      } else if (!memcmp(argptr2, "ith-reference", 14)) {
	if ((recode_modifier & RECODE_TYPEMASK) != RECODE_LGEN) {
	  logprint("Error: --with-reference must be used with --recode lgen.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --with-reference flag deprecated.  Use '--recode lgen-ref' instead.\n");
	recode_modifier += RECODE_LGEN_REF - RECODE_LGEN;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-covar", 11)) {
	if (calculation_type & (CALC_MAKE_BED | CALC_RECODE)) {
	  sprintf(logbuf, "Error: --write-covar cannot be used with --make-bed or --recode.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --write-covar cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_WRITE_COVAR;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-cluster", 13)) {
	if (!cluster.fname) {
	  sprintf(logbuf, "Error: --write-cluster must be used with --within.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "omit-unassigned")) {
	    sprintf(logbuf, "Error: Invalid --write-cluster parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_WRITE_CLUSTER_OMIT_UNASSIGNED;
	}
        calculation_type |= CALC_WRITE_CLUSTER;
      } else if (!memcmp(argptr2, "rite-set", 9)) {
	if (!set_info.fname) {
	  sprintf(logbuf, "Error: --write-set must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_WRITE_SET;
	set_info.modifier |= SET_WRITE_LIST;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-set-r2", 11)) {
        if (!set_info.fname) {
	  sprintf(logbuf, "Error: --write-set-r2 must be used with --set/--make-set.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        set_info.modifier |= SET_R2_WRITE;
        logprint("Note: --write-set-r2 deprecated.  Use '--set-r2 write'.\n");
        goto main_param_zero;
      } else if (!memcmp(argptr2, "ith-freqs", 10)) {
        ld_info.modifier |= LD_WITH_FREQS;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'x':
      if (!memcmp(argptr2, "chr-model", 10)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --xchr-model must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE)) {
	  sprintf(logbuf, "Error: --xchr-model cannot be used with --%s %s.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_GENOTYPIC)? "genotypic" : ((glm_modifier & GLM_HETHOM)? "hethom" : ((glm_modifier & GLM_DOMINANT)? "dominant" : "recessive")), errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((argv[cur_arg + 1][1] != '\0') || (argv[cur_arg + 1][0] < '0') || (argv[cur_arg + 1][0] > '3')) {
	  sprintf(logbuf, "Error: Invalid --xchr-model parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	glm_xchr_model = (uint32_t)(argv[cur_arg + 1][0] - '0');
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'z':
      if (!memcmp(argptr2, "ero-cluster", 12)) {
	if (!cluster.fname) {
	  sprintf(logbuf, "Error: --zero-cluster must be used with --within.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((calculation_type != CALC_MAKE_BED) || (geno_thresh != 1.0) || (mind_thresh != 1.0) || (hwe_thresh != 0.0) || (min_maf != 0.0) || (max_maf != 0.5)) {
	  // prevent old pipelines from silently breaking
	  sprintf(logbuf, "Error: --zero-cluster must now be used with --make-bed, no other output\ncommands, and no genotype-based filters.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cluster.zerofname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ero-cms", 8)) {
        misc_flags |= MISC_ZERO_CMS;
        goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    default:
      goto main_ret_INVALID_CMDLINE_2;

    main_param_zero:
      if (param_ct) {
        sprintf(logbuf, "Error: --%s doesn't accept parameters.%s", argptr, errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
  } while ((++cur_flag) < flag_ct);
  if (!outname_end) {
    outname_end = &(outname[5]);
  }

  // command-line restrictions which don't play well with alphabetical order
  if (load_rare) {
    if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
      if ((!(calculation_type & (CALC_REL_CUTOFF | CALC_UNRELATED_HERITABILITY))) || (calculation_type & (~(CALC_REL_CUTOFF | CALC_RELATIONSHIP | CALC_UNRELATED_HERITABILITY)))) {
	if (load_rare == LOAD_RARE_GRM) {
	  sprintf(logbuf, "Error: --grm-gz currently must be used with --rel-cutoff (possibly combined\nwith --make-grm-gz/--make-grm-bin) or --unrelated-heritability.%s", errstr_append);
	} else {
	  sprintf(logbuf, "Error: --grm-bin currently must be used with --rel-cutoff (possibly combined\nwith --make-grm-gz/--make-grm-bin) or --unrelated-heritability.%s", errstr_append);
	}
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
    if (!mperm_val) {
      if ((cnv_calc_type & CNV_INDIV_PERM) && (!cnv_indiv_mperms)) {
	sprintf(logbuf, "Error: --cnv-indiv-perm requires a permutation count.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if ((cnv_calc_type & CNV_TEST_REGION) && (!cnv_test_region_mperms)) {
	sprintf(logbuf, "Error: --cnv-test-region requires a permutation count.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
  }
  if ((cnv_intersect_filter_type & CNV_COUNT) && (!(cnv_calc_type & (CNV_INDIV_PERM | CNV_ENRICHMENT_TEST)))) {
    sprintf(logbuf, "Error: --cnv-count must be used with --cnv-indiv-perm or --cnv-enrichment-test.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (!phenoname) {
    if ((misc_flags & MISC_PRUNE) && (!(fam_cols & FAM_COL_6))) {
      sprintf(logbuf, "Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    } else if (pheno_modifier & PHENO_ALL) {
      sprintf(logbuf, "Error: --all-pheno must be used with --pheno.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  } else if (read_dists_fname && (!(calculation_type & (CALC_CLUSTER | CALC_IBS_TEST | CALC_GROUPDIST | CALC_NEIGHBOR | CALC_REGRESS_DISTANCE)))) {
    sprintf(logbuf, "Error: --read-dists cannot be used without --cluster, --ibs-test/--groupdist,\n--neighbour, or --regress-distance.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((cluster.ppc != 0.0) && (!read_genome_fname) && (calculation_type & (CALC_DISTANCE))) {
    sprintf(logbuf, "Error: --ppc cannot be used with --distance without --read-genome.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) && (((calculation_type & CALC_DISTANCE) && (dist_calc_type & (DISTANCE_1_MINUS_IBS | DISTANCE_ALCT))) || (calculation_type & CALC_PLINK1_DISTANCE_MATRIX)) && (!(calculation_type & CALC_GENOME))) {
    // actually allow this for now if --genome present, since it auto-clobbers
    // the wrong-unit distance matrix
    sprintf(logbuf, "Error: --cluster and --neighbour cannot be used with non-IBS distance matrix\ncalculations.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (matrix_flag_state == 1) {
    calculation_type |= CALC_PLINK1_IBS_MATRIX;
  }
  if (calculation_type & CALC_PLINK1_IBS_MATRIX) {
    if (exponent != 0.0) {
      logprint("Error: --ibs-matrix cannot be used with --distance-exp.\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (dist_calc_type & DISTANCE_IBS) {
      sprintf(logbuf, "Error: --ibs-matrix cannot be used with '--distance ibs'.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
    if (read_genome_fname && (cluster.ppc == 0.0)) {
      sprintf(logbuf, "Error: --read-genome is pointless with --ibs-matrix unless --ppc is also\npresent.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
    if (read_dists_fname) {
      sprintf(logbuf, "Error: --read-dists cannot be used with a distance matrix calculation.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
    if (parallel_tot > 1) {
      sprintf(logbuf, "Error: --parallel and --ibs-matrix cannot be used together.  Use --distance\ninstead.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  }
  if (update_map_modifier) {
    sprintf(logbuf, "Error: Deprecated parameter-free --update-%s cannot be used without\n--update-map.%s", (update_map_modifier == 1)? "chr" : "cm", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (((misc_flags & MISC_MERGEX) || splitx_bound2 || update_chr) && (((load_rare == LOAD_RARE_CNV) && (cnv_calc_type != CNV_WRITE)) || ((!load_rare) && (calculation_type != CALC_MAKE_BED)))) {
    sprintf(logbuf, "Error: --merge-x/--split-x/--update-chr must be used with --%s and no\nother commands.%s", load_rare? "cnv-write" : "make-bed", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (flip_subset_fname && (load_rare || (calculation_type != CALC_MAKE_BED) || (min_maf != 0.0) || (max_maf != 0.5) || (hwe_thresh != 0.0))) {
    sprintf(logbuf, "Error: --flip-subset must be used with --flip, --make-bed, and no other\ncommands or MAF-based filters.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (calculation_type & CALC_RECODE) {
    if (recode_modifier & (RECODE_23 | RECODE_BEAGLE)) {
      if (chrom_info.species != SPECIES_HUMAN) {
	sprintf(logbuf, "Error: --recode %s can only be used on human data.\n", (recode_modifier & RECODE_23)? "23" : "beagle");
	goto main_ret_INVALID_CMDLINE_3;
      }
      if ((recode_modifier & RECODE_23) && ((misc_flags & (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS)) == MISC_ALLOW_EXTRA_CHROMS)) {
	logprint("Error: --allow-extra-chr requires the '0' modifier when used with --recode 23.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if ((recode_modifier & RECODE_BEAGLE) && ((misc_flags & (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS)) == (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS))) {
        logprint("Error: --allow-extra-chr cannot have the '0' modifier when used with\n--recode beagle.\n");
	goto main_ret_INVALID_CMDLINE;
      }
    } else if ((recode_modifier & RECODE_BIMBAM) && ((misc_flags & (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS)) == MISC_ALLOW_EXTRA_CHROMS)) {
      logprint("Error: --allow-extra-chr requires the '0' modifier when used with\n--recode bimbam.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (sex_missing_pheno & MUST_HAVE_SEX) {
    if (load_rare & LOAD_RARE_CNV) {
      if (!(cnv_calc_type & CNV_WRITE)) {
        logprint("Error: --must-have-sex must be used with --cnv-write.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    } else {
      if (!(calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE))) {
        logprint("Error: --must-have-sex must be used with --make-bed/--recode/--write-covar.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
  }
  if (misc_flags & MISC_IMPUTE_SEX) {
    if (!(calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE))) {
      sprintf(logbuf, "Error: --impute-sex must be used with --make-bed/--recode/--write-covar.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    } else if (calculation_type & (~(CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE | CALC_SEXCHECK))) {
      sprintf(logbuf, "Error: --impute-sex cannot be used with any commands other than\n--make-bed/--recode/--write-covar.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  }
  if (cluster.qmatch_fname && (!cluster.qt_fname)) {
    sprintf(logbuf, "Error: --qt must be used with --qmatch.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (!cluster.fname) {
    if (mwithin_col && (!loop_assoc_fname)) {
      sprintf(logbuf, "Error: --mwithin must be used with --within.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    } else if (cluster.keep_fname) {
      logprint("Error: --keep-clusters must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.keep_flattened) {
      logprint("Error: --keep-cluster-names must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.remove_fname) {
      logprint("Error: --remove-clusters must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.remove_flattened) {
      logprint("Error: --remove-cluster-names must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (!set_info.fname) {
    if (set_info.modifier) {
      if (set_info.modifier & SET_COMPLEMENTS) {
	if (set_info.merged_set_name) {
	  logprint("Error: --make-set-complement-all must be used with --set/--make-set.\n");
	} else {
	  logprint("Error: --complement-sets must be used with --set/--make-set.\n");
	}
      } else { // only remaining possibility for now
        logprint("Error: --gene-all must be used with --set/--make-set.\n");
      }
      goto main_ret_INVALID_CMDLINE;
    } else if (set_info.genekeep_flattened) {
      logprint("Error: --gene must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (model_modifier & MODEL_SET_TEST) {
      logprint("Error: --assoc/--model set-test must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (glm_modifier & GLM_SET_TEST) {
      logprint("Error: --linear/--logistic set-test must be used with --set/--make-set.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  } else {
    if (model_modifier & MODEL_SET_TEST) {
      if ((!(model_modifier & MODEL_PERM)) && (!model_mperm_val)) {
        sprintf(logbuf, "Error: --assoc/--model set-test requires permutation.%s", errstr_append);
        goto main_ret_INVALID_CMDLINE_3;
      } else if (model_modifier & MODEL_FISHER) {
	sprintf(logbuf, "Error: --assoc/--model set-test cannot be used with Fisher's exact test.%s", errstr_append);
        goto main_ret_INVALID_CMDLINE_3;
      } else if ((mtest_adjust & ADJUST_GC) || (adjust_lambda != 0.0)) {
        sprintf(logbuf, "Error: --adjust 'gc' modifier and --lambda do not make sense with\n--assoc/--model set-test.%s", errstr_append);
        goto main_ret_INVALID_CMDLINE_3;
      }
    }
    if (glm_modifier & GLM_SET_TEST) {
      if ((!(glm_modifier & GLM_PERM)) && (!glm_mperm_val)) {
        sprintf(logbuf, "Error: --linear/--logistic set-test requires permutation.%s", errstr_append);
        goto main_ret_INVALID_CMDLINE_3;
      } else if ((mtest_adjust & ADJUST_GC) || (adjust_lambda != 0.0)) {
        sprintf(logbuf, "Error: --adjust 'gc' modifier and --lambda do not make sense with\n--linear/--logistic set-test.%s", errstr_append);
        goto main_ret_INVALID_CMDLINE_3;
      }
      logprint("Error: --linear/--logistic set-test is currently under development.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto main_ret_1;
    }
  }
  if ((!(calculation_type & CALC_LD)) || ((calculation_type & CALC_LD) && (ld_info.modifier & (LD_MATRIX_SHAPEMASK | LD_INTER_CHR)))) {
    if ((ld_info.snpstr || ld_info.snps_rl.name_ct) && (!(ld_info.modifier & LD_INTER_CHR))) {
      if (calculation_type & CALC_LD) {
	sprintf(logbuf, "Error: --ld-snp/--ld-snps/--ld-snp-list cannot be used with the --r/--r2 matrix\noutput modifiers.%s", errstr_append);
      } else {
        sprintf(logbuf, "Error: --ld-snp/--ld-snps/--ld-snp-list must be used with --r/--r2.%s", errstr_append);
      }
      goto main_ret_INVALID_CMDLINE_3;
    } else if (ld_info.window_size != 10) {
      if (calculation_type & CALC_LD) {
	sprintf(logbuf, "Error: --ld-window flag cannot be used with the --r/--r2 'inter-chr' or matrix\noutput modifiers.%s", errstr_append);
      } else {
        sprintf(logbuf, "Error: --ld-window flag must be used with --r/--r2.%s", errstr_append);
      }
      goto main_ret_INVALID_CMDLINE_3;
    } else if (ld_info.window_bp != 1000000) {
      if (calculation_type & CALC_LD) {
	sprintf(logbuf, "Error: --ld-window-kb flag cannot be used with the --r/--r2 'inter-chr' or\nmatrix output modifiers.%s", errstr_append);
      } else {
        sprintf(logbuf, "Error: --ld-window-kb flag must be used with --r/--r2.%s", errstr_append);
      }
      goto main_ret_INVALID_CMDLINE_3;
    } else if ((ld_info.window_r2 != 0.2) && (!(ld_info.modifier & LD_INTER_CHR))) {
      if (!(ld_info.modifier & LD_R2)) {
        logprint("Error: --ld-window-r2 flag must be used with --r2.\n");
        goto main_ret_INVALID_CMDLINE;
      } else {
	sprintf(logbuf, "Error: --ld-window-r2 flag cannot be used with the --r2 matrix output modifiers.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
  }
  if ((ld_info.modifier & LD_DPRIME) && (!(calculation_type & CALC_LD))) {
    sprintf(logbuf, "Error: --D/--dprime must be used with --r/--r2.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((ld_info.modifier & LD_WITH_FREQS) && (!(ld_info.modifier & LD_DPRIME))) {
    sprintf(logbuf, "Error: --r/--r2 'with-freqs' modifier must be used with 'dprime'.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }

  // --from-bp/-kb/-mb without any --to/--to-bp/...: include to end of
  // chromosome
  if ((marker_pos_start != -1) && (!markername_to) && (marker_pos_end == -1)) {
    marker_pos_end = 0x7fffffff;
  }
  if (!chrom_flag_present) {
    fill_chrom_mask(&chrom_info);
  }
  if (((marker_pos_start != -1) && (!markername_to)) || ((marker_pos_end != -1) && (!markername_from))) {
    // require exactly one chromosome to be defined given --from-bp/--to-bp
    // without --from/--to
    uii = next_set(chrom_info.chrom_mask, 0, CHROM_MASK_INITIAL_WORDS * BITCT);
    if (uii == CHROM_MASK_INITIAL_WORDS * BITCT) {
      uii = 0;
    } else {
      uii = next_set(chrom_info.chrom_mask, uii + 1, CHROM_MASK_INITIAL_WORDS * BITCT);
    }
    if (((uii == CHROM_MASK_INITIAL_WORDS * BITCT) && chrom_info.incl_excl_name_stack) || ((uii != CHROM_MASK_INITIAL_WORDS * BITCT) && (uii || (!chrom_info.incl_excl_name_stack) || chrom_info.incl_excl_name_stack->next))) {
      sprintf(logbuf, "Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb require a single chromosome to be\nidentified (either explicitly with --chr, or implicitly with --from/--to).%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  }

  if (mperm_save) {
    uii = 0;
    if ((calculation_type & CALC_MODEL) && (model_modifier & MODEL_MPERM)) {
      uii++;
    }
    if ((calculation_type & CALC_GLM) && (glm_modifier & GLM_MPERM)) {
      uii++;
    }
    if ((calculation_type & CALC_TESTMISS) && (testmiss_modifier & TESTMISS_MPERM)) {
      uii++;
    }
    if (uii != 1) {
      // prevent one permutation test's values from clobbering another's
      logprint("Error: --mperm-save{-all} must be used with exactly one max(T) permutation\ntest.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (calculation_type & CALC_MODEL) {
    if (!(model_modifier & (MODEL_ASSOC | MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
      if (mtest_adjust) {
	sprintf(logbuf, "Error: In order to use --model with --adjust, you must include the 'trend',\n'trend-only', 'dom', or 'rec' modifier.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
    if (model_cell_ct == -1) {
      model_cell_ct = (model_modifier & MODEL_FISHER)? 0 : 5;
    }
    if ((model_modifier & (MODEL_PERM | MODEL_MPERM | MODEL_GENEDROP)) == MODEL_GENEDROP) {
      model_modifier |= MODEL_PERM;
    }
  }
  if ((homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE)) && (!(calculation_type & CALC_HOMOZYG))) {
    if (homozyg.overlap_min == 0.95) {
      sprintf(logbuf, "Error: --homozyg-group must be used with another --homozyg... flag.%s", errstr_append);
    } else {
      sprintf(logbuf, "Error: --homozyg-match must be used with another --homozyg... flag.%s", errstr_append);
    }
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((load_params & 0x380) == 0x100) {
    sprintf(logbuf, "Error: --gen cannot be used without --data or --sample.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((!calculation_type) && (!(load_rare & (LOAD_RARE_LGEN | LOAD_RARE_DUMMY | LOAD_RARE_SIMULATE | LOAD_RARE_TRANSPOSE_MASK | LOAD_RARE_23 | LOAD_RARE_CNV | LOAD_RARE_VCF | LOAD_RARE_BCF))) && (famname[0] || load_rare)) {
    goto main_ret_NULL_CALC;
  }
  uii = 0; // short batch job?
  if (!(load_params || load_rare)) {
    if (!epi_info.summary_merge_prefix) {
      sprintf(logbuf, "Error: No input dataset.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
    // no input dataset needed in this case
    uii = 1;
  }

  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  subst_argv = NULL;
  script_buf = NULL;
  rerun_buf = NULL;
  flag_buf = NULL;
  flag_map = NULL;
  if (!rseeds) {
    sfmt_init_gen_rand(&sfmt, (uint32_t)time(NULL));
  } else {
    if (rseed_ct == 1) {
      sfmt_init_gen_rand(&sfmt, rseeds[0]);
    } else {
      sfmt_init_by_array(&sfmt, rseeds, rseed_ct);
    }
    free(rseeds);
    rseeds = NULL;
  }
  // guarantee contiguous malloc space outside of main workspace
  bubble = (char*)malloc(NON_WKSPACE_MIN * sizeof(char));
  if (!bubble) {
    goto main_ret_NOMEM;
  }

  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  sztmp = sizeof(int64_t);
  sysctl(mib, 2, &llxx, &sztmp, NULL, 0);
  llxx /= 1048576;
#else
#if _WIN32
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#else
  llxx = ((uint64_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
  if (!llxx) {
    default_alloc_mb = WKSPACE_DEFAULT_MB;
  } else if (llxx < (WKSPACE_MIN_MB * 2)) {
    default_alloc_mb = WKSPACE_MIN_MB;
  } else {
    default_alloc_mb = llxx / 2;
  }
  if (!malloc_size_mb) {
    malloc_size_mb = default_alloc_mb;
  } else if (malloc_size_mb < WKSPACE_MIN_MB) {
    malloc_size_mb = WKSPACE_MIN_MB;
  }
#ifndef __LP64__
  if (malloc_size_mb > 2047) {
    malloc_size_mb = 2047;
  }
#endif
  if (llxx) {
    sprintf(logbuf, "%" PRId64 " MB RAM detected; reserving %" PRIdPTR " MB for main workspace.\n", llxx, malloc_size_mb);
  } else {
    sprintf(logbuf, "Failed to calculate system memory.  Attempting to reserve %" PRIdPTR " MB.\n", malloc_size_mb);
  }
  logprintb();
  wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  while (!wkspace_ua) {
    malloc_size_mb = (malloc_size_mb * 3) / 4;
    if (malloc_size_mb < WKSPACE_MIN_MB) {
      malloc_size_mb = WKSPACE_MIN_MB;
    }
    wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (wkspace_ua) {
      sprintf(logbuf, "Allocated %" PRIdPTR " MB successfully, after larger attempt(s) failed.\n", malloc_size_mb);
      logprintb();
    } else if (malloc_size_mb == WKSPACE_MIN_MB) {
      goto main_ret_NOMEM;
    }
  }
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((uintptr_t)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = (malloc_size_mb * 1048576 - (uintptr_t)(wkspace - wkspace_ua)) & (~(CACHELINE - ONELU));
  free(bubble);

  if (epi_info.summary_merge_prefix) {
    retval = epi_summary_merge(&epi_info, outname, outname_end);
    // quit if there's nothing else to do
    if (retval || uii) {
      goto main_ret_1;
    }
  }

  pigz_init(g_thread_ct);
  if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
    // --unrelated-heritability and --rel-cutoff batch mode special cases
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      retval = unrelated_herit_batch(load_rare & LOAD_RARE_GRM_BIN, pedname, phenoname, mpheno_col, phenoname_str, missing_pheno, (misc_flags / MISC_UNRELATED_HERITABILITY_STRICT) & 1, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr);
    } else {
#endif
      retval = rel_cutoff_batch(load_rare & LOAD_RARE_GRM_BIN, pedname, outname, outname_end, rel_cutoff, rel_calc_type);
#ifndef NOLAPACK
    }
#endif
    /*
  } else if (genname[0]) {
    if (calculation_type & (~(CALC_DISTANCE | CALC_REGRESS_DISTANCE))) {
      logprint("Error: Only --distance calculations are currently supported with --data.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
    } else {
      if (!missing_code) {
	missing_code = (char*)"NA";
      }
      retval = plink_dosage(calculation_type, dist_calc_type, genname, samplename, outname, outname_end, missing_code, exponent, (misc_flags / MISC_MAF_SUCC) & 1, regress_iters, regress_d, g_thread_ct, parallel_idx, parallel_tot);
    }
    */
  } else if (load_rare & LOAD_RARE_CNV) {
    retval = plink_cnv(outname, outname_end, pedname, mapname, famname, phenoname, keepname, removename, filtername, misc_flags, update_chr, update_cm, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, filtervals_flattened, filter_binary, cnv_calc_type, cnv_min_seglen, cnv_max_seglen, cnv_min_score, cnv_max_score, cnv_min_sites, cnv_max_sites, cnv_intersect_filter_type, cnv_intersect_filter_fname, cnv_subset_fname, cnv_overlap_type, cnv_overlap_val, cnv_freq_type, cnv_freq_val, cnv_freq_val2, cnv_test_window, segment_modifier, segment_spanning_fname, cnv_indiv_mperms, cnv_test_mperms, cnv_test_region_mperms, cnv_enrichment_test_mperms, marker_pos_start, marker_pos_end, &chrom_info);
  } else if (load_rare & LOAD_RARE_GVAR) {
    retval = plink_gvar(outname, outname_end, pedname, mapname, famname);
  } else {
    // famname[0] indicates binary vs. text
    // extractname, excludename, keepname, and removename indicate the presence
    // of their respective flags
    // filtername indicates existence of filter
    // freqname signals --read-freq

    if (load_rare || (!famname[0])) {
      sptr = outname_end;
      if (calculation_type && (!(misc_flags & MISC_KEEP_AUTOCONV))) {
        sptr = memcpyb(sptr, "-temporary", 11);
      }
      uii = (sptr - outname);
      if (load_rare == LOAD_RARE_LGEN) {
        retval = lgen_to_bed(pedname, outname, sptr, missing_pheno, misc_flags, lgen_modifier, lgen_reference_fname, &chrom_info);
      } else if (load_rare & LOAD_RARE_TRANSPOSE_MASK) {
        retval = transposed_to_bed(pedname, famname, outname, sptr, misc_flags, &chrom_info);
      } else if (load_rare & LOAD_RARE_VCF) {
	retval = vcf_to_bed(pedname, outname, sptr, missing_pheno, misc_flags, const_fid, id_delim, vcf_min_qual, vcf_filter_exceptions_flattened, &chrom_info);
      } else if (load_rare & LOAD_RARE_BCF) {
	retval = bcf_to_bed(pedname, outname, sptr, missing_pheno, misc_flags, const_fid, id_delim, vcf_min_qual, vcf_filter_exceptions_flattened, &chrom_info);
      } else if (load_rare == LOAD_RARE_23) {
        retval = bed_from_23(pedname, outname, sptr, modifier_23, fid_23, iid_23, (pheno_23 == INFINITY)? ((double)missing_pheno) : pheno_23, paternal_id_23, maternal_id_23, &chrom_info);
      } else if (load_rare & LOAD_RARE_DUMMY) {
	retval = generate_dummy(outname, sptr, dummy_flags, dummy_marker_ct, dummy_indiv_ct, dummy_missing_geno, dummy_missing_pheno);
      } else if (load_rare & LOAD_RARE_SIMULATE) {
	retval = simulate_dataset(outname, sptr, simulate_flags, simulate_fname, simulate_cases, simulate_controls, simulate_prevalence, simulate_qt_indivs, simulate_missing, simulate_label);
	free(simulate_fname);
	simulate_fname = NULL;
	if (simulate_label) {
	  free(simulate_label);
	  simulate_label = NULL;
	}
      } else if (load_params & 0x380) {
	retval = oxford_to_bed(pedname, mapname, outname, sptr, hard_call_threshold, missing_code, missing_pheno, misc_flags, &chrom_info);
      } else {
	if (load_params & 0x30) {
	  sprintf(logbuf, "Error: --bed and --bim cannot be used without --bfile or --fam.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = ped_to_bed(pedname, mapname, outname, sptr, fam_cols, (misc_flags / MISC_AFFECTION_01) & 1, missing_pheno, &chrom_info);
	fam_cols |= FAM_COL_1 | FAM_COL_34 | FAM_COL_5;
	if (!(fam_cols & FAM_COL_6)) {
          fam_cols |= FAM_COL_6;
	  missing_pheno = -9;
	}
      }
      if (retval || (!calculation_type)) {
	goto main_ret_2;
      }
      memcpy(memcpya(pedname, outname, uii), ".bed", 5);
      memcpy(memcpya(mapname, outname, uii), ".bim", 5);
      memcpy(memcpya(famname, outname, uii), ".fam", 5);
      if (calculation_type && (!(misc_flags & MISC_KEEP_AUTOCONV))) {
	if (push_ll_str(&file_delete_list, pedname) || push_ll_str(&file_delete_list, mapname) || push_ll_str(&file_delete_list, famname)) {
	  goto main_ret_NOMEM;
	}
      }
      *outname_end = '\0';
    }
    if (ibc_type == 3) { // todo: make this less ugly
      ibc_type = 0;
    } else if (!ibc_type) {
      ibc_type = 1;
    }
    retval = plink(outname, outname_end, pedname, mapname, famname, cm_map_fname, cm_map_chrname, phenoname, extractname, excludename, keepname, removename, keepfamname, removefamname, filtername, freqname, read_dists_fname, read_dists_id_fname, evecname, mergename1, mergename2, mergename3, makepheno_str, phenoname_str, a1alleles, a2alleles, recode_allele_name, covar_fname, update_alleles_fname, read_genome_fname, update_chr, update_cm, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, loop_assoc_fname, flip_fname, flip_subset_fname, filtervals_flattened, condition_mname, condition_fname, filter_attrib_fname, filter_attrib_liststr, filter_attrib_indiv_fname, filter_attrib_indiv_liststr, thin_keep_prob, min_bp_space, mfilter_col, filter_binary, fam_cols, missing_pheno, output_missing_pheno, mpheno_col, pheno_modifier, &chrom_info, &oblig_missing_info, check_sex_fthresh, check_sex_mthresh, exponent, min_maf, max_maf, geno_thresh, mind_thresh, hwe_thresh, rel_cutoff, tail_bottom, tail_top, misc_flags, calculation_type, rel_calc_type, dist_calc_type, groupdist_iters, groupdist_d, regress_iters, regress_d, regress_rel_iters, regress_rel_d, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr, ibc_type, parallel_idx, parallel_tot, splitx_bound1, splitx_bound2, ppc_gap, sex_missing_pheno, hwe_modifier, genome_modifier, genome_min_pi_hat, genome_max_pi_hat, &homozyg, &cluster, neighbor_n1, neighbor_n2, &set_info, &ld_info, &epi_info, &clump_info, regress_pcs_modifier, max_pcs, recode_modifier, allelexxxx, merge_type, indiv_sort, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, &snps_range_list, covar_modifier, &covar_range_list, write_covar_modifier, write_covar_dummy_max_categories, mwithin_col, model_modifier, (uint32_t)model_cell_ct, model_mperm_val, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, &parameters_range_list, &tests_range_list, ci_size, pfilter, mtest_adjust, adjust_lambda, gxe_mcovar, &aperm, mperm_save, ibs_test_perms, perm_batch_size, lasso_h2, lasso_minlambda, &lasso_select_covars_range_list, testmiss_modifier, testmiss_mperm_val, &file_delete_list);
  }
 main_ret_2:
  free(wkspace_ua);
  while (0) {
  main_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  main_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  main_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  main_ret_INVALID_CMDLINE_2:
    invalid_arg(argv[cur_arg]);
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_4:
    sprintf(logbuf, "Error: --%s conflicts with another input flag.%s", argptr, errstr_append);
  main_ret_INVALID_CMDLINE_3:
    logprintb();
  main_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_NULL_CALC:
    logprint(null_calc_str);
    fputs(cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
#ifdef STABLE_BUILD
    break;
  main_unstable_disabled:
    logprint("Error: This flag's implementation is unfinished or unstable.  If you wish to\ntest it, use the latest development build.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
#endif
  }
 main_ret_1:
  fclose_cond(scriptfile);
  dispmsg(retval);
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  free_cond(makepheno_str);
  free_cond(phenoname_str);
  free_cond(a1alleles);
  free_cond(a2alleles);
  free_cond(filtervals_flattened);
  free_cond(evecname);
  free_cond(filtername);
  free_cond(read_dists_fname);
  free_cond(read_dists_id_fname);
  free_cond(freqname);
  free_cond(extractname);
  free_cond(excludename);
  free_cond(keepname);
  free_cond(removename);
  free_cond(keepfamname);
  free_cond(removefamname);
  free_cond(cm_map_fname);
  free_cond(cm_map_chrname);
  free_cond(phenoname);
  free_cond(recode_allele_name);
  free_cond(markername_from);
  free_cond(markername_to);
  free_cond(markername_snp);
  free_range_list(&snps_range_list);
  free_range_list(&covar_range_list);
  free_range_list(&lasso_select_covars_range_list);
  free_range_list(&parameters_range_list);
  free_range_list(&tests_range_list);
  free_cond(lgen_reference_fname);
  free_cond(covar_fname);
  free_cond(update_alleles_fname);
  free_cond(update_chr);
  free_cond(update_cm);
  free_cond(update_map);
  free_cond(update_name);
  free_cond(update_ids_fname);
  free_cond(update_parents_fname);
  free_cond(update_sex_fname);
  free_cond(loop_assoc_fname);
  free_cond(flip_fname);
  free_cond(flip_subset_fname);
  free_cond(read_genome_fname);
  free_cond(rseeds);
  free_cond(simulate_fname);
  free_cond(simulate_label);
  free_cond(cnv_intersect_filter_fname);
  free_cond(cnv_subset_fname);
  free_cond(segment_spanning_fname);
  free_cond(fid_23);
  free_cond(iid_23);
  free_cond(paternal_id_23);
  free_cond(maternal_id_23);
  free_cond(condition_mname);
  free_cond(condition_fname);
  free_cond(filter_attrib_fname);
  free_cond(filter_attrib_liststr);
  free_cond(filter_attrib_indiv_fname);
  free_cond(filter_attrib_indiv_liststr);
  free_cond(const_fid);
  free_cond(vcf_filter_exceptions_flattened);

  oblig_missing_cleanup(&oblig_missing_info);
  cluster_cleanup(&cluster);
  set_cleanup(&set_info);
  ld_epi_cleanup(&ld_info, &epi_info, &clump_info);
  if (file_delete_list) {
    do {
      ll_str_ptr = file_delete_list->next;
      unlink(file_delete_list->ss);
      free(file_delete_list);
      file_delete_list = ll_str_ptr;
    } while (file_delete_list);
  }
  forget_extra_chrom_names(&chrom_info);
  if (chrom_info.incl_excl_name_stack) {
    do {
      ll_str_ptr = chrom_info.incl_excl_name_stack->next;
      free(chrom_info.incl_excl_name_stack);
      chrom_info.incl_excl_name_stack = ll_str_ptr;
    } while (chrom_info.incl_excl_name_stack);
  }
  if (logfile) {
    if (!g_log_failed) {
      logstr("\nEnd time: ");
      time(&rawtime);
      logstr(ctime(&rawtime));
      if (fclose(logfile)) {
	fputs("Error: Failed to finish writing to log.\n", stdout);
      }
    } else {
      fclose(logfile);
    }
    logfile = NULL;
  }
  if (misc_flags & MISC_GPLINK) {
    memcpy(outname_end, ".gplink", 8);
    logfile = fopen(outname, "w");
    if (logfile) { // can't do much if an error occurs here...
      putc(retval? '1' : '0', logfile);
      putc('\n', logfile);
      fclose(logfile);
    }
  }

  return retval;
}
