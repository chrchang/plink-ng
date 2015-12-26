#include "plink_common.h"

#include "plink_lasso.h"
#include "plink_matrix.h"

#define WARM_START_ITERS 1000
#define NLAMBDA 100
#define DELTA_THRESHOLD 0.0001

int32_t transpose_covar(uintptr_t sample_valid_ct, uintptr_t covar_ct, uintptr_t* covar_nm, double* readptr, double* writeptr_start, double sqrt_n_recip, double sample_valid_ct_recip, double sample_valid_ctm1d) {
  // this may migrate to plink_common.c...
  double sum = 0.0;
  double ssq = 0.0;
  double* writeptr = writeptr_start;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  double subtract_by;
  double multiply_by;
  double dxx;
  for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_valid_ct; sample_uidx++, sample_idx++) {
    next_set_ul_unsafe_ck(covar_nm, &sample_uidx);
    dxx = readptr[sample_uidx * covar_ct];
    sum += dxx;
    ssq += dxx * dxx;
    *writeptr++ = dxx;
  }
  if (ssq * ((double)sample_valid_ct) == sum * sum) {
    return -1;
  }
  subtract_by = sum * sample_valid_ct_recip;
  multiply_by = sqrt_n_recip * sqrt(sample_valid_ctm1d / (ssq - sum * subtract_by));
  writeptr = writeptr_start;
  for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
    *writeptr = ((*writeptr) - subtract_by) * multiply_by;
    writeptr++;
  }
  return 0;
}

int32_t lasso_bigmem(FILE* bedfile, uintptr_t bed_offset, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm2, double lasso_h2, double lasso_minlambda, uint32_t select_covars, uintptr_t* select_covars_bitfield, double* pheno_d_collapsed, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uint32_t hh_or_mt_exists, uintptr_t sample_valid_ct, uintptr_t* sample_include2, uintptr_t* sample_male_include2, uintptr_t* loadbuf_raw, uintptr_t* loadbuf_collapsed, double* rand_matrix, double* misc_arr, double* residuals, uintptr_t* polymorphic_markers, uintptr_t* polymorphic_marker_ct_ptr, uint64_t* iter_tot_ptr, double** xhat_ptr) {
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  double* data_arr = (double*)g_bigstack_base; // marker-major
  double sqrt_n_recip = sqrt(1.0 / ((double)((intptr_t)sample_valid_ct)));
  double lambda_max = 0.0;
  double err_cur = 0.0;
  uint64_t iter_tot = 0;
  uintptr_t sample_valid_ctl2 = (sample_valid_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t polymorphic_marker_ct = 0;
  uintptr_t unselected_covar_ct = 0;
  uintptr_t final_mask = get_final_mask(sample_valid_ct);
  double lambda_min = lasso_minlambda;
  uint32_t chrom_fo_idx = 0xffffffffU; // exploit overflow
  uint32_t chrom_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t min_ploidy_1 = 0;
  int32_t retval = 0;
  double cur_mapping[4];
  double* xhat;
  double* prod_matrix;
  double* dptr;
  double* dptr2;
  uintptr_t* ulptr_end_init;
  uintptr_t* ulptr_end;
  uintptr_t* active_set;
  uintptr_t* ulptr;
  double sige;
  double loghi;
  double loglo;
  double logdelta;
  double lambda;
  double xjold;
  double err_last;
  double dxx;
  double dyy;
  double zz;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  uintptr_t ulii;
  uintptr_t iter;
  uintptr_t col_ct;
  uintptr_t col_ctl;
  uintptr_t col_idx;
  uintptr_t col_nz_ct;
  uintptr_t col_to_z;
  uintptr_t col_uidx;
  uintptr_t covar_idx;
  uintptr_t sample_idx;
  uintptr_t sample_idx_stop;
  uintptr_t marker_idx;
  uint32_t lambi;
  uint32_t marker_uidx;
  uint32_t homrar_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homset_ct;
  uint32_t uii;

  cur_mapping[1] = 0; // missing
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto lasso_bigmem_ret_READ_FAIL;
  }
  ulptr_end_init = &(loadbuf_collapsed[sample_valid_ct / BITCT2]);
  fputs("--lasso: Populating data matrix...", stdout);
  fflush(stdout);
  if (covar_ct) {
    dxx = 1.0 / ((double)((intptr_t)sample_valid_ct));
    dyy = (double)((intptr_t)(sample_valid_ct - 1));
    if (!select_covars_bitfield) {
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	if (transpose_covar(sample_valid_ct, covar_ct, covar_nm, &(covar_d[covar_idx]), &(data_arr[covar_idx * sample_valid_ct]), sqrt_n_recip, dxx, dyy)) {
	  goto lasso_bigmem_ret_CONST_COVAR;
	}
      }
      if (!select_covars) {
        unselected_covar_ct = covar_ct;
      }
    } else {
      ulii = 0;
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	if (IS_SET(select_covars_bitfield, covar_idx)) {
	  continue;
	}
	if (transpose_covar(sample_valid_ct, covar_ct, covar_nm, &(covar_d[covar_idx]), &(data_arr[ulii * sample_valid_ct]), sqrt_n_recip, dxx, dyy)) {
          goto lasso_bigmem_ret_CONST_COVAR;
	}
	ulii++;
      }
      unselected_covar_ct = ulii;
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	if (!IS_SET(select_covars_bitfield, covar_idx)) {
	  continue;
	}
	if (transpose_covar(sample_valid_ct, covar_ct, covar_nm, &(covar_d[covar_idx]), &(data_arr[ulii * sample_valid_ct]), sqrt_n_recip, dxx, dyy)) {
          goto lasso_bigmem_ret_CONST_COVAR;
	}
	ulii++;
      }
    }
  }
  dptr = &(data_arr[covar_ct * sample_valid_ct]);
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto lasso_bigmem_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
    }
    min_ploidy_1 |= uii;
    if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_sample_ct, loadbuf_collapsed, sample_valid_ct, pheno_nm2, final_mask, IS_SET(marker_reverse, marker_uidx))) {
      goto lasso_bigmem_ret_READ_FAIL;
    }
    if (min_ploidy_1) {
      haploid_fix(hh_or_mt_exists, sample_include2, sample_male_include2, sample_valid_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
    }
    vec_3freq(sample_valid_ctl2, loadbuf_collapsed, sample_include2, &missing_ct, &het_ct, &homset_ct);
    uii = sample_valid_ct - missing_ct;
    homrar_ct = uii - het_ct - homset_ct;
    if (!(((!homrar_ct) && ((!het_ct) || (!homset_ct))) || ((!het_ct) && (!homset_ct)))) {
      // ok, not monomorphic.  standardize to zero mean, unit variance
      SET_BIT(polymorphic_markers, marker_uidx);
      dyy = (double)(2 * homrar_ct + het_ct); // sum
      dxx = dyy / ((double)((int32_t)uii)); // mean
      dyy = sqrt_n_recip * sqrt(((double)((int32_t)(uii - 1))) / (4 * ((double)((int32_t)homrar_ct)) + ((double)((int32_t)het_ct)) - dyy * dxx)); // 1/(stdev * sqrt(n))
      cur_mapping[0] = (2 - dxx) * dyy; // 2x minor allele
      cur_mapping[2] = (1 - dxx) * dyy; // 1x minor allele
      cur_mapping[3] = (-dxx) * dyy; // no copies of minor allele
      sample_idx = 0;
      sample_idx_stop = BITCT2;
      ulptr = loadbuf_collapsed;
      ulptr_end = ulptr_end_init;
      while (1) {
        while (ulptr < ulptr_end) {
          cur_word = *ulptr++;
          for (; sample_idx < sample_idx_stop; sample_idx++, cur_word >>= 2) {
            cur_genotype = cur_word & 3;
            *dptr++ = cur_mapping[cur_genotype];
	  }
          sample_idx_stop += BITCT2;
	}
        if (sample_idx == sample_valid_ct) {
	  break;
	}
        ulptr_end++;
        sample_idx_stop = sample_valid_ct;
      }
      polymorphic_marker_ct++;
    }
  }
  *polymorphic_marker_ct_ptr = polymorphic_marker_ct;
  if (!polymorphic_marker_ct) {
    putchar('\n');
    logerrprint("Warning: Skipping --lasso since no polymorphic loci are present.\n");
    return 0;
  }
  col_ct = covar_ct + polymorphic_marker_ct;
  col_ctl = (col_ct + (BITCT - 1)) / BITCT;
  bigstack_shrink_top(data_arr, col_ct * sample_valid_ct * sizeof(double));
  sige = sqrt(1.0 - lasso_h2 + 1.0 / ((double)((intptr_t)sample_valid_ct)));
  zz = sige * sqrt_n_recip;
  if (rand_matrix) {
    bigstack_alloc_d(WARM_START_ITERS * WARM_START_ITERS, &prod_matrix);
    fputs("\r--lasso: Initializing warm start matrix...", stdout);
    fflush(stdout);
    fill_double_zero(misc_arr, WARM_START_ITERS);
    for (col_idx = 0; col_idx < col_ct;) {
      ulii = col_idx + WARM_START_ITERS;
      if (ulii > col_ct) {
	ulii = col_ct;
      }
      // splitting this into square blocks reduces memory consumption without
      // slowing things down (may even be faster due to locality).
      col_major_matrix_multiply(WARM_START_ITERS, ulii - col_idx, sample_valid_ct, rand_matrix, &(data_arr[col_idx * sample_valid_ct]), prod_matrix);
      dptr = prod_matrix;
      for (; col_idx < ulii; col_idx++) {
	for (uii = 0; uii < WARM_START_ITERS; uii++) {
	  dxx = fabs(*dptr++);
	  if (dxx > misc_arr[uii]) {
	    misc_arr[uii] = dxx;
	  }
	}
      }
    }
    lambda_min = destructive_get_dmedian(misc_arr, WARM_START_ITERS) * zz;
    logstr("--lasso:");
    LOGPRINTF(" using min lambda = %g.\n", lambda_min);
    bigstack_reset(prod_matrix);
  }
  bigstack_alloc_d(col_ct, &xhat);
  bigstack_alloc_ul(col_ctl, &active_set);
  *xhat_ptr = xhat;
  dptr = data_arr;
  for (col_idx = 0; col_idx < col_ct; col_idx++) {
    dxx = 0;
    dptr2 = pheno_d_collapsed;
    for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
      dxx += (*dptr++) * (*dptr2++);
    }
    dxx = fabs(dxx);
    if (dxx > lambda_max) {
      lambda_max = dxx;
    }
  }
  if (lambda_min >= lambda_max) {
    logprint("\n");
    logerrprint("Error: min lambda >= max lambda.\n");
    goto lasso_bigmem_ret_INVALID_CMDLINE;
  }
  loghi = log(lambda_max);
  loglo = log(lambda_min);
  logdelta = (loghi - loglo) / (NLAMBDA - 1.0);
  dptr = data_arr;
  for (col_idx = 0; col_idx < col_ct; col_idx++) {
    dxx = 0.0;
    dptr2 = pheno_d_collapsed;
    for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
      dxx += (*dptr++) * (*dptr2++);
    }
    xhat[col_idx] = dxx;
  }
  fputs("\r--lasso: Executing coordinate descent... ", stdout);
  for (lambi = 0; lambi < NLAMBDA; lambi++) {
    if (lambi > 10) {
      fputs("\b\b\b", stdout);
    } else if (lambi) {
      fputs("\b\b", stdout);
    }
    printf("%u%%", lambi); // only works since NLAMBDA is 100 for now
    fflush(stdout);
    lambda = exp(loghi - logdelta * ((double)((int32_t)lambi)));
    memcpy(residuals, pheno_d_collapsed, sample_valid_ct * sizeof(double));
    for (col_idx = 0; col_idx < col_ct; col_idx++) {
      dptr = residuals;
      dptr2 = &(data_arr[col_idx * sample_valid_ct]);
      dxx = -xhat[col_idx];
      for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
        *dptr += (*dptr2++) * dxx;
        dptr++;
      }
    }
    iter = 0;
    fill_all_bits(active_set, col_ct);
    col_nz_ct = col_ct;
    while (1) {
      col_uidx = 0;
      col_to_z = 0;
      for (marker_idx = 0; marker_idx < col_nz_ct; marker_idx++, col_uidx++) {
        col_uidx = next_set_unsafe(active_set, col_uidx);
        xjold = xhat[col_uidx];
        dxx = xjold;
        dptr = &(data_arr[col_uidx * sample_valid_ct]);
        dptr2 = residuals;
        for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
          dxx += (*dptr++) * (*dptr2++);
	}
	if (col_uidx >= unselected_covar_ct) {
	  if (dxx > 0.0) {
	    dxx = MAXV(dxx - lambda, 0.0);
	  } else {
	    dxx = MINV(dxx + lambda, 0.0);
	  }
	}
        xhat[col_uidx] = dxx;
        if (dxx == 0.0) {
          CLEAR_BIT(active_set, col_uidx);
	  col_to_z++;
	}
        dptr = residuals;
        dptr2 = &(data_arr[col_uidx * sample_valid_ct]);
        dxx -= xjold;
        for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
          *dptr -= (*dptr2++) * dxx;
          dptr++;
	}
      }
      col_nz_ct -= col_to_z;
      err_last = err_cur;
      err_cur = 0.0;
      col_uidx = 0;
      for (marker_idx = 0; marker_idx < col_nz_ct; marker_idx++, col_uidx++) {
        col_uidx = next_set_unsafe(active_set, col_uidx);
        err_cur += fabs(xhat[col_uidx]);
      }
      err_cur *= lambda;
      dptr = residuals;
      for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
        dxx = *dptr++;
        err_cur += dxx * dxx;
      }
      if (iter++) {
        if (((1.0 - (MINV(err_last, err_cur) / MAXV(err_last, err_cur))) < DELTA_THRESHOLD) || (err_cur != err_cur)) {
	  iter_tot += iter;
          break;
	}
      }
    }
  }
  *iter_tot_ptr = iter_tot;
  while (0) {
  lasso_bigmem_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  lasso_bigmem_ret_CONST_COVAR:
    logerrprint("Error: --lasso covariate is constant.\n");
  lasso_bigmem_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

uint32_t load_and_normalize(FILE* bedfile, uintptr_t* loadbuf_raw, uintptr_t unfiltered_sample_ct, uintptr_t* loadbuf_collapsed, uintptr_t sample_valid_ct, uintptr_t* pheno_nm2, uintptr_t final_mask, uint32_t do_reverse, uint32_t min_ploidy_1, uint32_t hh_or_mt_exists, uintptr_t* sample_include2, uintptr_t* sample_male_include2, uint32_t is_x, uint32_t is_y, double sqrt_n_recip, double* data_window_ptr) {
  uintptr_t sample_valid_ctl2 = (sample_valid_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t sample_idx = 0;
  uintptr_t sample_idx_stop = BITCT2;
  uintptr_t* ulptr_end_init = &(loadbuf_collapsed[sample_valid_ct / BITCT2]);
  double cur_mapping[4];
  uintptr_t* ulptr = loadbuf_collapsed;
  uintptr_t* ulptr_end = ulptr_end_init;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  double dxx;
  double dyy;
  uint32_t homrar_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homset_ct;
  uint32_t uii;
  if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_sample_ct, loadbuf_collapsed, sample_valid_ct, pheno_nm2, final_mask, do_reverse)) {
    return 2; // read failure
  }
  if (min_ploidy_1) {
    haploid_fix(hh_or_mt_exists, sample_include2, sample_male_include2, sample_valid_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
  }
  vec_3freq(sample_valid_ctl2, loadbuf_collapsed, sample_include2, &missing_ct, &het_ct, &homset_ct);
  uii = sample_valid_ct - missing_ct;
  homrar_ct = uii - het_ct - homset_ct;
  if (((!homrar_ct) && ((!het_ct) || (!homset_ct))) || ((!het_ct) && (!homset_ct))) {
    return 1; // not polymorphic
  }
  dyy = (double)(2 * homrar_ct + het_ct); // sum
  dxx = dyy / ((double)((int32_t)uii)); // mean
  dyy = sqrt_n_recip * sqrt(((double)((int32_t)(uii - 1))) / (4 * ((double)((int32_t)homrar_ct)) + ((double)((int32_t)het_ct)) - dyy * dxx)); // 1/(stdev * sqrt(n))
  cur_mapping[0] = (2 - dxx) * dyy; // 2x minor allele
  cur_mapping[1] = 0.0;
  cur_mapping[2] = (1 - dxx) * dyy; // 1x minor allele
  cur_mapping[3] = (-dxx) * dyy; // no copies of minor allele
  while (1) {
    while (ulptr < ulptr_end) {
      cur_word = *ulptr++;
      for (; sample_idx < sample_idx_stop; sample_idx++, cur_word >>= 2) {
	cur_genotype = cur_word & 3;
	*data_window_ptr++ = cur_mapping[cur_genotype];
      }
      sample_idx_stop += BITCT2;
    }
    if (sample_idx == sample_valid_ct) {
      return 0;
    }
    ulptr_end++;
    sample_idx_stop = sample_valid_ct;
  }
}


int32_t lasso_smallmem(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm2, double lasso_h2, double lasso_minlambda, uint32_t select_covars, uintptr_t* select_covars_bitfield, double* pheno_d_collapsed, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uint32_t hh_or_mt_exists, uintptr_t sample_valid_ct, uintptr_t* sample_include2, uintptr_t* sample_male_include2, uintptr_t* loadbuf_raw, uintptr_t* loadbuf_collapsed, double* rand_matrix, double* misc_arr, double* residuals, uintptr_t* polymorphic_markers, uintptr_t* polymorphic_marker_ct_ptr, uint64_t* iter_tot_ptr, double** xhat_ptr) {
  // Instead of populating and normalizing data_arr before the coordinate
  // descent, we reload and renormalize the data every iteration.
  // Since (i) there's probably a larger number of samples involved, and (ii)
  // there's also more computational work to do per sample, multithreading is
  // more profitable here than in the bigmem case.
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  double* covar_data_arr = NULL;
  double sqrt_n_recip = sqrt(1.0 / ((double)((intptr_t)sample_valid_ct)));
  double lambda_max = 0.0;
  double err_cur = 0.0;
  double lambda_min = lasso_minlambda;
  uint64_t iter_tot = 0;
  uintptr_t unselected_covar_ct = 0;
  uintptr_t final_mask = get_final_mask(sample_valid_ct);
  uintptr_t polymorphic_marker_ct = 0;
  uint32_t chrom_fo_idx = 0xffffffffU; // exploit overflow
  uint32_t chrom_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t min_ploidy_1 = 0;
  uint32_t partial_marker_idx = 0;
  int32_t retval = 0;
  double* prod_matrix = NULL;
  double* data_window;
  double* xhat;
  double* dptr;
  double* dptr2;
  uintptr_t* active_set;
  double sige;
  double loghi;
  double loglo;
  double logdelta;
  double lambda;
  double xjold;
  double err_last;
  double dxx;
  double dyy;
  double zz;
  uintptr_t ulii;
  uintptr_t iter;
  uintptr_t col_ct;
  uintptr_t col_ctl;
  uintptr_t col_idx;
  uintptr_t col_nz_ct;
  uintptr_t col_to_z;
  uintptr_t col_uidx;
  uintptr_t covar_idx;
  uintptr_t sample_idx;
  uintptr_t marker_idx;
  uint32_t lambi;
  uint32_t marker_uidx;
  uint32_t last_marker_uidx;
  uint32_t uii;

  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto lasso_smallmem_ret_READ_FAIL;
  }
  fputs("Using memory-conserving LASSO implementation.\n", stdout);
  if (covar_ct) {
    if (bigstack_alloc_d(covar_ct * sample_valid_ct, &covar_data_arr)) {
      goto lasso_smallmem_ret_NOMEM;
    }
    dxx = 1.0 / ((double)((intptr_t)sample_valid_ct));
    dyy = (double)((intptr_t)(sample_valid_ct - 1));
    if (!select_covars_bitfield) {
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	if (transpose_covar(sample_valid_ct, covar_ct, covar_nm, &(covar_d[covar_idx]), &(covar_data_arr[covar_idx * sample_valid_ct]), sqrt_n_recip, dxx, dyy)) {
	  goto lasso_smallmem_ret_CONST_COVAR;
	}
      }
      if (!select_covars) {
	unselected_covar_ct = covar_ct;
      }
    } else {
      ulii = 0;
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	if (IS_SET(select_covars_bitfield, covar_idx)) {
	  continue;
	}
	if (transpose_covar(sample_valid_ct, covar_ct, covar_nm, &(covar_d[covar_idx]), &(covar_data_arr[ulii * sample_valid_ct]), sqrt_n_recip, dxx, dyy)) {
          goto lasso_smallmem_ret_CONST_COVAR;
	}
	ulii++;
      }
      unselected_covar_ct = ulii;
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	if (!IS_SET(select_covars_bitfield, covar_idx)) {
	  continue;
	}
	if (transpose_covar(sample_valid_ct, covar_ct, covar_nm, &(covar_d[covar_idx]), &(covar_data_arr[ulii * sample_valid_ct]), sqrt_n_recip, dxx, dyy)) {
          goto lasso_smallmem_ret_CONST_COVAR;
	}
	ulii++;
      }
    }
  }
  sige = sqrt(1.0 - lasso_h2 + 1.0 / ((double)((intptr_t)sample_valid_ct)));
  zz = sige * sqrt_n_recip;
  // put this on top of the permanent stack portion so we can shrink it when we
  // know the true column count
  if (bigstack_alloc_d(covar_ct + marker_ct, &xhat)) {
    goto lasso_smallmem_ret_NOMEM;
  }
  if (rand_matrix) {
    if (bigstack_alloc_d(sample_valid_ct * WARM_START_ITERS, &data_window) ||
        bigstack_alloc_d(WARM_START_ITERS * WARM_START_ITERS, &prod_matrix)) {
      goto lasso_smallmem_ret_NOMEM;
    }
    fputs("\r--lasso: Initializing warm start matrix...", stdout);
    fflush(stdout);
    fill_double_zero(misc_arr, WARM_START_ITERS);
  } else {
    if (bigstack_alloc_d(sample_valid_ct, &data_window)) {
      goto lasso_smallmem_ret_NOMEM;
    }
  }
  dptr = data_window;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    // count polymorphic markers here, compute xhat[] and lambda_max, perform
    // rand_matrix multiply if necessary
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	goto lasso_smallmem_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
    }
    uii = load_and_normalize(bedfile, loadbuf_raw, unfiltered_sample_ct, loadbuf_collapsed, sample_valid_ct, pheno_nm2, final_mask, IS_SET(marker_reverse, marker_uidx), min_ploidy_1, hh_or_mt_exists, sample_include2, sample_male_include2, is_x, is_y, sqrt_n_recip, dptr);
    if (uii == 2) {
      goto lasso_smallmem_ret_READ_FAIL;
    }
    if (uii == 1) {
      continue;
    }
    SET_BIT(polymorphic_markers, marker_uidx);
    dxx = 0.0;
    for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
      dxx += dptr[sample_idx] * pheno_d_collapsed[sample_idx];
    }
    xhat[covar_ct + polymorphic_marker_ct + partial_marker_idx] = dxx;
    partial_marker_idx++;
    dxx = fabs(dxx);
    if (dxx > lambda_max) {
      lambda_max = dxx;
    }
    if (partial_marker_idx == WARM_START_ITERS) {
      if (rand_matrix) {
	col_major_matrix_multiply(WARM_START_ITERS, WARM_START_ITERS, sample_valid_ct, rand_matrix, data_window, prod_matrix);
	dptr = prod_matrix;
	for (col_idx = 0; col_idx < WARM_START_ITERS; col_idx++) {
	  for (ulii = 0; ulii < WARM_START_ITERS; ulii++) {
	    dxx = fabs(*dptr++);
	    if (dxx > misc_arr[ulii]) {
	      misc_arr[ulii] = dxx;
	    }
	  }
	}
      }
      polymorphic_marker_ct += WARM_START_ITERS;
      partial_marker_idx = 0;
      dptr = data_window;
    } else if (!rand_matrix) {
      dptr = data_window;
    } else {
      dptr = &(dptr[sample_valid_ct]);
    }
  }
  if (!(polymorphic_marker_ct + partial_marker_idx)) {
    logerrprint("Warning: Skipping --lasso since no polymorphic markers are present.\n");
    return 0;
  }
  if (rand_matrix) {
    if (partial_marker_idx) {
      col_major_matrix_multiply(WARM_START_ITERS, partial_marker_idx, sample_valid_ct, rand_matrix, data_window, prod_matrix);
      dptr = prod_matrix;
      for (col_idx = 0; col_idx < partial_marker_idx; col_idx++) {
	for (ulii = 0; ulii < WARM_START_ITERS; ulii++) {
	  dxx = fabs(*dptr++);
	  if (dxx > misc_arr[ulii]) {
	    misc_arr[ulii] = dxx;
	  }
	}
      }
    }
    lambda_min = destructive_get_dmedian(misc_arr, WARM_START_ITERS) * zz;
    logstr("--lasso:");
    LOGPRINTF(" using min lambda = %g.\n", lambda_min);
  }
  polymorphic_marker_ct += partial_marker_idx;
  *polymorphic_marker_ct_ptr = polymorphic_marker_ct;
  col_ct = covar_ct + polymorphic_marker_ct;
  bigstack_reset(data_window);
  bigstack_shrink_top(xhat, col_ct * sizeof(double));
  col_ctl = (col_ct + (BITCT - 1)) / BITCT;
  *xhat_ptr = xhat;
  if (lambda_min >= lambda_max) {
    logprint("\n");
    logerrprint("Error: min lambda >= max lambda.\n");
    goto lasso_smallmem_ret_INVALID_CMDLINE;
  }
  if (bigstack_alloc_ul(col_ctl, &active_set) ||
      bigstack_alloc_d(sample_valid_ct, &data_window)) {
    goto lasso_smallmem_ret_NOMEM;
  }
  loghi = log(lambda_max);
  loglo = log(lambda_min);
  logdelta = (loghi - loglo) / (NLAMBDA - 1.0);
  fputs("\r--lasso: Executing coordinate descent... ", stdout);
  for (lambi = 0; lambi < NLAMBDA; lambi++) {
    if (lambi > 10) {
      fputs("\b\b\b", stdout);
    } else if (lambi) {
      fputs("\b\b", stdout);
    }
    printf("%u%%", lambi); // only works since NLAMBDA is 100 for now
    fflush(stdout);
    lambda = exp(loghi - logdelta * ((double)((int32_t)lambi)));
    memcpy(residuals, pheno_d_collapsed, sample_valid_ct * sizeof(double));
    for (col_idx = 0; col_idx < covar_ct; col_idx++) {
      dptr = residuals;
      dptr2 = &(covar_data_arr[col_idx * sample_valid_ct]);
      dxx = -xhat[col_idx];
      for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
        *dptr += (*dptr2++) * dxx;
        dptr++;
      }
    }
    if (fseeko(bedfile, bed_offset, SEEK_SET)) {
      goto lasso_smallmem_ret_READ_FAIL;
    }
    for (marker_uidx = 0, marker_idx = 0, chrom_end = 0, chrom_fo_idx = 0xffffffffU; marker_idx < polymorphic_marker_ct; marker_uidx++, marker_idx++) {
      if (!IS_SET(polymorphic_markers, marker_uidx)) {
	marker_uidx = next_set_unsafe(polymorphic_markers, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto lasso_smallmem_ret_READ_FAIL;
	}
      }
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
      }
      if (load_and_normalize(bedfile, loadbuf_raw, unfiltered_sample_ct, loadbuf_collapsed, sample_valid_ct, pheno_nm2, final_mask, IS_SET(marker_reverse, marker_uidx), min_ploidy_1, hh_or_mt_exists, sample_include2, sample_male_include2, is_x, is_y, sqrt_n_recip, data_window)) {
	goto lasso_smallmem_ret_READ_FAIL;
      }
      dptr = residuals;
      dptr2 = data_window;
      dxx = -xhat[marker_idx + covar_ct];
      for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
	*dptr += (*dptr2++) * dxx;
	dptr++;
      }
    }
    iter = 0;
    fill_all_bits(active_set, col_ct);
    col_nz_ct = col_ct;
    while (1) {
      col_uidx = 0;
      col_to_z = 0;
      for (marker_uidx = 0, last_marker_uidx = 0xfffffffeU, chrom_end = 0, chrom_fo_idx = 0xffffffffU, marker_idx = 0; marker_idx < col_nz_ct; col_uidx++) {
	if (col_uidx < covar_ct) {
	  if (!(IS_SET(active_set, col_uidx))) {
	    continue;
	  }
	  dptr = &(covar_data_arr[col_uidx * sample_valid_ct]);
	} else {
          marker_uidx = next_set_unsafe(polymorphic_markers, marker_uidx);
	  if (!(IS_SET(active_set, col_uidx))) {
	    marker_uidx++;
	    continue;
	  }
	  if (marker_uidx != last_marker_uidx + 1) {
	    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	      goto lasso_smallmem_ret_READ_FAIL;
	    }
	  }
	  if (marker_uidx >= chrom_end) {
	    chrom_fo_idx++;
	    refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
	  }
	  if (load_and_normalize(bedfile, loadbuf_raw, unfiltered_sample_ct, loadbuf_collapsed, sample_valid_ct, pheno_nm2, final_mask, IS_SET(marker_reverse, marker_uidx), min_ploidy_1, hh_or_mt_exists, sample_include2, sample_male_include2, is_x, is_y, sqrt_n_recip, data_window)) {
	    goto lasso_smallmem_ret_READ_FAIL;
	  }
	  dptr = data_window;
	}
	marker_idx++;
        xjold = xhat[col_uidx];
        dxx = xjold;
        dptr2 = residuals;
        for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
          dxx += (*dptr++) * (*dptr2++);
	}
	if (col_uidx >= unselected_covar_ct) {
	  if (dxx > 0.0) {
	    dxx = MAXV(dxx - lambda, 0.0);
	  } else {
	    dxx = MINV(dxx + lambda, 0.0);
	  }
	}
        xhat[col_uidx] = dxx;
        if (dxx == 0.0) {
          CLEAR_BIT(active_set, col_uidx);
	  col_to_z++;
	}
        dptr = residuals;
	if (col_uidx < covar_ct) {
          dptr2 = &(covar_data_arr[col_uidx * sample_valid_ct]);
	} else {
	  dptr2 = data_window;
	  last_marker_uidx = marker_uidx;
	  marker_uidx++;
	}
        dxx -= xjold;
        for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
          *dptr -= (*dptr2++) * dxx;
          dptr++;
	}
      }
      col_nz_ct -= col_to_z;
      err_last = err_cur;
      err_cur = 0.0;
      for (col_uidx = 0, marker_idx = 0; marker_idx < col_nz_ct; marker_idx++, col_uidx++) {
        col_uidx = next_set_unsafe(active_set, col_uidx);
        err_cur += fabs(xhat[col_uidx]);
      }
      err_cur *= lambda;
      dptr = residuals;
      for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
        dxx = *dptr++;
        err_cur += dxx * dxx;
      }
      if (iter++) {
        if (((1.0 - (MINV(err_last, err_cur) / MAXV(err_last, err_cur))) < DELTA_THRESHOLD) || (err_cur != err_cur)) {
	  iter_tot += iter;
          break;
	}
      }
    }
  }
  *iter_tot_ptr = iter_tot;
  while (0) {
  lasso_smallmem_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  lasso_smallmem_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  lasso_smallmem_ret_CONST_COVAR:
    logerrprint("Error: --lasso covariate is constant.\n");
  lasso_smallmem_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

int32_t lasso(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t pheno_nm_ct, double lasso_h2, double lasso_minlambda, Range_list* select_covars_range_list_ptr, uint64_t misc_flags, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_male, uint32_t hh_or_mt_exists) {
  // Coordinate descent LASSO.  Based on a MATLAB script by Shashaank
  // Vattikuti.
  // Not yet multithreaded.  (Main loop is fairly tightly coupled, so getting
  // a performance benefit will be a bit tricky.)
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_sample_ctv2 = 2 * ((unfiltered_sample_ct + (BITCT - 1)) / BITCT);
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t polymorphic_marker_ct = 0;
  uint64_t iter_tot = 0;
  double* xhat = NULL;
  double* rand_matrix = NULL;
  double* misc_arr = NULL;
  uintptr_t* sample_male_include2 = NULL;
  uintptr_t* select_covars_bitfield = NULL;
  char* wptr_start = NULL;
  uint32_t report_zeroes = (misc_flags / MISC_LASSO_REPORT_ZEROES) & 1;
  uint32_t select_covars = (misc_flags / MISC_LASSO_SELECT_COVARS) & 1;
  uint32_t chrom_fo_idx = 0xffffffffU; // exploit overflow
  uint32_t chrom_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t min_ploidy_1 = 0;
  int32_t retval = 0;
  double* pheno_d_collapsed;
  double* residuals;
  double* dptr;
  double* dptr2;
  double* dptr3;
  uintptr_t* pheno_nm2;
  uintptr_t* sample_include2;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_collapsed;
  uintptr_t* polymorphic_markers;
  char* wptr;
  char* wptr2;
  double sqrt_n_recip;
  double dxx;
  double dyy;
  double dzz;
  uint64_t ullii;
  uintptr_t sample_valid_ct;
  uintptr_t sample_valid_ctv2;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t sample_uidx;
  uintptr_t sample_uidx_stop;
  uintptr_t sample_idx;
  uint32_t marker_uidx;
  uint32_t uii;
  if (!covar_ct) {
    sample_valid_ct = pheno_nm_ct;
  } else {
    sample_valid_ct = popcount_longs(covar_nm, (pheno_nm_ct + (BITCT - 1)) / BITCT);
  }
  if (sample_valid_ct < 2) {
    if (pheno_nm_ct < 2) {
      logerrprint("Warning: Skipping --lasso since less than two phenotypes are present.\n");
    } else {
      logerrprint("Warning: Skipping --lasso since too many samples have missing covariates.\n");
    }
    goto lasso_ret_1;
  }
  if (sample_valid_ct == pheno_nm_ct) {
    pheno_nm2 = pheno_nm;
  } else {
    if (bigstack_calloc_ul(unfiltered_sample_ctl, &pheno_nm2)) {
      goto lasso_ret_NOMEM;
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < pheno_nm_ct; sample_uidx++, sample_idx++) {
      next_set_ul_unsafe_ck(pheno_nm, &sample_uidx);
      if (IS_SET(covar_nm, sample_idx)) {
        SET_BIT(pheno_nm2, sample_uidx);
      }
    }
  }
  sample_valid_ctv2 = 2 * ((sample_valid_ct + (BITCT - 1)) / BITCT);
  sqrt_n_recip = sqrt(1.0 / ((double)((intptr_t)sample_valid_ct)));
  if (bigstack_alloc_ul(sample_valid_ctv2, &sample_include2) ||
      bigstack_alloc_ul(unfiltered_sample_ctv2, &loadbuf_raw) ||
      bigstack_alloc_ul(sample_valid_ctv2, &loadbuf_collapsed) ||
      bigstack_alloc_ul(unfiltered_marker_ctl, &polymorphic_markers) ||
      bigstack_alloc_d(sample_valid_ct, &pheno_d_collapsed) ||
      bigstack_alloc_d(sample_valid_ct, &residuals)) {
    goto lasso_ret_NOMEM;
  }
  if (lasso_minlambda == -1) {
    if (bigstack_alloc_d(sample_valid_ct * WARM_START_ITERS, &rand_matrix) ||
        bigstack_alloc_d(WARM_START_ITERS, &misc_arr)) {
      goto lasso_ret_NOMEM;
    }
  }
  dxx = 0; // sum
  dyy = 0; // ssq
  dptr = pheno_d_collapsed;
  sample_uidx = 0;
  sample_idx = 0;
  if (pheno_d) {
    do {
      sample_uidx = next_set_ul_unsafe(pheno_nm2, sample_uidx);
      sample_uidx_stop = next_unset_ul(pheno_nm2, sample_uidx, unfiltered_sample_ct);
      sample_idx += sample_uidx_stop - sample_uidx;
      dptr2 = &(pheno_d[sample_uidx]);
      sample_uidx = sample_uidx_stop;
      dptr3 = &(pheno_d[sample_uidx_stop]);
      do {
	dzz = *dptr2++;
	*dptr++ = dzz;
	dxx += dzz;
	dyy += dzz * dzz;
      } while (dptr2 < dptr3);
    } while (sample_idx < sample_valid_ct);
  } else {
    do {
      sample_uidx = next_set_ul_unsafe(pheno_nm2, sample_uidx);
      sample_uidx_stop = next_unset_ul(pheno_nm2, sample_uidx, unfiltered_sample_ct);
      sample_idx += sample_uidx_stop - sample_uidx;
      do {
        dzz = ((int32_t)IS_SET(pheno_c, sample_uidx));
	*dptr++ = dzz;
        dxx += dzz;
      } while ((++sample_uidx) < sample_uidx_stop);
    } while (sample_idx < sample_valid_ct);
    dyy = dxx;
  }
  if (dyy * ((double)((intptr_t)sample_valid_ct)) == dxx * dxx) {
    logerrprint("Warning: Skipping --lasso since phenotype is constant.\n");
    goto lasso_ret_1;
  }
  dzz = dxx / ((double)((intptr_t)sample_valid_ct)); // mean
  dyy = sqrt_n_recip * sqrt(((double)((intptr_t)(sample_valid_ct - 1))) / (dyy - dxx * dzz)); // 1/(stdev * sqrt(n))
  dptr = pheno_d_collapsed;
  for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
    *dptr = ((*dptr) - dzz) * dyy;
    dptr++;
  }
  fill_quatervec_55(sample_include2, sample_valid_ct);
  fill_ulong_zero(polymorphic_markers, unfiltered_marker_ctl);
  if ((chrom_info_ptr->mt_code != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->mt_code)) {
    hh_or_mt_exists |= NXMHH_EXISTS;
  }
  if (alloc_collapsed_haploid_filters(unfiltered_sample_ct, sample_valid_ct, hh_or_mt_exists, 1, pheno_nm2, sex_male, &sample_include2, &sample_male_include2)) {
    goto lasso_ret_NOMEM;
  }
  if (select_covars && select_covars_range_list_ptr->name_ct) {
    if (!covar_ct) {
      logerrprint("Error: No covariates loaded for --lasso-select-covars.\n");
      goto lasso_ret_INVALID_CMDLINE;
    }
    retval = string_range_list_to_bitfield_alloc(covar_names, covar_ct, max_covar_name_len, select_covars_range_list_ptr, &select_covars_bitfield, "lasso-select-covars", "--covar file");
    if (retval) {
      goto lasso_ret_1;
    }
  }
  uii = marker_ct + covar_ct;
  // maximum size of remaining allocations, for memory hog mode (col_ct is at
  // most marker_ct + covar_ct):
  // 1. data_arr: col_ct * sample_valid_ct * sizeof(double)
  // plus
  //   2a. xhat: col_ct * sizeof(double)
  //   2b. active_set: col_ctl * sizeof(intptr_t)
  // or
  //   3. prod_matrix: WARM_START_ITERS * WARM_START_ITERS * sizeof(double)
  // (whichever is larger)
  ullii = CACHEALIGN(((uint64_t)uii) * sizeof(double)) + CACHEALIGN(((uii + 7) / 8));
  // assumes WARM_START_ITERS is even
  if (rand_matrix) {
    uljj = (sample_valid_ct * WARM_START_ITERS) - 1;
    for (ulii = 0; ulii < uljj; ulii += 2) {
      rand_matrix[ulii] = rand_normal(&(rand_matrix[ulii + 1]));
    }
    if (ullii < CACHEALIGN(WARM_START_ITERS * WARM_START_ITERS * sizeof(double))) {
      ullii = CACHEALIGN(WARM_START_ITERS * WARM_START_ITERS * sizeof(double));
    }
  }
  ullii += CACHEALIGN(((uint64_t)uii) * sample_valid_ct * sizeof(double));
  // if (0) {
  if (ullii <= g_bigstack_left) {
    retval = lasso_bigmem(bedfile, bed_offset, marker_exclude, marker_ct, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, pheno_nm2, lasso_h2, lasso_minlambda, select_covars, select_covars_bitfield, pheno_d_collapsed, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, hh_or_mt_exists, sample_valid_ct, sample_include2, sample_male_include2, loadbuf_raw, loadbuf_collapsed, rand_matrix, misc_arr, residuals, polymorphic_markers, &polymorphic_marker_ct, &iter_tot, &xhat);
  } else {
    retval = lasso_smallmem(threads, bedfile, bed_offset, marker_exclude, marker_ct, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, pheno_nm2, lasso_h2, lasso_minlambda, select_covars, select_covars_bitfield, pheno_d_collapsed, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, hh_or_mt_exists, sample_valid_ct, sample_include2, sample_male_include2, loadbuf_raw, loadbuf_collapsed, rand_matrix, misc_arr, residuals, polymorphic_markers, &polymorphic_marker_ct, &iter_tot, &xhat);
    // retval = RET_NOMEM;
  }
  if (retval || (!polymorphic_marker_ct)) {
    goto lasso_ret_1;
  }
  memcpy(outname_end, ".lasso", 7);
  if (fopen_checked(outname, "w", &outfile)) {
    goto lasso_ret_OPEN_FAIL;
  }
  if (fputs_checked("CHR\tSNP\tA1\tEFFECT\n", outfile)) {
    goto lasso_ret_WRITE_FAIL;
  }
  tbuf[MAXLINELEN] = '\t';
  if (select_covars) {
    if (select_covars_bitfield) {
      marker_idx = covar_ct - popcount_longs(select_covars_bitfield, (covar_ct + (BITCT - 1)) / BITCT);
    } else {
      marker_idx = 0;
    }
    for (marker_uidx = 0; marker_idx < covar_ct; marker_uidx++, marker_idx++) {
      if (select_covars_bitfield) {
        next_set_unsafe_ck(select_covars_bitfield, &marker_uidx);
      }
      dxx = xhat[marker_idx];
      if ((!report_zeroes) && (dxx == 0)) {
	continue;
      }
      wptr = memcpya(tbuf, "COV\t", 4);
      wptr = strcpyax(wptr, &(covar_names[marker_uidx * max_covar_name_len]), '\t');
      wptr = memcpyl3a(wptr, "NA\t");
      wptr = double_g_writex(wptr, dxx, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto lasso_ret_WRITE_FAIL;
      }
    }
  }
  for (marker_uidx = 0, marker_idx = 0, marker_idx2 = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (marker_uidx >= chrom_end) {
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &uii, &min_ploidy_1);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = chrom_name_write(tbuf, chrom_info_ptr, uii);
      *wptr_start++ = '\t';
    }
    wptr = strcpyax(wptr_start, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
    if (IS_SET(polymorphic_markers, marker_uidx)) {
      dxx = xhat[covar_ct + (marker_idx2++)];
      if ((!report_zeroes) && (dxx == 0)) {
	continue;
      }
      wptr2 = double_g_writex(&(tbuf[MAXLINELEN + 1]), dxx, '\n');
    } else {
      if (!report_zeroes) {
	continue;
      }
      wptr2 = memcpyl3a(&(tbuf[MAXLINELEN + 1]), "NA\n");
    }
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto lasso_ret_WRITE_FAIL;
    }
    fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
    if (fwrite_checked(&(tbuf[MAXLINELEN]), (uintptr_t)(wptr2 - (&(tbuf[MAXLINELEN]))), outfile)) {
      goto lasso_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto lasso_ret_WRITE_FAIL;
  }
  putchar('\r');
  LOGPRINTFWW("--lasso report written to %s. Total iterations: %" PRIu64 ".\n", outname, iter_tot);

  while (0) {
  lasso_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  lasso_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  lasso_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  lasso_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 lasso_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}
