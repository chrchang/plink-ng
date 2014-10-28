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
  double* data_arr = (double*)wkspace_base; // marker-major
  double sqrt_n_recip = sqrt(1.0 / ((double)((intptr_t)sample_valid_ct)));
  double lambda_max = 0.0;
  double err_cur = 0.0;
  uint64_t iter_tot = 0;
  uintptr_t sample_valid_ctl2 = (sample_valid_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t polymorphic_marker_ct = 0;
  uintptr_t unselected_covar_ct = 0;
  uintptr_t final_mask = get_final_mask(sample_valid_ct);
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
  double lambda_min;
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
    logprint("Warning: Skipping --lasso since no polymorphic loci are present.\n");
    return 0;
  }
  col_ct = covar_ct + polymorphic_marker_ct;
  col_ctl = (col_ct + (BITCT - 1)) / BITCT;
  wkspace_shrink_top(data_arr, col_ct * sample_valid_ct * sizeof(double));
  xhat = (double*)wkspace_alloc(col_ct * sizeof(double));
  active_set = (uintptr_t*)wkspace_alloc(col_ctl * sizeof(intptr_t));
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
  sige = sqrt(1.0 - lasso_h2 + 1.0 / ((double)((intptr_t)sample_valid_ct)));
  zz = sige * sqrt_n_recip;
  if (rand_matrix) {
    prod_matrix = (double*)wkspace_alloc(col_ct * WARM_START_ITERS * sizeof(double));
    fputs("\r--lasso: Initializing warm start matrix...", stdout);
    fflush(stdout);
    col_major_matrix_multiply(WARM_START_ITERS, col_ct, sample_valid_ct, rand_matrix, data_arr, prod_matrix);
    for (ulii = 0; ulii < WARM_START_ITERS; ulii++) {
      dxx = 0.0;
      dptr = &(prod_matrix[ulii]);
      for (col_idx = 0; col_idx < col_ct; col_idx++) {
	dyy = fabs(dptr[col_idx * WARM_START_ITERS]);
	if (dyy > dxx) {
	  dxx = dyy;
	}
      }
      misc_arr[ulii] = dxx * zz;
    }
    lambda_min = destructive_get_dmedian(misc_arr, WARM_START_ITERS);
    logstr("--lasso:");
    LOGPRINTF(" using min lambda = %g.\n", lambda_min);
  } else {
    lambda_min = lasso_minlambda;
    if (lasso_minlambda >= lambda_max) {
      logprint("\nError: min lambda >= max lambda.\n");
      goto lasso_bigmem_ret_INVALID_CMDLINE;
    }
  }
  loghi = log(lambda_max);
  loglo = log(lambda_min);
  logdelta = (loghi - loglo) / (NLAMBDA - 1.0);
  for (col_idx = 0; col_idx < col_ct; col_idx++) {
    dxx = 0.0;
    dptr = &(data_arr[col_idx * sample_valid_ct]);
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
    logprint("Error: --lasso covariate is constant.\n");
  lasso_bigmem_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

/*
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
  uint64_t iter_tot = 0;
  uintptr_t sample_valid_ctl2 = (sample_valid_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t polymorphic_marker_ct = 0;
  uintptr_t final_mask = get_final_mask(sample_valid_ct);
  uint32_t chrom_fo_idx = 0xffffffffU; // exploit overflow
  uint32_t chrom_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t min_ploidy_1 = 0;
  int32_t retval = 0;
  double cur_mapping[4];
  double* data_window;
  double* xhat;
  double* prod_matrix;
  double* dptr;
  double* dptr2;
  uintptr_t* ulptr_end_init;
  uintptr_t* ulptr_end;
  uintptr_t* active_set;
  uintptr_t* ulptr;
  double sige;
  double lambda_min;
  double loghi;
  double loglo;
  double logdelta;
  double lambda;
  double xjold;
  double err_last;
  double dxx;
  double dyy;
  double dzz;
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
  uintptr_t sample_uidx;
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
    goto lasso_smallmem_ret_READ_FAIL;
  }
  ulptr_end_init = &(loadbuf_collapsed[sample_valid_ct / BITCT2]);
  fputs("Using memory-conserving LASSO implementation.\n", stdout);
  if (covar_ct) {
    if (wkspace_alloc_d_checked()) {
      goto lasso_smallmem_ret_NOMEM;
    }
    dptr = covar_data_arr;
    for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
      dxx = 0; // sum
      dyy = 0; // ssq
      dptr2 = &(covar_d[covar_idx]);
      for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_valid_ct; sample_uidx++, sample_idx++) {
	next_set_ul_unsafe_ck(covar_nm, &sample_uidx);
        dzz = dptr2[sample_uidx * covar_ct];
        dxx += dzz;
	dyy += dzz * dzz;
        *dptr++ = dzz;
      }
      if (dyy * ((double)sample_valid_ct) == dxx * dxx) {
        logprint("Error: --lasso covariate is constant.\n");
        goto lasso_smallmem_ret_INVALID_CMDLINE;
      }
      dzz = dxx / ((double)((intptr_t)sample_valid_ct));
      dyy = sqrt_n_recip * sqrt(((double)((intptr_t)(sample_valid_ct - 1))) / (dyy - dxx * dzz));
      dptr = &(data_arr[covar_idx * sample_valid_ct]);
      for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
        *dptr = ((*dptr) - dzz) * dyy;
        dptr++;
      }
    }
  }
  dptr = &(data_arr[covar_ct * sample_valid_ct]);
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
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
    if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_sample_ct, loadbuf_collapsed, sample_valid_ct, pheno_nm2, final_mask, IS_SET(marker_reverse, marker_uidx))) {
      goto lasso_smallmem_ret_READ_FAIL;
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
    logprint("Warning: Skipping --lasso since no polymorphic markers are present.\n");
    return 0;
  }
  col_ct = covar_ct + polymorphic_marker_ct;
  col_ctl = (col_ct + (BITCT - 1)) / BITCT;
  wkspace_shrink_top(data_arr, col_ct * sample_valid_ct * sizeof(double));
  xhat = (double*)wkspace_alloc(col_ct * sizeof(double));
  active_set = (uintptr_t*)wkspace_alloc(col_ctl * sizeof(intptr_t));
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
  sige = sqrt(1.0 - lasso_h2 + 1.0 / ((double)((intptr_t)sample_valid_ct)));
  zz = sige * sqrt_n_recip;
  if (rand_matrix) {
    prod_matrix = (double*)wkspace_alloc(col_ct * WARM_START_ITERS * sizeof(double));
    fputs("\r--lasso: Initializing warm start matrix...", stdout);
    fflush(stdout);
    col_major_matrix_multiply(WARM_START_ITERS, col_ct, sample_valid_ct, rand_matrix, data_arr, prod_matrix);
    for (ulii = 0; ulii < WARM_START_ITERS; ulii++) {
      dxx = 0.0;
      dptr = &(prod_matrix[ulii]);
      for (col_idx = 0; col_idx < col_ct; col_idx++) {
	dyy = fabs(dptr[col_idx * WARM_START_ITERS]);
	if (dyy > dxx) {
	  dxx = dyy;
	}
      }
      misc_arr[ulii] = dxx * zz;
    }
    lambda_min = destructive_get_dmedian(misc_arr, WARM_START_ITERS);
    logstr("--lasso:");
    LOGPRINTF(" using min lambda = %g.\n", lambda_min);
  } else {
    lambda_min = lasso_minlambda;
    if (lasso_minlambda >= lambda_max) {
      logprint("\nError: min lambda >= max lambda.\n");
      goto lasso_smallmem_ret_INVALID_CMDLINE;
    }
  }
  loghi = log(lambda_max);
  loglo = log(lambda_min);
  logdelta = (loghi - loglo) / (NLAMBDA - 1.0);
  for (col_idx = 0; col_idx < col_ct; col_idx++) {
    dxx = 0.0;
    dptr = &(data_arr[col_idx * sample_valid_ct]);
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
	if (col_uidx >= covar_ct) {
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
  lasso_smallmem_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  lasso_smallmem_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  lasso_smallmem_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}
*/

int32_t lasso(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t pheno_nm_ct, double lasso_h2, double lasso_minlambda, Range_list* select_covars_range_list_ptr, uint64_t misc_flags, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_male, uint32_t hh_or_mt_exists) {
  // Coordinate descent LASSO.  Based on a MATLAB script by Shashaank
  // Vattikuti.
  // Not yet multithreaded.  (Main loop is fairly tightly coupled, so getting
  // a performance benefit will be a bit tricky.)
  unsigned char* wkspace_mark = wkspace_base;
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
      logprint("Warning: Skipping --lasso since less than two phenotypes are present.\n");
    } else {
      logprint("Warning: Skipping --lasso since too many samples have missing covariates.\n");
    }
    goto lasso_ret_1;
  }
  if (sample_valid_ct == pheno_nm_ct) {
    pheno_nm2 = pheno_nm;
  } else {
    if (wkspace_alloc_ul_checked(&pheno_nm2, unfiltered_sample_ctl * sizeof(intptr_t))) {
      goto lasso_ret_NOMEM;
    }
    fill_ulong_zero(pheno_nm2, unfiltered_sample_ctl);
    for (sample_uidx = 0, sample_idx = 0; sample_idx < pheno_nm_ct; sample_uidx++, sample_idx++) {
      next_set_ul_unsafe_ck(pheno_nm, &sample_uidx);
      if (IS_SET(covar_nm, sample_idx)) {
        SET_BIT(pheno_nm2, sample_uidx);
      }
    }
  }
  sample_valid_ctv2 = 2 * ((sample_valid_ct + (BITCT - 1)) / BITCT);
  sqrt_n_recip = sqrt(1.0 / ((double)((intptr_t)sample_valid_ct)));
  if (wkspace_alloc_ul_checked(&sample_include2, sample_valid_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_sample_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf_collapsed, sample_valid_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&polymorphic_markers, unfiltered_marker_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&pheno_d_collapsed, sample_valid_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&residuals, sample_valid_ct * sizeof(double))) {
    goto lasso_ret_NOMEM;
  }
  if (lasso_minlambda == -1) {
    if (wkspace_alloc_d_checked(&rand_matrix, sample_valid_ct * WARM_START_ITERS * sizeof(double)) ||
        wkspace_alloc_d_checked(&misc_arr, WARM_START_ITERS * sizeof(double))) {
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
    logprint("Warning: Skipping --lasso since phenotype is constant.\n");
    goto lasso_ret_1;
  }
  dzz = dxx / ((double)((intptr_t)sample_valid_ct)); // mean
  dyy = sqrt_n_recip * sqrt(((double)((intptr_t)(sample_valid_ct - 1))) / (dyy - dxx * dzz)); // 1/(stdev * sqrt(n))
  dptr = pheno_d_collapsed;
  for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
    *dptr = ((*dptr) - dzz) * dyy;
    dptr++;
  }
  fill_vec_55(sample_include2, sample_valid_ct);
  fill_ulong_zero(polymorphic_markers, unfiltered_marker_ctl);
  if ((chrom_info_ptr->mt_code != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->mt_code)) {
    hh_or_mt_exists |= NXMHH_EXISTS;
  }
  if (alloc_collapsed_haploid_filters(unfiltered_sample_ct, sample_valid_ct, hh_or_mt_exists, 1, pheno_nm2, sex_male, &sample_include2, &sample_male_include2)) {
    goto lasso_ret_NOMEM;
  }
  if (select_covars && select_covars_range_list_ptr->name_ct) {
    if (!covar_ct) {
      logprint("Error: No covariates loaded for --lasso-select-covars.\n");
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
  // 2. prod_matrix: col_ct * WARM_START_ITERS * sizeof(double)
  // 3. xhat: col_ct * sizeof(double)
  // 4. active_set: col_ctl * sizeof(intptr_t)
  ullii = CACHEALIGN(((uint64_t)uii) * sample_valid_ct * sizeof(double)) + CACHEALIGN(((uint64_t)uii) * sizeof(double)) + CACHEALIGN(((uii + 7) / 8));
  // assumes WARM_START_ITERS is even
  if (rand_matrix) {
    uljj = (sample_valid_ct * WARM_START_ITERS) - 1;
    for (ulii = 0; ulii < uljj; ulii += 2) {
      rand_matrix[ulii] = rand_normal(&(rand_matrix[ulii + 1]));
    }
    ullii += CACHEALIGN(((uint64_t)uii) * WARM_START_ITERS * sizeof(double)) + CACHEALIGN(((uint64_t)uii) * sizeof(double)) + CACHEALIGN(((uii + 7) / 8));
  }
  if (ullii <= wkspace_left) {
    retval = lasso_bigmem(bedfile, bed_offset, marker_exclude, marker_ct, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, pheno_nm2, lasso_h2, lasso_minlambda, select_covars, select_covars_bitfield, pheno_d_collapsed, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, hh_or_mt_exists, sample_valid_ct, sample_include2, sample_male_include2, loadbuf_raw, loadbuf_collapsed, rand_matrix, misc_arr, residuals, polymorphic_markers, &polymorphic_marker_ct, &iter_tot, &xhat);
  } else {
    // retval = lasso_smallmem(threads, bedfile, bed_offset, marker_exclude, marker_ct, marker_reverse, chrom_info_ptr, unfiltered_sample_ct, pheno_nm2, lasso_h2, lasso_minlambda, select_covars, select_covars_bitfield, pheno_d_collapsed, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, hh_or_mt_exists, sample_valid_ct, sample_include2, sample_male_include2, loadbuf_raw, loadbuf_collapsed, rand_matrix, misc_arr, residuals, polymorphic_markers, &polymorphic_marker_ct, &iter_tot, &xhat);
    retval = RET_NOMEM;
  }
  if (retval) {
    goto lasso_ret_1;
  }
  memcpy(outname_end, ".lasso", 7);
  if (fopen_checked(&outfile, outname, "w")) {
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
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
