#include "wdist_lasso.h"
#include "wdist_matrix.h"

#define WARM_START_ITERS 1000
#define NLAMBDA 100
#define DELTA_THRESHOLD 0.0001

int32_t lasso(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t pheno_nm_ct, double lasso_h2, uint32_t report_zeroes, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* sex_male, uint32_t hh_exists) {
  // Coordinate descent LASSO.  Based on a MATLAB script by Shashaank
  // Vattikuti.
  // Not yet multithreaded.  (Main loop is fairly tightly coupled, so getting
  // a performance benefit will be a bit tricky.)
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctv2 = 2 * ((unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t pheno_nm_ctl2 = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t pheno_nm_ctv2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t polymorphic_marker_ct = 0;
  uint64_t iter_tot = 0;
  double lambda_max = 0.0;
  double err_cur = 0.0;
  uintptr_t* indiv_male_include2 = NULL;
  char* wptr_start = NULL;
  uint32_t chrom_fo_idx = 0xffffffffU; // exploit overflow
  uint32_t chrom_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t is_haploid = 0;
  int32_t retval = 0;
  double cur_mapping[4];
  double* data_arr;
  double* pheno_d_collapsed;
  double* rand_matrix;
  double* misc_arr;
  double* prod_matrix;
  double* xhat;
  double* residuals;
  double* dptr;
  double* dptr2;
  double* dptr3;
  uintptr_t* indiv_include2;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_collapsed;
  uintptr_t* polymorphic_markers;
  uintptr_t* ulptr;
  uintptr_t* ulptr_end;
  uintptr_t* ulptr_end_init;
  uintptr_t* active_set;
  char* wptr;
  double sqrt_n_recip;
  double sige;
  double dxx;
  double dyy;
  double dzz;
  double zz;
  double lambda_min;
  double loghi;
  double loglo;
  double logdelta;
  double lambda;
  double xjold;
  double err_last;
  uintptr_t polymorphic_marker_ctl;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t iter;
  uintptr_t col_nz_ct;
  uintptr_t col_to_z;
  uintptr_t col_uidx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t indiv_idx_stop;
  uint32_t lambi;
  uint32_t marker_uidx;
  uint32_t homrar_ct;
  uint32_t missing_ct;
  uint32_t het_ct;
  uint32_t homset_ct;
  uint32_t uii;
  if (pheno_nm_ct < 2) {
    logprint("Warning: Skipping --lasso since less than two phenotypes are present.\n");
    goto lasso_ret_1;
  }
  sqrt_n_recip = sqrt(1.0 / ((double)((intptr_t)pheno_nm_ct)));
  if (wkspace_alloc_ul_checked(&indiv_include2, pheno_nm_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf_collapsed, pheno_nm_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&polymorphic_markers, unfiltered_marker_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&pheno_d_collapsed, pheno_nm_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&rand_matrix, pheno_nm_ct * WARM_START_ITERS * sizeof(double)) ||
      wkspace_alloc_d_checked(&misc_arr, WARM_START_ITERS * sizeof(double)) ||
      wkspace_alloc_d_checked(&residuals, pheno_nm_ct * sizeof(double))) {
    goto lasso_ret_NOMEM;
  }
  dxx = 0; // sum
  dyy = 0; // ssq
  dptr = pheno_d_collapsed;
  indiv_uidx = 0;
  indiv_idx = 0;
  do {
    indiv_uidx = next_set_ul_unsafe(pheno_nm, indiv_uidx);
    indiv_uidx_stop = next_unset_ul(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    dptr2 = &(pheno_d[indiv_uidx]);
    indiv_uidx = indiv_uidx_stop;
    dptr3 = &(pheno_d[indiv_uidx_stop]);
    do {
      dzz = *dptr2++;
      *dptr++ = dzz;
      dxx += dzz;
      dyy += dzz * dzz;
    } while (dptr2 < dptr3);
  } while (indiv_idx < pheno_nm_ct);
  if (dyy * ((double)((intptr_t)pheno_nm_ct)) == dxx * dxx) {
    logprint("Warning: Skipping --lasso since phenotype is constant.\n");
    goto lasso_ret_1;
  }
  dzz = dxx / ((double)((intptr_t)pheno_nm_ct)); // mean
  dyy = sqrt_n_recip * sqrt(((double)((intptr_t)(pheno_nm_ct - 1))) / (dyy - dxx * dzz)); // 1/(stdev * sqrt(n))
  dptr = pheno_d_collapsed;
  for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
    *dptr = ((*dptr) - dzz) * dyy;
    dptr++;
  }
  fill_vec_55(indiv_include2, pheno_nm_ct);
  fill_ulong_zero(polymorphic_markers, unfiltered_marker_ctl);
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, pheno_nm_ct, hh_exists, 1, pheno_nm, sex_male, &indiv_include2, &indiv_male_include2)) {
    goto lasso_ret_NOMEM;
  }
  // this is the big one.  todo: write a fallback implementation which does not
  // require the entire matrix to be in memory at once.
  if (wkspace_alloc_d_checked(&data_arr, marker_ct * pheno_nm_ct * sizeof(double))) {
    goto lasso_ret_NOMEM;
  }

  cur_mapping[1] = 0; // missing
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto lasso_ret_READ_FAIL;
  }
  dptr = data_arr;
  ulptr_end_init = &(loadbuf_collapsed[pheno_nm_ct / BITCT2]);
  printf("--lasso: Populating data matrix...");
  fflush(stdout);
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto lasso_ret_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
    }
    if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_collapsed, pheno_nm_ct, pheno_nm, IS_SET(marker_reverse, marker_uidx))) {
      goto lasso_ret_READ_FAIL;
    }
    if (is_haploid) {
      haploid_fix(hh_exists, indiv_include2, indiv_male_include2, pheno_nm_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
    }
    vec_3freq(pheno_nm_ctl2, loadbuf_collapsed, indiv_include2, &missing_ct, &het_ct, &homset_ct);
    uii = pheno_nm_ct - missing_ct;
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
      indiv_idx = 0;
      indiv_idx_stop = BITCT2;
      ulptr = loadbuf_collapsed;
      ulptr_end = ulptr_end_init;
      while (1) {
        while (ulptr < ulptr_end) {
          cur_word = *ulptr++;
          for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
            cur_genotype = cur_word & 3;
            *dptr++ = cur_mapping[cur_genotype];
	  }
          indiv_idx_stop += BITCT2;
	}
        if (indiv_idx == pheno_nm_ct) {
	  break;
	}
        ulptr_end++;
        indiv_idx_stop = pheno_nm_ct;
      }
      polymorphic_marker_ct++;
    }
  }
  if (!polymorphic_marker_ct) {
    logprint("Warning: Skipping --lasso since no polymorphic markers are present.\n");
    goto lasso_ret_1;
  }
  polymorphic_marker_ctl = (polymorphic_marker_ct + (BITCT - 1)) / BITCT;
  wkspace_reset((unsigned char*)data_arr);
  data_arr = (double*)wkspace_alloc(polymorphic_marker_ct * pheno_nm_ct * sizeof(double));
  if (wkspace_alloc_d_checked(&prod_matrix, polymorphic_marker_ct * WARM_START_ITERS * sizeof(double)) ||
      wkspace_alloc_d_checked(&xhat, polymorphic_marker_ct * sizeof(double)) ||
      wkspace_alloc_ul_checked(&active_set, polymorphic_marker_ctl * sizeof(intptr_t))) {
    goto lasso_ret_NOMEM;
  }
  dptr = data_arr;
  for (marker_idx = 0; marker_idx < polymorphic_marker_ct; marker_idx++) {
    dxx = 0;
    dptr2 = pheno_d_collapsed;
    for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
      dxx += (*dptr++) * (*dptr2++);
    }
    dxx = fabs(dxx);
    if (dxx > lambda_max) {
      lambda_max = dxx;
    }
  }
  sige = sqrt(1.0 - lasso_h2);
  zz = (sige + sqrt_n_recip) * sqrt_n_recip;
  uljj = (pheno_nm_ct * WARM_START_ITERS) - 1;
  printf("\r--lasso: Initializing warm start matrix...");
  fflush(stdout);
  // assumes WARM_START_ITERS is even
  for (ulii = 0; ulii < uljj; ulii += 2) {
    rand_matrix[ulii] = rand_normal(&(rand_matrix[ulii + 1]));
  }
  col_major_matrix_multiply(WARM_START_ITERS, polymorphic_marker_ct, pheno_nm_ct, rand_matrix, data_arr, prod_matrix);
  for (ulii = 0; ulii < WARM_START_ITERS; ulii++) {
    dxx = 0.0;
    dptr = &(prod_matrix[ulii]);
    for (marker_idx = 0; marker_idx < polymorphic_marker_ct; marker_idx++) {
      dyy = fabs(dptr[marker_idx * WARM_START_ITERS]);
      if (dyy > dxx) {
        dxx = dyy;
      }
    }
    misc_arr[ulii] = dxx * zz;
  }
  lambda_min = destructive_get_dmedian(misc_arr, WARM_START_ITERS);
  loghi = log(lambda_max);
  loglo = log(lambda_min);
  logdelta = (loghi - loglo) / (NLAMBDA - 1.0);
  for (marker_idx = 0; marker_idx < polymorphic_marker_ct; marker_idx++) {
    dxx = 0.0;
    dptr = &(data_arr[marker_idx * pheno_nm_ct]);
    dptr2 = pheno_d_collapsed;
    for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
      dxx += (*dptr++) * (*dptr2++);
    }
    xhat[marker_idx] = dxx;
  }
  fflush(stdout);
  for (lambi = 0; lambi < NLAMBDA; lambi++) {
    printf("\r--lasso: Executing coordinate descent... %u%%", lambi);
    fflush(stdout);
    lambda = exp(loghi - logdelta * ((double)((int32_t)lambi)));
    memcpy(residuals, pheno_d_collapsed, pheno_nm_ct * sizeof(double));
    for (marker_idx = 0; marker_idx < polymorphic_marker_ct; marker_idx++) {
      dptr = residuals;
      dptr2 = &(data_arr[marker_idx * pheno_nm_ct]);
      dxx = -xhat[marker_idx];
      for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
        *dptr += (*dptr2++) * dxx;
        dptr++;
      }
    }
    iter = 0;
    fill_ulong_one(active_set, polymorphic_marker_ctl);
    ulii = polymorphic_marker_ct & (BITCT - 1);
    if (ulii) {
      active_set[polymorphic_marker_ctl - 1] = (~ZEROLU) >> (BITCT - ulii);
    }
    col_nz_ct = polymorphic_marker_ct;
    while (1) {
      col_uidx = 0;
      col_to_z = 0;
      for (marker_idx = 0; marker_idx < col_nz_ct; marker_idx++, col_uidx++) {
        col_uidx = next_set_unsafe(active_set, col_uidx);
        xjold = xhat[col_uidx];
        dxx = xjold;
        dptr = &(data_arr[col_uidx * pheno_nm_ct]);
        dptr2 = residuals;
        for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
          dxx += (*dptr++) * (*dptr2++);
	}
        if (dxx > 0.0) {
	  dxx = MAXV(dxx - lambda, 0.0);
	} else {
          dxx = MINV(dxx + lambda, 0.0);
	}
        xhat[col_uidx] = dxx;
        if (dxx == 0.0) {
          CLEAR_BIT(active_set, col_uidx);
	  col_to_z++;
	}
        dptr = residuals;
        dptr2 = &(data_arr[col_uidx * pheno_nm_ct]);
        dxx -= xjold;
        for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
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
      for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
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
  memcpy(outname_end, ".lasso", 7);
  if (fopen_checked(&outfile, outname, "w")) {
    goto lasso_ret_OPEN_FAIL;
  }
  if (fputs_checked("CHR\tSNP\tA1\tEFFECT\n", outfile)) {
    goto lasso_ret_WRITE_FAIL;
  }
  chrom_fo_idx = 0xffffffffU;
  chrom_end = 0;
  for (marker_uidx = 0, marker_idx = 0, marker_idx2 = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (marker_uidx >= chrom_end) {
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = chrom_name_write(tbuf, chrom_info_ptr, uii, zero_extra_chroms);
      *wptr_start++ = '\t';
    }
    wptr = strcpyax(wptr_start, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
    if (max_marker_allele_len == 1) {
      *wptr++ = marker_alleles[2 * marker_uidx];
    } else {
      wptr = strcpya(wptr, &(marker_alleles[2 * marker_uidx * max_marker_allele_len]));
    }
    *wptr++ = '\t';
    if (IS_SET(polymorphic_markers, marker_uidx)) {
      dxx = xhat[marker_idx2++];
      if ((!report_zeroes) && (dxx == 0)) {
	continue;
      }
      wptr = double_g_writex(wptr, dxx, '\n');
    } else {
      if (!report_zeroes) {
	continue;
      }
      wptr = memcpyl3a(wptr, "NA\n");
    }
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto lasso_ret_WRITE_FAIL;
    }    
  }
  if (fclose_null(&outfile)) {
    goto lasso_ret_WRITE_FAIL;
  }
  putchar('\r');
  sprintf(logbuf, "--lasso report written to %s.  Total iterations: %llu.\n", outname, iter_tot);
  logprintb();

  while (0) {
  lasso_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  lasso_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  lasso_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  lasso_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 lasso_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
