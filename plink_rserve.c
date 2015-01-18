// Rconnection.cc is part of a continuously updated package and uses
// C++-specific language constructs, so we omit this module when compiling with
// gcc instead of g++.
#if defined __cplusplus && !defined _WIN32

#include "plink_common.h"
#include "plink_cluster.h"

#define MAIN
#define SOCK_ERRORS
#include "sisocks.h"
#include "Rconnection.h"

#define RPLUGIN_BLOCK_SIZE 100

int32_t rserve_call(char* rplugin_fname, uint32_t rplugin_port, uint32_t rplugin_debug, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uintptr_t* pheno_c, double* pheno_d, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uintptr_t covar_ct, double* covar_d, char* outname, char* outname_end) {
  // See PLINK 1.07 r.cpp.
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  char* wkspace_end = (char*)(&(wkspace_base[wkspace_left]));
  int32_t* geno_int_buf = NULL;
  Rinteger* r_n = NULL;
  Rinteger* r_s = NULL;
  Rdouble* r_p = NULL;
  Rdouble* r_cov = NULL;
  Rconnection* rc = NULL;
  char* chrom_name_ptr = NULL;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = (unfiltered_sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t pheno_nm_ctl2 = (pheno_nm_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t line_idx = 0;
  uintptr_t final_mask = get_final_mask(pheno_nm_ct);
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  uint32_t chrom_end = 0;
  uint32_t chrom_name_len = 0;
  uint32_t pct = 0;
  int32_t retval = 0;
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_LEN];
  uintptr_t marker_uidx;
  uintptr_t marker_uidx_base;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* ulptr;
  Rinteger* r_l;
  Rinteger* r_g;
  Rdouble* r_data;
  double* pheno_d2;
  double* dptr;
  char* inbuf_start;
  char* inbuf_end;
  char* bufptr;
  int32_t* sample_to_cluster;
  int32_t* iptr;
  double dxx;
  uint32_t marker_idx_base;
  uint32_t block_size;
  uint32_t block_offset;
  uint32_t loop_end;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t cur_data_len;
  uint32_t uii;
  int32_t ii;
  if (fopen_checked(&infile, rplugin_fname, "r")) {
    goto rserve_call_ret_OPEN_FAIL;
  }
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_sample_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, pheno_nm_ctl2 * RPLUGIN_BLOCK_SIZE * sizeof(intptr_t))) {
    goto rserve_call_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  for (ulii = 1; ulii <= RPLUGIN_BLOCK_SIZE; ulii++) {
    loadbuf[ulii * pheno_nm_ctl2 - 1] = 0;
  }
  inbuf_start = (char*)wkspace_base;
  inbuf_end = inbuf_start;
  while (1) {
    if ((uintptr_t)(wkspace_end - inbuf_start) < MAXLINELEN) {
      goto rserve_call_ret_NOMEM;
    }
    inbuf_end[MAXLINELEN - 1] = ' ';
    if (!fgets(inbuf_end, MAXLINELEN, infile)) {
      break;
    }
    line_idx++;
    if (!(inbuf_end[MAXLINELEN - 1])) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --R file is pathologically long.\n", line_idx);
      goto rserve_call_ret_INVALID_FORMAT_2;
    }
    uii = strlen(inbuf_end);
    // standardize line ending
    if (uii) {
      if (inbuf_end[uii - 1] == '\n') {
	if ((uii >= 2) && (inbuf_end[uii - 2] == '\r')) {
	  inbuf_end[uii - 2] = '\n';
	  inbuf_end = &(inbuf_end[uii - 1]);
	} else {
	  inbuf_end = &(inbuf_end[uii]);
	}
      } else {
	inbuf_end[uii] = '\n';
	inbuf_end = &(inbuf_end[uii + 1]);
      }
    } else {
      *inbuf_end++ = '\n';
    }
  }
  if (fclose_null(&infile)) {
    goto rserve_call_ret_READ_FAIL;
  }
  if (inbuf_end == inbuf_start) {
    logprint("Error: Empty --R file.\n");
    goto rserve_call_ret_INVALID_FORMAT;
  }
  *inbuf_end = '\0';
  wkspace_alloc(1 + ((uintptr_t)(inbuf_end - inbuf_start)));
  if (pheno_c) {
    if (wkspace_alloc_d_checked(&pheno_d2, pheno_nm_ct * sizeof(double))) {
      goto rserve_call_ret_NOMEM;
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < pheno_nm_ct; sample_uidx++, sample_idx++) {
      next_set_unsafe_ck(pheno_nm, &sample_uidx);
      pheno_d2[sample_idx] = (double)((int32_t)(1 + is_set(pheno_c, sample_uidx)));
    }
  } else if (pheno_nm_ct != unfiltered_sample_ct) {
    pheno_d2 = (double*)alloc_and_init_collapsed_arr_incl((char*)pheno_d, sizeof(double), unfiltered_sample_ct, pheno_nm, pheno_nm_ct, 1);
    if (!pheno_d2) {
      goto rserve_call_ret_NOMEM;
    }
  } else {
    pheno_d2 = pheno_d;
  }
  if (cluster_ct) {
    if (wkspace_alloc_i_checked(&sample_to_cluster, unfiltered_sample_ct * sizeof(int32_t))) {
      goto rserve_call_ret_NOMEM;
    }
    fill_int_one(sample_to_cluster, pheno_nm_ct);
    fill_unfiltered_sample_to_cluster(unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, (uint32_t*)sample_to_cluster);
    inplace_collapse_uint32_incl((uint32_t*)sample_to_cluster, unfiltered_sample_ct, pheno_nm, pheno_nm_ct);
    wkspace_shrink_top(sample_to_cluster, pheno_nm_ct * sizeof(int32_t));
  } else {
    if (wkspace_alloc_i_checked(&sample_to_cluster, pheno_nm_ct * sizeof(int32_t))) {
      goto rserve_call_ret_NOMEM;
    }
    fill_int_zero(sample_to_cluster, pheno_nm_ct);
  }
  if (!rplugin_debug) {
    if (wkspace_alloc_i_checked(&geno_int_buf, RPLUGIN_BLOCK_SIZE * pheno_nm_ct * sizeof(int32_t))) {
      goto rserve_call_ret_NOMEM;
    }
    rc = new Rconnection("127.0.0.1", rplugin_port);
    ii = rc->connect();
    if (ii) {
      sockerrorchecks(tbuf, 128, -1);
      LOGPRINTFWW("Error: Unable to connect (code %d: %s).\n", ii, tbuf);
      goto rserve_call_ret_NETWORK;
    }
    rc->eval("options(echo=F)");
    
    r_s = new Rinteger(sample_to_cluster, pheno_nm_ct);
    r_p = new Rdouble(pheno_d2, pheno_nm_ct);
    r_n = new Rinteger((int32_t*)(&pheno_nm_ct), 1);
    rc->assign("n", r_n);
    rc->assign("PHENO", r_p);
    rc->assign("CLUSTER", r_s);
    rc->eval("CLUSTER[CLUSTER==-1] <- NA");
    if (covar_ct) {
      r_cov = new Rdouble(covar_d, pheno_nm_ct * covar_ct);
      rc->assign("c", r_cov);
      rc->eval("COVAR<-matrix(c,nrow=n,byrow=T)");
    } else {
      rc->eval("COVAR<-NA");
    }
    memcpy(outname_end, ".auto.R", 8);
  } else {
    memcpy(outname_end, ".debug.R", 9);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto rserve_call_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("--R%s: writing to %s ... ", rplugin_debug? " debug" : "", outname);
  fputs("0%", stdout);
  fflush(stdout);
  loop_end = marker_ct / 100;
  if (rplugin_debug) {
    bufptr = memcpya(tbuf, "n <- ", 5);
    bufptr = uint32_write(bufptr, pheno_nm_ct);
    bufptr = memcpya(bufptr, "\nPHENO <- c( ", 13);
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      goto rserve_call_ret_WRITE_FAIL;
    }
    for (sample_idx = 0; sample_idx < pheno_nm_ct - 1; sample_idx++) {
      bufptr = double_g_write(tbuf, pheno_d2[sample_idx]);
      bufptr = memcpya(bufptr, ", ", 2);
      fwrite(tbuf, 1, (uintptr_t)(bufptr - tbuf), outfile);
    }
    bufptr = double_g_write(tbuf, pheno_d2[sample_idx]);
    bufptr = memcpya(bufptr, " ) \n", 4);
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      goto rserve_call_ret_WRITE_FAIL;
    }
    if (covar_ct) {
      fputs("c <- c( ", outfile);
      uljj = pheno_nm_ct * covar_ct - 1;
      for (ulii = 0; ulii < uljj; ulii++) {
	bufptr = double_g_write(tbuf, covar_d[ulii]);
	bufptr = memcpya(bufptr, ", ", 2);
	fwrite(tbuf, 1, (uintptr_t)(bufptr - tbuf), outfile);
      }
      bufptr = double_g_write(tbuf, covar_d[ulii]);
      fwrite(tbuf, 1, (uintptr_t)(bufptr - tbuf), outfile);
      fputs(" ) \nCOVAR <- matrix( c , nrow = n , byrow=T)\n", outfile);
    } else {
      fputs("COVAR <- matrix( NA , nrow = n , ncol = 0 , byrow = T)\n", outfile);
    }
    fputs("CLUSTER <- c( ", outfile);
    for (sample_idx = 0; sample_idx < pheno_nm_ct - 1; sample_idx++) {
      bufptr = int32_write(tbuf, sample_to_cluster[sample_idx]);
      bufptr = memcpya(bufptr, ", ", 2);
      fwrite(tbuf, 1, (uintptr_t)(bufptr - tbuf), outfile);
    }
    bufptr = int32_write(tbuf, sample_to_cluster[sample_idx]);
    bufptr = memcpya(bufptr, " ) \n", 4);
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      goto rserve_call_ret_WRITE_FAIL;
    }
  }
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto rserve_call_ret_READ_FAIL;
  }
  for (marker_idx_base = 0; marker_idx_base < marker_ct;) {
    block_size = marker_ct - marker_idx_base;
    if (block_size > RPLUGIN_BLOCK_SIZE) {
      block_size = RPLUGIN_BLOCK_SIZE;
    }
    ulptr = loadbuf;
    marker_uidx_base = marker_uidx;
    for (block_offset = 0; block_offset < block_size; marker_uidx++, block_offset++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto rserve_call_ret_READ_FAIL;
	}
      }
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_sample_ct, ulptr, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx))) {
	goto rserve_call_ret_READ_FAIL;
      }
      // 0 -> 3
      // 1 -> 0
      // 2 -> 2
      // 3 -> 1
      // i.e.
      //   bottom bit unset -> top bit set
      //   top bit same as bottom bit -> bottom bit set
      for (uii = 0; uii < pheno_nm_ctl2; uii++) {
        ulii = *ulptr;
	uljj = (~ulii) & FIVEMASK;
	ulii = (ulii >> 1) & FIVEMASK;
	*ulptr++ = (uljj << 1) | (ulii ^ uljj);
      }
      // overflow bits on last word don't matter
    }
    if (!rplugin_debug) {
      // transpose, subtract 1
      iptr = geno_int_buf;
      for (sample_idx = 0; sample_idx < pheno_nm_ct; sample_idx++) {
	ulptr = &(loadbuf[sample_idx / BITCT2]);
	uii = 2 * (sample_idx & (BITCT2 - 1));
	for (block_offset = 0; block_offset < block_size; block_offset++, ulptr = &(ulptr[pheno_nm_ctl2])) {
	  ulii = ((*ulptr) >> uii) & 3;
	  *iptr++ = ((int32_t)((uint32_t)ulii)) - 1;
	}
      }
      r_l = new Rinteger((int32_t*)(&block_size), 1);
      r_g = new Rinteger(geno_int_buf, block_size * pheno_nm_ct);
      rc->assign("l", r_l);
      rc->assign("g", r_g);
      rc->eval("GENO<-matrix(g,nrow=n,byrow=T)");
      rc->eval("GENO[GENO==-1] <- NA");
      delete r_l;
      delete r_g;
      rc->eval(inbuf_start);
      r_data = (Rdouble*)rc->eval("Rplink(PHENO,GENO,CLUSTER,COVAR)");
      if (r_data) {
	dptr = r_data->doubleArray();
	for (marker_uidx = marker_uidx_base, block_offset = 0; block_offset < block_size; marker_uidx++, block_offset++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  if (marker_uidx >= chrom_end) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1];
	    } while (marker_uidx >= chrom_end);
	    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    chrom_name_ptr = chrom_name_buf5w4write(chrom_name_buf, chrom_info_ptr, uii, &chrom_name_len);
	  }
	  bufptr = memcpyax(tbuf, chrom_name_ptr, chrom_name_len, ' ');
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = uint32_writew10x(bufptr, marker_pos[marker_uidx], ' ');
	  if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	    goto rserve_call_ret_WRITE_FAIL;
	  }
	  fputs_w4(marker_allele_ptrs[2 * marker_uidx], outfile);
	  tbuf[0] = ' ';
	  bufptr = &(tbuf[1]);
	  cur_data_len = (int32_t)(*dptr++);
	  for (uii = 0; uii < cur_data_len; uii++) {
	    dxx = *dptr++;
	    if (realnum(dxx)) {
	      bufptr = double_g_write(bufptr, dxx);
	    } else {
	      bufptr = memcpya(bufptr, "NA", 2);
	    }
	    *bufptr++ = '\t';
	    if (bufptr > &(tbuf[MAXLINELEN])) {
	      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
		goto rserve_call_ret_WRITE_FAIL;
	      }
	      bufptr = tbuf;
	    }
	  }
	  *bufptr++ = '\n';
	  if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	    goto rserve_call_ret_WRITE_FAIL;
	  }
	}
	delete r_data;
      } else {
	for (marker_uidx = marker_uidx_base, block_offset = 0; block_offset < block_size; marker_uidx++, block_offset++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  if (marker_uidx >= chrom_end) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1];
	    } while (marker_uidx >= chrom_end);
	    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    chrom_name_ptr = chrom_name_buf5w4write(chrom_name_buf, chrom_info_ptr, uii, &chrom_name_len);
	  }
	  bufptr = memcpyax(tbuf, chrom_name_ptr, chrom_name_len, ' ');
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = uint32_writew10x(bufptr, marker_pos[marker_uidx], ' ');
	  if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	    goto rserve_call_ret_WRITE_FAIL;
	  }
	  fputs_w4(marker_allele_ptrs[2 * marker_uidx], outfile);
	  fputs(" NA\n", outfile);
	}
      }
    } else {
      bufptr = memcpya(tbuf, "l <- ", 5);
      bufptr = uint32_write(bufptr, block_size);
      bufptr = memcpya(bufptr, "\ng <- c( ", 9);
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	goto rserve_call_ret_WRITE_FAIL;
      }
      block_offset = 0;
      sample_idx = 0;
      while (1) {
        bufptr = tbuf;
	ulptr = &(loadbuf[sample_idx / BITCT2]);
	uii = 2 * (sample_idx & (BITCT2 - 1));
	for (block_offset = 0; block_offset < block_size; block_offset++, ulptr = &(ulptr[pheno_nm_ctl2])) {
	  ulii = ((*ulptr) >> uii) & 3;
	  if (!ulii) {
	    bufptr = memcpya(bufptr, "-1", 2);
	  } else {
	    // '/' = ascii 47
	    *bufptr++ = (unsigned char)(47 + ulii);
	  }
	  bufptr = memcpya(bufptr, ", ", 2);
	}
	if (++sample_idx == pheno_nm_ct) {
	  break;
	}
	if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	  goto rserve_call_ret_WRITE_FAIL;
	}
      }
      bufptr = memcpya(&(bufptr[-2]), " ) \nGENO <- matrix( g , nrow = n ,byrow=T)\nGENO[GENO == -1 ] <- NA \n\n\n", 70);
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	goto rserve_call_ret_WRITE_FAIL;
      }
      if (fwrite_checked(inbuf_start, inbuf_end - inbuf_start, outfile)) {
	goto rserve_call_ret_WRITE_FAIL;
      }
      putc('\n', outfile);
    }
    marker_idx_base += block_size;
    if (marker_idx_base >= loop_end) {
      if (marker_idx_base < marker_ct) {
	if (pct >= 10) {
	  putchar('\b');
	}
	pct = (marker_idx_base * 100LLU) / marker_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = (((uint64_t)pct + 1LLU) * marker_ct) / 100;
      }
    }
  }

  if (pct >= 10) {
    putchar('\b');
  }
  fputs("\b\b", stdout);
  logprint("done.\n");
  while (0) {
  rserve_call_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  rserve_call_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  rserve_call_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  rserve_call_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  rserve_call_ret_INVALID_FORMAT_2:
    logprintb();
  rserve_call_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  rserve_call_ret_NETWORK:
    retval = RET_NETWORK;
    break;
  }
  fclose_cond(infile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  delete r_n;
  delete r_p;
  delete r_s;
  delete r_cov;
  return retval;
}

#endif // __cplusplus, !_WIN32
