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

int32_t rserve_call(char* rplugin_fname, uint32_t rplugin_port, uint32_t rplugin_debug, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uintptr_t* pheno_c, double* pheno_d, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uintptr_t covar_ct, double* covar_d, char* outname, char* outname_end) {
  // See PLINK 1.07 r.cpp.
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  char* wkspace_end = (char*)(&(wkspace_base[wkspace_left]));
  char* inbuf_start = (char*)wkspace_base;
  char* inbuf_ptr = inbuf_start;
  Rinteger* r_n = NULL;
  Rinteger* r_s = NULL;
  Rdouble* r_p = NULL;
  Rdouble* r_cov = NULL;
  Rconnection* rc = NULL;
  uintptr_t line_idx = 0;
  int32_t retval = 0;
  uintptr_t ulii;
  uintptr_t uljj;
  double* pheno_d2;
  char* bufptr;
  int32_t* sample_to_cluster;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t uii;
  int32_t ii;
  if (fopen_checked(&infile, rplugin_fname, "r")) {
    goto rserve_call_ret_OPEN_FAIL;
  }
  while (1) {
    if ((uintptr_t)(wkspace_end - inbuf_start) < MAXLINELEN) {
      goto rserve_call_ret_NOMEM;
    }
    inbuf_ptr[MAXLINELEN - 1] = ' ';
    if (!fgets(inbuf_ptr, MAXLINELEN, infile)) {
      break;
    }
    line_idx++;
    if (!(inbuf_ptr[MAXLINELEN - 1])) {
      sprintf(logbuf, "Error: Line %" PRIuPTR " of --R file is pathologically long.\n", line_idx);
      goto rserve_call_ret_INVALID_FORMAT_2;
    }
    uii = strlen(inbuf_ptr);
    // standardize line ending
    if (uii) {
      if (inbuf_ptr[uii - 1] == '\n') {
	if ((uii >= 2) && (inbuf_ptr[uii - 2] == '\r')) {
	  inbuf_ptr[uii - 2] = '\n';
	  inbuf_ptr = &(inbuf_ptr[uii - 1]);
	}
      } else {
	inbuf_ptr[uii] = '\n';
	inbuf_ptr = &(inbuf_ptr[uii + 1]);
      }
    } else {
      *inbuf_ptr++ = '\n';
    }
  }
  if (fclose_null(&infile)) {
    goto rserve_call_ret_READ_FAIL;
  }
  if (inbuf_ptr == inbuf_start) {
    logprint("Error: Empty --R file.\n");
    goto rserve_call_ret_INVALID_FORMAT;
  }
  wkspace_alloc((uintptr_t)(inbuf_ptr - inbuf_start));
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
  LOGPRINTFWW5("--R%s debug: writing to %s ... ", rplugin_debug? " debug" : "", outname);
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
