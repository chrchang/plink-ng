// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Rconnection.cc is part of a continuously updated package and uses
// C++-specific language constructs, so we omit this module when compiling with
// gcc instead of g++.
#if defined __cplusplus && !defined _WIN32

#define MAIN
#define SOCK_ERRORS
#include "sisocks.h"
#include "Rconnection.h"

#include "plink_common.h"
#include "plink_cluster.h"

#define RPLUGIN_BLOCK_SIZE 100

int32_t rserve_call(char* rplugin_fname, char* rplugin_host_or_socket, int32_t rplugin_port, uint32_t rplugin_debug, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uintptr_t* pheno_c, double* pheno_d, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uintptr_t covar_ct, double* covar_d, char* outname, char* outname_end) {
  // See PLINK 1.07 r.cpp.
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  int32_t* geno_int_buf = nullptr;
  Rinteger* r_n = nullptr;
  Rinteger* r_s = nullptr;
  Rdouble* r_p = nullptr;
  Rdouble* r_cov = nullptr;
  Rconnection* rc = nullptr;
  char* chrom_name_ptr = nullptr;
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
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_SLEN];
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
  Rexp* r_eval_dummy;  // bugfix (11 Sep 2018): must delete rc->eval return-val
  double* pheno_d2;
  double* covar_d2 = nullptr;
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
  if (!pheno_nm_ct) {
    // bugfix (11 Sep 2018): this was segfaulting instead of printing an error
    // message
    logerrprint("Error: --R only processes samples with nonmissing phenotype values, and all\nphenotype values are missing.\n");
    goto rserve_call_ret_1;
  }
  if (fopen_checked(rplugin_fname, "r", &infile)) {
    goto rserve_call_ret_OPEN_FAIL;
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(pheno_nm_ctl2 * RPLUGIN_BLOCK_SIZE, &loadbuf)) {
    goto rserve_call_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  for (ulii = 1; ulii <= RPLUGIN_BLOCK_SIZE; ulii++) {
    loadbuf[ulii * pheno_nm_ctl2 - 1] = 0;
  }
  inbuf_start = (char*)g_bigstack_base;
  inbuf_end = inbuf_start;
  while (1) {
    if (((uintptr_t)g_bigstack_end) - ((uintptr_t)inbuf_start) < MAXLINELEN) {
      goto rserve_call_ret_NOMEM;
    }
    inbuf_end[MAXLINELEN - 1] = ' ';
    if (!fgets(inbuf_end, MAXLINELEN, infile)) {
      break;
    }
    line_idx++;
    if (!(inbuf_end[MAXLINELEN - 1])) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --R file is pathologically long.\n", line_idx);
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
    logerrprint("Error: Empty --R file.\n");
    goto rserve_call_ret_INVALID_FORMAT;
  }
  *inbuf_end = '\0';
  bigstack_alloc(1 + ((uintptr_t)(inbuf_end - inbuf_start)));
  if (pheno_c) {
    if (bigstack_alloc_d(pheno_nm_ct, &pheno_d2)) {
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
  if (covar_ct) {
    if (pheno_nm_ct == unfiltered_sample_ct) {
      covar_d2 = covar_d;
    } else {
      if (bigstack_alloc_d(pheno_nm_ct * covar_ct, &covar_d2)) {
	goto rserve_call_ret_NOMEM;
      }
      for (sample_uidx = 0, sample_idx = 0; sample_idx < pheno_nm_ct; ++sample_uidx, ++sample_idx) {
	next_set_unsafe_ck(pheno_nm, &sample_uidx);
	memcpy(&(covar_d2[sample_idx * covar_ct]), &(covar_d[sample_uidx * covar_ct]), covar_ct * sizeof(double));
      }
    }
  }
  if (cluster_ct) {
    if (bigstack_alloc_i(unfiltered_sample_ct, &sample_to_cluster)) {
      goto rserve_call_ret_NOMEM;
    }
    fill_int_one(pheno_nm_ct, sample_to_cluster);
    fill_unfiltered_sample_to_cluster(unfiltered_sample_ct, cluster_ct, cluster_map, cluster_starts, (uint32_t*)sample_to_cluster);
    inplace_collapse_uint32_incl((uint32_t*)sample_to_cluster, unfiltered_sample_ct, pheno_nm, pheno_nm_ct);
    bigstack_shrink_top(sample_to_cluster, pheno_nm_ct * sizeof(int32_t));
  } else {
    if (bigstack_calloc_i(pheno_nm_ct, &sample_to_cluster)) {
      goto rserve_call_ret_NOMEM;
    }
  }
  if (!rplugin_debug) {
    if (bigstack_alloc_i(RPLUGIN_BLOCK_SIZE * ((uintptr_t)pheno_nm_ct), &geno_int_buf)) {
      goto rserve_call_ret_NOMEM;
    }
    if (rplugin_port == -2) {
      rplugin_port = 6311;
    }
    rc = new Rconnection(rplugin_host_or_socket, rplugin_port);
    ii = rc->connect();
    if (ii) {
      sockerrorchecks(g_textbuf, 128, -1);
      LOGERRPRINTFWW("Error: Unable to connect (code %d: %s).\n", ii, g_textbuf);
      goto rserve_call_ret_NETWORK;
    }
    rc->eval("options(echo=F)");

    r_s = new Rinteger(sample_to_cluster, pheno_nm_ct);
    r_p = new Rdouble(pheno_d2, pheno_nm_ct);
    r_n = new Rinteger((int32_t*)(&pheno_nm_ct), 1);
    rc->assign("n", r_n);
    rc->assign("PHENO", r_p);
    rc->assign("CLUSTER", r_s);
    r_eval_dummy = rc->eval("CLUSTER[CLUSTER==-1] <- NA");
    if (r_eval_dummy) {
      delete r_eval_dummy;
    }
    if (covar_ct) {
      r_cov = new Rdouble(covar_d2, pheno_nm_ct * covar_ct);
      rc->assign("c", r_cov);
      r_eval_dummy = rc->eval("COVAR<-matrix(c,nrow=n,byrow=T)");
    } else {
      r_eval_dummy = rc->eval("COVAR<-NA");
    }
    if (r_eval_dummy) {
      delete r_eval_dummy;
    }
    memcpy(outname_end, ".auto.R", 8);
  } else {
    memcpy(outname_end, ".debug.R", 9);
  }
  if (fopen_checked(outname, "w", &outfile)) {
    goto rserve_call_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("--R%s: writing to %s ... ", rplugin_debug? " debug" : "", outname);
  fputs("0%", stdout);
  fflush(stdout);
  loop_end = marker_ct / 100;
  if (rplugin_debug) {
    bufptr = memcpya(g_textbuf, "n <- ", 5);
    bufptr = uint32toa(pheno_nm_ct, bufptr);
    bufptr = memcpya(bufptr, "\nPHENO <- c( ", 13);
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
      goto rserve_call_ret_WRITE_FAIL;
    }
    for (sample_idx = 0; sample_idx < pheno_nm_ct - 1; sample_idx++) {
      bufptr = dtoa_g(pheno_d2[sample_idx], g_textbuf);
      bufptr = memcpya(bufptr, ", ", 2);
      fwrite(g_textbuf, 1, (uintptr_t)(bufptr - g_textbuf), outfile);
    }
    bufptr = dtoa_g(pheno_d2[sample_idx], g_textbuf);
    bufptr = memcpya(bufptr, " ) \n", 4);
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
      goto rserve_call_ret_WRITE_FAIL;
    }
    if (covar_ct) {
      fputs("c <- c( ", outfile);
      uljj = pheno_nm_ct * covar_ct - 1;
      for (ulii = 0; ulii < uljj; ulii++) {
	bufptr = dtoa_g(covar_d2[ulii], g_textbuf);
	bufptr = memcpya(bufptr, ", ", 2);
	fwrite(g_textbuf, 1, (uintptr_t)(bufptr - g_textbuf), outfile);
      }
      bufptr = dtoa_g(covar_d2[ulii], g_textbuf);
      fwrite(g_textbuf, 1, (uintptr_t)(bufptr - g_textbuf), outfile);
      fputs(" ) \nCOVAR <- matrix( c , nrow = n , byrow=T)\n", outfile);
    } else {
      // old code (this might be better?  but --R backward compatibility is
      // more important than --R-debug...):
      // fputs("COVAR <- matrix( NA , nrow = n , ncol = 0 , byrow = T)\n", outfile);
      fputs("COVAR<-NA\n", outfile);
    }
    fputs("CLUSTER <- c( ", outfile);
    for (sample_idx = 0; sample_idx < pheno_nm_ct - 1; sample_idx++) {
      bufptr = int32toa(sample_to_cluster[sample_idx], g_textbuf);
      bufptr = memcpya(bufptr, ", ", 2);
      fwrite(g_textbuf, 1, (uintptr_t)(bufptr - g_textbuf), outfile);
    }
    bufptr = int32toa(sample_to_cluster[sample_idx], g_textbuf);
    bufptr = memcpya(bufptr, " ) \nCLUSTER[CLUSTER==-1] <- NA\n", 31);
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
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
      if (load_and_collapse_incl(unfiltered_sample_ct, pheno_nm_ct, pheno_nm, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, ulptr)) {
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
      r_eval_dummy = rc->eval("GENO<-matrix(g,nrow=n,byrow=T)");
      if (r_eval_dummy) {
        delete r_eval_dummy;
      }
      r_eval_dummy = rc->eval("GENO[GENO==-1] <- NA");
      if (r_eval_dummy) {
        delete r_eval_dummy;
      }
      delete r_l;
      delete r_g;
      r_eval_dummy = rc->eval(inbuf_start);
      if (r_eval_dummy) {
        delete r_eval_dummy;
      }
      r_data = (Rdouble*)rc->eval("Rplink(PHENO,GENO,CLUSTER,COVAR)");
      if (r_data) {
	dptr = r_data->doubleArray();
	for (marker_uidx = marker_uidx_base, block_offset = 0; block_offset < block_size; marker_uidx++, block_offset++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  if (marker_uidx >= chrom_end) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1];
	    } while (marker_uidx >= chrom_end);
	    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, uii, &chrom_name_len, chrom_name_buf);
	  }
	  bufptr = memcpyax(g_textbuf, chrom_name_ptr, chrom_name_len, ' ');
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = uint32toa_w10x(marker_pos[marker_uidx], ' ', bufptr);
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto rserve_call_ret_WRITE_FAIL;
	  }
	  fputs_w4(marker_allele_ptrs[2 * marker_uidx], outfile);
	  g_textbuf[0] = ' ';
	  bufptr = &(g_textbuf[1]);
	  cur_data_len = (int32_t)(*dptr++);
	  for (uii = 0; uii < cur_data_len; uii++) {
	    dxx = *dptr++;
	    if (realnum(dxx)) {
	      bufptr = dtoa_g(dxx, bufptr);
	    } else {
	      bufptr = memcpya(bufptr, "NA", 2);
	    }
	    *bufptr++ = '\t';
	    if (bufptr > &(g_textbuf[MAXLINELEN])) {
	      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
		goto rserve_call_ret_WRITE_FAIL;
	      }
	      bufptr = g_textbuf;
	    }
	  }
	  *bufptr++ = '\n';
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto rserve_call_ret_WRITE_FAIL;
	  }
	}
	delete r_data;
      } else {
	for (marker_uidx = marker_uidx_base, block_offset = 0; block_offset < block_size; marker_uidx++, block_offset++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  if (marker_uidx >= chrom_end) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1];
	    } while (marker_uidx >= chrom_end);
	    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, uii, &chrom_name_len, chrom_name_buf);
	  }
	  bufptr = memcpyax(g_textbuf, chrom_name_ptr, chrom_name_len, ' ');
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = uint32toa_w10x(marker_pos[marker_uidx], ' ', bufptr);
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto rserve_call_ret_WRITE_FAIL;
	  }
	  fputs_w4(marker_allele_ptrs[2 * marker_uidx], outfile);
	  fputs(" NA\n", outfile);
	}
      }
    } else {
      bufptr = memcpya(g_textbuf, "l <- ", 5);
      bufptr = uint32toa(block_size, bufptr);
      bufptr = memcpya(bufptr, "\ng <- c( ", 9);
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto rserve_call_ret_WRITE_FAIL;
      }
      sample_idx = 0;
      while (1) {
        bufptr = g_textbuf;
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
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	  goto rserve_call_ret_WRITE_FAIL;
	}
      }
      bufptr = memcpya(&(bufptr[-2]), " ) \nGENO <- matrix( g , nrow = n ,byrow=T)\nGENO[GENO == -1 ] <- NA \n\n\n", 70);
      if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	goto rserve_call_ret_WRITE_FAIL;
      }
      if (fwrite_checked(inbuf_start, inbuf_end - inbuf_start, outfile)) {
	goto rserve_call_ret_WRITE_FAIL;
      }
      putc_unlocked('\n', outfile);
    }
    marker_idx_base += block_size;
    if (marker_idx_base >= loop_end) {
      if (marker_idx_base < marker_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = (marker_idx_base * 100LLU) / marker_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = (((uint64_t)pct + 1LLU) * marker_ct) / 100;
      }
    }
  }

  if (pct >= 10) {
    putc_unlocked('\b', stdout);
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
    logerrprintb();
  rserve_call_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  rserve_call_ret_NETWORK:
    retval = RET_NETWORK;
    break;
  }
 rserve_call_ret_1:
  fclose_cond(infile);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  delete r_n;
  delete r_p;
  delete r_s;
  delete r_cov;
  return retval;
}

#endif // __cplusplus, !_WIN32
