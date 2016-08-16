#include "plink_common.h"

#include <ctype.h>
#include <time.h>
// no more mmap() dependency
// #include <fcntl.h>

#ifndef _WIN32
// #include <sys/mman.h>
#include <unistd.h>
#endif

// #include <sys/stat.h>
#include <sys/types.h>

#define D_EPSILON 0.000244140625

int32_t bgen_to_gen(char* bgenname, char* out_genname, uint32_t snpid_chr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* infile = nullptr;
  FILE* out_genfile = nullptr;
  uintptr_t mc_ct = 0;
  uintptr_t max_mc_len = 0;
  uintptr_t line_idx = 0;
  char* loadbuf = nullptr;
  char* sorted_mc = nullptr;

  uint32_t sample_ct = 0;
  uint32_t col_ct = 3;
  uint32_t is_binary_pheno = 0;
  uint32_t bgen_hardthresh = 0;
  uint32_t marker_ct = 0;
  int32_t retval = 0;
  uint32_t uint_arr[5];
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* wptr;
  uintptr_t* writebuf;
  uintptr_t* ulptr;
  uint16_t* bgen_probs;
  uint16_t* usptr;
  uintptr_t loadbuf_size;
  uintptr_t slen;
  uintptr_t cur_word;
  uintptr_t ulii;
  uintptr_t uljj;
  double dxx;
  double dyy;
  double dzz;
  double drand;
  uLongf zlib_ulongf;
  uint32_t missing_pheno_len;
  uint32_t marker_uidx;
  uint32_t sample_ct4;
  uint32_t sample_ctl2;
  uint32_t sample_idx;
  uint32_t col_idx;
  uint32_t shiftval;
  uint32_t bgen_compressed;
  uint32_t bgen_multichar_alleles;
  uint32_t identical_alleles;
  uint32_t ujj;
  uint32_t ukk;
  int32_t ii;
  uint16_t usii;
  uint16_t usjj;
  uint16_t uskk;
  char cc;
  char cc2;
  {
    if (fopen_checked(out_genname, FOPEN_WB, &out_genfile)) {
      goto bgen_to_gen_ret_OPEN_FAIL;
    }

    uint32_t uint_arr[5];
    if (fopen_checked(bgenname, FOPEN_RB, &in_bgenfile)) {
      goto bgen_to_gen_ret_OPEN_FAIL;
    }
    if (fread(uint_arr, 1, 20, in_bgenfile) < 20) {
      goto bgen_to_gen_ret_READ_FAIL;
    }
    if (uint_arr[1] > uint_arr[0]) {
      logerrprint("Error: Invalid .bgen header.\n");
      goto bgen_to_gen_ret_INVALID_FORMAT;
    }
    const uint32_t raw_marker_ct = uint_arr[2];
    if (!raw_marker_ct) {
      logerrprint("Error: .bgen file contains no variants.\n");
      goto bgen_to_gen_ret_INVALID_FORMAT;
    }
    const uint32_t sample_ct = uint_arr[3];
    if (uint_arr[4] && (uint_arr[4] != 0x6e656762)) {
      logerrprint("Error: Invalid .bgen magic number.\n");
      goto bgen_to_gen_ret_INVALID_FORMAT;
    }
    if (fseeko(in_bgenfile, uint_arr[1], SEEK_SET)) {
      goto bgen_to_gen_ret_READ_FAIL;
    }
    uint32_t uii;
    if (fread(&uii, 1, 4, in_bgenfile) < 4) {
      goto bgen_to_gen_ret_READ_FAIL;
    }
    if (uii & (~5)) {
      uii = (uii >> 2) & 15;
      if (uii == 2) {
	logerrprint("Error: BGEN v1.2 input is not yet supported.  Use gen-convert or a similar tool\nto downcode to BGEN v1.1.\n");
      } else if (uii > 2) {
	logerrprint("Error: Unrecognized BGEN version.  Use gen-convert or a similar tool to\ndowncode to BGEN v1.1.\n");
      } else {
	logerrprint("Error: Unrecognized flags in .bgen header.\n");
      }
      goto bgen_to_gen_ret_INVALID_FORMAT;
    }

    // supports BGEN v1.0 and v1.1.
    uint16_t* bgen_probs = (uint16_t*)bigstack_alloc(6LU * sample_ct);
    if (!bgen_probs) {
      goto bgen_to_gen_ret_NOMEM;
    }
    char* loadbuf = (char*)g_bigstack_base;
    uintptr_t loadbuf_size = bigstack_left();
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size < 3 * 65536) {
      goto bgen_to_gen_ret_NOMEM;
    }
    if (fseeko(in_bgenfile, 4 + uint_arr[0], SEEK_SET)) {
      goto bgen_to_gen_ret_READ_FAIL;
    }
    bgen_compressed = uii & 1;
    bgen_multichar_alleles = (uii >> 2) & 1;
    memcpyl3(g_textbuf, " 0 ");
    for (marker_uidx = 0; marker_uidx < raw_marker_ct; marker_uidx++) {
      if (fread(&uii, 1, 4, in_bgenfile) < 4) {
	goto bgen_to_gen_ret_READ_FAIL;
      }
      if (uii != sample_ct) {
	logerrprint("Error: Unexpected number of samples specified in SNP block header.\n");
	goto bgen_to_gen_ret_INVALID_FORMAT;
      }
      if (bgen_multichar_alleles) {
	// v1.1
	if (fread(&usii, 1, 2, in_bgenfile) < 2) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (!snpid_chr) {
	  if (fseeko(in_bgenfile, usii, SEEK_CUR)) {
	    goto bgen_to_gen_ret_READ_FAIL;
	  }
	  bufptr = loadbuf;
	} else {
	  if (!usii) {
	    logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
	    goto bgen_to_gen_ret_INVALID_FORMAT;
	  }
	  if (fread(loadbuf, 1, usii, in_bgenfile) < usii) {
	    goto bgen_to_gen_ret_READ_FAIL;
	  }
	  loadbuf[usii] = '\0';
	  bufptr = &(loadbuf[usii + 1]);
	}
	if (fread(&usjj, 1, 2, in_bgenfile) < 2) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (!usjj) {
	  logerrprint("Error: Length-0 rsID in .bgen file.\n");
	  goto bgen_to_gen_ret_INVALID_FORMAT;
	}
	if (fread(bufptr, 1, usjj, in_bgenfile) < usjj) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	bufptr2 = &(bufptr[usjj]);
	if (fread(&uskk, 1, 2, in_bgenfile) < 2) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (!snpid_chr) {
	  if (!uskk) {
	    logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
	    goto bgen_to_gen_ret_INVALID_FORMAT;
	  }
	  usii = uskk;
	  if (fread(bufptr2, 1, usii, in_bgenfile) < usii) {
	    goto bgen_to_gen_ret_READ_FAIL;
	  }
	  if ((usii == 2) && (!memcmp(bufptr2, "NA", 2))) {
	    // convert 'NA' to 0
	    usii = 1;
	    memcpy(bufptr2, "0", 2);
	  } else {
	    bufptr2[usii] = '\0';
	  }
	} else {
	  if (fseeko(in_bgenfile, uskk, SEEK_CUR)) {
	    goto bgen_to_gen_ret_READ_FAIL;
	  }
	  bufptr2 = loadbuf;
	}
	if (fread(uint_arr, 1, 8, in_bgenfile) < 8) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (!uint_arr[1]) {
	  logerrprint("Error: Length-0 allele ID in .bgen file.\n");
	  goto bgen_to_gen_ret_INVALID_FORMAT;
	}
	// bufptr2 (chromosome code) is already zero-terminated, with known
	// length usii
	int32_t cur_chrom_code;
	retval = get_or_add_chrom_code(bufptr2, ".bgen file", 0, usii, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
	if (retval) {
	  goto bgen_to_gen_ret_1;
	}
	if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	  // skip rest of current SNP
	  if (fseeko(in_bgenfile, uint_arr[1], SEEK_CUR)) {
	    goto bgen_to_gen_ret_READ_FAIL;
	  }
	  if (fread(&uii, 1, 4, in_bgenfile) < 4) {
	    goto bgen_to_gen_ret_READ_FAIL;
	  }
	  if (bgen_compressed) {
	    if (fseeko(in_bgenfile, uii, SEEK_CUR)) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	    if (fread(&uii, 1, 4, in_bgenfile) < 4) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	    if (fseeko(in_bgenfile, uii, SEEK_CUR)) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	  } else {
	    if (fseeko(in_bgenfile, uii + ((uint64_t)sample_ct) * 6, SEEK_CUR)) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	  }
	  continue;
	}
	fputs(bufptr2, outfile_bim);
	if (putc_checked(' ', outfile_bim)) {
	  goto bgen_to_gen_ret_WRITE_FAIL;
	}
	fwrite(bufptr, 1, usjj, outfile_bim);
	bufptr = uint32toa_x(uint_arr[0], ' ', &(g_textbuf[3]));
	fwrite(g_textbuf, 1, bufptr - g_textbuf, outfile_bim);

	// halve the limit since there are two alleles
	// (may want to enforce NON_BIGSTACK_MIN allele length limit?)
	if (uint_arr[1] >= loadbuf_size / 2) {
	  if (loadbuf_size < MAXLINEBUFLEN) {
	    goto bgen_to_gen_ret_NOMEM;
	  }
	  logerrprint("Error: Excessively long allele in .bgen file.\n");
	  goto bgen_to_gen_ret_INVALID_FORMAT;
	}
	if (fread(loadbuf, 1, uint_arr[1], in_bgenfile) < uint_arr[1]) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	loadbuf[uint_arr[1]] = ' ';
	if (fread(&uii, 1, 4, in_bgenfile) < 4) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (uii >= loadbuf_size / 2) {
	  if (loadbuf_size < MAXLINEBUFLEN) {
	    goto bgen_to_gen_ret_NOMEM;
	  }
	  logerrprint("Error: Excessively long allele in .bgen file.\n");
	  goto bgen_to_gen_ret_INVALID_FORMAT;
	}
	bufptr = &(loadbuf[uint_arr[1] + 1]);
	if (fread(bufptr, 1, uii, in_bgenfile) < uii) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	bufptr[uii] = '\n';
	identical_alleles = (uii == uint_arr[1]) && (!memcmp(loadbuf, bufptr, uii));
	if (!identical_alleles) {
	  if (fwrite_checked(loadbuf, uint_arr[1] + uii + 2, outfile_bim)) {
	    goto bgen_to_gen_ret_WRITE_FAIL;
	  }
	} else {
	  fputs("0 ", outfile_bim);
	  if (fwrite_checked(bufptr, uii + 1, outfile_bim)) {
	    goto bgen_to_gen_ret_WRITE_FAIL;
	  }
	}
      } else {
	// v1.0
	uii = 0;
	if (fread(&uii, 1, 1, in_bgenfile) < 1) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (fread(loadbuf, 1, 2 * uii + 9, in_bgenfile) < (2 * uii + 9)) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	// save marker ID length since we might clobber it
	ukk = (unsigned char)(loadbuf[uii + 1]);
	int32_t cur_chrom_code;
	if (!snpid_chr) {
	  cur_chrom_code = ((unsigned char)(loadbuf[2 * uii + 2]));
	  if (cur_chrom_code > 24) {
	    if (cur_chrom_code == 255) {
	      // unknown
	      cur_chrom_code = 0;
	    } else if (cur_chrom_code > 252) {
	      // XY or MT
	      cur_chrom_code = cur_chrom_code - 228;
	    } else {
	      logerrprint("Error: Invalid chromosome code in BGEN v1.0 file.\n");
	      goto bgen_to_gen_ret_INVALID_FORMAT;
	    }
	  }
	  uint32toa_x((uint32_t)cur_chrom_code, '\0', loadbuf);
	  bufptr = loadbuf;
	} else {
	  ujj = (unsigned char)loadbuf[0];
	  bufptr = &(loadbuf[1]);
	  if ((ujj == 2) && (!memcmp(bufptr, "NA", 2))) {
	    *bufptr = '0';
	    ujj = 1;
	  }
	  bufptr[ujj] = '\0';
	  retval = get_or_add_chrom_code(bufptr, ".bgen file", 0, ujj, allow_extra_chroms, chrom_info_ptr, &cur_chrom_code);
	  if (retval) {
	    goto bgen_to_gen_ret_1;
	  }
	}
	if (!is_set(chrom_info_ptr->chrom_mask, cur_chrom_code)) {
	  if (bgen_compressed) {
	    if (fread(&uii, 1, 4, in_bgenfile) < 4) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	    if (fseeko(in_bgenfile, uii, SEEK_CUR)) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	  } else {
	    if (fseeko(in_bgenfile, ((uint64_t)sample_ct) * 6, SEEK_CUR)) {
	      goto bgen_to_gen_ret_READ_FAIL;
	    }
	  }
	  continue;
	}
	fputs(bufptr, outfile_bim);
	if (putc_checked(' ', outfile_bim)) {
	  goto bgen_to_gen_ret_WRITE_FAIL;
	}
	fwrite(&(loadbuf[uii + 2]), 1, ukk, outfile_bim);
	memcpy(&ujj, &(loadbuf[2 * uii + 3]), 4);
	bufptr = uint32toa_x(ujj, ' ', &(g_textbuf[3]));
	identical_alleles = (loadbuf[2 * uii + 7] == loadbuf[2 * uii + 8]);
	if (!identical_alleles) {
	  *bufptr++ = loadbuf[2 * uii + 7];
	} else {
	  *bufptr++ = '0';
	}
	*bufptr++ = ' ';
	*bufptr++ = loadbuf[2 * uii + 8];
	*bufptr++ = '\n';
	if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile_bim)) {
	  goto bgen_to_gen_ret_WRITE_FAIL;
	}
      }
      if (bgen_compressed) {
	if (fread(&uii, 1, 4, in_bgenfile) < 4) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	if (uii > loadbuf_size) {
	  if (loadbuf_size < MAXLINEBUFLEN) {
	    goto bgen_to_gen_ret_NOMEM;
	  }
	  logerrprint("Error: Excessively long compressed SNP block in .bgen file.\n");
	  goto bgen_to_gen_ret_INVALID_FORMAT;
	}
	if (fread(loadbuf, 1, uii, in_bgenfile) < uii) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
	zlib_ulongf = 6 * sample_ct;
	if (uncompress((Bytef*)bgen_probs, &zlib_ulongf, (Bytef*)loadbuf, uii) != Z_OK) {
	  logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
	  goto bgen_to_gen_ret_INVALID_FORMAT;
	}
      } else {
	if (fread(bgen_probs, 1, 6 * sample_ct, in_bgenfile) < 6 * sample_ct) {
	  goto bgen_to_gen_ret_READ_FAIL;
	}
      }
      cur_word = 0;
      shiftval = 0;
      ulptr = writebuf;
      usptr = bgen_probs;
      if (!is_randomized) {
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, usptr = &(usptr[3])) {
	  if (usptr[2] >= bgen_hardthresh) {
	    ulii = 3;
	  } else if (usptr[1] >= bgen_hardthresh) {
	    ulii = 2;
	  } else if (usptr[0] >= bgen_hardthresh) {
	    ulii = 0;
	  } else {
	    ulii = 1;
	  }
	  cur_word |= ulii << shiftval;
	  shiftval += 2;
	  if (shiftval == BITCT) {
	    *ulptr++ = cur_word;
	    cur_word = 0;
	    shiftval = 0;
	  }
	}
      } else {
	uii = 0;
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++, usptr = &(usptr[3])) {
	  // fast handling of common cases
	  ukk = usptr[2];
	  if (ukk >= 32768) {
	    ulii = 3;
	  } else if (usptr[1] >= 32768) {
	    ulii = 2;
	  } else if (usptr[0] >= 32768) {
	    ulii = 0;
	  } else {
	    while (1) {
	      uii >>= 16;
	      if (!uii) {
		uii = sfmt_genrand_uint32(&g_sfmt) | 0x80000000U;
	      }
	      ujj = uii & 32767;
	      if (ujj < ukk) {
		ulii = 3;
		break;
	      } else {
		ukk += usptr[1];
		if (ujj < ukk) {
		  ulii = 2;
		  break;
		} else {
		  ukk += usptr[0];
		  if (ujj < ukk) {
		    ulii = 0;
		    break;
		  } else if (ukk < 32766) {
		    ulii = 1;
		    break;
		  } else {
		    ukk = usptr[2];
		  }
		}
	      }
	    }
	  }
	  cur_word |= ulii << shiftval;
	  shiftval += 2;
	  if (shiftval == BITCT) {
	    *ulptr++ = cur_word;
	    cur_word = 0;
	    shiftval = 0;
	  }
	}
      }
      if (shiftval) {
	*ulptr++ = cur_word;
      }
      if (identical_alleles) {
	for (ulptr = writebuf; ulptr < (&(writebuf[sample_ctl2])); ulptr++) {
	  ulii = *ulptr;
	  *ulptr = ((~ulii) << 1) | ulii | FIVEMASK;
	}
	if (sample_ct % 4) {
	  writebuf[sample_ctl2 - 1] &= (ONELU << (2 * (sample_ct % BITCT2))) - ONELU;
	}
      }
      if (fwrite_checked(writebuf, sample_ct4, outfile)) {
	goto bgen_to_gen_ret_WRITE_FAIL;
      }
      marker_ct++;
      if (!(marker_ct % 1000)) {
	if (marker_ct == marker_uidx + 1) {
	  printf("\r--bgen: %uk variants converted.", marker_ct / 1000);
	} else {
	  printf("\r--bgen: %uk variants converted (out of %u).", marker_ct / 1000, marker_uidx + 1);
	}
	fflush(stdout);
      }
    }
    if (fclose_null(&in_bgenfile)) {
      goto bgen_to_gen_ret_READ_FAIL;
    }
    if (fclose_null(&out_genfile)) {
      goto bgen_to_gen_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    LOGPRINTFWW("%s.bed + %s.bim + %s.fam written.\n", outname, outname, outname);
  }
  while (0) {
  bgen_to_gen_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  bgen_to_gen_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  bgen_to_gen_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  bgen_to_gen_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  bgen_to_gen_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file has fewer tokens than expected.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  bgen_to_gen_ret_SAMPLE_LONG_LINE:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file is pathologically long.\n", line_idx);
    retval = RET_INVALID_FORMAT;
    break;
  bgen_to_gen_ret_INVALID_SAMPLE_HEADER_2:
    logerrprint("Error: Invalid second header line in .sample file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  bgen_to_gen_ret_INVALID_SAMPLE_HEADER_1:
    logerrprint("Error: Invalid first header line in .sample file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  bgen_to_gen_ret_INVALID_FORMAT_2:
    logerrprintb();
  bgen_to_gen_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  bgen_to_gen_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 bgen_to_gen_ret_1:
  fclose_cond(in_bgenfile);
  fclose_cond(out_genfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

int32_t main(int32_t argc, char** argv) {
  if (argc != 3) {
    fputs("Usage: bgen_to_gen [input .bgen] [output .gen]", stdout);
  }
}
