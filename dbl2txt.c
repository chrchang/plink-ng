#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RET_NOMEM 1
#define RET_OPEN_FAIL 2
#define RET_READ_FAIL 3
#define RET_WRITE_FAIL 4
#define RET_INVALID_CMDLINE 5
#define RET_INVALID_FORMAT 6

static const char digit2_table[] = {
  "0001020304050607080910111213141516171819"
  "2021222324252627282930313233343536373839"
  "4041424344454647484950515253545556575859"
  "6061626364656667686970717273747576777879"
  "8081828384858687888990919293949596979899"};

static inline char* memcpya(char* target, const void* source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(target[ct]);
}

static inline char* uint32_write1p5(char* start, uint32_t quotient, uint32_t remainder) {
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 1000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  remainder -= 1000 * quotient;
  if (remainder) {
    quotient = remainder / 10;
    start += 2;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 10 * quotient;
    if (remainder) {
      start[2] = '0' + remainder;
      return &(start[3]);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline void uint32_write4(char* start, uint32_t uii) {
  // Write exactly four digits (padding with zeroes if necessary); useful for
  // e.g. floating point encoders.
  uint32_t quotient = uii / 100;
  uii -= 100 * quotient;
  memcpy(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
}

static inline void uint32_write6(char* start, uint32_t uii) {
  uint32_t quotient = uii / 10000;
  uint32_write4(memcpya(start, &(digit2_table[quotient * 2]), 2), uii - 10000 * quotient);
}

static inline char* uint32_write6trunc(char* start, uint32_t uii) {
  uint32_t quotient = uii / 10000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  uii -= 10000 * quotient;
  if (uii) {
    quotient = uii / 100;
    start += 2;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
      memcpy(start, &(digit2_table[uii * 2]), 2);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

char* double_write6(char* start, double dxx) {
  // 6 sig fig number, 0.999995 <= dxx < 999999.5
  // Just hardcoding all six cases, in the absence of a better approach...
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.99995) {
    if (dxx < 9.999995) {
      dxx += 0.000005;
      quotient = (int32_t)dxx;
      return uint32_write1p5(start, quotient, ((int32_t)(dxx * 100000)) - (quotient * 100000));
    }
    dxx += 0.00005;
    uii = (int32_t)dxx;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    uii = ((int32_t)(dxx * 10000)) - (uii * 10000);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    quotient = uii / 100;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
    double_write6_pretail:
      memcpy(start, &(digit2_table[uii * 2]), 2);
    }
  double_write6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 9999.995) {
    if (dxx < 999.9995) {
      dxx += 0.0005;
      uii = (int32_t)dxx;
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      uii = ((int32_t)(dxx * 1000)) - (uii * 1000);
      if (!uii) {
	return start;
      }
      *start++ = '.';
      quotient = uii / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      uii -= quotient * 10;
      if (!uii) {
        goto double_write6_tail;
      }
      start[2] = '0' + uii;
      return &(start[3]);
    }
    dxx += 0.005;
    uii = (int32_t)dxx;
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uii = ((int32_t)(dxx * 100)) - (uii * 100);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    goto double_write6_pretail;
  } else if (dxx < 99999.95) {
    dxx += 0.05;
    uii = (int32_t)dxx;
    quotient = uii / 10000;
    *start = '0' + quotient;
    remainder = uii - 10000 * quotient;
    quotient = remainder / 100;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    remainder = remainder - 100 * quotient;
    start = memcpya(start, &(digit2_table[remainder * 2]), 2);
    uii = ((int32_t)(dxx * 10)) - (uii * 10);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    *start = '0' + uii;
    return &(start[1]);
  } else {
    uint32_write6(start, (int32_t)(dxx + 0.5));
    return &(start[6]);
  }
}

char* double_g_write(char* start, double dxx) {
  uint32_t xp10 = 0;
  uint32_t uii;
  if (dxx != dxx) {
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.999995e-5) {
    // 6 sig fig exponential notation, small
    if (dxx < 9.999995e-16) {
      if (dxx < 9.999995e-128) {
	if (dxx == 0.0) {
	  *start = '0';
	  return &(start[1]);
	} else if (dxx < 9.999995e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.999995e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.999995e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.999995e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.999995e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.999995e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.999995e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.999995e-1) {
      dxx *= 10;
      xp10++;
    }
    dxx += 0.000005;
    uii = (int32_t)dxx;
    start = memcpya(uint32_write1p5(start, uii, ((int32_t)(dxx * 100000)) - (uii * 100000)), "e-", 2);
    if (xp10 >= 100) {
      uii = xp10 / 100;
      *start++ = '0' + uii;
      xp10 -= 100 * uii;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 999999.5) {
    // 6 sig fig exponential notation, large
    if (dxx >= 9.999995e15) {
      if (dxx >= 9.999995e127) {
	if (dxx == INFINITY) {
	  *((uint32_t*)start) = *((uint32_t*)"inf");
	  return &(start[3]);
	} else if (dxx >= 9.999995e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.999995e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.999995e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.999995e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.999995e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.999995e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.999995e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.999995e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    dxx += 0.000005;
    uii = (int32_t)dxx;
    start = memcpya(uint32_write1p5(start, uii, ((int32_t)(dxx * 100000)) - (uii * 100000)), "e+", 2);
    if (xp10 >= 100) {
      uii = xp10 / 100;
      *start++ = '0' + uii;
      xp10 -= 100 * uii;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 0.9999995) {
    return double_write6(start, dxx);
  } else {
    // 6 sig fig decimal, no less than ~0.0001
    start = memcpya(start, "0.", 2);
    if (dxx < 9.999995e-3) {
      dxx *= 100;
      start = memcpya(start, "00", 2);
    }
    if (dxx < 9.999995e-2) {
      dxx *= 10;
      *start++ = '0';
    }
    return uint32_write6trunc(start, (int32_t)((dxx * 1000000) + 0.5));
  }
}

char* float_write6(char* start, float dxx) {
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.99995) {
    if (dxx < 9.999995) {
      dxx += 0.000005;
      quotient = (int32_t)dxx;
      return uint32_write1p5(start, quotient, ((int32_t)(dxx * 100000)) - (quotient * 100000));
    }
    dxx += 0.00005;
    uii = (int32_t)dxx;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    uii = ((int32_t)(dxx * 10000)) - (uii * 10000);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    quotient = uii / 100;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
    float_write6_pretail:
      memcpy(start, &(digit2_table[uii * 2]), 2);
    }
  float_write6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 9999.995) {
    if (dxx < 999.9995) {
      dxx += 0.0005;
      uii = (int32_t)dxx;
      quotient = uii / 100;
      *start = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
      uii = ((int32_t)(dxx * 1000)) - (uii * 1000);
      if (!uii) {
	return start;
      }
      *start++ = '.';
      quotient = uii / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      uii -= quotient * 10;
      if (!uii) {
        goto float_write6_tail;
      }
      start[2] = '0' + uii;
      return &(start[3]);
    }
    dxx += 0.005;
    uii = (int32_t)dxx;
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uii = ((int32_t)(dxx * 100)) - (uii * 100);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    goto float_write6_pretail;
  } else if (dxx < 99999.95) {
    dxx += 0.05;
    uii = (int32_t)dxx;
    quotient = uii / 10000;
    *start = '0' + quotient;
    remainder = uii - 10000 * quotient;
    quotient = remainder / 100;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    remainder = remainder - 100 * quotient;
    start = memcpya(start, &(digit2_table[remainder * 2]), 2);
    uii = ((int32_t)(dxx * 10)) - (uii * 10);
    if (!uii) {
      return start;
    }
    *start = '.';
    start[1] = '0' + uii;
    return &(start[2]);
  } else {
    uint32_write6(start, (int32_t)(dxx + 0.5));
    return &(start[6]);
  }
}

char* float_g_write(char* start, float dxx) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  if (dxx != dxx) {
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.999995e-5) {
    if (dxx < 9.999995e-16) {
      if (dxx == 0.0) {
	*start = '0';
	return &(start[1]);
      } else if (dxx < 9.999995e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      } else {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.999995e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.999995e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.999995e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.999995e-1) {
      dxx *= 10;
      xp10++;
    }
    dxx += 0.000005;
    quotient = (int32_t)dxx;
    return memcpya(memcpya(uint32_write1p5(start, quotient, ((int32_t)(dxx * 100000)) - (quotient * 100000)), "e-", 2), &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 999999.5) {
    if (dxx >= 9.999995e15) {
      if (dxx == INFINITY) {
	*((uint32_t*)start) = *((uint32_t*)"inf");
	return &(start[3]);
      } else if (dxx >= 9.999995e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      } else {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.999995e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.999995e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.999995e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.999995e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    dxx += 0.000005;
    quotient = (int32_t)dxx;
    return memcpya(memcpya(uint32_write1p5(start, quotient, ((int32_t)(dxx * 100000)) - (quotient * 100000)), "e+", 2), &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 0.9999995) {
    return float_write6(start, dxx);
  } else {
    // 6 sig fig decimal, no less than ~0.0001
    start = memcpya(start, "0.", 2);
    if (dxx < 9.999995e-3) {
      dxx *= 100;
      start = memcpya(start, "00", 2);
    }
    if (dxx < 9.999995e-2) {
      dxx *= 10;
      *start++ = '0';
    }
    return uint32_write6trunc(start, (int32_t)((dxx * 1000000) + 0.5));
  }
}

void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

int32_t main(int32_t argc, char** argv) {
  FILE* infile = NULL;
  FILE* outfile = NULL;
  uint32_t is_triangle = 0;
  uint32_t is_float = 0;
  int32_t retval = 0;
  double* readbuf = NULL;
  float* readbuf_f = NULL;
  char* writebuf = NULL;
  double* dptr;
  float* fptr;
  char* wptr;
  uint64_t file_size;
  uintptr_t row_idx;
  uintptr_t col_idx;
  uintptr_t row_ct;
  uintptr_t col_ct;
  if ((argc != 3) && (argc != 4)) {
    goto main_ret_INVALID_CMDLINE;
  }
  if (argc == 4) {
    if (strcmp("--float", argv[1]) && strcmp("-f", argv[1])) {
      goto main_ret_INVALID_CMDLINE;
    }
    is_float = 1;
  }
  infile = fopen(argv[argc - 2], "rb");
  if (!infile) {
    printf("Failed to open %s.\n", argv[argc - 2]);
    goto main_ret_OPEN_FAIL;
  }
  if (fseeko(infile, 0, SEEK_END)) {
    goto main_ret_READ_FAIL;
  }
  file_size = ftello(infile);
  if (!file_size) {
    printf("Error: %s is empty.\n", argv[argc - 2]);
    goto main_ret_INVALID_FORMAT;
  }
  row_ct = (int32_t)(sqrt((double)((int64_t)(file_size / sizeof(double)))) + 0.0000000001);
  if (is_float || (((uint64_t)row_ct) * row_ct * sizeof(double) != file_size)) {
    row_ct = (int32_t)sqrt((double)((int64_t)(file_size / (sizeof(double) / 2))));
    if (is_float || (((uint64_t)row_ct) * (row_ct + 1) * (sizeof(double) / 2) != file_size)) {
      is_float = 1;
      row_ct = (int32_t)(sqrt((double)((int64_t)(file_size / sizeof(float)))) + 0.0000000001);
      if ((((uint64_t)row_ct) * row_ct * sizeof(float) != file_size)) {
	row_ct = (int32_t)sqrt((double)((int64_t)(file_size / (sizeof(float) / 2))));
	if ((((uint64_t)row_ct) * (row_ct + 1) * (sizeof(float) / 2) != file_size)) {
	  printf("Error: %s's size indicates it is not a binary floating-point matrix.\n", argv[argc - 2]);
	  goto main_ret_INVALID_FORMAT;
	}
	is_triangle = 1;
      }
    } else {
      is_triangle = 1;
    }
  }
  col_ct = row_ct;
  rewind(infile);
  outfile = fopen(argv[argc - 1], "w");
  if (!outfile) {
    printf("Failed to open %s.\n", argv[argc - 1]);
    goto main_ret_OPEN_FAIL;
  }
  if (!is_float) {
    readbuf = (double*)malloc(row_ct * sizeof(double));
    if (!readbuf) {
      goto main_ret_NOMEM;
    }
  } else {
    readbuf_f = (float*)malloc(row_ct * sizeof(float));
    if (!readbuf_f) {
      goto main_ret_NOMEM;
    }
  }
  writebuf = (char*)malloc(row_ct * 16);
  if (!writebuf) {
    goto main_ret_NOMEM;
  }
  for (row_idx = 0; row_idx < row_ct; row_idx++) {
    if (is_triangle) {
      col_ct = row_idx + 1;
    }
    if (!is_float) {
      if (fread(readbuf, 1, col_ct * sizeof(double), infile) < col_ct * sizeof(double)) {
	goto main_ret_READ_FAIL;
      }
      dptr = readbuf;
      wptr = writebuf;
      for (col_idx = 0; col_idx < col_ct; col_idx++) {
	wptr = double_g_write(wptr, *dptr++);
	*wptr++ = '\t';
      }
    } else {
      if (fread(readbuf_f, 1, col_ct * sizeof(float), infile) < col_ct * sizeof(float)) {
	goto main_ret_READ_FAIL;
      }
      fptr = readbuf_f;
      wptr = writebuf;
      for (col_idx = 0; col_idx < col_ct; col_idx++) {
	wptr = float_g_write(wptr, *fptr++);
	*wptr++ = '\t';
      }
    }
    wptr[-1] = '\n';
    fwrite(writebuf, 1, wptr - writebuf, outfile);
    if (ferror(outfile)) {
      goto main_ret_WRITE_FAIL;
    }
  }
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
  main_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  main_ret_INVALID_CMDLINE:
    fputs("Usage: dbl2text <--float> [input filename] [output filename]\n", stdout);
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  fclose_cond(outfile);
  if (readbuf) {
    free(readbuf);
  }
  if (readbuf_f) {
    free(readbuf_f);
  }
  if (writebuf) {
    free(writebuf);
  }
  return retval;
}
