// dose2plink, memory-efficient C implementation
// Copyright (C) 2014 Christopher Chang; based on a Perl script by Sarah
// Medland (http://genepi.qimr.edu.au/staff/sarahMe/dose2plink ).

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <locale.h>
#include <unistd.h>

#ifdef _WIN32
  #ifndef _WIN64
    #define WINVER 0x0500
  #endif
#else // Unix
  #include <sys/stat.h>
#endif

#ifdef _WIN32
  #define PRId64 "I64d"
  #define PRIu64 "I64u"
  #define fseeko fseeko64
  #define ftello ftello64
  #include <windows.h>
#else
  #ifdef __cplusplus
    #ifndef PRId64
      #define PRId64 "lld"
    #endif
  #endif
#endif

#ifdef __APPLE__
  typedef unsigned long long uint64_t;
  typedef long long int64_t;
#else
  #define uint64_t unsigned long long
  #define int64_t long long
#endif

#ifdef _WIN64
  #define __LP64__
#else
  #ifndef __LP64__
    #ifndef uintptr_t
      #define uintptr_t unsigned long
    #endif
    #ifndef intptr_t
      #define intptr_t long
    #endif
  #endif
#endif

#ifdef __LP64__
  #define ZEROLU 0LLU
  #define ONELU 1LLU
  #ifdef _WIN32
    #ifndef PRIuPTR
      #define PRIuPTR PRIu64
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR PRId64
    #endif
  #else
    #ifndef PRIuPTR
      #define PRIuPTR "lu"
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR "ld"
    #endif
  #endif
#else
  #define ZEROLU 0LU
  #define ONELU 1LU
  #ifndef PRIuPTR
    #define PRIuPTR "lu"
  #endif
  #ifndef PRIdPTR
    #define PRIdPTR "ld"
  #endif
#endif

#include "../zlib-1.2.11/zlib.h"

#ifdef __APPLE__
  #include <sys/sysctl.h>
#endif

#define RET_HELP 1
#define RET_NOMEM 2
#define RET_OPEN_FAIL 3
#define RET_INVALID_CMDLINE 4
#define RET_READ_FAIL 5
#define RET_INVALID_FORMAT 6
#define RET_WRITE_FAIL 7

const char errstr_fopen[] = "Error: Failed to open %s.\n";

unsigned char* wkspace;
unsigned char* wkspace_base;
uintptr_t wkspace_left;

void wkspace_reset(void* new_base) {
  uintptr_t freed_bytes = wkspace_base - (unsigned char*)new_base;
  wkspace_base = (unsigned char*)new_base;
  wkspace_left += freed_bytes;
}

#define WKSPACE_MIN_MB 129
#define WKSPACE_DEFAULT_MB 2048

#define CACHELINE 64

#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - ONELU)))

// defend against very long indels.
#define MINLINEBUFLEN 67108864

#define MAXLINEBUFLEN 0x7fffffc0

void disp_usage(FILE* stream) {
  fputs(
"Usage: dose2plink [flags...]\n\n"
"  -d, --dose [fname]   : Specify full .dose/.mldose filename.  Required.\n"
"  -i, --info [fname]   : Specify full .info/.mlinfo filename.  Required.\n"
"  -o, --out [prefix]   : Set output filename prefix (default 'plink_dosage').\n"
"  -m, --memory [val]   : Set size, in MB, of initial workspace malloc attempt.\n"
"  --no-gz (or '-gz 0') : Turn off zipping of the output .pdat file.\n"
, stream);
}

void dispmsg(int32_t retval) {
  switch (retval) {
  case RET_NOMEM:
    fputs("\nError: Out of memory.\n", stderr);
    break;
  case RET_READ_FAIL:
    fputs("\nError: File read failure.\n", stderr);
    break;
  case RET_WRITE_FAIL:
    fputs("\nError: File write failure.\n", stderr);
    break;
  }
}

void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    fprintf(stderr, errstr_fopen, fname);
    return -1;
  }
  return 0;
}

static inline int32_t putc_checked(int32_t ii, FILE* outfile) {
  putc(ii, outfile);
  return ferror(outfile);
}

static inline int32_t fputs_checked(const char* ss, FILE* outfile) {
  fputs(ss, outfile);
  return ferror(outfile);
}

int32_t fwrite_checked(const void* buf, size_t len, FILE* outfile) {
  while (len > 0x7ffe0000) {
    // OS X can't perform >2GB writes
    fwrite(buf, 1, 0x7ffe0000, outfile);
    buf = &(((unsigned char*)buf)[0x7ffe0000]);
    len -= 0x7ffe0000;
  }
  fwrite(buf, 1, len, outfile);
  return ferror(outfile);
}

void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

int32_t fclose_null(FILE** fptr_ptr) {
  int32_t ii = ferror(*fptr_ptr);
  int32_t jj = fclose(*fptr_ptr);
  *fptr_ptr = NULL;
  return ii || jj;
}

int32_t gzopen_checked(gzFile* target_ptr, const char* fname, const char* mode) {
  *target_ptr = gzopen(fname, mode);
  if (!(*target_ptr)) {
    fprintf(stderr, errstr_fopen, fname);
    return -1;
  }
  return 0;
}

static inline void gzclose_cond(gzFile gz_file) {
  if (gz_file) {
    gzclose(gz_file);
  }
}

unsigned char* wkspace_alloc(uintptr_t size) {
  unsigned char* retval;
  if (wkspace_left < size) {
    return NULL;
  }
  size = CACHEALIGN(size);
  retval = wkspace_base;
  wkspace_base += size;
  wkspace_left -= size;
  return retval;
}

static inline int32_t wkspace_alloc_c_checked(char** dc_ptr, uintptr_t size) {
  *dc_ptr = (char*)wkspace_alloc(size);
  return !(*dc_ptr);
}

static inline int32_t is_eoln_kns(unsigned char ucc) {
  return (ucc < 32);
}

static inline int32_t no_more_tokens_kns(char* sptr) {
  return ((!sptr) || is_eoln_kns(*sptr));
}

static inline char* skip_initial_spaces(char* sptr) {
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

static inline int32_t is_space_or_eoln(unsigned char ucc) {
  return (ucc <= 32);
}

int32_t strcmp_se(char* s_read, const char* s_const, uint32_t len) {
  return memcmp(s_read, s_const, len) || (!is_space_or_eoln(s_read[len]));
}

char* next_token(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while ((*sptr != ' ') && (*sptr != '\t')) {
    if (!(*sptr)) {
      return NULL;
    }
    sptr++;
  }
  return skip_initial_spaces(sptr);
}

char* next_token_mult(char* sptr, uint32_t ct) {
  if (!sptr) {
    return NULL;
  }
  do {
    while ((*sptr != ' ') && (*sptr != '\t')) {
      if (!(*sptr)) {
	return NULL;
      }
      sptr++;
    }
    sptr = skip_initial_spaces(sptr);
  } while (--ct);
  return sptr;
}

// assumes we are currently in a token
static inline char* token_endnn(char* sptr) {
  while (!is_space_or_eoln(*(++sptr)));
  return sptr;
}

static inline char* memcpya(char* target, const void* source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(target[ct]);
}

uint32_t scan_posintptr(char* ss, uintptr_t* valp) {
  // Reads an integer in [1, 2^BITCT - 1].  Assumes first character is
  // nonspace. 
  uintptr_t val = (uint32_t)((unsigned char)*ss) - 48;
  uintptr_t cur_digit;
  if (val < 10) {
    while (1) {
    scan_posintptr_main_loop:
      cur_digit = (uint32_t)((unsigned char)(*(++ss))) - 48;
      if (cur_digit >= 10) {
	if (val) {
	  *valp = val;
	  return 0;
	}
	return 1;
      }
      if ((val >= (~ZEROLU) / 10) && ((val > (~ZEROLU) / 10) || (cur_digit > (~ZEROLU) % 10))) {
	return 1;
      }
      val = val * 10 + cur_digit;
    }
  } else if (val == 0xfffffffbU) {
    val = (uint32_t)((unsigned char)(*(++ss))) - 48;
    if (val < 10) {
      goto scan_posintptr_main_loop;
    }
  }
  return 1;
}

static inline uint32_t scan_double(char* ss, double* valp) {
  char* ss2;
  *valp = strtod(ss, &ss2);
  return (ss == ss2);
}

static const char digit2_table[] = {
  "0001020304050607080910111213141516171819"
  "2021222324252627282930313233343536373839"
  "4041424344454647484950515253545556575859"
  "6061626364656667686970717273747576777879"
  "8081828384858687888990919293949596979899"};

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

void uint32_write4(char* start, uint32_t uii) {
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

static const double banker_round8[] = {0.499999995, 0.500000005};

static inline uint32_t double_bround(double dxx, const double* banker_round) {
  uint32_t result = (int32_t)dxx;
  return result + (int32_t)((dxx - ((int32_t)result)) + banker_round[result & 1]);
}

// These are separate functions so the compiler can optimize the integer
// divisions.
static inline void double_bround1(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10;
  *remainderp = remainder - (*quotientp) * 10; 
}

static inline void double_bround2(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100;
  *remainderp = remainder - (*quotientp) * 100; 
}

static inline void double_bround3(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000;
  *remainderp = remainder - (*quotientp) * 1000; 
}

static inline void double_bround4(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10000;
  *remainderp = remainder - (*quotientp) * 10000; 
}

static inline void double_bround5(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100000;
  *remainderp = remainder - (*quotientp) * 100000; 
}

char* double_write6(char* start, double dxx) {
  // 6 sig fig number, 0.999995 <= dxx < 999999.5
  // Just hardcoding all six cases, in the absence of a better approach...
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.999949999999) {
    if (dxx < 9.9999949999999) {
      double_bround5(dxx, banker_round8, &quotient, &remainder);
      return uint32_write1p5(start, quotient, remainder);
    }
    double_bround4(dxx, banker_round8, &quotient, &remainder);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    double_write6_pretail:
      memcpy(start, &(digit2_table[remainder * 2]), 2);
    }
  double_write6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 9999.9949999999) {
    if (dxx < 999.99949999999) {
      double_bround3(dxx, banker_round8, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto double_write6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    double_bround2(dxx, banker_round8, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto double_write6_pretail;
  } else if (dxx < 99999.949999999) {
    double_bround1(dxx, banker_round8, &uii, &remainder);
    quotient = uii / 10000;
    *start = '0' + quotient;
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    uii = uii - 100 * quotient;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    *start = '0' + remainder;
    return &(start[1]);
  } else {
    uint32_write6(start, double_bround(dxx, banker_round8));
    return &(start[6]);
  }
}

char* double_g_write(char* start, double dxx) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9999949999999e-5) {
    // 6 sig fig exponential notation, small
    if (dxx < 9.9999949999999e-16) {
      if (dxx < 9.9999949999999e-128) {
	if (dxx == 0.0) {
	  *start = '0';
	  return &(start[1]);
	} else if (dxx < 9.9999949999999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9999949999999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9999949999999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9999949999999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9999949999999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9999949999999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999949999999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999949999999e-1) {
      dxx *= 10;
      xp10++;
    }
    double_bround5(dxx, banker_round8, &quotient, &remainder);
    start = memcpya(uint32_write1p5(start, quotient, remainder), "e-", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 999999.49999999) {
    // 6 sig fig exponential notation, large
    if (dxx >= 9.9999949999999e15) {
      if (dxx >= 9.9999949999999e127) {
	if (dxx == INFINITY) {
	  *((uint32_t*)start) = *((uint32_t*)"inf");
	  return &(start[3]);
	} else if (dxx >= 9.9999949999999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9999949999999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9999949999999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9999949999999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9999949999999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9999949999999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999949999999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999949999999e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    double_bround5(dxx, banker_round8, &quotient, &remainder);
    start = memcpya(uint32_write1p5(start, quotient, remainder), "e+", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 0.99999949999999) {
    return double_write6(start, dxx);
  } else {
    // 6 sig fig decimal, no less than ~0.0001
    start = memcpya(start, "0.", 2);
    if (dxx < 9.9999949999999e-3) {
      dxx *= 100;
      start = memcpya(start, "00", 2);
    }
    if (dxx < 9.9999949999999e-2) {
      dxx *= 10;
      *start++ = '0';
    }
    return uint32_write6trunc(start, double_bround(dxx * 1000000, banker_round8));
  }
}

int32_t main(int32_t argc, char** argv) {
  gzFile dosefile = NULL;
  gzFile infofile = NULL;
  gzFile gz_outfile_pdat = NULL;
  FILE* outfile_pdat = NULL;
  FILE* outfile_pfam = NULL;
  char* outname = NULL;
  char* outname_end = NULL;
  unsigned char* wkspace_ua = NULL;
  uintptr_t dose_line_idx = 0;
  uintptr_t info_line_idx = 0;
  uintptr_t marker_ct = 0;
  uintptr_t sample_ct = 0;
  uintptr_t dosagebuf_width = 1;
  uintptr_t next_loadbuf_size = MINLINEBUFLEN - 1;
  intptr_t malloc_size_mb = 0;
  uint32_t dosename_param_idx = 0;
  uint32_t infoname_param_idx = 0;
  uint32_t gzip_pdat = 2;
  int32_t retval = 0;
  char numbuf[16];
#ifdef __APPLE__
  int32_t mib[2];
  size_t sztmp;
#else
#ifdef _WIN32
  MEMORYSTATUSEX memstatus;
#endif
#endif
  intptr_t default_alloc_mb;
  uintptr_t loadbuf_size;
  uintptr_t dosagebuf_size;
  uintptr_t sample_idx;
  uintptr_t marker_idx;
  uintptr_t marker_idx_start;
  uintptr_t ulii;
  double* dosagebuf;
  double* dptr;
  char* param_ptr;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  int64_t llxx;
  uint32_t param_idx;
  uint32_t pass_idx;
  uint32_t uii;
  unsigned char ucc;
  unsigned char ucc2;
  if (argc == 1) {
    goto main_ret_HELP;
  }
  for (param_idx = 1; param_idx < (uint32_t)argc; param_idx++) {
    if ((!strcmp(argv[param_idx], "--help")) || (!strcmp(argv[param_idx], "-help")) || (!strcmp(argv[param_idx], "-?")) || (!strcmp(argv[param_idx], "-h"))) {
      goto main_ret_HELP;
    }
  }

  if (argc > 11) {
    fputs("Error: Too many parameters.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  }
  for (param_idx = 1; param_idx < (uint32_t)argc; param_idx++) {
    if (argv[param_idx][0] != '-') {
      fputs("Error: Invalid parameter sequence.\n\n", stderr);
      goto main_ret_INVALID_CMDLINE_2;
    }
    param_ptr = &(argv[param_idx][1]);
    if (*param_ptr == '-') {
      // allow both single- and double-dash
      param_ptr++;
    }
    if ((!strcmp(param_ptr, "dose")) || (!strcmp(param_ptr, "d"))) {
      if (dosename_param_idx) {
        fputs("Error: Multiple instances of --dose.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      if ((uint32_t)argc == ++param_idx) {
	fputs("Error: Missing --dose parameter.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      dosename_param_idx = param_idx;
    } else if ((!strcmp(param_ptr, "info")) || (!strcmp(param_ptr, "i"))) {
      if (infoname_param_idx) {
        fputs("Error: Multiple instances of --info.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      if ((uint32_t)argc == ++param_idx) {
	fputs("Error: Missing --info parameter.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      infoname_param_idx = param_idx;
    } else if ((!strcmp(param_ptr, "out")) || (!strcmp(param_ptr, "o"))) {
      if (outname) {
        fputs("Error: Multiple instances of --out.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      if ((uint32_t)argc == ++param_idx) {
	fputs("Error: Missing --out parameter.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      uii = strlen(argv[param_idx]);
      outname = (char*)malloc(uii + 9);
      if (!outname) {
	goto main_ret_NOMEM;
      }
      memcpy(outname, argv[param_idx], uii);
      outname_end = &(outname[uii]);
    } else if (!strcmp(param_ptr, "no-gz")) {
      if (gzip_pdat != 2) {
        fputs("Error: Multiple instances of --{no-}gz.\n", stderr);
        goto main_ret_INVALID_CMDLINE_2;
      }
      gzip_pdat = 0;
    } else if (!strcmp(param_ptr, "gz")) {
      if ((uint32_t)argc == ++param_idx) {
	fputs("Error: Missing -gz parameter.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      if (gzip_pdat != 2) {
        fputs("Error: Multiple instances of --{no-}gz.\n", stderr);
        goto main_ret_INVALID_CMDLINE_2;
      }
      if (!strcmp(argv[param_idx], "0")) {
        gzip_pdat = 0;
      } else if (!strcmp(argv[param_idx], "1")) {
	gzip_pdat = 1;
      } else {
	fputs("Error: '-gz' must be followed by '0' or '1'.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
    } else if ((!strcmp(param_ptr, "memory")) || (!strcmp(param_ptr, "m"))) {
      if ((uint32_t)argc == ++param_idx) {
        fputs("Error: Missing --memory parameter.\n\n", stderr);
        goto main_ret_INVALID_CMDLINE_2;
      }
      if (scan_posintptr(argv[param_idx], (uintptr_t*)&malloc_size_mb)) {
	fprintf(stderr, "Error: Invalid --memory parameter '%s'.\n\n", argv[param_idx]);
        goto main_ret_INVALID_CMDLINE_2;
      }
      if (malloc_size_mb < WKSPACE_MIN_MB) {
	fprintf(stderr, "Error: Invalid --memory parameter '%s' (minimum %u).\n\n", argv[param_idx], WKSPACE_MIN_MB);
	goto main_ret_INVALID_CMDLINE_2;
      }
#ifndef __LP64__
      if (malloc_size_mb > 2047) {
        fputs("Error: --memory parameter too large for 32-bit build (max 2047).\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
#endif
    } else {
      fprintf(stderr, "Error: Invalid flag '%s'.\n\n", argv[param_idx]);
      goto main_ret_INVALID_CMDLINE_2;
    }
  }
  if (!dosename_param_idx) {
    fputs("Error: --dose is required.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  } else if (!infoname_param_idx) {
    fputs("Error: --info is required.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  }
  if (!outname) {
    outname = (char*)malloc(21);
    if (!outname) {
      goto main_ret_NOMEM;
    }
    memcpy(outname, "plink-dosage", 12);
    outname_end = &(outname[12]);
  }

#ifdef __APPLE__
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  sztmp = sizeof(int64_t);
  sysctl(mib, 2, &llxx, &sztmp, NULL, 0);
  llxx /= 1048576;
#else
#ifdef _WIN32
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
    printf("%" PRId64 " MB RAM detected; reserving %" PRIdPTR " MB for main workspace.\n", llxx, malloc_size_mb);
  } else {
    printf("Failed to calculate system memory.  Attempting to reserve %" PRIdPTR " MB.\n", malloc_size_mb);
  }
  wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576);
  while (!wkspace_ua) {
    malloc_size_mb = (malloc_size_mb * 3) / 4;
    if (malloc_size_mb < WKSPACE_MIN_MB) {
      malloc_size_mb = WKSPACE_MIN_MB;
    }
    wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576);
    if (wkspace_ua) {
      printf("Allocated %" PRIdPTR " MB successfully, after larger attempts failed.\n", malloc_size_mb);
    } else if (malloc_size_mb == WKSPACE_MIN_MB) {
      goto main_ret_NOMEM;
    }
  }
  wkspace = (unsigned char*)CACHEALIGN((uintptr_t)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = (malloc_size_mb * 1048576 - (uintptr_t)(wkspace - wkspace_ua)) & (~(CACHELINE - ONELU));

  // standardize strtod() behavior
  setlocale(LC_NUMERIC, "C");

  // okay, boilerplate done, real program can start

  if (gzopen_checked(&infofile, argv[infoname_param_idx], "rb")) {
    goto main_ret_OPEN_FAIL;
  }
  if (gzbuffer(infofile, 131072)) {
    goto main_ret_NOMEM;
  }
  if (gzopen_checked(&dosefile, argv[dosename_param_idx], "rb")) {
    goto main_ret_OPEN_FAIL;
  }
  if (gzbuffer(dosefile, 131072)) {
    goto main_ret_NOMEM;
  }
  memcpy(outname_end, ".pfam", 6);
  if (fopen_checked(&outfile_pfam, outname, "w")) {
    goto main_ret_OPEN_FAIL;
  }
  if (gzip_pdat) {
    memcpy(outname_end, ".pdat.gz", 9);
    if (gzopen_checked(&gz_outfile_pdat, outname, "wb")) {
      goto main_ret_OPEN_FAIL;
    }
    if (gzputs(gz_outfile_pdat, "SNP\tA1\tA2") == -1) {
      goto main_ret_WRITE_FAIL;
    }
  } else {
    memcpy(outname_end, ".pdat", 6);
    if (fopen_checked(&outfile_pdat, outname, "w")) {
      goto main_ret_OPEN_FAIL;
    }
    if (fputs_checked("SNP\tA1\tA2", outfile_pdat)) {
      goto main_ret_WRITE_FAIL;
    }
  }
  // Line loading buffer = min(half of workspace, 2 GB).
  // Rest of workspace stores dosages from bottom up.  FIDs and IIDs are not
  // actually kept in memory; instead, they're written immediately to both the
  // .pfam and .pdat files upon reading during the first pass.
  loadbuf_size = CACHEALIGN(wkspace_left / 2);
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size < MINLINEBUFLEN) {
    goto main_ret_NOMEM;
  }
  if (wkspace_alloc_c_checked(&loadbuf, loadbuf_size)) {
    goto main_ret_NOMEM;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  dosagebuf = (double*)wkspace_base;
  dosagebuf_size = wkspace_left / sizeof(double);

  // First pass: start by assuming entire file is loadable.  If we run out of
  // space, discard the last half of each loaded line; continue halving until
  // we reach the end of the .dose file (or we determine we have insufficient
  // memory for even a single .pdat line).
  // Subsequent passes, if necessary: we now know the number of samples and the
  // maximum .dose line length, so we can make efficient use of memory.

  // pass 1, read .dose and write .pfam
  while (1) {
    dose_line_idx++;
    if (!gzgets(dosefile, loadbuf, loadbuf_size)) {
      if (!gzeof(dosefile)) {
        goto main_ret_READ_FAIL;
      }
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        fprintf(stderr, "Error: Line %" PRIuPTR " of --dose file is pathologically long.\n", dose_line_idx);
        goto main_ret_INVALID_FORMAT;
      }
      goto main_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      goto main_update_max_line_1;
    }
    bufptr2 = bufptr;
    ucc = bufptr2[0];
    ucc2 = *(++bufptr2);
    while ((ucc != '-') || (ucc2 != '>')) {
      if (is_space_or_eoln(ucc2)) {
	fprintf(stderr, "Error: No '->' in first token of line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	goto main_ret_INVALID_FORMAT;
      }
      ucc = ucc2;
      ucc2 = *(++bufptr2);
    }
    if (bufptr2 == &(bufptr[1])) {
      fprintf(stderr, "Error: No FID in line %" PRIuPTR " of --dose file.\n", dose_line_idx);
      goto main_ret_INVALID_FORMAT;
    }
    bufptr3 = ++bufptr2;
    // may as well be paranoid
    ucc2 = *bufptr3;
    if (is_space_or_eoln(ucc2)) {
      fprintf(stderr, "Error: No IID in line %" PRIuPTR " of --dose file.\n", dose_line_idx);
      goto main_ret_INVALID_FORMAT;
    }
    do {
      ucc = ucc2;
      ucc2 = *(++bufptr3);
      if ((ucc == '-') && (ucc2 == '>')) {
	fprintf(stderr, "Error: Multiple '->'s in first token of line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	goto main_ret_INVALID_FORMAT;
      }
    } while (!is_space_or_eoln(ucc2));

    bufptr2[-2] = '\0';
    if (fputs_checked(bufptr, outfile_pfam)) {
      goto main_ret_WRITE_FAIL;
    }
    bufptr2[-1] = '\t';
    *bufptr3 = '\0';
    fputs(&(bufptr2[-1]), outfile_pfam);
    fputs("\t0\t0\t1\t-9\n", outfile_pfam);
    if (gzip_pdat) {
      if (gzputc(gz_outfile_pdat, '\t') == -1) {
	goto main_ret_WRITE_FAIL;
      }
      if (gzputs(gz_outfile_pdat, bufptr) == -1) {
	goto main_ret_WRITE_FAIL;
      }
      if (gzputs(gz_outfile_pdat, &(bufptr2[-1])) == -1) {
	goto main_ret_WRITE_FAIL;
      }
    } else {
      putc('\t', outfile_pdat);
      fputs(bufptr, outfile_pdat);
      fputs(&(bufptr2[-1]), outfile_pdat);
    }
    *bufptr3 = ucc2;
    bufptr = next_token(skip_initial_spaces(bufptr3));
    if (no_more_tokens_kns(bufptr)) {
      fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
      goto main_ret_INVALID_FORMAT;
    }
    if (!sample_ct) {
      do {
	if (marker_ct == dosagebuf_size) {
	  goto main_ret_NOMEM;
	}
	if (scan_double(bufptr, &(dosagebuf[marker_ct]))) {
	  fprintf(stderr, "Error: Invalid dosage value on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
	marker_ct++;
	bufptr = skip_initial_spaces(token_endnn(bufptr));
      } while (!is_eoln_kns(*bufptr));
      dosagebuf_width = marker_ct;
    } else {
      if ((sample_ct + 1) * dosagebuf_width > dosagebuf_size) {
	if (dosagebuf_width == 1) {
	  goto main_ret_NOMEM;
	}
	ulii = dosagebuf_width / 2;
        for (sample_idx = 1; sample_idx < sample_ct; sample_idx++) {
	  memcpy(&(dosagebuf[sample_idx * ulii]), &(dosagebuf[sample_idx * dosagebuf_width]), ulii * sizeof(double));
	}
	dosagebuf_width = ulii;
      }
      dptr = &(dosagebuf[sample_ct * dosagebuf_width]);
      for (marker_idx = 0; marker_idx < dosagebuf_width; marker_idx++) {
        if (is_eoln_kns(*bufptr)) {
          fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
        if (scan_double(bufptr, dptr)) {
	  fprintf(stderr, "Error: Invalid dosage value on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
	dptr++;
        bufptr = skip_initial_spaces(token_endnn(bufptr));
      }
      if (dosagebuf_width < marker_ct) {
	bufptr = next_token_mult(bufptr, marker_ct - dosagebuf_width - 1);
	if (no_more_tokens_kns(bufptr)) {
	  fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
        bufptr = skip_initial_spaces(token_endnn(bufptr));
      }
      if (!is_eoln_kns(*bufptr)) {
	fprintf(stderr, "Error: More tokens than expected on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	goto main_ret_INVALID_FORMAT;
      }
    }
    sample_ct++;
  main_update_max_line_1:
    ulii = strlen(bufptr) + (uintptr_t)(bufptr - loadbuf);
    if (ulii >= next_loadbuf_size) {
      next_loadbuf_size = ulii + 1;
    }
  }
  if (fclose_null(&outfile_pfam)) {
    goto main_ret_WRITE_FAIL;
  }

  numbuf[0] = '\t';
  fputs("Starting dose2plink ("
#ifdef __LP64__
"64"
#else
"32"
#endif
	"-bit).\n", stdout);

  // pass 1, read .info and write .pdat
  for (marker_idx = ~ZEROLU; marker_idx != dosagebuf_width;) {
    info_line_idx++;
    if (!gzgets(infofile, loadbuf, loadbuf_size)) {
      if (!gzeof(infofile)) {
	goto main_ret_READ_FAIL;
      }
      fputs("Error: --info file has fewer lines than expected.\n", stderr);
      goto main_ret_INVALID_FORMAT;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        fprintf(stderr, "Error: Line %" PRIuPTR " of --info file is pathologically long.\n", info_line_idx);
        goto main_ret_INVALID_FORMAT;
      }
      goto main_ret_NOMEM;
    }

    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      goto main_update_max_line_2;
    }
    if (marker_idx == ~ZEROLU) {
      if (strcmp_se(bufptr, "SNP", 3)) {
	fputs("Error: First field of --info file must be 'SNP'.\n", stderr);
	goto main_ret_INVALID_FORMAT;
      }
      bufptr = skip_initial_spaces(&(bufptr[3]));
      if (strcmp_se(bufptr, "Al1", 3)) {
	fputs("Error: Second field of --info file must be 'Al1'.\n", stderr);
	goto main_ret_INVALID_FORMAT;
      }
      bufptr = skip_initial_spaces(&(bufptr[3]));
      if (strcmp_se(bufptr, "Al2", 3)) {
	fputs("Error: Third field of --info file must be 'Al2'.\n", stderr);
	goto main_ret_INVALID_FORMAT;
      }
      if (gzip_pdat) {
	if (gzputc(gz_outfile_pdat, '\n') == -1) {
	  goto main_ret_WRITE_FAIL;
	}
      } else {
	if (putc_checked('\n', outfile_pdat)) {
	  goto main_ret_WRITE_FAIL;
	}
      }
      bufptr = &(bufptr[3]);
    } else {
      bufptr2 = token_endnn(bufptr);
      dptr = &(dosagebuf[marker_idx]);
      if (gzip_pdat) {
	if (!gzwrite(gz_outfile_pdat, bufptr, bufptr2 - bufptr)) {
	  goto main_ret_WRITE_FAIL;
	}
	bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
	  fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
	bufptr2 = token_endnn(bufptr);
	*(--bufptr) = '\t';
	if (!gzwrite(gz_outfile_pdat, bufptr, bufptr2 - bufptr)) {
	  goto main_ret_WRITE_FAIL;
	}
	bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
	  fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
	bufptr2 = token_endnn(bufptr);
	*(--bufptr) = '\t';
	if (!gzwrite(gz_outfile_pdat, bufptr, bufptr2 - bufptr)) {
	  goto main_ret_WRITE_FAIL;
	}
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	  *double_g_write(&(numbuf[1]), dptr[sample_idx * dosagebuf_width]) = '\0';
	  if (gzputs(gz_outfile_pdat, numbuf) == -1) {
	    goto main_ret_WRITE_FAIL;
	  }
	}
	if (gzputc(gz_outfile_pdat, '\n') == -1) {
	  goto main_ret_WRITE_FAIL;
	}
      } else {
	if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile_pdat)) {
	  goto main_ret_WRITE_FAIL;
	}
	bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
	  fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
	bufptr2 = token_endnn(bufptr);
	*(--bufptr) = '\t';
	if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile_pdat)) {
	  goto main_ret_WRITE_FAIL;
	}
	bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
	  fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	  goto main_ret_INVALID_FORMAT;
	}
	bufptr2 = token_endnn(bufptr);
	*(--bufptr) = '\t';
	if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile_pdat)) {
	  goto main_ret_WRITE_FAIL;
	}
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	  *double_g_write(&(numbuf[1]), dptr[sample_idx * dosagebuf_width]) = '\0';
	  fputs(numbuf, outfile_pdat);
	}
	if (putc_checked('\n', outfile_pdat)) {
	  goto main_ret_WRITE_FAIL;
	}
      }
    }
    marker_idx++;
  main_update_max_line_2:
    ulii = strlen(bufptr) + (uintptr_t)(bufptr - loadbuf);
    if (ulii >= next_loadbuf_size) {
      next_loadbuf_size = ulii + 1;
    }
  }

  if (dosagebuf_width < marker_ct) {
    pass_idx = 0;
    marker_idx_start = dosagebuf_width;
    wkspace_reset(loadbuf);
    loadbuf_size = CACHEALIGN(next_loadbuf_size);
    loadbuf = (char*)wkspace_alloc(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    dosagebuf_width = wkspace_left / (sizeof(double) * sample_ct);
    dosagebuf = (double*)wkspace_base;
    do {
      pass_idx++;
      printf("Pass %u complete (%" PRIuPTR "/%" PRIuPTR ").\n", pass_idx, marker_idx_start, marker_ct);
      if (marker_ct - marker_idx_start < dosagebuf_width) {
	dosagebuf_width = marker_ct - marker_idx_start;
      }
      dose_line_idx = 0;
      gzrewind(dosefile);
      // later passes: read .dose
      dptr = dosagebuf;
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
        do {
	  dose_line_idx++;
	  if (!gzgets(dosefile, loadbuf, loadbuf_size)) {
            goto main_ret_READ_FAIL;
	  }
	  bufptr = skip_initial_spaces(loadbuf);
	} while (is_eoln_kns(*bufptr));
        bufptr = next_token_mult(bufptr, 2 + marker_idx_start);
	for (marker_idx = 0; marker_idx < dosagebuf_width; marker_idx++) {
          if (scan_double(bufptr, dptr)) {
	    fprintf(stderr, "Error: Invalid dosage value on line %" PRIuPTR " of --dose file.\n", dose_line_idx);
	    goto main_ret_INVALID_FORMAT;
	  }
	  dptr++;
	  bufptr = skip_initial_spaces(token_endnn(bufptr));
	}
      }

      // later passes: read .info and write .pdat
      for (marker_idx = 0; marker_idx < dosagebuf_width;) {
	// marker_idx is actually an offset here
	info_line_idx++;
	if (!gzgets(infofile, loadbuf, loadbuf_size)) {
	  if (!gzeof(infofile)) {
	    goto main_ret_READ_FAIL;
	  }
          fputs("Error: --info file has fewer lines than expected.\n", stderr);
	  goto main_ret_INVALID_FORMAT;
	}
	if (!loadbuf[loadbuf_size - 1]) {
	  if (loadbuf_size == MAXLINEBUFLEN) {
	    fprintf(stderr, "Error: Line %" PRIuPTR " of --info file is pathologically long.\n", info_line_idx);
	    goto main_ret_INVALID_FORMAT;
	  }
	  goto main_ret_NOMEM;
	}
	bufptr = skip_initial_spaces(loadbuf);
	if (is_eoln_kns(*bufptr)) {
	  continue;
	}
	// er, this belongs in a function
	bufptr2 = token_endnn(bufptr);
	dptr = &(dosagebuf[marker_idx]);
	if (gzip_pdat) {
	  if (!gzwrite(gz_outfile_pdat, bufptr, bufptr2 - bufptr)) {
	    goto main_ret_WRITE_FAIL;
	  }
	  bufptr = skip_initial_spaces(bufptr2);
	  if (is_eoln_kns(*bufptr)) {
	    fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	    goto main_ret_INVALID_FORMAT;
	  }
	  bufptr2 = token_endnn(bufptr);
	  *(--bufptr) = '\t';
	  if (!gzwrite(gz_outfile_pdat, bufptr, bufptr2 - bufptr)) {
	    goto main_ret_WRITE_FAIL;
	  }
	  bufptr = skip_initial_spaces(bufptr2);
	  if (is_eoln_kns(*bufptr)) {
	    fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	    goto main_ret_INVALID_FORMAT;
	  }
	  bufptr2 = token_endnn(bufptr);
	  *(--bufptr) = '\t';
	  if (!gzwrite(gz_outfile_pdat, bufptr, bufptr2 - bufptr)) {
	    goto main_ret_WRITE_FAIL;
	  }
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    *double_g_write(&(numbuf[1]), dptr[sample_idx * dosagebuf_width]) = '\0';
	    if (gzputs(gz_outfile_pdat, numbuf) == -1) {
	      goto main_ret_WRITE_FAIL;
	    }
	  }
	  if (gzputc(gz_outfile_pdat, '\n') == -1) {
	    goto main_ret_WRITE_FAIL;
	  }
	} else {
	  if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile_pdat)) {
	    goto main_ret_WRITE_FAIL;
	  }
	  bufptr = skip_initial_spaces(bufptr2);
	  if (is_eoln_kns(*bufptr)) {
	    fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	    goto main_ret_INVALID_FORMAT;
	  }
	  bufptr2 = token_endnn(bufptr);
	  *(--bufptr) = '\t';
	  if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile_pdat)) {
	    goto main_ret_WRITE_FAIL;
	  }
	  bufptr = skip_initial_spaces(bufptr2);
	  if (is_eoln_kns(*bufptr)) {
	    fprintf(stderr, "Error: Fewer tokens than expected on line %" PRIuPTR " of --info file.\n", info_line_idx);
	    goto main_ret_INVALID_FORMAT;
	  }
	  bufptr2 = token_endnn(bufptr);
	  *(--bufptr) = '\t';
	  if (fwrite_checked(bufptr, bufptr2 - bufptr, outfile_pdat)) {
	    goto main_ret_WRITE_FAIL;
	  }
	  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	    *double_g_write(&(numbuf[1]), dptr[sample_idx * dosagebuf_width]) = '\0';
	    fputs(numbuf, outfile_pdat);
	  }
	  if (putc_checked('\n', outfile_pdat)) {
	    goto main_ret_WRITE_FAIL;
	  }
	}
	marker_idx++;
      }
      marker_idx_start += dosagebuf_width;
    } while (marker_idx_start < marker_ct);
  }
  if (gzip_pdat) {
    if (gzclose(gz_outfile_pdat) != Z_OK) {
      gz_outfile_pdat = NULL;
      goto main_ret_WRITE_FAIL;
    }
    gz_outfile_pdat = NULL;
  } else {
    if (fclose_null(&outfile_pdat)) {
      goto main_ret_WRITE_FAIL;
    }
  }
  printf("%s written (%" PRIuPTR " markers, %" PRIuPTR " samples).\n", outname, marker_ct, sample_ct);

  while (0) {
  main_ret_HELP:
    fputs(
"dose2plink "
#ifdef __LP64__
"64"
#else
"32"
#endif
"-bit (9 Nov 2014)  Christopher Chang (chrchang@alumni.caltech.edu)\n"
"Converts a MaCH/Minimac dosage file to a format usable by PLINK 1 --dosage.\n\n"
"Based on a Perl script by Sarah Medland; see\n"
"  http://www.genepi.qimr.edu.au/staff/sarahMe/dose2plink.html\n"
"for original source code and documentation.  This implementation is\n"
"GPLv3-licensed due to inclusion of some code from PLINK 1.9; if that's a\n"
"problem, Sarah's Perl script is unencumbered by copyleft.\n\n"
"Warning: all samples are assumed to be from males; you may need to postprocess\n"
"or replace the output .pfam file to correct the resulting errors.\n\n"
, stdout);
    disp_usage(stdout);
    retval = RET_HELP;
    break;
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
  main_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  main_ret_INVALID_CMDLINE_2:
    disp_usage(stderr);
    retval = RET_INVALID_CMDLINE;
    break;
  }
  gzclose_cond(dosefile);
  gzclose_cond(infofile);
  gzclose_cond(gz_outfile_pdat);
  fclose_cond(outfile_pdat);
  fclose_cond(outfile_pfam);
  free_cond(wkspace_ua);
  free_cond(outname);
  dispmsg(retval);
  return retval;
}
