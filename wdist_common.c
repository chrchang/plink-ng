#include "wdist_common.h"
#include "pigz.h"

const char errstr_fopen[] = "Error: Failed to open %s.\n";
const char errstr_append[] = "\nFor more information, try 'wdist --help [flag name]' or 'wdist --help | more'.\n";
const char errstr_thread_create[] = "\nError: Failed to create thread.\n";

char tbuf[MAXLINELEN * 4 + 256];

sfmt_t sfmt;

FILE* logfile = NULL;
char logbuf[MAXLINELEN]; // safe sprintf buffer, if one is needed
int32_t debug_on = 0;
int32_t log_failed = 0;
uintptr_t g_indiv_ct;
uint32_t g_thread_ct;

void logstr(const char* ss) {
  if (!debug_on) {
    if (fprintf(logfile, "%s", ss) < 0) {
      printf("\nWarning: Logging failure on:\n%s\nFurther logging will not be attempted in this run.\n", ss);
      log_failed = 1;
    }
  } else {
    if (log_failed) {
      printf("%s", ss);
      fflush(stdout);
    } else {
      if (fprintf(logfile, "%s", ss) < 0) {
        printf("\nError: Debug logging failure.  Dumping to standard output:\n%s", ss);
	log_failed = 1;
      } else {
	fflush(logfile);
      }
    }
  }
}

void logprint(const char* ss) {
  logstr(ss);
  fputs(ss, stdout);
}

void logprintb() {
  logstr(logbuf);
  fputs(logbuf, stdout);
}

int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    sprintf(logbuf, errstr_fopen, fname);
    logprintb();
    return -1;
  }
  return 0;
}

int32_t gzopen_checked(gzFile* target_ptr, const char* fname, const char* mode) {
  *target_ptr = gzopen(fname, mode);
  if (!(*target_ptr)) {
    sprintf(logbuf, errstr_fopen, fname);
    logprintb();
    return -1;
  }
  return 0;
}

// manually managed, very large stack
unsigned char* wkspace_base;
uintptr_t wkspace_left;

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

void wkspace_reset(void* new_base) {
  uintptr_t freed_bytes = wkspace_base - (unsigned char*)new_base;
  wkspace_base = (unsigned char*)new_base;
  wkspace_left += freed_bytes;
}

int32_t atoiz(char* ss, int32_t* sval) {
  int32_t ii = atoi(ss);
  ii = atoi(ss);
  if ((ii < 1) && ((*ss != '0') || (ss[1] != '\0'))) {
    return -1;
  }
  *sval = ii;
  return 0;
}

uint32_t strtoui32(char* ss, uint32_t* valp) {
  // one way to wrap strtoul() in a way that handles "0" and "4294967295"
  // nicely.  not high-performance; we probably want to write our own
  // string-to-number conversion routines at some point.
  uintptr_t ulii = strtoul(ss, NULL, 10);
  if (!ulii) {
    if (!memcmp(ss, "0", 2)) {
      *valp = 0;
      return 0;
    }
    return 1;
#ifdef __LP64__
  } else if (ulii > 4294967295LLU) {
#else
  } else if (ulii == ULONG_MAX) {
    if (!memcmp(ss, "4294967295", 11)) {
      *valp = 4294967295U;
      return 0;
    }
#endif
    return 1;
  } else {
    *valp = ulii;
    return 0;
  }
}

int32_t get_next_noncomment(FILE* fptr, char** lptr_ptr) {
  char* lptr;
  do {
    if (!fgets(tbuf, MAXLINELEN, fptr)) {
      return -1;
    }
    lptr = skip_initial_spaces(tbuf);
  } while (is_eoln_or_comment(*lptr));
  *lptr_ptr = lptr;
  return 0;
}

int32_t get_next_noncomment_excl(FILE* fptr, char** lptr_ptr, uintptr_t* marker_exclude, uintptr_t* marker_uidx_ptr) {
  while (!get_next_noncomment(fptr, lptr_ptr)) {
    if (!is_set(marker_exclude, *marker_uidx_ptr)) {
      return 0;
    }
    *marker_uidx_ptr += 1;
  }
  return -1;
}

char* item_end(char* sptr) {
  char cc;
  if (!sptr) {
    return NULL;
  }
  cc = *sptr;
  while (!is_space_or_eoln(cc)) {
    cc = *(++sptr);
  }
  return cc? sptr : NULL;
}

char* item_endl(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
  }
  return sptr;
}

void get_top_two(uint32_t* uint_arr, uint32_t uia_size, uint32_t* top_idx_ptr, uint32_t* second_idx_ptr) {
  uint32_t cur_idx = 2;
  uint32_t top_idx;
  uint32_t top_val;
  uint32_t second_idx;
  uint32_t second_val;
  uint32_t cur_val;
  if (uint_arr[1] > uint_arr[0]) {
    top_idx = 1;
  } else {
    top_idx = 0;
  }
  second_idx = 1 ^ top_idx;
  top_val = uint_arr[top_idx];
  second_val = uint_arr[second_idx];
  do {
    cur_val = uint_arr[cur_idx];
    if (cur_val > second_val) {
      if (cur_val > top_val) {
	second_val = top_val;
	second_idx = top_idx;
	top_val = cur_val;
	top_idx = cur_idx;
      } else {
	second_val = cur_val;
	second_idx = cur_idx;
      }
    }
  } while (++cur_idx < uia_size);
  *top_idx_ptr = top_idx;
  *second_idx_ptr = second_idx;
}

int32_t intlen(int32_t num) {
  int32_t retval;
  if (num < 0) {
    num = -num;
    retval = 2;
  } else {
    retval = 1;
  }
  while (num > 9) {
    num /= 10;
    retval++;
  }
  return retval;
}

int32_t strlen_se(char* ss) {
  int32_t val = 0;
  while (!is_space_or_eoln(*ss++)) {
    val++;
  }
  return val;
}

int32_t strcmp_se(char* s_read, const char* s_const, int32_t len) {
  return memcmp(s_read, s_const, len) || (!is_space_or_eoln(s_read[len]));
}

char* next_item(char* sptr) {
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

char* next_item_mult(char* sptr, uint32_t ct) {
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

// number-to-string encoders

static const char digit2_table[] = {
  "0001020304050607080910111213141516171819"
  "2021222324252627282930313233343536373839"
  "4041424344454647484950515253545556575859"
  "6061626364656667686970717273747576777879"
  "8081828384858687888990919293949596979899"};

char* uint32_write(char* start, uint32_t uii) {
  // Memory-efficient fast integer writer.  (You can do a bit better sometimes
  // by using a larger lookup table, but on average I doubt that pays off.)
  //
  // Originally the arguments were in the other order (was trying to follow
  // Google's "inputs first, than outputs" coding style guidelines), but then I
  // realized that chained invocation of this function is much easier to read
  // if I make the target buffer the first argument.
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      *start = '0' + uii;
      return &(start[1]);
    } else if (uii >= 100) {
      quotient = uii / 100;
      *start++ = '0' + quotient;
      uii -= quotient * 100;
    }
    return memcpya(start, &(digit2_table[uii * 2]), 2);
  } else if (uii < 10000000) {
    if (uii >= 100000) {
      if (uii < 1000000) {
	goto uint32_write_6;
      }
      quotient = uii / 1000000;
      *start++ = '0' + quotient;
      goto uint32_write_6b;
    } else if (uii < 10000) {
      goto uint32_write_4;
    }
    quotient = uii / 10000;
    *start++ = '0' + quotient;
  } else {
    if (uii >= 100000000) {
      quotient = uii / 100000000;
      if (uii >= 1000000000) {
	start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      } else {
	*start++ = '0' + quotient;
      }
      uii -= 100000000 * quotient;
    }
    quotient = uii / 1000000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  uint32_write_6b:
    uii -= 1000000 * quotient;
  uint32_write_6:
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 10000 * quotient;
 uint32_write_4:
  quotient = uii / 100;
  uii -= 100 * quotient;
  return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
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

static inline void uint32_write8(char* start, uint32_t uii) {
  uint32_t quotient = uii / 1000000;
  uint32_write6(memcpya(start, &(digit2_table[quotient * 2]), 2), uii - 1000000 * quotient);
}

char* int64_write(char* start, int64_t llii) {
  int64_t top_digits;
  uint32_t bottom_eight;
  uint32_t middle_eight;
  if (llii < 0) {
    if (llii < -9223372036854775807LL) {
      // special case, can't be represented positive
      return memcpya(start, "-9223372036854775808", 20);
    }
    *start++ = '-';
    llii = -llii;
  }
  if (llii <= 4294967295LL) {
    return uint32_write(start, (uint32_t)llii);
  }
  top_digits = llii / 100000000LL;
  bottom_eight = (uint32_t)(llii - (top_digits * 100000000));
  if (top_digits <= 4294967295LL) {
    start = uint32_write(start, (uint32_t)top_digits);
    uint32_write8(start, bottom_eight);
    return &(start[8]);
  }
  llii = top_digits / 100000000LL;
  middle_eight = (uint32_t)(top_digits - (llii * 100000000));
  start = uint32_write(start, (uint32_t)llii);
  uint32_write8(start, middle_eight);
  uint32_write8(&(start[8]), bottom_eight);
  return &(start[16]);
}

char* uint32_writew7(char* start, uint32_t uii) {
  // Minimum field width 7.
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      memset(start, 32, 6);
      start[6] = '0' + uii;
      return &(start[7]);
    } else if (uii < 100) {
      memset(start, 32, 5);
    } else {
      memset(start, 32, 4);
      quotient = uii / 100;
      start[4] = '0' + quotient;
      uii -= quotient * 100;
    }
    return memcpya(&(start[5]), &(digit2_table[uii * 2]), 2);
  } else if (uii < 10000000) {
    if (uii >= 100000) {
      if (uii >= 1000000) {
	quotient = uii / 1000000;
	*start++ = '0' + quotient;
	goto uint32_writew7_6b;
      }
      *start++ = ' ';
      goto uint32_writew7_6;
    } else if (uii < 100000) {
      start = memseta(start, 32, 2);
      quotient = uii / 10000;
      *start++ = '0' + quotient;
    } else {
      start = memseta(start, 32, 3);
      goto uint32_writew7_4;
    }
  } else {
    if (uii >= 100000000) {
      quotient = uii / 100000000;
      if (uii >= 1000000000) {
	start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      } else {
	*start++ = '0' + quotient;
      }
      uii -= 100000000 * quotient;
    }
    quotient = uii / 1000000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  uint32_writew7_6b:
    uii -= 1000000 * quotient;
  uint32_writew7_6:
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 10000 * quotient;
 uint32_writew7_4:
  quotient = uii / 100;
  uii -= 100 * quotient;
  return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
}

char* uint32_writew8(char* start, uint32_t uii) {
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      memset(start, 32, 7);
      start[7] = '0' + uii;
      return &(start[8]);
    } else if (uii < 100) {
      memset(start, 32, 6);
    } else {
      memset(start, 32, 5);
      quotient = uii / 100;
      start[5] = '0' + quotient;
      uii -= quotient * 100;
    }
    return memcpya(&(start[6]), &(digit2_table[uii * 2]), 2);
  } else if (uii < 10000000) {
    if (uii >= 100000) {
      if (uii < 1000000) {
	start = memseta(start, 32, 2);
	goto uint32_writew8_6;
      }
      quotient = uii / 1000000;
      *start = ' ';
      start[1] = '0' + quotient;
      start += 2;
      goto uint32_writew8_6b;
    } else if (uii < 10000) {
      start = memseta(start, 32, 4);
      goto uint32_writew8_4;
    }
    memset(start, 32, 3);
    quotient = uii / 10000;
    start[3] = '0' + quotient;
    start += 4;
  } else {
    if (uii >= 100000000) {
      quotient = uii / 100000000;
      if (uii >= 1000000000) {
	start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      } else {
	*start++ = '0' + quotient;
      }
      uii -= 100000000 * quotient;
    }
    quotient = uii / 1000000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  uint32_writew8_6b:
    uii -= 1000000 * quotient;
  uint32_writew8_6:
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 10000 * quotient;
 uint32_writew8_4:
  quotient = uii / 100;
  uii -= 100 * quotient;
  return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
}

char* uint32_writew10(char* start, uint32_t uii) {
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      memset(start, 32, 9);
      start[9] = '0' + uii;
      return &(start[10]);
    } else if (uii < 100) {
      memset(start, 32, 8);
    } else {
      memset(start, 32, 7);
      quotient = uii / 100;
      start[7] = '0' + quotient;
      uii -= quotient * 100;
    }
    return memcpya(&(start[8]), &(digit2_table[uii * 2]), 2);
  } else if (uii < 10000000) {
    if (uii >= 100000) {
      if (uii < 1000000) {
	start = memseta(start, 32, 4);
	goto uint32_writew10_6;
      }
      quotient = uii / 1000000;
      memset(start, 32, 3);
      start[3] = '0' + quotient;
      start += 4;
      goto uint32_writew10_6b;
    } else if (uii < 10000) {
      start = memseta(start, 32, 6);
      goto uint32_writew10_4;
    }
    memset(start, 32, 5);
    quotient = uii / 10000;
    start[5] = '0' + quotient;
    start += 6;
  } else {
    if (uii >= 100000000) {
      quotient = uii / 100000000;
      if (uii >= 1000000000) {
	memcpy(start, &(digit2_table[quotient * 2]), 2);
      } else {
	*start = ' ';
	start[1] = '0' + quotient;
      }
      uii -= 100000000 * quotient;
    } else {
      memset(start, 32, 2);
    }
    quotient = uii / 1000000;
    memcpy(&(start[2]), &(digit2_table[quotient * 2]), 2);
    start += 4;
  uint32_writew10_6b:
    uii -= 1000000 * quotient;
  uint32_writew10_6:
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 10000 * quotient;
 uint32_writew10_4:
  quotient = uii / 100;
  uii -= 100 * quotient;
  return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
}

static inline char* uint32_write4trunc(char* start, uint32_t uii) {
  // Given 0 < uii < 10000, writes uii without *trailing* zeroes.  (I.e. this
  // is for floating-point encoder use.)
  uint32_t quotient = uii / 100;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  uii -= 100 * quotient;
  if (uii) {
    start += 2;
    memcpy(start, &(digit2_table[uii * 2]), 2);
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
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

static inline char* uint32_write1p3(char* start, uint32_t quotient, uint32_t remainder) {
  // quotient = (int32_t)dxx;
  // remainder = ((int32_t)(dxx * 1000)) - (quotient * 1000);
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 10;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  remainder -= 10 * quotient;
  if (remainder) {
    start[2] = '0' + remainder;
    return &(start[3]);
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
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

char* double_write4(char* start, double dxx) {
  // 4 sig fig number, 0.9995 <= dxx < 9999.5
  uint32_t uii;
  uint32_t quotient;
  if (dxx < 99.995) {
    if (dxx < 9.9995) {
      dxx += 0.0005;
      quotient = (int32_t)dxx;
      return uint32_write1p3(start, quotient, ((int32_t)(dxx * 1000)) - (quotient * 1000));
    }
    dxx += 0.005;
    uii = (int32_t)dxx;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    uii = ((int32_t)(dxx * 100)) - (uii * 100);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    memcpy(start, &(digit2_table[uii * 2]), 2);
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 999.95) {
    dxx += 0.05;
    uii = (int32_t)dxx;
    quotient = uii / 100;
    *start = '0' + quotient;
    quotient = uii - 100 * quotient;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    uii = ((int32_t)(dxx * 10)) - (uii * 10);
    if (!uii) {
      return start;
    }
    *start = '.';
    start[1] = '0' + uii;
    return &(start[2]);
  } else {
    uint32_write4(start, (int32_t)(dxx + 0.5));
    return &(start[4]);
  }
}

char* double_e_write(char* start, double dxx) {
  uint32_t xp10 = 0;
  uint32_t uii;
  char sign;
  if (dxx != dxx) {
    // do this first to avoid generating exception
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx >= 9.9999995e-1) {
    if (dxx >= 9.9999995e7) {
      if (dxx >= 9.9999995e127) {
	if (dxx == INFINITY) {
	  *((uint32_t*)start) = *((uint32_t*)"inf");
	  return &(start[3]);
	} else if (dxx >= 9.9999995e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9999995e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9999995e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9999995e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
      if (dxx >= 9.9999995e7) {
	dxx *= 1.0e-8;
	xp10 |= 8;
      }
    }
    if (dxx >= 9.9999995e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999995e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999995) {
      dxx *= 1.0e-1;
      xp10++;
    }
    sign = '+';
  } else {
    if (dxx < 9.9999995e-8) {
      // general case
      if (dxx < 9.9999995e-128) {
	if (dxx == 0.0) {
	  return memcpya(start, "0.000000e+00", 12);
	}
	if (dxx < 9.9999995e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9999995e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9999995e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9999995e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
      if (dxx < 9.9999995e-8) {
	dxx *= 100000000;
	xp10 |= 8;
      }
    }
    if (dxx < 9.999995e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999995e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999995e-1) {
      dxx *= 10;
      xp10++;
    }
    sign = '-';
  }
  dxx += 0.0000005;
  uii = (int32_t)dxx;
  *start++ = '0' + uii;
  *start++ = '.';
  uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
  uint32_write6(start, uii);
  start += 6;
  *start++ = 'e';
  *start++ = sign;
  if (xp10 >= 100) {
    uii = xp10 / 100;
    *start++ = '0' + uii;
    xp10 -= uii * 100;
  }
  return memcpya(start, &(digit2_table[xp10 * 2]), 2);
}

char* float_e_write(char* start, float dxx) {
  uint32_t xp10 = 0;
  uint32_t uii;
  char sign;
  if (dxx != dxx) {
    // do this first to avoid generating exception
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx >= 9.9999995e-1) {
    if (dxx >= 9.9999995e15) {
      if (dxx == INFINITY) {
	*((uint32_t*)start) = *((uint32_t*)"inf");
	return &(start[3]);
      } else if (dxx >= 9.9999995e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      } else {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9999995e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9999995e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999995e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999995) {
      dxx *= 1.0e-1;
      xp10++;
    }
    sign = '+';
  } else {
    if (dxx < 9.9999995e-16) {
      if (dxx == 0.0) {
	return memcpya(start, "0.000000e+00", 12);
      } else if (dxx < 9.9999995e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      } else {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9999995e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.999995e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999995e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999995e-1) {
      dxx *= 10;
      xp10++;
    }
    sign = '-';
  }
  dxx += 0.0000005;
  uii = (int32_t)dxx;
  *start++ = '0' + uii;
  *start++ = '.';
  uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
  uint32_write6(start, uii);
  start += 6;
  *start++ = 'e';
  *start++ = sign;
  return memcpya(start, &(digit2_table[xp10 * 2]), 2);
}

char* double_f_writew6(char* start, double dxx) {
  int64_t llii;
  uint32_t uii;
  if (dxx != dxx) {
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 9.9999995) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.9999995) {
	goto double_f_writew6_10;
      }
    } else {
      *start++ = ' ';
    }
    dxx += 0.0000005;
    uii = ((int32_t)dxx);
    *start++ = '0' + uii;
    uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
  double_f_writew6_dec:
    *start++ = '.';
    uint32_write6(start, uii);
    return &(start[6]);
  }
 double_f_writew6_10:
  dxx += 0.0000005;
#ifndef __LP64__
  if (dxx < 2147.375) {
    uii = (int32_t)dxx;
    start = uint32_write(start, uii);
    uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
    goto double_f_writew6_dec;
  }
#endif
  // 2 ^ 63 int64_t max, divided by 1000000
  if (dxx <= 9223372036854.75) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (uint32_t)(((int64_t)(dxx * 1000000)) - (llii * 1000000));
    goto double_f_writew6_dec;
  } else if (dxx < 9223372036854775808.0) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (int32_t)((dxx - ((double)llii)) * 1000000);
    goto double_f_writew6_dec;
  }
  if (dxx == INFINITY) {
    *((uint32_t*)start) = *((uint32_t*)"inf");
    return &(start[3]);
  }
  // don't worry about optimizing %f on huge-ass finite numbers for now, since
  // it should be irrelevant for PLINK
  start += sprintf(start, "%.6f", dxx);
  return start;
}

char* double_f_writew74(char* start, double dxx) {
  int64_t llii;
  uint32_t uii;
  uint32_t quotient;
  if (dxx != dxx) {
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 9.99995) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.99995) {
	goto double_f_writew74_10;
      }
    } else {
      *start++ = ' ';
    }
    dxx += 0.00005;
    uii = ((int32_t)dxx);
    *start++ = '0' + uii;
    uii = ((int32_t)(dxx * 10000)) - (uii * 10000);
  double_f_writew74_dec:
    *start++ = '.';
    quotient = uii / 100;
    uii -= 100 * quotient;
    return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
  }
 double_f_writew74_10:
  dxx += 0.00005;
#ifndef __LP64__
  if (dxx < 214748.25) {
    uii = (int32_t)dxx;
    start = uint32_write(start, uii);
    uii = ((int32_t)(dxx * 10000)) - (uii * 10000);
    goto double_f_writew74_dec;
  }
#endif
  // 2 ^ 63 int64_t max, divided by 10000
  if (dxx <= 922337203685477.5) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (uint32_t)(((int64_t)(dxx * 10000)) - (llii * 10000));
    goto double_f_writew74_dec;
  } else if (dxx < 9223372036854775808.0) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (int32_t)((dxx - ((double)llii)) * 10000);
    goto double_f_writew74_dec;
  }
  if (dxx == INFINITY) {
    *((uint32_t*)start) = *((uint32_t*)"inf");
    return &(start[3]);
  }
  start += sprintf(start, "%.4f", dxx);
  return start;
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

char* double_g_writewx4(char* start, double dxx, uint32_t min_width) {
  // assumes min_width >= 4.
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t uii;
  if (dxx != dxx) {
    memcpy(memseta(start, 32, min_width - 4), " nan", 4);
    return &(start[min_width]);
  } else if (dxx < 0) {
    *wpos++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9995e-5) {
    // 4 sig fig exponential notation, small
    if (dxx < 9.9995e-16) {
      if (dxx < 9.9995e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.9995e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9995e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9995e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9995e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9995e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9995e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9995e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9995e-1) {
      dxx *= 10;
      xp10++;
    }
    dxx += 0.0005;
    uii = (int32_t)dxx;
    wpos = uint32_write1p3(wpos, uii, ((int32_t)(dxx * 1000)) - (uii * 1000));
    uii = wpos - wbuf;
    if (xp10 >= 100) {
      if (uii < min_width - 5) {
	memcpy(memseta(start, 32, min_width - 5 - uii), wbuf, uii);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, uii);
      }
      uii = xp10 / 100;
      start = memcpyax(start, "e-", 2, '0' + uii);
      xp10 -= 100 * uii;
    } else {
      if (uii < min_width - 4) {
	memcpy(memseta(start, 32, min_width - 4 - uii), wbuf, uii);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, uii);
      }
      start = memcpya(start, "e-", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 9999.5) {
    // 4 sig fig exponential notation, large
    if (dxx >= 9.9995e15) {
      if (dxx >= 9.9995e127) {
	if (dxx == INFINITY) {
	  start = memseta(start, 32, min_width - 4);
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.9995e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9995e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9995e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9995e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9995e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9995e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9995e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9995e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    dxx += 0.0005;
    uii = (int32_t)dxx;
    wpos = uint32_write1p3(wpos, uii, ((int32_t)(dxx * 1000)) - (uii * 1000));
    uii = wpos - wbuf;
    if (xp10 >= 100) {
      if (uii < min_width - 5) {
	memcpy(memseta(start, 32, min_width - 5 - uii), wbuf, uii);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, uii);
      }
      uii = xp10 / 100;
      start = memcpyax(start, "e+", 2, '0' + uii);
      xp10 -= 100 * uii;
    } else {
      if (uii < min_width - 4) {
	memcpy(memseta(start, 32, min_width - 4 - uii), wbuf, uii);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, uii);
      }
      start = memcpya(start, "e+", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else {
    if (dxx >= 0.99995) {
      wpos = double_write4(wpos, dxx);
    } else {
      // 4 sig fig decimal, no less than ~0.0001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.9995e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.9995e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uint32_write4trunc(wpos, (int32_t)((dxx * 10000) + 0.5));
    }
    uii = wpos - wbuf;
    if (uii < min_width) {
      memcpy(memseta(start, 32, min_width - uii), wbuf, uii);
      return &(start[min_width]);
    } else {
      return memcpya(start, wbuf, uii);
    }
  }
}

char* chrom_print_human(char* buf, uint32_t num) {
  uint32_t n10;
  if (num < 10) {
    *buf = '0' + num;
    return &(buf[1]);
  } else if (num < 23) {
    n10 = num / 10;
    *buf = '0' + n10;
    buf[1] = '0' + (num - 10 * n10);
    return &(buf[2]);
  } else if (num < 25) {
    // X is 24th letter of alphabet, and 23rd chromosome
    *buf = 'A' + num;
    return &(buf[1]);
  } else if (num == 25) {
    memcpy(buf, "XY", 2);
    return &(buf[2]);
  } else {
    memcpy(buf, "MT", 2);
    return &(buf[2]);
  }
}

void chrom_print_human_terminate(char* buf, uint32_t num) {
  char* ss = chrom_print_human(buf, num);
  *ss = '\0';
}

void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* pre_shiftp, uint32_t* post_shiftp, uint32_t* incrp) {
  // Enables fast integer division by a constant not known until runtime.  See
  // http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html .
  // Assumes divisor is not zero, of course.
  uint32_t down_multiplier = 0;
  uint32_t down_exponent = 0;
  uint32_t has_magic_down = 0;
  uint32_t quotient;
  uint32_t remainder;
  uint32_t ceil_log_2_d;
  uint32_t exponent;
  uint32_t uii;
  if (divisor & (divisor - 1)) {
    quotient = 2147483648U / divisor;
    remainder = 2147483648U - (quotient * divisor);
    ceil_log_2_d = 32 - __builtin_clz(divisor);
    for (exponent = 0; ; exponent++) {
      if (remainder >= divisor - remainder) {
        quotient = quotient * 2 + 1;
	remainder = remainder * 2 - divisor;
      } else {
	quotient = quotient * 2;
	remainder = remainder * 2;
      }
      if ((exponent >= ceil_log_2_d) || (divisor - remainder) <= (1U << exponent)) {
	break;
      }
      if ((!has_magic_down) && (remainder <= (1U << exponent))) {
	has_magic_down = 1;
	down_multiplier = quotient;
	down_exponent = exponent;
      }
    }
    if (exponent < ceil_log_2_d) {
      *multp = quotient + 1;
      *pre_shiftp = 0;
      *post_shiftp = 32 + exponent;
      *incrp = 0;
      return;
    } else if (divisor & 1) {
      *multp = down_multiplier;
      *pre_shiftp = 0;
      *post_shiftp = 32 + down_exponent;
      *incrp = 1;
      return;
    } else {
      *pre_shiftp = __builtin_ctz(divisor);
      magic_num(divisor >> (*pre_shiftp), multp, &uii, post_shiftp, incrp);
      return;
    }
  } else {
    // power of 2
    *multp = 1;
    *pre_shiftp = 0;
    *post_shiftp = __builtin_ctz(divisor);
    *incrp = 0;
  }
}

void set_bit(uintptr_t* bit_arr, uint32_t loc, uintptr_t* bit_set_ct_ptr) {
  uint32_t maj = loc / BITCT;
  uintptr_t min = ONELU << (loc % BITCT);
  if (!(bit_arr[maj] & min)) {
    bit_arr[maj] |= min;
    *bit_set_ct_ptr += 1;
  }
}

void set_bit_sub(uintptr_t* bit_arr, uint32_t loc, uintptr_t* bit_unset_ct_ptr) {
  uint32_t maj = loc / BITCT;
  uintptr_t min = ONELU << (loc % BITCT);
  if (!(bit_arr[maj] & min)) {
    bit_arr[maj] |= min;
    *bit_unset_ct_ptr -= 1;
  }
}

void clear_bit(uintptr_t* exclude_arr, uint32_t loc, uintptr_t* include_ct_ptr) {
  uint32_t maj = loc / BITCT;
  uintptr_t min = ONELU << (loc % BITCT);
  if (exclude_arr[maj] & min) {
    exclude_arr[maj] -= min;
    *include_ct_ptr += 1;
  }
}

// unsafe if you don't know there's another included marker or person remaining
int32_t next_non_set_unsafe(uintptr_t* exclude_arr, uint32_t loc) {
  uint32_t idx = loc / BITCT;
  uintptr_t ulii;
  exclude_arr = &(exclude_arr[idx]);
  ulii = (~(*exclude_arr)) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    idx++;
  } while (*(++exclude_arr) == ~ZEROLU);
  return (idx * BITCT) + CTZLU(~(*exclude_arr));
}

int32_t next_non_set(uintptr_t* exclude_arr, uint32_t loc, uint32_t ceil) {
  // safe version.  ceil >= 1.
  uint32_t idx = loc / BITCT;
  uint32_t max_idx;
  uintptr_t ulii;
  exclude_arr = &(exclude_arr[idx]);
  ulii = (~(*exclude_arr)) >> (loc % BITCT);
  if (ulii) {
    return MINV(loc + CTZLU(ulii), ceil);
  }
  max_idx = (ceil - 1) / BITCT;
  do {
    if ((++idx) > max_idx) {
      return ceil;
    }
  } while (*(++exclude_arr) == ~ZEROLU);
  return MINV((idx * BITCT) + CTZLU(~(*exclude_arr)), ceil);
}

int32_t next_set_unsafe(uintptr_t* include_arr, uint32_t loc) {
  uint32_t idx = loc / BITCT;
  uintptr_t ulii;
  include_arr = &(include_arr[idx]);
  ulii = (*include_arr) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    idx++;
  } while (*(++include_arr) == 0);
  return (idx * BITCT) + CTZLU(*include_arr);
}

void fill_vec_55(uintptr_t* vec, uint32_t ct) {
  uint32_t ctl = 2 * ((ct + (BITCT - 1)) / BITCT);
  uint32_t rem = ct & (BITCT - 1);
  uintptr_t* second_to_last = &(vec[ctl - 2]);
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* vecp = (__m128i*)vec;
  __m128i* vec_end = (__m128i*)(&(vec[ctl]));
  do {
    *vecp++ = m1;
  } while (vecp < vec_end);
#else
  uintptr_t* vec_end = &(vec[ctl]);
  do {
    *vec++ = FIVEMASK;
  } while (vec < vec_end);
#endif
  if (rem > BITCT2) {
    second_to_last[1] &= (~ZEROLU) >> ((BITCT - rem) * 2);
  } else if (rem) {
    *second_to_last &= (~ZEROLU) >> ((BITCT2 - rem) * 2);
    second_to_last[1] = 0;
  }
}

const char acgtarr[] = "ACGT";

void indiv_delim_convert(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char oldc, char newc) {
  // assumes there is exactly one delimiter to convert per name
  uintptr_t indiv_uidx = 0;
  uintptr_t indiv_idx = 0;
  char* nptr;
  for (; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    nptr = (char*)memchr(&(person_ids[indiv_uidx * max_person_id_len]), (unsigned char)oldc, max_person_id_len);
    *nptr = newc;
    indiv_uidx++;
  }
}

// human: 22, X, Y, XY, MT
// cow: 29, X, Y
// dog: 38, X, Y, XY
// horse: 31, X, Y
// mouse: 19, X, Y
// rice: 12 (haploid, not supported for now)
// sheep: 26, X, Y
// const uint64_t species_def_chrom_mask[] = {0x027fffffLLU, 0x3fffffffLLU, 0x27fffffffffLLU, 0xffffffffLLU, 0x000fffffLLU, 0LLU, 0x07ffffffLLU};
// const uint64_t species_def_chrom_mask[] = {0x07ffffffLLU, 0xffffffffLLU, 0x3ffffffffffLLU, 0x3ffffffffLLU, 0x003fffffLLU, 0LLU, 0x1fffffffLLU};
const uint64_t species_autosome_mask[] = {0x007ffffeLLU, 0x3ffffffeLLU, 0x7ffffffffeLLU, 0xfffffffeLLU, 0x000ffffeLLU, 0LLU, 0x07fffffeLLU};
// const uint64_t species_valid_chrom_mask[] = {0x3c0007ffffffLLU, 0xc00ffffffffLLU, 0x1fffffffffffLLU, 0xc00ffffffffLLU, 0xc00000fffffLLU, 0LLU, 0xc001fffffffLLU};
const uint64_t species_valid_chrom_mask[] = {0x07ffffffLLU, 0xffffffffLLU, 0x3ffffffffffLLU, 0x3ffffffffLLU, 0x003fffffLLU, 0LLU, 0x1fffffffLLU};
// const uint64_t species_valid_chrom_mask[] = {0x3c0007ffffffLLU, 0xc00ffffffffLLU, 0x1fffffffffffLLU, 0xc00ffffffffLLU, 0xc00000fffffLLU, 0x00001fffLLU, 0xc001fffffffLLU}
// const char species_regchrom_ct_p1[] = {23, 30, 39, 40, 20, 13, 27};
const char species_x_code[] = {23, 30, 39, 32, 20, -1, 27};
const char species_y_code[] = {24, 31, 40, 33, 21, -1, 28};
const char species_xy_code[] = {25, -1, 41, -1, -1, -1, -1};
const char species_mt_code[] = {26, -1, -1, -1, -1, -1, -1};
const char species_max_code[] = {26, 31, 41, 33, 21, 12, 28};
const uint64_t species_haploid_mask[] = {0x05800000LLU, 0xc0000000LLU, 0x18000000000LLU, 0x300000000LLU, 0x00300000LLU, 0x00001fffLLU, 0x18000000LLU};
char species_singulars[][7] = {"person", "animal", "animal", "animal", "animal", "plant", "animal"};
char species_plurals[][8] = {"people", "animals", "animals", "animals", "animals", "plants", "animals"};

char* species_singular = NULL;
char* species_plural = NULL;

int32_t marker_code_raw(char* sptr) {
  // any character <= ' ' is considered a terminator
  int32_t ii;
  if (sptr[1] > ' ') {
    if (sptr[2] > ' ') {
      return -1;
    }
    if ((sptr[0] == 'X') || (sptr[0] == 'x')) {
      if ((sptr[1] == 'Y') || (sptr[1] == 'y')) {
	return CHROM_XY;
      }
      return -1;
    }
    if ((sptr[0] == 'M') || (sptr[0] == 'm')) {
      if ((sptr[1] == 'T') || (sptr[1] == 't')) {
	return CHROM_MT;
      }
      return -1;
    }
    if ((sptr[0] >= '0') && (sptr[0] <= '9')) {
      if ((sptr[1] >= '0') && (sptr[1] <= '9')) {
        ii = ((sptr[0] - '0') * 10 + (sptr[1] - '0'));
	if (ii < MAX_POSSIBLE_CHROM) {
	  return ii;
	} else {
	  return -1;
	}
      } else {
	return -1;
      }
    }
  } else if ((sptr[0] >= '0') && (sptr[0] <= '9')) {
    return (sptr[0] - '0');
  } else if ((sptr[0] == 'X') || (sptr[0] == 'x')) {
    return CHROM_X;
  } else if ((sptr[0] == 'Y') || (sptr[0] == 'y')) {
    return CHROM_Y;
  } else if ((sptr[0] == 'M') || (sptr[0] == 'm')) {
    return CHROM_MT;
  }
  return -1;
}

int32_t marker_code(uint32_t species, char* sptr) {
  // does not require string to be null-terminated, and does not perform
  // exhaustive error-checking
  int32_t ii = marker_code_raw(sptr);
  if (ii >= MAX_POSSIBLE_CHROM) {
    switch (ii) {
    case CHROM_X:
      ii = species_x_code[species];
      break;
    case CHROM_Y:
      ii = species_y_code[species];
      break;
    case CHROM_XY:
      ii = species_xy_code[species];
      break;
    case CHROM_MT:
      ii = species_mt_code[species];
    }
  } else if ((ii == -1) || (!(species_valid_chrom_mask[species] & (1LLU << ii)))) {
    return -1;
  }
  return ii;
}

int32_t marker_code2(uint32_t species, char* sptr, uint32_t slen) {
  char* s_end = &(sptr[slen]);
  char tmpc = *s_end;
  int32_t retval;
  *s_end = ' ';
  retval = marker_code(species, sptr);
  *s_end = tmpc;
  return retval;
}

int32_t get_marker_chrom(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx) {
  uint32_t* marker_binsearch = chrom_info_ptr->chrom_file_order_marker_idx;
  int32_t chrom_min = 0;
  int32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t chrom_cur;
  while (chrom_ct - chrom_min > 1) {
    chrom_cur = (chrom_ct + chrom_min) / 2;
    if (marker_binsearch[chrom_cur] > marker_uidx) {
      chrom_ct = chrom_cur;
    } else {
      chrom_min = chrom_cur;
    }
  }
  return chrom_info_ptr->chrom_file_order[chrom_min];
}

void refresh_chrom_info(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx, uint32_t set_hh_missing, uint32_t is_all_nonmale, uint32_t* chrom_end_ptr, uint32_t* chrom_fo_idx_ptr, uint32_t* is_x_ptr, uint32_t* is_haploid_ptr) {
  uint32_t species = chrom_info_ptr->species;
  int32_t chrom_idx;
  *chrom_end_ptr = chrom_info_ptr->chrom_file_order_marker_idx[(*chrom_fo_idx_ptr) + 1];
  while (marker_uidx >= (*chrom_end_ptr)) {
    *chrom_end_ptr = chrom_info_ptr->chrom_file_order_marker_idx[(++(*chrom_fo_idx_ptr)) + 1];
  }
  chrom_idx = chrom_info_ptr->chrom_file_order[*chrom_fo_idx_ptr];
  *is_x_ptr = set_hh_missing && (chrom_idx == species_x_code[species]);
  *is_haploid_ptr = set_hh_missing && ((species_haploid_mask[species] >> chrom_idx) & 1LLU);
  if (is_all_nonmale) {
    *is_haploid_ptr = (*is_haploid_ptr) && (!(*is_x_ptr));
  }
}

// WDIST's natural sort uses the following logic:
// - All alphabetic characters act as if they are capitalized, except for
// tiebreaking purposes (where ASCII is used).
// - Numbers are compared by magnitude, with the exception of...
// - Numbers with leading zero(es).  If you're putting extraneous zeroes in
// front of IDs, we assume they're there to force particular items to be sorted
// earlier, rather than just appearing at random.  So, unlike many natural sort
// implementations, we sort 00200 < 021 < 20: all numbers with n leading zeroes
// are sorted before all numbers with (n-1) leading zeroes; magnitude only
// applies if the leading zero counts match.  This handles e.g. subbasement
// room numbering properly.
//
// This won't always do what you want if your IDs have variable-length decimals
// in them (e.g. it yields 0.99 < 0.101); if you don't want to fall back on
// ASCII sort, enforce a fixed number of digits after the decimal point.  Also
// note that ASCII sort is outright better for e.g. numbers represented in
// hexadecimal or base 36.  In principle, it's possible to reliably autodetect
// some of these cases (especially hexadecimal numbers beginning with "0x"),
// but that'll never be perfect so we just let the user toggle the sort method.
int32_t strcmp_natural_scan_forward(const char* s1, const char* s2) {
  // assumes s1 and s2 currently point to the middle of a mismatching number,
  // where s1 < s2.
  char c1;
  char c2;
  do {
    c1 = *(++s1);
    c2 = *(++s2);
    if (is_not_digit(c1)) {
      return -1;
    }
  } while (is_digit(c2));
  return 1;
}

// We have the following major states:
//   0 (initial): strings perfectly match so far, last char (if any) is
//                nonnumeric.
//   1: strings perfectly match so far, last char is numeric.
//   2: strings match except for capitalization, last char is nonnumeric.
//   3: strings match except for capitalization, last char is numeric.
// strcmp_natural_tiebroken() expresses the logic for states 2 and 3, while
// strcmp_natural_uncasted() handles states 0 and 1.
int32_t strcmp_natural_tiebroken(const char* s1, const char* s2) {
  // assumes ties should be broken in favor of s2.
  char c1 = *(++s1);
  char c2 = *(++s2);
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 2
  strcmp_natural_tiebroken_state_2:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	c1 -= 32;
      }
      if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
      }
      if (c1 < c2) {
	return -1;
      } else if (c1 > c2) {
	return 1;
      }
    } else if (!c1) {
      return -1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  }
  if (is_not_nzdigit(c1) || is_not_nzdigit(c2)) {
    return (c1 < c2)? -1 : 1;
  }
  do {
    // state 3
    if (c1 != c2) {
      if (is_digit(c2)) {
	if (c1 < c2) {
	  return strcmp_natural_scan_forward(s1, s2);
	} else {
	  return -strcmp_natural_scan_forward(s2, s1);
	}
      }
      return 1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  } while (is_digit(c1));
  if (is_digit(c2)) {
    return -1;
  }
  // skip the while (is_not_digit...) check
  goto strcmp_natural_tiebroken_state_2;
}

static inline int32_t strcmp_natural_uncasted(const char* s1, const char* s2) {
  char c1 = *s1;
  char c2 = *s2;
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 0
  strcmp_natural_uncasted_state_0:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	if (c2 + 32 == c1) {
	  return -strcmp_natural_tiebroken(s2, s1);
	} else if ((c2 < 'a') || (c2 > 'z')) {
	  c1 -= 32;
	}
      } else if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
	if (c1 == c2) {
	  return strcmp_natural_tiebroken(s1, s2);
	}
      }
      return (c1 < c2)? -1 : 1;
    } else if (!c1) {
      return 0;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  }
  if (is_not_nzdigit(c1) || is_not_nzdigit(c2)) {
    return (c1 < c2)? -1 : 1;
  }
  do {
    // state 1
    if (c1 != c2) {
      if (is_digit(c2)) {
	if (c1 < c2) {
	  return strcmp_natural_scan_forward(s1, s2);
	} else {
	  return -strcmp_natural_scan_forward(s2, s1);
	}
      }
      return 1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  } while (is_digit(c1));
  if (is_digit(c2)) {
    return -1;
  }
  goto strcmp_natural_uncasted_state_0;
}

int32_t strcmp_natural(const void* s1, const void* s2) {
  return strcmp_natural_uncasted((char*)s1, (char*)s2);
}

int32_t strcmp_deref(const void* s1, const void* s2) {
  return strcmp(*(char**)s1, *(char**)s2);
}

int32_t strcmp_natural_deref(const void* s1, const void* s2) {
  return strcmp_natural_uncasted(*(char**)s1, *(char**)s2);
}

int32_t is_missing(char* bufptr, int32_t missing_pheno, int32_t missing_pheno_len, int32_t affection_01) {
  if ((atoi(bufptr) == missing_pheno) && is_space_or_eoln(bufptr[missing_pheno_len])) {
    return 1;
  } else if ((!affection_01) && (*bufptr == '0') && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

int32_t eval_affection(char* bufptr, int32_t missing_pheno, int32_t missing_pheno_len, int32_t affection_01) {
  if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
    return 1;
  } else if (((*bufptr == '0') || (*bufptr == '1') || ((*bufptr == '2') && (!affection_01))) && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

int32_t triangle_divide(int64_t cur_prod, int32_t modif) {
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod, and neither term in the product is negative.  (Note the
  // lack of a divide by two; cur_prod should also be double its "true" value
  // as a result.)
  int64_t vv;
  if (cur_prod == 0) {
    if (modif < 0) {
      return -modif;
    } else {
      return 0;
    }
  }
  vv = (int64_t)sqrt((double)cur_prod);
  while ((vv - 1) * (vv + modif - 1) >= cur_prod) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod) {
    vv++;
  }
  return vv;
}

void parallel_bounds(int32_t ct, int32_t start, int32_t parallel_idx, int32_t parallel_tot, int32_t* bound_start_ptr, int32_t* bound_end_ptr) {
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = (int64_t)ct * (ct + modif);
  *bound_start_ptr = triangle_divide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = triangle_divide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// set align to 1 for no alignment
void triangle_fill(uint32_t* target_arr, int32_t ct, int32_t pieces, int32_t parallel_idx, int32_t parallel_tot, int32_t start, int32_t align) {
  int64_t ct_tr;
  int64_t cur_prod;
  int32_t modif = 1 - start * 2;
  int32_t cur_piece = 1;
  int32_t lbound;
  int32_t ubound;
  int32_t ii;
  int32_t align_m1;
  parallel_bounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = (int64_t)lbound * (lbound + modif);
  ct_tr = ((int64_t)ubound * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide(cur_prod, modif);
    ii = (lbound - start) & align_m1;
    if ((ii) && (ii != align_m1)) {
      lbound = start + ((lbound - start) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (lbound > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

int32_t write_ids(char* outname, uint32_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len) {
  FILE* outfile;
  uint32_t uii;
  if (fopen_checked(&outfile, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
    if (!is_set(indiv_exclude, uii) && (fprintf(outfile, "%s\n", &(person_ids[uii * max_person_id_len])) < 0)) {
      return RET_WRITE_FAIL;
    }
  }
  if (fclose(outfile)) {
    return RET_WRITE_FAIL;
  }
  return 0;
}

int32_t distance_d_write_ids(char* outname, char* outname_end, int32_t dist_calc_type, uint32_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len) {
  int32_t retval;
  if (dist_calc_type & DISTANCE_ALCT) {
    strcpy(outname_end, ".dist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
  }
  return 0;
}

int32_t relationship_req(uint64_t calculation_type) {
  return (calculation_type & (CALC_RELATIONSHIP | CALC_UNRELATED_HERITABILITY | CALC_REL_CUTOFF | CALC_REGRESS_REL));
}

int32_t distance_req(uint64_t calculation_type) {
  return ((calculation_type & CALC_DISTANCE) || ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))));
}

int32_t double_cmp(const void* aa, const void* bb) {
  double cc = *((const double*)aa) - *((const double*)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

int32_t double_cmp_deref(const void* aa, const void* bb) {
  double cc = **((const double**)aa) - **((const double**)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

#ifndef __cplusplus
int32_t llcmp(const void* aa, const void* bb) {
  int64_t diff = *((const int64_t*)aa) - *((const int64_t*)bb);
  if (diff > 0) {
    return 1;
  } else if (diff < 0) {
    return -1;
  } else {
    return 0;
  }
}
#endif

// alas, qsort_r not available on some Linux distributions

// Normally use qsort_ext(), but this version is necessary before wkspace has
// been allocated.
void qsort_ext2(char* main_arr, int32_t arr_length, int32_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int32_t secondary_item_len, char* proxy_arr, int32_t proxy_len) {
  int32_t ii;
  for (ii = 0; ii < arr_length; ii++) {
    *(char**)(&(proxy_arr[ii * proxy_len])) = &(main_arr[ii * item_length]);
    memcpy(&(proxy_arr[ii * proxy_len + sizeof(void*)]), &(secondary_arr[ii * secondary_item_len]), secondary_item_len);
  }
  qsort(proxy_arr, arr_length, proxy_len, comparator_deref);
  for (ii = 0; ii < arr_length; ii++) {
    memcpy(&(secondary_arr[ii * secondary_item_len]), &(proxy_arr[ii * proxy_len + sizeof(void*)]), secondary_item_len);
    memcpy(&(proxy_arr[ii * proxy_len]), *(char**)(&(proxy_arr[ii * proxy_len])), item_length);
  }
  for (ii = 0; ii < arr_length; ii++) {
    memcpy(&(main_arr[ii * item_length]), &(proxy_arr[ii * proxy_len]), item_length);
  }
}

// This actually tends to be faster than just sorting an array of indices,
// because of memory locality issues.
int32_t qsort_ext(char* main_arr, int32_t arr_length, int32_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int32_t secondary_item_len) {
  // main_arr = packed array of equal-length items to sort
  // arr_length = number of items
  // item_length = byte count of each main_arr item
  // comparator_deref = returns positive if *first > *second, 0 if equal,
  //                    negative if *first < *second.  Note the extra
  //                    dereference.
  // secondary_arr = packed array of fixed-length records associated with the
  //                 main_arr items, to be resorted in the same way.  (e.g.
  //                 if one is building an index, this could start as a sorted
  //                 0..(n-1) sequence of integers; then, post-sort, this would
  //                 be a lookup table for the original position of each
  //                 main_arr item.)
  // secondary_item_len = byte count of each secondary_arr item
  int32_t proxy_len = secondary_item_len + sizeof(void*);
  unsigned char* wkspace_mark = wkspace_base;
  char* proxy_arr;
  if (!arr_length) {
    return 0;
  }
  if (proxy_len < item_length) {
    proxy_len = item_length;
  }
  if (wkspace_alloc_c_checked(&proxy_arr, arr_length * proxy_len)) {
    return -1;
  }
  qsort_ext2(main_arr, arr_length, item_length, comparator_deref, secondary_arr, secondary_item_len, proxy_arr, proxy_len);
  wkspace_reset(wkspace_mark);
  return 0;
}

uint32_t doublearr_greater_than(double* sorted_dbl_arr, uint32_t arr_length, double dxx) {
  // assumes arr_length is nonzero.
  // dxx guaranteed to be larger than sorted_dbl_arr[min_idx - 1] if it exists,
  // but NOT necessarily sorted_dbl_arr[min_idx].
  int32_t min_idx = 0;
  // similarly, dxx guaranteed to be no greater than sorted_dbl_arr[max_idx +
  // 1] if it exists, but not necessarily sorted_dbl_arr[max_idx].
  // Signed integer since it could become -1.
  int32_t max_idx = arr_length - 1;
  uint32_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uint32_t)min_idx) + ((uint32_t)max_idx)) / 2;
    if (dxx > sorted_dbl_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (dxx > sorted_dbl_arr[((uint32_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

int32_t bsearch_str(char* id_buf, char* lptr, intptr_t max_id_len, int32_t min_idx, int32_t max_idx) {
  int32_t mid_idx;
  int32_t ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp(id_buf, &(lptr[mid_idx * max_id_len]));
  if (ii) {
    if (ii < 0) {
      return bsearch_str(id_buf, lptr, max_id_len, min_idx, mid_idx - 1);
    } else {
      return bsearch_str(id_buf, lptr, max_id_len, mid_idx + 1, max_idx);
    }
  } else {
    return mid_idx;
  }
}

int32_t bsearch_str_natural(char* id_buf, char* lptr, intptr_t max_id_len, int32_t min_idx, int32_t max_idx) {
  int32_t mid_idx;
  int32_t ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp_natural(id_buf, &(lptr[mid_idx * max_id_len]));
  if (ii) {
    if (ii < 0) {
      return bsearch_str_natural(id_buf, lptr, max_id_len, min_idx, mid_idx - 1);
    } else {
      return bsearch_str_natural(id_buf, lptr, max_id_len, mid_idx + 1, max_idx);
    }
  } else {
    return mid_idx;
  }
}

void fill_idbuf_fam_indiv(char* idbuf, char* fam_indiv, char fillchar) {
  char* iend_ptr = item_endnn(fam_indiv);
  uint32_t slen = (iend_ptr - fam_indiv);
  uint32_t slen2;
  memcpy(idbuf, fam_indiv, slen);
  idbuf[slen] = fillchar;
  fam_indiv = skip_initial_spaces(iend_ptr);
  iend_ptr = item_endnn(fam_indiv);
  slen2 = (iend_ptr - fam_indiv);
  memcpyx(&(idbuf[slen + 1]), fam_indiv, slen2, '\0');
}

int32_t bsearch_fam_indiv(char* id_buf, char* lptr, intptr_t max_id_len, int32_t filter_line_ct, char* fam_id, char* indiv_id) {
  // id_buf = workspace
  // lptr = packed, sorted list of ID strings to search over
  // fam_id and indiv_id are considered terminated by any space/eoln character
  int32_t ii;
  int32_t jj;
  if (!filter_line_ct) {
    return -1;
  }
  ii = strlen_se(fam_id);
  jj = strlen_se(indiv_id);
  if (ii + jj + 2 > max_id_len) {
    return -1;
  }
  memcpyx(memcpyax(id_buf, fam_id, ii, '\t'), indiv_id, jj, '\0');
  return bsearch_str(id_buf, lptr, max_id_len, 0, filter_line_ct - 1);
}

int32_t distance_open(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, char* outname, char* outname_end, const char* varsuffix, const char* mode, int32_t dist_calc_type, int32_t parallel_idx, int32_t parallel_tot) {
  if (dist_calc_type & DISTANCE_ALCT) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".dist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".dist%s", varsuffix);
    }
    strcpy(tbuf, outname_end);
    if (fopen_checked(outfile_ptr, outname, mode)) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mibs%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mibs%s", varsuffix);
    }
    strcpy(&(tbuf[MAX_POST_EXT]), outname_end);
    if (fopen_checked(outfile2_ptr, outname, mode)) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mdist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mdist%s", varsuffix);
    }
    strcpy(&(tbuf[MAX_POST_EXT * 2]), outname_end);
    if (fopen_checked(outfile3_ptr, outname, mode)) {
      return 1;
    }
  }
  return 0;
}

void distance_print_done(int32_t format_code, char* outname, char* outname_end) {
  putchar('\r');
  if (!format_code) {
    strcpy(outname_end, tbuf);
    sprintf(logbuf, "Distances (allele counts) written to %s.\n", outname);
  } else if (format_code == 1) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT]));
    sprintf(logbuf, "IBS matrix written to %s.\n", outname);
  } else if (format_code == 2) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT * 2]));
    sprintf(logbuf, "Distances (proportions) written to %s.\n", outname);
  }
  logprintb();
}

void bitfield_andnot(uintptr_t* vv, uintptr_t* exclude_vec, uintptr_t ct) {
  // vv := vv ANDNOT exclude_vec
  // on 64-bit systems, assumes vv and exclude_vec are 16-byte aligned
  // note that this is the reverse of the _mm_andnot() operand order
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)vv;
  __m128i* ev128 = (__m128i*)exclude_vec;
  __m128i* vv128_end = &(vv128[ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_andnot_si128(*ev128++, *vv128);
    vv128++;
  }
  if (ct & 1) {
    ct--;
    vv[ct] &= ~(exclude_vec[ct]);
  }
#else
  uintptr_t* vec_end = &(vv[ct]);
  do {
    *vv++ &= ~(*exclude_vec++);
  } while (vv < vec_end);
#endif
}

#ifdef __LP64__
// Basic SSE2 implementation of Lauradoux/Walisch popcount.
static inline uintptr_t popcount_vecs(__m128i* vptr, uintptr_t ct) {
  // popcounts vptr[0..(ct-1)].  Assumes ct is a multiple of 3 (0 ok).
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  __m128i* vend;
  __m128i count1, count2, half1, half2;
  __uni16 acc;

  while (ct >= 30) {
    ct -= 30;
    acc.vi = _mm_setzero_si128();
    vend = &(vptr[30]);
  popcount_vecs_main_loop:
    do {
      count1 = *vptr++;
      count2 = *vptr++;
      half1 = *vptr++;
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1
      // count2 store a partial bitcount covering themselves AND another bit
      // from elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    acc.vi = _mm_setzero_si128();
    vend = &(vptr[ct]);
    ct = 0;
    goto popcount_vecs_main_loop;
  }
  return tot;
}

static inline uintptr_t popcount_vecs_exclude(__m128i* vptr, __m128i* exclude_ptr, uintptr_t ct) {
  // popcounts vptr ANDNOT exclude_ptr[0..(ct-1)].  ct is a multiple of 3.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  __m128i* vend;
  __m128i count1, count2, half1, half2;
  __uni16 acc;

  while (ct >= 30) {
    ct -= 30;
    acc.vi = _mm_setzero_si128();
    vend = &(vptr[30]);
  popcount_vecs_exclude_main_loop:
    do {
      // nots the FIRST value
      count1 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
      count2 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
      half1 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    acc.vi = _mm_setzero_si128();
    vend = &(vptr[ct]);
    ct = 0;
    goto popcount_vecs_exclude_main_loop;
  }
  return tot;
}
#endif

uintptr_t popcount_longs(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  // given an aligned long array lptr[], this efficiently popcounts
  // lptr[start_idx..(end_idx - 1)].
  uintptr_t tot = 0;
  uintptr_t* lptr_end = &(lptr[end_idx]);
#ifdef __LP64__
  uintptr_t six_ct;
  __m128i* vptr;
  if (start_idx == end_idx) {
    return 0;
  }
  if (start_idx & 1) {
    tot = popcount_long(lptr[start_idx++]);
  }
  vptr = (__m128i*)(&(lptr[start_idx]));
  six_ct = (end_idx - start_idx) / 6;
  tot += popcount_vecs(vptr, six_ct * 3);
  lptr = &(lptr[start_idx + six_ct * 6]);
#else
  // The humble 16-bit lookup table actually beats
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
  // on my development machine by a hair.
  // However, if we take the hint from Lauradoux/Walisch and postpone the
  // multiply and right shift, this is no longer true.  Ah well.
  uintptr_t* lptr_six_end;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr = &(lptr[start_idx]);
  lptr_six_end = &(lptr[6 * ((end_idx - start_idx) / 6)]);
  while (lptr < lptr_six_end) {
    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount_long(*lptr++);
  }
  return tot;
}

uintptr_t popcount_chars(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  // given a CHAR array c[] starting at the aligned address of lptr, this
  // efficiently popcounts c[start_idx..(end_idx - 1)].
  uint32_t extra_ct = (start_idx % (BITCT / 8));
  uintptr_t si8l = start_idx / (BITCT / 8);
  uintptr_t ei8l = end_idx / (BITCT / 8);
  uintptr_t tot = 0;
  if (si8l == ei8l) {
    return popcount_long(lptr[si8l] & ((ONELU << ((end_idx % (BITCT / 8)) * 8)) - (ONELU << (extra_ct * 8))));
  } else {
    if (extra_ct) {
      tot = popcount_long(lptr[si8l++] >> (extra_ct * 8));
    }
    tot += popcount_longs(lptr, si8l, ei8l);
    extra_ct = end_idx % (BITCT / 8);
    if (extra_ct) {
      tot += popcount_long(lptr[ei8l] & ((ONELU << (extra_ct * 8)) - ONELU));
    }
    return tot;
  }
}

uintptr_t popcount_longs_exclude(uintptr_t* lptr, uintptr_t* exclude_arr, uintptr_t end_idx) {
  // popcounts lptr ANDNOT exclude_arr[0..(end_idx-1)].
  // N.B. on 64-bit systems, assumes lptr and exclude_arr are 16-byte aligned.
  uintptr_t tot = 0;
  uintptr_t* lptr_end = &(lptr[end_idx]);
#ifdef __LP64__
  uintptr_t six_ct = end_idx / 6;
  tot += popcount_vecs_exclude((__m128i*)lptr, (__m128i*)exclude_arr, six_ct * 3);
  lptr = &(lptr[six_ct * 6]);
  exclude_arr = &(exclude_arr[six_ct * 6]);
#else
  uintptr_t* lptr_six_end;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr_six_end = &(lptr[end_idx - (end_idx % 6)]);
  while (lptr < lptr_six_end) {
    loader = (*lptr++) & (~(*exclude_arr++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = (*lptr++) & (~(*exclude_arr++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount_long((*lptr++) & (~(*exclude_arr++)));
  }
  return tot;
}

#ifdef __LP64__
void count_2freq_dbl_60v(__m128i* vptr, __m128i* vend, __m128i* mask1vp, __m128i* mask2vp, uint32_t* ct1abp, uint32_t* ct1cp, uint32_t* ct2abp, uint32_t* ct2cp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct1_ab;
  __m128i to_ct_abtmp;
  __m128i to_ct1_c;
  __m128i to_ct2_ab;
  __m128i to_ct2_c;
  __uni16 acc1_ab;
  __uni16 acc1_c;
  __uni16 acc2_ab;
  __uni16 acc2_c;

  acc1_ab.vi = _mm_setzero_si128();
  acc1_c.vi = _mm_setzero_si128();
  acc2_ab.vi = _mm_setzero_si128();
  acc2_c.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = *mask1vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct1_ab = _mm_add_epi64(loader3, loader2);
    to_ct1_c = _mm_andnot_si128(loader3, loader2);
    loader2 = *mask2vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct2_ab = _mm_add_epi64(loader3, loader2);
    to_ct2_c = _mm_andnot_si128(loader3, loader2);
    to_ct1_ab = _mm_add_epi64(_mm_and_si128(to_ct1_ab, m2), _mm_and_si128(_mm_srli_epi64(to_ct1_ab, 2), m2));
    to_ct2_ab = _mm_add_epi64(_mm_and_si128(to_ct2_ab, m2), _mm_and_si128(_mm_srli_epi64(to_ct2_ab, 2), m2));

    loader = *vptr++;
    loader2 = *mask1vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct1_c = _mm_add_epi64(to_ct1_c, _mm_andnot_si128(loader3, loader2));
    to_ct1_ab = _mm_add_epi64(to_ct1_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));
    loader2 = *mask2vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct2_c = _mm_add_epi64(to_ct2_c, _mm_andnot_si128(loader3, loader2));
    to_ct2_ab = _mm_add_epi64(to_ct2_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));

    loader = *vptr++;
    loader2 = *mask1vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct1_c = _mm_add_epi64(to_ct1_c, _mm_andnot_si128(loader3, loader2));
    to_ct1_ab = _mm_add_epi64(to_ct1_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));
    loader2 = *mask2vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct2_c = _mm_add_epi64(to_ct2_c, _mm_andnot_si128(loader3, loader2));
    to_ct2_ab = _mm_add_epi64(to_ct2_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));

    to_ct1_c = _mm_add_epi64(_mm_and_si128(to_ct1_c, m2), _mm_and_si128(_mm_srli_epi64(to_ct1_c, 2), m2));
    to_ct2_c = _mm_add_epi64(_mm_and_si128(to_ct2_c, m2), _mm_and_si128(_mm_srli_epi64(to_ct2_c, 2), m2));

    acc1_ab.vi = _mm_add_epi64(acc1_ab.vi, _mm_add_epi64(_mm_and_si128(to_ct1_ab, m4), _mm_and_si128(_mm_srli_epi64(to_ct1_ab, 4), m4)));
    acc1_c.vi = _mm_add_epi64(acc1_c.vi, _mm_add_epi64(_mm_and_si128(to_ct1_c, m4), _mm_and_si128(_mm_srli_epi64(to_ct1_c, 4), m4)));
    acc2_ab.vi = _mm_add_epi64(acc2_ab.vi, _mm_add_epi64(_mm_and_si128(to_ct2_ab, m4), _mm_and_si128(_mm_srli_epi64(to_ct2_ab, 4), m4)));
    acc2_c.vi = _mm_add_epi64(acc2_c.vi, _mm_add_epi64(_mm_and_si128(to_ct2_c, m4), _mm_and_si128(_mm_srli_epi64(to_ct2_c, 4), m4)));
  } while (vptr < vend);
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  acc1_ab.vi = _mm_add_epi64(_mm_and_si128(acc1_ab.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1_ab.vi, 8), m8));
  acc1_c.vi = _mm_and_si128(_mm_add_epi64(acc1_c.vi, _mm_srli_epi64(acc1_c.vi, 8)), m8);
  acc2_ab.vi = _mm_add_epi64(_mm_and_si128(acc2_ab.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2_ab.vi, 8), m8));
  acc2_c.vi = _mm_and_si128(_mm_add_epi64(acc2_c.vi, _mm_srli_epi64(acc2_c.vi, 8)), m8);
  *ct1abp += ((acc1_ab.u8[0] + acc1_ab.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct1cp += ((acc1_c.u8[0] + acc1_c.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct2abp += ((acc2_ab.u8[0] + acc2_ab.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct2cp += ((acc2_c.u8[0] + acc2_c.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_3freq_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* even_ctp, uint32_t* odd_ctp, uint32_t* homset_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i even1;
  __m128i odd1;
  __m128i homset1;
  __m128i even2;
  __m128i odd2;
  __m128i homset2;
  __uni16 acc_even;
  __uni16 acc_odd;
  __uni16 acc_homset;

  acc_even.vi = _mm_setzero_si128();
  acc_odd.vi = _mm_setzero_si128();
  acc_homset.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = *maskvp++;
    odd1 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even1 = _mm_and_si128(loader2, loader);
    homset1 = _mm_and_si128(odd1, loader);
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even1 = _mm_add_epi64(even1, _mm_and_si128(loader2, loader));
    odd1 = _mm_add_epi64(odd1, loader3);
    homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even1 = _mm_add_epi64(even1, _mm_and_si128(loader2, loader));
    odd1 = _mm_add_epi64(odd1, loader3);
    homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));

    even1 = _mm_add_epi64(_mm_and_si128(even1, m2), _mm_and_si128(_mm_srli_epi64(even1, 2), m2));
    odd1 = _mm_add_epi64(_mm_and_si128(odd1, m2), _mm_and_si128(_mm_srli_epi64(odd1, 2), m2));
    homset1 = _mm_add_epi64(_mm_and_si128(homset1, m2), _mm_and_si128(_mm_srli_epi64(homset1, 2), m2));

    loader = *vptr++;
    loader2 = *maskvp++;
    odd2 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even2 = _mm_and_si128(loader2, loader);
    homset2 = _mm_and_si128(odd2, loader);
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even2 = _mm_add_epi64(even2, _mm_and_si128(loader2, loader));
    odd2 = _mm_add_epi64(odd2, loader3);
    homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even2 = _mm_add_epi64(even2, _mm_and_si128(loader2, loader));
    odd2 = _mm_add_epi64(odd2, loader3);
    homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));

    even1 = _mm_add_epi64(even1, _mm_add_epi64(_mm_and_si128(even2, m2), _mm_and_si128(_mm_srli_epi64(even2, 2), m2)));
    odd1 = _mm_add_epi64(odd1, _mm_add_epi64(_mm_and_si128(odd2, m2), _mm_and_si128(_mm_srli_epi64(odd2, 2), m2)));
    homset1 = _mm_add_epi64(homset1, _mm_add_epi64(_mm_and_si128(homset2, m2), _mm_and_si128(_mm_srli_epi64(homset2, 2), m2)));

    acc_even.vi = _mm_add_epi64(acc_even.vi, _mm_add_epi64(_mm_and_si128(even1, m4), _mm_and_si128(_mm_srli_epi64(even1, 4), m4)));
    acc_odd.vi = _mm_add_epi64(acc_odd.vi, _mm_add_epi64(_mm_and_si128(odd1, m4), _mm_and_si128(_mm_srli_epi64(odd1, 4), m4)));
    acc_homset.vi = _mm_add_epi64(acc_homset.vi, _mm_add_epi64(_mm_and_si128(homset1, m4), _mm_and_si128(_mm_srli_epi64(homset1, 4), m4)));
  } while (vptr < vend);
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  acc_even.vi = _mm_add_epi64(_mm_and_si128(acc_even.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_even.vi, 8), m8));
  acc_odd.vi = _mm_add_epi64(_mm_and_si128(acc_odd.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_odd.vi, 8), m8));
  acc_homset.vi = _mm_add_epi64(_mm_and_si128(acc_homset.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_homset.vi, 8), m8));
  *even_ctp += ((acc_even.u8[0] + acc_even.u8[1]) * 0x1000100010001LLU) >> 48;
  *odd_ctp += ((acc_odd.u8[0] + acc_odd.u8[1]) * 0x1000100010001LLU) >> 48;
  *homset_ctp += ((acc_homset.u8[0] + acc_homset.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_3freq_xx_120v(__m128i* vptr, __m128i* vend, __m128i* include_vec, __m128i* male_vec, uint32_t* missing_ctp, uint32_t* odd_ctp, uint32_t* homset_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i loader4;
  __m128i missing1;
  __m128i odd1;
  __m128i homset1;
  __m128i missing2;
  __m128i odd2;
  __m128i homset2;
  __uni16 accm;
  __uni16 acco;
  __uni16 acch;

  accm.vi = _mm_setzero_si128();
  acco.vi = _mm_setzero_si128();
  acch.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missing1 = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader)));
    odd1 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    homset1 = _mm_and_si128(odd1, loader);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missing1 = _mm_add_epi64(missing1, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader))));
    loader3 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    odd1 = _mm_add_epi64(odd1, loader3);
    homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missing1 = _mm_add_epi64(missing1, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader))));
    loader3 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    odd1 = _mm_add_epi64(odd1, loader3);
    homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));

    missing1 = _mm_add_epi64(_mm_and_si128(missing1, m2), _mm_and_si128(_mm_srli_epi64(missing1, 2), m2));
    odd1 = _mm_add_epi64(_mm_and_si128(odd1, m2), _mm_and_si128(_mm_srli_epi64(odd1, 2), m2));
    homset1 = _mm_add_epi64(_mm_and_si128(homset1, m2), _mm_and_si128(_mm_srli_epi64(homset1, 2), m2));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missing2 = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader)));
    odd2 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    homset2 = _mm_and_si128(odd2, loader);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missing2 = _mm_add_epi64(missing2, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader))));
    loader3 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    odd2 = _mm_add_epi64(odd2, loader3);
    homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missing2 = _mm_add_epi64(missing2, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader))));
    loader3 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    odd2 = _mm_add_epi64(odd2, loader3);
    homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));

    missing1 = _mm_add_epi64(missing1, _mm_add_epi64(_mm_and_si128(missing2, m2), _mm_and_si128(_mm_srli_epi64(missing2, 2), m2)));
    odd1 = _mm_add_epi64(odd1, _mm_add_epi64(_mm_and_si128(odd2, m2), _mm_and_si128(_mm_srli_epi64(odd2, 2), m2)));
    homset1 = _mm_add_epi64(homset1, _mm_add_epi64(_mm_and_si128(homset2, m2), _mm_and_si128(_mm_srli_epi64(homset2, 2), m2)));

    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(missing1, m4), _mm_and_si128(_mm_srli_epi64(missing1, 4), m4)));
    acco.vi = _mm_add_epi64(acco.vi, _mm_add_epi64(_mm_and_si128(odd1, m4), _mm_and_si128(_mm_srli_epi64(odd1, 4), m4)));
    acch.vi = _mm_add_epi64(acch.vi, _mm_add_epi64(_mm_and_si128(homset1, m4), _mm_and_si128(_mm_srli_epi64(homset1, 4), m4)));
  } while (vptr < vend);
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  acco.vi = _mm_add_epi64(_mm_and_si128(acco.vi, m8), _mm_and_si128(_mm_srli_epi64(acco.vi, 8), m8));
  acch.vi = _mm_add_epi64(_mm_and_si128(acch.vi, m8), _mm_and_si128(_mm_srli_epi64(acch.vi, 8), m8));
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
  *odd_ctp += ((acco.u8[0] + acco.u8[1]) * 0x1000100010001LLU) >> 48;
  *homset_ctp += ((acch.u8[0] + acch.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
void count_2freq_dbl_6(uintptr_t* lptr, uintptr_t* mask1p, uintptr_t* mask2p, uint32_t* ct1abp, uint32_t* ct1cp, uint32_t* ct2abp, uint32_t* ct2cp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = *mask1p++;
  uintptr_t loader3 = (loader >> 1) & loader2;
  uintptr_t to_ct1_ab;
  uintptr_t to_ct1_c;
  uintptr_t to_ct2_ab;
  uintptr_t to_ct2_c;
  uintptr_t to_ct_abtmp;
  uintptr_t partial1_ab;
  uintptr_t partial1_c;
  uintptr_t partial2_ab;
  uintptr_t partial2_c;
  loader2 &= loader;
  to_ct1_ab = loader2 + loader3;
  to_ct1_c = loader2 & (~loader3);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct2_ab = loader2 + loader3;
  to_ct2_c = loader2 & (~loader3);

  to_ct1_ab = (to_ct1_ab & 0x33333333) + ((to_ct1_ab >> 2) & 0x33333333);
  to_ct2_ab = (to_ct2_ab & 0x33333333) + ((to_ct2_ab >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  partial1_ab = (to_ct1_ab & 0x0f0f0f0f) + ((to_ct1_ab >> 4) & 0x0f0f0f0f);
  partial1_c = (to_ct1_c & 0x33333333) + ((to_ct1_c >> 2) & 0x33333333);
  partial2_ab = (to_ct2_ab & 0x0f0f0f0f) + ((to_ct2_ab >> 4) & 0x0f0f0f0f);
  partial2_c = (to_ct2_c & 0x33333333) + ((to_ct2_c >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct1_ab = loader2 + loader3;
  to_ct1_c = loader2 & (~loader3);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct2_ab = loader2 + loader3;
  to_ct2_c = loader2 & (~loader3);

  to_ct1_ab = (to_ct1_ab & 0x33333333) + ((to_ct1_ab >> 2) & 0x33333333);
  to_ct2_ab = (to_ct2_ab & 0x33333333) + ((to_ct2_ab >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  partial1_ab += (to_ct1_ab & 0x0f0f0f0f) + ((to_ct1_ab >> 4) & 0x0f0f0f0f);
  partial1_c += (to_ct1_c & 0x33333333) + ((to_ct1_c >> 2) & 0x33333333);
  partial2_ab += (to_ct2_ab & 0x0f0f0f0f) + ((to_ct2_ab >> 4) & 0x0f0f0f0f);
  partial2_c += (to_ct2_c & 0x33333333) + ((to_ct2_c >> 2) & 0x33333333);

  partial1_c = (partial1_c & 0x0f0f0f0f) + ((partial1_c >> 4) & 0x0f0f0f0f);
  partial2_c = (partial2_c & 0x0f0f0f0f) + ((partial2_c >> 4) & 0x0f0f0f0f);

  *ct1abp += (partial1_ab * 0x01010101) >> 24;
  *ct1cp += (partial1_c * 0x01010101) >> 24;
  *ct2abp += (partial2_ab * 0x01010101) >> 24;
  *ct2cp += (partial2_c * 0x01010101) >> 24;
}

void count_3freq_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ctap, uint32_t* ctbp, uint32_t* ctcp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = *maskp++;
  uint32_t to_ct_a1 = loader & loader2;
  uint32_t to_ct_b1 = (loader >> 1) & loader2;
  uint32_t to_ct_c1 = loader & to_ct_b1;
  uintptr_t loader3;
  uint32_t to_ct_a2;
  uint32_t to_ct_b2;
  uint32_t to_ct_c2;
  uintptr_t partial_a;
  uintptr_t partial_b;
  uintptr_t partial_c;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;

  loader = *lptr++;
  loader2 = *maskp++;
  to_ct_a2 = loader & loader2;
  to_ct_b2 = (loader >> 1) & loader2;
  to_ct_c2 = loader & to_ct_b2;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;

  to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
  to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
  partial_a = (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
  to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
  to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
  partial_b = (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
  to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
  to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
  partial_c = (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = *maskp++;
  to_ct_a1 = loader & loader2;
  to_ct_b1 = (loader >> 1) & loader2;
  to_ct_c1 = loader & to_ct_b1;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;

  loader = *lptr++;
  loader2 = *maskp++;
  to_ct_a2 = loader & loader2;
  to_ct_b2 = (loader >> 1) & loader2;
  to_ct_c2 = loader & to_ct_b2;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;
  loader = *lptr;
  loader2 = *maskp;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;

  to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
  to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
  partial_a += (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
  to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
  to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
  partial_b += (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
  to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
  to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
  partial_c += (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

  *ctap += (partial_a * 0x01010101) >> 24;
  *ctbp += (partial_b * 0x01010101) >> 24;
  *ctcp += (partial_c * 0x01010101) >> 24;
}

void count_3freq_xx_12(uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* missing_ctp, uint32_t* odd_ctp, uint32_t* homset_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t loader4 = *male_vec++;
  uintptr_t missing1 = loader3 & (loader4 | (loader & (~loader2)));
  uintptr_t odd1 = loader3 & loader2 & (~loader4);
  uintptr_t homset1 = odd1 & loader;
  uintptr_t missing2;
  uintptr_t odd2;
  uintptr_t homset2;
  uintptr_t accm;
  uintptr_t acco;
  uintptr_t acch;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing1 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd1 += loader3;
  homset1 += loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing1 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd1 += loader3;
  homset1 += loader3 & loader;

  accm = (missing1 & 0x33333333) + ((missing1 >> 2) & 0x33333333);
  acco = (odd1 & 0x33333333) + ((odd1 >> 2) & 0x33333333);
  acch = (homset1 & 0x33333333) + ((homset1 >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing2 = loader3 & (loader4 | (loader & (~loader2)));
  odd2 = loader3 & loader2 & (~loader4);
  homset2 = odd2 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing2 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd2 += loader3;
  homset2 += loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing2 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd2 += loader3;
  homset2 += loader3 & loader;

  accm += (missing2 & 0x33333333) + ((missing2 >> 2) & 0x33333333);
  acco += (odd2 & 0x33333333) + ((odd2 >> 2) & 0x33333333);
  acch += (homset2 & 0x33333333) + ((homset2 >> 2) & 0x33333333);
  accm = (accm & 0x0f0f0f0f) + ((accm >> 4) & 0x0f0f0f0f);
  acco = (acco & 0x0f0f0f0f) + ((acco >> 4) & 0x0f0f0f0f);
  acch = (acch & 0x0f0f0f0f) + ((acch >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing1 = loader3 & (loader4 | (loader & (~loader2)));
  odd1 = loader3 & loader2 & (~loader4);
  homset1 = odd1 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing1 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd1 += loader3;
  homset1 += loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing1 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd1 += loader3;
  homset1 += loader3 & loader;

  missing1 = (missing1 & 0x33333333) + ((missing1 >> 2) & 0x33333333);
  odd1 = (odd1 & 0x33333333) + ((odd1 >> 2) & 0x33333333);
  homset1 = (homset1 & 0x33333333) + ((homset1 >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing2 = loader3 & (loader4 | (loader & (~loader2)));
  odd2 = loader3 & loader2 & (~loader4);
  homset2 = odd2 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing2 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd2 += loader3;
  homset2 += loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missing2 += loader3 & (loader4 | (loader & (~loader2)));
  loader3 &= loader2 & (~loader4);
  odd2 += loader3;
  homset2 += loader3 & loader;

  missing1 += (missing2 & 0x33333333) + ((missing2 >> 2) & 0x33333333);
  odd1 += (odd2 & 0x33333333) + ((odd2 >> 2) & 0x33333333);
  homset1 += (homset2 & 0x33333333) + ((homset2 >> 2) & 0x33333333);
  accm += (missing1 & 0x0f0f0f0f) + ((missing1 >> 4) & 0x0f0f0f0f);
  acco += (odd1 & 0x0f0f0f0f) + ((odd1 >> 4) & 0x0f0f0f0f);
  acch += (homset1 & 0x0f0f0f0f) + ((homset1 >> 4) & 0x0f0f0f0f);

  *missing_ctp += (accm * 0x01010101) >> 24;
  *odd_ctp += (acco * 0x01010101) >> 24;
  *homset_ctp += (acch * 0x01010101) >> 24;
}
#endif

#ifdef __LP64__
void count_set_freq_60v(__m128i* vptr, __m128i* vend, __m128i* include_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i odds;
  __m128i evens;
  __m128i missings;
  __uni16 acc;
  __uni16 accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    odds = _mm_and_si128(loader2, loader3);
    evens = _mm_and_si128(odds, loader);
    missings = _mm_and_si128(loader, _mm_andnot_si128(loader2, loader3));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    odds = _mm_add_epi64(odds, _mm_and_si128(loader2, loader3));
    loader3 = _mm_and_si128(loader, loader3);
    evens = _mm_add_epi64(evens, _mm_and_si128(loader2, loader3));
    missings = _mm_add_epi64(missings, _mm_andnot_si128(loader2, loader3));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    odds = _mm_add_epi64(odds, _mm_and_si128(loader2, loader3));
    loader3 = _mm_and_si128(loader, loader3);
    evens = _mm_add_epi64(evens, _mm_and_si128(loader2, loader3));
    missings = _mm_add_epi64(missings, _mm_andnot_si128(loader2, loader3));

    odds = _mm_add_epi64(_mm_and_si128(odds, m2), _mm_and_si128(_mm_srli_epi64(odds, 2), m2));
    missings = _mm_add_epi64(_mm_and_si128(missings, m2), _mm_and_si128(_mm_srli_epi64(missings, 2), m2));
    odds = _mm_add_epi64(odds, _mm_add_epi64(_mm_and_si128(evens, m2), _mm_and_si128(_mm_srli_epi64(evens, 2), m2)));

    // each 4-bit value here <= 6, so safe to add before m4 mask
    accm.vi = _mm_add_epi64(accm.vi, _mm_and_si128(_mm_add_epi64(missings, _mm_srli_epi64(missings, 4)), m4));

    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(odds, m4), _mm_and_si128(_mm_srli_epi64(odds, 4), m4)));
  } while (vptr < vend);
  // and each 8-bit value here <= 120
  accm.vi = _mm_and_si128(_mm_add_epi64(accm.vi, _mm_srli_epi64(accm.vi, 8)), m8);

  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_hap_120v(__m128i* vptr, __m128i* vend, __m128i* include_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __uni16 acc;
  __uni16 accm;
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i partial;
  __m128i partialm;
  __m128i partial2;
  __m128i partial2m;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_and_si128(loader3, _mm_and_si128(loader, loader2));
    partialm = _mm_and_si128(loader3, _mm_xor_si128(loader, loader2));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partial2 = _mm_add_epi64(_mm_and_si128(partial, m2), _mm_and_si128(_mm_srli_epi64(partial, 2), m2));
    partial2m = _mm_add_epi64(_mm_and_si128(partialm, m2), _mm_and_si128(_mm_srli_epi64(partialm, 2), m2));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_and_si128(loader3, _mm_and_si128(loader, loader2));
    partialm = _mm_and_si128(loader3, _mm_xor_si128(loader, loader2));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partial2 = _mm_add_epi64(partial2, _mm_add_epi64(_mm_and_si128(partial, m2), _mm_and_si128(_mm_srli_epi64(partial, 2), m2)));
    partial2m = _mm_add_epi64(partial2m, _mm_add_epi64(_mm_and_si128(partialm, m2), _mm_and_si128(_mm_srli_epi64(partialm, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(partial2, m4), _mm_and_si128(_mm_srli_epi64(partial2, 4), m4)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(partial2m, m4), _mm_and_si128(_mm_srli_epi64(partial2m, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_x_60v(__m128i* vptr, __m128i* vend, __m128i* include_vec, __m128i* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i loader4;
  __m128i set_odds;
  __m128i set_evens;
  __m128i missings_nm;
  __m128i missings_m;
  __m128i males;
  __uni16 acc;
  __uni16 accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = _mm_andnot_si128(*male_vec, loader3);
    set_evens = _mm_and_si128(loader, loader4); // subtract missings_nm later
    set_odds = _mm_and_si128(loader2, loader4);
    missings_nm = _mm_andnot_si128(loader2, set_evens);
    males = _mm_and_si128(loader3, *male_vec++);
    set_evens = _mm_or_si128(set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), males));
    missings_m = _mm_and_si128(_mm_xor_si128(loader, loader2), males);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = _mm_andnot_si128(*male_vec, loader3);
    set_odds = _mm_add_epi64(set_odds, _mm_and_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader, loader4);
    set_evens = _mm_add_epi64(set_evens, loader4);
    missings_nm = _mm_add_epi64(missings_nm, _mm_andnot_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader3, *male_vec++);
    set_evens = _mm_add_epi64(set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), loader4));
    missings_m = _mm_add_epi64(missings_m, _mm_and_si128(_mm_xor_si128(loader, loader2), loader4));
    males = _mm_add_epi64(males, loader4);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = _mm_andnot_si128(*male_vec, loader3);
    set_odds = _mm_add_epi64(set_odds, _mm_and_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader, loader4);
    set_evens = _mm_add_epi64(set_evens, loader4);
    missings_nm = _mm_add_epi64(missings_nm, _mm_andnot_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader3, *male_vec++);
    set_evens = _mm_add_epi64(set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), loader4));
    missings_m = _mm_add_epi64(missings_m, _mm_and_si128(_mm_xor_si128(loader, loader2), loader4));
    males = _mm_add_epi64(males, loader4);

    set_evens = _mm_sub_epi64(set_evens, missings_nm);
    missings_nm = _mm_slli_epi64(_mm_add_epi64(_mm_and_si128(missings_nm, m2), _mm_and_si128(_mm_srli_epi64(missings_nm, 2), m2)), 1);
    set_odds = _mm_add_epi64(_mm_and_si128(set_odds, m2), _mm_and_si128(_mm_srli_epi64(set_odds, 2), m2));
    missings_nm = _mm_add_epi64(missings_nm, _mm_add_epi64(_mm_and_si128(missings_m, m2), _mm_and_si128(_mm_srli_epi64(missings_m, 2), m2)));
    set_odds = _mm_add_epi64(set_odds, _mm_add_epi64(_mm_and_si128(set_evens, m2), _mm_and_si128(_mm_srli_epi64(set_evens, 2), m2)));
    missings_nm = _mm_add_epi64(missings_nm, _mm_add_epi64(_mm_and_si128(males, m2), _mm_and_si128(_mm_srli_epi64(males, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(set_odds, m4), _mm_and_si128(_mm_srli_epi64(set_odds, 4), m4)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(missings_nm, m4), _mm_and_si128(_mm_srli_epi64(missings_nm, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_xx_60v(__m128i* vptr, __m128i* vend, __m128i* include_vec, __m128i* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i loader4;
  __m128i set_odds;
  __m128i set_evens;
  __m128i missings;
  __uni16 acc;
  __uni16 accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missings = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader)));
    set_odds = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    set_evens = _mm_and_si128(set_odds, loader);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missings = _mm_add_epi64(missings, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader))));
    loader3 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    set_odds = _mm_add_epi64(set_odds, loader3);
    set_evens = _mm_add_epi64(set_evens, _mm_and_si128(loader3, loader));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    missings = _mm_add_epi64(missings, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_andnot_si128(loader2, loader))));
    loader3 = _mm_and_si128(loader3, _mm_andnot_si128(loader4, loader2));
    set_odds = _mm_add_epi64(set_odds, loader3);
    set_evens = _mm_add_epi64(set_evens, _mm_and_si128(loader3, loader));

    set_odds = _mm_add_epi64(_mm_and_si128(set_odds, m2), _mm_and_si128(_mm_srli_epi64(set_odds, 2), m2));
    missings = _mm_add_epi64(_mm_and_si128(missings, m2), _mm_and_si128(_mm_srli_epi64(missings, 2), m2));
    set_odds = _mm_add_epi64(set_odds, _mm_add_epi64(_mm_and_si128(set_evens, m2), _mm_and_si128(_mm_srli_epi64(set_evens, 2), m2)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_and_si128(_mm_add_epi64(missings, _mm_srli_epi64(missings, 4)), m4));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(set_odds, m4), _mm_and_si128(_mm_srli_epi64(set_odds, 4), m4)));
  } while (vptr < vend);
  accm.vi = _mm_and_si128(_mm_add_epi64(accm.vi, _mm_srli_epi64(accm.vi, 8)), m8);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_y_120v(__m128i* vptr, __m128i* vend, __m128i* include_vec, __m128i* nonmale_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i loader4;
  __m128i sets1;
  __m128i missings1;
  __m128i sets2;
  __m128i missings2;
  __uni16 acc;
  __uni16 accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets1 = _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2));
    missings1 = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2)));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets1 = _mm_add_epi64(sets1, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings1 = _mm_add_epi64(missings1, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets1 = _mm_add_epi64(sets1, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings1 = _mm_add_epi64(missings1, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));
    sets1 = _mm_add_epi64(_mm_and_si128(sets1, m2), _mm_and_si128(_mm_srli_epi64(sets1, 2), m2));
    missings1 = _mm_add_epi64(_mm_and_si128(missings1, m2), _mm_and_si128(_mm_srli_epi64(missings1, 2), m2));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets2 = _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2));
    missings2 = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2)));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets2 = _mm_add_epi64(sets2, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings2 = _mm_add_epi64(missings2, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets2 = _mm_add_epi64(sets2, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings2 = _mm_add_epi64(missings2, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));
    sets1 = _mm_add_epi64(sets1, _mm_add_epi64(_mm_and_si128(sets2, m2), _mm_and_si128(_mm_srli_epi64(sets2, 2), m2)));
    missings1 = _mm_add_epi64(missings1, _mm_add_epi64(_mm_and_si128(missings2, m2), _mm_and_si128(_mm_srli_epi64(missings2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sets1, m4), _mm_and_si128(_mm_srli_epi64(sets1, 4), m4)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(missings1, m4), _mm_and_si128(_mm_srli_epi64(missings1, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
void count_set_freq_6(uintptr_t* lptr, uintptr_t* include_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t odds = loader2 & loader3;
  uintptr_t evens = odds & loader;
  uintptr_t missings = (~loader2) & loader3 & loader;
  uintptr_t acc;
  uintptr_t accm;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  odds = (odds & 0x33333333) + ((odds >> 2) & 0x33333333);
  odds += (evens & 0x33333333) + ((evens >> 2) & 0x33333333);
  accm = (missings & 0x33333333) + ((missings >> 2) & 0x33333333);
  acc = (odds & 0x0f0f0f0f) + ((odds >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds = loader2 & loader3;
  evens = odds & loader;
  missings = (~loader2) & loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  odds = (odds & 0x33333333) + ((odds >> 2) & 0x33333333);
  accm += (missings & 0x33333333) + ((missings >> 2) & 0x33333333);
  odds += (evens & 0x33333333) + ((evens >> 2) & 0x33333333);
  accm = (accm & 0x0f0f0f0f) + ((accm >> 4) & 0x0f0f0f0f);
  acc += (odds & 0x0f0f0f0f) + ((odds >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

void count_set_freq_hap_12(uintptr_t* lptr, uintptr_t* include_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t partial = loader & loader2 & loader3;
  uintptr_t partialm = (loader ^ loader2) & loader3;
  uintptr_t partial2;
  uintptr_t partial2m;
  uintptr_t acc;
  uintptr_t accm;
  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 = (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m = (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial = loader & loader2 & loader3;
  partialm = (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 += (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m += (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);
  acc = (partial2 & 0x0f0f0f0f) + ((partial2 >> 4) & 0x0f0f0f0f);
  accm = (partial2m & 0x0f0f0f0f) + ((partial2m >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial = loader & loader2 & loader3;
  partialm = (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 = (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m = (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial = loader & loader2 & loader3;
  partialm = (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 += (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m += (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);
  acc += (partial2 & 0x0f0f0f0f) + ((partial2 >> 4) & 0x0f0f0f0f);
  accm += (partial2m & 0x0f0f0f0f) + ((partial2m >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

void count_set_freq_x_6(uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t loader4 = loader3 & (~(*male_vec));
  uintptr_t set_odds = loader2 & loader4;
  uintptr_t set_evens = loader & loader4;
  uintptr_t missings_nm = set_evens & (~loader2);
  uintptr_t missings_m;
  uintptr_t males;
  uintptr_t acc;
  uintptr_t accm;
  males = loader3 & (*male_vec++);
  set_evens |= loader & loader2 & males;
  missings_m = (loader ^ loader2) & males;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  set_evens -= missings_nm;
  set_odds = (set_odds & 0x33333333) + ((set_odds >> 2) & 0x33333333);
  set_odds += (set_evens & 0x33333333) + ((set_evens >> 2) & 0x33333333);
  missings_nm = ((missings_nm & 0x33333333) + ((missings_nm >> 2) & 0x33333333)) * 2;
  missings_nm += (missings_m & 0x33333333) + ((missings_m >> 2) & 0x33333333);
  missings_nm += (males & 0x33333333) + ((males >> 2) & 0x33333333);
  acc = (set_odds & 0x0f0f0f0f) + ((set_odds >> 4) & 0x0f0f0f0f);
  accm = (missings_nm & 0x0f0f0f0f) + ((missings_nm >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds = loader2 & loader4;
  set_evens = loader & loader4;
  missings_nm = set_evens & (~loader2);
  males = loader3 & (*male_vec++);
  set_evens |= loader & loader2 & males;
  missings_m = (loader ^ loader2) & males;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  set_evens -= missings_nm;
  set_odds = (set_odds & 0x33333333) + ((set_odds >> 2) & 0x33333333);
  set_odds += (set_evens & 0x33333333) + ((set_evens >> 2) & 0x33333333);
  missings_nm = ((missings_nm & 0x33333333) + ((missings_nm >> 2) & 0x33333333)) * 2;
  missings_nm += (missings_m & 0x33333333) + ((missings_m >> 2) & 0x33333333);
  missings_nm += (males & 0x33333333) + ((males >> 2) & 0x33333333);
  acc += (set_odds & 0x0f0f0f0f) + ((set_odds >> 4) & 0x0f0f0f0f);
  accm += (missings_nm & 0x0f0f0f0f) + ((missings_nm >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

void count_set_freq_xx_6(uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t loader4 = *male_vec++;
  uintptr_t missings = loader3 & ((loader & (~loader2)) | loader4);
  uintptr_t set_odds = loader3 & loader2 & (~loader4);
  uintptr_t set_evens = loader & set_odds;
  uintptr_t acc;
  uintptr_t accm;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missings += loader3 & ((loader & (~loader2)) | loader4);
  loader3 &= loader2 & (~loader4);
  set_odds += loader3;
  set_evens += loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missings += loader3 & ((loader & (~loader2)) | loader4);
  loader3 &= loader2 & (~loader4);
  set_odds += loader3;
  set_evens += loader3 & loader;

  set_odds = (set_odds & 0x33333333) + ((set_odds >> 2) & 0x33333333);
  missings = (missings & 0x33333333) + ((missings >> 2) & 0x33333333);
  set_odds += (set_evens & 0x33333333) + ((set_evens >> 2) & 0x33333333);
  accm = (missings + (missings >> 4)) & 0x0f0f0f0f;
  acc = (set_odds & 0x0f0f0f0f) + ((set_odds >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missings = loader3 & ((loader & (~loader2)) | loader4);
  set_odds = loader3 & loader2 & (~loader4);
  set_evens = loader & set_odds;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missings += loader3 & ((loader & (~loader2)) | loader4);
  loader3 &= loader2 & (~loader4);
  set_odds += loader3;
  set_evens += loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *male_vec++;
  missings += loader3 & ((loader & (~loader2)) | loader4);
  loader3 &= loader2 & (~loader4);
  set_odds += loader3;
  set_evens += loader3 & loader;

  set_odds = (set_odds & 0x33333333) + ((set_odds >> 2) & 0x33333333);
  missings = (missings & 0x33333333) + ((missings >> 2) & 0x33333333);
  set_odds += (set_evens & 0x33333333) + ((set_evens >> 2) & 0x33333333);
  accm += (missings + (missings >> 4)) & 0x0f0f0f0f;
  acc += (set_odds & 0x0f0f0f0f) + ((set_odds >> 4) & 0x0f0f0f0f);
  *missing_ctp += (accm * 0x01010101) >> 24;
  *set_ctp += (acc * 0x01010101) >> 24;
}

void count_set_freq_y_12(uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* nonmale_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t loader4 = *nonmale_vec++;
  uintptr_t sets1 = loader3 & loader & loader2 & (~loader4);
  uintptr_t missings1 = loader3 & (loader4 | (loader ^ loader2));
  uintptr_t sets2;
  uintptr_t missings2;
  uintptr_t acc;
  uintptr_t accm;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));
  sets1 = (sets1 & 0x33333333) + ((sets1 >> 2) & 0x33333333);
  missings1 = (missings1 & 0x33333333) + ((missings1 >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 = loader3 & loader & loader2 & (~loader4);
  missings2 = loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));
  sets1 += (sets2 & 0x33333333) + ((sets2 >> 2) & 0x33333333);
  missings1 += (missings2 & 0x33333333) + ((missings2 >> 2) & 0x33333333);
  acc = (sets1 & 0x0f0f0f0f) + ((sets1 >> 4) & 0x0f0f0f0f);
  accm = (missings1 & 0x0f0f0f0f) + ((missings1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 = loader3 & loader & loader2 & (~loader4);
  missings1 = loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));
  sets1 = (sets1 & 0x33333333) + ((sets1 >> 2) & 0x33333333);
  missings1 = (missings1 & 0x33333333) + ((missings1 >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 = loader3 & loader & loader2 & (~loader4);
  missings2 = loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));
  sets1 += (sets2 & 0x33333333) + ((sets2 >> 2) & 0x33333333);
  missings1 += (missings2 & 0x33333333) + ((missings2 >> 2) & 0x33333333);
  acc += (sets1 & 0x0f0f0f0f) + ((sets1 >> 4) & 0x0f0f0f0f);
  accm += (missings1 & 0x0f0f0f0f) + ((missings1 >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}
#endif

void vec_set_freq(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  // Assuming include_vec describes e.g. cases, and an autosomal marker, this
  // counts the number of case set alleles loaded in lptr[], as well as the
  // number of cases with missing genotype info.
  // See single_marker_freqs_and_hwe() for discussion.
  // missing count: popcount2(genotype & (~(genotype >> 1)) & 0x5555...)
  // set allele count: popcount(genotype) - missing count
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t missing_incr;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_6x_end;
  indiv_ctl2 -= indiv_ctl2 % 6;
  while (indiv_ctl2 >= 60) {
    cur_decr = 60;
  vec_set_freq_loop:
    lptr_6x_end = &(lptr[cur_decr]);
    count_set_freq_60v((__m128i*)lptr, (__m128i*)lptr_6x_end, (__m128i*)include_vec, &acc, &accm);
    lptr = lptr_6x_end;
    include_vec = &(include_vec[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto vec_set_freq_loop;
  }
#else
  uintptr_t* lptr_six_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 6)]);
  while (lptr < lptr_six_end) {
    count_set_freq_6(lptr, include_vec, &acc, &accm);
    lptr = &(lptr[6]);
    include_vec = &(include_vec[6]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *include_vec++;
    missing_incr = popcount2_long(loader & (~(loader >> 1)) & loader2);
    accm += missing_incr;
    acc += popcount_long(loader & (loader2 * 3)) - missing_incr;
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void vec_set_freq_x(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  // diploid counting for nonmales, haploid counting for males
  // missing_ct := male_obs + male_missing + 2 * female_missing
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
  uintptr_t missing_incr;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_6x_end;
  indiv_ctl2 -= indiv_ctl2 % 6;
  while (indiv_ctl2 >= 60) {
    cur_decr = 60;
  vec_set_freq_loop:
    lptr_6x_end = &(lptr[cur_decr]);
    count_set_freq_x_60v((__m128i*)lptr, (__m128i*)lptr_6x_end, (__m128i*)include_vec, (__m128i*)male_vec, &acc, &accm);
    lptr = lptr_6x_end;
    include_vec = &(include_vec[cur_decr]);
    male_vec = &(male_vec[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto vec_set_freq_loop;
  }
#else
  uintptr_t* lptr_six_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 6)]);
  while (lptr < lptr_six_end) {
    count_set_freq_x_6(lptr, include_vec, male_vec, &acc, &accm);
    lptr = &(lptr[6]);
    include_vec = &(include_vec[6]);
    male_vec = &(male_vec[6]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = loader >> 1;
    loader3 = *include_vec++;
    loader4 = loader3 & (~(*male_vec));
    missing_incr = popcount2_long(loader & (~loader2) & loader4);
    accm += 2 * missing_incr;
    acc += popcount_long(loader & (loader4 * 3)) - missing_incr;

    loader4 = loader3 & (*male_vec++);
    acc += popcount2_long(loader & loader2 & loader4);
    accm += popcount_long(((loader ^ loader2) & loader4) | (loader4 << 1));
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void vec_set_freq_y(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* nonmale_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  indiv_ctl2 -= indiv_ctl2 % 12;
  while (indiv_ctl2 >= 120) {
    cur_decr = 120;
  vec_set_freq_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_set_freq_y_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)include_vec, (__m128i*)nonmale_vec, &acc, &accm);
    lptr = lptr_12x_end;
    include_vec = &(include_vec[cur_decr]);
    nonmale_vec = &(nonmale_vec[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto vec_set_freq_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 12)]);
  while (lptr < lptr_twelve_end) {
    count_set_freq_y_12(lptr, include_vec, nonmale_vec, &acc, &accm);
    lptr = &(lptr[12]);
    include_vec = &(include_vec[12]);
    nonmale_vec = &(nonmale_vec[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = loader >> 1;
    loader3 = *include_vec++;
    loader4 = *nonmale_vec++;
    acc += popcount2_long(loader & loader2 & loader3 & (~loader4));
    accm += popcount2_long(loader3 & ((loader ^ loader2) | loader4));
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void vec_set_freq_xx(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp) {
  // diploid counting for nonmales, males count as missing
  // set_ct = popcount(genotype & (((genotype >> 1) & nonmale & 0x555...) * 3))
  // missing_ct = popcount2(((genotype & (~genotype >> 1)) | male) & 0x5555...)
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_6x_end;
  indiv_ctl2 -= indiv_ctl2 % 6;
  while (indiv_ctl2 >= 60) {
    cur_decr = 60;
  vec_set_freq_loop:
    lptr_6x_end = &(lptr[cur_decr]);
    count_set_freq_xx_60v((__m128i*)lptr, (__m128i*)lptr_6x_end, (__m128i*)include_vec, (__m128i*)male_vec, &acc, &accm);
    lptr = lptr_6x_end;
    include_vec = &(include_vec[cur_decr]);
    male_vec = &(male_vec[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto vec_set_freq_loop;
  }
#else
  uintptr_t* lptr_six_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 6)]);
  while (lptr < lptr_six_end) {
    count_set_freq_xx_6(lptr, include_vec, male_vec, &acc, &accm);
    lptr = &(lptr[6]);
    include_vec = &(include_vec[6]);
    male_vec = &(male_vec[6]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = loader >> 1;
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    acc += popcount_long(loader & (3 * (loader2 & (~loader4) & loader3)));
    accm += popcount2_long(((loader & (~loader2)) | loader4) & loader3);
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void vec_3freq(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uint32_t* missing_ctp, uint32_t* het_ctp, uint32_t* homset_ctp) {
  // generic routine for getting all counts.
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uint32_t acc_even = 0;
  uint32_t acc_odd = 0;
  uint32_t acc_and = 0;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  indiv_ctl2 -= indiv_ctl2 % 12;
  while (indiv_ctl2 >= 120) {
    cur_decr = 120;
  vec_homset_freq_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)include_vec, &acc_even, &acc_odd, &acc_and);
    lptr = lptr_12x_end;
    include_vec = &(include_vec[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto vec_homset_freq_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 12)]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, include_vec, &acc_even, &acc_odd, &acc_and);
    lptr = &(lptr[12]);
    include_vec = &(include_vec[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *include_vec++;
    loader3 = loader2 & (loader >> 1);
    acc_even += popcount2_long(loader & loader2);
    acc_odd += popcount2_long(loader3);
    acc_and += popcount2_long(loader & loader3);
  }
  *missing_ctp = acc_even - acc_and;
  *het_ctp = acc_odd - acc_and;
  *homset_ctp = acc_and;
}

void vec_3freq_xx(uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* missing_ctp, uint32_t* het_ctp, uint32_t* homset_ctp) {
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
  uint32_t accm = 0;
  uint32_t acco = 0;
  uint32_t acch = 0;
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  indiv_ctl2 -= indiv_ctl2 % 12;
  while (indiv_ctl2 >= 120) {
    cur_decr = 120;
  vec_homset_freq_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_xx_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)include_vec, (__m128i*)male_vec, &accm, &acco, &acch);
    lptr = lptr_12x_end;
    include_vec = &(include_vec[cur_decr]);
    male_vec = &(male_vec[cur_decr]);
    indiv_ctl2 -= cur_decr;
  }
  if (indiv_ctl2) {
    cur_decr = indiv_ctl2;
    goto vec_homset_freq_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[indiv_ctl2 - (indiv_ctl2 % 12)]);
  while (lptr < lptr_twelve_end) {
    count_3freq_xx_12(lptr, include_vec, male_vec, &accm, &acco, &acch);
    lptr = &(lptr[12]);
    include_vec = &(include_vec[12]);
    male_vec = &(male_vec[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = loader >> 1;
    loader3 = *include_vec++;
    loader4 = *male_vec++;
    accm += popcount2_long(loader3 & (loader4 | (loader & (~loader2))));
    loader3 &= (~loader4) & loader2;
    acco += popcount2_long(loader3);
    acch += popcount2_long(loader3 & loader);
  }
  *missing_ctp = accm;
  *het_ctp = acco - acch;
  *homset_ctp = acch;
}

uint32_t count_chrom_markers(Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uintptr_t* marker_exclude) {
  uint32_t min_idx;
  uint32_t max_idx;
  uint32_t min_idxl;
  uint32_t min_idxlr;
  uint32_t max_idxl;
  uint32_t max_idxlr;
  uint32_t ct;
  if (!((chrom_info_ptr->chrom_mask >> chrom_idx) & 1LLU)) {
    return 0;
  }
  min_idx = chrom_info_ptr->chrom_start[chrom_idx];
  max_idx = chrom_info_ptr->chrom_end[chrom_idx];
  min_idxl = min_idx / BITCT;
  min_idxlr = min_idx & (BITCT - 1);
  max_idxl = max_idx / BITCT;
  max_idxlr = max_idx & (BITCT - 1);
  if (min_idxl == max_idxl) {
    return max_idx - min_idx - popcount_long(marker_exclude[min_idxl] & ((ONELU << max_idxlr) - (ONELU << min_idxlr)));
  } else {
    ct = 0;
    if (min_idxlr) {
      ct = popcount_long(marker_exclude[min_idxl++] >> min_idxlr);
    }
    if (max_idxl > min_idxl) {
      ct += popcount_longs(marker_exclude, min_idxl, max_idxl);
    }
    if (max_idxlr) {
      ct += popcount_long(marker_exclude[max_idxl] & ((ONELU << max_idxlr) - ONELU));
    }
    return max_idx - min_idx - ct;
  }
}

uint32_t count_non_autosomal_markers(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, uint32_t count_x) {
  uint32_t species = chrom_info_ptr->species;
  uint64_t hapmask = species_haploid_mask[species];
  uint32_t max_chrom = species_max_code[species];
  uint32_t ct = 0;
  uint32_t cur_chrom = 0;
  int32_t x_code = species_x_code[species];
  for (; cur_chrom <= max_chrom; cur_chrom++) {
    if ((hapmask >> cur_chrom) & 1LLU) {
      if (count_x || (cur_chrom != (uint32_t)x_code)) {
	ct += count_chrom_markers(chrom_info_ptr, cur_chrom, marker_exclude);
      }
    }
  }
  return ct;
}

uint32_t block_load_autosomal(FILE* bedfile, int32_t bed_offset, uintptr_t* marker_exclude, uint32_t marker_ct_autosomal, uint32_t block_max_size, uintptr_t unfiltered_indiv_ct4, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_weights, unsigned char* readbuf, uint32_t* chrom_fo_idx_ptr, uintptr_t* marker_uidx_ptr, uintptr_t* marker_idx_ptr, uint32_t* block_size_ptr, double* set_allele_freq_buf, float* set_allele_freq_buf_fl, uint32_t* wtbuf) {
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uintptr_t marker_idx = *marker_idx_ptr;
  uint32_t chrom_fo_idx = *chrom_fo_idx_ptr;
  uint32_t chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
  uint32_t markers_read = 0;
  uint32_t is_x;
  uint32_t is_haploid;

  if (block_max_size > marker_ct_autosomal - marker_idx) {
    block_max_size = marker_ct_autosomal - marker_idx;
  }
  while (markers_read < block_max_size) {
    if (is_set(marker_exclude, marker_uidx)) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	return RET_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      while (1) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	if (!is_haploid) {
	  break;
	}
	marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  return RET_READ_FAIL;
	}
      }
    }
    if (fread(&(readbuf[markers_read * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
      return RET_READ_FAIL;
    }
    if (set_allele_freq_buf) {
      set_allele_freq_buf[markers_read] = set_allele_freqs[marker_uidx];
    } else if (set_allele_freq_buf_fl) {
      set_allele_freq_buf_fl[markers_read] = (float)set_allele_freqs[marker_uidx];
    }
    if (wtbuf) {
      wtbuf[markers_read] = marker_weights[marker_idx];
    }
    markers_read++;
    marker_idx++;
    marker_uidx++;
  }

  *chrom_fo_idx_ptr = chrom_fo_idx;
  *marker_uidx_ptr = marker_uidx;
  *marker_idx_ptr = marker_idx;
  *block_size_ptr = markers_read;
  return 0;
}

void exclude_to_vec_include(uintptr_t unfiltered_indiv_ct, uintptr_t* include_arr, uintptr_t* exclude_arr) {
  uint32_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = *exclude_arr++;
    ulkk = FIVEMASK;
    ulmm = FIVEMASK;
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *include_arr++ = ulkk;
    *include_arr++ = ulmm;
  } while (--unfiltered_indiv_ctl);
  ulii = unfiltered_indiv_ct & (BITCT - 1);
  if (ulii) {
    include_arr--;
    if (ulii < BITCT2) {
      *include_arr-- = 0;
    } else {
      ulii -= BITCT2;
    }
    *include_arr &= (ONELU << (ulii * 2)) - ONELU;
  }
}

void vec_init_invert(uintptr_t entry_ct, uintptr_t* target_arr, uintptr_t* source_arr) {
  // initializes a half-bitfield as the inverse of another
  // assumes target_arr and source_arr are 16-byte aligned, and vec_entry_ct is
  // even and positive.
  uint32_t vec_wsize = 2 * ((entry_ct + (BITCT - 1)) / BITCT);
  uintptr_t* second_to_last = &(target_arr[vec_wsize - 2]);
  uint32_t rem = entry_ct & (BITCT - 1);
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* tptr = (__m128i*)target_arr;
  __m128i* sptr = (__m128i*)source_arr;
  __m128i* tptr_end = (__m128i*)(&(target_arr[vec_wsize]));
  do {
    *tptr++ = _mm_andnot_si128(*sptr++, m1);
  } while (tptr < tptr_end);
#else
  uintptr_t* tptr_end = &(target_arr[vec_wsize]);
  do {
    *target_arr++ = FIVEMASK & (~(*source_arr++));
  } while (target_arr < tptr_end);
#endif
  if (rem > BITCT2) {
    second_to_last[1] &= (~ZEROLU) >> ((BITCT - rem) * 2);
  } else if (rem) {
    *second_to_last &= (~ZEROLU) >> ((BITCT2 - rem) * 2);
    second_to_last[1] = 0;
  }
}

void vec_init_andnot(uintptr_t vec_wsize, uintptr_t* target_arr, uintptr_t* source_arr, uintptr_t* exclude_arr) {
  // initializes a half-bitfield as source_arr ANDNOT exclude_arr
#ifdef __LP64__
  __m128i* tptr = (__m128i*)target_arr;
  __m128i* sptr = (__m128i*)source_arr;
  __m128i* xptr = (__m128i*)exclude_arr;
  __m128i* tptr_end = (__m128i*)(&(target_arr[vec_wsize]));
  do {
    *tptr++ = _mm_andnot_si128(*xptr++, *sptr++);
  } while (tptr < tptr_end);
#else
  uintptr_t* tptr_end = &(target_arr[vec_wsize]);
  do {
    *target_arr++ = (*source_arr++) & (~(*exclude_arr++));
  } while (target_arr < tptr_end);
#endif
}

void vec_include_mask_in(uintptr_t unfiltered_indiv_ct, uintptr_t* include_arr, uintptr_t* mask_arr) {
  uint32_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = ~(*mask_arr++);
    ulkk = *include_arr;
    ulmm = include_arr[1];
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *include_arr++ = ulkk;
    *include_arr++ = ulmm;
  } while (--unfiltered_indiv_ctl);
}

void vec_include_mask_out(uintptr_t unfiltered_indiv_ct, uintptr_t* include_arr, uintptr_t* mask_arr) {
  uint32_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = *mask_arr++;
    ulkk = *include_arr;
    ulmm = include_arr[1];
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *include_arr++ = ulkk;
    *include_arr++ = ulmm;
  } while (--unfiltered_indiv_ctl);
}

void vec_include_mask_out_intersect(uintptr_t unfiltered_indiv_ct, uintptr_t* include_arr, uintptr_t* mask_arr, uintptr_t* mask2_arr) {
  uint32_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = (*mask_arr++) & (*mask2_arr++);
    ulkk = *include_arr;
    ulmm = include_arr[1];
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *include_arr++ = ulkk;
    *include_arr++ = ulmm;
  } while (--unfiltered_indiv_ctl);
}

void hh_reset(unsigned char* loadbuf, uintptr_t* indiv_include2, uintptr_t unfiltered_indiv_ct) {
  uintptr_t indiv_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_indiv_ct + 3) / 4]);
  unsigned char* iicp;
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_indiv_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  uint32_t* indiv_include2_alias32;
  __m128i* loadbuf_alias;
  __m128i* iivp;
  __m128i vii;
  __m128i vjj;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    iivp = (__m128i*)indiv_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 64;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      vii = *loadbuf_alias;
      vjj = _mm_and_si128(_mm_andnot_si128(vii, _mm_srli_epi64(vii, 1)), *iivp++);
      *loadbuf_alias++ = _mm_sub_epi64(vii, vjj);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    iicp = (unsigned char*)iivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    indiv_include2_alias32 = (uint32_t*)indiv_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / BITCT2;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *loadbuf_alias32;
      ujj = ((uii >> 1) & (~uii)) & (*indiv_include2_alias32++);
      *loadbuf_alias32++ = uii - ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    iicp = (unsigned char*)indiv_include2_alias32;
  } else {
    iicp = (unsigned char*)indiv_include2;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / BITCT2;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *loadbuf_alias32;
      ujj = ((uii >> 1) & (~uii)) & (*indiv_include2++);
      *loadbuf_alias32++ = uii - ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  iicp = (unsigned char*)indiv_include2;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = ((ucc >> 1) & (~ucc)) & (*iicp++);
    *loadbuf++ = ucc - ucc2;
  }
}

static uint32_t g_pct;
static uintptr_t g_dw_indiv1idx;
static uintptr_t g_dw_indiv2idx;
static uintptr_t g_dw_min_indiv;
static uintptr_t g_dw_max_indiv1idx;
static uint64_t g_dw_start_offset;
static uint64_t g_dw_hundredth;
static double* g_dw_dists;
static double* g_dw_dist_ptr;
static unsigned char* g_dw_membuf;
static double g_dw_half_marker_ct_recip;

uint32_t distance_d_write_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx + 1 < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_dw_dist_ptr++, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_tri_emitn_ret;
      }
    }
    if (g_dw_indiv2idx + 1 == g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_dw_dist_ptr++, '\n');
    }
    if ((((uint64_t)g_dw_indiv1idx) * (g_dw_indiv1idx + 1) / 2 - g_dw_start_offset) >= g_dw_hundredth * g_pct) {
      g_pct = (((uint64_t)g_dw_indiv1idx) * (g_dw_indiv1idx + 1) / 2 - g_dw_start_offset) / g_dw_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv1idx++;
    g_dw_indiv2idx = 0;
  }
 distance_d_write_tri_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t ulii;
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_dw_dist_ptr++, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_sq0_emitn_ret;
      }
    }
    ulii = 1 + (((uintptr_t)(readbuf_end - sptr_cur)) / 2);
    if (ulii < (g_indiv_ct - g_dw_indiv2idx)) {
      sptr_cur = memcpya(sptr_cur, &(g_dw_membuf[1]), 2 * ulii);
      g_dw_indiv2idx += ulii;
      goto distance_d_write_sq0_emitn_ret;
    }
    ulii = g_indiv_ct - g_dw_indiv2idx;
    sptr_cur = memcpyax(sptr_cur, &(g_dw_membuf[1]), 2 * ulii - 1, '\n');
    g_dw_indiv1idx++;
    if ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_dw_max_indiv1idx - g_dw_min_indiv)) {
      g_pct = ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU) / (g_dw_max_indiv1idx - g_dw_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_sq0_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_dw_dist_ptr++, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_sq_emitn_ret;
      }
    }
    if (g_dw_indiv2idx == g_dw_indiv1idx) {
      *sptr_cur++ = '0';
      g_dw_indiv2idx++;
    }
    while (g_dw_indiv2idx < g_indiv_ct) {
      *sptr_cur = '\t';
      sptr_cur = double_g_write(&(sptr_cur[1]), g_dw_dists[((g_dw_indiv2idx * (g_dw_indiv2idx - 1)) / 2) + g_dw_indiv1idx]);
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    g_dw_indiv1idx++;
    if ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_dw_max_indiv1idx - g_dw_min_indiv)) {
      g_pct = ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU) / (g_dw_max_indiv1idx - g_dw_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_sq_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_ibs_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, 1.0 - (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_tri_emitn_ret;
      }
    }
    if (g_dw_indiv2idx == g_dw_indiv1idx) {
      sptr_cur = memcpya(sptr_cur, "1\n", 2);
    }
    g_dw_indiv1idx++;
    if ((((uint64_t)g_dw_indiv1idx) * (g_dw_indiv1idx + 1) / 2 - g_dw_start_offset) >= g_dw_hundredth * g_pct) {
      g_pct = (((uint64_t)g_dw_indiv1idx) * (g_dw_indiv1idx + 1) / 2 - g_dw_start_offset) / g_dw_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_ibs_tri_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_ibs_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t ulii;
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, 1.0 - (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq0_emitn_ret;
      }
    }
    if (g_dw_indiv2idx == g_dw_indiv1idx) {
      *sptr_cur++ = '1';
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq0_emitn_ret;
      }
    }
    ulii = (1 + ((uintptr_t)(readbuf_end - sptr_cur))) / 2;
    if (ulii < (g_indiv_ct - g_dw_indiv2idx)) {
      sptr_cur = memcpya(sptr_cur, g_dw_membuf, 2 * ulii);
      g_dw_indiv2idx += ulii;
      goto distance_d_write_ibs_sq0_emitn_ret;
    }
    ulii = g_indiv_ct - g_dw_indiv2idx;
    sptr_cur = memcpyax(sptr_cur, g_dw_membuf, 2 * ulii, '\n');
    g_dw_indiv1idx++;
    if ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_dw_max_indiv1idx - g_dw_min_indiv)) {
      g_pct = ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU) / (g_dw_max_indiv1idx - g_dw_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_ibs_sq0_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_ibs_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, 1.0 - (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq_emitn_ret;
      }
    }
    if (g_dw_indiv2idx == g_dw_indiv1idx) {
      *sptr_cur++ = '1';
      g_dw_indiv2idx++;
    }
    while (g_dw_indiv2idx < g_indiv_ct) {
      *sptr_cur = '\t';
      sptr_cur = double_g_write(&(sptr_cur[1]), 1.0 - (g_dw_dists[((g_dw_indiv2idx * (g_dw_indiv2idx - 1)) / 2) + g_dw_indiv1idx]) * g_dw_half_marker_ct_recip);
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    g_dw_indiv1idx++;
    if ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_dw_max_indiv1idx - g_dw_min_indiv)) {
      g_pct = ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU) / (g_dw_max_indiv1idx - g_dw_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_ibs_sq_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_1mibs_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx + 1 < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_tri_emitn_ret;
      }
    }
    if (g_dw_indiv2idx + 1 == g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\n');
    }
    if ((((uint64_t)g_dw_indiv1idx) * (g_dw_indiv1idx + 1) / 2 - g_dw_start_offset) >= g_dw_hundredth * g_pct) {
      g_pct = (((uint64_t)g_dw_indiv1idx) * (g_dw_indiv1idx + 1) / 2 - g_dw_start_offset) / g_dw_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv1idx++;
    g_dw_indiv2idx = 0;
  }
 distance_d_write_1mibs_tri_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_1mibs_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t ulii;
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_sq0_emitn_ret;
      }
    }
    ulii = 1 + (((uintptr_t)(readbuf_end - sptr_cur)) / 2);
    if (ulii < (g_indiv_ct - g_dw_indiv2idx)) {
      sptr_cur = memcpya(sptr_cur, &(g_dw_membuf[1]), 2 * ulii);
      g_dw_indiv2idx += ulii;
      goto distance_d_write_1mibs_sq0_emitn_ret;
    }
    ulii = g_indiv_ct - g_dw_indiv2idx;
    sptr_cur = memcpyax(sptr_cur, &(g_dw_membuf[1]), 2 * ulii - 1, '\n');
    g_dw_indiv1idx++;
    if ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_dw_max_indiv1idx - g_dw_min_indiv)) {
      g_pct = ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU) / (g_dw_max_indiv1idx - g_dw_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_1mibs_sq0_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_1mibs_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_dw_indiv1idx < g_dw_max_indiv1idx) {
    while (g_dw_indiv2idx < g_dw_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, (*g_dw_dist_ptr++) * g_dw_half_marker_ct_recip, '\t');
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_sq_emitn_ret;
      }
    }
    if (g_dw_indiv2idx == g_dw_indiv1idx) {
      *sptr_cur++ = '0';
      g_dw_indiv2idx++;
    }
    while (g_dw_indiv2idx < g_indiv_ct) {
      *sptr_cur = '\t';
      sptr_cur = double_g_write(&(sptr_cur[1]), (g_dw_dists[((g_dw_indiv2idx * (g_dw_indiv2idx - 1)) / 2) + g_dw_indiv1idx]) * g_dw_half_marker_ct_recip);
      g_dw_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    g_dw_indiv1idx++;
    if ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_dw_max_indiv1idx - g_dw_min_indiv)) {
      g_pct = ((g_dw_indiv1idx - g_dw_min_indiv) * 100LLU) / (g_dw_max_indiv1idx - g_dw_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_dw_indiv2idx = 0;
  }
 distance_d_write_1mibs_sq_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t distance_d_write(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, int32_t dist_calc_type, char* outname, char* outname_end, double* dists, double half_marker_ct_recip, uint32_t indiv_ct, int32_t first_indiv_idx, int32_t end_indiv_idx, int32_t parallel_idx, int32_t parallel_tot, unsigned char* membuf) {
  // membuf assumed to be of at least size indiv_ct * 8.
  int32_t shape = dist_calc_type & DISTANCE_SHAPEMASK;
  int32_t write_alcts = dist_calc_type & DISTANCE_ALCT;
  int32_t write_ibs_matrix = dist_calc_type & DISTANCE_IBS;
  int32_t write_1mibs_matrix = dist_calc_type & DISTANCE_1_MINUS_IBS;
  int32_t retval = 0;
  double dxx;
  double dyy;
  double* dist_ptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t* glptr;
  uint32_t uii;
  uintptr_t indiv_idx_ct;
  int32_t ii;
  int32_t jj;
  char* cptr;
  g_dw_start_offset = ((int64_t)first_indiv_idx * (first_indiv_idx - 1)) / 2;
  indiv_idx_ct = (uintptr_t)(((int64_t)end_indiv_idx * (end_indiv_idx - 1)) / 2 - g_dw_start_offset);
  g_dw_hundredth = 1 + (indiv_idx_ct / 100);
  if (first_indiv_idx == 1) {
    first_indiv_idx = 0;
  }
  if (shape == DISTANCE_SQ0) {
    cptr = (char*)(&ulii);
    for (uii = 0; uii < sizeof(intptr_t); uii += 2) {
      cptr[uii] = '\t';
      cptr[uii + 1] = '0';
    }
    ii = (indiv_ct * 2 + sizeof(intptr_t) - 2) / sizeof(intptr_t);
    glptr = (uintptr_t*)membuf;
    for (jj = 0; jj < ii; jj++) {
      *glptr++ = ulii;
    }
    g_dw_membuf = membuf;
  }
  g_pct = 1;
  if (dist_calc_type & DISTANCE_BIN) {
    if (distance_open(outfile_ptr, outfile2_ptr, outfile3_ptr, outname, outname_end, ".bin", "wb", dist_calc_type, parallel_idx, parallel_tot)) {
      goto distance_d_write_ret_OPEN_FAIL;
    }
    if (shape == DISTANCE_TRI) {
      if (write_alcts) {
	fputs("Writing...", stdout);
	fflush(stdout);
	if (fwrite_checkedz(dists, indiv_idx_ct * sizeof(double), *outfile_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	ulii = 0;
	do {
	  uljj = (indiv_idx_ct * g_pct) / 100L;
	  for (; ulii < uljj; ulii++) {
	    dxx = 1.0 - (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  printf("\rWriting... %d%%", g_pct++);
	  fflush(stdout);
	} while (g_pct <= 100);
	distance_print_done(1, outname, outname_end);
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	ulii = 0;
	do {
	  uljj = (indiv_idx_ct * g_pct) / 100L;
	  for (; ulii < uljj; ulii++) {
	    dxx = (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  printf("\rWriting... %d%%", g_pct++);
	  fflush(stdout);
	} while (g_pct <= 100);
	distance_print_done(2, outname, outname_end);
      }
    } else {
      if (shape == DISTANCE_SQ0) {
	fill_double_zero((double*)membuf, indiv_ct);
      }
      if (write_alcts) {
	dxx = 0.0;
	dist_ptr = dists;
	for (ii = first_indiv_idx; ii < end_indiv_idx; ii++) {
	  if (fwrite_checkedz(dist_ptr, ii * sizeof(double), *outfile_ptr)) {
	    goto distance_d_write_ret_WRITE_FAIL;
	  }
	  dist_ptr = &(dist_ptr[ii]);
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checkedz(membuf, (indiv_ct - ii) * sizeof(double), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix, no need to handle parallel case
	    if (fwrite_checked(&dxx, sizeof(double), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      if (fwrite_checked(&(dists[(ulii * (ulii - 1)) / 2 + ii]), sizeof(double), *outfile_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (int64_t)g_pct * (end_indiv_idx - first_indiv_idx)) {
	    g_pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
	g_pct = 1;
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	dyy = 1.0;
	membuf[1] = '1';
	for (ii = first_indiv_idx; ii < end_indiv_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    dxx = 1.0 - (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checkedz(membuf, (indiv_ct - ii) * sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&dyy, sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      dxx = 1.0 - dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (int64_t)g_pct * (end_indiv_idx - first_indiv_idx)) {
	    g_pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	membuf[1] = '0';
	if (fclose_null(outfile2_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(1, outname, outname_end);
	g_pct = 1;
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	dyy = 0.0;
	for (ii = first_indiv_idx; ii < end_indiv_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    dxx = (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checkedz(membuf, (indiv_ct - ii) * sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&dyy, sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	      dxx = dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_indiv_idx) * 100LL >= (int64_t)g_pct * (end_indiv_idx - first_indiv_idx)) {
	    g_pct = ((ii - first_indiv_idx) * 100LL) / (end_indiv_idx - first_indiv_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile3_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(2, outname, outname_end);
      }
    }
  } else {
    g_dw_min_indiv = first_indiv_idx;
    g_dw_max_indiv1idx = end_indiv_idx;
    g_dw_dists = dists;
    g_dw_half_marker_ct_recip = half_marker_ct_recip;
    if (write_alcts) {
      g_dw_indiv1idx = first_indiv_idx;
      g_dw_indiv2idx = 0;
      g_dw_dist_ptr = dists;
      if (dist_calc_type & DISTANCE_GZ) {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".dist.%u.gz", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".dist.gz");
	}
	if (shape == DISTANCE_SQ) {
	  parallel_compress(outname, distance_d_write_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  parallel_compress(outname, distance_d_write_sq0_emitn);
	} else {
	  parallel_compress(outname, distance_d_write_tri_emitn);
	}
      } else {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".dist.%u", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".dist");
	}
	if (shape == DISTANCE_SQ) {
	  retval = write_uncompressed(outname, distance_d_write_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  retval = write_uncompressed(outname, distance_d_write_sq0_emitn);
	} else {
	  retval = write_uncompressed(outname, distance_d_write_tri_emitn);
	}
	if (retval) {
	  goto distance_d_write_ret_1;
	}
      }
      sprintf(logbuf, "\rDistances (allele counts) written to %s.\n", outname);
      logprintb();
      g_pct = 1;
    }
    if (write_1mibs_matrix) {
      g_dw_indiv1idx = first_indiv_idx;
      g_dw_indiv2idx = 0;
      g_dw_dist_ptr = dists;
      if (dist_calc_type & DISTANCE_GZ) {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mdist.%u.gz", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mdist.gz");
	}
	if (shape == DISTANCE_SQ) {
	  parallel_compress(outname, distance_d_write_1mibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  parallel_compress(outname, distance_d_write_1mibs_sq0_emitn);
	} else {
	  parallel_compress(outname, distance_d_write_1mibs_tri_emitn);
	}
      } else {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mdist.%u", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mdist");
	}
	if (shape == DISTANCE_SQ) {
	  retval = write_uncompressed(outname, distance_d_write_1mibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  retval = write_uncompressed(outname, distance_d_write_1mibs_sq0_emitn);
	} else {
	  retval = write_uncompressed(outname, distance_d_write_1mibs_tri_emitn);
	}
	if (retval) {
	  goto distance_d_write_ret_1;
	}
      }
      sprintf(logbuf, "\rDistances (proportions) written to %s.\n", outname);
      logprintb();
      g_pct = 1;
    }
    if (write_ibs_matrix) {
      g_dw_indiv1idx = first_indiv_idx;
      g_dw_indiv2idx = 0;
      g_dw_dist_ptr = dists;
      g_dw_start_offset = ((int64_t)first_indiv_idx * (first_indiv_idx + 1)) / 2;
      g_dw_hundredth = 1 + ((((int64_t)end_indiv_idx * (end_indiv_idx + 1)) / 2 - g_dw_start_offset) / 100);
      if (dist_calc_type & DISTANCE_GZ) {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mibs.%u.gz", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mibs.gz");
	}
	if (shape == DISTANCE_SQ) {
	  parallel_compress(outname, distance_d_write_ibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  parallel_compress(outname, distance_d_write_ibs_sq0_emitn);
	} else {
	  parallel_compress(outname, distance_d_write_ibs_tri_emitn);
	}
      } else {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mibs.%u", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mibs");
	}
	if (shape == DISTANCE_SQ) {
	  retval = write_uncompressed(outname, distance_d_write_ibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  retval = write_uncompressed(outname, distance_d_write_ibs_sq0_emitn);
	} else {
	  retval = write_uncompressed(outname, distance_d_write_ibs_tri_emitn);
	}
	if (retval) {
	  goto distance_d_write_ret_1;
	}
      }
      sprintf(logbuf, "\rIBS matrix written to %s.\n", outname);
      logprintb();
    }
  }
  while (0) {
  distance_d_write_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  distance_d_write_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 distance_d_write_ret_1:
  return retval;
}

void collapse_arr(char* item_arr, int32_t fixed_item_len, uintptr_t* exclude_arr, int32_t exclude_arr_size) {
  // collapses array of fixed-length items
  int32_t ii = 0;
  int32_t jj;
  while ((ii < exclude_arr_size) && (!is_set(exclude_arr, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < exclude_arr_size) {
    if (!is_set(exclude_arr, ii)) {
      memcpy(&(item_arr[(jj++) * fixed_item_len]), &(item_arr[ii * fixed_item_len]), fixed_item_len);
    }
  }
}

void collapse_bitarr(uintptr_t* bitarr, uintptr_t* exclude_arr, uint32_t orig_ct) {
  uint32_t uii = 0;
  uint32_t ujj;
  while ((uii < orig_ct) && (!is_set(exclude_arr, uii))) {
    uii++;
  }
  ujj = uii;
  while (++uii < orig_ct) {
    if (!is_set(exclude_arr, uii)) {
      if (is_set(bitarr, uii)) {
        // set bit jj
        bitarr[ujj / BITCT] |= (ONELU << (ujj % BITCT));
      } else {
	bitarr[ujj / BITCT] &= (~(ONELU << (ujj % BITCT)));
      }
      ujj++;
    }
  }
}

void collapse_bitarr_incl(uintptr_t* bitarr, uintptr_t* include_arr, uint32_t orig_ct) {
  uint32_t uii = 0;
  uint32_t ujj;
  while ((uii < orig_ct) && is_set(include_arr, uii)) {
    uii++;
  }
  ujj = uii;
  while (++uii < orig_ct) {
    if (is_set(include_arr, uii)) {
      if (is_set(bitarr, uii)) {
        // set bit jj
        bitarr[ujj / BITCT] |= (ONELU << (ujj % BITCT));
      } else {
	bitarr[ujj / BITCT] &= (~(ONELU << (ujj % BITCT)));
      }
      ujj++;
    }
  }
}

double rand_unif(void) {
  return (sfmt_genrand_uint32(&sfmt) + 0.5) * RECIP_2_32;
}

// implementation used in PLINK stats.cpp
double normdist(double zz) {
  double sqrt2pi = 2.50662827463;
  double t0;
  double z1;
  double p0;
  t0 = 1 / (1 + 0.2316419 * fabs(zz));
  z1 = exp(-0.5 * zz * zz) / sqrt2pi;
  p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
 return zz >= 0 ? 1 - p0 : p0;
}

double rand_normal(double* secondval_ptr) {
  // N(0, 1)
  double dxx = sqrt(-2 * log(rand_unif()));
  double dyy = 2 * PI * rand_unif();
  *secondval_ptr = dxx * cos(dyy);
  return dxx * sin(dyy);
}

void pick_d(unsigned char* cbuf, uint32_t ct, uint32_t dd) {
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  memset(cbuf, 0, ct);
  ukk = (uint32_t)(4294967296LLU % ct);
  for (uii = 0; uii < dd; uii++) {
    do {
      do {
        ujj = sfmt_genrand_uint32(&sfmt);
      } while (ujj < ukk);
      ujj %= ct;
    } while (cbuf[ujj]);
    cbuf[ujj] = 1;
  }
}

void pick_d_small(unsigned char* tmp_cbuf, int32_t* ibuf, uint32_t ct, uint32_t dd) {
  uint32_t uii;
  pick_d(tmp_cbuf, ct, dd);
  for (uii = 0; uii < ct; uii++) {
    if (tmp_cbuf[uii]) {
      *ibuf++ = uii;
    }
  }
  *ibuf = ct;
}

void print_pheno_stdev(double* pheno_d, uint32_t indiv_ct) {
  double reg_tot_x = 0.0;
  double reg_tot_xx = 0.0;
  double dxx;
  uint32_t uii;
  for (uii = 0; uii < indiv_ct; uii++) {
    dxx = pheno_d[uii];
    reg_tot_x += dxx;
    reg_tot_xx += dxx * dxx;
  }
  sprintf(logbuf, "Phenotype stdev: %g\n", sqrt((reg_tot_xx - reg_tot_x * reg_tot_x / indiv_ct) / (indiv_ct - 1)));
  logprintb();
}

uint32_t set_default_jackknife_d(uint32_t ct) {
  uint32_t dd = (uint32_t)pow((double)ct, 0.600000000001);
  sprintf(logbuf, "Setting d=%u for jackknife.\n", dd);
  logprintb();
  return dd;
}

void generate_perm1_interleaved(uint32_t tot_ct, uint32_t set_ct, uintptr_t perm_idx, uintptr_t perm_ct, uintptr_t* perm_buf) {
  uintptr_t tot_ctl = (tot_ct + (BITCT - 1)) / BITCT;
  uintptr_t tot_rem = tot_ct & (BITCT - 1);
  uint32_t tot_quotient = (uint32_t)(4294967296LLU / tot_ct);
  uint32_t upper_bound = tot_ct * tot_quotient;
  uintptr_t uljj = perm_ct - perm_idx;
  uint32_t totq_preshift;
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uintptr_t* pbptr;
  uint32_t num_set;
  uint32_t urand;
  uintptr_t ulii;
  // seeing as how we're gonna divide by the same number a billion times or so,
  // it just might be worth optimizing that division...
  magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
  if (set_ct * 2 < tot_ct) {
    for (ulii = 0; ulii < tot_ctl; ulii++) {
      fill_ulong_zero(&(perm_buf[perm_idx + (ulii * perm_ct)]), uljj);
    }
    for (; perm_idx < perm_ct; perm_idx++) {
      pbptr = &(perm_buf[perm_idx]);
      for (num_set = 0; num_set < set_ct; num_set++) {
	do {
	  do {
	    urand = sfmt_genrand_uint32(&sfmt);
	  } while (urand >= upper_bound);
	  // this is identical to ulii = urand / tot_quotient
	  ulii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  uljj = ulii / BITCT;
	  ulii &= (BITCT - 1);
	} while ((pbptr[uljj * perm_ct] >> ulii) & 1);
	pbptr[uljj * perm_ct] |= (ONELU << ulii);
      }
    }
  } else {
    for (ulii = 0; ulii < tot_ctl; ulii++) {
      fill_ulong_one(&(perm_buf[perm_idx + (ulii * perm_ct)]), uljj);
    }
    // "set" has reversed meaning here
    set_ct = tot_ct - set_ct;
    for (; perm_idx < perm_ct; perm_idx++) {
      pbptr = &(perm_buf[perm_idx]);
      for (num_set = 0; num_set < set_ct; num_set++) {
	do {
	  do {
	    urand = sfmt_genrand_uint32(&sfmt);
	  } while (urand >= upper_bound);
	  ulii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  uljj = ulii / BITCT;
	  ulii &= (BITCT - 1);
	} while (!((pbptr[uljj * perm_ct] >> ulii) & 1));
	pbptr[uljj * perm_ct] &= ~(ONELU << ulii);
      }
    }
    if (tot_rem) {
      uljj = (~ZEROLU) >> (BITCT - tot_rem);
      pbptr = &(perm_buf[(tot_ctl - 1) * perm_ct + perm_idx]);
      for (ulii = perm_idx; ulii < perm_ct; ulii++) {
	*pbptr &= uljj;
	pbptr++;
      }
    }
  }
}

void join_threads(pthread_t* threads, uint32_t ctp1) {
  if (!(--ctp1)) {
    return;
  }
#if _WIN32
  WaitForMultipleObjects(ctp1, threads, 1, INFINITE);
#else
  uint32_t uii;
  for (uii = 0; uii < ctp1; uii++) {
    pthread_join(threads[uii], NULL);
  }
#endif
}

#if _WIN32
int32_t spawn_threads(pthread_t* threads, unsigned (__stdcall *start_routine)(void*), uintptr_t ct)
#else
int32_t spawn_threads(pthread_t* threads, void* (*start_routine)(void*), uintptr_t ct)
#endif
{
  uintptr_t ulii;
  if (ct == 1) {
    return 0;
  }
  for (ulii = 1; ulii < ct; ulii++) {
#if _WIN32
    threads[ulii - 1] = (HANDLE)_beginthreadex(NULL, 4096, start_routine, (void*)ulii, 0, NULL);
    if (!threads[ulii - 1]) {
      join_threads(threads, ulii);
      return -1;
    }
#else
    if (pthread_create(&(threads[ulii - 1]), NULL, start_routine, (void*)ulii)) {
      join_threads(threads, ulii);
      return -1;
    }
#endif
  }
  return 0;
}

// ----- multithread globals -----
static double* g_pheno_d;
static uintptr_t g_jackknife_iters;
static uint32_t g_jackknife_d;
static double g_reg_tot_xy;
static double g_reg_tot_x;
static double g_reg_tot_y;
static double g_reg_tot_xx;
static double g_reg_tot_yy;
static double* g_jackknife_precomp;
static double* g_dists;
static double g_calc_result[4][MAX_THREADS_P1];
static unsigned char* g_generic_buf;

// double regress_jack(int32_t* ibuf) {
double regress_jack(int32_t* ibuf, double* ret2_ptr) {
  int32_t* iptr = ibuf;
  int32_t* jptr = &(ibuf[g_jackknife_d]);
  uint32_t uii;
  int32_t jj;
  int32_t kk;
  double* dptr;
  double* dptr2;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double neg_tot_yy = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  while (iptr < jptr) {
    dptr2 = &(g_jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_DIST]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  iptr = ibuf;
  for (uii = 1; uii < g_jackknife_d; uii++) {
    jj = *(++iptr);
    dxx1 = g_pheno_d[jj];
    jptr = ibuf;
    dptr = &(g_dists[((intptr_t)jj * (jj - 1)) / 2]);
    while (jptr < iptr) {
      kk = *jptr++;
      dxx = (dxx1 + g_pheno_d[kk]) * 0.5;
      dyy = dptr[kk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = g_reg_tot_y - neg_tot_y;
  dyy = g_indiv_ct - g_jackknife_d;
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((g_reg_tot_xy - neg_tot_xy) - dxx * (g_reg_tot_x - neg_tot_x) / dyy) / ((g_reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = g_reg_tot_x - neg_tot_x;
  return ((g_reg_tot_xy - neg_tot_xy) - dxx * (g_reg_tot_y - neg_tot_y) / dyy) / ((g_reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

THREAD_RET_TYPE regress_jack_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t* ibuf = (int32_t*)(&(g_generic_buf[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t))]));
  unsigned char* cbuf = &(g_generic_buf[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t)) + (g_jackknife_d + 1) * sizeof(int32_t)]);
  uint64_t ulii;
  uint64_t uljj = g_jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  double dxx;
  double ret2;
  for (ulii = 0; ulii < g_jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, g_indiv_ct, g_jackknife_d);
    dxx = regress_jack(ibuf, &ret2);
    // dxx = regress_jack(ibuf);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / g_jackknife_iters;
      printf("\r%" PRIu64 "%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * g_jackknife_iters) / 100;
    }
  }
  g_calc_result[0][tidx] = sum;
  g_calc_result[1][tidx] = sum_sq;
  g_calc_result[2][tidx] = sum2;
  g_calc_result[3][tidx] = sum2_sq;
  THREAD_RETURN;
}

int32_t regress_distance(uint64_t calculation_type, double* dists_local, double* pheno_d_local, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uint32_t thread_ct, uintptr_t regress_iters, uint32_t regress_d) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t ulii;
  uint32_t uii;
  double* dist_ptr;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  double* dptr5;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  double dvv;
  double duu;
  pthread_t threads[MAX_THREADS];

  g_dists = dists_local;
  g_pheno_d = pheno_d_local;

  // beta = (mean(xy) - mean(x)*mean(y)) / (mean(x^2) - mean(x)^2)
  if (unfiltered_indiv_ct != g_indiv_ct) {
    // destructive!  make copy in the future
    collapse_arr((char*)g_pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
  }
  if (!(calculation_type & CALC_REGRESS_REL)) {
    print_pheno_stdev(g_pheno_d, g_indiv_ct);
  }
  ulii = g_indiv_ct;
  ulii = ulii * (ulii - 1) / 2;
  g_reg_tot_xy = 0.0;
  g_reg_tot_x = 0.0;
  g_reg_tot_y = 0.0;
  g_reg_tot_xx = 0.0;
  g_reg_tot_yy = 0.0;
  dptr4 = g_dists;
  dist_ptr = g_pheno_d;
  // Linear regression slope is a function of sum(xy), sum(x), sum(y),
  // sum(x^2), and n.  To speed up the jackknife calculation, we precompute
  // (i) the global xy, x, y, x^2, and
  // (ii) the xy, x, y, x^2 for each row.
  // Then for each delete-d jackknife iteration, we take the global sums,
  // subtract the partial row sums corresponding to the deleted individuals,
  // and then add back the elements in the intersection of two deletions.
  g_jackknife_precomp = (double*)wkspace_alloc(g_indiv_ct * JACKKNIFE_VALS_DIST * sizeof(double));
  if (!g_jackknife_precomp) {
    return RET_NOMEM;
  }
  fill_double_zero(g_jackknife_precomp, g_indiv_ct * JACKKNIFE_VALS_DIST);
  for (uii = 1; uii < g_indiv_ct; uii++) {
    dzz = *(++dist_ptr);
    dptr2 = g_pheno_d;
    dptr3 = &(g_jackknife_precomp[uii * JACKKNIFE_VALS_DIST]);
    dptr5 = g_jackknife_precomp;
    do {
      dxx = (dzz + *dptr2++) * 0.5;
      dyy = (*dptr4++);
      dww = dxx * dyy;
      dvv = dxx * dxx;
      duu = dyy * dyy;
      g_reg_tot_xy += dww;
      *dptr3 += dww;
      *dptr5 += dww;
      dptr5++;
      g_reg_tot_x += dxx;
      dptr3[1] += dxx;
      *dptr5 += dxx;
      dptr5++;
      g_reg_tot_y += dyy;
      dptr3[2] += dyy;
      *dptr5 += dyy;
      dptr5++;
      g_reg_tot_xx += dvv;
      dptr3[3] += dvv;
      *dptr5 += dvv;
      dptr5++;
      g_reg_tot_yy += duu;
      dptr3[4] += duu;
      *dptr5 += duu;
      dptr5++;
    } while (dptr2 < dist_ptr);
  }

  dxx = ulii;
  sprintf(logbuf, "Regression slope (y = genomic distance, x = avg phenotype): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y / dxx) / (g_reg_tot_xx - g_reg_tot_x * g_reg_tot_x / dxx));
  logprintb();
  sprintf(logbuf, "Regression slope (y = avg phenotype, x = genomic distance): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y / dxx) / (g_reg_tot_yy - g_reg_tot_y * g_reg_tot_y / dxx));
  logprintb();

  g_jackknife_iters = (regress_iters + thread_ct - 1) / thread_ct;
  if (regress_d) {
    g_jackknife_d = regress_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(g_indiv_ct);
  }
  g_generic_buf = wkspace_alloc(thread_ct * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t)));
  if (!g_generic_buf) {
    return RET_NOMEM;
  }
  if (spawn_threads(threads, &regress_jack_thread, thread_ct)) {
    logprint(errstr_thread_create);
    return RET_THREAD_CREATE_FAIL;
  }
  ulii = 0;
  regress_jack_thread((void*)ulii);
  dyy = g_calc_result[0][0]; // sum
  dzz = g_calc_result[1][0]; // sum of squares
  dww = g_calc_result[2][0]; // reverse regression sum
  dvv = g_calc_result[3][0]; // reverse regression sum of squares
  join_threads(threads, thread_ct);
  for (uii = 0; uii < thread_ct - 1; uii++) {
    dyy += g_calc_result[0][uii + 1];
    dzz += g_calc_result[1][uii + 1];
    dww += g_calc_result[2][uii + 1];
    dvv += g_calc_result[3][uii + 1];
  }
  regress_iters = g_jackknife_iters * thread_ct;
  putchar('\r');
  sprintf(logbuf, "Jackknife s.e.: %g\n", sqrt((g_indiv_ct / ((double)g_jackknife_d)) * (dzz - dyy * dyy / regress_iters) / (regress_iters - 1)));
  logprintb();
  sprintf(logbuf, "Jackknife s.e. (y = avg phenotype): %g\n", sqrt((g_indiv_ct / ((double)g_jackknife_d)) * (dvv - dww * dww / regress_iters) / (regress_iters - 1)));
  logprintb();
  wkspace_reset(wkspace_mark);
  return 0;
}
