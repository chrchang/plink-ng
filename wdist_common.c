#include "wdist_common.h"
#include "pigz.h"

const char errstr_fopen[] = "Error: Failed to open %s.\n";
const char errstr_append[] = "\nFor more information, try '" PROG_NAME_STR " --help [flag name]' or '" PROG_NAME_STR " --help | more'.\n";
const char errstr_thread_create[] = "\nError: Failed to create thread.\n";
const char cmdline_format_str[] = "\n  " PROG_NAME_STR " [input flag(s)...] {command flag(s)...} {other flag(s)...}\n  " PROG_NAME_STR " --help {flag name(s)...}\n\n";
const char errstr_phenotype_format[] = "Error: Improperly formatted phenotype file.\n";

char tbuf[MAXLINELEN * 4 + 256];

sfmt_t sfmt;

FILE* logfile = NULL;
char logbuf[MAXLINELEN]; // safe sprintf buffer, if one is needed
int32_t debug_on = 0;
int32_t log_failed = 0;
uintptr_t g_indiv_ct;
uint32_t g_thread_ct;

uint32_t push_ll_str(Ll_str** ll_stack_ptr, const char* ss) {
  uint32_t slen = strlen(ss);
  Ll_str* new_ll_str = (Ll_str*)malloc(sizeof(Ll_str) + slen + 1);
  if (!new_ll_str) {
    return 1;
  }
  new_ll_str->next = *ll_stack_ptr;
  memcpy(new_ll_str->ss, ss, slen + 1);
  *ll_stack_ptr = new_ll_str;
  return 0;
}

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

uint32_t match_upper(char* ss, const char* fixed_str) {
  // Returns whether uppercased ss matches nonempty fixed_str.  Assumes
  // fixed_str contains nothing but letters and a null terminator.
  char cc = *fixed_str++;
  do {
    if ((((unsigned char)(*ss++)) & 0xdf) != ((unsigned char)cc)) {
      return 0;
    }
    cc = *fixed_str++;
  } while (cc);
  return !(*ss);
}

uint32_t match_upper_nt(char* ss, const char* fixed_str, uint32_t ct) {
  do {
    if ((((unsigned char)(*ss++)) & 0xdf) != ((unsigned char)(*fixed_str++))) {
      return 0;
    }
  } while (--ct);
  return 1;
}

int32_t atoiz(char* ss, int32_t* sval) {
  // accepts nonnegative integers
  int32_t ii = atoi(ss);
  if ((ii < 1) && ((*ss != '0') || (ss[1] != '\0'))) {
    return -1;
  }
  *sval = ii;
  return 0;
}

int32_t atoiz2(char* ss, int32_t* sval) {
  // version of atoiz which does not require the number to be null-terminated
  int32_t ii = atoi(ss);
  if ((ii < 1) && ((*ss != '0') || (ss[1] > ' '))) {
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
  } else if (ulii > 0xffffffffLLU) {
#else
  } else if (ulii == ULONG_MAX) {
    if (!memcmp(ss, "4294967295", 11)) {
      *valp = 0xffffffffU;
      return 0;
    }
#endif
    return 1;
  } else {
    *valp = ulii;
    return 0;
  }
}

uint32_t scan_two_doubles(char* ss, double* val1p, double* val2p) {
  char* ss2;
  *val1p = strtod(ss, &ss2);
  if (ss == ss2) {
    return 1;
  }
  ss = skip_initial_spaces(ss2);
  *val2p = strtod(ss, &ss2);
  return (ss == ss2)? 1 : 0;
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
    if (!is_set_ul(marker_exclude, *marker_uidx_ptr)) {
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

void get_top_two(uint32_t* uint_arr, uintptr_t uia_size, uintptr_t* top_idx_ptr, uintptr_t* second_idx_ptr) {
  uintptr_t cur_idx = 2;
  uintptr_t top_idx;
  uint32_t top_val;
  uintptr_t second_idx;
  uint32_t second_val;
  uintptr_t cur_val;
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

int32_t strcmp_se(char* s_read, const char* s_const, uint32_t len) {
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

uint32_t count_tokens(char* bufptr) {
  uint32_t token_ct = 0;
  bufptr = skip_initial_spaces(bufptr);
  while (!is_eoln_kns(*bufptr)) {
    token_ct++;
    bufptr = skip_initial_spaces(item_endnn(bufptr));
  }
  return token_ct;
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

char* int32_write(char* start, int32_t ii) {
  if (ii < 0) {
    if (ii < -2147483647) {
      return memcpya(start, "-2147483648", 11);
    }
    *start++ = '-';
    ii = -ii;
  }
  return uint32_write(start, (uint32_t)ii);
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
  if (llii <= 0xffffffffLL) {
    return uint32_write(start, (uint32_t)llii);
  }
  top_digits = llii / 100000000LL;
  bottom_eight = (uint32_t)(llii - (top_digits * 100000000));
  if (top_digits <= 0xffffffffLL) {
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

char* uint32_writew6(char* start, uint32_t uii) {
  // Minimum field width 6.
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      memset(start, 32, 5);
      start[5] = '0' + uii;
      return &(start[6]);
    } else if (uii < 100) {
      memset(start, 32, 4);
    } else {
      memset(start, 32, 3);
      quotient = uii / 100;
      start[3] = '0' + quotient;
      uii -= quotient * 100;
    }
    return memcpya(&(start[4]), &(digit2_table[uii * 2]), 2);
  } else if (uii < 10000000) {
    if (uii >= 100000) {
      if (uii >= 1000000) {
	quotient = uii / 1000000;
	*start++ = '0' + quotient;
	goto uint32_writew6_6b;
      }
      goto uint32_writew6_6;
    } else if (uii < 100000) {
      *start++ = ' ';
      quotient = uii / 10000;
      *start++ = '0' + quotient;
    } else {
      start = memseta(start, 32, 2);
      goto uint32_writew6_4;
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
  uint32_writew6_6b:
    uii -= 1000000 * quotient;
  uint32_writew6_6:
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 10000 * quotient;
 uint32_writew6_4:
  quotient = uii / 100;
  uii -= 100 * quotient;
  return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[uii * 2]), 2);
}

char* uint32_writew7(char* start, uint32_t uii) {
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

static inline char* uint32_write2trunc(char* start, uint32_t uii) {
  // Given 0 < uii < 100, writes uii without *trailing* zeroes.  (I.e. this is
  // for floating-point encoder use.)
  memcpy(start, &(digit2_table[uii * 2]), 2);
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uint32_write4trunc(char* start, uint32_t uii) {
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

static inline char* uint32_write8trunc(char* start, uint32_t uii) {
  uint32_t quotient = uii / 1000000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  uii -= 1000000 * quotient;
  if (uii) {
    quotient = uii / 10000;
    start += 2;
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
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uint32_write1p1(char* start, uint32_t quotient, uint32_t remainder) {
  start[0] = '0' + quotient;
  if (!remainder) {
    return &(start[1]);
  }
  start[1] = '.';
  start[2] = '0' + remainder;
  return &(start[3]);
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

static inline char* uint32_write1p7(char* start, uint32_t quotient, uint32_t remainder) {
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 100000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  remainder -= 100000 * quotient;
  if (remainder) {
    quotient = remainder / 1000;
    start += 2;
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

char* double_write2(char* start, double dxx) {
  // 2 sig fig number, 0.95 <= dxx < 99.5
  uint32_t quotient;
  if (dxx < 9.95) {
    dxx += 0.05;
    quotient = (int32_t)dxx;
    return uint32_write1p1(start, quotient, ((int32_t)(dxx * 10)) - (quotient * 10));
  }
  return memcpya(start, &(digit2_table[((uint32_t)((int32_t)(dxx + 0.5))) * 2]), 2);
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

char* double_write8(char* start, double dxx) {
  // 8 sig fig number, 0.99999995 <= dxx < 99999999.5
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.9999995) {
    if (dxx < 9.99999995) {
      dxx += 0.00000005;
      quotient = (int32_t)dxx;
      return uint32_write1p7(start, quotient, ((int32_t)(dxx * 10000000)) - (quotient * 10000000));
    }
    dxx += 0.0000005;
    uii = (int32_t)dxx;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    quotient = uii / 10000;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 10000 * quotient;
    if (uii) {
      start += 2;
    double_write8_pretail4:
      quotient = uii / 100;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      uii -= 100 * quotient;
      if (uii) {
	start += 2;
      double_write8_pretail2:
        memcpy(start, &(digit2_table[uii * 2]), 2);
      }
    }
  double_write8_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 9999.99995) {
    if (dxx < 999.999995) {
      dxx += 0.000005;
      uii = (int32_t)dxx;
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      uii = ((int32_t)(dxx * 100000)) - (uii * 100000);
      if (!uii) {
	return start;
      }
      *start++ = '.';
      quotient = uii / 1000;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      uii -= quotient * 1000;
      if (!uii) {
        goto double_write8_tail;
      }
      start += 2;
    double_write8_pretail3:
      quotient = uii / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      uii -= quotient * 10;
      if (!uii) {
	goto double_write8_tail;
      }
      start[2] = '0' + uii;
      return &(start[3]);
    }
    dxx += 0.00005;
    uii = (int32_t)dxx;
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uii = ((int32_t)(dxx * 10000)) - (uii * 10000);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    goto double_write8_pretail4;
  } else if (dxx < 999999.995) {
    if (dxx < 99999.9995) {
      dxx += 0.0005;
      uii = (int32_t)dxx;
      quotient = uii / 10000;
      *start = '0' + quotient;
      remainder = uii - 10000 * quotient;
      quotient = remainder / 100;
      start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
      remainder = remainder - 100 * quotient;
      start = memcpya(start, &(digit2_table[remainder * 2]), 2);
      uii = ((int32_t)(dxx * 1000)) - (uii * 1000);
      if (!uii) {
	return start;
      }
      *start++ = '.';
      goto double_write8_pretail3;
    }
    dxx += 0.005;
    uii = (int32_t)dxx;
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    remainder = uii - 10000 * quotient;
    quotient = remainder / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    remainder = remainder - 100 * quotient;
    start = memcpya(start, &(digit2_table[remainder * 2]), 2);
    uii = ((int32_t)(dxx * 100)) - (uii * 100);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    goto double_write8_pretail2;
  } else if (dxx < 9999999.95) {
    dxx += 0.05;
    uii = (int32_t)dxx;
    quotient = uii / 1000000;
    *start = '0' + quotient;
    remainder = uii - 1000000 * quotient;
    quotient = remainder / 10000;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    remainder = uii - 10000 * quotient;
    quotient = remainder / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    remainder = remainder - 100 * quotient;
    start = memcpya(start, &(digit2_table[remainder * 2]), 2);
    uii = ((int32_t)(dxx * 10)) - (uii * 10);
    if (!uii) {
      return start;
    }
    *start++ = '.';
    start[1] = '0' + uii;
    return &(start[2]);
  } else {
    uint32_write8(start, (int32_t)(dxx + 0.5));
    return &(start[8]);
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

char* double_f_writew3(char* start, double dxx) {
  int64_t llii;
  uint32_t uii;
  uint32_t quotient;
  if (dxx != dxx) {
    *((uint32_t*)start) = *((uint32_t*)"nan");
    return &(start[3]);
  } else if (dxx < 9.9995) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.9995) {
        goto double_f_writew3_10;
      }
    }
    dxx += 0.0005;
    uii = ((int32_t)dxx);
    *start++ = '0' + uii;
    uii = ((int32_t)(dxx * 1000)) - (uii * 1000);
  double_f_writew3_dec:
    *start++ = '.';
    quotient = uii / 100;
    uii -= 100 * quotient;
    *start++ = '0' + quotient;
    return memcpya(start, &(digit2_table[uii * 2]), 2);
  }
 double_f_writew3_10:
  dxx += 0.0005;
#ifndef __LP64__
  if (dxx <= 2147483.625) {
    uii = (int32_t)dxx;
    start = uint32_write(start, uii);
    uii = ((int32_t)(dxx * 1000)) - (uii * 1000);
    goto double_f_writew3_dec;
  }
#endif
  // 2 ^ 63 int64_t max, divided by 1000
  if (dxx <= 9223372036854775.75) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (uint32_t)(((int64_t)(dxx * 1000)) - (llii * 1000));
    goto double_f_writew3_dec;
  } else if (dxx < 9223372036854775808.0) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    // actually, since there are just 53 bits of precision, this should always
    // be 0, but whatever...
    uii = (int32_t)((dxx - ((double)llii)) * 1000);
    goto double_f_writew3_dec;
  }
  if (dxx == INFINITY) {
    *((uint32_t*)start) = *((uint32_t*)"inf");
    return &(start[3]);
  }
  // don't worry about optimizing %f on huge-ass finite numbers for now, since
  // it should be irrelevant for PLINK
  start += sprintf(start, "%.3f", dxx);
  return start;
}

char* double_f_writew96(char* start, double dxx) {
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
	goto double_f_writew96_10;
      }
    } else {
      *start++ = ' ';
    }
    dxx += 0.0000005;
    uii = ((int32_t)dxx);
    *start++ = '0' + uii;
    uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
  double_f_writew96_dec:
    *start++ = '.';
    uint32_write6(start, uii);
    return &(start[6]);
  }
 double_f_writew96_10:
  dxx += 0.0000005;
#ifndef __LP64__
  if (dxx < 2147.375) {
    uii = (int32_t)dxx;
    start = uint32_write(start, uii);
    uii = ((int32_t)(dxx * 1000000)) - (uii * 1000000);
    goto double_f_writew96_dec;
  }
#endif
  // 2 ^ 63 int64_t max, divided by 1000000
  if (dxx <= 9223372036854.75) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (uint32_t)(((int64_t)(dxx * 1000000)) - (llii * 1000000));
    goto double_f_writew96_dec;
  } else if (dxx < 9223372036854775808.0) {
    llii = (int64_t)dxx;
    start = int64_write(start, llii);
    uii = (int32_t)((dxx - ((double)llii)) * 1000000);
    goto double_f_writew96_dec;
  }
  if (dxx == INFINITY) {
    *((uint32_t*)start) = *((uint32_t*)"inf");
    return &(start[3]);
  }
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

char* double_g_writewx2(char* start, double dxx, uint32_t min_width) {
  // assumes min_width >= 5.
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
  if (dxx < 9.95e-5) {
    // 2 sig fig exponential notation, small
    if (dxx < 9.95e-16) {
      if (dxx < 9.95e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.95e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.95e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.95e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.95e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.95e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.95e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.95e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.95e-1) {
      dxx *= 10;
      xp10++;
    }
    dxx += 0.05;
    uii = (int32_t)dxx;
    wpos = uint32_write1p1(wpos, uii, ((int32_t)(dxx * 10)) - (uii * 10));
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
  } else if (dxx >= 99.5) {
    // 2 sig fig exponential notation, large
    if (dxx >= 9.95e15) {
      if (dxx >= 9.95e127) {
	if (dxx == INFINITY) {
	  start = memseta(start, 32, min_width - 4);
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.95e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.95e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.95e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.95e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.95e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.95e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.95e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.95e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    dxx += 0.05;
    uii = (int32_t)dxx;
    wpos = uint32_write1p1(wpos, uii, ((int32_t)(dxx * 10)) - (uii * 10));
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
    if (dxx >= 0.995) {
      wpos = double_write2(wpos, dxx);
    } else {
      // 2 sig fig decimal, no less than ~0.0001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.95e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.95e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uint32_write2trunc(wpos, (int32_t)((dxx * 100) + 0.5));
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

char* double_g_writewx4(char* start, double dxx, uint32_t min_width) {
  // assumes min_width >= 5.
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

char* double_g_writewx8(char* start, double dxx, uint32_t min_width) {
  // assumes min_width >= 8.
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
  if (dxx < 9.99999995e-5) {
    // 8 sig fig exponential notation, small
    if (dxx < 9.99999995e-16) {
      if (dxx < 9.99999995e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.99999995e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.99999995e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.99999995e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.99999995e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.99999995e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.99999995e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.99999995e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.99999995e-1) {
      dxx *= 10;
      xp10++;
    }
    dxx += 0.00000005;
    uii = (int32_t)dxx;
    wpos = uint32_write1p7(wpos, uii, ((int32_t)(dxx * 10000000)) - (uii * 10000000));
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
  } else if (dxx >= 99999999.5) {
    // 8 sig fig exponential notation, large
    if (dxx >= 9.99999995e15) {
      if (dxx >= 9.99999995e127) {
	if (dxx == INFINITY) {
	  start = memseta(start, 32, min_width - 4);
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.99999995e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.99999995e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.99999995e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.99999995e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.99999995e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.99999995e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.99999995e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.99999995e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    dxx += 0.00000005;
    uii = (int32_t)dxx;
    wpos = uint32_write1p7(wpos, uii, ((int32_t)(dxx * 10000000)) - (uii * 10000000));
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
    if (dxx >= 0.999999995) {
      wpos = double_write8(wpos, dxx);
    } else {
      // 8 sig fig decimal, no less than ~0.0001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.99999995e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.99999995e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uint32_write8trunc(wpos, (int32_t)((dxx * 100000000) + 0.5));
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
  } else if (num > 26) {
    // --allow-extra-chr 0
    *buf = '0';
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
    quotient = 0x80000000U / divisor;
    remainder = 0x80000000U - (quotient * divisor);
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

void fill_bits(uintptr_t* bit_arr, uintptr_t loc_start, uintptr_t len) {
  uintptr_t maj_start = loc_start / BITCT;
  uintptr_t maj_end = (loc_start + len) / BITCT;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bit_arr[maj_start] |= (ONELU << ((loc_start + len) % BITCT)) - (ONELU << (loc_start % BITCT));
  } else {
    bit_arr[maj_start] |= ~((ONELU << (loc_start % BITCT)) - ONELU);
    fill_ulong_one(&(bit_arr[maj_start + 1]), maj_end - maj_start - 1);
    minor = (loc_start + len) % BITCT;
    if (minor) {
      bit_arr[maj_end] |= (ONELU << minor) - ONELU;
    }
  }
}

void clear_bits(uintptr_t* bit_arr, uintptr_t loc_start, uintptr_t len) {
  uintptr_t maj_start = loc_start / BITCT;
  uintptr_t maj_end = (loc_start + len) / BITCT;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bit_arr[maj_start] &= ~((ONELU << ((loc_start + len) % BITCT)) - (ONELU << (loc_start % BITCT)));
  } else {
    bit_arr[maj_start] &= ((ONELU << (loc_start % BITCT)) - ONELU);
    fill_ulong_zero(&(bit_arr[maj_start + 1]), maj_end - maj_start - 1);
    minor = (loc_start + len) % BITCT;
    if (minor) {
      bit_arr[maj_end] &= ~((ONELU << minor) - ONELU);
    }
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

uint32_t next_unset_unsafe(uintptr_t* bit_arr, uint32_t loc) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (~(*bit_arr_ptr)) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bit_arr_ptr);
  } while (ulii == ~ZEROLU);
  return ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(~ulii);
}

#ifdef __LP64__
uintptr_t next_unset_ul_unsafe(uintptr_t* bit_arr, uintptr_t loc) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (~(*bit_arr_ptr)) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bit_arr_ptr);
  } while (ulii == ~ZEROLU);
  return (((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(~ulii));
}
#endif

uint32_t next_unset(uintptr_t* bit_arr, uint32_t loc, uint32_t ceil) {
  // safe version.  ceil >= 1.
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (~(*bit_arr_ptr)) >> (loc % BITCT);
  uintptr_t* bit_arr_last;
  if (ulii) {
    loc += CTZLU(ulii);
    return MINV(loc, ceil);
  }
  bit_arr_last = &(bit_arr[(ceil - 1) / BITCT]);
  do {
    if (bit_arr_ptr == bit_arr_last) {
      return ceil;
    }
    ulii = *(++bit_arr_ptr);
  } while (ulii == ~ZEROLU);
  loc = ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(~ulii);
  return MINV(loc, ceil);
}

#ifdef __LP64__
uintptr_t next_unset_ul(uintptr_t* bit_arr, uintptr_t loc, uintptr_t ceil) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (~(*bit_arr_ptr)) >> (loc % BITCT);
  uintptr_t* bit_arr_last;
  if (ulii) {
    ulii = loc + CTZLU(ulii);
    return MINV(ulii, ceil);
  }
  bit_arr_last = &(bit_arr[(ceil - 1) / BITCT]);
  do {
    if (bit_arr_ptr == bit_arr_last) {
      return ceil;
    }
    ulii = *(++bit_arr_ptr);
  } while (ulii == ~ZEROLU);
  ulii = ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(~ulii);
  return MINV(ulii, ceil);
}
#endif

uint32_t next_set_unsafe(uintptr_t* bit_arr, uint32_t loc) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (*bit_arr_ptr) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bit_arr_ptr);
  } while (!ulii);
  return ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(ulii);
}

#ifdef __LP64__
uintptr_t next_set_ul_unsafe(uintptr_t* bit_arr, uintptr_t loc) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (*bit_arr_ptr) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bit_arr_ptr);
  } while (!ulii);
  return ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(ulii);
}
#endif

uint32_t next_set(uintptr_t* bit_arr, uint32_t loc, uint32_t ceil) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (*bit_arr_ptr) >> (loc % BITCT);
  uintptr_t* bit_arr_last;
  uint32_t rval;
  if (ulii) {
    rval = loc + CTZLU(ulii);
    return MINV(rval, ceil);
  }
  bit_arr_last = &(bit_arr[(ceil - 1) / BITCT]);
  do {
    if (bit_arr_ptr == bit_arr_last) {
      return ceil;
    }
    ulii = *(++bit_arr_ptr);
  } while (!ulii);
  rval = ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(ulii);
  return MINV(rval, ceil);
}

#ifdef __LP64__
uintptr_t next_set_ul(uintptr_t* bit_arr, uintptr_t loc, uintptr_t ceil) {
  uintptr_t* bit_arr_ptr = &(bit_arr[loc / BITCT]);
  uintptr_t ulii = (*bit_arr_ptr) >> (loc % BITCT);
  uintptr_t* bit_arr_last;
  if (ulii) {
    ulii = loc + CTZLU(ulii);
    return MINV(ulii, ceil);
  }
  bit_arr_last = &(bit_arr[(ceil - 1) / BITCT]);
  do {
    if (bit_arr_ptr == bit_arr_last) {
      return ceil;
    }
    ulii = *(++bit_arr_ptr);
  } while (!ulii);
  ulii = ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + CTZLU(ulii);
  return MINV(ulii, ceil);
}
#endif

int32_t last_set_bit(uintptr_t* bit_arr, uint32_t word_ct) {
  uintptr_t* bit_arr_ptr = &(bit_arr[word_ct]);
  uintptr_t ulii;
  do {
    ulii = *(--bit_arr_ptr);
    if (ulii) {
      return ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + BITCT - 1 - CLZLU(ulii);
    }
  } while (bit_arr_ptr > bit_arr);
  return -1;
}

int32_t last_clear_bit(uintptr_t* bit_arr, uint32_t ceil) {
  uintptr_t* bit_arr_ptr = &(bit_arr[ceil / BITCT]);
  uintptr_t ulii;
  ceil = ceil % BITCT;
  if (ceil) {
    ulii = (~(*bit_arr_ptr)) & ((ONELU << ceil) - ONELU);
    if (ulii) {
      return (ceil | (BITCT - ONELU)) - CLZLU(ulii);
    }
  }
  while (bit_arr_ptr > bit_arr) {
    ulii = ~(*(--bit_arr_ptr));
    if (ulii) {
      return ((uintptr_t)(bit_arr_ptr - bit_arr)) * BITCT + BITCT - 1 - CLZLU(ulii);
    }
  }
  return -1;
}

void fill_idx_to_uidx(uintptr_t* exclude_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx) {
  uint32_t* idx_to_uidx_end = &(idx_to_uidx[item_ct]);
  uint32_t item_uidx = 0;
  uint32_t* cur_stop;
  do {
    item_uidx = next_unset_unsafe(exclude_arr, item_uidx);
    if (&(idx_to_uidx[unfiltered_item_ct - item_uidx]) < idx_to_uidx_end) {
      cur_stop = &(idx_to_uidx[next_set_unsafe(exclude_arr, item_uidx)]);
    } else {
      cur_stop = idx_to_uidx_end;
    }
    do {
      *idx_to_uidx++ = item_uidx++;
    } while (idx_to_uidx < cur_stop);
  } while (idx_to_uidx < idx_to_uidx_end);
}

void fill_uidx_to_idx(uintptr_t* exclude_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx) {
  uint32_t exclude_ct = unfiltered_item_ct - item_ct;
  uint32_t excluded_so_far = 0;
  uint32_t item_idx = 0;
  uint32_t* uidx_to_idx_ptr;
  uint32_t item_uidx;
  uint32_t item_idx_stop;
  do {
    item_uidx = next_unset_unsafe(exclude_arr, item_idx + excluded_so_far);
    uidx_to_idx_ptr = &(uidx_to_idx[item_uidx]);
    excluded_so_far = item_uidx - item_idx;
    if (exclude_ct > excluded_so_far) {
      item_idx_stop = next_set_unsafe(exclude_arr, item_uidx) - excluded_so_far;
    } else {
      item_idx_stop = item_ct;
    }
    do {
      *uidx_to_idx_ptr++ = item_idx++;
    } while (item_idx < item_idx_stop);
  } while (item_idx < item_ct);
}

void fill_uidx_to_idx_incl(uintptr_t* include_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx) {
  uint32_t exclude_ct = unfiltered_item_ct - item_ct;
  uint32_t excluded_so_far = 0;
  uint32_t item_idx = 0;
  uint32_t* uidx_to_idx_ptr;
  uint32_t item_uidx;
  uint32_t item_idx_stop;
  do {
    item_uidx = next_set_unsafe(include_arr, item_idx + excluded_so_far);
    uidx_to_idx_ptr = &(uidx_to_idx[item_uidx]);
    excluded_so_far = item_uidx - item_idx;
    if (exclude_ct > excluded_so_far) {
      item_idx_stop = next_unset_unsafe(include_arr, item_uidx) - excluded_so_far;
    } else {
      item_idx_stop = item_ct;
    }
    do {
      *uidx_to_idx_ptr++ = item_idx++;
    } while (item_idx < item_idx_stop);
  } while (item_idx < item_ct);
}

void fill_uidx_to_idx_incl(uintptr_t* include_arr, uint32_t item_ct, uint32_t* uidx_to_idx) {
  uintptr_t item_uidx = 0;
  uint32_t item_idx;
  for (item_idx = 0; item_idx < item_ct; item_idx++) {
    item_uidx = next_set_unsafe(include_arr, item_uidx);
    uidx_to_idx[item_uidx++] = item_idx;
  }
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

// global since species_str() may be called by functions which don't actually
// care about Chrom_info
const char* g_species_singular = NULL;
const char* g_species_plural = NULL;

char* chrom_name_write(char* buf, Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t zero_extra_chroms) {
  // assumes chrom_idx is valid
  if (chrom_idx <= chrom_info_ptr->max_code) {
    return uint32_write(buf, chrom_idx);
  } else if (zero_extra_chroms) {
    *buf++ = '0';
    return buf;
  } else {
    return strcpya(buf, chrom_info_ptr->nonstd_names[chrom_idx]);
  }
}

void forget_extra_chrom_names(Chrom_info* chrom_info_ptr) {
  uint32_t name_ct = chrom_info_ptr->name_ct;
  char** nonstd_names;
  uint32_t chrom_name_idx;
  // guard against init_species() not being called yet
  if (name_ct) {
    nonstd_names = &(chrom_info_ptr->nonstd_names[chrom_info_ptr->max_code + 1]);
    for (chrom_name_idx = 0; chrom_name_idx < name_ct; chrom_name_idx++) {
      free(nonstd_names[chrom_name_idx]);
      nonstd_names[chrom_name_idx] = NULL;
    }
    chrom_info_ptr->name_ct = 0;
  }
}

uint32_t haploid_chrom_present(Chrom_info* chrom_info_ptr) {
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  uint32_t uii;
  for (uii = 0; uii < CHROM_MASK_INITIAL_WORDS; uii++) {
    if (chrom_mask[uii] & haploid_mask[uii]) {
      return 1;
    }
  }
  return 0;
}

int32_t marker_code_raw(char* sptr) {
  // any character <= ' ' is considered a terminator
  int32_t ii;
  if (*sptr == 'c') {
    if ((sptr[1] == 'h') && (sptr[2] == 'r')) {
      sptr = &(sptr[3]);
    } else {
      return -1;
    }
  }
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

int32_t marker_code(Chrom_info* chrom_info_ptr, char* sptr) {
  // does not require string to be null-terminated, and does not perform
  // exhaustive error-checking
  int32_t ii = marker_code_raw(sptr);
  char** nonstd_names = chrom_info_ptr->nonstd_names;
  uint32_t uii;
  uint32_t ujj;
  uint32_t slen;
  uint32_t slen2;
  if (ii >= MAX_POSSIBLE_CHROM) {
    switch (ii) {
    case CHROM_X:
      ii = chrom_info_ptr->x_code;
      break;
    case CHROM_Y:
      ii = chrom_info_ptr->y_code;
      break;
    case CHROM_XY:
      ii = chrom_info_ptr->xy_code;
      break;
    case CHROM_MT:
      ii = chrom_info_ptr->mt_code;
    }
  } else if (ii == -1) {
    ujj = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
    slen = strlen_se(sptr);
    for (uii = chrom_info_ptr->max_code + 1; uii < ujj; uii++) {
      slen2 = strlen(nonstd_names[uii]);
      if ((slen == slen2) && (!memcmp(nonstd_names[uii], sptr, slen))) {
	return (int32_t)uii;
      }
    }
    return -1;
  } else if (((uint32_t)ii) > chrom_info_ptr->max_code) {
    return -1;
  }
  return ii;
}

int32_t marker_code2(Chrom_info* chrom_info_ptr, char* sptr, uint32_t slen) {
  // when the chromosome name doesn't end with a space
  char* s_end = &(sptr[slen]);
  char tmpc = *s_end;
  int32_t retval;
  *s_end = ' ';
  retval = marker_code(chrom_info_ptr, sptr);
  *s_end = tmpc;
  return retval;
}

uint32_t get_marker_chrom_fo_idx(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx) {
  uint32_t* marker_binsearch = chrom_info_ptr->chrom_file_order_marker_idx;
  uint32_t chrom_fo_min = 0;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t chrom_fo_cur;
  while (chrom_ct - chrom_fo_min > 1) {
    chrom_fo_cur = (chrom_ct + chrom_fo_min) / 2;
    if (marker_binsearch[chrom_fo_cur] > marker_uidx) {
      chrom_ct = chrom_fo_cur;
    } else {
      chrom_fo_min = chrom_fo_cur;
    }
  }
  return chrom_fo_min;
}

int32_t resolve_or_add_chrom_name(Chrom_info* chrom_info_ptr, char* bufptr, int32_t* chrom_idx_ptr) {
  uint32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  char** nonstd_names = chrom_info_ptr->nonstd_names;
  uint32_t slen = strlen_se(bufptr);
  Ll_str* name_stack_ptr = chrom_info_ptr->incl_excl_name_stack;
  uint32_t in_name_stack = 0;
  uint32_t chrom_idx;
  uint32_t slen2;
  for (chrom_idx = chrom_info_ptr->max_code + 1; chrom_idx < chrom_code_end; chrom_idx++) {
    slen2 = strlen(nonstd_names[chrom_idx]);
    if ((slen == slen2) && (!memcmp(bufptr, nonstd_names[chrom_idx], slen))) {
      *chrom_idx_ptr = (int32_t)chrom_idx;
      return 0;
    }
  }
  if (chrom_code_end == MAX_POSSIBLE_CHROM) {
    logprint("Error: Too many distinct nonstandard chromosome names.\n");
    return RET_INVALID_FORMAT;
  }
  nonstd_names[chrom_code_end] = (char*)malloc(slen + 1);
  if (!nonstd_names[chrom_code_end]) {
    return RET_NOMEM;
  }
  while (name_stack_ptr) {
    slen2 = strlen(name_stack_ptr->ss);
    if ((slen == slen2) && (!memcmp(bufptr, name_stack_ptr->ss, slen))) {
      in_name_stack = 1;
      break;
    }
    name_stack_ptr = name_stack_ptr->next;
  }
  if ((in_name_stack && chrom_info_ptr->is_include_stack) || ((!in_name_stack) && (!chrom_info_ptr->is_include_stack))) {
    SET_BIT(chrom_info_ptr->chrom_mask, chrom_code_end);
  }
  memcpy(nonstd_names[chrom_code_end], bufptr, slen);
  nonstd_names[chrom_code_end][slen] = '\0';
  *chrom_idx_ptr = (int32_t)chrom_code_end;
  chrom_info_ptr->name_ct += 1;
  return 0;
}

void refresh_chrom_info(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx, uint32_t allow_x_haploid, uint32_t is_all_nonmale, uint32_t* chrom_end_ptr, uint32_t* chrom_fo_idx_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* is_haploid_ptr) {
  int32_t chrom_idx;
  *chrom_end_ptr = chrom_info_ptr->chrom_file_order_marker_idx[(*chrom_fo_idx_ptr) + 1];
  while (marker_uidx >= (*chrom_end_ptr)) {
    *chrom_end_ptr = chrom_info_ptr->chrom_file_order_marker_idx[(++(*chrom_fo_idx_ptr)) + 1];
  }
  chrom_idx = chrom_info_ptr->chrom_file_order[*chrom_fo_idx_ptr];
  *is_x_ptr = allow_x_haploid && (chrom_idx == chrom_info_ptr->x_code);
  *is_y_ptr = (chrom_idx == chrom_info_ptr->y_code);
  *is_haploid_ptr = allow_x_haploid && is_set(chrom_info_ptr->haploid_mask, chrom_idx);
  if (is_all_nonmale) {
    *is_haploid_ptr = (*is_haploid_ptr) && (!(*is_x_ptr));
  }
}

int32_t single_chrom_start(Chrom_info* chrom_info_ptr, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude) {
  // Assumes there is at least one marker, and there are no split chromosomes.
  // Returns first marker_uidx in chromosome if there is only one, or -1 if
  // there's more than one chromosome.
  uint32_t first_marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  uint32_t last_marker_chrom = get_marker_chrom(chrom_info_ptr, last_clear_bit(marker_exclude, unfiltered_marker_ct));
  if (get_marker_chrom(chrom_info_ptr, first_marker_uidx) == last_marker_chrom) {
    return first_marker_uidx;
  }
  return -1;
}

int32_t strcmp_casted(const void* s1, const void* s2) {
  return strcmp((char*)s1, (char*)s2);
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

int32_t get_uidx_from_unsorted(char* idstr, uintptr_t* exclude_arr, uintptr_t id_ct, char* unsorted_ids, uintptr_t max_id_len) {
  uintptr_t id_uidx = 0;
  uintptr_t slen_p1 = strlen(idstr) + 1;
  uintptr_t id_idx;
  if (slen_p1 > max_id_len) {
    return -1;
  }
  for (id_idx = 0; id_idx < id_ct; id_idx++) {
    id_uidx = next_non_set_unsafe(exclude_arr, id_uidx);
    if (!memcmp(idstr, &(unsorted_ids[id_uidx * max_id_len]), slen_p1)) {
      return (int32_t)((uint32_t)id_uidx);
    }
    id_uidx++;
  }
  return -1;
}

char* scan_for_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len) {
  uintptr_t id_idx;
  id_ct--;
  for (id_idx = 0; id_idx < id_ct; id_idx++) {
    if (!strcmp(&(sorted_ids[id_idx * max_id_len]), &(sorted_ids[(id_idx + 1) * max_id_len]))) {
      return &(sorted_ids[id_idx * max_id_len]);
    }
  }
  return NULL;
}

int32_t is_missing_pheno(char* bufptr, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t affection_01) {
  if ((atoi(bufptr) == missing_pheno) && is_space_or_eoln(bufptr[missing_pheno_len])) {
    return 1;
  } else if ((!affection_01) && (*bufptr == '0') && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

int32_t eval_affection(char* bufptr, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t affection_01) {
  if (is_missing_pheno(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
    return 1;
  } else if (((*bufptr == '0') || (*bufptr == '1') || ((*bufptr == '2') && (!affection_01))) && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

uint32_t triangle_divide(int64_t cur_prod, int32_t modif) {
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

void parallel_bounds(uint32_t ct, int32_t start, uint32_t parallel_idx, uint32_t parallel_tot, int32_t* bound_start_ptr, int32_t* bound_end_ptr) {
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = ((int64_t)ct) * (ct + modif);
  *bound_start_ptr = triangle_divide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = triangle_divide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// set align to 1 for no alignment
void triangle_fill(uint32_t* target_arr, uint32_t ct, uint32_t pieces, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t align) {
  int32_t modif = 1 - start * 2;
  uint32_t cur_piece = 1;
  int64_t ct_tr;
  int64_t cur_prod;
  int32_t lbound;
  int32_t ubound;
  uint32_t uii;
  uint32_t align_m1;
  parallel_bounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = ((int64_t)lbound) * (lbound + modif);
  ct_tr = (((int64_t)ubound) * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide(cur_prod, modif);
    uii = (lbound - ((int32_t)start)) & align_m1;
    if ((uii) && (uii != align_m1)) {
      lbound = start + ((lbound - ((int32_t)start)) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (((uint32_t)lbound) > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

int32_t write_ids(char* outname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len) {
  FILE* outfile;
  uintptr_t ulii;
  if (fopen_checked(&outfile, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  for (ulii = 0; ulii < unfiltered_indiv_ct; ulii++) {
    if (!IS_SET(indiv_exclude, ulii) && (fprintf(outfile, "%s\n", &(person_ids[ulii * max_person_id_len])) < 0)) {
      return RET_WRITE_FAIL;
    }
  }
  if (fclose(outfile)) {
    return RET_WRITE_FAIL;
  }
  return 0;
}

int32_t distance_d_write_ids(char* outname, char* outname_end, uint32_t dist_calc_type, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len) {
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

int32_t distance_req(uint64_t calculation_type, char* read_dists_fname) {
  return ((calculation_type & CALC_DISTANCE) || ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!read_dists_fname) && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))));
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

int32_t double_cmp_decr(const void* aa, const void* bb) {
  double cc = *((const double*)aa) - *((const double*)bb);
  if (cc < 0.0) {
    return 1;
  } else if (cc > 0.0) {
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

int32_t char_cmp_deref(const void* aa, const void* bb) {
  return (int32_t)(**((const char**)aa) - **((const char**)bb));
}

int32_t intcmp(const void* aa, const void* bb) {
  return *((const int32_t*)aa) - *((const int32_t*)bb);
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
void qsort_ext2(char* main_arr, intptr_t arr_length, intptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, intptr_t secondary_item_len, char* proxy_arr, intptr_t proxy_len) {
  intptr_t lii;
  for (lii = 0; lii < arr_length; lii++) {
    *(char**)(&(proxy_arr[lii * proxy_len])) = &(main_arr[lii * item_length]);
    memcpy(&(proxy_arr[lii * proxy_len + sizeof(void*)]), &(secondary_arr[lii * secondary_item_len]), secondary_item_len);
  }
  qsort(proxy_arr, arr_length, proxy_len, comparator_deref);
  for (lii = 0; lii < arr_length; lii++) {
    memcpy(&(secondary_arr[lii * secondary_item_len]), &(proxy_arr[lii * proxy_len + sizeof(void*)]), secondary_item_len);
    memcpy(&(proxy_arr[lii * proxy_len]), *(char**)(&(proxy_arr[lii * proxy_len])), item_length);
  }
  for (lii = 0; lii < arr_length; lii++) {
    memcpy(&(main_arr[lii * item_length]), &(proxy_arr[lii * proxy_len]), item_length);
  }
}

// This actually tends to be faster than just sorting an array of indices,
// because of memory locality issues.
int32_t qsort_ext(char* main_arr, intptr_t arr_length, intptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, intptr_t secondary_item_len) {
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
  intptr_t proxy_len = secondary_item_len + sizeof(void*);
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

int32_t sort_item_ids_noalloc(char* sorted_ids, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*)) {
  // Stores a lexicographically sorted list of IDs in sorted_ids and the raw
  // positions of the corresponding markers/indivs in *id_map_ptr.  Does not
  // include excluded markers/indivs in the list.
  // Assumes sorted_ids and id_map have been allocated; use the sort_item_ids()
  // wrapper if they haven't been.
  // Note that this DOES still perform a "stack" allocation (in the qsort_ext()
  // call).
  uint32_t uii = 0;
  char* dup_id;
  char* tptr;
  uint32_t ujj;
  if (!item_ct) {
    return 0;
  }
  if (!collapse_idxs) {
    for (ujj = 0; ujj < item_ct; ujj++) {
      uii = next_non_set_unsafe(exclude_arr, uii);
      memcpy(&(sorted_ids[ujj * max_id_len]), &(item_ids[uii * max_id_len]), max_id_len);
      id_map[ujj] = uii++;
    }
  } else {
    for (ujj = 0; ujj < item_ct; ujj++) {
      uii = next_non_set_unsafe(exclude_arr, uii);
      memcpy(&(sorted_ids[ujj * max_id_len]), &(item_ids[uii * max_id_len]), max_id_len);
      id_map[ujj] = ujj;
      uii++;
    }
  }
  if (qsort_ext(sorted_ids, item_ct, max_id_len, comparator_deref, (char*)id_map, sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (!allow_dups) {
    dup_id = scan_for_duplicate_ids(sorted_ids, item_ct, max_id_len);
    if (dup_id) {
      tptr = strchr(dup_id, '\t');
      if (tptr) {
        *tptr = ' ';
      }
      sprintf(logbuf, "Error: Duplicate ID %s.\n", dup_id);
      logprintb();
      return RET_INVALID_FORMAT;
    }
  }
  return 0;
}

int32_t sort_item_ids(char** sorted_ids_ptr, uint32_t** id_map_ptr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t exclude_ct, char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*)) {
  uintptr_t item_ct = unfiltered_ct - exclude_ct;
  // id_map on bottom because --indiv-sort frees *sorted_ids_ptr
  if (wkspace_alloc_ui_checked(id_map_ptr, item_ct * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(sorted_ids_ptr, item_ct * max_id_len)) {
    return RET_NOMEM;
  }
  return sort_item_ids_noalloc(*sorted_ids_ptr, *id_map_ptr, unfiltered_ct, exclude_arr, item_ct, item_ids, max_id_len, allow_dups, collapse_idxs, comparator_deref);
}

uintptr_t uint64arr_greater_than(uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii) {
  // assumes arr_length is nonzero, and sorted_uint64_arr is in nondecreasing
  // order.
  // ullii guaranteed to be larger than sorted_uint64_arr[min_idx - 1] if it
  // exists, but NOT necessarily sorted_uint64_arr[min_idx].
  intptr_t min_idx = 0;
  // similarly, ullii guaranteed to be no greater than
  // sorted_uint64_arr[max_idx + 1] if it exists, but not necessarily
  // sorted_uint64_arr[max_idx].  Signed integer since it could become -1.
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (ullii > sorted_uint64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (ullii > sorted_uint64_arr[((uintptr_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

uintptr_t doublearr_greater_than(double* sorted_dbl_arr, uintptr_t arr_length, double dxx) {
  // returns number of items in sorted_dbl_arr which dxx is greater than.
  // assumes array is nonempty and sorted in nondecreasing order
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (dxx > sorted_dbl_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (dxx > sorted_dbl_arr[((uintptr_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

uintptr_t nonincr_doublearr_leq_stride(double* nonincr_dbl_arr, uintptr_t arr_length, uintptr_t stride, double dxx) {
  // assumes relevant elements of array are sorted in nonincreasing order
  // instead, and they are spaced stride units apart
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (dxx <= nonincr_dbl_arr[mid_idx * stride]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (dxx <= nonincr_dbl_arr[((uintptr_t)min_idx) * stride]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

void update_neighbor(uintptr_t indiv_ct, uint32_t neighbor_n2, uintptr_t indiv_idx1, uintptr_t indiv_idx2, double cur_ibs, double* neighbor_quantiles, uint32_t* neighbor_qindices) {
  uintptr_t exceed_ct;
  uintptr_t cur_write;
  exceed_ct = nonincr_doublearr_leq_stride(&(neighbor_quantiles[indiv_idx1]), neighbor_n2, indiv_ct, cur_ibs);
  if (exceed_ct < neighbor_n2) {
    for (cur_write = neighbor_n2 - 1; cur_write > exceed_ct; cur_write--) {
      neighbor_quantiles[cur_write * indiv_ct + indiv_idx1] = neighbor_quantiles[(cur_write - 1) * indiv_ct + indiv_idx1];
      neighbor_qindices[cur_write * indiv_ct + indiv_idx1] = neighbor_qindices[(cur_write - 1) * indiv_ct + indiv_idx1];
    }
    neighbor_quantiles[(exceed_ct * indiv_ct) + indiv_idx1] = cur_ibs;
    neighbor_qindices[(exceed_ct * indiv_ct) + indiv_idx1] = indiv_idx2;
  }
  exceed_ct = nonincr_doublearr_leq_stride(&(neighbor_quantiles[indiv_idx2]), neighbor_n2, indiv_ct, cur_ibs);
  if (exceed_ct < neighbor_n2) {
    for (cur_write = neighbor_n2 - 1; cur_write > exceed_ct; cur_write--) {
      neighbor_quantiles[cur_write * indiv_ct + indiv_idx2] = neighbor_quantiles[(cur_write - 1) * indiv_ct + indiv_idx2];
      neighbor_qindices[cur_write * indiv_ct + indiv_idx2] = neighbor_qindices[(cur_write - 1) * indiv_ct + indiv_idx2];
    }
    neighbor_quantiles[(exceed_ct * indiv_ct) + indiv_idx2] = cur_ibs;
    neighbor_qindices[(exceed_ct * indiv_ct) + indiv_idx2] = indiv_idx1;
  }
}

uintptr_t bsearch_str_lb(char* lptr, uintptr_t arr_length, uintptr_t max_id_len, char* id_buf) {
  // assumes nonempty array
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (strcmp(id_buf, &(lptr[mid_idx * max_id_len])) > 0) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (strcmp(id_buf, &(lptr[((uintptr_t)min_idx) * max_id_len])) > 0) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

int32_t bsearch_str(char* id_buf, char* lptr, uintptr_t max_id_len, intptr_t min_idx, intptr_t max_idx) {
  intptr_t mid_idx;
  int32_t ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp(id_buf, &(lptr[((uintptr_t)mid_idx) * max_id_len]));
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

int32_t bsearch_str_natural(char* id_buf, char* lptr, uintptr_t max_id_len, intptr_t min_idx, intptr_t max_idx) {
  intptr_t mid_idx;
  int32_t ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp_natural(id_buf, &(lptr[((uintptr_t)mid_idx) * max_id_len]));
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

int32_t bsearch_fam_indiv(char* id_buf, char* lptr, uintptr_t max_id_len, uint32_t filter_line_ct, char* fam_id, char* indiv_id) {
  // id_buf = workspace
  // lptr = packed, sorted list of ID strings to search over
  // fam_id and indiv_id are considered terminated by any space/eoln character
  uint32_t uii;
  uint32_t ujj;
  if (!filter_line_ct) {
    return -1;
  }
  uii = strlen_se(fam_id);
  ujj = strlen_se(indiv_id);
  if (uii + ujj + 2 > max_id_len) {
    return -1;
  }
  memcpyx(memcpyax(id_buf, fam_id, uii, '\t'), indiv_id, ujj, '\0');
  return bsearch_str(id_buf, lptr, max_id_len, 0, filter_line_ct - 1);
}

void bsearch_fam(char* id_buf, char* lptr, uintptr_t max_id_len, uint32_t filter_line_ct, char* fam_id, uint32_t* first_idx_ptr, uint32_t* last_idx_ptr) {
  uint32_t slen;
  uint32_t fidx;
  uint32_t loff;
  if (!filter_line_ct) {
    goto bsearch_fam_ret_null;
  }
  slen = strlen_se(fam_id);
  if (slen + 3 > max_id_len) {
    goto bsearch_fam_ret_null;
  }
  memcpy(id_buf, fam_id, slen);
  id_buf[slen] = '\t';
  id_buf[slen + 1] = '\0';
  fidx = bsearch_str_lb(lptr, filter_line_ct, max_id_len, id_buf);
  if (fidx == filter_line_ct) {
    goto bsearch_fam_ret_null;
  }
  id_buf[slen] = ' ';
  loff = bsearch_str_lb(&(lptr[fidx * max_id_len]), filter_line_ct - fidx, max_id_len, id_buf);
  if (!loff) {
    goto bsearch_fam_ret_null;
  }
  *first_idx_ptr = fidx;
  *last_idx_ptr = fidx + loff;
  return;
 bsearch_fam_ret_null:
  *first_idx_ptr = 0;
  *last_idx_ptr = 0;
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

void bitfield_and(uintptr_t* vv, uintptr_t* include_vec, uintptr_t word_ct) {
  // vv := vv AND include_vec
  // on 64-bit systems, assumes vv and include_vec are 16-byte aligned
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)vv;
  __m128i* iv128 = (__m128i*)include_vec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_and_si128(*iv128++, *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    vv[word_ct] &= include_vec[word_ct];
  }
#else
  uintptr_t* vec_end = &(vv[word_ct]);
  do {
    *vv++ &= *include_vec++;
  } while (vv < vec_end);
#endif
}

void bitfield_andnot(uintptr_t* vv, uintptr_t* exclude_vec, uintptr_t word_ct) {
  // vv := vv ANDNOT exclude_vec
  // on 64-bit systems, assumes vv and exclude_vec are 16-byte aligned
  // note that this is the reverse of the _mm_andnot() operand order
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)vv;
  __m128i* ev128 = (__m128i*)exclude_vec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_andnot_si128(*ev128++, *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    vv[word_ct] &= ~(exclude_vec[word_ct]);
  }
#else
  uintptr_t* vec_end = &(vv[word_ct]);
  do {
    *vv++ &= ~(*exclude_vec++);
  } while (vv < vec_end);
#endif
}

void bitfield_andnot_reversed_args(uintptr_t* vv, uintptr_t* include_vec, uintptr_t word_ct) {
  // vv := (~vv) AND include_vec
  // on 64-bit systems, assumes vv and exclude_vec are 16-byte aligned
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)vv;
  __m128i* iv128 = (__m128i*)include_vec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_andnot_si128(*vv128, *iv128++);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    vv[word_ct] = (~vv[word_ct]) & include_vec[word_ct];
  }
#else
  uintptr_t* vec_end = &(vv[word_ct]);
  do {
    *vv = (~(*vv)) & (*include_vec++);
    vv++;
  } while (vv < vec_end);
#endif
}

void bitfield_ornot(uintptr_t* vv, uintptr_t* inverted_or_vec, uintptr_t word_ct) {
  // vv := vv OR (~inverted_or_vec)
  // on 64-bit systems, assumes vv and inverted_or_vec are 16-byte aligned
#ifdef __LP64__
#ifdef _WIN32
  const __m128i all1 = {-1LL, -1LL};
#else
  const __m128i all1 = {0xffffffffffffffffLLU, 0xffffffffffffffffLLU};
#endif
  __m128i* vv128 = (__m128i*)vv;
  __m128i* ev128 = (__m128i*)inverted_or_vec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_or_si128(_mm_xor_si128(*ev128++, all1), *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    vv[word_ct] |= ~(inverted_or_vec[word_ct]);
  }
#else
  uintptr_t* vec_end = &(vv[word_ct]);
  do {
    *vv++ |= ~(*inverted_or_vec++);
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

uintptr_t popcount_bit_idx(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / BITCT;
  uintptr_t start_idxlr = start_idx & (BITCT - 1);
  uintptr_t end_idxl = end_idx / BITCT;
  uintptr_t end_idxlr = end_idx & (BITCT - 1);
  uintptr_t ct = 0;
  if (start_idxl == end_idxl) {
    return popcount_long(lptr[start_idxl] & ((ONELU << end_idxlr) - (ONELU << start_idxlr)));
  }
  if (start_idxlr) {
    ct = popcount_long(lptr[start_idxl++] >> start_idxlr);
  }
  if (end_idxl > start_idxl) {
    ct += popcount_longs(lptr, start_idxl, end_idxl);
  }
  if (end_idxlr) {
    ct += popcount_long(lptr[end_idxl] & ((ONELU << end_idxlr) - ONELU));
  }
  return ct;
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

void vertical_bitct_subtract(uintptr_t* bit_arr, uint32_t item_ct, uint32_t* sum_arr) {
  // assumes trailing bits are zeroed out
  uintptr_t cur_word;
  uint32_t idx_offset;
  uint32_t last_set_bit;
  for (idx_offset = 0; idx_offset < item_ct; idx_offset += BITCT) {
    cur_word = *bit_arr++;
    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      sum_arr[idx_offset + last_set_bit] -= 1;
      cur_word &= cur_word - ONELU;
    }
  }
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

uint32_t numeric_range_list_to_bitfield(Range_list* range_list_ptr, uint32_t item_ct, uintptr_t* bitfield, uint32_t offset, uint32_t ignore_overflow) {
  char* names = range_list_ptr->names;
  unsigned char* starts_range = range_list_ptr->starts_range;
  uint32_t name_ct = range_list_ptr->name_ct;
  uint32_t name_max_len = range_list_ptr->name_max_len;
  uint32_t idx_max = item_ct + offset;
  uint32_t name_idx;
  uint32_t idx1;
  uint32_t idx2;
  for (name_idx = 0; name_idx < name_ct; name_idx++) {
    idx1 = atoi(&(names[name_idx * name_max_len]));
    if (idx1 >= idx_max) {
      if (ignore_overflow) {
	continue;
      }
      return 1;
    }
    if (starts_range[name_idx]) {
      name_idx++;
      idx2 = atoi(&(names[name_idx * name_max_len]));
      if (idx2 >= idx_max) {
	if (!ignore_overflow) {
	  return 1;
	}
        idx2 = idx_max - 1;
      }
      fill_bits(bitfield, idx1 - offset, (idx2 - idx1) + 1);
    } else {
      set_bit(bitfield, idx1 - offset);
    }
  }
  return 0;
}

int32_t string_range_list_to_bitfield(char* header_line, uint32_t item_ct, Range_list* range_list_ptr, char* sorted_ids, uint32_t* id_map, int32_t* seen_idxs, const char* range_list_flag, const char* file_descrip, uintptr_t* bitfield) {
  uintptr_t max_id_len = range_list_ptr->name_max_len;
  uint32_t name_ct = range_list_ptr->name_ct;
  intptr_t name_ct_m1 = (uintptr_t)(name_ct - 1);
  uint32_t item_idx = 0;
  int32_t retval = 0;
  char* bufptr;
  char cc;
  int32_t ii;
  while (1) {
    bufptr = item_endnn(header_line);
    cc = *bufptr;
    *bufptr = '\0';
    ii = bsearch_str(header_line, sorted_ids, max_id_len, 0, name_ct_m1);
    *bufptr = cc;
    if (ii != -1) {
      if (seen_idxs[(uint32_t)ii] != -1) {
        sprintf(logbuf, "Error: Duplicate --%s token in %s.\n", range_list_flag, file_descrip);
        goto string_range_list_to_bitfield_ret_INVALID_FORMAT;
      }
      seen_idxs[(uint32_t)ii] = item_idx;
      if (ii && range_list_ptr->starts_range[(uint32_t)(ii - 1)]) {
        if (seen_idxs[ii - 1] == -1) {
          sprintf(logbuf, "Error: Second element of --%s range appears before first element in\n%s.\n", range_list_flag, file_descrip);
          goto string_range_list_to_bitfield_ret_INVALID_CMDLINE;
	}
	fill_bits(bitfield, seen_idxs[ii - 1], (item_idx - seen_idxs[ii - 1]) + 1);
      } else if (!(range_list_ptr->starts_range[(uint32_t)ii])) {
	SET_BIT(bitfield, item_idx);
      }
    }
    if (++item_idx == item_ct) {
      break;
    }
    header_line = skip_initial_spaces(&(bufptr[1]));
  }
  for (item_idx = 0; item_idx < name_ct; item_idx++) {
    if (seen_idxs[item_idx] == -1) {
      goto string_range_list_to_bitfield_ret_INVALID_CMDLINE_2;
    }
  }
  while (0) {
  string_range_list_to_bitfield_ret_INVALID_CMDLINE_2:
    sprintf(logbuf, "Error: Missing --%s token in %s.\n", range_list_flag, file_descrip);
  string_range_list_to_bitfield_ret_INVALID_CMDLINE:
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  string_range_list_to_bitfield_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

uint32_t count_non_autosomal_markers(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, uint32_t count_x) {
  // for backward compatibility, unplaced markers are considered to be
  // autosomal here
  uint32_t ct = 0;
  int32_t x_code = chrom_info_ptr->x_code;
  int32_t y_code = chrom_info_ptr->y_code;
  int32_t mt_code = chrom_info_ptr->mt_code;
  if (count_x && (x_code != -1)) {
    ct += count_chrom_markers(chrom_info_ptr, x_code, marker_exclude);
  }
  if (y_code != -1) {
    ct += count_chrom_markers(chrom_info_ptr, y_code, marker_exclude);
  }
  if (mt_code != -1) {
    ct += count_chrom_markers(chrom_info_ptr, mt_code, marker_exclude);
  }
  return ct;
}

uint32_t get_max_chrom_size(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, uint32_t* last_chrom_fo_idx_ptr) {
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t max_chrom_size = 0;
  uint32_t last_chrom_fo_idx = 0;
  uint32_t chrom_fo_idx;
  uint32_t cur_chrom_size;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    cur_chrom_size = count_chrom_markers(chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], marker_exclude);
    if (cur_chrom_size) {
      last_chrom_fo_idx = chrom_fo_idx;
      if (cur_chrom_size > max_chrom_size) {
        max_chrom_size = cur_chrom_size;
      }
    }
  }
  if (last_chrom_fo_idx_ptr) {
    *last_chrom_fo_idx_ptr = last_chrom_fo_idx;
  }
  return max_chrom_size;
}

void count_genders(uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uint32_t* male_ct_ptr, uint32_t* female_ct_ptr, uint32_t* unk_ct_ptr) {
  uint32_t male_ct = 0;
  uint32_t female_ct = 0;
  uint32_t unk_ct = 0;
  uint32_t unfiltered_indiv_ctld = unfiltered_indiv_ct / BITCT;
  uint32_t unfiltered_indiv_ct_rem = unfiltered_indiv_ct & (BITCT - 1);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t indiv_bidx;
  for (indiv_bidx = 0; indiv_bidx < unfiltered_indiv_ctld; indiv_bidx++) {
    ulii = ~(*indiv_exclude++);
  count_genders_last_loop:
    uljj = *sex_nm++;
    unk_ct += popcount_long(ulii & (~uljj));
    ulii &= uljj;
    uljj = *sex_male++;
    male_ct += popcount_long(ulii & uljj);
    female_ct += popcount_long(ulii & (~uljj));
  }
  if (unfiltered_indiv_ct_rem) {
    ulii = (~(*indiv_exclude)) & ((ONELU << unfiltered_indiv_ct_rem) - ONELU);
    unfiltered_indiv_ct_rem = 0;
    goto count_genders_last_loop;
  }
  *male_ct_ptr = male_ct;
  *female_ct_ptr = female_ct;
  *unk_ct_ptr = unk_ct;
}

uint32_t load_and_collapse(FILE* bedfile, uintptr_t* rawbuf, uint32_t unfiltered_indiv_ct, uintptr_t* mainbuf, uint32_t indiv_ct, uintptr_t* indiv_exclude) {
  uintptr_t cur_write = 0;
  uint32_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uint32_t indiv_uidx = 0;
  uint32_t indiv_idx_low = 0;
  uint32_t indiv_idx;
  if (unfiltered_indiv_ct == indiv_ct) {
    rawbuf = mainbuf;
  }
  if (fread(rawbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
    return RET_READ_FAIL;
  }
  if (unfiltered_indiv_ct == indiv_ct) {
    return 0;
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    cur_write |= ((rawbuf[indiv_uidx / BITCT2] >> (2 * (indiv_uidx % BITCT2))) & (3 * ONELU)) << (indiv_idx_low * 2);
    if (++indiv_idx_low == BITCT2) {
      *mainbuf++ = cur_write;
      cur_write = 0;
      indiv_idx_low = 0;
    }
    indiv_uidx++;
  }
  if (indiv_idx_low) {
    *mainbuf = cur_write;
  }
  return 0;
}

void collapse_copy_2bitarr_incl(uintptr_t* rawbuf, uintptr_t* mainbuf, uint32_t indiv_ct, uintptr_t* indiv_include) {
  uintptr_t cur_write = 0;
  uintptr_t indiv_uidx = 0;
  uint32_t indiv_idx_low = 0;
  uint32_t indiv_idx;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_set_unsafe(indiv_include, indiv_uidx);
    cur_write |= ((rawbuf[indiv_uidx / BITCT2] >> (2 * (indiv_uidx % BITCT2))) & (3 * ONELU)) << (indiv_idx_low * 2);
    if (++indiv_idx_low == BITCT2) {
      *mainbuf++ = cur_write;
      cur_write = 0;
      indiv_idx_low = 0;
    }
    indiv_uidx++;
  }
  if (indiv_idx_low) {
    *mainbuf = cur_write;
  }
}

uint32_t load_and_collapse_incl(FILE* bedfile, uintptr_t* rawbuf, uint32_t unfiltered_indiv_ct, uintptr_t* mainbuf, uint32_t indiv_ct, uintptr_t* indiv_include) {
  uint32_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  if (unfiltered_indiv_ct == indiv_ct) {
    rawbuf = mainbuf;
  }
  if (fread(rawbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
    return RET_READ_FAIL;
  }
  if (unfiltered_indiv_ct != indiv_ct) {
    collapse_copy_2bitarr_incl(rawbuf, mainbuf, indiv_ct, indiv_include);
  }
  return 0;
}

uint32_t block_load_autosomal(FILE* bedfile, int32_t bed_offset, uintptr_t* marker_exclude, uint32_t marker_ct_autosomal, uint32_t block_max_size, uintptr_t unfiltered_indiv_ct4, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_weights, unsigned char* readbuf, uint32_t* chrom_fo_idx_ptr, uintptr_t* marker_uidx_ptr, uintptr_t* marker_idx_ptr, uint32_t* block_size_ptr, double* set_allele_freq_buf, float* set_allele_freq_buf_fl, uint32_t* wtbuf) {
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uintptr_t marker_idx = *marker_idx_ptr;
  uint32_t chrom_fo_idx = *chrom_fo_idx_ptr;
  uint32_t chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
  uint32_t markers_read = 0;
  uint32_t autosome_ct = chrom_info_ptr->autosome_ct;
  uint32_t xy_code = (uint32_t)chrom_info_ptr->xy_code;
  uint32_t max_code = chrom_info_ptr->max_code;
  uint32_t cur_chrom;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;

  if (block_max_size > marker_ct_autosomal - marker_idx) {
    block_max_size = marker_ct_autosomal - marker_idx;
  }
  while (markers_read < block_max_size) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	return RET_READ_FAIL;
      }
    }
    if (marker_uidx >= chrom_end) {
      while (1) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	cur_chrom = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if ((cur_chrom <= autosome_ct) || (cur_chrom == xy_code) || (cur_chrom > max_code)) {
	  // for now, unplaced chromosomes are all "autosomal"
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

void vec_include_init(uintptr_t unfiltered_indiv_ct, uintptr_t* new_include2, uintptr_t* old_include) {
  uint32_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = ~(*old_include++);
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
    *new_include2++ = ulkk;
    *new_include2++ = ulmm;
  } while (--unfiltered_indiv_ctl);
  ulii = unfiltered_indiv_ct & (BITCT - 1);
  if (ulii) {
    new_include2--;
    if (ulii < BITCT2) {
      *new_include2-- = 0;
    } else {
      ulii -= BITCT2;
    }
    *new_include2 &= (ONELU << (ulii * 2)) - ONELU;
  }
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
  // Initializes a half-bitfield as the inverse of another.  Assumes target_arr
  // and source_arr are doubleword-aligned.
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

void hh_reset_y(unsigned char* loadbuf, uintptr_t* indiv_include2, uintptr_t* indiv_male_include2, uintptr_t unfiltered_indiv_ct) {
  uintptr_t indiv_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_indiv_ct + 3) / 4]);
  unsigned char* iicp;
  unsigned char* imicp;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char ucc3;
  uintptr_t unfiltered_indiv_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  uint32_t* indiv_include2_alias32;
  uint32_t* indiv_male_include2_alias32;
  __m128i* loadbuf_alias;
  __m128i* iivp;
  __m128i* imivp;
  __m128i vii;
  __m128i vjj;
  __m128i vkk;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    iivp = (__m128i*)indiv_include2;
    imivp = (__m128i*)indiv_male_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 64;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      // indiv_include2 & ~indiv_male_include2: force to 01
      // indiv_male_include2: convert 10 to 01, keep everything else
      vii = *imivp++;
      vjj = *iivp++;
      vkk = _mm_and_si128(*loadbuf_alias, _mm_or_si128(vii, _mm_slli_epi64(vii, 1)));
      *loadbuf_alias++ = _mm_or_si128(_mm_andnot_si128(vii, vjj), _mm_sub_epi64(vkk, _mm_and_si128(_mm_andnot_si128(vkk, _mm_srli_epi64(vkk, 1)), m1)));
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    iicp = (unsigned char*)iivp;
    imicp = (unsigned char*)imivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    indiv_include2_alias32 = (uint32_t*)indiv_include2;
    indiv_male_include2_alias32 = (uint32_t*)indiv_male_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 16;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *indiv_male_include2_alias32++;
      ujj = *indiv_include2_alias32++;
      ukk = (*loadbuf_alias32) & (uii * 3);
      *loadbuf_alias32++ = ((~uii) & ujj) | (ukk - ((~ukk) & (ukk >> 1) & 0x55555555));
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    iicp = (unsigned char*)indiv_include2_alias32;
    imicp = (unsigned char*)indiv_male_include2_alias32;
  } else {
    iicp = (unsigned char*)indiv_include2;
    imicp = (unsigned char*)indiv_male_include2;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 16;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *indiv_male_include2++;
      ujj = *indiv_include2++;
      ukk = (*loadbuf_alias32) & (uii * 3);
      *loadbuf_alias32++ = ((~uii) & ujj) | (ukk - ((~ukk) & (ukk >> 1) & 0x55555555));
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  iicp = (unsigned char*)indiv_include2;
  imicp = (unsigned char*)indiv_male_include2;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *imicp++;
    ucc2 = *iicp++;
    ucc3 = (*loadbuf) & (ucc * 3);
    *loadbuf++ = ((~ucc) & ucc2) | (ucc3 - ((~ucc3) & (ucc3 >> 1) & 0x55));
  }
}

/*
void force_unset_missing(unsigned char* loadbuf, uintptr_t* indiv_male_include2, uintptr_t unfiltered_indiv_ct) {
  uintptr_t indiv_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_indiv_ct + 3) / 4]);
  unsigned char* imicp;
  unsigned char ucc;
  unsigned char ucc3;
  uintptr_t unfiltered_indiv_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ukk;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  uint32_t* indiv_male_include2_alias32;
  __m128i* loadbuf_alias;
  __m128i* imivp;
  __m128i vii;
  __m128i vkk;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    imivp = (__m128i*)indiv_male_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 64;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      vii = *imivp++;
      vkk = _mm_and_si128(*loadbuf_alias, _mm_or_si128(vii, _mm_slli_epi64(vii, 1)));
      *loadbuf_alias++ = _mm_or_si128(_mm_andnot_si128(vii, m1), vkk);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    imicp = (unsigned char*)imivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    indiv_male_include2_alias32 = (uint32_t*)indiv_male_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 16;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *indiv_male_include2_alias32++;
      ukk = (*loadbuf_alias32) & (uii * 3);
      *loadbuf_alias32++ = ((~uii) & 0x55555555) | ukk;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    imicp = (unsigned char*)indiv_male_include2_alias32;
  } else {
    imicp = (unsigned char*)indiv_male_include2;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 16;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *indiv_male_include2++;
      ukk = (*loadbuf_alias32) & (uii * 3);
      *loadbuf_alias32++ = ((~uii) & 0x55555555) | ukk;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  imicp = (unsigned char*)indiv_male_include2;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *imicp++;
    ucc3 = (*loadbuf) & (ucc * 3);
    *loadbuf++ = ((~ucc) & 0x55) | ucc3;
  }
}
*/

void haploid_fix_multiple(uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, uint32_t xmhh_exists, uint32_t nxmhh_exists, uintptr_t* indiv_include2, uintptr_t* indiv_male_include2, uintptr_t unfiltered_indiv_ct, uintptr_t byte_ct_per_marker, unsigned char* loadbuf) {
  uintptr_t marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx_start);
  uint32_t chrom_fo_idx = get_marker_chrom_fo_idx(chrom_info_ptr, marker_uidx);
  uintptr_t marker_idx = 0;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uintptr_t chrom_end;
  uintptr_t marker_idx_chrom_end;

  while (marker_idx < marker_ct) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    is_x = (chrom_info_ptr->x_code == (int32_t)chrom_idx);
    is_y = (chrom_info_ptr->y_code == (int32_t)chrom_idx);
    is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
    marker_idx_chrom_end = marker_idx + chrom_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, chrom_end);
    if (marker_idx_chrom_end > marker_ct) {
      marker_idx_chrom_end = marker_ct;
    }
    if (is_haploid) {
      if (is_x) {
	if (xmhh_exists) {
	  for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	    hh_reset(&(loadbuf[marker_idx * byte_ct_per_marker]), indiv_male_include2, unfiltered_indiv_ct);
	  }
	}
      } else if (is_y) {
        for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
          hh_reset_y(&(loadbuf[marker_idx * byte_ct_per_marker]), indiv_include2, indiv_male_include2, unfiltered_indiv_ct);
	}
	/*
	for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	  force_unset_missing(&(loadbuf[marker_idx * byte_ct_per_marker]), indiv_male_include2, unfiltered_indiv_ct);
	}
	*/
      } else if (nxmhh_exists) {
	for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	  hh_reset(&(loadbuf[marker_idx * byte_ct_per_marker]), indiv_include2, unfiltered_indiv_ct);
	}
      }
    }
    marker_idx = marker_idx_chrom_end;
    chrom_fo_idx++;
  }
}


void reverse_loadbuf(unsigned char* loadbuf, uintptr_t unfiltered_indiv_ct) {
  uintptr_t indiv_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_indiv_ct + 3) / 4]);
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_indiv_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* loadbuf_alias;
  __m128i vii;
  __m128i vjj;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 64;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      vii = *loadbuf_alias;
      // we want to exchange 00 and 11, and leave 01/10 untouched.  So make
      // vjj := 11 iff vii is 00/11, and vjj := 00 otherwise; then xor.
      vjj = _mm_andnot_si128(_mm_xor_si128(vii, _mm_srli_epi64(vii, 1)), m1);
      vjj = _mm_or_si128(vjj, _mm_slli_epi64(vjj, 1));
      *loadbuf_alias++ = _mm_xor_si128(vii, vjj);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / BITCT2;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *loadbuf_alias32;
      ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
      ujj *= 3;
      *loadbuf_alias32++ = uii ^ ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / BITCT2;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *loadbuf_alias32;
      ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
      ujj *= 3;
      *loadbuf_alias32++ = uii ^ ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = 0x55 & (~(ucc ^ (ucc >> 1)));
    ucc2 *= 3;
    *loadbuf++ = ucc ^ ucc2;
  }
  uii = unfiltered_indiv_ct & 3;
  if (uii) {
    loadbuf[-1] &= (0xff >> (8 - 2 * uii));
  }
}

void force_missing(unsigned char* loadbuf, uintptr_t* force_missing_include2, uintptr_t unfiltered_indiv_ct) {
  uintptr_t indiv_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_indiv_ct + 3) / 4]);
  unsigned char* fmicp;
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_indiv_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  uint32_t* force_missing_include2_alias32;
  __m128i* loadbuf_alias;
  __m128i* fmivp;
  __m128i vii;
  __m128i vjj;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    fmivp = (__m128i*)force_missing_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / 64;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      vii = *loadbuf_alias;
      vjj = *fmivp++;
      vii = _mm_or_si128(vii, vjj);
      vjj = _mm_slli_epi64(vjj, 1);
      *loadbuf_alias++ = _mm_andnot_si128(vjj, vii);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    fmicp = (unsigned char*)fmivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    force_missing_include2_alias32 = (uint32_t*)force_missing_include2;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / BITCT2;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *loadbuf_alias32;
      ujj = *force_missing_include2_alias32++;
      uii |= ujj;
      ujj <<= 1;
      *loadbuf_alias32++ = uii & (~ujj);
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    fmicp = (unsigned char*)force_missing_include2_alias32;
  } else {
    fmicp = (unsigned char*)force_missing_include2;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_indiv_ctd = unfiltered_indiv_ct / BITCT2;
    for (; indiv_bidx < unfiltered_indiv_ctd; indiv_bidx++) {
      uii = *loadbuf_alias32;
      ujj = *force_missing_include2++;
      uii |= ujj;
      ujj <<= 1;
      *loadbuf_alias32++ = uii & (~ujj);
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  fmicp = (unsigned char*)force_missing_include2;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = *fmicp++;
    ucc |= ucc2;
    ucc2 <<= 1;
    *loadbuf++ = ucc & (~ucc2);
  }
}

int32_t open_and_size_string_list(char* fname, FILE** infile_ptr, uintptr_t* list_len_ptr, uintptr_t* max_str_len_ptr) {
  // assumes file is not open yet, and tbuf is safe to clobber
  uint32_t max_len = 0;
  uintptr_t list_len = 0;
  int32_t retval = 0;
  char* bufptr;
  uint32_t cur_len;
  if (fopen_checked(infile_ptr, fname, "r")) {
    goto open_and_size_string_list_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, *infile_ptr)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Pathologically long line in %s.\n", fname);
      logprintb();
      goto open_and_size_string_list_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    // don't complain about more than one entry on a line for now
    list_len++;
    cur_len = strlen_se(bufptr);
    if (cur_len >= max_len) {
      max_len = cur_len + 1;
    }
  }
  if (!feof(*infile_ptr)) {
    goto open_and_size_string_list_ret_READ_FAIL;
  }
  *list_len_ptr = list_len;
  *max_str_len_ptr = max_len;
  while (0) {
  open_and_size_string_list_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  open_and_size_string_list_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  open_and_size_string_list_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int32_t load_string_list(FILE** infile_ptr, uintptr_t max_str_len, char* str_list) {
  // assumes file is open (probably by open_and_size_string_list), and tbuf is
  // safe to clobber
  int32_t retval = 0;
  char* bufptr;
  uint32_t cur_len;
  rewind(*infile_ptr);
  while (fgets(tbuf, MAXLINELEN, *infile_ptr)) {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    cur_len = strlen_se(bufptr);
    memcpy(str_list, bufptr, cur_len);
    str_list[cur_len] = '\0';
    str_list = &(str_list[max_str_len]);
  }
  if (!feof(*infile_ptr)) {
    goto load_string_list_ret_READ_FAIL;
  }
  while (0) {
  load_string_list_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  return retval;
}

int32_t open_and_skip_first_lines(FILE** infile_ptr, char* fname, char* loadbuf, uintptr_t loadbuf_size, uint32_t lines_to_skip) {
  loadbuf[loadbuf_size - 1] = ' ';
  if (fopen_checked(infile_ptr, fname, "r")) {
    return RET_OPEN_FAIL;
  }
  while (lines_to_skip) {
    if (!fgets(loadbuf, loadbuf_size, *infile_ptr)) {
      if (feof(*infile_ptr)) {
	sprintf(logbuf, "Error: Fewer lines than expected in %s.\n", fname);
	logprintb();
	return RET_INVALID_FORMAT;
      } else {
	return RET_READ_FAIL;
      }
    }
    if (!(loadbuf[loadbuf_size - 1])) {
      return RET_NOMEM;
    }
    lines_to_skip--;
  }
  return 0;
}

int32_t load_to_first_token(FILE* infile, uintptr_t loadbuf_size, char comment_char, const char* file_descrip, char* loadbuf, char** bufptr_ptr) {
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      if ((loadbuf_size == MAXLINELEN) || (loadbuf_size == MAXLINEBUFLEN)) {
	sprintf(logbuf, "Error: Pathologically long line in %s.", file_descrip);
	logprintb();
	return RET_INVALID_FORMAT;
      } else {
	return RET_NOMEM;
      }
    }
    *bufptr_ptr = skip_initial_spaces(loadbuf);
    if (!is_eoln_kns(**bufptr_ptr)) {
      if ((**bufptr_ptr) != comment_char) {
        return 0;
      }
    }
  }
  if (!feof(infile)) {
    return RET_READ_FAIL;
  }
  sprintf(logbuf, "Error: Empty %s.", file_descrip);
  logprintb();
  return RET_INVALID_FORMAT;
}

int32_t open_and_load_to_first_token(FILE** infile_ptr, char* fname, uintptr_t loadbuf_size, char comment_char, const char* file_descrip, char* loadbuf, char** bufptr_ptr) {
  loadbuf[loadbuf_size - 1] = ' ';
  if (fopen_checked(infile_ptr, fname, "r")) {
    return RET_OPEN_FAIL;
  }
  return load_to_first_token(*infile_ptr, loadbuf_size, comment_char, file_descrip, loadbuf, bufptr_ptr);
}

int32_t scan_max_strlen(char* fname, uint32_t colnum, uint32_t colnum2, uint32_t headerskip, char skipchar, uintptr_t* max_str_len_ptr, uintptr_t* max_str2_len_ptr) {
  // colnum and colnum2 are 1-based indices.  If colnum2 is zero, only colnum
  // is scanned.
  // Includes terminating null in lengths.
  FILE* infile = NULL;
  char* loadbuf = (char*)wkspace_base;
  uintptr_t loadbuf_size = wkspace_left;
  uintptr_t max_str_len = *max_str_len_ptr;
  uintptr_t max_str2_len = 0;
  uint32_t colmin;
  uint32_t coldiff;
  char* str1_ptr;
  char* str2_ptr;
  char cc;
  uintptr_t cur_str_len;
  int32_t retval;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto scan_max_strlen_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, fname, loadbuf, loadbuf_size, headerskip);
  if (retval) {
    goto scan_max_strlen_ret_1;
  }
  if (colnum < colnum2) {
    max_str2_len = *max_str2_len_ptr;
    colmin = colnum - 1;
    coldiff = colnum2 - colnum;
  } else if (colnum2) {
    max_str2_len = max_str_len;
    max_str_len = *max_str2_len_ptr;
    colmin = colnum2 - 1;
    coldiff = colnum - colnum2;
  } else {
    colmin = colnum - 1;
    coldiff = 0;
    colnum2 = 0xffffffffU;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      goto scan_max_strlen_ret_NOMEM;
    }
    str1_ptr = skip_initial_spaces(loadbuf);
    cc = *str1_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colmin) {
      str1_ptr = next_item_mult(str1_ptr, colmin);
    }
    if (coldiff) {
      str2_ptr = next_item_mult(str1_ptr, coldiff);
    } else {
      str2_ptr = str1_ptr;
    }
    if (no_more_items_kns(str2_ptr)) {
      // probably want option for letting this slide in the future
      sprintf(logbuf, "Error: Fewer tokens than expected in %s line.\n", fname);
      goto scan_max_strlen_ret_INVALID_FORMAT;
    }
    cur_str_len = strlen_se(str1_ptr);
    if (cur_str_len >= max_str_len) {
      max_str_len = cur_str_len + 1;
    }
    if (coldiff) {
      cur_str_len = strlen_se(str2_ptr);
      if (cur_str_len >= max_str2_len) {
	max_str2_len = cur_str_len + 1;
      }
    }
  }
  if (!feof(infile)) {
    goto scan_max_strlen_ret_READ_FAIL;
  }
  if (colnum < colnum2) {
    *max_str_len_ptr = max_str_len;
    if (coldiff) {
      *max_str2_len_ptr = max_str2_len;
    }
  } else {
    *max_str_len_ptr = max_str2_len;
    *max_str2_len_ptr = max_str_len;
  }
  while (0) {
  scan_max_strlen_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  scan_max_strlen_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  scan_max_strlen_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 scan_max_strlen_ret_1:
  fclose_cond(infile);
  return retval;
}

int32_t scan_max_fam_indiv_strlen(char* fname, uint32_t colnum, uintptr_t* max_person_id_len_ptr) {
  // colnum is a 1-based index with the FID column number; IID column is
  // assumed to follow.
  // Includes terminating null in lengths.
  FILE* infile = NULL;
  char* loadbuf = (char*)wkspace_base;
  uintptr_t loadbuf_size = wkspace_left;
  uintptr_t max_person_id_len = *max_person_id_len_ptr;
  char* bufptr;
  char* bufptr2;
  uintptr_t cur_person_id_len;
  int32_t retval;
  colnum--;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto scan_max_fam_indiv_strlen_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, fname, loadbuf, loadbuf_size, 0);
  if (retval) {
    goto scan_max_fam_indiv_strlen_ret_1;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      goto scan_max_fam_indiv_strlen_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (colnum) {
      bufptr = next_item_mult(bufptr, colnum);
    }
    bufptr2 = next_item(bufptr);
    if (no_more_items_kns(bufptr2)) {
      sprintf(logbuf, "Error: Fewer tokens than expected in %s line.\n", fname);
      goto scan_max_fam_indiv_strlen_ret_INVALID_FORMAT;
    }
    cur_person_id_len = strlen_se(bufptr) + strlen_se(bufptr2) + 2;
    if (cur_person_id_len > max_person_id_len) {
      max_person_id_len = cur_person_id_len;
    }
  }
  if (!feof(infile)) {
    goto scan_max_fam_indiv_strlen_ret_READ_FAIL;
  }
  *max_person_id_len_ptr = max_person_id_len;
  while (0) {
  scan_max_fam_indiv_strlen_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  scan_max_fam_indiv_strlen_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  scan_max_fam_indiv_strlen_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 scan_max_fam_indiv_strlen_ret_1:
  fclose_cond(infile);
  return retval;
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
  // collapses array of fixed-length items, based on exclusion bitarray
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

void collapse_copy_bitarr(uint32_t orig_ct, uintptr_t* bit_arr, uintptr_t* exclude_arr, uint32_t filtered_ct, uintptr_t* output_arr) {
  uintptr_t ulii = 0;
  uint32_t item_uidx = 0;
  uint32_t write_bit = 0;
  uint32_t item_idx;
  fill_ulong_zero(output_arr, ((filtered_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
  for (item_idx = 0; item_idx < filtered_ct; item_idx++) {
    item_uidx = next_non_set_unsafe(exclude_arr, item_uidx);
    ulii |= ((bit_arr[item_uidx / BITCT] >> (item_uidx % BITCT)) & ONELU) << write_bit;
    if (++write_bit == BITCT) {
      *output_arr++ = ulii;
      ulii = 0;
      write_bit = 0;
    }
    item_uidx++;
  }
  if (filtered_ct % BITCT) {
    *output_arr = ulii;
  }
}

void collapse_copy_bitarr_incl(uint32_t orig_ct, uintptr_t* bit_arr, uintptr_t* include_arr, uint32_t filtered_ct, uintptr_t* output_arr) {
  uintptr_t ulii = 0;
  uintptr_t item_uidx = 0;
  uint32_t write_bit = 0;
  uint32_t item_idx;
  fill_ulong_zero(output_arr, ((filtered_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
  for (item_idx = 0; item_idx < filtered_ct; item_idx++) {
    item_uidx = next_set_unsafe(include_arr, item_uidx);
    ulii |= ((bit_arr[item_uidx / BITCT] >> (item_uidx % BITCT)) & ONELU) << write_bit;
    if (++write_bit == BITCT) {
      *output_arr++ = ulii;
      ulii = 0;
      write_bit = 0;
    }
    item_uidx++;
  }
  if (filtered_ct % BITCT) {
    *output_arr = ulii;
  }
}

void collapse_copy_bitarr_to_vec_incl(uint32_t orig_ct, uintptr_t* bit_arr, uintptr_t* include_arr, uint32_t filtered_ct, uintptr_t* output_vec) {
  uintptr_t ulii = 0;
  uintptr_t item_uidx = 0;
  uint32_t write_bit = 0;
  uint32_t item_idx;
  fill_ulong_zero(output_vec, 2 * ((filtered_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
  for (item_idx = 0; item_idx < filtered_ct; item_idx++) {
    item_uidx = next_set_unsafe(include_arr, item_uidx);
    ulii |= ((bit_arr[item_uidx / BITCT] >> (item_uidx % BITCT)) & ONELU) << (write_bit * 2);
    if (++write_bit == BITCT2) {
      *output_vec++ = ulii;
      ulii = 0;
      write_bit = 0;
    }
    item_uidx++;
  }
  if (filtered_ct % BITCT2) {
    *output_vec++ = ulii;
  }
  if ((filtered_ct + (BITCT2 - 1)) & BITCT2) {
    *output_vec = 0;
  }
}

uint32_t collapse_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len, uint32_t* id_starts) {
  // Collapses array of sorted IDs to remove duplicates, and writes
  // pre-collapse positions to id_starts (so e.g. duplication count of any
  // individual ID can be determined via subtraction).
  // Assumes id_ct is positive.
  uintptr_t read_idx;
  uintptr_t write_idx;
  id_starts[0] = 0;
  for (read_idx = 1; read_idx < id_ct; read_idx++) {
    if (!strcmp(&(sorted_ids[(read_idx - 1) * max_id_len]), &(sorted_ids[read_idx * max_id_len]))) {
      break;
    }
    id_starts[read_idx] = read_idx;
  }
  write_idx = read_idx;
  while (++read_idx < id_ct) {
    if (strcmp(&(sorted_ids[(write_idx - 1) * max_id_len]), &(sorted_ids[read_idx * max_id_len]))) {
      strcpy(&(sorted_ids[write_idx * max_id_len]), &(sorted_ids[read_idx * max_id_len]));
      id_starts[write_idx++] = read_idx;
    }
  }
  return write_idx;
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
  ukk = (uint32_t)(0x100000000LLU % ct);
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

void pick_d_small(unsigned char* tmp_cbuf, uint32_t* uibuf, uint32_t ct, uint32_t dd) {
  uint32_t uii;
  pick_d(tmp_cbuf, ct, dd);
  for (uii = 0; uii < ct; uii++) {
    if (tmp_cbuf[uii]) {
      *uibuf++ = uii;
    }
  }
  *uibuf = ct;
}

void init_sfmt64_from_sfmt32(sfmt_t* sfmt32, sfmt_t* sfmt64) {
  // sfmt_genrand_uint64() is not supposed to be called after
  // sfmt_genrand_uint32() is called on the same generator.  To work around
  // this, we initialize a new sfmt64 generator with this function when
  // necessary, and stick to genrand_uint32() calls with the main generator.
  uint32_t init_arr[4];
  uint32_t uii;
  for (uii = 0; uii < 4; uii++) {
    init_arr[uii] = sfmt_genrand_uint32(sfmt32);
  }
  sfmt_init_by_array(sfmt64, init_arr, 4);
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
  uint32_t tot_quotient = (uint32_t)(0x100000000LLU / tot_ct);
  uint32_t upper_bound = tot_ct * tot_quotient - 1;
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
	  } while (urand > upper_bound);
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
	  } while (urand > upper_bound);
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

void cluster_dist_divide(uintptr_t indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_starts, double* cluster_sdistances) {
  uintptr_t tcoord;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  double dxx;
  for (ulii = 0; ulii < cluster_ct; ulii++) {
    uii = cluster_starts[ulii + 1] - cluster_starts[ulii];
    if (uii > 1) {
      dxx = 1.0 / ((double)((int32_t)uii));
      uljj = (ulii * (ulii + 1)) / 2;
      for (tcoord = (ulii * (ulii - 1)) / 2; tcoord < uljj; tcoord++) {
	cluster_sdistances[tcoord] *= dxx;
      }
      for (uljj = ulii + 1; uljj < indiv_ct; uljj++) {
	cluster_sdistances[tri_coord_no_diag(ulii, uljj)] *= dxx;
      }
    }
  }
}

void cluster_dist_multiply(uintptr_t indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_starts, double* cluster_sdistances) {
  uintptr_t tcoord;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  double dxx;
  for (ulii = 0; ulii < cluster_ct; ulii++) {
    uii = cluster_starts[ulii + 1] - cluster_starts[ulii];
    if (uii > 1) {
      dxx = ((double)((int32_t)uii));
      uljj = (ulii * (ulii + 1)) / 2;
      for (tcoord = (ulii * (ulii - 1)) / 2; tcoord < uljj; tcoord++) {
	cluster_sdistances[tcoord] *= dxx;
      }
      for (uljj = ulii + 1; uljj < indiv_ct; uljj++) {
	cluster_sdistances[tri_coord_no_diag(ulii, uljj)] *= dxx;
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

// double regress_jack(uint32_t* uibuf) {
double regress_jack(uint32_t* uibuf, double* ret2_ptr) {
  uint32_t* uiptr = uibuf;
  uint32_t* uiptr2 = &(uibuf[g_jackknife_d]);
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
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
  while (uiptr < uiptr2) {
    dptr2 = &(g_jackknife_precomp[(*uiptr++) * JACKKNIFE_VALS_DIST]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  uiptr = uibuf;
  for (uii = 1; uii < g_jackknife_d; uii++) {
    ujj = *(++uiptr);
    dxx1 = g_pheno_d[ujj];
    uiptr2 = uibuf;
    dptr = &(g_dists[(((uintptr_t)ujj) * (ujj - 1)) / 2]);
    while (uiptr2 < uiptr) {
      ukk = *uiptr2++;
      dxx = (dxx1 + g_pheno_d[ukk]) * 0.5;
      dyy = dptr[ukk];
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
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t* uibuf = (uint32_t*)(&(g_generic_buf[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t))]));
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
    pick_d_small(cbuf, uibuf, g_indiv_ct, g_jackknife_d);
    dxx = regress_jack(uibuf, &ret2);
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
