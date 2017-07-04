// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_common.h"

#include <unistd.h> // sysconf()

#ifdef __APPLE__
  // needed for sysctl() call
  #include <sys/sysctl.h>
#endif

#include <time.h> // cleanup_logfile()

#ifdef __cplusplus
namespace plink2 {
#endif

const char g_errstr_fopen[] = "Error: Failed to open %s.\n";

char g_textbuf[kTextbufSize];

// now initialized by init_bigstack
const char* g_one_char_strs = nullptr;
// If one-base indels become sufficiently common, might want to predefine
// g_two_char_strs[], and update allele string construction/destruction
// accordingly.  (Though that should either be programmatically initialized, or
// only cover a subset of the space; 192k is a lot to increase the binary image
// size for a single simple table.)

const char* g_input_missing_geno_ptr = nullptr; // in addition to '.'
const char* g_output_missing_geno_ptr = nullptr; // now '.'

FILE* g_logfile = nullptr;

char g_logbuf[kMaxMediumLine * 2];

uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
uint32_t g_stderr_written_to = 0;

boolerr_t push_llstr(const char* ss, ll_str_t** ll_stack_ptr) {
  uintptr_t blen = strlen(ss) + 1;
  ll_str_t* new_llstr;
  if (pgl_malloc(sizeof(ll_str_t) + blen, &new_llstr)) {
    return 1;
  }
  new_llstr->next = *ll_stack_ptr;
  memcpy(new_llstr->ss, ss, blen);
  *ll_stack_ptr = new_llstr;
  return 0;
}

void logstr(const char* ss) {
  if (!g_debug_on) {
    fputs(ss, g_logfile);
    if (ferror(g_logfile)) {
      putchar('\n');
      fflush(stdout);
      fprintf(stderr, "Warning: Logging failure on:\n%s\nFurther logging will not be attempted in this run.\n", ss);
      g_log_failed = 1;
    }
  } else {
    if (g_log_failed) {
      fflush(stdout);
      fputs(ss, stderr);
    } else {
      fputs(ss, g_logfile);
      if (ferror(g_logfile)) {
	putchar('\n');
	fflush(stdout);
        fprintf(stderr, "Error: Debug logging failure.  Dumping to stderr:\n%s", ss);
	g_log_failed = 1;
	g_stderr_written_to = 1;
      } else {
	fflush(g_logfile);
      }
    }
  }
}

void logprint(const char* ss) {
  logstr(ss);
  fputs(ss, stdout);
}

void logerrprint(const char* ss) {
  logstr(ss);
  fflush(stdout);
  fputs(ss, stderr);
  g_stderr_written_to = 1;
}

void logprintb() {
  logstr(g_logbuf);
  fputs(g_logbuf, stdout);
}

void logerrprintb() {
  logstr(g_logbuf);
  fflush(stdout);
  fputs(g_logbuf, stderr);
  g_stderr_written_to = 1;
}

void wordwrap(uint32_t suffix_len, char* ss) {
  // Input: A null-terminated string with no intermediate newlines.  If
  //        suffix_len is zero, there should be a terminating \n; otherwise,
  //        the last character should be a space.  The allocation the string is
  //        part of must include at least ~80 bytes past the string end.
  // Effect: Spaces are replaced with newlines in a manner that plays well with
  //         80 column terminal windows.  (Multi-space blocks are never
  //         collapsed.)
  char* token_start = ss;
  char* line_end = &(ss[79]);
  char* token_end;
  while (1) {
    while (*token_start == ' ') {
      ++token_start;
    }
    if (token_start > line_end) {
      do {
	*line_end = '\n';
	line_end = &(line_end[80]);
      } while (token_start > line_end);
    }
    token_end = strchr(token_start, ' ');
    if (!token_end) {
      if (&(token_start[79]) == line_end) {
	return;
      }
      token_end = (char*)rawmemchr(token_start, '\0');
      if (!suffix_len) {
	if (token_end <= &(line_end[1])) {
	  // okay if end-of-string is one past the end, because function
	  // assumes last character is \n in suffix_len == 0 case
	  assert(token_end[-1] == '\n');
	  return;
	}
      } else {
        if (&(token_end[suffix_len]) <= line_end) {
	  return;
	}
	// because of terminal space assumption, token_start actually points
	// to the end of the string
	assert(token_start[-1] == ' ');
      }
      token_start[-1] = '\n';
      return;
    }
    if (token_end > line_end) {
      if (&(token_start[79]) != line_end) {
	token_start[-1] = '\n';
        line_end = &(token_start[79]);
	if (token_end > line_end) {
	  // single really long token, can't do anything beyond putting it on
	  // its own line
          *token_end = '\n';
	  line_end = &(token_end[80]);
	}
      } else {
	// single really long token, *and* previous token was either
	// nonexistent or long
	*token_end = '\n';
	line_end = &(token_end[80]);
      }
    }
    token_start = &(token_end[1]);
  }
}

void wordwrapb(uint32_t suffix_len) {
  wordwrap(suffix_len, g_logbuf);
}


boolerr_t fopen_checked(const char* fname, const char* mode, FILE** target_ptr) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    logprint("\n");
    LOGERRPRINTFWW(g_errstr_fopen, fname);
    return 1;
  }
  return 0;
}

interr_t fwrite_flush2(char* buf_flush, FILE* outfile, char** write_iter_ptr) {
  char* buf = &(buf_flush[-((int32_t)kMaxMediumLine)]);
  char* buf_end = *write_iter_ptr;
  *write_iter_ptr = buf;
  return fwrite_checked(buf, (uintptr_t)(buf_end - buf), outfile);
}


uint32_t int_slen(int32_t num) {
  int32_t slen = 1;
  uint32_t absnum;
  if (num < 0) {
    absnum = -num;
    ++slen;
  } else {
    absnum = num;
  }
  while (absnum > 99) {
    // division by a constant is faster for unsigned ints
    absnum /= 100;
    slen += 2;
  }
  if (absnum > 9) {
    ++slen;
  }
  return slen;
}

int32_t strcmp_se(const char* s_read, const char* s_const, uint32_t s_const_len) {
  return memcmp(s_read, s_const, s_const_len) || (!is_space_or_eoln(s_read[s_const_len]));
}

int32_t strcmp_casted(const void* s1, const void* s2) {
  return strcmp((const char*)s1, (const char*)s2);
}

// PLINK 2's natural sort uses the following logic:
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
int32_t strcmp_natural_scan_forward(const unsigned char* s1, const unsigned char* s2) {
  // assumes s1 and s2 currently point to the middle of a mismatching number,
  // where s1 < s2.
  unsigned char c1;
  unsigned char c2;
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
int32_t strcmp_natural_tiebroken(const unsigned char* s1, const unsigned char* s2) {
  // assumes ties should be broken in favor of s2.
  unsigned char c1 = *(++s1);
  unsigned char c2 = *(++s2);
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
      }
      if (c1 > c2) {
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
	}
	return -strcmp_natural_scan_forward(s2, s1);
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

static inline int32_t strcmp_natural_uncasted(const unsigned char* s1, const unsigned char* s2) {
  unsigned char c1 = *s1;
  unsigned char c2 = *s2;
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 0
  strcmp_natural_uncasted_state_0:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	if (c2 + 32 == c1) {
	  return -strcmp_natural_tiebroken(s2, s1);
	}
	if ((c2 < 'a') || (c2 > 'z')) {
	  c1 -= 32;
	}
      } else if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
	if (c1 == c2) {
	  return strcmp_natural_tiebroken(s1, s2);
	}
      }
      return (c1 < c2)? -1 : 1;
    }
    if (!c1) {
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
	}
	return -strcmp_natural_scan_forward(s2, s1);
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
  return strcmp_natural_uncasted((const unsigned char*)s1, (const unsigned char*)s2);
}

int32_t strcmp_deref(const void* s1, const void* s2) {
  // const_cast
  return strcmp(*(char**)((uintptr_t)s1), *(char**)((uintptr_t)s2));
}

int32_t strcmp_natural_deref(const void* s1, const void* s2) {
  // const_cast
  return strcmp_natural_uncasted(*(unsigned char**)((uintptr_t)s1), *(unsigned char**)((uintptr_t)s2));
}

int32_t float_cmp(const void* aa, const void* bb) {
  const float fxx = *((const float*)aa);
  const float fyy = *((const float*)bb);
  if (fxx < fyy) {
    return -1;
  }
  return (fxx > fyy);
}

int32_t double_cmp(const void* aa, const void* bb) {
  const double dxx = *((const double*)aa);
  const double dyy = *((const double*)bb);
  if (dxx < dyy) {
    return -1;
  }
  return (dxx > dyy);
}

int32_t double_cmp_decr(const void* aa, const void* bb) {
  const double dxx = *((const double*)aa);
  const double dyy = *((const double*)bb);
  if (dxx > dyy) {
    return -1;
  }
  return (dxx < dyy);
}

int32_t intcmp(const void* aa, const void* bb) {
  return *((const int32_t*)aa) - *((const int32_t*)bb);
}

int32_t uint64cmp(const void* aa, const void* bb) {
  const uint64_t ullaa = *((const uint64_t*)aa);
  const uint64_t ullbb = *((const uint64_t*)bb);
  if (ullaa < ullbb) {
    return -1;
  }
  return (ullaa > ullbb);
}

#ifndef __cplusplus
int32_t uint64cmp_decr(const void* aa, const void* bb) {
  const uint64_t ullaa = *((const uint64_t*)aa);
  const uint64_t ullbb = *((const uint64_t*)bb);
  if (ullaa > ullbb) {
    return -1;
  }
  return (ullaa < ullbb);
}
#endif

#ifdef __cplusplus
float destructive_get_fmedian(uintptr_t len, float* unsorted_arr) {
  if (!len) {
    return 0.0;
  }
  const uintptr_t len_d2 = len / 2;
  std::nth_element(unsorted_arr, &(unsorted_arr[len_d2]), &(unsorted_arr[len]));
  const float median_upper = unsorted_arr[len_d2];
  if (len % 2) {
    return median_upper;
  }
  return (get_fmax(len_d2, unsorted_arr) + median_upper) * 0.5f;
}

double destructive_get_dmedian(uintptr_t len, double* unsorted_arr) {
  if (!len) {
    return 0.0;
  }
  const uintptr_t len_d2 = len / 2;
  std::nth_element(unsorted_arr, &(unsorted_arr[len_d2]), &(unsorted_arr[len]));
  const double median_upper = unsorted_arr[len_d2];
  if (len % 2) {
    return median_upper;
  }
  return (get_dmax(len_d2, unsorted_arr) + median_upper) * 0.5;
}
#else
// these will probably be used in __cplusplus case too
float get_fmedian(const float* sorted_arr, uintptr_t len) {
  if (!len) {
    return 0.0f;
  }
  if (len % 2) {
    return sorted_arr[len / 2];
  }
  return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5f;
}

double get_dmedian(const double* sorted_arr, uintptr_t len) {
  if (!len) {
    return 0.0;
  }
  if (len % 2) {
    return sorted_arr[len / 2];
  }
  return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5;
}

float destructive_get_fmedian(uintptr_t len, float* unsorted_arr) {
  // no, I'm not gonna bother reimplementing introselect just for folks who
  // insist on compiling this as pure C instead of C++
  qsort(unsorted_arr, len, sizeof(float), float_cmp);
  return get_fmedian(unsorted_arr, len);
}

double destructive_get_dmedian(uintptr_t len, double* unsorted_arr) {
  qsort(unsorted_arr, len, sizeof(double), double_cmp);
  return get_dmedian(unsorted_arr, len);
}
#endif

// alas, qsort_r not available on some Linux distributions

// note that this can be expected to have size 16 bytes, not 12, on 64-bit
// systems
typedef struct str_sort_indexed_deref_struct {
  const char* strptr;
  uint32_t orig_idx;
#ifdef __cplusplus
  bool operator<(const struct str_sort_indexed_deref_struct& rhs) const {
    return (strcmp(strptr, rhs.strptr) < 0);
  }
#endif
} str_sort_indexed_deref_t;

#ifdef __cplusplus
typedef struct strbuf36_ui_struct {
  char strbuf[36];
  uint32_t orig_idx;
  bool operator<(const struct strbuf36_ui_struct& rhs) const {
    return (strcmp_natural_uncasted((const unsigned char*)strbuf, (const unsigned char*)(rhs.strbuf)) < 0);
  }
} Strbuf36_ui;

typedef struct strbuf60_ui_struct {
  char strbuf[60];
  uint32_t orig_idx;
  bool operator<(const struct strbuf60_ui_struct& rhs) const {
    return (strcmp_natural_uncasted((const unsigned char*)strbuf, (const unsigned char*)(rhs.strbuf)) < 0);
  }
} Strbuf60_ui;

static_assert(sizeof(Strbuf36_ui) == 40, "Strbuf36_ui is not laid out as expected.");
static_assert(offsetof(Strbuf36_ui, orig_idx) == 36, "Strbuf36_ui is not laid out as expected.");
static_assert(sizeof(Strbuf60_ui) == 64, "Strbuf60_ui is not laid out as expected.");
static_assert(offsetof(Strbuf60_ui, orig_idx) == 60, "Strbuf60_ui is not laid out as expected.");
uintptr_t get_strboxsort_wentry_blen(uintptr_t max_str_blen) {
  if (max_str_blen <= 36) {
    return sizeof(Strbuf36_ui);
  }
  if (max_str_blen <= 60) {
    return sizeof(Strbuf60_ui);
  }
  return max_str_blen;
}

typedef struct str_nsort_indexed_deref_struct {
  char* strptr;
  uint32_t orig_idx;
  bool operator<(const struct str_nsort_indexed_deref_struct& rhs) const {
    return (strcmp_natural_uncasted((unsigned char*)strptr, (unsigned char*)(rhs.strptr)) < 0);
  }
} str_nsort_indexed_deref_t;
#else
uintptr_t get_strboxsort_wentry_blen(uintptr_t max_str_blen) {
  return MAXV(max_str_blen, sizeof(str_sort_indexed_deref_t));
}
#endif

// assumed that sort_wkspace has size >= str_ct *
// max(sizeof(str_sort_indexed_deref_t), max_str_blen)
void sort_strbox_indexed2_fallback(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace) {
  str_sort_indexed_deref_t* wkspace_alias = (str_sort_indexed_deref_t*)sort_wkspace;
  for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    wkspace_alias[str_idx].strptr = &(strbox[str_idx * max_str_blen]);
    wkspace_alias[str_idx].orig_idx = id_map[str_idx];
  }
  if (!use_nsort) {
#ifdef __cplusplus
    std::sort(wkspace_alias, &(wkspace_alias[str_ct]));
#else
    qsort(wkspace_alias, str_ct, sizeof(str_sort_indexed_deref_t), strcmp_deref);
#endif
  } else {
#ifdef __cplusplus
    str_nsort_indexed_deref_t* wkspace_alias2 = (str_nsort_indexed_deref_t*)wkspace_alias;
    std::sort(wkspace_alias2, &(wkspace_alias2[str_ct]));
#else
    qsort(wkspace_alias, str_ct, sizeof(str_sort_indexed_deref_t), strcmp_natural_deref);
#endif
  }
  for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    id_map[str_idx] = wkspace_alias[str_idx].orig_idx;
  }
#ifndef __cplusplus
  if (max_str_blen < sizeof(str_sort_indexed_deref_t)) {
    // actually better to use non-deref sort here, but just get this working
    // properly for now
    for (uint32_t new_idx = 0; new_idx < str_ct; ++new_idx) {
      const char* strptr = wkspace_alias[new_idx].strptr;
      strcpy(&(((char*)wkspace_alias)[new_idx * max_str_blen]), strptr);
    }
  } else {
#endif
    // bugfix: need to handle id_map[str_idx] != str_idx
    uint32_t new_idx = str_ct;
    do {
      --new_idx;
      const char* strptr = wkspace_alias[new_idx].strptr;
      strcpy(&(((char*)wkspace_alias)[new_idx * max_str_blen]), strptr);
    } while (new_idx);
#ifndef __cplusplus
  }
#endif
  memcpy(strbox, wkspace_alias, str_ct * max_str_blen);
}

#ifdef __cplusplus
typedef struct word_cmp40b_struct {
  uintptr_t words[40 / kBytesPerWord];
  bool operator<(const struct word_cmp40b_struct& rhs) const {
    uint32_t idx = 0;
    do {
      const uintptr_t cur_word = words[idx];
      const uintptr_t rhs_word = rhs.words[idx];
      if (cur_word != rhs_word) {
	// could pre-reverse the strings?
        const uintptr_t xor_word = cur_word ^ rhs_word;
	const uint32_t lshift = (kBitsPerWord - 8) - (CTZLU(xor_word) & (kBitsPerWord - 8));
	return (cur_word << lshift) < (rhs_word << lshift);
      }
    } while (++idx < (40 / kBytesPerWord));
    return false;
  }
} word_cmp40b_t;

typedef struct word_cmp64b_struct {
  uintptr_t words[64 / kBytesPerWord];
  bool operator<(const struct word_cmp64b_struct& rhs) const {
    uint32_t idx = 0;
    do {
      const uintptr_t cur_word = words[idx];
      const uintptr_t rhs_word = rhs.words[idx];
      if (cur_word != rhs_word) {
        const uintptr_t xor_word = cur_word ^ rhs_word;
	const uint32_t lshift = (kBitsPerWord - 8) - (CTZLU(xor_word) & (kBitsPerWord - 8));
	return (cur_word << lshift) < (rhs_word << lshift);
      }
    } while (++idx < (64 / kBytesPerWord));
    return false;
  }
} word_cmp64b_t;

static_assert(sizeof(word_cmp40b_t) == 40, "word_cmp40b_t does not have the expected size.");
static_assert(sizeof(word_cmp64b_t) == 64, "word_cmp64b_t does not have the expected size.");

void sort_strbox_40b_finish(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, Strbuf36_ui* filled_wkspace, char* sorted_strbox, uint32_t* id_map) {
  if (!use_nsort) {
    word_cmp40b_t* wkspace_alias = (word_cmp40b_t*)filled_wkspace;
    std::sort(wkspace_alias, &(wkspace_alias[str_ct]));
  } else {
    std::sort(filled_wkspace, &(filled_wkspace[str_ct]));
  }
  for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    strcpy(&(sorted_strbox[str_idx * max_str_blen]), filled_wkspace[str_idx].strbuf);
    id_map[str_idx] = filled_wkspace[str_idx].orig_idx;
  }
}

void sort_strbox_64b_finish(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, Strbuf60_ui* filled_wkspace, char* sorted_strbox, uint32_t* id_map) {
  if (!use_nsort) {
    word_cmp64b_t* wkspace_alias = (word_cmp64b_t*)filled_wkspace;
    std::sort(wkspace_alias, &(wkspace_alias[str_ct]));
  } else {
    std::sort(filled_wkspace, &(filled_wkspace[str_ct]));
  }
  for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    strcpy(&(sorted_strbox[str_idx * max_str_blen]), filled_wkspace[str_idx].strbuf);
    id_map[str_idx] = filled_wkspace[str_idx].orig_idx;
  }
}

// Normally use sort_strbox_indexed(), but this version is necessary before
// g_bigstack has been allocated.
void sort_strbox_indexed2(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace) {
  if (max_str_blen <= 36) {
    Strbuf36_ui* wkspace_alias = (Strbuf36_ui*)sort_wkspace;
    for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
      const char* cur_str = &(strbox[str_idx * max_str_blen]);
      strcpy(wkspace_alias[str_idx].strbuf, cur_str);
      wkspace_alias[str_idx].orig_idx = id_map[str_idx];
    }
    sort_strbox_40b_finish(str_ct, max_str_blen, use_nsort, wkspace_alias, strbox, id_map);
    return;
  }
  if (max_str_blen <= 60) {
    Strbuf60_ui* wkspace_alias = (Strbuf60_ui*)sort_wkspace;
    for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
      const char* cur_str = &(strbox[str_idx * max_str_blen]);
      strcpy(wkspace_alias[str_idx].strbuf, cur_str);
      wkspace_alias[str_idx].orig_idx = id_map[str_idx];
    }
    sort_strbox_64b_finish(str_ct, max_str_blen, use_nsort, wkspace_alias, strbox, id_map);
    return;
  }
  sort_strbox_indexed2_fallback(str_ct, max_str_blen, use_nsort, strbox, id_map, sort_wkspace);
}
#endif

boolerr_t sort_strbox_indexed(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map) {
  if (str_ct < 2) {
    return 0;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  const uintptr_t wkspace_entry_blen = get_strboxsort_wentry_blen(max_str_blen);
  unsigned char* sort_wkspace;
  if (bigstack_alloc_uc(str_ct * wkspace_entry_blen, &sort_wkspace)) {
    return 1;
  }
  sort_strbox_indexed2(str_ct, max_str_blen, use_nsort, strbox, id_map, sort_wkspace);
  bigstack_reset(bigstack_mark);
  return 0;
}

boolerr_t sort_strbox_indexed_malloc(uintptr_t str_ct, uintptr_t max_str_blen, char* strbox, uint32_t* id_map) {
  if (str_ct < 2) {
    return 0;
  }
  const uintptr_t wkspace_entry_blen = get_strboxsort_wentry_blen(max_str_blen);
  unsigned char* sort_wkspace;
  if (pgl_malloc(str_ct * wkspace_entry_blen, &sort_wkspace)) {
    return 1;
  }
  sort_strbox_indexed2(str_ct, max_str_blen, 0, strbox, id_map, sort_wkspace);
  free(sort_wkspace);
  return 0;
}


uint32_t copy_and_dedup_sorted_strptrs_to_strbox(char** sorted_strptrs, uintptr_t str_ct, uintptr_t max_str_blen, char* strbox) {
  if (!str_ct) {
    return 0;
  }
  char** sorted_strptrs_iter = sorted_strptrs;
  char** sorted_strptrs_end = &(sorted_strptrs[str_ct]);
  uintptr_t write_idx = 0;
  uint32_t prev_slen = 0xffffffffU;
  char* prev_str = nullptr;
  do {
    char* cur_str = *sorted_strptrs_iter++;
    const uint32_t cur_slen = strlen(cur_str);
    if ((cur_slen != prev_slen) || memcmp(cur_str, prev_str, prev_slen)) {
      memcpy(&(strbox[write_idx * max_str_blen]), cur_str, cur_slen + 1);
      ++write_idx;
      prev_str = cur_str;
    }
  } while (sorted_strptrs_iter != sorted_strptrs_end);
  return write_idx;
}

uint32_t uint32arr_greater_than(const uint32_t* sorted_uint32_arr, uint32_t arr_length, uint32_t uii) {
  // (strangely, this seems to be equal to or better than std::lower_bound with
  // -O2 optimization, but can become much slower with -O3?)
  
  // assumes arr_length is nonzero, and sorted_uint32_arr is in nondecreasing
  // order.  (useful for searching marker_pos.)
  // also assumes arr_length < 2^31.
  // uii guaranteed to be larger than sorted_uint32_arr[min_idx - 1] if it
  // exists, but NOT necessarily sorted_uint32_arr[min_idx].
  int32_t min_idx = 0;
  // similarly, uii guaranteed to be no greater than
  // sorted_uint32_arr[max_idx + 1] if it exists, but not necessarily
  // sorted_uint32_arr[max_idx].  Signed integer since it could become -1, and
  // min_idx in turn is signed so comparisons are safe.
  int32_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uint32_t mid_idx = (((uint32_t)min_idx) + ((uint32_t)max_idx)) / 2;
    if (uii > sorted_uint32_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (uii > sorted_uint32_arr[((uint32_t)min_idx)]);
}

uintptr_t uint64arr_greater_than(const uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (ullii > sorted_uint64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (ullii > sorted_uint64_arr[((uintptr_t)min_idx)]);
}

uintptr_t doublearr_greater_than(const double* sorted_dbl_arr, uintptr_t arr_length, double dxx) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (dxx > sorted_dbl_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (dxx > sorted_dbl_arr[((uintptr_t)min_idx)]);
}

uintptr_t uint64arr_geq(const uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (ullii >= sorted_uint64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (ullii >= sorted_uint64_arr[((uintptr_t)min_idx)]);
}

uint32_t param_count(char** argv, uint32_t argc, uint32_t flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any nonnumeric parameter beginning with "-" as
  // optional.
  ++flag_idx;
  uint32_t cur_idx = flag_idx;
  while ((cur_idx < argc) && (!is_flag(argv[cur_idx]))) {
    ++cur_idx;
  }
  return cur_idx - flag_idx;
}

boolerr_t enforce_param_ct_range(const char* flag_name, uint32_t param_ct, uint32_t min_ct, uint32_t max_ct) {
  if (param_ct > max_ct) {
    if (max_ct > min_ct) {
      sprintf(g_logbuf, "Error: %s accepts at most %u parameter%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    } else {
      sprintf(g_logbuf, "Error: %s only accepts %u parameter%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    }
    return 1;
  }
  if (param_ct >= min_ct) {
    return 0;
  }
  if (min_ct == 1) {
    sprintf(g_logbuf, "Error: Missing %s parameter.\n", flag_name);
  } else {
    sprintf(g_logbuf, "Error: %s requires %s%u parameters.\n", flag_name, (min_ct < max_ct)? "at least " : "", min_ct);
  }
  return 1;
}

pglerr_t sort_cmdline_flags(uint32_t max_flag_blen, uint32_t flag_ct, char* flag_buf, uint32_t* flag_map) {
  // Assumes flag_ct is the number of flag (as opposed to value) parameters,
  // flag_buf[] points to a rectangular char* array (width max_flag_blen) of
  // flag names with leading dash(es) stripped, and flag_map[] maps flag_buf[]
  // entries to argv[] entries.
  // Lexicographically sorts flag_buf (updating flag_map in the process), and
  // then checks for duplicates.
  // Okay for flag_buf to contain entries with spaces (plink 1.9's alias
  // resolution takes advantage of this).
  assert(flag_ct); // this must be skipped if there are no flags at all
  if (sort_strbox_indexed_malloc(flag_ct, max_flag_blen, flag_buf, flag_map)) {
    return kPglRetNomem;
  }
  uint32_t prev_flag_len = strlen_se(flag_buf);
  char* prev_flag_ptr = flag_buf;
  for (uint32_t cur_flag_idx = 1; cur_flag_idx < flag_ct; ++cur_flag_idx) {
    char* cur_flag_ptr = &(prev_flag_ptr[max_flag_blen]);
    const uint32_t cur_flag_len = strlen_se(cur_flag_ptr);
    if ((prev_flag_len == cur_flag_len) && (!memcmp(prev_flag_ptr, cur_flag_ptr, cur_flag_len))) {
      cur_flag_ptr[cur_flag_len] = '\0'; // just in case of aliases
      fflush(stdout);
      fprintf(stderr, "Error: Duplicate --%s flag.\n", cur_flag_ptr);
      // g_stderr_written_to = 1;
      return kPglRetInvalidCmdline;
    }
    prev_flag_ptr = cur_flag_ptr;
    prev_flag_len = cur_flag_len;
  }
  return kPglRetSuccess;
}

pglerr_t init_logfile(uint32_t always_stderr, char* outname, char* outname_end) {
  strcpy(outname_end, ".log");
  g_logfile = fopen(outname, "w");
  if (!g_logfile) {
    fflush(stdout);
    fprintf(stderr, "Error: Failed to open %s for logging.\n", outname);
    // g_stderr_written_to = 1;
    return kPglRetOpenFail;
  }
  fprintf(always_stderr? stderr : stdout, "Logging to %s.\n", outname);
  return kPglRetSuccess;
}

boolerr_t cleanup_logfile(uint32_t print_end_time) {
  char* write_iter = strcpya(g_logbuf, "End time: ");
  time_t rawtime;
  time(&rawtime);
  write_iter = strcpya0(write_iter, ctime(&rawtime)); // has trailing \n
  if (print_end_time) {
    fputs(g_logbuf, stdout);
  }
  boolerr_t ret_boolerr = 0;
  if (g_logfile) {
    if (!g_log_failed) {
      logstr("\n");
      logstr(g_logbuf);
      if (fclose(g_logfile)) {
	fflush(stdout);
	fputs("Error: Failed to finish writing to log.\n", stderr);
	ret_boolerr = 1;
      }
    } else {
      fclose(g_logfile);
    }
    g_logfile = nullptr;
  }
  return ret_boolerr;
}

// manually managed, very large stack
unsigned char* g_bigstack_base = nullptr;
unsigned char* g_bigstack_end = nullptr;

uintptr_t detect_mb() {
  int64_t llxx;
  // return zero if detection failed
  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  int32_t mib[2];
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  size_t sztmp = sizeof(int64_t);
  sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
  llxx /= 1048576;
#else
#ifdef _WIN32
  MEMORYSTATUSEX memstatus;
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#else
  llxx = ((uint64_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
  return llxx;
}

uintptr_t get_default_alloc_mb() {
  const uintptr_t total_mb = detect_mb();
  if (!total_mb) {
    return kBigstackDefaultMb;
  }
  if (total_mb < (kBigstackMinMb * 2)) {
    return kBigstackMinMb;
  }
  return (total_mb / 2);
}

pglerr_t init_bigstack(uintptr_t malloc_size_mb, uintptr_t* malloc_mb_final_ptr, unsigned char** bigstack_ua_ptr) {
  // guarantee contiguous malloc space outside of main workspace
  unsigned char* bubble;
  if (pgl_malloc(kNonBigstackMin, &bubble)) {
    return kPglRetNomem;
  }
  assert(malloc_size_mb >= kBigstackMinMb);
#ifndef __LP64__
  assert(malloc_size_mb <= 2047);
#endif
  // don't use pgl_malloc here since we don't automatically want to set
  // g_failed_alloc_attempt_size on failure
  unsigned char* bigstack_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  // this is thwarted by overcommit, but still better than nothing...
  while (!bigstack_ua) {
    malloc_size_mb = (malloc_size_mb * 3) / 4;
    if (malloc_size_mb < kBigstackMinMb) {
      malloc_size_mb = kBigstackMinMb;
    }
    bigstack_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if ((!bigstack_ua) && (malloc_size_mb == kBigstackMinMb)) {
      // switch to "goto cleanup" pattern if any more exit points are needed
      g_failed_alloc_attempt_size = kBigstackMinMb * 1048576;
      free(bubble);
      return kPglRetNomem;
    }
  }
  // force 64-byte align to make cache line sensitivity work
  unsigned char* bigstack_initial_base = (unsigned char*)round_up_pow2((uintptr_t)bigstack_ua, kCacheline);
  g_bigstack_base = bigstack_initial_base;
  // last 512 bytes now reserved for g_one_char_strs
  g_bigstack_end = &(bigstack_initial_base[round_down_pow2(malloc_size_mb * 1048576 - 512 - (uintptr_t)(bigstack_initial_base - bigstack_ua), kCacheline)]);
  free(bubble);
  uintptr_t* one_char_iter = (uintptr_t*)g_bigstack_end;
#ifdef __LP64__
  // assumes little-endian
  uintptr_t cur_word = 0x3000200010000LLU;
  for (uint32_t uii = 0; uii < 64; ++uii) {
    *one_char_iter++ = cur_word;
    cur_word += 0x4000400040004LLU;
  }
#else
  uintptr_t cur_word = 0x10000;
  for (uint32_t uii = 0; uii < 128; ++uii) {
    *one_char_iter++ = cur_word;
    cur_word += 0x20002;
  }
#endif
  g_one_char_strs = (const char*)g_bigstack_end;

  // plink2 doesn't actually need these here, but short programs using
  // plink2_common benefit from this
  g_input_missing_geno_ptr = (const char*)(&(g_one_char_strs[96]));
  g_output_missing_geno_ptr = (const char*)(&(g_one_char_strs[92]));
  
  *malloc_mb_final_ptr = malloc_size_mb;
  *bigstack_ua_ptr = bigstack_ua;
  return kPglRetSuccess;
}


/*
uint32_t match_upper(const char* ss, const char* fixed_str) {
  char cc = *fixed_str++;
  do {
    if ((((unsigned char)(*ss++)) & 0xdf) != ((unsigned char)cc)) {
      return 0;
    }
    cc = *fixed_str++;
  } while (cc);
  return !(*ss);
}
*/

uint32_t match_upper_counted(const char* ss, const char* fixed_str, uint32_t ct) {
  for (uint32_t uii = 0; uii < ct; ++uii) {
    if ((((unsigned char)ss[uii]) & 0xdf) != ((unsigned char)fixed_str[uii])) {
      return 0;
    }
  }
  return 1;
}

/*
void str_toupper(char* ss) {
  while (1) {
    const uint32_t uii = (unsigned char)(*ss);
    if (!uii) {
      return;
    }
    if (((uint32_t)(uii - 97)) < 26) {
      // 'a' has ASCII code 97
      *ss = uii - 32;
    }
    ++ss;
  }
}

void buf_toupper(uint32_t slen, char* ss) {
  for (uint32_t pos = 0; pos < slen; ++pos) {
    const uint32_t uii = (unsigned char)(ss[pos]);
    if (((uint32_t)(uii - 97)) < 26) {
      ss[pos] = uii - 32;
    }
  }
}

void strcpy_toupper(char* target, const char* source) {
  while (1) {
    uint32_t uii = (unsigned char)(*source++);
    if (!uii) {
      return;
    }
    if (((uint32_t)(uii - 97)) < 26) {
      uii -= 32;
    }
    *target++ = uii;
  }
}
*/

uint32_t is_alphanumeric(const char* ss) {
  while (1) {
    uint32_t uii = (unsigned char)(*ss++);
    if (!uii) {
      return 1;
    }
    if (((uii - 48) > 9) && (((uii & 0xffffffdfU) - 65) > 25)) {
      return 0;
    }
  }
}

boolerr_t scan_posintptr(const char* ss, uintptr_t* valp) {
  // Reads an integer in [1, 2^kBitsPerWord - 1].  Assumes first character is
  // nonspace.
  assert(((unsigned char)ss[0]) > 32);
  uintptr_t val = (uintptr_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
#ifdef __LP64__
    if (val != 0xfffffffffffffffbLLU) {
      return 1;
    }
#else
    if (val != 0xfffffffbU) {
      return 1;
    }
#endif
    val = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (!val) {
    val = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
#ifdef __LP64__
  // limit is 20 digits, we've already read one
  const char* ss_limit = &(ss[20]);
#else
  const char* ss_limit = &(ss[10]);
#endif
  while (1) {
    const uintptr_t cur_digit = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    const uintptr_t cur_digit2 = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (ss == ss_limit) {
      if ((cur_digit2 < 10) || ((val >= (~k0LU) / 10) && ((val > (~k0LU) / 10) || (cur_digit > (~k0LU) % 10)))) {
	return 1;
      }
      *valp = val * 10 + cur_digit;
      return 0;
    }
    if (cur_digit2 >= 10) {
      *valp = val * 10 + cur_digit;
      return 0;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
  }
}

#ifdef __LP64__
static inline boolerr_t scanadv_uint_capped_finish(uint64_t cap, char** ss_ptr, uint32_t* valp) {
  unsigned char* ss = (unsigned char*)(*ss_ptr);
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = (uint64_t)(*ss++) - 48;
    if (cur_digit >= 10) {
      break;
    }
    // val = val * 10 + cur_digit;
    const uint64_t cur_digit2 = (uint64_t)(*ss++) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (val > cap) {
	return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (val > cap) {
      return 1;
    }
  }
  *valp = val;
  *ss_ptr = (char*)(&(ss[-1]));
  return 0;
}

boolerr_t scanadv_posint_capped(uint64_t cap, char** ss_ptr, uint32_t* valp) {
  unsigned char* ss = (unsigned char*)(*ss_ptr);
  *valp = (uint32_t)(*ss++) - 48;
  if (*valp >= 10) {
    if (*valp != 0xfffffffbU) {
      return 1;
    }
    *valp = (uint32_t)(*ss++) - 48;
    if (*valp >= 10) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = (uint32_t)(*ss++) - 48;
    if ((*valp) >= 10) {
      return 1;
    }
  }
  *ss_ptr = (char*)ss;
  return scanadv_uint_capped_finish(cap, ss_ptr, valp);
}

boolerr_t scanadv_uint_capped(uint64_t cap, char** ss_ptr, uint32_t* valp) {
  unsigned char* ss = (unsigned char*)(*ss_ptr);
  *valp = (uint32_t)(*ss++) - 48;
  if (*valp >= 10) {
    if (*valp != 0xfffffffbU) {
      // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
      if ((*valp != 0xfffffffdU) || (*ss != '0')) {
	return 1;
      }
      // accept "-0", "-00", etc.
      while (*(++ss) == '0');
      *valp = 0;
      *ss_ptr = (char*)ss;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
    }
    // accept leading '+'
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (*valp >= 10) {
      return 1;
    }
  }
  *ss_ptr = (char*)ss;
  return scanadv_uint_capped_finish(cap, ss_ptr, valp);
}
#else
boolerr_t scanadv_posint_capped32(uint32_t cap_div_10, uint32_t cap_mod_10, char** ss_ptr, uint32_t* valp) {
  unsigned char* ss = (unsigned char*)ss_ptr;
  uint32_t val = (uint32_t)(*ss++) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      return 1;
    }
    val = (uint32_t)(*ss++) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (!val) {
    val = (uint32_t)(*ss++);
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)(*ss++) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      *ss_ptr = (char*)(&(ss[-1]));
      return 0;
    }
    if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

boolerr_t scanadv_uint_capped32(uint32_t cap_div_10, uint32_t cap_mod_10, char** ss_ptr, uint32_t* valp) {
  unsigned char* ss = (unsigned char*)ss_ptr;
  uint32_t val = (uint32_t)(*ss++) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      if ((val != 0xfffffffd) || (*ss != '0')) {
	return 1;
      }
      while (*(++ss) == '0');
      *valp = 0;
      *ss_ptr = (char*)ss;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)(*ss++) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      *ss_ptr = (char*)(&(ss[-1]));
      return 0;
    }
    if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}
#endif

static const double kPositivePow10[16] = {1, 1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10, 1.0e11, 1.0e12, 1.0e13, 1.0e14, 1.0e15};
static const double kPositivePowTen16[16] = {1, 1.0e16, 1.0e32, 1.0e48, 1.0e64, 1.0e80, 1.0e96, 1.0e112, 1.0e128, 1.0e144, 1.0e160, 1.0e176, 1.0e192, 1.0e208, 1.0e224, 1.0e240};
static const double kNegativePow10[16] = {1, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14, 1.0e-15};
static const double kNegativePowTen16[8] = {1, 1.0e-16, 1.0e-32, 1.0e-48, 1.0e-64, 1.0e-80, 1.0e-96, 1.0e-112};

char* scanadv_double(char* ss, double* valp) {
  // requires first character to be nonspace (to succeed; it fails without
  //   segfaulting on space/eoln/null)
  // don't care about hexadecimal
  // ok to lose last ~2 bits of precision
  // ok if this yields incorrect results on >1GB strings
  // fail on nan/infinity/overflow instead of usual strtod behavior
  uint32_t cur_char_code = (unsigned char)(*ss);
  const uint32_t is_negative = (cur_char_code == 45);
  if (is_negative || (cur_char_code == 43)) {
    cur_char_code = (unsigned char)(*(++ss));
  }
  uint32_t cur_digit = cur_char_code - 48;
  int32_t e10 = 0;
  char* dot_ptr;
  int64_t digits;
#ifdef __LP64__
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
      if (cur_digit >= 10) {
	if (cur_digit == 0xfffffffeU) {
	  dot_ptr = ss;
	  goto scanadv_double_parse_decimal;
	}
	goto scanadv_double_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
      // this check should work differently in 32-bit version
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    // (could keep ~19 instead, but if we're systematically losing the last two
    // bits of precision anyway...)
    char* last_sig_fig_ptr = ss;
    do {
      cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
    } while (cur_digit < 10);
    e10 = (int32_t)((uint32_t)((uintptr_t)(ss - last_sig_fig_ptr))) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
	cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
      } while (cur_digit < 10);
    }
    goto scanadv_double_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = ss;
  cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits = cur_digit;
 scanadv_double_parse_decimal:
  while (1) {
    cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - (int32_t)((uint32_t)((uintptr_t)(ss - dot_ptr)));
      break;
    }
    digits = digits * 10 + cur_digit;
    if (digits >= 10000000000000000LL) {
      e10 = -(int32_t)((uint32_t)((uintptr_t)(ss - dot_ptr)));
      do {
	cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
      } while (cur_digit < 10);
      break;
    }
  }
 scanadv_double_parse_exponent:
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = (unsigned char)(*(++ss));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = (unsigned char)(*(++ss));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 107374182) {
	// may as well guard against exponent overflow
	if (!exp_is_negative) {
	  return nullptr;
	}
	*valp = 0;
	do {
	  cur_digit = ((unsigned char)(*(++ss))) - 48;
	} while (cur_digit < 10);
	return ss;
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ((unsigned char)(*(++ss))) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#else // not __LP64__
  int32_t digits_short;
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits_short = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
      if (cur_digit >= 10) {
	if (cur_digit == 0xfffffffeU) {
	  dot_ptr = ss;
	  goto scanadv_double_parse_decimal;
	}
	digits = digits_short;
	goto scanadv_double_parse_exponent;
      }
      digits_short = digits_short * 10 + cur_digit;
    } while (digits_short < 100000000);
    digits = digits_short;
    do {
      cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
      if (cur_digit >= 10) {
	if (cur_digit == 0xfffffffeU) {
	  dot_ptr = ss;
	  goto scanadv_double_parse_decimal_long;
	}
	goto scanadv_double_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    char* last_sig_fig_ptr = ss;
    do {
      cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
    } while (cur_digit < 10);
    e10 = (int32_t)((uint32_t)((uintptr_t)(ss - last_sig_fig_ptr))) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
	cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
      } while (cur_digit < 10);
    }
    goto scanadv_double_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = ss;
  cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits_short = cur_digit;
 scanadv_double_parse_decimal:
  while (1) {
    cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - (int32_t)((uint32_t)((uintptr_t)(ss - dot_ptr)));
      digits = digits_short;
      break;
    }
    digits_short = digits_short * 10 + cur_digit;
    if (digits_short >= 100000000) {
      digits = digits_short;
    scanadv_double_parse_decimal_long:
      while (1) {
	cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
	if (cur_digit >= 10) {
	  e10 = 1 - (int32_t)((uint32_t)((uintptr_t)(ss - dot_ptr)));
	  goto scanadv_double_parse_exponent;
	}
	digits = digits * 10 + cur_digit;
	if (digits >= 10000000000000000LL) {
	  e10 = -(int32_t)((uint32_t)((uintptr_t)(ss - dot_ptr)));
	  do {
	    cur_digit = ((uint32_t)((unsigned char)(*(++ss)))) - 48;
	  } while (cur_digit < 10);
	  goto scanadv_double_parse_exponent;
	}
      }
    }
  }
 scanadv_double_parse_exponent:
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = (unsigned char)(*(++ss));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = (unsigned char)(*(++ss));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 107374182) {
	// may as well guard against exponent overflow
	if (!exp_is_negative) {
	  return nullptr;
	}
	*valp = 0;
	do {
	  cur_digit = ((unsigned char)(*(++ss))) - 48;
	} while (cur_digit < 10);
	return ss;
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ((unsigned char)(*(++ss))) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#endif
  if (digits == 0) {
    *valp = 0;
    return ss;
  }
  if (is_negative) {
    digits = -digits;
  }
  double dxx = (double)digits;
  if (e10) {
    if (e10 < 0) {
      uint32_t pos_exp = (uint32_t)(-e10);
      dxx *= kNegativePow10[pos_exp & 15];
      pos_exp /= 16;
      if (pos_exp) {
	dxx *= kNegativePowTen16[pos_exp & 7];
	if (pos_exp > 7) {
	  if (pos_exp > 23) {
	    dxx = 0;
	  } else if (pos_exp > 15) {
	    dxx *= 1.0e-256;
	  } else {
	    dxx *= 1.0e-128;
	  }
	}
      }
    } else {
      uint32_t pos_exp = (uint32_t)e10;
      dxx *= kPositivePow10[pos_exp & 15];
      pos_exp /= 16;
      if (pos_exp) {
	dxx *= kPositivePowTen16[pos_exp & 15];
	if (pos_exp > 15) {
	  // overflow check
	  // last digits are "54" instead of "57" since that's the threshold
	  // beyond which multiply-by-1e256 overflows
	  if ((pos_exp > 31) || (dxx > 1.7976931348623154e52)) {
	    return nullptr;
	  }
	  dxx *= 1.0e256;
	}
      }
    }
  }
  *valp = dxx;
  return ss;
}

void get_top_two_ui(const uint32_t* __restrict uint_arr, uintptr_t uia_size, uintptr_t* __restrict top_idx_ptr, uintptr_t* __restrict second_idx_ptr) {
  assert(uia_size > 1);
  uintptr_t top_idx = (uint_arr[1] > uint_arr[0])? 1 : 0;
  uintptr_t second_idx = 1 ^ top_idx;
  uint32_t top_val = uint_arr[top_idx];
  uint32_t second_val = uint_arr[second_idx];
  uintptr_t cur_idx;
  uintptr_t cur_val;
  for (cur_idx = 2; cur_idx < uia_size; ++cur_idx) {
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
  }
  *top_idx_ptr = top_idx;
  *second_idx_ptr = second_idx;
}

char* comma_or_space_next_token_mult(char* sptr, uint32_t ct, uint32_t comma_delim) {
  assert(ct);
  if (!comma_delim) {
    return next_token_mult(sptr, ct);
  }
  if (!sptr) {
    return nullptr;
  }
  // assumes initial spaces in current token have been skipped
  // ok if we're at the end of the token
  unsigned char ucc = *sptr;
  assert(ucc != ' ');
  while (1) {
    // avoid strchr to keep "ASCII code < 32 == newline" consistent
    // (tab handling is quirky right now--permitted at the beginning of a
    // token, but treated as newline later--but it should never appear so
    // no point in e.g. adding an extra parameter to skip_initial_spaces();
    // just need to make sure the quirky behavior is consistent.)
    if (ucc < 32) {
      return nullptr;
    }
    if (ucc == ',') {
      do {
	ucc = (unsigned char)(*(++sptr));
      } while ((ucc == ' ') || (ucc == '\t'));
      if (!(--ct)) {
	return sptr;
      }
      continue;
    }
    ucc = (unsigned char)(*(++sptr));
  }
}

uint32_t count_tokens(const char* bufptr) {
  uint32_t token_ct = 0;
  // skip_initial_spaces/token_endnn spelled out due to const qualifier
  while ((*bufptr == ' ') || (*bufptr == '\t')) {
    ++bufptr;
  }
  while (!is_eoln_kns(*bufptr)) {
    ++token_ct;
    while (!is_space_or_eoln(*(++bufptr)));
    while ((*bufptr == ' ') || (*bufptr == '\t')) {
      ++bufptr;
    }
  }
  return token_ct;
}

/*
uint32_t comma_or_space_count_tokens(const char* bufptr, uint32_t comma_delim) {
  if (comma_delim) {
    // assumes nonempty line (treats trailing empty string as a token).
    uint32_t token_ct = 1;
    unsigned char ucc = (unsigned char)(*bufptr++);
    while (1) {
      if (ucc < 32) {
	return token_ct;
      }
      if (ucc == ',') {
	// spelled out due to const qualifier
	do {
	  ucc = (unsigned char)(*bufptr++);
	} while ((ucc == ' ') || (ucc == '\t'));
	token_ct++;
	continue;
      }
      ucc = (unsigned char)(*bufptr++);
    }
  }
  return count_tokens(bufptr);
}

uint32_t count_and_measure_multistr(const char* multistr, uintptr_t* max_blen_ptr) {
  // assumes multistr is nonempty
  assert(multistr[0]);
  uint32_t ct = 0;
  uintptr_t max_blen = *max_blen_ptr;
  do {
    const uintptr_t blen = strlen(multistr) + 1;
    if (blen > max_blen) {
      max_blen = blen;
    }
    multistr = &(multistr[blen]);
    ++ct;
  } while (*multistr);
  *max_blen_ptr = max_blen;
  return ct;
}
*/

boolerr_t count_and_measure_multistr_reverse_alloc(char* multistr, uintptr_t max_str_ct, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr, char*** strptr_arrp) {
  // assumes multistr is nonempty
  assert(multistr[0]);
  uintptr_t ct = 0;
  uintptr_t max_blen = *max_blen_ptr;
  char** strptr_arr_iter = *strptr_arrp;
  do {
    if (++ct > max_str_ct) {
      return 1;
    }
    const uintptr_t blen = strlen(multistr) + 1;
    if (blen > max_blen) {
      max_blen = blen;
    }
    *(--strptr_arr_iter) = multistr;
    multistr = &(multistr[blen]);
  } while (*multistr);
  *str_ct_ptr = ct;
  *max_blen_ptr = max_blen;
  *strptr_arrp = strptr_arr_iter;
  return 0;
}

// number-to-string encoders

static const uint16_t kDigitPair[100] = {
  0x3030, 0x3130, 0x3230, 0x3330, 0x3430, 0x3530, 0x3630, 0x3730, 0x3830, 0x3930,
  0x3031, 0x3131, 0x3231, 0x3331, 0x3431, 0x3531, 0x3631, 0x3731, 0x3831, 0x3931,
  0x3032, 0x3132, 0x3232, 0x3332, 0x3432, 0x3532, 0x3632, 0x3732, 0x3832, 0x3932,
  0x3033, 0x3133, 0x3233, 0x3333, 0x3433, 0x3533, 0x3633, 0x3733, 0x3833, 0x3933,
  0x3034, 0x3134, 0x3234, 0x3334, 0x3434, 0x3534, 0x3634, 0x3734, 0x3834, 0x3934,
  0x3035, 0x3135, 0x3235, 0x3335, 0x3435, 0x3535, 0x3635, 0x3735, 0x3835, 0x3935,
  0x3036, 0x3136, 0x3236, 0x3336, 0x3436, 0x3536, 0x3636, 0x3736, 0x3836, 0x3936,
  0x3037, 0x3137, 0x3237, 0x3337, 0x3437, 0x3537, 0x3637, 0x3737, 0x3837, 0x3937,
  0x3038, 0x3138, 0x3238, 0x3338, 0x3438, 0x3538, 0x3638, 0x3738, 0x3838, 0x3938,
  0x3039, 0x3139, 0x3239, 0x3339, 0x3439, 0x3539, 0x3639, 0x3739, 0x3839, 0x3939};

char* uint32toa(uint32_t uii, char* start) {
  // Memory-efficient fast integer writer.  (You can do a bit better sometimes
  // by using a larger lookup table, but on average I doubt that pays off.)
  // Returns a pointer to the end of the integer (not null-terminated).
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      *start++ = '0' + uii;
      return start;
    }
    if (uii < 100) {
      goto uint32toa_2;
    }
    quotient = uii / 100;
    *start++ = '0' + quotient;
  } else {
    if (uii < 10000000) {
      if (uii >= 100000) {
	if (uii < 1000000) {
	  goto uint32toa_6;
	}
	quotient = uii / 1000000;
	*start++ = '0' + quotient;
	goto uint32toa_6b;
      }
      if (uii < 10000) {
	goto uint32toa_4;
      }
      quotient = uii / 10000;
      *start++ = '0' + quotient;
    } else {
      if (uii >= 100000000) {
	quotient = uii / 100000000;
	if (uii >= 1000000000) {
	  start = memcpya(start, &(kDigitPair[quotient]), 2);
	} else {
	  *start++ = '0' + quotient;
	}
	uii -= 100000000 * quotient;
      }
      quotient = uii / 1000000;
      start = memcpya(start, &(kDigitPair[quotient]), 2);
    uint32toa_6b:
      uii -= 1000000 * quotient;
    uint32toa_6:
      quotient = uii / 10000;
      start = memcpya(start, &(kDigitPair[quotient]), 2);
    }
    uii -= 10000 * quotient;
  uint32toa_4:
    // could make a uitoa_z4() call here, but that's slightly slower
    quotient = uii / 100;
    start = memcpya(start, &(kDigitPair[quotient]), 2);
  }
  uii -= 100 * quotient;
 uint32toa_2:
  return memcpya(start, &(kDigitPair[uii]), 2);
}

char* int32toa(int32_t ii, char* start) {
  uint32_t uii = ii;
  if (ii < 0) {
    // -INT_MIN is undefined, but negating the unsigned int equivalent works
    *start++ = '-';
    uii = -uii;
  }
  return uint32toa(uii, start);
}

char* uitoa_z4(uint32_t uii, char* start) {
  uint32_t quotient = uii / 100;
  assert(quotient < 100);
  uii -= 100 * quotient;
  start = memcpya(start, &(kDigitPair[quotient]), 2);
  return memcpya(start, &(kDigitPair[uii]), 2);
}

char* uitoa_z6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  start = memcpya(start, &(kDigitPair[quotient]), 2);
  return uitoa_z4(uii - 10000 * quotient, start);
}

char* uitoa_z8(uint32_t uii, char* start) {
  uint32_t quotient = uii / 1000000;
  start = memcpya(start, &(kDigitPair[quotient]), 2);
  return uitoa_z6(uii - 1000000 * quotient, start);
}

char* int64toa(int64_t llii, char* start) {
  uint64_t ullii = llii;
  uint64_t top_digits;
  uint32_t bottom_eight;
  uint32_t middle_eight;
  if (llii < 0) {
    *start++ = '-';
    ullii = -ullii;
  }
  if (ullii <= 0xffffffffLLU) {
    return uint32toa((uint32_t)ullii, start);
  }
  top_digits = ullii / 100000000;
  bottom_eight = (uint32_t)(ullii - (top_digits * 100000000));
  if (top_digits <= 0xffffffffLLU) {
    start = uint32toa((uint32_t)top_digits, start);
    return uitoa_z8(bottom_eight, start);
  }
  ullii = top_digits / 100000000;
  middle_eight = (uint32_t)(top_digits - (ullii * 100000000));
  start = uint32toa((uint32_t)ullii, start);
  start = uitoa_z8(middle_eight, start);
  return uitoa_z8(bottom_eight, start);
}


char* uitoa_trunc4(uint32_t uii, char* start) {
  uint32_t quotient = uii / 100;
  memcpy(start, &(kDigitPair[quotient]), 2);
  uii -= 100 * quotient;
  if (uii) {
    start += 2;
    memcpy(start, &(kDigitPair[uii]), 2);
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  memcpy(start, &(kDigitPair[quotient]), 2);
  uii -= 10000 * quotient;
  if (uii) {
    quotient = uii / 100;
    start += 2;
    memcpy(start, &(kDigitPair[quotient]), 2);
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
      memcpy(start, &(kDigitPair[uii]), 2);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc8(uint32_t uii, char* start) {
  uint32_t quotient = uii / 1000000;
  memcpy(start, &(kDigitPair[quotient]), 2);
  uii -= 1000000 * quotient;
  if (uii) {
    quotient = uii / 10000;
    start += 2;
    memcpy(start, &(kDigitPair[quotient]), 2);
    uii -= 10000 * quotient;
    if (uii) {
      quotient = uii / 100;
      start += 2;
      memcpy(start, &(kDigitPair[quotient]), 2);
      uii -= 100 * quotient;
      if (uii) {
	start += 2;
	memcpy(start, &(kDigitPair[uii]), 2);
      }
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* rtoa_p5(uint32_t remainder, char* start) {
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  uint32_t quotient = remainder / 1000;
  memcpy(start, &(kDigitPair[quotient]), 2);
  remainder -= 1000 * quotient;
  if (remainder) {
    quotient = remainder / 10;
    start += 2;
    memcpy(start, &(kDigitPair[quotient]), 2);
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

static inline char* qrtoa_1p5(uint32_t quotient, uint32_t remainder, char* start) {
  *start++ = '0' + quotient;
  return rtoa_p5(remainder, start);
}

static inline char* qrtoa_1p7(uint32_t quotient, uint32_t remainder, char* start) {
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 100000;
  memcpy(start, &(kDigitPair[quotient]), 2);
  remainder -= 100000 * quotient;
  if (remainder) {
    quotient = remainder / 1000;
    start += 2;
    memcpy(start, &(kDigitPair[quotient]), 2);
    remainder -= 1000 * quotient;
    if (remainder) {
      quotient = remainder / 10;
      start += 2;
      memcpy(start, &(kDigitPair[quotient]), 2);
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

// Okay, time to do banker's rounding when printing doubles.  14 digits of
// precision are used in judging equality to 0.5 (actual precision of doubles
// is 15-17 digits); the intention is to capture all directly loaded or exactly
// computed edge cases (so enough tolerance is needed to survive the internal
// multiplications by powers of 10, etc.), while rounding a negligible number
// of honest-to-god 0.4999999s up and 0.5000001s down.
// To avoid inadvertent printing of an extra digit, there's a deliberate gap
// between the 99.9994999...-type bounds and the largest numbers that would
// actually round down.
static const double kBankerRound6[] = {0.4999995, 0.5000005};
static const double kBankerRound8[] = {0.499999995, 0.500000005};

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

static inline void double_bround6(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000000;
  *remainderp = remainder - (*quotientp) * 1000000; 
}

static inline void double_bround7(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10000000;
  *remainderp = remainder - (*quotientp) * 10000000; 
}

char* dtoa_so6(double dxx, char* start) {
  // 6 sig fig number, 0.999995 <= dxx < 999999.5
  // 'so' = "significand only"
  // Just hardcoding all six cases, in the absence of a better approach...
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.999949999999) {
    if (dxx < 9.9999949999999) {
      double_bround5(dxx, kBankerRound8, &quotient, &remainder);
      return qrtoa_1p5(quotient, remainder, start);
    }
    double_bround4(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy(start, &(kDigitPair[quotient]), 2);
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so6_pretail:
      memcpy(start, &(kDigitPair[remainder]), 2);
    }
  dtoa_so6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  }
  if (dxx < 9999.9949999999) {
    if (dxx < 999.99949999999) {
      double_bround3(dxx, kBankerRound8, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(kDigitPair[quotient]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto dtoa_so6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    double_bround2(dxx, kBankerRound8, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so6_pretail;
  }
  if (dxx >= 99999.949999999) {
    return uitoa_z6(double_bround(dxx, kBankerRound8), start);
  }
  double_bround1(dxx, kBankerRound8, &uii, &remainder);
  quotient = uii / 10000;
  *start = '0' + quotient;
  uii -= 10000 * quotient;
  quotient = uii / 100;
  start = memcpya(&(start[1]), &(kDigitPair[quotient]), 2);
  uii = uii - 100 * quotient;
  start = memcpya(start, &(kDigitPair[uii]), 2);
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  *start = '0' + remainder;
  return &(start[1]);
}

char* dtoa_so8(double dxx, char* start) {
  // 8 sig fig number, 0.99999995 <= dxx < 99999999.5
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.999999499999) {
    if (dxx < 9.9999999499999) {
      double_bround7(dxx, kBankerRound6, &quotient, &remainder);
      return qrtoa_1p7(quotient, remainder, start);
    }
    double_bround6(dxx, kBankerRound6, &quotient, &remainder);
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 10000;
    memcpy(start, &(kDigitPair[quotient]), 2);
    remainder -= 10000 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so8_pretail4:
      quotient = remainder / 100;
      memcpy(start, &(kDigitPair[quotient]), 2);
      remainder -= 100 * quotient;
      if (remainder) {
	start += 2;
      dtoa_so8_pretail2:
        memcpy(start, &(kDigitPair[remainder]), 2);
      }
    }
  dtoa_so8_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  }
  if (dxx < 9999.9999499999) {
    if (dxx < 999.99999499999) {
      double_bround5(dxx, kBankerRound6, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(kDigitPair[quotient]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 1000;
      memcpy(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 1000;
      if (!remainder) {
        goto dtoa_so8_tail;
      }
      start += 2;
    dtoa_so8_pretail3:
      quotient = remainder / 10;
      memcpy(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
	goto dtoa_so8_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    double_bround4(dxx, kBankerRound6, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so8_pretail4;
  }
  if (dxx < 999999.99499999) {
    if (dxx < 99999.999499999) {
      double_bround3(dxx, kBankerRound6, &uii, &remainder);
      quotient = uii / 10000;
      *start = '0' + quotient;
      uii -= 10000 * quotient;
      quotient = uii / 100;
      start = memcpya(&(start[1]), &(kDigitPair[quotient]), 2);
      uii -= 100 * quotient;
      start = memcpya(start, &(kDigitPair[uii]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      goto dtoa_so8_pretail3;
    }
    double_bround2(dxx, kBankerRound6, &uii, &remainder);
    quotient = uii / 10000;
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    uii -= 100 * quotient;
    start = memcpya(start, &(kDigitPair[uii]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so8_pretail2;
  }
  if (dxx >= 9999999.9499999) {
    return uitoa_z8(double_bround(dxx, kBankerRound6), start);
  }
  double_bround1(dxx, kBankerRound6, &uii, &remainder);
  quotient = uii / 1000000;
  *start = '0' + quotient;
  uii -= 1000000 * quotient;
  quotient = uii / 10000;
  start = memcpya(&(start[1]), &(kDigitPair[quotient]), 2);
  uii -= 10000 * quotient;
  quotient = uii / 100;
  start = memcpya(start, &(kDigitPair[quotient]), 2);
  uii -= 100 * quotient;
  start = memcpya(start, &(kDigitPair[uii]), 2);
  if (!remainder) {
    return start;
  }
  *start = '.';
  start[1] = '0' + remainder;
  return &(start[2]);
}

char* dtoa_g(double dxx, char* start) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpyl3a(start, "nan");
  }
  if (dxx < 0) {
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
	}
	if (dxx < 9.9999949999999e-256) {
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
      ++xp10;
    }
    double_bround5(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya(qrtoa_1p5(quotient, remainder, start), "e-", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 999999.49999999) {
    // 6 sig fig exponential notation, large
    if (dxx >= 9.9999949999999e15) {
      if (dxx >= 9.9999949999999e127) {
	if (dxx > DBL_MAX) {
	  return memcpyl3a(start, "inf");
	}
	if (dxx >= 9.9999949999999e255) {
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
    double_bround5(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya(qrtoa_1p5(quotient, remainder, start), "e+", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 0.99999949999999) {
    return dtoa_so6(dxx, start);
  }
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
  return uitoa_trunc6(double_bround(dxx * 1000000, kBankerRound8), start);
}

char* dtoa_g_p8(double dxx, char* start) {
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpyl3a(start, "nan");
  }
  if (dxx < 0) {
    *wpos++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9999999499999e-5) {
    // 8 sig fig exponential notation, small
    if (dxx < 9.9999999499999e-16) {
      if (dxx < 9.9999999499999e-128) {
	if (dxx == 0.0) {
	  *start = '0';
	  return &(start[1]);
        }
	if (dxx < 9.9999999499999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9999999499999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9999999499999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9999999499999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9999999499999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9999999499999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999999499999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999999499999e-1) {
      dxx *= 10;
      ++xp10;
    }
    double_bround7(dxx, kBankerRound6, &quotient, &remainder);
    wpos = qrtoa_1p7(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      start = memcpya(start, wbuf, remainder);
      quotient = xp10 / 100;
      start = memcpyax(start, "e-", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      start = memcpya(start, wbuf, remainder);
      start = memcpya(start, "e-", 2);
    }
    return memcpya(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 99999999.499999) {
    // 8 sig fig exponential notation, large
    if (dxx >= 9.9999999499999e15) {
      if (dxx >= 9.9999999499999e127) {
	if (dxx > DBL_MAX) {
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  }
	  return memcpya(start, "-inf", 4);
	}
	if (dxx >= 9.9999999499999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9999999499999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9999999499999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9999999499999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9999999499999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9999999499999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999999499999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999999499999e0) {
      dxx *= 1.0e-1;
      ++xp10;
    }
    double_bround7(dxx, kBankerRound6, &quotient, &remainder);
    wpos = qrtoa_1p7(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      start = memcpya(start, wbuf, remainder);
      quotient = xp10 / 100;
      start = memcpyax(start, "e+", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      start = memcpya(start, wbuf, remainder);
      start = memcpya(start, "e+", 2);
    }
    return memcpya(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 0.99999999499999) {
    wpos = dtoa_so8(dxx, wpos);
  } else {
    // 8 sig fig decimal, no less than ~0.0001
    wpos = memcpya(wpos, "0.", 2);
    if (dxx < 9.9999999499999e-3) {
      dxx *= 100;
      wpos = memcpya(wpos, "00", 2);
    }
    if (dxx < 9.9999999499999e-2) {
      dxx *= 10;
      *wpos++ = '0';
    }
    wpos = uitoa_trunc8(double_bround(dxx * 100000000, kBankerRound6), wpos);
  }
  remainder = wpos - wbuf;
  return memcpya(start, wbuf, remainder);
}

static_assert(kDosageMax == 32768, "print_dosage() needs to be updated.");
char* print_dosage(uint64_t dosage, char* start) {
  // 3 digit precision seems like the best compromise between accuracy and
  // avoidance of rounding ugliness
  // (Rounding ugliness is not actually hidden for e.g. 1000 Genomes phase 1,
  // since there are lots of 0.05 and 0.1 dosages which all get rounded in the
  // same direction; oh well.)

  // +16 since we need to round .99951 up to 1
  const uint64_t dosage_p16 = dosage + 16;
  start = uint32toa(dosage_p16 / kDosageMax, start);
  const uint32_t remainder_p16 = ((uint32_t)dosage_p16) & (kDosageMax - 1);
  if (remainder_p16 < 33) {
    return start;
  }
  // (1000 * remainder + 16384) / 32768
  //   1/16 = .0625 -> print 0.062
  //   3/16 = .1875 -> print 0.188
  //   5/16 = .3125 -> print 0.312
  // const uint32_t three_decimal_places = ((125 * remainder + 2048) / 4096) - ((remainder % 8192) == 2048);
  const uint32_t three_decimal_places = ((125 * remainder_p16 + 48) / 4096) - ((remainder_p16 % 8192) == 4048);
  // three_decimal_places guaranteed to be nonzero here
  *start++ = '.';
  const uint32_t first_decimal_place = three_decimal_places / 100;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_two_digits = three_decimal_places - first_decimal_place * 100;
  if (last_two_digits) {
    memcpy(start, &(kDigitPair[last_two_digits]), 2);
    return &(start[1 + (start[1] != '0')]);
  }
  return start;
}
/*
char* dtoa_f_p5_clipped(double dxx, char* start) {
  if (dxx != dxx) {
    return memcpyl3a(start, "nan");
  }
  if (dxx < 0.0) {
    // note that "-0" will be printed for very small negative numbers; do we
    // want this?
    *start++ = '-';
    dxx = -dxx;
  }
#ifdef __LP64__
  if (dxx < 4294967295.999994) {
    // We could use different levels of banker's rounding for different-size
    // quotients, but that's overkill for now; revisit after basic dosage
    // support is working.
    dxx *= 100000;
    uint64_t remainder = (int64_t)dxx;
    remainder += (int64_t)((dxx - ((int64_t)remainder)) + kBankerRound6[remainder & 1]);
    uint64_t quotient = remainder / 100000;
    remainder = remainder - quotient * 100000;
    start = uint32toa(quotient, start);
    return rtoa_p5(remainder, start);
  }
#else
  if (dxx < 2147483647.999994) {
    // avoid 64-bit integer math in 32-bit build.
    // (todo: a bit of benchmarking)
    const uintptr_t quotient = (intptr_t)dxx;
    const double remainder_d = (dxx - ((intptr_t)quotient)) * 100000;
    const uint32_t remainder_d_trunc = (int32_t)remainder_d;
    const uint32_t remainder = (int32_t)(remainder_d + kBankerRound6[remainder_d_trunc & 1]);
    start = uint32toa(quotient, start);
    return rtoa_p5(remainder, start);
  }
#endif
  if (dxx == INFINITY) {
    return memcpyl3a(start, "inf");
  }
  // just punt larger numbers to glibc for now, this isn't a bottleneck
  start += sprintf(start, "%.5f", dxx);
  // .5f doesn't strip trailing zeroes, do that manually
  for (uint32_t uii = 0; uii < 5; ++uii) {
    if (start[-1] != '0') {
      return start;
    }
    --start;
  }
  return &(start[-1]); // strip the decimal point
}
*/


// Briefly had banker's rounding for floats, but then I realized that the only
// float-printing function calls are --make-grm related, they all request 6-7
// digits of precision, and at that point it's impossible to distinguish exact
// 0.5-matches in the remainder.  So we just have generic rounding functions
// here, with similar interfaces to the double-rounding functions to minimize
// the need for separate reasoning about this code.
static inline uint32_t float_round(float fxx) {
  return (uint32_t)((int32_t)(fxx + 0.5));
}

static inline void float_round1(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 10);
  *quotientp = remainder / 10;
  *remainderp = remainder - (*quotientp) * 10;
}

static inline void float_round2(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 100);
  *quotientp = remainder / 100;
  *remainderp = remainder - (*quotientp) * 100;
}

static inline void float_round3(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 1000);
  *quotientp = remainder / 1000;
  *remainderp = remainder - (*quotientp) * 1000;
}

static inline void float_round4(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 10000);
  *quotientp = remainder / 10000;
  *remainderp = remainder - (*quotientp) * 10000;
}

static inline void float_round5(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 100000);
  *quotientp = remainder / 100000;
  *remainderp = remainder - (*quotientp) * 100000;
}

char* ftoa_so6(float fxx, char* start) {
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  // difference between consecutive floats near 10 can be as large as
  // 10 * 2^{-23}, which is just under 1.2e-6.  So, to avoid printing an extra
  // digit, we have to set this bound to be robust to an addition error of size
  // 6e-7.
  // (possible todo: just brute-force test this on all <2^32 possible floats
  // and look for a better threshold)
  if (fxx < 99.999944) {
    if (fxx < 9.9999944) {
      float_round5(fxx, &quotient, &remainder);
      return qrtoa_1p5(quotient, remainder, start);
    }
    float_round4(fxx, &quotient, &remainder);
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy(start, &(kDigitPair[quotient]), 2);
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    ftoa_so6_pretail:
      memcpy(start, &(kDigitPair[remainder]), 2);
    }
  ftoa_so6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  }
  if (fxx < 9999.9944) {
    if (fxx < 999.99944) {
      float_round3(fxx, &uii, &remainder);
      quotient = uii / 100;
      *start = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(&(start[1]), &(kDigitPair[quotient]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto ftoa_so6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    float_round2(fxx, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto ftoa_so6_pretail;
  }
  if (fxx >= 99999.944) {
    return uitoa_z6(float_round(fxx), start);
  }
  float_round1(fxx, &uii, &remainder);
  quotient = uii / 10000;
  *start = '0' + quotient;
  uii -= 10000 * quotient;
  quotient = uii / 100;
  start = memcpya(&(start[1]), &(kDigitPair[quotient]), 2);
  uii = uii - 100 * quotient;
  start = memcpya(start, &(kDigitPair[uii]), 2);
  if (!remainder) {
    return start;
  }
  *start = '.';
  start[1] = '0' + remainder;
  return &(start[2]);
}

char* ftoa_g(float fxx, char* start) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  if (fxx != fxx) {
    return memcpyl3a(start, "nan");
  }
  if (fxx < 0) {
    *start++ = '-';
    fxx = -fxx;
  }
  if (fxx < 9.9999944e-5) {
    if (fxx < 9.9999944e-16) {
      if (fxx == 0.0) {
	*start = '0';
	return &(start[1]);
      }
      if (fxx < 9.9999944e-32) {
	fxx *= 1.0e32;
	xp10 |= 32;
      } else {
	fxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (fxx < 9.9999944e-8) {
      fxx *= 100000000;
      xp10 |= 8;
    }
    if (fxx < 9.9999944e-4) {
      fxx *= 10000;
      xp10 |= 4;
    }
    if (fxx < 9.9999944e-2) {
      fxx *= 100;
      xp10 |= 2;
    }
    if (fxx < 9.9999944e-1) {
      fxx *= 10;
      ++xp10;
    }
    float_round5(fxx, &quotient, &remainder);
    return memcpya(memcpya(qrtoa_1p5(quotient, remainder, start), "e-", 2), &(kDigitPair[xp10]), 2);
  }
  if (fxx >= 999999.44) {
    if (fxx >= 9.9999944e15) {
      if (fxx > FLT_MAX) {
	return memcpyl3a(start, "inf");
      }
      if (fxx >= 9.9999944e31) {
	fxx *= 1.0e-32;
	xp10 |= 32;
      } else {
	fxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (fxx >= 9.9999944e7) {
      fxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (fxx >= 9.9999944e3) {
      fxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (fxx >= 9.9999944e1) {
      fxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (fxx >= 9.9999944e0) {
      fxx *= 1.0e-1;
      ++xp10;
    }
    float_round5(fxx, &quotient, &remainder);
    return memcpya(memcpya(qrtoa_1p5(quotient, remainder, start), "e+", 2), &(kDigitPair[xp10]), 2);
  }
  if (fxx >= 0.99999944) {
    return ftoa_so6(fxx, start);
  }
  // 6 sig fig decimal, no less than ~0.0001
  start = memcpya(start, "0.", 2);
  if (fxx < 9.9999944e-3) {
    fxx *= 100;
    start = memcpya(start, "00", 2);
  }
  if (fxx < 9.9999944e-2) {
    fxx *= 10;
    *start++ = '0';
  }
  return uitoa_trunc6(float_round(fxx * 1000000), start);
}


/*
void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp) {
  // Enables fast integer division by a constant not known until runtime.  See
  // http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html .
  // Assumes divisor is not zero, of course.
  // May want to populate a struct instead.
  // (May not need this any more?)
  assert(divisor);
  if (!(divisor & (divisor - 1))) {
    // power of 2
    *multp = 1;
    *pre_shiftp = 0;
    *post_shiftp = __builtin_ctz(divisor);
    *incrp = 0;
    return;
  }
  uint32_t quotient = 0x80000000U / divisor;
  uint32_t remainder = 0x80000000U - (quotient * divisor);
  const uint32_t ceil_log_2_d = 32 - __builtin_clz(divisor);
  uint32_t down_multiplier = 0;
  uint32_t down_exponent = 0;
  uint32_t has_magic_down = 0;
  uint32_t exponent;
  for (exponent = 0; ; ++exponent) {
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
  }
  if (divisor & 1) {
    *multp = down_multiplier;
    *pre_shiftp = 0;
    *post_shiftp = 32 + down_exponent;
    *incrp = 1;
    return;
  }
  *pre_shiftp = __builtin_ctz(divisor);
  uint32_t dummy;
  magic_num(divisor >> (*pre_shiftp), multp, &dummy, post_shiftp, incrp);
}
*/


void fill_bits_nz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr) {
  assert(end_idx > start_idx);
  uintptr_t maj_start = start_idx / kBitsPerWord;
  uintptr_t maj_end = end_idx / kBitsPerWord;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] |= (k1LU << (end_idx % kBitsPerWord)) - (k1LU << (start_idx % kBitsPerWord));
  } else {
    bitarr[maj_start] |= ~((k1LU << (start_idx % kBitsPerWord)) - k1LU);
    fill_ulong_one(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = end_idx % kBitsPerWord;
    if (minor) {
      bitarr[maj_end] |= (k1LU << minor) - k1LU;
    }
  }
}

void clear_bits_nz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr) {
  assert(end_idx > start_idx);
  uintptr_t maj_start = start_idx / kBitsPerWord;
  uintptr_t maj_end = end_idx / kBitsPerWord;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] &= ~((k1LU << (end_idx % kBitsPerWord)) - (k1LU << (start_idx % kBitsPerWord)));
  } else {
    bitarr[maj_start] &= ((k1LU << (start_idx % kBitsPerWord)) - k1LU);
    fill_ulong_zero(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = end_idx % kBitsPerWord;
    if (minor) {
      bitarr[maj_end] &= ~((k1LU << minor) - k1LU);
    }
  }
}

#ifdef __LP64__
uintptr_t next_set_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % kBitsPerWord);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  return ((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + CTZLU(ulii);
}
#endif

uint32_t next_unset(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil) {
  assert(ceil >= 1);
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % kBitsPerWord);
  if (ulii) {
    loc += CTZLU(ulii);
    return MINV(loc, ceil);
  }
  const uintptr_t* bitarr_last = &(bitarr[(ceil - 1) / kBitsPerWord]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (ulii == ~k0LU);
  loc = ((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + CTZLU(~ulii);
  return MINV(loc, ceil);
}


boolerr_t bigstack_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = (unsigned char*)bigstack_alloc(ct);
  if (!(*uc_arr_ptr)) {
    return 1;
  }
  memset(*uc_arr_ptr, 0, ct);
  return 0;
}

boolerr_t bigstack_calloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = (double*)bigstack_alloc(ct * sizeof(double));
  if (!(*d_arr_ptr)) {
    return 1;
  }
  fill_double_zero(ct, *d_arr_ptr);
  return 0;
}

boolerr_t bigstack_calloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = (float*)bigstack_alloc(ct * sizeof(float));
  if (!(*f_arr_ptr)) {
    return 1;
  }
  fill_float_zero(ct, *f_arr_ptr);
  return 0;
}

boolerr_t bigstack_calloc_usi(uintptr_t ct, uint16_t** usi_arr_ptr) {
  *usi_arr_ptr = (uint16_t*)bigstack_alloc(ct * sizeof(int16_t));
  if (!(*usi_arr_ptr)) {
    return 1;
  }
  memset(*usi_arr_ptr, 0, ct * sizeof(int16_t));
  return 0;
}

boolerr_t bigstack_calloc_ui(uintptr_t ct, uint32_t** ui_arr_ptr) {
  *ui_arr_ptr = (uint32_t*)bigstack_alloc(ct * sizeof(int32_t));
  if (!(*ui_arr_ptr)) {
    return 1;
  }
  fill_uint_zero(ct, *ui_arr_ptr);
  return 0;
}

boolerr_t bigstack_calloc_ul(uintptr_t ct, uintptr_t** ul_arr_ptr) {
  *ul_arr_ptr = (uintptr_t*)bigstack_alloc(ct * sizeof(intptr_t));
  if (!(*ul_arr_ptr)) {
    return 1;
  }
  fill_ulong_zero(ct, *ul_arr_ptr);
  return 0;
}

boolerr_t bigstack_calloc_ull(uintptr_t ct, uint64_t** ull_arr_ptr) {
  *ull_arr_ptr = (uint64_t*)bigstack_alloc(ct * sizeof(int64_t));
  if (!(*ull_arr_ptr)) {
    return 1;
  }
  fill_ull_zero(ct, *ull_arr_ptr);
  return 0;
}

boolerr_t bigstack_end_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = (unsigned char*)bigstack_end_alloc(ct);
  if (!(*uc_arr_ptr)) {
    return 1;
  }
  memset(*uc_arr_ptr, 0, ct);
  return 0;
}

boolerr_t bigstack_end_calloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = (double*)bigstack_end_alloc(ct * sizeof(double));
  if (!(*d_arr_ptr)) {
    return 1;
  }
  fill_double_zero(ct, *d_arr_ptr);
  return 0;
}

boolerr_t bigstack_end_calloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = (float*)bigstack_end_alloc(ct * sizeof(float));
  if (!(*f_arr_ptr)) {
    return 1;
  }
  fill_float_zero(ct, *f_arr_ptr);
  return 0;
}

boolerr_t bigstack_end_calloc_ui(uintptr_t ct, uint32_t** ui_arr_ptr) {
  *ui_arr_ptr = (uint32_t*)bigstack_end_alloc(ct * sizeof(int32_t));
  if (!(*ui_arr_ptr)) {
    return 1;
  }
  fill_uint_zero(ct, *ui_arr_ptr);
  return 0;
}

boolerr_t bigstack_end_calloc_ul(uintptr_t ct, uintptr_t** ul_arr_ptr) {
  *ul_arr_ptr = (uintptr_t*)bigstack_end_alloc(ct * sizeof(intptr_t));
  if (!(*ul_arr_ptr)) {
    return 1;
  }
  fill_ulong_zero(ct, *ul_arr_ptr);
  return 0;
}

boolerr_t bigstack_end_calloc_ull(uintptr_t ct, uint64_t** ull_arr_ptr) {
  *ull_arr_ptr = (uint64_t*)bigstack_end_alloc(ct * sizeof(int64_t));
  if (!(*ull_arr_ptr)) {
    return 1;
  }
  fill_ull_zero(ct, *ull_arr_ptr);
  return 0;
}


void bitarr_invert(uintptr_t bit_ct, uintptr_t* bitarr) {
  uintptr_t* bitarr_stop = &(bitarr[bit_ct / kBitsPerWord]);
  while (bitarr < bitarr_stop) {
    *bitarr = ~(*bitarr);
    ++bitarr;
  }
  const uint32_t trailing_bit_ct = bit_ct % kBitsPerWord;
  if (trailing_bit_ct) {
    *bitarr = (~(*bitarr)) & ((k1LU << trailing_bit_ct) - k1LU);
  }
}

void bitarr_invert_copy(const uintptr_t* __restrict source_bitarr, uintptr_t bit_ct, uintptr_t* __restrict target_bitarr) {
  const uintptr_t* source_bitarr_stop = &(source_bitarr[bit_ct / kBitsPerWord]);
  while (source_bitarr < source_bitarr_stop) {
    *target_bitarr++ = ~(*source_bitarr++);
  }
  const uint32_t trailing_bit_ct = bit_ct % kBitsPerWord;
  if (trailing_bit_ct) {
    *target_bitarr = (~(*source_bitarr)) & ((k1LU << trailing_bit_ct) - k1LU);
  }
}

void bitvec_and_copy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec) {
#ifdef __LP64__
  vul_t* target_bitvvec = (vul_t*)target_bitvec;
  const vul_t* source1_bitvvec = (const vul_t*)source1_bitvec;
  const vul_t* source2_bitvvec = (const vul_t*)source2_bitvec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  for (uintptr_t ulii = 0; ulii < full_vec_ct; ++ulii) {
    target_bitvvec[ulii] = source1_bitvvec[ulii] & source2_bitvvec[ulii];
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = source1_bitvec[base_idx] & source2_bitvec[base_idx];
    target_bitvec[base_idx + 1] = source1_bitvec[base_idx + 1] & source2_bitvec[base_idx + 1];
  }
  #endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = source1_bitvec[word_ct - 1] & source2_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    target_bitvec[widx] = source1_bitvec[widx] & source2_bitvec[widx];
  }
#endif
}

void bitvec_andnot_copy(const uintptr_t* __restrict source_bitvec, const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec) {
  // target_bitvec := source_bitvec AND (~exclude_bitvec)
#ifdef __LP64__
  vul_t* target_bitvvec = (vul_t*)target_bitvec;
  const vul_t* source_bitvvec = (const vul_t*)source_bitvec;
  const vul_t* exclude_bitvvec = (const vul_t*)exclude_bitvec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  for (uintptr_t ulii = 0; ulii < full_vec_ct; ++ulii) {
    target_bitvvec[ulii] = source_bitvvec[ulii] & (~exclude_bitvvec[ulii]);
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = source_bitvec[base_idx] & (~exclude_bitvec[base_idx]);
    target_bitvec[base_idx + 1] = source_bitvec[base_idx + 1] & (~exclude_bitvec[base_idx + 1]);
  }
  #endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = source_bitvec[word_ct - 1] & (~exclude_bitvec[word_ct - 1]);
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    target_bitvec[widx] = source_bitvec[widx] & (~exclude_bitvec[widx]);
  }
#endif
}

void bitvec_or(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec) {
  // main_bitvec := main_bitvec OR arg_bitvec
#ifdef __LP64__
  vul_t* main_bitvvec_iter = (vul_t*)main_bitvec;
  const vul_t* arg_bitvvec_iter = (const vul_t*)arg_bitvec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] |= arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] |= arg_bitvec[base_idx + 1]
  }
  #endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] |= arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] |= arg_bitvec[widx];
  }
#endif
}

void bitvec_andnot2(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := (~main_bitvec) AND include_bitvec
  // this corresponds _mm_andnot() operand order
#ifdef __LP64__
  vul_t* main_bitvvec_iter = (vul_t*)main_bitvec;
  const vul_t* include_bitvvec_iter = (const vul_t*)include_bitvec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = (~(*main_bitvvec_iter)) & (*include_bitvvec_iter++);
    ++main_bitvvec_iter;
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] = (~main_bitvec[base_idx]) & include_bitvec[base_idx];
    main_bitvec[base_idx + 1] = (~main_bitvec[base_idx + 1]) & include_bitvec[base_idx + 1];
  }
  #endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] = (~main_bitvec[word_ct - 1]) & include_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] = (~main_bitvec[widx]) & include_bitvec[widx];
  }
#endif
}

void set_het_missing(uintptr_t word_ct, uintptr_t* genovec) {
  // 01 -> 11, nothing else changes
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  vul_t* geno_vvec_iter = (vul_t*)genovec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    const vul_t cur_geno_vword = *geno_vvec_iter;
    const vul_t cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  if (full_vec_ct & 2) {
    vul_t cur_geno_vword = *geno_vvec_iter;
    vul_t cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    vul_t cur_geno_vword = *geno_vvec_iter;
    vul_t cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    uintptr_t geno_word = genovec[base_idx];
    genovec[base_idx] = geno_word | ((geno_word & kMask5555) << 1);
    geno_word = genovec[base_idx + 1];
    genovec[base_idx + 1] = geno_word | ((geno_word & kMask5555) << 1);
  }
  #endif
  if (word_ct & 1) {
    const uintptr_t geno_word = genovec[word_ct - 1];
    genovec[word_ct - 1] = geno_word | ((geno_word & kMask5555) << 1);
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    genovec[widx] = geno_word | ((geno_word & kMask5555) << 1);
  }
#endif
}

void genoarr_to_nonmissing(const uintptr_t* genoarr, uint32_t sample_ct, uintptr_t* nonmissing_bitarr) {
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  const uintptr_t* genoarr_iter = genoarr;
  halfword_t* nonmissing_bitarr_iter = (halfword_t*)nonmissing_bitarr;
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    uintptr_t ww = ~(*genoarr_iter++);
    ww = (ww | (ww >> 1)) & kMask5555;
    *nonmissing_bitarr_iter++ = pack_word_to_halfword(ww);
  }
  // zero trailing bits up to word boundary, in a way that doesn't create
  // aliasing issues
  // (if zeroing is needed up to vector boundary, that's the caller's
  // responsibility)
  const uint32_t trail_ct = sample_ct % kBitsPerWordD2;
  if (trail_ct) {
    nonmissing_bitarr_iter[-1] &= (1U << trail_ct) - 1;
  }
  if (sample_ctl2 % 2) {
    *nonmissing_bitarr_iter = 0;
  }
}

uint32_t genoarr_count_missing_notsubset_unsafe(const uintptr_t* genoarr, const uintptr_t* exclude_mask, uint32_t sample_ct) {
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  const uintptr_t* genoarr_iter = genoarr;
  const halfword_t* exclude_alias_iter = (halfword_t*)exclude_mask;
  uint32_t missing_ct = 0;
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    uintptr_t ww = *genoarr_iter++;
    ww = ww & (ww >> 1);
    const uint32_t include_mask = ~(*exclude_alias_iter++);
    missing_ct += popcount01_long(ww & unpack_halfword_to_word(include_mask));
  }
  return missing_ct;
}


int32_t get_variant_uidx_without_htable(const char* idstr, char** variant_ids, const uintptr_t* variant_include, uint32_t variant_ct) {
  const uint32_t id_blen = strlen(idstr) + 1;
  uint32_t variant_uidx = 0;
  int32_t retval = -1;
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    if (!memcmp(idstr, variant_ids[variant_uidx], id_blen)) {
      if (retval != -1) {
	// duplicate
	return -2;
      }
      retval = (int32_t)variant_uidx;
    }
  }
  return retval;
}

// assumes 0b11 == missing
// retired in favor of genoarr_to_nonmissing followed by a for loop?
/*
void copy_when_nonmissing(const uintptr_t* loadbuf, const void* source, uintptr_t elem_size, uintptr_t unfiltered_sample_ct, uintptr_t missing_ct, void* dest) {
  // tried hardcoding elem_size == 4 and 8; that was only ~4% faster, so may
  // as well use a single general-purpose function.
  if (!missing_ct) {
    memcpy(dest, source, unfiltered_sample_ct * elem_size);
    return;
  }
  const uintptr_t* loadbuf_iter = loadbuf;
  const uintptr_t* loadbuf_end = &(loadbuf[QUATERCT_TO_WORDCT(unfiltered_sample_ct)]);
  const unsigned char* source_alias = (const unsigned char*)source;
  char* dest_iter = (char*)dest;
  uintptr_t copy_start_idx = 0;
  uintptr_t sample_idx_offset = 0;
  do {
    uintptr_t cur_word = *loadbuf_iter++;
    cur_word = cur_word & (cur_word >> 1) & kMask5555;
    while (cur_word) {
      const uintptr_t new_missing_idx = sample_idx_offset + (CTZLU(cur_word) / 2);
      if (new_missing_idx != copy_start_idx) {
        const uintptr_t diff = new_missing_idx - copy_start_idx;
	dest_iter = memcpya(dest_iter, &(source_alias[copy_start_idx * elem_size]), diff * elem_size);
      }
      copy_start_idx = new_missing_idx + 1;
      cur_word &= cur_word - 1;
    }
    sample_idx_offset += kBitsPerWordD2;
  } while (loadbuf_iter < loadbuf_end);
  const uintptr_t diff = unfiltered_sample_ct - copy_start_idx;
  if (diff) {
    memcpy(dest_iter, &(source_alias[copy_start_idx * elem_size]), diff * elem_size);
  }
}
*/


// MurmurHash3, from
// https://code.google.com/p/smhasher/source/browse/trunk/MurmurHash3.cpp
static inline uint32_t rotl32(uint32_t x, int8_t r) {
  return (x << r) | (x >> (32 - r));
}

static inline uint32_t getblock32(const uint32_t* p, int i) {
  return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

static inline uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

uint32_t murmurhash3_32(const void* key, uint32_t len) {
  const uint8_t* data = (const uint8_t*)key;
  const int32_t nblocks = len / 4;

  uint32_t h1 = 0;
  // uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  const uint32_t* blocks = (const uint32_t*)(data + nblocks*4);

  int32_t i;
  uint32_t k1;
  for(i = -nblocks; i; i++) {
      k1 = getblock32(blocks,i);

      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;
   
      h1 ^= k1;
      h1 = rotl32(h1,13);
      h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t* tail = (const uint8_t*)(data + nblocks*4);

  k1 = 0;

  switch(len & 3) {
    case 3:
      k1 ^= tail[2] << 16;
      // fall through
    case 2:
      k1 ^= tail[1] << 8;
      // fall through
    case 1:
      k1 ^= tail[0];
      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;
      h1 ^= k1;
  }

  //----------
  // finalization

  h1 ^= len;

  return fmix32(h1);
}


uint32_t is_composite6(uintptr_t num) {
  // assumes num is congruent to 1 or 5 mod 6.
  // can speed this up by ~50% by hardcoding avoidance of multiples of 5/7,
  // but this isn't currently a bottleneck so I'll keep this simple
  assert((num % 6 == 1) || (num % 6 == 5));
  uintptr_t divisor = 5;
  while (divisor * divisor <= num) {
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 2;
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 4;
  }
  return 0;
}

uintptr_t geqprime(uintptr_t floor) {
  // assumes floor is odd and greater than 1.  Returns 5 if floor = 3,
  // otherwise returns the first prime >= floor.
  assert((floor % 2) && (floor > 1));
  uintptr_t ulii = floor % 3;
  if (!ulii) {
    floor += 2;
  } else if (ulii == 1) {
    goto geqprime_1mod6;
  }
  while (is_composite6(floor)) {
    floor += 2;
  geqprime_1mod6:
    if (!is_composite6(floor)) {
      return floor;
    }
    floor += 4;
  }
  return floor;
}

uintptr_t leqprime(uintptr_t ceil) {
  // assumes ceil is odd and greater than 4.  Returns the first prime <= ceil.
  assert((ceil % 2) && (ceil > 4));
  uintptr_t ulii = ceil % 3;
  if (!ulii) {
    ceil -= 2;
  } else if (ulii == 2) {
    goto leqprime_5mod6;
  }
  while (is_composite6(ceil)) {
    ceil -= 2;
  leqprime_5mod6:
    if (!is_composite6(ceil)) {
      return ceil;
    }
    ceil -= 4;
  }
  return ceil;
}

boolerr_t htable_good_size_alloc(uint32_t item_ct, uintptr_t bytes_avail, uint32_t** htable_ptr, uint32_t* htable_size_ptr) {
  bytes_avail &= (~(kCacheline - k1LU));
  uint32_t htable_size = get_htable_fast_size(item_ct);
  if (htable_size > bytes_avail / sizeof(int32_t)) {
    if (!bytes_avail) {
      return 1;
    }
    htable_size = leqprime((bytes_avail / sizeof(int32_t)) - 1);
    if (htable_size < item_ct * 2) {
      return 1;
    }
  }
  *htable_ptr = (uint32_t*)bigstack_alloc_raw_rd(htable_size * sizeof(int32_t));
  *htable_size_ptr = htable_size;
  return 0;
}

uint32_t populate_strbox_htable(const char* strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable) {
  // may want subset_mask parameter later
  fill_uint_one(str_htable_size, str_htable);
  const char* strbox_iter = strbox;
  for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    const uint32_t slen = strlen(strbox_iter);
    uint32_t hashval = hashceil(strbox_iter, slen, str_htable_size);
    // previously used quadratic probing, but turns out that that isn't
    // meaningfully better than linear probing.
    // uint32_t next_incr = 1;
    while (1) {
      const uint32_t cur_htable_entry = str_htable[hashval];
      if (cur_htable_entry == 0xffffffffU) {
	str_htable[hashval] = str_idx;
	break;
      }
      if (!memcmp(strbox_iter, &(strbox[cur_htable_entry * max_str_blen]), slen + 1)) {
	// guaranteed to be positive
	return str_idx;
      }
      if (++hashval == str_htable_size) {
	hashval = 0;
      }
      /*
      // defend against overflow
      const uint32_t top_diff = str_htable_size - hashval;
      if (top_diff > next_incr) {
	hashval += next_incr;
      } else {
	hashval = next_incr - top_diff;
      }
      next_incr += 2;
      */
    }
    strbox_iter = &(strbox_iter[max_str_blen]);
  }
  return 0;
}

// could merge this with non-subset case, but this isn't much code
/*
uint32_t populate_strbox_subset_htable(const uintptr_t* __restrict subset_mask, const char* strbox, uintptr_t raw_str_ct, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable) {
  // may want subset_mask parameter later
  fill_uint_one(str_htable_size, str_htable);
  uintptr_t str_uidx = 0;
  for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx, ++str_uidx) {
    next_set_ul_unsafe_ck(subset_mask, &str_uidx);
    const char* cur_str = &(strbox[str_uidx * max_str_blen]);
    const uint32_t slen = strlen(cur_str);
    uint32_t hashval = hashceil(cur_str, slen, str_htable_size);
    while (1) {
      const uint32_t cur_htable_entry = str_htable[hashval];
      if (cur_htable_entry == 0xffffffffU) {
	str_htable[hashval] = str_uidx;
	break;
      }
      if (!memcmp(cur_str, &(strbox[cur_htable_entry * max_str_blen]), slen + 1)) {
	// guaranteed to be positive
	return str_uidx;
      }
      if (++hashval == str_htable_size) {
	hashval = 0;
      }
    }
  }
  return 0;
}
*/

uint32_t id_htable_find(const char* cur_id, char** item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size) {
  // returns 0xffffffffU on failure
  uint32_t hashval = hashceil(cur_id, cur_id_slen, id_htable_size);
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == 0xffffffffU) || (!strcmp(cur_id, item_ids[cur_htable_idval]))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

// assumes cur_id_slen < max_str_blen.
// requires cur_id to be null-terminated.
uint32_t strbox_htable_find(const char* cur_id, const char* strbox, const uint32_t* id_htable, uintptr_t max_str_blen, uint32_t cur_id_slen, uint32_t id_htable_size) {
  uint32_t hashval = hashceil(cur_id, cur_id_slen, id_htable_size);
  const uint32_t cur_id_blen = cur_id_slen + 1;
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == 0xffffffffU) || (!memcmp(cur_id, &(strbox[cur_htable_idval * max_str_blen]), cur_id_blen))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t variant_id_dupflag_htable_find(const char* idbuf, char** variant_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen) {
  // assumes duplicate variant IDs are flagged, but full variant_uidx linked
  // lists are not stored
  // idbuf does not need to be null-terminated (note that this is currently
  // achieved in a way that forces variant_ids[] entries to not be too close
  // to the end of bigstack, otherwise memcmp behavior is potentially
  // undefined)
  // returns 0xffffffffU on failure, value with bit 31 set on duplicate
  if (cur_id_slen > max_id_slen) {
    return 0xffffffffU;
  }
  uint32_t hashval = hashceil(idbuf, cur_id_slen, id_htable_size);
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == 0xffffffffU) || ((!memcmp(idbuf, variant_ids[cur_htable_idval & 0x7fffffff], cur_id_slen)) && (!variant_ids[cur_htable_idval & 0x7fffffff][cur_id_slen]))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t variant_id_dup_htable_find(const char* idbuf, char** variant_ids, const uint32_t* id_htable, const uint32_t* htable_dup_base, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen, uint32_t* llidx_ptr) {
  // Permits duplicate entries.  Similar to plink 1.9
  // extract_exclude_process_token().
  // - Returns 0xffffffffU on failure (llidx currently unset in that case),
  //   otherwise returns the index of the first match (which will have the
  //   highest index, due to how the linked list is constructed)
  // - Sets second_llidx to 0xffffffffU if not a duplicate, otherwise it's the
  //   position in htable_dup_base[] of the next {variant_uidx, next_llidx}
  //   linked list entry.
  // - idbuf does not need to be null-terminated.
  if (cur_id_slen > max_id_slen) {
    return 0xffffffffU;
  }
  uint32_t hashval = hashceil(idbuf, cur_id_slen, id_htable_size);
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    const uint32_t cur_dup = cur_htable_idval >> 31;
    uint32_t cur_llidx;
    uint32_t variant_uidx;
    if (cur_dup) {
      // 0xffffffffU empty-entry code has high bit set, so only need to check
      // here
      if (cur_htable_idval == 0xffffffffU) {
	return 0xffffffffU;
      }
      cur_llidx = cur_htable_idval << 1;
      variant_uidx = htable_dup_base[cur_llidx];
    } else {
      cur_llidx = 0xffffffffU;
      variant_uidx = cur_htable_idval;
    }
    const char* sptr = variant_ids[variant_uidx];
    if ((!memcmp(idbuf, sptr, cur_id_slen)) && (!sptr[cur_id_slen])) {
      *llidx_ptr = cur_llidx;
      return variant_uidx;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

char* scan_for_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_blen) {
  --id_ct;
  for (uintptr_t id_idx = 0; id_idx < id_ct; ++id_idx) {
    if (!strcmp(&(sorted_ids[id_idx * max_id_blen]), &(sorted_ids[(id_idx + 1) * max_id_blen]))) {
      return &(sorted_ids[id_idx * max_id_blen]);
    }
  }
  return nullptr;
}

uint32_t collapse_duplicate_ids(uintptr_t id_ct, uintptr_t max_id_blen, char* sorted_ids, uint32_t* id_starts) {
  // Collapses array of sorted IDs to remove duplicates, and writes
  // pre-collapse positions to id_starts (so e.g. duplication count of any
  // sample ID can be determined via subtraction) if it isn't nullptr.
  // Returns id_ct of collapsed array.
  if (!id_ct) {
    return 0;
  }
  uintptr_t read_idx = 1;
  uintptr_t write_idx;
  if (id_starts) {
    id_starts[0] = 0;
    for (; read_idx < id_ct; ++read_idx) {
      if (!strcmp(&(sorted_ids[(read_idx - 1) * max_id_blen]), &(sorted_ids[read_idx * max_id_blen]))) {
	break;
      }
      id_starts[read_idx] = read_idx;
    }
    write_idx = read_idx;
    while (++read_idx < id_ct) {
      if (strcmp(&(sorted_ids[(write_idx - 1) * max_id_blen]), &(sorted_ids[read_idx * max_id_blen]))) {
	strcpy(&(sorted_ids[write_idx * max_id_blen]), &(sorted_ids[read_idx * max_id_blen]));
	id_starts[write_idx++] = read_idx;
      }
    }
  } else {
    for (; read_idx < id_ct; ++read_idx) {
      if (!strcmp(&(sorted_ids[(read_idx - 1) * max_id_blen]), &(sorted_ids[read_idx * max_id_blen]))) {
	break;
      }
    }
    write_idx = read_idx;
    while (++read_idx < id_ct) {
      if (strcmp(&(sorted_ids[(write_idx - 1) * max_id_blen]), &(sorted_ids[read_idx * max_id_blen]))) {
	strcpy(&(sorted_ids[write_idx * max_id_blen]), &(sorted_ids[read_idx * max_id_blen]));
	++write_idx;
      }
    }
  }
  return write_idx;
}

pglerr_t copy_sort_strbox_subset_noalloc(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t collapse_idxs, uint32_t use_nsort, char* __restrict sorted_strbox, uint32_t* __restrict id_map) {
  // Stores a lexicographically sorted list of IDs in sorted_strbox and the raw
  // positions of the corresponding markers/samples in *id_map_ptr.  Does not
  // include excluded markers/samples in the list.
  // Assumes sorted_strbox and id_map have been allocated; use the
  // copy_sort_strbox_subset() wrapper if they haven't been.
  // Note that this DOES still perform a "stack" allocation (in the qsort_ext()
  // call).
  if (!str_ct) {
    return kPglRetSuccess;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
#ifdef __cplusplus
    if (max_str_blen <= 60) {
      uintptr_t wkspace_entry_blen = (max_str_blen > 36)? sizeof(Strbuf60_ui) : sizeof(Strbuf36_ui);
      char* sort_wkspace;
      if (bigstack_alloc_c(str_ct * wkspace_entry_blen, &sort_wkspace)) {
	goto copy_sort_strbox_subset_noalloc_ret_NOMEM;
      }
      uint32_t str_uidx = 0;
      const uint32_t wkspace_entry_blen_m4 = wkspace_entry_blen - 4;
      char* sort_wkspace_iter = sort_wkspace;
      for (uint32_t str_idx = 0; str_idx < str_ct; ++str_idx, ++str_uidx) {
	next_set_unsafe_ck(subset_mask, &str_uidx);
	strcpy(sort_wkspace_iter, &(orig_strbox[str_uidx * max_str_blen]));
	sort_wkspace_iter = &(sort_wkspace_iter[wkspace_entry_blen_m4]);
	if (collapse_idxs) {
	  *((uint32_t*)sort_wkspace_iter) = str_idx;
	} else {
	  *((uint32_t*)sort_wkspace_iter) = str_uidx;
	}
	sort_wkspace_iter = &(sort_wkspace_iter[sizeof(int32_t)]);
      }
      if (wkspace_entry_blen == 40) {
	sort_strbox_40b_finish(str_ct, max_str_blen, use_nsort, (Strbuf36_ui*)sort_wkspace, sorted_strbox, id_map);
      } else {
	sort_strbox_64b_finish(str_ct, max_str_blen, use_nsort, (Strbuf60_ui*)sort_wkspace, sorted_strbox, id_map);
      }
    } else {
#endif
      str_sort_indexed_deref_t* sort_wkspace = (str_sort_indexed_deref_t*)bigstack_alloc(str_ct * sizeof(str_sort_indexed_deref_t));
      if (!sort_wkspace) {
	goto copy_sort_strbox_subset_noalloc_ret_NOMEM;
      }
      uint32_t str_uidx = 0;
      for (uint32_t str_idx = 0; str_idx < str_ct; ++str_idx, ++str_uidx) {
	next_set_unsafe_ck(subset_mask, &str_uidx);
	sort_wkspace[str_idx].strptr = (const char*)(&(orig_strbox[str_uidx * max_str_blen]));
	if (collapse_idxs) {
	  sort_wkspace[str_idx].orig_idx = str_idx;
	} else {
	  sort_wkspace[str_idx].orig_idx = str_uidx;
	}
      }
      if (!use_nsort) {
#ifdef __cplusplus
        std::sort(sort_wkspace, &(sort_wkspace[str_ct]));
#else
        qsort(sort_wkspace, str_ct, sizeof(str_sort_indexed_deref_t), strcmp_deref);
#endif
      } else {
#ifdef __cplusplus
	str_nsort_indexed_deref_t* wkspace_alias = (str_nsort_indexed_deref_t*)sort_wkspace;
	std::sort(wkspace_alias, &(wkspace_alias[str_ct]));
#else
	qsort(sort_wkspace, str_ct, sizeof(str_sort_indexed_deref_t), strcmp_natural_deref);
#endif
      }
      for (uintptr_t str_idx = 0; str_idx < str_ct; ++str_idx) {
	strcpy(&(sorted_strbox[str_idx * max_str_blen]), sort_wkspace[str_idx].strptr);
	id_map[str_idx] = sort_wkspace[str_idx].orig_idx;
      }
#ifdef __cplusplus
    }
#endif
    if (!allow_dups) {
      char* dup_id = scan_for_duplicate_ids(sorted_strbox, str_ct, max_str_blen);
      if (dup_id) {
	char* tptr = dup_id;
	while (1) {
	  tptr = strchr(tptr, '\t');
	  if (!tptr) {
	    break;
	  }
	  *tptr++ = ' ';
	}
	LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", dup_id);
	goto copy_sort_strbox_subset_noalloc_ret_MALFORMED_INPUT;
      }
    }
  }
  while (0) {
  copy_sort_strbox_subset_noalloc_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  copy_sort_strbox_subset_noalloc_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t copy_sort_strbox_subset(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t collapse_idxs, uint32_t use_nsort, char** sorted_strbox_ptr, uint32_t** id_map_ptr) {
  // id_map on bottom because --indiv-sort frees *sorted_strbox_ptr
  if (bigstack_alloc_ui(str_ct, id_map_ptr) ||
      bigstack_alloc_c(str_ct * max_str_blen, sorted_strbox_ptr)) {
    return kPglRetNomem;
  }
  return copy_sort_strbox_subset_noalloc(subset_mask, orig_strbox, str_ct, max_str_blen, allow_dups, collapse_idxs, use_nsort, *sorted_strbox_ptr, *id_map_ptr);
}

int32_t bsearch_str(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx) {
  // does not assume null-terminated idbuf, or nonempty array.
  if (cur_id_slen >= max_id_blen) {
    return -1;
  }
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    const int32_t ii = memcmp(idbuf, &(sorted_strbox[mid_idx * max_id_blen]), cur_id_slen);
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if ((ii < 0) || sorted_strbox[mid_idx * max_id_blen + cur_id_slen]) {
      end_idx = mid_idx;
    } else {
      return ((uint32_t)mid_idx);
    }
  }
  return -1;
}

int32_t bsearch_str_natural(const char* idbuf, const char* sorted_strbox, uintptr_t max_id_blen, uintptr_t end_idx) {
  // unlike bsearch_str(), caller is responsible for slen >= max_id_blen check
  // if appropriate here
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    const int32_t ii = strcmp_natural(idbuf, &(sorted_strbox[mid_idx * max_id_blen]));
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if (ii < 0) {
      end_idx = mid_idx;
    } else {
      return ((uint32_t)mid_idx);
    }
  }
  return -1;
}

uintptr_t bsearch_str_lb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx) {
  // returns number of elements in sorted_strbox[] less than idbuf.
  if (cur_id_slen > max_id_blen) {
    cur_id_slen = max_id_blen;
  }
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    if (memcmp(idbuf, &(sorted_strbox[mid_idx * max_id_blen]), cur_id_slen) > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}


uint32_t sid_col_required(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier) {
  // note that MAYBESID and SID can both be set
  if (maybe_modifier & 2) {
    return 1;
  }
  if (sids && (maybe_modifier & 1)) {
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      if (memcmp(&(sids[sample_uidx * max_sid_blen]), "0", 2)) {
	return 1;
      }
    }
  }
  return 0;
}

// sample_augid_map_ptr == nullptr ok
pglerr_t augid_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t** sample_augid_map_ptr, char** sample_augids_ptr, uintptr_t* max_sample_augid_blen_ptr) {
  if (!sids) {
    max_sid_blen = 2;
  }
  const uintptr_t max_sample_augid_blen = max_sample_id_blen + max_sid_blen;
  *max_sample_augid_blen_ptr = max_sample_augid_blen;
  uint32_t* sample_augid_map = nullptr;
  if (sample_augid_map_ptr) {
    if (bigstack_alloc_ui(sample_ct, sample_augid_map_ptr)) {
      return kPglRetNomem;
    }
    sample_augid_map = *sample_augid_map_ptr;
  }
  if (bigstack_alloc_c(max_sample_augid_blen * sample_ct, sample_augids_ptr)) {
    return kPglRetNomem;
  }
  char* sample_augids_iter = *sample_augids_ptr;
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include, &sample_uidx);
    char* write_iter = strcpyax(sample_augids_iter, &(sample_ids[sample_uidx * max_sample_id_blen]), '\t');
    if (sids) {
      strcpy(write_iter, &(sids[sample_uidx * max_sid_blen]));
    } else {
      strcpy(write_iter, "0");
    }
    sample_augids_iter = &(sample_augids_iter[max_sample_augid_blen]);
    if (sample_augid_map) {
      sample_augid_map[sample_idx] = sample_uidx;
    }
  }
  return kPglRetSuccess;
}

pglerr_t sorted_xidbox_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, xid_mode_t xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr) {
  if (!(xid_mode & kfXidModeFlagSid)) {
    // two fields
    *max_xid_blen_ptr = max_sample_id_blen;
    return copy_sort_strbox_subset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 0, 0, use_nsort, sorted_xidbox_ptr, xid_map_ptr);
  }
  // three fields
  if (augid_init_alloc(sample_include, sample_ids, sids, sample_ct, max_sample_id_blen, max_sid_blen, xid_map_ptr, sorted_xidbox_ptr, max_xid_blen_ptr)) {
    return kPglRetNomem;
  }
  if (sort_strbox_indexed(sample_ct, *max_xid_blen_ptr, use_nsort, *sorted_xidbox_ptr, *xid_map_ptr)) {
    return kPglRetNomem;
  }
  char* dup_id = scan_for_duplicate_ids(*sorted_xidbox_ptr, sample_ct, *max_xid_blen_ptr);
  if (dup_id) {
    char* tptr = (char*)rawmemchr(dup_id, '\t');
    *tptr = ' ';
    tptr = (char*)rawmemchr(&(tptr[1]), '\t');
    *tptr = ' ';
    LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", dup_id);
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

boolerr_t sorted_xidbox_read_find(const char* __restrict sorted_xidbox, const uint32_t* __restrict xid_map, uintptr_t max_xid_blen, uintptr_t end_idx, uint32_t comma_delim, xid_mode_t xid_mode, char** read_pp, uint32_t* sample_uidx_ptr, char* __restrict idbuf) {
  // idbuf = workspace
  // sorted_xidbox = packed, sorted list of ID strings to search over.
  //
  // input *read_pp must point to beginning of FID; this is a change from plink
  // 1.9.
  //
  // *read_pp is now set to point to the end of the last parsed token instead
  // of the beginning of the next; this is another change from plink 1.9.
  //
  // returns 1 on missing token *or* if the sample ID is not present.  cases
  // can be distinguished by checking whether *read_pp == nullptr.
  char* first_token_start = *read_pp;
  uintptr_t blen_sid = 0;
  char* token_iter;
  char* iid_ptr;
  char* sid_ptr = nullptr;
  uintptr_t slen_fid;
  uintptr_t slen_iid;
  if (comma_delim) {
    token_iter = first_token_start;
    unsigned char ucc = (unsigned char)(*token_iter);
    while (ucc != ',') {
      if (ucc < 32) {
	if (!(xid_mode & kfXidModeFlagOneTokenOk)) {
	  *read_pp = nullptr;
	  return 1;
	}
	slen_fid = (uintptr_t)(token_iter - first_token_start);
	goto sorted_xidbox_read_find_comma_single_token;
      }
      ucc = (unsigned char)(*(++token_iter));
    }
    slen_fid = (uintptr_t)(token_iter - first_token_start);
    if (xid_mode & kfXidModeFlagNeverFid) {
    sorted_xidbox_read_find_comma_single_token:
      iid_ptr = first_token_start;
      slen_iid = slen_fid;
    } else {
      do {
	ucc = (unsigned char)(*(++token_iter));
      } while ((ucc == ' ') || (ucc == '\t'));
      iid_ptr = token_iter;
      while ((ucc >= 32) && (ucc != ',')) {
	ucc = (unsigned char)(*(++token_iter));
      }
      slen_iid = (uintptr_t)(token_iter - iid_ptr);
    }
    // token_iter now points to comma/eoln at end of IID
    if (xid_mode & kfXidModeFlagSid) {
      if (*token_iter != ',') {
	return 1;
      }
      do {
	ucc = (unsigned char)(*(++token_iter));
      } while ((ucc == ' ') || (ucc == '\t'));
      sid_ptr = token_iter;
      while ((ucc >= 32) && (ucc != ',')) {
	ucc = (unsigned char)(*(++token_iter));
      }
      blen_sid = 1 + (uintptr_t)(token_iter - sid_ptr);
      if (token_iter == sid_ptr) {
	// special case: treat missing SID as '0'
	blen_sid = 2;
	// const_cast, since token_endnn doesn't return const pointer
        // function is too long for me to be comfortable just turning off
        // -Wcast-qual...
        sid_ptr = (char*)((uintptr_t)(&(g_one_char_strs[96])));
      }
    }
  } else {
    assert(!is_eoln_kns(*first_token_start));
    token_iter = token_endnn(first_token_start);
    slen_fid = (uintptr_t)(token_iter - first_token_start);
    if (xid_mode & kfXidModeFlagNeverFid) {
    sorted_xidbox_read_find_space_single_token:
      iid_ptr = first_token_start;
      slen_iid = slen_fid;
    } else {
      token_iter = skip_initial_spaces(token_iter);
      if (is_eoln_kns(*token_iter)) {
	if (!(xid_mode & kfXidModeFlagOneTokenOk)) {
	  *read_pp = nullptr;
	  return 1;
	}
	// need to backtrack
	token_iter = &(first_token_start[slen_fid]);
	goto sorted_xidbox_read_find_space_single_token;
      }
      iid_ptr = token_iter;
      token_iter = token_endnn(iid_ptr);
      slen_iid = (uintptr_t)(token_iter - iid_ptr);
    }
    // token_iter now points to space/eoln at end of IID
    if (xid_mode & kfXidModeFlagSid) {
      token_iter = skip_initial_spaces(token_iter);
      if (is_eoln_kns(*token_iter)) {
	*read_pp = nullptr;
	return 1;
      }
      sid_ptr = token_iter;
      token_iter = token_endnn(sid_ptr);
      blen_sid = 1 + (uintptr_t)(token_iter - sid_ptr);
    }
  }
  *read_pp = token_iter;
  uintptr_t slen_final = slen_fid + slen_iid + blen_sid + 1;
  if (slen_final >= max_xid_blen) {
    // avoid buffer overflow
    return 1;
  }
  char* idbuf_end = memcpya(memcpyax(idbuf, first_token_start, slen_fid, '\t'), iid_ptr, slen_iid);
  if (blen_sid) {
    *idbuf_end++ = '\t';
    memcpy(idbuf_end, sid_ptr, blen_sid - 1);
  }
  return sorted_idbox_find(idbuf, sorted_xidbox, xid_map, slen_final, max_xid_blen, end_idx, sample_uidx_ptr);
}

void init_range_list(range_list_t* range_list_ptr) {
  range_list_ptr->names = nullptr;
  range_list_ptr->starts_range = nullptr;
  range_list_ptr->name_ct = 0;
  range_list_ptr->name_max_blen = 0;
}

void cleanup_range_list(range_list_t* range_list_ptr) {
  free_cond(range_list_ptr->names);
  // starts_range now uses same allocation
}

boolerr_t numeric_range_list_to_bitarr(const range_list_t* range_list_ptr, uint32_t bitarr_size, uint32_t offset, uint32_t ignore_overflow, uintptr_t* bitarr) {
  // bitarr assumed to be initialized (but not necessarily zero-initialized)
  const char* names = range_list_ptr->names;
  const unsigned char* starts_range = range_list_ptr->starts_range;
  const uint32_t name_ct = range_list_ptr->name_ct;
  const uint32_t name_max_blen = range_list_ptr->name_max_blen;
  const uint32_t idx_max = bitarr_size + offset;
  for (uint32_t name_idx = 0; name_idx < name_ct; ++name_idx) {
    uint32_t idx1;
    if (scan_uint_capped(&(names[name_idx * name_max_blen]), idx_max, &idx1)) {
      if (ignore_overflow) {
	continue;
      }
      return 1;
    }
    if (starts_range[name_idx]) {
      ++name_idx;
      uint32_t idx2;
      if (scan_uint_capped(&(names[name_idx * name_max_blen]), idx_max, &idx2)) {
	if (!ignore_overflow) {
	  return 1;
	}
        idx2 = idx_max - 1;
      }
      fill_bits_nz(idx1 - offset, (idx2 - offset) + 1, bitarr);
    } else {
      set_bit(idx1 - offset, bitarr);
    }
  }
  return 0;
}

pglerr_t string_range_list_to_bitarr(char* header_line, const range_list_t* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t* bitarr, int32_t* __restrict seen_idxs) {
  // bitarr assumed to be zero-initialized
  // if fixed_len is zero, header_line is assumed to be a list of
  // space-delimited unequal-length names
  assert(token_ct);
  assert(!popcount_longs(bitarr, BITCT_TO_WORDCT(token_ct)));
  pglerr_t reterr = kPglRetSuccess;
  {
    char* header_line_iter = header_line;
    const uintptr_t name_ct = range_list_ptr->name_ct;
    const uintptr_t max_id_blen = range_list_ptr->name_max_blen;
    uint32_t item_idx = 0;
    while (1) {
      char* token_end = comma_or_space_token_end(header_line_iter, comma_delim);
      uint32_t cmdline_pos;
      if (!sorted_idbox_find(header_line_iter, sorted_ids, id_map, (uintptr_t)(token_end - header_line_iter), max_id_blen, name_ct, &cmdline_pos)) {
	if (seen_idxs[cmdline_pos] != -1) {
	  sprintf(g_logbuf, "Error: Duplicate --%s token in %s.\n", range_list_flag, file_descrip);
	  goto string_range_list_to_bitarr_ret_MALFORMED_INPUT_2;
	}
	seen_idxs[cmdline_pos] = item_idx;
	if (cmdline_pos && range_list_ptr->starts_range[cmdline_pos - 1]) {
	  if (seen_idxs[cmdline_pos - 1] == -1) {
	    LOGPREPRINTFWW("Error: Second element of --%s range appears before first element in %s.\n", range_list_flag, file_descrip);
	    goto string_range_list_to_bitarr_ret_INVALID_CMDLINE_2;
	  }
	  fill_bits_nz(seen_idxs[cmdline_pos - 1], item_idx + 1, bitarr);
	} else if (!(range_list_ptr->starts_range[cmdline_pos])) {
	  SET_BIT(item_idx, bitarr);
	}
      }
      if (++item_idx == token_ct) {
	break;
      }
      if (fixed_len) {
	header_line_iter = &(header_line_iter[fixed_len]);
      } else {
	header_line_iter = skip_initial_spaces(&(token_end[1]));
      }
    }
    for (uint32_t cmdline_pos = 0; cmdline_pos < name_ct; ++cmdline_pos) {
      if (seen_idxs[cmdline_pos] == -1) {
	goto string_range_list_to_bitarr_ret_INVALID_CMDLINE_3;
      }
    }
  }
  while (0) {
  string_range_list_to_bitarr_ret_INVALID_CMDLINE_3:
    sprintf(g_logbuf, "Error: Missing --%s token in %s.\n", range_list_flag, file_descrip);
  string_range_list_to_bitarr_ret_INVALID_CMDLINE_2:
    logerrprintb();
    reterr = kPglRetInvalidCmdline;
    break;
  string_range_list_to_bitarr_ret_MALFORMED_INPUT_2:
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  }
  return reterr;
}

pglerr_t string_range_list_to_bitarr_alloc(char* header_line, const range_list_t* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t** bitarr_ptr) {
  // wrapper for string_range_list_to_bitarr which allocates the bitfield and
  // temporary buffers on the heap
  uintptr_t token_ctl = BITCT_TO_WORDCT(token_ct);
  uintptr_t name_ct = range_list_ptr->name_ct;
  int32_t* seen_idxs;
  char* sorted_ids;
  uint32_t* id_map;
  if (bigstack_calloc_ul(token_ctl, bitarr_ptr) ||
      bigstack_alloc_i(name_ct, &seen_idxs)) {
    return kPglRetNomem;
  }
  // kludge to use copy_sort_strbox_subset()
  fill_all_bits(name_ct, (uintptr_t*)seen_idxs);
  if (copy_sort_strbox_subset((uintptr_t*)seen_idxs, range_list_ptr->names, name_ct, range_list_ptr->name_max_blen, 0, 0, 0, &sorted_ids, &id_map)) {
    return kPglRetNomem;
  }
  fill_int_one(name_ct, seen_idxs);
  pglerr_t reterr = string_range_list_to_bitarr(header_line, range_list_ptr, sorted_ids, id_map, range_list_flag, file_descrip, token_ct, fixed_len, comma_delim, *bitarr_ptr, seen_idxs);
  bigstack_reset(seen_idxs);
  return reterr;
}


const char g_xymt_log_names[][5] = {"chrX", "chrY", "XY", "chrM", "PAR1", "PAR2"};

static_assert(!(kChrRawEnd % kBytesPerVec), "kChrRawEnd must be a multiple of kBytesPerVec.");
pglerr_t init_chr_info(chr_info_t* cip) {
  // "constructor".  initializes with maximum capacity.  doesn't use bigstack.
  // chr_mask, haploid_mask: bits
  // chr_file_order, chr_idx_to_foidx: int32s
  // chr_fo_vidx_start: int32s, with an extra trailing element
  // nonstd_names: intptr_ts
  // nonstd_id_htable: kChrHtableSize int32s

  // this assumes kChrRawEnd is divisible by kBytesPerVec
  const uintptr_t vecs_required = 2 * BITCT_TO_VECCT(kChrRawEnd) + 3 * (kChrRawEnd / kInt32PerVec) + 1 + (kChrRawEnd / kWordsPerVec) + INT32CT_TO_VECCT(kChrHtableSize);

  // needed for proper cleanup
  cip->name_ct = 0;
  cip->incl_excl_name_stack = nullptr;
  if (vecaligned_malloc(vecs_required * kBytesPerVec, &(cip->chr_mask))) {
    return kPglRetNomem;
  }
  uintptr_t* alloc_iter = &(cip->chr_mask[BITCT_TO_VECCT(kChrRawEnd) * kWordsPerVec]);
  cip->haploid_mask = alloc_iter;
  alloc_iter = &(alloc_iter[BITCT_TO_VECCT(kChrRawEnd) * kWordsPerVec]);
  cip->chr_file_order = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[(kChrRawEnd / kInt32PerVec) * kWordsPerVec]);
  cip->chr_fo_vidx_start = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[((kChrRawEnd / kInt32PerVec) + 1) * kWordsPerVec]);
  cip->chr_idx_to_foidx = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[(kChrRawEnd / kInt32PerVec) * kWordsPerVec]);
  cip->nonstd_names = (char**)alloc_iter;
  alloc_iter = &(alloc_iter[kChrRawEnd]);
  cip->nonstd_id_htable = (uint32_t*)alloc_iter;
  // alloc_iter = &(alloc_iter[((kChrHtableSize + (kInt32PerVec - 1)) / kInt32PerVec) * kWordsPerVec]);
  // fill_uint_one(kChrHtableSize, cip->nonstd_id_htable);

  fill_ulong_zero(kChrMaskWords, cip->chr_mask);
  fill_ulong_zero(kChrExcludeWords, cip->chr_exclude);

  // this is a change from plink 1.x.  MT > M since the former matches Ensembl,
  // while the latter doesn't match any major resource.  no "chr" to reduce
  // file sizes and reduce the impact of this change.
  cip->output_encoding = kfChrOutputMT;
  
  cip->zero_extra_chrs = 0;
  cip->is_include_stack = 0;
  cip->chrset_source = kChrsetSourceDefault;
  cip->autosome_ct = 22;
  for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
    cip->xymt_codes[xymt_idx] = 23 + xymt_idx;
  }
  cip->haploid_mask[0] = 0x1800000;
  fill_ulong_zero(kChrMaskWords - 1, &(cip->haploid_mask[1]));
  return kPglRetSuccess;
}

// explicit plink 1.07 species (now initialized by command line parser):
// human: 22, X, Y, XY, MT, PAR1, PAR2 (PAR1/PAR2 added, XY deprecated in plink
//   2.0)
// cow: 29, X, Y, MT
// dog: 38, X, Y, XY, MT
// horse: 31, X, Y
// mouse: 19, X, Y
// rice: 12
// sheep: 26, X, Y

// must be safe to call this twice.
void finalize_chrset(misc_flags_t misc_flags, chr_info_t* cip) {
  uint32_t autosome_ct = cip->autosome_ct;
  uint32_t max_code = autosome_ct;
  for (uint32_t xymt_idx_p1 = kChrOffsetCt; xymt_idx_p1; --xymt_idx_p1) {
    if (cip->xymt_codes[xymt_idx_p1 - 1] >= 0) {
      max_code = autosome_ct + xymt_idx_p1;
      break;
    }
  }
  
  // could initialize haploid_mask bits (after the first) here, instead of
  // earlier...
  
  cip->max_numeric_code = MINV(max_code, autosome_ct + 4);
  cip->max_code = max_code;
  uintptr_t* chr_mask = cip->chr_mask;
  uintptr_t last_chr_mask_word = chr_mask[kChrMaskWords - 1];
  int32_t* xymt_codes = cip->xymt_codes;
  if (last_chr_mask_word) {
    // avoids repeating some work if this is called twice
    chr_mask[kChrMaskWords - 1] = 0;

    uint32_t xymt_include = last_chr_mask_word >> (kBitsPerWord - kChrOffsetCt);
    do {
      const uint32_t xymt_idx = __builtin_ctz(xymt_include);
      const int32_t cur_chr_code = xymt_codes[xymt_idx];
      if (cur_chr_code >= 0) {
	set_bit(cur_chr_code, chr_mask);
      }
      xymt_include &= xymt_include - 1;
    } while (xymt_include);
  } else if (are_all_words_zero(chr_mask, kChrExcludeWords) && (!cip->is_include_stack)) {
    // init_default_chr_mask()
    fill_all_bits(cip->autosome_ct + 1, chr_mask);
    for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
      const int32_t cur_chr_code = cip->xymt_codes[xymt_idx];
      if (cur_chr_code >= 0) {
	set_bit(cur_chr_code, chr_mask);
      }
    }
  } else if (misc_flags & (kfMiscAutosomePar | kfMiscAutosomeOnly)) {
    fill_bits_nz(1, cip->autosome_ct + 1, chr_mask);
    clear_bits_nz(cip->autosome_ct + 1, kChrExcludeWords * kBitsPerWord, chr_mask);
    if (misc_flags & kfMiscAutosomePar) {
      int32_t par_chr_code = cip->xymt_codes[kChrOffsetXY];
      if (par_chr_code >= 0) {
	set_bit(par_chr_code, chr_mask);
      }
      par_chr_code = cip->xymt_codes[kChrOffsetPAR1];
      if (par_chr_code >= 0) {
	set_bit(par_chr_code, chr_mask);
      }
      par_chr_code = cip->xymt_codes[kChrOffsetPAR2];
      if (par_chr_code >= 0) {
	set_bit(par_chr_code, chr_mask);
      }
    }
  }
  
  uintptr_t* chr_exclude = cip->chr_exclude;
  uintptr_t last_chr_exclude_word = chr_exclude[kChrExcludeWords - 1];
  uint32_t xymt_exclude = last_chr_exclude_word >> (kBitsPerWord - kChrOffsetCt);
  last_chr_exclude_word &= (k1LU << (kBitsPerWord - kChrOffsetCt)) - k1LU;
  for (uint32_t widx = 0; widx < kChrExcludeWords - 1; ++widx) {
    chr_mask[widx] &= ~chr_exclude[widx];
  }
  chr_mask[kChrExcludeWords - 1] &= ~last_chr_exclude_word;
  if (xymt_exclude) {
    do {
      const uint32_t xymt_idx = __builtin_ctz(xymt_exclude);
      const int32_t cur_chr_code = xymt_codes[xymt_idx];
      if (cur_chr_code >= 0) {
	clear_bit(cur_chr_code, chr_mask);
      }
      xymt_exclude &= xymt_exclude - 1;
    } while (xymt_exclude);
  }
  fill_uint_one(max_code + 1, cip->chr_idx_to_foidx);
}

void forget_extra_chr_names(uint32_t reinitialize, chr_info_t* cip) {
  const uint32_t name_ct = cip->name_ct;
  if (name_ct) {
    char** nonstd_names = cip->nonstd_names;
    const uint32_t chr_idx_last = cip->max_code + name_ct;
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx <= chr_idx_last; ++chr_idx) {
      free(nonstd_names[chr_idx]);
      nonstd_names[chr_idx] = nullptr;
    }
    if (reinitialize) {
      // fill_uint_one(kChrHtableSize, cip->nonstd_id_htable);
      cip->name_ct = 0;
    }
  }
}

// not currently called.  might want to do so in the future.
pglerr_t finalize_chr_info(chr_info_t* cip) {
  const uint32_t chr_ct = cip->chr_ct;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = cip->max_code + 1 + name_ct;
  const uint32_t chr_code_bitvec_ct = BITCT_TO_VECCT(chr_code_end);
  const uint32_t chr_ct_int32vec_ct = INT32CT_TO_VECCT(chr_ct);
  const uint32_t chr_ct_p1_int32vec_ct = 1 + (chr_ct / kInt32PerVec);
  const uint32_t chr_code_end_int32vec_ct = INT32CT_TO_VECCT(chr_code_end);
  const uint32_t chr_code_end_wordvec_ct = WORDCT_TO_VECCT(chr_code_end);
  uint32_t final_vecs_required = 2 * chr_code_bitvec_ct + chr_ct_int32vec_ct + chr_ct_p1_int32vec_ct + chr_code_end_int32vec_ct;
  if (name_ct) {
    final_vecs_required += chr_code_end_wordvec_ct + INT32CT_TO_VECCT(kChrHtableSize);
  }
  uintptr_t* new_alloc;
  if (vecaligned_malloc(final_vecs_required * kBytesPerVec, &new_alloc)) {
    return kPglRetNomem;
  }
  uintptr_t* old_alloc = cip->chr_mask;
  uintptr_t* new_alloc_iter = new_alloc;

  memcpy(new_alloc_iter, cip->chr_mask, chr_code_bitvec_ct * kBytesPerVec);
  cip->chr_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_code_bitvec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->haploid_mask, chr_code_bitvec_ct * kBytesPerVec);
  cip->haploid_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_code_bitvec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_file_order, chr_ct_int32vec_ct * kBytesPerVec);
  cip->chr_file_order = (uint32_t*)new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_ct_int32vec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_fo_vidx_start, chr_ct_p1_int32vec_ct * kBytesPerVec);
  cip->chr_fo_vidx_start = (uint32_t*)new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_ct_p1_int32vec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_idx_to_foidx, chr_code_end_int32vec_ct * kBytesPerVec);
  cip->chr_idx_to_foidx = (uint32_t*)new_alloc_iter;

  if (!name_ct) {
    cip->nonstd_names = nullptr;
    cip->nonstd_id_htable = nullptr;
  } else {
    new_alloc_iter = &(new_alloc_iter[chr_code_end_int32vec_ct * kWordsPerVec]);

    memcpy(new_alloc_iter, cip->nonstd_names, chr_code_end_wordvec_ct * kBytesPerVec);
    cip->nonstd_names = (char**)new_alloc_iter;
    new_alloc_iter = &(new_alloc_iter[chr_code_end_wordvec_ct * kWordsPerVec]);

    memcpy(new_alloc_iter, cip->nonstd_id_htable, kChrHtableSize * sizeof(int32_t));
    cip->nonstd_id_htable = (uint32_t*)new_alloc_iter;
  }
  vecaligned_free(old_alloc);
  return kPglRetSuccess;
}

void cleanup_chr_info(chr_info_t* cip) {
  if (cip->chr_mask) {
    forget_extra_chr_names(0, cip);
    vecaligned_free(cip->chr_mask);
    cip->chr_mask = nullptr;
  }
  ll_str_t* llstr_ptr = cip->incl_excl_name_stack;
  while (llstr_ptr) {
    ll_str_t* next_ptr = llstr_ptr->next;
    free(llstr_ptr);
    llstr_ptr = next_ptr;
  }
  cip->incl_excl_name_stack = nullptr;
}

char* chr_name_std(const chr_info_t* cip, uint32_t chr_idx, char* buf) {
  if (chr_idx > cip->max_numeric_code) {
    // this is encoding-independent.  users who require all numbers should use
    // 25 == XY instead.
    // this code will probably need to be changed later if we add more standard
    // nonnumeric codes.
    memcpyl3(buf, "PAR");
    buf[3] = '0' + (((int32_t)chr_idx) - cip->max_numeric_code);
    return &(buf[4]);
  }
  const uint32_t output_encoding = cip->output_encoding;
  if (output_encoding & (kfChrOutputPrefix | kfChrOutput0M)) {
    if (output_encoding == kfChrOutput0M) {
      // force two chars
      if (chr_idx <= cip->autosome_ct) {
	buf = memcpya(buf, &(kDigitPair[chr_idx]), 2);
      } else if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetY]) {
	buf = strcpya(buf, "XY");
      } else {
	*buf++ = '0';
	if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetX]) {
	  *buf++ = 'X';
	} else {
	  // assumes only X/Y/XY/MT defined
	  *buf++ = ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetY])? 'Y' : 'M';
	}
      }
      return buf;
    }
    buf = memcpyl3a(buf, "chr");
  }
  if ((!(output_encoding & (kfChrOutputM | kfChrOutputMT))) || (chr_idx <= cip->autosome_ct)) {
    return uint32toa(chr_idx, buf);
  }
  if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetX]) {
    *buf++ = 'X';
  } else if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetY]) {
    *buf++ = 'Y';
  } else if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetXY]) {
    buf = strcpya(buf, "XY");
  } else {
    *buf++ = 'M';
    if (output_encoding & kfChrOutputMT) {
      *buf++ = 'T';
    }
  }
  return buf;
}

char* chr_name_write(const chr_info_t* cip, uint32_t chr_idx, char* buf) {
  // assumes chr_idx is valid
  if (!chr_idx) {
    *buf++ = '0';
    return buf;
  }
  if (chr_idx <= cip->max_code) {
    return chr_name_std(cip, chr_idx, buf);
  }
  if (cip->zero_extra_chrs) {
    *buf++ = '0';
    return buf;
  }
  return strcpya(buf, cip->nonstd_names[chr_idx]);
}

uint32_t get_max_chr_slen(const chr_info_t* cip) {
  // does not include trailing null
  // can be overestimate
  // if more functions start calling this, it should just be built into
  // load_bim() instead
  if (cip->zero_extra_chrs) {
    return 3 + kMaxChrTextnum;
  }
  const uint32_t chr_ct = cip->chr_ct;
  const uint32_t max_code = cip->max_code;
  uint32_t max_chr_slen = 3 + kMaxChrTextnum;
  for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    if (!is_set(cip->chr_mask, chr_idx)) {
      continue;
    }
    if (chr_idx > max_code) {
      const uint32_t name_slen = strlen(cip->nonstd_names[chr_idx]);
      if (name_slen > max_chr_slen) {
	max_chr_slen = name_slen;
      }
    }
  }
  return max_chr_slen;
}

uint32_t haploid_chr_present(const chr_info_t* cip) {
  const uintptr_t* chr_mask = cip->chr_mask;
  const uintptr_t* haploid_mask = cip->haploid_mask;
  // since we don't load haploid vs. diploid info from ##contig header lines,
  // this is sufficient
  for (uint32_t widx = 0; widx < kChrExcludeWords; ++widx) {
    if (chr_mask[widx] & haploid_mask[widx]) {
      return 1;
    }
  }
  return 0;
}

static inline int32_t single_cap_letter_chrom(uint32_t cap_letter) {
  if (cap_letter == 'X') {
    return kChrRawX;
  }
  if (cap_letter == 'Y') {
    return kChrRawY;
  }
  if (cap_letter == 'M') {
    return kChrRawMT;
  }
  return -1;
}

static_assert(kMaxChrTextnumSlen == 2, "get_chr_code_raw() must be updated.");
int32_t get_chr_code_raw(const char* sptr) {
  // any character <= ' ' is considered a terminator
  // note that char arithmetic tends to be compiled to int32 operations, so we
  // mostly work with ints here
  uint32_t first_char_code = (unsigned char)sptr[0];
  uint32_t first_char_toi;
  if (first_char_code < 58) {
  get_chr_code_raw_digits:
    first_char_toi = first_char_code - '0';
    if (first_char_toi < 10) {
      const uint32_t second_char_code = (unsigned char)sptr[1];
      if (second_char_code <= ' ') {
	return first_char_toi;
      }
      if (((unsigned char)sptr[2]) <= ' ') {
	const uint32_t second_char_toi = second_char_code - '0';
	if (second_char_toi < 10) {
	  return first_char_toi * 10 + second_char_toi;
	}
	if (!first_char_toi) {
	  // accept '0X', '0Y', '0M' emitted by Oxford software
	  return single_cap_letter_chrom(second_char_code & 0xdf);
	}
      }
    }
    return -1;
  }
  first_char_code &= 0xdf;
  uint32_t second_char_code = (unsigned char)sptr[1];
  if (first_char_code == 'P') {
    // chrPAR1 *not* supported; has to be PAR1 by itself.
    // can't do uint16_t compare of multiple characters, since we could be
    // dealing with a length-1 null-terminated string; that IS faster when it's
    // safe, though
    if (((second_char_code & 0xdf) == 'A') && ((((unsigned char)sptr[2]) & 0xdf) == 'R')) {
      const uint32_t par_idx_m1 = ((unsigned char)sptr[3]) - '1';
      if ((par_idx_m1 < 2) && (((unsigned char)sptr[4]) <= ' ')) {
	return kChrRawPAR1 + par_idx_m1;
      }
    }
    return -1;
  }
  if (first_char_code == 'C') {
    if (((second_char_code & 0xdf) != 'H') || ((((unsigned char)sptr[2]) & 0xdf) != 'R')) {
      return -1;
    }
    sptr = &(sptr[3]);
    first_char_code = (unsigned char)sptr[0];
    if (first_char_code < 58) {
      goto get_chr_code_raw_digits;
    }
    first_char_code &= 0xdf;
    second_char_code = (unsigned char)sptr[1];
  }
  if (second_char_code <= ' ') {
    return single_cap_letter_chrom(first_char_code);
  }
  if (((unsigned char)sptr[2]) <= ' ') {
    second_char_code &= 0xdf;
    if ((first_char_code == 'X') && (second_char_code == 'Y')) {
      return kChrRawXY;
    } else if ((first_char_code == 'M') && (second_char_code == 'T')) {
      return kChrRawMT;
    }
  }
  return -1;
}

int32_t get_chr_code(const char* chr_name, const chr_info_t* cip, uint32_t name_slen) {
  // requires chr_name to be null-terminated
  // in practice, name_slen will usually already be known, may as well avoid
  // redundant strlen() calls even though this uglifies the interface
  // does not perform exhaustive error-checking
  // -1 = --allow-extra-chr ok, -2 = total fail
  int32_t chr_code_raw = get_chr_code_raw(chr_name);
  if (((uint32_t)chr_code_raw) <= cip->max_code) {
    return chr_code_raw;
  }
  if (chr_code_raw != -1) {
    if (chr_code_raw >= ((int32_t)kMaxContigs)) {
      return cip->xymt_codes[chr_code_raw - kMaxContigs];
    }
    return -2;
  }
  if (!cip->name_ct) {
    return -1;
  }
  // 0xffffffffU gets casted to -1
  return (int32_t)id_htable_find(chr_name, cip->nonstd_names, cip->nonstd_id_htable, name_slen, kChrHtableSize);
}

int32_t get_chr_code_counted(const chr_info_t* cip, uint32_t name_slen, char* chr_name) {
  // when the chromosome name isn't null-terminated
  char* s_end = &(chr_name[name_slen]);
  const char tmpc = *s_end;
  *s_end = '\0';
  const int32_t chr_code = get_chr_code(chr_name, cip, name_slen);
  *s_end = tmpc;
  return chr_code;
}

void chr_error(const char* chr_name, const char* file_descrip, const chr_info_t* cip, uintptr_t line_idx, int32_t error_code) {
  // assumes chr_name is null-terminated
  const int32_t raw_code = get_chr_code_raw(chr_name);
  logprint("\n");
  if (line_idx) {
    LOGERRPRINTFWW("Error: Invalid chromosome code '%s' on line %" PRIuPTR " of %s.\n", chr_name, line_idx, file_descrip);
  } else {
    LOGERRPRINTFWW("Error: Invalid chromosome code '%s' in %s.\n", chr_name, file_descrip);
  }
  if ((raw_code > ((int32_t)cip->max_code)) && ((raw_code <= (int32_t)(kMaxChrTextnum + kChrOffsetCt)) || (raw_code >= ((int32_t)kMaxContigs)))) {
    if (cip->chrset_source == kChrsetSourceDefault) {
      logerrprint("(This is disallowed for humans.  Check if the problem is with your data, or if\nyou forgot to define a different chromosome set with e.g. --chr-set.).\n");
    } else if (cip->chrset_source == kChrsetSourceCmdline) {
      logerrprint("(This is disallowed by your command-line flags.)\n");
    } else {
      // kChrsetSourceFile
      logerrprint("(This is disallowed by the file's own ##chrSet header line.)\n");
    }
    // maybe want to print message(s) depending on whether chromosome set was
    // defined on the command line or by the input file?
  } else if (error_code == -1) {
    logerrprint("(Use --allow-extra-chr to force it to be accepted.)\n");
  }
}

pglerr_t try_to_add_chr_name(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chrs, int32_t* chr_idx_ptr, chr_info_t* cip) {
  // assumes chr_name is either nonstandard (i.e. not "2", "chr2", "chrX",
  // etc.), or a rejected xymt.
  // requires chr_name to be null-terminated
  // assumes chr_idx currently has the return value of get_chr_code()
  if ((!allow_extra_chrs) || ((*chr_idx_ptr) == -2)) {
    chr_error(chr_name, file_descrip, cip, line_idx, *chr_idx_ptr);
    return kPglRetMalformedInput;
  }

  // quasi-bugfix: remove redundant hash table check
  
  if (chr_name[0] == '#') {
    // redundant with some of the comment-skipping loaders, but this isn't
    // performance-critical
    logprint("\n");
    logerrprint("Error: Chromosome/contig names may not begin with '#'.\n");
    return kPglRetMalformedInput;
  }
  if (name_slen > kMaxIdSlen) {
    logprint("\n");
    if (line_idx) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has an excessively long chromosome/contig name. (The " PROG_NAME_STR " limit is " MAX_ID_SLEN_STR " characters.)\n", line_idx, file_descrip);
    } else {
      LOGERRPRINTFWW("Error: Excessively long chromosome/contig name in %s. (The " PROG_NAME_STR " limit is " MAX_ID_SLEN_STR " characters.)\n", file_descrip);
    }
    return kPglRetMalformedInput;
  }
  const uint32_t max_code_p1 = cip->max_code + 1;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = max_code_p1 + name_ct;
  if (chr_code_end == kMaxContigs) {
    logprint("\n");
    logerrprint("Error: Too many distinct nonstandard chromosome/contig names.\n");
    return kPglRetMalformedInput;
  }
  if (!name_ct) {
    // lazy initialization
    fill_uint_one(kChrHtableSize, cip->nonstd_id_htable);
  }
  char** nonstd_names = cip->nonstd_names;
  if (pgl_malloc(name_slen + 1, &(nonstd_names[chr_code_end]))) {
    return kPglRetNomem;
  }
  ll_str_t* name_stack_ptr = cip->incl_excl_name_stack;
  uint32_t in_name_stack = 0;
  while (name_stack_ptr) {
    // there shouldn't be many of these, so sorting is unimportant
    if (!strcmp(chr_name, name_stack_ptr->ss)) {
      in_name_stack = 1;
      break;
    }
    name_stack_ptr = name_stack_ptr->next;
  }
  if ((in_name_stack && cip->is_include_stack) || ((!in_name_stack) && (!cip->is_include_stack))) {
    SET_BIT(chr_code_end, cip->chr_mask);
    if (cip->haploid_mask[0] & 1) {
      SET_BIT(chr_code_end, cip->haploid_mask);
    }
  }
  memcpy(nonstd_names[chr_code_end], chr_name, name_slen + 1);
  *chr_idx_ptr = (int32_t)chr_code_end;
  cip->name_ct = name_ct + 1;
  uint32_t* id_htable = cip->nonstd_id_htable;
  uint32_t hashval = hashceil(chr_name, name_slen, kChrHtableSize);
  while (1) {
    if (id_htable[hashval] == 0xffffffffU) {
      id_htable[hashval] = chr_code_end;
      return kPglRetSuccess;
    }
    if (++hashval == kChrHtableSize) {
      hashval = 0;
    }
  }
}

uintptr_t popcount_bit_idx(const uintptr_t* bitvec, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / kBitsPerWord;
  const uintptr_t start_idxlr = start_idx & (kBitsPerWord - 1);
  const uintptr_t end_idxl = end_idx / kBitsPerWord;
  const uintptr_t end_idxlr = end_idx & (kBitsPerWord - 1);
  uintptr_t ct = 0;
  if (start_idxl == end_idxl) {
    return popcount_long(bitvec[start_idxl] & ((k1LU << end_idxlr) - (k1LU << start_idxlr)));
  }
  if (start_idxlr) {
    ct = popcount_long(bitvec[start_idxl++] >> start_idxlr);
  }
  if (end_idxl > start_idxl) {
    ct += popcount_longs_nzbase(bitvec, start_idxl, end_idxl);
  }
  if (end_idxlr) {
    ct += popcount_long(bitvec[end_idxl] & ((k1LU << end_idxlr) - k1LU));
  }
  return ct;
}

static inline uintptr_t popcount_vecs_intersect(const vul_t* __restrict vvec1_iter, const vul_t* __restrict vvec2_iter, uintptr_t vec_ct) {
  // popcounts vvec1 AND vvec2[0..(ct-1)].  ct is a multiple of 3.
  assert(!(vec_ct % 3));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);
  const vul_t m8 = VCONST_UL(kMask00FF);
  uintptr_t tot = 0;
  while (1) {
    univec_t acc;
    acc.vi = vul_setzero();
    const vul_t* vvec1_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
	return tot;
      }
      vvec1_stop = &(vvec1_iter[vec_ct]);
      vec_ct = 0;
    } else {
      vvec1_stop = &(vvec1_iter[30]);
      vec_ct -= 30;
    }
    do {
      vul_t count1 = (*vvec1_iter++) & (*vvec2_iter++);
      vul_t count2 = (*vvec1_iter++) & (*vvec2_iter++);
      vul_t half1 = (*vvec1_iter++) & (*vvec2_iter++);
      const vul_t half2 = vul_rshift(half1, 1) & m1;
      half1 = half1 & m1;
      count1 = count1 - (vul_rshift(count1, 1) & m1);
      count2 = count2 - (vul_rshift(count2, 1) & m1);
      count1 = count1 + half1;
      count2 = count2 + half2;
      count1 = (count1 & m2) + (vul_rshift(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vul_rshift(count2, 2) & m2);
      acc.vi = acc.vi + (count1 & m4) + (vul_rshift(count1, 4) & m4);
    } while (vvec1_iter < vvec1_stop);
    acc.vi = (acc.vi & m8) + (vul_rshift(acc.vi, 8) & m8);
    tot += univec_hsum_16bit(acc);
  }
}

uintptr_t popcount_longs_intersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct) {
  uintptr_t tot = 0;
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  const uintptr_t trivec_ct = word_ct / (3 * kWordsPerVec);
  tot += popcount_vecs_intersect((const vul_t*)bitvec1_iter, (const vul_t*)bitvec2_iter, trivec_ct * 3);
  bitvec1_iter = &(bitvec1_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2_iter = &(bitvec2_iter[trivec_ct * (3 * kWordsPerVec)]);
  while (bitvec1_iter < bitvec1_end) {
    tot += popcount_long((*bitvec1_iter++) & (*bitvec2_iter++));
  }
  return tot;
}

/*
uintptr_t count_11_vecs(const vul_t* geno_vvec, uintptr_t vec_ct) {
  // Counts number of aligned 11s in vptr[0..(vec_ct-1)].  Assumes vec_ct is a
  // multiple of 6 (0 ok).
  assert(!(vec_ct % 6));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);
  const vul_t m8 = VCONST_UL(kMask00FF);
  const vul_t* geno_vvec_iter = geno_vvec;
  const vul_t* geno_vvec_end = &(geno_vvec[vec_ct]);
  uintptr_t tot = 0;

  while (1) {
    const vul_t* geno_vvec_stop = &(geno_vvec_iter[60]);

    univec_t acc;
    acc.vi = vul_setzero();
    
    if (geno_vvec_stop > geno_vvec_end) {
      if (geno_vvec_iter == geno_vvec_end) {
	return tot;
      }
      geno_vvec_stop = geno_vvec_end;
    }
    do {
      vul_t cur_geno_vword = *geno_vvec_iter++;
      vul_t count1 = cur_geno_vword & m1;
      count1 = count1 & vul_rshift(cur_geno_vword, 1);
      
      cur_geno_vword = *geno_vvec_iter++;
      vul_t cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count1 = count1 + cur_11;

      cur_geno_vword = *geno_vvec_iter++;
      cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count1 = count1 + cur_11;
      count1 = (count1 & m2) + (vul_rshift(count1, 2) & m2);

      cur_geno_vword = *geno_vvec_iter++;
      vul_t count2 = cur_geno_vword & m1;
      count2 = count2 & vul_rshift(cur_geno_vword, 1);
      
      cur_geno_vword = *geno_vvec_iter++;
      vul_t cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count2 = count2 + cur_11;

      cur_geno_vword = *geno_vvec_iter++;
      cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count2 = count2 + cur_11;
      count1 = count1 + (count2 & m2) + (vul_rshift(count2, 2) & m2);

      acc.vi = acc.vi + (count1 & m4) + (vul_rshift(count1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    acc.vi = (acc.vi & m8) + (vul_rshift(acc.vi, 8) & m8);
    tot += univec_hsum_16bit(acc);
  }
}

uintptr_t count_11_longs(const uintptr_t* genovec, uintptr_t word_ct) {
  uintptr_t tot = 0;
  if (word_ct >= (6 * kWordsPerVec)) {
    assert(IS_VEC_ALIGNED(genovec));
    const uintptr_t remainder = word_ct % (6 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = count_11_vecs((const vul_t*)genovec, main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    genovec = &(genovec[main_block_word_ct]);
  }
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx < word_ct; ++trailing_word_idx) {
    const uintptr_t cur_geno_word = genovec[trailing_word_idx];
    tot += popcount01_long(cur_geno_word & (cur_geno_word >> 1) & kMask5555);
  }
}
*/

uint32_t are_all_bits_zero(const uintptr_t* bitarr, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / kBitsPerWord;
  const uintptr_t start_idxlr = start_idx & (kBitsPerWord - 1);
  const uintptr_t end_idxl = end_idx / kBitsPerWord;
  const uintptr_t end_idxlr = end_idx & (kBitsPerWord - 1);
  if (start_idxl == end_idxl) {
    return !(bitarr[start_idxl] & ((k1LU << end_idxlr) - (k1LU << start_idxlr)));
  }
  if (start_idxlr && (bitarr[start_idxl++] >> start_idxlr)) {
    return 0;
  }
  for (; start_idxl < end_idxl; ++start_idxl) {
    if (bitarr[start_idxl]) {
      return 0;
    }
  }
  if (!end_idxlr) {
    return 1;
  }
  return !(bitarr[end_idxl] & ((k1LU << end_idxlr) - k1LU));
}

void copy_bitarr_range(const uintptr_t* __restrict src_bitarr, uintptr_t src_start_bitidx, uintptr_t target_start_bitidx, uintptr_t len, uintptr_t* __restrict target_bitarr) {
  // assumes len is positive, and relevant bits of target_bitarr are zero
  const uintptr_t* src_bitarr_iter = &(src_bitarr[src_start_bitidx / kBitsPerWord]);
  uint32_t src_rshift = src_start_bitidx % kBitsPerWord;
  uintptr_t* target_bitarr_iter = &(target_bitarr[target_start_bitidx / kBitsPerWord]);
  uint32_t target_initial_lshift = target_start_bitidx % kBitsPerWord;
  uintptr_t cur_src_word;
  if (target_initial_lshift) {
    const uint32_t initial_copy_bitct = kBitsPerWord - target_initial_lshift;
    if (len <= initial_copy_bitct) {
      goto copy_bitarr_range_last_partial_word;
    }
    cur_src_word = (*src_bitarr_iter) >> src_rshift;
    if (src_rshift >= target_initial_lshift) {
      ++src_bitarr_iter;
      cur_src_word |= (*src_bitarr_iter) << (kBitsPerWord - src_rshift);
    }
    *target_bitarr_iter++ |= cur_src_word << target_initial_lshift;
    src_rshift = (src_rshift + initial_copy_bitct) % kBitsPerWord;
    len -= initial_copy_bitct;
  }
  {
    const uintptr_t fullword_ct = len / kBitsPerWord;
    if (!src_rshift) {
      memcpy(target_bitarr_iter, src_bitarr_iter, fullword_ct * sizeof(intptr_t));
      target_bitarr_iter = &(target_bitarr_iter[fullword_ct]);
      src_bitarr_iter = &(src_bitarr_iter[fullword_ct]);
    } else {
      const uint32_t src_lshift = kBitsPerWord - src_rshift;
      cur_src_word = *src_bitarr_iter;
      for (uintptr_t widx = 0; widx < fullword_ct; ++widx) {
	const uintptr_t next_src_word = *(++src_bitarr_iter);
	*target_bitarr_iter++ = (cur_src_word >> src_rshift) | (next_src_word << src_lshift);
	cur_src_word = next_src_word;
      }
    }
  }
  len %= kBitsPerWord;
  if (len) {
    target_initial_lshift = 0;
  copy_bitarr_range_last_partial_word:
    cur_src_word = (*src_bitarr_iter) >> src_rshift;
    if (len + src_rshift > kBitsPerWord) {
      cur_src_word |= src_bitarr_iter[1] << (kBitsPerWord - src_rshift);
    }
    *target_bitarr_iter |= (cur_src_word & ((~k0LU) >> (kBitsPerWord - ((uint32_t)len)))) << target_initial_lshift;
  }
}

void interleaved_mask_zero(const uintptr_t* __restrict interleaved_mask, uintptr_t vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = vec_ct / 2;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t* interleaved_mask_iter = (const vul_t*)interleaved_mask;
  vul_t* genovvec_iter = (vul_t*)genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const vul_t mask_vvec = *interleaved_mask_iter++;
    vul_t mask_first = mask_vvec & m1;
    mask_first = mask_first | vul_lshift(mask_first, 1);
    vul_t mask_second = (~m1) & mask_vvec;
    mask_second = mask_second | vul_rshift(mask_second, 1);
    *genovvec_iter = (*genovvec_iter) & mask_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) & mask_second;
    ++genovvec_iter;
  }
  if (vec_ct & 1) {
    vul_t mask_first = *interleaved_mask_iter;
    mask_first = mask_first | vul_lshift(mask_first, 1);
    *genovvec_iter = (*genovvec_iter) & mask_first;
  }
#else
  const uintptr_t* interleaved_mask_iter = interleaved_mask;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const uintptr_t mask_word = *interleaved_mask_iter++;
    *genovec_iter &= (mask_word & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter &= ((mask_word >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (vec_ct & 1) {
    const uintptr_t mask_word = *interleaved_mask_iter;
    *genovec_iter &= mask_word * 3;
  }
#endif
}

void interleaved_set_missing(const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = vec_ct / 2;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t* interleaved_set_iter = (const vul_t*)interleaved_set;
  vul_t* genovvec_iter = (vul_t*)genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const vul_t set_vvec = *interleaved_set_iter++;
    vul_t set_first = set_vvec & m1;
    set_first = set_first | vul_lshift(set_first, 1);
    vul_t set_second = (~m1) & set_vvec;
    set_second = set_second | vul_rshift(set_second, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) | set_second;
    ++genovvec_iter;
  }
  if (vec_ct & 1) {
    vul_t set_first = *interleaved_set_iter;
    set_first = set_first | vul_lshift(set_first, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
  }
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    *genovec_iter |= (set_word & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter |= ((set_word >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (vec_ct & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    *genovec_iter |= set_word * 3;
  }
#endif
}

void set_male_het_missing(const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec) {
  const uint32_t twovec_ct = vec_ct / 2;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t* sex_male_interleaved_iter = (const vul_t*)sex_male_interleaved;
  vul_t* genovvec_iter = (vul_t*)genovec;
  for (uint32_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const vul_t sex_male_vvec = *sex_male_interleaved_iter++;
    // we wish to bitwise-or with (sex_male_quatervec_01 & genovec) << 1
    const vul_t sex_male_first = sex_male_vvec & m1;
    const vul_t sex_male_second_shifted = (~m1) & sex_male_vvec;
    vul_t cur_geno_vword = *genovvec_iter;
    
    const vul_t missing_male_vword = sex_male_first & cur_geno_vword;
    
    *genovvec_iter++ = cur_geno_vword | vul_lshift(missing_male_vword, 1);
    cur_geno_vword = *genovvec_iter;
    *genovvec_iter++ = cur_geno_vword | (sex_male_second_shifted & vul_lshift(cur_geno_vword, 1));
  }
  if (vec_ct & 1) {
    const vul_t sex_male_first = (*sex_male_interleaved_iter) & m1;
    const vul_t cur_geno_vword = *genovvec_iter;
    const vul_t missing_male_vword = sex_male_first & cur_geno_vword;
    *genovvec_iter = cur_geno_vword | vul_lshift(missing_male_vword, 1);
  }
#else
  const uintptr_t* sex_male_interleaved_iter = sex_male_interleaved;
  uintptr_t* genovec_iter = genovec;
  for (uint32_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const uintptr_t sex_male_word = *sex_male_interleaved_iter++;
    uintptr_t cur_geno_word = *genovec_iter;
    *genovec_iter++ = cur_geno_word | ((sex_male_word & kMask5555 & cur_geno_word) << 1);
    cur_geno_word = *genovec_iter;
    *genovec_iter++ = cur_geno_word | (sex_male_word & kMaskAAAA & (cur_geno_word << 1));
  }
  if (vec_ct & 1) {
    const uintptr_t sex_male_word = *sex_male_interleaved_iter;
    uintptr_t cur_geno_word = *genovec_iter;
    *genovec_iter = cur_geno_word | ((sex_male_word & kMask5555 & cur_geno_word) << 1);
  }
#endif
}

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
//
// Similar to pgr_detect_genovec_hets_unsafe(). 
void mask_genovec_hets_unsafe(const uintptr_t* __restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr) {
  halfword_t* bitarr_alias = (halfword_t*)bitarr;
  for (uint32_t widx = 0; widx < raw_sample_ctl2; ++widx) {
    const uintptr_t cur_word = genovec[widx];
    uintptr_t ww = (~(cur_word >> 1)) & cur_word & kMask5555; // low 1, high 0
    bitarr_alias[widx] &= pack_word_to_halfword(ww);
  }
}

/*
uint32_t chr_window_max(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bp, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max) {
  if (cur_window_max >= ct_max) {
    return ct_max;
  }
  const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  uint32_t variant_uidx = next_set(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
  const uint32_t variant_ct = popcount_bit_idx(variant_include, variant_uidx, chr_end);
  if (variant_ct <= cur_window_max) {
    return cur_window_max;
  }
  uint32_t window_idx_first = 0;
  uint32_t window_uidx_first = variant_uidx;
  uint32_t window_bp_first = variant_bp[variant_uidx];
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_uidx, ++variant_idx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    uint32_t variant_bp_thresh = variant_bp[variant_uidx];
    if (variant_bp_thresh < bp_max) {
      variant_bp_thresh = 0;
    } else {
      variant_bp_thresh -= bp_max;
    }
    if (variant_bp_thresh > window_bp_first) {
      do {
        ++window_uidx_first;
        next_set_unsafe_ck(variant_include, &window_uidx_first);
        window_bp_first = variant_bp[window_uidx_first];
        ++window_idx_first;
      } while (variant_bp_thresh > window_bp_first);
    } else if (variant_idx - window_idx_first == cur_window_max) {
      if (++cur_window_max == ct_max) {
	return cur_window_max;
      }
    }
  }
  return cur_window_max;
}
*/

// advances forward_ct set bits; forward_ct must be positive.  (stays put if
// forward_ct == 1 and current bit is set.  may want to tweak this interface,
// easy to introduce off-by-one bugs...)
// In usual 64-bit case, also assumes bitvec is 16-byte aligned and the end of
// the trailing 16-byte block can be safely read from.
uintptr_t jump_forward_set_unsafe(const uintptr_t* bitvec, uintptr_t cur_pos, uintptr_t forward_ct) {
  assert(forward_ct);
  uintptr_t widx = cur_pos / kBitsPerWord;
  uintptr_t ulii = cur_pos % kBitsPerWord;
  const uintptr_t* bptr = &(bitvec[widx]);
  uintptr_t uljj;
  uintptr_t ulkk;
#ifdef __LP64__
  const vul_t* vptr;
  assert(IS_VEC_ALIGNED(bitvec));
#endif
  if (ulii) {
    uljj = (*bptr) >> ulii;
    ulkk = popcount_long(uljj);
    if (ulkk >= forward_ct) {
    jump_forward_set_unsafe_finish:
      ulkk = CTZLU(uljj);
      while (--forward_ct) {
        uljj &= uljj - 1;
        ulkk = CTZLU(uljj);
      }
      return widx * kBitsPerWord + ulii + ulkk;
    }
    forward_ct -= ulkk;
    ++widx;
    ++bptr;
  }
  ulii = 0;
#ifdef __LP64__
  while (widx & (kWordsPerVec - k1LU)) {
    uljj = *bptr;
    ulkk = popcount_long(uljj);
    if (ulkk >= forward_ct) {
      goto jump_forward_set_unsafe_finish;
    }
    forward_ct -= ulkk;
    ++widx;
    ++bptr;
  }
  vptr = (const vul_t*)bptr;
  while (forward_ct > kBitsPerWord * (3 * kWordsPerVec)) {
    uljj = ((forward_ct - 1) / (kBitsPerWord * (3 * kWordsPerVec))) * 3;
    ulkk = popcount_vecs(vptr, uljj);
    vptr = &(vptr[uljj]);
    forward_ct -= ulkk;
  }
  bptr = (const uintptr_t*)vptr;
  while (forward_ct > kBitsPerWord) {
    forward_ct -= popcount_long(*bptr++);
  }
#else
  while (forward_ct > kBitsPerWord) {
    uljj = (forward_ct - 1) / kBitsPerWord;
    ulkk = popcount_longs(bptr, uljj);
    bptr = &(bptr[uljj]);
    forward_ct -= ulkk;
  }
#endif
  while (1) {
    uljj = *bptr;
    ulkk = popcount_long(uljj);
    if (ulkk >= forward_ct) {
      widx = (uintptr_t)(bptr - bitvec);
      goto jump_forward_set_unsafe_finish;
    }
    forward_ct -= ulkk;
    ++bptr;
  }
}

void compute_uidx_start_partition(const uintptr_t* variant_include, uint64_t variant_ct, uint32_t thread_ct, uint32_t first_variant_uidx, uint32_t* variant_uidx_starts) {
  assert(variant_ct);
  uint32_t cur_variant_uidx_start = next_set_unsafe(variant_include, first_variant_uidx);
  uint32_t cur_variant_idx_start = 0;
  variant_uidx_starts[0] = cur_variant_uidx_start;
  for (uint32_t tidx = 1; tidx < thread_ct; ++tidx) {
    const uint32_t new_variant_idx_start = (tidx * variant_ct) / thread_ct;
    if (new_variant_idx_start != cur_variant_idx_start) {
      cur_variant_uidx_start = jump_forward_set_unsafe(variant_include, cur_variant_uidx_start + 1, new_variant_idx_start - cur_variant_idx_start);
      cur_variant_idx_start = new_variant_idx_start;
    }
    variant_uidx_starts[tidx] = cur_variant_uidx_start;
  }
}

uint32_t not_only_xymt(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t raw_variant_ct, uint32_t xymt_offset) {
  const uint32_t xymt_code = (uint32_t)cip->xymt_codes[xymt_offset];
  const uint32_t cur_chr_fo_idx = cip->chr_idx_to_foidx[xymt_code];
  const uint32_t chr_start = cip->chr_fo_vidx_start[cur_chr_fo_idx];
  if (chr_start) {
    const uint32_t first_uidx = next_set_unsafe(variant_include, 0);
    if (first_uidx < chr_start) {
      return 1;
    }
  }
  const uint32_t chr_end = cip->chr_fo_vidx_start[cur_chr_fo_idx + 1];
  return (chr_end < raw_variant_ct) && (next_set(variant_include, chr_end, raw_variant_ct) != raw_variant_ct);
}

uint32_t count_non_autosomal_variants(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t count_x, uint32_t count_mt) {
  // for backward compatibility, unplaced markers are considered to be
  // autosomal here
  uint32_t ct = 0;
  if (count_x) {
    int32_t x_code;
    if (xymt_exists(cip, kChrOffsetX, &x_code)) {
      ct += count_chr_variants_unsafe(variant_include, cip, x_code);
    }
  }
  int32_t y_code;
  if (xymt_exists(cip, kChrOffsetY, &y_code)) {
    ct += count_chr_variants_unsafe(variant_include, cip, y_code);
  }
  if (count_mt) {
    int32_t mt_code;
    if (xymt_exists(cip, kChrOffsetMT, &mt_code)) {
      ct += count_chr_variants_unsafe(variant_include, cip, mt_code);
    }
  }
  return ct;
}

pglerr_t conditional_allocate_non_autosomal_variants(const chr_info_t* cip, const char* calc_descrip, uint32_t raw_variant_ct, uintptr_t** variant_include_ptr, uint32_t* variant_ct_ptr) {
  const uint32_t non_autosomal_variant_ct = count_non_autosomal_variants(*variant_include_ptr, cip, 1, 1);
  if (!non_autosomal_variant_ct) {
    return kPglRetSuccess;
  }
  LOGPRINTF("Excluding %u variant%s on non-autosomes from %s.\n", non_autosomal_variant_ct, (non_autosomal_variant_ct == 1)? "" : "s", calc_descrip);
  *variant_ct_ptr -= non_autosomal_variant_ct;
  if (!(*variant_ct_ptr)) {
    // this may not always be an error condition
    LOGERRPRINTF("Error: No variants remaining for %s.\n", calc_descrip);
    return kPglRetInconsistentInput;
  }
  const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
  uintptr_t* working_variant_include;
  if (bigstack_alloc_ul(raw_variant_ctl, &working_variant_include)) {
    return kPglRetNomem;
  }
  memcpy(working_variant_include, *variant_include_ptr, raw_variant_ctl * sizeof(intptr_t));
  int32_t x_code;
  if (xymt_exists(cip, kChrOffsetX, &x_code)) {
    uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
    clear_bits_nz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], working_variant_include);
  }
  int32_t y_code;
  if (xymt_exists(cip, kChrOffsetX, &y_code)) {
    uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)y_code];
    clear_bits_nz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], working_variant_include);
  }
  int32_t mt_code;
  if (xymt_exists(cip, kChrOffsetX, &mt_code)) {
    uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)mt_code];
    clear_bits_nz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], working_variant_include);
  }
  *variant_include_ptr = working_variant_include;
  return kPglRetSuccess;
}

void fill_subset_chr_fo_vidx_start(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t* subset_chr_fo_vidx_start) {
  const uint32_t chr_ct = cip->chr_ct;
  subset_chr_fo_vidx_start[0] = 0;
  uint32_t variant_uidx = 0;
  uint32_t variant_idx = 0;
  for (uint32_t chr_fo_idx = 1; chr_fo_idx <= chr_ct; ++chr_fo_idx) {
    const uint32_t chr_end_variant_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
    variant_idx += popcount_bit_idx(variant_include, variant_uidx, chr_end_variant_uidx);
    subset_chr_fo_vidx_start[chr_fo_idx] = variant_idx;
    variant_uidx = chr_end_variant_uidx;
  }
}

boolerr_t allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr) {
  char* newptr;
  if (allele_slen == 1) {
    // const_cast
    newptr = (char*)((uintptr_t)(&(g_one_char_strs[((unsigned char)(*newval)) * 2])));
  } else {
    char* new_alloc;
    if (pgl_malloc(allele_slen + 1, &new_alloc)) {
      return 1;
    }
    memcpyx(new_alloc, newval, allele_slen, '\0');
    newptr = new_alloc;
  }
  *allele_ptr = newptr;
  return 0;
}

boolerr_t allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr) {
  char* newptr;
  if (allele_slen == 1) {
    // const_cast
    newptr = (char*)((uintptr_t)(&(g_one_char_strs[((unsigned char)(*newval)) * 2])));
  } else {
    char* new_alloc;
    if (pgl_malloc(allele_slen + 1, &new_alloc)) {
      return 1;
    }
    memcpyx(new_alloc, newval, allele_slen, '\0');
    newptr = new_alloc;
  }
  const uintptr_t bigstack_end_addr = (uintptr_t)g_bigstack_end;
  const uintptr_t maxdiff = ((uintptr_t)(&(g_one_char_strs[512]))) - bigstack_end_addr;
  // take advantage of unsigned wraparound
  if ((((uintptr_t)(*allele_ptr)) - bigstack_end_addr) >= maxdiff) {
    free(*allele_ptr);
  }
  *allele_ptr = newptr;
  return 0;
}

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, char** allele_storage) {
  // Now doesn't improperly free bigstack allocations (as long as they aren't
  // past g_bigstack_end), and doesn't need to be called at all most of the
  // time.
  
  // An alternative representation: have a separate bitarray which indicates
  // whether the allele_storage[] element should be interpreted as a heap
  // pointer or an in-place zero-terminated string (i.e. string length can be
  // up to 7 on 64-bit systems).  I expect that to be more efficient for new
  // datasets, but let's get the simple (and 1.9-codebase-compatible)
  // implementation working first, and then benchmark the fancier code later.
  if (allele_storage && (max_allele_slen > 1)) {
    const uintptr_t bigstack_end_addr = (uintptr_t)g_bigstack_end;
    const uintptr_t maxdiff = ((uintptr_t)(&(g_one_char_strs[512]))) - bigstack_end_addr;
    for (uintptr_t idx = 0; idx < allele_storage_entry_ct; ++idx) {
      char* cur_entry = allele_storage[idx];
      assert(cur_entry);
      // take advantage of unsigned wraparound
      if ((((uintptr_t)cur_entry) - bigstack_end_addr) >= maxdiff) {
	free(cur_entry);
      }
    }
  }
}

char g_missing_catname[kMaxMissingPhenostrBlen];
char g_output_missing_pheno[kMaxMissingPhenostrBlen];
char g_legacy_output_missing_pheno[kMaxMissingPhenostrBlen];

void init_pheno() {
  strcpy(g_missing_catname, "NONE");
  strcpy(g_output_missing_pheno, "NA");
  strcpy(g_legacy_output_missing_pheno, "-9");
}

uint32_t is_categorical_phenostr(const char* phenostr) {
  uint32_t first_char_code = (unsigned char)(*phenostr++);
  // allow leading +/-
  if ((first_char_code == 43) || (first_char_code == 45)) {
    first_char_code = (unsigned char)(*phenostr++);
  }
  if (((first_char_code - 48) < 10) || (first_char_code == 44) || (first_char_code < 32)) {
    // the last two conditions are for detecting CSV empty strings
    return 0;
  }
  if (first_char_code == 46) {
    // decimal point.  classify based on whether next character is a digit.
    const uint32_t second_char_code = (unsigned char)phenostr[0];
    return ((second_char_code - 48) >= 10);
  }
  // allow any capitalization of "NA"/"nan", but not "inf"
  if ((first_char_code & 0xdf) != 78) {
    return 1;
  }
  const uint32_t second_char_code = (unsigned char)phenostr[0];
  if ((second_char_code & 0xdf) != 65) {
    return 1;
  }
  const uint32_t third_char_code = (unsigned char)phenostr[1];
  if ((third_char_code & 0xdf) == 78) {
    return (((unsigned char)phenostr[2]) > ' ');
  }
  return (third_char_code > 32);
}

uint32_t is_categorical_phenostr_nocsv(const char* phenostr) {
  uint32_t first_char_code = (unsigned char)(*phenostr++);
  // allow leading +/-
  if ((first_char_code == 43) || (first_char_code == 45)) {
    first_char_code = (unsigned char)(*phenostr++);
  }
  if ((first_char_code - 48) < 10) {
    return 0;
  }
  if (first_char_code == 46) {
    // decimal point.  classify based on whether next character is a digit.
    const uint32_t second_char_code = (unsigned char)phenostr[0];
    return ((second_char_code - 48) >= 10);
  }
  // allow any capitalization of "NA"/"nan", but not "inf"
  if ((first_char_code & 0xdf) != 78) {
    return 1;
  }
  const uint32_t second_char_code = (unsigned char)phenostr[0];
  if ((second_char_code & 0xdf) != 65) {
    return 1;
  }
  const uint32_t third_char_code = (unsigned char)phenostr[1];
  if ((third_char_code & 0xdf) == 78) {
    return (((unsigned char)phenostr[2]) > ' ');
  }
  return (third_char_code > 32);
}

uint32_t first_cc_or_qt_pheno_idx(const pheno_col_t* pheno_cols, uint32_t pheno_ct) {
  for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
    if (pheno_cols[pheno_idx].type_code < kPhenoDtypeCat) {
      return pheno_idx;
    }
  }
  return 0xffffffffU;
}

uint32_t is_const_covar(const pheno_col_t* covar_col, const uintptr_t* sample_include, uint32_t sample_ct) {
  if (sample_ct < 2) {
    return 1;
  }
  uint32_t sample_uidx = next_set_unsafe(sample_include, 0);
  if (covar_col->type_code == kPhenoDtypeQt) {
    const double* covar_vals = covar_col->data.qt;
    const double first_covar_val = covar_vals[sample_uidx];
    for (uint32_t sample_idx = 1; sample_idx < sample_ct; ++sample_idx) {
      ++sample_uidx;
      next_set_unsafe_ck(sample_include, &sample_uidx);
      if (covar_vals[sample_uidx] != first_covar_val) {
	return 0;
      }
    }
    return 1;
  }
  assert(covar_col->type_code == kPhenoDtypeCat);
  const uint32_t* covar_vals = covar_col->data.cat;
  const uint32_t first_covar_val = covar_vals[sample_uidx];
  for (uint32_t sample_idx = 1; sample_idx < sample_ct; ++sample_idx) {
    ++sample_uidx;
    next_set_unsafe_ck(sample_include, &sample_uidx);
    if (covar_vals[sample_uidx] != first_covar_val) {
      return 0;
    }
  }
  return 1;
}

uint32_t identify_remaining_cats(const uintptr_t* sample_include, const pheno_col_t* covar_col, uint32_t sample_ct, uintptr_t* cat_covar_wkspace) {
  // assumes covar_col->type_code == kPhenoTypeCat
  const uint32_t nonnull_cat_ct = covar_col->nonnull_category_ct;
  const uint32_t* covar_vals = covar_col->data.cat;
  const uint32_t word_ct = 1 + (nonnull_cat_ct / kBitsPerWord);
  fill_ulong_zero(word_ct, cat_covar_wkspace);
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include, &sample_uidx);
    set_bit(covar_vals[sample_uidx], cat_covar_wkspace);
  }
  return popcount_longs(cat_covar_wkspace, word_ct);
}

void cleanup_pheno_cols(uint32_t pheno_ct, pheno_col_t* pheno_cols) {
  if (pheno_cols) {
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      vecaligned_free_cond(pheno_cols[pheno_idx].nonmiss);
    }
    free(pheno_cols);
  }
}

boolerr_t parse_next_range(char** argv, uint32_t param_ct, char range_delim, uint32_t* cur_param_idx_ptr, char** cur_arg_pptr, char** range_start_ptr, uint32_t* rs_len_ptr, char** range_end_ptr, uint32_t* re_len_ptr) {
  // Starts reading from argv[cur_param_idx][cur_pos].  If a valid range is
  // next, range_start + rs_len + range_end + re_len are updated.  If only a
  // single item is next, range_end is set to nullptr and range_start + rs_len
  // are updated.  If there are no items left, range_start is set to nullptr.
  // If the input is not well-formed, -1 is returned instead of 0.
  uint32_t cur_param_idx = *cur_param_idx_ptr;
  if (cur_param_idx > param_ct) {
    *cur_arg_pptr = nullptr;
    return 0;
  }
  char* cur_arg_ptr = *cur_arg_pptr;
  while (1) {
    char cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
	*range_start_ptr = nullptr;
	return 0;
      }
      cur_arg_ptr = argv[cur_param_idx];
      cc = *cur_arg_ptr;
    }
    if (cc == range_delim) {
      return 1;
    }
    if (cc != ',') {
      break;
    }
    ++cur_arg_ptr;
  }
  *range_start_ptr = cur_arg_ptr;
  char cc;
  do {
    cc = *(++cur_arg_ptr);
    if ((!cc) || (cc == ',')) {
      *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
      *cur_arg_pptr = cur_arg_ptr;
      *range_end_ptr = nullptr;
      return 0;
    }
  } while (cc != range_delim);
  *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
  cc = *(++cur_arg_ptr);
  if ((!cc) || (cc == ',') || (cc == range_delim)) {
    return 1;
  }
  *range_end_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if (cc == range_delim) {
      return 1;
    }
  } while (cc && (cc != ','));
  *re_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_end_ptr));
  *cur_arg_pptr = cur_arg_ptr;
  return 0;
}

pglerr_t parse_chr_ranges(const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t allow_extra_chrs, uint32_t xymt_subtract, char range_delim, char** argv, chr_info_t* cip, uintptr_t* chr_mask) {
  pglerr_t reterr = kPglRetSuccess;
  {
    char* cur_arg_ptr = argv[1];
    char* range_end = nullptr;
    uint32_t cur_param_idx = 1;
    uint32_t rs_len = 0;
    uint32_t re_len = 0;
    while (1) {
      char* range_start;
      if (parse_next_range(argv, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", flagname_p, argv[cur_param_idx]);
	goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
      }
      if (!range_start) {
	break;
      }
      const char cc = range_start[rs_len];
      range_start[rs_len] = '\0';
      int32_t chr_code_start = get_chr_code_raw(range_start);
      if (chr_code_start < 0) {
	if (!allow_extra_chrs) {
	  sprintf(g_logbuf, "Error: Invalid --%s chromosome code '%s'.\n", flagname_p, range_start);
	  goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
	}
	if (range_end) {
	  goto parse_chr_ranges_ret_INVALID_CMDLINE_NONSTD;
	}
        if (push_llstr(range_start, &(cip->incl_excl_name_stack))) {
	  goto parse_chr_ranges_ret_NOMEM;
	}
      } else {
	if (chr_code_start >= ((int32_t)kMaxContigs)) {
	  chr_code_start -= xymt_subtract;
	}
	if (range_end) {
	  const char cc2 = range_end[re_len];
	  range_end[re_len] = '\0';
	  int32_t chr_code_end = get_chr_code_raw(range_end);
	  if (chr_code_end < 0) {
	    if (!allow_extra_chrs) {
	      sprintf(g_logbuf, "Error: Invalid --%s chromosome code '%s'.\n", flagname_p, range_end);
	      goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
	    }
	    goto parse_chr_ranges_ret_INVALID_CMDLINE_NONSTD;
	  }
	  if (chr_code_end >= ((int32_t)kMaxContigs)) {
	    // prohibit stuff like "--chr par1-par2", "--chr x-y", "--chr x-26"
	    sprintf(g_logbuf, "Error: --%s chromosome code '%s' cannot be the end of a range.\n", flagname_p, range_end);
	    goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
	  }
	  if (chr_code_end <= chr_code_start) {
	    sprintf(g_logbuf, "Error: --%s chromosome code '%s' is not greater than '%s'.\n", flagname_p, range_end, range_start);
	    goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
	  }
	  range_end[re_len] = cc2;
	  fill_bits_nz(chr_code_start, chr_code_end + 1, chr_mask);
	} else {
          set_bit(chr_code_start, chr_mask);
	}
      }
      range_start[rs_len] = cc;
    }
    // no compelling reason to prohibit "--not-chr ,"
  }
  while (0) {
  parse_chr_ranges_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  parse_chr_ranges_ret_INVALID_CMDLINE_NONSTD:
    logerrprint("Error: Chromosome ranges cannot include nonstandard names.\n");
    reterr = kPglRetInvalidCmdline;
    break;
  parse_chr_ranges_ret_INVALID_CMDLINE_WWA:
    wordwrapb(0);
    logerrprintb();
    logerrprint(errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return reterr;
}

pglerr_t parse_name_ranges(char** argv, const char* errstr_append, uint32_t param_ct, uint32_t require_posint, char range_delim, range_list_t* range_list_ptr) {
  uint32_t name_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t name_max_blen = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  char* cur_name_str;
  char* dup_check;
  unsigned char* cur_name_starts_range;
  uint32_t last_val;
  uint32_t cur_val;
  // two passes.  first pass: count parameters, determine name_max_blen;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(argv, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argv[0], argv[cur_param_idx]);
	logerrprint(errstr_append);
        return kPglRetInvalidCmdline;
      }
      if (!range_start) {
	break;
      }
      ++name_ct;
      if (rs_len > name_max_blen) {
	name_max_blen = rs_len; // does NOT include trailing null yet
      }
      if (range_end) {
	++name_ct;
	if (re_len > name_max_blen) {
	  name_max_blen = re_len;
	}
      }
    }
  }
  if (!name_ct) {
    LOGERRPRINTF("Error: %s requires at least one value.\n%s", argv[0], errstr_append);
    return kPglRetInvalidCmdline;
  }
  range_list_ptr->name_max_blen = ++name_max_blen;
  range_list_ptr->name_ct = name_ct;
  if (pgl_malloc(name_ct * (((uintptr_t)name_max_blen) + 1), &range_list_ptr->names)) {
    return kPglRetNomem;
  }
  range_list_ptr->starts_range = (unsigned char*)(&(range_list_ptr->names[name_ct * ((uintptr_t)name_max_blen)]));
  cur_name_str = range_list_ptr->names;
  cur_name_starts_range = range_list_ptr->starts_range;
  cur_param_idx = 1;
  cur_arg_ptr = argv[1];
  while (1) {
    // second pass; this can't fail since we already validated
    parse_next_range(argv, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      if (require_posint) {
	last_val = 0;
	for (cur_param_idx = 0; cur_param_idx < name_ct; ++cur_param_idx) {
	  cur_name_str = &(range_list_ptr->names[cur_param_idx * ((uintptr_t)name_max_blen)]);
	  dup_check = cur_name_str; // actually a numeric check
	  do {
	    if (is_not_digit(*dup_check)) {
	      LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	      return kPglRetInvalidCmdline;
	    }
	  } while (*(++dup_check));
	  if (scan_posint_defcap(cur_name_str, &cur_val)) {
	    LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	    return kPglRetInvalidCmdline;
	  }
	  if (range_list_ptr->starts_range[cur_param_idx]) {
	    last_val = cur_val;
	  } else {
	    if (cur_val <= last_val) {
	      LOGERRPRINTFWW("Error: Invalid %s range '%s-%s'.\n", argv[0], &(range_list_ptr->names[(cur_param_idx - 1) * name_max_blen]), cur_name_str);
	      return kPglRetInvalidCmdline;
	    }
	    last_val = 0;
	  }
	}
      }
      return kPglRetSuccess;
    }
    memcpyx(cur_name_str, range_start, rs_len, 0);
    dup_check = range_list_ptr->names;
    while (dup_check < cur_name_str) {
      if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	LOGERRPRINTFWW("Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	return kPglRetInvalidCmdline;
      }
      dup_check = &(dup_check[name_max_blen]);
    }
    cur_name_str = &(cur_name_str[name_max_blen]);
    if (range_end) {
      *cur_name_starts_range++ = 1;
      memcpyx(cur_name_str, range_end, re_len, 0);
      dup_check = range_list_ptr->names;
      while (dup_check < cur_name_str) {
	if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	  LOGERRPRINTFWW("Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	  return kPglRetInvalidCmdline;
	}
        dup_check = &(dup_check[name_max_blen]);
      }
      cur_name_str = &(cur_name_str[name_max_blen]);
      *cur_name_starts_range++ = 0;
    } else {
      *cur_name_starts_range++ = 0;
    }
  }
}

void join_threads(uint32_t ctp1, pthread_t* threads) {
  if (!(--ctp1)) {
    return;
  }
#ifdef _WIN32
  WaitForMultipleObjects(ctp1, threads, 1, INFINITE);
  for (uint32_t uii = 0; uii < ctp1; ++uii) {
    // fix handle leak?
    CloseHandle(threads[uii]);
  }
#else
  for (uint32_t uii = 0; uii < ctp1; ++uii) {
    pthread_join(threads[uii], nullptr);
  }
#endif
}

boolerr_t spawn_threads(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, pthread_t* threads) {
  uintptr_t ulii;
  if (ct == 1) {
    return 0;
  }
  for (ulii = 1; ulii < ct; ++ulii) {
#ifdef _WIN32
    threads[ulii - 1] = (HANDLE)_beginthreadex(nullptr, 4096, start_routine, (void*)ulii, 0, nullptr);
    if (!threads[ulii - 1]) {
      join_threads(ulii, threads);
      return 1;
    }
#else
    if (pthread_create(&(threads[ulii - 1]), nullptr, start_routine, (void*)ulii)) {
      join_threads(ulii, threads);
      return 1;
    }
#endif
  }
  return 0;
}

// Main plink 2.0 threading framework:
// * On all operating systems, g_is_last_thread_block indicates whether all
//   threads should terminate upon completion of the current block.
// * On Linux and OS X, if we aren't dealing with the final block,
//   spawn_threads2z() also reinitializes g_thread_active_ct.
// * On Linux and OS X, spawn_threads2z() checks if g_thread_mutex_initialized
//   is set.  If not, it, it is set, g_thread_sync_mutex,
//   g_thread_cur_block_done_condvar and g_thread_start_next_condvar are
//   initialized, then threads are launched.
//   If it has, pthread_cond_broadcast() acts on g_thread_start_next_condvar.
// * On Windows, spawn_threads2z() checks if g_thread_mutex_initialized is set.
//   If it has not, it, along with g_thread_start_next_event[] and
//   g_thread_cur_block_done_events[], are initialized, then the threads are
//   launched.  If it has, SetEvent() acts on g_thread_start_next_event[].
//   (It used to act on only one event; then I realized that safely dealing
//   with a manual-reset event could be a pain if the first thread finishes
//   before the last one wakes up...)
// * Thread functions are expected to be of the form
//     THREAD_FUNC_DECL function_name(void* arg) {
//       uintptr_t tidx = (uintptr_t)arg;
//       ...
//       while (1) {
//         ... // process current block
//         if (g_is_last_thread_block) {
//           THREAD_RETURN;
//         }
//         THREAD_BLOCK_FINISH(tidx);
//       }
//     }
// * On Linux and OS X, THREAD_BLOCK_FINISH() acquires a mutex, decrements
//   g_thread_active_ct, calls pthread_cond_signal() on
//   g_thread_cur_block_done_condvar iff g_thread_active_ct is now zero, then
//   unconditionally calls pthread_cond_wait on g_thread_start_next_condvar and
//   the mutex.
// * On Windows, THREAD_BLOCK_FINISH() calls SetEvent() on
//   g_thread_cur_block_done_events[tidx], then waits on
//   g_thread_start_next_event[tidx].
// * If the termination variable is set, join_threads2z() waits for all threads
//   to complete, then cleans up all multithreading objects.  Otherwise, on
//   Linux and OS X, it acquires the mutex and calls pthread_cond_wait() on
//   g_thread_cur_block_done_condvar and the mutex; and on Windows, it calls
//   WaitForMultipleObjects() on g_thread_cur_block_done_events[].
//   WaitForMultipleObjects has a 64 object limit, and for now it doesn't seem
//   too important to use a for loop to handle more objects?... well, we can
//   add that if anyone wants it, but for now the Windows thread limit is 64.

uintptr_t g_thread_spawn_ct;
uint32_t g_is_last_thread_block = 0;
#ifdef _WIN32
HANDLE g_thread_start_next_event[kMaxThreads];
HANDLE g_thread_cur_block_done_events[kMaxThreads];
#else
pthread_attr_t g_smallstack_thread_attr;
static pthread_mutex_t g_thread_sync_mutex;
static pthread_cond_t g_thread_cur_block_done_condvar;
static pthread_cond_t g_thread_start_next_condvar;
static uint32_t g_thread_active_ct;

void THREAD_BLOCK_FINISH(__attribute__((unused)) uintptr_t tidx) {
  const uintptr_t initial_spawn_ct = g_thread_spawn_ct;
  pthread_mutex_lock(&g_thread_sync_mutex);
  if (!(--g_thread_active_ct)) {
    pthread_cond_signal(&g_thread_cur_block_done_condvar);
  }
  while (g_thread_spawn_ct == initial_spawn_ct) {
    // spurious wakeup guard
    pthread_cond_wait(&g_thread_start_next_condvar, &g_thread_sync_mutex);
  }
  pthread_mutex_unlock(&g_thread_sync_mutex);
}
#endif
static uint32_t g_thread_mutex_initialized = 0;

void join_threads2z(uint32_t ct, uint32_t is_last_block, pthread_t* threads) {
#ifdef _WIN32
  if (!is_last_block) {
    WaitForMultipleObjects(ct, g_thread_cur_block_done_events, 1, INFINITE);
  } else {
    WaitForMultipleObjects(ct, threads, 1, INFINITE);
    for (uint32_t uii = 0; uii < ct; ++uii) {
      // fix handle leak?
      CloseHandle(threads[uii]);
      
      CloseHandle(g_thread_start_next_event[uii]);
      CloseHandle(g_thread_cur_block_done_events[uii]);
    }
    g_thread_mutex_initialized = 0;
  }
#else
  if (!is_last_block) {
    pthread_mutex_lock(&g_thread_sync_mutex);
    while (g_thread_active_ct) {
      pthread_cond_wait(&g_thread_cur_block_done_condvar, &g_thread_sync_mutex);
    }
    // keep mutex until next block loaded
  } else {
    for (uint32_t uii = 0; uii < ct; ++uii) {
      pthread_join(threads[uii], nullptr);
    }
    // slightly inefficient if there are multiple multithreaded commands being
    // run, but if different commands require different numbers of threads,
    // optimizing this sort of thing away could introduce bugs...
    pthread_mutex_destroy(&g_thread_sync_mutex);
    pthread_cond_destroy(&g_thread_cur_block_done_condvar);
    pthread_cond_destroy(&g_thread_start_next_condvar);
    g_thread_mutex_initialized = 0;
  }
#endif
}

boolerr_t spawn_threads2z(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, uint32_t is_last_block, pthread_t* threads) {
  // start_routine() might need this
  if (g_is_last_thread_block != is_last_block) {
    // might save us an unnecessary memory write that confuses the cache
    // coherency logic?
    g_is_last_thread_block = is_last_block;
  }
#ifdef _WIN32
  if (!g_thread_mutex_initialized) {
    g_thread_spawn_ct = 0;
    g_thread_mutex_initialized = 1;
    for (uintptr_t ulii = 0; ulii < ct; ++ulii) {
      g_thread_start_next_event[ulii] = CreateEvent(nullptr, FALSE, FALSE, nullptr);
      g_thread_cur_block_done_events[ulii] = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    }
    for (uintptr_t ulii = 0; ulii < ct; ++ulii) {
      threads[ulii] = (HANDLE)_beginthreadex(nullptr, kDefaultThreadStack, start_routine, (void*)ulii, 0, nullptr);
      if (!threads[ulii]) {
	if (ulii) {
	  join_threads2z(ulii, is_last_block, threads);
	  if (!is_last_block) {
	    // fix handle leak?
	    for (uintptr_t uljj = 0; uljj < ulii; ++uljj) {
	      CloseHandle(threads[uljj]);
	    }
	  }
	}
	if ((!is_last_block) || (!ulii)) {
	  for (uint32_t uii = 0; uii < ct; ++uii) {
	    CloseHandle(g_thread_start_next_event[uii]);
	    CloseHandle(g_thread_cur_block_done_events[uii]);
	  }
	  g_thread_mutex_initialized = 0;
	}
	return 1;
      }
    }
  } else {
    g_thread_spawn_ct++;
    for (uintptr_t ulii = 0; ulii < ct; ++ulii) {
      SetEvent(g_thread_start_next_event[ulii]);
    }
  }
#else
  if (!is_last_block) {
    g_thread_active_ct = ct;
  }
  if (!g_thread_mutex_initialized) {
    g_thread_spawn_ct = 0; // tidx 0 may need to know modulus
    g_thread_mutex_initialized = 1;
    if (pthread_mutex_init(&g_thread_sync_mutex, nullptr) ||
        pthread_cond_init(&g_thread_cur_block_done_condvar, nullptr) ||
        pthread_cond_init(&g_thread_start_next_condvar, nullptr)) {
      return 1;
    }
    for (uintptr_t ulii = 0; ulii < ct; ++ulii) {
      if (pthread_create(&(threads[ulii]), nullptr, start_routine, (void*)ulii)) {
	if (ulii) {
	  join_threads2z(ulii, is_last_block, threads);
	  if (!is_last_block) {
	    // not worth the trouble of demanding that all callers handle
	    // pthread_create() failure cleanly
	    // (in contrast, error_cleanup_threads2z is relevant when an input
	    // .pgen is malformed, which could happen a lot)
	    for (uintptr_t uljj = 0; uljj < ulii; ++uljj) {
	      pthread_cancel(threads[uljj]);
	    }
	  }
	}
	if ((!is_last_block) || (!ulii)) {
	  pthread_mutex_destroy(&g_thread_sync_mutex);
	  pthread_cond_destroy(&g_thread_cur_block_done_condvar);
	  pthread_cond_destroy(&g_thread_start_next_condvar);
	  g_thread_mutex_initialized = 0;
	}
	return 1;
      }
    }
  } else {
    ++g_thread_spawn_ct;
    // still holding mutex
    pthread_mutex_unlock(&g_thread_sync_mutex);
    pthread_cond_broadcast(&g_thread_start_next_condvar);
  }
#endif
  return 0;
}

void error_cleanup_threads2z(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, pthread_t* threads) {
  if (!spawn_threads2z(start_routine, ct, 1, threads)) {
    join_threads2z(ct, 1, threads);
  }
}


// multithread globals
static const uintptr_t* g_subset_mask = nullptr;
static char** g_item_ids;
static uint32_t* g_id_htable = nullptr;

// currently by item_idx, not item_uidx
static uint32_t* g_item_id_hashes = nullptr;
static uint32_t g_item_ct = 0;
static uint32_t g_id_htable_size = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_item_uidx_starts[16];

THREAD_FUNC_DECL calc_id_hash_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t* subset_mask = g_subset_mask;
  char** item_ids = g_item_ids;
  uint32_t* item_id_hashes = g_item_id_hashes;
  const uint32_t id_htable_size = g_id_htable_size;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t fill_start = round_down_pow2((id_htable_size * ((uint64_t)tidx)) / calc_thread_ct, kInt32PerCacheline);
  uint32_t fill_end;
  if (tidx + 1 < calc_thread_ct) {
    fill_end = round_down_pow2((id_htable_size * (((uint64_t)tidx) + 1)) / calc_thread_ct, kInt32PerCacheline);
  } else {
    fill_end = id_htable_size;
  }
  fill_uint_one(fill_end - fill_start, &(g_id_htable[fill_start]));

  const uint32_t item_ct = g_item_ct;
  const uint32_t item_idx_end = (item_ct * (((uint64_t)tidx) + 1)) / calc_thread_ct;
  uint32_t item_uidx = g_item_uidx_starts[tidx];
  for (uint32_t item_idx = (item_ct * ((uint64_t)tidx)) / calc_thread_ct; item_idx < item_idx_end; ++item_idx, ++item_uidx) {
    next_set_unsafe_ck(subset_mask, &item_uidx);
    const char* sptr = item_ids[item_uidx];
    const uint32_t slen = strlen(sptr);
    item_id_hashes[item_idx] = hashceil(sptr, slen, id_htable_size);
  }
  THREAD_RETURN;
}

pglerr_t populate_id_htable_mt(const uintptr_t* subset_mask, char** item_ids, uintptr_t item_ct, uint32_t store_all_dups, uint32_t id_htable_size, uint32_t thread_ct, uint32_t* id_htable) {
  // Change from plink 1.9: if store_all_dups is false, we don't error out on
  // the first encountered duplicate ID; instead, we just flag it in the hash
  // table.  So if '.' is the only duplicate ID, and it never appears in a
  // variant ID list, plink2 never complains.
  //
  // When store_all_dups is true, additional linked lists are allocated past
  // the end of id_htable to track all raw indexes of duplicate names.
  if (!item_ct) {
    return kPglRetSuccess;
  }
  unsigned char* bigstack_end_mark = g_bigstack_end;
  pglerr_t reterr = kPglRetSuccess;
  {
    // this seems to be a sweet spot
    if (thread_ct > 16) {
      thread_ct = 16;
    }
    if (thread_ct > item_ct / 65536) {
      thread_ct = item_ct / 65536;
      if (!thread_ct) {
	thread_ct = 1;
      }
    }
    if (bigstack_end_alloc_ui(item_ct, &g_item_id_hashes)) {
      goto populate_id_htable_mt_ret_NOMEM;
    }
    g_subset_mask = subset_mask;
    g_item_ids = item_ids;
    g_id_htable = id_htable;
    g_item_ct = item_ct;
    g_id_htable_size = id_htable_size;
    g_calc_thread_ct = thread_ct;
    pthread_t threads[16];
    {
      uint32_t item_uidx = next_set_unsafe(subset_mask, 0);
      uint32_t item_idx = 0;
      g_item_uidx_starts[0] = item_uidx;
      for (uintptr_t tidx = 1; tidx < thread_ct; ++tidx) {
	const uint32_t item_idx_new = (item_ct * ((uint64_t)tidx)) / thread_ct;
	item_uidx = jump_forward_set_unsafe(subset_mask, item_uidx + 1, item_idx_new - item_idx);
	g_item_uidx_starts[tidx] = item_uidx;
	item_idx = item_idx_new;
      }
    }
    if (spawn_threads(calc_id_hash_thread, thread_ct, threads)) {
      goto populate_id_htable_mt_ret_THREAD_CREATE_FAIL;
    }
    calc_id_hash_thread((void*)0);
    join_threads(thread_ct, threads);
    // could also partial-sort and actually fill the hash table in a
    // multithreaded manner, but I'll postpone that for now since it's tricky
    // to make that work with duplicate ID handling, and it also is a
    // substantially smaller bottleneck than hash value computation.
    uint32_t item_uidx = 0;
    if (!store_all_dups) {
      for (uint32_t item_idx = 0; item_idx < item_ct; ++item_uidx, ++item_idx) {
	next_set_unsafe_ck(subset_mask, &item_uidx);
	uint32_t hashval = g_item_id_hashes[item_idx];
	uint32_t cur_htable_entry = id_htable[hashval];
	if (cur_htable_entry == 0xffffffffU) {
	  id_htable[hashval] = item_uidx;
	} else {
	  const char* sptr = item_ids[item_uidx];
	  while (1) {
	    // could also use memcmp, guaranteed to be safe due to where
	    // variant IDs are allocated
	    if (!strcmp(sptr, item_ids[cur_htable_entry & 0x7fffffff])) {
	      if (!(cur_htable_entry >> 31)) {
	        id_htable[hashval] |= 0x80000000U;
	      }
	      break;
	    }
	    if (++hashval == id_htable_size) {
	      hashval = 0;
	    }
	    cur_htable_entry = id_htable[hashval];
	    if (cur_htable_entry == 0xffffffffU) {
	      id_htable[hashval] = item_uidx;
	      break;
	    }
	  }
	}
      }
    } else {
      const uintptr_t cur_bigstack_left = bigstack_left();
      uint32_t max_extra_alloc_m4;
#ifdef __LP64__
      if (cur_bigstack_left >= 0x400000000LLU) {
	// this can never be hit
	max_extra_alloc_m4 = 0xfffffffaU;
      } else {
#endif
	if (cur_bigstack_left < 4 * sizeof(int32_t)) {
	  goto populate_id_htable_mt_ret_NOMEM;
	}
	max_extra_alloc_m4 = (cur_bigstack_left / sizeof(int32_t)) - 4;
#ifdef __LP64__
      }
#endif
      uint32_t extra_alloc = 0;
      uint32_t prev_llidx = 0;
      // needs to be synced with extract_exclude_flag_norange()
      // multithread this?
      uint32_t* htable_dup_base = (uint32_t*)g_bigstack_base;
      for (uint32_t item_idx = 0; item_idx < item_ct; ++item_uidx, ++item_idx) {
	next_set_unsafe_ck(subset_mask, &item_uidx);
	uint32_t hashval = g_item_id_hashes[item_idx];
	uint32_t cur_htable_entry = id_htable[hashval];
	if (cur_htable_entry == 0xffffffffU) {
	  id_htable[hashval] = item_uidx;
	} else {
	  const char* sptr = item_ids[item_uidx];
	  while (1) {
	    const uint32_t cur_dup = cur_htable_entry >> 31;
	    uint32_t prev_uidx;
	    if (cur_dup) {
	      prev_llidx = cur_htable_entry * 2;
	      prev_uidx = htable_dup_base[prev_llidx];
	    } else {
	      prev_uidx = cur_htable_entry;
	    }
	    if (!strcmp(sptr, item_ids[prev_uidx])) {
	      if (extra_alloc > max_extra_alloc_m4) {
		goto populate_id_htable_mt_ret_NOMEM;
	      }
	      // point to linked list entry instead
	      if (!cur_dup) {
		htable_dup_base[extra_alloc] = cur_htable_entry;
		htable_dup_base[extra_alloc + 1] = 0xffffffffU; // list end
		prev_llidx = extra_alloc;
		extra_alloc += 2;
	      }
	      htable_dup_base[extra_alloc] = item_uidx;
	      htable_dup_base[extra_alloc + 1] = prev_llidx;
	      id_htable[hashval] = 0x80000000U | (extra_alloc >> 1);
	      extra_alloc += 2;
	      break; // bugfix
	    }
	    if (++hashval == id_htable_size) {
	      hashval = 0;
	    }
	    cur_htable_entry = id_htable[hashval];
	    if (cur_htable_entry == 0xffffffffU) {
	      id_htable[hashval] = item_uidx;
	      break;
	    }
	  }
	}
      }
      if (extra_alloc) {
	// bugfix: forgot to align this
	bigstack_alloc_raw_rd(extra_alloc * sizeof(int32_t));
      }
    }
  }
  while (0) {
  populate_id_htable_mt_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  populate_id_htable_mt_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  bigstack_end_reset(bigstack_end_mark);
  return reterr;
}

pglerr_t alloc_and_populate_id_htable_mt(const uintptr_t* subset_mask, char** item_ids, uintptr_t item_ct, uint32_t max_thread_ct, uint32_t** id_htable_ptr, uint32_t** htable_dup_base_ptr, uint32_t* id_htable_size_ptr) {
  uint32_t id_htable_size = get_htable_fast_size(item_ct);
  // 4 bytes per variant for hash buffer
  // if store_all_dups, up to 8 bytes per variant in extra_alloc for duplicate
  //   tracking
  const uint32_t store_all_dups = (htable_dup_base_ptr != nullptr);
  const uintptr_t nonhtable_alloc = round_up_pow2(item_ct * sizeof(int32_t), kCacheline) + store_all_dups * round_up_pow2(item_ct * 2 * sizeof(int32_t), kCacheline);
  uintptr_t max_bytes = round_down_pow2(bigstack_left(), kCacheline);
  // force max_bytes >= 5 so leqprime() doesn't fail
  if (nonhtable_alloc + (item_ct + 6) * sizeof(int32_t) > max_bytes) {
    return kPglRetNomem;
  }
  max_bytes -= nonhtable_alloc;
  if (id_htable_size * sizeof(int32_t) > max_bytes) {
    id_htable_size = leqprime((max_bytes / sizeof(int32_t)) - 1);
    const uint32_t min_htable_size = get_htable_min_size(item_ct);
    if (id_htable_size < min_htable_size) {
      id_htable_size = min_htable_size;
    }
  }
  *id_htable_ptr = (uint32_t*)bigstack_alloc_raw_rd(id_htable_size * sizeof(int32_t));
  if (store_all_dups) {
    *htable_dup_base_ptr = &((*id_htable_ptr)[round_up_pow2(id_htable_size, kInt32PerCacheline)]);
  }
  *id_htable_size_ptr = id_htable_size;
  return populate_id_htable_mt(subset_mask, item_ids, item_ct, store_all_dups, id_htable_size, max_thread_ct, *id_htable_ptr);
}

pglerr_t multithread_load_init(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, pgen_file_info_t* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** dosage_present_ptr, dosage_t*** dosage_val_bufs_ptr, uint32_t* read_block_size_ptr, unsigned char** main_loadbufs, pthread_t** threads_ptr, pgen_reader_t*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr) {
  uintptr_t cachelines_avail = bigstack_left() / kCacheline;
  uint32_t read_block_size = kPglVblockSize;
  uint64_t multiread_cacheline_ct;
  while (1) {
    multiread_cacheline_ct = pgfi_multiread_get_cacheline_req(variant_include, pgfip, variant_ct, read_block_size);
    // limit each raw load buffer to 1/4 of remaining workspace
    // if there's an additional per-variant allocation, put it in the same bin
    // as the load buffers
    if ((multiread_cacheline_ct + (((uint64_t)per_variant_xalloc_byte_ct) * read_block_size) / kCacheline) * 4 <= cachelines_avail) {
      break;
    }
    // lots of callers require read_block_size to be either raw_variant_ct or a
    // multiple of kBitsPerVec
#ifdef __LP64__
    if (read_block_size <= kBitsPerVec) {
      return kPglRetNomem;
    }
#else
    if (read_block_size <= kCacheline) {
      return kPglRetNomem;
    }
#endif
    read_block_size /= 2;
  }
#ifndef __LP64__
  if (multiread_cacheline_ct > (kMaxBytesPerIO / kCacheline)) {
    return kPglRetNomem;
  }
#endif
  main_loadbufs[0] = bigstack_alloc_raw(multiread_cacheline_ct * kCacheline);
  main_loadbufs[1] = bigstack_alloc_raw(multiread_cacheline_ct * kCacheline);
  pgfip->block_base = main_loadbufs[0];
  *read_block_size_ptr = read_block_size;
  cachelines_avail -= 2 * (multiread_cacheline_ct + (((uint64_t)per_variant_xalloc_byte_ct) * read_block_size) / kCacheline);
  // reduce calc_thread_ct if necessary
  uint32_t calc_thread_ct = *calc_thread_ct_ptr;
  if (calc_thread_ct > read_block_size) {
    calc_thread_ct = read_block_size;
    *calc_thread_ct_ptr = calc_thread_ct;
  }

  // pgr_pps, threads_ptr, read_variant_uidx_starts_ptr, (*pgr_pps)[tidx],
  //   pgr_alloc; deliberately a slight overestimate
  const uintptr_t pgr_struct_alloc = round_up_pow2(sizeof(pgen_reader_t), kCacheline);
  uintptr_t thread_alloc_cacheline_ct = 1 + 1 + 1 + (pgr_struct_alloc / kCacheline) + pgr_alloc_cacheline_ct + thread_xalloc_cacheline_ct;

  const uint32_t sample_ctcl2 = QUATERCT_TO_CLCT(sample_ct);
  const uint32_t sample_ctcl = BITCT_TO_CLCT(sample_ct);
  
  // todo: increase in multiallelic case
  const uintptr_t dosage_vals_cl = DIV_UP(sample_ct, (kCacheline / sizeof(dosage_t)));
  if (genovecs_ptr) {
    thread_alloc_cacheline_ct += 1 + sample_ctcl2;
    if (dosage_present_ptr) {
      assert(dosage_val_bufs_ptr);
      thread_alloc_cacheline_ct += 2 + sample_ctcl + dosage_vals_cl;
    }
  }
  if (thread_alloc_cacheline_ct * calc_thread_ct > cachelines_avail) {
    if (thread_alloc_cacheline_ct > cachelines_avail) {
      return kPglRetNomem;
    }
    calc_thread_ct = cachelines_avail / thread_alloc_cacheline_ct;
    *calc_thread_ct_ptr = calc_thread_ct;
  }

  const uint32_t array_of_ptrs_alloc = round_up_pow2(calc_thread_ct * sizeof(intptr_t), kCacheline);
  *pgr_pps = (pgen_reader_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
  *threads_ptr = (pthread_t*)bigstack_alloc_raw(array_of_ptrs_alloc);
  *read_variant_uidx_starts_ptr = (uint32_t*)bigstack_alloc_raw_rd(calc_thread_ct * sizeof(int32_t));
  for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
    (*pgr_pps)[tidx] = (pgen_reader_t*)bigstack_alloc_raw(pgr_struct_alloc);
    // pgr_preinit(g_pgr_ptrs[tidx]);
    unsigned char* pgr_alloc = bigstack_alloc_raw(pgr_alloc_cacheline_ct * kCacheline);

    // shouldn't be possible for this to fail
    pgr_init(nullptr, 0, pgfip, (*pgr_pps)[tidx], pgr_alloc);
  }
  if (genovecs_ptr) {
    *genovecs_ptr = (uintptr_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
    if (dosage_present_ptr) {
      *dosage_present_ptr = (uintptr_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
      *dosage_val_bufs_ptr = (dosage_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
    }
    const uintptr_t genovec_alloc = sample_ctcl2 * kCacheline;
    const uintptr_t dosage_present_alloc = sample_ctcl * kCacheline;
    const uintptr_t dosage_vals_alloc = dosage_vals_cl * kCacheline;
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      (*genovecs_ptr)[tidx] = (uintptr_t*)bigstack_alloc_raw(genovec_alloc);
      if (dosage_present_ptr) {
	(*dosage_present_ptr)[tidx] = (uintptr_t*)bigstack_alloc_raw(dosage_present_alloc);
	(*dosage_val_bufs_ptr)[tidx] = (dosage_t*)bigstack_alloc_raw(dosage_vals_alloc);
      }
    }
  }
  return kPglRetSuccess;
}

pglerr_t write_sample_ids(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* outname, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen) {
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto write_sample_ids_ret_OPEN_FAIL;
    }
    char* textbuf = g_textbuf;
    char* write_iter = textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    uintptr_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_ul_unsafe_ck(sample_include, &sample_uidx);
      write_iter = strcpya(write_iter, &(sample_ids[sample_uidx * max_sample_id_blen]));
      if (sids) {
	*write_iter++ = '\t';
	write_iter = strcpya(write_iter, &(sids[sample_uidx * max_sid_blen]));
      }
      append_binary_eoln(&write_iter);
      if (write_iter >= textbuf_flush) {
	if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	  goto write_sample_ids_ret_WRITE_FAIL;
	}
	write_iter = textbuf;
      }
    }
    if (write_iter > textbuf) {
      if (fwrite_checked(textbuf, write_iter - textbuf, outfile)) {
	goto write_sample_ids_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto write_sample_ids_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_sample_ids_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_sample_ids_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
