// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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


#include "plink2_cmdline.h"

#include <unistd.h> // sysconf()

#ifdef __APPLE__
  // needed for sysctl() call
  #include <sys/sysctl.h>
#endif

#include <time.h> // cleanup_logfile()
#include <unistd.h> // getcwd(), gethostname(), sysconf()

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

char g_logbuf[kLogbufSize];

uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
uint32_t g_stderr_written_to = 0;

void logstr(const char* ss) {
  if (!g_debug_on) {
    fputs(ss, g_logfile);
    if (ferror_unlocked(g_logfile)) {
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
      if (ferror_unlocked(g_logfile)) {
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
      token_end = strnul(token_start);
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

boolerr_t fwrite_flush2(char* buf_flush, FILE* outfile, char** write_iter_ptr) {
  char* buf = &(buf_flush[-((int32_t)kMaxMediumLine)]);
  char* buf_end = *write_iter_ptr;
  *write_iter_ptr = buf;
  return fwrite_checked(buf, (uintptr_t)(buf_end - buf), outfile);
}

boolerr_t fclose_flush_null(char* buf_flush, char* write_iter, FILE** outfile_ptr) {
  char* buf = &(buf_flush[-((int32_t)kMaxMediumLine)]);
  if (write_iter != buf) {
    if (fwrite_checked(buf, (uintptr_t)(write_iter - buf), *outfile_ptr)) {
      return 1;
    }
  }
  return fclose_null(outfile_ptr);
}


static const uint32_t kPow10[] =
{1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

#ifdef USE_AVX2
static const unsigned char kLzUintSlenBase[] =
{9, 9, 9, 8,
 8, 8, 7, 7,
 7, 6, 6, 6,
 6, 5, 5, 5,
 4, 4, 4, 3,
 3, 3, 3, 2,
 2, 2, 1, 1,
 1, 0, 0, 0,
 1};  // uint_slen(0) needs to be 1, not zero

uint32_t uint_slen(uint32_t num) {
  const uint32_t lz_ct = _lzcnt_u32(num);
  const uint32_t slen_base = kLzUintSlenBase[lz_ct];
  return slen_base + (num >= kPow10[slen_base]);
}
#else
// could also use something like ((32 - lz_ct) * 77) >> 8, since 77/256 is a
// sufficiently good approximation of ln(2)/ln(10), but that's a bit slower and
// this table doesn't take much space
static const unsigned char kLzUintSlenBase[] =
{9, 9, 9, 8,
 8, 8, 7, 7,
 7, 6, 6, 6,
 6, 5, 5, 5,
 4, 4, 4, 3,
 3, 3, 3, 2,
 2, 2, 1, 1,
 1};

uint32_t uint_slen(uint32_t num) {
  // tried divide-by-10 and divide-by-100 loops, they were slower
  // also tried a hardcoded binary tree, it was better but still slower

  // __builtin_clz(0) is undefined
  if (num < 10) {
    return 1;
  }
  const uint32_t lz_ct = __builtin_clz(num);
  const uint32_t slen_base = kLzUintSlenBase[lz_ct];
  return slen_base + (num >= kPow10[slen_base]);
}
#endif

/*
int32_t strcmp_se(const char* s_read, const char* s_const, uint32_t s_const_slen) {
  return memcmp(s_read, s_const, s_const_slen) || (!is_space_or_eoln(s_read[s_const_slen]));
}
*/

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
  const char* strptr;
  uint32_t orig_idx;
  bool operator<(const struct str_nsort_indexed_deref_struct& rhs) const {
    return (strcmp_natural_uncasted((const unsigned char*)strptr, (const unsigned char*)(rhs.strptr)) < 0);
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


uint32_t copy_and_dedup_sorted_strptrs_to_strbox(const char* const* sorted_strptrs, uintptr_t str_ct, uintptr_t max_str_blen, char* strbox) {
  if (!str_ct) {
    return 0;
  }
  const char* const* sorted_strptrs_iter = sorted_strptrs;
  const char* const* sorted_strptrs_end = &(sorted_strptrs[str_ct]);
  uintptr_t write_idx = 0;
  uint32_t prev_slen = UINT32_MAX;
  const char* prev_str = nullptr;
  do {
    const char* cur_str = *sorted_strptrs_iter++;
    const uint32_t cur_slen = strlen(cur_str);
    if ((cur_slen != prev_slen) || memcmp(cur_str, prev_str, prev_slen)) {
      memcpy(&(strbox[write_idx * max_str_blen]), cur_str, cur_slen + 1);
      ++write_idx;
      prev_str = cur_str;
    }
  } while (sorted_strptrs_iter != sorted_strptrs_end);
  return write_idx;
}


void strptr_arr_sort_main(uintptr_t str_ct, uint32_t use_nsort, str_sort_indexed_deref_t* wkspace_alias) {
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
}

/*
boolerr_t strptr_arr_indexed_sort(const char* const* unsorted_strptrs, uint32_t str_ct, uint32_t use_nsort, uint32_t* id_map) {
  if (str_ct < 2) {
    if (str_ct) {
      id_map[0] = 0;
    }
    return 0;
  }
  if (bigstack_left() < str_ct * sizeof(str_sort_indexed_deref_t)) {
    return 1;
  }
  str_sort_indexed_deref_t* wkspace_alias = (str_sort_indexed_deref_t*)g_bigstack_base;
  for (uint32_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    wkspace_alias[str_idx].strptr = unsorted_strptrs[str_idx];
    wkspace_alias[str_idx].orig_idx = str_idx;
  }
  strptr_arr_sort_main(str_ct, use_nsort, wkspace_alias);
  for (uint32_t str_idx = 0; str_idx < str_ct; ++str_idx) {
    id_map[str_idx] = wkspace_alias[str_idx].orig_idx;
  }
  bigstack_reset(wkspace_alias);
  return 0;
}
*/


uint32_t uint32arr_greater_than(const uint32_t* sorted_uint32_arr, uint32_t arr_length, uint32_t uii) {
  // (strangely, this seems to be equal to or better than std::lower_bound with
  // -O2 optimization, but can become much slower with -O3?)

  // assumes arr_length is nonzero, and sorted_uint32_arr is in nondecreasing
  // order.  (useful for searching variant_bps[].)
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

uint32_t param_count(const char* const* argvc, uint32_t argc, uint32_t flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any nonnumeric parameter beginning with "-" as
  // optional.
  ++flag_idx;
  uint32_t cur_idx = flag_idx;
  while ((cur_idx < argc) && (!is_flag(argvc[cur_idx]))) {
    ++cur_idx;
  }
  return cur_idx - flag_idx;
}

boolerr_t enforce_param_ct_range(const char* flag_name, uint32_t param_ct, uint32_t min_ct, uint32_t max_ct) {
  if (param_ct > max_ct) {
    if (max_ct > min_ct) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s accepts at most %u parameter%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "Error: %s only accepts %u parameter%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    }
    return 1;
  }
  if (param_ct >= min_ct) {
    return 0;
  }
  if (min_ct == 1) {
    snprintf(g_logbuf, kLogbufSize, "Error: Missing %s parameter.\n", flag_name);
  } else {
    snprintf(g_logbuf, kLogbufSize, "Error: %s requires %s%u parameters.\n", flag_name, (min_ct < max_ct)? "at least " : "", min_ct);
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

/*
boolerr_t push_llstr_counted(const char* ss, uint32_t slen, ll_str_t** ll_stack_ptr) {
  ll_str_t* new_llstr;
  if (pgl_malloc(sizeof(ll_str_t) + slen + 1, &new_llstr)) {
    return 1;
  }
  new_llstr->next = *ll_stack_ptr;
  memcpy(new_llstr->ss, ss, slen);
  new_llstr->ss[slen] = '\0';
  *ll_stack_ptr = new_llstr;
  return 0;
}

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
static inline boolerr_t scanadv_uint_capped_finish(uint64_t cap, const char** ss_ptr, uint32_t* valp) {
  const unsigned char* ss = (const unsigned char*)(*ss_ptr);
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
  *ss_ptr = (const char*)(&(ss[-1]));
  return 0;
}

boolerr_t scanadv_posint_capped(uint64_t cap, const char** ss_ptr, uint32_t* valp) {
  const unsigned char* ss = (const unsigned char*)(*ss_ptr);
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
  *ss_ptr = (const char*)ss;
  return scanadv_uint_capped_finish(cap, ss_ptr, valp);
}

boolerr_t scanadv_uint_capped(uint64_t cap, const char** ss_ptr, uint32_t* valp) {
  const unsigned char* ss = (const unsigned char*)(*ss_ptr);
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
      *ss_ptr = (const char*)ss;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
    }
    // accept leading '+'
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (*valp >= 10) {
      return 1;
    }
  }
  *ss_ptr = (const char*)ss;
  return scanadv_uint_capped_finish(cap, ss_ptr, valp);
}
#else
boolerr_t scanadv_posint_capped32(uint32_t cap_div_10, uint32_t cap_mod_10, const char** ss_ptr, uint32_t* valp) {
  const unsigned char* ss = (const unsigned char*)ss_ptr;
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
      *ss_ptr = (const char*)(&(ss[-1]));
      return 0;
    }
    if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

boolerr_t scanadv_uint_capped32(uint32_t cap_div_10, uint32_t cap_mod_10, const char** ss_ptr, uint32_t* valp) {
  const unsigned char* ss = (const unsigned char*)ss_ptr;
  uint32_t val = (uint32_t)(*ss++) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      if ((val != 0xfffffffd) || (*ss != '0')) {
        return 1;
      }
      while (*(++ss) == '0');
      *valp = 0;
      *ss_ptr = (const char*)ss;
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
      *ss_ptr = (const char*)(&(ss[-1]));
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

CXXCONST_CP scanadv_double(const char* ss, double* valp) {
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
  const char* dot_ptr;
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
    const char* last_sig_fig_ptr = ss;
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
        return (CXXCONST_CP)ss;
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
    const char* last_sig_fig_ptr = ss;
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
        return (CXXCONST_CP)ss;
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
    return (CXXCONST_CP)ss;
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
  return (CXXCONST_CP)ss;
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

CXXCONST_CP comma_or_space_next_token_mult(const char* sptr, uint32_t ct, uint32_t comma_delim) {
  assert(ct);
  if (!comma_delim) {
    return (CXXCONST_CP)next_token_mult(sptr, ct);
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
        return (CXXCONST_CP)sptr;
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
*/

uint32_t count_and_measure_multistr(const char* multistr, uintptr_t* max_blen_ptr) {
  uint32_t ct = 0;
  uintptr_t max_blen = *max_blen_ptr;
  while (*multistr) {
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

boolerr_t count_and_measure_multistr_reverse_alloc(const char* multistr, uintptr_t max_str_ct, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr, const char*** strptr_arrp) {
  // assumes multistr is nonempty
  assert(multistr[0]);
  uintptr_t ct = 0;
  uintptr_t max_blen = *max_blen_ptr;
  const char** strptr_arr_iter = *strptr_arrp;
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

boolerr_t multistr_to_strbox_dedup_arena_alloc(unsigned char* arena_top, const char* multistr, unsigned char** arena_bottom_ptr, char** sorted_strbox_ptr, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr) {
  const char** strptr_arr = (const char**)arena_top;
  uintptr_t max_str_blen = 0;
  uint32_t str_ct;
  if (count_and_measure_multistr_reverse_alloc(multistr, (arena_top - (*arena_bottom_ptr)) / sizeof(intptr_t), &str_ct, &max_str_blen, &strptr_arr)) {
    return 1;
  }
  const uintptr_t strbox_byte_ct = round_up_pow2(str_ct * max_str_blen, kCacheline);
  if ((((uintptr_t)strptr_arr) - ((uintptr_t)(*arena_bottom_ptr))) < strbox_byte_ct) {
    return 1;
  }
  strptr_arr_sort(str_ct, TO_CONSTCPP(strptr_arr));
  *sorted_strbox_ptr = (char*)(*arena_bottom_ptr);
  str_ct = copy_and_dedup_sorted_strptrs_to_strbox(strptr_arr, str_ct, max_str_blen, *sorted_strbox_ptr);
  *arena_bottom_ptr += round_up_pow2(str_ct * max_str_blen, kCacheline);
  *str_ct_ptr = str_ct;
  *max_blen_ptr = max_str_blen;
  return 0;
}

// number-to-string encoders

const uint16_t kDigitPair[100] = {
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

char* uitoa_z5(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  *start++ = '0' + quotient;
  return uitoa_z4(uii - 10000 * quotient, start);
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

// "prob" means that the number is guaranteed to be in [0, 1].
// no leading space is printed.  trailing zeroes (/decimal point) are erased
//   iff there is equality to ~13 decimal places.
char* dtoa_f_probp6_spaced(double dxx, char* start) {
  double dxx_10_6 = dxx * 1000000;
  const uint32_t dec_digits = double_bround(dxx_10_6, kBankerRound8);
  *start++ = '0' + (dec_digits == 1000000);
  *start++ = '.';
  start = uitoa_z6(dec_digits, start);
  if (fabs(dxx_10_6 - (double)((int32_t)dec_digits)) >= 0.00000005) {
    return start;
  }
  trailing_zeroes_to_spaces(start);
  return start;
}

char* dtoa_f_probp6_clipped(double dxx, char* start) {
  double dxx_10_6 = dxx * 1000000;
  const uint32_t dec_digits = double_bround(dxx_10_6, kBankerRound8);
  *start++ = '0' + (dec_digits == 1000000);
  *start++ = '.';
  start = uitoa_z6(dec_digits, start);
  if (fabs(dxx_10_6 - (double)((int32_t)dec_digits)) >= 0.00000005) {
    return start;
  }
  return clip_trailing_zeroes(start);
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
CSINLINE uint32_t float_round(float fxx) {
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


void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp) {
  // Enables fast integer division by a constant not known until runtime.  See
  // http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html .
  // Assumes divisor is not zero, of course.
  // May want to populate a struct instead.
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
    bitarr[maj_start] = bzhi(bitarr[maj_start], start_idx % kBitsPerWord);
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

// floor permitted to be -1, though not smaller than that.
int32_t prev_set(const uintptr_t* bitarr, uint32_t loc, int32_t floor) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uint32_t remainder = loc % kBitsPerWord;
  uintptr_t ulii;
  if (remainder) {
    ulii = bzhi(*bitarr_ptr, remainder);
    if (ulii) {
      const uint32_t set_bit_loc = (loc | (kBitsPerWord - 1)) - CLZLU(ulii);
      return MAXV(((int32_t)set_bit_loc), floor);
    }
  }
  const uintptr_t* bitarr_last = &(bitarr[((uint32_t)(floor + 1)) / kBitsPerWord]);
  do {
    if (bitarr_ptr <= bitarr_last) {
      return floor;
    }
    ulii = *(--bitarr_ptr);
  } while (!ulii);
  const uint32_t set_bit_loc = (uint32_t)(((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + kBitsPerWord - 1 - CLZLU(ulii));
  return MAXV(((int32_t)set_bit_loc), floor);
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
    *bitarr = bzhi(~(*bitarr), trailing_bit_ct);
  }
}

void bitarr_invert_copy(const uintptr_t* __restrict source_bitarr, uintptr_t bit_ct, uintptr_t* __restrict target_bitarr) {
  const uintptr_t* source_bitarr_stop = &(source_bitarr[bit_ct / kBitsPerWord]);
  while (source_bitarr < source_bitarr_stop) {
    *target_bitarr++ = ~(*source_bitarr++);
  }
  const uint32_t trailing_bit_ct = bit_ct % kBitsPerWord;
  if (trailing_bit_ct) {
    *target_bitarr = bzhi(~(*source_bitarr), trailing_bit_ct);
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
    main_bitvec[base_idx + 1] |= arg_bitvec[base_idx + 1];
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

/*
void bitvec_ornot(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec) {
  // main_bitvec := main_bitvec OR (~arg_bitvec)
#ifdef __LP64__
  vul_t* main_bitvvec_iter = (vul_t*)main_bitvec;
  const vul_t* arg_bitvvec_iter = (const vul_t*)arg_bitvec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    // todo: verify the compiler isn't dumb here
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] |= ~arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] |= ~arg_bitvec[base_idx + 1]
  }
  #endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] |= ~arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] |= ~arg_bitvec[widx];
  }
#endif
}
*/


int32_t get_variant_uidx_without_htable(const char* idstr, const char* const* variant_ids, const uintptr_t* variant_include, uint32_t variant_ct) {
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


// MurmurHash3, from
// https://code.google.com/p/smhasher/source/browse/trunk/MurmurHash3.cpp
CSINLINE uint32_t rotl32(uint32_t x, int8_t r) {
  return (x << r) | (x >> (32 - r));
}

static inline uint32_t getblock32(const uint32_t* p, int i) {
  return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

CSINLINE2 uint32_t fmix32(uint32_t h) {
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


/*
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
*/

boolerr_t htable_good_size_alloc(uint32_t item_ct, uintptr_t bytes_avail, uint32_t** htable_ptr, uint32_t* htable_size_ptr) {
  bytes_avail &= (~(kCacheline - k1LU));
  uint32_t htable_size = get_htable_fast_size(item_ct);
  if (htable_size > bytes_avail / sizeof(int32_t)) {
    if (!bytes_avail) {
      return 1;
    }
    htable_size = bytes_avail / sizeof(int32_t);
    // htable_size = leqprime((bytes_avail / sizeof(int32_t)) - 1);
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
      if (cur_htable_entry == UINT32_MAX) {
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
      if (cur_htable_entry == UINT32_MAX) {
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

uint32_t id_htable_find(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size) {
  // returns UINT32_MAX on failure
  uint32_t hashval = hashceil(cur_id, cur_id_slen, id_htable_size);
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == UINT32_MAX) || (!strcmp(cur_id, item_ids[cur_htable_idval]))) {
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
    if ((cur_htable_idval == UINT32_MAX) || (!memcmp(cur_id, &(strbox[cur_htable_idval * max_str_blen]), cur_id_blen))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t variant_id_dupflag_htable_find(const char* idbuf, const char* const* variant_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen) {
  // assumes duplicate variant IDs are flagged, but full variant_uidx linked
  // lists are not stored
  // idbuf does not need to be null-terminated (note that this is currently
  // achieved in a way that forces variant_ids[] entries to not be too close
  // to the end of bigstack, otherwise memcmp behavior is potentially
  // undefined)
  // returns UINT32_MAX on failure, value with bit 31 set on duplicate
  if (cur_id_slen > max_id_slen) {
    return UINT32_MAX;
  }
  uint32_t hashval = hashceil(idbuf, cur_id_slen, id_htable_size);
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == UINT32_MAX) || ((!memcmp(idbuf, variant_ids[cur_htable_idval & 0x7fffffff], cur_id_slen)) && (!variant_ids[cur_htable_idval & 0x7fffffff][cur_id_slen]))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t variant_id_dup_htable_find(const char* idbuf, const char* const* variant_ids, const uint32_t* id_htable, const uint32_t* htable_dup_base, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen, uint32_t* llidx_ptr) {
  // Permits duplicate entries.  Similar to plink 1.9
  // extract_exclude_process_token().
  // - Returns UINT32_MAX on failure (llidx currently unset in that case),
  //   otherwise returns the index of the first match (which will have the
  //   highest index, due to how the linked list is constructed)
  // - Sets llidx to UINT32_MAX if not a duplicate, otherwise it's the
  //   position in htable_dup_base[] of the next {variant_uidx, next_llidx}
  //   linked list entry.
  // - idbuf does not need to be null-terminated.
  if (cur_id_slen > max_id_slen) {
    return UINT32_MAX;
  }
  uint32_t hashval = hashceil(idbuf, cur_id_slen, id_htable_size);
  while (1) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    const uint32_t cur_dup = cur_htable_idval >> 31;
    uint32_t cur_llidx;
    uint32_t variant_uidx;
    if (cur_dup) {
      // UINT32_MAX empty-entry code has high bit set, so only need to check
      // here
      if (cur_htable_idval == UINT32_MAX) {
        return UINT32_MAX;
      }
      cur_llidx = cur_htable_idval << 1;
      variant_uidx = htable_dup_base[cur_llidx];
    } else {
      cur_llidx = UINT32_MAX;
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

CXXCONST_CP scan_for_duplicate_ids(const char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_blen) {
  --id_ct;
  for (uintptr_t id_idx = 0; id_idx < id_ct; ++id_idx) {
    if (!strcmp(&(sorted_ids[id_idx * max_id_blen]), &(sorted_ids[(id_idx + 1) * max_id_blen]))) {
      return (CXXCONST_CP)(&(sorted_ids[id_idx * max_id_blen]));
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
  // Note that this DOES still perform a "stack" allocation.
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

uintptr_t fwdsearch_str_lb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx, uintptr_t cur_idx) {
  uintptr_t next_incr = 1;
  uintptr_t start_idx = cur_idx;
  while (cur_idx < end_idx) {
    if (memcmp(idbuf, &(sorted_strbox[cur_idx * max_id_blen]), cur_id_slen) <= 0) {
      end_idx = cur_idx;
      break;
    }
    start_idx = cur_idx + 1;
    cur_idx += next_incr;
    next_incr *= 2;
  }
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
    if (idx1 < offset) {
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

pglerr_t string_range_list_to_bitarr(const char* header_line, const range_list_t* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t* bitarr, int32_t* __restrict seen_idxs) {
  // bitarr assumed to be zero-initialized
  // if fixed_len is zero, header_line is assumed to be a list of
  // space-delimited unequal-length names
  assert(token_ct);
  assert(!popcount_longs(bitarr, BITCT_TO_WORDCT(token_ct)));
  pglerr_t reterr = kPglRetSuccess;
  {
    const char* header_line_iter = header_line;
    const uintptr_t name_ct = range_list_ptr->name_ct;
    const uintptr_t max_id_blen = range_list_ptr->name_max_blen;
    uint32_t item_idx = 0;
    while (1) {
      const char* token_end = comma_or_space_token_end(header_line_iter, comma_delim);
      uint32_t cmdline_pos;
      if (!sorted_idbox_find(header_line_iter, sorted_ids, id_map, (uintptr_t)(token_end - header_line_iter), max_id_blen, name_ct, &cmdline_pos)) {
        if (seen_idxs[cmdline_pos] != -1) {
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate --%s token in %s.\n", range_list_flag, file_descrip);
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
    snprintf(g_logbuf, kLogbufSize, "Error: Missing --%s token in %s.\n", range_list_flag, file_descrip);
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

pglerr_t string_range_list_to_bitarr_alloc(const char* header_line, const range_list_t* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t** bitarr_ptr) {
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
    ct += popcount_long(bzhi(bitvec[end_idxl], end_idxlr));
  }
  return ct;
}

#ifdef USE_AVX2
uintptr_t popcount_avx2_intersect(const vul_t* __restrict vvec1_iter, const vul_t* __restrict vvec2_iter, uintptr_t vec_ct) {
  // See popcnt_avx2() in libpopcnt.  vec_ct must be a multiple of 16.
  vul_t cnt = vul_setzero();
  vul_t ones = vul_setzero();
  vul_t twos = vul_setzero();
  vul_t fours = vul_setzero();
  vul_t eights = vul_setzero();
  for (uintptr_t vec_idx = 0; vec_idx < vec_ct; vec_idx += 16) {
    vul_t twos_a = csa256(vvec1_iter[vec_idx + 0] & vvec2_iter[vec_idx + 0], vvec1_iter[vec_idx + 1] & vvec2_iter[vec_idx + 1], &ones);
    vul_t twos_b = csa256(vvec1_iter[vec_idx + 2] & vvec2_iter[vec_idx + 2], vvec1_iter[vec_idx + 3] & vvec2_iter[vec_idx + 3], &ones);
    vul_t fours_a = csa256(twos_a, twos_b, &twos);

    twos_a = csa256(vvec1_iter[vec_idx + 4] & vvec2_iter[vec_idx + 4], vvec1_iter[vec_idx + 5] & vvec2_iter[vec_idx + 5], &ones);
    twos_b = csa256(vvec1_iter[vec_idx + 6] & vvec2_iter[vec_idx + 6], vvec1_iter[vec_idx + 7] & vvec2_iter[vec_idx + 7], &ones);
    vul_t fours_b = csa256(twos_a, twos_b, &twos);
    const vul_t eights_a = csa256(fours_a, fours_b, &fours);

    twos_a = csa256(vvec1_iter[vec_idx + 8] & vvec2_iter[vec_idx + 8], vvec1_iter[vec_idx + 9] & vvec2_iter[vec_idx + 9], &ones);
    twos_b = csa256(vvec1_iter[vec_idx + 10] & vvec2_iter[vec_idx + 10], vvec1_iter[vec_idx + 11] & vvec2_iter[vec_idx + 11], &ones);
    fours_a = csa256(twos_a, twos_b, &twos);

    twos_a = csa256(vvec1_iter[vec_idx + 12] & vvec2_iter[vec_idx + 12], vvec1_iter[vec_idx + 13] & vvec2_iter[vec_idx + 13], &ones);
    twos_b = csa256(vvec1_iter[vec_idx + 14] & vvec2_iter[vec_idx + 14], vvec1_iter[vec_idx + 15] & vvec2_iter[vec_idx + 15], &ones);
    fours_b = csa256(twos_a, twos_b, &twos);
    const vul_t eights_b = csa256(fours_a, fours_b, &fours);
    const vul_t sixteens = csa256(eights_a, eights_b, &eights);
    cnt = cnt + popcount_avx2_single(sixteens);
  }
  cnt = vul_lshift(cnt, 4);
  cnt = cnt + vul_lshift(popcount_avx2_single(eights), 3);
  cnt = cnt + vul_lshift(popcount_avx2_single(fours), 2);
  cnt = cnt + vul_lshift(popcount_avx2_single(twos), 1);
  cnt = cnt + popcount_avx2_single(ones);
  return hsum64(cnt);
}

uintptr_t popcount_longs_intersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct) {
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  const uintptr_t block_ct = word_ct / (16 * kWordsPerVec);
  uintptr_t tot = 0;
  if (block_ct) {
    tot = popcount_avx2_intersect((const vul_t*)bitvec1_iter, (const vul_t*)bitvec2_iter, block_ct * 16);
    bitvec1_iter = &(bitvec1_iter[block_ct * (16 * kWordsPerVec)]);
    bitvec2_iter = &(bitvec2_iter[block_ct * (16 * kWordsPerVec)]);
  }
  while (bitvec1_iter < bitvec1_end) {
    tot += popcount_long((*bitvec1_iter++) & (*bitvec2_iter++));
  }
  return tot;
}
#else // !USE_AVX2
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
#endif // !USE_AVX2

#ifdef USE_SSE42
void popcount_longs_intersect_3val(const uintptr_t* __restrict bitvec1, const uintptr_t* __restrict bitvec2, uint32_t word_ct, uint32_t* __restrict popcount1_ptr, uint32_t* __restrict popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr) {
  uint32_t ct1 = 0;
  uint32_t ct2 = 0;
  uint32_t ct3 = 0;
  for (uint32_t widx = 0; widx < word_ct; ++widx) {
    const uintptr_t word1 = bitvec1[widx];
    const uintptr_t word2 = bitvec2[widx];
    ct1 += popcount_long(word1);
    ct2 += popcount_long(word2);
    ct3 += popcount_long(word1 & word2);
  }
  *popcount1_ptr = ct1;
  *popcount2_ptr = ct2;
  *popcount_intersect_ptr = ct3;
}
#else
static inline void popcount_vecs_intersect_3val(const vul_t* __restrict vvec1_iter, const vul_t* __restrict vvec2_iter, uint32_t vec_ct, uint32_t* __restrict popcount1_ptr, uint32_t* popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr) {
  // ct must be a multiple of 3.
  assert(!(vec_ct % 3));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);
  uint32_t ct1 = 0;
  uint32_t ct2 = 0;
  uint32_t ct3 = 0;
  while (1) {
    univec_t acc1;
    univec_t acc2;
    univec_t acc3;
    acc1.vi = vul_setzero();
    acc2.vi = vul_setzero();
    acc3.vi = vul_setzero();
    const vul_t* vvec1_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
        *popcount1_ptr = ct1;
        *popcount2_ptr = ct2;
        *popcount_intersect_ptr = ct3;
        return;
      }
      vvec1_stop = &(vvec1_iter[vec_ct]);
      vec_ct = 0;
    } else {
      vvec1_stop = &(vvec1_iter[30]);
      vec_ct -= 30;
    }
    do {
      vul_t count1a = *vvec1_iter++;
      vul_t count2a = *vvec2_iter++;
      vul_t count3a = count1a & count2a;
      vul_t count1b = *vvec1_iter++;
      vul_t count2b = *vvec2_iter++;
      vul_t count3b = count1b & count2b;
      vul_t half1a = *vvec1_iter++;
      vul_t half2a = *vvec2_iter++;
      const vul_t half1b = vul_rshift(half1a, 1) & m1;
      const vul_t half2b = vul_rshift(half2a, 1) & m1;
      half1a = half1a & m1;
      half2a = half2a & m1;

      count1a = count1a - (vul_rshift(count1a, 1) & m1);
      count2a = count2a - (vul_rshift(count2a, 1) & m1);
      count3a = count3a - (vul_rshift(count3a, 1) & m1);
      count1b = count1b - (vul_rshift(count1b, 1) & m1);
      count2b = count2b - (vul_rshift(count2b, 1) & m1);
      count3b = count3b - (vul_rshift(count3b, 1) & m1);
      count1a = count1a + half1a;
      count2a = count2a + half2a;
      count3a = count3a + (half1a & half2a);
      count1b = count1b + half1b;
      count2b = count2b + half2b;
      count3b = count3b + (half1b & half2b);

      count1a = (count1a & m2) + (vul_rshift(count1a, 2) & m2);
      count2a = (count2a & m2) + (vul_rshift(count2a, 2) & m2);
      count3a = (count3a & m2) + (vul_rshift(count3a, 2) & m2);
      count1a = count1a + (count1b & m2) + (vul_rshift(count1b, 2) & m2);
      count2a = count2a + (count2b & m2) + (vul_rshift(count2b, 2) & m2);
      count3a = count3a + (count3b & m2) + (vul_rshift(count3b, 2) & m2);
      acc1.vi = acc1.vi + (count1a & m4) + (vul_rshift(count1a, 4) & m4);
      acc2.vi = acc2.vi + (count2a & m4) + (vul_rshift(count2a, 4) & m4);
      acc3.vi = acc3.vi + (count3a & m4) + (vul_rshift(count3a, 4) & m4);
    } while (vvec1_iter < vvec1_stop);
    const vul_t m8 = VCONST_UL(kMask00FF);
    acc1.vi = (acc1.vi & m8) + (vul_rshift(acc1.vi, 8) & m8);
    acc2.vi = (acc2.vi & m8) + (vul_rshift(acc2.vi, 8) & m8);
    acc3.vi = (acc3.vi & m8) + (vul_rshift(acc3.vi, 8) & m8);
    ct1 += univec_hsum_16bit(acc1);
    ct2 += univec_hsum_16bit(acc2);
    ct3 += univec_hsum_16bit(acc3);
  }
}

void popcount_longs_intersect_3val(const uintptr_t* __restrict bitvec1, const uintptr_t* __restrict bitvec2, uint32_t word_ct, uint32_t* __restrict popcount1_ptr, uint32_t* __restrict popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr) {
  const uint32_t trivec_ct = word_ct / (3 * kWordsPerVec);
  uint32_t ct1;
  uint32_t ct2;
  uint32_t ct3;
  popcount_vecs_intersect_3val((const vul_t*)bitvec1, (const vul_t*)bitvec2, trivec_ct * 3, &ct1, &ct2, &ct3);
  const uint32_t words_consumed = trivec_ct * (3 * kWordsPerVec);
  bitvec1 = &(bitvec1[words_consumed]);
  bitvec2 = &(bitvec2[words_consumed]);
  const uint32_t remainder = word_ct - words_consumed;
  for (uint32_t widx = 0; widx < remainder; ++widx) {
    const uintptr_t word1 = bitvec1[widx];
    const uintptr_t word2 = bitvec2[widx];
    ct1 += popcount_long(word1);
    ct2 += popcount_long(word2);
    ct3 += popcount_long(word1 & word2);
  }
  *popcount1_ptr = ct1;
  *popcount2_ptr = ct2;
  *popcount_intersect_ptr = ct3;
}
#endif

uint32_t are_all_bits_zero(const uintptr_t* bitarr, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / kBitsPerWord;
  const uintptr_t start_idxlr = start_idx & (kBitsPerWord - 1);
  const uintptr_t end_idxl = end_idx / kBitsPerWord;
  const uintptr_t end_idxlr = end_idx & (kBitsPerWord - 1);
  if (start_idxl == end_idxl) {
    return !(bitarr[start_idxl] & ((k1LU << end_idxlr) - (k1LU << start_idxlr)));
  }
  if (start_idxlr) {
    if (bitarr[start_idxl++] >> start_idxlr) {
      return 0;
    }
  }
  for (; start_idxl < end_idxl; ++start_idxl) {
    if (bitarr[start_idxl]) {
      return 0;
    }
  }
  if (!end_idxlr) {
    return 1;
  }
  return !bzhi(bitarr[end_idxl], end_idxlr);
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
      while (--forward_ct) {
        uljj &= uljj - 1;
      }
      ulkk = CTZLU(uljj);
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
#ifdef USE_AVX2
  while (forward_ct > kBitsPerWord * (16 * kWordsPerVec)) {
    uljj = ((forward_ct - 1) / (kBitsPerWord * (16 * kWordsPerVec))) * 16;
    ulkk = popcount_avx2(vptr, uljj);
    vptr = &(vptr[uljj]);
    forward_ct -= ulkk;
  }
#else
  while (forward_ct > kBitsPerWord * (3 * kWordsPerVec)) {
    uljj = ((forward_ct - 1) / (kBitsPerWord * (3 * kWordsPerVec))) * 3;
    ulkk = popcount_vecs_old(vptr, uljj);
    vptr = &(vptr[uljj]);
    forward_ct -= ulkk;
  }
#endif
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

void compute_partition_aligned(const uintptr_t* variant_include, uint32_t orig_thread_ct, uint32_t first_variant_uidx, uint32_t cur_variant_idx, uint32_t cur_variant_ct, uint32_t alignment, uint32_t* variant_uidx_starts, uint32_t* vidx_starts) {
  // Minimize size of the largest chunk, under the condition that all
  // intermediate variant_idx values are divisible by alignment.
  //
  // There are three possibilities:
  // 1. The straightforward solution one gets from rounding cur_variant_idx
  //    down, and (cur_variant_idx + variant_ct) up, is optimal.  Call this
  //    chunk size C.
  // 2. C - leading_idx_ct is ideal.
  // 3. C - trailing_idx_ct is ideal.

  assert(cur_variant_ct);
  const uint32_t variant_idx_end = cur_variant_idx + cur_variant_ct;
  const uint32_t log2_align = __builtin_ctz(alignment);
  const uint32_t leading_idx_ct = cur_variant_idx % alignment;
  const uint32_t trailing_idx_ct = (-variant_idx_end) % alignment;
  const uint32_t block_ct = (cur_variant_ct + leading_idx_ct + trailing_idx_ct) >> log2_align;
  uint32_t cur_variant_uidx_start = next_set_unsafe(variant_include, first_variant_uidx);
  variant_uidx_starts[0] = cur_variant_uidx_start;
  vidx_starts[0] = cur_variant_idx;
  const uint32_t thread_ct = MINV(orig_thread_ct, block_ct);
  if (thread_ct > 1) {
    const uint32_t std_blocks_per_thread = 1 + (block_ct - 1) / thread_ct;
    // Possibilities 2 and 3 are only live if
    //   block_ct == (std_blocks_per_thread - 1) * thread_ct + 1, or
    //   block_ct == (std_blocks_per_thread - 1) * thread_ct + 2,
    // In the first situation, the best solution is
    //    min(possibility 2, possibility 3).
    // In the second situation, the best solution is
    //    max(possibility 2, possibility 3).
    uint32_t central_variant_ct = std_blocks_per_thread * alignment;
    uint32_t first_block_variant_ct = central_variant_ct - leading_idx_ct;
    const uint32_t remainder_m1 = block_ct - (std_blocks_per_thread - 1) * thread_ct - 1;
    if (remainder_m1 <= 1) {
      central_variant_ct -= alignment;
      if ((!remainder_m1) && (leading_idx_ct < trailing_idx_ct)) {
        first_block_variant_ct -= alignment;
      }
    }
    cur_variant_uidx_start = jump_forward_set_unsafe(variant_include, cur_variant_uidx_start + 1, first_block_variant_ct);
    cur_variant_idx += first_block_variant_ct;
    variant_uidx_starts[1] = cur_variant_uidx_start;
    vidx_starts[1] = cur_variant_idx;
    for (uint32_t tidx = 2; tidx < thread_ct; ++tidx) {
      cur_variant_uidx_start = jump_forward_set_unsafe(variant_include, cur_variant_uidx_start + 1, central_variant_ct);
      cur_variant_idx += central_variant_ct;
      // bugfix (14 Nov 2017): this decrement was in the wrong place
      if (tidx == remainder_m1) {
        central_variant_ct -= alignment;
      }
      variant_uidx_starts[tidx] = cur_variant_uidx_start;
      vidx_starts[tidx] = cur_variant_idx;
    }
  }
  if (thread_ct < orig_thread_ct) {
    uint32_t last_vidx_ct = variant_idx_end - vidx_starts[thread_ct - 1];
    cur_variant_uidx_start = jump_forward_set_unsafe(variant_include, cur_variant_uidx_start + 1, last_vidx_ct);
    for (uint32_t tidx = thread_ct; tidx < orig_thread_ct; ++tidx) {
      variant_uidx_starts[tidx] = cur_variant_uidx_start;
    }
    for (uint32_t tidx = thread_ct; tidx < orig_thread_ct; ++tidx) {
      vidx_starts[tidx] = variant_idx_end;
    }
  }
  vidx_starts[orig_thread_ct] = variant_idx_end;
}

boolerr_t parse_next_range(const char* const* argvc, uint32_t param_ct, char range_delim, uint32_t* cur_param_idx_ptr, const char** cur_arg_pptr, const char** range_start_ptr, uint32_t* rs_len_ptr, const char** range_end_ptr, uint32_t* re_len_ptr) {
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
  const char* cur_arg_ptr = *cur_arg_pptr;
  while (1) {
    char cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
        *range_start_ptr = nullptr;
        return 0;
      }
      cur_arg_ptr = argvc[cur_param_idx];
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
  if (((*rs_len_ptr) > kMaxIdSlen) || (!cc) || (cc == ',') || (cc == range_delim)) {
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
  if ((*re_len_ptr) > kMaxIdSlen) {
    return 1;
  }
  *cur_arg_pptr = cur_arg_ptr;
  return 0;
}

pglerr_t parse_name_ranges(const char* const* argvc, const char* errstr_append, uint32_t param_ct, uint32_t require_posint, char range_delim, range_list_t* range_list_ptr) {
  uint32_t name_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t name_max_blen = 0;
  const char* cur_arg_ptr;
  const char* range_start;
  uint32_t rs_len;
  const char* range_end;
  uint32_t re_len;
  unsigned char* cur_name_starts_range;
  uint32_t last_val;
  uint32_t cur_val;
  // two passes.  first pass: count parameters, determine name_max_blen;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argvc[1];
    while (1) {
      if (parse_next_range(argvc, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
        LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argvc[0], argvc[cur_param_idx]);
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
    LOGERRPRINTF("Error: %s requires at least one value.\n%s", argvc[0], errstr_append);
    return kPglRetInvalidCmdline;
  }
  range_list_ptr->name_max_blen = ++name_max_blen;
  range_list_ptr->name_ct = name_ct;
  if (pgl_malloc(name_ct * (((uintptr_t)name_max_blen) + 1), &range_list_ptr->names)) {
    return kPglRetNomem;
  }
  range_list_ptr->starts_range = (unsigned char*)(&(range_list_ptr->names[name_ct * ((uintptr_t)name_max_blen)]));
  char* cur_name_str = range_list_ptr->names;
  cur_name_starts_range = range_list_ptr->starts_range;
  cur_param_idx = 1;
  cur_arg_ptr = argvc[1];
  while (1) {
    // second pass; this can't fail since we already validated
    parse_next_range(argvc, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      if (require_posint) {
        last_val = 0;
        for (cur_param_idx = 0; cur_param_idx < name_ct; ++cur_param_idx) {
          cur_name_str = &(range_list_ptr->names[cur_param_idx * ((uintptr_t)name_max_blen)]);
          const char* dup_check = cur_name_str; // actually a numeric check
          do {
            if (is_not_digit(*dup_check)) {
              LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argvc[0], cur_name_str);
              return kPglRetInvalidCmdline;
            }
          } while (*(++dup_check));
          if (scan_posint_defcap(cur_name_str, &cur_val)) {
            LOGERRPRINTFWW("Error: Invalid %s parameter '%s'.\n", argvc[0], cur_name_str);
            return kPglRetInvalidCmdline;
          }
          if (range_list_ptr->starts_range[cur_param_idx]) {
            last_val = cur_val;
          } else {
            if (cur_val <= last_val) {
              LOGERRPRINTFWW("Error: Invalid %s range '%s-%s'.\n", argvc[0], &(range_list_ptr->names[(cur_param_idx - 1) * name_max_blen]), cur_name_str);
              return kPglRetInvalidCmdline;
            }
            last_val = 0;
          }
        }
      }
      return kPglRetSuccess;
    }
    memcpyx(cur_name_str, range_start, rs_len, 0);
    const char* dup_check = range_list_ptr->names;
    while (dup_check < cur_name_str) {
      if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
        LOGERRPRINTFWW("Error: Duplicate %s parameter '%s'.\n", argvc[0], cur_name_str);
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
          LOGERRPRINTFWW("Error: Duplicate %s parameter '%s'.\n", argvc[0], cur_name_str);
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


uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c, double* solutions) {
  // Additional research into numerical stability may be in order here...
  double a2 = coef_a * coef_a;
  double qq = (a2 - 3 * coef_b) * (1.0 / 9.0);
  double rr = (2 * a2 * coef_a - 9 * coef_a * coef_b + 27 * coef_c) * (1.0 / 54.0);
  double r2 = rr * rr;
  double q3 = qq * qq * qq;
  double adiv3 = coef_a * (1.0 / 3.0);
  double sq;
  double dxx;
  if (r2 < q3) {
    // three real roots
    sq = sqrt(qq);
    dxx = acos(rr / (qq * sq)) * (1.0 / 3.0);
    sq *= -2;
    solutions[0] = sq * cos(dxx) - adiv3;
    solutions[1] = sq * cos(dxx + (2.0 * kPi / 3.0)) - adiv3;
    solutions[2] = sq * cos(dxx - (2.0 * kPi / 3.0)) - adiv3;
    // now sort and check for within-epsilon equality
    if (solutions[0] > solutions[1]) {
      dxx = solutions[0];
      solutions[0] = solutions[1];
      if (dxx > solutions[2]) {
        solutions[1] = solutions[2];
        solutions[2] = dxx;
      } else {
        solutions[1] = dxx;
      }
      if (solutions[0] > solutions[1]) {
        dxx = solutions[0];
        solutions[0] = solutions[1];
        solutions[1] = dxx;
      }
    } else if (solutions[1] > solutions[2]) {
      dxx = solutions[1];
      solutions[1] = solutions[2];
      solutions[2] = dxx;
    }
    if (solutions[1] - solutions[0] < kEpsilon) {
      solutions[1] = solutions[2];
      return 2 - (solutions[1] - solutions[0] < kEpsilon);
    }
    return 3 - (solutions[2] - solutions[1] < kEpsilon);
  }
  dxx = -pow(fabs(rr) + sqrt(r2 - q3), 1.0 / 3.0);
  if (dxx == 0.0) {
    solutions[0] = -adiv3;
    return 1;
  }
  if (rr < 0.0) {
    dxx = -dxx;
  }
  sq = qq / dxx;
  solutions[0] = dxx + sq - adiv3;
  // use of regular epsilon here has actually burned us
  if (fabs(dxx - sq) >= (kEpsilon * 8)) {
    return 1;
  }
  if (dxx >= 0.0) {
    solutions[1] = solutions[0];
    solutions[0] = -dxx - adiv3;
  } else {
    solutions[1] = -dxx - adiv3;
  }
  return 2;
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

#ifndef _WIN32
pthread_attr_t g_smallstack_thread_attr;
#endif

boolerr_t spawn_threads(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, pthread_t* threads) {
  uintptr_t ulii;
  if (ct == 1) {
    return 0;
  }
  for (ulii = 1; ulii < ct; ++ulii) {
#ifdef _WIN32
    threads[ulii - 1] = (HANDLE)_beginthreadex(nullptr, kDefaultThreadStack, start_routine, (void*)ulii, 0, nullptr);
    if (!threads[ulii - 1]) {
      join_threads(ulii, threads);
      return 1;
    }
#else
    if (pthread_create(&(threads[ulii - 1]), &g_smallstack_thread_attr, start_routine, (void*)ulii)) {
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
            for (uintptr_t uljj = 0; uljj < ulii; ++uljj) {
              TerminateThread(threads[uljj], 0);
            }
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
      if (pthread_create(&(threads[ulii]), &g_smallstack_thread_attr, start_routine, (void*)ulii)) {
        if (ulii) {
          if (is_last_block) {
            join_threads2z(ulii, 1, threads);
          } else {
            const uintptr_t unstarted_thread_ct = ct - ulii;
            pthread_mutex_lock(&g_thread_sync_mutex);
            // minor bugfix (21 Aug 2017): join_threads2z hangs forever if not
            //   last block, since cur_block_done_condvar is only signaled
            //   when g_thread_active_ct decreased to zero
            g_thread_active_ct -= unstarted_thread_ct;
            while (g_thread_active_ct) {
              pthread_cond_wait(&g_thread_cur_block_done_condvar, &g_thread_sync_mutex);
            }
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
static const char* const* g_item_ids;
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
  const char* const* item_ids = g_item_ids;
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

pglerr_t populate_id_htable_mt(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t store_all_dups, uint32_t id_htable_size, uint32_t thread_ct, uint32_t* id_htable) {
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
        if (cur_htable_entry == UINT32_MAX) {
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
            if (cur_htable_entry == UINT32_MAX) {
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
        if (cur_htable_entry == UINT32_MAX) {
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
                htable_dup_base[extra_alloc + 1] = UINT32_MAX; // list end
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
            if (cur_htable_entry == UINT32_MAX) {
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

pglerr_t alloc_and_populate_id_htable_mt(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t max_thread_ct, uint32_t** id_htable_ptr, uint32_t** htable_dup_base_ptr, uint32_t* id_htable_size_ptr) {
  uint32_t id_htable_size = get_htable_fast_size(item_ct);
  // 4 bytes per variant for hash buffer
  // if store_all_dups, up to 8 bytes per variant in extra_alloc for duplicate
  //   tracking
  const uint32_t store_all_dups = (htable_dup_base_ptr != nullptr);
  const uintptr_t nonhtable_alloc = round_up_pow2(item_ct * sizeof(int32_t), kCacheline) + store_all_dups * round_up_pow2(item_ct * 2 * sizeof(int32_t), kCacheline);
  uintptr_t max_bytes = round_down_pow2(bigstack_left(), kCacheline);
  // force max_bytes >= 5 so leqprime() doesn't fail
  // (not actually relevant any more, but whatever)
  if (nonhtable_alloc + (item_ct + 6) * sizeof(int32_t) > max_bytes) {
    return kPglRetNomem;
  }
  max_bytes -= nonhtable_alloc;
  if (id_htable_size * sizeof(int32_t) > max_bytes) {
    id_htable_size = max_bytes / sizeof(int32_t);
    // id_htable_size = leqprime((max_bytes / sizeof(int32_t)) - 1);
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


uint32_t edit1_match(const char* s1, const char* s2, uint32_t len1, uint32_t len2) {
  // permit one difference of the following forms:
  // - inserted/deleted character
  // - replaced character
  // - adjacent pair of swapped characters
  uint32_t diff_found = 0;
  uint32_t pos = 0;
  if (len1 == len2) {
    while (pos < len1) {
      if (s1[pos] != s2[pos]) {
        if (diff_found) {
          if ((diff_found == 2) || (s1[pos] != s2[pos - 1]) || (s1[pos - 1] != s2[pos])) {
            return 0;
          }
        }
        ++diff_found;
      }
      ++pos;
    }
  } else if (len1 == len2 - 1) {
    do {
      if (s1[pos - diff_found] != s2[pos]) {
        if (diff_found) {
          return 0;
        }
        ++diff_found;
      }
      ++pos;
    } while (pos < len2);
  } else if (len1 == len2 + 1) {
    do {
      if (s1[pos] != s2[pos - diff_found]) {
        if (diff_found) {
          return 0;
        }
        ++diff_found;
      }
      ++pos;
    } while (pos < len1);
  } else {
    return 0;
  }
  return 1;
}

CONSTU31(kMaxEqualHelpParams, 64);

void help_print(const char* cur_params, help_ctrl_t* help_ctrl_ptr, uint32_t postprint_newline, const char* payload) {
  if (help_ctrl_ptr->param_ct) {
    strcpy(g_textbuf, cur_params);
    uint32_t cur_param_ct = 1;
    char* cur_param_start[kMaxEqualHelpParams];
    cur_param_start[0] = g_textbuf;
    char* textbuf_iter = strchr(g_textbuf, '\t');
    while (textbuf_iter) {
      *textbuf_iter++ = '\0';
      cur_param_start[cur_param_ct++] = textbuf_iter;
      textbuf_iter = strchr(textbuf_iter, '\t');
    }
    if (help_ctrl_ptr->iters_left) {
      const uint32_t orig_unmatched_ct = help_ctrl_ptr->unmatched_ct;
      if (help_ctrl_ptr->unmatched_ct) {
        uint32_t arg_uidx = 0;
        if (help_ctrl_ptr->iters_left == 2) {
          for (uint32_t arg_idx = 0; arg_idx < orig_unmatched_ct; ++arg_idx, ++arg_uidx) {
            arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
            for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
              if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
                SET_BIT(arg_uidx, help_ctrl_ptr->perfect_match_arr);
                SET_BIT(arg_uidx, help_ctrl_ptr->prefix_match_arr);
                SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
                help_ctrl_ptr->unmatched_ct -= 1;
                break;
              }
            }
          }
        } else {
          uint32_t cur_param_slens[kMaxEqualHelpParams];
          for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
            cur_param_slens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
          }
          for (uint32_t arg_idx = 0; arg_idx < orig_unmatched_ct; ++arg_idx, ++arg_uidx) {
            arg_uidx = next_unset_unsafe(help_ctrl_ptr->all_match_arr, arg_uidx);
            const uint32_t slen = help_ctrl_ptr->param_slens[arg_uidx];
            for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
              if (cur_param_slens[cur_param_idx] > slen) {
                if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], slen)) {
                  SET_BIT(arg_uidx, help_ctrl_ptr->prefix_match_arr);
                  SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
                  help_ctrl_ptr->unmatched_ct -= 1;
                  break;
                }
              }
            }
          }
        }
      }
    } else {
      uint32_t cur_param_slens[kMaxEqualHelpParams];
      for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
        cur_param_slens[cur_param_idx] = strlen(cur_param_start[cur_param_idx]);
      }
      uint32_t print_this = 0;
      for (uint32_t arg_uidx = 0; arg_uidx < help_ctrl_ptr->param_ct; ++arg_uidx) {
        if (IS_SET(help_ctrl_ptr->prefix_match_arr, arg_uidx)) {
          if (!print_this) {
            if (IS_SET(help_ctrl_ptr->perfect_match_arr, arg_uidx)) {
              for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
                if (!strcmp(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx])) {
                  print_this = 1;
                  break;
                }
              }
            } else {
              const uint32_t slen = help_ctrl_ptr->param_slens[arg_uidx];
              for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
                if (cur_param_slens[cur_param_idx] > slen) {
                  if (!memcmp(help_ctrl_ptr->argv[arg_uidx], cur_param_start[cur_param_idx], slen)) {
                    print_this = 1;
                    break;
                  }
                }
              }
            }
          }
        } else {
          for (uint32_t cur_param_idx = 0; cur_param_idx < cur_param_ct; ++cur_param_idx) {
            if (edit1_match(cur_param_start[cur_param_idx], help_ctrl_ptr->argv[arg_uidx], cur_param_slens[cur_param_idx], help_ctrl_ptr->param_slens[arg_uidx])) {
              print_this = 1;
              if (!IS_SET(help_ctrl_ptr->all_match_arr, arg_uidx)) {
                SET_BIT(arg_uidx, help_ctrl_ptr->all_match_arr);
                help_ctrl_ptr->unmatched_ct -= 1;
              }
              break;
            }
          }
        }
      }
      if (print_this) {
        const uint32_t payload_slen = strlen(payload);
        const char* payload_end;
        if (payload[payload_slen - 2] == '\n') {
          payload_end = &(payload[payload_slen - 1]);
        } else {
          payload_end = &(payload[payload_slen]);
        }
        if (help_ctrl_ptr->preprint_newline) {
          putc_unlocked('\n', stdout);
        }
        help_ctrl_ptr->preprint_newline = postprint_newline;
        const char* payload_iter = payload;
        do {
          const char* line_end = (const char*)rawmemchr(payload_iter, '\n') + 1;
          uint32_t line_slen = (uint32_t)(line_end - payload_iter);
          if (line_slen > 2) {
            payload_iter = &(payload_iter[2]);
            line_slen -= 2;
          }
          memcpyx(g_textbuf, payload_iter, line_slen, 0);
          fputs(g_textbuf, stdout);
          payload_iter = line_end;
        } while (payload_iter < payload_end);
      }
    }
  } else {
    fputs(payload, stdout);
  }
}


void plink2_cmdline_meta_preinit(plink2_cmdline_meta_t* pcmp) {
  pcmp->subst_argv = nullptr;
  pcmp->script_buf = nullptr;
  pcmp->rerun_buf = nullptr;
  pcmp->flag_buf = nullptr;
  pcmp->flag_map = nullptr;
}

const char errstr_nomem[] = "Error: Out of memory.  The --memory flag may be helpful.\n";
const char errstr_write[] = "Error: File write failure.\n";
const char errstr_read[] = "Error: File read failure.\n";
const char errstr_thread_create[] = "Error: Failed to create thread.\n";

// assumes logfile is open
void disp_exit_msg(pglerr_t reterr) {
  if (reterr) {
    if (reterr == kPglRetNomem) {
      logprint("\n");
      logerrprint(errstr_nomem);
      if (g_failed_alloc_attempt_size) {
        LOGERRPRINTF("Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
      }
    } else if (reterr == kPglRetReadFail) {
      logprint("\n");
      logerrprint(errstr_read);
    } else if (reterr == kPglRetWriteFail) {
      logprint("\n");
      logerrprint(errstr_write);
    } else if (reterr == kPglRetThreadCreateFail) {
      logprint("\n");
      logerrprint(errstr_thread_create);
    }
  }
}

// useful when there's e.g. a filename and an optional modifier, and we want to
// permit either parmeter ordering
boolerr_t check_extra_param(const char* const* argvc, const char* permitted_modif, uint32_t* other_idx_ptr) {
  const uint32_t idx_base = *other_idx_ptr;
  if (!strcmp(argvc[idx_base], permitted_modif)) {
    *other_idx_ptr = idx_base + 1;
  } else if (strcmp(argvc[idx_base + 1], permitted_modif)) {
    LOGERRPRINTF("Error: Invalid %s parameter sequence.\n", argvc[0]);
    return 1;
  }
  return 0;
}

char extract_char_param(const char* ss) {
  // maps c, 'c', and "c" to c, and anything else to the null char.  This is
  // intended to support e.g. always using '#' to designate a # parameter
  // without worrying about differences between shells.
  const char cc = ss[0];
  if (((cc == '\'') || (cc == '"')) && (ss[1]) && (ss[2] == cc) && (!ss[3])) {
    return ss[1];
  }
  if (cc && (!ss[1])) {
    return cc;
  }
  return '\0';
}

pglerr_t cmdline_alloc_string(const char* source, const char* flag_name, uint32_t max_slen, char** sbuf_ptr) {
  const uint32_t slen = strlen(source);
  if (slen > max_slen) {
    LOGERRPRINTF("Error: %s parameter too long.\n", flag_name);
    return kPglRetInvalidCmdline;
  }
  const uint32_t blen = slen + 1;
  if (pgl_malloc(blen, sbuf_ptr)) {
    return kPglRetNomem;
  }
  memcpy(*sbuf_ptr, source, blen);
  return kPglRetSuccess;
}

pglerr_t alloc_fname(const char* source, const char* flagname_p, uint32_t extra_size, char** fnbuf_ptr) {
  const uint32_t blen = strlen(source) + 1;
  if (blen > (kPglFnamesize - extra_size)) {
    LOGERRPRINTF("Error: --%s filename too long.\n", flagname_p);
    return kPglRetOpenFail;
  }
  if (pgl_malloc(blen + extra_size, fnbuf_ptr)) {
    return kPglRetNomem;
  }
  memcpy(*fnbuf_ptr, source, blen);
  return kPglRetSuccess;
}

pglerr_t alloc_and_flatten(const char* const* sources, uint32_t param_ct, uint32_t max_blen, char** flattened_buf_ptr) {
  uintptr_t tot_blen = 1;
  for (uint32_t param_idx = 0; param_idx < param_ct; ++param_idx) {
    const uint32_t cur_blen = 1 + strlen(sources[param_idx]);
    if (cur_blen > max_blen) {
      return kPglRetInvalidCmdline;
    }
    tot_blen += cur_blen;
  }
  char* buf_iter;
  if (pgl_malloc(tot_blen, &buf_iter)) {
    return kPglRetNomem;
  }
  *flattened_buf_ptr = buf_iter;
  for (uint32_t param_idx = 0; param_idx < param_ct; ++param_idx) {
    buf_iter = strcpyax(buf_iter, sources[param_idx], '\0');
  }
  *buf_iter = '\0';
  return kPglRetSuccess;
}

pglerr_t rerun(const char* ver_str, const char* ver_str2, const char* prog_name_str, uint32_t rerun_argv_pos, uint32_t rerun_parameter_present, int32_t* argc_ptr, uint32_t* first_arg_idx_ptr, char*** argv_ptr, char*** subst_argv_ptr, char** rerun_buf_ptr) {
  // caller is responsible for freeing rerun_buf

  // ok, requiring zlib/zstd here is totally not worth it
  FILE* rerunfile = nullptr;

  char** subst_argv2 = nullptr;
  char** argv = *argv_ptr;
  char* rerun_fname = rerun_parameter_present? argv[rerun_argv_pos + 1] : g_textbuf;
  uintptr_t line_idx = 1;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!rerun_parameter_present) {
      char* write_iter = strcpya(rerun_fname, prog_name_str);
      strcpy(write_iter, ".log");
    }
    rerunfile = fopen(rerun_fname, FOPEN_RB);
    if (!rerunfile) {
      goto rerun_ret_OPEN_FAIL;
    }
    char* textbuf = g_textbuf;
    textbuf[kMaxMediumLine - 1] = ' ';
    if (!fgets(textbuf, kMaxMediumLine, rerunfile)) {
      fputs(ver_str, stdout);
      fputs(ver_str2, stdout);
      fputs("Error: Empty log file for --rerun.\n", stderr);
      goto rerun_ret_MALFORMED_INPUT;
    }
    if (!textbuf[kMaxMediumLine - 1]) {
      goto rerun_ret_LONG_LINE;
    }
    if (!fgets(textbuf, kMaxMediumLine, rerunfile)) {
      fputs(ver_str, stdout);
      fputs(ver_str2, stdout);
      fputs("Error: Only one line in --rerun log file.\n", stderr);
      goto rerun_ret_MALFORMED_INPUT;
    }
    ++line_idx;
    if (!textbuf[kMaxMediumLine - 1]) {
      goto rerun_ret_LONG_LINE;
    }
    // don't bother supporting "xx arguments: --aa bb --cc --dd" format
    while ((!str_startswith2(textbuf, "Options in effect:")) || (textbuf[strlen("Options in effect:")] >= ' ')) {
      ++line_idx;
      if (!fgets(textbuf, kMaxMediumLine, rerunfile)) {
        fputs(ver_str, stdout);
        fputs(ver_str2, stdout);
        fputs("Error: Invalid log file for --rerun.\n", stderr);
        goto rerun_ret_MALFORMED_INPUT;
      }
    }
    char* all_args_write_iter = textbuf;
    char* textbuf_limit = &(textbuf[kMaxMediumLine]);
    uint32_t loaded_arg_ct = 0;
    // We load each of the option lines in sequence into textbuf, always
    // overwriting the previous line's newline.  (Note that textbuf[] has
    // size > 2 * kMaxMediumLine; this lets us avoid additional
    // dynamic memory allocation as long as we impose the constraint that all
    // lines combined add up to less than kMaxMediumLine.)
    while (1) {
      all_args_write_iter[kMaxMediumLine - 1] = ' ';
      if (!fgets(all_args_write_iter, kMaxMediumLine, rerunfile)) {
        break;
      }
      ++line_idx;
      if (!all_args_write_iter[kMaxMediumLine - 1]) {
        goto rerun_ret_LONG_LINE;
      }
      char* arg_iter = skip_initial_spaces(all_args_write_iter);
      if (is_eoln_kns(*arg_iter)) {
        *all_args_write_iter = '\0';
        break;
      }
      char* token_end;
      do {
        token_end = token_endnn(arg_iter);
        ++loaded_arg_ct;
        arg_iter = skip_initial_spaces(token_end);
      } while (!is_eoln_kns(*arg_iter));
      all_args_write_iter = token_end;
      if (all_args_write_iter >= textbuf_limit) {
        fputs(ver_str, stdout);
        fputs(ver_str2, stdout);
        fputs("Error: --rerun argument sequence too long.\n", stderr);
        goto rerun_ret_MALFORMED_INPUT;
      }
    }
    fclose_null(&rerunfile);
    const uint32_t line_byte_ct = 1 + (uintptr_t)(all_args_write_iter - textbuf);
    char* rerun_buf;
    if (pgl_malloc(line_byte_ct, &rerun_buf)) {
      goto rerun_ret_NOMEM;
    }
    *rerun_buf_ptr = rerun_buf;
    memcpy(rerun_buf, textbuf, line_byte_ct);
    const uint32_t argc = (uint32_t)(*argc_ptr);
    const uint32_t first_arg_idx = *first_arg_idx_ptr;
    char* rerun_first_token = skip_initial_spaces(rerun_buf);
    const char* arg_iter = rerun_first_token;
    // now use textbuf as a lame bitfield
    memset(textbuf, 1, loaded_arg_ct);
    uint32_t loaded_arg_idx = 0;
    uint32_t duplicate_arg_ct = 0;
    do {
      if (no_more_tokens_kns(arg_iter)) {
        fputs(ver_str, stdout);
        fputs(ver_str2, stdout);
        fputs("Error: Line 2 of --rerun log file has fewer tokens than expected.\n", stderr);
        goto rerun_ret_MALFORMED_INPUT;
      }
      const char* flagname_p = is_flag_start(arg_iter);
      if (flagname_p) {
        const uint32_t slen = strlen_se(flagname_p);
        uint32_t cmdline_arg_idx = first_arg_idx;
        for (; cmdline_arg_idx < argc; cmdline_arg_idx++) {
          const char* later_flagname_p = is_flag_start(argv[cmdline_arg_idx]);
          if (later_flagname_p) {
            const uint32_t slen2 = strlen(later_flagname_p);
            if ((slen == slen2) && (!memcmp(flagname_p, later_flagname_p, slen))) {
              cmdline_arg_idx = UINT32_MAX;
              break;
            }
          }
        }
        if (cmdline_arg_idx == UINT32_MAX) {
          // matching flag, override --rerun
          do {
            ++duplicate_arg_ct;
            textbuf[loaded_arg_idx++] = 0;
            if (loaded_arg_idx == loaded_arg_ct) {
              break;
            }
            arg_iter = next_token(arg_iter);
          } while (!is_flag(arg_iter));
        } else {
          ++loaded_arg_idx;
          arg_iter = next_token(arg_iter);
        }
      } else {
        ++loaded_arg_idx;
        arg_iter = next_token(arg_iter);
      }
    } while (loaded_arg_idx < loaded_arg_ct);
    if (pgl_malloc((argc + loaded_arg_ct - duplicate_arg_ct - rerun_parameter_present - 1 - first_arg_idx) * sizeof(intptr_t), &subst_argv2)) {
      goto rerun_ret_NOMEM;
    }
    uint32_t new_arg_idx = rerun_argv_pos - first_arg_idx;
    memcpy(subst_argv2, &(argv[first_arg_idx]), new_arg_idx * sizeof(intptr_t));
    char* arg_nullterminate_iter = rerun_first_token;
    for (loaded_arg_idx = 0; loaded_arg_idx < loaded_arg_ct; ++loaded_arg_idx) {
      arg_nullterminate_iter = skip_initial_spaces(arg_nullterminate_iter);
      char* token_end = token_endnn(arg_nullterminate_iter);
      if (textbuf[loaded_arg_idx]) {
        subst_argv2[new_arg_idx++] = arg_nullterminate_iter;
        *token_end = '\0';
      }
      arg_nullterminate_iter = &(token_end[1]);
    }
    const uint32_t final_copy_start_idx = rerun_argv_pos + rerun_parameter_present + 1;
    memcpy(&(subst_argv2[new_arg_idx]), &(argv[final_copy_start_idx]), (argc - final_copy_start_idx) * sizeof(intptr_t));
    *first_arg_idx_ptr = 0;
    *argc_ptr = new_arg_idx + argc - final_copy_start_idx;
    if (*subst_argv_ptr) {
      free(*subst_argv_ptr);
    }
    *subst_argv_ptr = subst_argv2;
    *argv_ptr = subst_argv2;
    subst_argv2 = nullptr;
  }
  while (0) {
  rerun_ret_NOMEM:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    reterr = kPglRetNomem;
    break;
  rerun_ret_OPEN_FAIL:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    fprintf(stderr, g_errstr_fopen, rerun_fname);
    reterr = kPglRetOpenFail;
    break;
  rerun_ret_LONG_LINE:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    fprintf(stderr, "Error: Line %" PRIuPTR " of --rerun log file is pathologically long.\n", line_idx);
  rerun_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  free_cond(subst_argv2);
  fclose_cond(rerunfile);
  return reterr;
}


// Handles --script, --rerun, --help, --version, and --silent.
// subst_argv, script_buf, and rerun_buf must be initialized to nullptr.
pglerr_t cmdline_parse_phase1(const char* ver_str, const char* ver_str2, const char* prog_name_str, const char* notestr_null_calc2, const char* cmdline_format_str, const char* errstr_append, uint32_t max_flag_blen, pglerr_t(* disp_help_fn)(uint32_t, const char* const*), int* argc_ptr, char*** argv_ptr, plink2_cmdline_meta_t* pcmp, uint32_t* first_arg_idx_ptr, uint32_t* flag_ct_ptr) {
  FILE* scriptfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    int argc = *argc_ptr;
    const char* const* argvc = TO_CONSTCPCONSTP(*argv_ptr);
    char** subst_argv = nullptr;
    uint32_t first_arg_idx = 1;
    for (uint32_t arg_idx = 1; arg_idx < (uint32_t)argc; ++arg_idx) {
      if ((!strcmp("-script", argvc[arg_idx])) || (!strcmp("--script", argvc[arg_idx]))) {
        const uint32_t param_ct = param_count(argvc, argc, arg_idx);
        if (enforce_param_ct_range(argvc[arg_idx], param_ct, 1, 1)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
        }
        for (uint32_t arg_idx2 = arg_idx + 2; arg_idx2 < (uint32_t)argc; ++arg_idx2) {
          if ((!strcmp("-script", argvc[arg_idx2])) || (!strcmp("--script", argvc[arg_idx2]))) {
            fputs(ver_str, stdout);
            fputs(ver_str2, stdout);
            fputs("Error: Multiple --script flags.  Merge the files into one.\n", stderr);
            fputs(errstr_append, stderr);
            goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
          }
        }
        // logging not yet active, so don't use fopen_checked()
        scriptfile = fopen(argvc[arg_idx + 1], FOPEN_RB);
        if (!scriptfile) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fprintf(stderr, g_errstr_fopen, argvc[arg_idx + 1]);
          goto cmdline_parse_phase1_ret_OPEN_FAIL;
        }
        if (fseeko(scriptfile, 0, SEEK_END)) {
          goto cmdline_parse_phase1_ret_READ_FAIL;
        }
        int64_t fsize = ftello(scriptfile);
        if (fsize < 0) {
          goto cmdline_parse_phase1_ret_READ_FAIL;
        }
        if (fsize > 0x7ffffffe) {
          // could actually happen if user enters parameters in the wrong
          // order, so may as well catch it and print a somewhat informative
          // error message
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs("Error: --script file too large.", stderr);
          goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
        }
        rewind(scriptfile);
        const uint32_t fsize_ui = (uint64_t)fsize;
        if (pgl_malloc(fsize_ui + 1, &pcmp->script_buf)) {
          goto cmdline_parse_phase1_ret_NOMEM;
        }
        char* script_buf = pcmp->script_buf;
        if (!fread_unlocked(script_buf, fsize_ui, 1, scriptfile)) {
          goto cmdline_parse_phase1_ret_READ_FAIL;
        }
        script_buf[fsize_ui] = '\0';
        fclose_null(&scriptfile);
        uint32_t num_script_params = 0;
        char* script_buf_iter = script_buf;
        uint32_t char_code;
        do {
          uint32_t char_code_m1;
          do {
            char_code_m1 = ((uint32_t)((unsigned char)(*script_buf_iter++))) - 1;
          } while (char_code_m1 < 32);
          if (char_code_m1 == UINT32_MAX) {
            break;
          }
          ++num_script_params;
          do {
            char_code = (uint32_t)((unsigned char)(*script_buf_iter++));
          } while (char_code > 32);
        } while (char_code);
        if (script_buf_iter != (&(script_buf[fsize_ui + 1]))) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs("Error: Null byte in --script file.\n", stderr);
          goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
        }
        // probable todo: detect duplicate flags in the same manner as --rerun
        const uint32_t new_param_ct = num_script_params + argc - 3;
        if (pgl_malloc(new_param_ct * sizeof(intptr_t), &subst_argv)) {
          goto cmdline_parse_phase1_ret_NOMEM;
        }
        pcmp->subst_argv = subst_argv;
        memcpy(subst_argv, &(argvc[1]), arg_idx * sizeof(intptr_t));
        const uint32_t load_param_idx_end = arg_idx + num_script_params;
        script_buf_iter = &(script_buf[-1]);
        for (uint32_t param_idx = arg_idx; param_idx < load_param_idx_end; ++param_idx) {
          while (((unsigned char)(*(++script_buf_iter))) <= 32);
          subst_argv[param_idx] = script_buf_iter;
          while (((unsigned char)(*(++script_buf_iter))) > 32);
          // could enforce some sort of length limit here
          *script_buf_iter = '\0';
        }
        memcpy(&(subst_argv[load_param_idx_end]), &(argvc[arg_idx + 2]), (argc - arg_idx - 2) * sizeof(intptr_t));
        argc = new_param_ct;
        *argc_ptr = argc;
        first_arg_idx = 0;
        argvc = TO_CONSTCPCONSTP(subst_argv);
        *argv_ptr = subst_argv;
        break;
      }
    }
    for (uint32_t arg_idx = first_arg_idx; arg_idx < (uint32_t)argc; ++arg_idx) {
      if ((!strcmp("-rerun", argvc[arg_idx])) || (!strcmp("--rerun", argvc[arg_idx]))) {
        const uint32_t param_ct = param_count(argvc, argc, arg_idx);
        if (enforce_param_ct_range(argvc[arg_idx], param_ct, 0, 1)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
        }
        for (uint32_t arg_idx2 = arg_idx + param_ct + 1; arg_idx2 < (uint32_t)argc; ++arg_idx2) {
          if ((!strcmp("-rerun", argvc[arg_idx2])) || (!strcmp("--rerun", argvc[arg_idx2]))) {
            fputs(ver_str, stdout);
            fputs(ver_str2, stdout);
            fputs("Error: Duplicate --rerun flag.\n", stderr);
            goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
          }
        }
        reterr = rerun(ver_str, ver_str2, prog_name_str, arg_idx, param_ct, &argc, &first_arg_idx, argv_ptr, &pcmp->subst_argv, &pcmp->rerun_buf);
        if (reterr) {
          goto cmdline_parse_phase1_ret_1;
        }
        *argc_ptr = argc;
        subst_argv = pcmp->subst_argv;
        argvc = TO_CONSTCPCONSTP(*argv_ptr);
        break;
      }
    }
    if ((first_arg_idx < (uint32_t)argc) && (!is_flag(argvc[first_arg_idx]))) {
      fputs("Error: First parameter must be a flag.\n", stderr);
      fputs(errstr_append, stderr);
      goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
    }
    uint32_t flag_ct = 0;
    uint32_t version_present = 0;
    uint32_t silent_present = 0;
    for (uint32_t arg_idx = first_arg_idx; arg_idx < (uint32_t)argc; ++arg_idx) {
      const char* flagname_p = is_flag_start(argvc[arg_idx]);
      if (flagname_p) {
        if (!strcmp("help", flagname_p)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          if ((!first_arg_idx) || (arg_idx != 1) || subst_argv) {
            fputs("--help present, ignoring other flags.\n", stdout);
          }
          if ((arg_idx == ((uint32_t)argc) - 1) && flag_ct) {
            // make "plink [valid flags/parameters] --help" work, and skip the
            // parameters
            const char** help_argv;
            if (pgl_malloc(flag_ct * sizeof(intptr_t), &help_argv)) {
              goto cmdline_parse_phase1_ret_NOMEM2;
            }
            uint32_t arg_idx2 = 0;
            for (uint32_t flag_idx = 0; flag_idx < flag_ct; ++flag_idx) {
              while (!is_flag_start(argvc[++arg_idx2]));
              help_argv[flag_idx] = argvc[arg_idx2];
            }
            reterr = disp_help_fn(flag_ct, help_argv);
            free(help_argv);
          } else {
            reterr = disp_help_fn(argc - arg_idx - 1, &(argvc[arg_idx + 1]));
          }
          if (!reterr) {
            reterr = kPglRetHelp;
          }
          goto cmdline_parse_phase1_ret_1;
        }
        if ((!strcmp("h", flagname_p)) || (!strcmp("?", flagname_p))) {
          // these just act like the no-parameter case
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          if ((!first_arg_idx) || (arg_idx != 1) || subst_argv) {
            printf("-%c present, ignoring other flags.\n", *flagname_p);
          }
          fputs(cmdline_format_str, stdout);
          fputs(notestr_null_calc2, stdout);
          reterr = kPglRetHelp;
          goto cmdline_parse_phase1_ret_1;
        }
        if (!strcmp("version", flagname_p)) {
          version_present = 1;
        } else if (!strcmp("silent", flagname_p)) {
          silent_present = 1;
        }
        if (strlen(flagname_p) >= max_flag_blen) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          // shouldn't be possible for this to overflow the buffer...
          snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized flag ('%s').\n", argvc[arg_idx]);
          wordwrapb(0);
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto cmdline_parse_phase1_ret_INVALID_CMDLINE;
        }
        ++flag_ct;
      }
    }
    if (version_present) {
      fputs(ver_str, stdout);
      putc_unlocked('\n', stdout);
      reterr = kPglRetHelp;
      goto cmdline_parse_phase1_ret_1;
    }
    if (silent_present) {
      if (!freopen("/dev/null", "w", stdout)) {
        fputs("Warning: --silent failed.", stderr);
        g_stderr_written_to = 1;
      }
    }
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    *first_arg_idx_ptr = first_arg_idx;
    *flag_ct_ptr = flag_ct;
  }
  while (0) {
  cmdline_parse_phase1_ret_NOMEM:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
  cmdline_parse_phase1_ret_NOMEM2:
    fputs(errstr_nomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  cmdline_parse_phase1_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  cmdline_parse_phase1_ret_READ_FAIL:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    fputs(errstr_read, stderr);
    reterr = kPglRetReadFail;
    break;
  cmdline_parse_phase1_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 cmdline_parse_phase1_ret_1:
  fclose_cond(scriptfile);
  return reterr;
}

// Assumes cmdline_parse_phase1() has completed, flag names have been copied to
// flag_buf/flag_map, aliases handled, and PROG_NAME_STR has been copied to
// outname (no null-terminator needed).  outname_end must be initialized to
// nullptr.
// This sorts the flag names so they're processed in a predictable order,
// handles --out if present, initializes the log, and determines the number of
// processors the OS wants us to think the machine has.
pglerr_t cmdline_parse_phase2(const char* ver_str, const char* errstr_append, const char* const* argvc, uint32_t prog_name_str_slen, uint32_t max_flag_blen, int32_t argc, uint32_t flag_ct, plink2_cmdline_meta_t* pcmp, char* outname, char** outname_end_ptr, int32_t* known_procs_ptr, uint32_t* max_thread_ct_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    char* flag_buf = pcmp->flag_buf;
    uint32_t* flag_map = pcmp->flag_map;
    reterr = sort_cmdline_flags(max_flag_blen, flag_ct, flag_buf, flag_map);
    if (reterr) {
      if (reterr == kPglRetNomem) {
        goto cmdline_parse_phase2_ret_NOMEM_NOLOG;
      }
      goto cmdline_parse_phase2_ret_1;
    }

    for (uint32_t cur_flag_idx = 0; cur_flag_idx < flag_ct; ++cur_flag_idx) {
      const int32_t memcmp_out_result = memcmp("out", &(flag_buf[cur_flag_idx * max_flag_blen]), 4);
      if (!memcmp_out_result) {
        const uint32_t arg_idx = flag_map[cur_flag_idx];
        const uint32_t param_ct = param_count(argvc, argc, arg_idx);
        if (enforce_param_ct_range(argvc[arg_idx], param_ct, 1, 1)) {
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto cmdline_parse_phase2_ret_INVALID_CMDLINE;
        }
        if (strlen(argvc[arg_idx + 1]) > (kPglFnamesize - kMaxOutfnameExtBlen)) {
          fflush(stdout);
          fputs("Error: --out parameter too long.\n", stderr);
          goto cmdline_parse_phase2_ret_OPEN_FAIL;
        }
        const uint32_t slen = strlen(argvc[arg_idx + 1]);
        memcpy(outname, argvc[arg_idx + 1], slen + 1);
        *outname_end_ptr = &(outname[slen]);
      }
      if (memcmp_out_result <= 0) {
        break;
      }
    }
    if (init_logfile(0, outname, (*outname_end_ptr)? (*outname_end_ptr) : &(outname[prog_name_str_slen]))) {
      goto cmdline_parse_phase2_ret_OPEN_FAIL;
    }
    logstr(ver_str);
    logstr("\n");
    logprint("Options in effect:\n");
    for (uint32_t cur_flag_idx = 0; cur_flag_idx < flag_ct; ++cur_flag_idx) {
      logprint("  --");
      logprint(&(flag_buf[cur_flag_idx * max_flag_blen]));
      uint32_t arg_idx = flag_map[cur_flag_idx] + 1;
      while ((arg_idx < (uint32_t)argc) && (!is_flag(argvc[arg_idx]))) {
        logprint(" ");
        logprint(argvc[arg_idx++]);
      }
      logprint("\n");
    }
    logprint("\n");

#ifdef _WIN32
    DWORD windows_dw = kTextbufSize;
    if (GetComputerName(g_textbuf, &windows_dw))
#else
    if (gethostname(g_textbuf, kTextbufSize) != -1)
#endif
    {
      logstr("Hostname: ");
      logstr(g_textbuf);
    }
    logstr("\nWorking directory: ");
    if (!getcwd(g_textbuf, kPglFnamesize)) {
      goto cmdline_parse_phase2_ret_READ_FAIL;
    }
    logstr(g_textbuf);
    logstr("\n");
    logprint("Start time: ");
    time_t rawtime;
    time(&rawtime);
    logprint(ctime(&rawtime));
    // ctime string always has a newline at the end
    logstr("\n");

#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    *max_thread_ct_ptr = sysinfo.dwNumberOfProcessors;
    *known_procs_ptr = *max_thread_ct_ptr;
#else
    *known_procs_ptr = sysconf(_SC_NPROCESSORS_ONLN);
    *max_thread_ct_ptr = ((*known_procs_ptr) == -1)? 1 : ((uint32_t)(*known_procs_ptr));
#endif
    // don't subtract 1 any more since, when max_thread_ct > 2, one of the
    // (virtual) cores will be dedicated to I/O and have lots of downtime.
    //
    // may make kMaxThreads a parameter later.
    if (*max_thread_ct_ptr > kMaxThreads) {
      *max_thread_ct_ptr = kMaxThreads;
    }
  }
  while (0) {
  cmdline_parse_phase2_ret_NOMEM_NOLOG:
    fputs(errstr_nomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  cmdline_parse_phase2_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  cmdline_parse_phase2_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    disp_exit_msg(reterr);
    break;
  cmdline_parse_phase2_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 cmdline_parse_phase2_ret_1:
  return reterr;
}

pglerr_t cmdline_parse_phase3(uintptr_t max_default_mb, uintptr_t malloc_size_mb, uint32_t memory_require, plink2_cmdline_meta_t* pcmp, unsigned char** bigstack_ua_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    if (pcmp->subst_argv) {
      free(pcmp->subst_argv);
      pcmp->subst_argv = nullptr;
    }
    if (pcmp->script_buf) {
      free(pcmp->script_buf);
      pcmp->script_buf = nullptr;
    }
    if (pcmp->rerun_buf) {
      free(pcmp->rerun_buf);
      pcmp->rerun_buf = nullptr;
    }
    if (pcmp->flag_buf) {
      free(pcmp->flag_buf);
      pcmp->flag_buf = nullptr;
    }
    if (pcmp->flag_map) {
      free(pcmp->flag_map);
      pcmp->flag_map = nullptr;
    }

    uint64_t total_mb = detect_mb();
    if (!malloc_size_mb) {
      if (!total_mb) {
        malloc_size_mb = max_default_mb? max_default_mb : kBigstackDefaultMb;
      } else if (total_mb < (kBigstackMinMb * 2)) {
        malloc_size_mb = kBigstackMinMb;
      } else {
        malloc_size_mb = total_mb / 2;
        if (max_default_mb && (malloc_size_mb > max_default_mb)) {
          malloc_size_mb = max_default_mb;
        }
      }
    }
    assert(malloc_size_mb >= kBigstackMinMb);
#ifndef __LP64__
    if (malloc_size_mb > kMalloc32bitMbMax) {
      malloc_size_mb = kMalloc32bitMbMax;
    }
#endif
    if (total_mb) {
      snprintf(g_logbuf, kLogbufSize, "%" PRIu64 " MB RAM detected; reserving %" PRIuPTR " MB for main workspace.\n", total_mb, malloc_size_mb);
    } else {
      snprintf(g_logbuf, kLogbufSize, "Failed to determine total system memory.  Attempting to reserve %" PRIuPTR " MB.\n", malloc_size_mb);
    }
    logprintb();
    uintptr_t malloc_mb_final;
    if (init_bigstack(malloc_size_mb, &malloc_mb_final, bigstack_ua_ptr)) {
      goto cmdline_parse_phase3_ret_NOMEM;
    }
    if (malloc_size_mb != malloc_mb_final) {
      if (memory_require) {
        goto cmdline_parse_phase3_ret_NOMEM;
      }
      LOGPRINTF("Allocated %" PRIuPTR " MB successfully, after larger attempt(s) failed.\n", malloc_mb_final);
    }

#ifndef _WIN32
    pthread_attr_init(&g_smallstack_thread_attr);
    pthread_attr_setstacksize(&g_smallstack_thread_attr, kDefaultThreadStack);
#endif
  }
  while (0) {
  cmdline_parse_phase3_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  return reterr;
}

void plink2_cmdline_meta_cleanup(plink2_cmdline_meta_t* pcmp) {
  free_cond(pcmp->subst_argv);
  free_cond(pcmp->script_buf);
  free_cond(pcmp->rerun_buf);
  free_cond(pcmp->flag_buf);
  free_cond(pcmp->flag_map);
}


void init_cmp_expr(cmp_expr_t* cmp_expr_ptr) {
  cmp_expr_ptr->pheno_name = nullptr;
}

void cleanup_cmp_expr(cmp_expr_t* cmp_expr_ptr) {
  free_cond(cmp_expr_ptr->pheno_name);
}

// may want CXXCONST_CP treatment later
const char* parse_next_binary_op(const char* expr_str, uint32_t expr_slen, const char** op_start_ptr, cmp_binary_op_t* binary_op_ptr) {
  // !=, <>: kCmpOperatorNoteq
  // <: kCmpOperatorLe
  // <=: kCmpOperatorLeq
  // =, ==: kCmpOperatorEq
  // >=: kCmpOperatorGeq
  // >: kCmpOperatorGe
  const char* next_eq = (char*)memchr(expr_str, '=', expr_slen);
  const char* next_lt = (char*)memchr(expr_str, '<', expr_slen);
  const char* next_gt = (char*)memchr(expr_str, '>', expr_slen);
  if (!next_eq) {
    if (!next_lt) {
      if (!next_gt) {
        return nullptr;
      }
      *op_start_ptr = next_gt;
      *binary_op_ptr = kCmpOperatorGe;
      return &(next_gt[1]);
    }
    if (next_gt == (&(next_lt[1]))) {
      *op_start_ptr = next_lt;
      *binary_op_ptr = kCmpOperatorNoteq;
      return &(next_lt[2]);
    }
    if ((!next_gt) || (next_gt > next_lt)) {
      *op_start_ptr = next_lt;
      *binary_op_ptr = kCmpOperatorLe;
      return &(next_lt[1]);
    }
    *op_start_ptr = next_gt;
    *binary_op_ptr = kCmpOperatorGe;
    return &(next_gt[1]);
  }
  if ((!next_lt) || (next_lt > next_eq)) {
    if ((!next_gt) || (next_gt > next_eq)) {
      if ((next_eq != expr_str) && (next_eq[-1] == '!')) {
        *op_start_ptr = &(next_eq[-1]);
        *binary_op_ptr = kCmpOperatorNoteq;
        return &(next_eq[1]);
      }
      *op_start_ptr = next_eq;
      *binary_op_ptr = kCmpOperatorEq;
      return (next_eq[1] == '=')? (&(next_eq[2])) : (&(next_eq[1]));
    }
    *op_start_ptr = next_gt;
    if (next_eq == (&(next_gt[1]))) {
      *binary_op_ptr = kCmpOperatorGeq;
      return &(next_gt[2]);
    }
    *binary_op_ptr = kCmpOperatorGe;
    return &(next_gt[1]);
  }
  if (next_gt == (&(next_lt[1]))) {
    *op_start_ptr = next_lt;
    *binary_op_ptr = kCmpOperatorNoteq;
    return &(next_lt[2]);
  }
  if ((!next_gt) || (next_gt > next_lt)) {
    *op_start_ptr = next_lt;
    if (next_eq == (&(next_lt[1]))) {
      *binary_op_ptr = kCmpOperatorLeq;
      return &(next_lt[2]);
    }
    *binary_op_ptr = kCmpOperatorLe;
    return &(next_lt[1]);
  }
  *op_start_ptr = next_gt;
  if (next_eq == (&(next_gt[1]))) {
    *binary_op_ptr = kCmpOperatorGeq;
    return &(next_gt[2]);
  }
  *binary_op_ptr = kCmpOperatorGe;
  return &(next_gt[1]);
}

pglerr_t validate_and_alloc_cmp_expr(const char* const* sources, const char* flag_name, uint32_t param_ct, cmp_expr_t* cmp_expr_ptr) {
  // Currently four use cases:
  //   [pheno/covar name] [operator] [pheno val]: regular comparison
  //   [pheno/covar name]: existence check
  //   [INFO key] [operator] [val]: regular comparison
  //   [INFO key]: existence check
  // Some key/value validation is deferred to load_pvar()/keep_remove_if(),
  // since the requirements are different (e.g. no semicolons in anything
  // INFO-related, categorical phenotypes can be assumed to not start with a
  // valid number).
  // May support or/and, parentheses later, but need to be careful to not slow
  // down load_pvar() too much in the no-INFO-filter case.
  pglerr_t reterr = kPglRetSuccess;
  {
    if ((param_ct != 1) && (param_ct != 3)) {
      goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
    }
    const char* pheno_name_start = sources[0];
    const char* pheno_val_start;
    uint32_t pheno_name_slen;
    uint32_t pheno_val_slen;
    if (param_ct == 3) {
      pheno_name_slen = strlen(pheno_name_start);
      const char* op_str = sources[1];
      uint32_t op_slen = strlen(op_str);
      // ok to have single/double quotes around operator
      if (op_slen > 2) {
        const char cc = op_str[0];
        if (((cc == '\'') || (cc == '"')) && (op_str[op_slen - 1] == cc)) {
          ++op_str;
          op_slen -= 2;
        }
      }
      const char* op_start;
      const char* op_end = parse_next_binary_op(op_str, op_slen, &op_start, &cmp_expr_ptr->binary_op);
      if ((!op_end) || (*op_end) || (op_start != op_str)) {
        goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
      }
      pheno_val_start = sources[2];
      pheno_val_slen = strlen(pheno_val_start);
    } else {
      // permit param_ct == 1 as long as tokens are unambiguous
      uint32_t expr_slen = strlen(pheno_name_start);
      const char* op_start;
      pheno_val_start = parse_next_binary_op(pheno_name_start, expr_slen, &op_start, &cmp_expr_ptr->binary_op);
      if ((!pheno_val_start) || (!(*pheno_val_start)) || (op_start == pheno_name_start)) {
        goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
      }
      pheno_name_slen = (uintptr_t)(op_start - pheno_name_start);
      // quasi-bugfix (13 Dec 2017): allow single argument to contain internal
      // spaces;
      //   --keep-if "PHENO1 > 1"
      // is more intuitive and consistent with usage of other command-line
      // tools than
      //   --keep-if PHENO1 '>' 1
      //
      // To prevent --rerun from breaking, if there's a space after the
      // operator, there must be a space before the operator as well, etc.
      if (*pheno_val_start == ' ') {
        if (pheno_name_start[pheno_name_slen - 1] != ' ') {
          goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
        }
        do {
          ++pheno_val_start;
        } while (*pheno_val_start == ' ');
        do {
          --pheno_name_slen;
          if (!pheno_name_slen) {
            goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
          }
        } while (pheno_name_start[pheno_name_slen - 1] == ' ');
      }
      pheno_val_slen = expr_slen - ((uintptr_t)(pheno_val_start - pheno_name_start));
    }
    if (memchr(pheno_name_start, ' ', pheno_name_slen) || memchr(pheno_val_start, ' ', pheno_val_slen)) {
      goto validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC;
    }
    if ((pheno_name_slen > kMaxIdSlen) || (pheno_val_slen > kMaxIdSlen)) {
      LOGERRPRINTF("Error: ID too long in %s expression.\n", flag_name);
      goto validate_and_alloc_cmp_expr_ret_INVALID_CMDLINE;
    }
    if ((cmp_expr_ptr->binary_op != kCmpOperatorNoteq) && (cmp_expr_ptr->binary_op != kCmpOperatorEq)) {
      double dxx;
      if (!scanadv_double(pheno_val_start, &dxx)) {
        LOGERRPRINTFWW("Error: Invalid %s value '%s' (finite number expected).\n", flag_name, pheno_val_start);
        goto validate_and_alloc_cmp_expr_ret_INVALID_CMDLINE;
      }
    }
    char* new_pheno_name_buf;
    if (pgl_malloc(2 + pheno_name_slen + pheno_val_slen, &new_pheno_name_buf)) {
      goto validate_and_alloc_cmp_expr_ret_NOMEM;
    }
    memcpyx(new_pheno_name_buf, pheno_name_start, pheno_name_slen, '\0');
    // pheno_val_start guaranteed to be null-terminated for now
    memcpy(&(new_pheno_name_buf[pheno_name_slen + 1]), pheno_val_start, pheno_val_slen + 1);
    cmp_expr_ptr->pheno_name = new_pheno_name_buf;
  }
  while (0) {
  validate_and_alloc_cmp_expr_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  validate_and_alloc_cmp_expr_ret_INVALID_EXPR_GENERIC:
    LOGERRPRINTF("Error: Invalid %s expression.\n", flag_name);
  validate_and_alloc_cmp_expr_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return reterr;
}

// See e.g. meta_analysis_open_and_read_header() in plink 1.9.
// Assumes at least one search term.
pglerr_t search_header_line(const char* header_line_iter, const char* const* search_multistrs, const char* flagname_p, uint32_t search_col_ct, uint32_t* found_col_ct_ptr, uint32_t* found_type_bitset_ptr, uint32_t* col_skips, uint32_t* col_types) {
  assert(search_col_ct <= 32);
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    uint32_t search_term_ct = 0;
    uintptr_t max_blen = 0;
    for (uint32_t search_col_idx = 0; search_col_idx < search_col_ct; ++search_col_idx) {
      const uint32_t cur_search_term_ct = count_and_measure_multistr(search_multistrs[search_col_idx], &max_blen);
      assert(cur_search_term_ct <= (1 << 26));
      search_term_ct += cur_search_term_ct;
    }
    char* merged_strbox;
    uint32_t* id_map;
    uint32_t* priority_vals;
    uint64_t* cols_and_types;
    if (bigstack_alloc_c(search_term_ct * max_blen, &merged_strbox) ||
        bigstack_alloc_ui(search_term_ct, &id_map) ||
        bigstack_alloc_ui(search_col_ct, &priority_vals) ||
        bigstack_alloc_ull(search_col_ct, &cols_and_types)) {
      goto search_header_line_ret_NOMEM;
    }
    uint32_t search_term_idx = 0;
    for (uint32_t search_col_idx = 0; search_col_idx < search_col_ct; ++search_col_idx) {
      const char* multistr_iter = search_multistrs[search_col_idx];
      uint32_t priority_idx = 0;  // 0 = highest priority
      while (*multistr_iter) {
        const uint32_t cur_blen = strlen(multistr_iter) + 1;
        memcpy(&(merged_strbox[search_term_idx * max_blen]), multistr_iter, cur_blen);
        id_map[search_term_idx] = search_col_idx + (priority_idx * 32);
        ++search_term_idx;
        ++priority_idx;
        multistr_iter = &(multistr_iter[cur_blen]);
      }
    }
    sort_strbox_indexed(search_term_ct, max_blen, 0, merged_strbox, id_map);
    assert(search_term_ct);
    const char* duplicate_search_term = scan_for_duplicate_ids(merged_strbox, search_term_ct, max_blen);
    if (duplicate_search_term) {
      LOGERRPRINTFWW("Error: Duplicate term '%s' in --%s column search order.\n", duplicate_search_term, flagname_p);
      goto search_header_line_ret_INVALID_CMDLINE;
    }

    fill_uint_one(search_col_ct, priority_vals);
    fill_ull_one(search_col_ct, cols_and_types);
    uint32_t col_idx = 0;
    while (1) {
      const char* token_end = token_endnn(header_line_iter);
      const uint32_t token_slen = (uintptr_t)(token_end - header_line_iter);
      int32_t ii = bsearch_str(header_line_iter, merged_strbox, token_slen, max_blen, search_term_ct);
      if (ii != -1) {
        const uint32_t cur_map_idx = id_map[(uint32_t)ii];
        const uint32_t search_col_idx = cur_map_idx & 31;
        const uint32_t priority_idx = cur_map_idx >> 5;
        if (priority_vals[search_col_idx] >= priority_idx) {
          if (priority_vals[search_col_idx] == priority_idx) {
            LOGERRPRINTFWW("Error: Duplicate column header '%s' in --%s file.\n", &(merged_strbox[max_blen * cur_map_idx]), flagname_p);
            goto search_header_line_ret_MALFORMED_INPUT;
          }
          priority_vals[search_col_idx] = priority_idx;
          cols_and_types[search_col_idx] = (((uint64_t)col_idx) << 32) | search_col_idx;
        }
      }
      header_line_iter = skip_initial_spaces(token_end);
      if (is_eoln_kns(*header_line_iter)) {
        break;
      }
      ++col_idx;
    }
    uint32_t found_type_bitset = 0;
    for (uint32_t search_col_idx = 0; search_col_idx < search_col_ct; ++search_col_idx) {
      if (priority_vals[search_col_idx] != UINT32_MAX) {
        found_type_bitset |= 1U << search_col_idx;
      }
    }
    const uint32_t found_col_ct = popcount_long(found_type_bitset);
    *found_col_ct_ptr = found_col_ct;
    *found_type_bitset_ptr = found_type_bitset;
    if (found_col_ct) {
#ifdef __cplusplus
      std::sort(cols_and_types, &(cols_and_types[search_col_ct]));
#else
      qsort(cols_and_types, search_col_ct, sizeof(int64_t), uint64cmp);
#endif
      uint32_t prev_col_idx = cols_and_types[0] >> 32;
      col_skips[0] = prev_col_idx;
      col_types[0] = (uint32_t)cols_and_types[0];
      for (uint32_t found_col_idx = 1; found_col_idx < found_col_ct; ++found_col_idx) {
        const uint64_t cur_col_and_type = cols_and_types[found_col_idx];
        const uint32_t cur_col_idx = cur_col_and_type >> 32;
        col_skips[found_col_idx] = cur_col_idx - prev_col_idx;
        col_types[found_col_idx] = (uint32_t)cur_col_and_type;
        prev_col_idx = cur_col_idx;
      }
    }
  }
  while (0) {
  search_header_line_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  search_header_line_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  search_header_line_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
