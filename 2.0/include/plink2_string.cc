// This library is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "plink2_string.h"

#include <assert.h>
#include <stddef.h>  // offsetof()
#include <stdlib.h>  // free()

#ifdef __cplusplus
namespace plink2 {
#endif

#if defined(__LP64__) && !defined(_GNU_SOURCE)
CXXCONST_VOIDP rawmemchr(const void* ss, int cc) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, ss);
  const VecI8* ss_viter = R_CAST(const VecI8*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecI8 vvec_all_needle = veci8_set1(cc);
  VecI8 cur_vvec = *ss_viter;
  VecI8 needle_match_vvec = (cur_vvec == vvec_all_needle);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, ss_viter);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
  uint32_t matching_bytes = veci8_movemask(needle_match_vvec);
  matching_bytes &= UINT32_MAX << leading_byte_ct;
  // This is typically short-range, so the Memrchr() double-vector strategy is
  // unlikely to be profitable.  (todo: experiment with header-inline)
  while (!matching_bytes) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    needle_match_vvec = (cur_vvec == vvec_all_needle);
    matching_bytes = veci8_movemask(needle_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzu32(matching_bytes);
#  else
  uint64_t matching_nybbles = arm_shrn4_i8(needle_match_vvec);
  matching_nybbles &= UINT64_MAX << (4 * leading_byte_ct);
  // This is typically short-range, so the Memrchr() double-vector strategy is
  // unlikely to be profitable.  (todo: experiment with header-inline)
  while (!matching_nybbles) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    needle_match_vvec = (cur_vvec == vvec_all_needle);
    matching_nybbles = arm_shrn4_i8(needle_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzw(matching_nybbles) / 4;
#  endif
  return R_CAST(CXXCONST_VOIDP, R_CAST(uintptr_t, ss_viter) + byte_offset_in_vec);
}
#endif

#ifdef __LP64__
CXXCONST_VOIDP rawmemchr2(const void* ss, unsigned char ucc1, unsigned char ucc2) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, ss);
  const VecI8* ss_viter = R_CAST(const VecI8*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecI8 vvec_all_ucc1 = veci8_set1(ucc1);
  const VecI8 vvec_all_ucc2 = veci8_set1(ucc2);
  VecI8 cur_vvec = *ss_viter;
  VecI8 ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
  VecI8 ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, ss_viter);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
  uint32_t matching_bytes = veci8_movemask(ucc1_match_vvec | ucc2_match_vvec);
  matching_bytes &= UINT32_MAX << leading_byte_ct;
  while (!matching_bytes) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
    ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
    matching_bytes = veci8_movemask(ucc1_match_vvec | ucc2_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzu32(matching_bytes);
#  else
  uint64_t matching_nybbles = arm_shrn4_i8(ucc1_match_vvec | ucc2_match_vvec);
  matching_nybbles &= UINT64_MAX << (4 * leading_byte_ct);
  while (!matching_nybbles) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
    ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
    matching_nybbles = arm_shrn4_i8(ucc1_match_vvec | ucc2_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzw(matching_nybbles) / 4;
#  endif
  return &(DowncastToXC(ss_viter)[byte_offset_in_vec]);
}

CXXCONST_VOIDP rawmemchr3(const void* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, ss);
  const VecI8* ss_viter = R_CAST(const VecI8*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecI8 vvec_all_ucc1 = veci8_set1(ucc1);
  const VecI8 vvec_all_ucc2 = veci8_set1(ucc2);
  const VecI8 vvec_all_ucc3 = veci8_set1(ucc3);
  VecI8 cur_vvec = *ss_viter;
  VecI8 ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
  VecI8 ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
  VecI8 ucc3_match_vvec = (cur_vvec == vvec_all_ucc3);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, ss_viter);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
  uint32_t matching_bytes = veci8_movemask(ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  matching_bytes &= UINT32_MAX << leading_byte_ct;
  while (!matching_bytes) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
    ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
    ucc3_match_vvec = (cur_vvec == vvec_all_ucc3);
    matching_bytes = veci8_movemask(ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzu32(matching_bytes);
#  else
  uint64_t matching_nybbles = arm_shrn4_i8(ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  matching_nybbles &= UINT64_MAX << (4 * leading_byte_ct);
  while (!matching_nybbles) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
    ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
    ucc3_match_vvec = (cur_vvec == vvec_all_ucc3);
    matching_nybbles = arm_shrn4_i8(ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzw(matching_nybbles) / 4;
#  endif
  return &(DowncastToXC(ss_viter)[byte_offset_in_vec]);
}

CXXCONST_CP strchrnul3(const char* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, ss);
  const VecI8* ss_viter = R_CAST(const VecI8*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecI8 vvec_all_zero = veci8_setzero();
  const VecI8 vvec_all_ucc1 = veci8_set1(ucc1);
  const VecI8 vvec_all_ucc2 = veci8_set1(ucc2);
  const VecI8 vvec_all_ucc3 = veci8_set1(ucc3);
  VecI8 cur_vvec = *ss_viter;
  VecI8 zero_match_vvec = (cur_vvec == vvec_all_zero);
  VecI8 ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
  VecI8 ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
  VecI8 ucc3_match_vvec = (cur_vvec == vvec_all_ucc3);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, ss_viter);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
  uint32_t matching_bytes = veci8_movemask(zero_match_vvec | ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  matching_bytes &= UINT32_MAX << leading_byte_ct;
  while (!matching_bytes) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    zero_match_vvec = (cur_vvec == vvec_all_zero);
    ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
    ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
    ucc3_match_vvec = (cur_vvec == vvec_all_ucc3);
    matching_bytes = veci8_movemask(zero_match_vvec | ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzu32(matching_bytes);
#  else
  uint64_t matching_nybbles = arm_shrn4_i8(zero_match_vvec | ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  matching_nybbles &= UINT64_MAX << (4 * leading_byte_ct);
  while (!matching_nybbles) {
    ++ss_viter;
    cur_vvec = *ss_viter;
    zero_match_vvec = (cur_vvec == vvec_all_zero);
    ucc1_match_vvec = (cur_vvec == vvec_all_ucc1);
    ucc2_match_vvec = (cur_vvec == vvec_all_ucc2);
    ucc3_match_vvec = (cur_vvec == vvec_all_ucc3);
    matching_nybbles = arm_shrn4_i8(zero_match_vvec | ucc1_match_vvec | ucc2_match_vvec | ucc3_match_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzw(matching_nybbles) / 4;
#  endif
  return &(DowncastToXC(ss_viter)[byte_offset_in_vec]);
}
#else
CXXCONST_VOIDP rawmemchr2(const void* ss, unsigned char ucc1, unsigned char ucc2) {
  for (const unsigned char* ss_iter = S_CAST(const unsigned char*, ss); ; ++ss_iter) {
    unsigned char ucc = *ss_iter;
    if ((ucc == ucc1) || (ucc == ucc2)) {
      return S_CAST(CXXCONST_VOIDP, ss_iter);
    }
  }
}

CXXCONST_VOIDP rawmemchr3(const void* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3) {
  for (const unsigned char* ss_iter = S_CAST(const unsigned char*, ss); ; ++ss_iter) {
    unsigned char ucc = *ss_iter;
    if ((ucc == ucc1) || (ucc == ucc2) || (ucc == ucc3)) {
      return S_CAST(CXXCONST_VOIDP, ss_iter);
    }
  }
}

CXXCONST_CP strchrnul3(const char* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3) {
  for (; ; ++ss) {
    unsigned char ucc = *ss;
    // at this point, checking against a precomputed 256-byte array may be
    // better?  test this.
    if ((!ucc) || (ucc == ucc1) || (ucc == ucc2) || (ucc == ucc3)) {
      return S_CAST(CXXCONST_CP, ss);
    }
  }
}
#endif

void WordWrap(uint32_t suffix_len, char* strbuf) {
  // Input: A null-terminated string with no intermediate newlines.  If
  //        suffix_len is zero, there should be a terminating \n; otherwise,
  //        the last character should be a space.  The allocation the string is
  //        part of must include at least ~80 bytes past the string end.
  // Effect: Spaces are replaced with newlines in a manner that plays well with
  //         80 column terminal windows.  (Multi-space blocks are never
  //         collapsed.)
  // Considered UTF-8 awareness, but then decided against it after reading
  //   https://denisbider.blogspot.com/2015/09/when-monospace-fonts-arent-unicode.html .
  char* token_start = strbuf;
  char* line_end = &(strbuf[79]);
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

void WordWrapMultiline(char* strbuf) {
  for (char* line_start = strbuf; *line_start != '\0'; ) {
    char* line_end = Strchrnul(line_start, '\n');
    if (*line_end == '\0') {
      line_end[0] = '\n';
      line_end[1] = '\0';
    }
    char* next_line_start = &(line_end[1]);
    // Temporarily clobber first character of next line so WordWrap()
    // conditions are satisfied.
    const char first_char_of_next_line = *next_line_start;
    *next_line_start = '\0';
    WordWrap(0, line_start);
    *next_line_start = first_char_of_next_line;
    line_start = next_line_start;
  }
}

// This implementation is from Kendall Willets.  See
//   https://lemire.me/blog/2021/06/03/computing-the-number-of-digits-of-an-integer-even-faster/

#ifdef USE_AVX2
// bugfix (14 Aug 2022): _lzcnt_u32(0) is 32.
static const uint64_t kLzcntUintSlenTable[] =
  {42949672960LL, 42949672960LL, 41949672960LL, 41949672960LL, 41949672960LL,
   38554705664LL, 38554705664LL, 38554705664LL, 34349738368LL, 34349738368LL,
   34349738368LL, 34349738368LL, 30063771072LL, 30063771072LL, 30063771072LL,
   25769703776LL, 25769703776LL, 25769703776LL, 21474826480LL, 21474826480LL,
   21474826480LL, 21474826480LL, 17179868184LL, 17179868184LL, 17179868184LL,
   12884901788LL, 12884901788LL, 12884901788LL,  8589934582LL,  8589934582LL,
    8589934582LL,  4294967296LL,  4294967296LL};

uint32_t UintSlen(uint32_t num) {
  const uint32_t lz_ct = _lzcnt_u32(num);
  return (num + kLzcntUintSlenTable[lz_ct]) >> 32;
}
#else
static const uint64_t kBsrUintSlenTable[] =
  {4294967296LL,  8589934582LL,  8589934582LL,  8589934582LL,  12884901788LL,
   12884901788LL, 12884901788LL, 17179868184LL, 17179868184LL, 17179868184LL,
   21474826480LL, 21474826480LL, 21474826480LL, 21474826480LL, 25769703776LL,
   25769703776LL, 25769703776LL, 30063771072LL, 30063771072LL, 30063771072LL,
   34349738368LL, 34349738368LL, 34349738368LL, 34349738368LL, 38554705664LL,
   38554705664LL, 38554705664LL, 41949672960LL, 41949672960LL, 41949672960LL,
   42949672960LL, 42949672960LL};

uint32_t UintSlen(uint32_t num) {
  // tried divide-by-10 and divide-by-100 loops, they were slower
  // also tried a hardcoded binary tree, it was better but still slower

  // bsru32(0) is undefined
  num |= 1;
  const uint32_t top_bit_pos = bsru32(num);
  return (num + kBsrUintSlenTable[top_bit_pos]) >> 32;
}
#endif

// May read (kBytesPerWord - 1) bytes past the end of each string.
// This can be quite a bit faster than stock strcmp on x86 for sorting, though
// benchmarking is necessary in general (compiler may perform its own strcmp
// optimizations).
int32_t strcmp_overread(const char* s1, const char* s2) {
#ifndef NO_UNALIGNED
  const uintptr_t* s1_alias = R_CAST(const uintptr_t*, s1);
  const uintptr_t* s2_alias = R_CAST(const uintptr_t*, s2);
  for (uintptr_t widx = 0; ; ++widx) {
    uintptr_t w1 = s1_alias[widx];
    const uintptr_t zcheck = DetectFirstZeroByte(w1);
    uintptr_t w2 = s2_alias[widx];
    if (zcheck) {
      // Mask out bytes past the known null.
      const uintptr_t mask = zcheck ^ (zcheck - 1);
      w1 &= mask;
      w2 &= mask;
      if (w1 == w2) {
        return 0;
      }
    } else if (w1 == w2) {
      continue;
    }
    // bugfix (30 Jun 2018): forgot to adhere to strcmp instead of std::sort
    // interface
#  ifdef __LP64__
    return (__builtin_bswap64(w1) < __builtin_bswap64(w2))? -1 : 1;
#  else
    return (__builtin_bswap32(w1) < __builtin_bswap32(w2))? -1 : 1;
#  endif
  }
#else // NO_UNALIGNED
  return strcmp(s1, s2);
#endif
}

int32_t strcmp_casted(const void* s1, const void* s2) {
  return strcmp(S_CAST(const char*, s1), S_CAST(const char*, s2));
}

int32_t strcmp_overread_casted(const void* s1, const void* s2) {
  return strcmp_overread(S_CAST(const char*, s1), S_CAST(const char*, s2));
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
int32_t strcmp_natural_scan_forward(const char* s1, const char* s2) {
  // assumes s1 and s2 currently point to the middle of a mismatching number,
  // where s1 < s2.
  char c1;
  char c2;
  do {
    c1 = *(++s1);
    c2 = *(++s2);
    if (IsNotDigit(c1)) {
      return -1;
    }
  } while (IsDigit(c2));
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
  unsigned char uc1 = *(++s1);
  unsigned char uc2 = *(++s2);
  while (IsNotNzdigit(uc1) && IsNotNzdigit(uc2)) {
    // state 2
  strcmp_natural_tiebroken_state_2:
    if (uc1 != uc2) {
      if ((uc1 >= 'a') && (uc1 <= 'z')) {
        uc1 -= 32;
      }
      if ((uc2 >= 'a') && (uc2 <= 'z')) {
        uc2 -= 32;
      }
      if (uc1 < uc2) {
        return -1;
      }
      if (uc1 > uc2) {
        return 1;
      }
    } else if (!uc1) {
      return -1;
    }
    uc1 = *(++s1);
    uc2 = *(++s2);
  }
  if (IsNotNzdigit(uc1) || IsNotNzdigit(uc2)) {
    return (uc1 < uc2)? -1 : 1;
  }
  do {
    // state 3
    if (uc1 != uc2) {
      if (IsDigit(uc2)) {
        if (uc1 < uc2) {
          return strcmp_natural_scan_forward(s1, s2);
        }
        return -strcmp_natural_scan_forward(s2, s1);
      }
      return 1;
    }
    uc1 = *(++s1);
    uc2 = *(++s2);
  } while (IsDigit(uc1));
  if (IsDigit(uc2)) {
    return -1;
  }
  // skip the while (is_not_digit...) check
  goto strcmp_natural_tiebroken_state_2;
}

int32_t strcmp_natural_uncasted(const char* s1, const char* s2) {
  unsigned char uc1 = *s1;
  unsigned char uc2 = *s2;
  while (IsNotNzdigit(uc1) && IsNotNzdigit(uc2)) {
    // state 0
  strcmp_natural_uncasted_state_0:
    if (uc1 != uc2) {
      if ((uc1 >= 'a') && (uc1 <= 'z')) {
        if (uc2 + 32 == uc1) {
          return -strcmp_natural_tiebroken(s2, s1);
        }
        if ((uc2 < 'a') || (uc2 > 'z')) {
          uc1 -= 32;
        }
      } else if ((uc2 >= 'a') && (uc2 <= 'z')) {
        uc2 -= 32;
        if (uc1 == uc2) {
          return strcmp_natural_tiebroken(s1, s2);
        }
      }
      return (uc1 < uc2)? -1 : 1;
    }
    if (!uc1) {
      return 0;
    }
    uc1 = *(++s1);
    uc2 = *(++s2);
  }
  if (IsNotNzdigit(uc1) || IsNotNzdigit(uc2)) {
    return (uc1 < uc2)? -1 : 1;
  }
  do {
    // state 1
    if (uc1 != uc2) {
      if (IsDigit(uc2)) {
        if (uc1 < uc2) {
          return strcmp_natural_scan_forward(s1, s2);
        }
        return -strcmp_natural_scan_forward(s2, s1);
      }
      return 1;
    }
    uc1 = *(++s1);
    uc2 = *(++s2);
  } while (IsDigit(uc1));
  if (IsDigit(uc2)) {
    return -1;
  }
  goto strcmp_natural_uncasted_state_0;
}

int32_t strcmp_natural(const void* s1, const void* s2) {
  return strcmp_natural_uncasted(S_CAST(const char*, s1), S_CAST(const char*, s2));
}

int32_t strcmp_deref(const void* s1, const void* s2) {
  return strcmp(*S_CAST(const char* const*, s1), *S_CAST(const char* const*, s2));
}

int32_t strcmp_overread_deref(const void* s1, const void* s2) {
  return strcmp_overread(*S_CAST(const char* const*, s1), *S_CAST(const char* const*, s2));
}

int32_t strcmp_natural_deref(const void* s1, const void* s2) {
  return strcmp_natural_uncasted(*S_CAST(const char* const*, s1), *S_CAST(const char* const*, s2));
}

#ifdef __cplusplus
static_assert(sizeof(Strbuf28Ui) == 32, "Strbuf28Ui is not laid out as expected.");
static_assert(offsetof(Strbuf28Ui, orig_idx) == 28, "Strbuf28Ui is not laid out as expected.");
static_assert(sizeof(Strbuf60Ui) == 64, "Strbuf60Ui is not laid out as expected.");
static_assert(offsetof(Strbuf60Ui, orig_idx) == 60, "Strbuf60Ui is not laid out as expected.");
uintptr_t GetStrboxsortWentryBlen(uintptr_t max_str_blen) {
  if (max_str_blen <= 28) {
    return sizeof(Strbuf28Ui);
  }
  if (max_str_blen <= 60) {
    return sizeof(Strbuf60Ui);
  }
  return max_str_blen;
}
#else
uintptr_t GetStrboxsortWentryBlen(uintptr_t max_str_blen) {
  return MAXV(max_str_blen, sizeof(StrSortIndexedDeref));
}
#endif

// Assumed that sort_wkspace has size >= str_ct *
// max(sizeof(StrSortIndexedDeref), max_str_blen).
// Must be ok to overread.
void SortStrboxIndexed2Fallback(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace) {
  StrSortIndexedDerefOverread* wkspace_alias = S_CAST(StrSortIndexedDerefOverread*, sort_wkspace);
  for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    wkspace_alias[str_idx].strptr = &(strbox[str_idx * max_str_blen]);
    wkspace_alias[str_idx].orig_idx = id_map[str_idx];
  }
  if (!use_nsort) {
    STD_SORT_PAR_UNSEQ(str_ct, strcmp_overread_deref, wkspace_alias);
  } else {
#ifdef __cplusplus
    StrNsortIndexedDeref* wkspace_alias2 = R_CAST(StrNsortIndexedDeref*, wkspace_alias);
    STD_SORT_PAR_UNSEQ(str_ct, nullptr, wkspace_alias2);
#else
    STD_SORT_PAR_UNSEQ(str_ct, strcmp_natural_deref, wkspace_alias);
#endif
  }
  for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    id_map[str_idx] = wkspace_alias[str_idx].orig_idx;
  }
#ifndef __cplusplus
  if (max_str_blen < sizeof(StrSortIndexedDeref)) {
    // actually better to use non-deref sort here, but just get this working
    // properly for now
    for (uint32_t new_idx = 0; new_idx != str_ct; ++new_idx) {
      const char* strptr = wkspace_alias[new_idx].strptr;
      // todo: check whether memcpy(., ., max_str_blen) tends to be better
      strcpy(&(DowncastToC(wkspace_alias)[new_idx * max_str_blen]), strptr);
    }
  } else {
#endif
    // bugfix: need to handle id_map[str_idx] != str_idx
    uint32_t new_idx = str_ct;
    do {
      --new_idx;
      const char* strptr = wkspace_alias[new_idx].strptr;
      strcpy(&(DowncastToC(wkspace_alias)[new_idx * max_str_blen]), strptr);
    } while (new_idx);
#ifndef __cplusplus
  }
#endif
  memcpy(strbox, wkspace_alias, str_ct * max_str_blen);
}

#ifdef __cplusplus
typedef struct WordCmp32bStruct {
  uintptr_t words[32 / kBytesPerWord];
  bool operator<(const struct WordCmp32bStruct& rhs) const {
    uint32_t idx = 0;
    do {
      const uintptr_t cur_word = words[idx];
      const uintptr_t rhs_word = rhs.words[idx];
      if (cur_word != rhs_word) {
#  ifdef __LP64__
        return __builtin_bswap64(cur_word) < __builtin_bswap64(rhs_word);
#  else
        return __builtin_bswap32(cur_word) < __builtin_bswap32(rhs_word);
#  endif
      }
    } while (++idx < (32 / kBytesPerWord));
    return false;
  }
} WordCmp32b;

typedef struct WordCmp64bStruct {
  uintptr_t words[64 / kBytesPerWord];
  bool operator<(const struct WordCmp64bStruct& rhs) const {
    uint32_t idx = 0;
    do {
      const uintptr_t cur_word = words[idx];
      const uintptr_t rhs_word = rhs.words[idx];
      if (cur_word != rhs_word) {
#  ifdef __LP64__
        return __builtin_bswap64(cur_word) < __builtin_bswap64(rhs_word);
#  else
        return __builtin_bswap32(cur_word) < __builtin_bswap32(rhs_word);
#  endif
      }
    } while (++idx < (64 / kBytesPerWord));
    return false;
  }
} WordCmp64b;

static_assert(sizeof(WordCmp32b) == 32, "WordCmp32b does not have the expected size.");
static_assert(sizeof(WordCmp64b) == 64, "WordCmp64b does not have the expected size.");

void SortStrbox32bFinish(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, Strbuf28Ui* filled_wkspace, char* sorted_strbox, uint32_t* id_map) {
  if (!use_nsort) {
    WordCmp32b* wkspace_alias = R_CAST(WordCmp32b*, filled_wkspace);
    STD_SORT_PAR_UNSEQ(str_ct, nullptr, wkspace_alias);
  } else {
    STD_SORT_PAR_UNSEQ(str_ct, nullptr, filled_wkspace);
  }
  for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    strcpy(&(sorted_strbox[str_idx * max_str_blen]), filled_wkspace[str_idx].strbuf);
    id_map[str_idx] = filled_wkspace[str_idx].orig_idx;
  }
}

void SortStrbox64bFinish(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, Strbuf60Ui* filled_wkspace, char* sorted_strbox, uint32_t* id_map) {
  if (!use_nsort) {
    WordCmp64b* wkspace_alias = R_CAST(WordCmp64b*, filled_wkspace);
    STD_SORT_PAR_UNSEQ(str_ct, nullptr, wkspace_alias);
  } else {
    STD_SORT_PAR_UNSEQ(str_ct, nullptr, filled_wkspace);
  }
  for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    strcpy(&(sorted_strbox[str_idx * max_str_blen]), filled_wkspace[str_idx].strbuf);
    id_map[str_idx] = filled_wkspace[str_idx].orig_idx;
  }
}

// Normally use SortStrboxIndexed(), but this version is necessary before
// g_bigstack has been allocated.
void SortStrboxIndexed2(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace) {
  if (max_str_blen <= 28) {
    Strbuf28Ui* wkspace_alias = S_CAST(Strbuf28Ui*, sort_wkspace);
    for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
      const char* cur_str = &(strbox[str_idx * max_str_blen]);
      strcpy(wkspace_alias[str_idx].strbuf, cur_str);
      wkspace_alias[str_idx].orig_idx = id_map[str_idx];
    }
    SortStrbox32bFinish(str_ct, max_str_blen, use_nsort, wkspace_alias, strbox, id_map);
    return;
  }
  if (max_str_blen <= 60) {
    Strbuf60Ui* wkspace_alias = S_CAST(Strbuf60Ui*, sort_wkspace);
    for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
      const char* cur_str = &(strbox[str_idx * max_str_blen]);
      strcpy(wkspace_alias[str_idx].strbuf, cur_str);
      wkspace_alias[str_idx].orig_idx = id_map[str_idx];
    }
    SortStrbox64bFinish(str_ct, max_str_blen, use_nsort, wkspace_alias, strbox, id_map);
    return;
  }
  SortStrboxIndexed2Fallback(str_ct, max_str_blen, use_nsort, strbox, id_map, sort_wkspace);
}
#endif  // __cplusplus

// Must be ok to overread.
BoolErr SortStrboxIndexedMalloc(uintptr_t str_ct, uintptr_t max_str_blen, char* strbox, uint32_t* id_map) {
  if (str_ct < 2) {
    return 0;
  }
  const uintptr_t wkspace_entry_blen = GetStrboxsortWentryBlen(max_str_blen);
  unsigned char* sort_wkspace;
  if (unlikely(pgl_malloc(str_ct * wkspace_entry_blen, &sort_wkspace))) {
    return 1;
  }
  SortStrboxIndexed2(str_ct, max_str_blen, 0, strbox, id_map, sort_wkspace);
  free(sort_wkspace);
  return 0;
}

uint32_t CopyAndDedupSortedStrptrsToStrbox(const char* const* sorted_strptrs, uintptr_t str_ct, uintptr_t max_str_blen, char* strbox) {
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
    if ((cur_slen != prev_slen) || (!memequal(cur_str, prev_str, prev_slen))) {
      memcpy(&(strbox[write_idx * max_str_blen]), cur_str, cur_slen + 1);
      ++write_idx;
      prev_str = cur_str;
    }
  } while (sorted_strptrs_iter != sorted_strptrs_end);
  return write_idx;
}


void StrptrArrSortMain(uintptr_t str_ct, uint32_t overread_ok, uint32_t use_nsort, StrSortIndexedDeref* wkspace_alias) {
  if (!use_nsort) {
    if (overread_ok) {
#ifdef __cplusplus
      StrSortIndexedDerefOverread* wkspace_alias2 = R_CAST(StrSortIndexedDerefOverread*, wkspace_alias);
      STD_SORT_PAR_UNSEQ(str_ct, nullptr, wkspace_alias2);
#else
      STD_SORT_PAR_UNSEQ(str_ct, strcmp_overread_deref, wkspace_alias);
#endif
    } else {
      STD_SORT_PAR_UNSEQ(str_ct, strcmp_deref, wkspace_alias);
    }
  } else {
#ifdef __cplusplus
    StrNsortIndexedDeref* wkspace_alias2 = R_CAST(StrNsortIndexedDeref*, wkspace_alias);
    STD_SORT_PAR_UNSEQ(str_ct, nullptr, wkspace_alias2);
#else
    STD_SORT_PAR_UNSEQ(str_ct, strcmp_natural_deref, wkspace_alias);
#endif
  }
}

void SortStrptrArrIndexed2(uint32_t str_ct, uint32_t leave_first_alone, uint32_t overread_ok, uint32_t use_nsort, const char** strptrs, uint32_t* new_to_old_idx, uint32_t* old_to_new_idx, void* wkspace) {
  const uint32_t str_sort_ct = str_ct - leave_first_alone;
  if (str_sort_ct < 2) {
    if (new_to_old_idx) {
      for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
        new_to_old_idx[str_idx] = str_idx;
      }
    }
    if (old_to_new_idx) {
      for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
        old_to_new_idx[str_idx] = str_idx;
      }
    }
    return;
  }
  StrSortIndexedDeref* wkspace_alias = S_CAST(StrSortIndexedDeref*, wkspace);
  const char** strptrs_to_sort = &(strptrs[leave_first_alone]);
  for (uint32_t str_idx = 0; str_idx != str_sort_ct; ++str_idx) {
    wkspace_alias[str_idx].strptr = strptrs_to_sort[str_idx];
    wkspace_alias[str_idx].orig_idx = str_idx + leave_first_alone;
  }
  StrptrArrSortMain(str_sort_ct, overread_ok, use_nsort, wkspace_alias);
  if (leave_first_alone) {
    if (new_to_old_idx) {
      new_to_old_idx[0] = 0;
      new_to_old_idx = &(new_to_old_idx[1]);
    }
    if (old_to_new_idx) {
      old_to_new_idx[0] = 0;
    }
  }
  for (uint32_t str_idx = 0; str_idx != str_sort_ct; ++str_idx) {
    strptrs_to_sort[str_idx] = wkspace_alias[str_idx].strptr;
    const uint32_t orig_idx = wkspace_alias[str_idx].orig_idx;
    if (new_to_old_idx) {
      new_to_old_idx[str_idx] = orig_idx;
    }
    if (old_to_new_idx) {
      old_to_new_idx[orig_idx] = str_idx + leave_first_alone;
    }
  }
}

/*
uint32_t match_upper(const char* str_iter, const char* fixed_str) {
  char cc = *fixed_str++;
  do {
    if ((((unsigned char)(*str_iter++)) & 0xdf) != ((unsigned char)cc)) {
      return 0;
    }
    cc = *fixed_str++;
  } while (cc);
  return !(*str_iter);
}
*/

/*
#ifdef __LP64__
CXXCONST_CP FirstPrecharFar(const char* str_iter, uint32_t char_code) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, str_iter);
  const VecUc* str_viter = R_CAST(const VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecUc vvec_add = vecuc_set1(128 - char_code);
  VecUc cur_vvec = *str_viter;
  VecUc non_prechar_vvec = vecuc_adds(cur_vvec, vvec_add);
  uint32_t matching_bytes = ~vecuc_movemask(non_prechar_vvec);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, str_viter);
  matching_bytes &= UINT32_MAX << leading_byte_ct;
  // This typically isn't *that* long-range, so the Memrchr()
  // double-vector strategy is unlikely to be profitable?
  while (!matching_bytes) {
    ++str_viter;
    cur_vvec = *str_viter;
    non_prechar_vvec = vecuc_adds(cur_vvec, vvec_add);
    matching_bytes = ~vecuc_movemask(non_prechar_vvec);
  }
  const uint32_t byte_offset_in_vec = ctzu32(matching_bytes);
  return &(DowncastToXC(str_viter)[byte_offset_in_vec]);
}
#endif
*/

uint32_t MatchUpperCounted(const char* str, const char* fixed_str, uint32_t ct) {
  for (uint32_t uii = 0; uii != ct; ++uii) {
    if ((ctou32(str[uii]) & 0xdf) != ctou32(fixed_str[uii])) {
      return 0;
    }
  }
  return 1;
}

// might want to make this std::array in the future?
static const char kToUpper[256] = {
  '\0', '\1', '\2', '\3', '\4', '\5', '\6', '\7',
  '\10', '\11', '\12', '\13', '\14', '\15', '\16', '\17',
  '\20', '\21', '\22', '\23', '\24', '\25', '\26', '\27',
  '\30', '\31', '\32', '\33', '\34', '\35', '\36', '\37',
  '\40', '\41', '\42', '\43', '\44', '\45', '\46', '\47',
  '\50', '\51', '\52', '\53', '\54', '\55', '\56', '\57',
  '\60', '\61', '\62', '\63', '\64', '\65', '\66', '\67',
  '\70', '\71', '\72', '\73', '\74', '\75', '\76', '\77',
  '\100', '\101', '\102', '\103', '\104', '\105', '\106', '\107',
  '\110', '\111', '\112', '\113', '\114', '\115', '\116', '\117',
  '\120', '\121', '\122', '\123', '\124', '\125', '\126', '\127',
  '\130', '\131', '\132', '\133', '\134', '\135', '\136', '\137',
  '\140', '\101', '\102', '\103', '\104', '\105', '\106', '\107',
  '\110', '\111', '\112', '\113', '\114', '\115', '\116', '\117',
  '\120', '\121', '\122', '\123', '\124', '\125', '\126', '\127',
  '\130', '\131', '\132', '\173', '\174', '\175', '\176', '\177',
  '\200', '\201', '\202', '\203', '\204', '\205', '\206', '\207',
  '\210', '\211', '\212', '\213', '\214', '\215', '\216', '\217',
  '\220', '\221', '\222', '\223', '\224', '\225', '\226', '\227',
  '\230', '\231', '\232', '\233', '\234', '\235', '\236', '\237',
  '\240', '\241', '\242', '\243', '\244', '\245', '\246', '\247',
  '\250', '\251', '\252', '\253', '\254', '\255', '\256', '\257',
  '\260', '\261', '\262', '\263', '\264', '\265', '\266', '\267',
  '\270', '\271', '\272', '\273', '\274', '\275', '\276', '\277',
  '\300', '\301', '\302', '\303', '\304', '\305', '\306', '\307',
  '\310', '\311', '\312', '\313', '\314', '\315', '\316', '\317',
  '\320', '\321', '\322', '\323', '\324', '\325', '\326', '\327',
  '\330', '\331', '\332', '\333', '\334', '\335', '\336', '\337',
  '\340', '\341', '\342', '\343', '\344', '\345', '\346', '\347',
  '\350', '\351', '\352', '\353', '\354', '\355', '\356', '\357',
  '\360', '\361', '\362', '\363', '\364', '\365', '\366', '\367',
  '\370', '\371', '\372', '\373', '\374', '\375', '\376', '\377'
};

uint32_t strcaseequal(const char* str1, const char* str2, uint32_t ct) {
  for (uint32_t uii = 0; uii != ct; ++uii) {
    if (kToUpper[ctou32(str1[uii])] != kToUpper[ctou32(str2[uii])]) {
      return 0;
    }
  }
  return 1;
}

/*
void str_toupper(char* str_iter) {
  while (1) {
    const uint32_t uii = (unsigned char)(*str_iter);
    if (!uii) {
      return;
    }
    if (((uint32_t)(uii - 97)) < 26) {
      // 'a' has ASCII code 97
      *str_iter = uii - 32;
    }
    ++str_iter;
  }
}

void buf_toupper(uint32_t slen, char* strbuf) {
  for (uint32_t pos = 0; pos != slen; ++pos) {
    const uint32_t uii = (unsigned char)(strbuf[pos]);
    if (((uint32_t)(uii - 97)) < 26) {
      strbuf[pos] = uii - 32;
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

char* memcpya_toupper(char* __restrict target, const char* __restrict source, uint32_t slen) {
  for (uint32_t pos = 0; pos != slen; ++pos) {
    uint32_t uii = ctou32(source[pos]);
    if ((uii - 97) < 26) {
      uii -= 32;
    }
    target[pos] = uii;
  }
  return &(target[slen]);
}
*/


uint32_t IsAlphanumeric(const char* str_iter) {
  for (; ; ++str_iter) {
    uint32_t uii = ctou32(*str_iter);
    if (!uii) {
      return 1;
    }
    if (((uii - 48) > 9) && (((uii & 0xffffffdfU) - 65) > 25)) {
      return 0;
    }
  }
}

BoolErr ScanPosintptr(const char* str_iter, uintptr_t* valp) {
  // Reads an integer in [1, 2^kBitsPerWord - 1].  Assumes first character is
  // nonspace.
  assert(ctow(str_iter[0]) > 32);
  uintptr_t val = ctow(*str_iter++) - 48;
  if (val >= 10) {
#ifdef __LP64__
    if (unlikely(val != 0xfffffffffffffffbLLU)) {
      return 1;
    }
#else
    if (unlikely(val != 0xfffffffbU)) {
      return 1;
    }
#endif
    val = ctow(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  while (!val) {
    val = ctow(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
#ifdef __LP64__
  // limit is 20 digits, we've already read one
  const char* str_limit = &(str_iter[20]);
#else
  const char* str_limit = &(str_iter[10]);
#endif
  while (1) {
    const uintptr_t cur_digit = ctow(*str_iter++) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    const uintptr_t cur_digit2 = ctow(*str_iter++) - 48;
    if (str_iter == str_limit) {
      if (unlikely((cur_digit2 < 10) || ((val >= (~k0LU) / 10) && ((val > (~k0LU) / 10) || (cur_digit > (~k0LU) % 10))))) {
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

// Requires (cap+1) * 100 < 2^64.
static inline BoolErr ScanmovU64CappedFinish(uint64_t cap, const char** str_iterp, uint64_t* valp) {
  const char* str_iter = *str_iterp;
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = ctou64(*str_iter++) - 48;
    if (cur_digit >= 10) {
      break;
    }
    // val = val * 10 + cur_digit;
    const uint64_t cur_digit2 = ctou64(*str_iter++) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (unlikely(val > cap)) {
        return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (unlikely(val > cap)) {
      return 1;
    }
  }
  *valp = val;
  *str_iterp = &(str_iter[-1]);
  return 0;
}

BoolErr ScanmovU64Capped(uint64_t cap, const char** str_iterp, uint64_t* valp) {
  const char* str_iter = *str_iterp;
  *valp = ctou32(*str_iter++) - 48;
  if (*valp >= 10) {
    if (unlikely(*valp != 0xfffffffbU)) {
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely((*valp) >= 10)) {
      return 1;
    }
  }
  *str_iterp = str_iter;
  return ScanmovU64CappedFinish(cap, str_iterp, valp);
}

#ifdef __LP64__
static inline BoolErr ScanmovUintCappedFinish(uint64_t cap, const char** str_iterp, uint32_t* valp) {
  const char* str_iter = *str_iterp;
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = ctou64(*str_iter++) - 48;
    if (cur_digit >= 10) {
      break;
    }
    // val = val * 10 + cur_digit;
    const uint64_t cur_digit2 = ctou64(*str_iter++) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (unlikely(val > cap)) {
        return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (unlikely(val > cap)) {
      return 1;
    }
  }
  *valp = val;
  *str_iterp = &(str_iter[-1]);
  return 0;
}

BoolErr ScanmovPosintCapped(uint64_t cap, const char** str_iterp, uint32_t* valp) {
  const char* str_iter = *str_iterp;
  *valp = ctou32(*str_iter++) - 48;
  if (*valp >= 10) {
    if (unlikely(*valp != 0xfffffffbU)) {
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely((*valp) >= 10)) {
      return 1;
    }
  }
  *str_iterp = str_iter;
  return ScanmovUintCappedFinish(cap, str_iterp, valp);
}

BoolErr ScanmovUintCapped(uint64_t cap, const char** str_iterp, uint32_t* valp) {
  const char* str_iter = *str_iterp;
  *valp = ctou32(*str_iter++) - 48;
  if (*valp >= 10) {
    if (*valp != 0xfffffffbU) {
      // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
      if (unlikely((*valp != 0xfffffffdU) || (*str_iter != '0'))) {
        return 1;
      }
      // accept "-0", "-00", etc.
      while (*(++str_iter) == '0');
      *valp = 0;
      *str_iterp = str_iter;
      return (ctou32(*str_iter) - 48) < 10;
    }
    // accept leading '+'
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  *str_iterp = str_iter;
  return ScanmovUintCappedFinish(cap, str_iterp, valp);
}

// 2^{-31} < floor <= 0 <= cap < 2^31
BoolErr ScanmovIntBounded(uint64_t abs_floor, uint64_t cap, const char** str_iterp, int32_t* valp) {
  const char* str_iter = *str_iterp;
  *valp = ctou32(*str_iter++) - 48;
  int32_t sign = 1;
  if (ctou32(*valp) >= 10) {
    if (*valp == -3) {
      sign = -1;
      cap = abs_floor;
    } else if (unlikely(*valp != -5)) {  // accept leading '+'
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  *str_iterp = str_iter;
  if (ScanmovUintCappedFinish(cap, str_iterp, I32ToU32(valp))) {
    return 1;
  }
  *valp *= sign;
  return 0;
}
#else
BoolErr ScanmovPosintCapped32(uint32_t cap_div_10, uint32_t cap_mod_10, const char** str_iterp, uint32_t* valp) {
  const char* str_iter = *str_iterp;
  uint32_t val = ctou32(*str_iter++) - 48;
  if (val >= 10) {
    if (unlikely(val != 0xfffffffbU)) {
      return 1;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  while (!val) {
    val = ctou32(*str_iter++);
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  for (; ; ++str_iter) {
    const uint32_t cur_digit = ctou32(*str_iter) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      *str_iterp = str_iter;
      return 0;
    }
    if (unlikely((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10)))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

BoolErr ScanmovUintCapped32(uint32_t cap_div_10, uint32_t cap_mod_10, const char** str_iterp, uint32_t* valp) {
  const char* str_iter = *str_iterp;
  uint32_t val = ctou32(*str_iter++) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      if (unlikely((val != 0xfffffffd) || (*str_iter != '0'))) {
        return 1;
      }
      while (*(++str_iter) == '0');
      *valp = 0;
      *str_iterp = str_iter;
      return (ctou32(*str_iter) - 48) < 10;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  for (; ; ++str_iter) {
    const uint32_t cur_digit = ctou32(*str_iter) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      *str_iterp = str_iter;
      return 0;
    }
    if (unlikely((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10)))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

BoolErr ScanmovIntBounded32(uint32_t abs_floor_div_10, uint32_t abs_floor_mod_10, uint32_t cap_div_10, uint32_t cap_mod_10, const char** str_iterp, int32_t* valp) {
  const char* str_iter = *str_iterp;
  uint32_t val = ctou32(*str_iter++) - 48;
  int32_t sign = 1;
  if (val >= 10) {
    if (val == 0xfffffffdU) {
      sign = -1;
      cap_div_10 = abs_floor_div_10;
      cap_mod_10 = abs_floor_mod_10;
    } else if (unlikely(val != 0xfffffffbU)) {
      return 1;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  for (; ; ++str_iter) {
    const uint32_t cur_digit = ctou32(*str_iter) - 48;
    if (cur_digit >= 10) {
      *valp = sign * S_CAST(int32_t, val);
      *str_iterp = str_iter;
      return 0;
    }
    if (unlikely((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10)))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}
#endif

static const double kPositivePow10[] = {1, 1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10, 1.0e11, 1.0e12, 1.0e13, 1.0e14, 1.0e15};
static const double kPositivePowTen16[] = {1, 1.0e16, 1.0e32, 1.0e48, 1.0e64, 1.0e80, 1.0e96, 1.0e112, 1.0e128, 1.0e144, 1.0e160, 1.0e176, 1.0e192, 1.0e208, 1.0e224, 1.0e240};
static const double kNegativePow10[] = {1, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14, 1.0e-15};
static const double kNegativePowTen16[] = {1, 1.0e-16, 1.0e-32, 1.0e-48, 1.0e-64, 1.0e-80, 1.0e-96, 1.0e-112};

CXXCONST_CP ScanadvDouble(const char* str_iter, double* valp) {
  // requires first character to be nonspace (to succeed; it fails without
  //   segfaulting on space/eoln/null)
  // don't care about hexadecimal
  // ok to lose last ~2 bits of precision
  // ok if behavior undefined on >1GB strings in 32-bit case, >2GB for 64-bit
  // fail on nan/infinity/overflow instead of usual strtod behavior
  uint32_t cur_char_code = ctou32(*str_iter);
  const uint32_t is_negative = (cur_char_code == 45);
  if (is_negative || (cur_char_code == 43)) {
    cur_char_code = ctou32(*(++str_iter));
  }
  uint32_t cur_digit = cur_char_code - 48;
  intptr_t e10 = 0;
  const char* dot_ptr;
  int64_t digits;
#ifdef __LP64__
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvDouble_parse_decimal;
        }
        goto ScanadvDouble_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    // (could keep ~19 instead, but if we're systematically losing the last two
    // bits of precision anyway...)
    const char* last_sig_fig_ptr = str_iter;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
    } while (cur_digit < 10);
    e10 = S_CAST(intptr_t, str_iter - last_sig_fig_ptr) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
    }
    goto ScanadvDouble_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = str_iter;
  cur_digit = ctou32(*(++str_iter)) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits = cur_digit;
 ScanadvDouble_parse_decimal:
  while (1) {
    cur_digit = ctou32(*(++str_iter)) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
      break;
    }
    digits = digits * 10 + cur_digit;
    if (digits >= 10000000000000000LL) {
      e10 = -S_CAST(intptr_t, str_iter - dot_ptr);
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
      break;
    }
  }
 ScanadvDouble_parse_exponent:
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = ctou32(*(++str_iter));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = ctou32(*(++str_iter));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 214748364) {
        // may as well guard against exponent overflow
        if (!exp_is_negative) {
          return nullptr;
        }
        *valp = 0;
        do {
          cur_digit = ctou32(*(++str_iter)) - 48;
        } while (cur_digit < 10);
        return S_CAST(CXXCONST_CP, str_iter);
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ctou32(*(++str_iter)) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#else  // not __LP64__
  int32_t digits_short;
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits_short = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvDouble_parse_decimal;
        }
        digits = digits_short;
        goto ScanadvDouble_parse_exponent;
      }
      digits_short = digits_short * 10 + cur_digit;
    } while (digits_short < 100000000);
    digits = digits_short;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvDouble_parse_decimal_long;
        }
        goto ScanadvDouble_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    const char* last_sig_fig_ptr = str_iter;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
    } while (cur_digit < 10);
    e10 = S_CAST(intptr_t, str_iter - last_sig_fig_ptr) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
    }
    goto ScanadvDouble_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = str_iter;
  cur_digit = ctou32(*(++str_iter)) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits_short = cur_digit;
 ScanadvDouble_parse_decimal:
  while (1) {
    cur_digit = ctou32(*(++str_iter)) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
      digits = digits_short;
      break;
    }
    digits_short = digits_short * 10 + cur_digit;
    if (digits_short >= 100000000) {
      digits = digits_short;
    ScanadvDouble_parse_decimal_long:
      while (1) {
        cur_digit = ctou32(*(++str_iter)) - 48;
        if (cur_digit >= 10) {
          e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
          goto ScanadvDouble_parse_exponent;
        }
        digits = digits * 10 + cur_digit;
        if (digits >= 10000000000000000LL) {
          e10 = -S_CAST(intptr_t, str_iter - dot_ptr);
          do {
            cur_digit = ctou32(*(++str_iter)) - 48;
          } while (cur_digit < 10);
          goto ScanadvDouble_parse_exponent;
        }
      }
    }
  }
 ScanadvDouble_parse_exponent:
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = ctou32(*(++str_iter));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = ctou32(*(++str_iter));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 107374182) {
        // may as well guard against exponent overflow
        // 2^30 instead of 2^31 limit since this gets added to e10
        if (!exp_is_negative) {
          return nullptr;
        }
        *valp = 0;
        do {
          cur_digit = ctou32(*(++str_iter)) - 48;
        } while (cur_digit < 10);
        return S_CAST(CXXCONST_CP, str_iter);
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ctou32(*(++str_iter)) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#endif
  if (digits == 0) {
    *valp = 0;
    return S_CAST(CXXCONST_CP, str_iter);
  }
  if (is_negative) {
    digits = -digits;
  }
  double dxx = S_CAST(double, digits);
  if (e10) {
    if (e10 < 0) {
      uint32_t pos_exp = -e10;
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
      uint32_t pos_exp = e10;
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
  return S_CAST(CXXCONST_CP, str_iter);
}

CXXCONST_CP ScanadvLn(const char* str_iter, double* ln_ptr) {
  // revised ScanadvDouble() which currently requires number to be nonnegative
  // returns -DBL_MAX on 0
  uint32_t cur_char_code = ctou32(*str_iter);
  const uint32_t is_negative = (cur_char_code == 45);
  if (is_negative || (cur_char_code == 43)) {
    cur_char_code = ctou32(*(++str_iter));
  }
  uint32_t cur_digit = cur_char_code - 48;
  intptr_t e10 = 0;
  const char* dot_ptr;
  int64_t digits;
#ifdef __LP64__
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvLn_parse_decimal;
        }
        goto ScanadvLn_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    // (could keep ~19 instead, but if we're systematically losing the last two
    // bits of precision anyway...)
    const char* last_sig_fig_ptr = str_iter;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
    } while (cur_digit < 10);
    e10 = S_CAST(intptr_t, str_iter - last_sig_fig_ptr) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
    }
    goto ScanadvLn_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = str_iter;
  cur_digit = ctou32(*(++str_iter)) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits = cur_digit;
 ScanadvLn_parse_decimal:
  while (1) {
    cur_digit = ctou32(*(++str_iter)) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
      break;
    }
    digits = digits * 10 + cur_digit;
    if (digits >= 10000000000000000LL) {
      e10 = -S_CAST(intptr_t, str_iter - dot_ptr);
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
      break;
    }
  }
 ScanadvLn_parse_exponent:
  if (is_negative && (digits != 0)) {
    return nullptr;
  }
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = ctou32(*(++str_iter));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = ctou32(*(++str_iter));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 214748364) {
        // may as well guard against exponent overflow
        if (!exp_is_negative) {
          return nullptr;
        }
        *ln_ptr = -DBL_MAX;
        do {
          cur_digit = ctou32(*(++str_iter)) - 48;
        } while (cur_digit < 10);
        return S_CAST(CXXCONST_CP, str_iter);
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ctou32(*(++str_iter)) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#else  // not __LP64__
  int32_t digits_short;
  if (cur_digit < 10) {
    // ok, we have at least one digit
    digits_short = cur_digit;
    // to check: best to skip leading zeroes and compare against 17 instead of
    // 10^16?
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvLn_parse_decimal;
        }
        digits = digits_short;
        goto ScanadvLn_parse_exponent;
      }
      digits_short = digits_short * 10 + cur_digit;
    } while (digits_short < 100000000);
    digits = digits_short;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
      if (cur_digit >= 10) {
        if (cur_digit == 0xfffffffeU) {
          dot_ptr = str_iter;
          goto ScanadvLn_parse_decimal_long;
        }
        goto ScanadvLn_parse_exponent;
      }
      digits = digits * 10 + cur_digit;
    } while (digits < 10000000000000000LL);
    // we have 17 significant digits; count the rest, but don't worry about
    // contents
    const char* last_sig_fig_ptr = str_iter;
    do {
      cur_digit = ctou32(*(++str_iter)) - 48;
    } while (cur_digit < 10);
    e10 = S_CAST(intptr_t, str_iter - last_sig_fig_ptr) - 1;
    if (cur_digit == 0xfffffffeU) {
      do {
        cur_digit = ctou32(*(++str_iter)) - 48;
      } while (cur_digit < 10);
    }
    goto ScanadvLn_parse_exponent;
  }
  if (cur_digit != 0xfffffffeU) {
    return nullptr;
  }
  // first (nonsign) character is dot, verify we have a digit after it
  dot_ptr = str_iter;
  cur_digit = ctou32(*(++str_iter)) - 48;
  if (cur_digit >= 10) {
    return nullptr;
  }
  digits_short = cur_digit;
 ScanadvLn_parse_decimal:
  while (1) {
    cur_digit = ctou32(*(++str_iter)) - 48;
    if (cur_digit >= 10) {
      e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
      digits = digits_short;
      break;
    }
    digits_short = digits_short * 10 + cur_digit;
    if (digits_short >= 100000000) {
      digits = digits_short;
    ScanadvLn_parse_decimal_long:
      while (1) {
        cur_digit = ctou32(*(++str_iter)) - 48;
        if (cur_digit >= 10) {
          e10 = 1 - S_CAST(intptr_t, str_iter - dot_ptr);
          goto ScanadvLn_parse_exponent;
        }
        digits = digits * 10 + cur_digit;
        if (digits >= 10000000000000000LL) {
          e10 = -S_CAST(intptr_t, str_iter - dot_ptr);
          do {
            cur_digit = ctou32(*(++str_iter)) - 48;
          } while (cur_digit < 10);
          goto ScanadvLn_parse_exponent;
        }
      }
    }
  }
 ScanadvLn_parse_exponent:
  if (is_negative && (digits != 0)) {
    return nullptr;
  }
  if ((cur_digit & 0xdf) == 21) { // 'E' - '0' is 21
    cur_char_code = ctou32(*(++str_iter));
    const uint32_t exp_is_negative = (cur_char_code == 45);
    if (exp_is_negative || (cur_char_code == 43)) {
      cur_char_code = ctou32(*(++str_iter));
    }
    cur_digit = cur_char_code - 48;
    int32_t cur_exp = 0;
    while (cur_digit < 10) {
      if (cur_exp >= 107374182) {
        // may as well guard against exponent overflow
        if (!exp_is_negative) {
          return nullptr;
        }
        *ln_ptr = -DBL_MAX;
        do {
          cur_digit = ctou32(*(++str_iter)) - 48;
        } while (cur_digit < 10);
        return S_CAST(CXXCONST_CP, str_iter);
      }
      cur_exp = cur_exp * 10 + cur_digit;
      cur_digit = ctou32(*(++str_iter)) - 48;
    }
    if (exp_is_negative) {
      cur_exp = -cur_exp;
    }
    e10 += cur_exp;
  }
#endif
  if (digits == 0) {
    *ln_ptr = -DBL_MAX;
    return S_CAST(CXXCONST_CP, str_iter);
  }
  double ln_val = log(S_CAST(double, digits));
  if (e10) {
    ln_val += e10 * kLn10;
  }
  *ln_ptr = ln_val;
  return S_CAST(CXXCONST_CP, str_iter);
}

BoolErr ScanPosintCappedx(const char* str_iter, uint64_t cap, uint32_t* valp) {
  double val;
  if ((!ScantokDouble(str_iter, &val)) || (val < 1.0) || (val > S_CAST(double, cap))) {
    return 1;
  }
  *valp = S_CAST(uint32_t, val);
  return (val != S_CAST(double, *valp));
}

BoolErr ScanUintCappedx(const char* str_iter, uint64_t cap, uint32_t* valp) {
  double val;
  if ((!ScantokDouble(str_iter, &val)) || (val < 0.0) || (val > S_CAST(double, cap))) {
    return 1;
  }
  *valp = S_CAST(uint32_t, val);
  return (val != S_CAST(double, *valp));
}

BoolErr ScanIntAbsBoundedx(const char* str_iter, int64_t bound, int32_t* valp) {
  const double bound_d = S_CAST(double, bound);
  double val;
  if ((!ScantokDouble(str_iter, &val)) || (val < -bound_d) || (val > bound_d)) {
    return 1;
  }
  *valp = S_CAST(int32_t, val);
  return (val != S_CAST(double, *valp));
}

BoolErr ScanPosintptrx(const char* str_iter, uintptr_t* valp) {
  double val;
  if ((!ScantokDouble(str_iter, &val)) || (val < 1.0) || (val > S_CAST(double, ~k0LU))) {
    return 1;
  }
  *valp = S_CAST(uintptr_t, val);
  return (val != S_CAST(double, *valp));
}

void GetTopTwoUi(const uint32_t* __restrict uint_arr, uintptr_t uia_size, uintptr_t* __restrict top_idx_ptr, uintptr_t* __restrict second_idx_ptr) {
  assert(uia_size > 1);
  uintptr_t top_idx = (uint_arr[1] > uint_arr[0])? 1 : 0;
  uintptr_t second_idx = 1 ^ top_idx;
  uint32_t top_val = uint_arr[top_idx];
  uint32_t second_val = uint_arr[second_idx];
  uintptr_t cur_idx;
  uintptr_t cur_val;
  for (cur_idx = 2; cur_idx != uia_size; ++cur_idx) {
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

#ifdef USE_AVX2
CXXCONST_CP NextTokenMultFar(const char* str_iter, uint32_t ct) {
  // assert(ct);
  if (!str_iter) {
    return nullptr;
  }
  uint32_t transition_bits_to_skip_p1 = ct * 2;
  const uintptr_t starting_addr = R_CAST(uintptr_t, str_iter);
  const VecUc* str_viter_start = R_CAST(const VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecUc vvec_all_tab = vecuc_set1(9);
  const VecUc vvec_all_space = vecuc_set1(32);
  VecUc cur_vvec = *str_viter_start;
  VecUc tabspace_vvec = (cur_vvec == vvec_all_tab) | (cur_vvec == vvec_all_space);
  // Underlying intrinsics are _mm256_cmpeq_... and _mm256_cmpgt_...
  // So we don't really want to use <= or !=.
  // (Er, there's also an signed vs. unsigned comparison impedence mismatch:
  // there is no unsigned version of _mm256_cmpgt_epi8.  So we don't really
  // want to use gcc-vector-extension < at all here when we can use a mask to
  // do the same job.)
  const VecUc vvec_all96 = vecuc_set1(96);
  VecUc post31_vvec = vecuc_adds(cur_vvec, vvec_all96);
  uint32_t delimiter_bytes = vecuc_movemask(tabspace_vvec);
  const uint32_t initial_nonterminating_bytes = vecuc_movemask(post31_vvec) | delimiter_bytes;
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, str_viter_start);
  uint32_t leading_mask = UINT32_MAX << leading_byte_ct;
  delimiter_bytes &= leading_mask;
  uint32_t terminating_bytes = leading_mask & (~initial_nonterminating_bytes);
  uint32_t prev_delim_highbit = 0;
  for (const VecUc* str_viter = str_viter_start; ; ) {
    // The number of 0->1 + 1->0 bit transitions in x is
    //   popcount(x ^ (x << 1)).
    uint32_t cur_transitions = S_CAST(Vec8thUint, delimiter_bytes ^ (delimiter_bytes << 1) ^ prev_delim_highbit);
    uint32_t cur_transition_ct = PopcountVec8thUint(cur_transitions);
    if (cur_transition_ct >= transition_bits_to_skip_p1) {
      cur_transitions = ClearBottomSetBits(transition_bits_to_skip_p1 - 1, cur_transitions);
      uint32_t byte_offset_in_vec = ctzu32(cur_transitions);
      if (terminating_bytes << (31 - byte_offset_in_vec)) {
        return nullptr;
      }
      return &(DowncastToXC(str_viter)[byte_offset_in_vec]);
    }
    if (terminating_bytes) {
      return nullptr;
    }
    transition_bits_to_skip_p1 -= cur_transition_ct;
    prev_delim_highbit = transition_bits_to_skip_p1 & 1;
    ++str_viter;
    cur_vvec = *str_viter;
    tabspace_vvec = (cur_vvec == vvec_all_tab) | (cur_vvec == vvec_all_space);
    post31_vvec = vecuc_adds(cur_vvec, vvec_all96);
    delimiter_bytes = vecuc_movemask(tabspace_vvec);
    terminating_bytes = ~(vecuc_movemask(post31_vvec) | delimiter_bytes);
  }
}

const char* TokenLexK0(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens) {
  const uint32_t relevant_col_ct_m1 = relevant_col_ct - 1;
  uint32_t relevant_col_idx = 0;
  uint32_t transition_bits_to_skip_p1 = col_skips[relevant_col_idx] * 2;
  const uintptr_t starting_addr = R_CAST(uintptr_t, str_iter);
  const VecUc* str_viter = R_CAST(const VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecUc vvec_all_tab = vecuc_set1(9);
  const VecUc vvec_all_space = vecuc_set1(32);
  VecUc cur_vvec = *str_viter;
  VecUc tabspace_vvec = (cur_vvec == vvec_all_tab) | (cur_vvec == vvec_all_space);
  const VecUc vvec_all96 = vecuc_set1(96);
  VecUc post31_vvec = vecuc_adds(cur_vvec, vvec_all96);
  uint32_t delimiter_bytes = vecuc_movemask(tabspace_vvec);
  const uint32_t initial_nonterminating_bytes = vecuc_movemask(post31_vvec) | delimiter_bytes;
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, str_viter);
  uint32_t leading_mask = UINT32_MAX << leading_byte_ct;
  delimiter_bytes &= leading_mask;
  uint32_t terminating_bytes = leading_mask & (~initial_nonterminating_bytes);
  uint32_t prev_delim_highbit = 0;
  if (!transition_bits_to_skip_p1) {
    // Special case: col_skips[0] can be zero.  When it is, str_iter must point
    // to the beginning of a token.
    prev_delim_highbit = -leading_mask;
    transition_bits_to_skip_p1 = 1;
  }
  while (1) {
    // The number of 0->1 + 1->0 bit transitions in x is
    //   popcount(x ^ (x << 1)).
    uint32_t cur_transitions = S_CAST(Vec8thUint, delimiter_bytes ^ (delimiter_bytes << 1) ^ prev_delim_highbit);
    uint32_t cur_transition_ct = PopcountVec8thUint(cur_transitions);
    while (cur_transition_ct >= transition_bits_to_skip_p1) {
      cur_transitions = ClearBottomSetBits(transition_bits_to_skip_p1 - 1, cur_transitions);
      uint32_t byte_offset_in_vec = ctzu32(cur_transitions);
      const char* token_start = &(DowncastKToC(str_viter)[byte_offset_in_vec]);
      const uint32_t cur_col_type = col_types[relevant_col_idx];
      token_ptrs[cur_col_type] = token_start;
      // Find end of token.
      if (relevant_col_idx == relevant_col_ct_m1) {
        if (terminating_bytes << (31 - byte_offset_in_vec)) {
          return nullptr;
        }
        // Handle this separately, since it's ok for the last token to be
        // terminated by something other than tab/space.
        // cur_transitions is a misnomer here; we just need the bottom bit to
        // tell us where the token end is.

        // 33..255 is a bit more annoying to check for efficiently than
        // 32..255.
        const VecUc vvec_all95 = vecuc_set1(95);
        cur_transitions &= cur_transitions - 1;
        cur_transitions |= terminating_bytes;
        while (!cur_transitions) {
          ++str_viter;
          cur_vvec = *str_viter;
          const VecUc postspace_vvec = vecuc_adds(cur_vvec, vvec_all95);
          cur_transitions = S_CAST(Vec8thUint, ~vecuc_movemask(postspace_vvec));
        }
        byte_offset_in_vec = ctzu32(cur_transitions);
        const char* token_end = &(DowncastKToC(str_viter)[byte_offset_in_vec]);
        token_slens[cur_col_type] = token_end - token_start;
        return token_end;
      }
      cur_transition_ct -= transition_bits_to_skip_p1;
      if (!cur_transition_ct) {
        // Token end is in a later vector.
        do {
          if (terminating_bytes) {
            return nullptr;
          }
          ++str_viter;
          cur_vvec = *str_viter;
          tabspace_vvec = (cur_vvec == vvec_all_tab) | (cur_vvec == vvec_all_space);
          post31_vvec = vecuc_adds(cur_vvec, vvec_all96);
          delimiter_bytes = vecuc_movemask(tabspace_vvec);
          terminating_bytes = ~(vecuc_movemask(post31_vvec) | delimiter_bytes);
        } while (!delimiter_bytes);
        cur_transitions = S_CAST(Vec8thUint, delimiter_bytes ^ (delimiter_bytes << 1));
        cur_transition_ct = PopcountVec8thUint(cur_transitions);
      } else {
        cur_transitions &= cur_transitions - 1;
      }
      byte_offset_in_vec = ctzu32(cur_transitions);
      const char* token_end = &(DowncastKToC(str_viter)[byte_offset_in_vec]);
      token_slens[cur_col_type] = token_end - token_start;
      ++relevant_col_idx;
      transition_bits_to_skip_p1 = 2 * col_skips[relevant_col_idx];
    }
    if (terminating_bytes) {
      return nullptr;
    }
    transition_bits_to_skip_p1 -= cur_transition_ct;
    prev_delim_highbit = transition_bits_to_skip_p1 & 1;
    ++str_viter;
    cur_vvec = *str_viter;
    tabspace_vvec = (cur_vvec == vvec_all_tab) | (cur_vvec == vvec_all_space);
    post31_vvec = vecuc_adds(cur_vvec, vvec_all96);
    delimiter_bytes = vecuc_movemask(tabspace_vvec);
    terminating_bytes = ~(vecuc_movemask(post31_vvec) | delimiter_bytes);
  }
}
#endif  // USE_AVX2

CXXCONST_CP NextCsvMult(const char* str_iter, uint32_t ct) {
  if (!str_iter) {
    return nullptr;
  }
  // assumes initial spaces in current token have been skipped
  // ok if we're at the end of the token
  unsigned char ucc = *str_iter;
  assert(ucc != ' ');
  while (ucc >= ' ') {
    // avoid strchr to keep "ASCII code < 32 == newline" consistent
    // (tab handling is quirky right now--permitted at the beginning of a
    // token, but treated as newline later--but it should never appear so
    // no point in e.g. adding an extra parameter to FirstNonTspace();
    // just need to make sure the quirky behavior is consistent.)
    if (ucc != ',') {
      ucc = *(++str_iter);
      continue;
    }
    do {
      ucc = *(++str_iter);
    } while ((ucc == ' ') || (ucc == '\t'));
    if (!(--ct)) {
      return S_CAST(CXXCONST_CP, str_iter);
    }
  }
  return nullptr;
}

#ifdef USE_AVX2
const char* CsvLexK(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens) {
  // Temporary; optimize later.
  for (uint32_t relevant_col_idx = 0; relevant_col_idx != relevant_col_ct; ++relevant_col_idx) {
    const uint32_t cur_col_type = col_types[relevant_col_idx];
    str_iter = NextCsvMult(str_iter, col_skips[relevant_col_idx]);
    if (!str_iter) {
      return nullptr;
    }
    token_ptrs[cur_col_type] = str_iter;
    const char* token_end = CsvFieldEnd(str_iter);
    token_slens[cur_col_type] = token_end - str_iter;
    str_iter = token_end;
  }
  return str_iter;
}
#endif

uint32_t CountTokens(const char* str_iter) {
  uint32_t token_ct = 0;
  str_iter = FirstNonTspace(str_iter);
  while (!IsEolnKns(*str_iter)) {
    ++token_ct;
    str_iter = FirstNonTspace(CurTokenEnd(str_iter));
  }
  return token_ct;
}

/*
uint32_t CommaOrSpaceCountTokens(const char* str_iter, uint32_t comma_delim) {
  if (comma_delim) {
    // assumes nonempty line (treats trailing empty string as a token).
    uint32_t token_ct = 1;
    unsigned char ucc = (unsigned char)(*str_iter++);
    while (1) {
      if (ucc < 32) {
        return token_ct;
      }
      if (ucc == ',') {
        // spelled out due to const qualifier
        do {
          ucc = (unsigned char)(*str_iter++);
        } while ((ucc == ' ') || (ucc == '\t'));
        token_ct++;
        continue;
      }
      ucc = (unsigned char)(*str_iter++);
    }
  }
  return CountTokens(str_iter);
}
*/

uint32_t CountAndMeasureMultistr(const char* multistr, uintptr_t* max_blen_ptr) {
  // Straightforward to accelerate this with movemask/shrn4 if it ever matters.
  // (Doesn't matter now since it's only used for command-line processing.)
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

// number-to-string encoders
// u32toa and i64toa moved to plink2_base

char* i32toa(int32_t ii, char* start) {
  uint32_t uii = ii;
  if (ii < 0) {
    // -INT_MIN is undefined, but negating the unsigned int equivalent works
    *start++ = '-';
    uii = -uii;
  }
  return u32toa(uii, start);
}

char* u32toa_z5(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  *start++ = '0' + quotient;
  return uitoa_z4(uii - 10000 * quotient, start);
}


char* u32toa_trunc4(uint32_t uii, char* start) {
  uint32_t quotient = uii / 100;
  memcpy_k(start, &(kDigitPair[quotient]), 2);
  uii -= 100 * quotient;
  if (uii) {
    start += 2;
    memcpy_k(start, &(kDigitPair[uii]), 2);
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  memcpy_k(start, &(kDigitPair[quotient]), 2);
  uii -= 10000 * quotient;
  if (uii) {
    quotient = uii / 100;
    start += 2;
    memcpy_k(start, &(kDigitPair[quotient]), 2);
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
      memcpy_k(start, &(kDigitPair[uii]), 2);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc8(uint32_t uii, char* start) {
  uint32_t quotient = uii / 1000000;
  memcpy_k(start, &(kDigitPair[quotient]), 2);
  uii -= 1000000 * quotient;
  if (uii) {
    quotient = uii / 10000;
    start += 2;
    memcpy_k(start, &(kDigitPair[quotient]), 2);
    uii -= 10000 * quotient;
    if (uii) {
      quotient = uii / 100;
      start += 2;
      memcpy_k(start, &(kDigitPair[quotient]), 2);
      uii -= 100 * quotient;
      if (uii) {
        start += 2;
        memcpy_k(start, &(kDigitPair[uii]), 2);
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
  memcpy_k(start, &(kDigitPair[quotient]), 2);
  remainder -= 1000 * quotient;
  if (remainder) {
    quotient = remainder / 10;
    start += 2;
    memcpy_k(start, &(kDigitPair[quotient]), 2);
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
  memcpy_k(start, &(kDigitPair[quotient]), 2);
  remainder -= 100000 * quotient;
  if (remainder) {
    quotient = remainder / 1000;
    start += 2;
    memcpy_k(start, &(kDigitPair[quotient]), 2);
    remainder -= 1000 * quotient;
    if (remainder) {
      quotient = remainder / 10;
      start += 2;
      memcpy_k(start, &(kDigitPair[quotient]), 2);
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
static const STD_ARRAY_DECL(double, 2, kBankerRound6) = STD_ARRAY_INIT_START() {0.4999995, 0.5000005} STD_ARRAY_INIT_END();
static const STD_ARRAY_DECL(double, 2, kBankerRound8) = STD_ARRAY_INIT_START() {0.499999995, 0.500000005} STD_ARRAY_INIT_END();

static inline uint32_t BankerRoundD(double dxx, STD_ARRAY_KREF(double, 2) banker_round) {
  uint32_t result = S_CAST(int32_t, dxx);
  return result + S_CAST(int32_t, (dxx - u31tod(result)) + banker_round[result & 1]);
}

// These are separate functions so the compiler can optimize the integer
// divisions.
static inline void BankerRoundD1(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10;
  *remainderp = remainder - (*quotientp) * 10;
}

static inline void BankerRoundD2(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100;
  *remainderp = remainder - (*quotientp) * 100;
}

static inline void BankerRoundD3(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000;
  *remainderp = remainder - (*quotientp) * 1000;
}

static inline void BankerRoundD4(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10000;
  *remainderp = remainder - (*quotientp) * 10000;
}

static inline void BankerRoundD5(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100000;
  *remainderp = remainder - (*quotientp) * 100000;
}

static inline void BankerRoundD6(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000000;
  *remainderp = remainder - (*quotientp) * 1000000;
}

static inline void BankerRoundD7(double dxx, STD_ARRAY_KREF(double, 2) banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000000;
  uint32_t remainder = S_CAST(int32_t, dxx);
  remainder += S_CAST(int32_t, (dxx - u31tod(remainder)) + banker_round[remainder & 1]);
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
      BankerRoundD5(dxx, kBankerRound8, &quotient, &remainder);
      return qrtoa_1p5(quotient, remainder, start);
    }
    BankerRoundD4(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy_k(start, &(kDigitPair[quotient]), 2);
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so6_pretail:
      memcpy_k(start, &(kDigitPair[remainder]), 2);
    }
  dtoa_so6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  }
  if (dxx < 9999.9949999999) {
    if (dxx < 999.99949999999) {
      BankerRoundD3(dxx, kBankerRound8, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya_k(start, &(kDigitPair[quotient]), 2);
      if (!remainder) {
        return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy_k(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto dtoa_so6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    BankerRoundD2(dxx, kBankerRound8, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so6_pretail;
  }
  if (dxx >= 99999.949999999) {
    return u32toa_z6(BankerRoundD(dxx, kBankerRound8), start);
  }
  BankerRoundD1(dxx, kBankerRound8, &uii, &remainder);
  quotient = uii / 10000;
  *start = '0' + quotient;
  uii -= 10000 * quotient;
  quotient = uii / 100;
  start = memcpya_k(&(start[1]), &(kDigitPair[quotient]), 2);
  uii = uii - 100 * quotient;
  start = memcpya_k(start, &(kDigitPair[uii]), 2);
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
      BankerRoundD7(dxx, kBankerRound6, &quotient, &remainder);
      return qrtoa_1p7(quotient, remainder, start);
    }
    BankerRoundD6(dxx, kBankerRound6, &quotient, &remainder);
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 10000;
    memcpy_k(start, &(kDigitPair[quotient]), 2);
    remainder -= 10000 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so8_pretail4:
      quotient = remainder / 100;
      memcpy_k(start, &(kDigitPair[quotient]), 2);
      remainder -= 100 * quotient;
      if (remainder) {
        start += 2;
      dtoa_so8_pretail2:
        memcpy_k(start, &(kDigitPair[remainder]), 2);
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
      BankerRoundD5(dxx, kBankerRound6, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya_k(start, &(kDigitPair[quotient]), 2);
      if (!remainder) {
        return start;
      }
      *start++ = '.';
      quotient = remainder / 1000;
      memcpy_k(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 1000;
      if (!remainder) {
        goto dtoa_so8_tail;
      }
      start += 2;
    dtoa_so8_pretail3:
      quotient = remainder / 10;
      memcpy_k(start, &(kDigitPair[quotient]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto dtoa_so8_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    BankerRoundD4(dxx, kBankerRound6, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so8_pretail4;
  }
  if (dxx < 999999.99499999) {
    if (dxx < 99999.999499999) {
      BankerRoundD3(dxx, kBankerRound6, &uii, &remainder);
      quotient = uii / 10000;
      *start = '0' + quotient;
      uii -= 10000 * quotient;
      quotient = uii / 100;
      start = memcpya_k(&(start[1]), &(kDigitPair[quotient]), 2);
      uii -= 100 * quotient;
      start = memcpya_k(start, &(kDigitPair[uii]), 2);
      if (!remainder) {
        return start;
      }
      *start++ = '.';
      goto dtoa_so8_pretail3;
    }
    BankerRoundD2(dxx, kBankerRound6, &uii, &remainder);
    quotient = uii / 10000;
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
    uii -= 100 * quotient;
    start = memcpya_k(start, &(kDigitPair[uii]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so8_pretail2;
  }
  if (dxx >= 9999999.9499999) {
    return uitoa_z8(BankerRoundD(dxx, kBankerRound6), start);
  }
  BankerRoundD1(dxx, kBankerRound6, &uii, &remainder);
  quotient = uii / 1000000;
  *start = '0' + quotient;
  uii -= 1000000 * quotient;
  quotient = uii / 10000;
  start = memcpya_k(&(start[1]), &(kDigitPair[quotient]), 2);
  uii -= 10000 * quotient;
  quotient = uii / 100;
  start = memcpya_k(start, &(kDigitPair[quotient]), 2);
  uii -= 100 * quotient;
  start = memcpya_k(start, &(kDigitPair[uii]), 2);
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
    return strcpya_k(start, "nan");
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
    BankerRoundD5(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya_k(qrtoa_1p5(quotient, remainder, start), "e-", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya_k(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 999999.49999999) {
    // 6 sig fig exponential notation, large
    if (dxx >= 9.9999949999999e15) {
      if (dxx >= 9.9999949999999e127) {
        if (dxx > DBL_MAX) {
          return strcpya_k(start, "inf");
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
    BankerRoundD5(dxx, kBankerRound8, &quotient, &remainder);
    start = memcpya_k(qrtoa_1p5(quotient, remainder, start), "e+", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya_k(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 0.99999949999999) {
    return dtoa_so6(dxx, start);
  }
  // 6 sig fig decimal, no less than ~0.0001
  start = memcpya_k(start, "0.", 2);
  if (dxx < 9.9999949999999e-3) {
    dxx *= 100;
    start = memcpya_k(start, "00", 2);
  }
  if (dxx < 9.9999949999999e-2) {
    dxx *= 10;
    *start++ = '0';
  }
  return uitoa_trunc6(BankerRoundD(dxx * 1000000, kBankerRound8), start);
}

char* dtoa_g_p8(double dxx, char* start) {
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return strcpya_k(start, "nan");
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
    BankerRoundD7(dxx, kBankerRound6, &quotient, &remainder);
    wpos = qrtoa_1p7(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      start = memcpya(start, wbuf, remainder);
      quotient = xp10 / 100;
      start = memcpyax_k(start, "e-", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      start = memcpya(start, wbuf, remainder);
      start = memcpya_k(start, "e-", 2);
    }
    return memcpya_k(start, &(kDigitPair[xp10]), 2);
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
    BankerRoundD7(dxx, kBankerRound6, &quotient, &remainder);
    wpos = qrtoa_1p7(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      start = memcpya(start, wbuf, remainder);
      quotient = xp10 / 100;
      start = memcpyax_k(start, "e+", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      start = memcpya(start, wbuf, remainder);
      start = memcpya_k(start, "e+", 2);
    }
    return memcpya_k(start, &(kDigitPair[xp10]), 2);
  }
  if (dxx >= 0.99999999499999) {
    wpos = dtoa_so8(dxx, wpos);
  } else {
    // 8 sig fig decimal, no less than ~0.0001
    wpos = memcpya_k(wpos, "0.", 2);
    if (dxx < 9.9999999499999e-3) {
      dxx *= 100;
      wpos = memcpya_k(wpos, "00", 2);
    }
    if (dxx < 9.9999999499999e-2) {
      dxx *= 10;
      *wpos++ = '0';
    }
    wpos = uitoa_trunc8(BankerRoundD(dxx * 100000000, kBankerRound6), wpos);
  }
  remainder = wpos - wbuf;
  return memcpya(start, wbuf, remainder);
}

// "prob" means that the number is guaranteed to be in [0, 1].
// no leading space is printed.  trailing zeroes (/decimal point) are erased
//   iff there is equality to ~13 decimal places.
char* dtoa_f_probp6_spaced(double dxx, char* start) {
  double dxx_10_6 = dxx * 1000000;
  const uint32_t dec_digits = BankerRoundD(dxx_10_6, kBankerRound8);
  *start++ = '0' + (dec_digits == 1000000);
  *start++ = '.';
  start = u32toa_z6(dec_digits, start);
  if (fabs(dxx_10_6 - u31tod(dec_digits)) >= 0.00000005) {
    return start;
  }
  TrailingZeroesToSpaces(start);
  return start;
}

char* dtoa_f_probp6_clipped(double dxx, char* start) {
  double dxx_10_6 = dxx * 1000000;
  const uint32_t dec_digits = BankerRoundD(dxx_10_6, kBankerRound8);
  *start++ = '0' + (dec_digits == 1000000);
  *start++ = '.';
  start = u32toa_z6(dec_digits, start);
  if (fabs(dxx_10_6 - u31tod(dec_digits)) >= 0.00000005) {
    return start;
  }
  return ClipTrailingZeroes(start);
}

/*
char* dtoa_f_p5_clipped(double dxx, char* start) {
  if (dxx != dxx) {
    return strcpya_k(start, "nan");
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
    start = u32toa(quotient, start);
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
    start = u32toa(quotient, start);
    return rtoa_p5(remainder, start);
  }
#endif
  if (dxx == INFINITY) {
    return strcpya_k(start, "inf");
  }
  // just punt larger numbers to glibc for now, this isn't a bottleneck
  start += sprintf(start, "%.5f", dxx);
  // .5f doesn't strip trailing zeroes, do that manually
  for (uint32_t uii = 0; uii != 5; ++uii) {
    if (start[-1] != '0') {
      return start;
    }
    --start;
  }
  return &(start[-1]);  // strip the decimal point
}
*/

// todo: benchmark the exponential-notation part of this vs. dtoa_g(); maybe
// dtoa_g() should actually call this (or at least the exponential-notation
// part, put into its own function) for the most extreme values?
char* lntoa_g(double ln_val, char* start) {
  // log(999999.49999999)
  if (ln_val < 13.81551005796414) {
    // log(9.9999949999999e-5)
    if (ln_val > -9.210340871976317) {
      // No exponential notation.

      // log(0.99999949999999)
      if (ln_val > -5.000001349509205e-7) {
        // may as well fast-path x=1; since the most common use-case for this
        // function is p-value printing, x=1 should happen a lot more than x>1.
        // log(1.0000050000001)
        if (ln_val < 4.999987599993995e-6) {
          *start++ = '1';
          return start;
        }
        return dtoa_so6(exp(ln_val), start);
      }
      double dxx = exp(ln_val);
      // 6 sig fig decimal, no less than ~0.0001
      start = memcpya_k(start, "0.", 2);
      if (dxx < 9.9999949999999e-3) {
        dxx *= 100;
        start = memcpya_k(start, "00", 2);
      }
      if (dxx < 9.9999949999999e-2) {
        dxx *= 10;
        *start++ = '0';
      }
      return uitoa_trunc6(BankerRoundD(dxx * 1000000, kBankerRound8), start);
    }
    // if exponent is in danger of overflowing int32, just print '0'
    if (ln_val < 0x7ffffffb * (-kLn10)) {
      *start++ = '0';
      return start;
    }
  } else {
    // if exponent is in danger of overflowing int32, just print 'inf'
    if (ln_val > 0x7ffffffb * kLn10) {
      return strcpya_k(start, "inf");
    }
  }
  int32_t xp10 = S_CAST(int32_t, (5.000001349509205e-7 + ln_val) * kRecipLn10);
  double mantissa = exp(xp10 * (-kLn10) + ln_val);
  // mantissa will usually be in [.9999995, 9.999995], but |ln_val| can be
  // larger than 2^32, and floating point errors in either direction are
  // definitely possible (<20 bits of precision).
  if (mantissa < 0.99999949999999) {
    mantissa *= 10;
    xp10 -= 1;
  } else if (mantissa > 9.9999949999999) {
    mantissa *= 0.1;
    xp10 += 1;
  }
  uint32_t quotient;
  uint32_t remainder;
  BankerRoundD5(mantissa, kBankerRound8, &quotient, &remainder);
  start = qrtoa_1p5(quotient, remainder, start);
  if (xp10 < 0) {
    start = memcpya_k(start, "e-", 2);
    if (xp10 > -10) {
      *start++ = '0';
    }
    return u32toa(-xp10, start);
  }
  start = memcpya_k(start, "e+", 2);
  if (xp10 < 10) {
    *start++ = '0';
  }
  return u32toa(xp10, start);
}

// Previously had specialized float-printing functions, but upon reflection it
// makes sense to just promote to double like printf does.


CXXCONST_CP FindSortedStrboxDuplicate(const char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_blen) {
  --id_ct;
  for (uintptr_t id_idx = 0; id_idx != id_ct; ++id_idx) {
    if (strequal_overread(&(sorted_ids[id_idx * max_id_blen]), &(sorted_ids[(id_idx + 1) * max_id_blen]))) {
      return S_CAST(CXXCONST_CP, &(sorted_ids[id_idx * max_id_blen]));
    }
  }
  return nullptr;
}


int32_t bsearch_strbox(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx) {
  // does not assume null-terminated idbuf, or nonempty array.
  if (cur_id_slen >= max_id_blen) {
    return -1;
  }
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    const int32_t ii = Memcmp(idbuf, &(sorted_strbox[mid_idx * max_id_blen]), cur_id_slen);
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if ((ii < 0) || sorted_strbox[mid_idx * max_id_blen + cur_id_slen]) {
      end_idx = mid_idx;
    } else {
      return S_CAST(uint32_t, mid_idx);
    }
  }
  return -1;
}

int32_t bsearch_strbox_natural(const char* idbuf, const char* sorted_strbox, uintptr_t max_id_blen, uintptr_t end_idx) {
  // unlike bsearch_strbox(), caller is responsible for slen >= max_id_blen
  // check if appropriate here
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    const int32_t ii = strcmp_natural(idbuf, &(sorted_strbox[mid_idx * max_id_blen]));
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if (ii < 0) {
      end_idx = mid_idx;
    } else {
      return S_CAST(uint32_t, mid_idx);
    }
  }
  return -1;
}

int32_t bsearch_strptr_overread(const char* idbuf, const char* const* sorted_strptrs, uintptr_t end_idx) {
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    const int32_t ii = strcmp_overread(idbuf, sorted_strptrs[mid_idx]);
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if (ii < 0) {
      end_idx = mid_idx;
    } else {
      return S_CAST(uint32_t, mid_idx);
    }
  }
  return -1;
}

int32_t bsearch_strptr_natural(const char* idbuf, const char* const* sorted_strptrs, uintptr_t end_idx) {
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    const int32_t ii = strcmp_natural(idbuf, sorted_strptrs[mid_idx]);
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if (ii < 0) {
      end_idx = mid_idx;
    } else {
      return S_CAST(uint32_t, mid_idx);
    }
  }
  return -1;
}

uintptr_t bsearch_strbox_lb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx) {
  // returns number of elements in sorted_strbox[] less than idbuf.
  if (cur_id_slen > max_id_blen) {
    cur_id_slen = max_id_blen;
  }
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    if (Memcmp(idbuf, &(sorted_strbox[mid_idx * max_id_blen]), cur_id_slen) > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}

uintptr_t bsearch_strbox_lb_natural(const char* idbuf, const char* nsorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx) {
  if (cur_id_slen > max_id_blen) {
    cur_id_slen = max_id_blen;
  }
  uintptr_t start_idx = 0;
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    if (strcmp_natural(idbuf, &(nsorted_strbox[mid_idx * max_id_blen])) > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}

uintptr_t ExpsearchStrLb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx, uintptr_t cur_idx) {
  uintptr_t next_incr = 1;
  uintptr_t start_idx = cur_idx;
  while (cur_idx < end_idx) {
    if (Memcmp(idbuf, &(sorted_strbox[cur_idx * max_id_blen]), cur_id_slen) <= 0) {
      end_idx = cur_idx;
      break;
    }
    start_idx = cur_idx + 1;
    cur_idx += next_incr;
    next_incr *= 2;
  }
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    if (Memcmp(idbuf, &(sorted_strbox[mid_idx * max_id_blen]), cur_id_slen) > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}

uintptr_t ExpsearchNsortStrLb(const char* idbuf, const char* nsorted_strbox, uintptr_t max_id_blen, uintptr_t end_idx, uintptr_t cur_idx) {
  uintptr_t next_incr = 1;
  uintptr_t start_idx = cur_idx;
  while (cur_idx < end_idx) {
    if (strcmp_natural(idbuf, &(nsorted_strbox[cur_idx * max_id_blen])) <= 0) {
      end_idx = cur_idx;
      break;
    }
    start_idx = cur_idx + 1;
    cur_idx += next_incr;
    next_incr *= 2;
  }
  while (start_idx < end_idx) {
    const uintptr_t mid_idx = (start_idx + end_idx) / 2;
    if (strcmp_natural(idbuf, &(nsorted_strbox[mid_idx * max_id_blen])) > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}

uint32_t IsInfStr(const char* ss, uint32_t slen, uint32_t* is_neg_ptr) {
  const char first_char = *ss;
  if (first_char == '-') {
    *is_neg_ptr = 1;
    ++ss;
    --slen;
  } else if (first_char == '+') {
    ++ss;
    --slen;
  }
  if (slen == 3) {
    uint32_t four_chars;
    memcpy(&four_chars, ss, 4);  // OVERREAD
    // assumes little-endian
    return ((four_chars & 0xdfdfdf) == 0x464e49);
  }
  if (slen != 8) {
    return 0;
  }
  uint64_t eight_chars;
  memcpy(&eight_chars, ss, 8);
  return ((eight_chars & 0xdfdfdfdfdfdfdfdfLLU) == 0x5954494e49464e49LLU);
}


#ifdef __LP64__
CXXCONST_CP Memrchr(const char* str_start, char needle, uintptr_t slen) {
  const VecI8 vvec_all_needle = veci8_set1(needle);
  const uintptr_t str_end_addr = R_CAST(uintptr_t, str_start) + slen;
  const uint32_t trailing_byte_ct = str_end_addr % kBytesPerVec;
  const VecI8* str_rev_viter = R_CAST(const VecI8*, RoundDownPow2(str_end_addr, kBytesPerVec));
  if (trailing_byte_ct) {
    // This is a GNU vector extension parallel-equality check, which gets
    // translated to e.g. _mm256_cmpeq_epi8().
    // As long as we're performing aligned reads, it's safe to read bytes
    // beyond str_end as long as they're in the same vector; we only risk
    // violating process read permissions if we cross a page boundary.
    // (For this reason, I don't bother with AVX unaligned reads.)
    const VecI8 match_vvec = (*str_rev_viter == vvec_all_needle);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    uint32_t matching_bytes = veci8_movemask(match_vvec);
    matching_bytes &= (1U << (trailing_byte_ct % kBytesPerVec)) - 1;
    if (str_start > DowncastKToC(str_rev_viter)) {
      const uint32_t leading_byte_ct = R_CAST(uintptr_t, str_start) % kBytesPerVec;
      matching_bytes &= -(1U << leading_byte_ct);
      // Special-case this, since main_loop_iter_ct below would underflow.
      if (!matching_bytes) {
        return nullptr;
      }
    }
    if (matching_bytes) {
      const uint32_t byte_offset_in_vec = bsru32(matching_bytes);
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  else
    uint64_t matching_nybbles = arm_shrn4_i8(match_vvec);
    matching_nybbles &= (1LLU << (4 * (trailing_byte_ct % kBytesPerVec))) - 1;
    if (str_start > DowncastKToC(str_rev_viter)) {
      const uint32_t leading_byte_ct = R_CAST(uintptr_t, str_start) % kBytesPerVec;
      matching_nybbles &= -(1LLU << (4 * leading_byte_ct));
      // Special-case this, since main_loop_iter_ct below would underflow.
      if (!matching_nybbles) {
        return nullptr;
      }
    }
    if (matching_nybbles) {
      const uint32_t byte_offset_in_vec = bsrw(matching_nybbles) / 4;
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  endif
  }
  const uintptr_t main_loop_iter_ct = (R_CAST(uintptr_t, str_rev_viter) - R_CAST(uintptr_t, str_start)) / (2 * kBytesPerVec);
  for (uintptr_t ulii = 0; ulii != main_loop_iter_ct; ++ulii) {
    // For long lines, looping over two vectors at a time is most efficient on
    // my Mac (also tried 1 and 4).
    --str_rev_viter;
    const VecI8 match_vvec1 = (*str_rev_viter == vvec_all_needle);
    --str_rev_viter;
    const VecI8 match_vvec0 = (*str_rev_viter == vvec_all_needle);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const uint32_t matching_bytes = veci8_movemask(match_vvec1 | match_vvec0);
    if (matching_bytes) {
      const uint32_t matching_bytes1 = veci8_movemask(match_vvec1);
      if (matching_bytes1) {
        const uint32_t byte_offset_in_vec = bsru32(matching_bytes1);
        return &(DowncastToXC(&(str_rev_viter[1]))[byte_offset_in_vec]);
      }
      const uint32_t byte_offset_in_vec = bsru32(matching_bytes);
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  else
    const uint64_t matching_nybbles = arm_shrn4_i8(match_vvec1 | match_vvec0);
    if (matching_nybbles) {
      const uint64_t matching_nybbles1 = arm_shrn4_i8(match_vvec1);
      if (matching_nybbles1) {
        const uint32_t byte_offset_in_vec = bsrw(matching_nybbles1) / 4;
        return &(DowncastToXC(&(str_rev_viter[1]))[byte_offset_in_vec]);
      }
      const uint32_t byte_offset_in_vec = bsrw(matching_nybbles) / 4;
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  endif
  }
  while (1) {
    uintptr_t remaining_byte_ct_underflow = R_CAST(uintptr_t, str_rev_viter) - R_CAST(uintptr_t, str_start);
    if (S_CAST(intptr_t, remaining_byte_ct_underflow) <= 0) {
      return nullptr;
    }
    --str_rev_viter;
    const VecI8 match_vvec = (*str_rev_viter == vvec_all_needle);
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const uint32_t matching_bytes = veci8_movemask(match_vvec);
    if (matching_bytes) {
      const uint32_t byte_offset_in_vec = bsru32(matching_bytes);
      if (byte_offset_in_vec + remaining_byte_ct_underflow < kBytesPerVec) {
        return nullptr;
      }
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  else
    const uint64_t matching_nybbles = arm_shrn4_i8(match_vvec);
    if (matching_nybbles) {
      const uint32_t byte_offset_in_vec = bsrw(matching_nybbles) / 4;
      if (byte_offset_in_vec + remaining_byte_ct_underflow < kBytesPerVec) {
        return nullptr;
      }
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  endif
  }
}

CXXCONST_CP LastSpaceOrEoln(const char* str_start, uintptr_t slen) {
  // See TokenLexK0().
  const VecUc vvec_all95 = vecuc_set1(95);
#  ifdef SIMDE_ARM_NEON_A32V8_NATIVE
  const VecUc vvec_all94 = vecuc_set1(94);
#  endif
  const uintptr_t str_end_addr = R_CAST(uintptr_t, str_start) + slen;
  const uint32_t trailing_byte_ct = str_end_addr % kBytesPerVec;
  const VecUc* str_rev_viter = R_CAST(const VecUc*, RoundDownPow2(str_end_addr, kBytesPerVec));
  if (trailing_byte_ct) {
    // As long as we're performing aligned reads, it's safe to read bytes
    // beyond str_end as long as they're in the same vector; we only risk
    // violating process read permissions if we cross a page boundary.
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const VecUc postspace_vvec = vecuc_adds(*str_rev_viter, vvec_all95);
    uint32_t nontoken_bytes = S_CAST(Vec8thUint, ~vecuc_movemask(postspace_vvec));
    nontoken_bytes &= (1U << (trailing_byte_ct % kBytesPerVec)) - 1;
    if (str_start > DowncastKToC(str_rev_viter)) {
      const uint32_t leading_byte_ct = R_CAST(uintptr_t, str_start) % kBytesPerVec;
      nontoken_bytes &= -(1U << leading_byte_ct);
      // Special-case this, since main_loop_iter_ct below would underflow.
      if (!nontoken_bytes) {
        return nullptr;
      }
    }
    if (nontoken_bytes) {
      const uint32_t byte_offset_in_vec = bsru32(nontoken_bytes);
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  else
    // bugfix (14 Apr 2025): postspace_vvec was not of vec0255 form.
    const VecUc pre33_vvec = vecuc_signed_cmpgt(vecuc_add(*str_rev_viter, vvec_all95), vvec_all94);
    uint64_t nontoken_nybbles = arm_shrn4_uc(pre33_vvec);
    nontoken_nybbles &= (1LLU << (4 * (trailing_byte_ct % kBytesPerVec))) - 1;
    if (str_start > DowncastKToC(str_rev_viter)) {
      const uint32_t leading_byte_ct = R_CAST(uintptr_t, str_start) % kBytesPerVec;
      nontoken_nybbles &= -(1LLU << (4 * leading_byte_ct));
      // Special-case this, since main_loop_iter_ct below would underflow.
      if (!nontoken_nybbles) {
        return nullptr;
      }
    }
    if (nontoken_nybbles) {
      const uint32_t byte_offset_in_vec = bsrw(nontoken_nybbles) / 4;
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  endif
  }
  const uintptr_t main_loop_iter_ct = (R_CAST(uintptr_t, str_rev_viter) - R_CAST(uintptr_t, str_start)) / (2 * kBytesPerVec);
  for (uintptr_t ulii = 0; ulii != main_loop_iter_ct; ++ulii) {
    --str_rev_viter;
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const VecUc postspace_vvec1 = vecuc_adds(*str_rev_viter, vvec_all95);
    --str_rev_viter;
    const VecUc postspace_vvec0 = vecuc_adds(*str_rev_viter, vvec_all95);
    const uint32_t nontoken_bytes = S_CAST(Vec8thUint, ~vecuc_movemask(postspace_vvec1 & postspace_vvec0));
    if (nontoken_bytes) {
      const uint32_t nontoken_bytes1 = S_CAST(Vec8thUint, ~vecuc_movemask(postspace_vvec1));
      if (nontoken_bytes1) {
        const uint32_t byte_offset_in_vec = bsru32(nontoken_bytes1);
        return &(DowncastToXC(&(str_rev_viter[1]))[byte_offset_in_vec]);
      }
      const uint32_t byte_offset_in_vec = bsru32(nontoken_bytes);
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  else
    const VecUc pre33_vvec1 = vecuc_signed_cmpgt(vecuc_add(*str_rev_viter, vvec_all95), vvec_all94);
    --str_rev_viter;
    const VecUc pre33_vvec0 = vecuc_signed_cmpgt(vecuc_add(*str_rev_viter, vvec_all95), vvec_all94);
    const uint64_t nontoken_nybbles = arm_shrn4_uc(pre33_vvec1 | pre33_vvec0);
    if (nontoken_nybbles) {
      const uint32_t nontoken_nybbles1 = arm_shrn4_uc(pre33_vvec1);
      if (nontoken_nybbles1) {
        const uint32_t byte_offset_in_vec = bsrw(nontoken_nybbles1) / 4;
        return &(DowncastToXC(&(str_rev_viter[1]))[byte_offset_in_vec]);
      }
      const uint32_t byte_offset_in_vec = bsrw(nontoken_nybbles) / 4;
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  endif
  }
  while (1) {
    uintptr_t remaining_byte_ct_underflow = R_CAST(uintptr_t, str_rev_viter) - R_CAST(uintptr_t, str_start);
    if (S_CAST(intptr_t, remaining_byte_ct_underflow) <= 0) {
      return nullptr;
    }
    --str_rev_viter;
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const VecUc postspace_vvec = vecuc_adds(*str_rev_viter, vvec_all95);
    const uint32_t nontoken_bytes = S_CAST(Vec8thUint, ~vecuc_movemask(postspace_vvec));
    if (nontoken_bytes) {
      const uint32_t byte_offset_in_vec = bsru32(nontoken_bytes);
      if (byte_offset_in_vec + remaining_byte_ct_underflow < kBytesPerVec) {
        return nullptr;
      }
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  else
    const VecUc pre33_vvec = vecuc_signed_cmpgt(vecuc_add(*str_rev_viter, vvec_all95), vvec_all94);
    const uint64_t nontoken_nybbles = arm_shrn4_uc(pre33_vvec);
    if (nontoken_nybbles) {
      const uint32_t byte_offset_in_vec = bsrw(nontoken_nybbles) / 4;
      if (byte_offset_in_vec + remaining_byte_ct_underflow < kBytesPerVec) {
        return nullptr;
      }
      return &(DowncastToXC(str_rev_viter)[byte_offset_in_vec]);
    }
#  endif
  }
}
#endif  // __LP64__

/*
void ReplaceAllInstances(char old_char, char new_char, uint32_t slen, char* dst) {
  // probable todo: add SIMD version
  while (1) {
    char* strchr_result = memchr(dst, old_char, slen);
    if (!strchr_result) {
      return;
    }
    *strchr_result++ = new_char;
    slen -= strchr_result - dst;
    dst = strchr_result;
  }
}
*/

void TabsToSpaces(char* ss_iter) {
  while (1) {
    ss_iter = strchr(ss_iter, '\t');
    if (!ss_iter) {
      return;
    }
    *ss_iter++ = ' ';
  }
}

BoolErr ReplaceCharAdvChecked(char old_char, char new_char, char** str_ptr) {
  char* token_end;
  for (char* str_iter = *str_ptr; ; str_iter = &(token_end[1])) {
    token_end = strchrnul2(str_iter, old_char, new_char);
    if ((*token_end) != old_char) {
      if (likely(!(*token_end))) {
        *str_ptr = token_end;
        return 0;
      }
      return 1;
    }
    *token_end = new_char;
  }
}

#ifdef __cplusplus
}  // namespace plink2
#endif
