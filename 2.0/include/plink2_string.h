#ifndef __PLINK2_STRING_H__
#define __PLINK2_STRING_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Standalone string-printing and parsing functions which neither make
// permanent memory allocations nor use g_bigstack for temporary allocations.

#include "plink2_base.h"

#include <math.h>  // fabs(), isfinite()
#include <stddef.h>  // offsetof()

#ifdef __cplusplus
#  include <algorithm>
#  if __cplusplus >= 201902L
#    include <execution>
#  endif
#  ifdef _WIN32
    // Windows C++11 <algorithm> resets these values :(
#    undef PRIu64
#    undef PRId64
#    define PRIu64 "I64u"
#    define PRId64 "I64d"
#    undef PRIuPTR
#    undef PRIdPTR
#    ifdef __LP64__
#      define PRIuPTR PRIu64
#      define PRIdPTR PRId64
#    else
#      if __cplusplus < 201103L
#        define PRIuPTR "lu"
#        define PRIdPTR "ld"
#      else
#        define PRIuPTR "u"
#        define PRIdPTR "d"
#      endif
#    endif
#  endif
#endif

#ifdef _WIN32
#  define EOLN_STR "\r\n"
#else
#  define EOLN_STR "\n"
#endif

// generic maximum line byte length, currently also used as a default I/O
// buffer size.  .ped/.vcf/etc. lines can of course be longer.
CONSTI32(kMaxMediumLine, 131072);

// apparently these aren't always defined in limits.h
#ifndef DBL_MAX
#  define DBL_MAX 1.7976931348623157e308
#endif
#ifndef FLT_MAX
#  define FLT_MAX S_CAST(float, 3.40282347e38)
#endif

// These are needed by plink2_stats.  May want to define these elsewhere, but
// they can't live in plink2_cmdline any more.
static const double kE = 2.7182818284590452;
static const double kPi = 3.1415926535897932;
static const double kSqrt2 = 1.4142135623730951;
static const double kRecipE = 0.36787944117144233;
static const double kLn2 = 0.6931471805599453;
static const double kRecip2m53 = 0.00000000000000011102230246251565404236316680908203125;

// floating point comparison-to-nonzero tolerance, currently 2^{-30}
static const double kEpsilon = 0.000000000931322574615478515625;
// less tolerant versions (2^{-35}, 2^{-44}) for some exact calculations
static const double kSmallEpsilon = 0.00000000000005684341886080801486968994140625;

// 2^{-21}, must be >= sqrt(kSmallEpsilon)
static const double kBigEpsilon = 0.000000476837158203125;

// 2^{-83} bias to give exact tests maximum ability to determine tiny p-values.
// (~2^{-53} is necessary to take advantage of denormalized small numbers, then
// allow tail sum to be up to 2^30.)
static const double kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;

#ifdef __cplusplus
#  define STD_SORT(ct, fallback_cmp, arr) std::sort(&((arr)[0]), (&((arr)[ct])))
#  if __cplusplus >= 201902L
// this should only be used for arrays of length >= variant_ct or sample_ct
// (sample_ct is cutting it close).
// macro should still be used in e.g. non-__cplusplus blocks, so that we have
// the option of falling back on a hand-coded parallel sort.
#    define STD_SORT_PAR_UNSEQ(ct, fallback_cmp, arr) std::sort(std::execution::par_unseq, &((arr)[0]), (&((arr)[ct])))
#  else
#    define STD_SORT_PAR_UNSEQ(ct, fallback_cmp, arr) std::sort(&((arr)[0]), (&((arr)[ct])))
#  endif
#else
#  define STD_SORT(ct, fallback_cmp, arr) qsort((arr), (ct), sizeof(*(arr)), (fallback_cmp))
#  define STD_SORT_PAR_UNSEQ(ct, fallback_cmp, arr) qsort((arr), (ct), sizeof(*(arr)), (fallback_cmp))
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

// A bunch of library functions operating on char*s don't modify the buffer
// themselves, but return a char* which should inherit the constness of the
// input parameter.  We want them to
//   1. check const-correctness when this is compiled as C++
//   2. still be valid C99
// and the method of achieving this should be minimally bug-prone.
//
// Current hack:
//   1. First declare the const char*-accepting version, with return type of
//      CXXCONST_CP.  Make all return statements include a C-style cast, even
//      when not strictly necessary (e.g. relaying the return value of another
//      function of this sort), unless nullptr is being returned.
//   2. Then declare the char*-accepting version in an immediately following
//      #ifdef __cplusplus block, which is simply a wrapper that uses
//      const_cast twice.
// This is kind of ugly, so I may not use this approach in non-library source
// code.
//
// There are related issues with const char** and const char* const*.
// Unfortunately, while C99 implicitly converts char* to const char*, it does
// not do char** -> const char* const*, so caller-side casts are required.
// (And char** -> const char** implicit conversion doesn't happen in C++ either
// for good reason, see http://c-faq.com/ansi/constmismatch.html .)  The 'good'
// news is that this removes the need for duplicate C++ function prototypes,
// but it's generally an even worse situation than the single-indirection case;
// these macro names are intentionally verbose to encourage assignment to the
// appropriate-qualified type as soon as possible.
#ifdef __cplusplus
#  define CXXCONST_CP const char*
#  define CXXCONST_VOIDP const void*
#  define TO_CONSTCPCONSTP(char_pp) (char_pp)
#else
#  define CXXCONST_CP char*
#  define CXXCONST_VOIDP void*
#  define TO_CONSTCPCONSTP(char_pp) ((const char* const*)(char_pp))
#endif

#ifdef _GNU_SOURCE
// There was some recent (2016) discussion on the gcc mailing list on strlen()
// vs. rawmemchr(., 0), where it was claimed that rawmemchr(., 0) could be
// compiled to something slower than &(.[strlen(.)]), rather than being at
// least as good.  However, this didn't happen when I tried to benchmark this,
// so I'll stick to the function that returns the right type (except when
// rawmemchr itself has to be emulated).
HEADER_INLINE CXXCONST_CP strnul(const char* str) {
  return S_CAST(CXXCONST_CP, rawmemchr(str, 0));
}

#  ifdef __cplusplus
HEADER_INLINE char* strnul(char* str) {
  return const_cast<char*>(strnul(const_cast<const char*>(str)));
}
#  endif

#else  // !_GNU_SOURCE

#  ifdef __LP64__
CXXCONST_VOIDP rawmemchr(const void* ss, int cc);

HEADER_INLINE CXXCONST_CP strnul(const char* str) {
  return S_CAST(CXXCONST_CP, rawmemchr(str, 0));
}
#  else  // !_GNU_SOURCE, !__LP64__
HEADER_INLINE CXXCONST_VOIDP rawmemchr(const void* ss, int cc) {
  return S_CAST(CXXCONST_VOIDP, memchr(ss, cc, 0x80000000U - kBytesPerVec));
}

HEADER_INLINE CXXCONST_CP strnul(const char* str) {
  return S_CAST(CXXCONST_CP, &(str[strlen(str)]));
}
#  endif

#  ifdef __cplusplus
HEADER_INLINE void* rawmemchr(void* ss, int cc) {
  return const_cast<void*>(rawmemchr(const_cast<const void*>(ss), cc));
}

HEADER_INLINE char* strnul(char* str) {
  return const_cast<char*>(strnul(const_cast<const char*>(str)));
}
#  endif

#endif  // !_GNU_SOURCE

// See also AdvToDelimOrEnd later below, which is an obvious alternative memchr
// interface.

// ReadLineStream emits lines which are *not* null-terminated, but are
// guaranteed to have trailing '\n's.
CXXCONST_VOIDP rawmemchr2(const void* ss, unsigned char ucc1, unsigned char ucc2);

HEADER_INLINE CXXCONST_CP strchrnul_n(const char* ss, unsigned char ucc1) {
  return S_CAST(CXXCONST_CP, rawmemchr2(ss, ucc1, '\n'));
}

CXXCONST_VOIDP rawmemchr3(const void* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3);

HEADER_INLINE CXXCONST_CP strchrnul2(const char* ss, unsigned char ucc1, unsigned char ucc2) {
  return S_CAST(CXXCONST_CP, rawmemchr3(ss, ucc1, ucc2, '\0'));
}

HEADER_INLINE CXXCONST_CP strchrnul2_n(const char* ss, unsigned char ucc1, unsigned char ucc2) {
  return S_CAST(CXXCONST_CP, rawmemchr3(ss, ucc1, ucc2, '\n'));
}

CXXCONST_CP strchrnul3(const char* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3);

#ifdef __cplusplus
HEADER_INLINE void* rawmemchr2(void* ss, unsigned char ucc1, unsigned char ucc2) {
  return const_cast<void*>(rawmemchr2(const_cast<const void*>(ss), ucc1, ucc2));
}

HEADER_INLINE char* strchrnul_n(char* ss, unsigned char ucc1) {
  return const_cast<char*>(strchrnul_n(const_cast<const char*>(ss), ucc1));
}

HEADER_INLINE void* rawmemchr3(void* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3) {
  return const_cast<void*>(rawmemchr3(const_cast<const void*>(ss), ucc1, ucc2, ucc3));
}

HEADER_INLINE char* strchrnul2(char* ss, unsigned char ucc1, unsigned char ucc2) {
  return const_cast<char*>(strchrnul2(const_cast<const char*>(ss), ucc1, ucc2));
}

HEADER_INLINE char* strchrnul2_n(char* ss, unsigned char ucc1, unsigned char ucc2) {
  return const_cast<char*>(strchrnul2_n(const_cast<const char*>(ss), ucc1, ucc2));
}

HEADER_INLINE char* strchrnul3(char* ss, unsigned char ucc1, unsigned char ucc2, unsigned char ucc3) {
  return const_cast<char*>(strchrnul3(const_cast<const char*>(ss), ucc1, ucc2, ucc3));
}
#endif

#ifndef _GNU_SOURCE
#  ifdef __LP64__
HEADER_INLINE CXXCONST_CP strchrnul(const char* str, int needle) {
  return S_CAST(CXXCONST_CP, rawmemchr2(str, 0, needle));
}
#  else
HEADER_INLINE CXXCONST_CP strchrnul(const char* str, int cc) {
  const char* strchr_result = strchr(str, cc);
  if (strchr_result) {
    return S_CAST(CXXCONST_CP, strchr_result);
  }
  return S_CAST(CXXCONST_CP, strnul(str));
}
#  endif

#  ifdef __cplusplus
HEADER_INLINE char* strchrnul(char* ss, int needle) {
  return const_cast<char*>(strchrnul(const_cast<const char*>(ss), needle));
}
#  endif
#endif

// These return 1 at eoln.
HEADER_INLINE uint32_t strchrnul_n_mov(unsigned char ucc1, const char** ss_ptr) {
  const char* ss_next = strchrnul_n(*ss_ptr, ucc1);
  *ss_ptr = ss_next;
  return (*ss_next != ucc1);
}

HEADER_INLINE uint32_t incr_strchrnul_n_mov(unsigned char ucc1, const char** ss_ptr) {
  const char* ss_next = strchrnul_n(&((*ss_ptr)[1]), ucc1);
  *ss_ptr = ss_next;
  return (*ss_next != ucc1);
}

// input for WordWrap should have no intermediate '\n's.  If suffix_len is 0,
// there should be a terminating \n.
void WordWrap(uint32_t suffix_len, char* strbuf);

uint32_t UintSlen(uint32_t num);

HEADER_INLINE uint32_t IntSlen(int32_t num) {
  // see abs_i32()
  const uint32_t neg_sign_bit = -(S_CAST(uint32_t, num) >> 31);
  return UintSlen((S_CAST(uint32_t, num) ^ neg_sign_bit) - neg_sign_bit) - neg_sign_bit;
}

// memcpya(), memseta() defined in plink2_base.h

HEADER_INLINE char* memcpyax(void* __restrict dst, const void* __restrict src, uint32_t ct, char extra_char) {
  memcpy(dst, src, ct);
  S_CAST(char*, dst)[ct] = extra_char;
  return &(S_CAST(char*, dst)[ct + 1]);
}

HEADER_INLINE void memcpyx(void* __restrict dst, const void* __restrict src, uint32_t ct, char extra_char) {
  memcpy(dst, src, ct);
  S_CAST(char*, dst)[ct] = extra_char;
}

HEADER_INLINE char* strcpya(char* __restrict dst, const void* __restrict src) {
  const uintptr_t slen = strlen(S_CAST(const char*, src));
  return memcpya(dst, src, slen);
}

HEADER_INLINE char* strcpyax(char* __restrict dst, const void* __restrict src, char extra_char) {
  const uintptr_t slen = strlen(S_CAST(const char*, src));
  memcpy(dst, src, slen);
  dst[slen] = extra_char;
  return &(dst[slen + 1]);
}

// MinGW support for stpcpy is a mess, so I'll use a capitalized name to route
// around the problem.
#if defined(_GNU_SOURCE) || defined(__APPLE__) || (_POSIX_C_SOURCE >= 200809L)
HEADER_INLINE char* Stpcpy(char* __restrict dst, const char* __restrict src) {
  return stpcpy(dst, src);
}
#else
HEADER_INLINE char* Stpcpy(char* __restrict dst, const char* __restrict src) {
  uintptr_t slen = strlen(src);
  memcpy(dst, src, slen + 1);
  return &(dst[slen]);
}
#endif

#if defined(__LP64__) && (__cplusplus >= 201103L)
template <uint32_t N> char* MemcpyaxK(void* __restrict dst, const void* __restrict src, char extra_char) {
  memcpyo_k(dst, src, N);
  S_CAST(char*, dst)[N] = extra_char;
  return &(S_CAST(char*, dst)[N + 1]);
}

#  define memcpyax_k(dst, src, ct, extra_char) plink2::MemcpyaxK<ct>(dst, src, extra_char)

template <uint32_t N> int32_t StrequalK(const char* s1, const char* k_s2, uint32_t s1_slen) {
  return (s1_slen == N) && memequal_k(s1, k_s2, N);
};

constexpr uint32_t CompileTimeSlen(const char* k_str) {
  return k_str[0]? (1 + CompileTimeSlen(&(k_str[1]))) : 0;
}

// can also use sizeof(k_s2) - 1, but that's less safe
#  define strequal_k(s1, k_s2, s1_slen) plink2::StrequalK<plink2::CompileTimeSlen(k_s2)>(s1, k_s2, s1_slen)

#  define strequal_k_unsafe(s1, k_s2) memequal_k(s1, k_s2, 1 + plink2::CompileTimeSlen(k_s2))

#  define strcpy_k(dst, src) plink2::MemcpyKImpl<plink2::CompileTimeSlen(src) + 1>::MemcpyK(dst, src);

#  define strcpya_k(dst, src) plink2::MemcpyaoK<plink2::CompileTimeSlen(src)>(dst, src);
#else  // !(defined(__LP64__) && (__cplusplus >= 201103L))
HEADER_INLINE char* memcpyax_k(void* __restrict dst, const void* __restrict src, uint32_t ct, char extra_char) {
  return memcpyax(dst, src, ct, extra_char);
}

HEADER_INLINE int32_t strequal_k(const char* s1, const char* k_s2, uint32_t s1_slen) {
  // any sane compiler should compute s2_slen at compile-time if k_s2 is a
  // constant string
  const uint32_t s2_slen = strlen(k_s2);
  return (s1_slen == s2_slen) && memequal(s1, k_s2, s2_slen);
}

// Can use this when it's always safe to read first (1 + strlen(k_s2)) bytes of
// s1.
HEADER_INLINE int32_t strequal_k_unsafe(const char* s1, const char* k_s2) {
  const uint32_t s2_blen = 1 + strlen(k_s2);
  return memequal(s1, k_s2, s2_blen);
}

HEADER_INLINE void strcpy_k(char* __restrict dst, const void* __restrict src) {
  strcpy(dst, S_CAST(const char*, src));
}

HEADER_INLINE char* strcpya_k(char* __restrict dst, const void* __restrict src) {
  return strcpya(dst, src);
}
#endif

#if defined(__cplusplus)
#  if __cplusplus >= 201103L
HEADER_INLINE bool isfinite_f(float fxx) {
  using namespace std;
  return isfinite(fxx);
}
#  else
#    ifdef isfinite
#      define isfinite_f isfinite
#    else
HEADER_INLINE bool isfinite_f(float fxx) {
  return (fxx == fxx) && (fxx != INFINITY) && (fxx != -INFINITY);
}
#    endif
#  endif
#else
#  define isfinite_f isfinite
#endif

HEADER_INLINE int32_t IsSpaceOrEoln(unsigned char ucc) {
  return (ucc <= 32);
}

// Assumes it's safe to read first 1 + strlen(s_const) bytes of s_read, i.e.
// this is ALWAYS 'unsafe'.
// Differs from strequal_k_unsafe() since strings are not considered equal when
// s_read[strlen(s_const)] isn't a token-ender.
HEADER_INLINE int32_t tokequal_k(const char* s_read, const char* s_const) {
  const uint32_t s_const_slen = strlen(s_const);
  return memequal(s_read, s_const, s_const_slen) && IsSpaceOrEoln(s_read[s_const_slen]);
}

// s_prefix must be strictly contained.
HEADER_INLINE int32_t StrStartsWith(const char* s_read, const char* s_prefix_const, uint32_t s_read_slen) {
  const uint32_t s_const_slen = strlen(s_prefix_const);
  return (s_read_slen > s_const_slen) && memequal(s_read, s_prefix_const, s_const_slen);
}

// permits s_read and s_prefix to be equal.
HEADER_INLINE int32_t StrStartsWith0(const char* s_read, const char* s_prefix_const, uint32_t s_read_slen) {
  const uint32_t s_const_slen = strlen(s_prefix_const);
  return (s_read_slen >= s_const_slen) && memequal(s_read, s_prefix_const, s_const_slen);
}

// Can use this when it's always safe to read first strlen(s_prefix_const)
// bytes of s_read.
HEADER_INLINE int32_t StrStartsWithUnsafe(const char* s_read, const char* s_prefix_const) {
  const uint32_t s_const_slen = strlen(s_prefix_const);
  return memequal(s_read, s_prefix_const, s_const_slen);
}

// s_suffix must be strictly contained.
HEADER_INLINE int32_t StrEndsWith(const char* s_read, const char* s_suffix_const, uint32_t s_read_slen) {
  const uint32_t s_const_slen = strlen(s_suffix_const);
  return (s_read_slen > s_const_slen) && memequal(&(s_read[s_read_slen - s_const_slen]), s_suffix_const, s_const_slen);
}

// These are likely to be revised to take const void* parameters, and moved to
// plink2_base.
// This requires len >= 4.
uintptr_t FirstUnequal4(const char* s1, const char* s2, uintptr_t slen);

HEADER_INLINE uintptr_t FirstUnequal(const char* s1, const char* s2, uintptr_t slen) {
  // Returns position of first mismatch, or slen if none was found.
  if (slen >= 4) {
    return FirstUnequal4(s1, s2, slen);
  }
  for (uintptr_t pos = 0; pos != slen; ++pos) {
    if (s1[pos] != s2[pos]) {
      return pos;
    }
  }
  return slen;
}

// May read (kBytesPerWord - 1) bytes past the end of each string.
HEADER_INLINE int32_t strequal_overread(const char* s1, const char* s2) {
  const uintptr_t* s1_alias = R_CAST(const uintptr_t*, s1);
  const uintptr_t* s2_alias = R_CAST(const uintptr_t*, s2);
  for (uintptr_t widx = 0; ; ++widx) {
    const uintptr_t w1 = s1_alias[widx];
    const uintptr_t zcheck = DetectFirstZeroByte(w1);
    const uintptr_t w2 = s2_alias[widx];
    const uintptr_t xor_word = w1 ^ w2;
    if (zcheck) {
      // Mask out bytes past the known null.
      const uintptr_t mask = zcheck ^ (zcheck - 1);
      return (xor_word & mask)? 0 : 1;
    }
    if (xor_word) {
      return 0;
    }
  }
}

int32_t strcmp_overread(const char* s1, const char* s2);

// Support for sorting arrays of strings, represented either as an array of
// const char*s, or a single [# of strings] x [max byte width] array of chars
// (suitable for old-school qsort(), referred to as a 'strbox' here).
// The SortStrboxIndexed functions automatically construct an array of const
// char*s and sort that when the max byte width is large.
int32_t strcmp_casted(const void* s1, const void* s2);

int32_t strcmp_overread_casted(const void* s1, const void* s2);


int32_t strcmp_natural(const void* s1, const void* s2);


int32_t strcmp_deref(const void* s1, const void* s2);

int32_t strcmp_overread_deref(const void* s1, const void* s2);

int32_t strcmp_natural_deref(const void* s1, const void* s2);


int32_t strcmp_natural_uncasted(const char* s1, const char* s2);

#ifdef __cplusplus
typedef struct Strbuf28UiStruct {
  char strbuf[28];
  uint32_t orig_idx;
  bool operator<(const struct Strbuf28UiStruct& rhs) const {
    return (strcmp_natural_uncasted(strbuf, rhs.strbuf) < 0);
  }
} Strbuf28Ui;

typedef struct Strbuf60UiStruct {
  char strbuf[60];
  uint32_t orig_idx;
  bool operator<(const struct Strbuf60UiStruct& rhs) const {
    return (strcmp_natural_uncasted(strbuf, rhs.strbuf) < 0);
  }
} Strbuf60Ui;
#endif

uintptr_t GetStrboxsortWentryBlen(uintptr_t max_str_blen);

#ifdef __cplusplus
typedef struct StrSortDerefStruct {
  const char* strptr;

  bool operator<(const struct StrSortDerefStruct& rhs) const {
    return (strcmp(strptr, rhs.strptr) < 0);
  }
} StrSortDeref;

HEADER_INLINE bool strcmp_overread_lt(const char* s1, const char* s2) {
  const uintptr_t* s1_alias = R_CAST(const uintptr_t*, s1);
  const uintptr_t* s2_alias = R_CAST(const uintptr_t*, s2);
  for (uintptr_t widx = 0; ; ++widx) {
    uintptr_t w1 = s1_alias[widx];
    const uintptr_t zcheck = DetectFirstZeroByte(w1);
    uintptr_t w2 = s2_alias[widx];
    if (zcheck) {
      // Mask out bytes past the known null.
      // Note that we can't safely include the garbage bytes past the null in
      // the comparison even if they aren't being changed, because they may be
      // uninitialized, and if they're also past a page boundary the OS may
      // return different values for consecutive queries on the same address!
      // See e.g. "CppCon 2016: Nicholas Ormrod 'The strange details of
      // std::string at Facebook'"
      // (https://www.youtube.com/watch?v=kPR8h4-qZdk ) starting at ~22:15.
      const uintptr_t mask = zcheck ^ (zcheck - 1);
      w1 &= mask;
      w2 &= mask;
      if (w1 == w2) {
        return false;
      }
      goto strcmp_overread_lt_finish;
    }
    if (w1 != w2) {
    strcmp_overread_lt_finish:
#  ifdef __LP64__
      return __builtin_bswap64(w1) < __builtin_bswap64(w2);
#  else
      return __builtin_bswap32(w1) < __builtin_bswap32(w2);
#  endif
    }
  }
}

typedef struct StrSortDerefOverreadStruct {
  // Must be safe to read up to (kBytesPerWord - 1) bytes past the end of these
  // strings.  Enough of a speed advantage to be worth using whenever possible,
  // though.
  const char* strptr;

  bool operator<(const struct StrSortDerefOverreadStruct& rhs) const {
    return strcmp_overread_lt(strptr, rhs.strptr);
  }
} StrSortDerefOverread;

typedef struct StrNsortDerefStruct {
  const char* strptr;
  bool operator<(const struct StrNsortDerefStruct& rhs) const {
    return (strcmp_natural(strptr, rhs.strptr) < 0);
  }
} StrNsortDeref;

typedef struct StrNsortIndexedDerefStruct {
  const char* strptr;
  uint32_t orig_idx;
  bool operator<(const struct StrNsortIndexedDerefStruct& rhs) const {
    return (strcmp_natural_uncasted(strptr, rhs.strptr) < 0);
  }
} StrNsortIndexedDeref;

HEADER_INLINE void StrptrArrSort(uintptr_t ct, const char** strptr_arr) {
  std::sort(R_CAST(StrSortDeref*, strptr_arr), &(R_CAST(StrSortDeref*, strptr_arr)[ct]));
}

HEADER_INLINE void StrptrArrSortOverread(uintptr_t ct, const char** strptr_arr) {
  std::sort(R_CAST(StrSortDerefOverread*, strptr_arr), &(R_CAST(StrSortDerefOverread*, strptr_arr)[ct]));
}

HEADER_INLINE void StrptrArrNsort(uintptr_t ct, const char** strptr_arr) {
  std::sort(R_CAST(StrNsortDeref*, strptr_arr), &(R_CAST(StrNsortDeref*, strptr_arr)[ct]));
}

// need to expose these for plink2_cmdline bigstack-allocating
// SortStrboxIndexed()'s use
void SortStrbox32bFinish(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, Strbuf28Ui* filled_wkspace, char* sorted_strbox, uint32_t* id_map);

void SortStrbox64bFinish(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, Strbuf60Ui* filled_wkspace, char* sorted_strbox, uint32_t* id_map);

// Must be ok to overread.
void SortStrboxIndexed2(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace);
#else  // !__cplusplus
HEADER_INLINE void StrptrArrSort(uintptr_t ct, const char** strptr_arr) {
  qsort(strptr_arr, ct, sizeof(intptr_t), strcmp_deref);
}

HEADER_INLINE void StrptrArrSortOverread(uintptr_t ct, const char** strptr_arr) {
  qsort(strptr_arr, ct, sizeof(intptr_t), strcmp_overread_deref);
}

HEADER_INLINE void StrptrArrNsort(uintptr_t ct, const char** strptr_arr) {
  qsort(strptr_arr, ct, sizeof(intptr_t), strcmp_natural_deref);
}

// Must be ok to overread.
void SortStrboxIndexed2Fallback(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace);

HEADER_INLINE void SortStrboxIndexed2(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace) {
  SortStrboxIndexed2Fallback(str_ct, max_str_blen, use_nsort, strbox, id_map, sort_wkspace);
}
#endif

// Uses malloc instead of bigstack.
// Must be ok to overread strbox.
BoolErr SortStrboxIndexedMalloc(uintptr_t str_ct, uintptr_t max_str_blen, char* strbox, uint32_t* id_map);

// Returns dedup'd strbox entry count.
uint32_t CopyAndDedupSortedStrptrsToStrbox(const char* const* sorted_strptrs, uintptr_t str_ct, uintptr_t max_str_blen, char* strbox);

// note that this can be expected to have size 16 bytes, not 12, on 64-bit
// systems
typedef struct StrSortIndexedDerefStruct {
  const char* strptr;
  uint32_t orig_idx;
#ifdef __cplusplus
  bool operator<(const struct StrSortIndexedDerefStruct& rhs) const {
    return (strcmp(strptr, rhs.strptr) < 0);
  }
#endif
} StrSortIndexedDeref;

typedef struct StrSortIndexedDerefOverreadStruct {
  // must be safe to read up to (kBytesPerWord - 1) bytes past the end of these
  // strings.
  const char* strptr;

  uint32_t orig_idx;
#ifdef __cplusplus
  bool operator<(const struct StrSortIndexedDerefOverreadStruct& rhs) const {
    return strcmp_overread_lt(strptr, rhs.strptr);
  }
#endif
} StrSortIndexedDerefOverread;

void StrptrArrSortMain(uintptr_t str_ct, uint32_t overread_ok, uint32_t use_nsort, StrSortIndexedDeref* wkspace_alias);


HEADER_INLINE int32_t IsLetter(unsigned char ucc) {
  return (((ucc & 192) == 64) && (((ucc - 1) & 31) < 26));
}

// if we need the digit value, better to use (unsigned char)cc - '0'...
HEADER_INLINE int32_t IsDigit(unsigned char ucc) {
  return (ucc <= '9') && (ucc >= '0');
}

HEADER_INLINE int32_t IsNotDigit(unsigned char ucc) {
  return (ucc > '9') || (ucc < '0');
}

HEADER_INLINE int32_t IsNotNzdigit(unsigned char ucc) {
  return (ucc > '9') || (ucc <= '0');
}

// May as well treat all chars < 32, except tab, as eoln...
// kns = "known non-space" (where tab counts as a space)
// This is of course identical to IsSpaceOrEoln(), but intent should be
// clearer and we can insert a debug-assert that we aren't at a space/tab.
HEADER_INLINE int32_t IsEolnKns(unsigned char ucc) {
  // could assert ucc is not a space/tab?
  return (ucc <= 32);
}

HEADER_INLINE int32_t IsEolnOrCommentKns(unsigned char ucc) {
  return (ucc < 32) || (ucc == '#');
}

HEADER_INLINE int32_t NoMoreTokensKns(const char* str) {
  return ((!str) || IsEolnKns(*str));
}

HEADER_INLINE CXXCONST_CP FirstNonChar(const char* str_iter, char cc) {
  while (*str_iter == cc) {
    ++str_iter;
  }
  return S_CAST(CXXCONST_CP, str_iter);
}

#ifdef __cplusplus
HEADER_INLINE char* FirstNonChar(char* str_iter, char cc) {
  return const_cast<char*>(FirstNonChar(const_cast<const char*>(str_iter), cc));
}
#endif

HEADER_INLINE CXXCONST_CP FirstNonTspace(const char* str_iter) {
  while ((*str_iter == ' ') || (*str_iter == '\t')) {
    ++str_iter;
  }
  return S_CAST(CXXCONST_CP, str_iter);
}

#ifdef __cplusplus
HEADER_INLINE char* FirstNonTspace(char* str_iter) {
  return const_cast<char*>(FirstNonTspace(const_cast<const char*>(str_iter)));
}
#endif

HEADER_INLINE CXXCONST_CP FirstPostspaceBounded(const char* str_iter, const char* str_end) {
  for (; str_iter != str_end; ++str_iter) {
    if (ctou32(*str_iter) > ' ') {
      break;
    }
  }
  return S_CAST(CXXCONST_CP, str_iter);
}

#ifdef __cplusplus
HEADER_INLINE char* FirstPostspaceBounded(char* str_iter, char* str_end) {
  return const_cast<char*>(FirstPostspaceBounded(const_cast<const char*>(str_iter), const_cast<const char*>(str_end)));
}
#endif


// See also (93) on TAOCP vol 4a, pp. 153.  Todo: benchmark a
// FirstPrecharUnsafe() function which uses that (unsafe because it reads up to
// 7 characters past buffer end).
HEADER_INLINE CXXCONST_CP FirstPrechar(const char* str_iter, uint32_t char_code) {
  while (ctou32(*str_iter) >= char_code) {
    ++str_iter;
  }
  return S_CAST(CXXCONST_CP, str_iter);
}

// It is worth distinguishing between incremental parsing functions where it's
// rarely necessary to scan more than 40 characters or so, and scans which are
// likely to be long-range.  The former is frequently best handled with a
// simple loop even when AVX2 movemask is available; it may not even be worth
// inserting a conditional to select between the two when length is known.  The
// latter benefits greatly from movemask.
//
// The following standard library and plink2-library scanning functions can be
// trusted to use movemask:
//   strlen, strchr, memchr
//   rawmemchr, strchrnul
//   rawmemchr2, rawmemchr3, strnul, strchrnul_n, strchrnul2, strchrnul3,
//     strchrnul_n_mov, incr_strchrnul_n_mov
//   NextTokenMultFar
//   AdvToNthDelimChecked, AdvToNthDelim, AdvToDelimOrEnd, Memrchr,
//     LastSpaceOrEoln

/*
#ifdef __LP64__
// Requires char_code <= 128.
CXXCONST_CP FirstPrecharFar(const char* str_iter, uint32_t char_code);
#else
HEADER_INLINE CXXCONST_CP FirstPrecharFar(const char* str_iter, uint32_t char_code) {
  return FirstPrechar(str_iter, char_code);
}
#endif
*/

HEADER_INLINE CXXCONST_CP FirstPrespace(const char* str_iter) {
  return S_CAST(CXXCONST_CP, FirstPrechar(str_iter, ' '));
}

HEADER_INLINE CXXCONST_CP NextPrespace(const char* str_iter) {
  return S_CAST(CXXCONST_CP, FirstPrechar(&(str_iter[1]), ' '));
}

HEADER_INLINE CXXCONST_CP FirstSpaceOrEoln(const char* str_iter) {
  return S_CAST(CXXCONST_CP, FirstPrechar(str_iter, 33));
}

// assumes we are currently in a token -- UNSAFE OTHERWISE
HEADER_INLINE CXXCONST_CP CurTokenEnd(const char* str_iter) {
  // assert(ctou32(*str_iter) > 32);
  return S_CAST(CXXCONST_CP, FirstPrechar(&(str_iter[1]), 33));
}

#ifdef __cplusplus
HEADER_INLINE char* FirstPrechar(char* str_iter, uint32_t char_code) {
  return const_cast<char*>(FirstPrechar(const_cast<const char*>(str_iter), char_code));
}

/*
HEADER_INLINE char* FirstPrecharFar(char* str_iter, uint32_t char_code) {
  return const_cast<char*>(FirstPrecharFar(const_cast<const char*>(str_iter), char_code));
}
*/

HEADER_INLINE char* FirstPrespace(char* str_iter) {
  return const_cast<char*>(FirstPrespace(const_cast<const char*>(str_iter)));
}

HEADER_INLINE char* NextPrespace(char* str_iter) {
  return const_cast<char*>(NextPrespace(const_cast<const char*>(str_iter)));
}

HEADER_INLINE char* FirstSpaceOrEoln(char* str_iter) {
  return const_cast<char*>(FirstSpaceOrEoln(const_cast<const char*>(str_iter)));
}

HEADER_INLINE char* CurTokenEnd(char* str_iter) {
  return const_cast<char*>(CurTokenEnd(const_cast<const char*>(str_iter)));
}
#endif

HEADER_INLINE CXXCONST_CP CsvFieldEnd(const char* token_iter) {
  unsigned char ucc = *token_iter;
  while ((ucc >= ' ') && (ucc != ',')) {
    ucc = *(++token_iter);
  }
  return S_CAST(CXXCONST_CP, token_iter);
}

// length-zero tokens and non-leading spaces are permitted in the
// comma-delimiter case
HEADER_INLINE CXXCONST_CP CommaOrTspaceTokenEnd(const char* token_iter, uint32_t comma_delim) {
  if (comma_delim) {
    return CsvFieldEnd(token_iter);
  }
  return S_CAST(CXXCONST_CP, CurTokenEnd(token_iter));
}

#ifdef __cplusplus
HEADER_INLINE char* CsvFieldEnd(char* token_iter) {
  return const_cast<char*>(CsvFieldEnd(const_cast<const char*>(token_iter)));
}

HEADER_INLINE char* CommaOrTspaceTokenEnd(char* token_iter, uint32_t comma_delim) {
  return const_cast<char*>(CommaOrTspaceTokenEnd(const_cast<const char*>(token_iter), comma_delim));
}
#endif

HEADER_INLINE CXXCONST_CP CommaOrTspaceFirstToken(const char* token_end_iter, uint32_t comma_delim) {
  // assumes token_end_iter is non-null, returns nullptr if there are no more
  // tokens
  // assert(token_end_iter);
  if (comma_delim) {
    if ((*token_end_iter) != ',') {
      return nullptr;
    }
    return S_CAST(CXXCONST_CP, FirstNonTspace(&(token_end_iter[1])));
  }
  const char* str = FirstNonTspace(token_end_iter);
  return IsEolnKns(*str)? nullptr : S_CAST(CXXCONST_CP, str);
}

#ifdef __cplusplus
HEADER_INLINE char* CommaOrTspaceFirstToken(char* token_end_iter, uint32_t comma_delim) {
  return const_cast<char*>(CommaOrTspaceFirstToken(const_cast<const char*>(token_end_iter), comma_delim));
}
#endif


// Returns whether uppercased str_iter matches nonempty fixed_str.  Assumes
// fixed_str contains nothing but letters and a null terminator.
// uint32_t match_upper(const char* str_iter, const char* fixed_str);

uint32_t MatchUpperCounted(const char* str, const char* fixed_str, uint32_t ct);

HEADER_INLINE uint32_t MatchUpperKLen(const char* str, const char* fixed_str, uint32_t str_slen) {
  const uint32_t fixed_slen = strlen(fixed_str);
  if (str_slen != fixed_slen) {
    return 0;
  }
  return MatchUpperCounted(str, fixed_str, fixed_slen);
}

HEADER_INLINE uint32_t MatchUpperK(const char* str, const char* fixed_str) {
  return MatchUpperCounted(str, fixed_str, strlen(fixed_str));
}

uint32_t strcaseequal(const char* str1, const char* str2, uint32_t ct);

/*
void str_toupper(char* str_iter);

void buf_toupper(uint32_t slen, char* strbuf);

void strcpy_toupper(char* target, const char* source);

char* memcpya_toupper(char* __restrict target, const char* __restrict source, uint32_t slen);
*/

uint32_t IsAlphanumeric(const char* str_iter);

// ScanPosintCapped(), ScanUintCapped(), ScanIntAbsBounded(), ScanInt32(),
// ScanPosintDefcap(), ScanUintDefcap(), ScanIntAbsDefcap(), ScanUintIcap() in
// plink2_base

BoolErr ScanPosintptr(const char* str_iter, uintptr_t* valp);

#ifdef __LP64__
BoolErr ScanmovPosintCapped(uint64_t cap, const char** str_iterp, uint32_t* valp);

BoolErr ScanmovUintCapped(uint64_t cap, const char** str_iterp, uint32_t* valp);

// 2^{-31} < -abs_floor <= 0 <= cap < 2^31
BoolErr ScanmovIntBounded(uint64_t abs_floor, uint64_t cap, const char** str_iterp, int32_t* valp);
#else
BoolErr ScanmovPosintCapped32(uint32_t cap_div_10, uint32_t cap_mod_10, const char** str_iterp, uint32_t* valp);

BoolErr ScanmovUintCapped32(uint32_t cap_div_10, uint32_t cap_mod_10, const char** str_iterp, uint32_t* valp);

BoolErr ScanmovIntBounded32(uint32_t abs_floor_div_10, uint32_t abs_floor_mod_10, uint32_t cap_div_10, uint32_t cap_mod_10, const char** str_iterp, int32_t* valp);

HEADER_INLINE BoolErr ScanmovPosintCapped(uint32_t cap, const char** str_iterp, uint32_t* valp) {
 return ScanmovPosintCapped32(cap / 10, cap % 10, str_iterp, valp);
}

HEADER_INLINE BoolErr ScanmovUintCapped(uint32_t cap, const char** str_iterp, uint32_t* valp) {
 return ScanmovUintCapped32(cap / 10, cap % 10, str_iterp, valp);
}

HEADER_INLINE BoolErr ScanmovIntBounded(uint32_t abs_floor, uint32_t cap, const char** str_iterp, int32_t* valp) {
  return ScanmovIntBounded32(abs_floor / 10, abs_floor % 10, cap / 10, cap % 10, str_iterp, valp);
}
#endif

HEADER_INLINE BoolErr ScanmovUintDefcap(const char** str_iterp, uint32_t* valp) {
  return ScanmovUintCapped(0x7ffffffe, str_iterp, valp);
}

// This has different semantics from ScanmovPosintCapped, etc. since integer
// readers don't take much code (so it's fine to have a bunch of similar
// functions, optimized for slightly different use cases), but we only want one
// core floating point reader.
// (update, 3 Feb 2018: renamed the integer readers above to start with
// scanmov_ instead of scanadv_, to reflect the interface difference between
// returning a pointer and moving the input pointer forward.)
CXXCONST_CP ScanadvDouble(const char* str_iter, double* valp);

// Thin wrapper that verifies the token ends with whitespace/eoln; should be
// used whenever comma/semicolon/etc. delimiters are not ok.
HEADER_INLINE CXXCONST_CP ScantokDouble(const char* str_iter, double* valp) {
  CXXCONST_CP parsed_end = ScanadvDouble(str_iter, valp);
  if ((!parsed_end) || (!IsSpaceOrEoln(*parsed_end))) {
    return nullptr;
  }
  return parsed_end;
}

#ifdef __cplusplus
HEADER_INLINE char* ScanadvDouble(char* str_iter, double* valp) {
  return const_cast<char*>(ScanadvDouble(const_cast<const char*>(str_iter), valp));
}

HEADER_INLINE char* ScantokDouble(char* str_iter, double* valp) {
  return const_cast<char*>(ScantokDouble(const_cast<const char*>(str_iter), valp));
}
#endif

// remove unlikely() if any caller ever tries to reparse the string as
// something else in error case (this is a valid ScanadvDouble() use case)
HEADER_INLINE BoolErr ScanFloat(const char* ss, float* valp) {
  double dxx;
  if (unlikely(!ScantokDouble(ss, &dxx))) {
    return 1;
  }
  if (unlikely(fabs(dxx) > 3.4028235677973362e38)) {
    return 1;
  }
  *valp = S_CAST(float, dxx);
  return 0;
}

CXXCONST_CP ScanadvLn(const char* str_iter, double* ln_ptr);

HEADER_INLINE CXXCONST_CP ScantokLn(const char* str_iter, double* ln_ptr) {
  CXXCONST_CP parsed_end = ScanadvLn(str_iter, ln_ptr);
  if ((!parsed_end) || (!IsSpaceOrEoln(*parsed_end))) {
    return nullptr;
  }
  return parsed_end;
}

#ifdef __cplusplus
HEADER_INLINE char* ScanadvLn(char* str_iter, double* ln_ptr) {
  return const_cast<char*>(ScanadvLn(const_cast<const char*>(str_iter), ln_ptr));
}

HEADER_INLINE char* ScantokLn(char* str_iter, double* ln_ptr) {
  return const_cast<char*>(ScantokLn(const_cast<const char*>(str_iter), ln_ptr));
}
#endif

// These provide the same interface as ScanPosintCapped(), etc., but there are
// two differences in behavior to make these more suitable for parsing of
// command-line parameters:
// - The strings are initially parsed as floating-point, and then (if an
//   integer is expected) the function errors out if the result isn't an exact
//   integer.  This allows exponential notation to be used.
// - Unlike atoi()/ScanPosintCapped(), the function errors out if parsing stops
//   at non-whitespace.
// The performance cost of this behavior is relatively high: these functions
// shouldn't be used for internal file-reading loops.
BoolErr ScanPosintCappedx(const char* str_iter, uint64_t cap, uint32_t* valp);

BoolErr ScanUintCappedx(const char* str_iter, uint64_t cap, uint32_t* valp);

BoolErr ScanIntAbsBoundedx(const char* str_iter, int64_t bound, int32_t* valp);

HEADER_INLINE BoolErr ScanInt32x(const char* str, int32_t* valp) {
  return ScanIntAbsBoundedx(str, 0x7fffffff, valp);
}

HEADER_INLINE BoolErr ScanPosintDefcapx(const char* str, uint32_t* valp) {
  return ScanPosintCappedx(str, 0x7ffffffe, valp);
}

HEADER_INLINE BoolErr ScanUintDefcapx(const char* str, uint32_t* valp) {
  return ScanUintCappedx(str, 0x7ffffffe, valp);
}

BoolErr ScanPosintptrx(const char* str_iter, uintptr_t* valp);


HEADER_INLINE void AppendBinaryEoln(char** dst_ptr) {
#ifdef _WIN32
  (*dst_ptr)[0] = '\r';
  (*dst_ptr)[1] = '\n';
  *dst_ptr += 2;
#else
  **dst_ptr = '\n';
  *dst_ptr += 1;
#endif
}

HEADER_INLINE void DecrAppendBinaryEoln(char** dst_ptr) {
#ifdef _WIN32
  (*dst_ptr)[-1] = '\r';
  (*dst_ptr)[0] = '\n';
  *dst_ptr += 1;
#else
  (*dst_ptr)[-1] = '\n';
#endif
}

void GetTopTwoUi(const uint32_t* __restrict uint_arr, uintptr_t uia_size, uintptr_t* __restrict top_idx_ptr, uintptr_t* __restrict second_idx_ptr);

// safer than CurTokenEnd(), since it handles length zero
// "se" = space/eoln treated as terminators
HEADER_INLINE uintptr_t strlen_se(const char* ss) {
  const char* ss2 = ss;
  while (!IsSpaceOrEoln(*ss2)) {
    ss2++;
  }
  return ss2 - ss;
}

// just an alias for rawmemchr which doesn't require a subsequent static-cast.
HEADER_INLINE CXXCONST_CP AdvToDelim(const char* str_iter, char delim) {
  return S_CAST(CXXCONST_CP, rawmemchr(str_iter, delim));
}

#ifdef __cplusplus
HEADER_INLINE char* AdvToDelim(char* str_iter, char delim) {
  return const_cast<char*>(AdvToDelim(const_cast<const char*>(str_iter), delim));
}
#endif

HEADER_INLINE CXXCONST_CP AdvPastDelim(const char* str_iter, char delim) {
  return &(AdvToDelim(str_iter, delim)[1]);
}

#ifdef __cplusplus
HEADER_INLINE char* AdvPastDelim(char* str_iter, char delim) {
  return &(AdvToDelim(str_iter, delim)[1]);
}
#endif


#ifdef __LP64__
// This is a major VCF-parsing bottleneck, and inlining it makes a big
// difference.

// ct must be nonzero
HEADER_INLINE CXXCONST_CP AdvToNthDelimChecked(const char* str_iter, const char* str_end, uint32_t ct, char delim) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, str_iter);
  const uintptr_t ending_addr = R_CAST(uintptr_t, str_end);
  VecUc* str_viter = R_CAST(VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecUc vvec_all_delim = vecuc_set1(delim);
  VecUc cur_vvec = *str_viter;
  VecUc delim_vvec = (cur_vvec == vvec_all_delim);
  uint32_t delimiter_bytes = vecuc_movemask(delim_vvec);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, str_viter);
  const uint32_t leading_mask = UINT32_MAX << leading_byte_ct;
  delimiter_bytes &= leading_mask;
  for (uint32_t remaining_delim_ct = ct; ; ) {
    const uint32_t cur_delim_ct = PopcountVec8thUint(delimiter_bytes);
    if (cur_delim_ct >= remaining_delim_ct) {
      delimiter_bytes = ClearBottomSetBits(remaining_delim_ct - 1, delimiter_bytes);
      const uint32_t byte_offset_in_vec = ctzu32(delimiter_bytes);
      const uintptr_t result_addr = R_CAST(uintptr_t, str_viter) + byte_offset_in_vec;
      if (result_addr >= ending_addr) {
        return nullptr;
      }
      return R_CAST(CXXCONST_CP, result_addr);
    }
    remaining_delim_ct -= cur_delim_ct;
    ++str_viter;
    if (R_CAST(uintptr_t, str_viter) >= ending_addr) {
      return nullptr;
    }
    cur_vvec = *str_viter;
    delim_vvec = (cur_vvec == vvec_all_delim);
    delimiter_bytes = vecuc_movemask(delim_vvec);
  }
}

HEADER_INLINE CXXCONST_CP AdvToNthDelim(const char* str_iter, uint32_t ct, char delim) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, str_iter);
  VecUc* str_viter = R_CAST(VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
  const VecUc vvec_all_delim = vecuc_set1(delim);
  VecUc cur_vvec = *str_viter;
  VecUc delim_vvec = (cur_vvec == vvec_all_delim);
  uint32_t delimiter_bytes = vecuc_movemask(delim_vvec);
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, str_viter);
  const uint32_t leading_mask = UINT32_MAX << leading_byte_ct;
  delimiter_bytes &= leading_mask;
  for (uint32_t remaining_delim_ct = ct; ; ) {
    const uint32_t cur_delim_ct = PopcountVec8thUint(delimiter_bytes);
    if (cur_delim_ct >= remaining_delim_ct) {
      delimiter_bytes = ClearBottomSetBits(remaining_delim_ct - 1, delimiter_bytes);
      const uint32_t byte_offset_in_vec = ctzu32(delimiter_bytes);
      const uintptr_t result_addr = R_CAST(uintptr_t, str_viter) + byte_offset_in_vec;
      return R_CAST(CXXCONST_CP, result_addr);
    }
    remaining_delim_ct -= cur_delim_ct;
    ++str_viter;
    cur_vvec = *str_viter;
    delim_vvec = (cur_vvec == vvec_all_delim);
    delimiter_bytes = vecuc_movemask(delim_vvec);
  }
}
#else  // !__LP64__
HEADER_INLINE CXXCONST_CP AdvToNthDelimChecked(const char* str_iter, const char* str_end, uint32_t ct, char delim) {
  for (uint32_t remaining_delim_ct = ct; ; ) {
    const char* next_delim = S_CAST(const char*, memchr(str_iter, delim, str_end - str_iter));
    if (!next_delim) {
      return nullptr;
    }
    if (!(--remaining_delim_ct)) {
      return S_CAST(CXXCONST_CP, next_delim);
    }
    str_iter = &(next_delim[1]);
  }
}

HEADER_INLINE CXXCONST_CP AdvToNthDelim(const char* str_iter, uint32_t ct, char delim) {
  for (uint32_t remaining_delim_ct = ct; ; ) {
    const char* next_delim = AdvToDelim(str_iter, delim);
    if (!(--remaining_delim_ct)) {
      return S_CAST(CXXCONST_CP, next_delim);
    }
    str_iter = &(next_delim[1]);
  }
}
#endif

#ifdef __cplusplus
HEADER_INLINE char* AdvToNthDelimChecked(char* str_iter, char* str_end, uint32_t ct, char delim) {
  return const_cast<char*>(AdvToNthDelimChecked(const_cast<const char*>(str_iter), const_cast<const char*>(str_end), ct, delim));
}

HEADER_INLINE char* AdvToNthDelim(char* str_iter, uint32_t ct, char delim) {
  return const_cast<char*>(AdvToNthDelim(const_cast<const char*>(str_iter), ct, delim));
}
#endif

// ok if str_iter is at end of current token
HEADER_INLINE CXXCONST_CP NextToken(const char* str_iter) {
  if (!str_iter) {
    return nullptr;
  }
  unsigned char ucc = *str_iter;
  while (ucc > ' ') {
    ucc = *(++str_iter);
  }
  while ((ucc == ' ') || (ucc == '\t')) {
    ucc = *(++str_iter);
  }
  return (ucc > 32)? S_CAST(CXXCONST_CP, str_iter) : nullptr;
}

#ifdef __cplusplus
HEADER_INLINE char* NextToken(char* str_iter) {
  return const_cast<char*>(NextToken(const_cast<const char*>(str_iter)));
}
#endif

HEADER_INLINE CXXCONST_CP NextTokenMult(const char* str_iter, uint32_t ct) {
  // assert(ct);
  if (!str_iter) {
    return nullptr;
  }
  unsigned char ucc = *str_iter;
  do {
    while (ucc > 32) {
      ucc = *(++str_iter);
    }
    while ((ucc == ' ') || (ucc == '\t')) {
      ucc = *(++str_iter);
    }
    if (ucc <= 32) {
      return nullptr;
    }
  } while (--ct);
  return S_CAST(CXXCONST_CP, str_iter);
}

#ifdef USE_AVX2
// todo: determine minimum ct where this pays off.
CXXCONST_CP NextTokenMultFar(const char* str_iter, uint32_t ct);
#else
HEADER_INLINE CXXCONST_CP NextTokenMultFar(const char* str_iter, uint32_t ct) {
  return NextTokenMult(str_iter, ct);
}
#endif

#ifdef __cplusplus
HEADER_INLINE char* NextTokenMult(char* str_iter, uint32_t ct) {
  return const_cast<char*>(NextTokenMult(const_cast<const char*>(str_iter), ct));
}

HEADER_INLINE char* NextTokenMultFar(char* str_iter, uint32_t ct) {
  return const_cast<char*>(NextTokenMultFar(const_cast<const char*>(str_iter), ct));
}
#endif

HEADER_INLINE CXXCONST_CP NextTokenMult0(const char* str, uint32_t ct) {
  // tried replacing this with ternary operator, but that actually seemed to
  // slow things down a bit under gcc 4.2.1 (tail call optimization issue?).
  // todo: recheck this under newer gcc/clang.
  if (ct) {
    return S_CAST(CXXCONST_CP, NextTokenMult(str, ct));
  }
  return S_CAST(CXXCONST_CP, str);
}

#ifdef __cplusplus
HEADER_INLINE char* NextTokenMult0(char* str, uint32_t ct) {
  return const_cast<char*>(NextTokenMult0(const_cast<const char*>(str), ct));
}
#endif

#ifdef USE_AVX2
const char* TokenLexK0(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens);

HEADER_INLINE const char* TokenLexK(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens) {
  return TokenLexK0(str_iter, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens);
}
#else
// assumes str_iter != nullptr
// returns nullptr on missing token, otherwise returns pointer to end of last
//   lexed token
HEADER_INLINE const char* TokenLexK(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens) {
  for (uint32_t relevant_col_idx = 0; relevant_col_idx != relevant_col_ct; ++relevant_col_idx) {
    const uint32_t cur_col_type = col_types[relevant_col_idx];
    str_iter = NextTokenMult(str_iter, col_skips[relevant_col_idx]);
    if (!str_iter) {
      return nullptr;
    }
    token_ptrs[cur_col_type] = str_iter;
    const char* token_end = CurTokenEnd(str_iter);
    token_slens[cur_col_type] = token_end - str_iter;
    str_iter = token_end;
  }
  return str_iter;
}

HEADER_INLINE const char* TokenLexK0(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens) {
  if (!col_skips[0]) {
    const uint32_t cur_col_type = col_types[0];
    const char* first_token_end = CurTokenEnd(str_iter);
    token_ptrs[cur_col_type] = str_iter;
    token_slens[cur_col_type] = first_token_end - str_iter;
    str_iter = first_token_end;
    ++col_types;
    ++col_skips;
    --relevant_col_ct;
  }
  return TokenLexK(str_iter, col_types, col_skips, relevant_col_ct, token_ptrs, token_slens);
}
#endif

HEADER_INLINE char* TokenLex(char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, char** token_ptrs, uint32_t* token_slens) {
  return K_CAST(char*, TokenLexK(str_iter, col_types, col_skips, relevant_col_ct, K_CAST(const char**, token_ptrs), token_slens));
}

HEADER_INLINE char* TokenLex0(char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, char** token_ptrs, uint32_t* token_slens) {
  return K_CAST(char*, TokenLexK0(str_iter, col_types, col_skips, relevant_col_ct, K_CAST(const char**, token_ptrs), token_slens));
}

// ct must be positive for these functions.
CXXCONST_CP NextCsvMult(const char* str_iter, uint32_t ct);

HEADER_INLINE CXXCONST_CP CommaOrTspaceNextTokenMult(const char* str_iter, uint32_t ct, uint32_t comma_delim) {
  if (!comma_delim) {
    return NextTokenMult(str_iter, ct);
  }
  return NextCsvMult(str_iter, ct);
}

#ifdef __cplusplus
HEADER_INLINE char* NextCsvMult(char* str_iter, uint32_t ct) {
  return const_cast<char*>(NextCsvMult(const_cast<const char*>(str_iter), ct));
}

HEADER_INLINE char* CommaOrTspaceNextTokenMult(char* str_iter, uint32_t ct, uint32_t comma_delim) {
  return const_cast<char*>(CommaOrTspaceNextTokenMult(const_cast<const char*>(str_iter), ct, comma_delim));
}
#endif

#ifdef USE_AVX2
const char* CsvLexK(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens);
#else
HEADER_INLINE const char* CsvLexK(const char* str_iter, const uint32_t* col_types, const uint32_t* col_skips, uint32_t relevant_col_ct, const char** token_ptrs, uint32_t* token_slens) {
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

// todo: movemask version of this
uint32_t CountTokens(const char* str_iter);

// uint32_t CommaOrSpaceCountTokens(const char* str_iter, uint32_t comma_delim);

// empty multistr ok
uint32_t CountAndMeasureMultistr(const char* multistr, uintptr_t* max_blen_ptr);

extern const uint16_t kDigitPair[];

char* u32toa(uint32_t uii, char* start);

char* i32toa(int32_t ii, char* start);

char* u32toa_z5(uint32_t uii, char* start);

char* i64toa(int64_t llii, char* start);

#ifdef __LP64__
// really just for printing line numbers
// must be less than 2^63
HEADER_INLINE char* wtoa(uintptr_t ulii, char* start) {
  return i64toa(ulii, start);
}
#else
HEADER_INLINE char* wtoa(uintptr_t ulii, char* start) {
  return u32toa(ulii, start);
}
#endif

char* u32toa_trunc4(uint32_t uii, char* start);

// Write-buffer-sizing constants, to support replacement of dtoa_g() with
// higher-precision variants like dtoa_g_p8() or full-blown Ryu.

// -0.000123456
// -1.23456e-38
CONSTI32(kMaxFloatGSlen, 12);

CONSTI32(kMaxDoubleGSlen, 13);

// 1.23456e-2147483647
CONSTI32(kMaxLnGSlen, 19);

char* dtoa_g(double dxx, char* start);

// We try to avoid micromanaging floating point printing and just use %g
// everywhere, but occasionally we explicitly need more precision.
//
// dtoa_g_p8 provides generic 8-digit precision (instead of %g's 6-digit
// default), while print_dosage in plink2_common provides up to 3 places after
// the decimal point when dealing with dosages (which are internally
// represented as 32768ths).
// (may want to replace _p8 with _p10 for perfect int32 handling.)
char* dtoa_g_p8(double dxx, char* start);

static const double kLn10 = 2.3025850929940457;
static const double kRecipLn10 = 0.43429448190325176;
static const double kLnNormalMin = -708.3964185322641;

char* lntoa_g(double ln_val, char* start);

HEADER_INLINE void TrailingZeroesToSpaces(char* start) {
  --start;
  while (*start == '0') {
    *start-- = ' ';
  }
  if (*start == '.') {
    *start = ' ';
  }
}

HEADER_INLINE char* ClipTrailingZeroes(char* start) {
  char cc;
  do {
    cc = *(--start);
  } while (cc == '0');
  return &(start[(cc != '.')]);
}

// "prob" means that the number is guaranteed to be in [0, 1].
// no leading space is printed.  trailing zeroes (/decimal point) are erased
//   iff there is equality to ~13 decimal places.
char* dtoa_f_probp6_spaced(double dxx, char* start);

char* dtoa_f_probp6_clipped(double dxx, char* start);

// char* dtoa_f_p5_clipped(double dxx, char* start);

// dedicated ftoa_g() discontinued
HEADER_INLINE char* ftoa_g(float fxx, char* start) {
  return dtoa_g(S_CAST(double, fxx), start);
}

HEADER_INLINE char* u32toa_x(uint32_t uii, char extra_char, char* start) {
  char* penult = u32toa(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* i32toa_x(int32_t ii, char extra_char, char* start) {
  char* penult = i32toa(ii, start);
  *penult = extra_char;
  return &(penult[1]);
}


// overread must be ok.
CXXCONST_CP ScanForDuplicateIds(const char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_blen);

#ifdef __cplusplus
HEADER_INLINE char* ScanForDuplicateIds(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_blen) {
  return const_cast<char*>(ScanForDuplicateIds(const_cast<const char*>(sorted_ids), id_ct, max_id_blen));
}
#endif

// Collapses array of sorted IDs to remove duplicates, and writes pre-collapse
// positions to id_starts (so e.g. duplication count of any sample ID can be
// determined via subtraction) if it isn't nullptr.
// Overread must be ok.
// Returns id_ct of collapsed array.
uint32_t CollapseDuplicateIds(uintptr_t id_ct, uintptr_t max_id_blen, char* sorted_ids, uint32_t* id_starts);


// returns position of string, or -1 if not found.
int32_t bsearch_str(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx);

// requires null-terminated string
int32_t bsearch_str_natural(const char* idbuf, const char* sorted_strbox, uintptr_t max_id_blen, uintptr_t end_idx);


// returns number of elements in sorted_strbox[] less than idbuf.
uintptr_t bsearch_str_lb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx);

// same result as bsearch_str_lb(), but checks against [cur_idx],
// [cur_idx + 1], [cur_idx + 3], [cur_idx + 7], etc. before finishing with a
// binary search, and assumes cur_id_slen <= max_id_blen and end_idx > 0.
uintptr_t ExpsearchStrLb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx, uintptr_t cur_idx);

// Null-terminated string required.
uintptr_t ExpsearchNsortStrLb(const char* idbuf, const char* nsorted_strbox, uintptr_t max_id_blen, uintptr_t end_idx, uintptr_t cur_idx);

// this is frequently preferable to bsearch_str(), since it's way too easy to
// forget to convert the sorted-stringbox index to the final index
// sample_id_map == nullptr is permitted; in this case id will be an index into
// the sorted array
HEADER_INLINE BoolErr SortedIdboxFind(const char* idbuf, const char* sorted_idbox, const uint32_t* id_map, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx, uint32_t* id_ptr) {
  const int32_t ii = bsearch_str(idbuf, sorted_idbox, cur_id_slen, max_id_blen, end_idx);
  if (ii == -1) {
    return 1;
  }
  *id_ptr = id_map? id_map[S_CAST(uint32_t, ii)] : S_CAST(uint32_t, ii);
  return 0;
}

#ifdef __arm__
#  error "Unaligned accesses in IsNanStr()."
#endif
// This returns 1 on any capitalization of 'na' or 'nan', 0 otherwise.
// todo: check whether there's actually any point to the uint16_t type-pun
HEADER_INLINE uint32_t IsNanStr(const char* ss, uint32_t slen) {
  if ((slen > 3) || (slen == 1)) {
    return 0;
  }
  if (!slen) {
    return 1;
  }
  const uint32_t first_two_chars_code = R_CAST(const uint16_t*, ss)[0];
  // assumes little-endian
  if ((first_two_chars_code & 0xdfdf) != 0x414e) {
    return 0;
  }
  return (slen == 2) || ((ctou32(ss[2]) & 0xdf) == 78);
}

// Matches "inf"/"infinity", any capitalization, can have sign in front.
// Assumes is_neg zero-initialized.
// Assumes one-char overread is ok.
uint32_t IsInfStr(const char* ss, uint32_t slen, uint32_t* is_neg_ptr);


HEADER_INLINE CXXCONST_CP AdvToDelimOrEnd(const char* str_iter, const char* str_end, char delim) {
  CXXCONST_CP memchr_result = S_CAST(CXXCONST_CP, memchr(str_iter, delim, str_end - str_iter));
  if (memchr_result) {
    return memchr_result;
  }
  return S_CAST(CXXCONST_CP, str_end);
}

// tried memchr_likely_short(), not worth it

#ifdef __cplusplus
HEADER_INLINE char* AdvToDelimOrEnd(char* str_iter, char* str_end, char needle) {
  return const_cast<char*>(AdvToDelimOrEnd(const_cast<const char*>(str_iter), const_cast<const char*>(str_end), needle));
}
#endif

// memrchr() not available on some platforms.  This implementation also
// includes a tweak which trades off a bit of performance around length 35 for
// substantially better performance on the longer lines often seen in e.g. VCF
// files, so we don't normally use the base implementation when it's available,
// hence the initial capital letter.
#ifdef __LP64__
CXXCONST_CP Memrchr(const char* str_start, char needle, uintptr_t slen);

CXXCONST_CP LastSpaceOrEoln(const char* str_start, uintptr_t slen);
#else  // !__LP64__
HEADER_INLINE CXXCONST_CP Memrchr(const char* str_start, char needle, uintptr_t slen) {
#  ifdef _GNU_SOURCE
  return S_CAST(CXXCONST_CP, memrchr(str_start, ctou32(needle), slen));
#  else  // !_GNU_SOURCE
  // Could check one word at a time for not-that-small slen.
  for (uintptr_t pos = slen; pos; ) {
    if (str_start[--pos] == needle) {
      return S_CAST(CXXCONST_CP, &(str_start[pos]));
    }
  }
  return nullptr;
#  endif  // !_GNU_SOURCE
}

HEADER_INLINE CXXCONST_CP LastSpaceOrEoln(const char* str_start, uintptr_t slen) {
  for (uintptr_t pos = slen; pos; ) {
    if (ctou32(str_start[--pos]) <= 32) {
      return S_CAST(CXXCONST_CP, &(str_start[pos]));
    }
  }
  return nullptr;
}
#endif  // !__LP64__

#ifdef __cplusplus
HEADER_INLINE char* Memrchr(char* str_start, char needle, uintptr_t slen) {
  return const_cast<char*>(Memrchr(const_cast<const char*>(str_start), needle, slen));
}

HEADER_INLINE char* LastSpaceOrEoln(char* str_start, uintptr_t slen) {
  return const_cast<char*>(LastSpaceOrEoln(const_cast<const char*>(str_start), slen));
}
#endif

// void ReplaceAllInstances(char old_char, char new_char, uint32_t slen, char* dst);

void TabsToSpaces(char* ss_iter);

// Errors out if new_char is already present.
BoolErr ReplaceCharAdvChecked(char old_char, char new_char, char** str_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_STRING_H__
