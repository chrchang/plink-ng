#ifndef __PLINK2_BASE_H__
#define __PLINK2_BASE_H__

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


// Low-level C99/C++03/C++11 library covering basic I/O, bitarrays, and
// Windows/OS X/Linux portability.  We try to benefit from as much C++ type
// safety as we can without either breaking compatibility with C-only codebases
// or making extension of pgenlib/plink2 code more difficult than the old
// type-unsafe style.
//
// Parameter conventions:
// - Input parameters, then in/out, then pure outputs, then temporary buffers.
//   Reference-style input parameters tend to go in the very front, to make it
//   more obvious that they aren't in/out.
// - "bitarr" indicates a word-aligned, packed array of bits, while "bitvec"
//   indicates vector-alignment in 64-bit builds.  ("vector" always means SIMD
//   inputs/outputs here; C++ std::vector is not used in this codebase.)
// - Most pointers are stationary; moving pointers have an _iter suffix.


// 10000 * major + 100 * minor + patch
// Exception to CONSTU31, since we want the preprocessor to have access to this
// value.  Named with all caps as a consequence.
#define PLINK2_BASE_VERNUM 100


#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
  #define __STDC_FORMAT_MACROS 1
#endif
#include <inttypes.h>
#include <limits.h> // CHAR_BIT, PATH_MAX

// #define NDEBUG
#include <assert.h>

#ifdef _WIN32
  // needed for MEMORYSTATUSEX
  #ifndef _WIN64
    #define WINVER 0x0500
  #else
    #define __LP64__
  #endif
  #include <windows.h>
#endif

#ifdef __LP64__
  #ifndef __SSE2__
    // todo: remove this requirement, the 32-bit vul_t-using code does most of
    // what we need
    #error "64-bit builds currently require SSE2.  Try producing a 32-bit build instead."
  #endif
  #include <emmintrin.h>
  #ifdef __SSE4_2__
    #define USE_SSE42
    #ifdef __AVX2__
      #include <immintrin.h>
      #define USE_AVX2
    #endif
  #endif
#endif


// done with #includes, can start C++ namespace
#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef __cplusplus
  #define HEADER_INLINE inline
  #if __cplusplus >= 201103L
    #define HEADER_CINLINE constexpr
    #define CSINLINE static constexpr
    #if __cplusplus > 201103L
      #define HEADER_CINLINE2 constexpr
      #define CSINLINE2 static constexpr
    #else
      #define HEADER_CINLINE2 inline
      #define CSINLINE2 static inline
    #endif
  #else
    #define HEADER_CINLINE inline
    #define HEADER_CINLINE2 inline
    #define CSINLINE static inline
    #define CSINLINE2 static inline
  #endif
  #if __cplusplus <= 199711L
    // this may be defined anyway, at least on OS X
    #ifndef static_assert
      // todo: check other cases
      #define static_assert(cond, msg)
    #endif
  #endif
#else
  #define HEADER_INLINE static inline
  #define HEADER_CINLINE static inline
  #define HEADER_CINLINE2 static inline
  #define CSINLINE static inline
  #define CSINLINE2 static inline
  // _Static_assert() should work in gcc 4.6+
  #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 6)
    #if defined(__APPLE__) && defined(__has_feature) && defined(__has_extension)
      // clang
      #if __has_feature(c_static_assert) || __has_extension(c_static_assert)
        #define static_assert _Static_assert
      #else
        #define static_assert(cond, msg)
      #endif
    #else
      #define static_assert(cond, msg)
    #endif
  #else
    #define static_assert _Static_assert
  #endif
#endif

#define __maybe_unused __attribute__((unused))

// Error return types.  All of these evaluate to true on error and false on
// success, but otherwise they have slightly different semantics:
// * pglerr_t is the general-purpose enum.  Unlike an enum, implicit conversion
//   *to* int, not just from int, is prevented by the C++11 compiler (and the
//   C++11-compiler-validated code still works under C99).  (To achieve this
//   additional safety, we engage in a bit of code duplication which would be
//   unreasonable for flagsets.)
//   Explicit cast to uint32_t, but not int32_t, is supported, to reflect the
//   fact that all error codes are positive.
// * boolerr_t allows implicit conversion from int, but conversion back to
//   uint32_t requires an explicit cast.  (It should always be 0/1-valued, but
//   this isn't enforced by the compiler.)
// * interr_t allows implicit conversion from int, but conversion back to
//   int32_t requires an explicit cast.  It mainly serves as a holding pen for
//   C standard library error return values, which can be negative.
#if __cplusplus >= 201103L
struct pglerr_t {
  enum class ec
#else
typedef enum
#endif
  {
  kPglRetSuccess,
  kPglRetSkipped,
  kPglRetNomem,
  kPglRetOpenFail,
  kPglRetReadFail,
  kPglRetWriteFail,
  // MalformedInput should be returned on low-level file format violations,
  // while InconsistentInput should be returned for higher-level logical
  // problems like mismatched files.
  kPglRetMalformedInput,
  kPglRetInconsistentInput,
  kPglRetInvalidCmdline,
  kPglRetHelp,
  kPglRetThreadCreateFail,
  kPglRetNetworkFail,
  kPglRetSampleMajorBed = 32,
  kPglRetWarningErrcode = 61,
  kPglRetImproperFunctionCall = 62,
  kPglRetNotYetSupported = 63,
  kPglRetLongLine = 126,
  kPglRetEmptyFile = 127}
#if __cplusplus >= 201103L
  ;

  pglerr_t() {}
  
  pglerr_t(const pglerr_t& source) : value_(source.value_) {}

  pglerr_t(ec source) : value_(source) {}

  operator ec() const {
    return value_;
  }

  explicit operator uint32_t() const {
    return static_cast<uint32_t>(value_);
  }
  
  explicit operator bool() const {
    return (static_cast<uint32_t>(value_) != 0);
  }

private:
  ec value_;
};

const pglerr_t kPglRetSuccess = pglerr_t::ec::kPglRetSuccess;
const pglerr_t kPglRetSkipped = pglerr_t::ec::kPglRetSkipped;
const pglerr_t kPglRetNomem = pglerr_t::ec::kPglRetNomem;
const pglerr_t kPglRetOpenFail = pglerr_t::ec::kPglRetOpenFail;
const pglerr_t kPglRetReadFail = pglerr_t::ec::kPglRetReadFail;
const pglerr_t kPglRetWriteFail = pglerr_t::ec::kPglRetWriteFail;
const pglerr_t kPglRetMalformedInput = pglerr_t::ec::kPglRetMalformedInput;
const pglerr_t kPglRetInconsistentInput = pglerr_t::ec::kPglRetInconsistentInput;
const pglerr_t kPglRetInvalidCmdline = pglerr_t::ec::kPglRetInvalidCmdline;
const pglerr_t kPglRetHelp = pglerr_t::ec::kPglRetHelp;
const pglerr_t kPglRetThreadCreateFail = pglerr_t::ec::kPglRetThreadCreateFail;
const pglerr_t kPglRetNetworkFail = pglerr_t::ec::kPglRetNetworkFail;
const pglerr_t kPglRetSampleMajorBed = pglerr_t::ec::kPglRetSampleMajorBed;
const pglerr_t kPglRetWarningErrcode = pglerr_t::ec::kPglRetWarningErrcode;
const pglerr_t kPglRetImproperFunctionCall = pglerr_t::ec::kPglRetImproperFunctionCall;
const pglerr_t kPglRetNotYetSupported = pglerr_t::ec::kPglRetNotYetSupported;
const pglerr_t kPglRetLongLine = pglerr_t::ec::kPglRetLongLine;
const pglerr_t kPglRetEmptyFile = pglerr_t::ec::kPglRetEmptyFile;
#else
  pglerr_t;
#endif

#if __cplusplus >= 201103L
// allow efficient arithmetic on these, but force them to require explicit
// int32_t/uint32_t casts; only permit implicit assignment from
// int32_t/uint32_t by default.
// built-in bool type does too many things we don't want...

// expected to be integer-valued, but not necessarily 0/1 or positive
struct interr_t {
  interr_t() {}
  
  interr_t(int32_t source) : value_(source) {}

  explicit operator int32_t() const {
    return static_cast<int32_t>(value_);
  }
  
  explicit operator bool() const {
    return (value_ != 0);
  }
  
private:
  int32_t value_;
};

// expected to be 0/1-valued
struct boolerr_t {
  boolerr_t() {}
  
  boolerr_t(uint32_t source) : value_(source) {}

  explicit operator uint32_t() const {
    return static_cast<uint32_t>(value_);
  }
  
  explicit operator bool() const {
    return (value_ != 0);
  }
  
private:
  uint32_t value_;
};
#else
  typedef int32_t interr_t;
  typedef uint32_t boolerr_t;
#endif

// make this work on 32-bit as well as 64-bit systems, across
// Windows/OS X/Linux
// (todo: clean this up a bit.  it's inherently a baling-wire-and-duct-tape
// sort of thing, though...)
#ifdef _WIN32
  // must compile with -std=gnu++11, not c++11, on 32-bit Windows since
  // otherwise fseeko64 not defined...
  #define fseeko fseeko64
  #define ftello ftello64
  #define FOPEN_RB "rb"
  #define FOPEN_WB "wb"
  #define FOPEN_AB "ab"
  #ifdef __LP64__
    #define getc_unlocked _fgetc_nolock
    #define putc_unlocked _fputc_nolock
  #else
    #define getc_unlocked getc
    #define putc_unlocked putc
  #endif
  #if __cplusplus < 201103L
    #define uint64_t unsigned long long
    #define int64_t long long
  #endif
#else
  #define FOPEN_RB "r"
  #define FOPEN_WB "w"
  #define FOPEN_AB "a"
#endif

#ifdef _WIN32
  #define PRId64 "I64d"
  #define PRIu64 "I64u"
#else
  #ifdef __cplusplus
    #ifndef PRId64
      #define PRId64 "lld"
    #endif
  #endif
#endif

#ifdef _WIN64
  #define CTZLU __builtin_ctzll
  #define CLZLU __builtin_clzll
#else
  #define CTZLU __builtin_ctzl
  #define CLZLU __builtin_clzl
  #ifndef __LP64__
    // needed to prevent GCC 6 build failure
    #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
      #if (__cplusplus < 201103L) && !defined(__APPLE__)
	#ifndef uintptr_t
	  #define uintptr_t unsigned long
	#endif
	#ifndef intptr_t
	  #define intptr_t long
	#endif
      #endif
    #endif
  #endif
#endif

#ifdef __LP64__
  #ifdef _WIN32 // i.e. Win64

    #undef PRIuPTR
    #undef PRIdPTR
    #define PRIuPTR PRIu64
    #define PRIdPTR PRId64
    #define PRIxPTR2 "016I64x"

  #else // not _WIN32

    #ifndef PRIuPTR
      #define PRIuPTR "lu"
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR "ld"
    #endif
    #define PRIxPTR2 "016lx"

  #endif // Win64

#else // not __LP64__

  // without this, we get ridiculous warning spew...
  // not 100% sure this is the right cutoff, but this has been tested on 4.7
  // and 4.8 build machines, so it plausibly is.
  #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8) && (__cplusplus < 201103L)
    #undef PRIuPTR
    #undef PRIdPTR
    #define PRIuPTR "lu"
    #define PRIdPTR "ld"
  #endif
  
  #define PRIxPTR2 "08lx"

#endif

#ifndef HAVE_NULLPTR
  #ifndef __cplusplus
    #define nullptr NULL
  #else
    #if __cplusplus <= 199711L
      #ifndef nullptr
        #define nullptr NULL
      #endif
    #endif
  #endif
#endif

// Checked a bunch of alternatives to #define constants.  For integer constants
// less than 2^31, enum {} avoids macro expansion issues that actually matter,
// and that more than cancels out any tiny increase in binary size due to
// additional debugger information (which has value, anyway).  However, we
// don't want to use this under C++ due to enumeral/non-enumeral conditional
// expression warnings, so this isn't one-size-fits-all; and plain old const
// int has all the functionality we want under C++ (including internal linkage,
// so it's fine to define them in header files).  Thus we wrap the
// implementation in a macro.
//
// Otherwise, the macro expansion thing is still annoying but we suck it up due
// to the need for too much duplicate C vs. C++ code ("initializer element is
// not constant" when using const [type] in C99...)
//
// We start most global library-specific numeric constant names here with
// "kPgl", which should have a vanishingly small chance of colliding with
// anything in C99.  Note that stuff like kBytesPerWord is not considered
// library-specific, so it's exempt from having "Pgl" in the name.  Also, the
// few string literals here are of the FOPEN_WB sort, which have similar usage
// patterns to e.g. PRIuPTR which shouldn't be renamed, so those remain
// all-caps.
#ifdef __cplusplus
  #define CONSTU31(name, expr) const uint32_t name = (expr)
#else
  #define CONSTU31(name, expr) enum {name = (expr)}
#endif

// useful because of its bitwise complement: ~k0LU is a word with all 1 bits,
// while ~0 is always 32 1 bits.
// LLU is used over ULL for searchability (no conflict with NULL).
static const uintptr_t k0LU = (uintptr_t)0;

// mainly useful for bitshifts: (k1LU << 32) works in 64-bit builds, while
// (1 << 32) is undefined.  also used as a quicker-to-type way of casting
// numbers/expressions to uintptr_t (via multiplication).
static const uintptr_t k1LU = (uintptr_t)1;


#ifdef __LP64__
  #ifdef USE_AVX2
    CONSTU31(kBytesPerVec, 32);
    CONSTU31(kBytesPerFVec, 32);
    // bleah, have to define these here, vector_size doesn't see enum values
    typedef uintptr_t vul_t __attribute__ ((vector_size (32)));
    typedef float vf_t __attribute__ ((vector_size (32)));
    typedef short vs_t __attribute__ ((vector_size (32)));
    typedef char vc_t __attibute__ ((vector_size (32)));
  #else
    CONSTU31(kBytesPerVec, 16);
    CONSTU31(kBytesPerFVec, 16);
    typedef uintptr_t vul_t __attribute__ ((vector_size (16)));
    typedef float vf_t __attribute__ ((vector_size (16)));
    typedef short vs_t __attribute__ ((vector_size (16)));
    typedef char vc_t __attribute__ ((vector_size (16)));
  #endif
  CONSTU31(kBitsPerWord, 64);
  CONSTU31(kBitsPerWordLog2, 6);

  typedef uint32_t halfword_t;
  typedef uint16_t quarterword_t;

  #ifdef USE_AVX2
    #define VCONST_UL(xx) {xx, xx, xx, xx}
    #define VCONST_S(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
    #define VCONST_C(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
    #define vul_setzero() (vul_t)_mm256_setzero_si256()
    #define vul_rshift(vv, ct) ((vul_t)_mm256_srli_epi64((__m256i)(vv), ct))
    #define vul_lshift(vv, ct) ((vul_t)_mm256_slli_epi64((__m256i)(vv), ct))
  #else
    #define VCONST_UL(xx) {xx, xx}
    #define VCONST_S(xx) {xx, xx, xx, xx, xx, xx, xx, xx}
    #define VCONST_C(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
    // vv = VCONST_UL(k0LU) doesn't work (only ok for initialization)
    #define vul_setzero() (vul_t)_mm_setzero_si128()
    // "vv >> ct" doesn't work, and Scientific Linux gcc 4.4 might not optimize
    // VCONST_UL shift properly (todo: test this)
    #define vul_rshift(vv, ct) ((vul_t)_mm_srli_epi64((__m128i)(vv), ct))
    #define vul_lshift(vv, ct) ((vul_t)_mm_slli_epi64((__m128i)(vv), ct))
  #endif
#else // not __LP64__
  CONSTU31(kBytesPerVec, 4);
  CONSTU31(kBytesPerFVec, 4);
  CONSTU31(kBitsPerWord, 32);
  CONSTU31(kBitsPerWordLog2, 5);

  typedef uint16_t halfword_t;
  typedef uint8_t quarterword_t;

  typedef uintptr_t vul_t;
  typedef float vf_t;
  // vs_t and vc_t aren't worth the trouble of scaling down to 32-bit

  #define VCONST_UL(xx) (xx)
  #define vul_setzero() k0LU
  #define vul_rshift(vv, ct) ((vv) >> (ct))
  #define vul_lshift(vv, ct) ((vv) << (ct))
#endif

static const uintptr_t kMask5555 = (~((uintptr_t)0)) / 3;
static const uintptr_t kMaskAAAA = ((~((uintptr_t)0)) / 3) * 2;
static const uintptr_t kMask3333 = (~((uintptr_t)0)) / 5;
static const uintptr_t kMask1111 = (~((uintptr_t)0)) / 15;
static const uintptr_t kMask0F0F = (~((uintptr_t)0)) / 17;
static const uintptr_t kMask0101 = (~((uintptr_t)0)) / 255;
static const uintptr_t kMask00FF = (~((uintptr_t)0)) / 257;
static const uintptr_t kMask0001 = (~((uintptr_t)0)) / 65535;
static const uintptr_t kMask0000FFFF = (~((uintptr_t)0)) / 65537;
static const uintptr_t kMask00000001 = (~((uintptr_t)0)) / 4294967295U;

static const uintptr_t kMask000000FF = (~((uintptr_t)0)) / 16843009;
static const uintptr_t kMask000F = (~((uintptr_t)0)) / 4369;
static const uintptr_t kMask0303 = (~((uintptr_t)0)) / 85;

CONSTU31(kBitsPerVec, kBytesPerVec * CHAR_BIT);
CONSTU31(kQuatersPerVec, kBytesPerVec * 4);

CONSTU31(kBitsPerWordD2, kBitsPerWord / 2);
CONSTU31(kBitsPerWordD4, kBitsPerWord / 4);

// number of bytes in a word
CONSTU31(kBytesPerWord, kBitsPerWord / CHAR_BIT);

static_assert(CHAR_BIT == 8, "plink2_base requires CHAR_BIT == 8.");
static_assert(sizeof(int32_t) == 4, "plink2_base requires sizeof(int32_t) == 4.");
static_assert(sizeof(int64_t) == 8, "plink2_base requires sizeof(int64_t) == 8.");
static_assert(sizeof(intptr_t) == kBytesPerWord, "plink2_base requires sizeof(intptr_t) == kBytesPerWord.");

CONSTU31(kWordsPerVec, kBytesPerVec / kBytesPerWord);
CONSTU31(kInt32PerVec, kBytesPerVec / 4);

CONSTU31(kCacheline, 64);

CONSTU31(kBitsPerCacheline, kCacheline * CHAR_BIT);
CONSTU31(kQuatersPerCacheline, kCacheline * 4);
CONSTU31(kInt32PerCacheline, kCacheline / sizeof(int32_t));
CONSTU31(kInt64PerCacheline, kCacheline / sizeof(int64_t));
CONSTU31(kWordsPerCacheline, kCacheline / kBytesPerWord);
CONSTU31(kDoublesPerCacheline, kCacheline / sizeof(double));
CONSTU31(kVecsPerCacheline, kCacheline / kBytesPerVec);

// could use ioctl, etc. to dynamically determine this later, and pass it as a
// parameter to e.g. pgfi_multiread
CONSTU31(kDiskBlockSize, 4096);

// unsafe to fread or fwrite more bytes than this on e.g. OS X
CONSTU31(kMaxBytesPerIO, 0x7ffff000);


// note that this is NOT foolproof: see e.g.
// http://insanecoding.blogspot.com/2007/11/pathmax-simply-isnt.html .  (This
// is why I haven't bothered with OS-based #ifdefs here.)  But it should be
// good enough in practice.  And PATH_MAX itself is still relevant due to use
// of realpath().
CONSTU31(kPglFnamesize, 4096);
#if defined(PATH_MAX) && !defined(_WIN32)
static_assert(kPglFnamesize >= PATH_MAX, "plink2_base assumes PATH_MAX <= 4096.  (Safe to increase kPglFnamesize to address this, up to 131072.)");
#endif


typedef union {
  vul_t vi;

  // not actually 8 bytes in 32-bit builds
  uintptr_t u8[kBitsPerVec / kBitsPerWord];

  uint32_t u4[kBytesPerVec / sizeof(int32_t)];
} univec_t;

typedef union {
  vf_t vf;
  float f4[kBytesPerFVec / sizeof(float)];
} univecf_t;

#ifdef __LP64__
// Still need this for access to _mm_mul{lo,hi}_epi16?
typedef union {
  __m128i vi;
  uintptr_t u8[2];
} univec16_t;

// Needed for hand-optimized logistic regression.
typedef union {
  __m128 vf;
  float f4[4];
} univec16f_t;

HEADER_CINLINE uintptr_t univec16_hsum_32bit(univec16_t uv) {
  return ((uv.u8[0] + uv.u8[1]) * kMask00000001) >> 32;
}
#endif

// sum must fit in 16 bits
HEADER_CINLINE uintptr_t univec_hsum_16bit(univec_t uv) {
#ifdef __LP64__
  #ifdef USE_AVX2
  return ((uv.u8[0] + uv.u8[1] + uv.u8[2] + uv.u8[3]) * kMask0001) >> 48;
  #else
  return ((uv.u8[0] + uv.u8[1]) * kMask0001) >> 48;
  #endif
#else
  return (uv.u8[0] * kMask0001) >> 16;
#endif
}

// sum must fit in 32 bits
HEADER_CINLINE uintptr_t univec_hsum_32bit(univec_t uv) {
#ifdef __LP64__
  #ifdef USE_AVX2
  return ((uv.u8[0] + uv.u8[1] + uv.u8[2] + uv.u8[3]) * kMask00000001) >> 32;
  #else
  return ((uv.u8[0] + uv.u8[1]) * kMask00000001) >> 32;
  #endif
#else
  return uv.u8[0];
#endif
}

HEADER_CINLINE2 uintptr_t unpack_halfword_to_word(uintptr_t hw) {
#ifdef __LP64__
  hw = (hw | (hw << 16)) & kMask0000FFFF;
#endif
  hw = (hw | (hw << 8)) & kMask00FF;
  hw = (hw | (hw << 4)) & kMask0F0F;
  hw = (hw | (hw << 2)) & kMask3333;
  return ((hw | (hw << 1)) & kMask5555);
}

HEADER_CINLINE2 halfword_t pack_word_to_halfword(uintptr_t ww) {
  // assumes only even bits of ww can be set
  ww = (ww | (ww >> 1)) & kMask3333;
  ww = (ww | (ww >> 2)) & kMask0F0F;
  ww = (ww | (ww >> 4)) & kMask00FF;
#ifdef __LP64__
  ww = (ww | (ww >> 8)) & kMask0000FFFF;
#endif
  return (halfword_t)(ww | (ww >> kBitsPerWordD4));
}

// alignment must be a power of 2
// tried splitting out round_down_pow2_ui() and _up_pow2_ui() functions, no
// practical difference
HEADER_CINLINE uintptr_t round_down_pow2(uintptr_t val, uintptr_t alignment) {
  return val & (~(alignment - 1));
}

HEADER_CINLINE uint64_t round_down_pow2_ull(uint64_t val, uint64_t alignment) {
  return val & (~(alignment - 1));
}

HEADER_CINLINE uintptr_t round_up_pow2(uintptr_t val, uintptr_t alignment) {
  return (val + alignment - 1) & (~(alignment - 1));
}


// this is best when the divisor is constant (so (divisor - 1) can be
// collapsed), and handles val == 0 properly.  if the divisor isn't constant
// and val is guaranteed to be nonzero, go with explicit
// "1 + (val - 1) / divisor".
//
// Thought about conditional use of constexpr here, but that has annoying
// integer-widening effects.  Unless we split the use cases into DIV_UP,
// DIVL_UP, and DIV64_UP; this may be worth doing at some point.
// Note that this fails if (val + divisor - 1) overflows the widest integer
// type on the left.
#define DIV_UP(val, divisor) (((val) + (divisor) - 1) / (divisor))

// "NZ" means nonzero in two ways:
// * result is in [1, modulus], not [0, modulus - 1]
// * val should not be zero (though this expression still works if val is zero
//   and modulus is a hardcoded power of 2)
#define MOD_NZ(val, modulus) (1 + (((val) - 1) % (modulus)))

HEADER_CINLINE2 uint32_t abs_int32(int32_t ii) {
  const uint32_t neg_sign_bit = -(((uint32_t)ii) >> 31);
  return (((uint32_t)ii) ^ neg_sign_bit) - neg_sign_bit;
}

extern uintptr_t g_failed_alloc_attempt_size;

#if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 7) && !defined(__APPLE__)
// putting this in the header file caused a bunch of gcc 4.4 strict-aliasing
// warnings, while not doing so seems to inhibit some malloc-related compiler
// optimizations, bleah
// compromise: header-inline iff gcc version >= 4.7 (might not be the right
// cutoff?)
boolerr_t pgl_malloc(uintptr_t size, void* pp);
#else
HEADER_INLINE boolerr_t pgl_malloc(uintptr_t size, void* pp) {
  *((unsigned char**)pp) = (unsigned char*)malloc(size);
  if (*((unsigned char**)pp)) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}
#endif

// This must be used for all fwrite() calls where len could be >= 2^31, since
// OS X raw fwrite() doesn't work in that case.
static_assert(sizeof(size_t) == sizeof(intptr_t), "plink2_base assumes size_t and intptr_t are synonymous.");
interr_t fwrite_checked(const void* buf, uintptr_t len, FILE* outfile);

interr_t fread_checked2(void* buf, uintptr_t len, FILE* infile, uintptr_t* bytes_read_ptr);

HEADER_INLINE boolerr_t fread_checked(void* buf, uintptr_t len, FILE* infile) {
  uintptr_t bytes_read;
  if (fread_checked2(buf, len, infile, &bytes_read)) {
    return 1;
  }
  return (bytes_read != len);
}

HEADER_INLINE boolerr_t fclose_null(FILE** fptr_ptr) {
  int32_t ii = ferror(*fptr_ptr);
  int32_t jj = fclose(*fptr_ptr);
  *fptr_ptr = nullptr;
  return ii || jj;
}


#ifdef __LP64__
// Reads an integer in [1, cap].
// * Errors out unless first character is a digit, or is '+' followed by a
//   digit.  Initial whitespace is not permitted.
// * Like atoi(), this considereds the number to be terminated by *any*
//   nondigit character.  E.g. "1000genomes" is treated as a valid instance of
//   1000 rather than a nonnumeric token, and "98.6" is treated as 98.  (See
//   scanadv_posint_capped(), scanadv_uint_capped(), etc. in plink2_common if
//   you want strtol-like semantics, where a pointer to the end of the string
//   is returned.)
// * Errors out on overflow.  This may be the biggest advantage over atoi().
boolerr_t scan_posint_capped(const char* ss, uint64_t cap, uint32_t* valp);

// [0, cap]
boolerr_t scan_uint_capped(const char* ss, uint64_t cap, uint32_t* valp);

// [-bound, bound]
boolerr_t scan_int_abs_bounded(const char* ss, uint64_t bound, int32_t* valp);
#else // not __LP64__
// Need to be more careful in 32-bit case due to overflow.
// A funny-looking div_10/mod_10 interface is used since the cap will usually
// be a constant, and we want the integer division/modulus to occur at compile
// time.
boolerr_t scan_posint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp);

boolerr_t scan_uint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp);

boolerr_t scan_int_abs_bounded32(const char* ss, uint32_t bound_div_10, uint32_t bound_mod_10, int32_t* valp);

HEADER_INLINE boolerr_t scan_posint_capped(const char* ss, uint32_t cap, uint32_t* valp) {
  return scan_posint_capped32(ss, cap / 10, cap % 10, valp);
}

HEADER_INLINE boolerr_t scan_uint_capped(const char* ss, uint32_t cap, uint32_t* valp) {
  return scan_uint_capped32(ss, cap / 10, cap % 10, valp);
}

HEADER_INLINE boolerr_t scan_int_abs_bounded(const char* ss, uint32_t bound, int32_t* valp) {
  return scan_int_abs_bounded32(ss, bound / 10, bound % 10, valp);
}
#endif


// intentionally rejects -2^31 for now
HEADER_INLINE boolerr_t scan_int32(const char* ss, int32_t* valp) {
  return scan_int_abs_bounded(ss, 0x7fffffff, valp);
}

// default cap = 0x7ffffffe
HEADER_INLINE boolerr_t scan_posint_defcap(const char* ss, uint32_t* valp) {
  return scan_posint_capped(ss, 0x7ffffffe, valp);
}

HEADER_INLINE boolerr_t scan_uint_defcap(const char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, 0x7ffffffe, valp);
}

HEADER_INLINE boolerr_t scan_int_abs_defcap(const char* ss, int32_t* valp) {
  return scan_int_abs_bounded(ss, 0x7ffffffe, valp);
}

HEADER_INLINE boolerr_t scan_uint_icap(const char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, 0x7fffffff, valp);
}

HEADER_INLINE unsigned char* memseta(void* target, unsigned char val, uintptr_t ct) {
  memset(target, val, ct);
  return &(((unsigned char*)target)[ct]);
}

HEADER_INLINE char* memcpya(void* __restrict target, const void* __restrict source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(((char*)target)[ct]);
}

#define BITCT_TO_VECCT(val) DIV_UP(val, kBitsPerVec)
#define BITCT_TO_WORDCT(val) DIV_UP(val, kBitsPerWord)
#define BITCT_TO_ALIGNED_WORDCT(val) (kWordsPerVec * BITCT_TO_VECCT(val))
#define BITCT_TO_CLCT(val) DIV_UP(val, kBitsPerCacheline)

// more verbose than (val + 3) / 4, but may as well make semantic meaning
// obvious; any explicit DIV_UP(val, 4) expressions should have a different
// meaning
// (not needed for bitct -> bytect, DIV_UP(val, CHAR_BIT) is clear enough)
#define QUATERCT_TO_BYTECT(val) DIV_UP(val, 4)

#define QUATERCT_TO_VECCT(val) DIV_UP(val, kQuatersPerVec)
#define QUATERCT_TO_WORDCT(val) DIV_UP(val, kBitsPerWordD2)
#define QUATERCT_TO_ALIGNED_WORDCT(val) (kWordsPerVec * QUATERCT_TO_VECCT(val))
#define QUATERCT_TO_CLCT(val) DIV_UP(val, kQuatersPerCacheline)


#define INT32CT_TO_VECCT(val) DIV_UP(val, kInt32PerVec)
#define INT32CT_TO_CLCT(val) DIV_UP(val, kInt32PerCacheline)

#define WORDCT_TO_VECCT(val) DIV_UP(val, kWordsPerVec)
#define WORDCT_TO_CLCT(val) DIV_UP(val, kWordsPerCacheline)

#ifdef __LP64__
  #define INT64CT_TO_VECCT(val) DIV_UP(val, kBytesPerVec / 8)
#else
  #define INT64CT_TO_VECCT(val) ((val) * 2)
#endif
#define INT64CT_TO_CLCT(val) DIV_UP(val, kInt64PerCacheline)
#define DBLCT_TO_VECCT INT64CT_TO_VECCT

#define VECCT_TO_CLCT(val) DIV_UP(val, kVecsPerCacheline)

// C++11 standard guarantees std::min and std::max return leftmost minimum in
// case of equality; best to adhere to that
// We don't actually use std::min/max since casting one argument when comparing
// e.g. a uint32_t with a uintptr_t is pointlessly verbose
#define MAXV(aa, bb) (((bb) > (aa))? (bb) : (aa))
#define MINV(aa, bb) (((bb) < (aa))? (bb) : (aa))

#define GET_QUATERARR_ENTRY(ulptr, idx) (((ulptr)[(idx) / kBitsPerWordD2] >> (2 * ((idx) % kBitsPerWordD2))) & 3)
#define ASSIGN_QUATERARR_ENTRY(idx, newval, ulptr) (ulptr)[(idx) / kBitsPerWordD2] = ((ulptr)[(idx) / kBitsPerWordD2] & (~((3 * k1LU) << (2 * ((idx) % kBitsPerWordD2))))) | (((uintptr_t)(newval)) << (2 * ((idx) % kBitsPerWordD2)))
// todo: check if ASSIGN_QUATERARR_ENTRY optimizes newval=0 out
#define CLEAR_QUATERARR_ENTRY(idx, ulptr) (ulptr)[(idx) / kBitsPerWordD2] &= ~((3 * k1LU) << (idx % kBitsPerWordD2))

#define GET_HEXADECARR_ENTRY(ulptr, idx) (((ulptr)[(idx) / kBitsPerWordD4] >> (4 * ((idx) % kBitsPerWordD4))) & 15)


// don't use pglerr_t here since there's only one failure mode, it's
// obvious what it is, and stacking multiple aligned_mallocs in a single
// if-statement is useful.
boolerr_t aligned_malloc(uintptr_t size, uintptr_t alignment, void* aligned_pp);

// ok for ct == 0
void fill_all_bits(uintptr_t ct, uintptr_t* bitarr);

void bitvec_and(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void bitvec_andnot(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

uint32_t next_set_unsafe(const uintptr_t* bitarr, uint32_t loc);

uint32_t next_unset_unsafe(const uintptr_t* bitarr, uint32_t loc);

// uint32_t next_nonmissing_unsafe(const uintptr_t* genoarr, uint32_t loc);

uint32_t next_set(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil);

uint32_t prev_set_unsafe(const uintptr_t* bitarr, uint32_t loc);

HEADER_INLINE uint32_t are_all_words_zero(const uintptr_t* word_arr, uintptr_t word_ct) {
  while (word_ct--) {
    if (*word_arr++) {
      return 0;
    }
  }
  return 1;
}

HEADER_INLINE uint32_t are_all_bits_one(const uintptr_t* bitarr, uintptr_t bit_ct) {
  const uintptr_t fullword_ct = bit_ct / kBitsPerWord;
  for (uintptr_t widx = 0; widx < fullword_ct; ++widx) {
    if (~(bitarr[widx])) {
      return 0;
    }
  }
  const uint32_t trailing_bit_ct = bit_ct % kBitsPerWord;
  return (!trailing_bit_ct) || ((~(bitarr[fullword_ct])) << (kBitsPerWord - trailing_bit_ct));
}

#ifdef USE_SSE42
HEADER_CINLINE uint32_t popcount2_long(uintptr_t val) {
  return __builtin_popcountll(val) + __builtin_popcountll(val & kMaskAAAA);
}
#else
HEADER_CINLINE2 uint32_t popcount2_long(uintptr_t val) {
  val = (val & kMask3333) + ((val >> 2) & kMask3333);
  return (((val + (val >> 4)) & kMask0F0F) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

  // the simple version, good enough for all non-time-critical stuff
  // (without SSE4.2, popcount_longs() tends to be >3x as fast on arrays.
  // with SSE4.2, there's no noticeable difference.)
#ifdef USE_SSE42
HEADER_CINLINE uint32_t popcount_long(uintptr_t val) {
  return __builtin_popcountll(val);
}
#else
HEADER_CINLINE2 uint32_t popcount_long(uintptr_t val) {
  // sadly, this is still faster than the clang implementation of the intrinsic
  // as of 2016
  return popcount2_long(val - ((val >> 1) & kMask5555));
}
#endif

#ifdef USE_SSE42
HEADER_CINLINE uint32_t popcount_2_longs(uintptr_t val0, uintptr_t val1) {
  return __builtin_popcountll(val0) + __builtin_popcountll(val1);
}
#else
HEADER_CINLINE2 uint32_t popcount_2_longs(uintptr_t val0, uintptr_t val1) {
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  const uintptr_t four_bit = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  // up to 16 values in 0..12; sum fits in 8 bits
  return (((four_bit & kMask0F0F) + ((four_bit >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

#ifndef __LP64__
HEADER_CINLINE2 uint32_t popcount_4_longs(uintptr_t val0, uintptr_t val1, uintptr_t val2, uintptr_t val3) {
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  val2 -= (val2 >> 1) & kMask5555;
  val3 -= (val3 >> 1) & kMask5555;
  const uintptr_t four_bit_0 = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  const uintptr_t four_bit_1 = (val2 & kMask3333) + ((val2 >> 2) & kMask3333) + (val3 & kMask3333) + ((val3 >> 2) & kMask3333);
  return (((four_bit_0 & kMask0F0F) + ((four_bit_0 >> 4) & kMask0F0F) + (four_bit_1 & kMask0F0F) + ((four_bit_1 >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

// assumes vec_ct is a multiple of 3
uintptr_t popcount_vecs(const vul_t* bit_vvec, uintptr_t vec_ct);

#define IS_VEC_ALIGNED(addr) (!(((uintptr_t)(addr)) % kBytesPerVec))

HEADER_INLINE uintptr_t popcount_longs(const uintptr_t* bitvec, uintptr_t word_ct) {
  // Efficiently popcounts bitvec[0..(word_ct - 1)].  In the 64-bit case,
  // bitvec[] must be 16-byte aligned.
  // The popcount_longs_nzbase() wrapper takes care of starting from a later
  // index.
  // No need for a separate USE_SSE42 implementation, there's no noticeable
  // speed difference.
  uintptr_t tot = 0;
  if (word_ct >= (3 * kWordsPerVec)) {
    assert(IS_VEC_ALIGNED(bitvec));
    const uintptr_t remainder = word_ct % (3 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = popcount_vecs((const vul_t*)bitvec, main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    bitvec = &(bitvec[main_block_word_ct]);
  }
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx < word_ct; ++trailing_word_idx) {
    tot += popcount_long(bitvec[trailing_word_idx]);
  }
  return tot;
}

// these don't read past the end of bitarr
uintptr_t popcount_bytes(const unsigned char* bitarr, uintptr_t byte_ct);
uintptr_t popcount_bytes_masked(const unsigned char* bitarr, const uintptr_t* mask_arr, uintptr_t byte_ct);

// requires positive word_ct
// stay agnostic a bit longer re: word_ct := DIV_UP(entry_ct, kBitsPerWord)
// vs. word_ct := 1 + (entry_ct / kBitsPerWord)
// (this is a source of bugs, though; interface should probably be changed to
// use entry_ct once multiallelic/dosage implementation is done)
void fill_cumulative_popcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts);

void uidxs_to_idxs(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, const uint32_t idx_list_len, uint32_t* idx_list);


HEADER_INLINE boolerr_t vecaligned_malloc(uintptr_t size, void* aligned_pp) {
#if defined(__APPLE__) || !defined(__LP64__)
  const boolerr_t ret_boolerr = pgl_malloc(size, aligned_pp);
  assert(IS_VEC_ALIGNED(*((uintptr_t*)aligned_pp)));
  return ret_boolerr;
#else
  return aligned_malloc(size, kBytesPerVec, aligned_pp);
#endif
}

HEADER_INLINE boolerr_t cachealigned_malloc(uintptr_t size, void* aligned_pp) {
  return aligned_malloc(size, kCacheline, aligned_pp);
}

HEADER_INLINE void aligned_free(void* aligned_ptr) {
  free((uintptr_t*)(((uintptr_t*)aligned_ptr)[-1]));
}

HEADER_INLINE void aligned_free_cond(void* aligned_ptr) {
  if (aligned_ptr) {
    free((uintptr_t*)(((uintptr_t*)aligned_ptr)[-1]));
  }
}

// C spec is a bit broken here
HEADER_INLINE void free_const(const void* memptr) {
  // const_cast
  free((void*)((uintptr_t)memptr));
}
 
HEADER_INLINE void free_cond(const void* memptr) {
  if (memptr) {
    free_const(memptr);
  }
}

#if defined(__APPLE__) || !defined(__LP64__)
HEADER_INLINE void vecaligned_free(void* aligned_ptr) {
  free(aligned_ptr);
}

HEADER_INLINE void vecaligned_free_cond(void* aligned_ptr) {
  free_cond(aligned_ptr);
}
#else
HEADER_INLINE void vecaligned_free(void* aligned_ptr) {
  aligned_free(aligned_ptr);
}

HEADER_INLINE void vecaligned_free_cond(void* aligned_ptr) {
  aligned_free_cond(aligned_ptr);
}
#endif

// now compiling with gcc >= 4.4 (or clang equivalent) on all platforms, so
// safe to use memset everywhere
HEADER_INLINE void fill_uint_zero(uintptr_t entry_ct, uint32_t* uiarr) {
  memset(uiarr, 0, entry_ct * sizeof(int32_t));
}

HEADER_INLINE void fill_ulong_zero(uintptr_t entry_ct, uintptr_t* ularr) {
  memset(ularr, 0, entry_ct * sizeof(intptr_t));
}

HEADER_INLINE void fill_ull_zero(uintptr_t entry_ct, uint64_t* ullarr) {
  memset(ullarr, 0, entry_ct * sizeof(int64_t));
}

HEADER_INLINE void fill_ulong_one(uintptr_t entry_ct, uintptr_t* ularr) {
  for (uintptr_t idx = 0; idx < entry_ct; ++idx) {
    ularr[idx] = ~k0LU;
  }
}

#define IS_SET(ulptr, idx) (((ulptr)[(idx) / kBitsPerWord] >> ((idx) % kBitsPerWord)) & 1)

#define SET_BIT(idx, arr) ((arr)[(idx) / kBitsPerWord] |= k1LU << ((idx) % kBitsPerWord))

#define CLEAR_BIT(idx, arr) ((arr)[(idx) / kBitsPerWord] &= ~(k1LU << ((idx) % kBitsPerWord)))

HEADER_INLINE void assign_bit(uintptr_t idx, uintptr_t newbit, uintptr_t* arr) {
  const uintptr_t inv_mask = k1LU << (idx % kBitsPerWord);
  uintptr_t* cur_word_ptr = &(arr[idx / kBitsPerWord]);
  *cur_word_ptr = ((*cur_word_ptr) & (~inv_mask)) | (inv_mask * newbit);
}

HEADER_INLINE void next_set_unsafe_ck(const uintptr_t* __restrict bitarr, uint32_t* __restrict loc_ptr) {
  if (!IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_set_unsafe(bitarr, *loc_ptr);
  }
}

HEADER_INLINE void next_unset_unsafe_ck(const uintptr_t* __restrict bitarr, uint32_t* __restrict loc_ptr) {
  if (IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_unset_unsafe(bitarr, *loc_ptr);
  }
}

// todo: test this against extracting a nonmissing bitarr first
/*
HEADER_INLINE void next_nonmissing_unsafe_ck(const uintptr_t* __restrict genoarr, uint32_t* __restrict loc_ptr) {
  if (GET_QUATERARR_ENTRY(genoarr, *loc_ptr) == 3) {
    *loc_ptr = next_nonmissing_unsafe(genoarr, *loc_ptr);
  }
}
*/

void copy_bitarr_subset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uintptr_t* __restrict output_bitarr);

// Equivalent to popcount_bit_idx(subset_mask, 0, raw_idx).
HEADER_INLINE uint32_t raw_to_subsetted_pos(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, uint32_t raw_idx) {
  // this should be much better than keeping a uidx_to_idx array!
  const uint32_t raw_widx = raw_idx / kBitsPerWord;
  return subset_cumulative_popcounts[raw_widx] + popcount_long(subset_mask[raw_widx] & ((k1LU << (raw_idx % kBitsPerWord)) - k1LU));
}

HEADER_INLINE void zero_trailing_bits(uintptr_t bit_ct, uintptr_t* bitarr) {
  uintptr_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    bitarr[bit_ct / kBitsPerWord] &= (k1LU << trail_ct) - k1LU;
  }
}

HEADER_CINLINE uint32_t bytes_to_represent_ui(uint32_t uii) {
  return (4 - (__builtin_clz(uii) / CHAR_BIT));
}


// transpose_quaterblock(), which is more plink-specific, is in
// pgenlib_internal
CONSTU31(kPglBitTransposeBatch, kBitsPerCacheline);
CONSTU31(kPglBitTransposeWords, kWordsPerCacheline);
CONSTU31(kPglBitTransposeBufbytes, (kPglBitTransposeBatch * kPglBitTransposeBatch) / (CHAR_BIT / 2));
CONSTU31(kPglBitTransposeBufwords, kPglBitTransposeBufbytes / kBytesPerWord);
// up to 512x512; vecaligned_buf must have size 64k
// write_iter must be allocated up to at least
//   round_up_pow2(write_batch_size, 2) rows
void transpose_bitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, vul_t* vecaligned_buf);


// Flagset conventions:
// * Each 32-bit and 64-bit flagset has its own type, which is guaranteed to be
//   the appropriate width.  (Todo: verify that bit 31 works properly in 32-bit
//   case.)
// * Constant flag names start with "kf[CamelCase description]", followed by a
//   description that shouldn't suck too badly.  The zero flagset is always
//   named kf[CamelCase description]0.
// * The type name is always of the form [snake_case description]_flags_t.
// * To gain the desired level of type-checking under C++11 without pointless
//   verbosity, &, |, ^, ~, &=, |=, and ^= operations are defined; [my_flags_t
//   variable] |= [another my_flags_t variable] & [a my_flags_t constant] works
//   without an explicit cast.  (Defining "struct my_flags_t" separately from
//   the enum global-scope-constants container is necessary to make |= work
//   without a cast.  inline is needed due to duplicate operator definitions
//   across multiple files.)
// * To slightly reduce the chance of breakage under C99/C++03, the enum is
//   nameless; the flagset type is just a uint32_t/uint64_t alias.  This is
//   because the C99 and C++03 specs do not provide enough control over the
//   enum base type to make it safe for the enum to serve as the flagset type.
// * Implicit conversion to int is not prevented for now, since I'm trying to
//   keep pglerr_t-style code duplication to a minimum.
#if __cplusplus >= 201103L

  // could avoid the typedef here, but that leads to a bit more verbosity.
  #define FLAGSET_DEF_START() typedef enum : uint32_t {
  #define FLAGSET_DEF_END(tname) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator|(tname ## _PLINK2_BASE_DO_NOT_USE__ aa, tname ## _PLINK2_BASE_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(static_cast<uint32_t>(aa) | static_cast<uint32_t>(bb)); \
} \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator&(tname ## _PLINK2_BASE_DO_NOT_USE__ aa, tname ## _PLINK2_BASE_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(static_cast<uint32_t>(aa) & static_cast<uint32_t>(bb)); \
} \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator^(tname ## _PLINK2_BASE_DO_NOT_USE__ aa, tname ## _PLINK2_BASE_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(static_cast<uint32_t>(aa) ^ static_cast<uint32_t>(bb)); \
} \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator~(tname ## _PLINK2_BASE_DO_NOT_USE__ aa) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(~static_cast<uint32_t>(aa)); \
} \
  \
struct tname { \
  tname() {} \
  \
  tname(const tname& source) : value_(source.value_) {} \
  \
  tname(const tname ## _PLINK2_BASE_DO_NOT_USE__ source) : value_(static_cast<uint32_t>(source)) {} \
  \
  explicit tname(uint32_t source) : value_(source) {} \
  \
  operator tname ## _PLINK2_BASE_DO_NOT_USE__() const { \
    return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(value_); \
  } \
  \
  tname& operator|=(const tname ## _PLINK2_BASE_DO_NOT_USE__ rhs) { \
    value_ |= rhs; \
    return *this; \
  } \
  \
  tname& operator&=(const tname ## _PLINK2_BASE_DO_NOT_USE__ rhs) { \
    value_ &= rhs; \
    return *this; \
  } \
  \
  tname& operator^=(const tname ## _PLINK2_BASE_DO_NOT_USE__ rhs) { \
    value_ ^= rhs; \
    return *this; \
  } \
  \
private: \
  uint32_t value_; \
}

  #define FLAGSET64_DEF_START() typedef enum : uint64_t {
  #define FLAGSET64_DEF_END(tname) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator|(tname ## _PLINK2_BASE_DO_NOT_USE__ aa, tname ## _PLINK2_BASE_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(static_cast<uint64_t>(aa) | static_cast<uint64_t>(bb)); \
} \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator&(tname ## _PLINK2_BASE_DO_NOT_USE__ aa, tname ## _PLINK2_BASE_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(static_cast<uint64_t>(aa) & static_cast<uint64_t>(bb)); \
} \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator^(tname ## _PLINK2_BASE_DO_NOT_USE__ aa, tname ## _PLINK2_BASE_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(static_cast<uint64_t>(aa) ^ static_cast<uint64_t>(bb)); \
} \
  \
inline tname ## _PLINK2_BASE_DO_NOT_USE__ operator~(tname ## _PLINK2_BASE_DO_NOT_USE__ aa) { \
  return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(~static_cast<uint64_t>(aa)); \
} \
  \
struct tname { \
  tname() {} \
  \
  tname(const tname& source) : value_(source.value_) {} \
  \
  tname(const tname ## _PLINK2_BASE_DO_NOT_USE__ source) : value_(static_cast<uint64_t>(source)) {} \
  \
  explicit tname(uint64_t source) : value_(source) {} \
  \
  operator tname ## _PLINK2_BASE_DO_NOT_USE__() const { \
    return static_cast<tname ## _PLINK2_BASE_DO_NOT_USE__>(value_); \
  } \
  \
  tname& operator|=(const tname ## _PLINK2_BASE_DO_NOT_USE__ rhs) { \
    value_ |= rhs; \
    return *this; \
  } \
  \
  tname& operator&=(const tname ## _PLINK2_BASE_DO_NOT_USE__ rhs) { \
    value_ &= rhs; \
    return *this; \
  } \
  \
  tname& operator^=(const tname ## _PLINK2_BASE_DO_NOT_USE__ rhs) { \
    value_ ^= rhs; \
    return *this; \
  } \
  \
private: \
  uint64_t value_; \
}

  #define ENUM_U31_DEF_START() typedef enum : uint32_t {
  #define ENUM_U31_DEF_END(tname) } tname

#else

  #define FLAGSET_DEF_START() enum {
  #define FLAGSET_DEF_END(tname) } ; \
typedef uint32_t tname

  // don't use a nameless enum here, since we want to be able to static_assert
  // the enum size.
  // best to artificially add an element to the end for now to force width to
  // 64-bit, otherwise gcc actually shrinks it even when the constants are
  // defined with LLU.
  #define FLAGSET64_DEF_START() typedef enum {
  #define FLAGSET64_DEF_END(tname) , \
  tname ## PLINK2_BASE_DO_NOT_USE__ALL_64_SET__ = ~(0LLU) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
static_assert(sizeof(tname ## _PLINK2_BASE_DO_NOT_USE__) == 8, "64-bit flagset constants are not actually uint64_ts."); \
typedef uint64_t tname

  #define ENUM_U31_DEF_START() typedef enum {
  #define ENUM_U31_DEF_END(tname) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
typedef uint32_t tname

#endif

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PGENLIB_INTERNAL_H__
