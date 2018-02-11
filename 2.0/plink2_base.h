#ifndef __PLINK2_BASE_H__
#define __PLINK2_BASE_H__

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
//
// Type-choice guidelines:
// - Integers are unsigned by default, signed only when necessary.
//   It's necessary to choose one or the other to avoid drowning in a sea of
//   casts and unexpected behavior.  Each choice has its own surprising
//   pitfalls that the developer had better be aware of; and I definitely do
//   not take the position that unsigned is the better default *in all C and/or
//   C++ code*.  However, for this codebase, the extremely high frequency of
//   bitwise operations makes unsigned-default the only sane choice.
//   Some consequences of this choice:
//   - All pointer differences that are part of a larger arithmetic or
//     comparison expression are explicitly casted to uintptr_t.
//   - Since uint64_t -> double conversion is frequently slower than int64_t ->
//     double conversion, u63tod() should be used when the integer is known to
//     be less than 2^63.  If we also know it's less than 2^31, u31tod() can
//     provide a performance improvement on Win32.
// - Integers that can be >= 2^32 in some of the largest existing datasets, but
//   are usually smaller, should be defined as uintptr_t, to strike a good
//   balance between 32-bit performance and 64-bit scaling.  Exhaustive
//   overflow checking in the 32-bit build is a non-goal; but I do aim for very
//   high statistical reliability, by inserting checks whenever it occurs to me
//   that overflow is especially likely (e.g. when multiplying two potentially
//   large 32-bit numbers).
// - Bitarrays and 'quaterarrays' (packed arrays of 2-bit elements, such as a
//   row of a plink 1.x .bed file) are usually uintptr_t*, to encourage
//   word-at-a-time iteration without requiring vector-alignment.  Quite a few
//   low-level library functions cast them to VecUL*; as mentioned above, the
//   affected function parameter names should end in 'vec' to document the
//   alignment requirements.
// - A buffer/iterator expected to contain only UTF-8 text should be char*.
//   unsigned char* should be reserved for byte-array buffers and iterators
//   which are expected to interact with some non-text bytes, and generic
//   memory-location pointers which will be subject to pointer arithmetic.
//   (Note that this creates some clutter in low-level parsing code: since the
//   signedness of char is platform-dependent, it becomes necessary to use e.g.
//   ctou32() a fair bit.)
// - unsigned char is an acceptable argument type for functions intended to
//   process a single text character, thanks to the automatic cast; it's just
//   unsigned char* that should be avoided.
// - void* return values should be restricted to generic pointers which are
//   *not* expected to be subject to pointer arithmetic.  void* as input
//   parameter type should only be used when there are at least two equally
//   valid input types, NOT counting VecUL*.


// The -Wshorten-64-to-32 diagnostic forces the code to be cluttered with
// meaningless uintptr_t -> uint32_t static casts (number known to be < 2^32,
// just stored in a uintptr_t because there's no speed penalty and we generally
// want to think in terms of word-based operations).  The code is more readable
// if S_CAST(uint32_t, [potentially wider value]) is reserved for situations
// where a higher bit may actually be set.  This pragma can always be commented
// out on the few occasions where inappropriate silent truncation is suspected.
#ifdef __APPLE__
  // todo: explicitly detect clang vs. gcc
#  pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

// 10000 * major + 100 * minor + patch
// Exception to CONSTU31, since we want the preprocessor to have access to this
// value.  Named with all caps as a consequence.
#define PLINK2_BASE_VERNUM 300


#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#  define __STDC_FORMAT_MACROS 1
#endif
#include <inttypes.h>
#include <limits.h> // CHAR_BIT, PATH_MAX

// #define NDEBUG
#include <assert.h>

#ifdef _WIN32
  // needed for MEMORYSTATUSEX
#  ifndef _WIN64
#    define WINVER 0x0500
#  else
#    define __LP64__
#  endif
#  include <windows.h>
#endif

#ifdef __LP64__
#  ifndef __SSE2__
    // todo: remove this requirement, the 32-bit VecUL-using code does most of
    // what we need
#    error "64-bit builds currently require SSE2.  Try producing a 32-bit build instead."
#  endif
#  include <emmintrin.h>
#  ifdef __SSE4_2__
#    define USE_SSE42
#    ifdef __AVX2__
#      include <immintrin.h>
#      ifndef __BMI__
#        error "AVX2 builds require -mbmi as well."
#      endif
#      ifndef __BMI2__
#        error "AVX2 builds require -mbmi2 as well."
#      endif
#      ifndef __LZCNT__
#        error "AVX2 builds require -mlzcnt as well."
#      endif
#      define USE_AVX2
#    endif
#  endif
#endif

#ifndef UINT32_MAX
  // can theoretically be undefined in C++03
#  define UINT32_MAX 0xffffffffU
#endif


// done with #includes, can start C++ namespace
#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef __cplusplus
#  define HEADER_INLINE inline
// Previously went on a wild constexpr spree, but now these are mostly unused.
// Reserve for cases where (i) there's a clear constant-initialization use case
// for an imaginable downstream program (I'm looking at you, DivUp() and
// RoundUpPow2()...), or (ii) it allows a useful static_assert to be inserted
// for a hardcoded constant.
#  if __cplusplus >= 201103L
#    define HEADER_CINLINE constexpr
#    define CSINLINE static constexpr
#    if __cplusplus > 201103L
#      define HEADER_CINLINE2 constexpr
#      define CSINLINE2 static constexpr
#    else
#      define HEADER_CINLINE2 inline
#      define CSINLINE2 static inline
#    endif
#  else
#    define HEADER_CINLINE inline
#    define HEADER_CINLINE2 inline
#    define CSINLINE static inline
#    define CSINLINE2 static inline
#  endif
#  if __cplusplus <= 199711L
    // this may be defined anyway, at least on OS X
#    ifndef static_assert
      // todo: check other cases
#      define static_assert(cond, msg)
#    endif
#  endif
#else
#  define HEADER_INLINE static inline
#  define HEADER_CINLINE static inline
#  define HEADER_CINLINE2 static inline
#  define CSINLINE static inline
#  define CSINLINE2 static inline
  // _Static_assert() should work in gcc 4.6+
#  if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 6)
#    if defined(__APPLE__) && defined(__has_feature) && defined(__has_extension)
      // clang
#      if __has_feature(c_static_assert) || __has_extension(c_static_assert)
#        define static_assert _Static_assert
#      else
#        define static_assert(cond, msg)
#      endif
#    else
#      define static_assert(cond, msg)
#    endif
#  else
#    define static_assert _Static_assert
#  endif
#endif

#define __maybe_unused __attribute__((unused))

#ifdef __cplusplus
#  define K_CAST(type, val) (const_cast<type>(val))
#  define R_CAST(type, val) (reinterpret_cast<type>(val))
#  define S_CAST(type, val) (static_cast<type>(val))
#else
#  define K_CAST(type, val) ((type)(val))
#  define R_CAST(type, val) ((type)(val))
#  define S_CAST(type, val) ((type)(val))
#endif

HEADER_INLINE double u31tod(uint32_t uii) {
  const int32_t ii = uii;
  assert(ii >= 0);
  return S_CAST(double, ii);
}

HEADER_INLINE double u63tod(uint64_t ullii) {
  const int64_t llii = ullii;
  assert(llii >= 0);
  return S_CAST(double, llii);
}

HEADER_INLINE float u31tof(uint32_t uii) {
  const int32_t ii = uii;
  assert(ii >= 0);
  return S_CAST(float, ii);
}

HEADER_INLINE uint32_t ctou32(char cc) {
  return S_CAST(unsigned char, cc);
}

HEADER_INLINE uintptr_t ctoul(char cc) {
  return S_CAST(unsigned char, cc);
}

HEADER_INLINE uint64_t ctou64(char cc) {
  return S_CAST(unsigned char, cc);
}

// Error return types.  All of these evaluate to true on error and false on
// success, but otherwise they have slightly different semantics:
// * PglErr is the general-purpose enum.  Unlike an enum, implicit conversion
//   *to* int, not just from int, is prevented by the C++11 compiler (and the
//   C++11-compiler-validated code still works under C99).  (To achieve this
//   additional safety, we engage in a bit of code duplication which would be
//   unreasonable for flagsets.)
//   Explicit cast to uint32_t, but not int32_t, is supported, to reflect the
//   fact that all error codes are positive.
// * BoolErr allows implicit conversion from int, but conversion back to
//   uint32_t requires an explicit cast.  (It should always be 0/1-valued, but
//   this isn't enforced by the compiler.)
// * IntErr allows implicit conversion from int, but conversion back to
//   int32_t requires an explicit cast.  It mainly serves as a holding pen for
//   C standard library error return values, which can be negative.
#if __cplusplus >= 201103L
struct PglErr {
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

  PglErr() {}

  PglErr(const PglErr& source) : value_(source.value_) {}

  PglErr(ec source) : value_(source) {}

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

const PglErr kPglRetSuccess = PglErr::ec::kPglRetSuccess;
const PglErr kPglRetSkipped = PglErr::ec::kPglRetSkipped;
const PglErr kPglRetNomem = PglErr::ec::kPglRetNomem;
const PglErr kPglRetOpenFail = PglErr::ec::kPglRetOpenFail;
const PglErr kPglRetReadFail = PglErr::ec::kPglRetReadFail;
const PglErr kPglRetWriteFail = PglErr::ec::kPglRetWriteFail;
const PglErr kPglRetMalformedInput = PglErr::ec::kPglRetMalformedInput;
const PglErr kPglRetInconsistentInput = PglErr::ec::kPglRetInconsistentInput;
const PglErr kPglRetInvalidCmdline = PglErr::ec::kPglRetInvalidCmdline;
const PglErr kPglRetHelp = PglErr::ec::kPglRetHelp;
const PglErr kPglRetThreadCreateFail = PglErr::ec::kPglRetThreadCreateFail;
const PglErr kPglRetNetworkFail = PglErr::ec::kPglRetNetworkFail;
const PglErr kPglRetSampleMajorBed = PglErr::ec::kPglRetSampleMajorBed;
const PglErr kPglRetWarningErrcode = PglErr::ec::kPglRetWarningErrcode;
const PglErr kPglRetImproperFunctionCall = PglErr::ec::kPglRetImproperFunctionCall;
const PglErr kPglRetNotYetSupported = PglErr::ec::kPglRetNotYetSupported;
const PglErr kPglRetLongLine = PglErr::ec::kPglRetLongLine;
const PglErr kPglRetEmptyFile = PglErr::ec::kPglRetEmptyFile;
#else
  PglErr;
#endif

#if __cplusplus >= 201103L
// allow efficient arithmetic on these, but force them to require explicit
// int32_t/uint32_t casts; only permit implicit assignment from
// int32_t/uint32_t by default.
// built-in bool type does too many things we don't want...

// expected to be integer-valued, but not necessarily 0/1 or positive
struct IntErr {
  IntErr() {}

  IntErr(int32_t source) : value_(source) {}

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
struct BoolErr {
  BoolErr() {}

  BoolErr(uint32_t source) : value_(source) {}

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
typedef int32_t IntErr;
typedef uint32_t BoolErr;
#endif

// make this work on 32-bit as well as 64-bit systems, across
// Windows/OS X/Linux
// (todo: clean this up a bit.  it's inherently a baling-wire-and-duct-tape
// sort of thing, though...)
#ifdef _WIN32
  // must compile with -std=gnu++11, not c++11, on 32-bit Windows since
  // otherwise fseeko64 not defined...
#  define fseeko fseeko64
#  define ftello ftello64
#  define FOPEN_RB "rb"
#  define FOPEN_WB "wb"
#  define FOPEN_AB "ab"
#  define ferror_unlocked ferror
#  ifdef __LP64__
#    define getc_unlocked _fgetc_nolock
#    define putc_unlocked _fputc_nolock
    // todo: find mingw-w64 build which properly links _fread_nolock, and
    // conditional-compile
#    define fread_unlocked fread
#    define fwrite_unlocked fwrite
#  else
#    define getc_unlocked getc
#    define putc_unlocked putc
#    define fread_unlocked fread
#    define fwrite_unlocked fwrite
#  endif
#  if __cplusplus < 201103L
#    define uint64_t unsigned long long
#    define int64_t long long
#  endif
#else
#  define FOPEN_RB "r"
#  define FOPEN_WB "w"
#  define FOPEN_AB "a"
#  ifdef __APPLE__
#    define fread_unlocked fread
#    define fwrite_unlocked fwrite
#  endif
#endif

#ifdef _WIN32
#  undef PRId64
#  undef PRIu64
#  define PRId64 "I64d"
#  define PRIu64 "I64u"
#else
#  ifdef __cplusplus
#    ifndef PRId64
#      define PRId64 "lld"
#    endif
#  endif
#endif

#ifdef _WIN64
#  define CTZLU __builtin_ctzll
#  define CLZLU __builtin_clzll
#else
#  define CTZLU __builtin_ctzl
#  define CLZLU __builtin_clzl
#  ifndef __LP64__
    // needed to prevent GCC 6 build failure
#    if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
#      if (__cplusplus < 201103L) && !defined(__APPLE__)
#        ifndef uintptr_t
#          define uintptr_t unsigned long
#        endif
#        ifndef intptr_t
#          define intptr_t long
#        endif
#      endif
#    endif
#  endif
#endif

#ifdef __LP64__
#  ifdef _WIN32 // i.e. Win64

#    undef PRIuPTR
#    undef PRIdPTR
#    define PRIuPTR PRIu64
#    define PRIdPTR PRId64
#    define PRIxPTR2 "016I64x"

#  else // not _WIN32

#    ifndef PRIuPTR
#      define PRIuPTR "lu"
#    endif
#    ifndef PRIdPTR
#      define PRIdPTR "ld"
#    endif
#    define PRIxPTR2 "016lx"

#  endif // Win64

#else // not __LP64__

  // without this, we get ridiculous warning spew...
  // not 100% sure this is the right cutoff, but this has been tested on 4.7
  // and 4.8 build machines, so it plausibly is.
#  if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8) && (__cplusplus < 201103L)
#    undef PRIuPTR
#    undef PRIdPTR
#    define PRIuPTR "lu"
#    define PRIdPTR "ld"
#  endif

#  define PRIxPTR2 "08lx"

#endif

#ifndef HAVE_NULLPTR
#  ifndef __cplusplus
#    define nullptr NULL
#  else
#    if __cplusplus <= 199711L
#      ifndef nullptr
#        define nullptr NULL
#      endif
#    endif
#  endif
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
#  define CONSTU31(name, expr) const uint32_t name = (expr)
#else
#  define CONSTU31(name, expr) enum {name = (expr)}
#endif

// useful because of its bitwise complement: ~k0LU is a word with all 1 bits,
// while ~0 is always 32 1 bits.
// LLU is used over ULL for searchability (no conflict with NULL).
static const uintptr_t k0LU = S_CAST(uintptr_t, 0);

// mainly useful for bitshifts: (k1LU << 32) works in 64-bit builds, while
// (1 << 32) is undefined.  also used as a quicker-to-type way of casting
// numbers/expressions to uintptr_t (via multiplication).
static const uintptr_t k1LU = S_CAST(uintptr_t, 1);


#ifdef __LP64__
#  ifdef USE_AVX2
CONSTU31(kBytesPerVec, 32);

// 16 still seems to noticeably outperform 32 on my Mac test machine, and
// is about equal on my Linux test machine, but libraries (and my AVX2
// code...) for the latter should improve over time.  Define FVEC_32 once
// it's time to switch over.
// #define FVEC_32

// bleah, have to define these here, vector_size doesn't see enum values
typedef uintptr_t VecUL __attribute__ ((vector_size (32)));
typedef short VecS __attribute__ ((vector_size (32)));
typedef char VecC __attribute__ ((vector_size (32)));
#  else
CONSTU31(kBytesPerVec, 16);
typedef uintptr_t VecUL __attribute__ ((vector_size (16)));
typedef short VecS __attribute__ ((vector_size (16)));
typedef char VecC __attribute__ ((vector_size (16)));
#  endif
CONSTU31(kBitsPerWord, 64);
CONSTU31(kBitsPerWordLog2, 6);

typedef uint32_t Halfword;
typedef uint16_t Quarterword;

#  ifdef USE_AVX2
// todo: check if _mm256_set1_... makes a difference, and if yes, which
// direction
#    define VCONST_UL(xx) {xx, xx, xx, xx}
#    define VCONST_S(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_C(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
#    define vecul_setzero() R_CAST(VecUL, _mm256_setzero_si256())
#    define vecul_srli(vv, ct) R_CAST(VecUL, _mm256_srli_epi64(R_CAST(__m256i, vv), ct))
#    define vecul_slli(vv, ct) R_CAST(VecUL, _mm256_slli_epi64(R_CAST(__m256i, vv), ct))
#  else
#    define VCONST_UL(xx) {xx, xx}
#    define VCONST_S(xx) {xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_C(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
// vv = VCONST_UL(k0LU) doesn't work (only ok for initialization)
#    define vecul_setzero() R_CAST(VecUL, _mm_setzero_si128())
// "vv >> ct" doesn't work, and Scientific Linux gcc 4.4 might not optimize
// VCONST_UL shift properly (todo: test this)
#    define vecul_srli(vv, ct) R_CAST(VecUL, _mm_srli_epi64(R_CAST(__m128i, vv), ct))
#    define vecul_slli(vv, ct) R_CAST(VecUL, _mm_slli_epi64(R_CAST(__m128i, vv), ct))
#  endif

#  ifdef FVEC_32
#    ifndef __FMA__
#      error "32-byte-float-vector builds require FMA3 as well."
#    endif
CONSTU31(kBytesPerFVec, 32);
typedef float VecF __attribute__ ((vector_size (32)));
#    define VCONST_F(xx) {xx, xx, xx, xx, xx, xx, xx, xx}
#    define vecf_setzero() R_CAST(VecF, _mm256_setzero_ps())
#  else
CONSTU31(kBytesPerFVec, 16);
typedef float VecF __attribute__ ((vector_size (16)));
#    define VCONST_F(xx) {xx, xx, xx, xx}
#    define vecf_setzero() R_CAST(VecF, _mm_setzero_ps())
#  endif
#else // not __LP64__
CONSTU31(kBytesPerVec, 4);
CONSTU31(kBytesPerFVec, 4);
CONSTU31(kBitsPerWord, 32);
CONSTU31(kBitsPerWordLog2, 5);

typedef uint16_t Halfword;
typedef uint8_t Quarterword;

typedef uintptr_t VecUL;
typedef float VecF;
// VecS and VecC aren't worth the trouble of scaling down to 32-bit

#  define VCONST_UL(xx) (xx)
#  define vecul_setzero() k0LU
#  define vecul_srli(vv, ct) ((vv) >> (ct))
#  define vecul_slli(vv, ct) ((vv) << (ct))
#endif

// Unfortunately, we need to spell out S_CAST(uintptr_t, 0) instead of just
// typing k0LU in C99.
static const uintptr_t kMask5555 = (~S_CAST(uintptr_t, 0)) / 3;
static const uintptr_t kMaskAAAA = ((~S_CAST(uintptr_t, 0)) / 3) * 2;
static const uintptr_t kMask3333 = (~S_CAST(uintptr_t, 0)) / 5;
static const uintptr_t kMask1111 = (~S_CAST(uintptr_t, 0)) / 15;
static const uintptr_t kMask0F0F = (~S_CAST(uintptr_t, 0)) / 17;
static const uintptr_t kMask0101 = (~S_CAST(uintptr_t, 0)) / 255;
static const uintptr_t kMask00FF = (~S_CAST(uintptr_t, 0)) / 257;
static const uintptr_t kMask0001 = (~S_CAST(uintptr_t, 0)) / 65535;
static const uintptr_t kMask0000FFFF = (~S_CAST(uintptr_t, 0)) / 65537;
static const uintptr_t kMask00000001 = (~S_CAST(uintptr_t, 0)) / 4294967295U;

static const uintptr_t kMask000000FF = (~S_CAST(uintptr_t, 0)) / 16843009;
static const uintptr_t kMask000F = (~S_CAST(uintptr_t, 0)) / 4369;
static const uintptr_t kMask0303 = (~S_CAST(uintptr_t, 0)) / 85;

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

CONSTU31(kFloatPerFVec, kBytesPerFVec / 4);

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
  VecUL vi;

  // not actually 8 bytes in 32-bit builds, probably want to rename
  uintptr_t u8[kWordsPerVec];

  uint32_t u4[kInt32PerVec];
} UniVec;

typedef union {
  VecF vf;
  float f4[kFloatPerFVec];
} UniVecF;

// sum must fit in 16 bits
HEADER_INLINE uintptr_t UniVecHsum16(UniVec uv) {
#ifdef __LP64__
#  ifdef USE_AVX2
  return ((uv.u8[0] + uv.u8[1] + uv.u8[2] + uv.u8[3]) * kMask0001) >> 48;
#  else
  return ((uv.u8[0] + uv.u8[1]) * kMask0001) >> 48;
#  endif
#else
  return (uv.u8[0] * kMask0001) >> 16;
#endif
}

// sum must fit in 32 bits
HEADER_INLINE uintptr_t UniVecHsum32(UniVec uv) {
#ifdef __LP64__
#  ifdef USE_AVX2
  return ((uv.u8[0] + uv.u8[1] + uv.u8[2] + uv.u8[3]) * kMask00000001) >> 32;
#  else
  return ((uv.u8[0] + uv.u8[1]) * kMask00000001) >> 32;
#  endif
#else
  return uv.u8[0];
#endif
}

HEADER_INLINE float VecFHsum(VecF vecf) {
  UniVecF uvf;
  uvf.vf = vecf;
#ifdef __LP64__
#  ifdef FVEC_32
  // tested various uses of _mm256_hadd_ps, couldn't get them to be faster
  return uvf.f4[0] + uvf.f4[1] + uvf.f4[2] + uvf.f4[3] + uvf.f4[4] + uvf.f4[5] + uvf.f4[6] + uvf.f4[7];
#  else
  return uvf.f4[0] + uvf.f4[1] + uvf.f4[2] + uvf.f4[3];
#  endif
#else
  return uvf.f4[0];
#endif
}

#ifdef USE_AVX2
HEADER_INLINE uintptr_t UnpackHalfwordToWord(uintptr_t hw) {
  return _pdep_u64(hw, kMask5555);
}

HEADER_INLINE Halfword PackWordToHalfword(uintptr_t ww) {
  // Assumes only even bits of ww can be set.
  return _pext_u64(ww, kMask5555);
}
#else // !USE_AVX2
HEADER_INLINE uintptr_t UnpackHalfwordToWord(uintptr_t hw) {
#  ifdef __LP64__
  hw = (hw | (hw << 16)) & kMask0000FFFF;
#  endif
  hw = (hw | (hw << 8)) & kMask00FF;
  hw = (hw | (hw << 4)) & kMask0F0F;
  hw = (hw | (hw << 2)) & kMask3333;
  return ((hw | (hw << 1)) & kMask5555);
}

HEADER_INLINE Halfword PackWordToHalfword(uintptr_t ww) {
  // assumes only even bits of ww can be set
  ww = (ww | (ww >> 1)) & kMask3333;
  ww = (ww | (ww >> 2)) & kMask0F0F;
  ww = (ww | (ww >> 4)) & kMask00FF;
#  ifdef __LP64__
  ww = (ww | (ww >> 8)) & kMask0000FFFF;
#  endif
  return S_CAST(Halfword, ww | (ww >> kBitsPerWordD4));
}
#endif // !USE_AVX2

// alignment must be a power of 2
// tried splitting out RoundDownPow2U32() and RoundUpPow2U32() functions, no
// practical difference
HEADER_CINLINE uintptr_t RoundDownPow2(uintptr_t val, uintptr_t alignment) {
  return val & (~(alignment - 1));
}

HEADER_CINLINE uint64_t RoundDownPow2U64(uint64_t val, uint64_t alignment) {
  return val & (~(alignment - 1));
}

HEADER_CINLINE uintptr_t RoundUpPow2(uintptr_t val, uintptr_t alignment) {
  return (val + alignment - 1) & (~(alignment - 1));
}


// This is best when the divisor is constant (so (divisor - 1) can be
// collapsed), and handles val == 0 properly.  If the divisor isn't constant
// and val is guaranteed to be nonzero, go with explicit
// "1 + (val - 1) / divisor".
//
// Note that this fails if (val + divisor - 1) overflows the widest integer
// type on the left.
//
// Since forced-uint32_t RoundDownPow2 was pointless, it stands to reason that
// the same applies to DivUp.  With that said, we may as well make divisor a
// uint32_t just in case this ever gets used on a not-known-at-compile-time
// divisor, since 64/64 can be slower than 64/32.
HEADER_CINLINE uintptr_t DivUp(uintptr_t val, uint32_t divisor) {
  return (val + divisor - 1) / divisor;
}

HEADER_CINLINE uint64_t DivUpU64(uint64_t val, uint32_t divisor) {
  return (val + divisor - 1) / divisor;
}

// "Nz" means nonzero in two ways:
// * result is in [1, modulus], not [0, modulus - 1]
// * val should not be zero (though this expression still works if val is zero
//   and modulus is a hardcoded power of 2)
HEADER_INLINE uint32_t ModNz(uintptr_t val, uint32_t modulus) {
  return (1 + ((val - 1) % modulus));
}

HEADER_INLINE uint32_t abs_i32(int32_t ii) {
  const uint32_t neg_sign_bit = -(S_CAST(uint32_t, ii) >> 31);
  return (S_CAST(uint32_t, ii) ^ neg_sign_bit) - neg_sign_bit;
}

extern uintptr_t g_failed_alloc_attempt_size;

#if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 7) && !defined(__APPLE__)
// putting this in the header file caused a bunch of gcc 4.4 strict-aliasing
// warnings, while not doing so seems to inhibit some malloc-related compiler
// optimizations, bleah
// compromise: header-inline iff gcc version >= 4.7 (might not be the right
// cutoff?)
BoolErr pgl_malloc(uintptr_t size, void* pp);
#else
// Unfortunately, defining the second parameter to be of type void** doesn't do
// the right thing.
HEADER_INLINE BoolErr pgl_malloc(uintptr_t size, void* pp) {
  *S_CAST(unsigned char**, pp) = S_CAST(unsigned char*, malloc(size));
  if (*S_CAST(unsigned char**, pp)) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}
#endif

// This must be used for all fwrite() calls where len could be >= 2^31, since
// OS X raw fwrite() doesn't work in that case.
static_assert(sizeof(size_t) == sizeof(intptr_t), "plink2_base assumes size_t and intptr_t are synonymous.");
BoolErr fwrite_checked(const void* buf, uintptr_t len, FILE* outfile);

// Only use this if loading < len bytes is not an error.
// IntErr fread_checked2(void* buf, uintptr_t len, FILE* infile, uintptr_t* bytes_read_ptr);

BoolErr fread_checked(void* buf, uintptr_t len, FILE* infile);

HEADER_INLINE BoolErr fclose_null(FILE** fptr_ptr) {
  int32_t ii = ferror_unlocked(*fptr_ptr);
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
//   ScanmovPosintCapped(), ScanmovUintCapped(), etc. in plink2_common if
//   you want strtol-like semantics, where the pointer is advanced.)
// * Errors out on overflow.  This may be the biggest advantage over atoi().
BoolErr ScanPosintCapped(const char* str_iter, uint64_t cap, uint32_t* valp);

// [0, cap]
BoolErr ScanUintCapped(const char* str_iter, uint64_t cap, uint32_t* valp);

// [-bound, bound]
BoolErr ScanIntAbsBounded(const char* str_iter, uint64_t bound, int32_t* valp);
#else // not __LP64__
// Need to be more careful in 32-bit case due to overflow.
// A funny-looking div_10/mod_10 interface is used since the cap will usually
// be a constant, and we want the integer division/modulus to occur at compile
// time.
BoolErr ScanPosintCapped32(const char* str_iter, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp);

BoolErr ScanUintCapped32(const char* str_iter, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp);

BoolErr ScanIntAbsBounded32(const char* str_iter, uint32_t bound_div_10, uint32_t bound_mod_10, int32_t* valp);

HEADER_INLINE BoolErr ScanPosintCapped(const char* str, uint32_t cap, uint32_t* valp) {
  return ScanPosintCapped32(str, cap / 10, cap % 10, valp);
}

HEADER_INLINE BoolErr ScanUintCapped(const char* str, uint32_t cap, uint32_t* valp) {
  return ScanUintCapped32(str, cap / 10, cap % 10, valp);
}

HEADER_INLINE BoolErr ScanIntAbsBounded(const char* str, uint32_t bound, int32_t* valp) {
  return ScanIntAbsBounded32(str, bound / 10, bound % 10, valp);
}
#endif


// intentionally rejects -2^31 for now
// (that's a reason why this doesn't have the shorter name 'ScanI32')
HEADER_INLINE BoolErr ScanInt32(const char* str, int32_t* valp) {
  return ScanIntAbsBounded(str, 0x7fffffff, valp);
}

// default cap = 0x7ffffffe
HEADER_INLINE BoolErr ScanPosintDefcap(const char* str, uint32_t* valp) {
  return ScanPosintCapped(str, 0x7ffffffe, valp);
}

HEADER_INLINE BoolErr ScanUintDefcap(const char* str, uint32_t* valp) {
  return ScanUintCapped(str, 0x7ffffffe, valp);
}

HEADER_INLINE BoolErr ScanIntAbsDefcap(const char* str, int32_t* valp) {
  return ScanIntAbsBounded(str, 0x7ffffffe, valp);
}

HEADER_INLINE BoolErr ScanUintIcap(const char* str, uint32_t* valp) {
  return ScanUintCapped(str, 0x7fffffff, valp);
}


// memcpya() tends to be used to copy known-length text strings, while
// memseta() has more mixed usage but char* type is also at least as common as
// unsigned char*; append comes up less when working with raw byte arrays.  So
// give these char* return types.
HEADER_INLINE char* memseta(void* target, unsigned char val, uintptr_t ct) {
  memset(target, val, ct);
  return &(S_CAST(char*, target)[ct]);
}

HEADER_INLINE char* memcpya(void* __restrict target, const void* __restrict source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(S_CAST(char*, target)[ct]);
}

HEADER_CINLINE uintptr_t BitCtToVecCt(uintptr_t val) {
  return DivUp(val, kBitsPerVec);
}

HEADER_CINLINE uintptr_t BitCtToWordCt(uintptr_t val) {
  return DivUp(val, kBitsPerWord);
}

HEADER_CINLINE uintptr_t BitCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * BitCtToVecCt(val);
}

HEADER_CINLINE uintptr_t BitCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kBitsPerCacheline);
}

HEADER_CINLINE uintptr_t Int32CtToVecCt(uintptr_t val) {
  return DivUp(val, kInt32PerVec);
}

HEADER_CINLINE uintptr_t Int32CtToCachelineCt(uintptr_t val) {
  return DivUp(val, kInt32PerCacheline);
}

HEADER_CINLINE uintptr_t WordCtToVecCt(uintptr_t val) {
  return DivUp(val, kWordsPerVec);
}

HEADER_CINLINE uintptr_t WordCtToCachelineCtU64(uintptr_t val) {
  return DivUpU64(val, kWordsPerCacheline);
}

#ifdef __LP64__
HEADER_CINLINE uintptr_t Int64CtToVecCt(uintptr_t val) {
  return DivUp(val, kBytesPerVec / 8);
}
#else
HEADER_CINLINE uintptr_t Int64CtToVecCt(uintptr_t val) {
  return val * 2;
}
#endif

HEADER_CINLINE uintptr_t Int64CtToCachelineCt(uintptr_t val) {
  return DivUp(val, kInt64PerCacheline);
}

HEADER_CINLINE uintptr_t DblCtToVecCt(uintptr_t val) {
  return Int64CtToVecCt(val);
}

HEADER_CINLINE uintptr_t VecCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kVecsPerCacheline);
}

HEADER_CINLINE uintptr_t VecCtToCachelineCtU64(uintptr_t val) {
  return DivUpU64(val, kVecsPerCacheline);
}

// C++11 standard guarantees std::min and std::max return leftmost minimum in
// case of equality; best to adhere to that
// We don't actually use std::min/max since casting one argument when comparing
// e.g. a uint32_t with a uintptr_t is pointlessly verbose.  Compiler will
// still warn against comparison of signed with unsigned.
#define MAXV(aa, bb) (((bb) > (aa))? (bb) : (aa))
#define MINV(aa, bb) (((bb) < (aa))? (bb) : (aa))


// don't use PglErr here since there's only one failure mode, it's
// obvious what it is, and stacking multiple aligned_mallocs in a single
// if-statement is useful.
BoolErr aligned_malloc(uintptr_t size, uintptr_t alignment, void* aligned_pp);

// ok for ct == 0
void SetAllBits(uintptr_t ct, uintptr_t* bitarr);

void BitvecAnd(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void BitvecAndNot(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

uintptr_t FindFirst1BitFrom(const uintptr_t* bitarr, uintptr_t loc);

uintptr_t FindFirst0BitFrom(const uintptr_t* bitarr, uintptr_t loc);

// uintptr_t NextNonmissingUnsafe(const uintptr_t* genoarr, uintptr_t loc);

uint32_t FindFirst1BitFromBounded(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil);

uint32_t FindLast1BitBefore(const uintptr_t* bitarr, uint32_t loc);

HEADER_INLINE uint32_t AllWordsAreZero(const uintptr_t* word_arr, uintptr_t word_ct) {
  while (word_ct--) {
    if (*word_arr++) {
      return 0;
    }
  }
  return 1;
}

HEADER_INLINE uint32_t AllBitsAreOne(const uintptr_t* bitarr, uintptr_t bit_ct) {
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
HEADER_CINLINE uint32_t QuatersumWord(uintptr_t val) {
  return __builtin_popcountll(val) + __builtin_popcountll(val & kMaskAAAA);
}
#else
HEADER_CINLINE2 uint32_t QuatersumWord(uintptr_t val) {
  val = (val & kMask3333) + ((val >> 2) & kMask3333);
  return (((val + (val >> 4)) & kMask0F0F) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

// the simple version, good enough for all non-time-critical stuff
// (without SSE4.2, PopcountWords() tends to be >3x as fast on arrays.  with
// SSE4.2 but no AVX2, there's no noticeable difference.  with AVX2,
// PopcountWords() gains another factor of 1.5-2x.)
#ifdef USE_SSE42
HEADER_CINLINE uint32_t PopcountWord(uintptr_t val) {
  return __builtin_popcountll(val);
}
#else
HEADER_CINLINE2 uint32_t PopcountWord(uintptr_t val) {
  // Sadly, this was still faster than the clang implementation of the
  // intrinsic as of 2016.
  return QuatersumWord(val - ((val >> 1) & kMask5555));
}
#endif

#ifdef USE_SSE42
HEADER_CINLINE uint32_t Popcount2Words(uintptr_t val0, uintptr_t val1) {
  return __builtin_popcountll(val0) + __builtin_popcountll(val1);
}
#else
HEADER_CINLINE2 uint32_t Popcount2Words(uintptr_t val0, uintptr_t val1) {
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  const uintptr_t four_bit = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  // up to 16 values in 0..12; sum fits in 8 bits
  return (((four_bit & kMask0F0F) + ((four_bit >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

#ifndef __LP64__
HEADER_CINLINE2 uint32_t Popcount4Words(uintptr_t val0, uintptr_t val1, uintptr_t val2, uintptr_t val3) {
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  val2 -= (val2 >> 1) & kMask5555;
  val3 -= (val3 >> 1) & kMask5555;
  const uintptr_t four_bit_0 = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  const uintptr_t four_bit_1 = (val2 & kMask3333) + ((val2 >> 2) & kMask3333) + (val3 & kMask3333) + ((val3 >> 2) & kMask3333);
  return (((four_bit_0 & kMask0F0F) + ((four_bit_0 >> 4) & kMask0F0F) + (four_bit_1 & kMask0F0F) + ((four_bit_1 >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

HEADER_INLINE uint32_t IsVecAligned(const void* ptr) {
  return !(R_CAST(uintptr_t, ptr) % kBytesPerVec);
}

// Updated PopcountWords() code is based on
// https://github.com/kimwalisch/libpopcnt .  libpopcnt license text follows.

/*
 * libpopcnt.h - C/C++ library for counting the number of 1 bits (bit
 * population count) in an array as quickly as possible using
 * specialized CPU instructions i.e. POPCNT, AVX2, AVX512, NEON.
 *
 * Copyright (c) 2016 - 2017, Kim Walisch
 * Copyright (c) 2016 - 2017, Wojciech Mula
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef USE_AVX2
// 'Csa' = carry, save, add
HEADER_INLINE VecUL Csa256(VecUL bb, VecUL cc, VecUL* lp) {
  const VecUL aa = *lp;
  const VecUL uu = aa ^ bb;
  *lp = uu ^ cc;
  return (aa & bb) | (uu & cc);
}

HEADER_INLINE VecUL PopcountVecAvx2(VecUL vv) {
  const __m256i vi = R_CAST(__m256i, vv);
  __m256i lookup1 = _mm256_setr_epi8(
                                     4, 5, 5, 6, 5, 6, 6, 7,
                                     5, 6, 6, 7, 6, 7, 7, 8,
                                     4, 5, 5, 6, 5, 6, 6, 7,
                                     5, 6, 6, 7, 6, 7, 7, 8
                                     );

  __m256i lookup2 = _mm256_setr_epi8(
                                     4, 3, 3, 2, 3, 2, 2, 1,
                                     3, 2, 2, 1, 2, 1, 1, 0,
                                     4, 3, 3, 2, 3, 2, 2, 1,
                                     3, 2, 2, 1, 2, 1, 1, 0
                                     );

  __m256i low_mask = _mm256_set1_epi8(0x0f);
  __m256i lo = _mm256_and_si256(vi, low_mask);
  __m256i hi = _mm256_and_si256(_mm256_srli_epi16(vi, 4), low_mask);
  __m256i popcnt1 = _mm256_shuffle_epi8(lookup1, lo);
  __m256i popcnt2 = _mm256_shuffle_epi8(lookup2, hi);

  return R_CAST(VecUL, _mm256_sad_epu8(popcnt1, popcnt2));
}

HEADER_INLINE uint64_t Hsum64(VecUL vv) {
  UniVec vu;
  vu.vi = vv;
  return vu.u8[0] + vu.u8[1] + vu.u8[2] + vu.u8[3];
  // return _mm256_extract_epi64((__m256i)vv, 0) + _mm256_extract_epi64((__m256i)vv, 1) + _mm256_extract_epi64((__m256i)vv, 2) + _mm256_extract_epi64((__m256i)vv, 3);
}

// assumes vec_ct is a multiple of 16
uintptr_t PopcountVecsAvx2(const VecUL* bit_vvec, uintptr_t vec_ct);

HEADER_INLINE uintptr_t PopcountWords(const uintptr_t* bitvec, uintptr_t word_ct) {
  // Efficiently popcounts bitvec[0..(word_ct - 1)].  In the 64-bit case,
  // bitvec[] must be 16-byte aligned.
  // The PopcountWordsNzbase() wrapper takes care of starting from a later
  // index.
  uintptr_t tot = 0;
  if (word_ct >= (16 * kWordsPerVec)) {
    assert(IsVecAligned(bitvec));
    const uintptr_t remainder = word_ct % (16 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = PopcountVecsAvx2(R_CAST(const VecUL*, bitvec), main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    bitvec = &(bitvec[main_block_word_ct]);
  }
  // todo: check if libpopcnt manual-4x-unroll makes a difference on any test
  // machine (I'd prefer to trust the compiler to take care of that...)
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx < word_ct; ++trailing_word_idx) {
    tot += PopcountWord(bitvec[trailing_word_idx]);
  }
  return tot;
}
#else // !USE_AVX2
// assumes vec_ct is a multiple of 3
uintptr_t PopcountVecsNoSse42(const VecUL* bit_vvec, uintptr_t vec_ct);

HEADER_INLINE uintptr_t PopcountWords(const uintptr_t* bitvec, uintptr_t word_ct) {
  uintptr_t tot = 0;
#  ifndef USE_SSE42
  if (word_ct >= (3 * kWordsPerVec)) {
    assert(IsVecAligned(bitvec));
    const uintptr_t remainder = word_ct % (3 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = PopcountVecsNoSse42(R_CAST(const VecUL*, bitvec), main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    bitvec = &(bitvec[main_block_word_ct]);
  }
#  endif
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx < word_ct; ++trailing_word_idx) {
    tot += PopcountWord(bitvec[trailing_word_idx]);
  }
  return tot;
}
#endif // !USE_AVX2

// Turns out memcpy(&cur_word, bytearr, ct) can't be trusted to be fast when ct
// isn't known at compile time.
//
// ct must be less than sizeof(intptr_t).  ct == 0 handled correctly, albeit
// inefficiently.
HEADER_INLINE uintptr_t ProperSubwordLoad(const void* bytearr, uint32_t ct) {
  const unsigned char* bytearr_iter = S_CAST(const unsigned char*, bytearr);
  bytearr_iter = &(bytearr_iter[ct]);
  uintptr_t cur_word = 0;
  if (ct & 1) {
    cur_word = *(--bytearr_iter);
  }
#ifdef __LP64__
  if (ct & 2) {
    cur_word <<= 16;
    bytearr_iter = &(bytearr_iter[-2]);
    cur_word |= *R_CAST(const uint16_t*, bytearr_iter);
  }
  if (ct & 4) {
    cur_word <<= 32;
    cur_word |= *S_CAST(const uint32_t*, bytearr);
  }
#else
  if (ct & 2) {
    cur_word <<= 16;
    cur_word |= *S_CAST(const uint16_t*, bytearr);
  }
#endif
  return cur_word;
}

HEADER_INLINE uintptr_t SubwordLoad(const void* bytearr, uint32_t ct) {
  if (ct == kBytesPerWord) {
    return *S_CAST(const uintptr_t*, bytearr);
  }
  return ProperSubwordLoad(bytearr, ct);
}

// ct must be in 1..4.
HEADER_INLINE uint32_t SubUintLoad(const void* bytearr, uint32_t ct) {
  if (ct & 1) {
    const unsigned char* bytearr_iter = S_CAST(const unsigned char*, bytearr);
    uint32_t cur_uint = *bytearr_iter;
    if (ct == 3) {
      ++bytearr_iter;
      cur_uint |= S_CAST(uint32_t, *R_CAST(const uint16_t*, bytearr_iter)) << 8;
    }
    return cur_uint;
  }
  if (ct == 2) {
    return *S_CAST(const uint16_t*, bytearr);
  }
  return *S_CAST(const uint32_t*, bytearr);
}

// tried making this non-inline, loop took more than 50% longer
HEADER_INLINE void ProperSubwordStore(uintptr_t cur_word, uint32_t byte_ct, void* target) {
  unsigned char* target_iter = S_CAST(unsigned char*, target);
  if (byte_ct & 1) {
    *target_iter++ = cur_word;
    cur_word >>= 8;
  }
#ifdef __LP64__
  if (byte_ct & 2) {
    *R_CAST(uint16_t*, target_iter) = cur_word;
    cur_word >>= 16;
    target_iter = &(target_iter[2]);
  }
  if (byte_ct & 4) {
    *R_CAST(uint32_t*, target_iter) = S_CAST(uint32_t, cur_word);
  }
#else
  if (byte_ct & 2) {
    *R_CAST(uint16_t*, target_iter) = cur_word;
  }
#endif
}

HEADER_INLINE void ProperSubwordStoreMov(uintptr_t cur_word, uint32_t byte_ct, unsigned char** targetp) {
  ProperSubwordStore(cur_word, byte_ct, *targetp);
  *targetp += byte_ct;
}

HEADER_INLINE void SubwordStore(uintptr_t cur_word, uint32_t byte_ct, void* target) {
  if (byte_ct == kBytesPerWord) {
    *S_CAST(uintptr_t*, target) = cur_word;
    return;
  }
  ProperSubwordStore(cur_word, byte_ct, target);
}

HEADER_INLINE void SubwordStoreMov(uintptr_t cur_word, uint32_t byte_ct, unsigned char** targetp) {
  SubwordStore(cur_word, byte_ct, *targetp);
  *targetp += byte_ct;
}

// byte_ct must be in 1..4.
HEADER_INLINE void SubUintStore(uint32_t cur_uint, uint32_t byte_ct, void* target) {
  if (byte_ct & 1) {
    unsigned char* target_iter = S_CAST(unsigned char*, target);
    *target_iter = cur_uint;
    if (byte_ct == 3) {
      ++target_iter;
      *R_CAST(uint16_t*, target_iter) = cur_uint >> 8;
    }
    return;
  }
  if (byte_ct == 2) {
    *S_CAST(uint16_t*, target) = cur_uint;
    return;
  }
  *S_CAST(uint32_t*, target) = cur_uint;
  return;
}

HEADER_INLINE void SubUintStoreMov(uint32_t cur_uint, uint32_t byte_ct, unsigned char** targetp) {
  SubUintStore(cur_uint, byte_ct, *targetp);
  *targetp += byte_ct;
}

// requires positive word_ct
// stay agnostic a bit longer re: word_ct := DIV_UP(entry_ct, kBitsPerWord)
// vs. word_ct := 1 + (entry_ct / kBitsPerWord)
// (this is a source of bugs, though; interface should probably be changed to
// use entry_ct once multiallelic/dosage implementation is done)
void FillCumulativePopcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts);

// If idx_list is a list of valid unfiltered indexes, this converts them
// in-place to corresponding filtered indexes.
void UidxsToIdxs(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, const uint32_t idx_list_len, uint32_t* idx_list);


HEADER_INLINE BoolErr vecaligned_malloc(uintptr_t size, void* aligned_pp) {
#ifdef USE_AVX2
  return aligned_malloc(size, kBytesPerVec, aligned_pp);
#else
#  if defined(__APPLE__) || !defined(__LP64__)
  const BoolErr ret_boolerr = pgl_malloc(size, aligned_pp);
  assert(IsVecAligned(*S_CAST(uintptr_t**, aligned_pp)));
  return ret_boolerr;
#  else
  return aligned_malloc(size, kBytesPerVec, aligned_pp);
#  endif
#endif
}

HEADER_INLINE BoolErr cachealigned_malloc(uintptr_t size, void* aligned_pp) {
  return aligned_malloc(size, kCacheline, aligned_pp);
}

HEADER_INLINE void aligned_free(void* aligned_ptr) {
  free(R_CAST(void*, S_CAST(uintptr_t*, aligned_ptr)[-1]));
}

HEADER_INLINE void aligned_free_cond(void* aligned_ptr) {
  if (aligned_ptr) {
    free(R_CAST(void*, S_CAST(uintptr_t*, aligned_ptr)[-1]));
  }
}

// C spec is slightly broken here
HEADER_INLINE void free_const(const void* memptr) {
  free(K_CAST(void*, memptr));
}

HEADER_INLINE void free_cond(const void* memptr) {
  if (memptr) {
    free_const(memptr);
  }
}

#ifdef USE_AVX2
HEADER_INLINE void vecaligned_free(void* aligned_ptr) {
  aligned_free(aligned_ptr);
}

HEADER_INLINE void vecaligned_free_cond(void* aligned_ptr) {
  aligned_free_cond(aligned_ptr);
}
#else
#  if defined(__APPLE__) || !defined(__LP64__)
HEADER_INLINE void vecaligned_free(void* aligned_ptr) {
  free(aligned_ptr);
}

HEADER_INLINE void vecaligned_free_cond(void* aligned_ptr) {
  free_cond(aligned_ptr);
}
#  else
HEADER_INLINE void vecaligned_free(void* aligned_ptr) {
  aligned_free(aligned_ptr);
}

HEADER_INLINE void vecaligned_free_cond(void* aligned_ptr) {
  aligned_free_cond(aligned_ptr);
}
#  endif
#endif

// now compiling with gcc >= 4.4 (or clang equivalent) on all platforms, so
// safe to use memset everywhere
HEADER_INLINE void ZeroUiArr(uintptr_t entry_ct, uint32_t* uiarr) {
  memset(uiarr, 0, entry_ct * sizeof(int32_t));
}

HEADER_INLINE void ZeroUlArr(uintptr_t entry_ct, uintptr_t* ularr) {
  memset(ularr, 0, entry_ct * sizeof(intptr_t));
}

HEADER_INLINE void ZeroU64Arr(uintptr_t entry_ct, uint64_t* ullarr) {
  memset(ullarr, 0, entry_ct * sizeof(int64_t));
}

HEADER_INLINE void SetAllUlArr(uintptr_t entry_ct, uintptr_t* ularr) {
  // todo: test this against memset(, 255, ) and manually vectorized loop
  for (uintptr_t idx = 0; idx < entry_ct; ++idx) {
    ularr[idx] = ~k0LU;
  }
}

HEADER_INLINE uintptr_t IsSet(const uintptr_t* bitarr, uintptr_t idx) {
  return (bitarr[idx / kBitsPerWord] >> (idx % kBitsPerWord)) & 1;
}

HEADER_INLINE void SetBit(uintptr_t idx, uintptr_t* bitarr) {
  bitarr[idx / kBitsPerWord] |= k1LU << (idx % kBitsPerWord);
}

HEADER_INLINE void ClearBit(uintptr_t idx, uintptr_t* bitarr) {
  bitarr[idx / kBitsPerWord] &= ~(k1LU << (idx % kBitsPerWord));
}

HEADER_INLINE void AssignBit(uintptr_t idx, uintptr_t newbit, uintptr_t* bitarr) {
  const uintptr_t inv_mask = k1LU << (idx % kBitsPerWord);
  uintptr_t* cur_word_ptr = &(bitarr[idx / kBitsPerWord]);
  *cur_word_ptr = ((*cur_word_ptr) & (~inv_mask)) | (inv_mask * newbit);
}

HEADER_INLINE void FindFirst1BitFromU32(const uintptr_t* __restrict bitarr, uint32_t* __restrict loc_ptr) {
  if (!IsSet(bitarr, *loc_ptr)) {
    *loc_ptr = FindFirst1BitFrom(bitarr, *loc_ptr);
  }
}

HEADER_INLINE void FindFirst0BitFromU32(const uintptr_t* __restrict bitarr, uint32_t* __restrict loc_ptr) {
  if (IsSet(bitarr, *loc_ptr)) {
    *loc_ptr = FindFirst0BitFrom(bitarr, *loc_ptr);
  }
}

// todo: test this against extracting a nonmissing bitarr first
/*
HEADER_INLINE void NextNonmissingUnsafeCk32(const uintptr_t* __restrict genoarr, uint32_t* __restrict loc_ptr) {
  if (GetQuaterarrEntry(genoarr, *loc_ptr) == 3) {
    *loc_ptr = NextNonmissingUnsafe(genoarr, *loc_ptr);
  }
}
*/

// tried _bzhi_u64() in AVX2 case, it was actually worse on my Mac (more
// opaque to compiler?)
// This is undefined if idx == kBitsPerWord.
HEADER_INLINE uintptr_t bzhi(uintptr_t ww, uint32_t idx) {
  return ww & ((k1LU << idx) - k1LU);
}

// This is undefined if idx == 0.
HEADER_INLINE uintptr_t bzhi_max(uintptr_t ww, uint32_t idx) {
  return ww & ((~k0LU) >> (kBitsPerWord - idx));
}

// Don't bother defining blsr(), compiler will automatically use the
// instruction under -mbmi and regular code is more readable.

// Equivalent to popcount_bit_idx(subset_mask, 0, raw_idx).
HEADER_INLINE uint32_t RawToSubsettedPos(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, uint32_t raw_idx) {
  // this should be much better than keeping a uidx_to_idx array!
  // (update: there are more compact indexes, but postpone for now, this is
  // is nice and simple and gets us most of what we need.)
  const uint32_t raw_widx = raw_idx / kBitsPerWord;
  return subset_cumulative_popcounts[raw_widx] + PopcountWord(bzhi(subset_mask[raw_widx], raw_idx % kBitsPerWord));
}

HEADER_INLINE void ZeroTrailingBits(uintptr_t bit_ct, uintptr_t* bitarr) {
  const uint32_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    bitarr[bit_ct / kBitsPerWord] = bzhi(bitarr[bit_ct / kBitsPerWord], trail_ct);
  }
}

// requires nonzero uii
HEADER_CINLINE uint32_t BytesToRepresentUi(uint32_t uii) {
  return (4 - (__builtin_clz(uii) / CHAR_BIT));
}

#ifdef __LP64__
HEADER_INLINE void ZeroTrailingWords(uint32_t word_ct, uintptr_t* bitvec) {
  const uint32_t remainder = word_ct % kWordsPerVec;
  if (remainder) {
    ZeroUlArr(kWordsPerVec - remainder, &(bitvec[word_ct]));
  }
}
#else
HEADER_INLINE void ZeroTrailingWords(__maybe_unused uint32_t word_ct, __maybe_unused uintptr_t* bitvec) {
}
#endif

void CopyBitarrSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uintptr_t* __restrict output_bitarr);

// expand_size + read_start_bit must be positive.
void ExpandBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, uint32_t word_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target);

// equivalent to calling ExpandBytearr() followed by CopyBitarrSubset()
void ExpandThenSubsetBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, const uintptr_t* __restrict subset_mask, uint32_t expand_size, uint32_t subset_size, uint32_t read_start_bit, uintptr_t* __restrict target);

// mid_popcount must be positive
void ExpandBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, uint32_t word_ct, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target);

// mid_popcount must be positive
// if mid_start_bit == 1, mid_popcount should not include that bit
void ExpandThenSubsetBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target);

// these don't read past the end of bitarr
uintptr_t PopcountBytes(const unsigned char* bitarr, uintptr_t byte_ct);
uintptr_t PopcountBytesMasked(const unsigned char* bitarr, const uintptr_t* mask_arr, uintptr_t byte_ct);


// transpose_quaterblock(), which is more plink-specific, is in
// pgenlib_internal
CONSTU31(kPglBitTransposeBatch, kBitsPerCacheline);
CONSTU31(kPglBitTransposeWords, kWordsPerCacheline);
// * Up to 512x512; vecaligned_buf must have size 32k (64k in 32-bit case)
//   (er, can reduce buffer size to 512 bytes in 64-bit case...)
// * write_iter must be allocated up to at least
//   RoundUpPow2(write_batch_size, 2) rows
// * We use pointers with different types to read from and write to buf0/buf1,
//   so defining the base type as unsigned char* is theoretically necessary to
//   avoid breaking strict-aliasing rules, while the restrict qualifiers should
//   tell the compiler it doesn't need to be paranoid about writes to one of
//   the buffers screwing with reads from the other.
#ifdef __LP64__
CONSTU31(kPglBitTransposeBufbytes, kPglBitTransposeBatch);
void TransposeBitblockInternal(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, void* buf0);

HEADER_INLINE void TransposeBitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecUL* vecaligned_buf) {
  TransposeBitblockInternal(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, vecaligned_buf);
}

#else // !__LP64__
CONSTU31(kPglBitTransposeBufbytes, (kPglBitTransposeBatch * kPglBitTransposeBatch) / (CHAR_BIT / 2));
void TransposeBitblockInternal(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecUL* __restrict buf0, VecUL* __restrict buf1);

// If this ever needs to be called on an input byte array, read_iter could be
// changed to const void*; in that case, read_ul_stride should be changed to a
// byte count.
HEADER_INLINE void TransposeBitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecUL* vecaligned_buf) {
  TransposeBitblockInternal(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, vecaligned_buf, &(vecaligned_buf[kPglBitTransposeBufbytes / (2 * kBytesPerWord)]));
}
#endif

CONSTU31(kPglBitTransposeBufwords, kPglBitTransposeBufbytes / kBytesPerWord);

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
//   keep PglErr-style code duplication to a minimum.
#if __cplusplus >= 201103L

  // could avoid the typedef here, but that leads to a bit more verbosity.
#  define FLAGSET_DEF_START() typedef enum : uint32_t {
#  define FLAGSET_DEF_END(tname) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
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

#  define FLAGSET64_DEF_START() typedef enum : uint64_t {
#  define FLAGSET64_DEF_END(tname) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
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

#  define ENUM_U31_DEF_START() typedef enum : uint32_t {
#  define ENUM_U31_DEF_END(tname) } tname

#else

#  define FLAGSET_DEF_START() enum {
#  define FLAGSET_DEF_END(tname) } ; \
typedef uint32_t tname

  // don't use a nameless enum here, since we want to be able to static_assert
  // the enum size.
  // best to artificially add an element to the end for now to force width to
  // 64-bit, otherwise gcc actually shrinks it even when the constants are
  // defined with LLU.
#  define FLAGSET64_DEF_START() typedef enum {
#  define FLAGSET64_DEF_END(tname) , \
  tname ## PLINK2_BASE_DO_NOT_USE__ALL_64_SET__ = ~(0LLU) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
static_assert(sizeof(tname ## _PLINK2_BASE_DO_NOT_USE__) == 8, "64-bit flagset constants are not actually uint64_ts."); \
typedef uint64_t tname

#  define ENUM_U31_DEF_START() typedef enum {
#  define ENUM_U31_DEF_END(tname) } tname ## _PLINK2_BASE_DO_NOT_USE__ ; \
typedef uint32_t tname

#endif

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PGENLIB_INTERNAL_H__
