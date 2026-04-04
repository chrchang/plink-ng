#ifndef __PLINK2_BASE_H__
#define __PLINK2_BASE_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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


// Low-level C99/C++03/C++11 library covering basic I/O, SWAR, and
// Windows/macOS/Linux portability.  We try to benefit from as much C++ type
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
//   - Actual SIMD code has been moved into plink2_simd, to make this a lighter
//     dependency.
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
// - Bitarrays and 'nyparrays' (packed arrays of 2-bit elements, such as a row
//   of a plink 1.x .bed file) are usually uintptr_t*, to encourage
//   word-at-a-time iteration without requiring vector-alignment.  Quite a few
//   low-level library functions cast them to VecW*.  As mentioned above, the
//   affected function parameter names must end in 'vec' when this creates an
//   alignment requirement.
// - A buffer/iterator expected to contain only UTF-8 text should be char*.
//   unsigned char* should be reserved for byte-array buffers and iterators
//   which are expected to interact with some non-text bytes, and generic
//   memory-location pointers which will be subject to pointer arithmetic.
//   (Note that this creates some clutter in low-level parsing code: since the
//   signedness of char is platform-dependent, it becomes necessary to use e.g.
//   ctou32() a fair bit.)
// - unsigned char is an acceptable argument type for functions intended to
//   process a single text character, thanks to implicit char -> unsigned char
//   conversion; it's just unsigned char* that should be avoided.
// - void* return values should be restricted to generic pointers which are
//   *not* expected to be subject to pointer arithmetic.  void* as input
//   parameter type should only be used when there are at least two equally
//   valid input types, NOT counting VecW*.


// gcc 8.3.0 has been miscompiling the vector code in the ParseOnebitUnsafe()
// function in pgenlib_read.cc for the last several years.  gcc 8.4 does not
// have this problem, and neither does any other gcc major version I've tested
// to date.
#if !defined(__clang__) && defined(__LP64__)
#  if (__GNUC__ == 8) && (__GNUC_MINOR__ < 4)
#    error "gcc 8.3 is known to have a miscompilation bug that was fixed in 8.4."
#  endif
#endif

#ifndef __clang__
#  if (__GNUC__ < 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ < 9))
// May eventually add MSVC support to gain access to MKL on Windows, but can't
// justify doing that before all major features are implemented.
// We're currently building with a version of libdeflate requiring 4.9+; may as
// well drop the remaining bits of our own code requiring an earlier version.
#    error "gcc 4.9+ or clang equivalent required."
#  endif
#endif

// The -Wshorten-64-to-32 diagnostic forces the code to be cluttered with
// meaningless uintptr_t -> uint32_t static casts (number known to be < 2^32,
// just stored in a uintptr_t because there's no speed penalty and we generally
// want to think in terms of word-based operations).  The code is more readable
// if S_CAST(uint32_t, <potentially wider value>) is reserved for situations
// where a higher bit may actually be set.  This pragma can always be commented
// out on the few occasions where inappropriate silent truncation is suspected.
#ifdef __clang__
#  pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

// 10000 * major + 100 * minor + patch
// Exception to CONSTI32, since we want the preprocessor to have access
// to this value.  Named with all caps as a consequence.
#define PLINK2_BASE_VERNUM 900

// We now try to adhere to include-what-you-use in simple cases.  However,
// we don't want to repeat either platform-specific ifdefs, or stuff like
// "#define _FILE_OFFSET_BITS 64" before "#include <stdio.h>", so we use
// "IWYU pragma: export" for those cases.
// TODO: unfortunately, there is no straightforward way to tell IWYU to stop
// looking behind the likes of <stdio.h> or "../simde/x86/sse2.h".  I'm
// manually ignoring the associated false-positives for now, but this should be
// systematized.

#define _FILE_OFFSET_BITS 64
#ifdef __USE_MINGW_ANSI_STDIO
#  undef __USE_MINGW_ANSI_STDIO
#endif
#define __USE_MINGW_ANSI_STDIO 1

#include <assert.h>
#ifndef __STDC_FORMAT_MACROS
#  define __STDC_FORMAT_MACROS 1
#endif
#include <inttypes.h>  // IWYU pragma: export
#include <limits.h>  // CHAR_BIT, PATH_MAX
#include <stddef.h>  // offsetof()
#include <stdio.h>  // IWYU pragma: export
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
  // needed for EnterCriticalSection, etc.
#  ifndef _WIN64
#    define WINVER 0x0501
#  else
#    define __LP64__
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>  // IWYU pragma: export
#endif

#ifdef __cplusplus
#  include <algorithm>  // IWYU pragma: export
#  if __cplusplus >= 201103L
#    include <array>  // IWYU pragma: export
// Disable for pgenlibr for now, adding -ltbb might be enough but I want to
// test that outside pgenlibr build first.
#    if __cplusplus >= 201902L && defined(__GNUC__) && !defined(NO_UNALIGNED) && !defined(__clang__)
#      include <execution>  // IWYU pragma: export
#    endif
#  endif
#endif

#ifdef __x86_64__
#  include <xmmintrin.h>  // IWYU pragma: export
#endif

// We define USE_SSE42 and USE_AVX2 here, instead of plink2_simd, since
// x86-64-v2 includes efficient single-word popcount, and x86-64-v3 has a
// bit-manipulation side.
#ifdef __SSE4_2__
// probable todo: rename this to something like USE_X86_64_V2 (and similar
// change for USE_AVX2)
#  define USE_SSE42
#  ifdef __AVX2__
#    if defined(__BMI__) && defined(__BMI2__) && defined(__LZCNT__)
#      include <immintrin.h>
#      define USE_AVX2
#    else
// Graceful downgrade, in case -march=native misfires on a VM.  See
// https://github.com/chrchang/plink-ng/issues/155 .
#      warning "AVX2 builds require -mbmi, -mbmi2, and -mlzcnt as well.  Downgrading to SSE4.2 build."
#      undef USE_AVX2
#    endif
#  endif
#endif

// done with #includes, can start C++ namespace...
#ifdef __cplusplus
namespace plink2 {
#endif

// ...though a bunch of symbols remain to be #defined; try to reduce the number
// over time.

#ifndef UINT32_MAX
  // can theoretically be undefined in C++03
#  define UINT32_MAX 0xffffffffU
#endif

#define UINT32_MAXM1 0xfffffffeU

#ifndef UINT64_MAX
#  define UINT64_MAX 0xffffffffffffffffLLU
#endif

#ifdef __cplusplus
#  define HEADER_INLINE inline
// Previously went on a wild constexpr spree, but now these are mostly unused.
// Reserve for cases where (i) there's a clear constant-initialization use case
// for an imaginable downstream program (I'm looking at you, DivUp() and
// RoundUpPow2()...), or (ii) it allows a useful static_assert to be inserted
// for a hardcoded constant.
#  if __cplusplus >= 201103L
#    define HAS_CONSTEXPR
#    define PREFER_CONSTEXPR constexpr
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
#    define PREFER_CONSTEXPR const
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
#  define PREFER_CONSTEXPR const
#  define HEADER_INLINE static inline
#  define HEADER_CINLINE static inline
#  define HEADER_CINLINE2 static inline
#  define CSINLINE static inline
#  define CSINLINE2 static inline
#  define static_assert _Static_assert
#endif

#define __maybe_unused __attribute__((unused))

// Rule of thumb: Use these macros if, and only if, the condition would always
// trigger exit-from-program.  As a side effect, this makes it more
// straightforward, if still tedious, to make global changes to error-handling
// strategy (always dump backtrace and exit immediately?), though provision
// must still be made for sometimes-error-sometimes-not return paths which
// don't get an unlikely annotation.
#ifndef likely
#  define likely(expr) __builtin_expect(!!(expr), 1)
#  define unlikely(expr) __builtin_expect(!!(expr), 0)
#endif

#ifdef __cplusplus
#  define K_CAST(type, val) (const_cast<type>(val))
#  define R_CAST(type, val) (reinterpret_cast<type>(val))
#  define S_CAST(type, val) (static_cast<type>(val))
#else
#  define K_CAST(type, val) ((type)(val))
#  define R_CAST(type, val) ((type)(val))
#  define S_CAST(type, val) ((type)(val))
#endif

// (from Linux kernel)
// container_of - cast a member of a structure out to the containing structure
// @ptr: the pointer to the member.
// @type: the type of the container struct this is embedded in.
// @member: the name of the member within the struct.
#define container_of(ptr, type, member) \
  (R_CAST(type*, R_CAST(char*, ptr) - offsetof(type, member)))

// original macro doesn't work in C++ when e.g. ptr is a const char*, and the
// quick workaround of casting away the const is unsafe.
#define const_container_of(ptr, type, member) \
  (R_CAST(const type*, R_CAST(const char*, ptr) - offsetof(type, member)))

HEADER_INLINE double u31tod(uint32_t uii) {
  const int32_t ii = uii;
  assert(ii >= 0);
  return S_CAST(double, ii);
}

HEADER_INLINE double swtod(intptr_t lii) {
  return S_CAST(double, lii);
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

HEADER_INLINE uintptr_t ctow(char cc) {
  return S_CAST(unsigned char, cc);
}

HEADER_INLINE uint64_t ctou64(char cc) {
  return S_CAST(unsigned char, cc);
}

HEADER_INLINE unsigned char* CToUc(char* pp) {
  return R_CAST(unsigned char*, pp);
}

// Error return types.  All of these evaluate to true on error and false on
// success, but otherwise they have slightly different semantics:
// * PglErr is the general-purpose enum.  Unlike an enum, implicit conversion
//   *to* int, not just from int, is prevented by the C++11 compiler (and the
//   C++11-compiler-validated code still works under C99).  (To achieve this
//   additional safety, we engage in a bit of code duplication which would be
//   unreasonable for flagsets.)
//   (Previously, explicit cast to uint32_t, but not int32_t, was supported, to
//   reflect the fact that all error codes are positive.  This was deemed
//   silly.)
// * BoolErr allows implicit conversion from int, but conversion back to
//   uint32_t requires an explicit cast.  (It should always be 0/1-valued, but
//   this isn't enforced by the compiler.)
// * IntErr allows implicit conversion from int, but conversion back to
//   int32_t requires an explicit cast.  It mainly serves as a holding pen for
//   C standard library error return values, which can be negative.
#if __cplusplus >= 201103L && !defined(NO_CPP11_TYPE_ENFORCEMENT)
// Cython doesn't support these, and the previous workaround of forcing C++98
// is breaking down.  So allow this C++11 feature to be selectively disabled.
#  define CPP11_TYPE_ENFORCEMENT
#endif

#ifdef CPP11_TYPE_ENFORCEMENT
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
  // problems like mismatched files (generally solvable by fixing the command
  // line), and DegenerateData for properly-formatted-and-matched files that
  // yields degenerate computational results due to e.g. divide by zero or
  // insufficient rank.
  kPglRetMalformedInput,
  kPglRetInconsistentInput,
  kPglRetInvalidCmdline,
  kPglRetThreadCreateFail,
  kPglRetNetworkFail,
  kPglRetVarRecordTooLarge,
  kPglRetUnsupportedInstructions,
  kPglRetDegenerateData,
  kPglRetDecompressFail, // also distinguish this from MalformedInput
  kPglRetRewindFail,
  kPglRetGpuFail,
  kPglRetInternalUse1,
  kPglRetSampleMajorBed = 32,
  kPglRetNomemCustomMsg = 59,
  kPglRetInternalError = 60,
  kPglRetWarningErrcode = 61,
  kPglRetImproperFunctionCall = 62,
  kPglRetNotYetSupported = 63,

  // These are only for internal use.  If any of these reach the top level
  // instead of being handled or converted to another error code, that's a bug,
  // and plink2 prints a message to that effect.
  kPglRetRetry = 124,
  kPglRetHelp = 125,
  kPglRetLongLine = 126,
  kPglRetEof = 127}
#ifdef CPP11_TYPE_ENFORCEMENT
  ;

  PglErr() {}

  PglErr(const PglErr& source) : value_(source.value_) {}

  PglErr(ec source) : value_(source) {}

  // Allow explicit conversion from uint64_t, and NOT uint32_t, to this error
  // type, to support reproducible multithreaded error reporting (where
  // multiple threads may atomically attempt to modify a single uint64_t with
  // the error code in the low 32 bits and a priority number in the high bits).
  explicit PglErr(uint64_t source) : value_(static_cast<ec>(source)) {}

  PglErr& operator=(const PglErr&) = default;

  operator ec() const {
    return value_;
  }

  explicit operator uint32_t() const {
    return static_cast<uint32_t>(value_);
  }

  explicit operator int32_t() const {
    return static_cast<int32_t>(value_);
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
const PglErr kPglRetRetry = PglErr::ec::kPglRetRetry;
const PglErr kPglRetHelp = PglErr::ec::kPglRetHelp;
const PglErr kPglRetThreadCreateFail = PglErr::ec::kPglRetThreadCreateFail;
const PglErr kPglRetNetworkFail = PglErr::ec::kPglRetNetworkFail;
const PglErr kPglRetVarRecordTooLarge = PglErr::ec::kPglRetVarRecordTooLarge;
const PglErr kPglRetUnsupportedInstructions = PglErr::ec::kPglRetUnsupportedInstructions;
const PglErr kPglRetDegenerateData = PglErr::ec::kPglRetDegenerateData;
const PglErr kPglRetDecompressFail = PglErr::ec::kPglRetDecompressFail;
const PglErr kPglRetRewindFail = PglErr::ec::kPglRetRewindFail;
const PglErr kPglRetGpuFail = PglErr::ec::kPglRetGpuFail;
const PglErr kPglRetInternalUse1 = PglErr::ec::kPglRetInternalUse1;
const PglErr kPglRetSampleMajorBed = PglErr::ec::kPglRetSampleMajorBed;
const PglErr kPglRetWarningErrcode = PglErr::ec::kPglRetWarningErrcode;
const PglErr kPglRetNomemCustomMsg = PglErr::ec::kPglRetNomemCustomMsg;
const PglErr kPglRetInternalError = PglErr::ec::kPglRetInternalError;
const PglErr kPglRetImproperFunctionCall = PglErr::ec::kPglRetImproperFunctionCall;
const PglErr kPglRetNotYetSupported = PglErr::ec::kPglRetNotYetSupported;
const PglErr kPglRetLongLine = PglErr::ec::kPglRetLongLine;
const PglErr kPglRetEof = PglErr::ec::kPglRetEof;
#else
  PglErr;
#endif

#ifdef CPP11_TYPE_ENFORCEMENT
// allow efficient arithmetic on these, but force them to require explicit
// int32_t/uint32_t casts; only permit implicit assignment from
// int32_t/uint32_t by default.
// built-in bool type does too many things we don't want...

// expected to be integer-valued, but not necessarily 0/1 or positive
struct IntErr {
  IntErr() {}

  IntErr(int32_t source) : value_(source) {}

  explicit operator int32_t() const {
    return value_;
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
    return value_;
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

extern const char kErrprintfFopen[];
extern const char kErrprintfFread[];
extern const char kErrprintfRewind[];
extern const char kErrstrNomem[];
extern const char kErrstrWrite[];
extern const char kErrprintfDecompress[];

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
#  define feof_unlocked feof
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
#else  // Linux or OS X
#  define FOPEN_RB "r"
#  define FOPEN_WB "w"
#  define FOPEN_AB "a"
#  if defined(__APPLE__) || defined(__FreeBSD__) || defined(__NetBSD__)
#    define fread_unlocked fread
#    define fwrite_unlocked fwrite
#  endif
#  if defined(__NetBSD__)
#    define ferror_unlocked ferror
#    define feof_unlocked feof
#  endif
#endif

#ifdef __cplusplus
#  ifndef PRId64
#    define PRId64 "lld"
#  endif
#endif

// These are useful for defending against base-pointer integer overflow on bad
// input.
HEADER_INLINE BoolErr PtrAddCk(const unsigned char* end, intptr_t incr, const unsigned char** basep) {
  *basep += incr;
  return unlikely((end - (*basep)) < 0);
}

// 'W' for writable
HEADER_INLINE BoolErr PtrWSubCk(unsigned char* base, intptr_t decr, unsigned char** endp) {
  *endp -= decr;
  return unlikely(((*endp) - base) < 0);
}

HEADER_INLINE BoolErr PtrCheck(const void* end, const void* base, intptr_t req) {
  const unsigned char* end_uc = S_CAST(const unsigned char*, end);
  const unsigned char* base_uc = S_CAST(const unsigned char*, base);
  return unlikely((end_uc - base_uc) < req);
}

// We want this to return an uint32_t, not an int32_t.
HEADER_INLINE uint32_t ctzu32(uint32_t uii) {
  return __builtin_ctz(uii);
}

// this should always compile down to bsr.
HEADER_INLINE uint32_t bsru32(uint32_t uii) {
  return 31 - __builtin_clz(uii);
}

#ifdef _WIN64
HEADER_INLINE uint32_t ctzw(unsigned long long ullii) {
  return __builtin_ctzll(ullii);
}

HEADER_INLINE uint32_t bsrw(unsigned long long ullii) {
  // Note that this actually compiles to a single instruction on x86; it's
  // naked __builtin_clzll which requires an additional subtraction.
  return 63 - __builtin_clzll(ullii);
}
#else
HEADER_INLINE uint32_t ctzw(unsigned long ulii) {
  return __builtin_ctzl(ulii);
}

HEADER_INLINE uint32_t bsrw(unsigned long ulii) {
  return (8 * sizeof(intptr_t) - 1) - __builtin_clzl(ulii);
}
#endif

#ifdef __LP64__
#  ifndef PRIuPTR
#    define PRIuPTR "lu"
#  endif
#  ifndef PRIdPTR
#    define PRIdPTR "ld"
#  endif
#  define PRIxPTR2 "016lx"
#else  // not __LP64__

#  ifdef _WIN32
#    ifndef PRIuPTR
#      define PRIuPTR "u"
#    endif
#    ifndef PRIdPTR
#      define PRIdPTR "d"
#    endif
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
// in [-2^31, 2^31), enum {} avoids macro expansion issues that actually
// matter, and that more than cancels out any tiny increase in binary size due
// to additional debugger information (which has value, anyway).  However, we
// don't want to use this under C++ due to enumeral/non-enumeral conditional
// expression warnings, so this isn't one-size-fits-all; and plain old const
// int has all the functionality we want under C++ (including internal linkage,
// so it's fine to define them in header files).  Thus we wrap the
// implementation in a macro.
//
// Otherwise, the macro expansion thing is still annoying but we suck it up due
// to the need for too much duplicate C vs. C++ code ("initializer element is
// not constant" when using const <type> in C99...)
//
// We start most plink2- and pgenlib-specific numeric constant names here with
// "kPgl", which should have a vanishingly small chance of colliding with
// anything in C99.  Note that stuff like kBytesPerWord is not considered
// library-specific, so it's exempt from having "Pgl" in the name.  Also, the
// few string literals here are of the FOPEN_WB sort, which have similar usage
// patterns to e.g. PRIuPTR which shouldn't be renamed, so those remain
// all-caps.
//
// (Update, May 2018: CONSTU31 was renamed to CONSTI32 and changed to type
// int32_t, to prevent C vs. C++ differences.  This almost never makes a
// difference, since when int32_t and uint32_t are mixed in the same
// expression, the former gets converted to unsigned.  However, unpleasant
// surprises are occasionally possible when mixing these constants with
// uint16_t or unsigned char values, since then the unsigned values are
// promoted to int32_t.  Also, this essentially forces use of -Wno-sign-compare
// when using gcc 4.4.
//
// Biggest thing to watch out for is mixing of Halfword with these constants in
// 32-bit builds.  Dosage and Vec8thUint are also relevant.)
#ifdef __cplusplus
#  define CONSTI32(name, expr) const int32_t name = (expr)
#else
#  define CONSTI32(name, expr) enum {name = (expr)}
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
CONSTI32(kBitsPerWord, 64);
CONSTI32(kBitsPerWordLog2, 6);

typedef uint32_t Halfword;
typedef uint16_t Quarterword;
#else
CONSTI32(kBitsPerWord, 32);
CONSTI32(kBitsPerWordLog2, 5);

typedef uint16_t Halfword;
typedef uint8_t Quarterword;
#endif

HEADER_INLINE uint32_t* I32ToU32(int32_t* pp) {
  return R_CAST(uint32_t*, pp);
}

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

CONSTI32(kBitsPerWordD2, kBitsPerWord / 2);
CONSTI32(kBitsPerWordD4, kBitsPerWord / 4);

// number of bytes in a word
CONSTI32(kBytesPerWord, kBitsPerWord / CHAR_BIT);

CONSTI32(kInt16PerWord, kBytesPerWord / 2);

static_assert(CHAR_BIT == 8, "plink2_base requires CHAR_BIT == 8.");
static_assert(sizeof(int8_t) == 1, "plink2_base requires sizeof(int8_t) == 1.");
static_assert(sizeof(int16_t) == 2, "plink2_base requires sizeof(int16_t) == 2.");
static_assert(sizeof(int32_t) == 4, "plink2_base requires sizeof(int32_t) == 4.");
static_assert(sizeof(int) >= 4, "plink2_base requires sizeof(int) >= 4.");
static_assert(sizeof(intptr_t) == kBytesPerWord, "plink2_base requires sizeof(intptr_t) == kBytesPerWord.");
static_assert(sizeof(int64_t) == 8, "plink2_base requires sizeof(int64_t) == 8.");

#if defined(__APPLE__) && defined(__LP64__) && !defined(__x86_64__)
#  define CACHELINE128
CONSTI32(kCacheline, 128);
#else
#  define CACHELINE64
CONSTI32(kCacheline, 64);
#endif

CONSTI32(kBitsPerCacheline, kCacheline * CHAR_BIT);
CONSTI32(kNypsPerCacheline, kCacheline * 4);
CONSTI32(kInt16PerCacheline, kCacheline / sizeof(int16_t));
CONSTI32(kInt32PerCacheline, kCacheline / sizeof(int32_t));
CONSTI32(kInt64PerCacheline, kCacheline / sizeof(int64_t));
CONSTI32(kWordsPerCacheline, kCacheline / kBytesPerWord);
CONSTI32(kDoublesPerCacheline, kCacheline / sizeof(double));

// could use ioctl, etc. to dynamically determine this later, and pass it as a
// parameter to e.g. PgfiMultiread()
CONSTI32(kDiskBlockSize, 4096);

CONSTI32(kPglFwriteBlockSize, 131072);

// unsafe to fread or fwrite more bytes than this on e.g. OS X
CONSTI32(kMaxBytesPerIO, 0x7ffff000);

// Maximum size of "dynamically" allocated line load buffer.  (This is the
// limit that applies to .vcf and similar files.)  Inconvenient to go higher
// since fgets() takes a int32_t size argument.
#if defined(__APPLE__) && defined(__LP64__) && !defined(__x86_64__)
CONSTI32(kMaxLongLine, 0x7fffff80);
#else
CONSTI32(kMaxLongLine, 0x7fffffc0);
#endif
static_assert(!(kMaxLongLine % kCacheline), "kMaxLongLine must be a multiple of kCacheline.");

#ifdef __APPLE__
// macOS is limited to 256?
CONSTI32(kMaxOpenFiles, 252);
#else
// Can't assume more than 512 are allowed on Windows, with current compilation
// settings.
CONSTI32(kMaxOpenFiles, 504);
#endif

// note that this is NOT foolproof: see e.g.
// http://insanecoding.blogspot.com/2007/11/pathmax-simply-isnt.html .  (This
// is why I haven't bothered with OS-based #ifdefs here.)  But it should be
// good enough in practice.  And PATH_MAX itself is still relevant due to use
// of realpath().
CONSTI32(kPglFnamesize, 4096);
#if defined(PATH_MAX) && !defined(_WIN32)
static_assert(kPglFnamesize >= PATH_MAX, "plink2_base assumes PATH_MAX <= 4096.  (Safe to increase kPglFnamesize to address this, up to 131072.)");
#endif

// safe errstr_buf size for PgenInitPhase{1,2}(), PgrValidate(),
// BitmapReaderInitPhase{1,2}()
CONSTI32(kPglErrstrBufBlen, kPglFnamesize + 256);

// shared between .pgen and plink-bitmap formats
// currently must be power of 2, and multiple of (kBitsPerWord / 2)
CONSTI32(kPglDifflistGroupSize, 64);

#ifdef CPP11_TYPE_ENFORCEMENT
// Main application of std::array in this codebase is enforcing length when
// passing references between functions.  Conversely, if the array type has
// different lengths in different functions (e.g. col_skips[]/col_types[]), we
// actively want to avoid &arr[0] clutter.
// When neither applies, it doesn't really matter whether we use this or not;
// I normally won't use it unless it plausibly makes sense to pass
// fixed-length-array-references in the future.
#  define STD_ARRAY_DECL(tt, nn, vv) std::array<tt, nn> vv
#  define STD_ARRAY_REF(tt, nn) std::array<tt, nn>&

// necessary if tt is a pointer type, otherwise optional
#  define STD_ARRAY_KREF(tt, nn) const std::array<tt, nn>&

#  define STD_ARRAY_COPY(src, nn, dst) static_assert(sizeof(dst) == sizeof((dst)[0]) * nn, "Invalid STD_ARRAY_COPY() invocation."); (dst) = (src)

#  define STD_ARRAY_PTR_TYPE(tt, nn) std::array<tt, nn>*
#  define STD_ARRAY_PTR_DECL(tt, nn, vv) std::array<tt, nn>* vv

// argh, need double-braces for C++11 std::array and single-braces for C
#  define STD_ARRAY_INIT_START() {
#  define STD_ARRAY_INIT_END() }

template <class T, std::size_t N> void STD_ARRAY_FILL0(std::array<T, N>& arr) {
  arr.fill(0);
}

// plain STD_ARRAY_FILL0() can't be used on array-references due to fallback
// code.
// this macro ensures that we *only* use it with uint32_t array-references
#  define STD_ARRAY_REF_FILL0(ct, aref) static_assert(ct * sizeof(aref[0]) == sizeof(aref), "invalid STD_ARRAY_REF_FILL0() invocation"); aref.fill(0)

// reasonable to default to applying one of the following two restrictions to
// pointer- and flexible-array-member-containing structs, though I haven't done
// so exhaustively.
#  define NONCOPYABLE(struct_name) \
  struct_name() = default; \
  struct_name(const struct_name&) = delete; \
  struct_name& operator=(const struct_name&) = delete

#  define MOVABLE_BUT_NONCOPYABLE(struct_name) \
  struct_name() = default; \
  struct_name(const struct_name&) = delete; \
  struct_name& operator=(const struct_name&) = delete; \
  struct_name(struct_name&&) = default; \
  struct_name& operator=(struct_name&&) = default

#else
#  define STD_ARRAY_DECL(tt, nn, vv) tt vv[nn]
#  define STD_ARRAY_REF(tt, nn) tt* const
#  define STD_ARRAY_KREF(tt, nn) tt const* const
#  define STD_ARRAY_COPY(src, nn, dst) memcpy(dst, src, nn * sizeof(dst[0]));
#  define STD_ARRAY_PTR_TYPE(tt, nn) tt(*)[nn]
#  define STD_ARRAY_PTR_DECL(tt, nn, vv) tt(*vv)[nn]
#  define STD_ARRAY_INIT_START()
#  define STD_ARRAY_INIT_END()
#  define STD_ARRAY_FILL0(arr) memset(arr, 0, sizeof(arr))
#  define STD_ARRAY_REF_FILL0(ct, aref) memset(aref, 0, ct * sizeof(*aref))

#  define NONCOPYABLE(struct_name)
#  define MOVABLE_BUT_NONCOPYABLE(struct_name)
#endif

#ifdef USE_AVX2
HEADER_INLINE uintptr_t UnpackHalfwordToWord(uintptr_t hw) {
  return _pdep_u64(hw, kMask5555);
}

HEADER_INLINE uintptr_t UnpackHalfwordToWordShift1(uintptr_t hw) {
  return _pdep_u64(hw, kMaskAAAA);
}

HEADER_INLINE Halfword PackWordToHalfword(uintptr_t ww) {
  return _pext_u64(ww, kMask5555);
}

HEADER_INLINE Halfword PackWordToHalfwordMask5555(uintptr_t ww) {
  return _pext_u64(ww, kMask5555);
}

HEADER_INLINE Halfword PackWordToHalfwordMaskAAAA(uintptr_t ww) {
  return _pext_u64(ww, kMaskAAAA);
}

HEADER_INLINE uintptr_t PackTwo5555Mask(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = PackWordToHalfwordMask5555(lo);
  const uintptr_t hi_packed = _pext_u64(hi, kMask5555);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}

HEADER_INLINE Halfword Pack3333(uintptr_t ww) {
  return _pext_u64(ww, kMask3333);
}

HEADER_INLINE Halfword Pack3333Mask(uintptr_t ww) {
  return _pext_u64(ww, kMask3333);
}

HEADER_INLINE uintptr_t PackTwo3333Mask(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = Pack3333Mask(lo);
  const uintptr_t hi_packed = _pext_u64(hi, kMask3333);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}

HEADER_INLINE uintptr_t Unpack0F0F(uintptr_t hw) {
  return _pdep_u64(hw, kMask0F0F);
}

HEADER_INLINE Halfword Pack0F0F(uintptr_t ww) {
  return _pext_u64(ww, kMask0F0F);
}

HEADER_INLINE Halfword Pack0F0FMask(uintptr_t ww) {
  return _pext_u64(ww, kMask0F0F);
}

HEADER_INLINE uintptr_t PackTwo0F0F(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = _pext_u64(lo, kMask0F0F);
  const uintptr_t hi_packed = _pext_u64(hi, kMask0F0F);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}

HEADER_INLINE uintptr_t PackTwo0F0FMask(uintptr_t lo, uintptr_t hi) {
  return PackTwo0F0F(lo, hi);
}

HEADER_INLINE uintptr_t Unpack0303(uintptr_t qw) {
  return _pdep_u64(qw, kMask0303);
}

HEADER_INLINE Quarterword Pack0303(uintptr_t ww) {
  return _pext_u64(ww, kMask0303);
}

HEADER_INLINE Quarterword Pack0303Mask(uintptr_t ww) {
  return _pext_u64(ww, kMask0303);
}
#else  // !USE_AVX2
HEADER_INLINE uintptr_t UnpackHalfwordToWord(uintptr_t hw) {
#  ifdef __LP64__
  hw = (hw | (hw << 16)) & kMask0000FFFF;
#  endif
  hw = (hw | (hw << 8)) & kMask00FF;
  hw = (hw | (hw << 4)) & kMask0F0F;
  hw = (hw | (hw << 2)) & kMask3333;
  return ((hw | (hw << 1)) & kMask5555);
}

HEADER_INLINE uintptr_t UnpackHalfwordToWordShift1(uintptr_t hw) {
  return UnpackHalfwordToWord(hw) << 1;
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

HEADER_INLINE Halfword PackWordToHalfwordMask5555(uintptr_t ww) {
  return PackWordToHalfword(ww & kMask5555);
}

HEADER_INLINE Halfword PackWordToHalfwordMaskAAAA(uintptr_t ww) {
  return PackWordToHalfword((ww >> 1) & kMask5555);
}

HEADER_INLINE uintptr_t PackTwo5555Mask(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = PackWordToHalfwordMask5555(lo);

  // Avoid Halfword cast here.
  hi = hi & kMask5555;
  hi = (hi | (hi >> 1)) & kMask3333;
  hi = (hi | (hi >> 2)) & kMask0F0F;
  hi = (hi | (hi >> 4)) & kMask00FF;
#  ifdef __LP64__
  hi = (hi | (hi >> 8)) & kMask0000FFFF;
#  endif
  const uintptr_t hi_packed = hi | (hi >> kBitsPerWordD4);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}

HEADER_INLINE Halfword Pack3333(uintptr_t ww) {
  ww = (ww | (ww >> 2)) & kMask0F0F;
  ww = (ww | (ww >> 4)) & kMask00FF;
#  ifdef __LP64__
  ww = (ww | (ww >> 8)) & kMask0000FFFF;
#  endif
  return S_CAST(Halfword, ww | (ww >> kBitsPerWordD4));
}

HEADER_INLINE Halfword Pack3333Mask(uintptr_t ww) {
  return Pack3333(ww & kMask3333);
}

HEADER_INLINE uintptr_t PackTwo3333Mask(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = Pack3333Mask(lo);

  // Avoid Halfword cast here.
  hi = hi & kMask3333;
  hi = (hi | (hi >> 2)) & kMask0F0F;
  hi = (hi | (hi >> 4)) & kMask00FF;
#  ifdef __LP64__
  hi = (hi | (hi >> 8)) & kMask0000FFFF;
#  endif
  const uintptr_t hi_packed = hi | (hi >> kBitsPerWordD4);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}

HEADER_INLINE uintptr_t Unpack0F0F(uintptr_t hw) {
#  ifdef __LP64__
  hw = (hw | (hw << 16)) & kMask0000FFFF;
#  endif
  hw = (hw | (hw << 8)) & kMask00FF;
  return ((hw | (hw << 4)) & kMask0F0F);
}

HEADER_INLINE Halfword Pack0F0F(uintptr_t ww) {
  ww = (ww | (ww >> 4)) & kMask00FF;
#  ifdef __LP64__
  ww = (ww | (ww >> 8)) & kMask0000FFFF;
#  endif
  return S_CAST(Halfword, ww | (ww >> kBitsPerWordD4));
}

HEADER_INLINE Halfword Pack0F0FMask(uintptr_t ww) {
  return Pack0F0F(ww & kMask0F0F);
}

HEADER_INLINE uintptr_t PackTwo0F0F(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = Pack0F0F(lo);

  // Avoid Halfword cast here.
  hi = (hi | (hi >> 4)) & kMask00FF;
#  ifdef __LP64__
  hi = (hi | (hi >> 8)) & kMask0000FFFF;
#  endif
  const uintptr_t hi_packed = hi | (hi >> kBitsPerWordD4);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}

HEADER_INLINE uintptr_t PackTwo0F0FMask(uintptr_t lo, uintptr_t hi) {
  const Halfword lo_packed = Pack0F0FMask(lo);

  // Avoid Halfword cast here.
  hi = hi & kMask0F0F;
  hi = (hi | (hi >> 4)) & kMask00FF;
#  ifdef __LP64__
  hi = (hi | (hi >> 8)) & kMask0000FFFF;
#  endif
  const uintptr_t hi_packed = hi | (hi >> kBitsPerWordD4);
  return lo_packed | (hi_packed << kBitsPerWordD2);
}


HEADER_INLINE uintptr_t Unpack0303(uintptr_t qw) {
  // ................................................fedcba9876543210
#  ifdef __LP64__
  qw = (qw | (qw << 24)) & kMask000000FF;
  // ........................fedcba98........................76543210
#  endif
  qw = qw | (qw << 12);
  // ............fedcba98....fedcba98............76543210....76543210

  qw = qw | (qw << 6);
  // ......fedcbaXXdcbaXXdcbaXXdcba98......765432XX5432XX5432XX543210

  return (qw & kMask0303);
  // ......fe......dc......ba......98......76......54......32......10
}

HEADER_INLINE Quarterword Pack0303(uintptr_t ww) {
  // ......fe......dc......ba......98......76......54......32......10

  ww = ww | (ww >> 6);
  // ......fe....fedc....dcba....ba98....9876....7654....5432....3210

  ww = ww | (ww >> 12);
  // ......fe....fedc..fedcbafedcba98dcba9876ba9876549876543276543210

#  ifdef __LP64__
  ww = ww & kMask000000FF;
  // ........................fedcba98........................76543210

  return S_CAST(Quarterword, ww | (ww >> 24));
#  else
  return S_CAST(Quarterword, ww);
#  endif
}

HEADER_INLINE uintptr_t Pack0303Mask(uintptr_t ww) {
  return Pack0303(ww & kMask0303);
}
#endif  // !USE_AVX2

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

HEADER_CINLINE uint64_t RoundUpPow2U64(uint64_t val, uint64_t alignment) {
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

// No need for ModNzU64 in practice, since high bits don't affect result when
// modulus is a power of 2.

// Equivalent to (static_cast<int32_t>(uii) < 0).  Most frequently used on
// possibly-error chromosome indexes.
HEADER_INLINE uint32_t IsI32Neg(uint32_t uii) {
  return uii >> 31;
}

HEADER_INLINE uint32_t abs_i32(int32_t ii) {
  // Arithmetic right shift.  0xffffffffU when ii is negative, 0 otherwise.
  const uint32_t neg_sign_bit = S_CAST(uint32_t, ii >> 31);

  return (S_CAST(uint32_t, ii) ^ neg_sign_bit) - neg_sign_bit;
}

// Returns first number in {2^0, 2^1, ..., 2^{kBitsPerWord - 1}) >= posint.
// Assumes posint in [1, 2^{kBitsPerWord - 1}].
HEADER_INLINE uintptr_t CeilPow2(uintptr_t posint) {
  return k1LU << bsrw(2 * posint - 1);
}

extern uint64_t g_failed_alloc_attempt_size;
// with NDEBUG undefined, may want to define a bunch of macros so that line
// number is printed as well; see e.g.
//   https://stackoverflow.com/questions/15884793/how-to-get-the-name-or-file-and-line-of-caller-method

// Unfortunately, defining the second parameter to be of type void** doesn't do
// the right thing.
// update (18 Feb 2022): looks like inlined pgl_malloc is not compiled as
// intended by gcc 11, due to new ipa-modref pass?  Open to suggestions on how
// to fix this; maybe it's now necessary to define type-specific malloc
// wrappers for clean compilation, ugh...
// update (22 Nov 2025): now explicitly disabling ipa-modref instead of making
// this non-inline, since the latter doesn't hold up under LTO.
#if defined(__cplusplus)
template <class T> HEADER_INLINE BoolErr pgl_malloc(uintptr_t size, T** pp) {
  *pp = S_CAST(T*, malloc(size));
  if (likely(*pp)) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}

#  if (__GNUC__ >= 11) && !defined(__clang__)
#    pragma GCC push_options
#    pragma GCC optimize("-fno-ipa-modref")
#  endif
HEADER_INLINE BoolErr pgl_malloc(uintptr_t size, uintptr_t* addr_ptr) {
  *addr_ptr = R_CAST(uintptr_t, malloc(size));
  if (likely(*addr_ptr)) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}
#  if (__GNUC__ >= 11) && !defined(__clang__)
#    pragma GCC pop_options
#  endif
#else
#  if (__GNUC__ >= 11) && !defined(__clang__)
#    pragma GCC push_options
#    pragma GCC optimize("-fno-ipa-modref")
#  endif
HEADER_INLINE BoolErr pgl_malloc(uintptr_t size, void* pp) {
  *S_CAST(unsigned char**, pp) = S_CAST(unsigned char*, malloc(size));
  if (likely(*S_CAST(unsigned char**, pp))) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}
#  if (__GNUC__ >= 11) && !defined(__clang__)
#    pragma GCC pop_options
#  endif
#endif

// This must be used for all fwrite() calls where len could be >= 2^31, since
// OS X raw fwrite() doesn't work in that case.
static_assert(sizeof(size_t) == sizeof(intptr_t), "plink2_base assumes size_t and intptr_t are synonymous.");
BoolErr fwrite_checked(const void* buf, uintptr_t len, FILE* outfile);

HEADER_INLINE IntErr putc_checked(int32_t ii, FILE* outfile) {
  putc_unlocked(ii, outfile);
  return ferror_unlocked(outfile);
}

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
//   ScanmovPosintCapped(), ScanmovUintCapped(), etc. in plink2_string if
//   you want strtol-like semantics, where the pointer is moved.)
// * Errors out on overflow.  This may be the biggest advantage over atoi().
BoolErr ScanPosintCapped(const char* str_iter, uint64_t cap, uint32_t* valp);

// [0, cap]
BoolErr ScanUintCapped(const char* str_iter, uint64_t cap, uint32_t* valp);

// [-bound, bound]
BoolErr ScanIntAbsBounded(const char* str_iter, uint64_t bound, int32_t* valp);
#else  // not __LP64__
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
// give the shortest-name forms char* return types.
HEADER_INLINE char* memseta(void* target, unsigned char val, uintptr_t ct) {
  memset(target, val, ct);
  return &(S_CAST(char*, target)[ct]);
}

HEADER_INLINE unsigned char* memsetua(void* target, unsigned char val, uintptr_t ct) {
  memset(target, val, ct);
  return &(S_CAST(unsigned char*, target)[ct]);
}

HEADER_CINLINE uintptr_t BitCtToWordCt(uintptr_t val) {
  return DivUp(val, kBitsPerWord);
}

HEADER_CINLINE uintptr_t BitCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kBitsPerCacheline);
}

HEADER_CINLINE uintptr_t Int32CtToCachelineCt(uintptr_t val) {
  return DivUp(val, kInt32PerCacheline);
}

HEADER_CINLINE uintptr_t WordCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kWordsPerCacheline);
}

HEADER_CINLINE uint64_t WordCtToCachelineCtU64(uint64_t val) {
  return DivUpU64(val, kWordsPerCacheline);
}

HEADER_CINLINE uintptr_t Int64CtToCachelineCt(uintptr_t val) {
  return DivUp(val, kInt64PerCacheline);
}

;;;
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

#ifdef USE_SSE42
HEADER_CINLINE uint32_t NypsumWord(uintptr_t val) {
  return __builtin_popcountll(val) + __builtin_popcountll(val & kMaskAAAA);
}
#else
HEADER_CINLINE2 uint32_t NypsumWord(uintptr_t val) {
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
  // Sadly, this was still faster than the LLVM implementation of the intrinsic
  // as of 2016.
  return NypsumWord(val - ((val >> 1) & kMask5555));
}
#endif

#ifdef USE_SSE42
HEADER_INLINE uint32_t Popcount2Words(uintptr_t val0, uintptr_t val1) {
  return __builtin_popcountll(val0) + __builtin_popcountll(val1);
}
#else
HEADER_INLINE uint32_t Popcount2Words(uintptr_t val0, uintptr_t val1) {
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  const uintptr_t four_bit = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  // up to 16 values in 0..12; sum fits in 8 bits
  return (((four_bit & kMask0F0F) + ((four_bit >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

#ifndef __LP64__
HEADER_INLINE uint32_t Popcount4Words(uintptr_t val0, uintptr_t val1, uintptr_t val2, uintptr_t val3) {
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  val2 -= (val2 >> 1) & kMask5555;
  val3 -= (val3 >> 1) & kMask5555;
  const uintptr_t four_bit_0 = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  const uintptr_t four_bit_1 = (val2 & kMask3333) + ((val2 >> 2) & kMask3333) + (val3 & kMask3333) + ((val3 >> 2) & kMask3333);
  return (((four_bit_0 & kMask0F0F) + ((four_bit_0 >> 4) & kMask0F0F) + (four_bit_1 & kMask0F0F) + ((four_bit_1 >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
}
#endif

#ifdef USE_SSE42
HEADER_INLINE uint32_t PopcountHW(uint32_t val) {
  return __builtin_popcount(val);
}

HEADER_INLINE uint32_t PopcountByte(uint32_t val) {
  return __builtin_popcount(val);
}
#else
#  ifdef __LP64__
HEADER_INLINE uint32_t PopcountHW(uint32_t val) {
  val = val - ((val >> 1) & 0x55555555);
  val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
  val = (val + (val >> 4)) & 0x0f0f0f0f;
  return (val * 0x1010101) >> 24;
}
#  else
HEADER_INLINE uint32_t PopcountHW(uint32_t val) {
  val = val - ((val >> 1) & 0x5555);
  val = (val & 0x3333) + ((val >> 2) & 0x3333);
  val = (val + (val >> 4)) & 0x0f0f;
  return (val + (val >> 8)) & 0xff;
}
#  endif

HEADER_INLINE uint32_t PopcountByte(uint32_t val) {
  val = val - ((val >> 1) & 0x55);
  val = (val & 0x33) + ((val >> 2) & 0x33);
  return (val + (val >> 4)) & 0xf;
}
#endif

// Downcasts don't risk alignment issues.
HEADER_INLINE unsigned char* DowncastToUc(void* pp) {
  return S_CAST(unsigned char*, pp);
}

HEADER_INLINE char* DowncastToC(void* pp) {
  return S_CAST(char*, pp);
}

HEADER_INLINE const unsigned char* DowncastKToUc(const void* pp) {
  return S_CAST(const unsigned char*, pp);
}

HEADER_INLINE const char* DowncastKToC(const void* pp) {
  return S_CAST(const char*, pp);
}

HEADER_INLINE Halfword* DowncastWToHW(uintptr_t* pp) {
  return R_CAST(Halfword*, pp);
}

HEADER_INLINE uint32_t* DowncastWToU32(uintptr_t* pp) {
  return R_CAST(uint32_t*, pp);
}

HEADER_INLINE uint16_t* DowncastWToU16(uintptr_t* pp) {
  return R_CAST(uint16_t*, pp);
}

HEADER_INLINE int16_t* DowncastWToI16(uintptr_t* pp) {
  return R_CAST(int16_t*, pp);
}

HEADER_INLINE uint16_t* DowncastU32ToU16(uint32_t* pp) {
  return R_CAST(uint16_t*, pp);
}

HEADER_INLINE const Halfword* DowncastKWToHW(const uintptr_t* pp) {
  return R_CAST(const Halfword*, pp);
}

HEADER_INLINE const uint32_t* DowncastKWToU32(const uintptr_t* pp) {
  return R_CAST(const uint32_t*, pp);
}

HEADER_INLINE const uint16_t* DowncastKWToU16(const uintptr_t* pp) {
  return R_CAST(const uint16_t*, pp);
}

HEADER_INLINE const uint16_t* DowncastKU64ToU16(const uint64_t* pp) {
  return R_CAST(const uint16_t*, pp);
}


HEADER_INLINE uint32_t AlignKToW(const void* prestart, const uintptr_t** result_ptr) {
  const unsigned char* prestart_uc = S_CAST(const unsigned char*, prestart);
  const uint32_t lead_byte_ct = (-R_CAST(uintptr_t, prestart_uc)) % kBytesPerWord;
  *result_ptr = R_CAST(const uintptr_t*, &(prestart_uc[lead_byte_ct]));
  return lead_byte_ct;
}


// Turns out memcpy(&cur_word, bytearr, ct) can't be trusted to be fast when ct
// isn't known at compile time.
//
// ct must be less than sizeof(intptr_t).  ct == 0 handled correctly, albeit
// inefficiently.
#ifndef NO_UNALIGNED
HEADER_INLINE uintptr_t ProperSubwordLoad(const void* bytearr, uint32_t ct) {
  const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
#  ifdef __LP64__
  if (ct >= 4) {
    const uint32_t remainder = ct - 4;
    bytearr_uc = &(bytearr_uc[remainder]);
    uintptr_t cur_word = *R_CAST(const uint32_t*, bytearr_uc);
    if (remainder) {
      cur_word <<= remainder * CHAR_BIT;
      cur_word |= *S_CAST(const uint32_t*, bytearr);
    }
    return cur_word;
  }
#  endif
  if (ct >= 2) {
    const uint32_t remainder = ct & 1;
    uintptr_t cur_word = *R_CAST(const uint16_t*, &(bytearr_uc[remainder]));
    if (remainder) {
      cur_word <<= 8;
      cur_word |= bytearr_uc[0];
    }
    return cur_word;
  }
  return ct? bytearr_uc[0] : 0;
}

HEADER_INLINE uintptr_t SubwordLoad(const void* bytearr, uint32_t ct) {
  if (ct == S_CAST(uint32_t, kBytesPerWord)) {
    return *S_CAST(const uintptr_t*, bytearr);
  }
  return ProperSubwordLoad(bytearr, ct);
}

// ct must be in 1..4.
HEADER_INLINE uint32_t SubU32Load(const void* bytearr, uint32_t ct) {
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

// ct must be 1 or 2.
HEADER_INLINE uint16_t SubU16Load(const void* bytearr, uint32_t ct) {
  if (ct == 1) {
    const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
    return bytearr_uc[0];
  }
  return *S_CAST(const uint16_t*, bytearr);
}

// ct must be in 1..sizeof(Halfword).
#ifdef __LP64__
HEADER_INLINE Halfword SubHWLoad(const void* bytearr, uint32_t ct) {
  return SubU32Load(bytearr, ct);
}
#else
HEADER_INLINE Halfword SubHWLoad(const void* bytearr, uint32_t ct) {
  return SubU16Load(bytearr, ct);
}
#endif

// tried making this non-inline, loop took more than 50% longer
HEADER_INLINE void ProperSubwordStore(uintptr_t cur_word, uint32_t byte_ct, void* target) {
  unsigned char* target_iter = S_CAST(unsigned char*, target);
#  ifdef __LP64__
  if (byte_ct >= 4) {
    *R_CAST(uint32_t*, target_iter) = cur_word;
    if (byte_ct == 4) {
      return;
    }
    const uint32_t remainder = byte_ct - 4;
    target_iter = &(target_iter[remainder]);
    cur_word >>= remainder * CHAR_BIT;
    *R_CAST(uint32_t*, target_iter) = cur_word;
    return;
  }
#  endif
  if (byte_ct & 1) {
    *target_iter++ = cur_word;
    cur_word >>= 8;
  }
  if (byte_ct & 2) {
    *R_CAST(uint16_t*, target_iter) = cur_word;
  }
}

HEADER_INLINE void SubwordStore(uintptr_t cur_word, uint32_t byte_ct, void* target) {
  if (byte_ct == kBytesPerWord) {
    *S_CAST(uintptr_t*, target) = cur_word;
    return;
  }
  ProperSubwordStore(cur_word, byte_ct, target);
}

// byte_ct must be in 1..4.
HEADER_INLINE void SubU32Store(uint32_t cur_uint, uint32_t byte_ct, void* target) {
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
#else  // NO_UNALIGNED
HEADER_INLINE uintptr_t ProperSubwordLoad(const void* bytearr, uint32_t ct) {
  uintptr_t cur_word = 0;
  memcpy(&cur_word, bytearr, ct);
  return cur_word;
}

HEADER_INLINE uintptr_t SubwordLoad(const void* bytearr, uint32_t ct) {
  uintptr_t cur_word = 0;
  memcpy(&cur_word, bytearr, ct);
  return cur_word;
}

HEADER_INLINE uint32_t SubU32Load(const void* bytearr, uint32_t ct) {
  uint32_t cur_uint = 0;
  memcpy(&cur_uint, bytearr, ct);
  return cur_uint;
}

HEADER_INLINE Halfword SubHWLoad(const void* bytearr, uint32_t ct) {
  Halfword cur_uint = 0;
  memcpy(&cur_uint, bytearr, ct);
  return cur_uint;
}

HEADER_INLINE void ProperSubwordStore(uintptr_t cur_word, uint32_t byte_ct, void* target) {
  memcpy(target, &cur_word, byte_ct);
}

HEADER_INLINE void SubwordStore(uintptr_t cur_word, uint32_t byte_ct, void* target) {
  memcpy(target, &cur_word, byte_ct);
}

HEADER_INLINE void SubU32Store(uint32_t cur_uint, uint32_t byte_ct, void* target) {
  memcpy(target, &cur_uint, byte_ct);
}
#endif // NO_UNALIGNED

HEADER_INLINE uint64_t SubU64Load(const void* bytearr, uint32_t ct) {
#ifdef __LP64__
  return SubwordLoad(bytearr, ct);
#else
  uint64_t cur_u64 = 0;
  memcpy(&cur_u64, bytearr, ct);
  return cur_u64;
#endif
}

HEADER_INLINE void ProperSubwordStoreMov(uintptr_t cur_word, uint32_t byte_ct, unsigned char** targetp) {
  ProperSubwordStore(cur_word, byte_ct, *targetp);
  *targetp += byte_ct;
}

HEADER_INLINE void SubwordStoreMov(uintptr_t cur_word, uint32_t byte_ct, unsigned char** targetp) {
  SubwordStore(cur_word, byte_ct, *targetp);
  *targetp += byte_ct;
}

HEADER_INLINE void SubU32StoreMov(uint32_t cur_uint, uint32_t byte_ct, unsigned char** targetp) {
  SubU32Store(cur_uint, byte_ct, *targetp);
  *targetp += byte_ct;
}

#ifndef NO_UNALIGNED
#  ifdef __LP64__
HEADER_INLINE void SubU64StoreMov(uint64_t cur_u64, uint32_t byte_ct, unsigned char** targetp) {
  return SubwordStoreMov(cur_u64, byte_ct, targetp);
}
#  else
HEADER_INLINE void SubU64StoreMov(uint64_t cur_u64, uint32_t byte_ct, unsigned char** targetp) {
  if (byte_ct > 4) {
    *R_CAST(uint32_t*, *targetp) = cur_u64;
    *targetp += 4;
    byte_ct -= 4;
    cur_u64 >>= 32;
  }
  return SubU32StoreMov(cur_u64, byte_ct, targetp);
}
#  endif
#else
HEADER_INLINE void SubU64StoreMov(uint64_t cur_u64, uint32_t byte_ct, unsigned char** targetp) {
  memcpy(*targetp, &cur_u64, byte_ct);
  *targetp += byte_ct;
}
#endif


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

// Common pattern is to initialize a data structure, then use it without
// mutating it, then free it.
// Yes, in principle you should have a non-const pointer to the data structure
// to free it.  But my experience is that that makes too many beneficial uses
// of const less ergonomic, without offering enough in return.  So we cheat a
// bit, in a way that can be systematically cleaned up later if necessary.
HEADER_INLINE void free_const(const void* memptr) {
  free(K_CAST(void*, memptr));
}

HEADER_INLINE void free_cond(const void* memptr) {
  if (memptr) {
    free_const(memptr);
  }
}


HEADER_INLINE char* memcpya(void* __restrict target, const void* __restrict source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(S_CAST(char*, target)[ct]);
}

HEADER_INLINE unsigned char* memcpyua(void* __restrict target, const void* __restrict source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(S_CAST(unsigned char*, target)[ct]);
}

// Tried beating memcpy for usually-small strings not known to have length <=
// 8, gave up.

HEADER_INLINE char* strcpya(char* __restrict dst, const void* __restrict src) {
  const uintptr_t slen = strlen(S_CAST(const char*, src));
  return memcpya(dst, src, slen);
}

#if __cplusplus >= 201103L
constexpr uint32_t CompileTimeSlen(const char* k_str) {
  return k_str[0]? (1 + CompileTimeSlen(&(k_str[1]))) : 0;
}
#endif

// A few more string-rendering functions that would normally live in
// plink2_string, to work around PRI{d,u}PTR and PRI{d,u}64 warning on CRAN
// Windows builds.

extern const uint16_t kDigitPair[];

#ifndef NO_UNALIGNED
// Unfortunately, memcpy[a]_k() implementation has components which force it to
// be in plink2_simd, and at least when I checked a few years ago, compilers
// were still prone to generating bad code for known-size-2 (though not 1, 4,
// 8) memcpys.
// (obvious todo: recheck this for current set of production compilers, then
// recheck when updating build machine(s), until we're able to delete all this
// code.)
HEADER_INLINE void memcpy_k2(void* __restrict dst, const void* __restrict src) {
  *S_CAST(uint16_t*, dst) = *S_CAST(const uint16_t*, src);
}

#else
HEADER_INLINE void memcpy_k2(void* __restrict dst, const void* __restrict src) {
  memcpy(dst, src, 2);
}
#endif

HEADER_INLINE void memcpy_kw(void* __restrict dst, const void* __restrict src) {
  memcpy(dst, src, kBytesPerWord);
}

#ifdef __LP64__
HEADER_INLINE void memcpy_khw(void* __restrict dst, const void* __restrict src) {
  memcpy(dst, src, 4);
}

HEADER_INLINE void memcpy_kqw(void* __restrict dst, const void* __restrict src) {
  memcpy_k2(dst, src);
}
#else
HEADER_INLINE void memcpy_khw(void* __restrict dst, const void* __restrict src) {
  memcpy_k2(dst, src);
}

HEADER_INLINE void memcpy_kqw(void* __restrict dst, const void* __restrict src) {
  memcpy(dst, src, 1);
}
#endif

HEADER_INLINE char* memcpya_k2(void* __restrict dst, const void* __restrict src) {
  memcpy_k2(dst, src);
  char* dst_c = S_CAST(char*, dst);
  return &(dst_c[2]);
}

/*
HEADER_INLINE char* memcpya_kw(void* __restrict dst, const void* __restrict src) {
  memcpy_kw(dst, src);
  char* dst_c = S_CAST(char*, dst);
  return &(dst_c[kBytesPerWord]);
}
*/


char* u32toa(uint32_t uii, char* start);

HEADER_INLINE char* uitoa_z4(uint32_t uii, char* start) {
  uint32_t quotient = uii / 100;
  assert(quotient < 100);
  uii -= 100 * quotient;
  start = memcpya_k2(start, &(kDigitPair[quotient]));
  return memcpya_k2(start, &(kDigitPair[uii]));
}

HEADER_INLINE char* u32toa_z6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  start = memcpya_k2(start, &(kDigitPair[quotient]));
  return uitoa_z4(uii - 10000 * quotient, start);
}

HEADER_INLINE char* uitoa_z8(uint32_t uii, char* start) {
  uint32_t quotient = uii / 1000000;
  start = memcpya_k2(start, &(kDigitPair[quotient]));
  return u32toa_z6(uii - 1000000 * quotient, start);
}

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

HEADER_INLINE void CopyFromUnalignedW(uintptr_t* dst, const unsigned char* src) {
  memcpy_kw(dst, src);
}

HEADER_INLINE void CopyFromUnalignedU32(uint32_t* dst, const unsigned char* src) {
  memcpy(dst, src, sizeof(int32_t));
}

HEADER_INLINE void CopyFromUnalignedU16(uint16_t* dst, const unsigned char* src) {
  memcpy_k2(dst, src);
}

HEADER_INLINE void CopyFromUnalignedI16(int16_t* dst, const unsigned char* src) {
  memcpy_k2(dst, src);
}

HEADER_INLINE void CopyFromUnalignedF(float* dst, const unsigned char* src) {
  memcpy(dst, src, sizeof(float));
}

// [u]int16_t (and unsigned char) arithmetic is slower than uint32_t/uint64_t
// arithmetic on some CPUs.
// todo: check whether uint32_t vs. uintptr_t return value matters
HEADER_INLINE uint32_t CopyFromUnalignedU16ZX(const unsigned char* src) {
#ifndef NO_UNALIGNED
  return *R_CAST(const uint16_t*, src);
#else
  uint16_t cur_u16;
  memcpy_k2(&cur_u16, src);
  return cur_u16;
#endif
}

HEADER_INLINE int32_t CopyFromUnalignedI16ZX(const unsigned char* src) {
#ifndef NO_UNALIGNED
  return *R_CAST(const int16_t*, src);
#else
  int16_t cur_i16;
  memcpy_k2(&cur_i16, src);
  return cur_i16;
#endif
}

HEADER_INLINE void CopyToUnalignedW(unsigned char* dst, const uintptr_t* src) {
  memcpy_kw(dst, src);
}

HEADER_INLINE void CopyToUnalignedU64(unsigned char* dst, const uint64_t* src) {
  memcpy(dst, src, sizeof(int64_t));
}

HEADER_INLINE void CopyToUnalignedU32(unsigned char* dst, const uint32_t* src) {
  memcpy(dst, src, sizeof(int32_t));
}

HEADER_INLINE void CopyToUnalignedF(unsigned char* dst, const float* src) {
  memcpy(dst, src, sizeof(float));
}

HEADER_INLINE void CopyFromUnalignedOffsetW(uintptr_t* dst, const unsigned char* src, uintptr_t offset) {
  memcpy_kw(dst, &(src[offset * kBytesPerWord]));
}

HEADER_INLINE void CopyFromUnalignedOffsetHW(Halfword* dst, const unsigned char* src, uintptr_t offset) {
  memcpy_khw(dst, &(src[offset * sizeof(Halfword)]));
}

HEADER_INLINE void CopyFromUnalignedOffsetQW(Quarterword* dst, const unsigned char* src, uintptr_t offset) {
  memcpy_kqw(dst, &(src[offset * sizeof(Quarterword)]));
}

HEADER_INLINE void CopyFromUnalignedOffsetU64(uint64_t* dst, const unsigned char* src, uintptr_t offset) {
  memcpy(dst, &(src[offset * sizeof(int64_t)]), sizeof(int64_t));
}

HEADER_INLINE void CopyFromUnalignedOffsetU32(uint32_t* dst, const unsigned char* src, uintptr_t offset) {
  memcpy(dst, &(src[offset * sizeof(int32_t)]), sizeof(int32_t));
}

HEADER_INLINE void CopyFromUnalignedOffsetI32(int32_t* dst, const unsigned char* src, uintptr_t offset) {
  memcpy(dst, &(src[offset * sizeof(int32_t)]), sizeof(int32_t));
}

HEADER_INLINE void CopyFromUnalignedOffsetU16(uint16_t* dst, const unsigned char* src, uintptr_t offset) {
  memcpy_k2(dst, &(src[offset * sizeof(int16_t)]));
}

HEADER_INLINE void CopyFromUnalignedOffsetI16(int16_t* dst, const unsigned char* src, uintptr_t offset) {
  memcpy_k2(dst, &(src[offset * sizeof(int16_t)]));
}

HEADER_INLINE void CopyFromUnalignedOffsetF(float* dst, const unsigned char* src, uintptr_t offset) {
  memcpy(dst, &(src[offset * sizeof(float)]), sizeof(float));
}

HEADER_INLINE uint32_t CopyFromUnalignedOffsetU16ZX(const unsigned char* src, uintptr_t offset) {
#ifndef NO_UNALIGNED
  return R_CAST(const uint16_t*, src)[offset];
#else
  uint16_t cur_u16;
  memcpy_k2(&cur_u16, &(src[offset * sizeof(int16_t)]));
  return cur_u16;
#endif
}

HEADER_INLINE int32_t CopyFromUnalignedOffsetI16ZX(const unsigned char* src, uintptr_t offset) {
#ifndef NO_UNALIGNED
  return R_CAST(const int16_t*, src)[offset];
#else
  int16_t cur_i16;
  memcpy_k2(&cur_i16, &(src[offset * sizeof(int16_t)]));
  return cur_i16;
#endif
}

HEADER_INLINE void CopyToUnalignedOffsetW(unsigned char* dst, const uintptr_t* src, uintptr_t offset) {
  memcpy(&(dst[offset * sizeof(intptr_t)]), src, sizeof(intptr_t));
}

HEADER_INLINE void CopyToUnalignedOffsetHW(unsigned char* dst, const Halfword* src, uintptr_t offset) {
  memcpy_khw(&(dst[offset * sizeof(Halfword)]), src);
}

HEADER_INLINE void CopyToUnalignedOffsetU64(unsigned char* dst, const uint64_t* src, uintptr_t offset) {
  memcpy(&(dst[offset * sizeof(int64_t)]), src, sizeof(int64_t));
}

HEADER_INLINE void CopyToUnalignedOffsetU32(unsigned char* dst, const uint32_t* src, uintptr_t offset) {
  memcpy(&(dst[offset * sizeof(int32_t)]), src, sizeof(int32_t));
}

HEADER_INLINE void CopyToUnalignedOffsetU16(unsigned char* dst, const uint16_t* src, uintptr_t offset) {
  memcpy_k2(&(dst[offset * sizeof(int16_t)]), src);
}

HEADER_INLINE void CopyToUnalignedOffsetI16(unsigned char* dst, const int16_t* src, uintptr_t offset) {
  memcpy_k2(&(dst[offset * sizeof(int16_t)]), src);
}

HEADER_INLINE void CopyToUnalignedOffsetF(unsigned char* dst, const float* src, uintptr_t offset) {
  memcpy(&(dst[offset * sizeof(float)]), src, sizeof(float));
}

HEADER_INLINE void CopyFromUnalignedIncrW(uintptr_t* dst, const unsigned char** srcp) {
  memcpy(dst, *srcp, kBytesPerWord);
  *srcp += kBytesPerWord;
}

HEADER_INLINE void CopyFromUnalignedIncrU32(uint32_t* dst, const unsigned char** srcp) {
  memcpy(dst, *srcp, sizeof(int32_t));
  *srcp += sizeof(int32_t);
}

HEADER_INLINE void CopyFromUnalignedIncrU16(uint16_t* dst, const unsigned char** srcp) {
  memcpy_k2(dst, *srcp);
  *srcp += sizeof(int16_t);
}

HEADER_INLINE void CopyFromUnalignedIncrI16(int16_t* dst, const unsigned char** srcp) {
  memcpy_k2(dst, *srcp);
  *srcp += sizeof(int16_t);
}

HEADER_INLINE uint32_t CopyFromUnalignedIncrU16ZX(const unsigned char** srcp) {
#ifndef NO_UNALIGNED
  const uint32_t result = *R_CAST(const uint16_t*, *srcp);
#else
  uint16_t result;
  memcpy_k2(&result, *srcp);
#endif
  *srcp += sizeof(int16_t);
  return result;
}

HEADER_INLINE int32_t CopyFromUnalignedIncrI16ZX(const unsigned char** srcp) {
#ifndef NO_UNALIGNED
  const int32_t result = *R_CAST(const int16_t*, *srcp);
#else
  int16_t result;
  memcpy_k2(&result, *srcp);
#endif
  *srcp += sizeof(int16_t);
  return result;
}

// no need to spell out 'CopyToUnalignedIncr' here; and we append constants
// often enough to justify no src pointer

HEADER_INLINE void AppendW(uintptr_t ulii, unsigned char** targetp) {
  memcpy(*targetp, &ulii, kBytesPerWord);
  *targetp += kBytesPerWord;
}

HEADER_INLINE void AppendU64(uint64_t ullii, unsigned char** targetp) {
  memcpy(*targetp, &ullii, sizeof(int64_t));
  *targetp += sizeof(int64_t);
}

HEADER_INLINE void AppendU32(uint32_t uii, unsigned char** targetp) {
  memcpy(*targetp, &uii, sizeof(int32_t));
  *targetp += sizeof(int32_t);
}

HEADER_INLINE void AppendU16(uint16_t usii, unsigned char** targetp) {
  memcpy_k2(*targetp, &usii);
  *targetp += sizeof(int16_t);
}

HEADER_INLINE void CAppendW(uintptr_t ulii, char** targetp) {
  memcpy(*targetp, &ulii, kBytesPerWord);
  *targetp += kBytesPerWord;
}

HEADER_INLINE void CAppendU32(uint32_t uii, char** targetp) {
  memcpy(*targetp, &uii, sizeof(int32_t));
  *targetp += sizeof(int32_t);
}

HEADER_INLINE void CAppendU16(uint16_t usii, char** targetp) {
  memcpy_k2(*targetp, &usii);
  *targetp += sizeof(int16_t);
}


// now compiling with gcc >= 4.4 (or clang equivalent) on all platforms, so
// safe to use memset everywhere
HEADER_INLINE void ZeroU32Arr(uintptr_t entry_ct, uint32_t* u32arr) {
  memset(u32arr, 0, entry_ct * sizeof(int32_t));
}

HEADER_INLINE void ZeroWArr(uintptr_t entry_ct, uintptr_t* warr) {
  memset(warr, 0, entry_ct * sizeof(intptr_t));
}

HEADER_INLINE void ZeroU64Arr(uintptr_t entry_ct, uint64_t* u64arr) {
  memset(u64arr, 0, entry_ct * sizeof(int64_t));
}

HEADER_INLINE void ZeroPtrArr(uintptr_t entry_ct, void* pp) {
  memset(pp, 0, entry_ct * sizeof(intptr_t));
}

HEADER_INLINE void ZeroHwArr(uintptr_t entry_ct, Halfword* hwarr) {
  memset(hwarr, 0, entry_ct * sizeof(Halfword));
}

HEADER_INLINE void ZeroFArr(uintptr_t entry_ct, float* farr) {
  for (uintptr_t ulii = 0; ulii != entry_ct; ulii++) {
    *farr++ = 0.0;
  }
}

HEADER_INLINE void ZeroDArr(uintptr_t entry_ct, double* darr) {
  for (uintptr_t ulii = 0; ulii != entry_ct; ulii++) {
    *darr++ = 0.0;
  }
}

HEADER_INLINE void SetAllWArr(uintptr_t entry_ct, uintptr_t* warr) {
  // todo: test this against vecset()
  for (uintptr_t idx = 0; idx != entry_ct; ++idx) {
    warr[idx] = ~k0LU;
  }
}

HEADER_INLINE void SetAllU32Arr(uintptr_t entry_ct, uint32_t* u32arr) {
  for (uintptr_t ulii = 0; ulii != entry_ct; ulii++) {
    *u32arr++ = ~0U;
  }
}

uint32_t MaxElementU32(const uint32_t* u32arr, uintptr_t entry_ct);

double MaxElementD(const double* darr, uintptr_t entry_ct);

double MinElementD(const double* darr, uintptr_t entry_ct);


// tried _bzhi_u64() in AVX2 case, it was actually worse on my Mac (more
// opaque to compiler?)
// todo: check gcc behavior since it may be different: see
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=82298 .
//
// This is undefined if idx == kBitsPerWord.
HEADER_INLINE uintptr_t bzhi(uintptr_t ww, uint32_t idx) {
  return ww & ((k1LU << idx) - k1LU);
}

// This is undefined if idx == 0.
HEADER_INLINE uintptr_t bzhi_max(uintptr_t ww, uint32_t idx) {
  return ww & ((~k0LU) >> (kBitsPerWord - idx));
}

// Don't bother defining blsr(), compiler should automatically use the
// instruction under -mbmi and regular code is more readable?  (again, should
// verify this is true for gcc)

HEADER_INLINE uint32_t BytesToRepresentNzU32(uint32_t uii) {
  return 1 + (bsru32(uii) / CHAR_BIT);
}

// todo: make sure these are efficient for small ct
HEADER_INLINE void u16set(void* dst, uint16_t usii, uintptr_t ct) {
  uint16_t* dst_u16 = S_CAST(uint16_t*, dst);
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst_u16[ulii] = usii;
  }
}

HEADER_INLINE char* u16setsa(char* dst, uint16_t usii, uintptr_t ct) {
  u16set(dst, usii, ct);
  return &(dst[ct * 2]);
}

HEADER_INLINE uintptr_t ClearBottomSetBits(uint32_t ct, uintptr_t ulii) {
#ifdef USE_AVX2
  return _pdep_u64((~k0LU) << ct, ulii);
#else
  for (uint32_t uii = 0; uii != ct; ++uii) {
    ulii &= ulii - 1;
  }
  return ulii;
#endif
}

HEADER_INLINE uint32_t WordBitIdxToUidx(uintptr_t ulii, uint32_t bit_idx) {
  return ctzw(ClearBottomSetBits(bit_idx, ulii));
}

CONSTI32(kNybblesPerWord, 2 * kBytesPerWord);
CONSTI32(kNybblesPerCacheline, 2 * kCacheline);

HEADER_CINLINE uintptr_t NybbleCtToByteCt(uintptr_t val) {
  return DivUp(val, 2);
}

HEADER_CINLINE uintptr_t NybbleCtToWordCt(uintptr_t val) {
  return DivUp(val, kNybblesPerWord);
}

HEADER_INLINE uintptr_t GetNybbleArrEntry(const uintptr_t* nybblearr, uint32_t idx) {
  return (nybblearr[idx / kBitsPerWordD4] >> (4 * (idx % kBitsPerWordD4))) & 15;
}

// Returns zero when ww has no zero bytes, and a word where the lowest set bit
// is at position 8x + 7 when the first zero byte is [8x .. 8x+7].
HEADER_INLINE uintptr_t DetectFirstZeroByte(uintptr_t ww) {
  return (ww - kMask0101) & (~ww) & (kMask0101 * 0x80);
}

// From TAOCP 4a, 7.1.3, (91).
// Position 8x + 7 is always set iff byte x is zero.  All other bits are always
// zero.
HEADER_INLINE uintptr_t DetectAllZeroBytes(uintptr_t ww) {
  return (kMask0101 * 0x80) & (~(ww | ((ww | (kMask0101 * 0x80)) - kMask0101)));
}

HEADER_INLINE uintptr_t DetectAllZeroNybbles(uintptr_t ww) {
  return (kMask1111 * 8) & (~(ww | ((ww | (kMask1111 * 8)) - kMask1111)));
}


HEADER_INLINE void* arena_alloc_raw(uintptr_t size, unsigned char** arena_bottom_ptr) {
  assert(!(size % kCacheline));
  unsigned char* alloc_ptr = *arena_bottom_ptr;
  *arena_bottom_ptr = &(alloc_ptr[size]);
  return alloc_ptr;
}

HEADER_INLINE void* arena_alloc_raw_rd(uintptr_t size, unsigned char** arena_bottom_ptr) {
  unsigned char* alloc_ptr = *arena_bottom_ptr;
  *arena_bottom_ptr = &(alloc_ptr[RoundUpPow2(size, kCacheline)]);
  return alloc_ptr;
}

// A VINT is a sequence of bytes where each byte stores just 7 bits of an
// an integer, and the high bit is set when the integer has more nonzero bits.
// See e.g.
//   https://developers.google.com/protocol-buffers/docs/encoding#varints
// (Note that protocol buffers used "group varints" at one point, but then
// abandoned them.  I suspect they'd be simultaneously slower and less
// compact here.)

HEADER_INLINE unsigned char* Vint32Append(uint32_t uii, unsigned char* buf) {
  while (uii > 127) {
    *buf++ = (uii & 127) + 128;
    uii >>= 7;
  }
  *buf++ = uii;
  return buf;
}

// Returns 0x80000000U on read-past-end instead of UINT32_MAX so overflow check
// works properly in 32-bit build.  Named "GetVint31" to make it more obvious
// that a 2^31 return value can't be legitimate.
HEADER_INLINE uint32_t GetVint31(const unsigned char* buf_end, const unsigned char** buf_iterp) {
  if (likely(buf_end > (*buf_iterp))) {
    uint32_t vint32 = *((*buf_iterp)++);
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (likely(buf_end > (*buf_iterp))) {
      uint32_t uii = *((*buf_iterp)++);
      vint32 |= (uii & 127) << shift;
      if (uii <= 127) {
        return vint32;
      }
      shift += 7;
      // currently don't check for shift >= 32 (that's what ValidateVint31()
      // is for).
    }
  }
  return 0x80000000U;
}

// Number of bytes required to encode ulii as a varint.
HEADER_INLINE uint32_t VintBytect(uintptr_t ulii) {
  if (ulii < 128) {
    // bsrw(0) is undefined.
    return 1;
  }
  return 1 + (bsrw(ulii) / 7);
}

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
#ifdef CPP11_TYPE_ENFORCEMENT

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
  tname& operator=(const tname& rhs) = default; \
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
  tname& operator=(const tname& rhs) = default; \
  \
private: \
  uint64_t value_; \
}

#  define ENUM_U31_DEF_START() typedef enum : uint32_t {
#  define ENUM_U31_DEF_END(tname) } tname

#else  // !CPP11_TYPE_ENFORCEMENT

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

// This supports private struct members in code that still compiles as C.
//
// Internal code should access these members with GET_PRIVATE(), and define a
// pair of public C++-only GET_PRIVATE_...() member functions (one const and
// one non-const) that each return a reference to the member; see plink2_thread
// for examples.  In addition, .cc API code should define a small number of
// standard GET_PRIVATE() accessors at the top of the file, and practically all
// private-member access should occur through those file-scope accessors; this
// keeps the surface area under control.
//
// (Tried to define a DECLARE_PRIVATE(typ, member) macro as well, but didn't
// get that to work.  This is already pretty painless if it's restricted to key
// APIs, though.)
//
// probable todo: see if the intended effect can be achieved in a simpler
// manner with well-chosen explicit and implicit type conversions.
#ifdef __cplusplus
#  define GET_PRIVATE(par, member) (par).GET_PRIVATE_ ## member()
#else
#  define GET_PRIVATE(par, member) (par).member
#endif

int32_t u32cmp(const void* aa, const void* bb);

#ifdef __cplusplus
#  define STD_SORT(ct, fallback_cmp, arr) std::sort(&((arr)[0]), (&((arr)[ct])))
#  if __cplusplus >= 201902L && defined(__GNUC__) && !defined(NO_UNALIGNED) && !defined(__clang__)
// Keep in sync with "#include <execution>" condition.

// This should only be used for arrays of length >= variant_ct or sample_ct
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

HEADER_INLINE void swap_u32(uint32_t* uiip, uint32_t* ujjp) {
  const uint32_t swaptmp = *uiip;
  *uiip = *ujjp;
  *ujjp = swaptmp;
}

HEADER_INLINE void swap_i32(int32_t* iip, int32_t* jjp) {
  const int32_t swaptmp = *iip;
  *iip = *jjp;
  *jjp = swaptmp;
}

HEADER_INLINE void swap_f32(float* fxxp, float* fyyp) {
  const float swaptmp = *fxxp;
  *fxxp = *fyyp;
  *fyyp = swaptmp;
}

HEADER_INLINE void swap_f64(double* dxxp, double* dyyp) {
  const double swaptmp = *dxxp;
  *dxxp = *dyyp;
  *dyyp = swaptmp;
}

HEADER_INLINE void swap_kcp(const char** s1p, const char** s2p) {
  const char* swaptmp = *s1p;
  *s1p = *s2p;
  *s2p = swaptmp;
}

HEADER_INLINE void swap_cp(char** s1p, char** s2p) {
  char* swaptmp = *s1p;
  *s1p = *s2p;
  *s2p = swaptmp;
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_BASE_H__
