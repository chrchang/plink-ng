#ifndef __PGENLIB_INTERNAL_H__
#define __PGENLIB_INTERNAL_H__

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


// Low-level C99/C++03/C++11 library for reading .pgen (PLINK 2.0 binary) files
// (designed to produce good lowest-common-denominator binaries across
// Windows/OS X/Linux).  We try to benefit from as much C++ type safety as we
// can without either breaking compatibility with C-only codebases or making
// extension of pgenlib/plink2 code more difficult than the old type-unsafe
// style.
//
// File format design:
// - With the header loaded, it is possible to efficiently access a variant by
//   its index.  Since records can now be variable-length, this sometimes
//   requires storage of record lengths.
// - Due to the power of LD-based compression, we permit a variant record to
//   just store a list of differences from an earlier, fully stored variant.
//   However, only short-range dependence is permitted; sequential processing
//   of the file only requires caching of the most recent explicitly stored
//   variant.
// - Like the plink1 format, this is balanced for relatively easy reading and
//   writing; in particular, the mode-0x10/0x11 header is not read-optimized,
//   it passes up some obvious compression opportunities which would make it
//   more difficult to write e.g. an efficient file merger.  This isn't a big
//   deal if we don't have a huge number of one-sample .pgen files sharing a
//   single .bim file (or equivalent).  (If they don't share the same .bim
//   file, .bim overhead > .pgen overhead.)  If we ever do, we can define an
//   additional mode to handle that case more efficiently.
// - Building blocks are arrays of 1-bit, 2-bit, 4-bit, 1-byte, 2-byte, 3-byte,
//   and 4-byte values.  3/5/6/7(/9...)-bit values don't play well with
//   bitwise operations, and when it's important, there's usually a natural way
//   to split them into power-of-2-bit components.
//   (unsigned integers known to be smaller than 2^24, but not known to be
//   smaller than 2^16, are stored as 3-byte values on disk and "decompressed"
//   to uint32_t during loading.)
// - Missing value is usually all-1s.  (Only exceptions right now: plink1
//   backward compatibility mode; presence/absence of rare alts for variants
//   with >2 alt alleles is an array of 1-bit values, where absence = 0; and
//   presence/absence of phasing info is similar.)  Time to move away from 01
//   nonsense.
// - Individual variant records are prohibited from being >= 4GB, to reduce
//   integer overflow issues.  (This may be reduced to 2GB later, but I'll
//   attempt to handle the 2-4GB range properly for now since it's conceivable
//   for multiallelic records in very large datasets to reach that size.)
// - (later todo: include stuff like file creation command in .bim successor
//   header; that doesn't really belong in a binary file.)

// Parameter conventions:
// - Input parameters, then in/out, then pure outputs, then temporary buffers.
//   Reference-style input parameters tend to go in the very front, to make it
//   more obvious that they aren't in/out.
// - "bitarr" indicates a word-aligned, packed array of bits, while "bitvec"
//   indicates vector-alignment in 64-bit builds.  ("vector" always means SIMD
//   inputs/outputs here; C++ std::vector is not used in this codebase.)
// - "quaterarr" indicates a word-aligned, packed array of 2-bit values, while
//   "quatervec" is the vector-aligned equivalent.  Similarly, "hexadecarr"
//   marks the much rarer case of a packed array of 4-bit values, etc.
// - "quatervec_01" indicates a packed, vector-aligned array of 2-bit values
//   where each value is zero or one.  This data structure was used quite a bit
//   bit by plink 1.9 for operating on a subset of a 2-bit-genotype array.
// - "genovec" indicates a quatervec containing genotype information.
// - "interleaved_vec" is the plink 2.0 replacement for quatervec_01: we
//   basically stack pairs of adjacent vectors on top of each other and unpack
//   on the fly, since that tends to be faster than having to access twice as
//   much memory.
// - Most pointers are stationary; moving pointers have an _iter suffix.


// 10000 * major + 100 * minor + patch
// Exception to CONSTU31, since we want the preprocessor to have access to this
// value.  Named with all caps as a consequence.
#define PGENLIB_INTERNAL_VERNUM 602


#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS 1
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
  #if __cplusplus <= 199711L
    #define static_assert(cond, msg)
  #endif
#else
  #define HEADER_INLINE static inline
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

// other configuration-ish values needed by plink2_common subset
typedef unsigned char alt_allele_ct_t;
// don't use CONSTU31 for this since it may need the 32nd bit in the future
#define kPglMaxAltAlleleCt ((uint32_t)((alt_allele_ct_t)(-2)))

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
  #else
    CONSTU31(kBytesPerVec, 16);
    CONSTU31(kBytesPerFVec, 16);
    typedef uintptr_t vul_t __attribute__ ((vector_size (16)));
    typedef float vf_t __attribute__ ((vector_size (16)));
  #endif
  CONSTU31(kBitsPerWord, 64);
  CONSTU31(kBitsPerWordLog2, 6);

  typedef uint32_t halfword_t;
  typedef uint16_t quarterword_t;

  #ifdef USE_AVX2
    #define VCONST_UL(xx) {xx, xx, xx, xx}
    #define vul_setzero() (vul_t)_mm256_setzero_si256()
    #define vul_rshift(vv, ct) ((vul_t)_mm256_srli_epi64((__m256i)(vv), ct))
    #define vul_lshift(vv, ct) ((vul_t)_mm256_slli_epi64((__m256i)(vv), ct))
  #else
    #define VCONST_UL(xx) {xx, xx}
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

static_assert(CHAR_BIT == 8, "pgenlib_internal requires CHAR_BIT == 8.");
static_assert(sizeof(int32_t) == 4, "pgenlib_internal requires sizeof(int32_t) == 4.");
static_assert(sizeof(int64_t) == 8, "pgenlib_internal requires sizeof(int64_t) == 8.");
static_assert(sizeof(intptr_t) == kBytesPerWord, "pgenlib_internal requires sizeof(intptr_t) == kBytesPerWord.");

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
static_assert(kPglFnamesize >= PATH_MAX, "pgenlib_internal assumes PATH_MAX <= 4096.  (Safe to increase kPglFnamesize to address this, up to 131072.)");
#endif

// safe errstr_buf size for pgen_init_phase{1,2}()
CONSTU31(kPglErrstrBufBlen, kPglFnamesize + 256);


typedef union {
  vul_t vi;
  uintptr_t u8[kBitsPerVec / kBitsPerWord];
  uint32_t u4[kBytesPerVec / sizeof(int32_t)];
} univec_t;

typedef union {
  vf_t vf;
  float f4[kBytesPerFVec / sizeof(float)];
} univecf_t;

// sum must fit in 16 bits
HEADER_INLINE uintptr_t univec_hsum_16bit(univec_t uv) {
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
HEADER_INLINE uintptr_t univec_hsum_32bit(univec_t uv) {
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

HEADER_INLINE uintptr_t unpack_halfword_to_word(uintptr_t hw) {
#ifdef __LP64__
  hw = (hw | (hw << 16)) & kMask0000FFFF;
#endif
  hw = (hw | (hw << 8)) & kMask00FF;
  hw = (hw | (hw << 4)) & kMask0F0F;
  hw = (hw | (hw << 2)) & kMask3333;
  return ((hw | (hw << 1)) & kMask5555);
}

HEADER_INLINE halfword_t pack_word_to_halfword(uintptr_t ww) {
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
HEADER_INLINE uintptr_t round_down_pow2(uintptr_t val, uintptr_t alignment) {
  const uintptr_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return val & (~alignment_m1);
}

HEADER_INLINE uint64_t round_down_pow2_ull(uint64_t val, uint64_t alignment) {
  const uint64_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return val & (~alignment_m1);
}

HEADER_INLINE uintptr_t round_up_pow2(uintptr_t val, uintptr_t alignment) {
  const uintptr_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
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

HEADER_INLINE uint32_t abs_int32(int32_t ii) {
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
static_assert(sizeof(size_t) == sizeof(intptr_t), "pgenlib_internal assumes size_t and intptr_t are synonymous.");
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
//   1000 rather than a nonnumeric token, and "98.6" is treated as 98.  (May
//   want to write another function with strtol-like semantics, returning a
//   pointer to the end of the string so the caller can decide whether to error
//   out on a nonspace terminator.  Don't need the base parameter, though...)
// * Errors out on overflow.
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

HEADER_INLINE uint32_t popcount2_long(uintptr_t val) {
#ifdef USE_SSE42
  return __builtin_popcountll(val) + __builtin_popcountll(val & kMaskAAAA);
#else
  val = (val & kMask3333) + ((val >> 2) & kMask3333);
  return (((val + (val >> 4)) & kMask0F0F) * kMask0101) >> (kBitsPerWord - 8);
#endif
}

HEADER_INLINE uint32_t popcount_long(uintptr_t val) {
  // the simple version, good enough for all non-time-critical stuff
  // (without SSE4.2, popcount_longs() tends to be >3x as fast on arrays.
  // with SSE4.2, there's no noticeable difference.)
#ifdef USE_SSE42
  return __builtin_popcountll(val);
#else
  // sadly, this is still faster than the clang implementation of the intrinsic
  // as of 2016
  return popcount2_long(val - ((val >> 1) & kMask5555));
#endif
}

#ifdef USE_SSE42
HEADER_INLINE uint32_t popcount01_long(uintptr_t val) {
  return popcount_long(val);
}
#else
HEADER_INLINE uint32_t popcount01_long(uintptr_t val) {
  return popcount2_long(val);
}
#endif

HEADER_INLINE uint32_t popcount_2_longs(uintptr_t val0, uintptr_t val1) {
  // the simple version, good enough for all non-time-critical stuff
  // (without SSE4.2, popcount_longs() tends to be >3x as fast on arrays.
  // with SSE4.2, there's no noticeable difference.)
#ifdef USE_SSE42
  return __builtin_popcountll(val0) + __builtin_popcountll(val1);
#else
  val0 -= (val0 >> 1) & kMask5555;
  val1 -= (val1 >> 1) & kMask5555;
  const uintptr_t four_bit = (val0 & kMask3333) + ((val0 >> 2) & kMask3333) + (val1 & kMask3333) + ((val1 >> 2) & kMask3333);
  // up to 16 values in 0..12; sum fits in 8 bits
  return (((four_bit & kMask0F0F) + ((four_bit >> 4) & kMask0F0F)) * kMask0101) >> (kBitsPerWord - 8);
#endif
}

#ifndef __LP64__
HEADER_INLINE uint32_t popcount_4_longs(uintptr_t val0, uintptr_t val1, uintptr_t val2, uintptr_t val3) {
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

// assumes subset_mask has trailing zeroes up to the next vector boundary
void fill_interleaved_mask_vec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec);

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

HEADER_INLINE void copy_quaterarr(const uintptr_t* __restrict source_quaterarr, uint32_t quaterarr_entry_ct, uintptr_t* __restrict target_quaterarr) {
  memcpy(target_quaterarr, source_quaterarr, QUATERCT_TO_WORDCT(quaterarr_entry_ct) * kBytesPerWord);
}

void copy_bitarr_subset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uintptr_t* __restrict output_bitarr);

// may want bit past the end of subset_mask (i.e. position
// raw_quaterarr_entry_ct) to always be allocated and unset.  This removes the
// need for some explicit end-of-bitarray checks.
void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_quaterarr);

HEADER_INLINE uint32_t raw_to_subsetted_pos(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, uint32_t raw_idx) {
  // this should be much better than keeping a uidx_to_idx array!
  uint32_t raw_widx = raw_idx / kBitsPerWord;
  return subset_cumulative_popcounts[raw_widx] + popcount_long(subset_mask[raw_widx] & ((k1LU << (raw_idx % kBitsPerWord)) - k1LU));
}

void genovec_count_freqs_unsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* counts);

void genovec_count_subset_freqs(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts);

// slower genovec_count_subset_freqs() which does not require
// sample_include_interleaved_vec to be precomputed (and incidentally doesn't
// require vector alignment)
void genoarr_count_subset_freqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts);

void genoarr_count_subset_intersect_freqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, uint32_t* genocounts);

void genovec_invert_unsafe(uint32_t sample_ct, uintptr_t* genovec);

HEADER_INLINE uintptr_t invert_geno_word_unsafe(uintptr_t geno_word) {
  return (geno_word ^ ((~(geno_word << 1)) & kMaskAAAA));
}

HEADER_INLINE void zero_trailing_bits(uintptr_t bit_ct, uintptr_t* bitarr) {
  uintptr_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    bitarr[bit_ct / kBitsPerWord] &= (k1LU << trail_ct) - k1LU;
  }
}

// too easy to forget to multiply by 2
HEADER_INLINE void zero_trailing_quaters(uintptr_t quater_ct, uintptr_t* bitarr) {
  zero_trailing_bits(quater_ct * 2, bitarr);
}


HEADER_INLINE uint32_t bytes_to_represent_ui(uint32_t uii) {
  return (4 - (__builtin_clz(uii) / CHAR_BIT));
}

// A VINT is a sequence of bytes where each byte stores just 7 bits of an
// an integer, and the high bit is set when the integer has more nonzero bits.
// See e.g.
//   https://developers.google.com/protocol-buffers/docs/encoding#varints
// (Note that protocol buffers used "group varints" at one point, but then
// abandoned them.  I suspect they'd be simultaneously slower and less
// compact here.)

HEADER_INLINE unsigned char* vint32_append(uint32_t uii, unsigned char* buf) {
  while (uii > 127) {
    *buf++ = (uii & 127) + 128;
    uii >>= 7;
  }
  *buf++ = uii;
  return buf;
}

// Returns 0x80000000U instead of 0xffffffffU so overflow check works properly
// in 32-bit build.  Named "get_vint31" to make it more obvious that a 2^31
// return value can't be legitimate.
HEADER_INLINE uint32_t get_vint31(const unsigned char* buf_end, const unsigned char** bufpp) {
  if (buf_end > (*bufpp)) {
    uint32_t vint32 = *((*bufpp)++);
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (buf_end > (*bufpp)) {
      uint32_t uii = *((*bufpp)++);
      vint32 |= (uii & 127) << shift;
      if (uii <= 127) {
	return vint32;
      }
      shift += 7;
      // currently don't check for shift >= 32 (that's what validate_vint31()
      // is for).
    }
  }
  return 0x80000000U;
}

// Input must be validated, or bufp must be >= 5 characters before the end of
// the read buffer.
// todo: check if this has enough of a speed advantage over get_vint31() to
// justify using this in the main loops and catching SIGSEGV.  (update: using
// this over get_vint31() provides a ~3% speed advantage for
// load-and-recompress on the big test dataset.)
HEADER_INLINE uint32_t get_vint31_unsafe(const unsigned char** bufpp) {
  uint32_t vint32 = *(*bufpp)++;
  if (vint32 <= 127) {
    return vint32;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift < 32; shift += 7) {
    uint32_t uii = *(*bufpp)++;
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      return vint32;
    }
  }
  return 0x80000000U;
}

// Does not update buf_ptr.
HEADER_INLINE uint32_t peek_vint31(const unsigned char* buf_ptr, const unsigned char* buf_end) {
  if (buf_end > buf_ptr) {
    uint32_t vint32 = *buf_ptr++;
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (buf_end > buf_ptr) {
      uint32_t uii = *buf_ptr++;
      vint32 |= (uii & 127) << shift;
      if (uii <= 127) {
	return vint32;
      }
      shift += 7;
    }
  }
  return 0x80000000U;
}

/*
HEADER_INLINE uint32_t fget_vint31(FILE* ff) {
  // Can't be used when multiple threads are reading from ff.
  uint32_t vint32 = getc_unlocked(ff);
  if (vint32 <= 127) {
    return vint32;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift < 32; shift += 7) {
    uint32_t uii = getc_unlocked(ff);
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      return vint32;
    }
  }
  return 0x80000000U;
}

HEADER_INLINE void fput_vint31(uint32_t uii, FILE* ff) {
  // caller's responsibility to periodically check ferror
  while (uii > 127) {
    putc_unlocked((uii & 127) + 128, ff);
    uii >>= 7;
  }
  putc_unlocked(uii, ff);
}
*/

// main batch size
CONSTU31(kPglQuaterTransposeBatch, kQuatersPerCacheline);

// word width of each matrix row
CONSTU31(kPglQuaterTransposeWords, kWordsPerCacheline);

CONSTU31(kPglQuaterTransposeBufbytes, (kPglQuaterTransposeBatch * kPglQuaterTransposeBatch) / 2);
CONSTU31(kPglQuaterTransposeBufwords, kPglQuaterTransposeBufbytes / kBytesPerWord);
// up to 256x256; vecaligned_buf must have size 32k
// write_iter must be allocated up to at least
//   round_up_pow2(write_batch_size, 2) rows
void transpose_quaterblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, vul_t* vecaligned_buf);

CONSTU31(kPglBitTransposeBatch, kBitsPerCacheline);
CONSTU31(kPglBitTransposeWords, kWordsPerCacheline);
CONSTU31(kPglBitTransposeBufbytes, (kPglBitTransposeBatch * kPglBitTransposeBatch) / (CHAR_BIT / 2));
CONSTU31(kPglBitTransposeBufwords, kPglBitTransposeBufbytes / kBytesPerWord);
// up to 512x512; vecaligned_buf must have size 64k
// write_iter must be allocated up to at least
//   round_up_pow2(write_batch_size, 2) rows
void transpose_bitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, vul_t* vecaligned_buf);

// replaces each x with (32768 - x)
// okay for dosage_vals to be nullptr if dosage_ct == 0
void biallelic_dosage16_invert(uint32_t dosage_ct, uint16_t* dosage_vals);

void genovec_to_missingness_unsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict missingness);

// ----- end plink2_common subset -----

// other configuration-ish values
// this part of the specification is set in stone.

CONSTU31(kPglVblockSize, 65536);

// currently chosen so that it plus kPglFwriteBlockSize is < 2^32
static const uint32_t kPglMaxBytesPerVariant = 0xfffdffc0U;
// CONSTU31(kPglMaxBytesPerDataTrack, 0x7ffff000);
// static_assert(kMaxBytesPerIO >= (int32_t)kPglMaxBytesPerDataTrack, "pgenlib_internal assumes a single variant data track always fits in one fread/fwrite operation.");

// mmap is a horrible idea for 32-bit builds, and as long as we have non-mmap
// code we may as well not worry about Win64 CreateFileMapping.

// also, OS X mmap implementation seems to be crappy for large sequentially
// accessed files, compared to Linux.

// possible todo: SIGBUS handling?  do we ever want to try to recover from an
// I/O error?
#if defined(_WIN32) || !defined(__LP64__)
  #define NO_MMAP
#endif

// currently must be power of 2, and multiple of (kBitsPerWord / 2)
CONSTU31(kPglDifflistGroupSize, 64);

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
  #define FLAGSET_DEF_END(tname) } tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ ; \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator|(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa, tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(static_cast<uint32_t>(aa) | static_cast<uint32_t>(bb)); \
} \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator&(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa, tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(static_cast<uint32_t>(aa) & static_cast<uint32_t>(bb)); \
} \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator^(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa, tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(static_cast<uint32_t>(aa) ^ static_cast<uint32_t>(bb)); \
} \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator~(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(~static_cast<uint32_t>(aa)); \
} \
  \
struct tname { \
  tname() {} \
  \
  tname(const tname& source) : value_(source.value_) {} \
  \
  tname(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ source) : value_(static_cast<uint32_t>(source)) {} \
  \
  explicit tname(uint32_t source) : value_(source) {} \
  \
  operator tname ## _PGENLIB_INTERNAL_DO_NOT_USE__() const { \
    return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(value_); \
  } \
  \
  tname& operator|=(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ rhs) { \
    value_ |= rhs; \
    return *this; \
  } \
  \
  tname& operator&=(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ rhs) { \
    value_ &= rhs; \
    return *this; \
  } \
  \
  tname& operator^=(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ rhs) { \
    value_ ^= rhs; \
    return *this; \
  } \
  \
private: \
  uint32_t value_; \
}

  #define FLAGSET64_DEF_START() typedef enum : uint64_t {
  #define FLAGSET64_DEF_END(tname) } tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ ; \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator|(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa, tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(static_cast<uint64_t>(aa) | static_cast<uint64_t>(bb)); \
} \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator&(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa, tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(static_cast<uint64_t>(aa) & static_cast<uint64_t>(bb)); \
} \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator^(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa, tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ bb) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(static_cast<uint64_t>(aa) ^ static_cast<uint64_t>(bb)); \
} \
  \
inline tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ operator~(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ aa) { \
  return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(~static_cast<uint64_t>(aa)); \
} \
  \
struct tname { \
  tname() {} \
  \
  tname(const tname& source) : value_(source.value_) {} \
  \
  tname(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ source) : value_(static_cast<uint64_t>(source)) {} \
  \
  explicit tname(uint64_t source) : value_(source) {} \
  \
  operator tname ## _PGENLIB_INTERNAL_DO_NOT_USE__() const { \
    return static_cast<tname ## _PGENLIB_INTERNAL_DO_NOT_USE__>(value_); \
  } \
  \
  tname& operator|=(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ rhs) { \
    value_ |= rhs; \
    return *this; \
  } \
  \
  tname& operator&=(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ rhs) { \
    value_ &= rhs; \
    return *this; \
  } \
  \
  tname& operator^=(const tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ rhs) { \
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
  tname ## PGENLIB_INTERNAL_DO_NOT_USE__ALL_64_SET__ = ~(0LLU) } tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ ; \
static_assert(sizeof(tname ## _PGENLIB_INTERNAL_DO_NOT_USE__) == 8, "64-bit flagset constants are not actually uint64_ts."); \
typedef uint64_t tname

  #define ENUM_U31_DEF_START() typedef enum {
  #define ENUM_U31_DEF_END(tname) } tname ## _PGENLIB_INTERNAL_DO_NOT_USE__ ; \
typedef uint32_t tname

#endif

FLAGSET_DEF_START()
  kfPgenGlobal0,
  kfPgenGlobalLdCompressionPresent = (1 << 0),
  kfPgenGlobalDifflistOrLdPresent = (1 << 1),
  kfPgenGlobalHardcallPhasePresent = (1 << 2),
  kfPgenGlobalDosagePresent = (1 << 3),
  kfPgenGlobalDosagePhasePresent = (1 << 4),
  kfPgenGlobalAllNonref = (1 << 5)
FLAGSET_DEF_END(pgen_global_flags_t);

FLAGSET_DEF_START()
  kfPgrLdcache0,
  kfPgrLdcacheQuater = (1 << 0),
  kfPgrLdcacheDifflist = (1 << 1),
  kfPgrLdcacheAllHets = (1 << 2),
  kfPgrLdcacheRefalt1Genocounts = (1 << 3)
FLAGSET_DEF_END(pgr_ldcache_flags_t);

// difflist/LD compression should never involve more than
//   raw_sample_ct / kPglMaxDifflistLenDivisor
// entries.  (however, returned difflists can have up to twice as many entries,
// when a variant is LD-compressed and the reference variant is
// difflist-compressed.)
CONSTU31(kPglMaxDifflistLenDivisor, 8);

// threshold for using a deltalist to represent a bitarray on disk (currently
// relevant for dosage data)
CONSTU31(kPglMaxDeltalistLenDivisor, 9);

// The actual format:
// 1. 2 magic bytes 0x6c 0x1b.
//
// 2. Mode byte.
//      0x01 = plink1 variant-major.
//      0x02 = plink2 basic variant-major.  variant/sample counts in header,
//             00 = hom ref, 01 = het, 10 = hom alt, 11 = missing.  (vrtype 0)
//      0x03 = plink2 basic unphased dosage (vrtype 0x40)
//      0x04 = plink2 basic phased dosage (vrtype 0xc0)
//      These are designed to be easy to write.  Note that the dosage formats
//      require hardcalls to be stored as well; however, you can just set them
//      to all-missing and then use
//        plink2 --hard-call-threshold [...] --make-pgen
//      to populate them.
//
//      0x10 = variable-type and/or variable-length records present.
//      0x11 = mode 0x10, but with phase set information at the end of the
//             file.
//      larger values, and 0x05..0x0f, reserved for now.
//
// 3. If not plink1-format,
//    a. 4-byte # of variants; call this M.
//    b. 4-byte # of samples, call this N.
//    c. Additional 1-byte header "control" value (pgen_header_ctrl_t).  May be
//       extended in the future.
//       bits 0-3: Indicates vrtype and variant record length storage widths.
//         If bit 3 is unset, bits 0-1 store (vrec_len_byte_ct - 1), while bit
//         2 is set iff phase or dosage info is present (requiring 8 bits
//         instead of 4 bits for vrtypes).
//         If bit 3 is set, a specialized encoding is used which combines the
//         two pieces of information (reducing the overhead for files with few
//         samples).  The following encodings are currently defined:
//         1000: No difflist/LD/onebit compression, 2 bit
//               (vrec_len - ceil(sample_ct / 4))) value.  vrtype is zero if
//               the entry is zero, and 8 (multiallelic) if the record has 1-3
//               extra bytes.  Designed for single-sample files sharing a
//               single .bim-like file (note that if they don't share a .bim,
//               .bim size will dominate), but it's usable whenever there's no
//               variant where >2 samples have a rare alternate allele
//               (assuming <16 alt alleles).
//         1001: No difflist/LD/onebit compression, 4 bit
//               (vrec_len - ceil(sample_ct / 4)) value.  vrtype is zero if the
//               entry is zero, and 8 if the record has 1-15 extra bytes.
//       bits 4-5: alt allele count storage (00 = unstored, 01-11 = bytes per
//                 count)
//       bits 6-7: nonref flags info (00 = unstored, 01 = all ref/alt, 10 =
//                 never ref/alt, 11 = explicitly stored)
//       Bits 0-5 do not apply to the fixed-length modes (currently 0x02-0x04)
//       and should be zeroed out in that case.
//
// 4. If mode 0x10/0x11,
//    a. Array of 8-byte fpos values for the first variant in each vblock.
//       (Note that this suggests a way to support in-place insertions: some
//       unused space can be left between the vblocks.)
//    b. Sequence of header blocks, each containing information about
//       kPglVblockSize variants (except the last may be shorter).  All values
//       are known-width, to allow e.g. plink2 --make-pgen/--pmerge to compress
//       all variant records first, then fseek to the beginning of the output
//       file and write the header.
//         i. array of 4-bit or 1-byte vrtypes.
//        ii. array of variant record lengths (each occupying vrec_len_byte_ct
//            bytes, or 4 bits).
//       iii. if bits 4-5 of {3c} aren't 00, array of alt allele counts.
//        iv. nonref flags info, if explicitly stored
//      (this representation allows more efficient random access)
//    If mode 0x02-0x04, and nonref flags info explicitly stored, just that
//    bitarray.
//
// 5. The variant records.  See below for details.

// Difflist format (used for onebit, sparse variant, and LD compression):
//   a. [difflist_len VINT]
//   If difflist_len is zero, that's it.  Otherwise, the difflist is organized
//   into 64-element groups (the last group will usually be smaller), to make
//   extraction of e.g. a single sample less painful.  Note that with 20k
//   samples, a difflist is space-saving even with MAF 5%:
//     ~1/400 hom alt + ~38/400 het = (~39/400) * 20k
//                                  = ~1950 sample IDs.
//     that's 31 groups, requiring about 2 + 62 + 30 + 488 + 1919 = 2501 bytes
//     (can be slightly higher since a few ID deltas may be larger than 127);
//     uncompressed storage requires 5000 bytes.
//   b. [array of group start sample IDs, each of sample_id_byte_ct]
//   c. [array of 1-byte [delta segment lengths minus 63], with last entry
//      omitted]
//   d. [array of 2-bit replacement genotype values]
//   e. one "delta segment"/group: [array of [group size - 1] VINT values,
//      each indicating the difference between the current and previous sample
//      IDs; i.e. value is 1 for two consecutive samples]
//   f. [if multiallelic, array of appropriate-length disambiguation values]


// pgen_file_info_t and pgen_reader_t are the main exported "classes".
// Exported functions involving these data structure should all have
// "pgfi"/"pgr" in their names.

struct Pgen_file_info_struct {
  // ----- Header information, constant after initialization -----
  uint32_t raw_variant_ct;
  uint32_t raw_sample_ct;

  // 0 if variant records aren't all the same length.
  // If they are (e.g. PLINK 1 encoding; or vrtype bits 0-5 unset), we just
  // fseek to
  //   const_fpos_offset + const_vrec_width * ((uint64_t)variant_idx).
  uint64_t const_fpos_offset;

  uint32_t const_vrec_width;
  
  // see below.  positioned here instead of slightly later due to struct
  // packing behavior.
  uint32_t const_vrtype; // 256 for plink 1 encoding, 0xffffffffU for nonconst

  // size (raw_variant_ct + 1), so that the number of bytes of (zero-based)
  // variant n is var_fpos[n+1] - var_fpos[n].  nullptr if
  // const_vrec_width is nonzero.
  // It's not difficult to save some memory here (e.g. unless we're dealing
  // with >256 TB files, it's trivial to go from 8 bytes down to 6 bytes per
  // entry), but I doubt that's worth the trouble; let's worry about
  // O(mn)-or-worse stuff, and on-disk stuff, first.
  uint64_t* var_fpos;

  // representation type codes.
  //
  // bits 0-2:
  //   000 = Simple 2-bit encoding.
  //   100, 110, 111 = Simple difflist.  Low two bits store the base value.
  //         for multiallelic variants, if the base value is 0b11 (missing),
  //         auxiliary data track #1 only contains entries for explicitly
  //         listed 0b11 values, the rest are assumed to be actual missing
  //         data.  (101 should practically never happen--gross violation of
  //         Hardy-Weinberg equilibrium--so it's reserved for future use.)
  //   010 = Differences-from-earlier-variant encoding ("LD compression").  The
  //         last variant without this type of encoding is the base.
  //         To simplify random access logic, the first variant in each vblock
  //         is prohibited from using this encoding.
  //   011 = Inverted differences-from-earlier-variant encoding.  (This covers
  //         the case where a reference allele is "wrong".)  When decoding, the
  //         difflist should be processed first, then the entire genovec should
  //         be flipped.
  //   001 = 1-bit + difflist representation.  Suppose most calls are
  //         hom ref or het (e.g. a 20% MAF variant with ~4% hom alt1, ~36%
  //         het ref/alt1, ~64% hom ref), then the main datatrack has just the
  //         low bits of the usual 2-bit codes.  This is followed by a difflist
  //         containing the hom alt1 and missing genotypes.
  //         The main datatrack is preceded by a single byte indicating what
  //         the two common values are: 2 low bits = [set value - unset value],
  //         next 2 bits = unset value (6 possibilities).  Top 4 bits are
  //         reserved.  When the set value is 3, it does NOT represent
  //         rarealts; those must be explicitly spelled out in the difflist.
  // bit 3: more than 1 alt allele?
  // bit 4: hardcall phased?  if yes, auxiliary data track #2 contains
  //        phasing information for heterozygous calls.
  //        The first *bit* of the track indicates whether an explicit
  //        "phasepresent" bitarray is stored.  If it's set, the next het_ct
  //        bits are 1-bit values, where 0 = no phasing info known, and 1 =
  //        phasing info present.  If it's unset, phasing info is present for
  //        every het call.
  //        This is followed by a "phaseinfo" bitarray, where 0 = unswapped,
  //        1 = swapped (e.g. "1|0" in VCF).
  //        This track is normally unpacked into fixed-size bitarrays when
  //        loaded, but a raw mode is also provided (which doesn't support
  //        subsetting).
  //        By default, entire chromosomes/contigs are assumed to be phased
  //        together.  (Todo: support contiguous phase sets.)
  //
  // bits 5-6:
  //   00 = no dosage data.
  //   01 = dosage list.  auxiliary data track #3 contains a delta-encoded list
  //        of sample IDs (like a difflist, but with no genotypes).  if dosage
  //        is unphased, track #5 contains a 16-bit (0..2^15; 65535 missing
  //        value is only permitted in unconditional-dosage case) value for
  //        each allele except the last alt; if it's phased, it uses the order
  //          [hap1 ref prob] [hap2 ref prob] [hap1 alt1 prob] ...
  //        where the values are in 0..2^14 (to minimize rounding headaches).
  //        Note that this and the other dosage modes are in ADDITION to
  //        hardcalls.  This increases filesize by up to 12.5%, but makes the
  //        reader substantially simpler; --hard-call-threshold logic is nicely
  //        compartmentalized.
  //   10 = unconditional dosage (just track #5).
  //   11 = dosage bitarray.  in this case, auxiliary data track #3 contains an
  //        array of 1-bit values indicating which samples have dosages.
  //   bgen 1.2 format no longer permits fractional missingness, so no good
  //   reason for us to support it.
  //   considered putting *all* dosage data at the end of the file (like I will
  //   do for phase set info); this could actually be worthwhile for
  //   unconditional dosages, but it doesn't work well when only some samples
  //   have dosage data.
  // bit 7: some dosages phased?  if yes, and dosages are not unconditionally
  //        present, auxiliary data track #4 is either a single zero byte
  //        (indicating that all dosages are phased), or a bitarray of length
  //        (dosage_ct + 1) where the first bit is set, and the other bits
  //        indicate whether phase info is present for that sample (unset = no
  //        phasing info)
  //        note that this is independent of bit 4; either can be set without
  //        the other.
  //
  // Representation of variable ploidy (MT) was considered, but rejected since
  // dosages should be at least as appropriate for MT.
  // Oxford/VCF-style storage of separate probabilities for every possible
  // genotype (e.g. P(AA), P(AB), P(BB) instead of just 2P(AA) + P(AB) and
  // 2P(BB) + P(AB)) is tentatively rejected due to (i) lack of relevance to
  // PLINK's analysis functions and (ii) high storage cost where we can afford
  // it least.  However, this is subject to reevaluation if (i) changes.
  //
  //
  // base pointer is null if mode is 0x01-0x04 (const_vrtype != 0xffffffffU).
  // if not nullptr, required to be length >=
  //   max(raw_variant_ct + 1, round_up_pow2(raw_variant_ct, kBytesPerWord))
  unsigned char* vrtypes;

  // alt allele counts.  if >1, auxiliary data track #1 disambiguates all the
  // "missing or rare alt" explicit hardcalls.  genotype representation is:
  //   low bits: smaller [1-based alt allele idx], 0 = ref
  //   high bits: larger [1-based alt allele idx]
  //   ...
  //   2 alts: 1-bit array with 0 = missing, 1 = nonmissing.  Then, for the
  //     nonmissing subset, 2 low bits.  The high bits entry is omitted because
  //     the value has to be alt2; optimize the common case!
  //   3 alts: 1-bit nonmissingness array,  2 low bits, 2 high bits
  //   4-15 alts: 1-bit nonmissingness array, then 4 low bits, 4 high bits
  //   16-255 alts: 1-bit nonmissingness array; then 8 low bits, 8 high bits
  //   (the following is also defined, but not implemented for now:
  //   256-4095 alts: 1-bit nonmissingness array; 12 low bits, 12 high bits
  //   4096-65535 alts: 1-bit nonmissingness array; 16 low bits, 16 high bits
  //   65536-16777215 alts: 1-bit nonmissingness array; 24 low bits, 24 high
  //   bits; the latter might be necessary in the most variable regions if we
  //   use tiles...)
  // This can be nullptr if all alt allele counts are 1.
  // (actually, we store the allele index offsets, so
  // (allele_idx_offsets[n+1] - allele_idx_offsets[n]) is the number of alleles
  // for variant n.  Otherwise, we'd need another data structure to support
  // fast allele name lookup.)
  uintptr_t* allele_idx_offsets;

  uintptr_t* nonref_flags;

  // If pgr.nonref_flags is nullptr and kfPgenGlobalAllNonref is unset, all
  // reference alleles are assumed to be correct.
  pgen_global_flags_t gflags;
  
  uint32_t max_alt_allele_ct;
  uint32_t max_dosage_alt_allele_ct;

  // * nullptr if using mmap
  // * if using per-variant fread(), this is non-null during Pgen_file_info
  //   initialization, but it's then "moved" to the first Pgen_reader and set
  //   to nullptr.
  FILE* shared_ff;
  
  const unsigned char* block_base; // nullptr if using per-variant fread()
  uint64_t block_offset; // 0 for mmap
#ifndef NO_MMAP
  uint64_t file_size;
#endif
};

typedef struct Pgen_file_info_struct pgen_file_info_t;

struct Pgen_reader_struct {
  // would like to make this const, but that makes initialization really
  // annoying in C99
  struct Pgen_file_info_struct fi;
  
  // ----- Mutable state -----
  // If we don't fseek, what's the next variant we'd read?  (Still relevant
  // with mmap due to how LD decompression is implemented.)
  uint32_t fp_vidx;

  // ** per-variant fread()-only **
  FILE* ff;
  unsigned char* fread_buf;
  // ** end per-variant fread()-only **
  
  // if LD compression is present, cache the last non-LD-compressed variant
  uint32_t ldbase_vidx;

  // flags indicating which base_variant buffers are populated
  pgr_ldcache_flags_t ldbase_stypes;
  
  uint32_t ldbase_difflist_len;

  // these should be treated as private after initial allocation.
  // not currently guaranteed to have trailing zeroes.
  uintptr_t* ldbase_genovec;
  uintptr_t* ldbase_raregeno;

  // when ldbase_difflist_ids[] is initialized, element [ldbase_difflist_len]
  // must be set to sample_ct.
  uint32_t* ldbase_difflist_sample_ids;

  uintptr_t* ldbase_all_hets;
  
  // common genotype can be looked up from vrtypes[]

  uint32_t ldbase_refalt1_genocounts[4];
  
  uintptr_t* workspace_vec; // must hold raw_sample_ct entries
  
  // currently must hold (raw_sample_ct / kPglMaxDifflistLenDivisor)
  // entries; may need to double the sizes later
  // some top-level interface functions use these, so several lower-level
  // functions cannot
  uintptr_t* workspace_raregeno_vec;
  uint32_t* workspace_difflist_sample_ids;

  // must hold (raw_sample_ct / kPglMaxDifflistLenDivisor) entries
  uintptr_t* workspace_raregeno_tmp_loadbuf;
  uint32_t* workspace_difflist_sample_ids_tmp;

  uintptr_t* workspace_aux1_nonmissing_vec;
  uintptr_t* workspace_aux1_code_vec;

  uintptr_t* workspace_all_hets;

  uint32_t* workspace_ambig_sample_ids;
  uint32_t workspace_ambig_id_ct;

  uintptr_t* workspace_dosage_present;  
  uintptr_t* workspace_dosage_phased;
  
  // phase set loading (mode 0x11) unimplemented for now; should be a sequence
  // of (sample ID, [uint32_t phase set begin, set end), [set begin, set end),
  // ...).
};

typedef struct Pgen_reader_struct pgen_reader_t;

// might want this value to be typed...
CONSTU31(kPglVrtypePlink1, 256);

HEADER_INLINE uint32_t get_pgfi_vrtype(const pgen_file_info_t* pgfip, uint32_t vidx) {
  if (pgfip->vrtypes) {
    return pgfip->vrtypes[vidx];
  }
  return pgfip->const_vrtype;
}

HEADER_INLINE uint64_t get_pgfi_fpos(const pgen_file_info_t* pgfip, uintptr_t vidx) {
  if (pgfip->var_fpos) {
    return pgfip->var_fpos[vidx];
  }
  return pgfip->const_fpos_offset + pgfip->const_vrec_width * ((uint64_t)vidx);
}

HEADER_INLINE uint32_t get_pgfi_vrec_width(const pgen_file_info_t* pgfip, uint32_t vidx) {
  if (pgfip->var_fpos) {
    return (uint32_t)(pgfip->var_fpos[vidx + 1] - pgfip->var_fpos[vidx]);
  }
  return pgfip->const_vrec_width;
}

HEADER_INLINE uint32_t pgfi_is_simple_format(const pgen_file_info_t* pgfip) {
  return (pgfip->const_vrtype != 0xffffffffU);
}

HEADER_INLINE uint32_t vrtype_difflist(uint32_t vrtype) {
  return (vrtype & 4);
}

HEADER_INLINE uint32_t vrtype_ld_compressed(uint32_t vrtype) {
  return (vrtype & 6) == 2;
}

HEADER_INLINE uint32_t vrtype_multiallelic(uint32_t vrtype) {
  return (vrtype & 8);
}

HEADER_INLINE uint32_t vrtype_hphase(uint32_t vrtype) {
  return (vrtype & 0x10);
}

HEADER_INLINE uint32_t vrtype_aux_tracks_present(uint32_t vrtype) {
  return (vrtype & 0x78);
}

HEADER_INLINE uint32_t vrtype_variable_width(uint32_t vrtype) {
  return (vrtype & 0x3e);
}

HEADER_INLINE uint32_t vrtype_dosage(uint32_t vrtype) {
  return (vrtype & 0x60);
}

HEADER_INLINE uintptr_t get_aux1_allele_bytect(uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct) {
  assert(alt_allele_ct >= 2);
  if (alt_allele_ct == 2) {
    return DIV_UP(aux1_nonmissing_ct, 4);
  }
  if (alt_allele_ct == 3) {
    return DIV_UP(aux1_nonmissing_ct, 2);
  }
  // one byte per entry for alt_allele_ct < 16, two bytes for 16..255
  return ((alt_allele_ct >= 16) + 1) * aux1_nonmissing_ct;
  // todo: alt_allele_ct > 255
}

// pgen_file_info_t initialization is split into two phases, to decouple
// plink2's arena allocator from this library.
//
// Phase 1: Open the .pgen; verify that the initial bytes are consistent with
//   the file format; load/verify sample and variant counts, initialize
//   pgfi.const_vrtype, pgfi.const_vrec_width, and pgfi.const_fpos_offset;
//   determine initial memory allocation requirement.  first_alloc_cacheline_ct
//   does not include allele counts and nonref flags, since it may be more
//   appropriate to allocate those arrays earlier (during loading of a
//   .bim-like file).
//
//   pgfi.var_fpos is set to nullptr if pgfi.const_vrec_width is nonzero.
//   pgfi.vrtypes/var_allele_cts are set to nullptr in the plink1-format case.
//
//   raw_sample_ct and raw_variant_ct should be 0xffffffffU if not previously
//   known.
//
// Intermission: Caller obtains a block of pgfi_alloc_cacheline_ct * 64 bytes,
//   64-byte aligned.  The cachealigned_malloc() function can be used for this
//   purpose.  If necessary, pgfi.allele_idx_offsets and pgfi.nonref_flags
//   should be pointed at already-loaded data, or allocated so they can be
//   loaded during phase 2.
//
// Phase 2: Initialize most pointers in the pgen_reader_t struct to appropriate
//   positions in first_alloc.  For modes 0x10-0x11, load pgfi.var_fpos and
//   pgfi.vrtypes, load/validate pgfi.allele_idx_offsets and pgfi.nonref_flags
//   if appropriate, and initialize pgfi.gflags, pgfi.max_alt_allele_ct, and
//   pgfi.max_dosage_alt_allele_ct.
//
// Finally, if block-fread mode is being used, pgfi.block_base must be
//   initialized to point to a memory large enough to handle the largest
//   pgfi_block_read() operation that will be attempted.
//   pgfi_blockload_get_cacheline_req() can be used to determine the necessary
//   buffer size.

// This type may change if we introduce a more read-optimized format in the
// future.  Right now it just tracks the presence/absence of two optional
// pieces of information: allele counts and nonref flags.
typedef uint32_t pgen_header_ctrl_t;

void pgfi_preinit(pgen_file_info_t* pgfip);

// There are three modes of operation:
// 1. mmaped file.  Appropriate for handling multiple queries across different
//    parts of the genome in parallel.  Suboptimal for whole-genome queries.
//    Doesn't currently run on Windows.
// 2. fread block-load.  Block-load operations are single-threaded, while
//    decompression/counting is multithreaded.  Appropriate for whole-genome
//    queries, since even with a SSD, reading from multiple parts of a file
//    simultaneously doesn't work well.
// 3. fread single-variant-at-a-time.  Simpler interface than block-load, and
//    doesn't share its inability to handle multiple queries at a time, but
//    less performant for CPU-heavy operations on the whole genome.
//
// To specify mode 1, pass in use_mmap == 1 here.
// To specify mode 2, pass in use_mmap == 0 here, and use_blockload == 1 during
//   phase2.
// To specify mode 3, pass in use_mmap == 0 here, and use_blockload == 0 during
//   phase2.
pglerr_t pgfi_init_phase1(const char* fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, uint32_t use_mmap, pgen_header_ctrl_t* header_ctrl_ptr, pgen_file_info_t* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf);

// If allele_cts_already_loaded is set, but they're present in the file,
// they'll be validated; similarly for nonref_flags_already_loaded.
pglerr_t pgfi_init_phase2(pgen_header_ctrl_t header_ctrl, uint32_t allele_cts_already_loaded, uint32_t nonref_flags_already_loaded, uint32_t use_blockload, uint32_t vblock_idx_start, uint32_t vidx_end, uint32_t* max_vrec_width_ptr, pgen_file_info_t* pgfip, unsigned char* pgfi_alloc, uintptr_t* pgr_alloc_cacheline_ct_ptr, char* errstr_buf);


uint64_t pgfi_multiread_get_cacheline_req(const uintptr_t* variant_include, const pgen_file_info_t* pgfip, uint32_t variant_ct, uint32_t block_size);

// variant_include can be nullptr; in that case, we simply load all the
// variants (load_variant_ct must be variant_uidx_end - variant_uidx_start).)
// IMPORTANT: pgfi.block_offset must be manually copied to each reader for now.
//   (todo: probably replace pgr.fi with a pointer.  when doing that, need to
//   ensure multiple per-variant readers still works.)
pglerr_t pgfi_multiread(const uintptr_t* variant_include, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t load_variant_ct, pgen_file_info_t* pgfip);


void pgr_preinit(pgen_reader_t* pgrp);

// Before pgr_init() is called, the caller must obtain a block of
// pgr_alloc_cacheline_ct * 64 bytes (this value is returned by
// pgfi_init_phase2), 64-byte aligned; this is the pgr_alloc parameter.
//
// There's also a modal usage difference:
//
// * Modes 1-2 (mmap, block-fread): There is one pgen_file_info_t per file
//   which doesn't belong to any reader.  After it's initialized, multiple
//   pgen_reader_ts can be based off of it.  When the pgen_file_info_t is
//   destroyed, those pgen_reader_ts are invalidated and should be destroyed if
//   that hasn't already happened.
//
//   fname parameter must be nullptr.
//
// * Mode 3 (per-variant fread): Destruction of the original pgen_file_info_t
//   struct does not invalidate any extant pgen_reader_t instances (at least
//   from pgenlib_internal's perspective).  Instead, destruction of the
//   corresponding memory block or allele_idx_offsets/nonref_flags invalidates
//   the associated pgen_reader_ts.
//
//   The only difference between the first reader and later readers of the same
//   file is that the first reader steals the shared_ff used to read the
//   header.
//
//   fname parameter must be non-null.

pglerr_t pgr_init(const char* fname, uint32_t max_vrec_width, pgen_file_info_t* pgfip, pgen_reader_t* pgrp, unsigned char* pgr_alloc);

// practically all these functions require genovec to be allocated up to
// vector, not word, boundary
void pgr_plink1_to_plink2_inplace_unsafe(uint32_t sample_ct, uintptr_t* genovec);

void pgr_plink2_to_plink1_inplace_unsafe(uint32_t sample_ct, uintptr_t* genovec);

void pgr_difflist_to_genovec_unsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec);

// This will normally extract only the genotype indexes corresponding to set
// bits in sample_include.  Set sample_ct == raw_sample_ct if you don't want
// any subsetting to occur (in this case sample_include is ignored, can be
// nullptr).
// Only the maintrack is loaded.  00 = hom ref, 01 = het ref/alt1,
// 10 = hom alt1, 11 = missing or anything else.
// If multiallelic_relevant is set, and the current variant is multiallelic,
// pgr.workspace_ambig_sample_ids and pgr.workspace_ambig_id_ct are updated.
// "unsafe": sample_ct cannot be zero.  Trailing bits of genovec are not zeroed
// out.
// Ok if genovec only has space for sample_ct values.
pglerr_t pgr_read_refalt1_genovec_subset_unsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec);

// Loads the specified variant as a difflist if that's more efficient, setting
// difflist_common_geno to the common genotype value in that case.  Otherwise,
// genovec is populated and difflist_common_geno is set to 0xffffffffU.
//
// Note that the returned difflist_len can be much larger than
// max_simple_difflist_len when the variant is LD-encoded; it's bounded by
//   2 * (raw_sample_ct / kPglMaxDifflistLenDivisor).
pglerr_t pgr_read_refalt1_difflist_or_genovec_subset_unsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr);

// This is necessary when changing sample_include, unless the new query is
// iterating from the first variant.  (Which can almost never be assumed in
// plink2 since variant_include[] may not include the first variant.)
HEADER_INLINE void pgr_clear_ld_cache(pgen_reader_t* pgrp) {
  pgrp->ldbase_stypes &= kfPgrLdcacheAllHets;

  // bugfix, ld_load_necessary() was otherwise claiming that reload wasn't
  // necessary in certain cases
  pgrp->ldbase_vidx = 0x80000000U;
}

pglerr_t pgr_get_refalt1_genotype_counts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uint32_t* genocounts);

// allele_idx is set to 0 for ref, 1 for alt1, 2 for alt2, etc.
// frequencies are computed on the fly.  ties are broken in favor of the
// lower-indexed allele.
// possible todo: also provide ..._common2_then_subset() function.
// better default than the functions above for machine learning/GWAS, since the
// reference allele is "wrong" sometimes.
pglerr_t pgr_read_genovec_subset_then_common2(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec, uint32_t* __restrict maj_allele_idx_ptr, uint32_t* __restrict second_allele_idx_ptr, uint32_t* __restrict allele_ct_buf);

// Loads a quatervec with counts of a single allele (allele_idx 0 corresponds
// to the reference allele, allele_idx 1 corresponds to alt1, etc.).  0b11 ==
// missing call.
// Note that calling this with allele_idx == 0 is similar to a plink1 load
// (except with missing == 0b11, of course).
// todo: provide a difflist interface once anyone wants it.
pglerr_t pgr_read_allele_countvec_subset_unsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, pgen_reader_t* pgrp, uintptr_t* __restrict allele_countvec);

// todo: add functions which directly support MAF-based queries.  Note that
// when the difflist representation is used, we can disqualify some low-MAF
// variants without actually loading the genotype data, since the size of the
// record puts an upper bound on the alt allele frequency.

// requires trailing bits of genovec to be zeroed out, AND does not update high
// bits of last word if raw_sample_ctl2 is odd.
void detect_genovec_hets_hw(const uintptr_t*__restrict genovec, uint32_t raw_sample_ctl2, halfword_t* __restrict all_hets_hw);

// requires trailing bits of genovec to be zeroed out.
HEADER_INLINE void pgr_detect_genovec_hets_unsafe(const uintptr_t*__restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict all_hets) {
  halfword_t* all_hets_alias = (halfword_t*)all_hets;
  detect_genovec_hets_hw(genovec, raw_sample_ctl2, all_hets_alias);
  if (raw_sample_ctl2 % 2) {
    all_hets_alias[raw_sample_ctl2] = 0;
  }
}

HEADER_INLINE void pgr_detect_genovec_hets(const uintptr_t* __restrict genovec, uint32_t raw_sample_ct, uintptr_t* __restrict all_hets) {
  detect_genovec_hets_hw(genovec, QUATERCT_TO_WORDCT(raw_sample_ct), (halfword_t*)all_hets);
  zero_trailing_bits(raw_sample_ct, all_hets);
}

// pglerr_t pgr_read_refalt1_genovec_hphase_raw_unsafe(uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phaseraw, uint32_t* phasepresent_ct_ptr);

pglerr_t pgr_read_refalt1_genovec_hphase_subset_unsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr);

// ok for both dosage_present and dosage_vals to be nullptr when no dosage data
// is present
pglerr_t pgr_read_refalt1_genovec_dosage16_subset_unsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr);

pglerr_t pgr_get_ref_nonref_genotype_counts_and_dosage16s(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, double* mach_r2_ptr, uint32_t* genocounts, uint64_t* all_dosages);

// ok for both dosage_present and dosage_vals to be nullptr when no dosage data
// is present
pglerr_t pgr_read_refalt1_genovec_hphase_dosage16_subset_unsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr);

// interface used by --make-pgen, just performs basic LD/difflist decompression
// (still needs multiallelic and dosage-phase extensions)
pglerr_t pgr_read_raw(uint32_t vidx, pgen_global_flags_t read_gflags, pgen_reader_t* pgrp, uintptr_t** loadbuf_iter_ptr, unsigned char* loaded_vrtype_ptr);

pglerr_t pgr_validate(pgen_reader_t* pgrp, char* errstr_buf);

// missingness bit is set iff hardcall is not present (even if dosage info *is*
// present)
pglerr_t pgr_read_missingness(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf);

// either missingness_hc (hardcall) or missingness_dosage must be non-null
pglerr_t pgr_read_missingness_multi(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* __restrict missingness_hc, uintptr_t* __restrict missingness_dosage, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf);


// failure = kPglRetReadFail
boolerr_t pgfi_cleanup(pgen_file_info_t* pgfip);

boolerr_t pgr_cleanup(pgen_reader_t* pgrp);


struct Pgen_writer_common_struct {
  uint32_t variant_ct;
  uint32_t sample_ct;
  pgen_global_flags_t phase_dosage_gflags; // subset of gflags

  // there should be a single copy of these arrays shared by all threads.
  // allele_idx_offsets is read-only.
  uint64_t* vblock_fpos;
  unsigned char* vrec_len_buf;
  uintptr_t* vrtype_buf;
  const uintptr_t* allele_idx_offsets;
  uintptr_t* explicit_nonref_flags; // usually nullptr

  // you can get a ~1-2% speedup by writing directly to genovec and swapping
  // it with ldbase_genovec when appropriate; don't think that's worth
  // supporting, given the messier API.
  // uintptr_t* genovec;
  
  uint32_t ldbase_genocounts[4];

  // should match ftello() return value in singlethreaded case, but be set to
  // zero in multithreaded case
  uint64_t vblock_fpos_offset;
  
  // these must hold sample_ct entries (could be fewer if not subsetting, but
  // let's play it safe)
  // genovec_invert_buf also used as phaseinfo and dphase_present temporary
  // storage
  uintptr_t* genovec_invert_buf;
  uintptr_t* ldbase_genovec;
  
  // these must hold 2 * (sample_ct / kPglMaxDifflistLenDivisor) entries
  uintptr_t* ldbase_raregeno;
  uint32_t* ldbase_difflist_sample_ids; // 1 extra entry, == sample_ct
  
  // this must fit 64k variants in multithreaded case
  unsigned char* fwrite_buf;
  unsigned char* fwrite_bufp;

  uint32_t ldbase_common_geno; // 0xffffffffU if ldbase_genovec present
  uint32_t ldbase_difflist_len;
  
  // I'll cache this for now
  uintptr_t vrec_len_byte_ct;
  
  uint32_t vidx;
};

typedef struct Pgen_writer_common_struct pgen_writer_common_t;

CONSTU31(kPglFwriteBlockSize, 131072);

// Given packed arrays of unphased biallelic genotypes in uncompressed plink2
// binary format (00 = hom ref, 01 = het ref/alt1, 10 = hom alt1, 11 =
// missing), {Single,Multi}threaded_pgen_writer performs difflist (sparse
// variant), one bit (mostly-two-value), and LD compression before writing to
// disk, and backfills the header at the end.  CPRA -> CPR merging is under
// development.
// The major difference between the two interfaces is that
// Multithreaded_pgen_writer forces you to process large blocks of variants at
// a time (64k per thread).  So Singlethreaded_pgen_writer is still worth using
// in some cases (memory is very limited, I/O is slow, no programmer time to
// spare for the additional complexity).

struct Singlethreaded_pgen_writer_struct {
  struct Pgen_writer_common_struct pwc;
  FILE* pgen_outfile;  
};

struct Multithreaded_pgen_writer_struct {
  FILE* pgen_outfile;
  uint32_t thread_ct;
  struct Pgen_writer_common_struct* pwcs[];
};

typedef struct Singlethreaded_pgen_writer_struct st_pgen_writer_t;
typedef struct Multithreaded_pgen_writer_struct mt_pgen_writer_t;

void spgw_preinit(st_pgen_writer_t* spgwp);

// nonref_flags_storage values:
//   0 = no info stored
//   1 = always trusted
//   2 = always untrusted
//   3 = use explicit_nonref_flags
pglerr_t spgw_init_phase1(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, pgen_global_flags_t phase_dosage_gflags, uint32_t nonref_flags_storage, st_pgen_writer_t* spgwp, uintptr_t* alloc_cacheline_ct_ptr, uint32_t* max_vrec_len_ptr);

void spgw_init_phase2(uint32_t max_vrec_len, st_pgen_writer_t* spgwp, unsigned char* spgw_alloc);

// moderately likely that there isn't enough memory to use the maximum number
// of threads, so this returns per-thread memory requirements before forcing
// the caller to specify thread count
// (eventually should write code which falls back on st_pgen_writer_t
// when there isn't enough memory for even a single 64k variant block, at least
// for the most commonly used plink 2.0 functions)
void mpgw_init_phase1(const uintptr_t* __restrict allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, pgen_global_flags_t phase_dosage_gflags, uintptr_t* alloc_base_cacheline_ct_ptr, uint64_t* alloc_per_thread_cacheline_ct_ptr, uint32_t* vrec_len_byte_ct_ptr, uint64_t* vblock_cacheline_ct_ptr);

pglerr_t mpgw_init_phase2(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, pgen_global_flags_t phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, mt_pgen_writer_t* mpgwp);


// trailing bits of genovec must be zeroed out
void pwc_append_biallelic_genovec(const uintptr_t* __restrict genovec, pgen_writer_common_t* pwcp);

pglerr_t spgw_append_biallelic_genovec(const uintptr_t* __restrict genovec, st_pgen_writer_t* spgwp);

// trailing bits of raregeno must be zeroed out
// all raregeno entries assumed to be unequal to difflist_common_geno; the
// difflist should be compacted first if this isn't true (might be possible
// with multiallelic projections?)
// difflist_len must be <= 2 * (sample_ct / kPglMaxDifflistLenDivisor);
// there's an assert checking this
void pwc_append_biallelic_difflist_limited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, pgen_writer_common_t* pwcp);

pglerr_t spgw_append_biallelic_difflist_limited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, st_pgen_writer_t* spgwp);

// trailing bits of refalt1_genovec must be zeroed out
// not implemented yet
pglerr_t spgw_append_multiallelic_counts(const uintptr_t** __restrict alt_countvecs);

// phasepresent == nullptr ok, that indicates that ALL heterozygous calls are
// phased.  Caller should use e.g. pwc_append_biallelic_genovec() if it's known
// in advance that no calls are phased.
// Ok for phaseinfo to have bits set at non-het calls, NOT currently okay for
//   phasepresent
// void pwc_append_biallelic_genovec_hphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, pgen_writer_common_t* pwcp);

// phasepresent == nullptr ok
// ok for trailing bits of phaseinfo to not be zeroed out, NOT currently ok for
//   phasepresent
pglerr_t spgw_append_biallelic_genovec_hphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, st_pgen_writer_t* spgwp);

// dosage_vals[] has length dosage_ct, not sample_ct
void pwc_append_biallelic_genovec_dosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, pgen_writer_common_t* pwcp);

pglerr_t spgw_append_biallelic_genovec_dosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, st_pgen_writer_t* spgwp);

void pwc_append_biallelic_genovec_hphase_dosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, pgen_writer_common_t* pwcp);

pglerr_t spgw_append_biallelic_genovec_hphase_dosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, st_pgen_writer_t* spgwp);

// dphase_present can be nullptr if dosage_ct == dphase_ct
// dosage_present cannot be null for nonzero dosage_ct
// dosage_vals[] has length dosage_ct + dphase_ct
// pglerr_t spgw_append_biallelic_genovec_dphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* dphase_present, const uint16_t* dosage_vals, uint32_t dosage_ct, uint32_t dphase_ct, st_pgen_writer_t* spgwp);


// Backfills header info, then closes the file.
pglerr_t spgw_finish(st_pgen_writer_t* spgwp);

// Last flush automatically backfills header info and closes the file.
// (caller should set mpgwp = nullptr after that)
pglerr_t  mpgw_flush(mt_pgen_writer_t* mpgwp);


// these close the file if open, but do not free any memory
// mpgw_cleanup() handles mpgwp == nullptr, since it shouldn't be allocated on
// the stack
boolerr_t spgw_cleanup(st_pgen_writer_t* spgwp);
boolerr_t mpgw_cleanup(mt_pgen_writer_t* mpgwp);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PGENLIB_INTERNAL_H__
