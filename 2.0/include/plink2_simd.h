#ifndef __PLINK2_SIMD_H__
#define __PLINK2_SIMD_H__

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

// SIMD-dependent primitives which have been moved out of plink2_base due to
// the size of the SIMDe dependency.

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "plink2_base.h"

#ifdef __LP64__
// Possible todo: working no-SSE2 fallback on 64-bit little-endian platforms
// unsupported by simde.  Can perform early test by compiling on arm64 without
// simde.
// (But not yet aware of a use case that matters.)
#  define USE_SSE2
#  ifdef __x86_64__
#    include <emmintrin.h>  // IWYU pragma: export
// xmmintrin.h moved to plink2_base since we want e.g. plink2_float to have
// access to _MM_SET_FLUSH_ZERO_MODE()
#  else
#    define SIMDE_ENABLE_NATIVE_ALIASES
// Since e.g. an old zstd system header breaks the build, and plink2 is
// expected to remain under active development for the next few years, we
// currently default to using vendored copies of zstd/libdeflate/simde, which
// are manually updated as necessary.
// To use system headers, define IGNORE_BUNDLED_{ZSTD,LIBDEFLATE,SIMDE}.
#    ifdef IGNORE_BUNDLED_SIMDE
#      include <simde/x86/sse2.h>  // IWYU pragma: export
#    else
#      include "../simde/x86/sse2.h"  // IWYU pragma: export
#    endif
#    ifdef SIMDE_ARM_NEON_A32V8_NATIVE
// For Apple M1, we effectively use SSE2 + constrained _mm_shuffle_epi8().
// - We don't want to use simde's emulated _mm_shuffle_epi8 since it has an
//   extra and-with-0x8f step that we never need.
//   In the event the and-with-0x8f is actually needed, we'll define
//   vec..._x86_shuffle8() helper functions.
// - M1 also doesn't have efficient word-popcount.
#      define USE_SHUFFLE8
#    endif
#  endif
// USE_SSE42 and USE_AVX2 defined in plink2_base, since both have non-SIMD
// components.
#  ifdef USE_SSE42
#    define USE_SHUFFLE8
#    include <smmintrin.h>
#  endif
#  define ALIGNV16 __attribute__ ((aligned (16)))
#else
#  define ALIGNV16
// As of Feb 2026, we are dropping support for 32-bit processors without SSE.
#  include <xmmintrin.h>  // IWYU pragma: export
#endif

// done with #includes, can start C++ namespace...
#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef USE_SSE2
#  ifdef USE_AVX2
CONSTI32(kBytesPerVec, 32);

// 16 still seems to noticeably outperform 32 on my Mac test machine, and
// is about equal on my Linux test machine, probably due to reduced clock
// frequency when 32-byte floating point vector operations are used (as in, ALL
// operations, sometimes on ALL cores, become slower when a single core
// performs a 32-byte fp vector operation).
// However, processor power management, numeric libraries, and my AVX2 code
// should improve over time.  There will probably come a time where switching
// to 32-byte fp is worthwhile.
#    define FVEC_32

// bleah, have to define these here, vector_size doesn't see enum values
typedef uintptr_t VecW __attribute__ ((vector_size (32)));
typedef uint32_t VecU32 __attribute__ ((vector_size (32)));
typedef int32_t VecI32 __attribute__ ((vector_size (32)));
typedef unsigned short VecU16 __attribute__ ((vector_size (32)));
typedef short VecI16 __attribute__ ((vector_size (32)));
// documentation says 'char', but int8_t works fine under gcc 4.4 and conveys
// intent better (char not guaranteed to be signed)
typedef int8_t VecI8 __attribute__ ((vector_size (32)));
typedef unsigned char VecUc __attribute__ ((vector_size (32)));

HEADER_INLINE VecW VecToW(__m256i vv) {
  return R_CAST(VecW, vv);
}

HEADER_INLINE VecU32 VecToU32(__m256i vv) {
  return R_CAST(VecU32, vv);
}

HEADER_INLINE VecI32 VecToI32(__m256i vv) {
  return R_CAST(VecI32, vv);
}

HEADER_INLINE VecU16 VecToU16(__m256i vv) {
  return R_CAST(VecU16, vv);
}

HEADER_INLINE VecI16 VecToI16(__m256i vv) {
  return R_CAST(VecI16, vv);
}

HEADER_INLINE VecUc VecToUc(__m256i vv) {
  return R_CAST(VecUc, vv);
}

HEADER_INLINE VecI8 VecToI8(__m256i vv) {
  return R_CAST(VecI8, vv);
}

HEADER_INLINE __m256i WToVec(VecW vv) {
  return R_CAST(__m256i, vv);
}

HEADER_INLINE __m256i U32ToVec(VecU32 vv) {
  return R_CAST(__m256i, vv);
}

HEADER_INLINE __m256i I32ToVec(VecI32 vv) {
  return R_CAST(__m256i, vv);
}

HEADER_INLINE __m256i U16ToVec(VecU16 vv) {
  return R_CAST(__m256i, vv);
}

HEADER_INLINE __m256i I16ToVec(VecI16 vv) {
  return R_CAST(__m256i, vv);
}

HEADER_INLINE __m256i UcToVec(VecUc vv) {
  return R_CAST(__m256i, vv);
}

HEADER_INLINE __m256i I8ToVec(VecI8 vv) {
  return R_CAST(__m256i, vv);
}

// _mm256_set1_... seems to have the same performance; could use that instead.
#    define VCONST_W(xx) {xx, xx, xx, xx}
#    define VCONST_S(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_C(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_UC VCONST_C

// vv = VCONST_W(k0LU) doesn't work (only ok for initialization)
HEADER_INLINE VecW vecw_setzero() {
  return VecToW(_mm256_setzero_si256());
}

HEADER_INLINE VecU32 vecu32_setzero() {
  return VecToU32(_mm256_setzero_si256());
}

HEADER_INLINE VecU16 vecu16_setzero() {
  return VecToU16(_mm256_setzero_si256());
}

HEADER_INLINE VecI16 veci16_setzero() {
  return VecToI16(_mm256_setzero_si256());
}

HEADER_INLINE VecUc vecuc_setzero() {
  return VecToUc(_mm256_setzero_si256());
}

HEADER_INLINE VecI8 veci8_setzero() {
  return VecToI8(_mm256_setzero_si256());
}

// "vv >> ct" doesn't work, and Scientific Linux gcc 4.4 might not optimize
// VCONST_W shift properly (todo: test this)
HEADER_INLINE VecW vecw_srli(VecW vv, uint32_t ct) {
  return VecToW(_mm256_srli_epi64(WToVec(vv), ct));
}

HEADER_INLINE VecW vecw_slli(VecW vv, uint32_t ct) {
  return VecToW(_mm256_slli_epi64(WToVec(vv), ct));
}

HEADER_INLINE VecU32 vecu32_srli(VecU32 vv, uint32_t ct) {
  return VecToU32(_mm256_srli_epi32(U32ToVec(vv), ct));
}

HEADER_INLINE VecU32 vecu32_slli(VecU32 vv, uint32_t ct) {
  return VecToU32(_mm256_slli_epi32(U32ToVec(vv), ct));
}

HEADER_INLINE VecU16 vecu16_srli(VecU16 vv, uint32_t ct) {
  return VecToU16(_mm256_srli_epi16(U16ToVec(vv), ct));
}

HEADER_INLINE VecU16 vecu16_slli(VecU16 vv, uint32_t ct) {
  return VecToU16(_mm256_slli_epi16(U16ToVec(vv), ct));
}

// Compiler still doesn't seem to be smart enough to use andnot properly.
HEADER_INLINE VecW vecw_and_notfirst(VecW excl, VecW main) {
  return VecToW(_mm256_andnot_si256(WToVec(excl), WToVec(main)));
}

HEADER_INLINE VecU32 vecu32_and_notfirst(VecU32 excl, VecU32 main) {
  return VecToU32(_mm256_andnot_si256(U32ToVec(excl), U32ToVec(main)));
}

HEADER_INLINE VecI32 veci32_and_notfirst(VecI32 excl, VecI32 main) {
  return VecToI32(_mm256_andnot_si256(I32ToVec(excl), I32ToVec(main)));
}

HEADER_INLINE VecU16 vecu16_and_notfirst(VecU16 excl, VecU16 main) {
  return VecToU16(_mm256_andnot_si256(U16ToVec(excl), U16ToVec(main)));
}

HEADER_INLINE VecI16 veci16_and_notfirst(VecI16 excl, VecI16 main) {
  return VecToI16(_mm256_andnot_si256(I16ToVec(excl), I16ToVec(main)));
}

HEADER_INLINE VecUc vecuc_and_notfirst(VecUc excl, VecUc main) {
  return VecToUc(_mm256_andnot_si256(UcToVec(excl), UcToVec(main)));
}

HEADER_INLINE VecI8 veci8_and_notfirst(VecI8 excl, VecI8 main) {
  return VecToI8(_mm256_andnot_si256(I8ToVec(excl), I8ToVec(main)));
}

HEADER_INLINE VecW vecw_set1(uintptr_t ulii) {
  return VecToW(_mm256_set1_epi64x(ulii));
}

HEADER_INLINE VecU32 vecu32_set1(uint32_t uii) {
  return VecToU32(_mm256_set1_epi32(uii));
}

HEADER_INLINE VecI32 veci32_set1(int32_t ii) {
  return VecToI32(_mm256_set1_epi32(ii));
}

HEADER_INLINE VecU16 vecu16_set1(unsigned short usi) {
  return VecToU16(_mm256_set1_epi16(usi));
}

HEADER_INLINE VecI16 veci16_set1(short si) {
  return VecToI16(_mm256_set1_epi16(si));
}

HEADER_INLINE VecUc vecuc_set1_epi16(unsigned short usi) {
  return VecToUc(_mm256_set1_epi16(usi));
}

HEADER_INLINE VecUc vecuc_set1(unsigned char ucc) {
  return VecToUc(_mm256_set1_epi8(ucc));
}

HEADER_INLINE VecI8 veci8_set1(char cc) {
  return VecToI8( _mm256_set1_epi8(cc));
}

HEADER_INLINE uint32_t vecw_movemask(VecW vv) {
  return _mm256_movemask_epi8(WToVec(vv));
}

HEADER_INLINE uint32_t vecu32_movemask(VecU32 vv) {
  return _mm256_movemask_epi8(U32ToVec(vv));
}

HEADER_INLINE uint32_t veci32_movemask(VecI32 vv) {
  return _mm256_movemask_epi8(I32ToVec(vv));
}

HEADER_INLINE uint32_t vecu16_movemask(VecU16 vv) {
  return _mm256_movemask_epi8(U16ToVec(vv));
}

HEADER_INLINE uint32_t veci16_movemask(VecI16 vv) {
  return _mm256_movemask_epi8(I16ToVec(vv));
}

HEADER_INLINE uint32_t veci8_movemask(VecI8 vv) {
  return _mm256_movemask_epi8(I8ToVec(vv));
}

HEADER_INLINE uint32_t vecuc_movemask(VecUc vv) {
  return _mm256_movemask_epi8(UcToVec(vv));
}

// Repeats elements in second lane in AVX2 case.
HEADER_INLINE VecW vecw_setr8(char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToW(_mm256_setr_epi8(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0, e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

HEADER_INLINE VecU16 vecu16_setr8(char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToU16(_mm256_setr_epi8(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0, e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

HEADER_INLINE VecUc vecuc_setr8(char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToUc(_mm256_setr_epi8(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0, e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

// Discards last 16 arguments in SSE2/SSE4.2 case.
HEADER_INLINE VecW vecw_setr8x(char e31, char e30, char e29, char e28, char e27, char e26, char e25, char e24, char e23, char e22, char e21, char e20, char e19, char e18, char e17, char e16, char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToW(_mm256_setr_epi8(e31, e30, e29, e28, e27, e26, e25, e24, e23, e22, e21, e20, e19, e18, e17, e16, e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

HEADER_INLINE VecUc vecuc_setr8x(char e31, char e30, char e29, char e28, char e27, char e26, char e25, char e24, char e23, char e22, char e21, char e20, char e19, char e18, char e17, char e16, char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToUc(_mm256_setr_epi8(e31, e30, e29, e28, e27, e26, e25, e24, e23, e22, e21, e20, e19, e18, e17, e16, e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

HEADER_INLINE VecW vecw_setr32(uint32_t e3, uint32_t e2, uint32_t e1, uint32_t e0) {
  return VecToW(_mm256_setr_epi32(e3, e2, e1, e0, e3, e2, e1, e0));
}


HEADER_INLINE VecW vecw_unpacklo8(VecW evens, VecW odds) {
  return VecToW(_mm256_unpacklo_epi8(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi8(VecW evens, VecW odds) {
  return VecToW(_mm256_unpackhi_epi8(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecI8 veci8_unpacklo8(VecI8 evens, VecI8 odds) {
  return VecToI8(_mm256_unpacklo_epi8(I8ToVec(evens), I8ToVec(odds)));
}

HEADER_INLINE VecI8 veci8_unpackhi8(VecI8 evens, VecI8 odds) {
  return VecToI8(_mm256_unpackhi_epi8(I8ToVec(evens), I8ToVec(odds)));
}

HEADER_INLINE VecUc vecuc_unpacklo8(VecUc evens, VecUc odds) {
  return VecToUc(_mm256_unpacklo_epi8(UcToVec(evens), UcToVec(odds)));
}

HEADER_INLINE VecUc vecuc_unpackhi8(VecUc evens, VecUc odds) {
  return VecToUc(_mm256_unpackhi_epi8(UcToVec(evens), UcToVec(odds)));
}

HEADER_INLINE VecW vecw_unpacklo16(VecW evens, VecW odds) {
  return VecToW(_mm256_unpacklo_epi16(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi16(VecW evens, VecW odds) {
  return VecToW(_mm256_unpackhi_epi16(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpacklo32(VecW evens, VecW odds) {
  return VecToW(_mm256_unpacklo_epi32(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi32(VecW evens, VecW odds) {
  return VecToW(_mm256_unpackhi_epi32(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_permute0xd8_if_avx2(VecW vv) {
  return VecToW(_mm256_permute4x64_epi64(WToVec(vv), 0xd8));
}

HEADER_INLINE VecI8 veci8_permute0xd8_if_avx2(VecI8 vv) {
  return VecToI8(_mm256_permute4x64_epi64(I8ToVec(vv), 0xd8));
}

HEADER_INLINE VecUc vecuc_permute0xd8_if_avx2(VecUc vv) {
  return VecToUc(_mm256_permute4x64_epi64(UcToVec(vv), 0xd8));
}

// Could have a single-src gather_even function, but that should wait until
// there is a clear SSE2 use case.
HEADER_INLINE VecW vecw_gather_even(VecW src_lo, VecW src_hi, VecW m8) {
  const VecW gathered_laneswapped = VecToW(_mm256_packus_epi16(WToVec(src_lo & m8), WToVec(src_hi & m8)));
  return vecw_permute0xd8_if_avx2(gathered_laneswapped);
}

HEADER_INLINE VecUc vecuc_gather_even(VecUc src_lo, VecUc src_hi, VecUc m8) {
  const VecUc gathered_laneswapped = VecToUc(_mm256_packus_epi16(UcToVec(src_lo & m8), UcToVec(src_hi & m8)));
  return vecuc_permute0xd8_if_avx2(gathered_laneswapped);
}

HEADER_INLINE VecUc vecuc_gather_odd(VecUc src_lo, VecUc src_hi) {
  const VecUc gathered_laneswapped = VecToUc(_mm256_packus_epi16(_mm256_srli_epi16(UcToVec(src_lo), 8), _mm256_srli_epi16(UcToVec(src_hi), 8)));
  return vecuc_permute0xd8_if_avx2(gathered_laneswapped);
}

HEADER_INLINE VecW vecw_shuffle8(VecW table, VecW indexes) {
  return VecToW(_mm256_shuffle_epi8(WToVec(table), WToVec(indexes)));
}

HEADER_INLINE VecU16 vecu16_shuffle8(VecU16 table, VecU16 indexes) {
  return VecToU16(_mm256_shuffle_epi8(U16ToVec(table), U16ToVec(indexes)));
}

HEADER_INLINE VecUc vecuc_shuffle8(VecUc table, VecUc indexes) {
  return VecToUc(_mm256_shuffle_epi8(UcToVec(table), UcToVec(indexes)));
}

HEADER_INLINE uintptr_t vecw_extract64_0(VecW vv) {
  return _mm256_extract_epi64(WToVec(vv), 0);
}

HEADER_INLINE uintptr_t vecw_extract64_1(VecW vv) {
  return _mm256_extract_epi64(WToVec(vv), 1);
}

// *** AVX2-only section ***
HEADER_INLINE uintptr_t vecw_extract64_2(VecW vv) {
  return _mm256_extract_epi64(WToVec(vv), 2);
}

HEADER_INLINE uintptr_t vecw_extract64_3(VecW vv) {
  return _mm256_extract_epi64(WToVec(vv), 3);
}

// todo: permute

// *** end AVX2-only section ***

#    define kVec8thUintMax UINT32_MAX

typedef uint16_t Vec16thUint;
typedef uint32_t Vec8thUint;
typedef uint64_t Vec4thUint;

HEADER_INLINE VecW vecw_load(const void* mem_addr) {
  return VecToW(_mm256_load_si256(S_CAST(const __m256i*, mem_addr)));
}

// There may be some value in adding a 4-consecutive-vector load function when
// addresses are expected to be unaligned: see
//   https://www.agner.org/optimize/blog/read.php?i=627&v=t

HEADER_INLINE VecW vecw_loadu(const void* mem_addr) {
  return VecToW(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE VecU32 vecu32_loadu(const void* mem_addr) {
  return VecToU32(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE VecI32 veci32_loadu(const void* mem_addr) {
  return VecToI32(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE VecU16 vecu16_loadu(const void* mem_addr) {
  return VecToU16(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE VecI16 veci16_loadu(const void* mem_addr) {
  return VecToI16(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE VecUc vecuc_loadu(const void* mem_addr) {
  return VecToUc(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE VecI8 veci8_loadu(const void* mem_addr) {
  return VecToI8(_mm256_loadu_si256(S_CAST(const __m256i*, mem_addr)));
}

HEADER_INLINE void vec_storeu(void* mem_addr, __m256i vv) {
  _mm256_storeu_si256(S_CAST(__m256i*, mem_addr), vv);
}

HEADER_INLINE VecI32 veci32_max(VecI32 v1, VecI32 v2) {
  return VecToI32(_mm256_max_epi32(I32ToVec(v1), I32ToVec(v2)));
}

HEADER_INLINE VecI16 veci16_max(VecI16 v1, VecI16 v2) {
  return VecToI16(_mm256_max_epi16(I16ToVec(v1), I16ToVec(v2)));
}

HEADER_INLINE VecW vecw_sad(VecW v1, VecW v2) {
  return VecToW(_mm256_sad_epu8(WToVec(v1), WToVec(v2)));
}

HEADER_INLINE VecUc vecuc_add(VecUc v1, VecUc v2) {
  return VecToUc(_mm256_add_epi8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecUc vecuc_adds(VecUc v1, VecUc v2) {
  return VecToUc(_mm256_adds_epu8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecUc vecuc_signed_cmpgt(VecUc v1, VecUc v2) {
  return VecToUc(_mm256_cmpgt_epi8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecU16 vecu16_min8(VecU16 v1, VecU16 v2) {
  return VecToU16(_mm256_min_epu8(U16ToVec(v1), U16ToVec(v2)));
}

HEADER_INLINE VecUc vecuc_min(VecUc v1, VecUc v2) {
  return VecToUc(_mm256_min_epu8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecW vecw_blendv(VecW aa, VecW bb, VecW mask) {
  return VecToW(_mm256_blendv_epi8(WToVec(aa), WToVec(bb), WToVec(mask)));
}

HEADER_INLINE VecU32 vecu32_blendv(VecU32 aa, VecU32 bb, VecU32 mask) {
  return VecToU32(_mm256_blendv_epi8(U32ToVec(aa), U32ToVec(bb), U32ToVec(mask)));
}

HEADER_INLINE VecU16 vecu16_blendv(VecU16 aa, VecU16 bb, VecU16 mask) {
  return VecToU16(_mm256_blendv_epi8(U16ToVec(aa), U16ToVec(bb), U16ToVec(mask)));
}

HEADER_INLINE VecUc vecuc_blendv(VecUc aa, VecUc bb, VecUc mask) {
  return VecToUc(_mm256_blendv_epi8(UcToVec(aa), UcToVec(bb), UcToVec(mask)));
}

#  else  // USE_SSE2, !USE_AVX2

CONSTI32(kBytesPerVec, 16);
typedef uintptr_t VecW __attribute__ ((vector_size (16)));
typedef uint32_t VecU32 __attribute ((vector_size (16)));
typedef int32_t VecI32 __attribute ((vector_size (16)));
typedef unsigned short VecU16 __attribute__ ((vector_size (16)));
typedef short VecI16 __attribute__ ((vector_size (16)));
typedef int8_t VecI8 __attribute__ ((vector_size (16)));
typedef unsigned char VecUc __attribute__ ((vector_size (16)));

HEADER_INLINE VecW VecToW(__m128i vv) {
  return R_CAST(VecW, vv);
}

HEADER_INLINE VecU32 VecToU32(__m128i vv) {
  return R_CAST(VecU32, vv);
}

HEADER_INLINE VecI32 VecToI32(__m128i vv) {
  return R_CAST(VecI32, vv);
}

HEADER_INLINE VecU16 VecToU16(__m128i vv) {
  return R_CAST(VecU16, vv);
}

HEADER_INLINE VecI16 VecToI16(__m128i vv) {
  return R_CAST(VecI16, vv);
}

HEADER_INLINE VecUc VecToUc(__m128i vv) {
  return R_CAST(VecUc, vv);
}

HEADER_INLINE VecI8 VecToI8(__m128i vv) {
  return R_CAST(VecI8, vv);
}

HEADER_INLINE __m128i WToVec(VecW vv) {
  return R_CAST(__m128i, vv);
}

HEADER_INLINE __m128i U32ToVec(VecU32 vv) {
  return R_CAST(__m128i, vv);
}

HEADER_INLINE __m128i I32ToVec(VecI32 vv) {
  return R_CAST(__m128i, vv);
}

HEADER_INLINE __m128i U16ToVec(VecU16 vv) {
  return R_CAST(__m128i, vv);
}

HEADER_INLINE __m128i I16ToVec(VecI16 vv) {
  return R_CAST(__m128i, vv);
}

HEADER_INLINE __m128i UcToVec(VecUc vv) {
  return R_CAST(__m128i, vv);
}

HEADER_INLINE __m128i I8ToVec(VecI8 vv) {
  return R_CAST(__m128i, vv);
}

#    define VCONST_W(xx) {xx, xx}
#    define VCONST_S(xx) {xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_C(xx) {xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_UC VCONST_C

HEADER_INLINE VecW vecw_setzero() {
  return VecToW(_mm_setzero_si128());
}

HEADER_INLINE VecU32 vecu32_setzero() {
  return VecToU32(_mm_setzero_si128());
}

HEADER_INLINE VecU16 vecu16_setzero() {
  return VecToU16(_mm_setzero_si128());
}

HEADER_INLINE VecI16 veci16_setzero() {
  return VecToI16(_mm_setzero_si128());
}

HEADER_INLINE VecUc vecuc_setzero() {
  return VecToUc(_mm_setzero_si128());
}

HEADER_INLINE VecI8 veci8_setzero() {
  return VecToI8(_mm_setzero_si128());
}

// simde is incompatible with defining these as inline functions
#    define vecw_srli(vv, ct) VecToW(_mm_srli_epi64(WToVec(vv), ct))

#    define vecw_slli(vv, ct) VecToW(_mm_slli_epi64(WToVec(vv), ct))

#    define vecu32_srli(vv, ct) VecToU32(_mm_srli_epi32(U32ToVec(vv), ct))

#    define vecu32_slli(vv, ct) VecToU32(_mm_slli_epi32(U32ToVec(vv), ct))

#    define vecu16_srli(vv, ct) VecToU16(_mm_srli_epi16(U16ToVec(vv), ct))

#    define vecu16_slli(vv, ct) VecToU16(_mm_slli_epi16(U16ToVec(vv), ct))

HEADER_INLINE VecW vecw_and_notfirst(VecW excl, VecW main) {
  return VecToW(_mm_andnot_si128(WToVec(excl), WToVec(main)));
}

HEADER_INLINE VecU32 vecu32_and_notfirst(VecU32 excl, VecU32 main) {
  return VecToU32(_mm_andnot_si128(U32ToVec(excl), U32ToVec(main)));
}

HEADER_INLINE VecI32 veci32_and_notfirst(VecI32 excl, VecI32 main) {
  return VecToI32(_mm_andnot_si128(I32ToVec(excl), I32ToVec(main)));
}

HEADER_INLINE VecU16 vecu16_and_notfirst(VecU16 excl, VecU16 main) {
  return VecToU16(_mm_andnot_si128(U16ToVec(excl), U16ToVec(main)));
}

HEADER_INLINE VecI16 veci16_and_notfirst(VecI16 excl, VecI16 main) {
  return VecToI16(_mm_andnot_si128(I16ToVec(excl), I16ToVec(main)));
}

HEADER_INLINE VecUc vecuc_and_notfirst(VecUc excl, VecUc main) {
  return VecToUc(_mm_andnot_si128(UcToVec(excl), UcToVec(main)));
}

HEADER_INLINE VecI8 veci8_and_notfirst(VecI8 excl, VecI8 main) {
  return VecToI8(_mm_andnot_si128(I8ToVec(excl), I8ToVec(main)));
}

HEADER_INLINE VecW vecw_set1(uintptr_t ulii) {
  return VecToW(_mm_set1_epi64x(ulii));
}

HEADER_INLINE VecU32 vecu32_set1(uint32_t uii) {
  return VecToU32(_mm_set1_epi32(uii));
}

HEADER_INLINE VecI32 veci32_set1(int32_t ii) {
  return VecToI32(_mm_set1_epi32(ii));
}

HEADER_INLINE VecU16 vecu16_set1(unsigned short usi) {
  return VecToU16(_mm_set1_epi16(usi));
}

HEADER_INLINE VecI16 veci16_set1(short si) {
  return VecToI16(_mm_set1_epi16(si));
}

HEADER_INLINE VecUc vecuc_set1_epi16(unsigned short usi) {
  return VecToUc(_mm_set1_epi16(usi));
}

HEADER_INLINE VecUc vecuc_set1(unsigned char ucc) {
  return VecToUc(_mm_set1_epi8(ucc));
}

HEADER_INLINE VecI8 veci8_set1(char cc) {
  return VecToI8(_mm_set1_epi8(cc));
}

HEADER_INLINE uint32_t vecw_movemask(VecW vv) {
  return _mm_movemask_epi8(WToVec(vv));
}

HEADER_INLINE uint32_t vecu32_movemask(VecU32 vv) {
  return _mm_movemask_epi8(U32ToVec(vv));
}

HEADER_INLINE uint32_t veci32_movemask(VecI32 vv) {
  return _mm_movemask_epi8(I32ToVec(vv));
}

HEADER_INLINE uint32_t vecu16_movemask(VecU16 vv) {
  return _mm_movemask_epi8(U16ToVec(vv));
}

HEADER_INLINE uint32_t veci16_movemask(VecI16 vv) {
  return _mm_movemask_epi8(I16ToVec(vv));
}

HEADER_INLINE uint32_t veci8_movemask(VecI8 vv) {
  return _mm_movemask_epi8(I8ToVec(vv));
}

HEADER_INLINE uint32_t vecuc_movemask(VecUc vv) {
  return _mm_movemask_epi8(UcToVec(vv));
}

CONSTI32(kVec8thUintMax, 65535);

// #    define kVec8thUintMax 65535

typedef unsigned char Vec16thUint;
typedef uint16_t Vec8thUint;
typedef uint32_t Vec4thUint;

HEADER_INLINE VecW vecw_load(const void* mem_addr) {
  return VecToW(_mm_load_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecW vecw_loadu(const void* mem_addr) {
  return VecToW(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecU32 vecu32_loadu(const void* mem_addr) {
  return VecToU32(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecI32 veci32_loadu(const void* mem_addr) {
  return VecToI32(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecU16 vecu16_loadu(const void* mem_addr) {
  return VecToU16(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecI16 veci16_loadu(const void* mem_addr) {
  return VecToI16(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecUc vecuc_loadu(const void* mem_addr) {
  return VecToUc(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE VecI8 veci8_loadu(const void* mem_addr) {
  return VecToI8(_mm_loadu_si128(S_CAST(const __m128i*, mem_addr)));
}

HEADER_INLINE void vec_storeu(void* mem_addr, __m128i vv) {
  _mm_storeu_si128(S_CAST(__m128i*, mem_addr), vv);
}

// Repeats arguments in AVX2 case.
HEADER_INLINE VecW vecw_setr8(char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToW(_mm_setr_epi8(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

HEADER_INLINE VecU16 vecu16_setr8(char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToU16(_mm_setr_epi8(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

HEADER_INLINE VecUc vecuc_setr8(char e15, char e14, char e13, char e12, char e11, char e10, char e9, char e8, char e7, char e6, char e5, char e4, char e3, char e2, char e1, char e0) {
  return VecToUc(_mm_setr_epi8(e15, e14, e13, e12, e11, e10, e9, e8, e7, e6, e5, e4, e3, e2, e1, e0));
}

// Discards last 16 arguments in SSE2/SSE4.2 case.
HEADER_INLINE VecW vecw_setr8x(
    char e31, char e30, char e29, char e28,
    char e27, char e26, char e25, char e24,
    char e23, char e22, char e21, char e20,
    char e19, char e18, char e17, char e16,
    __maybe_unused char e15, __maybe_unused char e14,
    __maybe_unused char e13, __maybe_unused char e12,
    __maybe_unused char e11, __maybe_unused char e10,
    __maybe_unused char e9, __maybe_unused char e8,
    __maybe_unused char e7, __maybe_unused char e6,
    __maybe_unused char e5, __maybe_unused char e4,
    __maybe_unused char e3, __maybe_unused char e2,
    __maybe_unused char e1, __maybe_unused char e0) {
  return VecToW(_mm_setr_epi8(e31, e30, e29, e28, e27, e26, e25, e24, e23, e22, e21, e20, e19, e18, e17, e16));
}

HEADER_INLINE VecUc vecuc_setr8x(
    char e31, char e30, char e29, char e28,
    char e27, char e26, char e25, char e24,
    char e23, char e22, char e21, char e20,
    char e19, char e18, char e17, char e16,
    __maybe_unused char e15, __maybe_unused char e14,
    __maybe_unused char e13, __maybe_unused char e12,
    __maybe_unused char e11, __maybe_unused char e10,
    __maybe_unused char e9, __maybe_unused char e8,
    __maybe_unused char e7, __maybe_unused char e6,
    __maybe_unused char e5, __maybe_unused char e4,
    __maybe_unused char e3, __maybe_unused char e2,
    __maybe_unused char e1, __maybe_unused char e0) {
  return VecToUc(_mm_setr_epi8(e31, e30, e29, e28, e27, e26, e25, e24, e23, e22, e21, e20, e19, e18, e17, e16));
}

HEADER_INLINE VecW vecw_setr32(uint32_t e3, uint32_t e2, uint32_t e1, uint32_t e0) {
  return VecToW(_mm_setr_epi32(e3, e2, e1, e0));
}

HEADER_INLINE VecW vecw_unpacklo8(VecW evens, VecW odds) {
  return VecToW(_mm_unpacklo_epi8(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi8(VecW evens, VecW odds) {
  return VecToW(_mm_unpackhi_epi8(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecI8 veci8_unpacklo8(VecI8 evens, VecI8 odds) {
  return VecToI8(_mm_unpacklo_epi8(I8ToVec(evens), I8ToVec(odds)));
}

HEADER_INLINE VecI8 veci8_unpackhi8(VecI8 evens, VecI8 odds) {
  return VecToI8(_mm_unpackhi_epi8(I8ToVec(evens), I8ToVec(odds)));
}

HEADER_INLINE VecUc vecuc_unpacklo8(VecUc evens, VecUc odds) {
  return VecToUc(_mm_unpacklo_epi8(UcToVec(evens), UcToVec(odds)));
}

HEADER_INLINE VecUc vecuc_unpackhi8(VecUc evens, VecUc odds) {
  return VecToUc(_mm_unpackhi_epi8(UcToVec(evens), UcToVec(odds)));
}

HEADER_INLINE VecW vecw_unpacklo16(VecW evens, VecW odds) {
  return VecToW(_mm_unpacklo_epi16(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi16(VecW evens, VecW odds) {
  return VecToW(_mm_unpackhi_epi16(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpacklo32(VecW evens, VecW odds) {
  return VecToW(_mm_unpacklo_epi32(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi32(VecW evens, VecW odds) {
  return VecToW(_mm_unpackhi_epi32(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpacklo64(VecW evens, VecW odds) {
  return VecToW(_mm_unpacklo_epi64(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_unpackhi64(VecW evens, VecW odds) {
  return VecToW(_mm_unpackhi_epi64(WToVec(evens), WToVec(odds)));
}

HEADER_INLINE VecW vecw_permute0xd8_if_avx2(VecW vv) {
  return vv;
}

HEADER_INLINE VecI8 veci8_permute0xd8_if_avx2(VecI8 vv) {
  return vv;
}

HEADER_INLINE VecUc vecuc_permute0xd8_if_avx2(VecUc vv) {
  return vv;
}

HEADER_INLINE VecW vecw_gather_even(VecW src_lo, VecW src_hi, VecW m8) {
  return VecToW(_mm_packus_epi16(WToVec(src_lo & m8), WToVec(src_hi & m8)));
}

HEADER_INLINE VecUc vecuc_gather_even(VecUc src_lo, VecUc src_hi, VecUc m8) {
  return VecToUc(_mm_packus_epi16(UcToVec(src_lo & m8), UcToVec(src_hi & m8)));
}

HEADER_INLINE VecUc vecuc_gather_odd(VecUc src_lo, VecUc src_hi) {
  return VecToUc(_mm_packus_epi16(_mm_srli_epi16(UcToVec(src_lo), 8), _mm_srli_epi16(UcToVec(src_hi), 8)));
}

#    ifdef USE_SHUFFLE8
#      ifdef SIMDE_ARM_NEON_A64V8_NATIVE
// See simde_mm_shuffle_epi8().
HEADER_INLINE __m128i _mm_shuffle_epi8(__m128i a, __m128i b) {
  SIMDE_ALIGN_TO_16 int8x16_t a_;
  SIMDE_ALIGN_TO_16 uint8x16_t b_;
  SIMDE_ALIGN_TO_16 int8x16_t r_;
  memcpy(&a_, &a, sizeof(a_));
  memcpy(&b_, &b, sizeof(b_));
  r_ = vqtbl1q_s8(a_, b_);
  __m128i r;
  memcpy(&r, &r_, sizeof(r));
  return r;
}
#      endif
HEADER_INLINE VecW vecw_shuffle8(VecW table, VecW indexes) {
  return VecToW(_mm_shuffle_epi8(WToVec(table), WToVec(indexes)));
}

HEADER_INLINE VecU16 vecu16_shuffle8(VecU16 table, VecU16 indexes) {
  return VecToU16(_mm_shuffle_epi8(U16ToVec(table), U16ToVec(indexes)));
}

HEADER_INLINE VecUc vecuc_shuffle8(VecUc table, VecUc indexes) {
  return VecToUc(_mm_shuffle_epi8(UcToVec(table), UcToVec(indexes)));
}
#    endif
#    ifdef USE_SSE42
HEADER_INLINE VecI32 veci32_max(VecI32 v1, VecI32 v2) {
  return VecToI32(_mm_max_epi32(I32ToVec(v1), I32ToVec(v2)));
}

HEADER_INLINE uintptr_t vecw_extract64_0(VecW vv) {
  return _mm_extract_epi64(WToVec(vv), 0);
}

HEADER_INLINE uintptr_t vecw_extract64_1(VecW vv) {
  return _mm_extract_epi64(WToVec(vv), 1);
}

HEADER_INLINE VecW vecw_blendv(VecW aa, VecW bb, VecW mask) {
  return VecToW(_mm_blendv_epi8(WToVec(aa), WToVec(bb), WToVec(mask)));
}

HEADER_INLINE VecU32 vecu32_blendv(VecU32 aa, VecU32 bb, VecU32 mask) {
  return VecToU32(_mm_blendv_epi8(U32ToVec(aa), U32ToVec(bb), U32ToVec(mask)));
}

HEADER_INLINE VecU16 vecu16_blendv(VecU16 aa, VecU16 bb, VecU16 mask) {
  return VecToU16(_mm_blendv_epi8(U16ToVec(aa), U16ToVec(bb), U16ToVec(mask)));
}

HEADER_INLINE VecUc vecuc_blendv(VecUc aa, VecUc bb, VecUc mask) {
  return VecToUc(_mm_blendv_epi8(UcToVec(aa), UcToVec(bb), UcToVec(mask)));
}
#    else // USE_SSE2, !USE_SSE42
HEADER_INLINE uintptr_t vecw_extract64_0(VecW vv) {
  return R_CAST(uintptr_t, _mm_movepi64_pi64(WToVec(vv)));
}

// compiler recognizes this on ARMv8
HEADER_INLINE uintptr_t vecw_extract64_1(VecW vv) {
  const __m128i v0 = _mm_srli_si128(WToVec(vv), 8);
  return R_CAST(uintptr_t, _mm_movepi64_pi64(v0));
}

// N.B. we do *not* enforce the low bits of each mask byte matching the high
// bit.
HEADER_INLINE VecW vecw_blendv(VecW aa, VecW bb, VecW mask) {
  return vecw_and_notfirst(mask, aa) | (mask & bb);
}

HEADER_INLINE VecU32 vecu32_blendv(VecU32 aa, VecU32 bb, VecU32 mask) {
  return vecu32_and_notfirst(mask, aa) | (mask & bb);
}

HEADER_INLINE VecU16 vecu16_blendv(VecU16 aa, VecU16 bb, VecU16 mask) {
  return vecu16_and_notfirst(mask, aa) | (mask & bb);
}

HEADER_INLINE VecUc vecuc_blendv(VecUc aa, VecUc bb, VecUc mask) {
  return vecuc_and_notfirst(mask, aa) | (mask & bb);
}
#    endif

HEADER_INLINE VecI16 veci16_max(VecI16 v1, VecI16 v2) {
  return VecToI16(_mm_max_epi16(I16ToVec(v1), I16ToVec(v2)));
}

HEADER_INLINE VecW vecw_sad(VecW v1, VecW v2) {
  return VecToW(_mm_sad_epu8(WToVec(v1), WToVec(v2)));
}

HEADER_INLINE VecUc vecuc_add(VecUc v1, VecUc v2) {
  return VecToUc(_mm_add_epi8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecUc vecuc_adds(VecUc v1, VecUc v2) {
  return VecToUc(_mm_adds_epu8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecUc vecuc_signed_cmpgt(VecUc v1, VecUc v2) {
  return VecToUc(_mm_cmpgt_epi8(UcToVec(v1), UcToVec(v2)));
}

HEADER_INLINE VecU16 vecu16_min8(VecU16 v1, VecU16 v2) {
  return VecToU16(_mm_min_epu8(U16ToVec(v1), U16ToVec(v2)));
}

HEADER_INLINE VecUc vecuc_min(VecUc v1, VecUc v2) {
  return VecToUc(_mm_min_epu8(UcToVec(v1), UcToVec(v2)));
}
#  endif  // USE_SSE2, !USE_AVX2

HEADER_INLINE void vecw_storeu(void* mem_addr, VecW vv) {
  vec_storeu(mem_addr, WToVec(vv));
}

HEADER_INLINE void vecu32_storeu(void* mem_addr, VecU32 vv) {
  vec_storeu(mem_addr, U32ToVec(vv));
}

HEADER_INLINE void veci32_storeu(void* mem_addr, VecI32 vv) {
  vec_storeu(mem_addr, I32ToVec(vv));
}

HEADER_INLINE void vecu16_storeu(void* mem_addr, VecU16 vv) {
  vec_storeu(mem_addr, U16ToVec(vv));
}

HEADER_INLINE void veci16_storeu(void* mem_addr, VecI16 vv) {
  vec_storeu(mem_addr, I16ToVec(vv));
}

HEADER_INLINE void vecuc_storeu(void* mem_addr, VecUc vv) {
  vec_storeu(mem_addr, UcToVec(vv));
}

HEADER_INLINE VecW vecw_bytesum(VecW src, VecW m0) {
  return vecw_sad(src, m0);
}

CONSTI32(kVec8thUintPerWord, sizeof(intptr_t) / sizeof(Vec8thUint));

#  ifdef FVEC_32

#    ifndef __FMA__
#      error "32-byte-float-vector builds require FMA3 as well."
#    endif

CONSTI32(kBytesPerFVec, 32);
CONSTI32(kBytesPerDVec, 32);
typedef float VecF __attribute__ ((vector_size (32)));
typedef double VecD __attribute__ ((vector_size (32)));

#    define VCONST_F(xx) {xx, xx, xx, xx, xx, xx, xx, xx}
#    define VCONST_D(xx) {xx, xx, xx, xx}

HEADER_INLINE VecF VecToF(__m256 xxv) {
  return R_CAST(VecF, xxv);
}

HEADER_INLINE VecD VecToD(__m256d xxv) {
  return R_CAST(VecD, xxv);
}

HEADER_INLINE __m256 FToVec(VecF xxv) {
  return R_CAST(__m256, xxv);
}

HEADER_INLINE __m256d DToVec(VecD xxv) {
  return R_CAST(__m256d, xxv);
}

HEADER_INLINE VecF vecf_setzero() {
  return VecToF(_mm256_setzero_ps());
}

HEADER_INLINE VecD vecd_setzero() {
  return VecToD(_mm256_setzero_pd());
}

#  else  // !FVEC_32

CONSTI32(kBytesPerFVec, 16);
CONSTI32(kBytesPerDVec, 16);
typedef float VecF __attribute__ ((vector_size (16)));
typedef double VecD __attribute__ ((vector_size (16)));

#    define VCONST_F(xx) {xx, xx, xx, xx}
#    define VCONST_D(xx) {xx, xx}

HEADER_INLINE VecF VecToF(__m128 xxv) {
  return R_CAST(VecF, xxv);
}

HEADER_INLINE VecD VecToD(__m128d xxv) {
  return R_CAST(VecD, xxv);
}

HEADER_INLINE __m128 FToVec(VecF xxv) {
  return R_CAST(__m128, xxv);
}

HEADER_INLINE __m128d DToVec(VecD xxv) {
  return R_CAST(__m128d, xxv);
}

HEADER_INLINE VecF vecf_setzero() {
  return VecToF(_mm_setzero_ps());
}

HEADER_INLINE VecD vecd_setzero() {
  return VecToD(_mm_setzero_pd());
}

#  endif  // !FVEC_32

HEADER_INLINE VecUc VecWToUc(VecW vv) {
  return R_CAST(VecUc, vv);
}

HEADER_INLINE VecW VecU16ToW(VecU16 vv) {
  return R_CAST(VecW, vv);
}

HEADER_INLINE VecW VecUcToW(VecUc vv) {
  return R_CAST(VecW, vv);
}

HEADER_INLINE void vecw_lo_and_hi_nybbles(VecW cur_vec, VecW m4, VecW* vec_lo_ptr, VecW* vec_hi_ptr) {
  // Assumes m4 is VCONST_W(kMask0F0F).
  // Returned vec_lo and vec_hi have top nybble of each byte zeroed out.
  cur_vec = vecw_permute0xd8_if_avx2(cur_vec);
  // AVX2:
  //   vec_even contains {0, 2, 4, ..., 14, 32, 34, ..., 46,
  //                      16, 18, ..., 30, 48, ... 62}
  //   vec_odd contains {1, 3, 5, ..., 15, 33, 35, ..., 47,
  //                     17, 19, ..., 31, 49, ..., 63}
  // SSE2:
  //   vec_even contains {0, 2, 4, ..., 30}
  //   vec_odd contains {1, 3, 5, ..., 31}
  const VecW vec_even = cur_vec & m4;
  const VecW vec_odd = vecw_srli(cur_vec, 4) & m4;

  // AVX2:
  //   vec_lo contains {0, 1, 2, ..., 31}
  //   vec_hi contains {32, 33, 34, ..., 63}
  // SSE2:
  //   vec_lo contains {0, 1, 2, ..., 15}
  //   vec_hi contains {16, 17, 18, ..., 31}
  *vec_lo_ptr = vecw_unpacklo8(vec_even, vec_odd);
  *vec_hi_ptr = vecw_unpackhi8(vec_even, vec_odd);
}
#else  // !USE_SSE2
#  ifdef __LP64__
CONSTI32(kBytesPerVec, 8);
#  else
CONSTI32(kBytesPerVec, 4);
#  endif
CONSTI32(kBytesPerFVec, 4);
CONSTI32(kBytesPerDVec, 8);

typedef uintptr_t VecW;
typedef uint32_t VecU32;
typedef float VecF;
typedef double VecD;
// VecI16 and VecI8 aren't worth the trouble of scaling down to 32-bit

#  define VCONST_W(xx) (xx)

HEADER_INLINE VecW vecw_setzero() {
  return k0LU;
}

HEADER_INLINE VecW vecw_srli(VecW vv, uint32_t ct) {
  return vv >> ct;
}

HEADER_INLINE VecW vecw_slli(VecW vv, uint32_t ct) {
  return vv << ct;
}

HEADER_INLINE VecW vecw_set1(uintptr_t ulii) {
  return ulii;
}

HEADER_INLINE VecW vecw_loadu(const void* mem_addr) {
  return *S_CAST(const VecW*, mem_addr);
}

#  ifdef __LP64__
HEADER_INLINE VecW vecw_bytesum(VecW src, __maybe_unused VecW m0) {
  src = (src & 0x00ff00ff00ff00ffLLU) + ((src >> 8) & 0x00ff00ff00ff00ffLLU);
  return (src * 0x1000100010001LLU) >> 48;
}
#  else
HEADER_INLINE VecW vecw_bytesum(VecW src, __maybe_unused VecW m0) {
  src = (src & 0x00ff00ff) + ((src >> 8) & 0x00ff00ff);
  return (src & 0xffff) + (src >> 16);
}
#  endif

HEADER_INLINE VecW vecw_and_notfirst(VecW excl, VecW main) {
  return (~excl) & main;
}

HEADER_INLINE VecU32 vecu32_and_notfirst(VecU32 excl, VecU32 main) {
  return (~excl) & main;
}
#endif  // !USE_SSE2

CONSTI32(kBitsPerVec, kBytesPerVec * CHAR_BIT);

// We now use Knuth's Nyp/Nybble vocabulary for 2-bit and 4-bit elements,
// respectively.
CONSTI32(kNypsPerVec, kBytesPerVec * 4);
CONSTI32(kNybblesPerVec, kBytesPerVec * 2);
CONSTI32(kWordsPerVec, kBytesPerVec / kBytesPerWord);
CONSTI32(kInt32PerVec, kBytesPerVec / 4);
CONSTI32(kInt16PerVec, kBytesPerVec / 2);

CONSTI32(kFloatPerFVec, kBytesPerFVec / 4);
CONSTI32(kDoublePerDVec, kBytesPerDVec / 8);

CONSTI32(kVecsPerCacheline, kCacheline / kBytesPerVec);

// debug
HEADER_INLINE void PrintVec(const void* vv) {
  const unsigned char* vv_alias = S_CAST(const unsigned char*, vv);
  for (uint32_t uii = 0; uii != kBytesPerVec; ++uii) {
    printf("%u ", vv_alias[uii]);
  }
  printf("\n");
}

HEADER_INLINE void PrintVecD(const VecD* vv_ptr, const char* preprint) {
  fputs(preprint, stdout);
  const double* vv_alias = R_CAST(const double*, vv_ptr);
  for (uint32_t uii = 0; uii != kDoublePerDVec; ++uii) {
    printf(" %g", vv_alias[uii]);
  }
  fputs("\n", stdout);
}

typedef union {
  VecW vw;

  STD_ARRAY_DECL(uintptr_t, kWordsPerVec, w);

  STD_ARRAY_DECL(uint32_t, kInt32PerVec, u32);
} UniVec;

typedef union {
  VecF vf;
  STD_ARRAY_DECL(float, kFloatPerFVec, f4);
} UniVecF;

typedef union {
  VecD vd;
  STD_ARRAY_DECL(double, kDoublePerDVec, d8);
} UniVecD;

// sum must fit in 16 bits
HEADER_INLINE uintptr_t UniVecHsum16(UniVec uv) {
#ifdef __LP64__
#  ifdef USE_AVX2
  return ((uv.w[0] + uv.w[1] + uv.w[2] + uv.w[3]) * kMask0001) >> 48;
#  else
  return ((uv.w[0] + uv.w[1]) * kMask0001) >> 48;
#  endif
#else
  return (uv.w[0] * kMask0001) >> 16;
#endif
}

// sum must fit in 32 bits
HEADER_INLINE uintptr_t UniVecHsum32(UniVec uv) {
#ifdef __LP64__
#  ifdef USE_AVX2
  return ((uv.w[0] + uv.w[1] + uv.w[2] + uv.w[3]) * kMask00000001) >> 32;
#  else
  return ((uv.w[0] + uv.w[1]) * kMask00000001) >> 32;
#  endif
#else
  return uv.w[0];
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

HEADER_INLINE double VecDHsum(VecD vecd) {
  UniVecD uvd;
  uvd.vd = vecd;
#ifdef __LP64__
#  ifdef FVEC_32
  return uvd.d8[0] + uvd.d8[1] + uvd.d8[2] + uvd.d8[3];
#  else
  return uvd.d8[0] + uvd.d8[1];
#  endif
#else
  return uvd.d8[0];
#endif
}

#ifdef USE_AVX2
HEADER_INLINE Vec4thUint UnpackVec8thUintTo4th(Vec8thUint hw) {
  return _pdep_u64(hw, kMask5555);
}

HEADER_INLINE Vec8thUint PackVec4thUintTo8th(Vec4thUint ww) {
  return _pext_u64(ww, kMask5555);
}

HEADER_INLINE Vec16thUint PackVec8thUintTo16th(Vec8thUint ww) {
  return _pext_u64(ww, kMask5555);
}

// See https://stackoverflow.com/questions/21622212/how-to-perform-the-inverse-of-mm256-movemask-epi8-vpmovmskb .
HEADER_INLINE VecUc InverseMovemaskFF(Vec8thUint mask) {
  __m256i vmask = _mm256_set1_epi32(mask);
  const __m256i byte_gather = _mm256_setr_epi64x(0, kMask0101, 2 * kMask0101, 3 * kMask0101);
  vmask = _mm256_shuffle_epi8(vmask, byte_gather);
  const __m256i bit_mask = _mm256_set1_epi64x(0x7fbfdfeff7fbfdfeLL);
  vmask = _mm256_or_si256(vmask, bit_mask);
  return R_CAST(VecUc, _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1)));
}

// If we're only interested in the even bits of mask.  No need to mask out odd
// bits before calling.
HEADER_INLINE VecUc InverseMovespreadmaskFF(Vec4thUint mask) {
  __m256i vmask = _mm256_set1_epi64x(mask);
  const __m256i byte_gather = _mm256_setr_epi32(0, 0x01010101, 0x02020202, 0x03030303, 0x04040404, 0x05050505, 0x06060606, 0x07070707);
  vmask = _mm256_shuffle_epi8(vmask, byte_gather);
  const __m256i bit_mask = _mm256_set1_epi32(0xbfeffbfeU);
  vmask = _mm256_or_si256(vmask, bit_mask);
  return R_CAST(VecUc, _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1)));
}
#else  // !USE_AVX2
#  ifdef USE_SSE2
HEADER_INLINE Vec4thUint UnpackVec8thUintTo4th(Vec8thUint hw) {
  hw = (hw | (hw << 8)) & 0x00ff00ffU;
  hw = (hw | (hw << 4)) & 0x0f0f0f0fU;
  hw = (hw | (hw << 2)) & 0x33333333U;
  return (hw | (hw << 1)) & 0x55555555U;
}

HEADER_INLINE Vec8thUint PackVec4thUintTo8th(Vec4thUint ww) {
  ww = (ww | (ww >> 1)) & kMask3333;
  ww = (ww | (ww >> 2)) & kMask0F0F;
  ww = (ww | (ww >> 4)) & kMask00FF;
  return S_CAST(Vec8thUint, ww | (ww >> 8));
}

HEADER_INLINE Vec16thUint PackVec8thUintTo16th(Vec8thUint ww) {
  ww = (ww | (ww >> 1)) & 0x3333;
  ww = (ww | (ww >> 2)) & 0x0f0f;
  return S_CAST(Vec16thUint, ww | (ww >> 4));
}

#    ifdef USE_SSE42
HEADER_INLINE VecUc InverseMovemaskFF(Vec8thUint mask) {
  __m128i vmask = _mm_set1_epi16(mask);
  const __m128i byte_gather = _mm_setr_epi32(0, 0, 0x01010101, 0x01010101);
  vmask = _mm_shuffle_epi8(vmask, byte_gather);
  const __m128i bit_mask = _mm_set1_epi64x(0x7fbfdfeff7fbfdfeLL);
  vmask = _mm_or_si128(vmask, bit_mask);
  return R_CAST(VecUc, _mm_cmpeq_epi8(vmask, _mm_set1_epi64x(-1)));
}

HEADER_INLINE VecUc InverseMovespreadmaskFF(Vec4thUint mask) {
  __m128i vmask = _mm_set1_epi32(mask);
  const __m128i byte_gather = _mm_setr_epi32(0, 0x01010101, 0x02020202, 0x03030303);
  vmask = _mm_shuffle_epi8(vmask, byte_gather);
  const __m128i bit_mask = _mm_set1_epi32(0xbfeffbfeU);
  vmask = _mm_or_si128(vmask, bit_mask);
  return R_CAST(VecUc, _mm_cmpeq_epi8(vmask, _mm_set1_epi64x(-1)));
}
#    endif

#  endif // !USE_SSE2
#endif  // !USE_AVX2

HEADER_CINLINE uintptr_t BitCtToVecCt(uintptr_t val) {
  return DivUp(val, kBitsPerVec);
}

HEADER_CINLINE uintptr_t BitCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * BitCtToVecCt(val);
}

HEADER_CINLINE uintptr_t Int32CtToVecCt(uintptr_t val) {
  return DivUp(val, kInt32PerVec);
}

HEADER_CINLINE uintptr_t WordCtToVecCt(uintptr_t val) {
  return DivUp(val, kWordsPerVec);
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

HEADER_CINLINE uintptr_t DblCtToVecCt(uintptr_t val) {
  return Int64CtToVecCt(val);
}

HEADER_CINLINE uintptr_t VecCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kVecsPerCacheline);
}

HEADER_CINLINE uint64_t VecCtToCachelineCtU64(uint64_t val) {
  return DivUpU64(val, kVecsPerCacheline);
}

#ifdef USE_SSE2
#  ifdef USE_SSE42
HEADER_INLINE uint32_t PopcountVec8thUint(uint32_t val) {
  return __builtin_popcount(val);
}
#  else
HEADER_INLINE uint32_t PopcountVec8thUint(uint32_t val) {
  // May as well exploit the fact that only the low 16 bits may be set.
  val = val - ((val >> 1) & 0x5555);
  val = (val & 0x3333) + ((val >> 2) & 0x3333);
  val = (val + (val >> 4)) & 0x0f0f;
  return (val + (val >> 8)) & 0xff;
}
#  endif
#endif

#ifdef USE_SSE2
// On ARM, emulated movemask isn't great.  But there are alternative
// instructions that efficiently perform what you usually want to do with
// movemasks:
//   https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon

#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
// vec0255 refers to a VecUc where all bytes are equal to 0 or 255.
// (possible todo: define this as its own type, with automatic downcast to
// VecUc but not the other way around.)
//
// Return value in nonzero case is architecture-dependent (though always
// nonzero).  So, best to only use this in other architecture-specific code;
// hence the leading underscore in the function name.
HEADER_INLINE uint64_t _vec0255_is_nonzero(VecUc vv) {
  return vecuc_movemask(vv);
}

HEADER_INLINE uint32_t vec0255_set_ct(VecUc vv) {
  return PopcountVec8thUint(vecuc_movemask(vv));
}

HEADER_INLINE uint32_t vec0255_is_all_set(VecUc vv) {
  return (vecuc_movemask(vv) == kVec8thUintMax);
}

HEADER_INLINE uint32_t vec0255u16_is_all_set(VecU16 vv) {
  return (vecu16_movemask(vv) == kVec8thUintMax);
}

HEADER_INLINE uint32_t m128is_are_equal(__m128i v1, __m128i v2) {
  return (_mm_movemask_epi8(_mm_cmpeq_epi8(v1, v2)) == 65535);
}
#  else
HEADER_INLINE uint64_t arm_shrn4_uc(VecUc vv) {
  uint16x8_t vv_;
  memcpy(&vv_, &vv, sizeof(vv_));
  return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vv_, 4)), 0);
}

HEADER_INLINE uint64_t arm_shrn4_i8(VecI8 vv) {
  uint16x8_t vv_;
  memcpy(&vv_, &vv, sizeof(vv_));
  return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vv_, 4)), 0);
}

HEADER_INLINE uint64_t arm_shrn4_u16(VecU16 vv) {
  uint16x8_t vv_;
  memcpy(&vv_, &vv, sizeof(vv_));
  return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vv_, 4)), 0);
}

HEADER_INLINE uint64_t arm_shrn4_m128i(__m128i vv) {
  uint16x8_t vv_;
  memcpy(&vv_, &vv, sizeof(vv_));
  return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vv_, 4)), 0);
}

HEADER_INLINE uint64_t _vec0255_is_nonzero(VecUc vv) {
  return arm_shrn4_uc(vv);
}

// set_bits4 must only have bits in 0, 4, ..., 60 set.
// move this out of ARM-only ifdef if we ever want this elsewhere.
HEADER_INLINE uint32_t popcount_bits4(uint64_t set_bits4) {
  // Branchlessly count the number of set bits, taking advantage of the limited
  // set of positions they can be in.
  // Multiplication by the magic constant kMask1111 usually puts the sum of
  // all 16 bits of interest in the high nybble of the result... except that
  // the nybble overflows when all 16 bits are set.  We work around this by
  // (i) multiplying by (kMask1111 >> 4) instead, which excludes the lowest bit
  //     from the high-nybble sum, and
  // (ii) then adding the lowest bit afterward.
  const uint32_t set_ct_excluding_lowest = (set_bits4 * (kMask1111 >> 4)) >> 60;
  return set_ct_excluding_lowest + (set_bits4 & 1);
}

// xx must have all nybbles equal to 0 or 15.
HEADER_INLINE uint32_t count_set_nybbles(uint64_t xx) {
  return popcount_bits4(xx & kMask1111);
}

HEADER_INLINE uint32_t vec0255_set_ct(VecUc vv) {
  return count_set_nybbles(arm_shrn4_uc(vv));
}

HEADER_INLINE uint32_t vec0255_is_all_set(VecUc vv) {
  return (arm_shrn4_uc(vv) == UINT64_MAX);
}

HEADER_INLINE uint32_t vec0255u16_is_all_set(VecU16 vv) {
  return (arm_shrn4_u16(vv) == UINT64_MAX);
}

HEADER_INLINE uint32_t m128is_are_equal(__m128i v1, __m128i v2) {
  const uint64_t set_nybbles = arm_shrn4_m128i(_mm_cmpeq_epi8(v1, v2));
  return (set_nybbles == UINT64_MAX);
}
#  endif

HEADER_INLINE uint32_t vecucs_are_equal(VecUc v1, VecUc v2) {
  return vec0255_is_all_set(v1 == v2);
}

HEADER_INLINE uint32_t vecu16s_are_equal(VecU16 v1, VecU16 v2) {
  return vec0255u16_is_all_set(v1 == v2);
}
#endif

// Downcasts don't risk alignment issues.
HEADER_INLINE uintptr_t* DowncastVecWToW(VecW* pp) {
  return R_CAST(uintptr_t*, pp);
}

HEADER_INLINE uint32_t* DowncastVecWToU32(VecW* pp) {
  return R_CAST(uint32_t*, pp);
}

#ifdef USE_SSE2
HEADER_INLINE Vec8thUint* DowncastWToV8(uintptr_t* pp) {
  return R_CAST(Vec8thUint*, pp);
}
#endif

HEADER_INLINE const uintptr_t* DowncastKVecWToW(const VecW* pp) {
  return R_CAST(const uintptr_t*, pp);
}

HEADER_INLINE const uint16_t* DowncastKVecWToU16(const VecW* pp) {
  return R_CAST(const uint16_t*, pp);
}


HEADER_INLINE uint32_t IsVecAligned(const void* ptr) {
  return !(R_CAST(uintptr_t, ptr) % kBytesPerVec);
}

#ifdef USE_SSE2
HEADER_INLINE void AlignWToVec(uintptr_t** pp) {
  const uintptr_t addr = R_CAST(uintptr_t, *pp);
  *pp = R_CAST(uintptr_t*, RoundUpPow2(addr, kBytesPerVec));
}
#else
HEADER_INLINE void AlignWToVec(__maybe_unused uintptr_t** pp) {
}
#endif

HEADER_INLINE void AlignKUcToVec(const unsigned char** pp) {
  const uintptr_t addr = R_CAST(uintptr_t, *pp);
  *pp = R_CAST(const unsigned char*, RoundUpPow2(addr, kBytesPerVec));
}

/*
HEADER_INLINE uint32_t AlignToVecW(void* prestart, VecW** result_ptr) {
  unsigned char* prestart_uc = S_CAST(unsigned char*, prestart);
  const uint32_t lead_byte_ct = (-R_CAST(uintptr_t, prestart_uc)) % kBytesPerVec;
  *result_ptr = R_CAST(VecW*, &(prestart_uc[lead_byte_ct]));
  return lead_byte_ct;
}
*/

HEADER_INLINE uint32_t AlignKToAW(const void* prestart, const uintptr_t** result_ptr) {
  const unsigned char* prestart_uc = S_CAST(const unsigned char*, prestart);
  const uint32_t lead_byte_ct = (-R_CAST(uintptr_t, prestart_uc)) % kBytesPerVec;
  *result_ptr = R_CAST(const uintptr_t*, &(prestart_uc[lead_byte_ct]));
  return lead_byte_ct;
}


HEADER_INLINE BoolErr vecaligned_malloc(uintptr_t size, void* aligned_pp) {
#ifdef USE_AVX2
  return aligned_malloc(size, kBytesPerVec, aligned_pp);
#else
#  if defined(__APPLE__) || !defined(__LP64__)
  const BoolErr ret_boolerr = pgl_malloc(size, S_CAST(uintptr_t*, aligned_pp));
  assert(IsVecAligned(*S_CAST(uintptr_t**, aligned_pp)));
  return ret_boolerr;
#  else
  return aligned_malloc(size, kBytesPerVec, aligned_pp);
#  endif
#endif
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


#if defined(USE_SSE2) && !defined(NO_UNALIGNED)
int32_t memequal(const void* m1, const void* m2, uintptr_t byte_ct);

// This is also better than the June 2018 OS X/LLVM stock implementation,
// especially for small values of ct.
// (gcc 7.1 and clang 6.0.0 should have better stock implementations;
// re-benchmark this once Linux build machine is upgraded to Ubuntu 18.04.)
int32_t Memcmp(const void* m1, const void* m2, uintptr_t ct);
#else
HEADER_INLINE int32_t memequal(const void* m1, const void* m2, uintptr_t byte_ct) {
  return !memcmp(m1, m2, byte_ct);
}

HEADER_INLINE int32_t Memcmp(const void* m1, const void* m2, uintptr_t ct) {
  return memcmp(m1, m2, ct);
}
#endif


#if defined(USE_SSE2) && defined(__cplusplus) && !defined(NO_UNALIGNED)
// See https://stackoverflow.com/questions/9510514/integer-range-based-template-specialisation .

template <bool> struct TRange;

// This makes MemequalKImpl<byte_ct> expand to
// MemequalKImpl<byte_ct, TRange<true> >.
// If a later single-parameter template defines the same thing, that takes
// precedence.
template <uint32_t N, typename = TRange<true> > struct MemequalKImpl {
  static int32_t MemequalK(const void* m1, const void* m2) {
    return memequal(m1, m2, N);
  }
};

template <> struct MemequalKImpl<1> {
  static int32_t MemequalK(const void* m1, const void* m2) {
    const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
    const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
    return (m1_uc[0] == m2_uc[0]);
  }
};

template <> struct MemequalKImpl<2> {
  static int32_t MemequalK(const void* m1, const void* m2) {
    return ((*R_CAST(const uint16_t*, m1)) == (*R_CAST(const uint16_t*, m2)));
  }
};

template <> struct MemequalKImpl<3> {
  static int32_t MemequalK(const void* m1, const void* m2) {
    const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
    const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
    return
      ((*R_CAST(const uint16_t*, m1)) == (*R_CAST(const uint16_t*, m2))) &&
      (m1_uc[2] == m2_uc[2]);
  }
};

template <> struct MemequalKImpl<4> {
  static int32_t MemequalK(const void* m1, const void* m2) {
    return ((*R_CAST(const uint32_t*, m1)) == (*R_CAST(const uint32_t*, m2)));
  }
};

template <uint32_t N> struct MemequalKImpl<N, TRange<(5 <= N) && (N <= 7)> > {
  static int32_t MemequalK(const void* m1, const void* m2) {
    const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
    const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
    return
      ((*R_CAST(const uint32_t*, m1)) == (*R_CAST(const uint32_t*, m2))) &&
      ((*R_CAST(const uint32_t*, &(m1_uc[N - 4]))) == (*R_CAST(const uint32_t*, &(m2_uc[N - 4]))));
  }
};

template <> struct MemequalKImpl<8> {
  static int32_t MemequalK(const void* m1, const void* m2) {
    return ((*R_CAST(const uint64_t*, m1)) == (*R_CAST(const uint64_t*, m2)));
  }
};

template <uint32_t N> struct MemequalKImpl<N, TRange<(9 <= N) && (N <= 15)> > {
  static int32_t MemequalK(const void* m1, const void* m2) {
    const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
    const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
    return
      ((*R_CAST(const uint64_t*, m1)) == (*R_CAST(const uint64_t*, m2))) &&
      ((*R_CAST(const uint64_t*, &(m1_uc[N - 8]))) == (*R_CAST(const uint64_t*, &(m2_uc[N - 8]))));
  }
};

template <> struct MemequalKImpl<16> {
  static int32_t MemequalK(const void* m1, const void* m2) {
    const __m128i v1 = _mm_loadu_si128(S_CAST(const __m128i*, m1));
    const __m128i v2 = _mm_loadu_si128(S_CAST(const __m128i*, m2));
    return m128is_are_equal(v1, v2);
  }
};

template <uint32_t N> struct MemequalKImpl<N, TRange<(17 <= N) && (N <= 24)> > {
  static int32_t MemequalK(const void* m1, const void* m2) {
    const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
    const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
    const __m128i v1 = _mm_loadu_si128(S_CAST(const __m128i*, m1));
    const __m128i v2 = _mm_loadu_si128(S_CAST(const __m128i*, m2));
    return
      m128is_are_equal(v1, v2) &&
      ((*R_CAST(const uint64_t*, &(m1_uc[N - 8]))) == (*R_CAST(const uint64_t*, &(m2_uc[N - 8]))));
  }
};

template <uint32_t N> struct MemequalKImpl<N, TRange<(25 <= N) && (N <= 31)> > {
  static int32_t MemequalK(const void* m1, const void* m2) {
    __m128i v1 = _mm_loadu_si128(S_CAST(const __m128i*, m1));
    __m128i v2 = _mm_loadu_si128(S_CAST(const __m128i*, m2));
    if (!m128is_are_equal(v1, v2)) {
      return 0;
    }
    const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
    const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
    v1 = _mm_loadu_si128(R_CAST(const __m128i*, &(m1_uc[N - 16])));
    v2 = _mm_loadu_si128(R_CAST(const __m128i*, &(m2_uc[N - 16])));
    return m128is_are_equal(v1, v2);
  }
};

#  define memequal_k(m1, m2, byte_ct) plink2::MemequalKImpl<byte_ct>::MemequalK(m1, m2)

template <uint32_t N, typename = TRange<true> > struct MemcpyKImpl {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    memcpy(dst, src, N);
  }
};

// Patch a bunch of cases where some commonly-used gcc and clang versions
// generate suboptimal code.  (Since this code is shamelessly x86-specific, we
// don't worry about the formal undefinedness of unaligned pointer dereferences
// here.)
// (todo: check if/when this has been fixed, and remove this bloat once all
// production build machines have sufficiently new compilers.)
template <> struct MemcpyKImpl<2> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    *S_CAST(uint16_t*, dst) = *S_CAST(const uint16_t*, src);
  }
};

template <> struct MemcpyKImpl<3> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    *S_CAST(uint16_t*, dst) = *S_CAST(const uint16_t*, src);
    dst_uc[2] = src_uc[2];
  }
};

template <> struct MemcpyKImpl<5> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    *S_CAST(uint32_t*, dst) = *S_CAST(const uint32_t*, src);
    dst_uc[4] = src_uc[4];
  }
};

template <> struct MemcpyKImpl<6> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    uint16_t* dst_u16 = S_CAST(uint16_t*, dst);
    const uint16_t* src_u16 = S_CAST(const uint16_t*, src);
    *S_CAST(uint32_t*, dst) = *S_CAST(const uint32_t*, src);
    dst_u16[2] = src_u16[2];
  }
};

template <> struct MemcpyKImpl<7> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    *S_CAST(uint32_t*, dst) = *S_CAST(const uint32_t*, src);
    *R_CAST(uint32_t*, &(dst_uc[3])) = *R_CAST(const uint32_t*, &(src_uc[3]));
  }
};

template <> struct MemcpyKImpl<9> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    *S_CAST(uint64_t*, dst) = *S_CAST(const uint64_t*, src);
    dst_uc[8] = src_uc[8];
  }
};

template <> struct MemcpyKImpl<10> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    uint16_t* dst_u16 = S_CAST(uint16_t*, dst);
    const uint16_t* src_u16 = S_CAST(const uint16_t*, src);
    *S_CAST(uint64_t*, dst) = *S_CAST(const uint64_t*, src);
    dst_u16[4] = src_u16[4];
  }
};

template <uint32_t N> struct MemcpyKImpl<N, TRange<(11 <= N) && (N <= 12)> > {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    *S_CAST(uint64_t*, dst) = *S_CAST(const uint64_t*, src);
    *R_CAST(uint32_t*, &(dst_uc[N - 4])) = *R_CAST(const uint32_t*, &(src_uc[N - 4]));
  }
};

template <uint32_t N> struct MemcpyKImpl<N, TRange<(13 <= N) && (N <= 15)> > {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    *S_CAST(uint64_t*, dst) = *S_CAST(const uint64_t*, src);
    *R_CAST(uint64_t*, &(dst_uc[N - 8])) = *R_CAST(const uint64_t*, &(src_uc[N - 8]));
  }
};

template <> struct MemcpyKImpl<17> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    const __m128i vv = _mm_loadu_si128(S_CAST(const __m128i*, src));
    _mm_storeu_si128(S_CAST(__m128i*, dst), vv);
    dst_uc[16] = src_uc[16];
  }
};

template <> struct MemcpyKImpl<18> {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    uint16_t* dst_u16 = S_CAST(uint16_t*, dst);
    const uint16_t* src_u16 = S_CAST(const uint16_t*, src);
    const __m128i vv = _mm_loadu_si128(S_CAST(const __m128i*, src));
    _mm_storeu_si128(S_CAST(__m128i*, dst), vv);
    dst_u16[8] = src_u16[8];
  }
};

template <uint32_t N> struct MemcpyKImpl<N, TRange<(19 <= N) && (N <= 20)> > {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    const __m128i vv = _mm_loadu_si128(S_CAST(const __m128i*, src));
    _mm_storeu_si128(S_CAST(__m128i*, dst), vv);
    *R_CAST(uint32_t*, &(dst_uc[N - 4])) = *R_CAST(const uint32_t*, &(src_uc[N - 4]));
  }
};

template <uint32_t N> struct MemcpyKImpl<N, TRange<(21 <= N) && (N <= 24)> > {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    const __m128i vv = _mm_loadu_si128(S_CAST(const __m128i*, src));
    _mm_storeu_si128(S_CAST(__m128i*, dst), vv);
    *R_CAST(uint64_t*, &(dst_uc[N - 8])) = *R_CAST(const uint64_t*, &(src_uc[N - 8]));
  }
};

template <uint32_t N> struct MemcpyKImpl<N, TRange<(25 <= N) && (N <= 31)> > {
  static void MemcpyK(void* __restrict dst, const void* __restrict src) {
    unsigned char* dst_uc = S_CAST(unsigned char*, dst);
    const unsigned char* src_uc = S_CAST(const unsigned char*, src);
    const __m128i v1 = _mm_loadu_si128(S_CAST(const __m128i*, src));
    const __m128i v2 = _mm_loadu_si128(R_CAST(const __m128i*, &(src_uc[N - 16])));
    _mm_storeu_si128(S_CAST(__m128i*, dst), v1);
    _mm_storeu_si128(R_CAST(__m128i*, &(dst_uc[N - 16])), v2);
  }
};

// Note that there's no difference between memcpy() and memcpy_k() for common
// 'well-behaved' sizes like 1, 4, 8, and 16.  It's the funny numbers in
// between, which often arise with constant strings, which this template is
// targeting.
#  define memcpy_k(dst, src, ct) plink2::MemcpyKImpl<ct>::MemcpyK(dst, src)

template <uint32_t N> char* MemcpyaK(void* __restrict dst, const void* __restrict src) {
  MemcpyKImpl<N>::MemcpyK(dst, src);
  char* dst_c = S_CAST(char*, dst);
  return &(dst_c[N]);
}

#  define memcpya_k(dst, src, ct) plink2::MemcpyaK<ct>(dst, src)
#  define memcpyua_k(dst, src, ct) CToUc(plink2::MemcpyaK<ct>(dst, src))

template <uint32_t N> struct MemcpyoKImpl {
  static void MemcpyoK(void* __restrict dst, const void* __restrict src) {
    MemcpyKImpl<N>::MemcpyK(dst, src);
  }
};

template <> struct MemcpyoKImpl<3> {
  static void MemcpyoK(void* __restrict dst, const void* __restrict src) {
    *S_CAST(uint32_t*, dst) = *S_CAST(const uint32_t*, src);
  }
};

template <> struct MemcpyoKImpl<7> {
  static void MemcpyoK(void* __restrict dst, const void* __restrict src) {
    *S_CAST(uint64_t*, dst) = *S_CAST(const uint64_t*, src);
  }
};

template <> struct MemcpyoKImpl<15> {
  static void MemcpyoK(void* __restrict dst, const void* __restrict src) {
    const __m128i vv = _mm_loadu_si128(S_CAST(const __m128i*, src));
    _mm_storeu_si128(S_CAST(__m128i*, dst), vv);
  }
};

// interestingly, __m256i copy does not seem to be better in 31 byte case

#  define memcpyo_k(dst, src, ct) plink2::MemcpyoKImpl<ct>::MemcpyoK(dst, src)

template <uint32_t N> char* MemcpyaoK(void* __restrict dst, const void* __restrict src) {
  MemcpyoKImpl<N>::MemcpyoK(dst, src);
  char* dst_c = S_CAST(char*, dst);
  return &(dst_c[N]);
}

#  define memcpyao_k(dst, src, ct) plink2::MemcpyaoK<ct>(dst, src)
#  define memcpyuao_k(dst, src, ct) DowncastToUc(plink2::MemcpyaoK<ct>(dst, src))

#  else  // !(defined(__LP64__) && defined(__cplusplus) && !defined(NO_UNALIGNED))

HEADER_INLINE int32_t memequal_k(const void* m1, const void* m2, uintptr_t ct) {
  return !memcmp(m1, m2, ct);
}

HEADER_INLINE void memcpy_k(void* __restrict dst, const void* __restrict src, uintptr_t ct) {
  memcpy(dst, src, ct);
}

HEADER_INLINE char* memcpya_k(void* __restrict dst, const void* __restrict src, uintptr_t ct) {
  return memcpya(dst, src, ct);
}

HEADER_INLINE unsigned char* memcpyua_k(void* __restrict dst, const void* __restrict src, uintptr_t ct) {
  return memcpyua(dst, src, ct);
}

HEADER_INLINE void memcpyo_k(void* __restrict dst, const void* __restrict src, uintptr_t ct) {
  memcpy(dst, src, ct);
}

HEADER_INLINE char* memcpyao_k(void* __restrict dst, const void* __restrict src, uintptr_t ct) {
  return memcpya(dst, src, ct);
}

HEADER_INLINE unsigned char* memcpyuao_k(void* __restrict dst, const void* __restrict src, uintptr_t ct) {
  return memcpyua(dst, src, ct);
}

#endif

#if defined(USE_SSE2) && (__cplusplus >= 201103L) && !defined(NO_UNALIGNED)

#  define strcpy_k(dst, src) plink2::MemcpyKImpl<plink2::CompileTimeSlen(src) + 1>::MemcpyK(dst, src);

#  define strcpya_k(dst, src) plink2::MemcpyaoK<plink2::CompileTimeSlen(src)>(dst, src);

#else

HEADER_INLINE void strcpy_k(char* __restrict dst, const void* __restrict src) {
  strcpy(dst, S_CAST(const char*, src));
}

HEADER_INLINE char* strcpya_k(char* __restrict dst, const void* __restrict src) {
  return strcpya(dst, src);
}

#endif

#ifdef USE_SSE2
HEADER_INLINE void CopyToUnalignedOffsetV8(unsigned char* dst, const Vec8thUint* src, uintptr_t offset) {
  memcpy_k(&(dst[offset * sizeof(Vec8thUint)]), src, sizeof(Vec8thUint));
}

HEADER_INLINE void CopyToUnalignedOffsetV16(unsigned char* dst, const Vec16thUint* src, uintptr_t offset) {
  memcpy_k(&(dst[offset * sizeof(Vec16thUint)]), src, sizeof(Vec16thUint));
}
#endif


// analogous to memset()
// this can be slightly slower if e.g. system supports AVX2 but non-AVX2 plink2
// build is in use; fine to pay that price given the small-array advantage for
// now.  Should revisit this after next build-machine Ubuntu upgrade, though.
HEADER_INLINE void vecset(void* target_vec, uintptr_t ww, uintptr_t vec_ct) {
  VecW* target_vec_iter = S_CAST(VecW*, target_vec);
#ifdef USE_SSE2
  const VecW payload = VCONST_W(ww);
  for (uintptr_t vec_idx = 0; vec_idx != vec_ct; ++vec_idx) {
    *target_vec_iter++ = payload;
  }
#else
  for (uintptr_t vec_idx = 0; vec_idx != vec_ct; ++vec_idx) {
    *target_vec_iter++ = ww;
  }
#endif
}

#if defined(USE_SSE2) && !defined(NO_UNALIGNED)
// This requires nbytes >= 4.
uintptr_t FirstUnequal4(const void* arr1, const void* arr2, uintptr_t nbytes);

HEADER_INLINE uintptr_t FirstUnequal(const void* arr1, const void* arr2, uintptr_t nbytes) {
  // Returns position of first byte mismatch, or nbytes if none was found.
  if (nbytes >= 4) {
    return FirstUnequal4(arr1, arr2, nbytes);
  }
  const char* s1 = S_CAST(const char*, arr1);
  const char* s2 = S_CAST(const char*, arr2);
  for (uintptr_t pos = 0; pos != nbytes; ++pos) {
    if (s1[pos] != s2[pos]) {
      return pos;
    }
  }
  return nbytes;
}
#else // !(defined(USE_SSE2) && !defined(NO_UNALIGNED))
// This requires nbytes >= kBytesPerWord.
uintptr_t FirstUnequalW(const void* arr1, const void* arr2, uintptr_t nbytes);

HEADER_INLINE uintptr_t FirstUnequal(const void* arr1, const void* arr2, uintptr_t nbytes) {
  // Returns position of first byte mismatch, or nbytes if none was found.
  if (nbytes >= kBytesPerWord) {
    return FirstUnequalW(arr1, arr2, nbytes);
  }
  const char* s1 = S_CAST(const char*, arr1);
  const char* s2 = S_CAST(const char*, arr2);
  for (uintptr_t pos = 0; pos != nbytes; ++pos) {
    if (s1[pos] != s2[pos]) {
      return pos;
    }
  }
  return nbytes;
}
#endif


HEADER_INLINE uintptr_t FirstUnequalFrom(const void* arr1, const void* arr2, uintptr_t start, uintptr_t nbytes) {
  const char* s1 = S_CAST(const char*, arr1);
  const char* s2 = S_CAST(const char*, arr2);
  return start + FirstUnequal(&(s1[start]), &(s2[start]), nbytes - start);
}


#ifdef __LP64__
uintptr_t CountVintsNonempty(const unsigned char* buf, const unsigned char* buf_end);

HEADER_INLINE uintptr_t CountVints(const unsigned char* buf, const unsigned char* buf_end) {
  if (buf == buf_end) {
    return 0;
  }
  return CountVintsNonempty(buf, buf_end);
}
#else
uintptr_t CountVints(const unsigned char* buf, const unsigned char* buf_end);

HEADER_INLINE uintptr_t CountVintsNonempty(const unsigned char* buf, const unsigned char* buf_end) {
  return CountVints(buf, buf_end);
}
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_SIMD_H__
