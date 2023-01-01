#ifndef __PLINK2_FMATH_H__
#define __PLINK2_FMATH_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2023 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Most of the functions are C ports of Shigeo Mitsunari's fast math library
// functions posted at https://github.com/herumi/fmath .  The original license
// is https://opensource.org/licenses/BSD-3-Clause .
// I've added comments to expd() which reflect my understanding of the
// algorithm; these are BSD-licensed.

#include "plink2_base.h"

#include <math.h>

#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef __LP64__
#  ifdef FVEC_32
extern const uint32_t kFloatExpLookupInt[] __attribute__((aligned(32)));
#  else
extern const uint32_t kFloatExpLookupInt[] __attribute__((aligned(16)));
#  endif

#  ifdef FVEC_32
HEADER_INLINE VecF fmath_exp_ps(VecF xxv) {
  __m256 xx = R_CAST(__m256, xxv);
  const __m256i mask7ff = _mm256_set1_epi32(0x7fffffff);
  // 88
  const __m256i max_x = _mm256_set1_epi32(0x42b00000);
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m256i min_x = _mm256_set1_epi64x(-0x3d4fffff3d500000LL);
  // 2^10 / log(2)
  const __m256i const_aa = _mm256_set1_epi32(0x44b8aa3b);
  // log(2) / 2^10
  const __m256i const_bb = _mm256_set1_epi32(0x3a317218);
  const __m256i f1 = _mm256_set1_epi32(0x3f800000);
  const __m256i mask_s = _mm256_set1_epi32(0x000003ff);
  const __m256i i127s = _mm256_set1_epi32(0x0001fc00);
  const __m256i limit = _mm256_castps_si256(_mm256_and_ps(xx, R_CAST(__m256, mask7ff)));
  const int32_t over = _mm256_movemask_epi8(_mm256_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm256_min_ps(xx, R_CAST(__m256, max_x));
    xx = _mm256_max_ps(xx, R_CAST(__m256, min_x));
  }
  const __m256i rr = _mm256_cvtps_epi32(_mm256_mul_ps(xx, R_CAST(__m256, const_aa)));
  __m256 tt = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(rr), R_CAST(__m256, const_bb), xx);
  tt = _mm256_add_ps(R_CAST(__m256, f1), tt);
  const __m256i v8 = _mm256_and_si256(rr, mask_s);
  __m256i u8 = _mm256_add_epi32(rr, i127s);
  u8 = _mm256_srli_epi32(u8, 10);
  u8 = _mm256_slli_epi32(u8, 23);
  __m256i ti = _mm256_i32gather_epi32(R_CAST(const int*, kFloatExpLookupInt), v8, 4);
  __m256 t0 = _mm256_castsi256_ps(ti);
  t0 = _mm256_or_ps(t0, _mm256_castsi256_ps(u8));
  return R_CAST(VecF, _mm256_mul_ps(tt, t0));
}
#  else  // !FVEC_32
static inline VecF fmath_exp_ps(VecF xxv) {
  const float* const kFloatExpLookup = R_CAST(const float*, kFloatExpLookupInt);
  __m128 xx = xxv;
  const __m128i mask7ff = _mm_set1_epi32(0x7fffffff);

  // 88
  const __m128i max_x = _mm_set1_epi32(0x42b00000);
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m128i min_x = _mm_set1_epi64x(-0x3d4fffff3d500000LL);
  // 2^10 / log(2)
  const __m128i const_aa = _mm_set1_epi32(0x44b8aa3b);
  // log(2) / 2^10
  const __m128i const_bb = _mm_set1_epi32(0x3a317218);

  const __m128i f1 = _mm_set1_epi32(0x3f800000);
  const __m128i mask_s = _mm_set1_epi64x(0x3ff000003ffLLU);
  const __m128i i127s = _mm_set1_epi32(0x0001fc00);
  const __m128i limit = _mm_castps_si128(_mm_and_ps(xx, R_CAST(__m128, mask7ff)));
  const int32_t over = _mm_movemask_epi8(_mm_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm_min_ps(xx, R_CAST(__m128, max_x));
    xx = _mm_max_ps(xx, R_CAST(__m128, min_x));
  }
  const __m128i rr = _mm_cvtps_epi32(_mm_mul_ps(xx, R_CAST(__m128, const_aa)));
  __m128 tt = _mm_sub_ps(xx, _mm_mul_ps(_mm_cvtepi32_ps(rr), R_CAST(__m128, const_bb)));
  tt = _mm_add_ps(tt, R_CAST(__m128, f1));
  const __m128i v4 = _mm_and_si128(rr, mask_s);
  __m128i u4 = _mm_add_epi32(rr, i127s);
  u4 = _mm_srli_epi32(u4, 10);
  u4 = _mm_slli_epi32(u4, 23);
  const uint32_t v0 = _mm_cvtsi128_si32(v4);
  // uint32_t v1 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(2)));
  // uint32_t v2 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(4)));
  // uint32_t v3 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(6)));
  // make this work with LLVM
  const uint32_t v1 = _mm_extract_epi16(R_CAST(__m128i, v4), 2);
  const uint32_t v2 = _mm_extract_epi16(R_CAST(__m128i, v4), 4);
  const uint32_t v3 = _mm_extract_epi16(R_CAST(__m128i, v4), 6);

  __m128 t0 = _mm_set_ss(kFloatExpLookup[v0]);
  __m128 t1 = _mm_set_ss(kFloatExpLookup[v1]);
  const __m128 t2 = _mm_set_ss(kFloatExpLookup[v2]);
  const __m128 t3 = _mm_set_ss(kFloatExpLookup[v3]);
  t1 = _mm_movelh_ps(t1, t3);
  t1 = _mm_castsi128_ps(_mm_slli_epi64(_mm_castps_si128(t1), 32));
  t0 = _mm_movelh_ps(t0, t2);
  t0 = _mm_or_ps(t0, t1);
  t0 = _mm_or_ps(t0, _mm_castsi128_ps(u4));
  tt = _mm_mul_ps(tt, t0);
  return R_CAST(VecF, tt);
}
#  endif

CONSTI32(kFmathExpdSbit, 11);
CONSTI32(kFmathExpdS, 1 << kFmathExpdSbit);
CONSTI32(kFmathExpdAdj, kFmathExpdS * 0x3ff);
static const double kFmathExpdC1 = 1.0;
// I assume Mitsunari adjusted these constants slightly (away from 1/6 and 3)
// to reduce empirical error.
static const double kFmathExpdC2 = 0.16666666685227835064;
static const double kFmathExpdC3 = 3.0000000027955394;
static const double kFmathExpdA = S_CAST(double, kFmathExpdS) / kLn2;
static const double kFmathExpdRa = kLn2 / S_CAST(double, kFmathExpdS);

#  ifdef FVEC_32
extern const uint64_t kExpdLookupInt[] __attribute__((aligned(32)));
#  else
extern const uint64_t kExpdLookupInt[] __attribute__((aligned(16)));
#  endif

typedef union {
  double d8;
  uint64_t u8;
} __uni8;

// Assumes x is within the expd() bounds.
HEADER_INLINE double expd_bounded(double x) {
  const uint64_t b = 3ULL << 51;
  // kFmathExpdA is 1 / log(2^{1/2048}), so y := x * kFmathExpdA satisfies
  //   e^x == (2^{1/2048})^y
  // Then, we set d := 1.5 * 2^52 + y.  Recall that IEEE-754 doubles have 52
  // mantissa bits, with an implicit leading '1'.  |y| is far smaller than
  // 2^51, so we know 2^52 < d < 2^53; thus, the exponent bits correspond to a
  // multiplier of 1.  The low mantissa bits of d correspond to y rounded to
  // the nearest (signed) integer; call this value yi.
  __uni8 di;
  di.d8 = x * kFmathExpdA + b;
  // Define z := yi * log(2^{1/2048}).
  // kExpdLookupInt is a 2048-entry table encoding the mantissa bits of
  // 2^{0/2048}, 2^{1/2048}, ..., 2^{2047/2048}.  From this value, we can
  // easily determine e^z.
  const uint64_t iax = kExpdLookupInt[di.u8 & (kFmathExpdS - 1)];

  // Calculate (z - x).  Under the default rounding mode, the absolute value
  // should not be larger than log(2^{1/2048}) / 2.
  const double t = (di.d8 - b) * kFmathExpdRa - x;
  const uint64_t u = ((di.u8 + kFmathExpdAdj) >> kFmathExpdSbit) << 52;
  // We're now ready to use the Taylor expansion of e^x about x=z.
  //   e^x = f(z) - f'(z)(z-x) + f''(z)((z-x)^2)/2 - f'''(z)((z-x)^3)/6 + ...
  // Conveniently, f(z) = f'(z) = f''(z) = f'''(z), so we can factor it out and
  // just multiply by it at the end.
  // We can stop the Taylor expansion after the cubic term since |((z-x)^4)/24|
  // is guaranteed to be smaller than 4e-17, and the (1 - (z-x) + ((z-x)^2)/2 -
  // ((z-x)^3)/6) expression it's being added to is guaranteed to be larger
  // than 0.5, so the quartic term would never affect even the least
  // significant bit.
  const double y = (kFmathExpdC3 - t) * (t * t) * kFmathExpdC2 - t + kFmathExpdC1;

  di.u8 = u | iax;
  return y * di.d8;
}

// Maximum error of this function is ~6 ulps; there's some rounding error in
// the multiply operations, etc.
HEADER_INLINE double expd(double x) {
  if (x <= -708.39641853226414) {
    // Changed from the original bound of -708.39641853226408, since that
    // produced inconsistent results between expd() and expd_v() in my testing.
    // The proper value of this constant depends on the floating-point rounding
    // mode.
    return 0;
  }
  if (x >= 709.78254366799815) {
    // This bound is also changed (from 709.78271289338397), though in this
    // case the change doesn't affect behavior.
    return S_CAST(double, INFINITY);
  }
  return expd_bounded(x);
}

// Requires vector-alignment.
HEADER_INLINE void expd_v(double* px, uintptr_t n) {
  const double b = S_CAST(double, 3ULL << 51);
#  ifdef FVEC_32
  const uintptr_t r = n & 3;
  n &= ~3;
  const __m256d mC1 = _mm256_set1_pd(kFmathExpdC1);
  const __m256d mC2 = _mm256_set1_pd(kFmathExpdC2);
  const __m256d mC3 = _mm256_set1_pd(kFmathExpdC3);
  const __m256d ma = _mm256_set1_pd(kFmathExpdA);
  const __m256d mra = _mm256_set1_pd(kFmathExpdRa);
  const __m256i madj = _mm256_set1_epi64x(kFmathExpdAdj);
  const __m256i maskSbit = _mm256_set1_epi64x(kFmathExpdS - 1);
  const __m256d expMax = _mm256_set1_pd(709.78254366799815);
  const __m256d expMin = _mm256_set1_pd(-708.39641853226414);
  for (size_t i = 0; i < n; i += 4) {
    __m256d x = _mm256_load_pd(px);
    x = _mm256_min_pd(x, expMax);
    x = _mm256_max_pd(x, expMin);

    __m256d d = _mm256_mul_pd(x, ma);
    d = _mm256_add_pd(d, _mm256_set1_pd(b));
    const __m256i adr = _mm256_and_si256(_mm256_castpd_si256(d), maskSbit);
    const __m256i iax = _mm256_i64gather_epi64(R_CAST(const long long*, kExpdLookupInt), adr, 8);
    const __m256d t = _mm256_sub_pd(_mm256_mul_pd(_mm256_sub_pd(d, _mm256_set1_pd(b)), mra), x);
    __m256i u = _mm256_castpd_si256(d);
    u = _mm256_add_epi64(u, madj);
    u = _mm256_srli_epi64(u, kFmathExpdSbit);
    u = _mm256_slli_epi64(u, 52);
    u = _mm256_or_si256(u, iax);
    __m256d y = _mm256_mul_pd(_mm256_sub_pd(mC3, t), _mm256_mul_pd(t, t));
    y = _mm256_mul_pd(y, mC2);
    y = _mm256_add_pd(_mm256_sub_pd(y, t), mC1);
    _mm256_store_pd(px, _mm256_mul_pd(y, _mm256_castsi256_pd(u)));
    px += 4;
  }
#  else
  const uintptr_t r = n & 1;
  n &= ~1;
  const __m128d mC1 = _mm_set1_pd(kFmathExpdC1);
  const __m128d mC2 = _mm_set1_pd(kFmathExpdC2);
  const __m128d mC3 = _mm_set1_pd(kFmathExpdC3);
  const __m128d ma = _mm_set1_pd(kFmathExpdA);
  const __m128d mra = _mm_set1_pd(kFmathExpdRa);
#    if defined(__x86_64__) || defined(_WIN64)
  const __m128i madj = _mm_set1_epi64x(kFmathExpdAdj);
#    else
  const __m128i madj = _mm_set_epi32(0, kFmathExpdAdj, 0, kFmathExpdAdj);
#    endif
  const __m128d expMax = _mm_set1_pd(709.78254366799815);
  const __m128d expMin = _mm_set1_pd(-708.39641853226414);
  for (size_t i = 0; i < n; i += 2) {
    __m128d x = _mm_load_pd(px);
    x = _mm_min_pd(x, expMax);
    x = _mm_max_pd(x, expMin);

    __m128d d = _mm_mul_pd(x, ma);
    d = _mm_add_pd(d, _mm_set1_pd(b));
    const int adr0 = _mm_cvtsi128_si32(_mm_castpd_si128(d)) & (kFmathExpdS - 1);
    const int adr1 = _mm_cvtsi128_si32(_mm_srli_si128(_mm_castpd_si128(d), 8)) & (kFmathExpdS - 1);

    const __m128i iaxL = _mm_castpd_si128(_mm_load_sd(R_CAST(const double*, &(kExpdLookupInt[adr0]))));
    __m128i iax = _mm_castpd_si128(_mm_load_sd(R_CAST(const double*, &(kExpdLookupInt[adr1]))));
    iax = _mm_unpacklo_epi64(iaxL, iax);

    const __m128d t = _mm_sub_pd(_mm_mul_pd(_mm_sub_pd(d, _mm_set1_pd(b)), mra), x);
    __m128i u = _mm_castpd_si128(d);
    u = _mm_add_epi64(u, madj);
    u = _mm_srli_epi64(u, kFmathExpdSbit);
    u = _mm_slli_epi64(u, 52);
    u = _mm_or_si128(u, iax);
    __m128d y = _mm_mul_pd(_mm_sub_pd(mC3, t), _mm_mul_pd(t, t));
    y = _mm_mul_pd(y, mC2);
    y = _mm_add_pd(_mm_sub_pd(y, t), mC1);
    _mm_store_pd(px, _mm_mul_pd(y, _mm_castsi128_pd(u)));
    px += 2;
  }
#  endif
  for (uintptr_t i = 0; i < r; ++i) {
    px[i] = expd(px[i]);
  }
}

// 1 / (1 + e^{-x}).
// Looked into approximating this function directly instead of going through
// expd_bounded(), but I couldn't find a better approach.
HEADER_INLINE double logistic(double x) {
  if (x <= -708.39641853226414) {
    // Even though expd_bounded(-x) returns a finite value down to -709.78...,
    // the reciprocal of that value is only representable as a denormal double,
    // and plink2 disables denormals.
    return 0.0;
  }
  // As soon as e^{-x} is smaller than ~2^{-53}, the result rounds up to 1.
  if (x >= 36.736800569677098) {
    return 1.0;
  }
  return 1.0 / (1.0 + expd_bounded(-x));
}

// In-place vector 1 / (1 + e^{-x}).
// Requires vector-alignment, and clobbers trailing entries.
HEADER_INLINE void logistic_v_unsafe(double* px, uintptr_t n) {
  const double b = S_CAST(double, 3ULL << 51);
#  ifdef FVEC_32
  const __m256d mC1 = _mm256_set1_pd(kFmathExpdC1);
  const __m256d mC2 = _mm256_set1_pd(kFmathExpdC2);
  const __m256d mC3 = _mm256_set1_pd(kFmathExpdC3);
  const __m256d mnega = _mm256_set1_pd(-kFmathExpdA);
  const __m256d mra = _mm256_set1_pd(kFmathExpdRa);
  const __m256i madj = _mm256_set1_epi64x(kFmathExpdAdj);
  const __m256i maskSbit = _mm256_set1_epi64x(kFmathExpdS - 1);
  const __m256d expMax = _mm256_set1_pd(36.736800569677098);
  const __m256d expMin = _mm256_set1_pd(-708.39641853226414);
  for (size_t i = 0; i < n; i += 4) {
    __m256d x = _mm256_load_pd(px);
    x = _mm256_min_pd(x, expMax);
    x = _mm256_max_pd(x, expMin);

    // We want e^{-x}, not e^x.
    __m256d d = _mm256_mul_pd(x, mnega);
    d = _mm256_add_pd(d, _mm256_set1_pd(b));
    const __m256i adr = _mm256_and_si256(_mm256_castpd_si256(d), maskSbit);
    const __m256i iax = _mm256_i64gather_epi64(R_CAST(const long long*, kExpdLookupInt), adr, 8);
    // Subtract -x -> add x.
    const __m256d t = _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(d, _mm256_set1_pd(b)), mra), x);
    __m256i u = _mm256_castpd_si256(d);
    u = _mm256_add_epi64(u, madj);
    u = _mm256_srli_epi64(u, kFmathExpdSbit);
    u = _mm256_slli_epi64(u, 52);
    u = _mm256_or_si256(u, iax);
    __m256d y = _mm256_mul_pd(_mm256_sub_pd(mC3, t), _mm256_mul_pd(t, t));
    y = _mm256_mul_pd(y, mC2);
    y = _mm256_add_pd(_mm256_sub_pd(y, t), mC1);
    const __m256d enegx = _mm256_mul_pd(y, _mm256_castsi256_pd(u));
    // Add 1, then take the reciprocal, to get the final result.
    _mm256_store_pd(px, _mm256_div_pd(mC1, _mm256_add_pd(mC1, enegx)));
    px += 4;
  }
#  else
  const __m128d mC1 = _mm_set1_pd(kFmathExpdC1);
  const __m128d mC2 = _mm_set1_pd(kFmathExpdC2);
  const __m128d mC3 = _mm_set1_pd(kFmathExpdC3);
  const __m128d mnega = _mm_set1_pd(-kFmathExpdA);
  const __m128d mra = _mm_set1_pd(kFmathExpdRa);
#    if defined(__x86_64__) || defined(_WIN64)
  const __m128i madj = _mm_set1_epi64x(kFmathExpdAdj);
#    else
  const __m128i madj = _mm_set_epi32(0, kFmathExpdAdj, 0, kFmathExpdAdj);
#    endif
  const __m128d expMax = _mm_set1_pd(36.736800569677098);
  const __m128d expMin = _mm_set1_pd(-708.39641853226414);
  for (size_t i = 0; i < n; i += 2) {
    __m128d x = _mm_load_pd(px);
    x = _mm_min_pd(x, expMax);
    x = _mm_max_pd(x, expMin);

    __m128d d = _mm_mul_pd(x, mnega);
    d = _mm_add_pd(d, _mm_set1_pd(b));
    const int adr0 = _mm_cvtsi128_si32(_mm_castpd_si128(d)) & (kFmathExpdS - 1);
    const int adr1 = _mm_cvtsi128_si32(_mm_srli_si128(_mm_castpd_si128(d), 8)) & (kFmathExpdS - 1);

    const __m128i iaxL = _mm_castpd_si128(_mm_load_sd(R_CAST(const double*, &(kExpdLookupInt[adr0]))));
    __m128i iax = _mm_castpd_si128(_mm_load_sd(R_CAST(const double*, &(kExpdLookupInt[adr1]))));
    iax = _mm_unpacklo_epi64(iaxL, iax);

    const __m128d t = _mm_add_pd(_mm_mul_pd(_mm_sub_pd(d, _mm_set1_pd(b)), mra), x);
    __m128i u = _mm_castpd_si128(d);
    u = _mm_add_epi64(u, madj);
    u = _mm_srli_epi64(u, kFmathExpdSbit);
    u = _mm_slli_epi64(u, 52);
    u = _mm_or_si128(u, iax);
    __m128d y = _mm_mul_pd(_mm_sub_pd(mC3, t), _mm_mul_pd(t, t));
    y = _mm_mul_pd(y, mC2);
    y = _mm_add_pd(_mm_sub_pd(y, t), mC1);
    const __m128d enegx = _mm_mul_pd(y, _mm_castsi128_pd(u));
    _mm_store_pd(px, _mm_div_pd(mC1, _mm_add_pd(mC1, enegx)));
    px += 2;
  }
#  endif
}

HEADER_INLINE void logistic_v(double* px, uintptr_t n) {
  const uintptr_t r = n % kDoublePerDVec;
  n = RoundDownPow2(n, kDoublePerDVec);
  logistic_v_unsafe(px, n);
  px = &(px[n]);
  for (uintptr_t i = 0; i < r; ++i) {
    px[i] = logistic(px[i]);
  }
}

#else  // !__LP64__
HEADER_INLINE double expd(double x) {
  return exp(x);
}

HEADER_INLINE void expd_v(double* px, uintptr_t n) {
  for (uintptr_t ulii = 0; ulii != n; ++ulii) {
    px[ulii] = exp(px[ulii]);
  }
}

HEADER_INLINE double logistic(double x) {
  return 1.0 / (1.0 + exp(-x));
}

HEADER_INLINE void logistic_v_unsafe(double* px, uintptr_t n) {
  for (uintptr_t ulii = 0; ulii != n; ++ulii) {
    px[ulii] = logistic(px[ulii]);
  }
}

HEADER_INLINE void logistic_v(double* px, uintptr_t n) {
  logistic_v_unsafe(px, n);
}
#endif



#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_FMATH_H__
