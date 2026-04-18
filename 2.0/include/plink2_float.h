#ifndef __PLINK2_FLOAT_H__
#define __PLINK2_FLOAT_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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


// Basic floating-point constants and portable definitions.

#include <math.h>

#include "plink2_base.h"

// apparently these aren't always defined in limits.h
#ifndef DBL_MAX
#  define DBL_MAX 1.7976931348623157e308
#endif
#ifndef FLT_MAX
#  define FLT_MAX S_CAST(float, 3.40282347e38)
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

static const double FLT_MAX_D = 3.4028235677973362e38;
// INFINITY is usually a float, but apparently not in w64-mingw32-g++.
static const float INFINITY_F = S_CAST(float, INFINITY);
static const double INFINITY_D = S_CAST(double, INFINITY);

// Notes on float64 ('double') precision:
// * 1 + 2^{-52} is the smallest float64 greater than 1, and 1 - 2^{-53} is the
//   largest float64 less than 1.  Since plink2 flushes denormals to zero on
//   the main target platforms, a 'ULP' ("unit in the last place") can be
//   assumed to be between 2^{-52} and 2^{-53} times the total value.
//   When manually updating an epsilon value that's ever added to 1, multiples
//   of 2^{-52} should be used.  If epsilon is only ever subtracted from 1, it
//   is ok to use multiples of 2^{-53}.
// * With default rounding behavior, the maximum rounding error from a basic
//   arithmetic operation is 0.5 ULP, i.e. ~2^{-53} times the result.
//   plink2 functions are written to assume default rounding, and to avoid
//   assumptions about denormal handling.
// * Unfortunately, accuracy of functions like exp and log is currently
//   platform-dependent, and it isn't realistically worthwhile to purge that
//   source of platform variation when we're continuing to use different linear
//   algebra libraries on different platforms.
//   The best survey I've found of the current state of the world is
//     Gladman B, Innocente V, Mather J, Zimmermann P, (2024) Accuracy of
//     Mathematical Functions in Single, Double, Double Extended, and Quadruple
//     Precision.  hal-03141101v7.
//     https://inria.hal.science/hal-03141101v7/document
//   It reports the following maximum observed errors (not an exhaustive check
//   due to the size of the search space):
//     acos: 1.53 ULP (CUDA 12.2.1)
//     cbrt: 1.53e22 ULP (!) (AMD 4.2); next-worst is GNU libc 2.40 at 3.67 ULP
//           AMD failure mode only applies to denormals.
//     cos/sin: Inf (!) (LLVM 18.1.8); next-worst is CUDA 12.2.1 at 1.52 ULP
//              LLVM failure mode only applies to very large arguments.
//     exp: 1.5 ULP (MSVC 2022)
//     expm1: 3.06 ULP (MSVC 2022)
//     log: 0.946 ULP (Newlib 4.4.0, OpenLibm 0.8.3, FreeBSD 14.1)
//     log1p: 1.74 ULP (ArmPL 24.04)
//     sqrt: 0.5 ULP
static const double k2m32 = 1.0 / (1LL << 32);
static const double k2m50 = 1.0 / (1LL << 50);
static const double k2m52 = 1.0 / (1LL << 52);
static const double k2m53 = 1.0 / (1LL << 53);
static const double k2p50 = 1.0 * (1LL << 50);
static const double k2p64 = 4.0 * (1LL << 62);
static const double k2p100 = k2p50 * k2p50;
static const double k2p200 = k2p100 * k2p100;
static const double k2p400 = k2p200 * k2p200;
static const double k2p800 = k2p400 * k2p400;
static const double kE = 2.7182818284590452;
static const double kLn2 = 0.6931471805599453;
static const double kLn10 = 2.3025850929940457;
static const double kLnNormalMin = -708.3964185322641;
static const double kLnSqrtPi = 0.5723649429247001;
static const double kLnSqrt2Pi = 0.91893853320467278056;
static const double kPi = 3.1415926535897932;
static const double kRecipE = 0.36787944117144233;
static const double kRecipLn10 = 0.43429448190325176;
static const double kSqrt2 = 1.4142135623730951;

// Some more negative powers of 2, mostly used as tolerances for floating-point
// approximate-equality checks.
static const double k2m21 = 1.0 / (1 << 21);
static const double k2m30 = 1.0 / (1 << 30);
static const double k2m35 = 1.0 / (1LL << 35);
static const double k2m44 = 1.0 / (1LL << 44);
static const double k2m60 = 1.0 / (1LL << 60);

static const double kBigEpsilon = k2m21;  // must be >= sqrt(kSmallEpsilon)
static const double kEpsilon = k2m30;
static const double kSmallEpsilon = k2m44;

#if defined(__cplusplus)
#  if __cplusplus >= 201103L
HEADER_INLINE bool isfinite_f(float fxx) {
  using namespace std;
  return isfinite(fxx);
}

HEADER_INLINE bool isfinite_d(double dxx) {
  using namespace std;
  return isfinite(dxx);
}
#  else
#    ifdef isfinite
#      define isfinite_f isfinite
#      define isfinite_d isfinite
#    else
HEADER_INLINE bool isfinite_f(float fxx) {
  return (fxx == fxx) && (fxx != INFINITY_F) && (fxx != -INFINITY_F);
}

HEADER_INLINE bool isfinite_d(double dxx) {
  return (dxx == dxx) && (dxx != INFINITY_D)) && (dxx != -INFINITY_D));
}
#    endif
#  endif
#else
#  define isfinite_f isfinite
#  define isfinite_d isfinite
#endif

// Floating-point environment setup.
//
// - plink2 is designed to work properly when denormals are flushed to zero, so
//   setting that flag is an easy performance win.
// - On 32-bit systems, the implementation of plink2_highprec's double-double
//   type assumes intermediate results are rounded by the processor to 64 bits,
//   not 80.  ...honestly it's fine to just drop 32-bit support at this point,
//   but passing "-msse -mfpmath=sse" (which doesn't work with chips older than
//   ~2000) to gcc is a simple fix.

void flush_denormals();

// We treat minor floating-point differences across platforms as inevitable,
// since we don't want to lock ourselves to a primitive linear algebra library
// implementation.  But we try to avoid introducing such differences without
// good reason.
//
// There is a fused-multiply-add (FMA) instruction on ARM and x86_64 v3 (AVX2),
// but not older x86.  This can speed up some critical loops by ~30% and/or
// improve accuracy when used properly.  That can be enough benefit to justify
// the cost of more minor differences between x86_64 v1 and v3 floating-point
// results.  The C language standard allows compilers to automatically
// "contract" (x*y)+z in a single expression to fma(x,y,z) when the target has
// efficient hardware support for the latter, and both clang 14+ and gcc were
// doing this under some of our build configurations.
//
// Unfortunately, this can break correctness of code like that in
// plink2_highprec, and even the latest gcc does not honor the #pragma that is
// supposed to be able to shield such code.  That leaves two options:
// 1. Compile files like include/plink2_highprec.cc with -ffp-contract=off, and
//    allow contraction elsewhere.
// 2. Compile everything with -ffp-contract=off, and then insert allowed-FMAs
//    by hand.
// I've decided to go with option 2 since I'm generally trying to innovate at
// least a little bit on floating-point accuracy, the amount of time-critical
// floating-point code to review is relatively small, and my preexisting code
// was not written with much awareness of how automatic contractions could
// endanger accuracy.
// (Simple counterintuitive example: (x1*y1)-(x2*y2) can be contracted to
// fma(x1, y1, -x2*y2), potentially yielding a nonzero result when x1=x2 and
// y1=y2.  See also
//   https://stackoverflow.com/questions/73985098/clang-14-0-0-floating-point-optimizations
// .)
#if defined (__FMA__) || defined(__ARM_FEATURE_FMA)
#  define USE_FMA
// Note that this typically isn't worthwhile when multiplying by a power-of-2
// constant (doesn't require a regular multiply op), and that careful attention
// must be paid to dependency chains (if third argument isn't naturally
// available by the time the first two are, FMA is likely to slow things down;
// and in general, using multiple accumulators is a lot more helpful than
// inserting FP contractions).
HEADER_INLINE double prefer_fma(double a, double b, double c) {
  return fma(a, b, c);
}

HEADER_INLINE float prefer_fmaf(float a, float b, float c) {
  return fmaf(a, b, c);
}
#else
HEADER_INLINE double prefer_fma(double a, double b, double c) {
  return a * b + c;
}

HEADER_INLINE float prefer_fmaf(float a, float b, float c) {
  return a * b + c;
}
#endif

// Efficient alternatives to ceil() for nonnegative numbers.

// limit is assumed to be a positive int32.
HEADER_INLINE double ceil_smalleps_limit32(double xx, double limit) {
  if (xx > limit) {
    return limit;
  }
  return S_CAST(int32_t, (xx + 1) * (1 - kSmallEpsilon));
}

HEADER_INLINE double ceil_smalleps(double xx) {
  return S_CAST(int64_t, (xx + 1) * (1 - kSmallEpsilon));
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_FLOAT_H__
