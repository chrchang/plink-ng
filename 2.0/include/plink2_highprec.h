#ifndef __PLINK2_HIGHPREC_H__
#define __PLINK2_HIGHPREC_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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

#include "plink2_base.h"

#ifdef IGNORE_BUNDLED_MINI_GMP
#  include <gmp.h>  // IWYU pragma: export
#endif
#include <math.h>

#ifndef IGNORE_BUNDLED_MINI_GMP
#  include "../mini-gmp/mini-gmp.h"  // IWYU pragma: export
#endif
#include "plink2_float.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Support for computations requiring more precision than double/int64.

// Portable "double-double" operations supporting high-accuracy log-likelihood
// calculations, based on a small subset of the QD library
// (https://github.com/BL-highprecision/QD ).  See LICENSE.QD for that
// library's BSD-3-Clause-LBNL license.

typedef struct dd_real_struct {
  double x[2];
} dd_real;

extern const dd_real _ddr_log2;

CONSTI32(_ddr_n_ln_fact, 256);
extern const dd_real _ddr_ln_fact[_ddr_n_ln_fact];

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

// Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|.
HEADER_INLINE double qd_quick_two_sum(double a, double b, double* errp) {
  const double s = a + b;
  *errp = b - (s - a);
  return s;
}

// Computes fl(a-b) and err(a-b).  Assumes |a| >= |b|.
/*
HEADER_INLINE double qd_quick_two_diff(double a, double b, double* errp) {
  const double s = a - b;
  *errp = (a - s) - b;
  return s;
}
*/

// Computes fl(a+b) and err(a+b).
HEADER_INLINE double qd_two_sum(double a, double b, double *errp) {
  const double s = a + b;
  const double bb = s - a;
  *errp = (a - (s - bb)) + (b - bb);
  return s;
}

// Computes fl(a-b) and err(a-b).
HEADER_INLINE double qd_two_diff(double a, double b, double* errp) {
  const double s = a - b;
  const double bb = s - a;
  *errp = (a - (s - bb)) - (b + bb);
  return s;
}

#ifndef USE_FMA
// Computes high word and lo word of a
HEADER_INLINE void qd_split(double a, double* hip, double* lop) {
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    const double temp = _QD_SPLITTER * a;
    *hip = temp - (temp - a);
    *lop = a - *hip;
    *hip *= 268435456.0;          // 2^28
    *lop *= 268435456.0;          // 2^28
  } else {
    const double temp = _QD_SPLITTER * a;
    *hip = temp - (temp - a);
    *lop = a - *hip;
  }
}
#endif

// Computes fl(a*b) and err(a*b).
HEADER_INLINE double qd_two_prod(double a, double b, double* errp) {
#ifdef USE_FMA
  const double p = a * b;
  *errp = fma(a, b, -p);
  return p;
#else
  double a_hi, a_lo, b_hi, b_lo;
  const double p = a * b;
  qd_split(a, &a_hi, &a_lo);
  qd_split(b, &b_hi, &b_lo);
  *errp = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

// Computes fl(a*a) and err(a*a).  Faster than the above method.
HEADER_INLINE double qd_two_sqr(double a, double *errp) {
#ifdef USE_FMA
  const double p = a * a;
  *errp = fma(a, a, -p);
  return p;
#else
  double hi, lo;
  double q = a * a;
  qd_split(a, &hi, &lo);
  *errp = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

HEADER_INLINE dd_real ddr_maked(const double a) {
  dd_real retval;
  retval.x[0] = a;
  retval.x[1] = 0.0;
  return retval;
}

HEADER_INLINE dd_real ddr_make(const double a, const double b) {
  dd_real retval;
  retval.x[0] = a;
  retval.x[1] = b;
  return retval;
}

HEADER_INLINE dd_real ddr_addd(const dd_real a, double b) {
  double s1, s2;
  s1 = qd_two_sum(a.x[0], b, &s2);
  s2 += a.x[1];
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}

// Only satisfies Cray-style error bound.
HEADER_INLINE dd_real ddr_sloppy_add(const dd_real a, const dd_real b) {
  double s, e;

  s = qd_two_sum(a.x[0], b.x[0], &e);
  e += a.x[1] + b.x[1];
  s = qd_quick_two_sum(s, e, &e);
  return ddr_make(s, e);
}

/*
HEADER_INLINE dd_real ddr_ieee_add(const dd_real a, const dd_real b) {
  double s1, s2, t1, t2;

  s1 = qd_two_sum(a.x[0], b.x[0], &s2);
  t1 = qd_two_sum(a.x[1], b.x[1], &t2);
  s2 += t1;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  s2 += t2;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}
*/

HEADER_INLINE dd_real ddr_negate(const dd_real a) {
  return ddr_make(-a.x[0], -a.x[1]);
}

// double-double - double
HEADER_INLINE dd_real ddr_subd(const dd_real a, double b) {
  double s1, s2;
  s1 = qd_two_diff(a.x[0], b, &s2);
  s2 += a.x[1];
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}

HEADER_INLINE dd_real ddr_sloppy_sub(const dd_real a, const dd_real b) {
  double s, e;
  s = qd_two_diff(a.x[0], b.x[0], &e);
  e += a.x[1];
  e -= b.x[1];
  s = qd_quick_two_sum(s, e, &e);
  return ddr_make(s, e);
}

/*
HEADER_INLINE dd_real ddr_ieee_sub(const dd_real a, const dd_real b) {
  double s1, s2, t1, t2;
  s1 = qd_two_diff(a.x[0], b.x[0], &s2);
  t1 = qd_two_diff(a.x[1], b.x[1], &t2);
  s2 += t1;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  s2 += t2;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}
*/

HEADER_INLINE dd_real ddr_ldexp(const dd_real a, int32_t expi) {
  return ddr_make(ldexp(a.x[0], expi), ldexp(a.x[1], expi));
}

// double-double * double, where double is a power of 2.
HEADER_INLINE dd_real ddr_mul_pwr2(const dd_real a, double b) {
  return ddr_make(a.x[0] * b, a.x[1] * b);
}

HEADER_INLINE dd_real ddr_muld(const dd_real a, double b) {
  double p1, p2;

  p1 = qd_two_prod(a.x[0], b, &p2);
  p2 += a.x[1] * b;
  p1 = qd_quick_two_sum(p1, p2, &p2);
  return ddr_make(p1, p2);
}

HEADER_INLINE dd_real ddr_mul(const dd_real a, const dd_real b) {
  double p1, p2;

  p1 = qd_two_prod(a.x[0], b.x[0], &p2);
  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
  p1 = qd_quick_two_sum(p1, p2, &p2);
  return ddr_make(p1, p2);
}

HEADER_INLINE dd_real ddr_divd(const dd_real a, double b) {

  double q1, q2;
  double p1, p2;
  double s, e;
  dd_real r;

  q1 = a.x[0] / b;  // approximate quotient

  p1 = qd_two_prod(q1, b, &p2);
  s = qd_two_diff(a.x[0], p1, &e);
  e += a.x[1];
  e -= p2;

  // get next approximation
  q2 = (s + e) / b;

  // renormalize
  r.x[0] = qd_quick_two_sum(q1, q2, &r.x[1]);

  return r;
}

/*
HEADER_INLINE dd_real ddr_sloppy_div(const dd_real a, const dd_real b) {
  double s1, s2;
  double q1, q2;
  dd_real r;

  q1 = a.x[0] / b.x[0];  // approximate quotient

  r = ddr_muld(b, q1);
  s1 = qd_two_diff(a.x[0], r.x[0], &s2);
  s2 -= r.x[1];
  s2 += a.x[1];

  // get next approximation
  q2 = (s1 + s2) / b.x[0];

  // renormalize
  r.x[0] = qd_quick_two_sum(q1, q2, &r.x[1]);
  return r;
}
*/

HEADER_INLINE dd_real ddr_sqr(const dd_real a) {
  double p1, p2;
  double s1, s2;
  p1 = qd_two_sqr(a.x[0], &p2);
  p2 += 2.0 * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  s1 = qd_quick_two_sum(p1, p2, &s2);
  return ddr_make(s1, s2);
}

HEADER_INLINE int32_t ddr_is_zero(const dd_real a) {
  return (a.x[0] == 0.0);
}

HEADER_INLINE int32_t ddr_is_one(const dd_real a) {
  return (a.x[0] == 1.0) && (a.x[1] == 0.0);
}


// Cray error bound is fine for our log-factorial calculation and many other
// applications, but there are a few scenarios (see e.g.
//   https://people.eecs.berkeley.edu/~demmel/cs267/lecture21/lecture21.html
// ) where the slower ieee_add's "don't make catastrophic cancellation any
// worse than it has to be" behavior matters.
HEADER_INLINE dd_real ddr_add(const dd_real a, const dd_real b) {
  return ddr_sloppy_add(a, b);
}

HEADER_INLINE dd_real ddr_sub(const dd_real a, const dd_real b) {
  return ddr_sloppy_sub(a, b);
}


dd_real ddr_exp(const dd_real a);

HEADER_INLINE dd_real ddr_add3(const dd_real a, const dd_real b, const dd_real c) {
  return ddr_add(ddr_add(a, b), c);
}

HEADER_INLINE dd_real ddr_add4(const dd_real a, const dd_real b, const dd_real c, const dd_real d) {
  return ddr_add(ddr_add(ddr_add(a, b), c), d);
}

HEADER_INLINE dd_real ddr_add5(const dd_real a, const dd_real b, const dd_real c, const dd_real d, const dd_real e) {
  return ddr_add(ddr_add(ddr_add(ddr_add(a, b), c), d), e);
}

dd_real ddr_lfact(double xx);

HEADER_INLINE dd_real ddr_add_lfacts(const double a, const double b) {
  return ddr_add(ddr_lfact(a), ddr_lfact(b));
}

HEADER_INLINE dd_real ddr_add3_lfacts(const double a, const double b, const double c) {
  return ddr_add3(ddr_lfact(a), ddr_lfact(b), ddr_lfact(c));
}

HEADER_INLINE dd_real ddr_add4_lfacts(const double a, const double b, const double c, const double d) {
  return ddr_add4(ddr_lfact(a), ddr_lfact(b), ddr_lfact(c), ddr_lfact(d));
}

HEADER_INLINE dd_real ddr_add5_lfacts(const double a, const double b, const double c, const double d, const double e) {
  return ddr_add5(ddr_lfact(a), ddr_lfact(b), ddr_lfact(c), ddr_lfact(d), ddr_lfact(e));
}

// Preconditions:
// - numer_factorial_args[] and denom_factorial_args[] are
//   not-necessarily-sorted lists of length ffac_ct, describing a quotient of
//   factorial-products.  (If one list is longer than the other, just pad the
//   other with zeroes.)
// - pow2 is a power of 2 to multiply by at the end.
//
// This function errors out iff memory allocation fails.
//
// Postconditions on success:
// - *cmp_resultp is set to a positive value if the fraction > 1, a negative
//   value if the fraction < 1, and zero if it's exactly 1.
// - *dbl_ptr is the double representation of the fraction, error limited to
//   1-2 ulps.
// - numer_factorial_args[] and denom_factorial_args[] are sorted in
//   nondecreasing order.
//
// possible todo: accept precomputed ddr_lfact table
BoolErr CompareFactorialProducts(uint32_t ffac_ct, int64_t pow2, int64_t numer_pow2, uint32_t* numer_factorial_args, uint32_t* denom_factorial_args, dd_real* numer_ddr_ptr, mp_limb_t** gmp_wkspacep, uintptr_t* gmp_wkspace_limb_ctp, intptr_t* cmp_resultp, double* dbl_ptr);

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_HIGHPREC_H__
