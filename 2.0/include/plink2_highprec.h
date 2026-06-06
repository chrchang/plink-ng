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

#include <math.h>

#include "plink2_float.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Portable "double-double" and "triple-double" operations, primarily
// supporting high-accuracy log-likelihood calculations, based on a subset of
// the QD library (https://github.com/BL-highprecision/QD ).  See LICENSE.QD
// for that library's BSD-3-Clause-LBNL license.

// Initially, I ported the qd_real type.  However, from the benchmark results
// at e.g.
//   https://perso.ens-lyon.fr/gilles.villard/hplll/fpbench.html
// (and I had similar results) it is clear that dd_real is very
// well-implemented, but the approach doesn't scale well to qd_reals (GNU MPFR
// becomes clearly better at that precision).  Fortunately, it is
// straightforward to modify the qd_real functions to work with td_reals.


typedef struct td_real_struct {
  double x[3];
} td_real;

extern const td_real _tdr_log2;
extern const td_real _tdr_log05;
static const double _tdr_eps = 1.0947644252537633e-47; // = 2^-156

typedef struct dd_real_struct {
  double x[2];
} dd_real;

HEADER_INLINE void swap_tdr(td_real* ap, td_real* bp) {
  const td_real swaptmp = *ap;
  *ap = *bp;
  *bp = swaptmp;
}

HEADER_INLINE void swap_ddr(dd_real* ap, dd_real* bp) {
  const dd_real swaptmp = *ap;
  *ap = *bp;
  *bp = swaptmp;
}

extern const dd_real _ddr_e;
extern const dd_real _ddr_log2;
extern const dd_real _ddr_log05;
extern const dd_real _ddr_64log2;

CONSTI32(_tdr_n_ln_fact, 256);
extern const td_real _tdr_ln_fact[_tdr_n_ln_fact];

// Fast float64-accuracy version.
// Assumes xx is a nonnegative integer.
double Lfact(double xx);

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

HEADER_CINLINE dd_real ddr_maked(const double a) {
  const dd_real retval = {{a, 0.0}};
  return retval;
}

HEADER_CINLINE dd_real ddr_makei(const int64_t a) {
  const double hi = S_CAST(double, a);
  const dd_real retval = {{hi, S_CAST(double, a - S_CAST(int64_t, hi))}};
  return retval;
}

HEADER_CINLINE dd_real ddr_makeu64(const uint64_t a) {
  const double hi = S_CAST(double, a);
  const dd_real retval = {{hi, S_CAST(double, S_CAST(int64_t, a - S_CAST(uint64_t, hi)))}};
  return retval;
}

HEADER_CINLINE dd_real ddr_make_td(const td_real a) {
  const dd_real retval = {{a.x[0], a.x[1]}};
  return retval;
}

HEADER_CINLINE dd_real ddr_make(const double a, const double b) {
  const dd_real retval = {{a, b}};
  return retval;
}

HEADER_INLINE dd_real ddr_add2d(double a, double b) {
  double s2;
  double s1 = qd_two_sum(a, b, &s2);
  return ddr_make(s1, s2);
}

HEADER_INLINE dd_real ddr_addd(const dd_real a, double b) {
  double s2;
  double s1 = qd_two_sum(a.x[0], b, &s2);
  s2 += a.x[1];
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}

// Only satisfies Cray-style error bound.
HEADER_INLINE dd_real ddr_sloppy_add(const dd_real a, const dd_real b) {
  double e;
  double s = qd_two_sum(a.x[0], b.x[0], &e);
  e += a.x[1] + b.x[1];
  s = qd_quick_two_sum(s, e, &e);
  return ddr_make(s, e);
}

/*
HEADER_INLINE dd_real ddr_ieee_add(const dd_real a, const dd_real b) {
  double s2;
  double t2;
  double s1 = qd_two_sum(a.x[0], b.x[0], &s2);
  double t1 = qd_two_sum(a.x[1], b.x[1], &t2);
  s2 += t1;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  s2 += t2;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}
*/

HEADER_CINLINE dd_real ddr_negate(const dd_real a) {
  return ddr_make(-a.x[0], -a.x[1]);
}

// double-double - double
HEADER_INLINE dd_real ddr_subd(const dd_real a, double b) {
  double s2;
  double s1 = qd_two_diff(a.x[0], b, &s2);
  s2 += a.x[1];
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}

HEADER_INLINE dd_real ddr_sloppy_sub(const dd_real a, const dd_real b) {
  double e;
  double s = qd_two_diff(a.x[0], b.x[0], &e);
  e += a.x[1];
  e -= b.x[1];
  s = qd_quick_two_sum(s, e, &e);
  return ddr_make(s, e);
}

HEADER_INLINE dd_real ddr_ieee_sub(const dd_real a, const dd_real b) {
  double s2;
  double t2;
  double s1 = qd_two_diff(a.x[0], b.x[0], &s2);
  double t1 = qd_two_diff(a.x[1], b.x[1], &t2);
  s2 += t1;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  s2 += t2;
  s1 = qd_quick_two_sum(s1, s2, &s2);
  return ddr_make(s1, s2);
}

HEADER_INLINE dd_real ddr_ldexp(const dd_real a, int32_t expi) {
  return ddr_make(ldexp(a.x[0], expi), ldexp(a.x[1], expi));
}

// double-double * double, where double is a power of 2.
HEADER_CINLINE dd_real ddr_mul_pwr2(const dd_real a, double b) {
  return ddr_make(a.x[0] * b, a.x[1] * b);
}

// This is (i) very efficient when FMA is available, and (ii) is error-free
// when a and b are e.g. integers < 2^53.
HEADER_INLINE dd_real ddr_mul2d(double a, double b) {
  double p2;
  double p1 = qd_two_prod(a, b, &p2);
  return ddr_make(p1, p2);
}

HEADER_INLINE dd_real ddr_muld(const dd_real a, double b) {
  double p2;
  double p1 = qd_two_prod(a.x[0], b, &p2);
  p2 += a.x[1] * b;
  p1 = qd_quick_two_sum(p1, p2, &p2);
  return ddr_make(p1, p2);
}

HEADER_INLINE dd_real ddr_mul(const dd_real a, const dd_real b) {
  double p2;
  double p1 = qd_two_prod(a.x[0], b.x[0], &p2);
  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
  p1 = qd_quick_two_sum(p1, p2, &p2);
  return ddr_make(p1, p2);
}

HEADER_INLINE dd_real ddr_divd(const dd_real a, double b) {

  double p2;
  double e;
  dd_real r;

  const double q1 = a.x[0] / b;  // approximate quotient

  const double p1 = qd_two_prod(q1, b, &p2);
  double s = qd_two_diff(a.x[0], p1, &e);
  e += a.x[1];
  e -= p2;

  // get next approximation
  const double q2 = (s + e) / b;

  // renormalize
  r.x[0] = qd_quick_two_sum(q1, q2, &r.x[1]);

  return r;
}

/*
HEADER_INLINE dd_real ddr_sloppy_div(const dd_real a, const dd_real b) {
  double s2;
  dd_real r;

  const double q1 = a.x[0] / b.x[0];  // approximate quotient

  r = ddr_muld(b, q1);
  const double s1 = qd_two_diff(a.x[0], r.x[0], &s2);
  s2 -= r.x[1];
  s2 += a.x[1];

  // get next approximation
  const double q2 = (s1 + s2) / b.x[0];

  // renormalize
  r.x[0] = qd_quick_two_sum(q1, q2, &r.x[1]);
  return r;
}
*/

static inline dd_real ddr_accurate_div(const dd_real a, const dd_real b) {
  double q1 = a.x[0] / b.x[0];  // approximate quotient

  dd_real r = ddr_ieee_sub(a, ddr_muld(b, q1));

  double q2 = r.x[0] / b.x[0];
  r = ddr_ieee_sub(r, ddr_muld(b, q2));

  const double q3 = r.x[0] / b.x[0];

  q1 = qd_quick_two_sum(q1, q2, &q2);
  return ddr_addd(ddr_make(q1, q2), q3);
}

HEADER_INLINE dd_real ddr_sqr(const dd_real a) {
  double p2;
  double s2;
  const double p1 = qd_two_sqr(a.x[0], &p2);
  p2 += 2.0 * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  const double s1 = qd_quick_two_sum(p1, p2, &s2);
  return ddr_make(s1, s2);
}

HEADER_CINLINE int32_t ddr_is_zero(const dd_real a) {
  return (a.x[0] == 0.0);
}

HEADER_CINLINE int32_t ddr_is(const dd_real a, double b) {
  return (a.x[0] == b) && (a.x[1] == 0.0);
}


HEADER_CINLINE int32_t ddr_ltd(const dd_real a, const double b) {
  return (a.x[0] < b) || ((a.x[0] == b) && (a.x[1] < 0));
}

HEADER_CINLINE int32_t ddr_lt(const dd_real a, const dd_real b) {
  return (a.x[0] < b.x[0]) || ((a.x[0] == b.x[0]) && (a.x[1] < b.x[1]));
}

HEADER_CINLINE int32_t ddr_leqd(const dd_real a, const double b) {
  return (a.x[0] < b) || ((a.x[0] == b) && (a.x[1] <= 0));
}

HEADER_CINLINE int32_t ddr_leq(const dd_real a, const dd_real b) {
  return (a.x[0] < b.x[0]) || ((a.x[0] == b.x[0]) && (a.x[1] <= b.x[1]));
}

HEADER_CINLINE int32_t ddr_geqd(const dd_real a, const double b) {
  return (a.x[0] > b) || ((a.x[0] == b) && (a.x[1] >= 0));
}

HEADER_CINLINE int32_t ddr_geq(const dd_real a, const dd_real b) {
  return (a.x[0] > b.x[0]) || ((a.x[0] == b.x[0]) && (a.x[1] >= b.x[1]));
}

HEADER_CINLINE int32_t ddr_gtd(const dd_real a, const double b) {
  return (a.x[0] > b) || ((a.x[0] == b) && (a.x[1] > 0));
}

HEADER_CINLINE int32_t ddr_gt(const dd_real a, const dd_real b) {
  return (a.x[0] > b.x[0]) || ((a.x[0] == b.x[0]) && (a.x[1] > b.x[1]));
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

dd_real ddr_log(const dd_real a);

dd_real ddr_expm1(const dd_real a);

dd_real ddr_log1p(const dd_real a);


// Try to put smaller-magnitude values first (or just use ddr_sort_and_add()).
HEADER_INLINE dd_real ddr_add3(const dd_real a, const dd_real b, const dd_real c) {
  return ddr_add(ddr_add(a, b), c);
}

HEADER_INLINE dd_real ddr_add4(const dd_real a, const dd_real b, const dd_real c, const dd_real d) {
  return ddr_add(ddr_add(ddr_add(a, b), c), d);
}

HEADER_INLINE dd_real ddr_add5(const dd_real a, const dd_real b, const dd_real c, const dd_real d, const dd_real e) {
  return ddr_add(ddr_add(ddr_add(ddr_add(a, b), c), d), e);
}

// Assumes xx is a nonnegative integer < 2^52.
dd_real ddr_lfact(double xx);

HEADER_INLINE dd_real ddr_add_lfacts(const double a, const double b) {
  return ddr_add(ddr_lfact(a), ddr_lfact(b));
}

HEADER_INLINE dd_real ddr_add3_lfacts(const double a, const double b1, const double b2) {
  // Common case seems to be that there's a known-nonmaximal argument caller
  // can put first, but order of the other two varies.
  const dd_real lfact_a_ddr = ddr_lfact(a);
  const dd_real lfact_b1_ddr = ddr_lfact(b1);
  const dd_real lfact_b2_ddr = ddr_lfact(b2);
  if (b1 < b2) {
    return ddr_add3(lfact_a_ddr, lfact_b1_ddr, lfact_b2_ddr);
  } else {
    return ddr_add3(lfact_a_ddr, lfact_b2_ddr, lfact_b1_ddr);
  }
}

// ddr_add{4,5}_lfacts() removed to encourage sort-and-add for those cases.

// Sorts ddrs in place from least to greatest magnitude, and then adds.
// Assumes ct positive.
dd_real ddr_sort_and_add(uint32_t ct, dd_real* ddrs);

HEADER_INLINE dd_real ddr_sort_and_add3(dd_real a_ddr, dd_real b_ddr, dd_real c_ddr) {
  dd_real args[3];
  args[0] = a_ddr;
  args[1] = b_ddr;
  args[2] = c_ddr;
  return ddr_sort_and_add(3, args);
}

dd_real ddr_sort_and_add_lfacts(uint32_t ct, double* args);

HEADER_INLINE dd_real ddr_sort_and_add_3_lfacts(double a, double b, double c) {
  double args[3];
  args[0] = a;
  args[1] = b;
  args[2] = c;
  return ddr_sort_and_add_lfacts(3, args);
}

HEADER_INLINE dd_real ddr_sort_and_add_4_lfacts(double a, double b, double c, double d) {
  double args[4];
  args[0] = a;
  args[1] = b;
  args[2] = c;
  args[3] = d;
  return ddr_sort_and_add_lfacts(4, args);
}

HEADER_INLINE dd_real ddr_sort_and_add_5_lfacts(double a, double b, double c, double d, double e) {
  double args[5];
  args[0] = a;
  args[1] = b;
  args[2] = c;
  args[3] = d;
  args[4] = e;
  return ddr_sort_and_add_lfacts(5, args);
}

// lnprob_ddr <= 0, mult < 2^52, exp(lnprob_ddr) * mult < 0.5 or so.
// Avoids intermediate underflow when exp(lnprob_ddr) would underflow, but
// final result does not underflow.  Doesn't try to preserve last bit of
// accuracy, but does try to avoid throwing away up to 10 bits when logp is
// false and return value is positive and <= e^{-512}.
HEADER_INLINE double join_log_and_nonlog(dd_real lnprob_ddr, double mult, uint32_t logp) {
  if (logp) {
    return lnprob_ddr.x[0] + log(mult);
  }
  return mult * k2m64 * ddr_exp(ddr_add(lnprob_ddr, _ddr_64log2)).x[0];
}


HEADER_INLINE void qd_renorm3(double* c0p, double* c1p, double* c2p) {
  if ((*c0p == INFINITY_D) || (*c0p == -INFINITY_D)) {
    return;
  }

  double s2 = 0.0;

  double s0 = qd_quick_two_sum(*c1p, *c2p, c2p);
  *c0p = qd_quick_two_sum(*c0p, s0, c1p);

  s0 = *c0p;
  double s1 = *c1p;
  if (s1 != 0.0) {
    s1 = qd_quick_two_sum(s1, *c2p, &s2);
  } else {
    s0 = qd_quick_two_sum(s0, *c2p, &s1);
  }

  *c0p = s0;
  *c1p = s1;
  *c2p = s2;
}

HEADER_INLINE void qd_renorm4(double* c0p, double* c1p, double* c2p, double* c3p) {
  if ((*c0p == INFINITY_D) || (*c0p == -INFINITY_D)) {
    return;
  }

  double s2 = 0.0;

  double s0 = qd_quick_two_sum(*c2p, *c3p, c3p);
  s0 = qd_quick_two_sum(*c1p, s0, c2p);
  *c0p = qd_quick_two_sum(*c0p, s0, c1p);

  s0 = *c0p;
  double s1 = *c1p;

  if (s1 != 0.0) {
    s1 = qd_quick_two_sum(s1, *c2p, &s2);
    if (s2 != 0.0) {
      s2 += *c3p;
    } else {
      s1 = qd_quick_two_sum(s1, *c3p, &s2);
    }
  } else {
    s0 = qd_quick_two_sum(s0, *c2p, &s1);
    if (s1 != 0.0) {
      s1 = qd_quick_two_sum(s1, *c3p, &s2);
    } else {
      s0 = qd_quick_two_sum(s0, *c3p, &s1);
    }
  }

  *c0p = s0;
  *c1p = s1;
  *c2p = s2;
}

HEADER_INLINE void tdr_renorm(td_real* xp) {
  qd_renorm3(&(xp->x[0]), &(xp->x[1]), &(xp->x[2]));
}

HEADER_INLINE void qd_three_sum(double* ap, double* bp, double* cp) {
  double t2;
  double t3;
  double t1 = qd_two_sum(*ap, *bp, &t2);
  *ap = qd_two_sum(*cp, t1, &t3);
  *bp = qd_two_sum(t2, t3, cp);
}

HEADER_INLINE void qd_three_sum2(double* ap, double* bp, double c) {
  double t2;
  double t3;
  const double t1 = qd_two_sum(*ap, *bp, &t2);
  *ap = qd_two_sum(c, t1, &t3);
  *bp = t2 + t3;
}

HEADER_CINLINE td_real tdr_make1(const double a) {
  const td_real retval = {{a, 0.0, 0.0}};
  return retval;
}

HEADER_CINLINE td_real tdr_makeu64(const uint64_t a) {
  const double hi = S_CAST(double, a);
  const td_real retval = {{hi, S_CAST(double, S_CAST(int64_t, a - S_CAST(uint64_t, hi))), 0.0}};
  return retval;
}

HEADER_CINLINE td_real tdr_make2(const double a, const double b) {
  const td_real retval = {{a, b, 0.0}};
  return retval;
}

HEADER_CINLINE td_real tdr_make_dd(const dd_real a) {
  const td_real retval = {{a.x[0], a.x[1], 0.0}};
  return retval;
}

HEADER_CINLINE td_real tdr_make(double a, double b, double c) {
  const td_real retval = {{a, b, c}};
  return retval;
}

HEADER_INLINE td_real tdr_addd(const td_real a, double b) {
  double e;
  double c0 = qd_two_sum(a.x[0], b, &e);
  double c1 = qd_two_sum(a.x[1], e, &e);
  double c2 = qd_two_sum(a.x[2], e, &e);

  qd_renorm4(&c0, &c1, &c2, &e);

  return tdr_make(c0, c1, c2);
}

HEADER_INLINE td_real tdr_add_dd(const td_real a, const dd_real b) {
  double t0;
  double t1;

  double s0 = qd_two_sum(a.x[0], b.x[0], &t0);
  double s1 = qd_two_sum(a.x[1], b.x[1], &t1);

  s1 = qd_two_sum(s1, t0, &t0);

  double s2 = a.x[2];
  qd_three_sum2(&s2, &t0, t1);

  qd_renorm4(&s0, &s1, &s2, &t0);
  return tdr_make(s0, s1, s2);
}

// s = qd_quick_three_accum(a, b, c) adds c to the dd-pair (a, b).
// If the result does not fit in two doubles, then the sum is
// output into s and (a,b) contains the remainder.  Otherwise
// s is zero and (a,b) contains the sum.
HEADER_INLINE double qd_quick_three_accum(double* ap, double* bp, double c) {
  double s = qd_two_sum(*bp, c, bp);
  s = qd_two_sum(*ap, s, ap);

  const uint32_t za = (*ap != 0.0);
  const uint32_t zb = (*bp != 0.0);

  if (za && zb) {
    return s;
  }

  if (!zb) {
    *bp = *ap;
    *ap = s;
  } else {
    *ap = s;
  }

  return 0.0;
}

// todo: check if sloppy-add is good enough for all our use cases

HEADER_INLINE td_real tdr_ieee_add(const td_real a, const td_real b) {
  int32_t i = 0;
  int32_t j = 0;
  int32_t k = 0;
  // (u,v) = double-length accumulator
  double u;
  double v;
  if (fabs(a.x[i]) > fabs(b.x[j])) {
    u = a.x[i++];
  } else {
    u = b.x[j++];
  }
  if (fabs(a.x[i]) > fabs(b.x[j])) {
    v = a.x[i++];
  } else {
    v = b.x[j++];
  }

  u = qd_quick_two_sum(u, v, &v);

  double x[3];
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  do {
    if ((i >= 3) && (j >= 3)) {
      x[k] = u;
      if (k < 2) {
        x[++k] = v;
      }
      break;
    }

    double t;
    if (i >= 3) {
      t = b.x[j++];
    } else if (j >= 3) {
      t = a.x[i++];
    } else if (fabs(a.x[i]) > fabs(b.x[j])) {
      t = a.x[i++];
    } else {
      t = b.x[j++];
    }

    const double s = qd_quick_three_accum(&u, &v, t);

    if (s != 0.0) {
      x[k++] = s;
    }
  } while (k < 3);

  // add the rest.
  for (k = i; k < 3; k++) {
    x[2] += a.x[k];
  }
  for (k = j; k < 3; k++) {
    x[2] += b.x[k];
  }

  qd_renorm3(&(x[0]), &(x[1]), &(x[2]));
  return tdr_make(x[0], x[1], x[2]);
}

// tdr_sloppy_add() deliberately omitted for now

HEADER_CINLINE td_real tdr_negate(const td_real a) {
  return tdr_make(-a.x[0], -a.x[1], -a.x[2]);
}

HEADER_INLINE td_real tdr_subd(const td_real a, double b) {
  return tdr_addd(a, -b);
}

HEADER_INLINE td_real tdr_sub_dd(const td_real a, const dd_real b) {
  return tdr_add_dd(a, ddr_negate(b));
}

HEADER_INLINE td_real tdr_ldexp(const td_real a, int32_t expi) {
  return tdr_make(ldexp(a.x[0], expi), ldexp(a.x[1], expi), ldexp(a.x[2], expi));
}

HEADER_CINLINE td_real tdr_mul_pwr2(const td_real a, double b) {
  return tdr_make(a.x[0] * b, a.x[1] * b, a.x[2] * b);
}

HEADER_INLINE td_real tdr_muld(const td_real a, double b) {
  // simplified qdr_muld() that effectively assumes the original a.x[2] is 0;
  // the original a.x[3] becomes a.x[2] here
  double q0;
  double q1;
  double s2;
  const double p0 = qd_two_prod(a.x[0], b, &q0);
  const double p1 = qd_two_prod(a.x[1], b, &q1);
  // p2 = 0, q2 = 0
  const double p3 = a.x[2] * b;

  double s0 = p0;

  double s1 = qd_two_sum(q0, p1, &s2);

  s2 = qd_two_sum(s2, q1, &q1);

  double s3 = q1 + p3;

  qd_renorm4(&s0, &s1, &s2, &s3);
  return tdr_make(s0, s1, s2);
}

// todo: check if sloppy-mul is good enough for all our use cases

td_real tdr_accurate_mul(const td_real a, const td_real b);

td_real tdr_divd(const td_real a, double b);

td_real tdr_accurate_div(const td_real a, const td_real b);

td_real tdr_sqr(const td_real a);

HEADER_CINLINE int32_t tdr_is_zero(const td_real a) {
  return (a.x[0] == 0.0);
}

HEADER_CINLINE int32_t tdr_is(const td_real a, double b) {
  // a.x[1] == 0.0 should guarantee a.x[2] == 0.0.
  return (a.x[0] == b) && (a.x[1] == 0.0);
}

HEADER_CINLINE int32_t tdr_ltd(const td_real a, const double b) {
  return (a.x[0] < b) || ((a.x[0] == b) && (a.x[1] < 0));
}

HEADER_CINLINE int32_t tdr_lt(const td_real a, const td_real b) {
  if (a.x[0] != b.x[0]) {
    return (a.x[0] < b.x[0]);
  }
  if (a.x[1] != b.x[1]) {
    return (a.x[1] < b.x[1]);
  }
  return (a.x[2] < b.x[2]);
}

HEADER_CINLINE int32_t tdr_leqd(const td_real a, const double b) {
  return (a.x[0] < b) || ((a.x[0] == b) && (a.x[1] <= 0));
}

HEADER_CINLINE int32_t tdr_leq(const td_real a, const td_real b) {
  if (a.x[0] != b.x[0]) {
    return (a.x[0] < b.x[0]);
  }
  if (a.x[1] != b.x[1]) {
    return (a.x[1] < b.x[1]);
  }
  return (a.x[2] <= b.x[2]);
}

HEADER_CINLINE int32_t tdr_geqd(const td_real a, const double b) {
  return (a.x[0] > b) || ((a.x[0] == b) && (a.x[1] >= 0));
}

HEADER_CINLINE int32_t tdr_geq(const td_real a, const td_real b) {
  if (a.x[0] != b.x[0]) {
    return (a.x[0] > b.x[0]);
  }
  if (a.x[1] != b.x[1]) {
    return (a.x[1] > b.x[1]);
  }
  return (a.x[2] >= b.x[2]);
}

HEADER_CINLINE int32_t tdr_gtd(const td_real a, const double b) {
  return (a.x[0] > b) || ((a.x[0] == b) && (a.x[1] > 0));
}

HEADER_CINLINE int32_t tdr_gt(const td_real a, const td_real b) {
  if (a.x[0] != b.x[0]) {
    return (a.x[0] > b.x[0]);
  }
  if (a.x[1] != b.x[1]) {
    return (a.x[1] > b.x[1]);
  }
  return (a.x[2] > b.x[2]);
}


HEADER_INLINE td_real tdr_add(const td_real a, const td_real b) {
  return tdr_ieee_add(a, b);
}

HEADER_INLINE td_real tdr_sub(const td_real a, const td_real b) {
  return tdr_add(a, tdr_negate(b));
}

HEADER_INLINE td_real tdr_mul(const td_real a, const td_real b) {
  return tdr_accurate_mul(a, b);
}

td_real tdr_exp(const td_real a);

td_real tdr_log(const td_real a);

td_real tdr_log1p(const td_real a);

td_real tdr_lfact(double xx);

HEADER_INLINE td_real tdr_add_lfacts(const double a, const double b) {
  return tdr_add(tdr_lfact(a), tdr_lfact(b));
}

td_real tdr_sort_and_add(uint32_t ct, td_real* tdrs);

HEADER_INLINE td_real tdr_sort_and_add3(td_real a_tdr, td_real b_tdr, td_real c_tdr) {
  td_real args[3];
  args[0] = a_tdr;
  args[1] = b_tdr;
  args[2] = c_tdr;
  return tdr_sort_and_add(3, args);
}


HEADER_INLINE dd_real ddr_log_2arg(const dd_real p, const dd_real q) {
  if (q.x[0] < 0.03125) {
    return ddr_log1p(ddr_negate(q));
  }
  return ddr_log(p);
}


// Preconditions:
// - numer_factorial_args[] and denom_factorial_args[] are
//   not-necessarily-sorted lists of length ffac_ct, describing a quotient of
//   factorial-products.  (If one list is longer than the other, just pad the
//   other with zeroes.)  All entries < 2^52.
// - odds_ratio^odds_ratio_pow is an exponential term to multiply the quotient
//   by at the end.
// - starting_lnprobv_tdr can either be fully initialized to
//     log(numer_odds_ratio_pow / (numer_factorial_args[0]! ...
//                                 numer_factorial_args[ffac_ct-1]!)),
//   have the first two terms initialized and x[2] set to DBL_MAX (to indicate
//   that it has only been computed to dd_real precision), or have x[0]
//   initialized to DBL_MAX to indicate that the calculation has not happened
//   at all.  In the latter two cases, it may be updated to higher accuracy.
// - Similarly, ln_odds_ratio_tdr can either be fully initialized to
//   log(odds_ratio), or it can have x[2] or x[0] initialized to DBL_MAX.
//
// Postconditions on success:
// - Return value is positive if the fraction > 1, negative if the fraction <
//   1, and zero if it's exactly 1.
// - *dbl_ptr is the double representation of the fraction, error limited to
//   1-2 ULPs.
// - numer_factorial_args[] and denom_factorial_args[] are sorted in
//   nondecreasing order.
//
// This could take a precomputed tdr_lfact table as an additional pair of
// parameters, but I don't think that makes much of a difference.
intptr_t CompareFactorialProducts(uint32_t ffac_ct, td_real odds_ratio_tdr, int64_t odds_ratio_pow, int64_t numer_odds_ratio_pow, uint64_t* numer_factorial_args, uint64_t* denom_factorial_args, td_real* starting_lnprobv_tdr_ptr, td_real* ln_odds_ratio_tdr_ptr, double* dbl_ptr);

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_HIGHPREC_H__
