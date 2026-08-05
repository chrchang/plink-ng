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

#include "plink2_float.h"

#ifdef __APPLE__
#  include <fenv.h>
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

void flush_denormals() {
#if defined(__x86_64__) || defined(__i386__)
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#else
#  ifdef __APPLE__
  fesetenv(FE_DFL_DISABLE_DENORMS_ENV);
#  else
  // generic arm64
  uint64_t fpcr;
  asm volatile("mrs %0, fpcr" : "=r"(fpcr));
  fpcr |= (1LLU << 24);
  asm volatile("msr fpcr, %0" : : "r"(fpcr));
#  endif
#endif
}

double poly_eval(const double* coefs, uint32_t degree, double xx) {
  // possible todo: benchmark this threshold
  if (degree < 3) {
    double total = coefs[degree];
    while (degree) {
      total = prefer_fma(total, xx, coefs[--degree]);
    }
    return total;
  }
  // this should be a good balance between function size and asymptotic speed?
  const double x2 = xx * xx;
  double even_terms = prefer_fma(coefs[degree], x2, coefs[degree - 2]);
  double odd_terms = prefer_fma(coefs[degree - 1], x2, coefs[degree - 3]);
  degree -= 3;
  while (degree >= 2) {
    degree -= 2;
    even_terms = prefer_fma(even_terms, x2, coefs[degree + 1]);
    odd_terms = prefer_fma(odd_terms, x2, coefs[degree]);
  }
  if (degree == 1) {
    return prefer_fma(even_terms, x2, coefs[0]) + odd_terms * xx;
  }
  return prefer_fma(even_terms, xx, odd_terms);
}

double ratfun_eval_smallx(const double* numer_coefs, const double* denom_coefs, uint32_t degree, double xx) {
  if (degree < 3) {
    double numer = numer_coefs[degree];
    double denom = denom_coefs[degree];
    while (degree) {
      --degree;
      numer = prefer_fma(numer, xx, numer_coefs[degree]);
      denom = prefer_fma(denom, xx, denom_coefs[degree]);
    }
    return numer / denom;
  }
  const double x2 = xx * xx;
  double numer_even_terms = prefer_fma(numer_coefs[degree], x2, numer_coefs[degree - 2]);
  double denom_even_terms = prefer_fma(denom_coefs[degree], x2, denom_coefs[degree - 2]);
  double numer_odd_terms = prefer_fma(numer_coefs[degree - 1], x2, numer_coefs[degree - 3]);
  double denom_odd_terms = prefer_fma(denom_coefs[degree - 1], x2, denom_coefs[degree - 3]);
  degree -= 3;
  while (degree >= 2) {
    --degree;
    numer_even_terms = prefer_fma(numer_even_terms, x2, numer_coefs[degree]);
    denom_even_terms = prefer_fma(denom_even_terms, x2, denom_coefs[degree]);
    --degree;
    numer_odd_terms = prefer_fma(numer_odd_terms, x2, numer_coefs[degree]);
    denom_odd_terms = prefer_fma(denom_odd_terms, x2, denom_coefs[degree]);
  }
  if (degree == 1) {
    return (prefer_fma(numer_even_terms, x2, numer_coefs[0]) + numer_odd_terms * xx) / (prefer_fma(denom_even_terms, x2, denom_coefs[0]) + denom_odd_terms * xx);
  }
  return prefer_fma(numer_even_terms, xx, numer_odd_terms) / prefer_fma(denom_even_terms, xx, denom_odd_terms);
}

double ratfun_eval_largex(const double* numer_coefs, const double* denom_coefs, uint32_t degree, double xx) {
  const double invx = 1.0 / xx;
  if (degree < 3) {
    double numer = numer_coefs[0];
    double denom = denom_coefs[0];
    for (uint32_t uii = 1; uii <= degree; ++uii) {
      numer = prefer_fma(numer, invx, numer_coefs[uii]);
      denom = prefer_fma(denom, invx, denom_coefs[uii]);
    }
    return numer / denom;
  }
  const double invx2 = invx * invx;
  double numer_even_terms = prefer_fma(numer_coefs[0], invx2, numer_coefs[2]);
  double denom_even_terms = prefer_fma(denom_coefs[0], invx2, denom_coefs[2]);
  double numer_odd_terms = prefer_fma(numer_coefs[1], invx2, numer_coefs[3]);
  double denom_odd_terms = prefer_fma(denom_coefs[1], invx2, denom_coefs[3]);
  for (uint32_t uii = 4; uii < degree; uii += 2) {
    numer_even_terms = prefer_fma(numer_even_terms, invx2, numer_coefs[uii]);
    denom_even_terms = prefer_fma(denom_even_terms, invx2, denom_coefs[uii]);
    numer_odd_terms = prefer_fma(numer_odd_terms, invx2, numer_coefs[uii + 1]);
    denom_odd_terms = prefer_fma(denom_odd_terms, invx2, denom_coefs[uii + 1]);
  }
  if (degree % 2 == 0) {
    return (prefer_fma(numer_even_terms, invx2, numer_coefs[degree]) + numer_odd_terms * invx) / (prefer_fma(denom_even_terms, invx2, denom_coefs[degree]) + denom_odd_terms * invx);
  }
  return prefer_fma(numer_even_terms, invx, numer_odd_terms) / prefer_fma(denom_even_terms, invx, denom_odd_terms);
}

#ifdef __cplusplus
}
#endif
