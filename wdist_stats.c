#include "wdist_common.h"
#include "ipmpar.h"
#include "dcdflib.h"

double chiprob_p(double xx, double df) {
  int st = 0;
  int ww = 1;
  double bnd = 1;
  double pp;
  double qq;
  cdfchi(&ww, &pp, &qq, &xx, &df, &st, &bnd);
  if (st) {
    return -9;
  }
  return qq;
}

double inverse_chiprob(double qq, double df) {
  double pp = 1 - qq;
  int32_t st = 0;
  int32_t ww = 2;
  double bnd = 1;
  double xx;

  if (qq >= 1.0) {
    return 0;
  }
  cdfchi(&ww, &pp, &qq, &xx, &df, &st, &bnd);
  if (st != 0) {
    return -9;
  }
  return xx;
}

// Inverse normal distribution

//
// Lower tail quantile for standard normal distribution function.
//
// This function returns an approximation of the inverse cumulative
// standard normal distribution function.  I.e., given P, it returns
// an approximation to the X satisfying P = Pr{Z <= X} where Z is a
// random variable from the standard normal distribution.
//
// The algorithm uses a minimax approximation by rational functions
// and the result has a relative error whose absolute value is less
// than 1.15e-9.
//
// Author:      Peter J. Acklam
// Time-stamp:  2002-06-09 18:45:44 +0200
// E-mail:      jacklam@math.uio.no
// WWW URL:     http://www.math.uio.no/~jacklam
//
// C implementation adapted from Peter's Perl version

// Coefficients in rational approximations.

static const double ivn_a[] =
  {
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  };

static const double ivn_b[] =
  {
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01
  };

static const double ivn_c[] =
  {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
     2.938163982698783e+00
  };

static const double ivn_d[] =
  {
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
  };

#define IVN_LOW 0.02425
#define IVN_HIGH 0.97575

double ltqnorm(double p) {
  // assumes 0 < p < 1
  double q, r;

  if (p < IVN_LOW) {
    // Rational approximation for lower region
    q = sqrt(-2*log(p));
    return (((((ivn_c[0]*q+ivn_c[1])*q+ivn_c[2])*q+ivn_c[3])*q+ivn_c[4])*q+ivn_c[5]) /
      ((((ivn_d[0]*q+ivn_d[1])*q+ivn_d[2])*q+ivn_d[3])*q+1);
  } else if (p > IVN_HIGH) {
    // Rational approximation for upper region
    q  = sqrt(-2*log(1-p));
    return -(((((ivn_c[0]*q+ivn_c[1])*q+ivn_c[2])*q+ivn_c[3])*q+ivn_c[4])*q+ivn_c[5]) /
      ((((ivn_d[0]*q+ivn_d[1])*q+ivn_d[2])*q+ivn_d[3])*q+1);
  } else {
    // Rational approximation for central region
    q = p - 0.5;
    r = q*q;
    return (((((ivn_a[0]*r+ivn_a[1])*r+ivn_a[2])*r+ivn_a[3])*r+ivn_a[4])*r+ivn_a[5])*q /
      (((((ivn_b[0]*r+ivn_b[1])*r+ivn_b[2])*r+ivn_b[3])*r+ivn_b[4])*r+1);
  }
}

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22) {
  // Basic 2x2 Fisher exact test p-value calculation.
  double tprob = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double cur_prob = tprob;
  double cprob = 0;
  uint32_t uii;
  double cur11;
  double cur12;
  double cur21;
  double cur22;
  double preaddp;
  // Ensure we are left of the distribution center, m11 <= m22, and m12 <= m21.
  if (m12 > m21) {
    uii = m12;
    m12 = m21;
    m21 = uii;
  }
  if (m11 > m22) {
    uii = m11;
    m11 = m22;
    m22 = uii;
  }
  if ((((uint64_t)m11) * m22) > (((uint64_t)m12) * m21)) {
    uii = m11;
    m11 = m12;
    m12 = uii;
    uii = m21;
    m21 = m22;
    m22 = uii;
  }
  cur11 = m11;
  cur12 = m12;
  cur21 = m21;
  cur22 = m22;
  while (cur12 > 0.5) {
    cur11 += 1;
    cur22 += 1;
    cur_prob *= (cur12 * cur21) / (cur11 * cur22);
    cur12 -= 1;
    cur21 -= 1;
    if (cur_prob == INFINITY) {
      return 0;
    }
    if (cur_prob < EXACT_TEST_BIAS) {
      tprob += cur_prob;
      break;
    }
    cprob += cur_prob;
  }
  if (cprob == 0) {
    return 1;
  }
  if (cur12 > 0.5) {
    do {
      cur11 += 1;
      cur22 += 1;
      cur_prob *= (cur12 * cur21) / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
      if (tprob <= preaddp) {
	break;
      }
    } while (cur12 > 0.5);
  }
  if (m11) {
    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;
    cur_prob = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
    do {
      cur12 += 1;
      cur21 += 1;
      cur_prob *= (cur11 * cur22) / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
      if (tprob <= preaddp) {
	return preaddp / (cprob + preaddp);
      }
    } while (cur11 > 0.5);
  }
  return tprob / (cprob + tprob);
}

double fisher22_tail_pval(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t right_offset, double tot_prob, double right_prob, double tail_sum, uint32_t new_m11) {
  // Given that the left (w.r.t. m11) reference contingency table has
  // likelihood 1/tot_prob, the contingency table with m11 increased by
  // right_offset has likelihood right_prob/tot_prob, and the tails (up to but
  // not including the two references) sum to tail_sum/tot_prob, this
  // calculates the p-value of the given m11 (which must be on one tail).
  double left_prob = 1.0;
  double dxx = new_m11;
  double left11 = m11;
  double left12 = m12;
  double left21 = m21;
  double left22 = m22;
  double right11 = m11 + right_offset;
  double right12 = m12 - right_offset;
  double right21 = m21 - right_offset;
  double right22 = m22 + right_offset;
  double psum;
  if (new_m11 < m11) {
    dxx += 0.5; // unnecessary (53 vs. 32 bits precision), but whatever
    do {
      left12 += 1;
      left21 += 1;
      left_prob *= left11 * left22 / (left12 * left21);
      if (left_prob == 0) {
	return 0;
      }
      left11 -= 1;
      left22 -= 1;
    } while (left11 > dxx);
    psum = left_prob;
    dxx = left_prob * (1 + SMALLISH_EPSILON);
    while (right_prob > dxx) {
      right11 += 1;
      right22 += 1;
      right_prob *= right12 * right21 / (right11 * right22);
      right12 -= 1;
      right21 -= 1;
    }
    psum += right_prob;
  } else {
    dxx -= 0.5;
    do {
      right11 += 1;
      right22 += 1;
      right_prob *= right12 * right21 / (right11 * right22);
      if (right_prob == 0) {
	return 0;
      }
      right12 -= 1;
      right21 -= 1;
    } while (right11 < dxx);
    psum = right_prob;
    dxx = right_prob * (1 + SMALLISH_EPSILON);
    while (left_prob > dxx) {
      left12 += 1;
      left21 += 1;
      left_prob *= left11 * left22 / (left12 * left21);
      left11 -= 1;
      left22 -= 1;
    }
    psum += left_prob;
  }
  if (left11 * left22 > 0.5) {
    do {
      left12 += 1;
      left21 += 1;
      left_prob *= left11 * left22 / (left12 * left21);
      left11 -= 1;
      left22 -= 1;
      dxx = psum;
      psum += left_prob;
    } while (psum > dxx);
  }
  if (right12 * right21 > 0.5) {
    do {
      right11 += 1;
      right22 += 1;
      right_prob *= right12 * right21 / (right11 * right22);
      right12 -= 1;
      right21 -= 1;
      dxx = psum;
      psum += right_prob;
    } while (psum > dxx);
  }
  return psum / tot_prob;
}

void fisher22_precomp_pval_bounds(double pval, uint32_t row1_sum, uint32_t col1_sum, uint32_t total, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiep, double* tot_probp, double* right_probp, double* tail_sump) {
  // Treating m11 as the only variable, this returns the minimum and (maximum -
  // 1) values of m11 which are less extreme than the observed result.  Sets
  // the low bit of *tiep if m11 = minimum - 1 results in a tie (2^{-35}
  // tolerance), and the second-lowest bit if m11 = maximum does.  Finally,
  // returns the values necessary for invoking fisher22_tail_pval() with
  //   m11 := (m11_min - (tie & 1)) and
  //   right_offset := (m11_max - 1) + (tie / 2) - m11.
  //
  // Algorithm:
  // 1. Determine center.
  // 2. Sum unscaled probabilities in both directions to precision limit.
  // 3. Proceed further outwards to (pval * unscaled_psum) precision limit,
  //    fill in the remaining return values.
  double tot_prob = 1.0 / EXACT_TEST_BIAS;
  double left_prob = tot_prob;
  double right_prob = tot_prob;
  intptr_t m11_offset = 0;
  double tail_prob = 0;
  double dxx;
  double left11;
  double left12;
  double left21;
  double left22;
  double right11;
  double right12;
  double right21;
  double right22;
  double cur_prob;
  double cur11;
  double cur12;
  double cur21;
  double cur22;
  double threshold;
  intptr_t lii;
  uint32_t uii;
  if (!total) {
    // hardcode this to avoid divide-by-zero
    *m11_minp = 0;
    *m11_maxp = 0;
    *tiep = 2;
    // no need to initialize the other return values, they should never be used
    return;
  } else {
    *tiep = 0;
    if (pval == 0) {
      if (total >= row1_sum + col1_sum) {
	*m11_minp = 0;
	*m11_maxp = MINV(row1_sum, col1_sum) + 1;
      } else {
	*m11_minp = row1_sum + col1_sum - total;
	*m11_maxp = total - MAXV(row1_sum, col1_sum) + 1;
      }
      return;
    }
  }
  // Center must be adjacent to the x which satisfies
  //   (m11 + x)(m22 + x) = (m12 - x)(m21 - x), so
  //   x = (m12 * m21 - m11 * m22) / (m11 + m12 + m21 + m22)
  if (total >= row1_sum + col1_sum) {
    // m11 = 0;
    // m12 = row1_sum;
    // m21 = col1_sum;
    // m22 = total - row1_sum - col1_sum;
    lii = (((uint64_t)row1_sum) * col1_sum) / total;
    left11 = lii;
    left12 = row1_sum - lii;
    left21 = col1_sum - lii;
    left22 = (total - row1_sum - col1_sum) + lii;
  } else {
    // m11 = row1_sum + col1_sum - total;
    // m12 = row1_sum - m11;
    // m21 = col1_sum - m11;
    // m22 = 0;
    lii = (((uint64_t)(total - row1_sum)) * (total - col1_sum)) / total;
    // Force m11 <= m22 for internal calculation, then adjust at end.
    m11_offset = row1_sum + col1_sum - total;
    left11 = lii;
    left12 = total - col1_sum - lii;
    left21 = total - row1_sum - lii;
    left22 = m11_offset + lii;
  }
  // We rounded x down.  Should we have rounded up instead?
  if ((left11 + 1) * (left22 + 1) < left12 * left21) {
    left11 += 1;
    left12 -= 1;
    left21 -= 1;
    left22 += 1;
  }
  // Safe to force m12 <= m21.
  if (left12 > left21) {
    dxx = left12;
    left12 = left21;
    left21 = dxx;
  }
  // Sum right side to limit, then left.
  right11 = left11;
  right12 = left12;
  right21 = left21;
  right22 = left22;
  do {
    if (right12 < 0.5) {
      break;
    }
    right11 += 1;
    right22 += 1;
    right_prob *= (right12 * right21) / (right11 * right22);
    right12 -= 1;
    right21 -= 1;
    dxx = tot_prob;
    tot_prob += right_prob;
  } while (tot_prob > dxx);
  do {
    if (left11 < 0.5) {
      break;
    }
    left12 += 1;
    left21 += 1;
    left_prob *= (left11 * left22) / (left12 * left21);
    left11 -= 1;
    left22 -= 1;
    dxx = tot_prob;
    tot_prob += left_prob;
  } while (tot_prob > dxx);
  // Now traverse the tails to p-value precision limit.
  // Upper bound for tail sum, if current element c is included:
  //   (c + cr + cr^2 + ...) + (c + cs + cs^2 + ...)
  // = c(1/(1 - r) + 1/(1 - s))
  // Compare this to pval * tot_prob.
  // I.e. compare c to pval * tot_prob * (1-r)(1-s) / (2-r-s)
  dxx = 1 - (left11 * left22) / ((left12 + 1) * (left21 + 1));
  threshold = 1 - (right12 * right21) / ((right11 + 1) * (right22 + 1));
  threshold = pval * tot_prob * dxx * threshold / (dxx + threshold);
  while (left11 > 0.5) {
    if (left_prob < threshold) {
      tail_prob = left_prob;
      cur11 = left11;
      cur12 = left12;
      cur21 = left21;
      cur22 = left22;
      cur_prob = left_prob;
      do {
	cur12 += 1;
	cur21 += 1;
	cur_prob *= (cur11 * cur22) / (cur12 * cur21);
	cur11 -= 1;
	cur22 -= 1;
	dxx = tail_prob;
	tail_prob += cur_prob;
      } while (dxx < tail_prob);
      left11 += 1;
      left22 += 1;
      left_prob *= (left12 * left21) / (left11 * left22);
      left12 -= 1;
      left21 -= 1;
      break;
    }
    left12 += 1;
    left21 += 1;
    left_prob *= (left11 * left22) / (left12 * left21);
    left11 -= 1;
    left22 -= 1;
  }
  while (right12 > 0.5) {
    if (right_prob < threshold) {
      tail_prob += right_prob;
      cur11 = right11;
      cur12 = right12;
      cur21 = right21;
      cur22 = right22;
      cur_prob = right_prob;
      do {
	cur11 += 1;
	cur22 += 1;
	cur_prob *= (cur12 * cur21) / (cur11 * cur22);
	cur12 -= 1;
	cur21 -= 1;
	dxx = tail_prob;
	tail_prob += cur_prob;
      } while (dxx < tail_prob);
      right12 += 1;
      right21 += 1;
      right_prob *= (right11 * right22) / (right12 * right21);
      right11 -= 1;
      right22 -= 1;
      break;
    }
    right11 += 1;
    right22 += 1;
    right_prob *= (right12 * right21) / (right11 * right22);
    right12 -= 1;
    right21 -= 1;
  }
  dxx = pval * tot_prob * (1 - SMALLISH_EPSILON);
  threshold = pval * tot_prob * (1 + SMALLISH_EPSILON);
  while (1) {
    if (left_prob < right_prob * (1 - SMALLISH_EPSILON)) {
      if (tail_prob + left_prob > threshold) {
	break;
      }
      tail_prob += left_prob;
      uii = 1;
    } else if (right_prob < left_prob * (1 - SMALLISH_EPSILON)) {
      if (tail_prob + right_prob > threshold) {
	break;
      }
      tail_prob += right_prob;
      uii = 2;
    } else {
      if (tail_prob + left_prob + right_prob > threshold) {
	break;
      }
      tail_prob += left_prob + right_prob;
      uii = 3;
    }
    if (tail_prob > dxx) {
      *tiep = uii;
      break;
    }
    // if more speed is necessary, we could use a buffer to save all unscaled
    // probabilities during the initial outward traversal.
    if (uii & 1) {
      left11 += 1;
      left22 += 1;
      left_prob *= (left12 * left21) / (left11 * left22);
      left12 -= 1;
      left21 -= 1;
    }
    if (uii & 2) {
      right12 += 1;
      right21 += 1;
      right_prob *= (right11 * right22) / (right12 * right21);
      right11 -= 1;
      right22 -= 1;
    }
  }
  *m11_minp = m11_offset + ((intptr_t)left11);
  *m11_maxp = m11_offset + ((intptr_t)right11) + 1;
  dxx = 1.0 / left_prob;
  *tot_probp = tot_prob * dxx;
  *right_probp = right_prob * dxx;
  *tail_sump = tail_prob * dxx;
}

int32_t fisher23_tailsum(double* base_probp, double* saved12p, double* saved13p, double* saved22p, double* saved23p, double *totalp, uint32_t right_side) {
  double total = 0;
  double cur_prob = *base_probp;
  double tmp12 = *saved12p;
  double tmp13 = *saved13p;
  double tmp22 = *saved22p;
  double tmp23 = *saved23p;
  double tmps12;
  double tmps13;
  double tmps22;
  double tmps23;
  double prev_prob;
  // identify beginning of tail
  if (right_side) {
    if (cur_prob > EXACT_TEST_BIAS) {
      prev_prob = tmp13 * tmp22;
      while (prev_prob > 0.5) {
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= prev_prob / (tmp12 * tmp23);
	tmp13 -= 1;
	tmp22 -= 1;
	if (cur_prob <= EXACT_TEST_BIAS) {
	  break;
	}
	prev_prob = tmp13 * tmp22;
      }
      *base_probp = cur_prob;
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
    } else {
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
      while (1) {
	prev_prob = cur_prob;
	tmp13 += 1;
	tmp22 += 1;
	cur_prob *= (tmp12 * tmp23) / (tmp13 * tmp22);
	if (cur_prob < prev_prob) {
	  return 1;
	}
	tmp12 -= 1;
	tmp23 -= 1;
	// throw in extra (1 - SMALL_EPSILON) multiplier to prevent rounding
	// errors from causing this to keep going when the left-side test
	// stopped
	if (cur_prob > (1 - SMALL_EPSILON) * EXACT_TEST_BIAS) {
	  break;
	}
	total += cur_prob;
      }
      prev_prob = cur_prob;
      cur_prob = *base_probp;
      *base_probp = prev_prob;
    }
  } else {
    if (cur_prob > EXACT_TEST_BIAS) {
      prev_prob = tmp12 * tmp23;
      while (prev_prob > 0.5) {
	tmp13 += 1;
	tmp22 += 1;
	cur_prob *= prev_prob / (tmp13 * tmp22);
	tmp12 -= 1;
	tmp23 -= 1;
	if (cur_prob <= EXACT_TEST_BIAS) {
	  break;
	}
	prev_prob = tmp12 * tmp23;
      }
      *base_probp = cur_prob;
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
    } else {
      tmps12 = tmp12;
      tmps13 = tmp13;
      tmps22 = tmp22;
      tmps23 = tmp23;
      while (1) {
	prev_prob = cur_prob;
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= (tmp13 * tmp22) / (tmp12 * tmp23);
	if (cur_prob < prev_prob) {
	  return 1;
	}
	tmp13 -= 1;
	tmp22 -= 1;
	if (cur_prob > EXACT_TEST_BIAS) {
	  break;
	}
	total += cur_prob;
      }
      prev_prob = cur_prob;
      cur_prob = *base_probp;
      *base_probp = prev_prob;
    }
  }
  *saved12p = tmp12;
  *saved13p = tmp13;
  *saved22p = tmp22;
  *saved23p = tmp23;
  if (cur_prob > EXACT_TEST_BIAS) {
    *totalp = 0;
    return 0;
  }
  // sum tail to floating point precision limit
  if (right_side) {
    prev_prob = total;
    total += cur_prob;
    while (total > prev_prob) {
      tmps12 += 1;
      tmps23 += 1;
      cur_prob *= (tmps13 * tmps22) / (tmps12 * tmps23);
      tmps13 -= 1;
      tmps22 -= 1;
      prev_prob = total;
      total += cur_prob;
    }
  } else {
    prev_prob = total;
    total += cur_prob;
    while (total > prev_prob) {
      tmps13 += 1;
      tmps22 += 1;
      cur_prob *= (tmps12 * tmps23) / (tmps13 * tmps22);
      tmps12 -= 1;
      tmps23 -= 1;
      prev_prob = total;
      total += cur_prob;
    }
  }
  *totalp = total;
  return 0;
}

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23) {
  // 2x3 Fisher-Freeman-Halton exact test p-value calculation.
  // The number of tables involved here is still small enough that the network
  // algorithm (and the improved variants thereof that I've seen) are
  // suboptimal; a 2-dimensional version of the SNPHWE2 strategy has higher
  // performance.
  // 2x4, 2x5, and 3x3 should also be practical with this method, but beyond
  // that I doubt it's worth the trouble.
  // Complexity of approach is O(n^{df/2}), where n is number of observations.
  double cur_prob = (1 - SMALLISH_EPSILON) * EXACT_TEST_BIAS;
  double tprob = cur_prob;
  double cprob = 0;
  uint32_t dir = 0; // 0 = forwards, 1 = backwards
  double dyy = 0;
  double base_probl;
  double base_probr;
  double orig_base_probl;
  double orig_base_probr;
  double orig_row_prob;
  double row_prob;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double cur11;
  double cur21;
  double savedl12;
  double savedl13;
  double savedl22;
  double savedl23;
  double savedr12;
  double savedr13;
  double savedr22;
  double savedr23;
  double orig_savedl12;
  double orig_savedl13;
  double orig_savedl22;
  double orig_savedl23;
  double orig_savedr12;
  double orig_savedr13;
  double orig_savedr22;
  double orig_savedr23;
  double tmp12;
  double tmp13;
  double tmp22;
  double tmp23;
  double dxx;
  double preaddp;
  // Ensure m11 + m21 <= m12 + m22 <= m13 + m23.
  uii = m11 + m21;
  ujj = m12 + m22;
  if (uii > ujj) {
    ukk = m11;
    m11 = m12;
    m12 = ukk;
    ukk = m21;
    m21 = m22;
    m22 = ukk;
    ukk = uii;
    uii = ujj;
    ujj = ukk;
  }
  ukk = m13 + m23;
  if (ujj > ukk) {
    ujj = ukk;
    ukk = m12;
    m12 = m13;
    m13 = ukk;
    ukk = m22;
    m22 = m23;
    m23 = ukk;
  }
  if (uii > ujj) {
    ukk = m11;
    m11 = m12;
    m12 = ukk;
    ukk = m21;
    m21 = m22;
    m22 = ukk;
  }
  // Ensure majority of probability mass is in front of m11.
  if ((((uint64_t)m11) * (m22 + m23)) > (((uint64_t)m21) * (m12 + m13))) {
    ukk = m11;
    m11 = m21;
    m21 = ukk;
    ukk = m12;
    m12 = m22;
    m22 = ukk;
    ukk = m13;
    m13 = m23;
    m23 = ukk;
  }
  if ((((uint64_t)m12) * m23) > (((uint64_t)m13) * m22)) {
    base_probr = cur_prob;
    savedr12 = m12;
    savedr13 = m13;
    savedr22 = m22;
    savedr23 = m23;
    tmp12 = savedr12;
    tmp13 = savedr13;
    tmp22 = savedr22;
    tmp23 = savedr23;
    // m12 and m23 must be nonzero
    dxx = tmp12 * tmp23;
    do {
      tmp13 += 1;
      tmp22 += 1;
      cur_prob *= dxx / (tmp13 * tmp22);
      tmp12 -= 1;
      tmp23 -= 1;
      if (cur_prob <= EXACT_TEST_BIAS) {
	tprob += cur_prob;
	break;
      }
      cprob += cur_prob;
      if (cprob == INFINITY) {
	return 0;
      }
      dxx = tmp12 * tmp23;
      // must enforce tmp12 >= 0 and tmp23 >= 0 since we're saving these
    } while (dxx > 0.5);
    savedl12 = tmp12;
    savedl13 = tmp13;
    savedl22 = tmp22;
    savedl23 = tmp23;
    base_probl = cur_prob;
    do {
      tmp13 += 1;
      tmp22 += 1;
      cur_prob *= (tmp12 * tmp23) / (tmp13 * tmp22);
      tmp12 -= 1;
      tmp23 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
    } while (tprob > preaddp);
    tmp12 = savedr12;
    tmp13 = savedr13;
    tmp22 = savedr22;
    tmp23 = savedr23;
    cur_prob = base_probr;
    do {
      tmp12 += 1;
      tmp23 += 1;
      cur_prob *= (tmp13 * tmp22) / (tmp12 * tmp23);
      tmp13 -= 1;
      tmp22 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
    } while (tprob > preaddp);
  } else {
    base_probl = cur_prob;
    savedl12 = m12;
    savedl13 = m13;
    savedl22 = m22;
    savedl23 = m23;
    if (!((((uint64_t)m12) * m23) + (((uint64_t)m13) * m22))) {
      base_probr = cur_prob;
      savedr12 = savedl12;
      savedr13 = savedl13;
      savedr22 = savedl22;
      savedr23 = savedl23;
    } else {
      tmp12 = savedl12;
      tmp13 = savedl13;
      tmp22 = savedl22;
      tmp23 = savedl23;
      dxx = tmp13 * tmp22;
      do {
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= dxx / (tmp12 * tmp23);
	tmp13 -= 1;
	tmp22 -= 1;
	if (cur_prob <= EXACT_TEST_BIAS) {
	  tprob += cur_prob;
	  break;
	}
	cprob += cur_prob;
	if (cprob == INFINITY) {
	  return 0;
	}
	dxx = tmp13 * tmp22;
      } while (dxx > 0.5);
      savedr12 = tmp12;
      savedr13 = tmp13;
      savedr22 = tmp22;
      savedr23 = tmp23;
      base_probr = cur_prob;
      do {
	tmp12 += 1;
	tmp23 += 1;
	cur_prob *= (tmp13 * tmp22) / (tmp12 * tmp23);
	tmp13 -= 1;
	tmp22 -= 1;
	preaddp = tprob;
	tprob += cur_prob;
      } while (tprob > preaddp);
      tmp12 = savedl12;
      tmp13 = savedl13;
      tmp22 = savedl22;
      tmp23 = savedl23;
      cur_prob = base_probl;
      do {
	tmp13 += 1;
	tmp22 += 1;
	cur_prob *= (tmp12 * tmp23) / (tmp13 * tmp22);
	tmp12 -= 1;
	tmp23 -= 1;
	preaddp = tprob;
	tprob += cur_prob;
      } while (tprob > preaddp);
    }
  }
  row_prob = tprob + cprob;
  orig_base_probl = base_probl;
  orig_base_probr = base_probr;
  orig_row_prob = row_prob;
  orig_savedl12 = savedl12;
  orig_savedl13 = savedl13;
  orig_savedl22 = savedl22;
  orig_savedl23 = savedl23;
  orig_savedr12 = savedr12;
  orig_savedr13 = savedr13;
  orig_savedr22 = savedr22;
  orig_savedr23 = savedr23;
  for (; dir < 2; dir++) {
    cur11 = m11;
    cur21 = m21;
    if (dir) {
      base_probl = orig_base_probl;
      base_probr = orig_base_probr;
      row_prob = orig_row_prob;
      savedl12 = orig_savedl12;
      savedl13 = orig_savedl13;
      savedl22 = orig_savedl22;
      savedl23 = orig_savedl23;
      savedr12 = orig_savedr12;
      savedr13 = orig_savedr13;
      savedr22 = orig_savedr22;
      savedr23 = orig_savedr23;
      ukk = m11;
      if (ukk > m22 + m23) {
	ukk = m22 + m23;
      }
    } else {
      ukk = m21;
      if (ukk > m12 + m13) {
	ukk = m12 + m13;
      }
    }
    ukk++;
    while (--ukk) {
      if (dir) {
	cur21 += 1;
	if (savedl23) {
	  savedl13 += 1;
	  row_prob *= (cur11 * (savedl22 + savedl23)) / (cur21 * (savedl12 + savedl13));
	  base_probl *= (cur11 * savedl23) / (cur21 * savedl13);
	  savedl23 -= 1;
	} else {
	  savedl12 += 1;
	  row_prob *= (cur11 * (savedl22 + savedl23)) / (cur21 * (savedl12 + savedl13));
	  base_probl *= (cur11 * savedl22) / (cur21 * savedl12);
	  savedl22 -= 1;
	}
	cur11 -= 1;
      } else {
	cur11 += 1;
	if (savedl12) {
	  savedl22 += 1;
	  row_prob *= (cur21 * (savedl12 + savedl13)) / (cur11 * (savedl22 + savedl23));
	  base_probl *= (cur21 * savedl12) / (cur11 * savedl22);
	  savedl12 -= 1;
	} else {
	  savedl23 += 1;
	  row_prob *= (cur21 * (savedl12 + savedl13)) / (cur11 * (savedl22 + savedl23));
	  base_probl *= (cur21 * savedl13) / (cur11 * savedl23);
	  savedl13 -= 1;
	}
	cur21 -= 1;
      }
      if (fisher23_tailsum(&base_probl, &savedl12, &savedl13, &savedl22, &savedl23, &dxx, 0)) {
	break;
      }
      tprob += dxx;
      if (dir) {
	if (savedr22) {
	  savedr12 += 1;
	  base_probr *= ((cur11 + 1) * savedr22) / (cur21 * savedr12);
	  savedr22 -= 1;
	} else {
	  savedr13 += 1;
	  base_probr *= ((cur11 + 1) * savedr23) / (cur21 * savedr13);
	  savedr23 -= 1;
	}
      } else {
	if (savedr13) {
	  savedr23 += 1;
	  base_probr *= ((cur21 + 1) * savedr13) / (cur11 * savedr23);
	  savedr13 -= 1;
	} else {
	  savedr22 += 1;
	  base_probr *= ((cur21 + 1) * savedr12) / (cur11 * savedr22);
	  savedr12 -= 1;
	}
      }
      fisher23_tailsum(&base_probr, &savedr12, &savedr13, &savedr22, &savedr23, &dyy, 1);
      tprob += dyy;
      cprob += row_prob - dxx - dyy;
      if (cprob == INFINITY) {
	return 0;
      }
    }
    if (!ukk) {
      continue;
    }
    savedl12 += savedl13;
    savedl22 += savedl23;
    if (dir) {
      while (1) {
	preaddp = tprob;
	tprob += row_prob;
	if (tprob <= preaddp) {
	  break;
	}
	cur21 += 1;
	savedl12 += 1;
	row_prob *= (cur11 * savedl22) / (cur21 * savedl12);
	cur11 -= 1;
	savedl22 -= 1;
      }
    } else {
      while (1) {
	preaddp = tprob;
	tprob += row_prob;
	if (tprob <= preaddp) {
	  break;
	}
	cur11 += 1;
	savedl22 += 1;
	row_prob *= (cur21 * savedl12) / (cur11 * savedl22);
	cur21 -= 1;
	savedl12 -= 1;
      }
    }
  }
  return tprob / (tprob + cprob);
}

void chi22_precomp_coeffs(intptr_t row1_sum, intptr_t col1_sum, intptr_t total, double* expm11p, double* recip_sump) {
  // chisq = (m11 - expm11)^2 * recip_sum
  // (see discussion for chi22_precomp_val_bounds() below.)
  //
  // expm11 = row1_sum * col1_sum / total
  // expm12 = row1_sum * col2_sum / total, etc.
  // recip_sum = 1 / expm11 + 1 / expm12 + 1 / expm21 + 1 / expm22
  // = total * (1 / (row1_sum * col1_sum) + 1 / (row1_sum * col2_sum) +
  //            1 / (row2_sum * col1_sum) + 1 / (row2_sum * col2_sum))
  // = total^3 / (row1_sum * col1_sum * row2_sum * col2_sum)
  intptr_t row2_sum = total - row1_sum;
  intptr_t col2_sum = total - col1_sum;
  double dxx;
  if (row1_sum && row2_sum && col1_sum && col2_sum) {
    dxx = total;
    *expm11p = ((double)(((int64_t)row1_sum) * col1_sum)) / ((double)total);
    *recip_sump = dxx * dxx * dxx / (((double)row1_sum) * ((double)row2_sum) * ((double)col1_sum) * ((double)col2_sum));
  } else {
    // since an entire row or column is zero, either m11 or m22 is zero
    // row1_sum + col1_sum - total = m11 - m22
    if (row1_sum + col1_sum < total) {
      *expm11p = 0;
    } else {
      *expm11p = row1_sum + col1_sum - total;
    }
    *recip_sump = 0;
  }
}

void chi22_precomp_val_bounds(double chisq, intptr_t row1_sum, intptr_t col1_sum, intptr_t total, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiep) {
  // Treating m11 as the only variable, this returns the minimum and (maximum -
  // 1) values of m11 which produce smaller chisq statistics than given.  Also
  // sets the low bit of *tiep if m11 = minimum - 1 results in a tie (2^{-21}
  // tolerance), and the second-lowest bit if m11 = maximum does.
  double expm11;
  double recip_sum;
  double cur11;
  double dxx;
  intptr_t ceil11;
  intptr_t lii;
  *tiep = 0;
  chi22_precomp_coeffs(row1_sum, col1_sum, total, &expm11, &recip_sum);
  if (recip_sum == 0) {
    // sum-0 row or column, no freedom at all
    *m11_minp = expm11;
    *m11_maxp = expm11;
    if (chisq == 0) {
      *tiep = 2;
    }
    return;
  }

  // double cur_stat = (cur11 - exp11) * (cur11 - exp11) * recipx11 + (cur12 - exp12) * (cur12 - exp12) * recipx12 + (cur21 - exp21) * (cur21 - exp21) * recipx21 + (cur22 - exp22) * (cur22 - exp22) * recipx22;
  // However, we have
  //   cur11 - exp11 = -(cur12 - exp12) = -(cur21 - exp21) = cur22 - exp22
  // So the chisq statistic reduces to
  //   (cur11 - exp11)^2 * (recipx11 + recipx12 + recipx21 + recipx22).

  ceil11 = MINV(row1_sum, col1_sum);
  // chisq = (cur11 - expm11)^2 * recip_sum
  // -> expm11 +/- sqrt(chisq / recip_sum) = cur11
  recip_sum = sqrt(chisq / recip_sum);
  cur11 = expm11 - recip_sum;
  dxx = cur11 + 1 - BIG_EPSILON;
  if (dxx < 0) {
    *m11_minp = 0;
  } else {
    lii = (intptr_t)dxx;
    *m11_minp = lii;
    if (lii == (intptr_t)(cur11 + BIG_EPSILON)) {
      *tiep = 1;
    }
  }
  cur11 = expm11 + recip_sum;
  if (cur11 > ceil11 + BIG_EPSILON) {
    *m11_maxp = ceil11 + 1;
  } else {
    dxx = cur11 + 1 - BIG_EPSILON;
    lii = (intptr_t)dxx;
    *m11_maxp = lii;
    if (lii == (intptr_t)(cur11 + BIG_EPSILON)) {
      *tiep |= 2;
    }
  }
}
