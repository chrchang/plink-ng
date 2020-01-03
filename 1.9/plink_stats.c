// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

#include "plink_stats.h"
#include "ipmpar.h"
#include "dcdflib.h"

// 2^{-40} for now, since 2^{-44} was too small on real data
#define FISHER_EPSILON 0.0000000000009094947017729282379150390625

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

double calc_tprob(double tt, double df) {
  int32_t st = 0;
  int32_t ww = 1;
  double bnd = 1;
  double pp;
  double qq;
  if (!realnum(tt)) {
    return -9;
  }
  tt = fabs(tt);
  cdft(&ww, &pp, &qq, &tt, &df, &st, &bnd);
  if (st != 0) {
    return -9;
  }
  return 2 * qq;
}

double inverse_tprob(double dbl_qq, double df) {
  double qq = dbl_qq * 0.5;
  double pp = 1 - qq;
  int32_t st = 0;
  int32_t ww = 2;
  double bnd = 1;
  double tt;
  cdft(&ww, &pp, &qq, &tt, &df, &st, &bnd);
  if (st != 0) {
    return -9;
  }
  return tt;
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

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp) {
  // This function implements an exact SNP test of Hardy-Weinberg
  // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
  // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
  // Equilibrium. American Journal of Human Genetics. 76: 000 - 000.
  //
  // The original version was written by Jan Wigginton.
  //
  // This version was written by Christopher Chang.  It contains the following
  // improvements over the original SNPHWE():
  // - Proper handling of >64k genotypes.  Previously, there was a potential
  //   integer overflow.
  // - Detection and efficient handling of floating point overflow and
  //   underflow.  E.g. instead of summing a tail all the way down, the loop
  //   stops once the latest increment underflows the partial sum's 53-bit
  //   precision; this results in a large speedup when max heterozygote count
  //   >1k.
  // - No malloc() call.  It's only necessary to keep track of a few partial
  //   sums.
  // - Support for the mid-p variant of this test.  See Graffelman J, Moreno V
  //   (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium.
  //
  // Note that the SNPHWE_t() function below is a lot more efficient for
  // testing against a p-value inclusion threshold.  SNPHWE2() should only be
  // used if you need the actual p-value.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  int64_t rare_copies = 2LL * obs_homr + obs_hets;
  int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  int32_t tie_ct = 1;
  double curr_hets_t2 = obs_hets;
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;
  double tailp = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;
  double preaddp;
  if (!genotypes2) {
    if (midp) {
      return 0.5;
    } else {
      return 1;
    }
  }

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    while (curr_hets_t2 > 1.5) {
      // het_probs[curr_hets] = 1
      // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
	  tie_ct++;
	}
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (curr_hets_t2 > 1.5) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_homr_t1 > 0.5) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      curr_hets_t1 += 2;
      curr_homr_t1 -= 1;
      curr_homc_t1 -= 1;
    }
  } else {
    // tail 1 = lower
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
          tie_ct++;
	}
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_hets_t1 > 1.5) {
      curr_homr_t1 += 1;
      curr_homc_t1 += 1;
      lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      curr_hets_t1 -= 2;
    }
  }
  if (!midp) {
    return tailp / (tailp + centerp);
  } else {
    return (tailp - ((1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (tailp + centerp);
  }
}

int32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh) {
  // Threshold-test-only version of SNPHWE2() which is usually able to exit
  // from the calculation earlier.  Returns 0 if these counts are close enough
  // to Hardy-Weinberg equilibrium, 1 otherwise.
  //
  // Suppose, for definiteness, that the number of observed hets is no less
  // than expectation.  (Same ideas apply for the other case.)  We proceed as
  // follows:
  // - Sum the *relative* likelihoods of more likely smaller het counts.
  // - Determine the minimum tail mass to pass the threshold.
  // - The majority of the time, the tail boundary elements are enough to pass
  //   the threshold; we never need to sum the remainder of the tails.
  // - And in the case of disequilibrium, we will often be able to immediately
  //   determine that the tail sum cannot possibly pass the threshold, just by
  //   looking at the tail boundary elements and using a geometric series to
  //   upper-bound the tail sums.
  // - Only when neither of these conditions hold do we start traveling down
  //   the tails.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  int64_t rare_copies = 2LL * obs_homr + obs_hets;
  int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  double curr_hets_t2 = obs_hets; // tail 2
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;

  // Subtract epsilon from initial probability mass, so that we can compare to
  // 1 when determining tail vs. center membership without floating point error
  // biting us in the ass
  double tailp1 = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp1;
  double tailp2 = 0;
  double tail1_ceil;
  double tail2_ceil;
  double lastp1;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;

  // Initially, if center sum reaches this, the test can immediately fail.
  // Once center is summed, this is recalculated, and when tail sum has reached
  // this, we've passed.
  double exit_thresh;
  double exit_threshx;
  double ratio;
  double preaddp;
  if (!genotypes2) {
    return 0;
  }

  // Convert thresh into reverse odds ratio.
  thresh = (1 - thresh) / thresh;

  // Expected het count:
  //   2 * rarefreq * (1 - rarefreq) * genotypes
  // = 2 * (rare_copies / (2 * genotypes)) * (1 - rarefreq) * genotypes
  // = rare_copies * (1 - (rare_copies / (2 * genotypes)))
  // = (rare_copies * (2 * genotypes - rare_copies)) / (2 * genotypes)
  //
  // The computational identity is
  //   P(nhets == n) := P(nhets == n+2) * (n+2) * (n+1) /
  //                    (4 * homr(n) * homc(n))
  // where homr() and homc() are the number of homozygous rares/commons needed
  // to maintain the same allele frequencies.
  // This probability is always decreasing when proceeding away from the
  // expected het count.

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    if (obs_hets < 2) {
      return 0;
    }

    // An initial upper bound on the tail sum is useful, since it lets us
    // report test failure before summing the entire center.  We use the
    // trivial bound of 1 + floor(rare_copies / 2): that's the total number
    // of possible het counts, and the relative probability for each count must
    // be <= 1 if it's in the tail.
    exit_thresh = (1 + (rare_copies / 2)) * thresh * EXACT_TEST_BIAS;

    // het_probs[curr_hets] = 1
    // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1) / (4 * (curr_homr + 1) * (curr_homc + 1))
    do {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp2 = lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp > exit_thresh) {
	return 1;
      }
    } while (curr_hets_t2 > 1.5);
    exit_thresh = centerp / thresh;
    if (tailp1 + tailp2 >= exit_thresh) {
      return 0;
    }
    // c + cr + cr^2 + ... = c/(1-r), which is an upper bound for the tail sum
    ratio = (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * (curr_homr_t2 + 1) * (curr_homc_t2 + 1));
    tail2_ceil = tailp2 / (1 - ratio);
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    // ratio for the other tail
    lastp1 = (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
    tail1_ceil = tailp1 / (1 - lastp1);
    if (tail1_ceil + tail2_ceil < exit_thresh) {
      return 1;
    }
    lastp1 *= tailp1;
    tailp1 += lastp1;

    if (obs_homr > 1) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      exit_threshx = exit_thresh - tailp2;
      do {
	curr_hets_t1 += 2;
	curr_homr_t1 -= 1;
	curr_homc_t1 -= 1;
	lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
	preaddp = tailp1;
	tailp1 += lastp1;
	if (tailp1 > exit_threshx) {
	  return 0;
	}
	if (tailp1 <= preaddp) {
	  break;
	}
      } while (curr_homr_t1 > 1.5);
    }
    if (tailp1 + tail2_ceil < exit_thresh) {
      return 1;
    }
    exit_threshx = exit_thresh - tailp1;
    while (curr_hets_t2 > 1) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= exit_threshx) {
	return 0;
      }
      if (tailp2 <= preaddp) {
	return 1;
      }
      curr_hets_t2 -= 2;
    }
    return 1;
  } else {
    // tail 1 = lower
    if (!obs_homr) {
      return 0;
    }
    exit_thresh = (1 + (rare_copies / 2)) * thresh * EXACT_TEST_BIAS;
    do {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp2 = lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp > exit_thresh) {
	return 1;
      }
    } while (curr_homr_t2 > 0.5);
    exit_thresh = centerp / thresh;
    if (tailp1 + tailp2 >= exit_thresh) {
      return 0;
    }
    ratio = (4 * curr_homr_t2 * curr_homc_t2) / ((curr_hets_t2 + 2) * (curr_hets_t2 + 1));
    tail2_ceil = tailp2 / (1 - ratio);
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr + 1;
    curr_homc_t1 = obs_homc + 1;
    lastp1 = (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
    tail1_ceil = tailp1 / (1 - lastp1);
    lastp1 *= tailp1;
    tailp1 += lastp1;

    if (tail1_ceil + tail2_ceil < exit_thresh) {
      return 1;
    }
    if (obs_hets >= 4) {
      exit_threshx = exit_thresh - tailp2;
      do {
	curr_hets_t1 -= 2;
	curr_homr_t1 += 1;
	curr_homc_t1 += 1;
	lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
	preaddp = tailp1;
        tailp1 += lastp1;
	if (tailp1 > exit_threshx) {
	  return 0;
	}
	if (tailp1 <= preaddp) {
	  break;
	}
      } while (curr_hets_t1 > 3.5);
    }
    if (tailp1 + tail2_ceil < exit_thresh) {
      return 1;
    }
    exit_threshx = exit_thresh - tailp1;
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= exit_threshx) {
	return 0;
      }
      if (tailp2 <= preaddp) {
	return 1;
      }
    }
    return 1;
  }
}

int32_t SNPHWE_midp_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh) {
  // Mid-p version of SNPHWE_t().  (There are enough fiddly differences that I
  // think it's better for this to be a separate function.)  Assumes threshold
  // is smaller than 0.5.
  intptr_t obs_homc;
  intptr_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  int64_t rare_copies = 2LL * obs_homr + obs_hets;
  int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  double curr_hets_t2 = obs_hets; // tail 2
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;
  double tailp1 = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5;
  double centerp = tailp1;
  double lastp2 = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double tailp2 = 0;
  double tail1_ceil;
  double tail2_ceil;
  double lastp1;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;
  double exit_thresh;
  double exit_threshx;
  double ratio;
  double preaddp;
  if (!genotypes2) {
    return 0;
  }
  thresh = (1 - thresh) / thresh;
  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    if (obs_hets < 2) {
      return 0;
    }
    exit_thresh = (1 + (rare_copies / 2)) * thresh * EXACT_TEST_BIAS;
    do {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
	  // tie with original contingency table, apply mid-p correction here
	  // too
          tailp2 = tailp1;
          centerp += tailp1;
	} else {
	  tailp2 = lastp2;
	}
	break;
      }
      centerp += lastp2;
      if (centerp > exit_thresh) {
	return 1;
      }
    } while (curr_hets_t2 > 1.5);
    exit_thresh = centerp / thresh;
    if (tailp1 + tailp2 >= exit_thresh) {
      return 0;
    }
    ratio = (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * (curr_homr_t2 + 1) * (curr_homc_t2 + 1));
    // this needs to work in both the tie and no-tie cases
    tail2_ceil = tailp2 + lastp2 * ratio / (1 - ratio);
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    lastp1 = (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
    // always a tie here
    tail1_ceil = tailp1 * 2 / (1 - lastp1) - tailp1;
    if (tail1_ceil + tail2_ceil < exit_thresh) {
      return 1;
    }
    lastp1 *= tailp1 * 2;
    tailp1 += lastp1;

    if (obs_homr > 1) {
      exit_threshx = exit_thresh - tailp2;
      do {
	curr_hets_t1 += 2;
	curr_homr_t1 -= 1;
	curr_homc_t1 -= 1;
	lastp1 *= (4 * curr_homr_t1 * curr_homc_t1) / (curr_hets_t1 * (curr_hets_t1 - 1));
	preaddp = tailp1;
	tailp1 += lastp1;
	if (tailp1 > exit_threshx) {
	  return 0;
	}
	if (tailp1 <= preaddp) {
	  break;
	}
      } while (curr_homr_t1 > 1.5);
    }
    if (tailp1 + tail2_ceil < exit_thresh) {
      return 1;
    }
    exit_threshx = exit_thresh - tailp1;
    while (curr_hets_t2 > 1) {
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= exit_threshx) {
	return 0;
      }
      if (tailp2 <= preaddp) {
	return 1;
      }
      curr_hets_t2 -= 2;
    }
    return 1;
  } else {
    if (!obs_homr) {
      return 0;
    }
    exit_thresh = (1 + (rare_copies / 2)) * thresh * EXACT_TEST_BIAS;
    do {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
          tailp2 = tailp1;
          centerp += tailp1;
	} else {
	  tailp2 = lastp2;
	}
	break;
      }
      centerp += lastp2;
      if (centerp > exit_thresh) {
	return 1;
      }
    } while (curr_homr_t2 > 0.5);
    exit_thresh = centerp / thresh;
    if (tailp1 + tailp2 >= exit_thresh) {
      return 0;
    }
    ratio = (4 * curr_homr_t2 * curr_homc_t2) / ((curr_hets_t2 + 2) * (curr_hets_t2 + 1));
    tail2_ceil = tailp2 + lastp2 * ratio / (1 - ratio);
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr + 1;
    curr_homc_t1 = obs_homc + 1;
    lastp1 = (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
    tail1_ceil = 2 * tailp1 / (1 - lastp1) - tailp1;
    lastp1 *= 2 * tailp1;
    tailp1 += lastp1;

    if (tail1_ceil + tail2_ceil < exit_thresh) {
      return 1;
    }
    if (obs_hets >= 4) {
      exit_threshx = exit_thresh - tailp2;
      do {
	curr_hets_t1 -= 2;
	curr_homr_t1 += 1;
	curr_homc_t1 += 1;
	lastp1 *= (curr_hets_t1 * (curr_hets_t1 - 1)) / (4 * curr_homr_t1 * curr_homc_t1);
	preaddp = tailp1;
        tailp1 += lastp1;
	if (tailp1 > exit_threshx) {
	  return 0;
	}
	if (tailp1 <= preaddp) {
	  break;
	}
      } while (curr_hets_t1 > 3.5);
    }
    if (tailp1 + tail2_ceil < exit_thresh) {
      return 1;
    }
    exit_threshx = exit_thresh - tailp1;
    while (curr_homr_t2 > 0.5) {
      curr_hets_t2 += 2;
      lastp2 *= (4 * curr_homr_t2 * curr_homc_t2) / (curr_hets_t2 * (curr_hets_t2 - 1));
      curr_homr_t2 -= 1;
      curr_homc_t2 -= 1;
      preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= exit_threshx) {
	return 0;
      }
      if (tailp2 <= preaddp) {
	return 1;
      }
    }
    return 1;
  }
}

double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp) {
  // Basic 2x2 Fisher exact test p-value calculation.
  double tprob = (1 - FISHER_EPSILON) * EXACT_TEST_BIAS;
  double cur_prob = tprob;
  double cprob = 0;
  int32_t tie_ct = 1;
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
      if (cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS) {
        tie_ct++;
      }
      tprob += cur_prob;
      break;
    }
    cprob += cur_prob;
  }
  if ((cprob == 0) && (!midp)) {
    return 1;
  }
  while (cur12 > 0.5) {
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
  }
  if (m11) {
    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;
    cur_prob = (1 - FISHER_EPSILON) * EXACT_TEST_BIAS;
    do {
      cur12 += 1;
      cur21 += 1;
      cur_prob *= (cur11 * cur22) / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
      preaddp = tprob;
      tprob += cur_prob;
      if (tprob <= preaddp) {
        if (!midp) {
	  return preaddp / (cprob + preaddp);
	} else {
          return (preaddp - ((1 - FISHER_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (cprob + preaddp);
	}
      }
    } while (cur11 > 0.5);
  }
  if (!midp) {
    return tprob / (cprob + tprob);
  } else {
    return (tprob - ((1 - FISHER_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (cprob + tprob);
  }
}

double fisher22_tail_pval(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, int32_t right_offset, double tot_prob_recip, double right_prob, uint32_t midp, uint32_t new_m11) {
  // Given that the left (w.r.t. m11) reference contingency table has
  // likelihood 1/tot_prob, the contingency table with m11 increased by
  // right_offset has likelihood right_prob/tot_prob, and the tails (up to but
  // not including the two references) sum to tail_sum/tot_prob, this
  // calculates the p-value of the given m11 (which must be on one tail).
  double left_prob = 1.0;
  double dxx = ((intptr_t)new_m11);
  double cur11;
  double cur12;
  double cur21;
  double cur22;
  double psum;
  double thresh;
  if (new_m11 < m11) {
    cur11 = ((intptr_t)m11);
    cur12 = ((intptr_t)m12);
    cur21 = ((intptr_t)m21);
    cur22 = ((intptr_t)m22);
    dxx += 0.5; // unnecessary (53 vs. 32 bits precision), but whatever
    do {
      cur12 += 1;
      cur21 += 1;
      left_prob *= cur11 * cur22 / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
    } while (cur11 > dxx);
    if (left_prob == 0) {
      return 0;
    }
    if (!midp) {
      psum = left_prob;
    } else {
      psum = left_prob * 0.5;
    }
    thresh = left_prob * (1 + FISHER_EPSILON);
    do {
      if (cur11 < 0.5) {
	break;
      }
      cur12 += 1;
      cur21 += 1;
      left_prob *= cur11 * cur22 / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
      dxx = psum;
      psum += left_prob;
    } while (psum > dxx);
    cur11 = ((intptr_t)(m11 + right_offset));
    cur12 = ((intptr_t)(m12 - right_offset));
    cur21 = ((intptr_t)(m21 - right_offset));
    cur22 = ((intptr_t)(m22 + right_offset));
    while (right_prob > thresh) {
      cur11 += 1;
      cur22 += 1;
      right_prob *= cur12 * cur21 / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
    }
    if (right_prob > 0) {
      if (midp && (right_prob < thresh * (1 - 2 * FISHER_EPSILON))) {
	psum += right_prob * 0.5;
      } else {
        psum += right_prob;
      }
      do {
	cur11 += 1;
	cur22 += 1;
	right_prob *= cur12 * cur21 / (cur11 * cur22);
	cur12 -= 1;
	cur21 -= 1;
	dxx = psum;
	psum += right_prob;
      } while (psum > dxx);
    }
  } else {
    dxx -= 0.5;
    cur11 = ((intptr_t)(m11 + right_offset));
    cur12 = ((intptr_t)(m12 - right_offset));
    cur21 = ((intptr_t)(m21 - right_offset));
    cur22 = ((intptr_t)(m22 + right_offset));
    do {
      cur11 += 1;
      cur22 += 1;
      right_prob *= cur12 * cur21 / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
    } while (cur11 < dxx);
    if (right_prob == 0) {
      return 0;
    }
    if (!midp) {
      psum = right_prob;
    } else {
      psum = right_prob * 0.5;
    }
    thresh = right_prob * (1 + FISHER_EPSILON);
    do {
      if (cur12 < 0.5) {
	break;
      }
      cur11 += 1;
      cur22 += 1;
      right_prob *= cur12 * cur21 / (cur11 * cur22);
      cur12 -= 1;
      cur21 -= 1;
      dxx = psum;
      psum += right_prob;
    } while (psum > dxx);
    cur11 = ((intptr_t)m11);
    cur12 = ((intptr_t)m12);
    cur21 = ((intptr_t)m21);
    cur22 = ((intptr_t)m22);
    while (left_prob > thresh) {
      cur12 += 1;
      cur21 += 1;
      left_prob *= cur11 * cur22 / (cur12 * cur21);
      cur11 -= 1;
      cur22 -= 1;
    }
    if (left_prob > 0) {
      if (midp && (left_prob < thresh * (1 - 2 * FISHER_EPSILON))) {
	psum += left_prob * 0.5;
      } else {
        psum += left_prob;
      }
      do {
	cur12 += 1;
	cur21 += 1;
	left_prob *= cur11 * cur22 / (cur12 * cur21);
	cur11 -= 1;
	cur22 -= 1;
	dxx = psum;
	psum += left_prob;
      } while (psum > dxx);
    }
  }
  return psum * tot_prob_recip;
}

void fisher22_precomp_pval_bounds(double pval, uint32_t midp, uint32_t row1_sum, uint32_t col1_sum, uint32_t total, uint32_t* bounds, double* tprobs) {
  // bounds[0] = m11 min
  // bounds[1] = m11 (max + 1)
  // bounds[2] = m11 min after including ties
  // bounds[3] = m11 (max + 1) after including ties
  // Treating m11 as the only variable, this returns the minimum and (maximum +
  // 1) values of m11 which are less extreme than the observed result, and
  // notes ties (2^{-40} tolerance).  Also, returns the values necessary for
  // invoking fisher22_tail_pval() with
  //   m11 := bounds[2] and
  //   right_offset := bounds[3] - bounds[2] - 1
  // in tprobs[0], [1], and [2] (if tprobs is not NULL).
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
  double cmult = midp? 0.5 : 1.0;
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
    bounds[0] = 0;
    bounds[1] = 0;
    bounds[2] = 0;
    bounds[3] = 1;
    // no need to initialize the other return values, they should never be used
    return;
  } else {
    if (pval == 0) {
      if (total >= row1_sum + col1_sum) {
	bounds[0] = 0;
	bounds[1] = MINV(row1_sum, col1_sum) + 1;
      } else {
	bounds[0] = row1_sum + col1_sum - total;
	bounds[1] = total - MAXV(row1_sum, col1_sum) + 1;
      }
      bounds[2] = bounds[0];
      bounds[3] = bounds[1];
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
  dxx = pval * tot_prob * (1 - FISHER_EPSILON / 2);
  threshold = pval * tot_prob * (1 + FISHER_EPSILON / 2);
  lii = 0;
  while (1) {
    if (left_prob < right_prob * (1 - FISHER_EPSILON / 2)) {
      cur_prob = tail_prob + left_prob * cmult;
      if (cur_prob > threshold) {
	break;
      }
      tail_prob += left_prob;
      uii = 1;
    } else if (right_prob < left_prob * (1 - FISHER_EPSILON / 2)) {
      cur_prob = tail_prob + right_prob * cmult;
      if (cur_prob > threshold) {
	break;
      }
      tail_prob += right_prob;
      uii = 2;
    } else {
      cur_prob = tail_prob + (left_prob + right_prob) * cmult;
      if (cur_prob > threshold) {
	if (left11 == right11) {
	  cur_prob = tail_prob + left_prob * cmult;
	  // center: left and right refer to same table.  subcases:
	  // 1. cur_prob > threshold: center table has less extreme pval.
	  //    lii = 0, both intervals size 1
	  // 2. dxx < cur_prob < threshold: center table has equal pval.
	  //    lii = 1, less-than interval size 0 but leq interval size 1
	  // 3. cur_prob < dxx: even centermost table has more extreme pval
	  //    (only possible with mid-p adj).
	  //    lii = 0, we increment left11 so both intervals size 0
	  if (cur_prob < threshold) {
	    if (cur_prob > dxx) {
	      lii = 1;
	    } else {
	      left11++;
	      left22++;
	      left_prob *= (left12 * left21) / (left11 * left22);
	    }
	  }
	}
	break;
      }
      tail_prob += left_prob + right_prob;
      uii = 3;
    }
    if (cur_prob > dxx) {
      lii = uii;
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
  bounds[2] = m11_offset + ((intptr_t)left11);
  bounds[3] = m11_offset + ((intptr_t)right11) + 1;
  bounds[0] = bounds[2] + (lii & 1);
  bounds[1] = bounds[3] - (lii >> 1);
  if (!tprobs) {
    return;
  }
  dxx = 1.0 / left_prob;
  tprobs[0] = left_prob / tot_prob;
  tprobs[1] = right_prob * dxx;
  /*
  if (lii & 1) {
    tail_prob -= left_prob;
  }
  if (lii >> 1) {
    tail_prob -= right_prob;
  }
  tprobs[2] = tail_prob * dxx;
  */
}

int32_t fisher23_tailsum(double* base_probp, double* saved12p, double* saved13p, double* saved22p, double* saved23p, double *totalp, uint32_t* tie_ctp, uint32_t right_side) {
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
	if (cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS) {
	  // throw in extra (1 - SMALL_EPSILON) multiplier to prevent rounding
	  // errors from causing this to keep going when the left-side test
	  // stopped
	  if (cur_prob > (1 - SMALL_EPSILON) * EXACT_TEST_BIAS) {
	    break;
	  }
          *tie_ctp += 1;
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
	if (cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS) {
	  if (cur_prob > EXACT_TEST_BIAS) {
	    break;
	  }
          *tie_ctp += 1;
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
  if (cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS) {
    if (cur_prob > EXACT_TEST_BIAS) {
      // even most extreme table on this side is too probable
      *totalp = 0;
      return 0;
    }
    *tie_ctp += 1;
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

double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23, uint32_t midp) {
  // 2x3 Fisher-Freeman-Halton exact test p-value calculation.
  // The number of tables involved here is still small enough that the network
  // algorithm (and the improved variants thereof that I've seen) are
  // suboptimal; a 2-dimensional version of the SNPHWE2 strategy has higher
  // performance.
  // 2x4, 2x5, and 3x3 should also be practical with this method, but beyond
  // that I doubt it's worth the trouble.
  // Complexity of approach is O(n^{df/2}), where n is number of observations.
  double cur_prob = (1 - FISHER_EPSILON) * EXACT_TEST_BIAS;
  double tprob = cur_prob;
  double cprob = 0;
  double dyy = 0;
  uint32_t tie_ct = 1;
  uint32_t dir = 0; // 0 = forwards, 1 = backwards
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
	if (cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS) {
          tie_ct++;
	}
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
          if (cur_prob > (1 - 2 * FISHER_EPSILON) * EXACT_TEST_BIAS) {
            tie_ct++;
	  }
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
      if (fisher23_tailsum(&base_probl, &savedl12, &savedl13, &savedl22, &savedl23, &dxx, &tie_ct, 0)) {
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
      fisher23_tailsum(&base_probr, &savedr12, &savedr13, &savedr22, &savedr23, &dyy, &tie_ct, 1);
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
  if (!midp) {
    return tprob / (tprob + cprob);
  } else {
    return (tprob - ((1 - FISHER_EPSILON) * EXACT_TEST_BIAS * 0.5) * ((int32_t)tie_ct)) / (tprob + cprob);
  }
}

void chi22_get_coeffs(intptr_t row1_sum, intptr_t col1_sum, intptr_t total, double* expm11p, double* recip_sump) {
  // chisq = (m11 - expm11)^2 * recip_sum
  // (see discussion for chi22_precomp_val_bounds() below.)
  //
  // expm11 = row1_sum * col1_sum / total
  // expm12 = row1_sum * col2_sum / total, etc.
  // recip_sum = 1 / expm11 + 1 / expm12 + 1 / expm21 + 1 / expm22
  // = total * (1 / (row1_sum * col1_sum) + 1 / (row1_sum * col2_sum) +
  //            1 / (row2_sum * col1_sum) + 1 / (row2_sum * col2_sum))
  // = total^3 / (row1_sum * col1_sum * row2_sum * col2_sum)
  double m11_numer = ((uint64_t)row1_sum) * ((uint64_t)col1_sum);
  double denom = m11_numer * (((uint64_t)(total - row1_sum)) * ((uint64_t)(total - col1_sum)));
  double dxx;
  if (denom != 0) {
    dxx = total;
    *expm11p = m11_numer / dxx;
    *recip_sump = dxx * dxx * dxx / denom;
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

double chi22_eval(intptr_t m11, intptr_t row1_sum, intptr_t col1_sum, intptr_t total) {
  double expm11_numer = ((uint64_t)row1_sum) * ((uint64_t)col1_sum);
  double denom = expm11_numer * (((uint64_t)(total - row1_sum)) * ((uint64_t)(total - col1_sum)));
  double dxx;
  double dyy;
  if (denom != 0) {
    dxx = total;
    dyy = m11 * dxx - expm11_numer; // total * (m11 - expm11)
    return (dyy * dyy * dxx) / denom;
  } else {
    return 0;
  }
}

double chi22_evalx(intptr_t m11, intptr_t row1_sum, intptr_t col1_sum, intptr_t total) {
  // PLINK emulation.  returns -9 instead of 0 if row1_sum, row2_sum, col1_sum,
  // or col2_sum is zero, for identical "NA" reporting.
  double expm11_numer = ((uint64_t)row1_sum) * ((uint64_t)col1_sum);
  double denom = expm11_numer * (((uint64_t)(total - row1_sum)) * ((uint64_t)(total - col1_sum)));
  double dxx;
  double dyy;
  if (denom != 0) {
    dxx = total;
    dyy = m11 * dxx - expm11_numer; // total * (m11 - expm11)
    return (dyy * dyy * dxx) / denom;
  } else {
    return -9;
  }
}

void chi22_precomp_val_bounds(double chisq, intptr_t row1_sum, intptr_t col1_sum, intptr_t total, uint32_t* bounds, double* coeffs) {
  // Treating m11 as the only variable, this returns the minimum and (maximum +
  // 1) values of m11 which produce smaller chisq statistics than given in
  // bounds[0] and bounds[1] respectively, and smaller-or-equal interval
  // bounds in bounds[2] and bounds[3].
  double expm11;
  double recip_sum;
  double cur11;
  double dxx;
  intptr_t ceil11;
  intptr_t lii;
  chi22_get_coeffs(row1_sum, col1_sum, total, &expm11, &recip_sum);
  if (coeffs) {
    coeffs[0] = expm11;
    coeffs[1] = recip_sum;
  }
  if (recip_sum == 0) {
    // sum-0 row or column, no freedom at all
    bounds[0] = (intptr_t)expm11;
    bounds[1] = bounds[0];
    bounds[2] = bounds[0];
    if (chisq == 0) {
      bounds[3] = bounds[0] + 1;
    } else {
      bounds[3] = bounds[0];
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
    bounds[0] = 0;
    bounds[2] = 0;
  } else {
    lii = (intptr_t)dxx;
    bounds[2] = lii;
    if (lii == (intptr_t)(cur11 + BIG_EPSILON)) {
      bounds[0] = lii + 1;
    } else {
      bounds[0] = lii;
    }
  }
  cur11 = expm11 + recip_sum;
  if (cur11 > ceil11 + BIG_EPSILON) {
    bounds[1] = ceil11 + 1;
    bounds[3] = bounds[1];
  } else {
    dxx = cur11 + 1 - BIG_EPSILON;
    lii = (intptr_t)dxx;
    bounds[1] = lii;
    if (lii == (intptr_t)(cur11 + BIG_EPSILON)) {
      bounds[3] = lii + 1;
    } else {
      bounds[3] = lii;
    }
  }
}

double chi23_eval(intptr_t m11, intptr_t m12, intptr_t row1_sum, intptr_t col1_sum, intptr_t col2_sum, intptr_t total) {
  // assumes no sum-zero row
  intptr_t m13 = row1_sum - m11 - m12;
  intptr_t col3_sum = total - col1_sum - col2_sum;
  double col1_sumd;
  double col2_sumd;
  double col3_sumd;
  double tot_recip;
  double dxx;
  double expect;
  double delta;
  double chisq;
  col1_sumd = col1_sum;
  col2_sumd = col2_sum;
  col3_sumd = col3_sum;
  tot_recip = 1.0 / ((double)total);
  dxx = row1_sum * tot_recip;
  expect = dxx * col1_sumd;
  delta = m11 - expect;
  chisq = delta * delta / expect;
  expect = dxx * col2_sumd;
  delta = m12 - expect;
  chisq += delta * delta / expect;
  expect = dxx * col3_sumd;
  delta = m13 - expect;
  chisq += delta * delta / expect;
  dxx = (total - row1_sum) * tot_recip;
  expect = dxx * col1_sumd;
  delta = (col1_sum - m11) - expect;
  chisq += delta * delta / expect;
  expect = dxx * col2_sumd;
  delta = (col2_sum - m12) - expect;
  chisq += delta * delta / expect;
  expect = dxx * col3_sumd;
  delta = (col3_sum - m13) - expect;
  chisq += delta * delta / expect;
  if (chisq < (SMALL_EPSILON * SMALL_EPSILON)) {
    return 0;
  }
  return chisq;
}

void chi23_evalx(intptr_t m11, intptr_t m12, intptr_t m13, intptr_t m21, intptr_t m22, intptr_t m23, double* chip, uint32_t* dfp) {
  // Slightly different from PLINK calculation, since it detects lone nonzero
  // columns.
  intptr_t row1_sum = m11 + m12 + m13;
  intptr_t row2_sum = m21 + m22 + m23;
  intptr_t col1_sum = m11 + m21;
  intptr_t col2_sum = m12 + m22;
  intptr_t col3_sum = m13 + m23;
  intptr_t total;
  double col1_sumd;
  double col2_sumd;
  double col3_sumd;
  double tot_recip;
  double dxx;
  double expect;
  double delta;
  double chisq;
  if ((!row1_sum) || (!row2_sum)) {
    *chip = -9;
    *dfp = 0;
    return;
  }
  total = row1_sum + row2_sum;
  if (!col1_sum) {
    *chip = chi22_evalx(m12, row1_sum, col2_sum, total);
    if (*chip != -9) {
      *dfp = 1;
    } else {
      *dfp = 0;
    }
    return;
  } else if ((!col2_sum) || (!col3_sum)) {
    *chip = chi22_evalx(m11, row1_sum, col1_sum, total);
    if (*chip != -9) {
      *dfp = 1;
    } else {
      *dfp = 0;
    }
    return;
  }
  col1_sumd = col1_sum;
  col2_sumd = col2_sum;
  col3_sumd = col3_sum;
  tot_recip = 1.0 / ((double)total);
  dxx = row1_sum * tot_recip;
  expect = dxx * col1_sumd;
  delta = m11 - expect;
  chisq = delta * delta / expect;
  expect = dxx * col2_sumd;
  delta = m12 - expect;
  chisq += delta * delta / expect;
  expect = dxx * col3_sumd;
  delta = m13 - expect;
  chisq += delta * delta / expect;
  dxx = row2_sum * tot_recip;
  expect = dxx * col1_sumd;
  delta = m21 - expect;
  chisq += delta * delta / expect;
  expect = dxx * col2_sumd;
  delta = m22 - expect;
  chisq += delta * delta / expect;
  expect = dxx * col3_sumd;
  delta = m23 - expect;
  chisq += delta * delta / expect;
  if (chisq < (SMALL_EPSILON * SMALL_EPSILON)) {
    chisq = 0;
  }
  *chip = chisq;
  *dfp = 2;
}

double ca_trend_eval(intptr_t case_dom_ct, intptr_t case_ct, intptr_t het_ct, intptr_t homdom_ct, intptr_t total) {
  // case_dom_ct is an allele count (2 * homa2 + het), while other inputs are
  // observation counts.
  //
  // If case_missing_ct is fixed,
  //   row1_sum = case ct
  //   col1_sum = A2 ct
  //   case_ct * ctrl_ct * REC_ct * DOM_ct
  //   REC_ct = 2 * obs_11 + obs_12
  //   DOM_ct = 2 * obs_22 + obs_12
  //   CA = (obs_U / obs_T) * (case REC ct) - (obs_A / obs_T) * (ctrl DOM ct)
  //      = (case A2) * (obs_U / obs_T) - (obs_A / obs_T) * (DOM ct - case DOM)
  //      = (case A2) * (obs_U / obs_T) + (case DOM) * (obs_A / obs_T) - DOM*(A/T)
  //      = (case A2 ct) - total A2 ct * (A/T)
  //   CAT = CA * obs_T
  //   varCA_recip = obs_T * obs_T * obs_T /
  //     (obs_A * obs_U * (obs_T * (obs_12 + 4 * obs_22) - DOMct * DOMct))
  //   trend statistic = CAT * CAT * [varCA_recip / obs_T^2]
  double dom_ct = het_ct + 2 * homdom_ct;
  double totald = total;
  double case_ctd = case_ct;
  double cat = case_dom_ct * totald - dom_ct * case_ctd;
  double dxx = totald * (het_ct + 4 * ((int64_t)homdom_ct)) - dom_ct * dom_ct;

  // This should never be called with dxx == 0 (which happens when two columns
  // are all-zero).  Use ca_trend_evalx() to check for that.
  dxx *= case_ctd * (totald - case_ctd);
  return cat * cat * totald / dxx;
}

double ca_trend_evalx(intptr_t case_dom_ct, intptr_t case_ct, intptr_t het_ct, intptr_t homdom_ct, intptr_t total) {
  double dom_ct = het_ct + 2 * homdom_ct;
  double totald = total;
  double case_ctd = case_ct;
  double cat = case_dom_ct * totald - dom_ct * case_ctd;
  double dxx = totald * (het_ct + 4 * ((int64_t)homdom_ct)) - dom_ct * dom_ct;
  if (dxx != 0) {
    dxx *= case_ctd * (totald - case_ctd);
    return cat * cat * totald / dxx;
  } else {
    return -9;
  }
}

void ca_trend_precomp_val_bounds(double chisq, intptr_t case_ct, intptr_t het_ct, intptr_t homdom_ct, intptr_t total, uint32_t* bounds, double* coeffs) {
  // If case_missing_ct is fixed,
  //   row1_sum = case ct
  //   col1_sum = DOM ct
  //   case_ct * ctrl_ct * REC_ct * DOM_ct
  //   REC_ct = 2 * obs_11 + obs_12
  //   DOM_ct = 2 * obs_22 + obs_12
  //   CA = (obs_U / obs_T) * (case DOM ct) - (obs_A / obs_T) * (ctrl DOM ct)
  //      = (case DOM) * (obs_U / obs_T) - (obs_A / obs_T) * (DOM ct - case DOM)
  //      = (case DOM) * (obs_U / obs_T) + (case DOM) * (obs_A / obs_T) - DOM*(A/T)
  //      = (case DOM ct) - total DOM ct * (A/T)
  //   varCA_recip = obs_T * obs_T * obs_T /
  //     (obs_A * obs_U * (obs_T * (obs_12 + 4 * obs_22) - DOMct * DOMct))
  //   trend statistic = CA * CA * varCA_recip
  intptr_t dom_ct = het_ct + 2 * homdom_ct;
  double dom_ctd = dom_ct;
  double totald = total;
  double case_ctd = case_ct;
  double tot_recip = 1.0 / totald;
  double expm11 = dom_ctd * case_ctd * tot_recip;
  double dxx = case_ctd * (totald - case_ctd) * (totald * (het_ct + 4 * ((int64_t)homdom_ct)) - dom_ctd * dom_ctd);
  double varca_recip;
  double cur11;
  intptr_t ceil11;
  intptr_t lii;
  if (dxx == 0) {
    // bounds/coefficients should never be referenced in this case
    return;
  }
  varca_recip = totald * totald * totald / dxx;
  if (coeffs) {
    coeffs[0] = expm11;
    coeffs[1] = varca_recip;
  }

  // statistic: (cur11 - expm11)^2 * varca_recip
  ceil11 = case_ct * 2;
  if (dom_ct < ceil11) {
    ceil11 = dom_ct;
  }
  // chisq = (cur11 - expm11)^2 * varca_recip
  // -> expm11 +/- sqrt(chisq / varca_recip) = cur11
  varca_recip = sqrt(chisq / varca_recip);
  cur11 = expm11 - varca_recip;
  dxx = cur11 + 1 - BIG_EPSILON;
  if (dxx < 0) {
    bounds[0] = 0;
    bounds[2] = 0;
  } else {
    lii = (intptr_t)dxx;
    bounds[2] = lii;
    if (lii == (intptr_t)(cur11 + BIG_EPSILON)) {
      bounds[0] = lii + 1;
    } else {
      bounds[0] = lii;
    }
  }
  cur11 = expm11 + varca_recip;
  if (cur11 > ceil11 + BIG_EPSILON) {
    bounds[1] = ceil11 + 1;
    bounds[3] = bounds[1];
  } else {
    dxx = cur11 + 1 - BIG_EPSILON;
    lii = (intptr_t)dxx;
    bounds[1] = lii;
    if (lii == (intptr_t)(cur11 + BIG_EPSILON)) {
      bounds[3] = lii + 1;
    } else {
      bounds[3] = lii;
    }
  }
}

uint32_t linear_hypothesis_chisq(uintptr_t constraint_ct, uintptr_t param_ct, double* constraints_con_major, double* coef, double* cov_matrix, double* param_df_buf, double* param_df_buf2, double* df_df_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* df_buf, double* chisq_ptr) {
  // See PLINK model.cpp Model::linearHypothesis().
  //
  // outer = df_buf
  // inner = df_df_buf
  // tmp = param_df_buf
  // mi_buf only needs to be of length constraint_ct
  //
  // Since no PLINK function ever calls this with nonzero h[] values, this just
  // takes a df parameter for now; it's trivial to switch to the more general
  // interface later.
  double* dptr = constraints_con_major;
  double* dptr2;
  uintptr_t constraint_idx;
  uintptr_t constraint_idx2;
  uintptr_t param_idx;
  double dxx;
  double dyy;
  for (constraint_idx = 0; constraint_idx < constraint_ct; constraint_idx++) {
    dxx = 0;
    dptr2 = coef;
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      dxx += (*dptr++) * (*dptr2++);
    }
    df_buf[constraint_idx] = dxx;
  }
  // temporarily set param_df_buf2[][] to H-transpose
  transpose_copy(constraint_ct, param_ct, constraints_con_major, param_df_buf2);
  col_major_matrix_multiply(constraint_ct, param_ct, param_ct, param_df_buf2, cov_matrix, param_df_buf);
  // tmp[][] is now param-major
  col_major_matrix_multiply(constraint_ct, constraint_ct, param_ct, param_df_buf, constraints_con_major, df_df_buf);

  if (invert_matrix((uint32_t)constraint_ct, df_df_buf, mi_buf, param_df_buf2)) {
    return 1;
  }
  dxx = 0; // result
  dptr = df_df_buf;
  for (constraint_idx = 0; constraint_idx < constraint_ct; constraint_idx++) {
    dyy = 0; // tmp2[c]
    dptr2 = df_buf;
    for (constraint_idx2 = 0; constraint_idx2 < constraint_ct; constraint_idx2++) {
      dyy += (*dptr++) * (*dptr2++);
    }
    dxx += dyy * df_buf[constraint_idx];
  }
  *chisq_ptr = dxx;
  return 0;
}

double binom_2sided(uint32_t succ, uint32_t obs, uint32_t midp) {
  // straightforward to generalize this to any success probability
  double cur_succ_t2 = (int32_t)succ;
  double cur_fail_t2 = (int32_t)(obs - succ);
  double tailp = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;
  int32_t tie_ct = 1;
  // double rate_mult_incr = rate / (1 - rate);
  // double rate_mult_decr = (1 - rate) / rate;
  double cur_succ_t1;
  double cur_fail_t1;
  double preaddp;
  if (!obs) {
    return midp? 0.5 : 1;
  }
  if (obs < succ * 2) {
  // if (obs * rate < succ) {
    while (cur_succ_t2 > 0.5) {
      cur_fail_t2 += 1;
      lastp2 *= cur_succ_t2 / cur_fail_t2;
      // lastp2 *= rate_mult_decr * cur_succ_t2 / cur_fail_t2;
      cur_succ_t2 -= 1;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
	  tie_ct++;
	}
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (cur_succ_t2 > 0.5) {
      cur_fail_t2 += 1;
      lastp2 *= cur_succ_t2 / cur_fail_t2;
      // lastp2 *= rate_mult_decr * cur_succ_t2 / cur_fail_t2;
      cur_succ_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    cur_succ_t1 = (int32_t)(succ + 1);
    cur_fail_t1 = (int32_t)(obs - succ);
    while (cur_fail_t1 > 0.5) {
      lastp1 *= cur_fail_t1 / cur_succ_t1;
      // lastp1 *= rate_mult_incr * cur_fail_t1 / cur_succ_t1;
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      cur_succ_t1 += 1;
      cur_fail_t1 -= 1;
    }
  } else {
    while (cur_fail_t2 > 0.5) {
      cur_succ_t2++;
      lastp2 *= cur_fail_t2 / cur_succ_t2;
      // lastp2 *= rate_mult_incr * cur_fail_t2 / cur_succ_t2;
      cur_fail_t2--;
      if (lastp2 < EXACT_TEST_BIAS) {
	if (lastp2 > (1 - 2 * SMALL_EPSILON) * EXACT_TEST_BIAS) {
	  tie_ct++;
	}
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if ((centerp == 0) && (!midp)) {
      return 1;
    }
    while (cur_fail_t2 > 0.5) {
      cur_succ_t2 += 1;
      lastp2 *= cur_fail_t2 / cur_succ_t2;
      // lastp2 *= rate_mult_incr * cur_fail_t2 / cur_succ_t2;
      cur_fail_t2 -= 1;
      preaddp = tailp;
      tailp += lastp2;
      if (tailp <= preaddp) {
	break;
      }
    }
    cur_succ_t1 = (int32_t)succ;
    cur_fail_t1 = (int32_t)(obs - succ);
    while (cur_succ_t1 > 0.5) {
      cur_fail_t1 += 1;
      lastp1 *= cur_succ_t1 / cur_fail_t1;
      // lastp1 *= rate_mult_decr * cur_succ_t1 / cur_fail_t1;
      preaddp = tailp;
      tailp += lastp1;
      if (tailp <= preaddp) {
	break;
      }
      cur_succ_t1 -= 1;
    }
  }
  if (!midp) {
    return tailp / (tailp + centerp);
  } else {
    return (tailp - ((1 - SMALL_EPSILON) * EXACT_TEST_BIAS * 0.5) * tie_ct) / (tailp + centerp);
  }
}
