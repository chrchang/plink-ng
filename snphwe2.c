// Uncomment this if you want to compile a standalone program.
#define TEST_BUILD

#include <stdint.h>
#include <inttypes.h>

// for INFINITY
#include <math.h>

#ifdef TEST_BUILD
// just for printf() in sample program
#include <stdio.h>
#include <stdlib.h>
#endif

#define EPSILON 0.0000000000001

// 53-bit double precision limit
#define DOUBLE_PREC_LIMIT 0.00000000000000011102230246251565404236316680908203125
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2) {
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
  // integer overflow.
  // - Detection and efficient handling of floating point overflow and
  // underflow.  E.g. instead of summing a tail all the way down, the loop
  // stops once the latest increment underflows the partial sum's 53-bit
  // precision; this results in a large speedup when max heterozygote count
  // >1k.
  // - No malloc() call.  It's only necessary to keep track of a few partial
  // sums.
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
  double curr_hets_t2 = obs_hets;
  double curr_homr_t2 = obs_homr;
  double curr_homc_t2 = obs_homc;
  double tailp = (1 - EPSILON) * EXACT_TEST_BIAS;
  double centerp = 0;
  double lastp2 = tailp;
  double lastp1 = tailp;
  double curr_hets_t1;
  double curr_homr_t1;
  double curr_homc_t1;
  double preaddp;
  if (!genotypes2) {
    return 1;
  }

  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper
    while (curr_hets_t2 > 1.5) {
      // het_probs[curr_hets] = 1
      // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets
      curr_homr_t2 += 1;
      curr_homc_t2 += 1;
      lastp2 *= (curr_hets_t2 * (curr_hets_t2 - 1)) / (4 * curr_homr_t2 * curr_homc_t2);
      curr_hets_t2 -= 2;
      if (lastp2 < EXACT_TEST_BIAS) {
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if (centerp == 0) {
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
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 0;
      }
    }
    if (centerp == 0) {
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
  return tailp / (tailp + centerp);
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
  // the threshold; we never need to sum the remainder of the tails.
  // - And in the case of disequilibrium, we will often be able to immediately
  // determine that the tail sum cannot possibly pass the threshold, just by
  // looking at the tail boundary elements and using a geometric series to
  // upper-bound the tail sums.
  // - Only when neither of these conditions hold do we start traveling down
  // the tails.
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
  double tailp1 = (1 - EPSILON) * EXACT_TEST_BIAS;
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
    // report test failure before summing the entire center.
    // The immediate tail is easy: if the next element out on the tail is r,
    // then 1 + r + r^2 + ... = 1 / (1-r) works.
    // For the far tail, we currently use the trivial bound of
    //   1 + floor[het_exp_floor / 2]
    // (each far tail element must be no greater than 1 and there are at
    // most that many of them).  This bound could be improved, but it might not
    // be worth the precomputational effort.
    // ...and as long as we're using such a weak bound for the far tail,
    // there's no point to carefully calculating the near tail since we're very
    // unlikely to use the partial result before exiting from the function.
    exit_thresh = rare_copies * thresh * EXACT_TEST_BIAS;

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
    exit_thresh = rare_copies * thresh * EXACT_TEST_BIAS;
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

#ifdef TEST_BUILD
int main(int argc, char** argv) {
  FILE* test_file;
  int32_t hets;
  int32_t homs1;
  int32_t homs2;
  double p_value;
  char name[100];
  if ((argc == 4) || (argc == 5)) {
    hets = atoi(argv[1]);
    homs1 = atoi(argv[2]);
    homs2 = atoi(argv[3]);
    if ((hets < 0) || (homs1 < 0) || (homs2 < 0)) {
      printf("Error: Negative parameter.\n");
      return 3;
    }
    if (argc == 4) {
      p_value = SNPHWE2(hets, homs1, homs2);
      printf("P-value: %lg\n", p_value);
    } else {
      if (sscanf(argv[4], "%lg", &p_value) != 1) {
	printf("Error: Invalid p-value '%s'.\n", argv[4]);
	return 3;
      }
      if (SNPHWE_t(hets, homs1, homs2, p_value)) {
	printf("Test failed (p-value below threshold).\n");
      } else {
	printf("Test passed (p-value above threshold).\n");
      }
    }
    return 0;
  } else if (argc != 2) {
    printf(
"SNPHWE2 demo    http://www.sph.umich.edu/csg/abecasis/Exact/\n\n"
"  snphwe2 [het count] [hom count 1] [hom count 2] {threshold}\n"
"  snphwe2 [filename]\n\n"
"If a filename is provided, the file is expected to contain marker names in the\n"
"first column, heterozygote counts in the second, and homozygote counts in the\n"
"third and fourth.\n");
    return 1;
  }
  test_file = fopen(argv[1], "r");
  if (!test_file) {
    printf("Error: Unable to open file.\n");
    return 2;
  }
  while (!feof(test_file)) {
    fscanf(test_file, "%s %d %d %d\n", name, &hets, &homs1, &homs2);
    p_value = SNPHWE2(hets, homs1, homs2);
    printf("P-value for %s: %lg\n", name, p_value);
  }
  return 0;
}
#endif
