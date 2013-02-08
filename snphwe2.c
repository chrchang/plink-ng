// Uncomment this if you want to compile a standalone program.
// #define TEST_BUILD

#include <stdint.h>
#include <inttypes.h>

// for INFINITY
#include <math.h>

#ifdef TEST_BUILD
// just for printf() in sample program
#include <stdio.h>
#endif

#define EPSILON 0.0000000000001

// 53-bit double precision limit
#define DOUBLE_PREC_LIMIT 0.00000000000000011102230246251565404236316680908203125

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
  uintptr_t obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  uintptr_t obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
  uintptr_t rare_copies = 2 * obs_homr + obs_hets;
  uint32_t genotypes = obs_hets + obs_homc + obs_homr;
  uintptr_t genotypes2 = genotypes * 2;
  uintptr_t curr_hets_t2 = obs_hets;
  uintptr_t curr_homr_t2 = obs_homr;
  uintptr_t curr_homc_t2 = obs_homc;
  double tailp = 1;
  double centerp = 0;
  double lastp2 = 1;
  double lastp1 = 1;
  uintptr_t curr_hets_t1;
  uintptr_t curr_homr_t1;
  uintptr_t curr_homc_t1;
  uint32_t het_exp_floor;
  double tail_stop;
  if (!genotypes) {
    return 1;
  }

  het_exp_floor = (((uint64_t)rare_copies) * (genotypes2 - rare_copies)) / genotypes2;
  if (curr_hets_t2 > het_exp_floor) {
    // tail 1 = upper
    while (curr_hets_t2 > 1) {
      // het_probs[curr_hets] = 1
      // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets
      lastp2 *= ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1))) / ((double)((4LLU * (++curr_homr_t2)) * (++curr_homc_t2)));
      curr_hets_t2 -= 2;
      if (lastp2 < 1 + EPSILON) {
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
    tail_stop = tailp * DOUBLE_PREC_LIMIT;
    while (curr_hets_t2 > 1) {
      lastp2 *= ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1))) / ((double)((4LLU * (++curr_homr_t2)) * (++curr_homc_t2)));
      curr_hets_t2 -= 2;
      if (lastp2 < tail_stop) {
	break;
      }
      tailp += lastp2;
    }
    tail_stop = tailp * DOUBLE_PREC_LIMIT;
    curr_hets_t1 = obs_hets + 2;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_homr_t1) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      lastp1 *= ((double)((4LLU * (curr_homr_t1--)) * (curr_homc_t1--))) / ((double)(((uint64_t)curr_hets_t1) * (curr_hets_t1 - 1)));
      if (lastp1 < tail_stop) {
	break;
      }
      tailp += lastp1;
      curr_hets_t1 += 2;
    }
  } else {
    // tail 1 = lower
    while (curr_homr_t2) {
      curr_hets_t2 += 2;
      lastp2 *= ((double)((4LLU * (curr_homr_t2--)) * (curr_homc_t2--))) / ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1)));
      if (lastp2 < 1 + EPSILON) {
	tailp += lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 1;
      }
    }
    if (centerp == 0) {
      return 1;
    }
    tail_stop = tailp * DOUBLE_PREC_LIMIT;
    while (curr_homr_t2) {
      curr_hets_t2 += 2;
      lastp2 *= ((double)((4LLU * (curr_homr_t2--)) * (curr_homc_t2--))) / ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1)));
      if (lastp2 < tail_stop) {
	break;
      }
      tailp += lastp2;
    }
    tail_stop = tailp * DOUBLE_PREC_LIMIT;
    curr_hets_t1 = obs_hets;
    curr_homr_t1 = obs_homr;
    curr_homc_t1 = obs_homc;
    while (curr_hets_t1 > 1) {
      lastp1 *= ((double)(((uint64_t)curr_hets_t1) * (curr_hets_t1 - 1))) / ((double)((4LLU * (++curr_homr_t1)) * (++curr_homc_t1)));
      if (lastp1 < tail_stop) {
	break;
      }
      tailp += lastp1;
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
  uintptr_t obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  uintptr_t obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
  uintptr_t rare_copies = 2 * obs_homr + obs_hets;
  uint32_t genotypes = obs_hets + obs_homc + obs_homr;
  uintptr_t genotypes2 = genotypes * 2;
  uintptr_t curr_hets_t2 = obs_hets; // tail 2
  uintptr_t curr_homr_t2 = obs_homr;
  uintptr_t curr_homc_t2 = obs_homc;
  double tailp1 = 1;
  double centerp = 0;
  double lastp2 = 1;
  double tailp2 = 0;
  double tail1_ceil;
  double tail2_ceil;
  double lastp1;
  uintptr_t curr_hets_t1;
  uintptr_t curr_homr_t1;
  uintptr_t curr_homc_t1;
  double tail_thresh; // as soon as tail sum reaches this, we've passed
  double tail_threshx;
  double tail_stop;
  uint32_t het_exp_floor;
  double ratio;
  if (!genotypes) {
    return 0;
  }

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
  het_exp_floor = (((uint64_t)rare_copies) * (genotypes2 - rare_copies)) / genotypes2;
  if (curr_hets_t2 > het_exp_floor) {
    // tail 1 = upper
    if (curr_hets_t2 < 2) {
      return 0;
    }
    // het_probs[curr_hets] = 1
    // het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1) / (4 * (curr_homr + 1) * (curr_homc + 1))
    do {
      lastp2 *= ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1))) / ((double)((4LLU * (++curr_homr_t2)) * (++curr_homc_t2)));
      curr_hets_t2 -= 2;
      if (lastp2 < 1 + EPSILON) {
	tailp2 = lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 1;
      }
    } while (curr_hets_t2 > 1);
    tail_thresh = centerp * thresh / (1 - thresh);
    if (1 + tailp2 >= tail_thresh) {
      return 0;
    }
    // c + cr + cr^2 + ... = c/(1-r), which is an upper bound for the tail sum
    ratio = ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1))) / ((double)((4LLU * (curr_homr_t2 + 1)) * (curr_homc_t2 + 1)));
    tail2_ceil = tailp2 / (1 - ratio);
    if (obs_homr) {
      // het_probs[curr_hets + 2] = het_probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
      curr_hets_t1 = obs_hets + 2;
      curr_homr_t1 = obs_homr;
      curr_homc_t1 = obs_homc;
      lastp1 = ((double)((4LLU * (curr_homr_t1--)) * (curr_homc_t1--))) / ((double)(((uint64_t)curr_hets_t1) * (curr_hets_t1 - 1)));
      tail1_ceil = 1 / (1 - lastp1);
      if (tail1_ceil + tail2_ceil < tail_thresh) {
	return 1;
      }
      tail_threshx = tail_thresh - tailp2;
      tail_stop = (1 + tailp2) * DOUBLE_PREC_LIMIT;
      while (1) {
        tailp1 += lastp1;
	if (tailp1 >= tail_threshx) {
	  return 0;
	}
	if (!curr_homr_t1) {
	  break;
	}
	curr_hets_t1 += 2;
	lastp1 *= ((double)((4LLU * (curr_homr_t1--)) * (curr_homc_t1--))) / ((double)(((uint64_t)curr_hets_t1) * (curr_hets_t1 - 1)));
	if (lastp1 < tail_stop) {
	  break;
	}
      }
    }
    if (tailp1 + tail2_ceil < tail_thresh) {
      return 1;
    }
    tail_threshx = tail_thresh - tailp1;
    tail_stop = (tailp1 + tailp2) * DOUBLE_PREC_LIMIT;
    while (curr_hets_t2 > 1) {
      lastp2 *= ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1))) / ((double)((4LLU * (++curr_homr_t2)) * (++curr_homc_t2)));
      if (lastp2 < tail_stop) {
	return 1;
      }
      tailp2 += lastp2;
      if (tailp2 >= tail_threshx) {
	return 0;
      }
      curr_hets_t2 -= 2;
    }
    return 1;
  } else {
    // tail 1 = lower
    if (!curr_homr_t2) {
      return 0;
    }
    do {
      curr_hets_t2 += 2;
      lastp2 *= ((double)((4LLU * (curr_homr_t2--)) * (curr_homc_t2--))) / ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1)));
      if (lastp2 < 1 + EPSILON) {
	tailp2 = lastp2;
	break;
      }
      centerp += lastp2;
      if (centerp == INFINITY) {
	return 1;
      }
    } while (curr_homr_t2);
    tail_thresh = centerp * thresh / (1 - thresh);
    if (1 + tailp2 >= tail_thresh) {
      return 0;
    }
    ratio = ((double)((4LLU * curr_homr_t2) * curr_homc_t2)) / ((double)(((uint64_t)(curr_hets_t2 + 2)) * (curr_hets_t2 + 1)));
    tail2_ceil = tailp2 / (1 - ratio);
    if (obs_hets >= 2) {
      curr_hets_t1 = obs_hets;
      curr_homr_t1 = obs_homr;
      curr_homc_t1 = obs_homc;
      lastp1 = ((double)(((uint64_t)curr_hets_t1) * (curr_hets_t1 - 1))) / ((double)((4LLU * (++curr_homr_t1)) * (++curr_homc_t1)));
      curr_hets_t1 -= 2;
      tail1_ceil = tailp1 / (1 - lastp1);
      if (tail1_ceil + tail2_ceil < tail_thresh) {
	return 1;
      }
      tail_threshx = tail_thresh - tailp2;
      tail_stop = (1 + tailp2) * DOUBLE_PREC_LIMIT;
      while (1) {
        tailp1 += lastp1;
        if (tailp1 >= tail_threshx) {
	  return 0;
        }
	if (curr_hets_t1 < 2) {
	  break;
	}
	lastp1 *= ((double)(((uint64_t)curr_hets_t1) * (curr_hets_t1 - 1))) / ((double)((4LLU * (++curr_homr_t1)) * (++curr_homc_t1)));
	if (lastp1 < tail_stop) {
	  break;
	}
	curr_hets_t1 -= 2;
      }
    }
    if (tailp1 + tail2_ceil < tail_thresh) {
      return 1;
    }
    tail_threshx = tail_thresh - tailp1;
    tail_stop = (tailp1 + tailp2) * DOUBLE_PREC_LIMIT;
    while (curr_homr_t2) {
      curr_hets_t2 += 2;
      lastp2 *= ((double)((4LLU * (curr_homr_t2--)) * (curr_homc_t2--))) / ((double)(((uint64_t)curr_hets_t2) * (curr_hets_t2 - 1)));
      if (lastp2 < tail_stop) {
	return 1;
      }
      tailp2 += lastp2;
      if (tailp2 >= tail_threshx) {
	return 0;
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
  if (argc != 2) {
    printf(
"SNPHWE2 demo    http://www.sph.umich.edu/csg/abecasis/Exact/\n\n"
"  snphwe2 [filename]\n\n"
"The file is expected to contain marker names in the first column, heterozygote\n"
"counts in the second, and homozygote counts in the third and fourth.\n");
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
