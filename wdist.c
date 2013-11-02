// WDIST weighted genomic distance calculator
// Copyright (C) 2013  Christopher Chang  chrchang@alumni.caltech.edu

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Uncomment "#define NOLAPACK" in wdist_common.h to build without LAPACK.

#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include "wdist_common.h"

#ifdef __APPLE__
  #include <sys/types.h>
  #include <sys/sysctl.h>
#endif

#include "wdist_assoc.h"
#include "wdist_calc.h"
#include "wdist_cnv.h"
#include "wdist_data.h"
#include "wdist_dosage.h"
#include "wdist_help.h"
#include "wdist_homozyg.h"
#include "wdist_lasso.h"
#include "wdist_stats.h"
#include "pigz.h"

// default jackknife iterations
#define ITERS_DEFAULT 100000
#define DEFAULT_PPC_GAP 500000
#define MAX_PCS_DEFAULT 20

#define DEFAULT_IBS_TEST_PERMS 100000

#define LOAD_RARE_GRM 1
#define LOAD_RARE_GRM_BIN 2
#define LOAD_RARE_LGEN 4
#define LOAD_RARE_TRANSPOSE 8
#define LOAD_RARE_TPED 0x10
#define LOAD_RARE_TFAM 0x20
#define LOAD_RARE_TRANSPOSE_MASK (LOAD_RARE_TRANSPOSE | LOAD_RARE_TPED | LOAD_RARE_TFAM)
#define LOAD_RARE_DUMMY 0x40
#define LOAD_RARE_SIMULATE 0x80
#define LOAD_RARE_CNV 0x100
#define LOAD_RARE_GVAR 0x200
#define LOAD_RARE_23 0x400

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^31 - 1
#define PARALLEL_MAX 32768

const char ver_str[] =
#ifdef PLINK_BUILD
  "PLINK v1.50b1"
#else
  #ifdef STABLE_BUILD
  "WDIST v0.22.7"
  #else
  "WDIST v0.23.0p"
  #endif
#endif
#ifdef NOLAPACK
  "NL"
#endif
#ifdef __LP64__
  " 64-bit"
#else
  " 32-bit"
#endif
  " (2 Nov 2013)";
const char ver_str2[] =
#ifdef PLINK_BUILD
  "               [final website TBD]\n"
  "(C) 2005-2013 Shaun Purcell, Christopher Chang    GNU GPL version 3\n";
#else
  "    https://www.cog-genomics.org/wdist\n"
  "(C) 2013 Christopher Chang, GNU General Public License version 3\n";
#endif
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
const char errstr_filter_format[] = "Error: Improperly formatted filter file.\n";
const char errstr_freq_format[] = "Error: Improperly formatted frequency file.\n";
const char notestr_null_calc[] = "Note: No output requested.  Exiting.\n";
#ifdef STABLE_BUILD
const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --merge-list, --write-snplist, --freqx,\n--missing, --hardy, --ibc, --distance, --genome, --homozyg, --cluster,\n--neighbour, --model, --gxe, --logistic, --lasso, --indep, --make-rel,\n--make-grm-gz, --rel-cutoff, --regress-distance, and --ibs-test.\n\n'" PROG_NAME_STR " --help | more' describes all functions (warning: long).\n";
#else
  #ifndef NOLAPACK
const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --merge-list, --write-snplist, --freqx,\n--missing, --hardy, --ibc, --distance, --genome, --homozyg, --cluster,\n--neighbour, --model, --gxe, --logistic, --lasso, --indep, --make-rel,\n--make-grm-gz, --rel-cutoff, --regress-pcs, --regress-distance, --ibs-test,\nand --unrelated-heritability.\n\n'" PROG_NAME_STR " --help | more' describes all functions (warning: long).\n";
  #else
const char notestr_null_calc2[] = "Commands include --make-bed, --recode, --merge-list, --write-snplist, --freqx,\n--missing, --hardy, --ibc, --distance, --genome, --homozyg, --cluster,\n--neighbour, --model, --gxe, --logistic, --lasso, --indep, --make-rel,\n--make-grm-gz, --rel-cutoff, --regress-pcs, --regress-distance, and\n--ibs-test.\n\n'" PROG_NAME_STR " --help | more' describes all functions (warning: long).\n";
  #endif
#endif

intptr_t malloc_size_mb = 0;

unsigned char* wkspace;

void dispmsg(int32_t retval) {
  switch (retval) {
  case RET_NOMEM:
    logprint("\nError: Out of memory.  Try the --memory and/or --parallel flags.\n");
    break;
  case RET_WRITE_FAIL:
    logprint("\nError: File write failure.\n");
    break;
  case RET_READ_FAIL:
    logprint("\nError: File read failure.\n");
    break;
  }
}

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
  double tailp = (1 - SMALL_EPSILON) * EXACT_TEST_BIAS;
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

// back to our regular program

inline char* read_next_unsafe(char* target, char* source) {
  // assumes space- or tab-termination
  while ((*source != ' ') && (*source != '\t')) {
    *target++ = *source++;
  }
  return target;
}

inline char* read_next_unsafe_upd(char* target, char** source_ptr) {
  // assumes space- or tab-termination
  while ((**source_ptr != ' ') && (**source_ptr != '\t')) {
    *target++ = **source_ptr;
    *source_ptr += 1;
  }
  return target;
}

inline char* read_next_upd(char* target, char** source_ptr) {
  while (!is_space_or_eoln(**source_ptr)) {
    *target++ = **source_ptr;
    *source_ptr += 1;
  }
  return target;
}

inline int32_t is_contained(char* id_buf, char* lptr, int32_t max_id_len, int32_t filter_line_ct, char* fam_id, char* indiv_id) {
  return (bsearch_fam_indiv(id_buf, lptr, max_id_len, filter_line_ct, fam_id, indiv_id) != -1);
}

void random_thin_markers(double thin_keep_prob, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t removed_ct = 0;
  uint32_t uint32_thresh = (uint32_t)(thin_keep_prob * 4294967296.0 + 0.5);
  uint32_t marker_uidx_stop;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      if (sfmt_genrand_uint32(&sfmt) >= uint32_thresh) {
	SET_BIT(marker_exclude, marker_uidx);
	removed_ct++;
      }
    } while (++marker_uidx < marker_uidx_stop);
  }
  sprintf(logbuf, "--thin: %u marker%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", marker_ct - removed_ct);
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
}


double calc_wt_mean(double exponent, int32_t lhi, int32_t lli, int32_t hhi) {
  double lcount = (double)lli + ((double)lhi * 0.5);
  int64_t tot = lhi + lli + hhi;
  double dtot = (double)tot;
  int64_t subcount = lli; // avoid 32-bit integer overflow
  double weight;
  if ((!lhi) && ((!lli) || (!hhi))) {
    return 0.0;
  }
  weight = pow(2 * lcount * (dtot - lcount) / (dtot * dtot), -exponent);
  subcount = lhi * (subcount + hhi) + 2 * subcount * hhi;
  return (subcount * weight * 2) / (double)(tot * tot);
}

double calc_wt_mean_maf(double exponent, double maf) {
  // assume Hardy-Weinberg equilibrium
  // homozygote frequencies: maf^2, (1-maf)^2
  // heterozygote frequency: 2maf(1-maf)
  double ll_freq = maf * maf;
  double lh_freq = 2 * maf * (1.0 - maf);
  double hh_freq = (1.0 - maf) * (1.0 - maf);
  double weight;
  if (lh_freq == 0.0) {
    return 0.0;
  }
  weight = pow(lh_freq, -exponent);
  return (lh_freq * (ll_freq + lh_freq) + 2 * ll_freq * hh_freq) * weight;
}

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  unsigned char* wkspace_mark;
  unsigned char* wkspace_mark2;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  int32_t oo;
  uint32_t uii;
  uint32_t ujj;
  uint32_t initial_family_blocks;
  uintptr_t indiv_uidx;
  int64_t llii;
  char* family_ids;
  char* cur_person_id;
  char* last_family_id = NULL;
  char* cur_family_id;
  uint32_t* family_sizes;
  uint32_t* uiptr;
  uint32_t* uiptr2 = NULL;
  uint32_t fidx;
  int32_t family_size;
  uint32_t* remaining_indiv_idxs;
  int32_t* remaining_indiv_parent_idxs; // -1 = no parent (or nonshared)
  uint32_t remaining_indiv_ct;
  uint32_t indiv_idx_write;
  uintptr_t max_family_id_len = 0;
  char* indiv_ids;
  uintptr_t max_indiv_id_len = 0;
  int32_t max_pm_id_len;
  uint32_t family_id_ct;
  uint32_t* fis_ptr;
  char* stray_parent_ids;
  int32_t stray_parent_ct;
  uintptr_t* processed_indivs;
  int32_t founder_ct;
  int32_t max_family_nf = 0;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctlm = unfiltered_indiv_ctl * BITCT;
  uint32_t* complete_indiv_idxs;
  uintptr_t complete_indiv_idx_ct;
  double* rs_ptr;
  double* rel_writer;
  double dxx;
  double* tmp_rel_space = NULL;
  double* tmp_rel_writer = NULL;

  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    ujj = strlen_se(&(person_ids[indiv_uidx * max_person_id_len])) + 1;
    if (ujj > max_family_id_len) {
      max_family_id_len = ujj;
    }
    ujj = strlen(&(person_ids[indiv_uidx * max_person_id_len + ujj]));
    if (ujj >= max_indiv_id_len) {
      max_indiv_id_len = ujj + 1;
    }
  }
  if (max_paternal_id_len > max_maternal_id_len) {
    max_pm_id_len = max_paternal_id_len;
  } else {
    max_pm_id_len = max_maternal_id_len;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_space), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_rel_nf_idxs), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_idxs), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&family_ids, unfiltered_indiv_ct * max_family_id_len) ||
      wkspace_alloc_ui_checked(&family_sizes, unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }

  // copy all the items over in order, then qsort, then eliminate duplicates
  // and count family sizes.
  cur_family_id = family_ids;
  cur_person_id = person_ids;
  uiptr = family_sizes;
  *uiptr = 1;
  jj = strlen_se(cur_person_id);
  memcpyx(cur_family_id, cur_person_id, jj, 0);
  for (indiv_uidx = 1; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    cur_person_id = &(cur_person_id[max_person_id_len]);
    mm = strlen_se(cur_person_id);
    if ((jj != mm) || memcmp(cur_family_id, cur_person_id, mm)) {
      cur_family_id = &(cur_family_id[max_family_id_len]);
      memcpyx(cur_family_id, cur_person_id, mm, 0);
      jj = mm;
      *(++uiptr) = 1;
    } else {
      *uiptr += 1;
    }
  }
  initial_family_blocks = 1 + (uint32_t)(uiptr - family_sizes);
  if (qsort_ext(family_ids, initial_family_blocks, max_family_id_len, strcmp_deref, (char*)family_sizes, sizeof(int32_t))) {
    return RET_NOMEM;
  }

  last_family_id = family_ids;
  cur_family_id = &(family_ids[max_family_id_len]);
  family_id_ct = 1;
  uii = 1; // read idx
  if (initial_family_blocks != 1) {
    uiptr = family_sizes;
    while (strcmp(cur_family_id, last_family_id)) {
      family_id_ct++;
      uiptr++;
      last_family_id = cur_family_id;
      cur_family_id = &(cur_family_id[max_family_id_len]);
      uii++;
      if (uii == initial_family_blocks) {
	break;
      }
    }
    if (uii < initial_family_blocks) {
      uiptr2 = uiptr; // family_sizes read pointer
      *uiptr += *(++uiptr2);
      uii++;
      cur_family_id = &(cur_family_id[max_family_id_len]); // read pointer
      while (uii < initial_family_blocks) {
	while (!strcmp(cur_family_id, last_family_id)) {
	  *uiptr += *(++uiptr2);
	  uii++;
	  if (uii == initial_family_blocks) {
	    break;
	  }
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
	if (uii < initial_family_blocks) {
	  *(++uiptr) = *(++uiptr2);
	  last_family_id = &(last_family_id[max_family_id_len]);
	  strcpy(last_family_id, cur_family_id);
	  family_id_ct++;
	  uii++;
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
      }
    }
  }

  uiptr = family_sizes;
  if (family_id_ct < unfiltered_indiv_ct) {
    wkspace_reset(family_ids);
    family_ids = (char*)wkspace_alloc(family_id_ct * max_family_id_len);
    family_sizes = (uint32_t*)wkspace_alloc(family_id_ct * sizeof(int32_t));
    for (uii = 0; uii < family_id_ct; uii++) {
      family_sizes[uii] = *uiptr++;
    }
  }
  pri_ptr->family_ids = family_ids;
  pri_ptr->family_id_ct = family_id_ct;
  pri_ptr->max_family_id_len = max_family_id_len;
  pri_ptr->family_sizes = family_sizes;

  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_offsets), (family_id_ct + 1) * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_rel_space_offsets), (family_id_ct + 1) * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_founder_cts), family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_int_zero((int32_t*)(pri_ptr->family_founder_cts), family_id_ct);

  ii = 0; // running family_info offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_info_offsets[fidx] = ii;
    ii += family_size;
  }

  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  memcpy(uiptr, pri_ptr->family_info_offsets, family_id_ct * sizeof(int32_t));

  // Fill family_idxs, family_founder_cts, and founder portion of
  // family_rel_nf_idxs.
  cur_person_id = person_ids;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    jj = strlen_se(cur_person_id);
    memcpyx(tbuf, cur_person_id, jj, 0);
    kk = bsearch_str(tbuf, family_ids, max_family_id_len, 0, family_id_ct - 1);
    pri_ptr->family_idxs[indiv_uidx] = kk;
    if (IS_SET(founder_info, indiv_uidx)) {
      pri_ptr->family_founder_cts[kk] += 1;
      pri_ptr->family_rel_nf_idxs[indiv_uidx] = uiptr[kk];
      uiptr[kk] += 1;
    }
    cur_person_id = &(cur_person_id[max_person_id_len]);
  }
  wkspace_reset(uiptr);
  jj = 0; // running rel_space offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_rel_space_offsets[fidx] = jj;
    kk = pri_ptr->family_founder_cts[fidx];
    if (family_size - kk > max_family_nf) {
      max_family_nf = family_size - kk;
    }
    // No need to explicitly store the (kk * (kk - 1)) / 2 founder-founder
    // relationships.
    jj += ((int64_t)family_size * (family_size - 1) - (int64_t)kk * (kk - 1)) / 2;
  }
  // make it safe to determine size of blocks by subtracting from the next
  // offset, even if we're at the last family
  pri_ptr->family_info_offsets[family_id_ct] = unfiltered_indiv_ct;
  pri_ptr->family_rel_space_offsets[family_id_ct] = jj;
  if (wkspace_alloc_d_checked(&(pri_ptr->rel_space), jj * sizeof(double))) {
    return RET_NOMEM;
  }

  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  // populate family_info_space
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    uiptr[fidx] = pri_ptr->family_info_offsets[fidx];
  }
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    fidx = pri_ptr->family_idxs[indiv_uidx];
    pri_ptr->family_info_space[uiptr[fidx]] = indiv_uidx;
    uiptr[fidx] += 1;
  }
  wkspace_reset(wkspace_mark);

  if (wkspace_alloc_ul_checked(&processed_indivs, (unfiltered_indiv_ctl + (max_family_nf + (BITCT2 - 1)) / BITCT2) * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_one(&(processed_indivs[unfiltered_indiv_ctl]), (max_family_nf + (BITCT2 - 1)) / BITCT2);

  wkspace_mark2 = wkspace_base;
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = family_sizes[fidx];
    founder_ct = pri_ptr->family_founder_cts[fidx];
    remaining_indiv_ct = family_size - founder_ct;
    stray_parent_ct = 0;
    if (remaining_indiv_ct) {
      memcpy(processed_indivs, founder_info, unfiltered_indiv_ctl * sizeof(intptr_t));
      if (wkspace_alloc_ui_checked(&complete_indiv_idxs, family_size * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&remaining_indiv_idxs, remaining_indiv_ct * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&indiv_ids, family_size * max_indiv_id_len) ||
          wkspace_alloc_i_checked(&remaining_indiv_parent_idxs, remaining_indiv_ct * 2 * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&stray_parent_ids, remaining_indiv_ct * 2 * max_pm_id_len)) {
	return RET_NOMEM;
      }
      ii = pri_ptr->family_info_offsets[fidx];
      fis_ptr = &(pri_ptr->family_info_space[ii]);
      rs_ptr = &(pri_ptr->rel_space[pri_ptr->family_rel_space_offsets[fidx]]);
      rel_writer = rs_ptr;
      cur_person_id = indiv_ids;
      for (ii = 0; ii < family_size; ii++) {
	kk = fis_ptr[ii];
	jj = strlen_se(&(person_ids[kk * max_person_id_len])) + 1;
	strcpy(cur_person_id, &(person_ids[kk * max_person_id_len + jj]));
	cur_person_id = &(cur_person_id[max_indiv_id_len]);
      }

      // Compile list of non-founder family member indices, and identify
      // parents who are referred to by individual ID but are NOT in the
      // dataset.
      ii = 0; // family_info_space index
      complete_indiv_idx_ct = 0;
      cur_person_id = stray_parent_ids;
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	while (IS_SET(founder_info, fis_ptr[ii])) {
	  complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[ii];
	  ii++;
	}
	kk = fis_ptr[ii++];

	// does not track sex for now
	if (memcmp("0", &(paternal_ids[kk * max_paternal_id_len]), 2)) {
	  mm = bsearch_str(&(paternal_ids[kk * max_paternal_id_len]), indiv_ids, max_indiv_id_len, 0, family_size - 1);
	  if (mm == -1) {
	    strcpy(cur_person_id, &(paternal_ids[kk * max_paternal_id_len]));
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2] = -2;
	  } else {
            remaining_indiv_parent_idxs[uii * 2] = fis_ptr[mm];
	  }
	} else {
          remaining_indiv_parent_idxs[uii * 2] = -1;
	}
	if (memcmp("0", &(maternal_ids[kk * max_maternal_id_len]), 2)) {
	  mm = bsearch_str(&(maternal_ids[kk * max_maternal_id_len]), indiv_ids, max_indiv_id_len, 0, family_size - 1);
	  if (mm == -1) {
	    strcpy(cur_person_id, &(maternal_ids[kk * max_maternal_id_len]));
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2 + 1] = -2;
	  } else {
	    remaining_indiv_parent_idxs[uii * 2 + 1] = fis_ptr[mm];
	  }
	} else {
	  remaining_indiv_parent_idxs[uii * 2 + 1] = -1;
	}
        remaining_indiv_idxs[uii] = kk;
      }
      while (ii < family_size) {
	complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[ii];
	ii++;
      }
      qsort(stray_parent_ids, stray_parent_ct, max_pm_id_len, strcmp_casted);
      cur_person_id = stray_parent_ids;
      ii = 0; // read idx
      jj = 0; // write idx

      // Now filter out all such parents who aren't referenced at least twice.
      while (ii < stray_parent_ct - 1) {
        if (strcmp(&(stray_parent_ids[ii * max_pm_id_len]), &(stray_parent_ids[(ii + 1) * max_pm_id_len]))) {
	  ii++;
	  continue;
	}
	ii++;
	strcpy(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len]));
	do {
	  ii++;
        } while (!(strcmp(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len])) || (ii >= stray_parent_ct - 1)));
        cur_person_id = &(cur_person_id[max_pm_id_len]);
	jj++;
      }
      stray_parent_ct = jj;

      // Now allocate temporary relatedness table between nonfounders and
      // stray parents with multiple references.
      if (stray_parent_ct) {
        if (wkspace_alloc_d_checked(&tmp_rel_space, (family_size - founder_ct) * stray_parent_ct * sizeof(double))) {
	  return RET_NOMEM;
        }
	tmp_rel_writer = tmp_rel_space;
      }

      // Now fill in remainder of remaining_indiv_parent_idxs.
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	jj = remaining_indiv_idxs[uii];
	if (remaining_indiv_parent_idxs[uii * 2] == -2) {
	  kk = bsearch_str(&(paternal_ids[jj * max_paternal_id_len]), stray_parent_ids, max_pm_id_len, 0, stray_parent_ct - 1);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2] = kk;
	}
	if (remaining_indiv_parent_idxs[uii * 2 + 1] == -2) {
	  kk = bsearch_str(&(maternal_ids[jj * max_maternal_id_len]), stray_parent_ids, max_pm_id_len, 0, stray_parent_ct - 1);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2 + 1] = kk;
	}
      }
      llii = (int64_t)founder_ct * (founder_ct - 1);
      while (remaining_indiv_ct) {
	indiv_idx_write = 0;
	for (uii = 0; uii < remaining_indiv_ct; uii++) {
	  kk = remaining_indiv_parent_idxs[uii * 2];
	  mm = remaining_indiv_parent_idxs[uii * 2 + 1];
	  jj = remaining_indiv_idxs[uii];
	  if (((kk == -1) || is_set(processed_indivs, kk)) && ((mm == -1) || is_set(processed_indivs, mm))) {
	    for (nn = 0; nn < founder_ct; nn++) {
	      // relationship between kk and nnth founder
	      if ((kk >= (int32_t)unfiltered_indiv_ct) || (kk == -1)) {
		dxx = 0.0;
	      } else if (is_set(founder_info, kk)) {
		if (kk == (int32_t)complete_indiv_idxs[nn]) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
                dxx = 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
	      }
	      if (is_set(founder_info, mm)) {
		if (mm == (int32_t)complete_indiv_idxs[nn]) {
		  dxx += 0.5;
		}
	      } else if ((mm != -1) && (mm < (int32_t)unfiltered_indiv_ct)) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		dxx += 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
	      }
	      *rel_writer++ = dxx;
	    }
	    for (; nn < (int32_t)complete_indiv_idx_ct; nn++) {
	      if (kk == -1) {
		dxx = 0.0;
	      } else if (kk >= (int32_t)unfiltered_indiv_ct) {
		dxx = 0.5 * tmp_rel_space[(nn - founder_ct) * stray_parent_ct + kk - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, kk)) {
                dxx = 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + pri_ptr->family_rel_nf_idxs[kk]];
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
		if (oo == nn) {
		  dxx = 0.5;
		} else if (oo < nn) {
		  dxx = 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + oo];
		} else {
		  dxx = 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
		}
	      }
	      if (mm >= (int32_t)unfiltered_indiv_ct) {
		dxx += 0.5 * tmp_rel_space[(nn - founder_ct) * stray_parent_ct + mm - unfiltered_indiv_ctlm];
	      } else if (is_set(founder_info, mm)) {
		dxx += 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + pri_ptr->family_rel_nf_idxs[mm]];
	      } else if (mm != -1) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		if (oo == nn) {
		  dxx += 0.5;
		} else if (oo < nn) {
		  dxx += 0.5 * rs_ptr[((int64_t)nn * (nn - 1) - llii) / 2 + oo];
		} else {
		  dxx += 0.5 * rs_ptr[((int64_t)oo * (oo - 1) - llii) / 2 + nn];
		}
	      }
	      *rel_writer++ = dxx;
	    }
	    for (nn = 0; nn < stray_parent_ct; nn++) {
	      if (kk >= (int32_t)unfiltered_indiv_ct) {
		if (kk == nn + (int32_t)unfiltered_indiv_ctlm) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else if (kk == -1) {
                dxx = 0.0;
	      } else {
		oo = pri_ptr->family_rel_nf_idxs[kk];
		if (oo < founder_ct) {
		  dxx = 0.0;
		} else {
                  dxx = 0.5 * tmp_rel_space[(oo - founder_ct) * stray_parent_ct + nn];
		}
	      }
	      if (mm >= (int32_t)unfiltered_indiv_ct) {
		if (mm == nn + (int32_t)unfiltered_indiv_ctlm) {
		  dxx += 0.5;
		}
	      } else if (mm != -1) {
		oo = pri_ptr->family_rel_nf_idxs[mm];
		if (oo >= founder_ct) {
		  dxx += 0.5 * tmp_rel_space[(oo - founder_ct) * stray_parent_ct + nn];
		}
	      }
	      *tmp_rel_writer++ = dxx;
	    }
	    pri_ptr->family_rel_nf_idxs[jj] = complete_indiv_idx_ct;
	    complete_indiv_idxs[complete_indiv_idx_ct++] = jj;
	    set_bit(processed_indivs, jj);
	  } else {
            remaining_indiv_parent_idxs[indiv_idx_write * 2] = kk;
	    remaining_indiv_parent_idxs[indiv_idx_write * 2 + 1] = mm;
	    remaining_indiv_idxs[indiv_idx_write++] = jj;
	  }
	}
	if (indiv_idx_write == remaining_indiv_ct) {
	  logprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
	  return RET_INVALID_FORMAT;
	}
	remaining_indiv_ct = indiv_idx_write;
      }
      wkspace_reset(wkspace_mark2);
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t make_founders(uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t require_two, uintptr_t* indiv_exclude, uintptr_t* founder_info) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uint32_t new_founder_ct = 0;
  int32_t retval = 0;
  char* sorted_ids;
  char* id_buf;
  char* wptr;
  char* pat_ptr;
  char* mat_ptr;
  uintptr_t* nf_bitarr;
  uintptr_t indiv_uidx;
  uint32_t fam_len_p1;
  uint32_t missing_parent_ct;
  uint32_t cur_len;
  int32_t ii;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&nf_bitarr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto make_founders_ret_NOMEM;
  }
  bitfield_exclude_to_include(indiv_exclude, nf_bitarr, unfiltered_indiv_ct);
  bitfield_andnot(nf_bitarr, founder_info, unfiltered_indiv_ctl);
  indiv_uidx = next_set(nf_bitarr, 0, unfiltered_indiv_ct);
  if (indiv_uidx == unfiltered_indiv_ct) {
    logprint("Note: Skipping --make-founders since there are no nonfounders.\n");
    goto make_founders_ret_1;
  }
  sorted_ids = alloc_and_init_collapsed_arr(person_ids, max_person_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_ct, 0);
  if (!sorted_ids) {
    goto make_founders_ret_NOMEM;
  }
  qsort(sorted_ids, indiv_ct, max_person_id_len, strcmp_casted);
  do {
    pat_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    fam_len_p1 = strlen_se(pat_ptr) + 1;
    wptr = memcpya(id_buf, pat_ptr, fam_len_p1);
    missing_parent_ct = 0;
    pat_ptr = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
    cur_len = strlen(pat_ptr);
    if (cur_len + fam_len_p1 >= max_person_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, pat_ptr, cur_len + 1);
      ii = bsearch_str(id_buf, sorted_ids, max_person_id_len, 0, indiv_ct - 1);
      if (ii == -1) {
	missing_parent_ct++;
      }
    }
    mat_ptr = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
    cur_len = strlen(mat_ptr);
    if (cur_len + fam_len_p1 >= max_person_id_len) {
      missing_parent_ct++;
    } else {
      memcpy(wptr, mat_ptr, cur_len + 1);
      ii = bsearch_str(id_buf, sorted_ids, max_person_id_len, 0, indiv_ct - 1);
      if (ii == -1) {
	missing_parent_ct++;
      }
    }
    if (missing_parent_ct > require_two) {
      SET_BIT(founder_info, indiv_uidx);
      memcpy(pat_ptr, "0", 2);
      memcpy(mat_ptr, "0", 2);
      new_founder_ct++;
    }
    indiv_uidx++;
    next_set_ul_ck(nf_bitarr, &indiv_uidx, unfiltered_indiv_ct);
  } while (indiv_uidx < unfiltered_indiv_ct);
  sprintf(logbuf, "--make-founders: %u individual%s affected.\n", new_founder_ct, (new_founder_ct == 1)? "" : "s");
  logprintb();
  while (0) {
  make_founders_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
 make_founders_ret_1:
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_indiv_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr) {
  uint32_t mp_strlen = strlen(makepheno_str);
  uint32_t makepheno_all = ((mp_strlen == 1) && (makepheno_str[0] == '*'));
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* pheno_c = *pheno_c_ptr;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  int32_t ii;
  uint32_t person_idx;
  uint32_t tmp_len;
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  if (!pheno_c) {
    pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
    *pheno_c_ptr = pheno_c;
  }
  if (makepheno_all) {
    fill_ulong_one(pheno_nm, unfiltered_indiv_ctl);
    zero_trailing_bits(pheno_nm, unfiltered_indiv_ct);
  }
  tbuf[MAXLINELEN - 1] = ' ';  
  while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    bufptr0 = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    bufptr = next_item(bufptr0);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_phenotype_format);
      return RET_INVALID_FORMAT;
    }
    ii = bsearch_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, unfiltered_indiv_ct, bufptr0, bufptr);
    if (ii != -1) {
      person_idx = id_map[(uint32_t)ii];
      if (makepheno_all) {
	SET_BIT(pheno_c, person_idx);
      } else {
	SET_BIT(pheno_nm, person_idx);
	bufptr = next_item(bufptr);
        tmp_len = strlen_se(bufptr);
	if ((tmp_len == mp_strlen) && (!memcmp(bufptr, makepheno_str, mp_strlen))) {
	  SET_BIT(pheno_c, person_idx);
	}
      }
    }
  }
  if (!feof(phenofile)) {
    return RET_READ_FAIL;
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t convert_tail_pheno(uint32_t unfiltered_indiv_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod) {
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  uint32_t indiv_uidx;
  uint32_t indiv_uidx_stop;
  double dxx;
  if (!(*pheno_d_ptr)) {
    logprint("Error: --tail-pheno requires scalar phenotype data.\n");
    return RET_INVALID_FORMAT;
  }
  indiv_uidx = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (!pheno_c) {
    pheno_c = (uintptr_t*)malloc(indiv_uidx * sizeof(intptr_t));
    if (!pheno_c) {
      return RET_NOMEM;
    }
    *pheno_c_ptr = pheno_c;
  }
  fill_ulong_zero(pheno_c, indiv_uidx);
  indiv_uidx = 0;
  do {
    indiv_uidx = next_set(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    indiv_uidx_stop = next_unset(pheno_nm, indiv_uidx, unfiltered_indiv_ct);
    for (; indiv_uidx < indiv_uidx_stop; indiv_uidx++) {
      dxx = pheno_d[indiv_uidx];
      if (dxx > tail_bottom) {
        if (dxx > tail_top) {
          SET_BIT(pheno_c, indiv_uidx);
        } else {
	  CLEAR_BIT(pheno_nm, indiv_uidx);
        }
      }
    }
  } while (indiv_uidx_stop < unfiltered_indiv_ct);
  free(pheno_d);
  *pheno_d_ptr = NULL;
  return 0;
}

const unsigned char acgt_reverse_arr[] = "1B2DEF3HIJKLMNOPQRS4";
const unsigned char acgt_arr[] = "ACGT";
// g_one_char_strs offsets (double)
const unsigned char acgt_reverse_arr1[] = "b\204d\210\212\214f\220\222\224\226\230\232\234\236\240\242\244\246h";
const unsigned char acgt_arr1[] = "\202\206\216\250";
// const unsigned char acgt_reverse_arr1[] = "\"D$HJL&PRTVXZ\\^`bdf(";
// const unsigned char acgt_arr1[] = "BFNh";

static inline unsigned char conditional_convert(unsigned char diff, unsigned char ucc2_max, const unsigned char* convert_arr, unsigned char ucc) {
  unsigned char ucc2 = ucc - diff;
  return (ucc2 < ucc2_max)? convert_arr[ucc2] : ucc;
}

static inline void conditional_convert1(unsigned char diff, unsigned char ucc2_max, const unsigned char* convert_arr, char** allele_ptr) {
  unsigned char ucc2 = ((unsigned char)(**allele_ptr)) - diff;
  if (ucc2 < ucc2_max) {
    *allele_ptr = (char*)(&(g_one_char_strs[convert_arr[ucc2]]));
  }
}

void allelexxxx_recode(uint32_t allelexxxx, char** marker_allele_ptrs, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t recode_multichar = allelexxxx & ALLELE_RECODE_MULTICHAR;
  const unsigned char* convert_arr;
  const unsigned char* convert_arr1;
  char** map_ptr;
  char** map_ptr_stop;
  char* cptr;
  uint32_t marker_uidx_stop;
  unsigned char diff;
  unsigned char ucc2_max;
  unsigned char ucc;
  if (allelexxxx & ALLELE_RECODE_ACGT) {
    diff = 49;
    ucc2_max = 4;
    convert_arr = acgt_arr;
    convert_arr1 = acgt_arr1;
  } else {
    diff = 65;
    ucc2_max = 20;
    convert_arr = acgt_reverse_arr;
    convert_arr1 = acgt_reverse_arr1;
  }
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    map_ptr = &(marker_allele_ptrs[marker_uidx * 2]);
    map_ptr_stop = &(marker_allele_ptrs[marker_uidx_stop * 2]);
    marker_uidx = marker_uidx_stop;
    if (recode_multichar) {
      do {
	cptr = *map_ptr;
	if (!cptr[1]) {
	  conditional_convert1(diff, ucc2_max, convert_arr1, map_ptr);
	} else {
	  ucc = *cptr;
	  do {
	    *cptr = conditional_convert(diff, ucc2_max, convert_arr, ucc);
	    ucc = *(++cptr);
	  } while (ucc);
	}
      } while (++map_ptr < map_ptr_stop);
    } else {
      do {
        if (!(map_ptr[0][1])) {
          conditional_convert1(diff, ucc2_max, convert_arr1, map_ptr);
	}
      } while (++map_ptr < map_ptr_stop);
    }
  }
}

int32_t filter_indivs_file(char* filtername, char* sorted_person_ids, uintptr_t sorted_ids_len, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* filtervals_flattened, uint32_t mfilter_col) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t include_ct = 0;
  uintptr_t max_filterval_len = 0;
  uint32_t filterval_ct = 0;
  int32_t retval = 0;
  char* sorted_filtervals;
  uintptr_t* indiv_exclude_new;
  char* id_buf;
  char* bufptr;
  uint32_t filterval_idx;
  uint32_t slen;
  int32_t person_idx;
  bufptr = filtervals_flattened;
  do {
    filterval_ct++;
    slen = strlen(bufptr);
    if (slen >= max_filterval_len) {
      max_filterval_len = slen + 1;
    }
    bufptr = &(bufptr[slen + 1]);
  } while (*bufptr);
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&indiv_exclude_new, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_c_checked(&sorted_filtervals, filterval_ct * max_filterval_len)) {
    goto filter_indivs_file_ret_NOMEM;
  }
  fill_ulong_one(indiv_exclude_new, unfiltered_indiv_ctl);
  bufptr = filtervals_flattened;
  for (filterval_idx = 0; filterval_idx < filterval_ct; filterval_idx++) {
    slen = strlen(bufptr) + 1;
    memcpy(&(sorted_filtervals[filterval_idx * max_filterval_len]), bufptr, slen);
    bufptr = &(bufptr[slen]);
  }
  qsort(sorted_filtervals, filterval_ct, max_filterval_len, strcmp_casted);

  if (fopen_checked(&infile, filtername, "r")) {
    goto filter_indivs_file_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Excessively long line in --keep/--remove file (max %u chars).\n", MAXLINELEN - 3);
      goto filter_indivs_file_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      goto filter_indivs_file_ret_INVALID_FORMAT;
    }
    person_idx = bsearch_fam_indiv(id_buf, sorted_person_ids, max_person_id_len, sorted_ids_len, tbuf, bufptr);
    if (person_idx != -1) {
      person_idx = id_map[(uint32_t)person_idx];
      if (!is_set(indiv_exclude, person_idx)) {
	bufptr = next_item_mult(bufptr, mfilter_col);
	if (no_more_items_kns(bufptr)) {
	  goto filter_indivs_file_ret_INVALID_FORMAT;
	}
	bufptr[strlen_se(bufptr)] = '\0';
	if (bsearch_str(bufptr, sorted_filtervals, max_filterval_len, 0, filterval_ct - 1) != -1) {
	  if (is_set(indiv_exclude_new, person_idx)) {
	    clear_bit(indiv_exclude_new, person_idx);
	    include_ct++;
	  }
	}
      }
    }
  }
  if (!feof(infile)) {
    goto filter_indivs_file_ret_READ_FAIL;
  }
  memcpy(indiv_exclude, indiv_exclude_new, unfiltered_indiv_ctl * sizeof(intptr_t));
  *indiv_exclude_ct_ptr = unfiltered_indiv_ct - include_ct;

  while (0) {
  filter_indivs_file_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_indivs_file_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  filter_indivs_file_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_indivs_file_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  filter_indivs_file_ret_INVALID_FORMAT:
    logprint("Error: Too few columns in --filter file line.\n");
    retval = RET_INVALID_FORMAT;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

void filter_indivs_bitfields(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* orfield, int32_t orfield_flip, uintptr_t* ornot) {
  // indiv_exclude := indiv_exclude | orfield | (~ornot) if !orfield_flip
  //               := indiv_exclude | (~orfield) | (~ornot) otherwise
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t* ieptr = indiv_exclude;
  uintptr_t* ieend = &(indiv_exclude[unfiltered_indiv_ctl]);
  if (orfield_flip) {
    if (ornot) {
      do {
	*ieptr |= (~(*orfield++)) | (~(*ornot++));
      } while (++ieptr < ieend);
    } else {
      do {
	*ieptr |= ~(*orfield++);
      } while (++ieptr < ieend);
    }
  } else {
    if (ornot) {
      do {
	*ieptr |= (*orfield++) | (~(*ornot++));
      } while (++ieptr < ieend);
    } else {
      do {
	*ieptr |= (*orfield++) | (~(*ornot++));
      } while (++ieptr < ieend);
    }
  }
  zero_trailing_bits(indiv_exclude, unfiltered_indiv_ct);
  *indiv_exclude_ct_ptr = popcount_longs(indiv_exclude, 0, unfiltered_indiv_ctl);
}

void calc_plink_maxfid(uint32_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uint32_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t* plink_maxfid_ptr, uint32_t* plink_maxiid_ptr) {
  uintptr_t plink_maxfid = 4;
  uintptr_t plink_maxiid = 4;
  uint32_t indiv_uidx = 0;
  uint32_t indivs_done = 0;
  char* cptr;
  char* cptr2;
  char* cptr_end;
  uintptr_t slen;
  uint32_t indiv_uidx_stop;
  // imitate PLINK behavior (see Plink::prettyPrintLengths() in helper.cpp), to
  // simplify testing and avoid randomly breaking existing scripts
  do {
    indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
    indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
    indivs_done += indiv_uidx_stop - indiv_uidx;
    cptr = &(person_ids[indiv_uidx * max_person_id_len]);
    cptr_end = &(person_ids[indiv_uidx_stop * max_person_id_len]);
    indiv_uidx = indiv_uidx_stop;
    do {
      cptr2 = (char*)memchr(cptr, '\t', max_person_id_len);
      slen = (uintptr_t)(cptr2 - cptr);
      if (slen > plink_maxfid) {
	plink_maxfid = slen + 2;
      }
      slen = strlen(&(cptr2[1]));
      if (slen > plink_maxiid) {
        plink_maxiid = slen + 2;
      }
      cptr = &(cptr[max_person_id_len]);
    } while (cptr < cptr_end);
  } while (indivs_done < indiv_ct);
  *plink_maxfid_ptr = plink_maxfid;
  *plink_maxiid_ptr = plink_maxiid;
}

uint32_t calc_plink_maxsnp(uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len) {
  uintptr_t plink_maxsnp = 4;
  uintptr_t max_marker_id_len_m1 = max_marker_id_len - 1;
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  char* cptr;
  char* cptr_end;
  uintptr_t slen;
  uint32_t marker_uidx_stop;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
    cptr_end = &(marker_ids[marker_uidx_stop * max_marker_id_len]);
    marker_uidx = marker_uidx_stop;
    do {
      slen = strlen(cptr);
      if (slen > plink_maxsnp) {
	plink_maxsnp = slen + 2;
	if (plink_maxsnp >= max_marker_id_len_m1) {
	  return plink_maxsnp;
	}
      }
      cptr = &(cptr[max_marker_id_len]);
    } while (cptr < cptr_end);
  }
  return plink_maxsnp;
}

int32_t mind_filter(FILE* bedfile, uintptr_t bed_offset, double mind_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uint32_t mind_int_thresh = (int32_t)(mind_thresh * ((int32_t)marker_ct) + EPSILON);
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ct2l = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t marker_idx = 0;
  uint32_t indiv_exclude_ct = *indiv_exclude_ct_ptr;
  uint32_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  uint32_t indiv_uidx = 0;
  uint32_t indiv_idx = 0;
  uint32_t removed_ct = 0;
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* lptr;
  uint32_t* missing_cts;
  uintptr_t marker_uidx;
  uint32_t indiv_uidx_stop;
  uint32_t uii;
  uint32_t ujj;
  uintptr_t ulii;

  if (wkspace_alloc_ui_checked(&missing_cts, unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ct2l * sizeof(intptr_t))) {
    goto mind_filter_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ct2l - 1] = 0;
  fill_uint_zero(missing_cts, unfiltered_indiv_ct);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto mind_filter_ret_READ_FAIL;
  }
  ujj = unfiltered_indiv_ct2l * BITCT2;
  for (marker_uidx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto mind_filter_ret_READ_FAIL;
      }
    }
    if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
      goto mind_filter_ret_READ_FAIL;
    }
    lptr = loadbuf;
    for (uii = 0; uii < ujj; uii += BITCT2) {
      ulii = *lptr++;
      ulii = (ulii & FIVEMASK) & ((~ulii) >> 1);
      // now ulii has single bit set only at missing positions
      while (ulii) {
	missing_cts[uii + CTZLU(ulii) / 2] += 1;
	ulii &= ulii - 1;
      }
    }
  }
  do {
    indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
    indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    do {
      if (missing_cts[indiv_uidx] > mind_int_thresh) {
        SET_BIT(indiv_exclude, indiv_uidx);
	removed_ct++;
      }
    } while (++indiv_uidx < indiv_uidx_stop);
  } while (indiv_idx < indiv_ct);
  *indiv_exclude_ct_ptr += removed_ct;
  sprintf(logbuf, "%u %s removed due to missing genotype data (--mind).\n", removed_ct, species_str(removed_ct));
  logprintb();
  while (0) {
  mind_filter_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  mind_filter_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

#ifdef __LP64__
void freq_hwe_haploid_count_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct_nm1;
  __m128i to_ct_hmaj1;
  __m128i to_ct_nm2;
  __m128i to_ct_hmaj2;
  __uni16 acc_nm;
  __uni16 acc_hmaj;

  acc_nm.vi = _mm_setzero_si128();
  acc_hmaj.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj1 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(_mm_and_si128(to_ct_nm1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 2), m2));
    to_ct_hmaj1 = _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 2), m2));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj2 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_add_epi64(_mm_and_si128(to_ct_nm2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm2, 2), m2)));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_add_epi64(_mm_and_si128(to_ct_hmaj2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj2, 2), m2)));

    acc_nm.vi = _mm_add_epi64(acc_nm.vi, _mm_add_epi64(_mm_and_si128(to_ct_nm1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 4), m4)));
    acc_hmaj.vi = _mm_add_epi64(acc_hmaj.vi, _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 4), m4)));
  } while (vptr < vend);
  acc_nm.vi = _mm_add_epi64(_mm_and_si128(acc_nm.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_nm.vi, 8), m8));
  acc_hmaj.vi = _mm_add_epi64(_mm_and_si128(acc_hmaj.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_hmaj.vi, 8), m8));
  *ct_nmp += ((acc_nm.u8[0] + acc_nm.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct_hmajp += ((acc_hmaj.u8[0] + acc_hmaj.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
void freq_hwe_haploid_count_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader3 = loader >> 1;
  uintptr_t loader2 = loader ^ (~loader3);
  uint32_t to_ct_nm1;
  uint32_t to_ct_hmaj1;
  uint32_t to_ct_nm2;
  uint32_t to_ct_hmaj2;
  uintptr_t partial_nm;
  uintptr_t partial_hmaj;
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm = (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj = (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm += (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj += (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  *ct_nmp += (partial_nm * 0x01010101) >> 24;
  *ct_hmajp += (partial_hmaj * 0x01010101) >> 24;
}
#endif

static inline void single_marker_freqs_and_hwe(uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* founder_include2, uintptr_t* founder_ctrl_include2, uintptr_t* founder_case_include2, uintptr_t indiv_ct, uint32_t* ll_ctp, uint32_t* lh_ctp, uint32_t* hh_ctp, uint32_t indiv_f_ct, uint32_t* ll_ctfp, uint32_t* lh_ctfp, uint32_t* hh_ctfp, int32_t hwe_needed, uintptr_t indiv_f_ctl_ct, uint32_t* ll_hwep, uint32_t* lh_hwep, uint32_t* hh_hwep, int32_t hardy_needed, uintptr_t indiv_f_case_ct, uint32_t* ll_case_hwep, uint32_t* lh_case_hwep, uint32_t* hh_case_hwep) {
  // This is best understood from the bottom third up (which is the order it
  // was written).  It's way overkill for just determining genotype
  // frequencies, but a ruthlessly optimized version is needed for e.g.
  // permutation testing so we may as well get it working here.
  //
  // The idea, which underlies the IBS and LD pruners as well, is to obtain
  // multiple popcounts at once.  In the case of PLINK frequency/HWE
  // evaluation, we need 6-9 numbers:
  // homozyg minor, het, homozyg major counts for all individuals
  // homozyg minor, het, homozyg major counts for all founders
  // sometimes, homozyg minor, het, homozyg major counts for all ctrl founders
  //
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = indiv_ct - homozyg minor ct - het ct
  //
  // Thus, with the appropriate indiv_ct and these three popcounts, we can
  // calculate a set of genotype counts.  We perform the
  // not-that-exploitative version of these calculations in the bottom third of
  // this function, to deal with the remainder that doesn't fit into the
  // 12-word block size of the main loops.
  //
  // The middle third is a 12-word block popcount for 32-bit platforms (see
  // popcount_longs() in wdist_common.c; this is 12 words instead of 6 since
  // odd bits of the popcount targets are guaranteed to be zero, delaying
  // overflow).  It could be improved a bit, but we care more about reliability
  // than blistering speed for the 32-bit build (since it's used to check the
  // results of the actually blistering 64-bit code...).  Sure, hardware
  // popcount is significantly faster, but what x86-32 machine has hardware
  // popcount??...
  //
  // The top third is the portable Lauradoux/Walisch loop.
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_c = 0;
  uint32_t tot_a_f = 0;
  uint32_t tot_b_f = 0;
  uint32_t tot_c_f = 0;
  uint32_t tot_a_hwe = 0;
  uint32_t tot_b_hwe = 0;
  uint32_t tot_c_hwe = 0;
  uint32_t tot_a_chwe = 0;
  uint32_t tot_b_chwe = 0;
  uint32_t tot_c_chwe = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_c);
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_needed) {
      count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[cur_decr]);
      if (hardy_needed) {
	count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[cur_decr]);
      }
    }
    lptr = lptr_12x_end;
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_c);
    count_3freq_12(lptr, founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_needed) {
      count_3freq_12(lptr, founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[12]);
      if (hardy_needed) {
	count_3freq_12(lptr, founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[12]);
      }
    }
    lptr = &(lptr[12]);
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *indiv_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    // N.B. because of the construction of indiv_include2, only even-numbered
    // bits can be present here.  So popcount2_long is safe.
    tot_a += popcount2_long(loader2);
    tot_b += popcount2_long(loader3);
    tot_c += popcount2_long(loader & loader3);
    loader2 = *founder_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_a_f += popcount2_long(loader2);
    tot_b_f += popcount2_long(loader3);
    tot_c_f += popcount2_long(loader & loader3);
    if (hwe_needed) {
      loader2 = *founder_ctrl_include2++;
      loader3 = (loader >> 1) & loader2;
      loader2 &= loader;
      tot_a_hwe += popcount2_long(loader2);
      tot_b_hwe += popcount2_long(loader3);
      tot_c_hwe += popcount2_long(loader & loader3);
      if (hardy_needed) {
	loader2 = *founder_case_include2++;
	loader3 = (loader >> 1) & loader2;
	loader2 &= loader;
	tot_a_chwe += popcount2_long(loader2);
	tot_b_chwe += popcount2_long(loader3);
	tot_c_chwe += popcount2_long(loader & loader3);
      }
    }
  }
  *hh_ctp = tot_c;
  *lh_ctp = tot_b - tot_c;
  *ll_ctp = indiv_ct - tot_a - *lh_ctp;
  *hh_ctfp = tot_c_f;
  *lh_ctfp = tot_b_f - tot_c_f;
  *ll_ctfp = indiv_f_ct - tot_a_f - *lh_ctfp;
  if (hwe_needed) {
    *hh_hwep = tot_c_hwe;
    *lh_hwep = tot_b_hwe - tot_c_hwe;
    *ll_hwep = indiv_f_ctl_ct - tot_a_hwe - *lh_hwep;
    if (hardy_needed) {
      *hh_case_hwep = tot_c_chwe;
      *lh_case_hwep = tot_b_chwe - tot_c_chwe;
      *ll_case_hwep = indiv_f_case_ct - tot_a_chwe - *lh_case_hwep;
    }
  }
}

static inline uint32_t nonmissing_present_diff(uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* indiv_male_include2) {
  // possible todo: write entries to .ynm file, using same format as .hh
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  do {
    loader = *lptr++;
    loader2 = (*indiv_include2++) & (*indiv_male_include2++);
    // when really bored, check if compiler translates this into andnot
    // operations
    if ((~((~(loader >> 1)) & loader)) & loader2) {
      return 1;
    }
  } while (lptr < lptr_end);
  return 0;
}

static inline void haploid_single_marker_freqs(uintptr_t unfiltered_indiv_ct, uintptr_t unfiltered_indiv_ctl2, uintptr_t* lptr, uintptr_t* indiv_include2, uintptr_t* founder_include2, uintptr_t indiv_ct, uint32_t* ll_ctp, uint32_t* hh_ctp, uint32_t indiv_f_ct, uint32_t* ll_ctfp, uint32_t* hh_ctfp, uint32_t* hethap_incr_ptr) {
  // Here, we interpret heterozygotes as missing.
  // Nonmissing: (genotype ^ (~(genotype >> 1))) & 0x5555...
  // Homozygote major: (genotype & (genotype >> 1)) & 0x5555...
  uint32_t tot_a = 0;
  uint32_t tot_b = 0;
  uint32_t tot_hmaj = 0;
  uint32_t tot_nm_f = 0;
  uint32_t tot_hmaj_f = 0;
  uintptr_t* lptr_end = &(lptr[unfiltered_indiv_ctl2]);
  uint32_t hethap_incr;
  uint32_t tot_nm;
  uint32_t uii;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  uintptr_t* lptr_12x_end;
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
  // Given a buffer with PLINK binary genotypes for a single marker, let
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  // Then,
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  //               = indiv_ct - homozyg minor ct - het ct
    count_3freq_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = lptr_12x_end;
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    count_3freq_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_hmaj);
    freq_hwe_haploid_count_12(lptr, founder_include2, &tot_nm_f, &tot_hmaj_f);
    lptr = &(lptr[12]);
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  tot_nm = 2 * tot_hmaj + indiv_ct - tot_a - tot_b;
  hethap_incr = tot_b - tot_hmaj;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader3 = loader >> 1;
    loader4 = *indiv_include2++;
    // different from tot_nm_f because of +indiv_ct above
    tot_nm -= popcount2_long(loader & loader4);
    hethap_incr += popcount2_long(loader3 & (~loader) & loader4);
    loader2 = loader ^ (~loader3); // nonmissing?
    loader &= loader3; // homozyg A2?
    uii = popcount2_long(loader & loader4);
    tot_nm += 2 * uii - popcount2_long(loader3 & loader4);
    // tot_nm += popcount2_long(loader2 & loader4);
    tot_hmaj += uii;
    loader3 = *founder_include2++;
    tot_nm_f += popcount2_long(loader2 & loader3);
    tot_hmaj_f += popcount2_long(loader & loader3);
  }
  *hh_ctp = tot_hmaj;
  *ll_ctp = tot_nm - tot_hmaj;
  *hh_ctfp = tot_hmaj_f;
  *ll_ctfp = tot_nm_f - tot_hmaj_f;
  *hethap_incr_ptr = hethap_incr;
}

int32_t calc_freqs_and_hwe(FILE* bedfile, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_exclude_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* founder_info, int32_t nonfounders, int32_t maf_succ, double* set_allele_freqs, uintptr_t** marker_reverse_ptr, uint32_t** marker_allele_cts_ptr, uintptr_t bed_offset, uint32_t hwe_needed, uint32_t hwe_all, uint32_t hardy_needed, uintptr_t* pheno_nm, uintptr_t* pheno_c, int32_t** hwe_lls_ptr, int32_t** hwe_lhs_ptr, int32_t** hwe_hhs_ptr, int32_t** hwe_ll_cases_ptr, int32_t** hwe_lh_cases_ptr, int32_t** hwe_hh_cases_ptr, int32_t** hwe_ll_allfs_ptr, int32_t** hwe_lh_allfs_ptr, int32_t** hwe_hh_allfs_ptr, int32_t** hwe_hapl_allfs_ptr, int32_t** hwe_haph_allfs_ptr, uint32_t* indiv_male_ct_ptr, uint32_t* indiv_f_ct_ptr, uint32_t* indiv_f_male_ct_ptr, uint32_t wt_needed, unsigned char** marker_weights_base_ptr, double** marker_weights_ptr, double exponent, Chrom_info* chrom_info_ptr, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t is_split_chrom, uint32_t* hh_exists_ptr) {
  FILE* hhfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_indiv_ctv2 = 2 * unfiltered_indiv_ctl;
  int32_t retval = 0;
  uint32_t pct = 1;
  uint32_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  double indiv_ct_recip = 1.0 / ((double)((int32_t)indiv_ct));
  // sum of nonmissing rates over all markers
  // rate is in [0, 1] for each marker, so sum is in [0, marker_ct].
  double nonmissing_rate_tot = 0.0;
  // track this to defend against Y chromosome/0 males pathological case
  uintptr_t nonmissing_rate_tot_max = marker_ct;
  uint32_t indiv_f_ct = indiv_ct;
  uintptr_t indiv_f_ctrl_ct = indiv_ct;
  uintptr_t indiv_f_case_ct = indiv_ct;
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t ll_hwe = 0;
  uint32_t lh_hwe = 0;
  uint32_t hh_hwe = 0;
  uint32_t ll_case_hwe = 0;
  uint32_t lh_case_hwe = 0;
  uint32_t hh_case_hwe = 0;
  uint32_t cur_chrom_idx = 0;
  uint32_t nonmissing_nonmale_y = 0;
  int32_t ii = chrom_info_ptr->chrom_file_order[0];
  uint32_t is_haploid = is_set(chrom_info_ptr->haploid_mask, ii);
  uint32_t next_chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[1];
  uint32_t is_x = (ii == chrom_info_ptr->x_code);
  uint32_t is_y = (ii == chrom_info_ptr->y_code);
  uint32_t ll_ct = 0;
  uint32_t lh_ct = 0;
  uint32_t hh_ct = 0;
  uint32_t ll_ctf = 0;
  uint32_t lh_ctf = 0;
  uint32_t hh_ctf = 0;
  uint32_t ukk = 0;
  int32_t* hwe_lls = NULL;
  int32_t* hwe_lhs = NULL;
  int32_t* hwe_hhs = NULL;
  int32_t* hwe_ll_cases = NULL;
  int32_t* hwe_lh_cases = NULL;
  int32_t* hwe_hh_cases = NULL;
  int32_t* hwe_ll_allfs = NULL;
  int32_t* hwe_lh_allfs = NULL;
  int32_t* hwe_hh_allfs = NULL;
  uintptr_t* indiv_nonmale_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* founder_nonmale_include2 = NULL;
  uintptr_t* founder_ctrl_nonmale_include2 = NULL;
  uintptr_t* founder_male_include2 = NULL;
  uintptr_t* founder_case_include2 = NULL;
  uintptr_t* founder_case_nonmale_include2 = NULL;
  double* marker_weights = NULL;
  uint32_t indiv_nonmale_ct = 0;
  uint32_t indiv_f_nonmale_ct = 0;
  uint32_t indiv_f_ctl_nonmale_ct = 0;
  uint32_t indiv_f_case_nonmale_ct = 0;
  uint64_t hethap_ct = 0;
  double male_ct_recip = 0;
  uint32_t indiv_male_ct;
  uint32_t indiv_f_male_ct;
  uint32_t hethap_incr;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uintptr_t* loadbuf;
  uintptr_t* indiv_include2;
  uintptr_t* founder_include2;
  uintptr_t* founder_ctrl_include2;
  uint32_t nonmales_needed;
  uint32_t males_needed;
  uintptr_t* tmp_indiv_excl_mask;
  uintptr_t* tmp_indiv_excl_mask2;
  uintptr_t loop_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t* marker_allele_cts;
  double maf;
  uii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(marker_reverse_ptr, uii * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  fill_ulong_zero(*marker_reverse_ptr, uii);
  if (wt_needed) {
    if (wkspace_alloc_uc_checked(marker_weights_base_ptr, CACHELINE) ||
        wkspace_alloc_d_checked(marker_weights_ptr, unfiltered_marker_ct * sizeof(double))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    marker_weights = *marker_weights_ptr;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      marker_weights[marker_uidx] = -1.0;
    }
  }

  if (!hwe_needed) {
    *hwe_lls_ptr = (int32_t*)wkspace_base;
  } else {
    if (wkspace_alloc_i_checked(&hwe_lls, unfiltered_marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_i_checked(&hwe_lhs, unfiltered_marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_i_checked(&hwe_hhs, unfiltered_marker_ct * sizeof(int32_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    *hwe_lls_ptr = hwe_lls;
    *hwe_lhs_ptr = hwe_lhs;
    *hwe_hhs_ptr = hwe_hhs;
    if (hardy_needed) {
      if (wkspace_alloc_i_checked(&hwe_ll_cases, unfiltered_marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&hwe_lh_cases, unfiltered_marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&hwe_hh_cases, unfiltered_marker_ct * sizeof(int32_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
    }
    *hwe_ll_cases_ptr = hwe_ll_cases;
    *hwe_lh_cases_ptr = hwe_lh_cases;
    *hwe_hh_cases_ptr = hwe_hh_cases;
  }
  if (wkspace_alloc_i_checked(&hwe_ll_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_lh_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_hh_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_hapl_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_i_checked(&hwe_haph_allfs, unfiltered_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&marker_allele_cts, 2 * unfiltered_marker_ct * sizeof(int32_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  *hwe_ll_allfs_ptr = hwe_ll_allfs;
  *hwe_lh_allfs_ptr = hwe_lh_allfs;
  *hwe_hh_allfs_ptr = hwe_hh_allfs;
  *hwe_hapl_allfs_ptr = hwe_hapl_allfs;
  *hwe_haph_allfs_ptr = hwe_haph_allfs;

  if ((!pheno_c) || is_split_chrom) {
    hwe_all = 1;
  }

  fill_int_zero(hwe_ll_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_lh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hh_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_hapl_allfs, unfiltered_marker_ct);
  fill_int_zero(hwe_haph_allfs, unfiltered_marker_ct);
  fill_uint_zero(marker_allele_cts, 2 * unfiltered_marker_ct);
  *marker_allele_cts_ptr = marker_allele_cts;
  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ctv2 - 2] = 0;
  loadbuf[unfiltered_indiv_ctv2 - 1] = 0;
  exclude_to_vec_include(unfiltered_indiv_ct, indiv_include2, indiv_exclude);
  ii = chrom_info_ptr->x_code;
  nonmales_needed = (!is_split_chrom) && (ii != -1) && is_set(chrom_info_ptr->chrom_mask, ii);
  ii = chrom_info_ptr->y_code;
  males_needed = nonmales_needed || ((!is_split_chrom) && (ii != -1) && is_set(chrom_info_ptr->chrom_mask, ii));
  if (wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  memcpy(indiv_male_include2, indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_nm);
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_male);
  indiv_male_ct = popcount_longs(indiv_male_include2, 0, unfiltered_indiv_ctv2);
  if (indiv_male_ct) {
    male_ct_recip = 1.0 / ((double)((int32_t)indiv_male_ct));
  }
  *indiv_male_ct_ptr = indiv_male_ct;
  indiv_f_male_ct = indiv_male_ct;
  if (males_needed) {
    founder_male_include2 = indiv_male_include2;
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(indiv_nonmale_include2, indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      vec_include_mask_out_intersect(unfiltered_indiv_ct, indiv_nonmale_include2, sex_nm, sex_male);
      indiv_nonmale_ct = popcount_longs(indiv_nonmale_include2, 0, unfiltered_indiv_ctv2);
      founder_nonmale_include2 = indiv_nonmale_include2;
      founder_ctrl_nonmale_include2 = indiv_nonmale_include2;
    }
  }
  founder_include2 = indiv_include2;
  if (wkspace_alloc_ul_checked(&tmp_indiv_excl_mask, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto calc_freqs_and_hwe_ret_NOMEM;
  }
  if (!nonfounders) {
    if (wkspace_alloc_ul_checked(&founder_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    for (uii = 0; uii < unfiltered_indiv_ctl; uii++) {
      tmp_indiv_excl_mask[uii] = indiv_exclude[uii] | (~founder_info[uii]);
    }
    zero_trailing_bits(tmp_indiv_excl_mask, unfiltered_indiv_ct);
    exclude_to_vec_include(unfiltered_indiv_ct, founder_include2, tmp_indiv_excl_mask);
    if (males_needed) {
      if (wkspace_alloc_ul_checked(&founder_male_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_male_include2, indiv_male_include2, unfiltered_indiv_ctl * 2 * sizeof(intptr_t));
      vec_include_mask_in(unfiltered_indiv_ct, founder_male_include2, founder_info);
      indiv_f_male_ct = popcount_longs(founder_male_include2, 0, unfiltered_indiv_ctv2);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
	vec_include_mask_in(unfiltered_indiv_ct, founder_nonmale_include2, founder_info);
	indiv_f_nonmale_ct = popcount_longs(founder_nonmale_include2, 0, unfiltered_indiv_ctv2);
      }
    }
    founder_ctrl_include2 = founder_include2;
    uii = popcount_longs_exclude(founder_info, indiv_exclude, unfiltered_indiv_ctl);
    indiv_f_ct = uii;
    indiv_f_ctrl_ct = indiv_f_ct;
  } else {
    for (uii = 0; uii < unfiltered_indiv_ctl; uii++) {
      tmp_indiv_excl_mask[uii] = indiv_exclude[uii];
    }
    founder_ctrl_include2 = founder_include2;
  }

  if (!hwe_all) {
    if (wkspace_alloc_ul_checked(&founder_ctrl_include2, unfiltered_indiv_ctv2 *  sizeof(intptr_t)) ||
	wkspace_alloc_ul_checked(&tmp_indiv_excl_mask2, unfiltered_indiv_ctl * sizeof(intptr_t))) {
      goto calc_freqs_and_hwe_ret_NOMEM;
    }
    indiv_uidx = 0;
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctl; indiv_uidx++) {
      tmp_indiv_excl_mask2[indiv_uidx] = tmp_indiv_excl_mask[indiv_uidx] | (~(pheno_nm[indiv_uidx])) | pheno_c[indiv_uidx];
    }
    zero_trailing_bits(tmp_indiv_excl_mask2, unfiltered_indiv_ct);
    // tmp_indiv_excl_mask is now set for each indiv who is excluded, or a
    // nonfounder, or is noncontrol.
    indiv_f_ctrl_ct = unfiltered_indiv_ct - popcount_longs(tmp_indiv_excl_mask2, 0, unfiltered_indiv_ctl);
    exclude_to_vec_include(unfiltered_indiv_ct, founder_ctrl_include2, tmp_indiv_excl_mask2);
    if (nonmales_needed) {
      if (wkspace_alloc_ul_checked(&founder_ctrl_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      memcpy(founder_ctrl_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      vec_include_mask_out(unfiltered_indiv_ct, founder_ctrl_nonmale_include2, tmp_indiv_excl_mask2);
      indiv_f_ctl_nonmale_ct = popcount_longs(founder_ctrl_nonmale_include2, 0, unfiltered_indiv_ctv2);
    }
    if (hardy_needed) {
      if (wkspace_alloc_ul_checked(&founder_case_include2, unfiltered_indiv_ctv2 *  sizeof(intptr_t))) {
	goto calc_freqs_and_hwe_ret_NOMEM;
      }
      indiv_uidx = 0;
      for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctl; indiv_uidx++) {
	tmp_indiv_excl_mask[indiv_uidx] |= (~(pheno_nm[indiv_uidx])) | (~pheno_c[indiv_uidx]);
      }
      zero_trailing_bits(tmp_indiv_excl_mask, unfiltered_indiv_ct);
      indiv_f_case_ct = unfiltered_indiv_ct - popcount_longs(tmp_indiv_excl_mask, 0, unfiltered_indiv_ctl);
      exclude_to_vec_include(unfiltered_indiv_ct, founder_case_include2, tmp_indiv_excl_mask);
      if (nonmales_needed) {
	if (wkspace_alloc_ul_checked(&founder_case_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
	  goto calc_freqs_and_hwe_ret_NOMEM;
	}
	memcpy(founder_case_nonmale_include2, indiv_nonmale_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
	vec_include_mask_out(unfiltered_indiv_ct, founder_case_nonmale_include2, tmp_indiv_excl_mask);
	indiv_f_ctl_nonmale_ct = popcount_longs(founder_case_nonmale_include2, 0, unfiltered_indiv_ctv2);
      }
    }
  }

  *indiv_f_ct_ptr = indiv_f_ct;
  *indiv_f_male_ct_ptr = indiv_f_male_ct;
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_freqs_and_hwe_ret_READ_FAIL;
  }
  marker_uidx = 0;
  marker_idx = 0;
  logprint("Calculating allele frequencies...");
  fputs(" 0%", stdout);
  fflush(stdout);
  if (is_split_chrom) {
    // only set is_haploid if all chromosomes are haploid
    is_haploid = (chrom_info_ptr->chrom_mask[0]) & 1;
    is_x = 0;
    is_y = 0;
    next_chrom_start = unfiltered_marker_ct;
  }
  for (; pct <= 100; pct++) {
    loop_end = ((uint64_t)pct * marker_ct) / 100LU;
    for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_freqs_and_hwe_ret_READ_FAIL;
	}
      }
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto calc_freqs_and_hwe_ret_READ_FAIL;
      }
      if (marker_uidx >= next_chrom_start) {
	do {
	  next_chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[(++cur_chrom_idx) + 1];
	} while (marker_uidx >= next_chrom_start);
	ii = chrom_info_ptr->chrom_file_order[cur_chrom_idx];
	is_haploid = is_set(chrom_info_ptr->haploid_mask, ii);
	is_x = (ii == chrom_info_ptr->x_code);
	is_y = (ii == chrom_info_ptr->y_code);
      }
      if (!is_haploid) {
	single_marker_freqs_and_hwe(unfiltered_indiv_ctv2, loadbuf, indiv_include2, founder_include2, founder_ctrl_include2, founder_case_include2, indiv_ct, &ll_ct, &lh_ct, &hh_ct, indiv_f_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, indiv_f_ctrl_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, indiv_f_case_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	hwe_ll_allfs[marker_uidx] = ll_ctf;
	hwe_lh_allfs[marker_uidx] = lh_ctf;
	hwe_hh_allfs[marker_uidx] = hh_ctf;
	marker_allele_cts[2 * marker_uidx] += 2 * ll_ct + lh_ct;
	marker_allele_cts[2 * marker_uidx + 1] += 2 * hh_ct + lh_ct;
	uii = 2 * (ll_ctf + lh_ctf + hh_ctf + maf_succ);
	nonmissing_rate_tot += (ll_ct + lh_ct + hh_ct) * indiv_ct_recip;
	if (!uii) {
	  // avoid 0/0 division
	  set_allele_freqs[marker_uidx] = 0.5;
	} else {
	  set_allele_freqs[marker_uidx] = ((double)(2 * hh_ctf + lh_ctf + maf_succ)) / ((double)uii);
	}
	if (hwe_needed) {
	  hwe_lls[marker_uidx] = ll_hwe;
	  hwe_lhs[marker_uidx] = lh_hwe;
	  hwe_hhs[marker_uidx] = hh_hwe;
	  if (hardy_needed) {
	    hwe_ll_cases[marker_uidx] = ll_case_hwe;
	    hwe_lh_cases[marker_uidx] = lh_case_hwe;
	    hwe_hh_cases[marker_uidx] = hh_case_hwe;
	  }
	}
      } else {
	uii = 0;
	ujj = 0;
	if (is_x || is_y) {
	  if (is_x) {
	    single_marker_freqs_and_hwe(unfiltered_indiv_ctv2, loadbuf, indiv_nonmale_include2, founder_nonmale_include2, founder_ctrl_nonmale_include2, founder_case_nonmale_include2, indiv_nonmale_ct, &ll_ct, &lh_ct, &hh_ct, indiv_f_nonmale_ct, &ll_ctf, &lh_ctf, &hh_ctf, hwe_needed, indiv_f_ctl_nonmale_ct, &ll_hwe, &lh_hwe, &hh_hwe, hardy_needed, indiv_f_case_nonmale_ct, &ll_case_hwe, &lh_case_hwe, &hh_case_hwe);
	    hwe_ll_allfs[marker_uidx] = ll_ctf;
	    hwe_lh_allfs[marker_uidx] = lh_ctf;
	    hwe_hh_allfs[marker_uidx] = hh_ctf;
	    marker_allele_cts[2 * marker_uidx] += 2 * ll_ct + lh_ct;
	    marker_allele_cts[2 * marker_uidx + 1] += 2 * hh_ct + lh_ct;
	    uii = 2 * (ll_ctf + lh_ctf + hh_ctf);
	    ujj = 2 * hh_ctf + lh_ctf;
	    ukk = ll_ct + lh_ct + hh_ct;
	    if (hwe_needed) {
	      hwe_lls[marker_uidx] = ll_hwe;
	      hwe_lhs[marker_uidx] = lh_hwe;
	      hwe_hhs[marker_uidx] = hh_hwe;
	      if (hardy_needed) {
		hwe_ll_cases[marker_uidx] = ll_case_hwe;
		hwe_lh_cases[marker_uidx] = lh_case_hwe;
		hwe_hh_cases[marker_uidx] = hh_case_hwe;
	      }
	    }
	  } else if (!nonmissing_nonmale_y) {
	    nonmissing_nonmale_y = nonmissing_present_diff(unfiltered_indiv_ctv2, loadbuf, indiv_include2, indiv_male_include2);
	  }
	  haploid_single_marker_freqs(unfiltered_indiv_ct, unfiltered_indiv_ctv2, loadbuf, indiv_male_include2, founder_male_include2, indiv_male_ct, &ll_ct, &hh_ct, indiv_f_male_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	  if (is_x) {
	    nonmissing_rate_tot += ((int32_t)(ll_ct + hh_ct + ukk)) * indiv_ct_recip;
	  } else if (indiv_male_ct) {
	    nonmissing_rate_tot += ((int32_t)(ll_ct + hh_ct)) * male_ct_recip;
	  } else {
	    nonmissing_rate_tot_max -= 1;
	  }
	} else {
	  haploid_single_marker_freqs(unfiltered_indiv_ct, unfiltered_indiv_ctv2, loadbuf, indiv_include2, founder_include2, indiv_ct, &ll_ct, &hh_ct, indiv_f_ct, &ll_ctf, &hh_ctf, &hethap_incr);
	  nonmissing_rate_tot += ((int32_t)(ll_ct + hh_ct)) * indiv_ct_recip;
	}
	if (hethap_incr) {
	  if (!hhfile) {
	    memcpy(outname_end, ".hh", 4);
	    if (fopen_checked(&hhfile, outname, "w")) {
	      goto calc_freqs_and_hwe_ret_OPEN_FAIL;
	    }
	  }
	  if (is_x) {
	    *hh_exists_ptr |= XMHH_EXISTS;
	  } else if (is_y) {
	    *hh_exists_ptr |= Y_FIX_NEEDED;
	  } else {
	    *hh_exists_ptr |= NXMHH_EXISTS;
	  }
	  if (is_x || is_y) {
	    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctv2; indiv_uidx++) {
	      ulii = loadbuf[indiv_uidx];
	      ulii = (ulii >> 1) & (~ulii) & indiv_male_include2[indiv_uidx];
	      while (ulii) {
		ukk = indiv_uidx * BITCT2 + CTZLU(ulii) / 2;
		fputs(&(person_ids[ukk * max_person_id_len]), hhfile);
		putc('\t', hhfile);
		fputs(&(marker_ids[marker_uidx * max_marker_id_len]), hhfile);
		putc('\n', hhfile);
		ulii &= ulii - ONELU;
	      }
	    }
	  } else {
	    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ctv2; indiv_uidx++) {
	      ulii = loadbuf[indiv_uidx];
	      ulii = (ulii >> 1) & (~ulii) & indiv_include2[indiv_uidx];
	      while (ulii) {
		ukk = indiv_uidx * BITCT2 + CTZLU(ulii) / 2;
		fputs(&(person_ids[ukk * max_person_id_len]), hhfile);
		putc('\t', hhfile);
		fputs(&(marker_ids[marker_uidx * max_marker_id_len]), hhfile);
		putc('\n', hhfile);
		ulii &= ulii - ONELU;
	      }
	    }
	  }
	  if (ferror(hhfile)) {
	    goto calc_freqs_and_hwe_ret_WRITE_FAIL;
	  }
	  hethap_ct += hethap_incr;
	}
	hwe_hapl_allfs[marker_uidx] = ll_ctf;
	hwe_haph_allfs[marker_uidx] = hh_ctf;
	marker_allele_cts[2 * marker_uidx] += ll_ct;
	marker_allele_cts[2 * marker_uidx + 1] += hh_ct;
	uii += ll_ctf + hh_ctf + 2 * maf_succ;
	ujj += hh_ctf + maf_succ;
	if (!uii) {
	  maf = 0.5;
	} else {
	  maf = ((double)ujj) / ((double)uii);
	}
	set_allele_freqs[marker_uidx] = maf;
	if (wt_needed) {
	  marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, maf);
	}
      }
    }
    if (pct < 100) {
      if (pct > 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");
  if (hethap_ct) {
    *outname_end = '\0';
    sprintf(logbuf, "Warning: %" PRIu64 " het. haploid genotype%s present (see %s.hh).\n", hethap_ct, (hethap_ct == 1LLU)? "" : "s", outname);
    logprintb();
  }
  if (nonmissing_nonmale_y) {
    logprint("Warning: Nonmissing nonmale Y chromosome genotype(s) present.\n");
    *hh_exists_ptr |= Y_FIX_NEEDED;
  }
  if (nonmissing_rate_tot <= 0.9999995 * ((double)((intptr_t)nonmissing_rate_tot_max))) {
    sprintf(logbuf, "Total genotyping rate %sis %g.\n", indiv_exclude_ct? "in remaining individuals " : "", nonmissing_rate_tot / ((double)((intptr_t)nonmissing_rate_tot_max)));
    logprintb();
  }
  while (0) {
  calc_freqs_and_hwe_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_freqs_and_hwe_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_freqs_and_hwe_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_freqs_and_hwe_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(hhfile);
  return retval;
}

// aptr1 = minor, aptr2 = major
int32_t load_one_freq(uint32_t alen1, const char* aptr1, uint32_t alen2, const char* aptr2, double maf, double* set_allele_freq_ptr, char* mastr1, char* mastr2, char missing_geno) {
  uint32_t malen1 = strlen(mastr1);
  uint32_t malen2 = strlen(mastr2);
  uint32_t uii;
  const char* aptr;
  if (maf > 0.5) {
    aptr = aptr2;
    uii = alen2;
    aptr2 = aptr1;
    alen2 = alen1;
    aptr1 = aptr;
    alen1 = uii;
    maf = 1.0 - maf;
  }
  if ((malen1 == alen1) && (!memcmp(mastr1, aptr1, alen1))) {
    if ((malen2 == alen2) && (!memcmp(mastr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0 - maf;
    } else {
      return -1;
    }
  } else if ((malen2 == alen1) && (!memcmp(mastr2, aptr1, alen1))) {
    if ((malen1 == alen2) && (!memcmp(mastr1, aptr2, alen2))) {
      *set_allele_freq_ptr = maf;
    } else {
      return -1;
    }
  } else if ((*aptr1 == missing_geno) && (alen1 == 1) && (maf == 0.0)) {
    if ((malen1 == alen2) && (!memcmp(mastr1, aptr2, alen2))) {
      *set_allele_freq_ptr = 0.0;
    } else if ((malen2 == alen2) && (!memcmp(mastr2, aptr2, alen2))) {
      *set_allele_freq_ptr = 1.0;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

int32_t read_external_freqs(char* freqname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs, uint32_t* marker_allele_cts, double* set_allele_freqs, uint32_t maf_succ, double exponent, uint32_t wt_needed, double* marker_weights) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* freqfile = NULL;
  uint32_t freq_counts = 0;
  uint32_t alen1 = 0;
  uint32_t alen2 = 0;
  char* aptr1 = NULL;
  char* aptr2 = NULL;
  int32_t retval = 0;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;
  char* loadbuf;
  char* sorted_ids;
  uint32_t* id_map;
  uintptr_t loadbuf_size;
  uint32_t chrom_idx;
  uint32_t marker_uidx;
  uint32_t uii;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  double maf;
  int32_t c_hom_a1;
  int32_t c_het;
  int32_t c_hom_a2;
  int32_t c_hap_a1;
  int32_t c_hap_a2;
  int32_t ii;
  if (fopen_checked(&freqfile, freqname, "r")) {
    goto read_external_freqs_ret_OPEN_FAIL;
  }
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto read_external_freqs_ret_1;
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto read_external_freqs_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  if (fgets(loadbuf, loadbuf_size, freqfile) == NULL) {
    logprint("Error: Empty --read-freq file.\n");
    goto read_external_freqs_ret_INVALID_FORMAT_2;
  }
  if (!memcmp(loadbuf, " CHR  ", 6)) {
    uii = strlen(loadbuf);
    if (loadbuf[uii - 2] == '0') { // --counts makes G0 the last column header
      freq_counts = 1;
    } else if (loadbuf[uii - 2] != 'S') { // NCHROBS
      goto read_external_freqs_ret_INVALID_FORMAT;
    }
    while (fgets(loadbuf, loadbuf_size, freqfile) != NULL) {
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
      ii = marker_code(chrom_info_ptr, bufptr);
      if (ii == -1) {
	goto read_external_freqs_ret_INVALID_FORMAT;
      }
      chrom_idx = ii;
      bufptr = next_item(bufptr); // now at beginning of marker name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
        if ((chrom_idx == get_marker_chrom(chrom_info_ptr, marker_uidx)) || (!chrom_idx) || (!get_marker_chrom(chrom_info_ptr, marker_uidx))) {
	  alen1 = strlen_se(bufptr2);
	  aptr1 = bufptr2;
	  bufptr2 = next_item(bufptr2);
	  if (no_more_items_kns(bufptr2)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  alen2 = strlen_se(bufptr2);
	  aptr2 = bufptr2;
	  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr2);
	  if (no_more_items_kns(bufptr)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  if (freq_counts) {
	    if (no_more_items_kns(next_item(bufptr))) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	    c_hom_a1 = atoi(bufptr);
	    c_hom_a2 = atoi(next_item(bufptr));
	    maf = ((double)c_hom_a1 + maf_succ) / ((double)(c_hom_a1 + c_hom_a2 + 2 * maf_succ));
	  } else {
	    if (scan_double(bufptr, &maf)) {
	      goto read_external_freqs_ret_INVALID_FORMAT;
	    }
	  }
	  if (load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	  }
        }
      }
    }
    if (freq_counts) {
      logprint(".frq.count file loaded.\n");
    } else {
      logprint(".frq file loaded.\n");
    }
  } else if (!memcmp(loadbuf, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)", 71)) {
    // changed from strcmp to avoid eoln problems
    // known --freqx format, v0.15.3 or later
    while (fgets(loadbuf, loadbuf_size, freqfile) != NULL) {
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      ii = marker_code(chrom_info_ptr, loadbuf);
      if (ii == -1) {
	goto read_external_freqs_ret_INVALID_FORMAT;
      }
      chrom_idx = ii;
      bufptr = next_item(loadbuf); // now at beginning of marker name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
        if ((chrom_idx == get_marker_chrom(chrom_info_ptr, marker_uidx)) || (!chrom_idx) || (!get_marker_chrom(chrom_info_ptr, marker_uidx))) {
	  alen1 = strlen_se(bufptr2);
	  aptr1 = bufptr2;
	  bufptr2 = next_item(bufptr2);
	  if (no_more_items_kns(bufptr2)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  alen2 = strlen_se(bufptr2);
	  aptr2 = bufptr2;
	  if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  bufptr = next_item(bufptr2);
	  bufptr2 = next_item(bufptr);
	  bufptr3 = next_item(bufptr2);
	  bufptr4 = next_item(bufptr3);
	  bufptr5 = next_item(bufptr4);
	  if (no_more_items_kns(bufptr5)) {
	    goto read_external_freqs_ret_INVALID_FORMAT;
	  }
	  c_hom_a1 = atoi(bufptr);
	  c_het = atoi(bufptr2);
	  c_hom_a2 = atoi(bufptr3);
	  c_hap_a1 = atoi(bufptr4);
	  c_hap_a2 = atoi(bufptr5);
	  maf = ((double)(c_hom_a1 * 2 + c_het + c_hap_a1 + maf_succ)) / ((double)(2 * (c_hom_a1 + c_het + c_hom_a2 + maf_succ) + c_hap_a1 + c_hap_a2));
	  if (load_one_freq(alen1, aptr1, alen2, aptr2, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	    goto read_external_freqs_ret_ALLELE_MISMATCH;
	  }
	  if (wt_needed) {
	    if (c_hap_a1 || c_hap_a2) {
	      marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	    } else {
	      marker_weights[marker_uidx] = calc_wt_mean(exponent, c_het, c_hom_a1, c_hom_a2);
	    }
	  }
        }
      }
    }
    logprint(".frqx file loaded.\n");
  } else {
    // Also support GCTA-style frequency files:
    // [marker ID]\t[reference allele]\t[frequency of reference allele]\n
    do {
      if (!loadbuf[loadbuf_size - 1]) {
	goto read_external_freqs_ret_TOO_LONG_LINE;
      }
      bufptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr = next_item(bufptr);
      if (!bufptr) {
        goto read_external_freqs_ret_INVALID_FORMAT;
      }
      read_next_terminate(loadbuf, loadbuf); // destructive read
      ii = bsearch_str(loadbuf, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        marker_uidx = id_map[(uint32_t)ii];
	alen1 = strlen_se(bufptr);
	aptr1 = bufptr;
        bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr)) {
          goto read_external_freqs_ret_INVALID_FORMAT;
	}
	if (scan_double(bufptr, &maf)) {
          goto read_external_freqs_ret_INVALID_FORMAT;
        }
	if (load_one_freq(1, missing_geno_ptr, alen1, aptr1, maf, &(set_allele_freqs[marker_uidx]), marker_allele_ptrs[marker_uidx * 2], marker_allele_ptrs[marker_uidx * 2 + 1], missing_geno)) {
	  goto read_external_freqs_ret_ALLELE_MISMATCH;
	}
	if (wt_needed) {
	  marker_weights[marker_uidx] = calc_wt_mean_maf(exponent, set_allele_freqs[marker_uidx]);
	}
      } else {
	// if there aren't exactly 3 columns, this isn't a GCTA .freq file
	bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr) || (!no_more_items_kns(next_item(bufptr)))) {
	  goto read_external_freqs_ret_INVALID_FORMAT;
	}
      }
    } while (fgets(loadbuf, loadbuf_size, freqfile) != NULL);
    logprint("GCTA-formatted .freq file loaded.\n");
  }
  while (0) {
  read_external_freqs_ret_TOO_LONG_LINE:
    if (loadbuf_size == MAXLINEBUFLEN) {
      logprint("Error: Pathologically long line in --freq[x] file.\n");
      retval = RET_INVALID_FORMAT;
      break;
    }
  read_external_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  read_external_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_external_freqs_ret_INVALID_FORMAT:
    logprint(errstr_freq_format);
  read_external_freqs_ret_INVALID_FORMAT_2:
    retval = RET_INVALID_FORMAT;
    break;
  read_external_freqs_ret_ALLELE_MISMATCH:
    sprintf(logbuf, "Error: Mismatch between .bim/.ped and --read-freq alleles at %s.\n", next_item(skip_initial_spaces(loadbuf)));
    logprintb();
    retval = RET_ALLELE_MISMATCH;
    break;
  }
 read_external_freqs_ret_1:
  fclose_cond(freqfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t write_freqs(char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t indiv_f_ct, uint32_t indiv_f_male_ct, uint64_t misc_flags, uintptr_t* marker_reverse) {
  FILE* outfile = NULL;
  uint32_t reverse = 0;
  uint32_t freq_counts = (misc_flags / MISC_FREQ_COUNTS) & 1;
  uint32_t freqx = (misc_flags / MISC_FREQX) & 1;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  char* tbuf2 = &(tbuf[MAXLINELEN]);
  int32_t retval = 0;
  char* minor_ptr;
  char* major_ptr;
  char* bufptr;
  uint32_t chrom_end;
  uint32_t marker_uidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t missing_ct;
  int32_t chrom_idx;
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_freqs_ret_OPEN_FAIL;
  }
  if (freqx) {
    if (fputs_checked("CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)\n", outfile)) {
      goto write_freqs_ret_WRITE_FAIL;
    }
  } else if (plink_maxsnp < 5) {
    if (freq_counts) {
      if (fputs_checked(" CHR  SNP   A1   A2     C1     C2     G0\n", outfile)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
      strcpy(tbuf, " %4s %4s %4s %6u %6u %6u\n");
    } else {
      if (fputs_checked(" CHR  SNP   A1   A2          MAF  NCHROBS\n", outfile)) {
        goto write_freqs_ret_WRITE_FAIL;
      }
      strcpy(tbuf, " %4s %4s %4s %12.4g %8d\n");
    }
  } else if (freq_counts) {
    sprintf(tbuf, " CHR %%%us   A1   A2     C1     C2     G0\n", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
    sprintf(tbuf, " %%%us %%4s %%4s %%6u %%6u %%6u\n", plink_maxsnp);
  } else {
    sprintf(tbuf, " CHR %%%us   A1   A2          MAF  NCHROBS\n", plink_maxsnp);
    fprintf(outfile, tbuf, "SNP");
    sprintf(tbuf, " %%%us %%4s %%4s %%12.4g %%8d\n", plink_maxsnp);
  }
  if (ferror(outfile)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->x_code);
    is_y = (chrom_idx == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end);
    while (marker_uidx < chrom_end) {
      reverse = IS_SET(marker_reverse, marker_uidx);
      major_ptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      minor_ptr = marker_allele_ptrs[marker_uidx * 2];
      if (freq_counts || freqx) {
	if (is_x) {
	  missing_ct = indiv_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx] + hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	} else if (is_haploid) {
	  if (is_y) {
	    missing_ct = indiv_f_male_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  } else {
	    missing_ct = indiv_f_ct - (hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
	  }
	} else {
	  missing_ct = indiv_f_ct - (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]);
	}
	if (freqx) {
	  bufptr = chrom_name_write(tbuf2, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms);
	  fwrite(tbuf2, 1, bufptr - tbuf2, outfile);
	  fprintf(outfile, "\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\n", &(marker_ids[marker_uidx * max_marker_id_len]), minor_ptr, major_ptr, reverse? hh_cts[marker_uidx] : ll_cts[marker_uidx], lh_cts[marker_uidx], reverse? ll_cts[marker_uidx] : hh_cts[marker_uidx], reverse? haph_cts[marker_uidx] : hapl_cts[marker_uidx], reverse? hapl_cts[marker_uidx] : haph_cts[marker_uidx], missing_ct);
	} else {
	  bufptr = width_force(4, tbuf2, chrom_name_write(tbuf2, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
	  fwrite(tbuf2, 1, bufptr - tbuf2, outfile);
	  fprintf(outfile, tbuf, &(marker_ids[marker_uidx * max_marker_id_len]), minor_ptr, major_ptr, 2 * ll_cts[marker_uidx] + lh_cts[marker_uidx] + hapl_cts[marker_uidx], 2 * hh_cts[marker_uidx] + lh_cts[marker_uidx] + haph_cts[marker_uidx], missing_ct);
	}
      } else {
	bufptr = width_force(4, tbuf2, chrom_name_write(tbuf2, chrom_info_ptr, get_marker_chrom(chrom_info_ptr, marker_uidx), zero_extra_chroms));
	fwrite(tbuf2, 1, bufptr - tbuf2, outfile);
	fprintf(outfile, tbuf, &(marker_ids[marker_uidx * max_marker_id_len]), minor_ptr, major_ptr, (reverse? set_allele_freqs[marker_uidx] : (1.0 - set_allele_freqs[marker_uidx])) + SMALL_EPSILON, 2 * (ll_cts[marker_uidx] + lh_cts[marker_uidx] + hh_cts[marker_uidx]) + hapl_cts[marker_uidx] + haph_cts[marker_uidx]);
      }
      if (ferror(outfile)) {
	goto write_freqs_ret_WRITE_FAIL;
      }
      marker_uidx = next_unset(marker_exclude, marker_uidx + 1, chrom_end);
    }
  }
  if (fclose_null(&outfile)) {
    goto write_freqs_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Allele frequencies written to %s.\n", outname);
  logprintb();
  while (0) {
  write_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t write_stratified_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uint32_t indiv_f_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uint32_t indiv_f_male_ct, uintptr_t* marker_reverse, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uint32_t* cur_cluster_map = cluster_map;
  uint32_t* cur_cluster_starts = cluster_starts;
  uint32_t* cluster_map_nonmale = NULL;
  uint32_t* cluster_starts_nonmale = NULL;
  uint32_t* cluster_map_male = NULL;
  uint32_t* cluster_starts_male = NULL;
  int32_t chrom_code_end = chrom_info_ptr->max_code + 1 + chrom_info_ptr->name_ct;
  uint32_t cslen = 10;
  int32_t retval = 0;
  uint32_t cur_cts[4];
  uintptr_t* readbuf;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t* uiptr3;
  char* csptr;
  char* col_2_start;
  char* wptr_start;
  char* wptr;
  char* sptr;
  uintptr_t chrom_end;
  uintptr_t marker_uidx;
  int32_t chrom_idx;
  uintptr_t clidx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t clmpos;
  uint32_t a1_obs;
  uint32_t tot_obs;
  uint32_t uii;
  if (wkspace_alloc_ul_checked(&readbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  if ((indiv_ct > indiv_f_ct) && (!nonfounders)) {
    if (wkspace_alloc_ui_checked(&cur_cluster_starts, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cur_cluster_map, indiv_f_ct * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cur_cluster_starts[0] = 0;
    uiptr = cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cluster_map[cluster_starts[clidx + 1]]);
      do {
	uii = *uiptr;
	if (IS_SET(founder_info, uii)) {
          cur_cluster_map[clmpos++] = uii;
	}
      } while ((++uiptr) < uiptr2);
      cur_cluster_starts[clidx + 1] = clmpos;
    }
  }
  chrom_idx = chrom_info_ptr->x_code;
  if ((chrom_idx != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
    if (wkspace_alloc_ui_checked(&cluster_starts_nonmale, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_map_nonmale, (indiv_f_ct - indiv_f_male_ct) * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cluster_starts_nonmale[0] = 0;
    uiptr = cur_cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
      while (uiptr < uiptr2) {
	uii = *uiptr++;
	if (!IS_SET(sex_male, uii)) {
          cluster_map_nonmale[clmpos++] = uii;
	}
      }
      cluster_starts_nonmale[clidx + 1] = clmpos;
    }
  }
  chrom_idx = chrom_info_ptr->y_code;
  if (cluster_map_nonmale || ((chrom_idx != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_idx))) {
    if (wkspace_alloc_ui_checked(&cluster_starts_male, (cluster_ct + 1) * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_map_male, indiv_f_male_ct * sizeof(int32_t))) {
      goto write_stratified_freqs_ret_NOMEM;
    }
    clmpos = 0;
    cluster_starts_male[0] = 0;
    uiptr = cur_cluster_map;
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
      while (uiptr < uiptr2) {
	uii = *uiptr++;
	if (IS_SET(sex_male, uii)) {
          cluster_map_male[clmpos++] = uii;
	}
      }
      cluster_starts_male[clidx + 1] = clmpos;
    }
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_stratified_freqs_ret_OPEN_FAIL;
  }
  sprintf(tbuf, " CHR %%%ds     CLST   A1   A2      MAF    MAC  NCHROBS\n", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  if (wkspace_alloc_c_checked(&csptr, 2 * max_marker_allele_len + 16)) {
    goto write_stratified_freqs_ret_NOMEM;
  }
  memset(csptr, 32, 10);
  for (chrom_idx = 0; chrom_idx < chrom_code_end; chrom_idx++) {
    if (!chrom_exists(chrom_info_ptr, chrom_idx)) {
      continue;
    }
    is_x = (chrom_idx == chrom_info_ptr->x_code);
    is_y = (chrom_idx == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
    marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_start[chrom_idx], chrom_end);
    if (marker_uidx >= chrom_end) {
      continue;
    }
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto write_stratified_freqs_ret_READ_FAIL;
    }
    col_2_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms));
    *col_2_start++ = ' ';
    do {
      sptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      uii = strlen(sptr);
      wptr_start = memseta(col_2_start, 32, plink_maxsnp - uii);
      wptr_start = memcpyax(wptr_start, sptr, uii, ' ');
      sptr = marker_allele_ptrs[marker_uidx * 2];
      wptr = fw_strcpy(4, sptr, &(csptr[1]));
      *wptr++ = ' ';
      sptr = marker_allele_ptrs[marker_uidx * 2 + 1];
      wptr = fw_strcpy(4, sptr, wptr);
      *wptr++ = ' ';
      cslen = (uintptr_t)(wptr - csptr);

      if (fread(readbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto write_stratified_freqs_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf((unsigned char*)readbuf, unfiltered_indiv_ct);
      }
      if (is_x) {
	uiptr = cluster_map_nonmale;
	uiptr2 = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr3 = &(cluster_map_nonmale[cluster_starts_nonmale[clidx + 1]]);
	  while (uiptr < uiptr3) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs = 2 * cur_cts[0] + cur_cts[2];
	  tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  fill_uint_zero(cur_cts, 4);
	  uiptr3 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr2 < uiptr3) {
	    uii = *uiptr2++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  a1_obs += cur_cts[0];
	  tot_obs += cur_cts[0] + cur_cts[3];
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else if (is_y) {
	uiptr = cluster_map_male;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr2 = &(cluster_map_male[cluster_starts_male[clidx + 1]]);
	  while (uiptr < uiptr2) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  if (is_haploid) {
	    a1_obs = cur_cts[0];
	    tot_obs = cur_cts[0] + cur_cts[3];
	  } else {
	    a1_obs = 2 * cur_cts[0] + cur_cts[2];
	    tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  }
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      } else {
        uiptr = cur_cluster_map;
	for (clidx = 0; clidx < cluster_ct; clidx++) {
	  wptr = fw_strcpy(8, &(cluster_ids[clidx * max_cluster_id_len]), wptr_start);
	  wptr = memcpyax(wptr, csptr, cslen, ' ');
	  fill_uint_zero(cur_cts, 4);
	  uiptr2 = &(cur_cluster_map[cur_cluster_starts[clidx + 1]]);
	  while (uiptr < uiptr2) {
	    uii = *uiptr++;
	    cur_cts[(readbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & (ONELU * 3)] += 1;
	  }
	  if (is_haploid) {
	    a1_obs = cur_cts[0];
	    tot_obs = cur_cts[0] + cur_cts[3];
	  } else {
	    a1_obs = 2 * cur_cts[0] + cur_cts[2];
	    tot_obs = 2 * (cur_cts[0] + cur_cts[2] + cur_cts[3]);
	  }
	  if (tot_obs) {
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)a1_obs)) / ((double)tot_obs), 8, ' ');
	    wptr = uint32_writew6x(wptr, a1_obs, ' ');
	    wptr = uint32_writew8(wptr, tot_obs);
	    wptr = memcpya(wptr, " \n", 2);
	  } else {
	    wptr = memcpya(wptr, "       0      0        0 \n", 26);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_stratified_freqs_ret_WRITE_FAIL;
	  }
	}
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto write_stratified_freqs_ret_READ_FAIL;
	}
      }
    } while (marker_uidx < chrom_end);
  }
  if (fclose_null(&outfile)) {
    goto write_stratified_freqs_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Cluster-stratified allele frequencies written to %s.\n", outname);
  logprintb();
  while (0) {
  write_stratified_freqs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_stratified_freqs_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_stratified_freqs_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_stratified_freqs_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t write_missingness_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uintptr_t* pheno_nm, uintptr_t* sex_male, uint32_t indiv_male_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ct2l = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_indiv_ctv2 = (unfiltered_indiv_ct2l + 1) & (~1);
  uintptr_t marker_ct_y = 0;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t* indiv_to_cluster = NULL;
  uint32_t* missing_ct_by_cluster = NULL;
  uint32_t* cluster_sizes = NULL;
  uint32_t* cluster_sizes_y = NULL;
  uint32_t indiv_uidx = 0;
  uint32_t indiv_idx = 0;
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* indiv_include2;
  uintptr_t* cur_nm;
  uintptr_t* lptr;
  uintptr_t* lptr2;
  uint32_t* missing_cts;
  uint32_t* cur_cluster_sizes;
  char* wptr;
  char* cptr;
  char* cptr2;
  uintptr_t marker_ct_nony;
  uintptr_t marker_uidx;
  uintptr_t ulii;
  uint32_t slen;
  uint32_t indiv_uidx_stop;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t cur_tot;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t clidx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  if (wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto write_missingness_reports_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ctv2 - 2] = 0;
  loadbuf[unfiltered_indiv_ctv2 - 1] = 0;
  exclude_to_vec_include(unfiltered_indiv_ct, indiv_include2, indiv_exclude);
  memcpy(indiv_male_include2, indiv_include2, unfiltered_indiv_ctv2 * sizeof(intptr_t));
  vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_male);
  if (chrom_info_ptr->y_code != -1) {
    marker_ct_y = count_chrom_markers(chrom_info_ptr, chrom_info_ptr->y_code, marker_exclude);
  }
  marker_ct_nony = marker_ct - marker_ct_y;
  fill_uint_zero(missing_cts, unfiltered_indiv_ct);
  ujj = unfiltered_indiv_ct2l * BITCT2;
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto write_missingness_reports_ret_READ_FAIL;
  }
  memcpy(outname_end, ".lmiss", 7);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  if (!cluster_ct) {
    sprintf(tbuf, " CHR %%%us   N_MISS   N_GENO   F_MISS\n", plink_maxsnp);
  } else {
    if (wkspace_alloc_ui_checked(&indiv_to_cluster, unfiltered_indiv_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&missing_ct_by_cluster, cluster_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_sizes, cluster_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_sizes_y, cluster_ct * sizeof(int32_t))) {
      goto write_missingness_reports_ret_NOMEM;
    }
    fill_uint_zero(indiv_to_cluster, unfiltered_indiv_ct);
    for (clidx = 0; clidx < cluster_ct; clidx++) {
      unn = cluster_starts[clidx + 1];
      ukk = clidx + 1;
      for (uii = cluster_starts[clidx]; uii < unn; uii++) {
	umm = cluster_map[uii];
	if (!IS_SET(indiv_exclude, umm)) {
          indiv_to_cluster[umm] = ukk;
	  cluster_sizes[clidx] += 1;
          if (IS_SET(sex_male, umm)) {
            cluster_sizes_y[clidx] += 1;
	  }
	}
      }
    }
    sprintf(tbuf, " CHR %%%us       CLST   N_MISS   N_GENO   N_CLST   F_MISS\n", plink_maxsnp);
  }
  fprintf(outfile, tbuf, "SNP");
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    is_x = (((int32_t)chrom_idx) == chrom_info_ptr->x_code);
    is_y = (((int32_t)chrom_idx) == chrom_info_ptr->y_code);
    is_haploid = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
    if (!is_y) {
      cur_nm = indiv_include2;
      cur_tot = indiv_ct;
      cur_cluster_sizes = cluster_sizes;
    } else {
      cur_nm = indiv_male_include2;
      cur_tot = indiv_male_ct;
      cur_cluster_sizes = cluster_sizes_y;
    }
    cptr = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms));
    *cptr++ = ' ';
    if (marker_uidx < chrom_end) {
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto write_missingness_reports_ret_READ_FAIL;
      }
      do {
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto write_missingness_reports_ret_READ_FAIL;
	}
        if (is_haploid) {
          haploid_fix(hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, is_x, is_y, (unsigned char*)loadbuf);
	}
	lptr = loadbuf;
	lptr2 = cur_nm;
	cptr2 = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr);
	*cptr2++ = ' ';
	if (!cluster_ct) {
          ukk = 0;
	  for (uii = 0; uii < ujj; uii += BITCT2) {
	    ulii = *lptr++;
	    ulii = ulii & ((~ulii) >> 1) & (*lptr2++);
	    while (ulii) {
	      missing_cts[uii + (CTZLU(ulii) / 2)] += 1;
	      ukk++;
	      ulii &= ulii - 1;
	    }
	  }
	  wptr = uint32_writew8x(cptr2, ukk, ' ');
          wptr = uint32_writew8x(wptr, cur_tot, ' ');
	  wptr = double_g_writewx4x(wptr, ((double)((int32_t)ukk)) / ((double)((int32_t)cur_tot)), 8, '\n');
          if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto write_missingness_reports_ret_WRITE_FAIL;
	  }
	} else {
	  fill_uint_zero(missing_ct_by_cluster, cluster_ct);
	  for (uii = 0; uii < ujj; uii += BITCT2) {
	    ulii = *lptr++;
	    ulii = ulii & ((~ulii) >> 1) & (*lptr2++);
	    while (ulii) {
	      ukk = uii + (CTZLU(ulii) / 2);
	      missing_cts[ukk] += 1;
	      ukk = indiv_to_cluster[ukk];
	      if (ukk) {
	        missing_ct_by_cluster[ukk - 1] += 1;
	      }
	      ulii &= ulii - 1;
	    }
	  }
	  for (clidx = 0; clidx < cluster_ct; clidx++) {
            wptr = fw_strcpy(10, &(cluster_ids[clidx * max_cluster_id_len]), cptr2);
	    *wptr++ = ' ';
	    uii = missing_ct_by_cluster[clidx];
            wptr = uint32_writew8x(wptr, uii, ' ');
	    umm = cur_cluster_sizes[clidx];
	    wptr = uint32_writew8x(wptr, umm, ' ');
	    wptr = uint32_writew8x(wptr, umm, ' ');
            wptr = double_g_writewx4x(wptr, ((double)((int32_t)uii)) / ((double)((int32_t)umm)), 8, '\n');
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	}
        marker_uidx++;
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	  if (marker_uidx < chrom_end) {
	    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	      goto write_missingness_reports_ret_WRITE_FAIL;
	    }
	  }
	}
      } while (marker_uidx < chrom_end);
    }
  }
  if (fclose_null(&outfile)) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  outname_end[1] = 'i';
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us MISS_PHENO   N_MISS   N_GENO   F_MISS\n", plink_maxfid, plink_maxiid);
  fprintf(outfile, tbuf, "FID", "IID");
  do {
    indiv_uidx = next_unset_unsafe(indiv_exclude, indiv_uidx);
    indiv_uidx_stop = next_set(indiv_exclude, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    do {
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      cptr2 = (char*)memchr(cptr, '\t', max_person_id_len);
      slen = (uintptr_t)(cptr2 - cptr);
      wptr = memseta(tbuf, 32, plink_maxfid - slen);
      wptr = memcpyax(wptr, cptr, slen, ' ');
      wptr = fw_strcpy(plink_maxiid, &(cptr2[1]), wptr);
      wptr = memseta(wptr, 32, 10);
      *wptr++ = 'Y' - (is_set(pheno_nm, indiv_uidx) * 11);
      *wptr++ = ' ';
      uii = missing_cts[indiv_uidx];
      wptr = uint32_writew8x(wptr, uii, ' ');
      ujj = marker_ct_nony + (is_set(sex_male, indiv_uidx) * marker_ct_y);
      wptr = uint32_writew8x(wptr, ujj, ' ');
      wptr = double_g_writewx4x(wptr, ((double)((int32_t)uii)) / ((double)((int32_t)ujj)), 8, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto write_missingness_reports_ret_WRITE_FAIL;
      }
    } while (++indiv_uidx < indiv_uidx_stop);
  } while (indiv_idx < indiv_ct);
  if (fclose_null(&outfile)) {
    goto write_missingness_reports_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  sprintf(logbuf, "--missing: Individual missing data report written to %s.imiss, and\nmarker-based %smissing data report written to %s.lmiss.\n", outname, cluster_ct? "cluster-stratified " : "", outname);
  logprintb();
  while (0) {
  write_missingness_reports_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_missingness_reports_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  write_missingness_reports_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

uintptr_t binary_geno_filter(double geno_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t indiv_ct, uintptr_t male_ct, uint32_t* marker_allele_cts, Chrom_info* chrom_info_ptr) {
  uint32_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uint32_t orig_exclude_ct = marker_exclude_ct;
  uint32_t geno_int_thresh;
  uint32_t marker_uidx;
  uint32_t chrom_end;
  uint32_t chrom_fo_idx;
  int32_t chrom_idx;
  uint32_t cur_ct;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    if (is_set(chrom_info_ptr->haploid_mask, chrom_idx)) {
      if (chrom_idx == chrom_info_ptr->x_code) {
	cur_ct = 2 * indiv_ct - male_ct;
      } else if (chrom_idx == chrom_info_ptr->y_code) {
	cur_ct = male_ct;
      } else {
	cur_ct = indiv_ct;
      }
    } else {
      cur_ct = 2 * indiv_ct;
    }
    geno_int_thresh = cur_ct - (int32_t)(geno_thresh * cur_ct);
    marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    while (marker_uidx < chrom_end) {
      if ((marker_allele_cts[2 * marker_uidx] + marker_allele_cts[2 * marker_uidx + 1]) < geno_int_thresh) {
	SET_BIT(marker_exclude, marker_uidx);
	marker_exclude_ct++;
      }
      marker_uidx++;
      next_unset_ck(marker_exclude, &marker_uidx, chrom_end);
    }
  }
  *marker_exclude_ct_ptr = marker_exclude_ct;
  return (marker_exclude_ct - orig_exclude_ct);
}

void calc_marker_reverse_bin(uintptr_t* marker_reverse, uintptr_t* marker_exclude, uint32_t unfiltered_marker_ct, uint32_t marker_ct, double* set_allele_freqs) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  do {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    for (; marker_uidx < marker_uidx_stop; marker_uidx++) {
      if (set_allele_freqs[marker_uidx] < 0.5) {
	SET_BIT(marker_reverse, marker_uidx);
      }
    }
  } while (markers_done < marker_ct);
}

int32_t hardy_report_write_line(FILE* outfile, char* prefix_buf, uint32_t prefix_len, uint32_t reverse, uint32_t ll_ct, uint32_t lh_ct, uint32_t hh_ct, char* midbuf_ptr, double pvalue) {
  char wbuf[48];
  char* cptr;
  uint32_t denom;
  double drecip;
  double minor_freq;
  fwrite(prefix_buf, 1, prefix_len, outfile);
  if (reverse) {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, hh_ct, '/'), lh_ct, '/'), ll_ct);
  } else {
    cptr = uint32_write(uint32_writex(uint32_writex(wbuf, ll_ct, '/'), lh_ct, '/'), hh_ct);
  }
  cptr = fw_strcpyn(20, cptr - wbuf, wbuf, midbuf_ptr);
  *cptr++ = ' ';
  denom = (ll_ct + lh_ct + hh_ct) * 2;
  if (denom) {
    drecip = 1.0 / ((double)denom);
    minor_freq = (2 * ll_ct + lh_ct) * drecip;
    cptr = double_g_writewx4x(double_g_writewx4x(double_g_writewx4x(cptr, (lh_ct * 2) * drecip, 8, ' '), minor_freq * (2 * hh_ct + lh_ct) * drecip * 2, 8, ' '), pvalue, 12, '\n');
  } else {
    cptr = memcpya(cptr, "     nan      nan           NA\n", 31);
  }
  return fwrite_checked(midbuf_ptr, (cptr - midbuf_ptr), outfile);
}

int32_t hardy_report(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, int32_t* hwe_ll_cases, int32_t* hwe_lh_cases, int32_t* hwe_hh_cases, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, uint32_t pheno_nm_ct, uintptr_t* pheno_c, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr) {
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  int32_t retval = 0;
  uint32_t pct = 0;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t prefix_len;
  uint32_t loop_end;
  uint32_t uii;
  uint32_t report_type;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t reverse;
  double* p_values;
  char* writebuf;
  char* cptr0;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5;
  if (pheno_nm_ct) {
    report_type = pheno_c? 0 : 1;
  } else {
    report_type = 2;
  }
  uii = report_type? 1 : 3;
  if (wkspace_alloc_d_checked(&p_values, uii * marker_ct * sizeof(double)) ||
      wkspace_alloc_c_checked(&writebuf, 2 * max_marker_allele_len + MAXLINELEN)) {
    goto hardy_report_ret_NOMEM;
  }

  // todo: multithread?
  if (report_type) {
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      p_values[marker_uidx] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx]);
    }
  } else {
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      p_values[marker_uidx * 3] = SNPHWE2(hwe_lh_allfs[marker_uidx], hwe_ll_allfs[marker_uidx], hwe_hh_allfs[marker_uidx]);
      p_values[marker_uidx * 3 + 1] = SNPHWE2(hwe_lh_cases[marker_uidx], hwe_ll_cases[marker_uidx], hwe_hh_cases[marker_uidx]);
      p_values[marker_uidx * 3 + 2] = SNPHWE2(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx]);
    }
  }
  marker_uidx = 0;
  marker_idx = 0;

  memcpy(outname_end, ".hwe", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto hardy_report_ret_OPEN_FAIL;
  }
  sprintf(logbuf, "Writing Hardy-Weinberg report to %s...", outname);
  logprintb();
  fputs(" 0%", stdout);
  sprintf(writebuf, " CHR %%%us     TEST   A1   A2                 GENO   O(HET)   E(HET)            P \n", plink_maxsnp);
  fprintf(outfile, writebuf, "SNP");
 
  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
  *cptr0++ = ' ';
  cptr = &(cptr0[10 + plink_maxsnp]);
  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
  if (report_type) {
    if (report_type == 1) {
      memcpy(&(cptr0[plink_maxsnp]), "  ALL(QT)           ", 20);
    } else {
      memcpy(&(cptr0[plink_maxsnp]), "  ALL(NP)           ", 20);
    }
    cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
	  *cptr0++ = ' ';
	  cptr = &(cptr0[10 + plink_maxsnp]);
	  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
	  if (report_type == 1) {
	    memcpy(&(cptr0[plink_maxsnp]), "  ALL(QT)           ", 20);
	  } else {
	    memcpy(&(cptr0[plink_maxsnp]), "  ALL(NP)           ", 20);
	  }
	  cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
	}
	fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr);
	reverse = IS_SET(marker_reverse, marker_uidx);
	cptr3 = marker_allele_ptrs[2 * marker_uidx];
	cptr4 = marker_allele_ptrs[2 * marker_uidx + 1];
	cptr5 = fw_strcpy(4, cptr3, cptr);
	*cptr5 = ' ';
	cptr5 = fw_strcpy(4, cptr4, &(cptr5[1]));
	*cptr5 = ' ';
	prefix_len = 1 + (cptr5 - writebuf);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[marker_uidx])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else {
    memset(&(cptr0[plink_maxsnp]), 32, 20);
    cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
    for (; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100LLU;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	  cptr0 = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, chrom_info_ptr->chrom_file_order[chrom_fo_idx], zero_extra_chroms));
          memset(&(cptr0[plink_maxsnp]), 32, 20);
	  cptr = &(cptr0[10 + plink_maxsnp]);
	  cptr2 = &(cptr[18 + 2 * max_marker_allele_len]);
	  prefix_len = 10 + ((uintptr_t)(cptr - writebuf));
	}
	fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), cptr0);
	memcpy(&(cptr0[4 + plink_maxsnp]), "  ALL", 5);
	reverse = IS_SET(marker_reverse, marker_uidx);
	cptr3 = marker_allele_ptrs[2 * marker_uidx];
	cptr4 = marker_allele_ptrs[2 * marker_uidx + 1];
	cptr5 = fw_strcpy(4, cptr3, cptr);
	*cptr5 = ' ';
	cptr5 = fw_strcpy(4, cptr4, &(cptr5[1]));
	*cptr5 = ' ';
	prefix_len = 1 + (cptr5 - writebuf);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_allfs[marker_uidx], hwe_lh_allfs[marker_uidx], hwe_hh_allfs[marker_uidx], cptr2, p_values[3 * marker_uidx])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[7 + plink_maxsnp]), "FF", 2);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_ll_cases[marker_uidx], hwe_lh_cases[marker_uidx], hwe_hh_cases[marker_uidx], cptr2, p_values[3 * marker_uidx + 1])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}

	memcpy(&(cptr0[4 + plink_maxsnp]), "UN", 2);
	if (hardy_report_write_line(outfile, writebuf, prefix_len, reverse, hwe_lls[marker_uidx], hwe_lhs[marker_uidx], hwe_hhs[marker_uidx], cptr2, p_values[3 * marker_uidx + 2])) {
	  goto hardy_report_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  }
  fputs("\b\b\b\b", stdout);
  logprint(" done.\n");

  while (0) {
  hardy_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  hardy_report_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  hardy_report_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}


void enforce_hwe_threshold(double hwe_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_all, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t markers_done;
  hwe_thresh += EPSILON;
  if (hwe_all) {
    hwe_lhs = hwe_lh_allfs;
    hwe_lls = hwe_ll_allfs;
    hwe_hhs = hwe_hh_allfs;
  }
  for (markers_done = 0; markers_done < marker_ct; marker_uidx++, markers_done++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (SNPHWE_t(hwe_lhs[marker_uidx], hwe_lls[marker_uidx], hwe_hhs[marker_uidx], hwe_thresh)) {
      SET_BIT(marker_exclude, marker_uidx);
      removed_ct++;
    }
  }
  sprintf(logbuf, "%u marker%s removed due to Hardy-Weinberg exact test (--hwe).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
}

void enforce_maf_threshold(double min_maf, double max_maf, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, double* set_allele_freqs) {
  uint32_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uint32_t marker_uidx = 0;
  uint32_t removed_ct = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  double dxx;
  min_maf -= EPSILON;
  max_maf += EPSILON;
  while (markers_done < marker_ct) {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      dxx = get_maf(set_allele_freqs[marker_uidx]);
      if ((dxx < min_maf) || (dxx > max_maf)) {
        SET_BIT(marker_exclude, marker_uidx);
        removed_ct++;
      }
    } while (++marker_uidx < marker_uidx_stop);
  }
  sprintf(logbuf, "%u marker%s removed due to MAF threshold(s) (--maf/--max-maf).\n", removed_ct, (removed_ct == 1)? "" : "s");
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
}

void enforce_min_bp_space(int32_t min_bp_space, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, Chrom_info* chrom_info_ptr) {
  uint32_t marker_ct = unfiltered_marker_ct - (uint32_t)(*marker_exclude_ct_ptr);
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t removed_ct = 0;
  uint32_t chrom_end = 0;
  uint32_t marker_uidx = next_unset(marker_exclude, 0, unfiltered_marker_ct);
  uint32_t chrom_fo_idx_p1 = 0;
  uint32_t marker_uidx_stop;
  int32_t last_pos;
  int32_t cur_pos;
  for (chrom_fo_idx_p1 = 1; chrom_fo_idx_p1 <= chrom_ct; chrom_fo_idx_p1++) {
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx_p1];
    if (marker_uidx >= chrom_end) {
      continue;
    }
    last_pos = -2147483647;
    do {
      marker_uidx_stop = next_set(marker_exclude, marker_uidx, chrom_end);
      do {
        cur_pos = marker_pos[marker_uidx];
        if (cur_pos < last_pos + min_bp_space) {
          SET_BIT(marker_exclude, marker_uidx);
	  removed_ct++;
	} else {
	  last_pos = cur_pos;
	}
      } while (++marker_uidx < marker_uidx_stop);
      marker_uidx = next_unset(marker_exclude, marker_uidx, unfiltered_marker_ct);
    } while (marker_uidx < chrom_end);
  }
  sprintf(logbuf, "--bp-space: %u marker%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", marker_ct - removed_ct);
  logprintb();
  *marker_exclude_ct_ptr += removed_ct;
}

void calc_marker_weights(double exponent, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, double* marker_weights) {
  uint32_t marker_uidx = 0;
  uint32_t markers_done = 0;
  uint32_t marker_uidx_stop;
  do {
    marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
    marker_uidx_stop = next_set(marker_exclude, marker_uidx, unfiltered_marker_ct);
    markers_done += marker_uidx_stop - marker_uidx;
    do {
      if (marker_weights[marker_uidx] < 0.0) {
	marker_weights[marker_uidx] = calc_wt_mean(exponent, lh_cts[marker_uidx], ll_cts[marker_uidx], hh_cts[marker_uidx]);
      }
    } while (++marker_uidx < marker_uidx_stop);
  } while (markers_done < marker_ct);
}

int32_t load_ax_alleles(Two_col_params* axalleles, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t is_a2) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = axalleles->skipchar;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  uint32_t colid_first = (axalleles->colid < axalleles->colx);
  uintptr_t marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t slen;
  char* sorted_marker_ids;
  char* colid_ptr;
  char* colx_ptr;
  uint32_t* marker_id_map;
  uint32_t colmin;
  uint32_t coldiff;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  char cc;
  int32_t retval;
  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto load_ax_alleles_ret_1;
  }
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto load_ax_alleles_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_ax_alleles_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  retval = open_and_skip_first_lines(&infile, axalleles->fname, loadbuf, loadbuf_size, axalleles->skip);
  if (colid_first) {
    colmin = axalleles->colid - 1;
    coldiff = axalleles->colx - axalleles->colid;
  } else {
    colmin = axalleles->colx - 1;
    coldiff = axalleles->colid - axalleles->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Pathologically long line in %s.\n", axalleles->fname);
        logprintb();
	goto load_ax_alleles_ret_INVALID_FORMAT;
      } else {
        goto load_ax_alleles_ret_NOMEM;
      }
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      if (colmin) {
        colid_ptr = next_item_mult(colid_ptr, colmin);
      }
      colx_ptr = next_item_mult(colid_ptr, coldiff);
    } else {
      colx_ptr = colid_ptr;
      if (colmin) {
	colx_ptr = next_item_mult(colx_ptr, colmin);
      }
      colid_ptr = next_item_mult(colx_ptr, coldiff);
    }
    colid_ptr[strlen_se(colid_ptr)] = '\0';
    sorted_idx = bsearch_str(colid_ptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate marker %s in --a%c-allele file.\n", colid_ptr, is_a2? '2' : '1');
      logprintb();
      goto load_ax_alleles_ret_INVALID_FORMAT;
    }
    SET_BIT(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    slen = strlen_se(colx_ptr);
    colx_ptr[slen] = '\0';
    if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + is_a2])) {
      CLEAR_BIT(marker_reverse, marker_uidx);
    } else if (!strcmp(colx_ptr, marker_allele_ptrs[marker_uidx * 2 + 1 - is_a2])) {
      SET_BIT(marker_reverse, marker_uidx);
    } else if (marker_allele_ptrs[marker_uidx * 2 + is_a2] == missing_geno_ptr) {
      if (allele_reset(&(marker_allele_ptrs[marker_uidx * 2 + is_a2]), colx_ptr, slen)) {
	goto load_ax_alleles_ret_NOMEM;
      }
      if (slen >= max_marker_allele_len) {
	max_marker_allele_len = slen + 1;
      }
      CLEAR_BIT(marker_reverse, marker_uidx);
    } else {
      sprintf(logbuf, "Warning: Impossible A%c allele assignment for marker %s.\n", is_a2? '2' : '1', colid_ptr);
      logprintb();
    }
  }
  if (!feof(infile)) {
    goto load_ax_alleles_ret_READ_FAIL;
  }
  *max_marker_allele_len_ptr = max_marker_allele_len;
  while (0) {
  load_ax_alleles_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_ax_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_ax_alleles_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_ax_alleles_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

void swap_reversed_marker_alleles(uintptr_t unfiltered_marker_ct, uintptr_t* marker_reverse, char** marker_allele_ptrs) {
  uintptr_t marker_uidx = 0;
  char* swap_ptr;
  while (1) {
    next_set_ul_ck(marker_reverse, &marker_uidx, unfiltered_marker_ct);
    if (marker_uidx == unfiltered_marker_ct) {
      return;
    }
    swap_ptr = marker_allele_ptrs[marker_uidx * 2];
    marker_allele_ptrs[marker_uidx * 2] = marker_allele_ptrs[marker_uidx * 2 + 1];
    marker_allele_ptrs[marker_uidx * 2 + 1] = swap_ptr;
    marker_uidx++;
  }
}

int32_t write_snplist(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t list_23_indels) {
  FILE* outfile;
  uintptr_t marker_uidx = 0;
  uintptr_t markers_done = 0;
  int32_t retval = 0;
  const char* a0ptr = g_missing_geno_ptr;
  const char* adptr = &(g_one_char_strs[136]); // "D"
  const char* aiptr = &(g_one_char_strs[146]); // "I"
  char* a1ptr;
  char* a2ptr;
  char* cptr;
  char* cptr_end;
  uintptr_t marker_uidx_stop;
  if (!list_23_indels) {
    memcpy(outname_end, ".snplist", 9);
  } else {
    memcpy(outname_end, ".indel", 7);
  }
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_snplist_ret_OPEN_FAIL;
  }
  if (!list_23_indels) {
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      markers_done += marker_uidx_stop - marker_uidx;
      cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      cptr_end = &(marker_ids[marker_uidx_stop * max_marker_id_len]);
      marker_uidx = marker_uidx_stop;
      do {
	fputs(cptr, outfile);
	if (putc_checked('\n', outfile)) {
          goto write_snplist_ret_WRITE_FAIL;
	}
        cptr = &(cptr[max_marker_id_len]);
      } while (cptr < cptr_end);
    } while (markers_done < marker_ct);
  } else {
    for (; markers_done < marker_ct; marker_uidx++, markers_done++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      a1ptr = marker_allele_ptrs[2 * marker_uidx];
      a2ptr = marker_allele_ptrs[2 * marker_uidx + 1];
      if ((a1ptr != adptr) && (a1ptr != aiptr)) {
        if ((a1ptr != a0ptr) || ((a2ptr != adptr) && (a2ptr != aiptr))) {
	  continue;
	}
      } else if ((a2ptr != adptr) && (a2ptr != aiptr) && (a2ptr != a0ptr)) {
	continue;
      }
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      if (putc_checked('\n', outfile)) {
	goto write_snplist_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto write_snplist_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "List of %smarker IDs written to %s.\n", list_23_indels? "indel " : "" , outname);
  logprintb();
  while (0) {
  write_snplist_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_snplist_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

static inline uint32_t are_marker_pos_needed(uint64_t calculation_type, uint32_t min_bp_space, uint32_t genome_skip_write) {
  return (calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_GENOME | CALC_HOMOZYG | CALC_LD_PRUNE | CALC_REGRESS_PCS | CALC_MODEL | CALC_GLM)) || min_bp_space || genome_skip_write;
}

static inline uint32_t are_marker_cms_needed(uint64_t calculation_type, Two_col_params* update_cm) {
  if (calculation_type & (CALC_MAKE_BED | CALC_RECODE)) {
    if (update_cm) {
      return MARKER_CMS_FORCED;
    } else {
      return MARKER_CMS_OPTIONAL;
    }
  }
  return 0;
}

static inline uint32_t are_marker_alleles_needed(uint64_t calculation_type, char* freqname, Homozyg_info* homozyg_ptr, Two_col_params* a1alleles, Two_col_params* a2alleles) {
  return (freqname || (calculation_type & (CALC_FREQ | CALC_HARDY | CALC_MAKE_BED | CALC_RECODE | CALC_REGRESS_PCS | CALC_MODEL | CALC_GLM | CALC_LASSO | CALC_LIST_23_INDELS)) || ((calculation_type & CALC_HOMOZYG) && (homozyg_ptr->modifier & HOMOZYG_GROUP_VERBOSE)) || a1alleles || a2alleles);
}

inline int32_t relationship_or_ibc_req(uint64_t calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

inline int32_t distance_wt_req(uint64_t calculation_type, char* read_dists_fname, uint32_t dist_calc_type) {
  return (((calculation_type & CALC_DISTANCE) || ((!read_dists_fname) && ((calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))))) && (!(dist_calc_type & DISTANCE_FLAT_MISSING)));
}

int32_t wdist(char* outname, char* outname_end, char* pedname, char* mapname, char* famname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* keepfamname, char* removefamname, char* filtername, char* freqname, char* read_dists_fname, char* read_dists_id_fname, char* evecname, char* mergename1, char* mergename2, char* mergename3, char* makepheno_str, char* phenoname_str, Two_col_params* a1alleles, Two_col_params* a2alleles, char* recode_allele_name, char* covar_fname, char* set_fname, char* subset_fname, char* update_alleles_fname, char* read_genome_fname, Two_col_params* update_chr, Two_col_params* update_cm, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* loop_assoc_fname, char* flip_fname, char* flip_subset_fname, char* filtervals_flattened, char* condition_mname, char* condition_fname, double thin_keep_prob, uint32_t min_bp_space, uint32_t mfilter_col, uint32_t filter_binary, uint32_t fam_cols, int32_t missing_pheno, char* output_missing_pheno, uint32_t mpheno_col, uint32_t pheno_modifier, Chrom_info* chrom_info_ptr, double exponent, double min_maf, double max_maf, double geno_thresh, double mind_thresh, double hwe_thresh, double rel_cutoff, double tail_bottom, double tail_top, uint64_t misc_flags, uint64_t calculation_type, uint32_t rel_calc_type, uint32_t dist_calc_type, uintptr_t groupdist_iters, uint32_t groupdist_d, uintptr_t regress_iters, uint32_t regress_d, uintptr_t regress_rel_iters, uint32_t regress_rel_d, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_covr, int32_t ibc_type, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t ppc_gap, uint32_t sex_missing_pheno, uint32_t genome_modifier, double genome_min_pi_hat, double genome_max_pi_hat, Homozyg_info* homozyg_ptr, Cluster_info* cluster_ptr, uint32_t neighbor_n1, uint32_t neighbor_n2, uint32_t ld_window_size, uint32_t ld_window_kb, uint32_t ld_window_incr, double ld_last_param, uint32_t regress_pcs_modifier, uint32_t max_pcs, uint32_t recode_modifier, uint32_t allelexxxx, uint32_t merge_type, uint32_t indiv_sort, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* snps_range_list_ptr, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uint32_t mwithin_col, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uint32_t gxe_mcovar, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t mperm_save, uint32_t ibs_test_perms, uint32_t perm_batch_size, double lasso_h2, double lasso_minlambda, Ll_str** file_delete_list_ptr) {
  FILE* bedfile = NULL;
  FILE* famfile = NULL;
  FILE* phenofile = NULL;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t* marker_exclude = NULL;
  uintptr_t max_marker_id_len = 0;
  // set_allele_freqs = frequency of allele corresponding to set bits in .bed
  //   (i.e. A2), or frequency of MAJOR allele in middle of text loading.
  double* set_allele_freqs = NULL;
  uintptr_t unfiltered_indiv_ct = 0;
  uintptr_t* indiv_exclude = NULL;
  uintptr_t indiv_exclude_ct = 0;
  uint32_t* indiv_sort_map = NULL;
  uintptr_t* founder_info = NULL;
  uintptr_t* sex_nm = NULL;
  uintptr_t* sex_male = NULL;
  uint32_t genome_skip_write = (cluster_ptr->ppc != 0.0) && (!(calculation_type & CALC_GENOME)) && (!read_genome_fname);
  uint32_t marker_pos_needed = are_marker_pos_needed(calculation_type, min_bp_space, genome_skip_write);
  uint32_t marker_cms_needed = are_marker_cms_needed(calculation_type, update_cm);
  uint32_t marker_alleles_needed = are_marker_alleles_needed(calculation_type, freqname, homozyg_ptr, a1alleles, a2alleles);
  uint32_t zero_extra_chroms = (misc_flags / MISC_ZERO_EXTRA_CHROMS) & 1;
  uint32_t uii = 0;
  int64_t llxx = 0;
  uint32_t nonfounders = (misc_flags / MISC_NONFOUNDERS) & 1;
  uint32_t pheno_all = pheno_modifier & PHENO_ALL;
  char* marker_ids = NULL;
  double* marker_cms = NULL;
  unsigned char* marker_weights_base = NULL;
  // marker_allele_ptrs[2 * i] is id of A1 (usually minor) allele at marker i
  // marker_allele_ptrs[2 * i + 1] is id of A2 allele
  // Single-character allele names point to g_one_char_strs[]; otherwise
  // string allocation occurs on the heap.
  char** marker_allele_ptrs = NULL;
  uintptr_t max_marker_allele_len = 2; // includes trailing null
  uintptr_t* marker_reverse = NULL;
  int32_t retval = 0;
  uint32_t map_is_unsorted = 0;
  uint32_t map_cols = 3;
  uint32_t affection = 0;
  uintptr_t* pheno_nm = NULL;
  uintptr_t* pheno_nosex_exclude = NULL; // --make-bed/--recode only
  uintptr_t* orig_pheno_nm = NULL; // --all-pheno + --pheno-merge
  uintptr_t* pheno_c = NULL;
  uintptr_t* orig_pheno_c = NULL;
  double* pheno_d = NULL;
  double* orig_pheno_d = NULL;
  double* marker_weights = NULL;
  uint32_t marker_weight_sum = 0;
  uint32_t* marker_weights_i = NULL;
  char* person_ids = NULL;
  uintptr_t max_person_id_len = 4;
  char* paternal_ids = NULL;
  uintptr_t max_paternal_id_len = 2;
  char* maternal_ids = NULL;
  uintptr_t max_maternal_id_len = 2;
  unsigned char* wkspace_mark = NULL;
  uintptr_t cluster_ct = 0;
  uint32_t* cluster_map = NULL; // unfiltered indiv IDs
  // index for cluster_map, length (cluster_ct + 1)
  // cluster_starts[n+1] - cluster_starts[n] = length of cluster n (0-based)
  uint32_t* cluster_starts = NULL;
  char* cluster_ids = NULL;
  uintptr_t max_cluster_id_len = 2;
  double* mds_plot_dmatrix_copy = NULL;
  uintptr_t* cluster_merge_prevented = NULL;
  double* cluster_sorted_ibs = NULL;
  char* cptr = NULL;
  uint64_t dists_alloc = 0;
  uintptr_t marker_exclude_ct = 0;
  char* pid_list = NULL;
  char* id_list = NULL;
  double missing_phenod = (double)missing_pheno;
  double ci_zt = 0.0;
  uint32_t missing_pheno_len = intlen(missing_pheno);
  uint32_t wt_needed = distance_wt_req(calculation_type, read_dists_fname, dist_calc_type);
  uintptr_t bed_offset = 3;
  uint32_t* marker_pos = NULL;
  uint32_t hh_exists = 0;
  uint32_t pheno_ctrl_ct = 0;
  uintptr_t covar_ct = 0;
  char* covar_names = NULL;
  uintptr_t max_covar_name_len = 0;
  uintptr_t* covar_nm = NULL;
  double* covar_d = NULL;
  uintptr_t* gxe_covar_nm = NULL;
  uintptr_t* gxe_covar_c = NULL;
  uint32_t plink_maxfid = 0;
  uint32_t plink_maxiid = 0;
  unsigned char* wkspace_mark2 = NULL;
  unsigned char* wkspace_mark_precluster = NULL;
  unsigned char* wkspace_mark_postcluster = NULL;
  pthread_t threads[MAX_THREADS];
  uintptr_t unfiltered_indiv_ct4;
  uintptr_t unfiltered_indiv_ctl;
  uint32_t* marker_allele_cts;
  uint32_t* uiptr;
  double* dptr;
  double* dptr2;
  double* rel_ibc;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t ujj;
  uint32_t ukk;
  double dxx;
  char* outname_end2;
  uintptr_t marker_ct;
  uint32_t plink_maxsnp;
  int32_t ii;
  int64_t llyy;
  int32_t* hwe_lls;
  int32_t* hwe_lhs;
  int32_t* hwe_hhs;
  int32_t* hwe_ll_cases;
  int32_t* hwe_lh_cases;
  int32_t* hwe_hh_cases;
  int32_t* hwe_ll_allfs;
  int32_t* hwe_lh_allfs;
  int32_t* hwe_hh_allfs;
  int32_t* hwe_hapl_allfs;
  int32_t* hwe_haph_allfs;
  uint32_t indiv_male_ct;
  uint32_t indiv_f_ct;
  uint32_t indiv_f_male_ct;
  uint32_t pheno_nm_ct;
  Pedigree_rel_info pri;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx_stop;
  uintptr_t marker_idx;
  uint32_t gender_unk_ct;

  if (update_cm && (!marker_cms_needed)) {
    logprint("Error: --update-cm results would never be used.  (Did you forget --make-bed?)\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if (update_map && (!marker_pos_needed)) {
    logprint("Error: --update-map results would never be used.  (Did you forget --make-bed?)\n");
    goto wdist_ret_INVALID_CMDLINE;
  }
  if (ci_size != 0.0) {
    ci_zt = ltqnorm(1 - (1 - ci_size) / 2);
  }
  if (rel_calc_type & REL_CALC_COV) {
    ibc_type = -1;
  }

  if (calculation_type & CALC_MAKE_BED) {
#if _WIN32
    uii = GetFullPathName(pedname, FNAMESIZE, tbuf, NULL);
    if ((!uii) || (uii > FNAMESIZE))
#else
    if (!realpath(pedname, tbuf))
#endif
    {
      sprintf(logbuf, "Error: Failed to open %s.\n", pedname);
      logprintb();
      goto wdist_ret_OPEN_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    // if file doesn't exist, realpath returns NULL on Linux instead of what
    // the path would be.
#if _WIN32
    uii = GetFullPathName(outname, FNAMESIZE, &(tbuf[FNAMESIZE + 64]), NULL);
    if (uii && (uii <= FNAMESIZE) && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64]))))
#else
    cptr = realpath(outname, &(tbuf[FNAMESIZE + 64]));
    if (cptr && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64]))))
#endif
    {
      logprint("Note: --make-bed input and output filenames match.  Appending '~' to input\nfilenames.\n");
      uii = strlen(pedname);
      memcpy(tbuf, pedname, uii + 1);
      memcpy(&(pedname[uii]), "~", 2);
      if (rename(tbuf, pedname)) {
	logprint("Error: Failed to append '~' to input .bed filename.\n");
	goto wdist_ret_OPEN_FAIL;
      }
      uii = strlen(mapname);
      memcpy(tbuf, mapname, uii + 1);
      memcpy(&(mapname[uii]), "~", 2);
      if (rename(tbuf, mapname)) {
	logprint("Error: Failed to append '~' to input .bim filename.\n");
	goto wdist_ret_OPEN_FAIL;
      }
      uii = strlen(famname);
      memcpy(tbuf, famname, uii + 1);
      memcpy(&(famname[uii]), "~", 2);
      if (rename(tbuf, famname)) {
	logprint("Error: Failed to append '~' to input .fam filename.\n");
	goto wdist_ret_OPEN_FAIL;
      }
    }
  }

  if (calculation_type & CALC_MERGE) {
    if (!(((fam_cols & FAM_COL_13456) == FAM_COL_13456) && (!(misc_flags & MISC_AFFECTION_01)) && (missing_pheno == -9))) {
      logprint("Error: --merge/--bmerge/--merge-list cannot be used with an irregularly\nformatted reference fileset (--no-fid, --no-parents, --no-sex, --no-pheno,\n--1).  Use --make-bed first.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    // Only append -merge to the filename stem if --make-bed or --recode lgen
    // is specified.
    ulii = bed_suffix_conflict(calculation_type, recode_modifier);
    if (ulii) {
      memcpy(outname_end, "-merge", 7);
    }
    retval = merge_datasets(pedname, mapname, famname, outname, ulii? &(outname_end[6]) : outname_end, mergename1, mergename2, mergename3, calculation_type, merge_type, indiv_sort, misc_flags, chrom_info_ptr);
    if (retval || (!(calculation_type & (~CALC_MERGE)))) {
      goto wdist_ret_1;
    }
    uljj = (uintptr_t)(outname_end - outname) + (ulii? 6 : 0);
    memcpy(memcpya(pedname, outname, uljj), ".bed", 5);
    memcpy(memcpya(famname, pedname, uljj), ".fam", 5);
    memcpy(memcpya(mapname, pedname, uljj), ".bim", 5);
    if ((calculation_type & CALC_MAKE_BED) && ulii) {
      if (push_ll_str(file_delete_list_ptr, pedname) || push_ll_str(file_delete_list_ptr, famname) || push_ll_str(file_delete_list_ptr, mapname)) {
	goto wdist_ret_NOMEM;
      }
    }
  }

  if (fopen_checked(&bedfile, pedname, "rb")) {
    goto wdist_ret_OPEN_FAIL;
  }
  if (fopen_checked(&famfile, famname, "rb")) {
    goto wdist_ret_OPEN_FAIL;
  }
  // load .bim, count markers, filter chromosomes
  if (update_name) {
    ulii = 0;
    retval = scan_max_strlen(update_name->fname, update_name->colid, update_name->colx, update_name->skip, update_name->skipchar, &max_marker_id_len, &ulii);
    if (retval) {
      goto wdist_ret_1;
    }
    if (ulii > 80) {
      // only warn on long new marker ID, since if there's a long old marker ID
      // and no long new one, it's reasonable to infer that the user is fixing
      // the problem, so we shouldn't spam them.
      logprint("Warning: Unusually long new marker ID(s) in --update-name file.  Double-check\nyour file and command-line parameters, and consider changing your naming\nscheme if you encounter memory problems.\n");
    }
    if (ulii > max_marker_id_len) {
      max_marker_id_len = ulii;
    }
  }
  if (!marker_alleles_needed) {
    allelexxxx = 0;
  }
  retval = load_bim(mapname, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &marker_exclude, &set_allele_freqs, &marker_allele_ptrs, &max_marker_allele_len, &marker_ids, (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1, chrom_info_ptr, &marker_cms, &marker_pos, freqname, calculation_type, recode_modifier, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, (misc_flags / MISC_EXCLUDE_MARKERNAME_SNP) & 1, snps_range_list_ptr, &map_is_unsorted, marker_pos_needed, marker_cms_needed, marker_alleles_needed, "bim", ((calculation_type == CALC_MAKE_BED) && (mind_thresh == 1.0) && (geno_thresh == 1.0) && (!update_map) && freqname)? NULL : "make-bed");
  if (retval) {
    goto wdist_ret_1;
  }

  // load .fam, count indivs
  uii = fam_cols & FAM_COL_6;
  if (uii && phenoname) {
    uii = (pheno_modifier & PHENO_MERGE) && (!makepheno_str);
  }
  if (update_ids_fname) {
    ulii = 0;
    retval = scan_max_fam_indiv_strlen(update_ids_fname, 3, &max_person_id_len);
    if (retval) {
      goto wdist_ret_1;
    }
  } else if (update_parents_fname) {
    retval = scan_max_strlen(update_parents_fname, 3, 4, 0, '\0', &max_paternal_id_len, &max_maternal_id_len);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  retval = load_fam(famfile, MAXLINELEN, fam_cols, uii, missing_pheno, missing_pheno_len, (misc_flags / MISC_AFFECTION_01) & 1, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &indiv_exclude);
  if (retval) {
    goto wdist_ret_1;
  }

  unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;

  if (misc_flags & MISC_MAKE_FOUNDERS_FIRST) {
    if (make_founders(unfiltered_indiv_ct, unfiltered_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, indiv_exclude, founder_info)) {
      goto wdist_ret_NOMEM;
    }
  }

  if (pheno_modifier & PHENO_MERGE) {
    if ((!pheno_c) && (!pheno_d)) {
      pheno_modifier &= ~PHENO_MERGE; // nothing to merge
    } else if (pheno_all) {
      orig_pheno_nm = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      if (!orig_pheno_nm) {
	goto wdist_ret_NOMEM;
      }
      memcpy(orig_pheno_nm, pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
      if (pheno_c) {
	orig_pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
	if (!orig_pheno_c) {
	  goto wdist_ret_NOMEM;
	}
	memcpy(orig_pheno_c, pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t));
      } else if (pheno_d) {
	orig_pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
	if (!orig_pheno_d) {
	  goto wdist_ret_NOMEM;
	}
	memcpy(orig_pheno_d, pheno_d, unfiltered_indiv_ct * sizeof(double));
      }
    }
  }
  count_genders(sex_nm, sex_male, unfiltered_indiv_ct, indiv_exclude, &uii, &ujj, &gender_unk_ct);
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (gender_unk_ct) {
    sprintf(logbuf, "%" PRIuPTR " marker%s and %" PRIuPTR " %s (%d male%s, %d female%s, %u ambiguous) loaded.\n", marker_ct, (marker_ct == 1)? "" : "s", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s", gender_unk_ct);
  } else {
    sprintf(logbuf, "%" PRIuPTR " marker%s and %" PRIuPTR " %s (%d male%s, %d female%s) loaded.\n", marker_ct, (marker_ct == 1)? "" : "s", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s");
  }
  logprintb();

  if (phenoname && fopen_checked(&phenofile, phenoname, "r")) {
    goto wdist_ret_OPEN_FAIL;
  }

  if (phenofile || update_ids_fname || update_parents_fname || update_sex_fname || (misc_flags & MISC_TAIL_PHENO)) {
    wkspace_mark = wkspace_base;
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, 0, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto wdist_ret_1;
    }

    if (makepheno_str) {
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_indiv_ct, cptr, max_person_id_len, uiptr, pheno_nm, &pheno_c);
      if (retval) {
	goto wdist_ret_1;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_indiv_ct, 0, cptr, max_person_id_len, uiptr, missing_pheno, missing_pheno_len, (misc_flags / MISC_AFFECTION_01) & 1, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d);
      if (retval) {
	if (retval == LOAD_PHENO_LAST_COL) {
	  logprint(errstr_phenotype_format);
	  logprint("Fewer tokens than expected in line.\n");
	  retval = RET_INVALID_FORMAT;
	  wkspace_reset(wkspace_mark);
	}
	goto wdist_ret_1;
      }
    }
    if (misc_flags & MISC_TAIL_PHENO) {
      retval = convert_tail_pheno(unfiltered_indiv_ct, pheno_nm, &pheno_c, &pheno_d, tail_bottom, tail_top, missing_phenod);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if ((calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST)) && (!pheno_c)) {
    logprint("Error: --ibs-test and --groupdist calculations require a case/control\nphenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_REGRESS_DISTANCE) && (!pheno_d)) {
    logprint("Error: --regress-distance calculation requires a scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_UNRELATED_HERITABILITY) && (!pheno_d)) {
    logprint("Error: --unrelated-heritability requires a scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & (CALC_REGRESS_PCS | CALC_REGRESS_PCS_DISTANCE)) && (!pheno_d)) {
    sprintf(logbuf, "Error: --regress-pcs%s requires a scalar phenotype.\n", (calculation_type & CALC_REGRESS_PCS_DISTANCE)? "-distance" : "");
    goto wdist_ret_INVALID_CMDLINE_2;
  } else if ((calculation_type & CALC_MODEL) && (!(model_modifier & MODEL_ASSOC)) && (!pheno_c)) {
    logprint("Error: --model requires a case/control phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_GXE) && (!pheno_d)) {
    logprint("Error: --gxe requires a scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_LASSO) && (!pheno_d)) {
    logprint("Error: --lasso requires a scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_CMH) && (!pheno_c)) {
    logprint("Error: --mh and --mh2 require a case/control phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_HOMOG) && (!pheno_c)) {
    logprint("Error: --homog requires a case/control phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_CLUSTER) && (!pheno_c)) {
    if (cluster_ptr->modifier & CLUSTER_CC) {
      logprint("Error: --cc requires a case/control phenotype.\n");
      goto wdist_ret_INVALID_CMDLINE;
    } else if ((cluster_ptr->max_cases != 0xffffffffU) || (cluster_ptr->max_ctrls != 0xffffffffU)) {
      logprint("Error: --mcc requires a case/control phenotype.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
  } else if ((calculation_type & CALC_GLM) && (!(pheno_modifier & PHENO_ALL))) {
    if (glm_modifier & GLM_LOGISTIC) {
      if (!pheno_c) {
	logprint("Error: --logistic without --all-pheno requires a case/control phenotype.\n");
	goto wdist_ret_INVALID_CMDLINE;
      }
    } else if (!pheno_d) {
      logprint("Error: --linear without --all-pheno requires a scalar phenotype.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
  }

  uii = update_cm || update_map || update_name || (marker_alleles_needed && (update_alleles_fname || (flip_fname && (!flip_subset_fname))));
  if (uii || extractname || excludename) {
    wkspace_mark = wkspace_base;
    // only permit duplicate marker IDs for --extract/--exclude
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, !uii, 0, strcmp_deref);
    if (retval) {
      goto wdist_ret_1;
    }
    // length of sorted list is NOT necessarily equal to unfiltered_marker_ct -
    // marker_exclude_ct for --exclude, since marker_exclude_ct may first
    // change from --update-map or --extract
    ulii = unfiltered_marker_ct - marker_exclude_ct;

    if (update_cm) {
      retval = update_marker_cms(update_cm, cptr, ulii, max_marker_id_len, uiptr, marker_cms);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (update_map) {
      retval = update_marker_pos(update_map, cptr, ulii, max_marker_id_len, uiptr, marker_exclude, &marker_exclude_ct, marker_pos, &map_is_unsorted, chrom_info_ptr);
    }
    if (update_name) {
      retval = update_marker_names(update_name, cptr, ulii, max_marker_id_len, uiptr, marker_ids);
      if (retval) {
	goto wdist_ret_1;
      }
      if (update_alleles_fname || (marker_alleles_needed && flip_fname && (!flip_subset_fname)) || extractname || excludename) {
	wkspace_reset(wkspace_mark);
	retval = sort_item_ids(&cptr, &uiptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
	if (retval) {
	  goto wdist_ret_1;
	}
	ulii = unfiltered_marker_ct - marker_exclude_ct;
      }
    }
    if (marker_alleles_needed) {
      if (update_alleles_fname) {
        retval = update_marker_alleles(update_alleles_fname, cptr, ulii, max_marker_id_len, uiptr, marker_allele_ptrs, &max_marker_allele_len, outname, outname_end);
        if (retval) {
	  goto wdist_ret_1;
        }
      }
      if (flip_fname && (!flip_subset_fname)) {
        retval = flip_strand(flip_fname, cptr, ulii, max_marker_id_len, uiptr, marker_allele_ptrs);
        if (retval) {
	  goto wdist_ret_1;
        }
      }
    }
    if (extractname) {
      retval = include_or_exclude(extractname, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
      if (retval) {
        goto wdist_ret_1;
      }
    }
    if (excludename) {
      retval = include_or_exclude(excludename, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1);
      if (retval) {
        goto wdist_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }

  if (allelexxxx) {
    allelexxxx_recode(allelexxxx, marker_allele_ptrs, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct);
  }

  if (update_ids_fname || update_parents_fname || update_sex_fname || keepname || keepfamname || removename || removefamname || filtername) {
    wkspace_mark = wkspace_base;
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto wdist_ret_1;
    }
    ulii = unfiltered_indiv_ct - indiv_exclude_ct;
    if (update_ids_fname) {
      retval = update_indiv_ids(update_ids_fname, cptr, ulii, max_person_id_len, uiptr, person_ids);
      if (retval) {
	goto wdist_ret_1;
      }
    } else {
      if (update_parents_fname) {
	retval = update_indiv_parents(update_parents_fname, cptr, ulii, max_person_id_len, uiptr, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
	if (retval) {
	  goto wdist_ret_1;
	}
      }
      if (update_sex_fname) {
        retval = update_indiv_sexes(update_sex_fname, cptr, ulii, max_person_id_len, uiptr, sex_nm, sex_male);
	if (retval) {
	  goto wdist_ret_1;
	}
      }
    }
    if (keepfamname) {
      retval = include_or_exclude(keepfamname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 6);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (keepname) {
      retval = include_or_exclude(keepname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 2);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (removefamname) {
      retval = include_or_exclude(removefamname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 7);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (removename) {
      retval = include_or_exclude(removename, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 3);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_indivs_file(filtername, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, filtervals_flattened, mfilter_col);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }
  if (thin_keep_prob != 1.0) {
    random_thin_markers(thin_keep_prob, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
  }

  if (gender_unk_ct && (calculation_type & (CALC_MAKE_BED | CALC_RECODE)) && (!sex_missing_pheno)) {
    // postponed to here to avoid e.g. --update-sex not working on missing-sex
    // individuals because --allow-no-sex was not specified
    pheno_nosex_exclude = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
    if (!pheno_nosex_exclude) {
      goto wdist_ret_NOMEM;
    }
    memcpy(pheno_nosex_exclude, pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
    bitfield_andnot(pheno_nosex_exclude, indiv_exclude, unfiltered_indiv_ctl);
    bitfield_andnot(pheno_nosex_exclude, sex_nm, unfiltered_indiv_ctl);
  }
  if (!(sex_missing_pheno & ALLOW_NO_SEX)) {
    bitfield_and(pheno_nm, sex_nm, unfiltered_indiv_ctl);
  }
  if (misc_flags & MISC_PRUNE) {
    bitfield_ornot(indiv_exclude, pheno_nm, unfiltered_indiv_ctl);
    zero_trailing_bits(indiv_exclude, unfiltered_indiv_ct);
    indiv_exclude_ct = popcount_longs(indiv_exclude, 0, unfiltered_indiv_ctl);
  }

  if (filter_binary & (FILTER_BINARY_CASES | FILTER_BINARY_CONTROLS)) {
    if (!pheno_c) {
      logprint("Error: --filter-cases/--filter-controls requires a case/control phenotype.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    ii = indiv_exclude_ct;
    // fcc == 1: exclude all zeroes in pheno_c
    // fcc == 2: exclude all ones in pheno_c
    // -> flip on fcc == 1
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, pheno_c, (filter_binary & FILTER_BINARY_CASES)? 1 : 0, pheno_nm);
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to case/control status (--filter-%s).\n", ii, species_str(ii), (filter_binary & FILTER_BINARY_CASES)? "cases" : "controls");
    logprintb();
  }
  if (filter_binary & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_male, (filter_binary & FILTER_BINARY_MALES)? 1 : 0, sex_nm);
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to gender filter (--filter-%s).\n", ii, species_str(ii), (filter_binary & FILTER_BINARY_MALES)? "males" : "females");
    logprintb();
  }
  if (filter_binary & (FILTER_BINARY_FOUNDERS | FILTER_BINARY_NONFOUNDERS)) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, founder_info, (filter_binary & FILTER_BINARY_FOUNDERS)? 1 : 0, NULL);
    ii = indiv_exclude_ct - ii;
    sprintf(logbuf, "%d %s removed due to founder status (--filter-%s).\n", ii, species_str(ii), (filter_binary & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
    logprintb();
  }

  if (fseeko(bedfile, 0, SEEK_END)) {
    goto wdist_ret_READ_FAIL;
  }
  llxx = ftello(bedfile);
  llyy = llxx - ((uint64_t)unfiltered_indiv_ct4) * unfiltered_marker_ct;
  rewind(bedfile);
  if (llyy == 3LL) {
    // v1.00 or later
    if (fread(tbuf, 1, 3, bedfile) < 3) {
      goto wdist_ret_READ_FAIL;
    }
    if (memcmp(tbuf, "l\x1b\x01", 3)) {
      if (memcmp(tbuf, "l\x1b", 3)) {
	if ((tbuf[0] == '#') || (!memcmp(tbuf, "chr", 3))) {
          logprint("Error: Invalid header bytes in PLINK .bed file.  (Is this a UCSC Genome Browser\nBED file instead?)\n");
	} else {
	  logprint("Error: Invalid header bytes in PLINK .bed file.\n");
	}
	goto wdist_ret_INVALID_FORMAT;
      }
      bed_offset = 2;
    }
  } else if (llyy == 1LL) {
    // v0.99
    if (fread(tbuf, 1, 1, bedfile) != 1) {
      goto wdist_ret_READ_FAIL;
    }
    if (*tbuf == '\x01') {
      bed_offset = 1;
    } else if (*tbuf == '\0') {
      bed_offset = 2;
    } else {
      logprint("Error: Invalid header bytes in .bed file.\n");
      goto wdist_ret_INVALID_FORMAT;
    }
  } else if (llyy != 0LL) {
    sprintf(logbuf, "Error: Invalid .bed file size (expected %" PRIu64 " bytes).\n", 3LLU + ((uint64_t)unfiltered_indiv_ct4) * unfiltered_marker_ct);
    goto wdist_ret_INVALID_FORMAT_2;
  } else {
    // pre-0.99, no magic number, indiv-major
    bed_offset = 2;
  }
  if (bed_offset == 2) {
    strcpy(outname_end, ".bed.tmp"); // not really temporary
    logprint("Individual-major .bed file detected.  Transposing to SNP-major form.\n");
    fclose(bedfile);
    retval = indiv_major_to_snp_major(pedname, outname, unfiltered_marker_ct);
    if (retval) {
      goto wdist_ret_1;
    }
    strcpy(pedname, outname);
    if (fopen_checked(&bedfile, pedname, "rb")) {
      goto wdist_ret_OPEN_FAIL;
    }
    bed_offset = 3;
  }
  if (mind_thresh < 1.0) {
    retval = mind_filter(bedfile, bed_offset, mind_thresh, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct);
    if (retval) {
      goto wdist_ret_1;
    }
  }
  if (cluster_ptr->fname) {
    retval = load_clusters(cluster_ptr->fname, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, cluster_ptr->keep_fname, cluster_ptr->keep_flattened, cluster_ptr->remove_fname, cluster_ptr->remove_flattened);
    if (retval) {
      goto wdist_ret_1;
    }
  }
  // er, obvious todo: convert this to a local variable
  g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  if (!g_indiv_ct) {
    sprintf(logbuf, "Error: No %s pass QC.\n", g_species_plural);
    goto wdist_ret_INVALID_FORMAT_2;
  }
  if ((g_indiv_ct == 1) && (relationship_or_ibc_req(calculation_type) || distance_req(calculation_type, read_dists_fname) || (calculation_type & (CALC_GENOME | CALC_CLUSTER | CALC_NEIGHBOR)))) {
    sprintf(logbuf, "Error: More than 1 %s required for pairwise analysis.\n", g_species_singular);
    goto wdist_ret_INVALID_FORMAT_2;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  bitfield_andnot(pheno_nm, indiv_exclude, unfiltered_indiv_ctl);
  if (pheno_c) {
    bitfield_and(pheno_c, pheno_nm, unfiltered_indiv_ctl);
  }
  bitfield_andnot(founder_info, indiv_exclude, unfiltered_indiv_ctl);
  pheno_nm_ct = popcount_longs(pheno_nm, 0, unfiltered_indiv_ctl);
  if (!pheno_nm_ct) {
    logprint("Note: No phenotypes present.\n");
    misc_flags |= MISC_HWE_ALL;
  } else if (pheno_c) {
    pheno_ctrl_ct = popcount_longs_exclude(pheno_nm, pheno_c, unfiltered_indiv_ctl);
    if (pheno_nm_ct != g_indiv_ct) {
      sprintf(logbuf, "%u case%s, %u control%s, and %" PRIuPTR " missing phenotype%s present.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "" : "s", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "" : "s", g_indiv_ct - pheno_nm_ct, (g_indiv_ct - pheno_nm_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "%u case%s and %u control%s present.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "" : "s", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "" : "s");
    }
    logprintb();
    if (!pheno_ctrl_ct) {
      misc_flags |= MISC_HWE_ALL;
    }
  } else {
    if (pheno_nm_ct != g_indiv_ct) {
      sprintf(logbuf, "%u quantitative phenotype%s present (%" PRIuPTR " missing).\n", pheno_nm_ct, (pheno_nm_ct == 1)? "" : "s", g_indiv_ct - pheno_nm_ct);
    } else {
      sprintf(logbuf, "%u quantitative phenotype%s present.\n", pheno_nm_ct, (pheno_nm_ct == 1)? "" : "s");
    }
    logprintb();
    misc_flags |= MISC_HWE_ALL;
  }

  if ((parallel_tot > 1) && (parallel_tot > g_indiv_ct / 2)) {
    sprintf(logbuf, "Error: Too many --parallel jobs (maximum %" PRIuPTR "/2 = %" PRIuPTR ").\n", g_indiv_ct, g_indiv_ct / 2);
    goto wdist_ret_INVALID_CMDLINE_2;
  }
  if (g_thread_ct > 1) {
    if ((calculation_type & (CALC_RELATIONSHIP | CALC_IBC | CALC_GDISTANCE_MASK | CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE | CALC_GENOME | CALC_REGRESS_REL | CALC_UNRELATED_HERITABILITY | CALC_LASSO)) || ((calculation_type & CALC_MODEL) && (model_modifier & (MODEL_PERM | MODEL_MPERM))) || ((calculation_type & CALC_GLM) && (glm_modifier & (GLM_PERM | GLM_MPERM))) || ((calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) && (!read_genome_fname) && ((cluster_ptr->ppc != 0.0) || (!read_dists_fname)))) {
      sprintf(logbuf, "Using %d threads (change this with --threads).\n", g_thread_ct);
      logprintb();
    } else {
      logprint("Using 1 thread (no multithreaded calculations invoked).\n");
    }
  }

  if ((calculation_type & (CALC_MISSING_REPORT | CALC_GENOME | CALC_HOMOZYG)) || cluster_ptr->mds_dim_ct) {
    calc_plink_maxfid(unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, &plink_maxfid, &plink_maxiid);
  }
  plink_maxsnp = calc_plink_maxsnp(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len);

  if (indiv_sort & (INDIV_SORT_NATURAL | INDIV_SORT_ASCII)) {
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, (indiv_sort & INDIV_SORT_NATURAL)? strcmp_natural_deref : strcmp_deref);
    if (retval) {
      goto wdist_ret_1;
    }
    indiv_sort_map = uiptr;
    wkspace_reset((unsigned char*)cptr);
  }

  if ((misc_flags & (MISC_MAKE_FOUNDERS | MISC_MAKE_FOUNDERS_FIRST)) == MISC_MAKE_FOUNDERS) {
    if (make_founders(unfiltered_indiv_ct, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, indiv_exclude, founder_info)) {
      goto wdist_ret_NOMEM;
    }
  }
  if ((calculation_type & CALC_GENOME) || genome_skip_write) {
    retval = populate_pedigree_rel_info(&pri, unfiltered_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & CALC_WRITE_CLUSTER) {
    retval = write_clusters(outname, outname_end, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, (misc_flags / MISC_WRITE_CLUSTER_OMIT_UNASSIGNED) & 1, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER))))) {
      goto wdist_ret_1;
    }
  }
  if (covar_fname) {
    // update this as more covariate-referencing commands are added
    if (!(calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_WRITE_COVAR | CALC_GXE | CALC_GLM | CALC_LASSO))) {
      logprint("Note: Ignoring --covar since no commands reference the covariates.\n");
    } else {
      // if only --gxe, ignore --covar-name/--covar-number
      uii = (calculation_type & (CALC_MAKE_BED | CALC_RECODE | CALC_WRITE_COVAR | CALC_GLM | CALC_LASSO))? 1 : 0;
      retval = load_covars(covar_fname, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, missing_phenod, uii? covar_modifier : (covar_modifier & COVAR_KEEP_PHENO_ON_MISSING_COV), uii? covar_range_list_ptr : NULL, gxe_mcovar, &covar_ct, &covar_names, &max_covar_name_len, pheno_nm, &pheno_nm_ct, &covar_nm, &covar_d, &gxe_covar_nm, &gxe_covar_c);
      if (retval) {
	goto wdist_ret_1;
      }
    }
  }

  retval = calc_freqs_and_hwe(bedfile, outname, outname_end, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, founder_info, nonfounders, (misc_flags / MISC_MAF_SUCC) & 1, set_allele_freqs, &marker_reverse, &marker_allele_cts, bed_offset, (hwe_thresh > 0.0) || (calculation_type & CALC_HARDY), (misc_flags / MISC_HWE_ALL) & 1, (pheno_nm_ct && pheno_c)? (calculation_type & CALC_HARDY) : 0, pheno_nm, pheno_nm_ct? pheno_c : NULL, &hwe_lls, &hwe_lhs, &hwe_hhs, &hwe_ll_cases, &hwe_lh_cases, &hwe_hh_cases, &hwe_ll_allfs, &hwe_lh_allfs, &hwe_hh_allfs, &hwe_hapl_allfs, &hwe_haph_allfs, &indiv_male_ct, &indiv_f_ct, &indiv_f_male_ct, wt_needed, &marker_weights_base, &marker_weights, exponent, chrom_info_ptr, sex_nm, sex_male, map_is_unsorted & UNSORTED_SPLIT_CHROM, &hh_exists);
  if (retval) {
    goto wdist_ret_1;
  }

  if (freqname) {
    retval = read_external_freqs(freqname, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, chrom_info_ptr, marker_allele_ptrs, marker_allele_cts, set_allele_freqs, (misc_flags / MISC_MAF_SUCC) & 1, exponent, wt_needed, marker_weights);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (!(misc_flags & MISC_KEEP_ALLELE_ORDER)) {
    calc_marker_reverse_bin(marker_reverse, marker_exclude, unfiltered_marker_ct, unfiltered_marker_ct - marker_exclude_ct, set_allele_freqs);
  }
  if (a1alleles || a2alleles) {
    retval = load_ax_alleles(a1alleles? a1alleles : a2alleles, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_allele_ptrs, &max_marker_allele_len, marker_reverse, marker_ids, max_marker_id_len, a2alleles? 1 : 0);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (marker_allele_ptrs) {
    swap_reversed_marker_alleles(unfiltered_marker_ct, marker_reverse, marker_allele_ptrs);
  }

  // contrary to the PLINK flowchart, --freq effectively resolves before
  // --geno.
  if (calculation_type & CALC_FREQ) {
    if (cluster_ct && (!(misc_flags & MISC_FREQX))) {
      memcpy(outname_end, ".frq.strat", 11);
      retval = write_stratified_freqs(bedfile, bed_offset, outname, plink_maxsnp, unfiltered_marker_ct, marker_exclude, zero_extra_chroms, chrom_info_ptr, marker_ids, max_marker_id_len, marker_allele_ptrs, max_marker_allele_len, unfiltered_indiv_ct, g_indiv_ct, indiv_f_ct, founder_info, nonfounders, sex_male, indiv_f_male_ct, marker_reverse, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len);
    } else {
      if (misc_flags & MISC_FREQX) {
	memcpy(outname_end, ".frqx", 6);
      } else if (misc_flags & MISC_FREQ_COUNTS) {
	memcpy(outname_end, ".frq.count", 11);
      } else {
	memcpy(outname_end, ".frq", 5);
      }
      retval = write_freqs(outname, plink_maxsnp, unfiltered_marker_ct, marker_exclude, set_allele_freqs, zero_extra_chroms, chrom_info_ptr, marker_ids, max_marker_id_len, marker_allele_ptrs, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, hwe_hapl_allfs, hwe_haph_allfs, indiv_f_ct, indiv_f_male_ct, misc_flags, marker_reverse);
    }
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ))))) {
      goto wdist_ret_1;
    }
  }
  if (calculation_type & CALC_MISSING_REPORT) {
    retval = write_missingness_reports(bedfile, bed_offset, outname, outname_end, plink_maxfid, plink_maxiid, plink_maxsnp, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, zero_extra_chroms, chrom_info_ptr, marker_ids, max_marker_id_len, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, pheno_nm, sex_male, indiv_male_ct, person_ids, max_person_id_len, cluster_ct, cluster_map, cluster_starts, cluster_ids, max_cluster_id_len, hh_exists);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_MISSING_REPORT))))) {
      goto wdist_ret_1;
    }
  }
  if (geno_thresh < 1.0) {
    ulii = binary_geno_filter(geno_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, g_indiv_ct, indiv_male_ct, marker_allele_cts, chrom_info_ptr);
    sprintf(logbuf, "%" PRIuPTR " marker%s removed due to missing genotype data (--geno).\n", ulii, (ulii == 1)? "" : "s");
    logprintb();
  }
  wkspace_reset(marker_allele_cts);
  marker_allele_cts = NULL;
  if (calculation_type & CALC_HARDY) {
    retval = hardy_report(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, hwe_lls, hwe_lhs, hwe_hhs, hwe_ll_cases, hwe_lh_cases, hwe_hh_cases, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, pheno_nm_ct, pheno_c, zero_extra_chroms, chrom_info_ptr);
    if (retval || (!(calculation_type & (~(CALC_MERGE | CALC_WRITE_CLUSTER | CALC_FREQ | CALC_HARDY))))) {
      goto wdist_ret_1;
    }
  }
  if (hwe_thresh > 0.0) {
    enforce_hwe_threshold(hwe_thresh, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, hwe_lls, hwe_lhs, hwe_hhs, (misc_flags / MISC_HWE_ALL) & 1, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs);
  }
  if ((min_maf != 0.0) || (max_maf != 0.5)) {
    enforce_maf_threshold(min_maf, max_maf, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, set_allele_freqs);
  }
  if (min_bp_space) {
    if (map_is_unsorted & (UNSORTED_BP | UNSORTED_SPLIT_CHROM)) {
      logprint("Error: --bp-space requires a sorted .bim file.\n");
      goto wdist_ret_INVALID_FORMAT;
    }
    enforce_min_bp_space(min_bp_space, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, chrom_info_ptr);
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (!marker_ct) {
    logprint("Error: All markers fail QC.\n");
    goto wdist_ret_INVALID_FORMAT;
  }

  if (wt_needed) {
    calc_marker_weights(exponent, unfiltered_marker_ct, marker_exclude, marker_ct, hwe_ll_allfs, hwe_lh_allfs, hwe_hh_allfs, marker_weights);
  }
  wkspace_reset(hwe_lls);

  if (wt_needed) {
    // N.B. marker_weights now on top of stack
    wkspace_reset(marker_weights_base);
    // normalize included marker weights to add to just under 2^32.  (switch to
    // 2^64 if/when 32-bit performance becomes less important than accuracy on
    // 50+ million marker sets.)
    dxx = 0.0;
    marker_uidx = 0;
    marker_idx = 0;
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      marker_idx += marker_uidx_stop - marker_uidx;
      dptr = &(marker_weights[marker_uidx]);
      dptr2 = &(marker_weights[marker_uidx_stop]);
      marker_uidx = marker_uidx_stop;
      do {
        dxx += *dptr++;
      } while (dptr < dptr2);
    } while (marker_idx < marker_ct);
    // subtract marker_ct to guard against marker_weight_sum overflow from
    // rounding
    dxx = (4294967296.0 - ((double)((intptr_t)marker_ct))) / dxx;
    marker_weights_i = (uint32_t*)wkspace_alloc(marker_idx * sizeof(int32_t));
    marker_uidx = 0;
    marker_idx = 0;
    uiptr = marker_weights_i;
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      marker_uidx_stop = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
      marker_idx += marker_uidx_stop - marker_uidx;
      dptr = &(marker_weights[marker_uidx]);
      dptr2 = &(marker_weights[marker_uidx_stop]);
      marker_uidx = marker_uidx_stop;
      do {
        uii = (uint32_t)((*dptr++) * dxx + 0.5);
        marker_weight_sum += uii;
        *uiptr++ = uii;
      } while (dptr < dptr2);
    } while (marker_idx < marker_ct);
  }

  sprintf(logbuf, "%" PRIuPTR " marker%s and %" PRIuPTR " %s pass filters and QC%s.\n", marker_ct, (marker_ct == 1)? "" : "s", g_indiv_ct, species_str(g_indiv_ct), (calculation_type & CALC_REL_CUTOFF)? " (before --rel-cutoff)": "");
  logprintb();

  if (calculation_type & CALC_WRITE_SNPLIST) {
    retval = write_snplist(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, NULL, 0);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & CALC_LIST_23_INDELS) {
    retval = write_snplist(outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, 1);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & (CALC_WRITE_COVAR | CALC_MAKE_BED | CALC_RECODE)) {
    if (pheno_nosex_exclude) {
      bitfield_andnot_reversed_args(pheno_nosex_exclude, pheno_nm, unfiltered_indiv_ctl);
    }
    if (covar_fname) {
      retval = write_covars(outname, outname_end, write_covar_modifier, write_covar_dummy_max_categories, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nosex_exclude? pheno_nosex_exclude : pheno_nm, pheno_c, pheno_d, missing_phenod, output_missing_pheno, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (calculation_type & CALC_MAKE_BED) {
      retval = make_bed(bedfile, bed_offset, mapname, map_cols, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, marker_reverse, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nosex_exclude? pheno_nosex_exclude : pheno_nm, pheno_c, pheno_d, output_missing_pheno, map_is_unsorted, indiv_sort_map, misc_flags, update_chr, flip_subset_fname, hh_exists, chrom_info_ptr);
      if (retval) {
        goto wdist_ret_1;
      }
    }
    if (calculation_type & CALC_RECODE) {
      retval = recode(recode_modifier, bedfile, bed_offset, famfile, outname, outname_end, recode_allele_name, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, marker_ids, max_marker_id_len, marker_cms, marker_allele_ptrs, max_marker_allele_len, marker_pos, marker_reverse, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nosex_exclude? pheno_nosex_exclude : pheno_nm, pheno_c, pheno_d, output_missing_pheno, misc_flags, hh_exists, chrom_info_ptr);
      if (retval) {
        goto wdist_ret_1;
      }
    }
    if (pheno_nosex_exclude) {
      free(pheno_nosex_exclude);
      pheno_nosex_exclude = NULL;
    }
  }

  if (calculation_type & CALC_HOMOZYG) {
    if (map_is_unsorted & UNSORTED_BP) {
      logprint("Error: Run-of-homozygosity scanning requires a sorted .map/.bim.  Retry this\ncommand after using --make-bed to sort your data.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    retval = calc_homozyg(homozyg_ptr, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, zero_extra_chroms, chrom_info_ptr, marker_pos, g_indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, outname, outname_end, pheno_nm, pheno_c, pheno_d, missing_pheno, sex_male);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & CALC_LD_PRUNE) {
    if (map_is_unsorted & UNSORTED_BP) {
      logprint("Error: LD-based marker pruning requires a sorted .map/.bim.  Retry this command\nafter using --make-bed to sort your data.\n");
      goto wdist_ret_INVALID_CMDLINE;
    }
    retval = ld_prune(bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, chrom_info_ptr, set_allele_freqs, marker_pos, unfiltered_indiv_ct, founder_info, sex_male, ld_window_size, ld_window_kb, ld_window_incr, ld_last_param, outname, outname_end, misc_flags, hh_exists);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & CALC_REGRESS_PCS) {
    retval = calc_regress_pcs(evecname, regress_pcs_modifier, max_pcs, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, marker_reverse, marker_ids, max_marker_id_len, marker_allele_ptrs, zero_extra_chroms, chrom_info_ptr, marker_pos, g_indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, sex_nm, sex_male, pheno_d, missing_phenod, outname, outname_end, hh_exists);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  // possibly no more need for marker_ids/marker_alleles, conditional unload to
  // clear space for IBS matrix, etc.?  (probably want to initially load at far
  // end of stack to make this workable...)

  if (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) {
    wkspace_mark_postcluster = wkspace_base;
    ulii = (g_indiv_ct * (g_indiv_ct - 1)) >> 1;
    if (cluster_ptr->mds_dim_ct) {
#ifndef __LP64__
      // catch 32-bit intptr_t overflow
      if (g_indiv_ct > 23169) {
        goto wdist_ret_NOMEM;
      }
#endif
      if (((!read_dists_fname) && (!read_genome_fname)) || (cluster_ptr->modifier & CLUSTER_MISSING)) {
	if ((!(cluster_ptr->modifier & CLUSTER_MDS)) || (!cluster_ct)) {
          if (wkspace_alloc_d_checked(&mds_plot_dmatrix_copy, ulii * sizeof(double))) {
            goto wdist_ret_NOMEM;
          }
	} else {
	  ulii = cluster_ct + g_indiv_ct - cluster_starts[cluster_ct];
          if (wkspace_alloc_d_checked(&mds_plot_dmatrix_copy, (ulii * (ulii - 1)) * (sizeof(double) / 2))) {
            goto wdist_ret_NOMEM;
          }
	}
      }
    }

    if (cluster_ct) {
      ulii = cluster_ct + g_indiv_ct - cluster_starts[cluster_ct];
#ifndef __LP64__
      if (ulii > 23169) {
	goto wdist_ret_NOMEM;
      }
#endif
      ulii = (ulii * (ulii - 1)) >> 1;
#ifndef __LP64__
    } else if (g_indiv_ct > 23169) {
      goto wdist_ret_NOMEM;
#endif
    }
    if (wkspace_alloc_ul_checked(&cluster_merge_prevented, ((ulii + (BITCT - 1)) / BITCT) * sizeof(intptr_t))) {
      goto wdist_ret_NOMEM;
    }
    if (cluster_ct || (calculation_type & CALC_GENOME) || genome_skip_write) {
      if (wkspace_alloc_d_checked(&cluster_sorted_ibs, ulii * sizeof(double))) {
	goto wdist_ret_NOMEM;
      }
      if (cluster_ptr->modifier & CLUSTER_GROUP_AVG) {
        fill_double_zero(cluster_sorted_ibs, ulii);
      } else {
	for (uljj = 0; uljj < ulii; uljj++) {
	  cluster_sorted_ibs[uljj] = 1.0;
	}
      }
    }
    wkspace_mark_precluster = wkspace_base;
  }

  wkspace_mark2 = wkspace_base;

  if (relationship_or_ibc_req(calculation_type)) {
    if (rel_calc_type & REL_CALC_SINGLE_PREC) {
      retval = calc_rel_f(threads, parallel_idx, parallel_tot, calculation_type, rel_calc_type, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, ibc_type, (float)rel_cutoff, set_allele_freqs, chrom_info_ptr);
    } else {
      retval = calc_rel(threads, parallel_idx, parallel_tot, calculation_type, rel_calc_type, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, ibc_type, rel_cutoff, set_allele_freqs, &rel_ibc, chrom_info_ptr);
    }
    if (retval) {
      goto wdist_ret_1;
    }

    if (calculation_type & CALC_REGRESS_REL) {
      retval = regress_rel_main(unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, regress_rel_iters, regress_rel_d, threads, pheno_d);
      if (retval) {
	goto wdist_ret_1;
      }
    }
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      retval = calc_unrelated_herit(calculation_type, ibc_type, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, pheno_d, rel_ibc, (misc_flags / MISC_UNRELATED_HERITABILITY_STRICT) & 1, unrelated_herit_covg, unrelated_herit_covr, unrelated_herit_tol);
      if (retval) {
	goto wdist_ret_1;
      }
    }
#endif
    wkspace_reset(g_indiv_missing_unwt);
    g_indiv_missing_unwt = NULL;
    g_missing_dbl_excluded = NULL;
  }

  if (calculation_type & CALC_REGRESS_PCS_DISTANCE) {
    logprint("Error: --regress-pcs-distance has not yet been written.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto wdist_ret_1;
  } else if (distance_req(calculation_type, read_dists_fname)) {
    retval = calc_distance(threads, parallel_idx, parallel_tot, bedfile, bed_offset, outname, outname_end, calculation_type, dist_calc_type, marker_exclude, marker_ct, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len, chrom_info_ptr, wt_needed, marker_weight_sum, marker_weights_i, exponent);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (read_dists_fname && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))) {
    // use delayed and specialized load for --cluster/--neighbour, since PPC
    // values may be needed, and user may want to process a distance matrix too
    // large to be loaded in memory by doing some pre-clustering
    dists_alloc = (g_indiv_ct * (g_indiv_ct - 1)) * (sizeof(double) / 2);
    if (wkspace_alloc_d_checked(&g_dists, dists_alloc)) {
      goto wdist_ret_NOMEM;
    }
    retval = read_dists(read_dists_fname, read_dists_id_fname, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, person_ids, max_person_id_len, 0, NULL, NULL, 0, 0, g_dists, 0, NULL, NULL);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & CALC_IBS_TEST) {
    retval = ibs_test_calc(threads, read_dists_fname, unfiltered_indiv_ct, indiv_exclude, ibs_test_perms, pheno_nm_ct, pheno_ctrl_ct, pheno_nm, pheno_c);
    if (retval) {
      goto wdist_ret_1;
    }
  }
  if (calculation_type & CALC_GROUPDIST) {
    retval = groupdist_calc(threads, unfiltered_indiv_ct, indiv_exclude, groupdist_iters, groupdist_d, pheno_nm_ct, pheno_ctrl_ct, pheno_nm, pheno_c);
    if (retval) {
      goto wdist_ret_1;
    }
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    retval = regress_distance(calculation_type, g_dists, pheno_d, unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, g_thread_ct, regress_iters, regress_d);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (read_dists_fname && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))) {
    wkspace_reset((unsigned char*)g_dists);
    g_dists = NULL;
  }

  if ((calculation_type & CALC_GENOME) || genome_skip_write) {
    wkspace_reset(wkspace_mark2);
    g_dists = NULL;
    retval = calc_genome(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, marker_pos, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info, parallel_idx, parallel_tot, outname, outname_end, nonfounders, calculation_type, genome_modifier, ppc_gap, genome_min_pi_hat, genome_max_pi_hat, pheno_nm, pheno_c, pri, genome_skip_write);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) {
    retval = calc_cluster_neighbor(threads, bedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, chrom_info_ptr, set_allele_freqs, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, read_dists_fname, read_dists_id_fname, read_genome_fname, outname, outname_end, calculation_type, cluster_ct, cluster_map, cluster_starts, cluster_ptr, missing_pheno, neighbor_n1, neighbor_n2, ppc_gap, pheno_c, mds_plot_dmatrix_copy, cluster_merge_prevented, cluster_sorted_ibs, wkspace_mark_precluster, wkspace_mark_postcluster);
    if (retval) {
      goto wdist_ret_1;
    }
  }

  if (calculation_type & (CALC_MODEL | CALC_GXE | CALC_GLM | CALC_LASSO | CALC_CMH | CALC_HOMOG)) {
    if ((!pheno_all) && (!loop_assoc_fname)) {
      outname_end2 = outname_end;
      goto wdist_skip_all_pheno;
    }
    uii = 0; // phenotype/cluster number
    if (loop_assoc_fname) {
      retval = load_clusters(loop_assoc_fname, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, NULL, NULL, NULL, NULL);
      if (retval) {
	goto wdist_ret_1;
      }
      *outname_end = '.';
      if (pheno_d) {
	free(pheno_d);
	pheno_d = NULL;
      }
      if (!pheno_c) {
	pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      }
    } else {
      wkspace_mark = wkspace_base;
      retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
      if (retval) {
	goto wdist_ret_1;
      }
      memcpy(outname_end, ".P", 2);
    }
    do {
      if (loop_assoc_fname) {
	if (uii == cluster_ct) {
	  break;
	}
	outname_end2 = strcpya(&(outname_end[1]), &(cluster_ids[uii * max_cluster_id_len]));
	fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
	ukk = cluster_starts[uii + 1];
	for (ujj = cluster_starts[uii]; ujj < ukk; ujj++) {
	  SET_BIT(pheno_c, cluster_map[ujj]);
	}
	uii++;
      } else {
	// --all-pheno
	if (pheno_modifier & PHENO_MERGE) {
	  memcpy(pheno_nm, orig_pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
	  if (orig_pheno_c) {
	    if (!pheno_c) {
	      free(pheno_d);
	      pheno_d = NULL;
	      pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
	    }
	    memcpy(pheno_c, orig_pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t));
	  } else {
	    memcpy(pheno_d, orig_pheno_d, unfiltered_indiv_ct * sizeof(double));
	  }
	} else {
	  fill_ulong_zero(pheno_nm, unfiltered_indiv_ctl);
	  if (pheno_c) {
	    free(pheno_c);
	    pheno_c = NULL;
	  }
	  if (pheno_d) {
	    free(pheno_d);
	    pheno_d = NULL;
	  }
	}
	uii++;
      wdist_skip_empty_pheno:
	rewind(phenofile);
	retval = load_pheno(phenofile, unfiltered_indiv_ct, indiv_exclude_ct, cptr, max_person_id_len, uiptr, missing_pheno, missing_pheno_len, (misc_flags / MISC_AFFECTION_01) & 1, uii, NULL, pheno_nm, &pheno_c, &pheno_d);
	if (retval == LOAD_PHENO_LAST_COL) {
	  wkspace_reset(wkspace_mark);
	  break;
	} else if (retval) {
	  goto wdist_ret_1;
	}
	bitfield_andnot(pheno_nm, indiv_exclude, unfiltered_indiv_ctl);
	if (gender_unk_ct && (!(sex_missing_pheno & ALLOW_NO_SEX))) {
	  bitfield_and(pheno_nm, sex_nm, unfiltered_indiv_ctl);
	}
	pheno_nm_ct = popcount_longs(pheno_nm, 0, unfiltered_indiv_ctl);
	if (!pheno_nm_ct) {
	  goto wdist_skip_empty_pheno;
	}
	outname_end2 = uint32_write(&(outname_end[2]), uii);
      }
      *outname_end2 = '\0';
    wdist_skip_all_pheno:
      if (calculation_type & CALC_MODEL) {
	if (pheno_d) {
	  retval = qassoc(threads, bedfile, bed_offset, outname, outname_end2, model_modifier, model_mperm_val, pfilter, mtest_adjust, adjust_lambda, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_allele_ptrs, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, cluster_ct, cluster_map, cluster_starts, aperm_min, aperm_max, aperm_alpha, aperm_beta, aperm_init_interval, aperm_interval_slope, mperm_save, pheno_nm_ct, pheno_nm, pheno_d, sex_male, hh_exists, perm_batch_size);
	} else {
	  retval = model_assoc(threads, bedfile, bed_offset, outname, outname_end2, model_modifier, model_cell_ct, model_mperm_val, ci_size, ci_zt, pfilter, mtest_adjust, adjust_lambda, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_len, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, cluster_ct, cluster_map, loop_assoc_fname? NULL : cluster_starts, aperm_min, aperm_max, aperm_alpha, aperm_beta, aperm_init_interval, aperm_interval_slope, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, sex_male);
	}
	if (retval) {
	  goto wdist_ret_1;
	}
      }
      if (calculation_type & CALC_GLM) {
	if (!(glm_modifier & GLM_NO_SNP)) {
          retval = glm_assoc(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, mtest_adjust, adjust_lambda, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, marker_allele_ptrs, max_marker_allele_len, marker_reverse, zero_extra_chroms, condition_mname, condition_fname, chrom_info_ptr, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, cluster_ct, cluster_map, cluster_starts, aperm_min, aperm_max, aperm_alpha, aperm_beta, aperm_init_interval, aperm_interval_slope, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, perm_batch_size);
	} else {
	  retval = glm_assoc_nosnp(threads, bedfile, bed_offset, outname, outname_end2, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, parameters_range_list_ptr, tests_range_list_ptr, ci_size, ci_zt, pfilter, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, cluster_ct, cluster_map, cluster_starts, mperm_save, pheno_nm_ct, pheno_nm, pheno_c, pheno_d, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, perm_batch_size);
	}
	if (retval) {
	  goto wdist_ret_1;
	}
      }
      // if case/control phenotype loaded with --all-pheno, skip --gxe
      if ((calculation_type & CALC_GXE) && pheno_d) {
	retval = gxe_assoc(bedfile, bed_offset, outname, outname_end, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, g_indiv_ct, indiv_exclude, pheno_nm, pheno_d, gxe_covar_nm, gxe_covar_c, sex_male, hh_exists);
	if (retval) {
	  goto wdist_ret_1;
	}
      }
      if ((calculation_type & CALC_LASSO) && pheno_d) {
	retval = lasso(threads, bedfile, bed_offset, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, marker_reverse, zero_extra_chroms, chrom_info_ptr, unfiltered_indiv_ct, pheno_nm_ct, lasso_h2, lasso_minlambda, (misc_flags / MISC_LASSO_REPORT_ZEROES) & 1, pheno_nm, pheno_d, covar_ct, covar_nm, covar_d, sex_male, hh_exists);
        if (retval) {
	  goto wdist_ret_1;
	}
      }
      if ((calculation_type & CALC_CMH) && pheno_c) {
        retval = assoc_cmh();
        if (retval) {
          goto wdist_ret_1;
	}
      }
      if ((calculation_type & CALC_HOMOG) && pheno_c) {
	retval = assoc_homog();
        if (retval) {
          goto wdist_ret_1;
	}
      }
    } while (pheno_all || loop_assoc_fname);
  }
  while (0) {
  wdist_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  wdist_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  wdist_ret_INVALID_FORMAT_2:
    logprintb();
  wdist_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  wdist_ret_INVALID_CMDLINE_2:
    logprintb();
  wdist_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  wdist_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
 wdist_ret_1:
  free_cond(pheno_nosex_exclude);
  free_cond(orig_pheno_d);
  free_cond(orig_pheno_c);
  free_cond(orig_pheno_nm);
  free_cond(pheno_d);
  free_cond(pheno_c);
  free_cond(id_list);
  free_cond(pid_list);
  fclose_cond(phenofile);
  fclose_cond(famfile);
  fclose_cond(bedfile);
  if (marker_allele_ptrs && (max_marker_allele_len > 2)) {
    ulii = unfiltered_marker_ct * 2;
    for (marker_uidx = 0; marker_uidx < ulii; marker_uidx++) {
      cptr = marker_allele_ptrs[marker_uidx];
      if ((cptr < g_one_char_strs) || (cptr >= &(g_one_char_strs[512]))) {
	free(cptr);
      }
    }
  }
  return retval;
}

// output-missing-phenotype + terminating null
#define MAX_FLAG_LEN 25

static inline int32_t is_flag(char* param) {
  char cc = param[1];
  return ((*param == '-') && ((cc > '9') || (cc < '0'))); 
}

static inline char* is_flag_start(char* param) {
  char cc = param[1];
  if ((*param == '-') && ((cc > '9') || (cc < '0'))) {
    return (cc == '-')? (&(param[2])) : (&(param[1]));
  }
  return NULL;
}

uint32_t param_count(int32_t argc, char** argv, int32_t flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any parameter not beginning with "--" as optional.
  int32_t opt_params = 0;
  int32_t cur_idx = flag_idx + 1;
  while (cur_idx < argc) {
    if (is_flag(argv[cur_idx])) {
      break;
    }
    opt_params++;
    cur_idx++;
  }
  return opt_params;
}

int32_t enforce_param_ct_range(uint32_t param_ct, char* flag_name, uint32_t min_ct, uint32_t max_ct) {
  if (param_ct > max_ct) {
    if (max_ct > min_ct) {
      sprintf(logbuf, "Error: %s accepts at most %d parameter%s.%s", flag_name, max_ct, (max_ct == 1)? "" : "s", errstr_append);
    } else {
      sprintf(logbuf, "Error: %s only accepts %d parameter%s.%s", flag_name, max_ct, (max_ct == 1)? "" : "s", errstr_append);
    }
    return -1;
  } else if (param_ct < min_ct) {
    if (min_ct == 1) {
      sprintf(logbuf, "Error: Missing %s parameter.%s", flag_name, errstr_append);
    } else {
      sprintf(logbuf, "Error: %s requires %s%d parameters.%s", flag_name, (min_ct < max_ct)? "at least " : "", min_ct, errstr_append);
    }
    return -1;
  }
  return 0;
}

int32_t parse_next_range(uint32_t param_ct, char range_delim, char** argv, uint32_t* cur_param_idx_ptr, char** cur_arg_pptr, char** range_start_ptr, uint32_t* rs_len_ptr, char** range_end_ptr, uint32_t* re_len_ptr) {
  // Starts reading from argv[cur_param_idx][cur_pos].  If a valid range is
  // next, range_start + rs_len + range_end + re_len are updated.  If only a
  // single item is next, range_end is set to NULL and range_start + rs_len are
  // updated.  If there are no items left, range_start is set to NULL.  If
  // the input is not well-formed, -1 is returned instead of 0.
  uint32_t cur_param_idx = *cur_param_idx_ptr;
  char* cur_arg_ptr = *cur_arg_pptr;
  char cc;
  if (cur_param_idx > param_ct) {
    *cur_arg_pptr = NULL;
    return 0;
  }
  while (1) {
    cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
	*range_start_ptr = NULL;
	return 0;
      }
      cur_arg_ptr = argv[cur_param_idx];
      cc = *cur_arg_ptr;
    }
    if (cc == range_delim) {
      return -1;
    }
    if (cc != ',') {
      break;
    }
    cur_arg_ptr++;
  }
  *range_start_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if ((!cc) || (cc == ',')) {
      *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
      *cur_arg_pptr = cur_arg_ptr;
      *range_end_ptr = NULL;
      return 0;
    }
  } while (cc != range_delim);
  *rs_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_start_ptr));
  cc = *(++cur_arg_ptr);
  if ((!cc) || (cc == ',') || (cc == range_delim)) {
    return -1;
  }
  *range_end_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if (cc == range_delim) {
      return -1;
    }
  } while (cc && (cc != ','));
  *re_len_ptr = (uintptr_t)(cur_arg_ptr - (*range_end_ptr));
  *cur_arg_pptr = cur_arg_ptr;
  return 0;
}

int32_t parse_chrom_ranges(uint32_t param_ct, char range_delim, char** argv, uintptr_t* chrom_mask, Chrom_info* chrom_info_ptr, uint32_t allow_extra_chroms, char* cur_flag_str) {
  uint32_t argct = 0;
  uint32_t cur_param_idx = 1;
  int32_t retval = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  int32_t marker_code_start;
  int32_t marker_code_end;
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", cur_flag_str, argv[cur_param_idx], errstr_append);
	goto parse_chrom_ranges_ret_INVALID_CMDLINE;
      }
      if (!range_start) {
	break;
      }
      marker_code_start = marker_code2(chrom_info_ptr, range_start, rs_len);
      if (marker_code_start == -1) {
	range_start[rs_len] = '\0';
	if (!allow_extra_chroms) {
	  sprintf(logbuf, "Error: Invalid --%s chromosome code '%s'.%s", cur_flag_str, range_start, errstr_append);
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE;
	} else if (range_end) {
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE_2;
	}
        if (push_ll_str(&(chrom_info_ptr->incl_excl_name_stack), range_start)) {
	  goto parse_chrom_ranges_ret_NOMEM;
	}
      } else if (range_end) {
        marker_code_end = marker_code2(chrom_info_ptr, range_end, re_len);
	if (marker_code_end == -1) {
	  if (!allow_extra_chroms) {
	    range_end[re_len] = '\0';
	    sprintf(logbuf, "Error: Invalid --%s chromosome code '%s'.%s", cur_flag_str, range_end, errstr_append);
	    goto parse_chrom_ranges_ret_INVALID_CMDLINE;
	  } else {
	    goto parse_chrom_ranges_ret_INVALID_CMDLINE_2;
	  }
	}
        if (marker_code_end <= marker_code_start) {
	  range_start[rs_len] = '\0';
	  range_end[re_len] = '\0';
	  sprintf(logbuf, "Error: --%s chromosome code '%s' is not greater than '%s'.%s", cur_flag_str, range_end, range_start, errstr_append);
	  goto parse_chrom_ranges_ret_INVALID_CMDLINE;
	}
	fill_bits(chrom_mask, marker_code_start, marker_code_end + 1 - marker_code_start);
      } else {
        SET_BIT(chrom_mask, marker_code_start);
      }
      argct++;
    }
  }
  if (!argct) {
    sprintf(logbuf, "Error: --%s requires at least one value.%s", cur_flag_str, errstr_append);
    logprintb();
    return -1;
  }
  while (0) {
  parse_chrom_ranges_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  parse_chrom_ranges_ret_INVALID_CMDLINE_2:
    logprint("Error: Chromosome ranges cannot include nonstandard names.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  parse_chrom_ranges_ret_INVALID_CMDLINE:
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

int32_t parse_name_ranges(uint32_t param_ct, char range_delim, char** argv, Range_list* range_list_ptr, uint32_t require_posint) {
  uint32_t name_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t name_max_len = 0;
  char* cur_arg_ptr;
  char* range_start;
  uint32_t rs_len;
  char* range_end;
  uint32_t re_len;
  char* cur_name_str;
  char* dup_check;
  unsigned char* cur_name_starts_range;
  int32_t last_val;
  int32_t cur_val;
  // two passes.  first pass: count parameters, determine name_max_len;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argv[1];
    while (1) {
      if (parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
	sprintf(logbuf, "Error: Invalid %s parameter '%s'.%s", argv[0], argv[cur_param_idx], errstr_append);
        logprintb();
        return RET_INVALID_CMDLINE;
      }
      if (!range_start) {
	break;
      }
      name_ct++;
      if (rs_len > name_max_len) {
	name_max_len = rs_len; // does NOT include trailing null yet
      }
      if (range_end) {
	name_ct++;
	if (re_len > name_max_len) {
	  name_max_len = re_len;
	}
      }
    }
  }
  if (!name_ct) {
    sprintf(logbuf, "Error: %s requires at least one value.%s", argv[0], errstr_append);
    logprintb();
    return RET_INVALID_CMDLINE;
  }
  range_list_ptr->name_max_len = ++name_max_len;
  range_list_ptr->name_ct = name_ct;
  range_list_ptr->names = (char*)malloc(name_ct * name_max_len * sizeof(char));
  if (!range_list_ptr->names) {
    return RET_NOMEM;
  }
  range_list_ptr->starts_range = (unsigned char*)malloc(name_ct * sizeof(char));
  if (!range_list_ptr->starts_range) {
    return RET_NOMEM;
  }
  cur_name_str = range_list_ptr->names;
  cur_name_starts_range = range_list_ptr->starts_range;
  cur_param_idx = 1;
  cur_arg_ptr = argv[1];
  while (1) {
    parse_next_range(param_ct, range_delim, argv, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      if (require_posint) {
	last_val = 0;
	for (cur_param_idx = 0; cur_param_idx < name_ct; cur_param_idx++) {
	  cur_name_str = &(range_list_ptr->names[cur_param_idx * name_max_len]);
	  dup_check = cur_name_str; // actually a numeric check
	  do {
	    if (is_not_digit(*dup_check)) {
	      sprintf(logbuf, "Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	      logprintb();
	      return RET_INVALID_CMDLINE;
	    }
	  } while (*(++dup_check));
	  cur_val = atoi(cur_name_str);
	  if (cur_val < 1) {
	    sprintf(logbuf, "Error: Invalid %s parameter '%s'.\n", argv[0], cur_name_str);
	    logprintb();
	    return RET_INVALID_CMDLINE;
	  }
	  if (range_list_ptr->starts_range[cur_param_idx]) {
	    last_val = cur_val;
	  } else {
	    if (cur_val <= last_val) {
	      sprintf(logbuf, "Error: Invalid %s range '%s-%s'.\n", argv[0], &(range_list_ptr->names[(cur_param_idx - 1) * name_max_len]), cur_name_str);
	      logprintb();
	      return RET_INVALID_CMDLINE;
	    }
	    last_val = 0;
	  }
	}
      }
      return 0;
    }
    memcpyx(cur_name_str, range_start, rs_len, 0);
    dup_check = range_list_ptr->names;
    while (dup_check < cur_name_str) {
      if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	sprintf(logbuf, "Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	logprintb();
	return RET_INVALID_CMDLINE;
      }
      dup_check = &(dup_check[name_max_len]);
    }
    cur_name_str = &(cur_name_str[name_max_len]);
    if (range_end) {
      *cur_name_starts_range++ = 1;
      memcpyx(cur_name_str, range_end, re_len, 0);
      dup_check = range_list_ptr->names;
      while (dup_check < cur_name_str) {
	if (!memcmp(dup_check, cur_name_str, rs_len + 1)) {
	  sprintf(logbuf, "Error: Duplicate %s parameter '%s'.\n", argv[0], cur_name_str);
	  logprintb();
	  return RET_INVALID_CMDLINE;
	}
        dup_check = &(dup_check[name_max_len]);
      }
      cur_name_str = &(cur_name_str[name_max_len]);
      *cur_name_starts_range++ = 0;
    } else {
      *cur_name_starts_range++ = 0;
    }
  }
}

void invalid_arg(char* argv) {
  sprintf(logbuf, "Error: Unrecognized flag ('%s').%s%s", argv, (argv[0] == '-')? "" : "  All flags must be preceded by 1-2 dashes.", errstr_append);
}

void print_ver() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

char extract_char_param(char* ss) {
  // maps c, 'c', and "c" to c, and anything else to the null char.  This is
  // intended to support e.g. always using '#' to designate a # parameter
  // without worrying about differences between shells.
  char cc = ss[0];
  if (((cc == '\'') || (cc == '"')) && (ss[1]) && (ss[2] == cc) && (!ss[3])) {
    return ss[1];
  } else if (cc && (!ss[1])) {
    return cc;
  } else {
    return '\0';
  }
}

int32_t alloc_string(char** sbuf, char* source) {
  uint32_t slen = strlen(source) + 1;
  *sbuf = (char*)malloc(slen * sizeof(char));
  if (!(*sbuf)) {
    return -1;
  }
  memcpy(*sbuf, source, slen);
  return 0;
}

int32_t alloc_fname(char** fnbuf, char* source, char* argptr, uint32_t extra_size) {
  uint32_t slen = strlen(source) + 1;
  if (slen > (FNAMESIZE - extra_size)) {
    sprintf(logbuf, "Error: --%s filename too long.\n", argptr);
    logprintb();
    return RET_OPEN_FAIL;
  }
  *fnbuf = (char*)malloc((slen + extra_size) * sizeof(char));
  if (!(*fnbuf)) {
    return RET_NOMEM;
  }
  memcpy(*fnbuf, source, slen);
  return 0;
}

int32_t alloc_and_flatten(char** flattened_buf_ptr, char** sources, uint32_t param_ct) {
  uint32_t totlen = 1;
  char* bufptr;
  uint32_t param_idx;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    totlen += 1 + strlen(sources[param_idx]);
  }
  bufptr = (char*)malloc(totlen);
  if (!bufptr) {
    return RET_NOMEM;
  }
  *flattened_buf_ptr = bufptr;
  for (param_idx = 0; param_idx < param_ct; param_idx++) {
    bufptr = strcpyax(bufptr, sources[param_idx], '\0');
  }
  *bufptr = '\0';
  return 0;
}

int32_t alloc_2col(Two_col_params** tcbuf, char** params_ptr, char* argptr, uint32_t param_ct) {
  uint32_t slen = strlen(*params_ptr) + 1;
  int32_t ii;
  char cc;
  if (slen > FNAMESIZE) {
    sprintf(logbuf, "Error: --%s filename too long.\n", argptr);
    logprintb();
    return RET_OPEN_FAIL;
  }
  *tcbuf = (Two_col_params*)malloc(sizeof(Two_col_params) + slen);
  if (!(*tcbuf)) {
    return RET_NOMEM;
  }
  memcpy((*tcbuf)->fname, params_ptr[0], slen);
  (*tcbuf)->skip = 0;
  (*tcbuf)->skipchar = '\0';
  if (param_ct > 1) {
    ii = atoi(params_ptr[1]);
    if (ii < 1) {
      sprintf(logbuf, "Error: Invalid --%s column number.\n", argptr);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    (*tcbuf)->colx = ii;
    if (param_ct > 2) {
      ii = atoi(params_ptr[2]);
      if (ii < 1) {
	sprintf(logbuf, "Error: Invalid --%s marker ID column number.\n", argptr);
	logprintb();
	return RET_INVALID_FORMAT;
      }
      (*tcbuf)->colid = ii;
      if (param_ct == 4) {
	cc = params_ptr[3][0];
	if ((cc < '0') || (cc > '9')) {
	  cc = extract_char_param(params_ptr[3]);
	  if (!cc) {
            goto alloc_2col_invalid_skip;
	  }
	  (*tcbuf)->skipchar = cc;
	} else {
	  if (atoiz(params_ptr[3], &ii)) {
	  alloc_2col_invalid_skip:
	    sprintf(logbuf, "Error: Invalid --%s skip parameter.  This needs to either be a\nsingle character (usually '#') which, when present at the start of a line,\nindicates it should be skipped; or the number of initial lines to skip.  (Note\nthat in shells such as bash, '#' is a special character that must be\nsurrounded by single- or double-quotes to be parsed correctly.)\n", argptr);
	    logprintb();
	    return RET_INVALID_FORMAT;
	  }
	  (*tcbuf)->skip = ii;
	}
      }
    } else {
      (*tcbuf)->colid = 1;
    }
    if ((*tcbuf)->colx == (*tcbuf)->colid) {
      sprintf(logbuf, "Error: Column numbers for --%s cannot be equal.%s", argptr, errstr_append);
      logprintb();
      return RET_INVALID_FORMAT;
    }
  } else {
    (*tcbuf)->colx = 2;
    (*tcbuf)->colid = 1;
  }
  return 0;
}

int32_t flag_match(const char* to_match, uint32_t* cur_flag_ptr, uint32_t flag_ct, char* flag_buf) {
  int32_t ii;
  while (*cur_flag_ptr < flag_ct) {
    ii = strcmp(to_match, &(flag_buf[(*cur_flag_ptr) * MAX_FLAG_LEN]));
    if (ii < 0) {
      return 0;
    }
    *cur_flag_ptr += 1;
    if (!ii) {
      flag_buf[((*cur_flag_ptr) - 1) * MAX_FLAG_LEN] = '\0';
      return 1;
    }
  }
  return 0;
}

uint32_t species_flag(uint32_t* species_code_ptr, uint32_t new_code) {
  if (*species_code_ptr) {
    logprint("Error: Multiple chromosome set flags.\n");
    return 1;
  }
  *species_code_ptr = new_code;
  return 0;
}

// these need global scope to stay around on all systems
const char species_singular_constants[][7] = {"person", "cow", "dog", "horse", "mouse", "plant", "sheep", "sample"};
const char species_plural_constants[][8] = {"people", "cows", "dogs", "horses", "mice", "plants", "sheep", "samples"};

int32_t init_delim_and_species(uint32_t flag_ct, char* flag_buf, uint32_t* flag_map, int32_t argc, char** argv, char* range_delim_ptr, Chrom_info* chrom_info_ptr) {
  // human: 22, X, Y, XY, MT
  // cow: 29, X, Y
  // dog: 38, X, Y, XY
  // horse: 31, X, Y
  // mouse: 19, X, Y
  // rice: 12
  // sheep: 26, X, Y
  const int32_t species_x_code[] = {23, 30, 39, 32, 20, -1, 27};
  const int32_t species_y_code[] = {24, 31, 40, 33, 21, -1, 28};
  const int32_t species_xy_code[] = {25, -1, 41, -1, -1, -1, -1};
  const int32_t species_mt_code[] = {26, -1, -1, -1, -1, -1, -1};
  const uint32_t species_max_code[] = {26, 31, 41, 33, 21, 12, 28};
  uint32_t species_code = SPECIES_HUMAN;
  uint32_t flag_idx = 0;
  uint32_t retval = 0;
  int32_t cur_arg;
  uint32_t param_ct;
  int32_t ii;
  uint32_t param_idx;
  fill_ulong_zero(chrom_info_ptr->haploid_mask, CHROM_MASK_WORDS);
  fill_ulong_zero(chrom_info_ptr->chrom_mask, CHROM_MASK_WORDS);
  if (flag_match("autosome-num", &flag_idx, flag_ct, flag_buf)) {
    species_code = SPECIES_UNKNOWN;
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    ii = atoi(argv[cur_arg + 1]);
    if ((ii < 1) || (ii > 59)) {
      sprintf(logbuf, "Error: Invalid --autosome-num parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    chrom_info_ptr->x_code = ii + 1;
    chrom_info_ptr->y_code = -1;
    chrom_info_ptr->xy_code = -1;
    chrom_info_ptr->mt_code = -1;
    chrom_info_ptr->max_code = ii + 1;
    chrom_info_ptr->autosome_ct = ii;
    set_bit(chrom_info_ptr->haploid_mask, ii + 1);
  }
  if (flag_match("chr-set", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_UNKNOWN)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 5)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    ii = atoi(argv[cur_arg + 1]);
    if ((!ii) || (ii > 59) || (ii < -59)) {
      sprintf(logbuf, "Error: Invalid --chr-set parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    if (ii < 0) {
      if (param_ct > 1) {
	sprintf(logbuf, "Error: --chr-set does not accept multiple parameters in haploid mode.%s", errstr_append);
	goto init_delim_and_species_ret_INVALID_CMDLINE_2;
      }
      ii = -ii;
      chrom_info_ptr->autosome_ct = ii;
      chrom_info_ptr->x_code = -1;
      chrom_info_ptr->y_code = -1;
      chrom_info_ptr->xy_code = -1;
      chrom_info_ptr->mt_code = -1;
      chrom_info_ptr->max_code = ii;
      fill_bits(chrom_info_ptr->haploid_mask, 0, ii + 1);
    } else {
      chrom_info_ptr->autosome_ct = ii;
      chrom_info_ptr->x_code = ii + 1;
      chrom_info_ptr->y_code = ii + 2;
      chrom_info_ptr->xy_code = ii + 3;
      chrom_info_ptr->mt_code = ii + 4;
      set_bit(chrom_info_ptr->haploid_mask, ii + 1);
      set_bit(chrom_info_ptr->haploid_mask, ii + 2);
      set_bit(chrom_info_ptr->haploid_mask, ii + 4);
      for (param_idx = 2; param_idx <= param_ct; param_idx++) {
	if (!strcmp(argv[cur_arg + param_idx], "no-x")) {
	  chrom_info_ptr->x_code = -1;
	  clear_bit(chrom_info_ptr->haploid_mask, ii + 1);
	} else if (!strcmp(argv[cur_arg + param_idx], "no-y")) {
	  chrom_info_ptr->y_code = -1;
	  clear_bit(chrom_info_ptr->haploid_mask, ii + 2);
	} else if (!strcmp(argv[cur_arg + param_idx], "no-xy")) {
	  chrom_info_ptr->xy_code = -1;
	} else if (!strcmp(argv[cur_arg + param_idx], "no-mt")) {
	  chrom_info_ptr->mt_code = -1;
	  clear_bit(chrom_info_ptr->haploid_mask, ii + 4);
	} else {
	  sprintf(logbuf, "Error: Invalid --chr-set parameter '%s'.%s", argv[cur_arg + param_idx], errstr_append);
	  goto init_delim_and_species_ret_INVALID_CMDLINE_2;
	}
      }
      if (chrom_info_ptr->mt_code != -1) {
	chrom_info_ptr->max_code = ii + 4;
      } else if (chrom_info_ptr->xy_code != -1) {
	chrom_info_ptr->max_code = ii + 3;
      } else if (chrom_info_ptr->y_code != -1) {
	chrom_info_ptr->max_code = ii + 2;
      } else if (chrom_info_ptr->x_code != -1) {
	chrom_info_ptr->max_code = ii + 1;
      } else {
	chrom_info_ptr->max_code = ii;
      }
    }
  }
  if (flag_match("cow", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_COW)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --cow doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("d", &flag_idx, flag_ct, flag_buf)) {
    // moved here to support --covar-name + --d
    cur_arg = flag_map[flag_idx - 1];
    param_ct = param_count(argc, argv, cur_arg);
    if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
    *range_delim_ptr = extract_char_param(argv[cur_arg + 1]);
    if (!(*range_delim_ptr)) {
      sprintf(logbuf, "Error: --d parameter too long (must be a single character).%s", errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    } else if ((*range_delim_ptr == '-') || (*range_delim_ptr == ',')) {
      sprintf(logbuf, "Error: --d parameter cannot be '-' or ','.%s", errstr_append);
      goto init_delim_and_species_ret_INVALID_CMDLINE_2;
    }
  }
  if (flag_match("dog", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_DOG)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --dog doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("horse", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_HORSE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --horse doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("mouse", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_MOUSE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --mouse doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("rice", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_RICE)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --rice doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  if (flag_match("sheep", &flag_idx, flag_ct, flag_buf)) {
    if (species_flag(&species_code, SPECIES_SHEEP)) {
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
    if (param_count(argc, argv, flag_map[flag_idx - 1])) {
      logprint("Error: --sheep doesn't accept parameters.\n");
      goto init_delim_and_species_ret_INVALID_CMDLINE;
    }
  }
  chrom_info_ptr->species = species_code;
  chrom_info_ptr->is_include_stack = 0;
  if (species_code != SPECIES_UNKNOWN) {
    chrom_info_ptr->x_code = species_x_code[species_code];
    chrom_info_ptr->y_code = species_y_code[species_code];
    chrom_info_ptr->xy_code = species_xy_code[species_code];
    chrom_info_ptr->mt_code = species_mt_code[species_code];
    chrom_info_ptr->max_code = species_max_code[species_code];
  }
  g_species_singular = species_singular_constants[species_code];
  g_species_plural = species_plural_constants[species_code];
  switch (species_code) {
  case SPECIES_HUMAN:
    chrom_info_ptr->autosome_ct = 22;
    chrom_info_ptr->haploid_mask[0] = 0x5800000;
    break;
  case SPECIES_COW:
    chrom_info_ptr->autosome_ct = 29;
    chrom_info_ptr->haploid_mask[0] = 0xc0000000LU;
    break;
  case SPECIES_DOG:
    chrom_info_ptr->autosome_ct = 38;
#ifdef __LP64__
    chrom_info_ptr->haploid_mask[0] = 0x18000000000LLU;
#else
    chrom_info_ptr->haploid_mask[1] = 0x180;
#endif
    break;
  case SPECIES_HORSE:
    chrom_info_ptr->autosome_ct = 31;
#ifdef __LP64__
    chrom_info_ptr->haploid_mask[0] = 0x300000000LLU;
#else
    chrom_info_ptr->haploid_mask[1] = 3;
#endif
    break;
  case SPECIES_MOUSE:
    chrom_info_ptr->autosome_ct = 19;
    chrom_info_ptr->haploid_mask[0] = 0x300000;
    break;
  case SPECIES_RICE:
    chrom_info_ptr->autosome_ct = 12;
    chrom_info_ptr->haploid_mask[0] = 0x1fff;
    break;
  case SPECIES_SHEEP:
    chrom_info_ptr->autosome_ct = 26;
    chrom_info_ptr->haploid_mask[0] = 0x18000000;
    break;
  }
  while (0) {
  init_delim_and_species_ret_INVALID_CMDLINE_2:
    logprintb();
  init_delim_and_species_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

void fill_chrom_mask(Chrom_info* chrom_info_ptr) {
  if (chrom_info_ptr->species != SPECIES_UNKNOWN) {
    fill_bits(chrom_info_ptr->chrom_mask, 0, chrom_info_ptr->max_code + 1);
  } else {
    fill_bits(chrom_info_ptr->chrom_mask, 0, chrom_info_ptr->autosome_ct + 1);
    // --chr-set support
    if (chrom_info_ptr->x_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->x_code);
    }
    if (chrom_info_ptr->y_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->y_code);
    }
    if (chrom_info_ptr->xy_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->xy_code);
    }
    if (chrom_info_ptr->mt_code != -1) {
      set_bit(chrom_info_ptr->chrom_mask, chrom_info_ptr->mt_code);
    }
  }
}

int32_t recode_type_set(uint32_t* recode_modifier_ptr, uint32_t cur_code) {
  if (*recode_modifier_ptr & (RECODE_TYPEMASK - cur_code)) {
    sprintf(logbuf, "Error: Conflicting --recode modifiers.%s", errstr_append);
    return -1;
  }
  *recode_modifier_ptr |= cur_code;
  return 0;
}

void range_list_init(Range_list* range_list_ptr) {
  range_list_ptr->names = NULL;
  range_list_ptr->starts_range = NULL;
  range_list_ptr->name_ct = 0;
  range_list_ptr->name_max_len = 0;
}

void free_range_list(Range_list* range_list_ptr) {
  free_cond(range_list_ptr->names);
  free_cond(range_list_ptr->starts_range);
}

int32_t main(int32_t argc, char** argv) {
  char* outname_end = NULL;
  char** subst_argv = NULL;
  char* script_buf = NULL;
  char* rerun_buf = NULL;
  char* flag_buf = NULL;
  uint32_t* flag_map = NULL;
  char* makepheno_str = NULL;
  char* phenoname_str = NULL;
  Two_col_params* a1alleles = NULL;
  Two_col_params* a2alleles = NULL;
  char* filtervals_flattened = NULL;
  char* evecname = NULL;
  char* filtername = NULL;
  char* read_dists_fname = NULL;
  char* read_dists_id_fname = NULL;
  char* freqname = NULL;
  char* extractname = NULL;
  char* excludename = NULL;
  char* keepname = NULL;
  char* removename = NULL;
  char* keepfamname = NULL;
  char* removefamname = NULL;
  char* phenoname = NULL;
  char* recode_allele_name = NULL;
  char* lgen_reference_fname = NULL;
  char* covar_fname = NULL;
  char* set_fname = NULL;
  char* subset_fname = NULL;
  char* update_alleles_fname = NULL;
  Two_col_params* update_chr = NULL;
  Two_col_params* update_cm = NULL;
  Two_col_params* update_map = NULL;
  Two_col_params* update_name = NULL;
  char* update_ids_fname = NULL;
  char* update_parents_fname = NULL;
  char* update_sex_fname = NULL;
  char* loop_assoc_fname = NULL;
  char* flip_fname = NULL;
  char* flip_subset_fname = NULL;
  char* read_genome_fname = NULL;
  char* condition_mname = NULL;
  char* condition_fname = NULL;
  int32_t retval = 0;
  uint32_t load_params = 0; // describes what file parameters have been provided
  uint32_t load_rare = 0;
  uint32_t fam_cols = FAM_COL_13456;
  uint32_t mpheno_col = 0;
  uint32_t mwithin_col = 0;
  uint64_t misc_flags = 0;
  double thin_keep_prob = 1.0;
  uint32_t min_bp_space = 0;
  double exponent = 0.0;
  double min_maf = 0.0;
  double max_maf = 0.5;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  double rel_cutoff = 0.025;
  uint32_t cur_arg = 1;
  uint64_t calculation_type = 0;
  uint32_t rel_calc_type = 0;
  uint32_t dist_calc_type = 0;
  uint32_t mfilter_col = 0;
  uint32_t pheno_modifier = 0;
  int32_t missing_pheno = -9;
  uintptr_t groupdist_iters = ITERS_DEFAULT;
  uint32_t groupdist_d = 0;
  uintptr_t regress_iters = ITERS_DEFAULT;
  uint32_t regress_d = 0;
  uintptr_t regress_rel_iters = ITERS_DEFAULT;
  uint32_t regress_rel_d = 0;
  double unrelated_herit_tol = 0.0000001;
  double unrelated_herit_covg = 0.45;
  double unrelated_herit_covr = 0.55;
  int32_t ibc_type = 0; // -1 for cov
  uint32_t parallel_idx = 0;
  uint32_t parallel_tot = 1;
  uint32_t sex_missing_pheno = 0;
  uint32_t write_covar_modifier = 0;
  uint32_t write_covar_dummy_max_categories = 49;
  uint32_t model_modifier = 0;
  int32_t model_cell_ct = -1;
  uint32_t gxe_mcovar = 0;
  uint32_t glm_modifier = 0;
  double glm_vif_thresh = 50.0;
  uint32_t glm_xchr_model = 1;
  uint32_t ppc_gap = DEFAULT_PPC_GAP;
  uint32_t* rseeds = NULL;
  uint32_t rseed_ct = 0;
  uint32_t genome_modifier = 0;
  double genome_min_pi_hat = -1.0;
  double genome_max_pi_hat = 1.0;
  FILE* scriptfile = NULL;
  uint32_t ld_window_size = 0;
  uint32_t ld_window_incr = 0;
  double ld_last_param = 0.0;
  uint32_t ld_window_kb = 0;
  uint32_t filter_binary = 0;
  uint32_t regress_pcs_modifier = 0;
  uint32_t max_pcs = MAX_PCS_DEFAULT;
  uint32_t recode_modifier = 0;
  uint32_t allelexxxx = 0;
  uint32_t merge_type = 0;
  uint32_t indiv_sort = 0;
  uint32_t cur_flag = 0;
  uint32_t flag_ct = 0;
  uint32_t dummy_marker_ct = 0;
  uint32_t dummy_indiv_ct = 0;
  uint32_t dummy_flags = 0;
  double dummy_missing_geno = 0.0;
  double dummy_missing_pheno = 0.0;
  char* simulate_fname = NULL;
  uint32_t simulate_flags = 0;
  uint32_t simulate_cases = 1000;
  uint32_t simulate_controls = 1000;
  double simulate_prevalence = 0.01;
  char* simulate_label = NULL;
  double simulate_missing = 0.0;
  uint32_t simulate_qt_indivs = 1000;
  char* markername_from = NULL;
  char* markername_to = NULL;
  char* markername_snp = NULL;
  uint32_t snp_window_size = 0;
  int32_t marker_pos_start = -1;
  int32_t marker_pos_end = -1;
  uint32_t lgen_modifier = 0;
  uint32_t covar_modifier = 0;
  uint32_t update_map_modifier = 0;
  uint32_t model_mperm_val = 0;
  uint32_t glm_mperm_val = 0;
  uint32_t mperm_save = 0;
  uint32_t mperm_val = 0;
  uint32_t aperm_min = 6;
  uint32_t aperm_max = 1000000;
  double aperm_alpha = 0;
  double aperm_beta = 0.0001;
  double aperm_init_interval = 1;
  double aperm_interval_slope = 0.001;
  double ci_size = 0.0;
  double pfilter = 1.0;
  uint32_t perm_batch_size = 0;
  uint32_t mtest_adjust = 0;
  double adjust_lambda = 0.0;
  uint32_t ibs_test_perms = DEFAULT_IBS_TEST_PERMS;
  uint32_t neighbor_n1 = 0;
  uint32_t neighbor_n2 = 0;
  uint32_t cnv_calc_type = 0;
  uint32_t cnv_indiv_mperms = 0;
  uint32_t cnv_test_mperms = 0;
  uint32_t cnv_test_region_mperms = 0;
  uint32_t cnv_enrichment_test_mperms = 0;
  uint32_t cnv_min_seglen = 0;
  uint32_t cnv_max_seglen = 0xffffffffU;
  double cnv_min_score = -INFINITY;
  double cnv_max_score = INFINITY;
  uint32_t cnv_min_sites = 0;
  uint32_t cnv_max_sites = 0xffffffffU;
  uint32_t cnv_intersect_filter_type = 0;
  char* cnv_intersect_filter_fname = NULL;
  char* cnv_subset_fname = NULL;
  uint32_t cnv_overlap_type = 0;
  double cnv_overlap_val = 0.0;
  uint32_t cnv_freq_type = 0;
  uint32_t cnv_freq_val = 0;
  double cnv_freq_val2 = 0.0;
  uint32_t cnv_test_window = 0;
  uint32_t segment_modifier = 0;
  double tail_bottom = 0.0;
  double tail_top = 0.0;
  double lasso_h2 = 0.0;
  double lasso_minlambda = -INFINITY;
  char* segment_spanning_fname = NULL;
  char* missing_code = NULL;
  char range_delim = '-';
  uint32_t modifier_23 = 0;
  double pheno_23 = INFINITY;
  char* fid_23 = NULL;
  char* iid_23 = NULL;
  char* paternal_id_23 = NULL;
  char* maternal_id_23 = NULL;
  char* convert_xy_23 = NULL;
  Ll_str* file_delete_list = NULL;
  uint32_t chrom_flag_present = 0;
  uintptr_t chrom_exclude[CHROM_MASK_INITIAL_WORDS];
  char outname[FNAMESIZE];
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char genname[FNAMESIZE];
  char samplename[FNAMESIZE];
  char mergename1[FNAMESIZE];
  char mergename2[FNAMESIZE];
  char mergename3[FNAMESIZE];
  char output_missing_pheno[32];
#ifdef __APPLE__
  int32_t mib[2];
  size_t sztmp;
#endif
  unsigned char* wkspace_ua;
  char** subst_argv2;
  uint32_t param_ct;
  time_t rawtime;
  char* argptr;
  char* sptr;
  char* bubble;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t num_params;
  int32_t in_param;
  Chrom_info chrom_info;
  Homozyg_info homozyg;
  Cluster_info cluster;
  Range_list snps_range_list;
  Range_list covar_range_list;
  Range_list parameters_range_list;
  Range_list tests_range_list;
  char* argptr2;
  char* flagptr;
  double dxx;
  char cc;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  intptr_t default_alloc_mb;
  int64_t llxx;
  Ll_str* ll_str_ptr;
#if _WIN32
  SYSTEM_INFO sysinfo;
  MEMORYSTATUSEX memstatus;
  DWORD windows_dw; // why the f*** does uint32_t not work?
#endif
  homozyg_init(&homozyg);
  cluster_init(&cluster);
  range_list_init(&snps_range_list);
  range_list_init(&covar_range_list);
  range_list_init(&parameters_range_list);
  range_list_init(&tests_range_list);

  chrom_info.name_ct = 0;
  chrom_info.incl_excl_name_stack = NULL;
  for (uii = 1; uii < (uint32_t)argc; uii++) {
    if ((!strcmp("-script", argv[uii])) || (!strcmp("--script", argv[uii]))) {
      ujj = param_count(argc, argv, uii);
      if (enforce_param_ct_range(ujj, argv[uii], 1, 1)) {
	print_ver();
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      for (ujj = uii + 2; ujj < (uint32_t)argc; ujj++) {
	if ((!strcmp("-script", argv[ujj])) || (!strcmp("--script", argv[ujj]))) {
	  print_ver();
	  printf("Error: Multiple --script flags.  Merge the files into one.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      // logging not yet active, so don't use fopen_checked()
      scriptfile = fopen(argv[uii + 1], "rb");
      if (!scriptfile) {
	print_ver();
	printf(errstr_fopen, argv[uii + 1]);
	goto main_ret_OPEN_FAIL;
      }
      if (fseeko(scriptfile, 0, SEEK_END)) {
	print_ver();
	goto main_ret_READ_FAIL;
      }
      llxx = ftello(scriptfile);
      if (llxx == -1) {
	print_ver();
	goto main_ret_READ_FAIL;
      } else if (llxx > 0x7fffffff) {
	// could actually happen if user enters parameters in the wrong order,
	// so may as well catch it and print a somewhat informative error msg
	print_ver();
        fputs("Error: --script file too large.\n", stdout);
        goto main_ret_NOMEM;
      }
      rewind(scriptfile);
      ujj = (uint32_t)((uint64_t)llxx);
      script_buf = (char*)malloc(ujj);
      if (!script_buf) {
	print_ver();
	goto main_ret_NOMEM;
      }
      ukk = fread(script_buf, 1, ujj, scriptfile);
      if (ukk < ujj) {
	print_ver();
	goto main_ret_READ_FAIL;
      }
      fclose_null(&scriptfile);
      num_params = 0;
      in_param = 0;
      for (ukk = 0; ukk < ujj; ukk++) {
	if (is_space_or_eoln(script_buf[ukk])) {
	  in_param = 0;
	} else if (!in_param) {
	  num_params++;
	  in_param = 1;
	}
      }
      subst_argv = (char**)malloc((num_params + argc - 3) * sizeof(char*));
      num_params = 0;
      in_param = 0;
      for (ukk = 1; ukk < uii; ukk++) {
        subst_argv[num_params++] = argv[ukk];
      }
      for (ukk = 0; ukk < ujj; ukk++) {
	if (is_space_or_eoln(script_buf[ukk])) {
	  if (in_param) {
	    script_buf[ukk] = '\0';
	    in_param = 0;
	  }
	} else if (!in_param) {
	  subst_argv[num_params++] = &(script_buf[ukk]);
	  in_param = 1;
	}
      }
      for (ujj = uii + 2; ujj < (uint32_t)argc; ujj++) {
	subst_argv[num_params++] = argv[ujj];
      }
      argc = num_params;
      cur_arg = 0;
      argv = subst_argv;
    }
  }
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    if ((!strcmp("-rerun", argv[uii])) || (!strcmp("--rerun", argv[uii]))) {
      ujj = param_count(argc, argv, uii);
      if (enforce_param_ct_range(ujj, argv[uii], 0, 1)) {
	print_ver();
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      for (ukk = uii + ujj + 1; ukk < (uint32_t)argc; ukk++) {
	if ((!strcmp("-rerun", argv[ukk])) || (!strcmp("--rerun", argv[ukk]))) {
	  print_ver();
	  fputs("Error: Duplicate --rerun flag.\n", stdout);
	  goto main_ret_INVALID_CMDLINE;
	}
      }
      if (ujj) {
	scriptfile = fopen(argv[uii + 1], "r");
      } else {
        scriptfile = fopen(PROG_NAME_STR ".log", "r");
      }
      if (!scriptfile) {
	print_ver();
	goto main_ret_OPEN_FAIL;
      }
      if (!fgets(tbuf, MAXLINELEN, scriptfile)) {
	print_ver();
	fputs("Error: Empty log file for --rerun.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      if (!fgets(tbuf, MAXLINELEN, scriptfile)) {
	print_ver();
	fputs("Error: Only one line in --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      fclose_null(&scriptfile);
      kk = atoi(tbuf);
      if ((kk < 1) || (kk > MAXLINELEN)) {
	print_ver();
	fputs("Error: Improperly formatted --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      ukk = strlen(tbuf) + 1;
      if (ukk == MAXLINELEN) {
	print_ver();
	fputs("Error: Second line too long in --rerun log file.\n", stdout);
	goto main_ret_INVALID_FORMAT;
      }
      rerun_buf = (char*)malloc(ukk);
      memcpy(rerun_buf, tbuf, ukk);

      memset(tbuf, 1, (uint32_t)kk);
      sptr = next_item_mult(rerun_buf, 2);
      umm = 0;
      ukk = 0;
      do {
	if (no_more_items_kns(sptr)) {
	  print_ver();
	  fputs("Error: Improperly formatted --rerun log file.\n", stdout);
	  goto main_ret_INVALID_FORMAT;
	}
	argptr = is_flag_start(sptr);
	if (argptr) {
          for (unn = cur_arg; unn < (uint32_t)argc; unn++) {
	    argptr2 = is_flag_start(argv[unn]);
	    if (argptr2) {
	      if (!strcmp(argptr, argptr2)) {
		unn = 0xffffffffU;
		break;
	      }
	    }
	  }
          if (unn == 0xffffffffU) {
	    // matching flag, override --rerun
            do {
	      ukk++;
	      tbuf[umm++] = 0;
	      if (umm == (uint32_t)kk) {
		break;
	      }
	      sptr = next_item(sptr);
	    } while (!is_flag(sptr));
	  } else {
	    umm++;
	    sptr = next_item(sptr);
	  }
	} else {
	  umm++;
          sptr = next_item(sptr);
	}
      } while (umm < (uint32_t)kk);
      subst_argv2 = (char**)malloc((argc + kk - ukk - ujj - 1 - cur_arg) * sizeof(char*));
      if (!subst_argv2) {
	print_ver();
	goto main_ret_NOMEM;
      }
      unn = 0;
      for (umm = cur_arg; umm < uii; umm++) {
	subst_argv2[unn++] = argv[umm];
      }
      sptr = next_item_mult(rerun_buf, 2);
      for (umm = 0; umm < (uint32_t)kk; umm++) {
        if (tbuf[umm]) {
	  ukk = strlen_se(sptr);
	  subst_argv2[unn++] = sptr;
	  sptr[ukk] = '\0';
	  if (umm != ((uint32_t)kk) - 1) {
	    sptr = skip_initial_spaces(&(sptr[ukk + 1]));
	  }
	} else {
	  sptr = next_item(sptr);
	}
      }
      for (umm = uii + ujj + 1; umm < (uint32_t)argc; umm++) {
	subst_argv2[unn++] = argv[umm];
      }
      cur_arg = 0;
      argc = unn;
      if (subst_argv) {
	free(subst_argv);
      }
      subst_argv = subst_argv2;
      argv = subst_argv2;
      subst_argv2 = NULL;
    }
  }
  if ((cur_arg < (uint32_t)argc) && (!is_flag(argv[cur_arg]))) {
    print_ver();
    printf("Error: First parameter must be a flag.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE;
  }
  flag_ct = 0;
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    argptr = is_flag_start(argv[uii]);
    if (argptr) {
      if (!strcmp("help", argptr)) {
	print_ver();
	if ((cur_arg != 1) || (uii != 1) || subst_argv) {
	  fputs("--help present, ignoring other flags.\n", stdout);
	}
	retval = disp_help(argc - uii - 1, &(argv[uii + 1]));
	goto main_ret_1;
      }
      if (strlen(argptr) >= MAX_FLAG_LEN) {
	print_ver();
	invalid_arg(argv[uii]);
	fputs(logbuf, stdout);
        goto main_ret_INVALID_CMDLINE;
      }
      flag_ct++;
    }
  }
  if (!flag_ct) {
    print_ver();
    fputs(notestr_null_calc, stdout);
    fputs(cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
    goto main_ret_1;
  }
  flag_buf = (char*)malloc(flag_ct * MAX_FLAG_LEN * sizeof(char));
  flag_map = (uint32_t*)malloc(flag_ct * sizeof(int32_t));
  if ((!flag_buf) || (!flag_map)) {
    print_ver();
    goto main_ret_NOMEM;
  }
  flagptr = flag_buf;
  umm = 0; // parameter count increase due to aliases
  for (uii = cur_arg; uii < (uint32_t)argc; uii++) {
    argptr = is_flag_start(argv[uii]);
    if (argptr) {
      ukk = strlen(argptr) + 1;
      // handle aliases now, so sorting will have the desired effects
      switch (*argptr) {
      case 'Z':
	if (!strcmp(argptr, "Z-genome")) {
	  memcpy(flagptr, "genome gz", 10);
	  umm++;
	  break;
	}
	goto main_flag_copy;
      case 'a':
	if ((ukk == 11) && (!memcmp(argptr, "allele", 6))) {
	  if (match_upper(&(argptr[6]), "ACGT")) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	} else if ((ukk == 12) && (!memcmp(argptr, "allele-", 7))) {
          if (!memcmp(&(argptr[7]), "1234", 4)) {
	    memcpy(flagptr, "allele1234", 11);
	    break;
	  } else if (match_upper(&(argptr[7]), "ACGT")) {
	    memcpy(flagptr, "alleleACGT", 11);
	    break;
	  }
	}
	goto main_flag_copy;
      case 'c':
        if (!strcmp(argptr, "chr-excl")) {
          logprint("Note: --chr-excl flag has been renamed to --not-chr.\n");
	  memcpy(flagptr, "not-chr", 8);
	  break;
	} else if (!strcmp(argptr, "cmh")) {
	  memcpy(flagptr, "mh", 3);
	  break;
	}
	goto main_flag_copy;
      case 'e':
	if (!strcmp(argptr, "extract-snp")) {
	  memcpy(flagptr, "snp", 4);
	  break;
	}
	goto main_flag_copy;
      case 'f':
	if (!strcmp(argptr, "frqx")) {
	  memcpy(flagptr, "freqx", 6);
	  break;
	}
	goto main_flag_copy;
      case 'g':
	if (!strcmp(argptr, "grm-cutoff")) {
          memcpy(flagptr, "rel-cutoff", 11);
	  break;
	}
	goto main_flag_copy;
      case 'k':
	if (!memcmp(argptr, "k", 2)) {
	  memcpy(flagptr, "K", 2);
	  break;
	}
	goto main_flag_copy;
      case 'l':
	if (!strcmp(argptr, "list")) {
	  memcpy(flagptr, "recode list", 12);
	  printf("Note: --list flag deprecated.  Use '--recode list' instead.\n");
	  recode_modifier |= RECODE_LIST;
	  misc_flags |= MISC_SET_HH_MISSING;
	  break;
	} else if (!strcmp(argptr, "load-dists")) {
          memcpy(flagptr, "read-dists", 11);
          printf("Note: --load-dists flag has been renamed to --read-dists.\n");
          break;
	}
	goto main_flag_copy;
      case 'm':
	if (!strcmp(argptr, "missing_code")) {
	  memcpy(flagptr, "missing-code", 13);
	  break;
	} else if (!strcmp(argptr, "mh1")) {
	  memcpy(flagptr, "mh", 3);
	  break;
	}
	goto main_flag_copy;
      case 'n':
	if (!strcmp(argptr, "neighbor")) {
	  memcpy(flagptr, "neighbour", 10);
	  break;
	} else if (!strcmp(argptr, "num_threads")) {
	  memcpy(flagptr, "threads", 8);
	  break;
	}
	goto main_flag_copy;
      case 'r':
	if ((ukk >= 8) && (!memcmp(argptr, "recode", 6))) {
	  ujj = 0; // alias match?
	  argptr2 = &(argptr[6]);
          switch (ukk) {
	  case 8:
            if (tolower(*argptr2) == 'a') {
	      memcpy(flagptr, "recode A", 9);
	      recode_modifier |= RECODE_A;
	      ujj = 1;
	    }
	    break;
	  case 9:
	    if (!memcmp(argptr2, "12", 2)) {
              memcpy(flagptr, "recode 12", 10);
	      recode_modifier |= RECODE_12;
	      ujj = 1;
            } else if (match_upper(argptr2, "AD")) {
              memcpy(flagptr, "recode AD", 10);
	      recode_modifier |= RECODE_AD;
	      ujj = 1;
	    } else if (match_upper(argptr2, "HV")) {
	      memcpy(flagptr, "recode HV-1chr", 15);
	      recode_modifier |= RECODE_HV_1CHR;
              printf("Note: --recodeHV flag deprecated.  Use '--recode HV' or '--recode HV-1chr'.\n");
	      ujj = 2;
	    }
	    break;
	  case 11:
	    if (!memcmp(argptr2, "-vcf", 4)) {
	      memcpy(flagptr, "recode vcf", 11);
	      recode_modifier |= RECODE_VCF;
	      ujj = 1;
	    }
	    break;
          case 12:
            if (!memcmp(argptr2, "-lgen", 5)) {
              memcpy(flagptr, "recode lgen", 12);
	      recode_modifier |= RECODE_LGEN;
	      // backwards compatibility
	      misc_flags |= MISC_SET_HH_MISSING;
              ujj = 1;
	    }
	    break;
	  case 13:
	    if (!memcmp(argptr2, "-rlist", 6)) {
	      memcpy(flagptr, "recode rlist", 13);
	      recode_modifier |= RECODE_RLIST;
	      misc_flags |= MISC_SET_HH_MISSING;
	      ujj = 1;
	    }
	    break;
	  case 14:
	    if (!memcmp(argptr2, "-beagle", 7)) {
	      memcpy(flagptr, "recode beagle", 14);
	      recode_modifier |= RECODE_BEAGLE;
	      ujj = 1;
	    } else if (!memcmp(argptr2, "-bimbam", 7)) {
	      memcpy(flagptr, "recode bimbam-1chr", 19);
	      recode_modifier |= RECODE_BIMBAM_1CHR;
	      misc_flags |= MISC_SET_HH_MISSING;
	      printf("Note: --recode-bimbam flag deprecated.  Use '--recode bimbam' or\n'--recode bimbam-1chr'.\n");
	      ujj = 2;
	    }
	    break;
	  case 17:
	    if (!memcmp(argptr2, "-fastphase", 10)) {
	      memcpy(flagptr, "recode fastphase-1chr", 22);
	      recode_modifier |= RECODE_FASTPHASE_1CHR;
	      misc_flags |= MISC_SET_HH_MISSING;
	      printf("Note: --recode-fastphase flag deprecated.  Use '--recode fastphase' or\n'--recode fastphase-1chr'.\n");
	      ujj = 2;
	    } else if (!memcmp(argptr2, "-structure", 10)) {
	      memcpy(flagptr, "recode structure", 17);
	      recode_modifier |= RECODE_STRUCTURE;
	      misc_flags |= MISC_SET_HH_MISSING;
	      ujj = 1;
	    }
	    break;
	  }
	  if (ujj) {
	    if (ujj == 1) {
	      printf("Note: --%s flag deprecated.  Use '%s ...'.\n", argptr, flagptr);
	    }
	    umm++;
	    break;
	  }
	} else if (!strcmp(argptr, "reference-allele")) {
	  memcpy(flagptr, "a1-allele", 10);
	  break;
	}
	goto main_flag_copy;
      case 't':
        if (!strcmp(argptr, "thread-num")) {
	  memcpy(flagptr, "threads", 8);
	  break;
	}
	goto main_flag_copy;
      case 'u':
	if (!strcmp(argptr, "update-freq")) {
	  memcpy(flagptr, "read-freq", 10);
	  break;
	} else if (!strcmp(argptr, "update-ref-allele")) {
	  // GCTA alias
	  memcpy(flagptr, "a1-allele", 10);
	  break;
	}
	// fall through
      default:
      main_flag_copy:
	memcpy(flagptr, argptr, ukk);
      }
      flagptr = &(flagptr[MAX_FLAG_LEN]);
      flag_map[cur_flag++] = uii;
    }
  }
  sptr = (char*)malloc(flag_ct * MAX_FLAG_LEN);
  if (!sptr) {
    print_ver();
    goto main_ret_NOMEM;
  }
  qsort_ext2(flag_buf, flag_ct, MAX_FLAG_LEN, strcmp_deref, (char*)flag_map, sizeof(int32_t), sptr, MAX_FLAG_LEN);
  free(sptr);
  ujj = strlen_se(flag_buf);
  for (cur_flag = 1; cur_flag < flag_ct; cur_flag++) {
    ukk = strlen_se(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    if ((ujj == ukk) && (!memcmp(&(flag_buf[(cur_flag - 1) * MAX_FLAG_LEN]), &(flag_buf[cur_flag * MAX_FLAG_LEN]), ukk))) {
      flag_buf[cur_flag * MAX_FLAG_LEN + ukk] = '\0'; // just in case of aliases
      print_ver();
      printf("Error: Duplicate --%s flag.\n", &(flag_buf[cur_flag * MAX_FLAG_LEN]));
      goto main_ret_INVALID_CMDLINE;
    }
    ujj = ukk;
  }

  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    if (!memcmp("silent", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 7)) {
      freopen("/dev/null", "w", stdout);
      break;
    }
  }
  print_ver();
  uii = 5;
  memcpy(outname, PROG_NAME_STR, 6);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    ii = memcmp("out", &(flag_buf[cur_flag * MAX_FLAG_LEN]), 4);
    if (!ii) {
      ujj = flag_map[cur_flag];
      ukk = param_count(argc, argv, ujj);
      if (enforce_param_ct_range(ukk, argv[ujj], 1, 1)) {
	fputs(logbuf, stdout);
	goto main_ret_INVALID_CMDLINE;
      }
      if (strlen(argv[ujj + 1]) > (FNAMESIZE - MAX_POST_EXT)) {
	fputs("Error: --out parameter too long.\n", stdout);
	goto main_ret_OPEN_FAIL;
      }
      uii = strlen(argv[ujj + 1]);
      memcpy(outname, argv[ujj + 1], uii + 1);
      outname_end = &(outname[uii]);
    }
    if (ii <= 0) {
      break;
    }
  }
  memcpy(&(outname[uii]), ".log", 5);
  logfile = fopen(outname, "w");
  if (!logfile) {
    printf("Error: Failed to open %s.  Try ", outname);
    if (!memcmp(outname, PROG_NAME_STR, 6)) {
      printf("using --out.\n");
    } else {
      printf("changing the --out parameter.\n");
    }
    goto main_ret_OPEN_FAIL;
  }
  printf("Logging to %s.\n", outname);
  outname[uii] = '\0';

  logstr(ver_str);
  sprintf(logbuf, "\n%d argument%s:", argc + umm - cur_arg, (argc + umm - cur_arg == 1)? "" : "s");
  logstr(logbuf);
  for (cur_flag = 0; cur_flag < flag_ct; cur_flag++) {
    logstr(" --");
    logstr(&(flag_buf[cur_flag * MAX_FLAG_LEN]));
    ii = flag_map[cur_flag] + 1;
    while ((ii < argc) && (!is_flag(argv[ii]))) {
      logstr(" ");
      logstr(argv[ii++]);
    }
  }
#if _WIN32
  windows_dw = 4 * MAXLINELEN + 256;
  if (GetComputerName(tbuf, &windows_dw))
#else
  if (gethostname(tbuf, 4 * MAXLINELEN + 256) != -1)
#endif
  {
    logstr("\nHostname: ");
    logstr(tbuf);
  }
  logstr("\nWorking directory: ");
  getcwd(tbuf, FNAMESIZE);
  logstr(tbuf);
  logstr("\nStart time: ");
  time(&rawtime);
  logstr(ctime(&rawtime));
  logstr("\n");

#if _WIN32
  GetSystemInfo(&sysinfo);
  g_thread_ct = sysinfo.dwNumberOfProcessors;
#else
  ii = sysconf(_SC_NPROCESSORS_ONLN);
  if (ii == -1) {
    g_thread_ct = 1;
  } else {
    g_thread_ct = ii;
  }
#endif
  if (g_thread_ct > 8) {
    if (g_thread_ct > MAX_THREADS) {
      g_thread_ct = MAX_THREADS;
    } else {
      g_thread_ct--;
    }
  }
  memcpy(mapname, PROG_NAME_STR ".map", 10);
  memcpy(pedname, PROG_NAME_STR ".ped", 10);
  famname[0] = '\0';
  genname[0] = '\0';
  samplename[0] = '\0';
  memcpyl3(output_missing_pheno, "-9");
  // stuff that must be processed before regular alphabetical loop
  retval = init_delim_and_species(flag_ct, flag_buf, flag_map, argc, argv, &range_delim, &chrom_info);
  if (retval) {
    goto main_ret_1;
  }
  fill_ulong_zero(chrom_exclude, CHROM_MASK_INITIAL_WORDS);
  cur_flag = 0;
  do {
    argptr = &(flag_buf[cur_flag * MAX_FLAG_LEN]);
    if (!(*argptr)) {
      // preprocessed
      continue;
    }
    argptr2 = &(argptr[1]);
    cur_arg = flag_map[cur_flag];
    param_ct = param_count(argc, argv, cur_arg);
    switch (*argptr) {
    case '1':
      if (*argptr2 == '\0') {
	misc_flags |= MISC_AFFECTION_01;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case '2':
      if (!memcmp(argptr2, "3file", 6)) {
	if (chrom_info.species != SPECIES_HUMAN) {
	  logprint("Error: --23file cannot be used with a nonhuman species flag.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 7)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = strlen(argv[cur_arg + 1]);
	if (ii > FNAMESIZE - 1) {
	  logprint("Error: --23file filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
        memcpy(pedname, argv[cur_arg + 1], ii + 1);
	if (param_ct > 1) {
	  if (strchr(argv[cur_arg + 2], ' ')) {
	    logprint("Error: Space present in --23file family ID.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (alloc_string(&fid_23, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	  if (param_ct > 2) {
	    if (strchr(argv[cur_arg + 3], ' ')) {
	      logprint("Error: Space present in --23file individual ID.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (alloc_string(&iid_23, argv[cur_arg + 3])) {
	      goto main_ret_NOMEM;
	    }
	    if (param_ct > 3) {
	      cc = extract_char_param(argv[cur_arg + 4]);
	      if ((cc == 'M') || (cc == 'm') || (cc == '1')) {
		modifier_23 |= M23_MALE;
	      } else if ((cc == 'F') || (cc == 'f') || (cc == '2')) {
		modifier_23 |= M23_FEMALE;
	      } else if (cc == '0') {
		modifier_23 |= M23_FORCE_MISSING_SEX;
	      } else if ((cc != 'I') && (cc != 'i')) {
		logprint("Error: Invalid --23file sex parameter (M or 1 = male, F or 2 = female,\nI = infer from data, 0 = force missing).\n");
		goto main_ret_INVALID_CMDLINE;
	      }
	      if (param_ct > 4) {
		if (scan_double(argv[cur_arg + 5], &pheno_23)) {
		  sprintf(logbuf, "Error: Invalid --23file phenotype '%s'.%s", argv[cur_arg + 5], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
		if (param_ct > 5) {
		  if (strchr(argv[cur_arg + 6], ' ')) {
		    logprint("Error: Space present in --23file paternal ID.\n");
		    goto main_ret_INVALID_CMDLINE;
		  }
		  if (alloc_string(&paternal_id_23, argv[cur_arg + 6])) {
		    goto main_ret_NOMEM;
		  }
		  if (param_ct > 6) {
		    if (strchr(argv[cur_arg + 7], ' ')) {
		      logprint("Error: Space present in --23file maternal ID.\n");
		      goto main_ret_INVALID_CMDLINE;
		    }
		    if (alloc_string(&maternal_id_23, argv[cur_arg + 7])) {
		      goto main_ret_NOMEM;
		    }
		  }
		}
	      }
	    }
	  }
	}
	load_rare = LOAD_RARE_23;
      } else if (!memcmp(argptr2, "3file-convert-xy", 17)) {
	if (load_rare != LOAD_RARE_23) {
	  logprint("Error: --23file-convert-xy must be used with --23file.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  retval = alloc_fname(&convert_xy_23, argv[cur_arg + 1], argptr, 0);
	  if (retval) {
	    goto main_ret_1;
	  }
	} else if (modifier_23 & M23_FEMALE) {
	  sprintf(logbuf, "Error: --23file-convert-xy requires a parameter when used on a nonmale genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (!(modifier_23 & M23_SEX)) {
	  logprint("Note: Inferring male sex because --23file-convert-xy had no parameter.\n");
	  modifier_23 |= M23_MALE;
	}
        modifier_23 |= M23_CONVERT_XY;
      } else if (!memcmp(argptr2, "3file-make-xylist", 18)) {
	if (load_rare != LOAD_RARE_23) {
	  logprint("Error: --23file-make-xylist must be used with --23file.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (modifier_23 & M23_FEMALE) {
	  sprintf(logbuf, "Error: --23file-make-xylist cannot be used on a female genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (convert_xy_23) {
	  logprint("Error: --23file-make-xylist cannot be used with a --23file-convert-xy file.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!(modifier_23 & M23_SEX)) {
	  logprint("Note: Inferring male sex due to use of --23file-make-xylist.\n");
	  modifier_23 |= M23_MALE;
	}
	modifier_23 |= M23_MAKE_XYLIST;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'K':
      if (*argptr2 == '\0') {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
          sprintf(logbuf, "Error: Invalid --K cluster count '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.min_ct = ii;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'a':
      if (!memcmp(argptr2, "utosome", 8)) {
	fill_bits(chrom_info.chrom_mask, 1, chrom_info.autosome_ct);
	chrom_info.is_include_stack = 1;
	chrom_flag_present = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "utosome-xy", 11)) {
	if (chrom_flag_present) {
          logprint("Error: --autosome-xy cannot be used with --autosome.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (chrom_info.xy_code == -1) {
	  sprintf(logbuf, "Error: --autosome-xy used with a species lacking an XY region.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	fill_bits(chrom_info.chrom_mask, 1, chrom_info.autosome_ct);
	set_bit(chrom_info.chrom_mask, chrom_info.xy_code);
	chrom_info.is_include_stack = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llow-extra-chr", 15)) {
	if (load_rare == LOAD_RARE_23) {
	  logprint("Error: --allow-extra-chr cannot currently be used with --23file.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  if (memcmp("0", argv[cur_arg + 1], 2)) {
            sprintf(logbuf, "Error: Invalid --allow-extra-chr parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_ZERO_EXTRA_CHROMS;
	}
	misc_flags |= MISC_ALLOW_EXTRA_CHROMS;
      } else if (!memcmp(argptr2, "llow-no-sex", 12)) {
        sex_missing_pheno |= ALLOW_NO_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ll", 3)) {
	logprint("Note: --all flag has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "llele1234", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  if (strcmp("multichar", argv[cur_arg + 1])) {
	    sprintf(logbuf, "Error: Invalid --allele1234 parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  allelexxxx = ALLELE_RECODE_MULTICHAR;
	} else {
	  allelexxxx = ALLELE_RECODE;
	}
      } else if (!memcmp(argptr2, "lleleACGT", 9)) {
	if (allelexxxx) {
	  logprint("Error: --allele1234 and --alleleACGT cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  if (strcmp("multichar", argv[cur_arg + 1])) {
	    sprintf(logbuf, "Error: Invalid --alleleACGT parameter '%s'.%s\n", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  allelexxxx = ALLELE_RECODE_ACGT | ALLELE_RECODE_MULTICHAR;
	} else {
	  allelexxxx = ALLELE_RECODE_ACGT;
	}
      } else if (!memcmp(argptr2, "llele-count", 12)) {
	lgen_modifier |= LGEN_ALLELE_COUNT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ll-pheno", 9)) {
	pheno_modifier |= PHENO_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ssoc", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 6)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "counts")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' modifier does not make sense with 'counts'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_ASSOC_COUNTS;
	  } else if (!strcmp(argv[cur_arg + uii], "fisher")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with 'fisher'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: --assoc 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    if (model_modifier & MODEL_QMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with 'genedrop'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "p2")) {
	    model_modifier |= MODEL_ASSOC_P2;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (model_modifier & MODEL_PERM) {
	      sprintf(logbuf, "Error: --assoc 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --assoc 'mperm' modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    kk = atoi(&(argv[cur_arg + uii][6]));
	    if (kk < 1) {
	      sprintf(logbuf, "Error: Invalid --assoc mperm parameter '%s'.%s", &(argv[cur_arg + uii][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_mperm_val = (uint32_t)kk;
	    model_modifier |= MODEL_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "qt-means")) {
	    if (model_modifier & MODEL_DMASK) {
	      sprintf(logbuf, "Error: --assoc 'qt-means' does not make sense with a case/control-specific\nmodifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_QT_MEANS;
	  } else if (!strcmp(argv[cur_arg + uii], "lin")) {
	    if (model_modifier & MODEL_DMASK) {
	      sprintf(logbuf, "Error: --assoc 'lin' does not make sense with a case/control-specific modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_LIN;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logprint("Error: Improper --assoc mperm syntax.  (Use '--assoc mperm=[value]'.)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    sprintf(logbuf, "Error: Invalid --assoc parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	model_modifier |= MODEL_ASSOC;
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "djust", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mtest_adjust = 1;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "gc")) {
	    mtest_adjust |= ADJUST_GC;
	  } else if (!strcmp(argv[cur_arg + uii], "log10")) {
	    mtest_adjust |= ADJUST_LOG10;
	  } else if (!strcmp(argv[cur_arg + uii], "qq-plot")) {
	    mtest_adjust |= ADJUST_QQ;
	  } else {
	    sprintf(logbuf, "Error: Invalid --adjust parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 6, 6)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --aperm min permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	aperm_min = ii + 1;
	ii = atoi(argv[cur_arg + 2]);
	if ((ii <= ((int32_t)aperm_min)) || (ii > APERM_MAX)) {
	  sprintf(logbuf, "Error: Invalid --aperm max permutation count '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	aperm_max = ii;
	if (scan_double(argv[cur_arg + 3], &aperm_alpha)) {
	  sprintf(logbuf, "Error: Invalid --aperm alpha threshold '%s'.%s", argv[cur_arg + 3], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 4], &aperm_beta)) {
	  sprintf(logbuf, "Error: Invalid --aperm beta '%s'.%s", argv[cur_arg + 4], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 5], &aperm_init_interval)) {
	  sprintf(logbuf, "Error: Invalid --aperm initial pruning interval '%s'.%s", argv[cur_arg + 5], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((aperm_init_interval < 1) || (aperm_init_interval > 1000000)) {
	  sprintf(logbuf, "Error: Invalid --aperm initial pruning interval '%s'.%s", argv[cur_arg + 5], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 6], &aperm_interval_slope)) {
	  sprintf(logbuf, "Error: Invalid --aperm pruning interval slope '%s'.%s", argv[cur_arg + 6], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((aperm_interval_slope < 0) || (aperm_interval_slope > 1)) {
	  sprintf(logbuf, "Error: Invalid --aperm pruning interval slope '%s'.%s", argv[cur_arg + 6], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "1-allele", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_2col(&a1alleles, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "2-allele", 9)) {
	if (a1alleles) {
	  logprint("Error: --a2-allele cannot be used with --a1-allele.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_2col(&a2alleles, &(argv[cur_arg + 1]), argptr, param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'b':
      if (!memcmp(argptr2, "file", 5)) {
	load_params |= 8;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 5)) {
	    logprint("Error: --bfile parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  sptr = (char*)PROG_NAME_STR;
	}
	if (!(load_params & 16)) {
	  memcpy(strcpya(pedname, sptr), ".bed", 5);
	}
	memcpy(strcpya(mapname, sptr), ".bim", 5);
	memcpy(strcpya(famname, sptr), ".fam", 5);
      } else if (!memcmp(argptr2, "ed", 3)) {
	load_params |= 16;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --bed parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "im", 3)) {
	load_params |= 32;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --bim parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "merge", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 2) {
	  sprintf(logbuf, "Error: --bmerge must have exactly 1 or 3 parameters.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (param_ct == 3) {
	  if (++jj > FNAMESIZE) {
	    logprint("Error: --bmerge .bed filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --bmerge .bim filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	  jj = strlen(argv[cur_arg + 3]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --bmerge .fam filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename3, argv[cur_arg + 3], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logprint("Error: --bmerge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(mergename1, argv[cur_arg + 1], jj), ".bed", 5);
	  memcpy(memcpya(mergename2, argv[cur_arg + 1], jj), ".bim", 5);
	  memcpy(memcpya(mergename3, argv[cur_arg + 1], jj), ".fam", 5);
	}
	calculation_type |= CALC_MERGE;
	merge_type |= MERGE_BINARY;
      } else if (!memcmp(argptr2, "p-space", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --bp-space minimum bp distance '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        min_bp_space = ii;
      } else if (!memcmp(argptr2, "eta", 4)) {
	logprint("Note: --beta flag deprecated.  Use e.g. '--logistic beta'.\n");
	glm_modifier |= GLM_BETA;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "d", 2)) {
	logprint("Note: --bd flag deprecated.  Use '--mh bd'.\n");
	calculation_type |= CALC_CMH;
	misc_flags |= MISC_CMH_BD;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'c':
      if (!memcmp(argptr2, "hr", 3)) {
	if (chrom_flag_present) {
	  sprintf(logbuf, "Error: --chr cannot be used with --autosome[-xy].%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = parse_chrom_ranges(param_ct, '-', &(argv[cur_arg]), chrom_info.chrom_mask, &chrom_info, (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1, argptr);
	if (retval) {
	  goto main_ret_1;
	}
	chrom_info.is_include_stack = 1;
	chrom_flag_present = 1;
      } else if (!memcmp(argptr2, "ompound-genotypes", 18)) {
	logprint("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nand .lgen files are optional if all alleles are single-character).\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ompress", 8)) {
	logprint("Error: --compress flag retired.  Use e.g. 'gzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ounts", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  if (model_modifier & MODEL_QMASK) {
	    sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with --counts.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  logprint("Note: --counts flag deprecated.  Use '--assoc counts' instead.\n");
          model_modifier |= MODEL_ASSOC_COUNTS;
	} else {
	  logprint("Note: --counts flag deprecated.  Use '--freq counts' or --freqx instead.\n");
	}
	misc_flags |= MISC_FREQ_COUNTS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ovar", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "keep-pheno-on-missing-cov")) {
	    uii = 2;
	  } else if (strcmp(argv[cur_arg + 2], "keep-pheno-on-missing-cov")) {
	    sprintf(logbuf, "Error: Invalid --covar parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
            goto main_ret_INVALID_CMDLINE_3;
	  }
          covar_modifier |= COVAR_KEEP_PHENO_ON_MISSING_COV;
	}
	retval = alloc_fname(&covar_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ovar-name", 10)) {
	if (!covar_fname) {
	  logprint("Error: --covar-name must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &covar_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
	covar_modifier |= COVAR_NAME;
      } else if (!memcmp(argptr2, "ovar-number", 12)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (covar_modifier & COVAR_NAME) {
	  logprint("Error: --covar-number cannot be used with --covar-name.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (!covar_fname) {
	  logprint("Error: --covar-number must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &covar_range_list, 1);
	if (retval) {
	  goto main_ret_1;
	}
	covar_modifier |= COVAR_NUMBER;
      } else if (!memcmp(argptr2, "ell", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &model_cell_ct)) {
	  sprintf(logbuf, "Error: Invalid --cell parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "i", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ci parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < 0.01) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --ci confidence interval size s must satisfy 0.01 <= s < 1.%s", errstr_append);
	}
	ci_size = dxx;
      } else if (!memcmp(argptr2, "luster", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cc")) {
            cluster.modifier |= CLUSTER_CC;
	  } else if (!strcmp(argv[cur_arg + uii], "group-avg")) {
	    if (cluster.modifier & CLUSTER_OLD_TIEBREAKS) {
              sprintf(logbuf, "Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    cluster.modifier |= CLUSTER_GROUP_AVG;
	  } else if (!strcmp(argv[cur_arg + uii], "missing")) {
	    cluster.modifier |= CLUSTER_MISSING;
	  } else if (!strcmp(argv[cur_arg + uii], "only2")) {
	    cluster.modifier |= CLUSTER_ONLY2;
	  } else if (!strcmp(argv[cur_arg + uii], "old-tiebreaks")) {
	    if (cluster.modifier & CLUSTER_GROUP_AVG) {
              sprintf(logbuf, "Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.%s", errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    }
	    cluster.modifier |= CLUSTER_OLD_TIEBREAKS;
	  } else {
            sprintf(logbuf, "Error: Invalid --cluster parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        calculation_type |= CALC_CLUSTER;
      } else if (!memcmp(argptr2, "c", 2)) {
        logprint("Note: --cc flag deprecated.  Use '--cluster cc'.\n");
        cluster.modifier |= CLUSTER_CC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "luster-missing", 15)) {
	if (calculation_type & CALC_CLUSTER) {
	  sprintf(logbuf, "Error: --cluster-missing cannot be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --cluster-missing flag deprecated.  Use '--cluster missing'.\n");
        calculation_type |= CALC_CLUSTER;
        cluster.modifier |= CLUSTER_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "file", 5)) {
        UNSTABLE;
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	sptr = argv[cur_arg + 1];
	uii = strlen(sptr);
	if (uii > (FNAMESIZE - 9)) {
	  logprint("Error: --cfile parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(memcpya(pedname, sptr, uii), ".cnv", 5);
	memcpy(memcpya(famname, sptr, uii), ".fam", 5);
	memcpy(memcpya(mapname, sptr, uii), ".cnv.map", 9);
	load_rare = LOAD_RARE_CNV;
      } else if (!memcmp(argptr2, "nv-count", 9)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_COUNT;
      } else if (!memcmp(argptr2, "nv-del", 7)) {
	UNSTABLE;
	cnv_calc_type |= CNV_DEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-disrupt", 11)) {
	UNSTABLE;
	cnv_overlap_type = CNV_DISRUPT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-dup", 7)) {
	UNSTABLE;
	if (cnv_calc_type & CNV_DEL) {
	  sprintf(logbuf, "Error: --cnv-dup cannot be used with --cnv-del.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_calc_type |= CNV_DUP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-enrichment-test", 19)) {
	UNSTABLE;
	if (!cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-enrichment-test must be used with --cnv-count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --cnv-enrichment-test permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_enrichment_test_mperms = ii;
	}
	cnv_calc_type |= CNV_ENRICHMENT_TEST;
      } else if (!memcmp(argptr2, "nv-exclude", 11)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-exclude cannot be used with --cnv-count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_EXCLUDE;
      } else if (!memcmp(argptr2, "nv-exclude-off-by-1", 20)) {
	UNSTABLE;
        cnv_calc_type |= CNV_EXCLUDE_OFF_BY_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "nv-freq-exclude-above", 22)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-exclude-above parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_ABOVE;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-exclude-below", 22)) {
	UNSTABLE;
	if (cnv_freq_type) {
	  logprint("Error: --cnv-freq-exclude-below cannot be used with --cnv-freq-exclude-above.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-exclude-below parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_BELOW;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-exclude-exact", 22)) {
	UNSTABLE;
	if (cnv_freq_type) {
	  logprint("Error: --cnv-freq-exclude-exact cannot be used with\n--cnv-freq-exclude-above/-below.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-exclude-exact parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_EXCLUDE_EXACT;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-include-exact", 22)) {
	UNSTABLE;
	if (cnv_freq_type) {
	  logprint("Error: --cnv-freq-include-exact cannot be used with\n--cnv-freq-exclude-above/-below/-exact.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-freq-include-exact parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_freq_type = CNV_FREQ_INCLUDE_EXACT;
	cnv_freq_val = ii;
      } else if (!memcmp(argptr2, "nv-freq-method2", 16)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &cnv_freq_val2) || (cnv_freq_val2 < 0) || (cnv_freq_val2 > 1)) {
	    sprintf(logbuf, "Error: Invalid --cnv-freq-method2 parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	cnv_freq_type |= CNV_FREQ_METHOD2;
	if (cnv_freq_val2 == 0) {
	  // allow >= comparison to be used
	  cnv_freq_val2 = SMALLISH_EPSILON;
	}
      } else if (!memcmp(argptr2, "nv-freq-overlap", 16)) {
	UNSTABLE;
	if (!(cnv_freq_type & CNV_FREQ_FILTER)) {
	  logprint("Error: --cnv-freq-overlap must be used with --cnv-freq-include-exact or\n--cnv-freq-exclude-above/-below/-exact.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_freq_type & CNV_FREQ_METHOD2) {
	  logprint("Error: --cnv-freq-overlap cannot be used with --cnv-freq-method2.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &cnv_freq_val2) || (cnv_freq_val2 < 0) || (cnv_freq_val2 > 1)) {
	    sprintf(logbuf, "Error: Invalid --cnv-freq-overlap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (cnv_freq_val2 == 0) {
	  cnv_freq_val2 = SMALLISH_EPSILON;
	}
	cnv_freq_type |= CNV_FREQ_OVERLAP;
      } else if (!memcmp(argptr2, "nv-indiv-perm", 14)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --cnv-indiv-perm permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_indiv_mperms = ii;
	}
	cnv_calc_type |= CNV_INDIV_PERM;
      } else if (!memcmp(argptr2, "nv-intersect", 13)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-intersect cannot be used with --cnv-count/--cnv-exclude.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_intersect_filter_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	cnv_intersect_filter_type = CNV_INTERSECT;
      } else if (!memcmp(argptr2, "nv-kb", 6)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx > 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --cnv-kb size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_min_seglen = (int32_t)(dxx * 1000 + EPSILON);
      } else if (!memcmp(argptr2, "nv-list", 8)) {
	UNSTABLE;
	if ((load_rare & (~LOAD_RARE_CNV)) || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	strcpya(pedname, argv[cur_arg + 1]);
	load_rare = LOAD_RARE_CNV;
      } else if (!memcmp(argptr2, "nv-make-map", 12)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-make-map cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "short")) {
            sprintf(logbuf, "Error: Invalid --cnv-make-map parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_calc_type |= CNV_MAKE_MAP;
	} else {
	  cnv_calc_type |= CNV_MAKE_MAP | CNV_MAKE_MAP_LONG;
	}
      } else if (!memcmp(argptr2, "nv-max-kb", 10)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-max-kb cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx > 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --cnv-max-kb size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_max_seglen = (int32_t)(dxx * 1000 + EPSILON);
	if (cnv_min_seglen > cnv_max_seglen) {
	  logprint("Error: --cnv-max-kb value cannot be smaller than --cnv-kb value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-max-score", 13)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-max-score cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_max_score)) {
	  sprintf(logbuf, "Error: Invalid --cnv-max-score value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "nv-max-sites", 13)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-max-sites cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], (int32_t*)(&cnv_max_sites))) {
	  sprintf(logbuf, "Error: Invalid --cnv-max-sites parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "nv-overlap", 11)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type == CNV_DISRUPT) {
	  sprintf(logbuf, "Error: --cnv-overlap cannot be used with --cnv-disrupt.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val < 0) || (cnv_overlap_val > 1))  {
	  sprintf(logbuf, "Error: Invalid --cnv-overlap value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_overlap_val != 0) {
	  // ignore --cnv-overlap 0
	  cnv_overlap_type = CNV_OVERLAP;
	}
	if ((cnv_freq_type & CNV_FREQ_FILTER) && (!(cnv_freq_type & (CNV_FREQ_OVERLAP | CNV_FREQ_METHOD2)))) {
	  logprint("Note: --cnv-overlap + --cnv-freq-... deprecated.  Use --cnv-freq-overlap.\n");
	  if (cnv_overlap_val != 0) {
	    cnv_freq_type |= CNV_FREQ_OVERLAP;
	    cnv_freq_val2 = cnv_overlap_val;
	  }
	}
      } else if (!memcmp(argptr2, "nv-region-overlap", 18)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-region-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type) {
	  sprintf(logbuf, "Error: --cnv-region-overlap cannot be used with --cnv-overlap/-disrupt.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val <= 0) || (cnv_overlap_val > 1))  {
	  sprintf(logbuf, "Error: Invalid --cnv-region-overlap value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_overlap_type = CNV_OVERLAP_REGION;
      } else if (!memcmp(argptr2, "nv-score", 9)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-score cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_min_score)) {
	  sprintf(logbuf, "Error: Invalid --cnv-score value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_min_score > cnv_max_score) {
	  logprint("Error: --cnv-score value cannot be greater than --cnv-max-score value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-sites", 9)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-sites cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], (int32_t*)(&cnv_min_sites))) {
	  sprintf(logbuf, "Error: Invalid --cnv-sites parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (cnv_min_sites > cnv_max_sites) {
	  logprint("Error: --cnv-sites value cannot be greater than --cnv-max-sites value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
      } else if (!memcmp(argptr2, "nv-subset", 10)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-subset cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!cnv_intersect_filter_type) {
	  sprintf(logbuf, "Error: --cnv-subset must be used with --cnv-intersect/-exclude/-count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cnv_subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "nv-test", 8)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-test cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 2) {
	  if (!strcmp(argv[cur_arg + 1], "1sided")) {
	    uii = 2;
	    cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
	  } else if (!strcmp(argv[cur_arg + 1], "2sided")) {
	    uii = 2;
	    cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
	  } else {
	    uii = 1;
	    if (!strcmp(argv[cur_arg + 2], "1sided")) {
	      cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
	    } else if (!strcmp(argv[cur_arg + 2], "2sided")) {
	      cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
	    } else {
	      sprintf(logbuf, "Error: Invalid --cnv-test parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	} else {
	  uii = 1;
	}
	ii = atoi(argv[cur_arg + uii]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --cnv-test permutation count '%s'.%s", argv[cur_arg + uii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_test_mperms = ii;
	cnv_calc_type |= CNV_TEST;
      } else if (!memcmp(argptr2, "nv-test-1sided", 15)) {
	UNSTABLE;
	if (cnv_calc_type & CNV_TEST_FORCE_2SIDED) {
	  logprint("Error: --cnv-test cannot be both 1-sided and 2-sided at the same time.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --cnv-test-1sided flag deprecated.  Use '--cnv-test 1sided'.\n");
	cnv_calc_type |= CNV_TEST_FORCE_1SIDED;
      } else if (!memcmp(argptr2, "nv-test-2sided", 15)) {
	UNSTABLE;
	if (cnv_calc_type & CNV_TEST_FORCE_1SIDED) {
	  logprint("Error: --cnv-test cannot be both 1-sided and 2-sided at the same time.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --cnv-test-2sided flag deprecated.  Use '--cnv-test 2sided'.\n");
	cnv_calc_type |= CNV_TEST_FORCE_2SIDED;
      } else if (!memcmp(argptr2, "nv-test-region", 15)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-test-region cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --cnv-test-region permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_test_region_mperms = ii;
	}
	cnv_calc_type |= CNV_TEST_REGION;
      } else if (!memcmp(argptr2, "nv-test-window", 15)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-test-window cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001)) {
	  sprintf(logbuf, "Error: Invalid --cnv-test-window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dxx *= 1000;
	if (dxx > 2147483647) {
	  cnv_test_window = 0x7fffffff;
	} else {
	  cnv_test_window = (int32_t)(dxx + EPSILON);
	}
      } else if (!memcmp(argptr2, "nv-union-overlap", 17)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-union-overlap cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (cnv_overlap_type) {
	  sprintf(logbuf, "Error: --cnv-union-overlap cannot be used with --cnv-[region-]overlap/-disrupt.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &cnv_overlap_val) || (cnv_overlap_val <= 0) || (cnv_overlap_val > 1)) {
	  sprintf(logbuf, "Error: Invalid --cnv-union-overlap value '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cnv_overlap_type = CNV_OVERLAP_UNION;
      } else if (!memcmp(argptr2, "nv-write", 9)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-write cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "freq")) {
            sprintf(logbuf, "Error: Invalid --cnv-write parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (!(cnv_freq_val & CNV_FREQ_METHOD2)) {
	    sprintf(logbuf, "Error: --cnv-write 'freq' modifier must be used with --cnv-freq-method2.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  cnv_calc_type |= CNV_WRITE_FREQ;
	}
	cnv_calc_type |= CNV_WRITE;
      } else if (!memcmp(argptr2, "nv-write-freq", 14)) {
	UNSTABLE;
	if (!(load_rare & LOAD_RARE_CNV)) {
	  logprint("Error: --cnv-write freq cannot be used without a .cnv fileset.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (!(cnv_freq_val & CNV_FREQ_METHOD2)) {
	  sprintf(logbuf, "Error: --cnv-write 'freq' modifier must be used with --cnv-freq-method2.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (!(cnv_calc_type & CNV_WRITE)) {
	  sprintf(logbuf, "Error: --cnv-write-freq must be used with --cnv-write.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --cnv-write-freq flag deprecated.  Use '--cnv-write freq'.\n");
	cnv_calc_type |= CNV_WRITE_FREQ;
      } else if (!memcmp(argptr2, "onsensus-match", 15)) {
        logprint("Note: --consensus-match flag deprecated.  Use '--homozyg consensus-match'.\n");
	homozyg.modifier |= HOMOZYG_CONSENSUS_MATCH;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ondition", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if (!strcmp("dominant", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_DOMINANT;
	  } else if (!strcmp("recessive", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_RECESSIVE;
	  } else {
	    uii = 2;
            if (!strcmp("dominant", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_DOMINANT;
	    } else if (!strcmp("recessive", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_RECESSIVE;
	    } else {
	      sprintf(logbuf, "Error: Invalid --condition parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	if (alloc_string(&condition_mname, argv[cur_arg + uii])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ondition-list", 14)) {
	if (condition_mname) {
	  logprint("Error: --condition-list cannot be used with --condition.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 2) {
	  if (!strcmp("dominant", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_DOMINANT;
	  } else if (!strcmp("recessive", argv[cur_arg + 2])) {
	    glm_modifier |= GLM_CONDITION_RECESSIVE;
	  } else {
	    uii = 2;
            if (!strcmp("dominant", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_DOMINANT;
	    } else if (!strcmp("recessive", argv[cur_arg + 1])) {
	      glm_modifier |= GLM_CONDITION_RECESSIVE;
	    } else {
	      sprintf(logbuf, "Error: Invalid --condition-list parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  }
	}
	retval = alloc_fname(&condition_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'd':
      if (!memcmp(argptr2, "ata", 4)) {
	if (load_params & 0xff) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x80;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 5)) {
	    logprint("Error: --data parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	} else {
	  sptr = (char*)PROG_NAME_STR;
	}
	if (!(load_params & 0x100)) {
	  memcpy(strcpya(genname, sptr), ".gen", 5);
	}
	if (!(load_params & 0x200)) {
	  memcpy(strcpya(samplename, sptr), ".sample", 8);
	}
      } else if (!memcmp(argptr2, "ebug", 5)) {
	debug_on = 1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ecompress", 10)) {
	logprint("Error: --decompress flag retired.  Use e.g. 'gunzip [filename]'.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "istance", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 7)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (dist_calc_type & DISTANCE_BIN) {
	      sprintf(logbuf, "Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (dist_calc_type & DISTANCE_GZ) {
	      sprintf(logbuf, "Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "ibs")) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logprint("Error: Duplicate --distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "1-ibs")) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logprint("Error: Duplicate --distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "allele-ct")) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logprint("Error: Duplicate --distance 'allele-ct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!strcmp(argv[cur_arg + uii], "3d")) {
	    dist_calc_type |= DISTANCE_3D;
	  } else if (!strcmp(argv[cur_arg + uii], "flat-missing")) {
	    dist_calc_type |= DISTANCE_FLAT_MISSING;
	  } else {
	    sprintf(logbuf, "Error: Invalid --distance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(dist_calc_type & DISTANCE_TYPEMASK)) {
	  dist_calc_type |= DISTANCE_ALCT;
	}
	calculation_type |= CALC_DISTANCE;
      } else if (!memcmp(argptr2, "istance-matrix", 15)) {
	if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	  sprintf(logbuf, "Error: --distance-matrix flag cannot be used with '--distance 1-ibs'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_PLINK_DISTANCE_MATRIX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ummy", 5)) {
	if (load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 6)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy individual count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_indiv_ct = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --dummy marker count.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dummy_marker_ct = ii;
        for (uii = 3; uii <= param_ct; uii++) {
	  if (match_upper(argv[cur_arg + uii], "ACGT")) {
	    if (dummy_flags & (DUMMY_1234 | DUMMY_12)) {
	      sprintf(logbuf, "Error: --dummy 'acgt' modifier cannot be used with '1234' or '12'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_ACGT;
	  } else if (!strcmp(argv[cur_arg + uii], "1234")) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_12)) {
	      sprintf(logbuf, "Error: --dummy '1234' modifier cannot be used with 'acgt' or '12'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_1234;
	  } else if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (dummy_flags & (DUMMY_ACGT | DUMMY_1234)) {
	      sprintf(logbuf, "Error: --dummy '12' modifier cannot be used with 'acgt' or '1234'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            dummy_flags |= DUMMY_12;
	  } else if (!strcmp(argv[cur_arg + uii], "scalar-pheno")) {
	    dummy_flags |= DUMMY_SCALAR_PHENO;
	  } else {
	    if ((dummy_flags & DUMMY_MISSING_PHENO) || scan_double(argv[cur_arg + uii], &dxx) || (dxx < 0.0) || (dxx > 1.0)) {
	      sprintf(logbuf, "Error: Invalid --dummy parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (dummy_flags & DUMMY_MISSING_GENO) {
	      dummy_missing_pheno = dxx;
	      dummy_flags |= DUMMY_MISSING_PHENO;
	    } else {
	      dummy_missing_geno = dxx;
	      dummy_flags |= DUMMY_MISSING_GENO;
	    }
	  }
	}
	load_rare = LOAD_RARE_DUMMY;
      } else if (!memcmp(argptr2, "ummy-coding", 12)) {
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --dummy-coding cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "no-round")) {
	    uii = 2;
	    write_covar_modifier |= WRITE_COVAR_DUMMY_NO_ROUND;
	  } else {
	    if (param_ct == 2) {
              if (strcmp(argv[cur_arg + 2], "no-round")) {
		sprintf(logbuf, "Error: Invalid --dummy-coding parameter sequence.%s", errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      write_covar_modifier |= WRITE_COVAR_DUMMY_NO_ROUND;
	    }
	    uii = 1;
	  }
	  if (uii <= param_ct) {
            ii = atoi(argv[cur_arg + uii]);
	    if (ii < 3) {
	      sprintf(logbuf, "Error: Invalid --dummy-coding max categories parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            write_covar_dummy_max_categories = ii;
	  }
	}
	write_covar_modifier |= WRITE_COVAR_DUMMY;
      } else if (!memcmp(argptr2, "ominant", 8)) {
	logprint("Note: --dominant flag deprecated.  Use e.g. '--linear dominant' (and\n'--condition-list [filename] dominant' to change covariate coding).\n");
	glm_modifier |= GLM_DOMINANT | GLM_CONDITION_DOMINANT;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'e':
      if (!memcmp(argptr2, "xtract", 7)) {
	if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "--extract cannot be used with a .cnv fileset.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&extractname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "xclude", 7)) {
	if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "--exclude cannot be used with a .cnv fileset.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&excludename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "xclude-snp", 11)) {
	if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --exclude-snp cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb\nand --to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
        misc_flags |= MISC_EXCLUDE_MARKERNAME_SNP;
      } else if (!memcmp(argptr2, "xclude-snps", 12)) {
	if (markername_snp) {
	  sprintf(logbuf, "Error: --exclude-snps cannot be used with --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --exclude-snps cannot be used with a .cnv fileset.  Use\n--from-bp/-kb/-mb and --to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &snps_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
        misc_flags |= MISC_EXCLUDE_MARKERNAME_SNP;
      } else if (!memcmp(argptr2, "xclude-before-extract", 22)) {
        logprint("Note: --exclude-before-extract has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "xponent", 8)) {
	if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
	  logprint("Error: --exponent cannot be used with --distance-matrix.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &exponent)) {
	  sprintf(logbuf, "Error: Invalid --exponent parameter '%s'.\n", argv[cur_arg + 1]);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'f':
      if (!memcmp(argptr2, "ile", 4)) {
	if (load_params & 0x3f9) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 1;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 5)) {
	    logprint("Error: --file parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  if (!(load_params & 2)) {
	    memcpy(strcpya(pedname, argv[cur_arg + 1]), ".ped", 5);
	  }
	  if (!(load_params & 4)) {
	    memcpy(strcpya(mapname, argv[cur_arg + 1]), ".map", 5);
	  }
	}
      } else if (!memcmp(argptr2, "am", 3)) {
	if (load_params & 0x3c7) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 64;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --fam parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(famname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "ilter", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&filtername, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        retval = alloc_and_flatten(&filtervals_flattened, &(argv[cur_arg + 2]), param_ct - 1);
	if (retval) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "ilter-cases", 12)) {
	filter_binary |= FILTER_BINARY_CASES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-controls", 15)) {
	if (filter_binary & FILTER_BINARY_CASES) {
	  logprint("Error: --filter-cases and --filter-controls cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_binary |= FILTER_BINARY_CONTROLS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-females", 14)) {
	filter_binary |= FILTER_BINARY_FEMALES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-males", 12)) {
	if (filter_binary & FILTER_BINARY_FEMALES) {
	  logprint("Error: --filter-males and --filter-females cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_binary |= FILTER_BINARY_MALES;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-founders", 15)) {
	filter_binary |= FILTER_BINARY_FOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ilter-nonfounders", 18)) {
	if (filter_binary & FILTER_BINARY_FOUNDERS) {
	  logprint("Error: --filter-founders and --filter-nonfounders cannot be used together.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	filter_binary |= FILTER_BINARY_NONFOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "req", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "counts")) {
            sprintf(logbuf, "Error: Invalid --freq parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  misc_flags |= MISC_FREQ_COUNTS;
	}
	calculation_type |= CALC_FREQ;
	if (misc_flags & MISC_FREQ_COUNTS) {
	  // --keep-allele-order also set for backward compatibility
	  misc_flags |= MISC_KEEP_ALLELE_ORDER;
	}
      } else if (!memcmp(argptr2, "reqx", 5)) {
	if (calculation_type & CALC_FREQ) {
	  sprintf(logbuf, "Error: --freqx cannot be used with --freq.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_FREQ;
	misc_flags |= MISC_FREQX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rom", 4)) {
	if (chrom_flag_present) {
	  sprintf(logbuf, "Error: --from cannot be used with --autosome[-xy] or --[not-]chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --from cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb\ninstead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_from, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if ((!memcmp(argptr2, "rom-bp", 7)) || (!memcmp(argptr2, "rom-kb", 7)) || (!memcmp(argptr2, "rom-mb", 7))) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --from-bp/-kb/-mb cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argptr2[4];
	if (cc == 'b') {
	  if (atoiz(argv[cur_arg + 1], &marker_pos_start)) {
	    sprintf(logbuf, "Error: Invalid --from-bp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  if (marker_pos_start != -1) {
	    logprint("Error: Multiple --from-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (scan_double(argv[cur_arg + 1], &dxx)) {
	    sprintf(logbuf, "Error: Invalid --from-kb/-mb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    marker_pos_start = 0;
	  } else if (dxx > 2147483647) {
	    marker_pos_start = 0x7fffffff;
	  } else {
	    marker_pos_start = (int32_t)(dxx + EPSILON);
	  }
	}
      } else if (!memcmp(argptr2, "isher", 6)) {
	if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --fisher cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --fisher flag deprecated.  Use '--assoc fisher' or '--model fisher'.\n");
	model_modifier |= MODEL_ASSOC | MODEL_FISHER | MODEL_ASSOC_FDEPR;
	calculation_type |= CALC_MODEL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "id", 3)) {
        logprint("Note: --fid flag deprecated.  Use '--recode vcf-fid'.\n");
	recode_modifier |= RECODE_FID;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lip", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&flip_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "lip-subset", 11)) {
	UNSTABLE;
	if (!flip_fname) {
          sprintf(logbuf, "Error: --flip-subset must be used with --flip.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (allelexxxx) {
	  // fix for this is too messy to be worthwhile
	  logprint("Error: --flip-subset cannot be used with --allele1234/ACGT.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&flip_subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'g':
      if (!memcmp(argptr2, "en", 3)) {
	if (load_params & 0x17f) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x100;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --gen parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(genname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "eno", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &geno_thresh)) {
	    sprintf(logbuf, "Error: Invalid --geno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((geno_thresh < 0.0) || (geno_thresh > 1.0)) {
	    sprintf(logbuf, "Error: Invalid --geno parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  geno_thresh = 0.1;
	}
      } else if ((!memcmp(argptr2, "enome", 6)) || (!memcmp(argptr2, "enome gz", 9))) {
	if (argptr2[5] == ' ') {
	  kk = 1;
	  genome_modifier |= GENOME_OUTPUT_GZ;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5 - kk)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "gz")) {
            genome_modifier |= GENOME_OUTPUT_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "rel-check")) {
            genome_modifier |= GENOME_REL_CHECK;
	  } else if (!strcmp(argv[cur_arg + uii], "full")) {
	    genome_modifier |= GENOME_OUTPUT_FULL;
	  } else if (!strcmp(argv[cur_arg + uii], "unbounded")) {
	    genome_modifier |= GENOME_IBD_UNBOUNDED;
	  } else if (!strcmp(argv[cur_arg + uii], "nudge")) {
            genome_modifier |= GENOME_NUDGE;
	  } else {
	    sprintf(logbuf, "Error: Invalid --genome parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_GENOME;
      } else if (!memcmp(argptr2, "enome-full", 11)) {
	if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --genome-full must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --genome-full flag deprecated.  Use '--genome full'.\n");
	genome_modifier |= GENOME_OUTPUT_FULL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "roupdist", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  groupdist_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((groupdist_iters < 2) || (groupdist_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --groupdist jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    ii = atoi(argv[cur_arg + 2]);
	    if (ii <= 0) {
	      sprintf(logbuf, "Error: Invalid --groupdist jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    groupdist_d = ii;
	  }
	}
	calculation_type |= CALC_GROUPDIST;
      } else if (!memcmp(argptr2, "rm", 3)) {
	logprint("Error: --grm has been retired due to inconsistent meaning across GCTA versions.\nUse --grm-gz or --grm-bin.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "rm-gz", 6)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 8)) {
	    logprint("Error: --grm-gz parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
        load_rare = LOAD_RARE_GRM;
      } else if (!memcmp(argptr2, "rm-bin", 7)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (param_ct) {
	  sptr = argv[cur_arg + 1];
	  if (strlen(sptr) > (FNAMESIZE - 11)) {
	    logprint("Error: --grm-bin parameter too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, sptr);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
        load_rare = LOAD_RARE_GRM_BIN;
      } else if (!memcmp(argptr2, "xe", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  logprint("Error: --gxe must be used with --covar.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --gxe parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  gxe_mcovar = ii;
	} else {
	  gxe_mcovar = 1;
	}
	calculation_type |= CALC_GXE;
      } else if (!memcmp(argptr2, "enedrop", 8)) {
	if (model_modifier & MODEL_QMASK) {
	  sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with --genedrop.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --genedrop flag deprecated.  Use e.g. '--model genedrop'.\n");
	model_modifier |= MODEL_GENEDROP;
	glm_modifier |= GLM_GENEDROP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "c", 2)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --gc must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --gc flag deprecated.  Use '--adjust gc'.\n");
	mtest_adjust |= ADJUST_GC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "file", 5)) {
	UNSTABLE;
	if (load_rare || (load_params & (~64))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	sptr = argv[cur_arg + 1];
	uii = strlen(sptr);
	if (uii > (FNAMESIZE - 6)) {
	  logprint("Error: --gfile parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(memcpya(pedname, sptr, uii), ".gvar", 6);
	if (!(load_params & 64)) {
	  memcpy(memcpya(famname, sptr, uii), ".fam", 5);
	}
	memcpy(memcpya(mapname, sptr, uii), ".map", 5);
	load_rare = LOAD_RARE_GVAR;
      } else if (!memcmp(argptr2, "enome-lists", 12)) {
	logprint("Error: --genome-lists flag retired.  Use --parallel.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "enome-minimal", 14)) {
	logprint("Error: --genome-minimal flag retired.  Use '--genome gz'.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if ((!memcmp(argptr2, "roup-avg", 9)) || (!memcmp(argptr2, "roup-average", 13))) {
        if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --group-avg must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (cluster.modifier & CLUSTER_OLD_TIEBREAKS) {
	  sprintf(logbuf, "Error: --cluster 'group-avg' and 'old-tiebreaks' cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        sprintf(logbuf, "Note: --%s flag deprecated.  Use '--cluster group-avg'.\n", argptr);
	logprintb();
	cluster.modifier |= CLUSTER_GROUP_AVG;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "enotypic", 9)) {
	if (glm_modifier & GLM_DOMINANT) {
	  logprint("Error: --genotypic cannot be used with --dominant.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --genotypic flag deprecated.  Use e.g. '--linear genotypic'.\n");
	glm_modifier |= GLM_GENOTYPIC;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'h':
      if (!memcmp(argptr2, "we", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &hwe_thresh)) {
	    sprintf(logbuf, "Error: Invalid --hwe parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((hwe_thresh < 0.0) || (hwe_thresh >= 1.0)) {
	    sprintf(logbuf, "Error: Invalid --hwe parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  hwe_thresh = 0.001;
	}
      } else if (!memcmp(argptr2, "we-all", 7)) {
	misc_flags |= MISC_HWE_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et", 3)) {
	sprintf(logbuf, "Error: --het retired.  Use --ibc.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "ardy", 5)) {
	calculation_type |= CALC_HARDY;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "we2", 4)) {
	sprintf(logbuf, "Error: --hwe2 retired.  Use the --hwe exact test.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "ardy2", 6)) {
	sprintf(logbuf, "Error: --hardy2 retired.  Use the exact test-based --hardy report.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if (!memcmp(argptr2, "omozyg", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "group")) {
	    if (homozyg.modifier & HOMOZYG_GROUP_VERBOSE) {
	      logprint("Error: --homozyg 'group' and 'group-verbose' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    homozyg.modifier |= HOMOZYG_GROUP;
	  } else if (!strcmp(argv[cur_arg + uii], "group-verbose")) {
	    if (homozyg.modifier & HOMOZYG_GROUP) {
	      logprint("Error: --homozyg 'group' and 'group-verbose' modifiers cannot be used together.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    homozyg.modifier |= HOMOZYG_GROUP_VERBOSE;
	  } else if (!strcmp(argv[cur_arg + uii], "consensus-match")) {
	    homozyg.modifier |= HOMOZYG_CONSENSUS_MATCH;
	  } else if (!strcmp(argv[cur_arg + uii], "extend")) {
	    homozyg.modifier |= HOMOZYG_EXTEND;
	  } else if (!strcmp(argv[cur_arg + uii], "subtract-1-from-lengths")) {
            homozyg.modifier |= HOMOZYG_OLD_LENGTHS;
	  } else {
	    sprintf(logbuf, "Error: Invalid --homozyg parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_HOMOZYG;
      } else if (!memcmp(argptr2, "omozyg-snp", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --homozyg-snp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_HOMOZYG;
	homozyg.min_snp = ii;
      } else if (!memcmp(argptr2, "omozyg-kb", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < EPSILON) || (dxx >= (2147483.647 + EPSILON))) {
	  sprintf(logbuf, "Error: Invalid --homozyg-kb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_HOMOZYG;
	// round up
	homozyg.min_bases = 1 + (uint32_t)((int32_t)(dxx * 1000 - EPSILON));
      } else if (!memcmp(argptr2, "omozyg-density", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx >= 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-density parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.max_bases_per_snp = dxx * 1000 + EPSILON;
      } else if (!memcmp(argptr2, "omozyg-gap", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx) || (dxx < 0.001) || (dxx >= 2147483.647)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.max_gap = ((int32_t)(dxx * 1000 + EPSILON));
      } else if (!memcmp(argptr2, "omozyg-het", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-het parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ii && (homozyg.modifier & HOMOZYG_EXTEND)) {
	  sprintf(logbuf, "Error: --homozyg-het with a nonzero parameter cannot be used with --homozyg\nextend.  For fine-grained control over heterozygote frequency, use\n--homozyg-window-snp and --homozyg-window-het instead.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_HOMOZYG;
        homozyg.max_hets = ii;
      } else if (!memcmp(argptr2, "omozyg-window-snp", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-snp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.window_size = ii;
      } else if (!memcmp(argptr2, "omozyg-window-kb", 17)) {
        logprint("Error: --homozyg-window-kb flag provisionally retired, since it had no effect\nin PLINK 1.07.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "omozyg-window-het", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-het parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.window_max_hets = ii;
      } else if (!memcmp(argptr2, "omozyg-window-missing", 22)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-missing parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.window_max_missing = ii;
      } else if (!memcmp(argptr2, "omozyg-window-threshold", 24)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-window-threshold parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_HOMOZYG;
	homozyg.hit_threshold = dxx;
      } else if (!memcmp(argptr2, "omozyg-match", 13)) {
	if (!(homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE))) {
          homozyg.modifier |= HOMOZYG_GROUP;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx) || (dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(logbuf, "Error: Invalid --homozyg-match parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	homozyg.overlap_min = dxx;
      } else if (!memcmp(argptr2, "omozyg-group", 13)) {
	if (homozyg.modifier & HOMOZYG_GROUP_VERBOSE) {
	  logprint("Note: --homozyg-group deprecated, and superseded by --homozyg group-verbose.\n");
	} else {
	  logprint("Note: --homozyg-group flag deprecated.  Use '--homozyg group'.\n");
	  homozyg.modifier |= HOMOZYG_GROUP;
	}
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omozyg-verbose", 15)) {
	if (!(homozyg.modifier & HOMOZYG_GROUP)) {
	  logprint("Error: --homozyg-verbose must be used with --homozyg group.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --homozyg-verbose flag deprecated.  Use '--homozyg group-verbose'.\n");
	homozyg.modifier = (homozyg.modifier & (~HOMOZYG_GROUP)) | HOMOZYG_GROUP_VERBOSE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omozyg-include-missing", 23)) {
        logprint("Error: --homozyg-include-missing flag provisionally retired, since it had no\neffect in PLINK 1.07.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ethom", 6)) {
	if (!(glm_modifier & GLM_GENOTYPIC)) {
	  logprint("Error: --hethom must be used with --genotypic.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --hethom flag deprecated.  Use e.g. '--linear hethom' (and\n'--condition-list [filename] recessive' to change covariate coding).\n");
	glm_modifier |= GLM_HETHOM | GLM_CONDITION_RECESSIVE;
	glm_modifier -= GLM_GENOTYPIC;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ide-covar", 10)) {
	logprint("Note: --hide-covar flag deprecated.  Use e.g. '--linear hide-covar'.\n");
	glm_modifier |= GLM_HIDE_COVAR;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "omog", 5)) {
        calculation_type |= CALC_HOMOG;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'i':
      if (!memcmp(argptr2, "bc", 3)) {
	calculation_type |= CALC_IBC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ndep-pairwise", 14)) {
	if (calculation_type & CALC_LD_PRUNE) {
	  logprint("Error: --indep-pairwise cannot be used with --indep.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 3, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || ((ii == 1) && (param_ct == 3))) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_window_size = ii;
	if (param_ct == 4) {
	  if (!match_upper(argv[cur_arg + 2], "KB")) {
	    sprintf(logbuf, "Error: Invalid --indep-pairwise parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ld_window_kb = 1;
	} else {
	  jj = strlen(argv[cur_arg + 1]);
	  if ((jj > 2) && match_upper(&(argv[cur_arg + 1][jj - 2]), "KB")) {
	    ld_window_kb = 1;
	  }
	}
	ii = atoi(argv[cur_arg + param_ct - 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid increment '%s' for --indep-pairwise.%s", argv[cur_arg + param_ct - 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_window_incr = ii;
	if (scan_double(argv[cur_arg + param_ct], &ld_last_param) || (ld_last_param < 0.0) || (ld_last_param >= 1.0)) {
	  sprintf(logbuf, "Error: Invalid --indep-pairwise r^2 threshold '%s'.%s", argv[cur_arg + param_ct], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_LD_PRUNE;
	misc_flags |= MISC_LD_PRUNE_PAIRWISE;
      } else if (!memcmp(argptr2, "ndep", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 3, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || ((ii == 1) && (param_ct == 3))) {
	  sprintf(logbuf, "Error: Invalid --indep window size '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_window_size = ii;
	if (param_ct == 4) {
	  if (!match_upper(argv[cur_arg + 2], "KB")) {
	    sprintf(logbuf, "Error: Invalid --indep parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ld_window_kb = 1;
	} else {
	  jj = strlen(argv[cur_arg + 1]);
	  if ((jj > 2) && match_upper(&(argv[cur_arg + 1][jj - 2]), "KB")) {
	    ld_window_kb = 1;
	  }
	}
	ii = atoi(argv[cur_arg + param_ct - 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid increment '%s' for --indep.%s", argv[cur_arg + param_ct - 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ld_window_incr = ii;
	if (scan_double(argv[cur_arg + param_ct], &ld_last_param)) {
	  sprintf(logbuf, "Error: Invalid --indep VIF threshold '%s'.%s", argv[cur_arg + param_ct], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ld_last_param < 1.0) {
	  sprintf(logbuf, "Error: --indep VIF threshold '%s' too small (must be >= 1).%s", argv[cur_arg + param_ct], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_LD_PRUNE;
      } else if (!memcmp(argptr2, "ndiv-sort", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = (argv[cur_arg + 1][1] == '\0');
        if ((!strcmp(argv[cur_arg + 1], "none")) || ((argv[cur_arg + 1][0] == '0') && jj)) {
	  indiv_sort = INDIV_SORT_NONE;
	} else if ((!strcmp(argv[cur_arg + 1], "natural")) || ((tolower(argv[cur_arg + 1][0]) == 'n') && jj)) {
	  indiv_sort = INDIV_SORT_NATURAL;
	} else if ((!strcmp(argv[cur_arg + 1], "ascii")) || ((tolower(argv[cur_arg + 1][0]) == 'a') && jj)) {
	  indiv_sort = INDIV_SORT_ASCII;
	} else {
	  sprintf(logbuf, "Error: '%s' is not a valid mode for --indiv-sort.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "bs-test", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  ii = atoi(argv[cur_arg + 1]);
	  if (ii < 1) {
	    sprintf(logbuf, "Error: Invalid --ibs-test permutation count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (ii < MAX_THREADS * 2) {
	    sprintf(logbuf, "Error: --ibs test permutation count '%s' too small (min %u).%s", argv[cur_arg + 1], MAX_THREADS * 2, errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  ibs_test_perms = ii;
	}
	calculation_type |= CALC_IBS_TEST;
      } else if (!memcmp(argptr2, "id", 3)) {
        logprint("Note: --iid flag deprecated.  Use '--recode vcf-iid'.\n");
	recode_modifier |= RECODE_IID;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "bm", 3)) {
        if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --ibm must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ibm parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx > 1.0)) {
	  sprintf(logbuf, "Error: --ibm threshold must be in (0, 1].%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.min_ibm = dxx;
      } else if (!memcmp(argptr2, "mpossible", 10)) {
	logprint("Error: --impossible flag retired.  Use '--genome nudge', or explicitly validate\nZ0/Z1/Z2/PI_HAT in your script.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "nteraction", 11)) {
	logprint("Note: --interaction flag deprecated.  Use e.g. '--linear interaction'.\n");
	glm_modifier |= GLM_INTERACTION;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'k':
      if (!memcmp(argptr2, "eep", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&keepname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-fam", 8)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&keepfamname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-allele-order", 17)) {
	misc_flags |= MISC_KEEP_ALLELE_ORDER;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-before-remove", 18)) {
        logprint("Note: --keep-before-remove has no effect.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-autoconv", 13)) {
        misc_flags |= MISC_KEEP_AUTOCONV;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "eep-clusters", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&(cluster.keep_fname), argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eep-cluster-names", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_and_flatten(&(cluster.keep_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'l':
      if (!memcmp(argptr2, "file", 5)) {
	if (load_rare || load_params) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 6) {
	    logprint("Error: --lfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  strcpy(pedname, argv[cur_arg + 1]);
	} else {
	  memcpy(pedname, PROG_NAME_STR, 6);
	}
	load_rare = LOAD_RARE_LGEN;
      } else if (!memcmp(argptr2, "oop-assoc", 10)) {
	if (pheno_modifier & PHENO_ALL) {
	  sprintf(logbuf, "Error: --loop-assoc cannot be used with --all-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if ((strlen(argv[cur_arg + 1]) == 7) && (!memcmp(argv[cur_arg + 1], "keep-", 5)) && match_upper(&(argv[cur_arg + 1][5]), "NA")) {
	    uii = 2;
	  } else if ((strlen(argv[cur_arg + 2]) != 7) || memcmp(argv[cur_arg + 2], "keep-", 5) || (!match_upper(&(argv[cur_arg + 2][5]), "NA"))) {
            sprintf(logbuf, "Error: Invalid --loop-assoc parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_LOAD_CLUSTER_KEEP_NA;
	}
	retval = alloc_fname(&loop_assoc_fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "og10", 5)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --log10 must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --log10 flag deprecated.  Use '--adjust log10'.\n");
	mtest_adjust |= ADJUST_LOG10;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ambda", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --lambda must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &adjust_lambda)) {
	  sprintf(logbuf, "Error: Invalid --lambda parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (adjust_lambda < 1) {
	  logprint("Note: --lambda parameter set to 1.\n");
	  adjust_lambda = 1;
	}
	mtest_adjust |= ADJUST_LAMBDA;
      } else if (!memcmp(argptr2, "ist-23-indels", 14)) {
        calculation_type |= CALC_LIST_23_INDELS;
      } else if ((!memcmp(argptr2, "inear", 6)) || (!memcmp(argptr2, "ogistic", 8))) {
#ifndef NOLAPACK
        if (calculation_type & CALC_GLM) {
	  logprint("Error: --logistic cannot be used with --linear.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
	if (*argptr2 == 'o') {
	  glm_modifier |= GLM_LOGISTIC;
#ifdef NOLAPACK
	} else {
	  logprint("Error: --linear requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	  goto main_ret_INVALID_CMDLINE;
#endif
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 9)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (glm_modifier & GLM_MPERM) {
	      sprintf(logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            glm_modifier |= GLM_PERM;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
            if (glm_modifier & GLM_PERM) {
	      sprintf(logbuf, "Error: --%s 'mperm' and 'perm' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (glm_modifier & GLM_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --%s 'mperm' modifier.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ii = atoi(&(argv[cur_arg + uii][6]));
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --%s mperm parameter '%s'.%s", argptr, &(argv[cur_arg + uii][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            glm_mperm_val = (uint32_t)ii;
            glm_modifier |= GLM_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    glm_modifier |= GLM_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    glm_modifier |= GLM_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "genotypic")) {
	    if (glm_modifier & (GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_GENOTYPIC;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "hethom")) {
	    if (glm_modifier & (GLM_GENOTYPIC | GLM_DOMINANT | GLM_RECESSIVE)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_HETHOM;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "dominant")) {
	    if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_RECESSIVE)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_DOMINANT;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "recessive")) {
	    if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT)) {
	      sprintf(logbuf, "Error: Conflicting --%s parameters.\n", argptr);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_RECESSIVE;
	    glm_xchr_model = 0;
	  } else if (!strcmp(argv[cur_arg + uii], "no-snp")) {
	    if (mtest_adjust) {
	      sprintf(logbuf, "Error: --%s no-snp cannot be used with --adjust.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    // defer the rest of the check
	    glm_modifier |= GLM_NO_SNP;
	  } else if (!strcmp(argv[cur_arg + uii], "hide-covar")) {
	    glm_modifier |= GLM_HIDE_COVAR;
	  } else if (!strcmp(argv[cur_arg + uii], "sex")) {
	    if (glm_modifier & GLM_NO_X_SEX) {
	      sprintf(logbuf, "Error: --%s 'sex' and 'no-x-sex' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "no-x-sex")) {
	    if (glm_modifier & GLM_SEX) {
	      sprintf(logbuf, "Error: --%s 'sex' and 'no-x-sex' cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_NO_X_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "interaction")) {
	    glm_modifier |= GLM_INTERACTION;
	  } else if (!strcmp(argv[cur_arg + uii], "standard-beta")) {
	    if (glm_modifier & GLM_LOGISTIC) {
	      sprintf(logbuf, "Error: --logistic does not have a 'standard-beta' modifier.  (Did you mean\n--linear or 'beta'?)%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    glm_modifier |= GLM_STANDARD_BETA;
	  } else if (!strcmp(argv[cur_arg + uii], "beta")) {
	    glm_modifier |= GLM_BETA;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    sprintf(logbuf, "Error: Improper --%s mperm syntax.  (Use '--%s mperm=[value]'.)\n", argptr, argptr);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", argptr, argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if ((glm_modifier & GLM_NO_SNP) && (glm_modifier & GLM_NO_SNP_EXCL)) {
	  sprintf(logbuf, "Error: --%s 'no-snp' modifier conflicts with another modifier.%s", argptr, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_GLM;
      } else if (!memcmp(argptr2, "d-xchr", 7)) {
        if (!(calculation_type & CALC_LD_PRUNE)) {
          sprintf(logbuf, "Error: --ld-xchr must be used with --indep[-pairwise].%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '3') || (argv[cur_arg + 1][1] != '\0')) {
	  sprintf(logbuf, "Error: Invalid --ld-xchr parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (cc == '2') {
          misc_flags |= MISC_LD_IGNORE_X;
	} else if (cc == '3') {
	  misc_flags |= MISC_LD_WEIGHTED_X;
	}
      } else if (!memcmp(argptr2, "asso", 5)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
          goto main_ret_INVALID_CMDLINE_3;
	}
        uii = 1;
	ujj = 0;
	if (param_ct > 1) {
          if (!strcmp(argv[cur_arg + 1], "report-zeroes")) {
	    uii = 2;
            misc_flags |= MISC_LASSO_REPORT_ZEROES;
	  } else if (!strcmp(argv[cur_arg + 2], "report-zeroes")) {
	    if (param_ct == 3) {
	      ujj = 3;
	    }
            misc_flags |= MISC_LASSO_REPORT_ZEROES;
	  } else {
            ujj = 2;
	    if (param_ct == 3) {
	      if (strcmp(argv[cur_arg + 3], "report-zeroes")) {
	        sprintf(logbuf, "Error: Invalid --lasso parameter sequence.%s", errstr_append);
	        goto main_ret_INVALID_CMDLINE_3;
	      }
	      misc_flags |= MISC_LASSO_REPORT_ZEROES;
	    }
	  }
	}
	if (scan_double(argv[cur_arg + uii], &lasso_h2) || (lasso_h2 > 1) || (lasso_h2 <= 0)) {
	  sprintf(logbuf, "Error: Invalid --lasso heritability estimate '%s'.%s", argv[cur_arg + uii], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (ujj) {
	  if (scan_double(argv[cur_arg + ujj], &lasso_minlambda) || (lasso_minlambda <= 0)) {
	    sprintf(logbuf, "Error: Invalid --lasso minimum lambda '%s'.%s", argv[cur_arg + ujj], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        calculation_type |= CALC_LASSO;
      } else {
        goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'm':
      if (!memcmp(argptr2, "ap", 3)) {
	if ((load_params & 0x3fc) || (load_rare & (~(LOAD_RARE_CNV | LOAD_RARE_GVAR)))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 4;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --map parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(mapname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "issing-genotype", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((argv[cur_arg + 1][1] != '\0') || (((unsigned char)cc) <= ' ') || ((cc > '0') && (cc <= '4')) || (cc == 'A') || (cc == 'C') || (cc == 'G') || (cc == 'T')) {
	  sprintf(logbuf, "Error: Invalid --missing-genotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	g_missing_geno_ptr = &(g_one_char_strs[((unsigned char)cc) * 2]);
	g_output_missing_geno_ptr = g_missing_geno_ptr;
      } else if (!memcmp(argptr2, "issing-phenotype", 17)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	missing_pheno = atoi(argv[cur_arg + 1]);
	if ((missing_pheno == 0) || (missing_pheno == 1)) {
	  sprintf(logbuf, "Error: Invalid --missing-phenotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if ((!memcmp(argptr2, "issing-code", 12))) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  missing_code = argv[cur_arg + 1];
	} else {
	  missing_code = (char*)"";
	}
      } else if (!memcmp(argptr2, "ake-pheno", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (alloc_string(&makepheno_str, argv[cur_arg + 2])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "pheno", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mpheno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mpheno_col = ii;
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mfilter parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mfilter_col = ii;
      } else if (!memcmp(argptr2, "emory", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	malloc_size_mb = atoi(argv[cur_arg + 1]);
	if (malloc_size_mb < WKSPACE_MIN_MB) {
	  if (malloc_size_mb > 0) {
	    sprintf(logbuf, "Error: Invalid --memory parameter '%s' (minimum %u).%s", argv[cur_arg + 1], WKSPACE_MIN_MB, errstr_append);
	  } else {
	    sprintf(logbuf, "Error: Invalid --memory parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  }
	  goto main_ret_INVALID_CMDLINE_3;
	}
#ifndef __LP64__
	if (malloc_size_mb > 2944) {
	  logprint("Error: --memory parameter too large for 32-bit version (max 2944).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
#endif
      } else if (!memcmp(argptr2, "af", 3)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &min_maf)) {
	    sprintf(logbuf, "Error: Invalid --maf parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (min_maf <= 0.0) {
	    sprintf(logbuf, "Error: --maf parameter '%s' too small (must be > 0).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (min_maf > max_maf) {
	    sprintf(logbuf, "Error: --maf parameter '%s' too large (must be <= %g).%s", argv[cur_arg + 1], max_maf, errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  min_maf = 0.01;
	}
      } else if (!memcmp(argptr2, "ax-maf", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &max_maf)) {
	  sprintf(logbuf, "Error: Invalid --max-maf parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (max_maf < min_maf) {
	  sprintf(logbuf, "Error: --max-maf parameter '%s' too small (must be >= %g).%s", argv[cur_arg + 1], min_maf, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (max_maf >= 0.5) {
	  sprintf(logbuf, "Error: --max-maf parameter '%s' too large (must be < 0.5).%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ind", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &mind_thresh)) {
	    sprintf(logbuf, "Error: Invalid --mind parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((mind_thresh < 0.0) || (mind_thresh > 1.0)) {
	    sprintf(logbuf, "Error: Invalid --mind parameter '%s' (must be between 0 and 1).%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  mind_thresh = 0.1;
	}
      } else if (!memcmp(argptr2, "ake-grm", 8)) {
	logprint("Error: --make-grm has been retired due to inconsistent meaning across GCTA\nversions.  Use --make-grm-gz or --make-grm-bin.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ake-grm-gz", 11)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rel_calc_type |= REL_CALC_GZ | REL_CALC_GRM;
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-grm-gz 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm-gz 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + uii], "no-gz")) {
	    rel_calc_type &= ~REL_CALC_GZ;
	  } else if ((!strcmp(argv[cur_arg + uii], "ibc2")) || (!strcmp(argv[cur_arg + uii], "ibc3"))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-grm-gz 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-grm-gz '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + uii][3] - '0';
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    rel_calc_type |= REL_CALC_SINGLE_PREC;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-grm-gz parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-grm-bin", 12)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  sprintf(logbuf, "Error: --make-grm-bin cannot be used with --make-grm-gz.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rel_calc_type |= REL_CALC_GRM_BIN | REL_CALC_SINGLE_PREC;
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-grm-bin 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if ((!strcmp(argv[cur_arg + 1], "ibc2")) || (!strcmp(argv[cur_arg + 1], "ibc3"))) {
	    ibc_type = argv[cur_arg + 1][3] - '0';
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-grm-bin parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "ake-rel", 8)) {
	if (calculation_type & CALC_RELATIONSHIP) {
	  sprintf(logbuf, "Error: --make-rel cannot be used with --make-grm-gz/--make-grm-bin.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "cov")) {
	    if (calculation_type & CALC_IBC) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with --ibc flag.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_COV;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (rel_calc_type & REL_CALC_BIN) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (rel_calc_type & REL_CALC_GZ) {
	      sprintf(logbuf, "Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_TRI) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	      sprintf(logbuf, "Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ0) {
	      sprintf(logbuf, "Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    rel_calc_type |= REL_CALC_TRI;
	  } else if ((!strcmp(argv[cur_arg + uii], "ibc2")) || (!strcmp(argv[cur_arg + uii], "ibc3"))) {
	    if (rel_calc_type & REL_CALC_COV) {
	      sprintf(logbuf, "Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (ibc_type) {
	      sprintf(logbuf, "Error: --make-rel '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    ibc_type = argv[cur_arg + uii][3] - '0';
	  } else if (!strcmp(argv[cur_arg + uii], "single-prec")) {
	    rel_calc_type |= REL_CALC_SINGLE_PREC;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-rel parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (!(rel_calc_type & REL_CALC_SHAPEMASK)) {
	  rel_calc_type |= (rel_calc_type & REL_CALC_BIN)? REL_CALC_SQ : REL_CALC_TRI;
	}
	calculation_type |= CALC_RELATIONSHIP;
      } else if (!memcmp(argptr2, "atrix", 6)) {
	if (exponent != 0.0) {
	  sprintf(logbuf, "Error: --matrix cannot be used with --exponent.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (dist_calc_type & DISTANCE_IBS) {
	  sprintf(logbuf, "Error: --matrix cannot be used with '--distance ibs'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_PLINK_IBS_MATRIX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "af-succ", 8)) {
	misc_flags |= MISC_MAF_SUCC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ap3", 4)) {
	logprint("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ake-bed", 8)) {
        if (misc_flags & MISC_KEEP_AUTOCONV) {
	  sprintf(logbuf, "Error: --make-bed cannot be used with --keep-autoconv.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  sprintf(logbuf, "Error: --%s doesn't accept parameters.%s%s", argptr, ((param_ct == 1) && (!outname_end))? "  (Did you forget '--out'?)" : "", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	calculation_type |= CALC_MAKE_BED;
      } else if (!memcmp(argptr2, "erge", 5)) {
	if (calculation_type & CALC_MERGE) {
	  sprintf(logbuf, "Error: --merge cannot be used with --bmerge.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare & LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --merge does not currently support .cnv filesets.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (param_ct == 2) {
	  if (++jj > FNAMESIZE) {
	    logprint("Error: --merge .ped filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename1, argv[cur_arg + 1], jj);
	  jj = strlen(argv[cur_arg + 2]) + 1;
	  if (jj > FNAMESIZE) {
	    logprint("Error: --merge .map filename too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(mergename2, argv[cur_arg + 2], jj);
	} else {
	  if (jj > (FNAMESIZE - 5)) {
	    logprint("Error: --merge filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(mergename1, argv[cur_arg + 1], jj), ".ped", 5);
	  memcpy(memcpya(mergename2, argv[cur_arg + 1], jj), ".map", 5);
	}
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-list", 10)) {
	if (calculation_type & CALC_MERGE) {
	  logprint("Error: --merge-list cannot be used with --merge or --bmerge.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]) + 1;
	if (jj > FNAMESIZE) {
	  logprint("Error: --merge-list filename too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(mergename1, argv[cur_arg + 1], jj);
	merge_type |= MERGE_LIST;
	calculation_type |= CALC_MERGE;
      } else if (!memcmp(argptr2, "erge-mode", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((cc < '1') || (cc > '7') || (argv[cur_arg + 1][1] != '\0')) {
          sprintf(logbuf, "Error: Invalid --merge-mode parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((merge_type & MERGE_LIST) && (cc > '5')) {
	  sprintf(logbuf, "Error: --merge-mode 6-7 cannot be used with --merge-list.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        merge_type |= cc - '0';
      } else if (!memcmp(argptr2, "erge-equal-pos", 15)) {
	merge_type |= MERGE_EQUAL_POS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ust-have-sex", 13)) {
        sex_missing_pheno |= MUST_HAVE_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "covar", 6)) {
        if (!(calculation_type & CALC_GXE)) {
	  logprint("Error: --mcovar must be used with --covar and --gxe.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (gxe_mcovar > 1) {
	  logprint("Error: --mcovar cannot be used with a --gxe parameter.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mcovar parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --mcovar flag deprecated.  Use '--gxe [covariate index]'.\n");
	gxe_mcovar = ii;
      } else if (!memcmp(argptr2, "odel", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 5)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_ASSOC_FDEPR) {
	  model_modifier &= ~(MODEL_ASSOC | MODEL_ASSOC_FDEPR);
	} else if (model_modifier & MODEL_ASSOC) {
	  logprint("Error: --model cannot be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "fisher")) {
	    if (model_modifier & MODEL_TRENDONLY) {
	      sprintf(logbuf, "Error: --model 'fisher' and 'trend-only' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_FISHER;
	  } else if (!strcmp(argv[cur_arg + uii], "perm")) {
	    if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: --model 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PERM;
	  } else if (!strcmp(argv[cur_arg + uii], "genedrop")) {
	    model_modifier |= MODEL_GENEDROP;
	  } else if (!strcmp(argv[cur_arg + uii], "perm-count")) {
	    model_modifier |= MODEL_PERM_COUNT;
	  } else if (!strcmp(argv[cur_arg + uii], "dom")) {
	    if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PDOM;
	  } else if (!strcmp(argv[cur_arg + uii], "rec")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PREC;
	  } else if (!strcmp(argv[cur_arg + uii], "gen")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND | MODEL_TRENDONLY)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    if (mtest_adjust) {
	      sprintf(logbuf, "Error: --model perm-gen cannot be used with --adjust.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_modifier |= MODEL_PGEN;
	  } else if (!strcmp(argv[cur_arg + uii], "trend")) {
	    if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PTREND;
	  } else if (!strcmp(argv[cur_arg + uii], "trend-only")) {
	    if (model_modifier & (MODEL_FISHER | MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	      logprint("Error: Conflicting --model parameters.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
	  } else if ((strlen(argv[cur_arg + uii]) > 6) && (!memcmp(argv[cur_arg + uii], "mperm=", 6))) {
	    if (model_modifier & MODEL_PERM) {
	      sprintf(logbuf, "Error: --model 'mperm' and 'perm' cannot be used together.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (model_modifier & MODEL_MPERM) {
	      sprintf(logbuf, "Error: Duplicate --model 'mperm' modifier.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    kk = atoi(&(argv[cur_arg + uii][6]));
	    if (kk < 1) {
	      sprintf(logbuf, "Error: Invalid --model mperm parameter '%s'.%s", &(argv[cur_arg + uii][6]), errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    model_mperm_val = (uint32_t)kk;
	    model_modifier |= MODEL_MPERM;
	  } else if (!strcmp(argv[cur_arg + uii], "mperm")) {
	    logprint("Error: Improper --model mperm syntax.  (Use '--model mperm=[value]'.)\n");
	    goto main_ret_INVALID_CMDLINE;
	  } else {
	    sprintf(logbuf, "Error: Invalid --model parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_MODEL;
      } else if (!memcmp(argptr2, "odel-dom", 9)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-dom must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PREC | MODEL_PGEN | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-dom flag deprecated.  Use '--model dom'.\n");
	model_modifier |= MODEL_PDOM;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-gen", 9)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-gen must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-gen flag deprecated.  Use '--model gen'.\n");
	model_modifier |= MODEL_PGEN;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-rec", 9)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-rec must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PTREND)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-rec flag deprecated.  Use '--model rec'.\n");
	model_modifier |= MODEL_PREC;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "odel-trend", 11)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  logprint("Error: --model-trend must be used with --model.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & (MODEL_PDOM | MODEL_PGEN | MODEL_PREC)) {
	  logprint("Error: Conflicting --model parameters.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --model-trend flag deprecated.  Use '--model trend'.\n");
	model_modifier |= MODEL_PTREND;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "perm", 5)) {
	if (model_modifier & (MODEL_PERM | MODEL_MPERM)) {
	  sprintf(logbuf, "Error: --mperm cannot be used with --%s %sperm.%s", (model_modifier & MODEL_ASSOC)? "assoc" : "model", (model_modifier & MODEL_PERM)? "" : "m", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_PERM | GLM_MPERM)) {
	  sprintf(logbuf, "Error: --mperm cannot be used with --%s %sperm.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_PERM)? "" : "m", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mperm parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mperm_val = (uint32_t)ii;
	if (load_rare & LOAD_RARE_CNV) {
	  if ((cnv_calc_type & CNV_INDIV_PERM) && (!cnv_indiv_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-indiv-perm [perm. count]'.\n");
	    cnv_indiv_mperms = mperm_val;
	  } else if ((cnv_calc_type & CNV_TEST_REGION) && (!cnv_test_region_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-test-region [perm. count]'.\n");
	  } else if ((cnv_calc_type & CNV_ENRICHMENT_TEST) && (!cnv_enrichment_test_mperms)) {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-enrichment-test [perm. count]'.\n");
	  } else {
	    logprint("Note: --mperm flag deprecated.  Use e.g. '--cnv-test [permutation count]'.\n");
            if (!(cnv_calc_type & (CNV_INDIV_PERM | CNV_ENRICHMENT_TEST | CNV_TEST | CNV_TEST_REGION))) {
	      cnv_calc_type |= CNV_TEST;
	    }
	    cnv_test_mperms = mperm_val;
	  }
	  // if e.g. --cnv-test-region had a valid parameter, don't clobber it
	  if (!cnv_test_region_mperms) {
	    cnv_test_region_mperms = mperm_val;
	  }
	  if (!cnv_enrichment_test_mperms) {
	    cnv_enrichment_test_mperms = mperm_val;
	  }
	} else {
	  logprint("Note: --mperm flag deprecated.  Use e.g. '--model mperm=[value]'.\n");
	  model_mperm_val = mperm_val;
	  model_modifier |= MODEL_MPERM;
	  glm_mperm_val = mperm_val;
	  glm_modifier |= GLM_MPERM;
	}
      } else if (!memcmp(argptr2, "perm-save", 10)) {
	if (glm_modifier & GLM_NO_SNP) {
          sprintf(logbuf, "Error: --mperm-save cannot be used with --linear/--logistic no-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	mperm_save |= MPERM_DUMP_BEST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "perm-save-all", 14)) {
	mperm_save |= MPERM_DUMP_ALL;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "c", 2)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mc must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --mc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.max_size = ii;
      } else if (!memcmp(argptr2, "cc", 2)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mcc must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mcc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (((uint32_t)ii) > cluster.max_size) {
          logprint("Error: --mcc parameter exceeds --mc parameter.\n");
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.max_cases = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mcc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (((uint32_t)ii) > cluster.max_size) {
          logprint("Error: --mcc parameter exceeds --mc parameter.\n");
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.max_ctrls = ii;
      } else if (!memcmp(argptr2, "atch", 5)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --match must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cluster.match_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
	  if (alloc_string(&cluster.match_missing_str, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "atch-type", 10)) {
	if (!cluster.match_fname) {
	  sprintf(logbuf, "Error: --match-type must be used with --match.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&cluster.match_type_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ds-plot", 8)) {
#ifdef NOLAPACK
	// PLINK 1.07's SVD-based non-LAPACK implementation does not conform to
	// classical MDS, so we do not replicate it.
        logprint("Error: --mds-plot requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mds-plot must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cluster.mds_dim_ct = 0;
        for (uii = 1; uii <= param_ct; uii++) {
          if (!strcmp(argv[cur_arg + uii], "by-cluster")) {
	    cluster.modifier |= CLUSTER_MDS;
	  } else if (!strcmp(argv[cur_arg + uii], "eigvals")) {
	    cluster.modifier |= CLUSTER_MDS_EIGVALS;
	  } else {
	    ii = atoi(argv[cur_arg + uii]);
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --mds-plot parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
              goto main_ret_INVALID_CMDLINE_3;
	    } else if (cluster.mds_dim_ct) {
	      sprintf(logbuf, "Error: Invalid --mds-plot parameter sequence.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    cluster.mds_dim_ct = ii;
	  }
	}
#endif
      } else if (!memcmp(argptr2, "ds-cluster", 11)) {
	if (!(calculation_type & CALC_CLUSTER)) {
	  sprintf(logbuf, "Error: --mds-cluster must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --mds-cluster flag deprecated.  Use '--mds-plot by-cluster'.\n");
        cluster.modifier |= CLUSTER_MDS;
      } else if (!memcmp(argptr2, "within", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --mwithin parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        mwithin_col = ii;
      } else if (!memcmp(argptr2, "in", 3)) {
        if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --min must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --min parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < -1.0) || (dxx > 1.0)) {
          sprintf(logbuf, "Error: --min threshold must be between -1 and 1 inclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (dxx > genome_max_pi_hat) {
	  logprint("Error: --min value cannot be greater than --max value.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        genome_modifier |= GENOME_FILTER_PI_HAT;
	genome_min_pi_hat = dxx;
      } else if (!memcmp(argptr2, "ax", 3)) {
        if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --max must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --max parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx < -1.0) || (dxx > 1.0)) {
          sprintf(logbuf, "Error: --max threshold must be between -1 and 1 inclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	genome_modifier |= GENOME_FILTER_PI_HAT;
	genome_max_pi_hat = dxx;
      } else if (!memcmp(argptr2, "ake-founders", 13)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "require-2-missing")) {
	    misc_flags |= MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING;
	  } else if (!strcmp(argv[cur_arg + uii], "first")) {
	    misc_flags |= MISC_MAKE_FOUNDERS_FIRST;
	  } else {
	    sprintf(logbuf, "Error: Invalid --make-founders parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
        misc_flags |= MISC_MAKE_FOUNDERS;
      } else if (!memcmp(argptr2, "issing", 7)) {
	calculation_type |= CALC_MISSING_REPORT;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "h", 2)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "bd")) {
	    sprintf(logbuf, "Error: Invalid --mh parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  misc_flags |= MISC_CMH_BD;
	}
	calculation_type |= CALC_CMH;
	logprint("Error: --mh is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "h2", 3)) {
	if (calculation_type & CALC_CMH) {
	  logprint("Error: --mh2 cannot be used with --mh.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	calculation_type |= CALC_CMH;
        misc_flags |= MISC_CMH2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "lma", 4)) {
        logprint("Error: --mlma is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "lma-loco", 9)) {
        logprint("Error: --mlma-loco is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "lma-no-adj-covar", 17)) {
        logprint("Error: --mlma-no-adj-covar is not implemented yet.\n");
        goto main_ret_INVALID_CMDLINE;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'n':
      if (!memcmp(argptr2, "o-fid", 6)) {
	fam_cols &= ~FAM_COL_1;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-parents", 10)) {
	fam_cols &= ~FAM_COL_34;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-sex", 6)) {
	if (filter_binary & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
	  logprint("Error: --filter-males/--filter-females cannot be used with --no-sex.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	fam_cols &= ~FAM_COL_5;
	sex_missing_pheno |= ALLOW_NO_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-pheno", 8)) {
	fam_cols &= ~FAM_COL_6;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "onfounders", 11)) {
	misc_flags |= MISC_NONFOUNDERS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "eighbour", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --neighbour parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	neighbor_n1 = ii;
	ii = atoi(argv[cur_arg + 2]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --neighbour parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	neighbor_n2 = ii;
	if (neighbor_n2 < neighbor_n1) {
	  sprintf(logbuf, "Error: Second --neighbour parameter cannot be smaller than first parameter.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_NEIGHBOR;
      } else if (!memcmp(argptr2, "ot-chr", 7)) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --from cannot be used with --autosome[-xy] or --[not-]chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	// allowed:
	//   --allow-extra-chr --chr 5-22 bobs_chrom --not-chr 17
	// allowed:
	//   --allow-extra-chr --not-chr 12-17 bobs_chrom
	// does not make sense, disallowed:
	//   --allow-extra-chr --chr 5-22 --not-chr bobs_chrom

	// --allow-extra-chr present, --chr/--autosome[-xy] not present
	uii = ((misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1) && (!chrom_info.is_include_stack);
	retval = parse_chrom_ranges(param_ct, '-', &(argv[cur_arg]), chrom_exclude, &chrom_info, uii, argptr);
	if (retval) {
	  goto main_ret_1;
	}
	if (chrom_info.is_include_stack) {
	  fill_chrom_mask(&chrom_info);
	}
	for (uii = 0; uii < CHROM_MASK_INITIAL_WORDS; uii++) {
	  chrom_info.chrom_mask[uii] &= ~chrom_exclude[uii];
	}
	if (all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS) && ((!((misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1)) || (chrom_info.is_include_stack && (!chrom_info.incl_excl_name_stack)))) {
	  sprintf(logbuf, "Error: All chromosomes excluded.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	chrom_flag_present = 1;
      } else if (!memcmp(argptr2, "udge", 5)) {
        if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --nudge must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --nudge flag deprecated.  Use '--genome nudge'.\n");
        genome_modifier |= GENOME_NUDGE;
        goto main_param_zero;
      } else if (!memcmp(argptr2, "o-snp", 6)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --no-snp must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (mtest_adjust) {
	  sprintf(logbuf, "Error: --no-snp cannot be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (mperm_save & MPERM_DUMP_BEST) {
	  sprintf(logbuf, "Error: --no-snp cannot be used with --mperm-save.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((glm_modifier & (GLM_NO_SNP_EXCL - GLM_HETHOM - GLM_DOMINANT)) || ((glm_modifier & (GLM_HETHOM | GLM_DOMINANT)) && (!(glm_modifier & (GLM_CONDITION_DOMINANT | GLM_CONDITION_RECESSIVE))))) {
	  sprintf(logbuf, "Error: --no-snp conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --no-snp flag deprecated.  Use e.g. '--linear no-snp'.\n");
        glm_modifier |= GLM_NO_SNP;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "o-x-sex", 8)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --no-x-sex must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_NO_SNP | GLM_SEX)) {
	  sprintf(logbuf, "Error: --no-x-sex conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --no-x-sex flag deprecated.  Use e.g. '--linear no-x-sex'.\n");
	glm_modifier |= GLM_NO_X_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "oweb", 5)) {
        logprint("Note: --noweb has no effect since no web check is implemented yet.\n");
	goto main_param_zero;
      } else {
        goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'o':
      if (!memcmp(argptr2, "utput-missing-genotype", 23)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argv[cur_arg + 1][0];
	if ((argv[cur_arg + 1][1] != '\0') || (((unsigned char)cc) <= ' ')) {
	  sprintf(logbuf, "Error: Invalid --output-missing-genotype parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	g_output_missing_geno_ptr = &(g_one_char_strs[((unsigned char)cc) * 2]);
      } else if (!memcmp(argptr2, "utput-missing-phenotype", 24)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > 31) {
	  logprint("Error: --output-missing-phenotype string too long (max 31 chars).\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  logprint("Error: --output-missing-phenotype parameter currently must be numeric.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	memcpy(output_missing_pheno, argv[cur_arg + 1], jj + 1);
      } else if (memcmp(argptr2, "ut", 3)) {
	// --out is a special case due to logging
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'p':
      if (!memcmp(argptr2, "ed", 3)) {
	if ((load_params & 0x3fa) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 2;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --ped parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(pedname, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "heno", 5)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (makepheno_str) {
	  logprint("Error: --pheno and --make-pheno flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&phenoname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "heno-name", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (mpheno_col != 0) {
	  logprint("Error: --mpheno and --pheno-name flags cannot coexist.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (alloc_string(&phenoname_str, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "heno-merge", 11)) {
	pheno_modifier |= PHENO_MERGE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rune", 5)) {
	misc_flags |= MISC_PRUNE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "arallel", 8)) {
	if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	  sprintf(logbuf, "Error: --parallel cannot be used with '--distance square'.  Use '--distance\nsquare0' or plain --distance instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((dist_calc_type & DISTANCE_BIN) && (!(dist_calc_type & DISTANCE_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --parallel cannot be used with plain '--distance bin'.  Use '--distance\nbin square0' or '--distance bin triangle' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((rel_calc_type & REL_CALC_SHAPEMASK) == REL_CALC_SQ) {
	  sprintf(logbuf, "Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain '--make-rel' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((rel_calc_type & REL_CALC_BIN) && (!(rel_calc_type & REL_CALC_SHAPEMASK))) {
	  sprintf(logbuf, "Error: --parallel cannot be used with plain '--make-rel bin'.  Use '--make-rel\nbin square0' or '--make-rel bin triangle' instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) {
	  sprintf(logbuf, "Error: --parallel and --distance-matrix/--matrix cannot be used together.  Use\n--distance instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_GROUPDIST) {
	  sprintf(logbuf, "Error: --parallel and --groupdist cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_CLUSTER) {
	  sprintf(logbuf, "Error: --parallel and --cluster cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_NEIGHBOR) {
	  sprintf(logbuf, "Error: --parallel and --neighbour cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 2, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if ((ii < 1) || (ii > PARALLEL_MAX)) {
	  sprintf(logbuf, "Error: Invalid --parallel job index '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	parallel_idx = ii - 1; // internal 0..(n-1) indexing
	ii = atoi(argv[cur_arg + 2]);
	if ((ii < 2) || (ii > PARALLEL_MAX) || (((uint32_t)ii) < parallel_idx)) {
	  sprintf(logbuf, "Error: Invalid --parallel total job count '%s'.%s", argv[cur_arg + 2], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        parallel_tot = ii;
      } else if (!memcmp(argptr2, "pc-gap", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ppc-gap parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	dxx *= 1000;
	if (dxx < 0) {
	  ppc_gap = 0;
	} else if (dxx > 2147483647) {
	  ppc_gap = 0x7fffffff;
	} else {
	  ppc_gap = (int32_t)(dxx + EPSILON);
	}
      } else if (!memcmp(argptr2, "erm", 4)) {
	if ((model_modifier & MODEL_MPERM) && (calculation_type & CALC_MODEL)) {
	  sprintf(logbuf, "Error: --perm cannot be used with --%s mperm.%s", (model_modifier & MODEL_ASSOC)? "assoc" : "model", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((calculation_type & CALC_GLM) && (glm_modifier & (GLM_MPERM + GLM_NO_SNP))) {
	  sprintf(logbuf, "Error: --perm cannot be used with --%s %s.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_MPERM)? "mperm" : "no-snp", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (model_modifier & MODEL_MPERM) {
	  sprintf(logbuf, "Error: --perm cannot be used with --mperm.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	model_modifier |= MODEL_PERM;
        glm_modifier |= GLM_PERM;
	logprint("Note: --perm flag deprecated.  Use e.g. '--model perm'.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "erm-count", 10)) {
	model_modifier |= MODEL_PERM_COUNT;
	glm_modifier |= GLM_PERM_COUNT;
	logprint("Note: --perm-count flag deprecated.  Use e.g. '--model perm-count'.\n");
	goto main_param_zero;
      } else if (!memcmp(argptr2, "2", 2)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  logprint("Error: --p2 must be used with --assoc.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (model_modifier & MODEL_QMASK) {
	  sprintf(logbuf, "Error: --assoc 'qt-means'/'lin' does not make sense with --p2.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --p2 flag deprecated.  Use '--assoc p2 ...'.\n");
	model_modifier |= MODEL_ASSOC_P2;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "filter", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --pfilter parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --pfilter threshold must be between 0 and 1 exclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	pfilter = dxx;
      } else if (!memcmp(argptr2, "erm-batch-size", 1)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	perm_batch_size = atoi(argv[cur_arg + 1]);
	if ((perm_batch_size < 1) || (perm_batch_size > 0x7fffffff)) {
	  sprintf(logbuf, "Error: Invalid --perm-batch-size parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "pc", 3)) {
	if (!(calculation_type & (CALC_NEIGHBOR | CALC_CLUSTER))) {
          sprintf(logbuf, "Error: --ppc must be used with --cluster or --neigbour.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --ppc parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((dxx <= 0.0) || (dxx >= 1.0)) {
	  sprintf(logbuf, "Error: --ppc threshold must be between 0 and 1 exclusive.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        cluster.ppc = dxx;
      } else if (!memcmp(argptr2, "ool-size", 9)) {
	if (!(homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE))) {
          logprint("Error: --pool-size must be used with --homozyg group[-verbose].\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 2) {
	  sprintf(logbuf, "Error: Invalid --pool-size parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	homozyg.pool_size_min = ii;
      } else if (!memcmp(argptr2, "arameters", 10)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --parameters must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &parameters_range_list, 1);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'q':
      if (!memcmp(argptr2, "t-means", 8)) {
	if ((!(calculation_type & CALC_MODEL)) || (!(model_modifier & MODEL_ASSOC))) {
	  sprintf(logbuf, "Error: --qt-means must be used with --assoc.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_DMASK) {
	  sprintf(logbuf, "Error: --qt-means does not make sense with a case/control-specific --assoc\nmodifier.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --qt-means flag deprecated.  Use '--assoc qt-means ...'.\n");
	model_modifier |= MODEL_QT_MEANS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "q-plot", 7)) {
        if (!mtest_adjust) {
	  sprintf(logbuf, "Error: --qq-plot must be used with --adjust.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --qq-plot flag deprecated.  Use '--adjust qq-plot'.\n");
	mtest_adjust |= ADJUST_QQ;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "match", 6)) {
        if (!(calculation_type & CALC_CLUSTER)) {
          sprintf(logbuf, "Error: --qmatch must be used with --cluster.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&cluster.qmatch_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
	  if (alloc_string(&cluster.qmatch_missing_str, argv[cur_arg + 2])) {
	    goto main_ret_NOMEM;
	  }
	}
      } else if (!memcmp(argptr2, "t", 2)) {
        if (!cluster.qmatch_fname) {
          sprintf(logbuf, "Error: --qt must be used with --qmatch.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&cluster.qt_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'r':
      if (!memcmp(argptr2, "emove", 6)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&removename, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "emove-fam", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&removefamname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "emove-clusters", 15)) {
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&(cluster.remove_fname), argv[cur_arg + 1], argptr, 0);
        if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "emove-cluster-names", 20)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_and_flatten(&(cluster.remove_flattened), &(argv[cur_arg + 1]), param_ct);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "el-cutoff", 10)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel cannot be used with %s.  (Use a combination of\n--make-rel, --keep/--remove, and a filtering script.)%s", argptr, errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (scan_double(argv[cur_arg + 1], &rel_cutoff) || (rel_cutoff <= 0.0) || (rel_cutoff >= 1.0)) {
	    sprintf(logbuf, "Error: Invalid %s parameter '%s'.%s", argptr, argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_REL_CUTOFF;
      } else if (!memcmp(argptr2, "egress-distance", 16)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --regress-distance cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  regress_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((regress_iters < 2) || (regress_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --regress-distance jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    ii = atoi(argv[cur_arg + 2]);
	    if (ii <= 0) {
	      sprintf(logbuf, "Error: Invalid --regress-distance jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    regress_d = ii;
	  }
	}
	calculation_type |= CALC_REGRESS_DISTANCE;
      } else if (!memcmp(argptr2, "egress-rel", 11)) {
	if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --regress-rel flags cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (rel_calc_type & REL_CALC_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --regress-rel cannot currently be used with a single-precision\nrelationship matrix.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  regress_rel_iters = strtoul(argv[cur_arg + 1], NULL, 10);
	  if ((regress_rel_iters < 2) || (regress_rel_iters == ULONG_MAX)) {
	    sprintf(logbuf, "Error: Invalid --regress-rel jackknife iteration count '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if (param_ct == 2) {
	    ii = atoi(argv[cur_arg + 2]);
	    if (ii <= 0) {
	      sprintf(logbuf, "Error: Invalid --regress-rel jackknife delete parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    regress_rel_d = ii;
	  }
	}
	calculation_type |= CALC_REGRESS_REL;
      } else if (!memcmp(argptr2, "egress-pcs", 11)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 5)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&evecname, argv[cur_arg + 1], argptr, 9);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "normalize-pheno")) {
	    regress_pcs_modifier |= REGRESS_PCS_NORMALIZE_PHENO;
	  } else if (!strcmp(argv[cur_arg + uii], "sex-specific")) {
	    regress_pcs_modifier |= REGRESS_PCS_SEX_SPECIFIC;
	  } else if (!strcmp(argv[cur_arg + uii], "clip")) {
	    regress_pcs_modifier |= REGRESS_PCS_CLIP;
	  } else if ((max_pcs != MAX_PCS_DEFAULT) || (argv[cur_arg + uii][0] < '0') || (argv[cur_arg + uii][0] > '9')) {
	    sprintf(logbuf, "Error: Invalid --regress-pcs parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    ii = atoi(argv[cur_arg + uii]);
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --regress-pcs maximum principal component count '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    max_pcs = ii;
	  }
	}
	calculation_type |= CALC_REGRESS_PCS;
      } else if (!memcmp(argptr2, "egress-pcs-distance", 20)) {
	if (calculation_type & CALC_REGRESS_PCS) {
	  sprintf(logbuf, "Error: --regress-pcs-distance cannot be used with --regress-pcs.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (calculation_type & CALC_DISTANCE) {
	  sprintf(logbuf, "Error: --regress-pcs-distance cannot be used with --distance.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 11)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&evecname, argv[cur_arg + 1], argptr, 9);
	if (retval) {
	  goto main_ret_1;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "normalize-pheno")) {
	    regress_pcs_modifier |= REGRESS_PCS_NORMALIZE_PHENO;
	  } else if (!strcmp(argv[cur_arg + uii], "sex-specific")) {
	    regress_pcs_modifier |= REGRESS_PCS_SEX_SPECIFIC;
	  } else if (!strcmp(argv[cur_arg + uii], "square")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if (parallel_tot > 1) {
	      sprintf(logbuf, "Error: --parallel cannot be used with '--regress-pcs-distance square'.  Use\nthe square0 or triangle shape instead.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ;
	  } else if (!strcmp(argv[cur_arg + uii], "square0")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_TRI) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square0' and 'triangle' modifiers can't coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_SQ0;
	  } else if (!strcmp(argv[cur_arg + uii], "triangle")) {
	    if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    } else if ((dist_calc_type & DISTANCE_SHAPEMASK) == DISTANCE_SQ0) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'square0' and 'triangle' modifiers can't coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_TRI;
	  } else if (!strcmp(argv[cur_arg + uii], "gz")) {
	    if (dist_calc_type & DISTANCE_BIN) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_GZ;
	  } else if (!strcmp(argv[cur_arg + uii], "bin")) {
	    if (dist_calc_type & DISTANCE_GZ) {
	      sprintf(logbuf, "Error: --regress-pcs-distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    dist_calc_type |= DISTANCE_BIN;
	  } else if (!strcmp(argv[cur_arg + uii], "ibs")) {
	    if (dist_calc_type & DISTANCE_IBS) {
	      logprint("Error: Duplicate --regress-pcs-distance 'ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "1-ibs")) {
	    if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
	      logprint("Error: Duplicate --regress-pcs-distance '1-ibs' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_1_MINUS_IBS;
	  } else if (!strcmp(argv[cur_arg + uii], "allele-ct")) {
	    if (dist_calc_type & DISTANCE_ALCT) {
	      logprint("Error: Duplicate --regress-pcs-distance 'allele-ct' modifier.\n");
	      goto main_ret_INVALID_CMDLINE;
	    }
	    dist_calc_type |= DISTANCE_ALCT;
	  } else if (!strcmp(argv[cur_arg + uii], "3d")) {
	    dist_calc_type |= DISTANCE_3D;
	  } else if (!strcmp(argv[cur_arg + uii], "flat-missing")) {
	    dist_calc_type |= DISTANCE_FLAT_MISSING;
	  } else if ((max_pcs != MAX_PCS_DEFAULT) || (argv[cur_arg + uii][0] < '0') || (argv[cur_arg + uii][0] > '9')) {
	    sprintf(logbuf, "Error: Invalid --regress-pcs-distance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else {
	    ii = atoi(argv[cur_arg + uii]);
	    if (ii < 1) {
	      sprintf(logbuf, "Error: Invalid --regress-pcs-distance maximum PC count '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    max_pcs = ii;
	  }
	}
	if (!(dist_calc_type & DISTANCE_TYPEMASK)) {
	  dist_calc_type |= DISTANCE_ALCT;
	}
	calculation_type |= CALC_REGRESS_PCS_DISTANCE;
      } else if (!memcmp(argptr2, "ead-freq", 9)) {
	if (calculation_type & CALC_FREQ) {
	  sprintf(logbuf, "Error: --freq and --read-freq flags cannot coexist.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&freqname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "ecode", 6)) || (!memcmp(argptr2, "ecode ", 6))) {
	if (argptr2[5] == ' ') {
	  kk = 1;
	} else {
	  kk = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4 - kk)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (recode_modifier & (RECODE_A | RECODE_AD)) {
	      sprintf(logbuf, "Error: --recode '12' modifier cannot be used with 'A' or 'AD'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_12;
	  } else if (!strcmp(argv[cur_arg + uii], "compound-genotypes")) {
	    if (recode_modifier & RECODE_STRUCTURE) {
              logprint("Error: --recode 'compound-genotypes' modifier cannot be used with 'structure'.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    recode_modifier |= RECODE_COMPOUND;
	  } else if (!strcmp(argv[cur_arg + uii], "23")) {
	    if (recode_type_set(&recode_modifier, RECODE_23)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
          } else if ((!argv[cur_arg + uii][1]) && (tolower(argv[cur_arg + uii][0]) == 'a')) {
	    if (recode_type_set(&recode_modifier, RECODE_A)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if ((!argv[cur_arg + uii][2]) && match_upper(argv[cur_arg + uii], "AD")) {
	    if (recode_type_set(&recode_modifier, RECODE_AD)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (match_upper(argv[cur_arg + uii], "HV")) {
	    if (!argv[cur_arg + uii][2]) {
	      if (recode_type_set(&recode_modifier, RECODE_HV)) {
	        goto main_ret_INVALID_CMDLINE_3;
	      }
	    } else if (!strcmp(&(argv[cur_arg + uii][2]), "-1chr")) {
	      if (recode_type_set(&recode_modifier, RECODE_HV_1CHR)) {
	        goto main_ret_INVALID_CMDLINE_3;
	      }
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "tab")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TAB;
	  } else if (!strcmp(argv[cur_arg + uii], "tabx")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_TAB | RECODE_DELIMX;
	  } else if (!strcmp(argv[cur_arg + uii], "spacex")) {
	    if (recode_modifier & (RECODE_TAB | RECODE_DELIMX)) {
	      sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_DELIMX;
	  } else if (!strcmp(argv[cur_arg + uii], "beagle")) {
	    if (recode_type_set(&recode_modifier, RECODE_BEAGLE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "bimbam")) {
	    if (recode_type_set(&recode_modifier, RECODE_BIMBAM)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "bimbam-1chr")) {
	    if (recode_type_set(&recode_modifier, RECODE_BIMBAM_1CHR)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "fastphase")) {
	    if (recode_type_set(&recode_modifier, RECODE_FASTPHASE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "fastphase-1chr")) {
	    if (recode_type_set(&recode_modifier, RECODE_FASTPHASE_1CHR)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "lgen")) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "lgen-ref")) {
	    if (recode_type_set(&recode_modifier, RECODE_LGEN_REF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "list")) {
	    if (recode_type_set(&recode_modifier, RECODE_LIST)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "rlist")) {
	    if (recode_type_set(&recode_modifier, RECODE_RLIST)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "structure")) {
	    if (recode_modifier & RECODE_COMPOUND) {
              logprint("Error: --recode 'compound-genotypes' modifier cannot be used with 'structure'.\n");
              goto main_ret_INVALID_CMDLINE;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_STRUCTURE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "transpose")) {
	    if (recode_type_set(&recode_modifier, RECODE_TRANSPOSE)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "vcf")) {
	    if (recode_modifier & RECODE_VCF) {
	      sprintf(logbuf, "Error: Conflicting or redundant --recode modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	  } else if (!strcmp(argv[cur_arg + uii], "vcf-fid")) {
	    if (recode_modifier & (RECODE_VCF | RECODE_IID)) {
	      sprintf(logbuf, "Error: Conflicting or redundant --recode modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_FID;
	  } else if (!strcmp(argv[cur_arg + uii], "vcf-iid")) {
	    if (recode_modifier & (RECODE_VCF | RECODE_FID)) {
	      sprintf(logbuf, "Error: Conflicting or redundant --recode modifiers.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (recode_type_set(&recode_modifier, RECODE_VCF)) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    recode_modifier |= RECODE_IID;
	  } else {
	    sprintf(logbuf, "Error: Invalid --recode parameter '%s'.%s%s", argv[cur_arg + uii], ((uii == param_ct) && (!outname_end))? "  (Did you forget '--out'?)" : "", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	calculation_type |= CALC_RECODE;
      } else if (!memcmp(argptr2, "ecode-whap", 11)) {
        logprint("Error: --recode-whap flag retired since WHAP is no longer supported.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ecode-allele", 13)) {
	if (!(recode_modifier & (RECODE_A | RECODE_AD))) {
	  sprintf(logbuf, "Error: --recode-allele must be used with --recode A or --recode AD.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (recode_modifier & RECODE_12) {
	  sprintf(logbuf, "Error: --recode-allele cannot be used with --recode 12.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&recode_allele_name, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "eference", 9)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&lgen_reference_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	lgen_modifier |= LGEN_REFERENCE;
      } else if (!memcmp(argptr2, "ead-genome", 11)) {
	if (calculation_type & CALC_GENOME) {
          sprintf(logbuf, "Error: --read-genome cannot be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (!(calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR))) {
          sprintf(logbuf, "Error: --read-genome cannot be used without --cluster or --neighbour.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	} else if ((cluster.ppc == 0.0) && ((calculation_type & (CALC_DISTANCE | CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) || read_dists_fname)) {
          sprintf(logbuf, "Error: --read-genome is pointless with --distance/--[distance-]matrix or\n--read-dists unless --ppc is also present.%s", errstr_append);
          goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        retval = alloc_fname(&read_genome_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ead-genome-list", 19)) {
	logprint("Error: --read-genome-list flag retired.  Use --parallel + Unix cat instead.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ead-genome-minimal", 19)) {
	logprint("Error: --read-genome-minimal flag retired.  Use --genome gz + --read-genome\ninstead.");
	goto main_ret_INVALID_CMDLINE;
      } else if (!memcmp(argptr2, "ead-dists", 10)) {
	if (calculation_type & (CALC_DISTANCE | CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) {
	  sprintf(logbuf, "Error: --read-dists cannot be used with a distance matrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (cluster.modifier & CLUSTER_MISSING) {
          sprintf(logbuf, "Error: --read-dists cannot be used with '--cluster missing'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&read_dists_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
        if (param_ct == 2) {
          retval = alloc_fname(&read_dists_id_fname, argv[cur_arg + 2], argptr, 0);
          if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "el-check", 9)) {
        if (!(calculation_type & CALC_GENOME)) {
          sprintf(logbuf, "Error: --rel-check must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        logprint("Note: --rel-check flag deprecated.  Use '--genome rel-check'.\n");
        genome_modifier |= GENOME_REL_CHECK;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ecessive", 9)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --recessive must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT)) {
	  sprintf(logbuf, "Error: --recessive conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --recessive flag deprecated.  Use e.g. '--linear recessive' (and\n'--condition-list [filename] recessive' to change covariate coding).\n");
	glm_modifier |= GLM_RECESSIVE | GLM_CONDITION_RECESSIVE;
	glm_xchr_model = 0;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 's':
      if (!memcmp(argptr2, "ample", 6)) {
	if ((load_params & 0x27f) || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	load_params |= 0x200;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
	  logprint("Error: --sample parameter too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	strcpy(samplename, argv[cur_arg + 1]);
      } else if (!memcmp(argptr2, "eed", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 0x7fffffff)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	rseed_ct = param_ct;
	rseeds = (uint32_t*)malloc(param_ct * sizeof(int32_t));
	for (uii = 1; uii <= param_ct; uii++) {
	  if (strtoui32(argv[cur_arg + uii], &(rseeds[uii - 1]))) {
	    sprintf(logbuf, "Error: Invalid --seed parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
      } else if (!memcmp(argptr2, "np", 3)) {
        if (markername_from) {
	  sprintf(logbuf, "Error: --snp cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (marker_pos_start != -1) {
	  sprintf(logbuf, "Error: --snp cannot be used with --from-bp/-kb/-mb.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if ((!all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS)) || chrom_info.incl_excl_name_stack) {
	  sprintf(logbuf, "Error: --snp cannot be used with --autosome[-xy] or --[not-]chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
          sprintf(logbuf, "Error: --snp cannot be used with --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
          sprintf(logbuf, "Error: --snp cannot be used with --exclude-snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --snp cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb and\n--to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_snp, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "nps", 4)) {
	if (markername_from) {
	  sprintf(logbuf, "Error: --snps cannot be used with --from.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (marker_pos_start != -1) {
	  sprintf(logbuf, "Error: --snps cannot be used with --from-bp/-kb/-mb.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
	  sprintf(logbuf, "Error: --snps cannot be used with --snp or --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
	  sprintf(logbuf, "Error: --snps cannot be used with --exclude-snps.%s", errstr_append);
	} else if (load_rare == LOAD_RARE_CNV) {
	  sprintf(logbuf, "Error: --snps cannot be used with a .cnv fileset.  Use --from-bp/-kb/-mb and\n--to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	// mise well allow --snps + --autosome/--autosome-xy/--chr/--not-chr
	retval = parse_name_ranges(param_ct, range_delim, &(argv[cur_arg]), &snps_range_list, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "et-hh-missing", 14)) {
	misc_flags |= MISC_SET_HH_MISSING;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "et", 3)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&set_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ubset", 6)) {
	UNSTABLE;
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!set_fname) {
	  logprint("Error: --subset must be used with --set.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&subset_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if ((!memcmp(argptr2, "imulate", 8)) || (!memcmp(argptr2, "imulate-qt", 11))) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 3)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&simulate_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
	if (argptr2[7] == '-') {
	  simulate_flags |= SIMULATE_QT;
	}
	for (uii = 2; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "tags")) {
	    if (simulate_flags & SIMULATE_HAPS) {
	      sprintf(logbuf, "Error: --%s 'tags' and 'haps' modifiers cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    simulate_flags |= SIMULATE_TAGS;
	  } else if (!strcmp(argv[cur_arg + uii], "haps")) {
	    if (simulate_flags & SIMULATE_TAGS) {
	      sprintf(logbuf, "Error: --%s 'tags' and 'haps' modifiers cannot be used together.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    simulate_flags |= SIMULATE_HAPS;
	  } else if (match_upper(argv[cur_arg + uii], "ACGT")) {
	    if (simulate_flags & (SIMULATE_1234 | SIMULATE_12)) {
	      sprintf(logbuf, "Error: --%s 'acgt' modifier cannot be used with '1234' or '12'.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            simulate_flags |= SIMULATE_ACGT;
	  } else if (!strcmp(argv[cur_arg + uii], "1234")) {
	    if (simulate_flags & (SIMULATE_ACGT | SIMULATE_12)) {
	      sprintf(logbuf, "Error: --%s '1234' modifier cannot be used with 'acgt' or '12'.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            simulate_flags |= SIMULATE_1234;
	  } else if (!strcmp(argv[cur_arg + uii], "12")) {
	    if (simulate_flags & (SIMULATE_ACGT | SIMULATE_1234)) {
	      sprintf(logbuf, "Error: --%s '12' modifier cannot be used with 'acgt' or '1234'.%s", argptr, errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
            simulate_flags |= SIMULATE_12;
	  } else {
	    sprintf(logbuf, "Error: Invalid --%s parameter '%s'.%s", argptr, argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	load_rare = LOAD_RARE_SIMULATE;
      } else if (!memcmp(argptr2, "imulate-ncases", 15)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare != LOAD_RARE_SIMULATE) {
	  sprintf(logbuf, "Error: --simulate-ncases must be used with --simulate.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --simulate-ncases parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	simulate_cases = ii;
      } else if (!memcmp(argptr2, "imulate-ncontrols", 18)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare != LOAD_RARE_SIMULATE) {
	  sprintf(logbuf, "Error: --simulate-ncontrols must be used with --simulate.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (atoiz(argv[cur_arg + 1], &ii)) {
	  sprintf(logbuf, "Error: Invalid --simulate-ncontrols parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((!ii) && (!simulate_cases)) {
	  logprint("Error: '--simulate-ncases 0' cannot be used with '--simulate-ncontrols 0'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	simulate_controls = ii;
      } else if (!memcmp(argptr2, "imulate-prevalence", 19)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &simulate_prevalence) || (simulate_prevalence < 0) || (simulate_prevalence > 1)) {
	  sprintf(logbuf, "Error: Invalid --simulate-prevalence parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "imulate-label", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&simulate_label, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if (!memcmp(argptr2, "imulate-missing", 16)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &simulate_missing) || (simulate_missing < 0) || (simulate_missing > 1)) {
	  sprintf(logbuf, "Error: Invalid --simulate-missing parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "imulate-n", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare == LOAD_RARE_SIMULATE) {
	  sprintf(logbuf, "Error: --simulate-n must be used with --simulate-qt, not --simulate.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --simulate-n parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	simulate_qt_indivs = ii;
      } else if (!memcmp(argptr2, "imulate-haps", 13)) {
	if (simulate_flags & SIMULATE_TAGS) {
	  sprintf(logbuf, "Error: --simulate-tags cannot be used with --simulate-haps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --simulate-haps flag deprecated.  Use e.g. '--simulate haps'.\n");
	simulate_flags |= SIMULATE_HAPS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "imulate-tags", 13)) {
	if (simulate_flags & SIMULATE_HAPS) {
	  sprintf(logbuf, "Error: --simulate-tags cannot be used with --simulate-haps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --simulate-tags flag deprecated.  Use e.g. '--simulate tags'.\n");
	simulate_flags |= SIMULATE_TAGS;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ex", 3)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --sex must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & GLM_NO_X_SEX) {
	  sprintf(logbuf, "Error: --sex conflicts with a --%s modifier.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --sex flag deprecated.  Use e.g. '--linear sex'.\n");
	glm_modifier |= GLM_SEX;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "tandard-beta", 13)) {
	if ((!(calculation_type & CALC_GLM)) || (glm_modifier & GLM_LOGISTIC)) {
	  sprintf(logbuf, "Error: --standard-beta must be used wtih --linear.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --standard-beta flag deprecated.  Use '--linear standard-beta'.\n");
	glm_modifier |= GLM_STANDARD_BETA;
	goto main_param_zero;
      } else if (memcmp(argptr2, "ilent", 6)) {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 't':
      if (!memcmp(argptr2, "ail-pheno", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (makepheno_str) {
	  sprintf(logbuf, "Error: --tail-pheno cannot be used with --make-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &tail_bottom)) {
	  sprintf(logbuf, "Error: Invalid --tail-pheno parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct == 1) {
	  tail_top = tail_bottom;
	} else {
	  if (scan_double(argv[cur_arg + 2], &tail_top)) {
	    sprintf(logbuf, "Error: Invalid --tail-pheno parameter '%s'.%s", argv[cur_arg + 2], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (tail_bottom > tail_top) {
	  sprintf(logbuf, "Error: Ltop cannot be larger than Hbottom for --tail-pheno.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	misc_flags |= MISC_TAIL_PHENO;
      } else if (!memcmp(argptr2, "hreads", 7)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	ii = atoi(argv[cur_arg + 1]);
	if (ii < 1) {
	  sprintf(logbuf, "Error: Invalid --threads parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (ii > MAX_THREADS) {
	  sprintf(logbuf, "Note: Reducing --threads parameter to %d.  (If this is not large enough,\nrecompile with a larger MAX_THREADS setting.)\n", MAX_THREADS);
	  logprintb();
	  ii = MAX_THREADS;
	}
	g_thread_ct = ii;
      } else if (!memcmp(argptr2, "ab", 3)) {
	logprint("Note: --tab flag deprecated.  Use '--recode tab ...'.\n");
	if (recode_modifier & RECODE_DELIMX) {
	  sprintf(logbuf, "Error: Multiple --recode delimiter modifiers.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	recode_modifier |= RECODE_TAB;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "ranspose", 9)) {
	logprint("Note: --transpose flag deprecated.  Use '--recode transpose ...'.\n");
	if (recode_modifier & RECODE_LGEN) {
	  sprintf(logbuf, "Error: --recode 'transpose' and 'lgen' modifiers cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	recode_modifier |= RECODE_TRANSPOSE;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "fam", 4)) {
	if (load_params || load_rare) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logprint("Error: --tfam filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(famname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TFAM;
      } else if (!memcmp(argptr2, "file", 5)) {
	if (load_params || (load_rare & (~LOAD_RARE_TFAM))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  jj = strlen(argv[cur_arg + 1]);
	  if (jj > FNAMESIZE - 6) {
	    logprint("Error: --tfile filename prefix too long.\n");
	    goto main_ret_OPEN_FAIL;
	  }
	  memcpy(memcpya(pedname, argv[cur_arg + 1], jj), ".tped", 6);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(memcpya(famname, argv[cur_arg + 1], jj), ".tfam", 6);
	  }
	} else {
	  memcpy(pedname, PROG_NAME_STR ".tped", 11);
	  if (!(load_rare & LOAD_RARE_TFAM)) {
	    memcpy(famname, PROG_NAME_STR ".tfam", 11);
	  }
	}
	load_rare |= LOAD_RARE_TRANSPOSE;
      } else if (!memcmp(argptr2, "ped", 4)) {
	if (load_params || (load_rare & (~(LOAD_RARE_TRANSPOSE | LOAD_RARE_TFAM)))) {
	  goto main_ret_INVALID_CMDLINE_4;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	jj = strlen(argv[cur_arg + 1]);
	if (jj > FNAMESIZE - 1) {
	  logprint("Error: --tped filename prefix too long.\n");
	  goto main_ret_OPEN_FAIL;
	}
	memcpy(pedname, argv[cur_arg + 1], jj + 1);
	load_rare |= LOAD_RARE_TPED;
      } else if (!memcmp(argptr2, "o", 2)) {
	if ((!all_words_zero(chrom_info.chrom_mask, CHROM_MASK_INITIAL_WORDS)) || chrom_info.incl_excl_name_stack) {
	  sprintf(logbuf, "Error: --to cannot be used with --autosome[-xy] or --[not-]chr.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_snp) {
	  sprintf(logbuf, "Error: --to cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
	  sprintf(logbuf, "Error: --to cannot be used with --snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "Error: --to cannot be used with a .cnv fileset.  Use --to-bp/-kb/-mb instead.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (alloc_string(&markername_to, argv[cur_arg + 1])) {
	  goto main_ret_NOMEM;
	}
      } else if ((!memcmp(argptr2, "o-bp", 5)) || (!memcmp(argptr2, "o-kb", 5)) || (!memcmp(argptr2, "o-mb", 5))) {
	if (markername_snp && (!(misc_flags & MISC_EXCLUDE_MARKERNAME_SNP))) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (snps_range_list.names) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --snps.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (markername_to) {
	  sprintf(logbuf, "Error: --to-bp/-kb/-mb cannot be used with --to.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((!markername_from) && (marker_pos_start == -1)) {
	  marker_pos_start = 0;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	cc = argptr2[2];
	if (cc == 'b') {
	  if (atoiz(argv[cur_arg + 1], &ii)) {
	    sprintf(logbuf, "Error: Invalid --to-bp parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	} else {
	  if (marker_pos_end != -1) {
	    logprint("Error: Multiple --to-bp/-kb/-mb values.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  if (scan_double(argv[cur_arg + 1], &dxx)) {
	    sprintf(logbuf, "Error: Invalid --to-kb/-mb parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  dxx *= (cc == 'k')? 1000 : 1000000;
	  if (dxx < 0) {
	    ii = 0;
	  } else if (dxx > 2147483647) {
	    ii = 0x7fffffff;
	  } else {
	    ii = (int32_t)(dxx + EPSILON);
	  }
	}
	if (ii < marker_pos_start) {
	  marker_pos_end = marker_pos_start;
	  marker_pos_start = ii;
	} else {
	  marker_pos_end = ii;
	}
      } else if (!memcmp(argptr2, "rend", 5)) {
	if ((!(calculation_type & CALC_MODEL)) || (model_modifier & MODEL_ASSOC)) {
	  sprintf(logbuf, "Error: --trend must be used with --model.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (model_modifier & MODEL_FISHER) {
	  sprintf(logbuf, "Error: --trend cannot be used with --model fisher.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (model_modifier & (MODEL_PDOM | MODEL_PREC | MODEL_PGEN)) {
	  sprintf(logbuf, "Error: --trend cannot be used with --model dom/rec/gen.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --trend flag deprecated.  Use '--model trend ...'.\n");
	model_modifier |= MODEL_PTREND | MODEL_TRENDONLY;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "hin", 4)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &thin_keep_prob)) {
	  sprintf(logbuf, "Error: Invalid --thin marker retention probability '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (thin_keep_prob < (0.5 / 4294967296.0)) {
	  sprintf(logbuf, "Error: --thin marker retention probability too small.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (thin_keep_prob >= (4294967295.5 / 4294967296.0)) {
	  sprintf(logbuf, "Error: --thin marker retention probability too large.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else if (!memcmp(argptr2, "ests", 5)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --tests must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((param_ct == 1) && (!strcmp(argv[cur_arg + 1], "all"))) {
	  glm_modifier |= GLM_TEST_ALL;
	} else {
	  if (glm_modifier & GLM_TEST_ALL) {
	    sprintf(logbuf, "Error: --test-all cannot be used with --tests.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  retval = parse_name_ranges(param_ct, '-', &(argv[cur_arg]), &tests_range_list, 1);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "est-all", 8)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --test-all must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --test-all flag deprecated.  Use '--tests all'.\n");
	glm_modifier |= GLM_TEST_ALL;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'u':
      if (!memcmp(argptr2, "nrelated-heritability", 22)) {
#ifdef NOLAPACK
        logprint("Error: --unrelated-heritability requires " PROG_NAME_CAPS " to be built with LAPACK.\n");
	goto main_ret_INVALID_CMDLINE;
#else
	UNSTABLE;
	if (rel_calc_type & REL_CALC_COV) {
	  sprintf(logbuf, "Error: --unrelated-heritability flag cannot coexist with a covariance\nmatrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (parallel_tot > 1) {
	  sprintf(logbuf, "Error: --parallel and --unrelated-heritability cannot be used together.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (rel_calc_type & REL_CALC_SINGLE_PREC) {
	  sprintf(logbuf, "Error: --unrelated-heritability flag cannot be used with a single-precision\nrelationship matrix calculation.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
	  if (calculation_type & CALC_REL_CUTOFF) {
	    sprintf(logbuf, "Error: --unrelated-heritability + --grm-gz/--grm-bin cannot be used with\n--rel-cutoff.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  } else if (!phenoname) {
	    sprintf(logbuf, "Error: --unrelated-heritability + --grm-gz/--grm-bin requires --pheno as well.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (!strcmp(argv[cur_arg + 1], "strict")) {
	    misc_flags |= MISC_UNRELATED_HERITABILITY_STRICT;
	    uii = 2;
	  } else {
	    uii = 1;
	  }
	  if (param_ct >= uii) {
	    if (scan_double(argv[cur_arg + uii], &unrelated_herit_tol)) {
	      sprintf(logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (unrelated_herit_tol <= 0.0) {
	      sprintf(logbuf, "Error: Invalid --unrelated-heritability EM tolerance parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    if (param_ct > uii) {
	      if (scan_double(argv[cur_arg + uii + 1], &unrelated_herit_covg)) {
		sprintf(logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.%s", argv[cur_arg + uii + 1], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if ((unrelated_herit_covg <= 0.0) || (unrelated_herit_covg > 1.0)) {
		sprintf(logbuf, "Error: Invalid --unrelated-heritability genomic covariance prior '%s'.%s", argv[cur_arg + uii + 1], errstr_append);
		goto main_ret_INVALID_CMDLINE_3;
	      }
	      if (param_ct == uii + 2) {
		if (scan_double(argv[cur_arg + uii + 2], &unrelated_herit_covr)) {
		  sprintf(logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.%s", argv[cur_arg + uii + 2], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
		if ((unrelated_herit_covr <= 0.0) || (unrelated_herit_covr > 1.0)) {
		  sprintf(logbuf, "Error: Invalid --unrelated-heritability residual covariance prior '%s'.%s", argv[cur_arg + uii + 2], errstr_append);
		  goto main_ret_INVALID_CMDLINE_3;
		}
	      } else {
		unrelated_herit_covr = 1.0 - unrelated_herit_covg;
	      }
	    }
	  }
	}
	calculation_type |= CALC_UNRELATED_HERITABILITY;
#endif
      } else if (!memcmp(argptr2, "pdate-alleles", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&update_alleles_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "pdate-chr", 10)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-chr cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!param_ct) {
	  update_map_modifier = 1;
	} else {
	  retval = alloc_2col(&update_chr, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-cm", 9)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-cm cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_map_modifier) {
	  logprint("Error: --update-map 'cm' modifier cannot be used with 'chr'.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
        if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!param_ct) {
	  update_map_modifier = 2;
	} else {
	  retval = alloc_2col(&update_cm, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-ids", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	retval = alloc_fname(&update_ids_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "pdate-map", 10)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-map cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_map_modifier) {
	  if (param_ct != 1) {
	    sprintf(logbuf, "Error: Multi-parameter --update-map cannot be used with deprecated\nparameter-free --update-%s.\n", (update_map_modifier == 1)? "chr" : "cm");
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  retval = alloc_2col((update_map_modifier == 1)? (&update_chr) : (&update_cm), &(argv[cur_arg + 1]), argptr, 1);
	  if (retval) {
	    goto main_ret_1;
	  }
	  sprintf(logbuf, "Note: --update-map [filename] + parameter-free --update-%s deprecated.  Use\n--update-%s [filename] instead.\n", (update_map_modifier == 1)? "chr" : "cm", (update_map_modifier == 1)? "chr" : "cm");
	  logprintb();
	  update_map_modifier = 0;
	} else {
	  retval = alloc_2col(&update_map, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-name", 11)) {
	if (cnv_calc_type & CNV_MAKE_MAP) {
	  sprintf(logbuf, "--update-name cannot be used with --cnv-make-map.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_chr) {
	  logprint("Error: --update-name cannot be used with --update-chr.\n");
	  goto main_ret_INVALID_CMDLINE;
	} else if (update_cm) {
	  logprint("Error: --update-name cannot be used with --update-cm.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 4)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!param_ct) {
	  if (!update_map) {
	    sprintf(logbuf, "Error: Deprecated parameter-free --update-name cannot be used without\n--update-map.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	  if ((update_map->colx != 2) || (update_map->colid != 1) || (update_map->skip) || (update_map->skipchar)) {
	    logprint("Error: Multi-parameter --update-map cannot be used with deprecated\nparameter-free --update-name.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  logprint("Note: --update-map [filename] + parameter-free --update-name deprecated.  Use\n--update-name [filename] instead.\n");
	  update_name = update_map;
	  update_map = NULL;
	} else {
	  if (update_map) {
	    // no point in explaining the deprecated exception to this in the
	    // error message
	    logprint("Error: --update-name cannot be used with --update-map.\n");
	    goto main_ret_INVALID_CMDLINE;
	  }
	  retval = alloc_2col(&update_name, &(argv[cur_arg + 1]), argptr, param_ct);
	  if (retval) {
	    goto main_ret_1;
	  }
	}
      } else if (!memcmp(argptr2, "pdate-parents", 14)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_ids_fname) {
	  logprint("Error: --update-parents cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&update_parents_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "pdate-sex", 10)) {
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (update_ids_fname) {
	  logprint("Error: --update-sex cannot be used with --update-ids.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	retval = alloc_fname(&update_sex_fname, argv[cur_arg + 1], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "nbounded", 9)) {
	if (!(calculation_type & CALC_GENOME)) {
	  sprintf(logbuf, "Error: --unbounded must be used with --genome.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	logprint("Note: --unbounded flag deprecated.  Use '--genome unbounded'.\n");
	genome_modifier |= GENOME_IBD_UNBOUNDED;
	goto main_param_zero;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'v':
      if (!memcmp(argptr2, "if", 3)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --vif must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &glm_vif_thresh)) {
	  sprintf(logbuf, "Error: Invalid --linear/--logistic VIF threshold '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (glm_vif_thresh < 1.0) {
	  sprintf(logbuf, "Error: --linear/--logistic VIF threshold '%s' too small (must be >= 1).%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'w':
      if (!memcmp(argptr2, "rite-snplist", 13)) {
	calculation_type |= CALC_WRITE_SNPLIST;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "indow", 6)) {
        if (!markername_snp) {
	  sprintf(logbuf, "Error: --window must be used with --snp or --exclude-snp.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (scan_double(argv[cur_arg + 1], &dxx)) {
	  sprintf(logbuf, "Error: Invalid --window parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        dxx *= 500;
	if (dxx < 1) {
	  sprintf(logbuf, "Error: --window parameter '%s' too small.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (dxx > 2147483647) {
	  snp_window_size = 0x7fffffff;
	} else {
	  snp_window_size = (int32_t)(dxx + EPSILON);
	}
      } else if (!memcmp(argptr2, "ithin", 6)) {
        if (loop_assoc_fname) {
	  sprintf(logbuf, "Error: --within cannot be used with --loop-assoc.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((calculation_type & CALC_FREQ) && (misc_flags & MISC_FREQ_COUNTS)) {
	  sprintf(logbuf, "Error: --within cannot be used with '--freq counts'.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	uii = 1;
	if (param_ct == 2) {
	  if ((strlen(argv[cur_arg + 1]) == 7) && (!memcmp(argv[cur_arg + 1], "keep-", 5)) && match_upper(&(argv[cur_arg + 1][5]), "NA")) {
	    uii = 2;
	  } else if ((strlen(argv[cur_arg + 2]) != 7) || memcmp(argv[cur_arg + 2], "keep-", 5) || (!match_upper(&(argv[cur_arg + 2][5]), "NA"))) {
            sprintf(logbuf, "Error: Invalid --within parameter sequence.%s", errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_LOAD_CLUSTER_KEEP_NA;
	}
	retval = alloc_fname(&cluster.fname, argv[cur_arg + uii], argptr, 0);
	if (retval) {
	  goto main_ret_1;
	}
      } else if (!memcmp(argptr2, "ith-phenotype", 14)) {
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --with-phenotype cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 2)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	for (uii = 1; uii <= param_ct; uii++) {
	  if (!strcmp(argv[cur_arg + uii], "no-parents")) {
	    write_covar_modifier |= WRITE_COVAR_NO_PARENTS;
	  } else if (!strcmp(argv[cur_arg + uii], "no-sex")) {
	    if (write_covar_modifier & WRITE_COVAR_FEMALE_2) {
	      sprintf(logbuf, "Error: --with-phenotype 'female-2' modifier cannot be used with 'no-sex'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    write_covar_modifier |= WRITE_COVAR_NO_SEX;
	  } else if (!strcmp(argv[cur_arg + uii], "female-2")) {
	    if (write_covar_modifier & WRITE_COVAR_NO_SEX) {
	      sprintf(logbuf, "Error: --with-phenotype 'female-2' modifier cannot be used with 'no-sex'.%s", errstr_append);
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    write_covar_modifier |= WRITE_COVAR_FEMALE_2;
	  } else {
	    sprintf(logbuf, "Error: Invalid --with-phenotype parameter '%s'.%s", argv[cur_arg + uii], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
	}
	write_covar_modifier |= WRITE_COVAR_PHENO;
      } else if (!memcmp(argptr2, "ith-reference", 14)) {
	if ((recode_modifier & RECODE_TYPEMASK) != RECODE_LGEN) {
	  logprint("Error: --with-reference must be used with --recode lgen.\n");
	  goto main_ret_INVALID_CMDLINE;
	}
	logprint("Note: --with-reference flag deprecated.  Use '--recode lgen-ref' instead.\n");
	recode_modifier += RECODE_LGEN_REF - RECODE_LGEN;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-covar", 11)) {
	if (calculation_type & (CALC_MAKE_BED | CALC_RECODE)) {
	  sprintf(logbuf, "Error: --write-covar cannot be used with --make-bed or --recode.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (!covar_fname) {
	  sprintf(logbuf, "Error: --write-covar cannot be used without --covar.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
        calculation_type |= CALC_WRITE_COVAR;
	goto main_param_zero;
      } else if (!memcmp(argptr2, "rite-cluster", 13)) {
	if (!cluster.fname) {
	  sprintf(logbuf, "Error: --write-cluster must be used with --within.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 0, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (param_ct) {
	  if (strcmp(argv[cur_arg + 1], "omit-unassigned")) {
	    sprintf(logbuf, "Error: Invalid --write-cluster parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	    goto main_ret_INVALID_CMDLINE_3;
	  }
          misc_flags |= MISC_WRITE_CLUSTER_OMIT_UNASSIGNED;
	}
        calculation_type |= CALC_WRITE_CLUSTER;
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    case 'x':
      if (!memcmp(argptr2, "chr-model", 10)) {
	if (!(calculation_type & CALC_GLM)) {
	  sprintf(logbuf, "Error: --xchr-model must be used with --linear or --logistic.%s", errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	} else if (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE)) {
	  sprintf(logbuf, "Error: --xchr-model cannot be used with --%s %s.%s", (glm_modifier & GLM_LOGISTIC)? "logistic" : "linear", (glm_modifier & GLM_GENOTYPIC)? "genotypic" : ((glm_modifier & GLM_HETHOM)? "hethom" : ((glm_modifier & GLM_DOMINANT)? "dominant" : "recessive")), errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if (enforce_param_ct_range(param_ct, argv[cur_arg], 1, 1)) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
	if ((argv[cur_arg + 1][1] != '\0') || (argv[cur_arg + 1][0] < '0') || (argv[cur_arg + 1][0] > '3')) {
	  sprintf(logbuf, "Error: Invalid --xchr-model parameter '%s'.%s", argv[cur_arg + 1], errstr_append);
	  goto main_ret_INVALID_CMDLINE_3;
	}
	glm_xchr_model = (uint32_t)(argv[cur_arg + 1][0] - '0');
      } else {
	goto main_ret_INVALID_CMDLINE_2;
      }
      break;

    default:
      goto main_ret_INVALID_CMDLINE_2;

    main_param_zero:
      if (param_ct) {
        sprintf(logbuf, "Error: --%s doesn't accept parameters.%s", argptr, errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
  } while ((++cur_flag) < flag_ct);
  if (!outname_end) {
    outname_end = &(outname[5]);
  }

  // command-line restrictions which don't play well with alphabetical order
  if (load_rare) {
    if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
      if ((!(calculation_type & (CALC_REL_CUTOFF | CALC_UNRELATED_HERITABILITY))) || (calculation_type & (~(CALC_REL_CUTOFF | CALC_RELATIONSHIP | CALC_UNRELATED_HERITABILITY)))) {
	if (load_rare == LOAD_RARE_GRM) {
	  sprintf(logbuf, "Error: --grm-gz currently must be used with --rel-cutoff (possibly combined\nwith --make-grm-gz/--make-grm-bin) or --unrelated-heritability.%s", errstr_append);
	} else {
	  sprintf(logbuf, "Error: --grm-bin currently must be used with --rel-cutoff (possibly combined\nwith --make-grm-gz/--make-grm-bin) or --unrelated-heritability.%s", errstr_append);
	}
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
    if (!mperm_val) {
      if ((cnv_calc_type & CNV_INDIV_PERM) && (!cnv_indiv_mperms)) {
	sprintf(logbuf, "Error: --cnv-indiv-perm requires a permutation count.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      } else if ((cnv_calc_type & CNV_TEST_REGION) && (!cnv_test_region_mperms)) {
	sprintf(logbuf, "Error: --cnv-test-region requires a permutation count.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
  }
  if ((cnv_intersect_filter_type & CNV_COUNT) && (!(cnv_calc_type & (CNV_INDIV_PERM | CNV_ENRICHMENT_TEST)))) {
    sprintf(logbuf, "Error: --cnv-count must be used with --cnv-indiv-perm or --cnv-enrichment-test.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (!phenoname) {
    if ((misc_flags & MISC_PRUNE) && (!(fam_cols & FAM_COL_6))) {
      sprintf(logbuf, "Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    } else if (pheno_modifier & PHENO_ALL) {
      sprintf(logbuf, "Error: --all-pheno must be used with --pheno.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  } else if (read_dists_fname && (!(calculation_type & (CALC_CLUSTER | CALC_IBS_TEST | CALC_GROUPDIST | CALC_NEIGHBOR | CALC_REGRESS_DISTANCE)))) {
    sprintf(logbuf, "Error: --read-dists cannot be used without --cluster, --ibs-test/--groupdist,\n--neighbour, or --regress-distance.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((cluster.ppc != 0.0) && (!read_genome_fname) && (calculation_type & (CALC_DISTANCE))) {
    sprintf(logbuf, "Error: --ppc cannot be used with --distance without --read-genome.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if ((calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR)) && (((calculation_type & CALC_DISTANCE) && (dist_calc_type & (DISTANCE_1_MINUS_IBS | DISTANCE_ALCT))) || (calculation_type & CALC_PLINK_DISTANCE_MATRIX)) && (!(calculation_type & CALC_GENOME))) {
    // actually allow this for now if --genome present, since it auto-clobbers
    // the wrong-unit distance matrix
    sprintf(logbuf, "Error: --cluster and --neighbour cannot be used with non-IBS distance matrix\ncalculations.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (update_map_modifier) {
    sprintf(logbuf, "Error: Deprecated parameter-free --update-%s cannot be used without\n--update-map.%s", (update_map_modifier == 1)? "chr" : "cm", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (update_chr && (((load_rare == LOAD_RARE_CNV) && (cnv_calc_type != CNV_WRITE)) || ((!load_rare) && (calculation_type != CALC_MAKE_BED)))) {
    sprintf(logbuf, "Error: --update-chr/--update-map must be used with --make-bed and no other\ncommands.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (flip_subset_fname && (load_rare || (calculation_type != CALC_MAKE_BED) || (min_maf != 0.0) || (max_maf != 0.5) || (hwe_thresh != 0.0))) {
    sprintf(logbuf, "Error: --flip-subset must be used with --flip, --make-bed, and no other\ncommands or MAF-based filters.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (calculation_type & CALC_RECODE) {
    if (recode_modifier & (RECODE_23 | RECODE_BEAGLE)) {
      if (chrom_info.species != SPECIES_HUMAN) {
	sprintf(logbuf, "Error: --recode %s can only be used on human data.\n", (recode_modifier & RECODE_23)? "23" : "beagle");
	goto main_ret_INVALID_CMDLINE_3;
      }
      if ((recode_modifier & RECODE_23) && ((misc_flags & (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS)) == MISC_ALLOW_EXTRA_CHROMS)) {
	logprint("Error: --allow-extra-chr requires the '0' modifier when used with --recode 23.\n");
	goto main_ret_INVALID_CMDLINE;
      } else if ((recode_modifier & RECODE_BEAGLE) && ((misc_flags & (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS)) == (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS))) {
        logprint("Error: --allow-extra-chr cannot have the '0' modifier when used with\n--recode beagle.\n");
	goto main_ret_INVALID_CMDLINE;
      }
    } else if ((recode_modifier & RECODE_BIMBAM) && ((misc_flags & (MISC_ALLOW_EXTRA_CHROMS | MISC_ZERO_EXTRA_CHROMS)) == MISC_ALLOW_EXTRA_CHROMS)) {
      logprint("Error: --allow-extra-chr requires the '0' modifier when used with\n--recode bimbam.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (sex_missing_pheno & MUST_HAVE_SEX) {
    if (load_rare & LOAD_RARE_CNV) {
      if (!(cnv_calc_type & CNV_WRITE)) {
        logprint("Error: --must-have-sex must be used with --cnv-write.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    } else {
      if (!(calculation_type & (CALC_MAKE_BED | CALC_RECODE))) {
        logprint("Error: --must-have-sex must be used with --make-bed or --recode.\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
  }
  if (cluster.qmatch_fname && (!cluster.qt_fname)) {
    sprintf(logbuf, "Error: --qt must be used with --qmatch.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  if (!cluster.fname) {
    if (mwithin_col && (!loop_assoc_fname)) {
      sprintf(logbuf, "Error: --mwithin must be used with --within.%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    } else if (cluster.keep_fname) {
      sprintf(logbuf, "Error: --keep-clusters must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.keep_flattened) {
      sprintf(logbuf, "Error: --keep-cluster-names must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.remove_fname) {
      sprintf(logbuf, "Error: --remove-clusters must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    } else if (cluster.remove_flattened) {
      sprintf(logbuf, "Error: --remove-cluster-names must be used with --within.\n");
      goto main_ret_INVALID_CMDLINE;
    }
  }

  // --from-bp/-kb/-mb without any --to/--to-bp/...: include to end of
  // chromosome
  if ((marker_pos_start != -1) && (!markername_to) && (marker_pos_end == -1)) {
    marker_pos_end = 0x7fffffff;
  }
  if (!chrom_flag_present) {
    fill_chrom_mask(&chrom_info);
  }
  if (((marker_pos_start != -1) && (!markername_to)) || ((marker_pos_end != -1) && (!markername_from))) {
    // require exactly one chromosome to be defined given --from-bp/--to-bp
    // without --from/--to
    uii = next_set(chrom_info.chrom_mask, 0, CHROM_MASK_INITIAL_WORDS * BITCT);
    if (uii == CHROM_MASK_INITIAL_WORDS * BITCT) {
      uii = 0;
    } else {
      uii = next_set(chrom_info.chrom_mask, uii + 1, CHROM_MASK_INITIAL_WORDS * BITCT);
    }
    if (((uii == CHROM_MASK_INITIAL_WORDS * BITCT) && chrom_info.incl_excl_name_stack) || ((uii != CHROM_MASK_INITIAL_WORDS * BITCT) && (uii || (!chrom_info.incl_excl_name_stack) || chrom_info.incl_excl_name_stack->next))) {
      sprintf(logbuf, "Error: --from-bp/-kb/-mb and --to-bp/-kb/-mb require a single chromosome to be\nidentified (either explicitly with --chr, or implicitly with --from/--to).%s", errstr_append);
      goto main_ret_INVALID_CMDLINE_3;
    }
  }

  if (calculation_type & CALC_MODEL) {
    if (!(model_modifier & (MODEL_ASSOC | MODEL_PDOM | MODEL_PREC | MODEL_PTREND))) {
      if (mtest_adjust) {
	sprintf(logbuf, "Error: In order to use --model with --adjust, you must include the 'trend',\n'trend-only', 'dom', or 'rec' modifier.%s", errstr_append);
	goto main_ret_INVALID_CMDLINE_3;
      }
    }
    if (model_cell_ct == -1) {
      model_cell_ct = (model_modifier & MODEL_FISHER)? 0 : 5;
    }
    if ((model_modifier & (MODEL_PERM | MODEL_MPERM | MODEL_GENEDROP)) == MODEL_GENEDROP) {
      model_modifier |= MODEL_PERM;
    }
  }
  if ((homozyg.modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE)) && (!(calculation_type & CALC_HOMOZYG))) {
    if (homozyg.overlap_min == 0.95) {
      sprintf(logbuf, "Error: --homozyg-group must be used with another --homozyg... flag.%s", errstr_append);
    } else {
      sprintf(logbuf, "Error: --homozyg-match must be used with another --homozyg... flag.%s", errstr_append);
    }
    goto main_ret_INVALID_CMDLINE_3;
  }

  if ((!calculation_type) && (!(load_rare & (LOAD_RARE_LGEN | LOAD_RARE_DUMMY | LOAD_RARE_SIMULATE | LOAD_RARE_TRANSPOSE_MASK | LOAD_RARE_23 | LOAD_RARE_CNV))) && (famname[0] || load_rare)) {
    goto main_ret_NULL_CALC;
  }
  if (!(load_params || load_rare)) {
    sprintf(logbuf, "Error: No input dataset.%s", errstr_append);
    goto main_ret_INVALID_CMDLINE_3;
  }
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  subst_argv = NULL;
  script_buf = NULL;
  rerun_buf = NULL;
  flag_buf = NULL;
  flag_map = NULL;
  if (!rseeds) {
    sfmt_init_gen_rand(&sfmt, (uint32_t)time(NULL));
  } else {
    if (rseed_ct == 1) {
      sfmt_init_gen_rand(&sfmt, rseeds[0]);
    } else {
      sfmt_init_by_array(&sfmt, rseeds, rseed_ct);
    }
    free(rseeds);
    rseeds = NULL;
  }

  // guarantee contiguous malloc space outside of main workspace
  bubble = (char*)malloc(NON_WKSPACE_MIN * sizeof(char));
  if (!bubble) {
    goto main_ret_NOMEM;
  }

  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  sztmp = sizeof(int64_t);
  sysctl(mib, 2, &llxx, &sztmp, NULL, 0);
  llxx /= 1048576;
#else
#if _WIN32
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#else
  llxx = ((size_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
  if (!llxx) {
    default_alloc_mb = WKSPACE_DEFAULT_MB;
  } else if (llxx < (WKSPACE_MIN_MB * 2)) {
    default_alloc_mb = WKSPACE_MIN_MB;
  } else {
    default_alloc_mb = llxx / 2;
  }
  if (!malloc_size_mb) {
    malloc_size_mb = default_alloc_mb;
  } else if (malloc_size_mb < WKSPACE_MIN_MB) {
    malloc_size_mb = WKSPACE_MIN_MB;
  }
  if (llxx) {
    sprintf(logbuf, "%" PRId64 " MB RAM detected; reserving %" PRIdPTR " MB for main workspace.\n", llxx, malloc_size_mb);
  } else {
    sprintf(logbuf, "Failed to calculate system memory.  Attempting to reserve %" PRIdPTR " MB.\n", malloc_size_mb);
  }
  logprintb();
  wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  while (!wkspace_ua) {
    malloc_size_mb = (malloc_size_mb * 3) / 4;
    if (malloc_size_mb < WKSPACE_MIN_MB) {
      malloc_size_mb = WKSPACE_MIN_MB;
    }
    wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (wkspace_ua) {
      sprintf(logbuf, "Allocated %" PRIdPTR " MB successfully, after larger attempt(s) failed.\n", malloc_size_mb);
      logprintb();
    } else if (malloc_size_mb == WKSPACE_MIN_MB) {
      goto main_ret_NOMEM;
    }
  }
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((uintptr_t)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = malloc_size_mb * 1048576 - (uintptr_t)(wkspace - wkspace_ua);
  free(bubble);

  tbuf[MAXLINELEN - 6] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  pigz_init(g_thread_ct);
  if (load_rare & (LOAD_RARE_GRM | LOAD_RARE_GRM_BIN)) {
    // --unrelated-heritability and --rel-cutoff batch mode special cases
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      retval = unrelated_herit_batch(load_rare & LOAD_RARE_GRM_BIN, pedname, phenoname, mpheno_col, phenoname_str, missing_pheno, (misc_flags / MISC_UNRELATED_HERITABILITY_STRICT) & 1, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr);
    } else {
#endif
      retval = rel_cutoff_batch(load_rare & LOAD_RARE_GRM_BIN, pedname, outname, outname_end, rel_cutoff, rel_calc_type);
#ifndef NOLAPACK
    }
#endif
  } else if (genname[0]) {
    if (calculation_type & (~(CALC_DISTANCE | CALC_REGRESS_DISTANCE))) {
      logprint("Error: Only --distance calculations are currently supported with --data.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
    } else {
      if (!missing_code) {
	missing_code = (char*)"NA";
      }
      retval = wdist_dosage(calculation_type, dist_calc_type, genname, samplename, outname, outname_end, missing_code, exponent, (misc_flags / MISC_MAF_SUCC) & 1, regress_iters, regress_d, g_thread_ct, parallel_idx, parallel_tot);
    }
  } else if (load_rare & LOAD_RARE_CNV) {
    retval = wdist_cnv(outname, outname_end, pedname, mapname, famname, phenoname, keepname, removename, filtername, misc_flags, update_chr, update_cm, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, filtervals_flattened, filter_binary, cnv_calc_type, cnv_min_seglen, cnv_max_seglen, cnv_min_score, cnv_max_score, cnv_min_sites, cnv_max_sites, cnv_intersect_filter_type, cnv_intersect_filter_fname, cnv_subset_fname, cnv_overlap_type, cnv_overlap_val, cnv_freq_type, cnv_freq_val, cnv_freq_val2, cnv_test_window, segment_modifier, segment_spanning_fname, cnv_indiv_mperms, cnv_test_mperms, cnv_test_region_mperms, cnv_enrichment_test_mperms, marker_pos_start, marker_pos_end, &chrom_info);
  } else if (load_rare & LOAD_RARE_GVAR) {
    retval = wdist_gvar(outname, outname_end, pedname, mapname, famname);
  } else {
  // famname[0] indicates binary vs. text
  // extractname, excludename, keepname, and removename indicate the presence
  // of their respective flags
  // filtername indicates existence of filter
  // freqname signals --read-freq
    // if (load_rare) {
    if (load_rare || (!famname[0])) {
      sptr = outname_end;
      if (calculation_type && (!(misc_flags & MISC_KEEP_AUTOCONV))) {
        sptr = memcpyb(sptr, "-temporary", 11);
      }
      uii = (sptr - outname);
      if (load_rare == LOAD_RARE_LGEN) {
        retval = lgen_to_bed(pedname, outname, sptr, missing_pheno, misc_flags, lgen_modifier, lgen_reference_fname, &chrom_info);
      } else if (load_rare & LOAD_RARE_TRANSPOSE_MASK) {
        retval = transposed_to_bed(pedname, famname, outname, sptr, misc_flags, &chrom_info);
      } else if (load_rare == LOAD_RARE_23) {
        retval = bed_from_23(pedname, outname, sptr, modifier_23, fid_23, iid_23, (pheno_23 == INFINITY)? ((double)missing_pheno) : pheno_23, paternal_id_23, maternal_id_23, convert_xy_23, &chrom_info);
      } else if (load_rare & LOAD_RARE_DUMMY) {
	retval = generate_dummy(outname, sptr, dummy_flags, dummy_marker_ct, dummy_indiv_ct, dummy_missing_geno, dummy_missing_pheno);
      } else if (load_rare & LOAD_RARE_SIMULATE) {
	retval = simulate_dataset(outname, sptr, simulate_flags, simulate_fname, simulate_cases, simulate_controls, simulate_prevalence, simulate_qt_indivs, simulate_missing, simulate_label);
	free(simulate_fname);
	simulate_fname = NULL;
	if (simulate_label) {
	  free(simulate_label);
	  simulate_label = NULL;
	}
      } else {
        retval = ped_to_bed(pedname, mapname, outname, sptr, fam_cols, (misc_flags / MISC_AFFECTION_01) & 1, missing_pheno, &chrom_info);
	fam_cols |= FAM_COL_1 | FAM_COL_34 | FAM_COL_5;
	if (!(fam_cols & FAM_COL_6)) {
          fam_cols |= FAM_COL_6;
	  missing_pheno = -9;
	}
      }
      if (retval || (!calculation_type)) {
	goto main_ret_2;
      }
      memcpy(memcpya(pedname, outname, uii), ".bed", 5);
      memcpy(memcpya(mapname, outname, uii), ".bim", 5);
      memcpy(memcpya(famname, outname, uii), ".fam", 5);
      if (calculation_type && (!(misc_flags & MISC_KEEP_AUTOCONV))) {
	if (push_ll_str(&file_delete_list, pedname) || push_ll_str(&file_delete_list, mapname) || push_ll_str(&file_delete_list, famname)) {
	  goto main_ret_NOMEM;
	}
      }
      *outname_end = '\0';
    }
    if (ibc_type == 3) { // todo: make this less ugly
      ibc_type = 0;
    } else if (!ibc_type) {
      ibc_type = 1;
    }
    retval = wdist(outname, outname_end, pedname, mapname, famname, phenoname, extractname, excludename, keepname, removename, keepfamname, removefamname, filtername, freqname, read_dists_fname, read_dists_id_fname, evecname, mergename1, mergename2, mergename3, makepheno_str, phenoname_str, a1alleles, a2alleles, recode_allele_name, covar_fname, set_fname, subset_fname, update_alleles_fname, read_genome_fname, update_chr, update_cm, update_map, update_name, update_ids_fname, update_parents_fname, update_sex_fname, loop_assoc_fname, flip_fname, flip_subset_fname, filtervals_flattened, condition_mname, condition_fname, thin_keep_prob, min_bp_space, mfilter_col, filter_binary, fam_cols, missing_pheno, output_missing_pheno, mpheno_col, pheno_modifier, &chrom_info, exponent, min_maf, max_maf, geno_thresh, mind_thresh, hwe_thresh, rel_cutoff, tail_bottom, tail_top, misc_flags, calculation_type, rel_calc_type, dist_calc_type, groupdist_iters, groupdist_d, regress_iters, regress_d, regress_rel_iters, regress_rel_d, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr, ibc_type, parallel_idx, parallel_tot, ppc_gap, sex_missing_pheno, genome_modifier, genome_min_pi_hat, genome_max_pi_hat, &homozyg, &cluster, neighbor_n1, neighbor_n2, ld_window_size, ld_window_kb, ld_window_incr, ld_last_param, regress_pcs_modifier, max_pcs, recode_modifier, allelexxxx, merge_type, indiv_sort, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, &snps_range_list, covar_modifier, &covar_range_list, write_covar_modifier, write_covar_dummy_max_categories, mwithin_col, model_modifier, (uint32_t)model_cell_ct, model_mperm_val, glm_modifier, glm_vif_thresh, glm_xchr_model, glm_mperm_val, &parameters_range_list, &tests_range_list, ci_size, pfilter, mtest_adjust, adjust_lambda, gxe_mcovar, aperm_min, aperm_max, aperm_alpha, aperm_beta, aperm_init_interval, aperm_interval_slope, mperm_save, ibs_test_perms, perm_batch_size, lasso_h2, lasso_minlambda, &file_delete_list);
  }
 main_ret_2:
  free(wkspace_ua);
  while (0) {
  main_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  main_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  main_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  main_ret_INVALID_CMDLINE_2:
    invalid_arg(argv[cur_arg]);
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_4:
    sprintf(logbuf, "Error: --%s conflicts with another input flag.%s", argptr, errstr_append);
  main_ret_INVALID_CMDLINE_3:
    logprintb();
  main_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_NULL_CALC:
    logprint(notestr_null_calc);
    fputs(cmdline_format_str, stdout);
    fputs(notestr_null_calc2, stdout);
    retval = RET_NULL_CALC;
#ifdef STABLE_BUILD
    break;
  main_unstable_disabled:
    logprint("Error: This flag's implementation is unfinished or unstable.  If you wish to\ntest it, use the latest development build.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
#endif
  }
 main_ret_1:
  fclose_cond(scriptfile);
  dispmsg(retval);
  free_cond(subst_argv);
  free_cond(script_buf);
  free_cond(rerun_buf);
  free_cond(flag_buf);
  free_cond(flag_map);
  free_cond(makepheno_str);
  free_cond(phenoname_str);
  free_cond(a1alleles);
  free_cond(a2alleles);
  free_cond(filtervals_flattened);
  free_cond(evecname);
  free_cond(filtername);
  free_cond(read_dists_fname);
  free_cond(read_dists_id_fname);
  free_cond(freqname);
  free_cond(extractname);
  free_cond(excludename);
  free_cond(keepname);
  free_cond(removename);
  free_cond(keepfamname);
  free_cond(removefamname);
  free_cond(phenoname);
  free_cond(recode_allele_name);
  free_cond(markername_from);
  free_cond(markername_to);
  free_cond(markername_snp);
  free_range_list(&snps_range_list);
  free_range_list(&covar_range_list);
  free_range_list(&parameters_range_list);
  free_range_list(&tests_range_list);
  free_cond(lgen_reference_fname);
  free_cond(covar_fname);
  free_cond(set_fname);
  free_cond(subset_fname);
  free_cond(update_alleles_fname);
  free_cond(update_chr);
  free_cond(update_cm);
  free_cond(update_map);
  free_cond(update_name);
  free_cond(update_ids_fname);
  free_cond(update_parents_fname);
  free_cond(update_sex_fname);
  free_cond(loop_assoc_fname);
  free_cond(flip_fname);
  free_cond(flip_subset_fname);
  free_cond(read_genome_fname);
  free_cond(cluster.fname);
  free_cond(cluster.match_fname);
  free_cond(cluster.match_missing_str);
  free_cond(cluster.match_type_fname);
  free_cond(cluster.qmatch_fname);
  free_cond(cluster.qmatch_missing_str);
  free_cond(cluster.qt_fname);
  free_cond(cluster.keep_fname);
  free_cond(cluster.remove_fname);
  free_cond(cluster.keep_flattened);
  free_cond(cluster.remove_flattened);
  free_cond(rseeds);
  free_cond(simulate_fname);
  free_cond(simulate_label);
  free_cond(cnv_intersect_filter_fname);
  free_cond(cnv_subset_fname);
  free_cond(segment_spanning_fname);
  free_cond(fid_23);
  free_cond(iid_23);
  free_cond(paternal_id_23);
  free_cond(maternal_id_23);
  free_cond(convert_xy_23);
  free_cond(condition_mname);
  free_cond(condition_fname);
  if (file_delete_list) {
    do {
      ll_str_ptr = file_delete_list->next;
      unlink(file_delete_list->ss);
      free(file_delete_list);
      file_delete_list = ll_str_ptr;
    } while (file_delete_list);
  }
  forget_extra_chrom_names(&chrom_info);
  if (chrom_info.incl_excl_name_stack) {
    do {
      ll_str_ptr = chrom_info.incl_excl_name_stack->next;
      free(chrom_info.incl_excl_name_stack);
      chrom_info.incl_excl_name_stack = ll_str_ptr;
    } while (chrom_info.incl_excl_name_stack);
  }
  if (logfile) {
    if (!log_failed) {
      logstr("\nEnd time: ");
      time(&rawtime);
      logstr(ctime(&rawtime));
      if (fclose(logfile)) {
	fputs("Error: Failed to finish writing to log.\n", stdout);
      }
    } else {
      fclose(logfile);
    }
  }
  return retval;
}
