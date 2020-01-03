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


// The key ideas behind this calculator's design are:
//
// 1. Incremental processing of SNPs.  Each element A_{jk} of a distance or
// relationship matrix is a sum of M terms, one for each SNP, multiplied by a
// missingness correction at the end.  We can walk through the SNPs
// sequentially without keeping much in memory beyond partial sums;
// conveniently, this plays well with the decision made by the PLINK team a few
// years ago to switch to SNP-major binary files.
//
// 2. Multiplexing of markers using bitwise operations.  For instance, there
// are only seven possible ways SNP i can affect the relationship matrix entry
// between samples j and k:
//    a. j and k are both homozygous for the rare allele
//    b. one is homozygous rare, one is heterozygous
//    c. one is homozygous rare, one is homozygous common
//    d. both are heterozygous
//    e. one is heterozygous, one is homozygous common
//    f. both are homozygous common
//    g. one or both has missing genotype data
// Seven cases can be distinguished by three bits, so one can compose a
// sequence of bitwise operations that maps a pair of padded 2-bit PLINK
// genotypes to seven different 3-bit values according to this breakdown.
// On 64-bit machines, this allows integer operations to act on 20 markers
// simultaneously.  (There's space for a 21st, but we currently choose not to
// use it.)
//
// 3. Lookup tables describing the effect of 5-7 markers at a time on a
// distance or relationship, rather than just one.  For relationship matrices,
// idea #2 allows computation of a single 64-bit integer where bits 0-2
// describe the relationship on marker #0, bits 3-5 describe the relationship
// on marker #1, ..., all the way up to bits 57-59 describing the relationship
// on marker #19.  We then want to perform the update
//    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_19(x_19)
// where the x_i's are bit trios, and the f_i's map them to floating point
// terms.  We could do this with 20 table lookups and floating point additions.
// Or, we could structure this as
//    A_{jk} := A_{jk} + f_{0-4}(x_{0-4}) + ... + f_{15-19}(x_{15-19})
// where x_{0-4} denotes the lowest order *15* bits, and f_{0-4} maps them
// directly to f_0(x_0) + f_1(x_1) + f_2(x_2) + f_3(x_3) + f_4(x_4); similarly
// for f_{5-9}, f_{10-14}, and f_{15-19}.  This requires precomputation of four
// lookup tables of size 2^15, total size 1 MB (which fits comfortably in a
// typical L2/L3 cache these days), and licenses the use of four table lookups
// and adds instead of twenty.
//
// 4. Bitslicing algorithms for especially fast calculation of unweighted
// distances and SNP covariances.  Zero-exponent distance matrices and IBS
// matrices are special cases, reducing to Hamming weight computations plus a
// bit of dressing to handle missing markers.  Hamming weight computation on
// x86 CPUs has been studied extensively by others; a good reference as of this
// writing is
//    http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html .
// We use a variant of the Kim Walisch/Cedric Lauradoux bitslicing algorithm
// discussed there (with most 64-bit integer operations replaced by SSE2
// instructions), which runs several times faster than our corresponding
// nonzero exponent distance matrix computation.
//
// We also reduce SNP covariance calculation (used in LD-based pruning) to a
// few Hamming weights.  (This can also be done for covariances between
// samples, but only when there are no missing markers, so we do not include
// an implementation of that.)
//
// 5. Splitting the distance/relationship matrix into pieces of roughly equal
// size and assigning one thread to each piece.  This is an "embarrassingly
// parallel" problem with no need for careful synchronization.  Cluster
// computing is supported in essentially the same manner.
//
//
//
// In the end, we can't get around the O(MN^2) nature of these calculations,
// but we believe we have beaten down the leading constant by a large enough
// factor to meaningfully help researchers.

#include "plink_common.h"

#include "plink_calc.h"
#include "plink_cluster.h"
#include "plink_matrix.h"
#include "plink_misc.h"
#include "plink_stats.h"
#include "pigz.h"

// number of different types of jackknife values to precompute (x^2, x, y, xy)
#define JACKKNIFE_VALS_REL 5

// Number of snp-major .bed lines to read at once for distance calc if exponent
// is zero.  Currently assumed to be a multiple of 192, and no larger than
// 1920, by the popcount_..._multiword functions.  (The optimal value depends
// on both system-specific properties such as cache sizes, as well as the
// number of samples in the current calculation, so in principle it's best to
// select this value at runtime.  But 960 usually works well in practice in my
// experience.)
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)

// Must be multiple of 384, no larger than 3840.
#define GENOME_MULTIPLEX 1152
#define GENOME_MULTIPLEX2 (GENOME_MULTIPLEX * 2)

#define MAX_EM_ACCEL 100.0

void rel_init(Rel_info* relip) {
  relip->modifier = 0;
  relip->regress_rel_iters = ITERS_DEFAULT;
  relip->regress_rel_d = 0;
  relip->unrelated_herit_tol = 0.0000001;
  relip->unrelated_herit_covg = 0.45;
  relip->unrelated_herit_covr = 0.55;
  relip->cutoff = 0.025;
  relip->ibc_type = 0;
  relip->pc_ct = 20;
  relip->pca_cluster_names_flattened = nullptr;
  relip->pca_clusters_fname = nullptr;
}

void rel_cleanup(Rel_info* relip) {
  free_cond(relip->pca_cluster_names_flattened);
  free_cond(relip->pca_clusters_fname);
}

void update_rel_ibc(double* rel_ibc, uintptr_t* geno, double* set_allele_freqs, double* main_weights, int32_t ibc_type, uint32_t sample_ct, uint32_t window_size) {
  // first calculate weight array, then loop
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  double twt;
  double* wtptr;
  double mult = 1.0;
  uintptr_t ulii;
  double weights[BITCT * 12];
  double* weights1 = &(weights[64]);
  double* weights2 = &(weights[128]);
  double* weights3 = &(weights[256]);
  double* weights4 = &(weights[320]);
#ifdef __LP64__
  double* weights5 = &(weights[384]);
  double* weights6 = &(weights[448]);
  double* weights7 = &(weights[512]);
  double* weights8 = &(weights[640]);
  double* weights9 = &(weights[704]);
#endif
  double wtarr[BITCT2 * 5];
  double *wptr = weights;
  fill_double_zero(BITCT2 * 5, wtarr);
  for (uii = 0; uii < window_size; uii++) {
    if ((set_allele_freqs[uii] != 0.0) && (set_allele_freqs[uii] < (1.0 - EPSILON))) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[uii * 8] = 2;
          wtarr[uii * 8 + 2] = 2.0 - 1.0 / (2 * set_allele_freqs[uii] * (1.0 - set_allele_freqs[uii]));
          wtarr[uii * 8 + 3] = 2;
        } else {
          twt = 2 * set_allele_freqs[uii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - set_allele_freqs[uii]));
          }
          wtarr[uii * 8] = twt * twt * mult;
          wtarr[uii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[uii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - set_allele_freqs[uii];
        mult = 1 / (set_allele_freqs[uii] * twt);
        wtarr[uii * 8] = 1.0 + set_allele_freqs[uii] * set_allele_freqs[uii] * mult;
        wtarr[uii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * set_allele_freqs[uii];
          wtarr[uii * 8] = twt * twt;
          wtarr[uii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[uii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  // infinities are useful here for calling out inaccurate zero MAFs,
	  // but those markers should just be automatically filtered out
	  // instead?
	  wtarr[uii * 8 + 2] = INFINITY;
          if (set_allele_freqs[uii] == 0.0) {
            wtarr[uii * 8] = 0;
            wtarr[uii * 8 + 3] = INFINITY;
          } else {
            wtarr[uii * 8] = INFINITY;
            wtarr[uii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[uii * 8] = 1;
          wtarr[uii * 8 + 2] = -INFINITY;
          wtarr[uii * 8 + 3] = 1;
        }
      } else {
        if (set_allele_freqs[uii] == 0.0) {
          wtarr[uii * 8] = 1;
          wtarr[uii * 8 + 3] = INFINITY;
        } else {
          wtarr[uii * 8] = INFINITY;
          wtarr[uii * 8 + 3] = 1;
        }
      }
    }
    if (main_weights) {
      wtarr[uii * 8] *= main_weights[uii];
      wtarr[uii * 8 + 2] *= main_weights[uii];
      wtarr[uii * 8 + 3] *= main_weights[uii];
    }
  }
  for (ukk = 0; ukk < (BITCT * 5) / 32; ukk++) {
    wtptr = &(wtarr[16 * ukk]);
#ifdef __LP64__
    if ((ukk == 2) || (ukk == 7)) {
      for (uii = 0; uii < 8; uii++) {
	twt = wtptr[uii + 8];
	for (ujj = 0; ujj < 8; ujj++) {
	  *wptr++ = twt + wtptr[ujj];
	}
	wptr = &(wptr[8]);
      }
    } else {
      for (uii = 0; uii < 8; uii++) {
	twt = wtptr[uii + 8];
	for (ujj = 0; ujj < 8; ujj++) {
	  *wptr++ = twt + wtptr[ujj];
	}
      }
    }
#else
    if (ukk == 2) {
      for (uii = 0; uii < 8; uii++) {
	twt = wtptr[uii + 8];
	for (ujj = 0; ujj < 8; ujj++) {
	  *wptr++ = twt + wtptr[ujj];
	}
	wptr = &(wptr[8]);
      }
    } else {
      for (uii = 0; uii < 8; uii++) {
	twt = wtptr[uii + 8];
	for (ujj = 0; ujj < 8; ujj++) {
	  *wptr++ = twt + wtptr[ujj];
	}
      }
    }
#endif
  }
  for (umm = 0; umm < sample_ct; umm++) {
    ulii = *geno++;
#ifdef __LP64__
    *rel_ibc += weights9[ulii >> 57] + weights8[(ulii >> 51) & 63] + weights7[(ulii >> 44) & 127] + weights6[(ulii >> 38) & 63] + weights5[(ulii >> 32) & 63] + weights4[(ulii >> 25) & 63] + weights3[(ulii >> 19) & 63] + weights2[(ulii >> 12) & 127] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 25] + weights3[(ulii >> 19) & 63] + weights2[(ulii >> 12) & 127] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

void fill_subset_weights(double* subset_weights, double* main_weights) {
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  double wtarr[MULTIPLEX_DIST_EXP / 2];
  double* wt;
#ifdef __LP64__
  double twt[5];
  double twtf;
  __m128d* swpairs = (__m128d*)subset_weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
#else
  uint32_t upp;
  uint32_t uqq;
  double twt[7];
#endif
  memcpy(wtarr, main_weights, (MULTIPLEX_DIST_EXP / 2) * sizeof(double));
  for (uoo = 0; uoo < 2; uoo++) {
    wt = &(wtarr[7 * uoo]);
#ifdef __LP64__
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(wt[0] * 2, wt[0]);
#endif
    twt[0] = 0;
    for (uii = 0; uii < 4; uii += 1) {
      // turning the uii == 2 case into a memcpy doesn't actually seem to speed
      // this up
      if (uii & 1) {
	twt[0] += wt[6];
      }
      twt[1] = twt[0];
      for (ujj = 0; ujj < 4; ujj += 1) {
	if (ujj & 1) {
	  twt[1] += wt[5];
	}
	twt[2] = twt[1];
	for (ukk = 0; ukk < 4; ukk += 1) {
	  if (ukk & 1) {
	    twt[2] += wt[4];
	  }
	  twt[3] = twt[2];
	  for (umm = 0; umm < 4; umm++) {
	    if (umm & 1) {
	      twt[3] += wt[3];
	    }
	    twt[4] = twt[3];
	    for (unn = 0; unn < 4; unn++) {
	      if (unn & 1) {
		twt[4] += wt[2];
	      }
#ifdef __LP64__
	      twtf = twt[4];
	      vpen = _mm_set1_pd(twtf);
	      *swpairs++ = _mm_add_pd(vpen, vfinal1);
	      *swpairs++ = _mm_add_pd(vpen, vfinal2);
	      twtf += wt[1];
	      vpen = _mm_set1_pd(twtf);
	      *swpairs++ = _mm_add_pd(vpen, vfinal1);
	      *swpairs++ = _mm_add_pd(vpen, vfinal2);
	      *swpairs = *(swpairs - 2);
	      swpairs++;
	      *swpairs = *(swpairs - 2);
	      swpairs++;
	      vpen = _mm_set1_pd(twtf + wt[1]);
	      *swpairs++ = _mm_add_pd(vpen, vfinal1);
	      *swpairs++ = _mm_add_pd(vpen, vfinal2);
#else
              twt[5] = twt[4];
              for (upp = 0; upp < 4; upp++) {
                if (upp & 1) {
                  twt[5] += wt[1];
                }
                twt[6] = twt[5];
                for (uqq = 0; uqq < 4; uqq++) {
                  if (uqq & 1) {
                    twt[6] += wt[0];
                  }
                  *subset_weights++ = twt[6];
                }
              }
#endif
	    }
          }
	}
      }
    }
  }
#ifdef __LP64__
  for (uoo = 0; uoo < 3; uoo++) {
    wt = &(wtarr[14 + 6 * uoo]);
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(2 * wt[0], wt[0]);
    twt[0] = 0;
    for (uii = 0; uii < 4; uii += 1) {
      if (uii & 1) {
        twt[0] += wt[5];
      }
      twt[1] = twt[0];
      for (ujj = 0; ujj < 4; ujj += 1) {
        if (ujj & 1) {
          twt[1] += wt[4];
        }
        twt[2] = twt[1];
	for (ukk = 0; ukk < 4; ukk += 1) {
          if (ukk & 1) {
            twt[2] += wt[3];
          }
          twt[3] = twt[2];
          for (umm = 0; umm < 4; umm++) {
	    if (umm & 1) {
	      twt[3] += wt[2];
	    }
	    twtf = twt[3];
	    vpen = _mm_set1_pd(twtf);
	    *swpairs++ = _mm_add_pd(vpen, vfinal1);
	    *swpairs++ = _mm_add_pd(vpen, vfinal2);
	    twtf += wt[1];
	    vpen = _mm_set1_pd(twtf);
	    *swpairs++ = _mm_add_pd(vpen, vfinal1);
	    *swpairs++ = _mm_add_pd(vpen, vfinal2);
	    *swpairs = *(swpairs - 2);
	    swpairs++;
	    *swpairs = *(swpairs - 2);
	    swpairs++;
	    vpen = _mm_set1_pd(twtf + wt[1]);
	    *swpairs++ = _mm_add_pd(vpen, vfinal1);
	    *swpairs++ = _mm_add_pd(vpen, vfinal2);
          }
	}
      }
    }
  }
#endif
}

void fill_subset_weights_r(double* subset_weights, double* set_allele_freqs, double* main_weights, uint32_t var_std) {
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  // 20 markers to process in quintuplets, for 64-bit; 10, for 32-bit.
  // Each quintuplet of markers requires 40 wtarr entries, and induces
  // 2^15 writes to subset_weights[].
  double wtarr_raw[BITCT2 * 5 + 1];
  double* wtarr = wtarr_raw;
  double twt;
  double twt2;
  double twt3;
  double twt4;
  double* wtptr;
  double mean;
  double mean_m1;
  double mean_m2;
  double mult = 1.0;
  double aux;
#ifdef __LP64__
  __m128d* swpairs = (__m128d*)subset_weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
  __m128d vfinal3;
  __m128d vfinal4;
#else
  uint32_t uoo;
#endif
  if (((uintptr_t)wtarr) & 15) {
    // force 16-byte alignment; can't do this at compile-time since stack
    // pointer has no 16-byte align guarantee.
    // yes, this assumes doubles are 8 bytes.
    wtarr++;
  }
  for (uii = 0; uii < MULTIPLEX_REL / 3; uii += 1) {
    if (((set_allele_freqs[uii] != 0.0) && (set_allele_freqs[uii] < (1.0 - EPSILON))) || (!var_std)) {
      if (set_allele_freqs[uii] < 0.5) {
	mean = 2 * set_allele_freqs[uii];
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * (1.0 - set_allele_freqs[uii]));
        }
        aux = mean * mult;
	wtarr[uii * 8] = mean * aux;
        wtarr[uii * 8 + 1] = 0;
	wtarr[uii * 8 + 2] = mean_m1 * aux;
	wtarr[uii * 8 + 3] = mean_m2 * aux;
	wtarr[uii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[uii * 8 + 5] = mean_m2 * mean_m1 * mult;
	wtarr[uii * 8 + 6] = mean_m2 * mean_m2 * mult;
      } else {
	mean = 2 * (1.0 - set_allele_freqs[uii]);
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * set_allele_freqs[uii]);
        }
        aux = mean_m2 * mult;
	wtarr[uii * 8] = mean_m2 * aux;
        wtarr[uii * 8 + 1] = 0;
	wtarr[uii * 8 + 2] = mean_m1 * aux;
	wtarr[uii * 8 + 3] = mean * aux;
	wtarr[uii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[uii * 8 + 5] = mean_m1 * mean * mult;
	wtarr[uii * 8 + 6] = mean * mean * mult;
      }
    } else {
      if (set_allele_freqs[uii] == 0.0) {
        wtarr[uii * 8] = 0;
        wtarr[uii * 8 + 1] = 0;
        wtarr[uii * 8 + 2] = -1;
        wtarr[uii * 8 + 3] = -2;
        wtarr[uii * 8 + 4] = INFINITY;
        wtarr[uii * 8 + 5] = INFINITY;
        wtarr[uii * 8 + 6] = INFINITY;
      } else {
        wtarr[uii * 8] = INFINITY;
        wtarr[uii * 8 + 1] = 0;
        wtarr[uii * 8 + 2] = INFINITY;
        wtarr[uii * 8 + 3] = -2;
        wtarr[uii * 8 + 4] = INFINITY;
        wtarr[uii * 8 + 5] = -1;
        wtarr[uii * 8 + 6] = 0;
      }
    }
    if (main_weights) {
      for (ujj = 0; ujj < 7; ujj++) {
	wtarr[uii * 8 + ujj] *= main_weights[uii];
      }
    }
    wtarr[uii * 8 + 7] = 0;
  }
  for (unn = 0; unn < BITCT / 16; unn++) {
    wtptr = &(wtarr[40 * unn]);
#ifdef __LP64__
    vfinal1 = _mm_load_pd(wtptr);
    vfinal2 = _mm_load_pd(&(wtptr[2]));
    vfinal3 = _mm_load_pd(&(wtptr[4]));
    vfinal4 = _mm_load_pd(&(wtptr[6]));
#endif
    for (uii = 0; uii < 8; uii++) {
      twt = wtptr[uii + 32];
      for (ujj = 0; ujj < 8; ujj++) {
        twt2 = twt + wtptr[ujj + 24];
        for (ukk = 0; ukk < 8; ukk++) {
          twt3 = twt2 + wtptr[ukk + 16];
          for (umm = 0; umm < 8; umm++) {
            twt4 = twt3 + wtptr[umm + 8];
#ifdef __LP64__
            vpen = _mm_set1_pd(twt4);
            *swpairs++ = _mm_add_pd(vpen, vfinal1);
            *swpairs++ = _mm_add_pd(vpen, vfinal2);
            *swpairs++ = _mm_add_pd(vpen, vfinal3);
            *swpairs++ = _mm_add_pd(vpen, vfinal4);
#else
            for (uoo = 0; uoo < 8; uoo++) {
              *subset_weights++ = twt4 + wtptr[uoo];
            }
#endif
          }
        }
      }
    }
  }
}

void exclude_multi(uintptr_t* exclude_arr, int32_t* new_excl, uint32_t unfiltered_sample_ct, uintptr_t* exclude_ct_ptr) {
  uint32_t exclude_ct = *exclude_ct_ptr;
  int32_t* new_excl_end = &(new_excl[unfiltered_sample_ct - exclude_ct]);
  uint32_t sample_uidx = 0;
  uint32_t sample_uidx_stop;
  do {
    sample_uidx = next_unset_unsafe(exclude_arr, sample_uidx);
    sample_uidx_stop = next_set(exclude_arr, sample_uidx, unfiltered_sample_ct);
    do {
      if (*new_excl++ == -1) {
        SET_BIT(sample_uidx, exclude_arr);
	exclude_ct++;
      }
    } while (++sample_uidx < sample_uidx_stop);
  } while (new_excl < new_excl_end);
  *exclude_ct_ptr = exclude_ct;
}

void collapse_copy_phenod(double* target, double* pheno_d, uintptr_t* sample_exclude, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct) {
  uintptr_t sample_uidx = 0;
  double* target_end = &(target[sample_ct]);
  uintptr_t delta;
  do {
    sample_uidx = next_unset_ul_unsafe(sample_exclude, sample_uidx);
    delta = next_set_ul(sample_exclude, sample_uidx, unfiltered_sample_ct) - sample_uidx;
    memcpy(target, &(pheno_d[sample_uidx]), delta * sizeof(double));
    target = &(target[delta]);
    sample_uidx += delta;
  } while (target < target_end);
}

#ifdef __LP64__
// XOR + mask variants of vectorized Lauradoux/Walisch popcount.  (See
// popcount_vecs() in plink_common.c for basic documentation.)
// Note that the size of the popcounted buffer is a hardcoded constant
// (specifically, (MULTIPLEX_DIST / BITCT) * 16 bytes).  The current code
// assumes (MULTIPLEX_DIST / BITCT) is a multiple of 3, and no greater than 30.
static inline uint32_t popcount_xor_1mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** maskp) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i count1, count2, half1, half2;
  __univec acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  do {
    count1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    count2 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    half1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (xor2 < xor2_end);
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
  // are guaranteed to be <= 120, thus adding two together does not overflow
  // 255.
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  return ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
}

static inline uint32_t popcount_xor_2mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** mask1p, __m128i* mask2) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i count1, count2, half1, half2;
  __univec acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  do {
    count1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    count2 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    half1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (xor2 < xor2_end);
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  return ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
static inline uint32_t popcount_xor_1mask_multiword(uintptr_t** xor1p, uintptr_t* xor2, uintptr_t** maskp) {
  uintptr_t* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  uint32_t bit_count = 0;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  do {
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    bit_count += (tmp_stor * 0x01010101) >> 24;
  } while (xor2 < xor2_end);
  return bit_count;
}

static inline uint32_t popcount_xor_2mask_multiword(uintptr_t** xor1p, uintptr_t* xor2, uintptr_t** mask1p, uintptr_t* mask2) {
  uintptr_t* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  uint32_t bit_count = 0;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  do {
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    loader -= ((loader >> 1) & FIVEMASK);
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    ulii += (loader & 0x33333333) + ((loader >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    loader -= ((loader >> 1) & FIVEMASK);
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    ulii += (loader & 0x33333333) + ((loader >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    bit_count += (tmp_stor * 0x01010101) >> 24;
  } while (xor2 < xor2_end);
  return bit_count;
}
#endif

// ----- multithread globals -----
double* g_rel_dists = nullptr;
uint32_t* g_sample_missing_unwt = nullptr;
uint32_t* g_missing_dbl_excluded = nullptr;
double* g_dists = nullptr;

static uint32_t g_thread_start[MAX_THREADS_P1];
static int32_t* g_idists;
static uintptr_t* g_pheno_nm = nullptr;
static uintptr_t* g_pheno_c = nullptr;
static unsigned char* g_geno = nullptr;
static double* g_subset_weights;
static uint32_t* g_subset_weights_i;
static double g_reg_tot_xy;
static double g_reg_tot_x;
static double g_reg_tot_y;
static double g_reg_tot_xx;
static double g_reg_tot_yy;
static uint32_t g_ctrl_ct;
static uint32_t g_case_ct;
static uintptr_t g_sample_ct;
static uintptr_t g_jackknife_iters;
static uint32_t g_jackknife_d;
static double g_calc_result[MAX_THREADS_P1][9];
static uintptr_t* g_masks;
static uintptr_t* g_mmasks;
static uint32_t* g_missing_tot_weights;
static uint32_t* g_sample_missing;
static double* g_jackknife_precomp = nullptr;
static uint32_t* g_genome_main = nullptr;
static uintptr_t g_marker_window[GENOME_MULTIPLEX * 2];
static double* g_pheno_packed;

// interleaved: first word of *every* permutation, then second word, etc.
static uintptr_t* g_perm_rows;

static uintptr_t* g_perm_col_buf;
static double* g_ibs_test_partial_sums;
static double* g_perm_results;
static uintptr_t g_perm_ct;
static double g_half_marker_ct_recip;
static uint32_t g_load_dists;
static unsigned char* g_generic_buf;

void ibs_test_init_col_buf(uintptr_t row_idx, uintptr_t perm_ct, uintptr_t* perm_rows, uintptr_t* perm_col_buf) {
  uintptr_t perm_idx = 0;
  uintptr_t block_size = BITCT;
  uintptr_t rshift_ct = row_idx & (BITCT - 1);
  uintptr_t* gpr_ptr;
  uintptr_t ulii;
  uintptr_t block_pos;
  gpr_ptr = &(perm_rows[(row_idx / BITCT) * perm_ct]);
  do {
    if (perm_idx + BITCT > perm_ct) {
      block_size = perm_ct - perm_idx;
    }
    ulii = 0;
    block_pos = 0;
    do {
      ulii |= (((*gpr_ptr++) >> rshift_ct) & ONELU) << block_pos;
    } while (++block_pos < block_size);
    perm_idx += block_size;
    *perm_col_buf++ = ulii;
  } while (perm_idx < perm_ct);
}

double fill_psbuf(uintptr_t block_size, double half_marker_ct_recip, uint32_t load_dists, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* dists, uintptr_t* col_uidxp, double* psbuf, double* ssq0p) {
  // also updates total sum and sums of squares
  double tot = 0.0;
  uintptr_t col_idx = 0;
  uintptr_t col_uidx = *col_uidxp;
  uint32_t sub_block_size = 8;
  uint32_t sub_block_idx;
  double subtot;
  uintptr_t ulii;
  double increment[8];
  double ssq[2];
  double dxx;
  ssq[0] = 0.0;
  ssq[1] = 0.0;
  do {
    if (col_idx + 8 > block_size) {
      sub_block_size = block_size - col_idx;
    }
    sub_block_idx = 0;
    subtot = 0.0;
    do {
      next_set_ul_unsafe_ck(pheno_nm, &col_uidx);
      if (load_dists) {
	dxx = dists[col_uidx];
      } else {
	dxx = 1.0 - dists[col_uidx] * half_marker_ct_recip;
      }
      increment[sub_block_idx] = subtot - dxx;
      subtot += dxx;
      ssq[IS_SET(pheno_c, col_uidx)] += dxx * dxx;
      col_uidx++;
    } while (++sub_block_idx < sub_block_size);
    tot += subtot;
    while (sub_block_idx < 8) {
      increment[sub_block_idx++] = subtot;
    }
    ulii = 0;
    dxx = subtot;
    goto fill_psbuf_loop_start;

    do {
      dxx += increment[CTZLU(++ulii)];
    fill_psbuf_loop_start:
      *psbuf++ = dxx;
    } while (ulii < 255);
    col_idx += sub_block_size;
  } while (col_idx < block_size);
  *col_uidxp = col_uidx;
  ssq0p[0] += ssq[0];
  ssq0p[1] += ssq[1];
  return tot;
}

void ibs_test_process_perms(uintptr_t* perm_row_start, uintptr_t perm_ct, uint32_t sub_block_ct, double block_tot, double* psbuf, uintptr_t* perm_col_buf, double* perm_results) {
  uintptr_t perm_idx = 0;
  uintptr_t block_size = BITCT;
  double dxx;
  uintptr_t block_pos;
  uintptr_t col_bits;
  uintptr_t ulii;
  uint32_t sub_block_idx;
  do {
    col_bits = *perm_col_buf++;
    if (perm_idx + BITCT > perm_ct) {
      block_size = perm_ct - perm_idx;
    }
    block_pos = 0;
    if (sub_block_ct == BITCT / 8) {
      do {
	sub_block_idx = 0;
	ulii = *perm_row_start++;
#ifdef __LP64__
	dxx = psbuf[(uint8_t)ulii] + psbuf[256 + ((uint8_t)(ulii >> 8))] + psbuf[512 + ((uint8_t)(ulii >> 16))] + psbuf[768 + ((uint8_t)(ulii >> 24))] + psbuf[1024 + ((uint8_t)(ulii >> 32))] + psbuf[1280 + ((uint8_t)(ulii >> 40))] + psbuf[1536 + ((uint8_t)(ulii >> 48))] + psbuf[1792 + (ulii >> 56)];
#else
        dxx = psbuf[(uint8_t)ulii] + psbuf[256 + ((uint8_t)(ulii >> 8))] + psbuf[512 + ((uint8_t)(ulii >> 16))] + psbuf[768 + (ulii >> 24)];
#endif
	if (col_bits & 1) {
	  perm_results[2 * block_pos + 1] += dxx;
	} else {
	  perm_results[2 * block_pos] += dxx;
	  perm_results[2 * block_pos + 1] += (block_tot - dxx);
	}
	col_bits >>= 1;
      } while (++block_pos < block_size);
    } else {
      do {
	sub_block_idx = 0;
	ulii = *perm_row_start++;
	dxx = psbuf[(uint8_t)ulii];
	while (++sub_block_idx < sub_block_ct) {
	  dxx += psbuf[256 * sub_block_idx + ((uint8_t)(ulii >> (8 * sub_block_idx)))];
	}
	if (col_bits & 1) {
	  perm_results[2 * block_pos + 1] += dxx;
	} else {
	  perm_results[2 * block_pos] += dxx;
	  perm_results[2 * block_pos + 1] += (block_tot - dxx);
	}
	col_bits >>= 1;
      } while (++block_pos < block_size);
    }
    perm_idx += block_size;
    perm_results = &(perm_results[2 * block_size]);
  } while (perm_idx < perm_ct);
}

void ibs_test_range(uint32_t tidx, uintptr_t* perm_col_buf, double* perm_results) {
  // (11-bit chunks were tested and found wanting.)

  // 256 possible bytes * (BITCT / 8) bytes per word
  double* psptr = &(g_ibs_test_partial_sums[tidx * (32 * BITCT)]);
  double dist_tot = 0.0;
  uintptr_t row_uidx = 0;
  uintptr_t pct = 0;
  uintptr_t pct_div = 1 + ((g_thread_start[1] * (g_thread_start[1] - 1)) / 100);
  uintptr_t pct_next = pct_div;
  uintptr_t perm_ct = g_perm_ct;
  uintptr_t row_bound1 = g_thread_start[tidx];
  uintptr_t row_bound2 = g_thread_start[tidx + 1];
  uintptr_t* pheno_nm = g_pheno_nm;
  uintptr_t* pheno_c = g_pheno_c;
  uintptr_t* perm_rows = g_perm_rows;
  double* dists = g_dists;
  double half_marker_ct_recip = g_half_marker_ct_recip;
  uint32_t load_dists = g_load_dists;
  double ssq[3];
  double* dptr;
  double block_tot;
  uintptr_t row_idx;
  uintptr_t col_idx;
  uintptr_t col_uidx;
  uintptr_t block_size;
  uintptr_t ulii;
  uint32_t row_set;
  ssq[0] = 0.0;
  ssq[1] = 0.0;
  ssq[2] = 0.0;
  for (row_idx = 0; row_idx < row_bound1; row_uidx++, row_idx++) {
    next_set_ul_unsafe_ck(pheno_nm, &row_uidx);
  }
  for (; row_idx < row_bound2; row_uidx++, row_idx++) {
    next_set_ul_unsafe_ck(pheno_nm, &row_uidx);
    dptr = &(dists[(row_uidx * (row_uidx - 1)) / 2]);
    col_idx = 0;
    col_uidx = 0;
    row_set = IS_SET(pheno_c, row_uidx);
    ibs_test_init_col_buf(row_idx, perm_ct, perm_rows, perm_col_buf);
    do {
      if (col_idx + BITCT > row_idx) {
	block_size = row_idx - col_idx;
      } else {
	block_size = BITCT;
      }
      block_tot = fill_psbuf(block_size, half_marker_ct_recip, load_dists, pheno_nm, pheno_c, dptr, &col_uidx, psptr, &(ssq[row_set]));
      dist_tot += block_tot;
      ibs_test_process_perms(&(perm_rows[(col_idx / BITCT) * perm_ct]), perm_ct, (block_size + 7) / 8, block_tot, psptr, perm_col_buf, perm_results);
      col_idx += block_size;
    } while (col_idx < row_idx);
    if (!tidx) {
      // technically should change other triangular pct loops to this as well,
      // to guard against int64 overflow with 400m+ people...
      ulii = row_idx * (row_idx + 1);
      if (ulii >= pct_next) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
	}
	pct = ulii / pct_div;
	printf("\b\b%" PRIuPTR "%%", pct);
	fflush(stdout);
	pct_next = pct_div * (pct + 1);
      }
    }
  }
  g_calc_result[tidx][0] = dist_tot;
  g_calc_result[tidx][1] = ssq[0];
  g_calc_result[tidx][2] = ssq[1];
  g_calc_result[tidx][3] = ssq[2];
}

THREAD_RET_TYPE ibs_test_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t perm_ctcl = (g_perm_ct + (CACHELINE * 8 - 1)) / (CACHELINE * 8);
  uintptr_t perm_ctcld = (g_perm_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  ibs_test_range((uint32_t)tidx, &(g_perm_col_buf[tidx * perm_ctcl * (CACHELINE / sizeof(intptr_t))]), &(g_perm_results[tidx * perm_ctcld * 2 * CACHELINE_DBL]));
  THREAD_RETURN;
}

void incr_dists_i(uint32_t* idists, uintptr_t* geno, uintptr_t* masks, uint32_t start_idx, uint32_t end_idx) {
#ifdef __LP64__
  __m128i* glptr;
  __m128i* glptr2;
  __m128i* mptr;
  __m128i* mcptr_start;
  uintptr_t* lptr;
#else
  uintptr_t* glptr;
  uintptr_t* glptr2;
  uintptr_t* mptr;
  uintptr_t* mcptr_start;
#endif
  uint32_t uii;
  int32_t jj;
  uintptr_t mask_fixed;
  for (uii = start_idx; uii < end_idx; uii++) {
    jj = uii * (MULTIPLEX_2DIST / BITCT);
#ifdef __LP64__
    glptr = (__m128i*)geno;
    glptr2 = (__m128i*)(&(geno[jj]));
    lptr = &(masks[jj]);
    mcptr_start = (__m128i*)lptr;
    mask_fixed = *lptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *lptr++;
    }
    mptr = (__m128i*)masks;
#else
    glptr = geno;
    glptr2 = &(geno[jj]);
    mcptr_start = &(masks[jj]);
    mptr = mcptr_start;
    mask_fixed = *mptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *mptr++;
    }
    mptr = masks;
#endif
    if (~mask_fixed) {
      while (glptr < glptr2) {
	*idists += popcount_xor_2mask_multiword(&glptr, glptr2, &mptr, mcptr_start);
	idists++;
      }
    } else {
      while (glptr < glptr2) {
	*idists += popcount_xor_1mask_multiword(&glptr, glptr2, &mptr);
	idists++;
      }
    }
  }
}

void incr_wt_dist_missing(uint32_t* mtw, uint32_t* subset_weights_i, uintptr_t* mmasks, uint32_t start_idx, uint32_t end_idx) {
  uintptr_t* glptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  for (uii = start_idx; uii < end_idx; uii++) {
    glptr = mmasks;
    ulii = mmasks[uii];
    if (ulii) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++) & ulii;
        while (uljj) {
          mtw[ujj] += subset_weights_i[CTZLU(uljj)];
          uljj &= uljj - 1;
        }
      }
    }
    mtw = &(mtw[uii]);
  }
}

void incr_dists_rm(uint32_t* idists, uintptr_t* mmasks, uint32_t start_idx, uint32_t end_idx) {
  // count missing intersection
  uintptr_t* mlptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  for (uii = start_idx; uii < end_idx; uii++) {
    mlptr = mmasks;
    ulii = mmasks[uii];
    if (ulii) {
      for (ujj = 0; ujj < uii; ujj++) {
        uljj = (*mlptr++) & ulii;
	if (uljj) {
	  idists[ujj] += popcount_long(uljj);
	}
      }
    }
    idists = &(idists[uii]);
  }
}

THREAD_RET_TYPE calc_ibs_thread(void* arg) {
  uintptr_t tidx = (intptr_t)arg;
  uintptr_t ulii = g_thread_start[tidx];
  uintptr_t uljj = g_thread_start[0];
  uintptr_t offset = (((uint64_t)ulii) * (ulii - 1) - ((uint64_t)uljj) * (uljj - 1)) / 2;
  uint32_t* weighted_missing_ptr = g_missing_tot_weights;
  uint32_t* flat_missing_ptr = g_missing_dbl_excluded;
  uint32_t* idists_ptr = (uint32_t*)(&(g_idists[offset]));
  uintptr_t* geno_ptr = (uintptr_t*)g_geno;
  uintptr_t* masks_ptr = g_masks;
  uintptr_t* mmasks_ptr = g_mmasks;
  uint32_t end_idx = g_thread_start[tidx + 1];
  uint32_t is_last_block;
  if (weighted_missing_ptr) {
    weighted_missing_ptr = &(weighted_missing_ptr[offset]);
  }
  if (flat_missing_ptr) {
    flat_missing_ptr = &(flat_missing_ptr[offset]);
  }
  while (1) {
    is_last_block = g_is_last_thread_block;
    if (weighted_missing_ptr) {
      // g_subset_weights_i moves around
      incr_wt_dist_missing(weighted_missing_ptr, g_subset_weights_i, mmasks_ptr, ulii, end_idx);
    }
    if (flat_missing_ptr) {
      incr_dists_rm(flat_missing_ptr, mmasks_ptr, ulii, end_idx);
    }
    if (is_last_block || ((g_thread_spawn_ct % (MULTIPLEX_DIST / BITCT)) == ((MULTIPLEX_DIST / BITCT) - 1))) {
      incr_dists_i(idists_ptr, geno_ptr, masks_ptr, ulii, end_idx);
    }
    if ((!tidx) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

void incr_genome(uint32_t* genome_main, uintptr_t* geno, uintptr_t* masks, uintptr_t sample_ct, uint32_t start_idx, uint32_t end_idx) {
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i xor_buf[GENOME_MULTIPLEX / BITCT];
  __m128i* xor_buf_end = &(xor_buf[GENOME_MULTIPLEX / BITCT]);
  __m128i* maskptr;
  __m128i* maskptr_fixed;
  __m128i* maskptr_fixed_tmp;
  __m128i* xor_ptr;
  __m128i* glptr_fixed_tmp;
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i count_ibs1;
  __m128i count_ibs0;
  __m128i count2_ibs1;
  __m128i count2_ibs0;
  __univec acc_ibs1;
  __univec acc_ibs0;
  uintptr_t* lptr;
  __m128i* glptr;
  __m128i* glptr_fixed;
  __m128i* glptr_end;
#else
  uintptr_t* glptr;
  uintptr_t* glptr_fixed;
  uintptr_t* glptr_end;
  uintptr_t* maskptr;
  uintptr_t* maskptr_fixed;
  uintptr_t* maskptr_fixed_tmp;
  uintptr_t* glptr_fixed_tmp;
  uintptr_t xor_buf[GENOME_MULTIPLEX2 / BITCT];
  uintptr_t* xor_buf_end = &(xor_buf[GENOME_MULTIPLEX2 / BITCT]);
  uintptr_t* xor_ptr;
  uint32_t bit_count_ibs1 = 0;
  uint32_t bit_count_ibs0 = 0;
  uintptr_t bitfield_ibs1;
  uintptr_t bitfield_ibs0;
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t tmp_stor_ibs1;
  uintptr_t tmp_stor_ibs0;
#endif
  uintptr_t* glptr_back;
  uintptr_t ibs_incr;
  uint32_t uii;
  uint32_t ujj;
  int32_t offset;
  uintptr_t uland;
  uintptr_t ulval;
  uintptr_t next_ppc_marker_hybrid;
  uintptr_t mask_fixed_test;
  uintptr_t* marker_window_ptr;
  int32_t lowct2 = g_ctrl_ct * 2;
  int32_t highct2 = g_case_ct * 2;
#ifdef __LP64__
  glptr_end = (__m128i*)(&(geno[sample_ct * (GENOME_MULTIPLEX2 / BITCT)]));
#else
  glptr_end = &(geno[sample_ct * (GENOME_MULTIPLEX2 / BITCT)]);
#endif
  for (uii = start_idx; uii < end_idx; uii++) {
    ujj = uii * (GENOME_MULTIPLEX2 / BITCT);
#ifdef __LP64__
    glptr_fixed = (__m128i*)(&(geno[ujj]));
    glptr = (__m128i*)(&(geno[ujj + (GENOME_MULTIPLEX2 / BITCT)]));
    lptr = &(masks[ujj]);
    maskptr = (__m128i*)(&(masks[ujj + (GENOME_MULTIPLEX2 / BITCT)]));
    maskptr_fixed = (__m128i*)lptr;
    mask_fixed_test = *lptr++;
    for (ujj = 0; ujj < GENOME_MULTIPLEX2 / BITCT - 1; ujj++) {
      mask_fixed_test &= *lptr++;
    }
#else
    glptr_fixed = &(geno[ujj]);
    glptr = &(geno[ujj + (GENOME_MULTIPLEX2 / BITCT)]);
    maskptr_fixed = &(masks[ujj]);
    maskptr = maskptr_fixed;
    mask_fixed_test = *maskptr++;
    for (ujj = 0; ujj < GENOME_MULTIPLEX2 / BITCT - 1; ujj++) {
      mask_fixed_test &= *maskptr++;
    }
#endif
    if (~mask_fixed_test) {
      while (glptr < glptr_end) {
	xor_ptr = xor_buf;
	glptr_back = (uintptr_t*)glptr;
	glptr_fixed_tmp = glptr_fixed;
	maskptr_fixed_tmp = maskptr_fixed;
#ifdef __LP64__
	acc_ibs1.vi = _mm_setzero_si128();
	acc_ibs0.vi = _mm_setzero_si128();
	do {
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
	  count_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count_ibs0;

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
	  count2_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count2_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count2_ibs0;

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

          count_ibs1 = _mm_add_epi64(_mm_and_si128(count_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count_ibs1, 2), m2));
          count_ibs0 = _mm_add_epi64(_mm_and_si128(count_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count_ibs0, 2), m2));
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_add_epi64(_mm_and_si128(count2_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs1, 2), m2)));
          count_ibs0 = _mm_add_epi64(count_ibs0, _mm_add_epi64(_mm_and_si128(count2_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs0, 2), m2)));
          acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_add_epi64(_mm_and_si128(count_ibs1, m4), _mm_and_si128(_mm_srli_epi64(count_ibs1, 4), m4)));
          acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_add_epi64(_mm_and_si128(count_ibs0, m4), _mm_and_si128(_mm_srli_epi64(count_ibs0, 4), m4)));
	} while (xor_ptr < xor_buf_end);
#if GENOME_MULTIPLEX > 1920
	acc_ibs1.vi = _mm_add_epi64(_mm_and_si128(acc_ibs1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs1.vi, 8), m8));
	acc_ibs0.vi = _mm_add_epi64(_mm_and_si128(acc_ibs0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs0.vi, 8), m8));
#else
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 8)), m8);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 8)), m8);
#endif
	*genome_main += ((acc_ibs1.u8[0] + acc_ibs1.u8[1]) * 0x1000100010001LLU) >> 48;
	genome_main++;
	*genome_main += ((acc_ibs0.u8[0] + acc_ibs0.u8[1]) * 0x1000100010001LLU) >> 48;
#else
        bit_count_ibs1 = 0;
	bit_count_ibs0 = 0;
	do {
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 = (loader ^ loader2) & FIVEMASK;
	  bitfield_ibs0 = (loader & loader2) & FIVEMASK;
	  *xor_ptr++ = bitfield_ibs0;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
          bitfield_ibs1 = (bitfield_ibs1 & 0x33333333) + ((bitfield_ibs1 >> 2) & 0x33333333);
	  bitfield_ibs0 = (bitfield_ibs0 & 0x33333333) + ((bitfield_ibs0 >> 2) & 0x33333333);
	  tmp_stor_ibs1 = (bitfield_ibs1 + (bitfield_ibs1 >> 4)) & 0x0f0f0f0f;
	  tmp_stor_ibs0 = (bitfield_ibs0 + (bitfield_ibs0 >> 4)) & 0x0f0f0f0f;

          bit_count_ibs1 += (tmp_stor_ibs1 * 0x01010101) >> 24;
	  bit_count_ibs0 += (tmp_stor_ibs0 * 0x01010101) >> 24;
	} while (xor_ptr < xor_buf_end);
	*genome_main += bit_count_ibs1;
	genome_main++;
	*genome_main += bit_count_ibs0;
#endif
	genome_main++;
	next_ppc_marker_hybrid = *genome_main - lowct2;
	if (next_ppc_marker_hybrid < GENOME_MULTIPLEX2) {
	  ibs_incr = 0; // hethet low-order, ibs0 high-order

          // This PPC test, rather than the IBS matrix, is now the limiting
	  // step of --genome if the markers aren't very dense.
	  //
          // I've taken a few "maintenance nightmare" liberties with this code
          // to speed it up, such as using a single lookup table that stores
          // values in two different forms (distinguished by the high bit of
          // the value), and using gotos for non-error-conditions, since the
          // loop is quite short.  In the long run, we want to implement
	  // support for better IBD estimation; see e.g. Browning B.L.,
	  // Browning S.R. "A Fast, Powerful Method for Detecting Identity by
	  // Descent", which discusses some weaknesses of PLINK --genome.
	  do {
	    offset = next_ppc_marker_hybrid / BITCT;
	    marker_window_ptr = &(g_marker_window[offset * BITCT]);
	    next_ppc_marker_hybrid = ~ZEROLU << (next_ppc_marker_hybrid & (BITCT - 1));
	  incr_genome_2mask_loop:
	    uland = glptr_back[offset] & (((uintptr_t*)glptr_fixed)[offset]);
	    // het is represented as 11, so
	    //   (uland & (uland << 1)) & 0xaaaaaaaaaaaaaaaa
	    // stores whether a particular marker is a hethet hit in the
	    // corresponding odd bit.
	    //
	    // homozygotes are represented as 01 and 10, so
	    //   (ulxor & (ulxor >> 1)) & 0x5555555555555555
	    // stores whether a particular marker is a hom1-hom2 hit in the
	    // corresponding even bit.  (het-missing pairs also set that bit,
	    // but the masking filters that out.)
	    //
	    // ~ZEROLU << xx masks out the bottom xx bits.
	    ulval = (((uland & (uland << 1)) & AAAAMASK) | (((uintptr_t*)xor_buf)[offset]));
	    do {
	      ulval &= next_ppc_marker_hybrid;
	      if (ulval) {
		ujj = CTZLU(ulval);
		next_ppc_marker_hybrid = marker_window_ptr[ujj];
		ibs_incr += (ONELU << ((ujj & 1) * BITCT2));
	      } else if (offset < ((GENOME_MULTIPLEX2 - BITCT) / BITCT)) {
		offset++;
		next_ppc_marker_hybrid = ~ZEROLU;
		marker_window_ptr = &(marker_window_ptr[BITCT]);
		goto incr_genome_2mask_loop;
	      } else {
		*genome_main = highct2;
		goto incr_genome_2mask_exit;
	      }
	    } while (next_ppc_marker_hybrid & (ONELU << (BITCT - 1)));
	  } while (next_ppc_marker_hybrid < GENOME_MULTIPLEX2);
	  *genome_main = next_ppc_marker_hybrid + lowct2;
        incr_genome_2mask_exit:
	  genome_main++;
          *genome_main += ibs_incr & ((~ZEROLU) >> BITCT2);
	  genome_main++;
	  *genome_main += ibs_incr >> BITCT2;
	  genome_main++;
	} else {
	  genome_main = &(genome_main[3]);
	}
      }
    } else {
      while (glptr < glptr_end) {
	xor_ptr = xor_buf;
	glptr_back = (uintptr_t*)glptr;
	glptr_fixed_tmp = glptr_fixed;
#ifdef __LP64__
	acc_ibs1.vi = _mm_setzero_si128();
	acc_ibs0.vi = _mm_setzero_si128();
	do {
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
	  count_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count_ibs0;
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
	  count2_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count2_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count2_ibs0;
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

          count_ibs1 = _mm_add_epi64(_mm_and_si128(count_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count_ibs1, 2), m2));
          count_ibs0 = _mm_add_epi64(_mm_and_si128(count_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count_ibs0, 2), m2));
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_add_epi64(_mm_and_si128(count2_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs1, 2), m2)));
          count_ibs0 = _mm_add_epi64(count_ibs0, _mm_add_epi64(_mm_and_si128(count2_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs0, 2), m2)));
          acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_add_epi64(_mm_and_si128(count_ibs1, m4), _mm_and_si128(_mm_srli_epi64(count_ibs1, 4), m4)));
          acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_add_epi64(_mm_and_si128(count_ibs0, m4), _mm_and_si128(_mm_srli_epi64(count_ibs0, 4), m4)));
	} while (xor_ptr < xor_buf_end);
#if GENOME_MULTIPLEX > 1920
	acc_ibs1.vi = _mm_add_epi64(_mm_and_si128(acc_ibs1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs1.vi, 8), m8));
	acc_ibs0.vi = _mm_add_epi64(_mm_and_si128(acc_ibs0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs0.vi, 8), m8));
#else
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 8)), m8);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 8)), m8);
#endif
        *genome_main += ((acc_ibs1.u8[0] + acc_ibs1.u8[1]) * 0x1000100010001LLU) >> 48;
	genome_main++;
        *genome_main += ((acc_ibs0.u8[0] + acc_ibs0.u8[1]) * 0x1000100010001LLU) >> 48;
#else
        bit_count_ibs1 = 0;
	bit_count_ibs0 = 0;
	do {
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 = (loader ^ loader2) & FIVEMASK;
	  bitfield_ibs0 = (loader & loader2) & FIVEMASK;
	  *xor_ptr++ = bitfield_ibs0;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
          bitfield_ibs1 = (bitfield_ibs1 & 0x33333333) + ((bitfield_ibs1 >> 2) & 0x33333333);
	  bitfield_ibs0 = (bitfield_ibs0 & 0x33333333) + ((bitfield_ibs0 >> 2) & 0x33333333);
	  tmp_stor_ibs1 = (bitfield_ibs1 + (bitfield_ibs1 >> 4)) & 0x0f0f0f0f;
	  tmp_stor_ibs0 = (bitfield_ibs0 + (bitfield_ibs0 >> 4)) & 0x0f0f0f0f;

          bit_count_ibs1 += (tmp_stor_ibs1 * 0x01010101) >> 24;
	  bit_count_ibs0 += (tmp_stor_ibs0 * 0x01010101) >> 24;
	} while (xor_ptr < xor_buf_end);
	*genome_main += bit_count_ibs1;
	genome_main++;
	*genome_main += bit_count_ibs0;
#endif
	genome_main++;
	next_ppc_marker_hybrid = *genome_main - lowct2;
	if (next_ppc_marker_hybrid < GENOME_MULTIPLEX2) {
	  ibs_incr = 0; // hethet low-order, ibs0 high-order
	  do {
	    offset = next_ppc_marker_hybrid / BITCT;
	    marker_window_ptr = &(g_marker_window[offset * BITCT]);
	    next_ppc_marker_hybrid = (~ZEROLU) << (next_ppc_marker_hybrid & (BITCT - 1));
	  incr_genome_1mask_loop:
	    uland = glptr_back[offset] & (((uintptr_t*)glptr_fixed)[offset]);
	    ulval = ((uland & (uland << 1)) & AAAAMASK) | (((uintptr_t*)xor_buf)[offset]);
	    do {
	      ulval &= next_ppc_marker_hybrid;
	      if (ulval) {
		ujj = CTZLU(ulval);
		next_ppc_marker_hybrid = marker_window_ptr[ujj];
		ibs_incr += (ONELU << ((ujj & 1) * BITCT2));
	      } else if (offset < ((GENOME_MULTIPLEX2 - BITCT) / BITCT)) {
		offset++;
		next_ppc_marker_hybrid = ~ZEROLU;
		marker_window_ptr = &(marker_window_ptr[BITCT]);
		goto incr_genome_1mask_loop;
	      } else {
		*genome_main = highct2;
		goto incr_genome_1mask_exit;
	      }
	    } while (next_ppc_marker_hybrid & (ONELU << (BITCT - 1)));
	  } while (next_ppc_marker_hybrid < GENOME_MULTIPLEX2);
	  *genome_main = next_ppc_marker_hybrid + lowct2;
	incr_genome_1mask_exit:
	  genome_main++;
	  *genome_main += ibs_incr & ((~ZEROLU) >> BITCT2);
	  genome_main++;
	  *genome_main += ibs_incr >> BITCT2;
	  genome_main++;
	} else {
	  genome_main = &(genome_main[3]);
	}
      }
    }
  }
}

void incr_dists_rm_inv(uint32_t* idists, uintptr_t* mmasks, uintptr_t sample_ct_m1, uint32_t start_idx, uint32_t end_idx) {
  // inverted loops for --genome --parallel
  uintptr_t* glptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  for (uii = start_idx; uii < end_idx; uii++) {
    ulii = mmasks[uii];
    ukk = sample_ct_m1 - uii;
    if (ulii) {
      glptr = &(mmasks[uii + 1]);
      for (ujj = 0; ujj < ukk; ujj++) {
        uljj = (*glptr++) & ulii;
	if (uljj) {
	  idists[ujj] += popcount_long(uljj);
	}
      }
    }
    idists = &(idists[ukk]);
  }
}

THREAD_RET_TYPE calc_genome_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t sample_ct_m1 = sample_ct - 1;
  uintptr_t ulii = g_thread_start[tidx];
  uintptr_t uljj = g_thread_start[0];
  // this is different from the regular offset because incr_dists_rm_inv() has
  // custom arithmetic
  uintptr_t offsetm = ((uint64_t)sample_ct) * (ulii - uljj) - ((((uint64_t)ulii) * (ulii + 1) - ((uint64_t)uljj) * (uljj + 1)) / 2);
  uintptr_t offset = offsetm * 5;
  uint32_t* missing_ptr = &(g_missing_dbl_excluded[offsetm]);
  uint32_t* genome_main_ptr = &(g_genome_main[offset]);
  uintptr_t* geno_ptr = (uintptr_t*)g_geno;
  uintptr_t* masks_ptr = g_masks;
  uintptr_t* mmasks_ptr = g_mmasks;
  uint32_t end_idx = g_thread_start[tidx + 1];
  uint32_t is_last_block;
  while (1) {
    is_last_block = g_is_last_thread_block;
    incr_dists_rm_inv(missing_ptr, mmasks_ptr, sample_ct_m1, ulii, end_idx);
    if (is_last_block || ((g_thread_spawn_ct % (GENOME_MULTIPLEX / BITCT)) == (GENOME_MULTIPLEX / BITCT) - 1)) {
      incr_genome(genome_main_ptr, geno_ptr, masks_ptr, sample_ct, ulii, end_idx);
    }
    if ((!tidx) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

void incr_dists(double* dists, uintptr_t* geno, uintptr_t* masks, double* weights, uint32_t start_idx, uint32_t end_idx) {
  uintptr_t* glptr;
  uintptr_t ulii;
  uintptr_t mask_fixed;
  uintptr_t uljj;
  uintptr_t* mptr;
  double* weights1 = &(weights[16384]);
#ifdef __LP64__
  double* weights2 = &(weights[32768]);
  double* weights3 = &(weights[36864]);
  double* weights4 = &(weights[40960]);
#endif
  uint32_t uii;
  uint32_t ujj;
  for (uii = start_idx; uii < end_idx; uii++) {
    glptr = geno;
    ulii = geno[uii];
    mptr = masks;
    mask_fixed = masks[uii];
#ifdef __LP64__
    if (mask_fixed == ~ZEROLU) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + weights[uljj & 16383];
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + weights[uljj & 16383];
	dists++;
      }
    }
#else
    if (mask_fixed == 0x0fffffff) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
	*dists += weights1[uljj >> 14] + weights[uljj & 16383];
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
	*dists += weights1[uljj >> 14] + weights[uljj & 16383];
	dists++;
      }
    }
#endif
  }
}

THREAD_RET_TYPE calc_wdist_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t ulii = g_thread_start[tidx];
  uintptr_t uljj = g_thread_start[0];
  uintptr_t offset = (((uint64_t)ulii) * (ulii - 1) - ((uint64_t)uljj) * (uljj - 1)) / 2;
  double* dists_ptr = &(g_dists[offset]);
  uintptr_t* geno_ptr = (uintptr_t*)g_geno;
  uintptr_t* masks_ptr = g_masks;
  uintptr_t* mmasks_ptr = g_mmasks;
  double* subset_weights_ptr = g_subset_weights;
  uint32_t* subset_weights_i_ptr = g_subset_weights_i;
  uint32_t* weighted_missing_ptr = &(g_missing_tot_weights[offset]);
  uint32_t end_idx = g_thread_start[tidx + 1];
  uint32_t is_last_block;
  while (1) {
    is_last_block = g_is_last_thread_block;
    incr_dists(dists_ptr, geno_ptr, masks_ptr, subset_weights_ptr, ulii, end_idx);
    if (is_last_block || (g_thread_spawn_ct & 1)) {
      // subset_weights_i is stationary here
      incr_wt_dist_missing(weighted_missing_ptr, subset_weights_i_ptr, mmasks_ptr, ulii, end_idx);
    }
    if ((!tidx) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

void incr_dists_r(double* dists, uintptr_t* geno, uintptr_t* masks, uint32_t tidx, double* weights) {
  uintptr_t* glptr;
  uintptr_t* maskptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t basemask;
  double* weights1 = &(weights[32768]);
#ifdef __LP64__
  double* weights2 = &(weights[65536]);
  double* weights3 = &(weights[98304]);
#endif
  uint32_t uii;
  uint32_t ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    glptr = geno;
    ulii = geno[uii];
    maskptr = masks;
    basemask = masks[uii];
    if (!basemask) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = ((*glptr++) + ulii) | (*maskptr++);
#ifdef __LP64__
	*dists += weights[(uint16_t)uljj] + weights1[(uint16_t)(uljj >> 16)] + weights2[(uint16_t)(uljj >> 32)] + weights3[uljj >> 48];
#else
	*dists += weights[(uint16_t)uljj] + weights1[uljj >> 16];
#endif
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
        uljj = ((*glptr++) + ulii) | ((*maskptr++) | basemask);
#ifdef __LP64__
	*dists += weights[(uint16_t)uljj] + weights1[(uint16_t)(uljj >> 16)] + weights2[(uint16_t)(uljj >> 32)] + weights3[uljj >> 48];
#else
	*dists += weights[(uint16_t)uljj] + weights1[uljj >> 16];
#endif
	dists++;
      }
    }
  }
}

THREAD_RET_TYPE calc_rel_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t ulii = g_thread_start[tidx];
  uintptr_t uljj = g_thread_start[0];
  uintptr_t offset = (((uint64_t)ulii) * (ulii - 1) - ((uint64_t)uljj) * (uljj - 1)) / 2;
  double* rel_ptr = &(g_rel_dists[offset]);
  uintptr_t* geno_ptr = (uintptr_t*)g_geno;
  uintptr_t* masks_ptr = g_masks;
  uintptr_t* mmasks_ptr = g_mmasks;
  uint32_t* missing_ptr = &(g_missing_dbl_excluded[offset]);
  double* subset_weights_ptr = g_subset_weights;
  uint32_t end_idx = g_thread_start[tidx + 1];
  uint32_t is_last_block;
  while (1) {
    is_last_block = g_is_last_thread_block;
    incr_dists_r(rel_ptr, geno_ptr, masks_ptr, (uint32_t)tidx, subset_weights_ptr);
    if (is_last_block || ((g_thread_spawn_ct % 3) == 2)) {
      incr_dists_rm(missing_ptr, mmasks_ptr, ulii, end_idx);
    }
    if ((!tidx) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE calc_wt_rel_thread(void* arg) {
  // this needs more work
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t ulii = g_thread_start[tidx];
  uintptr_t uljj = g_thread_start[0];
  uintptr_t offset = (((uint64_t)ulii) * (ulii - 1) - ((uint64_t)uljj) * (uljj - 1)) / 2;
  double* rel_ptr = &(g_rel_dists[offset]);
  uintptr_t* geno_ptr = (uintptr_t*)g_geno;
  uintptr_t* masks_ptr = g_masks;
  uintptr_t* mmasks_ptr = g_mmasks;
  uint32_t* missing_ptr = &(g_missing_dbl_excluded[offset]);
  double* subset_weights_ptr = g_subset_weights;
  uint32_t end_idx = g_thread_start[tidx + 1];
  uint32_t is_last_block;
  while (1) {
    is_last_block = g_is_last_thread_block;
    incr_dists_r(rel_ptr, geno_ptr, masks_ptr, (uint32_t)tidx, subset_weights_ptr);
    if (is_last_block || ((g_thread_spawn_ct % 3) == 2)) {
      incr_dists_rm(missing_ptr, mmasks_ptr, ulii, end_idx);
    }
    if ((!tidx) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

THREAD_RET_TYPE calc_missing_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t ulii = g_thread_start[tidx];
  uintptr_t uljj = g_thread_start[0];
  uintptr_t offset = (((uint64_t)ulii) * (ulii - 1) - ((uint64_t)uljj) * (uljj - 1)) / 2;
  uintptr_t* mmasks_ptr = g_mmasks;
  uint32_t* missing_ptr = &(g_missing_dbl_excluded[offset]);
  uint32_t end_idx = g_thread_start[tidx + 1];
  while (1) {
    incr_dists_rm(missing_ptr, mmasks_ptr, ulii, end_idx);
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

void groupdist_jack(uint32_t jackknife_d, uint32_t case_ct, uint32_t ctrl_ct, double reg_tot_x, double reg_tot_xy, double reg_tot_y, uint32_t* uibuf, double* jackknife_precomp, double* dists, uintptr_t* pheno_c, double* returns) {
  uint32_t* uiptr = uibuf;
  uint32_t* uiptr2 = &(uibuf[jackknife_d]);
  double neg_tot_uu = 0.0;
  double neg_tot_au = 0.0;
  double neg_tot_aa = 0.0;
  uint32_t neg_a = 0;
  uint32_t neg_u = 0;
  double* dptr;
  uint32_t* uiptr3;
  uintptr_t sample_idx;
  uint32_t uii;
  while (uiptr < uiptr2) {
    dptr = &(jackknife_precomp[(*uiptr++) * JACKKNIFE_VALS_GROUPDIST]);
    neg_tot_uu += *dptr++;
    neg_tot_au += *dptr++;
    neg_tot_aa += *dptr++;
  }
  uiptr = uibuf;
  while (uiptr < uiptr2) {
    sample_idx = *uiptr;
    uiptr3 = uibuf;
    dptr = &(dists[(sample_idx * (sample_idx - 1)) / 2]);
    if (IS_SET(pheno_c, sample_idx)) {
      neg_a++;
      while (uiptr3 < uiptr) {
	uii = *uiptr3++;
	if (IS_SET(pheno_c, uii)) {
	  neg_tot_aa -= dptr[uii];
	} else {
	  neg_tot_au -= dptr[uii];
	}
      }
    } else {
      neg_u++;
      while (uiptr3 < uiptr) {
	uii = *uiptr3++;
	if (IS_SET(pheno_c, uii)) {
	  neg_tot_au -= dptr[uii];
	} else {
	  neg_tot_uu -= dptr[uii];
	}
      }
    }
    uiptr++;
  }
  returns[0] = (reg_tot_x - neg_tot_aa) / (double)(((intptr_t)(case_ct - neg_a) * (case_ct - neg_a - 1)) / 2);
  returns[1] = (reg_tot_xy - neg_tot_au) / (double)((intptr_t)(case_ct - neg_a) * (ctrl_ct - neg_u));
  returns[2] = (reg_tot_y - neg_tot_uu) / (double)(((intptr_t)(ctrl_ct - neg_u) * (ctrl_ct - neg_u - 1)) / 2);
}

void small_remap(uint32_t* uibuf, uint32_t ct, uint32_t dd) {
  uint32_t* uibuf_end = &(uibuf[dd]);
  uint32_t missings = 0;
  uint32_t curpos = 0;
  do {
    if (!IS_SET(g_pheno_nm, curpos)) {
      missings++;
    } else if (*uibuf == curpos - missings) {
      *uibuf++ = curpos;
    }
    curpos++;
  } while (uibuf < uibuf_end);
}

void pick_d(unsigned char* cbuf, uint32_t ct, uint32_t dd, sfmt_t* sfmtp) {
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  memset(cbuf, 0, ct);
#ifdef __LP64__
  ukk = (uint32_t)(0x100000000LLU % ct);
#else
  ukk = 2 * (0x80000000U % ct);
#endif
  for (uii = 0; uii < dd; uii++) {
    do {
      do {
        ujj = sfmt_genrand_uint32(sfmtp);
      } while (ujj < ukk);
      ujj %= ct;
    } while (cbuf[ujj]);
    cbuf[ujj] = 1;
  }
}

void pick_d_small(unsigned char* tmp_cbuf, uint32_t* uibuf, uint32_t ct, uint32_t dd, sfmt_t* sfmtp) {
  uint32_t uii;
  pick_d(tmp_cbuf, ct, dd, sfmtp);
  for (uii = 0; uii < ct; uii++) {
    if (tmp_cbuf[uii]) {
      *uibuf++ = uii;
    }
  }
  *uibuf = ct;
}

THREAD_RET_TYPE groupdist_jack_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t sample_ct = g_sample_ct;
  uint32_t case_ct = g_case_ct;
  uint32_t ctrl_ct = g_ctrl_ct;
  uint32_t jackknife_d = g_jackknife_d;
  uint32_t* uibuf = (uint32_t*)(&(g_geno[tidx * round_up_pow2(case_ct + ctrl_ct + (jackknife_d + 1) * sizeof(int32_t), CACHELINE)]));
  unsigned char* cbuf = &(g_geno[tidx * round_up_pow2(case_ct + ctrl_ct + (jackknife_d + 1) * sizeof(int32_t), CACHELINE) + (jackknife_d + 1) * sizeof(int32_t)]);
  uintptr_t jackknife_iters = g_jackknife_iters;
  uintptr_t uljj = jackknife_iters / 100;
  double reg_tot_x = g_reg_tot_x;
  double reg_tot_xy = g_reg_tot_xy;
  double reg_tot_y = g_reg_tot_y;
  double* jackknife_precomp = g_jackknife_precomp;
  double* dists = g_dists;
  uintptr_t* pheno_c = g_pheno_c;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  double returns[3];
  double results[9];
  uintptr_t ulii;
  fill_double_zero(9, results);
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, uibuf, case_ct + ctrl_ct, jackknife_d, sfmtp);
    if (case_ct + ctrl_ct < sample_ct) {
      small_remap(uibuf, case_ct + ctrl_ct, jackknife_d);
    }
    groupdist_jack(jackknife_d, case_ct, ctrl_ct, reg_tot_x, reg_tot_xy, reg_tot_y, uibuf, jackknife_precomp, dists, pheno_c, returns);
    results[0] += returns[0];
    results[1] += returns[1];
    results[2] += returns[2];
    results[3] += returns[0] * returns[0];
    results[4] += returns[1] * returns[1];
    results[5] += returns[2] * returns[2];
    results[6] += returns[0] * returns[1];
    results[7] += returns[0] * returns[2];
    results[8] += returns[1] * returns[2];
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100LLU) / jackknife_iters;
      printf("\r%" PRIuPTR "%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1LLU) * jackknife_iters) / 100;
    }
  }
  // don't write until end, to avoid false sharing
  for (ulii = 0; ulii < 9; ulii++) {
    g_calc_result[tidx][ulii] = results[ulii];
  }
  THREAD_RETURN;
}

double regress_rel_jack(uint32_t jackknife_d, uint32_t sample_ct, double reg_tot_xy, double reg_tot_x, double reg_tot_y, double reg_tot_xx, double reg_tot_yy, uint32_t* uibuf, double* jackknife_precomp, double* rel_dists, double* pheno_packed, double* ret2_ptr) {
  uint32_t* uiptr = uibuf;
  uint32_t* uiptr2 = &(uibuf[jackknife_d]);
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double* dptr;
  double* dptr2;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double neg_tot_yy = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  while (uiptr < uiptr2) {
    dptr2 = &(jackknife_precomp[(*uiptr++) * JACKKNIFE_VALS_REL]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  uiptr = uibuf;
  for (uii = 1; uii < jackknife_d; uii++) {
    ujj = *(++uiptr);
    dxx1 = pheno_packed[ujj];
    uiptr2 = uibuf;
    dptr = &(rel_dists[((uintptr_t)ujj * (ujj - 1)) / 2]);
    while (uiptr2 < uiptr) {
      ukk = *uiptr2++;
      dxx = (dxx1 + pheno_packed[ukk]) * 0.5;
      dyy = dptr[ukk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = reg_tot_y - neg_tot_y;
  dyy = ((int32_t)sample_ct) - ((int32_t)jackknife_d);
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_x - neg_tot_x) / dyy) / ((reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = reg_tot_x - neg_tot_x;
  return ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_y - neg_tot_y) / dyy) / ((reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

THREAD_RET_TYPE regress_rel_jack_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t jackknife_iters = g_jackknife_iters;
  double reg_tot_xy = g_reg_tot_xy;
  double reg_tot_x = g_reg_tot_x;
  double reg_tot_y = g_reg_tot_y;
  double reg_tot_xx = g_reg_tot_xx;
  double reg_tot_yy = g_reg_tot_yy;
  uint32_t jackknife_d = g_jackknife_d;
  uint32_t* uibuf = (uint32_t*)(&(g_geno[tidx * round_up_pow2(sample_ct + (jackknife_d + 1) * sizeof(int32_t), CACHELINE)]));
  unsigned char* cbuf = &(g_geno[tidx * round_up_pow2(sample_ct + (jackknife_d + 1) * sizeof(int32_t), CACHELINE) + (jackknife_d + 1) * sizeof(int32_t)]);
  double* jackknife_precomp = g_jackknife_precomp;
  double* rel_dists = g_rel_dists;
  double* pheno_packed = g_pheno_packed;
  uintptr_t uljj = jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint64_t ulii;
  double dxx;
  double ret2;
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, uibuf, sample_ct, jackknife_d, sfmtp);
    dxx = regress_rel_jack(jackknife_d, sample_ct, reg_tot_xy, reg_tot_x, reg_tot_y, reg_tot_xx, reg_tot_yy, uibuf, jackknife_precomp, rel_dists, pheno_packed, &ret2);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100LLU) / jackknife_iters;
      printf("\r%" PRIuPTR "%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1LLU) * jackknife_iters) / 100;
    }
  }
  g_calc_result[tidx][0] = sum;
  g_calc_result[tidx][1] = sum_sq;
  g_calc_result[tidx][2] = sum2;
  g_calc_result[tidx][3] = sum2_sq;
  THREAD_RETURN;
}

void print_pheno_stdev(double* pheno_d, uint32_t sample_ct) {
  double reg_tot_x = 0.0;
  double reg_tot_xx = 0.0;
  double dxx;
  uint32_t uii;
  for (uii = 0; uii < sample_ct; uii++) {
    dxx = pheno_d[uii];
    reg_tot_x += dxx;
    reg_tot_xx += dxx * dxx;
  }
  LOGPRINTF("Phenotype stdev: %g\n", sqrt((reg_tot_xx - reg_tot_x * reg_tot_x / sample_ct) / (sample_ct - 1)));
}

uint32_t set_default_jackknife_d(uint32_t ct) {
  uint32_t dd = (uint32_t)pow((double)ct, 0.600000000001);
  LOGPRINTF("Setting d=%u for jackknife.\n", dd);
  return dd;
}

int32_t regress_rel_main(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Rel_info* relip, pthread_t* threads, double* pheno_d) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t regress_rel_iters = relip->regress_rel_iters;
  double reg_tot_xy = 0;
  double reg_tot_x = 0;
  double reg_tot_y = 0;
  double reg_tot_xx = 0;
  double reg_tot_yy = 0;
  uint32_t regress_rel_d = relip->regress_rel_d;
  double* rel_ptr;
  double* pheno_ptr;
  double* pheno_ptr2;
  double* jp_fixed_ptr;
  double* jp_moving_ptr;
  uint32_t uii;
  uintptr_t ulii;
  uintptr_t trimatrix_size;
  double trimatrix_size_recip;
  double half_avg_pheno;
  double dxx;
  double dyy;
  double dxxyy;
  double dxxsq;
  double dyysq;
  g_sample_ct = sample_ct;
  if (bigstack_alloc_d(sample_ct, &g_pheno_packed)) {
    return RET_NOMEM;
  }
  collapse_copy_phenod(g_pheno_packed, pheno_d, sample_exclude, unfiltered_sample_ct, sample_ct);
  print_pheno_stdev(g_pheno_packed, sample_ct);
  trimatrix_size = ((uintptr_t)sample_ct * (sample_ct - 1)) / 2;
  rel_ptr = g_rel_dists;
  pheno_ptr = g_pheno_packed;
  if (bigstack_calloc_d(sample_ct * JACKKNIFE_VALS_REL, &g_jackknife_precomp)) {
    return RET_NOMEM;
  }
  for (uii = 1; uii < sample_ct; uii++) {
    half_avg_pheno = *(++pheno_ptr);
    pheno_ptr2 = g_pheno_packed;
    jp_fixed_ptr = &(g_jackknife_precomp[uii * JACKKNIFE_VALS_REL]);
    jp_moving_ptr = g_jackknife_precomp;
    while (pheno_ptr2 < pheno_ptr) {
      dxx = (half_avg_pheno + (*pheno_ptr2++)) * 0.5;
      dyy = (*rel_ptr++);
      dxxyy = dxx * dyy;
      dxxsq = dxx * dxx;
      dyysq = dyy * dyy;
      reg_tot_xy += dxxyy;
      jp_fixed_ptr[0] += dxxyy;
      *jp_moving_ptr += dxxyy;
      jp_moving_ptr++;
      reg_tot_x += dxx;
      jp_fixed_ptr[1] += dxx;
      *jp_moving_ptr += dxx;
      jp_moving_ptr++;
      reg_tot_y += dyy;
      jp_fixed_ptr[2] += dyy;
      *jp_moving_ptr += dyy;
      jp_moving_ptr++;
      reg_tot_xx += dxxsq;
      jp_fixed_ptr[3] += dxxsq;
      *jp_moving_ptr += dxxsq;
      jp_moving_ptr++;
      reg_tot_yy += dyysq;
      jp_fixed_ptr[4] += dyysq;
      *jp_moving_ptr += dyysq;
      jp_moving_ptr++;
    }
  }
  g_reg_tot_xy = reg_tot_xy;
  g_reg_tot_x = reg_tot_x;
  g_reg_tot_y = reg_tot_y;
  g_reg_tot_xx = reg_tot_xx;
  g_reg_tot_yy = reg_tot_yy;
  trimatrix_size_recip = 1.0 / (double)trimatrix_size;
  LOGPRINTF("Regression slope (y = genomic relationship, x = avg phenotype): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y * trimatrix_size_recip) / (reg_tot_xx - reg_tot_x * reg_tot_x * trimatrix_size_recip));
  LOGPRINTF("                 (y = avg phenotype, x = genomic relationship): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y * trimatrix_size_recip) / (reg_tot_yy - reg_tot_y * reg_tot_y * trimatrix_size_recip));
  g_jackknife_iters = (regress_rel_iters + g_thread_ct - 1) / g_thread_ct;
  if (regress_rel_d) {
    g_jackknife_d = regress_rel_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(sample_ct);
  }
  if (bigstack_alloc_uc(g_thread_ct * round_up_pow2(sample_ct * (g_jackknife_d + 1) * sizeof(int32_t), CACHELINE), &g_geno)) {
    return RET_NOMEM;
  }
  if (bigstack_init_sfmtp(g_thread_ct)) {
    return RET_NOMEM;
  }
  if (spawn_threads(threads, &regress_rel_jack_thread, g_thread_ct)) {
    return RET_THREAD_CREATE_FAIL;
  }
  ulii = 0;
  regress_rel_jack_thread((void*)ulii);
  dxx = g_calc_result[0][0]; // relationship on pheno
  dxxsq = g_calc_result[0][1];

  dyy = g_calc_result[0][2]; // pheno on relationship
  dyysq = g_calc_result[0][3];

  join_threads(threads, g_thread_ct);
  for (uii = 0; uii < g_thread_ct - 1; uii++) {
    dxx += g_calc_result[uii + 1][0];
    dxxsq += g_calc_result[uii + 1][1];
    dyy += g_calc_result[uii + 1][2];
    dyysq += g_calc_result[uii + 1][3];
  }
  ulii = g_jackknife_iters * g_thread_ct;
  putc_unlocked('\r', stdout);
  LOGPRINTF("Jackknife s.e. (y = genomic relationship): %g\n", sqrt(((sample_ct - g_jackknife_d) / (double)g_jackknife_d) * (dxxsq - dxx * dxx / (double)ulii) / ((double)ulii - 1)));
  LOGPRINTF("               (y = phenotype): %g\n", sqrt(((sample_ct - g_jackknife_d) / (double)g_jackknife_d) * (dyysq - dyy * dyy / (double)ulii) / ((double)ulii - 1)));
  bigstack_reset(bigstack_mark);
  return 0;
}

// Replaces matrix[][] with mult_val * matrix[][] + add_val * I.
// Multithreading doesn't help here.
void matrix_const_mult_add(uint32_t sample_ct, double* matrix, double mult_val, double add_val) {
  uint32_t uii;
  uint32_t loop_end = sample_ct - 1;
  uint32_t ujj;
  double* dptr = matrix;
#ifdef __LP64__
  __m128d* vptr;
  __m128d v_mult_val = _mm_set1_pd(mult_val);
#endif
  for (uii = 0; uii < loop_end; uii++) {
    *dptr = (*dptr) * mult_val + add_val;
    dptr++;
#ifdef __LP64__
    if ((uintptr_t)dptr & 8) {
      *dptr *= mult_val;
      dptr++;
      ujj = 1;
    } else {
      ujj = 0;
    }
    vptr = (__m128d*)dptr;
    while (ujj < loop_end) {
      *vptr = _mm_mul_pd(*vptr, v_mult_val);
      vptr++;
      ujj += 2;
    }
    dptr = (double*)vptr;
    if (ujj < sample_ct) {
      *dptr *= mult_val;
      dptr++;
    }
#else
    for (ujj = 0; ujj < sample_ct; ujj++) {
      *dptr *= mult_val;
      dptr++;
    }
#endif
  }
  *dptr = (*dptr) * mult_val + add_val;
}

// sums[idx] = matrix[idx][1] + matrix[idx][2] + ....  Ideally, we can assume
// the matrix is symmetric and only reference the FORTRAN upper right, but for
// now we can't.
void matrix_row_sum_ur(uintptr_t sample_ct, double* sums, double* matrix) {
  uintptr_t sample_idx;
  double* dptr;
  double acc;
  double* sptr_end;
  double* sptr;
  fill_double_zero(sample_ct, sums);
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    dptr = &(matrix[sample_idx * sample_ct]);
    acc = 0.0;
    sptr_end = &(sums[sample_idx]);
    sptr = sums;
    while (sptr < sptr_end) {
      acc += *dptr;
      *sptr += *dptr++;
      sptr++;
    }
    *sptr += acc + *dptr;
  }
}

#ifndef NOLAPACK
static inline void collapse_copy_phenod_incl(double* target, double* pheno_d, uintptr_t* sample_include, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct) {
  uintptr_t sample_uidx = 0;
  double* target_end = &(target[sample_ct]);
  uintptr_t delta;
  do {
    sample_uidx = next_set_ul_unsafe(sample_include, sample_uidx);
    delta = next_unset_ul(sample_include, sample_uidx, unfiltered_sample_ct) - sample_uidx;
    memcpy(target, &(pheno_d[sample_uidx]), delta * sizeof(double));
    target = &(target[delta]);
    sample_uidx += delta;
  } while (target < target_end);
}

// one-trait REML via EM.
//
// wkbase is assumed to have space for three cache-aligned
// sample_ct * sample_ct double matrices plus three more rows.  The unpacked
// relationship matrix is stored in the SECOND slot.
//
// g_sample_ct currently must be set.
void reml_em_one_trait(double* wkbase, double* pheno, double* covg_ref, double* covr_ref, double tol, uint32_t is_strict) {
  double ll_change;
  uintptr_t sample_ct = g_sample_ct;
  int64_t mat_offset = sample_ct;
  double* rel_dists;
  MATRIX_INVERT_BUF1_TYPE* irow;
  __CLPK_integer lwork;
  double* row;
  double* row2;
  double* work;
  double* dptr;
  double* dptr2;
  double* matrix_pvg;
  __CLPK_integer sample_ct_li = sample_ct;
  double dxx;
  double dyy;
  double max_jump;
  double dzz;
  double dlg;
  double dle;
  double covg_cur_change = 1.0;
  double covr_cur_change = 1.0;
  double covg_last_change;
  double covr_last_change;
  double sample_ct_recip = 1 / ((double)((intptr_t)sample_ct));
  uintptr_t sample_idx;
  int32_t jj;
#ifdef _WIN32
  char blas_char;
  int32_t sample_ct_i32 = sample_ct;
#endif
  mat_offset = round_up_pow2_ull(mat_offset * mat_offset, CACHELINE);
  rel_dists = &(wkbase[mat_offset]);
  row = &(wkbase[mat_offset * 3]);
  irow = (MATRIX_INVERT_BUF1_TYPE*)row;
  row2 = &(row[sample_ct]);
  work = &(wkbase[mat_offset * 2]);
  lwork = mat_offset;
  matrix_pvg = work;
  if (!lwork) {
    lwork = CACHELINE_DBL;
  }
  fill_double_zero(mat_offset, matrix_pvg);
  fill_double_zero(sample_ct, row2);
  printf("      ");
  do {
    memcpy(wkbase, rel_dists, mat_offset * sizeof(double));
    matrix_const_mult_add(sample_ct, wkbase, *covg_ref, *covr_ref);
    invert_matrix(sample_ct_li, wkbase, irow, work);
    matrix_row_sum_ur(sample_ct, row, wkbase);
    dxx = 0.0;
    dptr = row;
    dptr2 = &(row[sample_ct]);
    while (dptr < dptr2) {
      dxx += *dptr++;
    }
    dxx = -1 / dxx;
#ifdef _WIN32
    jj = 1;
    dger_(&sample_ct_i32, &sample_ct_i32, &dxx, row, &jj, row, &jj, wkbase, &sample_ct_i32);
    // todo: test dsymm
    dyy = 1.0;
    dzz = 0.0;
    col_major_matrix_multiply(sample_ct_i32, sample_ct_i32, sample_ct_i32, wkbase, rel_dists, matrix_pvg);
    // blas_char = 'N';
    // dgemm_(&blas_char, &blas_char, &sample_ct_i32, &sample_ct_i32, &sample_ct_i32, &dyy, wkbase, &sample_ct_i32, rel_dists, &sample_ct_i32, &dzz, matrix_pvg, &sample_ct_i32);
    dlg = 0.0;
    dle = 0.0;
    jj = sample_ct + 1;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      dlg -= matrix_pvg[sample_idx * jj];
      dle -= wkbase[sample_idx * jj];
    }
    blas_char = 'U';
    jj = 1;
    dsymv_(&blas_char, &sample_ct_i32, &dyy, wkbase, &sample_ct_i32, pheno, &jj, &dzz, row2, &jj);
    dsymv_(&blas_char, &sample_ct_i32, &dyy, matrix_pvg, &sample_ct_i32, row2, &jj, &dzz, row, &jj);
    dlg += ddot_(&sample_ct_i32, pheno, &jj, row, &jj);
    dsymv_(&blas_char, &sample_ct_i32, &dyy, wkbase, &sample_ct_i32, row2, &jj, &dzz, row, &jj);
    dle += ddot_(&sample_ct_i32, pheno, &jj, row, &jj);
#else
    cblas_dger(CblasColMajor, sample_ct, sample_ct, dxx, row, 1, row, 1, wkbase, sample_ct);
    // unfortunately, cblas_dsymm is much worse than cblas_dgemm on OS X
    col_major_matrix_multiply(sample_ct, sample_ct, sample_ct, wkbase, rel_dists, matrix_pvg);
    // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, sample_ct, sample_ct, sample_ct, 1.0, wkbase, sample_ct, rel_dists, sample_ct, 0.0, matrix_pvg, sample_ct);
    dlg = 0.0;
    dle = 0.0;
    jj = sample_ct + 1;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      dlg -= matrix_pvg[sample_idx * jj];
      dle -= wkbase[sample_idx * jj];
    }
    cblas_dsymv(CblasColMajor, CblasUpper, sample_ct, 1.0, wkbase, sample_ct, pheno, 1, 0.0, row2, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, sample_ct, 1.0, matrix_pvg, sample_ct, row2, 1, 0.0, row, 1);
    dlg += cblas_ddot(sample_ct, pheno, 1, row, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, sample_ct, 1.0, wkbase, sample_ct, row2, 1, 0.0, row, 1);
    dle += cblas_ddot(sample_ct, pheno, 1, row, 1);
#endif
    covg_last_change = covg_cur_change;
    covr_last_change = covr_cur_change;
    covg_cur_change = (*covg_ref) * (*covg_ref) * dlg * sample_ct_recip;
    covr_cur_change = (*covr_ref) * (*covr_ref) * dle * sample_ct_recip;
    if (is_strict) {
      max_jump = 1.0;
    } else {
      // acceleration factor:
      // min(half covg distance to 0 or 1, covr distance to 0 or 1, pi/4 divided
      // by last angular change, 1.0 / (1 - ratio of last two step lengths),
      // MAX_EM_ACCEL)
      dxx = atan2(covg_last_change, covr_last_change) - atan2(covg_cur_change, covr_cur_change);
      if (dxx < 0.0) {
	dxx = -dxx;
      }
      if (dxx > PI) {
	dxx = 2 * PI - dxx;
      }
      dyy = sqrt((covg_cur_change * covg_cur_change + covr_cur_change * covr_cur_change) / (covg_last_change * covg_last_change + covr_last_change * covr_last_change));
      if (covg_cur_change < 0.0) {
	max_jump = *covg_ref * (-0.5) / covg_cur_change;
      } else {
	max_jump = (1.0 - *covg_ref) * 0.5 / covg_cur_change;
      }
      if (covr_cur_change < 0.0) {
	dzz = *covr_ref * (-0.5) / covr_cur_change;
      } else {
	dzz = (1.0 - *covr_ref) * 0.5 / covr_cur_change;
      }
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      dzz = (PI / 4) / dxx;
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      if (dyy < 1.0) {
	dzz = 1 / (1.0 - dyy);
      }
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      if (max_jump < 1.0) {
	max_jump = 1.0;
      } else if (max_jump > MAX_EM_ACCEL) {
	max_jump = MAX_EM_ACCEL;
      }
    }
    *covg_ref += covg_cur_change * max_jump;
    *covr_ref += covr_cur_change * max_jump;
    ll_change = (covg_cur_change * dlg) + (covr_cur_change * dle);
    printf("\b\b\b\b\b\b      \rcovg: %g  covr: %g  EM step log likelihood change: %g", *covg_ref, *covr_ref, ll_change);
    fflush(stdout);
  } while (ll_change > tol);
  putc_unlocked('\n', stdout);
  sprintf(g_logbuf, "covg: %g  covr: %g\n", *covg_ref, *covr_ref);
  logstr(g_logbuf);
}

void mean_zero_var_one_in_place(uint32_t sample_ct, double* pheno_d) {
  double sum = 0.0;
  double ssq = 0.0;
  double* dptr = pheno_d;
  double dxx;
  double sample_ctd;
  double mean;
  double stdev_recip;
  uint32_t sample_idx;
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    dxx = *dptr++;
    sum += dxx;
    ssq += dxx * dxx;
  }
  sample_ctd = (double)((int32_t)sample_ct);
  mean = sum / sample_ctd;
  // --unrelated-heritability is intended to be a straight port of Carson's
  // MATLAB script; it does not have "- 1" in the numerator here.
  stdev_recip = sqrt(sample_ctd / (ssq - sum * mean));
  dptr = pheno_d;
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    *dptr = ((*dptr) - mean) * stdev_recip;
    dptr++;
  }
}

int32_t calc_unrelated_herit(uint64_t calculation_type, Rel_info* relip, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, double* pheno_d, double* rel_ibc) {
  uint32_t is_strict = (relip->modifier / REL_UNRELATED_HERITABILITY_STRICT) & 1;
  int32_t ibc_type = relip->ibc_type;
  double unrelated_herit_covg = relip->unrelated_herit_covg;
  double unrelated_herit_covr = relip->unrelated_herit_covr;
  double unrelated_herit_tol = relip->unrelated_herit_tol;
  uintptr_t ulii;
  uintptr_t uljj;
  double* pheno_ptr;
  double* ibc_ptr;
  double* rel_base;
  g_sample_ct = sample_ct;
  g_missing_dbl_excluded = nullptr;
  ulii = sample_ct;
  ulii = round_up_pow2(ulii * ulii, CACHELINE_DBL);
  rel_base = &(g_rel_dists[ulii]);
  ulii = ulii * 3 + round_up_pow2(sample_ct, CACHELINE_DBL) * 3;
  // no bigstack_shrink_top here since this actually grows the allocation...
  bigstack_reset(g_rel_dists);
  if (bigstack_alloc_d(ulii, &g_rel_dists)) {
    return RET_NOMEM;
  }
  pheno_ptr = &(g_rel_dists[ulii - round_up_pow2(sample_ct, CACHELINE_DBL)]);
  collapse_copy_phenod(pheno_ptr, pheno_d, sample_exclude, unfiltered_sample_ct, sample_ct);
  mean_zero_var_one_in_place(sample_ct, pheno_ptr);
  if (calculation_type & CALC_IBC) {
    ibc_ptr = &(rel_ibc[ibc_type * sample_ct]);
  } else {
    ibc_ptr = rel_ibc;
  }
  for (ulii = 0; ulii < sample_ct; ulii++) {
    memcpy(&(rel_base[ulii * sample_ct]), &(g_rel_dists[(ulii * (ulii - 1)) / 2]), ulii * sizeof(double));
    rel_base[ulii * (sample_ct + 1)] = *ibc_ptr++;
    for (uljj = ulii + 1; uljj < sample_ct; uljj++) {
      rel_base[ulii * sample_ct + uljj] = g_rel_dists[(uljj * (uljj - 1)) / 2 + ulii];
    }
  }
  reml_em_one_trait(g_rel_dists, pheno_ptr, &unrelated_herit_covg, &unrelated_herit_covr, unrelated_herit_tol, is_strict);
  LOGPRINTF("h^2 estimate: %g\n", unrelated_herit_covg);
  return 0;
}

int32_t unrelated_herit_batch(uint32_t load_grm_bin, char* grmname, char* phenoname, uint32_t mpheno_col, char* phenoname_str, int32_t missing_pheno, Rel_info* relip) {
  char* grmname_end = (char*)memchr(grmname, 0, FNAMESIZE);
  FILE* infile = nullptr;
  FILE* grm_binfile = nullptr;
  gzFile grm_gzfile = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t max_sample_id_len = 4;
  uintptr_t unfiltered_sample_ct = 0;
  uintptr_t sample_uidx = 0;
  uintptr_t line_idx = 0;
  double unrelated_herit_tol = relip->unrelated_herit_tol;
  double unrelated_herit_covg = relip->unrelated_herit_covg;
  double unrelated_herit_covr = relip->unrelated_herit_covr;
  uint32_t is_strict = (relip->modifier / REL_UNRELATED_HERITABILITY_STRICT) & 1;
  uintptr_t* pheno_c = nullptr;
  double* pheno_d = nullptr;
  uintptr_t cur_sample_id_len;
  uintptr_t unfiltered_sample_ctl;
  uintptr_t* pheno_nm;
  uintptr_t pheno_nm_ct;
  uintptr_t sample_uidx2;
  char* sorted_ids;
  uint32_t* id_map;
  char* bufptr;
  char* bufptr2;
  double* matrix_wkbase;
  double* pheno_ptr;
  double* rel_base;
  double* row_ptr;
  uint64_t fpos;
  uintptr_t ulii;
  uintptr_t uljj;
  double dxx;
  float fxx;
  int32_t retval;
  // 1. load IDs
  // 2. load phenotypes and check for missing samples
  // 3. collapse phenotypes if necessary,  load (subset of) relationship matrix
  // 4. call reml_em_one_trait()
  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(grmname, "r", &infile)) {
    goto unrelated_herit_batch_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, grmname);
      goto unrelated_herit_batch_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = token_endnn(bufptr);
    cur_sample_id_len = bufptr2 - bufptr;
    bufptr2 = skip_initial_spaces(bufptr2);
    if (is_eoln_kns(*bufptr2)) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, grmname);
      goto unrelated_herit_batch_ret_INVALID_FORMAT_2;
    }
    cur_sample_id_len += strlen_se(bufptr2) + 2;
    if (cur_sample_id_len > max_sample_id_len) {
      max_sample_id_len = cur_sample_id_len;
    }
    unfiltered_sample_ct++;
  }
  if (!feof(infile)) {
    goto unrelated_herit_batch_ret_READ_FAIL;
  }
  if (unfiltered_sample_ct < 2) {
    logerrprint("Error: Less than two samples in .grm.id file.\n");
    goto unrelated_herit_batch_ret_INVALID_FORMAT;
  }
  rewind(infile);
  unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  if (bigstack_calloc_ul(unfiltered_sample_ctl, &pheno_nm) ||
      bigstack_end_alloc_c(unfiltered_sample_ct * max_sample_id_len, &sorted_ids) ||
      bigstack_end_alloc_ui(unfiltered_sample_ct, &id_map)) {
    goto unrelated_herit_batch_ret_NOMEM;
  }
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = token_endnn(bufptr);
    cur_sample_id_len = bufptr2 - bufptr;
    memcpy(&(sorted_ids[sample_uidx * max_sample_id_len]), bufptr, cur_sample_id_len);
    sorted_ids[sample_uidx * max_sample_id_len + (cur_sample_id_len++)] = '\t';
    bufptr2 = skip_initial_spaces(bufptr2);
    ulii = strlen_se(bufptr2);
    memcpy(&(sorted_ids[sample_uidx * max_sample_id_len + cur_sample_id_len]), bufptr2, ulii);
    sorted_ids[(sample_uidx++) * max_sample_id_len + cur_sample_id_len + ulii] = '\0';
  }
  if (!feof(infile)) {
    goto unrelated_herit_batch_ret_READ_FAIL;
  }
  for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
    id_map[sample_uidx] = sample_uidx;
  }
  if (qsort_ext(sorted_ids, unfiltered_sample_ct, max_sample_id_len, strcmp_deref, (char*)id_map, sizeof(int32_t))) {
    goto unrelated_herit_batch_ret_NOMEM;
  }

  fclose_null(&infile);
  if (fopen_checked(phenoname, "r", &infile)) {
    goto unrelated_herit_batch_ret_OPEN_FAIL;
  }
  retval = load_pheno(infile, unfiltered_sample_ct, 0, sorted_ids, max_sample_id_len, id_map, missing_pheno, 0, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d, nullptr, 0);
  bigstack_end_reset(bigstack_end_mark);
  fclose_null(&infile);
  if (retval) {
    goto unrelated_herit_batch_ret_1;
  }
  if (!pheno_d) {
    logerrprint("Error: --unrelated-heritability requires scalar phenotype.\n");
    goto unrelated_herit_batch_ret_INVALID_CMDLINE;
  }
  pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_sample_ctl);
  if (pheno_nm_ct < 2) {
    logerrprint("Error: Less than two phenotypes present.\n");
    goto unrelated_herit_batch_ret_INVALID_FORMAT;
  }
  ulii = round_up_pow2(pheno_nm_ct * pheno_nm_ct, CACHELINE_DBL);
  uljj = ulii * 3 + round_up_pow2(pheno_nm_ct, CACHELINE_DBL) * 3;
  if (bigstack_alloc_d(ulii, &matrix_wkbase)) {
    goto unrelated_herit_batch_ret_NOMEM;
  }
  g_sample_ct = pheno_nm_ct;
  pheno_ptr = &(matrix_wkbase[uljj - round_up_pow2(pheno_nm_ct, CACHELINE_DBL)]);
  collapse_copy_phenod_incl(pheno_ptr, pheno_d, pheno_nm, unfiltered_sample_ct, pheno_nm_ct);
  rel_base = &(matrix_wkbase[ulii]);
  mean_zero_var_one_in_place(pheno_nm_ct, pheno_ptr);
  sample_uidx = 0;
  if (load_grm_bin) {
    memcpy(grmname_end, ".grm.bin", 9);
    if (fopen_checked(grmname, FOPEN_RB, &grm_binfile)) {
      goto unrelated_herit_batch_ret_OPEN_FAIL;
    }
    if (fseeko(grm_binfile, 0, SEEK_END)) {
      goto unrelated_herit_batch_ret_READ_FAIL;
    }
    // n(n+1)/2 * 4 bytes per float
    fpos = ((uint64_t)unfiltered_sample_ct) * (unfiltered_sample_ct + 1) * 2;
    if (((int64_t)fpos) != ftello(grm_binfile)) {
      LOGPREPRINTFWW("Error: --grm-bin expects size of %s to be %" PRIu64 " bytes.\n", grmname, fpos);
      goto unrelated_herit_batch_ret_INVALID_FORMAT_2;
    }
    rewind(grm_binfile);
    for (ulii = 0; ulii < pheno_nm_ct; ulii++) {
      if (IS_SET(pheno_nm, sample_uidx)) {
	fpos = sample_uidx * (sample_uidx + 1) * (sizeof(float) / 2);
      } else {
        sample_uidx = next_set_ul_unsafe(pheno_nm, sample_uidx);
	fpos = sample_uidx * (sample_uidx + 1) * (sizeof(float) / 2);
        if (fseeko(grm_binfile, fpos, SEEK_SET)) {
          goto unrelated_herit_batch_ret_READ_FAIL;
	}
      }
      row_ptr = &(rel_base[ulii * pheno_nm_ct]);
      sample_uidx2 = 0;
      for (uljj = 0; uljj <= ulii; uljj++) {
        if (!IS_SET(pheno_nm, sample_uidx2)) {
          sample_uidx2 = next_set_ul_unsafe(pheno_nm, sample_uidx2);
          if (fseeko(grm_binfile, fpos + (sample_uidx2 * sizeof(float)), SEEK_SET)) {
            goto unrelated_herit_batch_ret_READ_FAIL;
	  }
	}
        if (fread(&fxx, sizeof(float), 1, grm_binfile) != sizeof(float)) {
	  goto unrelated_herit_batch_ret_READ_FAIL;
	}
        *row_ptr++ = (double)fxx;
        sample_uidx2++;
      }
      sample_uidx++;
    }
    fclose_null(&grm_binfile);
  } else {
    memcpy(grmname_end, ".grm.gz", 8);
    retval = gzopen_read_checked(grmname, &grm_gzfile);
    if (retval) {
      goto unrelated_herit_batch_ret_1;
    }
    ulii = 0;
    for (sample_uidx = 0; sample_uidx < pheno_nm_ct; sample_uidx++) {
      if (!IS_SET(pheno_nm, sample_uidx)) {
        for (sample_uidx2 = 0; sample_uidx2 <= sample_uidx; sample_uidx2++) {
          if (!gzgets(grm_gzfile, g_textbuf, MAXLINELEN)) {
	    goto unrelated_herit_batch_ret_READ_FAIL;
	  }
	  if (!g_textbuf[MAXLINELEN - 1]) {
	    goto unrelated_herit_batch_ret_INVALID_FORMAT_3;
	  }
	}
      } else {
	row_ptr = &(rel_base[ulii * pheno_nm_ct]);
	for (sample_uidx2 = 0; sample_uidx2 <= sample_uidx; sample_uidx2++) {
	  if (!gzgets(grm_gzfile, g_textbuf, MAXLINELEN)) {
	    goto unrelated_herit_batch_ret_READ_FAIL;
	  }
	  if (!g_textbuf[MAXLINELEN - 1]) {
	    goto unrelated_herit_batch_ret_INVALID_FORMAT_3;
	  }
	  if (IS_SET(pheno_nm, sample_uidx2)) {
	    bufptr = next_token_mult(g_textbuf, 3);
	    if (no_more_tokens_kns(bufptr)) {
	      goto unrelated_herit_batch_ret_INVALID_FORMAT_3;
	    }
	    if (scan_double(bufptr, &dxx)) {
	      goto unrelated_herit_batch_ret_INVALID_FORMAT_3;
	    }
            *row_ptr++ = dxx;
	  }
	}
        ulii++;
      }
    }
    gzclose(grm_gzfile);
    grm_gzfile = nullptr;
  }
  // fill in upper right
  for (ulii = 0; ulii < pheno_nm_ct; ulii++) {
    row_ptr = &(rel_base[ulii * pheno_nm_ct + ulii + 1]);
    for (uljj = ulii + 1; uljj < pheno_nm_ct; uljj++) {
      *row_ptr++ = rel_base[uljj * pheno_nm_ct + ulii];
    }
  }
  LOGPRINTF("--unrelated-heritability: %" PRIuPTR " phenotypes loaded.\n", pheno_nm_ct);
  reml_em_one_trait(matrix_wkbase, pheno_ptr, &unrelated_herit_covg, &unrelated_herit_covr, unrelated_herit_tol, is_strict);
  LOGPRINTF("h^2 estimate: %g\n", unrelated_herit_covg);
  while (0) {
  unrelated_herit_batch_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  unrelated_herit_batch_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  unrelated_herit_batch_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  unrelated_herit_batch_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  unrelated_herit_batch_ret_INVALID_FORMAT_3:
    logerrprint("Error: Invalid .grm.gz file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  unrelated_herit_batch_ret_INVALID_FORMAT_2:
    logerrprintb();
  unrelated_herit_batch_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 unrelated_herit_batch_ret_1:
  aligned_free_cond(pheno_c);
  fclose_cond(infile);
  fclose_cond(grm_binfile);
  gzclose_cond(grm_gzfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return retval;
}
#endif

int32_t ibs_test_calc(pthread_t* threads, char* read_dists_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t perm_ct, uintptr_t pheno_nm_ct, uintptr_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t pheno_nm_ctl = BITCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t perm_ctcl = 1 + (perm_ct / CACHELINE_BIT);
  uintptr_t perm_ctclm = perm_ctcl * (CACHELINE / sizeof(intptr_t));
  uintptr_t perm_ctcld = 1 + (perm_ct / CACHELINE_DBL);
  uintptr_t perm_ctcldm = perm_ctcld * CACHELINE_DBL;
  uintptr_t case_ct = pheno_nm_ct - pheno_ctrl_ct;
  uint32_t tidx = 1;
  int32_t retval = 0;
  int32_t perm_test[6];
  uintptr_t* perm_rows;
  double* perm_results;
  double ctrl_ctrl_ct;
  double ctrl_case_ct;
  double case_case_ct;
  double tot_sum;
  double ctrl_ctrl_tot;
  double ctrl_case_tot;
  double case_case_tot;
  double ctrl_ctrl_ssq;
  double ctrl_case_ssq;
  double case_case_ssq;
  double tot_mean;
  double ingroups_mean;
  double ctrl_ctrl_mean;
  double ctrl_case_mean;
  double case_case_mean;
  double ctrl_ctrl_var;
  double ctrl_case_var;
  double case_case_var;
  double ctrl_ctrl_tot1;
  double ctrl_case_tot1;
  double case_case_tot1;
  double case_case_minus_ctrl_ctrl;
  double case_case_minus_ctrl_case;
  double ctrl_ctrl_minus_ctrl_case;
  double between_ssq;
  double total_ssq;
  double perm_ct_recip;
  uintptr_t ulii;
  uintptr_t uljj = 0;
#ifdef __LP64__
  __m128d* rvptr1;
  __m128d* rvptr2;
#else
  double* rptr1;
  double* rptr2;
#endif
  uintptr_t perm_idx;
  g_load_dists = read_dists_fname? 1 : 0;
  g_sample_ct = sample_ct;
  perm_ct += 1; // first permutation = original config
  if (pheno_ctrl_ct < 2) {
    logerrprint("Warning: Skipping --ibs-test due to too few controls (minimum 2).\n");
    goto ibs_test_calc_ret_1;
  } else if (case_ct < 2) {
    logerrprint("Warning: Skipping --ibs-test due to too few cases (minimum 2).\n");
    goto ibs_test_calc_ret_1;
  }
  for (ulii = 0; ulii < 6; ulii++) {
    perm_test[ulii] = 0;
  }
  ctrl_ctrl_ct = (pheno_ctrl_ct * (pheno_ctrl_ct - 1)) / 2;
  ctrl_case_ct = pheno_ctrl_ct * case_ct;
  case_case_ct = (case_ct * (case_ct - 1)) / 2;
  g_perm_ct = perm_ct;
  // g_pheno_nm and g_pheno_c should be nullptr
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &g_pheno_nm) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, &g_pheno_c)) {
    goto ibs_test_calc_ret_NOMEM;
  }
  copy_bitarr_subset_excl(pheno_nm, sample_exclude, unfiltered_sample_ct, sample_ct, g_pheno_nm);
  copy_bitarr_subset_excl(pheno_c, sample_exclude, unfiltered_sample_ct, sample_ct, g_pheno_c);
  if (bigstack_alloc_d(g_thread_ct * 32 * BITCT, &g_ibs_test_partial_sums) ||
      bigstack_alloc_ul(perm_ct * pheno_nm_ctl, &perm_rows) ||
      bigstack_alloc_ul(perm_ctclm * g_thread_ct, &g_perm_col_buf) ||
      bigstack_calloc_d(2 * perm_ctcldm * g_thread_ct, &perm_results)) {
    goto ibs_test_calc_ret_NOMEM;
  }
  g_perm_results = perm_results;
  g_perm_rows = perm_rows;

  // first permutation = original
  copy_bitarr_subset(g_pheno_c, g_pheno_nm, unfiltered_sample_ct, pheno_nm_ct, perm_rows);
  for (ulii = pheno_nm_ctl - 1; ulii; ulii--) {
    perm_rows[ulii * perm_ct] = perm_rows[ulii];
  }

  printf("--ibs-test (%" PRIuPTR " permutations): [generating permutations]", perm_ct - 1);
  fflush(stdout);
  // minor todo: multithread this
  // less minor: cluster support
  generate_perm1_interleaved(pheno_nm_ct, case_ct, 1, perm_ct, perm_rows);
  fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                       \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b0%", stdout);
  fflush(stdout);
  for (ulii = 0; ulii < pheno_nm_ct; ulii++) {
    uljj += ((perm_rows[((ulii / BITCT) * perm_ct)] >> (ulii & (BITCT - 1))) & 1);
  }
  triangle_fill(pheno_nm_ct, g_thread_ct, 0, 1, 1, 1, g_thread_start);
  if (spawn_threads(threads, &ibs_test_thread, g_thread_ct)) {
    goto ibs_test_calc_ret_THREAD_CREATE_FAIL;
  }
  ulii = 0;
  ibs_test_thread((void*)ulii);
  join_threads(threads, g_thread_ct);

  tot_sum = g_calc_result[0][0];
  ctrl_ctrl_ssq = g_calc_result[0][1];
  ctrl_case_ssq = g_calc_result[0][2];
  case_case_ssq = g_calc_result[0][3];
  for (; tidx < g_thread_ct; tidx++) {
    tot_sum += g_calc_result[tidx][0];
    ctrl_ctrl_ssq += g_calc_result[tidx][1];
    ctrl_case_ssq += g_calc_result[tidx][2];
    case_case_ssq += g_calc_result[tidx][3];
#ifdef __LP64__
    rvptr1 = (__m128d*)perm_results;
    rvptr2 = (__m128d*)(&(perm_results[2 * perm_ctcldm * tidx]));
    for (perm_idx = 0; perm_idx < perm_ct; perm_idx++) {
      *rvptr1 = _mm_add_pd(*rvptr1, *rvptr2++);
      rvptr1++;
    }
#else
    rptr1 = perm_results;
    rptr2 = &(perm_results[2 * perm_ctcldm * tidx]);
    for (perm_idx = 0; perm_idx < perm_ct; perm_idx++) {
      *rptr1++ += *rptr2++;
      *rptr1++ += *rptr2++;
    }
#endif
  }
  ctrl_ctrl_tot = perm_results[0];
  ctrl_case_tot = perm_results[1];
  case_case_tot = tot_sum - ctrl_ctrl_tot - ctrl_case_tot;

  tot_mean = tot_sum / (ctrl_ctrl_ct + ctrl_case_ct + case_case_ct);
  ingroups_mean = (ctrl_ctrl_tot + case_case_tot) / (ctrl_ctrl_ct + case_case_ct);
  ctrl_ctrl_mean = ctrl_ctrl_tot / ctrl_ctrl_ct;
  ctrl_case_mean = ctrl_case_tot / ctrl_case_ct;
  case_case_mean = case_case_tot / case_case_ct;

  ctrl_ctrl_var = ctrl_ctrl_ssq - ctrl_ctrl_tot * ctrl_ctrl_mean;
  ctrl_case_var = ctrl_case_ssq - ctrl_case_tot * ctrl_case_mean;
  case_case_var = case_case_ssq - case_case_tot * case_case_mean;

  total_ssq = ctrl_ctrl_var + ctrl_case_var + case_case_var;
  between_ssq = ctrl_case_ct * (ctrl_case_mean - tot_mean) * (ctrl_case_mean - tot_mean) + (ctrl_ctrl_ct + case_case_ct) * (ingroups_mean - tot_mean) * (ingroups_mean - tot_mean);

  case_case_minus_ctrl_ctrl = case_case_tot - ctrl_ctrl_tot;
  case_case_minus_ctrl_case = case_case_tot - ctrl_case_tot;
  ctrl_ctrl_minus_ctrl_case = ctrl_ctrl_tot - ctrl_case_tot;

  for (ulii = 1; ulii < perm_ct; ulii++) {
    ctrl_ctrl_tot1 = perm_results[ulii * 2];
    ctrl_case_tot1 = perm_results[ulii * 2 + 1];
    case_case_tot1 = tot_sum - ctrl_ctrl_tot1 - ctrl_case_tot1;
    if (ctrl_case_tot1 < ctrl_case_tot) {
      perm_test[0] += 1;
    }
    if (case_case_tot1 - ctrl_ctrl_tot1 < case_case_minus_ctrl_ctrl) {
      perm_test[1] += 1;
    }
    if (case_case_tot1 < case_case_tot) {
      perm_test[2] += 1;
    }
    if (ctrl_ctrl_tot1 < ctrl_ctrl_tot) {
      perm_test[3] += 1;
    }
    if (case_case_tot1 - ctrl_case_tot1 < case_case_minus_ctrl_case) {
      perm_test[4] += 1;
    }
    if (ctrl_ctrl_tot1 - ctrl_case_tot1 < ctrl_ctrl_minus_ctrl_case) {
      perm_test[5] += 1;
    }
  }

  fputs("\r                                         \r", stdout);
  logprint("--ibs-test results:\n");
  LOGPRINTF("  Between-group IBS (mean, SD)   = %g, %g\n", ctrl_case_mean, sqrt(ctrl_case_var / (ctrl_case_ct - 1)));
  LOGPRINTF("  In-group (case) IBS (mean, SD) = %g, %g\n", case_case_mean, sqrt(case_case_var / (case_case_ct - 1)));
  LOGPRINTF("  In-group (ctrl) IBS (mean, SD) = %g, %g\n", ctrl_ctrl_mean, sqrt(ctrl_ctrl_var / (ctrl_ctrl_ct - 1)));
  LOGPRINTF("  Approximate proportion of variance between group = %g\n", between_ssq / total_ssq);
  perm_ct_recip = 1.0 / ((double)perm_ct);
  fputs("  IBS group-difference empirical p-values:\n", stdout);
  LOGPRINTF("     T1: Case/control less similar                p = %g\n", (perm_test[0]) * perm_ct_recip);
  LOGPRINTF("     T2: Case/control more similar                p = %g\n\n", (perm_ct - perm_test[0]) * perm_ct_recip);
  LOGPRINTF("     T3: Case/case less similar than ctrl/ctrl    p = %g\n", (perm_test[1]) * perm_ct_recip);
  LOGPRINTF("     T4: Case/case more similar than ctrl/ctrl    p = %g\n\n", (perm_ct - perm_test[1]) * perm_ct_recip);
  LOGPRINTF("     T5: Case/case less similar                   p = %g\n", (perm_test[2]) * perm_ct_recip);
  LOGPRINTF("     T6: Case/case more similar                   p = %g\n\n", (perm_ct - perm_test[2]) * perm_ct_recip);
  LOGPRINTF("     T7: Control/control less similar             p = %g\n", (perm_test[3]) * perm_ct_recip);
  LOGPRINTF("     T8: Control/control more similar             p = %g\n\n", (perm_ct - perm_test[3]) * perm_ct_recip);
  LOGPRINTF("     T9: Case/case less similar than case/ctrl    p = %g\n", (perm_test[4]) * perm_ct_recip);
  LOGPRINTF("    T10: Case/case more similar than case/ctrl    p = %g\n\n", (perm_ct - perm_test[4]) * perm_ct_recip);
  LOGPRINTF("    T11: Ctrl/ctrl less similar than case/ctrl    p = %g\n", (perm_test[5]) * perm_ct_recip);
  LOGPRINTF("    T12: Ctrl/ctrl more similar than case/ctrl    p = %g\n", (perm_ct - perm_test[5]) * perm_ct_recip);

  while (0) {
  ibs_test_calc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ibs_test_calc_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 ibs_test_calc_ret_1:
  bigstack_reset(bigstack_mark);
  g_pheno_nm = nullptr;
  g_pheno_c = nullptr;
  return retval;
}

int32_t groupdist_calc(pthread_t* threads, uint32_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t groupdist_iters, uint32_t groupdist_d, uint32_t pheno_nm_ct, uint32_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  double* dist_ptr = g_dists;
  double dhh_ssq = 0.0;
  double dhl_ssq = 0.0;
  double dll_ssq = 0.0;
  int32_t retval = 0;
  uintptr_t* pheno_nm_local;
  uintptr_t* pheno_c_local;
  double* jackknife_precomp;
  int32_t ll_size;
  int32_t lh_size;
  int32_t hh_size;
  double* ll_pool;
  double* lh_pool;
  double* hh_pool;
  double* ll_poolp;
  double* lh_poolp;
  double* hh_poolp;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t sample_idx;
  double ll_med;
  double lh_med;
  double hh_med;
  double dll_sd;
  double dhl_sd;
  double dhh_sd;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  uint32_t is_case;
  if (pheno_ctrl_ct < 2) {
    logerrprint("Warning: Skipping --groupdist due to too few controls (minimum 2).\n");
    goto groupdist_calc_ret_1;
  }
  g_case_ct = pheno_nm_ct - pheno_ctrl_ct;
  if (g_case_ct < 2) {
    logerrprint("Warning: Skipping --groupdist due to too few cases (minimum 2).\n");
    goto groupdist_calc_ret_1;
  }
  g_ctrl_ct = pheno_ctrl_ct;
  g_sample_ct = sample_ct;
  // g_pheno_nm and g_pheno_c should be nullptr
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &pheno_nm_local) ||
      bigstack_alloc_ul(unfiltered_sample_ctl, &pheno_c_local)) {
    goto groupdist_calc_ret_NOMEM;
  }
  g_pheno_nm = pheno_nm_local;
  g_pheno_c = pheno_c_local;
  copy_bitarr_subset_excl(pheno_nm, sample_exclude, unfiltered_sample_ct, sample_ct, pheno_nm_local);
  copy_bitarr_subset_excl(pheno_c, sample_exclude, unfiltered_sample_ct, sample_ct, pheno_c_local);
  ll_size = ((uintptr_t)g_ctrl_ct * (g_ctrl_ct - 1)) / 2;
  lh_size = g_ctrl_ct * g_case_ct;
  hh_size = ((uintptr_t)g_case_ct * (g_case_ct - 1)) / 2;
  g_reg_tot_y = 0.0;
  g_reg_tot_xy = 0.0;
  g_reg_tot_x = 0.0;
  if (groupdist_d) {
    g_jackknife_d = groupdist_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(g_case_ct + g_ctrl_ct);
  }
  if (bigstack_alloc_d(ll_size, &ll_pool) ||
      bigstack_alloc_d(lh_size, &lh_pool) ||
      bigstack_alloc_d(hh_size, &hh_pool) ||
      bigstack_alloc_uc(g_thread_ct * round_up_pow2(g_case_ct + g_ctrl_ct + (g_jackknife_d + 1) * sizeof(int32_t), CACHELINE), &g_geno)) {
    goto groupdist_calc_ret_NOMEM;
  }
  ll_poolp = ll_pool;
  lh_poolp = lh_pool;
  hh_poolp = hh_pool;
  for (sample_idx = 1; sample_idx < sample_ct; sample_idx++) {
    if (IS_SET(pheno_nm_local, sample_idx)) {
      if (IS_SET(pheno_c_local, sample_idx)) {
	for (uii = 0; uii < sample_idx; uii++) {
	  if (IS_SET(pheno_nm_local, uii)) {
	    dxx = *dist_ptr;
	    if (IS_SET(pheno_c_local, uii)) {
	      *hh_poolp++ = dxx;
	      g_reg_tot_x += dxx;
	      dhh_ssq += dxx * dxx;
	    } else {
	      *lh_poolp++ = dxx;
	      g_reg_tot_xy += dxx;
	      dhl_ssq += dxx * dxx;
	    }
	  }
	  dist_ptr++;
	}
      } else {
	for (uii = 0; uii < sample_idx; uii++) {
	  if (IS_SET(pheno_nm_local, uii)) {
	    dxx = *dist_ptr;
	    if (IS_SET(pheno_c_local, uii)) {
	      *lh_poolp++ = dxx;
	      g_reg_tot_xy += dxx;
	      dhl_ssq += dxx * dxx;
	    } else {
	      *ll_poolp++ = dxx;
	      g_reg_tot_y += dxx;
	      dll_ssq += dxx * dxx;
	    }
	  }
	  dist_ptr++;
	}
      }
    } else {
      dist_ptr += sample_idx;
    }
  }
  ll_med = destructive_get_dmedian(ll_size, ll_pool);
  lh_med = destructive_get_dmedian(lh_size, lh_pool);
  hh_med = destructive_get_dmedian(hh_size, hh_pool);
  logprint("Case/control distance analysis:\n");
  if (g_case_ct < 2) {
    dxx = 0.0;
    dhh_sd = 0.0;
  } else {
    dww = (double)(((uintptr_t)g_case_ct * (g_case_ct - 1)) / 2);
    dxx = g_reg_tot_x / dww;
    dhh_sd = sqrt((dhh_ssq / dww - dxx * dxx) / (dww - 1.0));
  }
  if (!(g_case_ct && g_ctrl_ct)) {
    dyy = 0.0;
    dhl_sd = 0.0;
  } else {
    dww = (double)((uintptr_t)g_case_ct * g_ctrl_ct);
    dyy = g_reg_tot_xy / dww;
    dhl_sd = sqrt((dhl_ssq / dww - dyy * dyy) / (dww - 1.0));
  }
  if (g_ctrl_ct < 2) {
    dzz = 0.0;
    dll_sd = 0.0;
  } else {
    dww = (double)(((uintptr_t)g_ctrl_ct * (g_ctrl_ct - 1)) / 2);
    dzz = g_reg_tot_y / dww;
    dll_sd = sqrt((dll_ssq / dww - dzz * dzz) / (dww - 1.0));
  }
  LOGPRINTF("  Mean (sd), median dists between 2x affected     : %g (%g), %g\n", dxx, dhh_sd, hh_med);
  LOGPRINTF("  Mean (sd), median dists between aff. and unaff. : %g (%g), %g\n", dyy, dhl_sd, lh_med);
  LOGPRINTF("  Mean (sd), median dists between 2x unaffected   : %g (%g), %g\n\n", dzz, dll_sd, ll_med);
  if (2 * g_jackknife_d >= (g_case_ct + g_ctrl_ct)) {
    logprint("Delete-d jackknife skipped because d is too large.\n");
  } else {
    if (bigstack_calloc_d(sample_ct * JACKKNIFE_VALS_GROUPDIST, &jackknife_precomp)) {
      goto groupdist_calc_ret_NOMEM;
    }
    g_jackknife_precomp = jackknife_precomp;
    if (bigstack_init_sfmtp(g_thread_ct)) {
      goto groupdist_calc_ret_NOMEM;
    }
    // to precompute:
    // tot_uu, tot_au, tot_aa
    dist_ptr = g_dists;
    for (sample_idx = 1; sample_idx < sample_ct; sample_idx++) {
      if (IS_SET(pheno_nm_local, sample_idx)) {
	uii = 0;
	is_case = IS_SET(pheno_c_local, sample_idx);
	dyy = 0;
	dzz = 0;
	do {
	  if (IS_SET(pheno_nm_local, uii)) {
	    dxx = dist_ptr[uii];
	    if (IS_SET(pheno_c_local, uii)) {
	      jackknife_precomp[(uii * JACKKNIFE_VALS_GROUPDIST) + is_case + 1] += dxx;
	      dzz += dxx;
	    } else {
	      jackknife_precomp[(uii * JACKKNIFE_VALS_GROUPDIST) + is_case] += dxx;
	      dyy += dxx;
	    }
	  }
	} while ((++uii) < sample_idx);
	jackknife_precomp[(sample_idx * JACKKNIFE_VALS_GROUPDIST) + is_case] += dyy;
	jackknife_precomp[(sample_idx * JACKKNIFE_VALS_GROUPDIST) + is_case + 1] += dzz;
      }
      dist_ptr = &(dist_ptr[sample_idx]);
    }

    g_jackknife_iters = (groupdist_iters + g_thread_ct - 1) / g_thread_ct;
    if (spawn_threads(threads, &groupdist_jack_thread, g_thread_ct)) {
      goto groupdist_calc_ret_THREAD_CREATE_FAIL;
    }
    ulii = 0;
    groupdist_jack_thread((void*)ulii);
    join_threads(threads, g_thread_ct);
    for (uii = 1; uii < g_thread_ct; uii++) {
      for (ujj = 0; ujj < 9; ujj++) {
	g_calc_result[0][ujj] += g_calc_result[uii][ujj];
      }
    }
    dxx = 1.0 / ((double)(g_jackknife_iters * g_thread_ct));
    // avg[X-Y] = avg(X) - avg(Y)
    // Jackknife SE(X-Y) = sqrt(((obs_ct - jackknife_d) / jackknife_d) *
    //                          avg(([actual X-Y] - [avg X-Y])^2))
    //                   = sqrt(((obs_ct - jackknife_d) / jackknife_d) *
    //                          avg(([actual X-Y])^2 - [avg X-Y]^2))
    //                   = sqrt(((obs_ct - jackknife_d) / jackknife_d) *
    //                          (avg([actual X-Y]^2) - [avg X-Y]^2))
    //                   = sqrt(((obs_ct - jackknife_d) / jackknife_d) *
    //                          (avg(X^2) - avg(2XY) + avg(Y^2) - avg[X-Y]^2))
    for (uii = 0; uii < 9; uii++) {
      g_calc_result[0][uii] *= dxx;
    }
    putc_unlocked('\r', stdout);
    dxx = g_calc_result[0][0] - g_calc_result[0][1];
    LOGPRINTF("  AA mean - AU mean avg difference (s.e.): %g (%g)\n", dxx, sqrt(((g_case_ct + g_ctrl_ct - g_jackknife_d) / ((double)g_jackknife_d)) * (g_calc_result[0][3] + g_calc_result[0][4] - 2 * g_calc_result[0][6] - dxx * dxx)));
    dxx = g_calc_result[0][0] - g_calc_result[0][2];
    LOGPRINTF("  AA mean - UU mean avg difference (s.e.): %g (%g)\n", dxx, sqrt(((g_case_ct + g_ctrl_ct - g_jackknife_d) / ((double)g_jackknife_d)) * (g_calc_result[0][3] + g_calc_result[0][5] - 2 * g_calc_result[0][7] - dxx * dxx)));
    dxx = g_calc_result[0][1] - g_calc_result[0][2];
    LOGPRINTF("  AU mean - UU mean avg difference (s.e.): %g (%g)\n", dxx, sqrt(((g_case_ct + g_ctrl_ct - g_jackknife_d) / ((double)g_jackknife_d)) * (g_calc_result[0][4] + g_calc_result[0][5] - 2 * g_calc_result[0][8] - dxx * dxx)));
  }
  while (0) {
  groupdist_calc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  groupdist_calc_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 groupdist_calc_ret_1:
  bigstack_reset(bigstack_mark);
  g_pheno_nm = nullptr;
  g_pheno_c = nullptr;
  return retval;
}

void normalize_phenos(double* new_phenos, uint32_t sample_ct, uintptr_t* sample_exclude, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t sex_exclude) {
  uint32_t incl_males = sex_exclude & 1;
  uint32_t incl_females = sex_exclude & 2;
  uint32_t incl_unknown = sex_exclude & 4;
  uint32_t sample_uidx = 0xffffffffU; // deliberate overflow
  double pheno_tot = 0.0;
  double pheno_sq_tot = 0.0;
  uint32_t pheno_ct = 0;
  uint32_t sample_idx = 0xffffffffU;
  double dxx;
  double mean;
  double stdev_recip;

  while ((++sample_idx) < sample_ct) {
    sample_uidx++;
    next_unset_unsafe_ck(sample_exclude, &sample_uidx);
    if (IS_SET(sex_nm, sample_uidx)) {
      if (IS_SET(sex_male, sample_uidx)) {
	if (!incl_males) {
	  continue;
	}
      } else if (!incl_females) {
	continue;
      }
    } else if (!incl_unknown) {
      continue;
    }
    dxx = new_phenos[sample_idx];
    pheno_tot += dxx;
    pheno_sq_tot += dxx * dxx;
    pheno_ct++;
  }
  if (!pheno_ct) {
    return;
  }
  mean = pheno_tot / pheno_ct;
  if (pheno_ct == 1) {
    stdev_recip = 0;
  } else {
    stdev_recip = sqrt((double)(pheno_ct - 1) / (pheno_sq_tot - pheno_tot * mean));
  }
  sample_uidx = 0xffffffffU;
  sample_idx = 0xffffffffU;
  while ((++sample_idx) < sample_ct) {
    sample_uidx++;
    next_unset_unsafe_ck(sample_exclude, &sample_uidx);
    if (IS_SET(sex_nm, sample_uidx)) {
      if (IS_SET(sex_male, sample_uidx)) {
	if (!incl_males) {
	  continue;
	}
      } else if (!incl_females) {
	continue;
      }
    } else if (!incl_unknown) {
      continue;
    }
    new_phenos[sample_idx] = (new_phenos[sample_idx] - mean) * stdev_recip;
  }
}

static uintptr_t g_cg_sample1idx;
static uintptr_t g_cg_sample2idx;
static uintptr_t g_cg_sample1uidx;
static uintptr_t g_cg_sample2uidx;
static uintptr_t g_cg_tstc;
static uint32_t g_cg_marker_ct;
static uintptr_t g_cg_max_sample_id_len;
static uint32_t g_cg_max_sample_fid_len;
static uint32_t g_cg_max_sample_iid_len;
static uintptr_t g_cg_max_paternal_id_len;
static uintptr_t g_cg_max_maternal_id_len;
static uintptr_t* g_cg_sample_exclude;
static uintptr_t* g_cg_pheno_nm;
static uintptr_t* g_cg_pheno_c;
static uintptr_t* g_cg_founder_info;
static char* g_cg_sample_ids;
static char* g_cg_paternal_ids;
static char* g_cg_maternal_ids;
static uint32_t g_cg_genome_modifier;
static Pedigree_rel_info g_cg_pri;
static uint32_t g_cg_family_id_fixed;
static uint32_t g_cg_is_founder_fixed;
static uint32_t g_cg_rel_space_id_fixed;
static int64_t g_cg_llfct;
static char* g_cg_sptr_start;
static char* g_cg_fam1;
static char* g_cg_fam2;
static char* g_cg_sample1;
static char* g_cg_pat1;
static char* g_cg_mat1;
static int64_t g_cg_cur_line;
static int64_t g_cg_tot_lines;
static uint32_t g_pct;
static uintptr_t g_cg_gmcell;
static uintptr_t g_cg_mdecell;
static double g_cg_e00;
static double g_cg_e01;
static double g_cg_e02;
static double g_cg_e11;
static double g_cg_e12;
static double g_cg_min_pi_hat;
static double g_cg_max_pi_hat;

static uintptr_t g_dw_sample1idx;
static uintptr_t g_dw_sample2idx;
static uintptr_t g_dw_min_sample;
static uintptr_t g_dw_max_sample1idx;
static uint64_t g_dw_start_offset;
static uint64_t g_dw_hundredth;
static double* g_dw_dists;
static double* g_dw_dist_ptr;
static unsigned char* g_dw_membuf;
static double g_dw_half_marker_ct_recip;

int32_t write_ids(char* outname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, char* sample_ids, uintptr_t max_sample_id_len) {
  uintptr_t sample_uidx = 0;
  FILE* outfile;
  if (fopen_checked(outname, "w", &outfile)) {
    return RET_OPEN_FAIL;
  }
  while (1) {
    next_unset_ul_ck(sample_exclude, unfiltered_sample_ct, &sample_uidx);
    if (sample_uidx == unfiltered_sample_ct) {
      break;
    }
    fputs(&(sample_ids[sample_uidx * max_sample_id_len]), outfile);
    if (putc_checked('\n', outfile)) {
      return RET_WRITE_FAIL;
    }
    sample_uidx++;
  }
  if (fclose_null(&outfile)) {
    return RET_WRITE_FAIL;
  }
  return 0;
}

int32_t distance_d_write_ids(char* outname, char* outname_end, uint32_t dist_calc_type, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, char* sample_ids, uintptr_t max_sample_id_len) {
  int32_t retval;
  if (dist_calc_type & DISTANCE_ALCT) {
    strcpy(outname_end, ".dist.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      return retval;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      return retval;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      return retval;
    }
  }
  return 0;
}

int32_t distance_open(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, char* outname, char* outname_end, const char* varsuffix, const char* mode, int32_t dist_calc_type, int32_t parallel_idx, int32_t parallel_tot) {
  if (dist_calc_type & DISTANCE_ALCT) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".dist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".dist%s", varsuffix);
    }
    strcpy(g_textbuf, outname_end);
    if (fopen_checked(outname, mode, outfile_ptr)) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mibs%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mibs%s", varsuffix);
    }
    strcpy(&(g_textbuf[MAX_POST_EXT]), outname_end);
    if (fopen_checked(outname, mode, outfile2_ptr)) {
      return 1;
    }
  }
  if (dist_calc_type & DISTANCE_1_MINUS_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mdist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mdist%s", varsuffix);
    }
    strcpy(&(g_textbuf[MAX_POST_EXT * 2]), outname_end);
    if (fopen_checked(outname, mode, outfile3_ptr)) {
      return 1;
    }
  }
  return 0;
}

void distance_print_done(int32_t format_code, char* outname, char* outname_end) {
  putc_unlocked('\r', stdout);
  if (!format_code) {
    strcpy(outname_end, g_textbuf);
    sprintf(g_logbuf, "Distances (allele counts) written to %s .\n", outname);
  } else if (format_code == 1) {
    strcpy(outname_end, &(g_textbuf[MAX_POST_EXT]));
    sprintf(g_logbuf, "IBS matrix written to %s .\n", outname);
  } else if (format_code == 2) {
    strcpy(outname_end, &(g_textbuf[MAX_POST_EXT * 2]));
    sprintf(g_logbuf, "Distances (proportions) written to %s .\n", outname);
  }
  wordwrapb(0);
  logprintb();
}

uint32_t distance_d_write_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uint64_t start_offset = g_dw_start_offset;
  uint64_t hundredth = g_dw_hundredth;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  double* dist_ptr = g_dw_dist_ptr;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx + 1 < sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_tri_emitn_ret;
      }
    }
    if (sample2idx + 1 == sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\n', sptr_cur);
    }
    if ((((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) >= hundredth * pct) {
      pct = (((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) / hundredth;
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample1idx++;
    sample2idx = 0;
  }
 distance_d_write_tri_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  uintptr_t min_sample = g_dw_min_sample;
  double* dist_ptr = g_dw_dist_ptr;
  unsigned char* membuf1 = &(g_dw_membuf[1]);
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  uintptr_t ulii;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_sq0_emitn_ret;
      }
    }
    ulii = 1 + (((uintptr_t)(readbuf_end - sptr_cur)) / 2);
    if (ulii < (sample_ct - sample2idx)) {
      sptr_cur = memcpya(sptr_cur, membuf1, 2 * ulii);
      sample2idx += ulii;
      goto distance_d_write_sq0_emitn_ret;
    }
    ulii = sample_ct - sample2idx;
    sptr_cur = memcpyax(sptr_cur, membuf1, 2 * ulii - 1, '\n');
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_sq0_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  uintptr_t min_sample = g_dw_min_sample;
  double* dists = g_dw_dists;
  double* dist_ptr = g_dw_dist_ptr;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_sq_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      *sptr_cur++ = '0';
      sample2idx++;
    }
    while (sample2idx < sample_ct) {
      *sptr_cur++ = '\t';
      sptr_cur = dtoa_g(dists[((sample2idx * (sample2idx - 1)) / 2) + sample1idx], sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_sq_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_ibs_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uint64_t start_offset = g_dw_start_offset;
  uint64_t hundredth = g_dw_hundredth;
  double half_marker_ct_recip = g_dw_half_marker_ct_recip;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  double* dist_ptr = g_dw_dist_ptr;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(1.0 - (*dist_ptr++) * half_marker_ct_recip, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_tri_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      sptr_cur = memcpya(sptr_cur, "1\n", 2);
    }
    sample1idx++;
    if ((((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) >= hundredth * pct) {
      pct = (((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) / hundredth;
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_ibs_tri_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_ibs_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  double half_marker_ct_recip = g_dw_half_marker_ct_recip;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  uintptr_t min_sample = g_dw_min_sample;
  double* dist_ptr = g_dw_dist_ptr;
  unsigned char* membuf = g_dw_membuf;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  uintptr_t ulii;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(1.0 - (*dist_ptr++) * half_marker_ct_recip, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq0_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      *sptr_cur++ = '1';
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq0_emitn_ret;
      }
    }
    ulii = (1 + ((uintptr_t)(readbuf_end - sptr_cur))) / 2;
    if (ulii < (sample_ct - sample2idx)) {
      sptr_cur = memcpya(sptr_cur, membuf, 2 * ulii);
      sample2idx += ulii;
      goto distance_d_write_ibs_sq0_emitn_ret;
    }
    ulii = sample_ct - sample2idx;
    sptr_cur = memcpyax(sptr_cur, membuf, 2 * ulii, '\n');
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_ibs_sq0_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_ibs_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  double half_marker_ct_recip = g_dw_half_marker_ct_recip;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  uintptr_t min_sample = g_dw_min_sample;
  double* dists = g_dw_dists;
  double* dist_ptr = g_dw_dist_ptr;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(1.0 - (*dist_ptr++) * half_marker_ct_recip, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      *sptr_cur++ = '1';
      sample2idx++;
    }
    while (sample2idx < sample_ct) {
      *sptr_cur++ = '\t';
      sptr_cur = dtoa_g(1.0 - (dists[((sample2idx * (sample2idx - 1)) / 2) + sample1idx]) * half_marker_ct_recip, sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_ibs_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_ibs_sq_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_1mibs_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uint64_t start_offset = g_dw_start_offset;
  uint64_t hundredth = g_dw_hundredth;
  double half_marker_ct_recip = g_dw_half_marker_ct_recip;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  double* dist_ptr = g_dw_dist_ptr;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx + 1 < sample1idx) {
      sptr_cur = dtoa_gx((*dist_ptr++) * half_marker_ct_recip, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_tri_emitn_ret;
      }
    }
    if (sample2idx + 1 == sample1idx) {
      sptr_cur = dtoa_gx((*dist_ptr++) * half_marker_ct_recip, '\n', sptr_cur);
    }
    if ((((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) >= hundredth * pct) {
      pct = (((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) / hundredth;
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample1idx++;
    sample2idx = 0;
  }
 distance_d_write_1mibs_tri_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_1mibs_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  double half_marker_ct_recip = g_dw_half_marker_ct_recip;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  uintptr_t min_sample = g_dw_min_sample;
  double* dist_ptr = g_dw_dist_ptr;
  unsigned char* membuf1 = &(g_dw_membuf[1]);
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  uintptr_t ulii;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx((*dist_ptr++) * half_marker_ct_recip, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_sq0_emitn_ret;
      }
    }
    ulii = 1 + (((uintptr_t)(readbuf_end - sptr_cur)) / 2);
    if (ulii < (sample_ct - sample2idx)) {
      sptr_cur = memcpya(sptr_cur, membuf1, 2 * ulii);
      sample2idx += ulii;
      goto distance_d_write_1mibs_sq0_emitn_ret;
    }
    ulii = sample_ct - sample2idx;
    sptr_cur = memcpyax(sptr_cur, membuf1, 2 * ulii - 1, '\n');
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_1mibs_sq0_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t distance_d_write_1mibs_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  double half_marker_ct_recip = g_dw_half_marker_ct_recip;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_dw_max_sample1idx;
  uintptr_t min_sample = g_dw_min_sample;
  double* dists = g_dw_dists;
  double* dist_ptr = g_dw_dist_ptr;
  uintptr_t sample1idx = g_dw_sample1idx;
  uintptr_t sample2idx = g_dw_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx((*dist_ptr++) * half_marker_ct_recip, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_sq_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      *sptr_cur++ = '0';
      sample2idx++;
    }
    while (sample2idx < sample_ct) {
      *sptr_cur++ = '\t';
      sptr_cur = dtoa_g((dists[((sample2idx * (sample2idx - 1)) / 2) + sample1idx]) * half_marker_ct_recip, sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto distance_d_write_1mibs_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 distance_d_write_1mibs_sq_emitn_ret:
  g_dw_dist_ptr = dist_ptr;
  g_dw_sample1idx = sample1idx;
  g_dw_sample2idx = sample2idx;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t distance_d_write(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, int32_t dist_calc_type, char* outname, char* outname_end, double* dists, double half_marker_ct_recip, uint32_t sample_ct, int32_t first_sample_idx, int32_t end_sample_idx, int32_t parallel_idx, int32_t parallel_tot, unsigned char* membuf) {
  // membuf assumed to be of at least size sample_ct * 8.
  uint32_t bin4 = dist_calc_type & DISTANCE_BIN4;
  int32_t shape = dist_calc_type & DISTANCE_SHAPEMASK;
  int32_t write_alcts = dist_calc_type & DISTANCE_ALCT;
  int32_t write_ibs_matrix = dist_calc_type & DISTANCE_IBS;
  int32_t write_1mibs_matrix = dist_calc_type & DISTANCE_1_MINUS_IBS;
  int32_t retval = 0;
  unsigned char overflow_buf[262144];
  double dxx;
  double dyy;
  double* dist_ptr;
  uintptr_t* glptr;
  uintptr_t sample_idx_ct;
  uintptr_t ulii;
  uintptr_t uljj;
  float fxx;
  float fyy;
  uint32_t uii;
  int32_t ii;
  int32_t jj;
  g_sample_ct = sample_ct;
  g_dw_start_offset = ((int64_t)first_sample_idx * (first_sample_idx - 1)) / 2;
  sample_idx_ct = (uintptr_t)(((int64_t)end_sample_idx * (end_sample_idx - 1)) / 2 - g_dw_start_offset);
  g_dw_hundredth = 1 + (sample_idx_ct / 100);
  if (first_sample_idx == 1) {
    first_sample_idx = 0;
  }
  g_pct = 1;
  if (dist_calc_type & (DISTANCE_BIN | DISTANCE_BIN4)) {
    if (distance_open(outfile_ptr, outfile2_ptr, outfile3_ptr, outname, outname_end, ".bin", FOPEN_WB, dist_calc_type, parallel_idx, parallel_tot)) {
      goto distance_d_write_ret_OPEN_FAIL;
    }
    if (shape == DISTANCE_TRI) {
      if (write_alcts) {
	fputs("Writing...", stdout);
	fflush(stdout);
        if (!bin4) {
	  if (fwrite_checked(dists, sample_idx_ct * sizeof(double), *outfile_ptr)) {
	    goto distance_d_write_ret_WRITE_FAIL;
	  }
	} else {
	  for (ulii = 0; ulii < sample_idx_ct; ulii++) {
            fxx = dists[ulii];
	    if (fwrite_checked(&fxx, sizeof(float), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	}
	distance_print_done(0, outname, outname_end);
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	ulii = 0;
	do {
	  uljj = (((uint64_t)sample_idx_ct) * g_pct) / 100LL;
	  if (!bin4) {
	    for (; ulii < uljj; ulii++) {
	      dxx = 1.0 - (*dist_ptr++) * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  } else {
	    for (; ulii < uljj; ulii++) {
	      fxx = (float)(1.0 - (*dist_ptr++) * half_marker_ct_recip);
	      if (fwrite_checked(&fxx, sizeof(float), *outfile2_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  printf("\rWriting... %d%%", g_pct++);
	  fflush(stdout);
	} while (g_pct <= 100);
	distance_print_done(1, outname, outname_end);
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	ulii = 0;
	do {
	  uljj = (((uint64_t)sample_idx_ct) * g_pct) / 100LL;
	  if (!bin4) {
	    for (; ulii < uljj; ulii++) {
	      dxx = (*dist_ptr++) * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  } else {
	    for (; ulii < uljj; ulii++) {
	      fxx = (float)((*dist_ptr++) * half_marker_ct_recip);
	      if (fwrite_checked(&fxx, sizeof(float), *outfile3_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  printf("\rWriting... %d%%", g_pct++);
	  fflush(stdout);
	} while (g_pct <= 100);
	distance_print_done(2, outname, outname_end);
      }
    } else if (!bin4) {
      if (shape == DISTANCE_SQ0) {
	fill_double_zero(sample_ct, (double*)membuf);
      }
      if (write_alcts) {
	dxx = 0.0;
	dist_ptr = dists;
	for (ii = first_sample_idx; ii < end_sample_idx; ii++) {
	  if (fwrite_checked(dist_ptr, ii * sizeof(double), *outfile_ptr)) {
	    goto distance_d_write_ret_WRITE_FAIL;
	  }
	  dist_ptr = &(dist_ptr[(uint32_t)ii]);
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (sample_ct - ii) * sizeof(double), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix, no need to handle parallel case
	    if (fwrite_checked(&dxx, sizeof(double), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < sample_ct; ulii++) {
	      if (fwrite_checked(&(dists[(ulii * (ulii - 1)) / 2 + ii]), sizeof(double), *outfile_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_sample_idx) * 100LL >= (int64_t)g_pct * (end_sample_idx - first_sample_idx)) {
	    g_pct = ((ii - first_sample_idx) * 100LL) / (end_sample_idx - first_sample_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
	g_pct = 1;
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	dyy = 1.0;
	for (ii = first_sample_idx; ii < end_sample_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    dxx = 1.0 - (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (sample_ct - ii) * sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&dyy, sizeof(double), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < sample_ct; ulii++) {
	      dxx = 1.0 - dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile2_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_sample_idx) * 100LL >= (int64_t)g_pct * (end_sample_idx - first_sample_idx)) {
	    g_pct = ((ii - first_sample_idx) * 100LL) / (end_sample_idx - first_sample_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile2_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(1, outname, outname_end);
	g_pct = 1;
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	dyy = 0.0;
	for (ii = first_sample_idx; ii < end_sample_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    dxx = (*dist_ptr++) * half_marker_ct_recip;
	    if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (sample_ct - ii) * sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&dyy, sizeof(double), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < sample_ct; ulii++) {
	      dxx = dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip;
	      if (fwrite_checked(&dxx, sizeof(double), *outfile3_ptr)) {
		goto distance_d_write_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((ii - first_sample_idx) * 100LL >= (int64_t)g_pct * (end_sample_idx - first_sample_idx)) {
	    g_pct = ((ii - first_sample_idx) * 100LL) / (end_sample_idx - first_sample_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile3_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(2, outname, outname_end);
      }
    } else {
      if (shape == DISTANCE_SQ0) {
	fill_float_zero(sample_ct, (float*)membuf);
      }
      if (write_alcts) {
	fxx = 0.0;
	dist_ptr = dists;
	for (ii = first_sample_idx; ii < end_sample_idx; ii++) {
	  for (uii = 0; uii < (uint32_t)ii; uii++) {
	    fxx = (float)dist_ptr[uii];
	    fwrite(&fxx, 4, 1, *outfile_ptr);
	  }
	  dist_ptr = &(dist_ptr[(uint32_t)ii]);
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (sample_ct - ii) * sizeof(float), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix, no need to handle parallel case
	    if (fwrite_checked(&fxx, sizeof(float), *outfile_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < sample_ct; ulii++) {
	      fyy = (float)dists[(ulii * (ulii - 1)) / 2 + ii];
	      fwrite(&fyy, 4, 1, *outfile_ptr);
	    }
	  }
	  if ((ii - first_sample_idx) * 100LL >= (int64_t)g_pct * (end_sample_idx - first_sample_idx)) {
	    g_pct = ((ii - first_sample_idx) * 100LL) / (end_sample_idx - first_sample_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
	g_pct = 1;
      }
      if (write_ibs_matrix) {
	dist_ptr = dists;
	fyy = 1.0;
	for (ii = first_sample_idx; ii < end_sample_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    fxx = (float)(1.0 - (*dist_ptr++) * half_marker_ct_recip);
	    if (fwrite_checked(&fxx, sizeof(float), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (sample_ct - ii) * sizeof(float), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&fyy, sizeof(float), *outfile2_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < sample_ct; ulii++) {
	      fxx = (float)(1.0 - dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip);
	      fwrite(&fxx, 4, 1, *outfile2_ptr);
	    }
	  }
	  if ((ii - first_sample_idx) * 100LL >= (int64_t)g_pct * (end_sample_idx - first_sample_idx)) {
	    g_pct = ((ii - first_sample_idx) * 100LL) / (end_sample_idx - first_sample_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile2_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(1, outname, outname_end);
	g_pct = 1;
      }
      if (write_1mibs_matrix) {
	dist_ptr = dists;
	fyy = 0.0;
	for (ii = first_sample_idx; ii < end_sample_idx; ii++) {
	  for (jj = 0; jj < ii; jj++) {
	    fxx = (float)((*dist_ptr++) * half_marker_ct_recip);
	    if (fwrite_checked(&fxx, sizeof(float), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  }
	  if (shape == DISTANCE_SQ0) {
	    if (fwrite_checked(membuf, (sample_ct - ii) * sizeof(float), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	  } else {
	    // square matrix
	    if (fwrite_checked(&fyy, sizeof(float), *outfile3_ptr)) {
	      goto distance_d_write_ret_WRITE_FAIL;
	    }
	    for (ulii = ii + 1; ulii < sample_ct; ulii++) {
	      fxx = (float)(dists[(ulii * (ulii - 1)) / 2 + ii] * half_marker_ct_recip);
	      fwrite(&fxx, 4, 1, *outfile3_ptr);
	    }
	  }
	  if ((ii - first_sample_idx) * 100LL >= (int64_t)g_pct * (end_sample_idx - first_sample_idx)) {
	    g_pct = ((ii - first_sample_idx) * 100LL) / (end_sample_idx - first_sample_idx);
	    printf("\rWriting... %d%%", g_pct++);
	    fflush(stdout);
	  }
	}
	if (fclose_null(outfile3_ptr)) {
	  goto distance_d_write_ret_WRITE_FAIL;
	}
	distance_print_done(2, outname, outname_end);
      }
    }
  } else {
    if (shape == DISTANCE_SQ0) {
      // assume little-endian
#ifdef __LP64__
      ulii = 0x3009300930093009LLU;
#else
      ulii = 0x30093009;
#endif
      ii = (sample_ct * 2 + sizeof(intptr_t) - 2) / sizeof(intptr_t);
      glptr = (uintptr_t*)membuf;
      for (jj = 0; jj < ii; jj++) {
	*glptr++ = ulii;
      }
      g_dw_membuf = membuf;
    }
    g_dw_min_sample = first_sample_idx;
    g_dw_max_sample1idx = end_sample_idx;
    g_dw_dists = dists;
    g_dw_half_marker_ct_recip = half_marker_ct_recip;
    if (write_alcts) {
      g_dw_sample1idx = first_sample_idx;
      g_dw_sample2idx = 0;
      g_dw_dist_ptr = dists;
      if (dist_calc_type & DISTANCE_GZ) {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".dist.%u.gz", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".dist.gz");
	}
	if (shape == DISTANCE_SQ) {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_sq0_emitn);
	} else {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_tri_emitn);
	}
      } else {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".dist.%u", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".dist");
	}
	if (shape == DISTANCE_SQ) {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_sq0_emitn);
	} else {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_tri_emitn);
	}
	if (retval) {
	  goto distance_d_write_ret_1;
	}
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("Distances (allele counts) written to %s .\n", outname);
      g_pct = 1;
    }
    if (write_1mibs_matrix) {
      g_dw_sample1idx = first_sample_idx;
      g_dw_sample2idx = 0;
      g_dw_dist_ptr = dists;
      if (dist_calc_type & DISTANCE_GZ) {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mdist.%u.gz", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mdist.gz");
	}
	if (shape == DISTANCE_SQ) {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_1mibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_1mibs_sq0_emitn);
	} else {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_1mibs_tri_emitn);
	}
      } else {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mdist.%u", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mdist");
	}
	if (shape == DISTANCE_SQ) {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_1mibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_1mibs_sq0_emitn);
	} else {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_1mibs_tri_emitn);
	}
	if (retval) {
	  goto distance_d_write_ret_1;
	}
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("Distances (proportions) written to %s .\n", outname);
      g_pct = 1;
    }
    if (write_ibs_matrix) {
      g_dw_sample1idx = first_sample_idx;
      g_dw_sample2idx = 0;
      g_dw_dist_ptr = dists;
      g_dw_start_offset = ((int64_t)first_sample_idx * (first_sample_idx + 1)) / 2;
      g_dw_hundredth = 1 + ((((int64_t)end_sample_idx * (end_sample_idx + 1)) / 2 - g_dw_start_offset) / 100);
      if (dist_calc_type & DISTANCE_GZ) {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mibs.%u.gz", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mibs.gz");
	}
	if (shape == DISTANCE_SQ) {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_ibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_ibs_sq0_emitn);
	} else {
	  parallel_compress(outname, overflow_buf, 0, distance_d_write_ibs_tri_emitn);
	}
      } else {
	if (parallel_tot > 1) {
	  sprintf(outname_end, ".mibs.%u", parallel_idx + 1);
	} else {
	  sprintf(outname_end, ".mibs");
	}
	if (shape == DISTANCE_SQ) {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_ibs_sq_emitn);
	} else if (shape == DISTANCE_SQ0) {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_ibs_sq0_emitn);
	} else {
	  retval = write_uncompressed(outname, overflow_buf, 0, distance_d_write_ibs_tri_emitn);
	}
	if (retval) {
	  goto distance_d_write_ret_1;
	}
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("IBS matrix written to %s .\n", outname);
    }
  }
  while (0) {
  distance_d_write_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  distance_d_write_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 distance_d_write_ret_1:
  return retval;
}

uint32_t calc_genome_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t* sample_exclude = g_cg_sample_exclude;
  uint32_t* sample_missing_unwt = g_sample_missing_unwt;
  uint32_t* missing_dbl_excluded = g_missing_dbl_excluded;
  uint32_t* genome_main = g_genome_main;
  uintptr_t sample_ct = g_sample_ct;
  uint32_t is_rel_check = g_cg_genome_modifier & GENOME_REL_CHECK;
  uint32_t is_unbounded = g_cg_genome_modifier & GENOME_IBD_UNBOUNDED;
  uint32_t is_nudge = g_cg_genome_modifier & GENOME_NUDGE;
  uint32_t filter_pi_hat = g_cg_genome_modifier & GENOME_FILTER_PI_HAT;
  uint32_t output_full = g_cg_genome_modifier & GENOME_OUTPUT_FULL;
  uintptr_t sample1idx = g_cg_sample1idx;
  uintptr_t sample2idx = g_cg_sample2idx;
  uintptr_t sample1uidx = g_cg_sample1uidx;
  uintptr_t sample2uidx = g_cg_sample2uidx;
  uintptr_t tstc = g_cg_tstc;
  uint32_t marker_ct = g_cg_marker_ct;
  uintptr_t max_sample_id_len = g_cg_max_sample_id_len;
  uint32_t max_sample_fid_len = g_cg_max_sample_fid_len;
  uint32_t max_sample_iid_len = g_cg_max_sample_iid_len;
  uintptr_t max_paternal_id_len = g_cg_max_paternal_id_len;
  uintptr_t max_maternal_id_len = g_cg_max_maternal_id_len;
  uintptr_t* pheno_c = g_cg_pheno_c;
  uintptr_t* founder_info = g_cg_founder_info;
  char* sample_ids = g_cg_sample_ids;
  char* paternal_ids = g_cg_paternal_ids;
  char* maternal_ids = g_cg_maternal_ids;
  uint32_t family_id_fixed = g_cg_family_id_fixed;
  uint32_t is_founder_fixed = g_cg_is_founder_fixed;
  uint32_t rel_space_id_fixed = g_cg_rel_space_id_fixed;
  int64_t llfct = g_cg_llfct;
  char* sptr_start = g_cg_sptr_start;
  char* fam1 = g_cg_fam1;
  char* fam2 = g_cg_fam2;
  char* sample1 = g_cg_sample1;
  char* pat1 = g_cg_pat1;
  char* mat1 = g_cg_mat1;
  int64_t cur_line = g_cg_cur_line;
  int64_t tot_lines = g_cg_tot_lines;
  uint32_t pct = g_pct;
  uintptr_t gmcell = g_cg_gmcell;
  uintptr_t mdecell = g_cg_mdecell;
  double e00 = g_cg_e00;
  double e01 = g_cg_e01;
  double e02 = g_cg_e02;
  double e11 = g_cg_e11;
  double e12 = g_cg_e12;
  char* pat2;
  char* mat2;
  char* cptr;
  char* sample2;
  char* sptr_cur_start;
  uint32_t founder_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  int32_t nn;
  int32_t oo;
  double dxx;
  double dyy;
  double dxx1;
  double dxx2;
  if (!sample2idx) {
    // first line, if not 2nd or later part of parallel write
    if (!sample1idx) {
      sptr_cur += sprintf(sptr_cur, g_textbuf, " FID1", " IID1", " FID2", " IID2");
    }
    sample1idx = g_thread_start[0];
    sample2idx = sample1idx + 1;
    g_textbuf[0] = ' ';
  }
  while (sample1idx < tstc) {
    if (sample2idx == sample1idx + 1) {
      cptr = &(sample_ids[sample1uidx * max_sample_id_len]);
      uii = strlen_se(cptr);
      memcpyx(fam1, cptr, uii, '\0');
      sample1 = next_token(cptr);
      pat1 = &(paternal_ids[sample1uidx * max_paternal_id_len]);
      mat1 = &(maternal_ids[sample1uidx * max_maternal_id_len]);
      is_founder_fixed = IS_SET(founder_info, sample1uidx);
      rel_space_id_fixed = g_cg_pri.family_rel_nf_idxs[sample1uidx];
      family_id_fixed = g_cg_pri.family_idxs[sample1uidx];
      founder_ct = g_cg_pri.family_founder_cts[family_id_fixed];
      llfct = ((int64_t)founder_ct * (founder_ct - 1)) - 2 * g_cg_pri.family_rel_space_offsets[family_id_fixed];
      sptr_start = fw_strcpyn(max_sample_fid_len - 1, uii, fam1, &(g_textbuf[1]));
      *sptr_start++ = ' ';
      sptr_start = fw_strcpy(max_sample_iid_len - 1, sample1, sptr_start);
      *sptr_start++ = ' ';
    }
    while (sample2idx < sample_ct) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample2uidx);
      sptr_cur_start = sptr_cur;
      sptr_cur = memcpya(sptr_cur, g_textbuf, sptr_start - g_textbuf);
      cptr = &(sample_ids[sample2uidx * max_sample_id_len]);
      uii = strlen_se(cptr);
      memcpyx(fam2, cptr, uii, '\0');
      sample2 = skip_initial_spaces(&(cptr[uii + 1]));
      pat2 = &(paternal_ids[sample2uidx * max_paternal_id_len]);
      mat2 = &(maternal_ids[sample2uidx * max_maternal_id_len]);
      sptr_cur = fw_strcpyn(max_sample_fid_len - 1, uii, fam2, sptr_cur);
      *sptr_cur++ = ' ';
      sptr_cur = fw_strcpy(max_sample_iid_len - 1, sample2, sptr_cur);
      *sptr_cur++ = ' ';
      if (!strcmp(fam1, fam2)) {
	// FS = sibling
	// HS = half-sibling
	// PO = parent-offspring
	// OT = other
	while (1) {
	  if (!(is_founder_fixed || IS_SET(founder_info, sample2uidx))) {
	    if ((!strcmp(pat1, pat2)) && (!strcmp(mat1, mat2))) {
	      sptr_cur = memcpyl3a(sptr_cur, "FS ");
	      break;
	    } else if ((!strcmp(pat1, pat2)) || (!strcmp(mat1, mat2))) {
	      sptr_cur = memcpyl3a(sptr_cur, "HS ");
	      break;
	    }
	  }
	  if ((!strcmp(pat1, sample2)) || (!strcmp(mat1, sample2)) || (!strcmp(pat2, sample1)) || (!strcmp(mat2, sample1))) {
	    sptr_cur = memcpyl3a(sptr_cur, "PO ");
	    break;
	  }
	  sptr_cur = memcpyl3a(sptr_cur, "OT ");
	  break;
	}

	// insert relationship
	oo = IS_SET(founder_info, sample2uidx);
	if (is_founder_fixed && oo) {
	  dxx = 0.0;
	} else {
	  // rel_space[] is a sequence of beheaded triangles.
	  nn = g_cg_pri.family_rel_nf_idxs[sample2uidx];
	  if (is_founder_fixed || ((nn > (int32_t)rel_space_id_fixed) && (!oo))) {
	    dxx = g_cg_pri.rel_space[rel_space_id_fixed + ((int64_t)nn * (nn - 1) - llfct) / 2];
	  } else {
	    dxx = g_cg_pri.rel_space[nn + ((int64_t)rel_space_id_fixed * (rel_space_id_fixed - 1) - llfct) / 2];
	  }
	}
	sptr_cur = width_force(5, sptr_cur, dtoa_g(dxx, sptr_cur));
      } else if (!is_rel_check) {
	sptr_cur = memcpya(sptr_cur, "UN    NA", 8);
      } else {
        sptr_cur = sptr_cur_start;
	goto calc_genome_emitn_skip_line;
      }
      nn = marker_ct - sample_missing_unwt[sample1idx] - sample_missing_unwt[sample2idx] + missing_dbl_excluded[mdecell];
      oo = nn - genome_main[gmcell] - genome_main[gmcell + 1];
      dxx = ((double)((int32_t)genome_main[gmcell + 1])) / (e00 * nn);
      dyy = (((double)((int32_t)genome_main[gmcell])) - dxx * e01 * nn) / (e11 * nn);
      dxx1 = ((double)oo - nn * (dxx * e02 + dyy * e12)) / ((double)nn);
      if (!is_unbounded) {
	if (dxx > 1) {
	  dxx = 1;
	  dyy = 0;
	  dxx1 = 0;
	} else if (dyy > 1) {
	  dyy = 1;
	  dxx = 0;
	  dxx1 = 0;
	} else if (dxx1 > 1) {
	  dxx1 = 1;
	  dyy = 0;
	  dxx = 0;
	} else if (dxx < 0) {
	  dxx2 = 1.0 / (dyy + dxx1);
	  dyy *= dxx2;
	  dxx1 *= dxx2;
	  dxx = 0;
	}
	if (dyy < 0) {
	  dxx2 = 1.0 / (dxx + dxx1);
	  dxx *= dxx2;
	  dxx1 *= dxx2;
	  dyy = 0;
	}
	if (dxx1 < 0) {
	  dxx2 = 1.0 / (dxx + dyy);
	  dxx *= dxx2;
	  dyy *= dxx2;
	  dxx1 = 0;
	}
      }
      dxx2 = dyy * 0.5 + dxx1; // PI_HAT
      if (filter_pi_hat && ((dxx2 < g_cg_min_pi_hat) || (dxx2 > g_cg_max_pi_hat))) {
	sptr_cur = sptr_cur_start;
        goto calc_genome_emitn_skip_line;
      }
      if (is_nudge && (dxx2 * dxx2 < dxx1)) {
        dxx = (1 - dxx2) * (1 - dxx2);
	dyy = 2 * dxx2 * (1 - dxx2);
	dxx1 = dxx2 * dxx2;
      }
      *sptr_cur++ = ' ';
      sptr_cur = dtoa_f_w7p4x(dxx, ' ', sptr_cur);
      sptr_cur = dtoa_f_w7p4x(dyy, ' ', sptr_cur);
      sptr_cur = dtoa_f_w7p4x(dxx1, ' ', sptr_cur);
      sptr_cur = dtoa_f_w7p4(dxx2, sptr_cur);

      if (pheno_c) {
	uii = IS_SET(g_cg_pheno_nm, sample1uidx);
	ujj = IS_SET(g_cg_pheno_nm, sample2uidx);
	ukk = IS_SET(pheno_c, sample1uidx);
	umm = IS_SET(pheno_c, sample2uidx);
	if (((!uii) || (!ukk)) && ((!ujj) || (!umm))) {
	  memcpy(sptr_cur, "  -1 ", 5);
	} else if (uii && ujj && ukk && umm) {
	  memcpy(sptr_cur, "   1 ", 5);
	} else {
	  memcpy(sptr_cur, "   0 ", 5);
	}
      } else {
	memcpy(sptr_cur, "  NA ", 5);
      }
      sptr_cur += 5;
      dxx = (double)genome_main[gmcell + 4];
      dyy = (double)genome_main[gmcell + 3];
      dxx1 = 1.0 / ((double)(genome_main[gmcell + 4] + genome_main[gmcell + 3]));
      dxx2 = normdist((dxx * dxx1 - 0.666666) / (sqrt(0.2222222 * dxx1)));
      sptr_cur = dtoa_f_w9p6x(1.0 - (genome_main[gmcell] + 2 * genome_main[gmcell + 1]) / ((double)(2 * nn)), ' ', sptr_cur);
      if (dxx2 != dxx2) {
	sptr_cur = memcpya(sptr_cur, "     NA ", 8);
      } else {
	sptr_cur = dtoa_f_w7p4x(dxx2, ' ', sptr_cur);
      }
      if (genome_main[gmcell + 3]) {
	sptr_cur = dtoa_f_w7p4(dxx / dyy, sptr_cur);
      } else {
	sptr_cur = memcpya(sptr_cur, "     NA", 7);
      }
      if (output_full) {
	*sptr_cur++ = ' ';
	sptr_cur = uint32toa_w7x(genome_main[gmcell + 1], ' ', sptr_cur);
	sptr_cur = uint32toa_w7x(genome_main[gmcell], ' ', sptr_cur);
	sptr_cur = uint32toa_w7x(oo, ' ', sptr_cur);
	sptr_cur = dtoa_f_w7p4x(dyy, ' ', sptr_cur);
	sptr_cur = dtoa_f_w7p4(dxx, sptr_cur);
      }
      *sptr_cur++ = '\n';
    calc_genome_emitn_skip_line:
      gmcell += 5;
      mdecell++;
      sample2idx++;
      sample2uidx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_genome_emitn_ret;
      }
    }
    cur_line += sample_ct - sample1idx - 1;
    sample1idx++;
    sample1uidx++;
    next_unset_ul_unsafe_ck(sample_exclude, &sample1uidx);
    sample2idx = sample1idx + 1;
    sample2uidx = sample1uidx + 1;
  }
 calc_genome_emitn_ret:
  if (cur_line * 100 >= tot_lines * pct) {
    pct = (cur_line * 100) / tot_lines;
    printf("\rWriting... %u%%", pct++);
    fflush(stdout);
  }
  g_cg_sample1idx = sample1idx;
  g_cg_sample2idx = sample2idx;
  g_cg_sample1uidx = sample1uidx;
  g_cg_sample2uidx = sample2uidx;
  g_cg_family_id_fixed = family_id_fixed;
  g_cg_is_founder_fixed = is_founder_fixed;
  g_cg_rel_space_id_fixed = rel_space_id_fixed;
  g_cg_llfct = llfct;
  g_cg_sptr_start = sptr_start;
  g_cg_sample1 = sample1;
  g_cg_pat1 = pat1;
  g_cg_mat1 = mat1;
  g_cg_cur_line = cur_line;
  g_cg_tot_lines = tot_lines;
  g_pct = pct;
  g_cg_gmcell = gmcell;
  g_cg_mdecell = mdecell;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t calc_genome(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, double* set_allele_freqs, uint32_t* nchrobs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, int32_t nonfounders, uint64_t calculation_type, uint32_t genome_modifier, uint32_t ppc_gap, double min_pi_hat, double max_pi_hat, uintptr_t* pheno_nm, uintptr_t* pheno_c, Pedigree_rel_info pri, uint32_t skip_write) {
  FILE* outfile = nullptr;
  int32_t retval = 0;
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  unsigned char* loadbuf = nullptr; // from file
  int32_t ibd_prect = 0;
  int64_t cur_line = 0;
  double e00 = 0;
  double e01 = 0;
  double e02 = 0;
  double e11 = 0;
  double e12 = 0;
  uint32_t mp_lead_unfiltered_idx = 0;
  uint32_t mp_lead_idx = 0;
  uint32_t chrom_fo_idx = 0;
  uint32_t dist_thread_ct = g_thread_ct;
  char wbuf[16];
  int32_t missing_ct_buf[BITCT];
  double set_allele_freq_buf[GENOME_MULTIPLEX];
  uint32_t nchrobs_buf[GENOME_MULTIPLEX];
  unsigned char* overflow_buf;
  unsigned char* gptr;
  char* cptr;
  uintptr_t* geno;
  uintptr_t* masks;
  uintptr_t* mmasks;
  uintptr_t* glptr;
  uintptr_t* glptr2;
  uintptr_t* glptr3;
  uint32_t* missing_dbl_excluded;
  uint32_t* sample_missing_unwt;
  uint32_t* genome_main;
  uint32_t* giptr;
  uint32_t* giptr2;
  uint32_t* giptr3;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  uint32_t upp;
  uint32_t uqq;
  int32_t missing_ct_all;
  double dpp;
  double dqq;
  double dpp_sq;
  double dqq_sq;
  double dxx;
  double dyy;
  double dxx_recip;
  double dyy_recip;
  double dxx1;
  double dxx2;
  double dyy1;
  double dyy2;
  int64_t tot_cells;
  double num_alleles;
  double num_allelesf2;
  double num_allelesf3;
  uint32_t chrom_end;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t is_haploid;
  uintptr_t marker_uidx;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uint32_t marker_idx;
  uint32_t pct;

  if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logerrprint("Error: --genome cannot be used on haploid genomes.\n");
    goto calc_genome_ret_INVALID_CMDLINE;
  }
#ifdef __LP64__
  // 32-bit value in main loop can be up to
  //   sample_ct * (GENOME_MULTIPLEX2 / BITCT) - 1
  if (sample_ct > 119304647) {
    logerrprint("Error: --genome does not support > 119304647 samples.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto calc_genome_ret_1;
  }
#endif
  g_sample_ct = sample_ct;
  if (dist_thread_ct > sample_ct / 32) {
    dist_thread_ct = sample_ct / 32;
    if (!dist_thread_ct) {
      dist_thread_ct = 1;
    }
  }
  g_cg_max_sample_fid_len = plink_maxfid;
  g_cg_max_sample_iid_len = plink_maxiid;
  g_cg_min_pi_hat = min_pi_hat;
  g_cg_max_pi_hat = max_pi_hat;

  triangle_fill(sample_ct, dist_thread_ct, parallel_tot - parallel_idx - 1, parallel_tot, 1, 1, g_thread_start);
  // invert order, for --genome --parallel to naturally work
  for (uii = 0; uii <= dist_thread_ct / 2; uii++) {
    ujj = g_thread_start[uii];
    g_thread_start[uii] = sample_ct - g_thread_start[dist_thread_ct - uii];
    g_thread_start[dist_thread_ct - uii] = sample_ct - ujj;
  }

  if (!parallel_idx) {
    cur_line = 1;
  }
  g_cg_tstc = g_thread_start[dist_thread_ct];
  // f(tstc) - f(g_thread_start[0])
  // f(0) = 0
  // f(1) = sample_ct - 1
  // f(2) = 2sample_ct - 3
  // ...
  // f(n) = nsample_ct - n(n+1)/2
  tot_cells = (int64_t)sample_ct * (g_cg_tstc - g_thread_start[0]) - ((int64_t)g_cg_tstc * (g_cg_tstc + 1) - (int64_t)g_thread_start[0] * (g_thread_start[0] + 1)) / 2;
  g_cg_tot_lines = cur_line + tot_cells;
  if (bigstack_alloc_ui(tot_cells, &missing_dbl_excluded) ||
      bigstack_alloc_ui(sample_ct, &sample_missing_unwt) ||
      bigstack_alloc_ui(tot_cells * 5, &genome_main) ||
      bigstack_alloc_uc(GENOME_MULTIPLEX * unfiltered_sample_ct4, &loadbuf) ||
      bigstack_alloc_ul(sample_ct * (GENOME_MULTIPLEX / BITCT2), &geno) ||
      bigstack_alloc_ul(sample_ct * (GENOME_MULTIPLEX / BITCT2), &masks) ||
      bigstack_alloc_ul(sample_ct, &mmasks) ||
      bigstack_alloc_c(plink_maxfid + 1, &g_cg_fam1) ||
      bigstack_alloc_c(plink_maxfid + 1, &g_cg_fam2) ||
      bigstack_alloc_uc(262144, &overflow_buf)) {
    goto calc_genome_ret_NOMEM;
  }

  g_geno = (unsigned char*)geno;
  g_masks = masks;
  g_mmasks = mmasks;
  g_missing_dbl_excluded = missing_dbl_excluded;
  g_sample_missing_unwt = sample_missing_unwt;
  g_genome_main = genome_main;
  fill_uint_zero(tot_cells, missing_dbl_excluded);
  fill_uint_zero(sample_ct, sample_missing_unwt);
  fill_uint_zero(tot_cells * 5, genome_main);
  if (!IS_SET(marker_exclude, 0)) {
    if (fseeko(bedfile, bed_offset, SEEK_SET)) {
      retval = RET_READ_FAIL;
      goto calc_genome_ret_1;
    }
  }
  marker_uidx = 0; // raw marker index
  marker_idx = 0; // after excluding missing
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
  // subtract X/MT/haploid markers from marker_ct
  ukk = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
  if (ukk) {
    LOGPRINTF("Excluding %u variant%s on non-autosomes from IBD calculation.\n", ukk, (ukk == 1)? "" : "s");
    marker_ct -= ukk;
  }
  do {
    ukk = marker_ct - marker_idx;
    if (ukk > GENOME_MULTIPLEX) {
      ukk = GENOME_MULTIPLEX;
    }
    glptr2 = g_marker_window;
    for (ujj = 0; ujj < ukk; ujj++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  retval = RET_READ_FAIL;
	  goto calc_genome_ret_1;
	}
      }
      if (marker_uidx >= chrom_end) {
	while (1) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_mt, &is_haploid);
	  if ((!is_haploid) && (!is_mt)) {
	    break;
	  }
          marker_uidx = next_unset_ul_unsafe(marker_exclude, chrom_end);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    retval = RET_READ_FAIL;
	    goto calc_genome_ret_1;
	  }
	}
      }
      if (fread(&(loadbuf[ujj * unfiltered_sample_ct4]), 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4) {
	retval = RET_READ_FAIL;
	goto calc_genome_ret_1;
      }
      set_allele_freq_buf[ujj] = set_allele_freqs[marker_uidx];
      if (nchrobs) {
	nchrobs_buf[ujj] = nchrobs[marker_uidx];
      }
      // See comments in incr_genome(): the PPC test is time-critical and
      // we do a bit of unusual precomputation here to speed it up.
      //
      // Objective: Fill glptr[0] and glptr[1] with either
      // * a bitmask that excludes the correct number of markers, if the next
      //   eligible marker for the PPC test is within the same (BITCT / 2)
      //   marker window, or
      // * twice the offset of the next marker eligible for the PPC test,
      //   relative to the bottom of the currently loaded window,
      // because distinguishing between these two cases is effectively free.
      //
      // Then advance glptr two spaces.  The double storage eliminates a
      // divide-by-two in the inner loop at a low cost in cache space.
      if (mp_lead_unfiltered_idx < marker_uidx) {
	mp_lead_unfiltered_idx = marker_uidx;
	mp_lead_idx = marker_idx + ujj;
      }
      if (mp_lead_unfiltered_idx < chrom_end) {
	uii = ppc_gap + marker_pos[marker_uidx];
	if (marker_pos[mp_lead_unfiltered_idx] <= uii) {
	  do {
	    if (!IS_SET(marker_exclude, mp_lead_unfiltered_idx)) {
	      mp_lead_idx++;
	    }
	    mp_lead_unfiltered_idx++;
	  } while ((mp_lead_unfiltered_idx < chrom_end) && (IS_SET(marker_exclude, mp_lead_unfiltered_idx) || (marker_pos[mp_lead_unfiltered_idx] <= uii)));
	}
      }
      if (mp_lead_unfiltered_idx < unfiltered_marker_ct) {
	ulii = 2 * (mp_lead_idx - marker_idx);
	if (ulii < BITCT + (2 * (ujj & (~(BITCT2 - 1))))) {
	  ulii = (~ZEROLU) << (ulii & (BITCT - 1));
	}
      } else {
	ulii = 2 * (unfiltered_marker_ct + GENOME_MULTIPLEX);
      }

      *glptr2++ = ulii;
      *glptr2++ = ulii;
      marker_uidx++;
    }
    if (ukk < GENOME_MULTIPLEX) {
      memset(&(loadbuf[ukk * unfiltered_sample_ct4]), 0, (GENOME_MULTIPLEX - ukk) * unfiltered_sample_ct4);
      fill_ulong_zero(sample_ct * (GENOME_MULTIPLEX / BITCT2), geno);
      fill_ulong_zero(sample_ct * (GENOME_MULTIPLEX / BITCT2), masks);
      for (umm = ukk * 2; umm < GENOME_MULTIPLEX2; umm++) {
	*glptr2++ = GENOME_MULTIPLEX2;
      }
    }
    g_case_ct = marker_idx + ukk;
    is_last_block = (g_case_ct == marker_ct);
    for (ujj = 0; ujj < ukk; ujj += BITCT) {
      glptr = &(geno[ujj / BITCT2]);
      glptr2 = &(masks[ujj / BITCT2]);
      glptr3 = mmasks;
      giptr = sample_missing_unwt;
      sample_uidx = 0;
      fill_int_zero(BITCT, missing_ct_buf);
      missing_ct_all = 0;
      for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	// er, switch this to new loop structure
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	uoo = (sample_uidx % 4) * 2;
	ulii = 0;
	ulkk = 0;
        gptr = &(loadbuf[sample_uidx / 4 + ujj * unfiltered_sample_ct4]);
	uqq = (nonfounders || IS_SET(founder_info, sample_uidx));
	if (uqq) {
	  for (upp = 0; upp < BITCT2; upp++) {
	    uljj = (gptr[upp * unfiltered_sample_ct4] >> uoo) & 3;
	    ulii |= uljj << (upp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[upp] += 1;
	      ulkk |= ONELU << upp;
	      *giptr += 1;
	    }
	  }
	} else {
	  missing_ct_all++;
	  for (upp = 0; upp < BITCT2; upp++) {
	    uljj = (gptr[upp * unfiltered_sample_ct4] >> uoo) & 3;
	    ulii |= uljj << (upp * 2);
	    if (uljj == 1) {
	      ulkk |= ONELU << upp;
	      *giptr += 1;
	    }
	  }
	}
	ulii ^= FIVEMASK;
	*glptr++ = ulii;
	ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	*glptr2++ = ulii * 3;
	*glptr3 = ulkk;

	ulii = 0;
	ulkk = 0;
	gptr = &(loadbuf[sample_uidx / 4 + (ujj + BITCT2) * unfiltered_sample_ct4]);
	if (uqq) {
	  for (upp = 0; upp < BITCT2; upp++) {
	    uljj = (gptr[upp * unfiltered_sample_ct4] >> uoo) & 3;
	    ulii |= uljj << (upp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[upp + BITCT2] += 1;
	      ulkk |= ONELU << upp;
	      *giptr += 1;
	    }
	  }
	} else {
	  for (upp = 0; upp < BITCT2; upp++) {
	    uljj = (gptr[upp * unfiltered_sample_ct4] >> uoo) & 3;
	    ulii |= uljj << (upp * 2);
	    if (uljj == 1) {
	      ulkk |= ONELU << upp;
	      *giptr += 1;
	    }
	  }
	}
	ulii ^= FIVEMASK;
	*glptr = ulii;
	ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	*glptr2 = ulii * 3;
	*glptr3++ |= ulkk << BITCT2;
	glptr = &(glptr[(GENOME_MULTIPLEX2 / BITCT) - 1]);
	glptr2 = &(glptr2[(GENOME_MULTIPLEX2 / BITCT) - 1]);
	giptr++;
      }
      unn = ukk - ujj;
      if (unn > BITCT) {
	unn = BITCT;
      }
      for (umm = 0; umm < unn; umm++) {
	dpp = set_allele_freq_buf[ujj + umm];
	dqq = 1.0 - dpp;
	if ((!nchrobs) || (nchrobs_buf[ujj + umm] == 0xffffffffU)) {
	  num_alleles = (double)(2 * (sample_ct - missing_ct_buf[umm] - missing_ct_all));
	} else {
	  num_alleles = (double)nchrobs_buf[ujj + umm];
	}
	if ((num_alleles > 3) && (dpp > 0.0) && (dqq > 0.0)) {
	  // update e00, e01, e02, e11, e12, ibd_prect
	  // see Plink::preCalcGenomeIBD() in genome.cpp
	  num_allelesf2 = num_alleles * num_alleles / ((num_alleles - 1) * (num_alleles - 2));
	  num_allelesf3 = num_allelesf2 * num_alleles / (num_alleles - 3);
	  dxx = dpp * num_alleles;
	  dyy = dqq * num_alleles;
	  dxx_recip = 1.0 / dxx;
	  dyy_recip = 1.0 / dyy;
	  dpp_sq = dpp * dpp;
	  dqq_sq = dqq * dqq;
	  dxx1 = (dxx - 1) * dxx_recip;
	  dxx2 = dxx1 * (dxx - 2) * dxx_recip;
	  dyy1 = (dyy - 1) * dyy_recip;
	  dyy2 = dyy1 * (dyy - 2) * dyy_recip;
	  e00 += 2 * dpp_sq * dqq_sq * dxx1 * dyy1 * num_allelesf3;
	  e01 += 4 * dpp * dqq * num_allelesf3 * (dpp_sq * dxx2 + dqq_sq * dyy2);
	  e02 += num_allelesf3 * (dqq_sq * dqq_sq * dyy2 * (dyy - 3) * dyy_recip + dpp_sq * dpp_sq * dxx2 * (dxx - 3) * dxx_recip + 4 * dpp_sq * dqq_sq * dxx1 * dyy1);
	  e11 += 2 * dpp * dqq * num_allelesf2 * (dpp * dxx1 + dqq * dyy1);
	  e12 += num_allelesf2 * (dpp_sq * dpp * dxx2 + dqq_sq * dqq * dyy2 + dpp_sq * dqq * dxx1 + dpp * dqq_sq * dyy1);
	  ibd_prect++;
	}
      }
      uii = is_last_block && (ujj + BITCT >= ukk);
      if (spawn_threads2(threads, &calc_genome_thread, dist_thread_ct, uii)) {
	goto calc_genome_ret_THREAD_CREATE_FAIL;
      }
      ulii = 0;
      calc_genome_thread((void*)ulii);
      join_threads2(threads, dist_thread_ct, uii);
    }

    marker_idx = g_case_ct;
    g_ctrl_ct = g_case_ct;
    printf("\r%u markers complete.", marker_idx);
    fflush(stdout);
  } while (!is_last_block);
  fputs("\rIBD calculations complete.  \n", stdout);
  logstr("IBD calculations complete.\n");
  if (skip_write) {
    goto calc_genome_ret_1;
  }
  dxx = 1.0 / ((double)ibd_prect);
  g_cg_e00 = e00 * dxx;
  g_cg_e01 = e01 * dxx;
  g_cg_e02 = e02 * dxx;
  g_cg_e11 = e11 * dxx;
  g_cg_e12 = e12 * dxx;

  if (calculation_type & CALC_PLINK1_IBS_MATRIX) {
    strcpy(outname_end, ".mibs");
    if (fopen_checked(outname, "w", &outfile)) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    giptr = genome_main;
    giptr2 = missing_dbl_excluded;
    pct = 1;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      giptr3 = sample_missing_unwt;
      uii = marker_ct - giptr3[sample_idx];
      uljj = sample_idx - 1; // not referenced when sample_idx == 0
      for (ulii = 0; ulii < sample_idx; ulii++) {
        cptr = dtoa_gx(1.0 - ((double)(genome_main[uljj * 5] + 2 * genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + missing_dbl_excluded[uljj]))), ' ', wbuf);
        fwrite(wbuf, 1, cptr - wbuf, outfile);
	uljj += sample_ct - ulii - 2;
      }
      putc_unlocked('1', outfile);
      putc_unlocked(' ', outfile);
      giptr3++;
      for (ujj = sample_idx + 1; ujj < sample_ct; ujj++) {
	cptr = dtoa_gx(1.0 - ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++)))), ' ', wbuf);
	fwrite(wbuf, 1, cptr - wbuf, outfile);
	giptr = &(giptr[5]);
      }
      if (putc_checked('\n', outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (sample_idx * 100 >= (pct * sample_ct)) {
        pct = (sample_idx * 100) / sample_ct;
	printf("\rWriting... %d%%", pct++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
    *outname_end = '\0';
    LOGPRINTFWW("IBS matrix written to %s.mibs , and IDs written to %s.mibs.id .\n", outname, outname);
  }

  if (calculation_type & CALC_PLINK1_DISTANCE_MATRIX) {
    strcpy(outname_end, ".mdist");
    if (fopen_checked(outname, "w", &outfile)) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    giptr = genome_main;
    giptr2 = missing_dbl_excluded;
    ukk = 1;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      giptr3 = sample_missing_unwt;
      uii = marker_ct - giptr3[sample_idx];
      uljj = sample_idx - 1;
      for (ulii = 0; ulii < sample_idx; ulii++) {
	cptr = dtoa_gx(((double)(genome_main[uljj * 5] + 2 * genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + missing_dbl_excluded[uljj]))), ' ', wbuf);
	fwrite(wbuf, 1, cptr - wbuf, outfile);
	uljj += sample_ct - ulii - 2;
      }
      putc_unlocked('0', outfile);
      putc_unlocked(' ', outfile);
      giptr3++;
      for (ujj = sample_idx + 1; ujj < sample_ct; ujj++) {
	cptr = dtoa_gx(((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++)))), ' ', wbuf);
        fwrite(wbuf, 1, cptr - wbuf, outfile);
	giptr = &(giptr[5]);
      }
      if (putc_checked('\n', outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (sample_idx * 100 >= (ukk * sample_ct)) {
	ukk = (sample_idx * 100) / sample_ct;
	printf("\rWriting... %d%%", ukk++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
    *outname_end = '\0';
    LOGPRINTFWW("Distances (proportions) written to %s.mdist , and IDs written to %s.mdist.id .\n", outname, outname);
  }

  if (!parallel_idx) {
    if (genome_modifier & GENOME_OUTPUT_FULL) {
      sprintf(g_textbuf, "%%%us%%%us%%%us%%%us RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO    IBS0    IBS1    IBS2  HOMHOM  HETHET\n", plink_maxfid, plink_maxiid, plink_maxfid, plink_maxiid);
    } else {
      sprintf(g_textbuf, "%%%us%%%us%%%us%%%us RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO\n", plink_maxfid, plink_maxiid, plink_maxfid, plink_maxiid);
    }
  }
  g_pct = 1;
  g_cg_gmcell = 0;
  g_cg_mdecell = 0;
  g_cg_marker_ct = marker_ct;
  g_cg_max_sample_id_len = max_sample_id_len;
  g_cg_max_paternal_id_len = max_paternal_id_len;
  g_cg_max_maternal_id_len = max_maternal_id_len;
  g_cg_pheno_nm = pheno_nm;
  g_cg_pheno_c = pheno_c;
  g_cg_founder_info = founder_info;
  g_cg_sample_ids = sample_ids;
  g_cg_paternal_ids = paternal_ids;
  g_cg_maternal_ids = maternal_ids;
  g_cg_genome_modifier = genome_modifier;
  g_cg_pri = pri;
  g_cg_cur_line = 0;

  // signal to emitn() that this is the first line, and whether or not this is
  // a later part of a parallel write
  g_cg_sample1idx = parallel_idx;
  g_cg_sample2idx = 0;
  g_cg_sample_exclude = sample_exclude;
  g_cg_sample1uidx = jump_forward_unset_unsafe(sample_exclude, 0, g_thread_start[0] + 1);
  g_cg_sample2uidx = g_cg_sample1uidx + 1;
  if (genome_modifier & GENOME_OUTPUT_GZ) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome.gz");
    }
    parallel_compress(outname, overflow_buf, 0, calc_genome_emitn);
  } else {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome");
    }
    retval = write_uncompressed(outname, overflow_buf, 0, calc_genome_emitn);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }
  putc_unlocked('\r', stdout);
  LOGPRINTFWW("Finished writing %s .\n", outname);
  while (0) {
  calc_genome_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_genome_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_genome_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_genome_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  calc_genome_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_genome_ret_1:
  fclose_cond(outfile);
  if ((!retval) && (calculation_type & (CALC_CLUSTER | CALC_NEIGHBOR))) {
    bigstack_reset(loadbuf);
  } else {
    bigstack_reset(bigstack_mark);
  }
  return retval;
}

static inline void rel_cut_arr_dec(int32_t* rel_ct_arr_elem, uint32_t* exactly_one_rel_ct_ptr) {
  int32_t rcae = *rel_ct_arr_elem - 1;
  *rel_ct_arr_elem = rcae;
  if (rcae < 2) {
    if (rcae) {
      *exactly_one_rel_ct_ptr += 1;
    } else {
      *exactly_one_rel_ct_ptr -= 1;
    }
  }
}

int32_t do_rel_cutoff(uint64_t calculation_type, double rel_cutoff, double* rel_ibc, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* outname, char* outname_end, uintptr_t unfiltered_sample_ct, char* sample_ids, uintptr_t max_sample_id_len) {
  uint32_t samples_excluded = 0;
  uint32_t exactly_one_rel_ct = 0;
  uintptr_t sample_ct = unfiltered_sample_ct - (*sample_exclude_ct_ptr);
  unsigned char* bigstack_mark = g_bigstack_base;
  double* rel_dists = g_rel_dists;
  double* dist_ptr = rel_dists;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  uint32_t* giptr;
  uint32_t* giptr2;
  // number of too-close relations, -1 if excluded
  int32_t* rel_ct_arr;
  uintptr_t sample_idx;
  uint32_t uii;
  uintptr_t ulii;
  int32_t kk;
  int32_t mm;
  int32_t retval;

  // Algorithm:
  // - Whenever there is at least one sample with exactly one remaining
  //   too-close relation, prune the other side of that relationship, because
  //   doing so is never suboptimal.
  // - Otherwise, there's no efficient rule that is always optimal (assuming P
  //   != NP, anyway), so we use a simple heuristic: prune the first sample
  //   with the largest number of remaining too-close relationships.

  if (bigstack_calloc_i(sample_ct, &rel_ct_arr)) {
    return RET_NOMEM;
  }
  for (sample_idx = 1; sample_idx < sample_ct; sample_idx++) {
    for (uii = 0; uii < sample_idx; uii++) {
      if (*dist_ptr++ > rel_cutoff) {
	rel_ct_arr[sample_idx] += 1;
	rel_ct_arr[uii] += 1;
      }
    }
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    if (rel_ct_arr[sample_idx] == 1) {
      exactly_one_rel_ct++;
    }
  }
  while (1) {
    kk = 0;
    if (exactly_one_rel_ct) {
      // there is at least one sample with exactly one too-close relation left,
      // find the first one
      while (rel_ct_arr[kk] != 1) {
	kk++;
      }
      // and now find the identity of the other side
      dist_ptr = &(rel_dists[((intptr_t)kk * (kk - 1)) / 2]);
      for (mm = 0; mm < kk; mm++) {
	if (*dist_ptr > rel_cutoff) {
	  *dist_ptr = 0.0;
	  break;
	}
	dist_ptr++;
      }
      if (mm == kk) {
	do {
	  mm++;
	  dist_ptr = &(rel_dists[tri_coord_no_diag((uint32_t)kk, (uint32_t)mm)]);
	} while (*dist_ptr <= rel_cutoff);
	*dist_ptr = 0.0;
      }
      rel_ct_arr[kk] = 0;
      exactly_one_rel_ct--;
      if (rel_ct_arr[mm] == 1) {
        // speed up the easy case
	exactly_one_rel_ct--;
	rel_ct_arr[mm] = -1;
	samples_excluded++;
	continue;
      }
    } else {
      // find identity of first sample with maximum number of remaining
      // too-close relations
      // kk is highest too-close pair count so far
      mm = -1; // associated sample index
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	if (rel_ct_arr[sample_idx] > kk) {
	  kk = rel_ct_arr[sample_idx];
	  mm = sample_idx;
	}
      }
      // no too-close relations left at all, we're done
      if (mm == -1) {
	break;
      }
    }
    dist_ptr = &(rel_dists[((intptr_t)mm * (mm - 1)) / 2]);
    for (kk = 0; kk < mm; kk++) {
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[kk]), &exactly_one_rel_ct);
      }
      dist_ptr++;
    }
    for (ulii = mm + 1; ulii < sample_ct; ulii++) {
      dist_ptr = &(rel_dists[tri_coord_no_diag((uint32_t)mm, ulii)]);
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[ulii]), &exactly_one_rel_ct);
      }
    }
    rel_ct_arr[mm] = -1;
    samples_excluded++;
  }
  exclude_multi(sample_exclude, rel_ct_arr, unfiltered_sample_ct, sample_exclude_ct_ptr);
  if (samples_excluded) {
    dist_ptr = rel_dists; // write
    dptr2 = rel_dists; // read
    dptr3 = rel_ibc; // write
    dptr4 = rel_ibc; // read
    giptr = g_missing_dbl_excluded; // write
    giptr2 = g_missing_dbl_excluded; // read
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      if (rel_ct_arr[sample_idx] != -1) {
	if (calculation_type & CALC_IBC) {
	  dptr3[sample_ct] = dptr4[sample_ct];
	  dptr3[sample_ct * 2] = dptr4[sample_ct * 2];
	}
	*dptr3 = *dptr4++;
	dptr3++;
	for (uii = 0; uii < sample_idx; uii++) {
	  if (rel_ct_arr[uii] != -1) {
	    *dist_ptr = *dptr2++;
	    dist_ptr++;
	    *giptr = *giptr2++;
	    giptr++;
	  } else {
	    dptr2++;
	    giptr2++;
	  }
	}
      } else {
	dptr4++;
	dptr2 = &(dptr2[sample_idx]);
	giptr2 = &(giptr2[sample_idx]);
      }
    }
    sample_ct -= samples_excluded;
    if (calculation_type & CALC_IBC) {
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	*dptr3++ = *dptr4++;
      }
      dptr4 = &(dptr4[samples_excluded]);
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	*dptr3++ = *dptr4++;
      }
    }
    giptr = g_sample_missing_unwt;
    giptr2 = g_sample_missing_unwt;
    for (sample_idx = 0; sample_idx < sample_ct + samples_excluded; sample_idx++) {
      if (rel_ct_arr[sample_idx] != -1) {
	*giptr = *giptr2++;
	giptr++;
      } else {
	giptr2++;
      }
    }
  }
  LOGPRINTF("%u %s excluded by --rel-cutoff.\n", samples_excluded, species_str(samples_excluded));
  if (!(calculation_type & (CALC_RELATIONSHIP | CALC_GDISTANCE_MASK))) {
    strcpy(outname_end, ".rel.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      return retval;
    }
    LOGPRINTFWW("Remaining sample IDs written to %s .\n", outname);
  }
  bigstack_reset(bigstack_mark);
  return 0;
}

uint32_t g_rcb_row;
uint32_t g_rcb_col;
uint32_t g_rcb_new_row;
uint32_t g_rcb_new_col;
uint32_t g_rcb_sample_ct;
uint64_t g_rcb_progress;
uint64_t g_rcb_hundredth;
FILE* g_rcb_in_binfile;
FILE* g_rcb_in_bin_nfile;
gzFile g_rcb_cur_gzfile;
int32_t* g_rcb_rel_ct_arr;

uint32_t rel_cutoff_batch_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  gzFile cur_gzfile = g_rcb_cur_gzfile;
  uint32_t sample_ct = g_rcb_sample_ct;
  uint32_t row = g_rcb_row;
  uint32_t col = g_rcb_col;
  uint32_t new_row = g_rcb_new_row;
  uint32_t new_col = g_rcb_new_col;
  uint64_t progress = g_rcb_progress;
  uint64_t hundredth = g_rcb_hundredth;
  uint32_t pct = g_pct;
  int32_t* rel_ct_arr = g_rcb_rel_ct_arr;
  char wbuf[16];
  char* cptr;
  uint32_t wbuf_ct;
  uint32_t uii;
  while (row < sample_ct) {
    if (rel_ct_arr[row] == -1) {
      for (uii = 0; uii <= row; uii++) {
	gzgets(cur_gzfile, g_textbuf, MAXLINELEN);
      }
    } else {
      cptr = uint32toa_x(new_row, '\t', wbuf);
      wbuf_ct = (uintptr_t)(cptr - wbuf);
      while (col <= row) {
        gzgets(cur_gzfile, g_textbuf, MAXLINELEN);
	if (rel_ct_arr[col++] != -1) {
	  cptr = next_token_mult(g_textbuf, 2);
          uii = strlen(cptr);
	  sptr_cur = memcpya(sptr_cur, wbuf, wbuf_ct);
	  sptr_cur = uint32toa_x(++new_col, '\t', sptr_cur);
          sptr_cur = memcpya(sptr_cur, cptr, uii);
	  if (sptr_cur >= readbuf_end) {
	    goto rel_cutoff_batch_emitn_ret;
	  }
	}
      }
      new_row++;
      new_col = 0;
    }
    row++;
    progress += row;
    if (progress >= pct * hundredth) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      pct = 1 + (progress / hundredth);
      printf("\b\b%u%%", pct - 1);
      fflush(stdout);
    }
    col = 0;
  }
 rel_cutoff_batch_emitn_ret:
  g_rcb_row = row;
  g_rcb_col = col;
  g_rcb_new_row = new_row;
  g_rcb_new_col = new_col;
  g_rcb_progress = progress;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t rel_cutoff_batch_rbin_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  FILE* in_binfile = g_rcb_in_binfile;
  FILE* in_bin_nfile = g_rcb_in_bin_nfile;
  uint32_t sample_ct = g_rcb_sample_ct;
  uint32_t row = g_rcb_row;
  uint32_t col = g_rcb_col;
  uint32_t new_row = g_rcb_new_row;
  uint32_t new_col = g_rcb_new_col;
  uint64_t progress = g_rcb_progress;
  uint64_t hundredth = g_rcb_hundredth;
  uint32_t pct = g_pct;
  int32_t* rel_ct_arr = g_rcb_rel_ct_arr;
  char wbuf[16];
  char* cptr;
  uint32_t wbuf_ct;
  uint32_t uii;
  float fxx;
  float fyy;
  while (row < sample_ct) {
    if (rel_ct_arr[row] == -1) {
      fseeko(in_binfile, (row + 1) * sizeof(float), SEEK_CUR);
      fseeko(in_bin_nfile, (row + 1) * sizeof(float), SEEK_CUR);
    } else {
      cptr = uint32toa_x(new_row, '\t', wbuf);
      wbuf_ct = (uintptr_t)(cptr - wbuf);
      while (col <= row) {
	if (rel_ct_arr[col] == -1) {
	  uii = col;
	  while (++col <= row) {
	    if (rel_ct_arr[col] != -1) {
	      break;
	    }
	  }
	  fseeko(in_binfile, (col - uii) * sizeof(float), SEEK_CUR);
	  fseeko(in_bin_nfile, (col - uii) * sizeof(float), SEEK_CUR);
	  if (col > row) {
	    break;
	  }
	}
	sptr_cur = memcpya(sptr_cur, wbuf, wbuf_ct);
	sptr_cur = uint32toa_x(++new_col, '\t', sptr_cur);
	if ((fread(&fxx, sizeof(float), 1, in_bin_nfile) != sizeof(float)) || (fread(&fyy, sizeof(float), 1, in_binfile) != sizeof(float))) {
	  // can't use return code here
	  putc_unlocked('\n', stdout);
	  fflush(stdout);
	  fputs("Error: File read failure.\n", stderr);
	  exit(RET_READ_FAIL);
	}
	sptr_cur = uint32toa_x((int32_t)fxx, '\t', sptr_cur);
	sptr_cur = ftoa_ex(fyy, '\n', sptr_cur);
	col++;
	if (sptr_cur >= readbuf_end) {
	  goto rel_cutoff_batch_rbin_emitn_ret;
	}
      }
      new_row++;
      new_col = 0;
    }
    row++;
    progress += row;
    if (progress >= pct * hundredth) {
      if (pct > 10) {
	putc_unlocked('\b', stdout);
      }
      pct = 1 + (progress / hundredth);
      printf("\b\b%u%%", pct - 1);
      fflush(stdout);
    }
    col = 0;
  }
 rel_cutoff_batch_rbin_emitn_ret:
  g_rcb_row = row;
  g_rcb_col = col;
  g_rcb_new_row = new_row;
  g_rcb_new_col = new_col;
  g_rcb_progress = progress;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t rel_cutoff_batch(uint32_t load_grm_bin, char* grmname, char* outname, char* outname_end, Rel_info* relip) {
  // Specialized --rel-cutoff usable on larger files.
  char* grmname_end = (char*)memchr(grmname, 0, FNAMESIZE);
  uintptr_t sample_ct = 0;
  uintptr_t line_idx = 0;
  double rel_cutoff = relip->cutoff;
  FILE* idfile = nullptr;
  FILE* outfile = nullptr;
  FILE* out_bin_nfile = nullptr;
  FILE* in_binfile = nullptr;
  FILE* in_bin_nfile = nullptr;
  gzFile cur_gzfile = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  uint32_t samples_excluded = 0;
  uint32_t exactly_one_rel_ct = 0;
  uint32_t rel_calc_type = relip->modifier & REL_CALC_MASK;
  uintptr_t* compact_rel_table;
  uintptr_t* rtptr;
  unsigned char* overflow_buf;
  char* bufptr;
  uint64_t ullii;
  uint64_t ulljj;
  uint64_t progress;
  uint64_t hundredth;
  uintptr_t tot_words;
  uintptr_t words_left;
  uintptr_t wl_floor;
  uintptr_t cur_word;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint32_t inword_idx;
  uint32_t inword_bound;
  uint32_t uii;
  uint32_t row;
  uint32_t col;
  uint32_t new_row;
  uint32_t pct;
  uintptr_t sample_idx;
  int32_t* rel_ct_arr;
  float rel_cutoff_f;
  float fxx;
  float fyy;
  double dxx;
  int32_t retval;
  int32_t kk;
  int32_t mm;
  int32_t cur_prune;
  ulii = (uintptr_t)(grmname_end - grmname);
  if ((ulii == (uintptr_t)(outname_end - outname)) && (!memcmp(grmname, outname, ulii))) {
    LOGERRPRINTF("Error: --rel-cutoff input and output ID filenames cannot match.%s\n", strcmp(outname, PROG_NAME_STR)? "" : "  (Use --out.)");
    goto rel_cutoff_batch_ret_INVALID_CMDLINE;
  }
  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(grmname, "r", &idfile)) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, idfile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, grmname);
      goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
    }
    if (is_eoln_kns(*(skip_initial_spaces(g_textbuf)))) {
      continue;
    }
    sample_ct++;
  }
  if (!feof(idfile)) {
    goto rel_cutoff_batch_ret_READ_FAIL;
  }
  fclose_null(&idfile);
  ullii = sample_ct;
  ullii = BITCT_TO_WORDCT((ullii * (ullii - 1)) / 2);
#ifndef __LP64__
  if (ullii >= 0x20000000) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
#endif
  tot_words = ullii;
  if (bigstack_calloc_ul(tot_words, &compact_rel_table) ||
      bigstack_calloc_i(sample_ct, &rel_ct_arr) ||
      bigstack_alloc_uc(262144, &overflow_buf)) {
    goto rel_cutoff_batch_ret_NOMEM;
  }

  fputs("Reading... 0%", stdout);
  fflush(stdout);
  words_left = tot_words;
  rtptr = compact_rel_table;
  row = 0;
  col = 0;
  if (load_grm_bin) {
    memcpy(grmname_end, ".grm.bin", 9);
    if (fopen_checked(grmname, FOPEN_RB, &in_binfile)) {
      goto rel_cutoff_batch_ret_OPEN_FAIL;
    }
    rel_cutoff_f = (float)rel_cutoff;
    for (pct = 1; pct <= 100; pct++) {
      wl_floor = (((uint64_t)tot_words) * (100 - pct)) / 100;
      while (words_left > wl_floor) {
	cur_word = 0;
	if (--words_left) {
	  inword_bound = BITCT;
	} else {
	  // only sample_ct mod (BITCT * 2) matters for remainder
	  uii = sample_ct & (BITCT * 2 - 1);
	  inword_bound = ((uii * (uii - 1)) / 2) & (BITCT - 1);
	  if (!inword_bound) {
	    inword_bound = BITCT;
	  }
	}
	for (inword_idx = 0; inword_idx < inword_bound; inword_idx++) {
	  if (!fread(&fxx, 4, 1, in_binfile)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	  if (row == col) {
	    row++;
	    col = 0;
	    if (!fread(&fxx, 4, 1, in_binfile)) {
	      goto rel_cutoff_batch_ret_READ_FAIL;
	    }
	  }
	  if (fxx > rel_cutoff_f) {
	    rel_ct_arr[row] += 1;
	    rel_ct_arr[col] += 1;
	    cur_word |= (ONELU << inword_idx);
	  }
	  col++;
	}
	*rtptr++ = cur_word;
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    fclose_null(&in_binfile);
  } else {
    memcpy(grmname_end, ".grm.gz", 8);
    retval = gzopen_read_checked(grmname, &cur_gzfile);
    if (retval) {
      goto rel_cutoff_batch_ret_1;
    }
    for (pct = 1; pct <= 100; pct++) {
      wl_floor = (((uint64_t)tot_words) * (100 - pct)) / 100;
      while (words_left > wl_floor) {
	cur_word = 0;
	if (--words_left) {
	  inword_bound = BITCT;
	} else {
	  // only sample_ct mod (BITCT * 2) matters for remainder
	  uii = sample_ct & (BITCT * 2 - 1);
	  inword_bound = ((uii * (uii - 1)) / 2) & (BITCT - 1);
	  if (!inword_bound) {
	    inword_bound = BITCT;
	  }
	}
	for (inword_idx = 0; inword_idx < inword_bound; inword_idx++) {
	  if (!gzgets(cur_gzfile, g_textbuf, MAXLINELEN)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	  if (row == col) {
	    row++;
	    col = 0;
	    if (!gzgets(cur_gzfile, g_textbuf, MAXLINELEN)) {
	      goto rel_cutoff_batch_ret_READ_FAIL;
	    }
	  }
	  bufptr = next_token_mult(g_textbuf, 3);
	  if (no_more_tokens_kns(bufptr)) {
	    goto rel_cutoff_batch_ret_INVALID_FORMAT_GENERIC;
	  }
	  if (scan_double(bufptr, &dxx)) {
	    goto rel_cutoff_batch_ret_INVALID_FORMAT_GENERIC;
	  }
	  if (dxx > rel_cutoff) {
	    rel_ct_arr[row] += 1;
	    rel_ct_arr[col] += 1;
	    cur_word |= (ONELU << inword_idx);
	  }
	  col++;
	}
	*rtptr++ = cur_word;
      }
      if (pct < 100) {
	if (pct > 10) {
	  putc_unlocked('\b', stdout);
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (!gzgets(cur_gzfile, g_textbuf, MAXLINELEN)) {
      goto rel_cutoff_batch_ret_READ_FAIL;
    }
    if (gzgets(cur_gzfile, g_textbuf, MAXLINELEN)) {
      goto rel_cutoff_batch_ret_INVALID_FORMAT_GENERIC;
    }
    gzclose(cur_gzfile);
    cur_gzfile = nullptr;
  }
  putc_unlocked('\r', stdout);
  LOGPRINTFWW("%s read complete.  Pruning.\n", grmname);

  // would prefer to just call do_rel_cutoff(), but unfortunately that
  // interferes with the intended "handle extra-large datasets" mission of this
  // function, which necessitates a compact bit representation of the
  // relationship matrix... fortunately, the algorithm is pretty simple.
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    if (rel_ct_arr[sample_idx] == 1) {
      exactly_one_rel_ct++;
    }
  }

  while (1) {
    kk = 0;
    cur_prune = -1;
    if (exactly_one_rel_ct) {
      while (rel_ct_arr[kk] != 1) {
	kk++;
      }
      ullii = (((int64_t)kk) * (kk - 1)) / 2;
      ulii = ullii / BITCT;
      rtptr = &(compact_rel_table[ulii]);
      inword_idx = ullii & (BITCT - 1);
      ulljj = ullii + kk;
      uljj = ulljj / BITCT;
      inword_bound = ulljj & (BITCT - 1);
      if (uljj == ulii) {
	uljj = (ONELU << inword_bound) - (ONELU << inword_idx);
        ulkk = (*rtptr) & uljj;
	if (ulkk) {
	  *rtptr &= ~uljj;
	  cur_prune = CTZLU(ulkk) - inword_idx;
	}
      } else {
        ulkk = (*rtptr) & (~((ONELU << inword_idx) - ONELU));
	if (ulkk) {
	  *rtptr &= (ONELU << inword_idx) - ONELU;
          cur_prune = CTZLU(ulkk) - inword_idx;
	} else {
	  col = BITCT - inword_idx;
          row = col + (uljj - ulii - 1) * BITCT;
          while (col < row) {
            ulkk = *(++rtptr);
            if (ulkk) {
	      *rtptr = 0;
              cur_prune = CTZLU(ulkk) + col;
	      break;
	    }
	    col += BITCT;
	  }
	  if (cur_prune == -1) {
            ulkk = (*(++rtptr)) & ((ONELU << inword_bound) - ONELU);
            if (ulkk) {
	      *rtptr &= (~((ONELU << inword_bound) - ONELU));
	      cur_prune = CTZLU(ulkk) + col;
	    }
	  }
	}
      }
      if (cur_prune == -1) {
        mm = kk + 1;
        while (1) {
          ullii = ((((int64_t)mm) * (mm - 1)) / 2) + kk;
	  ulii = ullii / BITCT;
	  inword_idx = ullii & (BITCT - 1);
	  if (compact_rel_table[ulii] & (ONELU << inword_idx)) {
	    compact_rel_table[ulii] &= ~(ONELU << inword_idx);
	    cur_prune = mm;
	    break;
	  }
          mm++;
	}
      }
      rel_ct_arr[kk] = 0;
      exactly_one_rel_ct--;
      if (rel_ct_arr[cur_prune] == 1) {
	exactly_one_rel_ct--;
	rel_ct_arr[cur_prune] = -1;
	samples_excluded++;
	continue;
      }
    } else {
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
	if (rel_ct_arr[sample_idx] > kk) {
	  kk = rel_ct_arr[sample_idx];
	  cur_prune = sample_idx;
	}
      }
      if (cur_prune == -1) {
	break;
      }
    }
    // zero out cur_prune row/column, update other array entries
    ullii = (((int64_t)cur_prune) * (cur_prune - 1)) / 2;
    ulii = ullii / BITCT;
    rtptr = &(compact_rel_table[ulii]);
    inword_idx = ullii & (BITCT - 1);
    ulljj = ullii + cur_prune;
    uljj = ulljj / BITCT;
    inword_bound = ulljj & (BITCT - 1);
    if (uljj == ulii) {
      uljj = (ONELU << inword_bound) - (ONELU << inword_idx);
      ulkk = (*rtptr) & uljj;
      if (ulkk) {
	do {
	  uii = CTZLU(ulkk) - inword_idx;
	  rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
	  ulkk &= ulkk - 1;
	} while (ulkk);
	*rtptr &= ~uljj;
      }
    } else {
      ulkk = (*rtptr) & (~((ONELU << inword_idx) - ONELU));
      if (ulkk) {
	do {
	  uii = CTZLU(ulkk) - inword_idx;
	  rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
	  ulkk &= ulkk - 1;
	} while (ulkk);
	*rtptr &= (ONELU << inword_idx) - ONELU;
      }
      col = BITCT - inword_idx;
      row = col + (uljj - ulii - 1) * BITCT;
      while (col < row) {
	ulkk = *(++rtptr);
	if (ulkk) {
	  do {
	    uii = CTZLU(ulkk) + col;
	    rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
	    ulkk &= ulkk - 1;
	  } while (ulkk);
	  *rtptr = 0;
	}
	col += BITCT;
      }
      ulkk = (*(++rtptr)) & ((ONELU << inword_bound) - ONELU);
      if (ulkk) {
	do {
	  uii = CTZLU(ulkk) + col;
	  rel_cut_arr_dec(&(rel_ct_arr[uii]), &exactly_one_rel_ct);
          ulkk &= ulkk - 1;
	} while (ulkk);
	*rtptr &= (~((ONELU << inword_bound) - ONELU));
      }
    }

    for (uljj = cur_prune + 1; uljj < sample_ct; uljj++) {
      ullii = ((((uint64_t)uljj) * (uljj - 1)) / 2) + cur_prune;
      ulii = ullii / BITCT;
      rtptr = &(compact_rel_table[ulii]);
      inword_idx = ullii & (BITCT - 1);
      if ((*rtptr) & (ONELU << inword_idx)) {
        rel_cut_arr_dec(&(rel_ct_arr[uljj]), &exactly_one_rel_ct);
        *rtptr &= ~(ONELU << inword_idx);
      }
    }
    rel_ct_arr[cur_prune] = -1;
    samples_excluded++;
  }

  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(grmname, "r", &idfile)) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }

  memcpy(outname_end, ".grm.id", 8);
  if (fopen_checked(outname, "w", &outfile)) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }

  for (sample_idx = 0; sample_idx < sample_ct;) {
    if (fgets(g_textbuf, MAXLINELEN, idfile) == nullptr) {
      goto rel_cutoff_batch_ret_READ_FAIL;
    }
    if (is_eoln_kns(*(skip_initial_spaces(g_textbuf)))) {
      continue;
    }
    if (rel_ct_arr[sample_idx] != -1) {
      if (fputs_checked(g_textbuf, outfile)) {
	goto rel_cutoff_batch_ret_WRITE_FAIL;
      }
    }
    sample_idx++;
  }

  fclose_null(&idfile);
  fclose_null(&outfile);

  LOGPRINTF("%u %s excluded by --rel-cutoff.\n", samples_excluded, species_str(samples_excluded));
  LOGPRINTFWW("Remaining sample IDs written to %s .\n", outname);
  if (rel_calc_type & (REL_CALC_GRM | REL_CALC_GRM_BIN)) {
    if (load_grm_bin) {
      memcpy(grmname_end, ".grm.bin", 9);
      if (fopen_checked(grmname, FOPEN_RB, &in_binfile)) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      g_rcb_in_binfile = in_binfile;
      memcpy(grmname_end, ".grm.N.bin", 11);
      if (fopen_checked(grmname, FOPEN_RB, &in_bin_nfile)) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      g_rcb_in_bin_nfile = in_bin_nfile;
    } else {
      memcpy(grmname_end, ".grm.gz", 8);
      retval = gzopen_read_checked(grmname, &cur_gzfile);
      if (retval) {
	goto rel_cutoff_batch_ret_1;
      }
      g_rcb_cur_gzfile = cur_gzfile;
    }
    fputs("Rewriting matrix... 0%", stdout);
    fflush(stdout);
    if (rel_calc_type & REL_CALC_GRM) {
      g_pct = 1;
      g_rcb_row = 0;
      g_rcb_col = 0;
      g_rcb_new_row = 1;
      g_rcb_new_col = 0;
      g_rcb_progress = 0;
      g_rcb_hundredth = 1 + ((((uint64_t)sample_ct) * (sample_ct - 1)) / 200);
      g_rcb_sample_ct = sample_ct;
      g_rcb_rel_ct_arr = rel_ct_arr;
      if (load_grm_bin) {
	if (rel_calc_type & REL_CALC_GZ) {
	  memcpy(outname_end, ".grm.gz", 8);
	  parallel_compress(outname, overflow_buf, 0, rel_cutoff_batch_rbin_emitn);
	} else {
	  memcpy(outname_end, ".grm", 5);
	  retval = write_uncompressed(outname, overflow_buf, 0, rel_cutoff_batch_rbin_emitn);
	  if (retval) {
	    goto rel_cutoff_batch_ret_1;
	  }
	}
      } else {
	if (rel_calc_type & REL_CALC_GZ) {
	  memcpy(outname_end, ".grm.gz", 8);
	  parallel_compress(outname, overflow_buf, 0, rel_cutoff_batch_emitn);
	} else {
	  memcpy(outname_end, ".grm", 5);
	  retval = write_uncompressed(outname, overflow_buf, 0, rel_cutoff_batch_emitn);
	  if (retval) {
	    goto rel_cutoff_batch_ret_1;
	  }
	}
      }
    } else {
      pct = 1;
      row = 0;
      col = 0;
      new_row = 1;
      progress = 0;
      hundredth = 1 + ((((uint64_t)sample_ct) * (sample_ct - 1)) / 200);
      memcpy(outname_end, ".grm.N.bin", 11);
      if (fopen_checked(outname, FOPEN_WB, &out_bin_nfile)) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      memcpy(outname_end, ".grm.bin", 9);
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      while (row < sample_ct) {
	if (rel_ct_arr[row] == -1) {
	  if (load_grm_bin) {
	    fseeko(in_binfile, (row + 1) * sizeof(float), SEEK_CUR);
	    fseeko(in_bin_nfile, (row + 1) * sizeof(float), SEEK_CUR);
	  } else {
	    for (uii = 0; uii <= row; uii++) {
	      gzgets(cur_gzfile, g_textbuf, MAXLINELEN);
	    }
	  }
	} else {
	  if (load_grm_bin) {
	    while (col <= row) {
	      if (rel_ct_arr[col] == -1) {
		uii = col;
		while (++col <= row) {
		  if (rel_ct_arr[col] != -1) {
		    break;
		  }
		}
		fseeko(in_bin_nfile, (col - uii) * sizeof(float), SEEK_CUR);
		fseeko(in_bin_nfile, (col - uii) * sizeof(float), SEEK_CUR);
		if (col > row) {
		  break;
		}
	      }
	      if ((fread(&fxx, sizeof(float), 1, in_bin_nfile) != sizeof(float)) || (fread(&fyy, sizeof(float), 1, in_binfile) != sizeof(float))) {
		goto rel_cutoff_batch_ret_READ_FAIL;
	      }
	      fwrite(&fxx, 4, 1, out_bin_nfile);
	      fwrite(&fyy, 4, 1, outfile);
	      col++;
	    }
	  } else {
	    while (col <= row) {
	      gzgets(cur_gzfile, g_textbuf, MAXLINELEN);
	      if (rel_ct_arr[col++] != -1) {
		bufptr = next_token_mult(g_textbuf, 2);
		if (scan_float(bufptr, &fxx)) {
		  goto rel_cutoff_batch_ret_INVALID_FORMAT_GENERIC;
		}
		bufptr = next_token(bufptr);
		if (scan_float(bufptr, &fyy)) {
		  goto rel_cutoff_batch_ret_INVALID_FORMAT_GENERIC;
		}
		fwrite(&fxx, 4, 1, out_bin_nfile);
		fwrite(&fyy, 4, 1, outfile);
	      }
	    }
	  }
	  new_row++;
	}
	row++;
	progress += row;
	if (progress >= pct * hundredth) {
	  if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
	  pct = 1 + (progress / hundredth);
	  printf("\b\b%u%%", pct - 1);
	  fflush(stdout);
	}
	col = 0;
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTFWW("Pruned relationship matrix written to %s .\n", outname);
  }
  retval = 0;
  while (0) {
  rel_cutoff_batch_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  rel_cutoff_batch_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  rel_cutoff_batch_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  rel_cutoff_batch_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  rel_cutoff_batch_ret_INVALID_FORMAT_GENERIC:
    logprint("\n");
    logerrprint("Error: Improperly formatted .grm.gz file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  rel_cutoff_batch_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  rel_cutoff_batch_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 rel_cutoff_batch_ret_1:
  fclose_cond(idfile);
  fclose_cond(in_binfile);
  fclose_cond(in_bin_nfile);
  fclose_cond(outfile);
  fclose_cond(out_bin_nfile);
  gzclose_cond(cur_gzfile);
  bigstack_reset(bigstack_mark);
  return retval;
}

static uint32_t g_cr_marker_ct;
static uintptr_t g_cr_sample1idx;
static uintptr_t g_cr_sample2idx;
static uintptr_t g_cr_min_sample;
static uintptr_t g_cr_max_sample1idx;
static uint64_t g_cr_start_offset;
static uint64_t g_cr_hundredth;
static uint32_t* g_cr_mdeptr;
static double* g_cr_dist_ptr;
static double* g_cr_ibc_ptr;

uint32_t calc_rel_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uint64_t start_offset = g_cr_start_offset;
  uint64_t hundredth = g_cr_hundredth;
  double* dist_ptr = g_cr_dist_ptr;
  double* ibc_ptr = g_cr_ibc_ptr;
  uintptr_t max_sample1idx = g_cr_max_sample1idx;
  uintptr_t sample1idx = g_cr_sample1idx;
  uintptr_t sample2idx = g_cr_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_tri_emitn_ret;
      }
    }
    sptr_cur = dtoa_gx(*ibc_ptr++, '\n', sptr_cur);
    sample1idx++;
    if ((((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) >= hundredth * pct) {
      pct = (((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) / hundredth;
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 calc_rel_tri_emitn_ret:
  g_cr_sample1idx = sample1idx;
  g_cr_sample2idx = sample2idx;
  g_cr_dist_ptr = dist_ptr;
  g_cr_ibc_ptr = ibc_ptr;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  unsigned char* zbuf = g_geno;
  double* dist_ptr = g_cr_dist_ptr;
  double* ibc_ptr = g_cr_ibc_ptr;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_cr_max_sample1idx;
  uintptr_t min_sample = g_cr_min_sample;
  uintptr_t sample1idx = g_cr_sample1idx;
  uintptr_t sample2idx = g_cr_sample2idx;
  uint32_t pct = g_pct;
  uintptr_t ulii;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq0_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      sptr_cur = dtoa_g(*ibc_ptr++, sptr_cur);
      sample2idx++;
    }
    if (sptr_cur >= readbuf_end) {
      goto calc_rel_sq0_emitn_ret;
    } else {
      ulii = (1 + (uintptr_t)(readbuf_end - sptr_cur)) / 2;
      if (ulii < (sample_ct - sample2idx)) {
	sptr_cur = memcpya(sptr_cur, zbuf, 2 * ulii);
	sample2idx += ulii;
	goto calc_rel_sq0_emitn_ret;
      }
      ulii = 2 * (sample_ct - sample2idx);
      sptr_cur = memcpya(sptr_cur, zbuf, ulii);
    }
    *sptr_cur++ = '\n';
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 calc_rel_sq0_emitn_ret:
  g_cr_sample1idx = sample1idx;
  g_cr_sample2idx = sample2idx;
  g_cr_dist_ptr = dist_ptr;
  g_cr_ibc_ptr = ibc_ptr;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  double* rel_dists = g_rel_dists;
  double* dist_ptr = g_cr_dist_ptr;
  double* ibc_ptr = g_cr_ibc_ptr;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t max_sample1idx = g_cr_max_sample1idx;
  uintptr_t min_sample = g_cr_min_sample;
  uintptr_t sample1idx = g_cr_sample1idx;
  uintptr_t sample2idx = g_cr_sample2idx;
  uint32_t pct = g_pct;
  while (sample1idx < max_sample1idx) {
    while (sample2idx < sample1idx) {
      sptr_cur = dtoa_gx(*dist_ptr++, '\t', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq_emitn_ret;
      }
    }
    if (sample2idx == sample1idx) {
      sptr_cur = dtoa_g(*ibc_ptr++, sptr_cur);
      sample2idx++;
    }
    while (sample2idx < sample_ct) {
      *sptr_cur++ = '\t';
      sptr_cur = dtoa_g(rel_dists[tri_coord_no_diag(sample1idx, sample2idx)], sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    sample1idx++;
    if ((sample1idx - min_sample) * 100LLU >= ((uint64_t)pct) * (max_sample1idx - min_sample)) {
      pct = ((sample1idx - min_sample) * 100LLU) / (max_sample1idx - min_sample);
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 calc_rel_sq_emitn_ret:
  g_cr_sample1idx = sample1idx;
  g_cr_sample2idx = sample2idx;
  g_cr_dist_ptr = dist_ptr;
  g_cr_ibc_ptr = ibc_ptr;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_grm_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uint64_t start_offset = g_cr_start_offset;
  uint64_t hundredth = g_cr_hundredth;
  uint32_t* sample_missing_unwt = g_sample_missing_unwt;
  uint32_t* mdeptr = g_cr_mdeptr;
  double* dist_ptr = g_cr_dist_ptr;
  double* ibc_ptr = g_cr_ibc_ptr;
  uintptr_t max_sample1idx = g_cr_max_sample1idx;
  uintptr_t sample1idx = g_cr_sample1idx;
  uintptr_t sample2idx = g_cr_sample2idx;
  uint32_t marker_ct = g_cr_marker_ct;
  uint32_t pct = g_pct;
  char wbuf[16];
  char* wbuf_end;
  uint32_t wbuf_len;
  uint32_t uii;
  while (sample1idx < max_sample1idx) {
    uii = marker_ct - sample_missing_unwt[sample1idx];
    wbuf_end = uint32toa_x(sample1idx + 1, '\t', wbuf);
    wbuf_len = (uintptr_t)(wbuf_end - wbuf);
    while (sample2idx < sample1idx) {
      sptr_cur = memcpya(sptr_cur, wbuf, wbuf_len);
      sptr_cur = uint32toa_x(sample2idx + 1, '\t', sptr_cur);
      sptr_cur = uint32toa_x((uii - sample_missing_unwt[sample2idx]) + (*mdeptr++), '\t', sptr_cur);
      sptr_cur = dtoa_ex(*dist_ptr++, '\n', sptr_cur);
      sample2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_grm_emitn_ret;
      }
    }
    sptr_cur = memcpya(sptr_cur, wbuf, wbuf_len);
    sptr_cur = uint32toa_x(++sample1idx, '\t', sptr_cur);
    sptr_cur = uint32toa_x(uii, '\t', sptr_cur);
    sptr_cur = dtoa_ex(*ibc_ptr++, '\n', sptr_cur);
    if ((((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) >= hundredth * pct) {
      pct = (((uint64_t)sample1idx) * (sample1idx + 1) / 2 - start_offset) / hundredth;
      printf("\rWriting... %u%%", pct++);
      fflush(stdout);
    }
    sample2idx = 0;
  }
 calc_rel_grm_emitn_ret:
  g_cr_sample1idx = sample1idx;
  g_cr_sample2idx = sample2idx;
  g_cr_dist_ptr = dist_ptr;
  g_cr_ibc_ptr = ibc_ptr;
  g_cr_mdeptr = mdeptr;
  g_pct = pct;
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t block_load(FILE* bedfile, int32_t bed_offset, uintptr_t* marker_exclude, uint32_t marker_ct, uint32_t block_max_size, uintptr_t unfiltered_sample_ct4, unsigned char* readbuf, uintptr_t* marker_uidx_ptr, uintptr_t* marker_idx_ptr, uint32_t* block_size_ptr) {
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uintptr_t marker_idx = *marker_idx_ptr;
  uint32_t markers_read = 0;
  if (block_max_size > marker_ct - marker_idx) {
    block_max_size = marker_ct - marker_idx;
  }
  while (markers_read < block_max_size) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	return RET_READ_FAIL;
      }
    }
    if (fread(&(readbuf[markers_read * unfiltered_sample_ct4]), 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4) {
      return RET_READ_FAIL;
    }
    markers_read++;
    marker_idx++;
    marker_uidx++;
  }

  *marker_uidx_ptr = marker_uidx;
  *marker_idx_ptr = marker_idx;
  *block_size_ptr = markers_read;
  return 0;
}

void copy_set_allele_freqs(uintptr_t marker_uidx, uintptr_t* marker_exclude, uint32_t block_max_size, uintptr_t marker_idx, uint32_t marker_ct, uintptr_t* marker_reverse, double* set_allele_freqs, double* set_allele_freq_buf) {
  uint32_t markers_read = 0;
  if (block_max_size > marker_ct - marker_idx) {
    block_max_size = marker_ct - marker_idx;
  }
  while (markers_read < block_max_size) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    if ((!marker_reverse) || (!IS_SET(marker_reverse, marker_uidx))) {
      set_allele_freq_buf[markers_read] = set_allele_freqs[marker_uidx];
    } else {
      set_allele_freq_buf[markers_read] = 1.0 - set_allele_freqs[marker_uidx];
    }
    markers_read++;
    marker_idx++;
    marker_uidx++;
  }
}

int32_t load_distance_wts(char* distance_wts_fname, uintptr_t unfiltered_marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t noheader, uint32_t conditional_alloc_exclude, uintptr_t** marker_exclude_ptr, uint32_t* marker_ct_ptr, double** main_weights_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* infile = nullptr;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t line_idx = 0;

  // special case: weight-0 assignment effectively doesn't exist, but we still
  // want to check for repeated IDs there.
  uint32_t zcount = 0;

  int32_t retval = 0;
  uintptr_t* marker_include;
  double* main_weights_tmp;
  double* dptr;
  char* bufptr;
  uint32_t* marker_id_htable;
  double dxx;
  uint32_t marker_id_htable_size;
  uint32_t marker_uidx;
  uint32_t marker_idx;
  uint32_t idlen;
  uint32_t marker_ct;
  if (bigstack_end_calloc_ul(unfiltered_marker_ctl, &marker_include) ||
      bigstack_end_alloc_d(unfiltered_marker_ct, &main_weights_tmp)) {
    goto load_distance_wts_ret_NOMEM;
  }
  retval = alloc_and_populate_id_htable(unfiltered_marker_ct, *marker_exclude_ptr, *marker_ct_ptr, marker_ids, max_marker_id_len, 0, &marker_id_htable_size, &marker_id_htable);
  if (retval) {
    goto load_distance_wts_ret_1;
  }
  if (fopen_checked(distance_wts_fname, "r", &infile)) {
    goto load_distance_wts_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, infile)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, distance_wts_fname);
      goto load_distance_wts_ret_INVALID_FORMAT_2;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (!noheader) {
      noheader = 1;
      continue;
    }
    // variant ID in first column, weight in second
    idlen = strlen_se(bufptr);
    marker_uidx = id_htable_find(bufptr, idlen, marker_id_htable, marker_id_htable_size, marker_ids, max_marker_id_len);
    if (marker_uidx == 0xffffffffU) {
      continue;
    }
    if (is_set(marker_include, marker_uidx)) {
      bufptr[idlen] = '\0';
      LOGPREPRINTFWW("Error: Duplicate variant ID '%s' in --distance-wts file.\n", bufptr);
      goto load_distance_wts_ret_INVALID_FORMAT_2;
    }
    set_bit(marker_uidx, marker_include);
    bufptr = skip_initial_spaces(&(bufptr[idlen]));
    if (is_eoln_kns(*bufptr)) {
      sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --distance-wts file has fewer tokens than expected.\n", line_idx);
      goto load_distance_wts_ret_INVALID_FORMAT_2;
    }
    if (scan_double(bufptr, &dxx)) {
      goto load_distance_wts_ret_INVALID_WEIGHT;
    }
    if (!((dxx >= 0.0) && (dxx != INFINITY))) {
      goto load_distance_wts_ret_INVALID_WEIGHT;
    }
    if (dxx == 0.0) {
      zcount++;
    }
    main_weights_tmp[marker_uidx] = dxx;
  }
  if (!feof(infile)) {
    goto load_distance_wts_ret_READ_FAIL;
  }
  bigstack_reset(bigstack_mark);
  marker_ct = popcount_longs(marker_include, unfiltered_marker_ctl) - zcount;
  if (!marker_ct) {
    logerrprint("Error: No valid nonzero entries in --distance-wts file.\n");
    goto load_distance_wts_ret_INVALID_FORMAT;
  }
  if ((marker_ct != (*marker_ct_ptr))) {
    if (conditional_alloc_exclude) {
      if (bigstack_alloc_ul(unfiltered_marker_ctl, marker_exclude_ptr)) {
	goto load_distance_wts_ret_NOMEM;
      }
    }
    bitarr_invert_copy(marker_include, unfiltered_marker_ct, *marker_exclude_ptr);
    *marker_ct_ptr = marker_ct;
  }
  if (bigstack_alloc_d(marker_ct, main_weights_ptr)) {
    goto load_distance_wts_ret_NOMEM;
  }
  dptr = *main_weights_ptr;
  *marker_ct_ptr = marker_ct;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++) {
    next_set_unsafe_ck(marker_include, &marker_uidx);
    dxx = main_weights_tmp[marker_uidx];
    if (dxx != 0.0) {
      *dptr++ = dxx;
      marker_idx++;
    }
  }
  while (0) {
  load_distance_wts_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_distance_wts_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_distance_wts_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_distance_wts_ret_INVALID_WEIGHT:
    sprintf(g_logbuf, "Error: Invalid weight on line %" PRIuPTR " of --distance-wts file.\n", line_idx);
  load_distance_wts_ret_INVALID_FORMAT_2:
    logerrprintb();
  load_distance_wts_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 load_distance_wts_ret_1:
  bigstack_end_reset(bigstack_end_mark);
  fclose_cond(infile);
  return retval;
}

int32_t calc_rel(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, uint64_t calculation_type, Rel_info* relip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* distance_wts_fname, uint32_t distance_wts_noheader, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* sample_ids, uintptr_t max_sample_id_len, double* set_allele_freqs, double** rel_ibc_ptr, Chrom_info* chrom_info_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t sample_ct = unfiltered_sample_ct - (*sample_exclude_ct_ptr);
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  FILE* outfile = nullptr;
  FILE* out_bin_nfile = nullptr;
  uintptr_t* marker_exclude = marker_exclude_orig;
  uint32_t rel_calc_type = relip->modifier & REL_CALC_MASK;
  int32_t ibc_type = relip->ibc_type;
  int32_t retval = 0;
  uint32_t dist_thread_ct = g_thread_ct;
  uint32_t rel_req = relationship_req(calculation_type);
  uint32_t all_missing_warning = 0;
  int64_t llxx = 0;
  double rel_cutoff = relip->cutoff;
  double* dist_ptr = nullptr;
  double* dptr3 = nullptr;
  double* dptr4 = nullptr;
  double* rel_dists = nullptr;
  double* main_weights = nullptr;
  double* main_weights_ptr = nullptr;
  double* dptr2;
  double set_allele_freq_buf[MULTIPLEX_DIST];
  uint64_t start_offset;
  uint64_t hundredth;
  unsigned char* overflow_buf;
  char* wptr;
  char* fam_id;
  char* sample_id;
  uintptr_t* geno;
  uintptr_t* masks;
  uintptr_t* mmasks;
  double* subset_weights;
  double* rel_ibc;
  uint32_t* mdeptr;
  uint32_t* sample_missing_unwt;
  uint32_t cur_markers_loaded;
  uint32_t win_marker_idx;
  uint32_t is_last_block;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  float fxx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t rel_shape;
  uint32_t min_sample;
  uint32_t max_parallel_sample;
  uint32_t pct;
  unsigned char* gptr;
  unsigned char* gptr2;
  uint32_t* giptr;
  uint32_t* giptr2;
  uintptr_t* glptr2;
  if (distance_wts_fname) {
    logerrprint("Error: --make-{rel,grm-gz,grm-bin} + --distance-wts is currently under\ndevelopment.\n");
    goto calc_rel_ret_1;
  }

  // timing results on the NIH 512-core machine suggest that it's
  // counterproductive to make thread count exceed about n/64
  if (dist_thread_ct > sample_ct / 64) {
    dist_thread_ct = sample_ct / 64;
    if (!dist_thread_ct)  {
      dist_thread_ct = 1;
    }
  }
  // currently must be bottom allocation, since plink() will free it
  if (bigstack_calloc_ui(sample_ct, &sample_missing_unwt)) {
    goto calc_rel_ret_NOMEM;
  }
  g_sample_missing_unwt = sample_missing_unwt;
  g_sample_ct = sample_ct;
  if (dist_thread_ct > sample_ct / 2) {
    dist_thread_ct = sample_ct / 2;
  }
  triangle_fill(sample_ct, dist_thread_ct, parallel_idx, parallel_tot, 1, 1, g_thread_start);
  if (calculation_type & CALC_IBC) {
    uii = sample_ct * 3;
  } else {
    uii = sample_ct;
  }
  if (bigstack_calloc_d(uii, rel_ibc_ptr)) {
    goto calc_rel_ret_NOMEM;
  }
  rel_ibc = *rel_ibc_ptr;
  if (rel_req) {
    llxx = g_thread_start[dist_thread_ct];
    llxx = ((llxx * (llxx - 1)) - (int64_t)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
    if (!(calculation_type & (CALC_PCA | CALC_UNRELATED_HERITABILITY))) {
      // if the memory isn't needed for CALC_PCA or
      // CALC_UNRELATED_HERITABILITY, positioning the missingness matrix here
      // will let us avoid recalculating it if --distance-matrix or --matrix is
      // requested
      if (bigstack_calloc_ui(llxx, &g_missing_dbl_excluded)) {
	goto calc_rel_ret_NOMEM;
      }
    }
    if (bigstack_calloc_d(llxx, &rel_dists)) {
      goto calc_rel_ret_NOMEM;
    }
    g_rel_dists = rel_dists;
  }
  bigstack_mark = g_bigstack_base;
  // bigstack allocations after this point are freed normally
  if (rel_req && (!g_missing_dbl_excluded)) {
    if (bigstack_calloc_ui(llxx, &g_missing_dbl_excluded)) {
      goto calc_rel_ret_NOMEM;
    }
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_rel_ret_READ_FAIL;
  }
  if (bigstack_alloc_ul(sample_ct, &geno) ||
      bigstack_alloc_ul(sample_ct, &mmasks) ||
      bigstack_alloc_uc(MULTIPLEX_REL * unfiltered_sample_ct4, &gptr) ||
      bigstack_alloc_ul(sample_ct, &masks) ||
      bigstack_alloc_d(2048 * BITCT, &subset_weights) ||
      bigstack_alloc_uc(262144, &overflow_buf)) {
    goto calc_rel_ret_NOMEM;
  }
  g_geno = (unsigned char*)geno;
  g_masks = masks;
  g_mmasks = mmasks;
  g_subset_weights = subset_weights;

  // Exclude markers on non-autosomal chromosomes for now.
  retval = conditional_allocate_non_autosomal_markers(chrom_info_ptr, unfiltered_marker_ct, marker_exclude_orig, marker_ct, 1, 1, "relationship matrix calc", &marker_exclude, &uii);
  if (retval) {
    goto calc_rel_ret_1;
  }
  marker_ct -= uii;

  if (distance_wts_fname) {
    retval = load_distance_wts(distance_wts_fname, unfiltered_marker_ct, marker_ids, max_marker_id_len, distance_wts_noheader, (marker_exclude == marker_exclude_orig), &marker_exclude, &marker_ct, &main_weights);
    if (retval) {
      goto calc_rel_ret_1;
    }
  }

  // See comments at the beginning of this file, and those in the main
  // CALC_DISTANCE loop.  The main difference between this calculation and
  // the (nonzero exponent) distance calculation is that we have to pad
  // each marker to 3 bits and use + instead of XOR to distinguish the
  // cases.
  do {
    copy_set_allele_freqs(marker_uidx, marker_exclude, MULTIPLEX_REL, marker_idx, marker_ct, marker_reverse, set_allele_freqs, set_allele_freq_buf);
    if (main_weights) {
      main_weights_ptr = &(main_weights[marker_idx]);
    }
    retval = block_load(bedfile, bed_offset, marker_exclude, marker_ct, MULTIPLEX_REL, unfiltered_sample_ct4, gptr, &marker_uidx, &marker_idx, &cur_markers_loaded);
    if (retval) {
      goto calc_rel_ret_1;
    }
    if (cur_markers_loaded < MULTIPLEX_REL) {
      memset(&(gptr[cur_markers_loaded * unfiltered_sample_ct4]), 0, (MULTIPLEX_REL - cur_markers_loaded) * unfiltered_sample_ct4);
      fill_double_zero(MULTIPLEX_REL - cur_markers_loaded, &(set_allele_freq_buf[cur_markers_loaded]));
    }
    fill_ulong_zero(sample_ct, mmasks);

    is_last_block = (marker_idx == marker_ct);
    for (win_marker_idx = 0; win_marker_idx < cur_markers_loaded; win_marker_idx += MULTIPLEX_REL / 3) {
      fill_ulong_zero(sample_ct, masks);
      sample_idx = 0;
      glptr2 = geno;
      for (sample_uidx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
	next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
	ulii = 0;
	gptr2 = &(gptr[sample_uidx / 4 + win_marker_idx * unfiltered_sample_ct4]);
	uii = (sample_uidx % 4) * 2;
	umm = 0;
	unn = 0;
	for (ukk = 0; ukk < (BITCT / 16); ukk++) {
	  for (ujj = 0; ujj < 5; ujj++) {
	    uljj = (gptr2[umm * unfiltered_sample_ct4] >> uii) & 3;
	    if (uljj == 1) {
	      masks[sample_idx] |= (7 * ONELU) << unn;
	      mmasks[sample_idx] |= ONELU << (win_marker_idx + umm);
	      sample_missing_unwt[sample_idx] += 1;
	    }
	    ulii |= uljj << unn;
	    umm++;
	    unn += 3;
	  }
	  unn++;
	}
	*glptr2++ = ulii;
      }
      ujj = is_last_block && (cur_markers_loaded - win_marker_idx <= MULTIPLEX_REL / 3);
      if (!ujj) {
	ukk = MULTIPLEX_REL / 3;
      } else {
	ukk = cur_markers_loaded - win_marker_idx;
      }
      if (calculation_type & CALC_IBC) {
	for (uii = 0; uii < 3; uii++) {
	  update_rel_ibc(&(rel_ibc[uii * sample_ct]), geno, &(set_allele_freq_buf[win_marker_idx]), main_weights_ptr? (&(main_weights_ptr[win_marker_idx])) : nullptr, uii, sample_ct, ukk);
	}
      } else {
	update_rel_ibc(rel_ibc, geno, &(set_allele_freq_buf[win_marker_idx]), main_weights_ptr? (&(main_weights_ptr[win_marker_idx])) : nullptr, ibc_type, sample_ct, ukk);
      }
      if (rel_req) {
	fill_subset_weights_r(subset_weights, &(set_allele_freq_buf[win_marker_idx]), main_weights_ptr? (&(main_weights_ptr[win_marker_idx])) : nullptr, (ibc_type != -1));
	ulii = 0;
	if (!main_weights_ptr) {
	  if (spawn_threads2(threads, &calc_rel_thread, dist_thread_ct, ujj)) {
	    goto calc_rel_ret_THREAD_CREATE_FAIL;
	  }
	  calc_rel_thread((void*)ulii);
	} else {
	  if (spawn_threads2(threads, &calc_wt_rel_thread, dist_thread_ct, ujj)) {
	    goto calc_rel_ret_THREAD_CREATE_FAIL;
	  }
	  calc_wt_rel_thread((void*)ulii);
	}
	join_threads2(threads, dist_thread_ct, ujj);
      }
    }
    printf("\r%" PRIuPTR " markers complete.", marker_idx);
    fflush(stdout);
  } while (!is_last_block);
  if (rel_req) {
    putc_unlocked('\r', stdout);
    logprint("Relationship matrix calculation complete.\n");
    dist_ptr = rel_dists;
  } else {
    putc_unlocked('\n', stdout);
  }
  dptr2 = rel_ibc;
  if (calculation_type & CALC_IBC) {
    dptr3 = &(rel_ibc[sample_ct]);
    dptr4 = &(rel_ibc[sample_ct * 2]);
  }
  giptr2 = g_missing_dbl_excluded;
  min_sample = g_thread_start[0];
  max_parallel_sample = g_thread_start[dist_thread_ct];
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    uii = marker_ct - sample_missing_unwt[sample_idx];
    if ((!all_missing_warning) && (!uii)) {
      all_missing_warning = 1;
    }
    if ((sample_idx >= min_sample) && (sample_idx < max_parallel_sample)) {
      if (rel_req) {
	giptr = sample_missing_unwt;
	for (ujj = 0; ujj < sample_idx; ujj++) {
	  *dist_ptr /= uii - (*giptr++) + (*giptr2++);
	  dist_ptr++;
	}
      }
    }
    if (calculation_type & CALC_IBC) {
      *dptr2 /= uii;
      dptr2++;
      *dptr3 /= uii;
      dptr3++;
      *dptr4 /= uii;
      dptr4++;
    } else {
      *dptr2 /= uii;
      dptr2++;
    }
  }
  if (calculation_type & CALC_REL_CUTOFF) {
    retval = do_rel_cutoff(calculation_type, rel_cutoff, rel_ibc, sample_exclude, sample_exclude_ct_ptr, outname, outname_end, unfiltered_sample_ct, sample_ids, max_sample_id_len);
    if (retval) {
      goto calc_rel_ret_1;
    }
    sample_ct = unfiltered_sample_ct - *sample_exclude_ct_ptr;
    g_sample_ct = sample_ct;
  }

  if (calculation_type & CALC_IBC) {
    strcpy(outname_end, ".ibc");
    if (fopen_checked(outname, "w", &outfile)) {
      goto calc_rel_ret_OPEN_FAIL;
    }
    dptr2 = rel_ibc;
    dptr3 = &(rel_ibc[sample_ct]);
    dptr4 = &(rel_ibc[sample_ct * 2]);
    if (fputs_checked("FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3\n", outfile)) {
      goto calc_rel_ret_WRITE_FAIL;
    }
    char* textbuf = g_textbuf;
    sample_uidx = 0;
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      fam_id = &(sample_ids[sample_uidx * max_sample_id_len]);
      sample_id = (char*)memchr(fam_id, '\t', max_sample_id_len);
      wptr = memcpyax(textbuf, fam_id, (uintptr_t)(sample_id - fam_id), '\t');
      wptr = strcpyax(wptr, &(sample_id[1]), '\t');
      wptr = uint32toa_x(marker_ct - sample_missing_unwt[sample_idx], '\t', wptr);
      wptr = dtoa_gx(*dptr3++ - 1.0, '\t', wptr);
      wptr = dtoa_gx(*dptr4++ - 1.0, '\t', wptr);
      wptr = dtoa_gx(*dptr2++ - 1.0, '\n', wptr);

      if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
	goto calc_rel_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_rel_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("%s written.\n", outname);
  }
  if (calculation_type & CALC_RELATIONSHIP) {
    rel_shape = rel_calc_type & REL_CALC_SHAPEMASK;
    if (parallel_tot == 1) {
      // nasty --rel-cutoff bug
      max_parallel_sample = sample_ct;
    } else {
      // can't run --rel-cutoff with --parallel, so this is safe
      max_parallel_sample = g_thread_start[dist_thread_ct];
    }
    min_sample = g_thread_start[0];
    if (min_sample == 1) {
      min_sample = 0;
    }
    if (calculation_type & CALC_IBC) {
      dptr2 = &(rel_ibc[((uint32_t)ibc_type) * sample_ct + min_sample]);
    } else {
      dptr2 = &(rel_ibc[min_sample]);
    }
    start_offset = ((uint64_t)min_sample * (min_sample - 1)) / 2;
    hundredth = 1 + (((((uint64_t)max_parallel_sample * (max_parallel_sample + 1)) / 2) - start_offset) / 100);
    if (rel_calc_type & REL_CALC_BIN) {
      pct = 1;
      if (rel_shape == REL_CALC_SQ0) {
	fill_double_zero(sample_ct - 1, (double*)geno);
      }
      strcpy(outname_end, ".rel.bin");
      if (parallel_tot > 1) {
	sprintf(&(outname_end[8]), ".%u", parallel_idx + 1);
      }
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto calc_rel_ret_OPEN_FAIL;
      }
      for (sample_idx = min_sample; sample_idx < max_parallel_sample; sample_idx++) {
	if (fwrite_checked(&(rel_dists[((int64_t)sample_idx * (sample_idx - 1)) / 2 - start_offset]), sample_idx * sizeof(double), outfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	if (fwrite_checked(dptr2++, sizeof(double), outfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	if (rel_shape == REL_CALC_TRI) {
	  if ((((uint64_t)sample_idx + 1) * (sample_idx + 2) / 2 - start_offset) >= hundredth * pct) {
	    pct = (((uint64_t)sample_idx + 1) * (sample_idx + 2) / 2 - start_offset) / hundredth;
	    printf("\rWriting... %u%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (rel_shape == REL_CALC_SQ0) {
	    if (fwrite_checked(geno, (sample_ct - sample_idx - 1) * sizeof(double), outfile)) {
	      goto calc_rel_ret_WRITE_FAIL;
	    }
	  } else {
	    for (uii = sample_idx + 1; uii < sample_ct; uii++) {
	      if (fwrite_checked(&(rel_dists[((uintptr_t)uii * (uii - 1) / 2) + sample_idx - start_offset]), sizeof(double), outfile)) {
		goto calc_rel_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((sample_idx + 1 - min_sample) * 100 >= pct * (max_parallel_sample - min_sample)) {
	    pct = ((sample_idx + 1 - min_sample) * 100) / (max_parallel_sample - min_sample);
	    printf("\rWriting... %u%%", pct++);
	    fflush(stdout);
	  }
	}
      }
      if (fclose_null(&outfile)) {
	goto calc_rel_ret_WRITE_FAIL;
      }
    } else if (rel_calc_type & REL_CALC_GRM_BIN) {
      pct = 1;
      memcpy(outname_end, ".grm.N.bin", 11);
      if (parallel_tot > 1) {
	outname_end[10] = '.';
	uint32toa_x(parallel_idx + 1, '\0', &(outname_end[11]));
      }
      if (fopen_checked(outname, FOPEN_WB, &out_bin_nfile)) {
	goto calc_rel_ret_OPEN_FAIL;
      }
      memcpy(outname_end, ".grm.bin", 9);
      if (parallel_tot > 1) {
	outname_end[8] = '.';
	uint32toa_x(parallel_idx + 1, '\0', &(outname_end[9]));
      }
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto calc_rel_ret_OPEN_FAIL;
      }
      mdeptr = g_missing_dbl_excluded;
      for (sample_idx = min_sample; sample_idx < max_parallel_sample; sample_idx++) {
	dptr3 = &(rel_dists[((int64_t)sample_idx * (sample_idx - 1)) / 2 - start_offset]);
	for (ulii = 0; ulii < sample_idx; ulii++) {
	  fxx = (float)(*dptr3++);
	  fwrite(&fxx, sizeof(float), 1, outfile);
	}
	fxx = (float)(*dptr2++);
	if (fwrite_checked(&fxx, sizeof(float), outfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	uii = marker_ct - sample_missing_unwt[sample_idx];
	for (ujj = 0; ujj < sample_idx; ujj++) {
	  fxx = (float)((int32_t)(uii - sample_missing_unwt[ujj] + (*mdeptr++)));
	  fwrite(&fxx, sizeof(float), 1, out_bin_nfile);
	}
	fxx = (float)((int32_t)uii);
	if (fwrite_checked(&fxx, sizeof(float), out_bin_nfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	if ((((uint64_t)sample_idx + 1) * (sample_idx + 2) / 2 - start_offset) >= hundredth * pct) {
	  pct = (((uint64_t)sample_idx + 1) * (sample_idx + 2) / 2 - start_offset) / hundredth;
	  printf("\rWriting... %u%%", pct++);
	  fflush(stdout);
	}
      }
    } else if (rel_calc_type & REL_CALC_BIN4) {
      // need to downcode all doubles to floats
      pct = 1;
      if (rel_shape == REL_CALC_SQ0) {
	fill_float_zero(sample_ct - 1, (float*)geno);
      }
      // make this .rel.bin4?
      strcpy(outname_end, ".rel.bin");
      if (parallel_tot > 1) {
	sprintf(&(outname_end[8]), ".%u", parallel_idx + 1);
      }
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
	goto calc_rel_ret_OPEN_FAIL;
      }
      for (sample_idx = min_sample; sample_idx < max_parallel_sample; sample_idx++) {
	dptr3 = &(rel_dists[((int64_t)sample_idx * (sample_idx - 1)) / 2 - start_offset]);
	for (ulii = 0; ulii < sample_idx; ulii++) {
	  fxx = (float)(*dptr3++);
	  fwrite(&fxx, sizeof(float), 1, outfile);
	}
	fxx = (float)(*dptr2++);
	if (fwrite_checked(&fxx, sizeof(float), outfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	if (rel_shape == REL_CALC_TRI) {
	  if ((((uint64_t)sample_idx + 1) * (sample_idx + 2) / 2 - start_offset) >= hundredth * pct) {
	    pct = (((uint64_t)sample_idx + 1) * (sample_idx + 2) / 2 - start_offset) / hundredth;
	    printf("\rWriting... %u%%", pct++);
	    fflush(stdout);
	  }
	} else {
	  if (rel_shape == REL_CALC_SQ0) {
	    if (fwrite_checked(geno, (sample_ct - sample_idx - 1) * sizeof(float), outfile)) {
	      goto calc_rel_ret_WRITE_FAIL;
	    }
	  } else {
	    for (ulii = sample_idx + 1; ulii < sample_ct; ulii++) {
	      fxx = (float)(rel_dists[(ulii * (ulii - 1) / 2) + sample_idx - start_offset]);
	      fwrite(&fxx, 4, 1, outfile);
	    }
	  }
	  if ((sample_idx + 1 - min_sample) * 100 >= pct * (max_parallel_sample - min_sample)) {
	    pct = ((sample_idx + 1 - min_sample) * 100) / (max_parallel_sample - min_sample);
	    printf("\rWriting... %u%%", pct++);
	    fflush(stdout);
	  }
	}
      }
      if (fclose_null(&outfile)) {
	goto calc_rel_ret_WRITE_FAIL;
      }
    } else {
      g_cr_marker_ct = marker_ct;
      g_pct = 1;
      g_cr_start_offset = start_offset;
      g_cr_hundredth = hundredth;
      g_cr_sample1idx = min_sample;
      g_cr_sample2idx = 0;
      g_cr_max_sample1idx = max_parallel_sample;
      g_cr_mdeptr = g_missing_dbl_excluded;
      g_cr_dist_ptr = rel_dists;
      g_cr_ibc_ptr = dptr2;
      if (rel_calc_type & REL_CALC_GRM) {
	if (rel_calc_type & REL_CALC_GZ) {
	  if (parallel_tot > 1) {
	    sprintf(outname_end, ".grm.%u.gz", parallel_idx + 1);
	  } else {
	    strcpy(outname_end, ".grm.gz");
	  }
	  parallel_compress(outname, overflow_buf, 0, calc_rel_grm_emitn);
	} else {
	  strcpy(outname_end, ".grm");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[4]), ".%u", parallel_idx + 1);
	  }
	  retval = write_uncompressed(outname, overflow_buf, 0, calc_rel_grm_emitn);
	  if (retval) {
	    goto calc_rel_ret_1;
	  }
	}
      } else {
	if (rel_calc_type & REL_CALC_GZ) {
	  if (parallel_tot > 1) {
	    sprintf(outname_end, ".rel.%u.gz", parallel_idx + 1);
	  } else {
	    strcpy(outname_end, ".rel.gz");
	  }
	} else {
	  strcpy(outname_end, ".rel");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[4]), ".%u", parallel_idx + 1);
	  }
	}
	if (rel_shape == REL_CALC_TRI) {
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, overflow_buf, 0, calc_rel_tri_emitn);
	  } else {
	    retval = write_uncompressed(outname, overflow_buf, 0, calc_rel_tri_emitn);
	    if (retval) {
	      goto calc_rel_ret_1;
	    }
	  }
	} else if (rel_shape == REL_CALC_SQ0) {
	  // if we wanted to port to big-endian...
	  // cptr2 = (char*)(&ulii);
	  // for (uii = 0; uii < sizeof(intptr_t); uii += 2) {
	  //   cptr2[uii] = '\t';
	  //   cptr2[uii + 1] = '0';
	  // }
#ifdef __LP64__
	  ulii = 0x3009300930093009LLU;
#else
	  ulii = 0x30093009LU;
#endif
	  uii = (sample_ct * 2 + sizeof(intptr_t) - 4) / sizeof(intptr_t);
	  glptr2 = geno;
	  for (ujj = 0; ujj < uii; ujj++) {
	    *glptr2++ = ulii;
	  }
	  g_cr_min_sample = min_sample;
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, overflow_buf, 0, calc_rel_sq0_emitn);
	  } else {
	    retval = write_uncompressed(outname, overflow_buf, 0, calc_rel_sq0_emitn);
	    if (retval) {
	      goto calc_rel_ret_1;
	    }
	  }
	} else {
	  g_cr_min_sample = min_sample;
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, overflow_buf, 0, calc_rel_sq_emitn);
	  } else {
	    retval = write_uncompressed(outname, overflow_buf, 0, calc_rel_sq_emitn);
	    if (retval) {
	      goto calc_rel_ret_1;
	    }
	  }
	}
      }
    }
    putc_unlocked('\r', stdout);
    if (!parallel_idx) {
      wptr = strcpya(g_logbuf, "Relationship matrix ");
      if (parallel_tot > 1) {
	wptr = strcpya(wptr, "component ");
      }
      wptr = strcpya(wptr, "written to ");
      wptr = strcpya(wptr, outname);
      strcpy(&(outname_end[4]), ".id");
      retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
      if (retval) {
	goto calc_rel_ret_1;
      }
      sprintf(wptr, " , and IDs written to %s .\n", outname);
    } else {
      sprintf(g_logbuf, "Relationship matrix component written to %s .\n", outname);
    }
    wordwrapb(0);
    logprintb();
  }
  if (all_missing_warning) {
    // note that it's still possible for the relationship matrix to have nans
    // even if this isn't true (empty intersection, or inaccurate 0 MAF), but
    // at least this catches the common case
    if (calculation_type & (CALC_REGRESS_REL | CALC_PCA | CALC_UNRELATED_HERITABILITY)) {
      logerrprint("Error: Sample(s) present with no genotype data.  Use --mind to filter them out.\n");
      retval = RET_INVALID_FORMAT;
    } else {
      logerrprint("Warning: Sample(s) present with no genotype data.  Use --mind to filter them\nout.\n");
    }
  }
  while (0) {
  calc_rel_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_rel_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_rel_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_rel_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_rel_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_rel_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(out_bin_nfile);
  return retval;
}

#ifndef NOLAPACK
int32_t calc_pca(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, Rel_info* relip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* pca_sample_exclude, uintptr_t pca_sample_ct, char* sample_ids, uintptr_t max_sample_id_len, double* set_allele_freqs, Chrom_info* chrom_info_ptr, double* rel_ibc) {
  // similar to mds_plot() in plink_cluster.c
  // todo: provide a randomized approximation algorithm as well.  (This can
  // wait, though; far more important to implement stuff that doesn't already
  // exist.  EIGENSOFT is not *that* hard to use.)
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t pca_sample_ctl2 = QUATERCT_TO_WORDCT(pca_sample_ct);
  uintptr_t proj_sample_ct = sample_ct - pca_sample_ct;
  uintptr_t proj_sample_ctl2 = QUATERCT_TO_WORDCT(proj_sample_ct);
  uintptr_t final_mask = get_final_mask(pca_sample_ct);
  double nz = 0.0;
  double zz = -1.0;
  uintptr_t* loadbuf_proj = nullptr;
  uintptr_t* sample_exclude_proj = nullptr;
  double* proj_sample_loadings = nullptr;
  double* proj_allhom_wts = nullptr;
  uint32_t* proj_missing_cts = nullptr;
  uint32_t write_headers = relip->modifier & REL_PCA_HEADER;
  uint32_t pc_ct = relip->pc_ct;
  uint32_t var_wts = relip->modifier & REL_PCA_VAR_WTS;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t ibc_type = relip->ibc_type;
  int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  int32_t retval = 0;
  __CLPK_integer info = 0;
  __CLPK_integer lwork = -1;
  __CLPK_integer liwork = -1;
  char jobz = 'V';
  char range = 'I';
  char uplo = 'U';
  char delimiter = (relip->modifier & REL_PCA_TABS)? '\t' : ' ';
  double var_wt_incr[4];
  double* rel_dists;
  double* main_matrix;
  double* out_w;
  double* out_z;
  double* work;
  double* sample_loadings;
  double* cur_var_wts;
  double* pc_sums;
  double* eigval_inv_sqrts;
  double* dptr;
  double* dptr2;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* ulptr;
  char* cptr;
  char* wptr_start;
  char* wptr;
  double optim_lwork;
  double dxx;
  double dyy;
  uintptr_t chrom_end;
  uintptr_t marker_uidx;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t ulii;
  __CLPK_integer mdim;
  __CLPK_integer i1;
  __CLPK_integer i2;
  __CLPK_integer out_m;
  __CLPK_integer ldz;
  __CLPK_integer optim_liwork;
  __CLPK_integer* iwork;
  __CLPK_integer* isuppz;
  uint32_t chrom_fo_idx;
  uint32_t pc_idx;
  uint32_t uii;
  uint32_t ujj;
  int32_t chrom_idx;
  // calc_rel() already verified that diploid data is present
  g_missing_dbl_excluded = nullptr;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
  if ((pc_ct > pca_sample_ct) || (pc_ct > marker_ct)) {
    if (pca_sample_ct <= marker_ct) {
      pc_ct = pca_sample_ct;
      sprintf(g_logbuf, "Warning: calculating %u PCs, since there are only %u samples.\n", pc_ct, pc_ct);
    } else {
      pc_ct = marker_ct;
      sprintf(g_logbuf, "Warning: calculating %u PCs, since there are only %u autosomal markers.\n", pc_ct, pc_ct);
    }
    logerrprintb();
  }
  ulii = (pca_sample_ct * (pca_sample_ct + 1)) / 2;
  if (ulii < sample_ct * pc_ct) {
    ulii = sample_ct * pc_ct;
  }
  rel_dists = g_rel_dists;
  bigstack_reset(rel_dists);
  if (bigstack_alloc_d(ulii, &rel_dists) ||
      bigstack_alloc_d(pca_sample_ct * pca_sample_ct, &main_matrix)) {
    goto calc_pca_ret_NOMEM;
  }
  ulii = ((pca_sample_ct - 1) * (pca_sample_ct - 2)) >> 1;
  for (sample_idx = pca_sample_ct - 1; sample_idx;) {
    memcpy(&(main_matrix[sample_idx * pca_sample_ct]), &(rel_dists[ulii]), sample_idx * sizeof(double));
    ulii -= (--sample_idx);
  }
  if (calculation_type & CALC_IBC) {
    dptr = &(rel_ibc[ibc_type * pca_sample_ct]);
  } else {
    dptr = rel_ibc;
  }
  ulii = pca_sample_ct + 1;
  for (sample_idx = 0; sample_idx < pca_sample_ct; sample_idx++) {
    main_matrix[sample_idx * ulii] = *dptr++;
  }
  ulii = pca_sample_ct;
  fputs("[extracting eigenvalues and eigenvectors]", stdout);
  fflush(stdout);
  mdim = ulii;
  i2 = mdim;
  i1 = i2 + 1 - pc_ct;
  if (bigstack_calloc_d(pc_ct, &out_w) ||
      bigstack_calloc_d(pc_ct * ulii, &out_z)) {
    goto calc_pca_ret_NOMEM;
  }
  isuppz = (__CLPK_integer*)bigstack_alloc(2 * pc_ct * sizeof(__CLPK_integer));
  if (!isuppz) {
    goto calc_pca_ret_NOMEM;
  }
  fill_int_zero(2 * pc_ct * (sizeof(__CLPK_integer) / sizeof(int32_t)), (int32_t*)isuppz);
  ldz = mdim;

  dsyevr_(&jobz, &range, &uplo, &mdim, main_matrix, &mdim, &nz, &nz, &i1, &i2, &zz, &out_m, out_w, out_z, &ldz, isuppz, &optim_lwork, &lwork, &optim_liwork, &liwork, &info);
  lwork = (int32_t)optim_lwork;
  if (bigstack_calloc_d(lwork, &work)) {
    goto calc_pca_ret_NOMEM;
  }
  liwork = optim_liwork;
  iwork = (__CLPK_integer*)bigstack_alloc(liwork * sizeof(__CLPK_integer));
  if (!iwork) {
    goto calc_pca_ret_NOMEM;
  }
  fill_int_zero(liwork * (sizeof(__CLPK_integer) / sizeof(int32_t)), (int32_t*)iwork);
  dsyevr_(&jobz, &range, &uplo, &mdim, main_matrix, &mdim, &nz, &nz, &i1, &i2, &zz, &out_m, out_w, out_z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
  if ((info != 0) || (out_w[0] != out_w[0])) {
    logprint("\n");
    logerrprint("Error: Failed to extract eigenvector(s) from GRM.\n");
    goto calc_pca_ret_DEGENERATE_DATA;
  }

  // * out_w[0..(pc_ct-1)] contains eigenvalues
  // * out_z[(ii*ulii)..(ii*ulii + ulii - 1)] is eigenvector corresponding to
  //   out_w[ii]
  sample_loadings = rel_dists; // PC-major
  dptr = sample_loadings;
  // transpose out_z to sample-major.  eigenvals/vecs still "backwards"
  for (sample_idx = 0; sample_idx < pca_sample_ct; sample_idx++) {
    dptr2 = &(out_z[sample_idx]);
    for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
      *dptr++ = dptr2[pc_idx * pca_sample_ct];
    }
  }
  bigstack_reset(out_z);
  if (var_wts || proj_sample_ct) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
        bigstack_alloc_ul(pca_sample_ctl2, &loadbuf) ||
        bigstack_alloc_d(pc_ct, &cur_var_wts) ||
        bigstack_calloc_d(pc_ct, &pc_sums) ||
        bigstack_alloc_d(pc_ct, &eigval_inv_sqrts)) {
      goto calc_pca_ret_NOMEM;
    }
    loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
    loadbuf[pca_sample_ctl2 - 1] = 0;
    dptr = sample_loadings;
    for (sample_idx = 0; sample_idx < pca_sample_ct; sample_idx++) {
      for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
	pc_sums[pc_idx] += *dptr++;
      }
    }
    for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
      eigval_inv_sqrts[pc_idx] = sqrt(1 / out_w[pc_idx]);
    }
    if (proj_sample_ct) {
      if (bigstack_alloc_ul(unfiltered_sample_ctl, &sample_exclude_proj) ||
          bigstack_alloc_ul(proj_sample_ctl2, &loadbuf_proj) ||
          bigstack_calloc_d(proj_sample_ct * pc_ct, &proj_sample_loadings) ||
	  bigstack_calloc_d(pc_ct, &proj_allhom_wts) ||
          bigstack_calloc_ui(proj_sample_ct, &proj_missing_cts)) {
	goto calc_pca_ret_NOMEM;
      }
      memcpy(sample_exclude_proj, sample_exclude, unfiltered_sample_ctl * sizeof(intptr_t));
      bitvec_ornot(pca_sample_exclude, unfiltered_sample_ctl, sample_exclude_proj);
      zero_trailing_bits(unfiltered_sample_ct, sample_exclude_proj);
      loadbuf_proj[proj_sample_ctl2 - 1] = 0;
    }
    char* textbuf_mid = &(g_textbuf[MAXLINELEN]);
    if (var_wts) {
      memcpy(outname_end, ".eigenvec.var", 14);
      if (fopen_checked(outname, "w", &outfile)) {
	goto calc_pca_ret_OPEN_FAIL;
      }
      if (write_headers) {
	wptr = memcpyl3a(g_textbuf, "CHR");
	*wptr++ = delimiter;
	wptr = memcpyl3a(wptr, "VAR");
	// design fix: include allele codes
	*wptr++ = delimiter;
	wptr = strcpya(wptr, "A1");
	*wptr++ = delimiter;
	wptr = strcpya(wptr, "A2");
	for (pc_idx = 1; pc_idx <= pc_ct; pc_idx++) {
	  *wptr++ = delimiter;
	  wptr = memcpya(wptr, "PC", 2);
	  wptr = uint32toa(pc_idx, wptr);
	}
	*wptr++ = '\n';
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto calc_pca_ret_WRITE_FAIL;
	}
      }
    }
    var_wt_incr[0] = 0;
    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if (is_set(chrom_info_ptr->haploid_mask, chrom_idx) || (chrom_idx == mt_code)) {
	continue;
      }
      marker_uidx = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
      wptr_start = chrom_name_write(chrom_info_ptr, chrom_idx, g_textbuf);
      *wptr_start++ = delimiter;
      if (marker_uidx < chrom_end) {
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto calc_pca_ret_READ_FAIL;
	}
	do {
	  if (IS_SET(marker_exclude, marker_uidx)) {
	    marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	    if (marker_uidx == chrom_end) {
	      break;
	    }
	    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	      goto calc_pca_ret_READ_FAIL;
	    }
	  }
	  // if projecting, loadbuf_raw contains raw data, so we can just
	  // follow up with copy_quaterarr_nonempty_subset_excl() and
	  // reverse_loadbuf()
	  if (load_and_collapse(unfiltered_sample_ct, pca_sample_ct, pca_sample_exclude, final_mask, IS_SET(marker_reverse, marker_uidx), bedfile, loadbuf_raw, loadbuf)) {
	    goto calc_pca_ret_READ_FAIL;
	  }
	  // Variant weight matrix = X^T * S * D^{-1/2}, where X^T is the
	  // variance-standardized genotype matrix, S is the sample weight
	  // matrix, and D is a diagonal eigenvalue matrix.
	  fill_double_zero(pc_ct, cur_var_wts);
	  dxx = set_allele_freqs[marker_uidx];
	  dyy = sqrt(1 / (2 * dxx * (1.0 - dxx)));
	  ulptr = loadbuf;

	  var_wt_incr[1] = dyy; // het
	  var_wt_incr[2] = 2 * (1.0 - dxx) * dyy; // missing
	  var_wt_incr[3] = 2 * dyy; // hom minor
	  sample_uidx = 0; // using this as an offset here
	  do {
	    ulii = ~(*ulptr++);
	    if (sample_uidx + BITCT2 > pca_sample_ct) {
	      ulii &= (ONELU << ((pca_sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
	    }
	    while (ulii) {
	      uii = CTZLU(ulii) & (BITCT - 2);
	      ujj = (ulii >> uii) & 3;
	      sample_idx = sample_uidx + (uii / 2);
	      dxx = var_wt_incr[ujj];
	      dptr = cur_var_wts;
	      dptr2 = &(sample_loadings[sample_idx * pc_ct]);
	      for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
		*dptr += dxx * (*dptr2++);
		dptr++;
	      }
	      ulii &= ~((3 * ONELU) << uii);
	    }
	    sample_uidx += BITCT2;
	  } while (sample_uidx < pca_sample_ct);
	  for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
	    cur_var_wts[pc_idx] -= pc_sums[pc_idx] * var_wt_incr[2];
	    cur_var_wts[pc_idx] *= eigval_inv_sqrts[pc_idx];
	  }
	  if (proj_sample_ct) {
	    copy_quaterarr_nonempty_subset_excl(loadbuf_raw, sample_exclude_proj, unfiltered_sample_ct, proj_sample_ct, loadbuf_proj);
	    if (IS_SET(marker_reverse, marker_uidx)) {
	      reverse_loadbuf(proj_sample_ct, (unsigned char*)loadbuf_proj);
	    }
	    ulptr = loadbuf_proj;
	    sample_uidx = 0;
	    do {
	      ulii = ~(*ulptr++);
	      if (sample_uidx + BITCT2 > proj_sample_ct) {
		ulii &= (ONELU << ((proj_sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
	      }
	      while (ulii) {
		uii = CTZLU(ulii) & (BITCT - 2);
                ujj = (ulii >> uii) & 3;
		sample_idx = sample_uidx + (uii / 2);
		dxx = var_wt_incr[ujj];
		dptr = &(proj_sample_loadings[sample_idx * pc_ct]);
		dptr2 = cur_var_wts;
	        for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
		  *dptr += dxx * cur_var_wts[pc_idx];
		  dptr++;
		}
		if (ujj == 2) {
		  proj_missing_cts[sample_idx] += 1;
		}
                ulii &= ~((3 * ONELU) << uii);
	      }
	      sample_uidx += BITCT2;
	    } while (sample_uidx < proj_sample_ct);
	    dptr = proj_allhom_wts;
	    dxx = var_wt_incr[2];
	    for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
	      *dptr += dxx * cur_var_wts[pc_idx];
	      dptr++;
	    }
	  }
	  if (var_wts) {
	    wptr = strcpya(wptr_start, &(marker_ids[marker_uidx * max_marker_id_len]));
	    *wptr++ = delimiter;
	    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	      goto calc_pca_ret_WRITE_FAIL;
	    }
	    fputs(marker_allele_ptrs[marker_uidx * 2], outfile);
	    putc_unlocked(delimiter, outfile);
	    if (fputs_checked(marker_allele_ptrs[marker_uidx * 2 + 1], outfile)) {
	      goto calc_pca_ret_WRITE_FAIL;
	    }
	    wptr = textbuf_mid;
	    dptr = &(cur_var_wts[pc_ct]);
	    dptr2 = cur_var_wts;
	    do {
	      *wptr++ = delimiter;
	      wptr = dtoa_g(*(--dptr), wptr);
	    } while (dptr > dptr2);
	    *wptr++ = '\n';
	    if (fwrite_checked(textbuf_mid, wptr - textbuf_mid, outfile)) {
	      goto calc_pca_ret_WRITE_FAIL;
	    }
	  }
	} while (++marker_uidx < chrom_end);
      }
    }
    if (var_wts) {
      if (fclose_null(&outfile)) {
        goto calc_pca_ret_WRITE_FAIL;
      }
    }
    if (proj_sample_ct) {
      dptr = proj_sample_loadings;
      for (sample_idx = 0; sample_idx < proj_sample_ct; sample_idx++) {
	dxx = 1.0 / ((double)((int32_t)(marker_ct - proj_missing_cts[sample_idx])));
	for (pc_idx = 0; pc_idx < pc_ct; pc_idx++) {
	  *dptr = ((*dptr) - proj_allhom_wts[pc_idx]) * eigval_inv_sqrts[pc_idx] * dxx;
	  dptr++;
	}
      }
      // now fill final sample_loadings from back to front
      sample_uidx = unfiltered_sample_ct - 1;
      sample_idx = sample_ct - 1;
      ulii = pca_sample_ct;
      while (1) {
	if (IS_SET(sample_exclude, sample_uidx)) {
	  sample_uidx = prev_unset_unsafe(sample_exclude, sample_uidx);
	}
        if (IS_SET(sample_exclude_proj, sample_uidx)) {
	  memcpy(&(sample_loadings[sample_idx * pc_ct]), &(sample_loadings[(--ulii) * pc_ct]), pc_ct * sizeof(double));
	} else {
	  memcpy(&(sample_loadings[sample_idx * pc_ct]), &(proj_sample_loadings[(--proj_sample_ct) * pc_ct]), pc_ct * sizeof(double));
	  if (!proj_sample_ct) {
	    break;
	  }
	}
        sample_uidx--;
	sample_idx--;
      }
    }
  }

  memcpy(outname_end, ".eigenval", 10);
  if (fopen_checked(outname, "w", &outfile)) {
    goto calc_pca_ret_OPEN_FAIL;
  }
  wptr = g_textbuf;
  pc_idx = pc_ct;
  do {
    pc_idx--;
    wptr = dtoa_gx(out_w[pc_idx], '\n', wptr);
  } while (pc_idx);
  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
    goto calc_pca_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile)) {
    goto calc_pca_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".eigenvec", 10);
  if (fopen_checked(outname, "w", &outfile)) {
    goto calc_pca_ret_OPEN_FAIL;
  }
  if (write_headers) {
    wptr = memcpyl3a(g_textbuf, "FID");
    *wptr++ = delimiter;
    wptr = memcpyl3a(wptr, "IID");
    for (pc_idx = 1; pc_idx <= pc_ct; pc_idx++) {
      *wptr++ = delimiter;
      wptr = memcpya(wptr, "PC", 2);
      wptr = uint32toa(pc_idx, wptr);
    }
    *wptr++ = '\n';
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto calc_pca_ret_WRITE_FAIL;
    }
  }
  for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_idx++, sample_uidx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    if (delimiter == '\t') {
      wptr = strcpya(g_textbuf, cptr);
    } else {
      wptr_start = (char*)memchr(cptr, '\t', max_sample_id_len);
      wptr = memcpya(g_textbuf, cptr, wptr_start - cptr);
      *wptr++ = ' ';
      wptr = strcpya(wptr, &(wptr_start[1]));
    }
    dptr = &(sample_loadings[(sample_idx + 1) * pc_ct]);
    dptr2 = &(sample_loadings[sample_idx * pc_ct]);
    do {
      *wptr++ = delimiter;
      wptr = dtoa_g(*(--dptr), wptr);
    } while (dptr > dptr2);
    *wptr++ = '\n';
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
      goto calc_pca_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto calc_pca_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  putc_unlocked('\r', stdout);
  if (var_wts) {
    LOGPRINTFWW("--pca: Results saved to %s.eigenval , %s.eigenvec , and %s.eigenvec.var .\n", outname, outname, outname);
  } else {
    LOGPRINTFWW("--pca: Results saved to %s.eigenval and %s.eigenvec .\n", outname, outname);
  }
  while (0) {
  calc_pca_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_pca_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_pca_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_pca_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_pca_ret_DEGENERATE_DATA:
    retval = RET_DEGENERATE_DATA;
    break;
  }
  return retval;
}
#endif

int32_t calc_ibm(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uint32_t marker_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Chrom_info* chrom_info_ptr) {
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t dist_thread_ct = g_thread_ct;
  int32_t retval = 0;
  uintptr_t* marker_exclude = marker_exclude_orig;
  uint32_t* giptr = nullptr;
  unsigned char* bigstack_mark;
  unsigned char* bedbuf;
  unsigned char* gptr;
  uint32_t* sample_missing_unwt;
  uintptr_t* mmasks;
  uintptr_t sample_uidx;
  uintptr_t ulii;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uintptr_t* glptr;
  int64_t llxx;
  g_sample_ct = sample_ct;
  if (dist_thread_ct > sample_ct / 32) {
    dist_thread_ct = sample_ct / 32;
    if (!dist_thread_ct)  {
      dist_thread_ct = 1;
    }
  }
  triangle_fill(sample_ct, dist_thread_ct, 0, 1, 1, 1, g_thread_start);
  llxx = g_thread_start[dist_thread_ct];
  llxx = (llxx * (llxx - 1)) / 2;
  if (bigstack_calloc_ui(llxx, &g_missing_dbl_excluded) ||
      bigstack_calloc_ui(sample_ct, &sample_missing_unwt)) {
    goto calc_ibm_ret_NOMEM;
  }
  g_sample_missing_unwt = sample_missing_unwt;
  bigstack_mark = g_bigstack_base;

  if (bigstack_alloc_ul(sample_ct, &mmasks) ||
      bigstack_alloc_uc(MULTIPLEX_DIST * unfiltered_sample_ct4, &bedbuf)) {
    goto calc_ibm_ret_NOMEM;
  }
  g_mmasks = mmasks;
  fseeko(bedfile, bed_offset, SEEK_SET);
  retval = conditional_allocate_non_autosomal_markers(chrom_info_ptr, unfiltered_marker_ct, marker_exclude_orig, marker_ct, 1, 1, "IBM calculation", &marker_exclude, &uii);
  if (retval) {
    goto calc_ibm_ret_1;
  }
  marker_ct -= uii;
  do {
    retval = block_load(bedfile, bed_offset, marker_exclude, marker_ct, MULTIPLEX_DIST, unfiltered_sample_ct4, bedbuf, &marker_uidx, &marker_idx, &ujj);
    if (retval) {
      goto calc_ibm_ret_1;
    }
    if (ujj < MULTIPLEX_DIST) {
      memset(&(bedbuf[ujj * unfiltered_sample_ct4]), 0, (MULTIPLEX_DIST - ujj) * unfiltered_sample_ct4);
    }
    is_last_block = (marker_idx == marker_ct);
    for (ukk = 0; ukk < ujj; ukk += BITCT) {
      glptr = mmasks;
      giptr = sample_missing_unwt;
      for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
	if (!IS_SET(sample_exclude, sample_uidx)) {
	  unn = (sample_uidx % 4) * 2;
	  ulii = 0;
	  gptr = &(bedbuf[sample_uidx / 4 + ukk * unfiltered_sample_ct4]);
	  for (umm = 0; umm < BITCT2; umm++) {
	    if (((gptr[umm * unfiltered_sample_ct4] >> unn) & 3) == 1) {
	      ulii |= ONELU << umm;
	      *giptr += 1;
	    }
	  }
	  *glptr = ulii;
	  ulii = 0;
	  gptr = &(bedbuf[sample_uidx / 4 + (ukk + BITCT2) * unfiltered_sample_ct4]);
	  for (umm = 0; umm < BITCT2; umm++) {
	    if (((gptr[umm * unfiltered_sample_ct4] >> unn) & 3) == 1) {
	      ulii |= ONELU << umm;
	      *giptr += 1;
	    }
	  }
	  *glptr++ |= ulii << BITCT2;
	  giptr++;
	}
      }

      uii = is_last_block && (ukk + BITCT >= ujj);
      if (spawn_threads2(threads, &calc_missing_thread, dist_thread_ct, uii)) {
	goto calc_ibm_ret_THREAD_CREATE_FAIL;
      }
      ulii = 0;
      calc_missing_thread((void*)ulii);
      join_threads2(threads, dist_thread_ct, uii);
    }

    printf("\r%" PRIuPTR " markers complete.", marker_idx);
    fflush(stdout);
  } while (!is_last_block);
  putc_unlocked('\r', stdout);
  bigstack_reset(bigstack_mark);
  while (0) {
  calc_ibm_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_ibm_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_ibm_ret_1:
  // caller will free memory if there was an error
  return retval;
}

int32_t calc_distance(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* read_dists_fname, char* distance_wts_fname, double distance_exp, uint64_t calculation_type, uint32_t dist_calc_type, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, Chrom_info* chrom_info_ptr) {
  // if calculation_type == 0, this must perform the basic unweighted
  // computation and not write to disk.
  FILE* outfile = nullptr;
  FILE* outfile2 = nullptr;
  FILE* outfile3 = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uint64_t dists_alloc = 0;
  uint32_t missing_wt_needed = ((calculation_type & CALC_DISTANCE) || ((!read_dists_fname) && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE)))) && (!(dist_calc_type & DISTANCE_FLAT_MISSING));
  uint32_t unwt_needed = 0;
  uint32_t marker_weight_sum = 0;
  int32_t retval = 0;
  uintptr_t* marker_exclude = marker_exclude_orig;
  uint32_t* dist_missing_wts_i = nullptr;
  uint32_t* sample_missing = nullptr;
  uint32_t* sample_missing_unwt = nullptr;
  uint32_t* giptr = nullptr;
  uint32_t* giptr2 = nullptr;
  char* writebuf = nullptr;
  double* main_weights = nullptr;
  double* subset_weights = nullptr;
  uint32_t dist_thread_ct = g_thread_ct;
  double set_allele_freq_buf[MULTIPLEX_DIST];
  uint32_t wtbuf[MULTIPLEX_DIST];
  uintptr_t* geno;
  uintptr_t* masks;
  uintptr_t* mmasks;
  char* wptr;
  unsigned char* bigstack_mark;
  unsigned char* bedbuf;
  unsigned char* gptr;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t pct;
  uint32_t tstc;
  int32_t* iptr;
  uint32_t* giptr3;
  uintptr_t* glptr;
  uintptr_t* glptr2;
  uintptr_t* glptr3;
  double* dist_missing_wts;
  double* dptr2;
  double marker_weight_sum_d;
  double dxx;
  double dyy;
  uint32_t multiplex;
  int64_t llxx;
  g_sample_ct = sample_ct;
  if (dist_thread_ct > sample_ct / 32) {
    dist_thread_ct = sample_ct / 32;
    if (!dist_thread_ct)  {
      dist_thread_ct = 1;
    }
  }
  triangle_fill(sample_ct, dist_thread_ct, parallel_idx, parallel_tot, 1, 1, g_thread_start);
  llxx = g_thread_start[dist_thread_ct];
  llxx = ((llxx * (llxx - 1)) - (int64_t)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
  dists_alloc = llxx * sizeof(double);
#ifndef __LP64__
  if (dists_alloc > 0x7fffffff) {
    goto calc_distance_ret_NOMEM;
  }
#endif
  if ((calculation_type & (CALC_PLINK1_DISTANCE_MATRIX | CALC_PLINK1_IBS_MATRIX)) || (dist_calc_type & DISTANCE_FLAT_MISSING)) {
    if (bigstack_calloc_ui(llxx, &g_missing_dbl_excluded) ||
        bigstack_calloc_ui(sample_ct, &sample_missing_unwt)) {
      goto calc_distance_ret_NOMEM;
    }
    g_sample_missing_unwt = sample_missing_unwt;
    unwt_needed = 1;
  } else {
    // defensive
    g_missing_dbl_excluded = nullptr;
  }
  // Additional + CACHELINE is to fix aliasing bug that shows up with -O2 in
  // some cases.
  if (bigstack_alloc_d(dists_alloc + CACHELINE_DBL, &g_dists)) {
    goto calc_distance_ret_NOMEM;
  }
  // stack allocations before this point must be freed by the caller.
  bigstack_mark = g_bigstack_base;
  if (missing_wt_needed) {
    if (bigstack_calloc_ui(llxx, &g_missing_tot_weights) ||
        bigstack_calloc_ui(sample_ct, &sample_missing)) {
      goto calc_distance_ret_NOMEM;
    }
    g_sample_missing = sample_missing;
  } else {
    g_missing_tot_weights = nullptr;
  }

  ujj = distance_wts_fname || (distance_exp != 0.0); // special weights?
  if (!ujj) {
    g_idists = (int32_t*)((char*)bigstack_mark - round_up_pow2(llxx * sizeof(int32_t), CACHELINE));
    fill_int_zero(llxx, g_idists);
  } else {
    fill_double_zero(llxx, g_dists);
  }

  retval = conditional_allocate_non_autosomal_markers(chrom_info_ptr, unfiltered_marker_ct, marker_exclude_orig, marker_ct, 1, 1, "distance matrix calc", &marker_exclude, &uii);
  if (retval) {
    goto calc_distance_ret_1;
  }
  marker_ct -= uii;

  if (distance_wts_fname) {
    retval = load_distance_wts(distance_wts_fname, unfiltered_marker_ct, marker_ids, max_marker_id_len, dist_calc_type & DISTANCE_WTS_NOHEADER, (marker_exclude == marker_exclude_orig), &marker_exclude, &marker_ct, &main_weights);
    if (retval) {
      goto calc_distance_ret_1;
    }
  }

  // stack allocations past this point are freed BEFORE results are written.
  if (bigstack_alloc_ul(ujj? sample_ct : (sample_ct * (MULTIPLEX_2DIST / BITCT)), &masks) ||
      bigstack_alloc_ul(sample_ct, &mmasks)) {
    goto calc_distance_ret_NOMEM;
  }

  // Load or compute nonuniform marker weighting scheme.
  if (distance_exp != 0.0) {
    if (bigstack_alloc_d(marker_ct, &main_weights)) {
      goto calc_distance_ret_NOMEM;
    }
    for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      dxx = set_allele_freqs[marker_uidx];
      dyy = 2 * dxx * (1.0 - dxx);
      if (dyy != 0.0) {
	dyy = pow(dyy, -distance_exp);
      }
      main_weights[marker_idx] = dyy;
    }
  }
  // Now compute missing observation weights.  (Note that these are usually not
  // the same as the raw marker weights: for instance, a missing observation at
  // a MAF-0 marker has no weight at all.)
  if (missing_wt_needed) {
    // hack: overwrite dist_missing_wts while populating dist_missing_wts_i.
    // CACHELINE padding added to reduce risk of an aliasing problem.
    if (bigstack_alloc_ui(CACHELINE_INT32, &dist_missing_wts_i) ||
        bigstack_alloc_d(marker_ct, &dist_missing_wts)) {
      goto calc_distance_ret_NOMEM;
    }
    dyy = 0.0; // raw weight sum
    for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      // assume HWE, compute expected contribution to distance statistic:
      //   expected minor allele obs: 2 * maf
      //   P(0 copies) = (1 - maf) * (1 - maf)
      //   P(1 copy)   = 2 * maf * (1 - maf)
      //   P(2 copies) = maf * maf
      //   frequency of distance-1 pairs:
      //       freq[0-1 pair] + freq[1-2 pair]
      //     =   2 * (1 - maf) * (1 - maf) * 2 * maf * (1 - maf)
      //       + 2 * 2 * maf * (1 - maf) * maf * maf
      //     =   4 * maf * (1 - maf) * (maf * maf + (1 - maf) * (1 - maf))
      //         4 * maf * (1 - maf) * (2 * maf * maf - 2 * maf + 1)
      //   frequency of distance-2 pairs:
      //     2 * (1 - maf) * (1 - maf) * maf * maf
      //   expected distance:
      //       4 * maf * (1 - maf) * (2 * maf * maf - 2 * maf + 1
      //                              + maf * (1 - maf))
      //     = 4 * maf * (1 - maf) * (maf * maf - maf + 1)
      //     constant factor doesn't matter here
      dxx = set_allele_freqs[marker_uidx];
      if ((dxx != 0.0) && (dxx != 1.0)) {
	dxx = dxx * (1.0 - dxx) * (dxx * dxx - dxx + 1);
	if (main_weights) {
	  dxx *= main_weights[marker_idx];
	}
      }
      dist_missing_wts[marker_idx] = dxx;
      dyy += dxx;
    }

    // now normalize to sum to just under 2^32.  (switch to 2^64 if/when 32-bit
    // performance becomes less important than accuracy on 50+ million marker
    // sets.)
    // subtract marker_ct to guard against rounding-driven overflow
    dyy = (4294967296.0 - ((double)((intptr_t)marker_ct))) / dyy;
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      uii = (uint32_t)(dist_missing_wts[marker_idx] * dyy + 0.5);
      marker_weight_sum += uii;
      dist_missing_wts_i[marker_idx] = uii;
    }
  }
  marker_weight_sum_d = (double)marker_weight_sum;

  if (!main_weights) {
    multiplex = MULTIPLEX_DIST;
    bigstack_alloc_ul(sample_ct * (MULTIPLEX_2DIST / BITCT), &geno);
  } else {
    multiplex = MULTIPLEX_DIST_EXP;
    bigstack_alloc_ul(sample_ct, &geno);
  }
  if (!geno) {
    goto calc_distance_ret_NOMEM;
  }
  g_geno = (unsigned char*)geno;
  g_masks = masks;
  g_mmasks = mmasks;

  if (bigstack_alloc_uc(multiplex * unfiltered_sample_ct4, &bedbuf)) {
    goto calc_distance_ret_NOMEM;
  }
  if (main_weights) {
#ifdef __LP64__
    if (bigstack_alloc_d(45056, &subset_weights)) {
      goto calc_distance_ret_NOMEM;
    }
#else
    if (bigstack_alloc_d(32768, &subset_weights)) {
      goto calc_distance_ret_NOMEM;
    }
    g_subset_weights_i = wtbuf;
#endif
    g_subset_weights = subset_weights;
  }
  fseeko(bedfile, bed_offset, SEEK_SET);
  marker_uidx = 0;
  marker_idx = 0;
  do {
    for (ujj = 0; ujj < multiplex; ujj++) {
      set_allele_freq_buf[ujj] = 0.5;
    }
    fill_int_zero(multiplex, (int32_t*)wtbuf);

    // For each pair (g_j, g_k) of 2-bit PLINK genotypes, we perform the
    // following operations:
    //
    // 1. XOR each genotype with 01.  This shuffles the genotype
    // representation to:
    //    00 = missing (important for simplifying step 2)
    //    01 = homozygote 1
    //    10 = homozygote 2
    //    11 = heterozygote
    //
    // 2. Next, compute
    //    mask_i := ((g_i | (g_i >> 1)) & 01) * 11
    // which is 00 whenever g_i is missing, and 11 otherwise.
    //
    // 3. Then, (g_j ^ g_k) & (mask_j & mask_k) distinguishes the
    // possible distances between the genotypes:
    //    - It's equal to zero iff either g_j == g_k or one/both is
    //      missing.  It's fine for these cases to overlap since either
    //      way we do not want to increment the numerator of our final
    //      distance.  (We can handle the effect of missingness on the
    //      denominator outside the main loop.)
    //    - It's equal to 01 or 10 iff neither is missing and exactly
    //      one is a heterozygote.
    //    - It's equal to 11 iff one is homozygote rare and the other is
    //      homozygote common.
    //
    // 4. Finally, we perform the update
    //    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_31(x_31)
    // in the nonzero exponent case, or
    //    A_{jk} := A_{jk} + f(x_0) + f(x_1) + ... + f(x_959)
    // in the zero exponent case.
    //
    // For nonzero exponents, we structure the update as
    //    A_{jk} := A_{jk} + f_{0-6}(x_{0-6}) + f_{7-13}(x_{7-13}) +
    //              f_{14-19}(x_{14-19}) + f_{20-25}(x_{20-25}) +
    //              f_{26-31}(x_{26-31})
    // which requires 352 KB of table space.  (This is a conservative
    // choice; the 2 MB 8-8-8-8 table would work better on some newer
    // systems.)
    //
    // See the comments at the beginning of this file for discussion of
    // the zero exponent special case.

    copy_set_allele_freqs(marker_uidx, marker_exclude, multiplex, marker_idx, marker_ct, nullptr, set_allele_freqs, set_allele_freq_buf);
    if (missing_wt_needed) {
      uii = marker_ct - marker_idx;
      if (uii > multiplex) {
	uii = multiplex;
      }
      memcpy(wtbuf, &(dist_missing_wts_i[marker_idx]), uii * sizeof(int32_t));
    }
    retval = block_load(bedfile, bed_offset, marker_exclude, marker_ct, multiplex, unfiltered_sample_ct4, bedbuf, &marker_uidx, &marker_idx, &ujj);
    if (retval) {
      goto calc_distance_ret_1;
    }
    if (ujj < multiplex) {
      memset(&(bedbuf[ujj * unfiltered_sample_ct4]), 0, (multiplex - ujj) * unfiltered_sample_ct4);
      if (!main_weights) {
	fill_ulong_zero(sample_ct * (MULTIPLEX_2DIST / BITCT), geno);
	fill_ulong_zero(sample_ct * (MULTIPLEX_2DIST / BITCT), masks);
      } else {
	fill_ulong_zero(sample_ct, geno);
	fill_ulong_zero(sample_ct, masks);
      }
    }
    is_last_block = (marker_idx == marker_ct);
    if (!main_weights) {
      for (ukk = 0; ukk < ujj; ukk += BITCT) {
	glptr = &(geno[ukk / BITCT2]);
	glptr2 = &(masks[ukk / BITCT2]);
	glptr3 = mmasks;
	if (missing_wt_needed) {
	  giptr = sample_missing;
	}
	if (unwt_needed) {
	  giptr2 = sample_missing_unwt;
	}
	for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
	  if (!IS_SET(sample_exclude, sample_uidx)) {
	    unn = (sample_uidx % 4) * 2;
	    ulii = 0;
	    ulkk = 0;
	    gptr = &(bedbuf[sample_uidx / 4 + ukk * unfiltered_sample_ct4]);
	    for (umm = 0; umm < BITCT2; umm++) {
	      uljj = (gptr[umm * unfiltered_sample_ct4] >> unn) & 3;
	      ulii |= uljj << (umm * 2);
	      if (uljj == 1) {
		ulkk |= ONELU << umm;
		if (missing_wt_needed) {
		  *giptr += wtbuf[umm + ukk];
		}
		if (unwt_needed) {
		  *giptr2 += 1;
		}
	      }
	    }
	    // use xor to convert representation to 0 = missing,
	    // 1 or 2 = homozygote, 3 = heterozygote
	    ulii ^= FIVEMASK;
	    *glptr++ = ulii;
	    ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	    *glptr2++ = ulii * 3;
	    *glptr3 = ulkk;
	    ulii = 0;
	    ulkk = 0;
	    gptr = &(bedbuf[sample_uidx / 4 + (ukk + BITCT2) * unfiltered_sample_ct4]);
	    for (umm = 0; umm < BITCT2; umm++) {
	      uljj = (gptr[umm * unfiltered_sample_ct4] >> unn) & 3;
	      ulii |= uljj << (umm * 2);
	      if (uljj == 1) {
		ulkk |= ONELU << umm;
		if (missing_wt_needed) {
		  *giptr += wtbuf[umm + ukk + BITCT2];
		}
		if (unwt_needed) {
		  *giptr2 += 1;
		}
	      }
	    }
	    ulii ^= FIVEMASK;
	    *glptr = ulii;
	    ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	    *glptr2 = ulii * 3;
	    *glptr3++ |= ulkk << BITCT2;
	    glptr = &(glptr[(MULTIPLEX_2DIST / BITCT) - 1]);
	    glptr2 = &(glptr2[(MULTIPLEX_2DIST / BITCT) - 1]);
	    if (missing_wt_needed) {
	      giptr++;
	    }
	    if (unwt_needed) {
	      giptr2++;
	    }
	  }
	}

	if (missing_wt_needed) {
	  g_subset_weights_i = &(wtbuf[ukk]);
	}
	uii = is_last_block && (ukk + BITCT >= ujj);
	if (spawn_threads2(threads, &calc_ibs_thread, dist_thread_ct, uii)) {
	  goto calc_distance_ret_THREAD_CREATE_FAIL;
	}
	ulii = 0;
	calc_ibs_thread((void*)ulii);
	join_threads2(threads, dist_thread_ct, uii);
      }
    } else {
      fill_ulong_zero(sample_ct, mmasks);
      for (ukk = 0; ukk < ujj; ukk += MULTIPLEX_DIST_EXP / 2) {
	glptr = geno;
	glptr2 = masks;
	glptr3 = mmasks;
	giptr3 = sample_missing;
	for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
	  if (!IS_SET(sample_exclude, sample_uidx)) {
	    unn = (sample_uidx % 4) * 2;
	    ulii = 0;
	    ulkk = 0;
	    gptr = &(bedbuf[sample_uidx / 4 + ukk * unfiltered_sample_ct4]);
	    for (umm = 0; umm < MULTIPLEX_DIST_EXP / 2; umm++) {
	      uljj = (gptr[umm * unfiltered_sample_ct4] >> unn) & 3;
	      ulii |= uljj << (umm * 2);
	      if (uljj == 1) {
		ulkk |= ONELU << umm;
		*giptr3 += wtbuf[umm + ukk];
	      }
	    }
#ifdef __LP64__
	    ulii ^= FIVEMASK;
	    *glptr++ = ulii;
	    ulii = (ulii | (ulii >> 1)) & FIVEMASK;
#else
	    // note that FIVEMASK does NOT work here because the mask is
	    // only 28 bits
	    ulii ^= 0x05555555;
	    *glptr++ = ulii;
	    ulii = (ulii | (ulii >> 1)) & 0x05555555;
#endif
	    *glptr2++ = ulii * 3;
	    *glptr3++ |= ulkk << ukk;
	    giptr3++;
	  }
	}
	fill_subset_weights(subset_weights, &(main_weights[ukk]));
	uii = is_last_block && (ukk + (MULTIPLEX_DIST_EXP / 3) >= ujj);
	if (spawn_threads2(threads, &calc_wdist_thread, dist_thread_ct, uii)) {
	  goto calc_distance_ret_THREAD_CREATE_FAIL;
	}
	ulii = 0;
	calc_wdist_thread((void*)ulii);
	join_threads2(threads, dist_thread_ct, uii);
      }
    }
    printf("\r%" PRIuPTR " markers complete.", marker_idx);
    fflush(stdout);
  } while (!is_last_block);
  putc_unlocked('\r', stdout);
  logprint("Distance matrix calculation complete.\n");
  bigstack_reset(masks);
  if (calculation_type & (CALC_PLINK1_DISTANCE_MATRIX | CALC_PLINK1_IBS_MATRIX)) {
    if (bigstack_alloc_c(16 * sample_ct, &writebuf)) {
      goto calc_distance_ret_NOMEM;
    }
  }
  if (calculation_type & CALC_PLINK1_DISTANCE_MATRIX) {
    strcpy(outname_end, ".mdist");
    if (fopen_checked(outname, "w", &outfile)) {
      goto calc_distance_ret_OPEN_FAIL;
    }
    iptr = g_idists;
    giptr = g_missing_dbl_excluded;
    pct = 1;
    // parallel_tot must be 1 for --distance-matrix
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      giptr2 = sample_missing_unwt;
      uii = marker_ct - giptr2[sample_idx];
      wptr = writebuf;
      for (ujj = 0; ujj < sample_idx; ujj++) {
	wptr = dtoa_gx(((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++))), ' ', wptr);
      }
      wptr = memcpya(wptr, "0 ", 2);
      giptr2++;
      for (ulii = sample_idx + 1; ulii < sample_ct; ulii++) {
	uljj = tri_coord_no_diag(sample_idx, ulii);
	wptr = dtoa_gx(((double)g_idists[uljj]) / (2 * (uii - (*giptr2++) + g_missing_dbl_excluded[uljj])), ' ', wptr);
      }
      *wptr++ = '\n';
      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	goto calc_distance_ret_WRITE_FAIL;
      }
      if (sample_idx * 100LLU >= ((uint64_t)pct * sample_ct)) {
	pct = (sample_idx * 100LLU) / sample_ct;
	printf("\rWriting... %u%%", pct++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_distance_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    if (!parallel_idx) {
      wptr = strcpya(g_logbuf, "Distances (proportions) written to ");
      wptr = strcpya(wptr, outname);
      strcpy(outname_end, ".mdist.id");
      retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
      if (retval) {
	goto calc_distance_ret_1;
      }
      sprintf(wptr, " , and IDs to %s .\n", outname);
    } else {
      sprintf(g_logbuf, "Distances (proportions) written to %s .\n", outname);
    }
    wordwrapb(0);
    logprintb();
  }
  if (calculation_type & CALC_PLINK1_IBS_MATRIX) {
    strcpy(outname_end, ".mibs");
    if (fopen_checked(outname, "w", &outfile)) {
      goto calc_distance_ret_OPEN_FAIL;
    }
    iptr = g_idists;
    giptr = g_missing_dbl_excluded;
    pct = 1;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      giptr2 = sample_missing_unwt;
      uii = marker_ct - giptr2[sample_idx];
      wptr = writebuf;
      for (ujj = 0; ujj < sample_idx; ujj++) {
	wptr = dtoa_gx(1.0 - (((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++)))), ' ', wptr);
      }
      wptr = memcpya(wptr, "1 ", 2);
      giptr2++;
      for (ulii = sample_idx + 1; ulii < sample_ct; ulii++) {
	uljj = tri_coord_no_diag(sample_idx, ulii);
	wptr = dtoa_gx(1.0 - (((double)g_idists[uljj]) / (2 * (uii - (*giptr2++) + g_missing_dbl_excluded[uljj]))), ' ', wptr);
      }
      *wptr++ = '\n';
      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	goto calc_distance_ret_WRITE_FAIL;
      }
      if (sample_idx * 100 >= (pct * sample_ct)) {
	pct = (sample_idx * 100) / sample_ct;
	printf("\rWriting... %u%%", pct++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_distance_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
    if (retval) {
      goto calc_distance_ret_1;
    }
    outname_end[5] = '\0';
    LOGPRINTFWW("IBS matrix written to %s , and IDs to %s.id .\n", outname, outname);
  } while (!is_last_block);
  tstc = g_thread_start[dist_thread_ct];
  if (missing_wt_needed) {
    giptr = g_missing_tot_weights;
    dptr2 = g_dists;
    if (!main_weights) {
      iptr = g_idists;
      for (sample_idx = g_thread_start[0]; sample_idx < tstc; sample_idx++) {
	giptr2 = sample_missing;
	uii = giptr2[sample_idx];
	for (ujj = 0; ujj < sample_idx; ujj++) {
	  *dptr2++ = (marker_weight_sum_d / ((marker_weight_sum - uii - (*giptr2++)) + (*giptr++))) * (*iptr++);
	}
      }
    } else {
      for (sample_idx = g_thread_start[0]; sample_idx < tstc; sample_idx++) {
	giptr2 = sample_missing;
	uii = giptr2[sample_idx];
	for (ujj = 0; ujj < sample_idx; ujj++) {
	  *dptr2 *= (marker_weight_sum_d / ((marker_weight_sum - uii - (*giptr2++)) + (*giptr++)));
	  dptr2++;
	}
      }
    }
  } else if (dist_calc_type & DISTANCE_FLAT_MISSING) {
    dptr2 = g_dists;
    giptr = g_missing_dbl_excluded;
    if (!main_weights) {
      iptr = g_idists;
      if (dist_calc_type & DISTANCE_CLUSTER) {
	// save as IBS
        for (sample_idx = g_thread_start[0]; sample_idx < tstc; sample_idx++) {
	  giptr2 = sample_missing_unwt;
	  uii = marker_ct - giptr2[sample_idx];
	  for (ujj = 0; ujj < sample_idx; ujj++) {
	    *dptr2++ = 1.0 - (((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++))));
	  }
	}
      } else {
	for (sample_idx = g_thread_start[0]; sample_idx < tstc; sample_idx++) {
	  giptr2 = sample_missing_unwt;
	  uii = marker_ct - giptr2[sample_idx];
	  for (ujj = 0; ujj < sample_idx; ujj++) {
	    *dptr2++ = (((double)marker_ct) / (uii - (*giptr2++) + (*giptr++))) * (*iptr++);
	  }
	}
      }
    } else {
      for (sample_idx = g_thread_start[0]; sample_idx < tstc; sample_idx++) {
	giptr2 = sample_missing_unwt;
	uii = marker_ct - giptr2[sample_idx];
	for (ujj = 0; ujj < sample_idx; ujj++) {
	  *dptr2 *= ((double)marker_ct) / (uii - (*giptr2++) + (*giptr++));
	  dptr2++;
	}
      }
    }
  }

  if (calculation_type & (CALC_DISTANCE | CALC_IBS_TEST)) {
    if ((distance_exp == 0.0) || (!(dist_calc_type & (DISTANCE_IBS | DISTANCE_1_MINUS_IBS)))) {
      g_half_marker_ct_recip = 0.5 / (double)marker_ct;
    } else {
      dyy = 0.0;
      marker_uidx = 0;
      for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	dxx = set_allele_freqs[marker_uidx];
	if ((dxx > 0.0) && (dxx < 1.0)) {
	  dyy += pow(2 * dxx * (1.0 - dxx), -distance_exp);
	} else {
	  dyy += 1.0;
	}
      }
      g_half_marker_ct_recip = 0.5 / dyy;
    }
    if (calculation_type & CALC_DISTANCE) {
      if (!parallel_idx) {
	retval = distance_d_write_ids(outname, outname_end, dist_calc_type, unfiltered_sample_ct, sample_exclude, sample_ids, max_sample_id_len);
	if (retval) {
	  goto calc_distance_ret_1;
	}
        LOGPRINTFWW("IDs written to %s .\n", outname);
      }
      retval = distance_d_write(&outfile, &outfile2, &outfile3, dist_calc_type, outname, outname_end, g_dists, g_half_marker_ct_recip, sample_ct, g_thread_start[0], g_thread_start[dist_thread_ct], parallel_idx, parallel_tot, (unsigned char*)geno);
      if (retval) {
	goto calc_distance_ret_1;
      }
    }
  }
  bigstack_reset(bigstack_mark);

  while (0) {
  calc_distance_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_distance_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_distance_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_distance_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_distance_ret_1:
  fclose_cond(outfile);
  fclose_cond(outfile2);
  fclose_cond(outfile3);
  return retval;
}

void cluster_dist_divide(uintptr_t sample_ct, uintptr_t cluster_ct, uint32_t* cluster_starts, double* cluster_sdistances) {
  uintptr_t tcoord;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  double dxx;
  for (ulii = 0; ulii < cluster_ct; ulii++) {
    uii = cluster_starts[ulii + 1] - cluster_starts[ulii];
    if (uii > 1) {
      dxx = 1.0 / ((double)((int32_t)uii));
      uljj = (ulii * (ulii + 1)) / 2;
      for (tcoord = (ulii * (ulii - 1)) / 2; tcoord < uljj; tcoord++) {
	cluster_sdistances[tcoord] *= dxx;
      }
      for (uljj = ulii + 1; uljj < sample_ct; uljj++) {
	cluster_sdistances[tri_coord_no_diag(ulii, uljj)] *= dxx;
      }
    }
  }
}

void cluster_dist_multiply(uintptr_t sample_ct, uintptr_t cluster_ct, uint32_t* cluster_starts, double* cluster_sdistances) {
  uintptr_t tcoord;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  double dxx;
  for (ulii = 0; ulii < cluster_ct; ulii++) {
    uii = cluster_starts[ulii + 1] - cluster_starts[ulii];
    if (uii > 1) {
      dxx = ((double)((int32_t)uii));
      uljj = (ulii * (ulii + 1)) / 2;
      for (tcoord = (ulii * (ulii - 1)) / 2; tcoord < uljj; tcoord++) {
	cluster_sdistances[tcoord] *= dxx;
      }
      for (uljj = ulii + 1; uljj < sample_ct; uljj++) {
	cluster_sdistances[tri_coord_no_diag(ulii, uljj)] *= dxx;
      }
    }
  }
}

int32_t calc_cluster_neighbor(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, char* read_dists_fname, char* read_dists_id_fname, char* read_genome_fname, char* outname, char* outname_end, uint64_t calculation_type, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Cluster_info* cp, int32_t missing_pheno, uint32_t neighbor_n1, uint32_t neighbor_n2, uint32_t ppc_gap, uintptr_t* pheno_c, double* mds_plot_dmatrix_copy, uintptr_t* cluster_merge_prevented, double* cluster_sorted_ibs, unsigned char* bigstack_mark_precluster, unsigned char* bigstack_mark_postcluster) {
  // --cluster and --neighbour.  They are handled by the same function because
  // they initially process the distance matrix/PPC test results in roughly the
  // same way, but we have removed the PLINK 1.07 requirement that --cluster be
  // invoked in every --neighbour run.
  //
  // This sorted list/heap --cluster implementation is O(n^2 log n) time,
  // O(n^2) space.  While this is a major improvement over PLINK 1.07's O(n^3)
  // time, it is possible to do even better in some scenarios.  In particular,
  // if there are no clustering restrictions, Defays' CLINK algorithm (see e.g.
  // http://comjnl.oxfordjournals.org/content/20/4/364.full.pdf ) requires only
  // O(n) space, so it's an excellent complement to --distance/--genome +
  // --parallel on very large datasets (e.g. 500k samples), and should be added
  // as a special case in the future.
  FILE* outfile = nullptr;
  uint32_t* cluster_sorted_ibs_indices = nullptr;
#ifdef __LP64__
  // uint64_t* cluster_sorted_ibs_indices_big = nullptr;
#endif
  uint32_t* sample_to_cluster = nullptr;
  double* neighbor_quantiles = nullptr;
  double* neighbor_quantile_means = nullptr;
  double* neighbor_quantile_stdev_recips = nullptr;
  uint32_t* neighbor_qindices = nullptr;
  uint32_t* ppc_fail_counts = nullptr;
  uint32_t* cur_cluster_sizes = nullptr;
  uint32_t* cur_cluster_case_cts = nullptr;
  uint32_t* cur_cluster_remap = nullptr;
  uint32_t* cluster_index = nullptr;
  uintptr_t* collapsed_pheno_c = nullptr;
  uint32_t* sample_idx_to_uidx = nullptr;
  uint32_t* late_clidx_to_sample_uidx = nullptr;
  uintptr_t* ibs_ties = nullptr;
  uint32_t* genome_main = g_genome_main;
  double* dptr = nullptr;
  double min_ppc = cp->ppc;
  double min_ibm = cp->min_ibm;
  double min_zx = 0.0;
  uint32_t ibm_constraint = (min_ibm != 0.0);
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  double sample_ct_recip = 1.0 / ((double)((intptr_t)sample_ct));
  uint32_t do_neighbor = (calculation_type / CALC_NEIGHBOR) & 1;
  uint32_t is_group_avg = (cp->modifier / CLUSTER_GROUP_AVG) & 1;
  uint32_t is_mds_cluster = cp->modifier & CLUSTER_MDS;
  uint32_t is_old_tiebreaks = cp->modifier & CLUSTER_OLD_TIEBREAKS;
  uint32_t mds_fill_nonclust = (!is_mds_cluster) && mds_plot_dmatrix_copy;
  uint32_t cluster_missing = cp->modifier & CLUSTER_MISSING;
  uint32_t use_genome_dists = (!g_dists) && (!read_dists_fname) && (!cluster_missing);
  uint32_t cur_cluster_ct = sample_ct;
  uintptr_t neighbor_row_ct = neighbor_n2 - neighbor_n1 + 1;
  uint32_t mc_warning = 0;
  uint32_t mcc_warning = 0;
  uint32_t ppc_warning = 0;
  uint32_t ibm_warning = 0;
  uintptr_t tcoord = 0;
  uintptr_t clidx1 = 0;
  uintptr_t clidx2 = 0;
  uint32_t case_ct = 0xffffffffU;
  uint32_t ctrl_ct = 0xffffffffU;
  double dyy = 0.0;
  uint32_t pct = 1;
  int32_t retval = 0;
  uintptr_t initial_triangle_size;
  uintptr_t heap_size;
  double* dptr2;
  uint32_t* sample_missing_unwt;
  uint32_t* missing_dbl_excluded;
  uint32_t* sample_missing_ptr;
  uint32_t* dbl_exclude_ptr;
  uint32_t* merge_sequence;
  uint32_t* uiptr;
  uint32_t* uiptr2;
  uint32_t merge_ct;
  char* wptr_start;
  char* fam_id;
  char* sample_id;
  char* wptr;
  double dxx;
  double dxx1;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t sample_idx1;
  uintptr_t sample_idx2;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  g_sample_ct = sample_ct;

  if (sample_ct <= cp->min_ct) {
    // todo: check if anyone uses --mds-plot without caring about clustering;
    // if yes, we should stop forcing --mds-plot to be used with --cluster.
    logerrprint("Error: --K parameter too large (>= sample count).\n");
    goto calc_cluster_neighbor_ret_INVALID_CMDLINE;
  }
  if (cluster_ct) {
    cur_cluster_ct += cluster_ct - cluster_starts[cluster_ct];
    if (cur_cluster_ct <= cp->min_ct) {
      logerrprint("Error: --K parameter too large (>= initial cluster count).\n");
      goto calc_cluster_neighbor_ret_INVALID_CMDLINE;
    }
    sample_to_cluster = (uint32_t*)malloc(sample_ct * sizeof(int32_t));
    if (!sample_to_cluster) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    if (cluster_starts[cluster_ct] < sample_ct) {
      late_clidx_to_sample_uidx = (uint32_t*)malloc((sample_ct - cluster_starts[cluster_ct]) * sizeof(int32_t));
      if (!late_clidx_to_sample_uidx) {
        goto calc_cluster_neighbor_ret_NOMEM;
      }
    }
    retval = fill_sample_to_cluster(unfiltered_sample_ct, sample_exclude, sample_ct, cluster_ct, cluster_map, cluster_starts, sample_to_cluster, late_clidx_to_sample_uidx);
    if (retval) {
      goto calc_cluster_neighbor_ret_1;
    }
  }
  initial_triangle_size = (((uintptr_t)cur_cluster_ct) * (cur_cluster_ct - 1)) >> 1;
  heap_size = initial_triangle_size;
  if (!cluster_sorted_ibs) {
    // --genome calculation not run, no --within
    if (g_dists) {
      cluster_sorted_ibs = g_dists;
    } else {
      if (bigstack_calloc_d(initial_triangle_size, &cluster_sorted_ibs)) {
        goto calc_cluster_neighbor_ret_NOMEM;
      }
    }
    bigstack_mark_precluster = g_bigstack_base;
  }
  if (do_neighbor) {
    if (neighbor_n2 >= sample_ct) {
      logerrprint("Error: Second --neighbour parameter too large (>= population size).\n");
      goto calc_cluster_neighbor_ret_INVALID_CMDLINE;
    }
    ulii = neighbor_n2 * sample_ct;
    neighbor_quantiles = (double*)malloc(ulii * sizeof(double));
    if (!neighbor_quantiles) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    neighbor_qindices = (uint32_t*)malloc(ulii * sizeof(int32_t));
    if (!neighbor_qindices) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    neighbor_quantile_means = (double*)malloc(neighbor_row_ct * sizeof(double));
    if (!neighbor_quantile_means) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    neighbor_quantile_stdev_recips = (double*)malloc(neighbor_row_ct * sizeof(double));
    if (!neighbor_quantile_stdev_recips) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    fill_double_zero(ulii, neighbor_quantiles);
  }
  fill_ulong_zero(BITCT_TO_WORDCT(initial_triangle_size), cluster_merge_prevented);
  if ((min_ppc != 0.0) || genome_main || read_genome_fname) {
    if (do_neighbor && (min_ppc != 0.0)) {
      ppc_fail_counts = (uint32_t*)malloc(sample_ct * sizeof(int32_t));
      if (!ppc_fail_counts) {
        goto calc_cluster_neighbor_ret_NOMEM;
      }
      fill_uint_zero(sample_ct, ppc_fail_counts);
    }
    if (read_genome_fname) {
      retval = read_genome(read_genome_fname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, cluster_merge_prevented, use_genome_dists? cluster_sorted_ibs : nullptr, neighbor_n2, neighbor_quantiles, neighbor_qindices, ppc_fail_counts, min_ppc, cluster_ct && (!is_group_avg), cluster_ct, cluster_starts, sample_to_cluster);
      if (retval) {
	goto calc_cluster_neighbor_ret_1;
      }
    } else {
      // N.B. the first --genome "row" is actually the longest, while the
      // reverse is the case for every other internal distance matrix.
      // This changes coordinate mapping for g_missing_dbl_excluded too!

      // premultiply Z threshold to remove a multiply from the inner loop
      if (min_ppc != 0.0) {
        min_zx = ltqnorm(min_ppc) * sqrt(0.2222222);
      }
      ulii = 0;
      dbl_exclude_ptr = g_missing_dbl_excluded;
      sample_missing_unwt = g_sample_missing_unwt;
      for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
        sample_missing_ptr = &(sample_missing_unwt[sample_idx1]);
	uii = marker_ct - (*sample_missing_ptr++);
	for (sample_idx2 = sample_idx1 + 1; sample_idx2 < sample_ct; sample_idx2++) {
	  if (cluster_ct) {
            clidx1 = sample_to_cluster[sample_idx1];
            clidx2 = sample_to_cluster[sample_idx2];
            if (clidx1 < clidx2) {
	      tcoord = tri_coord_no_diag(clidx1, clidx2);
	    } else {
              tcoord = tri_coord_no_diag(clidx2, clidx1);
	    }
	  }
          if (use_genome_dists) {
	    dxx = 1.0 - (((double)(genome_main[ulii] + 2 * genome_main[ulii + 1])) / ((double)(2 * (uii - (*sample_missing_ptr++) + (*dbl_exclude_ptr++)))));
	    if (do_neighbor) {
	      update_neighbor(sample_ct, neighbor_n2, sample_idx1, sample_idx2, dxx, neighbor_quantiles, neighbor_qindices);
	    }
	    if (cluster_ct) {
	      if (clidx1 != clidx2) {
		if (!is_group_avg) {
		  if (dxx < cluster_sorted_ibs[tcoord]) {
                    cluster_sorted_ibs[tcoord] = dxx;
		  }
		} else {
                  cluster_sorted_ibs[tcoord] += dxx;
		}
	      }
	    } else {
	      cluster_sorted_ibs[tri_coord_no_diag(sample_idx1, sample_idx2)] = dxx;
	    }
	    if (mds_plot_dmatrix_copy) {
              mds_plot_dmatrix_copy[tri_coord_no_diag(sample_idx1, sample_idx2)] = dxx;
	    }
	  }
	  if (min_ppc != 0.0) {
	    dxx = (double)genome_main[ulii + 4];
	    dxx1 = 1.0 / ((double)(genome_main[ulii + 4] + genome_main[ulii + 3]));
	    if (((dxx * dxx1 - 0.666666) / sqrt(dxx1)) < min_zx) {
	      if (do_neighbor) {
		ppc_fail_counts[sample_idx1] += 1;
		ppc_fail_counts[sample_idx2] += 1;
	      }
	      if (!cluster_ct) {
		set_bit_ul(tri_coord_no_diag(sample_idx1, sample_idx2), cluster_merge_prevented);
	      } else {
		if (clidx1 != clidx2) {
		  SET_BIT(tcoord, cluster_merge_prevented);
		} else if (!ppc_warning) {
		  logerrprint("Warning: Initial cluster assignment violates PPC test constraint.\n");
		  ppc_warning = 1;
		}
	      }
	    }
	  }
	  ulii += 5;
	}
      }
      bigstack_reset(genome_main);
    }
  } else if (((!cluster_missing) || do_neighbor) && (!read_dists_fname)) {
    // calculate entire distance matrix, or use already-calculated matrix in
    // memory
    if (!g_dists) {
      retval = calc_distance(threads, 0, 1, bedfile, bed_offset, outname, outname_end, nullptr, nullptr, 0.0, 0, DISTANCE_FLAT_MISSING | DISTANCE_CLUSTER, unfiltered_marker_ct, marker_exclude, marker_ct, nullptr, 0, set_allele_freqs, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, chrom_info_ptr);
      if (retval) {
        goto calc_cluster_neighbor_ret_1;
      }
      if (!cluster_ct) {
        cluster_sorted_ibs = g_dists;
      }
    }
    if (do_neighbor) {
      dptr = g_dists;
      for (sample_idx1 = 1; sample_idx1 < sample_ct; sample_idx1++) {
	for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	  update_neighbor(sample_ct, neighbor_n2, sample_idx1, sample_idx2, *dptr++, neighbor_quantiles, neighbor_qindices);
	}
      }
    }
    if (mds_fill_nonclust) {
      memcpy(mds_plot_dmatrix_copy, g_dists, (sample_ct * (sample_ct - 1)) * (sizeof(double) / 2));
    }
    if (cluster_ct) {
      dptr = g_dists;
      if (!is_group_avg) {
	for (sample_idx1 = 1; sample_idx1 < sample_ct; sample_idx1++) {
	  clidx1 = sample_to_cluster[sample_idx1];
	  for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	    clidx2 = sample_to_cluster[sample_idx2];
	    if (clidx1 != clidx2) {
	      dxx = *dptr++;
	      if (clidx1 < clidx2) {
	        dptr2 = &(cluster_sorted_ibs[tri_coord_no_diag(clidx1, clidx2)]);
	      } else {
	        dptr2 = &(cluster_sorted_ibs[tri_coord_no_diag(clidx2, clidx1)]);
	      }
	      if (dxx < (*dptr2)) {
		*dptr2 = dxx;
	      }
	    } else {
	      dptr++;
	    }
	  }
	}
      } else {
	for (sample_idx1 = 1; sample_idx1 < sample_ct; sample_idx1++) {
	  clidx1 = sample_to_cluster[sample_idx1];
	  for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	    clidx2 = sample_to_cluster[sample_idx2];
	    if (clidx1 < clidx2) {
	      cluster_sorted_ibs[tri_coord_no_diag(clidx1, clidx2)] += *dptr++;
	    } else if (clidx1 > clidx2) {
	      cluster_sorted_ibs[tri_coord_no_diag(clidx2, clidx1)] += *dptr++;
	    } else {
	      dptr++;
	    }
	  }
	}
      }
      bigstack_reset(g_dists);
    }
  }
  if (read_dists_fname) {
    retval = read_dists(read_dists_fname, read_dists_id_fname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, cluster_ct, cluster_starts, sample_to_cluster, (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))? 1 : 0, cluster_ct && (!is_group_avg), cluster_sorted_ibs, neighbor_n2, neighbor_quantiles, neighbor_qindices);
    if (retval) {
      goto calc_cluster_neighbor_ret_1;
    }
  }
  sample_idx_to_uidx = (uint32_t*)malloc(sample_ct * sizeof(int32_t));
  if (!sample_idx_to_uidx) {
    goto calc_cluster_neighbor_ret_NOMEM;
  }
  fill_idx_to_uidx(sample_exclude, unfiltered_sample_ct, sample_ct, sample_idx_to_uidx);
  if (do_neighbor) {
    memcpy(outname_end, ".nearest", 9);
    if (fopen_checked(outname, "w", &outfile)) {
      goto calc_cluster_neighbor_ret_OPEN_FAIL;
    }
    if (fputs_checked("         FID          IID     NN      MIN_DST            Z         FID2         IID2 ", outfile)) {
      goto calc_cluster_neighbor_ret_WRITE_FAIL;
    }
    dptr = &(neighbor_quantiles[(neighbor_n1 - 1) * sample_ct]);
    for (sample_idx1 = neighbor_n1 - 1; sample_idx1 < neighbor_n2; sample_idx1++) {
      dxx = 0.0; // sum
      dxx1 = 0.0; // ssq
      for (sample_idx2 = 0; sample_idx2 < sample_ct; sample_idx2++) {
        dyy = *dptr++;
        dxx += dyy;
        dxx1 += dyy * dyy;
      }
      dyy = dxx * sample_ct_recip;
      neighbor_quantile_means[sample_idx1 + 1 - neighbor_n1] = dyy;
      neighbor_quantile_stdev_recips[sample_idx1 + 1 - neighbor_n1] = sqrt(((double)((intptr_t)(sample_ct - 1))) / (dxx1 - dxx * dyy));
    }
    if (min_ppc != 0.0) {
      fputs("   PROP_DIFF ", outfile);
    }
    putc_unlocked('\n', outfile);
    dxx1 = 1.0 / ((double)((intptr_t)(sample_ct - 1)));
    for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
      fam_id = &(sample_ids[sample_idx_to_uidx[sample_idx1] * max_sample_id_len]);
      sample_id = (char*)memchr(fam_id, '\t', max_sample_id_len);
      wptr_start = fw_strcpyn(12, (uint32_t)(sample_id - fam_id), fam_id, g_textbuf);
      *wptr_start++ = ' ';
      wptr_start = fw_strcpy(12, &(sample_id[1]), wptr_start);
      *wptr_start++ = ' ';
      if (min_ppc != 0.0) {
        dyy = ((double)((int32_t)ppc_fail_counts[sample_idx1])) * dxx1;
      }
      for (ulii = 0; ulii < neighbor_row_ct; ulii++) {
        wptr = uint32toa_w6x(ulii + neighbor_n1, ' ', wptr_start);
	sample_idx2 = sample_idx1 + ulii * sample_ct;
        dxx = neighbor_quantiles[sample_idx2];
	wptr = dtoa_g_wxp4x(dxx, 12, ' ', wptr);
        wptr = dtoa_g_wxp4x((dxx - neighbor_quantile_means[ulii]) * neighbor_quantile_stdev_recips[ulii], 12, ' ', wptr);
	sample_idx2 = neighbor_qindices[sample_idx2];
        fam_id = &(sample_ids[sample_idx_to_uidx[sample_idx2] * max_sample_id_len]);
        sample_id = (char*)memchr(fam_id, '\t', max_sample_id_len);
        wptr = fw_strcpyn(12, (uint32_t)(sample_id - fam_id), fam_id, wptr);
        *wptr++ = ' ';
        wptr = fw_strcpy(12, &(sample_id[1]), wptr);
        *wptr++ = ' ';
	if (min_ppc != 0.0) {
          wptr = dtoa_g_wxp4x(dyy, 12, ' ', wptr);
	}
        *wptr++ = '\n';
        if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto calc_cluster_neighbor_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_cluster_neighbor_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--neighbour report written to %s .\n", outname);
    if (!(calculation_type & CALC_CLUSTER)) {
      goto calc_cluster_neighbor_ret_1;
    }
  }
  if (cluster_missing || ibm_constraint) {
    if (!g_missing_dbl_excluded) {
      retval = calc_ibm(threads, bedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_ct, unfiltered_sample_ct, sample_exclude, sample_ct, chrom_info_ptr);
      if (retval) {
        goto calc_cluster_neighbor_ret_1;
      }
      if (!cluster_missing) {
        logprint("IBM matrix calculated.");
        fputs("      ", stdout);
        logprint("\n");
      }
      // N.B. this is still used until the end of the block
      bigstack_reset(g_missing_dbl_excluded);
    }
    dxx1 = 1.0 / ((double)((int32_t)marker_ct));
    if (cluster_missing) {
      memcpy(outname_end, ".mdist.missing", 15);
      if (fopen_checked(outname, "w", &outfile)) {
	goto calc_cluster_neighbor_ret_OPEN_FAIL;
      }
      fputs("Writing IBM matrix... 0%    \b\b\b\b", stdout);
      fflush(stdout);
    }
    sample_missing_unwt = g_sample_missing_unwt;
    missing_dbl_excluded = g_missing_dbl_excluded;
    for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
      sample_missing_ptr = sample_missing_unwt;
      uii = sample_missing_ptr[sample_idx1];
      if (!cluster_ct) {
	dptr = &(cluster_sorted_ibs[(sample_idx1 * (sample_idx1 - 1)) >> 1]);
      } else if (mds_fill_nonclust) {
	dptr = &(mds_plot_dmatrix_copy[(sample_idx1 * (sample_idx1 - 1)) >> 1]);
      }
      ulii = sample_ct - 2;
      if (!genome_main) {
	dbl_exclude_ptr = &(missing_dbl_excluded[(sample_idx1 * (sample_idx1 - 1)) >> 1]);
	if (!cluster_ct) {
	  for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	    dxx = 1.0 - ((double)((int32_t)(uii + (*sample_missing_ptr++) - 2 * (*dbl_exclude_ptr++)))) * dxx1;
	    if (cluster_missing) {
	      *dptr++ = dxx;
	      wptr = dtoa_gx(dxx, ' ', g_textbuf);
	      fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	    }
	    if (dxx < min_ibm) {
	      set_bit_ul(tri_coord_no_diag(sample_idx2, sample_idx1), cluster_merge_prevented);
	    }
	  }
	} else {
	  clidx1 = sample_to_cluster[sample_idx1];
	  for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	    dxx = 1.0 - ((double)((int32_t)(uii + (*sample_missing_ptr++) - 2 * (*dbl_exclude_ptr++)))) * dxx1;
	    clidx2 = sample_to_cluster[sample_idx2];
	    if (cluster_missing) {
	      if (mds_fill_nonclust) {
		*dptr++ = dxx;
	      }
	      if (clidx1 != clidx2) {
		if (clidx1 < clidx2) {
		  dptr2 = &(cluster_sorted_ibs[tri_coord_no_diag(clidx1, clidx2)]);
		} else {
		  dptr2 = &(cluster_sorted_ibs[tri_coord_no_diag(clidx2, clidx1)]);
		}
		if (!is_group_avg) {
		  if (dxx < (*dptr2)) {
		    *dptr2 = dxx;
		  }
		} else {
		  *dptr2 += dxx;
		}
	      }
	      wptr = dtoa_gx(dxx, ' ', g_textbuf);
	      fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	    }
	    if (dxx < min_ibm) {
	      if (clidx1 < clidx2) {
		set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	      } else if (clidx1 > clidx2) {
		set_bit_ul(tri_coord_no_diag(clidx2, clidx1), cluster_merge_prevented);
	      } else {
		ibm_warning = 1;
	      }
	    }
	  }
	}
      } else {
	dbl_exclude_ptr = &(missing_dbl_excluded[sample_idx1 - 1]);
	if (!cluster_ct) {
	  for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	    dxx = 1.0 - ((double)((int32_t)(uii + (*sample_missing_ptr++) - 2 * (*dbl_exclude_ptr)))) * dxx1;
	    dbl_exclude_ptr = &(dbl_exclude_ptr[ulii - sample_idx2]);
	    if (cluster_missing) {
	      *dptr++ = dxx;
	      wptr = dtoa_gx(dxx, ' ', g_textbuf);
	      fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	    }
	    if (dxx < min_ibm) {
	      set_bit_ul(tri_coord_no_diag(sample_idx2, sample_idx1), cluster_merge_prevented);
	    }
	  }
	} else {
	  clidx1 = sample_to_cluster[sample_idx1];
	  for (sample_idx2 = 0; sample_idx2 < sample_idx1; sample_idx2++) {
	    dxx = 1.0 - ((double)((int32_t)(uii + (*sample_missing_ptr++) - 2 * (*dbl_exclude_ptr)))) * dxx1;
	    dbl_exclude_ptr = &(dbl_exclude_ptr[ulii - sample_idx2]);
	    clidx2 = sample_to_cluster[sample_idx2];
	    if (cluster_missing) {
	      if (mds_fill_nonclust) {
		*dptr++ = dxx;
	      }
	      if (clidx1 != clidx2) {
		if (clidx1 < clidx2) {
		  dptr2 = &(cluster_sorted_ibs[tri_coord_no_diag(clidx1, clidx2)]);
		} else {
		  dptr2 = &(cluster_sorted_ibs[tri_coord_no_diag(clidx2, clidx1)]);
		}
		if (!is_group_avg) {
		  if (dxx < (*dptr2)) {
		    *dptr2 = dxx;
		  }
		} else {
		  *dptr2 += dxx;
		}
	      }
	      wptr = dtoa_gx(dxx, ' ', g_textbuf);
	      fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	    }
	    if (dxx < min_ibm) {
	      if (clidx1 < clidx2) {
		set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	      } else if (clidx1 > clidx2) {
		set_bit_ul(tri_coord_no_diag(clidx2, clidx1), cluster_merge_prevented);
	      } else {
		ibm_warning = 1;
	      }
	    }
	  }
	}
      }
      if (cluster_missing) {
	if (putc_checked('1', outfile)) {
	  goto calc_cluster_neighbor_ret_WRITE_FAIL;
	}
	putc_unlocked(' ', outfile);
	if (!genome_main) {
	  for (sample_idx2 = sample_idx1 + 1; sample_idx2 < sample_ct; sample_idx2++) {
	    dxx = 1.0 - ((double)((int32_t)(uii + (*(++sample_missing_ptr)) - 2 * missing_dbl_excluded[((sample_idx2 * (sample_idx2 - 1)) >> 1) + sample_idx1]))) * dxx1;
	    wptr = dtoa_gx(dxx, ' ', g_textbuf);
	    fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	  }
	} else {
	  // f(0) = 0
	  // f(1) = sample_ct - 1
	  // f(2) = 2 * sample_ct - 3
	  // ...
          // f(n) = n * sample_ct - n(n+1)/2
	  dbl_exclude_ptr = &(missing_dbl_excluded[sample_ct * sample_idx1 - ((sample_idx1 * (sample_idx1 + 1)) >> 1)]);
	  for (sample_idx2 = sample_idx1 + 1; sample_idx2 < sample_ct; sample_idx2++) {
	    dxx = 1.0 - ((double)((int32_t)(uii + (*(++sample_missing_ptr)) - 2 * (*dbl_exclude_ptr++)))) * dxx1;
	    wptr = dtoa_gx(dxx, ' ', g_textbuf);
	    fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	  }
	}
	if (putc_checked('\n', outfile)) {
	  goto calc_cluster_neighbor_ret_WRITE_FAIL;
	}
        if ((sample_idx1 + 1) * 100LLU >= sample_ct * ((uint64_t)pct)) {
          if (pct > 10) {
	    putc_unlocked('\b', stdout);
	  }
	  pct = ((sample_idx1 + 1) * 100LLU) / sample_ct;
          printf("\b\b%u%%", pct++);
	  fflush(stdout);
	}
      }
    }
    if (cluster_missing) {
      if (fclose_null(&outfile)) {
	goto calc_cluster_neighbor_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      LOGPRINTFWW("IBM matrix written to %s .\n", outname);
      if (ibm_warning) {
	logerrprint("Warning: Initial cluster assignment violates IBM constraint.\n");
      }
      if (mds_fill_nonclust && (!cluster_ct)) {
	memcpy(mds_plot_dmatrix_copy, cluster_sorted_ibs, initial_triangle_size * sizeof(double));
      }
    }
  }
  if (is_group_avg) {
    cluster_dist_divide(sample_ct, cluster_ct, cluster_starts, cluster_sorted_ibs);
  }
  if (mds_plot_dmatrix_copy && is_mds_cluster) {
    memcpy(mds_plot_dmatrix_copy, cluster_sorted_ibs, initial_triangle_size * sizeof(double));
  }
  if (pheno_c) {
    case_ct = popcount_longs(pheno_c, unfiltered_sample_ctl);
    ctrl_ct = sample_ct - case_ct;
    if ((cp->modifier & CLUSTER_CC) || (cp->max_cases < case_ct) || (cp->max_ctrls < ctrl_ct)) {
      if (bigstack_alloc_ul(sample_ctl, &collapsed_pheno_c)) {
	goto calc_cluster_neighbor_ret_NOMEM;
      }
      cur_cluster_case_cts = (uint32_t*)malloc(cur_cluster_ct * sizeof(int32_t));
      if (!cur_cluster_case_cts) {
	goto calc_cluster_neighbor_ret_NOMEM;
      }
      copy_bitarr_subset_excl(pheno_c, sample_exclude, unfiltered_sample_ct, sample_ct, collapsed_pheno_c);
      fill_uint_zero(cur_cluster_ct, cur_cluster_case_cts);
      if (!cluster_ct) {
	for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
	  if (IS_SET(collapsed_pheno_c, sample_idx1)) {
	    cur_cluster_case_cts[sample_idx1] += 1;
	  }
	}
      } else {
	for (sample_idx1 = 0; sample_idx1 < sample_ct; sample_idx1++) {
	  if (IS_SET(collapsed_pheno_c, sample_idx1)) {
	    cur_cluster_case_cts[sample_to_cluster[sample_idx1]] += 1;
	  }
	}
      }
      bigstack_reset(collapsed_pheno_c);
    }
  }
  if (cur_cluster_case_cts || is_group_avg || (cp->max_size < sample_ct)) {
    cur_cluster_sizes = (uint32_t*)malloc(cur_cluster_ct * sizeof(int32_t));
    if (!cur_cluster_sizes) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    uii = cp->max_size;
    for (clidx1 = 0; clidx1 < cluster_ct; clidx1++) {
      ujj = cluster_starts[clidx1 + 1] - cluster_starts[clidx1];
      cur_cluster_sizes[clidx1] = ujj;
      if (ujj >= uii) {
	if (ujj != uii) {
	  mc_warning = 1;
	}
	fill_bits((clidx1 * (clidx1 - 1)) >> 1, clidx1, cluster_merge_prevented);
	for (clidx2 = clidx1 + 1; clidx2 < cur_cluster_ct; clidx2++) {
	  set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	}
      } else if (ujj > uii / 2) {
	ujj = uii - ujj;
	for (clidx2 = 0; clidx2 < clidx1; clidx2++) {
	  if (cluster_starts[clidx2 + 1] - cluster_starts[clidx2] > ujj) {
	    set_bit_ul(tri_coord_no_diag(clidx2, clidx1), cluster_merge_prevented);
	  }
	}
	for (clidx2 = clidx1 + 1; clidx2 < cluster_ct; clidx2++) {
	  if (cluster_starts[clidx2 + 1] - cluster_starts[clidx2] > ujj) {
	    set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	  }
	}
      }
    }
    if (mc_warning) {
      logerrprint("Warning: Initial cluster assignment violates --mc restriction.\n");
    }
    for (clidx1 = cluster_ct; clidx1 < cur_cluster_ct; clidx1++) {
      cur_cluster_sizes[clidx1] = 1;
    }
    if ((cp->max_cases < case_ct) || (cp->max_ctrls < ctrl_ct)) {
      uii = cp->max_cases;
      if (uii < case_ct) {
	for (clidx1 = 0; clidx1 < cluster_ct; clidx1++) {
	  ujj = cur_cluster_case_cts[clidx1];
	  if (ujj > uii) {
	    mcc_warning = 1;
	    fill_bits((clidx1 * (clidx1 - 1)) >> 1, clidx1, cluster_merge_prevented);
	    for (clidx2 = clidx1 + 1; clidx2 < cur_cluster_ct; clidx2++) {
	      set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	    }
	  } else if (ujj > uii / 2) {
	    ujj = uii - ujj;
	    for (clidx2 = 0; clidx2 < clidx1; clidx2++) {
	      if (cur_cluster_case_cts[clidx2] > ujj) {
		set_bit_ul(tri_coord_no_diag(clidx2, clidx1), cluster_merge_prevented);
	      }
	    }
	    for (clidx2 = clidx1 + 1; clidx2 < cluster_ct; clidx2++) {
	      if (cur_cluster_case_cts[clidx2] > ujj) {
		set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	      }
	    }
	    if (!ujj) {
	      for (clidx2 = cluster_ct; clidx2 < cur_cluster_ct; clidx2++) {
                if (cur_cluster_case_cts[clidx2]) {
		  set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
		}
	      }
	    }
	  }
	}
	if (uii == 1) {
	  for (clidx1 = cluster_ct; clidx1 < cur_cluster_ct; clidx1++) {
	    if (cur_cluster_case_cts[clidx1]) {
	      for (clidx2 = clidx1 + 1; clidx2 < cur_cluster_ct; clidx2++) {
		if (cur_cluster_case_cts[clidx2]) {
		  set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
		}
	      }
	    }
	  }
	}
      }
      uii = cp->max_ctrls;
      if (uii < ctrl_ct) {
	for (clidx1 = 0; clidx1 < cluster_ct; clidx1++) {
	  ujj = cur_cluster_sizes[clidx1] - cur_cluster_case_cts[clidx1];
	  if (ujj > uii) {
	    mcc_warning = 1;
	    fill_bits((clidx1 * (clidx1 - 1)) >> 1, clidx1, cluster_merge_prevented);
	    for (clidx2 = clidx1 + 1; clidx2 < cur_cluster_ct; clidx2++) {
	      set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	    }
	  } else if (ujj > uii / 2) {
	    ujj = uii - ujj;
	    for (clidx2 = 0; clidx2 < clidx1; clidx2++) {
	      if (cur_cluster_sizes[clidx2] - cur_cluster_case_cts[clidx2] > ujj) {
		set_bit_ul(tri_coord_no_diag(clidx2, clidx1), cluster_merge_prevented);
	      }
	    }
	    for (clidx2 = clidx1 + 1; clidx2 < cluster_ct; clidx2++) {
	      if (cur_cluster_sizes[clidx2] - cur_cluster_case_cts[clidx2] > ujj) {
		set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
	      }
	    }
	    if (!ujj) {
	      for (clidx2 = cluster_ct; clidx2 < cur_cluster_ct; clidx2++) {
                if (!cur_cluster_case_cts[clidx2]) {
		  set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
		}
	      }
	    }
	  }
	}
	if (uii == 1) {
	  for (clidx1 = cluster_ct; clidx1 < cur_cluster_ct; clidx1++) {
	    if (!cur_cluster_case_cts[clidx1]) {
	      for (clidx2 = clidx1 + 1; clidx2 < cur_cluster_ct; clidx2++) {
		if (!cur_cluster_case_cts[clidx2]) {
		  set_bit_ul(tri_coord_no_diag(clidx1, clidx2), cluster_merge_prevented);
		}
	      }
	    }
	  }
	}
      }
      if (mcc_warning) {
        logerrprint("Warning: Initial cluster assignment violates --mcc restriction.\n");
      }
    }
  }

  bigstack_reset(bigstack_mark_precluster);
  if (cp->match_fname || cp->qmatch_fname) {
    retval = cluster_enforce_match(cp, missing_pheno, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, cluster_ct, cluster_starts, sample_to_cluster, cluster_merge_prevented);
    if (retval) {
      goto calc_cluster_neighbor_ret_1;
    }
  }

  tcoord = next_set_ul(cluster_merge_prevented, 0, initial_triangle_size);
  logprint("Clustering...");
  printf(" [sorting IB%c values]", cluster_missing? 'M' : 'S');
  fflush(stdout);
#ifdef __LP64__
  if (cur_cluster_ct <= 65536) {
#endif
    // Objective: Produce a list of inter-cluster IBS values sorted in
    // nonincreasing order, where we can look up the distance endpoints.
    // We could do this with qsort_ext(), but that could use twice as much
    // memory as necessary, which is bad since this is a major memory
    // bottleneck.  (E.g. with 20k samples and no pre-clustering, we're talking
    // about 20k * 10k * 12 = ~2.4GB minimum, which already overflows my dev
    // machine's 2GB default workspace allocation; so I cannot in good
    // conscience double this requirement to 4.8GB.)
    //
    // So we perform a special in-place qsort_ext():
    // 0. Ensure cluster_sorted_ibs is on top of the stack.  (This is why we
    //    use malloc more than bigstack_alloc here.)
    // 1. Convert it to an array of 12-byte structs (first 8 bytes = original
    //    IBS value, last 4 bytes = cluster indices) by allocating ~50% more
    //    memory, and copying values back-to-front.  (If there are >65536
    //    initial clusters, we will use 16-byte structs instead.)
    // 2. Remove distances corresponding to prohibited cluster merges from this
    //    array (to speed up the sort, and possibly reduce step 4's memory
    //    requirement), and fill in cluster indices.
    // 3. Sort.
    // 4. If is_group_avg, un-interleave the array-of-structs into two arrays;
    //    otherwise extract just the cluster indices.
    if (tcoord == initial_triangle_size) {
      umm = cur_cluster_ct;
    } else {
      // f(0) = 1
      // f(1) = f(2) = 2
      // f(3) = f(4) = f(5) = 3... (triangle_divide() with different rounding)
#ifdef __LP64__
      umm = (int32_t)sqrt((intptr_t)(tcoord * 2));
#else
      umm = (int32_t)sqrt(2 * ((double)((intptr_t)tcoord)));
#endif
      if (tcoord >= (umm * (umm + 1)) / 2) {
	umm++;
      }
      heap_size -= popcount_longs_nzbase(cluster_merge_prevented, tcoord / BITCT, BITCT_TO_WORDCT(initial_triangle_size));
    }
    if (!heap_size) {
      logerrprint("Error: No cluster merges possible.\n");
      goto calc_cluster_neighbor_ret_INVALID_CMDLINE;
    }
    bigstack_reset(cluster_sorted_ibs);
    if (bigstack_alloc_ui((3 + is_group_avg) * heap_size + CACHELINE_INT32, &cluster_sorted_ibs_indices)) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    ulii = initial_triangle_size - 1;
    do {
      cluster_sorted_ibs_indices[ulii * 3 + CACHELINE_INT32] = cluster_sorted_ibs_indices[ulii * 2];
      cluster_sorted_ibs_indices[ulii * 3 + CACHELINE_INT32 + 1] = cluster_sorted_ibs_indices[ulii * 2 + 1];
    } while (ulii--);

    uiptr = &(cluster_sorted_ibs_indices[CACHELINE_INT32 + 2]);
    for (uii = 1; uii < umm; uii++) {
      ukk = uii << 16;
      for (ujj = 0; ujj < uii; ujj++) {
	*uiptr = ukk | ujj;
	uiptr = &(uiptr[3]);
      }
    }
    if (uii < cur_cluster_ct) {
      umm = tcoord - ((uii * (uii - 1)) / 2);
      ukk = uii << 16;
      for (ujj = 0; ujj < umm; ujj++) {
        *uiptr = ukk | ujj;
	uiptr = &(uiptr[3]);
      }
      umm = cur_cluster_ct;
      uiptr = &(uiptr[-2]);
      do {
        uiptr2 = &(cluster_sorted_ibs_indices[CACHELINE_INT32 + ((uintptr_t)((uii * (uii - 1)) / 2)) * 3]);
	ukk = uii << 16;
	do {
	  if (!IS_SET(cluster_merge_prevented, tcoord)) {
	    *uiptr++ = uiptr2[3 * ujj];
	    *uiptr++ = uiptr2[3 * ujj + 1];
            *uiptr++ = ukk | ujj;
	  }
	  tcoord++;
	} while ((++ujj) < uii);
	ujj = 0;
      } while ((++uii) < umm);
    }
    qsort(&(cluster_sorted_ibs_indices[CACHELINE_INT32]), heap_size, 12, double_cmp_decr);
    if (!is_group_avg) {
      if (is_old_tiebreaks) {
	ulii = 1 + (heap_size / BITCT);
	if (bigstack_calloc_ul(ulii, &ibs_ties)) {
	  goto calc_cluster_neighbor_ret_NOMEM;
	}
	// copy this earlier after cluster_index allocated?
	uljj = heap_size - 1;
	for (ulii = 0; ulii < uljj; ulii++) {
	  if ((cluster_sorted_ibs_indices[ulii * 3 + CACHELINE_INT32] == cluster_sorted_ibs_indices[ulii * 3 + 3 + CACHELINE_INT32]) && (cluster_sorted_ibs_indices[ulii * 3 + 1 + CACHELINE_INT32] == cluster_sorted_ibs_indices[ulii * 3 + 4 + CACHELINE_INT32])) {
	    SET_BIT(ulii, ibs_ties);
	  }
	}
      }
      for (ulii = 0; ulii < heap_size; ulii++) {
	cluster_sorted_ibs_indices[ulii] = cluster_sorted_ibs_indices[CACHELINE_INT32 + 2 + (ulii * 3)];
      }
      if (ulii < cur_cluster_ct) {
	// this guarantees write_cluster_solution() has enough space
	ulii = cur_cluster_ct;
      }
      bigstack_reset((&(cluster_sorted_ibs_indices[round_up_pow2(ulii, CACHELINE_INT32)])));
    } else {
      uiptr = &(cluster_sorted_ibs_indices[CACHELINE_INT32 + 3 * heap_size]);
      for (ulii = 0; ulii < heap_size; ulii++) {
        *uiptr++ = cluster_sorted_ibs_indices[CACHELINE_INT32 + 2 + (ulii * 3)];
      }
      uiptr = cluster_sorted_ibs_indices;
      for (ulii = 0; ulii < heap_size; ulii++) {
        *uiptr++ = cluster_sorted_ibs_indices[CACHELINE_INT32 + ulii * 3];
        *uiptr++ = cluster_sorted_ibs_indices[CACHELINE_INT32 + ulii * 3 + 1];
      }
      bigstack_shrink_top(cluster_sorted_ibs, heap_size * sizeof(double));
      memcpy(g_bigstack_base, &(cluster_sorted_ibs_indices[CACHELINE_INT32 + 3 * heap_size]), heap_size * sizeof(int32_t));
      ulii = heap_size;
      if (ulii < cur_cluster_ct) {
        ulii = cur_cluster_ct;
      }
      bigstack_alloc_ui(ulii, &cluster_sorted_ibs_indices);
    }
    if (bigstack_calloc_ui(initial_triangle_size, &cluster_index)) {
      goto calc_cluster_neighbor_ret_NOMEM;
    }
    ujj = heap_size;
    if (!is_group_avg) {
      for (uii = 0; uii < ujj; uii++) {
	ukk = cluster_sorted_ibs_indices[uii];
	cluster_index[tri_coord_no_diag_32(ukk & 65535, ukk >> 16)] = uii;
      }
    } else {
      for (uii = 0; uii < ujj; uii++) {
	ukk = cluster_sorted_ibs_indices[uii];
	cluster_index[tri_coord_no_diag_32(ukk & 65535, ukk >> 16)] = uii + 1;
      }
    }
#ifdef __LP64__
  } else {
    logerrprint("Error: --cluster cannot handle >65536 initial clusters yet.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto calc_cluster_neighbor_ret_1;
  }
#endif

  if (bigstack_alloc_ui(2 * (sample_ct - cp->min_ct), &merge_sequence)) {
    goto calc_cluster_neighbor_ret_NOMEM;
  }
  cur_cluster_remap = (uint32_t*)malloc(cur_cluster_ct * sizeof(int32_t));
  if (!cur_cluster_remap) {
    goto calc_cluster_neighbor_ret_NOMEM;
  }
  for (clidx1 = 0; clidx1 < cur_cluster_ct; clidx1++) {
    cur_cluster_remap[clidx1] = clidx1;
  }

  if (!is_group_avg) {
    merge_ct = cluster_main(cur_cluster_ct, cluster_merge_prevented, heap_size, cluster_sorted_ibs_indices, cluster_index, cur_cluster_sizes, sample_ct, cur_cluster_case_cts, case_ct, ctrl_ct, cur_cluster_remap, cp, ibs_ties, merge_sequence);
  } else {
    merge_ct = cluster_group_avg_main(cur_cluster_ct, cluster_merge_prevented, heap_size + 1, &(cluster_sorted_ibs[-1]), &(cluster_sorted_ibs_indices[-1]), cluster_index, cur_cluster_sizes, sample_ct, cur_cluster_case_cts, case_ct, ctrl_ct, cur_cluster_remap, cp, merge_sequence);
  }
  fputs("\rClustering...", stdout);
  logprint(" done.");
  fputs("                        ", stdout);
  logprint("\n");
  retval = write_cluster_solution(outname, outname_end, sample_to_cluster, sample_ct, cluster_map, cluster_starts, late_clidx_to_sample_uidx, cluster_ct, cur_cluster_ct, sample_ids, max_sample_id_len, pheno_c, sample_idx_to_uidx, cp, cur_cluster_remap, cluster_sorted_ibs_indices, merge_ct, merge_sequence);
  if (retval) {
    goto calc_cluster_neighbor_ret_1;
  }

#ifndef NOLAPACK
  if (cp->mds_dim_ct) {
    if (!mds_plot_dmatrix_copy) {
      // --read-dists or --read-genome, and not cluster_missing
      bigstack_reset(bigstack_mark_postcluster);
      if (!is_mds_cluster) {
	ulii = sample_ct;
      } else {
	ulii = cur_cluster_ct;
      }
      if (bigstack_alloc_d((ulii * (ulii - 1)) / 2, &mds_plot_dmatrix_copy)) {
	goto calc_cluster_neighbor_ret_NOMEM;
      }
      if (is_mds_cluster || (!read_dists_fname)) {
	fill_double_zero((ulii * (ulii - 1)) / 2, mds_plot_dmatrix_copy);
      }
      if (read_dists_fname) {
	retval = read_dists(read_dists_fname, read_dists_id_fname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, is_mds_cluster? cluster_ct : 0, is_mds_cluster? cluster_starts : nullptr, is_mds_cluster? sample_to_cluster : nullptr, 2, 0, mds_plot_dmatrix_copy, 0, nullptr, nullptr);
      } else {
	retval = read_genome(read_genome_fname, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, nullptr, mds_plot_dmatrix_copy, 0, nullptr, nullptr, nullptr, 0.0, 0, is_mds_cluster? cluster_ct : 0, is_mds_cluster? cluster_starts : nullptr, is_mds_cluster? sample_to_cluster : nullptr);
      }
      if (retval) {
	goto calc_cluster_neighbor_ret_1;
      }
    } else {
      bigstack_reset(cluster_merge_prevented);
      if (is_mds_cluster) {
        cluster_dist_multiply(sample_ct, cluster_ct, cluster_starts, mds_plot_dmatrix_copy);
      }
    }
    if (cp->modifier & CLUSTER_MDS_EIGENDECOMP) {
      retval = mds_plot_eigendecomp(outname, outname_end, sample_exclude, sample_ct, sample_idx_to_uidx, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, cur_cluster_ct, merge_ct, sample_to_cluster, cur_cluster_remap, cp->mds_dim_ct, is_mds_cluster, cp->modifier & CLUSTER_MDS_EIGVALS, mds_plot_dmatrix_copy);
    } else {
      retval = mds_plot(outname, outname_end, sample_exclude, sample_ct, sample_idx_to_uidx, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, cur_cluster_ct, merge_ct, sample_to_cluster, cur_cluster_remap, cp->mds_dim_ct, is_mds_cluster, cp->modifier & CLUSTER_MDS_EIGVALS, mds_plot_dmatrix_copy);
    }
    if (retval) {
      goto calc_cluster_neighbor_ret_1;
    }
  }
#endif

  while (0) {
  calc_cluster_neighbor_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_cluster_neighbor_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_cluster_neighbor_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_cluster_neighbor_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 calc_cluster_neighbor_ret_1:
  bigstack_reset(bigstack_mark_postcluster);
  free_cond(sample_to_cluster);
  free_cond(neighbor_quantiles);
  free_cond(neighbor_qindices);
  free_cond(neighbor_quantile_means);
  free_cond(neighbor_quantile_stdev_recips);
  free_cond(ppc_fail_counts);
  free_cond(cur_cluster_case_cts);
  free_cond(cur_cluster_sizes);
  free_cond(sample_idx_to_uidx);
  free_cond(late_clidx_to_sample_uidx);
  free_cond(cur_cluster_remap);
  fclose_cond(outfile);
  return retval;
}

double regress_jack(uint32_t jackknife_d, uint32_t sample_ct, double reg_tot_xy, double reg_tot_x, double reg_tot_y, double reg_tot_xx, double reg_tot_yy, uint32_t* uibuf, double* jackknife_precomp, double* dists, double* pheno_packed, double* ret2_ptr) {
  uint32_t* uiptr = uibuf;
  uint32_t* uiptr2 = &(uibuf[jackknife_d]);
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double* dptr;
  double* dptr2;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double neg_tot_yy = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  while (uiptr < uiptr2) {
    dptr2 = &(jackknife_precomp[(*uiptr++) * JACKKNIFE_VALS_DIST]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  uiptr = uibuf;
  for (uii = 1; uii < jackknife_d; uii++) {
    ujj = *(++uiptr);
    dxx1 = pheno_packed[ujj];
    uiptr2 = uibuf;
    dptr = &(dists[(((uintptr_t)ujj) * (ujj - 1)) / 2]);
    while (uiptr2 < uiptr) {
      ukk = *uiptr2++;
      dxx = (dxx1 + pheno_packed[ukk]) * 0.5;
      dyy = dptr[ukk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = reg_tot_y - neg_tot_y;
  dyy = ((int32_t)sample_ct) - ((int32_t)jackknife_d);
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_x - neg_tot_x) / dyy) / ((reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = reg_tot_x - neg_tot_x;
  return ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_y - neg_tot_y) / dyy) / ((reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

THREAD_RET_TYPE regress_jack_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t jackknife_iters = g_jackknife_iters;
  uintptr_t sample_ct = g_sample_ct;
  uint32_t jackknife_d = g_jackknife_d;
  uint32_t* uibuf = (uint32_t*)(&(g_generic_buf[tidx * round_up_pow2(sample_ct + (jackknife_d + 1) * sizeof(int32_t), CACHELINE)]));
  unsigned char* cbuf = &(g_generic_buf[tidx * round_up_pow2(sample_ct + (jackknife_d + 1) * sizeof(int32_t), CACHELINE) + (jackknife_d + 1) * sizeof(int32_t)]);
  uintptr_t uljj = jackknife_iters / 100;
  double* jackknife_precomp = g_jackknife_precomp;
  double* dists = g_dists;
  double* pheno_packed = g_pheno_packed;
  double reg_tot_xy = g_reg_tot_xy;
  double reg_tot_x = g_reg_tot_x;
  double reg_tot_y = g_reg_tot_y;
  double reg_tot_xx = g_reg_tot_xx;
  double reg_tot_yy = g_reg_tot_yy;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uintptr_t ulii;
  double dxx;
  double ret2;
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, uibuf, sample_ct, jackknife_d, sfmtp);
    dxx = regress_jack(jackknife_d, sample_ct, reg_tot_xy, reg_tot_x, reg_tot_y, reg_tot_xx, reg_tot_yy, uibuf, jackknife_precomp, dists, pheno_packed, &ret2);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100LLU) / jackknife_iters;
      printf("\r%" PRIuPTR "%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1LLU) * jackknife_iters) / 100;
    }
  }
  g_calc_result[tidx][0] = sum;
  g_calc_result[tidx][1] = sum_sq;
  g_calc_result[tidx][2] = sum2;
  g_calc_result[tidx][3] = sum2_sq;
  THREAD_RETURN;
}

int32_t regress_distance(pthread_t* threads, uint64_t calculation_type, double* pheno_d, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t thread_ct, uintptr_t regress_iters, uint32_t regress_d) {
  unsigned char* bigstack_mark = g_bigstack_base;
  double reg_tot_xy = 0;
  double reg_tot_x = 0;
  double reg_tot_y = 0;
  double reg_tot_xx = 0;
  double reg_tot_yy = 0;
  int32_t retval = 0;
  uintptr_t ulii;
  uint32_t uii;
  double* pheno_packed;
  double* jackknife_precomp;
  double* dist_ptr;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  double* dptr5;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  double dvv;
  double duu;

  g_sample_ct = sample_ct;

  // beta = (mean(xy) - mean(x)*mean(y)) / (mean(x^2) - mean(x)^2)
  pheno_packed = (double*)alloc_and_init_collapsed_arr((char*)pheno_d, sizeof(double), unfiltered_sample_ct, sample_exclude, sample_ct, 1);
  g_pheno_packed = pheno_packed;
  if (!(calculation_type & CALC_REGRESS_REL)) {
    print_pheno_stdev(pheno_packed, sample_ct);
  }
  ulii = sample_ct;
  ulii = ulii * (ulii - 1) / 2;
  dptr4 = g_dists;
  dist_ptr = pheno_packed;
  // Linear regression slope is a function of sum(xy), sum(x), sum(y),
  // sum(x^2), and n.  To speed up the jackknife calculation, we precompute
  // (i) the global xy, x, y, x^2, and
  // (ii) the xy, x, y, x^2 for each row.
  // Then for each delete-d jackknife iteration, we take the global sums,
  // subtract the partial row sums corresponding to the deleted samples, and
  // then add back the elements in the intersection of two deletions.
  if (bigstack_calloc_d(sample_ct * JACKKNIFE_VALS_DIST, &jackknife_precomp)) {
    goto regress_distance_ret_NOMEM;
  }
  g_jackknife_precomp = jackknife_precomp;
  if (bigstack_init_sfmtp(g_thread_ct)) {
    goto regress_distance_ret_NOMEM;
  }
  for (uii = 1; uii < sample_ct; uii++) {
    dzz = *(++dist_ptr);
    dptr2 = pheno_packed;
    dptr3 = &(jackknife_precomp[uii * JACKKNIFE_VALS_DIST]);
    dptr5 = jackknife_precomp;
    do {
      dxx = (dzz + *dptr2++) * 0.5;
      dyy = (*dptr4++);
      dww = dxx * dyy;
      dvv = dxx * dxx;
      duu = dyy * dyy;
      reg_tot_xy += dww;
      *dptr3 += dww;
      *dptr5 += dww;
      dptr5++;
      reg_tot_x += dxx;
      dptr3[1] += dxx;
      *dptr5 += dxx;
      dptr5++;
      reg_tot_y += dyy;
      dptr3[2] += dyy;
      *dptr5 += dyy;
      dptr5++;
      reg_tot_xx += dvv;
      dptr3[3] += dvv;
      *dptr5 += dvv;
      dptr5++;
      reg_tot_yy += duu;
      dptr3[4] += duu;
      *dptr5 += duu;
      dptr5++;
    } while (dptr2 < dist_ptr);
  }

  g_reg_tot_xy = reg_tot_xy;
  g_reg_tot_x = reg_tot_x;
  g_reg_tot_y = reg_tot_y;
  g_reg_tot_xx = reg_tot_xx;
  g_reg_tot_yy = reg_tot_yy;

  dxx = ulii;
  LOGPRINTF("Regression slope (y = genomic distance, x = avg phenotype): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y / dxx) / (reg_tot_xx - reg_tot_x * reg_tot_x / dxx));
  LOGPRINTF("Regression slope (y = avg phenotype, x = genomic distance): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y / dxx) / (reg_tot_yy - reg_tot_y * reg_tot_y / dxx));

  g_jackknife_iters = (regress_iters + thread_ct - 1) / thread_ct;
  if (regress_d) {
    g_jackknife_d = regress_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(sample_ct);
  }
  if (bigstack_alloc_uc(thread_ct * round_up_pow2(sample_ct + (g_jackknife_d + 1) * sizeof(int32_t), CACHELINE), &g_generic_buf)) {
    goto regress_distance_ret_NOMEM;
  }
  if (spawn_threads(threads, &regress_jack_thread, thread_ct)) {
    goto regress_distance_ret_THREAD_CREATE_FAIL;
  }
  ulii = 0;
  regress_jack_thread((void*)ulii);
  dyy = g_calc_result[0][0]; // sum
  dzz = g_calc_result[0][1]; // sum of squares
  dww = g_calc_result[0][2]; // reverse regression sum
  dvv = g_calc_result[0][3]; // reverse regression sum of squares
  join_threads(threads, thread_ct);
  for (uii = 0; uii < thread_ct - 1; uii++) {
    dyy += g_calc_result[uii + 1][0];
    dzz += g_calc_result[uii + 1][1];
    dww += g_calc_result[uii + 1][2];
    dvv += g_calc_result[uii + 1][3];
  }
  regress_iters = g_jackknife_iters * thread_ct;
  putc_unlocked('\r', stdout);
  LOGPRINTF("Jackknife s.e.: %g\n", sqrt(((sample_ct - g_jackknife_d) / ((double)g_jackknife_d)) * (dzz - dyy * dyy / regress_iters) / (regress_iters - 1)));
  LOGPRINTF("Jackknife s.e. (y = avg phenotype): %g\n", sqrt(((sample_ct - g_jackknife_d) / ((double)g_jackknife_d)) * (dvv - dww * dww / regress_iters) / (regress_iters - 1)));
  while (0) {
  regress_distance_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  regress_distance_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
  }
  bigstack_reset(bigstack_mark);
  return retval;
}
