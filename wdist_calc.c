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
// between individuals j and k:
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
// typical L2 cache these days), and licenses the use of four table lookups and
// adds instead of twenty.
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
// individuals, but only when there are no missing markers, so WDIST does not
// include an implementation of that.)
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

#include "wdist_common.h"
#include "pigz.h"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#ifdef NOLAPACK
#define MATRIX_INVERT_BUF1_TYPE double
#define __CLPK_integer int
#else // not NOLAPACK
#ifndef __APPLE__
#ifdef __cplusplus
extern "C" {
#endif
  typedef double __CLPK_doublereal;
#ifdef _WIN32
#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_STRUCTURE
#include "lapack/lapacke/include/lapacke.h"
  typedef int32_t __CLPK_integer;

  void dger_(int* m, int* n, double* alpha, double* x, int* incx, double* y,
             int* incy, double* a, int* lda);

  void dgemm_(char* transa, char* transb, int* m, int* n, int* k,
              double* alpha, double* a, int* lda, double* b, int* ldb,
              double* beta, double* c, int* ldc);

  void dsymv_(char* uplo, int* n, double* alpha, double* a, int* lda,
              double* x, int* incx, double* beta, double* y, int* incy);

  double ddot_(int* n, double* dx, int* incx, double* dy, int* incy);

  void dgetrf_(__CLPK_integer* m, __CLPK_integer* n,
               __CLPK_doublereal* a, __CLPK_integer* lda,
               __CLPK_integer* ipiv, __CLPK_integer* info);

#else // not _WIN32
#include <cblas.h>
#ifdef __LP64__
  typedef int32_t __CLPK_integer;
#else
  typedef long int __CLPK_integer;
#endif
  int dgetrf_(__CLPK_integer* m, __CLPK_integer* n,
              __CLPK_doublereal* a, __CLPK_integer* lda,
              __CLPK_integer* ipiv, __CLPK_integer* info);

  int dgetri_(__CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_integer* ipiv,
              __CLPK_doublereal* work, __CLPK_integer* lwork,
              __CLPK_integer* info);

  double dlange_(char* norm, __CLPK_integer* m, __CLPK_integer* n,
                 __CLPK_doublereal* a, __CLPK_integer* lda,
                 __CLPK_doublereal* work);

  int dgecon_(char* norm, __CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_doublereal* anorm,
              __CLPK_doublereal* rcond, __CLPK_doublereal* work,
              __CLPK_integer* iwork, __CLPK_integer* info);
#endif

  void xerbla_(void) {} // fix static linking error
#ifdef __cplusplus
}
#endif // __cplusplus
#endif // __APPLE__
#define MATRIX_INVERT_BUF1_TYPE __CLPK_integer
#endif // NOLAPACK

#define MATRIX_SINGULAR_RCOND 1e-14

// number of different types of jackknife values to precompute (x^2, x, y, xy)
#define JACKKNIFE_VALS_REL 5

// Number of snp-major .bed lines to read at once for distance calc if exponent
// is zero.  Currently assumed to be a multiple of 192, and no larger than
// 1920, by the popcount_..._multiword functions.  (The optimal value depends
// on both system-specific properties such as cache sizes, as well as the
// number of individuals in the current calculation, so in principle it's best
// to select this value at runtime.  But 960 usually works well in practice in
// my experience.)
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)

#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)

// Must be multiple of 384, no larger than 3840.
#define GENOME_MULTIPLEX 1152
#define GENOME_MULTIPLEX2 (GENOME_MULTIPLEX * 2)

#define MAX_EM_ACCEL 100.0

void update_rel_ibc(double* rel_ibc, uintptr_t* geno, double* set_allele_freqs, int32_t ibc_type, uint32_t indiv_ct) {
  // first calculate weight array, then loop
  int32_t ii;
  int32_t jj;
  int32_t kk;
  uint32_t uii;
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
  fill_double_zero(wtarr, BITCT2 * 5);
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if ((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[ii * 8] = 2;
          wtarr[ii * 8 + 2] = 2.0 - 1.0 / (2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]));
          wtarr[ii * 8 + 3] = 2;
        } else {
          twt = 2 * set_allele_freqs[ii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - set_allele_freqs[ii]));
          }
          wtarr[ii * 8] = twt * twt * mult;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - set_allele_freqs[ii];
        mult = 1 / (set_allele_freqs[ii] * twt);
        wtarr[ii * 8] = 1.0 + set_allele_freqs[ii] * set_allele_freqs[ii] * mult;
        wtarr[ii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * set_allele_freqs[ii];
          wtarr[ii * 8] = twt * twt;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  wtarr[ii * 8 + 2] = INFINITY;
          if (set_allele_freqs[ii] == 0.0) {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          } else {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (set_allele_freqs[ii] == 0.0) {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        } else {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      }
    }
  }
  for (kk = 0; kk < (BITCT * 5) / 32; kk++) {
    wtptr = &(wtarr[16 * kk]);
#ifdef __LP64__
    if ((kk == 2) || (kk == 7)) {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
	wptr = &(wptr[8]);
      }
    } else {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
      }
    }
#else
    if (kk == 2) {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
	wptr = &(wptr[8]);
      }
    } else {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
      }
    }
#endif
  }
  for (uii = 0; uii < indiv_ct; uii++) {
    ulii = *geno++;
#ifdef __LP64__
    *rel_ibc += weights9[ulii >> 57] + weights8[(ulii >> 51) & 63] + weights7[(ulii >> 44) & 127] + weights6[(ulii >> 38) & 63] + weights5[(ulii >> 32) & 63] + weights4[(ulii >> 25) & 63] + weights3[(ulii >> 19) & 63] + weights2[(ulii >> 12) & 127] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 25] + weights3[(ulii >> 19) & 63] + weights2[(ulii >> 12) & 127] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

void update_rel_f_ibc(float* rel_ibc, uintptr_t* geno, float* set_allele_freqs, int32_t ibc_type, uint32_t indiv_ct) {
  // first calculate weight array, then loop
  int32_t ii;
  int32_t jj;
  int32_t kk;
  uint32_t uii;
  float twt;
  float* wtptr;
  float mult = 1.0;
  uintptr_t ulii;
  float weights[BITCT * 12];
  float* weights1 = &(weights[64]);
  float* weights2 = &(weights[128]);
  float* weights3 = &(weights[256]);
  float* weights4 = &(weights[320]);
#ifdef __LP64__
  float* weights5 = &(weights[384]);
  float* weights6 = &(weights[448]);
  float* weights7 = &(weights[512]);
  float* weights8 = &(weights[640]);
  float* weights9 = &(weights[704]);
#endif
  float wtarr[BITCT2 * 5];
  float *wptr = weights;
  fill_float_zero(wtarr, BITCT2 * 5);
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if ((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[ii * 8] = 2;
          wtarr[ii * 8 + 2] = 2.0 - 1.0 / (2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]));
          wtarr[ii * 8 + 3] = 2;
        } else {
          twt = 2 * set_allele_freqs[ii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - set_allele_freqs[ii]));
          }
          wtarr[ii * 8] = twt * twt * mult;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - set_allele_freqs[ii];
        mult = 1 / (set_allele_freqs[ii] * twt);
        wtarr[ii * 8] = 1.0 + set_allele_freqs[ii] * set_allele_freqs[ii] * mult;
        wtarr[ii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * set_allele_freqs[ii];
          wtarr[ii * 8] = twt * twt;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  wtarr[ii * 8 + 2] = INFINITY;
          if (set_allele_freqs[ii] == 0.0) {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          } else {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (set_allele_freqs[ii] == 0.0) {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        } else {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      }
    }
  }
  for (kk = 0; kk < (BITCT * 5) / 32; kk++) {
    wtptr = &(wtarr[16 * kk]);
#ifdef __LP64__
    if ((kk == 2) || (kk == 7)) {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
	wptr = &(wptr[8]);
      }
    } else {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
      }
    }
#else
    if (kk == 2) {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
	wptr = &(wptr[8]);
      }
    } else {
      for (ii = 0; ii < 8; ii++) {
	twt = wtptr[ii + 8];
	for (jj = 0; jj < 8; jj++) {
	  *wptr++ = twt + wtptr[jj];
	}
      }
    }
#endif
  }
  for (uii = 0; uii < indiv_ct; uii++) {
    ulii = *geno++;
#ifdef __LP64__
    *rel_ibc += weights9[ulii >> 57] + weights8[(ulii >> 51) & 63] + weights7[(ulii >> 44) & 127] + weights6[(ulii >> 38) & 63] + weights5[(ulii >> 32) & 63] + weights4[(ulii >> 25) & 63] + weights3[(ulii >> 19) & 63] + weights2[(ulii >> 12) & 127] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 25] + weights3[(ulii >> 19) & 63] + weights2[(ulii >> 12) & 127] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

void fill_weights(double* weights, double* set_allele_freqs, double exponent) {
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  int32_t oo;
  double wtarr[MULTIPLEX_DIST_EXP / 2];
  double* wt;
#ifdef __LP64__
  double twt[5];
  double twtf;
  __m128d* wpairs = (__m128d*)weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
#else
  int32_t pp;
  int32_t qq;
  double twt[7];
#endif
  for (ii = 0; ii < MULTIPLEX_DIST_EXP / 2; ii++) {
    wtarr[ii] = pow(2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]), -exponent);
  }
  for (oo = 0; oo < 2; oo++) {
    wt = &(wtarr[7 * oo]);
#ifdef __LP64__
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(wt[0] * 2, wt[0]);
#endif
    twt[0] = 0;
    for (ii = 0; ii < 4; ii += 1) {
      // bizarrely, turning the ii == 2 case into a memcpy doesn't actually
      // seem to speed this up
      if (ii & 1) {
	twt[0] += wt[6];
      }
      twt[1] = twt[0];
      for (jj = 0; jj < 4; jj += 1) {
	if (jj & 1) {
	  twt[1] += wt[5];
	}
	twt[2] = twt[1];
	for (kk = 0; kk < 4; kk += 1) {
	  if (kk & 1) {
	    twt[2] += wt[4];
	  }
	  twt[3] = twt[2];
	  for (mm = 0; mm < 4; mm++) {
	    if (mm & 1) {
	      twt[3] += wt[3];
	    }
	    twt[4] = twt[3];
	    for (nn = 0; nn < 4; nn++) {
	      if (nn & 1) {
		twt[4] += wt[2];
	      }
#ifdef __LP64__
	      twtf = twt[4];
	      vpen = _mm_set1_pd(twtf);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
	      twtf += wt[1];
	      vpen = _mm_set1_pd(twtf);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
	      *wpairs = *(wpairs - 2);
	      wpairs++;
	      *wpairs = *(wpairs - 2);
	      wpairs++;
	      vpen = _mm_set1_pd(twtf + wt[1]);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
#else
              twt[5] = twt[4];
              for (pp = 0; pp < 4; pp++) {
                if (pp & 1) {
                  twt[5] += wt[1];
                }
                twt[6] = twt[5];
                for (qq = 0; qq < 4; qq++) {
                  if (qq & 1) {
                    twt[6] += wt[0];
                  }
                  *weights++ = twt[6];
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
  for (oo = 0; oo < 3; oo++) {
    wt = &(wtarr[14 + 6 * oo]);
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(2 * wt[0], wt[0]);
    twt[0] = 0;
    for (ii = 0; ii < 4; ii += 1) {
      if (ii & 1) {
        twt[0] += wt[5];
      }
      twt[1] = twt[0];
      for (jj = 0; jj < 4; jj += 1) {
        if (jj & 1) {
          twt[1] += wt[4];
        }
        twt[2] = twt[1];
	for (kk = 0; kk < 4; kk += 1) {
          if (kk & 1) {
            twt[2] += wt[3];
          }
          twt[3] = twt[2];
          for (mm = 0; mm < 4; mm++) {
	    if (mm & 1) {
	      twt[3] += wt[2];
	    }
	    twtf = twt[3];
	    vpen = _mm_set1_pd(twtf);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
	    twtf += wt[1];
	    vpen = _mm_set1_pd(twtf);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
	    *wpairs = *(wpairs - 2);
	    wpairs++;
	    *wpairs = *(wpairs - 2);
	    wpairs++;
	    vpen = _mm_set1_pd(twtf + wt[1]);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
          }
	}
      }
    }
  }
#endif
}

void fill_weights_r(double* weights, double* set_allele_freqs, int32_t var_std) {
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  // 20 markers to process in quintuplets, for 64-bit; 10, for 32-bit.
  // Each quintuplet of markers requires 40 wtarr entries, and induces
  // 2^15 writes to weights[].
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
  __m128d* wpairs = (__m128d*)weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
  __m128d vfinal3;
  __m128d vfinal4;
#else
  int32_t oo;
#endif
  if (((uintptr_t)wtarr) & 15) {
    // force 16-byte alignment; can't do this at compile-time since stack
    // pointer has no 16-byte align guarantee.
    // yes, this assumes doubles are 8 bytes.
    wtarr++;
  }
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if (((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) || (!var_std)) {
      if (set_allele_freqs[ii] < 0.5) {
	mean = 2 * set_allele_freqs[ii];
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * (1.0 - set_allele_freqs[ii]));
        }
        aux = mean * mult;
	wtarr[ii * 8] = mean * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean_m2 * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m2 * mean_m1 * mult;
	wtarr[ii * 8 + 6] = mean_m2 * mean_m2 * mult;
      } else {
	mean = 2 * (1.0 - set_allele_freqs[ii]);
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * set_allele_freqs[ii]);
        }
        aux = mean_m2 * mult;
	wtarr[ii * 8] = mean_m2 * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m1 * mean * mult;
	wtarr[ii * 8 + 6] = mean * mean * mult;
      }
    } else {
      if (set_allele_freqs[ii] == 0.0) {
        wtarr[ii * 8] = 0;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = -1;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = INFINITY;
        wtarr[ii * 8 + 6] = INFINITY;
      } else {
        wtarr[ii * 8] = INFINITY;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = INFINITY;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = -1;
        wtarr[ii * 8 + 6] = 0;
      }
    }
    wtarr[ii * 8 + 7] = 0;
  }
  for (nn = 0; nn < BITCT / 16; nn++) {
    wtptr = &(wtarr[40 * nn]);
#ifdef __LP64__
    vfinal1 = _mm_load_pd(wtptr);
    vfinal2 = _mm_load_pd(&(wtptr[2]));
    vfinal3 = _mm_load_pd(&(wtptr[4]));
    vfinal4 = _mm_load_pd(&(wtptr[6]));
#endif
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 32];
      for (jj = 0; jj < 8; jj++) {
        twt2 = twt + wtptr[jj + 24];
        for (kk = 0; kk < 8; kk++) {
          twt3 = twt2 + wtptr[kk + 16];
          for (mm = 0; mm < 8; mm++) {
            twt4 = twt3 + wtptr[mm + 8];
#ifdef __LP64__
            vpen = _mm_set1_pd(twt4);
            *wpairs++ = _mm_add_pd(vpen, vfinal1);
            *wpairs++ = _mm_add_pd(vpen, vfinal2);
            *wpairs++ = _mm_add_pd(vpen, vfinal3);
            *wpairs++ = _mm_add_pd(vpen, vfinal4);
#else
            for (oo = 0; oo < 8; oo++) {
              *weights++ = twt4 + wtptr[oo];
            }
#endif
          }
        }
      }
    }
  }
}

void fill_weights_r_f(float* weights_f, float* set_allele_freqs_f, int32_t var_std) {
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  // 20 markers to process in quintuplets, for 64-bit; 10, for 32-bit.
  // Each quintuplet of markers requires 40 wtarr entries, and induces
  // 2^15 writes to weights_f[].
  float wtarr_raw[BITCT2 * 5 + 3];
  float* wtarr = wtarr_raw;
  float twt;
  float twt2;
  float twt3;
  float twt4;
  float* wtptr;
  float mean;
  float mean_m1;
  float mean_m2;
  float mult = 1.0;
  float aux;
#ifdef __LP64__
  __m128* wquads = (__m128*)weights_f;
  __m128 vpen;
  __m128 vfinal1;
  __m128 vfinal2;
#else
  int32_t oo;
#endif
  ii = (((uintptr_t)wtarr) & 15);
  if (ii) {
    // force 16-byte alignment; can't do this at compile-time since stack
    // pointer has no 16-byte align guarantee.
    // yes, this assumes floats are 4 bytes.
    wtarr = &(wtarr[4 - (ii / 4)]);
  }
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if (((set_allele_freqs_f[ii] != 0.0) && (set_allele_freqs_f[ii] < (1.0 - EPSILON))) || (!var_std)) {
      if (set_allele_freqs_f[ii] < 0.5) {
	mean = 2 * set_allele_freqs_f[ii];
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * (1.0 - set_allele_freqs_f[ii]));
        }
        aux = mean * mult;
	wtarr[ii * 8] = mean * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean_m2 * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m2 * mean_m1 * mult;
	wtarr[ii * 8 + 6] = mean_m2 * mean_m2 * mult;
      } else {
	mean = 2 * (1.0 - set_allele_freqs_f[ii]);
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * set_allele_freqs_f[ii]);
        }
        aux = mean_m2 * mult;
	wtarr[ii * 8] = mean_m2 * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m1 * mean * mult;
	wtarr[ii * 8 + 6] = mean * mean * mult;
      }
    } else {
      if (set_allele_freqs_f[ii] == 0.0) {
        wtarr[ii * 8] = 0;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = -1;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = INFINITY;
        wtarr[ii * 8 + 6] = INFINITY;
      } else {
        wtarr[ii * 8] = INFINITY;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = INFINITY;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = -1;
        wtarr[ii * 8 + 6] = 0;
      }
    }
    wtarr[ii * 8 + 7] = 0;
  }
  for (nn = 0; nn < BITCT / 16; nn++) {
    wtptr = &(wtarr[40 * nn]);
#ifdef __LP64__
    vfinal1 = _mm_load_ps(wtptr);
    vfinal2 = _mm_load_ps(&(wtptr[4]));
#endif
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 32];
      for (jj = 0; jj < 8; jj++) {
        twt2 = twt + wtptr[jj + 24];
        for (kk = 0; kk < 8; kk++) {
          twt3 = twt2 + wtptr[kk + 16];
          for (mm = 0; mm < 8; mm++) {
            twt4 = twt3 + wtptr[mm + 8];
#ifdef __LP64__
            vpen = _mm_set1_ps(twt4);
            *wquads++ = _mm_add_ps(vpen, vfinal1);
            *wquads++ = _mm_add_ps(vpen, vfinal2);
#else
            for (oo = 0; oo < 8; oo++) {
              *weights_f++ = twt4 + wtptr[oo];
            }
#endif
          }
        }
      }
    }
  }
}

/*
void collapse_chrom_marker_idxs(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, int32_t unfiltered_marker_ct) {
  uint32_t* chrom_fo = chrom_info_ptr->chrom_file_order;
  uint32_t* chrom_fo_midxs = chrom_info_ptr->chrom_file_order_marker_idx;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t midx = 0;
  int32_t new_midx = 0;
  int32_t chrom_end_midx;
  uint32_t chrom_fo_idx;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_fo_midxs[chrom_fo_idx] = new_midx;
    chrom_info_ptr->chrom_start[chrom_fo[chrom_fo_idx]] = new_midx;
    chrom_end_midx = chrom_fo_midxs[chrom_fo_idx + 1];
    for (; midx < chrom_end_midx; midx++) {
      if (!is_set(marker_exclude, midx)) {
	new_midx++;
      }
    }
    // todo: collapse when chromosome completely eliminated
    chrom_info_ptr->chrom_end[chrom_fo[chrom_fo_idx]] = new_midx;
  }
  chrom_fo_midxs[chrom_fo_idx] = new_midx;
}
*/

  /*
// (using this as a dumping ground for old code that is likely enough to be
// useful later that I don't want to be forced to dredge it from git)

    collapse_arr(marker_alleles, 2, marker_exclude, unfiltered_marker_ct);
    sscanf(output_missing_pheno, "%lg", &missing_phenod);
    // if this becomes much more of a maintenance nightmare, consider exiting
    // function and reloading from .bed from scratch
    if (g_pheno_c) {
      collapse_arr(g_pheno_c, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    } else if (g_pheno_d) {
      collapse_arr((char*)g_pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
    }
    if (sex_info) {
      collapse_arr((char*)sex_info, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    }
    collapse_bitarr(marker_reverse, marker_exclude, unfiltered_marker_ct);
    if ((calculation_type & CALC_WRITE_SNPLIST) || ((calculation_type & CALC_RECODE) && (recode_modifier & RECODE_LGEN))) {
      collapse_arr(marker_ids, max_marker_id_len, marker_exclude, unfiltered_marker_ct);
    }
    if (calculation_type & (CALC_WRITE_SNPLIST | CALC_GENOME | CALC_LD_PRUNE)) {
      collapse_chrom_marker_idxs(chrom_info_ptr, marker_exclude, unfiltered_marker_ct);
      if (calculation_type & (CALC_GENOME | CALC_LD_PRUNE)) {
	collapse_arr((char*)marker_pos, sizeof(int32_t), marker_exclude, unfiltered_marker_ct);
      }
    }
    collapse_arr(person_ids, max_person_id_len, indiv_exclude, unfiltered_indiv_ct);
    if (fam_col_34) {
      collapse_arr(paternal_ids, max_paternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_arr(maternal_ids, max_maternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_bitarr(founder_info, indiv_exclude, unfiltered_indiv_ct);
    }
    if (wt_needed) {
      collapse_arr((char*)g_marker_weights, sizeof(double), marker_exclude, unfiltered_marker_ct);
    }
    collapse_arr((char*)set_allele_freqs, sizeof(double), marker_exclude, unfiltered_marker_ct);
    unfiltered_marker_ct -= marker_exclude_ct;
    marker_exclude_ct = 0;
    fill_ulong_zero(marker_exclude, (unfiltered_marker_ct + (BITCT - 1)) / BITCT);
    unfiltered_indiv_ct -= indiv_exclude_ct;
    unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
    indiv_exclude_ct = 0;
    fill_ulong_zero(indiv_exclude, (unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  }
  if (!allow_no_sex) {
    ii = indiv_exclude_ct;
    exclude_ambiguous_sex(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_info);
    ii = indiv_exclude_ct - ii;
    if (ii) {
      sprintf(logbuf, "%d individual%s with unknown sex removed (prevent with --allow-no-sex).\n", ii, (ii == 1)? "" : "s");
      logprintb();
    }
    g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
    if (!g_indiv_ct) {
      logprint("Error: No people remaining.\n");
      goto wdist_ret_INVALID_FORMAT;
    }

inline int32_t flexclose_null(FILE** outfile_ptr, gzFile* gz_outfile_ptr) {
  int32_t ii;
  if (*outfile_ptr) {
    return fclose_null(outfile_ptr);
  } else {
    ii = gzclose(*gz_outfile_ptr);
    *gz_outfile_ptr = NULL;
    return (ii != Z_OK);
  }
}

  */

inline double SQR(const double a) {
  return a * a;
}

#ifdef __cplusplus
inline double SIGN(const double &a, const double &b) {
  // PLINK helper.h SIGN() template specialized to doubles.
  return (b >= 0)? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#else
inline double SIGN(const double a, const double b) {
  // PLINK helper.h SIGN() template specialized to doubles.
  return (b >= 0)? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#endif

double pythag(const double a, const double b) {
  // PLINK stats.cpp pythag().
  double absa,absb;
 
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

int32_t svdcmp_c(int32_t m, double* a, double* w, double* v) {
  // C port of PLINK stats.cpp svdcmp().
  // Note that this function is NOT thread-safe, due to the buffer allocated
  // from the workspace stack.  Pass in a preallocated buffer if that's not
  // okay.
  unsigned char* wkspace_mark = wkspace_base;
  int32_t n = m;
  int32_t flag;
  int32_t l = 0; // suppress compile warning
  int32_t i,its,j,jj,k,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double volatile temp;
  double* rv1;
  if (wkspace_alloc_d_checked(&rv1, m * sizeof(double))) {
    return -1;
  }

  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k * m + i]);
      if (scale != 0.0) {
	for (k=i;k<m;k++) {
	  a[k * m + i] /= scale;
	  s += a[k * m + i]*a[k * m + i];
	}
	f=a[i * m + i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i * m + i]=f-g;
	for (j=l-1;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a[k * m + i]*a[k * m + j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k * m + j] += f*a[k * m + i];
	}
	for (k=i;k<m;k++) a[k * m + i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k<n;k++) scale += fabs(a[i * m + k]);
      if (scale != 0.0) {
	for (k=l-1;k<n;k++) {
	  a[i * m + k] /= scale;
	  s += a[i * m + k]*a[i * m + k];
	}
	f=a[i * m + l-1];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i * m + l-1]=f-g;
	for (k=l-1;k<n;k++) rv1[k]=a[i * m + k]/h;
	for (j=l-1;j<m;j++) {
	  for (s=0.0,k=l-1;k<n;k++) s += a[j * m + k]*a[i * m + k];
	  for (k=l-1;k<n;k++) a[j * m + k] += s*rv1[k];
	}
	for (k=l-1;k<n;k++) a[i * m + k] *= scale;
      }
    }
    anorm=MAXV(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
	for (j=l;j<n;j++)
	  v[j * m + i]=(a[i * m + j]/a[i * m + l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i * m + k]*v[k * m + j];
	  for (k=l;k<n;k++) v[k * m + j] += s*v[k * m + i];
	}
      }
      for (j=l;j<n;j++) v[i * m + j]=v[j * m + i]=0.0;
    }
    v[i * m + i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MINV(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i * m + j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k * m + i]*a[k * m + j];
	f=(s/a[i * m + i])*g;
	for (k=i;k<m;k++) a[k * m + j] += f*a[k * m + i];
      }
      for (j=i;j<m;j++) a[j * m + i] *= g;
    } else for (j=i;j<m;j++) a[j * m + i]=0.0;
    ++a[i * m + i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	temp=fabs(rv1[l])+anorm;
	if (temp == anorm) {
	  flag=0;
	  break;
	}
	temp=fabs(w[nm])+anorm;
	if (temp == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<k+1;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  temp = fabs(f)+anorm;
	  if (temp == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j * m + nm];
	    z=a[j * m + i];
	    a[j * m + nm]=y*c+z*s;
	    a[j * m + i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j * m + k] = -v[j * m + k];
	}
	break;
      }
      if (its == 29) 
	return 0; // cannot converge: multi-collinearity?
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj * m + j];
	  z=v[jj * m + i];
	  v[jj * m + j]=x*c+z*s;
	  v[jj * m + i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj * m + j];
	  z=a[jj * m + i];
	  a[jj * m + j]=y*c+z*s;
	  a[jj * m + i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  wkspace_reset(wkspace_mark);
  return 1;
}

int32_t invert_matrix_trunc_singular_svd(__CLPK_integer dim, double* matrix, double* dbl_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim) {
  // for no truncation, call this with min_dim == dim
  const double eps = 1e-24;
  unsigned char* wkspace_mark = wkspace_base;
  double* matrix_copy;
  double* small_matrix = NULL;
  uint32_t min_dim_u = min_dim;
  int32_t flag;
  uint32_t cur_dim;
  uint32_t max_dim; // known singular
  uint32_t uii;
  int32_t i;
  int32_t j;
  int32_t k;

  if (wkspace_alloc_d_checked(&matrix_copy, dim * dim * sizeof(double))) {
    return -1;
  }
  memcpy(matrix_copy, matrix, dim * dim * sizeof(double));
  flag = svdcmp_c(dim, matrix, dbl_1d_buf, dbl_2d_buf);

  if (!flag) {
    max_dim = dim;
    if (dim > min_dim + 1) {
      if (wkspace_alloc_d_checked(&small_matrix, (dim - 1) * (dim - 1) * sizeof(double))) {
        return -1;
      }
    }
    while (max_dim > min_dim_u + 1) {
      cur_dim = (max_dim + min_dim_u) / 2;
      for (uii = 0; uii < cur_dim; uii++) {
        memcpy(&(small_matrix[uii * cur_dim]), &(matrix_copy[uii * dim]), cur_dim * sizeof(double));
      }
      flag = svdcmp_c(cur_dim, small_matrix, dbl_1d_buf, dbl_2d_buf);
      if (flag) {
        max_dim = cur_dim;
      } else {
        min_dim_u = cur_dim;
      }
    }
    wkspace_reset(wkspace_mark);
    return min_dim_u;
  }
  // Look for singular values
  double wmax = 0;
  for (i=0; i<dim; i++) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i<dim; i++) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : 1/dbl_1d_buf[i];
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  // [nxn].[t(v)] 
  for (i=0; i<dim; i++) {
    fill_double_zero(dbl_1d_buf, dim);
    for (j=0; j<dim; j++) {
      for (k=0; k<dim; k++) {
	dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j < dim; j++) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

#ifdef NOLAPACK
int32_t invert_matrix(int32_t dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* dbl_1d_buf, double* dbl_2d_buf) {
  // C port of PLINK stats.cpp's svd_inverse() function.

  // w -> dbl_1d_buf
  // v -> dbl_2d_buf
  const double eps = 1e-24;
  int32_t i;
  int32_t j;
  int32_t k;
  if (svdcmp_c(dim, matrix, dbl_1d_buf, dbl_2d_buf) == -1) {
    return -1;
  }

  // Look for singular values
  double wmax = 0;
  for (i=0; i<dim; i++) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i<dim; i++) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : 1/dbl_1d_buf[i];
  }
  
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  // [nxn].[t(v)] 
  for (i=0; i<dim; i++) {
    fill_double_zero(dbl_1d_buf, dim);
    for (j=0; j<dim; j++) {
      for (k=0; k<dim; k++) {
	dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j < dim; j++) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  return 0;
}

int32_t invert_matrix_trunc_singular(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim) {
  return invert_matrix_trunc_singular_svd(dim, matrix, int_1d_buf, dbl_2d_buf, min_dim);
}
#else
int32_t invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf) {
  // dgetrf_/dgetri_ is more efficient than dpotrf_/dpotri_ on OS X.
  __CLPK_integer lwork = dim * dim;
  __CLPK_integer info;
  dgetrf_(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  dgetri_(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  // todo: check if there are any major error conditions, return -1 on them
  return 0;
}

int32_t invert_matrix_trunc_singular(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim) {
  // Variant of invert_matrix which checks if the matrix is singular.  If it
  // is, this determines the minimum number of rows/columns which can be
  // removed from the bottom right to make the matrix nonsingular, knowing that
  // the upper left min_dim * min_dim matrix is nonsingular, and returns the
  // zero-based index of the first removed row/column.
  __CLPK_integer lwork = dim * dim;
  char cc = '1';
  unsigned char* wkspace_mark = wkspace_base;
  double* matrix_copy;
  __CLPK_integer info;
  double norm;
  double rcond;
  norm = dlange_(&cc, &dim, &dim, matrix, &dim, dbl_2d_buf);
  if (wkspace_alloc_d_checked(&matrix_copy, dim * dim * sizeof(double))) {
    return -1;
  }
  memcpy(matrix_copy, matrix, dim * dim * sizeof(double));
  dgetrf_(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  if (info > 0) {
    rcond = 0.0;
  } else {
    dgecon_(&cc, &dim, matrix, &dim, &norm, &rcond, dbl_2d_buf, &(int_1d_buf[dim]), &info);
  }
  if (rcond < MATRIX_SINGULAR_RCOND) {
    memcpy(matrix, matrix_copy, dim * dim * sizeof(double));
    wkspace_reset(wkspace_mark);
    return invert_matrix_trunc_singular_svd(dim, matrix, (double*)int_1d_buf, dbl_2d_buf, min_dim);
  } else {
    dgetri_(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  }
  wkspace_reset(wkspace_mark);
  return 0;
}
#endif

int32_t get_chrom_end(Chrom_info* chrom_info_ptr, uintptr_t marker_idx) {
  return chrom_info_ptr->chrom_end[get_marker_chrom(chrom_info_ptr, marker_idx)];
}

void exclude_multi(uintptr_t* exclude_arr, int32_t* new_excl, uintptr_t indiv_ct, uintptr_t* exclude_ct_ptr) {
  uint32_t uii;
  int32_t true_loc = 0;
  for (uii = 0; uii < indiv_ct; uii++) {
    true_loc = next_non_set_unsafe(exclude_arr, true_loc);
    if (new_excl[uii] == -1) {
      set_bit(exclude_arr, true_loc, exclude_ct_ptr);
    }
    true_loc++;
  }
}

void collapse_copy_phenod(double *target, double* pheno_d, uintptr_t* indiv_exclude, uintptr_t indiv_ct) {
  int32_t ii = 0;
  double* target_end = &(target[indiv_ct]);
  while (target < target_end) {
    ii = next_non_set_unsafe(indiv_exclude, ii);
    *target++ = pheno_d[ii++];
  }
}

#ifdef __LP64__
// XOR + mask variants of vectorized Lauradoux/Walisch popcount.  (See
// popcount_vecs() in wdist_common.c for basic documentation.)
// Note that the size of the popcounted buffer is a hardcoded constant
// (specifically, (MULTIPLEX_DIST / BITCT) * 16 bytes).  The current code
// assumes (MULTIPLEX_DIST / BITCT) is a multiple of 3, and no greater than 30.
static inline uint32_t popcount_xor_1mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** maskp) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
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
  __uni16 acc;
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

static inline void ld_dot_prod(__m128i* vec1, __m128i* vec2, __m128i* mask1, __m128i* mask2, int32_t* return_vals, int32_t iters) {
  // Main routine for computation of \sum_i^M (x_i - \mu_x)(y_i - \mu_y), where
  // x_i, y_i \in \{-1, 0, 1\}, but there are missing values.
  //
  //
  // We decompose this sum into
  //   \sum_i x_iy_i - \mu_y\sum_i x_i - \mu_x\sum_i y_i +
  //   (M - # missing)\mu_x\mu_y.
  // *Without* missing values, this can be handled very cleanly.  The last
  // three terms can all be precomputed, and \sum_i x_iy_i can be handled in a
  // manner very similar to bitwise Hamming distance.  This is several times as
  // fast as the lookup tables used for relationship matrices.
  //
  // Unfortunately, when missing values are present, \mu_y\sum_i x_i and
  // \mu_x\sum_i y_i must be handled in the main loop, and this removes much of
  // the speed advantage.  So the best applications of the underlying ternary
  // dot product algorithm used here lie elsewhere.  Nevertheless, it is still
  // faster, so we use it.
  //
  //
  // Input:
  // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
  // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
  //   nonmissing).
  // * return_vals provides space for return values.
  // * iters is the number of 48-byte windows to process, anywhere from 1 to 10
  //   inclusive.  (No, this is not the interface you'd use for a
  //   general-purpose library.)
  //
  // This function performs the update
  //   return_vals[0] += (-N) + \sum_i x_iy_i
  //   return_vals[1] += N_y + \sum_i x_i
  //   return_vals[2] += N_x + \sum_i y_i
  // where N is the number of individuals processed after applying the
  // missingness masks indicated by the subscripts.  The calculation currently
  // proceeds as follows:
  //
  // 1. N + \sum_i x_i = popcount_variant(vec1 & mask2)
  // The "variant" suffix refers to starting with two-bit integers instead of
  // one-bit integers in our summing process, so we get to skip a few
  // operations.  (Once all researchers are using machines with fast hardware
  // popcount, a slightly different implementation would be better.)
  //
  // 2. zcheck := (vec1 | vec2) & 0x5555...
  // Detects whether at least one member of the pair has a 0/missing value.
  //
  // 3. popcount_variant(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
  // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i dot
  // product.
  //
  // MULTIPLEX_LD sets of values are handled per function call.  If fewer
  // values are present, it is currently safe to zero out the ends of all
  // input vectors.

  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader1;
  __m128i loader2;
  __m128i sum1;
  __m128i sum2;
  __m128i sum12;
  __m128i tmp_sum1;
  __m128i tmp_sum2;
  __m128i tmp_sum12;
  __uni16 acc;
  __uni16 acc1;
  __uni16 acc2;
  acc.vi = _mm_setzero_si128();
  acc1.vi = _mm_setzero_si128();
  acc2.vi = _mm_setzero_si128();
  do {
    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    sum1 = _mm_and_si128(sum1, loader1);
    sum2 = _mm_and_si128(sum2, loader2);
    // use andnot to eliminate need for 0xaaaa... to occupy an xmm register
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
    sum12 = _mm_or_si128(sum12, loader1);

    // sum1, sum2, and sum12 now store the (biased) two-bit sums of
    // interest
    sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2), _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
    sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2), _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
    sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
    tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
    tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
    loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
    tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

    sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
    sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
    sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

    acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
    acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
  } while (--iters);
  // moved down since we're out of xmm registers
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
#if MULTIPLEX_LD > 960
  acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
  acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
  acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
  return_vals[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
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

static inline void ld_dot_prod(uintptr_t* vec1, uintptr_t* vec2, uintptr_t* mask1, uintptr_t* mask2, int32_t* return_vals, int32_t iters) {
  uint32_t final_sum1 = 0;
  uint32_t final_sum2 = 0;
  uint32_t final_sum12 = 0;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t sum1;
  uintptr_t sum2;
  uintptr_t sum12;
  uintptr_t tmp_sum1;
  uintptr_t tmp_sum2;
  uintptr_t tmp_sum12;
  do {
    // (The important part of the header comment on the 64-bit version is
    // copied below.)
    //
    // Input:
    // * vec1 and vec2 are encoded -1 -> 00, 0/missing -> 01, 1 -> 10.
    // * mask1 and mask2 mask out missing values (i.e. 00 for missing, 11 for
    //   nonmissing).
    // * return_vals provides space for return values.
    // * iters is the number of 12-byte windows to process, anywhere from 1 to
    //   40 inclusive.  (No, this is not the interface you'd use for a
    //   general-purpose library.)  [32- and 64-bit differ here.]
    //
    // This function performs the update
    //   return_vals[0] += (-N) + \sum_i x_iy_i
    //   return_vals[1] += N_y + \sum_i x_i
    //   return_vals[2] += N_x + \sum_i y_i
    // where N is the number of individuals processed after applying the
    // missingness masks indicated by the subscripts.  The calculation
    // currently proceeds as follows:
    //
    // 1. N + \sum_i x_i = popcount_variant(vec1 & mask2)
    // The "variant" suffix refers to starting with two-bit integers instead of
    // one-bit integers in our summing process, so we get to skip a few
    // operations.  (Once all reserachers are using machines with fast hardware
    // popcount, a slightly different implementation would be better.)
    //
    // 2. zcheck := (vec1 | vec2) & 0x5555...
    // Detects whether at least one member of the pair has a 0/missing value.
    //
    // 3. popcount_variant(((vec1 ^ vec2) & (0xaaaa... - zcheck)) | zcheck)
    // Subtracting this *from* a bias will give us our desired \sum_i x_iy_i
    // dot product.

    loader1 = *vec1++;
    loader2 = *vec2++;
    sum1 = *mask2++;
    sum2 = *mask1++;
    sum12 = (loader1 | loader2) & FIVEMASK;

    sum1 = sum1 & loader1;
    sum2 = sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
    sum12 = sum12 | loader1;

    sum1 = (sum1 & 0x33333333) + ((sum1 >> 2) & 0x33333333);
    sum2 = (sum2 & 0x33333333) + ((sum2 >> 2) & 0x33333333);
    sum12 = (sum12 & 0x33333333) + ((sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;

    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    loader1 = *vec1++;
    loader2 = *vec2++;
    tmp_sum1 = *mask2++;
    tmp_sum2 = *mask1++;
    tmp_sum12 = (loader1 | loader2) & FIVEMASK;

    tmp_sum1 = tmp_sum1 & loader1;
    tmp_sum2 = tmp_sum2 & loader2;
    loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
    tmp_sum12 = tmp_sum12 | loader1;

    sum1 += (tmp_sum1 & 0x33333333) + ((tmp_sum1 >> 2) & 0x33333333);
    sum2 += (tmp_sum2 & 0x33333333) + ((tmp_sum2 >> 2) & 0x33333333);
    sum12 += (tmp_sum12 & 0x33333333) + ((tmp_sum12 >> 2) & 0x33333333);

    sum1 = (sum1 & 0x0f0f0f0f) + ((sum1 >> 4) & 0x0f0f0f0f);
    sum2 = (sum2 & 0x0f0f0f0f) + ((sum2 >> 4) & 0x0f0f0f0f);
    sum12 = (sum12 & 0x0f0f0f0f) + ((sum12 >> 4) & 0x0f0f0f0f);

    // technically could do the multiply-and-shift only once every two rounds
    final_sum1 += (sum1 * 0x01010101) >> 24;
    final_sum2 += (sum2 * 0x01010101) >> 24;
    final_sum12 += (sum12 * 0x01010101) >> 24;
  } while (--iters);
  return_vals[0] -= final_sum12;
  return_vals[1] += final_sum1;
  return_vals[2] += final_sum2;
}
#endif

// ----- multithread globals -----
double* g_rel_dists = NULL;
uint32_t* g_indiv_missing_unwt = NULL;
uint32_t* g_missing_dbl_excluded = NULL;
double* g_dists = NULL;

static uint32_t g_thread_start[MAX_THREADS_P1];
static float* g_rel_f_dists = NULL;
static int32_t* g_idists;
static uintptr_t* g_pheno_nm = NULL;
static uintptr_t* g_pheno_c = NULL;
static unsigned char* g_geno = NULL;
// see incr_dists()
#if (2048 * BITCT) < 45056
#error "Insufficient initial size for g_weights[]."
#endif
static double g_weights[2048 * BITCT];
static float* g_weights_f = (float*)g_weights;
static uint32_t* g_weights_i = (uint32_t*)g_weights;
static double g_reg_tot_xy;
static double g_reg_tot_x;
static double g_reg_tot_y;
static double g_reg_tot_xx;
static double g_reg_tot_yy;
static uint32_t g_ctrl_ct;
static uint32_t g_case_ct;
static uintptr_t g_jackknife_iters;
static uint32_t g_jackknife_d;
static double g_calc_result[MAX_THREADS][9];
static uintptr_t* g_masks;
static uintptr_t* g_mmasks;
static uint32_t* g_missing_tot_weights;
static uint32_t* g_indiv_missing;
static double* g_jackknife_precomp = NULL;
static uint32_t* g_genome_main;
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

void ibs_test_init_col_buf(uintptr_t row_idx, uintptr_t* perm_col_buf) {
  uintptr_t perm_idx = 0;
  uintptr_t block_size = BITCT;
  uintptr_t rshift_ct = row_idx & (BITCT - 1);
  uintptr_t* gpr_ptr;
  uintptr_t ulii;
  uintptr_t block_pos;
  gpr_ptr = &(g_perm_rows[(row_idx / BITCT) * g_perm_ct]);
  do {
    if (perm_idx + BITCT > g_perm_ct) {
      block_size = g_perm_ct - perm_idx;
    }
    ulii = 0;
    block_pos = 0;
    do {
      ulii |= (((*gpr_ptr++) >> rshift_ct) & ONELU) << block_pos;
    } while (++block_pos < block_size);
    perm_idx += block_size;
    *perm_col_buf++ = ulii;
  } while (perm_idx < g_perm_ct);
}

double fill_psbuf(uintptr_t block_size, double* dists, uintptr_t* col_uidxp, double* psbuf, double* ssq0p) {
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
      col_uidx = next_set_unsafe(g_pheno_nm, col_uidx);
      if (g_load_dists) {
	dxx = dists[col_uidx];
      } else {
	dxx = 1.0 - dists[col_uidx] * g_half_marker_ct_recip;
      }
      increment[sub_block_idx] = subtot - dxx;
      subtot += dxx;
      ssq[is_set(g_pheno_c, col_uidx)] += dxx * dxx;
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

void ibs_test_process_perms(uintptr_t* perm_row_start, uint32_t sub_block_ct, double block_tot, double* psbuf, uintptr_t* perm_col_buf, double* perm_results) {
  uintptr_t perm_idx = 0;
  uintptr_t block_size = BITCT;
  double dxx;
  uintptr_t block_pos;
  uintptr_t col_bits;
  uintptr_t ulii;
  uint32_t sub_block_idx;
  do {
    col_bits = *perm_col_buf++;
    if (perm_idx + BITCT > g_perm_ct) {
      block_size = g_perm_ct - perm_idx;
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
  } while (perm_idx < g_perm_ct);
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
  for (row_idx = 0; row_idx < g_thread_start[tidx]; row_idx++) {
    row_uidx = next_set_unsafe(g_pheno_nm, row_uidx);
    row_uidx++;
  }
  for (; row_idx < g_thread_start[tidx + 1]; row_idx++) {
    row_uidx = next_set_unsafe(g_pheno_nm, row_uidx);
    dptr = &(g_dists[(row_uidx * (row_uidx - 1)) / 2]);
    col_idx = 0;
    col_uidx = 0;
    row_set = is_set(g_pheno_c, row_uidx);
    ibs_test_init_col_buf(row_idx, perm_col_buf);
    do {
      if (col_idx + BITCT > row_idx) {
	block_size = row_idx - col_idx;
      } else {
	block_size = BITCT;
      }
      block_tot = fill_psbuf(block_size, dptr, &col_uidx, psptr, &(ssq[row_set]));
      dist_tot += block_tot;
      ibs_test_process_perms(&(g_perm_rows[(col_idx / BITCT) * g_perm_ct]), (block_size + 7) / 8, block_tot, psptr, perm_col_buf, perm_results);
      col_idx += block_size;
    } while (col_idx < row_idx);
    if (!tidx) {
      // technically should change other triangular pct loops to this as well,
      // to guard against int64 overflow with 400m+ people...
      ulii = row_idx * (row_idx + 1);
      if (ulii >= pct_next) {
	if (pct >= 10) {
	  putchar('\b');
	}
	pct = ulii / pct_div;
	printf("\b\b%" PRIuPTR "%%", pct);
	fflush(stdout);
	pct_next = pct_div * (pct + 1);
      }
    }
    row_uidx++;
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
  ibs_test_range((uint32_t)tidx, &(g_perm_col_buf[tidx * perm_ctcl * (CACHELINE / sizeof(intptr_t))]), &(g_perm_results[2 * tidx * perm_ctcld * CACHELINE_DBL]));
  THREAD_RETURN;
}

void incr_dists_i(int32_t* idists, uintptr_t* geno, int32_t tidx) {
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
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    jj = uii * (MULTIPLEX_2DIST / BITCT);
#ifdef __LP64__
    glptr = (__m128i*)geno;
    glptr2 = (__m128i*)(&(geno[jj]));
    lptr = &(g_masks[jj]);
    mcptr_start = (__m128i*)lptr;
    mask_fixed = *lptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *lptr++;
    }
    mptr = (__m128i*)g_masks;
#else
    glptr = geno;
    glptr2 = &(geno[jj]);
    mcptr_start = &(g_masks[jj]);
    mptr = mcptr_start;
    mask_fixed = *mptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *mptr++;
    }
    mptr = g_masks;
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

THREAD_RET_TYPE calc_idist_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_dists_i(&(g_idists[((int64_t)ii * (ii - 1) - (int64_t)jj * (jj - 1)) / 2]), (uintptr_t*)g_geno, (int)tidx);
  THREAD_RETURN;
}

void incr_genome(uint32_t* genome_main, uintptr_t* geno, int32_t tidx) {
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
  __uni16 acc_ibs1;
  __uni16 acc_ibs0;
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
  glptr_end = (__m128i*)(&(geno[g_indiv_ct * (GENOME_MULTIPLEX2 / BITCT)]));
#else
  glptr_end = &(geno[g_indiv_ct * (GENOME_MULTIPLEX2 / BITCT)]);
#endif
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    ujj = uii * (GENOME_MULTIPLEX2 / BITCT);
#ifdef __LP64__
    glptr_fixed = (__m128i*)(&(geno[ujj]));
    glptr = (__m128i*)(&(geno[ujj + (GENOME_MULTIPLEX2 / BITCT)]));
    lptr = &(g_masks[ujj]);
    maskptr = (__m128i*)(&(g_masks[ujj + (GENOME_MULTIPLEX2 / BITCT)]));
    maskptr_fixed = (__m128i*)lptr;
    mask_fixed_test = *lptr++;
    for (ujj = 0; ujj < GENOME_MULTIPLEX2 / BITCT - 1; ujj++) {
      mask_fixed_test &= *lptr++;
    }
#else
    glptr_fixed = &(geno[ujj]);
    glptr = &(geno[ujj + (GENOME_MULTIPLEX2 / BITCT)]);
    maskptr_fixed = &(g_masks[ujj]);
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

          // This PPC test is now the limiting step of --genome, not the IBS
	  // matrix.
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
	    next_ppc_marker_hybrid = ~ZEROLU << (next_ppc_marker_hybrid & (BITCT - 1));
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

THREAD_RET_TYPE calc_genome_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_genome(&(g_genome_main[((int64_t)g_indiv_ct * (ii - jj) - ((int64_t)ii * (ii + 1) - (int64_t)jj * (jj + 1)) / 2) * 5]), (uintptr_t*)g_geno, (int)tidx);
  THREAD_RETURN;
}

void incr_dists(double* dists, uintptr_t* geno, int32_t tidx) {
  uintptr_t* glptr;
  uintptr_t ulii;
  uintptr_t mask_fixed;
  uintptr_t uljj;
  uintptr_t* mptr;
  double* weights1 = &(g_weights[16384]);
#ifdef __LP64__
  double* weights2 = &(g_weights[32768]);
  double* weights3 = &(g_weights[36864]);
  double* weights4 = &(g_weights[40960]);
#endif
  uint32_t uii;
  uint32_t ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    glptr = geno;
    ulii = geno[uii];
    mptr = g_masks;
    mask_fixed = g_masks[uii];
#ifdef __LP64__
    if (mask_fixed == ~ZEROLU) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + g_weights[uljj & 16383];
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + g_weights[uljj & 16383];
	dists++;
      }
    }
#else
    if (mask_fixed == 0x0fffffff) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
	*dists += weights1[uljj >> 14] + g_weights[uljj & 16383];
	dists++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
	*dists += weights1[uljj >> 14] + g_weights[uljj & 16383];
	dists++;
      }
    }
#endif
  }
}

THREAD_RET_TYPE calc_dist_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_dists(&(g_dists[((int64_t)ii * (ii - 1) - (int64_t)jj * (jj - 1)) / 2]), (uintptr_t*)g_geno, (int)tidx);
  THREAD_RETURN;
}

void incr_wt_dist_missing(uint32_t* mtw, int32_t tidx) {
  uintptr_t* glptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    glptr = g_mmasks;
    ulii = g_mmasks[uii];
    if (ulii) {
      for (ujj = 0; ujj < uii; ujj++) {
	uljj = (*glptr++) & ulii;
        while (uljj) {
          mtw[ujj] += g_weights_i[CTZLU(uljj)];
          uljj &= uljj - 1;
        }
      }
    }
    mtw = &(mtw[uii]);
  }
}

THREAD_RET_TYPE calc_distm_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_wt_dist_missing(&(g_missing_tot_weights[((int64_t)ii * (ii - 1) - (int64_t)jj * (jj - 1)) / 2]), (int)tidx);
  THREAD_RETURN;
}

void incr_dists_r(double* dists, uintptr_t* geno, uintptr_t* masks, int32_t tidx, double* weights) {
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
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_dists_r(&(g_rel_dists[((int64_t)ii * (ii - 1) - (int64_t)jj * (jj - 1)) / 2]), (uintptr_t*)g_geno, g_masks, (int32_t)tidx, g_weights);
  THREAD_RETURN;
}

void incr_dists_r_f(float* dists_f, uintptr_t* geno, uintptr_t* masks, int32_t tidx, float* weights_f) {
  uintptr_t* glptr;
  uintptr_t* maskptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t basemask;
  float* weights1 = &(weights_f[32768]);
#ifdef __LP64__
  float* weights2 = &(weights_f[65536]);
  float* weights3 = &(weights_f[98304]);
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
	*dists_f += weights_f[(uint16_t)uljj] + weights1[(uint16_t)(uljj >> 16)] + weights2[(uint16_t)(uljj >> 32)] + weights3[uljj >> 48];
#else
	*dists_f += weights_f[(uint16_t)uljj] + weights1[uljj >> 16];
#endif
	dists_f++;
      }
    } else {
      for (ujj = 0; ujj < uii; ujj++) {
        uljj = ((*glptr++) + ulii) | ((*maskptr++) | basemask);
#ifdef __LP64__
	*dists_f += weights_f[(uint16_t)uljj] + weights1[(uint16_t)(uljj >> 16)] + weights2[(uint16_t)(uljj >> 32)] + weights3[uljj >> 48];
#else
	*dists_f += weights_f[(uint16_t)uljj] + weights1[uljj >> 16];
#endif
	dists_f++;
      }
    }
  }
}

THREAD_RET_TYPE calc_rel_f_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_dists_r_f(&(g_rel_f_dists[((int64_t)ii * (ii - 1) - (int64_t)jj * (jj - 1)) / 2]), (uintptr_t*)g_geno, g_masks, (int)tidx, g_weights_f);
  THREAD_RETURN;
}

void incr_dists_rm(uint32_t* idists, int32_t tidx, uint32_t* thread_start) {
  // count missing intersection, optimized for sparsity
  uintptr_t* mlptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  for (uii = thread_start[tidx]; uii < thread_start[tidx + 1]; uii++) {
    mlptr = g_mmasks;
    ulii = g_mmasks[uii];
    if (ulii) {
      for (ujj = 0; ujj < uii; ujj++) {
        uljj = (*mlptr++) & ulii;
        while (uljj) {
          idists[ujj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[uii]);
  }
}

THREAD_RET_TYPE calc_missing_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  incr_dists_rm(&(g_missing_dbl_excluded[((int64_t)ii * (ii - 1) - (int64_t)jj * (jj - 1)) / 2]), (int)tidx, g_thread_start);
  THREAD_RETURN;
}

void incr_dists_rm_inv(uint32_t* idists, int32_t tidx) {
  // inverted loops for --genome --parallel
  uintptr_t* glptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t indiv_ct_m1 = g_indiv_ct - 1;
  uint32_t uii;
  uint32_t ujj;
  for (uii = g_thread_start[tidx]; uii < g_thread_start[tidx + 1]; uii++) {
    ulii = g_mmasks[uii];
    if (ulii) {
      glptr = &(g_mmasks[uii + 1]);
      // ujj is deliberately biased down by 1
      for (ujj = uii; ujj < indiv_ct_m1; ujj++) {
        uljj = (*glptr++) & ulii;
        while (uljj) {
          idists[ujj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[indiv_ct_m1 - uii - 1]);
  }
}

THREAD_RET_TYPE calc_genomem_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t ii = g_thread_start[tidx];
  int32_t jj = g_thread_start[0];
  // f(0) = 0
  // f(1) = ic - 2
  // f(2) = 2ic - 5
  // f(3) = 3ic - 9
  // ...
  // f(n) = nic - (n+1)(n+2)/2 + 1
  incr_dists_rm_inv(&(g_missing_dbl_excluded[(int64_t)g_indiv_ct * (ii - jj) - ((int64_t)(ii + 1) * (ii + 2) - (int64_t)(jj + 1) * (jj + 2)) / 2]), (int)tidx);
  THREAD_RETURN;
}

void groupdist_jack(int32_t* ibuf, double* returns) {
  int32_t* iptr = ibuf;
  int32_t* jptr = &(ibuf[g_jackknife_d]);
  double neg_tot_uu = 0.0;
  double neg_tot_au = 0.0;
  double neg_tot_aa = 0.0;
  uint32_t neg_a = 0;
  uint32_t neg_u = 0;
  double* dptr;
  int32_t* iptr2;
  uintptr_t indiv_idx;
  uint32_t ii2;
  while (iptr < jptr) {
    dptr = &(g_jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_GROUPDIST]);
    neg_tot_uu += *dptr++;
    neg_tot_au += *dptr++;
    neg_tot_aa += *dptr++;
  }
  iptr = ibuf;
  while (iptr < jptr) {
    indiv_idx = *iptr;
    iptr2 = ibuf;
    dptr = &(g_dists[(indiv_idx * (indiv_idx - 1)) / 2]);
    if (is_set(g_pheno_c, indiv_idx)) {
      neg_a++;
      while (iptr2 < iptr) {
	ii2 = *iptr2++;
	if (is_set(g_pheno_c, ii2)) {
	  neg_tot_aa -= dptr[ii2];
	} else {
	  neg_tot_au -= dptr[ii2];
	}
      }
    } else {
      neg_u++;
      while (iptr2 < iptr) {
	ii2 = *iptr2++;
	if (is_set(g_pheno_c, ii2)) {
	  neg_tot_au -= dptr[ii2];
	} else {
	  neg_tot_uu -= dptr[ii2];
	}
      }
    }
    iptr++;
  }
  returns[0] = (g_reg_tot_x - neg_tot_aa) / (double)(((intptr_t)(g_case_ct - neg_a) * (g_case_ct - neg_a - 1)) / 2);
  returns[1] = (g_reg_tot_xy - neg_tot_au) / (double)((intptr_t)(g_case_ct - neg_a) * (g_ctrl_ct - neg_u));
  returns[2] = (g_reg_tot_y - neg_tot_uu) / (double)(((intptr_t)(g_ctrl_ct - neg_u) * (g_ctrl_ct - neg_u - 1)) / 2);
}

void small_remap(int32_t* ibuf, uint32_t ct, uint32_t dd) {
  int32_t* ibuf_end = &(ibuf[dd]);
  int32_t missings = 0;
  int32_t curpos = 0;
  do {
    if (!is_set(g_pheno_nm, curpos)) {
      missings++;
    } else if (*ibuf == curpos - missings) {
      *ibuf++ = curpos;
    }
    curpos++;
  } while (ibuf < ibuf_end);
}

THREAD_RET_TYPE groupdist_jack_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t* ibuf = (int32_t*)(&(g_geno[tidx * CACHEALIGN(g_case_ct + g_ctrl_ct + (g_jackknife_d + 1) * sizeof(int32_t))]));
  unsigned char* cbuf = &(g_geno[tidx * CACHEALIGN(g_case_ct + g_ctrl_ct + (g_jackknife_d + 1) * sizeof(int32_t)) + (g_jackknife_d + 1) * sizeof(int32_t)]);
  uint64_t ullii;
  uint64_t ulljj = g_jackknife_iters / 100;
  double returns[3];
  double results[9];
  double new_old_diff[3];
  fill_double_zero(results, 9);
  for (ullii = 0; ullii < g_jackknife_iters; ullii++) {
    pick_d_small(cbuf, ibuf, g_case_ct + g_ctrl_ct, g_jackknife_d);
    if (g_case_ct + g_ctrl_ct < g_indiv_ct) {
      small_remap(ibuf, g_case_ct + g_ctrl_ct, g_jackknife_d);
    }
    groupdist_jack(ibuf, returns);
    if (ullii > 0) {
      new_old_diff[0] = returns[0] - results[0];
      new_old_diff[1] = returns[1] - results[1];
      new_old_diff[2] = returns[2] - results[2];
      results[0] += new_old_diff[0] / (ullii + 1); // AA mean
      results[1] += new_old_diff[1] / (ullii + 1); // AU mean
      results[2] += new_old_diff[2] / (ullii + 1); // UU mean
      results[3] += (returns[0] - results[0]) * new_old_diff[0]; // AA var
      results[4] += (returns[1] - results[1]) * new_old_diff[1]; // AU var
      results[5] += (returns[2] - results[2]) * new_old_diff[2]; // UU var
      results[6] += (returns[0] - results[0]) * new_old_diff[1]; // AA-AU cov
      results[7] += (returns[0] - results[0]) * new_old_diff[2]; // AA-UU cov
      results[8] += (returns[1] - results[1]) * new_old_diff[2]; // AU-UU cov
    } else {
      results[0] += returns[0];
      results[1] += returns[1];
      results[2] += returns[2];
    }
    if ((!tidx) && (ullii >= ulljj)) {
      ulljj = (ullii * 100) / g_jackknife_iters;
      printf("\r%" PRIu64 "%%", ulljj);
      fflush(stdout);
      ulljj = ((ulljj + 1) * g_jackknife_iters) / 100;
    }
  }
  // don't write until end, to avoid false sharing
  for (ullii = 0; ullii < 9; ullii++) {
    g_calc_result[tidx][ullii] = results[ullii];
  }
  THREAD_RETURN;
}

double regress_rel_jack(int32_t* ibuf, double* ret2_ptr) {
  int32_t* iptr = ibuf;
  int32_t* jptr = &(ibuf[g_jackknife_d]);
  uint32_t uii;
  int32_t jj;
  int32_t kk;
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
  while (iptr < jptr) {
    dptr2 = &(g_jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_REL]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  iptr = ibuf;
  for (uii = 1; uii < g_jackknife_d; uii++) {
    jj = *(++iptr);
    dxx1 = g_pheno_packed[jj];
    jptr = ibuf;
    dptr = &(g_rel_dists[((intptr_t)jj * (jj - 1)) / 2]);
    while (jptr < iptr) {
      kk = *jptr++;
      dxx = (dxx1 + g_pheno_packed[kk]) * 0.5;
      dyy = dptr[kk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = g_reg_tot_y - neg_tot_y;
  dyy = g_indiv_ct - g_jackknife_d;
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((g_reg_tot_xy - neg_tot_xy) - dxx * (g_reg_tot_x - neg_tot_x) / dyy) / ((g_reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = g_reg_tot_x - neg_tot_x;
  return ((g_reg_tot_xy - neg_tot_xy) - dxx * (g_reg_tot_y - neg_tot_y) / dyy) / ((g_reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

THREAD_RET_TYPE regress_rel_jack_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  int32_t* ibuf = (int32_t*)(&(g_geno[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t))]));
  unsigned char* cbuf = &(g_geno[tidx * CACHEALIGN(g_indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t)) + (g_jackknife_d + 1) * sizeof(int32_t)]);
  uint64_t ullii;
  uint64_t ulljj = g_jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  double dxx;
  double ret2;
  for (ullii = 0; ullii < g_jackknife_iters; ullii++) {
    pick_d_small(cbuf, ibuf, g_indiv_ct, g_jackknife_d);
    dxx = regress_rel_jack(ibuf, &ret2);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ullii >= ulljj)) {
      ulljj = (ullii * 100) / g_jackknife_iters;
      printf("\r%" PRIu64 "%%", ulljj);
      fflush(stdout);
      ulljj = ((ulljj + 1) * g_jackknife_iters) / 100;
    }
  }
  g_calc_result[tidx][0] = sum;
  g_calc_result[tidx][1] = sum_sq;
  g_calc_result[tidx][2] = sum2;
  g_calc_result[tidx][3] = sum2_sq;
  THREAD_RETURN;
}

int32_t regress_rel_main(uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t regress_rel_iters, int32_t regress_rel_d, pthread_t* threads, double* pheno_d) {
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
  g_pheno_packed = (double*)wkspace_alloc(indiv_ct * sizeof(double));
  if (!g_pheno_packed) {
    return RET_NOMEM;
  }
  collapse_copy_phenod(g_pheno_packed, pheno_d, indiv_exclude, indiv_ct);
  print_pheno_stdev(g_pheno_packed, indiv_ct);
  trimatrix_size = ((uintptr_t)indiv_ct * (indiv_ct - 1)) / 2;
  g_reg_tot_xy = 0.0;
  g_reg_tot_x = 0.0;
  g_reg_tot_y = 0.0;
  g_reg_tot_xx = 0.0;
  g_reg_tot_yy = 0.0;
  rel_ptr = g_rel_dists;
  pheno_ptr = g_pheno_packed;
  g_jackknife_precomp = (double*)wkspace_alloc(indiv_ct * JACKKNIFE_VALS_REL * sizeof(double));
  if (!g_jackknife_precomp) {
    return RET_NOMEM;
  }
  fill_double_zero(g_jackknife_precomp, indiv_ct * JACKKNIFE_VALS_REL);
  for (uii = 1; uii < indiv_ct; uii++) {
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
      g_reg_tot_xy += dxxyy;
      jp_fixed_ptr[0] += dxxyy;
      *jp_moving_ptr += dxxyy;
      jp_moving_ptr++;
      g_reg_tot_x += dxx;
      jp_fixed_ptr[1] += dxx;
      *jp_moving_ptr += dxx;
      jp_moving_ptr++;
      g_reg_tot_y += dyy;
      jp_fixed_ptr[2] += dyy;
      *jp_moving_ptr += dyy;
      jp_moving_ptr++;
      g_reg_tot_xx += dxxsq;
      jp_fixed_ptr[3] += dxxsq;
      *jp_moving_ptr += dxxsq;
      jp_moving_ptr++;
      g_reg_tot_yy += dyysq;
      jp_fixed_ptr[4] += dyysq;
      *jp_moving_ptr += dyysq;
      jp_moving_ptr++;
    }
  }
  trimatrix_size_recip = 1.0 / (double)trimatrix_size;
  sprintf(logbuf, "Regression slope (y = genomic relationship, x = avg phenotype): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y * trimatrix_size_recip) / (g_reg_tot_xx - g_reg_tot_x * g_reg_tot_x * trimatrix_size_recip));
  logprintb();
  sprintf(logbuf, "                 (y = avg phenotype, x = genomic relationship): %g\n", (g_reg_tot_xy - g_reg_tot_x * g_reg_tot_y * trimatrix_size_recip) / (g_reg_tot_yy - g_reg_tot_y * g_reg_tot_y * trimatrix_size_recip));
  logprintb();
  g_jackknife_iters = (regress_rel_iters + g_thread_ct - 1) / g_thread_ct;
  if (regress_rel_d) {
    g_jackknife_d = regress_rel_d;
  } else {
    g_jackknife_d = set_default_jackknife_d(indiv_ct);
  }
  g_geno = wkspace_alloc(g_thread_ct * CACHEALIGN(indiv_ct + (g_jackknife_d + 1) * sizeof(int32_t)));
  if (!g_geno) {
    return RET_NOMEM;
  }
  if (spawn_threads(threads, &regress_rel_jack_thread, g_thread_ct)) {
    logprint(errstr_thread_create);
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
  putchar('\r');
  sprintf(logbuf, "Jackknife s.e. (y = genomic relationship): %g\n", sqrt((indiv_ct / (double)g_jackknife_d) * (dxxsq - dxx * dxx / (double)ulii) / ((double)ulii - 1)));
  logprintb();
  sprintf(logbuf, "               (y = phenotype): %g\n", sqrt((indiv_ct / (double)g_jackknife_d) * (dyysq - dyy * dyy / (double)ulii) / ((double)ulii - 1)));
  logprintb();
  return 0;
}

// Replaces matrix[][] with mult_val * matrix[][] + add_val * I.
// Multithreading doesn't help here.
void matrix_const_mult_add(double* matrix, double mult_val, double add_val) {
  uint32_t uii;
  uint32_t loop_end = g_indiv_ct - 1;
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
    if (ujj < g_indiv_ct) {
      *dptr *= mult_val;
      dptr++;
    }
#else
    for (ujj = 0; ujj < g_indiv_ct; ujj++) {
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
void matrix_row_sum_ur(double* sums, double* matrix) {
  uintptr_t indiv_idx;
  double* dptr;
  double acc;
  double* sptr_end;
  double* sptr;
  fill_double_zero(sums, g_indiv_ct);
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    dptr = &(matrix[indiv_idx * g_indiv_ct]);
    acc = 0.0;
    sptr_end = &(sums[indiv_idx]);
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
// one-trait REML via EM.
//
// wkbase is assumed to have space for three cache-aligned
// indiv_ct * indiv_ct double matrices plus three more rows.  The unpacked
// relationship matrix is stored in the SECOND slot.
void reml_em_one_trait(double* wkbase, double* pheno, double* covg_ref, double* covr_ref, double tol, int32_t strict) {
  double ll_change;
  int64_t mat_offset = g_indiv_ct;
  double* rel_dists;
  MATRIX_INVERT_BUF1_TYPE* irow;
  __CLPK_integer lwork;
  double* row;
  double* row2;
  double* work;
  double* dptr;
  double* dptr2;
  double* matrix_pvg;
  __CLPK_integer indiv_ct_li = g_indiv_ct;
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
  double indiv_ct_d = 1 / (double)g_indiv_ct;
  uintptr_t indiv_idx;
  int32_t jj;
#if _WIN32
  char blas_char;
  int32_t indiv_ct_i32 = g_indiv_ct;
#endif
  mat_offset = CACHEALIGN_DBL(mat_offset * mat_offset);
  rel_dists = &(wkbase[mat_offset]);
  row = &(wkbase[mat_offset * 3]);
  irow = (MATRIX_INVERT_BUF1_TYPE*)row;
  row2 = &(row[g_indiv_ct]);
  work = &(wkbase[mat_offset * 2]);
  lwork = mat_offset;
  matrix_pvg = work;
  if (!lwork) {
    lwork = CACHELINE_DBL;
  }
  fill_double_zero(matrix_pvg, mat_offset);
  fill_double_zero(row2, g_indiv_ct);
  do {
    memcpy(wkbase, rel_dists, mat_offset * sizeof(double));
    matrix_const_mult_add(wkbase, *covg_ref, *covr_ref);
    invert_matrix(indiv_ct_li, wkbase, irow, work);
    matrix_row_sum_ur(row, wkbase);
    dxx = 0.0;
    dptr = row;
    dptr2 = &(row[g_indiv_ct]);
    while (dptr < dptr2) {
      dxx += *dptr++;
    }
    dxx = -1 / dxx;
#if _WIN32
    jj = 1;
    dger_(&indiv_ct_i32, &indiv_ct_i32, &dxx, row, &jj, row, &jj, wkbase, &indiv_ct_i32);
    // todo: test dsymm
    dyy = 1.0;
    dzz = 0.0;
    blas_char = 'N';
    dgemm_(&blas_char, &blas_char, &indiv_ct_i32, &indiv_ct_i32, &indiv_ct_i32, &dyy, wkbase, &indiv_ct_i32, rel_dists, &indiv_ct_i32, &dzz, matrix_pvg, &indiv_ct_i32);
    dlg = 0.0;
    dle = 0.0;
    jj = g_indiv_ct + 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      dlg -= matrix_pvg[indiv_idx * jj];
      dle -= wkbase[indiv_idx * jj];
    }
    blas_char = 'U';
    jj = 1;
    dsymv_(&blas_char, &indiv_ct_i32, &dyy, wkbase, &indiv_ct_i32, pheno, &jj, &dzz, row2, &jj);
    dsymv_(&blas_char, &indiv_ct_i32, &dyy, matrix_pvg, &indiv_ct_i32, row2, &jj, &dzz, row, &jj);
    dlg += ddot_(&indiv_ct_i32, pheno, &jj, row, &jj);
    dsymv_(&blas_char, &indiv_ct_i32, &dyy, wkbase, &indiv_ct_i32, row2, &jj, &dzz, row, &jj);
    dle += ddot_(&indiv_ct_i32, pheno, &jj, row, &jj);
#else
    cblas_dger(CblasColMajor, g_indiv_ct, g_indiv_ct, dxx, row, 1, row, 1, wkbase, g_indiv_ct);
    // unfortunately, cblas_dsymm is much worse than cblas_dgemm on OS X
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, g_indiv_ct, g_indiv_ct, g_indiv_ct, 1.0, wkbase, g_indiv_ct, rel_dists, g_indiv_ct, 0.0, matrix_pvg, g_indiv_ct);
    dlg = 0.0;
    dle = 0.0;
    jj = g_indiv_ct + 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      dlg -= matrix_pvg[indiv_idx * jj];
      dle -= wkbase[indiv_idx * jj];
    }
    cblas_dsymv(CblasColMajor, CblasUpper, g_indiv_ct, 1.0, wkbase, g_indiv_ct, pheno, 1, 0.0, row2, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, g_indiv_ct, 1.0, matrix_pvg, g_indiv_ct, row2, 1, 0.0, row, 1);
    dlg += cblas_ddot(g_indiv_ct, pheno, 1, row, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, g_indiv_ct, 1.0, wkbase, g_indiv_ct, row2, 1, 0.0, row, 1);
    dle += cblas_ddot(g_indiv_ct, pheno, 1, row, 1);
#endif
    covg_last_change = covg_cur_change;
    covr_last_change = covr_cur_change;
    covg_cur_change = (*covg_ref) * (*covg_ref) * dlg * indiv_ct_d;
    covr_cur_change = (*covr_ref) * (*covr_ref) * dle * indiv_ct_d;
    if (strict) {
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
  putchar('\n');
  sprintf(logbuf, "covg: %g  covr: %g\n", *covg_ref, *covr_ref);
  logstr(logbuf);
}

int32_t calc_unrelated_herit(uint64_t calculation_type, int32_t ibc_type, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, double* pheno_d, double* rel_ibc, double unrelated_herit_covg, double unrelated_herit_covr, double unrelated_herit_tol) {
  double dxx = 0.0;
  double dyy = 0.0;
  double dzz;
  uintptr_t ulii;
  uintptr_t uljj;
  double* dist_ptr;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  g_missing_dbl_excluded = NULL;
  ulii = g_indiv_ct;
  ulii = CACHEALIGN_DBL(ulii * ulii);
  dptr4 = &(g_rel_dists[ulii]);
  ulii = ulii * 3 + CACHEALIGN_DBL(g_indiv_ct) * 3;
  wkspace_reset(g_rel_dists);
  g_rel_dists = (double*)wkspace_alloc(ulii * sizeof(double));
  if (!g_rel_dists) {
    return RET_NOMEM;
  }
  dptr2 = &(g_rel_dists[ulii - CACHEALIGN_DBL(g_indiv_ct)]);
  collapse_copy_phenod(dptr2, pheno_d, indiv_exclude, unfiltered_indiv_ct);
  dptr3 = dptr2;
  dist_ptr = &(dptr2[g_indiv_ct]);
  while (dptr3 < dist_ptr) {
    dzz = *dptr3++;
    dxx += dzz;
    dyy += dzz * dzz;
  }
  dxx /= (double)g_indiv_ct;
  dxx = 1 / sqrt((dyy / (double)g_indiv_ct) - dxx * dxx);
  dptr3 = dptr2;
  while (dptr3 < dist_ptr) {
    *dptr3 *= dxx;
    dptr3++;
  }
  if (calculation_type & CALC_IBC) {
    dptr3 = &(rel_ibc[ibc_type * g_indiv_ct]);
  } else {
    dptr3 = rel_ibc;
  }
  for (ulii = 0; ulii < g_indiv_ct; ulii++) {
    memcpy(&(dptr4[ulii * g_indiv_ct]), &(g_rel_dists[(ulii * (ulii - 1)) / 2]), ulii * sizeof(double));
    dptr4[ulii * (g_indiv_ct + 1)] = *dptr3++;
    for (uljj = ulii + 1; uljj < g_indiv_ct; uljj++) {
      dptr4[ulii * g_indiv_ct + uljj] = g_rel_dists[(uljj * (uljj - 1)) / 2 + ulii];
    }
  }
  reml_em_one_trait(g_rel_dists, dptr2, &unrelated_herit_covg, &unrelated_herit_covr, unrelated_herit_tol, calculation_type & CALC_UNRELATED_HERITABILITY_STRICT);
  sprintf(logbuf, "h^2 estimate: %g\n", unrelated_herit_covg);
  logprintb();
  return 0;
}
#endif

/*
double get_dmedian_i(int32_t* sorted_arr, int32_t len) {
  if (len) {
    if (len % 2) {
      return (double)sorted_arr[len / 2];
    } else {
      return ((double)sorted_arr[len / 2] + (double)sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
}
*/

double get_dmedian(double* sorted_arr, int32_t len) {
  if (len) {
    if (len % 2) {
      return sorted_arr[len / 2];
    } else {
      return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
}

int32_t ibs_test_calc(pthread_t* threads, uint64_t calculation_type, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t perm_ct, uintptr_t pheno_nm_ct, uintptr_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ctl  = (g_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t pheno_nm_ctl = (pheno_nm_ct + (BITCT - 1)) / BITCT;
  uintptr_t perm_ctcl = (perm_ct + (CACHELINE * 8)) / (CACHELINE * 8);
  uintptr_t perm_ctclm = perm_ctcl * (CACHELINE / sizeof(intptr_t));
  uintptr_t perm_ctcld = (perm_ct + CACHELINE_DBL) / CACHELINE_DBL;
  uintptr_t perm_ctcldm = perm_ctcld * CACHELINE_DBL;
  uintptr_t case_ct = pheno_nm_ct - pheno_ctrl_ct;
  uint32_t tidx = 1;
  int32_t retval = 0;
  int32_t perm_test[6];
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
  g_load_dists = (calculation_type & CALC_LOAD_DISTANCES)? 1 : 0;
  perm_ct += 1; // first permutation = original config
  if (pheno_ctrl_ct < 2) {
    logprint("Skipping --ibs-test: Too few controls (minimum 2).\n");
    goto ibs_test_calc_ret_1;
  } else if (case_ct < 2) {
    logprint("Skipping --ibs-test: Too few cases (minimum 2).\n");
    goto ibs_test_calc_ret_1;
  }
  for (ulii = 0; ulii < 6; ulii++) {
    perm_test[ulii] = 0;
  }
  ctrl_ctrl_ct = (pheno_ctrl_ct * (pheno_ctrl_ct - 1)) / 2;
  ctrl_case_ct = pheno_ctrl_ct * case_ct;
  case_case_ct = (case_ct * (case_ct - 1)) / 2;
  g_perm_ct = perm_ct;
  // g_pheno_nm and g_pheno_c should be NULL
  if (wkspace_alloc_ul_checked(&g_pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto ibs_test_calc_ret_NOMEM;
  }
  memcpy(g_pheno_nm, pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
  memcpy(g_pheno_c, pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t));
  collapse_bitarr(g_pheno_nm, indiv_exclude, unfiltered_indiv_ct);
  collapse_bitarr(g_pheno_c, indiv_exclude, unfiltered_indiv_ct);
  if (wkspace_alloc_d_checked(&g_ibs_test_partial_sums, g_thread_ct * 32 * BITCT * sizeof(double)) ||
      wkspace_alloc_ul_checked(&g_perm_rows, perm_ct * pheno_nm_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_perm_col_buf, perm_ctclm * sizeof(intptr_t) * g_thread_ct) ||
      wkspace_alloc_d_checked(&g_perm_results, 2 * perm_ctcldm * sizeof(double) * g_thread_ct)) {
    goto ibs_test_calc_ret_NOMEM;
  }
  fill_double_zero(g_perm_results, 2 * perm_ctcldm * g_thread_ct);

  // first permutation = original
  memcpy(g_perm_rows, g_pheno_c, indiv_ctl * sizeof(intptr_t));
  collapse_bitarr_incl(g_perm_rows, g_pheno_nm, g_indiv_ct);
  for (ulii = pheno_nm_ctl - 1; ulii; ulii--) {
    g_perm_rows[ulii * perm_ct] = g_perm_rows[ulii];
  }

  printf("--ibs-test (%" PRIuPTR " permutations): [generating permutations]", perm_ct - 1);
  fflush(stdout);
  // minor todo: multithread this
  generate_perm1_interleaved(pheno_nm_ct, case_ct, 1, perm_ct, g_perm_rows);
  fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                       \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b0%", stdout);
  fflush(stdout);
  for (ulii = 0; ulii < pheno_nm_ct; ulii++) {
    uljj += ((g_perm_rows[((ulii / BITCT) * perm_ct)] >> (ulii & (BITCT - 1))) & 1);
  }
  triangle_fill(g_thread_start, pheno_nm_ct, g_thread_ct, 0, 1, 1, 1);
  if (spawn_threads(threads, &ibs_test_thread, g_thread_ct)) {
    goto ibs_test_calc_ret_THREAD_CREATE_FAIL;
  }
  ulii = 0;
  fflush(stdout);
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
    rvptr1 = (__m128d*)g_perm_results;
    rvptr2 = (__m128d*)(&(g_perm_results[2 * perm_ctcldm * tidx]));
    for (perm_idx = 0; perm_idx < perm_ct; perm_idx++) {
      *rvptr1 = _mm_add_pd(*rvptr1, *rvptr2++);
      rvptr1++;
    }
#else
    rptr1 = g_perm_results;
    rptr2 = &(g_perm_results[2 * perm_ctcldm * tidx]);
    for (perm_idx = 0; perm_idx < perm_ct; perm_idx++) {
      *rptr1++ += *rptr2++;
      *rptr1++ += *rptr2++;
    }
#endif
  }
  ctrl_ctrl_tot = g_perm_results[0];
  ctrl_case_tot = g_perm_results[1];
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
    ctrl_ctrl_tot1 = g_perm_results[ulii * 2];
    ctrl_case_tot1 = g_perm_results[ulii * 2 + 1];
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
  sprintf(logbuf, "  Between-group IBS (mean, SD)   = %g, %g\n", ctrl_case_mean, sqrt(ctrl_case_var / (ctrl_case_ct - 1)));
  logprintb();
  sprintf(logbuf, "  In-group (case) IBS (mean, SD) = %g, %g\n", case_case_mean, sqrt(case_case_var / (case_case_ct - 1)));
  logprintb();
  sprintf(logbuf, "  In-group (ctrl) IBS (mean, SD) = %g, %g\n", ctrl_ctrl_mean, sqrt(ctrl_ctrl_var / (ctrl_ctrl_ct - 1)));
  logprintb();
  sprintf(logbuf, "  Approximate proportion of variance between group = %g\n", between_ssq / total_ssq);
  logprintb();
  perm_ct_recip = 1.0 / ((double)perm_ct);
  fputs("  IBS group-difference empirical p-values:\n", stdout);
  sprintf(logbuf, "     T1: Case/control less similar                p = %g\n", (perm_test[0]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T2: Case/control more similar                p = %g\n\n", (perm_ct - perm_test[0]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T3: Case/case less similar than ctrl/ctrl    p = %g\n", (perm_test[1]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T4: Case/case more similar than ctrl/ctrl    p = %g\n\n", (perm_ct - perm_test[1]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T5: Case/case less similar                   p = %g\n", (perm_test[2]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T6: Case/case more similar                   p = %g\n\n", (perm_ct - perm_test[2]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T7: Control/control less similar             p = %g\n", (perm_test[3]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T8: Control/control more similar             p = %g\n\n", (perm_ct - perm_test[3]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "     T9: Case/case less similar than case/ctrl    p = %g\n", (perm_test[4]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "    T10: Case/case more similar than case/ctrl    p = %g\n\n", (perm_ct - perm_test[4]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "    T11: Ctrl/ctrl less similar than case/ctrl    p = %g\n", (perm_test[5]) * perm_ct_recip);
  logprintb();
  sprintf(logbuf, "    T12: Ctrl/ctrl more similar than case/ctrl    p = %g\n", (perm_ct - perm_test[5]) * perm_ct_recip);
  logprintb();


  while (0) {
  ibs_test_calc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ibs_test_calc_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 ibs_test_calc_ret_1:
  wkspace_reset(wkspace_mark);
  g_pheno_nm = NULL;
  g_pheno_c = NULL;
  return retval;
}

int32_t groupdist_calc(pthread_t* threads, uint32_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t groupdist_iters, uint32_t groupdist_d, uint32_t pheno_nm_ct, uint32_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  double* dist_ptr = g_dists;
  double dhh_ssq = 0.0;
  double dhl_ssq = 0.0;
  double dll_ssq = 0.0;
  int32_t retval = 0;
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
  uint32_t indiv_idx;
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
    logprint("Skipping --groupdist: Too few controls (minimum 2).\n");
    goto groupdist_calc_ret_1;
  }
  g_case_ct = pheno_nm_ct - pheno_ctrl_ct;
  if (g_case_ct < 2) {
    logprint("Skipping --groupdist: Too few cases (minimum 2).\n");
    goto groupdist_calc_ret_1;
  }
  g_ctrl_ct = pheno_ctrl_ct;
  // g_pheno_nm and g_pheno_c should be NULL
  if (wkspace_alloc_ul_checked(&g_pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    goto groupdist_calc_ret_NOMEM;
  }
  memcpy(g_pheno_nm, pheno_nm, unfiltered_indiv_ctl * sizeof(intptr_t));
  memcpy(g_pheno_c, pheno_c, unfiltered_indiv_ctl * sizeof(intptr_t));
  collapse_bitarr(g_pheno_nm, indiv_exclude, unfiltered_indiv_ct);
  collapse_bitarr(g_pheno_c, indiv_exclude, unfiltered_indiv_ct);
  zero_trailing_bits(g_pheno_nm, g_indiv_ct);
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
  if (wkspace_alloc_d_checked(&ll_pool, ll_size * sizeof(double)) ||
      wkspace_alloc_d_checked(&lh_pool, lh_size * sizeof(double)) ||
      wkspace_alloc_d_checked(&hh_pool, hh_size * sizeof(double)) ||
      wkspace_alloc_uc_checked(&g_geno, g_thread_ct * CACHEALIGN(g_case_ct + g_ctrl_ct + (g_jackknife_d + 1) * sizeof(int32_t)))) {
    goto groupdist_calc_ret_NOMEM;
  }
  ll_poolp = ll_pool;
  lh_poolp = lh_pool;
  hh_poolp = hh_pool;
  for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
    if (is_set(g_pheno_nm, indiv_idx)) {
      if (is_set(g_pheno_c, indiv_idx)) {
	for (uii = 0; uii < indiv_idx; uii++) {
	  if (is_set(g_pheno_nm, uii)) {
	    dxx = *dist_ptr;
	    if (is_set(g_pheno_c, uii)) {
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
	for (uii = 0; uii < indiv_idx; uii++) {
	  if (is_set(g_pheno_nm, uii)) {
	    dxx = *dist_ptr;
	    if (is_set(g_pheno_c, uii)) {
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
      dist_ptr += indiv_idx;
    }
  }
#ifdef __cplusplus
  // std::sort is faster than qsort for basic types.  See e.g. Anders
  // Kaseorg's answer to
  // http://www.quora.com/Software-Engineering/Generally-how-much-faster-is-C-compared-to-C++
  std::sort(ll_pool, &(ll_pool[ll_size]));
  std::sort(lh_pool, &(lh_pool[lh_size]));
  std::sort(hh_pool, &(hh_pool[hh_size]));
#else
  qsort(ll_pool, ll_size, sizeof(double), double_cmp);
  qsort(lh_pool, lh_size, sizeof(double), double_cmp);
  qsort(hh_pool, hh_size, sizeof(double), double_cmp);
#endif
  ll_med = get_dmedian(ll_pool, ll_size);
  lh_med = get_dmedian(lh_pool, lh_size);
  hh_med = get_dmedian(hh_pool, hh_size);
  logprint("Case/control distance analysis:\n");
  if (g_case_ct < 2) {
    dxx = 0.0;
    dhh_sd = 0.0;
  } else {
    dww = (double)(((uintptr_t)g_case_ct * (g_case_ct - 1)) / 2);
    dxx = g_reg_tot_x / dww;
    dhh_sd = sqrt((dhh_ssq / dww - dxx * dxx) / (dww - 1.0));
  }
  if (!(g_case_ct * g_ctrl_ct)) {
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
  sprintf(logbuf, "  Mean (sd), median dists between 2x affected     : %g (%g), %g\n", dxx, dhh_sd, hh_med);
  logprintb();
  sprintf(logbuf, "  Mean (sd), median dists between aff. and unaff. : %g (%g), %g\n", dyy, dhl_sd, lh_med);
  logprintb();
  sprintf(logbuf, "  Mean (sd), median dists between 2x unaffected   : %g (%g), %g\n\n", dzz, dll_sd, ll_med);
  logprintb();
  if (2 * g_jackknife_d >= (g_case_ct + g_ctrl_ct)) {
    logprint("Delete-d jackknife skipped because d is too large.\n");
  } else {
    if (wkspace_alloc_d_checked(&g_jackknife_precomp, g_indiv_ct * JACKKNIFE_VALS_GROUPDIST)) {
      goto groupdist_calc_ret_NOMEM;
    }
    fill_double_zero(g_jackknife_precomp, g_indiv_ct * JACKKNIFE_VALS_GROUPDIST);
    // to precompute:
    // tot_uu, tot_au, tot_aa
    dist_ptr = g_dists;
    for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
      if (is_set(g_pheno_nm, indiv_idx)) {
	uii = 0;
	is_case = is_set(g_pheno_c, indiv_idx);
	dyy = 0;
	dzz = 0;
	do {
	  if (is_set(g_pheno_nm, uii)) {
	    dxx = dist_ptr[uii];
	    if (is_set(g_pheno_c, uii)) {
	      g_jackknife_precomp[(uii * JACKKNIFE_VALS_GROUPDIST) + is_case + 1] += dxx;
	      dzz += dxx;
	    } else {
	      g_jackknife_precomp[(uii * JACKKNIFE_VALS_GROUPDIST) + is_case] += dxx;
	      dyy += dxx;
	    }
	  }
	} while ((++uii) < indiv_idx);
	g_jackknife_precomp[(indiv_idx * JACKKNIFE_VALS_GROUPDIST) + is_case] += dyy;
	g_jackknife_precomp[(indiv_idx * JACKKNIFE_VALS_GROUPDIST) + is_case + 1] += dzz;
      }
      dist_ptr = &(dist_ptr[indiv_idx]);
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
    dxx = 1.0 / g_thread_ct;
    g_calc_result[0][0] *= dxx;
    g_calc_result[0][1] *= dxx;
    g_calc_result[0][2] *= dxx;
    dxx /= (g_jackknife_iters - 1) * g_thread_ct;
    for (uii = 3; uii < 9; uii++) {
      g_calc_result[0][uii] *= dxx;
    }
    putchar('\r');
    sprintf(logbuf, "  AA mean - AU mean avg difference (s.e.): %g (%g)\n", g_calc_result[0][0] - g_calc_result[0][1], sqrt(((g_case_ct + g_ctrl_ct) / ((double)g_jackknife_d)) * (g_calc_result[0][3] + g_calc_result[0][4] - 2 * g_calc_result[0][6])));
    logprintb();
    sprintf(logbuf, "  AA mean - UU mean avg difference (s.e.): %g (%g)\n", g_calc_result[0][0] - g_calc_result[0][2], sqrt(((g_case_ct + g_ctrl_ct) / ((double)g_jackknife_d)) * (g_calc_result[0][3] + g_calc_result[0][5] - 2 * g_calc_result[0][7])));
    logprintb();
    sprintf(logbuf, "  AU mean - UU mean avg difference (s.e.): %g (%g)\n", g_calc_result[0][1] - g_calc_result[0][2], sqrt(((g_case_ct + g_ctrl_ct) / ((double)g_jackknife_d)) * (g_calc_result[0][4] + g_calc_result[0][5] - 2 * g_calc_result[0][8])));
    logprintb();
  }
  while (0) {
  groupdist_calc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  groupdist_calc_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 groupdist_calc_ret_1:
  wkspace_reset(wkspace_mark);
  g_pheno_nm = NULL;
  g_pheno_c = NULL;
  return retval;
}

void normalize_phenos(double* new_phenos, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t sex_exclude) {
  uint32_t incl_males = sex_exclude & 1;
  uint32_t incl_females = sex_exclude & 2;
  uint32_t incl_unknown = sex_exclude & 4;
  uintptr_t indiv_uidx = ~ZEROLU;
  double pheno_tot = 0.0;
  double pheno_sq_tot = 0.0;
  uint32_t pheno_ct = 0;
  uintptr_t indiv_idx = ~ZEROLU;
  double dxx;
  double mean;
  double stdev_recip;

  while ((++indiv_idx) < indiv_ct) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx + 1);
    if (is_set(sex_nm, indiv_uidx)) {
      if (is_set(sex_male, indiv_uidx)) {
	if (!incl_males) {
	  continue;
	}
      } else if (!incl_females) {
	continue;
      }
    } else if (!incl_unknown) {
      continue;
    }
    dxx = new_phenos[indiv_idx];
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
  indiv_uidx = ~ZEROLU;
  indiv_idx = ~ZEROLU;
  while ((++indiv_idx) < indiv_ct) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx + 1);
    if (is_set(sex_nm, indiv_uidx)) {
      if (is_set(sex_male, indiv_uidx)) {
	if (!incl_males) {
	  continue;
	}
      } else if (!incl_females) {
	continue;
      }
    } else if (!incl_unknown) {
      continue;
    }
    new_phenos[indiv_idx] = (new_phenos[indiv_idx] - mean) * stdev_recip;
  }
}

int32_t ld_process_load(unsigned char* loadbuf, uintptr_t* geno_buf, uintptr_t* mask_buf, uintptr_t* missing_buf, double* marker_stdev_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uint32_t indiv_ctl, int32_t indiv_trail_ct, uint32_t is_haploid, uint32_t is_x, uintptr_t* sex_male) {
  uintptr_t unfiltered_idx = 0;
  uint32_t write_offset;
  int32_t sloop_max;
  int32_t write_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t new_geno;
  uintptr_t new_mask;
  uintptr_t new_missing;
  int32_t missing_ct = 0;
  int32_t sq_sum = 0;
  int32_t sum = -indiv_ct;
  double non_missing_recip;
  for (write_offset = 0; write_offset < indiv_ctl * BITCT; write_offset += BITCT) {
    sloop_max = indiv_ct - write_offset;
    if (sloop_max > BITCT2) {
      sloop_max = BITCT2;
    }
    new_geno = 0;
    new_mask = 0;
    new_missing = 0;
    // Nothing time-critical here, but may as well do it branchlessly.
    // Desired encodings:
    // new_geno: nonset homozygote -> 00
    //           het/missing       -> 01
    //           set homozygote    -> 10
    // Given PLINK encoding xx, this is (xx - (xx >> 1)).
    //
    // new_missing: autosomes: missing   -> 1
    //                         otherwise -> 0
    //              xx & ((xx >> 1) ^ 1) is one way to do this.
    //
    //              non-X haploid: het/missing -> 1
    //                             otherwise   -> 0
    //              (xx ^ (xx >> 1)) & 01 works.
    //
    //              X: [het + male]/missing -> 1
    //                 otherwise            -> 0
    //              screw it, just branch
    //
    // new_mask: missing   -> 00
    //           otherwise -> 11
    // (new_missing ^ 1) * 3 works.
    if (!is_haploid) {
      for (write_idx = 0; write_idx < sloop_max; write_idx++) {
	unfiltered_idx = next_non_set_unsafe(indiv_exclude, unfiltered_idx);
	ulii = (loadbuf[unfiltered_idx / 4] >> ((unfiltered_idx % 4) * 2)) & 3;

	uljj = ulii >> 1;
	new_geno |= (ulii - uljj) << (write_idx * 2);
	uljj = ulii & (uljj ^ 1);
	new_missing |= uljj << write_idx;
	new_mask |= ((uljj ^ 1) * 3) << (write_idx * 2);
	unfiltered_idx++;
      }
    } else if (!is_x) {
      for (write_idx = 0; write_idx < sloop_max; write_idx++) {
	unfiltered_idx = next_non_set_unsafe(indiv_exclude, unfiltered_idx);
	ulii = (loadbuf[unfiltered_idx / 4] >> ((unfiltered_idx % 4) * 2)) & 3;

	uljj = ulii >> 1;
	new_geno |= (ulii - uljj) << (write_idx * 2);
	uljj = (ulii ^ uljj) & 1;
	new_missing |= uljj << write_idx;
	new_mask |= ((uljj ^ 1) * 3) << (write_idx * 2);
	unfiltered_idx++;
      }
    } else {
      for (write_idx = 0; write_idx < sloop_max; write_idx++) {
	unfiltered_idx = next_non_set_unsafe(indiv_exclude, unfiltered_idx);
	ulii = (loadbuf[unfiltered_idx / 4] >> ((unfiltered_idx % 4) * 2)) & 3;

	uljj = ulii >> 1;
	new_geno |= (ulii - uljj) << (write_idx * 2);
	uljj = ((ulii == 1) || ((ulii == 2) && is_set(sex_male, unfiltered_idx)))? 1 : 0;
	new_missing |= uljj << write_idx;
	new_mask |= ((uljj ^ 1) * 3) << (write_idx * 2);
	unfiltered_idx++;
      }
    }
    *geno_buf++ = new_geno;
    *mask_buf++ = new_mask;
    // new_mask needed for proper post-ending handling
    sq_sum += popcount2_long((new_geno ^ FIVEMASK) & FIVEMASK & new_mask);
    sum += popcount2_long(new_geno);
    new_geno = 0;
    new_mask = 0;
    if (write_offset + BITCT2 <= indiv_ct) {
      // +0 hom1, +1 het or missing, +2 hom2
      sloop_max = indiv_ct - write_offset - BITCT2;
      if (sloop_max > BITCT2) {
	sloop_max = BITCT2;
      }
      for (write_idx = 0; write_idx < sloop_max; write_idx++) {
	unfiltered_idx = next_non_set_unsafe(indiv_exclude, unfiltered_idx);
	ulii = (loadbuf[unfiltered_idx / 4] >> ((unfiltered_idx % 4) * 2)) & 3;
	uljj = ulii >> 1;
	new_geno |= (ulii - uljj) << (write_idx * 2);
	uljj = ulii & (uljj ^ 1);
	new_missing |= uljj << (write_idx + BITCT2);
	new_mask |= ((uljj ^ 1) * 3) << (write_idx * 2);
	unfiltered_idx++;
      }
    }
    *geno_buf++ = new_geno;
    *mask_buf++ = new_mask;
    sq_sum += popcount2_long((new_geno ^ FIVEMASK) & FIVEMASK & new_mask);
    sum += popcount2_long(new_geno);
    *missing_buf++ = new_missing;
    missing_ct += popcount_long(new_missing);
  }
  fill_ulong_zero(geno_buf, indiv_trail_ct);
  fill_ulong_zero(mask_buf, indiv_trail_ct);
  non_missing_recip = 1.0 / (indiv_ct - missing_ct);
  *marker_stdev_ptr = non_missing_recip * sqrt(((int64_t)sq_sum) * (indiv_ct - missing_ct) - ((int64_t)sum) * sum);
  return missing_ct;
}

uint32_t sparse_intersection_ct(uintptr_t* sparse_buf1, uintptr_t* sparse_buf2, int32_t len) {
  int32_t ii;
  uint32_t ct = 0;
  uintptr_t ulii;
  for (ii = 0; ii < len; ii++) {
    ulii = (*sparse_buf1++) & (*sparse_buf2++);
    while (ulii) {
      ulii &= ulii - 1;
      ct++;
    }
  }
  return ct;
}

inline int32_t nz_chrom(Chrom_info* chrom_info_ptr, uintptr_t marker_idx) {
  return (marker_idx >= chrom_info_ptr->chrom_end[0]) || (marker_idx < chrom_info_ptr->chrom_start[0]);
}

uint32_t ld_prune_next_valid_chrom_start(uintptr_t* marker_exclude, uintptr_t cur_idx, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_marker_ct) {
  cur_idx = next_non_set_unsafe(marker_exclude, cur_idx);
  if (!nz_chrom(chrom_info_ptr, cur_idx)) {
    cur_idx = chrom_info_ptr->chrom_end[0];
    while (is_set(marker_exclude, cur_idx)) {
      if (cur_idx == unfiltered_marker_ct) {
	return cur_idx;
      }
      cur_idx++;
    }
  }
  return cur_idx;
}

void ld_prune_start_chrom(uint32_t ld_window_kb, uint32_t* cur_chrom_ptr, uint32_t* chrom_end_ptr, uint32_t window_unfiltered_start, uint32_t* live_indices, uint32_t* start_arr, uint32_t* window_unfiltered_end_ptr, uint32_t ld_window_size, int32_t* cur_window_size_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uint32_t* is_haploid_ptr, uint32_t* is_x_ptr) {
  uint32_t cur_chrom = get_marker_chrom(chrom_info_ptr, window_unfiltered_start);
  uint32_t window_unfiltered_end = window_unfiltered_start + 1;
  uint32_t chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
  uint32_t uii = 0;
  uint32_t species = chrom_info_ptr->species;
  uint32_t window_size;
  live_indices[0] = window_unfiltered_start;
  if (ld_window_kb) {
    window_size = 0;
    while ((window_unfiltered_start + window_size < chrom_end) && (marker_pos[window_unfiltered_start + window_size] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
      window_size++;
    }
  } else {
    window_size = ld_window_size;
  }
  for (uii = 1; uii < window_size; uii++) {
    while (is_set(marker_exclude, window_unfiltered_end)) {
      window_unfiltered_end++;
      if (window_unfiltered_end == chrom_end) {
	break;
      }
    }
    if (window_unfiltered_end == chrom_end) {
      break;
    }
    start_arr[uii - 1] = window_unfiltered_end;
    live_indices[uii] = window_unfiltered_end;
    window_unfiltered_end++;
  }
  *cur_window_size_ptr = (int)uii;
  start_arr[uii - 1] = window_unfiltered_end;
  *cur_chrom_ptr = cur_chrom;
  *chrom_end_ptr = chrom_end;
  *window_unfiltered_end_ptr = window_unfiltered_end;
  *is_haploid_ptr = (species_haploid_mask[species] >> cur_chrom) & 1LLU;
  *is_x_ptr = (species_x_code[species] == ((int32_t)cur_chrom))? 1 : 0;
}

int32_t calc_regress_pcs(char* evecname, int32_t regress_pcs_normalize_pheno, int32_t regress_pcs_sex_specific, int32_t regress_pcs_clip, int32_t max_pcs, FILE* pedfile, int32_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, double* pheno_d, double missing_phenod, FILE** outfile_ptr, char* outname, char* outname_end) {
  FILE* evecfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  int32_t retval = 0;
  int32_t pc_ct = 0;
  uint32_t pct = 1;
  int32_t is_eigenvec = 0; // GCTA .eigenvec format instead of SMARTPCA .evec?
  char wbuf[32];
  FILE* outfile;
  int32_t pc_ct_p1; // plus 1 to account for intercept
  double* pc_matrix;
  double* pc_orig_prod_sums; // pc_ct_p1 * pc_ct_p1, upper triangle filled
  double* pc_prod_sums; // (X'X)^{-1}
  double* x_prime_y; // X'Y
  double* beta_vec; // (X'X)^{-1}X'Y
  double* residual_vec;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  unsigned char* loadbuf;
  uint32_t* missing_cts;
  char* bufptr;
  char* id_buf;
  uint32_t uii;
  int32_t ii;
  int32_t jj;
  int32_t cur_missing;
  char* person_id_ptr;
  MATRIX_INVERT_BUF1_TYPE* inv_1d_buf;
  double* dbl_2d_buf;
  double dxx;
  if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_c_checked(&id_buf, max_person_id_len)) {
    return RET_NOMEM;
  }
  
  // try unaltered filename.  If that fails and the unaltered filename didn't
  // have an .evec or .eigenvec extension, then also try appending .evec and
  // appending .eigenvec.
  evecfile = fopen(evecname, "r");
  if (!evecfile) {
    ii = strlen(evecname);
    if (((ii >= 5) && (!memcmp(".evec", &(evecname[ii - 5]), 5))) || ((ii >= 9) && (!memcmp(".eigenvec", &(evecname[ii - 9]), 9)))) {
      return RET_OPEN_FAIL;
    }
    strcpy(&(evecname[ii]), ".evec");
    evecfile = fopen(evecname, "r");
    if (!evecfile) {
      strcpy(&(evecname[ii]), ".eigenvec");
      if (fopen_checked(&evecfile, evecname, "r")) {
	return RET_OPEN_FAIL;
      }
    }
  }

  tbuf[MAXLINELEN - 7] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  if (!fgets(tbuf, MAXLINELEN - 6, evecfile)) {
    if (feof(evecfile)) {
      goto calc_regress_pcs_ret_INVALID_FORMAT;
    } else {
      goto calc_regress_pcs_ret_READ_FAIL;
    }
  }
  if (!tbuf[MAXLINELEN - 7]) {
    logprint("Error: Excessively long line in .evec/.eigenvec file.\n");
    goto calc_regress_pcs_ret_INVALID_FORMAT2;
  }
  bufptr = skip_initial_spaces(tbuf);
  if (no_more_items_kns(bufptr)) {
    goto calc_regress_pcs_ret_INVALID_FORMAT;
  }
  if (memcmp(bufptr, "#eigvals:", 9)) {
    is_eigenvec = 1;
    bufptr = next_item(bufptr);
  }
  bufptr = next_item(bufptr);
  while ((!no_more_items_kns(bufptr)) && ((*bufptr == '-') || ((*bufptr >= '0') && (*bufptr <= '9')))) {
    pc_ct++;
    bufptr = next_item(bufptr);
  }
  if (!pc_ct) {
    goto calc_regress_pcs_ret_INVALID_FORMAT;
  }
  if (pc_ct > max_pcs) {
    sprintf(logbuf, "%svec format detected.  Regressing on %d PC%s (out of %d).\n", is_eigenvec? "GCTA .eigen" : "SMARTPCA .e", max_pcs, (max_pcs == 1)? "" : "s", pc_ct);
    pc_ct = max_pcs;
  } else {
    sprintf(logbuf, "%svec format detected.  Regressing on %d principal component%s.\n", is_eigenvec? "GCTA .eigen" : "SMARTPCA .e", pc_ct, (pc_ct == 1)? "" : "s");
  }
  logprintb();
  pc_ct_p1 = pc_ct + 1;
  if (wkspace_alloc_d_checked(&pc_matrix, pc_ct_p1 * indiv_ct * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&pc_orig_prod_sums, pc_ct_p1 * pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&pc_prod_sums, pc_ct_p1 * pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&x_prime_y, pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&beta_vec, pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&residual_vec, pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  inv_1d_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(pc_ct_p1 * sizeof(MATRIX_INVERT_BUF1_TYPE));
  if (!inv_1d_buf) {
    goto calc_regress_pcs_ret_NOMEM;
  }
  if (wkspace_alloc_d_checked(&dbl_2d_buf, pc_ct_p1 * pc_ct_p1 * sizeof(double))) {
    goto calc_regress_pcs_ret_NOMEM;
  }

  if (is_eigenvec) {
    indiv_idx = 0;
    while (1) {
      // todo: validate, and perhaps permute, family/indiv IDs
      bufptr = next_item_mult(skip_initial_spaces(tbuf), 2);
      for (ii = 0; ii < pc_ct; ii++) {
	if (no_more_items_kns(bufptr)) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	if (sscanf(bufptr, "%lg", &(pc_matrix[ii * indiv_ct + indiv_idx])) != 1) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	bufptr = next_item(bufptr);
      }
      pc_matrix[pc_ct * indiv_ct + indiv_idx] = 1.0; // intercept
      if (++indiv_idx >= indiv_ct) {
	break;
      }
      if (!fgets(tbuf, MAXLINELEN, evecfile)) {
	if (feof(evecfile)) {
	  sprintf(logbuf, "Error: Fewer %s in .eigenvec file than expected.\n", species_plural);
	  goto calc_regress_pcs_ret_INVALID_FORMAT_3;
	} else {
	  goto calc_regress_pcs_ret_READ_FAIL;
	}
      }
    }
  } else {
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      if (!fgets(tbuf, MAXLINELEN, evecfile)) {
	if (feof(evecfile)) {
	  sprintf(logbuf, "Error: Fewer %s in .evec file than expected.\n", species_plural);
	  goto calc_regress_pcs_ret_INVALID_FORMAT_3;
	} else {
	  goto calc_regress_pcs_ret_READ_FAIL;
	}
      }
      bufptr = next_item(skip_initial_spaces(tbuf));
      for (ii = 0; ii < pc_ct; ii++) {
	if (no_more_items_kns(bufptr)) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	if (sscanf(bufptr, "%lg", &(pc_matrix[ii * indiv_ct + indiv_idx])) != 1) {
	  goto calc_regress_pcs_ret_INVALID_FORMAT;
	}
	bufptr = next_item(bufptr);
      }
      pc_matrix[pc_ct * indiv_ct + indiv_idx] = 1.0;
    }
  }
  if (fgets(tbuf, MAXLINELEN, evecfile)) {
    if (!no_more_items_kns(skip_initial_spaces(tbuf))) {
      sprintf(logbuf, "Error: More %s in .e%svec file than expected.\n", species_plural, is_eigenvec? "igen" : "");
      goto calc_regress_pcs_ret_INVALID_FORMAT_3;
    }
  }
  fclose_null(&evecfile);

  // precalculate (X'X)
  fill_double_zero(pc_orig_prod_sums, pc_ct_p1 * pc_ct_p1);
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    for (ii = 0; ii < pc_ct_p1; ii++) {
      for (jj = ii; jj < pc_ct_p1; jj++) {
        pc_orig_prod_sums[ii * pc_ct_p1 + jj] += pc_matrix[ii * indiv_ct + indiv_idx] * pc_matrix[jj * indiv_ct + indiv_idx];
      }
    }
  }

  fill_uint_zero(missing_cts, indiv_ct);
  marker_uidx = 0;
  // .gen instead of .bgen because latter actually has lower precision(!) (15
  // bits instead of the ~20 you get from printf("%g", dxx)), and there's no
  // need for repeated random access.
  strcpy(outname_end, ".gen");
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  if (fseeko(pedfile, bed_offset, SEEK_SET)) {
    return RET_READ_FAIL;
  }
  outfile = *outfile_ptr;
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    if (is_set(marker_exclude, marker_uidx)) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      if (fseeko(pedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	return RET_READ_FAIL;
      }
    }
    if (fread(loadbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
      return RET_READ_FAIL;
    }
    bufptr = uint32_writex(tbuf, get_marker_chrom(chrom_info_ptr, marker_uidx), ' ');
    fwrite(tbuf, 1, bufptr - tbuf, outfile);
    fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
    tbuf[0] = ' ';
    bufptr = uint32_writex(&(tbuf[1]), marker_pos[marker_uidx], ' ');
    if (max_marker_allele_len == 1) {
      bufptr[0] = marker_alleles[2 * marker_uidx];
      bufptr[1] = ' ';
      bufptr[2] = marker_alleles[2 * marker_uidx + 1];
      if (fwrite_checked(tbuf, 3 + (uintptr_t)(bufptr - tbuf), outfile)) {
        return RET_WRITE_FAIL;
      }
    } else {
      fwrite(tbuf, 1, bufptr - tbuf, outfile);
      fputs(&(marker_alleles[2 * marker_uidx * max_marker_allele_len]), outfile);
      putc(' ', outfile);
      if (fputs(&(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]), outfile) == EOF) {
        return RET_WRITE_FAIL;
      }
    }
    memcpy(pc_prod_sums, pc_orig_prod_sums, pc_ct_p1 * pc_ct_p1 * sizeof(double));
    fill_double_zero(x_prime_y, pc_ct_p1);
    indiv_uidx = 0;
    cur_missing = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
      uii = ((loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3);
      if (uii == 1) {
	cur_missing++;
	missing_cts[indiv_idx] += 1;
	for (ii = 0; ii < pc_ct_p1; ii++) {
	  for (jj = ii; jj < pc_ct_p1; jj++) {
	    pc_prod_sums[ii * pc_ct_p1 + jj] -= pc_matrix[ii * indiv_ct + indiv_idx] * pc_matrix[jj * indiv_ct + indiv_idx];
	  }
	}
      } else {
	uii = uii - (uii >> 1);
        for (ii = 0; ii < pc_ct_p1; ii++) {
	  x_prime_y[ii] += pc_matrix[ii * indiv_ct + indiv_idx] * uii;
	}
      }
      indiv_uidx++;
    }
    // last-minute update of lower triangle
    for (ii = 1; ii < pc_ct_p1; ii++) {
      for (jj = 0; jj < ii; jj++) {
	pc_prod_sums[ii * pc_ct_p1 + jj] = pc_prod_sums[jj * pc_ct_p1 + ii];
      }
    }
    invert_matrix(pc_ct_p1, pc_prod_sums, inv_1d_buf, dbl_2d_buf);
    for (ii = 0; ii < pc_ct_p1; ii++) {
      dxx = 0.0;
      for (jj = 0; jj < pc_ct_p1; jj++) {
        dxx += pc_prod_sums[ii * pc_ct_p1 + jj] * x_prime_y[jj];
      }
      beta_vec[ii] = dxx;
    }
    indiv_uidx = 0;
    wbuf[0] = ' ';
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
      uii = ((loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3);
      if (uii == 1) {
	fputs(" 0 0 0", outfile);
      } else {
	uii = uii - (uii >> 1);
	dxx = 0.0;
        for (ii = 0; ii < pc_ct_p1; ii++) {
          dxx += pc_matrix[ii * indiv_ct + indiv_idx] * beta_vec[ii];
	}
	dxx = (double)uii - dxx;
	// now dxx is the residual, normally but not always in [0, 2]
	if (dxx < 1.0) {
	  if (dxx < 0.0) {
	    if (regress_pcs_clip) {
	      fputs(" 1 0 0", outfile);
	    } else {
	      bufptr = double_g_write(memcpyl3a(double_g_write(&(wbuf[1]), 1.0 - dxx * 0.5), " 0 "), dxx * 0.5);
	      fwrite(wbuf, 1, bufptr - wbuf, outfile);
	    }
	  } else {
	    bufptr = memcpya(double_g_write(double_g_writex(&(wbuf[1]), 1.0 - dxx, ' '), dxx), " 0", 2);
	    fwrite(wbuf, 1, bufptr - wbuf, outfile);
	  }
	} else {
	  if (dxx > 2.0) {
	    if (regress_pcs_clip) {
	      fputs(" 0 0 1", outfile);
	    } else {
	      bufptr = double_g_write(memcpyl3a(double_g_write(&(wbuf[1]), 1.0 - dxx * 0.5), " 0 "), dxx * 0.5);
	      fwrite(wbuf, 1, bufptr - wbuf, outfile);
	    }
	  } else {
            bufptr = double_g_write(double_g_writex(memcpya(&(wbuf[1]), "0 ", 2), 2.0 - dxx, ' '), dxx - 1.0);
	    fwrite(wbuf, 1, bufptr - wbuf, outfile);
	  }
	}
      }
      indiv_uidx++;
    }
    if (putc('\n', outfile) == EOF) {
      return RET_WRITE_FAIL;
    }
    if (marker_idx * 100LLU >= ((uint64_t)pct * marker_ct)) {
      pct = ((uint64_t)marker_idx * 100) / marker_ct;
      printf("\r%d%%", pct++);
      fflush(stdout);
    }
    marker_uidx++;
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  strcpy(outname_end, ".sample");
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  outfile = *outfile_ptr;
  if (fputs_checked("ID_1 ID_2 missing sex phenotype\n0 0 0 D P\n", outfile)) {
    return RET_WRITE_FAIL;
  }
  // regress phenotype
  fill_double_zero(x_prime_y, pc_ct_p1);
  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    dxx = pheno_d[indiv_uidx];
    for (ii = 0; ii < pc_ct_p1; ii++) {
      x_prime_y[ii] += pc_matrix[ii * indiv_ct + indiv_idx] * dxx;
    }
    indiv_uidx++;
  }
  for (ii = 1; ii < pc_ct_p1; ii++) {
    for (jj = 0; jj < ii; jj++) {
      pc_orig_prod_sums[ii * pc_ct_p1 + jj] = pc_orig_prod_sums[jj * pc_ct_p1 + ii];
    }
  }
  invert_matrix(pc_ct_p1, pc_orig_prod_sums, inv_1d_buf, dbl_2d_buf);
  for (ii = 0; ii < pc_ct_p1; ii++) {
    dxx = 0.0;
    for (jj = 0; jj < pc_ct_p1; jj++) {
      dxx += pc_orig_prod_sums[ii * pc_ct_p1 + jj] * x_prime_y[jj];
    }
    beta_vec[ii] = dxx;
  }

  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    dxx = 0.0;
    for (ii = 0; ii < pc_ct_p1; ii++) {
      dxx += pc_matrix[ii * indiv_ct + indiv_idx] * beta_vec[ii];
    }
    residual_vec[indiv_idx] = pheno_d[indiv_uidx] - dxx;
    indiv_uidx++;
  }

  if (regress_pcs_normalize_pheno) {
    if (regress_pcs_sex_specific) {
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_nm, sex_male, 1);
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_nm, sex_male, 2);
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_nm, sex_male, 4);
    } else {
      normalize_phenos(residual_vec, indiv_ct, unfiltered_indiv_ct, indiv_exclude, sex_nm, sex_male, 7);
    }
  }

  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    person_id_ptr = &(person_ids[indiv_uidx * max_person_id_len]);
    uii = strlen_se(person_id_ptr);
    // todo: adjust pheno_d, double-check missing gender behavior
    fwrite(person_id_ptr, 1, uii, outfile);
    putc(' ', outfile);
    fputs(&(person_id_ptr[uii + 1]), outfile);
    tbuf[0] = ' ';
    bufptr = double_g_writex(&(tbuf[1]), (double)missing_cts[indiv_uidx] / (double)marker_ct, ' ');
    *bufptr = sexchar(sex_nm, sex_male, indiv_uidx);
    bufptr[1] = ' ';
    bufptr = double_g_writex(&(bufptr[2]), residual_vec[indiv_idx], '\n');
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      return RET_WRITE_FAIL;
    }
    indiv_uidx++;
  }
  if (fclose_null(outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  *outname_end = '\0';
  putchar('\r');
  sprintf(logbuf, "Principal component regression residuals and %sphenotype Z-scores %s%s.gen and %s.sample.\n", regress_pcs_sex_specific? "sex-specific " : "", regress_pcs_sex_specific? "\nwritten to " : "written to\n", outname, outname);
  logprintb();
  wkspace_reset(wkspace_mark);
  while (0) {
  calc_regress_pcs_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_regress_pcs_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_regress_pcs_ret_INVALID_FORMAT_3:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  calc_regress_pcs_ret_INVALID_FORMAT:
    logprint("Error: Improperly formatted .evec file.\n");
  calc_regress_pcs_ret_INVALID_FORMAT2:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(evecfile);
  return retval;
}

static uintptr_t g_cg_indiv1idx;
static uintptr_t g_cg_indiv2idx;
static uintptr_t g_cg_tstc;
static uint32_t g_cg_marker_ct;
static uintptr_t g_cg_max_person_id_len;
static uint32_t g_cg_max_person_fid_len;
static uint32_t g_cg_max_person_iid_len;
static uintptr_t g_cg_max_paternal_id_len;
static uintptr_t g_cg_max_maternal_id_len;
static uintptr_t* g_cg_pheno_nm;
static uintptr_t* g_cg_pheno_c;
static uintptr_t* g_cg_founder_info;
static char* g_cg_person_ids;
static char* g_cg_paternal_ids;
static char* g_cg_maternal_ids;
static int32_t g_cg_genome_output_full;
static int32_t g_cg_genome_ibd_unbounded;
static Pedigree_rel_info g_cg_pri;
static uint32_t g_cg_family_id_fixed;
static int32_t g_cg_is_founder_fixed;
static uint32_t g_cg_rel_space_id_fixed;
static int64_t g_cg_llfct;
static char* g_cg_sptr_start;
static char* g_cg_fam1;
static char* g_cg_fam2;
static char* g_cg_indiv1;
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

uint32_t calc_genome_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  char* cptr;
  char* indiv2;
  char* pat2 = NULL;
  char* mat2 = NULL;
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
  if (!g_cg_indiv2idx) {
    // first line, if not 2nd or later part of parallel write
    if (!g_cg_indiv1idx) {
      sptr_cur += sprintf(sptr_cur, tbuf, " FID1", " IID1", " FID2", " IID2");
    }
    g_cg_indiv1idx = g_thread_start[0];
    g_cg_indiv2idx = g_cg_indiv1idx + 1;
    tbuf[0] = ' ';
  }
  do {
    if (g_cg_indiv2idx == g_cg_indiv1idx + 1) {
      cptr = &(g_cg_person_ids[g_cg_indiv1idx * g_cg_max_person_id_len]);
      uii = strlen_se(cptr);
      memcpyx(g_cg_fam1, cptr, uii, '\0');
      g_cg_indiv1 = next_item(cptr);
      if (g_cg_paternal_ids) {
	g_cg_pat1 = &(g_cg_paternal_ids[g_cg_indiv1idx * g_cg_max_paternal_id_len]);
	g_cg_mat1 = &(g_cg_maternal_ids[g_cg_indiv1idx * g_cg_max_maternal_id_len]);
	g_cg_is_founder_fixed = is_founder(g_cg_founder_info, g_cg_indiv1idx);
	g_cg_rel_space_id_fixed = g_cg_pri.family_rel_nf_idxs[g_cg_indiv1idx];
	g_cg_family_id_fixed = g_cg_pri.family_idxs[g_cg_indiv1idx];
	founder_ct = g_cg_pri.family_founder_cts[g_cg_family_id_fixed];
	g_cg_llfct = (int64_t)founder_ct * (founder_ct - 1);
      }
      g_cg_sptr_start = fw_strcpyn(g_cg_max_person_fid_len - 1, uii, g_cg_fam1, &(tbuf[1]));
      *g_cg_sptr_start++ = ' ';
      g_cg_sptr_start = fw_strcpy(g_cg_max_person_iid_len - 1, g_cg_indiv1, g_cg_sptr_start);
      *g_cg_sptr_start++ = ' ';
    }
    while (g_cg_indiv2idx < g_indiv_ct) {
      sptr_cur = memcpya(sptr_cur, tbuf, g_cg_sptr_start - tbuf);
      cptr = &(g_cg_person_ids[g_cg_indiv2idx * g_cg_max_person_id_len]);
      uii = strlen_se(cptr);
      memcpyx(g_cg_fam2, cptr, uii, '\0');
      indiv2 = skip_initial_spaces(&(cptr[uii + 1]));
      if (g_cg_paternal_ids) {
	pat2 = &(g_cg_paternal_ids[g_cg_indiv2idx * g_cg_max_paternal_id_len]);
	mat2 = &(g_cg_maternal_ids[g_cg_indiv2idx * g_cg_max_maternal_id_len]);
      }
      sptr_cur = fw_strcpyn(g_cg_max_person_fid_len - 1, uii, g_cg_fam2, sptr_cur);
      *sptr_cur++ = ' ';
      sptr_cur = fw_strcpy(g_cg_max_person_iid_len - 1, indiv2, sptr_cur);
      *sptr_cur++ = ' ';
      if (!strcmp(g_cg_fam1, g_cg_fam2)) {
	while (1) {
	  if (g_cg_paternal_ids) {
	    if (!(g_cg_is_founder_fixed || is_set(g_cg_founder_info, g_cg_indiv2idx))) {
	      if ((!strcmp(g_cg_pat1, pat2)) && (!strcmp(g_cg_mat1, mat2))) {
		sptr_cur = memcpyl3a(sptr_cur, "FS ");
		break;
	      } else if ((!strcmp(g_cg_pat1, pat2)) || (!strcmp(g_cg_mat1, mat2))) {
		sptr_cur = memcpyl3a(sptr_cur, "HS ");
		break;
	      }
	    }
	    if ((!strcmp(g_cg_pat1, indiv2)) || (!strcmp(g_cg_mat1, indiv2)) || (!strcmp(pat2, g_cg_indiv1)) || (!strcmp(mat2, g_cg_indiv1))) {
	      sptr_cur = memcpyl3a(sptr_cur, "PO ");
	      break;
	    }
	  }
	  sptr_cur = memcpyl3a(sptr_cur, "OT ");
	  break;
	}
	// insert relationship
	if (!g_cg_paternal_ids) {
	  sptr_cur = memcpya(sptr_cur, "    0", 5);
	} else {
	  oo = is_set(g_cg_founder_info, g_cg_indiv2idx);
	  if (g_cg_is_founder_fixed && oo) {
	    dxx = 0.0;
	  } else {
	    nn = g_cg_pri.family_rel_nf_idxs[g_cg_indiv2idx];
	    if (g_cg_is_founder_fixed || ((nn > (int32_t)g_cg_rel_space_id_fixed) && (!oo))) {
	      dxx = g_cg_pri.rel_space[g_cg_rel_space_id_fixed + ((int64_t)nn * (nn - 1) - g_cg_llfct) / 2];
	    } else {
	      dxx = g_cg_pri.rel_space[nn + ((int64_t)g_cg_rel_space_id_fixed * (g_cg_rel_space_id_fixed - 1) - g_cg_llfct) / 2];
	    }
	  }
	  sptr_cur = width_force(5, sptr_cur, double_g_write(sptr_cur, dxx));
	}
      } else {
	sptr_cur = memcpya(sptr_cur, "UN    NA", 8);
      }
      nn = g_cg_marker_ct - g_indiv_missing_unwt[g_cg_indiv1idx] - g_indiv_missing_unwt[g_cg_indiv2idx] + g_missing_dbl_excluded[g_cg_mdecell];
      oo = nn - g_genome_main[g_cg_gmcell] - g_genome_main[g_cg_gmcell + 1];
      dxx = (double)g_genome_main[g_cg_gmcell + 1] / (g_cg_e00 * nn);
      dyy = ((double)g_genome_main[g_cg_gmcell] - dxx * g_cg_e01 * nn) / (g_cg_e11 * nn);
      dxx1 = ((double)oo - nn * (dxx * g_cg_e02 + dyy * g_cg_e12)) / ((double)nn);
      if (!g_cg_genome_ibd_unbounded) {
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
	  dxx1 = 0;
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
      *sptr_cur++ = ' ';
      sptr_cur = double_f_writew74(double_f_writew74x(double_f_writew74x(double_f_writew74x(sptr_cur, dxx, ' '), dyy, ' '), dxx1, ' '), dyy * 0.5 + dxx1);

      if (g_cg_pheno_c) {
	uii = is_set(g_cg_pheno_nm, g_cg_indiv1idx);
	ujj = is_set(g_cg_pheno_nm, g_cg_indiv2idx);
	ukk = is_set(g_cg_pheno_c, g_cg_indiv1idx);
	umm = is_set(g_cg_pheno_c, g_cg_indiv2idx);
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
      dxx = (double)g_genome_main[g_cg_gmcell + 4];
      dyy = (double)g_genome_main[g_cg_gmcell + 3];
      dxx1 = 1.0 / ((double)(g_genome_main[g_cg_gmcell + 4] + g_genome_main[g_cg_gmcell + 3]));
      dxx2 = normdist((dxx * dxx1 - 0.666666) / (sqrt(0.2222222 * dxx1)));
      sptr_cur = double_f_writew74x(double_f_writew6x(sptr_cur, 1.0 - (g_genome_main[g_cg_gmcell] + 2 * g_genome_main[g_cg_gmcell + 1]) / ((double)(2 * nn)), ' '), dxx2, ' ');
      if (g_genome_main[g_cg_gmcell + 3]) {
	sptr_cur = double_f_writew74(sptr_cur, dxx / dyy);
      } else {
	sptr_cur = memcpya(sptr_cur, "     NA", 7);
      }
      if (g_cg_genome_output_full) {
	*sptr_cur = ' ';
	sptr_cur = double_f_writew74(double_f_writew74x(uint32_writew7x(uint32_writew7x(uint32_writew7x(&(sptr_cur[1]), g_genome_main[uii + 1], ' '), g_genome_main[g_cg_gmcell], ' '), oo, ' '), dyy, ' '), dxx);
      }
      *sptr_cur++ = '\n';
      g_cg_gmcell += 5;
      g_cg_mdecell++;
      g_cg_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_genome_emitn_ret;
      }
    }
    g_cg_cur_line += g_indiv_ct - g_cg_indiv1idx - 1;
    g_cg_indiv1idx++;
    g_cg_indiv2idx = g_cg_indiv1idx + 1;
  } while (g_cg_indiv1idx < g_cg_tstc);
 calc_genome_emitn_ret:
  if (g_cg_cur_line * 100 >= g_cg_tot_lines * g_pct) {
    g_pct = (g_cg_cur_line * 100) / g_cg_tot_lines;
    printf("\rWriting... %u%%", g_pct++);
    fflush(stdout);
  }
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t calc_genome(pthread_t* threads, FILE* bedfile, int32_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, double* set_allele_freqs, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info, int32_t parallel_idx, int32_t parallel_tot, char* outname, char* outname_end, int32_t nonfounders, uint64_t calculation_type, int32_t genome_output_gz, int32_t genome_output_full, int32_t genome_ibd_unbounded, int32_t ppc_gap, uintptr_t* pheno_nm, uintptr_t* pheno_c, Pedigree_rel_info pri) {
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  int32_t retval = 0;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  char wbuf[16];
  unsigned char* loadbuf; // from file
  unsigned char* gptr;
  char* cptr;
  char* cptr2;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  int32_t nn;
  int32_t oo;
  int32_t pp;
  int32_t qq;
  uintptr_t* glptr;
  uintptr_t* glptr2;
  uintptr_t* glptr3;
  uint32_t* giptr;
  uint32_t* giptr2;
  uint32_t* giptr3;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  int32_t missing_ct_buf[BITCT];
  int32_t missing_ct_all;
  double set_allele_freq_buf[GENOME_MULTIPLEX];
  int32_t ibd_prect = 0;
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
  int64_t cur_line = 0;
  int64_t tot_cells;
  double num_alleles;
  double num_allelesf2;
  double num_allelesf3;
  uint32_t mp_lead_unfiltered_idx = 0;
  uint32_t mp_lead_idx = 0;
  uint32_t chrom_fo_idx = 0;
  uint32_t chrom_end;
  uint32_t is_x;
  uint32_t is_haploid;
  uintptr_t marker_uidx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uint32_t pct;

  // imitate PLINK behavior (see Plink::prettyPrintLengths() in helper.cpp), to
  // avoid randomly breaking existing scripts
  g_cg_max_person_fid_len = 4;
  g_cg_max_person_iid_len = 4;

  triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_tot - parallel_idx - 1, parallel_tot, 1, 1);
  // invert order, for --genome --parallel to naturally work
  for (uii = 0; uii <= g_thread_ct / 2; uii++) {
    jj = g_thread_start[uii];
    g_thread_start[uii] = g_indiv_ct - g_thread_start[g_thread_ct - uii];
    g_thread_start[g_thread_ct - uii] = g_indiv_ct - jj;
  }

  if (!parallel_idx) {
    cur_line = 1;
  }
  g_cg_tstc = g_thread_start[g_thread_ct];
  // f(tstc) - f(g_thread_start[0])
  // f(0) = 0
  // f(1) = indiv_ct - 1
  // f(2) = 2indiv_ct - 3
  // ...
  // f(n) = nindiv_ct - n(n+1)/2
  tot_cells = (int64_t)g_indiv_ct * (g_cg_tstc - g_thread_start[0]) - ((int64_t)g_cg_tstc * (g_cg_tstc + 1) - (int64_t)g_thread_start[0] * (g_thread_start[0] + 1)) / 2;
  g_cg_tot_lines = cur_line + tot_cells;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (!is_set(indiv_exclude, indiv_uidx)) {
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      uii = strlen_se(cptr);
      cptr2 = next_item(cptr);
      if (uii > g_cg_max_person_fid_len) {
	g_cg_max_person_fid_len = uii + 2;
      }
      uii = strlen_se(cptr2);
      if (uii > g_cg_max_person_iid_len) {
        g_cg_max_person_iid_len = uii + 2;
      }
    }
  }
  if (wkspace_alloc_ui_checked(&g_missing_dbl_excluded, tot_cells * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&g_indiv_missing_unwt, g_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&g_genome_main, tot_cells * 5 * sizeof(int32_t)) ||
      wkspace_alloc_uc_checked(&loadbuf, GENOME_MULTIPLEX * unfiltered_indiv_ct4) ||
      wkspace_alloc_uc_checked(&g_geno, g_indiv_ct * (GENOME_MULTIPLEX / 4)) ||
      wkspace_alloc_ul_checked(&g_masks, g_indiv_ct * (GENOME_MULTIPLEX / 4)) ||
      wkspace_alloc_ul_checked(&g_mmasks, g_indiv_ct * sizeof(intptr_t)) ||
      wkspace_alloc_c_checked(&g_cg_fam1, g_cg_max_person_fid_len + 1) ||
      wkspace_alloc_c_checked(&g_cg_fam2, g_cg_max_person_fid_len + 1)) {
    goto calc_genome_ret_NOMEM;
  }

  fill_int_zero((int32_t*)g_missing_dbl_excluded, tot_cells);
  fill_int_zero((int32_t*)g_indiv_missing_unwt, g_indiv_ct);
  fill_int_zero((int32_t*)g_genome_main, tot_cells * 5);
  if (!is_set(marker_exclude, 0)) {
    if (fseeko(bedfile, bed_offset, SEEK_SET)) {
      retval = RET_READ_FAIL;
      goto calc_genome_ret_1;
    }
  }
  marker_uidx = 0; // raw marker index
  g_ctrl_ct = 0; // after excluding missing (abuse of variable name, should fix this later)
  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
  // subtract X/haploid markers from marker_ct
  ukk = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1);
  if (ukk) {
    sprintf(logbuf, "Excluding %u marker%s on non-autosomes from IBD calculation.\n", ukk, (ukk == 1)? "" : "s");
    logprintb();
    marker_ct -= ukk;
  }
  g_cg_e00 = 0.0;
  g_cg_e01 = 0.0;
  g_cg_e02 = 0.0;
  g_cg_e11 = 0.0;
  g_cg_e12 = 0.0;
  do {
    kk = marker_ct - g_ctrl_ct;
    if (kk > GENOME_MULTIPLEX) {
      kk = GENOME_MULTIPLEX;
    }
    glptr2 = g_marker_window;
    for (ujj = 0; ujj < (uint32_t)kk; ujj++) {
      if (is_set(marker_exclude, marker_uidx)) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  retval = RET_READ_FAIL;
	  goto calc_genome_ret_1;
	}
      }
      if (marker_uidx >= chrom_end) {
	while (1) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  if (!is_haploid) {
	    break;
	  }
          marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    retval = RET_READ_FAIL;
	    goto calc_genome_ret_1;
	  }
	}
      }
      if (fread(&(loadbuf[ujj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	retval = RET_READ_FAIL;
	goto calc_genome_ret_1;
      }
      set_allele_freq_buf[ujj] = set_allele_freqs[marker_uidx];
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
	mp_lead_idx = g_ctrl_ct + ujj;
      }
      if (mp_lead_unfiltered_idx < chrom_end) {
	uii = ppc_gap + marker_pos[marker_uidx];
	if (marker_pos[mp_lead_unfiltered_idx] <= uii) {
	  do {
	    if (!is_set(marker_exclude, mp_lead_unfiltered_idx)) {
	      mp_lead_idx++;
	    }
	    mp_lead_unfiltered_idx++;
	  } while ((mp_lead_unfiltered_idx < chrom_end) && (is_set(marker_exclude, mp_lead_unfiltered_idx) || (marker_pos[mp_lead_unfiltered_idx] <= uii)));
	}
      }
      if (mp_lead_unfiltered_idx < unfiltered_marker_ct) {
	ulii = 2 * (mp_lead_idx - g_ctrl_ct);
	if (ulii < BITCT + (2 * (ujj & (~(BITCT2 - 1))))) {
	  ulii = ~ZEROLU << (ulii & (BITCT - 1));
	}
      } else {
	ulii = 2 * (unfiltered_marker_ct + GENOME_MULTIPLEX);
      }

      *glptr2++ = ulii;
      *glptr2++ = ulii;
      marker_uidx++;
    }
    if (kk < GENOME_MULTIPLEX) {
      memset(&(loadbuf[kk * unfiltered_indiv_ct4]), 0, (GENOME_MULTIPLEX - kk) * unfiltered_indiv_ct4);
      fill_long_zero((intptr_t*)g_geno, g_indiv_ct * (GENOME_MULTIPLEX / BITCT2));
      fill_ulong_zero(g_masks, g_indiv_ct * (GENOME_MULTIPLEX / BITCT2));
      for (mm = kk * 2; mm < GENOME_MULTIPLEX2; mm++) {
	*glptr2++ = GENOME_MULTIPLEX2;
      }
    }
    g_case_ct = g_ctrl_ct + kk;
    for (jj = 0; jj < kk; jj += BITCT) {
      glptr = &(((uintptr_t*)g_geno)[jj / BITCT2]);
      glptr2 = &(g_masks[jj / BITCT2]);
      glptr3 = g_mmasks;
      giptr = g_indiv_missing_unwt;
      indiv_uidx = 0;
      fill_int_zero(missing_ct_buf, BITCT);
      missing_ct_all = 0;
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	oo = (indiv_uidx % 4) * 2;
	ulii = 0;
	ulkk = 0;
        gptr = &(loadbuf[indiv_uidx / 4 + jj * unfiltered_indiv_ct4]);
	qq = (nonfounders || is_founder(founder_info, indiv_uidx));
	if (qq) {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[pp] += 1;
	      ulkk |= ONELU << pp;
	      *giptr += 1;
	    }
	  }
	} else {
	  missing_ct_all++;
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      ulkk |= ONELU << pp;
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
	gptr = &(loadbuf[indiv_uidx / 4 + (jj + BITCT2) * unfiltered_indiv_ct4]);
	if (qq) {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[pp + BITCT2] += 1;
	      ulkk |= ONELU << pp;
	      *giptr += 1;
	    }
	  }
	} else {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      ulkk |= ONELU << pp;
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
	indiv_uidx++;
      }
      nn = kk - jj;
      if (nn > BITCT) {
	nn = BITCT;
      }
      for (mm = 0; mm < nn; mm++) {
	dpp = set_allele_freq_buf[jj + mm];
	dqq = 1.0 - dpp;
	num_alleles = (double)(2 * (g_indiv_ct - missing_ct_buf[mm] - missing_ct_all));
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
	  g_cg_e00 += 2 * dpp_sq * dqq_sq * dxx1 * dyy1 * num_allelesf3;
	  g_cg_e01 += 4 * dpp * dqq * num_allelesf3 * (dpp_sq * dxx2 + dqq_sq * dyy2);
	  g_cg_e02 += num_allelesf3 * (dqq_sq * dqq_sq * dyy2 * (dyy - 3) * dyy_recip + dpp_sq * dpp_sq * dxx2 * (dxx - 3) * dxx_recip + 4 * dpp_sq * dqq_sq * dxx1 * dyy1);
	  g_cg_e11 += 2 * dpp * dqq * num_allelesf2 * (dpp * dxx1 + dqq * dyy1);
	  g_cg_e12 += num_allelesf2 * (dpp_sq * dpp * dxx2 + dqq_sq * dqq * dyy2 + dpp_sq * dqq * dxx1 + dpp * dqq_sq * dyy1);
	  ibd_prect++;
	}
      }
      if (spawn_threads(threads, &calc_genomem_thread, g_thread_ct)) {
	goto calc_genome_ret_THREAD_CREATE_FAIL;
      }
      incr_dists_rm_inv(g_missing_dbl_excluded, 0);
      join_threads(threads, g_thread_ct);
    }

    if (spawn_threads(threads, &calc_genome_thread, g_thread_ct)) {
      goto calc_genome_ret_THREAD_CREATE_FAIL;
    }
    incr_genome(g_genome_main, (uintptr_t*)g_geno, 0);
    join_threads(threads, g_thread_ct);
    g_ctrl_ct = g_case_ct;
    printf("\r%d markers complete.", g_ctrl_ct);
    fflush(stdout);
  } while (g_ctrl_ct < marker_ct);
  fputs("\rIBD calculations complete.  \n", stdout);
  logstr("IBD calculations complete.\n");
  dxx = 1.0 / (double)ibd_prect;
  g_cg_e00 *= dxx;
  g_cg_e01 *= dxx;
  g_cg_e02 *= dxx;
  g_cg_e11 *= dxx;
  g_cg_e12 *= dxx;

  if (calculation_type & CALC_PLINK_IBS_MATRIX) {
    strcpy(outname_end, ".mibs");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    giptr = g_genome_main;
    giptr2 = g_missing_dbl_excluded;
    pct = 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      giptr3 = g_indiv_missing_unwt;
      uii = marker_ct - giptr3[indiv_idx];
      uljj = (int)indiv_idx - 1; // not referenced when indiv_idx == 0
      for (ulii = 0; ulii < indiv_idx; ulii++) {
        cptr = double_g_writex(wbuf, 1.0 - ((double)(g_genome_main[uljj * 5] + 2 * g_genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + g_missing_dbl_excluded[uljj]))), ' ');
        fwrite(wbuf, 1, cptr - wbuf, outfile);
	uljj += g_indiv_ct - ulii - 2;
      }
      putc('1', outfile);
      putc(' ', outfile);
      giptr3++;
      for (ujj = indiv_idx + 1; ujj < g_indiv_ct; ujj++) {
	cptr = double_g_writex(wbuf, 1.0 - ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++)))), ' ');
	fwrite(wbuf, 1, cptr - wbuf, outfile);
	giptr = &(giptr[5]);
      }
      if (putc('\n', outfile) == EOF) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (indiv_idx * 100 >= (pct * g_indiv_ct)) {
        pct = (indiv_idx * 100) / g_indiv_ct;
	printf("\rWriting... %d%%", pct++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    putchar('\r');
    sprintf(logbuf, "IBS matrix written to %s.\n", outname);
    logprintb();
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }

  if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
    strcpy(outname_end, ".mdist");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPEN_FAIL;
    }
    giptr = g_genome_main;
    giptr2 = g_missing_dbl_excluded;
    kk = 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      giptr3 = g_indiv_missing_unwt;
      uii = marker_ct - giptr3[indiv_idx];
      uljj = indiv_idx - 1;
      for (ulii = 0; ulii < indiv_idx; ulii++) {
	cptr = double_g_writex(wbuf, ((double)(g_genome_main[uljj * 5] + 2 * g_genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + g_missing_dbl_excluded[uljj]))), ' ');
	fwrite(wbuf, 1, cptr - wbuf, outfile);
	uljj += g_indiv_ct - ulii - 2;
      }
      putc('0', outfile);
      putc(' ', outfile);
      giptr3++;
      for (ujj = indiv_idx + 1; ujj < g_indiv_ct; ujj++) {
	cptr = double_g_writex(wbuf, ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++)))), ' ');
        fwrite(wbuf, 1, cptr - wbuf, outfile);
	giptr = &(giptr[5]);
      }
      if (putc('\n', outfile) == EOF) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (indiv_idx * 100 >= (kk * g_indiv_ct)) {
	kk = (indiv_idx * 100) / g_indiv_ct;
	printf("\rWriting... %d%%", kk++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    putchar('\r');
    sprintf(logbuf, "Distances (proportions) written to %s.\n", outname);
    logprintb();
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }

  if (!parallel_idx) {
    if (genome_output_full) {
      sprintf(tbuf, "%%%us%%%us%%%us%%%us RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO    IBS0    IBS1    IBS2  HOMHOM  HETHET\n", g_cg_max_person_fid_len, g_cg_max_person_iid_len, g_cg_max_person_fid_len, g_cg_max_person_iid_len);
    } else {
      sprintf(tbuf, "%%%us%%%us%%%us%%%us RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO\n", g_cg_max_person_fid_len, g_cg_max_person_iid_len, g_cg_max_person_fid_len, g_cg_max_person_iid_len);
    }
  }
  g_pct = 1;
  g_cg_gmcell = 0;
  g_cg_mdecell = 0;
  g_cg_marker_ct = marker_ct;
  g_cg_max_person_id_len = max_person_id_len;
  g_cg_max_paternal_id_len = max_paternal_id_len;
  g_cg_max_maternal_id_len = max_maternal_id_len;
  g_cg_pheno_nm = pheno_nm;
  g_cg_pheno_c = pheno_c;
  g_cg_founder_info = founder_info;
  g_cg_person_ids = person_ids;
  g_cg_paternal_ids = paternal_ids;
  g_cg_maternal_ids = maternal_ids;
  g_cg_genome_output_full = genome_output_full;
  g_cg_genome_ibd_unbounded = genome_ibd_unbounded;
  g_cg_pri = pri;
  g_cg_cur_line = 0;

  // signal to emitn() that this is the first line, and whether or not this is
  // a later part of a parallel write
  g_cg_indiv1idx = parallel_idx;
  g_cg_indiv2idx = 0;
  if (genome_output_gz) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome.gz");
    }
    parallel_compress(outname, calc_genome_emitn);
  } else {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome");
    }
    retval = write_uncompressed(outname, calc_genome_emitn);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }
  putchar('\r');
  sprintf(logbuf, "Finished writing %s.\n", outname);
  logprintb();
  wkspace_reset(wkspace_mark);
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
  calc_genome_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_genome_ret_1:
  gzclose_cond(gz_outfile);
  fclose_cond(outfile);
  return retval;
}

int32_t ld_prune(FILE* bedfile, int32_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_nm, uintptr_t* sex_male, int32_t ld_window_size, int32_t ld_window_kb, int32_t ld_window_incr, double ld_last_param, char* outname, char* outname_end, uint64_t calculation_type) {
  // todo: replace is_set with founder-sensitive check
  // for future consideration: chromosome-based multithread/parallel?
  FILE* outfile_in = NULL;
  FILE* outfile_out = NULL;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t indiv_ctl = (g_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t indiv_ct_mld = (g_indiv_ct + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
  int32_t indiv_ct_mld_m1 = (int)indiv_ct_mld - 1;
#ifdef __LP64__
  int32_t indiv_ct_mld_rem = (MULTIPLEX_LD / 192) - (indiv_ct_mld * MULTIPLEX_LD - g_indiv_ct) / 192;
#else
  int32_t indiv_ct_mld_rem = (MULTIPLEX_LD / 48) - (indiv_ct_mld * MULTIPLEX_LD - g_indiv_ct) / 48;
#endif
  uintptr_t indiv_ct_mld_long = indiv_ct_mld * (MULTIPLEX_LD / BITCT2);
  int32_t indiv_trail_ct = indiv_ct_mld_long - indiv_ctl * 2;
  int32_t retval = 0;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* geno = NULL;
  uintptr_t* pruned_arr;
  uint32_t* live_indices;
  uint32_t* start_arr;
  uint32_t marker_unfiltered_idx;
  uintptr_t marker_idx;
  int32_t pct;
  uint32_t pct_thresh;
  int32_t pairwise = ((calculation_type & CALC_LD_PRUNE_PAIRWISE) != 0);
  uint32_t window_unfiltered_start;
  uint32_t window_unfiltered_end;
  int32_t cur_window_size;
  int32_t old_window_size;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  uint32_t uii;
  uint32_t ujj;
  uint32_t cur_chrom;
  uint32_t chrom_end;
  uint32_t is_haploid;
  uint32_t is_x;
  double* marker_stdevs;
  unsigned char* loadbuf;
  uint32_t* missing_cts;
  uintptr_t window_max = 0;
  uintptr_t ulii;
  double dxx;
  double cov12;
  uint32_t fixed_non_missing_ct;
  uint32_t non_missing_ct;
  int32_t dp_result[3];
  double non_missing_recip;
#ifdef __LP64__
  __m128i* geno_fixed_vec_ptr;
  __m128i* geno_var_vec_ptr;
  __m128i* mask_fixed_vec_ptr;
  __m128i* mask_var_vec_ptr;
#else
  uintptr_t* geno_fixed_vec_ptr;
  uintptr_t* geno_var_vec_ptr;
  uintptr_t* mask_fixed_vec_ptr;
  uintptr_t* mask_var_vec_ptr;
#endif
  uintptr_t cur_exclude_ct;
  int32_t tot_exclude_ct = 0;
  int32_t prev_end;
  int32_t at_least_one_prune = 0;
  char* sptr;
  FILE* fptr;
  double* cov_matrix = NULL;
  double* new_cov_matrix = NULL;
  MATRIX_INVERT_BUF1_TYPE* irow = NULL;
  __CLPK_integer window_rem_li;
  __CLPK_integer old_window_rem_li;
  double* work = NULL;
  int32_t window_rem;
  uint32_t* idx_remap = NULL;
  double prune_ld_r2;

  if (ld_window_kb) {
    // determine maximum number of markers that may need to be loaded at once
    for (cur_chrom = 0; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
      if (chrom_exists(chrom_info_ptr, cur_chrom)) {
        uii = chrom_info_ptr->chrom_start[cur_chrom];
	chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
        do {
	  ujj = uii + 1;
	  while ((ujj < chrom_end) && (marker_pos[ujj] <= marker_pos[uii] + (1000 * ld_window_size))) {
	    ujj++;
	  }
          if (ujj - uii > window_max) {
	    window_max = ujj - uii;
	  }
	  uii++;
	} while (ujj < chrom_end);
      }
    }
  }
  if (pairwise) {
    prune_ld_r2 = sqrt(ld_last_param);
  } else {
    prune_ld_r2 = 0.999999;
  }

  window_unfiltered_start = ld_prune_next_valid_chrom_start(marker_exclude, 0, chrom_info_ptr, unfiltered_marker_ct);
  if (window_unfiltered_start == unfiltered_marker_ct) {
    if (pairwise) {
      logprint("Error: No valid markers for --indep-pairwise.\n");
    } else {
      logprint("Error: No valid markers for --indep.\n");
    }
    return RET_INVALID_FORMAT;
  }

  if (wkspace_alloc_ul_checked(&pruned_arr, unfiltered_marker_ctl * sizeof(intptr_t))) {
    goto ld_prune_ret_NOMEM;
  }

  memcpy(pruned_arr, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));

  if (!ld_window_kb) {
    window_max = ld_window_size;
  }
  ulii = window_max;
  if (wkspace_alloc_ui_checked(&live_indices, ulii * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&start_arr, ulii * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&marker_stdevs, ulii * sizeof(double)) ||
      wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4) ||
      wkspace_alloc_ul_checked(&geno, ulii * indiv_ct_mld_long * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_masks, ulii * indiv_ct_mld_long * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_mmasks, ulii * indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&missing_cts, g_indiv_ct * sizeof(int32_t))) {
    goto ld_prune_ret_NOMEM;
  }
  if (!pairwise) {
    if (wkspace_alloc_d_checked(&cov_matrix, window_max * window_max * sizeof(double)) ||
        wkspace_alloc_d_checked(&new_cov_matrix, window_max * window_max * sizeof(double)) ||
        wkspace_alloc_ui_checked(&idx_remap, window_max * sizeof(int32_t))) {
      goto ld_prune_ret_NOMEM;
    }

    irow = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(window_max * 2 * sizeof(MATRIX_INVERT_BUF1_TYPE));
    if (!irow) {
      goto ld_prune_ret_NOMEM;
    }

    if (window_max < 4) {
      ulii = 4;
    } else {
      ulii = window_max;
    }
    if (wkspace_alloc_d_checked(&work, ulii * window_max * sizeof(double))) {
      goto ld_prune_ret_NOMEM;
    }
  }
  do {
    prev_end = 0;
    ld_prune_start_chrom(ld_window_kb, &cur_chrom, &chrom_end, window_unfiltered_start, live_indices, start_arr, &window_unfiltered_end, ld_window_size, &cur_window_size, unfiltered_marker_ct, pruned_arr, chrom_info_ptr, marker_pos, &is_haploid, &is_x);
    old_window_size = 1;
    if (cur_window_size > 1) {
      for (ulii = 0; ulii < (uintptr_t)cur_window_size; ulii++) {
	if (fseeko(bedfile, bed_offset + (live_indices[ulii] * unfiltered_indiv_ct4), SEEK_SET)) {
	  goto ld_prune_ret_READ_FAIL;
	}
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto ld_prune_ret_READ_FAIL;
	}
        missing_cts[ulii] = ld_process_load(loadbuf, &(geno[ulii * indiv_ct_mld_long]), &(g_masks[ulii * indiv_ct_mld_long]), &(g_mmasks[ulii * indiv_ctl]), &(marker_stdevs[ulii]), unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, indiv_ctl, indiv_trail_ct, is_haploid, is_x, sex_male);
      }
    }
    pct = 1;
    pct_thresh = window_unfiltered_start + ((int64_t)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100;
    cur_exclude_ct = 0;
    while ((window_unfiltered_start < chrom_end) || (cur_window_size > 1)) {
      if (cur_window_size > 1) {
	for (ii = 0; ii < cur_window_size; ii++) {
	  if (marker_stdevs[ii] == 0.0) {
	    set_bit(pruned_arr, live_indices[ii], &cur_exclude_ct);
	  }
	}
	do {
	  at_least_one_prune = 0;
	  for (ii = 0; ii < cur_window_size - 1; ii++) {
	    if (is_set(pruned_arr, live_indices[ii])) {
	      continue;
	    }
	    fixed_non_missing_ct = g_indiv_ct - missing_cts[ii];
#ifdef __LP64__
	    geno_fixed_vec_ptr = (__m128i*)(&(geno[ii * indiv_ct_mld_long]));
	    mask_fixed_vec_ptr = (__m128i*)(&(g_masks[ii * indiv_ct_mld_long]));
#else
	    geno_fixed_vec_ptr = &(geno[ii * indiv_ct_mld_long]);
	    mask_fixed_vec_ptr = &(g_masks[ii * indiv_ct_mld_long]);
#endif
	    jj = ii + 1;
	    while (live_indices[jj] < start_arr[ii]) {
	      jj++;
	      if (jj == cur_window_size) {
		break;
	      }
	    }
	    for (; jj < cur_window_size; jj++) {
	      if (is_set(pruned_arr, live_indices[jj])) {
		continue;
	      }
#ifdef __LP64__
	      geno_var_vec_ptr = (__m128i*)(&(geno[jj * indiv_ct_mld_long]));
	      mask_var_vec_ptr = (__m128i*)(&(g_masks[jj * indiv_ct_mld_long]));
#else
	      geno_var_vec_ptr = &(geno[jj * indiv_ct_mld_long]);
	      mask_var_vec_ptr = &(g_masks[jj * indiv_ct_mld_long]);
#endif

	      non_missing_ct = fixed_non_missing_ct - missing_cts[jj] + sparse_intersection_ct(&(g_mmasks[ii * indiv_ctl]), &(g_mmasks[jj * indiv_ctl]), indiv_ctl);
	      dp_result[0] = g_indiv_ct;
	      // reversed from what I initially thought because I'm passing the
	      // jj-associated buffers before the ii-associated ones.
	      dp_result[1] = -fixed_non_missing_ct;
	      dp_result[2] = missing_cts[jj] - g_indiv_ct;
#ifdef __LP64__
	      for (kk = 0; kk < indiv_ct_mld_m1; kk++) {
		ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), dp_result, MULTIPLEX_LD / 192);
		geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT]);
		mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT]);
	      }
	      ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT)]), dp_result, indiv_ct_mld_rem);
	      geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT]);
	      mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT]);
#else
	      for (kk = 0; kk < indiv_ct_mld_m1; kk++) {
		ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), dp_result, MULTIPLEX_LD / 48);
		geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
		mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
	      }
	      ld_dot_prod(geno_var_vec_ptr, &(geno_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), mask_var_vec_ptr, &(mask_fixed_vec_ptr[kk * (MULTIPLEX_LD / BITCT2)]), dp_result, indiv_ct_mld_rem);
	      geno_var_vec_ptr = &(geno_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
	      mask_var_vec_ptr = &(mask_var_vec_ptr[MULTIPLEX_LD / BITCT2]);
#endif
	      non_missing_recip = 1.0 / ((double)non_missing_ct);
	      cov12 = non_missing_recip * (dp_result[0] - (non_missing_recip * dp_result[1]) * dp_result[2]);
	      // r-squared
	      dxx = cov12 / (marker_stdevs[ii] * marker_stdevs[jj]);
	      if (!pairwise) {
		cov_matrix[ii * window_max + jj] = dxx;
	      }
	      if (fabs(dxx) > prune_ld_r2) {
		at_least_one_prune = 1;
		// remove marker with lower MAF
		if (get_maf(set_allele_freqs[live_indices[ii]]) < get_maf(set_allele_freqs[live_indices[jj]])) {
		  set_bit(pruned_arr, live_indices[ii], &cur_exclude_ct);
		} else {
		  set_bit(pruned_arr, live_indices[jj], &cur_exclude_ct);
		  jj++;
		  while (jj < cur_window_size) {
		    if (!is_set(pruned_arr, live_indices[jj])) {
		      break;
		    }
		    jj++;
		  }
		  if (jj < cur_window_size) {
		    start_arr[ii] = live_indices[jj];
		  }
		}
		break;
	      }
	    }
	    if (jj == cur_window_size) {
	      start_arr[ii] = window_unfiltered_end;
	    }
	  }
	} while (at_least_one_prune);
	if (!pairwise) {
	  window_rem = 0;
	  old_window_rem_li = 0;
	  for (ii = 0; ii < old_window_size; ii++) {
	    if (is_set(pruned_arr, live_indices[ii])) {
	      continue;
	    }
            idx_remap[window_rem++] = ii;
	  }
	  old_window_rem_li = window_rem;
	  for (; ii < cur_window_size; ii++) {
	    if (is_set(pruned_arr, live_indices[ii])) {
	      continue;
	    }
            idx_remap[window_rem++] = ii;
	  }
	  while (window_rem > 1) {
	    new_cov_matrix[0] = 1.0;
	    for (ii = 1; ii < window_rem; ii++) {
	      kk = idx_remap[ii];
	      for (jj = 0; jj < ii; jj++) {
		dxx = cov_matrix[idx_remap[jj] * window_max + kk];
		new_cov_matrix[jj * window_rem + ii] = dxx;
		new_cov_matrix[ii * window_rem + jj] = dxx;
	      }
	      new_cov_matrix[ii * (window_rem + 1)] = 1.0;
	    }
	    window_rem_li = window_rem;
	    jj = invert_matrix_trunc_singular(window_rem_li, new_cov_matrix, irow, work, old_window_rem_li);
	    while (jj) {
	      if (jj == -1) {
		goto ld_prune_ret_NOMEM;
	      }
              set_bit(pruned_arr, live_indices[idx_remap[jj]], &cur_exclude_ct);
	      window_rem--;
	      for (ii = jj; ii < window_rem; ii++) {
		idx_remap[ii] = idx_remap[ii + 1];
	      }
	      new_cov_matrix[0] = 1.0;
	      for (ii = 1; ii < window_rem; ii++) {
		kk = idx_remap[ii];
		for (jj = 0; jj < ii; jj++) {
		  dxx = cov_matrix[idx_remap[jj] * window_max + kk];
		  new_cov_matrix[jj * window_rem + ii] = dxx;
		  new_cov_matrix[ii * window_rem + jj] = dxx;
		}
		new_cov_matrix[ii * (window_rem + 1)] = 1.0;
	      }
              window_rem_li = window_rem;
	      jj = invert_matrix_trunc_singular(window_rem_li, new_cov_matrix, irow, work, old_window_rem_li);
	    }
	    dxx = new_cov_matrix[0];
	    jj = 0;
	    for (ii = 1; ii < window_rem; ii++) {
              if (new_cov_matrix[ii * (window_rem + 1)] > dxx) {
		dxx = new_cov_matrix[ii * (window_rem + 1)];
		jj = ii;
	      }
	    }
	    if (dxx > ld_last_param) {
	      set_bit(pruned_arr, live_indices[idx_remap[jj]], &cur_exclude_ct);
	      window_rem--;
	      if (idx_remap[jj] < (uint32_t)old_window_size) {
		old_window_rem_li--;
	      }
	      for (ii = jj; ii < window_rem; ii++) {
                idx_remap[ii] = idx_remap[ii + 1];
	      }
	    } else {
	      // break out
	      window_rem = 1;
	    }
	  }
	}
      }
      for (ii = 0; ii < ld_window_incr; ii++) {
	while (is_set(marker_exclude, window_unfiltered_start)) {
	  if (window_unfiltered_start == chrom_end) {
	    break;
	  }
	  window_unfiltered_start++;
	}
	if (window_unfiltered_start == chrom_end) {
	  break;
	}
	window_unfiltered_start++;
      }
      if (window_unfiltered_start == chrom_end) {
	break;
      }
      if (window_unfiltered_start >= pct_thresh) {
	pct = (((int64_t)(window_unfiltered_start - chrom_info_ptr->chrom_start[cur_chrom])) * 100) / (chrom_end - chrom_info_ptr->chrom_start[cur_chrom]);
	printf("\r%d%%", pct++);
	fflush(stdout);
	pct_thresh = chrom_info_ptr->chrom_start[cur_chrom] + (((int64_t)pct * (chrom_end - chrom_info_ptr->chrom_start[cur_chrom])) / 100);
      }
      jj = 0;
      // copy back previously loaded/computed results
      while (live_indices[jj] < window_unfiltered_start) {
	jj++;
	if (jj == cur_window_size) {
	  break;
	}
      }
      for (ii = 0; jj < cur_window_size; jj++) {
	if (is_set(pruned_arr, live_indices[jj])) {
	  continue;
	}
	memcpy(&(geno[ii * indiv_ct_mld_long]), &(geno[jj * indiv_ct_mld_long]), indiv_ct_mld_long * sizeof(intptr_t));
	memcpy(&(g_masks[ii * indiv_ct_mld_long]), &(g_masks[jj * indiv_ct_mld_long]), indiv_ct_mld_long * sizeof(intptr_t));
	memcpy(&(g_mmasks[ii * indiv_ctl]), &(g_mmasks[jj * indiv_ctl]), indiv_ctl * sizeof(intptr_t));
	marker_stdevs[ii] = marker_stdevs[jj];
	live_indices[ii] = live_indices[jj];
	start_arr[ii] = start_arr[jj];
	missing_cts[ii] = missing_cts[jj];
	if (!pairwise) {
	  for (kk = 0; kk < ii; kk++) {
	    cov_matrix[kk * window_max + ii] = cov_matrix[idx_remap[kk] * window_max + jj];
	  }
	  idx_remap[ii] = jj;
	}
	ii++;
      }

      prev_end = ii;
      cur_window_size = ii;
      if (ld_window_kb) {
	jj = 0;
	while ((window_unfiltered_end + jj < chrom_end) && (marker_pos[window_unfiltered_end + jj] <= marker_pos[window_unfiltered_start] + (1000 * ld_window_size))) {
	  jj++;
	}
      } else {
	jj = ld_window_incr;
      }
      old_window_size = cur_window_size;
      for (ii = 0; ii < jj; ii++) {
	while ((window_unfiltered_end < chrom_end) && is_set(marker_exclude, window_unfiltered_end)) {
	  window_unfiltered_end++;
	}
	if (window_unfiltered_end < chrom_end) {
	  live_indices[cur_window_size] = window_unfiltered_end;
	  if (cur_window_size > prev_end) {
	    start_arr[cur_window_size - 1] = window_unfiltered_end;
	  }
	  if (fseeko(bedfile, bed_offset + (window_unfiltered_end * unfiltered_indiv_ct4), SEEK_SET)) {
	    goto ld_prune_ret_READ_FAIL;
	  }
	  if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	    goto ld_prune_ret_READ_FAIL;
	  }
	  missing_cts[cur_window_size] = ld_process_load(loadbuf, &(geno[cur_window_size * indiv_ct_mld_long]), &(g_masks[cur_window_size * indiv_ct_mld_long]), &(g_mmasks[cur_window_size * indiv_ctl]), &(marker_stdevs[cur_window_size]), unfiltered_indiv_ct, indiv_exclude, g_indiv_ct, indiv_ctl, indiv_trail_ct, is_haploid, is_x, sex_male);
	  cur_window_size++;
	  window_unfiltered_end++;
	}
      }
      if (cur_window_size > prev_end) {
	start_arr[cur_window_size] = window_unfiltered_end;
      }
    }
    ii = get_marker_chrom(chrom_info_ptr, window_unfiltered_start - 1);
    putchar('\r');
    sprintf(logbuf, "Pruned %" PRIuPTR " markers from chromosome %d, leaving %" PRIuPTR ".\n", cur_exclude_ct, ii, chrom_info_ptr->chrom_end[ii] - chrom_info_ptr->chrom_start[ii] - cur_exclude_ct);
    logprintb();
    tot_exclude_ct += cur_exclude_ct;

    // advance chromosomes as necessary
    window_unfiltered_start = ld_prune_next_valid_chrom_start(pruned_arr, window_unfiltered_start, chrom_info_ptr, unfiltered_marker_ct);
  } while (window_unfiltered_start < unfiltered_marker_ct);

  sprintf(logbuf, "Pruning complete.  %d of %d markers removed.\n", tot_exclude_ct, marker_ct);
  logprintb();
  strcpy(outname_end, ".prune.in");
  if (fopen_checked(&outfile_in, outname, "w")) {
    goto ld_prune_ret_OPEN_FAIL;
  }
  strcpy(outname_end, ".prune.out");
  if (fopen_checked(&outfile_out, outname, "w")) {
    goto ld_prune_ret_OPEN_FAIL;
  }
  marker_unfiltered_idx = 0;
  marker_idx = 0;
  pct = 1;
  ii = 0;
  for (cur_chrom = 1; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
    if (!(chrom_info_ptr->chrom_mask & (1LLU << cur_chrom))) {
      continue;
    }
    if (chrom_info_ptr->chrom_end[cur_chrom]) {
      ii += chrom_info_ptr->chrom_end[cur_chrom] - chrom_info_ptr->chrom_start[cur_chrom];
    }
  }
  pct_thresh = ((int64_t)pct * ii) / 100;
  for (cur_chrom = 1; cur_chrom < MAX_POSSIBLE_CHROM; cur_chrom++) {
    chrom_end = chrom_info_ptr->chrom_end[cur_chrom];
    if (!chrom_end) {
      continue;
    }
    marker_unfiltered_idx = chrom_info_ptr->chrom_start[cur_chrom];
    for (; marker_unfiltered_idx < chrom_end; marker_unfiltered_idx++) {
      if (!is_set(marker_exclude, marker_unfiltered_idx)) {
	sptr = &(marker_ids[marker_unfiltered_idx * max_marker_id_len]);
	fptr = is_set(pruned_arr, marker_unfiltered_idx)? outfile_out : outfile_in;
	fwrite(sptr, 1, strlen(sptr), fptr);
	if (putc('\n', fptr) == EOF) {
	  goto ld_prune_ret_WRITE_FAIL;
	}
      }
      marker_idx++;
      if (marker_idx == pct_thresh) {
	printf("\rWriting... %d%%", pct);
	fflush(stdout);
	pct = ((int64_t)marker_idx * 100) / ii + 1;
        pct_thresh = ((int64_t)pct * ii) / 100;
      }
    }
  }
  if (fclose_null(&outfile_in)) {
    goto ld_prune_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile_out)) {
    goto ld_prune_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  putchar('\r');
  sprintf(logbuf, "Marker lists written to %s.prune.in and %s.prune.out.\n", outname, outname);
  logprintb();

  while (0) {
  ld_prune_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ld_prune_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ld_prune_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ld_prune_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile_in);
  fclose_cond(outfile_out);
  wkspace_reset(wkspace_mark);
  return retval;
}

inline void rel_cut_arr_dec(int32_t* rel_ct_arr_elem, uint32_t* exactly_one_rel_ct_ptr) {
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

int32_t do_rel_cutoff(uint64_t calculation_type, double rel_cutoff, double* rel_ibc, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, char* person_ids, uintptr_t max_person_id_len) {
  int32_t indivs_excluded = 0;
  uint32_t exactly_one_rel_ct = 0;
  unsigned char* wkspace_mark = wkspace_base;
  double* dist_ptr = g_rel_dists;
  double* dptr2;
  double* dptr3;
  double* dptr4;
  uint32_t* giptr;
  uint32_t* giptr2;
  // number of too-close relations, -1 if excluded
  int32_t* rel_ct_arr;
  uintptr_t indiv_idx;
  uint32_t uii;
  uintptr_t ulii;
  int32_t kk;
  int32_t mm;
  int32_t retval;
  
  // Algorithm:
  // - Whenever there is at least one individual with exactly one
  // remaining too-close relation, prune the other side of that
  // relationship, because doing so is never suboptimal.
  // - Otherwise, there's no efficient rule that is always optimal
  // (assuming P != NP, anyway), so we use a simple heuristic: prune the
  // first individual with the largest number of remaining too-close
  // relationships.

  if (wkspace_alloc_i_checked(&rel_ct_arr, g_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_int_zero(rel_ct_arr, g_indiv_ct);
  for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
    for (uii = 0; uii < indiv_idx; uii++) {
      if (*dist_ptr++ > rel_cutoff) {
	rel_ct_arr[indiv_idx] += 1;
	rel_ct_arr[uii] += 1;
      }
    }
  }
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    if (rel_ct_arr[indiv_idx] == 1) {
      exactly_one_rel_ct++;
    }
  }
  while (1) {
    kk = 0;
    if (exactly_one_rel_ct) {
      // there is at least one individual with exactly one too-close
      // relation left, find the first one
      while (rel_ct_arr[kk] != 1) {
	kk++;
      }
      // and now find the identity of the other side
      dist_ptr = &(g_rel_dists[((intptr_t)kk * (kk - 1)) / 2]);
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
	  dist_ptr = &(g_rel_dists[((intptr_t)mm * (mm - 1)) / 2 + kk]);
	} while (*dist_ptr <= rel_cutoff);
	*dist_ptr = 0.0;
      }
      rel_ct_arr[kk] = 0;
      exactly_one_rel_ct--;
      if (rel_ct_arr[mm] == 1) {
        // speed up the easy case
	exactly_one_rel_ct--;
	rel_ct_arr[mm] = -1;
	indivs_excluded++;
	continue;
      }
    } else {
      // find identity of first individual with maximum number of
      // remaining too-close relations
      // kk is highest too-close pair count so far
      mm = -1; // associated individual index
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	if (rel_ct_arr[indiv_idx] > kk) {
	  kk = rel_ct_arr[indiv_idx];
	  mm = indiv_idx;
	}
      }
      // no too-close relations left at all, we're done
      if (mm == -1) {
	break;
      }
    }
    dist_ptr = &(g_rel_dists[((intptr_t)mm * (mm - 1)) / 2]);
    for (kk = 0; kk < mm; kk++) {
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[kk]), &exactly_one_rel_ct);
      }
      dist_ptr++;
    }
    for (ulii = mm + 1; ulii < g_indiv_ct; ulii++) {
      dist_ptr = &(g_rel_dists[(ulii * (ulii - 1)) / 2 + mm]);
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[ulii]), &exactly_one_rel_ct);
      }
    }
    rel_ct_arr[mm] = -1;
    indivs_excluded++;
  }
  exclude_multi(indiv_exclude, rel_ct_arr, g_indiv_ct, indiv_exclude_ct_ptr);
  if (indivs_excluded) {
    dist_ptr = g_rel_dists; // write
    dptr2 = g_rel_dists; // read
    dptr3 = rel_ibc; // write
    dptr4 = rel_ibc; // read
    giptr = g_missing_dbl_excluded; // write
    giptr2 = g_missing_dbl_excluded; // read
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      if (rel_ct_arr[indiv_idx] != -1) {
	if (calculation_type & CALC_IBC) {
	  dptr3[g_indiv_ct] = dptr4[g_indiv_ct];
	  dptr3[g_indiv_ct * 2] = dptr4[g_indiv_ct * 2];
	}
	*dptr3 = *dptr4++;
	dptr3++;
	for (uii = 0; uii < indiv_idx; uii++) {
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
	dptr2 = &(dptr2[indiv_idx]);
	giptr2 = &(giptr2[indiv_idx]);
      }
    }
    g_indiv_ct -= indivs_excluded;
    if (calculation_type & CALC_IBC) {
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	*dptr3++ = *dptr4++;
      }
      dptr4 = &(dptr4[indivs_excluded]);
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	*dptr3++ = *dptr4++;
      }
    }
    giptr = g_indiv_missing_unwt;
    giptr2 = g_indiv_missing_unwt;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct + indivs_excluded; indiv_idx++) {
      if (rel_ct_arr[indiv_idx] != -1) {
	*giptr = *giptr2++;
	giptr++;
      } else {
	giptr2++;
      }
    }
  }
  sprintf(logbuf, "%d %s excluded by --rel-cutoff.\n", indivs_excluded, species_str(indivs_excluded));
  logprintb();
  if (!(calculation_type & (CALC_RELATIONSHIP | CALC_GDISTANCE_MASK))) {
    strcpy(outname_end, ".rel.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
    sprintf(logbuf, "Remaining individual IDs written to %s.\n", outname);
    logprintb();
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

uint32_t g_rcb_row;
uint32_t g_rcb_col;
uint32_t g_rcb_new_row;
uint32_t g_rcb_new_col;
uint32_t g_rcb_indiv_ct;
uint64_t g_rcb_progress;
uint64_t g_rcb_hundredth;
FILE* g_rcb_in_binfile;
FILE* g_rcb_in_bin_nfile;
gzFile g_rcb_cur_gzfile;
int32_t* g_rcb_rel_ct_arr;

uint32_t rel_cutoff_batch_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  char wbuf[16];
  char* cptr;
  uint32_t wbuf_ct;
  uint32_t uii;
  while (g_rcb_row < g_rcb_indiv_ct) {
    if (g_rcb_rel_ct_arr[g_rcb_row] == -1) {
      for (uii = 0; uii <= g_rcb_row; uii++) {
	gzgets(g_rcb_cur_gzfile, tbuf, MAXLINELEN);
      }
    } else {
      cptr = uint32_writex(wbuf, g_rcb_new_row, '\t');
      wbuf_ct = (uintptr_t)(cptr - wbuf);
      while (g_rcb_col <= g_rcb_row) {
        gzgets(g_rcb_cur_gzfile, tbuf, MAXLINELEN);
	if (g_rcb_rel_ct_arr[g_rcb_col++] != -1) {
	  cptr = next_item_mult(tbuf, 2);
          uii = strlen(cptr);
          sptr_cur = memcpya(uint32_writex(memcpya(sptr_cur, wbuf, wbuf_ct), ++g_rcb_new_col, '\t'), cptr, uii);
	  if (sptr_cur >= readbuf_end) {
	    goto rel_cutoff_batch_emitn_ret;
	  }
	}
      }
      g_rcb_new_row++;
      g_rcb_new_col = 0;
    }
    g_rcb_row++;
    g_rcb_progress += g_rcb_row;
    if (g_rcb_progress >= g_pct * g_rcb_hundredth) {
      if (g_pct > 10) {
	putchar('\b');
      }
      g_pct = 1 + (g_rcb_progress / g_rcb_hundredth);
      printf("\b\b%u%%", g_pct - 1);
      fflush(stdout);
    }
    g_rcb_col = 0;
  }  
 rel_cutoff_batch_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t rel_cutoff_batch_rbin_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  char wbuf[16];
  char* cptr;
  uint32_t wbuf_ct;
  uint32_t uii;
  float fxx;

  while (g_rcb_row < g_rcb_indiv_ct) {
    if (g_rcb_rel_ct_arr[g_rcb_row] == -1) {
      fseeko(g_rcb_in_binfile, (g_rcb_row + 1) * sizeof(float), SEEK_CUR);
      fseeko(g_rcb_in_bin_nfile, (g_rcb_row + 1) * sizeof(float), SEEK_CUR);
    } else {
      cptr = uint32_writex(wbuf, g_rcb_new_row, '\t');
      wbuf_ct = (uintptr_t)(cptr - wbuf);
      while (g_rcb_col <= g_rcb_row) {
	if (g_rcb_rel_ct_arr[g_rcb_col] == -1) {
	  uii = g_rcb_col;
	  while (++g_rcb_col <= g_rcb_row) {
	    if (g_rcb_rel_ct_arr[g_rcb_col] != -1) {
	      break;
	    }
	  }
	  fseeko(g_rcb_in_binfile, (g_rcb_col - uii) * sizeof(float), SEEK_CUR);
	  fseeko(g_rcb_in_bin_nfile, (g_rcb_col - uii) * sizeof(float), SEEK_CUR);
	  if (g_rcb_col > g_rcb_row) {
	    break;
	  }
	}
	sptr_cur = uint32_writex(memcpya(sptr_cur, wbuf, wbuf_ct), ++g_rcb_new_col, '\t');
	fread(&fxx, 4, 1, g_rcb_in_bin_nfile);
	sptr_cur = uint32_writex(sptr_cur, (int32_t)fxx, '\t');
	fread(&fxx, 4, 1, g_rcb_in_binfile);
	sptr_cur = float_e_writex(sptr_cur, fxx, '\n');
	g_rcb_col++;
	if (sptr_cur >= readbuf_end) {
	  goto rel_cutoff_batch_rbin_emitn_ret;
	}
      }
      g_rcb_new_row++;
      g_rcb_new_col = 0;
    }
    g_rcb_row++;
    g_rcb_progress += g_rcb_row;
    if (g_rcb_progress >= g_pct * g_rcb_hundredth) {
      if (g_pct > 10) {
	putchar('\b');
      }
      g_pct = 1 + (g_rcb_progress / g_rcb_hundredth);
      printf("\b\b%u%%", g_pct - 1);
      fflush(stdout);
    }
    g_rcb_col = 0;
  }  
 rel_cutoff_batch_rbin_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t rel_cutoff_batch(uint32_t load_grm_bin, char* grmname, char* outname, char* outname_end, double rel_cutoff, int32_t rel_calc_type) {
  // Specialized --rel-cutoff usable on larger files.
  char* grmname_end = (char*)memchr(grmname, 0, FNAMESIZE);
  uintptr_t indiv_ct = 0;
  FILE* idfile = NULL;
  FILE* outfile = NULL;
  FILE* out_bin_nfile = NULL;
  gzFile cur_gzfile = NULL;
  gzFile gz_outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uint32_t indivs_excluded = 0;
  uint32_t exactly_one_rel_ct = 0;
  uintptr_t* compact_rel_table;
  uintptr_t* rtptr;
  char* bufptr;
  uint64_t ullii;
  uint64_t ulljj;
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
  uintptr_t indiv_idx;
  int32_t* rel_ct_arr;
  float rel_cutoff_f;
  float fxx;
  float fyy;
  double dxx;
  int32_t retval;
  int32_t kk;
  int32_t mm;
  int32_t cur_prune;
  g_rcb_in_binfile = NULL;
  g_rcb_in_bin_nfile = NULL;
  g_rcb_cur_gzfile = NULL;
  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(&idfile, grmname, "r")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, idfile)) {
    if (is_eoln_kns(*(skip_initial_spaces(tbuf)))) {
      continue;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in .grm.id file.\n");
      goto rel_cutoff_batch_ret_INVALID_FORMAT_1;
    }
    indiv_ct++;
  }
  if (!feof(idfile)) {
    goto rel_cutoff_batch_ret_READ_FAIL;
  }
  fclose_null(&idfile);
  ullii = indiv_ct;
  ullii = ((ullii * (ullii - 1)) / 2 + BITCT - 1) / BITCT;
#ifndef __LP64__
  if (ullii >= 0x20000000) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
#endif
  tot_words = ullii;
  if (wkspace_alloc_ul_checked(&compact_rel_table, tot_words * sizeof(intptr_t))) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
  fill_ulong_zero(compact_rel_table, tot_words);
  if (wkspace_alloc_i_checked(&rel_ct_arr, indiv_ct * sizeof(int32_t))) {
    goto rel_cutoff_batch_ret_NOMEM;
  }
  fill_int_zero(rel_ct_arr, indiv_ct);

  fputs("Reading... 0%", stdout);
  fflush(stdout);
  words_left = tot_words;
  rtptr = compact_rel_table;
  row = 0;
  col = 0;
  if (load_grm_bin) {
    memcpy(grmname_end, ".grm.bin", 9);
    if (fopen_checked(&g_rcb_in_binfile, grmname, "rb")) {
      goto rel_cutoff_batch_ret_OPEN_FAIL;
    }
    rel_cutoff_f = (float)rel_cutoff;
    for (g_pct = 1; g_pct <= 100; g_pct++) {
      wl_floor = (((uint64_t)tot_words) * (100 - g_pct)) / 100;
      while (words_left > wl_floor) {
	cur_word = 0;
	if (--words_left) {
	  inword_bound = BITCT;
	} else {
	  // only indiv_ct mod (BITCT * 2) matters for remainder
	  uii = indiv_ct & (BITCT * 2 - 1);
	  inword_bound = ((uii * (uii - 1)) / 2) & (BITCT - 1);
	  if (!inword_bound) {
	    inword_bound = BITCT;
	  }
	}
	for (inword_idx = 0; inword_idx < inword_bound; inword_idx++) {
	  if (!fread(&fxx, 4, 1, g_rcb_in_binfile)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	  if (row == col) {
	    row++;
	    col = 0;
	    if (!fread(&fxx, 4, 1, g_rcb_in_binfile)) {
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
    }
    fclose_null(&g_rcb_in_binfile);
  } else {
    memcpy(grmname_end, ".grm.gz", 8);
    if (gzopen_checked(&cur_gzfile, grmname, "rb")) {
      goto rel_cutoff_batch_ret_OPEN_FAIL;
    }
    if (gzbuffer(cur_gzfile, 131072)) {
      goto rel_cutoff_batch_ret_NOMEM;
    }
    for (g_pct = 1; g_pct <= 100; g_pct++) {
      wl_floor = (((uint64_t)tot_words) * (100 - g_pct)) / 100;
      while (words_left > wl_floor) {
	cur_word = 0;
	if (--words_left) {
	  inword_bound = BITCT;
	} else {
	  // only indiv_ct mod (BITCT * 2) matters for remainder
	  uii = indiv_ct & (BITCT * 2 - 1);
	  inword_bound = ((uii * (uii - 1)) / 2) & (BITCT - 1);
	  if (!inword_bound) {
	    inword_bound = BITCT;
	  }
	}
	for (inword_idx = 0; inword_idx < inword_bound; inword_idx++) {
	  if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
	    goto rel_cutoff_batch_ret_READ_FAIL;
	  }
	  if (row == col) {
	    row++;
	    col = 0;
	    if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
	      goto rel_cutoff_batch_ret_READ_FAIL;
	    }
	  }
	  bufptr = next_item_mult(tbuf, 3);
	  if (no_more_items_kns(bufptr)) {
	    goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
	  }
	  if (sscanf(bufptr, "%lg", &dxx) != 1) {
	    goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
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
      if (g_pct < 100) {
	if (g_pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", g_pct);
	fflush(stdout);
      }
    }
    if (!gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
      goto rel_cutoff_batch_ret_READ_FAIL;
    }
    if (gzgets(cur_gzfile, tbuf, MAXLINELEN)) {
      goto rel_cutoff_batch_ret_INVALID_FORMAT_2;
    }
    gzclose(cur_gzfile);
    cur_gzfile = NULL;
  }
  putchar('\r');
  sprintf(logbuf, "%s read complete.  Pruning.\n", grmname);
  logprintb();

  // would prefer to just call do_rel_cutoff(), but unfortunately that
  // interferes with the intended "handle extra-large datasets" mission of this
  // function, which necessitates a compact bit representation of the
  // relationship matrix... fortunately, the algorithm is pretty simple.
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if (rel_ct_arr[indiv_idx] == 1) {
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
	indivs_excluded++;
	continue;
      }
    } else {
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	if (rel_ct_arr[indiv_idx] > kk) {
	  kk = rel_ct_arr[indiv_idx];
	  cur_prune = indiv_idx;
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

    for (uljj = cur_prune + 1; uljj < indiv_ct; uljj++) {
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
    indivs_excluded++;
  }

  memcpy(grmname_end, ".grm.id", 8);
  if (fopen_checked(&idfile, grmname, "r")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }

  memcpy(outname_end, ".grm.id", 8);
  if (fopen_checked(&outfile, outname, "w")) {
    goto rel_cutoff_batch_ret_OPEN_FAIL;
  }

  for (indiv_idx = 0; indiv_idx < indiv_ct;) {
    if (fgets(tbuf, MAXLINELEN, idfile) == NULL) {
      goto rel_cutoff_batch_ret_READ_FAIL;
    }
    if (is_eoln_kns(*(skip_initial_spaces(tbuf)))) {
      continue;
    }
    if (rel_ct_arr[indiv_idx] != -1) {
      if (fputs(tbuf, outfile) == EOF) {
	goto rel_cutoff_batch_ret_WRITE_FAIL;
      }
    }
    indiv_idx++;
  }

  fclose_null(&idfile);
  fclose_null(&outfile);

  sprintf(logbuf, "%d %s excluded by --rel-cutoff.\n", indivs_excluded, species_str(indivs_excluded));
  logprintb();
  sprintf(logbuf, "Remaining individual IDs written to %s.\n", outname);
  logprintb();
  if (rel_calc_type & (REL_CALC_GRM | REL_CALC_GRM_BIN)) {
    g_pct = 1;
    g_rcb_row = 0;
    g_rcb_col = 0;
    g_rcb_new_row = 1;
    g_rcb_new_col = 0;
    g_rcb_progress = 0;
    g_rcb_hundredth = 1 + ((((uint64_t)indiv_ct) * (indiv_ct - 1)) / 100);
    if (load_grm_bin) {
      memcpy(grmname_end, ".grm.bin", 9);
      if (fopen_checked(&g_rcb_in_binfile, grmname, "rb")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      memcpy(grmname_end, ".grm.N.bin", 11);
      if (fopen_checked(&g_rcb_in_bin_nfile, grmname, "rb")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
    } else {
      memcpy(grmname_end, ".grm.gz", 8);
      if (gzopen_checked(&g_rcb_cur_gzfile, grmname, "rb")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      if (gzbuffer(g_rcb_cur_gzfile, 131072)) {
	goto rel_cutoff_batch_ret_NOMEM;
      }
    }
    fputs("Rewriting matrix... 0%", stdout);
    fflush(stdout);
    if (rel_calc_type & REL_CALC_GRM) {
      g_rcb_indiv_ct = indiv_ct;
      g_rcb_rel_ct_arr = rel_ct_arr;
      if (load_grm_bin) {
	if (rel_calc_type & REL_CALC_GZ) {
	  memcpy(outname_end, ".grm.gz", 8);
	  parallel_compress(outname, rel_cutoff_batch_rbin_emitn);
	} else {
	  memcpy(outname_end, ".grm", 5);
	  retval = write_uncompressed(outname, rel_cutoff_batch_rbin_emitn);
	  if (retval) {
	    goto rel_cutoff_batch_ret_1;
	  }
	}
      } else {
	if (rel_calc_type & REL_CALC_GZ) {
	  memcpy(outname_end, ".grm.gz", 8);
	  parallel_compress(outname, rel_cutoff_batch_emitn);
	} else {
	  memcpy(outname_end, ".grm", 5);
	  retval = write_uncompressed(outname, rel_cutoff_batch_emitn);
	  if (retval) {
	    goto rel_cutoff_batch_ret_1;
	  }
	}
      }
    } else {
      memcpy(outname_end, ".grm.N.bin", 11);
      if (fopen_checked(&out_bin_nfile, outname, "wb")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      memcpy(outname_end, ".grm.bin", 9);
      if (fopen_checked(&outfile, outname, "wb")) {
	goto rel_cutoff_batch_ret_OPEN_FAIL;
      }
      while (g_rcb_row < indiv_ct) {
	if (rel_ct_arr[g_rcb_row] == -1) {
	  if (load_grm_bin) {
	    fseeko(g_rcb_in_binfile, (g_rcb_row + 1) * sizeof(float), SEEK_CUR);
	    fseeko(g_rcb_in_bin_nfile, (g_rcb_row + 1) * sizeof(float), SEEK_CUR);
	  } else {
	    for (uii = 0; uii <= g_rcb_row; uii++) {
	      gzgets(g_rcb_cur_gzfile, tbuf, MAXLINELEN);
	    }
	  }
	} else {
	  if (load_grm_bin) {
	    while (g_rcb_col <= g_rcb_row) {
	      if (rel_ct_arr[g_rcb_col] == -1) {
		uii = g_rcb_col;
		while (++g_rcb_col <= g_rcb_row) {
		  if (rel_ct_arr[g_rcb_col] != -1) {
		    break;
		  }
		}
		fseeko(g_rcb_in_bin_nfile, (g_rcb_col - uii) * sizeof(float), SEEK_CUR);
		fseeko(g_rcb_in_bin_nfile, (g_rcb_col - uii) * sizeof(float), SEEK_CUR);
		if (g_rcb_col > g_rcb_row) {
		  break;
		}
	      }
	      fread(&fxx, 4, 1, g_rcb_in_bin_nfile);
	      fwrite(&fxx, 4, 1, out_bin_nfile);
	      fread(&fxx, 4, 1, g_rcb_in_binfile);
	      fwrite(&fxx, 4, 1, outfile);
	      g_rcb_col++;
	    }
	  } else {
	    while (g_rcb_col <= g_rcb_row) {
	      gzgets(g_rcb_cur_gzfile, tbuf, MAXLINELEN);
	      if (rel_ct_arr[g_rcb_col++] != -1) {
		bufptr = next_item_mult(tbuf, 2);
		sscanf(bufptr, "%g %e", &fxx, &fyy);
		fwrite(&fxx, 4, 1, out_bin_nfile);
		fwrite(&fyy, 4, 1, outfile);
	      }
	    }
	  }
	  g_rcb_new_row++;
	  g_rcb_new_col = 0;
	}
	g_rcb_row++;
	g_rcb_progress += g_rcb_row;
	if (g_rcb_progress >= g_pct * g_rcb_hundredth) {
	  if (g_pct > 10) {
	    putchar('\b');
	  }
	  g_pct = 1 + (g_rcb_progress / g_rcb_hundredth);
	  printf("\b\b%u%%", g_pct - 1);
	  fflush(stdout);
	}
	g_rcb_col = 0;
      }
    }
    putchar('\r');
    sprintf(logbuf, "Pruned relationship matrix written to %s.\n", outname);
    logprintb();
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
  rel_cutoff_batch_ret_INVALID_FORMAT_2:
    putchar('\n');
    logprint("Error: Improperly formatted .grm.gz file.\n");
  rel_cutoff_batch_ret_INVALID_FORMAT_1:
    retval = RET_INVALID_FORMAT;
    break;
  }
 rel_cutoff_batch_ret_1:
  fclose_cond(idfile);
  fclose_cond(g_rcb_in_binfile);
  fclose_cond(g_rcb_in_bin_nfile);
  fclose_cond(outfile);
  fclose_cond(out_bin_nfile);
  gzclose_cond(cur_gzfile);
  gzclose_cond(g_rcb_cur_gzfile);
  gzclose_cond(gz_outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t do_rel_cutoff_f(uint64_t calculation_type, float rel_cutoff, float* rel_ibc, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, char* person_ids, uintptr_t max_person_id_len) {
  int32_t indivs_excluded = 0;
  uint32_t exactly_one_rel_ct = 0;
  unsigned char* wkspace_mark = wkspace_base;
  float* dist_ptr = g_rel_f_dists;
  float* dptr2;
  float* dptr3;
  float* dptr4;
  uint32_t* giptr;
  uint32_t* giptr2;
  // number of too-close relations, -1 if excluded
  int32_t* rel_ct_arr;
  uintptr_t indiv_idx;
  uint32_t uii;
  uintptr_t ulii;
  int32_t kk;
  int32_t mm;
  int32_t retval;
  
  if (wkspace_alloc_i_checked(&rel_ct_arr, g_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_int_zero(rel_ct_arr, g_indiv_ct);
  for (indiv_idx = 1; indiv_idx < g_indiv_ct; indiv_idx++) {
    for (uii = 0; uii < indiv_idx; uii++) {
      if (*dist_ptr++ > rel_cutoff) {
	rel_ct_arr[indiv_idx] += 1;
	rel_ct_arr[uii] += 1;
      }
    }
  }
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    if (rel_ct_arr[indiv_idx] == 1) {
      exactly_one_rel_ct++;
    }
  }
  while (1) {
    kk = 0;
    if (exactly_one_rel_ct) {
      // there is at least one individual with exactly one too-close
      // relation left, find the first one
      while (rel_ct_arr[kk] != 1) {
	kk++;
      }
      // and now find the identity of the other side
      dist_ptr = &(g_rel_f_dists[((intptr_t)kk * (kk - 1)) / 2]);
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
	  dist_ptr = &(g_rel_f_dists[((intptr_t)mm * (mm - 1)) / 2 + kk]);
	} while (*dist_ptr <= rel_cutoff);
	*dist_ptr = 0.0;
      }
      rel_ct_arr[kk] = 0;
      exactly_one_rel_ct--;
      if (rel_ct_arr[mm] == 1) {
        // speed up the easy case
	exactly_one_rel_ct--;
	rel_ct_arr[mm] = -1;
	indivs_excluded++;
	continue;
      }
    } else {
      // find identity of first individual with maximum number of
      // remaining too-close relations
      // kk is highest too-close pair count so far
      mm = -1; // associated individual index
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	if (rel_ct_arr[indiv_idx] > kk) {
	  kk = rel_ct_arr[indiv_idx];
	  mm = indiv_idx;
	}
      }
      // no too-close relations left at all, we're done
      if (mm == -1) {
	break;
      }
    }
    dist_ptr = &(g_rel_f_dists[((intptr_t)mm * (mm - 1)) / 2]);
    for (kk = 0; kk < mm; kk++) {
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[kk]), &exactly_one_rel_ct);
      }
      dist_ptr++;
    }
    for (ulii = mm + 1; ulii < g_indiv_ct; ulii++) {
      dist_ptr = &(g_rel_f_dists[(ulii * (ulii - 1)) / 2 + mm]);
      if (*dist_ptr > rel_cutoff) {
	*dist_ptr = 0.0;
	rel_cut_arr_dec(&(rel_ct_arr[ulii]), &exactly_one_rel_ct);
      }
    }
    rel_ct_arr[mm] = -1;
    indivs_excluded++;
  }
  exclude_multi(indiv_exclude, rel_ct_arr, g_indiv_ct, indiv_exclude_ct_ptr);
  if (indivs_excluded) {
    dist_ptr = g_rel_f_dists; // write
    dptr2 = g_rel_f_dists; // read
    dptr3 = rel_ibc; // write
    dptr4 = rel_ibc; // read
    giptr = g_missing_dbl_excluded; // write
    giptr2 = g_missing_dbl_excluded; // read
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      if (rel_ct_arr[indiv_idx] != -1) {
	if (calculation_type & CALC_IBC) {
	  dptr3[g_indiv_ct] = dptr4[g_indiv_ct];
	  dptr3[g_indiv_ct * 2] = dptr4[g_indiv_ct * 2];
	}
	*dptr3 = *dptr4++;
	dptr3++;
	for (uii = 0; uii < indiv_idx; uii++) {
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
	dptr2 = &(dptr2[indiv_idx]);
	giptr2 = &(giptr2[indiv_idx]);
      }
    }
    g_indiv_ct -= indivs_excluded;
    if (calculation_type & CALC_IBC) {
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	*dptr3++ = *dptr4++;
      }
      dptr4 = &(dptr4[indivs_excluded]);
      for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
	*dptr3++ = *dptr4++;
      }
    }
    giptr = g_indiv_missing_unwt;
    giptr2 = g_indiv_missing_unwt;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct + indivs_excluded; indiv_idx++) {
      if (rel_ct_arr[indiv_idx] != -1) {
	*giptr = *giptr2++;
	giptr++;
      } else {
	giptr2++;
      }
    }
  }
  sprintf(logbuf, "%d %s excluded by --rel-cutoff.\n", indivs_excluded, species_str(indivs_excluded));
  logprintb();
  if (!(calculation_type & (CALC_RELATIONSHIP | CALC_GDISTANCE_MASK))) {
    strcpy(outname_end, ".rel.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      return retval;
    }
    sprintf(logbuf, "Remaining individual IDs written to %s.\n", outname);
    logprintb();
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

static uint32_t g_cr_marker_ct;
static uintptr_t g_cr_indiv1idx;
static uintptr_t g_cr_indiv2idx;
static uintptr_t g_cr_min_indiv;
static uintptr_t g_cr_max_indiv1idx;
static uint64_t g_cr_start_offset;
static uint64_t g_cr_hundredth;
static uint32_t* g_cr_mdeptr;
static double* g_cr_dist_ptr;
static double* g_cr_ibc_ptr;
static float* g_crf_dist_ptr;
static float* g_crf_ibc_ptr;

uint32_t calc_rel_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_cr_dist_ptr++, '\t');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_tri_emitn_ret;
      }
    }
    sptr_cur = double_g_writex(sptr_cur, *g_cr_ibc_ptr++, '\n');
    g_cr_indiv1idx++;
    if ((((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
      g_pct = (((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) / g_cr_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_tri_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t ulii;
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_cr_dist_ptr++, '\t');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq0_emitn_ret;
      }
    }
    if (g_cr_indiv2idx == g_cr_indiv1idx) {
      sptr_cur = double_g_write(sptr_cur, *g_cr_ibc_ptr++);
      g_cr_indiv2idx++;
    }
    if (sptr_cur >= readbuf_end) {
      goto calc_rel_sq0_emitn_ret;
    } else {
      ulii = (1 + (uintptr_t)(readbuf_end - sptr_cur)) / 2;
      if (ulii < (g_indiv_ct - g_cr_indiv2idx)) {
	sptr_cur = memcpya(sptr_cur, g_geno, 2 * ulii);
	g_cr_indiv2idx += ulii;
	goto calc_rel_sq0_emitn_ret;
      }
      ulii = 2 * (g_indiv_ct - g_cr_indiv2idx);
      sptr_cur = memcpya(sptr_cur, g_geno, ulii);
    }
    *sptr_cur++ = '\n';
    g_cr_indiv1idx++;
    if ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_cr_max_indiv1idx - g_cr_min_indiv)) {
      g_pct = ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU) / (g_cr_max_indiv1idx - g_cr_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_sq0_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = double_g_writex(sptr_cur, *g_cr_dist_ptr++, '\t');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq_emitn_ret;
      }
    }
    if (g_cr_indiv2idx == g_cr_indiv1idx) {
      sptr_cur = double_g_write(sptr_cur, *g_cr_ibc_ptr++);
      g_cr_indiv2idx++;
    }
    while (g_cr_indiv2idx < g_indiv_ct) {
      *sptr_cur = '\t';
      sptr_cur = double_g_write(&(sptr_cur[1]), g_rel_dists[((g_cr_indiv2idx * (g_cr_indiv2idx - 1)) / 2) + g_cr_indiv1idx]);
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    g_cr_indiv1idx++;
    if ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_cr_max_indiv1idx - g_cr_min_indiv)) {
      g_pct = ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU) / (g_cr_max_indiv1idx - g_cr_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_sq_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_grm_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  char wbuf[16];
  char* wbuf_end;
  uint32_t wbuf_len;
  uint32_t uii;
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    uii = g_cr_marker_ct - g_indiv_missing_unwt[g_cr_indiv1idx];
    wbuf_end = uint32_writex(wbuf, g_cr_indiv1idx + 1, '\t');
    wbuf_len = (uintptr_t)(wbuf_end - wbuf);
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = double_e_writex(uint32_writex(uint32_writex(memcpya(sptr_cur, wbuf, wbuf_len), g_cr_indiv2idx + 1, '\t'), (uii - g_indiv_missing_unwt[g_cr_indiv2idx]) + (*g_cr_mdeptr++), '\t'), *g_cr_dist_ptr++, '\n');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_grm_emitn_ret;
      }
    }
    sptr_cur = double_e_writex(uint32_writex(uint32_writex(memcpya(sptr_cur, wbuf, wbuf_len), ++g_cr_indiv1idx, '\t'), uii, '\t'), *g_cr_ibc_ptr++, '\n');
    if ((((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
      g_pct = (((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) / g_cr_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_grm_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t calc_rel(pthread_t* threads, int32_t parallel_idx, int32_t parallel_tot, uint64_t calculation_type, int32_t rel_calc_type, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* person_ids, uintptr_t max_person_id_len, int32_t var_std, int32_t ibc_type, double rel_cutoff, double* set_allele_freqs, double** rel_ibc_ptr, Chrom_info* chrom_info_ptr) {
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  int32_t retval = 0;
  int64_t llxx = 0;
  double* dist_ptr = NULL;
  double* dptr3 = NULL;
  double* dptr4 = NULL;
  uint32_t chrom_fo_idx = 0;
  double* dptr2;
  double set_allele_freq_buf[MULTIPLEX_DIST];
  char wbuf[96];
  char* wptr;
  uint32_t cur_markers_loaded;
  uint32_t win_marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  double* rel_ibc;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t rel_shape;
  uint32_t min_indiv;
  uint32_t max_parallel_indiv;
  unsigned char* wkspace_mark;
  unsigned char* gptr;
  unsigned char* gptr2;
  uint32_t* giptr;
  uint32_t* giptr2;
  uintptr_t* glptr2;
  if (wkspace_alloc_ui_checked(&g_indiv_missing_unwt, g_indiv_ct * sizeof(int32_t))) {
    goto calc_rel_ret_NOMEM;
  }
  fill_int_zero((int32_t*)g_indiv_missing_unwt, g_indiv_ct);
;
  triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_idx, parallel_tot, 1, 1);
  if (relationship_req(calculation_type)) {
    llxx = g_thread_start[g_thread_ct];
    llxx = ((llxx * (llxx - 1)) - (int64_t)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
    if (!(calculation_type & CALC_UNRELATED_HERITABILITY)) {
      // if the memory isn't needed for CALC_UNRELATED_HERITABILITY,
      // positioning the missingness matrix here will let us avoid
      // recalculating it if --distance-matrix or --matrix is requested
      if (wkspace_alloc_ui_checked(&g_missing_dbl_excluded, llxx * sizeof(int32_t))) {
	goto calc_rel_ret_NOMEM;
      }
      fill_int_zero((int32_t*)g_missing_dbl_excluded, llxx);
    }
    if (wkspace_alloc_d_checked(&g_rel_dists, llxx * sizeof(double))) {
      goto calc_rel_ret_NOMEM;
    }
    fill_double_zero(g_rel_dists, llxx);
  }
  if (calculation_type & CALC_IBC) {
    uii = g_indiv_ct * 3;
  } else {
    uii = g_indiv_ct;
  }
  if (wkspace_alloc_d_checked(rel_ibc_ptr, uii * sizeof(double))) {
    goto calc_rel_ret_NOMEM;
  }
  rel_ibc = *rel_ibc_ptr;
  fill_double_zero(rel_ibc, uii);
  wkspace_mark = wkspace_base;
  if (relationship_req(calculation_type) && (!g_missing_dbl_excluded)) {
    if (wkspace_alloc_ui_checked(&g_missing_dbl_excluded, llxx * sizeof(int32_t))) {
      goto calc_rel_ret_NOMEM;
    }
    fill_int_zero((int32_t*)g_missing_dbl_excluded, llxx);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_rel_ret_READ_FAIL;
  }
  if (wkspace_alloc_uc_checked(&g_geno, g_indiv_ct * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_mmasks, g_indiv_ct * sizeof(intptr_t)) ||
      wkspace_alloc_uc_checked(&gptr, MULTIPLEX_REL * unfiltered_indiv_ct4) ||
      wkspace_alloc_ul_checked(&g_masks, g_indiv_ct * sizeof(intptr_t))) {
    goto calc_rel_ret_NOMEM;
  }

  // Exclude markers on non-autosomal chromosomes for now.
  uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1);
  if (uii) {
    if (uii == marker_ct) {
      logprint("Error: No autosomal markers for relationship matrix calculation.\n");
      retval = RET_INVALID_CMDLINE;
      goto calc_rel_ret_1;
    }
    sprintf(logbuf, "Excluding %u marker%s on non-autosomes from relationship matrix calc.\n", uii, (uii == 1)? "" : "s");
    logprintb();
    marker_ct -= uii;
  }

  // See comments at the beginning of this file, and those in the main
  // CALC_DISTANCE loop.  The main difference between this calculation and
  // the (nonzero exponent) distance calculation is that we have to pad
  // each marker to 3 bits and use + instead of XOR to distinguish the
  // cases.
  do {
    retval = block_load_autosomal(bedfile, bed_offset, marker_exclude, marker_ct, MULTIPLEX_REL, unfiltered_indiv_ct4, chrom_info_ptr, set_allele_freqs, NULL, gptr, &chrom_fo_idx, &marker_uidx, &marker_idx, &cur_markers_loaded, set_allele_freq_buf, NULL, NULL);
    if (retval) {
      goto calc_rel_ret_1;
    }
    if (cur_markers_loaded < MULTIPLEX_REL) {
      memset(&(gptr[cur_markers_loaded * unfiltered_indiv_ct4]), 0, (MULTIPLEX_REL - cur_markers_loaded) * unfiltered_indiv_ct4);
      fill_double_zero(&(set_allele_freq_buf[cur_markers_loaded]), MULTIPLEX_REL - cur_markers_loaded);
    }
    fill_ulong_zero(g_mmasks, g_indiv_ct);

    for (win_marker_idx = 0; win_marker_idx < cur_markers_loaded; win_marker_idx += MULTIPLEX_REL / 3) {
      fill_ulong_zero(g_masks, g_indiv_ct);
      indiv_idx = 0;
      glptr2 = (uintptr_t*)g_geno;
      for (indiv_uidx = 0; indiv_idx < g_indiv_ct; indiv_uidx++) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	ulii = 0;
	gptr2 = &(gptr[indiv_uidx / 4 + win_marker_idx * unfiltered_indiv_ct4]);
	uii = (indiv_uidx % 4) * 2;
	umm = 0;
	unn = 0;
	for (ukk = 0; ukk < (BITCT / 16); ukk++) {
	  for (ujj = 0; ujj < 5; ujj++) {
	    uljj = (gptr2[umm * unfiltered_indiv_ct4] >> uii) & 3;
	    if (uljj == 1) {
	      g_masks[indiv_idx] |= (7 * ONELU) << unn;
	      g_mmasks[indiv_idx] |= ONELU << (win_marker_idx + umm);
	      g_indiv_missing_unwt[indiv_idx] += 1;
	    }
	    ulii |= uljj << unn;
	    umm++;
	    unn += 3;
	  }
	  unn++;
	}
	*glptr2++ = ulii;
	indiv_idx++;
      }
      if (calculation_type & CALC_IBC) {
	for (uii = 0; uii < 3; uii++) {
	  update_rel_ibc(&(rel_ibc[uii * g_indiv_ct]), (uintptr_t*)g_geno, &(set_allele_freq_buf[win_marker_idx]), uii, g_indiv_ct);
	}
      } else {
	update_rel_ibc(rel_ibc, (uintptr_t*)g_geno, &(set_allele_freq_buf[win_marker_idx]), ibc_type, g_indiv_ct);
      }
      if (relationship_req(calculation_type)) {
	fill_weights_r(g_weights, &(set_allele_freq_buf[win_marker_idx]), var_std);
	if (spawn_threads(threads, &calc_rel_thread, g_thread_ct)) {
	  goto calc_rel_ret_THREAD_CREATE_FAIL;
	}
	incr_dists_r(g_rel_dists, (uintptr_t*)g_geno, g_masks, 0, g_weights);
	join_threads(threads, g_thread_ct);
      }
    }
    if (relationship_req(calculation_type)) {
      if (spawn_threads(threads, &calc_missing_thread, g_thread_ct)) {
	goto calc_rel_ret_THREAD_CREATE_FAIL;
      }
      incr_dists_rm(g_missing_dbl_excluded, 0, g_thread_start);
      join_threads(threads, g_thread_ct);
    }
    printf("\r%" PRIuPTR " markers complete.", marker_idx);
    fflush(stdout);
  } while (marker_idx < marker_ct);
  if (relationship_req(calculation_type)) {
    putchar('\r');
    logprint("Relationship matrix calculation complete.\n");
    dist_ptr = g_rel_dists;
  } else {
    putchar('\n');
  }
  dptr2 = rel_ibc;
  if (calculation_type & CALC_IBC) {
    dptr3 = &(rel_ibc[g_indiv_ct]);
    dptr4 = &(rel_ibc[g_indiv_ct * 2]);
  }
  giptr2 = g_missing_dbl_excluded;
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    uii = marker_ct - g_indiv_missing_unwt[indiv_idx];
    if ((indiv_idx >= g_thread_start[0]) && (indiv_idx < g_thread_start[g_thread_ct])) {
      if (relationship_req(calculation_type)) {
	giptr = g_indiv_missing_unwt;
	for (ujj = 0; ujj < indiv_idx; ujj++) {
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
    retval = do_rel_cutoff(calculation_type, rel_cutoff, rel_ibc, indiv_exclude, indiv_exclude_ct_ptr, outname, outname_end, unfiltered_indiv_ct, person_ids, max_person_id_len);
    if (retval) {
      goto calc_rel_ret_1;
    }
  }

  if (calculation_type & CALC_IBC) {
    strcpy(outname_end, ".ibc");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_rel_ret_OPEN_FAIL;
    }
    dptr2 = rel_ibc;
    dptr3 = &(rel_ibc[g_indiv_ct]);
    dptr4 = &(rel_ibc[g_indiv_ct * 2]);
    if (fputs_checked("FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3\n", outfile)) {
      goto calc_rel_ret_WRITE_FAIL;
    }
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      wptr = double_g_writex(double_g_writex(double_g_writex(uint32_writex(uint32_writex(uint32_writex(wbuf, indiv_idx + 1, '\t'), indiv_idx + 1, '\t'), marker_ct - g_indiv_missing_unwt[indiv_idx], '\t'), *dptr3++ - 1.0, '\t'), *dptr4++ - 1.0, '\t'), *dptr2++ - 1.0, '\n');
      if (fwrite_checked(wbuf, wptr - wbuf, outfile)) {
	goto calc_rel_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_rel_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "%s written.\n", outname);
    logprintb();
  }
  if (calculation_type & CALC_RELATIONSHIP) {
    rel_shape = rel_calc_type & REL_CALC_SHAPEMASK;
    if (parallel_tot == 1) {
      // nasty --rel-cutoff bug
      max_parallel_indiv = g_indiv_ct;
    } else {
      // can't run --rel-cutoff with --parallel, so this is safe
      max_parallel_indiv = g_thread_start[g_thread_ct];
    }
    min_indiv = g_thread_start[0];
    if (min_indiv == 1) {
      min_indiv = 0;
    }
    if (calculation_type & CALC_IBC) {
      dptr2 = &(rel_ibc[ibc_type * g_indiv_ct + min_indiv]);
    } else {
      dptr2 = &(rel_ibc[min_indiv]);
    }
    g_cr_marker_ct = marker_ct;
    g_pct = 1;
    g_cr_start_offset = ((uint64_t)min_indiv * (min_indiv - 1)) / 2;
    g_cr_hundredth = 1 + (((((uint64_t)max_parallel_indiv * (max_parallel_indiv + 1)) / 2) - g_cr_start_offset) / 100);
    if (rel_calc_type & REL_CALC_BIN) {
      if (rel_shape == REL_CALC_SQ0) {
	fill_double_zero((double*)g_geno, g_indiv_ct - 1);
      }
      strcpy(outname_end, ".rel.bin");
      if (parallel_tot > 1) {
	sprintf(&(outname_end[8]), ".%u", parallel_idx + 1);
      }
      if (fopen_checked(&outfile, outname, "wb")) {
	goto calc_rel_ret_OPEN_FAIL;
      }
      for (indiv_idx = min_indiv; indiv_idx < max_parallel_indiv; indiv_idx++) {
	if (fwrite_checkedz(&(g_rel_dists[((int64_t)indiv_idx * (indiv_idx - 1)) / 2 - g_cr_start_offset]), indiv_idx * sizeof(double), outfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	if (fwrite_checked(dptr2++, sizeof(double), outfile)) {
	  goto calc_rel_ret_WRITE_FAIL;
	}
	if (rel_shape == REL_CALC_TRI) {
	  if ((((uint64_t)indiv_idx + 1) * (indiv_idx + 2) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
	    g_pct = (((uint64_t)indiv_idx + 1) * (indiv_idx + 2) / 2 - g_cr_start_offset) / g_cr_hundredth;
	    printf("\rWriting... %u%%", g_pct++);
	    fflush(stdout);
	  }
	} else {
	  if (rel_shape == REL_CALC_SQ0) {
	    if (fwrite_checkedz(g_geno, (g_indiv_ct - indiv_idx - 1) * sizeof(double), outfile)) {
	      goto calc_rel_ret_WRITE_FAIL;
	    }
	  } else {
	    for (uii = indiv_idx + 1; uii < g_indiv_ct; uii++) {
	      if (fwrite_checked(&(g_rel_dists[((uintptr_t)uii * (uii - 1) / 2) + indiv_idx - g_cr_start_offset]), sizeof(double), outfile)) {
		goto calc_rel_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((indiv_idx + 1 - min_indiv) * 100 >= g_pct * (max_parallel_indiv - min_indiv)) {
	    g_pct = ((indiv_idx + 1 - min_indiv) * 100) / (max_parallel_indiv - min_indiv);
	    printf("\rWriting... %u%%", g_pct++);
	    fflush(stdout);
	  }
	}
      }
      if (fclose_null(&outfile)) {
	goto calc_rel_ret_WRITE_FAIL;
      }
    } else {
      g_cr_indiv1idx = min_indiv;
      g_cr_indiv2idx = 0;
      g_cr_max_indiv1idx = max_parallel_indiv;
      g_cr_mdeptr = g_missing_dbl_excluded;
      g_cr_dist_ptr = g_rel_dists;
      g_cr_ibc_ptr = dptr2;
      if (rel_calc_type & REL_CALC_GRM) {
	if (rel_calc_type & REL_CALC_GZ) {
	  if (parallel_tot > 1) {
	    sprintf(outname_end, ".grm.%u.gz", parallel_idx + 1);
	  } else {
	    strcpy(outname_end, ".grm.gz");
	  }
	  parallel_compress(outname, calc_rel_grm_emitn);
	} else {
	  strcpy(outname_end, ".grm");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[4]), ".%u", parallel_idx + 1);
	  }
	  retval = write_uncompressed(outname, calc_rel_grm_emitn);
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
	    parallel_compress(outname, calc_rel_tri_emitn);
	  } else {
	    retval = write_uncompressed(outname, calc_rel_tri_emitn);
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
	  uii = (g_indiv_ct * 2 + sizeof(intptr_t) - 4) / sizeof(intptr_t);
	  glptr2 = (uintptr_t*)g_geno;
	  for (ujj = 0; ujj < uii; ujj++) {
	    *glptr2++ = ulii;
	  }
	  g_cr_min_indiv = min_indiv;
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, calc_rel_sq0_emitn);
	  } else {
	    retval = write_uncompressed(outname, calc_rel_sq0_emitn);
	    if (retval) {
	      goto calc_rel_ret_1;
	    }
	  }
	} else {
	  g_cr_min_indiv = min_indiv;
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, calc_rel_sq_emitn);
	  } else {
	    retval = write_uncompressed(outname, calc_rel_sq_emitn);
	    if (retval) {
	      goto calc_rel_ret_1;
	    }
	  }
	}
      }
    }
    putchar('\r');
    sprintf(logbuf, "Relationship matrix written to %s.\n", outname);
    logprintb();
    if (!parallel_idx) {
      strcpy(&(outname_end[4]), ".id");
      retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
      if (retval) {
	goto calc_rel_ret_1;
      }
    }
  }
  wkspace_reset(wkspace_mark);
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
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_rel_ret_1:
  fclose_cond(outfile);
  gzclose_cond(gz_outfile);
  return retval;
}

uint32_t calc_rel_f_tri_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = float_g_writex(sptr_cur, *g_crf_dist_ptr++, '\t');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_f_tri_emitn_ret;
      }
    }
    sptr_cur = float_g_writex(sptr_cur, *g_crf_ibc_ptr++, '\n');
    g_cr_indiv1idx++;
    if ((((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
      g_pct = (((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) / g_cr_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_f_tri_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_f_sq0_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  uintptr_t ulii;
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = float_g_writex(sptr_cur, *g_crf_dist_ptr++, '\t');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq0_emitn_ret;
      }
    }
    if (g_cr_indiv2idx == g_cr_indiv1idx) {
      sptr_cur = float_g_write(sptr_cur, *g_crf_ibc_ptr++);
      g_cr_indiv2idx++;
    }
    if (sptr_cur >= readbuf_end) {
      goto calc_rel_sq0_emitn_ret;
    } else {
      ulii = (1 + (uintptr_t)(readbuf_end - sptr_cur)) / 2;
      if (ulii < (g_indiv_ct - g_cr_indiv2idx)) {
	sptr_cur = memcpya(sptr_cur, g_geno, 2 * ulii);
	g_cr_indiv2idx += ulii;
	goto calc_rel_sq0_emitn_ret;
      }
      ulii = 2 * (g_indiv_ct - g_cr_indiv2idx);
      sptr_cur = memcpya(sptr_cur, g_geno, ulii);
    }
    *sptr_cur++ = '\n';
    g_cr_indiv1idx++;
    if ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_cr_max_indiv1idx - g_cr_min_indiv)) {
      g_pct = ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU) / (g_cr_max_indiv1idx - g_cr_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_sq0_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_f_sq_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = float_g_writex(sptr_cur, *g_crf_dist_ptr++, '\t');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq_emitn_ret;
      }
    }
    if (g_cr_indiv2idx == g_cr_indiv1idx) {
      sptr_cur = float_g_write(sptr_cur, *g_crf_ibc_ptr++);
      g_cr_indiv2idx++;
    }
    while (g_cr_indiv2idx < g_indiv_ct) {
      *sptr_cur++ = '\t';
      sptr_cur = float_g_write(sptr_cur, g_rel_f_dists[((g_cr_indiv2idx * (g_cr_indiv2idx - 1)) / 2) + g_cr_indiv1idx]);
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	goto calc_rel_sq_emitn_ret;
      }
    }
    *sptr_cur++ = '\n';
    g_cr_indiv1idx++;
    if ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU >= ((uint64_t)g_pct) * (g_cr_max_indiv1idx - g_cr_min_indiv)) {
      g_pct = ((g_cr_indiv1idx - g_cr_min_indiv) * 100LLU) / (g_cr_max_indiv1idx - g_cr_min_indiv);
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
 calc_rel_sq_emitn_ret:
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

uint32_t calc_rel_f_grm_emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  char* sptr_cur = (char*)(&(readbuf[overflow_ct]));
  char* readbuf_end = (char*)(&(readbuf[PIGZ_BLOCK_SIZE]));
  char wbuf[16];
  char* wbuf_end;
  uint32_t wbuf_len;
  uint32_t uii;
  while (g_cr_indiv1idx < g_cr_max_indiv1idx) {
    uii = g_cr_marker_ct - g_indiv_missing_unwt[g_cr_indiv1idx];
    wbuf_end = uint32_writex(wbuf, g_cr_indiv1idx + 1, '\t');
    wbuf_len = (uintptr_t)(wbuf_end - wbuf);
    while (g_cr_indiv2idx < g_cr_indiv1idx) {
      sptr_cur = float_e_writex(uint32_writex(uint32_writex(memcpya(sptr_cur, wbuf, wbuf_len), g_cr_indiv2idx + 1, '\t'), (uii - g_indiv_missing_unwt[g_cr_indiv2idx]) + (*g_cr_mdeptr++), '\t'), *g_crf_dist_ptr++, '\n');
      g_cr_indiv2idx++;
      if (sptr_cur >= readbuf_end) {
	return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
      }
    }
    sptr_cur = float_e_writex(uint32_writex(uint32_writex(memcpya(sptr_cur, wbuf, wbuf_len), ++g_cr_indiv1idx, '\t'), uii, '\t'), *g_crf_ibc_ptr++, '\n');
    if ((((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
      g_pct = (((uint64_t)g_cr_indiv1idx) * (g_cr_indiv1idx + 1) / 2 - g_cr_start_offset) / g_cr_hundredth;
      printf("\rWriting... %u%%", g_pct++);
      fflush(stdout);
    }
    g_cr_indiv2idx = 0;
  }
  return (uintptr_t)(((unsigned char*)sptr_cur) - readbuf);
}

int32_t calc_rel_f(pthread_t* threads, int32_t parallel_idx, int32_t parallel_tot, uint64_t calculation_type, int32_t rel_calc_type, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* person_ids, uintptr_t max_person_id_len, int32_t var_std, int32_t ibc_type, float rel_cutoff, double* set_allele_freqs, Chrom_info* chrom_info_ptr) {
  // N.B. ACTA may currently outperform this when compiled with ICC and run on
  // a heavily multicore 64-bit Linux system.  If ACTA ever gets to the point
  // where it wins when compiled with gcc, on both 32- and 64-bit systems with
  // as few as a single core, and on OS X/Windows as well as Linux, then it's
  // time to replace this with the ACTA implementation.
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  FILE* outfile = NULL;
  FILE* out_bin_nfile = NULL;
  gzFile gz_outfile = NULL;
  int32_t retval = 0;
  uint64_t ullxx = 0;
  float* dist_ptr = NULL;
  float* dptr3 = NULL;
  float* dptr4 = NULL;
  uint32_t chrom_fo_idx = 0;
  float* dptr2;
  float set_allele_freq_buf[MULTIPLEX_DIST];
  char wbuf[96];
  char* wptr;
  uint32_t cur_markers_loaded;
  uint32_t win_marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  float* rel_ibc;
  uintptr_t ulii;
  uintptr_t uljj;
  float fxx;
  uint32_t* mdeptr;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t rel_shape;
  uint32_t min_indiv;
  uint32_t max_parallel_indiv;
  unsigned char* wkspace_mark;
  unsigned char* gptr;
  unsigned char* gptr2;
  uint32_t* giptr;
  uint32_t* giptr2;
  uintptr_t* glptr2;
  if (wkspace_alloc_ui_checked(&g_indiv_missing_unwt, g_indiv_ct * sizeof(int32_t))) {
    goto calc_rel_f_ret_NOMEM;
  }
  fill_int_zero((int32_t*)g_indiv_missing_unwt, g_indiv_ct);
;
  triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_idx, parallel_tot, 1, 1);
  if (relationship_req(calculation_type)) {
    ullxx = g_thread_start[g_thread_ct];
    ullxx = ((ullxx * (ullxx - 1)) - (int64_t)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
    if (wkspace_alloc_ui_checked(&g_missing_dbl_excluded, ullxx * sizeof(int32_t)) ||
        wkspace_alloc_f_checked(&g_rel_f_dists, ullxx * sizeof(float))) {
      goto calc_rel_f_ret_NOMEM;
    }
    fill_int_zero((int32_t*)g_missing_dbl_excluded, ullxx);
    fill_float_zero(g_rel_f_dists, ullxx);
  }
  if (calculation_type & CALC_IBC) {
    uii = g_indiv_ct * 3;
  } else {
    uii = g_indiv_ct;
  }
  if (wkspace_alloc_f_checked(&rel_ibc, uii * sizeof(float))) {
    goto calc_rel_f_ret_NOMEM;
  }
  fill_float_zero(rel_ibc, uii);
  wkspace_mark = wkspace_base;

  if (relationship_req(calculation_type) && (!g_missing_dbl_excluded)) {
    if (wkspace_alloc_ui_checked(&g_missing_dbl_excluded, ullxx * sizeof(int32_t))) {
      goto calc_rel_f_ret_NOMEM;
    }
    fill_int_zero((int32_t*)g_missing_dbl_excluded, ullxx);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_rel_f_ret_READ_FAIL;
  }
  if (wkspace_alloc_uc_checked(&g_geno, g_indiv_ct * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_mmasks, g_indiv_ct * sizeof(intptr_t)) ||
      wkspace_alloc_uc_checked(&gptr, MULTIPLEX_REL * unfiltered_indiv_ct4) ||
      wkspace_alloc_ul_checked(&g_masks, g_indiv_ct * sizeof(intptr_t))) {
    goto calc_rel_f_ret_NOMEM;
  }

  uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1);
  if (uii) {
    if (uii == marker_ct) {
      logprint("Error: No autosomal markers for relationship matrix calculation.\n");
      retval = RET_INVALID_CMDLINE;
      goto calc_rel_f_ret_1;
    }
    sprintf(logbuf, "Excluding %u marker%s on non-autosomes from relationship matrix calc.\n", uii, (uii == 1)? "" : "s");
    logprintb();
    marker_ct -= uii;
  }

  do {
    retval = block_load_autosomal(bedfile, bed_offset, marker_exclude, marker_ct, MULTIPLEX_REL, unfiltered_indiv_ct4, chrom_info_ptr, set_allele_freqs, NULL, gptr, &chrom_fo_idx, &marker_uidx, &marker_idx, &cur_markers_loaded, NULL, set_allele_freq_buf, NULL);
    if (retval) {
      goto calc_rel_f_ret_1;
    }
    if (cur_markers_loaded < MULTIPLEX_REL) {
      memset(&(gptr[cur_markers_loaded * unfiltered_indiv_ct4]), 0, (MULTIPLEX_REL - cur_markers_loaded) * unfiltered_indiv_ct4);
      fill_float_zero(&(set_allele_freq_buf[cur_markers_loaded]), MULTIPLEX_REL - cur_markers_loaded);
    }
    fill_ulong_zero(g_mmasks, g_indiv_ct);

    for (win_marker_idx = 0; win_marker_idx < cur_markers_loaded; win_marker_idx += MULTIPLEX_REL / 3) {
      fill_ulong_zero(g_masks, g_indiv_ct);
      indiv_idx = 0;
      glptr2 = (uintptr_t*)g_geno;
      for (indiv_uidx = 0; indiv_idx < g_indiv_ct; indiv_uidx++) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	ulii = 0;
	gptr2 = &(gptr[indiv_uidx / 4 + win_marker_idx * unfiltered_indiv_ct4]);
	uii = (indiv_uidx % 4) * 2;
	umm = 0;
	unn = 0;
	for (ukk = 0; ukk < (BITCT / 16); ukk++) {
	  for (ujj = 0; ujj < 5; ujj++) {
	    uljj = (gptr2[umm * unfiltered_indiv_ct4] >> uii) & 3;
	    if (uljj == 1) {
	      g_masks[indiv_idx] |= (7 * ONELU) << unn;
	      g_mmasks[indiv_idx] |= ONELU << (win_marker_idx + umm);
	      g_indiv_missing_unwt[indiv_idx] += 1;
	    }
	    ulii |= uljj << unn;
	    umm++;
	    unn += 3;
	  }
	  unn++;
	}
	*glptr2++ = ulii;
	indiv_idx++;
      }
      if (calculation_type & CALC_IBC) {
	for (uii = 0; uii < 3; uii++) {
	  update_rel_f_ibc(&(rel_ibc[uii * g_indiv_ct]), (uintptr_t*)g_geno, &(set_allele_freq_buf[win_marker_idx]), uii, g_indiv_ct);
	}
      } else {
	update_rel_f_ibc(rel_ibc, (uintptr_t*)g_geno, &(set_allele_freq_buf[win_marker_idx]), ibc_type, g_indiv_ct);
      }
      if (relationship_req(calculation_type)) {
	fill_weights_r_f(g_weights_f, &(set_allele_freq_buf[win_marker_idx]), var_std);
	if (spawn_threads(threads, &calc_rel_f_thread, g_thread_ct)) {
	  goto calc_rel_f_ret_THREAD_CREATE_FAIL;
	}
	incr_dists_r_f(g_rel_f_dists, (uintptr_t*)g_geno, g_masks, 0, g_weights_f);
	join_threads(threads, g_thread_ct);
      }
    }
    if (relationship_req(calculation_type)) {
      if (spawn_threads(threads, &calc_missing_thread, g_thread_ct)) {
	goto calc_rel_f_ret_THREAD_CREATE_FAIL;
      }
      incr_dists_rm(g_missing_dbl_excluded, 0, g_thread_start);
      join_threads(threads, g_thread_ct);
    }
    printf("\r%" PRIuPTR " markers complete.", marker_idx);
    fflush(stdout);
  } while (marker_idx < marker_ct);
  if (relationship_req(calculation_type)) {
    putchar('\r');
    logprint("Single-precision relationship matrix calculation complete.\n");
    dist_ptr = g_rel_f_dists;
  } else {
    putchar('\n');
  }
  dptr2 = rel_ibc;
  if (calculation_type & CALC_IBC) {
    dptr3 = &(rel_ibc[g_indiv_ct]);
    dptr4 = &(rel_ibc[g_indiv_ct * 2]);
  }
  giptr2 = g_missing_dbl_excluded;
  for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
    uii = marker_ct - g_indiv_missing_unwt[indiv_idx];
    if ((indiv_idx >= g_thread_start[0]) && (indiv_idx < g_thread_start[g_thread_ct])) {
      if (relationship_req(calculation_type)) {
	giptr = g_indiv_missing_unwt;
	for (ujj = 0; ujj < indiv_idx; ujj++) {
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
    retval = do_rel_cutoff_f(calculation_type, rel_cutoff, rel_ibc, indiv_exclude, indiv_exclude_ct_ptr, outname, outname_end, unfiltered_indiv_ct, person_ids, max_person_id_len);
    if (retval) {
      goto calc_rel_f_ret_1;
    }
  }

  if (calculation_type & CALC_IBC) {
    strcpy(outname_end, ".ibc");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_rel_f_ret_OPEN_FAIL;
    }
    dptr2 = rel_ibc;
    dptr3 = &(rel_ibc[g_indiv_ct]);
    dptr4 = &(rel_ibc[g_indiv_ct * 2]);
    if (fputs_checked("FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3\n", outfile)) {
      goto calc_rel_f_ret_WRITE_FAIL;
    }
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      wptr = float_g_writex(float_g_writex(float_g_writex(uint32_writex(uint32_writex(uint32_writex(wbuf, indiv_idx + 1, '\t'), indiv_idx + 1, '\t'), marker_ct - g_indiv_missing_unwt[indiv_idx], '\t'), *dptr3++ - 1.0, '\t'), *dptr4++ - 1.0, '\t'), *dptr2++ - 1.0, '\n');
      if (fwrite_checked(wbuf, wptr - wbuf, outfile)) {
	goto calc_rel_f_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_rel_f_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "%s written.\n", outname);
    logprintb();
  }
  if (calculation_type & CALC_RELATIONSHIP) {
    rel_shape = rel_calc_type & REL_CALC_SHAPEMASK;
    if (parallel_tot == 1) {
      max_parallel_indiv = g_indiv_ct;
    } else {
      max_parallel_indiv = g_thread_start[g_thread_ct];
    }
    min_indiv = g_thread_start[0];
    if (min_indiv == 1) {
      min_indiv = 0;
    }
    if (calculation_type & CALC_IBC) {
      dptr2 = &(rel_ibc[ibc_type * g_indiv_ct + min_indiv]);
    } else {
      dptr2 = &(rel_ibc[min_indiv]);
    }
    g_cr_marker_ct = marker_ct;
    g_pct = 1;
    g_cr_start_offset = ((uint64_t)min_indiv * (min_indiv - 1)) / 2;
    g_cr_hundredth = 1 + (((((uint64_t)max_parallel_indiv * (max_parallel_indiv + 1)) / 2) - g_cr_start_offset) / 100);

    if (rel_calc_type & REL_CALC_BIN) {
      if (rel_shape == REL_CALC_SQ0) {
	fill_float_zero((float*)g_geno, g_indiv_ct - 1);
      }
      strcpy(outname_end, ".rel.bin");
      if (parallel_tot > 1) {
	sprintf(&(outname_end[8]), ".%u", parallel_idx + 1);
      }
      if (fopen_checked(&outfile, outname, "wb")) {
	goto calc_rel_f_ret_OPEN_FAIL;
      }
      for (indiv_idx = min_indiv; indiv_idx < max_parallel_indiv; indiv_idx++) {
	if (fwrite_checkedz(&(g_rel_f_dists[((int64_t)indiv_idx * (indiv_idx - 1)) / 2 - g_cr_start_offset]), indiv_idx * sizeof(float), outfile)) {
	  goto calc_rel_f_ret_WRITE_FAIL;
	}
	if (fwrite_checked(dptr2++, sizeof(float), outfile)) {
	  goto calc_rel_f_ret_WRITE_FAIL;
	}
	if (rel_shape == REL_CALC_TRI) {
	  if ((((uint64_t)indiv_idx + 1) * (indiv_idx + 2) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
	    g_pct = (((uint64_t)indiv_idx + 1) * (indiv_idx + 2) / 2 - g_cr_start_offset) / g_cr_hundredth;
	    printf("\rWriting... %u%%", g_pct++);
	    fflush(stdout);
	  }
	} else {
	  if (rel_shape == REL_CALC_SQ0) {
	    if (fwrite_checkedz(g_geno, (g_indiv_ct - indiv_idx - 1) * sizeof(float), outfile)) {
	      goto calc_rel_f_ret_WRITE_FAIL;
	    }
	  } else {
	    for (uii = indiv_idx + 1; uii < g_indiv_ct; uii++) {
	      if (fwrite_checked(&(g_rel_f_dists[((uintptr_t)uii * (uii - 1) / 2) + indiv_idx - g_cr_start_offset]), sizeof(float), outfile)) {
		goto calc_rel_f_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if ((indiv_idx + 1 - min_indiv) * 100 >= g_pct * (max_parallel_indiv - min_indiv)) {
	    g_pct = ((indiv_idx + 1 - min_indiv) * 100) / (max_parallel_indiv - min_indiv);
	    printf("\rWriting... %u%%", g_pct++);
	    fflush(stdout);
	  }
	}
      }
      if (fclose_null(&outfile)) {
	goto calc_rel_f_ret_WRITE_FAIL;
      }
    } else if (rel_calc_type & REL_CALC_GRM_BIN) {
      memcpy(outname_end, ".grm.N.bin", 11);
      if (parallel_tot > 1) {
	outname[10] = '.';
	uint32_writex(&(outname[11]), parallel_idx + 1, '\0');
      }
      if (fopen_checked(&out_bin_nfile, outname, "wb")) {
	goto calc_rel_f_ret_OPEN_FAIL;
      }
      memcpy(outname_end, ".grm.bin", 9);
      if (parallel_tot > 1) {
	outname[8] = '.';
	uint32_writex(&(outname[9]), parallel_idx + 1, '\0');
      }
      if (fopen_checked(&outfile, outname, "wb")) {
	goto calc_rel_f_ret_OPEN_FAIL;
      }
      mdeptr = g_missing_dbl_excluded;
      for (indiv_idx = min_indiv; indiv_idx < max_parallel_indiv; indiv_idx++) {
	if (fwrite_checkedz(&(g_rel_f_dists[((int64_t)indiv_idx * (indiv_idx - 1)) / 2 - g_cr_start_offset]), indiv_idx * sizeof(float), outfile)) {
	  goto calc_rel_f_ret_WRITE_FAIL;
	}
	if (fwrite_checked(dptr2++, sizeof(float), outfile)) {
	  goto calc_rel_f_ret_WRITE_FAIL;
	}
	uii = marker_ct - g_indiv_missing_unwt[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  fxx = (float)((int32_t)(uii - g_indiv_missing_unwt[ujj] + (*mdeptr++)));
	  fwrite(&fxx, 4, 1, out_bin_nfile);
	}
	fxx = (float)((int32_t)uii);
	if (fwrite_checked(&fxx, sizeof(float), out_bin_nfile)) {
	  goto calc_rel_f_ret_WRITE_FAIL;
	}
	if ((((uint64_t)indiv_idx + 1) * (indiv_idx + 2) / 2 - g_cr_start_offset) >= g_cr_hundredth * g_pct) {
	  g_pct = (((uint64_t)indiv_idx + 1) * (indiv_idx + 2) / 2 - g_cr_start_offset) / g_cr_hundredth;
	  printf("\rWriting... %u%%", g_pct++);
	  fflush(stdout);
	}
      }
    } else {
      g_cr_indiv1idx = min_indiv;
      g_cr_indiv2idx = 0;
      g_cr_max_indiv1idx = max_parallel_indiv;
      g_cr_mdeptr = g_missing_dbl_excluded;
      g_crf_dist_ptr = g_rel_f_dists;
      g_crf_ibc_ptr = dptr2;
      if (rel_calc_type & REL_CALC_GRM) {
	if (rel_calc_type & REL_CALC_GZ) {
	  if (parallel_tot > 1) {
	    sprintf(outname_end, ".grm.%u.gz", parallel_idx + 1);
	  } else {
	    strcpy(outname_end, ".grm.gz");
	  }
	  parallel_compress(outname, calc_rel_f_grm_emitn);
	} else {
	  strcpy(outname_end, ".grm");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[4]), ".%u", parallel_idx + 1);
	  }
	  retval = write_uncompressed(outname, calc_rel_f_grm_emitn);
	  if (retval) {
	    goto calc_rel_f_ret_1;
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
	    parallel_compress(outname, calc_rel_f_tri_emitn);
	  } else {
	    retval = write_uncompressed(outname, calc_rel_f_tri_emitn);
	    if (retval) {
	      goto calc_rel_f_ret_1;
	    }
	  }
	} else if (rel_shape == REL_CALC_SQ0) {
#ifdef __LP64__
	  ulii = 0x3009300930093009LLU;
#else
	  ulii = 0x30093009LU;
#endif
	  uii = (g_indiv_ct * 2 + sizeof(intptr_t) - 4) / sizeof(intptr_t);
	  glptr2 = (uintptr_t*)g_geno;
	  for (ujj = 0; ujj < uii; ujj++) {
	    *glptr2++ = ulii;
	  }
	  g_cr_min_indiv = min_indiv;
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, calc_rel_f_sq0_emitn);
	  } else {
	    retval = write_uncompressed(outname, calc_rel_f_sq0_emitn);
	    if (retval) {
	      goto calc_rel_f_ret_1;
	    }
	  }
	} else {
	  g_cr_min_indiv = min_indiv;
	  if (rel_calc_type & REL_CALC_GZ) {
	    parallel_compress(outname, calc_rel_f_sq_emitn);
	  } else {
	    retval = write_uncompressed(outname, calc_rel_f_sq_emitn);
	    if (retval) {
	      goto calc_rel_f_ret_1;
	    }
	  }
	}
      }
    }
    putchar('\r');
    sprintf(logbuf, "Relationship matrix written to %s.\n", outname);
    logprintb();
    if (!parallel_idx) {
      strcpy(&(outname_end[4]), ".id");
      retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
      if (retval) {
	goto calc_rel_f_ret_1;
      }
    }
  }
  wkspace_reset(wkspace_mark);
  while (0) {
  calc_rel_f_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_rel_f_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  calc_rel_f_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_rel_f_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_rel_f_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_rel_f_ret_1:
  fclose_cond(outfile);
  fclose_cond(out_bin_nfile);
  gzclose_cond(gz_outfile);
  return retval;
}

int32_t calc_distance(pthread_t* threads, int32_t parallel_idx, int32_t parallel_tot, FILE* bedfile, int32_t bed_offset, FILE** outfile_ptr, char* outname, char* outname_end, uint64_t calculation_type, int32_t dist_calc_type, int32_t distance_flat_missing, uintptr_t* marker_exclude, uint32_t marker_ct, double* set_allele_freqs, uintptr_t unfiltered_indiv_ct, uintptr_t unfiltered_indiv_ct4, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len, Chrom_info* chrom_info_ptr, int32_t wt_needed, uint32_t* marker_weights_i, int32_t exp0, double exponent) {
  FILE* outfile = NULL;
  FILE* outfile2 = NULL;
  FILE* outfile3 = NULL;
  uint64_t dists_alloc = 0;
  uint32_t unwt_needed = 0;
  uint32_t unwt_needed_full = 0;
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  uint32_t chrom_fo_idx = 0;
  int32_t retval = 0;
  uint32_t* giptr = NULL;
  uint32_t* giptr2 = NULL;
  double set_allele_freq_buf[MULTIPLEX_DIST];
  uint32_t wtbuf[MULTIPLEX_DIST];
  char wbuf[16];
  char* wptr;
  unsigned char* wkspace_mark;
  unsigned char* bedbuf;
  unsigned char* gptr;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
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
  double* dptr2;
  double dxx;
  uint32_t marker_ct_autosomal;
  uint32_t multiplex;
  uint32_t chrom_end;
  int64_t llxx;
  triangle_fill(g_thread_start, g_indiv_ct, g_thread_ct, parallel_idx, parallel_tot, 1, 1);
  llxx = g_thread_start[g_thread_ct];
  llxx = ((llxx * (llxx - 1)) - (int64_t)g_thread_start[0] * (g_thread_start[0] - 1)) / 2;
  dists_alloc = llxx * sizeof(double);
#ifndef __LP64__
  if (dists_alloc > 2147483647) {
    goto calc_distance_ret_NOMEM;
  }
#endif
  // Additional + CACHELINE is to fix aliasing bug that shows up with -O2 in
  // some cases.
  // The g_missing_dbl_excluded check is to avoid recalculation, if a
  // relationship matrix was already calculated, the g_missing_dbl_excluded
  // matrix was not overwritten, and no --rel-cutoff filtering was performed
  // to make the results obsolete.  It's unimportant to keep this
  // optimization around if it ever complicates maintenance.
  if (((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX )) || distance_flat_missing) && (!g_missing_dbl_excluded)) {
    g_missing_dbl_excluded = (uint32_t*)wkspace_alloc(llxx * sizeof(int32_t));
    if (!g_missing_dbl_excluded) {
      goto calc_distance_ret_NOMEM;
    }
    fill_int_zero((int32_t*)g_missing_dbl_excluded, llxx);
    unwt_needed = 1;
    if (!g_indiv_missing_unwt) {
      g_indiv_missing_unwt = (uint32_t*)wkspace_alloc(g_indiv_ct * sizeof(int32_t));
      if (!g_indiv_missing_unwt) {
	goto calc_distance_ret_NOMEM;
      }
      unwt_needed_full = 1;
      fill_int_zero((int32_t*)g_indiv_missing_unwt, g_indiv_ct);
    }
  }
  g_dists = (double*)wkspace_alloc(dists_alloc + CACHELINE);
  if (!g_dists) {
    goto calc_distance_ret_NOMEM;
  }
  wkspace_mark = wkspace_base;
  if (wt_needed) {
    g_missing_tot_weights = (uint32_t*)wkspace_alloc(llxx * sizeof(int32_t));
    if (!g_missing_tot_weights) {
      goto calc_distance_ret_NOMEM;
    }
    fill_int_zero((int32_t*)g_missing_tot_weights, llxx);
    g_indiv_missing = (uint32_t*)wkspace_alloc(g_indiv_ct * sizeof(int32_t));
    if (!g_indiv_missing) {
      goto calc_distance_ret_NOMEM;
    }
    fill_int_zero((int32_t*)g_indiv_missing, g_indiv_ct);
  }

  if (exp0) {
    g_idists = (int32_t*)((char*)wkspace_mark - CACHEALIGN(llxx * sizeof(int32_t)));
    fill_int_zero(g_idists, llxx);
    g_masks = (uintptr_t*)wkspace_alloc(g_indiv_ct * (MULTIPLEX_2DIST / 8));
  } else {
    fill_double_zero(g_dists, llxx);
    g_masks = (uintptr_t*)wkspace_alloc(g_indiv_ct * sizeof(intptr_t));
  }
  if (!g_masks) {
    goto calc_distance_ret_NOMEM;
  }
  g_mmasks = (uintptr_t*)wkspace_alloc(g_indiv_ct * sizeof(intptr_t));
  if (!g_mmasks) {
    goto calc_distance_ret_NOMEM;
  }

  if (exp0) {
    multiplex = MULTIPLEX_DIST;
    g_geno = wkspace_alloc(g_indiv_ct * (MULTIPLEX_2DIST / 8));
  } else {
    multiplex = MULTIPLEX_DIST_EXP;
    g_geno = wkspace_alloc(g_indiv_ct * sizeof(intptr_t));
  }
  if (!g_geno) {
    goto calc_distance_ret_NOMEM;
  }

  if (wkspace_alloc_uc_checked(&bedbuf, multiplex * unfiltered_indiv_ct4)) {
    goto calc_distance_ret_NOMEM;
  }
  fseeko(bedfile, bed_offset, SEEK_SET);
  uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1);
  marker_ct_autosomal = marker_ct - uii;
  if (uii) {
    sprintf(logbuf, "Excluding %u marker%s on non-autosomes from distance matrix calc.\n", uii, (uii == 1)? "" : "s");
    logprintb();
  }
  while (marker_idx < marker_ct_autosomal) {
    for (ujj = 0; ujj < multiplex; ujj++) {
      set_allele_freq_buf[ujj] = 0.5;
    }
    fill_int_zero((int32_t*)wtbuf, multiplex);

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

    retval = block_load_autosomal(bedfile, bed_offset, marker_exclude, marker_ct_autosomal, multiplex, unfiltered_indiv_ct4, chrom_info_ptr, set_allele_freqs, marker_weights_i, bedbuf, &chrom_fo_idx, &marker_uidx, &marker_idx, &ujj, set_allele_freq_buf, NULL, wt_needed? wtbuf : NULL);
    if (retval) {
      goto calc_distance_ret_1;
    }
    if (ujj < multiplex) {
      memset(&(bedbuf[ujj * unfiltered_indiv_ct4]), 0, (multiplex - ujj) * unfiltered_indiv_ct4);
      if (exp0) {
	fill_long_zero((intptr_t*)g_geno, g_indiv_ct * (MULTIPLEX_2DIST / BITCT));
	fill_ulong_zero(g_masks, g_indiv_ct * (MULTIPLEX_2DIST / BITCT));
      } else {
	fill_long_zero((intptr_t*)g_geno, g_indiv_ct);
	fill_ulong_zero(g_masks, g_indiv_ct);
      }
    }
    if (exp0) {
      for (ukk = 0; ukk < ujj; ukk += BITCT) {
	glptr = &(((uintptr_t*)g_geno)[ukk / BITCT2]);
	glptr2 = &(g_masks[ukk / BITCT2]);
	glptr3 = g_mmasks;
	if (wt_needed) {
	  giptr = g_indiv_missing;
	}
	if (unwt_needed_full) {
	  giptr2 = g_indiv_missing_unwt;
	}
	for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
	  if (!is_set(indiv_exclude, indiv_uidx)) {
	    unn = (indiv_uidx % 4) * 2;
	    ulii = 0;
	    ulkk = 0;
	    gptr = &(bedbuf[indiv_uidx / 4 + ukk * unfiltered_indiv_ct4]);
	    for (umm = 0; umm < BITCT2; umm++) {
	      uljj = (gptr[umm * unfiltered_indiv_ct4] >> unn) & 3;
	      ulii |= uljj << (umm * 2);
	      if (uljj == 1) {
		ulkk |= ONELU << umm;
		if (wt_needed) {
		  *giptr += wtbuf[umm + ukk];
		}
		if (unwt_needed_full) {
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
	    gptr = &(bedbuf[indiv_uidx / 4 + (ukk + BITCT2) * unfiltered_indiv_ct4]);
	    for (umm = 0; umm < BITCT2; umm++) {
	      uljj = (gptr[umm * unfiltered_indiv_ct4] >> unn) & 3;
	      ulii |= uljj << (umm * 2);
	      if (uljj == 1) {
		ulkk |= ONELU << umm;
		if (wt_needed) {
		  *giptr += wtbuf[umm + ukk + BITCT2];
		}
		if (unwt_needed_full) {
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
	    if (wt_needed) {
	      giptr++;
	    }
	    if (unwt_needed_full) {
	      giptr2++;
	    }
	  }
	}

	if (wt_needed) {
	  g_weights_i = &(wtbuf[ukk]);
	  if (spawn_threads(threads, &calc_distm_thread, g_thread_ct)) {
	    goto calc_distance_ret_THREAD_CREATE_FAIL;
	  }
	  incr_wt_dist_missing(g_missing_tot_weights, 0);
	  join_threads(threads, g_thread_ct);
	}
	if (unwt_needed) {
	  if (spawn_threads(threads, &calc_missing_thread, g_thread_ct)) {
	    goto calc_distance_ret_THREAD_CREATE_FAIL;
	  }
	  incr_dists_rm(g_missing_dbl_excluded, 0, g_thread_start);
	  join_threads(threads, g_thread_ct);
	}
      }
      if (spawn_threads(threads, &calc_idist_thread, g_thread_ct)) {
	goto calc_distance_ret_THREAD_CREATE_FAIL;
      }
      incr_dists_i(g_idists, (uintptr_t*)g_geno, 0);
      join_threads(threads, g_thread_ct);
    } else {
      fill_ulong_zero(g_mmasks, g_indiv_ct);
      for (ukk = 0; ukk < ujj; ukk += MULTIPLEX_DIST_EXP / 2) {
	glptr = (uintptr_t*)g_geno;
	glptr2 = g_masks;
	glptr3 = g_mmasks;
	giptr3 = g_indiv_missing;
	for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
	  if (!is_set(indiv_exclude, indiv_uidx)) {
	    unn = (indiv_uidx % 4) * 2;
	    ulii = 0;
	    ulkk = 0;
	    gptr = &(bedbuf[indiv_uidx / 4 + ukk * unfiltered_indiv_ct4]);
	    for (umm = 0; umm < MULTIPLEX_DIST_EXP / 2; umm++) {
	      uljj = (gptr[umm * unfiltered_indiv_ct4] >> unn) & 3;
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
	fill_weights(g_weights, &(set_allele_freq_buf[ukk]), exponent);
	if (spawn_threads(threads, &calc_dist_thread, g_thread_ct)) {
	  goto calc_distance_ret_THREAD_CREATE_FAIL;
	}
	incr_dists(g_dists, (uintptr_t*)g_geno, 0);
	join_threads(threads, g_thread_ct);
      }
      g_weights_i = wtbuf;
      if (spawn_threads(threads, &calc_distm_thread, g_thread_ct)) {
	goto calc_distance_ret_THREAD_CREATE_FAIL;
      }
      incr_wt_dist_missing(g_missing_tot_weights, 0);
      join_threads(threads, g_thread_ct);
    }
    printf("\r%" PRIuPTR " markers complete.", marker_idx);
    fflush(stdout);
  }
  putchar('\r');
  logprint("Distance matrix calculation complete.\n");
  if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
    strcpy(outname_end, ".mdist");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto calc_distance_ret_OPEN_FAIL;
    }
    outfile = *outfile_ptr;
    iptr = g_idists;
    giptr = g_missing_dbl_excluded;
    pct = 1;
    // parallel_tot must be 1 for --distance-matrix
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      giptr2 = g_indiv_missing_unwt;
      uii = marker_ct_autosomal - giptr2[indiv_idx];
      for (ujj = 0; ujj < indiv_idx; ujj++) {
	wptr = double_g_writex(wbuf, ((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++))), ' ');
	fwrite(wbuf, 1, wptr - wbuf, outfile);
      }
      putc('0', outfile);
      putc(' ', outfile);
      giptr2++;
      for (ulii = indiv_idx + 1; ulii < g_indiv_ct; ulii++) {
	uljj = ulii * (ulii - 1) / 2 + indiv_idx;
	wptr = double_g_writex(wbuf, ((double)g_idists[uljj]) / (2 * (uii - (*giptr2++) + g_missing_dbl_excluded[uljj])), ' ');
	fwrite(wbuf, 1, wptr - wbuf, outfile);
      }
      if (putc('\n', outfile) == EOF) {
	goto calc_distance_ret_WRITE_FAIL;
      }
      if (indiv_idx * 100LLU >= ((uint64_t)pct * g_indiv_ct)) {
	pct = (indiv_idx * 100LLU) / g_indiv_ct;
	printf("\rWriting... %u%%", pct++);
	fflush(stdout);
      }
    }
    outfile = NULL;
    if (fclose_null(outfile_ptr)) {
      goto calc_distance_ret_WRITE_FAIL;
    }
    putchar('\r');
    sprintf(logbuf, "Distances (proportions) written to %s.\n", outname);
    logprintb();
    if (!parallel_idx) {
      strcpy(outname_end, ".mdist.id");
      retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
      if (retval) {
	goto calc_distance_ret_1;
      }
    }
  }
  if (calculation_type & CALC_PLINK_IBS_MATRIX) {
    strcpy(outname_end, ".mibs");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto calc_distance_ret_OPEN_FAIL;
    }
    outfile = *outfile_ptr;
    iptr = g_idists;
    giptr = g_missing_dbl_excluded;
    pct = 1;
    for (indiv_idx = 0; indiv_idx < g_indiv_ct; indiv_idx++) {
      giptr2 = g_indiv_missing_unwt;
      uii = marker_ct_autosomal - giptr2[indiv_idx];
      for (ujj = 0; ujj < indiv_idx; ujj++) {
	wptr = double_g_writex(wbuf, 1.0 - (((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++)))), ' ');
	fwrite(wbuf, 1, wptr - wbuf, outfile);
      }
      putc('1', outfile);
      putc(' ', outfile);
      giptr2++;
      for (ulii = indiv_idx + 1; ulii < g_indiv_ct; ulii++) {
	uljj = (ulii * (ulii - 1)) / 2 + indiv_idx;
	wptr = double_g_writex(wbuf, 1.0 - (((double)g_idists[uljj]) / (2 * (uii - (*giptr2++) + g_missing_dbl_excluded[uljj]))), ' ');
	fwrite(wbuf, 1, wptr - wbuf, outfile);
      }
      if (putc('\n', outfile) == EOF) {
	goto calc_distance_ret_WRITE_FAIL;
      }
      if (indiv_idx * 100 >= (pct * g_indiv_ct)) {
	pct = (indiv_idx * 100) / g_indiv_ct;
	printf("\rWriting... %u%%", pct++);
	fflush(stdout);
      }
    }
    outfile = NULL;
    if (fclose_null(outfile_ptr)) {
      goto calc_distance_ret_WRITE_FAIL;
    }
    putchar('\r');
    sprintf(logbuf, "IBS matrix written to %s.\n", outname);
    logprintb();
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
    if (retval) {
      goto calc_distance_ret_1;
    }
  }
  tstc = g_thread_start[g_thread_ct];
  if (wt_needed) {
    giptr = g_missing_tot_weights;
    dptr2 = g_dists;
    if (exp0) {
      iptr = g_idists;
      for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	giptr2 = g_indiv_missing;
	uii = giptr2[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  *dptr2++ = (4294967295.0 / ((4294967295U - uii - (*giptr2++)) + (*giptr++))) * (*iptr++);
	}
      }
    } else {
      for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	giptr2 = g_indiv_missing;
	uii = giptr2[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  *dptr2 *= (4294967295.0 / ((4294967295U - uii - (*giptr2++)) + (*giptr++)));
	  dptr2++;
	}
      }
    }
  } else if (distance_flat_missing) {
    dptr2 = g_dists;
    giptr = g_missing_dbl_excluded;
    if (exp0) {
      iptr = g_idists;
      for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	giptr2 = g_indiv_missing_unwt;
	uii = marker_ct_autosomal - giptr2[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  *dptr2++ = (((double)marker_ct_autosomal) / (uii - (*giptr2++) + (*giptr++))) * (*iptr++);
	}
      }
    } else {
      for (indiv_idx = g_thread_start[0]; indiv_idx < tstc; indiv_idx++) {
	giptr2 = g_indiv_missing_unwt;
	uii = marker_ct_autosomal - giptr2[indiv_idx];
	for (ujj = 0; ujj < indiv_idx; ujj++) {
	  *dptr2 *= ((double)marker_ct_autosomal) / (uii - (*giptr2++) + (*giptr++));
	  dptr2++;
	}
      }
    }
  }

  if (calculation_type & (CALC_DISTANCE | CALC_IBS_TEST)) {
    if ((exponent == 0.0) || (!(dist_calc_type & (DISTANCE_IBS | DISTANCE_1_MINUS_IBS)))) {
      g_half_marker_ct_recip = 0.5 / (double)marker_ct_autosomal;
    } else {
      g_half_marker_ct_recip = 0.0;
      marker_uidx = 0;
      chrom_fo_idx = 0;
      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[1];
      for (marker_idx = 0; marker_idx < marker_ct_autosomal; marker_idx++) {
	marker_uidx = next_autosomal_unsafe(marker_exclude, marker_uidx, chrom_info_ptr, &chrom_end, &chrom_fo_idx);
	dxx = set_allele_freqs[marker_uidx];
	if ((dxx > 0.0) && (dxx < 1.0)) {
	  g_half_marker_ct_recip += pow(2 * dxx * (1.0 - dxx), -exponent);
	} else {
	  g_half_marker_ct_recip += 1.0;
	}
	marker_uidx++;
      }
      g_half_marker_ct_recip = 0.5 / g_half_marker_ct_recip;
    }
    if (calculation_type & CALC_DISTANCE) {
      if (!parallel_idx) {
	retval = distance_d_write_ids(outname, outname_end, dist_calc_type, unfiltered_indiv_ct, indiv_exclude, person_ids, max_person_id_len);
	if (retval) {
	  goto calc_distance_ret_1;
	}
      }
      retval = distance_d_write(&outfile, &outfile2, &outfile3, dist_calc_type, outname, outname_end, g_dists, g_half_marker_ct_recip, g_indiv_ct, g_thread_start[0], g_thread_start[g_thread_ct], parallel_idx, parallel_tot, g_geno);
      if (retval) {
	goto calc_distance_ret_1;
      }
    }
  }
  wkspace_reset(wkspace_mark);

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
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 calc_distance_ret_1:
  fclose_cond(outfile);
  fclose_cond(outfile2);
  fclose_cond(outfile3);
  return retval;
}
