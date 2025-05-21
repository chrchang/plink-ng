#ifndef __PLINK2_MATRIX_H__
#define __PLINK2_MATRIX_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2025 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Wrappers for frequent LAPACK calls (sometimes with no-LAPACK fallbacks).
// Now supports MKL backend.

// todo: allow this to take advantage of 64-bit integer LAPACK.  As of this
// writing, it's available on Amazon EC2 64-bit Linux instances, but I can't
// find it for Windows.  (And even if OS X vecLib adds it soon, we can't use it
// there anytime soon because static linking is not an option.)

#include "plink2_cmdline.h"

#ifdef NOLAPACK
typedef double MatrixInvertBuf1;
CONSTI32(kMatrixInvertBuf1ElemAlloc, 2 * sizeof(double));
CONSTI32(kMatrixInvertBuf1CheckedAlloc, 2 * sizeof(double));
#  define __CLPK_integer int

#else  // not NOLAPACK
#  ifdef __APPLE__
#    include <Accelerate/Accelerate.h>
#    define USE_CBLAS_XGEMM
#    ifdef ACCELERATE_NEW_LAPACK
typedef int32_t __CLPK_integer;
#    endif
#  elif defined(USE_AOCL)
#    define USE_CBLAS_XGEMM
#  endif

#  ifndef __APPLE__

#    ifdef __cplusplus
extern "C" {
#    endif
  typedef double __CLPK_doublereal;
#    ifdef __LP64__
#      ifdef LAPACK_ILP64
  typedef long long __CLPK_integer;
#      else
  typedef int32_t __CLPK_integer;
#      endif
#    else
#      ifdef LAPACK_ILP64
#        error "Invalid compile flags."
#      else
#        ifdef _WIN32
  // probably don't need this?
  typedef int32_t __CLPK_integer;
#        else
  typedef long int __CLPK_integer;
#        endif
#      endif
#    endif  // !__LP64__

#    ifdef _WIN32
  // openblas is easy enough to set up on Windows nowadays.
  // not worth the trouble of ripping out vector extensions, etc. just so we
  // can compile with Visual Studio and gain access to MKL
#      ifndef USE_OPENBLAS
#        error "Windows build currently requires OpenBLAS's LAPACK."
#      endif
#      define HAVE_LAPACK_CONFIG_H
#      define LAPACK_COMPLEX_STRUCTURE
#      include "lapacke.h"

  __CLPK_doublereal ddot_(__CLPK_integer* n, __CLPK_doublereal* dx,
                          __CLPK_integer* incx, __CLPK_doublereal* dy,
                          __CLPK_integer* incy);
  __CLPK_doublereal sdot_(__CLPK_integer* n, float* sx, __CLPK_integer* incx,
                          float* sy, __CLPK_integer* incy);
#    else  // Linux
#      ifdef USE_MKL
#        define USE_CBLAS_XGEMM
#        ifdef DYNAMIC_MKL
#          include <mkl_cblas.h>
#          include <mkl_lapack.h>
#        else
#          include "mkl_cblas.h"
#          include "mkl_lapack.h"
#        endif
static_assert(sizeof(MKL_INT) == 8, "Unexpected MKL_INT size.");
#      else
// If you want 64-bit index support, but not MKL (e.g. you're targeting an
// AMD processor), modify the Makefile to link to a LAPACK library recompiled
// with -fdefault-integer-8.

#        ifdef USE_CBLAS_XGEMM
#          include <cblas.h>
// LAPACK 3.9.1 changed function signatures for functions taking char*
// parameters; previous signatures actually corresponded to undefined behavior
// which manifested in real problems with gfortran 7-9.
#          ifdef LAPACK_FORTRAN_STRLEN_END
  int dpotri_(char* uplo, __CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_integer* info, size_t slen1);
#          else
  int dpotri_(char* uplo, __CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_integer* info);
#          endif
#        else
          // ARGH
          // cmake on Ubuntu 14 seems to require use of cblas_f77.h instead of
          // cblas.h.  Conversely, cblas_f77.h does not seem to be available on
          // the Scientific Linux ATLAS/LAPACK install, and right now that's my
          // only option for producing 32-bit static builds...
          // So.  Default include is cblas.h.  To play well with cmake + Ubuntu
          // 14 and 16 simultaneously, there is a CBLAS_F77_ON_OLD_GCC mode
          // which picks cblas_f77.h on Ubuntu 14 and cblas.h on 16.
#          ifdef FORCE_CBLAS_F77
#            include <cblas_f77.h>
#          elif !defined(CBLAS_F77_ON_OLD_GCC)
#            include <cblas.h>
#          else
#            if (__GNUC__ == 4)
#              include <cblas_f77.h>
#            else
#              if __has_include(<cblas.h>)
#                include <cblas.h>
#              else
#                include <cblas_f77.h>
#              endif
#            endif
#          endif
  __CLPK_doublereal ddot_(__CLPK_integer* n, __CLPK_doublereal* dx,
                          __CLPK_integer* incx, __CLPK_doublereal* dy,
                          __CLPK_integer* incy);

  __CLPK_doublereal sdot_(__CLPK_integer* n, float* sx, __CLPK_integer* incx,
                          float* sy, __CLPK_integer* incy);
  int dpotri_(char* uplo, __CLPK_integer* n, __CLPK_doublereal* a,
              __CLPK_integer* lda, __CLPK_integer* info);
#        endif
#      endif  // !USE_MKL
#      ifdef USE_CUDA
#        include "cuda/plink2_matrix_cuda.h"
#      endif
#    endif

  void xerbla_(void);
#    ifdef __cplusplus
} // extern "C"
#    endif

#  endif  // !__APPLE__

#  ifdef LAPACK_FORTRAN_STRLEN_END
#    define dpotri_wrap(a, b, c, d, e) dpotri_((a), (b), (c), (d), (e), 1)
#  else
#    define dpotri_wrap(a, b, c, d, e) dpotri_((a), (b), (c), (d), (e))
#  endif

typedef __CLPK_integer MatrixInvertBuf1;
// need to be careful about >= 2^32?
CONSTI32(kMatrixInvertBuf1ElemAlloc, sizeof(__CLPK_integer));
// invert_matrix_checked() usually requires a larger buffer
CONSTI32(kMatrixInvertBuf1CheckedAlloc, 2 * sizeof(__CLPK_integer));

#endif  // !NOLAPACK

#ifdef __cplusplus
namespace plink2 {
#endif

static const double kMatrixSingularRcond = 1e-14;

// Returns -1 if no inf/nan found.
// May move this to a more central location if there are other users.
intptr_t FirstInfOrNan(const double* vec, uintptr_t size);

uint32_t LowerTriangularFirstInfOrNan(const double* matrix, uintptr_t dim, uintptr_t* row_idx_ptr, uintptr_t* col_idx_ptr);

// Copies (C-order) lower-left to upper right.
void ReflectMatrix(uint32_t dim, double* matrix);

void ReflectStridedMatrix(uint32_t dim, uint32_t stride, double* matrix);

void ReflectFmatrix(uint32_t dim, uint32_t stride, float* matrix);

// If dim < stride, this zeroes out the trailing elements of each row.
void ReflectStridedMatrix0(uint32_t dim, uint32_t stride, double* matrix);

void ReflectFmatrix0(uint32_t dim, uint32_t stride, float* matrix);

// AddDVec and AddFVec execute main += arg, under the assumption that both
// vectors have been zero-padded.
// "ctav" is the float-length of both vectors, rounded up to a vector boundary.
// (With the current implementation, it's harmless to just pass in ct instead,
// but ctav is better practice since it indicates awareness of what this
// function is actually doing.)
HEADER_INLINE void AddDVec(const double* arg, uintptr_t ctav, double* main) {
  // LEA instruction supports multiply by 4 or 8, but *not* by 16 or 32, so
  // this awkward loop tends to generate slightly better code than the more
  // natural vector-based loop.
  for (uintptr_t ulii = 0; ulii < ctav; ulii += kDoublePerDVec) {
    *R_CAST(VecD*, &(main[ulii])) += *R_CAST(const VecD*, &(arg[ulii]));
  }
}

HEADER_INLINE void AddFVec(const float* arg, uintptr_t ctav, float* main) {
  for (uintptr_t ulii = 0; ulii < ctav; ulii += kFloatPerFVec) {
    *R_CAST(VecF*, &(main[ulii])) += *R_CAST(const VecF*, &(arg[ulii]));
  }
}

#ifdef __LP64__
// FillDVec sets all elements of dst to dxx, and zero-fills any trailing
// elements (w.r.t. vector boundaries).
// bugfix (21 Nov 2023): forgot about trailing zero-fill.
HEADER_INLINE void FillDVec(uintptr_t ct, double dxx, double* dst) {
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst[ulii] = dxx;
  }
  const uintptr_t remainder = ct % kDoublePerDVec;
  if (remainder) {
    const uintptr_t ctav = ct + kDoublePerDVec - remainder;
    for (uintptr_t ulii = ct; ulii != ctav; ++ulii) {
      dst[ulii] = 0.0;
    }
  }
}
#else
HEADER_INLINE void FillDVec(uintptr_t ct, double dxx, double* dst) {
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst[ulii] = dxx;
  }
}
#endif

#ifdef __LP64__
// FillFVec sets all elements of dst to fxx, and zero-fills any trailing
// elements (w.r.t. vector boundaries).
HEADER_INLINE void FillFVec(uintptr_t ct, float fxx, float* dst) {
  const uintptr_t fullvec_ct = ct / kFloatPerFVec;
  const VecF vfill = VCONST_F(fxx);
  VecF* dst_alias = R_CAST(VecF*, dst);
  for (uintptr_t ulii = 0; ulii != fullvec_ct; ++ulii) {
    dst_alias[ulii] = vfill;
  }
  const uintptr_t trailing_start_idx = fullvec_ct * kFloatPerFVec;
  if (trailing_start_idx != ct) {
    for (uintptr_t ulii = trailing_start_idx; ulii != ct; ++ulii) {
      dst[ulii] = fxx;
    }
    const uintptr_t trailing_stop_idx = trailing_start_idx + kFloatPerFVec;
    for (uintptr_t ulii = ct; ulii != trailing_stop_idx; ++ulii) {
      dst[ulii] = S_CAST(float, 0.0);
    }
  }
}
#else
HEADER_INLINE void FillFVec(uintptr_t ct, float fxx, float* dst) {
  for (uintptr_t ulii = 0; ulii != ct; ++ulii) {
    dst[ulii] = fxx;
  }
}
#endif

HEADER_INLINE double DotprodDShort(const double* vec1, const double* vec2, uint32_t ct) {
  double dotprod = 0.0;
  for (uint32_t uii = 0; uii != ct; ++uii) {
    dotprod += vec1[uii] * vec2[uii];
  }
  return dotprod;
}

HEADER_INLINE float DotprodFShort(const float* vec1, const float* vec2, uint32_t ct) {
  float dotprod = 0.0;
  for (uint32_t uii = 0; uii != ct; ++uii) {
    dotprod += vec1[uii] * vec2[uii];
  }
  return dotprod;
}

// todo: benchmark again after Spectre/Meltdown mitigation is deployed
CONSTI32(kDotprodDThresh, 17);
CONSTI32(kDotprodFThresh, 15);

#ifdef NOLAPACK
HEADER_INLINE double DotprodD(const double* vec1, const double* vec2, uint32_t ct) {
  return DotprodDShort(vec1, vec2, ct);
}

HEADER_INLINE double DotprodxD(const double* vec1, const double* vec2, uint32_t ct) {
  return DotprodDShort(vec1, vec2, ct);
}

HEADER_INLINE float DotprodF(const float* vec1, const float* vec2, uint32_t ct) {
  return DotprodFShort(vec1, vec2, ct);
}

BoolErr InvertMatrix(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf);

HEADER_INLINE BoolErr InvertMatrixChecked(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  return InvertMatrix(dim, matrix, dbl_1d_buf, dbl_2d_buf);
}

HEADER_INLINE BoolErr InvertSymmdefMatrix(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  ReflectMatrix(dim, matrix);
  return InvertMatrix(dim, matrix, dbl_1d_buf, dbl_2d_buf);
}

HEADER_INLINE BoolErr InvertSymmdefMatrixChecked(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  ReflectMatrix(dim, matrix);
  return InvertMatrix(dim, matrix, dbl_1d_buf, dbl_2d_buf);
}

// first half computes either LU or singular value decomposition, and
//   determinant
// second half actually inverts matrix, assuming 1d_buf and 2d_buf have results
//   from first half
BoolErr InvertStridedMatrixFirstHalf(int32_t dim, int32_t stride, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf);

HEADER_INLINE BoolErr InvertSymmdefMatrixFirstHalf(int32_t dim, int32_t stride, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  ReflectStridedMatrix(dim, stride, matrix);
  return InvertStridedMatrixFirstHalf(dim, stride, matrix, dbl_1d_buf, dbl_2d_buf);
}

BoolErr InvertFmatrixFirstHalf(int32_t dim, uint32_t stride, const float* matrix, double* half_inverted, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf);

HEADER_INLINE BoolErr InvertSymmdefFmatrixFirstHalf(int32_t dim, uint32_t stride, float* matrix, double* half_inverted, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  ReflectFmatrix(dim, stride, matrix);
  return InvertFmatrixFirstHalf(dim, stride, matrix, half_inverted, dbl_1d_buf, dbl_2d_buf);
}

HEADER_INLINE double HalfInvertedDet(__maybe_unused const double* half_inverted_iter, __maybe_unused const MatrixInvertBuf1* dbl_1d_buf, uint32_t dim) {
  // singular values in dbl_1d_buf
  double det_u = dbl_1d_buf[0];
  for (uint32_t uii = 1; uii != dim; ++uii) {
    det_u *= dbl_1d_buf[uii];
  }
  return fabs(det_u);
}

HEADER_INLINE double HalfSymmInvertedDet(__maybe_unused const double* half_inverted_iter, __maybe_unused const MatrixInvertBuf1* dbl_1d_buf, uint32_t dim, __maybe_unused uint32_t stride) {
  return HalfInvertedDet(half_inverted_iter, dbl_1d_buf, dim);
}

void InvertStridedMatrixSecondHalf(__CLPK_integer dim, __CLPK_integer stride, double* inverted_result, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf);

HEADER_INLINE void InvertSymmdefMatrixSecondHalf(__CLPK_integer dim, __CLPK_integer stride, double* inverted_result, __maybe_unused MatrixInvertBuf1* dbl_1d_buf, __maybe_unused double* dbl_2d_buf) {
  InvertStridedMatrixSecondHalf(dim, stride, inverted_result, dbl_1d_buf, dbl_2d_buf);
}

void InvertFmatrixSecondHalf(__CLPK_integer dim, uint32_t stride, double* half_inverted, float* inverted_result, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf);

HEADER_INLINE void InvertSymmdefFmatrixSecondHalf(__CLPK_integer dim, uint32_t stride, double* half_inverted, float* inverted_result, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  InvertFmatrixSecondHalf(dim, stride, half_inverted, inverted_result, dbl_1d_buf, dbl_2d_buf);
}

HEADER_INLINE BoolErr InvertSymmdefStridedMatrix(__CLPK_integer dim, __CLPK_integer stride, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  if (InvertSymmdefMatrixFirstHalf(dim, stride, matrix, dbl_1d_buf, dbl_2d_buf)) {
    return 1;
  }
  InvertSymmdefMatrixSecondHalf(dim, stride, matrix, dbl_1d_buf, dbl_2d_buf);
  return 0;
}

#else // !NOLAPACK
HEADER_INLINE double DotprodD(const double* vec1, const double* vec2, uint32_t ct) {
#  ifndef USE_CBLAS_XGEMM
  __CLPK_integer cti = ct;
  __CLPK_integer incxy = 1;
  return ddot_(&cti, K_CAST(double*, vec1), &incxy, K_CAST(double*, vec2), &incxy);
#  else
  return cblas_ddot(ct, vec1, 1, vec2, 1);
#  endif
}

HEADER_INLINE double DotprodxD(const double* vec1, const double* vec2, uint32_t ct) {
  if (ct > kDotprodDThresh) {
    // best threshold is machine-dependent; this is what I got on my MacBook
    return DotprodD(vec1, vec2, ct);
  }
  return DotprodDShort(vec1, vec2, ct);
}

// not worthwhile for ct < 16.
HEADER_INLINE float DotprodF(const float* vec1, const float* vec2, uint32_t ct) {
#  ifndef USE_CBLAS_XGEMM
  __CLPK_integer cti = ct;
  __CLPK_integer incxy = 1;
  return S_CAST(float, sdot_(&cti, K_CAST(float*, vec1), &incxy, K_CAST(float*, vec2), &incxy));
#  else
  return cblas_sdot(ct, vec1, 1, vec2, 1);
#  endif
}

// extra if-statement in DotprodxF() seems disproportionally expensive in
// test?... guess I won't have auto-branch for now.

BoolErr InvertMatrix(__CLPK_integer dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

BoolErr InvertMatrixChecked(__CLPK_integer dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);


// InvertSymmdef... functions only assume (C-order) lower left of matrix is
// filled, and only the lower left of the return matrix is valid.
BoolErr InvertSymmdefMatrix(__CLPK_integer dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

BoolErr InvertSymmdefStridedMatrix(__CLPK_integer dim, __CLPK_integer stride, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

// dbl_2d_buf must have room for at least max(dim, 3) * dim elements.
BoolErr InvertSymmdefMatrixChecked(__CLPK_integer dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

BoolErr InvertSymmdefMatrixFirstHalf(__CLPK_integer dim, __CLPK_integer stride, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);


BoolErr InvertFmatrixFirstHalf(__CLPK_integer dim, uint32_t stride, const float* matrix, double* half_inverted, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

BoolErr InvertSymmdefFmatrixFirstHalf(__CLPK_integer dim, uint32_t stride, float* matrix, double* half_inverted, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

/*
HEADER_INLINE double HalfInvertedDet(__maybe_unused const double* half_inverted_iter, __maybe_unused const MatrixInvertBuf1* int_1d_buf, uint32_t dim) {
  uint32_t dim_p1 = dim + 1;
  double det_u = *half_inverted_iter;
  for (uint32_t uii = 1; uii != dim; ++uii) {
    half_inverted_iter = &(half_inverted_iter[dim_p1]);
    det_u *= (*half_inverted_iter);
  }
  return fabs(det_u);
}
*/

HEADER_INLINE double HalfSymmInvertedDet(__maybe_unused const double* half_inverted_iter, __maybe_unused const MatrixInvertBuf1* int_1d_buf, uint32_t dim, uint32_t stride) {
  uint32_t stride_p1 = stride + 1;
  double sqrt_det_u = *half_inverted_iter;
  for (uint32_t uii = 1; uii != dim; ++uii) {
    half_inverted_iter = &(half_inverted_iter[stride_p1]);
    sqrt_det_u *= (*half_inverted_iter);
  }
  return sqrt_det_u * sqrt_det_u;
}

HEADER_INLINE void InvertSymmdefMatrixSecondHalf(__CLPK_integer dim, __CLPK_integer stride, double* matrix, __maybe_unused MatrixInvertBuf1* int_1d_buf, __maybe_unused double* dbl_2d_buf) {
  char uplo = 'U';
  __CLPK_integer info;
  dpotri_wrap(&uplo, &dim, matrix, &stride, &info);
}

void InvertFmatrixSecondHalf(__CLPK_integer dim, uint32_t stride, double* half_inverted, float* inverted_result, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);

void InvertSymmdefFmatrixSecondHalf(__CLPK_integer dim, uint32_t stride, double* half_inverted, float* inverted_result, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf);
#endif // !NOLAPACK

void ColMajorMatrixMultiply(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* outmatrix);

HEADER_INLINE void RowMajorMatrixMultiply(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* outmatrix) {
  return ColMajorMatrixMultiply(inmatrix2, inmatrix1, col2_ct, row1_ct, common_ct, outmatrix);
}

// this is essentially a full-blown dgemm wrapper, only missing the alpha
// parameter now
void ColMajorMatrixMultiplyStridedAddassign(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double beta, double* outmatrix);

HEADER_INLINE void ColMajorMatrixMultiplyStrided(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double* outmatrix) {
  ColMajorMatrixMultiplyStridedAddassign(inmatrix1, inmatrix2, row1_ct, stride1, col2_ct, stride2, common_ct, stride3, 0.0, outmatrix);
}

HEADER_INLINE void RowMajorMatrixMultiplyIncr(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* outmatrix) {
  return ColMajorMatrixMultiplyStridedAddassign(inmatrix2, inmatrix1, col2_ct, col2_ct, row1_ct, common_ct, common_ct, col2_ct, 1.0, outmatrix);
}

HEADER_INLINE void RowMajorMatrixMultiplyStrided(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double* outmatrix) {
  // stride1 should be close to common_ct
  // stride2 should be close to col2_ct
  // output matrix uses stride3, which should be close to col2_ct
  return ColMajorMatrixMultiplyStridedAddassign(inmatrix2, inmatrix1, col2_ct, stride2, row1_ct, stride1, common_ct, stride3, 0.0, outmatrix);
}

HEADER_INLINE void RowMajorMatrixMultiplyStridedIncr(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double* outmatrix) {
  return ColMajorMatrixMultiplyStridedAddassign(inmatrix2, inmatrix1, col2_ct, stride2, row1_ct, stride1, common_ct, stride3, 1.0, outmatrix);
}

// out^T := V^T * M
// out := M^T * V
void ColMajorVectorMatrixMultiplyStrided(const double* in_dvec1, const double* inmatrix2, __CLPK_integer common_ct, __CLPK_integer stride2, __CLPK_integer col2_ct, double* out_dvec);

// out := M * V
void ColMajorMatrixVectorMultiplyStrided(const double* inmatrix1, const double* in_dvec2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer common_ct, double* out_dvec);

void ColMajorFmatrixMultiplyStrided(const float* inmatrix1, const float* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, float* outmatrix);

HEADER_INLINE void RowMajorFmatrixMultiply(const float* inmatrix1, const float* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, float* outmatrix) {
  ColMajorFmatrixMultiplyStrided(inmatrix2, inmatrix1, col2_ct, col2_ct, row1_ct, common_ct, common_ct, col2_ct, outmatrix);
}

// out := M * V
void ColMajorFmatrixVectorMultiplyStrided(const float* inmatrix1, const float* in_fvec2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer common_ct, float* out_fvec);

// out^T := V^T * M
// out := M^T * V
void ColMajorFvectorMatrixMultiplyStrided(const float* in_fvec1, const float* inmatrix2, __CLPK_integer common_ct, __CLPK_integer stride2, __CLPK_integer col2_ct, float* out_fvec);

void MatrixTransposeCopy(const double* old_matrix, uint32_t old_maj, uint32_t new_maj, double* new_matrix_iter);

void FmatrixTransposeCopy(const float* old_matrix, uint32_t old_maj, uint32_t new_maj, uint32_t new_maj_max, float* new_matrix_iter);


// A(A^T), where A is row-major; result is dim x dim
// ONLY UPDATES LOWER TRIANGLE OF result[].
void MultiplySelfTranspose(const double* input_matrix, uint32_t dim, uint32_t col_ct, double* result);

void MultiplySelfTransposeStrided(const double* input_matrix, uint32_t dim, uint32_t col_ct, uint32_t stride, double* result);

void MultiplySelfTransposeStridedF(const float* input_matrix, uint32_t dim, uint32_t col_ct, uint32_t stride, float* result);


// (A^T)A
void TransposeMultiplySelfIncr(double* input_part, uint32_t dim, uint32_t partial_row_ct, double* result);

#ifndef NOLAPACK
BoolErr GetSvdRectLwork(uint32_t major_ct, uint32_t minor_ct, __CLPK_integer* lwork_ptr);

// currently a wrapper for dgesvd_().
IntErr SvdRect(uint32_t major_ct, uint32_t minor_ct, __CLPK_integer lwork, double* matrix, double* ss, unsigned char* svd_rect_wkspace);

BoolErr GetExtractEigvecsLworks(uint32_t dim, uint32_t pc_ct, __CLPK_integer* lwork_ptr, __CLPK_integer* liwork_ptr, uintptr_t* wkspace_byte_ct_ptr);

// currently a wrapper for dsyevr_().  Matrix is expected to be
// lower-triangular from C's perspective.
// reverse_eigvecs is eigenvector-major, but the vectors are in order of
// *increasing* eigenvalue.
BoolErr ExtractEigvecs(uint32_t dim, uint32_t pc_ct, __CLPK_integer lwork, __CLPK_integer liwork, double* matrix, double* eigvals, double* reverse_eigvecs, unsigned char* extract_eigvecs_wkspace);
#endif

// Computes inverse of
//   [ A   b^T ]
//   [ b   c   ]
// given precomputed A^{-1} (must be fully filled out, not just lower left).
// See e.g.
//   https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion .
// Only fills lower left of outmatrix.
// insert_idx specifies the zero-based row/column number of b/c in outmatrix.
BoolErr InvertRank1Symm(const double* a_inv, const double* bb, __CLPK_integer orig_dim, uint32_t insert_idx, double cc, double* __restrict outmatrix, double* __restrict ainv_b_buf);

// When you only need the diagonal from InvertRank1Symm().  insert_idx
// assumed to be orig_dim.
BoolErr InvertRank1SymmDiag(const double* a_inv, const double* bb, __CLPK_integer orig_dim, double cc, double* __restrict outdiag, double* __restrict ainv_b_buf);

//   [ A   B^T ]
//   [ B   D   ]
// B is row-major.
BoolErr InvertRank2Symm(const double* a_inv, const double* bb, __CLPK_integer orig_dim, __CLPK_integer b_stride, uint32_t insert_idx, double d11, double d12, double d22, double* __restrict outmatrix, double* __restrict b_ainv_buf, double* __restrict s_b_ainv_buf);

BoolErr InvertRank2SymmDiag(const double* a_inv, const double* bb, __CLPK_integer orig_dim, double d11, double d12, double d22, double* __restrict outdiag, double* __restrict b_ainv_buf, double* __restrict s_b_ainv_buf);

#ifdef NOLAPACK
BoolErr LinearRegressionInvMain(const double* xt_y_phenomaj, uint32_t predictor_ct, uint32_t pheno_ct, double* xtx_inv, double* fitted_coefs_phenomaj, MatrixInvertBuf1* mi_buf, double* dbl_2d_buf);
#else
BoolErr LinearRegressionInvMain(const double* xt_y_phenomaj, uint32_t predictor_ct, __CLPK_integer pheno_ct, double* xtx_inv, double* fitted_coefs_phenomaj);
#endif

// now assumes xtx_inv is predictors_pmaj * transpose on input
HEADER_INLINE BoolErr LinearRegressionInv(const double* pheno_d_pmaj, const double* predictors_pmaj, uint32_t predictor_ct, uint32_t sample_ct, uint32_t pheno_ct, double* xtx_inv, double* fitted_coefs_phenomaj, double* xt_y_phenomaj, __maybe_unused MatrixInvertBuf1* mi_buf, __maybe_unused double* dbl_2d_buf) {
  // MultiplySelfTranspose(predictors_pmaj, predictor_ct, sample_ct,
  //   xtx_inv);
  // categorical optimization possible here
  for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
    RowMajorMatrixMultiply(predictors_pmaj, &(pheno_d_pmaj[pheno_idx * sample_ct]), predictor_ct, 1, sample_ct, &(xt_y_phenomaj[pheno_idx * predictor_ct]));
  }
#ifdef NOLAPACK
  return LinearRegressionInvMain(xt_y_phenomaj, predictor_ct, pheno_ct, xtx_inv, fitted_coefs_phenomaj, mi_buf, dbl_2d_buf);
#else
  return LinearRegressionInvMain(xt_y_phenomaj, predictor_ct, pheno_ct, xtx_inv, fitted_coefs_phenomaj);
#endif
}

HEADER_INLINE BoolErr LinearRegressionDVec(const double* pheno_d, const double* predictors_pmaj, uint32_t predictor_ct, uint32_t sample_ct, double* fitted_coefs, double* xtx_inv_buf, double* xt_y_buf, __maybe_unused MatrixInvertBuf1* mi_buf, __maybe_unused double* dbl_2d_buf) {
  // categorical optimization possible here
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kDoublePerDVec);
  ColMajorVectorMatrixMultiplyStrided(pheno_d, predictors_pmaj, sample_ct, sample_ctav, predictor_ct, xt_y_buf);
  MultiplySelfTransposeStrided(predictors_pmaj, predictor_ct, sample_ct, sample_ctav, xtx_inv_buf);
#ifdef NOLAPACK
  return LinearRegressionInvMain(xt_y_buf, predictor_ct, 1, xtx_inv_buf, fitted_coefs, mi_buf, dbl_2d_buf);
#else
  return LinearRegressionInvMain(xt_y_buf, predictor_ct, 1, xtx_inv_buf, fitted_coefs);
#endif
}

// just for debugging
HEADER_INLINE void PrintVector(const double* vec, uintptr_t ct) {
  printf("%g", vec[0]);
  for (uintptr_t ulii = 1; ulii != ct; ++ulii) {
    printf(" %g", vec[ulii]);
  }
  printf("\n");
}

HEADER_INLINE void PrintFvector(const float* vec, uintptr_t ct) {
  printf("%g", S_CAST(double, vec[0]));
  for (uintptr_t ulii = 1; ulii != ct; ++ulii) {
    printf(" %g", S_CAST(double, vec[ulii]));
  }
  printf("\n");
}

HEADER_INLINE void PrintSymmMatrix(const double* matrix, uint32_t dim) {
  for (uint32_t uii = 0; uii != dim; ++uii) {
    printf("%g", matrix[uii * dim]);
    for (uint32_t ujj = 1; ujj <= uii; ++ujj) {
      printf(" %g", matrix[uii * dim + ujj]);
    }
    printf("\n");
  }
}

HEADER_INLINE void PrintSymmFmatrix(const float* matrix, uint32_t dim) {
  for (uint32_t uii = 0; uii != dim; ++uii) {
    printf("%g", S_CAST(double, matrix[uii * dim]));
    for (uint32_t ujj = 1; ujj <= uii; ++ujj) {
      printf(" %g", S_CAST(double, matrix[uii * dim + ujj]));
    }
    printf("\n");
  }
}

HEADER_INLINE void PrintMatrix(const double* matrix, uintptr_t row_ct, uintptr_t col_ct) {
  for (uintptr_t ulii = 0; ulii != row_ct; ++ulii) {
    printf("%g", matrix[ulii * col_ct]);
    for (uintptr_t uljj = 1; uljj != col_ct; ++uljj) {
      printf(" %g", matrix[ulii * col_ct + uljj]);
    }
    printf("\n");
  }
}

HEADER_INLINE void PrintStridedMatrix(const double* matrix, uintptr_t row_ct, uintptr_t col_ct, uintptr_t stride) {
  for (uintptr_t ulii = 0; ulii != row_ct; ++ulii) {
    printf("%g", matrix[ulii * stride]);
    for (uintptr_t uljj = 1; uljj != col_ct; ++uljj) {
      printf(" %g", matrix[ulii * stride + uljj]);
    }
    printf("\n");
  }
}

HEADER_INLINE void PrintFmatrix(const float* matrix, uintptr_t row_ct, uintptr_t col_ct, uintptr_t stride) {
  for (uintptr_t ulii = 0; ulii != row_ct; ++ulii) {
    printf("%g", S_CAST(double, matrix[ulii * stride]));
    for (uintptr_t uljj = 1; uljj != col_ct; ++uljj) {
      printf(" %g", S_CAST(double, matrix[ulii * stride + uljj]));
    }
    printf("\n");
  }
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MATRIX_H__
