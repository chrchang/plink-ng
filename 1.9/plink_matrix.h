#ifndef __PLINK_MATRIX_H__
#define __PLINK_MATRIX_H__

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


// Wrappers for frequent LAPACK calls (sometimes with no-LAPACK fallbacks).
// (Update, 11 Oct 2018: Backported PLINK 2.0's MKL support.)

#ifdef NOLAPACK

#  define MATRIX_INVERT_BUF1_TYPE double
#  define MATRIX_INVERT_BUF1_ELEM_ALLOC (2 * sizeof(double))
#  define MATRIX_INVERT_BUF1_CHECKED_ALLOC (2 * sizeof(double))
#  define __CLPK_integer int

#else // !NOLAPACK

#  ifdef __APPLE__
#    include <Accelerate/Accelerate.h>
#    define USE_CBLAS_XGEMM
#  endif

#  ifndef __APPLE__

#    ifdef __cplusplus
extern "C" {
#    endif

  typedef double __CLPK_doublereal;
#    if defined(__LP64__) || defined(_WIN32)
  typedef int32_t __CLPK_integer;
#    else
  typedef long int __CLPK_integer;
#    endif

#    ifdef _WIN32
  // openblas is easy enough to set up on Windows nowadays.
  // not worth the trouble of ripping out vector extensions, etc. just so we
  // can compile with Visual Studio and gain access to MKL.
  // (todo: upgrade from 0.2.19 to a later version, build setup will probably
  // need to change a bit)
#     define HAVE_LAPACK_CONFIG_H
#     define LAPACK_COMPLEX_STRUCTURE
#     include "lapacke.h"

  __CLPK_doublereal ddot_(__CLPK_integer* n, __CLPK_doublereal* dx,
                          __CLPK_integer* incx, __CLPK_doublereal* dy,
                          __CLPK_integer* incy);

  void dger_(int* m, int* n, double* alpha, double* x, int* incx, double* y,
             int* incy, double* a, int* lda);

  void dgemm_(char* transa, char* transb, __CLPK_integer* m, __CLPK_integer* n,
              __CLPK_integer* k, double* alpha, double* a, __CLPK_integer* lda,
              double* b, __CLPK_integer* ldb, double* beta, double* c,
              __CLPK_integer* ldc);

  void dsymv_(char* uplo, int* n, double* alpha, double* a, int* lda,
              double* x, int* incx, double* beta, double* y, int* incy);

  void dgetrf_(__CLPK_integer* m, __CLPK_integer* n,
               __CLPK_doublereal* a, __CLPK_integer* lda,
               __CLPK_integer* ipiv, __CLPK_integer* info);

  void sgemm_(char* transa, char* transb, __CLPK_integer* m, __CLPK_integer* n,
              __CLPK_integer* k, float* alpha, float* a, __CLPK_integer* lda,
              float* b, __CLPK_integer* ldb, float* beta, float* c,
              __CLPK_integer* ldc);

#    else  // Linux
#      ifdef USE_MKL
#        ifndef __LP64__
#          error "32-bit Linux build does not support Intel MKL."
#        endif
#        define USE_CBLAS_XGEMM
  // sizeof(MKL_INT) should be 4.
#        define MKL_LP64
#        ifdef DYNAMIC_MKL
#          include <mkl_service.h>
#          include <mkl_cblas.h>
#          include <mkl_lapack.h>
#        else
#          include "mkl_service.h"
#          include "mkl_cblas.h"
#          include "mkl_lapack.h"
#        endif
#      else
#        ifdef USE_CBLAS_XGEMM
#          include <cblas.h>
#        else
  // ARGH
  // cmake on Ubuntu 14 seems to require use of cblas_f77.h instead of cblas.h.
  // Conversely, cblas_f77.h does not seem to be available on the Scientific
  // Linux ATLAS/LAPACK install, and right now that's my only option for
  // producing 32-bit static builds...
  // So.  Default include is cblas.h.  To play well with cmake + Ubuntu 14 and
  // 16 simultaneously, there is a CBLAS_F77_ON_OLD_GCC mode which picks
  // cblas_f77.h on Ubuntu 14 and cblas.h on 16.
#          ifdef FORCE_CBLAS_F77
#            include <cblas_f77.h>
#          elif !defined(CBLAS_F77_ON_OLD_GCC)
#            include <cblas.h>
#          else
#            if (__GNUC__ <= 4)
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
#        endif

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

  void dgels_(char* trans, __CLPK_integer* m, __CLPK_integer* n,
              __CLPK_integer* nrhs, __CLPK_doublereal* a, __CLPK_integer* lda,
              __CLPK_doublereal* b, __CLPK_integer* ldb,
              __CLPK_doublereal* work, __CLPK_integer* lwork,
              __CLPK_integer* info);

  int dsyevr_(char* jobz, char* range, char* uplo, __CLPK_integer* n,
              __CLPK_doublereal* a, __CLPK_integer* lda, __CLPK_doublereal* vl,
              __CLPK_doublereal* vu, __CLPK_integer* il, __CLPK_integer* iu,
              __CLPK_doublereal* abstol, __CLPK_integer* m,
              __CLPK_doublereal* w, __CLPK_doublereal* z, __CLPK_integer* ldz,
              __CLPK_integer* isuppz, __CLPK_doublereal* work,
              __CLPK_integer* lwork, __CLPK_integer* iwork,
              __CLPK_integer* liwork, __CLPK_integer* info);

  void dgesdd_(char* jobs, __CLPK_integer* m, __CLPK_integer* n,
               __CLPK_doublereal* a, __CLPK_integer* lda, __CLPK_doublereal* s,
               __CLPK_doublereal* u, __CLPK_integer* ldu,
               __CLPK_doublereal* vt, __CLPK_integer* ldvt,
               __CLPK_doublereal* work, __CLPK_integer* lwork,
               __CLPK_integer* iwork, __CLPK_integer* info);

#        ifndef USE_CBLAS_XGEMM
  void dgemm_(char* transa, char* transb, __CLPK_integer* m, __CLPK_integer* n,
              __CLPK_integer* k, double* alpha, double* a, __CLPK_integer* lda,
              double* b, __CLPK_integer* ldb, double* beta, double* c,
              __CLPK_integer* ldc);

  void dsymv_(char* uplo, int* n, double* alpha, double* a, int* lda,
              double* x, int* incx, double* beta, double* y, int* incy);

  void sgemm_(char* transa, char* transb, __CLPK_integer* m, __CLPK_integer* n,
              __CLPK_integer* k, float* alpha, float* a, __CLPK_integer* lda,
              float* b, __CLPK_integer* ldb, float* beta, float* c,
              __CLPK_integer* ldc);
#        endif

#      endif  // !USE_MKL
#    endif


  void xerbla_(void);
#    ifdef __cplusplus
} // extern "C"
#    endif

#  endif  // !__APPLE__

#  define MATRIX_INVERT_BUF1_TYPE __CLPK_integer
#  define MATRIX_INVERT_BUF1_ELEM_ALLOC sizeof(__CLPK_integer)
// invert_matrix_checked() usually requires a larger buffer
#  define MATRIX_INVERT_BUF1_CHECKED_ALLOC (2 * sizeof(__CLPK_integer))

#endif  // !NOLAPACK

#define MATRIX_SINGULAR_RCOND 1e-14

#ifdef NOLAPACK
int32_t invert_matrix(int32_t dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* dbl_1d_buf, double* dbl_2d_buf);

#  define invert_matrix_checked invert_matrix
#else
int32_t invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf);

int32_t invert_matrix_checked(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf);
#endif

void col_major_matrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* inmatrix1, double* inmatrix2, double* outmatrix);

void col_major_fmatrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, float* inmatrix1, float* inmatrix2, float* outmatrix);

void transpose_copy(uintptr_t old_maj, uintptr_t new_maj, double* old_matrix, double* new_matrix);

void transpose_copy_float(uintptr_t old_maj, uintptr_t new_maj, uint32_t new_maj_max, float* old_matrix, float* new_matrix);

#endif // __PLINK_MATRIX_H__
