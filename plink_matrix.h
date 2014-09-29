#ifndef __PLINK_MATRIX_H__
#define __PLINK_MATRIX_H__

// Wrappers for frequent LAPACK calls (sometimes with no-LAPACK fallbacks).
// May want to make this comprehensive to make linking with Intel MKL practical
// in the future.

// todo: allow this to take advantage of 64-bit integer LAPACK.  As of this
// writing, it's available on Amazon EC2 64-bit Linux instances, but I can't
// find it for Windows.  (And even if OS X vecLib adds it soon, we can't use it
// there anytime soon because static linking is not an option.)

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

  void sgemm_(char* transa, char* transb, int* m, int* n, int* k,
              float* alpha, float* a, int* lda, float* b, int* ldb,
              float* beta, float* c, int* ldc);

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
#endif

  void xerbla_(void);
#ifdef __cplusplus
}
#endif // __cplusplus
#endif // __APPLE__
#define MATRIX_INVERT_BUF1_TYPE __CLPK_integer
#endif // NOLAPACK

#define MATRIX_SINGULAR_RCOND 1e-14

#ifdef NOLAPACK
int32_t invert_matrix(int32_t dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* dbl_1d_buf, double* dbl_2d_buf);

#define invert_matrix_checked invert_matrix
#else
int32_t invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf);

int32_t invert_matrix_checked(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf);
#endif

void col_major_matrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* inmatrix1, double* inmatrix2, double* outmatrix);

void col_major_fmatrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, float* inmatrix1, float* inmatrix2, float* outmatrix);

void transpose_copy(uintptr_t old_maj, uintptr_t new_maj, double* old_matrix, double* new_matrix);

void transpose_copy_float(uintptr_t old_maj, uintptr_t new_maj, uint32_t new_maj_max, float* old_matrix, float* new_matrix);

#endif // __PLINK_MATRIX_H__
