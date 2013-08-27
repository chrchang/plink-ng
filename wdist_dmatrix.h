#ifndef __WDIST_DMATRIX_H__
#define __WDIST_DMATRIX_H__

#include "wdist_common.h"

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
#else
int32_t invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf);
#endif

int32_t invert_matrix_trunc_singular(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim);

void col_major_matrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* inmatrix1, double* inmatrix2, double* outmatrix);

#endif // __WDIST_DMATRIX_H__
