#ifndef __PLINK2_MATRIX_H__
#define __PLINK2_MATRIX_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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

#include "plink2_common.h"

#ifdef NOLAPACK

  typedef double matrix_invert_buf1_t;
  CONSTU31(kMatrixInvertBuf1ElemAlloc, 2 * sizeof(double));
  CONSTU31(kMatrixInvertBuf1CheckedAlloc, 2 * sizeof(double));
  #define __CLPK_integer int
  typedef float matrix_finvert_buf1_t;
  CONSTU31(kMatrixFinvertBuf1CheckedAlloc, 2 * sizeof(float));

#else // not NOLAPACK
  #ifdef __APPLE__
    #include <Accelerate/Accelerate.h>
  #endif

  #ifndef __APPLE__

    #ifdef __cplusplus
extern "C" {
    #endif
  typedef double __CLPK_doublereal;
    #ifdef __LP64__
      #ifdef LAPACK_ILP64
  typedef long long __CLPK_integer;
      #else
  typedef int32_t __CLPK_integer;
      #endif
    #else
      #ifdef LAPACK_ILP64
	#error "Invalid compile flags."
      #else
	#ifdef _WIN32
  // probably don't need this?
  typedef int32_t __CLPK_integer;
	#else
  typedef long int __CLPK_integer;
        #endif
      #endif
    #endif // !__LP64__

  void xerbla_(void);
    #ifdef __cplusplus
}
    #endif

  #endif // !__APPLE__

  typedef __CLPK_integer matrix_invert_buf1_t;
  // need to be careful about >= 2^32?
  CONSTU31(kMatrixInvertBuf1ElemAlloc, sizeof(__CLPK_integer));
  // invert_matrix_checked() usually requires a larger buffer
  CONSTU31(kMatrixInvertBuf1CheckedAlloc, 2 * sizeof(__CLPK_integer));
  typedef __CLPK_integer matrix_finvert_buf1_t;
  CONSTU31(kMatrixFinvertBuf1CheckedAlloc, 2 * sizeof(__CLPK_integer));

#endif // NOLAPACK

#ifdef __cplusplus
namespace plink2 {
#endif

static const double kMatrixSingularRcond = 1e-14;

#ifdef NOLAPACK
boolerr_t invert_matrix(int32_t dim, double* matrix, matrix_invert_buf1_t* dbl_1d_buf, double* dbl_2d_buf);

HEADER_INLINE boolerr_t invert_matrix_checked(int32_t dim, double* matrix, matrix_invert_buf1_t* dbl_1d_buf, double* dbl_2d_buf) {
  return invert_matrix(dim, matrix, dbl_1d_buf, dbl_2d_buf);
}

// if we're using float32s instead of float64s, we care enough about low-level
// details that this split interface makes sense.
// first half computes either LU or singular value decomposition, and
//   determinant (could make latter optional)
// second half actually inverts matrix, assuming 1d_buf and 2d_buf have results
//   from first half
boolerr_t invert_fmatrix_first_half(int32_t dim, int32_t stride, float* matrix, float* absdet_ptr, matrix_finvert_buf1_t* flt_1d_buf, float* flt_2d_buf);

void invert_fmatrix_second_half(int32_t dim, int32_t stride, float* matrix, matrix_finvert_buf1_t* flt_1d_buf, float* flt_2d_buf);
#else
boolerr_t invert_matrix(__CLPK_integer dim, double* matrix, matrix_invert_buf1_t* int_1d_buf, double* dbl_2d_buf);

boolerr_t invert_matrix_checked(__CLPK_integer dim, double* matrix, matrix_invert_buf1_t* int_1d_buf, double* dbl_2d_buf);

boolerr_t invert_fmatrix_first_half(__CLPK_integer dim, __CLPK_integer stride, float* matrix, float* absdet_ptr, matrix_finvert_buf1_t* int_1d_buf, float* flt_2d_buf);

void invert_fmatrix_second_half(__CLPK_integer dim, __CLPK_integer stride, float* matrix, matrix_finvert_buf1_t* int_1d_buf, float* flt_2d_buf);
#endif

void col_major_matrix_multiply(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* outmatrix);

HEADER_INLINE void row_major_matrix_multiply(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* outmatrix) {
  return col_major_matrix_multiply(inmatrix2, inmatrix1, col2_ct, row1_ct, common_ct, outmatrix);
}

// this is essentially a full-blown dgemm wrapper, only missing the alpha
// parameter now
void col_major_matrix_multiply_strided_addassign(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double beta, double* outmatrix);

HEADER_INLINE void row_major_matrix_multiply_incr(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* outmatrix) {
  return col_major_matrix_multiply_strided_addassign(inmatrix2, inmatrix1, col2_ct, col2_ct, row1_ct, common_ct, common_ct, col2_ct, 1.0, outmatrix);
}

HEADER_INLINE void row_major_matrix_multiply_strided(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double* outmatrix) {
  // stride1 should be close to common_ct
  // stride2 should be close to col2_ct
  // output matrix uses stride3, which should be close to col2_ct
  return col_major_matrix_multiply_strided_addassign(inmatrix2, inmatrix1, col2_ct, stride2, row1_ct, stride1, common_ct, stride3, 0.0, outmatrix);
}

HEADER_INLINE void row_major_matrix_multiply_strided_incr(const double* inmatrix1, const double* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, double* outmatrix) {
  return col_major_matrix_multiply_strided_addassign(inmatrix2, inmatrix1, col2_ct, stride2, row1_ct, stride1, common_ct, stride3, 1.0, outmatrix);
}

void col_major_fmatrix_multiply_strided(const float* inmatrix1, const float* inmatrix2, __CLPK_integer row1_ct, __CLPK_integer stride1, __CLPK_integer col2_ct, __CLPK_integer stride2, __CLPK_integer common_ct, __CLPK_integer stride3, float* outmatrix);

void transpose_copy(const double* old_matrix, uint32_t old_maj, uint32_t new_maj, double* new_matrix_iter);

void transpose_copy_float(const float* old_matrix, uint32_t old_maj, uint32_t new_maj, uint32_t new_maj_max, float* new_matrix_iter);


// A(A^T), where A is row-major; result is dim x dim
// ONLY UPDATES LOWER TRIANGLE OF result[].
void multiply_self_transpose(double* input_matrix, uint32_t dim, uint32_t col_ct, double* result);

// (A^T)A
void transpose_multiply_self_incr(double* input_part, uint32_t dim, uint32_t partial_row_ct, double* result);

#ifndef NOLAPACK
boolerr_t get_svd_rect_lwork(uint32_t major_ct, uint32_t minor_ct, __CLPK_integer* lwork_ptr);

// currently a wrapper for dgesvd_().
boolerr_t svd_rect(uint32_t major_ct, uint32_t minor_ct, __CLPK_integer lwork, double* matrix, double* ss, unsigned char* svd_rect_wkspace);

boolerr_t get_extract_eigvecs_lworks(uint32_t dim, uint32_t pc_ct, __CLPK_integer* lwork_ptr, __CLPK_integer* liwork_ptr, uintptr_t* wkspace_byte_ct_ptr);

// currently a wrapper for dsyevr_().
// reverse_eigvecs is eigenvector-major, but the vectors are in order of
// *increasing* eigenvalue.
boolerr_t extract_eigvecs(uint32_t dim, uint32_t pc_ct, __CLPK_integer lwork, __CLPK_integer liwork, double* matrix, double* eigvals, double* reverse_eigvecs, unsigned char* extract_eigvecs_wkspace);
#endif

boolerr_t linear_regression_inv(const double* pheno_d, double* predictors_pmaj, uint32_t predictor_ct, uint32_t sample_ct, double* fitted_coefs, double* xtx_inv, double* xt_y, matrix_invert_buf1_t* mi_buf, double* dbl_2d_buf);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_MATRIX_H__
