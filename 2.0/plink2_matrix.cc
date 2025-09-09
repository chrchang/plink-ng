// This file is part of PLINK 2.0, copyright (C) 2005-2022 Shaun Purcell,
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

#include "plink2_matrix.h"

#include <assert.h>
#include <math.h>  // fabs(), sqrt()
#include <string.h>

#ifndef NOLAPACK
#  if defined(__APPLE__) || defined(USE_MKL)
#    define LAPACK_dgecon dgecon_
#    define LAPACK_dgesvd dgesvd_
#    define LAPACK_dgetrf dgetrf_
#    define LAPACK_dgetri dgetri_
#    define LAPACK_dlange dlange_
#    define LAPACK_dlansy dlansy_
#    define LAPACK_dpocon dpocon_
#    define LAPACK_dpotrf dpotrf_
#    define LAPACK_dpotrs dpotrs_
#    define LAPACK_dsyevr dsyevr_
#  endif
#endif  // !NOLAPACK


#ifdef __cplusplus
namespace plink2 {
#endif

intptr_t FirstInfOrNan(const double* vec, uintptr_t size) {
  // possible todo: vectorize this.  Main loop could bitwise-or results from
  // several vector-comparisons at a time, use e.g. movemask to check for at
  // least one hit, then backtrack to identify first location.
  const uint64_t* vec_alias = R_CAST(const uint64_t*, vec);
  for (uintptr_t ulii = 0; ulii != size; ++ulii) {
    if ((vec_alias[ulii] & (0x7ffLLU << 52)) == (0x7ffLLU << 52)) {
      return S_CAST(intptr_t, ulii);
    }
  }
  return -1;
}

uint32_t LowerTriangularFirstInfOrNan(const double* matrix, uintptr_t dim, uintptr_t* row_idx_ptr, uintptr_t* col_idx_ptr) {
  const uint64_t* row_start = R_CAST(const uint64_t*, matrix);
  for (uintptr_t row_idx = 0; row_idx != dim; ++row_idx) {
    for (uintptr_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      if ((row_start[col_idx] & (0x7ffLLU << 52)) == (0x7ffLLU << 52)) {
        *row_idx_ptr = row_idx;
        *col_idx_ptr = col_idx;
        return 1;
      }
    }
    row_start = &(row_start[dim]);
  }
  return 0;
}

void ReflectMatrix(uint32_t dim, double* matrix) {
  const uintptr_t dim_p1l = dim + 1;
  double* write_row = matrix;
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    double* read_col_iter = &(matrix[dim_p1l * row_idx + dim]);
    for (uint32_t col_idx = row_idx + 1; col_idx != dim; ++col_idx) {
      write_row[col_idx] = *read_col_iter;
      read_col_iter = &(read_col_iter[dim]);
    }
    write_row = &(write_row[dim]);
  }
}

void ReflectStridedMatrix(uint32_t dim, uint32_t stride, double* matrix) {
  const uintptr_t stride_p1l = stride + 1;
  double* write_row = matrix;
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    double* read_col_iter = &(matrix[stride_p1l * row_idx + stride]);
    for (uint32_t col_idx = row_idx + 1; col_idx != dim; ++col_idx) {
      write_row[col_idx] = *read_col_iter;
      read_col_iter = &(read_col_iter[stride]);
    }
    write_row = &(write_row[stride]);
  }
}

void ReflectFmatrix(uint32_t dim, uint32_t stride, float* matrix) {
  const uintptr_t stride_p1l = stride + 1;
  float* write_row = matrix;
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    float* read_col_iter = &(matrix[stride_p1l * row_idx + stride]);
    for (uint32_t col_idx = row_idx + 1; col_idx != dim; ++col_idx) {
      write_row[col_idx] = *read_col_iter;
      read_col_iter = &(read_col_iter[stride]);
    }
    write_row = &(write_row[stride]);
  }
}

void ReflectStridedMatrix0(uint32_t dim, uint32_t stride, double* matrix) {
  const uintptr_t stride_p1l = stride + 1;
  double* write_row = matrix;
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    double* read_col_iter = &(matrix[stride_p1l * row_idx + stride]);
    for (uint32_t col_idx = row_idx + 1; col_idx != dim; ++col_idx) {
      write_row[col_idx] = *read_col_iter;
      read_col_iter = &(read_col_iter[stride]);
    }
    ZeroDArr(stride - dim, &(write_row[dim]));
    write_row = &(write_row[stride]);
  }
}

void ReflectFmatrix0(uint32_t dim, uint32_t stride, float* matrix) {
  const uintptr_t stride_p1l = stride + 1;
  float* write_row = matrix;
  for (uint32_t row_idx = 0; row_idx != dim; ++row_idx) {
    float* read_col_iter = &(matrix[stride_p1l * row_idx + stride]);
    for (uint32_t col_idx = row_idx + 1; col_idx != dim; ++col_idx) {
      write_row[col_idx] = *read_col_iter;
      read_col_iter = &(read_col_iter[stride]);
    }
    ZeroFArr(stride - dim, &(write_row[dim]));
    write_row = &(write_row[stride]);
  }
}

static inline double Sqr(const double aa) {
  return aa * aa;
}

double Pythag(const double a, const double b) {
  // PLINK stats.cpp Pythag().
  double absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
    return absa*sqrt(1.0+Sqr(absb/absa));
  }
  return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+Sqr(absa/absb)));
}

#ifdef NOLAPACK
static inline double MultiplyBySgn(double aa, double bb) {
  // PLINK helper.h MultiplyBySgn() template specialized to doubles.
  return (bb >= 0)? (aa >= 0 ? aa : -aa) : (aa >= 0 ? -aa : aa);
}

uint32_t SvdcmpC(int32_t m, double* a, double* w, double* v) {
  // C port of PLINK stats.cpp svdcmp().
  //   m: (in) number of rows/columns
  //   a: (in/out) input matrix, becomes first part of decomposition
  //   w: (out) singular values
  //   v: (out) final part of decomposition
  //
  // Now thread-safe.
  double* rv1 = &(w[S_CAST(uint32_t, m)]);
  int32_t n = m;
  int32_t flag;
  int32_t l = 0;  // suppress compile warning
  int32_t i,its,j,jj,k,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double temp;

  g=scale=anorm=0.0;
  for (i=0;i!=n;++i) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k!=m;++k) scale += fabs(a[k * m + i]);
      if (scale != 0.0) {
        for (k=i;k!=m;++k) {
          a[k * m + i] /= scale;
          s += a[k * m + i]*a[k * m + i];
        }
        f=a[i * m + i];
        g = -MultiplyBySgn(sqrt(s),f);
        h=f*g-s;
        a[i * m + i]=f-g;
        for (j=l-1;j!=n;++j) {
          for (s=0.0,k=i;k!=m;++k) s += a[k * m + i]*a[k * m + j];
          f=s/h;
          for (k=i;k!=m;++k) a[k * m + j] += f*a[k * m + i];
        }
        for (k=i;k!=m;++k) a[k * m + i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k!=n;++k) scale += fabs(a[i * m + k]);
      if (scale != 0.0) {
        for (k=l-1;k!=n;++k) {
          a[i * m + k] /= scale;
          s += a[i * m + k]*a[i * m + k];
        }
        f=a[i * m + l-1];
        g = -MultiplyBySgn(sqrt(s),f);
        h=f*g-s;
        a[i * m + l-1]=f-g;
        for (k=l-1;k!=n;++k) rv1[k]=a[i * m + k]/h;
        for (j=l-1;j!=m;++j) {
          for (s=0.0,k=l-1;k!=n;++k) s += a[j * m + k]*a[i * m + k];
          for (k=l-1;k!=n;++k) a[j * m + k] += s*rv1[k];
        }
        for (k=l-1;k!=n;++k) a[i * m + k] *= scale;
      }
    }
    anorm=MAXV(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
        for (j=l;j!=n;++j)
          v[j * m + i]=(a[i * m + j]/a[i * m + l])/g;
        for (j=l;j!=n;++j) {
          for (s=0.0,k=l;k!=n;++k) s += a[i * m + k]*v[k * m + j];
          for (k=l;k!=n;++k) v[k * m + j] += s*v[k * m + i];
        }
      }
      for (j=l;j!=n;++j) v[i * m + j]=v[j * m + i]=0.0;
    }
    v[i * m + i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MINV(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j!=n;++j) a[i * m + j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j!=n;++j) {
        for (s=0.0,k=l;k<m;++k) s += a[k * m + i]*a[k * m + j];
        f=(s/a[i * m + i])*g;
        for (k=i;k!=m;++k) a[k * m + j] += f*a[k * m + i];
      }
      for (j=i;j!=m;++j) a[j * m + i] *= g;
    } else for (j=i;j!=m;++j) a[j * m + i]=0.0;
    ++a[i * m + i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its!=30;++its) {
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
        for (i=l;i!=k+1;++i) {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          temp = fabs(f)+anorm;
          if (temp == anorm) break;
          g=w[i];
          h=Pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=0;j!=m;++j) {
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
          for (j=0;j!=n;++j) v[j * m + k] = -v[j * m + k];
        }
        break;
      }
      if (its == 29)
        return 0;  // cannot converge: multi-collinearity?
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=Pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+MultiplyBySgn(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=Pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=0;jj!=n;++jj) {
          x=v[jj * m + j];
          z=v[jj * m + i];
          v[jj * m + j]=x*c+z*s;
          v[jj * m + i]=z*c-x*s;
        }
        z=Pythag(f,h);
        w[j]=z;
        if (z != 0.0) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=0;jj!=m;++jj) {
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
  return 1;
}

BoolErr InvertMatrix(int32_t dim, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  // C port of PLINK stats.cpp's svd_inverse() function.
  assert(dim > 0);

  // w -> dbl_1d_buf
  // v -> dbl_2d_buf
  const double eps = 1e-24;
  // not unlikely(), matrix inversion failure can be common
  if (!SvdcmpC(dim, matrix, dbl_1d_buf, dbl_2d_buf)) {
    return 1;
  }

  // Look for singular values
  double wmax = 0;
  int32_t i;
  for (i=0; i!=dim; ++i) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i!=dim; ++i) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : (1 / dbl_1d_buf[i]);
  }

  int32_t j;
  for (i=0; i!=dim; ++i) {
    for (j=0; j!=dim; ++j) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  int32_t k;
  // [nxn].[t(v)]
  for (i=0; i!=dim; ++i) {
    ZeroDArr(dim, dbl_1d_buf);
    for (j=0; j!=dim; ++j) {
      for (k=0; k!=dim; ++k) {
        dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j != dim; ++j) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  for (i=1; i!=dim; ++i) {
    for(j=0; j!=i; ++j) {
      const double tmp = matrix[i * dim + j];
      matrix[i * dim + j] = matrix[j * dim + i];
      matrix[j * dim + i] = tmp;
    }
  }
  return 0;
}

BoolErr InvertStridedMatrixFirstHalf(int32_t dim, int32_t stride, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  // Unfortunately, SvdcmpC() doesn't support a stride parameter.
  if (stride != dim) {
    const uint32_t dim_u = dim;
    const uintptr_t nbyte = dim_u * sizeof(double);
    const double* read_row = matrix;
    double* write_row = matrix;
    for (uint32_t row_idx = 1; row_idx != dim_u; ++row_idx) {
      read_row = &(read_row[S_CAST(uint32_t, stride)]);
      write_row = &(write_row[dim_u]);
      memmove(write_row, read_row, nbyte);
    }
  }

  return (!SvdcmpC(dim, matrix, dbl_1d_buf, dbl_2d_buf));
}

BoolErr InvertFmatrixFirstHalf(int32_t dim, uint32_t stride, const float* matrix, double* half_inverted, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  const float* read_row = matrix;
  double* write_row = half_inverted;
  for (uint32_t row_idx = 0; row_idx != S_CAST(uint32_t, dim); ++row_idx) {
    for (uint32_t col_idx = 0; col_idx != S_CAST(uint32_t, dim); ++col_idx) {
      write_row[col_idx] = S_CAST(double, read_row[col_idx]);
    }
    read_row = &(read_row[stride]);
    write_row = &(write_row[S_CAST(uint32_t, dim)]);
  }

  return (!SvdcmpC(dim, half_inverted, dbl_1d_buf, dbl_2d_buf));
}

void InvertStridedMatrixSecondHalf(lapack_int dim, lapack_int stride, double* matrix, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  // Look for singular values
  assert(dim > 0);
  const double eps = 1e-24;
  double wmax = 0;
  int32_t i;
  for (i=0; i!=dim; ++i) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i!=dim; ++i) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : (1 / dbl_1d_buf[i]);
  }

  int32_t j;
  for (i=0; i!=dim; ++i) {
    for (j=0; j!=dim; ++j) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  int32_t k;
  // [nxn].[t(v)]
  for (i=0; i!=dim; ++i) {
    ZeroDArr(dim, dbl_1d_buf);
    for (j=0; j!=dim; ++j) {
      for (k=0; k!=dim; ++k) {
        dbl_1d_buf[j] += matrix[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j != dim; ++j) {
      matrix[i * dim + j] = dbl_1d_buf[j];
    }
  }
  // transpose, then fix stride
  for (i=1; i!=dim; ++i) {
    for(j=0; j!=i; ++j) {
      const double tmp = matrix[i * dim + j];
      matrix[i * dim + j] = matrix[j * dim + i];
      matrix[j * dim + i] = tmp;
    }
  }
  if ((stride != dim) && (dim > 1)) {
    const uint32_t dim_m1 = dim - 1;
    const uintptr_t nbyte = S_CAST(uint32_t, dim) * sizeof(double);
    double* read_row = &(matrix[dim_m1 * dim]);
    double* write_row = &(matrix[dim_m1 * stride]);
    for (uint32_t bot_row_idx = 0; bot_row_idx != dim_m1; ++bot_row_idx) {
      memmove(write_row, read_row, nbyte);
      read_row -= dim;
      write_row -= stride;
    }
  }
}

void InvertFmatrixSecondHalf(lapack_int dim, uint32_t stride, double* half_inverted, float* inverted_result, MatrixInvertBuf1* dbl_1d_buf, double* dbl_2d_buf) {
  // Look for singular values
  assert(dim > 0);
  const double eps = 1e-24;
  double wmax = 0;
  int32_t i;
  for (i=0; i!=dim; ++i) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i!=dim; ++i) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : (1 / dbl_1d_buf[i]);
  }

  int32_t j;
  for (i=0; i!=dim; ++i) {
    for (j=0; j!=dim; ++j) {
      half_inverted[i * dim + j] = half_inverted[i * dim + j] * dbl_1d_buf[j];
    }
  }

  int32_t k;
  // [nxn].[t(v)]
  for (i=0; i!=dim; ++i) {
    ZeroDArr(dim, dbl_1d_buf);
    for (j=0; j!=dim; ++j) {
      for (k=0; k!=dim; ++k) {
        dbl_1d_buf[j] += half_inverted[i * dim + k] * dbl_2d_buf[j * dim + k];
      }
    }
    for (j = 0; j != dim; ++j) {
      half_inverted[i * dim + j] = dbl_1d_buf[j];
    }
  }
  inverted_result[0] = S_CAST(float, half_inverted[0]);
  for (i=1; i!=dim; ++i) {
    for(j=0; j!=i; ++j) {
      inverted_result[i * stride + j] = S_CAST(float, half_inverted[j * dim + i]);
      inverted_result[j * stride + i] = S_CAST(float, half_inverted[i * dim + j]);
    }
    inverted_result[i * stride + i] = S_CAST(float, half_inverted[i * dim + i]);
  }
}
#else  // !NOLAPACK
BoolErr InvertMatrix(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  // InvertSymmdefMatrix() is noticeably faster in the symmetric
  // positive-semidefinite case.
  lapack_int info;
  LAPACK_dgetrf(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  if (info) {
    return 1;
  }
  lapack_int lwork = dim * dim;
  LAPACK_dgetri(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  assert(info == 0);
  return 0;
}

BoolErr InvertMatrixChecked(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  // This used to fall back on PLINK 1.07's SVD-based implementation when the
  // rcond estimate was too small, but in practice that just slowed things down
  // without meaningfully improving inversion of nonsingular matrices.  So now
  // this just exits a bit earlier, while leaving the old "binary search for
  // the first row/column causing multicollinearity" logic to the caller.
  char cc = '1';
  double norm = LAPACK_dlange(&cc, &dim, &dim, matrix, &dim, dbl_2d_buf);
  lapack_int info;
  LAPACK_dgetrf(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dgecon(&cc, &dim, matrix, &dim, &norm, &rcond, dbl_2d_buf, &(int_1d_buf[S_CAST(uint32_t, dim)]), &info);
  if (rcond < kMatrixSingularRcond) {
    return 1;
  }
  lapack_int lwork = dim * dim;
  LAPACK_dgetri(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  return 0;
}

BoolErr InvertSymmdefMatrix(lapack_int dim, double* matrix, __maybe_unused MatrixInvertBuf1* int_1d_buf, __maybe_unused double* dbl_2d_buf) {
  char uplo = 'U';
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &dim, &info);
  if (info) {
    return 1;
  }
  LAPACK_dpotri(&uplo, &dim, matrix, &dim, &info);
  assert(info == 0);
  return 0;
}

BoolErr InvertSymmdefStridedMatrix(lapack_int dim, lapack_int stride, double* matrix, __maybe_unused MatrixInvertBuf1* int_1d_buf, __maybe_unused double* dbl_2d_buf) {
  char uplo = 'U';
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &stride, &info);
  if (info) {
    return 1;
  }
  LAPACK_dpotri(&uplo, &dim, matrix, &stride, &info);
  assert(info == 0);
  return 0;
}

BoolErr InvertSymmdefMatrixChecked(lapack_int dim, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  char cc = '1';
  char uplo = 'U';
  double norm = LAPACK_dlansy(&cc, &uplo, &dim, matrix, &dim, dbl_2d_buf);
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &dim, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dpocon(&uplo, &dim, matrix, &dim, &norm, &rcond, dbl_2d_buf, int_1d_buf, &info);
  if (rcond < kMatrixSingularRcond) {
    return 1;
  }
  LAPACK_dpotri(&uplo, &dim, matrix, &dim, &info);
  return 0;
}

BoolErr InvertSymmdefMatrixFirstHalf(lapack_int dim, lapack_int stride, double* matrix, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  char cc = '1';
  char uplo = 'U';
  double norm = LAPACK_dlansy(&cc, &uplo, &dim, matrix, &stride, dbl_2d_buf);
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, matrix, &stride, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dpocon(&uplo, &dim, matrix, &stride, &norm, &rcond, dbl_2d_buf, int_1d_buf, &info);
  return (rcond < kMatrixSingularRcond);
}

// quasi-bugfix (20 Sep 2017): give up on doing this with single-precision
// numbers.  Instead, convert to double-precision, then perform the usual
// inversion, then downcode back to single-precision.
BoolErr InvertFmatrixFirstHalf(lapack_int dim, uint32_t stride, const float* matrix, double* half_inverted, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  const float* read_row = matrix;
  double* write_row = half_inverted;
  for (uint32_t row_idx = 0; row_idx != S_CAST(uint32_t, dim); ++row_idx) {
    // could use _mm256_cvtps_pd() here
    for (uint32_t col_idx = 0; col_idx != S_CAST(uint32_t, dim); ++col_idx) {
      write_row[col_idx] = S_CAST(double, read_row[col_idx]);
    }
    read_row = &(read_row[stride]);
    write_row = &(write_row[S_CAST(uint32_t, dim)]);
  }

  char cc = '1';
  double norm = LAPACK_dlange(&cc, &dim, &dim, half_inverted, &dim, dbl_2d_buf);
  lapack_int info;
  LAPACK_dgetrf(&dim, &dim, half_inverted, &dim, int_1d_buf, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dgecon(&cc, &dim, half_inverted, &dim, &norm, &rcond, dbl_2d_buf, &(int_1d_buf[S_CAST(uint32_t, dim)]), &info);
  return (rcond < kMatrixSingularRcond);
}

BoolErr InvertSymmdefFmatrixFirstHalf(lapack_int dim, uint32_t stride, float* matrix, double* half_inverted, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  const float* read_row = matrix;
  double* write_row = half_inverted;
  for (uint32_t row_idx = 0; row_idx != S_CAST(uint32_t, dim); ++row_idx) {
    // could use _mm256_cvtps_pd() here
    for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      write_row[col_idx] = S_CAST(double, read_row[col_idx]);
    }
    read_row = &(read_row[stride]);
    write_row = &(write_row[S_CAST(uint32_t, dim)]);
  }

  char cc = '1';
  char uplo = 'U';
  double norm = LAPACK_dlansy(&cc, &uplo, &dim, half_inverted, &dim, dbl_2d_buf);
  lapack_int info;
  LAPACK_dpotrf(&uplo, &dim, half_inverted, &dim, &info);
  if (info > 0) {
    return 1;
  }
  double rcond;
  LAPACK_dpocon(&uplo, &dim, half_inverted, &dim, &norm, &rcond, dbl_2d_buf, int_1d_buf, &info);
  return (rcond < kMatrixSingularRcond);
}

void InvertFmatrixSecondHalf(lapack_int dim, uint32_t stride, double* half_inverted, float* inverted_result, MatrixInvertBuf1* int_1d_buf, double* dbl_2d_buf) {
  lapack_int lwork = dim * dim;
  lapack_int info;
  LAPACK_dgetri(&dim, half_inverted, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  const double* read_row = half_inverted;
  float* write_row = inverted_result;
  for (uint32_t row_idx = 0; row_idx != S_CAST(uint32_t, dim); ++row_idx) {
    // could use _mm256_cvtpd_ps() here
    for (uint32_t col_idx = 0; col_idx != S_CAST(uint32_t, dim); ++col_idx) {
      write_row[col_idx] = S_CAST(float, read_row[col_idx]);
    }
    read_row = &(read_row[S_CAST(uint32_t, dim)]);
    write_row = &(write_row[stride]);
  }
}

void InvertSymmdefFmatrixSecondHalf(lapack_int dim, uint32_t stride, double* half_inverted, float* inverted_result, __maybe_unused MatrixInvertBuf1* int_1d_buf, __maybe_unused double* dbl_2d_buf) {
  char uplo = 'U';
  lapack_int info;
  LAPACK_dpotri(&uplo, &dim, half_inverted, &dim, &info);
  const double* read_row = half_inverted;
  float* write_row = inverted_result;
  for (uint32_t row_idx = 0; row_idx != S_CAST(uint32_t, dim); ++row_idx) {
    // could use _mm256_cvtpd_ps() here
    for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
      write_row[col_idx] = S_CAST(float, read_row[col_idx]);
    }
    read_row = &(read_row[S_CAST(uint32_t, dim)]);
    write_row = &(write_row[stride]);
  }
}
#endif  // !NOLAPACK

void ColMajorMatrixMultiply(const double* inmatrix1, const double* inmatrix2, lapack_int row1_ct, lapack_int col2_ct, lapack_int common_ct, double* outmatrix) {
#ifdef NOLAPACK
  const uintptr_t row1_ct_l = row1_ct;
  const uintptr_t col2_ct_l = col2_ct;
  const uintptr_t common_ct_l = common_ct;
  // not optimized
  for (uintptr_t col_idx = 0; col_idx != col2_ct_l; ++col_idx) {
    for (uintptr_t row_idx = 0; row_idx != row1_ct_l; ++row_idx) {
      double cur_dotprod = 0.0;
      const double* dptr = &(inmatrix2[col_idx * common_ct]);
      for (uintptr_t com_idx = 0; com_idx != common_ct_l; ++com_idx) {
        cur_dotprod += (*dptr++) * inmatrix1[com_idx * row1_ct_l + row_idx];
      }
      *outmatrix++ = cur_dotprod;
    }
  }
#else
  // bugfix (30 Aug 2017): this fails on OS X when LDB > sqrt(2^31).
  // update: Windows does not have the same problem
  // update 2 (6 Sep 2017): the OS X failure seems to have been driven by 128k
  //   thread stack size; going up to the usual 512k appears to solve the
  //   problem
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row1_ct, col2_ct, common_ct, 1.0, inmatrix1, row1_ct, inmatrix2, common_ct, 0.0, outmatrix, row1_ct);
#endif  // !NOLAPACK
}

void ColMajorMatrixMultiplyStridedAddassign(const double* inmatrix1, const double* inmatrix2, lapack_int row1_ct, lapack_int stride1, lapack_int col2_ct, lapack_int stride2, lapack_int common_ct, lapack_int stride3, double beta, double* outmatrix) {
  // stride1 should be close to row1_ct
  // stride2 should be close to common_ct
  // output matrix uses stride3, which should be close to row1_ct
#ifdef NOLAPACK
  const uintptr_t row1_ct_l = row1_ct;
  const uintptr_t col2_ct_l = col2_ct;
  const uintptr_t common_ct_l = common_ct;
  // not optimized, no beta == 0 special case
  for (uintptr_t col_idx = 0; col_idx != col2_ct_l; ++col_idx) {
    double* outmatrix_row_iter = &(outmatrix[col_idx * stride3]);
    for (uintptr_t row_idx = 0; row_idx != row1_ct_l; ++row_idx) {
      double cur_entry = 0.0;
      const double* col2_iter = &(inmatrix2[col_idx * stride2]);
      for (uintptr_t com_idx = 0; com_idx != common_ct_l; com_idx++) {
        cur_entry += (*col2_iter++) * inmatrix1[com_idx * stride1 + row_idx];
      }
      *outmatrix_row_iter = (*outmatrix_row_iter) * beta + cur_entry;
      ++outmatrix_row_iter;
    }
  }
#else
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row1_ct, col2_ct, common_ct, 1.0, inmatrix1, stride1, inmatrix2, stride2, beta, outmatrix, stride3);
#endif  // !NOLAPACK
}

void ColMajorVectorMatrixMultiplyStrided(const double* in_dvec1, const double* inmatrix2, lapack_int common_ct, lapack_int stride2, lapack_int col2_ct, double* out_dvec) {
#ifdef NOLAPACK
  const uintptr_t col2_ct_l = col2_ct;
  const uintptr_t common_ct_l = common_ct;
  const uintptr_t stride_l = stride2;
  for (uintptr_t col_idx = 0; col_idx != col2_ct_l; ++col_idx) {
    double dxx = 0.0;
    const double* cur_col = &(inmatrix2[col_idx * stride_l]);
    for (uintptr_t row_idx = 0; row_idx != common_ct_l; ++row_idx) {
      dxx += in_dvec1[row_idx] * cur_col[row_idx];
    }
    out_dvec[col_idx] = dxx;
  }
#else
  cblas_dgemv(CblasColMajor, CblasTrans, common_ct, col2_ct, 1.0, inmatrix2, stride2, in_dvec1, 1, 0.0, out_dvec, 1);
#endif  // !NOLAPACK
}

void ColMajorMatrixVectorMultiplyStrided(const double* inmatrix1, const double* in_dvec2, lapack_int row1_ct, lapack_int stride1, lapack_int common_ct, double* out_dvec) {
#ifdef NOLAPACK
  ZeroDArr(row1_ct, out_dvec);
  const uintptr_t row1_ct_l = row1_ct;
  const uintptr_t common_ct_l = common_ct;
  const uintptr_t stride_l = stride1;
  const double* col_iter = inmatrix1;
  for (uintptr_t common_idx = 0; common_idx != common_ct_l; ++common_idx) {
    const double dxx = in_dvec2[common_idx];
    for (uintptr_t row1_idx = 0; row1_idx != row1_ct_l; ++row1_idx) {
      out_dvec[row1_idx] += col_iter[row1_idx] * dxx;
    }
    col_iter = &(col_iter[stride_l]);
  }
#else
  cblas_dgemv(CblasColMajor, CblasNoTrans, row1_ct, common_ct, 1.0, inmatrix1, stride1, in_dvec2, 1, 0.0, out_dvec, 1);
#endif  // !NOLAPACK
}

// er, should make this _addassign for consistency...
void ColMajorFmatrixMultiplyStrided(const float* inmatrix1, const float* inmatrix2, lapack_int row1_ct, lapack_int stride1, lapack_int col2_ct, lapack_int stride2, lapack_int common_ct, lapack_int stride3, float* outmatrix) {
#ifdef NOLAPACK
  const uintptr_t row1_ct_l = row1_ct;
  const uintptr_t col2_ct_l = col2_ct;
  const uintptr_t common_ct_l = common_ct;
  // not optimized
  for (uintptr_t col_idx = 0; col_idx != col2_ct_l; ++col_idx) {
    float* outmatrix_row_iter = &(outmatrix[col_idx * stride3]);
    for (uintptr_t row_idx = 0; row_idx != row1_ct_l; ++row_idx) {
      float cur_entry = 0.0;
      const float* col2_iter = &(inmatrix2[col_idx * stride2]);
      for (uintptr_t com_idx = 0; com_idx != common_ct_l; com_idx++) {
        cur_entry += (*col2_iter++) * inmatrix1[com_idx * stride1 + row_idx];
      }
      *outmatrix_row_iter++ = cur_entry;
    }
  }
#else
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row1_ct, col2_ct, common_ct, 1.0, inmatrix1, stride1, inmatrix2, stride2, 0.0, outmatrix, stride3);
#endif  // !NOLAPACK
}

void ColMajorFmatrixVectorMultiplyStrided(const float* inmatrix1, const float* in_fvec2, lapack_int row1_ct, lapack_int stride1, lapack_int common_ct, float* out_fvec) {
#ifdef NOLAPACK
  ZeroFArr(row1_ct, out_fvec);
  const uintptr_t row1_ct_l = row1_ct;
  const uintptr_t common_ct_l = common_ct;
  const uintptr_t stride_l = stride1;
  for (uintptr_t common_idx = 0; common_idx != common_ct_l; ++common_idx) {
    const float fxx = in_fvec2[common_idx];
    const float* col_iter = &(inmatrix1[common_idx * stride_l]);
    for (uintptr_t row1_idx = 0; row1_idx != row1_ct_l; ++row1_idx) {
      out_fvec[row1_idx] += (*col_iter++) * fxx;
    }
  }
#else
  cblas_sgemv(CblasColMajor, CblasNoTrans, row1_ct, common_ct, 1.0, inmatrix1, stride1, in_fvec2, 1, 0.0, out_fvec, 1);
#endif  // !NOLAPACK
}

void ColMajorFvectorMatrixMultiplyStrided(const float* in_fvec1, const float* inmatrix2, lapack_int common_ct, lapack_int stride2, lapack_int col2_ct, float* out_fvec) {
#ifdef NOLAPACK
  const uintptr_t col2_ct_l = col2_ct;
  const uintptr_t common_ct_l = common_ct;
  const uintptr_t stride_l = stride2;
  for (uintptr_t col_idx = 0; col_idx != col2_ct_l; ++col_idx) {
    float fxx = 0.0;
    const float* cur_col = &(inmatrix2[col_idx * stride_l]);
    for (uintptr_t row_idx = 0; row_idx != common_ct_l; ++row_idx) {
      fxx += in_fvec1[row_idx] * cur_col[row_idx];
    }
    out_fvec[col_idx] = fxx;
  }
#else
  cblas_sgemv(CblasColMajor, CblasTrans, common_ct, col2_ct, 1.0, inmatrix2, stride2, in_fvec1, 1, 0.0, out_fvec, 1);
#endif  // !NOLAPACK
}

// Briefly experimented with trying to speed this up, didn't make any progress.
// Tried vDSP_mtransD() on macOS and mkl_domatcopy() on Linux, neither of those
// helped, so it's likely there's no progress to be made.
void MatrixTransposeCopy(const double* old_matrix, uint32_t old_maj, uint32_t new_maj, double* new_matrix_iter) {
  for (uint32_t new_maj_idx = 0; new_maj_idx != new_maj; ++new_maj_idx) {
    const double* old_matrix_col_iter = &(old_matrix[new_maj_idx]);
    for (uint32_t old_maj_idx = 0; old_maj_idx != old_maj; ++old_maj_idx) {
      *new_matrix_iter++ = *old_matrix_col_iter;
      old_matrix_col_iter = &(old_matrix_col_iter[new_maj]);
    }
  }
}

void FmatrixTransposeCopy(const float* old_matrix, uint32_t old_maj, uint32_t new_maj, uint32_t new_maj_max, float* new_matrix_iter) {
  // new_maj = in-memory stride of old_matrix rows
  // new_maj_max = actual number of rows in new_matrix
  // (distinction is necessary for SSE alignment)
  for (uint32_t new_maj_idx = 0; new_maj_idx != new_maj_max; ++new_maj_idx) {
    const float* old_matrix_col_iter = &(old_matrix[new_maj_idx]);
    for (uint32_t old_maj_idx = 0; old_maj_idx != old_maj; ++old_maj_idx) {
      *new_matrix_iter++ = *old_matrix_col_iter;
      old_matrix_col_iter = &(old_matrix_col_iter[new_maj]);
    }
  }
}

/*
lapack_int qr_square_factor_float_get_lwork(uint32_t dim) {
  lapack_int dim_i = (lapack_int)dim;
  lapack_int sgeqrf_lwork = -1;
  float work[1];
  lapack_int info;
  float dummy;
  float dummy2;
  sgeqrf_(&dim_i, &dim_i, &dummy, &dim_i, &dummy2, work, &sgeqrf_lwork, &info);
  assert(info == 0);
  sgeqrf_lwork = (lapack_int)work[0];
  lapack_int sorgqr_lwork = -1;
  sorgqr_(&dim_i, &dim_i, &dim_i, &dummy, &dim_i, &dummy2, work, &sorgqr_lwork, &info);
  assert(info == 0);
  sorgqr_lwork = (lapack_int)work[0];
  return MAXV(sgeqrf_lwork, sorgqr_lwork);
}

BoolErr qr_square_factor_float(const float* input_matrix, uint32_t dim, uintptr_t stride, lapack_int lwork, float* qq, float* r_determinant_ptr, float* tau_buf, float* work_buf) {
  // only returns Q and, optionally, the product of R's diagonal entries (which
  // should be the determinant of the original matrix).
  // tau_buf should have space for dim entries
  assert(dim > 0);
  if (dim == stride) {
    memcpy(qq, input_matrix, dim * ((uintptr_t)dim) * sizeof(float));
  } else {
    for (uintptr_t col_idx = 0; col_idx != dim; ++col_idx) {
      memcpy(&(qq[col_idx * dim]), &(input_matrix[col_idx * stride]), dim * sizeof(float));
    }
  }
  lapack_int dim_i = (lapack_int)dim;
  lapack_int info;
  sgeqrf_(&dim_i, &dim_i, qq, &dim_i, tau_buf, work_buf, &lwork, &info);
  if (info != 0) {
    return 1;
  }
  if (r_determinant_ptr) {
    const uintptr_t dimp1 = dim + 1;
    float prod = qq[0];
    for (uintptr_t col_idx = 1; col_idx != dim; ++col_idx) {
      prod *= qq[col_idx * dimp1];
    }
    *r_determinant_ptr = prod;
  }
  sorgqr_(&dim_i, &dim_i, &dim_i, qq, &dim_i, tau_buf, work_buf, &lwork, &info);
  if (info != 0) {
    return 1;
  }
  return 0;
}
*/

// A(A^T), where A is row-major; result is dim x dim
// ONLY UPDATES LOWER TRIANGLE OF result[].
void MultiplySelfTranspose(const double* input_matrix, uint32_t dim, uint32_t col_ct, double* result) {
#ifdef NOLAPACK
  for (uintptr_t row1_idx = 0; row1_idx != dim; ++row1_idx) {
    const double* pred_row1 = &(input_matrix[row1_idx * col_ct]);
    double* result_row = &(result[row1_idx * dim]);
    for (uintptr_t row2_idx = 0; row2_idx <= row1_idx; ++row2_idx) {
      const double* pred_row2 = &(input_matrix[row2_idx * col_ct]);
      result_row[row2_idx] = DotprodD(pred_row1, pred_row2, col_ct);
    }
  }
#else
  // see ColMajorMatrixMultiply() remarks; same OS X issue here.
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, dim, col_ct, 1.0, input_matrix, col_ct, 0.0, result, dim);
#endif
}

void MultiplySelfTransposeStrided(const double* input_matrix, uint32_t dim, uint32_t col_ct, uint32_t stride, double* result) {
#ifdef NOLAPACK
  for (uintptr_t row1_idx = 0; row1_idx != dim; ++row1_idx) {
    const double* pred_row1 = &(input_matrix[row1_idx * stride]);
    double* result_row = &(result[row1_idx * dim]);
    for (uintptr_t row2_idx = 0; row2_idx <= row1_idx; ++row2_idx) {
      const double* pred_row2 = &(input_matrix[row2_idx * stride]);
      result_row[row2_idx] = DotprodD(pred_row1, pred_row2, col_ct);
    }
  }
#else
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, dim, col_ct, 1.0, input_matrix, stride, 0.0, result, dim);
#endif
}

void MultiplySelfTransposeStridedF(const float* input_matrix, uint32_t dim, uint32_t col_ct, uint32_t stride, float* result) {
#ifdef NOLAPACK
  for (uintptr_t row1_idx = 0; row1_idx != dim; ++row1_idx) {
    const float* pred_row1 = &(input_matrix[row1_idx * stride]);
    float* result_row = &(result[row1_idx * dim]);
    for (uintptr_t row2_idx = 0; row2_idx <= row1_idx; ++row2_idx) {
      const float* pred_row2 = &(input_matrix[row2_idx * stride]);
      result_row[row2_idx] = DotprodF(pred_row1, pred_row2, col_ct);
    }
  }
#else
  cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans, dim, col_ct, 1.0, input_matrix, stride, 0.0, result, dim);
#endif
}

void TransposeMultiplySelfIncr(double* input_part, uint32_t dim, uint32_t partial_row_ct, double* result) {
#ifdef NOLAPACK
  // friends do not let friends use this implementation
  const uintptr_t dim_l = dim;
  const uintptr_t row_ct_l = partial_row_ct;
  for (uintptr_t idx1 = 0; idx1 != dim_l; ++idx1) {
    const double* col1 = &(input_part[idx1]);
    double* write_iter = &(result[idx1 * dim_l]);
    for (uintptr_t idx2 = 0; idx2 <= idx1; ++idx2) {
      double cur_dotprod = *write_iter;
      const double* col2 = &(input_part[idx2]);
      for (uintptr_t row_idx = 0; row_idx != row_ct_l; ++row_idx) {
        cur_dotprod += col1[row_idx * dim_l] * col2[row_idx * dim_l];
      }
      *write_iter = cur_dotprod;
      ++write_iter;
    }
  }
#else
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, dim, partial_row_ct, 1.0, input_part, dim, 1.0, result, dim);
#endif  // !NOLAPACK
}

#ifndef NOLAPACK
BoolErr GetSvdRectLwork(uint32_t major_ct, uint32_t minor_ct, lapack_int* lwork_ptr) {
  char jobu = 'S';
  char jobvt = 'O';
  lapack_int tmp_m = minor_ct;
  lapack_int tmp_n = major_ct;
  lapack_int wkspace_size = -1;
  double wkspace_size_d;
  lapack_int info;
  LAPACK_dgesvd(&jobu, &jobvt, &tmp_m, &tmp_n, nullptr, &tmp_m, nullptr, nullptr, &tmp_m, nullptr, &tmp_m, &wkspace_size_d, &wkspace_size, &info);
#  ifdef LAPACK_ILP64
  if (unlikely(info)) {
    return 1;
  }
#  else
  if (unlikely(info || (wkspace_size_d > 2147483640.0))) {
    return 1;
  }
#  endif
  *lwork_ptr = RoundUpPow2(S_CAST(lapack_int, wkspace_size_d), kCacheline / sizeof(double));
  return 0;
}

IntErr SvdRect(uint32_t major_ct, uint32_t minor_ct, lapack_int lwork, double* matrix, double* ss, double* vv, unsigned char* svd_rect_wkspace) {
  double* work = R_CAST(double*, svd_rect_wkspace);
  char jobu = 'S';
  char jobvt = 'O';
  lapack_int tmp_m = minor_ct;
  lapack_int tmp_n = major_ct;
  lapack_int info;
  LAPACK_dgesvd(&jobu, &jobvt, &tmp_m, &tmp_n, matrix, &tmp_m, ss, vv, &tmp_m, nullptr, &tmp_m, work, &lwork, &info);
  return S_CAST(IntErr, info);
}

// dsyevr_ takes ~30% less time than dsyevd_ on OS X dev machine.  todo: retest
// for Linux 64-bit MKL.
BoolErr GetExtractEigvecsLworks(uint32_t dim, uint32_t pc_ct, lapack_int* lwork_ptr, lapack_int* liwork_ptr, uintptr_t* wkspace_byte_ct_ptr) {
#ifndef LAPACK_ILP64
  if (unlikely(dim > 46340)) {
    // (30 Dec 2024) maybe this still works on some systems despite matrix size
    // being unrepresentable by lapack_int?  but it doesn't on macOS:
    //   https://groups.google.com/g/plink2-users/c/oteXlRFHdgk
    return 1;
  }
#endif
  char jobz = 'V';
  char range = 'I';
  char uplo = 'U';
  lapack_int tmp_n = dim;
  lapack_int il = dim + 1 - pc_ct;
  lapack_int iu = dim;
  double abstol = -1.0;
  lapack_int lwork_dummy = -1;

  // defined to make macOS 13.3 happy.
  double ignored_vl = 0.0;  // gcc 14 warning
  double ignored_vu = 0.0;  // gcc 14 warning
  lapack_int ignored_m;

  double lwork_d;
  lapack_int liwork;
  lapack_int info;
  LAPACK_dsyevr(&jobz, &range, &uplo, &tmp_n, nullptr, &tmp_n, &ignored_vl, &ignored_vu, &il, &iu, &abstol, &ignored_m, nullptr, nullptr, &tmp_n, nullptr, &lwork_d, &lwork_dummy, &liwork, &lwork_dummy, &info);
#ifdef LAPACK_ILP64
  if (unlikely(info)) {
    return 1;
  }
#else
  if (unlikely(info || (lwork_d > 2147483648.0 - (kCacheline / sizeof(double))))) {
    return 1;
  }
#endif
  const lapack_int lwork = RoundUpPow2(S_CAST(lapack_int, lwork_d), kCacheline / sizeof(double));
  liwork = RoundUpPow2(liwork, kCacheline / sizeof(lapack_int));
  *lwork_ptr = lwork;
  *liwork_ptr = liwork;
  *wkspace_byte_ct_ptr = lwork * sizeof(double) + liwork * sizeof(lapack_int) + RoundUpPow2(2 * dim * sizeof(lapack_int), kCacheline);
  return 0;
}

BoolErr ExtractEigvecs(uint32_t dim, uint32_t pc_ct, lapack_int lwork, lapack_int liwork, double* matrix, double* eigvals, double* reverse_eigvecs, unsigned char* extract_eigvecs_wkspace) {
  char jobz = 'V';
  char range = 'I';
  char uplo = 'U';
  lapack_int tmp_n = dim;
  lapack_int il = dim + 1 - pc_ct;
  lapack_int iu = dim;
  double abstol = -1.0;
  lapack_int out_m;
  double* work = R_CAST(double*, extract_eigvecs_wkspace);
  lapack_int* iwork = R_CAST(lapack_int*, &(work[lwork]));
  lapack_int* isuppz = &(iwork[liwork]);
  lapack_int info;
  // vl and vu may actually be referenced in some implementations
  double dummy_d = 0.0;
  LAPACK_dsyevr(&jobz, &range, &uplo, &tmp_n, matrix, &tmp_n, &dummy_d, &dummy_d, &il, &iu, &abstol, &out_m, eigvals, reverse_eigvecs, &tmp_n, isuppz, work, &lwork, iwork, &liwork, &info);
  // bug workaround (11 Feb 2023): on macOS, when there are NaNs in the matrix,
  // info=0 is still returned but out_m=0.
  // More generally, with the current interface, we want to return an error
  // whenever the number of found eigenvalues is less than pc_ct.
  return (S_CAST(uint32_t, out_m) != pc_ct) || (info != 0);
}
#endif  // !NOLAPACK

BoolErr invert_rank1_symm_start(const double* a_inv, const double* bb, lapack_int orig_dim, double cc, double* __restrict ainv_b, double* k_recip_ptr) {
#ifdef NOLAPACK
  const uintptr_t orig_dim_l = orig_dim;
  const double* a_inv_iter = a_inv;
  for (uintptr_t ulii = 0; ulii != orig_dim_l; ++ulii) {
    ainv_b[ulii] = DotprodD(bb, a_inv_iter, orig_dim_l);
    a_inv_iter = &(a_inv_iter[orig_dim_l]);
  }
#else
  cblas_dgemv(CblasColMajor, CblasNoTrans, orig_dim, orig_dim, 1.0, a_inv, orig_dim, bb, 1, 0.0, ainv_b, 1);
#endif
  const double kk = cc - DotprodxD(bb, ainv_b, orig_dim);
  if (fabs(kk) < kMatrixSingularRcond) {
    return 1;
  }
  *k_recip_ptr = 1.0 / kk;
  return 0;
}

BoolErr InvertRank1Symm(const double* a_inv, const double* bb, lapack_int orig_dim, uint32_t insert_idx, double cc, double* __restrict outmatrix, double* __restrict ainv_b_buf) {
  double k_recip;
  if (invert_rank1_symm_start(a_inv, bb, orig_dim, cc, ainv_b_buf, &k_recip)) {
    return 1;
  }
  // [ A^{-1} + k_recip * outer_prod(ainv_b, ainv_b)   -k_recip * ainv_b ]
  // [         -k_recip * ainv_b                            k_recip      ]
  const uintptr_t orig_dim_l = orig_dim;
  const uintptr_t final_dim = orig_dim_l + 1;
  uintptr_t orig_row_idx = 0;
  const double* a_inv_row = a_inv;
  double* outmatrix_row = outmatrix;
  for (; orig_row_idx != insert_idx; ++orig_row_idx) {
    const double ainv_b_div_k = k_recip * ainv_b_buf[orig_row_idx];
    for (uintptr_t col_idx = 0; col_idx <= orig_row_idx; ++col_idx) {
      outmatrix_row[col_idx] = a_inv_row[col_idx] + ainv_b_div_k * ainv_b_buf[col_idx];
    }
    a_inv_row = &(a_inv_row[orig_dim_l]);
    outmatrix_row = &(outmatrix_row[final_dim]);
  }
  for (uintptr_t col_idx = 0; col_idx != insert_idx; ++col_idx) {
    outmatrix_row[col_idx] = -k_recip * ainv_b_buf[col_idx];
  }
  outmatrix_row[insert_idx] = k_recip;
  for (; orig_row_idx != orig_dim_l; ++orig_row_idx) {
    outmatrix_row = &(outmatrix_row[final_dim]);
    const double ainv_b_div_k = k_recip * ainv_b_buf[orig_row_idx];
    for (uintptr_t col_idx = 0; col_idx != insert_idx; ++col_idx) {
      outmatrix_row[col_idx] = a_inv_row[col_idx] + ainv_b_div_k * ainv_b_buf[col_idx];
    }
    outmatrix_row[insert_idx] = -ainv_b_div_k;
    double* outmatrix_write_base = &(outmatrix_row[1]);
    for (uintptr_t orig_col_idx = insert_idx; orig_col_idx <= orig_row_idx; ++orig_col_idx) {
      outmatrix_write_base[orig_col_idx] = a_inv_row[orig_col_idx] + ainv_b_div_k * ainv_b_buf[orig_col_idx];
    }
    a_inv_row = &(a_inv_row[orig_dim_l]);
  }
  return 0;
}

BoolErr InvertRank1SymmDiag(const double* a_inv, const double* bb, lapack_int orig_dim, double cc, double* __restrict outdiag, double* __restrict ainv_b_buf) {
  double k_recip;
  if (invert_rank1_symm_start(a_inv, bb, orig_dim, cc, ainv_b_buf, &k_recip)) {
    return 1;
  }
  const uintptr_t orig_dim_l = orig_dim;
  const uintptr_t orig_dim_p1 = orig_dim_l + 1;
  for (uintptr_t ulii = 0; ulii != orig_dim_l; ++ulii) {
    const double dxx = ainv_b_buf[ulii];
    outdiag[ulii] = a_inv[ulii * orig_dim_p1] + k_recip * dxx * dxx;
  }
  outdiag[orig_dim_l] = k_recip;
  return 0;
}

BoolErr InvertRank2SymmStart(const double* a_inv, const double* bb, lapack_int orig_dim, lapack_int b_stride, double d11, double d12, double d22, double* __restrict b_ainv, double* __restrict s_b_ainv, double* __restrict schur11_ptr, double* __restrict schur12_ptr, double* __restrict schur22_ptr) {
  const uintptr_t orig_dim_l = orig_dim;
  if (orig_dim) {
    // (have confirmed that dgemm is better than per-row multiplies even for
    // just 2 rows)
    RowMajorMatrixMultiplyStrided(bb, a_inv, 2, b_stride, orig_dim, orig_dim, orig_dim, orig_dim, b_ainv);

    // Schur complement = (D - B A^{-1} B^T)^{-1}
    if (orig_dim_l > kDotprodDThresh) {
      d11 -= DotprodD(bb, b_ainv, orig_dim);
      d12 -= DotprodD(bb, &(b_ainv[orig_dim_l]), orig_dim);
      d22 -= DotprodD(&(bb[b_stride]), &(b_ainv[orig_dim_l]), orig_dim);
    } else {
      d11 -= DotprodDShort(bb, b_ainv, orig_dim);
      d12 -= DotprodDShort(bb, &(b_ainv[orig_dim_l]), orig_dim);
      d22 -= DotprodDShort(&(bb[b_stride]), &(b_ainv[orig_dim_l]), orig_dim);
    }
  }

  // [ a b ]^{-1} = [ d  -b ] / (ad - b^2)
  // [ b d ]        [ -b a  ]
  const double det = d11 * d22 - d12 * d12;
  if (fabs(det) < kMatrixSingularRcond) {
    return 1;
  }
  const double det_recip = 1.0 / det;
  const double schur11 = d22 * det_recip;
  const double schur12 = -d12 * det_recip;
  const double schur22 = d11 * det_recip;
  for (uintptr_t col_idx = 0; col_idx != orig_dim_l; ++col_idx) {
    const double b_ainv_1 = b_ainv[col_idx];
    const double b_ainv_2 = b_ainv[col_idx + orig_dim_l];
    s_b_ainv[col_idx] = schur11 * b_ainv_1 + schur12 * b_ainv_2;
    s_b_ainv[col_idx + orig_dim_l] = schur12 * b_ainv_1 + schur22 * b_ainv_2;
  }
  *schur11_ptr = schur11;
  *schur12_ptr = schur12;
  *schur22_ptr = schur22;
  return 0;
}

BoolErr InvertRank2Symm(const double* a_inv, const double* bb, lapack_int orig_dim, lapack_int b_stride, uint32_t insert_idx, double d11, double d12, double d22, double* __restrict outmatrix, double* __restrict b_ainv_buf, double* __restrict s_b_ainv_buf) {
  double schur11;
  double schur12;
  double schur22;
  if (InvertRank2SymmStart(a_inv, bb, orig_dim, b_stride, d11, d12, d22, b_ainv_buf, s_b_ainv_buf, &schur11, &schur12, &schur22)) {
    return 1;
  }
  // [ A^{-1} + A{-1}B^T S B A{-1}   -A^{-1}B^T S ]
  // [         -S B A{-1}                  S      ]
  const uintptr_t orig_dim_l = orig_dim;
  const uintptr_t final_dim = orig_dim_l + 2;
  uintptr_t orig_row_idx = 0;
  const double* a_inv_row = a_inv;
  const double* b_ainv_row2 = &(b_ainv_buf[orig_dim_l]);
  const double* s_b_ainv_row2 = &(s_b_ainv_buf[orig_dim_l]);
  double* outmatrix_row = outmatrix;
  for (; orig_row_idx != insert_idx; ++orig_row_idx) {
    const double b_ainv_1 = b_ainv_buf[orig_row_idx];
    const double b_ainv_2 = b_ainv_row2[orig_row_idx];
    for (uintptr_t col_idx = 0; col_idx <= orig_row_idx; ++col_idx) {
      outmatrix_row[col_idx] = a_inv_row[col_idx] + b_ainv_1 * s_b_ainv_buf[col_idx] + b_ainv_2 * s_b_ainv_row2[col_idx];
    }
    a_inv_row = &(a_inv_row[orig_dim_l]);
    outmatrix_row = &(outmatrix_row[final_dim]);
  }
  for (uintptr_t col_idx = 0; col_idx != insert_idx; ++col_idx) {
    outmatrix_row[col_idx] = -s_b_ainv_buf[col_idx];
  }
  outmatrix_row[insert_idx] = schur11;
  outmatrix_row = &(outmatrix_row[final_dim]);
  for (uintptr_t col_idx = 0; col_idx != insert_idx; ++col_idx) {
    outmatrix_row[col_idx] = -s_b_ainv_row2[col_idx];
  }
  outmatrix_row[insert_idx] = schur12;
  outmatrix_row[insert_idx + 1] = schur22;
  for (; orig_row_idx != orig_dim_l; ++orig_row_idx) {
    outmatrix_row = &(outmatrix_row[final_dim]);
    const double b_ainv_1 = b_ainv_buf[orig_row_idx];
    const double b_ainv_2 = b_ainv_row2[orig_row_idx];
    for (uintptr_t col_idx = 0; col_idx != insert_idx; ++col_idx) {
      outmatrix_row[col_idx] = a_inv_row[col_idx] + b_ainv_1 * s_b_ainv_buf[col_idx] + b_ainv_2 * s_b_ainv_row2[col_idx];
    }
    outmatrix_row[insert_idx] = -s_b_ainv_buf[orig_row_idx];
    outmatrix_row[insert_idx + 1] = -s_b_ainv_row2[orig_row_idx];
    double* outmatrix_write_base = &(outmatrix_row[2]);
    for (uintptr_t orig_col_idx = insert_idx; orig_col_idx <= orig_row_idx; ++orig_col_idx) {
      outmatrix_write_base[orig_col_idx] = a_inv_row[orig_col_idx] + b_ainv_1 * s_b_ainv_buf[orig_col_idx] + b_ainv_2 * s_b_ainv_row2[orig_col_idx];
    }
    a_inv_row = &(a_inv_row[orig_dim_l]);
  }
  return 0;
}

BoolErr InvertRank2SymmDiag(const double* a_inv, const double* bb, lapack_int orig_dim, double d11, double d12, double d22, double* __restrict outdiag, double* __restrict b_ainv_buf, double* __restrict s_b_ainv_buf) {
  double schur11;
  double schur12;
  double schur22;
  if (InvertRank2SymmStart(a_inv, bb, orig_dim, orig_dim, d11, d12, d22, b_ainv_buf, s_b_ainv_buf, &schur11, &schur12, &schur22)) {
    return 1;
  }
  const uintptr_t orig_dim_l = orig_dim;
  const uintptr_t orig_dim_p1 = orig_dim_l + 1;
  const double* b_ainv_row2 = &(b_ainv_buf[orig_dim_l]);
  const double* s_b_ainv_row2 = &(s_b_ainv_buf[orig_dim_l]);
  for (uintptr_t ulii = 0; ulii != orig_dim_l; ++ulii) {
    outdiag[ulii] = a_inv[ulii * orig_dim_p1] + b_ainv_buf[ulii] * s_b_ainv_buf[ulii] + b_ainv_row2[ulii] * s_b_ainv_row2[ulii];
  }
  outdiag[orig_dim_l] = schur11;
  outdiag[orig_dim_l + 1] = schur22;
  return 0;
}

// now assumes xtx_inv is predictors_pmaj * transpose on input
#ifdef NOLAPACK
BoolErr LinearRegressionInvMain(const double* xt_y_phenomaj, uint32_t predictor_ct, uint32_t pheno_ct, double* xtx_inv, double* fitted_coefs_phenomaj, MatrixInvertBuf1* mi_buf, double* dbl_2d_buf) {
  if (InvertSymmdefMatrix(predictor_ct, xtx_inv, mi_buf, dbl_2d_buf)) {
    return 1;
  }
  RowMajorMatrixMultiply(xt_y_phenomaj, xtx_inv, pheno_ct, predictor_ct, predictor_ct, fitted_coefs_phenomaj);
  return 0;
}
#else

BoolErr LinearRegressionInvMain(const double* xt_y_phenomaj, uint32_t predictor_ct, lapack_int pheno_ct, double* xtx_inv, double* fitted_coefs_phenomaj) {
  char uplo = 'U';
  lapack_int tmp_n = predictor_ct;
  lapack_int info;
  LAPACK_dpotrf(&uplo, &tmp_n, xtx_inv, &tmp_n, &info);
  if (info) {
    return 1;
  }
  memcpy(fitted_coefs_phenomaj, xt_y_phenomaj, predictor_ct * pheno_ct * sizeof(double));
  LAPACK_dpotrs(&uplo, &tmp_n, &pheno_ct, xtx_inv, &tmp_n, fitted_coefs_phenomaj, &tmp_n, &info);
  assert(!info);
  LAPACK_dpotri(&uplo, &tmp_n, xtx_inv, &tmp_n, &info);
  return (info != 0);
}
#endif  // !NOLAPACK

#ifdef __cplusplus
}  // namespace plink2
#endif
