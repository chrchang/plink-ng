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


#include "plink_common.h"

#include "plink_matrix.h"

#ifndef NOLAPACK
#  ifndef __APPLE__
  void xerbla_(void) {} // fix static linking error
#  endif
#endif

static inline double SQR(const double a) {
  return a * a;
}

double pythag(const double a, const double b) {
  // PLINK stats.cpp pythag().
  double absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

#ifdef NOLAPACK
  #ifdef __cplusplus
static inline double SIGN(const double &a, const double &b) {
  // PLINK helper.h SIGN() template specialized to doubles.
  return (b >= 0)? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
  #else
static inline double SIGN(const double a, const double b) {
  // PLINK helper.h SIGN() template specialized to doubles.
  return (b >= 0)? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
  #endif

int32_t svdcmp_c(int32_t m, double* a, double* w, double* v) {
  // C port of PLINK stats.cpp svdcmp().
  // now thread-safe.
  double* rv1 = &(w[(uint32_t)m]);
  int32_t n = m;
  int32_t flag;
  int32_t l = 0; // suppress compile warning
  int32_t i,its,j,jj,k,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double temp;

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
  return 1;
}

int32_t invert_matrix(int32_t dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* dbl_1d_buf, double* dbl_2d_buf) {
  // C port of PLINK stats.cpp's svd_inverse() function.
  // Now thread-safe in NOLAPACK case.

  // w -> dbl_1d_buf
  // v -> dbl_2d_buf
  const double eps = 1e-24;
  int32_t i;
  int32_t j;
  int32_t k;
  i = svdcmp_c(dim, matrix, dbl_1d_buf, dbl_2d_buf);
  if (i == -1) {
    return -1;
  } else if (!i) {
    return 1;
  }

  // Look for singular values
  double wmax = 0;
  for (i=0; i<dim; i++) {
    wmax = dbl_1d_buf[i] > wmax ? dbl_1d_buf[i] : wmax;
  }
  double wmin = wmax * eps;
  for (i=0; i<dim; i++) {
    dbl_1d_buf[i] = dbl_1d_buf[i] < wmin ? 0 : (1 / dbl_1d_buf[i]);
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      matrix[i * dim + j] = matrix[i * dim + j] * dbl_1d_buf[j];
    }
  }

  // [nxn].[t(v)]
  for (i=0; i<dim; i++) {
    fill_double_zero(dim, dbl_1d_buf);
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
#else  // !NOLAPACK
int32_t invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf) {
  // dgetrf_/dgetri_ is more efficient than dpotrf_/dpotri_ on OS X.
  __CLPK_integer lwork = dim * dim;
  __CLPK_integer info;
  dgetrf_(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  dgetri_(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  if (info) {
    return 1;
  }
  return 0;
}

int32_t invert_matrix_checked(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf) {
  // This used to fall back on PLINK 1.07's SVD-based implementation when the
  // rcond estimate was too small, but in practice that just slowed things down
  // without meaningfully improving inversion of nonsingular matrices.  So now
  // this just exits a bit earlier, while leaving the old "binary search for
  // the first row/column causing multicollinearity" logic to the caller.
  __CLPK_integer lwork = dim * dim;
  char cc = '1';
  double norm = dlange_(&cc, &dim, &dim, matrix, &dim, dbl_2d_buf);
  __CLPK_integer info;
  double rcond;
  dgetrf_(&dim, &dim, matrix, &dim, int_1d_buf, &info);
  if (info > 0) {
    return 1;
  }
  dgecon_(&cc, &dim, matrix, &dim, &norm, &rcond, dbl_2d_buf, &(int_1d_buf[dim]), &info);
  if (rcond < MATRIX_SINGULAR_RCOND) {
    return 1;
  }
  dgetri_(&dim, matrix, &dim, int_1d_buf, dbl_2d_buf, &lwork, &info);
  return 0;
}
#endif  // !NOLAPACK

void col_major_matrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, double* inmatrix1, double* inmatrix2, double* outmatrix) {
#ifdef NOLAPACK
  uintptr_t row1_ct_l = row1_ct;
  uintptr_t col2_ct_l = col2_ct;
  uintptr_t common_ct_l = common_ct;
  uintptr_t row_idx;
  uintptr_t col_idx;
  uintptr_t com_idx;
  double* dptr;
  double dxx;
  // not optimized
  for (col_idx = 0; col_idx < col2_ct_l; col_idx++) {
    for (row_idx = 0; row_idx < row1_ct_l; row_idx++) {
      dxx = 0;
      dptr = &(inmatrix2[col_idx * common_ct]);
      for (com_idx = 0; com_idx < common_ct_l; com_idx++) {
        dxx += (*dptr++) * inmatrix1[com_idx * row1_ct_l + row_idx];
      }
      *outmatrix++ = dxx;
    }
  }
#else
#  ifndef USE_CBLAS_XGEMM
  char blas_char = 'N';
  double dyy = 1;
  double dzz = 0;
  dgemm_(&blas_char, &blas_char, &row1_ct, &col2_ct, &common_ct, &dyy, inmatrix1, &row1_ct, inmatrix2, &common_ct, &dzz, outmatrix, &row1_ct);
#  else
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row1_ct, col2_ct, common_ct, 1.0, inmatrix1, row1_ct, inmatrix2, common_ct, 0.0, outmatrix, row1_ct);
#  endif  // USE_CBLAS_XGEMM
#endif  // !NOLAPACK
}

void col_major_fmatrix_multiply(__CLPK_integer row1_ct, __CLPK_integer col2_ct, __CLPK_integer common_ct, float* inmatrix1, float* inmatrix2, float* outmatrix) {
#ifdef NOLAPACK
  uintptr_t row1_ct_l = row1_ct;
  uintptr_t col2_ct_l = col2_ct;
  uintptr_t common_ct_l = common_ct;
  uintptr_t row_idx;
  uintptr_t col_idx;
  uintptr_t com_idx;
  float* fptr;
  float fxx;
  // not optimized
  for (col_idx = 0; col_idx < col2_ct_l; col_idx++) {
    for (row_idx = 0; row_idx < row1_ct_l; row_idx++) {
      fxx = 0;
      fptr = &(inmatrix2[col_idx * common_ct]);
      for (com_idx = 0; com_idx < common_ct_l; com_idx++) {
        fxx += (*fptr++) * inmatrix1[com_idx * row1_ct_l + row_idx];
      }
      *outmatrix++ = fxx;
    }
  }
#else
#  ifndef USE_CBLAS_XGEMM
  char blas_char = 'N';
  float fyy = 1;
  float fzz = 0;
  sgemm_(&blas_char, &blas_char, &row1_ct, &col2_ct, &common_ct, &fyy, inmatrix1, &row1_ct, inmatrix2, &common_ct, &fzz, outmatrix, &row1_ct);
#  else
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row1_ct, col2_ct, common_ct, 1.0, inmatrix1, row1_ct, inmatrix2, common_ct, 0.0, outmatrix, row1_ct);
#  endif  // USE_CBLAS_XGEMM
#endif  // !NOLAPACK
}

// Todo: replace these with cache-oblivious, or at least -friendlier,
// algorithms.
void transpose_copy(uintptr_t old_maj, uintptr_t new_maj, double* old_matrix, double* new_matrix) {
  double* dptr;
  uintptr_t new_maj_idx;
  uintptr_t old_maj_idx;
  for (new_maj_idx = 0; new_maj_idx < new_maj; new_maj_idx++) {
    dptr = &(old_matrix[new_maj_idx]);
    for (old_maj_idx = 0; old_maj_idx < old_maj; old_maj_idx++) {
      *new_matrix++ = dptr[old_maj_idx * new_maj];
    }
  }
}

void transpose_copy_float(uintptr_t old_maj, uintptr_t new_maj, uint32_t new_maj_max, float* old_matrix, float* new_matrix) {
  // new_maj = in-memory stride of old_matrix rows
  // new_maj_max = actual number of rows in new_matrix
  // (distinction is necessary for SSE alignment)
  float* fptr;
  uintptr_t new_maj_idx;
  uintptr_t old_maj_idx;
  for (new_maj_idx = 0; new_maj_idx < new_maj_max; new_maj_idx++) {
    fptr = &(old_matrix[new_maj_idx]);
    for (old_maj_idx = 0; old_maj_idx < old_maj; old_maj_idx++) {
      *new_matrix++ = fptr[old_maj_idx * new_maj];
    }
  }
}
