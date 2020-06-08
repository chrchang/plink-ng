#ifndef __PLINK2_MATRIX_CUDA_H__
#define __PLINK2_MATRIX_CUDA_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

// CUDA-specific portion of plink2_matrix, to be compiled with nvcc instead of
// gcc/clang.

#ifdef __cplusplus
extern "C" {
#endif

// Avoid importing cuda_runtime.h and cublas_v2.h here, since it uses raw
// casts, etc.

int CudaGetDeviceCount();

int CudaSetDevice(int device_idx);

typedef struct CublasFmultiplierStruct {
  void* handle; // known to be a pointer
  float* dev_inmatrix1;
  float* dev_inmatrix2;
  float* dev_outmatrix;
  int row1_ct;
  int col2_ct;
  int common_ct;
  unsigned int in_borrowed;
  int device_idx;
} CublasFmultiplier;

void CublasFmultiplierPreinit(CublasFmultiplier* cfmp);

// These functions allow multiple CublasFmultipliers to share an input matrix
// loaded in device memory.
void CublasFmultiplierBorrowColMajor1(CublasFmultiplier* srcp, CublasFmultiplier* dstp);

void CublasFmultiplierBorrowColMajor2(CublasFmultiplier* srcp, CublasFmultiplier* dstp);

inline void CublasFmultiplierBorrowRowMajor1(CublasFmultiplier* srcp, CublasFmultiplier* dstp) {
  CublasFmultiplierBorrowColMajor2(srcp, dstp);
}

inline void CublasFmultiplierBorrowRowMajor2(CublasFmultiplier* srcp, CublasFmultiplier* dstp) {
  CublasFmultiplierBorrowColMajor1(srcp, dstp);
}

// Initializes device handle and allocates device memory.  (Borrows should
// happen first.)
// Returns 0 on success, 1 on failure.
int CublasFmultiplierColMajorInit(int row1_ct, int col2_ct, int common_ct, CublasFmultiplier* cfmp);

inline int CublasFmultiplierRowMajorInit(int row1_ct, int col2_ct, int common_ct, CublasFmultiplier* cfmp) {
  return CublasFmultiplierColMajorInit(col2_ct, row1_ct, common_ct, cfmp);
}

int CublasFmultiplierPreloadColMajor1(const float* inmatrix1, CublasFmultiplier* cfmp);

int CublasFmultiplierPreloadColMajor2(const float* inmatrix2, CublasFmultiplier* cfmp);

inline int CublasFmultiplierPreloadRowMajor1(const float* inmatrix1, CublasFmultiplier* cfmp) {
  return CublasFmultiplierPreloadColMajor2(inmatrix1, cfmp);
}

inline int CublasFmultiplierPreloadRowMajor2(const float* inmatrix2, CublasFmultiplier* cfmp) {
  return CublasFmultiplierPreloadColMajor1(inmatrix2, cfmp);
}

// Assumes CublasFmultiplyPreloadColMajor2() was previously called.
int CublasFmultiplyColMajor1(const float* inmatrix1, CublasFmultiplier* cfmp, float* outmatrix);

int CublasFmultiplyColMajor2(const float* inmatrix2, CublasFmultiplier* cfmp, float* outmatrix);

inline int CublasFmultiplyRowMajor1(const float* inmatrix1, CublasFmultiplier* cfmp, float* outmatrix) {
  return CublasFmultiplyColMajor2(inmatrix1, cfmp, outmatrix);
}

inline int CublasFmultiplyRowMajor2(const float* inmatrix2, CublasFmultiplier* cfmp, float* outmatrix) {
  return CublasFmultiplyColMajor1(inmatrix2, cfmp, outmatrix);
}

int CublasFmultiplierCleanup(CublasFmultiplier* cfmp);


// Only for testing purposes: inefficient host-memory reallocation, etc.
int CublasColMajorFmatrixMultiply(const float* inmatrix1, const float* inmatrix2, int row1_ct, int col2_ct, int common_ct, float* outmatrix);

inline int CublasRowMajorFmatrixMultiply(const float* inmatrix1, const float* inmatrix2, int row1_ct, int col2_ct, int common_ct, float* outmatrix) {
  return CublasColMajorFmatrixMultiply(inmatrix2, inmatrix1, col2_ct, row1_ct, common_ct, outmatrix);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif  // __PLINK2_MATRIX_CUDA_H__
