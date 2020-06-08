#include "cuda_runtime.h"
#include "cublas_v2.h"

#include "plink2_matrix_cuda.h"
#include <stdio.h>

int CudaGetDeviceCount() {
  int device_count;
  if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
    return 0;
  }
  return device_count;
}

int CudaSetDevice(int device_idx) {
  if (cudaSetDevice(device_idx) != cudaSuccess) {
    return 1;
  }
  return 0;
}

void CublasFmultiplierPreinit(CublasFmultiplier* cfmp) {
  cfmp->handle = NULL;
  cfmp->dev_inmatrix1 = NULL;
  cfmp->dev_inmatrix2 = NULL;
  cfmp->dev_outmatrix = NULL;
  cfmp->in_borrowed = 0;
  cfmp->device_idx = -1;
}

void CublasFmultiplierBorrowColMajor1(CublasFmultiplier* srcp, CublasFmultiplier* dstp) {
  dstp->dev_inmatrix1 = srcp->dev_inmatrix1;
  dstp->in_borrowed |= 1;
}

void CublasFmultiplierBorrowColMajor2(CublasFmultiplier* srcp, CublasFmultiplier* dstp) {
  dstp->dev_inmatrix2 = srcp->dev_inmatrix2;
  dstp->in_borrowed |= 2;
}

// Initializes device handle and allocates device memory.
// Returns 0 on success, 1 on failure.
int CublasFmultiplierColMajorInit(int row1_ct, int col2_ct, int common_ct, CublasFmultiplier* cfmp) {
  if ((cfmp->device_idx != -1) || (cfmp->handle != NULL)) {
    // sanity check: improper function call
    return 1;
  }
  if (cudaGetDevice(&cfmp->device_idx) != cudaSuccess) {
    cfmp->device_idx = -1;
    return 1;
  }
  cublasHandle_t* handlep = (cublasHandle_t*)(&cfmp->handle);
  cublasStatus_t stat = cublasCreate(handlep);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    return 1;
  }
  if (!(cfmp->in_borrowed & 1)) {
    cudaError_t cudaStat = cudaMalloc((void**)&cfmp->dev_inmatrix1, row1_ct * sizeof(float) * common_ct);
    if (cudaStat != cudaSuccess) {
      return 1;
    }
  }
  if (!(cfmp->in_borrowed & 2)) {
    cudaError_t cudaStat = cudaMalloc((void**)&cfmp->dev_inmatrix2, common_ct * sizeof(float) * col2_ct);
    if (cudaStat != cudaSuccess) {
      return 1;
    }
  }
  cudaError_t cudaStat = cudaMalloc((void**)&cfmp->dev_outmatrix, row1_ct * sizeof(float) * col2_ct);
  if (cudaStat != cudaSuccess) {
    return 1;
  }
  cfmp->row1_ct = row1_ct;
  cfmp->col2_ct = col2_ct;
  cfmp->common_ct = common_ct;
  return 0;
}

int CublasFmultiplierPreloadColMajor1(const float* inmatrix1, CublasFmultiplier* cfmp) {
  const int row1_ct = cfmp->row1_ct;
  cublasStatus_t stat = cublasSetMatrix(row1_ct, cfmp->common_ct, sizeof(float), inmatrix1, row1_ct, cfmp->dev_inmatrix1, row1_ct);
  return (stat != CUBLAS_STATUS_SUCCESS);
}

int CublasFmultiplierPreloadColMajor2(const float* inmatrix2, CublasFmultiplier* cfmp) {
  const int common_ct = cfmp->common_ct;
  cublasStatus_t stat = cublasSetMatrix(common_ct, cfmp->col2_ct, sizeof(float), inmatrix2, common_ct, cfmp->dev_inmatrix2, common_ct);
  return (stat != CUBLAS_STATUS_SUCCESS);
}

int CublasFmultiplyColMajor1(const float* inmatrix1, CublasFmultiplier* cfmp, float* outmatrix) {
  const int row1_ct = cfmp->row1_ct;
  const int col2_ct = cfmp->col2_ct;
  const int common_ct = cfmp->common_ct;
  cublasStatus_t stat = cublasSetMatrix(row1_ct, common_ct, sizeof(float), inmatrix1, row1_ct, cfmp->dev_inmatrix1, row1_ct);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    return 1;
  }
  float alpha = 1.0;
  float beta = 0.0;
  stat = cublasSgemm((cublasHandle_t)cfmp->handle, CUBLAS_OP_N, CUBLAS_OP_N, row1_ct, col2_ct, common_ct, &alpha, cfmp->dev_inmatrix1, row1_ct, cfmp->dev_inmatrix2, common_ct, &beta, cfmp->dev_outmatrix, row1_ct);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    return 1;
  }
  stat = cublasGetMatrix(row1_ct, col2_ct, sizeof(float), cfmp->dev_outmatrix, row1_ct, outmatrix, row1_ct);
  return (stat != CUBLAS_STATUS_SUCCESS);
}

int CublasFmultiplyColMajor2(const float* inmatrix2, CublasFmultiplier* cfmp, float* outmatrix) {
  const int row1_ct = cfmp->row1_ct;
  const int col2_ct = cfmp->col2_ct;
  const int common_ct = cfmp->common_ct;
  cublasStatus_t stat = cublasSetMatrix(common_ct, col2_ct, sizeof(float), inmatrix2, common_ct, cfmp->dev_inmatrix2, common_ct);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    return 1;
  }
  float alpha = 1.0;
  float beta = 0.0;
  stat = cublasSgemm((cublasHandle_t)cfmp->handle, CUBLAS_OP_N, CUBLAS_OP_N, row1_ct, col2_ct, common_ct, &alpha, cfmp->dev_inmatrix1, row1_ct, cfmp->dev_inmatrix2, common_ct, &beta, cfmp->dev_outmatrix, row1_ct);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    return 1;
  }
  stat = cublasGetMatrix(row1_ct, col2_ct, sizeof(float), cfmp->dev_outmatrix, row1_ct, outmatrix, row1_ct);
  return (stat != CUBLAS_STATUS_SUCCESS);
}

int CublasFmultiplierCleanup(CublasFmultiplier* cfmp) {
  if (cfmp->device_idx == -1) {
    return 0;
  }
  int retval = 0;
  if (cudaSetDevice(cfmp->device_idx) != cudaSuccess) {
    retval = 1;
  }
  cfmp->device_idx = -1;
  if (cudaFree(cfmp->dev_outmatrix) != cudaSuccess) {
    retval = 1;
  }
  cfmp->dev_outmatrix = NULL;
  if (!(cfmp->in_borrowed & 2)) {
    if (cudaFree(cfmp->dev_inmatrix2) != cudaSuccess) {
      retval = 1;
    }
  }
  cfmp->dev_inmatrix2 = NULL;
  if (!(cfmp->in_borrowed & 1)) {
    if (cudaFree(cfmp->dev_inmatrix1) != cudaSuccess) {
      retval = 1;
    }
  }
  cfmp->dev_inmatrix1 = NULL;
  cfmp->in_borrowed = 0;
  cublasDestroy((cublasHandle_t)cfmp->handle); // don't care about return value
  cfmp->handle = NULL;
  return retval;
}

// Only for testing purposes.
int CublasColMajorFmatrixMultiply(const float* inmatrix1, const float* inmatrix2, int row1_ct, int col2_ct, int common_ct, float* outmatrix) {
  CublasFmultiplier cfm;
  CublasFmultiplierPreinit(&cfm);
  int retval = 0;
  {
    if (CublasFmultiplierColMajorInit(row1_ct, col2_ct, common_ct, &cfm)) {
      goto CublasColMajorFmatrixMultiply_fail;
    }
    cublasStatus_t stat = cublasSetMatrix(row1_ct, common_ct, sizeof(float), inmatrix1, row1_ct, cfm.dev_inmatrix1, row1_ct);
    if (stat != CUBLAS_STATUS_SUCCESS) {
      goto CublasColMajorFmatrixMultiply_fail;
    }
    stat = cublasSetMatrix(common_ct, col2_ct, sizeof(float), inmatrix2, common_ct, cfm.dev_inmatrix2, common_ct);
    if (stat != CUBLAS_STATUS_SUCCESS) {
      goto CublasColMajorFmatrixMultiply_fail;
    }
    float alpha = 1.0;
    float beta = 0.0;
    stat = cublasSgemm((cublasHandle_t)cfm.handle, CUBLAS_OP_N, CUBLAS_OP_N, row1_ct, col2_ct, common_ct, &alpha, cfm.dev_inmatrix1, row1_ct, cfm.dev_inmatrix2, common_ct, &beta, cfm.dev_outmatrix, row1_ct);
    if (stat != CUBLAS_STATUS_SUCCESS) {
      goto CublasColMajorFmatrixMultiply_fail;
    }
    stat = cublasGetMatrix(row1_ct, col2_ct, sizeof(float), cfm.dev_outmatrix, row1_ct, outmatrix, row1_ct);
    if (stat != CUBLAS_STATUS_SUCCESS) {
      goto CublasColMajorFmatrixMultiply_fail;
    }
  }
  while (0) {
  CublasColMajorFmatrixMultiply_fail:
    retval = 1;
  }
  if (CublasFmultiplierCleanup(&cfm)) {
    retval = 1;
  }
  return retval;
}
