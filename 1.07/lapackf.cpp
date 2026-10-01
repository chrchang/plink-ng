


//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>


#include "plink.h"
#include "helper.h"

#ifdef WITH_LAPACK
#include "lapackf.h"

extern "C" int dgesdd_(char *jobz, int *m, int *n, double *a, 
		       int *lda, double *s, double *u, int *ldu, 
		       double *vt, int *ldvt, double *work, int *lwork, 
		       int *iwork, int *info);

extern "C" int dsyevx_( char * , char * , char * , int * ,
			double * , int * , double * , double * , int * ,
			int * , double * , int * , double * ,
			double * , int * , double * , int * ,
			int * , int * , int * ) ;
#endif


bool svd_lapack(int n, vector_t & A, vector_t & S, matrix_t & V)
{
  int m=n;
  vector_t U(n*n);
  vector_t tV(n*n);

  cout << "Using LAPACK SVD library function...\n";
    
#ifdef WITH_LAPACK
    
  int info=0;
    
  vector<int> iwork(8*m,0);    
  double optim_lwork;
  int lwork;
  lwork = -1;
  
  // Determine workspace needed
  dgesdd_("A", &m, &n, &A[0] , &m, &S[0], &U[0], &m, &tV[0], &n, &optim_lwork, &lwork, &iwork[0], &info);
  lwork = (int) optim_lwork;
  vector_t work( lwork, 0 );
  
  // Perform actual SVD
  dgesdd_("A", &m, &n, &A[0] , &m, &S[0], &U[0], &m, &tV[0], &n, &work[0], &lwork, &iwork[0], &info);
  
  // Copy and transpose V
  int k = 0;
  for( int i = 0; i < n; i++ )
    for( int j = 0; j < n; j++ )
      {
	V[j][i] = tV[k];
	++k;
      }
  
  return true;
  
#else
  
  // LAPACK support not compiled 
  return false;
  
#endif
  
}


bool eigen_lapack(int n, vector_t & A, vector_t & S, matrix_t & V)
{
    
  // Use eigenvalue decomposition instead of SVD
  // Get only the highest eigen-values, (par::cluster_mds_dim)
  
  int i1 = n - par::cluster_mds_dim + 1;
  int i2 = n;
  double z = -1;
  
  // Integer workspace size, 5N
  vector<int> iwork(5*n,0);    
  
  double optim_lwork;
  int lwork = -1;
  
  int out_m;
  vector_t out_w( par::cluster_mds_dim , 0 );
  vector_t out_z( n * par::cluster_mds_dim ,0 );
  
  int ldz = n;
  vector<int> ifail(n,0);
  int info=0;
  double nz = 0;
  
  // Get workspace
  
  dsyevx_("V" ,         // get eigenvalues and eigenvectors
	  "I" ,         // get interval of selected eigenvalues
	  "L" ,         // data stored as upper triangular
	  &n  ,         // order of matrix
	  &A[0] ,       // input matrix
	  &n ,          // LDA
	  &nz ,         // Vlower
	  &nz ,         // Vupper
	  &i1,          // from 1st ...
	  &i2,          // ... to nth eigenvalue
	  &z ,          // 0 for ABSTOL
	  &out_m,       // # of eigenvalues found
	  &out_w[0],    // first M entries contain sorted eigen-values
	  &out_z[0],    // array (can be mxm? nxn)
	  &ldz,         // make n at first
	  &optim_lwork, // Get optimal workspace 
	  &lwork,       // size of workspace
	  &iwork[0],    // int workspace
	  &ifail[0],    // output: failed to converge
	  &info );
  
  // Assign workspace
  
  lwork = (int) optim_lwork;
  vector_t work( lwork, 0 );

  dsyevx_("V" ,      // get eigenvalues and eigenvectors
	  "I" ,      // get interval of selected eigenvalues
	  "L" ,      // data stored as upper triangular
	  &n  ,      // order of matrix
	  &A[0] ,    // input matrix
	  &n ,       // LDA
	  &nz ,      // Vlower
	  &nz ,      // Vupper
	  &i1,       // from 1st ...
	  &i2,       // ... to nth eigenvalue
	  &z ,       // 0 for ABSTOL
	  &out_m,    // # of eigenvalues found
	  &out_w[0], // first M entries contain sorted eigen-values
	  &out_z[0], // array (can be mxm? nxn)
	  &ldz,      // make n at first
	  &work[0],  // Workspace
	  &lwork,    // size of workspace
	  &iwork[0], // int workspace
	  &ifail[0], // output: failed to converge
	  &info );      

  // Get eigenvalues, vectors
  for (int i=0; i< par::cluster_mds_dim; i++)
    S[i] = out_w[i];
  
  for (int i=0; i<n; i++)
    for (int j=0;j<par::cluster_mds_dim; j++)
      V[i][j] = out_z[ i + j*n ];
  
  return true;
  
}
