

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __LAPACK_FUNC_H__
#define __LAPACK_FUNC_H__

bool svd_lapack(int,vector_t & A, vector_t & S,  matrix_t & V);
bool eigen_lapack(int,vector_t & A, vector_t & S, matrix_t & V);

#endif
