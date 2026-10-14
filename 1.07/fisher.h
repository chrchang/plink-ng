

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2001    Robert Gentleman, Ross Ihaka 
 *                             and the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Memory Allocation (garbage collected) --- INCLUDING S compatibility ---
 */

#ifndef __FISHER_H__
#define __FISHER_H__

#include "plink.h"

#ifndef R_EXT_MEMORY_H_
#define R_EXT_MEMORY_H_

#ifdef  __cplusplus
extern "C" {
#endif

  char * vmaxget(void);
  void vmaxset(char*);

  void R_gc(void);

  char * R_alloc(long, int);
  char * S_alloc(long, int);
  char * S_realloc(char*, long, long, int);

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_MEMORY_H_ */


#ifndef R_EXT_BOOLEAN_H_
#define R_EXT_BOOLEAN_H_

#undef FALSE
#undef TRUE

#ifdef  __cplusplus
extern "C" {
#endif
  typedef enum { FALSE = 0, TRUE /*, MAYBE */ } Rboolean;

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_BOOLEAN_H_ */

#ifndef R_EXT_CONSTANTS_H_
#define R_EXT_CONSTANTS_H_

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

#define PI             M_PI
#define SINGLE_EPS     FLT_EPSILON
#define SINGLE_BASE    FLT_RADIX
#define SINGLE_XMIN    FLT_MIN
#define SINGLE_XMAX    FLT_MAX
#define DOUBLE_DIGITS  DBL_MANT_DIG
#define DOUBLE_EPS     DBL_EPSILON
#define DOUBLE_XMAX    DBL_MAX
#define DOUBLE_XMIN    DBL_MIN

#endif


// Fisher's exact test

void fexact(int *nrow, 
	    int *ncol, 
	    double *table, 
	    int *ldtabl,
	    double *expect, 
	    double *percnt, 
	    double *emin, 
	    double *prt,
	    double *pre, 
	    int *workspace);

double fisher(table_t t);

#endif
