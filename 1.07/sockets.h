

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __MODEL_H__
#define __MODEL_H__

#include<vector>
#include "plink.h"

using namespace std;

vector<string> socketConnection( Plink * P, 
				 string ip_addr , 
				 int port , 
				 string message ) ; 


#endif
