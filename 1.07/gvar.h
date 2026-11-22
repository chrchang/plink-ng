

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __GVAR_H__
#define __GVAR_H__

#include <vector>
#include <set>
#include "plink.h"

using namespace std;

class Variant : public Locus {
  
 public:

  int nallele;
  vector<string> alleles;
  vector_t freqs;
  map<string,int> acode;
  
  bool allelicVariation;
  bool copyNumberVariation;
  bool integerDosage;
  Variant() 
    { 
      nallele = 0;
      integerDosage = true;
      allelicVariation = copyNumberVariation = false;
    }
    
};

class GVariant {  
 public:

  GVariant()
    {
      missing = true;
      allele1 = allele2 = -1;
      dosage1 = dosage2 = 0;
    }
  bool missing;
  int allele1;
  int allele2;
  float dosage1;
  float dosage2;  
};

#endif
