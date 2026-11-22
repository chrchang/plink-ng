
//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2007 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <assert.h>

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"


#include "genogroup.h"
#include "haplowindow.h"


/////////////////////////////////////////////////////////////////////
// For a given window, collapse all genotypes into a unique groups,
// 'genoGroups' and perform subsequent EM on these entitities rather
// than on individuals


void HaploWindow::enumerateGenogroups()
{

  // Consider each individual

  for (int i=0; i < P->n ; i++) 
    {
            
      // Only phase non-missing founders
      
      if ( ! (P->sample[i]->founder && haplo->include[i]))
	continue;

      // Build a new multilocus genotype set
      
      MultiLocusGenotype * m = new MultiLocusGenotype;
      
      
      // Include sex here for X chr SNPs
      
      if ( haplo->X )
    	  m->g.push_back(P->sample[i]->sex);
      
      // Genotypes

      for (int s=0; s<ns; s++) 
	{

	  bool s1 = par::SNP_major ? 
	    P->SNP[ S[s] ]->one[i] :
	    P->sample[i]->one[ S[s] ];

	  bool s2 = par::SNP_major ? 
	    P->SNP[ S[s] ]->two[i] :
	    P->sample[i]->two[ S[s] ];
	  
	  m->g.push_back(s1);
	  m->g.push_back(s2);
	}


      // One individual, this individual

      m->count = 1;
      m->reference = i;
            

      // But have we already seen a similar genoGroup?

      set<MultiLocusGenotype*>::iterator im = genotypes.find(m);
      if (im == genotypes.end() ) 
	{
	  genoGroup[i] = m;
	  genotypes.insert( m );
	} 
      else 
	{
	  delete m;
	  (*im)->count++;
	  genoGroup[i] = *im;
	}

    } // Next individual
  

}

