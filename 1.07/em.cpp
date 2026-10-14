

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
#include <cassert>

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "genogroup.h"
#include "phase.h"
#include "haplowindow.h"

extern ofstream LOG;
using namespace std;




////////////////////////////////////////////////
// Original, single window EM algorithm, without
// genoGrouping -- this is now used for the meta-EM


void HaploPhase::performEM_original()
  {

  vector_t uc(nh, 0); // unambigous counts
  vector_t ac(nh, 0); // ambigous counts

  // Count numbers of unambigous haplotypes
  // as these stay constant throughout EM

  for (int i=0; i<P.n; i++)
    {
      if (P.sample[i]->founder && include[i])
	{
	  if (!ambig[i])
	    {
	      uc[hap1[i][0]]++;
	      if ( ! (haploid || (X && P.sample[i]->sex)))
		uc[hap2[i][0]]++;
	    }
	}
    }


  //////////////////  
  // Begin E-M

  double sampleLogLikelihood = 0;

  for (int j=0; j<= par::haplo_plem_meta_iter; j++)
    {


      ///////////
      // E-step

      for (int i=0; i<P.n; i++)
	{
	  if (P.sample[i]->founder && include[i])
	    {
	      if (ambig[i])
		{

		  double s=0;

		  // Haploid phases...
		  if (haploid || (X && P.sample[i]->sex))
		    {
		      for (int z=0; z<hap1[i].size(); z++)
			{
			  pp[i][z] = f[hap1[i][z]];
			  s += pp[i][z];
			}
		    }
		  else // ... or diploid
		    {
		      for (int z=0; z<hap1[i].size(); z++)
			{
			  pp[i][z] = f[hap1[i][z]] * f[hap2[i][z]];
			  if (hap1[i][z] != hap2[i][z])
			    pp[i][z] *= 2;
			  s += pp[i][z];
			}
		    }

		  for (int z=0; z<hap1[i].size(); z++)
		    pp[i][z] /= s;
		}
	    }
	}


      ///////////
      // M-step

      // unambiguous counts
      for (int h=0; h<nh; h++)
	f[h] = uc[h];
      
      // then add the fractional ones
      for (int i=0; i<P.n; i++)
	if (P.sample[i]->founder && include[i])
	  if (ambig[i])
	    {
	      
	      if (haploid || (X && P.sample[i]->sex))
		{
		  for (int z=0; z<hap1[i].size(); z++)
		    f[hap1[i][z]] += pp[i][z];
		}
	      else
		for (int z=0; z<hap1[i].size(); z++)
		  {
		    f[hap1[i][z]] += pp[i][z];
		    f[hap2[i][z]] += pp[i][z];
		  }
	    }
      
      // validN is the total number of *chromosomes*
      for (int h=0; h<nh; h++)
	f[h] /= (double)validN;



      ///////////////////////
      // Update likelihood 

      if (j % par::haplo_plem_meta_likelihood_iter == 0)
	{

	  double lnl = 0;
	  
	  for (int i=0; i<P.n; i++)
	    {
	      
	      if (P.sample[i]->founder && include[i])
		{
		  double lk = 0;
		  
		  if (haploid || (X && P.sample[i]->sex))
		    {
		      for (int z=0; z<hap1[i].size(); z++)
			lk += f[hap1[i][z]];
		    }
		  else
		    for (int z=0; z<hap1[i].size(); z++)
		      {
			lk += f[hap1[i][z]] * f[hap2[i][z]];
			if (hap1[i][z] != hap2[i][z])
			  lk += f[hap1[i][z]] * f[hap2[i][z]];
		      }
		  lnl -= log(lk);
		}
	    }
	  
	  if (j > 0 && 
	      sampleLogLikelihood - lnl < par::haplo_plem_meta_tol )
	    {
	      break;
	    }
	  sampleLogLikelihood = lnl;
	}
    }
}
