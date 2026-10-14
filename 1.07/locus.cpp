

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include "plink.h"
#include "options.h"

using namespace std;

vector<Z> Plink::calcLocusIBD(Individual * p1, Individual * p2, Z I)
{

  // ** TODO ** -- development code -- lots of possible speedups here.

  // 1) Remove I as an argument
  // 2) Now no need to calculate all probabilities -- just do the three required for that 
  //    particular locus


  // Store calculated P(M|Z) here

  vector<Z> ZL(nl);

  // Locus counter

  int l = 0;
  

  //////////////////////////////
  // All SNPs in the scan region

  for (int l2=par::run_start; l2<=par::run_end; l2++)
    {
      

      /////////////////////
      // Allele frequencies

      double p = locus[l2]->freq;
      double q = 1 - p;

  
      // Na = # alleles = 2N where N is number of individuals

      double Na = locus[l]->nm; 
      double x = p * Na;
      double y = q * Na;
      

      /////////////////////////////////////////////////
      // Assign P(M|Z) based on genotype for each SNP
      
      bool a1 = p1->one[l2];
      bool a2 = p1->two[l2];

      bool b1 = p2->one[l2];
      bool b2 = p2->two[l2];
      

      // Assign unit vector if either genotype is missing

      if ( ( a1 && (!a2) ) ||
	   ( b1 && (!b2) ) )
	{
 	  ZL[l].z0 = ZL[l].z1 = ZL[l].z2 = 1;
	  l++;
	  continue;
	}

      
      if ( a1 ) 
	{
	  if ( a2 ) 
	    {
	      if ( b1 ) 
		{
		  if ( b2 ) 
		    {
		      // aa, aa
		      
		      ZL[l].z0 = q*q*q*q    
			* ( (y-1)/y * (y-2)/y * (y-3)/y * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = q*q*q   * ( (y-1)/y * (y-2)/y *  Na/(Na-1) * Na/(Na-2));
		      ZL[l].z2 = q*q * ( (y-1)/y *  Na/(Na-1));
		      
		    }
		}
	      else
		{
		  if ( b2 ) 
		    {

		      // aa, Aa

		      ZL[l].z0 = 2*p*q*q*q  
			* ( (y-1)/y * (y-2)/y *  (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = p*q*q   *  ((y-1)/y  * Na/(Na-1) * Na/(Na-2));
		      ZL[l].z2 = 0;      
		      
		    }
		  else
		    {
		      
		      // aa, AA
		      
		      ZL[l].z0 = p*p*q*q    
			* ( (x-1)/x * (y-1)/y * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = 0;  
		      ZL[l].z2 = 0;
		      
		    }
		}
	    }
	}
      else
	{
	  if ( a2 ) 
	    {
	      if ( b1 ) 
		{
		  if ( b2 ) 
		    {

		      // Aa, aa
		      		      
		      ZL[l].z0 = 2*p*q*q*q  
			* ( (y-1)/y * (y-2)/y *  (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = p*q*q   *  ((y-1)/y  * Na/(Na-1) * Na/(Na-2));
		      ZL[l].z2 = 0;      
		      
		    }
		}
	      else
		{
		  if ( b2 ) 
		    {

		      // Aa, Aa

		      ZL[l].z0 = 4*p*p*q*q  
			* ( (x-1)/x * (y-1)/y * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = p*p*q   
			* ( (x-1)/x *  Na/(Na-1) * Na/(Na-2) ) 
			+ p*q*q 
			* ((y-1)/y  * Na/(Na-1) * Na/(Na-2));
		      ZL[l].z2 = 2*p*q *  Na/(Na-1) ;

		    }
		  else
		    {
		      
		      // Aa, AA

		      ZL[l].z0 = 2*p*p*p*q  
			* ( (x-1)/x * (x-2)/x *  (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = p*p*q  
			* ( (x-1)/x *  Na/(Na-1) * Na/(Na-2) );
		      ZL[l].z2 = 0;

		    }
		}
	    }
	  else
	    {
	      if ( b1 ) 
		{
		  if ( b2 ) 
		    {

		      // AA, aa
		      		      
		      ZL[l].z0 = p*p*q*q    
			* ( (x-1)/x * (y-1)/y * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = 0;  
		      ZL[l].z2 = 0;

		    }
		}
	      else
		{
		  if ( b2 ) 
		    {

		      // AA, Aa

		      ZL[l].z0 = 2*p*p*p*q  
			* ( (x-1)/x * (x-2)/x *  (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
		      ZL[l].z1 = p*p*q  
			* ( (x-1)/x *  Na/(Na-1) * Na/(Na-2) );
		      ZL[l].z2 = 0;

		    }
		  else
		    {

		      // AA, AA

		      ZL[l].z0 = p*p*p*p 
			* ( (x-1)/x * (x-2)/x * (x-3)/x * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)));
		      ZL[l].z1 = p*p*p  
			* ( (x-1)/x * (x-2)/x *  Na/(Na-1) * Na/(Na-2));
		      ZL[l].z2 = p*p 
			* ( (x-1)/x *  Na/(Na-1));

		    }
		}
	    }
	}
    

      /////////////////////////////////////
      // Fudge factor for genotyping error

        if ( ZL[l].z1 < 0.0001 ) 
  	{
  	  ZL[l].z1 = 0.0001;
  	  double S = ZL[l].z0 + ZL[l].z1 + ZL[l].z2;
  	  ZL[l].z0 /= S;
  	  ZL[l].z1 /= S;
  	  ZL[l].z2 /= S;
  	}      
      


      /////////////////
      // Next SNP

      l++;

    }

  return ZL;
}

      


//       ////////////////////////////////////////////
//       // 2. Allow for possible genotyping error
//       //    sum_{all possible true genotypes} P(observed G|true G) P(true G |IBD)

// //       double e = 0.005;
// //       double f = 0.001;

// //       double ER_AA_AA__AA_AA = 1- 2e - 2f - 2e*f - e*e - f*f;
// //       double ER_AA_AB__AA_AA = e;
// //       double ER_AA_BB__AA_AA = f;
	
// //       double ER_AB_AA__AA_AA = e;
// //       double ER_AB_AB__AA_AA = e*e;
// //       double ER_AB_BB__AA_AA = e*f;
      
// //       double ER_BB_AA__AA_AA = f;
// //       double ER_BB_AB__AA_AA = e*f;
// //       double ER_BB_BB__AA_AA = f*f;

      
// //       double ER_AA_AA__AA_AB = e;
// //       double ER_AA_AB__AA_AB = 1 - 3*e - 2*e*e - 2*e*f - f;
// //       double ER_AA_BB__AA_AB = e;

// //       double ER_AB_AA__AA_AB = e*e;
// //       double ER_AB_AB__AA_AB = e;
// //       double ER_AB_BB__AA_AB = e*e;
      
// //       double ER_BB_AA__AA_AB = f*e;
// //       double ER_BB_AB__AA_AB = f;
// //       double ER_BB_BB__AA_AB = f*e;

//       // Sum over all possible true genotypes P(Observed G | True G) P(True G |IBD = 1)

