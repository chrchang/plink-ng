

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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"
#include "stats.h"

void HaploPhase::haplotypicCC(map<int,int> & tests, int nt, bool display)
{

  vector<double> caseN(nt);
  vector<double> controlN(nt);
  
  // Consider each individual
  for (int i=0; i<P.n; i++)
    {

      Individual * person = P.sample[i];

      // Case?
      if ( ! person->missing )
	{
	  
	  if ( person->aff )
	    {
	      
	      for (int z = 0 ; z < hap1[i].size(); z++)
		{
		  
		  map<int,int>::iterator i1 = tests.find(hap1[i][z]);
		  map<int,int>::iterator i2 = tests.find(hap2[i][z]);

		  if ( i1 != tests.end() )
		    {
		      if (!ambig[i]) 
			caseN[i1->second]++;
		      else
			caseN[i1->second] += pp[i][z];		      
		    }
		  
		  if ( ! ( haploid || X && person->sex ) ) 
		    {
		      if ( i2 != tests.end() )
			{
			  if (!ambig[i]) 
			    caseN[i2->second]++;
			  else
			    caseN[i2->second] += pp[i][z];
			}
		    }

		}  
	    }
	  // Or control?
	  else 
	    {

	      for (int z = 0 ; z < hap1[i].size(); z++)
		{
		  
		  map<int,int>::iterator i1 = tests.find(hap1[i][z]);
		  map<int,int>::iterator i2 = tests.find(hap2[i][z]);
		  
		  if ( i1 != tests.end() )
		    {
		      if (!ambig[i]) 
			controlN[i1->second]++;
		      else
			controlN[i1->second] += pp[i][z];
		    }
		  
		  if ( ! ( haploid || X && person->sex ) ) 
		    {
		      if ( i2 != tests.end() )
			{
			  if (!ambig[i]) 
			    controlN[i2->second]++;
			  else
			    controlN[i2->second] += pp[i][z];
			}
		    }
		  
		  
		}  	  
	    }
	}
      
    } // next individual


  ////////////////////////////////////////
  // Display header haplotype information?
  
  if ( display )
    {
      HTEST << setw(10) << hname << " ";
  
      if (nt==2)
	{
	  
	  // Find test haplotype
	  int hh=0;
	  map<int,int>::iterator i1 = tests.begin();
	  while ( i1 != tests.end() )
	    {
	      if ( i1->second == 0 )
		hh = i1->first;
	      i1++;
	    }
	  
	  HTEST << setw(12) << haplotypeName(hh) << " ";
	  HTEST << setw(10) << caseN[0] / ( caseN[0] + caseN[1] ) << " "
		<< setw(10) << controlN[0] / ( controlN[0] + controlN[1] ) << " ";	
	  
	}
      else
	{
	  
	  HTEST << setw(12) << "OMNIBUS" << " "
		<< setw(10) << "NA" << " "
		<< setw(10) << "NA" << " ";
	}
    }
  

  ///////////////////////////////////////////////////
  // Standard chi-sq statistic: use for omnibus test
  
  vector<double> rowT(nt);
  double caseT = 0;
  double controlT = 0;
  
  for (int h=0; h<nt; h++)
    {
      rowT[h] = caseN[h] + controlN[h];
      caseT += caseN[h];
      controlT += controlN[h];       
    }
  
  double chi2 = 0;
  for (int h=0; h<nt; h++)
    {
      double exp = ( rowT[h] * caseT ) / (caseT + controlT);
      chi2 += ( ( caseN[h] - exp ) * ( caseN[h] - exp ) ) / exp ; 
      
      exp = ( rowT[h] * controlT ) / (caseT + controlT);
      chi2 += ( ( controlN[h] - exp ) * ( controlN[h] - exp ) ) / exp ;  
    }
  
  

  /////////////////////////////////////////
  // Display what we have? 

  if ( display )
    {
      
      if ( realnum(chi2) )
	{
	  HTEST << setw(10) << chi2 << " "
		<< setw(4) << nt-1 << " "
		<< setw(10) << chiprobP(chi2,nt-1) << " ";
	}
      else
	{
	  HTEST << setw(10) << "NA" << " "
		<< setw(4) << "NA" << " "
		<< setw(10) << "NA" << " ";
	}
      
      
      // Display OR and CI
      //       if (par::display_ci)
      // 	HTEST << setw(12) << string("L"+dbl2str(par::ci_level*100)) << " "
      // 	      << setw(12) << string("U"+dbl2str(par::ci_level*100)) << " ";
      
      
      for (int snps=0; snps<ns-1; snps++)
	HTEST << P.locus[S[snps]]->name << "|";
      
      HTEST << P.locus[S[ns-1]]->name << "\n";
    }


  /////////////////////////////////////////////////////////
  //
  // Adjust result based on empirical variance
  //
  /////////////////////////////////////////////////////////


  if ( useEmpiricalVariance )
    {
      
      set<int> hs;
      map<int,int>::iterator i1 = tests.begin();
      while ( i1 != tests.end() )
	{
	  if ( i1->second == 0 )
	    hs.insert( i1->first);
	  ++i1;
	}
      
      // This updates empiricalVariance in HaploPhase      
      calculateEmpiricalVariance(hs);
      
    }
  
  
  // Save chi-square statistic in temporary storage
  result = chi2;


  //////////////////////////////
  // Test based on proportions

  // Z = ( p1 - p2 ) / ( sqrt(  p*(1-p)*(1/n1+1/n2) ) ) 
  
  // we get p from HaploPhase (i.e. population frequency
  // Be careful of missing phenotypes here?

  // Haplotype-specific tests only

  if ( nt == 2 )
    {

      ///////////////////////////
      // Too much of a stretch? 
      
      if ( ratio < 0.01 ) 
	{
	  result = -1;
	  pvalue = -9;
	  odds = 1;
	  case_freq = 0;
	  control_freq = 0;
	  return;
	}
      
      double p2 = caseN[0] / ( caseN[0] + caseN[1] );
      double p1 = controlN[0] / ( controlN[0] + controlN[1] );

      // Based on asymptotic variance
      
      double n2 = caseN[0] + caseN[1];
      double n1 = controlN[0] + controlN[1];
      double p = ( n1 * p1 + n2 * p2 ) / ( n1 + n2 );
      
      // Instead of ( p2 - p1 ) / ( sqrt(  p*(1-p)*(1/n1+1/n2) ) )
      // use the empirical variance of dosage

      double Ze = ( p2 - p1 ) / ( sqrt( empiricalVariance * (1/n1+1/n2) ) );
      
      // Chi-squared statistic based on test of proportions, using
      // empirical estimate of variance

      result = chi2 = Ze * Ze; 
      case_freq = p2;
      control_freq = p1;

    }
  
  pvalue = chiprobP(chi2,nt-1);
  
  //////////////////////////////////
  // Return odds ratio

  if ( nt == 2 )
    {
      odds = ( caseN[0] * controlN[1] ) 
	/ ( controlN[0] * caseN[1] );
    }
  
}




///////////////////////////////////
// Multimarker test with weighting

void HaploPhase::haplotypicWeightedCC()
{

  vector_t weights;

  for (int i=0; i<nh; i++)
    {
      map<string,double>::iterator 
	whap = new_pred_weighted_allele[current].find( haplotypeName(i) );

      if ( whap != new_pred_weighted_allele[current].end() )
	{
	  weights.push_back( whap->second );
	}
      else
	{
	  weights.push_back( 0 );
	}
    }
  

   vector<double> caseN(2);
   vector<double> controlN(2);
  
   // Consider each individual
   for (int i=0; i<P.n; i++)
     {

       Individual * person = P.sample[i];

       // Case?
       if ( ! person->missing )
	 {
	   
	   if (person->aff )
	     {
	       
	       for (int z = 0 ; z < hap1[i].size(); z++)
		 {
		   
		   int h1 = hap1[i][z];
		   int h2 = hap2[i][z];
		   
		   if (!ambig[i]) 
		     {
		       caseN[0] += weights[h1];
		       caseN[1] += 1-weights[h1];
		     }
		   else
		     {
		       caseN[0] += weights[h1] * pp[i][z];
		       caseN[1] += (1-weights[h1]) * pp[i][z];
		     }

		  if ( ! ( haploid || X && person->sex ) ) 
		    {
		      if (!ambig[i]) 
			{
			  caseN[0] += weights[h2];
			  caseN[1] += 1-weights[h2];
			}
		      else
			{
			  caseN[0] += weights[h2] * pp[i][z];
			  caseN[1] += (1-weights[h2]) * pp[i][z];
			}
		    }
		 }
	     }
	   // Or control?
	   else 
	     {
	       for (int z = 0 ; z < hap1[i].size(); z++)
		 {
		   
		   int h1 = hap1[i][z];
		   int h2 = hap2[i][z];
		   
		   if (!ambig[i]) 
		     {
		       controlN[0] += weights[h1];
		       controlN[1] += 1-weights[h1];
		     }
		   else
		     {
		       controlN[0] += weights[h1] * pp[i][z];
		       controlN[1] += (1-weights[h1]) * pp[i][z];
		     }

		   if ( ! ( haploid || X && person->sex ) ) 
		     {
		       if (!ambig[i]) 
			 {
			   controlN[0] += weights[h2];
			   controlN[1] += 1-weights[h2];
			 }
		       else
			 {
			   controlN[0] += weights[h2] * pp[i][z];
			   controlN[1] += (1-weights[h2]) * pp[i][z];
			 }
		     }

		 }
	     }
	 }

     } // next individual

   

   HTEST << setw(10) << hname << " ";
   
   // set hh to integer
   
   HTEST << setw(12) << new_map[current]->allele1 << " ";
   
   HTEST << setw(10) << caseN[0] / ( caseN[0] + caseN[1] ) << " "
	 << setw(10) << controlN[0] / ( controlN[0] + controlN[1] ) << " ";
   
  
   vector<double> rowT(2);
   double caseT = 0;
   double controlT = 0;
  
   for (int h=0; h<2; h++)
     {
       rowT[h] = caseN[h] + controlN[h];
       caseT += caseN[h];
       controlT += controlN[h];       
     }

   double chi2 = 0;
   for (int h=0; h<2; h++)
     {
       double exp = ( rowT[h] * caseT ) / (caseT + controlT);
       chi2 += ( ( caseN[h] - exp ) * ( caseN[h] - exp ) ) / exp ; 
     
       exp = ( rowT[h] * controlT ) / (caseT + controlT);
       chi2 += ( ( controlN[h] - exp ) * ( controlN[h] - exp ) ) / exp ;  
     }

  
   if ( realnum(chi2) )
     {
       HTEST << setw(10) << chi2 << " "
	     << setw(4) << 1 << " "
	     << setw(10) << chiprobP(chi2,1) << " ";
     }
   else
     {
       HTEST << setw(10) << "NA" << " "
	     << setw(4) << "NA" << " "
	     << setw(10) << "NA" << " ";
     }
   
   for (int snps=0; snps<ns-1; snps++)
     HTEST << P.locus[S[snps]]->name << "|";
    
   HTEST << P.locus[S[ns-1]]->name << "\n";
   
}

