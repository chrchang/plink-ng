

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

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"
#include "stats.h"


void HaploPhase::haplotypicQTL(map<int,int> & tests, 
			       int nt, 
			       bool display_results )
{
  
  // Quantitative trait test based on a test vector; like TDT, assumes
  // only ever two groups, for now
  
  // No implementation of QTL omnibus test yet
  
  if (nt!=2) return;

  
  // Genotypic and phenotype mean, variance, covariance

  double genotypic_mean = 0;
  double genotypic_variance = 0;
  double qt_mean = 0;
  double qt_variance = 0;
  double covariance = 0;

  // number of individuals in analysis
  
  int numberIndividuals = 0;


  // For this, we will make the male X coding equivalent to
  // --xchr-model 1; except we will not add a covariate for 
  // sex; 

  //  Females:  0, 1, 2
  //  Males:    0,    1
        

  /////////////////////////////
  // Iterate over individuals
  
  for (int i = 0 ; i < P.n; i++)
    {
      
      if ( hap1[i].size() == 0 )
	continue;

      Individual * pperson = P.sample[i]->pperson;
      Individual * gperson = P.sample[i];

      if (!pperson->missing)
	{
	  
	  qt_mean += pperson->phenotype;
	  
	  
	  // Consider all possible phases
	  for (int z = 0 ; z < hap1[i].size(); z++)
	    {
	      
	      map<int,int>::iterator i1 = tests.find(hap1[i][z]);
	      map<int,int>::iterator i2 = tests.find(hap2[i][z]);
	      
	      // i1 and i2 should always point to a 0/1 variable; 
	      // but the coding is reversed (for god knows what reason)
	      // such as convention means the to-be-tested variant(s)
	      // have a 0; therefore reverse here.

	      int c1 = 1 - i1->second;
	      int c2 = 1 - i2->second;
	      
	      if ( i1 != tests.end() )
		{
		  if (!ambig[i]) 
		    genotypic_mean += c1;
		  else
		    genotypic_mean += c1 * pp[i][z];
		}
	      
	      if ( ! ( haploid || ( X && gperson->sex ) ) )
		{
		  if ( i2 != tests.end() )
		    {
		      if (!ambig[i]) 
			genotypic_mean += c2;
		      else
			genotypic_mean += c2 * pp[i][z];
		    }
		}
	    }
	  
	  numberIndividuals++;
	}

    } // Next individual

      
  qt_mean /= (double)numberIndividuals;
  genotypic_mean /= (double)numberIndividuals;      
      
  
  //////////////////////////////////
  // Iterate over individuals again
        
  for (int i=0; i< P.n; i++)
    {
      
      if ( hap1[i].size() == 0 )
	continue;
	  
      Individual * pperson = P.sample[i]->pperson;
      Individual * gperson = P.sample[i];
      
      if (!pperson->missing)
	{
	  	  
	  double g = 0;
	  
	  // Consider all possible phases
	  for (int z = 0 ; z < hap1[i].size(); z++)
	    {
	      
	      map<int,int>::iterator i1 = tests.find(hap1[i][z]);
	      map<int,int>::iterator i2 = tests.find(hap2[i][z]);
	      
	      // i1 and i2 should always point to a 0/1 variable; 
	      // but the coding is reversed (for god knows what reason)
	      // such as convention means the to-be-tested variant(s)
	      // have a 0; therefore reverse here.

	      int c1 = 1 - i1->second;
	      int c2 = 1 - i2->second;

	      if ( i1 != tests.end() )
		{
		  if (!ambig[i]) 
		    g += c1;
		  else
		    g += c1 * pp[i][z];
		}
	      
	      if ( ! ( haploid || ( X && gperson->sex ) ) )
		{
		  if ( i2 != tests.end() )
		    {
		      if (!ambig[i]) 
			g += c2;
		      else
			g += c2 * pp[i][z];
		    }
		}
	    }

	  qt_variance += (pperson->phenotype-qt_mean) 
	    * ( pperson->phenotype-qt_mean ) ;
	  
	  genotypic_variance += (g-genotypic_mean) 
	    * ( g-genotypic_mean ) ;
	  
	  covariance += ( pperson->phenotype - qt_mean ) 
	    * ( g - genotypic_mean ) ;
	  	  
	}
            

    } // Next individual
  

  // Statistics

  qt_variance /= (double)numberIndividuals - 1;
  genotypic_variance /= (double)numberIndividuals - 1;
  covariance /= (double)numberIndividuals - 1;
  
  // Test statistic
  
  double beta = covariance / genotypic_variance;
  
  double vbeta = ( qt_variance/genotypic_variance 
		   - (covariance * covariance ) 
		   / (genotypic_variance* genotypic_variance) 
		   ) / (numberIndividuals-2);
  
  double t = beta / sqrt(vbeta);

  double t_p = pT(t,numberIndividuals-2);  

  // Display results?

  if ( display_results ) 
	{
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && ( t_p > par::pfvalue || t_p < 0 ) ) 
	    goto skip_p2;
	  
	  double r2 =   (covariance * covariance ) 
	    / ( qt_variance * genotypic_variance ) ;
	  	  
	  
	  HTEST << setw(10) << hname << " ";
  
	  // Find test haplotype (assuming there is a single one; 
	  // otherwise we won't be in display mode, i.e. proxy
	  // association has it's own display)

	  int hh=0;
	  
	  map<int,int>::iterator i1 = tests.begin();
	  while ( i1 != tests.end() )
	    {
	      if ( i1->second == 0 )
		hh = i1->first;
	      i1++;
	    }
	  
	  HTEST << setw(12) << haplotypeName(hh) << " ";
	  
	  HTEST << setw(8) << numberIndividuals << " " 
		<< setw(10) << beta << " " 
		<< setw(10) << r2 << " " ;

	  if (t_p >= 0) 
 	    HTEST << setw(8) << t  << " " << setw(12) << t_p << " " ;
 	  else
 	    HTEST << setw(8) << "NA"  << " " << setw(12) << "NA"  << " " ;

	  // Display SNPs

	  for (int snps=0; snps<ns-1; snps++)
	    HTEST << P.locus[S[snps]]->name << "|";
	  
	  HTEST << P.locus[S[ns-1]]->name << "\n";
	}
  
 skip_p2:
  
  // Store chi-sq and regression coefficient

  result = t;
  pvalue = t_p;
  odds = beta;

}
