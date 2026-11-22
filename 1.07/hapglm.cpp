

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include "whap.h"
#include "helper.h"
#include "plink.h"
#include "options.h"
#include "perm.h"
#include "nlist.h"
#include "phase.h"
#include "model.h"
#include "linear.h"
#include "logistic.h"
#include "stats.h"




//////////////////////////////////////////////////////////////
// Implements --hap-logistic and --hap-linear functions
// Use framework provided by --chap/whap.cpp
// can perform either omnibus or haplotype specific tests


vector_t Plink::glmHaplotypeTest(bool print, Perm & perm)
{
  
  ///////////////////////////////////////////////
  //                                           //
  // Some basic setup first                    //
  //                                           //
  ///////////////////////////////////////////////


  // Use basic GLM function to fit linear and logistic 
  // models: although, let it know that there will not 
  // be a 'main' SNP 

  par::assoc_glm_without_main_snp = true;


  // Return a single result

  vector_t results;


  // Haplotypes at this position have already been phased 


  // Record the number of common haplotypes
  
  int nch = 0;  
  set<int> commonHaplotypes;
  for (int h=0; h < haplo->nh; h++)
    if ( haplo->f[h] >= par::min_hf )
      {
	++nch;
	commonHaplotypes.insert(h);
      }


//   if ( ! par::test_hap_GLM_omnibus ) 
//     haplo->HTEST << setw( haplo->ns + 1 ) << haplo->haplotypeName(0) << " ";

    
  if ( nch < 2 ) 
    {
      haplo->HTEST << setw(4) << haplo->ns << " " 
		   << setw(4) << nch << " " 
		   << setw(4) << locus[haplo->S[0]]->chr << " " 
		   << setw(12) << locus[haplo->S[0]]->bp << " " 
		   << setw(12) << locus[haplo->S[haplo->ns-1]]->bp << " " 
		   << setw(par::pp_maxsnp) << locus[haplo->S[0]]->name << " " 
		   << setw(par::pp_maxsnp) << locus[haplo->S[haplo->ns-1]]->name << " ";
      
      if ( ! par::test_hap_GLM_omnibus )
	{

	  if ( nch==1 )
	    haplo->HTEST << setw(12) << haplo->haplotypeName(0) << " ";
	  else
	    haplo->HTEST << setw(12) << "NA" << " ";

	  haplo->HTEST << setw(8) << "NA" << " " 
		       << setw(8) << "NA" << " " 
		       << setw(8) << "NA" << " " 
		       << setw(8) << "NA" << "\n";
	}
      else
	haplo->HTEST << setw(8) << "NA" << " " 
		     << setw(8) << "NA" << "\n";

      results.push_back( 0 );
      return results;
    }
  

  // Single SNP association
  
  if ( par::test_hap_GLM_omnibus )
    {

      haplo->HTEST << setw(4) << haplo->ns << " " 
		   << setw(4) << nch << " " 
		   << setw(4) << locus[haplo->S[0]]->chr << " " 
		   << setw(12) << locus[haplo->S[0]]->bp << " " 
		   << setw(12) << locus[haplo->S[haplo->ns-1]]->bp << " " 
		   << setw(par::pp_maxsnp) << locus[haplo->S[0]]->name << " " 
		   << setw(par::pp_maxsnp) << locus[haplo->S[haplo->ns-1]]->name << " ";


      // H-1 omnibus (H0 is ref.)

      haplo->sets.clear();

      set<int>::iterator i = commonHaplotypes.begin();

      // Skip first haplotype (this is reference)
      // All rare haplotypes will therefore be 
      // lumped in with the reference
      ++i;

      while ( i != commonHaplotypes.end() )
	{
	  haplo->sets.insert(*i);
	  ++i;
	}

      
      // Fit model
      glmAssoc(false,*pperm);


      // Report results
      haplo->result = model->isValid() ? model->getStatistic() : 0;
      haplo->pvalue = par::bt ? 
	chiprobP(haplo->result,1) : ((LinearModel*)model)->getPValue();


      // Calculate omnibus tests of H-1 terms
      // Assumes the terms are: e.g. for 4 haplotypes

      // 0 intercept

      // 1 haplotype 2 of H
      // 2 haplotype 3 of H
      // 3 haplotype 4 of H

      int df = nch-1;
      vector_t h;
      h.resize(df,0);
      matrix_t H;
      sizeMatrix(H,df,model->getNP());
      for (int j=0; j<df; j++)
	H[j][j+1] = 1;
      
      double chisq = model->isValid() ? model->linearHypothesis(H,h) : 0;
      double pvalue = chiprobP(chisq,df);
      
      if ( model->isValid() )
	{
	  haplo->HTEST << setw(8) << chisq << " " 
		       << setw(8) << pvalue << "\n";
	}
      else
	{
	  haplo->HTEST << setw(8) << "NA" << " " 
		       << setw(8) << "NA" << "\n";
	}

      // Clean up
      delete model;
      
      // Return 1-p, as will be different DF for different windows
      results.push_back( 1 - pvalue );      
      return results;
    }
  
  // Otherwise, we are performing H haplotype specific tests

  set<int>::iterator i = commonHaplotypes.begin();
  while ( i != commonHaplotypes.end() )
    {
      
      haplo->sets.clear();
      haplo->sets.insert(*i);
      
      // Fit model
      glmAssoc(false,*pperm);
      
      // Report results
      vector_t coef = model->getCoefs();

      // Note: the different direction of OR
      haplo->odds = par::bt ? exp(coef[1]) : coef[1];
      haplo->result = model->isValid() ? model->getStatistic() : 0;
      haplo->pvalue = par::bt ? chiprobP(haplo->result,1) : ((LinearModel*)model)->getPValue();
	  
      // Calculate omnibus tests of H-1 terms
	  
      haplo->HTEST << setw(4) << haplo->ns << " " 
		   << setw(4) << nch << " " 
		   << setw(4) << locus[haplo->S[0]]->chr << " " 
		   << setw(12) << locus[haplo->S[0]]->bp << " " 
		   << setw(12) << locus[haplo->S[haplo->ns-1]]->bp << " " 
		   << setw(par::pp_maxsnp) << locus[haplo->S[0]]->name << " " 
		   << setw(par::pp_maxsnp) << locus[haplo->S[haplo->ns-1]]->name << " ";
      
      haplo->HTEST << setw(12) << haplo->haplotypeName(*i) << " "
		   << setw(8) << haplo->f[*i] << " ";
      
      if ( model->isValid() )
	{
	  haplo->HTEST << setw(8) << haplo->odds << " " 
		       << setw(8) << haplo->result << " " 
		       << setw(8) << haplo->pvalue << "\n";
	}
      else
	{
	  haplo->HTEST << setw(8) << "NA" << " " 
		       << setw(8) << "NA" << " " 
		       << setw(8) << "NA" << "\n";
	}

      
      // Clean up
      delete model;

      // Return chi-sq (always 1df)
      results.push_back( haplo->result );

      // Next common haplotype
      ++i;
    }
  
  return results;
  
}

