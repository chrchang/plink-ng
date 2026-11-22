

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
#include <vector>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "stats.h"
#include "perm.h"



void Plink::perm_testGXE2(Perm & perm)
{

  // Assumes SNP-major mode

  if (!par::SNP_major) Ind2SNP();
  
  // This procedure is only for continuous traits
  if (par::bt)
    error("Can only use --gxe option with continuous phenotypes");
  
  // GxE test statistics
  vector<double> original;
  
  // Empirical p-valuess
  perm.setTests(nl_all);
  
  // Construct a binary covariate for GxE
  // Individuals who are missing for the 
  // covariate will already have been set 
  // to missing for the phenotype -- also 
  // allow for 0 to equal missing here 
  // (i.e. use affection status coding)

  for (int i=0; i<n; i++)
    {  
      if (sample[i]->covar == 0)
	sample[i]->missing = true;
      else if (sample[i]->covar == 2)
	sample[i]->bcovar = false;
      else 
	sample[i]->bcovar = true;
    }


  ////////////////////////////////
  // Set up permutation structure 

  perm.setPermClusters(*this);
  perm.originalOrder();
    
  
  ////////////////////////////////////
  // If we do perform permutation, 
  // check the permutation procedure here: 
  // i.e. pperson->bcovar or gperson->bcovar
  

  ////////////////////////////////////
  // Quantitative trait regression
  
  original = testQAssocGXE2(true,perm);
  


  ////////////////////////////////////
  // No permutations for now
  
  shutdown();

}     


/////////////////////////////////////////////
// Simple quantitative trait association test
// Assumes SNP-major mode


vector<double> Plink::testQAssocGXE2(bool print_results , 
				    Perm & perm )
{	
  
  vector<double> results(nl_all);

  ofstream ASC;
  if (print_results)
    
    {
      string f = par::output_file_name + ".qassoc.gxe";
      printLOG("Writing QT GxE association results to [ " + f + " ] \n");
      ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " " 
	    << setw(par::pp_maxsnp) << "SNP" << " " 
	    << setw(8) << "NMISS1" << " " 
	    << setw(10) << "BETA1" << " " 
	    << setw(10) << "SE1" << " " 
	    << setw(8) << "NMISS2" << " " 
	    << setw(10) << "BETA2" << " " 
	    << setw(10) << "SE2" << " " 
	    << setw(8) << "Z_GXE" << " " 
	    << setw(12) << "P_GXE" << " " 
	    << "\n";
	ASC.precision(4);
    }


  // Iterate over each locus

  vector<CSNP*>::iterator s = SNP.begin();
  int l = 0;
  while ( s != SNP.end() )
    {	
      
      // Skip possibly
      if (par::adaptive_perm && !perm.snp_test[l])
	{
	  // advance to next SNP
	  s++;
	  l++;
	  continue;
	}

      double g_mean1=0;
      double g_var1=0;
      double qt_mean1=0;
      double qt_var1=0;
      double qt_g_covar1=0;
      int nanal1 = 0;  

      double g_mean2=0;
      double g_var2=0;
      double qt_mean2=0;
      double qt_var2=0;
      double qt_g_covar2=0;
      int nanal2=0;
      
      
      ///////////////////////////////
      // Iterate over each individual

      vector<Individual*>::iterator person = sample.begin();
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      
      while ( person != sample.end() )
	{
	  // Permuted self
	  Individual * pperson = (*person)->pperson;
	  
	  // Genotype
	  bool s1 = *i1;
	  bool s2 = *i2;
	  
	  if (!pperson->missing)
	    {
	      if ( ! ( s1 && !s2) )     //   10 = missing 
		{
		  
		  if (pperson->bcovar)
		    qt_mean1 += pperson->phenotype;
		  else
		    qt_mean2 += pperson->phenotype;
		  
		  if ( (!s1) && (!s2) )      //   00 = hom(11)
		    {
		      if (pperson->bcovar) g_mean1+=2;
		      else g_mean2+=2;
		    }
		  else if ( (!s1) && s2)     //   01 = het(12)
		    {
		      if (pperson->bcovar) g_mean1++;
		      else g_mean2++;
		    }		  
		  
		  if (pperson->bcovar) nanal1++;
		  else nanal2++;
		}
	      
	    }
	  
	  // Next person
	  i1++;
	  i2++;
	  person++;
	}
      
      
      // Calculate mean
      qt_mean1 /= (double)nanal1;
      g_mean1 /= (double)nanal1;      
      
      qt_mean2 /= (double)nanal2;
      g_mean2 /= (double)nanal2;      


      // Iterate over individuals again
      person = sample.begin();
      i1 = (*s)->one.begin();
      i2 = (*s)->two.begin();
      
      while ( person != sample.end() )
	{
	  // Permuted self
	  Individual * pperson = (*person)->pperson;
	  
	  // Genotype
	  bool s1 = *i1;
	  bool s2 = *i2;

	  if (!pperson->missing)
	    {
	      if ( ! ( s1 && !s2) )      //   10 = missing 
		{
		  if (pperson->bcovar)
		    qt_var1 += (pperson->phenotype-qt_mean1) * ( pperson->phenotype-qt_mean1 ) ;
		  else
		    qt_var2 += (pperson->phenotype-qt_mean2) * ( pperson->phenotype-qt_mean2 ) ;

		  double g = 0;
		  
		  if ( (!s1) && (!s2) )      //   00 = hom(11)
		    g=2;
		  else if ( (!s1) && s2 )    //   01 = het(12)
		    g=1;
		  
		  if (pperson->bcovar)
		    {
		      g_var1 += (g-g_mean1) * ( g-g_mean1 ) ;
		      qt_g_covar1 += ( pperson->phenotype - qt_mean1 ) * ( g - g_mean1 ) ;
		    }
		  else
		    {
		      g_var2 += (g-g_mean2) * ( g-g_mean2 ) ;
		      qt_g_covar2 += ( pperson->phenotype - qt_mean2 ) * ( g - g_mean2 ) ;
		    }
		  
		}
	    }

	  // Next individual
	  i1++;
	  i2++;
	  person++;
	}
   

      qt_var1 /= (double)nanal1 - 1;
      g_var1 /= (double)nanal1 - 1;
      qt_g_covar1 /= (double)nanal1 - 1;
      
      qt_var2 /= (double)nanal2 - 1;
      g_var2 /= (double)nanal2 - 1;
      qt_g_covar2 /= (double)nanal2 - 1;

      double beta1 = qt_g_covar1 / g_var1;
      double vbeta1 = (qt_var1/g_var1 - (qt_g_covar1*qt_g_covar1)/(g_var1*g_var1) ) / (nanal1-2);

      double beta2 = qt_g_covar2 / g_var2;
      double vbeta2 = (qt_var2/g_var2 - (qt_g_covar2*qt_g_covar2)/(g_var2*g_var2) ) / (nanal2-2);

      double Z = (beta1-beta2) / sqrt( vbeta1 + vbeta2 ) ; 
      
      
      if (print_results)
	{
	  
	  ASC << setw(4) << locus[l]->chr  << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " ";
	  
	  if (realnum(Z)) 	   
	    {
	      ASC << setw(8) << nanal1 << " " 
		  << setw(10) << beta1 << " "
		  << setw(10) << sqrt(vbeta1) << " " 
		  << setw(8) << nanal2 << " " 
		  << setw(10) << beta2 << " " 
		  << setw(10) << sqrt(vbeta2) << " " 
		  << setw(8) << Z << " " 
		  << setw(12) << chiprobP(Z*Z,1) << "\n";
	    }
	  else
	    {
	      ASC << setw(8) << "NA" << " " 
		  << setw(10) << "NA" << " "
		  << setw(10) << "NA" << " " 
		  << setw(8) << "NA" << " " 
		  << setw(10) << "NA" << " " 
	 	  << setw(10) << "NA" << " " 
		  << setw(8) << "NA" << " " 
		  << setw(12) << "NA" << "\n";
	    }
	 	  
	}
      
      results[l] = Z;

      // Advance to next SNP
      s++;
      l++;
	
    }
  
  if (print_results)
    ASC.close();
  
  return results;
}

