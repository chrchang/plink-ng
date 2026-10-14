

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
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>

#include "plink.h"
#include "fisher.h"
#include "stats.h"
#include "helper.h"
#include "options.h"
#include "crandom.h"
#include "sets.h"
#include "perm.h"
#include "phase.h"

using namespace std;

extern ofstream LOG;


////////////////////////////////////////
// Standard 2x2 allelic association tests
// Genotypic, and quantitative trait association
// Permutation (within-cluster)

void Plink::calcAssociationWithPermutation(Perm & perm)
{
  
  
  // SNP-major mode analyses?

  if (par::assoc_glm)
    {
      if (par::SNP_major) 
	SNP2Ind();
    }
  else if (!par::SNP_major) 
    Ind2SNP();
  
  


  //////////////////////////////
  // Profile-based set test?

//   if ( par::set_score )
//     pS->profileTestInitialise();


  //////////////////////////////
  // LD-clump within each set?

  if ( par::set_test && par::set_r2 )
    {
      printLOG("Performing LD-based set test, with parameters:\n");
      printLOG("     r-squared  (--set-r2)   = " + dbl2str( par::set_r2_val ) + "\n" );
      printLOG("     p-value    (--set-p)    = " + dbl2str( chiprobP(par::set_chisq_threshold,1) ) + "\n" );
      printLOG("     max # SNPs (--set-max)  = " + int2str( par::set_max ) + "\n" );
	
      pS->makeLDSets();

    }

  //////////////////////////////
  // Step-wise set tests
  
  if ( par::set_step )
    {
      vector_t r = pS->fitStepwiseModel();
      shutdown();
    }

  // Basic association testing results
  vector<double> original;

  //////////////////////
  // Empirical p-values

  ////////////////////
  // Number of tests
  
  int ntests = nl_all;

  if ( par::set_test && ( par::set_r2 || par::set_score ) ) 
    ntests = pS->snpset.size();

  perm.setTests(ntests);


  // Case/control missingness test statistics
  vector<double> missing;

  // Observed marginals
  int aff;
  int unf;

  // Odds ratio
  vector<double> odds(nl_all);



  //////////////////////////////////////////
  // Working vectors for assoc_test_alt_perm

  vector<int> a1;
  vector<int> a2;
  vector<int> a0;

  // Expected values for the 2x2 test
  vector<double> exp_afffreq1;
  vector<double> exp_afffreq2;
  vector<double> exp_unffreq1;
  vector<double> exp_unffreq2;
  
  if (par::assoc_test_alt_perm)
    {
      a1.resize(nl_all);
      a2.resize(nl_all);
      a0.resize(nl_all);
      exp_afffreq1.resize(nl_all);
      exp_afffreq2.resize(nl_all);
      exp_unffreq1.resize(nl_all);
      exp_unffreq2.resize(nl_all);
    }


  ////////////////////////////////
  // Set up permutation structure 
  // (we need to perform this step
  //  whether or not we also 
  //  subsequently permute)

  perm.setPermClusters(*this);
  perm.originalOrder();
    


  ////////////////////////////////////////////////////////
  // Perform a test of missingness by case/control status
  
  if (par::test_missing && par::bt)
    original = testMiss(perm,true);


  //////////////////////////////////////////
  // Calculate original association results

  else if (par::bt)
    {
      
      if (par::full_model_assoc) 
	{
	  
	  ////////////////////////////////////
	  // Full model case/control test
	  
	  original = fullModelAssoc(true,perm);  
	}
      else if (par::assoc_glm)
	{
	  ///////////////////////////
	  // Logistic regression test
	  
	  original = glmAssoc(true,perm);

	}
      else if (par::CMH_test_1)
	{
	  
	  //////////////////////////////////////////
	  // 2 x 2 x K Cochran-Mantel-Haenszel test

	  original = calcMantelHaenszel_2x2xK(perm, true);

	}
      else if ( par::test_hap_GLM )
	{
	  
	  ///////////////////////////////////
	  // Haplotypic GLM tests (logistic)
	  string f = par::output_file_name + ".assoc.hap.logistic";
	  printLOG("Writing haplotype results to [ " + f + " ]\n");

	  haplo->HTEST.open(f.c_str(), ios::out);
	  haplo->HTEST.precision(3);	  
	  haplo->HTEST << setw(4) << "NSNP" << " " 
	       << setw(4) << "NHAP" << " " 
	       << setw(4) << "CHR" << " " 
	       << setw(12) << "BP1" << " " 
	       << setw(12) << "BP2" << " " 
	       << setw(par::pp_maxsnp) << "SNP1" << " " 
	       << setw(par::pp_maxsnp) << "SNP2" << " ";
	  
	  if ( ! par::test_hap_GLM_omnibus )
	    {
	      haplo->HTEST << setw(12) << "HAPLOTYPE" << " "
			   << setw(8) << "F" << " "
			   << setw(8) << "OR" << " ";
	    }
	  
	  haplo->HTEST << setw(8) << "STAT" << " " 
		       << setw(8) << "P" << "\n";

	  original = haplo->phaseAllHaplotypes(true,perm);

	  haplo->HTEST.close();

	}
      else
	{		
	  
	  //////////////////////////////////////
	  // Standard allelic case/control test
	  
	  original = testAssoc(aff,unf,
			       a1,a2,a0,
			       odds,	 
			       exp_afffreq1, exp_afffreq2,
			       exp_unffreq1, exp_unffreq2,
			       perm, 
			       true);
	  	  
	}	
  
    }
  else if (par::qt)
    {

      if (par::assoc_glm)
	{
	  ///////////////////////////
	  // Linear regression test
	  
	  original = glmAssoc(true,perm);
	}
      else if ( par::test_hap_GLM )
	{

	  ///////////////////////////
	  // Haplotypic GLM tests (linear)
	  string f = par::output_file_name + ".assoc.hap.linear";
	  printLOG("Writing haplotype results to [ " + f + " ]\n");
	  haplo->HTEST.open(f.c_str(), ios::out);
	  haplo->HTEST.precision(3);	  
	  haplo->HTEST << setw(4) << "NSNP" << " " 
	       << setw(4) << "NHAP" << " " 
	       << setw(4) << "CHR" << " " 
	       << setw(12) << "BP1" << " " 
	       << setw(12) << "BP2" << " " 
	       << setw(par::pp_maxsnp) << "SNP1" << " " 
	       << setw(par::pp_maxsnp) << "SNP2" << " ";

	  if ( ! par::test_hap_GLM_omnibus )
	    {
	      haplo->HTEST << setw(12) << "HAPLOTYPE" << " "
			   << setw(8) << "F" << " "
			   << setw(8) << "BETA" << " ";
	    }
	  
	  haplo->HTEST << setw(8) << "STAT" << " " 
		       << setw(8) << "P" << "\n";

	  original = haplo->phaseAllHaplotypes(true,perm);
	  haplo->HTEST.close();

	}
      else
	{
	  ////////////////////////////////////
	  // Quantitative trait regression
	  
	  original = testQAssoc(true,perm);
	}
      
    }
  

  // If we didn't know how many values to expect back, 
  // resize now (i.e. from haplotype tests)

  if ( par::test_hap_GLM )
    {
      ntests = original.size();      
      perm.setTests(ntests);
    }


  ////////////////////////////  
  // Display corrected p-values?
  
  string f0 = ".assoc";
  if (par::test_missing && par::bt)
    f0 = ".missing";
  else if (par::qt && !par::assoc_glm) 
    f0 = ".qassoc";
  else if (par::bt && par::assoc_glm)
    f0 += ".logistic";
  else if (par::qt && par::assoc_glm)
    f0 += ".linear";
  else if (par::CMH_test_1)
    f0 = ".cmh";
  else if (par::test_hap_GLM)
    f0 += par::bt ? ".hap.logistic" : ".hap.linear";
  else if (par::full_model_assoc) 
    {
      if (par::model_perm_best && par::permute)
	{
	  f0 = ".model.best"; 
	  printLOG("Using BEST of ALLELIC, DOM and REC for --model permutation\n");	  
	}
      else if (par::model_perm_gen && par::permute)
	{f0 = ".model.gen"; printLOG("Using GENO for --model permutation\n"); }
      else if (par::model_perm_dom)
	{f0 = ".model.dom"; printLOG("Using DOM for --model permutation/adjusted p-values\n"); }
      else if (par::model_perm_rec)
	{f0 = ".model.rec"; printLOG("Using REC for --model permutation/adjusted p-values\n"); }
      else if (par::model_perm_trend)
	{f0 = ".model.trend"; printLOG("Using CA-trend test  for --model permutation/adjusted p-values\n"); }
    }
  

  if (par::fisher_test )
    f0 += ".fisher";
  
  // Profile-based set-test?
  if ( par::set_test && par::set_r2 )
    f0 += ".set";
  else if ( par::set_test && par::set_score ) 
    f0 += ".set.score";

  // Assumes we have 1 df chi-square statistics
  
  if (par::multtest)
    {
      vector<double> obp(0);
      
      if ( par::fisher_test )
	{
	  for (int l=0; l<ntests;l++)
	    {
	      double chi = original[l] > 0 && realnum(original[l]) ? inverse_chiprob( 1-original[l] , 1) : -9;
	      obp.push_back( chi ) ;
	    }
	}
      else
	{
	  for (int l=0; l<ntests;l++)
	    obp.push_back(original[l]);
	}

      multcomp(obp,f0);
    }
  

  ////////////////////////////
  // If no permutation, then 
  // leave now

  if (!par::permute) 
    return;
 

  ////////////////////////
  // Score original sets

  vector<int> setsigsize;

  if (par::set_test) 
    {
      if ( par::set_r2 || par::set_score )
	{
	  // Score...
	  
	  // original = pS->profileTestScore();
	  
	  original = pS->fitLDSetTest(original,true);
	  
	  // ...and save # of significant SNPs	  
	  setsigsize.clear();
	  for (int i=0; i<pS->profileSNPs.size(); i++)
	    setsigsize.push_back( pS->s_min[i] );
	}
      else
	pS->cumulativeSetSum_WITHLABELS(*this,original);
    }


  /////////////////////////////
  // Ordered/rank permutation?

  if (par::mperm_rank)
    perm.setOriginalRanking(original);


  //////////////////////
  // Verbose dumping?

  if (par::mperm_save_all)
    printLOG("Dumping all permutation statistics to [ "
	     + par::output_file_name+".mperm.dump.all ]\n");
  else if (par::mperm_save_best)
    printLOG("Dumping best permutation statistics to [ "
	     + par::output_file_name+".mperm.dump.best ]\n");
  


  //////////////////////
  // Begin permutations

  bool finished = par::replicates == 0 ? true : false;

  while(!finished)
    {
      
      // Store permuted results
      
      vector<double> pr(ntests);

      
      if (par::perm_genedrop)
	{
	  if (par::perm_genedrop_and_swap)
	    perm.permuteInCluster();
	  perm.geneDrop();
	}
      else
	perm.permuteInCluster();

      if (par::test_missing)
	pr = testMiss(perm,false);
      else if ((!par::assoc_test_alt_perm) 
	       || par::qt 
	       || par::full_model_assoc 
	       || par::CMH_test_1
	       || par::assoc_glm)
	{
	  
    	  if (par::qt)
	    {
	      if (par::assoc_glm)
		pr = glmAssoc(false,perm);
	      else if ( par::test_hap_GLM )
		pr = haplo->phaseAllHaplotypes(false,perm);
	      else
		pr = testQAssoc(false,perm);
	      
	    }
	  else if (par::full_model_assoc)
	    pr = fullModelAssoc(false,perm);
	  else if (par::assoc_glm)
	    pr = glmAssoc(false,perm);
	  else if (par::CMH_test_1)
	    pr = calcMantelHaenszel_2x2xK(perm, false);
	  else if ( par::test_hap_GLM )
	    pr = haplo->phaseAllHaplotypes(false,perm);
	  else 
	    pr = testAssoc(aff,unf,
			   a1,a2,a0,
			   odds,	 
			   exp_afffreq1, exp_afffreq2,
			   exp_unffreq1, exp_unffreq2,
			   perm,
			   false);
	}
      else
	{
	  
	  /////////////////////////
	  // For binary traits only
	  
	  //  -------------
	  //  | A | B | E |  aff
	  //  ------------- 
	  //  | C | D | F |  unf
	  //  -------------
	  //   a1  a2  a0
	  
	  // a1 most likely to be common, followed by a2, then a0
	  // save aff+unf (do not alter by locus)
	  // and a1,a2,a0 marginals (which do alter by locus) 
	  // then we only need count A and B in each subsequent replicate:
	  // 
	  
	  int A, B, C, D, M;

	  ///////////////////////////////
	  // Iterate over SNPs
	  
	  vector<CSNP*>::iterator s = SNP.begin();
	  int l=0;
	  
	  while ( s != SNP.end() )
	    {	
	      
	      // In adaptive mode, possibly skip this test
	      if (par::adaptive_perm && (!perm.snp_test[l])) 
		{
		  s++;
		  l++;
		  continue;
		}
	      
	      
	      ///////////////// 	      
	      // clear counts
	      
	      D=M=0;
	      
	      ///////////////// 	      
	      // Autosomal or haploid?
	      
	      bool X=false, haploid=false;
	      if (par::chr_sex[locus[l]->chr]) X=true;
	      else if (par::chr_haploid[locus[l]->chr]) haploid=true;
	      
	      /////////////////////////////
	      // Iterate over individuals
	      
	      vector<bool>::iterator i1 = (*s)->one.begin();
	      vector<bool>::iterator i2 = (*s)->two.begin();
	      vector<Individual*>::iterator gperson = sample.begin();
	      
	      while ( gperson != sample.end() )
		{
		  
		  // Phenotype for this person (i.e. might be permuted)
		  Individual * pperson = (*gperson)->pperson;
		  
		  // SNP alleles
		  
		  bool s1 = *i1;
		  bool s2 = *i2;
		  
		  if ( ! pperson->missing )
		    {
		      if (! pperson->aff ) // unaffected 
			{	      
			  if ( haploid || ( X && (*gperson)->sex ) ) 
			    { 
			      if ( s2 )
				D++;               // (hemi, one count)
			      else if ( s1 ) M++;  // (missing, one count)
			    }
			  else
			    {
			      if ( s2 )
				{
				  if (!s1) D++;    // (het, one A count)
				  else D+=2;       // (hom, two B count)
				}
			      else if ( s1 ) M+=2; // (missing, two B count)
			    }
			}
		    }

		  // Next individual
		  gperson++;
		  i1++;
		  i2++;
		  
		}
	      
	      
	      // reconstruct rest of 2x2 table
	      C = unf - D - M;
	      A = a1[l] - C;
	      B = a2[l] - D;
	      
	      pr[l] = ( (A - exp_afffreq1[l]) * ( A - exp_afffreq1[l] ) ) / exp_afffreq1[l] 
		+ ( (C - exp_unffreq1[l]) * ( C - exp_unffreq1[l] ) ) / exp_unffreq1[l] 
		+ ( (B - exp_afffreq2[l]) * ( B - exp_afffreq2[l] ) ) / exp_afffreq2[l] 
		+ ( (D - exp_unffreq2[l]) * ( D - exp_unffreq2[l] ) ) / exp_unffreq2[l];
	      

	      // Next SNP
	      s++;
	      l++;

	    }
	}
      
      
      //////////////////////
      // Make sets?
      
      if (par::set_test) 
	{
	  if ( par::set_r2 )
	    pr = pS->fitLDSetTest(pr,false);
	  else if ( par::set_score ) 
	    pr = pS->profileTestScore();
	  else
	    pS->cumulativeSetSum_WITHOUTLABELS(pr,perm.current_reps()+1);
	}
  
          
      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);      

    } // next permutation
  
  if (!par::silent)
    cout << "\n\n";
  
  
  

  /////////////////////////////////////////////////////
  //                                                 //
  // Calculate SET-based empirical p-values          //
  //                                                 //
  /////////////////////////////////////////////////////

  if (par::set_test && ! (par::set_r2 || par::set_score) ) 
    {
      printLOG("Calculating empirical SET-based p-values\n");
      pS->empiricalSetPValues();
    }
  
    

  /////////////////////////////////////////////////////
  //                                                 //
  //  Display basic permutation results              //
  //                                                 //
  /////////////////////////////////////////////////////

  ofstream ASC;
  string f;

  if (par::adaptive_perm) f = par::output_file_name + f0 + ".perm";
  else f = par::output_file_name + f0 + ".mperm";

  ASC.open(f.c_str(),ios::out);
  ASC.precision(4);
  printLOG("Writing permutation association results to [ " + f + " ] \n");

  if ( par::test_hap_GLM )
    ASC << setw(10) << "TEST" << " ";
  else if ( ! ( par::set_score || par::set_r2 ) ) 
    {
      ASC << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp)<< "SNP" << " "; 
    }
  else
    {
      ASC << setw(12) << "SET" << " "
	  << setw(6) << "NSNP" << " "
	  << setw(6) << "NSIG" << " "
	  << setw(6) << "ISIG" << " ";
    }
  
  ASC << setw(12) << "EMP1" << " ";
  
  if ( !par::set_r2 ) 
    {
      if (par::adaptive_perm)
	ASC << setw(12)<< "NP" << " ";
      else if ( par::mperm_rank )
	ASC	<< setw(12)<< "EMP3" << " "
		<< setw(12)<< "RANK" << " ";
      else
	ASC	<< setw(12)<< "EMP2" << " ";
    }
  else 
    ASC << "SNPS";
  
  ASC << "\n";
  
  vector<double> pv(0);
  
  for (int l=0; l<ntests; l++)
    {	

      // Skip?, if filtering p-values
      if ( par::pfilter && perm.pvalue(l) > par::pfvalue ) 
	continue;
      
      if ( par::test_hap_GLM )
	{
	  ASC << setw(10) << ("T"+int2str(l)) << " ";
	}
      else if ( ! (par::set_score || par::set_r2 ) )
	{
	  ASC << setw(4)  << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " ";
	}
      else
	{
	  ASC << setw(12) << setname[l] << " "
	      << setw(6) << pS->snpset[l].size() << " "
	      << setw(6) << pS->numSig[l] << " " 
	      << setw(6) << pS->selectedSNPs[l].size() << " ";
	  
	}

      // All tests are 1 df
      double p = chiprobP(original[l],1);
      
      // ... except 2df genotypic test
      if ( par::model_perm_gen )
	p = chiprobP(original[l],2);
      
      if (par::multtest)
	pv.push_back(p);
      
      ASC << setw(12) << perm.pvalue(l) << " ";

      if ( ! par::set_r2 )
	{
	  if (par::adaptive_perm) 
	    ASC << setw(12) << perm.reps_done(l) << " ";
	  else if ( par::mperm_rank )
	    ASC << setw(12) << perm.max_pvalue(l) << " "
		<< setw(12) << perm.rank(l) << " ";
	  else
	    ASC << setw(12) << perm.max_pvalue(l) << " ";
	}
      else
	{
	  if ( pS->selectedSNPs[l].size() == 0 )
	    ASC << "NA";
	  else
	    for (int j=0; j<pS->selectedSNPs[l].size(); j++)
	      {
		ASC << locus[ snpset[l][pS->selectedSNPs[l][j]] ]->name;
		if ( j < pS->selectedSNPs[l].size() - 1 ) 
		  ASC << "|";
	      }
	  
	}

      ASC << "\n";
      
    }
  
  ASC.close();
  
  


  ////////////////////////////////////////////////////////
  //                                                    //
  // Display SET-based results (sum statistics)         //
  //                                                    //
  ////////////////////////////////////////////////////////


  if (par::set_test && ! (par::set_r2 || par::set_score) )
    {
      
      f = par::output_file_name + f0 + ".set";
      ASC.open(f.c_str(),ios::out);
      printLOG("Writing set-based association results to [ " + f + " ] \n");

      ASC.clear();
      
      ASC << setw(12) << "SET"  << " " 
	  << setw(6) << "NSNP" << " "
	  << setw(6) << "S"  << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " " 
	  << setw(12) << "T" << " " 
	  << setw(12) << "P_0" << " " 
	  << setw(12) << "P_1" << " " 
	  << setw(12) << "P_2" << " " 
	  << "\n";
      
      for (int i=0;i<pS->pv_set.size();i++)
	{
	  if ( pS->pv_set[i].size()>0 ) ASC << "\n";
	  for (int j=0;j<pS->pv_set[i].size();j++)
	    {
	      ASC << setw(12) << setname[i]  << " " 
		  << setw(6) << pS->snpset[i].size() << " "
		  << setw(6) 
		  << string("S"+int2str(j+1+pS->s_min[i])) << " " 
		  << setw(par::pp_maxsnp)
		  << pS->setsort[i][j]  << " " 
		  << setw(12)
		  << pS->stat_set[i][j][0] << " " 
		  << setw(12)
		  << pS->pv_set[i][j][0] << " " 
		  << setw(12)
		  << pS->pv_maxG_set[i][j]/(par::replicates+1)  << " " 
		  << setw(12)
		  << pS->pv_maxE_set[i][j]/(par::replicates+1)  << "\n"; 
	      
	    }
	}
	            
      ASC.close();
    }  


  // Call destructor to close PDUMP file
  perm.closeDUMP();

}


/////////////////////////////////////////////////
// Basic C/C association test

vector<double> Plink::testAssoc(int & aff,    
				int & unf,    
				vector<int> & a1,
				vector<int> & a2, 
				vector<int> & a0,
				vector<double> & odds,
				vector<double> & exp_afffreq1,
				vector<double> & exp_afffreq2,
				vector<double> & exp_unffreq1,
				vector<double> & exp_unffreq2, 
				Perm & perm,
				bool display)
{	
  
  ofstream ASC;

  if (display)
    {
      ////////////////////
      // Display results
      
      string f = par::output_file_name + ".assoc";
      if ( par::fisher_test ) 
	f += ".fisher";
      ASC.open(f.c_str(),ios::out);
      ASC.precision(4);
      printLOG("Writing main association results to [ " + f + " ] \n");
      
      ASC << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp)<< "SNP" << " " 
	  << setw(10) << "BP" << " "
	  << setw(4) << "A1" << " ";
      if ( par::assoc_counts )
	ASC << setw(8) << "C_A" << " " 
	    << setw(8) << "C_U" << " "; 
      else
	ASC << setw(8) << "F_A" << " " 
	    << setw(8) << "F_U" << " "; 

      ASC << setw(4) << "A2" << " ";
      if ( ! par::fisher_test ) 
	ASC << setw(12)<< "CHISQ" << " ";
      ASC << setw(12)<< "P" << " " 
	  << setw(12)<< "OR" << " " ;
      if (par::display_ci)
	ASC << setw(12)<< "SE" << " " 
	    << setw(12) << string("L"+dbl2str(par::ci_level*100)) << " "
	    << setw(12) << string("U"+dbl2str(par::ci_level*100)) << " ";

      ASC << "\n";
   
    }

  vector<double> original(nl_all);
  
  
  ///////////////////////////////
  // Iterate over SNPs
  
  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {	
      
      // In adaptive mode, possibly skip this test
      if (par::adaptive_perm && (!perm.snp_test[l])) 
	{
	  s++;
	  l++;
	  continue;
	}

      int A1 = 0, A2 = 0, A0 = 0;
      int U1 = 0, U2 = 0, U0 = 0;
      
      bool X=false, haploid=false;
      
      if (par::chr_sex[locus[l]->chr]) X=true;
      else if (par::chr_haploid[locus[l]->chr]) haploid=true;
      

      /////////////////////////////
      // Iterate over individuals
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator gperson = sample.begin();
 
      while ( gperson != sample.end() )
	{
	  
	  // Phenotype for this person (i.e. might be permuted)
	  Individual * pperson = (*gperson)->pperson;

	  // Is this individual missing?
	  if ( pperson->missing )
	    {
	      // Next person
	      gperson++;
	      i1++;
	      i2++;	      	      
	      continue;
	    }
	  
	  // SNP alleles	  
	  bool s1 = *i1;
	  bool s2 = *i2;
	  	       
	  // Type of marker
	  if ( haploid || ( X && (*gperson)->sex ) ) 
	    {
	      
	      /////////////////////////////////////
	      // Haploid marker (or male X)
	      
	      if (pperson->aff)  // if affected 
		{	      
		  if (!s1)
		    {
		      if (!s2)  
			A1++;
		      
		    }
		  else
		    {
		      if ( s2 ) 
			A2++;
		      else A0++; 
		    }
		}
	      else // unaffected if not missing
		{	      
		  if (!s1)
		    {
		      if (!s2)   
			U1++;
		    }
		  else
		    {
		      if ( s2)   
			U2++;
		      else U0++; 
		    }
		}
	    }
	  else
	    {
	      /////////////////////////////////////
	      // Autosomal marker
	      
	      if (pperson->aff) // if affected 
		{	      
		  if (!s1)
		    {
		      if (!s2)   
			A1+=2;
		      else
			{ A1++; A2++; }  
		    }
		  else
		    {
		      if ( s2 )  
			A2+=2;
		      else A0+=2;  
		    }
		}
	      else  // unaffected if not missing
		{	      
		  if (!s1)
		    {
		      if (!s2)  
			U1+=2;
		      else
			{ U1++; U2++; } 
		    }
		  else
		    {
		      if ( s2 )  
			U2+=2;
		      else U0+=2;  
		    }
		}
	      
	    }

	  // Next person
	  gperson++;
	  i1++;
	  i2++;
	}
      
      
      // Calculate standard association statistic
      
      // Total number of alleles
      double tot = A1+A2+U1+U2;
      
      double pvalue;

      // Total number of non-missing affecteds/unaffecteds
      if ( par::fisher_test ) 
	{
	  table_t t;
	  sizeTable(t,2,2);
	  t[0][0] = A1;
	  t[0][1] = U1;
	  t[1][0] = A2;
	  t[1][1] = U2;
	  pvalue = fisher(t);
	  original[l] = 1 - pvalue;
	}
      else if (!par::assoc_test_alt_perm)
	{
	  double Taff = A1+A2;
	  double Tunf = U1+U2;
	  
	  double Ta1 = A1+U1;
	  double Ta2 = A2+U2;
	  
	  double Texp_afffreq1 = (Taff*Ta1)/tot;
	  double Texp_unffreq1 = (Tunf*Ta1)/tot;
	  double Texp_afffreq2 = (Taff*Ta2)/tot;
	  double Texp_unffreq2 = (Tunf*Ta2)/tot;

 	  original[l] = ((A1 - Texp_afffreq1) * (A1 - Texp_afffreq1)) / Texp_afffreq1
 	    + ( (U1 - Texp_unffreq1) * ( U1 - Texp_unffreq1 ) ) / Texp_unffreq1 
 	    + ( (A2 - Texp_afffreq2) * ( A2 - Texp_afffreq2 ) ) / Texp_afffreq2 
 	    + ( (U2 - Texp_unffreq2) * ( U2 - Texp_unffreq2 ) ) / Texp_unffreq2 ;
	  
	}
      else // if (par::assoc_test_alt_perm)
	{

	  // Total number of non-missing affecteds/unaffecteds
	  aff = A1+A2;
	  unf = U1+U2;

	  a1[l] = A1+U1;
	  a2[l] = A2+U2;
	  a0[l] = A0+U0;

	  exp_afffreq1[l] = (aff*a1[l])/tot;
	  exp_unffreq1[l] = (unf*a1[l])/tot;
	  exp_afffreq2[l] = (aff*a2[l])/tot;
	  exp_unffreq2[l] = (unf*a2[l])/tot;

	  // Include missing alleles for final marginal values
	  aff += A0;
	  unf += U0;
	
	  original[l] = ((A1 - exp_afffreq1[l]) * (A1 - exp_afffreq1[l])) / exp_afffreq1[l]
	    + ( (U1 - exp_unffreq1[l]) * ( U1 - exp_unffreq1[l] ) ) / exp_unffreq1[l] 
	    + ( (A2 - exp_afffreq2[l]) * ( A2 - exp_afffreq2[l] ) ) / exp_afffreq2[l] 
	    + ( (U2 - exp_unffreq2[l]) * ( U2 - exp_unffreq2[l] ) ) / exp_unffreq2[l] ;

	}

      
      ////////////////////////////////////////////////////
      // Do we need to calculate odds ratio and p-value? 
      // (Either for display (of original data) or for 
      // other reasons (profile set-based test)
      
      if ( display || par::set_score )
	{
	  // Note: in set-score mode, use an adjusted form 
	  // of the odds-ratio, to ensure it is always valud
	  
	  if ( par::set_score )
	    odds[l] = (double)( (A1+0.5)*(U2+0.5) ) / (double)( (U1+0.5)*(A2+0.5) ) ;
	  else
	    {	      
	      // with v. large sample N, better to use: ad/bc = a/b * d/c

	      //odds[l] = (double)( A1*U2 ) / (double)( U1*A2 ) ;
	      
	      odds[l] = ( (double)A1 / (double)A2 ) * 
		( (double)U2 / (double)U1 ) ;
	    }
	  
	  if ( ! par::fisher_test ) 
	    pvalue = chiprobP(original[l],1);
	}
      
      if ( par::set_score )
	{
	  if ( pvalue <= par::set_score_p && pvalue >= 0 )
	    pS->profileTestSNPInformation( l, log(odds[l]) * -log10( pvalue ) );
	}

      if (display)
	{
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && ( pvalue > par::pfvalue || pvalue < 0 ) ) 
	    goto skip_p1;
	  
	  
	  // Now display results for this SNP
	  
	  ASC << setw(4)  << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 
	      << setw(10) << locus[l]->bp << " "
	      << setw(4)  << locus[l]->allele1 << " ";

	  if ( A1+A2 != 0 )
	    {
	      if ( par::assoc_counts )
		ASC << setw(8)  << A1 << " ";
	      else
		ASC << setw(8)  << (double)A1/(double)(A1+A2) << " ";
	    }
	  else
	    ASC << setw(8)  << "NA" << " ";
	  
	  if ( U1+U2 != 0 )
	    {
	      if ( par::assoc_counts )
		ASC << setw(8)  << U1 << " ";
	      else
		ASC << setw(8)  << (double)U1/(double)(U1+U2) << " ";
	    }
	  else
	    ASC << setw(8)  << "NA" << " ";
	  
	  ASC << setw(4)  << locus[l]->allele2 << " " ;
	  
	  if ( par::fisher_test ) 
	    {
	      if ( pvalue > -1 )
		ASC << setw(12) << pvalue << " ";
	      else
		ASC << setw(12) << "NA" << " ";
	    }
	  else
	    {
	      if ( pvalue > -1 )
		ASC << setw(12) << original[l]  << " " 
		    << setw(12) << pvalue << " ";
	      else
		ASC << setw(12) << "NA"  << " " 
		    << setw(12) << "NA" << " ";
	    }

	  double zero=0;
	  if (odds[l] != odds[l] || odds[l] == 1/zero || odds[l] == -1/zero )
	    {
	      ASC << setw(12) << "NA" << " ";
	      if (par::display_ci)
		ASC << setw(12) << "NA" << " "
		    << setw(12) << "NA" << " "
		    << setw(12) << "NA" << " ";	      
	    }
	  else
	    {
	      ASC << setw(12) << odds[l]  << " " ;
	      if (par::display_ci)
		{
		  double lOR = log(odds[l]);
		  double SE = sqrt(1/(double)A1 + 1/(double)A2 + 1/(double)U1 + 1/(double)U2);
		  double OR_lower = exp( lOR - par::ci_zt * SE );
		  double OR_upper = exp( lOR + par::ci_zt * SE );
		  
		  ASC << setw(12) << SE << " "
		      << setw(12) << OR_lower << " "
		      << setw(12) << OR_upper << " ";	      
		}
	    }
	  ASC << "\n";
	}
      
    skip_p1:

      // Next SNP
      s++;
      l++;

    }
  
  if (display) ASC.close();

  return original;
}



/////////////////////////////////////////////
// Full model association tests

vector<double> Plink::fullModelAssoc(bool print_results,
				     Perm & perm)
{
  
  
  if (print_results)
    printLOG("Full-model association tests, minimum genotype count: --cell " + 
	     int2str(par::min_geno_cell) + "\n");
  
  vector<double> results(nl_all);
  
  ofstream ASC;
  
  if (print_results)
    {
      string f = par::output_file_name + ".model";
      ASC.open(f.c_str(),ios::out);
      printLOG("Writing full model association results to [ " + f + " ] \n");

      ASC << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " " 
	  << setw(4) << "A1" << " "
	  << setw(4) << "A2" << " "
	  << setw(8) << "TEST" << " "
	  << setw(14) << "AFF" << " " 
	  << setw(14) << "UNAFF" << " ";
      if ( ! par::fisher_test )
	ASC << setw(12) << "CHISQ" << " " 
	    << setw(4) << "DF" << " ";
      ASC << setw(12) << "P" << "\n";
 
      ASC.precision(4);
  }
  

  ///////////////////////////////
  // Iterate over SNPs
  
  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {	
      
      // In adaptive mode, possibly skip this test
      if (par::adaptive_perm && (!perm.snp_test[l])) 
	{
	  s++;
	  l++;
	  continue;
	}
      
      int A11=0, A12=0, A22=0;
      int U11=0, U12=0, U22=0;
      

      //////////////////////// 
      // Autosomal or haploid?
      
      bool X=false, haploid=false;
      if (par::chr_sex[locus[l]->chr]) X=true;
      else if (par::chr_haploid[locus[l]->chr]) haploid=true;
      
      // Skip haploid markers
      if (haploid)
	{
	  s++;
	  l++;
	  continue;
	}
      

      /////////////////////////////
      // Iterate over individuals
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator gperson = sample.begin();
 
      while ( gperson != sample.end() )
	{
	  
	  // Phenotype for this person (i.e. might be permuted)
	  Individual * pperson = (*gperson)->pperson;
	  
	  // SNP alleles	  
	  bool s1 = *i1;
	  bool s2 = *i2;
	  
	  if ( ! pperson->missing  )      
	    {
	      
	      // Only consider diploid chromosomes
	      if ( ! ( X && (*gperson)->sex ) ) 
		{		  
		  
		  if ( pperson->aff )     // cases
		    {
		      if ( ! s1 )
			{
			  if ( ! s2 )     // Homozyg 00
			    A11++;
			  else            // Hetero  01
			    A12++;                     
			}
		      else if ( s2 )      // Homozyg 11
			A22++;
		    }
		  else
		    {
		      if ( !s1 )
			{
			  if ( !s2 )      // Homozyg 00
			    U11++;
			  else            // Hetero  01
			    U12++;                     
			}
		      else if ( s2 )      // Homozyg 11
			U22++;
		    }	      
		  
		  
		}
	    }

	      
	  // Next individual
	  gperson++;
	  i1++;
	  i2++;
	  
	}
      
      
      ///////////////////////////////////
      // Calculate association statistics
      
      double obs_A = A11 + A12 + A22;
      double obs_U = U11 + U12 + U22;
      double obs_T = obs_A + obs_U;
      
      double obs_1 = 2*(A11+U11) + A12 + U12;
      double obs_2 = 2*(A22+U22) + A12 + U12;
      
      double obs_11 = A11+U11;
      double obs_12 = A12+U12;
      double obs_22 = A22+U22;

      bool invalid = false;
      if (A11 < par::min_geno_cell || 
	  A12 < par::min_geno_cell || 
	  A22 < par::min_geno_cell) invalid = true;
      else if (U11 < par::min_geno_cell || 
	       U12 < par::min_geno_cell || 
	       U22 < par::min_geno_cell) invalid = true;

      
      if ( par::trend_only ) 
	invalid = true;
      
      ///////////////////////
      // Cochram-Armitage Trend test

      double CA = ( ( obs_U / obs_T * A12 ) - ( obs_A / obs_T * U12 ) )
	+ 2*( ( obs_U / obs_T * A22 ) - ( obs_A / obs_T * U22 ) ) ;
      
      double varCA = obs_A * obs_U 
	* ( ( obs_T * ( obs_12 + 4*obs_22 ) 
	      - ( obs_12+2*obs_22 ) * ( obs_12+2*obs_22 )  ) 
	    / (obs_T * obs_T * obs_T  )) ;
      
      double CA_chisq = (CA*CA) / varCA;
      double CA_p = chiprobP(CA_chisq,1);

      double mult_p, mult_chisq;
      
      ///////////////////////
      // Multiplicative model
      
      double obs_A1 = 2*A11 + A12;
      double obs_A2 = 2*A22 + A12;
      double obs_U1 = 2*U11 + U12;
      double obs_U2 = 2*U22 + U12;
      
      if ( par::fisher_test ) 
	{
	  table_t t;
	  sizeTable(t,2,2);
	  t[0][0] = (int)obs_A1;
	  t[1][0] = (int)obs_A2;
	  t[0][1] = (int)obs_U1;
	  t[1][1] = (int)obs_U2;
	  mult_p = fisher(t);
	      
	}
      else
	{
	  
	  double exp_A1 = (obs_A * obs_1 ) / obs_T;    // note 2's cancelled for obs_A and obs_T 
	  double exp_A2 = (obs_A * obs_2 ) / obs_T;    // which are counts of individuals, not 
	  double exp_U1 = (obs_U * obs_1 ) / obs_T;    // alleles
	  double exp_U2 = (obs_U * obs_2 ) / obs_T;
	  
	  mult_chisq =  ( ( obs_A1 - exp_A1 ) * ( obs_A1 - exp_A1 ) ) / exp_A1
	    + ( ( obs_A2 - exp_A2 ) * ( obs_A2 - exp_A2 ) ) / exp_A2
	    + ( ( obs_U1 - exp_U1 ) * ( obs_U1 - exp_U1 ) ) / exp_U1
	    + ( ( obs_U2 - exp_U2 ) * ( obs_U2 - exp_U2 ) ) / exp_U2;
	  
	  	  
	  ///////////////////////
	  // Multiplicative model
	  
	  mult_p = chiprobP(mult_chisq,1);	
	  
	}
      
      
      double gen_p, dom_p, rec_p;
      gen_p = dom_p = rec_p = -9;
      double dom_chisq, rec_chisq, gen_chisq;


      if (!invalid)
	{
	  
	  //////////////////////////////////////////////////////////////
	  // Standard chi-square test, or Fisher's exact
	  
	  if ( par::fisher_test )
	    {

	      ////////////
	      // General

	      table_t t;
	      sizeTable(t,3,2);
	      t[0][0] = A11;
	      t[1][0] = A12;
	      t[2][0] = A22;
	      t[0][1] = U11;
	      t[1][1] = U12;
	      t[2][1] = U22;
	      gen_p = fisher(t);


	      ////////////
	      // Dominant

	      sizeTable(t,2,2);
	      t[0][0] = A11+A12;
	      t[1][0] = A22;
	      t[0][1] = U11+U12;
	      t[1][1] = U22;
	      dom_p = fisher(t);

	      /////////////
	      // Recessive

	      sizeTable(t,2,2);
	      t[0][0] = A11;
	      t[1][0] = A12+A22;
	      t[0][1] = U11;
	      t[1][1] = U12+U22;
	      rec_p = fisher(t);


	    }
	  else
	    {

	      ///////////////////////
	      // General model
	  
	      double exp_A11 = (obs_A * obs_11 ) / obs_T;
	      double exp_A12 = (obs_A * obs_12 ) / obs_T;
	      double exp_A22 = (obs_A * obs_22 ) / obs_T;
	      double exp_U11 = (obs_U * obs_11 ) / obs_T;
	      double exp_U12 = (obs_U * obs_12 ) / obs_T;
	      double exp_U22 = (obs_U * obs_22 ) / obs_T;
	      
	      gen_chisq =  ( ( A11 - exp_A11 ) * ( A11 - exp_A11 ) ) / exp_A11
		+ ( ( A12 - exp_A12 ) * ( A12 - exp_A12 ) ) / exp_A12
		+ ( ( A22 - exp_A22 ) * ( A22 - exp_A22 ) ) / exp_A22
		+ ( ( U11 - exp_U11 ) * ( U11 - exp_U11 ) ) / exp_U11
		+ ( ( U12 - exp_U12 ) * ( U12 - exp_U12 ) ) / exp_U12
		+ ( ( U22 - exp_U22 ) * ( U22 - exp_U22 ) ) / exp_U22;
	      
	      
	      ///////////////////////
	      // Dominant (minor allele) (1) model 
	      
	      dom_chisq =  ( ( (A11+A12) - (exp_A11+exp_A12) ) * ( (A11+A12) - (exp_A11+exp_A12) ) ) / (exp_A11+exp_A12) 
		+ ( ( A22 - exp_A22 ) * ( A22 - exp_A22 ) ) / exp_A22
		+ ( ( (U11+U12) - (exp_U11+exp_U12) ) * ( (U11+U12) - (exp_U11+exp_U12) ) ) / (exp_U11+exp_U12) 
		+ ( ( U22 - exp_U22 ) * ( U22 - exp_U22 ) ) / exp_U22;
	      
	      
	      //////////////////////////////////////
	      // Recessive (minor allele) (1) model 
	      
	      rec_chisq =  ( ( (A22+A12) - (exp_A22+exp_A12) ) * ( (A22+A12) - (exp_A22+exp_A12) ) ) / (exp_A22+exp_A12) 
		+ ( ( A11 - exp_A11 ) * ( A11 - exp_A11 ) ) / exp_A11
		+ ( ( (U22+U12) - (exp_U22+exp_U12) ) * ( (U22+U12) - (exp_U22+exp_U12) ) ) / (exp_U22+exp_U12) 
		+ ( ( U11 - exp_U11 ) * ( U11 - exp_U11 ) ) / exp_U11;
	      
	      
	      //////////////////////////////////
	      // p-values and model comparisons 
	      
	      gen_p = chiprobP(gen_chisq,2);
	      dom_p = chiprobP(dom_chisq,1);
	      rec_p = chiprobP(rec_chisq,1);      
	      
	    }
	  
	}

      
      ////////////////////////////////////////////
      // Save best p-value for permutation test

      
      //////////////////////////
      // Save the desired result

      int best = 0 ;
      if (par::model_perm_best)
	{
	  
	  double best_p = mult_p;
	  
	  if (!invalid) 
	    {
	      // Skip general model (i.e. just compare ALLELIC, DOM, REC

	      //if (gen_p < best_p && gen_p >= 0 ) { best = 2; best_p = gen_p; }   // general
	      
	      if (dom_p < best_p && dom_p >= 0 ) { best_p = dom_p; }   // dom
	      if (rec_p < best_p && rec_p >= 0 ) { best_p = rec_p; }   // rec
	    }
	  
	  results[l] = 1-best_p;
	}
      else if ( par::model_perm_gen )
	results[l] = gen_p >= 0 ? 1-gen_p : -9 ;
      else if ( par::model_perm_dom )
	results[l] = dom_p >= 0 ? 1-dom_p : -9 ;
      else if ( par::model_perm_rec )
	results[l] = rec_p >= 0 ? 1-rec_p : -9;
      else if ( par::model_perm_trend )
	results[l] = CA_p >= 0 ? CA_chisq : -9;

      if (print_results)
	{
	  
	  /////////////
	  // Genotypic 
	  
	  if ( ! par::trend_only ) 
	    {
	      
	      ASC << setw(4) << locus[l]->chr << " " 
		  << setw(par::pp_maxsnp) << locus[l]->name << " " 
		  << setw(4) << locus[l]->allele1 << " "
		  << setw(4) << locus[l]->allele2 << " "
		  << setw(8) << "GENO" << " "
		  << setw(14) << int2str(A11)+"/"+int2str(A12)+"/"+int2str(A22)  << " " 
		  << setw(14) << int2str(U11)+"/"+int2str(U12)+"/"+int2str(U22) << " " ;
	      if (gen_p < -1) 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << "NA"  << " " 
			<< setw(4)  << "NA" << " ";
		  
		  ASC << setw(12) << "NA" << "\n" ;
		}
	      else 
		{
		  if ( ! par::fisher_test ) 
		    ASC << setw(12) << gen_chisq << " " 
			<< setw(4)  << "2" << " ";
		  
		  ASC << setw(12) << gen_p << "\n";
		}
	    }

	  
	  /////////////////
	  // CA trend test

	  ASC << setw(4) << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 
	      << setw(4) << locus[l]->allele1 << " "
	      << setw(4) << locus[l]->allele2 << " "
	      << setw(8) << "TREND" << " "
	      << setw(14) << int2str(A11*2+A12)+"/"+int2str(A12+A22*2) << " " 
	      << setw(14) << int2str(U11*2+U12)+"/"+int2str(U12+U22*2) << " "; 
	  if (CA_p < -1) 
	    {
	      if ( ! par::fisher_test )
		ASC << setw(12) << "NA" << " " 
		    << setw(4) << "NA" << " ";
	      ASC << setw(12) << "NA" << "\n" ;
	    }
	  else 
	    {
	      if ( ! par::fisher_test )
		ASC << setw(12) << CA_chisq << " " 
		    << setw(4) << "1" << " ";
	      ASC << setw(12) << CA_p << "\n" ;
	    }
	  

	  if ( ! par::trend_only )
	    {

	      /////////////
	      // Allelic 
	      
	      ASC << setw(4) << locus[l]->chr << " " 
		  << setw(par::pp_maxsnp) << locus[l]->name << " " 
		  << setw(4) << locus[l]->allele1 << " "
		  << setw(4) << locus[l]->allele2 << " "
		  << setw(8) << "ALLELIC" << " "
		  << setw(14) << int2str(A11*2+A12)+"/"+int2str(A12+A22*2) << " " 
		  << setw(14) << int2str(U11*2+U12)+"/"+int2str(U12+U22*2) << " ";
	      if (mult_p < -1) 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << "NA" << " " 
			<< setw(4) << "NA" << " ";
		  ASC << setw(12) << "NA" << "\n" ;
		}
	      else 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << mult_chisq << " " 
			<< setw(4) << "1" << " ";
		  ASC << setw(12) << mult_p << "\n" ;
		}
	      
	      /////////////
	      // Dominant 
	      
	      ASC << setw(4) << locus[l]->chr << " " 
		  << setw(par::pp_maxsnp) << locus[l]->name << " " 
		  << setw(4) << locus[l]->allele1 << " "
		  << setw(4) << locus[l]->allele2 << " "
		  << setw(8) << "DOM" << " "
		  << setw(14) << int2str(A11+A12)+"/"+int2str(A22) << " " 
		  << setw(14) << int2str(U11+U12)+"/"+int2str(U22) << " ";
	      if (dom_p < -1) 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << "NA" << " " 
			<< setw(4) << "NA" << " ";
		  ASC << setw(12) << "NA" << "\n" ;
		}
	      else 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << dom_chisq << " " 
			<< setw(4) << "1" << " ";
		  ASC << setw(12) << dom_p << "\n" ;
		}
	      
	      /////////////
	      // Recessive
	      
	      ASC << setw(4) << locus[l]->chr << " " 
		  << setw(par::pp_maxsnp) << locus[l]->name << " " 	    
		  << setw(4) << locus[l]->allele1 << " "
		  << setw(4) << locus[l]->allele2 << " "
		  << setw(8) << "REC" << " "
		  << setw(14) << int2str(A11)+"/"+int2str(A12+A22) << " " 
		  << setw(14) << int2str(U11)+"/"+int2str(U12+U22) << " ";
	      if (rec_p < -1) 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << "NA" << " " 
			<< setw(4) << "NA" << " ";
		  ASC << setw(12) << "NA" << "\n" ;
		}
	      else 
		{
		  if ( ! par::fisher_test )
		    ASC << setw(12) << rec_chisq << " " 
			<< setw(4) << "1" << " ";
		  ASC << setw(12) << rec_p << "\n" ;
		}
	    }
	}
      
      // Next SNP
      s++;
      l++;
	       
    }

  if (print_results) 
    {
      ASC.close();
    }

  return results;

}



/////////////////////////////////////////////
// Simple quantitative trait association test
// note: does not explicitly treat X/haploid
// markers differently: see --linear

vector<double> Plink::testQAssoc(bool print_results , 
				 Perm & perm )
{	
  
  vector<double> results(nl_all);
  
  if ( print_results && par::multtest )
    tcnt.resize(nl_all);
  
  ofstream ASC, QT_MEANS;

  if (print_results)    
    {
      string f = par::output_file_name + ".qassoc";
      printLOG("Writing QT association results to [ " + f + " ] \n");
      ASC.open(f.c_str(),ios::out);
      ASC << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " " 
	  << setw(10) << "BP" << " "
	  << setw(8) << "NMISS" << " " 
	  << setw(10) << "BETA" << " "
	  << setw(10) << "SE" << " "
	  << setw(10) << "R2" << " "
	  << setw(8) << "T" << " "
	  << setw(12) << "P" << " "
	  << "\n";
      ASC.precision(4);
      
      if ( par::qt_means )
	{
	  string f = par::output_file_name + ".qassoc.means";
	  printLOG("Writing QT genotypic means to [ " + f + " ] \n");
	  QT_MEANS.open(f.c_str(),ios::out);
	  QT_MEANS.precision(4);
	  QT_MEANS << setw(4) << "CHR" << " "
		   << setw(par::pp_maxsnp) << "SNP" << " "
		   << setw(6) << "VALUE" << " " 
		   << setw(8) << "G11" << " "
		   << setw(8) << "G12" << " "
		   << setw(8) << "G22" << "\n"; 
	}
      
    }
  

  ////////////////////////////
  // Iterate over each locus
  
  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {	
      
      // Skip possibly
      if (par::adaptive_perm && !perm.snp_test[l])
	{
	  l++;
	  s++;
	  continue;
	}

      double g_mean=0, g_var=0;
      double qt_mean=0, qt_var=0;
      double qt_g_covar=0;
      
      // number of individuals in analysis
      int nanal = 0;  

      /////////////////////////////
      // Iterate over individuals
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator gperson = sample.begin();
 
      while ( gperson != sample.end() )
	{
	  
	  // Phenotype for this person (i.e. might be permuted)

	  Individual * pperson = (*gperson)->pperson;
	  
	  // SNP alleles
	  
	  bool s1 = *i1;
	  bool s2 = *i2;

	  if (!pperson->missing)
	    {
	      if ( ! ( s1 && (!s2) ) )     //   10 = missing 
		{
		  qt_mean += pperson->phenotype;
		  
		  if (!s1)
		    {
		      if (!s2)      //   00 = hom(11)
			g_mean+=2;
		      else          //   01 = het(12)
			g_mean++;
		    }
		  nanal++;
		}		
	      
	    }
	  
	  
	  // Advance to the next person (for phenotype information
	  // and the two SNP alleles also)
	  
	  gperson++;
	  i1++;
	  i2++;

	}
      
      qt_mean /= (double)nanal;
      g_mean /= (double)nanal;      
      
      
      //////////////////////////////////
      // Iterate over individuals again
      
      i1 = (*s)->one.begin();
      i2 = (*s)->two.begin();
      gperson = sample.begin();
      
      while ( gperson != sample.end() )
	{
	  
	  // Phenotype for this person (i.e. might be permuted)
	  Individual * pperson = (*gperson)->pperson;
	  
	  // SNP alleles
	  bool s1 = *i1;
	  bool s2 = *i2;
	  
	  if (!pperson->missing)
	    {
	      if ( ! ( (s1) && (!s2) ) ) 
		{

		  qt_var += (pperson->phenotype-qt_mean) 
		    * ( pperson->phenotype-qt_mean ) ;
		  
		  double g = 0;
		  
		  if (!s1)
		    {
		      if (!s2)
			g=2;
		      else
			g=1;
		    }
		  
		  g_var += (g-g_mean) * ( g-g_mean ) ;
		  qt_g_covar += ( pperson->phenotype - qt_mean ) 
		    * ( g - g_mean ) ;
		  
		}
	    }

	  // Advance to the next person
	  gperson++;
	  i1++;
	  i2++;
	  
	}

      // Summary statistics

      qt_var /= (double)nanal - 1;
      g_var /= (double)nanal - 1;
      qt_g_covar /= (double)nanal - 1;
      
      // Test statistics
      double beta = qt_g_covar / g_var;
      double vbeta = ( qt_var/g_var - (qt_g_covar*qt_g_covar)/(g_var*g_var) ) / (nanal-2);
      double t = beta / sqrt(vbeta);
      double t_p;

      // Display results?
      
      if (print_results)
	{

// 	  double wald = (beta*beta) / ( vbeta ) ;
// 	  double wald_p = chiprobP(wald,1);

	  t_p = pT(t,nanal-2);
	  double r2 =   (qt_g_covar * qt_g_covar ) / ( qt_var * g_var ) ;
	  
//        double lrt = -nanal * log(1-r2);

	  // Skip?, if filtering p-values
	  
	  if ( par::pfilter && ( t_p > par::pfvalue || t_p < 0 ) ) 
	    goto skip_p2;
	  
	  ASC << setw(4) << locus[l]->chr  << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 
	      << setw(10) << locus[l]->bp << " "
	      << setw(8) << nanal << " ";
	  if ( ! realnum(beta) )
	    {
	      ASC << setw(10) << "NA" << " " 
		  << setw(10) << "NA" << " " 
		  << setw(10) << "NA" << " " ;
	    }
	  else
	    {
	      ASC << setw(10) << beta << " " 
		  << setw(10) << sqrt(vbeta) << " " 
		  << setw(10) << r2 << " " ;
	    }

 	  if (t_p >= 0) 
 	    ASC << setw(8) << t  << " " << setw(12) << t_p << " " ;
 	  else
 	    ASC << setw(8) << "NA"  << " " << setw(12) << "NA"  << " " ;
	  
	  ASC << "\n";	

	  if ( par::qt_means ) 
	    displayQTMeans(QT_MEANS, l);
	  
	}
      
      
    skip_p2:


      // Store chi-sq
      results[l] = t*t;

      // Store original p-value (for --adjust)
      if ( print_results && par::multtest )
	{
	  tcnt[l] = nanal-2;
	}

      // Next SNP
      s++;
      l++;

    }
  
  if (print_results)
    {
      ASC.close();
      if ( par::qt_means ) 
	QT_MEANS.close();
    }
  return results;
}




/////////////////////////////////////////////
// For a given SNP, calculate genotypic mean, 
// frequency and variance, and display 

void Plink::displayQTMeans(ofstream & QT_MEANS, int l)
{
  
  vector<CSNP*>::iterator s = SNP.begin()+ l;
  
  double g11=0, g12=0, g22=0;
  double x11=0, x12=0, x22=0;  
  double xx11=0, xx12=0, xx22=0;
  
  /////////////////////////////
    // Iterate over individuals
    
    vector<bool>::iterator i1 = (*s)->one.begin();
    vector<bool>::iterator i2 = (*s)->two.begin();
    vector<Individual*>::iterator person = sample.begin();
    
    while ( person != sample.end() )
      {
	
	// SNP alleles
	
	bool s1 = *i1;
	bool s2 = *i2;
	
	if ( ! (*person)->missing )
	  {
	    
	    if ( ! s1 ) 
	      {
		if ( ! s2 ) 
		  {
		    g11++;
		    x11 += (*person)->phenotype;
		    xx11 += (*person)->phenotype * (*person)->phenotype;
		  }
		else
		  {
		    g12++;
		    x12 += (*person)->phenotype;
		    xx12 += (*person)->phenotype * (*person)->phenotype;
		  }
	      }
	    else
	      {
		if ( s2 ) 
		  {
		    g22++;
		    x22 += (*person)->phenotype;
		    xx22 += (*person)->phenotype * (*person)->phenotype;		    
		  }
	      }
	  }
	
      	person++;
	i1++;
	i2++;
	
      }

    double nanal = g11 + g12 + g22;

    x11 /= g11;
    x12 /= g12;
    x22 /= g22;

    xx11 /= g11;
    xx12 /= g12;
    xx22 /= g22;

    double sd11 = g11>1 ? sqrt(xx11 - x11 * x11) * sqrt(g11/(g11-1)) : 0;
    double sd12 = g12>1 ? sqrt(xx12 - x12 * x12) * sqrt(g12/(g12-1)) : 0;
    double sd22 = g22>1 ? sqrt(xx22 - x22 * x22) * sqrt(g22/(g22-1)) : 0;
    
    string a1 = locus[l]->allele1;
    string a2 = locus[l]->allele2;
    if ( a1 == "" ) a1 = "*";    
    if ( a2 == "" ) a2 = "*";
        
    QT_MEANS << setw(4) << locus[l]->chr << " "
	     << setw(par::pp_maxsnp) << locus[l]->name << " "
	     << setw(6) << "GENO" << " "
	     << setw(8) << a1+"/"+a1 << " "
	     << setw(8) << a1+"/"+a2 << " "
	     << setw(8) << a2+"/"+a2 << "\n";
    
    QT_MEANS << setw(4) << locus[l]->chr << " "
	     << setw(par::pp_maxsnp) << locus[l]->name << " "
	     << setw(6) << "COUNTS" << " "
	     << setw(8) << g11 << " "
	     << setw(8) << g12 << " "
	     << setw(8) << g22 << "\n";

    g11 /= nanal;
    g12 /= nanal;
    g22 /= nanal;

    QT_MEANS << setw(4) << locus[l]->chr << " "
	     << setw(par::pp_maxsnp) << locus[l]->name << " "
	     << setw(6) << "FREQ" << " "
	     << setw(8) << g11 << " "
	     << setw(8) << g12 << " "
	     << setw(8) << g22 << "\n";

    QT_MEANS << setw(4) << locus[l]->chr << " "
	     << setw(par::pp_maxsnp) << locus[l]->name << " "
	     << setw(6) << "MEAN" << " ";
    if ( g11>0 ) 
      QT_MEANS << setw(8) << x11 << " ";
    else
      QT_MEANS << setw(8) << "NA" << " ";
    if ( g12>0 ) 
      QT_MEANS << setw(8) << x12 << " ";
    else
      QT_MEANS << setw(8) << "NA" << " ";
    if ( g22>0)
      QT_MEANS << setw(8) << x22 << "\n";
    else
      QT_MEANS << setw(8) << "NA" << "\n";

    QT_MEANS << setw(4) << locus[l]->chr << " "
	     << setw(par::pp_maxsnp) << locus[l]->name << " "
	     << setw(6) << "SD" << " ";
    if ( g11>0 ) 
      QT_MEANS << setw(8) << sd11 << " ";
    else
      QT_MEANS << setw(8) << "NA" << " ";
    if ( g12>0 ) 
      QT_MEANS << setw(8) << sd12 << " ";
    else
      QT_MEANS << setw(8) << "NA" << " ";
    if ( g22>0)
      QT_MEANS << setw(8) << sd22 << "\n";
    else
      QT_MEANS << setw(8) << "NA" << "\n";

}


//////////////////////////////////////////////////////////////
// Test difference in missingness between cases and controls

vector<double> Plink::testMiss(Perm & perm, bool display)
{	
  
  // Requires SNP-major mode
  if (!par::SNP_major) Ind2SNP();

  ofstream MIS;
  
  if (display)
    {
      
      string f = par::output_file_name + ".missing";
      MIS.open(f.c_str(),ios::out);
      MIS.precision(4);
      printLOG("Writing case/control missingness test to [ " + f + " ] \n");
      
      MIS << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp)<< "SNP" << " " 
	  << setw(12) << "F_MISS_A" << " "
	  << setw(12) << "F_MISS_U" << " " 
	  << setw(12)<< "P" << " " 
	  << "\n";
    }

  vector<double> missing(0);
  
  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {	
      
      // Skip possibly
      if (par::adaptive_perm && !perm.snp_test[l])
	{
	  l++;
	  s++;
	  continue;
	}
      
      int affmiss = 0;
      int affgeno = 0;
      int unfmiss = 0;
      int unfgeno = 0;
      

      ///////////////////////////////
      // Iterate over each individual
            
      vector<Individual*>::iterator gperson = sample.begin();
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      
      while ( gperson != sample.end() )
	{

	  Individual * pperson = (*gperson)->pperson;

	  // If we haven't excluded this individual
	  // on the basis of a cluster solution

	  if ( ! pperson->missing )
	    {
	      
	      // "1" allele count
	      
	      if ( pperson->aff) // affected 
		{	      
		  if ( *i1 &&  ! *i2 )
		    affmiss++;
		  else
		    affgeno++;
		}
	      else
		{
		  if ( *i1 && ! *i2 )
		    unfmiss++;
		  else
		    unfgeno++;
		}
	      
	    }
	  
	  // Next individual

	  gperson++;
	  i1++;
	  i2++;
	}
    
      
      // Calculate association statistic
      
      table_t t;
      sizeTable(t,2,2);
      t[0][0] = affmiss;
      t[0][1] = affgeno;
      t[1][0] = unfmiss;
      t[1][1] = unfgeno;
            
      double pvalue = fisher(t);
      
      // Record 1-p as empirical statistic

      if ( pvalue > -1 ) 
	missing.push_back(1-pvalue);
      else
	missing.push_back(1);

      if (display)
	{
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && pvalue > par::pfvalue ) 
	    {
	      // Next SNP
	      s++;
	      l++;
      	      continue;
	    }
	  
	  // Total number of affecteds/unaffecteds

	  int aff = affmiss + affgeno;
	  int unf = unfmiss + unfgeno;
	  
	  MIS << setw(4)  << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 
	      << setw(12) << affmiss / (double)aff << " "
	      << setw(12) << unfmiss / (double)unf << " ";

	  if ( pvalue > -1 ) 
	    MIS << setw(12) << pvalue << "\n";
	  else
	    MIS << setw(12) << "NA" << "\n";
	}


      // Next SNP

      s++;
      l++;

    }

  if (display)
    {
      MIS.close();
    }

  return missing;
}



void Plink::calcLDStatistics()
{
  
  ///////////////////////////////////////////
  // Calculate simple correlation (r or r^2) 
  // based on 3x3 genotype counts  
  
  // If --matrix option also specified, output as 
  // one matrix, otherwise use  
      
  ofstream LD;
  string f = par::output_file_name + ".ld";
  LD.open(f.c_str(),ios::out);
  
  printLOG("Writing LD statistics to [ " + f + " ] \n");
  if (!par::matrix)
    {
      LD << setw(6) << "CHR_A" << " "
	 << setw(12) << "BP_A" << " "
	 << setw(par::pp_maxsnp) << "SNP_A"  << " "
	 << setw(6) << "CHR_B"  << " " 
	 << setw(12) << "BP_B" << " "
	 << setw(par::pp_maxsnp) << "SNP_B" << " " ;
      if (par::disp_r1)
	LD << setw(12) << "R" << " " ;
      else
	LD << setw(12) << "R2" << " " ;
      LD << "\n";
    }


  set<int> ldAnchorSet;
  if ( par::ld_anchor_list )
    {
      checkFileExists( par::ld_SNP1_file );
      
      ifstream IN( par::ld_SNP1_file.c_str() , ios::in );
      
      map<string,int> mlocus;
      for (int l=0; l<nl_all; l++)
	mlocus.insert( make_pair( locus[l]->name , l ) );

      while ( !IN.eof() )
	{
	  string snp;
	  IN >> snp;
	  if (snp == "")
	    continue;
	  map<string,int>::iterator i = mlocus.find( snp );
	  if ( i != mlocus.end() )
	    ldAnchorSet.insert( i->second );
	}
      IN.close();
    }

  int ld_anchor_number = -1;
  if ( par::ld_anchor && ! par::ld_anchor_list )
    {      
      par::matrix = false;
      ld_anchor_number = getMarkerNumber(*this, par::ld_SNP1);
      if (ld_anchor_number == -1)
	error("--ld-snp {marker} not found");
    }
  
  
  int end = nl_all;
  // Do we need to go up to the last SNP?
  if ( (!par::matrix) && (!par::ld_anchor) ) end--;
  

  ///////////////////////////
  // First locus  
  
  for (int l1=0; l1<end; l1++)  
    {
      
      // Full matrix, or just lower-diagonal?
      
      int start=0;
      if ( (!par::matrix) && (!par::ld_anchor) ) start = l1+1;

      // Using an LD anchor?
      
      if ( par::ld_anchor )
	{
	  if ( par::ld_anchor_list )
	    {
	      if ( ldAnchorSet.find ( l1 ) == ldAnchorSet.end() )
		continue;
	    }
	  else
	    if ( l1 != ld_anchor_number )
	      continue;
	}
      
      // Second locus
	  
      for (int l2=start; l2<nl_all; l2++)      
	{
	  
	  // Outside of window?
	  if (par::disp_r_window)
	    {
      
	      if ( l2 - l1 >= par::disp_r_window_snp )
		continue;

	      if ( l1 - l2 >= par::disp_r_window_snp )
		continue;

	      if ( locus[l2]->chr != locus[l1]->chr ) 
		continue;

	      if ( locus[l2]->bp - locus[l1]->bp 
		   > par::disp_r_window_kb )
		continue;

	      if ( locus[l1]->bp - locus[l2]->bp 
		   > par::disp_r_window_kb )
		continue;
	    }


	  
	  ////////////////////////
	  // Calculate correlation
	  
	  double r = correlation2SNP(l1,l2,par::disp_r2,false);
	  
	  if (par::matrix)
	    {	  	  
	      LD << r << " ";
	    }
	  else
	    {
	      // Using a r^2 threshold? 
	      if ( par::disp_r1 ||
		   r >= par::disp_r_window_r2 )
		{
		  LD << setw(6) << locus[l1]->chr  << " " 
		     << setw(12) << locus[l1]->bp << " "
		     << setw(par::pp_maxsnp) << locus[l1]->name  << " " 
		     << setw(6) << locus[l2]->chr  << " " 
		     << setw(12) << locus[l2]->bp << " "
		     << setw(par::pp_maxsnp) << locus[l2]->name << " " 
		     << setw(12) << r << " " << "\n";		  
		}
	    }
	}
    
      if (par::matrix) LD << "\n";      
    }
  
  LD.close();
  

}



double Plink::correlation2SNP(int l1, int l2, bool squared, bool covariance, bool useflag)
{

  // Calculate simple correlation based on 0,1,2 allele counts
  // i.e. the 3x3 genotypic table, rather than the 2x2 haplotypic
  // table. 

  // r = cov(1,2) / sqrt( var(1).var(2) )	  
	  
  double X = 0;
  double X2 = 0;
  double Y = 0;
  double Y2 = 0;
  double XY = 0;
  double count = 0;
  
  bool haploid_snp1 = par::chr_haploid[locus[l1]->chr];
  bool X_snp1 = par::chr_sex[locus[l1]->chr];

  bool haploid_snp2 = par::chr_haploid[locus[l2]->chr];
  bool X_snp2 = par::chr_sex[locus[l2]->chr];
	  
  // Iterate over every individual
  // but only consider founders
  
  // Have a specicial case loop for when both SNPs
  // are autosomal (usual case)
  
  if ( haploid_snp1 || X_snp1 || 
       haploid_snp2 || X_snp2 )
    {
      
      // Allow for 1 or both SNPs to be non-autosomal
      
      for (int i=0; i<n; i++)
	{
	  
	  Individual * person = sample[i];
	  
	  // Only consider founders
	  if ( ! person->founder ) continue; 
	  
	  // Are we using a flag? 
	  if ( useflag && ! person->flag ) continue;



	  bool a1 = par::SNP_major ? SNP[l1]->one[i] : person->one[l1];
	  bool a2 = par::SNP_major ? SNP[l1]->two[i] : person->two[l1];
	  if ( a1 && (!a2) ) continue;
	  
	  bool b1 = par::SNP_major ? SNP[l2]->one[i] : person->one[l2];
	  bool b2 = par::SNP_major ? SNP[l2]->two[i] : person->two[l2];
	  if ( b1 && (!b2) ) continue;
	  
	  
	  // Only consider if non-missing at both loci
	  
	  if ( a1 && (!a2) ) continue;
	  if ( b1 && (!b2) ) continue;
	  
	  // Score individuals
	  
	  count++;
	  
	  int sx = 0, sy = 0;
	  
	  // Haploid/diploid marker 1 ?
	  
	  if ( haploid_snp1 || ( X_snp1 && person->sex ) )
	    {		
	      // Hemizygous "0"
	      if ( ! a1 ) sx=1;
	    }	
	  else
	    {
	      
	      // Score 2,1,0 for 00,01,11 genotypes
	      if ( ! a1 )
		{
		  if ( ! a2 ) 
		    sx=2;
		  else
		    sx=1;
		}
	    }
	  
	  // Haploid/diploid marker 2 ?
	  
	  if ( haploid_snp2 || ( X_snp2 && person->sex ) )
	    {		
	      // Hemizygous "0"
	      if ( ! b1 ) sy=1;
	    }	
	  else
	    {
	      
	      // Score 2,1,0 for 00,01,11 genotypes
	      if ( ! b1 )
		{
		  if ( ! b2 )
		    sy=2;
		  else      
		    sy=1;
		}
	    }
	  
	  X += sx;
	  Y += sy;
	  XY += sx*sy;
	  
	  // Sum squares
	  sx *= sx;
	  sy *= sy;
	  X2 += sx;
	  Y2 += sy;
	  
	  // consider next person;
	}	      
    }
  else
    {
      
      // Autosomal only version
      
      for (int i=0; i<n; i++)
	{
	  
	  Individual * person = sample[i];

	  // Only consider founders
	  if ( ! person->founder ) continue; 

	  // Are we using a flag? 
	  if ( useflag && ! person->flag ) continue;

	  
	  // Only consider if non-missing at both loci
	  
	  bool a1 = par::SNP_major ? SNP[l1]->one[i] : person->one[l1];
	  bool a2 = par::SNP_major ? SNP[l1]->two[i] : person->two[l1];
	  if ( a1 && (!a2) ) continue;
	  
	  bool b1 = par::SNP_major ? SNP[l2]->one[i] : person->one[l2];
	  bool b2 = par::SNP_major ? SNP[l2]->two[i] : person->two[l2];
	  if ( b1 && (!b2) ) continue;
	  

	  // Score individuals
	  
	  count++;
	  
	  int sx = 0, sy = 0;
	  
	  if ( ! a1 )
	    {
	      if ( ! a2 )
		sx=2;
	      else 
		sx=1;
	    }		 	 
	  
	  if ( ! b1 )
	    {
	      if ( ! b2 ) 
		sy=2;
	      else  
		sy=1;
	    }
	  
	  X += sx;
	  Y += sy;
	  XY += sx*sy;
	  sx *= sx;
	  sy *= sy;
	  X2 += sx;
	  Y2 += sy;
	  
	  // consider next person;
	}
    }
  
  
  // count refers to number of individuals
	  
  X /= count;
  X2 /= count;
  Y /= count;
  Y2 /= count;
  XY /= count;
  
  double var1 = X2 - X*X;
  double var2 = Y2 - Y*Y;
  double cov12 = XY - X*Y;
  double r;


  // Return either:
  //  covariance, correlation or correlation squared?
  
  if ( covariance ) 
    r = cov12;
  else
    {
      if ( squared ) 
	r = (cov12*cov12) / (var1*var2);
      else
	{
	  r = cov12 / (sqrt(var1)*sqrt(var2));
	  
	  // get sign of label: check minor allele assignment matches
	  
	  if ( ( (locus[l1]->allele1 > locus[l1]->allele2)	  
		 && (locus[l2]->allele2 < locus[l2]->allele2) )
	       ||
	       ( (locus[l1]->allele1 < locus[l1]->allele2)	  
		 && (locus[l2]->allele2 > locus[l2]->allele2) ) )
	    r *= -1;
	}
    }
  
  return r;
}
