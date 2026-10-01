

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream> 
#include <map>
#include <vector>
#include <set>
#include <cmath>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "crandom.h"
#include "sets.h"
#include "perm.h"
#include "fisher.h"
#include "helper.h"
#include "stats.h"


void Plink::perm_testTDT(Perm & perm)
{

  ////////////////////////////////////////
  // This function is the entry point for
  // both TDT and DFAM tests (this is the 
  // wrapper around the test functions for 
  // permutation, set-based tests, etc).



  //////////////////////////////////
  // Individual-major mode analysis
  
  if (par::SNP_major) 
    SNP2Ind();
  
  if ( ! par::bt ) 
    error("This analysis requires a binary disease phenotype");


  if ( par::set_r2 )
    {
      printLOG("Performing LD-based set test, with parameters:\n");
      printLOG("     r-squared  (--set-r2)   = " + dbl2str( par::set_r2_val ) + "\n" );
      printLOG("     p-value    (--set-p)    = " + dbl2str( chiprobP(par::set_chisq_threshold,1) ) + "\n" );
      printLOG("     max # SNPs (--set-max)  = " + int2str( par::set_max ) + "\n" );
	
      pS->makeLDSets();
    }
  
  if ( par::sibTDT_test )
    {

      ////////////////////////////////////////////////////////////////
      // Parse family and cluster sets, to ensure we do not count
      // anybody twice -- i.e. remove anybody who is in a family from
      // CMH-like analysis

      // Make a temporary set of Individuals*
      
      set<Individual*> plist;
      for (int i=0; i<n; i++)
	plist.insert(sample[i]);
  
      for (int k=0; k<nk; k++)
	{

	  for (int i=0; i<klist[k]->person.size(); i++)      
	    {

	    // Exclude non-singletons from list of people
	    
	    // Individual might have been filtered out, in which case
	    // we should revise klist[] in any case;
	    
	    if ( plist.find( klist[k]->person[i] ) == plist.end() )
	      {
		klist[k]->person.erase(klist[k]->person.begin() + i);
		i--;
	      }
	    else
	      {
		Individual * person = klist[k]->person[i];
		Family * fam = klist[k]->person[i]->family;
				
		if ( fam->parents || fam->sibship  ) 
		  {
		    klist[k]->person[i]->sol = -1;
		    klist[k]->person.erase(klist[k]->person.begin() + i);
		    i--;
		  }   

	      }
	    }
	}
    }
  

  ///////////////////////////////////////////
  // Calculate original results for true data

  vector<bool> dummy(family.size(),false);


  ///////////////////////////////////////////
  // Create cluster for permutation
  //  a) Permute within family
  //  b) Use any existing --within cluster scheme for unrelateds.
  
  // The preGeneDrop() function will blank sol for all individuals in
  // families; if a cluster has been loaded in, certain clusters will
  // be set to zero possibly. This is fine -- all we need to do is now
  // go through and add new clusters for each family. Start adding
  // from nk onwards (we just keep any zero-sized clusters in the
  // analysis, they will not harm anything). We only need to put
  // siblings in clusters (i.e. parents do not come into this)
    
    if ( par::sibTDT_test )
      for (int f=0; f<family.size(); f++)
	{

	  Family * fam = family[f];
	  
	  if ( fam->singleton )
	    continue;
	  
	  if ( fam->kid.size() < 2 )
	    continue;
	  
	  klist.push_back( new Cluster );
	  for (int c=0; c < fam->kid.size(); c++)
	    {
	      fam->kid[c]->sol = nk;
	      klist[nk]->person.push_back( fam->kid[c] );
	    }
	  nk++;
	}
    
    /////////////////////////////////
    // Determine the number of tests
    
    int ntests = nl_all;
    
    if ( par::set_r2 || par::set_score ) 
      ntests = pS->snpset.size();
    

    /////////////////////////////////
    // Empirical p-values
    
    perm.setTests(ntests);
    perm.setPermClusters(*this);
  
    string testname = ".tdt";
    if (par::sibTDT_test)
      testname = ".dfam";
    
  vector<double> original;
  
  if (par::sibTDT_test)
    original = testSibTDT(true, false, perm, dummy, dummy);
  else
    original = testTDT(true, false, perm, dummy, dummy);
  
  
  ////////////////////////////  
  // Display corrected p-values?
  
  if (par::multtest)
    {
      vector<double> obp(0);
      for (int l=0; l<nl_all;l++)
	obp.push_back(original[l]);
      multcomp(obp,testname);
    }

  ////////////////////////////////
  // If no permutations requested, 
  // we can finish here

  if (!par::permute) return;

  
  //////////////////////
  // Make sets?
  
  vector<int> setsigsize;
    
  if (par::set_test) 
    {
      if ( par::set_r2 )
	{
	  original = pS->fitLDSetTest(original,true);

	  // ...and save # of significant SNPs	  
	  setsigsize.clear();
	  for (int i=0; i<pS->profileSNPs.size(); i++)
	    setsigsize.push_back( pS->s_min[i] );

	}
      pS->cumulativeSetSum_WITHLABELS(*this,original);
    }

  //////////////////////
  // Begin permutations
  
  bool finished = false;

  while(!finished)
    {
      
      ///////////////////////////////////
      // Set up permutation list for TDT
      
      // Permutations are constant across family and markers
      // flipA/B[permutation][family]
      
      vector<bool> fA(family.size(),false);
      vector<bool> fP(family.size(),false);
      
      for (int f=0; f<family.size(); f++)
	{
	  if (CRandom::rand() < 0.5) fA[f] = true;
	  if (CRandom::rand() < 0.5) fP[f] = true;
	}
      
      // And also the label-swapping permutation for sibships 
      // and rest of sample

      perm.permuteInCluster();
      
      vector<double> pr;

      if (par::sibTDT_test)
	pr = testSibTDT(false, true, perm, fA, fP);
      else
	pr = testTDT(false, true, perm, fA, fP);
      
      //////////////////////
      // Make sets?

      if (par::set_test)
	{
	  if ( par::set_r2 )
	    pr = pS->fitLDSetTest(pr,false);
	  else
	    pS->cumulativeSetSum_WITHOUTLABELS(pr,perm.current_reps()+1);
	}

      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);

    } // next permutation

  if (!par::silent)
    cout << "\n\n";


  ///////////////////////////////////////////
  // Calculate SET-based empirical p-values
  
  if (par::set_test && ! (par::set_r2 || par::set_score) ) 
    {
      printLOG("Calculating empirical SET-based p-values\n");
      pS->empiricalSetPValues();
    }


  ////////////////////
  // Display results
  
  ofstream TDT;
  string f;

  if ( par::set_r2 ) 
    {
      
      if (par::adaptive_perm) f = par::output_file_name + testname + ".set.perm";
      else f = par::output_file_name + testname + ".set.mperm";
      
      TDT.open(f.c_str(),ios::out);
      TDT.precision(4);
      printLOG("Writing set-based results to [ " + f + " ] \n");
      
      TDT << setw(12) << "SET" << " "
	  << setw(6) << "NSNP" << " "
	  << setw(6) << "NSIG" << " "
	  << setw(6) << "ISIG" << " "
	  << setw(12)<< "STAT" << " " 
	  << setw(12) << "EMP1" << " "
	  << "SNPS" << "\n";
  
      vector<double> pv(0);
  
      for (int l=0; l<ntests; l++)
	{	
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && perm.pvalue(l) > par::pfvalue ) 
	    continue;
	  
	  TDT << setw(12) << setname[l] << " "
	      << setw(6) << pS->snpset[l].size() << " "
	      << setw(6) << pS->numSig[l] << " " 
	      << setw(6) << pS->selectedSNPs[l].size() << " ";
	  TDT << setw(12) << original[l]  << " " 
	      << setw(12) << perm.pvalue(l) << " ";

	  if ( pS->selectedSNPs[l].size() == 0 )
	    TDT << "NA";
	  else
	    for (int j=0; j<pS->selectedSNPs[l].size(); j++)
	      {
		TDT << locus[ snpset[l][pS->selectedSNPs[l][j]] ]->name;
		if ( j < pS->selectedSNPs[l].size() - 1 ) 
		  TDT << "|";
	      }
	  
	  TDT << "\n";
	}
    }
  else
    {
      // Standard empirical p-value reports

      string f;
      if (par::adaptive_perm) f = par::output_file_name + testname + ".perm";
      else f = par::output_file_name + testname + ".mperm";
      
      TDT.open(f.c_str(),ios::out);
      printLOG("Writing TDT permutation results to [ " + f + " ] \n"); 
      TDT.precision(4);
      
      TDT << setw(4) << "CHR" << " "
	  << setw(par::pp_maxsnp) << "SNP" << " ";
      
      if (par::perm_TDT_basic) TDT << setw(12) << "CHISQ_TDT" << " ";
      else if (par::perm_TDT_parent) TDT << setw(12) << "CHISQ_PAR" << " ";
      else TDT << setw(12) << "CHISQ_COM" << " ";
      
      TDT << setw(12) << "EMP1" << " ";
      if (par::adaptive_perm)
	TDT << setw(12) << "NP" << " " << "\n";  
      else
	TDT << setw(12) << "EMP2" << " " << "\n";  
      
      for (int l=0; l<nl_all; l++)
	{	
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && perm.pvalue(l) > par::pfvalue ) 
	    continue;
	  
	  TDT << setw(4) << locus[l]->chr << " "
	      << setw(par::pp_maxsnp) << locus[l]->name << " "; 
	  if (original[l] < -0.5)
	    TDT << setw(12) << "NA"  << " " 
		<< setw(12) << "NA"  << " " 
		<< setw(12) << "NA";
	  else
	    {
	      TDT << setw(12) << original[l] << " "
		  << setw(12) << perm.pvalue(l) << " ";
	      
	      if (par::adaptive_perm)
		TDT << setw(12) << perm.reps_done(l);
	      else
		TDT << setw(12) << perm.max_pvalue(l);
	    }
	  TDT << "\n";
	}
    }

  TDT.close();

  ////////////////////////////  
  // Display SET-based results
 
  if (par::set_test && ! par::set_r2 )
    {
      
      f = par::output_file_name + testname + ".set";
      TDT.open(f.c_str(),ios::out);
      printLOG("Writing set-based TDT results to [ " +f+ " ] \n");
      TDT.clear();

      // Header row
      TDT << setw(12) << "SET" << " "
	  << setw(6) << "S"  << " "
	  << setw(par::pp_maxsnp) << "SNP" << " "
	  << setw(12) << "T" << " "
	  << setw(12) << "P_0" << " "
	  << setw(12) << "P_1" << " "
	  << setw(12) << "P_2" << " "
	  << "\n";

      for (int i=0;i<pS->pv_set.size();i++)
	{
	  TDT << "\n";
	  for (int j=0;j<pS->pv_set[i].size();j++)
	    {
	      TDT << setw(12) << setname[i]  << " "
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
		  << pS->pv_maxE_set[i][j]/(par::replicates+1)  << " "
		  << "\n"; 
	    }
	}
      
      TDT.close();
    }


  
}



vector<double> Plink::testTDT(bool print_results, 
			      bool permute,
			      Perm & perm,
			      vector<bool> & flipA,
			      vector<bool> & flipP)
{
  
  
  // TDT and X chromosome: males are coded as homozygous i.e. father
  // should always be uninformative;
  // male child will always receive his X from father
  // females as usual
  
  ///////////////////////////
  // Vector to store results 

  vector<double> res(nl_all);
  double zt;

  ofstream TDT, MT;

  if (print_results)
    {
      string f = par::output_file_name + ".tdt";
      TDT.open(f.c_str(),ios::out);
      printLOG("Writing TDT results (asymptotic) to [ " + f + " ] \n");
      TDT << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " "
	  << setw(12) << "BP" << " "
	  << setw(3) << "A1" << " "
	  << setw(3) << "A2" << " "
	  << setw(6) << "T" << " "
	  << setw(6) << "U" << " "
	  << setw(12) << "OR" << " ";
      
      if (par::display_ci)
	TDT << setw(12) << string("L"+dbl2str(par::ci_level*100)) << " "
	    << setw(12) << string("U"+dbl2str(par::ci_level*100)) << " ";
      
      TDT << setw(12) << "CHISQ" << " "
	  << setw(12) << "P" << " ";
      
      if (par::discordant_parents)
	TDT << setw(12) << "A:U_PAR" << " "
	    << setw(12) << "CHISQ_PAR" << " "
	    << setw(12) << "P_PAR" << " "
	    << setw(12) << "CHISQ_COM" << " "
	    << setw(12) << "P_COM" << " ";
      
      TDT << "\n";
      
      if ( par::mating_tests )
	{
	  MT.open( (par::output_file_name + ".mt").c_str(), ios::out);
	  MT.precision(3);
	}

      if (par::display_ci)
	zt = ltqnorm( 1 - (1 - par::ci_level) / 2  ) ; 
  }


  ///////////////////////////////////
  // Perform analysis for each locus
  
  for (int l=0; l<nl_all; l++)
    {
                        
      // Adaptive permutation, skip this SNP?
      if (par::adaptive_perm && (!perm.snp_test[l])) 
	continue;

      // Transmission counts
      
      double t1 = 0;
      double t2 = 0;
      
      // Count over families

      for (int f=0; f<family.size(); f++)
	{

	  if ( ! family[f]->TDT ) continue;

	  int trA = 0;  // transmitted allele from first het parent
	  int unA = 0;  // untransmitted allele from first het parent
	  
	  int trB = 0;  // transmitted allele from second het parent
	  int unB = 0;  // untransmitted allele from second het parent
	  
	  Individual * pat = family[f]->pat;
	  Individual * mat = family[f]->mat;
	  vector<Individual *> kid = family[f]->kid;

	  bool pat1 = pat->one[l];
	  bool pat2 = pat->two[l];
 
	  bool mat1 = mat->one[l];
	  bool mat2 = mat->two[l];
	  
	  // We need two genotyped parents, with 
	  // at least one het

	  if ( pat1 == pat2 && 
	       mat1 == mat2 ) 
  	    continue;

	  if ( ( pat1 && !pat2 ) || 
	       ( mat1 && !mat2 ) ) 
	    continue;


	  // Consider all offspring in nuclear family
	  
	  for (int c=0; c<kid.size(); c++)
	    {

	      // Only consider affected children
	      if ( ! kid[c]->aff ) continue;
	      
	      bool kid1 = kid[c]->one[l];
	      bool kid2 = kid[c]->two[l];
	      
	      // Skip if offspring has missing genotype
	      if ( kid1 && !kid2 ) continue;
	      
	      // We've now established: no missing genotypes
	      // and at least one heterozygous parent

	      // Kid is 00

	      if ( (!kid1) && (!kid2) ) 
		{
		  if ( ( (!pat1) && pat2 ) && 
		       ( (!mat1) && mat2 ) )
		    { trA=1; unA=2; trB=1; unB=2; }
		  else 
		    { trA=1; unA=2; } 
		}
	      else if ( (!kid1) && kid2 )  // Kid is 01
		{
		  // het dad
		  if (pat1 != pat2 )
		    {
		      // het mum
		      if ( mat1 != mat2 )
			{ trA=1; trB=2; unA=2; unB=1; }
		      else if ( !mat1 ) 
			{ trA=2; unA=1; }
		      else { trA=1; unA=2; }
		    }
		  else if ( !pat1 ) 
		    {
		      trA=2; unA=1; 
		    }		    
		  else
		    {
		      trA=1; unA=2;
		    }
		}
	      else // kid is 1/1
		{
		  
		  if ( ( (!pat1) && pat2 ) && 
		       ( (!mat1) && mat2 ) )
		    { trA=2; unA=1; trB=2; unB=1; }
		  else 
		    { 
		      trA=2; unA=1;
		    }
		}
	      
	      // We have now populated trA (first transmission) 
	      // and possibly trB also 
	      
	      ////////////////////////////////////////
	      // Permutation? 50:50 flip (precomputed)
	      
	      if (permute) { 
		if (flipA[f]) 
		  { 
		    int t=trA; 
		    trA=unA; 
		    unA=t;  
		    
		    t=trB; 
		    trB=unB; 
		    unB=t;  
		  } 
	      }
	      
	      // Increment transmission counts
	      
	      if (trA==1) t1++;
	      if (trB==1) t1++;
	      if (trA==2) t2++;
	      if (trB==2) t2++; 
	      
	      if ( par::verbose)
		{
		  cout << "TDT\t" << locus[l]->name << " " 
		       << pat->fid << " : " 
		       << trA << " " 
		       << trB << "\n";
		}
	      
	    } // next offspring in family
	  
	}  // next nuclear family
      
      

      /////////////////////////////////////////////
      // Consider parental discordance information
      
      double p1 = 0; 
      double p2 = 0;
      double d1 = 0;
      double d2 = 0;

      if (par::discordant_parents)
	{
	  
	  // Count over families
	  
	  for (int f=0; f<family.size(); f++)
	    {

	      // Requires parental discordance...
	      if ( ! family[f]->discordant_parents ) 
		continue;
	      
	      Individual * pat = family[f]->pat;
	      Individual * mat = family[f]->mat;

	      bool pat1 = pat->one[l];
	      bool pat2 = pat->two[l];
	      
	      bool mat1 = mat->one[l];
	      bool mat2 = mat->two[l];
	      
	      
	      // ...and that both are genotyped
	      
	      if ( ( pat1 && !pat2 ) || 
		   ( mat1 && !mat2 ) ) 
		continue;
	      

	      ////////////////////////////////////////
	      // Permutation? 50:50 flip (precomputed)
	      
	      if (permute)
		{
		  if (flipP[f]) 
		    {
		      if (pat->aff) { pat->aff = false; mat->aff = true; }
		      else { pat->aff = true; mat->aff = false; }
		    }
		}
	      
	      // Get number of 'F' alleles that the affected parent has
	      // above the unaffected; this count is p1/d1
	      
	      //  excess T alleles in affected -> p1/d1
	      //  excess F alleles in unaffected -> p2/d2

	      if ( pat1 == mat1 &&
		   pat2 == mat2 ) 
		continue;  // d = 0;
	      else if ( pat->aff ) // affected pat 
		{
		  
		  if ( (!pat1) && (!pat2) ) // F/F
		    {
		      // mat will either be T/T or F/T
		      if ( mat1 ) d1++; // two extra T
		      else p1++; // one extra T
		    }
		  else if ( (!pat1 ) && pat2 ) // pat F/T
		    {
		      // mat either T/T or F/F
		      if ( mat1 ) 
			p1++; // one extra T
		      else
			p2++; // one less T
		    }
		  else // pat must be T/T
		    {
		      // mat will either be F/F or F/T
		      if ( ! mat2 ) d2++; // two less T
		      else p2++; // one less T
		    }
		} 
	      else // affected mat / score other direction
		{
		  if ( (!pat1) && (!pat2) ) // F/F
		    {
		      // mat will either be T/T or F/T
		      if ( mat1 ) d2++;
		      else p2++;
		    }
		  else if ( (!pat1 ) && pat2 ) // pat F/T
		    {
		      // mat either T/T or F/F
		      if ( mat1 ) 
			p2++;
		      else
			p1++;	      
		    }
		  else // pat must be T/T
		    {
		      // mat will either be F/F or F/T
		      if ( ! mat2 ) d1++;
		      else p1++;
		    }
		}
	    }
	  
	}
      

      ///////////////////////////////////
      // General family test

      if ( par::mating_tests ) 
	{
	  
	  table_t parenMT;
	  sizeTable(parenMT,3,3);

	  // Count over families
	  
	  for (int f=0; f<family.size(); f++)
	    {

	      Individual * pat = family[f]->pat;
	      Individual * mat = family[f]->mat;

	      if ( pat == NULL || mat == NULL ) 
		continue;

	      bool pat1 = pat->one[l];
	      bool pat2 = pat->two[l];
	      
	      bool mat1 = mat->one[l];
	      bool mat2 = mat->two[l];
	      
	      
	      // ...and that both are genotyped
	      
	      if ( ( pat1 && !pat2 ) || 
		   ( mat1 && !mat2 ) ) 
		continue;
	      

	      int i=0, j=0;
	      
	      if ( pat1 ) 
		++i;
	      if ( pat2 ) 
		++i;
	      
	      if ( mat1 ) 
		++j;
	      if ( mat2 ) 
		++j;
	      
	      ++parenMT[i][j];
	    }
	  
	  
	  double mean1 = 0, mean2 = 0; 
	  int total = 0;
	  for(int i=0; i<=2; i++)
	    for (int j=0; j<=2; j++)
	      {
		mean1 += parenMT[i][j] * i;
		mean2 += parenMT[i][j] * j;
		total += parenMT[i][j];
	      }
	  
	  mean1 /= (double)total;
	  mean2 /= (double)total;
	  
	  double var1 = 0, var2 = 0, covar = 0;
	  for(int i=0; i<=2; i++)
	    for (int j=0; j<=2; j++)
	      {
		var1 += ( i - mean1 )*(i-mean1)*parenMT[i][j];
		var2 += ( j - mean2 )*(j-mean2)*parenMT[i][j];
		covar += ( i - mean1 )*(j-mean2)*parenMT[i][j];
	      }
	  
	  var1 /= (double)total - 1.0;
	  var2 /= (double)total - 1.0;
	  covar /= (double)total - 1.0;
	  
	  double r = covar / sqrt( var1 * var2 );
	  	
	  //	  double t = fisher(parenMT);

	  double t = chiTable(parenMT);

	  MT << setw(4) << locus[l]->chr << " " 
	     << setw(par::pp_maxsnp) << locus[l]->name << " " 	
	     << setw(12) << locus[l]->bp << " "
	     << setw(8) << locus[l]->freq << " "
	     << setw(12) << t << " ";	  

	  t = symTable(parenMT);
	  MT << setw(12) << t << " ";

	  MT << setw(12) << r << " ";

	  for(int i=0; i<=2; i++)
	    for (int j=0; j<=2; j++)
		MT <<  setw(5) << parenMT[i][j] << " ";

	  MT << "\n";

	}



      
      /////////////////////////////
      // Finished counting: now compute
      // the statistics
      
      double tdt_chisq, par_chisq, com_chisq;
      tdt_chisq = par_chisq = com_chisq = -1;
      
      // Basic TDT test
      
      if (t1+t2 > 0)
	tdt_chisq = ((t1-t2)*(t1-t2))/(t1+t2);
      
      if (par::discordant_parents)
	{
	  // parenTDT 
	  
	  if ( p1+p2+d1+d2 > 0 ) 
	    par_chisq = (((p1+2*d1)-(p2+2*d2))*((p1+2*d1)-(p2+2*d2)))
	      /(p1+p2+4*(d1+d2));
	  
	  
	  // Combined test
	  
	  if ( t1+p1+4*d1+t2+p2+4*d2 > 0 )
	    com_chisq = ( ( (t1+p1+2*d1) - (t2+p2+2*d2) ) 
			  * ( (t1+p1+2*d1) - (t2+p2+2*d2) ) ) 
	      / ( t1+p1+4*d1+t2+p2+4*d2 ) ;
	}
      
      
      // Display asymptotic results
      
      if (print_results)
	{

	  double pvalue = chiprobP(tdt_chisq,1);
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && pvalue > par::pfvalue ) 
	    continue;

	  TDT.precision(4);
	  
	  TDT << setw(4) << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 	
	      << setw(12) << locus[l]->bp << " "
	      << setw(3) << locus[l]->allele1 << " "
	      << setw(3) << locus[l]->allele2 << " "
	      << setw(6) << t1 << " "
	      << setw(6) << t2 << " ";
     
	  // Odds ratio for T:U
	  double OR = t1 / t2;

	  if ( ! realnum(OR) ) 
	    {
	      TDT << setw(12) << "NA" << " ";
	      if (par::display_ci) 
		TDT << setw(12) << "NA" << " "
		    << setw(12) << "NA" << " ";
	    }
	  else
	    {
	      TDT << setw(12) << OR << " ";
	      
	      if (par::display_ci)
		{
		  double OR_lower = exp( log(OR) - zt * sqrt(1/t1+1/t2)) ;
		  double OR_upper = exp( log(OR) + zt * sqrt(1/t1+1/t2)) ;
		  
		  TDT << setw(12) << OR_lower << " "
		      << setw(12) << OR_upper << " ";
		}
	    }

	  if (tdt_chisq>=0)
	    TDT << setw(12) << tdt_chisq << " "
		<< setw(12) << chiprobP(tdt_chisq,1) << " ";
	  else
	    TDT << setw(12) << "NA" << " "
		<< setw(12) << "NA" << " "; 
	  
	  if (par::discordant_parents)
	    {
	      TDT << setw(12) 
		  <<  dbl2str(p1+2*d1)+":"+dbl2str(p2+2*d2) << " ";

	      if (par_chisq>=0)
		TDT << setw(12) << par_chisq  << " "
		    << setw(12) << chiprobP(par_chisq,1) << " ";
	      else
		TDT << setw(12) << "NA" << " "	
		    << setw(12) << "NA" << " "; 
	      
	      if (com_chisq>=0)
		TDT << setw(12) << com_chisq << " "
		    << setw(12) << chiprobP(com_chisq,1) << " ";
	      else
		TDT << setw(12) << "NA" << " "
		    << setw(12) << "NA" << " "; 
	    }
	  TDT << "\n";
	}


      ///////////////////////////////////////////
      // Choose which statistic for permutation

      if (par::perm_TDT_basic) res[l] = tdt_chisq;
      else if (par::perm_TDT_parent) res[l] = par_chisq;
      else res[l] = com_chisq;
    
    } // next locus


  //////////////////////////////
  // Close output file, if open

  if (print_results)
    {
      TDT.close();
      if ( par::mating_tests )
	MT.close();
    }

  ///////////////////////////////////////////
  // Return chosen statistic for permutation

  return res;

}

