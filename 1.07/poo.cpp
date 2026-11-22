

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
#include "stats.h"

////////////////////////////
// Parent-of-origin analysis

void Plink::perm_testTDT_POO(Perm & perm)
{
  
  //////////////////////////////////
  // Individual-major mode analysis
  
  if (par::SNP_major) SNP2Ind();
  

  ///////////////////////////////////////////
  // Calculate original results for true data

  vector<bool> dummy(family.size(),false);
  perm.setTests(nl_all);
  perm.setPermClusters(*this);

  vector<double> original = testTDT_POO(true, false, perm, dummy, dummy);

  
  ////////////////////////////  
  // Display corrected p-values?
  
  if (par::multtest)
    {
      vector<double> obp(0);
      for (int l=0; l<nl_all;l++)
	obp.push_back(original[l]);
      multcomp(obp,".tdt");
    }

  ////////////////////////////////
  // If no permutations requested, 
  // we can finish here

  if (!par::permute) return;

  
  //////////////////////
  // Make sets?
 
  if (par::set_test) pS->cumulativeSetSum_WITHLABELS(*this,original);



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
      vector<bool> fB(family.size(),false);
      
      for (int f=0; f<family.size(); f++)
	{
	  if (CRandom::rand() < 0.5) fA[f] = true;
	  if (CRandom::rand() < 0.5) fB[f] = true;
	}
      
      vector<double> pr = testTDT_POO(false, true, perm, fA, fB);
      
      //////////////////////
      // Make sets?

      if (par::set_test)
        pS->cumulativeSetSum_WITHOUTLABELS(pr,perm.current_reps()+1);

      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);

    } // next permutation

  cout << "\n\n";


  ///////////////////////////////////////////
  // Calculate SET-based empirical p-values
  
  if (par::set_test) 
    {
      printLOG("Calculating empirical SET-based p-values\n");
      pS->empiricalSetPValues();
    }


  ////////////////////
  // Display results
  
  ofstream TDT;
  string f;
  if (par::adaptive_perm) f = par::output_file_name + ".tdt.poo.perm";
  else f = par::output_file_name + ".tdt.poo.mperm";

  TDT.open(f.c_str(),ios::out);
  printLOG("Writing TDT parent-of-origin permutation results to [ " + f + " ] \n"); 
  TDT.precision(4);

  TDT << setw(4) << "CHR" << " "
      << setw(par::pp_maxsnp) << "SNP" << " ";
  
  TDT << setw(12) << "CHISQ_TDT" << " ";
  
  TDT << setw(12) << "EMP1" << " ";
  if (par::adaptive_perm)
    TDT << setw(12) << "NP" << " " << "\n";  
  else
    TDT << setw(12) << "EMP2" << " " << "\n";  

  for (int l=0; l<nl_all; l++)
    {	
      
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
  
  TDT.close();

  ////////////////////////////  
  // Display SET-based results
 
  if (par::set_test)
    {
      
      f = par::output_file_name + ".tdt.poo.set";
      TDT.open(f.c_str(),ios::out);
      printLOG("Writing set-based TDT parent-of-origin results to [ " +f+ " ] \n");
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
		  << string("S<"+int2str(j+1)) << " "
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



vector<double> Plink::testTDT_POO(bool print_results, 
				  bool permute,
				  Perm & perm,
				  vector<bool> & flipA,
				  vector<bool> & flipB)
{
  
  
  
  ///////////////////////////
  // Vector to store results 

  vector<double> res(nl_all);
  double zt;

  ofstream TDT;

  if (print_results)
    {
      string f = par::output_file_name + ".tdt.poo";
      TDT.open(f.c_str(),ios::out);
      printLOG("Writing TDT parent-of-origin results (asymptotic) to [ " + f + " ] \n");
      TDT << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " "
	  << setw(6) << "A1:A2" << " "
	  << setw(12) << "T:U_PAT" << " "
	  << setw(12) << "CHISQ_PAT" << " "
	  << setw(12) << "P_PAT" << " "
      	  << setw(12) << "T:U_MAT" << " "
	  << setw(12) << "CHISQ_MAT" << " "
	  << setw(12) << "P_MAT" << " "
	  << setw(12) << "Z_POO" << " "
	  << setw(12) << "P_POO" << " ";
      

//       if (par::display_ci)
// 	TDT << setw(12) << string("L"+int2str(int(par::ci_level*100))) << " "
// 	    << setw(12) << string("U"+int2str(int(par::ci_level*100))) << " ";
//       if (par::display_ci)
// 	zt = ltqnorm( 1 - (1 - par::ci_level) / 2  ) ; 

      
      TDT << "\n";
      
    }
  

  ///////////////////////////////////
  // Perform analysis for each locus
  
  for (int l=0; l<nl_all; l++)
    {

      // Adaptive permutation, skip this SNP?
      if (par::adaptive_perm && (!perm.snp_test[l])) 
	continue;
                        
      // Transmission counts
      
      double p1 = 0;
      double p2 = 0;
      
      double m1 = 0;
      double m2 = 0;

      // Count over families

      for (int f=0; f<family.size(); f++)
	{

	  if ( ! family[f]->TDT ) continue;
	
	  int trP = 0;  // transmitted allele from het father
	  int unP = 0;  // untransmitted allele from het father
	  
	  int trM = 0;  // transmitted allele from het mother
	  int unM = 0;  // untransmitted allele from het mother
	  	 
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

	      bool hhh = false; // flag for het X het => het 

	      // Kid is 00

	      if ( (!kid1) && (!kid2) ) 
		{
		  
		  // Paternal transmission?
		  if ( (!pat1) && pat2 ) 
		    {
		      trP=1; 
		      unP=2;
		    }
		  
		  // Maternal transmission?
		  if ( (!mat1) && mat2 )
		    {  
		      trM=1; 
		      unM=2; 
		    }

		}
	      else if ( (!kid1) && kid2 )  // Kid is 01
		{
		  
		  // Everybody heterozygous? 
		  if ( pat1 != pat2 && mat1 != mat2 )
		    hhh = true;
		  else
		    {
		      // het father
		      if ( pat1 != pat2 )
			{ 
			  // what did mother transmit?
			  if ( !mat1 )
			    {
			      trP=2; 
			      unP=1;
			    }
			  else
			    {
			      trP=1; 
			      unP=2;
			    }
			}
		      else
			{
			  // what did father transmit?
			  if ( !pat1 )
			    {
			      trM=2; 
			      unM=1;
			    }
			  else
			    {
			      trM=1; 
			      unM=2;
			    }
			}
		    }		
		}
	      else // kid is 1/1
		{
		  
		  // Paternal transmission?
		  if ( (!pat1) && pat2 ) 
		    {
		      trP=2; 
		      unP=1;
		    }
		  
		  // Maternal transmission?
		  if ( (!mat1) && mat2 )
		    {  
		      trM=2; 
		      unM=1; 
		    }		  
		  
		}
	      
	      
	      ///////////////
	      // Permutation?
	      
	      if (permute) { 
		
		// Determine whether to flip parental origin...
		if (par::perm_POO_poo)
		  {
		    if (flipA[f]) { int t=trP; trP=trM; trM=t;  } 
		    if (flipB[f]) { int t=unP; unP=unM; unM=t;  } 
		  }
		else // ... or allelic transmission
		  {
		    if (flipA[f]) { int t=trP; trP=unP; unP=t;  } 
		    if (flipB[f]) { int t=trM; trM=unM; unM=t;  } 
		  }

	      }

	      
	      // Increment transmission counts
	      
	      if (hhh)
		{
		  p1 += 0.5;
		  p2 += 0.5;
		  m1 += 0.5;
		  m2 += 0.5;
		}
	      else
		{
		  if (trP==1) p1++;
		  if (trM==1) m1++;
		  if (trP==2) p2++;
		  if (trM==2) m2++; 
		}

	    } // next offspring in family
	  
	}  // next nuclear family
      
      
      
      /////////////////////////////
      // Finished counting: now compute
      // the statistics
      
      double pat_chisq, mat_chisq, tot_chisq;
      pat_chisq = mat_chisq = tot_chisq = -1;
      
      // Basic TDT test
      
      if (p1+p2 > 0)
	pat_chisq = ((p1-p2)*(p1-p2))/(p1+p2);
      
      if (m1+m2 > 0)
	mat_chisq = ((m1-m2)*(m1-m2))/(m1+m2);
      
      double t1 = p1 + m1;
      double t2 = p2 + m2;

      if (t1+t2 > 0)
	tot_chisq = ((t1-t2)*(t1-t2))/(t1+t2);

      double pat_OR = p1 / p2;
      double pat_VOR = 1/p1 + 1/p2;
      double mat_OR = m1 / m2;
      double mat_VOR = 1/m1 + 1/m2;
      // Test of POO effect
      double z = ( log(pat_OR) - log(mat_OR) ) / sqrt( pat_VOR + mat_VOR );
      

      // Display asymptotic results
      
      if (print_results)
	{
	  TDT.precision(4);
	  
	  TDT << setw(4) << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 	
	      << setw(6) 
	      << string(locus[l]->allele1 + ":" + locus[l]->allele2) << " ";
	  
	  // Paternal transmissions
	  TDT << setw(12) << dbl2str(p1)+":"+dbl2str(p2) << " ";
	  	  
	  if (pat_chisq>=0)
	    TDT << setw(12) << pat_chisq << " "
		<< setw(12) << chiprobP(pat_chisq,1) << " ";
	  else
	    TDT << setw(12) << "NA" << " "
		<< setw(12) << "NA" << " "; 


	  // Maternal transmissions
	  TDT << setw(12) << dbl2str(m1)+":"+dbl2str(m2) << " ";
	  
	  if (mat_chisq>=0)
	    TDT << setw(12) << mat_chisq << " "
		<< setw(12) << chiprobP(mat_chisq,1) << " ";
	  else
	    TDT << setw(12) << "NA" << " "
		<< setw(12) << "NA" << " "; 


	  if ( realnum(z) )
	    TDT << setw(12) << z << " "
		<< setw(12) << normdist(-fabs(z)) * 2 << " ";
	  else
	    {
	      TDT << setw(12) << "NA" << " "
		  << setw(12) << "NA" << " ";
	    }

	  TDT << "\n";
	}


      ///////////////////////////////////////////
      // Choose which statistic for permutation
      
      if (par::perm_POO_poo) res[l] = realnum(z) ? z*z : -1 ;
      else if (par::perm_POO_pat) res[l] = pat_chisq;
      else if (par::perm_POO_mat) res[l] = mat_chisq;
      else if (par::perm_POO_best) res[l] = pat_chisq > mat_chisq ? pat_chisq : mat_chisq;
      
    } // next locus


  //////////////////////////////
  // Close output file, if open

  if (print_results)
    TDT.close();

  ///////////////////////////////////////////
  // Return chosen statistic for permutation

  return res;

}
