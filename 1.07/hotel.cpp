

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
#include <map>
#include <vector>
#include <cmath>

#include "plink.h"
#include "sets.h"
#include "options.h"
#include "helper.h"
#include "stats.h"
#include "perm.h"

using namespace std;

// Helper function

void calcHotelSetMeanVariance(vector<CSNP*> &,
			      vector<double> &,
			      vector<double> &,
			      vector<vector<double> > &,
			      vector<Individual*>&,
			      int,int);


void Plink::perm_testHotel(Perm & perm)
{
  
  if (!par::SNP_major)
    Ind2SNP();

  // Do not allow monomorphic alleles
  if (par::min_af==0)
    error("Cannot specify --maf 0 when using --T2; set --maf > 0");

  // Do not allow completely missing SNPs
  if (par::MAX_GENO_MISSING==1)
    error("Cannot specify --geno 1 when using --T2; set --geno < 1");

  // Are we using sets? If so, construct these now
  if (!par::set_test)
    error("You need to specify sets (--set option) with --T2");

  // Do not allow quantitative traits
  if (!par::bt)
    error("Cannot specify --T2 with quantitative traits");
  
 
  // Prune SET (0-sized sets, MAF==0 SNPs, etc) 
  pS->pruneSets(*this);

  int ns = snpset.size();
  


  ///////////////////////////////////////////
  // Count how many cases, how many controls

  int caseN = 0;
  int controlN = 0;

  for (int i=0; i < n; i++)
    if (!sample[i]->missing)
      if (sample[i]->aff)
	caseN++;
      else
	controlN++;

  if ( caseN == 0 || 
       controlN == 0 ) 
    error("No cases / no controls for T(2) test");
  
  // Multi-collinearity SNP pruning
  setFlags(true);
  pS->pruneMC(*this,true,par::vif_threshold);
      
  // Empirical p-values (1 per set)
  perm.setTests(ns);


  ////////////////////////////////
  // Set up permutation structure 
  // (we need to perform this step
  //  whether or not we also 
  //  subsequently permute)

  perm.setPermClusters(*this);
  perm.originalOrder();
  

  vector<double> original = calcHotel(true, 
				      perm,
				      *pS,
				      caseN, controlN);


  ////////////////////////////
  // If no permutation, then 
  // leave now

  if (!par::permute) return;


  //////////////////////
  // Begin permutations

  bool finished = false;
  while(!finished)
    {
      
      // Store permuted results
      
      vector<double> pr(ns);
      
      if (par::perm_genedrop)
	perm.geneDrop();
      else
	perm.permuteInCluster();

      
      pr = calcHotel(false, 
		     perm,
		     *pS,
		     caseN, controlN);
      
      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);

    } // next permutation
  
  if (!par::silent)
    cout << "\n\n";
  

  ////////////////////
  // Display results
  
  ofstream ASC;
  string f;
  if (par::adaptive_perm) f = par::output_file_name + ".T2.perm";
  else f = par::output_file_name + ".T2.mperm";

  ASC.open(f.c_str(),ios::out);
  ASC.precision(4);
  printLOG("Writing permutation T2 test results to [ " + f + " ] \n");

  ASC << setw(12) << "SET" << " " 
      << setw(4)<< "SIZE" << " " 
      << setw(12) << "EMP1" << " ";
  if (par::adaptive_perm)
    ASC << setw(12)<< "NP" << " ";
  else
    ASC	<< setw(12)<< "EMP2" << " ";
  ASC << "\n";
  

  for (int s=0; s<snpset.size(); s++)
    {	
      ASC << setw(12)  << setname[s] << " " 
	  << setw(4)   << snpset[s].size() << " ";
      
      ASC << setw(12) << perm.pvalue(s) << " ";
      
      if (par::adaptive_perm) 
	ASC << setw(12) << perm.reps_done(s) << " ";
      else
	ASC << setw(12) << perm.max_pvalue(s) << " ";
      
      ASC << "\n";
      
    }
  
 ASC.close();


}

vector<double> Plink::calcHotel(bool disp, 
				Perm & perm, 
				Set & S, 
				int ncase, 
				int ncontrol)
{

  ofstream ASC;
  if (disp)
    {
      
      string f = par::output_file_name + ".T2";
      
      ASC.open(f.c_str(),ios::out);
      ASC.precision(4);
      printLOG("Writing T2 test results to [ " + f + " ] \n");
      
      ASC << setw(12) << "SET" << " " 
	  << setw(4)<< "SIZE" << " " 
	  << setw(12)<< "T2" << " " 
	  << setw(12) << "DF1" << " "
	  << setw(12)<< "DF2" << " "
	  << setw(12)<< "P_HOTEL" << "\n";
    }

  // Number of SETs
  int ns = pS->snpset.size();
  vector<double> T2(ns,0);
  


  // Consider each SET
  for (int s=0; s<ns; s++)
    {
      
      // Adaptive permutation: skip possibly
      if (par::adaptive_perm && !perm.snp_test[s]) continue;
      
      // Consider each SNP, coded 1,0 and -1
      // Use mean-substitution for missing alleles
      
      int nss0 = snpset[s].size();
            
      //////////////////////////////////////////////
      // Create a vector of pointers for SNPs in set
      
      int nss = 0;
      vector<CSNP*> pSNP(0);
      for (int j=0; j<nss0; j++)
	{
	  // include this SNP? (MC considerations)
	  if ( pS->cur[s][j] ) 
	    {
	      // Add to set list
	      pSNP.push_back(  SNP[snpset[s][j]]  );
	      
	      // Increase the actual number of snps in set
	      nss++; 
	    }
	}
      

      vector<double> mean2(nss,0);    // Case mean
      vector<double> mean1(nss,0);    // Control mean
      vector<vector<double> > pooled; // Covariance matrix

      
      ///////////////////////////////////////
      // Calculate mean and variance (pooled)
      // after imputing missing SNPs
            
      calcHotelSetMeanVariance(pSNP,mean1,mean2,pooled,sample,ncase,ncontrol);
                 
  
      ///////////////////////////////
      // 2. Calculate test statistic
      
      for (int j1=0; j1<nss; j1++)
	for (int j2=j1; j2<nss; j2++)
	  {
	    pooled[j1][j2] /= (double)(ncase+ncontrol-2);

	    pooled[j1][j2] *= 1/(double)ncase + 1/(double)ncontrol;

	    if (j1!=j2)
	      pooled[j2][j1] = pooled[j1][j2];
	    
	  }
      
      // Get inverse of this matix
      bool flag = true;
      pooled = svd_inverse(pooled,flag);

      // Calculate T2 statistic
      vector<double> tmp(nss,0);
      for (int j1=0; j1<nss; j1++)
	for (int j2=0; j2<nss; j2++)
	  tmp[j1] += ( mean1[j2] - mean2[j2] ) * pooled[j1][j2];

      double stat_t2 = 0;
      for (int j1=0; j1<nss; j1++)
	stat_t2 += tmp[j1] * ( mean1[j1] - mean2[j1] );
      
      // T(2) is distributed as (n-1)p / (n-p)  F(p,n-p)
      
      // where n = number of individuals; p = number of variables
      // For two-sample T(2), replace n with n1+n2-1

      // Make statistic; test against F(nss, ncase+ncontrol-1 - nss) 
      
      stat_t2 /= 
	(
	 (double)((ncase+ncontrol-2)*nss)
	 / (double)(ncase+ncontrol-nss-1) 
	 );
      
      // Asymptotic p-value, use 1-p as test statistic
      // for permutation

      double T2p = pF(stat_t2,nss, ncase+ncontrol-nss-1);
      T2[s] = 1 - T2p;
      //      cout << "\nTest STAT = " << T2[s] << "\n";


      if (disp)
	{
	  ASC << setw(12) << setname[s] << " " 
	      << setw(4)  << snpset[s].size() << " " 
	      << setw(12) << stat_t2<< " " 
	      << setw(12) << nss << " "
	      << setw(12) << ncase+ncontrol-nss-1 << " "
	      << setw(12) << T2p << "\n";
	}
     	        
    }

  return T2;
}


void calcHotelSetMeanVariance(vector<CSNP*> & pSNP,
			      vector<double> & mean1,
			      vector<double> & mean2,
			      vector<vector<double> > & pooled,
			      vector<Individual*> & sample,
			      int ncase,
			      int ncontrol)
{

  int nss = mean1.size();

  vector<double> mean(nss,0);

  vector<int> cnt1(nss,0);
  vector<int> cnt2(nss,0);
  
  
  ////////////////////////////
  // Iterate over SNPs in SET
  
  vector<CSNP*>::iterator ps = pSNP.begin();
  int j=0;
  while ( ps != pSNP.end() )
    {
      

      ///////////////////////////
      // Iterate over individuals
      
      vector<Individual*>::iterator gperson = sample.begin();
      vector<bool>::iterator i1 = (*ps)->one.begin();
      vector<bool>::iterator i2 = (*ps)->two.begin();
      int i=0;
      
      while ( gperson != sample.end() )
	{
	  
	  // Permuted self
	  Individual * pperson = (*gperson)->pperson;
	  
	  // Affected individuals
	  if ( ! pperson->missing ) 
	    {
	      if (pperson->aff) 
		{
		  
		  if ( *i1 )
		    {		      
		      if ( *i2 ) // 11 homozygote
			{
			  mean[j]++;
			  cnt2[j]++;
			  mean2[j]++;			  
			}
		    }
		  else 
		    {
		      cnt2[j]++;
		      if ( ! *i2  ) // 00 homozygote
			{
			  mean[j]--;
			  mean2[j]--;			  
			}
		    }
		}
	      else 
		{
		  
		  if ( *i1 )
		    {		      
		      if ( *i2 ) // 11 homozygote
			{
			  mean[j]++;
			  cnt1[j]++;
			  mean1[j]++;			  
			}
		    }
		  else 
		    {
		      cnt1[j]++;
		      if ( ! *i2 ) // 00 homozygote
			{
			  mean[j]--;
			  mean1[j]--;
			}		      			
		    }
		}	  
	    }
	      
	  // Next individual
	  gperson++;
	  i1++;
	  i2++;
	  i++;
	  
	}
      
      // Next SNP in set
      ps++;
      j++;
      
    }
  
  
  
  // Having iterated over all individuals, we can now calculate the mean 
  // values, perform mean-substitution of missing data, and calculate the 
  // second order terms
  
  cout.precision(8);

  for (int j=0; j<nss; j++)
    {
      mean[j] /= (double)(cnt1[j]+cnt2[j]);
      mean1[j] /= (double)cnt1[j];
      mean2[j] /= (double)cnt2[j];
    }

  
  // Element of pooled covariance matrix S is 
  // S[x][y] = ( X[i] - mean1[i] ) ( X[j] - mean1[j] ) 
  
  pooled.resize(nss);
  for (int j=0; j<nss; j++)
    pooled[j].resize(nss,0);
  
  
  /////////////////////////////////////
  // Iterate over pairs of SNPs in SET
  
  // First SNP 
  vector<CSNP*>::iterator ps1 = pSNP.begin();
  int j1=0;
  while ( ps1 != pSNP.end() )
    {
      
      // Second SNP
      vector<CSNP*>::iterator ps2 = ps1;
      int j2=j1;
      while ( ps2 != pSNP.end() )
	{
	  
	  ///////////////////////////
	  // Iterate over individuals
	  
	  vector<Individual*>::iterator gperson = sample.begin();
	  vector<bool>::iterator i1_1 = (*ps1)->one.begin();
	  vector<bool>::iterator i2_1 = (*ps1)->two.begin();
	  
	  vector<bool>::iterator i1_2 = (*ps2)->one.begin();
	  vector<bool>::iterator i2_2 = (*ps2)->two.begin();
	  
	  while ( gperson != sample.end() )
	    {
	      
	      // Permuted self
	      Individual * pperson = (*gperson)->pperson;
	      
	      
	      // Set both values to sample mean
	      double v1 = mean[j1];
	      double v2 = mean[j2];
	      
	      // First SNP
	      if ( *i1_1 )
		{		      
		  if ( *i2_1 )    // 11 homozygote
		    {
		      v1 = 1;
		    }
		}
	      else 
		{
		  if ( ! *i2_1  ) // 00 homozygote
		    {
		      v1 = -1;
		    }
		  else
		    v1 = 0;       // 01 heterozygote
		}
	      
	      
	      // Second SNP
	      if ( *i1_2 )
		{		      
		  if ( *i2_2 )    // 11 homozygote
		    {
		      v2 = 1;
		    }
		}
	      else 
		{
		  if ( ! *i2_2  ) // 00 homozygote
		    {
		      v2 = -1;
		    }
		  else
		    v2 = 0;       // 01 heterozygote
		}
	      
	      
	      // Contribution to covariance term
	      if (! pperson->missing)
		{
		  if (pperson->aff) // affecteds
		    pooled[j1][j2] += ( v1 - mean2[j1] ) * ( v2 - mean2[j2] );
		  else // unaffecteds
		    pooled[j1][j2] += ( v1 - mean1[j1] ) * ( v2 - mean1[j2] );
		}

	      // Next individual
	      gperson++;
	      i1_1++;
	      i2_1++;
	      i1_2++;
	      i2_2++;
	    }
	  
	  // Next second SNP
	  ps2++;
	  j2++;
	}
      
      // Next first SNP
      ps1++;
      j1++;
      
    }
  

  // Make matrix symmetric
  for (int i=0; i<nss; i++)
    for (int j=i; j<nss; j++)
      pooled[j][i] = pooled[i][j];

  return;
}
