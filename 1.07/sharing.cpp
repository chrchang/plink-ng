

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

#include "plink.h"
#include "perm.h"
#include "options.h"
#include "helper.h"

void Plink::perm_sharingIBSTest(Perm & perm)
{

  // This is a SNP-major test
  if (!par::SNP_major) Ind2SNP();

  // Test statistic (set-based)
  vector<double> delta(snpset.size());
  
  // Empirical p-values
  perm.setTests(snpset.size());

  ////////////////////////////////
  // Fast binary affection coding

  if (!par::qt)
    affCoding(*this);


  ////////////////////////////////
  // Set up permutation structure 
  // (we need to perform this step
  //  whether or not we also 
  //  subsequently permute)

  perm.setPermClusters(*this);
  perm.originalOrder();


  /////////////////////
  // Create original
  
  delta = sharingIBSTest(perm);
  

  //////////////////////
  // Begin permutations
  
  bool finished = false;
  while(!finished)
    {
      
      perm.permuteInCluster();

      vector<double> pr = sharingIBSTest(perm);	  
      
      
      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,delta);
     
    } // next permutation

  if (!par::silent)
    cout << "\n\n";
  

  ////////////////////
  // Display results
  
  ofstream ASC;
  string f;
  if (par::adaptive_perm) f = par::output_file_name + ".sharing.perm";
  else f = par::output_file_name + ".sharing.mperm";
  
  ASC.open(f.c_str(),ios::out);
  ASC.precision(4);
  printLOG("Writing IBS sharing association results to [ " + f + " ] \n");
  
  ASC << setw(20) << "SET" << " " 
      << setw(12) << "EMP1" << " ";
  if (par::adaptive_perm)
    ASC << setw(12)<< "NP" << " ";
  else
    ASC	<< setw(12)<< "EMP2" << " ";
  ASC << "\n";
  
  
  for (int l=0; l<snpset.size(); l++)
    {

      ASC << setw(20) << setname[l] << " " 
	  << setw(12) << perm.pvalue(l) << " ";
      
      if (par::adaptive_perm) 
	ASC << setw(12) << perm.reps_done(l) << " ";
      else
	ASC << setw(12) << perm.max_pvalue(l) << " ";
      
      ASC << "\n";
    }

  ASC.close();
     
}

vector<double> Plink::sharingIBSTest(Perm & perm)
{

  
  
//   // number of rare (F) alleles shared
//   // FF FF  ->    
//   // FF FT  ->  
//   // FF TT  ->  -2 
  
//   // FT FF  ->  
//   // FT FT  ->  +1
//   // FT TT  ->  
  
//   // TT FF  ->  -2
//   // TT FT  -> 
//   // TT TT  -> 
  
//   // Number of sets
//   int ns = snpset.size();
  
//   // Test statistics
//   vector<double> delta(ns,0);
  
//   // Iterate over sets
//   for (int i=0; i<ns; i++)
//     {
      
//       // Temporary copy of test statistic
//       double d = 0;

//       // Consider each locus in set
      
//       for (int j=0; j< snpset[i].size(); j++)
// 	{
	  
// 	  CSNP * loc = SNP[snpset[i][j]];
	  
// 	  vector<bool>::iterator a1 = (*loc)->one.begin();
// 	  vector<bool>::iterator a2 = (*loc)->two.begin();
// 	  vector<Individual*>::iterator gperson1 = sample.begin();
// 	  int i1 = 0;

// 	  ////////////////
// 	  // Individual A 
	  
// 	  while ( i1 < n-1 )
// 	    {
	      
// 	      // Permuted self for first individual
// 	      Individual * pperson1 = (*gperson1)->pperson;
	      
// 	      // U-? -- first individual unaffected
// 	      if ( ! ( sample[perm_pheno[i1]]->aff )
// 	      {
		
// 		// Consider each locus: FF x TT -> -2
// 		for (int l=0; l<nl_all; l++)
// 		  {
// 		    // First member FF
// 		    if ( (!g1->one[l]) && (!g1->two[l]) )
// 		      {
// 			// Consider all other individuals
// 			for (int i2=i1+1; i2<n; i2++)
// 			  {
			    
// 			    // is FF x TT ?
// 			    if ( sample[perm_geno[i2]]->one[l] && 
// 				 sample[perm_geno[i2]]->two[l] ) 
// 			{
// 			  // only count discordant pairs
// 			  if ( sample[perm_geno[i2]]->aff )  
// 			    s1[l]-=2;
			  
// 			}
// 		    }
// 		}
// 	      else if ( (!g1->one[l]) && g1->two[l] )
// 		{
// 		  // ... if first member is FT
		  
// 		  // Consider all other individuals
// 		  for (int i2=i1+1; i2<n; i2++)
// 		    {
		      
// 		      Individual * g2 = sample[perm_geno[i2]];
		      
// 		      // FT x FT ->  +1
		      
// 		      if ( (!sample[perm_geno[i2]]->one[l]) && 
// 			   sample[perm_geno[i2]]->two[l] ) 
// 			{
// 			  // Discordant pair
// 			  if ( sample[perm_geno[i2]]->aff )
// 			    s1[l]++;
// 			}
// 		    }
// 		}
// 	      else if ( g1->one[l] && g1->two[l] )
// 		{
// 		  // ... if first member is TT
		  
// 		  // Consider all other individuals
// 		  for (int i2=i1+1; i2<n; i2++)
		    
// 		    {
// 		      if ( (!sample[perm_geno[i2]]->one[l]) && 
// 			   (!sample[perm_geno[i2]]->two[l]) ) 
// 			{
// 			  // Discordant pair
// 			  if ( sample[perm_geno[i2]]->aff )
// 			    s1[l]-=2;
// 			}
// 		    }
// 		}
// 	    } // next locus
	  
// 	}
//       else
// 	{
// 	  // otherwise, we know first individul is affected
// 	  // and so we must now scan for AU and AA pairs (s1, s2)
	  
// 	  // Consider each locus: FF x TT -> -2
// 	  for (int l=0; l<nl_all; l++)
// 	    {
// 	      // First member FF
// 	      if ( (!g1->one[l]) && (!g1->two[l]) )
// 		{
// 		  // Consider all other individuals
// 		  for (int i2=i1+1; i2<n; i2++)
// 		    {
		      
// 		      // is FF x TT ?
// 		      if ( sample[perm_geno[i2]]->one[l] && 
// 			   sample[perm_geno[i2]]->two[l] ) 
// 			{
// 			  // only count discordant pairs
// 			  if ( sample[perm_geno[i2]]->aff )  
// 			    s2[l]-=2; // conc aff
// 			  else
// 			    s1[l]-=2;  // disc
// 			}
// 		    }
// 		}
// 	      else if ( (!g1->one[l]) && g1->two[l] )
// 		{
// 		  // ... if first member is FT
		  
// 		  // Consider all other individuals
// 		  for (int i2=i1+1; i2<n; i2++)
// 		    {
// 		      // FT x FT ->  +1		      
// 		      if ( (!sample[perm_geno[i2]]->one[l]) && 
// 			   sample[perm_geno[i2]]->two[l] ) 
// 			{
// 			  if ( sample[perm_geno[i2]]->aff )
// 			    s2[l]++;
// 			  else 
// 			    s1[l]++;
// 			}
// 		    }
// 		}
// 	      else if ( g1->one[l] && g1->two[l] )
// 		{
// 		  // ... if first member is TT
		  
// 		  // Consider all other individuals
// 		  for (int i2=i1+1; i2<n; i2++)
		    
// 		    {
// 		      if ( (!sample[perm_geno[i2]]->one[l]) && 
// 			   (!sample[perm_geno[i2]]->two[l]) ) 
// 			{
// 			  if ( sample[perm_geno[i2]]->aff )
// 			    s2[l]-=2;
// 			  else
// 			    s1[l]-=2;
// 			}
// 		    }
// 		}
// 	    } // next locus
	  
	

// 	}
	  
      
//     } // next, first individual of pair

  vector<double> t(1);
  return t;
}




