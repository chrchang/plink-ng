

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

vector<double> Plink::testSibTDT(bool print_results, 
				 bool permute,
				 Perm & perm,
				 vector<bool> & flipA,
				 vector<bool> & flipP)
{
  
  ///////////////////////////
  // Vector to store results 

  vector<double> res(nl_all);

  ofstream TDT;
  
  if (print_results)
    {
      string f = par::output_file_name + ".dfam";
      TDT.open(f.c_str(),ios::out);
      printLOG("Writing DFAM results (asymptotic) to [ " + f + " ] \n");
      TDT << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " "
	  << setw(4) << "A1" << " "
	  << setw(4) << "A2" << " "
	  << setw(8) << "OBS" << " "
	  << setw(8) << "EXP" << " "
	  << setw(12) << "CHISQ" << " "
	  << setw(12) << "P" << " ";      
      TDT << "\n";
      
    }
  
  

  ///////////////////////////////////
  // Verbose display of pedigrees

  if ( par::dumpped && ! par::permute ) 
    {
      
      string str = par::output_file_name + ".pdump";
      printLOG("Dumping pedigree information to [ " + str + " ]\n");
      ofstream PD(str.c_str(),ios::out);

      // Of course, in practice, due to missing genotypes, the exact
      // family configurations may shift.

      for (int f=0; f<family.size(); f++)
	{
	  
          Family * fam = family[f];
	  
	  // Families with parents...
          if ( fam->parents && fam->TDT && par::dfam_tdt )
	    {
	      PD << "W/PARENTS\t" << fam->pat->fid << " : ";
              PD << fam->pat->iid << " x " << fam->mat->iid << " -> ";
              for ( int k=0; k<fam->kid.size(); k++)
                PD << fam->kid[k]->iid << " ";
              PD << "\n";

	    }
	  
	  // ...and those sibling without 2 parents
	  else if ( fam->sibship && par::dfam_sibs && fam->kid.size() > 1 )
	    {
	      PD << "SIBSHIP  \t" << fam->kid[0]->fid << " : ";
              for ( int k=0; k<fam->kid.size(); k++)
                PD << fam->kid[k]->iid << " ";
              PD << "\n";
	    }
	}
      
      // And unrelated clusters
      if ( par::dfam_unrelateds )
	for ( int k=0; k<nk; k++)
	  {
	    for (int c=0; c<klist[k]->person.size(); c++)
	      {
		Individual * person = klist[k]->person[c];
		
		if ( ! ( person->family->sibship || 
			 ( person->family->parents ) ) )
		  PD << "CLUSTER " << k << "\t"
		     << klist[k]->person[c]->fid << " : " 
		     << klist[k]->person[c]->iid << "\n";	  
	      }
	  }
      
      PD.close();
    }
  

  ///////////////////////////////////
  // Perform analysis for each locus
  
  for (int l=0; l<nl_all; l++)
    {
      
      // Adaptive permutation, skip this SNP?
      
      if (par::adaptive_perm && (!perm.snp_test[l])) 
      	continue;


      // Skip X/haploid markers for now
      
      if ( par::chr_sex[ locus[l]->chr ] ||
	   par::chr_haploid[ locus[l]->chr ] )
	{
	  continue;
	}
      
      // Allele T counts

      double numerator = 0;
      double denom = 0;
	  
      // Total counts

      double totalCount = 0;
      double totalExpected = 0;
      int totalInformative = 0;
      
      
      /////////////////////////
      // Count over families

      for (int f=0; f<family.size(); f++)
	{
	  
	  // A simple pedigree association test
	  
	  // 1. Break pedigrees into nuclear families
	  
	  // 2. Classify nuclear families into (a) those where both
	  // parents are genotyped (b) the others.
	  
	  // 3. For each type (a) family, obtain the count allele A
	  // among the affected children. and also its expected value
	  // and variance (under H0) given the genotypes of the two
	  // parents. These are given by the binomial distribution,
	  // e.g. if the parental genotypes are (AA, AB) and there are
	  // k affected children, then the the expected count of A is
	  // 1.5k, and its variance is 0.25k.
	  
	  // AA AB ->   AA 0.5    2
	  //            AB 0.5    1
	  //            E = 1.5 * K
	  //            V = 0.25 * K
	  
	  // BB AB ->   AB 0.5    1
	  //            BB 0.5    0
	  //            E = 0.5 * K
	  //            V = 0.25 * K
	  
	  // AB AB ->   AA 0.25   2
	  //            AB 0.50   1
	  //            BB 0.25   0	  
	  //            E = 0
	  //            V = 0.5 * K 
	  
	  
	  // 4. For each type (b) family, also obtain the count A
	  // among the affected children and its expected value and
	  // variance (under H0) given the genotypes of all the
	  // children. These are given by the hypergeometric
	  // distribution, e.g. if the sibship contain n A alleles m B
	  // alleles and there are N affected members, then the
	  // expected value of the count of A among the affected
	  // members would be 2N n/(n+m), and its variance is 2N
	  // nm/(n+m)^2 Note the factor of 2 is because each person
	  // has 2 alleles.
	  
	  // n "A" alleles, m "B" alleles
	  
	  // S = n | k affected siblings
	  // E = 2 k n/(n+m)
	  // V = 2 k nm/(n+m)^2
	  
	  // replacement issue?

	  // 5. An overall test statistic is (Sum of the Counts of A -
	  // Sum of expected counts )^2 / (Sum of variances)
	  


	  Family * fam = family[f];
	  
	  bool informative = false;
	  bool parents = true;
	  

	  // Type A: two genotyped parents and at least 1 affected individual
	  
	  if ( fam->parents && fam->TDT && par::dfam_tdt )
	    {
	      
	      Individual * pat = family[f]->pat;
	      Individual * mat = family[f]->mat;
	      	      
	      bool pat1 = pat->one[l];
	      bool pat2 = pat->two[l];
	      
	      bool mat1 = mat->one[l];
	      bool mat2 = mat->two[l];
	      
	      // We need two genotyped parents, with 
	      // at least one het
	      
	      if ( pat1 && (!pat2) || 
		   mat1 && (!mat2) ) 
		{
		  parents = false;
		  goto jump_to_sibships;
		}
	      
	      int heteroParents = 0;
	      bool homozygParent;
	      
	      if ( pat1 != pat2 )
		heteroParents++;
	      else
		homozygParent = pat1;
	      
	      if ( mat1 != mat2 ) 
		heteroParents++;
	      else
		homozygParent = mat1;
	      
	      if ( heteroParents == 0 )
		{
		  parents = false;
		  goto jump_to_sibships;
		}
	      
	      
	      // Consider all offspring in nuclear family
	      
	      double alleleCount = 0;
	      double childCount = 0;
	      
	      for (int c=0; c<family[f]->kid.size(); c++)
		{
		  
		  // Only consider affected children: based on true
		  // (not permuted) phenotype here: permutation works
		  // by flipping transmissions

		  if ( ! family[f]->kid[c]->aff ) continue;
		  
		  bool kid1 = family[f]->kid[c]->one[l];
		  bool kid2 = family[f]->kid[c]->two[l];
		  
		  // Skip if offspring has missing genotype
		  if ( kid1 && !kid2 ) continue;
		  
		  // We've now established: no missing genotypes
		  // and at least one heterozygous parent
		  
		  if ( permute )
		    {
		      if ( heteroParents == 1 )
			{
			  if ( ! homozygParent )
			    alleleCount++;
			  if ( flipA[f] )
			    alleleCount++;
			}
		      else // ...two heterozygous parents
			{
			  if ( flipA[f] )
			    alleleCount+=2;
			}
		    }
		  else
		    {
		      // No permtutation: standard counting
		      
		      if ( ! kid1 )
			alleleCount++;
		  
		      if ( ! kid2 )
			alleleCount++;		  
		    }

		  childCount++;
		  
		} // next offspring in family
	      

	      double expected = childCount;

	      if ( heteroParents == 1 )
		{
		  if ( ! homozygParent )
		    expected *= 1.5;
		  else
		    expected *= 0.5;
		}
	      
	      double variance = heteroParents * 0.25 * childCount;

 	      numerator += (alleleCount - expected);
  	      denom += variance;
	      	      
	      totalCount += alleleCount;
	      totalExpected += expected;
	      
	    }
	  
	  
	jump_to_sibships:
 
	  // Sibships, considering genotypes using the 
	  // multivariate hypergeometric distribution
	  
	  if ( ( fam->sibship || ( fam->parents && !parents ) ) 
	       && par::dfam_sibs )
	    {
	      
	      // Let the numbers of offspring with genotypes AA, AB
	      // and BB.
	      
	      // N=NAA+NAB+NBB
	      

	      // Let the numbers of affected offsrping with genotypes
	      // AA, AB and BB. be

	      // D=DAA+DAB+DBB
	      
	      // Then the observed count of A is DA = 2*DAA+DAB The
	      // expected value of DA = 2E(DAA) + E D(AB) = 2 NAA *
	      // (D/N) + NAB * (D/N). The variance of DA = 4 Var (DAA)
	      // + Var(DAB) + 4Cov(AA,DAB)
	      
	      // The variances and covariances from the multivariate
	      // hypergeometric distribution.
	      
	      double childCount = 0;
	      double affectedCount = 0;

	      double genotype1Count = 0; // Aa
	      double affectedGenotype1Count = 0;  // Aa

	      double genotype2Count = 0; // AA
	      double affectedGenotype2Count = 0;  // AA
	      
	      for (int c=0; c<family[f]->kid.size(); c++)
		{
		  
		  bool kid1 = family[f]->kid[c]->one[l];
		  bool kid2 = family[f]->kid[c]->two[l];
		  
		  // Skip if offspring has missing genotype
		  if ( kid1 && !kid2 ) continue;
		  
		  // Only consider affected children (possibly allowing
		  // for permutation)
		  

		  if ( family[f]->kid[c]->pperson->aff )
		    {
		      if ( ! kid1 )
			{
			  if ( ! kid2 )			    
			    affectedGenotype2Count++;			   
			  else
			    affectedGenotype1Count++;		 
			}
		      
		      
		      affectedCount++;
		    }
		  
		  
		  if ( !kid1 )
		    {
		      if ( !kid2 )			    
			genotype2Count++;
		      else 
			genotype1Count++;		  
		    }
		  
		  childCount++;
		  
		  
		} // next offspring in family
	      

	      if ( affectedCount > 0 && affectedCount != childCount )
		{
		  
		  double expectedGenotype2 = affectedCount 
		    * ( genotype2Count / childCount ) ; 
		  double expectedGenotype1 = affectedCount 
		    * ( genotype1Count / childCount ) ;
		  
		  double varianceGenotype2 = affectedCount 
		    * ( genotype2Count / childCount )
		    * ( 1 - genotype2Count / childCount )
		    * ( (childCount - affectedCount ) / ( childCount - 1 )  ) ;
		  
		  double varianceGenotype1 = affectedCount 
		    * ( genotype1Count / childCount )
		    * ( 1 - genotype1Count / childCount )
		    * ( (childCount - affectedCount ) / ( childCount - 1 )  ) ;
		  
		  double covarianceGenotype = 
		    - ( affectedCount * 
			( ( genotype2Count * genotype1Count ) 
			  / ( childCount * childCount ) )
			* ( ( (childCount - affectedCount ) 
			      / ( childCount - 1 ) ) ) );		  
		  
		  double affectedAlleleCount = 2 * affectedGenotype2Count 
		    + affectedGenotype1Count;
		  double expected = 2 * expectedGenotype2 
		    + expectedGenotype1;
		  double variance = 4 * varianceGenotype2 
		    + varianceGenotype1 + 4 * covarianceGenotype;

		  numerator += ( affectedAlleleCount - expected );
  		  denom += variance;
		  
 		  totalCount += affectedAlleleCount;
 		  totalExpected += expected;


//        		  cout << "SIB " 
//        		       << childCount << " "
//        		       << genotype1Count << " "
//        		       << genotype2Count << " "
//        		       << affectedCount << " VAR1,2,3= "
//     		       << varianceGenotype2 << " " 
//     		       << varianceGenotype1 << " " 
//     		       << covarianceGenotype << " " 
//        		       << affectedAlleleCount << "  [ "
//        		       << expected << " & "
//        		       << variance << "] \t"
//       		       << numerator << " / " 
//       		       << denom << " , " 
//        		       << totalCount << " "
//        		       << totalExpected << "\n";

		}
	      
	    }
	  

	} // Next nuclear family


      
      ///////////////////////////////////////
      // Now consider clusters of unrelateds
      
      // As sibling test, except allelic rather than genotypic
      // variance estimate (i.e. univariate rather than 
      // multivariate hypergeometic distribution, so equivalent
      // to CMH test

      if ( par::dfam_unrelateds ) 
	for ( int k=0; k<nk; k++)
	  {
	  
	    double affectedCount = 0;
	    double affectedAlleleCount = 0;
		  
	    double alleleCount = 0;
	    double childCount = 0;
	    
	    bool same = true;
	    
	    for (int c=0; c<klist[k]->person.size(); c++)
	      {
		
		Individual * person = klist[k]->person[c];
		
		// Skip any individuals who we've already 
		// analysed as a family
		if ( person->family->sibship ||
		     person->family->parents )
		  continue;
		
		bool s1 = person->one[l];
		bool s2 = person->two[l];
		
		// Skip if offspring has missing genotype
		if ( s1 && !s2 ) continue;
		
		// Are we seeing any genotypic discordance? 
		// Only consider families where we do (i.e. 
		// ignore (het, het) pairs, for example
		
		if ( c>0 && 
		     ( s1 != klist[k]->person[c-1]->one[l] || 
		       s2 != klist[k]->person[c-1]->two[l] ) )
		  same = false;
		
		// Only consider affected children
		if ( person->pperson->aff )
		  {
		    if ( ! s1 )
		      affectedAlleleCount++;
		    
		    if ( ! s2 )
		      affectedAlleleCount++;		  
		    
		    affectedCount++;
		  }
		
		if ( ! s1 )
		  alleleCount++;
		
		if ( ! s2 )
		  alleleCount++;		  
		
		childCount++;
		
	      } // next individual in cluster 
	    		  

	    // S = n | A affected individuals
	    // E = 2 k n/(n+m)
	    // V = 2 k nm/(n+m)^2
	    
	    if ( (!same) && childCount > 1 ) 
	      {
	      
	      double D = alleleCount;
	      double N = 2 * childCount;
	      double A = 2 * affectedCount;
	      
	      double expected = A * ( D / N ) ; 
	      double variance = A *(D/N)*(1-D/N)*((N-A)/(N-1));
	      
	      numerator += ( affectedAlleleCount - expected );
	      denom += variance;
	      
	      totalCount += affectedAlleleCount;
	      totalExpected += expected;
	      
	    }
	  
	} // Next cluster of unrelateds
      
      
      //////////////////////////////
      // Calculate DFAM statistic
      
      double chisq = numerator*numerator;
      chisq /= denom;
      
      
      //////////////////////////////
      // Display asymptotic results
      
      if (print_results)
	{
	  
	  double pvalue = chiprobP(chisq,1);
	  
	  // Skip?, if filtering p-values
	  if ( par::pfilter && pvalue > par::pfvalue ) 
	    continue;

	  TDT.precision(4);
	  
	  TDT << setw(4) << locus[l]->chr << " " 
	      << setw(par::pp_maxsnp) << locus[l]->name << " " 	
	      << setw(4) << locus[l]->allele1 << " "
	      << setw(4) << locus[l]->allele2 << " "
	      << setw(8) << totalCount << " " 
	      << setw(8) << totalExpected << " ";

	  if ( realnum(chisq) )
	    {
	      TDT << setw(12) << chisq << " " 
		  << setw(12) << pvalue << " ";
	    }
	  else
	    TDT << setw(12) << "NA" << " " 
		<< setw(12) << "NA" << " ";
	  
	  TDT << "\n";
	}
      

      /////////////////////////////////
      // Save statistic for permutation

      res[l] = chisq;

      
    } // next locus

  
  //////////////////////////////
  // Close output file, if open
  
  if (print_results)
    TDT.close();

  ///////////////////////////////////////////
  // Return chosen statistic for permutation
  
  return res;

}

