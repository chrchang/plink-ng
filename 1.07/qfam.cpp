

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
#include <map>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "linear.h"
#include "perm.h"
#include "crandom.h"
#include "stats.h"



void setCovariatesForSNP(Plink & P, int l)
{
  vector<Individual*>::iterator gperson = P.sample.begin();
  
  while ( gperson != P.sample.end() )
    {
	 
      // Assume a non-missing genotype
      (*gperson)->flag = true;
      
      bool s1 = (*gperson)->one[l];
      bool s2 = (*gperson)->two[l];	     
      
      if ( ! s1 ) 
	{
	  if ( ! s2 ) 
	    (*gperson)->T = 1;
	  else          
	    (*gperson)->T = 0;
	} 
      else
	{
	  if ( ! s2 ) 
	    (*gperson)->flag = false;
	  else 
	    (*gperson)->T = -1;
	}
      
      // Next individual
      gperson++;  	 
    }
  
}


void scoreBetween(Plink & P , int l)
{
  
  vector<Family*>::iterator f = P.family.begin();
  
  // Construct family B score for this SNP
  
  int fc=0;
  while ( f != P.family.end() )
    {
      
      if (par::verbose)
	{
	  
	  if ( (*f)->singleton )
	    {
	      cout << "SINGLETON(S)\t" 
		   << (*f)->kid[0]->fid << " : ";
	      for (int k=0; k < (*f)->kid.size() ;k++)
		cout << (*f)->kid[k]->iid << " ";
	      cout << "\n";
	    }
	  else if ( (*f)->sibship )
	    {
	      cout << "SIBSHIP  \t" << (*f)->kid[0]->fid << " : ";
	      for ( int k=0; k<(*f)->kid.size(); k++)
		cout << (*f)->kid[k]->iid << " ";
	      cout << "\n";
	    }
	  else if ( (*f)->parents )
	    {
	      cout << "W/ PARENTS\t" << (*f)->pat->fid << " : ";
	      cout << (*f)->pat->iid << " x " << (*f)->mat->iid << " -> ";
	      for ( int k=0; k<(*f)->kid.size(); k++)
		cout << (*f)->kid[k]->iid << " ";
	      cout << "\n";	     
	    }
	  else
	    cout << "UNDEFINED\t" 
		 << (*f)->pat->fid << " "
		 << (*f)->pat->iid << "\n";
	}
      
      
      double B = 0;
      bool Bset = false;
      
      // Include this entire family?
      (*f)->include = true;
      
      // Flag to indicate inclusion in parenQTDT (both parents genotyped)
      (*f)->discordant_parents = true;
      
      // Two theoretically genotyped parents?
      if ( (*f)->parents )
	{
	  // Two actually genotyped parents?
	  if ( (*f)->pat->flag && (*f)->mat->flag )
	    {
	      B = ( (*f)->pat->T + (*f)->mat->T ) * 0.5 ;
	      Bset = true;
	    }
	  else
	    (*f)->discordant_parents = false;
	}
      
      
      
      // Did this individual have parental genotype information to set B? If not...
      
      if ( !Bset ) 
	{
	  // Use sibling genotypes? This will default to one's own
	  // genotype (i.e. sibship of size 1 (singletons are coded 
	  // offspring here)
	     
	  // Number of sibling
	  int nsib = (*f)->kid.size();
	  
	  // Number of genotyped sibling
	     int nsib2 = nsib;
	     
	     for (int k=0; k<nsib; k++)
	       {
		 if ( (*f)->kid[k]->flag )
		   B += (*f)->kid[k]->T;
		 else
		   nsib2--;
	       }
	     
	     if (nsib2==0)
	       {
		 // No non-missing offspring in family
		 // so does not matter what we set here
		 (*f)->include = false;		 
		 B = -9;
	       }
	     else
	       B /= (double)nsib2;
	}
      
      // Store between family score
      (*f)->B = B;
      
      // Next family
      f++;
      fc++;
    }
}


void scoreBandW(Plink & P, int l , vector<bool> & include) 
{
  
  // Initially, everybody is included

  vector<Individual*>::iterator gperson = P.sample.begin();
  int i=0;
     
  while ( gperson != P.sample.end() )
    {
	 
      Individual * pperson = (*gperson)->pperson;

      if ( ! pperson->family ) 
	error("Internal problem: no family assigned for [ " 
	      + pperson->fid + " " + pperson->iid + " ]\n");

      // Valid phenotype...
      if ( ( ! pperson->missing ) )
	{
	  // ... and genotype?
	  if ( (*gperson)->flag )
	    {
	      
	      Family * f = (*gperson)->family;
	      
	      // Are we modelling parental phenotypes?
	      if ( par::QFAM_total ||              // total association test...
		   par::QFAM_between ||            // ...between association test...
		   ( par::QFAM_within2     
		     && f			
		     && f->discordant_parents ) || // ...parenQTDT...
		   ! (*gperson)->founder )         // ...or, not a founder
		{
		  
		  // Between-family component
		  (*gperson)->B = f ? f->B : 0;
		  
		  // Within-family component
		  (*gperson)->W = (*gperson)->T - (*gperson)->B;
		  
		}
	      else
		include[i] = false;	      
	    }
	  else include[i] = false;
	}
      else include[i] = false;
      
      // Next person
      gperson++;
      i++;
    }
     
}

//////////////////////////////////////////////////////////////////////
//
// For QFAM, and unlike all other tests, we use two different ways of
// permuting: either standard (all SNPs per replicate) or on a per SNP
// adaptive basis (i.e. all perms for a SNP; then move on to next
// SNP). This saves the work of constructing the family, etc, as we
// need to do each time.
//


void Plink::perm_testQTDT(Perm & perm)
{

  //////////////////////////////
  // Use individual-major coding
  
  if (par::SNP_major) 
    SNP2Ind();
  
  
  // for now, no covariates
  if ( par::clist_number > 0 ) 
    error("Cannot specify covariates with QFAM for now...\n");


  ////////////////////////////////////////////////
  // Specify special adaptive QFAM mode (i.e. one SNP
  // at a time)



  /////////////////////////////
  // Set up permutation indices
  
  vector<int> pbetween(family.size());
  vector<bool> pwithin(family.size(),false);
  for (int i=0; i < family.size(); i++)
    pbetween[i] = i;
  
  
  ///////////////
  // Output files

  string f = ".qfam";
  if (par::QFAM_within1) f += ".within";
  else if (par::QFAM_within2) f += ".parents";
  else if (par::QFAM_between) f += ".between";
  else if (par::QFAM_total) f += ".total";
  
  printLOG("Writing QFAM statistics to [ " + par::output_file_name + f + " ]\n");
  
  if (!par::permute) 
    printLOG("** Warning ** QFAM results require permutation to correct for family structure\n");
  else
    printLOG("Important: asymptotic p-values not necessarily corrected for family-structure:\n"
	     "           use empirical p-values for robust p-values from QFAM\n"
	     "           and consult the above file only for parameter estimates\n");
  

  ofstream QOUT((par::output_file_name+f).c_str(),ios::out); // dummy
  QOUT.precision(4);
  QOUT << setw(4) << "CHR" << " " 
       << setw(par::pp_maxsnp) << "SNP" << " " 
       << setw(10) << "BP" << " "
       << setw(4) << "A1" << " "
       << setw(10) << "TEST" << " "
       << setw(8) << "NIND" << " "
       << setw(10) << "BETA" << " ";
  if (par::display_ci)    
    QOUT << setw(8) << "SE" << " "
	 << setw(8) << "LOWER" << " "
	 << setw(8) << "UPPER" << " ";	    	      
  QOUT << setw(12) << "STAT" << " "
       << setw(12) << "P\n";


  //////////////////////
  // Familial clustering

  // C holds which family an individual belongs to 
  // (as element in the family[] array
  
  vector<int> C;
  map<Family*,int> famcnt;
  for (int f = 0 ; f < family.size() ; f++)
    famcnt.insert( make_pair( family[f] , f ) );      
  
  vector<Individual*>::iterator person = sample.begin(); 
  while ( person != sample.end() )
    {
      map<Family*,int>::iterator f = famcnt.find( (*person)->family );
      
      if ( f == famcnt.end() )
	error("Internal error in QFAM, allocating families to individuals...\n");
      else
	C.push_back( f->second );
      
      person++;  
    } 
  
  printLOG(int2str(family.size())+" nuclear families in analysis\n");
      
  if ( family.size()<2 )
    error("Halting: not enough nuclear families for this analysis\n");
  


  ////////////////////
  // Run original QFAM

  perm.setTests(nl_all);
  perm.setPermClusters(*this);

  // Force adaptive perm
  par::adaptive_perm = true;

  vector_t orig = calcQTDT(C, QOUT, false, perm, pbetween, pwithin);

  QOUT.close(); 



  ////////////////
  // Permutation

  if ( ! par::permute ) 
    return;
  
  // Adpative permutation will already have been conducted in original 
  // function call for QFAM (i.e. per-SNP adaptive permutation)

  if (!par::silent)
    cout << "\n\n";
  

  ////////////////////
  // Display results
  
  ofstream TDT;   

  f += ".perm";    
  TDT.open((par::output_file_name+f).c_str(),ios::out);
  printLOG("Writing QFAM permutation results to [ " 
	   + par::output_file_name + f + " ] \n"); 
  TDT.precision(4);
  
  TDT << setw(4) << "CHR" << " "
      << setw(par::pp_maxsnp) << "SNP" << " ";
  
  if (par::perm_TDT_basic) TDT << setw(12) << "STAT" << " ";
  
  TDT << setw(12) << "EMP1" << " ";
  TDT << setw(12) << "NP" << " " << "\n";  
  
  for (int l=0; l<nl_all; l++)
    {	
      
      TDT << setw(4) << locus[l]->chr << " "
	  << setw(par::pp_maxsnp) << locus[l]->name << " "; 
      
      if (orig[l] < -0.5)
	TDT << setw(12) << "NA"  << " " 
	    << setw(12) << "NA"  << " " 
	    << setw(12) << "NA";
      else
	{
	  TDT << setw(12) << orig[l] << " "
	      << setw(12) << perm.pvalue(l) << " "
	      << setw(12) << perm.reps_done(l);	  
	}
      TDT << "\n";
    }
  
  TDT.close();

  
  // Adjusted p-values, assumes 1-df chi-squares
  
  if (par::multtest)
    {
      
      vector<double> obp(0);
      for (int l=0; l<nl_all;l++)
	obp.push_back(inverse_chiprob(perm.pvalue(l),1));      
      
      multcomp(obp,f);
    }

  
   
}

vector_t Plink::calcQTDT(vector<int> & C,
			 ofstream & QOUT,
			 bool permuting, 
			 Perm & perm,
			 vector<int> & pbetween, 
			 vector<bool> & pwithin)
{
  

  /////////////////////////
  // Iterate over each SNP
  
  vector_t results(nl_all);
  
  for (int l=0; l<nl_all; l++)
    {     
      
      // Note: when using adaptive permutation in QFAM, we do not skip
      // a failed SNP here, as we permute on a per-SNP basis instead;
      // i.e. for this particular SNP we will perform enough
      // permutations to assess significance in this first instance of the 
      // call to calcQTDT().  
      
      // Skip X markers for now
      
      if (par::chr_sex[locus[l]->chr] || 
	  par::chr_haploid[locus[l]->chr])
	{
 	  results[l] = -1;
 	  continue;
	}
      
      if (par::verbose)
	cout << "\n ******************************************\n"
	     << "  LOCUS " << locus[l]->name << "\n\n";
      
      
      ////////////////////////////////////////////////////////////////
      // Create X vector that encodes the genotype for each individual
      // as 1,0,-1 (or -9 for missing)
      
      // Use the per-person 'flag' variable to indicate a non-missing genotype
      // at this SNP (i.e. for gperson)
      
      // Use 'covar' to store the X= 1,0,-1 codes for this SNP
      
      setCovariatesForSNP(*this,l);
      
      
      ///////////////////////////////////////
      // Score between and within components
	
      scoreBetween(*this,l);
	
      // Now, for each individual, set B and W 

      vector<bool> include(n,true);
      
      scoreBandW(*this,l,include);
     

      // Now we have created the family structure, B and W and flagged who is missing
      // in terms of genotype and phenotype
      
      // We can either proceed to return one value for this (in max(T) mode)
      // or to exhaust all permutations
      
     
      /////////////////////////
      // Prune out missing data (already done?)
     
      vector<Family*>::iterator f = family.begin();
      while ( f != family.end() ) 
        {
 	 if ( ! (*f)->include ) 
 	   {
 	     if ( (*f)->pat ) 
 	       (*f)->pat->flag = false;
	     
 	     if ( (*f)->mat ) 
 	       (*f)->mat->flag = false;
		 
 	     for ( int k = 0 ; k < (*f)->kid.size() ; k++) 
	       (*f)->kid[k]->flag = false;
	     
 	   }
 	 f++;
        }
     
     
      // Prune individuals
      for (int i=0; i<n; i++)
        if ( (!sample[i]->flag) || sample[i]->missing ) 
	  include[i] = false;
                    

     /////////////////////////
     // Optional display
     
     if (par::verbose)
       {
	 
	 for (int i=0; i<n; i++)
	   {
	     if ( include[i] ) 
	       cout << "INC\t";
	     else
	       cout << "EXC\t";
	     
	     cout << C[i] << "\t"
		  << sample[i]->fid << " " << sample[i]->iid << "\t"
		  << sample[i]->phenotype << "\t"
		  << genotype(*this,i,l) << " "
		  << sample[i]->T << " "
		  << sample[i]->B << " "
		  << sample[i]->W ;
	     cout << "\n";
	   }
	 cout << "\n\n";
       }
     
     
     
     ///////////////////////////////////
     // Form linear model
     
     Model * lm;
     LinearModel * m = new LinearModel(this);
     lm = m;
     
     // Copy pattern of missing data over, with 
     // some additional exclusions based on family 
     // structure
     
     lm->setMissing(include);
     
     // Add independent variables: T, B and/or W
     // and set the test parameter
     // (intercept is 0)
     
     // Covariates  Model
     // 0 Total
     // 1 Between
     // 2 Within
     
     // Model
     // 0      Intercept      Intercept
     // 1      Total          Between
     // 2      n/a            Within
     
     if (par::QFAM_total) 
       {
	 lm->label.push_back("TOT");
	 lm->testParameter = 1;	     
       }
     else if (par::QFAM_between)
       {
	 lm->label.push_back("BET");
	 //	 lm->label.push_back("WITH");
	 lm->testParameter = 1;	     
       }
     else if (par::QFAM_within1 || par::QFAM_within2) 
       {
	 //	 lm->label.push_back("BET");
	 lm->label.push_back("WITH");
	 lm->testParameter = 1;	     
       }
     
     // Build design matrix
     lm->buildDesignMatrix();
     
     // Fit linear model
     if ( par::QFAM_total && par::qt )
       lm->fitUnivariateLM();
     else
       lm->fitLM();

     // Check for multi-collinearity
     lm->validParameters();
     
     // Calculate Original Test statistic
     results[l] = lm->getStatistic();

     // Store,return and display this value?
     
     lm->displayResults(QOUT,locus[l]);

     
     ///////////////////
     // Now, permutation
     
     // 1) We have the complete, non-missing data: permute only this
     //    i.e. we do not need to worry about missing data; we are
     //    no longer controlling the correlation between SNPs, as we
     //    are permuting genotype, so we do not need to worry about this
     //    in any case.
     
     // 2) Keep the same Model in each case: directly re-state the X 
     //    variables in the design matrix, then re-fit model. This 
     //    will avoid the cost of building the model, pruning for missing
     //    data, etc, each iteration
     
     
     // Store original, and set up permutations
     // (i.e. return pperson to original order)

     perm.nextSNP();
     double original = results[l];
     

     ////////////////////////
     // Adaptive permutation
     
     ///////////////////////////////////////////////////
     // Set up permutation indices, specific to this SNP
     
     int tc = 0;

     while ( true ) 
       {
	 
	 // Permute between and within family components
	 
	 permute(pbetween);
	 
	 for (int i=0; i<family.size(); i++)
	   {
	     if (CRandom::rand() < 0.5) pwithin[i] = true;
	     else pwithin[i] = false;
	   }
	 
	 // Edit pbetween for this SNP, so that we keep missing 
	 // B components constant
	 
	 for (int f=0; f<pbetween.size(); f++)
	   {
	     
	     if ( 
		 // Permuted family is all missing
		 ( ! family[pbetween[f]]->include ) 
		 &&
		 // Recipient family is not...
		 family[f]->include )
	       {
		 // ... then swap 
		 
		 //   F  P(F)    -->
		 //   0  2       -->   0  2
		 //   1  0       -->   1  0
		 //   2  3*      -->   2  4
		 //   3* 4       -->   3* 3*
		 //   4  1       -->   4  1
		 //   ...
		 
		 // e.g. 3* is missing, so swap 3* and 4 in P(F), so 2
		 // and 4 end up together instead, 3* is invarint
		 
		 int missing_family = pbetween[f];
		 int swap_in_family = pbetween[pbetween[f]];
		 pbetween[missing_family] = missing_family;
		 pbetween[f] = swap_in_family;
		 
// 		 if (par::verbose)
// 		   {
// 		     cout << "FAM " << f << " (NOT MISS) has " << missing_family << " (MISS)\n";
// 		     cout << "FAM " << missing_family << " (MISS) has " << swap_in_family << " (?)\n";
// 		     cout << "SWAP MADE ..\n";
// 		     cout << "FAM " << f << "  has " << pbetween[f] << "\n";
// 		     cout << "FAM " << missing_family << " has " << pbetween[missing_family] << "\n\n";
// 		   }
		 
		 // And re-check this new pairing
		 f--;
		 
	       }		 
	   }
	 

// 	 if (par::verbose)
// 	   for (int f=0; f<pbetween.size(); f++)
// 	     {
// 	       if ( ! family[pbetween[f]]->include ) 
// 		 cout << " Permuted family is all missing " << f << "\t" << family[pbetween[f]]->kid[0]->fid << "\n";
// 	       if ( ! family[f]->include ) 
// 		 cout << " Recipient family is all missing " << f << "\t" << family[f]->kid[0]->fid << "\n";
// 	     }


//  	 if (true)
//  	   {
//  	     for (int i=0; i<n; i++)
//  	       {
//  		 cout << sample[i]->fid << "\t"
//  		      << include[i] << "\t"
//  		      << C[i] << "\t"
//  		      << pbetween[C[i]] << "\t"
//  		      << sample[i]->family->include << "\t"
//  		      << family[C[i]]->include << "\t"
//  		      << family[pbetween[C[i]]]->include << "\n";
//  	       }
//  	   }




	 //////////////////////////////////
	 // Reconstitute genotypes
	 // and fit back into LinearModel
	 
 	 int c=0;
 	 for (int i=0; i<n; i++)
 	   {
 	     if (include[i])
 	       {
 		 Family * pfam = family[ pbetween[C[i]] ];
 		 Individual * person = sample[i];
		 
 		 if ( par::QFAM_total )
 		   lm->X[c++][1] = pwithin[C[i]] ? pfam->B + person->W : pfam->B - person->W;
 		 else if ( par::QFAM_between )
 		   {
 		     lm->X[c++][1] = pfam->B;
 		   }
 		 else
 		   {
 		     lm->X[c++][1] = pwithin[C[i]] ? person->W : - person->W;
 		   }
		 
// 		 cout << "added " << person->fid << " " 
// 		      << person->iid << " " 
// 		      << lm->X[c-1][1] << "\n";


 	       }
 	   }
// 	 cout << "\n\n";



		 
	 ////////////////////////////////////
	 // Re-fit model
	 
	 if ( par::QFAM_total && par::qt )
	   lm->fitUnivariateLM();
	 else
	   lm->fitLM();
	 
	 // Check for multi-collinearity
	 lm->validParameters();
	 

	 // Calculate Original Test statistic; 

	 // Should not encounter this too much, but if not valid,
	 // count conservatively.

	 double r = lm->isValid() ? lm->getStatistic() : original + 1 ;
	 
// 	 cout << "Permutation ... \n";
// 	 if ( ! lm->isValid() )
// 	   cout << "NOT VALID>.. \n";

// 	 int c2 = 0;

// 	 for (int i=0; i<n; i++)
// 	   {
// 	     if ( include[i] ) 
// 	       cout << "INC\t";
// 	     else
// 	       cout << "EXC\t";
	     	     
// 	     cout << C[i] << "\t"
// 		  << sample[i]->fid << " " << sample[i]->iid << "\t"
// 		  << sample[i]->phenotype << "\t"
// 		  << genotype(*this,i,l) << " ";

// 	     if ( include[i] )
// 	       cout << lm->X[c2++][1] << " ";
// 	     else 
// 	       cout << "NA" << " ";
// 	     cout << "\n";
// 	   }
// 	 cout << "\n\n";






	 // Reset in case the previous model was not valid

	 lm->setValid();
	 
	 ////////////////////////////////////
	 // Test / update / are we finished ? 
	 
	 if ( perm.updateSNP( r , original , l ) )
	   {
	     if ( ! par::silent )
	       {
		 cout << "Adaptive permutation done for " 
		      << l+1 << " of " << nl_all << " SNPs            \r";
		 cout.flush(); 
	       }
	     break; // We are done for this SNP
	   }
	 
	 
       } // Next adaptive permutation

      
     // Clear up
     delete lm;
          
   } // Next SNP
 
 return results;
 
}


