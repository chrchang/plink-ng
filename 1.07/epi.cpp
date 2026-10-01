

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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "crandom.h"
#include "linear.h"
#include "logistic.h"
#include "stats.h"

extern ofstream LOG;

using namespace std;


////////////////////////////////////////
// Epistasis tests (no permutation)

void Plink::calcEpistasis()
{

  ///////////////////////////////////////////
  // SNP major mode or individual major mode?

  if (par::fast_epistasis)
    {
      if ( ! par::SNP_major ) 
	Ind2SNP();
    }
  else
    {
      if ( par::SNP_major ) 
	SNP2Ind();
    }


  //////////////////////////////////////////////
  // Set up results files

  ofstream EPI;
  string f = par::output_file_name;
  if (par::qt) 
    f += ".epi.qt";
  else 
    {
      if (par::epi_caseonly) f += ".epi.co";
      else f += ".epi.cc";
    }
  EPI.open(f.c_str(),ios::out);
  printLOG("Writing epistasis pairwise results to [ " + f + " ] \n");

  EPI.precision(4);

  if ( !par::epi_quickscan )
    {
      EPI << setw(4) << "CHR1" << " "
	  << setw(par::pp_maxsnp) << "SNP1" << " "
	  << setw(4) << "CHR2" << " "
	  << setw(par::pp_maxsnp) << "SNP2" << " ";
      
      if (!par::fast_epistasis)
	{
	  if (par::bt)
	    EPI << setw(12) << "OR_INT" << " ";
	  else
	    EPI << setw(12) << "BETA_INT" << " ";
	}
      
      
      EPI << setw(12) << "STAT" << " "
	  << setw(12) << "P" << " "
	  << "\n";
    }
  else
    EPI << setw(4) << "CHR1" << " "
	<< setw(par::pp_maxsnp) << "SNP1" << " "
	<< setw(4) << "CHR2" << " "
	<< setw(par::pp_maxsnp) << "SNP2" << " "
	<< setw(12) << "CHISQ" << " "
	<< "\n";
  
  
  ////////////////////////////////////////////////////////////////////
  // epi1 and epi2 thresholds were given in terms of 0.01 (two-sided)
  // calculate appropriate absolute Z scores
  
  printLOG("Threshold for displaying epistatic result (--epi1) : p <= "+dbl2str(par::epi_alpha1)+"\n");
  printLOG("Threshold for counting epistatic result (--epi2) : p <= "+dbl2str(par::epi_alpha2)+"\n");
  
  par::epi_alpha1 = fabs(ltqnorm(par::epi_alpha1 / 2));
  par::epi_alpha2 = fabs(ltqnorm(par::epi_alpha2 / 2));
  
  // Fast epistasis:  caae-only or case/control
  // Regression based test: case/control or quantitative trait

  // Take a list of SNPs, or all SNPs (vector<bool> epi1)
  // Test these against either themselves, or all SNPs (vector<bool> epi2)
  
  //  A     B
  //  ALL x ALL    skip e1>e2
  //  SET1 x ALL
  //  SET1 x SET1  skip e1>e2
  //  SET1 x SET2

  bool skip_symm = false;
 
  // Only output epistatic tests that have p < par::epi_alpha1;
  // Do not even attempt to save any epistatic results -- go straight to STDOUT
  
  // Also present summary results for all epi1 SNPs 
  // (i.e. average / proportion of significant epistatic tests
  //  at a certain alpha level, par::epi_alpha2


  vector<bool> sA(nl_all,false);
  vector<bool> sB(nl_all,false);
  
  // Are we using a test set? If so, construct now
  if (par::set_test) 
    {

      if (snpset.size()>2) 
	error("Can only specify one or two SETs when testing for epistasis");
      if (snpset.size()==0) 
	error("There are no valid sets specified");

      for (int e=0;e<snpset[0].size();e++)
	sA[snpset[0][e]] = true;
      
      // Has a second set been specified?
      
      if (snpset.size()==2)
	{
	  printLOG("SET1 x SET2 epistasis mode\n");
	  for (int e=0;e<snpset[1].size();e++)
	    sB[snpset[1][e]] = true;
	}
      else if (par::set_by_set) // Otherwise, has SET x SET flag been given? 
	{
	  printLOG("SET1 x SET1 epistasis mode\n");
	  skip_symm = true;
	  for (int e=0;e<snpset[0].size();e++)
	    sB[snpset[0][e]] = true;	  
	}
      else // All SNPs in second set
	{
	  printLOG("SET1 x ALL epistasis mode\n");
	  for (int e=0;e<nl_all;e++)
	    sB[e] = true;	  
	}
    }
  else
    {
      printLOG("ALL x ALL epistasis mode\n");
      skip_symm = true;
      for (int e=0;e<nl_all;e++)
	{
	  sA[e] = true;
	  sB[e] = true;	  
	}
    }
  

  // Use fast aff coding
  
  if (par::bt)
    affCoding(*this);
  
  
  // Count how many items in the SET1

  int epc = 0;
  for (vector<bool>::iterator e1 = sA.begin(); 
       e1 != sA.end();
       e1++)
    if (*e1) epc++;
  int epcc = 0;


  // Keep track of how many epistatic tests actually performed
  long int nepi = 0;

  vector<int> summary_sig(nl_all,0);  
  vector<int> summary_good(nl_all,0);  
  vector<double> best_score(nl_all,0);
  vector<int> best_partner(nl_all);
  

  //////////////////////////////////////////
  // Begin iterating over pairs : SET x SET
  
  for (int e1=0;e1<nl_all;e1++)
    {
      if (sA[e1]) 
	{
	  if (!par::silent)
	    {
	      cout << "Peforming tests of epistasis: group " 
		   << ++epcc << " of " << epc << "        \r";
	      cout.flush();
	    }
	  
	  for (int e2=0;e2<nl_all;e2++)
	    {

	      ///////////////////////////////////////////
	      // Skip this test under certain conditions

	      // The SNP not in the set
	      if (!sB[e2]) { cout << "skipping...\n"; continue; }

	      // We've already performed this test
	      if (e1>=e2 && skip_symm) continue;

	      // Same SNP 
	      if (e1==e2) continue;

	      // Skip X chromosome for now
	      if (par::chr_sex[locus[e1]->chr] || 
		  par::chr_sex[locus[e2]->chr] || 
		  par::chr_haploid[locus[e1]->chr] ||
		  par::chr_haploid[locus[e2]->chr]) 
		continue;
	      
	      // SNPs too close (case-only analysis)
	      if (par::epi_caseonly) 
		if ( locus[e1]->chr == locus[e2]->chr)
		  if ( fabs((double)(locus[e1]->bp - locus[e2]->bp)) 
		       < par::epi_caseonly_kb_gap*1000 )
		    continue;
	      


	      //////////////////////////////////
	      // Perform test of epistasis here
	      
	      if (par::bt && par::fast_epistasis)
		{
		  
		  double z;  // statistic from either method
		  
		  // Odds ratio test
		  // make two 2x2 tables
		  
		  int a11, a12, a21, a22;
		  int u11, u12, u21, u22;
		  a11=a12=a21=a22=0;
		  u11=u12=u21=u22=0;
		  
		  vector<bool>::iterator a1 = SNP[e1]->one.begin();
		  vector<bool>::iterator a2 = SNP[e1]->two.begin();
		  
		  vector<bool>::iterator b1 = SNP[e2]->one.begin();
		  vector<bool>::iterator b2 = SNP[e2]->two.begin();
		  
		  vector<Individual*>::iterator person = sample.begin();
		  
		  while ( person != sample.end() )
		    {
		      
		      if( (*person)->missing )
			{
			  // Next person
			  a1++;
			  a2++;
			  b1++;
			  b2++;
			  person++;
			  continue;
			}
		      
		      if ((*person)->aff) // if affected 
			{	      
			  
			  if ( ! *b1 )
			    {
			      if ( ! *b2  )  //   ??x00
				{
				  if ( ! *a1 )
				    {
				      if ( ! *a2 )
					a11+=4; // 00 x 00 
				      else
					{ a11+=2; a21+=2; } // 01 x 00
				    }
				  else if ( *a2 )
				    a21+=4;  // 11 x 00
				}
			      else  //   ??x01
				{
				  if ( ! *a1 )
				    {
				      if ( ! *a2 )
					{ a11+=2; a12+=2; }  // 00 x 01 
				      else
					{ a11++; a21++; a12++; a22++; } // 01x01
				    }
				  else if ( *a2 )
				    { a21+=2; a22+=2; }  // 11 x 01
				}
			    }
			  else if ( *b2 ) // ?? x 11 
			    {
			      
			      if ( ! *a1 )
				{
				  if ( ! *a2 )
				    a12+=4;    // 00 x 01 
				  else
				    { a12+=2; a22+=2; } // 01 x 01
				}
			      else if ( *a2 )
				a22+=4;  // 11 x 01
			      
			    }
			}     
		      
		      // Unaffecteds?
		      else if ( !par::epi_caseonly )  // unaffected 
			{
			      
			  if ( ! *b1 )
			    {
			      if ( ! *b2 )  //   ??x00
				{
				  if ( ! *a1 )
				    {
				      if ( ! *a2 )
					u11+=4;   // 00 x 00 
				      else
					{ u11+=2; u21+=2; }  // 01 x 00
				    }
				  else if ( *a2 )
				    u21+=4;  // 11 x 00
				}
			      else  //   ??x01
				{
				  if ( ! *a1 )
				    {
				      if ( ! *a2 )
					{ u11+=2; u12+=2; }  // 00 x 01 
				      else
					{ u11++; u21++; u12++; u22++; } // 01x01
				    }
				  else if ( *a2 )
				    { u21+=2; u22+=2; } // 11 x 01
				}
			    }
			  else if ( *b2 )  //  ?? x 11 
			    {
			      
			      if ( ! *a1 )
				{
				  if ( ! *a2 )
				    u12+=4;  // 00 x 01 
				  else
				    { u12+=2; u22+=2; }  // 01 x 01
				}
			      else if ( *a2 )
				u22+=4;   // 11 x 01
			      
			    }
			  
			}
		      
		      // Next person
		      a1++;
		      a2++;
		      b1++;
		      b2++;
		      person++;
		      
		    }
		  
		  
		  // Calculate log(OR) and SEs
		  
		  double or_aff, v_aff, or_unf, v_unf;
		  
		  or_aff = log( (double)(a11*a22)/ (double)(a12*a21) );
		  v_aff = 1/(double)a11 + 1/(double)a12 + 1/(double)a21 + 1/(double)a22;
		  
		  // Case-only z-score (if requested)
		  if (par::epi_caseonly)
		    z = fabs( or_aff / sqrt(v_aff) );
		  else // Standard case-control analysis 
		    {
		      or_unf = log( (double)(u11*u22)/ (double)(u12*u21) );
		      v_unf = 1/(double)u11 + 1/(double)u12 + 1/(double)u21 + 1/(double)u22;
		      z = fabs( (or_aff - or_unf) / sqrt ( v_aff + v_unf ) );
		    }
		  

		  //////////////////////////////
		  // --nop option in effect 
		  // Just output z score, if valid & above threshold
		  
		  if (par::epi_quickscan)
		    {
		      // Is this worth recording?
		      if ( realnum(z) )
			{
			  nepi++;
			  if (z >= par::epi_alpha1)
			    EPI << setw(4) << locus[e1]->chr << " " 
				<< setw(par::pp_maxsnp) << locus[e1]->name << " "
				<< setw(4) << locus[e2]->chr << " " 
				<< setw(par::pp_maxsnp) << locus[e2]->name << " "
				<< setw(12) << z*z << "\n";
			  EPI.flush();
			  continue;
			}
		    }
		  

		  /////////////////////////////////	
		  // More full parsing of results
		  
		  double zero = 0;

		  // Check this is a proper result
		  		  
		  if ( par::epi_filter && realnum(z) )
		    {  
		      		      
		      // One more test performed
		      nepi++;
		      
		      // Count as a good result
		      summary_good[e1]++;
		      if (sA[e2]) summary_good[e2]++;
		      
		      // Do we want to record this as part of the summary for the first set? 	  
		      if (z >= par::epi_alpha2)
			{
			  // first variable will always be in A set
			  summary_sig[e1]++;
			  // but the second may also be in A set
			  if (sA[e2]) summary_sig[e2]++;
			}
		      
		      // Is this result the best scrore yet for marker in set A?
		      if (z > best_score[e1]) 
			{
			  best_score[e1] = z;
			  best_partner[e1] = e2;
			}
		      
		      // The second marker might also be in set A 
		      if (sA[e2])
			{
			  if (z > best_score[e2])
			    {
			      best_score[e2] = z;
			      best_partner[e2] = e1;
			    }
			}
		      
		      // Is this worth recording?
		      
		      if (z >= par::epi_alpha1)
			{
			  EPI << setw(4) << locus[e1]->chr << " " 
			      << setw(par::pp_maxsnp) << locus[e1]->name << " "
			      << setw(4) << locus[e2]->chr << " "
			      << setw(par::pp_maxsnp) << locus[e2]->name  << " "
			      << setw(12) << z*z   << " "
			      << setw(12) << normdist(-z) * 2  << " "
			      << "\n";
			  EPI.flush();
			}
		      else
			continue; // skip to next pair (skip logistic test)
		      
		    }
		  else if (!par::epi_filter)
		    {
		      // Record all results here, whether NA or otherwise
		      EPI << setw(4) << locus[e1]->chr << " " 
			  << setw(par::pp_maxsnp) << locus[e1]->name << " "
			  << setw(4) << locus[e2]->chr << " "
			  << setw(par::pp_maxsnp) << locus[e2]->name  << " "
			  << setw(12) << z*z   << " "
			  << setw(12) << normdist(-z) * 2  << " "
			  << "\n";
		      EPI.flush();
		    }
		  else
		    continue; // if bad statistic for this test, do not try logistic
		  
		} // End of binary OR test
	      
	      
	      ///////////////////////////////////////////////
	      // Logistic or linear regression epistasis test
	      
	       if ( !par::fast_epistasis )
		 {

		   Model * lm;

		   if (par::bt)
		     {
		       LogisticModel * m = new LogisticModel(this);
		       lm = m;
		     }
		   else
		     {
		       LinearModel * m = new LinearModel(this);
		       lm = m;
		     }

		   // Set missing data

		   lm->setMissing();

		   
		   // Main effect of SNP 1

		   lm->addAdditiveSNP(e1); 
		   lm->label.push_back("ADD1");

		   
		   // Main effect of SNP 2

		   lm->addAdditiveSNP(e2); 
		   lm->label.push_back("ADD2");


		   // Epistasis

		   lm->addInteraction(1,2);
		   lm->label.push_back("EPI");	  

		   
		   // Build design matrix

		   lm->buildDesignMatrix();

		   
		   // Prune out any remaining missing individuals
		   // No longer needed

//		   lm->pruneY();


		   // Fit linear model

		   lm->fitLM();


		   // Did model fit okay?

		   lm->validParameters();


		   // Obtain estimates and statistic

		   lm->testParameter = 3; // interaction
		   vector_t b = lm->getCoefs();
		   double chisq = lm->getStatistic();
		   double pvalue = chiprobP(chisq,1);
		   double z = sqrt(chisq);

		   // Is this result worth displaying?
		   
		   if (lm->isValid())
		     {
		       
		       // One more valid test performed
		       nepi++;
		      
		       // Count as a good result

		       summary_good[e1]++;
		       if (sA[e2]) summary_good[e2]++;
		      
		       // Do we want to record this as part of the summary for the first set? 	  
		       if ( z >= par::epi_alpha2)
			 {
			   // first variable will always be in A set
			   summary_sig[e1]++;

			   // but the second may also be in A set
			   if (sA[e2]) summary_sig[e2]++;
			 }
		       
		       // Is this result the best scrore yet for marker in set A?

		       if (z > best_score[e1]) 
			 {
			   best_score[e1] = z;
			   best_partner[e1] = e2;
			 }
		       
		       // The second marker might also be in set A 

		       if (sA[e2])
			 {
			   if (z > best_score[e2])
			     {
			       best_score[e2] = z;
			       best_partner[e2] = e1;
			     }
			 }
		     }
		   
		   // Is this result worth displaying?
		   
		   if (  z >= par::epi_alpha1 )
		     {
		       EPI << setw(4) << locus[e1]->chr << " " 
			   << setw(par::pp_maxsnp) << locus[e1]->name << " "
			   << setw(4) << locus[e2]->chr << " "
			   << setw(par::pp_maxsnp) << locus[e2]->name  << " ";
		       if (lm->isValid())
			 {
			   if ( par::bt) 
			     EPI << setw(12) << exp(b[3])  << " "
				 << setw(12) << chisq  << " "
				 << setw(12) << pvalue << " "
				 << "\n";
			   else
			     EPI << setw(12) << b[3]  << " "
				 << setw(12) << chisq  << " "
				 << setw(12) << pvalue << " "
				 << "\n";
			 }
		       else
			 EPI << setw(12) << "NA"  << " "
			     << setw(12) << "NA"  << " "
			     << setw(12) << "NA" << " "
			     << "\n";

		       EPI.flush();
					  
		     }

		   // Clean up
		   delete lm;

		 }
	       

	    } // Next pair of SNPs
	}
      
    }  	  
  
  if (!par::silent)
    cout << "\n";  
  
  EPI.close();


  //////////////////////
  // Summary of results

  // Skip this for now
  
  if (true)
    {
      f += ".summary";
      EPI.open(f.c_str(),ios::out);
      EPI.clear();
      
      printLOG("Performed a total of "+int2str(nepi)+" valid SNPxSNP tests\n");
      
      printLOG("Writing epistasis summary results to [ " + f + " ] \n");
      
      EPI.precision(4);
      EPI << setw(4) << "CHR" << " "
	  << setw(par::pp_maxsnp) << "SNP" << " "
	  << setw(12) << "N_SIG" << " "
	  << setw(12) << "N_TOT" << " "
	  << setw(12) << "PROP" << " "
	  << setw(12) << "BEST_CHISQ"  << " "
	  << setw(4) << "BEST_CHR" << " " 
	  << setw(par::pp_maxsnp) << "BEST_SNP" << " "
	  << "\n";
      
      int c=0;
      for (int e1=0;e1<nl_all;e1++)
	{
	  if (sA[e1]) 
	    {
	      EPI << setw(4) << locus[e1]->chr << " "
		  << setw(par::pp_maxsnp) << locus[e1]->name << " "
		  << setw(12) << summary_sig[e1] << " "
		  << setw(12) << summary_good[e1] << " "
		  << setw(12) << (double)summary_sig[e1] / (double)summary_good[e1] << " "
		  << setw(12) << best_score[e1] * best_score[e1] << " "
		  << setw(4) << locus[best_partner[e1]]->chr  << " "
		  << setw(par::pp_maxsnp) << locus[best_partner[e1]]->name  << " "
		  << "\n";
	      
	    }
	}
      
      EPI.close();
    }


}
