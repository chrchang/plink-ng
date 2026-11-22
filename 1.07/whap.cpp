

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
//
// Conditional Haplotype tests (CHAP models)
//
//   A null and alternate model specified in terms of 
//           haplotypes (potentially grouped)
//           covariates
//           conditioning SNPs
//
//   Only focus on allelic, autosomal main effects right now
//
//
//   SNPs done with --condition {list}
//   Covariates with --covar / --covar-name / --covar-number
//   Haplotypes with --hap-snps
//
//   Covariates are always present under both alternate and null
//   SNPs can be dropped from the alternate (tested)
//   Haplotypes can be dropped (tested) and grouped


// A helper function
void displayHaploGroups(ofstream &, ChapModel &, HaploPhase *);

string ci(double coef, double se)
{
  string c = " (";
  if ( (!par::bt) || par::return_beta )
    {
      c += dbl2str( coef - par::ci_zt * se , 3 ) + "; " ;
      c += dbl2str( coef + par::ci_zt * se , 3 ) + " )" ;
    }
  else
    {
      c += dbl2str( exp( coef - par::ci_zt * se ), 3) + "; " ;
      c += dbl2str( exp( coef + par::ci_zt * se ), 3) + " )" ;
    }
  return c;
}


vector_t Plink::conditionalHaplotypeTest(bool print_results, Perm & perm)
{


  ///////////////////////////////////////////////
  //                                           //
  // Some basic setup first                    //
  //                                           //
  ///////////////////////////////////////////////

  Chap thisCModel(this, haplo);
  
  whap = & thisCModel;

  ChapModel alternateModel;
  ChapModel nullModel;

  
  // Use basic GLM function to fit linear and logistic 
  // models: although, let it know that there will not 
  // be a 'main' SNP 

  par::assoc_glm_without_main_snp = true;


  // Return a single result

  vector_t results(1);


  printLOG("Writing conditional haplotype tests to [ " 
	   + par::output_file_name + ".chap ]\n");

  ofstream CH;
  CH.open((par::output_file_name+".chap").c_str(),ios::out);
  CH.precision(3);
  
  CH << "+++ PLINK conditional haplotype test results +++ \n\n";
  

  ///////////////////////////////////////////////
  //                                           //
  // Phase haplotypes                          //
  //                                           //
  ///////////////////////////////////////////////

  haplo->phaseAllHaplotypes(false,*pperm);

  // Record the number of common haplotypes
  
  int nch = 0;  
  for (int h=0; h < haplo->nh; h++)
    if ( haplo->f[h] >= par::min_hf )
      ++nch;
    
  
  CH << haplo->ns << " SNPs, and " 
     << nch << " common haplotypes ( MHF >= " << par::min_hf << " ) "
     << "from " << haplo->nh << " possible\n\n";

  CH << setw(4) << "CHR" << " "
     << setw(12) << "BP" << " "
     << setw(12) << "SNP" << " "
     << setw(4) << "A1" << " "
     << setw(4) << "A2" << " "
     << setw(10) << "F" << "\n";

  for (int s=0; s< haplo->ns; s++)
    CH << setw(4) << locus[haplo->S[s]]->chr << " "
       << setw(12) << locus[haplo->S[s]]->bp << " "
       << setw(12) << locus[haplo->S[s]]->name << " "
       << setw(4) << locus[haplo->S[s]]->allele1 << " "
       << setw(4) << locus[haplo->S[s]]->allele2 << " "
       << setw(10) << locus[haplo->S[s]]->freq << "\n";
  CH << "\n";

  if ( nch == 0 ) 
    {
      results[0] = 0;
      CH << "Exiting... no common haplotypes to test\n";
      CH.close();
      return results;
    }
  
  
  ///////////////////////////////////////////////
  //                                           //
  // Build models                              //
  //                                           //
  ///////////////////////////////////////////////
 
  whap->setModels(alternateModel, nullModel);
  
  whap->build(alternateModel);

  whap->build(nullModel);



  ///////////////////////////////////////////////
  //                                           //
  // Display models                            //
  //                                           //
  ///////////////////////////////////////////////

  CH << "Haplogrouping: each {set} allowed a unique effect\n";
  CH << "Alternate model\n";
  displayHaploGroups(CH, alternateModel,haplo);

  CH << "Null model\n";
  displayHaploGroups(CH, nullModel,haplo);
  CH << "\n";
  
  if ( ! whap->isNested() )
    error("The null model is not nested in the alternate: please respecify");


  
  ///////////////////////////////////////////////
  //                                           //
  // Fit alternate model                       //
  //                                           //
  ///////////////////////////////////////////////

  whap->current = & alternateModel;
  
  glmAssoc(false,perm);
  
  Model * alternate = model;



  ///////////////////////////////////////////////
  //                                           //
  // Fit null model                            //
  //                                           //
  ///////////////////////////////////////////////

  // Ensure we have the same individuals
  vector<bool> missingInAlternate = alternate->getMissing();
  for (int i=0; i<n; i++)
    sample[i]->missing2 = missingInAlternate[i];

  whap->current = & nullModel;  
  
  glmAssoc(false,perm);
  
  Model * null = model;
  

  ///////////////////////////////////////////////
  //                                           //
  // Fit group-specific models                 //
  //                                           //
  ///////////////////////////////////////////////
  
  // If there is more than 1 group, these models test each group
  // against all others, in both the alternate and null

  vector_t alt_specific_pval;
  bool hasAltSpecifics = alternateModel.group.size() > 2 && par::chap_add_grp_specifics;
  
  if ( hasAltSpecifics ) 
    {

      ////////////////////////////////////////////
      // Create a model with no haplotype effects
  
      ChapModel simpleNullModel = nullModel;
      
      simpleNullModel.group.clear();

      set<int> t;	
      for (int h=0; h< haplo->nh; h++)
	if ( haplo->f[h] >= par::min_hf )
	  t.insert(h);
      simpleNullModel.group.push_back(t);      

      whap->current = & simpleNullModel;  
      glmAssoc(false,perm);
      Model * simplenull = model;

      
      
      ////////////////////////////////////////////
      // Run all haplotype-specific models

      

      for (int h=0; h< haplo->nh; h++)
	{
	  
	  if ( haplo->f[h] < par::min_hf ) 
	    continue;

	  ChapModel simpleAlternateModel = nullModel;
	  
	  simpleAlternateModel.group.clear();
	  simpleAlternateModel.group.resize(1);
	  
	  for (int h2=0; h2< haplo->nh; h2++)
	    {
	      
	      if ( haplo->f[h2] < par::min_hf ) 
		continue;
	      
	      if ( h2 == h )
		{
		  set<int> t;
		  t.insert(h);
		  simpleAlternateModel.group.push_back(t);
		}
	      else
		simpleAlternateModel.group[0].insert(h);	  

	    }
	  
      
	  ////////////////////////////////////////////
	  // Perform test, and model comparison
	  
	  whap->current = & simpleAlternateModel;  
	  glmAssoc(false,perm);
	  Model * simplealternate = model;
	  
	  
	  //////////////////////////
	  // Store test statistic  
	  
	  alt_specific_pval.push_back( modelComparisonPValue(simplealternate, 
							     simplenull) );
	  
	  
	  // Next haplo-group
	}
      
    }

  ///////////////////////////////////////////////
  //                                           //
  // Fit sub-null model                        //
  //                                           //
  ///////////////////////////////////////////////
  
  // If the null model contains >1 group, then 
  // perform the alternate:null comparisons 
  // separately for each sub-group, if the null
  // still contains fewer parameters than the 
  // alternate
  
  vector_t subnull_pval;

  int subnullModels = 0;

  // Worth doing? 
  if ( nullModel.group.size() > 1 && 
       alternateModel.group.size() - nullModel.group.size() > 1 )
    {

      for (int g = 0; g < nullModel.group.size(); g++)
	{
	  // For this null-group, do the haplotypes in the 
	  // alternate belong to >1 group? If so, perform a 
	  // separate test
	  
	  set<int>::iterator ih = nullModel.group[g].begin();
	  set<int> aGroup;
	  
	  while ( ih != nullModel.group[g].end() )
	    {
	      aGroup.insert( alternateModel.haploGroup.find(*ih)->second );	      	      
	      ++ih;
	      
	    }

	  if ( aGroup.size() > 1 ) 
	    {
	      
	      ++subnullModels;
	      
	      // Re-jig a new model, that is basically like the 
	      // alternate, except for this one group.
	      
	      ChapModel subnullModel = alternateModel;
	      
	      // Edit group only (haploGroup will not be used, 
	      // so ignore that for now...)
	      
	      set<int>::iterator ai = aGroup.begin();
	      

	      // Get first group: arbitrarily choose this group to
	      // be the merged-to group

	      int ng = *ai;
	      
	      // All other groups, merge with this first one
	      
	      ++ai;

	      // Iterate over each alternate-model haplogroup to be merged
	      while ( ai != aGroup.end() )
		{

		  // Iterate over all haplotypes in this haplogroup
		  set<int>::iterator s = alternateModel.group[ *ai ].begin();
		  while ( s != alternateModel.group[ *ai ].end() )
		    {
		      subnullModel.group[ ng ].insert( *s ); 
		      ++s;		      
		    }

		  subnullModel.group[ *ai ].clear();
		  ++ai;

		}
	     
	      // Now erase empty groups
	      for (int g=0; g<subnullModel.group.size(); g++)
		{
		  if ( subnullModel.group[g].size() == 0 )
		    {
		      subnullModel.group.erase( subnullModel.group.begin() + g );
		      --g;
		    }
		}

	      whap->current = & subnullModel;  
	      glmAssoc(false,perm);
	      Model * subnull = model;

	      // Store test statistic
	      
	      subnull_pval.push_back( modelComparisonPValue(alternate, subnull) );
	      
	    }
	  else
	    {
	      subnull_pval.push_back(-1);
	    }
	}          

    }



  ///////////////////////////////////////////////
  //                                           //
  // Report model comparisons                  //
  //                                           //
  ///////////////////////////////////////////////


  // Check that both models converged

  if ( ! ( alternate->isValid() && null->isValid() ) )
    error("Could not fit conditional haplotype models:\n   "
	  "collinearity issues from a badly-specified model\n");
      

  vector_t coeff1 = alternate->getCoefs();
  vector_t se1 = alternate->getSE();
  
  vector_t coeff0 = null->getCoefs();
  vector_t se0 = null->getSE();
  
  
  vector<string> label1 = alternate->label;
  vector<string> label0 = null->label;
   
 
  //////////////////////////////////
  // Convert to odds ratios?
  
  vector_t odds1 = coeff1;
  vector_t odds0 = coeff0;
  
  if ( par::bt )
    {
      for (int h=0; h<coeff1.size(); h++)
	odds1[h] = exp(coeff1[h]);
      for (int h=0; h<coeff0.size(); h++)
	odds0[h] = exp(coeff0[h]);
    }


  ////////////////////////////////////////////////////
  // List by groups; first null then alternate groups 

  string clabel_alternate = par::bt ? "OR(A)" : "BETA(A)";
  string clabel_null = par::bt ? "OR(N)" : "BETA(N)";

  string cunder_alternate = par::bt ? "-------" : "---------";
  string cunder_null = par::bt ? "-------" : "---------";
  
  int estimate_size = 12;
  if (par::display_ci) 
    {
      estimate_size += 16;
      cunder_alternate += "----------------";
      cunder_null += "----------------";
    }

  CH << setw( haplo->ns+5 ) << "HAPLO" << "   "
     << setw(10) << "FREQ" << " ";

  CH << setw( estimate_size ) << clabel_alternate << " ";

  if ( hasAltSpecifics )
    CH << setw(12) << "SPEC(A)" << " ";
  
  CH << setw( estimate_size ) << clabel_null << " ";


  if ( subnullModels > 1 )
    CH << setw(12) << "SUBNULL P" << " "; 
  CH << "\n";

  CH << setw( haplo->ns+5 ) << "-------" << "   "
     << setw(10) << "------" << " ";

  CH << setw( estimate_size ) << cunder_alternate << " ";

  if ( hasAltSpecifics )
    CH << setw(12) << "---------" << " ";
  
  CH << setw( estimate_size ) << cunder_null << " ";

  
  if ( subnullModels > 1 )
    CH << setw(12) << "-----------" << " "; 
  CH << "\n";
  
  for ( int g=0; g<nullModel.group.size(); g++ )
    {
      bool printed0 = false;
      for ( int g2=0; g2<alternateModel.group.size(); g2++ )
	{
	  set<int>::iterator ih = alternateModel.group[g2].begin();
	  
	  bool printed = false;
	  while ( ih != alternateModel.group[g2].end() )
	    {
	      
	      // Only display if this haplotype is in the
	      // null group
	      
	      if ( nullModel.group[g].find( *ih ) 
		   == nullModel.group[g].end() )
		{
		  ++ih;
		  continue;
		}
	      

	      CH << setw( haplo->ns+5 ) << haplo->haplotypeName( *ih ) << "   "
		 << setw(10) << haplo->f[ *ih ] << " ";
	      
	      if ( ! printed ) 
		{
		  
		  if ( g2==0 )
		    CH << setw( estimate_size ) << "(-ref-)" << " ";
		  else
		    {
		      int p = alternateModel.haploGroup.find(*ih)->second;
		      if ( realnum( odds1[p] ) ) 
			{
			  string r = dbl2str( odds1[p] , 4 );
			  if ( par::display_ci ) 
			    r += ci( coeff1[p] , se1[ p ]  );
			  CH << setw( estimate_size ) << r << " ";		  
			  
			}
		      else
			{
			  CH << setw( estimate_size ) << "NA" << " ";		  
			}
		    }

		  if ( hasAltSpecifics )
		    CH << setw(12) << alt_specific_pval[ g2 ] << " ";
		  
		  printed = true;
		}
	      else 
		{ 
		  CH << setw( estimate_size ) << "|   " << " ";
		  
		  if ( hasAltSpecifics )
		    CH << setw(12) << " " << " ";
		  
		}
	      
	      ///////////////////////////////
	      // Display corresponding null 
	      
	      if ( ! printed0 ) 
		{		  
		  if ( g==0 )
		    CH << setw( estimate_size ) << "(-ref-)" << " ";
		  else
		    {
		      int p = nullModel.haploGroup.find(*ih)->second;
		      if ( realnum( odds0[p] ) ) 
			{
			  string r = dbl2str( odds0[p] , 4 );
			  if ( par::display_ci ) 
			    r += ci( coeff0[p] , se0[ p ]  );
			  CH << setw( estimate_size ) << r << " ";		  			  
			}
		      else
			CH << setw( estimate_size ) << "NA" << " ";		  
		    }
		}
	      else 
		CH << setw( estimate_size ) << "|   " << " ";
	      
	      ////////////////////////////////
	      // Display corresponding subnull?
	      
	      if ( subnullModels > 1  ) 
		{
		  if ( ! printed0 )
		    {
		      if ( subnull_pval[g] < 0 ) 
			CH << setw(12) << "n/a" << " ";
		      else
			CH << setw(12) << subnull_pval[g] << " ";		  
		    }
		}
	      
	      printed0 = true;
	      
	      
	      CH << "\n";	  
	      ++ih;
	    }
	  
	  // Delimiter alternate groups, unless we 
	  // are also about to delimit a null group
	  
	}
      
      if ( g < nullModel.group.size() - 1 ) 
	CH << "\n";
      
    }
  
  CH << setw( haplo->ns+5 ) << "-------" << "   "
     << setw(10) << "------" << " ";

  
  CH << setw( estimate_size ) << cunder_alternate << " ";
  if ( hasAltSpecifics )
    CH << setw(12) << "---------" << " ";      
  CH << setw( estimate_size ) << cunder_null << " ";
    
  if ( subnullModels > 1 )
    CH << setw(12) << "-----------" << " "; 
  CH << "\n";
  

  /////////////////////////////////////////////////
  //  Display other covariates, conditioning SNPs
  
  // 0=intercept; 1 -> (H-1) haplotype-group coefficients; 
  // then conditioning SNPs; then other covariates

  int p1 = alternateModel.group.size();
  int p0 = nullModel.group.size();
  
  // Only an intercept, then need to add 1
  if ( p1 == 0 ) p1++;
  if ( p0 == 0 ) p0++;


  /////////////////////////////////////////////////////////
  // Conditioning SNPs: will always feature in alternate; 
  // may or may not feature in null
  
  if ( conditioner.size() > 0 )
    {
      CH << "\n";
      
      CH << setw( haplo->ns+5) << "SNPS" << " " 
	 << setw(12) << " " << " "
	 << setw( estimate_size ) << clabel_alternate << " "
	 << setw( estimate_size ) << clabel_null << "\n";
      
      CH << setw( haplo->ns+5) << "-----" << " " 
	 << setw(12) << " " << " "
	 << setw( estimate_size ) << cunder_alternate << " "
	 << setw( estimate_size ) << cunder_null << "\n";

      for (int s=0; s < conditioner.size(); s++)
	{

	  CH << setw( haplo->ns+5 ) << label1[p1] << " " 
	     << setw(12) << " " << " ";
	  if ( realnum( coeff1[ p1 ] ) )
	    {
	      string r = dbl2str( coeff1[ p1 ] , 4 );
	      if ( par::display_ci ) 
		r += ci( coeff1[ p1 ] , se1[ p1 ]  );
	      CH << setw( estimate_size ) << r << " ";		  
			  
	    }
	  else
	    CH << setw( estimate_size ) << "NA" << " ";

	  p1++;

	  if ( nullModel.masked_conditioning_snps[s] ) 
	    {
	      if ( realnum( coeff0[ p0 ] ) )
		{
		  string r = dbl2str( coeff0[ p0 ] , 4 );
		  if ( par::display_ci ) 
		    r += ci( coeff0[ p0 ] , se0[ p0 ]  );
		  CH << setw( estimate_size ) << r << " ";
		}
	      else 
		CH << setw( estimate_size ) << "NA" << "\n";
	      p0++;
	    }
	  else
	    CH << setw( estimate_size ) << "(dropped)" << "\n";
	  
	}

    }



  ///////////////////////////////////////////////////////
  // Other covariates: these will always feature in alternate 
  // and null

  if ( par::clist && par::clist_number > 0 ) 
    {

      CH << "\n";

      CH << setw( haplo->ns+5) << "COVAR" << " " 
	 << setw(12) << " " << " "
	 << setw( estimate_size ) << clabel_alternate << " "
	 << setw( estimate_size ) << clabel_null << "\n";

      CH << setw( haplo->ns+5) << "-----" << " " 
	 << setw(12) << " " << " "
	 << setw( estimate_size ) << cunder_alternate << " "
	 << setw( estimate_size ) << cunder_null << "\n";

      for (int s=0; s<par::clist_number; s++)
	{
	  CH << setw( haplo->ns+5) << label1[p1] << " " 
	     << setw(12) << " " << " ";

	  if ( ! par::display_ci )
	    {
	      CH << setw( estimate_size ) << coeff1[ p1 ] << " "
		 << setw( estimate_size ) << coeff0[ p0 ] << "\n";
	    }
	  else
	    {
	      string r = dbl2str( coeff1[ p1 ] , 4 ) + ci( coeff1[ p1 ] , se1[ p1 ]  );
	      CH << setw( estimate_size ) << r << " ";		  	      
	      r = dbl2str( coeff0[ p0 ] , 4 ) + ci( coeff0[ p0 ] , se0[ p0 ]  );
	      CH << setw( estimate_size ) << r << " ";		  	      
	    }
	  
	  // Advance to next coefficient, 
	  ++p1;
	  ++p0;

	}
    }

  if ( alternate->isSexInModel() )
    {

      // Have we already displayed a header for covariates?
      if ( ! ( par::clist && par::clist_number > 0 ) )
	{
	  CH << "\n";

	  CH << setw( haplo->ns+5) << "COVAR" << " " 
	     << setw(12) << " " << " "
	     << setw( estimate_size ) << clabel_alternate << " "
	     << setw( estimate_size ) << clabel_null << "\n";
	  
	  CH << setw( haplo->ns+5) << "-----" << " " 
	     << setw(12) << " " << " "
	     << setw(estimate_size ) << cunder_alternate << " "
	     << setw(estimate_size ) << cunder_null << "\n";
	}

      CH << setw( haplo->ns+5 ) << label1[p1] << " " 
	 << setw(12) << " " << " ";

      if ( ! par::display_ci )
	{
	  CH << setw(estimate_size) << coeff1[ p1 ] << " "
	     << setw(estimate_size) << coeff0[ p0 ] << "\n";
	}
      else
	{
	  string r = dbl2str( coeff1[ p1 ] , 4 ) + ci( coeff1[ p1 ] , se1[ p1 ]  );
	  CH << setw( estimate_size ) << r << " ";		  	      
	  r = dbl2str( coeff0[ p0 ] , 4 ) + ci( coeff0[ p0 ] , se0[ p0 ]  );
	  CH << setw( estimate_size ) << r << " ";		  	      
	}
      
      ++p1;
      ++p0;
    }
  
   if ( p1 != coeff1.size() || p0 != coeff0.size() )
     error("Internal error in whap.cpp -- p1,p0 do not align");
  
  /////////////////////////////////////////////////
  //  Display overall model comparison statistics


  CH << "\n"
       << "Model comparison test statistics:\n\n";
  
  CH << setw(25) << " " << " "
     << setw(10) << "Alternate" << " " 
     << setw(10) << "Null" << "\n";
  
  if ( par::bt )
    {
      LogisticModel * lalternate = (LogisticModel*)alternate;
      LogisticModel * lnull = (LogisticModel*)null;
      
      CH << setw(25) << "-2LL : " << " " 
	 << setw(10) << lalternate->getLnLk() << " "
	 << setw(10) << lnull->getLnLk() << " "
	 << "\n\n";
      
      if ( lalternate->getNP() - lnull->getNP()  == 0 ) 
	CH << setw(25) << "Likelihood ratio test: "
	   << " ( not a valid comparison: identical models, df = 0 )\n";	   
      else
	{
	  double lrt = lnull->getLnLk() - lalternate->getLnLk();
	  if ( lrt < 0 || !realnum(lrt) ) lrt = 0;
	  int df = lalternate->getNP() - lnull->getNP();
	  double pval = chiprobP( lrt, df);
	  
	  CH << setw(25) << "Likelihood ratio test: "
	     << "chi-square = "
	     << lrt
	     << "\n"
	     << setw(25) << " "
	     << "df = " 
	     << df
	     << "\n"
	     << setw(25) << " "
	     << "p = " ;
	  if ( pval < 0 || ! realnum(pval) ) 
	    CH << "NA";
	  else
	    CH << pval;
	  CH << "\n";
	}
      
    }
  else
    {
      // Quantitative traits

      LinearModel * lalternate = (LinearModel*)alternate;
      LinearModel * lnull = (LinearModel*)null;

      CH << setw(25) << "R-squared : " << " " 
	 << setw(10) << lalternate->calculateRSquared() << " "
	 << setw(10) << lnull->calculateRSquared() << " "
	 << "\n";
      CH << setw(25) << "Adjusted R-squared : " << " " 
	 << setw(10) << lalternate->calculateAdjustedRSquared() << " "
	 << setw(10) << lnull->calculateAdjustedRSquared() << " "
	 << "\n\n";

//       CH << setw(25) << "Mallow's C : " << " " 
// 	 << lalternate->calculateMallowC(lnull) 
// 	 << "\n";
      
      
      double F = lalternate->calculateFTest(lnull);
      
      if ( F < 0 || !realnum(F) ) F = 0;
      
      double pvalue = pF( F, 
			  alternate->getNP() - null->getNP(),
			  alternate->Ysize() - alternate->getNP() - 1 );
      
      string df = int2str(alternate->getNP() - null->getNP())
	+", "+int2str(alternate->Ysize() - alternate->getNP() - 1);

      CH << setw(25) << "F-statistic comparison : "
	 << "F = "
	 << F
	 << "\n"
	 << setw(25) << " "
	 << "df = " 
	 << df
	 << "\n"
	 << setw(25) << " "
	 << "p = ";
      if ( pvalue < 0 || !realnum(pvalue) ) 
	CH << "NA";
      else
	CH << pvalue;
      CH << "\n";
      
    }

  


  ///////////////////////////////////////////////
  //                                           //
  // We're done                                //
  //                                           //
  ///////////////////////////////////////////////


  CH.close();
  
  delete alternate;
  delete null;
  
  return results;
  
}



void Chap::determineTestType()
{  
  // REDUNDANT 
}


void Chap::build(ChapModel & model)
{
  
  bool isNull = (&model) == null;
  
  // Either:
  //  1) Grouping for alternate and/or null
  //  2) Specific SNPs for alternate and/or null
  //  3) Sole-variant framing
  //  4) Independent effects
  //  5) Haplotype-specific


  model.group.clear();

  bool useDefault = false;
  
  string modelDescription = isNull ? 
    par::chap_model0 : par::chap_model1;
	  

  if ( par::chap_specified_groups )
    {
      
      // Make comma as hash group-delimiter code
      modelDescription = searchAndReplace(modelDescription,","," # ");
      
      // Expand haplotype equality statements
      modelDescription = searchAndReplace(modelDescription,"="," ");
      
      // Tokenize
      vector<string> tok;
      string buf; 
      stringstream ss(modelDescription);
      while (ss >> buf)
	tok.push_back(buf);
      
      
      // Check the same haplotype isn't specified more than once
      set<string> hapin;
      for (int h=0; h<tok.size(); h++)
	{
	  if ( tok[h] == "#" )
	    continue;
	  if ( hapin.find( tok[h] ) != hapin.end() )
	    {
	      if ( ! par::silent )
		cout << "\n";
	      error("Symbol " + tok[h] + " appears more than once in haplotype list");
	    }
	  hapin.insert(tok[h]);
	}

      // We have a list of haplotypes which need 
      // parsing; allow for wildcards
      
      // * = all haplotypes in one group
      // % = all haplotypes in own group
      // # = group delimiter
      
      // Check: cannot have both * and %
      bool groupAll = false;
      bool ungroupAll = false;
      
      for (int i=0; i<tok.size(); i++)
	{
	  if ( tok[i] == "*" )
	    groupAll = true;
	  else if ( tok[i] == "%" )
	    ungroupAll = true;
	}
      
      if ( groupAll && ungroupAll ) 
	error("Cannot specify * and % on same haplotype model");
      
      if ( isNull && ! ungroupAll ) 
	groupAll = true;
      
      if ( (!isNull) && ! groupAll ) 
	ungroupAll = true;
      
      if ( groupAll || ungroupAll ) 
	{
	  
	  set<int> toAdd;
	  
	  for ( int h=0; h < H->nh; h++ )
	    {
	      // Is this haplotype already explicitly listed?
	      
	      string hname = H->haplotypeName( h );
	      bool listed = false;
	      
	      for (int i=0; i< tok.size(); i++)
		{
		  if ( tok[i] == hname ) 
		    {
		      listed = true;
		      break;
		    }
		}
	      
	      if ( ! listed ) toAdd.insert(h);
	      
	    }
	  
	  // Is there anything to add?
	  if ( toAdd.size() > 0 ) 
	    {
	      if ( groupAll ) 
		tok.push_back("#");
	      set<int>::iterator ih = toAdd.begin();
	      while ( ih != toAdd.end() )
		{
		  if ( ungroupAll ) 
		    tok.push_back("#");
		  tok.push_back( H->haplotypeName(*ih) );
		  ++ih;
		}
	    }
	}
      
      
      /////////////////////////////////////////////
      // Now work out the wild-card expanded list
      set<int> t;
      for ( int i = 0 ; i < tok.size() ; i++ ) 
	{
	  if ( tok[i] == "#" ) 
	    {
	      if ( t.size() > 0 ) 
		model.group.push_back( t );
	      t.clear();
	    }
	  else if ( tok[i] == "*" || tok[i] == "%" ) 
	    continue;
	  else
	    {
	      // find haplotype code the long way...
	      for (int h=0; h< H->nh; h++)
		{
		  if ( H->f[h] < par::min_hf ) 
		    continue;
		  
		  if ( tok[i] == H->haplotypeName(h) )
		    t.insert(h);
		}
	    }	      
	}
      
      if ( t.size() > 0 ) 
	model.group.push_back(t);
      
    }
  else if ( par::chap_specified_snps )
    {
      
      // Assume that modelDescription contains a list of SNPs
      
      map<string,int> mapping;
      for (int l=0; l<P->nl_all; l++)
	mapping.insert(make_pair( P->locus[l]->name,l));
      
      NList nl(P->nl_all);
      vector<int> snplist = nl.deparseStringList(modelDescription,&mapping);
      
      if (snplist.size() == 0 )
	useDefault = true;
      
      setSNPList(snplist, model);
      
    }
  else if ( par::chap_sole_variant && isNull )
    {
      
      // This could be a list of SNPs, or a list of haplotypes
      // If SNPs, in NList form (i.e. allowing for 
      // If haplotypes, just in common delimited form
      
      // Under the alternate, we do not do anything here (i.e. thus the
      // condition above, which means the default alternate coding 
      // will be used)
      
      // We may also have specified, with "--control-alleles", that only 
      // a subset of the implied groupings are constrained under the null
      
      modelDescription = par::chap_entity;
      
      // We need to determine whether haplotypes or SNPs specified
      
      map<string,int> mapping;
      for (int s=0; s<H->ns; s++)
	mapping.insert(make_pair( P->locus[H->S[s]]->name,s));
      for (int h=0; h<H->nh; h++)
	mapping.insert(make_pair( H->haplotypeName(h),H->ns + h));
      
      if ( mapping.size() != H->ns + H->nh ) 
	error("Problem, as some SNPs and haplotypes appear not to have unique names");

      NList nl(H->ns + H->nh);
      vector<int> lst = nl.deparseStringList(modelDescription,&mapping);
      
      if (lst.size() == 0 )
	useDefault = true;
      
      bool isSNP = false;
      bool isHAP = false;
      for (int i=0; i< lst.size(); i++)
	if ( lst[i] >= H->ns ) isHAP = true;
	else if ( lst[i] < H->ns ) isSNP = true;
      
      if ( isSNP && isHAP ) 
	error("Cannot specify both SNPs and hapltoypes for --chap-control");

      if ( isHAP && par::chap_sole_variant_specific_alleles ) 
	error("Can only specify --control-alleles when SNPs are listed for --control");
      
      if ( isSNP ) 
	{
	  // Convert to locus 0..nl_all coding
	  for (int i=0; i<lst.size(); i++)
	    lst[i] = H->S[lst[i]]; 
	  
	  setSNPList(lst, model);
	}
      else
	{
	  // For any haplotype found, make so that it has it's 
	  // own group

	  model.group.clear();
	  model.group.resize(1); // main null group

	  for (int h=0; h< H->nh; h++)
	    {
	      if ( H->f[h] < par::min_hf ) 
		continue;
	      
	      bool found = false;
	      
	      for (int i=0; i<lst.size(); i++)
		if ( lst[i] - H->ns == h )
		  {
		    set<int> t;
		    t.insert(h);
		    model.group.push_back(t);
		    found = true;
		  }
	      if ( ! found ) 
		model.group[0].insert(h);
	      
	    }
	}
    }
  else if ( par::chap_independent_effect && isNull )
    {

      // Set SNPs for all *except* the one(s) specified, under the null

      modelDescription = par::chap_entity;

      map<string,int> mapping;
      for (int s=0; s<H->ns; s++)
	mapping.insert(make_pair( P->locus[H->S[s]]->name,s));	  
	

      // Use NList to return negative complement of SNPs listed

      NList nl(H->ns,false);
      vector<int> snplist = nl.deparseStringList(modelDescription,&mapping);
      
      if (snplist.size() == 0 )
	useDefault = true;
      
      // Convert to locus 0..nl_all coding
      for (int i=0; i<snplist.size(); i++)
	snplist[i] = H->S[snplist[i]]; 
      
      setSNPList(snplist, model);
      
    }
  else if ( par::chap_haplotype_specific && ! isNull ) 
    {
      
      // Set a single haplotype under the alternative

      modelDescription = par::chap_entity;

      map<string,int> mapping;
      for (int h=0; h<H->nh; h++)
	mapping.insert(make_pair( H->haplotypeName(h),h));
      
      NList nl(H->nh);
      vector<int> lst = nl.deparseStringList(modelDescription,&mapping);
      
      if (lst.size() == 0 )
	useDefault = true;
      
      model.group.clear();
      model.group.resize(1); // main null group
      
      for (int h=0; h< H->nh; h++)
	{
	  
	  if ( H->f[h] < par::min_hf ) 
	    continue;
	  
	  bool found = false;
	  
	  for (int i=0; i<lst.size(); i++)
	    {
	      if ( lst[i] == h )
		{
		  set<int> t;
		  t.insert(h);
		  model.group.push_back(t);
		  found = true;
		}
	      if ( ! found ) 
		model.group[0].insert(h);
	    }
	  
	}
      
    }
  else
    {
      // If nothing else done by now...
      useDefault = true;
    }
  
  
  /////////////////////////////////////////////////
  // Use defaults? Default null model is all 
  // haplotypes in one group. 
  
  if ( useDefault )
    {
      
      model.group.clear();
      
      if ( isNull )
	{      
	  set<int> t;	
	  for (int h=0; h< H->nh; h++)
	    if ( H->f[h] >= par::min_hf )
	      t.insert(h);
	  model.group.push_back(t);      
	}
      else
	{
	  for (int h=0; h< H->nh; h++)
	    if ( H->f[h] >= par::min_hf )
	      {
		set<int> t;	
		t.insert(h);
		model.group.push_back(t);
	      }      
	}
    }

  

  ////////////////////////////////////////////////////
  // Also make haploGroup's also -- might not need these?

  model.haploGroup.clear();

  for (int g=0; g< model.group.size(); g++)
    {
      set<int>::iterator ih = model.group[g].begin();
      while ( ih != model.group[g].end() )
	{
	  
	  // cout << H->haplotypeName(*ih) << " " << g << "\n";

	  model.haploGroup.insert(make_pair( *ih, g));
	  
	  ++ih;
	}
    }

  
  /////////////////////////////////////////////
  // What about dropping conditioning SNPs?

  // By default, include all conditioning SNPs
  
  model.masked_conditioning_snps.clear();
  model.masked_conditioning_snps.resize(P->conditioner.size(),true);
      
  // *Unless* this is the null, and we have requested that some be
  // dropped

  if ( isNull && par::chap_drop_snps )
    {      
     
      // List of conditioning SNPs
      
      if ( P->conditioner.size() == 0 ) 
	error("No conditioning SNPs in the analysis:\n   cannot test " 
	      + par::chap_drop_snps_list );
      
      map<string,int> mapping;
      for (int c=0; c < P->conditioner.size(); c++)
	mapping.insert(make_pair( P->locus[P->conditioner[c]]->name,c));
      
      // List of SNPs to drop:
      
      NList nl( P->conditioner.size() );
      vector<int> lst = nl.deparseStringList(par::chap_drop_snps_list,&mapping);
      
      for (int l=0; l<lst.size(); l++)
	{
	  model.masked_conditioning_snps[l] = false;
	}
    }


}


void Chap::setModels(ChapModel & a, ChapModel & n)
{
  alternate = &a;
  null = &n;
}

bool Chap::isNested()
{

  // If two haplotypes have the same grouping in the alternate, they
  // must also have the same in the null
  
  for (int h1 = 0; h1 < H->nh-1 ; h1++)
    for (int h2 = h1+1; h2 < H->nh; h2++)
      {
	if ( alternate->haplotypesInSameGroup(h1,h2) )
	  if ( ! null->haplotypesInSameGroup(h1,h2) )
	    return false;
      }
  return true;
}


bool ChapModel::haplotypesInSameGroup(int h1, int h2)
{
  map<int,int>::iterator i1 = haploGroup.find(h1);
  map<int,int>::iterator i2 = haploGroup.find(h2);
  
  if ( i1 == haploGroup.end() ||
       i2 == haploGroup.end() )
    return true;

  if ( i1 == haploGroup.end() ||
       i2 == haploGroup.end() )
    error("Internal problem in ChapModel...");

  return ( i1->second == i2->second );
  
}


void Chap::setSNPList(vector<int> & snplist, ChapModel & model)
{

  boolvec_t snpmask(H->ns,false);
  for (int l=0;l<snplist.size(); l++)
    for ( int s=0; s < H->ns; s++)
      if ( H->S[s] == snplist[l] )
	snpmask[s] = true;		      
  
  map<int,int> subhaplotypes = H->makeSubHaplotypeSet(snpmask);
  
  model.group.clear();
  int cnt=0;
  map<int,int> added;
  for (int h = 0; h < H->nh; h++)
    {	      
      if ( H->f[h] >= par::min_hf )
	{		  
	  int g = subhaplotypes.find(h)->second;
	  map<int,int>::iterator gi = added.find(g);
	  if ( gi == added.end() )
	    {
	      added.insert(make_pair(g,cnt));
	      cnt++;
	      set<int> t;
	      t.insert(h);
	      model.group.push_back(t);
	    }
	  else
	    model.group[ added.find(g)->second ].insert(h);
	}
    }
}


void displayHaploGroups(ofstream & CH, ChapModel & m, HaploPhase * haplo)
{
  bool closeThenEnd = false;
  int cnt = 0;
  CH << "   ";
  for (int h=0; h < m.group.size(); h++)
    {

      CH << "{ ";
      set<int>::iterator ih = m.group[h].begin();

      int cnt2 = m.group[h].size();
      while ( ih != m.group[h].end() )
	{
	  if ( cnt>0 && cnt2 < m.group[h].size() ) 
	    CH << ", ";
	  CH << haplo->haplotypeName( *ih );
  
	  ++cnt;
	  --cnt2;

	  if ( cnt * haplo->ns > 40 ) 
	    {

	      if ( cnt2 == 0) 
		closeThenEnd = true;
	      else
		{
		  CH << "\n     ";
		  cnt = 0;
		}
	    }
  
	  ++ih;
	}
      // Close group      
      CH << " }  ";

      if ( closeThenEnd )
	{
	  CH << "\n   ";
	  closeThenEnd = false;
	  cnt=0;
	}
    }
  CH << "\n";

}



