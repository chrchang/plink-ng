

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

#include "linear.h"
#include "logistic.h"
#include "helper.h"
#include "plink.h"
#include "options.h"
#include "crandom.h"
#include "sets.h"
#include "perm.h"
#include "phase.h"
#include "whap.h"
#include "stats.h"

// Usage of Model

// Fit model                       LinearModel lm;
// Give pointer to PLINK           lm.setPlink(this);
// Set missing data                lm.setMissing();
// Set dependent (adds intercept)  lm.setDependent(Y);
// Addive effects, labels          lm.addAdditiveSNP(CSNP*); 
//                                 lm.label.push_back("ADD");
//	                           lm.addDominanceSNP(CSNP*);
// Sex effect?                     lm.addSexEffect();
//                                 lm.label.push_back("SEX");
// Covariates?                     lm.addCovariate(int);
//              	           lm.label.push_back("COV"+int2str(c+1));
// Interactions?       	           lm.addInteraction(int,int);
// Build design matrix             lm.buildDesignMatrix();
// Fit logistic model              lm.fitLM();
// Fit okay?                       lm.validParameters();
// Display results                 lm.displayResults();
// Get statistic                   lm.getStatistic();


vector_t Plink::glmAssoc(bool print_results, Perm & perm)
{
  
  // The model.cpp functions require a SNP-major structure, if SNP
  // data are being used.  There are some exceptions to this however, 
  // listed below

  if ( par::SNP_major && 
       ! ( par::epi_genebased 
	   || par::set_score 
	   || par::set_step 
	   || par::proxy_glm 
	   || par::dosage_assoc
	   || par::cnv_enrichment_test
	   || par::cnv_glm 
	   || par::score_test 
	   || par::rare_test 
	   || par::gvar ) )
    SNP2Ind();
  


  // Test all SNPs 1 at a time automatically, or is this 
  // a tailored single test?

  int ntests = par::assoc_glm_without_main_snp ? 1 : nl_all;

  vector<double> results(ntests);

  if ( print_results && par::qt && par::multtest )
    tcnt.resize(ntests);

  ofstream ASC;
  if (print_results)  
    {
      string f = par::output_file_name;
      if ( par::bt)
	{
	  f += ".assoc.logistic";
	  printLOG("Writing logistic model association results to [ " 
		   + f + " ] \n");
	}
      else
	{
	  f += ".assoc.linear";
	  printLOG("Writing linear model association results to [ " 
		   + f + " ] \n");
	}

      ASC.open(f.c_str(),ios::out);
      ASC << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " " 
	  << setw(10) << "BP" << " "
	  << setw(4) << "A1" << " "
	  << setw(10) << "TEST" << " "
	  << setw(8) << "NMISS" << " ";
      if ( par::bt && ! par::return_beta )
	ASC << setw(10) << "OR" << " ";
      else
	ASC << setw(10) << "BETA" << " ";
      
      if (par::display_ci)
	ASC << setw(8) << "SE" << " "
	    << setw(8) << string("L"+dbl2str(par::ci_level*100)) << " "
	    << setw(8) << string("U"+dbl2str(par::ci_level*100)) << " ";
      
      ASC << setw(12) << "STAT" << " " 
	  << setw(12) << "P" << " " 	  
	  << "\n";
      ASC.precision(4);
    }
  
  
  /////////////////////////////
  // Determine sex distribution

  int nmales = 0, nfemales = 0;
  for (int i=0; i<n; i++)
    if ( ! sample[i]->missing )
      {
	if ( sample[i]->sex )
	  nmales++;
	else
	  nfemales++;
      }
  
  bool variationInSex = nmales > 0 && nfemales > 0;

  
  //////////////////////////////////////////
  // Iterate over each locus, or just once

  for (int l=0; l<ntests; l++)
    {	

      // Skip possibly (in all-locus mode)
      
      if ( par::adaptive_perm && 
	   ( ! par::assoc_glm_without_main_snp ) && 
	   ( ! perm.snp_test[l]) )
	continue;
      

      //////////////////////////////////////////////////////////
      // X-chromosome, haploid?
      // xchr_model 0: skip non-autosomal SNPs
      
      bool X=false;
      bool automaticSex=false;

      if ( ! par::assoc_glm_without_main_snp )
	{
	  if ( par::xchr_model == 0 )
	    {
	      if ( par::chr_sex[locus[l]->chr] ||
		   par::chr_haploid[locus[l]->chr] )
		continue;
	    }
	  else 
	    if (par::chr_sex[locus[l]->chr]) 
	      X=true;
	}
      

      //////////////////////////////////////////////////////////
      // A new GLM
      
      Model * lm;

       
      //////////////////////////////////////////////////////////
      // Linear or logistic?

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
      
  
      //////////////////////////////////////////////////////////
      // A temporary fix

      if ( par::dosage_assoc || 
	   par::cnv_enrichment_test || 
	   par::cnv_glm ||
	   par::score_test || 
	   par::set_score || 
	   par::proxy_glm || 
	   par::gvar || 
	   par::rare_test ) 
 	lm->hasSNPs(false);


      //////////////////////////////////////////////////////////
      // Set missing data

      lm->setMissing();


      //////////////////////////////////////////////////////////
      // Set genetic model

      if ( par::glm_dominant )
	lm->setDominant();
      else if ( par::glm_recessive || par::twoDFmodel_hethom )
	lm->setRecessive();
      
      string mainEffect = "";
      
      bool genotypic = false;

      /////////////////////////////////////////////////
      // Main SNP
      
      if ( ! par::assoc_glm_without_main_snp ) 
	{
	  
	  genotypic = par::chr_haploid[locus[l]->chr] ? false : par::twoDFmodel ;
	  
	  // Models
	  //               AA    AB    BB
	  // Additive       0     1     2

	  // Dominant       0     1     1
	  // Recessive      0     0     1

	  // Genotypic(1)
	  // Additive       0     1     2
	  // Dom Dev.       0     1     0

	  // Genotypic(2)
	  // Homozygote     0     0     1
	  // Heterozygote   0     1     0

	  
	  ////////////////////////////////////////////////////////////
	  // An additive effect? (or single coded effect) of main SNP

	  if ( par::glm_recessive )
	    mainEffect = "REC";
	  else if ( par::glm_dominant ) 
	    mainEffect = "DOM";
	  else if ( par::twoDFmodel_hethom )
	    mainEffect = "HOM";
	  else
	    mainEffect = "ADD";
	  
	  lm->addAdditiveSNP(l); 
	  lm->label.push_back(mainEffect);
	  
	  
 	  //////////////////////////////////////////////////////////
 	  // Or a 2-df additive + dominance model?
	  
 	  if ( genotypic ) 
 	    {
 	      lm->addDominanceSNP(l);
	      
 	      if ( par::twoDFmodel_hethom )
 		lm->label.push_back("HET");
 	      else
 		lm->label.push_back("DOMDEV");
 	    }
	  
	}


      
      //////////////////////////////////////////////////////////
      // Haplotypes: WHAP test (grouped?)
      
      if ( par::chap_test )
	{
	  
	  // Use whap->group (a list of sets) to specify these, from
	  // the current model (either alternate or null)

	  // Start from second category (i.e. first is reference)
	  for (int h=1; h < whap->current->group.size(); h++)
	    {
	      lm->addHaplotypeDosage( whap->current->group[h] );	    
	      lm->label.push_back( "WHAP"+int2str(h+1) );
	    }
	}


      //////////////////////////////////////////////////////////
      // Haplotypes: proxy test
      
      if ( par::proxy_glm )
	{
	  
	  // Unlike WHAP tests, we now will only ever have two
	  // categories; and a single tested coefficient

	  set<int> t1 = haplo->makeSetFromMap(haplo->testSet);
	  lm->addHaplotypeDosage( t1 );
	  lm->label.push_back( "PROXY" );
	    
	}

      if ( par::test_hap_GLM )
	{
	  // Assume model specified in haplotype sets
	  // Either 1 versus all others, or H-1 versus  
	  // terms for omnibus

	  set<int>::iterator i = haplo->sets.begin();
	  while ( i != haplo->sets.end() )
	    {
	      set<int> t;
	      t.insert(*i);
	      lm->addHaplotypeDosage( t );
	      lm->label.push_back( haplo->haplotypeName( *i ) );
	      ++i;
	    }
	}


      
      //////////////////////////////////////////////////////////
      // Conditioning SNPs?
      // (might be X or autosomal, dealth with automatically)
      
      if (par::conditioning_snps)
	{
	  if ( par::chap_test ) 
	    {
	      for (int c=0; c<conditioner.size(); c++)
		{
		  if ( whap->current->masked_conditioning_snps[c] )
		    {
		      lm->addAdditiveSNP(conditioner[c]); 
		      lm->label.push_back(locus[conditioner[c]]->name);
		    }
		}
	    }
	  else
	    {
	      for (int c=0; c<conditioner.size(); c++)
		{
		  lm->addAdditiveSNP(conditioner[c]); 
		  lm->label.push_back(locus[conditioner[c]]->name);
		}
	    }
	}
      


      //////////////////////////////////////////////////////////      
      // Sex-covariate (necessary for X chromosome models, unless
      // explicitly told otherwise)
      
      if ( ( par::glm_sex_effect || ( X && !par::glm_no_auto_sex_effect ) )
	   && variationInSex )
	{
	  automaticSex = true;
	  lm->addSexEffect();
	  lm->label.push_back("SEX");	  
	}
      
  

      //////////////////////////////////////////////////////////
      // Covariates?

      if (par::clist)
	{
	  for (int c=0; c<par::clist_number; c++)
	    {
	      lm->addCovariate(c);
	      lm->label.push_back(clistname[c]);
	    }
	}


      //////////////////////////////////////////////////////////
      // Interactions

      // addInteraction() takes parameter numbers
      // i.e. not covariate codes
      
      // 0 intercept
      // 1 {A}
      //   {D}
      //   {conditioning SNPs}
      //   {sex efffect}
      //   {covariates}

      // Allow for interactions between conditioning SNPs, sex, covariates, etc
      	
      
      ////////////////////////////////////////
      // Basic SNP x covariate interaction? 
      
      // Currently -- do not allow interactions if no main effect 
      // SNP -- i.e. we need a recoding of things here.

      if ( par::simple_interaction && ! par::assoc_glm_without_main_snp )
	{
	  
	  // A, D and haplotypes by conditioning SNPs, sex, covariates
	  
	  int cindex = 2;
	  if ( genotypic )
	    cindex = 3;
	    	  
	  for (int c=0; c<conditioner.size(); c++)
	    {
	      lm->addInteraction(1,cindex);
	      lm->label.push_back(mainEffect+"xCSNP"+int2str(c+1));	  
	      
	      if ( genotypic )
		{
		  lm->addInteraction(2,cindex);
		  if ( par::twoDFmodel_hethom )
		    lm->label.push_back("HETxCSNP"+int2str(c+1));	  
		  else
		    lm->label.push_back("DOMDEVxCSNP"+int2str(c+1));	  
		}
	      
	      cindex++;
	    }

	  if ( automaticSex )
	    {
	      lm->addInteraction(1,cindex);
	      lm->label.push_back(mainEffect+"xSEX");	  
	      
	      if ( genotypic )
		{
		  
		  lm->addInteraction(2,cindex);
		  if ( par::twoDFmodel_hethom )
		    lm->label.push_back("HETxSEX");
		  else
		    lm->label.push_back("DOMDEVxSEX");
		}
	      
	      cindex++;
	    }
	  for (int c=0; c<par::clist_number; c++)
	    {
	      lm->addInteraction(1,cindex);
	      lm->label.push_back(mainEffect+"x"+clistname[c]);	  
	      
	      if ( genotypic )
		{
		  lm->addInteraction(2,cindex);
		  
		  if ( par::twoDFmodel_hethom )		  
		    lm->label.push_back("HETx"+clistname[c]); 
		  else
		    lm->label.push_back("DOMDEVx"+clistname[c]); 
		}
	      
	      cindex++;
	      
	    }
	}
      


      
      //////////////////////////////
      // Fancy X chromosome models
      
      if ( X && automaticSex && par::xchr_model > 2 )
	{
	  
	  // Interaction between allelic term and sex (i.e. 
	  // allow scale of male effect to vary)
	  
	  int sindex = 2;
	  if ( genotypic )
	    sindex++;
	  sindex += conditioner.size();
	  
	  lm->addInteraction(2,sindex);
	  lm->label.push_back("XxSEX");	  
	  
	  // xchr model 3 : test ADD + XxSEX
	  // xchr model 4 : test ADD + DOM + XxSEX
	}

      

      //////////////////////////////
      // Build design matrix

      lm->buildDesignMatrix();
      

      //////////////////////////////
      // Clusters specified?

      if ( par::include_cluster ) 
	{
	  lm->setCluster();
	}
      
      //////////////////////////////////////////////////
      // Fit linear or logistic model (Newton-Raphson)
      
      lm->fitLM();


      ////////////////////////////////////////
      // Check for multi-collinearity

      lm->validParameters();
      

      ////////////////////////////////////////
      // Obtain estimates and statistic

      if (print_results)
	lm->displayResults(ASC,locus[l]);
      


      ////////////////////////////////////////////////
      // Test linear hypothesis (multiple parameters)

      // Perform if:
      //   automatic 2df genotypic test  ( --genotypic )
      //     OR
      //   sex-tests            ( --xchr-model )
      //     OR
      //   test of everything   ( --test-all )
      //     OR
      //   user has specified user-defined test  ( --tests )
      
      if ( ( genotypic && ! par::glm_user_parameters ) 
	   || par::glm_user_test 
	   || par::test_full_model )
	{

	  vector_t h; // dim = number of fixes (to =0)
	  matrix_t H; // row = number of fixes; cols = np
	  int df;
	  string testname;
	  

	  ////////////////////////////////////////////////
	  // Joint test of all parameters

	  if (par::test_full_model) 
	    {
	      df = lm->getNP() - 1;
	      h.resize(df,0);
	      testname = "FULL_"+int2str(df)+"DF";
	      sizeMatrix(H,df,lm->getNP());
	      for (int i=0; i<df; i++)
		H[i][i+1] = 1;
	    }
	  
	  ////////////////////////////////////////////////
	  // Joint test of user-specified parameters
	  
	  else if (par::glm_user_test) 
	    {
	      df = par::test_list.size();
	      h.resize(df,0);
	      testname = "USER_"+int2str(df)+"DF";
	      sizeMatrix(H,df,lm->getNP());
	      for (int i=0; i<df; i++)
		if ( par::test_list[i]<lm->getNP() )
		  H[i][par::test_list[i]] = 1;
	      
	    }
	  
	  ////////////////////////////////////////////////
	  // Joint test of additive and dominant models
	  
	  else if ( genotypic )
	    {
	      testname = "GENO_2DF";
	      df = 2;
	      h.resize(2,0);
	      sizeMatrix(H,2,lm->getNP());
	      H[0][1] = H[1][2] = 1; 	  
	    }
	  
	  else if ( X && par::xchr_model == 3 )	  
	    {
	      testname = "XMOD_2DF";
	    }


	  ////////////////////////////////////////////////
	  // Joint test of all parameters

	  double chisq = lm->isValid() ? lm->linearHypothesis(H,h) : 0;
	  double pvalue = chiprobP(chisq,df);
	  
	  // If filtering p-values
	  if ( (!par::pfilter) || pvalue <= par::pfvalue ) 
	    {	 
	      
	      ASC << setw(4) << locus[l]->chr << " " 
		  << setw(par::pp_maxsnp) << locus[l]->name << " " 
		  << setw(10) << locus[l]->bp << " "		
		  << setw(4) << locus[l]->allele1 << " " 
		  << setw(10) << testname << " "
		  << setw(8) << lm->Ysize() << " " 
		  << setw(10) << "NA" << " ";
	      
	      if (par::display_ci)
		ASC << setw(8) << "NA" << " " 
		    << setw(8) << "NA" << " "
		    << setw(8) << "NA" << " ";
	      
	      if (lm->isValid() && realnum(chisq) )
		ASC << setw(12) << chisq << " " 
		    << setw(12) << pvalue << "\n"; 
	      else
		ASC << setw(12) << "NA" << " " 
		    << setw(12) << "NA" << "\n"; 
	    }


 	}
    
      
      ////////////////////////////////////////
      // Store statistic (1 df chisq), and p-value
      // if need be ( based on value of testParameter )

      if ( ! par::assoc_glm_without_main_snp )
	results[l] = lm->getStatistic();


      if ( par::qt && print_results && par::multtest )
	tcnt[l] = lm->Ysize() - lm->getNP();
      


      //////////////////////////////////////////////
      // Clear up linear model, if no longer needed
      
      if ( par::chap_test || 
	   par::test_hap_GLM || 
	   par::set_step || 
	   par::set_score || 
	   par::proxy_glm || 
	   par::dosage_assoc || 
	   par::cnv_enrichment_test ||	   
	   par::cnv_glm ||
	   par::score_test ||
	   par::gvar || 	   
	   par::rare_test ) 
	{
	  // Responsibility to clear up in parent routine
	  model = lm; 
	}
      else
	{
	  delete lm;	  
	}

      // Flush output buffer
      ASC.flush();


      // Next SNP
    }

  
 
  
  if (print_results)
    ASC.close();
  
  return results;

}

