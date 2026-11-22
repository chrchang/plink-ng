

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
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "model.h"
#include "sets.h"

extern ofstream LOG;
extern Plink * PP;

using namespace std;

void scoreRanges(int, 
		 vector<int> &,
		 map<int, set<Range*> > &,
		 map<Range*,int> &,
		 ofstream &);


void Plink::scoreIndividuals()
{

  map<string,int> mlocus;
  for(int l=0; l<nl_all;l++)
    mlocus.insert(make_pair(locus[l]->name,l));      
  
  string suffix = "";

  map<int,double> qscore;
  vector<double2> qthresh;
  vector<string> qlabel;

  if ( par::score_risk_on_qrange )
    {
      
      checkFileExists( par::score_qfile );
      checkFileExists( par::score_qrange_file );
      
      printLOG("Reading quantitative scores from [ " + par::score_qfile + " ]\n");
      printLOG("Reading score ranges from [ " + par::score_qrange_file + " ]\n");
      
      ifstream Q1( par::score_qfile.c_str() , ios::in );
      while ( ! Q1.eof() )
	{
	  string snp;
	  string str_score;
	  double score;	  
	  Q1 >> snp >> str_score;	  
	  if ( ! from_string<double>( score , str_score , std::dec ) )
	    continue;
	  if ( snp == "" )
	    continue;
	  
	  map<string,int>::iterator i1 = mlocus.find( snp );
	  if ( i1 != mlocus.end() )
	    qscore.insert( make_pair( i1->second , score ) );
	}
      Q1.close();
      printLOG("Read q-scores for " + int2str( qscore.size() ) + " SNPs\n");      
      
      Q1.open( par::score_qrange_file.c_str() , ios::in );
      while ( ! Q1.eof() )
	{
	  // Expect: name, lower, upper
	  string label;
	  double lower, upper;
	  Q1 >> label >> lower >> upper;
	  if ( label == "" )
	    continue;
	  double2 d2(lower,upper);
	  qthresh.push_back( d2 );
	  qlabel.push_back( label );
	}
      Q1.close();
      printLOG("Read " + int2str( qthresh.size() ) + " thresholds to apply\n");            

    }


  printLOG("Reading set of predictors from [ " + par::score_risk_file + " ]\n");
  checkFileExists(par::score_risk_file);

  
  // Main loop in which we consider either multiple takes from the score file, or we 
  // just run through once

  int qcnt = 0;
  
  while(1) 
    {

      string suffix = "";
      double2 th;

      if ( par::score_risk_on_qrange )
	{
	  suffix = "." + qlabel[qcnt];
	  th = qthresh[ qcnt ];
	  printLOG("Thresholding group " 
		   + qlabel[qcnt] 
		   + " on ( " 
		   + dbl2str(th.p1) 
		   + " -- " + dbl2str(th.p2) 
		   + " )\n");
	}

      ifstream PROFIN;
      PROFIN.open( par::score_risk_file.c_str(), ios::in );
      
      string problems;
      
      map<int,double> scores;
      map<int,bool> allele1;
      
      int cnt1 = 0, cnt2 = 0, cnt2b = 0, cnt3 = 0;
      
      while ( ! PROFIN.eof() )
	{
	  
	  // Format assumed: SNP allele score
	  
	  string snp, allele, sscore;
	  double score;
	  
	  PROFIN >> snp >> allele >> sscore;
	  
	  if ( sscore=="" )
	    continue;
	  

	  if ( ! from_string<double>( score, sscore , std::dec))
	    {
	      problems += "BADVAL\t" + snp + "\n";
	      continue;
	    }

	  ++cnt1;
	  
	  map<string,int>::iterator ilocus = mlocus.find(snp);
	  
	  // SNP not found
	  if ( ilocus == mlocus.end() )
	    {
	      problems += "NOSNP\t" + snp + "\n";
	      continue;
	    }


	  int l = ilocus->second;

	  ++cnt2;
	  
	  // Purposely not include this SNP base on Q-range?
	  // Are we only looking at subsets of SNPs?
	  if ( par::score_risk_on_qrange )
	    {
	      map<int,double>::iterator i1 = qscore.find( l ); 
	      if ( i1 == qscore.end() ) 
		continue;
	      
	      double sc = i1->second;

	      if ( sc < th.p1 || sc > th.p2 )
		continue;	      
	    }
	  
	  ++cnt2b;
	  
	  // Allele found?
	  if ( allele == locus[l]->allele1 )
	    {
	      scores.insert(make_pair(l,score));
	      allele1.insert(make_pair(l,false));
	      ++cnt3;
	    }
	  else if ( allele == locus[l]->allele2 )
	    {
	      scores.insert(make_pair(l,score));
	      allele1.insert(make_pair(l,true));
	      ++cnt3;
	    }
	  else
	    problems += "NOALLELE\t" 
	      + snp + " " + allele 
	      + " vs " + locus[l]->allele1 
	      + " " + locus[l]->allele2 + "\n";
	  
	}
      
      PROFIN.close();
      
      if ( par::score_risk_on_qrange )
	printLOG("Read " + int2str(cnt1) 
		 + " predictors; " + int2str(cnt2)      
		 + " mapped to SNPs; " + int2str(cnt2b)
		 + " selected; " + int2str(cnt3) 
		 + " to alleles\n");
      else
	printLOG("Read " + int2str(cnt1) 
		 + " predictors; " + int2str(cnt2)      
		 + " mapped to SNPs; " + int2str(cnt3) 
		 + " to alleles\n");
      
      if ( problems != "" )
	{ 
	  printLOG("Writing problem SNPs in predictor to [ " 
		   + par::output_file_name + suffix + ".nopred ]\n");
	  ofstream O1;
	  O1.open( (par::output_file_name + suffix + ".nopred").c_str() , ios::out );
	  O1 << problems ;
	  O1.close();
	  problems = "";
	}
      
      
      

  ////////////////////////////////
  // Calculate for each individual
  
  printLOG("Writing profiles to [ " + par::output_file_name + suffix + ".profile ]\n");


  ////////////////////////////////////////////
  // First, perform this for all SNPs

  vector_t profile;

  matrix_t set_profile;
  if ( par::profile_sets )
    pS->initialiseSetMapping();

  vector<int> cnt;
  vector<int> acount;
  
  calculateProfile(scores,allele1,profile,set_profile,cnt,acount);



  ///////////////////////////////
  // Report for all individuals

  ofstream PROFOUT;
  string f = par::output_file_name + suffix + ".profile";
  PROFOUT.open( f.c_str(), ios::out );

  PROFOUT << setw(par::pp_maxfid) << "FID" << " " 
	  << setw(par::pp_maxiid) << "IID" << " "
	  << setw(6) << "PHENO" << " " 
	  << setw(6) << "CNT" << " "
          << setw(6) << "CNT2" << " "
	  << setw(8) << "SCORE" << "\n";
  
  for ( int i=0; i<n; i++ )
    {
      Individual * person = sample[i];

      PROFOUT << setw(par::pp_maxfid) << person->fid << " " 
	      << setw(par::pp_maxiid) << person->iid << " "
	      << setw(6) << person->phenotype << " " 
	      << setw(6) << cnt[i] << " "
	      << setw(6) << acount[i] << " " 
	      << setw(8) << profile[i] << "\n";
    }
  
  PROFOUT.close();



  
  ///////////////////////////////////////////////
  // Test association with score and phenotype
  
  if ( par::score_test )
    {
      

      vector_t results;

      ///////////////////////////////////////////////////
      // Set up association model

      bool OLD_assoc_glm_without_main_snp = par::assoc_glm_without_main_snp;
      bool OLD_clist = par::clist;


      par::assoc_glm_without_main_snp = true;
      par::clist = true;
      par::clist_number = 0;


      int totalTerms = 1;

      par::clist_number += totalTerms;
      for (int i=0; i<n; i++)
	sample[i]->clist.resize( par::clist_number );
      
      // Fill in label forms
      int testTerm = par::clist_number - 1;

      clistname.resize( par::clist_number );
      clistname[ testTerm ] = "SCORE";

      ////////////////////////////////////////
      // Put relevant variables in clist slots
      
      for ( int i = 0; i < n; i++)
	sample[i]->clist[ testTerm ] = profile[i];
    
      
      ///////////////////////////////////////////
      // Perform association
      
      glmAssoc(false,*pperm);
            
      
      ///////////////////////////////////////////
      // Report results
      
      ofstream OUTF;
      
      bool valid = model->isValid();

      vector_t b = model->getCoefs();
      vector_t chisq(1,model->getStatistic());
      vector_t pval = model->getPVals();
  
      // NOTE: b includes intercept; pval doesn't
  
      double statistic = valid ? model->getStatistic() : 0;
      double pvalue = pval[ pval.size()-1 ];
      double beta = b[ b.size()-1 ];

      delete model;
  

      string f = par::output_file_name + suffix + ".profile.test";
      PROFOUT.open( f.c_str(), ios::out );

      PROFOUT << setw(12) << "TEST" << " "
	      << setw(12) << "BETA" << " "
	      << setw(12) << "P" << "\n";
      
      for (int c = 0; c < par::clist_number; c++)
	{
	  PROFOUT << setw(12) << clistname[c] << " "
		  << setw(12) << b[c+1] << " "
		  << setw(12) << pval[c] << "\n";
	}
      
      PROFOUT.close();

    }


  ////////////////////////////////////////////////////////
  // 
  // Recalculate scores for specified sets only
  

  if ( par::profile_sets )
    {
      ofstream PROFOUT;
      string f = par::output_file_name + suffix + ".profile.sets";
      PROFOUT.open( f.c_str(), ios::out );
            
      PROFOUT << setw(par::pp_maxfid) << "FID" << " " 
	      << setw(par::pp_maxiid) << "IID" << " "
	      << setw(6) << "PHENO" << " ";
      for (int s=0; s<snpset.size(); s++)
	PROFOUT << setname[s] << "\t";
      PROFOUT << "\n";
      
      for ( int i=0; i<n; i++ )
	{
	  Individual * person = sample[i];
	  
	  PROFOUT << setw(par::pp_maxfid) << person->fid << " " 
		  << setw(par::pp_maxiid) << person->iid << " "
		  << setw(6) << person->phenotype << " ";
	  // Consider each pathway
	  for (int s=0; s<snpset.size(); s++)
	    PROFOUT << set_profile[i][s] << "\t";
	  PROFOUT << "\n";
	}
      PROFOUT.close();
    }



  ////////////////////////////////////////////////////////
  // 
  // Wrap-up, or recalculate with a different set of SNPs?
  

  if ( ! par::score_risk_on_qrange )
    break;
  
  ++qcnt;

  if ( qcnt == qlabel.size() )
    break;
    }

}


void Plink::calculateProfile(map<int,double> & scores, 
			     map<int,bool> & allele1, 
			     vector_t & profile,
			     matrix_t & set_profile,
			     vector<int> & count,
			     vector<int> & acount )
{  
  

  // Generate a vector of scores, one for each individual, given then
  // scoring set (and allele direction) for a set of SNPs

  profile.resize(n,0);
  if ( par::profile_sets )
    sizeMatrix( set_profile , n , snpset.size() );
  count.resize(n,0);
  acount.resize(n,0);


  ///////////////////////////////////////
  // Do we want to score for genes also?

  map<string, set<Range> > ranges;  
  vector<string> rangeLabels;
  map<int,set<Range*> > snp2range;
  map<Range*,int> rangeCount;


  //////////////////////////////////////////
  // Consider each individual and calculate
  // the score
  
  for (int i=0; i<n; i++)
    {

      Individual * person = sample[i];
      
      map<int,double>::iterator i1 = scores.begin();
      map<int,bool>::iterator i2 = allele1.begin();
      
      double score = 0;
      vector_t set_score( snpset.size() , 0 );
      int cnt = 0;
      int cntActual = 0;
      int cntNamedAllele = 0;
	  
      vector<int> flaggedSNPs;
      
      while ( i1 != scores.end() )
	{
	  
	  int l = i1->first;
	  bool a1 = i2->second;

	  bool s1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
	  bool s2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];

	  bool missingGenotype = false;

	  double thisScore = 0;

	  /////////////////////////////////////////////
	  // Individual is missing this genotype
	  
	  // We with either skip, or impute mean

	  if ( s1 && ! s2 ) 
	    {
	      
	      if ( ! par::score_impute_expected ) 
		{
		  ++i1;
		  ++i2;
		  continue;
		}
	      
	      missingGenotype = true;
	      
	      if ( i2->second ) 
		thisScore = ( 1 - locus[l]->freq ) * i1->second;
	      else
		thisScore = locus[l]->freq * i1->second;
	      
	      if ( par::chr_haploid[ locus[l]->chr ] || 
		   ( par::chr_sex[ locus[l]->chr ] && person->sex ) )
		++cnt;
	      else
		{
		  cnt += 2;
		  thisScore *= 2;
		}
	    }


	  // Currently, just an allelic scoring: we could extend this 
	  // to genotypes, dominant/recessive models,
	  
	  bool sawNamedAllele = false;
	  
	  if ( ! missingGenotype ) 
	    { 
	      
	      if ( par::chr_haploid[ locus[l]->chr ] || 
		   ( par::chr_sex[ locus[l]->chr ] && person->sex ) )
		{
		  // A single copy
		  
		  if ( i2->second ) 
		    {
		      if ( s1 ) 
			{
			  thisScore = i1->second;
			  sawNamedAllele = true;
			  ++cntNamedAllele;
			}
		    }
		  else
		    {
		      if ( !s1 ) 
			{
			  thisScore = i1->second;
			  sawNamedAllele = true;
			  ++cntNamedAllele;
			}
		    }
		  
		  ++cnt;	      
		  ++cntActual;
		}
	      else // .. autosomal
		{
		  if ( i2->second ) 
		    {
		      if ( s1 ) 
			{
			  thisScore = i1->second;
			  sawNamedAllele = true;
			  ++cntNamedAllele;
			}
		      if ( s2 ) 
			{
			  thisScore += i1->second;
			  sawNamedAllele = true;
			  ++cntNamedAllele;
			}
		    }
		  else
		    {
		      if ( !s1 ) 
			{
			  thisScore = i1->second;
			  sawNamedAllele = true;
			  ++cntNamedAllele;
			}
		      if ( !s2 ) 
			{
			  thisScore += i1->second;	      
			  sawNamedAllele = true;
			  ++cntNamedAllele;
			}
		    }
		  cnt += 2;	      
		  cntActual +=2;
		}
	    }
	  

	  //////////////////////////////////////////
	  // Accumulate score
	  
	  score += thisScore;
	  
	  

	  //////////////////////////////////////////
	  // Score in pathways also?

	  if ( par::profile_sets )
	    {
	      map<int, set<int> >::iterator si = pS->setMapping.find(l);
	      set<int>::iterator i2 = si->second.begin();
	      while ( i2 != si->second.end() )
		{
		  set_score[ *i2 ] += thisScore;
		  ++i2;
		}
	    }
	  
	  ++i1;
	  ++i2;

	}
      
      
      // Get average per seen loci
      
      if ( cnt>0 ) 
	score /= (double)cnt;

      if ( par::profile_sets )
	{
	  for (int j=0; j<snpset.size(); j++)
	    set_score[j] /= (double)cnt;
	  set_profile[i] = set_score;	  
	}



      // Save for this individual (actual number)

      
      profile[i] = score;
      
      count[i] = cntActual;
      acount[i] = cntNamedAllele;

     
    } // Next individual
   
  return;
}


void scoreRanges(int i, 
		 vector<int> & f,		 
		 map<int, set<Range*> > & snp2range,
		 map<Range*,int> & rangeCount,
		 ofstream & ROUT)
{
  Individual * person = PP->sample[i];
  
  ROUT << setw(par::pp_maxfid) << person->fid << " "
       << setw(par::pp_maxiid) << person->iid << " "
       << setw(6) << person->phenotype << " ";
  
  // Get list of ranges that will be flagged

  set<Range*> mappedRanges;

  for (int l=0; l<f.size(); l++)
    {
      map<int, set<Range*> >::iterator ri = snp2range.find( f[l] );
      
      if ( ri == snp2range.end() )
	continue;

      set<Range*>::iterator si = ri->second.begin();
      while ( si != ri->second.end() )
	{
	  mappedRanges.insert( *si );
	  ++si;
	}
    }
  
  // We now have populated the set mappedRanges
  
  map<Range*,int>::iterator r = rangeCount.begin();

  while ( r != rangeCount.end() )
    {
      if ( mappedRanges.find(r->first) != mappedRanges.end() )
	{
	  ROUT << "1 ";
	  ++rangeCount[ r->first ];
	}
      else
	ROUT << "0 ";
      
      ++r;
      
      ROUT << "\n";
      
    }
}



// OLD VERSION WITH SCORE RANGES CODE IN PLACE

// void Plink::calculateProfile(map<int,double> & scores, 
// 			     map<int,bool> & allele1, 
// 			     vector_t & profile,
// 			     matrix_t & set_profile,
// 			     vector<int> & count,
// 			     vector<int> & acount )
// {  
  

//   // Generate a vector of scores, one for each individual, given then
//   // scoring set (and allele direction) for a set of SNPs

//   profile.resize(n,0);
//   if ( par::profile_sets )
//     sizeMatrix( set_profile , n , snpset.size() );
//   count.resize(n,0);
//   acount.resize(n,0);


//   ///////////////////////////////////////
//   // Do we want to score for genes also?

//   map<string, set<Range> > ranges;  
//   vector<string> rangeLabels;
//   map<int,set<Range*> > snp2range;
//   map<Range*,int> rangeCount;
//   ofstream ROUT;

//   if ( par::score_risk_ranges )
//     {
//       printLOG("Writing range scores to [ " 
// 	       + par::output_file_name 
// 	       + ".profile.ranges ]\n");

//       // Helper function to map ranges to SNPs
//       mapRangesToSNPs( par::score_risk_ranges_file,
// 		       ranges, 
// 		       snp2range );
      
//       map<string, set<Range> >::iterator r = ranges.begin();

//       while ( r != ranges.end() )
// 	{
// 	  set<Range> * theseRanges = &( r->second );
	  
// 	  set<Range>::iterator thisRangeSet = theseRanges->begin();
// 	  while ( thisRangeSet != theseRanges->end() )
// 	    {
// 	      rangeLabels.push_back( thisRangeSet->name );
// 	      ++thisRangeSet;
// 	    }	  
// 	  ++r;
// 	}
            
      
//       ROUT.open( ( par::output_file_name+".profile.ranges").c_str(), ios::out);
   
//       // Header
      
//       ROUT << setw(par::pp_maxfid) << "FID" << " " 
// 	   << setw(par::pp_maxiid) << "IID" << " "
// 	   << setw(6) << "PHENO" << " ";

//       for (int r=0; r<rangeLabels.size(); r++)
// 	ROUT << rangeLabels[r] << " ";
//       ROUT << "\n";
      	      	   
//     }
  

//   //////////////////////////////////////////
//   // Consider each individual and calculate
//   // the score
  
//   for (int i=0; i<n; i++)
//     {

//       Individual * person = sample[i];
      
//       map<int,double>::iterator i1 = scores.begin();
//       map<int,bool>::iterator i2 = allele1.begin();
      
//       double score = 0;
//       vector_t set_score( snpset.size() , 0 );
//       int cnt = 0;
//       int cntActual = 0;
//       int cntNamedAllele = 0;
	  
//       vector<int> flaggedSNPs;
      
//       while ( i1 != scores.end() )
// 	{
	  
// 	  int l = i1->first;
// 	  bool a1 = i2->second;

// 	  bool s1 = par::SNP_major ? SNP[l]->one[i] : person->one[l];
// 	  bool s2 = par::SNP_major ? SNP[l]->two[i] : person->two[l];

// 	  bool missingGenotype = false;

// 	  double thisScore = 0;

// 	  /////////////////////////////////////////////
// 	  // Individual is missing this genotype
	  
// 	  // We with either skip, or impute mean

// 	  if ( s1 && ! s2 ) 
// 	    {
	      
// 	      if ( ! par::score_impute_expected ) 
// 		{
// 		  ++i1;
// 		  ++i2;
// 		  continue;
// 		}
	      
// 	      missingGenotype = true;
	      
// 	      if ( i2->second ) 
// 		thisScore = ( 1 - locus[l]->freq ) * i1->second;
// 	      else
// 		thisScore = locus[l]->freq * i1->second;
	      
// 	      if ( par::chr_haploid[ locus[l]->chr ] || 
// 		   ( par::chr_sex[ locus[l]->chr ] && person->sex ) )
// 		++cnt;
// 	      else
// 		{
// 		  cnt += 2;
// 		  thisScore *= 2;
// 		}
// 	    }


// 	  // Currently, just an allelic scoring: we could extend this 
// 	  // to genotypes, dominant/recessive models,
	  
// 	  bool sawNamedAllele = false;
	  
// 	  if ( ! missingGenotype ) 
// 	    { 
	      
// 	      if ( par::chr_haploid[ locus[l]->chr ] || 
// 		   ( par::chr_sex[ locus[l]->chr ] && person->sex ) )
// 		{
// 		  // A single copy
		  
// 		  if ( i2->second ) 
// 		    {
// 		      if ( s1 ) 
// 			{
// 			  thisScore = i1->second;
// 			  sawNamedAllele = true;
// 			  ++cntNamedAllele;
// 			}
// 		    }
// 		  else
// 		    {
// 		      if ( !s1 ) 
// 			{
// 			  thisScore = i1->second;
// 			  sawNamedAllele = true;
// 			  ++cntNamedAllele;
// 			}
// 		    }
		  
// 		  ++cnt;	      
// 		  ++cntActual;
// 		}
// 	      else // .. autosomal
// 		{
// 		  if ( i2->second ) 
// 		    {
// 		      if ( s1 ) 
// 			{
// 			  thisScore = i1->second;
// 			  sawNamedAllele = true;
// 			  ++cntNamedAllele;
// 			}
// 		      if ( s2 ) 
// 			{
// 			  thisScore += i1->second;
// 			  sawNamedAllele = true;
// 			  ++cntNamedAllele;
// 			}
// 		    }
// 		  else
// 		    {
// 		      if ( !s1 ) 
// 			{
// 			  thisScore = i1->second;
// 			  sawNamedAllele = true;
// 			  ++cntNamedAllele;
// 			}
// 		      if ( !s2 ) 
// 			{
// 			  thisScore += i1->second;	      
// 			  sawNamedAllele = true;
// 			  ++cntNamedAllele;
// 			}
// 		    }
// 		  cnt += 2;	      
// 		  cntActual +=2;
// 		}
// 	    }
	  

// 	  //////////////////////////////////////////
// 	  // Accumulate score
	  
// 	  score += thisScore;
	  

// 	  //////////////////////////////////////////
// 	  // Pathway-specific scores?

// 	  if ( par::profile_sets )
// 	    {
// 	      for ( int j = 0 ; j < snpset.size(); j++ )
// 		{
// 		  // Is this SNP in this pathway? If so, 
// 		  // add to pathway-specific score

// 		  if ( 1 ) 
// 		    set_score[j] += thisScore;
// 		}
// 	    }


// 	  //////////////////////////////////////////
// 	  // Do we want to score a "yes/no" for a
// 	  // gene?

// 	  if ( par::score_risk_ranges && ! missingGenotype )
// 	    {

// 	      // Did we see at least one risk-increasing allele?
	      
// 	      if ( ( sawNamedAllele && i1->second > 0 ) ||
// 		   ( (!sawNamedAllele) && i1->second < 0 ) ) 
// 		{
// 		  flaggedSNPs.push_back(l);
// 		}
// 	    }
	  
// 	  ++i1;
// 	  ++i2;

// 	}
      
      
//       // Get average per seen loci
      
//       if ( cnt>0 ) 
// 	score /= (double)cnt;

//       // Save for this individual (actual number)

//       profile[i] = score;
//       set_profile[i] = set_score;
//       count[i] = cntActual;
//       acount[i] = cntNamedAllele;

//       // Score for ranges?
      
//       if ( par::score_risk_ranges )
// 	{
// 	  scoreRanges(i,flaggedSNPs,snp2range,rangeCount,ROUT);
// 	}

     
//     } // Next individual
  

//   if ( par::score_risk_ranges )
//     {
//       ROUT.close();
      
//       // Also ouput how many times each range was seen
      
//       printLOG("Writing range summary counts to [ "
// 	       + par::output_file_name 
// 	       + ".profile.ranges.summary ]\n");
      
//       ofstream ROUT2;
//       ROUT2.open( ( par::output_file_name+".profile.ranges.summary").c_str(), ios::out);
   
//       // Header
      
//       ROUT2 << setw(18) << "RANGE" << " " 
// 	    << setw(8) << "CNT" << "\n";
      
//       map<Range*,int>::iterator r = rangeCount.begin();
//       while ( r != rangeCount.end() )
// 	{
// 	  ROUT2 << setw(18) << r->first->name << " "
// 		<< setw(8) << r->second << "\n";
// 	  ++r;
// 	}
      
//       ROUT2.close();
//     }
  
 
//   return;
//}
