


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
#include <algorithm>

#include "plink.h"
#include "sets.h"
#include "options.h"
#include "helper.h"
#include "model.h"
#include "stats.h"
#include "phase.h"

extern Plink * PP;

Set::Set(vector<vector<int> > & ss) : snpset(ss)
{
  sizeSets();
}

void Set::sizeSets()
{

  cur.resize(snpset.size());
  for(int s=0;s<snpset.size();s++)
    cur[s].resize(snpset[s].size(),true);

  // Specific to SET-based tests

  if ( (par::assoc_test || par::TDT_test ) && par::set_test
       && !par::hotel)
    {
      s_min.resize(snpset.size());
      s_max.resize(snpset.size());
      stat_set.resize(snpset.size());
      pv_set.resize(snpset.size());
      pv_maxG_set.resize(snpset.size());
      pv_maxE_set.resize(snpset.size());

      for(int i=0;i<snpset.size();i++)
	{

	  // If no constraints given, then the
	  // number of tests == size of set

	  if (par::set_min==-1 ) s_min[i] = 0;
	  else if (par::set_min > snpset[i].size() )
	    s_min[i] = snpset[i].size();
	  else s_min[i] = par::set_min-1;
	  if (par::set_max==-1 || par::set_max > snpset[i].size() )
	    s_max[i] = snpset[i].size();
	  else s_max[i] = par::set_max;
	  if (s_min>s_max) s_min[i]=s_max[i];

	  int s = (s_max[i] - s_min[i]);

	  stat_set[i].resize(s);
	  pv_set[i].resize(s);
	  pv_maxG_set[i].resize(s);
	  pv_maxE_set[i].resize(s);

	  if ( ! par::set_score )
	    {
	      for (int j=0; j<s; j++)
		stat_set[i][j].resize(par::replicates+1);
	    }
	      for (int j=0; j<s; j++)
		pv_set[i][j].resize(par::replicates+1);

	}
    }

}


//////////////////////////////////////////////////////
//                                                  //
//    Remove 0-sized sets                           //
//                                                  //
//////////////////////////////////////////////////////

void Set::pruneSets(Plink & P)
{

  int pruned=0;

  for(int i=0;i<snpset.size();i++)
    {
      if (snpset[i].size() == 0)
	{
	  snpset.erase(snpset.begin()+i);
	  P.setname.erase(P.setname.begin()+i);
	  i--;
	  pruned++;
	}
    }

  P.printLOG(int2str(pruned)+" sets removed (0 valid SNPs)\n");

  // Resize all the other set arrays
  sizeSets();

}


//////////////////////////////////////////////////////
//                                                  //
//    Prune sets based on multi-collinearity        //
//                                                  //
//////////////////////////////////////////////////////

void Set::pruneMC(Plink & P, bool disp,double VIF_threshold)
{

  P.printLOG("Pruning sets based on variance inflation factor\n");

  for (int s=0; s<snpset.size(); s++)
    {

      if (!par::silent)
	cout << s+1 << " of " << snpset.size() << " sets pruned           \r";
      int nss = snpset[s].size();

      vector<double> mean;         // Sample mean
      vector<vector<double> > var; // Covariance matrix

      vector<int> nSNP(0);
      for (int j=0; j<nss; j++)
	{
	  nSNP.push_back( snpset[s][j] );
	}

      // Calculate covariance matrix (full sample)
      // (sizes and populates mean and var)
      // this routine uses the 'flag' variable

      var = calcSetCovarianceMatrix(nSNP);


      // Perform VIF pruning, setting filters (S.cur[][])

      vector<bool> p = vif_prune(var,VIF_threshold,nSNP);

      for (int i=0; i<nss; i++)
	if (!p[i]) cur[s][i]=false;

  }

  if (!par::silent)
    cout << "\n";

  if (disp)
    {
      ofstream SET1, SET2;
      string f = par::output_file_name + ".set.in";

      P.printLOG("Writing pruned-in set file to [ " + f + " ]\n");
      SET1.open(f.c_str(),ios::out);

      f = par::output_file_name + ".set.out";
      P.printLOG("Writing pruned-out set file to [ " + f + " ]\n");
      SET2.open(f.c_str(),ios::out);

      for (int s=0; s<snpset.size(); s++)
	{

	  int nss = snpset[s].size();

	  SET1 << P.setname[s] << "\n";
	  SET2 << P.setname[s] << "\n";

	  for (int j=0; j<nss; j++)
	    {
	      if (cur[s][j])
		SET1 << P.locus[snpset[s][j]]->name << "\n";
	      else
		SET2 << P.locus[snpset[s][j]]->name << "\n";
	    }

	  SET1 << "END\n\n";
	  SET2 << "END\n\n";
	}

      SET1.close();
      SET2.close();
    }
}



//////////////////////////////////////////////////////
//                                                  //
//    Remove SNPs not in any set                    //
//                                                  //
//////////////////////////////////////////////////////

void Set::dropNotSet(Plink & P)
{

  /////////////////////////////////////////////
  // Drop any SNPs that do not belong in a set

  vector<bool> drop(P.nl_all,true);
  for (int i=0;i<snpset.size();i++)
    {
      for (int j=0; j < snpset[i].size(); j++)
	{
	  drop[snpset[i][j]] = false;
	}
    }

  map<int,int> nmap;
  int cnt = 0;
  for (int l=0; l<P.nl_all; l++)
    if ( ! drop[l]  )
      nmap.insert( make_pair( l , cnt++ ) );


  P.deleteSNPs(drop);


  // We now need to update SNP codes

  for (int i=0;i<snpset.size();i++)
    for (int j=0; j < snpset[i].size(); j++)
      {
	int t = snpset[i][j];
	snpset[i][j] = nmap.find(t)->second;
      }

}



//////////////////////////////////////////////////////
//                                                  //
//    Create LD map within each set                 //
//                                                  //
//////////////////////////////////////////////////////

void Set::makeLDSets()
{

  ldSet.clear();
  ldSet.resize( snpset.size() );


  //////////////////////////////////////////////////////
  // If pre-calculated, we can read a .ldset file

  if ( par::set_r2_read )
    {
//       checkFileExists(par::set_r2_read_file);
//       PP->printLOG("Read LD set information from [ " + par::set_r2_read_file + " ]\n");

//       map<string,int> mlocus;
//       makeLocusMap(*PP,mlocus);

//       ifstream SIN;
//       SIN.open( par::set_r2_read_file.c_str() , ios::in );
//       while ( ! SIN.eof() )
// 	{

// 	  vector<string> l = tokenizeLine( SIN );

// 	  if ( SIN.eof() )
// 	    break;

// 	  if ( l.size() < 2 )
// 	    continue;

// 	  // SET ISNP PROXIES...
// 	  int nprox = l.size() - 2;

// 	  // Lookup SNP names
// 	  int isnp =

// 	  for ( int j = 0; j < nprox; j++)
// 	    {
// 	      int l1 = snpset[i][j];
// 	      int l2 = snpset[i][k];

// 	      double rsq = -1;

// 	      if ( par::set_r2_phase )
// 		rsq = PP->haplo->rsq(l1,l2);
// 	      else
// 		rsq = PP->correlation2SNP(l1,l2,true,false);

// 	      if ( rsq >= par::set_r2_val )
// 		{
// 		  ldSet[i][j].insert(k);
// 		  ldSet[i][k].insert(j);

// 	    }

// 	}
//      return;
    }


  //////////////////////////////////////////////////////
  // Otherwise, calculate LD based on raw genotype data

  for (int i=0;i<snpset.size();i++)
    {

      // Is this SNP pair in LD?

      ldSet[i].resize( snpset[i].size() );
      for (int j=0; j < snpset[i].size(); j++)
	ldSet[i][j].clear();

      // Consider each unique pair of SNPs
      // in this set

      for (int j=0; j < snpset[i].size(); j++)
 	{

	  for (int k=j+1; k < snpset[i].size(); k++)
	    {

	      int l1 = snpset[i][j];
	      int l2 = snpset[i][k];

	      double rsq = -1;

	      if ( PP->locus[l1]->chr == PP->locus[l2]->chr )
		{
		  if ( par::set_r2_phase )
		    rsq = PP->haplo->rsq(l1,l2);
		  else
		    rsq = PP->correlation2SNP(l1,l2,true,false);
		}

	      if ( rsq >= par::set_r2_val )
		{
		  ldSet[i][j].insert(k);
		  ldSet[i][k].insert(j);
		}
	    }

	}
    }


  // Output LD sets?
  if ( par::set_r2_write )
    {

      PP->printLOG("Writing LD sets to [ " + par::output_file_name + ".ldset ]\n");
      ofstream SOUT;
      SOUT.open( ( par::output_file_name + ".ldset").c_str() , ios::out);

      for (int i=0;i<snpset.size();i++)
	{

	  for (int j=0; j < snpset[i].size(); j++)
	    {
	      Locus * loc = PP->locus[ snpset[i][j] ];

	      set<int> & lset = ldSet[i][j];

	      if ( lset.size() > 0 )
		{
		  SOUT << PP->setname[i] << " ";
		  SOUT << loc->name << " ";
		  //SOUT << lset.size() << " ";

		  set<int>::iterator k = lset.begin();
		  while ( k != lset.end() )
		    {
		      int l = snpset[i][*k];
		      SOUT << PP->locus[l]->name << " ";
		      ++k;
		    }
		  SOUT << "\n";
		}
	    }
	}
      SOUT.close();
    }



}


//////////////////////////////////////////////////////
//                                                  //
//    Create map of SNP number of set codes         //
//                                                  //
//////////////////////////////////////////////////////

void Set::initialiseSetMapping()
{

  setMapping.clear();

  for (int i=0;i<snpset.size();i++)
    for (int j=0; j < snpset[i].size(); j++)
      {

	int l = snpset[i][j];

	map<int,set<int> >::iterator si = setMapping.find(l);

	// Either we haven't yet seen the SNP...
	if ( si == setMapping.end() )
	  {
	    set<int> t;
	    t.insert(i);
	    setMapping.insert(make_pair(l,t));
	  }
	else
	  {
	    // ... or we have
	    si->second.insert(i);
	  }

	// Next SNP
      }

}


//////////////////////////////////////////////////////
//                                                  //
//    Sum-statistic scoring (original)              //
//                                                  //
//////////////////////////////////////////////////////

void Set::cumulativeSetSum_WITHLABELS(Plink & P, vector<double> & original)
{

//   // Consider each set
//   for (int i=0;i<snpset.size();i++)
//     {

//       vector<SetSortedSNP> t;

//       // Gather set of all chi-sqs (map sorts them automatically)
//       for (int j=0; j < snpset[i].size(); j++)
// 	{
// 	  SetSortedSNP s;
// 	  s.chisq = original[snpset[i][j]];
// 	  s.name = P.locus[snpset[i][j]]->name;
// 	  s.locus = snpset[i][j];
// 	  t.push_back(s);
// 	}


//       // Sort t
//       sort(t.begin(),t.end());

//       // Store results for s_min through s_max
//       double s=0;
//       int j=0;
//       vector<string> t2;

//       for( vector<SetSortedSNP>::reverse_iterator p = t.rbegin(); p!=t.rend(); p++)
// 	{


// // 	  ////////////////////////////////
// // 	  // Using an r-sq threshold also?

// // 	  double max_r2 = 0;

// // 	  if ( par::set_r2 )
// // 	    {
// // 	      int l0 = p->locus;
// // 	      for (int l=0; l< inSet.size(); l++)
// // 		{
// // 		  double r = PP->haplo->rsq( l0, inSet[l] );
// // 		  if ( r > max_r2 )
// // 		    max_r2 = r;
// // 		}
// // 	    }

// 	  ////////////////////////////////
// 	  // Add this SNP to the set?

// // 	  if ( (!par::set_r2) ||
// // 	       max_r2 <= par::set_r2_val )
// // 	    {
// 	      s += p->chisq;
// 	      if (j>=s_min[i] && j<s_max[i])
// 		{
// 		  stat_set[i][j-s_min[i]][0] = s/(double)(j+1);
// 		  t2.push_back(p->name);
// //		  inSet.push_back(p->locus);
// 		}
// 	      j++;
// //	    }

// 	}

//       // And save
//       setsort.push_back(t2);

//    }
}


//////////////////////////////////////////////////////
//                                                  //
//    Sum-statistic scoring (permutation)           //
//                                                  //
//////////////////////////////////////////////////////

void Set::cumulativeSetSum_WITHOUTLABELS(vector<double> & perm, int p)
{

//   vector<double> t;

//   // Consider each set
//   for (int i=0;i<snpset.size();i++)
//     {

//       t.resize(0);

//       // Gather set of chi-sqs
//       for (int j=0;j<snpset[i].size();j++)
// 	t.push_back(perm[snpset[i][j]]);



//       /////////////////////////////
//       // Sort them

//       sort(t.begin(),t.end());



//       /////////////////////////////
//       // Store

// //      vector<int> inSet;

//       double s=0;
//       for (int j=0;j<s_max[i];j++)
// 	{


// 	  ////////////////////////////////
// 	  // Add this SNP to the set?

// // 	  double max_r2 = 0;

// // 	  if ( par::set_r2 )
// // 	    {
// // 	      int l0 = p->locus;
// // 	      for (int l=0; l< inSet.size(); l++)
// // 		{
// // 		  double r = PP->haplo->rsq( l0, inSet[l] );
// // 		  if ( r > max_r2 )
// // 		    max_r2 = r;
// // 		}
// // 	    }


// 	  ////////////////////////////////
// 	  // Add this SNP to the set?

// // 	  if ( (!par::set_r2) ||
// // 	       max_r2 <= par::set_r2_val )
// // 	    {

// 	      s += t[t.size()-1-j];
// 	      if (j>=s_min[i] && j<s_max[i])
// 		{
// 		  stat_set[i][j-s_min[i]][p] = s/(double)(j+1);
// 		}

// //	    }

// 	}


//    }
}


//////////////////////////////////////////////////////
//                                                  //
//    Sum-statistic empircal p-value calculation    //
//                                                  //
//////////////////////////////////////////////////////

void Set::empiricalSetPValues()
{

  int R = par::replicates;


  //////////////////////////////////////////////////
  // Basic p-values, for original and each replicate

  // For the j'th SNP of the i'th SET, calculate how many times
  // the other permutations exceed it (permutations 0 to R, where
  // 0 is the original result)

  for (int p0=0;p0<=R;p0++)   // index
    for (int p1=0;p1<=R;p1++) // all other perms (including self)
      for (int i=0;i<stat_set.size();i++)
	for (int j=0;j<stat_set[i].size();j++)
	  if (stat_set[i][j][p1] >= stat_set[i][j][p0] ) pv_set[i][j][p0]++;

  // Find best p-values per rep (overall, per set)

  for (int p=0;p<=R;p++)
    {

      double maxE_set = 1;
      vector<double> maxG_set(pv_set.size(),1);

      // Consider each score
      for (int i=0;i<pv_set.size();i++)
	for (int j=0;j<pv_set[i].size();j++)
	  {
	    // Make into p-value (will include self: i.e. N+1)
	    pv_set[i][j][p] /= R+1;

	    if (pv_set[i][j][p] < maxG_set[i]) maxG_set[i] = pv_set[i][j][p];
	    if (pv_set[i][j][p] < maxE_set) maxE_set = pv_set[i][j][p];
	  }

      // Score max values
      for (int i=0;i<pv_set.size();i++)
	for (int j=0;j<pv_set[i].size();j++)
	  {
	    if (maxG_set[i] <= pv_set[i][j][0]) pv_maxG_set[i][j]++;
	    if (maxE_set    <= pv_set[i][j][0]) pv_maxE_set[i][j]++;
	  }
    }

}




////////////////////////////////////////////////////////////////////
//                                                                //
//   Score-profile based test                                     //
//                                                                //
////////////////////////////////////////////////////////////////////

void Set::profileTestSNPInformation(int l, double odds)
{

  // If we are passed a SNP here, it is because it significant at the
  // specified par::set_score_p threshold

  // We have to ask: does this SNP belong to one or more sets?  If so,
  // store the SNP number and odds ratio, for each set (i.e.  build up
  // a profile to score; we do not need to save allele, as it is
  // always with reference to the minor one

  map<int, set<int> >::iterator si = setMapping.find(l);

  if ( si == setMapping.end() )
    {
      return;
    }

  set<int>::iterator li = si->second.begin();

  while ( li != si->second.end() )
    {
      profileSNPs[ *li ].push_back( l );
      profileScore[ *li ].push_back( odds );
      ++li;
    }

}




vector_t Set::profileTestScore()
{

  ///////////////////////////////////////////////////
  // For each set, calculate per-individual scores, then
  // regress this on the phenotype, then save a Wald
  // test statistic

  vector_t results;

  for (int i=0; i<snpset.size(); i++)
    {

      vector_t profile;
      vector<int> count;
      vector<int> acount;

      map<int,double> scores;
      map<int,bool> allele1;

      for (int j=0; j<profileSNPs[i].size(); j++)
	{
	  scores.insert(make_pair( profileSNPs[i][j], profileScore[i][j] ));
	  allele1.insert(make_pair( profileSNPs[i][j], false ));
	}

      // Record set size (# significant SNPs; use set_min to store this)

      s_min[i] = profileSNPs[i].size();


      ///////////////////////////////
      // Any significant SNPs?

      if ( scores.size() == 0 )
        {
          // Record a null score
          results.push_back( 0 );
          continue;
        }


      ////////////////////////////////
      // Calculate actual profile

      matrix_t dummy;
      PP->calculateProfile(scores, allele1, profile, dummy , count, acount);


      ///////////////////////////////////////////////
      // Save as the covariate, the mean score (i.e.
      // average by number of seen SNPs)

      for (int k=0; k < PP->n; k++)
	{

	  Individual * person = PP->sample[k];

	  if ( count[k] == 0 || person->flag )
	    person->missing = true;
	  else
	    {
	      person->clist[0] = profile[k] / (double)count[k];
	      person->missing = false;
	    }

	}


      ////////////////////////////////
      // Regress phenotype on profil

      PP->glmAssoc(false,*PP->pperm);


      //////////////////////////////////////////////
      // Reset original missing status

      vector<Individual*>::iterator it = PP->sample.begin();
      while ( it != PP->sample.end() )
	{
	  (*it)->missing = (*it)->flag;
	  ++it;
	}

      ////////////////////////////////////////////////
      // Save test statistic for permutation purposes

      double statistic = PP->model->getStatistic();

      PP->model->validParameters();

      if ( ! PP->model->isValid() )
	statistic = -1;

      results.push_back( statistic );


      ////////////////////////////////////////////////
      // Clear up GLM model

      delete PP->model;

    }


  // Finally, important to clear the profile scores now,
  // so that the next permutation starts from scratch

  profileSNPs.clear();
  profileScore.clear();

  profileSNPs.resize( snpset.size() );
  profileScore.resize( snpset.size() );

  return results;
}


void Set::profileTestInitialise()
{

  PP->printLOG("Initalising profile-based set test\n");

  // Set up the mapping to determine which set(s)
  // a given SNP is in

  initialiseSetMapping();

  // Clear the scores

  profileSNPs.clear();
  profileScore.clear();

  profileSNPs.resize( snpset.size() );
  profileScore.resize( snpset.size() );


  ///////////////////////////////////////////////////
  // Set-up for use of the Linear or Logistic Models

  par::assoc_glm_without_main_snp = true;

  if ( PP->clistname.size() > 0 )
    error("Cannot specify covariates with --set-score");


  //////////////////////////////////////////////
  // Use flag to store original missing status

  vector<Individual*>::iterator i = PP->sample.begin();
  while ( i != PP->sample.end() )
    {
      (*i)->flag = (*i)->missing;
      ++i;
    }

  /////////////////////////////////
  // Pretend we have covariates

  par::clist = true;
  par::clist_number = 1;

  PP->clistname.resize(1);
  PP->clistname[0] = "PROFILE";

  for (int i=0; i< PP->n; i++)
    {
      Individual * person = PP->sample[i];
      person->clist.resize(1);
    }

}



vector_t Set::fitLDSetTest( vector_t & singleSNP, bool save )
{

  int ns = snpset.size();

  vector_t results(ns,0);

  if ( save )
    {
      numSig.resize(ns,0);
      selectedSNPs.resize(ns);
    }


  ///////////////////////////////////////////
  // Down-weight under true model?

  if ( save && par::fix_lambda )
    {

      PP->printLOG("Downweighting observed statistics in set-test by a factor of "
		   + dbl2str( par::lambda ) + "\n");

      vector_t::iterator i = singleSNP.begin();
      while ( i != singleSNP.end() )
	{
	  *i = (*i) / par::lambda;
	  ++i;
	}
    }



  ///////////////////////////////////////////
  // Consider each set

  for (int i=0; i<snpset.size(); i++)
    {

      int nss = snpset[i].size();


      ///////////////////////////////////////////
      // Extract all single SNP statistics for
      // this set

      vector<SetSortedSNP> t(nss);

      for (int j=0;j<nss;j++)
	{
	  t[j].chisq = singleSNP[ snpset[i][j] ] ;
	  t[j].l = j;
	}


      ////////////////////////////////////////////
      // Sort them (by decreasing chisq statistic)

      sort(t.begin(),t.end());



      /////////////////////////////
      // Extract a score from this
      // set

      double score = 0;
      int inSet = 0;
      int isSig = 0;

      vector<int> selected(0);


      // Step through SNPs sequentially, adding to score

      for( vector<SetSortedSNP>::reverse_iterator p = t.rbegin(); p!=t.rend(); p++)
	{

	  // Is this score already too large?

	  if ( inSet == par::set_max )
	    {
	      break;
	    }

	  // Get SET-centric SNP code

	  int j = p->l;

	  // Are there any SNPs worth adding?

	  if ( p->chisq < par::set_chisq_threshold )
	    {
	      break;
	    }


	  // Record this SNP as significant
	  if ( save )
	    ++isSig;

	  // Is this SNP correlated to a SNP already in the list?

	  set<int> & ls = ldSet[i][j];
	  bool hasProxy = false;

	  for (int k=0; k<selected.size(); k++)
	    {
	      set<int>::iterator d = ls.find(selected[k]);
	      if ( d != ls.end() )
		{
		  hasProxy = true;
		  break;
		}
	    }

	  // Advance to next potential SNP

	  if ( hasProxy )
	    continue;

	  // Otherwise, add this SNP to the score

	  score += p->chisq;
	  ++inSet;
	  selected.push_back( j );

	}


      ///////////////////////////////////////////
      // Do we want to save anything here?

      if ( save )
	{
	  numSig[i] = isSig;
	  selectedSNPs[i] = selected;
	}


      ///////////////////////////////////////////
      // Statistic is the mean test statistic per
      // selected SNP

      results[i] = inSet>0 ? score/(double)inSet : 0 ;
    }

  return results;
}



vector_t Set::fitStepwiseModel()
{

  if ( par::SNP_major )
    PP->SNP2Ind();
  par::assoc_glm_without_main_snp = true;

  // We are using the conditioning SNPs list to swap
  // in and out the effects
  par::conditioning_snps = true;

  // Put a set of SNPs into the model

  // Allow
  //  Fixed covariates, as usual
  //  Handle all SNPs as conditioning SNPs
  //  but have a boolean vector that allows some to be fixed
  //

  vector_t results;

  for (int i=0; i<snpset.size(); i++)
    {
      cout << "considering SET " << PP->setname[i] << "\n";

      int ns = snpset[i].size();

      vector<bool> fixed(ns,false);
      vector<bool> inModel(ns,false);

      // Scan all SNPs not in the model, and add the best if above
      // threshold

      bool done = false;

      Model * bestModel = NULL;

      while ( ! done )
	{

	  int bestSNP = -1;
	  double lowestP = 1;

	  for (int j=0; j<ns; j++)
	    {

	      if ( fixed[j] || inModel[j] )
		continue;

	      PP->conditioner.clear();

	      // And now add this SNP
	      PP->conditioner.push_back( snpset[i][j] );

	      for (int j2=0; j2<ns; j2++)
		{
		  if ( fixed[j2] || inModel[j2] )
		    PP->conditioner.push_back( snpset[i][j2] );
		}

// 	      cout << "Testing model: ";
// 	      for (int j2=0; j2<PP->conditioner.size(); j2++)
//		cout << PP->locus[ PP->conditioner[j2] ]->name  << " ";

	      PP->glmAssoc(false,*PP->pperm);

	      // Conditioning test SNP will always be the first

	      // This function skips the intercept
	      vector_t pv = PP->model->getPVals();
	      double pval = pv[0];

	      if ( pval < lowestP && realnum(pval) )
		{

		  //		  cout << "Selecting this marker..." << pval << "\n";

		  // But do we really want to accept this based
		  // on the absolute threshold?

		  if ( pval < par::set_step_in )
		    {
		      if ( bestModel != NULL )
			delete bestModel;
		      bestModel = PP->model;
		    }
		  else
		    delete PP->model;

		  lowestP = pval;
		  bestSNP = j;
		}
	      else
		{
		  delete PP->model;
		}
	    }


	  if ( lowestP < par::set_step_in )
	    {
	      inModel[bestSNP] = true;
	    }
	  else
	    {
	      done = true;
	    }
	  // Do we need this still?
	  if ( bestSNP == -1 )
	    done = true;


	} // Conintue the stepwise procedure?


      // The final model is stored in bestModel

      // Or perhaps we did not find a model?
      if ( bestModel == NULL )
	continue;

      // Note: skips intercept
      vector_t pval = bestModel->getPVals();

      // But, annoyingly..., this includes intercept...
      //  hmm, should sort this out

      vector_t coef = bestModel->getCoefs();

      // Skip intercept
       for (int t = 1; t < bestModel->getNP(); t++)
	{
 	  cout << "fcoef   " << bestModel->label[t] << " " << coef[t] << "\t" << pval[t-1] << "\n";
 	}
      cout << "--------\n";

//       for (int j2=0; j2<ns; j2++)
//  	{
// 	  if ( inModel[j2] )
// 	    cout << PP->locus[ snpset[i][j2] ]->name << " (selected) \n";
// 	  if ( fixed[j2] )
// 	    cout << PP->locus[ snpset[i][j2] ]->name << " (fixed) \n";
// 	}
//       cout << "\n";
//       cout << "-----------\n";

      // Obtain full model p-valie
      //      results.push_back( bestModel->getStatistic() );


      // PLACE ALL THIS IN CONTEXT OF PERMTATION ALSO...

      if ( bestModel != NULL )
	delete bestModel;

      // Next set
    }

  return results;

}
