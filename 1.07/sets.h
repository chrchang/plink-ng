

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __SETS_H__
#define __SETS_H__

class SetSortedSNP
{
 public:
  double chisq;
  int l;
  bool operator< (const SetSortedSNP & s2) const
    {  return (chisq < s2.chisq); }  
};


class Set
{
  
 public:
  
  vector<int> s_min;
  vector<int> s_max;
  
  Set(vector<vector<int> > &);
  
  vector<vector<int> > & snpset;
  
  map<int,set<int> > setMapping;

  // Sum-statistics based set-based scores
  void empiricalSetPValues();
  void cumulativeSetSum_WITHOUTLABELS(vector<double>&,int);
  void cumulativeSetSum_WITHLABELS(Plink &, vector<double>&);
  
  
  // Profile-score based set-based tests
  void initialiseSetMapping();
  void profileTestSNPInformation(int,double);
  vector_t profileTestScore();
  void profileTestInitialise();
  
  vector<vector<int> > profileSNPs;
  vector<vector_t> profileScore;
  
  // Stepwise regression models
  vector_t fitStepwiseModel();

  // New, LD-aware single-statistic set test
  vector_t fitLDSetTest(vector_t&,bool);
  vector< vector<set<int> > > ldSet;
  vector<int> numSig;
  vector<vector<int> > selectedSNPs;

  // Helper functions  
  void sizeSets();
  void pruneSets(Plink&);
  void pruneMC(Plink &,bool,double);
  void dropNotSet(Plink &);
  void makeLDSets();

  // get better name ?? 
  vector<vector<string> > setsort;
  
  // All statistics
  vector<vector<vector<double> > > stat_set;
  
  // Empirical set-based p-values (p0, p1 and p2)
  vector<vector<vector<double> > > pv_set;
  vector<vector<double> > pv_maxG_set;
  vector<vector<double> > pv_maxE_set;

  // Include or drop this SNP (multi-collinearity)
  vector<vector<bool> > cur;
  
};

#endif
