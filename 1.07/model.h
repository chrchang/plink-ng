

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __MODEL_H__
#define __MODEL_H__

#include<vector>
#include "plink.h"

using namespace std;

class Model {
  
 public:
  
  Model();
  virtual ~Model() { };

  virtual void setDependent() = 0;
  virtual void pruneY() = 0;
  virtual void fitLM() = 0;
  virtual vector_t getCoefs() = 0;
  virtual vector_t getVar() = 0;
  virtual vector_t getSE() = 0;
  virtual vector_t getPVals() = 0;
  virtual void displayResults(ofstream &, Locus *) = 0;
  virtual void fitUnivariateLM() = 0;
  

  void setMissing();
  vector<bool> getMissing();
  void setMissing(vector<bool>&);
  void yokeMissing(Model *);
  void setHaploid();
  void setX();
  void setDominant();
  void setRecessive();
  void hasSNPs(bool);
  void addAdditiveSNP(int);
  void addDominanceSNP(int);
  void addHaplotypeDosage(set<int>&);
  void addSexEffect();
  bool isSexInModel();
  void addCovariate(int);
  void addInteraction(int,int);
  void buildDesignMatrix();  
  bool checkVIF();
  vector<bool> validParameters();
  bool isValid() { return all_valid; }
  double getStatistic();
  //  double getPValue();
  double linearHypothesis(matrix_t &, vector_t &);
  int Ysize() { return nind; }
  int getNP() { return np; } 
  void setValid() { all_valid = true; }
  vector<string> label;
  vector<int> order;
  vector<int> type;
  int testParameter;

  void noCluster();
  void setCluster();
  virtual void HuberWhite() = 0;
  
  // Independent variables (can be directly manipulated...)
  vector<vector<double> > X;

 protected:
 
  Plink * P;

  // Missing flag
  vector<bool> miss;

  int nind;
  int np;  // Main effects + interaction + intercept

  bool has_snps;

  vector<bool> xchr;
  vector<bool> haploid;

  bool sex_effect;
  
  vector<bool> valid;
  bool all_valid;

  vector_t coef;  // beta
  matrix_t S;     // Sigma


  // Term types
  
  enum terms { INTERCEPT, 
	       ADDITIVE, 
	       DOMDEV, 
	       HAPLOTYPE, 
	       SEX, 
	       COVARIATE, 
	       INTERACTION,
               QFAM };

  
  double buildIntercept();
  double buildAdditive(Individual *, int);
  double buildDominance(Individual *, int);
  double buildHaplotype(int, int);
  double buildSex(Individual *);
  double buildCovariate(Individual *, int);
  double buildInteraction(Individual *, int, vector_t &);
  double buildQFAM(Individual *);
  
  bool skip;

  // List of additive SNP effects
  // assuming SNP major mode

  vector<int> additive;

  int mAA;
  int mAB;
  int mBB;

  double mA, mB;
  
  // List of dominance deviation SNP effects

  vector<int> dominance;

  // List of covariates (clist)

  vector<int> covariate;

  // List of pairwise interactions
  // ( indexing previously specified components, 1,2,..)
  
  vector<int2> interaction;

  // List of sets of haplotypes
  
  vector<set<int> > haplotype;

  // Clustering information
  bool cluster;
  vector<int> clst;
  int nc;

};


#endif
