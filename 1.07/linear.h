

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __LINEAR_H__
#define __LINEAR_H__

#include<vector>
#include "plink.h"
#include "model.h"

using namespace std;

class LinearModel : public Model {

 public:

  LinearModel(Plink *);
  ~LinearModel() { };

  void setDependent();

  void fitLM();
  void fitUnivariateLM();
  void pruneY();
  void standardise();

  void reset();
  vector_t getCoefs();
  vector_t getVar();
  vector_t getSE();
  vector_t getPVals();
  double getPValue();
  void HuberWhite();

  void displayResults(ofstream &, Locus *);

  double calculateRSS(); 
  double calculateRSquared();
  double calculateAdjustedRSquared();
  double calculateMallowC(LinearModel *);
  double calculateFTest(LinearModel *);

 private:

  vector_t Y;
  vector<int> C;

  vector<double> se;
  double chisq;

  vector<double> sig;
  vector<double> w;
  vector<vector<double> > u;
  vector<vector<double> > v;

  double varY;
  double meanY;


  double RSS;

  void function(const int i, vector<double> & p );  
  void setVariance();
};


#endif
