

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __LOGISTIC_H__
#define __LOGISTIC_H__

#include<vector>
#include "plink.h"
#include "model.h"

using namespace std;

class LogisticModel : public Model {

 public:

  LogisticModel(Plink *);
  ~LogisticModel() { };

  void setDependent();
  
  void fitLM();
  void fitUnivariateLM() { };
  void reset();
  void pruneY();
  vector_t getCoefs();
  vector_t getVar();
  vector_t getSE();
  void displayResults(ofstream &, Locus *);
  double getLnLk();
  vector_t getPVals();
  double getPValue();
  void HuberWhite();

  
 private:
 
  vector_t p; 
  vector<int> Y;
  vector_t V; // diagonal p(1-p)
  double chisq;

};


#endif
