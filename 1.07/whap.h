

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2007 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __WHAP_H__
#define __WHAP_H__

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "plink.h"

class Plink;
class HaploPhase;

using namespace std;


class ChapModel {
  
 public:
  

  string model;
  
  map<int,int> haploGroup;
  vector<bool> masked_conditioning_snps;
  
  vector_t coef;
  vector_t se;
  vector_t p;
  vector<string> label;

  vector< set<int> > group;

  double lnLk;
  int df;
  double rsq;
  
  void buildFromSNPs(string);
  void buildFromGroups(string);
  bool haplotypesInSameGroup(int,int);
};



class Chap {

 public:
  

  ChapModel * alternate;
  ChapModel * null;
  
  ChapModel * current;

  Plink * P;
  HaploPhase * H;

  Chap(Plink * p_, HaploPhase * h_)
    { 
      P = p_;
      H = h_;
    }
  
  void determineTestType();
  
  void setModels(ChapModel &, ChapModel &);
  void build(ChapModel &);
  void setSNPList(vector<int> &, ChapModel &);

  
  bool isNested();


};

#endif
