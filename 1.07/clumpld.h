

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2007 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#ifndef CLUMP_LD_H_
#define CLUMP_LD_H_

#include <string>
#include <vector>
#include <map>

#include "plink.h"
#include "helper.h"


using namespace std;

class ResultTrio
{
public:
  
  double p;  // p-value
  string annot; // Annotation
  int f;     // which results file?
  string s;  // SNP name
  
  bool operator< (const ResultTrio & p2) const
  {
    return ( p < p2.p );
  }
  
};

class ClumpPair {
 public:
  string snp;
  int f;

  bool operator< (const ClumpPair & p2) const
  {
    if ( snp == p2.snp ) return ( f < p2.f );
    return ( snp < p2.snp );
  }

};

class ClumpResults {
 public:
  double p;
  string annot;
  
  bool operator< (const ClumpResults & p2) const
  {
    return ( p < p2.p );
  }
};

class clump_LD
{
public:

  Plink * P;
  HaploPhase * hp;

  //will be user defined

  double pval_cutoff; 
  double ld_distance;
  double second_pval_cutoff;
  float r2_cutoff;
  
  map<string, bool> clumped;
  vector<string> snps;
  vector<double> pvals;
  map<ClumpPair, ClumpResults>  assoc_results;

  vector<string> filename;
  
  // constructer
  clump_LD(Plink*,HaploPhase*,
	   double, double, double, float);
  
  // accessors
  void set_pval( double );
  void set_second_pval( double );
  void set_ld( double );
  void set_r2( double );
  
  // methods
  vector<ResultTrio> read_assoc_file(string);
  void clump();
  string allelePairs(int,int);
};


#endif /*CLUMP_LD_H_*/
