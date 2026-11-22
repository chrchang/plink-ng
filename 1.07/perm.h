

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __PERM_H__
#define __PERM_H__

#include "options.h"

using namespace std;

class Perm { 
  
  int t;                 // Number of tests
  long int replicates;   // Basic number of replicates
  long int performed;    // Number of replicates actually performed

  bool count;            // Output counts, not p-values

  ofstream PDUMP;        // Verbose permutation dump
  bool dump_best;        // Record just best permutation result
  bool dump_all;         // Record all permutation results
   
  vector<int> order;     // For --rank, record order of original
  vector<int> reorder;     // For --rank, record order of original, reverse mapping

  /////////////////////////////
  // Gene-dropping permutation

  bool genedrop;
  map<Individual*,int> idmap; 

  ///////////////////////////////////////////
  // Standard phenotype-swapping permutation

  // Parameters for adaptive permutation
  bool adaptive;
  
  int min;               // Minimum number of permutations
  double zt;             // SD CI range (based on par::adaptive_alpha)
  int interval;          // Prune tests every (I+N*I2) permutations

  // Main storage
  vector<int> R;         // number of successes
  vector<int> maxR;      // number of genome-wide successes
  vector<long int> N;    // number of trials (adaptive)
  
  // Cluster information
  vector< vector<int> > s;
  int ns;                // number of clusters
  
  Plink & P;             // reference to Plink class

 public:

  Perm(Plink &); 
  
  void closeDUMP() 
    {
      if (dump_all || dump_best)
	PDUMP.close();
    }
  
  // For basic permutation
  vector<int> pheno;     // label-swapped phenotype
  vector<int> geno;     // label-swapped phenotype

  vector<bool> test;             // whether to stop with this test
  vector<bool> snp_test;         // whether to skip these SNPs in 
                                 // a set-based test


  void                setTests(int x);
  void                setAdaptiveSetSNPs(int x);

  void                originalOrder();
  void                setOriginalRanking(vector_t&);
  bool                finished();
  bool                update(vector<double>&, vector<double>&);
  bool                updateSNP(double,double,int);
  void                nextSNP();
  vector<double> &    report();
  int                 current_reps() { return performed; }
  int                 reps_done(int);
  double              pvalue(int);
  double              max_pvalue(int);
  int                 rank(int);
  void                permuteInCluster();
  void                setPermClusters(Plink &);

  void                preGeneDrop();
  void                geneDrop();
  void                dropAlleles(Plink &,
				  Individual*,
				  int,
				  int,
				  vector<bool>&,
				  vector<bool>&,
				  vector<bool>&,
				  map<Individual*,int> &);
  
  
};


#endif
