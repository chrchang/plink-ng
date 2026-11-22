
//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __HAPWINDOW_H_
#define __HAPWINDOW_H__

class HaploPhase;
class Plink;
class MultiLocusGenotype;

class HaploWindow
{

 public:

  int ns; // Number of SNPs in haplotype
  int nh; // Number of possible haplotypes
  int np; // Number of phases, diploid

  // Parent 'region', PLINK
  HaploPhase * haplo;
  Plink * P;

  // Start and stop positions (relative to region)
  int start, stop;

  // Haplotype frequencies
  vector_t f;

  // Window haplotype codes
  vector<vector<bool> > hap;

  // Stub codes for each haplotype (for quick lookup)
  vector<int> leftStub;
  vector<int> rightStub;

  // Lookup table for haplotype number given SNPs
  map<vector<bool>,int> hapmap;

  // List of SNP numbers
  intvec_t S;

  // Posterior probabilities, per individual
  matrix_t pp;

  // Haplotype phases, per individual
  table_t hap1;
  table_t hap2;

  // Ambiguous for this window?
  boolvec_t ambig;

  // Unamiguous haplotype counts
  vector_t uc;

  // Store count, reference individual
  set<MultiLocusGenotype*> genotypes;

  // Store which genoGroup a person belongs to
  vector<MultiLocusGenotype*> genoGroup;

  // Finished with this window?
  bool converged;
  bool left_passed;
  bool right_passed;


  // Convergence
  vector<bool> zero;
  double sampleLogLikelihood;
  int iter;

  ///////////////////////////////
  // Functions

  HaploWindow(HaploPhase *, Plink *);

  ~HaploWindow();

  void expandGenogroups();
  void enumerateGenogroups();
  void pruneGenogroups(double t=par::haplo_plem_window_prune_phase);
  void enumerateHaplotypes(intvec_t &);
  void setStubCodes();
  void performEM();
  void enumeratePhase(int);
  void prunePhase(int,double t=par::haplo_plem_window_prune_phase);
  void reportPhase();
  string haplotypeName(int i);

  // Get overlap frequency from a window
  vector_t leftStubFrequency();
  vector_t rightStubFrequency();

  void tallyUnambiguousCounts();

};

#endif
