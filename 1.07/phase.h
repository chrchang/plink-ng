

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __HAPPHASE_H_
#define __HAPPHASE_H__

class HaploWindow;
class Plink;

class FamilyTransmissions
{
public:
  int pt, pu, mt, mu;  

  bool operator< (const FamilyTransmissions & b) const
  {
    if ( pt < b.pt )
      return true;
    if ( pt > b.pt )
      return false;

    if ( pu < b.pu )
      return true;
    if ( pu > b.pu )
      return false;

    if ( mt < b.mt )
      return true;
    if ( mt > b.mt )
      return false;

    if ( mu < b.mu )
      return true;
    if ( mu > b.mu )
      return false;

    return false;
  }

};


class HaploPhase
{
 public:

  Plink & P;

  int ns; // Number of SNPs in haplotype (region)
  int nw; // Numbner of windows in region
  int actual_nw; // If a subset if analysed
  int nh; // Number of possible haplotypes
  int nt; // Number of downcoded haplotypes
  int nsh; // Number of possible stub haplotypes
  int np; // Number of phases, diploid
  int haploid_np; // As above, haploid (do we need this?)

  string hname; // Name of haplotype locus
  int test_hap; // To be imputed haplotype 

  bool X; // Sex chromosome code
  bool haploid; // Haploid chromosome code

  int cnt_f; // Number of founders to be phased
  int current; // Number of current haplotype being tested
  
  bool reference_only; // Only consider reference panel

  //////////////////////////////////////
  // Lists of SNP sets (regions)

  // List of SNPs in haplotypes
  vector<vector<int> > new_pred_locus;

  // List of 'tag' haplotypes
  vector<string> new_pred_allele;

  // List of weighted multimarker predictors
  vector<map<string,double> > new_pred_weighted_allele;

  
  //////////////////////////////////////
  // Region-wide haplotype information

  // Coding for each haplotype (and HaploWindow coding)
  vector<vector<bool> > hap;
  vector<vector<int> > hapi;

  // List of predictor SNP numbers
  intvec_t S;

  // Estimated haplotype frequencies
  vector_t f;

  // Individual posterior probabilities
  matrix_t pp;
	
  vector<vector<int> > hap1;
  vector<vector<int> > hap2;

  // Lookup table for haplotype number given SNPs
  map<vector<bool>,int> hapmapb;
  map<vector<int>,int > hapmap;
       	
	
  // Phase markers and frequencies
  vector<int> ph_hap1;
  vector<int> ph_hap2;
  vector<double> ph_freq;

  vector<int> haploid_ph_hap1;
  vector<double> haploid_ph_freq;


  // Whether individual has ambiguous phase for region
  // i.e. hap1[i].size() == 1 
  vector<bool> ambig;

  // Should we skip this person?
  vector<bool> include;

  // Downcoding?
  bool subhaplotypes;
  map<int,int> downcoding;
 
  ////////////////////////////////////////
  // For EM, region is split into windows

  vector<HaploWindow*> windows;
  
  int startWindow;
  int finishWindow;

  ////////////////////////////////////////
  // For EM, region is split into windows

  void enumeratePhasedWindows(int);
  bool makeWaplotype(vector<int> &, vector<int> &);

  ////////////////////////////////////////
  // (Regional) haplotype association 

  // Transmission/untransmission counts
  vector<double> trans;
  vector<double> untrans;
  vector<map<FamilyTransmissions,double> > phasemap;

  ////////////////////////////////////////
  // (Regional) haplotype imputation (ML)

  vector<Locus*> new_map;
  vector<Locus*> actual_map;
  vector<vector<bool> > new_one;
  vector<vector<bool> > new_two;


  ////////////////////////////////////////
  // In segment-tracking-mode, individuals

  int p1;
  int p2;
  bool homozyg;

  // Output files: haplotype frequencies
  ofstream HFRQ;
  ofstream HTEST;
  ofstream HIMPUTE;
  ofstream HPHASE;
  ofstream VPHASE;
  ofstream WGT;

  // Temporary storage for chi-sqs from haplotype tests
  // and odds ratio (haplotype specific tests)
  double result;
  double pvalue;
  double odds;
  double case_freq;    // also T:U 
  double control_freq;

  HaploPhase(Plink & P_) :
    P(P_)
    {
      ambig.resize(P.n, false);
      include.resize(P.n, true);
      pp.resize(P.n);
      hap1.resize(P.n);
      hap2.resize(P.n);
      X=haploid = false;
      subhaplotypes = false;
      useEmpiricalVariance = true;
      reference_only = false;
      nonfounders = false;
    }
	
  // Read list of tests/tags
  void readTagFile();

  // Make sliding window list of tests
  void makeSlidingWindow(string);

  // Set a specific set of SNPs based on a command line
  void setSpecificSNPs(string);

  // Display haplotype frequencies
  void calculateHaplotypeFrequencies();

  // Track shared haplotypes
  void trackSharedHaplotypes();
  void trackThisSegment();
  vector_t trackedIBS;
  vector<int> trackedN;

  // Make test set for haplotype tests
  map<int,int> makeTestSet(boolvec_t &, boolvec_t &);

  // Make subhaplotype identity set
  map<int,int> makeSubHaplotypeSet(boolvec_t &);

  // Return subhaplotype name, formatted
  string getSubHaplotypeName(boolvec_t &, boolvec_t &, int);

  // Perform haplotype tests
  vector_t performHaplotypeTests(bool,Perm&);

  // Impute all haplotypes
  void imputeAllHaplotypes();

  // Display haplotype phases
  void calculatetHaplotypePhases();

  // Verbose displays
  void verboseDisplayWindows(int i, bool use_ref = true );

  // Return of dosage for a single or set of haplotypes
  double dosage(int i, set<int> & h);
  set<int> makeSetFromMap(map<int,int> & h);

  void reset()
    {
      ns = nh = np = 0;
      test_hap = -1;
      for (int i=0; i<hap.size(); i++)
	hap[i].clear();
      hap.clear();
      hapmap.clear();
      S.clear();
      f.clear();
      validN = 0;
      pp.resize(P.n);
      hap1.resize(P.n);
      hap2.resize(P.n);
      ambig.resize(P.n);
      include.resize(P.n);
      subhaplotypes = false;
      downcoding.clear();
      //calculateDp = false;
      //nonfounders = false;
      reference_only = false;

      for (int i=0; i<P.n; i++)
	{
	  pp[i].clear();
	  hap1[i].clear();
	  hap2[i].clear();
	  ambig[i] = false;
	  include[i] = true;
	}
    }
	
  void name(string n)
    {
      hname = n;
    }

  string haplotypeName(int i);

  int nHap()
    {
      return nh;
    }

  double testHaplotypeFreq()
    {
      if (test_hap>=0)
	return f[test_hap];
      else
	return -1;
    }

  // Main routine to driver phasing (also set up for assoc/perm testing)
  vector_t phaseAllHaplotypes(bool, Perm&);

  // Given list of SNP numbers, set up all possible haplotypes (hap)
  void enumerateHaplotypes(vector<int>&);

  // Give list of haplotypes in each phase
  void enumerateAllPhases();

  // Set test haplotype (query hap with allele string)
  void setTestHaplotype(string);

  // Return possible haplotype list
  vector<string> returnHaplotypes(vector<int>&);

  // Don't include individuals missing too much information
  void includeIndividuals( int i );
	
  // Construct and determine inclusion for nonfounder
  void validateNonfounder(int, vector<bool> &, vector<bool> &);
  
  // Determine possible haplotype phases for an individual
  void enumerateNonfounderPhase(int,
				vector<bool>&,vector<bool>&,
				int,int,
				int,int,
				vector<int>&,vector<int>&);

  // Possible offspring haplotypes given parents? (NOT USED)
  bool consistentNonfounderPhaseGivenParents(int,int,int,int,int,int,int);

  // Possible offspring haplotypes given offspring genotypes?
  bool consistentNonfounderPhaseGivenGenotypes(vector<bool>&,vector<bool>&,int,int);

  // Possible offspring haplotypes given offspring genotypes, exception for X?
  bool consistentNonfounderMalePhaseGivenXGenotypes(vector<bool> &,vector<bool> & s2,int); 


  // Score transmissions
  void transmissionCount(int, map<FamilyTransmissions,double> & );
  void scoreTransmissions(int,int,int,int,int,int,vector<int>&,vector<int>&);

  // Get rid of unlikely phases
  void prunePhase(int);

  // Check for unlikely genotypes
  void queryGenotype(int);
  void queryThisGenotype(int,int,int,vector_t&);

  // Phase non-founders, given we have haplotype frequencies
  // and score/rescore for TDT
  void phaseAndScoreNonfounder(int);

  // Use offspring information to help resolve parental phase
  void resolveWithKids(int);

  // E-M algorithm
  void performEM_original();
  void performAlternEM();

  // Report haplotype phase for an individual
  void reportPhase();

  // Report haplotype phase for an individual, alternate format
  void reportPhaseWideFormat();

  // Report haplotype frequencues
  void reportHaplotypeFrequencies();

  // Impute most likely genotype given a set of windows
  void mainImputation();

  // Helper function in stitching together windows in imputation mode
  void updateForImputation();

  void imputeThisHaplotype(int);
  double imputeHaplotypes(int, bool&, bool&);
  vector_t imputeGenotype(int, int);

  /////////////////////////////////////////////
  // Convenience functions to report LD, freqs

  double rsq(int, int);
  double dprime(int, int);
  double rsq_internal(int, int);
  double rsq_internal(boolvec_t &, boolvec_t &, boolvec_t &, boolvec_t &);
  double freq(boolvec_t &, boolvec_t &);
  
  bool calculateDp;

  //////////////////////////////////////
  // Haplotype-based association tests

  void haplotypicCC(map<int,int> &, int, bool);
  void haplotypicWeightedCC();

  void haplotypicTDT(map<int,int> &, int, bool);
  void haplotypicWeightedTDT();

  void haplotypicQTL(map<int,int> &, int, bool);

  map<int,int> testSet;
  set<int> sets;  // as above, slightly diff. specification

  int validN; // Number of non-missing founders

  // Perform non-founder fill-in phasing?
  bool nonfounders;


  // Empirical variance of dosage

  bool useEmpiricalVariance;
  void calculateEmpiricalVariance(int);
  void calculateEmpiricalVariance(set<int>&);
  set<int> returnHaplotypeSet(boolvec_t &, boolvec_t &); 
  double empiricalVariance;
  double ratio;
  
  // TDT empirical variance stores
  vector_t transmissionX;
  vector_t transmissionX2;
  double transmissionTotal;

};



#endif
