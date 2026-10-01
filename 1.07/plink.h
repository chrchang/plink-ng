

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////



#ifndef __PLINK_H__
#define __PLINK_H__

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <new>

#include "zed.h"

class CArgs;
class Perm;
class Set;
class Family;
class Cluster;
class HaploPhase;
class Locus;
class WMLocus;
class Individual;
class CSNP;
class Model;
class Chap;
class Variant;
class GVariant;

using namespace std;

typedef vector<vector<int> > table_t;
typedef vector<vector<double> > matrix_t;
typedef vector<double> vector_t;
typedef vector<bool> boolvec_t;
typedef vector<vector<bool> > boolmatrix_t;
typedef vector<int> intvec_t;
typedef vector<vector<double> > fmatrix_t;
typedef vector<float> floatvec_t;

typedef vector<Individual*>::iterator iIndividual;
typedef vector<Locus*>::iterator iLocus;
typedef vector<CSNP*>::iterator iSNP;
typedef vector<bool>::iterator iAllele;

class int2 {
 public:
  int p1;
  int p2;
  int2() { p1=p2=0; }
  int2(int a, int b) { p1=a; p2=b; }
  bool operator< (const int2 & b) const
    {
      return (p1 < b.p1 || (p1 == b.p1 && p2 < b.p2) );
    }
  bool operator== (const int2 & b) const
    {
      return (p1 == b.p1 && p2 == b.p2 );
    }

};

class double2 {
 public:
  double p1;
  double p2;
  double2() { p1=p2=0; }
  double2(double a, double b) { p1=a; p2=b; }
  bool operator< (const double2 & b) const
    {
      return (p1 < b.p1 || (p1 == b.p1 && p2 < b.p2) );
    }
  bool operator== (const double2 & b) const
    {
      return (p1 == b.p1 && p2 == b.p2 );
    }

};

class Pair2
{
public:

  double p;
  int l;

  bool operator< (const Pair2 & p2) const
  {
    return ( p < p2.p );
  }
};

class indivPair
{
public:

  Individual * p1;
  Individual * p2;

  bool operator< (const indivPair & b) const
  {
    return (p1 < b.p1 || (p1 == b.p1 && p2 < b.p2) );
  }
};


class Range {
public:
  int chr;
  int start;
  int stop;
  string name;
  int group;
  static map<string,int> groupNames;

  Range() { }
  Range(int p1, int p2, int p3, string p4)
    {
      chr = p1;
      start = p2;
      stop = p3;
      name = p4;
    }

  bool operator< (const Range & b) const
  {
    if ( chr < b.chr )
      return true;
    else if ( chr > b.chr )
      return false;
    if ( start < b.start )
      return true;
    else if ( start > b.start )
      return false;
    return stop < b.stop;
  }

  bool operator== (const Range & b) const
  {
    return chr == b.chr && start == b.start && stop == b.stop;
  }

};


class WMLocus {
 public:
  WMLocus() { chr=0; name=""; }

  void reset()
    {
      allele.clear();
      weight.clear();
    }

  int chr;
  string name;

  vector<string> allele;
  vector<double> weight;
};


class Individual {
 public:

  Individual() {
    fid=iid=pat=mat="";
    ip=im=-1;
    sex=false; phenotype=-9;
    sexcode="";
    aff=false;
    covar=-9;
    bcovar=false;
    clist.resize(0);
    clistMissing.resize(0);
    plist.resize(0);
    plistMissing.resize(0);
    missing=false;
    missing2=false;
    flag=true;
    one.resize(0);
    two.resize(0);
    sol=0;
    founder=true;
    pp=pm=NULL;
    family=NULL;
    kids.resize(0);
    pperson=this;
    T=W=B=0;
    gvar.resize(0);
  }

  string fid;
  string iid;

  // Parental codes
  string pat;
  string mat;

  // Pointers to parents
  Individual * pp;
  Individual * pm;

  // Parent slot number
  int ip;
  int im;

  // Relatedness functions
  int countMeioses(Individual*);

  // Permuted self
  Individual * pperson;

  // Children (pointers, slot numbers)
  vector<Individual*> kids;
  vector<int> ikids;

  bool sex;
  string sexcode;
  double phenotype;
  bool aff;
  double covar;
  bool bcovar;

  vector_t clist; // multiple covariates
  vector<bool> clistMissing;

  vector_t plist; // multiple phenotypes
  vector<bool> plistMissing;

  bool missing;
  bool missing2;
  bool flag;

  int sol;
  bool founder;
  Family * family;

  // SNP data
  vector<bool> one; // Person-major mode genotypes
  vector<bool> two;

  vector<bool>::iterator i1;
  vector<bool>::iterator i2;

  // Generic variant data
  vector<GVariant*> gvar;



  // Weighted, multi-allelic single marker

  WMLocus wmlocus;

  // For QFAM, within and total scores (temporary variables)
  double T;
  double B;
  double W;
};


// Main genotype storage, ordered by SNP
class CSNP
{
 public:

  vector<bool> one; // SNP-major mode genotypes
  vector<bool> two;

};

class Cluster
{
 public:
  vector<Individual*> person;
};

class Family
{
 public:

  Family()
    {
      include = false;
      parents = false;
      discordant_parents = false;
      singleton = false;
      sibship = false;
      TDT = false;
      pat = mat = NULL;
      kid.clear();
    }

  void copy(const Family & rhs)
    {
      include = rhs.include;
      pat = rhs.pat;
      mat = rhs.mat;
      kid.clear();
      for (unsigned int c=0; c<rhs.kid.size(); c++)
	kid.push_back(rhs.kid[c]);
    }

  Family & operator = (const Family & rhs)
  {
    copy(rhs);
    return *this;
  }

  bool include;
  bool parents;
  bool sibship;
  bool discordant_parents;
  bool singleton;
  bool TDT;

  Individual * pat;
  Individual * mat;
  vector<Individual *> kid;

  // Between-family genotypic score
  double B;
};



class Locus {
 public:
  Locus() { chr=0; name=""; allele1=""; allele2=""; freq=0; pos=0; bp=0; nm=0; }

  int chr;
  string name;
  string allele1;
  string allele2;

  double freq;     // of allele1
  double pos;      // cM map positions
  int bp;          // base-pair position
  int nm;          // number of non-missing alleles

  // Copy constructor
  Locus(const Locus& h1) { copy(h1); }
  Locus & operator= (const Locus & h1) { copy(h1); return *this; }

  void copy(const Locus &h1)
    {
      chr = h1.chr;
      name = h1.name;
      allele1 = h1.allele1;
      allele2 = h1.allele2;
      freq = h1.freq;
      pos = h1.pos;
      bp = h1.bp;
      nm = h1.nm;
    }

  bool operator< (const Locus & p2) const
    {
      return (chr < p2.chr || (chr == p2.chr && bp < p2.bp) );
    }

  bool operator== (const Locus & p2) const
    {
      return ( name == p2.name );
    }


};

namespace std
{
  template<>
    class less<Locus*> {
    public:
    bool operator()(Locus const* p1, Locus const* p2) const
      {

	// Locus comparison based first on distance,
	// but then pointers in case we have a degenerate map
	// file (i.e. so we can sort on position, but so that
	// set<Locus*> still works

	if(!p1)
	  return true;
	if(!p2)
	  return false;

	if (p1->chr < p2->chr)
	  return true;

	if (p1->chr > p2->chr)
	  return false;

	if (p1->bp < p2->bp)
	  return true;

	return false;

      }
  };
};

class MainExitException {
 public:
  void exitMessage()
    {
      cout << "Exception: ending...\n";
    }
};



class Z {
 public:
  Z() { z0=z1=z2=0;}

  double z0;
  double z1;
  double z2;
};


class CInfo {
 public:
    int lstart;
    int lstop;
    int bpstart;
    int bpstop;
};

class Segment {
 public:
  Segment()
    {
      start = finish = 0;
      p1 = p2 = NULL;
      count = baseline = freq = type = sites = 0;
      score = 0.0;
    }

  int start; // based on map position [0..nl_all]
  int finish;
  Individual * p1;
  Individual * p2;

  // Generics

  int count;
  int baseline;

  double weightedCount;
  double weightedBaseline;

  int freq;
  int type;
  int sites;
  double score;

  // Just base for CNVs for now (i.e. only consider p1)

  bool operator< (const Segment & b) const
  {
    if ( start < b.start ) return true;
    if ( start > b.start ) return false;
    if ( finish < b.finish ) return true;
    if ( finish > b.finish ) return false;
    if ( p1 < b.p1 ) return true;
    return false;
  }

  bool operator== (const Segment & b) const
  {
    return ( start == b.start &&
	     finish == b.finish &&
	     p1 == b.p1 );
  }

};


class ZZ {
 public:
  ZZ() { z00=z01=z02=0;
        z10=z11=z12=0;
        z20=z21=z22=0; }
  double z00;
  double z10;
  double z20;
  double z01;
  double z11;
  double z21;
  double z02;
  double z12;
  double z22;

};


class Plink
{
 public:
  Plink() {
    sample.resize(0);
    locus.resize(0);
    phenotype.resize(0);
    clistname.resize(0);
    plistname.resize(0);
    m1.resize(0);
    m2.resize(0);
    pos.resize(0);
    pihat_G.resize(0);
    warnings=false;
    n=0;
    nl_all=0;
    ngvar=0;
    cnt_f=npheno=nl=0;
    nk=1;
    kname.resize(1,"0");
    phenotype_name = "";
    scaffold.clear();
  }

  // Genotype/phenotype per individual file
  vector<Individual*>            sample;

  // SNP information (ordered by SNP/individual)
  vector<CSNP*>                   SNP;

  // Locus information
  vector<Locus*>                 locus;

  // Marker scaffold
  map<int,CInfo> scaffold;

  // Genetic variant information
  vector<Variant*>               gvar;

  // Family data
  vector<Family*>                family;

  // Phenotype names
  string                         phenoLabel;
  vector<string>                 clistname;
  vector<string>                 plistname;

  // number of individuals, pairs
  int n;       // total number of individuals
  int cnt_f;   // number of founders
  int npheno;  // number of individuals with informative phenotypes
  int np;
  int nl_all;  // all loci
  int ngvar;   // generic variants (non-SNP)
  int nl;      // test loci
  int nk;      // number of clusters

  string phenotype_name;

  // Generic output file
  ofstream OUTFILE;
  ZOutput ZOUTFILE;

  // Were any warnings set?
  bool warnings;

  // Cluster names
  map<string,int> kmap;
  vector<string> kname;
  vector<Cluster*> klist;

  // List of oblig-missing SNP/clusters
  set<int2> oblig_missing;

  // Conditioning SNPs, and mask
  vector<int> conditioner;
  vector<bool> conditioner_mask;

  // Skip the pair if not informative
  vector<bool> skip_pair;

  // String for current multipoint pair IDs
  string pairid;

  // Singlepoint locus-specific IBD, singlepoint
  // for each pair (one at a time)
  vector<Z> Zlocus;

  // Store genome-wide IBD only for informative pairs
  vector<Z> saved_IBDg;

  // Multipoint map variables
  vector<int> m1;    // left flanking marker
  vector<int> m2;    // right flanking marker
  vector<double> pos; // relative position between

  // Final matrices: pihats and squared differences, cross-products
  vector< vector<double> > pihat; // row=marker, col=pair

  // Segments; CNVs, ROHs, IBD segments
  vector<Segment> segment;        // condensed IBD segment record
  set<Range> geneList;            // Genic/regional intersection range list
  map<Range,set<Segment> > gene2segment;  // Map of segments per gene

  vector<int> indivSegmentGroup;  // allelic group for each segment(per-ind)
  vector<double> pihat_G;         // global pi-hat
  set<int2> related;              // set of T pairs above threshold
  vector<double> phenotype;       // SD or CP
  vector<int> pair1;     // First member of pair
  vector<int> pair2;     // Second member of pair
  vector<int> in_anal;   // Unique'd list of inds in regression

  // Variances, means
  double m_phenotype;
  double v_phenotype;
  double prev_bt;
  vector<double> m_pihat;
  vector<double> v_pihat;

  // Expected frequencies for IBS|IBD
  double E00, E10, E20;
  double E01, E11, E21;
  double E02, E12, E22;

  // T matrix elements
  double T00, T01, T02;
  double T10, T11, T12;
  double T20, T21, T22;

  // Storage for genome-wide max(r^2) values
  vector<double> maxr2;

  // Association SETs
  vector<string> setname;
  vector<vector<int> > snpset;

  // Storage for original results
  vector<vector<double> > original;

  // IBS matrix and cluster variables
  vector<vector<double> > mdist;  // IBS metric
  double pv;                      // temporary holder of p-value
  double dst;                     // temporary holder of IBS
  double pvIBS0;                  // holder for IBS0 pvalue count
  double pvIBS2het;               // holder for IBS2 het/het count

  // Epistasis tests
  vector<bool> epi1;
  vector<bool> epi2;

  // Lists of individuals
  set<Individual*> gset1;
  set<Individual*> gset2;

  // Working variables for merge_mode >=6
  long int diff_overlap;
  long int diff_nonmissing_overlap;
  long int diff_concordant_overlap;

  // Association test p-value storage: # inds
  vector_t tcnt;

  // Cache for LD values in proxy-windows
  map<int2,double> proxyLD;

  // Segmental test help variables
  map<indivPair,int> segmentCount;
  map<indivPair,double> segmentLength;

  map<indivPair,double> segmentCount2;
  map<indivPair,double> segmentCount2Baseline;

  // Expected overlap
  vector_t expectedOverlap;
  vector_t expectedOverlapBaseline;

  // Pointer to permutation class
  Perm * pperm;


  ///////////////////////
  // Functions

  // Input/output functions

  void           readData();
  void           readDataLongFormat();
  void           readFamFile(string);
  void           readMapFile(string,vector<bool>&,vector<int>&,int&);
  void           readTransposedData();
  void           readGenericVariantData();
  void           outputGenericVariantFile();
  void           convertGenericVariantData();
  void           updateMapFile();
  void           updateFamFile();
  void           updateAlleles();
  void           readStdIn();
  void           mergeData();
  bool           reconcileMerge(int,int,string,string,bool,bool,ofstream&,map<string,int>&);
  void           mergeBinaryData();
  void           mergeList();
  void           dummyLoader();
  void           simulateSNPs();
  void           simulateSNPs_QT();

  bool           readPhenoFile();
  bool           readMultiplePhenoFile();
  bool           readCovariateFile();
  bool           readCovListFile();
  bool           readClusterFile(bool verbose=true);
  void           readConditioningList();
  void           readBinData();
  void           readSet();
  void           prettyPrintLengths();
  void           printLOG(string);
  void           outputSetFile();
  void           setAssocSummary();


  void           Ind2SNP();
  void           SNP2Ind();

  // Summary statistic / data cleaning functions

  void           filterSNPs();
  void           processGVAR();

  void           calcStratifiedAlleleFreqs();
  void           hardyWeinbergCheck();
  double         calcInbreeding(Individual *,int,int,ofstream&);
  void           sexCheck();
  void           calcFst();

  void           findAllHomozygousRuns(Perm &);
  void           findHomoRuns(Individual *,ofstream&);
  void           findHomoWindow(Individual *,ofstream&);
  void           summariseHomoRuns();
  void           findIBSRuns(Individual *,Individual *,ofstream&);
  void           findMissRuns(Individual *,ofstream&);

  void           groupSegmentsSpanning(int);
  void           displaySegmentsLong();
  void           displaySegmentsBED();

  // CNV segment functions
  void           setUpForCNVList();
  void           readCNVList();
  void           processCNVList();
  vector_t       glmCNVBurdenModel(Perm &, bool);

  // Helper functions
  bool           missingGenotype(int,int);
  bool           obligMissing(int,int);

  void           outputPermedPhenotypes(Perm &);

  void           countCNVPerRegion(vector<int>&,vector<int>&);
  void           initialiseGeneCountAssociation(Perm &);

  // Family-based functions

  void           parseTrios();
  void           makeFounders();
  void           makeMissingParents();
  void           linkRelateds(map<Individual*,int> &,
			      map<string,Individual*> &);

  void           checkMendel();
  void           pseudoCaseControl();

  vector<double> testTDT(bool,bool,
			 Perm &,
			 vector<bool> &,
			 vector<bool> & );
  void           perm_testTDT(Perm &);

  vector<double> testSibTDT(bool,bool,
			    Perm &,
			    vector<bool> &,
			    vector<bool> & );

  void           perm_testQTDT(Perm &);
  vector<double> calcQTDT(vector<int> &,
			  ofstream&,
			  bool,
			  Perm &,
			  vector<int> &,
			  vector<bool> &);


  vector<double> testTDT_POO(bool,bool,
			     Perm &,
			     vector<bool> &,
			     vector<bool> & );
  void           perm_testTDT_POO(Perm &);

  // IBS sharing test statistics

  vector<double> sharingIBSTest(Perm &);
  void           perm_sharingIBSTest(Perm &);



  ////////////////////////////////////
  // Main pointers to other classes

  // Haplotype phasing/testing

  HaploPhase *   haplo;

  // GLM models

  Model * model;

  // Conditional haplotype tests (WHAP)

  Chap * whap;

  // Set-based functions

  Set * pS;


  // PLINK Functions
  int            readInformative();
  int            calcInformative();
  void           writeInformative();
  void           displayGenomeWideInfo();
  void           testGenomeIBDByCovariate(Perm &);
  void           permutationIBSTest(Perm &);
  void           displayGMULTI(Individual *, Individual *, int, ofstream &);
  void           preCalcGenomeIBD();
  void           preCalcMultiPoint();
  void           preCalcSinglePoint();
  void           preCalcPhenotypes();

  Z              calcGenomeIBS(Individual *, Individual *);
  void           calcGenomeIBM(Individual *, Individual *);
  Z              calcGenomeIBD(Individual *, Individual *, Z);


  vector<Z>      calcLocusIBD(Individual *, Individual *, Z);
  vector<double> calcMultiPoint(vector<Z> &, Z, ofstream &);
  vector<double> calcSinglePoint(vector<Z> &, Z);


  short          calcPhenotypes(vector<double> &, Individual *p1, Individual *p2);
  void           calcRegression(int);
  vector<double> doRegression(int,vector<double>&);
  void           preCalcRegression_PHENO(vector<double>&);
  void           preCalcRegression_PIHAT();

  // Association tests

  void           calcAssociationWithPermutation(Perm&);
  void           calcAssociationWithBootstrap();

  void           perm_testGXE2(Perm &);
  vector<double> testQAssocGXE2(bool,Perm &);

  void           calcGXE(Perm&);

  void           perm_testHotel(Perm &);
  vector<double> calcHotel(bool, Perm &, Set &,int,int);


  void           calcMH();
  void           calcHomog();
  vector<double> calcMantelHaenszel_2x2xK(Perm &, bool);
  vector<double> calcMantelHaenszel_ORD(vector<int>&,vector<int>&,vector<int>&);
  vector<double> calcMantelHaenszel_IxJxK(vector<int>&,vector<int>&,vector<int>&);

  void           calcLDStatistics();
  void           calcPairwiseLD();
  double         correlation2SNP(int,int,bool,bool,bool useFlag=false);
  void           pruneLD();
  void           calcFlipScan();
  void           setReferenceAllele();

  map<Range,vector<int> > mkBlks(int, int );


  void           setFlagToCase();
  void           setFlagToControl();

  void           calcEpistasis();
  void           driverSCREEPI();

  vector<double> testMiss(Perm &,bool);
  void           performMisHapTests();

  void           proxyWrapper();
  void           performProxyTests(int);

  void           scoreIndividuals();
  void           calculateProfile(map<int,double> &, map<int,bool> &, vector_t &, matrix_t &, vector<int> &,vector<int> &);

  vector<double> testAssoc(int &, int &,
			   vector<int> &, vector<int> &, vector<int> &,
			   vector<double> &,
			   vector<double> &,vector<double> &,
			   vector<double> &,vector<double> &,
			   Perm &,
			   bool);

  vector<double> testQAssoc(bool, Perm &);
  vector<double> fullModelAssoc(bool, Perm &);
  void           displayQTMeans(ofstream &, int l);
  vector_t       glmAssoc(bool, Perm &);

  vector_t       conditionalHaplotypeTest(bool, Perm &);
  vector_t       glmHaplotypeTest(bool, Perm &);

  void           multcomp(vector<double>&,string);

  void           buildT(double,bool,double,double);
  void           setMarkerRange();

  void           buildCluster();
  void           generateMDS();
  void           groupGenome();
  void           summaryIBD();
  void           findSegments(int,int,vector_t &,ofstream &);
  void           summaryIBDsegments(Perm & perm);
  void           summaryIBSsegments(Perm & perm);
  void           indivSegmentSummary();
  void           indivSegmentSummaryCalc(map<indivPair,int>&, map<indivPair,double>&,bool,bool);
  void           readSegmentFile(ifstream &);
  void           readSegmentFileMinimal(ifstream &);
  void           readHomozygSegmentFile(ifstream &);
  void           segmentPermutationTest(Perm &,bool,string,vector<int>&,vector<int>&,vector<int>&);

  void           segmentIndividualTest(Perm &);
  vector_t       perm_segmentIndividualTest(Perm&,bool,int,int,map<Individual*,int>&);
  void           homozygousSegmentPermutationTest(Perm &,string,vector<int>&,vector<int>&);
  void           validateSegments();
  void           positionPermuteSegments();

  void           runTestCNVwithQT(Perm &);
  vector_t       testCNVwithQT(double,int,int,vector_t&,vector_t&,vector_t&);

  void           runTestCNVwithGLM(Perm &);
  vector_t       testCNVwithGLM(bool, Perm &, vector<int> & );

  void           displayGenomePV();

  void           extractExcludeSet(bool);
  void           removeIndividuals(bool);
  void           keep2SetsForGenome();

  void           filterQualSNPs();
  void           filterQualGenotypes();

  void           makePhenotype();

  void           filterOnCovariate();
  void           filterOnCase();
  void           filterOnControl();
  void           filterOnMale();
  void           filterOnFemale();
  void           filterOnFounder();
  void           filterOnNonFounder();

  void           attribFilterSNP();
  void           attribFilterInd();


  void           zeroOnCluster();
  void           setObligMissing();

  int            deleteSNPs(vector<bool>&);
  int            deleteSNPs(set<string>&);
  int            deleteSNPs(set<Locus*>&);
  int            deleteIndividuals(vector<bool>&);
  int            deleteIndividuals(set<Individual*>&);
  void           thinSNPs();

  int            keepSNPs(set<string>&);
  int            keepSNPs(set<Locus*>&);
  int            keepIndividuals(set<Individual*>&);

  void           flipStrand();
  void           alleleRecoding();

  void           display_recoded_PEDFILE();
  void           display_recoded_PEDFILE_transpose();
  void           display_recoded_PEDFILE_AD();
  void           display_recoded_LONG();
  void           display_recoded_MUTLIST();

  void           output_fastphase_format();
  void           output_bimbam_format();
  void           output_structure_format();

  void           display_listByAllele();
  void           display_twolocus();
  void           display_pairList();
  void           display_indivReport();
  void           write_BITFILE();
  void           write_covariates();
  void           write_clusters();
  void           write_snplist();
  bool           openBinaryFile(string,ifstream&);
  void           setTable();
  void           writeSetFile();
  void           tagMode();

  void           processDosageFile();

  void           displayGeneReport();
  void           annotateFile();
  void           metaAnalysis();

  void           webcheck(CArgs &);
  void           lookup();
  void           lookup2();
  void           Rfunc();
  void           cleanUp();


  // Misc help functions

  void setFlags(bool f)
    {
      vector<Individual*>::iterator person = sample.begin();
      while ( person != sample.end() )
	{
	  (*person)->flag = f;
	  person++;
	}
    }


  // Additional functions
  void permTestRareDistribution(Perm &);
  void elfBaseline();
  void displayRareRange();
  vector_t testRareDistribution(Perm &,bool,map<Range,int2> & ranges);

};

#endif
