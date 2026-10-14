

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __HELPER_H__
#define __HELPER_H__

#include <string>
#include <vector>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "plink.h"
#include "options.h"

template<class T>
inline const T SQR(const T a) {return a*a;}
 
template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}
 
template<class T>
inline const T MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}
 
template<class T>
inline const T SIGN(const T &a, const T &b)
        {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
 
template<class T>
inline void SWAP(T &a, T &b)
        {T dum=a; a=b; b=dum;}


class Plink;
class Individual;
class CSNP;

using namespace std;

void sizeMatrix(matrix_t &, int,int);
void sizefMatrix(fmatrix_t &,int,int);
void sizeTable(table_t & , int, int);

void NoMem();

vector<bool> nvec_bool();

class CArgs
{
 public:
  CArgs(int,char**);

  int count()
    { return n; }
  
  bool any()
    { return n > 1 ? true : false; }
  
  void fromScript(string);
  void fromPriorLog(string);
  bool find(string);
  string value(string);
  int value_int(string);
  double value_double(string);
  long unsigned int value_lui(string);
  void check_unused_options(Plink &);
  bool parseOptions(string,string);

  vector<string> value(string,int);
  vector<string> varValue(string);
  
  vector<string> a;
  
 private:
  int n;
  vector<bool> parsed;
  vector<bool> option;
  vector<string> root_command;
  vector<string> original;
  map<string,string> optionLabel;
};

vector<string> parse2str(string);
vector<int> parse2int(string);
string searchAndReplace(string,string,string);
vector<string> commaParse(string);



template <class T>
bool from_string(T& t, 
		 const std::string& s, 
		 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

string display(vector<string> &);
string displayLine(vector<string> &);

void error(string);
void shutdown();
void checkDupes(Plink&);
bool readString(FILE *,string &);
void summaryBasics(Plink&);

string relType(Individual *, Individual *);

typedef vector<vector<double> > matrix_t;
typedef vector<double> vector_t;

void display(matrix_t &); 
void display(vector_t &); 
void display(vector<int> &); 

double genotypingRate(Plink &, int);
bool identicalSNPs(Plink *, int, int);
vector<string> listPossibleHaplotypes(Plink &, vector<int>);

void geno2matrix(vector<int> & snps, matrix_t &, boolmatrix_t &,bool);
 
int getInt(string,string);
long unsigned int getLongUnsignedInt(string,string);
double getDouble(string,string);

void permute(vector<long int>&);
void permute(vector<int>&);

vector<double> FDR_BH(vector<double>&);

void affCoding(Plink &);
void removeMissingPhenotypes(Plink & );

string genotype(Plink &, int i, int l);
string genotype(Plink & P, Individual *, int);
string genotypeToFile(Plink &, int i, int l);

int getChromosomeCode(string);
string chromosomeName(int);
int getMarkerChromosome(Plink &,string);
int getMarkerNumber(Plink &,string);
string leftWindowEdge(Plink & P, int bp, int chr);
string rightWindowEdge(Plink & P, int bp, int chr);
vector<int> getChromosomeMarkerRange(Plink &, int);
bool seeChromosome(Plink &,int);
vector<int> getChromosomeRange(Plink &);
vector<int> getWindowRange(Plink &P, int);

vector<vector<int> > two_locus_table(int,int);

void makePersonMap(Plink&,map<string,Individual*>&);
void makeLocusMap(Plink&,map<string,int>&);

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);

string int2str(int);
string dbl2str(double,int prc = -1);
string dbl2str_fixed(double, int prc = -1);
string longint2str(long int);
std::string sw(std::string s , int n);
std::string sw(double d , int n);
std::string sw(double d , int f, int n);
std::string sw(int i , int n);

std::string itoa(int, int);


void checkFileExists(string);
void checkFileExists(vector<string>);
bool doesFileExist(string);
bool compressed(string);
vector<string> tokenizeLine(ifstream&);
vector<string> tokenizeLine(string);
vector<string> tokenizeLine(ifstream &,string);

void defineDogChromosomes();
void defineMouseChromosomes();
void defineCowChromosomes();
void defineSheepChromosomes();
void defineHorseChromosomes();
void defineHumanChromosomes();
void defineRiceChromosomes();

vector<bool> vif_prune(vector<vector<double> > , double threshold,vector<int>&);
vector<vector<double> > calcSetCovarianceMatrix(vector<int> & nSNP);

void smoother(Plink & P, 
	      vector_t & input,
	      int n, 
	      vector_t & output1,
	      vector_t & output2,
	      vector<int> & count);

map<string, set<Range> > readRange(string);

double modelComparisonPValue(Model * alternate, Model * null);

set<Range*> rangeIntersect(Range & r1, map<string, set<Range> > & ranges);
set<Range*> mapRanges2SNP(int l, map<string, set<Range> > & ranges);
int2 mapSNPs2Range(Plink & P, const Range * range);
void makeScaffold(Plink & P);

void mapRangesToSNPs(string,
		     map<string, set<Range> > & ranges, 		     
		     map<int,set<Range*> > & snp2range);


map<string, set<Range> > filterRanges(map<string, set<Range> > & ranges, string filename);

#endif
