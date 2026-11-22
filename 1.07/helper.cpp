

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <ctime>

#include "helper.h"
#include "crandom.h"
#include "options.h"
#include "plink.h"
#include "perm.h"
#include "stats.h"
#include "nlist.h"
#include "model.h"
#include "logistic.h"
#include "linear.h"

#define FPMIN 1.0e-30

extern ofstream LOG;
extern Plink * PP;

vector<bool> nvec_bool()
{ 
  vector<bool> t(0);
  return t;
}

inline bool is_rare(Locus * loc) {
  if (loc->freq < 0 
      || loc->freq < par::min_af 
      || (1-loc->freq) < par::min_af)
    return true;
  else
    return false;
};

string display(vector<string> & s)
{
  string t = "";
  for (int i=0; i<s.size(); i++)
    t += s[i] + "\n";
  return t;
}

string displayLine(vector<string> & s)
{
  string t = "";
  for (int i=0; i<s.size(); i++)
    t += s[i] + "  ";
  t += "\n";
  return t;
}

void display(matrix_t & m) 
{
  if (par::silent) return;
  cout << "\n";
  for (int i=0; i< m.size(); i++)
    {
      cout << i << ")\t";
      for (int j=0; j<m[i].size(); j++)
	cout << m[i][j] << " ";	
      cout << "\n";
    }
  cout << "\n";
}

void display(vector_t & m) 
{
  if (par::silent) return;
  cout << "\n";
  for (int i=0; i< m.size(); i++)
    cout << i << ")\t" << m[i] << "\n";
  cout << "\n";
  cout << "\n";
}

void display(vector<int> & m) 
{
  if (par::silent) return;
  cout << "\n";
  for (int i=0; i< m.size(); i++)
    cout << i << ")\t" << m[i] << "\n";
  cout << "\n";
  cout << "\n";
}


CArgs::CArgs(int _n , char *argv[] )
{
  n = _n;
  a.push_back(argv[0]);
  parsed.resize(n,false);
  option.resize(n,false);
  for (int i=1 ; i < n ; i++ )
    a.push_back(argv[i]);
  original = a;
  
  // Valid option labels
  optionLabel.insert(make_pair("--lookup-gene2","LOOKUP"));
  optionLabel.insert(make_pair("--lookup-gene","LOOKUP"));

  optionLabel.insert(make_pair("--id-replace","IDHELP"));
  optionLabel.insert(make_pair("--id-match","IDHELP"));
  
  optionLabel.insert(make_pair("--meta-analysis","META"));
  optionLabel.insert(make_pair("--annotate","ANNOT"));
  optionLabel.insert(make_pair("--dosage","DOSAGE"));

}  



void CArgs::fromScript(string f)
{
  checkFileExists(f);
  ifstream INP(f.c_str());
  string buff;
  
  while( INP >> buff )
    {
      a.push_back(buff);
      parsed.push_back(false);
      option.push_back(false);
      n++;
      if (INP.eof()) continue;      
    }
  INP.close();
  original = a;
}

void CArgs::fromPriorLog(string f)
{

  checkFileExists(f);
  ifstream INP(f.c_str());
  string buff;

  // Points to note:
  //  1) check whether any flag has already been given:
  //    if it has, then do not add the flag (i.e. allow
  //    the possibility of over-writing values).  This
  //    will automatically take care of the multiple rerun
  //    issue, and will allow for a different naming, via
  //    the --out command which would be the most common use

  // Read log file up to where commands start

  while (1)
    {
      vector<string> tokens = tokenizeLine(INP);
      
      if ( tokens.size() == 3 &&
	   tokens[0] == "Options" && 
	   tokens[2] == "effect:" )
	break;

      if (INP.eof()) break;      
	    
    }

  while( 1 ) 
    {
      vector<string> tokens = tokenizeLine(INP);
      
      if ( tokens.size()==0 )
	break;
      
      // Are we over-riding this command / do we 
      // want to skip it? 

      if ( find(tokens[0]) )
	continue;
      
      for (int t=0; t<tokens.size(); t++)
	{
	  a.push_back(tokens[t]);
	  parsed.push_back(false);
	  option.push_back(false);
	  n++;
	}
      if (INP.eof()) break;      
    }

  INP.close();
  original = a;
}

bool CArgs::find(string s)
{  
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s) { parsed[i]=true; return true; }
  return false;
}  

string CArgs::value(string s)
{
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s && (i+1 < n) ) { parsed[i+1]=true; return a[i+1]; }
  error("Missing an argument for "+s);
} 


int CArgs::value_int(string s)
{
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s && (i+1 < n) ) 
      { 
	int x = getInt(a[i+1],s);
	parsed[i+1]=true; 
	return x; 
      }
  error("Missing an argument for "+s); 
} 

long unsigned int CArgs::value_lui(string s)
{
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s && (i+1 < n) ) 
      { 
	long unsigned int x = getLongUnsignedInt(a[i+1],s);
	parsed[i+1]=true; 
	return x; 
      }
  error("Missing an argument for "+s); 
} 

double CArgs::value_double(string s)
{
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s && (i+1 < n) ) 
      { 
	double x = getDouble(a[i+1],s);
	parsed[i+1]=true; 
	return x; 
      }
  error("Missing an argument for "+s);
} 

vector<string> CArgs::value(string s, int c)
{
  vector<string> r(0);
  
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s && (i+1 < n) ) 
      {
	for (int j=1;j<=c;j++) 
	  {
	    if ( (i+j) < a.size() )
	      {
		parsed[i+j]=true;
		r.push_back(a[i+j]);
	      }
	    else
	      error("Not enough arguments given for option: "+s+" ");
	  } 
      }

  if (r.size() != c) error("Not enough arguments given for option: "+s+" ");
  return r;
} 

vector<string> CArgs::varValue(string s)
{  
  vector<string> r(0);
  for (int i=0 ; i < n ; i++ )
    if (a[i] == s && (i+1 < n) ) 
      {
	for (int j=i+1;j<n;j++) 
	  {	    
	    // Stop at next command, a "+" sign, or end of command list
	    if ( a[j].substr(0,2)=="--" ) 
	      return r;
	    if ( a[j] == "+" )
	      return r;
	    r.push_back(a[j]);
	    parsed[j]=true;	    
	  }
      }  
  return r;
} 

vector<string> parse2str(string s)
{
  vector<string> y;
  string t="";

  for (int i=0 ; i < s.length() ; i++)
    if (s[i] == ',' || i == s.length()-1 )
      {
	if (i == s.length()-1) t += s[i];
	y.push_back(t);
	t = "";
      }
    else
      t += s[i];

  return y;
}

vector<int> parse2int(string s)
{
  vector<string> v = parse2str(s);
  vector<int> y;
  for (int i=0; i<v.size(); i++)
    y.push_back( atoi( v[i].c_str() ) );
  return y;
}

bool CArgs::parseOptions(string cmd, string opt)
{
  map<string,string>::iterator i = optionLabel.find(cmd);
  if ( i == optionLabel.end() )
    return false;
    
  // Get handle (creating it if needed)
  OptionSet * o = par::opt.addOption( i->second );
    
  NList n2(0);
  n2.setRangeChar(" ");
  n2.setDelimiter("=");

  vector<string> tok2 = n2.deparseStringList( opt );
  
  // Can only have 1 or fewer '=' signs
  if ( tok2.size() > 2 ) 
    return false;
  //    error("Problem with option " + opt + "\n");
  
  // Single flag?
  if ( tok2.size() == 1 || ( tok2.size()==2 && tok2[1]=="" ) )
    {
      vector<string> dummy;
      o->val.insert( make_pair( opt , dummy ) );
    }

  // Key-value(s) pair?
  if ( tok2.size() == 2 )
    {
      NList n2(0);
      n2.setRangeChar(" ");
      n2.setDelimiter(",");
      vector<string> r = n2.deparseStringList( tok2[1] );
      o->val.insert( make_pair( tok2[0] , r ) );
    }

  return true;
}


int getInt(string s, string a)
{
  int x;
  if(from_string<int>(x,s,std::dec))
    {
      return x;
    }
  else
    {
      error("Not valid integer argument for : "+a+" [ "+s+" ]");
    }
}

long unsigned int getLongUnsignedInt(string s, string a)
{
  long unsigned int x;
  if(from_string<long unsigned int>(x,s,std::dec))
    {
      return x;
    }
  else
    {
      error("Not valid integer argument for : "+a+" [ "+s+" ]");
    }
}

double getDouble(string s, string a)
{
  double x;
  if(from_string<double>(x,s,std::dec))
    {
      return x;
    }
  else
    {
      error("Not valid numeric argument for : "+a+" [ "+s+" ]");
    }

}


void CArgs::check_unused_options(Plink & P)
{

  // Any unused options get added as options

  string cmd = "";

  bool okay=true;

  for (int i=1; i<n; i++)  // argv[0] is command line
    {
      
      // Is this a command? 

      if ( a[i].substr(0,2)=="--" )
	cmd = a[i];
      
      // Is this a "+" symbol (end of var-args, before options)
      if ( a[i] == "+" )
	{
	  option[i] = true;
	  continue;
	}

      if (!parsed[i]) 
	{
	  if ( cmd == "" ) 
	    {
	      P.printLOG("** Unused command line option: " + a[i] + " \n");
	      okay=false;
	    }
	  else
	    {
	      // add as option to this command
	      if ( parseOptions( cmd , a[i] ) )
		option[i] = true;
	      else
		{
		  P.printLOG("** Unused command line option: " + a[i] + " \n");
		  okay=false;
		}
	    }
	}
    }

  if (!okay) error("Problem parsing the command line arguments.");


  // Otherwise, all is okay, but print a record of the options used to CERR
  // All options always start with "--"
  P.printLOG("Options in effect:");
  for (int i=1; i<original.size(); i++)
    {
      if (original[i].substr(0,2)=="--") P.printLOG("\n\t"+original[i]) ;
      else if ( option[i] ) 
	{
	  if ( i>0 && original[i-1] == "+" ) 
	    P.printLOG(" " + original[i]);
	  else
	    P.printLOG("\n\t  " + original[i]);	    
	}
      else P.printLOG(" "+original[i]); 
    }
  P.printLOG("\n\n");
}


void checkDupes(Plink & P)
{
  
  set<string> people;
  vector<Individual*>::iterator person = P.sample.begin();
  string errmsg;
  while (person != P.sample.end() )
    {
      string str = (*person)->fid + "_" + (*person)->iid;
      if ( people.find(str) != people.end() )
	errmsg += "Duplicate individual found: [ " + (*person)->fid + " " + (*person)->iid + " ]\n";
      else
	people.insert(str);
      person++;
    }

  if (errmsg.size()>0)
    {
      P.printLOG(" *** WARNING *** DUPLICATE INDIVIDUAL IDS FOUND *** \n");
      P.printLOG(errmsg+"\n");
    }
  people.clear();

  set<string> markers;
  vector<Locus*>::iterator loc = P.locus.begin();
  errmsg="";
  while (loc != P.locus.end() )
    {
      if ( markers.find( (*loc)->name ) != markers.end() )
	errmsg += "Duplicate marker name found: [ " + (*loc)->name + " ]\n";
      else
	markers.insert( (*loc)->name );
      loc++;
    }
  
  if (errmsg.size()>0)
    {
      P.printLOG(" *** WARNING *** DUPLICATE MARKERS FOUND *** \n");
      P.printLOG(errmsg+"\n");
    }
  
  markers.clear();
  
}

void error(string msg)
{
  cerr << "\nERROR: " << msg << "\n";
  LOG << "\nERROR: " << msg << "\n";  
  LOG.close();

  if (par::gplink)
    {
      ofstream GP((par::output_file_name+".gplink").c_str(),ios::out);
      GP << "1\n";
      GP.close();
    }
  
  PP->cleanUp();

  exit(1);
}

void shutdown()
{

  time_t curr=time(0);
  string tdstamp = ctime(&curr);

  if (!par::silent) cout << "\nAnalysis finished: " + tdstamp +"\n";
  LOG << "\nAnalysis finished: " + tdstamp +"\n";

  if (PP->warnings)
    {
      if (!par::silent) 
	cout << "*** One or more WARNINGS were issued (see this LOG file) ***\n";
      LOG << "*** One or more WARNINGS were issued (see this LOG file) ***\n";
    }

  LOG.close();

  if (par::gplink)
    {
      ofstream GP((par::output_file_name+".gplink").c_str(),ios::out);
      GP << "0\n";
      GP.close();
    }
  
  PP->cleanUp();

  exit(0);
}


void affCoding(Plink & P)
{
  // Create affection coding
  for (int i=0; i<P.n; i++)  
    {      
      Individual * person = P.sample[i];     
      if ( person->phenotype == 2 && ! person->missing ) 
	person->aff = true;      
      else
	person->aff = false;
    }
}


void summaryBasics(Plink & P)
{
  
  if (par::bt) 
    {
      
      int ncase = 0;
      int ncontrol = 0;
      int nmissing = 0;
      
      for (int i=0; i<P.sample.size(); i++)
	if ( P.sample[i]->missing ) 
	  nmissing++;
	else if ( P.sample[i]->phenotype == 1 )
	  ncontrol++;
	else if ( P.sample[i]->phenotype == 2 )
	  ncase++;

      P.printLOG("After filtering, " + int2str(ncase)+" cases, "
		 +int2str(ncontrol)+" controls and "
		 +int2str(nmissing)+" missing\n");
    }
  else
    {
      int nmissing = 0;
      
      for (int i=0; i<P.sample.size(); i++)
	if ( P.sample[i]->missing ) 
	  nmissing++;

      P.printLOG("After filtering, " + int2str(P.sample.size()-nmissing)+" individuals with non-missing status\n");
    }


  // Display sex counts
  int nmale = 0;
  int nfemale = 0;
  int nambig = 0;

  for (int i=0; i<P.n; i++)
    {
      Individual * person = P.sample[i];
      
      if (person->sexcode=="1")
	nmale++;
      else if (person->sexcode=="2")
	nfemale++;
      else 
	nambig++;
    }

  P.printLOG("After filtering, "+int2str(nmale)+" males, "+int2str(nfemale)
	     +" females, and "+int2str(nambig)+" of unspecified sex\n");
  
}


#define MISSING1(i,l) ( P.SNP[l]->one[i] && ( ! P.SNP[l]->two[i] ) ) 
#define MISSING2(i,l) ( P.sample[i]->one[l] && ( ! P.sample[i]->two[l] ) ) 

double genotypingRate(Plink & P, int l)
{

  // Because we do not store genotyping rate, provide this 
  // convenience function; 

  // Do *not* distinguish between obligatory and non-obligatory
  // missingness -- for the purpose of proxy-assoc, it is all the
  // same.

  int m = 0;

  if ( par::SNP_major ) 
    {
      for (int i=0;i<P.sample.size();i++)
	if ( MISSING1(i,l) ) m++;	 
    }  
  else
    {
      for (int i=0;i<P.sample.size();i++)
	if ( MISSING2(i,l) ) m++;	 
    }  
  
  return (double)m / (double)P.n ;
  
}


bool identicalSNPs(Plink * P, int l1, int l2)
{
  if ( par::SNP_major ) 
    {
      for (int i=0; i<P->n; i++)
	{
	  if ( P->SNP[l1]->one[i] != P->SNP[l2]->one[i] )
	    return false;
	  if ( P->SNP[l1]->two[i] != P->SNP[l2]->two[i] )
	    return false;
	}  
      return true;
    }
  else
    {
      for (int i=0; i<P->n; i++)
	{
	  if ( P->sample[i]->one[l1] != P->sample[i]->one[l2] ) 
	    return false;
	  if ( P->sample[i]->two[l1] != P->sample[i]->two[l2] ) 
	    return false;
	}  
      return true;
    }
}


vector<string> listPossibleHaplotypes(Plink & P, vector<int> S)
{
  
  vector<string> str;  
  unsigned int h=0;
  int ns = S.size();
  int nh = (int)pow((double)2,ns);
  vector<vector<bool> > hap;

  while(h<nh)
    {
      vector<bool> m1;
      
      unsigned int p=1;
      for (int s=0;s<ns;s++)
	{
	  if ( h & p ) m1.push_back(false);
	  else m1.push_back(true);
	  p <<= 1;
	}
      
      // Add to list of haplotypes
      hap.push_back(m1);
      
      // Consider next haplotype
      h++;     

   }

  for (int h=0; h<hap.size(); h++)     
    {
      string hstr;
      for (int s=0;s<ns;s++)
        if (!hap[h][s]) hstr += P.locus[S[s]]->allele1;
        else if (P.locus[S[s]]->allele2=="") hstr += "0";
	else hstr += P.locus[S[s]]->allele2;
      str.push_back(hstr);
    }
  return str;

}


bool readString(FILE * fp, string & s)
{
  
  bool done = false;
  s="";
  while (1)
    {
      char ch = fgetc(fp);
      if ( ch==' ' || ch == '\t' )
	{
	  if (done)
	    return true;
	}
      else if ( ch=='\n' || ch=='\r' || feof(fp) ) 
	{
	  if (done)
	    return true;	  
	  else
	    return false;
	}
      else
	{
	  s += ch;
	  done = true;
	}
    }
}


void removeMissingPhenotypes(Plink & P)
{
  vector<bool> del(P.sample.size());
  for (int i=0; i<P.sample.size(); i++)
    del[i] = P.sample[i]->missing;
  
  int n_removed = P.deleteIndividuals(del);

  if ( n_removed > 0 )
    P.printLOG(int2str(n_removed)+
	       " individuals removed because of missing phenotypes\n");  
}


void geno2matrix(vector<int> & snps, matrix_t & g, boolmatrix_t & m, bool dom)
{
  // return a S x N matrix coded 0,1,2, where m is missing (1=yes,0=no)
  m.clear();
  sizeMatrix(g,PP->n,snps.size());
  m.resize(PP->n);
  for (int p = 0 ; p < PP->n ; p++)
    m[p].resize(snps.size());
  
  for (int s = 0 ; s < snps.size() ; s++)
    {
      for (int p = 0 ; p < PP->n ; p++)
	{
	  bool s1 = par::SNP_major ? PP->SNP[snps[s]]->one[p] : PP->sample[p]->one[snps[s]];
	  bool s2 = par::SNP_major ? PP->SNP[snps[s]]->two[p] : PP->sample[p]->two[snps[s]];

	  if      ( (!s1) && (!s2) )  
	    {
	      ++g[p][s] = dom ? 1 : 2;
	    }
	  else if ( (!s1) && s2 ) 
	    {
	      ++g[p][s] = 1;
	    }
	  else if (  s1   && s2 ) 
	    {
	      ++g[p][s] = 0;
	    }
	  else
	    {
	      m[p][s] = true;
	    }
	}
    }
}



string genotypeToFile(Plink & P, int i, int l)
{
  
  // Return a genotype in suitable text format to be 
  // written to most output files
  
  string a1 = par::recode_12 ? "1" : P.locus[l]->allele1;
  string a2 = par::recode_12 ? "2" : P.locus[l]->allele2;

  bool s1 = par::SNP_major ? P.SNP[l]->one[i] : P.sample[i]->one[l];
  bool s2 = par::SNP_major ? P.SNP[l]->two[i] : P.sample[i]->two[l];

  if      ( (!s1) && (!s2) )  
    return par::recode_delimit+a1+par::recode_indelimit+a1;
  else if ( (!s1) && s2 ) 
    return par::recode_delimit+a1+par::recode_indelimit+a2;
  else if (  s1   && s2 ) 
    return par::recode_delimit+a2+par::recode_indelimit+a2;
  else 
    return par::recode_delimit + par::out_missing_genotype
      + par::recode_indelimit+par::out_missing_genotype;

  return "?";
}


string genotype(Plink & P, int i, int l)
{
  string delimit = "/";
  string g;
  Locus * loc = P.locus[l];

  if (par::SNP_major)
    {
      CSNP * s = P.SNP[l];
      if      ( (!s->one[i]) && (!s->two[i]) )  
	g = loc->allele1 + delimit + loc->allele1;
      else if ( (!s->one[i]) && s->two[i]) 
	g = loc->allele1 + delimit + loc->allele2;
      else if (  s->one[i]   && s->two[i]) 
	g = loc->allele2 + delimit + loc->allele2;
      else g = par::missing_genotype + delimit + par::missing_genotype;
    }
  else
    {
      Individual * person = P.sample[i];

      if      ( (!person->one[l]) && (!person->two[l]) )  
	g = loc->allele1 + delimit + loc->allele1;
      else if ( (!person->one[l]) && person->two[l]) 
	g = loc->allele1 + delimit + loc->allele2;
      else if (  person->one[l]   && person->two[l]) 
	g = loc->allele2 + delimit + loc->allele2;
      else g = par::missing_genotype + delimit + par::missing_genotype;
    }
  return g;
}


string genotype(Plink & P, Individual * person, int l)
{
  string delimit = "/";
  string g;
  Locus * loc = P.locus[l];
  
  if      ( (!person->one[l]) && (!person->two[l]) )  
    g = loc->allele1 + delimit + loc->allele1;
  else if ( (!person->one[l]) && person->two[l]) 
    g = loc->allele1 + delimit + loc->allele2;
  else if (  person->one[l]   && person->two[l]) 
    g = loc->allele2 + delimit + loc->allele2;
  else g = par::missing_genotype + delimit + par::missing_genotype;
  
  return g;
}

void permute(vector<long int> &a)
{
    // generate random permutation of 0..n-1
    // where n is a.size();

    const long int n = a.size( );
    
    for( long int i = 0; i < n; i++ )
        a[ i ] = i;
    
    for( long int j = 1; j < n; j++ )
    {
        long int pos = CRandom::rand(j+1);
        long int tmp = a[ j ];

        a[ j ] = a[ pos ];
        a[ pos ] = tmp;
    }
}


void permute(vector<int> &a)
{
  
  // Generate a random permutation of 0 
  // to n-1 where n is a.size();
  
  const int n = a.size( );
  
  for( int i = 0; i < n; i++ )
    a[ i ] = i;
  
  for( int j = 1; j < n; j++ )
    {
      int pos = CRandom::rand(j+1);
      int tmp = a[ j ];
      
      a[ j ] = a[ pos ];
      a[ pos ] = tmp;
    }
}



int getChromosomeCode(string chr)
{
  map<string,int>::iterator cc = par::chr_map.find(chr);
  return cc == par::chr_map.end() ? 0 : cc->second;
}

string chromosomeName(int c)
{
  if ( c < 0 || c >= par::chr_code.size() ) 
    return "0";
  return par::chr_code[c];
}

int getMarkerChromosome(Plink & P, string m)
{
  for (int i=0;i<P.locus.size();i++)
    if (P.locus[i]->name==m) return P.locus[i]->chr;
  return -1;
}

int getMarkerNumber(Plink & P, string m)
{
  for (int i=0;i<P.locus.size();i++)
    if (P.locus[i]->name==m) return i;
  return -1;  
}

bool seeChromosome(Plink & P, int c)
{
  for (int l=0;l<P.locus.size();l++)
    {
      if (P.locus[l]->chr==c) 
	return true;
      else if (P.locus[l]->chr>c)
	return false;
    }
  return false;  
}

vector<int> getChromosomeMarkerRange(Plink & P, int c)
{
  vector<int> m(2);
  m[0] = -1; // first
  m[1] = -1; // last
  
  for (int i=0;i<P.locus.size();i++)
    {
      if (P.locus[i]->chr==c) 
	{
	  if (i<=m[0] || m[0]==-1) m[0]=i;
	  if (i>=m[1] || m[1]==-1) m[1]=i;	  
	}      
    }
  return m;  
}

vector<int> getWindowRange(Plink &P, int s)
{
  vector<int> m(2);
  m[0] = s; // first SNP
  m[1] = s; // last SNP
  
  // move backwards
  int x=s;
  int chr=P.locus[s]->chr;
  int bp=P.locus[s]->bp;
  int win = (int)(par::window * 1000); // half window size in bases
  int nl = P.locus.size() - 1;

  while ( 1 ) 
    {
      if ( x==0 ) break;
      // Move one position, until on different chromosome, or outside window
      x--;
      if ( P.locus[x]->chr != chr ) { x++; break; }
      if ( bp - P.locus[x]->bp > win ) { x++; break; }
    }
  m[0]=x;
  x=s;
  while ( 1 ) 
    {
      if ( x== nl ) break;
      x++;
      if ( P.locus[x]->chr != chr ) { x--; break; }
      if ( P.locus[x]->bp - bp > win ) { x--; break; }
    }
  m[1]=x;
  return m;
}

vector<int> getChromosomeRange(Plink & P)
{
  vector<int> m(2);
  m[0] = -1; // first chromosome
  m[1] = -1; // last chromosome
  
  for (int i=0;i<P.locus.size();i++)
    {
      if (P.locus[i]->chr<=m[0] || m[0]==-1) m[0]=P.locus[i]->chr;
      if (P.locus[i]->chr>=m[1] || m[1]==-1) m[1]=P.locus[i]->chr;	  
    }
  return m;  
}

std::string int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string longint2str(long int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string dbl2str(double n, int prc)
{
  std::ostringstream s2;
  if ( prc > 0 ) 
    s2.precision(prc);
  s2 << n;
  return s2.str();
}


std::string dbl2str_fixed(double n, int prc)
{
  std::ostringstream s2;
  s2 << setiosflags( ios::fixed );
  if ( prc > 0 ) 
    s2.precision(prc);
  s2 << n;
  return s2.str();
}
std::string sw(std::string s , int n) 
{
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

std::string sw(double d , int n) 
{
  std::string s = realnum(d) ? dbl2str(d) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

std::string sw(double d , int f, int n) 
{  
  std::string s = realnum(d) ? ( f < 0 ? dbl2str(d,-f) : dbl2str_fixed(d,f) ) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

std::string sw(int i , int n) 
{
  std::string s = realnum(i) ? int2str(i) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

void NoMem()
{
  cerr << "*****************************************************\n"
       << "* FATAL ERROR    Exhausted system memory            *\n"
       << "*                                                   *\n"
       << "* You need a smaller dataset or a bigger computer...*\n"
       << "*                                                   *\n"
       << "* Forced exit now...                                *\n"
       << "*****************************************************\n\n";
  exit(1);
}


std::string itoa(int value, int base) {
	enum { kMaxDigits = 35 };
	std::string buf;
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.
	
	// check that the base if valid
	if (base < 2 || base > 16) return buf;
	
	int quotient = value;
	
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	
	// Append the negative sign for base 10
	if ( value < 0 && base == 10) buf += '-';
	
	std::reverse( buf.begin(), buf.end() );
	return buf;	
}

void checkFileExists(string f)
{

  ifstream inp;
  
  inp.open(f.c_str(), ifstream::in);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      string msg = "No file [ " + f + " ] exists.";
      error(msg);
    }
  inp.close();
  return;

}

void checkFileExists(vector<string> f)
{
  for (int k=0; k<f.size(); k++)
    checkFileExists(f[k]);
  return;
}

bool doesFileExist(string f)
{

  ifstream inp;
  
  inp.open(f.c_str(), ifstream::in);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      return false;
    }
  inp.close();
  return true;

}

bool compressed(string s)
{
  // Does this end in ".gz" or ".Z" ? 
  int l = s.size();
  if ( l > 3 && s.substr(l-3,3) == ".gz" ) return true;
  if ( l > 2 && s.substr(l-2,2) == ".Z" ) return true;
  return false;
}


vector<string> tokenizeLine(ifstream & F1)
{
  char cline[par::MAX_LINE_LENGTH];
  F1.getline(cline,par::MAX_LINE_LENGTH,'\n');
  string sline = cline;
  string buf; 
  stringstream ss(sline); 
  vector<string> tokens; 
  while (ss >> buf)
    tokens.push_back(buf);
  return tokens;
}

vector<string> tokenizeLine(string sline)
{
  string buf; 
  stringstream ss(sline); 
  vector<string> tokens; 
  while (ss >> buf)
    tokens.push_back(buf);
  return tokens;
}

vector<string> tokenizeLine(ifstream & F1,string d)
{
  char cline[par::MAX_LINE_LENGTH];
  F1.getline(cline,par::MAX_LINE_LENGTH,'\n');
  string sline = cline;
  NList nl(0);
  nl.setDelimiter(d);
  nl.setRangeChar(" ");
  return nl.deparseStringList( sline );
}

bool Plink::obligMissing(int i, int l)
{
  int2 p;
  p.p1 = l;
  p.p2 = sample[i]->sol; 
  return ( oblig_missing.find(p) != oblig_missing.end() );
}


bool Plink::missingGenotype(int i, int l)
{  
  bool s1 = par::SNP_major ? SNP[l]->one[i] : sample[i]->one[l];
  bool s2 = par::SNP_major ? SNP[l]->two[i] : sample[i]->two[l];
  return ( s1 && ! s2 );
}


void Plink::prettyPrintLengths()
{

  par::pp_maxfid = 4;
  par::pp_maxiid = 4;
  par::pp_maxsnp = 4;
  
  for (int i=0;i<n;i++)
    {
      if (sample[i]->fid.length() > par::pp_maxfid)
	par::pp_maxfid = sample[i]->fid.length() + 2;

      if (sample[i]->iid.length() > par::pp_maxiid)
	par::pp_maxiid = sample[i]->iid.length() + 2;
    }	  
  

  for (int l=0;l<nl_all;l++)
    if (locus[l]->name.length() > par::pp_maxsnp)
      par::pp_maxsnp = locus[l]->name.length() + 2;
  
}


vector<bool> vif_prune(vector<vector<double> > m , double threshold , vector<int> & varcode )
{
  
  // Number of variables
  int p = m.size();

  vector<bool> cur(p,true);
  
  // This only is needed if we have 2+ SNPs
  if (p<2) { 
    return cur;
  }

  vector<vector<double> > r = m;
  
  // Make 'm' a correlation matrix
  for (int i=0; i<p; i++)
    for (int j=0; j<p; j++)
      r[i][j] = m[i][j] / sqrt(m[i][i] * m[j][j]);

  // Number of excluded items
  int it = 0;
  
  
  // Any SNPs with zero variance should be automatically excluded
  for (int i=0; i<p; i++)
    if ( r[i][i] == 0 || !realnum(r[i][i]) )
      {
	cur[i] = false;
	it++;
      }

  
  // For any pair of perfectly correlated SNPs, exclude 1

  while(1)
    {
      
      bool done = true;
      for (int i=0; i<p-1; i++)
	{
	  if (cur[i]){
	    for (int j=i+1;j<p; j++)
	      {
		if (cur[j]) {
		  if ( fabs(r[i][j]) > par::prune_ld_r2 )
		    {
		      
		      if ( par::prune_ld_pairwise_maf )
			{
			  // Remove SNP with lower MAF
			  if ( PP->locus[ varcode[i] ]->freq < PP->locus[ varcode[j] ]->freq ) 
			    cur[i] = false;
			  else
			    cur[j] = false;			  
			  it++;
			  done = false;
			  break;			  
			}
		      else
			{
			  // Just remove first
			  cur[i] = false;
			  it++;
			  done = false;
			  break;			  
			}
		    }
		}
	      }
	  }
	}
      if (done) break;
    }
  
  // Skip VIF calculation?
  if (par::prune_ld_pairwise)
    return cur;

  // Calculate r^2 for each element versus all others
  // considering only the current non-pruned elements
  
  while (1)
    {
      
      // Build correlation matrix all included items
      vector<vector<double> > u;
      
      for (int i=0;i<p;i++)
	{
	  if ( cur[i] )
	    {
	      vector<double> mt;
	      
	      for (int j=0;j<p;j++)
		if ( cur[j] )
		  mt.push_back(r[i][j]);
	      
	      u.push_back(mt);
	    }
	}
      

      // Check enough markers left  
      if (u.size()<2) break;
      
      if (par::verbose)
	{
	  cout <<  " about to invert\n";
	  cout.precision(12);
	  for (int i=0; i<u.size(); i++)
	    {
	      for (int j=0; j<u[i].size(); j++)
		cout << u[i][j] << " ";
	      cout << "\n";
	      
	    }
	}
      
      // Get inverse
      bool flag = true;
      u = svd_inverse(u,flag);

      
      // Calculate VIFs
      double maxVIF = 0;
      int maxI;
      int cnt=0;
      for (int i=0;i<p;i++)
	if ( cur[i] )
	  {
	    
	    // r^2 = 1 - 1/x where x is diagonal element of inverted
	    // correlation matrix
	    // As VIF = 1 / ( 1 - r^2 ) , implies VIF = x

	    double vif = u[cnt][cnt];
	    
	    if ( maxVIF < vif )
	      {
		maxVIF = vif;
		maxI = i;
	      }
	    
	    cnt++;
	  }

      // How are we doing?
      if ( maxVIF > threshold ) 
	{
	  // exclude this item
	  cur[maxI] = false;	  
	}
      else
	{
	  break;
	}
      
      // Increase count of removed items
      it++;
      
      // Down to a single item or worse?
      if (it==p-1) break;
      
    }

  return cur;

}



vector<vector<double> > calcSetCovarianceMatrix(vector<int> & nSNP)
{
  
  int nss = nSNP.size();
  vector<vector<double> > var( nss );

  if ( nss == 0 )
    return var;
  
  for ( int i = 0; i < nss; i++ )
    {
      var[i].resize( nss );
    }
  
  // Use helper function to calculate correlation coefficient 
  // between two SNPs (that allows for haploid,diploid nature)
  
  // second flag 'true' indicates to return covariance term, not
  // correlation

  for (int i=0; i<nss; i++)
    for (int j=i; j<nss; j++)
      {
	var[i][j] = var[j][i] 
	  = PP->correlation2SNP( nSNP[i], nSNP[j] , false, true );	  
      }
  
  return var;
}



string leftWindowEdge(Plink & P, int chr, int bp)
{
  // Get nearest SNP 
  Locus * marker = NULL;
  int distance = -1;
    
  vector<Locus*>::iterator loc = P.locus.begin();
  while ( loc != P.locus.end() )
    {
      if ( (*loc)->chr == chr) 
	if ( (*loc)->bp >= bp )
	  if ( (*loc)->bp - bp < distance || ! marker )
	    {
	      distance =  (*loc)->bp - bp;
	      marker = *loc;
	    }
      loc++;
    }

  if (!marker) 
    error("Could not place marker for left window edge");

  return marker->name;  

}


string rightWindowEdge(Plink & P, int chr, int bp)
{

  // Get nearest SNP 
  Locus * marker = NULL;
  int distance = -1;
  
  vector<Locus*>::iterator loc = P.locus.begin();
  while ( loc != P.locus.end() )
    {
      if ( (*loc)->chr == chr) 
	if ( (*loc)->bp <= bp )
	  if ( bp - (*loc)->bp < distance || ! marker )
	    {
	      distance =  bp - (*loc)->bp;
	      marker = *loc;
	    }
      loc++;
    }

  if (!marker) 
    error("Could not place marker for right window edge");

  return marker->name;  
}


void Plink::setMarkerRange()
{

  // If chromosome code >0, implies a specific chromosome
  if (par::run_chr>0 && !par::position_window)
    {
      // Get first and last markers on this chromosome
      vector<int> m = getChromosomeMarkerRange((*this),par::run_chr);
      if(m[0]==-1 || m[1]==-1) error("--chr {chromosome} not found:"+int2str(par::run_chr));
      par::run_start = m[0];
      par::run_end = m[1];
    }
  else if (par::position_window)
    {
      // Physical position specified (chromosome and range)
      par::m1 = leftWindowEdge(*this, par::run_chr, par::from_window);
      par::m2 = rightWindowEdge(*this, par::run_chr, par::to_window);
      par::run_start = getMarkerNumber( (*this), par::m1 );
      par::run_end = getMarkerNumber( (*this), par::m2 );	
    }
  else
    {
      // Two SNPs specified (or a SNP and a range)

      // If a specific range on one chromosome is specified
      par::run_start = getMarkerNumber((*this),par::m1);
      par::run_end = getMarkerNumber((*this),par::m2);
      
      if (par::run_start==-1) error("--from {marker} not found");
      if (par::run_end==-1) error("--to {marker} not found");    

      // Do we require a window around a specific SNP?
      if ( par::run_start == par::run_end && par::window > 0 )
	{
	  vector<int> m = getWindowRange( *this , par::run_start ); 
	  par::run_start = m[0];
	  par::run_end = m[1];
	}

  
      if (getMarkerChromosome((*this),par::m1) != getMarkerChromosome((*this),par::m2))
	{
	  string msg = "--from {marker} and --to {marker} must lie on same chromosome";
	  msg += "\nwhereas these lie on chromosomes "+int2str(getMarkerChromosome((*this),par::m1));	
	  msg += " and "+int2str(getMarkerChromosome((*this),par::m2));
	  error(msg);
	}
    }

  // Get order right
  if (par::run_start > par::run_end)
    {
      int tmp = par::run_start;
      par::run_start = par::run_end;
      par::run_end = tmp;
    }

  int ccode = locus[par::run_start]->chr;
  printLOG("Scan region on chromosome " + int2str(ccode)
	   + " from [ " + locus[par::run_start]->name 
	   + " ] to [ " + locus[par::run_end]->name 
	   + " ]\n");
  
}



void defineHorseChromosomes()
{

  // 31 autosomes + X + Y + etc

  par::chr_haploid.resize(31 + 2 + 1 );
  par::chr_sex.resize(31 + 2 + 1 );
  par::chr_Y.resize(31 + 2 + 1 );
  par::chr_code.resize(31 + 2 + 1 );

  for (int i=0; i<=31; i++)
    {
      par::chr_haploid[i] = par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] = int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

  par::chr_sex[32] = true;
  par::chr_haploid[32] = false;
  par::chr_Y[32] = false;
  par::chr_code[32] = "X";
  par::chr_map.insert( make_pair("X",32) );
  par::chr_map.insert( make_pair("x",32) );
  par::chr_map.insert( make_pair("32",32) );

  par::chr_sex[33] = false;
  par::chr_haploid[33] = true;
  par::chr_Y[33] = true;
  par::chr_code[33] = "Y";
  par::chr_map.insert( make_pair("Y",33) );
  par::chr_map.insert( make_pair("y",33) );
  par::chr_map.insert( make_pair("33",33) );
}


void defineSheepChromosomes()
{

  // 2n = 54 
  // 26 autosomes + X + Y + etc

  par::chr_haploid.resize(26 + 2 + 1 );
  par::chr_sex.resize(26 + 2 + 1 );
  par::chr_Y.resize(26 + 2 + 1 );
  par::chr_code.resize(26 + 2 + 1 );

  for (int i=0; i<=26; i++)
    {
      par::chr_haploid[i] = par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] = int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

  par::chr_sex[27] = true;
  par::chr_haploid[27] = false;
  par::chr_Y[27] = false;
  par::chr_code[27] = "X";
  par::chr_map.insert( make_pair("X",27) );
  par::chr_map.insert( make_pair("x",27) );
  par::chr_map.insert( make_pair("27",27) );

  par::chr_sex[28] = false;
  par::chr_haploid[28] = true;
  par::chr_Y[28] = true;
  par::chr_code[28] = "Y";
  par::chr_map.insert( make_pair("Y",28) );
  par::chr_map.insert( make_pair("y",28) );
  par::chr_map.insert( make_pair("28",28) );
}

void defineRiceChromosomes()
{

  // 12 haploid chromosomes (+ 0 dummy code)

  par::chr_haploid.resize(12 + 1 );
  par::chr_sex.resize(12 + 1);
  par::chr_Y.resize(12 + 1);
  par::chr_code.resize(12 + 1);

  for (int i=0; i<=12; i++)
    {
      par::chr_haploid[i] = true;
      par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] = int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

}

void defineDogChromosomes()
{

  // 38 autosomes + X + Y + XY + 0 missing code

  par::chr_haploid.resize(38 + 3 + 1 );
  par::chr_sex.resize(38 + 3 + 1 );
  par::chr_Y.resize(38 + 3 + 1 );
  par::chr_code.resize(38 + 3 + 1 );

  for (int i=0; i<=38; i++)
    {
      par::chr_haploid[i] = par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] = int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

  par::chr_sex[39] = true;
  par::chr_haploid[39] = false;
  par::chr_Y[39] = false;
  par::chr_code[39] = "X";
  par::chr_map.insert( make_pair("X",39) );
  par::chr_map.insert( make_pair("x",39) );
  par::chr_map.insert( make_pair("39",39) );

  par::chr_sex[40] = false;
  par::chr_haploid[40] = true;
  par::chr_Y[40] = true;
  par::chr_code[40] = "Y";
  par::chr_map.insert( make_pair("Y",40) );
  par::chr_map.insert( make_pair("y",40) );
  par::chr_map.insert( make_pair("40",40) );

  par::chr_sex[41] = false;
  par::chr_haploid[41] = false;
  par::chr_Y[41] = false;
  par::chr_code[41] = "XY";
  par::chr_map.insert( make_pair("XY",41) );
  par::chr_map.insert( make_pair("xy",41) );
  par::chr_map.insert( make_pair("41",41) );
  
}

void defineMouseChromosomes()
{

  // 19 autosomes + X + Y + 0 missing code

  par::chr_haploid.resize(19 + 2 + 1 );
  par::chr_sex.resize(19 + 2 + 1 );
  par::chr_Y.resize(19 + 2 + 1 );
  par::chr_code.resize(19 + 2 + 1 );

  for (int i=0; i<=19; i++)
    {
      par::chr_haploid[i] = par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] = int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

  par::chr_sex[20] = true;
  par::chr_haploid[20] = false;
  par::chr_Y[20] = false;
  par::chr_code[20] = "X";
  par::chr_map.insert( make_pair("X",20) );
  par::chr_map.insert( make_pair("x",20) );
  par::chr_map.insert( make_pair("20",20) );

  par::chr_sex[21] = false;
  par::chr_haploid[21] = true;
  par::chr_Y[21] = true;
  par::chr_code[21] = "Y";
  par::chr_map.insert( make_pair("Y",21) );
  par::chr_map.insert( make_pair("y",21) );
  par::chr_map.insert( make_pair("21",21) );
  
}


void defineCowChromosomes()
{

  // 29 autosomes + X + Y + 0 missing code

  par::chr_haploid.resize(29 + 2 + 1 );
  par::chr_sex.resize(29 + 2 + 1 );
  par::chr_Y.resize(29 + 2 + 1 );
  par::chr_code.resize(29 + 2 + 1 );

  for (int i=0; i<=29; i++)
    {
      par::chr_haploid[i] = par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] =  int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

  par::chr_sex[30] = true;
  par::chr_haploid[30] = false;
  par::chr_Y[30] = false;
  par::chr_code[30] = "X";
  par::chr_map.insert( make_pair("X",30) );  
  par::chr_map.insert( make_pair("x",30) );
  par::chr_map.insert( make_pair("30",30) );
  
  par::chr_sex[31] = false;
  par::chr_haploid[31] = true;
  par::chr_Y[31] = true;
  par::chr_code[31] = "Y";
  par::chr_map.insert( make_pair("Y",31) );
  par::chr_map.insert( make_pair("y",31) );
  par::chr_map.insert( make_pair("31",31) );
}

void defineHumanChromosomes()
{

  // 22 autosomes + X + Y + XY + M + 0 missing code

  par::chr_haploid.resize( 22 + 2 + 2 + 1 );
  par::chr_sex.resize(22 + 2 + 2 + 1 );
  par::chr_Y.resize(22 + 2 + 2 + 1 );
  par::chr_code.resize(22 + 2 + 2 + 1 );
  
  for (int i=0; i<=22; i++)
    {
      par::chr_haploid[i] = par::chr_sex[i] = par::chr_Y[i] = false;
      par::chr_code[i] =  int2str(i);
      par::chr_map.insert( make_pair( int2str(i), i ) );
    }

  // X chromosome
  par::chr_sex[23] = true;
  par::chr_haploid[23] = false;
  par::chr_Y[23] = false;
  par::chr_code[23] = "X";
  par::chr_map.insert( make_pair("X",23) );
  par::chr_map.insert( make_pair("x",23) );
  par::chr_map.insert( make_pair("23",23) );

  // Y chromosome
  par::chr_sex[24] = false;
  par::chr_haploid[24] = true;
  par::chr_Y[24] = true;
  par::chr_code[24] = "Y";
  par::chr_map.insert( make_pair("Y",24) );
  par::chr_map.insert( make_pair("y",24) );
  par::chr_map.insert( make_pair("24",24) );

  // XY chromosome
  par::chr_sex[25] = false;
  par::chr_haploid[25] = false;
  par::chr_Y[25] = false;
  par::chr_code[25] = "XY";
  par::chr_map.insert( make_pair("XY",25) );
  par::chr_map.insert( make_pair("xy",25) );
  par::chr_map.insert( make_pair("25",25) );

  // MT chromosome
  par::chr_sex[26] = false;
  par::chr_haploid[26] = true;
  par::chr_Y[26] = false;  
  par::chr_code[26] = "MT";
  par::chr_map.insert( make_pair("MT",26) );
  par::chr_map.insert( make_pair("mt",26) );
  par::chr_map.insert( make_pair("M",26) );
  par::chr_map.insert( make_pair("m",26) );
  par::chr_map.insert( make_pair("26",26) );
}

void sizeMatrix(matrix_t & m , int r, int c)
{
  m.clear();
  m.resize(r);
  for (int i=0; i<r; i++)
    m[i].resize(c,0);
}

void sizefMatrix(fmatrix_t & m , int r, int c)
{
  m.clear();
  m.resize(r);
  for (int i=0; i<r; i++)
    m[i].resize(c,0);
}

void sizeTable(table_t & m , int r, int c)
{
  m.clear();
  m.resize(r);
  for (int i=0; i<r; i++)
    m[i].resize(c,0);
}



/*
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
  
  if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;
  
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) 
    {
      error("Internal error: negative count in HWE test: "
	    +int2str(obs_hets)+" "
	    +int2str(obs_hom1)+" "
	    +int2str(obs_hom2));
    }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes   = obs_hets + obs_homc + obs_homr;

  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if (het_probs == NULL) 
    error("Internal error: SNP-HWE: Unable to allocate array" );
  
  int i;
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] = 0.0;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  
  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;
  
  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
    {
      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
      sum += het_probs[curr_hets - 2];

      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
      curr_homr++;
      curr_homc++;
    }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
    {
      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
	/((curr_hets + 2.0) * (curr_hets + 1.0));
      sum += het_probs[curr_hets + 2];
      
      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
      curr_homr--;
      curr_homc--;
    }
  
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
   double p_hi = het_probs[obs_hets];
   for (i = obs_hets + 1; i <= rare_copies; i++)
     p_hi += het_probs[i];
   
   double p_lo = het_probs[obs_hets];
   for (i = obs_hets - 1; i >= 0; i--)
      p_lo += het_probs[i];
   
   double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
  */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++)
    {
      if (het_probs[i] > het_probs[obs_hets])
	continue;
      p_hwe += het_probs[i];
    }
   
  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}




// Convert dataset from Individual-major format to SNP-major

void Plink::Ind2SNP()
{
  printLOG("Converting data to SNP-major format\n");

  SNP.clear();
  
  vector<Individual*>::iterator person = sample.begin();
  
  // Initialise SNP positions per person
  while ( person != sample.end() )
    {
      (*person)->i1 = (*person)->one.end()-1;
      (*person)->i2 = (*person)->two.end()-1;
      person++;
    }
  
  // Copy, per SNP
  int l = 0;
  while ( l < nl_all  )
    {
 
      CSNP * newlocus = new CSNP;
      
      person = sample.begin();
       
      while ( person != sample.end() )
	{
	  // Add genotype to SNP-major storage
	  newlocus->one.push_back( *((*person)->i1) );
	  newlocus->two.push_back( *((*person)->i2) );
	  
	  // Shift one SNP back
	  (*person)->i1--;
	  (*person)->i2--;

	  // Remove individual-major storage
	  (*person)->one.pop_back();
	  (*person)->two.pop_back();
	  
	  // Advance to next person
	  person++;
	}
            
      // And add this new SNP to the main list
      SNP.push_back(newlocus);
      
      // Next SNP
      l++;
    }
 
  // We finally need to reverse the order of these 
  reverse(SNP.begin(), SNP.end());

  par::SNP_major = true;
}

// Convert dataset from SNP-major format to Individual-major
void Plink::SNP2Ind()
{
  
  printLOG("Converting data to Individual-major format\n");
    
  vector<Individual*>::iterator person = sample.begin();
  
  // Make sure these containers are empty
  while ( person != sample.end() )
    {
      (*person)->one.clear();
      (*person)->two.clear();
      person++;
    }
  
  ///////////////////////////////
  // Iterate over SNPs
  
  vector<CSNP*>::iterator s = SNP.begin();
 
  while ( s != SNP.end() )
    {
            
      /////////////////////////////
      // Iterate over individuals
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator person = sample.begin();
      
      while ( person != sample.end() )
	{
	  
	  // Add SNP alleles
	  (*person)->one.push_back(*i1);
	  (*person)->two.push_back(*i2);
	  
 	  // Shift one SNP back
 	  i1++;
 	  i2++;
	  	  
 	  // Advance to next person
 	  person++;
 	}
      
      // For this SNP, remove SNP-major storage completely
      delete (*s);
      
      // Next SNP
      s++;      
    }
  
  SNP.clear();
  par::SNP_major = false;
}


int Plink::deleteSNPs(set<string> & mset)
{
  vector<bool> b(nl_all);
  for (int l=0; l<nl_all; l++)
    if ( mset.find( locus[l]->name ) != mset.end() )
      b[l] = true;
    else
      b[l] = false;
  return deleteSNPs(b);
}

int Plink::deleteSNPs(set<Locus*> & mset)
{
  vector<bool> b(nl_all);
  for (int l=0; l<nl_all; l++)
    if ( mset.find( locus[l] ) != mset.end() )
      b[l] = true;
    else
      b[l] = false;
  return deleteSNPs(b);
}

int Plink::deleteIndividuals(set<Individual*> & pset)
{
  vector<bool> b(n);
  for (int i=0; i<n; i++)
    if ( pset.find( sample[i] ) != pset.end() )
      b[i] = true;
    else
      b[i] = false;
  return deleteIndividuals(b);
}



int Plink::keepSNPs(set<string> & mset)
{
  vector<bool> b(nl_all);
  for (int l=0; l<nl_all; l++)
    if ( mset.find( locus[l]->name ) != mset.end() )
      b[l] = false;
    else
      b[l] = true;
  return deleteSNPs(b);
}

int Plink::keepSNPs(set<Locus*> & mset)
{
  vector<bool> b(nl_all);
  for (int l=0; l<nl_all; l++)
    if ( mset.find( locus[l] ) != mset.end() )
      b[l] = false;
    else
      b[l] = true;
  return deleteSNPs(b);
}

int Plink::keepIndividuals(set<Individual*> & pset)
{
  vector<bool> b(n);
  for (int i=0; i<n; i++)
    if ( pset.find( sample[i] ) != pset.end() )
      b[i] = false;
    else
      b[i] = true;
  return deleteIndividuals(b);
}



int Plink::deleteSNPs(vector<bool> & del)
{
  
   
   // Remove SNPs that have T in vector 'del'
   // We expect this vector to be same length
   // current number of SNPs
   
   int original = SNP.size();


   //////////////////////////////////////////////
   // SNP-major mode

   if (par::SNP_major)
     {
       
       vector<CSNP*>::iterator s1 = SNP.begin();
       vector<CSNP*>::iterator s2 = SNP.begin();
       
       vector<bool>::iterator d = del.begin();
       int i = 0;
       while ( s1 != SNP.end() ) 
	 {
	   //cout << "i = " << i++ << "\n";
	   // Keep this SNP?
	   if ( ! *d ) 
	     {
	       *s2 = *s1;
	       s2++;
	     }
	   else
	     {
	       delete (*s1);
	     }

	   s1++;
	   d++;
	 }
       
       // Then erase the remaining SNPs
       SNP.erase(s2,SNP.end());
       
     }
  else
    {
      // Individual major-mode
      // copy remaining SNPs to a new list
  
      vector<Individual*>::iterator person = sample.begin();
      while ( person != sample.end() )
	{
	  
	  // Set individual iterator at start of SNP list
	  vector<bool>::iterator one1 = (*person)->one.begin();
	  vector<bool>::iterator two1 = (*person)->two.begin();

	  vector<bool>::iterator one2 = (*person)->one.begin();
	  vector<bool>::iterator two2 = (*person)->two.begin();

	  vector<bool>::iterator d = del.begin();

	  
	  while ( one1 != (*person)->one.end() )
	    {
	      
	      // Keep this marker
	      if ( ! * d ) 
		{  
		  *one2 = *one1;
		  *two2 = *two1;
		  
		  // Advance next saved SNP
		  one2++;
		  two2++;
		}

	      // Advance to next to-be-checked SNP
	      one1++;
	      two1++;
	      d++;
	    }
	  
	  // Then erase the remaining SNPs
	  (*person)->one.erase(one2,(*person)->one.end());
	  (*person)->two.erase(two2,(*person)->two.end());
	  
	  // Next individual
	  person++;	  
	}
    }

   
   ///////////////////////////////////////////////////
   // Second, remove deleted SNPs from the locus list

   vector<Locus*>::iterator loc1 = locus.begin();
   vector<Locus*>::iterator loc2 = locus.begin();
   vector<bool>::iterator d = del.begin();

   while ( loc1 != locus.end() )
     {    
       
       // Should we keep this SNP?
       if ( ! *d )
	 {
	   *loc2 = *loc1;
	   loc2++;
	 }
       else
	 delete (*loc1);

       loc1++;
       d++;
       
     }
       

   // Then erase the remaining SNPs
   // and the storage
   locus.erase(loc2,locus.end());

   // Keep track of the number of SNPs remaining
   nl_all = locus.size();
   
   // Return the number of SNPs we chucked
   return original - nl_all;

}




int Plink::deleteIndividuals(vector<bool> & del)
{


  // Remove SNPs that have T in vector 'del'
  // We expect this vector to be same length
  // current number of SNPs
  
  int original = sample.size();
  
  //////////////////////////////////////////////
  // SNP-major mode
  
  if (par::SNP_major)
    {
      
      // Erase genotype data (SNP-major order)
      // Consider each SNP in outer loop
      vector<CSNP*>::iterator s = SNP.begin();
      while ( s != SNP.end() )
	{ 
	  
	  vector<bool>::iterator one1 = (*s)->one.begin();
	  vector<bool>::iterator two1 = (*s)->two.begin();
	  
	  vector<bool>::iterator one2 = (*s)->one.begin();
	  vector<bool>::iterator two2 = (*s)->two.begin();
	  
	  vector<bool>::iterator d = del.begin();

	  // Consider each person
	  while ( one1 != (*s)->one.end() )
	    {
	      
	      // Keep this individual?
	      if ( ! *d ) 
		{
		  (*one2) = (*one1);
		  (*two2) = (*two1);

		  one2++;
		  two2++;
		}
	      
	      one1++;
	      two1++;
	      d++;
	    }
	 	  
	  // Erase old storage
	  (*s)->one.erase(one2,(*s)->one.end());
	  (*s)->two.erase(two2,(*s)->two.end());

	  // Next SNP
	  s++;
	}
    }

  /////////////////////////////////////////////
  // Whether SNP-major or individual-major
  // we still need to take care of the sample[]
  
  vector<Individual*>::iterator person1 = sample.begin();
  vector<Individual*>::iterator person2 = sample.begin();
  
  vector<bool>::iterator d = del.begin();
  
  while ( person1 != sample.end() )
    {
      // Keep this person
      if ( ! * d ) 
	{  
	  *person2 = *person1;
	  
	  // Advance next saved SNP
	  person2++;
	}
      else
	{
	  // Free storage
	  delete (*person1);	  
	}
      
      // Advance to next to-be-checked SNP
      person1++;
      d++;
    }
  
  // Erase pointers at end of vector
  sample.erase(person2,sample.end());
  
  // Adjust sample statistics
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);  
  
  return original - n; 
}






void Plink::filterOnCovariate()
{

  printLOG("Filtering individuals based on [ "+par::filter_filename+" ]\n");
  printLOG("Filtering criterion is [ " + par::filter_value 
	   + " ] for cluster " + int2str(par::mult_filter) +"\n");

  // Expand filter criteria to allow a list 
  string tmp = par::range_delimiter;
  par::range_delimiter = " ";
  NList tlist(0);
  vector<string> filters = tlist.deparseStringList( par::filter_value );  
  par::range_delimiter = tmp;
  

  // Swap q-match filename as covariate file

  string tmp_covar_file = par::include_cluster_filename;
  int tmp_mult_covar = par::mult_clst;

  par::include_cluster_filename = par::filter_filename;
  par::mult_clst = par::mult_filter;
  
  if (!readClusterFile())
    error("Problem reading filter file [ " + par::filter_filename + " ]\n");
  
  // Put back the original covariate specificiation  
  par::include_cluster_filename = tmp_covar_file;
  par::mult_clst = tmp_mult_covar;
  

  // Screen-based on covariate
  // vector to record which individuals to be deleted

  vector<bool> indel(sample.size(),false);
  int n_removed1=0;
  for (int i=0; i<sample.size(); i++)
    {

      bool removeThisSample = true;

      int thisK = sample[i]->sol;
      
      for (int j=0; j<filters.size(); j++)
	{	     
	  map<string,int>::iterator k = kmap.find( filters[j] );
	  
	  if ( k != kmap.end() && k->second == thisK )
	    {
	      removeThisSample = false;
	      break;
	    }
	}
      
      if ( removeThisSample )
	{
	  indel[i] = true;
	  n_removed1++;
	}
    }
  
  // And now remove these individuals, so that 
  // SNP-based statistics are calculated with 
  // these samples already excluded
  
  int n_removed = deleteIndividuals(indel);
  
  if (n_removed != n_removed1)
    error("Internal problem in filterOnCovariate, please contact SMP\n");
  printLOG(int2str(n_removed)+" individuals removed based on filter\n");

  // Remove these as clusters now
  for (int i=0; i<sample.size(); i++)
    {
      sample[i]->sol = 0;
    }
  
  nk=1;
  kname.resize(0);
  kmap.clear();
  
}



void Plink::filterOnCase()
{
  printLOG("Filtering cases only ... ");

  // Remove controls, missing
  vector<bool> indel(sample.size(),false);
  for (int i=0; i<sample.size(); i++)
    if ( (!sample[i]->aff) || sample[i]->missing ) 
      indel[i] = true;

  int n_original = sample.size();
  int n_removed = deleteIndividuals(indel);
  printLOG(int2str(n_original-n_removed)+" individuals remaining\n");
  
  // Reset number of individuals
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);
  
}

void Plink::filterOnControl()
{
  printLOG("Filtering controls only ... ");

  // Remove cases, missing
  vector<bool> indel(sample.size(),false);
  for (int i=0; i<sample.size(); i++)
    if ( sample[i]->aff || sample[i]->missing ) 
      indel[i] = true;

  int n_original = sample.size();  
  int n_removed = deleteIndividuals(indel);
  printLOG(int2str(n_original-n_removed)+" individuals remaining\n");
  
  // Reset number of individuals
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);
  
}

void Plink::filterOnMale()
{
  printLOG("Filtering males only ... ");
  
  vector<bool> indel(sample.size(),false);
  for (int i=0; i<sample.size(); i++)
    if ( sample[i]->sexcode != "1" )
      indel[i] = true;

  int n_original = sample.size();  
  int n_removed = deleteIndividuals(indel);
  printLOG(int2str(n_original-n_removed)+" individuals remaining\n");
  
  // Reset number of individuals
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);
  
}

void Plink::filterOnFemale()
{
  printLOG("Filtering females only ... ");
  
  vector<bool> indel(sample.size(),false);
  for (int i=0; i<sample.size(); i++)
    if ( sample[i]->sexcode != "2" )
      indel[i] = true;
  
  int n_original = sample.size();
  int n_removed = deleteIndividuals(indel);
  printLOG(int2str(n_original-n_removed)+" individuals remaining\n");
  
  // Reset number of individuals
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);
  
}

void Plink::filterOnFounder()
{
  printLOG("Filtering founders only ... ");
  
  vector<bool> indel(sample.size(),false);
  for (int i=0; i<sample.size(); i++)
    if ( ! sample[i]->founder )
      indel[i] = true;
  
  int n_original = sample.size();
  int n_removed = deleteIndividuals(indel);
  printLOG(int2str(n_original-n_removed)+" individuals remaining\n");
  
  // Reset number of individuals
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);
  
}

void Plink::filterOnNonFounder()
{
  printLOG("Filtering nonfounders only ... ");
  
  vector<bool> indel(sample.size(),false);
  for (int i=0; i<sample.size(); i++)
    if ( sample[i]->founder )
      indel[i] = true;
  
  int n_original = sample.size();
  int n_removed = deleteIndividuals(indel);
  printLOG(int2str(n_original-n_removed)+" individuals remaining\n");
  
  // Reset number of individuals
  n = sample.size();
  np = (int)((double)(n*(n-1))/(double)2);
  
}


void Plink::attribFilterSNP()
{
  string tmp = par::range_delimiter;
  par::range_delimiter = " ";

  printLOG("Filtering markers, from [ " + par::snp_attrib_file + " ] "
	   + "criterion: " + par::snp_attrib_value + "\n");
  checkFileExists( par::snp_attrib_file );  
  ifstream IN1( par::snp_attrib_file.c_str() , ios::in );

  NList nl(0);
  vector<string> vlist = nl.deparseStringList( par::snp_attrib_value );
  map<string,bool> vset;
  
  bool posMatch = false;
  bool negMatch = false;

  for (int i=0; i<vlist.size(); i++)
    {
      if ( vlist[i].substr(0,1) == "-" )
	{
	  vset.insert( make_pair( vlist[i].substr(1), false ) );
	  negMatch = true;
	}
      else
	{
	  vset.insert( make_pair( vlist[i], true ) );
	  posMatch = true;
	}
    }

  set<Locus*> mset;
  map<string,Locus*> mlocus;
  for (int l=0; l<nl_all; l++)
    mlocus.insert(make_pair(locus[l]->name,locus[l]));
  
  while ( ! IN1.eof() )
    {
      vector<string> tok = tokenizeLine( IN1 );
      if ( tok.size() == 0 )
	continue;

      map<string,Locus*>::iterator i = mlocus.find( tok[0] );
      if ( i == mlocus.end() )
	continue;
      
      bool match = false; // T if at least 1 positive match
      bool exclude = false; // T if at least 1 negative match
      
      // Logical OR for matching
      for (int j=1; j<tok.size(); j++)
	{
	  
	  // Is this token relevant
	  map<string,bool>::iterator k = vset.find( tok[j] );

	  // No
	  if ( k == vset.end() )
	    continue;

	  // Yes
	  if ( k->second ) 
	    match = true;
	  else
	    exclude = true;
	}

      
      // Keep this SNP?
      if ( ( match || (!posMatch) ) &&
	   ( (!exclude) || (!negMatch) ) )  
	mset.insert( i->second );

    }
  IN1.close();
  int rem = keepSNPs( mset );
  printLOG("Removed " + int2str(rem) + " SNPs based on this\n");

  par::range_delimiter = tmp;
  return;

}



void Plink::attribFilterInd()
{
  string tmp = par::range_delimiter;
  par::range_delimiter = " ";

  printLOG("Filtering individuals, from [ " + par::ind_attrib_file + " ] "
	   + "criterion: " + par::ind_attrib_value + "\n");
  checkFileExists( par::ind_attrib_file );  
  ifstream IN1( par::ind_attrib_file.c_str() , ios::in );
  
  NList nl(0);
  vector<string> vlist = nl.deparseStringList( par::ind_attrib_value );
  map<string,bool> vset;
  bool posMatch = false;
  bool negMatch = false;

  for (int i=0; i<vlist.size(); i++)
    {
      if ( vlist[i].substr(0,1) == "-" )
	{
	  vset.insert( make_pair( vlist[i].substr(1), false ) );
	  negMatch = true;
	}
      else
	{
	  vset.insert( make_pair( vlist[i], true ) );
	  posMatch = true;
	}
    }

  set<Individual*> mset;
  map<string,Individual*> mpeople;
  for (int i=0; i<n; i++)
    mpeople.insert(make_pair(sample[i]->fid + "_" + sample[i]->iid,sample[i]));
  
  while ( ! IN1.eof() )
    {
      vector<string> tok = tokenizeLine( IN1 );
      if ( tok.size() < 2  )
	continue;

      map<string,Individual*>::iterator i = 
	mpeople.find( tok[0] + "_" + tok[1] );

      if ( i == mpeople.end() )
	continue;

      bool match = false; // T if at least 1 positive match
      bool exclude = false; // T if at least 1 negative match

      for (int j=1; j<tok.size(); j++)
	{
	  map<string,bool>::iterator k = vset.find( tok[j] );
	  
	  if ( k == vset.end() )
	    continue;
	  if ( k->second ) 
	    match = true;
	  else
	    exclude = true;
	}

      // Keep these people
      if ( ( match || (!posMatch) ) &&
	   ( (!exclude) || (!negMatch) ) )  
	mset.insert( i->second );
    }
  IN1.close();
  int rem = keepIndividuals( mset );
  printLOG("Removed " + int2str(rem) + " individuals based on this\n");

  par::range_delimiter = tmp;
  return;
}


void Plink::dummyLoader()
{

  // Create dummy dataset full of heterozygotes

  int L = par::dummy_nsnp;
  int N = par::dummy_nind;

  for (int l=0;l<L;l++)
    {
      Locus * loc = new Locus;
      loc->name = "snp"+int2str(l);
      loc->chr = 1;
      loc->allele1 = "A";
      loc->allele2 = "B";
      loc->bp = l ; 
      loc->pos = 0;

      locus.push_back(loc);

      CSNP * newset = new CSNP;

      newset->one.resize(N);
      newset->two.resize(N);

      for ( int i = 0 ; i < N ; i++ ) 
	{
	  
	  int g = 0;
	  if (CRandom::rand() > 0.5)
	    g++;
	  if (CRandom::rand() > 0.5)
	    g++;
	  
	  if ( g == 0 ) 
	    {
	      newset->one[i] = false;
	      newset->two[i] = false;
	    }
	  else if ( g == 1 ) 
	    {
	      newset->one[i] = false;
	      newset->two[i] = true;
	    }
	  else
	    {
	      newset->one[i] = true;
	      newset->two[i] = true;
	    }
	}

      SNP.push_back(newset);
    }
  
  for (int i=0;i<N;i++)
    {
      Individual * person = new Individual;
      person->fid = person->iid = "per"+int2str(i);
      person->missing = false;
      person->pat = "0";
      person->mat = "0";

      if ( CRandom::rand() > 0.5 ) 
	person->phenotype = 1;
      else
	person->phenotype = 2;

      person->sex = false;
      person->sexcode = "2";
      sample.push_back(person);      
    }

}


void Plink::alleleRecoding()
{
  
  vector<Locus*>::iterator loc = locus.begin();

  while( loc != locus.end() )
  {
  
    if ( par::recode_1234 )
	{
        if ( (*loc)->allele1 == "A" || (*loc)->allele1 == "a" )
             (*loc)->allele1 = "1";
        if ( (*loc)->allele1 == "C" || (*loc)->allele1 == "c" )
             (*loc)->allele1 = "2";
        if ( (*loc)->allele1 == "G" || (*loc)->allele1 == "g" )
             (*loc)->allele1 = "3";
        if ( (*loc)->allele1 == "T" || (*loc)->allele1 == "t" )
             (*loc)->allele1 = "4";

        if ( (*loc)->allele2 == "A" || (*loc)->allele1 == "a" )
             (*loc)->allele2 = "1";
        if ( (*loc)->allele2 == "C" || (*loc)->allele1 == "c" )
             (*loc)->allele2 = "2";
        if ( (*loc)->allele2 == "G" || (*loc)->allele1 == "g" )
             (*loc)->allele2 = "3";
        if ( (*loc)->allele2 == "T" || (*loc)->allele1 == "t" )
             (*loc)->allele2 = "4";

      }
     else if ( par::recode_ACGT )
      {
        if ( (*loc)->allele1 == "1" )
             (*loc)->allele1 = "A";
        if ( (*loc)->allele1 == "2" )
             (*loc)->allele1 = "C";
        if ( (*loc)->allele1 == "3" )
             (*loc)->allele1 = "G";
        if ( (*loc)->allele1 == "4" )
             (*loc)->allele1 = "T";

        if ( (*loc)->allele2 == "1" )
             (*loc)->allele2 = "A";
        if ( (*loc)->allele2 == "2" )
             (*loc)->allele2 = "C";
        if ( (*loc)->allele2 == "3" )
             (*loc)->allele2 = "G";
        if ( (*loc)->allele2 == "4" )
             (*loc)->allele2 = "T";
      }

    loc++;
  }

}

vector<string> commaParse(string s)
{
  NList tlist(0);
  tlist.setRangeChar(" ");
  return tlist.deparseStringList( s );
}

string searchAndReplace(string str, 
			string searchString, 
			string replaceString)
{
  string::size_type pos = 0;
  while ( (pos = str.find(searchString, pos)) != string::npos ) 
    {
      str.replace( pos, searchString.size(), replaceString );
      pos++;
    }
  return str;
}


void makePersonMap(Plink &P, map<string,Individual*> & uid)
{
  for (int i=0; i<P.sample.size(); i++)
    uid.insert(make_pair(P.sample[i]->fid+"_"+P.sample[i]->iid,P.sample[i]));
}

void makeLocusMap(Plink &P, map<string,int> & mlocus)
{
  for (int l=0; l<P.nl_all; l++)
    mlocus.insert(make_pair(P.locus[l]->name,l));
}



void smoother(Plink & P, 
	      vector_t & input,
	      int n, 
	      vector_t & output1,
	      vector_t & output2,
	      vector<int> & count)
{

  // Take a vector a numbers, 0..nl_all, and 
  // smooth, respecting chromosome boundaries, 
  // based one par::seg_window_kb and par::seg_window_step
  
  if ( input.size() != P.nl_all )
    error("Problem in smoother()\n");

  
  output1.resize( P.nl_all );
  output2.resize( P.nl_all );
  count.resize( P.nl_all );


  for (int l=0; l<P.nl_all; l++)
    {

      double x1 = input[l];
      double c1 = n - x1;
  
      int involved = 1;
      
      // pick this SNP, move forwards and backwards

      Locus * loc1 = P.locus[l];
      
      
      int l2 = l;
      while ( true ) 
	{
	  l2--;
	  if ( l2 < 0 ) break;

	  Locus * loc2 = P.locus[l2];
	  
	  if ( loc2->chr != loc1->chr ) break;
	  if ( loc1->bp - loc2->bp > par::seg_test_window_bp  ) break;
	  x1 += input[l2];
	  c1 += n - input[l2];
	  involved++;
	}
      
      l2 = l;
      while ( true ) 
	{
	  l2++;
	  if ( l2 == P.nl_all ) break;

	  Locus * loc2 = P.locus[l2];

	  if ( loc2->chr != loc1->chr ) break;
	  if ( loc2->bp - loc1->bp > par::seg_test_window_bp ) break;
	  x1 += input[l2];
	  c1 += n - input[l2];
	  involved++;
	}
      
      output1[l] = x1;
      output2[l] = c1;
      count[l] = involved;
    }

  return;
}


map<string, set<Range> > readRange(string filename)
{
  

  // Format: CHR BP1 BP2 (NAME) (GROUP)

  // If same named range read twice, take largest range
  
  checkFileExists(filename);    
  ifstream IN(filename.c_str(),ios::in);
  IN.clear();
  PP->printLOG("Reading list of ranges from [ " + filename + " ]\n");
  PP->printLOG("Allowing a " + int2str( par::make_set_border/1000 ) 
	       + " kb window around each range\n");

  map<string, set<Range> > ranges;
  int rcount = 0;

  // Track number of ranges
  while (!IN.eof())
    {
      Range r;      
      char cline[par::MAX_LINE_LENGTH];
      IN.getline(cline,par::MAX_LINE_LENGTH,'\n');
      string sline = cline;
      if (sline=="") continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);

      if ( tokens.size() < 4 ) 
	error("Problem with line:\n" + sline );

      string chr = tokens[0];

      if ( ! from_string<int>( r.start  , tokens[1], std::dec ) )
	error("Problem with position : " + tokens[1] );

      if ( ! from_string<int>( r.stop  , tokens[2], std::dec ) )
	error("Problem with position : " + tokens[2] );
      
      // Add any specified border region 

      r.start -= par::make_set_border;
      r.stop += par::make_set_border;
      if (r.start < 0 ) r.start = 0;
      
      if ( chr == "" ) continue;
      if ( r.start > r.stop ) continue;
      r.chr = getChromosomeCode(chr);

      // Assign a name for this range/set
      r.name = tokens[3];
      
      bool hasGroupLabel =  tokens.size() >= 5;
	
      // Assign a group?
      if ( par::make_set_ignore_group || ! hasGroupLabel )
	r.group = -1;
      else
	{
	  map<string,int>::iterator i = Range::groupNames.find( tokens[4] );
	  if ( i != Range::groupNames.end() )
	    r.group = i->second;
	  else
	    {
	      int t = Range::groupNames.size();
	      Range::groupNames.insert( make_pair( tokens[4] , t ) );
	      r.group = t;
	    }
	}

      // Have we already seen this range, or is it new?
      
      set<Range> * s;

      string fullName = r.name;
      if ( hasGroupLabel )
	fullName += "_" + int2str(r.group);
      
      if ( ranges.find( fullName ) == ranges.end() )
	{
	  // Never seen this set before
	  set<Range> tmp;
	  tmp.insert(r);
	  ranges.insert(make_pair( fullName, tmp ));
	}
      else
	{
	  map<string, set<Range> >::iterator ri = ranges.find( fullName );

	  // Add this new range to the existing set
	  ri->second.insert( r );
	}
      
      ++rcount;
    }
  
  IN.close();
  
  PP->printLOG("Added " + int2str( rcount )
	       + " distinct ranges to " 
	       + int2str(ranges.size() ) 
	       + " sets\n");

  return ranges;
}


double modelComparisonPValue(Model * alternate, Model * null)
{
  
  if ( par::bt )
    {
      LogisticModel * lalternate = (LogisticModel*)alternate;
      LogisticModel * lnull = (LogisticModel*)null;
      
      return chiprobP( lnull->getLnLk() - lalternate->getLnLk() , 
			lalternate->getNP() - lnull->getNP() ) ;
      
    }
  else
    {
      
      LinearModel * lalternate = (LinearModel*)alternate;
      LinearModel * lnull = (LinearModel*)null;
      
      double F = lalternate->calculateFTest(lnull);
      
      return pF( F, 
		 alternate->getNP() - null->getNP(),
		 alternate->Ysize() - alternate->getNP() - 1 ) ;
    }

  return -1;
}


set<Range*> rangeIntersect(Range & r1, map<string, set<Range> > & ranges)
{
  
  // Return all ranges from 'ranges' that intersect with 'r1'  
  set<Range*> intersected;
  int bp1 = r1.start;
  int bp2 = r1.stop;
  int chr = r1.chr;
  
  map<string, set<Range> >::iterator si = ranges.begin();
  
  while ( si != ranges.end() )
    {
      // string rname = si->first;
      
      set<Range>::iterator ri = si->second.begin();
      while ( ri != si->second.end() )
	{
	  if ( ri->chr != chr )
	    {
	      ++ri;
	      continue;
	    }
	  
	  if ( ri->start <= bp2 && 
	       ri->stop >= bp1 )
	    {
	      intersected.insert( (Range*)&(*ri) );
	      // Only enter this name once
	      break;
	    }
	  ++ri;
	}
      ++si;
    }
  
  return intersected;
}

set<Range*> mapRanges2SNP(int l, map<string, set<Range> > & ranges)
{
    // For a single SNP, return a list of pointers to all 
    // ranges that span it
    
    set<Range*> r;
    // { TODO } 
    return r;
}

int2 mapSNPs2Range(Plink & P, const Range * range)
{

  // For a single range, return the start and stop SNPs that fall within
  // that range
    
  int chr = range->chr;
  
  // Cannot map chromosome
  if ( P.scaffold.find( chr ) == P.scaffold.end() )
    return int2(-1,-1);
    
  int lstart = P.scaffold[chr].lstart;
  int lstop = P.scaffold[chr].lstop;
  int bpstart = P.scaffold[chr].bpstart;
  int bpstop = P.scaffold[chr].bpstop;
  
  // Assume a roughly uniform SNP spacing, to get 
  // good gueses at where to start looking for this
  // range
  
  double prop_start = ( range->start - bpstart ) / (double)( bpstop - bpstart );
  double prop_stop = ( range->stop - bpstart ) / (double)( bpstop - bpstart );

  int guess_start = int ( lstart + prop_start * ( lstop-lstart ) );
  int guess_stop = int ( lstart + prop_stop * ( lstop-lstart ) );
    
  if ( guess_start < lstart ) 
    guess_start = lstart;
  if ( guess_stop < lstart ) 
    guess_stop = lstart;
  
  if ( guess_start > lstop ) 
    guess_start = lstop;
  if ( guess_stop > lstop ) 
    guess_stop = lstop;
  
  ////////////////////
  // Adjust start 
  
  while (1)
    {
      
      if ( P.locus[guess_start]->bp == range->start )
	break;
      
      if ( P.locus[guess_start]->bp > range->start ) 
	{
	  if ( guess_start == lstart || 
	       P.locus[guess_start-1]->bp < range->start )
	    break;
	  --guess_start;
	}
      else 
	{
	  if ( guess_start == lstop )
	    break;
	  
	  // else next SNP right must be first in range
	  if ( P.locus[guess_start+1]->bp >= range->start )
	    {
	      ++guess_start;
	      break;   
	    }
	  ++guess_start;	  
	}
    }
  
  
  ////////////////////
  // Adjust stop 
  
  while (1)
    {

      if ( P.locus[guess_stop]->bp == range->stop ) 
	break;

      if ( P.locus[guess_stop]->bp > range->stop ) 
	{

	  if ( guess_stop == lstart )
	    break;
	  
	  if ( P.locus[guess_stop-1]->bp <= range->stop )
	    {
	      --guess_stop;
	      break;
	    }

	  --guess_stop;
	}
      else 
	{
	  if ( guess_stop == lstop || 
	       P.locus[guess_stop+1]->bp > range->stop )
	    break;   
	  ++guess_stop;
	  
	}
    }

  

  if ( P.locus[guess_start]->bp > range->stop ||
       P.locus[guess_stop]->bp < range->start )
    return int2(-1,-1);

  return int2(guess_start,guess_stop);
}

    
void makeScaffold(Plink & P)
{
    int last = -1;
    
    P.scaffold.clear();

    int thisChromosome = P.locus[0]->chr;
    int nextChromosome;
    int lastSNP = P.nl_all-1;

    for (int l=0; l< P.nl_all; l++)
      {
	
	int chr = P.locus[l]->chr;
	
	// Have we seen this chromosome before? If not, 
	// add this to the list.
	
	if ( P.scaffold.find( chr ) == P.scaffold.end() )
	  {
	    CInfo ci;
	    ci.lstart = l;
	    ci.bpstart = P.locus[l]->bp;	    
	    P.scaffold.insert( make_pair( chr, ci ));	    
	  }
	
	// Is this the end of this chromosome? If so, 
	// also record that.

	if ( l == lastSNP || 
	     chr != P.locus[l+1]->chr )
	  {
	    map<int,CInfo>::iterator i = P.scaffold.find( chr );
	    (i->second).lstop = l;
	    (i->second).bpstop = P.locus[l]->bp;
	  }
	
      }
    
}



  
void mapRangesToSNPs(string filename,
		     map<string, set<Range> > & ranges, 		     
		     map<int,set<Range*> > & snp2range )
{
  
  // Read list of ranges
  ranges = readRange(filename);
  
  // Consider each range
  map<string, set<Range> >::iterator r = ranges.begin();
  int rc = 0;
  

  while ( r != ranges.end() )
    {
      
      set<Range> * theseRanges = &( r->second );
      // Assume a single, unqiue named range (and so set-size=1)
      
      set<Range>::iterator thisRangeSet = theseRanges->begin();	  

      const Range * thisRange = &(*thisRangeSet);
      
      // Assign SNPs
      int2 srange = mapSNPs2Range(*PP, thisRange);
      
      if ( srange.p1 != -1 ) 
	{
	  for ( int l= srange.p1; l<= srange.p2; l++)
	    {
	      if ( snp2range.find(l) == snp2range.end() )
		{
		  set<Range*> t;
		  t.insert( (Range*)thisRange );
		  snp2range.insert(make_pair( l, t ) );
		}
	      else
		{
		  set<Range*> * rp = &(snp2range.find(l)->second);
		  rp->insert((Range*)thisRange);
		}
	    }
	  ++rc;
	}
     
      // Next range/gene
      ++r;
      
    }
  
  // Save the range labels in the correct order

//   r = ranges.begin();
//   while ( r != ranges.end() )
//   {
//     set<Range> * theseRanges = &( r->second );
//     // Assume a single, unqiue named range (and so set-size=1)
//     set<Range>::iterator thisRangeSet = theseRanges->begin();
//     const Range * thisRange = &(*thisRangeSet);
//     rangeLabels.push_back( thisRange->name );
//     ++r;
//   }


  
  PP->printLOG(int2str(rc)+" ranges with at least 1 marker\n");
  PP->printLOG("Assigned " 
	      + int2str( snp2range.size() ) 
	      + " SNPs to at least 1 range\n");
  
}



void Plink::setFlagToCase()
{
  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      if ( person->missing )
	person->flag = false;
      else if ( person->aff )
	person->flag = true;
      else
	person->flag = false;
    }
}

void Plink::setFlagToControl()
{
  for (int i=0;i<n;i++)
    {
      Individual * person = sample[i];
      if ( person->missing )
	person->flag = false;
      else if ( person->aff )
	person->flag = false;
      else
	person->flag = true;
    }
}

string relType(Individual *a , Individual * b)
{
  // UN unrelated
  // FS sibling
  // HS half-sibling
  // PO parent-offspring
  // OT other

  if ( a->fid != b->fid )
    return "UN";
  
  if ( ! ( a->founder || b->founder ) )
    {
      if ( a->pat == b->pat && a->mat == b->mat ) 
	return "FS";
      else if ( a->pat == b->pat || a->mat == b->mat ) 
	return "HS";
    }

  if ( a->pat == b->iid || a->mat == b->iid ||
       b->pat == a->iid || b->mat == a->iid )
    return "PO";

  return "OT";

}

void Plink::outputPermedPhenotypes(Perm & perm)
{


  // Dummy test statistic
  perm.setTests(1);

  vector<double> pr(1,0);
  vector<double> original(1,0);

  perm.setPermClusters(*this);
  perm.originalOrder();

  matrix_t permphe;

  sizeMatrix(permphe, n, par::replicates );
  bool finished = false;
  int j = 0;
  while(!finished)
    {
      
      if (par::perm_genedrop)
	{
	  if (par::perm_genedrop_and_swap)
	    perm.permuteInCluster();
	  perm.geneDrop();
	}
      else
	perm.permuteInCluster();

      for (int i=0; i<n; i++) 
	permphe[i][j] = sample[i]->pperson->phenotype;
      
      finished = perm.update(pr,original);      
      ++j;
    } // next permutation
  
  ofstream PPHE;
  printLOG("Writing permuted phenotype file to [ " + par::output_file_name + ".pphe ]\n");
  PPHE.open( ( par::output_file_name + ".pphe").c_str() , ios::out );
  for (int i = 0 ; i < n ; i++ )
    {
      PPHE << sample[i]->fid << "\t"
	   << sample[i]->iid << "\t";
      for (int j=0; j<par::replicates; j++)
	PPHE << permphe[i][j] << "\t";
      PPHE << "\n";
    }
  PPHE.close();

}

vector<vector<int> > two_locus_table(int l1, int l2)
{
  

  //   0 1 2 tot
  // 0 a b c d 
  // 1 e f g h
  // 2 i g k l
  // M
  // tot

  // i.e. so t[4][4] contains # non-missing individuals
  
  vector< vector<int> > t(5);
  for (int i=0; i<5; i++)
    t[i].resize(5,0);
  
  for (int i=0; i<PP->n; i++)
    {
      Individual * person = PP->sample[i];
      if ( person->missing || ! person->founder ) 
	continue;

      bool a1 = par::SNP_major ? PP->SNP[l1]->one[i] : person->one[l1];
      bool a2 = par::SNP_major ? PP->SNP[l1]->two[i] : person->two[l1];

      bool b1 = par::SNP_major ? PP->SNP[l2]->one[i] : person->one[l2];
      bool b2 = par::SNP_major ? PP->SNP[l2]->two[i] : person->two[l2];
      
      if ( ! a1 ) 
	{
	  if ( ! a2 ) 
	    {
	      if ( ! b1 ) 
		{
		  if ( ! b2 ) 
		    ++t[0][0];
		  else
		    ++t[0][1];
		}
	      else
		{
		  if ( ! b2 ) 
		    ++t[0][3];
		  else
		    ++t[0][2];
		}
	    }
	  else
	    {
	      if ( ! b1 ) 
		{
		  if ( ! b2 ) 
		    ++t[1][0];
		  else
		    ++t[1][1];
		}
	      else
		{
		  if ( ! b2 ) 
		    ++t[1][3];
		  else
		    ++t[1][2];
		}
	    }
	}
      else
	{
	  if ( ! a2 ) 
	    {
	      if ( ! b1 ) 
		{
		  if ( ! b2 ) 
		    ++t[3][0];
		  else
		    ++t[3][1];
		}
	      else
		{
		  if ( ! b2 ) 
		    ++t[3][3];
		  else
		    ++t[3][2];
		}
	    }
	  else
	    {
	      if ( ! b1 ) 
		{
		  if ( ! b2 ) 
		    ++t[2][0];
		  else
		    ++t[2][1];
		}
	      else
		{
		  if ( ! b2 ) 
		    ++t[2][3];
		  else
		    ++t[2][2];
		}
	    }
	}
    }
  
  // Row and col totals
  for (int i = 0; i<4; i++)
    for (int j = 0; j<4; j++)
      {
	t[i][4] += t[i][j];
	t[4][j] += t[i][j];
	if ( i<3 && j<3 ) 
	  t[4][4] += t[i][j];
      }
  return t;
}


map<string, set<Range> > filterRanges(map<string, set<Range> > & ranges, string filename)
{
  set<string> isubset;
  set<string> inotfound;
  
  checkFileExists( filename );
  PP->printLOG("Reading gene subset list from [ " + filename + " ]\n");
  ifstream IN(filename.c_str(), ios::in);
  
  while ( ! IN.eof() )
    {
      string gname;
      IN >> gname;
      if ( gname=="" )
	continue;
      isubset.insert(gname);
    }

  
  // Copy over extracted set to here
  map<string, set<Range> > newRanges;
  set<string>::iterator i = isubset.begin();
  while ( i != isubset.end() )
    {

      map<string, set<Range> >::iterator rf = ranges.find( *i );
      if ( rf == ranges.end() )
	{
	  inotfound.insert( *i );
	  ++i;
	  continue;
	}
      
      newRanges.insert( make_pair( *i , rf->second ) );
      
      ++i;
    }

  PP->printLOG("Extracted " + int2str( newRanges.size() ) + " ranges from this list\n");
  if ( inotfound.size() > 0 )
    {
      PP->printLOG("Was unable to find " + int2str( inotfound.size() ) + " ranges\n");
      PP->printLOG("Writing this list of not-found genes to [ " + par::output_file_name + ".notfound ]\n");
      ofstream O2;
      O2.open( (par::output_file_name+".notfound").c_str() , ios::out);
      set<string>::iterator i1 = inotfound.begin();
      while ( i1 != inotfound.end() )
	{
	  O2 << *i1 << "\n";
	  ++i1;
	}
      O2.close();
    }
  
  return newRanges;
}
