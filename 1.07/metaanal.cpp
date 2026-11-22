

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

#include "stats.h"
#include "options.h"
#include "plink.h"
#include "helper.h"
#include "zed.h"
#include "nlist.h"

using namespace std;

extern Plink * PP;

class SInfo {
public:
  SInfo(double d, double se) : d(d), se(se) { }
  double d;
  double se;
};


class Alleles {
  
public:

  string snp;
  string a1;
  string a2;
  int chr;
  int bp;

  Alleles(string name) : snp(name) 
  {
    a1=a2="";  
    chr = 0;
    bp = 0;
  }

  Alleles(string name, int chr, int bp) : snp(name), chr(chr), bp(bp) 
  {
    a1=a2="";  
  }

  Alleles(string name, int chr, int bp, string a1, string a2) 
    : snp(name),  chr(chr), bp(bp), a1(a1), a2(a2) { }

  bool matches(string b1, string b2) const
  {
    return ( a1 == b1 && a2 == b2 ) || 
      ( a1 == b2 && a2 == b1 ) ;
  }
  
  bool swap(string b1) const
  { return a1 != b1; }
  
  bool operator< (const Alleles & b) const
  {
    // Base first on position, then name
    if ( chr < b.chr ) return true;
    if ( chr > b.chr ) return false;
    if (  bp < b.bp ) return true;
    if (  bp > b.bp ) return false;
    return snp < b.snp;

  }
  
  bool operator== (const Alleles & b) const
  {
    return (chr == b.chr && bp == b.bp && snp == b.snp );
  }
  
};

typedef map<Alleles,map<int,SInfo> > mymap_t;

void Plink::metaAnalysis()
{
  // Read in 2+ files; match SNPs on IDs
  // Find OR,SE and perform meta-analysis
  
  // Information indexed by SNP: each SNP
  // has 1+ files containing that SNP
  // and the needed info
  
  map<Alleles,map<int,SInfo> > store;
    
  printLOG("Performing meta-analysis of " 
	   + int2str( par::meta_files.size() )
	   + " files\n");


  OptionSet * opt = par::opt.getOptions("META");

  bool outputStudyEffects = opt->isSet("study");
  bool usePositions = ! opt->isSet("no-map");
  bool quantTrait = opt->isSet("qt");
  bool allelicInfo = ! opt->isSet("no-allele");
  //bool sensitivity = opt->isSet("drop1");
  bool reportAll = opt->isSet("report-all");
  bool logOR = opt->isSet("logscale");

  // Are sample sizes given?
  bool sampleN = opt->isSet("n");
  vector<double> n1;
  vector<double> n2;
  
  // Use sample N as weights?
  bool nWeights = opt->isSet("n-weight");
  
  bool uWeights = opt->isSet("weight");
  
  if ( nWeights && ! sampleN ) 
    error("Must give sample N's with n-weight option");
  if ( nWeights && uWeights )
    error("Cannot specify both n-weight and weight");

  // By default, look for SE field
  string weightField = "SE";
  if ( uWeights )
    {
      weightField = opt->getValue("weight");
      printLOG("Setting user-defined weight field to " + weightField + "\n");
    }
  
  if ( sampleN )
    {
      // Get sample-N information
      vector<string> ns = opt->getValues("n");
      
      // n=200,300,100,300
      // OR n=200/200,100/150,403,405 etc
      
      if ( ns.size() != par::meta_files.size() )
	error( int2str( ns.size()) 
	       + " sample Ns given, but " 
	       + int2str( par::meta_files.size()) 
	       + " files listed\n");
      
      n1.resize(ns.size());
      n2.resize(ns.size());
      
      for (int f=0; f<ns.size(); f++)
	{
	  if ( quantTrait )
	    {
	      if ( ! from_string<double>( n1[f] , ns[f] , std::dec ) )
		error("Problem with sample N = " + ns[f] );
	    }
	  else
	    {
	      NList l2(0);
	      l2.setRangeChar(" ");
	      l2.setDelimiter("/");
	      vector<string> tok2 = l2.deparseStringList( ns[f] );
	      if ( tok2.size() != 2 )
		error("Problem with sample N = " + display(tok2) + "       Expecting A/U format" );
	      if ( ! from_string<double>( n2[f] , tok2[0] , std::dec ) )
		error("Problem with sample N = " + tok2[0] );
	      if ( ! from_string<double>( n1[f] , tok2[1] , std::dec ) )
		error("Problem with sample N = " + tok2[1] ); 
	    }
	}
      
      double totN1 = 0, totN2 = 0;
      for (int f=0; f<ns.size(); f++)
	{
	  totN1 += n1[f]; 
	  totN2 += n2[f];
	}

      if ( quantTrait ) 
	printLOG("In total, " + dbl2str(totN1) + " individuals\n");
      else
	printLOG("In total, " 
		 + dbl2str(totN2) + " cases and " 
		 + dbl2str(totN1) + " controls\n");
    }
      

  // Are we filtering on a given chromosome or set of SNPs?

  if ( usePositions && par::run_chr > 0 ) 
    printLOG("Processing only chromosome " + int2str(par::run_chr) + "\n");

  if ( par::run_chr > 0 && ! usePositions ) 
    par::run_chr == -1;

  set<string> mset;
  if ( par::extract_set )
    {
      printLOG("Processing only SNPs in [ " + par::extract_file + " ]\n");
      ZInput z2( par::extract_file , compressed(par::extract_file) );
      while( ! z2.endOfFile() )
	{
	  vector<string> tok = z2.tokenizeLine();
	  for (int i=0; i<tok.size(); i++)
	    mset.insert( tok[i] );
	}
      z2.close();
    }
  
  if ( allelicInfo && ! usePositions )
    {
      printLOG("Note: no-map also implies no-allele\n");
      allelicInfo = false;
    }

  // If QT, assume beta and SE from linear regression
  // Otherwise (default), assume OR and SE(log(OR)) as per PLINK output
  //  (and thus take log of OR before combining)

  string dlabel = quantTrait ? "BETA" : "OR";
 
  int twoOrMore = 0;
  int rejected = 0;
  string rdet = "";

  for (int f = 0 ; f < par::meta_files.size() ; f++ )
    {

      checkFileExists( par::meta_files[f] );
      ZInput zin( par::meta_files[f]  , compressed( par::meta_files[f] ) );
      
      printLOG("Reading results from [ " + par::meta_files[f] + " ] ");

      int snps = 0;
      
      // Open a single results file
      // Read first (header) row
      
      string header = zin.readLine();
      vector<string> tokens = tokenizeLine( header );
      
      // Find appropriate columns to filter      
      int snp_column = -1;
      int pval_column = -1;
      int d_column = -1;
      int se_column = -1;
      int chr_column = -1;
      int bp_column = -1;
      int a1_column = -1;
      int a2_column = -1;

      for (int i=0; i<tokens.size(); i++)
	{	  
	  if ( tokens[i] == "SNP" )
	    snp_column = i;	        
	  if ( tokens[i] == "CHR" )
	    chr_column = i;
	  if ( tokens[i] == "BP" )
	    bp_column = i;
	  if ( tokens[i] == "A1" )
	    a1_column = i;
	  if ( tokens[i] == "A2" )
	    a2_column = i;
	  
	  if ( tokens[i] == dlabel )
	    d_column = i;	        
	  if ( tokens[i] == weightField )
	    se_column = i;	        	  
	  if ( tokens[i] == "P" )
	    pval_column = i;	        
	}
      
      if ( snp_column == -1 || 
	   d_column == -1 || 
	   se_column == -1 )
	error("Need fields SNP, OR and SE");
      
      if ( usePositions  && ( chr_column == -1 || bp_column == -1 ) )
	error("No CHR/BP fields in " + par::meta_files[f] );

      if ( allelicInfo && a1_column == -1 ) 
	error("Need A1 (A2) field, or specify no-allele option");
      
      int fsize = tokens.size();
      
      while ( ! zin.endOfFile() )
	{

	  vector<string> tokens = zin.tokenizeLine( );

	  if ( tokens.size() != fsize )
	    continue;

	  string snp = tokens[ snp_column ];

	  // Ignore this SNP?

	  if ( par::extract_set 
	       && mset.find(snp) == mset.end() )
	    continue;

	  double d, se;
	  bool okay = true;
	  
	  // A potential SNP to add

	  Alleles a(snp);

	  if ( usePositions )
	    {
	      a.chr = getChromosomeCode( tokens[ chr_column ] );

	      if ( a.chr == 0 ) 
		{
		  rdet += par::meta_files[f] + "\t" + snp + "\tBAD_CHR\n";
		  okay = false;
		}
	      
	      if ( par::run_chr > 0 && par::run_chr != a.chr )
		continue;
	      
	      if ( ! from_string<int>( a.bp , tokens[ bp_column ] , std::dec ) )
		{
		  rdet += par::meta_files[f] + "\t" + snp + "\tBAD_BP\n";
		  okay = false;
		}

	    }
	  
	  

	  if ( allelicInfo )
	    {
	      a.a1 = tokens[ a1_column ];
	      if ( a2_column != -1 )
		a.a2 = tokens[ a2_column ];

	      // We can only include polymorphic alleles
	      if ( a.a1 == par::missing_genotype )
		{
		  rdet += par::meta_files[f] + "\t" + snp + "\tMISSING_A1\n";
		  okay = false;
		}
	      if ( a2_column != -1 && a.a2 == par::missing_genotype )
		{
		  rdet += par::meta_files[f] + "\t" + snp + "\tMISSING_A2\n";
		  okay = false;
		}
	    }
	  
	  if ( ! from_string<double>( d , tokens[ d_column ] , std::dec ) )
	    {
	      rdet += par::meta_files[f] + "\t" + snp + "\tBAD_ES\n";
	      okay = false;
	    }

	  if ( ! from_string<double>( se , tokens[ se_column ] , std::dec ) )
	    {
	      rdet += par::meta_files[f] + "\t" + snp + "\tBAD_SE\n";
	      okay = false;
	    }

	  // Check alleles?
	  if ( allelicInfo )
	    {

	      mymap_t::iterator k = store.find( a );

	      if ( k != store.end() )
		{
		  if ( k->first.matches(a.a1,a.a2) )
		    {
		      if ( k->first.swap(a.a1) )
			{
			  // Need to swap effect direction			      
			  if ( quantTrait || logOR )
			    d = -d;
			  else
			    d = 1/d;			      
			}
		    }
		      else
			{
			  rdet += par::meta_files[f] + "\t" + snp + "\tALLELE_MISMATCH\n";
			  okay = false;
			}	     
		}
	    }
	

	  
	  if ( ! okay ) 
	    {
	      ++rejected;
	      continue;
	    }
	  
	  // If OR, take log unless told it is already as log
	  if ( ! ( quantTrait || logOR ) )
	    d = log(d);
	  
	  SInfo s(d,se);
	  
	  mymap_t::iterator i = store.find( a );
	  
	  if ( i == store.end() )
	    {
	      map<int,SInfo> t;
	      t.insert(make_pair(f,s));
	      store.insert(make_pair(a,t));
	    }
	  else
	    {
	      if ( i->second.size() == 1 )
		++twoOrMore;
	      i->second.insert(make_pair(f,s));
	    }
	  

	  ++snps;
	  
	}
      zin.close();
      printLOG( " with " + int2str( snps ) + " read\n");
    }
  
  printLOG(int2str(store.size()) + " unique SNPs, " 
	   + int2str(twoOrMore) + " in two or more files\n");
  
  if ( rejected > 0 )
    {
      printLOG("Rejected " + int2str( rejected ) 
	       + " SNPs, writing details to [ " 
	       + par::output_file_name 
	       + ".prob ]\n");
      ZOutput z2( par::output_file_name + ".prob" , false );
      z2.write(rdet);
      z2.close();
    }

  // NOTES
  //  No weights
  //  Assume OR,SE, so LOG taken -- add BETA, etc
  //  Do not check for allelic discrepancy
  
  // Perform meta-analysis

  printLOG("Writing meta-analysis results to [ " 
	   + par::output_file_name 
	   + ".meta ]\n");
  
	     
  ZOutput zout( par::output_file_name + ".meta" , false );

  if ( usePositions )
    zout << sw("CHR",4) << sw("BP",12) ;

  zout << sw( "SNP" , 15 );
  if ( allelicInfo )
    {
      zout << sw( "A1" , 4 )
	   << sw( "A2" , 4 );
    }
  zout << sw( "N" , 4 )
       << sw( "P" , 12 )
       << sw( "P(R)" , 12 )
       << sw( "OR" , 8 )
       << sw( "OR(R)" , 8 )
       << sw( "Q" , 8 )
       << sw( "I" , 8 );

  if( outputStudyEffects )
    for (int f=0; f<par::meta_files.size(); f++)
      zout << sw( "F"+int2str(f) , 8 );

  zout << "\n";
  

  mymap_t::iterator i = store.begin();
  
  
  while ( i != store.end() ) 
    {
      
      int n = i->second.size();
      
      if ( n > 1 ) 
	{

	  // CHR/BP positions?
	  if ( usePositions )
	    zout << sw( i->first.chr, 4) << sw( i->first.bp, 12) ;
	  
	  // SNP name
	  zout << sw( i->first.snp , 15 );
	  
	  // Allele codes?
	  if ( allelicInfo )
	    {
	      zout << sw( i->first.a1 , 4 );
	      if ( i->first.a2 == "" )
		zout << sw( "?" , 4 );
	      else
		zout << sw( i->first.a2 , 4 );
	    }
	  
	  // Number of studies it appears in
	  zout << sw( n , 4 );
	  
	  // vector of OR, SE
	  vector_t d;
	  vector_t se;
	  map<int,SInfo>::iterator j = i->second.begin();
	  while ( j != i->second.end() )
	    {
	      d.push_back( j->second.d );
	      se.push_back( j->second.se );
	      ++j;
	    }
	  
	  // Caclculate weights
	  double denom = 0, denom2 =0, numer = 0;
	  vector_t w(n);
	  vector_t w_random(n);
	  vector_t vars(n);
	  
	  for (int k=0; k<n; k++)
	    {
	      vars[k] = se[k]*se[k];
	      w[k] = 1/vars[k];
	      numer += w[k] * d[k];
	      denom += w[k];
	      denom2 += w[k] * w[k];
	    }

	  double fixedSumm = numer / denom;

	  double Q = 0;
	  for (int k=0; k<n; k++)
	    {
	      double t = d[k]-fixedSumm;
	      t *= t;
	      t /= vars[k];
	      Q += t;
	    }
	  
	  double tau2 = (Q-(n-1)) / ( denom - denom2/denom );
	  if ( tau2 < 0 ) tau2 = 0;      
	  
	  // Random effects model?
	  for (int k=0; k<n; k++)
	    w_random[k] = 1.0/( vars[k] + tau2 );
	  
	  numer = 0;
	  denom = 0;
	  for (int k=0; k<n; k++)
	    {
	      numer += w[k] * d[k];
	      denom += w[k];
	    }
	  double summ = numer / denom;
	  
	  numer = 0;
	  denom = 0;
	  for (int k=0; k<n; k++)
	    {
	      numer += w_random[k] * d[k];
	      denom += w_random[k];
	    }
	  double summ_random = numer / denom;
	  

	  // Variance: fixed effects

	  double varsum = 0;
	  numer = denom = 0;
	  for (int k=0; k<n; k++)
	    {
	      numer += w[k] * w[k] * vars[k];
	      denom += w[k];
	    }
	  varsum = numer / ( denom * denom ); 
	  
	  // Random effects
	  numer = denom = 0;
	  for (int k=0; k<n; k++)
	    {
	      numer += w_random[k] * w_random[k] * ( vars[k] + tau2 );
	      denom += w_random[k];
	    }

	  double varsum_random = numer / ( denom * denom ); 
	  double summtest = summ / sqrt(varsum);
	  double summtest_random = summ_random / sqrt(varsum_random);

	  double p1 = chiprobP( summtest*summtest , 1 );
	  double pR = chiprobP( summtest_random*summtest_random , 1 );
	  double pQ = chiprobP( Q , n-1 );
	  double I = 100 * ( ( Q-(n-1) ) / Q );
	  if ( I < 0 ) I=0;
	  if ( I > 100 ) I= 100;
	  

	  /////////////////////////////////////////
	  // Output
	  
	  if ( ! quantTrait )
	    {
	      summ = exp(summ);
	      summ_random = exp(summ_random);
	    }

	  // Convert -9 to NaN, so sw() handles printing

	  double z = 0;
	  if ( p1 < 0 ) p1 = 1/z;
	  if ( pR < 0 ) pR = 1/z;
	  if ( pQ < 0 ) pQ = 1/z;

	  zout << sw( p1 , -4, 12 )
	       << sw( pR , -4, 12 )
	       << sw( summ , 4, 8 )
	       << sw( summ_random , 4, 8 )
	       << sw( pQ , 4, 8 )
	       << sw( I ,2 , 8 );
	  
	  if( outputStudyEffects )
	    for (int f=0; f<par::meta_files.size(); f++)
	      {
		map<int,SInfo>::iterator k = i->second.find(f);
		if ( k != i->second.end() )
		  zout << sw( exp(k->second.d) , 4, 8 );
		else
		  zout << sw( "NA" , 8 );
	      }
	  
	  zout << "\n";
	}
      else if ( reportAll )	
	{

	  // Display an entry for even single study/zero study 
	  // SNPs
	  
	  if ( usePositions )
	    zout << sw( i->first.chr, 4) << sw( i->first.bp, 12) ;
	  
	  // SNP name
	  zout << sw( i->first.snp , 15 );

	  if ( allelicInfo )
	    {
	      zout << sw( i->first.a1 , 4 );
	      if ( i->first.a2 == "" )
		zout << sw( "?" , 4 );
	      else
		zout << sw( i->first.a2 , 4 );
	    }

	  // Number of studies it appears in
	  zout << sw( n , 4 );
	  
	  zout << sw( "NA" , 12 )
	       << sw( "NA" , 12 )
	       << sw( "NA" , 8 )
	       << sw( "NA" , 8 )
	       << sw( "NA" , 8 )
	       << sw( "NA" , 8 );
	  
	  if( outputStudyEffects )
	    for (int f=0; f<par::meta_files.size(); f++)
	      {
		map<int,SInfo>::iterator k = i->second.find(f);
		if ( k != i->second.end() )
		  zout << sw( exp(k->second.d) , 4, 8 );
		else
		  zout << sw( "NA" , 8 );
	      }
	  
	  zout << "\n";

	}

	  
      // Next SNP

      ++i;
      

    }
  
  zout.close();
  shutdown();

}
  
