

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
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
#include <algorithm>
#include <cmath>

#include "options.h"
#include "plink.h"
#include "helper.h"
#include "sets.h"
#include "stats.h"

extern Plink * PP;

class SortedResult
{
 public:

  double chisq;
  double p;
  int l;

  bool operator< (const SortedResult & s2) const
    {  return (p < s2.p); }  

};


void Plink::setAssocSummary()
{


  // For each set, 

  // 1) get largest unselected test statistic * Ne
  // 2) attenuate all other Ne factor of r^2 to SNP selected in 1
  // 3) repeat to (1) until all SNPs selected 
  // 4) calculate sum test statistic * Ne 



  /////////////////////////////////////
  // Output file

  ofstream SSUM;
  SSUM.open( (par::output_file_name + ".set.summary").c_str() , ios::out);
  SSUM.precision(4);

  

  ////////////////////////////////////
  // Open a single results file
  
  ifstream RESIN;
  checkFileExists( par::set_screen_resultfile );
  RESIN.open( par::set_screen_resultfile.c_str() , ios::in );
  
  // Read first (header) row
  vector<string> tokens = tokenizeLine( RESIN );

  int chisq_column = -1;
  int snp_column = -1;
  int p_column = -1;
  int cols = tokens.size();

  for (int i=0; i<tokens.size(); i++)
    {
      if ( tokens[i] == "CHISQ" )
	chisq_column = i;
      if ( tokens[i] == "SNP" )
	snp_column = i;	        
      if ( tokens[i] == "P" )
	p_column = i;
    }
  
  if ( ( chisq_column == -1 && p_column == -1 ) || snp_column == -1 ) 
    error("Need a file with SNP and CHISQ or P fields\n");
  
  
  
  // -1 code means not observed in results file
  vector_t stat(nl_all,-1);
  vector_t p(nl_all,-1);

  
  map<string,int> mlocus;
  makeLocusMap(*this,mlocus);
  set<int> slist;

  bool convert = chisq_column == -1 || p_column == -1 ;
  bool from_p = chisq_column == -1;

  while ( !RESIN.eof() )
    {
      vector<string> tokens = tokenizeLine( RESIN );
      if ( tokens.size() != cols )
	continue;

      string snp = tokens[snp_column];
      
      map<string,int>::iterator i1 = mlocus.find( snp );
      if ( i1 == mlocus.end() )
	continue;

      double x2, pv;

      if ( ! convert )
	{
	  if ( ! from_string<double>( x2, tokens[chisq_column] , std::dec))
	    continue;
	  
	  if ( ! from_string<double>( pv, tokens[p_column] , std::dec))
	    continue;
	}
      else
	{
	  if ( from_p )
	    {
	      if ( ! from_string<double>( pv, tokens[p_column] , std::dec))
		continue;
	      x2 = inverse_chiprob( pv, 1 );
	    }
	  else
	    {
	      if ( ! from_string<double>( x2, tokens[chisq_column] , std::dec))
		continue;
	      pv = chiprobP(x2,1);
	    }
	}

      
      // Bounds on possible p-values

      if ( pv == 1 ) 
	pv = 1-1e-10;
      
      if ( pv == 0 ) 
	pv = 1e-10;
      
      stat[ i1->second ] = x2;
      p[ i1->second ] = pv;
      slist.insert( i1->second );
    }

  RESIN.close();

  printLOG("Read results for " + int2str( slist.size() ) + " SNPs\n");
  printLOG("Writing set summary statistics to [ " + par::output_file_name + ".set.summary ]\n");


  //////////////////////////////////
  // Header row

  SSUM << setw(22) << "SET" << " " 
       << setw(6) << "NSNP" << " " 
       << setw(32) << "POS" << " "
       << setw(8) << "KB" << " "
       << setw(12) << "P1" << " "
       << setw(12) << "P2" << "\n";

      

  ofstream SVERB;
  if ( par::verbose )
    SVERB.open( (par::output_file_name+".set.summary.verbose").c_str(), ios::out );


  // Substract from all other SNPs this value times the r^2 with the other SNP


  /////////////////////////////////////
  // Score each Set
 
  for ( int j = 0; j < pS->snpset.size(); j++ )
    {


      // Gather and sort test statistics (chi-sq)

      vector<SortedResult> t;

      for (int i=0; i < snpset[j].size(); i++)
 	{
	  
	  if ( p[snpset[j][i]] > 0 && p[snpset[j][i]] <= 1 )
	    {
	      SortedResult s;
	      s.chisq = stat[snpset[j][i]];
	      s.l = snpset[j][i];
	      s.p = p[snpset[j][i]];
	      t.push_back(s);
	    }
 	}	  


      int ns = t.size();

      if ( ns == 0 ) 
	{
	  SSUM << setw(22) << setname[j] << " "
	       << setw(6) << 0 << " " 	       
	       << setw(32) << "NA" << " " 
	       << setw(8) << "NA" << " " 
	       << setw(12) << "NA" << " "
	       << setw(12) << "NA" << "\n";

	  continue;
	}


      // Sort statistics (in descrending order)
      
      sort(t.begin(),t.end());
      


      //////////////////////////////////////////////
      //                                          //
      // Makambi (2003), Delongchamp et al (2006) //
      //                                          //
      //////////////////////////////////////////////


      // Makambi; 
      // assume uniform weighting
      
      vector_t w(ns, 1.0/(double)ns);
      double var = 0;
      for (int l=0; l<ns; l++)
	var += w[l] * w[l];
      var *= 4;
      
      // Delongchamp

      double denom = ns;


      //////////////////////////////////////////////
      //                                          //
      //  Loop over pairs of SNPs                 // 
      //                                          //
      //////////////////////////////////////////////
      
      
      for (int l1=0; l1 < ns-1; l1++)
	{
	  int chr = locus[l1]->chr;
	  int bp = locus[l1]->bp;

	  if ( ! par::silent )
	    cout << "Set " << j << ", SNP " << l1 << " of " << ns << "             \r";
	  
	  for (int l2=l1+1; l2<ns; l2++)
	    {
	      
	      // Assume r^2=0 for effectively unlinked SNPs
	      if ( locus[l2]->chr != chr )
		continue;
	      if ( abs( bp - locus[l2]->bp ) > 1000000 ) 
		continue;
	      
	      double r = correlation2SNP( t[l1].l , t[l2].l, false, false);
	      
	      if ( ! realnum(r) ) 
		r = 0;
	      r = abs(r);
	      
	      // Makambi
	      var += 2 * w[l1] * w[l2] * ( 3.25 * r + 0.75 * r*r ) ;	    
	      // Delongchamp
	      denom += 2 * r;
	      	      
	    }
	}
      
	  
      // Makambi
      double df = 8.0 / var;
      if ( df < 2 ) df = 2;
      double score = 0;
      for (int l=0; l<ns; l++)
	score += w[l] *  -2 * log( chiprobP( t[l].chisq , 1 ) ); 
      score = df * ( score / 2.0 );
      double pv = chiprobP( score , df );
      
      
      // Delongchamp
      double numerator = 0;
      for (int l=0; l<ns; l++)
	numerator += ltqnorm( 1 - t[l].p );
      denom = sqrt(denom);
      double stat2 = numerator / denom;
      double pv2 = 1-normdist(stat2);
      

      
      //////////////////////////////////////////////
      //                                          //
      //    Output                                //
      //                                          //
      //////////////////////////////////////////////

      
      SSUM << setw(22) << setname[j] << " "
	   << setw(6) << ns << " ";

      int chr = locus[ t[0].l ]->chr;
      int minBP = locus[ t[0].l ]->bp;
      int maxBP = minBP;
      bool diffchr = false;

      for (int l=1;l<ns;l++)
	{
	  if ( locus[ t[l].l ]->chr != chr )
	    {
	      diffchr = true;
	      break;
	    }

	  if ( locus[ t[l].l ]->bp > maxBP ) 
	    maxBP = locus[ t[l].l ]->bp;
	  else if ( locus[ t[l].l ]->bp < minBP )
	    minBP = locus[ t[l].l ]->bp;
	  
	}
      
      if ( diffchr )
	{
	  SSUM << setw(32) << "NA" << " " 
	       << setw(8) << "NA" << " ";
	}
      else
	{
	  string pstr = "chr" 
	    + int2str(chr) + ":" 
	    + int2str( minBP) + ".." 
	    + int2str( maxBP );
	  SSUM << setw(32) << pstr << " " 
	       << setw(8) << (maxBP-minBP)/1000.0 << " ";
	}
	
      SSUM << setw(12) << pv << " "
	   << setw(12) << pv2 << "\n";
      

    }
  
  if ( ! par::silent )
    cout << "                                                           \n";


  SSUM.close();

  if ( par::verbose )
    SVERB.close();

}

