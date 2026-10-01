

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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <map>

#include "plink.h"
#include "stats.h"
#include "helper.h"
#include "options.h"

using namespace std;

void Plink::filterQualSNPs()
{
  
  // Remove entire SNPs if they do not pass the qual score
  // threshold

  vector<bool> del( nl_all, false );

  /////////////////////////////
  // Look-up table by SNP name

  map<string,int> mlocus;
  map<string,int>::iterator ilocus;
  
  for ( int l = 0 ; l < nl_all ; l++ )
    mlocus.insert(make_pair( locus[l]->name, l ) );
  
  checkFileExists( par::snp_qual_file );

  printLOG("Reading SNP quality scores from [ " + par::snp_qual_file + " ]\n");
  ifstream I( par::snp_qual_file.c_str() , ios::in );

  int ndel = 0;
  int nfound = 0;

  while ( ! I.eof() )
    {
      string snp;
      double qual;

      I >> snp >> qual;
      
      if ( snp == "" )
	continue;

      ilocus = mlocus.find( snp );

      if ( ilocus != mlocus.end() )
	{
	  ++nfound;
	  
	  if ( qual < par::snp_qual_min || 
	       qual > par::snp_qual_max ) 
	    {
	      ++ndel;
	      del[ ilocus->second ] = true;	  
	    }
	}      
    }
  
  I.close();

  printLOG("Read quality scores for " + int2str(nfound) + " of " + int2str(nl_all) + " SNPs\n");
  printLOG("Removing " + int2str(ndel) + " SNPs based on quality scores\n");
  

  ////////////////////////////////////////
  // Remove selected loci from locus list, 
  
  deleteSNPs(del);
  
}


void Plink::filterQualGenotypes()
{
  
  
  // Only blank out genotypes if they do not meet the quality 
  // score threshold

  
  // ?To add: an automatic option to treat these as obligatory missing?
  
  
  ////////////////////////////
  // Look-up table by SNP name

  map<string,int> mpeople;
  map<string,int> mlocus;

  for ( int l = 0 ; l < nl_all ; l++ )
    mlocus.insert(make_pair( locus[l]->name, l ) );
  
  for ( int i = 0 ; i < n ; i++ )
    mpeople.insert(make_pair( sample[i]->fid+"_"+sample[i]->iid, i ) );
  
  checkFileExists( par::geno_qual_file );

  printLOG("Reading genotype quality scores from [ " + par::geno_qual_file + " ]\n");
  ifstream I( par::geno_qual_file.c_str() , ios::in );

  long int ndel = 0;
  long int nfound = 0;
  
  
  // Format 0) Q
  //        1) FID/IID  
  //        2) SNP
  //        3) QUAL 
  
  // But can wild card, either person or SNP
  
  // Q * rs12345 {list all qual scores for rs12345 for all people (order as file)
  // Q P1 I1 *   {list all qual scores for person P1 I1 (order as file) } 
  // Q * * { list all qual scores for each person, for each SNP } 


  printLOG("Acceptable genotype quality score range is " + dbl2str( par::geno_qual_min ) + 
	   " to " + dbl2str( par::geno_qual_max ) + "\n");

  string p, m;
   
  while (!I.eof())
    {
      
      // Expecting a "Q" 
      
      string q;
      I >> q;

      if ( q == "" ) 
	continue;

      if ( q != "Q" && q != "q" )
	error("Problem with file format: leading 'Q' not found\n");
      
      
      // Read person

      int pcode = -1;
      bool person_wildcard = false;
      string fid , iid;

      I >> fid;

      if ( fid == "*" ) 
	person_wildcard = true;
      else
	{
	  I >> iid;
	  string pstring = fid + "_" + iid;            

	  map<string,int>::iterator i = mpeople.find( pstring );

	  if ( i != mpeople.end() ) 
	    pcode = i->second;
	  
	}
      

      // Read SNP 
      
      int scode = -1;
      string sstring;
      bool snp_wildcard = false;
      I >> sstring;
      if ( sstring == "*" ) 
	snp_wildcard = true;
      else
	{
	  map<string,int>::iterator i = mlocus.find( sstring );
	  if ( i != mlocus.end() ) 
	    scode = i->second;
	}


      // Now read quality score

      int pstart = pcode;
      int pstop = pcode;
      if ( person_wildcard ) 
	{
	  pstart = 0;
	  pstop = n - 1;
	}

      int sstart = scode;
      int sstop = scode;      
      if ( snp_wildcard ) 
	{
	  sstart = 0;
	  sstop = nl_all - 1;
	}
      
      for ( int p = pstart ; p <= pstop ; p++ )
	for ( int s = sstart ; s <= sstop ; s++ )
	  {
	    // Read qual score
	    
	    double q;
	    I >> q;

	    if ( p >= 0 && s >= 0 )
	      { 
		if ( q < par::geno_qual_min || 
		     q > par::geno_qual_max )
		  {
		    ++ndel;
		    
		    // Assume SNP major
		    
		    // Set to missing
		    SNP[s]->one[p] = true;
		    SNP[s]->two[p] = false;		
		  }
	      }
	  }
    }
  
  

  I.close();
  
  printLOG(int2str(ndel) + " genotypes did not meet quality score, set to missing\n");
  
}



