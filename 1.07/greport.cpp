

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

#include "options.h"
#include "plink.h"
#include "helper.h"

extern Plink * PP;


void Plink::displayGeneReport()
{

  // Simply read in any generic results file and list of SNPs by
  // ranges (which may be subsetted).
  
  //   if ( false )
  //     readMapFile(par::mapfile,include,include_pos,nl_actual);

  ofstream GREP;

  GREP.open( (par::output_file_name + ".range.report").c_str() , ios::out);
  
  map<string, set<Range> > ranges;
  

  // Read list of ranges
  ranges = readRange( par::greport_gene_list );

  // Filter ranges 
  
  if ( par::greport_subset ) 
    ranges = filterRanges( ranges, par::greport_subset_file );



  // Open a single results file
  
  ifstream RESIN;
  RESIN.open( par::greport_results.c_str() , ios::in );
  
  // Read first (header) row

  char cline[par::MAX_LINE_LENGTH];
  RESIN.getline(cline,par::MAX_LINE_LENGTH,'\n');
  
  string sline = cline;
  if (sline=="") 
    error("Problem reading [ " + par::greport_results + " ]\n");

  string buf; 
  stringstream ss(sline); 
  vector<string> tokens; 
  while (ss >> buf)
    tokens.push_back(buf);

  int chr_column = -1;
  int bp_column = -1;
  int pval_column = -1;
  int snp_column = -1;

  for (int i=0; i<tokens.size(); i++)
    {
      if ( tokens[i] == "CHR" )
	chr_column = i;
      
      if ( tokens[i] == "BP" )
	bp_column = i;	  
      
      if ( tokens[i] == "SNP" )
	snp_column = i;	        

      if ( tokens[i] == "P" )
	pval_column = i;	        
    }


  // Do we have a list of SNPs to specifically extract?

  set<string> extractSNP;
  if ( par::extract_set )
    {
      if ( snp_column == -1 ) 
	error("Did not find a SNP field, so cannot use --extract");

      checkFileExists( par::extract_file );
      PP->printLOG("Only extracting SNPs listed in [ " + par::extract_file + " ]\n");
      ifstream IN(par::extract_file.c_str(), ios::in);
      
      while ( ! IN.eof() )
	{
	  string snpname;
	  IN >> snpname;
	  if ( snpname=="" )
	    continue;
	  extractSNP.insert(snpname);
	}
      IN.close();      
      PP->printLOG("Read " + int2str( extractSNP.size() ) + " SNPs to extract\n");
    }


  if ( chr_column < 0 || bp_column < 0 )
    error("Could not find CHR and BP fields in results file");
  
  map<Range*,vector<string> > annotatedResults;

  string headerline = sline;
  int cnt = 0;
  while ( ! RESIN.eof() )
    {

//       if ( ! par::silent ) 
// 	cout << "Processing results line " << ++cnt << "        \r";

      //      vector<string> tokens = tokenizeLine( RESIN ); 

      char cline[par::MAX_LINE_LENGTH];
      RESIN.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      string sline = cline;
      if (sline=="") 
	continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
    
      if ( tokens.size() <= chr_column ||
	   tokens.size() <= bp_column )
	continue;
      
      // Using a p-value-filtering field? 

      double pvalue = 0;
      if ( pval_column != -1 )
	{
	  if ( tokens.size() <= pval_column )
	    continue;
	  
	  if ( ! from_string<double>( pvalue, tokens[pval_column] , std::dec))
	    continue;
	  
	  if ( par::pfilter && pvalue > par::pfvalue ) 
	    continue;
	  
	}

      if ( par::extract_set ) 
	{
	  if ( tokens.size() <= snp_column )
	    continue;
	  
	  if ( extractSNP.find( tokens[snp_column] ) == extractSNP.end() )
	    continue;
	}

      int thisChr = -1;
      int thisBP = -1;

      if ( ! from_string<int>( thisChr, tokens[chr_column] , std::dec))
	continue;

      if ( ! from_string<int>( thisBP, tokens[bp_column] , std::dec))
	continue;
      
      // Do we need to store this? i.e. what ranges is it actually in?
      // This information is in snp2range

      Range r1(thisChr,thisBP,thisBP,"dummy"); 
      
      set<Range*> implicated = rangeIntersect(r1,ranges);
      set<Range*>::iterator ri = implicated.begin();
      while ( ri != implicated.end() )
	{

	  string distance = dbl2str(( thisBP - ((*ri)->start + par::make_set_border)) /1000.00 , 4 ) + "kb" ;
	  
	  if ( annotatedResults.find( *ri ) == annotatedResults.end() )
	    {
	      vector<string> t(2);
	      t[0] = distance;
	      t[1] = sline;
	      annotatedResults.insert(make_pair( (Range *)(*ri) , t ) );
	    }
	  else
	    {
	      vector<string> & v = annotatedResults.find( *ri )->second;
	      v.push_back(distance);
	      v.push_back(sline);
	    }

	  ++ri;
	}

      // Read next line of results

    }


      
  // Iterate through these -- they will be in genomic order, hopefully
  
  map<string, set<Range> >::iterator ri = ranges.begin();

  while ( ri != ranges.end() )
    {
      set<Range>::iterator si = ri->second.begin();
      

      while ( si != ri->second.end() )
	{
	  
	  bool displayed = false;
	  
	  	  
	  map<Range*,vector<string> >::iterator ari;
	  ari = annotatedResults.find( (Range *)&(*si) );
	  	 	  
	  if ( ari != annotatedResults.end() )
	    {
	      for (int l=0; l< ari->second.size(); l+=2)
		{
		  if ( ! displayed ) 
		    {
		      GREP << ri->first << " -- chr" 
			   << chromosomeName( si->chr ) << ":" 
			   << si->start << ".."
			   << si->stop << " ( "
			   << (si->stop - si->start ) / 1000.00 << "kb ) ";
		      if ( par::make_set_border > 0 ) 
			GREP << " including " << par::make_set_border/1000.00 << "kb border ";		      
		      GREP << "\n\n" 
			   << setw(12) << "DIST" << " "
			   << headerline << "\n";
		      displayed = true;
		    }
		  
		  GREP << setw(12) << ari->second[l] << " "
		       << ari->second[l+1] << "\n";
		}
	    }
	  
	  if ( ! displayed ) 
	    {
	      if ( par::greport_display_empty ) 
		{
		  GREP << ri->first << " -- chr" 
		       << chromosomeName( si->chr ) << ":" 
		       << si->start << ".."
		       << si->stop << " ( "
		       << (si->stop - si->start ) / 1000.00 << "kb ) ";
		  if ( par::make_set_border > 0 ) 
		    GREP << " including " << par::make_set_border/1000.00 << "kb border ";		      
		  GREP << "   { nothing to report }\n\n";
		}
	    }
	  else
	    GREP << "\n\n";
	  
	  ++si;
	}
            
      ++ri;
    }
	 
  RESIN.close(); 
  GREP.close();

  if ( ! par::silent ) 
    cout << "\n";

  printLOG("Writing per-range report to [ " 
	   + par::output_file_name 
	   + ".range.report ]\n");

  shutdown();

}



