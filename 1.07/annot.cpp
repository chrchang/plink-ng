

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
#include <algorithm>

#include "options.h"
#include "plink.h"
#include "helper.h"
#include "zed.h"

extern Plink * PP;

map<string, set<Range> > filterRanges(map<string, set<Range> > & ranges, string filename);

void Plink::annotateFile()
{
  
  // Simply read in any generic results file and list of SNPs by
  // ranges (which may be subsetted). Input could be compressed
  
  
  checkFileExists( par::annot_filename );
   ZInput zin( par::annot_filename , compressed( par::annot_filename ) );

   // If input is compressed, also make output compressed

   string f = par::output_file_name + ".annot";
   if ( compressed( par::annot_filename ) ) f += ".gz";
   ZOutput zout( f , compressed( par::annot_filename ) );

   printLOG("Reading input from [ " + par::annot_filename + " ]\n");
   printLOG("Writing annotated file to [ " + f + " ]\n");


   // Range information
   map<string, set<Range> > ranges;

   // SNP attribute information
   map<string,set<string> > attrib;

   // Only output rows that correspond to these filters
   set<string> snp_filter;
   map<string, set<Range> > range_filter;



   // Read list of ranges
   OptionSet * annot_opt = par::opt.getOptions("ANNOT");
   if ( annot_opt->isSet("ranges") )
     {
       ranges = readRange( annot_opt->getValue("ranges") );      
       // Filter to a subset of attribs?
       if ( annot_opt->isSet("subset") )
	 ranges = filterRanges( ranges, annot_opt->getValue("subset") );
     }


  // Read list of SNP annotations, which can be compressed
  // Just simple format: rs-number , 1 + space-delimited annotations
  
  if ( annot_opt->isSet("attrib") )
    {
      string fname = annot_opt->getValue("attrib");
      checkFileExists( fname );      
      ZInput ZIN1( fname , compressed(fname) );
      while ( ! ZIN1.endOfFile() )
	{
	  vector<string> tok = ZIN1.tokenizeLine();
	  for (int j=1; j<tok.size(); j++)
	    {	      
	      map<string,set<string> >::iterator i = attrib.find( tok[0] );
	      if ( i == attrib.end() )
		{
		  set<string> t;
		  t.insert( tok[j]) ;
		  attrib.insert(make_pair( tok[0] , t ) );
		}
	      else
		i->second.insert( tok[j] );
	    }
	}
      ZIN1.close();
      printLOG("Read attributes for " + int2str(attrib.size()) + " SNPs\n");
    }
  
  
  // Filters?
    
  if ( annot_opt->isSet("filter") )
    range_filter = readRange( annot_opt->getValue("filter") );

  if ( annot_opt->isSet("snps") )
    {
      checkFileExists( annot_opt->getValue("snps") );
      ifstream IN1(annot_opt->getValue("snps").c_str() , ios::in );
      while ( ! IN1.eof() )
	{
	  string s;
	  IN1 >> s;
	  if ( s == "" ) continue;
	  snp_filter.insert(s);
	}
      IN1.close();
    }

  bool hasRanges = ranges.size() > 0;
  bool hasSNPs = attrib.size() > 0;
  bool filterRanges = range_filter.size() > 0;
  bool filterSNPs = snp_filter.size() > 0;
  bool needPosition = hasRanges || filterRanges;

  // Default is to output all rows; however, if the 
  // 'annot-only' option is set, then only output a 
  // row that has at least some annotation

  bool output_all = annot_opt->isSet("prune") ? false : true ;
  
  if ( ! ( hasRanges || hasSNPs || filterRanges || filterSNPs ) )
    error("Nothing to do -- stopping");
  

  // Open a single results file
  // Read first (header) row

  string header = zin.readLine();
  vector<string> tokens = tokenizeLine( header );
  
  // Find appropriate columns to filter

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


  if ( needPosition
       && ( chr_column < 0 || bp_column < 0 ) ) 
    error("Could not find CHR and BP fields in results file");

  // Always print distance field 
  bool track_distance = hasRanges && annot_opt->isSet("distance");
  
  // Minimal range output (i.e. not (distkb)
  bool minimal = annot_opt->isSet("minimal");

  // Use "NA" or "." for missing fields
  string missingValue = annot_opt->isSet("NA") ? "NA" : ".";
  
  // block0/1 reporting
  // find all possible annotations, then write fields on 0/1s for each variant
  bool block01 = annot_opt->isSet("block");
  
  

  // Determine all unique annotations
  set<string> uniqFields;

  if ( block01 ) 
    {

      map<string,set<string> >::iterator j = attrib.begin();
      while ( j != attrib.end() )
	{
	  set<string>::iterator k = j->second.begin();
	  while ( k != j->second.end() )
	    {
	      if ( uniqFields.find( *k ) == uniqFields.end() )
		uniqFields.insert( *k );
	      ++k;
	    }
	  ++j;
	}

      map<string, set<Range> >::iterator i = ranges.begin();
      while ( i != ranges.end() )
	{
	  if ( uniqFields.find( i->first ) == uniqFields.end() )
	    uniqFields.insert( i->first );
	  ++i;
	}
      printLOG("Found " + int2str( uniqFields.size() ) 
	       + " unique annotations\n");
    }
  
  // Write header back out, with additional field
  
  if ( block01 ) 
    {
      zout << header;
      set<string>::iterator l = uniqFields.begin();
      while ( l != uniqFields.end() )
        {
          zout << " " << *l;
          ++l;
        }
      zout << "\n";
    }
  else
    {
      if ( track_distance )
	zout << header << sw("DIST",12) << sw("SGN",12) << " ANNOT\n";
      else
	zout << header << " ANNOT\n";
    }

  int cnt = 0, cnt2 = 0;

  while ( ! zin.endOfFile() )
    {

      // Get line of output
      
      string input = zin.readLine();
      vector<string> tokens = tokenizeLine(input);

      if ( tokens.size() == 0 ) continue;

      if ( needPosition )
	{ 
	  if ( tokens.size() <= chr_column ||
	       tokens.size() <= bp_column )
	    continue;
	}

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

      // Filtering on pre-specified SNP names?

      if ( filterSNPs ) 
	{
	  if ( tokens.size() <= snp_column )
	    continue;
	  
	  if ( snp_filter.find( tokens[snp_column] ) == snp_filter.end() )
	    continue;
	}


      int thisChr = -1;
      int thisBP = -1;
      
      if ( needPosition )
	{
	  if ( ! from_string<int>( thisChr, tokens[chr_column] , std::dec))
	    continue;
	  if ( ! from_string<int>( thisBP, tokens[bp_column] , std::dec))
	    continue;
	}

      Range r1(thisChr,thisBP,thisBP,"dummy"); 
      
      // Filtering on a set of ranges?
      
      if ( filterRanges )
	{
	  bool include = false;	  
	  set<Range*> implicated = rangeIntersect(r1,range_filter);
	  if ( implicated.size() == 0 )
	    continue;
	}
      
      
      // Annotation to build up, if any
      string annotation = "";
      
      // If we need to track what we see (for block01 output)
      set<string> x;

      // 1) Ranges
      
      int min_distance = 999999999;
      int sign = 0;

      // Do we need to store this? i.e. what ranges is it actually in?
      // This information is in snp2range
    
      
      // Does this point overlap with any ranges of interest?
      
      if ( hasRanges )
	{
	  set<Range*> implicated = rangeIntersect(r1,ranges);      
	  set<Range*>::iterator ri = implicated.begin();
	  while ( ri != implicated.end() )
	    {
	      
	      string distance = "0";
	      if( thisBP < (*ri)->start + par::make_set_border ) 
		{
		  distance = "-" 
		    + dbl2str(( ( (*ri)->start + par::make_set_border ) - thisBP ) / 1000.00 , 4 ) + "kb" ;
		  
		  if ( track_distance )
		    if ( ( (*ri)->start + par::make_set_border ) - thisBP < min_distance )
		      {
			min_distance = ( (*ri)->start + par::make_set_border ) - thisBP;
			sign = -1;
		      }
		}
	      else if ( thisBP > (*ri)->stop - par::make_set_border )
		{
		  distance = "+" 
		    + dbl2str( ( thisBP - ( (*ri)->stop - par::make_set_border ) ) / 1000.00 , 4 ) + "kb" ;
		  
		  if ( track_distance )
		    if ( thisBP - ( (*ri)->stop - par::make_set_border ) < min_distance ) 
		      {
			min_distance =  thisBP - ( (*ri)->stop - par::make_set_border );		  
			sign = 1;
		      }
		}
	      else
		{
		  min_distance = 0;
		  sign = 0;
		}
	      
	      if ( annotation == "" )
		annotation += (*ri)->name;
	      else
		annotation += "|" + (*ri)->name;
	      
	      if ( ! minimal ) 
		annotation += "(" + distance + ")";
	      
	      // Do we need to track this?
	      if ( block01 )
		x.insert( (*ri)->name );
	      
	      ++ri;
	    }
	}

      
      // 2) Attributes
      
      if ( hasSNPs ) 
	{

	  map<string,set<string> >::iterator i = attrib.find(tokens[snp_column]);
	  if ( i != attrib.end() )
	    {
	      set<string>::iterator j = i->second.begin();
	      while ( j != i->second.end() )
		{
		  if (annotation=="" )
		    annotation += *j;
		  else
		    annotation += "|" + *j;
		  
		  // Do we need to track this?
		  if ( block01 )
		    x.insert( *j );		  
		  
		  ++j;
		}

	    }
	  
	}
      
      // Output this row (or possibly not)
      
      if ( block01 )
	{
	  zout << input << " ";
	  set<string>::iterator l = uniqFields.begin();
	  while ( l != uniqFields.end() )
	    {
	      if ( x.find( *l ) != x.end() ) 
		zout << " 1";
	      else 
		zout << " 0";
	      ++l;
	    }
	  zout << "\n";	  
	}
      else if ( annotation != "" )
	{
	  ++cnt2;

	  if ( track_distance ) 
	    {
	      zout << input << sw( min_distance / 1000.0 , 12 );
	      if ( sign == -1 ) zout << sw("-",4);
	      else if ( sign == 1 ) zout << sw("+",4);
	      else if ( sign == 0 ) zout << sw(missingValue,4);
	      zout << " " << annotation << "\n";
	    }
	  else
	    zout << input << " " << annotation << "\n";
	  
	}
      else if ( output_all )
	{
	  if ( track_distance ) 
	    zout << input << sw("NA",12) << sw("NA",4) << " " << missingValue << "\n";
	  else 
	    zout << input << " " << missingValue << "\n";
	}

      ++cnt;

      // Read next line of results

    }



  printLOG("Processed " + int2str(cnt) + " rows");

  if ( !block01 )
    printLOG(", " + int2str(cnt2) + " of which were annotated");

  printLOG("\n");

  zin.close();
  zout.close();
  shutdown();

}



