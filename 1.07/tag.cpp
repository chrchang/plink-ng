

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

void Plink::tagMode()
{
  

  ////////////////////////////
  // Look-up table by SNP name
  
  map<string,int> mlocus;
  for (int l = 0 ; l < nl_all ; l++ )
    mlocus.insert(make_pair( locus[l]->name, l ));
 
  set<int> toTag;
  bool testAll = par::gettag_file == "all";
  
  if ( !testAll ) 
    {
      
      checkFileExists(par::gettag_file);
      printLOG("Reading SNPs to tag from [ " 
	       + par::gettag_file + " ] \n");
      
      ifstream IN2(par::gettag_file.c_str(),ios::in);
      
 
      int cnt = 0;
      while ( ! IN2.eof() )
	{
	  
	  string snp;
	  string code;
	  if ( par::gettag_mode1 )
	    IN2 >> snp;
	  else if ( par::gettag_mode2 )
	    {
	      vector<string> tokens = tokenizeLine(IN2);
	      if ( tokens.size() == 0 )
		continue;
	      if ( tokens.size() != 2 )
		error("Expected two columns per line\n");
	      snp = tokens[0];
	      code = tokens[1];
	      if ( code == "0" )
		continue;
	    }
      
	  if ( snp == "" ) 
	    continue;
	  
	  ++cnt;
	  
	  // Can we find this SNP?
	  
	  map<string,int>::iterator i = mlocus.find( snp );
	  if ( i == mlocus.end() )
	    continue;
      
	  // Add to list 
	  
	  toTag.insert(i->second);
	  
	}
      
      printLOG("Read " + int2str(cnt) 
	       + " SNPs to tag, of which " 
	       + int2str(toTag.size()) 
	       + " are unique and present\n");
      
      IN2.close();
    }
  else
    {
      printLOG("Setting to tag all " + int2str(nl_all) + " SNPs in dataset\n");
      for (int l=0;l<nl_all;l++)
	toTag.insert(l);
    }
    
    
  
  ///////////////////////
  // Find tags and output

  
  // input
  //  rs0001
  //  rs0004
  
  // output
  //  rs0001
  //  rs0002
  //  rs0003
  //  rs0004

  ofstream O2;
  if ( par::gettag_listall ) 
  {
    printLOG("Writing list of specific tags per SNP to [ " + par::output_file_name + ".tags.list ]\n");
    O2.open( ( par::output_file_name + ".tags.list").c_str() , ios::out );
    O2 << setw( par::pp_maxsnp ) << "SNP" << " " 
       << setw(4) << "CHR" << " " 
       << setw(10) << "BP" << " " 
       << setw(4) << "NTAG" << " "
       << setw(10) << "LEFT" << " " 
       << setw(10) << "RIGHT" << " "
       << setw(8) << "KBSPAN" << " "
       << "TAGS\n";
    
  }

  set<int>::iterator i = toTag.begin();
  set<int> tagged;

  while ( i != toTag.end() )
    {
      // SNP to tag
      int l = *i;
      tagged.insert(l);
      
      set<int> thisTagged;
      int dist_left = locus[l]->bp;
      int dist_right = locus[l]->bp;

      // Move forwards and backwards, within range, and 
      // add any with r^2 above threshold to tagged set
      
      int j = l - 1;
      int chr = locus[l]->chr;
      int pos = locus[l]->bp;

      while (1)
	{
	  
	  if ( j < 0 ) 
	    break;
	  if ( locus[j]->chr != chr )
	    break;
	  if ( pos - locus[j]->bp > par::gettag_kb )
	    break;

	  if ( ! par::gettag_listall )
	    {
	      if ( tagged.find(j) != tagged.end() )
		{
		  --j;
		  continue;
		}
	    }
	  
	  double rsq = correlation2SNP( l,j,true,false);
	  if ( realnum(rsq) && rsq >= par::gettag_r2 )
	    {
	      tagged.insert(j);
	      if ( par::gettag_listall ) 
		{
		  thisTagged.insert(j);
		  if ( locus[j]->bp < dist_left )
		    dist_left = locus[j]->bp;
		}
	    }
	  --j;
	}
      

      // Now move right

      j = l+1;
      while (1)
	{
	  
	  if ( j == nl_all ) 
	    break;
	  if ( locus[j]->chr != chr )
	    break;
	  if ( locus[j]->bp - pos > par::gettag_kb )
	    break;
	  if ( ! par::gettag_listall )
	    {
	      if ( tagged.find(j) != tagged.end() )
		{
		  ++j;
		  continue;
		}
	    }

	  double rsq = correlation2SNP( l,j,true,false);
	  if ( realnum(rsq) && rsq >= par::gettag_r2 )
	    {
	      tagged.insert(j);
	      if ( par::gettag_listall )
		{
		  thisTagged.insert(j);
		  if ( locus[j]->bp > dist_right )
		    dist_right = locus[j]->bp;
		}
	    }
	  ++j;
	}
      
      if ( par::gettag_listall )
	{

	  O2 << setw( par::pp_maxsnp ) << locus[l]->name << " " 
	     << setw(4) << locus[l]->chr << " " 
	     << setw(10) << locus[l]->bp << " " 
	     << setw(4) << thisTagged.size() << " "
	     << setw(10) << dist_left << " " 
	     << setw(10) << dist_right << " "
	     << setw(8) << (dist_right - dist_left ) / 1000.0 << " ";

	  set<int>::iterator s = thisTagged.begin();
	  bool first = true;
	  while ( s != thisTagged.end() )
	    {
	      if ( first ) 
		{
		  first = false;		  
		}
	      else
		O2 << "|";
	      O2 << locus[*s]->name;
	      ++s;
	    }
	  if ( first ) 
	    O2 << "NONE";
	  O2 << "\n";
	}
      
      // Consider next SNP to tag
      ++i;
    }
  
  if ( par::gettag_listall ) 
    O2.close();
  

  if ( ! testAll )
    {
      printLOG("In total, added " 
	       + int2str( tagged.size() - toTag.size() )
	       + " tag SNPs\n");
      
      ofstream O1;
      printLOG("Writing tag list to [ " + par::output_file_name + ".tags ]\n");
      O1.open( (par::output_file_name+".tags").c_str(), ios::out);

      // Mode 1 : just write list
      if ( par::gettag_mode1 )
	{
	  set<int>::iterator l = tagged.begin();
	  while ( l != tagged.end() )
	    {
	      O1 << locus[*l]->name << "\n";
	      ++l;
	    }
	}
      
      // Mode 2 : Write all SNPs, with 0/1 code 
      if ( par::gettag_mode2 )
	{
	  for (int l = 0; l < nl_all ; l++ )
	    {
	      O1 << locus[l]->name << "\t";
	      if ( tagged.find(l) != tagged.end() )
		O1 << "1\n";
	      else
		O1 << "0\n";	  
	    }
	}
      
      O1.close();
    }

}
