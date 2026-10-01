

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
#include <cmath>
#include <vector>
#include <map>
#include <cstdlib>

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"
#include "nlist.h"
#include "stats.h"

extern ofstream LOG;

using namespace std;

// Format for haplotype file: format 1:
// SNP_ID  CHR  CM  BP  A1 A2   HAP   SNP1  SNP2 ... 

// i.e. length of Pred_allele indicate how many SNPs to expect:
// rs10001 5 0 10203  A G   TTC  rs001 rs002 rs003

// Alternatively: using wild cards: format 2:
// will name haplotype "H1_TTC_", "H1_CTT_", "H2_AA_", etc
// * rs001 rs002 rs003
// * rs002 rs004

// Alternatively: using wild cards: format 3:
// will name haplotype "MYHAP1_1_TTC_", "MYHAP1_1_CTT_", "GENEB_2_AA_", etc
// ** MYHAP1 rs001 rs002 rs003
// ** GENEB rs002 rs004

// Alternatively: using --whap weighted multimarker tests
// rs10001 5 0 10203  A G   rs001 rs002 rs003 / TTC / CCT 0.9 / TCT 0.1 
// i.e. "/" separator used, if number omitted, then assume PP=1

void HaploPhase::readTagFile()
{

  P.printLOG("\n");

  ///////////////////////////////////////////////
  // Haplotype inference is a SNP-major function
  
  if (!par::SNP_major) P.Ind2SNP();


  ////////////////////////
  // Lookup table for SNPs

  map<string,int> mlocus;
  map<string,int>::iterator ilocus;
  
  for (int l=0;l<P.nl_all;l++)
    mlocus.insert(make_pair(P.locus[l]->name,l));


  ///////////////////
  // Affection coding

  if (par::bt) affCoding(P);

  
  ////////////////////////////////
  // Read list of tags/haplotypes

  checkFileExists(par::tagfile);
  ifstream TAG(par::tagfile.c_str(), ios::in);
  TAG.clear();
  
  string f2 = par::output_file_name + ".mishap";
  ofstream MISHAP(f2.c_str(), ios::out);
  bool all_okay = true;

  // Count of new haplotypes we want to infer
  int hc=1; 
  
  while(!TAG.eof())
    {
      char c[500000];
      TAG.getline(c,500000,'\n');
      
      string l = c;

      // Catch blank lines or DOS carriage-returns
      if (l=="" || l=="\r") continue;

      // Tokenize line
      string buf; 
      stringstream ss(l); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      // whitepsace line?
      if (tokens.size() == 0)
	continue;
      
      // tokens[0]  predicted SNP rs#
      // tokens[1]  predicted SNP chromosome
      // tokens[2]  predicted SNP Morgan position
      // tokens[3]  predicted SNP base pair position
      // tokens[4]  predicted SNP allele 1
      // tokens[5]  predicted SNP allele 2
      // tokens[6]  tag allele (-> allele 1)
      // tokens[7+] predictor rs#(s)


      // If we see one or more wildcard specifications, then
      // automatically set phase_all_haps to be true;

      if ( tokens[0] == "*" )
	{ 
	  // Wildcard format

	  if ( tokens.size() == 1 )
	    error("Problem with " + 
		  par::tagfile + " line\n" + l 
		  + "\n: must have atleast one SNP list\n");

	  par::phase_hap_all = true; 
	}
      else if ( tokens[0] == "**")
	{ 
	  // Wildcard2 format
	  
	  if (tokens.size() == 2 )
	    error("Problem with " 
		  + par::tagfile + " line\n" + l 
		  + "\n: must have atleast one SNP list\n");

	  par::phase_hap_all = true; 

	}
      else if ( ! par::weighted_mm )
	{
	  // Standard format

	  if (tokens.size() < 8 ) 
	    {
	      string e = "Problem with " + 
		par::tagfile + " line\n" 
		+ l + 
		"\n (expecting at least 8 items, or to start with */** wildcard)\n";
	      error(e);	  
	    }
	} 
      else
	{
	  // Weighted multi-marker format
	  
	  
	  if (tokens.size() < 9 ) 
	    {
	      string e = "Problem with " + par::tagfile + " line\n" 
		+ l 
		+ "\n (expecting at least 9 items for --whap format file)\n";
	      error(e);	  
	    }


	}
 
      if (tokens[0].substr(tokens[0].size()-1) == "_" ) 
	error("Cannot use '_' in tag/haplotype name: reserved for wildcards\n");


      int len; 		        // length of haplotype
      vector<int> locusList;	// list of predictor #s
      
    
      // Is this particular line a wildcard? okay?
      bool wildcard = tokens[0] == "*" || tokens[0] == "**" ? true : false ;
      string wildname = "H";
      
      // Take the name from the second position? (** wildcard?) 
      if ( tokens[0] == "**" )
	{
	  wildname = tokens[1]+"_";
	  tokens.erase(tokens.begin()+1);
	}

      bool okay = true;




      /////////////////////////////////
      // Fully-specified haplotype 
      
      if (!wildcard) 
	{
	  
	  int offset = 7;
	  if ( par::weighted_mm ) offset = 6;

	  if ( par::weighted_mm )

	    {
	      // Find first "/" separator
	      int sep = 6;
	      while ( tokens[++sep] != "/" ) { }
	      len = sep - 6;
	    }
	  else
	    {
	      len = tokens[6].length();
	      
	      if (len != tokens.size() - offset )
		{
		  string e = "Problem with " + 
		    par::tagfile + " line\n" + l + "\n";
		  error(e);
		}      
	    }

	  
	  //////////////////////
	  // Lookup locus name
	  
	  for (int i=0; i<len; i++)
	    {
	      ilocus = mlocus.find(tokens[i+offset]);
	      if (ilocus != mlocus.end())
		{
		  locusList.push_back(ilocus->second);	  	  
		}
	      else 
		{
		  MISHAP << "NOSNP\t" << tokens[0] 
			 << "\t" << tokens[i+offset] << "\n";
		  okay = false;
		}
	    }
	  

	  /////////////////////////////////
	  // Check specified alleles exist	  

	  if (okay)
	    {

	      // Just one haplotype to check for non-weighted test version
	      
	      if ( ! par::weighted_mm )
		{
		  for (int s=0;s<len;s++)
		    if ( ! ( P.locus[locusList[s]]->allele1 
			     == tokens[6].substr(s,1) 
			     || P.locus[locusList[s]]->allele2 
			     == tokens[6].substr(s,1) ) )
		      {
			MISHAP << "NOALLELE\t" << tokens[0] 
			       << "\t" << P.locus[locusList[s]]->name
			       << "\t" << tokens[6] << "\n";
			okay = false;
		      }
		}

	      // Otherwise, we must parse through and check each
	      else
		{
		  
		  int allele = 5 + len + 2;

		  for ( int i = allele ; i < tokens.size() ; i++ ) 
		    {
		      // Is this a haplotype?
		      if ( i == allele || tokens[i-1] == "/" )
			{
			  for (int s=0;s<len;s++)
			    if ( ! ( P.locus[locusList[s]]->allele1 
				     == tokens[i].substr(s,1) 
				     || P.locus[locusList[s]]->allele2 
				     == tokens[i].substr(s,1) ) )
			      {
				MISHAP << "NOALLELE\t" << tokens[0] 
				       << "\t" << P.locus[locusList[s]]->name
				       << "\t" << tokens[i] << "\n";
				okay = false;
			      }
			}
		    }
		}

	    }

	  
	}
      else
	{
	  /////////////////////////////////
	  // Wildcard selection

          len = tokens.size()-1;

          // Lookup locus name
          for (int i=0; i<len; i++)
	    {
	      ilocus = mlocus.find(tokens[i+1]);
	      if (ilocus != mlocus.end())
		{
		  locusList.push_back(ilocus->second);	  	  
		}
	      else 
		{
		  MISHAP << "NOSNP\t" << tokens[i+1] << "\n";
		  okay = false;
		}
	    }
	}
      


      //////////////////////////////////////
      // Should we try to add this haploype
           
      if (!okay) 
	{
	  all_okay = false;
	  continue;
	}


      ////////////////////////////////////////////////////////
      // Check that all predictors are on the same chromosome
      
      if (!wildcard)  
	{    
	  for (int ck=0; ck<len; ck++)
	    if (P.locus[locusList[ck]]->chr != atoi(tokens[1].c_str()))
	      {
		MISHAP << "DIFF_CHR\t" << P.locus[locusList[ck]]->name 
		       << "\t" << tokens[0] << "\n";	  
		okay = false;
	      }
	}
      else
	{	
	  for (int ck=0; ck<len-1; ck++)
	    if (P.locus[locusList[ck]]->chr != P.locus[locusList[ck+1]]->chr)
	      {
		MISHAP << "DIFF_CHR\t" << P.locus[locusList[ck]]->name 
		       << "\t" << P.locus[locusList[ck+1]]->name << "\n";
		okay = false;
	      }
	}
      
      if (!okay) 
	{
	  all_okay = false;
	  continue;
	}
      
      
      ///////////////////////////////////////
      // Standard approach -- only one entry
	
      if ( ( !wildcard) && (!par::weighted_mm) )
	{
	  
	  // Add numbers for predictors to list
	  new_pred_locus.push_back(locusList);
	  
	  // Add which allele to look for (corresponding to allele1)
	  new_pred_allele.push_back(tokens[6]);
	  
	  // Make new entry in MAP file
	  Locus * loc = new Locus;
	  loc->name = tokens[0];
	  loc->chr = getChromosomeCode( tokens[1] );
	  loc->pos = atof(tokens[2].c_str());
	  loc->bp = atoi(tokens[3].c_str());  
	  loc->allele1 = tokens[4];
	  loc->allele2 = tokens[5];
	  
	  // Add this new locus to the list
	  new_map.push_back(loc);
	}
      else if ( ( !wildcard ) && par::weighted_mm )
	{

	  ///////////////////////////////////////
	  // Weighted MM test -- only one entry

  	  // Add numbers for predictors to list
	  new_pred_locus.push_back(locusList);
	  
	  // Add which allele(s) to look for (corresponding to allele 1)
	  // and then the corresponding weights
	  
	  map<string,double> whap;

	  int index = 5 + len + 2;
	  while ( index < tokens.size() ) 
	    {
	      // Read haplotype, then optionally a weight
	      
	      if ( tokens[index].size() != len )
		error("Problem with " + par::tagfile + " line\n" + l + "\n");
	      
	      // Read weight, or advance to next haplotype?
	      if ( index == tokens.size() - 1 )
		{
		  whap.insert(make_pair( tokens[index] , 1 ) );
		  break;
		}
	      else if ( tokens[index+1] == "/" )
		{
		  whap.insert(make_pair( tokens[index] , 1 ) );
		  index += 2;
		}
	      else
		{
		  double w = atof(tokens[index+1].c_str());
		  if ( ( !realnum(w) ) ||
		       w < 0 || w > 1 ) 
		    error("Problem with specified weight in line:\n"+l+"\n");
		  whap.insert(make_pair( tokens[index] , w ) );
		  index += 3;
		}
	    }

	  
	  new_pred_weighted_allele.push_back(whap);
	  
	  // Make new entry in MAP file
	  Locus * loc = new Locus;
	  loc->name = tokens[0];
	  loc->chr = getChromosomeCode( tokens[1] );
	  loc->pos = atof(tokens[2].c_str());
	  loc->bp = atoi(tokens[3].c_str());  
	  loc->allele1 = tokens[4];
	  loc->allele2 = tokens[5];
	  
	  // Add this new locus to the list
	  new_map.push_back(loc);


	  
	}
      else
	{

	  ////////////////////////////////////////////
	  // Wildcard approach -- just a single entry

	  // Add numbers for predictors to list
	  new_pred_locus.push_back(locusList);
	      
	  // Put in a dummy allele code (we ignore this...)
	  string hstr="";
	  for (int s=0;s<S.size();s++)
	    hstr += P.locus[locusList[s]]->allele1;
		  
	  new_pred_allele.push_back(hstr);
		
	  // Make new entry in MAP file
	  
	  Locus * loc = new Locus;
	  loc->name = wildname +int2str(hc)+"_DUMMY_"+hstr+"_";
	  loc->chr = P.locus[locusList[0]]->chr;
	  loc->pos = 0;
	  loc->bp = hc; 
	  loc->allele1 = "1";
	  loc->allele2 = "2";
	  
	  // Add this new locus to the list
	  new_map.push_back(loc);
	  
	   	  
	}

      // Increment new haplotype count
      hc++;
      
      // Read next in TAG file

    }

  TAG.close();
  MISHAP.close();

  if (!all_okay)
    P.printLOG("Warning: misspecified haplotypes found: listed in [ " + f2 + " ]\n");

  
  // End of reading haplotype list -- did we encounter any problems?

  P.printLOG("Read " + int2str(new_map.size()) + 
	   " haplotypes from [ " + par::tagfile + " ]\n");  

}



void HaploPhase::makeSlidingWindow(string winspec)
{
  
  P.printLOG("\n");
  
  ///////////////////////////////////////////////
  // Haplotype inference is a SNP-major function
  
  if (!par::SNP_major) 
    P.Ind2SNP();

  
  /////////////////////////////////////
  //

  NList nl(0);
  vector<string> tok = nl.deparseStringList( winspec ); 
  vector<int> spec; // size
  vector<int> spec2; // step
  for (int i=0; i<tok.size(); i++)
    {
      if ( tok[i].find("+") == string::npos )
	{
	  int t;
	  if ( ! from_string<int>( t, tok[i], std::dec ) )
	    error("Problem with specification of haplotype sliding window");
	  spec.push_back(t);
	  spec2.push_back(1);
	}
      else
	{
	  string u1 = tok[i].substr(0,tok[i].find("+"));
	  string u2 = tok[i].substr(tok[i].find("+")+1);
		      
	  int t;
	  if ( ! from_string<int>( t, u1, std::dec ) )
	    error("Problem with specification of haplotype sliding window");
	  spec.push_back(t);
	  if ( ! from_string<int>( t, u2 , std::dec ) )
	    error("Problem with specification of haplotype sliding window");
	  spec2.push_back(t);
	}
    }

  int w=1;
      
  for (int i=0; i<spec.size(); i++)
    {
      
      int winsize = spec[i];
      int winstep = spec2[i];
      
      //////////////////////////////////
      // Set up haplotype entries
      
      int start = 0;
      
      while ( 1 )
	{
	  
	  // Beyond last SNP?
	  
	  if ( start >= P.nl_all ) 
	    break;
	  
	  // Make a window, as large as possible
	  vector<int> snps;
	  
	  // Make sure it is restricted to one chromosome
	  int chr = P.locus[start]->chr;
	  
	  bool fail = false;
	  bool newChromosome = false;
	  int actualStop = start;
	  
	  // Add SNPs to window
	  
	  for (int s = start; s < start + winsize; s++)
	    {
	      
	      // No more SNPs left 
	      if ( s == P.nl_all  ) 
		{
		  fail = true;
		  break;
		}
	      
	      // Next chromosome?
	      if ( P.locus[s]->chr != chr )
		{	      
		  newChromosome = true;
		  actualStop = s;
		  break;
		}
	      
	      snps.push_back(s);	  
	      actualStop = s;
	    }
	  
	  // Have we come up to the end of this chromosome in an exact number?
	  if ( actualStop == P.nl_all-1 )
	    {
	      fail = true;
	    }
	  else if ( ( ! newChromosome ) && P.locus[actualStop+1]->chr != chr )
	    {
	      newChromosome = true;
	      actualStop++;
	    }
	  
	  
	  // Finished constructing this particular window: 
	  // do we have anything to add?
	  
	  if ( snps.size() == 0 )
	    {
	      if ( newChromosome )
		start = actualStop;
	      else
		start += winstep;
	      continue;
	    }

	  Locus * tmploc = new Locus;
	  tmploc->name = "WIN"+int2str(w++);
	  tmploc->chr = P.locus[start]->chr;
	  tmploc->pos = P.locus[start]->pos;
	  tmploc->bp = P.locus[start]->bp;
	  tmploc->allele1 = "1";
	  tmploc->allele2 = "2";
	  new_pred_locus.push_back(snps);
	  
	  new_map.push_back(tmploc);
	  new_pred_allele.push_back("");
	  
	  
	  // Advance window
	  
	  if ( fail )
	    break;
	  
	  if ( newChromosome )
	    start = actualStop;
	  else
	    start += winstep;
	  
	}
           
    }
  
  P.printLOG("Created " + int2str(w-1) + " sliding windows\n");

}



void HaploPhase::setSpecificSNPs(string snps)
{

  map<string,int> mapping;
  for (int l=0; l<P.nl_all;l++)
    mapping.insert(make_pair(P.locus[l]->name,l));
  
  NList nl(P.nl_all,true);   
  vector<int> snplist = nl.deparseStringList(snps,&mapping);
  
  if (snplist.size() == 0 )
    return;
  
  int start = snplist[0];

  new_pred_locus.push_back(snplist);
  
  Locus * tmploc = new Locus;
  tmploc->name = "WIN1";
  tmploc->chr = P.locus[start]->chr;
  tmploc->pos = P.locus[start]->pos;
  tmploc->bp = P.locus[start]->bp;
  tmploc->allele1 = "1";
  tmploc->allele2 = "2";	      
  new_map.push_back(tmploc);
  
  new_pred_allele.push_back("");		      
  
  return;
}
