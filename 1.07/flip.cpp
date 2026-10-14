

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

void Plink::flipStrand()
{

  ////////////////////////////
  // Look-up table by SNP name

  map<string,Locus*> mlocus;
  map<string,Locus*>::iterator ilocus;
  map<string,int> mnlocus;
  map<string,int>::iterator inlocus;

  vector<Locus*>::iterator loc = locus.begin();
  while ( loc != locus.end() )
    {
      mlocus.insert(make_pair( (*loc)->name, *loc) );
      loc++;
    }  


  ////////////////////////////
  // Look-up table by SNP name

  map<string,Individual*> mpeople;
  
  
  //////////////////////////////////////////////////////////
  // Performing this for all individuals, or just a subset?
  set<Individual*> pflip;
  
  if ( par::flip_subset )
    {

      if ( ! par::SNP_major ) 
	Ind2SNP();
      
      for (int i=0; i<n; i++)
	{
	  Individual * person = sample[i];
	  string id = person->fid + "_" + person->iid;
	  mpeople.insert( make_pair( id , person ) );
	}

      for (int l=0; l<nl_all; l++)
	{
	  mnlocus.insert( make_pair( locus[l]->name , l ) );
	}
      
      checkFileExists(par::flip_subset_file);  
      printLOG("Reading individuals to flip strand for [ " 
	       + par::flip_subset_file + " ] \n");
      ifstream IN2(par::flip_subset_file.c_str(),ios::in);
      int pcount1 = 0;
      int pcount2 = 0;

      while ( ! IN2.eof() )
	{
	  
	  string fid, iid;
	  IN2 >> fid >> iid;
	  string id = fid + "_" + iid;
	  if ( fid == "" ) 
	    continue;
	  ++pcount1;
	  map<string,Individual*>::iterator ipeople = mpeople.find( id );
	  if ( ipeople == mpeople.end() )
	    continue;
	  pflip.insert( ipeople->second );
	  ++pcount2;
	}
      printLOG("Read " + int2str(pcount1) + " individuals, of whom " 
	       + int2str(pcount2) + " were found, to flip\n");
      IN2.close();
    }
  

  ///////////////////////
  // Read in SNPs to flip

  checkFileExists(par::flip_file);  
  printLOG("Reading SNPs to flip strand from [ " 
	   + par::flip_file + " ] \n");
  ifstream INFILE(par::flip_file.c_str(),ios::in);
  INFILE.clear();
  int counter = 0;

  while (!INFILE.eof())
    {
      string m;
      INFILE >> m;
      if (m=="") continue;

      if ( par::flip_subset )
	{
	  inlocus = mnlocus.find(m);
	  if (inlocus != mnlocus.end() )
	    {
	      ++counter;
	      for (int i=0; i<n; i++)
		{
		  Individual * person = sample[i];
		  if ( pflip.find( person ) != pflip.end() )
		    {
		      int l = inlocus->second;
		      if ( SNP[ l ]->one[ i ] == 
			   SNP[ l ]->two[ i ] )
			{
			  SNP[ l ]->one[ i ] = ! SNP[ l ]->one[ i ];
			  SNP[ l ]->two[ i ] = ! SNP[ l ]->two[ i ];
			}
		    }		  
		}
	    }
	}
      else
	{
	  
	  ilocus = mlocus.find(m);
	  if (ilocus != mlocus.end())
	    {
	      counter++;
	      
	      // Flip strand
	      
	      Locus * loc = ilocus->second;
	      
	      if ( loc->allele1 == "A" ) loc->allele1 = "T";
	      else if ( loc->allele1 == "C" ) loc->allele1 = "G";
	      else if ( loc->allele1 == "G" ) loc->allele1 = "C";
	      else if ( loc->allele1 == "T" ) loc->allele1 = "A";
	      else if ( loc->allele1 == "1" ) loc->allele1 = "4";
	      else if ( loc->allele1 == "2" ) loc->allele1 = "3";
	      else if ( loc->allele1 == "3" ) loc->allele1 = "2";
	      else if ( loc->allele1 == "4" ) loc->allele1 = "1";
	      
	      if ( loc->allele2 == "A" ) loc->allele2 = "T";
	      else if ( loc->allele2 == "C" ) loc->allele2 = "G";
	      else if ( loc->allele2 == "G" ) loc->allele2 = "C";
	      else if ( loc->allele2 == "T" ) loc->allele2 = "A";
	      else if ( loc->allele2 == "1" ) loc->allele2 = "4";
	      else if ( loc->allele2 == "2" ) loc->allele2 = "3";
	      else if ( loc->allele2 == "3" ) loc->allele2 = "2";
	      else if ( loc->allele2 == "4" ) loc->allele2 = "1";
	      
	    }
	}
      // Next SNP
    }

  INFILE.close();
  
  printLOG("Flipped strand of " + int2str(counter) + " SNPs\n");
  
}



void Plink::calcFlipScan()
{
  
  ///////////////////////////////////////////
  // Screen for SNPs that appear to have been
  // flipped in one dataset versus another, 
  // based on patterns of LD.  Just use a simple
  // sliding window, keeping track of the number
  // and magnitude of concordant versus discordant
  // LD pairs

  
  

  ofstream OUT1;
  string f = par::output_file_name + ".flipscan";
  OUT1.open(f.c_str(),ios::out);
  OUT1.precision(3);

  printLOG("Writing FS statistics to [ " + f + " ] \n");

  
  OUT1 << setw(6) << "CHR" << " "
     << setw(par::pp_maxsnp) << "SNP"  << " "
     << setw(12) << "BP" << " "
     << setw(4) << "A1" << " "
     << setw(4) << "A2" << " "
     << setw(8) << "F"  << " "
     << setw(6) << "POS" << " "
     << setw(8) << "R_POS" << " "
     << setw(6) << "NEG" << " "
     << setw(8) << "R_NEG" << " "
     << "NEGSNPS\n";

  ofstream OUT1V;
  if ( par::flip_scan_verbose )
    {
      string f = par::output_file_name + ".flipscan.verbose";
      OUT1V.open(f.c_str(),ios::out);
      OUT1V.precision(3);
      printLOG("Writing FS verbose output to [ " + f + " ] \n");
      
      OUT1V << setw(6) << "CHR_INDX" << " "
	  << setw(par::pp_maxsnp) << "SNP_INDX"  << " "
	  << setw(12) << "BP_INDX" << " "
	  << setw(4) << "A1_INDX" << " "	
	  << setw(par::pp_maxsnp) << "SNP_PAIR"  << " "
	  << setw(12) << "BP_PAIR" << " "
	  << setw(4) << "A1_PAIR" << " "
	  << setw(8) << "R_A" << " "
	  << setw(8) << "R_U" << "\n";
	
    }

  
  ///////////////////////////
  // Index locus  
  
  for (int l1=0; l1<nl_all; l1++)  
    {

      int cntPlus = 0;
      int cntNeg = 0;

      double scorePlus = 0;
      double scoreNeg = 0;

      stringstream verbose_buffer;
      verbose_buffer.precision(3);
      
      string snplist = "";

      // Second locus
	  
      for (int l2=0; l2<nl_all; l2++)      
	{
	  
	  // No point in testing the same SNP

	  if ( l1 == l2 ) 
	    continue;
	  
	  // Outside of window?
	  
	  if ( l2 - l1 >= par::disp_r_window_snp )
	    continue;
	  
	  if ( l1 - l2 >= par::disp_r_window_snp )
	    continue;
	  
	  if ( locus[l2]->chr != locus[l1]->chr ) 
	    continue;
	  
	  if ( locus[l2]->bp - locus[l1]->bp 
	       > par::disp_r_window_kb )
	    continue;
	  
	  if ( locus[l1]->bp - locus[l2]->bp 
	       > par::disp_r_window_kb )
		continue;
	  
	  
	  //////////////////////////////////////
	  // Calculate correlation (un-squared)
	  // in cases and controls separately

	  setFlagToCase();
	  double rCase = correlation2SNP(l1,l2,false,false,true);
	  
	  setFlagToControl();
	  double rControl = correlation2SNP(l1,l2,false,false,true);

	  // Keep track of score
	  
	  bool sameDirection = true;
	  if ( ( rCase > 0 && rControl < 0 ) || 
	       ( rControl > 0 && rCase < 0 ) )
	    sameDirection = false;
	  
	  bool aboveThreshold = true;
	  
	  if ( ( rCase > - par::flip_scan_threshold && rCase < par::flip_scan_threshold ) &&
	       ( rControl > - par::flip_scan_threshold && rControl < par::flip_scan_threshold ) )
	    aboveThreshold = false;

	  if ( ! realnum( rCase ) ) 
	    aboveThreshold = false;

	  if ( ! realnum( rControl ) ) 
	    aboveThreshold = false;
	  
	  if ( aboveThreshold ) 
	    {
	      if ( sameDirection ) 
		{
		  ++cntPlus;
		  if ( rCase > 0 ) 
		    scorePlus += ( rCase + rControl ) / 2.0;
		  else
		    scorePlus += ( -rCase - rControl ) / 2.0;
		}
	      else
		{
		  ++cntNeg;
		  if ( rCase > 0 ) 
		    scoreNeg += ( rCase - rControl ) / 2.0;
		  else
		    scoreNeg += ( -rCase + rControl ) / 2.0;

		  if ( snplist == "" ) 
		    snplist = locus[l2]->name;
		  else
		    snplist += "|" + locus[l2]->name;
		  
		}
	   
	      if ( par::flip_scan_verbose )
		{
		  verbose_buffer << setw(6) << chromosomeName( locus[l1]->chr ) << " "
				 << setw(par::pp_maxsnp) << locus[l1]->name  << " "
				 << setw(12) << locus[l1]->bp << " "
				 << setw(4) << locus[l1]->allele1 << " "	
				 << setw(par::pp_maxsnp) << locus[l2]->name  << " "
				 << setw(12) << locus[l2]->bp << " "
				 << setw(4) << locus[l2]->allele1 << " "
				 << setw(8) << rCase << " "
				 << setw(8) << rControl << "\n";
		  
		}

	    }
	  
   
	 	  
	  
	  // Consider the paired SNP
	}
      
      
      // Calculate FS score for this index SNP

      
      // Display
      
      OUT1 << setw(6) << locus[l1]->chr  << " " 
	 << setw(par::pp_maxsnp) << locus[l1]->name  << " " 
	 << setw(12) << locus[l1]->bp << " "
	 << setw(4) << locus[l1]->allele1 << " "
	 << setw(4) << locus[l1]->allele2 << " "
	 << setw(8) << locus[l1]->freq  << " "
	 << setw(6) << cntPlus << " ";

      if ( cntPlus == 0 ) 
	OUT1 << setw(8) << "NA" << " ";
      else
	OUT1 << setw(8) << scorePlus/(double)cntPlus << " ";
	  
      OUT1 << setw(6) << cntNeg << " ";

      if ( cntNeg == 0 ) 
	OUT1 << setw(8) << "NA" << " "; 
      else
	OUT1 << setw(8) << scoreNeg/(double)cntNeg << " "; 
      
      OUT1 << snplist << "\n";

      if ( par::flip_scan_verbose )
	{
	  if ( cntNeg > 0 ) 
	    OUT1V << verbose_buffer.str();
	}


    } // Consider the next index SNP
  
  OUT1.close();
      
  if ( par::flip_scan_verbose )
    OUT1V.close();
  
}



void Plink::setReferenceAllele()
{

  ////////////////////////////
  // Look-up table by SNP name

  map<string,int> mnlocus;
  map<string,int>::iterator inlocus;

  for (int l=0; l<nl_all; l++)
    mnlocus.insert(make_pair( locus[l]->name, l) );
  

  ////////////////
  // Read in SNPs 

  checkFileExists( par::set_reference_allele_file ); 
  printLOG("Reading SNPs to set reference allele [ "
	   + par::set_reference_allele_file + " ] \n");
  ifstream INFILE(par::set_reference_allele_file.c_str(),ios::in);
  INFILE.clear();

  int counter_changed = 0;
  int counter_problem = 0;
  int counter_unchanged = 0;
  
  while (!INFILE.eof())
    {
      string m, ref;
      INFILE >> m >> ref;

      if (m=="") continue;
      
      inlocus = mnlocus.find(m);
      
      if (inlocus != mnlocus.end() )
	{
	  
	  int l = inlocus->second;
	  
	  // Check alleles
	  if ( ref == locus[l]->allele1 )
	    {
	      ++counter_unchanged;
	      continue;
	    }
	  
	  if ( ref != locus[l]->allele2 )
	    {
	      ++counter_problem;
	      continue;
	    }
	  

	  ++counter_changed;

	  // Swap alleles, recode genotypes internally

	  string tmp = locus[l]->allele2;
	  locus[l]->allele2 = locus[l]->allele1;
	  locus[l]->allele1 = tmp;

// 	  cout << "changed " << m << " ref= " << ref << " now = "  
// 	       << locus[l]->allele1 << " " << locus[l]->allele2 << "\n";

	  for (int i=0; i<n; i++)
	    {
	      
	      if ( SNP[ l ]->one[ i ] == 
		   SNP[ l ]->two[ i ] )
		{
		  SNP[ l ]->one[ i ] = ! SNP[ l ]->one[ i ];
		  SNP[ l ]->two[ i ] = ! SNP[ l ]->two[ i ];
		}
	      
	    }
	}
      
      
      // Next SNP
    }

  INFILE.close();
  
  printLOG("Set reference alleles for " + int2str(counter_changed+counter_unchanged) + " SNPs, ");
  printLOG(int2str(counter_changed) + " different from minor allele\n");
  if (counter_problem>0)
    printLOG("Also, " + int2str(counter_problem) + " couldn't be changed due to bad allele codes\n");
  
}

