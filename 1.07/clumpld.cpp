

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

#include "clumpld.h"
#include "phase.h"
#include "options.h"
#include "plink.h"


//////////////////////////////////////////////////
// Set user-defined parameters within constructor

clump_LD::clump_LD(Plink * pp,
		   HaploPhase * hp_,
		   double sig, 
		   double dist, 
		   double secondp, 
		   float r2c)
{
  P = pp;
  hp = hp_;
  pval_cutoff = sig;
  ld_distance = dist;
  second_pval_cutoff = secondp;
  r2_cutoff = r2c;  
  
}

/////////////
// helper functions

string returnFullRangeList(Range & r1, map<string, set<Range> > & ranges, bool verbose);


/////////////
// accessors

void clump_LD::set_pval( double sig ){
  pval_cutoff = sig;
}

void clump_LD::set_second_pval( double secondp ){
  second_pval_cutoff = secondp;
}

void clump_LD::set_ld( double dist ){
  ld_distance = dist;
}

void clump_LD::set_r2( double r2c ){
  r2_cutoff = r2c;
}


//////////////////////////////////////////////////
// read association file and pull out significant 
// results sort by pvalue

vector<ResultTrio> clump_LD::read_assoc_file(string fileList)
{
  
  vector<ResultTrio> sp;

  // We may be reading a single file, or more than one file
  fileList = searchAndReplace(fileList,","," ");
  
  // Tokenize
  filename.clear();
  string buf; 
  stringstream ss(fileList);
  while (ss >> buf)
    filename.push_back(buf);
  
  // Read each file
  
  for (int f=0; f<filename.size(); f++)
    {

      checkFileExists(filename[f]);
      ifstream RESIN( filename[f].c_str(), ios::in );

      int snp_column = -1;
      int pval_column = -1;
      vector<int> annot_field;

      // Get header row
      
      char cline[par::MAX_LINE_LENGTH];
      RESIN.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      string sline = cline;
      if (sline=="") 
	error("Problem reading [ " + par::clumpld_results + " ]\n");
      
      vector<string> tok_annot;

      if ( par::clumpld_annot )
	{
	  string afields = searchAndReplace(par::clumpld_annot_fields,","," ");
  	  string buf; 
	  stringstream ss(afields);
	  while (ss >> buf)
	    tok_annot.push_back(buf);	 
	}

      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      for (int i=0; i<tokens.size(); i++)
	{
	  if ( tokens[i] == "SNP" )
	    snp_column = i;
	  
	  if ( tokens[i] == par::clumpld_column )
	    pval_column = i;	  

	  if ( par::clumpld_annot ) 
	    for ( int f=0; f<tok_annot.size(); f++)
	      if ( tokens[i] == tok_annot[f] )
		annot_field.push_back(i);
	}
      
      
      //////////////////////////////////////
      // Did we find the necessary columns?
      
      if ( snp_column < 0 || pval_column < 0 )
	{
	  error("Could not find [" + par::clumpld_column 
		+ "] as given in --clump-field\n");
	}
      
      P->printLOG("Reading results for clumping from [ "
		  + filename[f] + " ]\n");
      P->printLOG("Extracting fields SNP and " 
		  + par::clumpld_column + "\n");

      while (!RESIN.eof())
	{
	  
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
	  
	  if ( tokens.size() <= snp_column ||
	       tokens.size() <= pval_column )
	    continue;
	  
	  ResultTrio pt;

	  // Keep track of which file this is from
	  pt.f = f+1;

	  if ( ! from_string<double>( pt.p, tokens[pval_column] , std::dec))
	    continue;
	  
	  pt.s = tokens[snp_column];
	  
	  // Create clumped vector, based just on SNP name
	  
	  if( pt.p <= pval_cutoff )
	    clumped.insert(make_pair(pt.s,false));
	  
	  ClumpPair cp;
	  cp.snp = pt.s;
	  cp.f = f+1;
	  
	  ClumpResults cr;
	  cr.p = pt.p;
	  cr.annot = "";
	  if ( par::clumpld_annot )
	    {
	      for (int f=0; f<annot_field.size(); f++)
		{
		  cr.annot += tokens[annot_field[f]];
		  if ( f < annot_field.size() - 1 )
		    cr.annot += ", ";
		}
	    }

	  assoc_results.insert( make_pair( cp, cr ));
	  sp.push_back(pt);
	  
	}
      
      RESIN.close();

      // Read in next results file
    }
  
  sort(sp.begin(), sp.end());
  return sp;
}


///////////////////////////////
// Perform clumping operation

void clump_LD::clump()
{
  

  P->printLOG("\nParameters for --clump:\n");
  P->printLOG("          p-value threshold for index SNPs = " 
	      + dbl2str(pval_cutoff) + "\n");

  P->printLOG("      Physical (kb) threshold for clumping = " 
	      + dbl2str(ld_distance/1000.0) + "\n");

  P->printLOG("     LD (r-squared) threshold for clumping = " 
	      + dbl2str(r2_cutoff) + "\n");

  P->printLOG("        p-value threshold for clumped SNPs = " 
	      + dbl2str(second_pval_cutoff) + "\n\n");


  if ( par::clumpld_annot )
    P->printLOG("Including annotation fields: " + par::clumpld_annot_fields +"\n");

  vector<ResultTrio> sp = read_assoc_file( par::clumpld_results );

  if ( par::clumpld_index1 )
    P->printLOG("Indexing only on [ " + filename[0] + " ]\n");
  else
    P->printLOG("Indexing on all files\n");

  if ( par::clumpld_only_show_replications )
    P->printLOG("Only showing cross-file clumps\n");

  if ( par::clumpld_only_show_replications_list )
    P->printLOG("Only showing non-index clumped SNPs\n");

  int zero, one, two, three, four;
  map<string,int> mlocus;
  map<int2,double> grouped_snps;

  ofstream CLMP;
  CLMP.open( (par::output_file_name + ".clumped").c_str() , ios::out);
  CLMP.precision(3);

  P->printLOG("Writing clumped results file to [ " + 
	      par::output_file_name + ".clumped ]\n");


  ofstream CLMP2;
  if ( par::clumpld_range_annotate && ! par::clumpld_verbose )
    {
      P->printLOG("Writing clumped ranges file to [ " + 
		  par::output_file_name + ".clumped.ranges ]\n");
      CLMP2.open( (par::output_file_name + ".clumped.ranges").c_str() , ios::out);      
      CLMP2 << setw(4) << "CHR" << " " 
	    << setw(par::pp_maxsnp) << "SNP" << " " 
	    << setw(10) << "P" << " "
	    << setw(6) << "N" << " "
	    << setw(28) << "POS" << " " 
	    << setw(10) << "KB" << " "
	    << "RANGES" << "\n";

    }

  ofstream BEST;
  if ( par::clumpld_best )
    {
      BEST.open( (par::output_file_name + ".clumped.best").c_str() , ios::out);
      BEST.precision(3);
      P->printLOG("Writing best per clump to [ " 
		  + par::output_file_name 
		  + ".clumped.best" + " ]\n");
      BEST << setw(par::pp_maxsnp) << "INDEX" << " "
	   << setw(par::pp_maxsnp) << "PSNP" << " "
	   << setw(6) << "RSQ" << " "
	   << setw(8) << "KB" << " "
	   << setw(8) << "P" << " "
	   << setw(8) << "ALLELES" << " "
	   << setw(8) << "F" << "\n";
    }




  //////////////////////////
  // Read a list of ranges?
  
  map<string, set<Range> > ranges;
  map<int,set<Range*> > snp2range;

  if ( par::clumpld_range_annotate )
  {

    // Helper function to map ranges to SNPs

    mapRangesToSNPs( par::clumpld_range_file,
		     ranges,		     
		     snp2range );
      
  }


  string vmessage;

  /////////////////
  // Output header

  if ( ! par::clumpld_verbose )
    CLMP << setw(4) << "CHR" << " "
	 << setw(4) << "F" << " "	 
	 << setw(par::pp_maxsnp) << "SNP" << " " 
	 << setw(10) << "BP" << " "
	 << setw(8) << "P" << " " 
	 << setw(8) << "TOTAL" << " "
	 << setw(6) << "NSIG" << " "
	 << setw(6) << "S05" << " "
	 << setw(6) << "S01" << " "
	 << setw(6) << "S001" << " "
	 << setw(6) << "S0001" << " "
	 << setw(6) << "SP2" << "\n";
  
 

  //////////////////
  // Build SNP map
  
  for (int l=0; l<P->nl_all; l++)
    mlocus.insert(make_pair(P->locus[l]->name, l));
  
  if (sp.size()==0)
    {
      P->printLOG("No significant results given --clump parameters\n");
      return;
    }
  
  
  
  ////////////////////////////////////////////////
  // Iterate through all association results by 
  // decreasing p-value until cutoff is reached

  int clumpCount = 0;

  for (int i = 0; i < sp.size(); i++)
    {

//       if ( ! par::silent)
// 	cout << "Tested " << i << " of " << sp.size() << " SNPs                         \r";

      zero=one=two=three=four=0;
      
      ////////////////////////////////////
      // End if p-value cutoff is reached
      
      if( sp[i].p > pval_cutoff )
	break;
      

      //////////////////////////////////
      // Skip already clumped this SNP
      // (unless running in "best" mode
      // where we always want to pull the 
      // best SNP
      
      if( clumped[sp[i].s] && ! par::clumpld_best )
	{
	  continue;
	}

      /////////////////////////////////////////////
      // Are we only indexing based on first file?

      if ( par::clumpld_index1 && sp[i].f > 1 ) 
	{
	  continue;
	}
      

      ///////////////////////////////////////////
      // Record dataset of index SNP

      int indexF = sp[i].f;


      ///////////////////////////////////////////
      // Make sure associated SNP is in SNP map
      
      map<string,int>::iterator ilocus;
      ilocus = mlocus.find(sp[i].s);
      int l = -1;
      if (ilocus != mlocus.end())
	{
	  l = ilocus->second;
	}
      else 
	{
	  
	  if ( !par::clumpld_verbose ) 
	    CLMP << setw(4) << "NA" << " "
		 << setw(4) << "NA" << " "
		 << setw(par::pp_maxsnp) << sp[i].s << " "
		 << setw(10) << "NA" << " "
		 << setw(10) << sp[i].p << " " 
		 << setw(8) << "NA" << " " 
		 << setw(6) << "NA" << " " 
		 << setw(6) << "NA" << " " 
		 << setw(6) << "NA" << " " 
		 << setw(6) << "NA" << " " 
		 << setw(6) << "NA" << " "
		 << setw(6) << "NA" << "\n";
	  else
	    vmessage += sp[i].s + " not found in dataset\n";
	  
	  if ( par::clumpld_best )
	    {
	      BEST << setw(par::pp_maxsnp) << sp[i].s << " "
		   << setw(par::pp_maxsnp) << "NF" << " "
		   << setw(6) << "NF" << " "
		   << setw(8) << "NF" << " "
		   << setw(8) << "NF" << " "
		   << setw(8) << "NF" << " "
		   << setw(8) << "NF" << " ";
		BEST << "\n";
	    }
	  

	    continue;
	}


      ///////////////////////////////
      // Check all SNPs in LD range
      
      // Set at reference SNP, l, and move out
      // left and right (l1,l2)

      int l1 = l; 
      int l2 = l;
      
      // If multiple files, allow for self comparison also
      if ( filename.size()>1 ) --l1;

      set<string> willClump;

      map<int,string> inPhaseAllele;

      while(1)
	{

	  bool failed1 = false, failed2 = false;
	  
	  // Expand outwards to physical limits/SNP sets
	  
	  // Moving right

	  if( l1 < P->locus.size()-1 )
	    {
	      // Advance a position
	      l1++; 

	      double r2a = -1;
	      
	      // Compute r-squared if this SNP is close
	      // enough to the index SNP, physically
	      
	      if( P->locus[l1]->chr == P->locus[l]->chr && 
		  P->locus[l1]->bp - P->locus[l]->bp < ld_distance )
		r2a = hp->rsq( l, l1 );
	      else
		failed1 = true;


	      //////////////////////////////////
	      // Skip already clumped this SNP
	      
	      if( par::clumpld_indep && clumped[ P->locus[l1]->name ] )
 		  continue;
	      
	      // If in LD with association result
	      
	      if( r2a > r2_cutoff )
		{
		  
		  // Record which alleles are correlated

		  inPhaseAllele.insert(make_pair(l1, allelePairs(l,l1)));


		  // Record that this SNP has been clumped

		  willClump.insert(P->locus[l1]->name);
		  
		  // Now look at the assocations (in the multiple files) 
		  // with this SNP
			      
		  for (int f=1; f<=filename.size(); f++)
		    {
		      
		      // Do not clump with self
		      if ( l == l1 && f == sp[i].f )
			continue;

		      // Are we requiring that only cross-file 
		      // (i.e. no index) SNPs are listed?
		      if ( par::clumpld_only_show_replications_list && f == indexF ) 
			continue;
		      

		      ClumpPair result;
		      result.snp = P->locus[l1]->name;
		      result.f = f;

		      // Result not found
		      if ( assoc_results.find(result) == assoc_results.end() )
			continue;
		      
		      double pval = assoc_results[result].p;
 
		      if( pval < second_pval_cutoff )
			{
			  int2 t(l1,f);
			  grouped_snps.insert(make_pair(t,r2a));
			}

		      if( pval < .0001 )
			four++;
		      else
			if( pval < .001 )
			  three++;
			else
			  if( pval < .01 )
			    two++;
			  else
			    if( pval < .05 )
			      one++;
			    else
			      zero++;
		    }
		}
	    }
	  else 
	    failed1 = true;
	      


	  //////////////
	  // Move left

	  if( l2 > 0 )
	    {
	      l2--;
	      
	      double r2b = -1;
	      

	      if( P->locus[l2]->chr == P->locus[l]->chr && 
		  P->locus[l]->bp - P->locus[l2]->bp < ld_distance )
		r2b = hp->rsq( l, l2 );
	      else
		failed2 = true;
	      

	      //////////////////////////////////
	      // Skip already clumped this SNP
	      
	      if( par::clumpld_indep && clumped[ P->locus[l2]->name ] )
  		  continue;


	      /////////////////////////////////////
	      // Does this SNP meet r-sq threshold?
	      
	      if( r2b > r2_cutoff )
		{
		  
		  // Find the allele in phase with rare allele for l
		  
		  inPhaseAllele.insert(make_pair(l2, allelePairs(l,l2)));


		  // Record that this SNP is clumped

		  willClump.insert(P->locus[l2]->name);
		  
		  // Now look at associations

		  for (int f=1; f<=filename.size(); f++)
		    {


		      // Are we requiring that only cross-file 
		      // (i.e. no index) SNPs are listed?
		      if ( par::clumpld_only_show_replications_list && f == indexF ) 
			continue;
		      
		      ClumpPair result;
		      result.snp = P->locus[l2]->name;
		      result.f = f;
		      
		      // Result not found
		      if ( assoc_results.find(result) == assoc_results.end() )
			continue;

		      double pval = assoc_results[result].p;
		      
		      if( pval < second_pval_cutoff )
			{
			  int2 t(l2,f);
			  grouped_snps.insert(make_pair(t,r2b));
			}

		      if( pval < .0001 )
			four++;
		      else
			if( pval < .001 )
			  three++;
			else
			  if( pval < .01 )
			    two++;
			  else
			    if( pval < .05 )
			      one++;
			    else
			      zero++;
		    }
		} 

	    }
	  else
	    failed2 = true;
	  
	  
	  // No point in looking further?

	  if( failed1 && failed2 ) 
	    break;	  
	
	}
      

      //////////////////////////////////////////////////////////
      // Report results

      // Are we only interested in cross-file clumpings?

      if ( par::clumpld_only_show_replications )
	{
	  bool seen_replication = false;
	  map<int2,double>::iterator gi = grouped_snps.begin();	  
	  int cnt=0;
	  int lastf;
	  while( gi != grouped_snps.end() )
	    {
	      if ( cnt > 0 ) 
		{
		  if ( gi->first.p2 != lastf ) 
		    seen_replication = true;		 
		  else if ( gi->first.p2 != indexF )
		    seen_replication =  true;
		}	      
	      lastf = gi->first.p2;
	      ++cnt;
	      ++gi;
	    }

	  // If no cross-file results, then do not report this 
	  // clump -- clear all flags for it

	  if ( ! par::clumpld_best )
	    if ( ! seen_replication ) 
	      {
		grouped_snps.clear();
		continue;
	      }
	}
    
      
      // Indicate that this SNP has now been clumped
      set<string>::iterator si = willClump.begin();
      while ( si != willClump.end() )
	{
	  clumped[ *si ] = true;
	  ++si;
	}

      
      
      int total = zero+one+two+three+four;
      
      if ( par::clumpld_verbose ) 
	{
	  // Repeat header
	  CLMP << "\n"
	       << setw(4) << "CHR" << " "
	       << setw(4) << "F" << " "
	       << setw(par::pp_maxsnp) << "SNP" << " " 
	       << setw(10) << "BP" << " "
	       << setw(10) << "P" << " " 
	       << setw(8) << "TOTAL" << " "
	       << setw(6) << "NSIG" << " "
	       << setw(6) << "S05" << " "
	       << setw(6) << "S01" << " "
	       << setw(6) << "S001" << " "
	       << setw(6) << "S0001" << "\n";	  
	}

      ClumpPair cp;
      cp.snp = P->locus[l]->name;
      cp.f = sp[i].f;
      

      CLMP << setw(4) << P->locus[l]->chr << " " 
	   << setw(4) << cp.f << " " 	   
	   << setw(par::pp_maxsnp) << P->locus[l]->name << " " 
	   << setw(10) << P->locus[l]->bp << " " 
	   << setw(10) << assoc_results[cp].p << " " 
	   << setw(8) << total << " " 
	   << setw(6) << zero << " " 
	   << setw(6) << one << " " 
	   << setw(6) << two << " " 
	   << setw(6) << three << " " 
	   << setw(6) << four << " ";
      

      /////////////////////////////////////
      // Convenient format that just gives
      // the single best SNP

      if ( par::clumpld_best )
	{
	  
	  // BEST, index SNP, best SNP (or NA), r^2, KB
	  
	  int    bestSNP    = -1;
	  double bestRsq    = -1;
	  double bestP      = -1;
	  string bestAllele = "NA";
	  int    bestF      = -1;
	  string bestAnnot  = "";
	  double bestKB     = 0;
	  bool   foundSelf  = false;

	  map<int2,double>::iterator gi = grouped_snps.begin();
	  while( gi != grouped_snps.end() )
	    {
	      
	      int l0 = gi->first.p1;
	      int f = gi->first.p2;


	      // The same SNP?
	      bool isSelf = false;
	      if ( l == l0 )
		{
		  foundSelf = true;
		  isSelf = true;
		}

	      ClumpPair cp;
	      cp.snp = P->locus[l0]->name;
	      cp.f = f;
	      
	      if ( par::clumpld_only_show_replications && f == indexF )
		{
		  ++gi;
		  continue;
		}
	      
	      
	      if ( assoc_results.find(cp) != assoc_results.end() )
		{
		  if ( isSelf || ( ( ! foundSelf ) && gi->second > bestRsq ) )
		    {
		      bestRsq = gi->second;
		      
		      bestSNP = l0;		      	      
		      bestF = f;
		      bestAllele = inPhaseAllele.find(l0)->second;
		      bestP = assoc_results[cp].p;
		      if ( par::clumpld_annot )
			bestAnnot = assoc_results[cp].annot;
		      bestKB = (double)(P->locus[l0]->bp - P->locus[l]->bp)/1000.0;
		    }
		}
	      
		
	      // Consider next proxy SNP
	      ++gi;
	    }

	  
	  // Report Best SNP
	  
	  if ( bestSNP > -1 ) 
	    {
	      BEST << setw(par::pp_maxsnp) << P->locus[l]->name << " "
		   << setw(par::pp_maxsnp) << P->locus[bestSNP]->name << " ";
	      if ( foundSelf ) 
		BEST << setw(6) << "*" << " ";
	      else
		BEST << setw(6) << bestRsq << " ";

	      BEST << setw(8) << bestKB << " "
		   << setw(8) << bestP << " "
		   << setw(8) << bestAllele << " "
		   << setw(8) << bestF << " ";
	      if ( par::clumpld_annot )
		BEST << bestAnnot;
	    }
	  else
	    {
	      BEST << setw(par::pp_maxsnp) << P->locus[l]->name << " "
		   << setw(par::pp_maxsnp) << "NA" << " "
		   << setw(6) << "NA" << " "
		   << setw(8) << "NA" << " "
		   << setw(8) << "NA" << " "
		   << setw(8) << "NA" << " "
		   << setw(8) << "NA" << " ";
	      
	    }
	  BEST << "\n";
	}

      

      //////////////////////////////////
      // verbose output

      if ( par::clumpld_verbose ) 
	{

	  int minBP = P->locus[l]->bp;
	  int maxBP = P->locus[l]->bp;
	  set<string> range_notes;
	  
	  // Now list one per line, sorted by distance
	  CLMP << "\n";

	  if ( grouped_snps.size() > 0 ) 
	    {
	      
	      CLMP << "\n" 
		   << setw(4) << " " << " " 		   
		   << setw(4) << " " << " " 		   
		   << setw(par::pp_maxsnp) << " " << " "
		   << setw(10) << "KB" << " "
		   << setw(8) << "RSQ" << " "
		   << setw(8) << "ALLELES" << " "
		   << setw(4) << "F" << " "		   
		   << setw(12) << "P" << " ";
	      if ( par::clumpld_annot )
		CLMP << setw(12) << "ANNOT" << "\n";
	      else
		CLMP << "\n";
	      
	      CLMP << setw(4) << "  (INDEX) " 
		   << setw(par::pp_maxsnp) << P->locus[l]->name << " ";

	      CLMP << setw(10) << (double)(P->locus[l]->bp - P->locus[l]->bp)/1000.0 << " ";

	      CLMP << setw(8) << "1.000" << " "
		   << setw(8) << P->locus[l]->allele1 << " "
		   << setw(4) << sp[i].f << " " 		       
		   << setw(12) << assoc_results[cp].p << " ";
	      if ( par::clumpld_annot )
		CLMP << setw(12) << assoc_results[cp].annot << "\n";
	      else
		CLMP << "\n";
	      
	      CLMP << "\n";


	      //////////////////////////////////////////////////
	      // Track if SNP info already listed in this group
	      
	      set<int> infoDisplayed;
	      
	      // Are we tracking ranges?
	      
	      if ( par::clumpld_range_annotate ) 
		{
		  map<int,set<Range*> >::iterator mi = snp2range.find(l);
		  if ( mi != snp2range.end() ) 
		    {
		      set<Range*>::iterator si = mi->second.begin();
		      while ( si != mi->second.end() )
			{
			  range_notes.insert( (*si)->name );
			  ++si;
			}
		    }
		}			  
	      

	      map<int2,double>::iterator gi = grouped_snps.begin();
	      
	      while( gi != grouped_snps.end() )
		{
		  
		  int l0 = gi->first.p1;
		  int f = gi->first.p2;

		  ClumpPair cp;
		  cp.snp = P->locus[l0]->name;
		  cp.f = f;
		  
  
		  if ( assoc_results.find(cp) != assoc_results.end() )
		    {
		      if ( infoDisplayed.find(l0) == infoDisplayed.end() )
			{
			  CLMP << setw(4) << " " << " "				   
			       << setw(4) << " " << " "				   
			       << setw(par::pp_maxsnp) << P->locus[l0]->name << " "
			       << setw(10) 
			       << (double)(P->locus[l0]->bp - P->locus[l]->bp)/1000.0 
			       << " ";
			  
			  CLMP << setw(8) << gi->second << " ";
			  
			  CLMP << setw(8) << inPhaseAllele.find(l0)->second << " "
			       << setw(4) << f << " " 		       
			       << setw(12) << assoc_results[cp].p << " ";
			  if ( par::clumpld_annot )
			    CLMP << setw(12) << assoc_results[cp].annot << "\n";
			  else
			    CLMP << "\n";
			  
			  // Note that we have now seen it
			  infoDisplayed.insert(l0);
			  
			  // Track overall range of this clump 
			  if ( P->locus[l0]->bp < minBP )
			    minBP = P->locus[l0]->bp;
			  if ( P->locus[l0]->bp > maxBP )
			    maxBP = P->locus[l0]->bp;


			  ////////////////////////////////
			  // Are we also tracking ranges?
			  
			  if ( par::clumpld_range_annotate ) 
			  {
			    map<int,set<Range*> >::iterator mi = snp2range.find(l0);
			    if ( mi != snp2range.end() ) 
			      {
				set<Range*>::iterator si = mi->second.begin();
				while ( si != mi->second.end() )
				  {
				    range_notes.insert( (*si)->name );
				    ++si;
				  }
			      }
			  }			  
			}
		      else
			{
			  CLMP << setw(4) << " " << " "
			       << setw(4) << " " << " "
			       << setw(par::pp_maxsnp) << " " << " "
			       << setw(10) << " " << " "
			       << setw(8) << " " << " "
			       << setw(8) << " " << " "
			       << setw(4) << f << " " 		       
			       << setw(12) << assoc_results[cp].p << " ";
			  if ( par::clumpld_annot )
			    CLMP << setw(12) << assoc_results[cp].annot << "\n";
			  else
			    CLMP << "\n";

			}
		    }
		  		  
		  gi++;
		}
	      
	      
	      // Output range of p2-passing SNPs in UCSC cut-and-paste friendly format
	      
	      
	      CLMP << "\n          RANGE: "
		   << "chr" << chromosomeName( P->locus[l]->chr ) << ":" 
		   << minBP << ".." << maxBP << "\n";
	      
	      CLMP << "           SPAN: "
		   << ( maxBP - minBP ) / 1000 << "kb\n";
	      
	    
	    }

	  ////////////////////////////////////
	  // Print ranges encountered here
	  
	  if ( par::clumpld_range_annotate )
	    {
	      
	      if ( grouped_snps.size() > 0 ) 
		{

		  set<string>::iterator ri = range_notes.begin();
		  bool first = true;
		  
		  CLMP << "     GENES w/SNPs: ";
		  while ( ri != range_notes.end() )
		    {
		      if ( first ) 
			{
			  CLMP << *ri;
			  first = false;
			}
		      else
			{
			  CLMP << "," << *ri;
			}
		      ++ri;
		    }
		  CLMP << "\n";
		}
	      
	      // Now list all genes in region
	      Range r1;
	      r1.start = minBP;
	      r1.stop = maxBP;
	      r1.chr = P->locus[l]->chr;
	      
	      if (  grouped_snps.size() == 0 )
		CLMP << "\n";

	      CLMP << "            GENES: ";
	      
	      
	      set<Range*> intRanges = rangeIntersect(r1,ranges);
	      set<Range*>::iterator ri2 = intRanges.begin();
	      int cnt = 0;
	      	     	      
	      while ( ri2 != intRanges.end() )
		{
		  if ( cnt == 0 )
		    {
		      CLMP << (*ri2)->name;
		    }
		  else if ( cnt % 8 == 0 )
		    {
		      CLMP << "\n                   " 
			   << (*ri2)->name;
		    }
		  else
		    {
		      CLMP << "," << (*ri2)->name;
		    }
		  ++ri2;
		  ++cnt;
		}
	      
	      CLMP << "\n";
	    } 
	  
	  
	}
      else
	{
	  
	  // Just list of SNP names
	  if( grouped_snps.size() == 0 )
	    CLMP << "NONE";
	  
	  map<int2,double>::iterator gi = grouped_snps.begin();
          int j = 0;
          while( gi != grouped_snps.end() )
            {
              CLMP << P->locus[ gi->first.p1 ]->name << "(" 
		   << gi->first.p2 << ")";
              if( j < grouped_snps.size()-1)
                CLMP << ",";
              j++;
              gi++;
            }

	  
	  // Non-verbose mode gene information 
	  // goes to a separate file (CLMP2)

	  if ( par::clumpld_range_annotate )
	    {

	      int minBP = P->locus[l]->bp;
	      int maxBP = P->locus[l]->bp;

	      map<int2,double>::iterator gi = grouped_snps.begin();
	      while( gi != grouped_snps.end() )
		{
		  int l0 = gi->first.p1;
		  if ( P->locus[l0]->bp < minBP )
		    minBP = P->locus[l0]->bp;
		  if ( P->locus[l0]->bp > maxBP )
		    maxBP = P->locus[l0]->bp;
		  ++gi;
		}
	      
	      Range r1;
	      r1.start = minBP;
	      r1.stop = maxBP;
	      r1.chr = P->locus[l]->chr;
	      string glist = returnFullRangeList(r1,ranges,false);
	      if ( glist == "" )
		glist = "(NONE)";
	      CLMP2 << setw(4) << P->locus[l]->chr << " " 
		    << setw(par::pp_maxsnp) << P->locus[l]->name << " " 
		    << setw(10) << assoc_results[cp].p << " "
		    << setw(6) << grouped_snps.size()+1 << " "
		    << setw(28) 
		    << ("chr"+chromosomeName(P->locus[l]->chr)+
			":"+int2str(minBP)+".."+int2str(maxBP)) << " "		
		    << setw(10) << ( maxBP - minBP ) / 1000.0 << " " 
		    << "["<<returnFullRangeList(r1,ranges,false) << "]\n";
	    }

	}
      
      
      if ( par::clumpld_verbose )
	CLMP << "\n"
	     << "------------------------------------------------------------------\n";

      CLMP << "\n";
      
      if ( ! par::silent)
	cout << "Formed " << ++clumpCount
	     << " clumps of the top " 
	     << i << " SNPs, so far            \r";
      
      grouped_snps.clear();


      // Next SNP
    }
  
  
  if ( ! par::silent)
    cout << "\n";
  
  CLMP << "\n" << vmessage << "\n";

  CLMP.close();

  if ( par::clumpld_range_annotate && ! par::clumpld_verbose )
    CLMP2.close();
  
  if ( par::clumpld_best )
    BEST.close();


}



string clump_LD::allelePairs(int l1, int l2)
{
  
  // The first SNP is always the index SNP
  // i.e. as we call rsq(l,l1) and rsq(l,l2) above

  // Identify common/common haplotype

  int ch = 0;
  for (int h=0; h<hp->nh; h++)
    if ( hp->haplotypeName(h) == P->locus[l1]->allele2 + P->locus[l2]->allele2 )
      ch = h;
  
  // Is D positive or negative?
  string s;

  if ( hp->f[ch] > (1 - P->locus[l1]->freq)*(1 - P->locus[l2]->freq) )
    s = P->locus[l1]->allele1 + P->locus[l2]->allele1 + "/" + P->locus[l1]->allele2 + P->locus[l2]->allele2;
  else
    s = P->locus[l1]->allele1 + P->locus[l2]->allele2 + "/" + P->locus[l1]->allele2 + P->locus[l2]->allele1;

  return s;
}



string returnFullRangeList(Range & r1, map<string, set<Range> > & ranges, bool verbose)
{
  string rlist = "";
  
  set<Range*> intRanges = rangeIntersect(r1,ranges);
  set<Range*>::iterator ri2 = intRanges.begin();
  int cnt = 0;
  
  while ( ri2 != intRanges.end() )
    {
      if ( cnt == 0 )
	{
	  rlist += (*ri2)->name;
	}
      else if ( verbose && cnt % 8 == 0 )
	{
	  rlist += "\n                   " + (*ri2)->name;
	}
      else
	{
	  rlist += "," + (*ri2)->name;
	}
      ++ri2;
      ++cnt;
    }
  return rlist;
}
