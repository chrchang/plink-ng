

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
#include <vector>
#include <map>
#include <iterator>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"
#include "cnv.h"
#include "crandom.h"
#include "model.h"

extern Plink * PP;

const double EPS_OVERLAP = 1e-6;


class sortedSegments 
{
public:
  Segment s;         
  bool operator< (const sortedSegments & b) const
  {
    if ( s.p1 < b.s.p1 ) return true;
    if ( s.p1 > b.s.p1 ) return false;
    if ( s.start < b.s.start ) return true;
    if ( s.start > b.s.start ) return false;
    if ( s.finish < b.s.finish ) return true;
    return false;
  }
};



double probOverlap(Segment * s, set<Range>& isection)
{

  // For this segement, calculate sum of all weights across all genes
  
  double dlen = PP->locus[s->finish]->bp - PP->locus[s->start]->bp + 1;

  set<Range>::iterator i = isection.begin();
  double wCount = 0;

  while ( i != isection.end() )
    {
      double glen = i->stop - i->start + 1;
      wCount += 1.0 / ( dlen + glen );		      
      ++i;
    }
  

  return wCount / 1000.00 ;

}

void CNVIndivReport::calculateResults()
{
  // Note -- T1, T2, T5 and T6 are based on *all* individuals
  // Tests T3, T4, T7 and T8 are based on all individuals with 1+ event
  
  if ( count == 0 )
    count = 1;
  
  t1 /= n;
  t2 /= n;
  t3 /= count;
  t4 /= count;

  t5 /= n;
  t6 /= n;
  t7 /= count;
  
  t8 /= n;


  //  t8 = t8 > 0 ? t5/t8 : 0 ;
  
  // Expected counts 

//   t9 /= n; // subset
//   t10 /= n; // baseline geneset

  //  t11 = t10>0 ? t9 / t10 : 0 ; 

  // Make main metric the ratio
  //  t8 = t11>0 ? t8 / t11 : 0 ;

}


void Plink::setUpForCNVList()
{

  if ( ! par::cnv_makemap ) 
    {

      ///////////////////////////////////////////////
      // .map file 

      checkFileExists(par::mapfile);

      printLOG("Reading marker information from [ " 
	       + par::mapfile + " ]\n");

      vector<bool> include;
      vector<int> include_pos(0);
      int nl_actual=0;
      
      readMapFile(par::mapfile,include,include_pos,nl_actual);


      ///////////////////////////////////////////////
      // .fam file 
      
      printLOG("Reading individual information from [ " 
	       + par::famfile + " ]\n");
      checkFileExists(par::famfile);
      readFamFile(par::famfile);  


      // Set some basics that we would have skipped by 
      // doing the above manually

      nl_all = locus.size();
      n = sample.size();
      prettyPrintLengths();

    }


}

void Plink::readCNVList()
{


  // These should be done already, but check

  nl_all = locus.size();
  n = sample.size();
  prettyPrintLengths();


  ///////////////////////////////////////////////
  // Intersect with a range file?
  
  set<Range> isection;
  set<Range> isection_baseline;
  map<Range,string> idescription;
  set<Range> iintersected;

  if ( par::cnv_intersect )
    {

      set<string> isubset;
      set<string> ifound;

      if ( par::cnv_intersect_subset )
	{
	  checkFileExists( par::cnv_intersect_subset_file );
	  printLOG("Reading intersection subset list from [ " 
		   + par::cnv_intersect_subset_file + " ]\n");
	  ifstream IN(par::cnv_intersect_subset_file.c_str(), ios::in);
	  
	  while ( ! IN.eof() )
	    {
	      string gname;
	      IN >> gname;
	      if ( gname=="" )
		continue;
	      isubset.insert(gname);
	    }

	  printLOG("Looking for subset of " 
		   + int2str( isubset.size() ) 
		   + " ranges\n");
	  IN.close();
	}
      

      checkFileExists( par::cnv_intersect_file );

      printLOG("Reading CNV intersection list from [ " 
	       + par::cnv_intersect_file + " ]\n");
      ifstream IN(par::cnv_intersect_file.c_str(), ios::in);
      
      while ( ! IN.eof() )
	{
	  // First three fields are CHR M1, M2
	  // Ignore rest of the line
	  
	  char cline[par::MAX_LINE_LENGTH];
	  IN.getline(cline,par::MAX_LINE_LENGTH,'\n');
	  string sline = cline;
	  if (sline=="") continue;
	  	  
	  string buf; 
	  stringstream ss(sline); 
	  vector<string> tokens; 
	  while (ss >> buf)
	    tokens.push_back(buf);
	  	  
	  if ( tokens.size() < 3 ) 
	    error("Problem with line:\n" + sline );

	  // If we are only extracting a subset of this list, then 
	  // we require a properly defined 4th column
	  
	  if ( par::cnv_intersect_subset )
	    {
	      if ( tokens.size() < 4 ) 
		error("Problem with line, no region name:\n" + sline );
	      
	      // Skip if this gene isn't on the list
	      if ( isubset.find( tokens[3] ) == isubset.end() )
		continue;
	      
	      // Add to list of found subset genes
	      ifound.insert( tokens[3] );
	    }

	  
	  // Otherwise load this in
	  
	  string chr = tokens[0];
	  string m1 = tokens[1];
	  string m2 = tokens[2];
	  
	  int p1,p2;
	  Range r;
	  if (chr=="") 
	    continue;
	  if ( ! from_string<int>( r.start,m1, std::dec ) )
	    error("Problem with position : " + m1 );
 	  if ( ! from_string<int>( r.stop,m2, std::dec ) )
	    error("Problem with position : " + m2 );

	  // Add any border
	  r.start -= par::cnv_region_border;
	  r.stop += par::cnv_region_border;

	  // Check for consistency
	  if ( r.start > r.stop ) 
	    error("Badly defined region:\n"+sline);

	  r.chr = getChromosomeCode(chr);
	  
	  isection.insert(r);
	  idescription.insert(make_pair(r,sline));

	  // Add region name to this one
	  r.name = tokens.size() >= 4 ? 
	    tokens[3] : "REGION-"+int2str(geneList.size()+1);
	  geneList.insert(r);   // A second, global copy

	}
      IN.close();


      

      if ( par::cnv_exclude ) 
	printLOG("Read " + int2str( isection.size() ) 
		 + " ranges to exclude from CNV list\n");
      else
	printLOG("Read " + int2str( isection.size() ) 
		 + " ranges to intersect with CNV list\n");
      

      if ( ifound.size() < isubset.size() )
	{
	  printLOG("Could not find " + int2str( isubset.size() - ifound.size() ) + " ranges\n");
	  printLOG("Writing this list to [ " + par::output_file_name + ".notfound ] \n");
	  ofstream O2;
	  O2.open( (par::output_file_name+".notfound").c_str() , ios::out);
	  set<string>::iterator i1 = isubset.begin();
	  while ( i1 != isubset.end() )
	    {
	      if ( ifound.find( *i1 ) == ifound.end() )
		O2 << *i1 << "\n";
	      ++i1;
	    }
	  O2.close();
	}

      
      ////////////////////////////////////////////////
      // Read in a second list of baseline ranges?

      if ( par::cnv_count_baseline )
	{
	  
	  checkFileExists( par::cnv_count_baseline_file );

	  printLOG("Reading CNV baseline count list from [ " 
		   + par::cnv_count_baseline_file + " ]\n");
	  ifstream IN1(par::cnv_count_baseline_file.c_str(), ios::in);
      
	  while ( ! IN1.eof() )
	    {
	      // First three fields are CHR M1, M2
	      // Ignore rest of the line
	      
	      vector<string> tokens = tokenizeLine( IN1 );
	      if ( tokens.size() == 0 )
		continue;
	      if ( tokens.size() < 3 ) 
		error("Problem with line:\n" + displayLine(tokens) );
	      
	      string chr = tokens[0];
	      string m1 = tokens[1];
	      string m2 = tokens[2];
	  
	      int p1,p2;
	      Range r;
	      if (chr=="") 
		continue;
	      if ( ! from_string<int>( r.start,m1, std::dec ) )
		error("Problem with position : " + m1 );
	      if ( ! from_string<int>( r.stop,m2, std::dec ) )
		error("Problem with position : " + m2 );

	      // Add any border
	      r.start -= par::cnv_region_border;
	      r.stop += par::cnv_region_border;
	      
	      // Check for consistency
	      if ( r.start > r.stop ) 
		error("Badly defined region:\n"+ displayLine( tokens ));

	      r.chr = getChromosomeCode(chr);
	  
	      isection_baseline.insert(r);
	      
	    }

	  IN1.close();

	  printLOG("Read " + int2str( isection_baseline.size() ) 
		   + " baseline count ranges from [ " 
		   + par::cnv_count_baseline_file + " ]\n");

	  
	}
      

    }




  ///////////////////////////////////////////////
  // .cnv file 

  printLOG("\nReading segment list (CNVs) from [ " 
	   + par::cnv_listname + " ]\n");
  
  checkFileExists( par::cnv_listname );

  ifstream IN;
  IN.open( par::cnv_listname.c_str() , ios::in );

  ofstream MOUT;
  if ( par::cnv_writelist )
    {      

      if ( par::output_file_name + ".cnv" == par::cnv_listname )
	error("CNV input and output file names cannot be the same, " 
	      + par::cnv_listname );
    }

  map<string,Individual*> uid;
  map<int2,int> mlocus;

  if ( ! par::cnv_makemap )
    {

      makePersonMap( *this, uid );

      for (int l=0; l<nl_all; l++)
	{
	  int2 p;
	  p.p1 = locus[l]->chr;
	  p.p2 = locus[l]->bp;
	  mlocus.insert(make_pair(p,l));
	}
    }


  map<string,Individual*>::iterator ii;
  map<int2,int>::iterator il1;
  map<int2,int>::iterator il2;
  
  int nseg=0;
  int n_mapped_to_person=0;
  int n_passed_filters=0;
  int n_intersects=0;
  int nall=0;

  set<int2> positions;


  while ( ! IN.eof() )
    {

      string fid, iid, chr, bp1, bp2;
      string type, scorestr, sitesstr;

      // FID IID CHR BP1 BP2 TYPE SCORE SITES
      
      IN >> fid >> iid 
	 >> chr >> bp1 >> bp2 
	 >> type >> scorestr >> sitesstr;
      
      
      if ( fid == "FID" || 
	   fid == "" )
	continue;

      nall++;

      // Lookup person
      if ( ! par::cnv_makemap )
	{
	  ii = uid.find( fid + "_" + iid );      
	  if ( ii == uid.end() )
	    continue;
	}
      
      ++n_mapped_to_person;

      // Correct type? 
      int t;
      if ( ! from_string<int>( t,type, std::dec ) )
	error("Problem with type specifier: " + type );
	    
      if ( par::cnv_del_only && t > 1 ) 
	continue;

      if ( par::cnv_dup_only && t < 3 ) 
	continue;

      
      int p1,p2;
      if ( ! from_string<int>( p1, bp1, std::dec ) )
	error("Problem with first position: " + bp1 );
      if ( ! from_string<int>( p2, bp2, std::dec ) )
	error("Problem with second position: " + bp2 );

      if ( p1 > p2 ) 
	error("Badly defined segment, " + bp1 + " > " + bp2);


      double score;
      int sites;
      
      if ( ! from_string<double>( score, scorestr, std::dec ) )
	error("Problem with score : " + scorestr );
	    
      if ( ! from_string<int>( sites, sitesstr, std::dec ) )
	error("Problem with sites : " + sitesstr );
      
      
      // Filters:

      double kb = (double)(p2 - p1) / 1000.0;

      if ( par::cnv_min_sites > 0 && par::cnv_min_sites > sites )
	continue;
      if ( par::cnv_max_sites > 0 && par::cnv_max_sites < sites )
	continue;

      if ( par::cnv_min_score > 0 && par::cnv_min_score > score )
	continue;
      if ( par::cnv_max_score > 0 && par::cnv_max_score < score )
	continue;

      if ( par::cnv_min_kb > 0 && par::cnv_min_kb > kb )
	continue;
      if ( par::cnv_max_kb > 0 && par::cnv_max_kb < kb )
	continue;

      ++n_passed_filters;


      ///////////////////////////////////////////////
      // Intersect with range as specified from file
      // (optionally counting overaps, instead)

      int segment_contains = 0;
      double weighted_segment_contains = 0;
      if ( par::cnv_intersect ) 
	{
	  if ( par::cnv_count ) 
	    {
	      segment_contains += 
		count_intersects(isection,getChromosomeCode(chr),p1,p2);
	      if ( par::cnv_weighted_gene_test )
		weighted_segment_contains += weighted_count_intersects(isection,getChromosomeCode(chr),p1,p2);
	    }
	  else
	    {
	      if ( par::cnv_exclude && 
		   intersects(isection,iintersected,getChromosomeCode(chr),p1,p2) )
		continue;
	      if ( (!par::cnv_exclude) && 
		   !intersects(isection,iintersected,getChromosomeCode(chr),p1,p2) )
		continue;
	    }	  
	  ++n_intersects;
	}


      ///////////////////////////////////////////////
      // Intersect with range as specified from a 
      // second file, to act as baseline count
      
      int segment_baseline = 0;
      double weighted_segment_baseline = 0;
      if ( par::cnv_count_baseline ) 
	{
	  segment_baseline += 
	    count_intersects(isection_baseline,getChromosomeCode(chr),p1,p2);
	  if ( par::cnv_weighted_gene_test )
	    weighted_segment_baseline += weighted_count_intersects(isection_baseline,getChromosomeCode(chr),p1,p2);
	}
      
      int2 p;
      p.p1 = getChromosomeCode( chr );
      
      if ( par::cnv_makemap )
	{

	  int2 p;
	  p.p1 = getChromosomeCode( chr );

	  // Start
	  p.p2 = p1;
	  positions.insert(p);

	  // End
	  p.p2 = p2;
	  positions.insert(p);
	  
	  // One position past end
	  p.p2++;
	  positions.insert(p);

	  
	}
      else
	{

	  // Can we map this to an exact marker?

	  p.p2 = p1;
	  il1 = mlocus.find( p );
	  
	  p.p2 = p2;
	  il2 = mlocus.find( p );

	  if ( il1 == mlocus.end() || 
	       il2 == mlocus.end() )
	    continue;

	  
	  // Seems okay, add segment to list

	  Segment s;
	  
	  s.start = il1->second;
	  s.finish = il2->second;
	  s.p1 = s.p2 = ii->second;
	  s.count = segment_contains;
	  s.baseline = par::cnv_count_baseline ? segment_baseline : 0 ;
	  s.weightedCount = weighted_segment_contains;
	  s.weightedBaseline = weighted_segment_baseline;
	  s.type = t;
	  s.score = score;
	  s.sites = sites;
	  segment.push_back(s);

	  

  	  if ( par::verbose ) 
  	    cout << "Adding " << s.count << " " << s.weightedCount << "\t"
  		 << " and  " << s.baseline << " " << s.weightedBaseline << "\tphe = "
		 << s.p1->phenotype << "\n";

	}
        
    }
  

  // Determine some measure of expected weighted count, given 
  // totality of CNVs and genes specified per individual
  
  expectedOverlap.resize(n,0);
  expectedOverlapBaseline.resize(n,0);

  if ( par::cnv_weighted_gene_test )
    {
      map<Individual*,int> mmap;
      for (int i=0; i<n; i++)
	mmap.insert(make_pair( sample[i],i ));

      vector<Segment>::iterator s0 = segment.begin();
      while ( s0 != segment.end() )
	{	  
	  int j = mmap.find(s0->p1)->second;
	  expectedOverlap[ j ] += probOverlap(&(*s0), isection );
	  	  
	  if ( par::cnv_count_baseline ) 
	    expectedOverlapBaseline[ j ] += probOverlap(&(*s0), isection_baseline );
	  
	  ++s0;
	}
      
    }

//   cout << "expected:\n";
//   display(expectedOverlap);
//   cout << "v2\n";
//   display(expectedOverlap);
  
  

  /////////////////////////////////////////////////////////////
  // Write-back any intersected segments (or, non-intersected 
  // segments, if in exclude mode)

  if ( par::cnv_intersect_writeback )
    {
      
      if ( par::cnv_exclude ) 
	printLOG("Writing back list to non-intersected regions to [ " 
		 + par::output_file_name + ".reg ]\n");
      else
	printLOG("Writing back list to intersected regions to [ " 
		 + par::output_file_name + ".reg ]\n");
      
      
      ofstream ROUT;
      ROUT.open( ( par::output_file_name+".reg").c_str(), ios::out );
      ROUT.precision(4);
      
      // Either a simple list of the intersected regions
      // or a verbose report that has both the region
      // and then a list of the CNVs in that region
      
      set<Range>::iterator ir = isection.begin();     
      while ( ir != isection.end() )
	{
	  set<Range>::iterator i = iintersected.find(*ir);  
	  bool writeback = par::cnv_exclude ? 
	    i == iintersected.end() :
	    i != iintersected.end() ; 	  
	  if ( writeback )
	    {	      
	      map<Range,string>::iterator is = idescription.find(*ir);	  

	      if ( par::cnv_intersect_writeback_verbose ) 		
		{
		  ROUT << "RANGE (+/- "  
		       << par::cnv_region_border/1000 
		       << "kb )  [ " 
		       << idescription[*ir] 
		       << " ]\n";

		  ROUT << setw(par::pp_maxfid) << "FID" << " "
		       << setw(par::pp_maxiid) << "IID" << " ";
		  ROUT << setw(8) << "PHE" << " ";
		  ROUT << setw(4) << "CHR" << " "
		       << setw(12) << "BP1" << " "   
		       << setw(12) << "BP2" << " "
		       << setw(6) << "TYPE" << " "
		       << setw(8) << "KB" << " "
		       << setw(8) << "OLAP" << " "
		       << setw(8) << "OLAP_U" << " "
		       << setw(8) << "OLAP_R" << "\n";

		    //	<< setw(12) << "SCORE" << " "
		    //  << setw(8) << "SITES" << "\n";
		  
		}
	      else
		ROUT << idescription[*ir] << "\n";
	      
	      // Get a list of CNVs that span this range
	      // and display in verbose mode
	      if ( par::cnv_intersect_writeback_verbose ) 		
		{
		  
		  Range trange = *ir;
		  set< Segment > segset = allSegmentsIntersecting( trange );
		  
		  set<Segment>::iterator s = segset.begin();
		  		  
		  while ( s != segset.end() )
		    {

		      Individual * person = s->p1;
		      Locus * loc1 = locus[s->start];
		      Locus * loc2 = locus[s->finish];
		      
		      ROUT << setw(par::pp_maxfid) << person->fid << " "
			   << setw(par::pp_maxiid) << person->iid << " ";
		      ROUT << setw(8) << person->phenotype << " ";
		      ROUT  << setw(4) << loc1->chr << " "
			    << setw(12) << loc1->bp << " "   
			    << setw(12) << loc2->bp << " ";
		      if ( s->type < 2 ) 
			ROUT << setw(6) << "DEL" << " ";
		      else
			ROUT << setw(6) << "DUP" << " ";


		      // Instead of SCORE and SITES, output kb length, and overlap stats
// 		      ROUT  << setw(12) << s->score << " "
// 			    << setw(8) << s->sites << "\n";
		      
		      // Start and stop of CNV
		      int p1 = loc1->bp;
		      int p2 = loc2->bp;
		      
		      // Start and stop of gene

		      double consensusStart = p1 > ir->start ? p1 : ir->start;
		      double consensusStop = p2 < ir->stop ? p2 : ir->stop;		      
		      double numerator = consensusStop - consensusStart + 1.0;
	
		      // 1,2,3 = default, union, region overlap

		      double denom1 = p2-p1+1.0;	      
		      
		      double unionStart = p1 < ir->start ? p1 : ir->start;
		      double unionStop = p2 > ir->stop ? p2 : ir->stop;
		      double denom2 = unionStop - unionStart + 1.0;

		      double denom3 = ir->stop-ir->start+1.0;
		            	      
		      double overlap1 = numerator / denom1 ;
		      double overlap2 = numerator / denom2 ;
		      double overlap3 = numerator / denom3 ;

		      ROUT << setw(8) << (loc2->bp - loc1->bp)/1000.0 << " "
			   << setw(8) << overlap1 << " "
			   << setw(8) << overlap2 << " "
			   << setw(8) << overlap3 << "\n";
		      
		      ++s;
		    }
		  ROUT << "\n";
		}
	    }
	  ++ir;	  
	}


      
      ROUT.close();
    }
  


  /////////////////////////////////////////////////////////////
  // Drop segments that are above or below a certain threshold

  int removed_freq_filter = 0;
  
  // Old approach was to define genomic regions 
  // and then filter based on these. 
  
  // New approach (below) is to focus on each event more directly, 
  // and ask how many others intersect with it.  Populate "freq"
  // for each event;  
  // Overlap definition wants to be 
  //  1) NOT "disrupt"
  //  2) Inforce union overlap so that frequency groups are 
  //     symmetric 
 
  // Clean up this code later (decide if we want to keep these old 
  // definitions

  if ( par::cnv_freq_method2 )
    {
      
      // Turn off potential write-back function for 
      // ranges; and manually set the intersect function
      // to work as needed here
      
      par::cnv_intersect_writeback = false;
      par::cnv_disrupt = false;
      par::cnv_union_overlap = true;
      par::cnv_region_overlap = false;
      
      for (int s1=0; s1<segment.size(); s1++)
	{
	  Segment * seg1 = &(segment[s1]);;
	  int chr = locus[ seg1->start ]->chr;
	  int bp1 = locus[ seg1->start ]->bp;
	  int bp2 = locus[ seg1->finish ]->bp;
	  
	  // Count overlap with self, once
	  ++seg1->freq;
	  
	  for (int s2=s1+1; s2<segment.size(); s2++)
	    {
	      Segment * seg2 = &(segment[s2]);
	      
	      // Skip if different chromosome
	      if ( locus[ seg2->start ]->chr != chr ) 
		continue;
	      
	      // Do these overlap?
	      
	      int t1 = locus[seg2->start ]->bp;
	      int t2 = locus[seg2->finish ]->bp;
	      
	      // This seg ends before the other starts?
	      if ( t2 < bp1 )
		continue;
	      
	      // Or does this seg start after the other ends?
	      if ( t1 > bp2 ) 
		continue;
	      
	      // Calculate overlap (union)
	      
	      double consensusStart = t1 > bp1 ? t1 : bp1;
	      double consensusStop = t2 < bp2 ? t2 : bp2;	      
	      double numerator = consensusStop - consensusStart + 1.0;
	      double unionStart = t1 < bp1 ? t1 : bp1;
	      double unionStop = t2 > bp2 ? t2 : bp2;
	      double denom = unionStop - unionStart + 1.0;
	      double overlap = numerator / denom ;
	      overlap += EPS_OVERLAP;
	      
	      if ( overlap > par::cnv_freq_method2_threshold )
		{
		  seg1->freq++;
		  seg2->freq++;
		}
	    }
	}
      
      
      
      if ( par::cnv_freq_include ) 
	{
	  printLOG("Filtering segments based on frequencies\n");	  
	  
	  vector<Segment>::iterator s = segment.begin();
	  
	  while ( s != segment.end() )
	    {
	      //cout << "seg count = " << s->freq << "\n";
	      
	      if ( par::cnv_freq_include_exact 
		   && s->freq != par::cnv_freq_include_cnt )
		{
		  s = segment.erase(s);
		  ++removed_freq_filter;
		}
	      else if ( par::cnv_freq_include_below 
			&& s->freq > par::cnv_freq_include_cnt )
		{
		  s = segment.erase(s);
		  ++removed_freq_filter;
		}
	      else if ( s->freq < par::cnv_freq_include_cnt )
		{
		  s = segment.erase(s);
		  ++removed_freq_filter;
		}
	      else
		s++;
	    }
	  
	  printLOG("Will remove " + int2str( removed_freq_filter ) 
		   + " CNVs based on frequency (after other filters)\n");
	  
	}
    }


  ////////////////////////////
  // Old frequency filter code
  
  if ( ! par::cnv_freq_method2 ) 
    {
      if ( par::cnv_freq_include || par::cnv_unique )
	{
	
	if ( par::cnv_unique )
	  printLOG("Filtering segments unique to cases/controls\n");
	
	if ( par::cnv_freq_include )
	  printLOG("Filtering segments based on frequencies\n");	  
	  
	  // Clear any existing segments (we are done with gene lists, etc)
	  // by now
	  
	  isection.clear();
	  
	  // 1) Find common regions
	  // 2) intersect or exclude as is fit
	  // 3) remove from list
	  
	  vector<int> caseCount = segmentCountCaseControls(this,true);
	  vector<int> controlCount = segmentCountCaseControls(this,false);
	  
	  // Determine regions to exclude
	  
	  bool inRegion = false;

	  Range r;
	  
	  for (int l=0; l<nl_all; l++)
	    {
	      
	      //////////////////////////////
	      // End of chromosome/genome?
	      
	      bool endOfChromosome = false;
	      
	      if ( l == nl_all-1 || locus[l+1]->chr != locus[l]->chr )
		{
		  endOfChromosome = true;
		}
	  
	      int count = caseCount[l] + controlCount[l];
	      bool uniq = ! ( caseCount[l] == 0 || controlCount[l] == 0 );
	      
// 	      cout << "dets " << locus[l]->chr << "\t" << locus[l]->bp << "\t" << uniq << "\t" << caseCount[l] << " : " 
// 		   << controlCount[l] << "\t";
	      
	      bool iregion;
	      
	      if ( par::cnv_freq_include_exact ) 
		{
		  // Define the regions we want to exclude

		  iregion = count != par::cnv_freq_include_cnt;

		}
	      else 
		{

		  // Define the regions we want to exclude

		  // Use inclusive thresholds for inclusion
		  // --cnv-freq-exclude-above X

		  if ( par::cnv_freq_include_below ) // inclusive thresholds
		    iregion = count > par::cnv_freq_include_cnt; 
		  else // --cnv-freq-exclude-below X
		    iregion = count < par::cnv_freq_include_cnt;  
		}
	      
	      //	      cout << iregion << "\t" << count << "\n";
	      
	      /////////////////////////////////////////////////////////
	      // Filter based on uniquenes to either cases or controls
	      
	      if ( inRegion )
		{
		  
		  // End of an iregion or unique region?
		  if ( ( ( par::cnv_unique && ! uniq ) || ( ! par::cnv_unique ) ) &&
		       ( ( par::cnv_freq_include && ! iregion ) || ( ! par::cnv_freq_include ) ) )
		    {
		      // Range goes up to the position just before this region
		      r.stop = locus[l]->bp-1;
		      if ( r.stop < 0 ) r.stop = 0;
		      isection.insert(r);
		      inRegion = false;
		    }
		  else if ( endOfChromosome ) // of we just have to stop anyway?
		    {
		      inRegion = false;
		      r.stop = locus[l]->bp;
		      isection.insert(r);
		      continue;
		    }
		}
	      // ...or, the start of a new region?
	      else 
		{
		  if ( ( par::cnv_unique && uniq ) || 
		       ( par::cnv_freq_include && iregion ) )
		    {
		      r.start = locus[l]->bp;
		      r.chr = locus[l]->chr;
		      inRegion = true;
		    } 
		  // But could also be the end of the chromosome
		  if ( inRegion && endOfChromosome )
		    {
		      inRegion = false;
		      r.stop = locus[l]->bp;
		      isection.insert(r);
		      continue;
		    }
		}
	    } // Next SNP
	  
	  
	  /////////////////////////////////////////////////////////////
	  // We now have a list of "iregion sections"  -- intersect 
	  // based on these
	  
	  set<Range>::iterator i2 = isection.begin();


//     	  cout << "List of bad regions, N = " << isection.size() << "\n";
//      	  while (i2 != isection.end() )
//      	    {
//      	      cout << "regs = " << i2->chr << "\t" << i2->start << " -> " << i2->stop << "\n";	
//      	      ++i2;
//      	    }
//    	  cout << "\n";


	  // Turn off potential write-back function for 
	  // ranges
	  
	  par::cnv_intersect_writeback = false;
	  
	  vector<Segment>::iterator s = segment.begin();

	  while ( s != segment.end() )
	    {
	      
	      bool doesIntersect = intersects(isection,
					      iintersected,
					      locus[s->start]->chr,
					      locus[s->start]->bp,
					      locus[s->finish]->bp);
	      
 	      if ( par::cnv_freq_include_exact_exclude )
 		{
		  doesIntersect = ! doesIntersect;
		}

	      bool willErase = doesIntersect;

// 	      bool willErase = par::cnv_freq_include_exact ? 
// 		! doesIntersect : doesIntersect ;
	      
	      // Remove this segment or keep?
	      if ( willErase )
		{
		  s = segment.erase(s);
		  ++removed_freq_filter;		  
		}
	      else
		++s;
	    }

	  
	  printLOG("Will remove " + int2str( removed_freq_filter ) 
		   + " CNVs based on frequency (after other filters)\n");
	  
	}
    }
  
  
  
  

  // Drop individuals for whom we do not see any segments of this
  // variety?

  if ( par::cnv_drop_no_segment )
    {
      vector<bool> indel(n,false);      
      indivSegmentSummaryCalc(segmentCount, segmentLength, true, true);
      
      for ( int i = 0; i < n; i++)
	{
	  indivPair t;
	  t.p1 = t.p2 = sample[i];
	  map<indivPair,int>::iterator ic = segmentCount.find(t);
	  map<indivPair,double>::iterator il = segmentLength.find(t);
	  	  
	  if ( ic != segmentCount.end() ) 
	    {
	      if ( ic->second == 0 )
		indel[i] = true;
	    }
	  else
	    indel[i] =true;
	}
      
      int n_removed = deleteIndividuals(indel);
      printLOG("Removed " + int2str(n_removed) + " individuals with <1 CNV\n");      
    }


  ///////////////////////////////////////////////////////////
  // Check that CNVs do not overlap at all (within a person)
    
  if ( par::cnv_check_overlap ) 
    {      
      set<sortedSegments> sorted;
      vector<Segment>::iterator s1 = segment.begin();
      while ( s1 != segment.end() )
	{
	  sortedSegments ss;
	  ss.s = *s1;
	  sorted.insert( ss );
	  ++s1;
	}
      

      // Look for overlap
      set<Segment*> olap;
      set<sortedSegments>::iterator s = sorted.begin();
      while ( s != sorted.end() )
	{
	  set<sortedSegments>::iterator snext = s;
	  ++snext;
	  if ( snext == sorted.end() )
	    break;
	  if ( s->s.p1 != snext->s.p1 )
	    {
	      ++s;
	      continue;
	    }
	  if ( snext->s.start <= s->s.finish ) 
	    {
	      olap.insert( (Segment*)&(s->s) );
	      olap.insert( (Segment*)&(snext->s) );
	    }

	  ++s;
	}

      if ( olap.size() > 0 )
	{
	  printLOG("Within-individual CNV overlap detected, involving " 
		   + int2str( olap.size() ) + " CNVs\n");
	  printLOG("Writing list to [ " + par::output_file_name + ".cnv.overlap ]\n");
	  ofstream O( ( par::output_file_name + ".cnv.overlap").c_str() , ios::out );

	  O << setw( par::pp_maxfid ) << "FID" << " " 
	    << setw( par::pp_maxfid ) << "IID" << " " 
	    << setw( 4 ) << "CHR" << " " 
	    << setw( 12 ) << "BP1" << " " 
	    << setw( 12 ) << "BP2" << "\n"; 

	  set<Segment*>::iterator i = olap.begin();
	  while ( i != olap.end() )
	    {
	      O << setw( par::pp_maxfid ) << (*i)->p1->fid << " " 
		<< setw( par::pp_maxfid ) << (*i)->p1->iid << " " 
		<< setw( 4 ) << locus[ (*i)->start ]->chr << " " 
		<< setw( 12 ) << locus[ (*i)->start ]->bp << " " 
		<< setw( 12 ) << locus[ (*i)->finish ]->bp << "\n"; 
	      ++i;
	    }	  
	  O.close();
	} 
      else
	{
	  printLOG("No overlapping samples found\n");
	}
    }

      
  ////////////////////////////////////////////////////////
  // Get full typeCount; for case/control data, split out 
  // by cases and controls
  
  map<int,int> typeCount;
  map<int,int> typeCaseCount;
  
  vector<Segment>::iterator s = segment.begin();
  while ( s != segment.end() )
    {
      
      //////////////////////////////////
      // Record in full 0/1, 3/4 space
      
      map<int,int> & myCount = par::bt && s->p1->aff ? 
	typeCaseCount : typeCount;
      
      map<int,int>::iterator it = myCount.find( s->type );
            
      if ( it == myCount.end() )
	{
	  myCount.insert(make_pair(s->type,1));
	}
      else
	{
	  ++(it->second);
	}

      // Next segment
      ++s;
    }



  //////////////////////////////
  // Make a new map file

  if ( ! par::cnv_makemap )
    {

      printLOG(int2str( n_mapped_to_person ) + " mapped to a person, of which " 
	       + int2str( n_passed_filters) + " passed filters\n");

      if ( par::cnv_intersect )
	{
	  if ( par::cnv_exclude )
	    printLOG( int2str( n_intersects ) 
		      + " kept after excluding specific regions\n");
	  else
	    printLOG( int2str( n_intersects ) 
		      + " intersected with one or more specified region\n");
	}
      
      int t = par::cnv_intersect ? n_intersects : n_passed_filters;
      t -= removed_freq_filter;
      if ( t - segment.size()  > 0  )
	printLOG( int2str( t - segment.size() ) 
		   + " did not map to any marker\n");

      printLOG(int2str( segment.size() ) 
	       + " of " + int2str(nall) 
	       + " mapped as valid segments\n");
      

      map<int,int>::iterator it1 = typeCaseCount.begin();
      map<int,int>::iterator it2 = typeCount.begin();
      set<int> obsCounts;

      while ( it1 != typeCaseCount.end() )
	{
	  obsCounts.insert( it1->first );
	  ++it1;
	}
      while ( it2 != typeCount.end() )
	{
	  obsCounts.insert( it2->first );
	  ++it2;
	}
      
      stringstream s2;
      
      if ( par::bt )
	s2 << setw(6) << "CopyN" << " " 
	   << setw(12) << "Case/Control" << "\n";
      else
	s2 << setw(6) << "CopyN" << " "
	   << setw(8) << "Count" << "\n";
      
      printLOG( s2.str() );
      s2.clear();
      
      set<int>::iterator it = obsCounts.begin();

      while ( it != obsCounts.end() )
	{
	  map<int,int>::iterator i1 = typeCaseCount.find( *it );
	  map<int,int>::iterator i2 = typeCount.find( *it );
	  
	  int n1 = 0, n2 = 0;

	  if ( i1 != typeCaseCount.end() )
	    n1 = i1->second;
	  
	  if ( i2 != typeCount.end() )
	    n2 = i2->second;

	  stringstream s;

	  if ( par::bt )
	    s << setw(6) << *it << " "
	      << setw(12) << (int2str(n1)+" / "+int2str(n2)) << "\n";
	  else
	    s << setw(6) << *it << " "
	      << setw(8) << n2 << "\n";
	  	  
	  printLOG( s.str() );

	  ++it;
	}
      printLOG("\n");
    }
  
  IN.close();
  


  ////////////////////////////////
  // Write a CNV list file? 

  if ( par::cnv_writelist )
    {

      printLOG("Writing new CNV list to [ " 
	       + par::output_file_name 
	       + ".cnv ]\n");
      
      MOUT.open( ( par::output_file_name + ".cnv").c_str() , ios::out );
      MOUT << setw(par::pp_maxfid) << "FID" << " "
	   << setw(par::pp_maxiid) << "IID" << " ";
      if (par::dump_covar_with_phenotype )
	MOUT << setw(8) << "PHE" << " ";
      MOUT << setw(4) << "CHR" << " "
	   << setw(12) << "BP1" << " "   
	   << setw(12) << "BP2" << " "
	   << setw(6) << "TYPE" << " "
	   << setw(12) << "SCORE" << " "
	   << setw(8) << "SITES" << " ";
      
      if ( par::cnv_write_freq )
	MOUT << setw(8) << "FREQ" << " ";
	  
      MOUT << "\n";
      
      vector<Segment>::iterator s = segment.begin();
      while ( s != segment.end() )
	{
	  
	  Individual * person = s->p1;
	  Locus * loc1 = locus[s->start];
	  Locus * loc2 = locus[s->finish];

	  MOUT << setw(par::pp_maxfid) << person->fid << " "
	       << setw(par::pp_maxiid) << person->iid << " ";
	  if (par::dump_covar_with_phenotype )
	    MOUT << setw(8) << person->phenotype << " ";
	  MOUT  << setw(4) << loc1->chr << " "
		<< setw(12) << loc1->bp << " "   
		<< setw(12) << loc2->bp << " "
		<< setw(6) << s->type << " "
		<< setw(12) << s->score << " "
		<< setw(8) << s->sites << " ";

	  if ( par::cnv_write_freq )
	    MOUT << setw(8) << s->freq << " ";   
	  
	  MOUT << "\n";
	 	  
	  ++s;
	}
      MOUT.close();

      printLOG("Writing new FAM file to [ "
               + par::output_file_name
               + ".fam ]\n");
      MOUT.open( ( par::output_file_name + ".fam").c_str() , ios::out );

      for (int i=0;i<n;i++)
	{
	  Individual * person = sample[i];
	  MOUT << person->fid << " "
	       << person->iid << " "
	       << person->pat << " "
	       << person->mat << " "
	       << person->sexcode << " ";
	  if (par::bt)
	    MOUT << (int)person->phenotype << "\n";
	  else
	    MOUT << person->phenotype << "\n";
	}      
      MOUT.close();
    }  


  /////////////////////////////////////
  // Collapse type (0/1 -> 2, 3/4 -> 2)

  s = segment.begin();
  while ( s != segment.end() )
    {
      s->type = s->type < 2 ? 1 : 2;
      ++s;
    }


  ////////////////////////////////
  // Write a new MAP file out?

  if ( par::cnv_makemap )
    {
      ofstream MOUT;
      printLOG("Writing new MAP file to [ " 
	       + par::output_file_name + ".cnv.map ]\n");
      MOUT.open( ( par::output_file_name + ".cnv.map").c_str() , ios::out );
      set<int2>::iterator ip = positions.begin();
      int nseg = 1;
      while ( ip != positions.end() )
	{
	  MOUT << ip->p1 << "\t"
	       << "p"+int2str(ip->p1)+"-"+int2str(ip->p2) << "\t"
	       << "0\t"
	       << ip->p2 << "\n"; 
	  ++ip;
	}
      MOUT.close();
      
      printLOG("Wrote " + int2str( positions.size() ) 
	       + " unique positions to file\n");
    }

 
}


void Plink::processCNVList()
{

  if (  par::cnv_writelist || par::cnv_makemap )
    return;

  // Per-individual summaries
  printLOG("Writing per-individual summary to [ " 
	   + par::output_file_name +
	   ".cnv.indiv ]\n");
  indivSegmentSummary();

  if ( par::cnv_enrichment_test )
    {
      glmCNVBurdenModel(*pperm,true);
    }


  // Display convenient segment view

  if (par::display_segment_long)
    displaySegmentsLong();


  // Display convenient BED/UCSC track

  if (par::display_cnv_track)
    displaySegmentsBED();


  // Find overlap in segments?

  if (par::segment_overlap)
    summariseHomoRuns();  
  


  /////////////////////////////////
  // Association/burden mapping
  

  ///////////////////////////////////////////////////////
  // Set up scaffold based on gene positions, if needed

  set<Range>::iterator g = geneList.begin();
  map<int,int> minpos;
  map<int,int> maxpos;
  while ( g != geneList.end() )
    {
      if ( minpos.find( g->chr ) == minpos.end() )
	{
	  minpos.insert(make_pair( g->chr, g->start ) );
	  maxpos.insert(make_pair( g->chr, g->stop ) );
	}
      else
	{
	  if ( minpos[g->chr] > g->start ) 
	    minpos[g->chr] = g->start;
	  if ( maxpos[g->chr] < g->stop )
	    maxpos[g->chr] = g->stop;
	}
      ++g;
    }

  map<int,int>::iterator i = minpos.begin();
  while ( i != minpos.end() )
    {
      CInfo cdet;
      cdet.bpstart = i->second;      
      cdet.bpstop = maxpos[i->first];
      scaffold.insert(make_pair(i->first,cdet));
      ++i;
    }


//   if ( par::cnv_glm )
//     runTestCNVwithGLM(*pperm);
 
  if ( par::qt ) 
    runTestCNVwithQT(*pperm);
  else
    {

      // Disease trait tests
      
      // Utilise existing tests for homozygous segments
      
      if ( par::seg_test_region )
      {
	printLOG("Writing positional summary to [ " 
		 + par::output_file_name 
		 + ".cnv.regional.summary ]\n");

	initialiseGeneCountAssociation(*pperm);
      }
      else
      {
	printLOG("Writing positional summary to [ " 
		 + par::output_file_name 
		 + ".cnv.summary ]\n");
	
	summaryIBSsegments(*pperm);
	
      }

    }
  
}


bool intersects(set<Range>& isection , 
		set<Range> & iintersects, 
		int chr, int p1, int p2)
{

  // Does this particular CNV intersect with atleast one range?
  // Potentially allowing for fractional overlaps)

  // Either consider all ranges and report back true when 
  // first intersection is seen
  
  // Or; consider all ranges, no matter what, keeping track
  // of which ranges have been intersected (i.e. if we have
  // overlapping ranges, for example)
  
  // Intersect: default -- if any of segment overlaps
  //      param 0.2     -- if at least20% of segment overlaps
  
  // Alternatively, "disrupt" mode simply asks whether the start 
  // or end of a CNV is within the region

  bool doesIntersect = false;
  
  set<Range>::iterator ir = isection.begin();
  
  while ( ir != isection.end() )
    {
      if ( ir->chr != chr )
	{
	  ++ir;
	  continue;
	}

      // Either use disrupt mode, or standard
      // intersect mode (which might include 
      // and overlap)

      if ( par::cnv_disrupt )
	{
	  
	  if (  ( p1 >= ir->start && p1 <= ir->stop ) ||	   
		( p2 >= ir->start && p2 <= ir->stop ) )
	    {
	      if ( par::cnv_intersect_writeback )
		{
		  doesIntersect = true;
		  iintersects.insert(*ir);
		}
	      else
		return true;
	    }
	}
      else // ... intersect or exclude mode
	{

	  if ( p1 <= ir->stop && p2 >= ir->start )
	    {

	      if ( par::cnv_overlap < 0 )
		{
		  
		  if ( par::cnv_intersect_writeback )
		    {
		      doesIntersect = true;
		      iintersects.insert(*ir);
		    }
		  else
		    return true;
		}
	      
	      /////////////////////////////
	      // Overlap-based comparison
		
	      // The CNV spans p1 to p2
	      // Region spans ir->start/stop

              // Denominator either of CNV itself, 
              // or is union

	      double consensusStart = p1 > ir->start ? p1 : ir->start;
	      double consensusStop = p2 < ir->stop ? p2 : ir->stop;
	      
	      double numerator = consensusStop - consensusStart + 1.0;

	      double denom;
	      if ( par::cnv_union_overlap )
		{
		  double unionStart = p1 < ir->start ? p1 : ir->start;
		  double unionStop = p2 > ir->stop ? p2 : ir->stop;
		  denom = unionStop - unionStart + 1.0;
		}
	      else if ( par::cnv_region_overlap )
		denom = ir->stop-ir->start+1.0;
	      else
		denom = p2-p1+1.0;

	      double overlap = numerator / denom ;
	      overlap += EPS_OVERLAP;

	      if ( overlap >= par::cnv_overlap )
		{
		  if ( par::cnv_intersect_writeback )
		    {
		      doesIntersect = true;
		      iintersects.insert(*ir);
		    }
		  else
		    return true;
		}
	    }
	}
      
      // Consider next range
      ++ir;
    }
  
  return doesIntersect;

}


int count_intersects(set<Range>& isection ,int chr,int p1,int p2)
{

  // How many ranges does this CNV intersect with?
  // Potentially allowing for fractional overlaps)

  // Intersect: default -- if any of segment overlaps
  //      param 0.2     -- if at least20% of segment overlaps

  int iCount = 0;

  set<Range>::iterator ir = isection.begin();
  while ( ir != isection.end() )
    {
      if ( ir->chr != chr )
	{
	  ++ir;
	  continue;
	}

      // Either use disrupt mode, or standard
      // intersect mode (which might include 
      // and overlap)

      if ( par::cnv_disrupt )
	{
	  
	  if (  ( p1 >= ir->start && p1 <= ir->stop ) ||	   
		( p2 >= ir->start && p2 <= ir->stop ) )
	    {
	      ++iCount;
	    }
	}
      else // ... intersect or exclude mode
	{	  
	  
	  if ( p1 < ir->stop && p2 > ir->start )
	    {
	      
	      if ( par::cnv_overlap < 0 ) 
		++iCount;
	      else
		{
		  // The CNV spans p1 to p2
		  double consensusStart = p1 > ir->start ? p1 : ir->start;
		  double consensusStop = p2 < ir->stop ? p2 : ir->stop;
		  double numerator = consensusStop - consensusStart + 1.0;
		  
		  double denom;
		  if ( par::cnv_union_overlap )
		    {
		      double unionStart = p1 < ir->start ? p1 : ir->start;
		      double unionStop = p2 > ir->stop ? p2 : ir->stop;
		      denom = unionStop - unionStart + 1.0;
		    }
		  else if ( par::cnv_region_overlap )
		    denom = ir->stop-ir->start+1.0;
		  else
		    denom = p2-p1+1.0;
		  
		  double overlap = numerator / denom ;
		  
		  overlap += EPS_OVERLAP;
		  
		  if ( overlap >= par::cnv_overlap )
		    ++iCount;
		}
	    }
	}
      
      ++ir;
    }
  return iCount;
}


double weighted_count_intersects(set<Range>& isection ,int chr,int p1,int p2)
{

  // For simple overlap or disruption statistics, calculate 
  // weighted score

  // Currently, no support for fractional overlaps

  double wCount = 0;
  int cnt = 0;

  set<Range>::iterator ir = isection.begin();
  while ( ir != isection.end() )
    {
      if ( ir->chr != chr )
	{
	  ++ir;
	  continue;
	}

      // Either use disrupt mode, or standard
      // intersect mode (which might include 
      // and overlap)

      if ( par::cnv_disrupt )
	{
	  
 	  if (  ( p1 >= ir->start && p1 <= ir->stop ) ||	   
 		( p2 >= ir->start && p2 <= ir->stop ) )
 	    {
	      // if ( D < G ) then W = D + G
	      // if ( D >= G ) then W = D + G - ( D - G + 1 ) 

	      double dlen = p2 - p1 + 1;
	      double glen = ir->stop - ir->start + 1;
	      double score = dlen < glen ? dlen + glen : dlen + glen - ( dlen-glen+1);
	      wCount += 1.0 / score;	      
	    }
	}
      else // ... intersect or exclude mode
	{	  
	  
	  if ( p1 < ir->stop && p2 > ir->start )
	    {
	      double dlen = p2 - p1 + 1;
	      double glen = ir->stop - ir->start + 1;
	      wCount += 1.0 / ( dlen + glen );		      


	      //Make 'weight' the intersection 
	      
// 	      double consensusStart = p1 > ir->start ? p1 : ir->start;
// 	      double consensusStop = p2 < ir->stop ? p2 : ir->stop;
// 	      wCount += consensusStop - consensusStart + 1.0;


// 	      cout.precision(8);
// 	      cout << "wCount = " << wCount << " from " << dlen << " " << glen << "\n";
// 	      if ( ++cnt == 2 ) cout << "TWO!\n";
	    }
	}
      
      ++ir;
    }

  // In 1 / KB units  
  wCount *= 1000;



  return wCount;
}



vector<int> segmentCountCaseControls(Plink * P, bool countCases)
{

  // Helper function to count segments (for CNVs, homozygous segments
  // only)

  vector<int> count(P->nl_all,0);
    
  vector<Segment>::iterator s = P->segment.begin();
  
  while ( s != P->segment.end() )
    {          
      if ( s->p1->pperson->aff == countCases ) 
	for (int l = s->start ; l <= s->finish; l++) 
	  count[l]++;
      s++;
    }
  return count;
}


set<Segment> allSegmentsIntersecting(Range & r)
{

  set<Segment> sset;
  set<Range> dummySet;    
  set<Range> rangeSet;
  rangeSet.insert(r);
    
  par::cnv_intersect_writeback = false;

  vector<Segment>::iterator s = PP->segment.begin();
  while ( s != PP->segment.end() )
    {
      
      bool testIntersection = intersects(rangeSet,
					 dummySet,
					 PP->locus[ s->start ]->chr,
					 PP->locus[ s->start ]->bp,
					 PP->locus[ s->finish ]->bp );
      
      // Record this CNV 
      if ( ( par::cnv_exclude && ! testIntersection ) ||
	   ( testIntersection && ! par::cnv_exclude ) )
	{
	  
	  sset.insert( *s );
	}
      
      ++s;
    }
  return sset;
}



void Plink::countCNVPerRegion(vector<int> & caseCount,
			      vector<int> & controlCount )
{

  // Return the number of case and control CNVs in each gene/region
  // These values are calculated first time, then stored in
  // gene2segment for subsequent lookup

  int nGenes = geneList.size();
  
  caseCount.clear();
  controlCount.clear();
  caseCount.resize(nGenes,0);
  controlCount.resize(nGenes,0);
  
  int gCount = 0;
  
  set<Range>::iterator g = geneList.begin();
  
  while ( g != geneList.end() )
    {
      
      Range tr = *g;

      // Either calculate, of if already done, just lookup the CNVs that
      // fall in this range (i.e. this function will be called many times
      // when using permutation)

      map<Range, set<Segment> >::iterator i = gene2segment.find( tr );
      if ( i == gene2segment.end() )
	{
	  set<Segment> s = allSegmentsIntersecting( tr );
	  gene2segment.insert(make_pair( tr, s ) );
	  i = gene2segment.find( tr );
	}
            
      set<Segment>::iterator si = i->second.begin();
      
      while ( si != i->second.end() )
	{
	  if ( si->p1->pperson->missing )
	    continue;
	  
	  if ( si->p1->pperson->aff )
	    ++caseCount[ gCount ];
	  else
	    ++controlCount[ gCount ];
	  
	  ++si;
	}

      ++gCount;
      ++g;
    }
  
  return;  

}


//////////////////////////////////////////////
// General segmental permutation test routine

void Plink::initialiseGeneCountAssociation( Perm & perm )
{

  printLOG("Performing region-based association mapping\n");
  int nt = geneList.size();
    
  vector<int> caseCount(nt,0);
  vector<int> controlCount(nt,0);


  // Count up number of CNVs in each gene/region
  
  countCNVPerRegion(caseCount,
		    controlCount);
  
  
  //////////////////////////////
  // And display

  string f = par::output_file_name + ".cnv.regional.summary";  

  ofstream SIBS;
  SIBS.open( f.c_str() , ios::out );
  
  SIBS << setw(4) << "CHR" << " " 
       << setw(16) << "REGION" << " " 
       << setw(12) << "BP1" << " "
       << setw(12) << "BP2" << " "
       << setw(8) << "AFF" << " " 
       << setw(8) << "UNAFF" << "\n"; 
  
  set<Range>::iterator ri = geneList.begin();
  int rCount = 0;
  while ( ri != geneList.end() )
    {      
      SIBS << setw(4) << ri->chr << " " 
	   << setw(16) << ri->name << " " 
	   << setw(12) << ri->start << " "
	   << setw(12) << ri->stop << " "
	   << setw(8) << caseCount[rCount] << " " 
	   << setw(8) << controlCount[rCount] << "\n";

      ++ri;
      ++rCount;
    }
  
  SIBS.close();
  
  
  /////////////////////
  // Permutation test?
  
  if (!par::permute) return;
  
  
  /////////////////////////////////////
  // Treat CNVs as homozygous segments
  
  par::homo_run = true;
 
  homozygousSegmentPermutationTest(perm,f,caseCount,controlCount);
  
  return;

}


void Plink::positionPermuteSegments()
{
  // Based on the global geneList set of genes, 
  // create a permuted version and rescore each segment count

  set<Range> permGeneList = geneList;
  
  // Permute on a within-chromosome basis
  //  cout << scaffold.size() << "\n";
  
  vector<int> offset( scaffold.size() , 0 );

  for (int c = 0 ; c < scaffold.size(); c++)
    {
      offset[c] = (int) ( CRandom::rand() * ( scaffold[c].bpstop - scaffold[c].bpstart ) );
    }    

  set<Range>::iterator g = permGeneList.begin();
  
  while ( g != permGeneList.end() )
    {
      
      Range * pg = (Range*)&(*g);

//       cout << "original =  " << pg->name << " : " << pg->start << " to " << pg->stop << " ";
//       cout << "offset = " << offset[pg->chr] << " ; ";
      
      pg->start += offset[pg->chr];
      pg->stop += offset[pg->chr];
      
      // Did we fall off the edge?

      if ( pg->stop >= scaffold[pg->chr].bpstop )
	{
	  //	  cout << " ** ";
	  int gsize = pg->stop - pg->start + 1; 
	  pg->start = scaffold[pg->chr].bpstart + ( pg->stop - scaffold[pg->chr].bpstop );
	  pg->stop  = pg->start + gsize;
	}           

      //      cout << "now " << pg->start << " to " << pg->stop << "\n";

      // Shuffle next gene
      ++g;
    }

  // i.e. update the "count" variable (simple # of genes intersected)

  if ( ! par::cnv_count ) 
    error("This function requires --cnv-count");

  // Update each segment
  
  vector<Segment>::iterator s = segment.begin();
  while ( s != segment.end() )
    {
      if ( ! par::cnv_weighted_gene_test )
	{
	  s->count = count_intersects(permGeneList , locus[s->start]->chr, locus[s->start]->bp, locus[s->finish]->bp );
	}
      else
	s->weightedCount = weighted_count_intersects(permGeneList , locus[s->start]->chr, locus[s->start]->bp, locus[s->finish]->bp );
      ++s;
    }

}




vector_t Plink::glmCNVBurdenModel(Perm & perm, bool print )
{
  vector_t results;

  string f = par::output_file_name + ".cnv.burden"; 
  
  if ( print )
    printLOG("Performing GLM-based CNV burden test, results in [ " + f + " ]\n");

  ///////////////////////////////////////////////////
  // Set up association model

  bool OLD_assoc_glm_without_main_snp = par::assoc_glm_without_main_snp;
  bool OLD_clist = par::clist;

  par::assoc_glm_without_main_snp = true;
  par::clist = true;

  // Model 0   no covariates
  //       1   CNV count
  //       2   CNV avg size
  //       3   CNV total kb
  //       4   CNV count + CNV avg size

  bool noCovar = false;

  if ( par::cnv_en_model == 0 )
    noCovar == true;
  

  /////////////////////////////////////////////////////
  // First calculate:
  //   What is mean CNV size for all observed CNVs?
  //   Do we see variation in CNV count in the dataset?

  double cnt = 0;
  double avg = 0;
  set<int> num;

  for ( int i = 0; i < n; i++)
    {
      indivPair t;
      t.p1 = t.p2 = sample[i];
      
      map<indivPair,int>::iterator ic = segmentCount.find(t);
      map<indivPair,double>::iterator il = segmentLength.find(t);
      
      // Insert correct # of terms here:
      if ( ic != segmentCount.end() )
	{
	  num.insert( ic->second );
	  cnt += ic->second;
	  avg += il->second;
	}
      else 
	num.insert( 0 );
      
    }
  avg /= cnt;

  if ( cnt == 0 ) 
    return results;

  // If equal burden, then do not include COUNT covariate

  bool equalBurden = num.size() == 1 ? true : false;

  // The number of additional terms we want to add go here (at the end of the clist)

  int totalTerms = 2;  // gene count + covariate
  if ( par::cnv_en_model == 4 && equalBurden ) par::cnv_en_model = 2;
  if ( par::cnv_en_model == 0 ) noCovar = true;
  if ( par::cnv_en_model == 1 && equalBurden ) noCovar = true;
  if ( par::cnv_en_model == 4 ) totalTerms = 3; 
  if ( noCovar ) totalTerms = 1;

  par::clist_number += totalTerms;
  for (int i=0; i<n; i++)
    sample[i]->clist.resize( par::clist_number );


  // If we want to use permutation, record this guy
  // And make this the first term

  int testTerm = par::clist_number - totalTerms;

  // Fill in label forms

  clistname.resize( par::clist_number );
  clistname[ testTerm ] = "GCNT";
  if ( ! noCovar )
    {
      if ( par::cnv_en_model == 1 )
	clistname[ testTerm + 1 ] = "NSEG";
      else if ( par::cnv_en_model == 2 )
	clistname[ testTerm + 1 ] = "AVGKB";
      else if ( par::cnv_en_model == 3 )
	clistname[ testTerm + 1 ] = "TOTKB";
      else if ( par::cnv_en_model == 4 )
	{
	  clistname[ testTerm + 1 ] = "NSEG";
	  clistname[ testTerm + 2 ] = "AVGKB";
	}
    }


  // Now these are part of Plink class
  // map<indivPair,int> segmentCount;
  // map<indivPair,double> segmentLength;
  
  // Generate basic summary file for all people
  indivSegmentSummaryCalc(segmentCount, segmentLength, true, true);
  

  if ( ! par::cnv_count )
    error("This function requires a --cnv-count set to be specified");


  ////////////////////////////////////////
  // Put relevant variables in clist slots
  
  for ( int i = 0; i < n; i++)
    {

      indivPair t;
      t.p1 = t.p2 = sample[i];

      map<indivPair,int>::iterator ic = segmentCount.find(t);
      map<indivPair,double>::iterator il = segmentLength.find(t);
      map<indivPair,double>::iterator ic2 = segmentCount2.find(t);
	  
      indivPair p = ic->first;
      
      // Insert correct # of terms here:
      if ( ic != segmentCount.end() )
	{
	  
	  double countCNV = ic->second;
	  double totalKB = il->second;
	  double avgKB = il->second / (double)ic->second;
	  double geneCount = ic2->second;

	  // Model = Y ~ gene-count + cnv-count + avg-cnv-size
	  
	  sample[i]->clist[ testTerm ] = geneCount;

	  if ( ! noCovar )
	    {
	      if ( par::cnv_en_model == 1 )
		sample[i]->clist[ testTerm + 1 ] = countCNV;
	      else if ( par::cnv_en_model == 2 )
		sample[i]->clist[ testTerm + 1 ] = avgKB;
	      else if ( par::cnv_en_model == 3 )
		sample[i]->clist[ testTerm + 1 ] = totalKB;
	      else if ( par::cnv_en_model == 4 )
		{
		  sample[i]->clist[ testTerm + 1 ] = countCNV;
		  sample[i]->clist[ testTerm + 2 ] = avgKB;
		}	     
	    }
	}
      else
	{	  
	  
	  // No CNVs at all

	  sample[i]->clist[ testTerm ] = 0;

	  if ( ! noCovar )
	    {
	      
	      if ( par::cnv_en_model == 1 )
		sample[i]->clist[ testTerm + 1 ] = 0;
	      else if ( par::cnv_en_model == 2 )
		sample[i]->clist[ testTerm + 1 ] = avg;
	      else if ( par::cnv_en_model == 3 )
		sample[i]->clist[ testTerm + 1 ] = 0;
	      else if ( par::cnv_en_model == 4 )
		{
		  sample[i]->clist[ testTerm + 1 ] = 0;
		  sample[i]->clist[ testTerm + 2 ] = avg;
		}	     
	    }
	}
      
    }  // Next individual
    

  //////////////////////////////////////////
  // Set up potential permutation test

  perm.setTests( totalTerms );
  perm.setPermClusters(*this);
  perm.originalOrder();
  vector_t original( totalTerms );


  ///////////////////////////////////////////
  // Perform association
  
  glmAssoc(false,*pperm);



  ///////////////////////////////////////////
  // Report results
  
  ofstream OUTF;
      
  bool valid = model->isValid();

  vector_t b = model->getCoefs();
  vector_t chisq(1,model->getStatistic());
  vector_t pval = model->getPVals();
  
  // NOTE: b includes intercept; pval doesn't
  
  double statistic = valid ? model->getStatistic() : 0;
  double pvalue = pval[ pval.size()-1 ];
  double beta = b[ b.size()-1 ];

  delete model;
  
  if ( print )
    {
      OUTF.open( f.c_str() , ios::out );

      OUTF << setw(12) << "TEST" << " "
	   << setw(12) << "BETA" << " "
	   << setw(12) << "P" << "\n";
      
      for (int c = 0; c < par::clist_number; c++)
	{
	  OUTF << setw(12) << clistname[c] << " "
	       << setw(12) << b[c+1] << " "
	       << setw(12) << pval[c] << "\n";
	}
      
      OUTF.close();
    }


  // Make test 1-sided, based on signed coef.

  for (int c = 0; c < par::clist_number; c++)
    original[c] = b[c+1];


  ////////////////////////////////////////////////////
  // Run permutations

  if ( par::permute ) 
    {    

      bool finished = false;
      int nv = 0;
      double z = 0;
      double unreal = 1/z;
      while(!finished)
	{      
	  perm.permuteInCluster();
	  glmAssoc(false,*pperm);
	  bool valid = model->isValid();
	  if ( ! valid ) { ++nv; }
	  vector_t b = model->getCoefs();
	  vector_t pr( par::clist_number );
	  for (int c = 0; c < par::clist_number; c++)
	    pr[c] = valid ? b[c+1] : unreal ;  
	  finished = perm.update(pr,original);      
	  delete model;
	}

      if (nv>0)
	cout << nv << " invalid models during permutation\n";
      
      if (!par::silent)
	cout << "\n\n";
      
      
      
      /////////////////////////////////////////////////
      // Display permuted results  
            
      f += ".mperm";
      printLOG("Writing CNV burden permutation results to [ "+f+" ]\n");
      ofstream FOUT;
      
      FOUT.open( f.c_str() , ios::out );
      FOUT.precision(4);
      
      FOUT << setw(8) << "TEST" << " "
	   << setw(12) << "EMP1" << "\n";

      
      for (int c = 0; c < par::clist_number; c++)
	{
	  FOUT << setw(8) << clistname[c] << " "
	       << setw(12) << perm.pvalue( c ) << "\n"; 
	       
	}           
      
      
      FOUT.close(); 
      
    }
  

  ///////////////////////////////////////
  // Some final tidying up

  par::assoc_glm_without_main_snp = OLD_assoc_glm_without_main_snp;
  par::clist = OLD_clist;  
  par::clist_number -= totalTerms;;
  clistname.resize( par::clist_number );
  for (int i=0; i<n; i++)
    sample[i]->clist.resize( par::clist_number );
  

  shutdown();
  return chisq;



}


