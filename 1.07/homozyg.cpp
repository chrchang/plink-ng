

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
#include <vector>
#include <map>
#include <iterator>
#include <cmath>

#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"
#include "stats.h"
#include "cnv.h"

extern Plink * PP;

namespace std
{
  template<>
  class less<Segment*> {
  public:
  bool operator()(Segment const* s1, Segment const* s2) const
  {
    if      ( s1->start  > s2->start ) return true;
    else if ( s1->start  < s2->start ) return false;
    else if ( s1->finish > s2->finish ) return true;
    else if ( s1->p1 < s2->p1 ) return true;
    else if ( s1->p1 > s2->p1 ) return false;
    else if ( s1->p2 < s2->p2 ) return true;
    else if ( s1->p2 > s2->p2 ) return false;
    else return false;
  }
  };
};


class Pool
{
public:
  set<Segment*,less<Segment*> > segs;
  vector<vector<int> > match;
  vector<int> matchcount;
  vector<int> group;
  vector<bool> index;
  int chr;
  int ng;
  int min;
  int max;
  int union_min;
  int union_max;
};


namespace std
{
  template<>
  class less<Pool*> {
  public:
  bool operator()(Pool const* p1, Pool const* p2) const
  {

    if      ( p1->segs.size() > p2->segs.size() ) return true;
    else if ( p1->segs.size() < p2->segs.size() ) return false;
    else
      {
	set<Segment*>::iterator s1 = p1->segs.begin();
	set<Segment*>::iterator s2 = p2->segs.begin();

	while ( s1 != p1->segs.end() )
	  {
	    if      ( (*s1)->start > (*s2)->start ) return true;
	    else if ( (*s1)->start < (*s2)->start ) return false;
	    else if ( (*s1)->finish > (*s2)->finish ) return true;
	    else if ( (*s1)->finish < (*s2)->finish ) return false;
	    s1++;
	    s2++;
	  }
      }
    return false;
  }
  };
};



vector_t compareCNVs(CNVIndivReport & a,
		     CNVIndivReport & b)
{
  vector_t res(8);

  // t1 = # events in cases
  // t2 = proportion of sample with 1+ event
  // t6 = proportion of sample with at least gene
  // t7 =

  if ( par::segment_test_1sided )
    {
      res[0] = a.t1 - b.t1;  // RATE
      res[1] = a.t2 - b.t2;  // PROP
      res[2] = a.t3 - b.t3;  // KBTOT
      res[3] = a.t4 - b.t4;  // KBAVG
      res[4] = a.t5 - b.t5;  // GRATE
      res[5] = a.t6 - b.t6;  // GPROP
      res[6] = a.t7 - b.t7;  // GRICH
      res[7] = a.t8 - b.t8;  // GRICH2
    }
  else
    {
      res[0] = fabs( a.t1 - b.t1 );
      res[1] = fabs( a.t2 - b.t2 );
      res[2] = fabs( a.t3 - b.t3 );
      res[3] = fabs( a.t4 - b.t4 );
      res[4] = fabs( a.t5 - b.t5 );
      res[5] = fabs( a.t6 - b.t6 );
      res[6] = fabs( a.t7 - b.t7 );
      res[7] = fabs( a.t8 - b.t8 );
    }

//    cout << "DETS: "
//         << a.count << " " << b.count << " -- G= "
//         << a.t5 << " " << b.t5 << " -- B= "
//         << a.t8 << " " << b.t8 << " -- EG= "
//         << a.t9 << " " << b.t9 << " EB= "
//         << a.t10 << " " << a.t10 << "\n";


  return res;

}



// Helper function
void summaryIndivSummaries(Plink * P,
			   int kmask,
			   map<indivPair,int> & segmentCount,
			   map<indivPair,double> & segmentLength,
			   CNVIndivReport & a,
			   CNVIndivReport & u,
			   vector_t & res)

{

  // Optionally only select individuals with kmask value,
  // unless kmask < 0

  // Return 2 (case/control) x object with
  //

  for ( int i = 0; i < P->n; i++)
    {

      Individual * person = P->sample[i];

      if ( kmask >= 0 )
	{
	  if ( person->sol != kmask )
	    continue;
	}

      indivPair t;
      t.p1 = t.p2 = P->sample[i];

      map<indivPair,int>::iterator ic = P->segmentCount.find(t);
      map<indivPair,double>::iterator il = P->segmentLength.find(t);
      map<indivPair,double>::iterator ic2 = P->segmentCount2.find(t);

      map<indivPair,double>::iterator ic2b;
      if ( par::cnv_count_baseline )
	ic2b = P->segmentCount2Baseline.find(t);

      //indivPair p = ic->first;

      if ( person->pperson->aff )
	++a.n;
      else
	++u.n;

      if ( ic == P->segmentCount.end() )
	continue;

      CNVIndivReport * pstat = person->pperson->aff ? &a : &u;


      // Basic CNV properties

      pstat->t1 += ic->second;
      pstat->t2++;
      pstat->t3 += il->second;
      pstat->t4 += il->second / (double)ic->second;


      // Geneset count statistics

      pstat->t5 += ic2->second;
      pstat->t6 += ic2->second > 0 ? 1 : 0;
      pstat->t7 += (double)ic2->second /  (double)il->second;
      pstat->t9 += PP->expectedOverlap[i];
      //      cout << " overlap = " << ic2->second << " of " << PP->expectedOverlap[i] << "\n";

      // Baseline/comparator geneset counts

      if ( par::cnv_count_baseline )
	{
	  if ( ic2b->second>0) pstat->count_baseline++;
	  pstat->t8 += ic2b->second;
	  pstat->t10 += PP->expectedOverlapBaseline[i];
	  //	  cout << "bline = " << ic2b->second << " " <<  PP->expectedOverlapBaseline[i] << "\n";
	}

      pstat->count++;


    } // Next individual


  // edit to make denom of t8 test total # of individuals
//   a.count_baseline = a.count;
//   u.count_baseline = u.count;


  // Save actual segment counts before we make the means

  a.segCount = (int)a.t1;
  u.segCount = (int)u.t1;


  ////////////////////////////////////////
  // Get means

  a.calculateResults();
  u.calculateResults();


  res = compareCNVs(a,u);

}


/////////////////////////////////////////////
// Entrypoint for all homozygosity run tests

void Plink::findAllHomozygousRuns(Perm & perm)
{

  if (par::SNP_major) SNP2Ind();

  string f = par::output_file_name + ".hom";

  // Calculate or read homozygous segments from a file?

  if ( ! par::read_segment_file )
    {

      ofstream HOM;
      HOM.open(f.c_str(),ios::out);
      HOM.precision(3);
      HOM.setf(ios::fixed);

      printLOG("\nWriting homozygosity-run information to [ "+f+" ] \n");
      printLOG("Run defined as: ");

      if (par::homo_run_kb)
	{
	  printLOG(int2str(par::homo_run_length_kb)+" kb ");
	  if (par::homo_run_snps) printLOG(", ");
	}

      if (par::homo_run_snps)
	printLOG(int2str(par::homo_run_length_snps)+" SNPs");
      printLOG("\n");

      //  printLOG("Allowing "+int2str(par::homo_run_het)+" hets per run\n");


      HOM << setw(par::pp_maxfid) << "FID" << " "
	  << setw(par::pp_maxiid) << "IID" << " "
	  << setw(8) << "PHE" << " "
	  << setw(4) << "CHR" << " "
	  << setw(par::pp_maxsnp) << "SNP1" << " "
	  << setw(par::pp_maxsnp) << "SNP2" << " "
	  << setw(12) << "POS1" << " "
	  << setw(12) << "POS2" << " "
	  << setw(10) << "KB" << " "
	  << setw(8) << "NSNP" << " "
	  << setw(8) << "DENSITY" << " "
	  << setw(8) << "PHOM" << " "
	  << setw(8) << "PHET" << "\n";


      printLOG("Homozygous segment criteria:\n");
      printLOG("  length (kb)       = " + int2str(par::homo_run_length_kb) + "\n");
      printLOG("  # SNPs (N)        = " + int2str(par::homo_run_length_snps) + "\n");
      printLOG("  density (kb/SNP)  = " + dbl2str(par::homo_run_density) + "\n");
      printLOG("  largest gap (kb)  = " + int2str(par::homo_run_gap) + "\n");


      // Rescale GAP to base-pairs

      par::homo_run_gap *= 1000;


      // Find segments for each individual

      for (int i1=0; i1<n; i1++)
	{
	  if (!par::silent)
	    cout << i1+1 << " of " << n << " individuals      \r";
	  //      findHomoRuns(sample[i1], HOM);
	  findHomoWindow(sample[i1], HOM);
	}

      if (!par::silent)
	cout << "\n\n";

      HOM.close();

    }
  else
    {
      // Read segments from a file

      ifstream SEG;
      SEG.open(par::read_segment_filename.c_str(),ios::in);
      printLOG("Reading homozygous-segment information from [ "
	       +par::read_segment_filename+" ]\n");
      checkFileExists(par::read_segment_filename);
      readHomozygSegmentFile(SEG);
      SEG.close();

    }



  ///////////////////////////////////
  // Display summary by individual

  printLOG("Writing segment summary to [ " + f +".indiv ]\n");
  indivSegmentSummary();

  ///////////////////////////////////
  // Display summary (i.e. per locus)

  printLOG("Writing segment summary to [ " + f +".summary ]\n\n");
  summaryIBSsegments(perm);


  //////////////////////////////////////
  // Look for shared homozygous regions
  // if --homo-match option given

  if (par::homo_summary_allelic_match)
    summariseHomoRuns();


}


// Homozygosity match, whole segments
bool segsMatch(Segment * s1, Segment * s2)
{

  if ( par::cnv_list )
    return true;

  int mismatch = 0;

  int start = s1->start > s2->start ? s1->start : s2->start;
  int finish = s1->finish < s2->finish ? s1->finish : s2->finish;

  // Assume individual major mode

  for (int l=start; l<= finish; l++)
    {
      // Only test homozygotes
      if ( s1->p1->one[l] == s1->p1->two[l] &&
	   s2->p1->one[l] == s2->p1->two[l] )
	if ( s1->p1->one[l] != s2->p1->one[l] )
	  mismatch++;
    }

  if ((double)mismatch/(double)(finish-start+1) > 1-par::fuzzy_homo ) return false;
  else return true;
}


// Homozygosity Match, supplying pool consensus start/finish
bool segsIBDMatchCON(Segment * s1, Segment * s2, int start, int finish)
{
  int match = 0;
  int valid = 0;

  // 1. Determine which allele is shared for each pair
  // 2. See whether this is the same allele between pairs

  // Assume individual major mode
    for (int l=start; l<= finish; l++)
    {

      // Ignore is 1+ missing or 2 hets

      bool a1 = s1->p1->one[l];
      bool a2 = s1->p1->two[l];

      bool b1 = s1->p2->one[l];
      bool b2 = s1->p2->two[l];

      bool c1 = s2->p1->one[l];
      bool c2 = s2->p1->two[l];

      bool d1 = s2->p2->one[l];
      bool d2 = s2->p2->two[l];

      bool allele1 = false;
      bool allele2 = false;

      // Any missing alleles?
      if ( (a1 && (!a2)) ||
	   (b1 && (!b2)) ||
	   (c1 && (!c2)) ||
	   (d1 && (!d2)) )
	continue;

      // Any double hets within each pair?
      if ( ((!a1) && a2) &&
	   ((!b1) && b2) )
	continue;

      if ( (!c1) && c2 &&
	   (!d1) && d2 )
	continue;

      // Any opposing homozygotes within pair?
      if ( (a1 == a2) &&
	   (b1 == b2) &&
	   (a1 != b1) )
	continue;

      if ( (c1 == c2) &&
	   (d1 == d2) &&
	   (c1 != d1) )
	continue;

      // Get the alleles for the pairs
      // (i.e. from the homozygote)

      if ( a1 == a2 ) allele1 = a1;
      else allele1 = b1;

      if ( c1 == c2 ) allele2 = c1;
      else allele2 = d1;

      // Do these match?
      if ( allele1 == allele2 )
	match++;

      valid++;
    }


  if ( valid<par::genome_test_min_snp ) return false;
  if ((double)match/(double)valid < par::fuzzy_homo ) return false;
  else return true;


}


 // Pairwise segmental match, whole segments
 bool segsIBDMatch(Segment * s1, Segment * s2)
{

  int match = 0;
  int valid = 0;

  int start = s1->start > s2->start ? s1->start : s2->start;
  int finish = s1->finish < s2->finish ? s1->finish : s2->finish;

   // 1. Determine which allele is shared for each pair
  // 2. See whether this is the same allele between pairs

  // Assume individual major mode

  for (int l=start; l<= finish; l++)
    {

      // Ignore is 1+ missing or 2 hets

      bool a1 = s1->p1->one[l];
      bool a2 = s1->p1->two[l];

      bool b1 = s1->p2->one[l];
      bool b2 = s1->p2->two[l];

      bool c1 = s2->p1->one[l];
      bool c2 = s2->p1->two[l];

      bool d1 = s2->p2->one[l];
      bool d2 = s2->p2->two[l];

      bool allele1 = false;
      bool allele2 = false;

      // Any missing alleles?
      if ( (a1 && (!a2)) ||
	   (b1 && (!b2)) ||
	   (c1 && (!c2)) ||
	   (d1 && (!d2)) )
	continue;

      // Any double hets within each pair?
      if ( ((!a1) && a2) &&
	   ((!b1) && b2) )
	continue;

      if ( (!c1) && c2 &&
	   (!d1) && d2 )
	continue;

      // Any opposing homozygotes within pair?
      if ( (a1 == a2) &&
	   (b1 == b2) &&
	   (a1 != b1) )
	continue;

      if ( (c1 == c2) &&
	   (d1 == d2) &&
	   (c1 != d1) )
	continue;

      // Get the alleles for the pairs
      // (i.e. from the homozygote)

      if ( a1 == a2 ) allele1 = a1;
      else allele1 = b1;

      if ( c1 == c2 ) allele2 = c1;
      else allele2 = d1;

      // Do these match?
      if ( allele1 == allele2 )
	match++;

      valid++;
    }


  if ( valid<par::genome_test_min_snp ) return false;
  if ((double)match/(double)valid < par::fuzzy_homo ) return false;
  else return true;
}


// Pairwise segmental match, supplying pool consensus start/finish
bool segsMatchCON(Segment * s1, Segment * s2, int start, int finish)
{

  if ( par::cnv_list )
    return true;

  int mismatch = 0;

  // Assume individual major mode

  for (int l=start; l<= finish; l++)
    {
      // Only test homozygotes
      if ( s1->p1->one[l] == s1->p1->two[l] &&
	   s2->p1->one[l] == s2->p1->two[l] )
	if ( s1->p1->one[l] != s2->p1->one[l] )
	  mismatch++;
    }

  if ((double)mismatch/(double)(finish-start+1) > 1-par::fuzzy_homo ) return false;
  else return true;
}


void displayPoolVerbose( Plink & P, Pool * pool , ofstream & OUT )
{

  // Figure out list of individuals, for the maximal region
  // i.e. union rather than intersection of the pool of segs

  set<Individual*> pset;
  vector<set<Individual*> > pgrpset( pool->ng );

  vector<Individual*> plist;
  vector<int> pstart;
  vector<int> pend;
  vector<int> pgrp;


  // Loop over each group in the pool
  for ( int g = 0; g < pool->ng; g++)
    {

      // Consider all segments in this pool
      set<Segment*>::iterator s = pool->segs.begin();
      int c2=0;
      while ( s != pool->segs.end() )
	{

	  // Not in group 'g' ?
	  if ( pool->group[c2] == g )
	    {

	      if ( pset.find( (*s)->p1 ) == pset.end() )
		{
		  pset.insert( (*s)->p1 );
		  plist.push_back( (*s)->p1 );
		  pstart.push_back( (*s)->start );
		  pend.push_back( (*s)->finish );
		  pgrp.push_back( g );
		}

	      if ( pset.find( (*s)->p2 ) == pset.end() )
		{
		  pset.insert( (*s)->p2 );
		  plist.push_back( (*s)->p2 );
		  pstart.push_back( (*s)->start );
		  pend.push_back( (*s)->finish );
		  pgrp.push_back( g );
		}

	      // Add to group-specfic unique individual list
	      if ( pgrpset[g].find( (*s)->p1 ) == pgrpset[g].end() )
		{
		  pgrpset[g].insert( (*s)->p1 );
		}

	      if ( pgrpset[g].find( (*s)->p2 ) == pgrpset[g].end() )
		{
		  pgrpset[g].insert( (*s)->p2 );
		}


	    }

	  // Next segment
	  s++;
	  c2++;
	  continue;

	}
    }

  // We now a have a unique list of individuals, in the same order in plist as
  // the grouping variable; now display, each row is a SNP

  // Header row (including IDs)

  OUT << setw(6) << " " << " "
      << setw(par::pp_maxfid) << "FID" << " "
      << setw(par::pp_maxiid) << "IID" << " "
      << setw(4) << "GRP" << " "
      << "\n";

  for (int i=0; i < plist.size(); i++)
    {
      OUT << setw(6) << int2str(i+1)+") " << " "
	  << setw(par::pp_maxfid) << plist[i]->fid << " "
	  << setw(par::pp_maxiid) << plist[i]->iid << "   ";

      // Display the group(s) this individual belongs to
      bool any = false;
      for ( int g = 1 ; g < pool->ng; g++ )
	{
	  if ( pgrpset[g].find(plist[i]) != pgrpset[g].end() )
	    {
	      if (!any) any = true;
	      else OUT << ", ";
	      OUT << g;
	    }
	}
      OUT << "\n";
    }

  OUT << "\n" << setw(par::pp_maxsnp) << "SNP" << " ";

  for (int i=0; i < plist.size(); i++)
    OUT << setw(5) << int2str(i+1)+" " << " ";
  OUT << "\n\n";


  ///////////////////////////////////////////////////////
  // Forcing a single pool? In this case, only display
  // the selected region, rather than the whole thing
  // (note: we have already displayed the union in the
  // *.overlap file, so okay to change this now

  if ( par::force_span )
    {
      pool->union_min = par::segment_snp1;
      pool->union_max = par::segment_snp2;
    }

  ///////////////////////////////////////////////////////
  // Display all SNPs

  for (int l = pool->union_min ; l <= pool->union_max ; l++)
    {
      if (! par::force_span ) if ( l == pool->min ) OUT << "\n";

      OUT << setw(par::pp_maxsnp) << P.locus[l]->name << " ";

      // Consider all individuals (in plist)
      for (int i=0; i < plist.size(); i++)
	{
	  if ( l >= pstart[i] && l <= pend[i] )
	    OUT << setw(5) << "["+genotype(P,plist[i],l)+"]" << " ";
	  else
	    OUT << setw(5) << " "+genotype(P,plist[i],l)+" " << " ";
	}
      OUT << "\n";

      if (! par::force_span ) if ( l == pool->max ) OUT << "\n";

    }

  OUT << "\n\n";


  ///////////////////////////////////////////////////////
  // Finally, find the consensus haplotype for each group

  // For example, pair within a group, for their overlap only
  // count number of greatest allele of the 4

  // Exactly the same principle applies whether looking at homozygous
  // or heterozygous runs; this will implicitly ignore double hets;
  // and provide a sensible way of dealing with missing data, hets in
  // homozygous runs, etc

  vector<string> ghap(pool->ng,"");

  for ( int g = 1; g < pool->ng; g++)
    {

      OUT << "Group " << g << "\n\n";

      for (int i=0; i < plist.size(); i++)
	{
	  // Is this individual in this group?
	  if ( pgrpset[g].find(plist[i]) != pgrpset[g].end() )
	    {
	      Individual * person = plist[i];
	      OUT << setw(6) << int2str(i+1)+") "
		  << setw(par::pp_maxfid) << person->fid << " "
		  << setw(par::pp_maxiid) << person->fid << " "
		  << setw(8) << person->phenotype << "\n";
	    }
	}

      OUT << "\n";
      OUT << "\n" << setw(par::pp_maxsnp) << "SNP" << " "
	  << setw(6) << " " << "  ";

      for (int i=0; i < plist.size(); i++)
	if ( pgrpset[g].find(plist[i]) != pgrpset[g].end() )
	  OUT << setw(5) << int2str(i+1)+" " << " ";
      OUT << "\n\n";


      // Consider each position
      for (int l = pool->union_min ; l <= pool->union_max ; l++)
	{

	  if (! par::force_span ) if ( l == pool->min ) OUT << "\n";
	  OUT << setw(par::pp_maxsnp) << P.locus[l]->name << " ";

	  // Loop over each group in the pool (start from 1..ng)
	  // (0 is unassigned)

	  // Keep track of most likely allele for this position for this group
	  int a1count = 0 ;
	  int a2count = 0 ;

	  // Consider all segments in this pool
	  set<Segment*>::iterator s = pool->segs.begin();
	  int c2=0;
	  while ( s != pool->segs.end() )
	    {

	      // Only count of this segment is spanning this particular
	      // position

	      if ( l >= (*s)->start && l <= (*s)->finish )
		{

		  // Not in group 'g' ?
		  if ( pool->group[c2] == g )
		    {

		      // Consider genotypes of p1, p2: get count of 4 alleles

		      bool a1 = (*s)->p1->one[l];
		      bool a2 = (*s)->p1->two[l];

		      bool b1 = (*s)->p2->one[l];
		      bool b2 = (*s)->p2->two[l];

		      // Homozygote?
		      if ( a1 == a2 )
			{
			  if ( a1 )
			    a2count++;
			  else
			    a1count++;
			}

		      if ( b1 == b2 )
			{
			  if ( b1 )
			    a2count++;
			  else
			    a1count++;
			}
		    }
		}

	      // Next segment
	      s++;
	      c2++;
	      continue;

	    } // Next segment in group


	  ///////////////////////////////////////////
	  // Display this group's consensus haplotype

	  string str = "?";
	  if ( a1count > a2count )
	    str = P.locus[l]->allele1;
	  else if ( a2count > a1count )
	    str = P.locus[l]->allele2;

	  OUT << setw(2) << str << " "
	      << setw(4) << " " << " ";

	  ghap[g] += str;

	  ///////////////////////////////////////////
	  // Consider all individuals (in plist)

	  for (int i=0; i < plist.size(); i++)
	    {
	      if ( pgrpset[g].find(plist[i]) != pgrpset[g].end() )
		//if ( pgrp[i] == g )
		{
		  if ( l >= pstart[i] && l <= pend[i] )
		    OUT << setw(5) << "["+genotype(P,plist[i],l)+"]" << " ";
		  else
		    OUT << setw(5) << " "+genotype(P,plist[i],l)+" " << " ";
		}
	    }

	  OUT << "\n";

	  if (! par::force_span ) if ( l == pool->max ) OUT << "\n";

	} // Next SNP

      OUT << "\n";

    } // Next group

  OUT << "\n\n";


  // Consider each position
  int lp=0;
  for (int l = pool->union_min ; l <= pool->union_max ; l++)
    {

      if (! par::force_span )  if ( l == pool->min ) OUT << "\n";
      OUT << setw(par::pp_maxsnp) << P.locus[l]->name << " ";

        for ( int g = 1; g < pool->ng; g++)
	  {
	    OUT << ghap[g][lp] << " ";
	  }

	lp++;

	OUT << "\n";
	if (! par::force_span )  if ( l == pool->max ) OUT << "\n";

    }


}


void Plink::groupSegmentsSpanning(int l)
{
  // Find list of segments spanning SNP 'l'
  // and return as a vector of ints -- this
  // is just a convenience wrapper around
  // summariseHomoRuns()

  par::force_span = true;
  par::segment_current_focal_snp = l;
  par::segment_silently_return_groups = true;
  bool old_silent = par::silent;
  par::silent = true;

  // A value of -1 means no segment

  indivSegmentGroup.resize(n);

  for (int i=0; i<n; i++)
    indivSegmentGroup[i] = -1;

  // This currently only works for homozygosity (i.e.  where each
  // individual can only belong to one group)

  // Populate Plink::indivSegmentGroup[]

  summariseHomoRuns();

  par::silent = old_silent;
}

void Plink::summariseHomoRuns()
{

  ofstream HOM;

  ////////////////////////////////
  // Output to screen, files, etc

  if ( ! par::segment_silently_return_groups )
    {


      if ( par::cnv_list )
	{
	  par::segment_overlap = false;
	  par::fuzzy_homo = 0;
	}

      // Codes:
      // CNVs -- par::cnv_list
      // IBD  -- par::segment_overlap
      // HBD  -- par::homo_run


      string f;

      if ( par::cnv_list )
	f = par::output_file_name + ".cnv.overlap";
      else if (par::segment_overlap)
	f = par::output_file_name + ".segment.overlap";
      else
	f = par::output_file_name + ".hom.overlap";

      HOM.open(f.c_str(),ios::out);
      HOM.precision(6);

      if ( par::cnv_list )
	printLOG("Writing CNV summary to [ "+f+" ] \n");
      else if (par::segment_overlap)
	printLOG("Writing IBD-run summary to [ "+f+" ] \n");
      else
	printLOG("Writing homozygosity-run summary to [ "+f+" ] \n");

      if (par::segment_verbose)
	printLOG("Writing IBD-run summary to [ "
		 + par::output_file_name
		 + ".segment.overlap.*.verbose"+" ] \n");
      else if (par::homozyg_verbose)
	printLOG("Writing homozygosity-run summary to [ "
		 + par::output_file_name
		 + ".hom.overlap.*.verbose"+" ] \n");

      if ( ! par::cnv_list)
	{
	  if (par::fuzzy_homo > 1e-6)
	    {
	      printLOG("Requiring allelic match, threshold "
		       + dbl2str(par::fuzzy_homo)+" identity\n");
	      printLOG("Requiring at least "
		       + int2str(par::genome_test_min_snp)
		       + " informative SNPs to match segments\n");
	    }
	  else
	    {
	      printLOG("Not requiring allelic match (overlapping only)\n");
	      par::fuzzy_homo = 0;
	    }
	}


      if (par::segment_overlap)
	HOM << setw(5) << "POOL" << " "
	    << setw(par::pp_maxfid) << "FID1" << " "
	    << setw(par::pp_maxiid) << "IID1" << " "
	    << setw(par::pp_maxfid) << "FID2" << " "
	    << setw(par::pp_maxiid) << "IID2" << " "
	    << setw(8) << "PHE" << " "
	    << setw(4) << "CHR" << " "
	    << setw(par::pp_maxsnp) << "SNP1" << " "
	    << setw(par::pp_maxsnp) << "SNP2" << " "
	    << setw(14) << "BP1" << " "
	    << setw(14) << "BP2" << " "
	    << setw(8) << "KB" << " "
	    << setw(8) << "NSNP" << " "
	    << setw(4) << "NSIM" << " "
	    << setw(6) << "GRP" << "\n";
      if ( par::cnv_list )
	HOM << setw(5) << "POOL" << " "
	    << setw(par::pp_maxfid) << "FID" << " "
	    << setw(par::pp_maxiid) << "IID" << " "
	    << setw(8) << "PHE" << " "
	    << setw(4) << "CHR" << " "
	    << setw(14) << "BP1" << " "
	    << setw(14) << "BP2" << " "
	    << setw(8) << "KB" << " "
	    << setw(6) << "TYPE" << " "
	    << setw(8) << "SCORE" << "\n";
      else
	HOM << setw(5) << "POOL" << " "
	    << setw(par::pp_maxfid) << "FID" << " "
	    << setw(par::pp_maxiid) << "IID" << " "
	    << setw(8) << "PHE" << " "
	    << setw(4) << "CHR" << " "
	    << setw(par::pp_maxsnp) << "SNP1" << " "
	    << setw(par::pp_maxsnp) << "SNP2" << " "
	    << setw(14) << "BP1" << " "
	    << setw(14) << "BP2" << " "
	    << setw(8) << "KB" << " "
	    << setw(8) << "NSNP" << " "
	    << setw(4) << "NSIM" << " "
	    << setw(6) << "GRP" << "\n";

    }



  // Make groups of segments
  // Group by maximally over-lapping region

  // Both irrespective of alleles
  // And considering alleles

  // i.e. all segments in a group must have at least one position
  // where they overlap Then group alleles within pool based on
  // allelic identity

  //   ----                  -----
  // ----                 ----------
  //   -----------------------
  //   |                     |




  // A unique list of overlapping segments
  set<Pool*,less<Pool*> > pools;


  /////////////////
  // 1. Make sets

  // Consider each position

  int start = 0;
  int finish = nl_all-1;

  if ( par::segment_silently_return_groups )
    start = finish = par::segment_current_focal_snp;
  else
    {
      if ( par::segment_m1 != "" )
	{
	  par::segment_snp1 = getMarkerNumber((*this),par::segment_m1);
	  if (par::segment_snp1==-1)
	    error("--segment-from {marker} not found");
	  start = par::segment_snp1;
	}

      if ( par::segment_m2 != "" )
	{
	  par::segment_snp2 = getMarkerNumber((*this),par::segment_m2);
	  if (par::segment_snp2==-1)
	    error("--segment-to {marker} not found");
	  finish = par::segment_snp2;
	}

      start = start > finish ? finish : start;
    }

  // A sequential position scan (i.e. to scan whole region,
  // potentially find non-overlapping pools, or force everything into
  // a single pool?


  if ( par::force_span )
    {

      Pool * thispool = new Pool;

      // Which segments contain this SNP?
      vector<Segment>::iterator s = segment.begin();
      while ( s != segment.end() )
	{

	  if ( s->start <= finish &&
	       s->finish >= start )
	    {
	      thispool->segs.insert( &(*s) );
	    }
	  s++;
	}


      // Add this pool to the overall list, if it is unique
      if ( thispool->segs.size() >= par::pool_size_min )
	{
	  pools.insert( thispool );
	}
      else
	{
	  delete thispool;

 	  if ( par::segment_silently_return_groups)
 	    return;

	  printLOG("No segments found in the spanned region, exiting now.\n");
	  return;
	}
    }
  else
    {
      for (int l=start; l<=finish; l++)
	{

	  if (!par::silent)
	    cout << "Considering position "
		 << l+1
		 << " of "
		 << nl_all
		 << "        \r";

	  Pool * thispool = new Pool;

	  // Which segments contain this SNP?
	  vector<Segment>::iterator s = segment.begin();
	  while ( s != segment.end() )
	    {

	      if ( s->start <= l &&
		   s->finish >= l )
		{
		  thispool->segs.insert( &(*s) );
		}
	      s++;
	    }


	  // Add this pool to the overall list, if it is unique
	  if ( thispool->segs.size() >= par::pool_size_min )
	    {
	      pools.insert( thispool );
	    }
	  else delete thispool;
	}

      if (!par::silent)
	cout << "\n";

      printLOG("Found "
	       +int2str(pools.size())
	       +" unique overlapping sets of segments\n");
    }


  //////////////////////////////////////////////////////
  // 2. Prune,  i.e.  drop AB if ABC and/or ABD exists?
  // Start at end (smallest set) and go upwards

  vector<bool> redundant(pools.size(),false);

  set<Pool*>::reverse_iterator p = pools.rbegin();
  int c=pools.size()-1;
  while ( p != pools.rend())
    {

      if (!par::silent)
	cout << c << " pools left to prune    \r";

      // Do we find all segments in this pool in a larger pool?

      set<Pool*>::iterator p2 = pools.begin();
      int c2=0;
      while ( p2 != pools.end() )
 	{

 	  // Already redundant
 	  if ( c2 >= c || redundant[c2] )
 	    {
 	      c2++;
 	      p2++;
 	      continue;
 	    }

 	  // Same pool
 	  if ( *p == *p2 )
 	    {
 	      c2++;
 	      p2++;
 	      continue;
 	    }

 	  set<Segment*>::iterator s = (*p)->segs.begin();
 	  bool embedded = true;

 	  while ( s != (*p)->segs.end() )
 	    {
 	      // Do we find this segment?

 	      set<Segment*>::iterator s2 = (*p2)->segs.begin();
 	      bool found = false;

 	      while ( s2 != (*p2)->segs.end() )
 		{
 		  if ( *s == *s2 )
 		    {
 		      found = true;
 		      break;
 		    }
 		  s2++;
 		}

 	      if ( ! found ) embedded = false;

 	      s++;
	    }

 	  if (embedded)
 	    {
 	      redundant[c] = true;
 	      break;
	    }

	  c2++;
 	  p2++;
 	}
      c--;
      p++;
    }


  if (!par::silent)
    cout << "\n";


  //////////////////////////////////////////////
  // 3. Determine allelic matching


  // Determine consensus region for each pool

  set<Pool*>::iterator pC = pools.begin();
  while ( pC != pools.end())
    {

      (*pC)->min=0;
      (*pC)->max=nl_all;

      (*pC)->union_min=nl_all;
      (*pC)->union_max=0;

      set<Segment*>::iterator s = (*pC)->segs.begin();
      while ( s != (*pC)->segs.end() )
	{
	  (*pC)->min = (*s)->start > (*pC)->min ? (*s)->start : (*pC)->min;
 	  (*pC)->max = (*s)->finish < (*pC)->max ? (*s)->finish : (*pC)->max;

	  (*pC)->union_min = (*s)->start < (*pC)->union_min ?
	    (*s)->start : (*pC)->union_min;
 	  (*pC)->union_max = (*s)->finish > (*pC)->union_max ?
	    (*s)->finish : (*pC)->union_max;

	  s++;

	}

      pC++;
    }

  // Find matches

  set<Pool*>::iterator pA = pools.begin();
  int cA=0;
  while ( pA != pools.end())
    {

      // Skip this pool?
      if ( redundant[cA] )
	{
 	  cA++;
 	  pA++;
 	  continue;
	}

      // For this pool, consider all segments in the pool in a pairwise
      // manner, and populate the match matrix

      set<Segment*>::iterator s1 = (*pA)->segs.begin();
      int c1=0;

      // Resize matching variables
      (*pA)->match.resize( (*pA)->segs.size() );
      for (int i=0; i< (*pA)->segs.size(); i++)
	(*pA)->match[i].clear();
      (*pA)->group.resize( (*pA)->segs.size() , 0 );

      (*pA)->matchcount.resize((*pA)->segs.size() , 0 );
      (*pA)->index.resize((*pA)->segs.size() , false );

      // Conisder each pair of segments
      while ( s1 != (*pA)->segs.end() )
	{

	  set<Segment*>::iterator s2 = (*pA)->segs.begin();
	  int c2=0;

	  while ( s2 != (*pA)->segs.end() )
	    {

	      if ( c2 >= c1 )
		{
		  s2++;
		  c2++;
		  continue;
		}

	      // Determine match function:
	      // Based on homozygosity or pairwise sharing
	      // Based on whole segments, or just pool consensus region

	      if (par::segment_overlap)
		{
		  // PAIRWISE SEGMENTAL MATCH

 		  if ( par::homo_run_consensus_match)
 		    {
 		      // Consensus match
 		      if ( segsIBDMatchCON( *s1, *s2, (*pA)->min, (*pA)->max ) )
 			{
 			  (*pA)->match[c1].push_back( c2 );
 			  (*pA)->match[c2].push_back( c1 );
 			  (*pA)->matchcount[c1]++;
 			  (*pA)->matchcount[c2]++;
 			}
 		    }
 		  else // else match whole segments (default)
 		    {
 		      if ( segsIBDMatch( *s1, *s2 ) )
 			{
 			  (*pA)->match[c1].push_back( c2 );
 			  (*pA)->match[c2].push_back( c1 );
 			  (*pA)->matchcount[c1]++;
 			  (*pA)->matchcount[c2]++;
 			}

 		    }
		}
	      else
		{
		  // HOMOZYGOSITY MATCH
		  // If a match, add pairs to lists
		  if ( par::homo_run_consensus_match)
		    {
		      // Consensus match
		      if ( segsMatchCON( *s1, *s2, (*pA)->min, (*pA)->max ) )
			{
			  (*pA)->match[c1].push_back( c2 );
			  (*pA)->match[c2].push_back( c1 );
			  (*pA)->matchcount[c1]++;
			  (*pA)->matchcount[c2]++;
			}
		    }
		  else // else match whole segments (default)
		    {
		      if ( segsMatch( *s1, *s2 ) )
			{
			  (*pA)->match[c1].push_back( c2 );
			  (*pA)->match[c2].push_back( c1 );
			  (*pA)->matchcount[c1]++;
			  (*pA)->matchcount[c2]++;
			}

		    }
		}

	      // Next segment (B)
	      s2++;
	      c2++;
	    }

	  // Next segment (A)
	  s1++;
	  c1++;
	}


      // Parse the list

      bool done = false;
      (*pA)->ng = 1;

      while ( ! done )
	{

	  // Find largest, ungrouped list
	  int maxlist = 0;
	  int maxlisti = -1;

	  for (int i=0; i < (*pA)->group.size() ; i++)
	    {
	      if ( (*pA)->group[i] == 0 )
		{
		  if ( (*pA)->match[i].size() >= maxlist )
		    {
		      maxlist = (*pA)->match[i].size();
		      maxlisti = i;
		    }
		}
	    }

	  // Set group for this, and matches
	  (*pA)->index[maxlisti] = true;
	  (*pA)->group[maxlisti] = (*pA)->ng;
	  for ( int j=0; j < maxlist; j++)
	    (*pA)->group[ (*pA)->match[maxlisti][j] ] = (*pA)->ng;

	  // Advance to next group
	  (*pA)->ng++;

	  // Are all segments grouped?
	  bool ungroup = false;

	  for (int i=0; i < (*pA)->group.size() ; i++)
	    {
	      if ( (*pA)->group[i] == 0 )
		{
		  ungroup = true;
		  //break;
		}
	    }

	  if ( ! ungroup ) done = true;

	}


      // Remove temporary storage
      (*pA)->match.clear();

      // Next pool
      pA++;
      cA++;
    }


  //////////////////////////////////////////
  // Populate indivSegmentGroup and return?

  if ( par::segment_silently_return_groups )
    {

      // Should only be a single pool -- but for now, let's ignore and
      // just use this code... check later...

      set<Pool*>::iterator p2 = pools.begin();
      c=0;
      while ( p2 != pools.end())
	{

	  if (redundant[c])
	    {
	      c++;
	      p2++;
	      continue;
	    }

	  // Loop over each group in the pool
	  for ( int g = 0; g < (*p2)->ng; g++)
	    {
	      // Consider all segments in this pool
	      set<Segment*>::iterator s = (*p2)->segs.begin();
	      int c2=0;
	      while ( s != (*p2)->segs.end() )
		{

		  // Not in group 'g' ?
		  if ( (*p2)->group[c2] != g )
		    {
		      s++;
		      c2++;
		      continue;
		    }

		  // Find person, given pointer.... yes, yes, I know
		  // this is a terrible way to do things... for now
		  // just assume homozygous segments (i.e. p1==p2)

		  for (int i=0; i<n; i++)
		    {
		      if ( indivSegmentGroup[i] == -1 &&
			   sample[i] == (*s)->p1 )
			indivSegmentGroup[i] = (*p2)->group[c2];
		    }

		  s++;
		  c2++;

		} // Next segment
	    } // Next group

	  c++;
	  p2++;

	} // Next pool


      return;
    }


  //////////////
  // 4. Display

  set<Pool*>::iterator p2 = pools.begin();
  c=0;
  while ( p2 != pools.end())
    {

      if (redundant[c])
	{
 	  c++;
  	  p2++;
  	  continue;
	}

      int ncase=0;
      int ncontrol=0;

      // Loop over each group in the pool
      for ( int g = 0; g < (*p2)->ng; g++)
	{

	  // Consider all segments in this pool
	  set<Segment*>::iterator s = (*p2)->segs.begin();
	  int c2=0;
	  while ( s != (*p2)->segs.end() )
	    {

	      // Not in group 'g' ?
	      if ( (*p2)->group[c2] != g )
		{
		  s++;
		  c2++;
		  continue;
		}

	      HOM << setw(5) << "S"+int2str(c+1) << " ";

	      if (par::segment_overlap)
		{

		  HOM << setw(par::pp_maxfid) << (*s)->p1->fid << " "
		      << setw(par::pp_maxiid) << (*s)->p1->iid << " "
		      << setw(par::pp_maxfid) << (*s)->p2->fid << " "
		      << setw(par::pp_maxiid) << (*s)->p2->iid << " ";
		  if (par::bt)
		    {
		      if ( (*s)->p1->missing || (*s)->p2->missing )
			HOM << setw(8) << "NA" << " ";
		      else if ( (!(*s)->p1->aff) && (!(*s)->p2->aff) )
			{
			  HOM << setw(8) << "-1" << " ";
			  ncontrol++;
			}
		      else if ( (*s)->p1->aff && (*s)->p2->aff )
			{
			  HOM << setw(8) << "1" << " ";
			  ncase++;
			}
		      else
			{
			  HOM << setw(8) << "0" << " ";
			  ncontrol++;
			}
		    }
		  else
		    HOM << setw(8) << "NA" << " ";
		}
	      else
		{
		  HOM << setw(par::pp_maxfid) << (*s)->p1->fid << " "
		      << setw(par::pp_maxiid) << (*s)->p1->iid << " ";

		  if (par::bt)
		    HOM << setw(8) << (*s)->p1->phenotype << " ";
		  else
		    {
		      HOM.precision(4);
		      HOM << setw(8) << (*s)->p1->phenotype << " ";
		      HOM.precision(8);
		    }

		  if ( (*s)->p1->aff )
		    ncase++;
		  else
		    ncontrol++;
		}

	      HOM << setw(4) << locus[(*s)->start]->chr << " ";

	      if ( ! par::cnv_list )
		HOM << setw(par::pp_maxsnp) << locus[(*s)->start]->name << " "
		    << setw(par::pp_maxsnp) << locus[(*s)->finish]->name << " ";

	      HOM << setw(14) << locus[(*s)->start]->bp << " "
		  << setw(14) << locus[(*s)->finish]->bp << " "
		  << setw(8) << (double)(locus[(*s)->finish]->bp
					 - locus[(*s)->start]->bp ) / 1000.0 << " ";

	      // CNV-specific info: type and score
	      if ( par::cnv_list )
		{
		  if ( (*s)->type == 1 )
		    HOM << setw(6) << "DEL" << " "
			<< setw(8) << (*s)->score << " ";
		  else
		    HOM << setw(6) << "DUP" << " "
			<< setw(8) << (*s)->score << " ";
		}

	      // Non-CNV specific info
	      if ( !par::cnv_list )
		{
		  HOM<< setw(8) << (*s)->finish - (*s)->start + 1 << " "
		     << setw(4) << (*p2)->matchcount[c2] << " ";
		  if ( (*p2)->index[c2] )
		    HOM << setw(6) << int2str((*p2)->group[c2])+"*" << " ";
		  else
		    HOM << setw(6) << int2str((*p2)->group[c2])+" " << " ";
		}

	      HOM << "\n";

	      s++;
	      c2++;

	    } // Next segment

	} // Next group

      // Consensus region
      if (! par::force_span )
	{
	  HOM << setw(5) << "S"+int2str(c+1) << " ";
	  HOM << setw(par::pp_maxfid) << "CON" << " "
	      << setw(par::pp_maxiid) <<  (*p2)->segs.size() << " ";
	  if (par::segment_overlap)
	    HOM << setw(par::pp_maxfid) << "NA" << " "
		<< setw(par::pp_maxiid) << "NA" << " ";
	  HOM << setw(8) << int2str(ncase)+":"+int2str(ncontrol) << " "
	      << setw(4) << locus[(*(*p2)->segs.begin())->start]->chr << " ";

	  if ( ! par::cnv_list )
	    HOM << setw(par::pp_maxsnp) << locus[(*p2)->min]->name << " "
		<< setw(par::pp_maxsnp) << locus[(*p2)->max]->name << " ";

	  HOM << setw(14) << locus[(*p2)->min]->bp << " "
	      << setw(14) << locus[(*p2)->max]->bp << " "
	      << setw(8) << (double)(locus[(*p2)->max]->bp
				     - locus[(*p2)->min]->bp ) / 1000.0 << " ";

	  if ( par::cnv_list )
	    {
	      HOM << setw(6) << "NA" << " "
		  << setw(8) << "NA" << " ";
	    }

	  if ( ! par::cnv_list )
	    {
	      HOM << setw(8) << (*p2)->max - (*p2)->min + 1;
	      HOM << setw(6) << "NA" << " "
		  << setw(6) << "NA" << " ";
	    }

	  HOM << "\n";

	}
      else
	{
	  HOM << setw(5) << "S"+int2str(c+1) << " ";
	  HOM << setw(par::pp_maxfid) << "FORCE" << " "
	      << setw(par::pp_maxiid) <<  (*p2)->segs.size() << " ";
	  if (par::segment_overlap)
	    HOM << setw(par::pp_maxfid) << "NA" << " "
		<< setw(par::pp_maxiid) << "NA" << " ";
	  HOM << setw(8) << int2str(ncase)+":"+int2str(ncontrol) << " "
	      << setw(4) << locus[(*(*p2)->segs.begin())->start]->chr << " ";

	  if ( ! par::cnv_list )
	    HOM << setw(par::pp_maxsnp) << locus[par::segment_snp1]->name << " "
		<< setw(par::pp_maxsnp) << locus[par::segment_snp2]->name << " ";

	  HOM << setw(14) << locus[par::segment_snp1]->bp << " "
	      << setw(14) << locus[par::segment_snp2]->bp << " "
	      << setw(8) << (double)(locus[par::segment_snp2]->bp
				     - locus[par::segment_snp1]->bp ) / 1000.0 << " ";
	  if ( ! par::cnv_list )
	    {
	      HOM << setw(8) << par::segment_snp2 - par::segment_snp1 + 1;
	      HOM << setw(6) << "NA" << " "
		  << setw(6) << "NA" << " ";
	    }
	  HOM<< "\n";

	}

      // Union region
      HOM << setw(5) << "S"+int2str(c+1) << " ";
      HOM << setw(par::pp_maxfid) << "UNION" << " "
	  << setw(par::pp_maxiid) <<  (*p2)->segs.size() << " ";
      if (par::segment_overlap)
	HOM << setw(par::pp_maxfid) << "NA" << " "
	    << setw(par::pp_maxiid) << "NA" << " ";
      HOM << setw(8) << int2str(ncase)+":"+int2str(ncontrol) << " "
	  << setw(4) << locus[(*(*p2)->segs.begin())->start]->chr << " ";

      if ( ! par::cnv_list )
	HOM << setw(par::pp_maxsnp) << locus[(*p2)->union_min]->name << " "
	    << setw(par::pp_maxsnp) << locus[(*p2)->union_max]->name << " ";

      HOM << setw(14) << locus[(*p2)->union_min]->bp << " "
	  << setw(14) << locus[(*p2)->union_max]->bp << " "
	  << setw(8) << (double)(locus[(*p2)->union_max]->bp
				 - locus[(*p2)->union_min]->bp ) / 1000.0 << " ";

      if ( par::cnv_list )
	{
	  HOM << setw(6) << "NA" << " "
	      << setw(8) << "NA" << " ";
	}

      if ( !par::cnv_list )
	HOM << setw(8) << (*p2)->union_max - (*p2)->union_min + 1
	    << setw(6) << "NA" << " "
	    << setw(6) << "NA" << " ";

      HOM << "\n\n";

      // Verbose mode? (not for CNV lists)

      if ( ( par::homozyg_verbose || par::segment_verbose )  && !par::cnv_list )
	{
	  string f;
	  if (par::segment_overlap)
	    f = par::output_file_name + ".segment.overlap.S"+int2str(c+1)+".verbose";
	  else
	    f = par::output_file_name + ".hom.overlap.S"+int2str(c+1)+".verbose";

	  ofstream VHOM(f.c_str(),ios::out);
	  VHOM.precision(2);
	  displayPoolVerbose( *this , *p2 , VHOM );
	  VHOM.close();
	}

      c++;
      p2++;

    } // Next pool


  HOM.close();
}




void Plink::findHomoRuns(Individual * person, ofstream & HOM)
{

  int l=0;
  int lasthom=0;
  int nmiss = 0;
  int nhet = 0;
  bool run = false;
  int start = 0;
  int end = 0;

  while ( l < nl_all )
    {

      // Skip haploid chromosomes / end any existing run
      if ( ( par::chr_sex[locus[l]->chr] && person->sex ) ||
	   par::chr_haploid[locus[l]->chr] )
	{
	  if (run)
	    {
	      end = l-1;
	      run = false;
	    }
	  else
	    {
	      l++;
	      continue;
	    }
	}

      // Outside of a run?
      if (!run)
	{
	  // A new run?
	  if (person->one[l] == person->two[l])
	    {
	      start = lasthom = l;
	      nmiss=0;
	      nhet=0;
	      run=true;
	    }
	}
      else // if already in a run, either end or increase length?
	{

	  if ( locus[l]->chr != locus[start]->chr ) // different chromosome?
	    {
	      end = l-1;
	      run = false;
	    }
	  else if ( l == (nl_all -1) ) // or end of all SNPs?
	    {
	      if ( person->one[l] == person->two[l] )
		lasthom=l;
	      end = lasthom;
	      run = false;
	    }
	  // found a het?
	  else if ( (!person->one[l]) && person->two[l])
	    {
	      if (nhet==par::homo_run_het)
		{
		  end = lasthom;
		  run = false;
		}
	      else
		nhet++;
	    }
	  else // ...continue run
	    {
	      lasthom=l;
	      if ( person->one[l] && (!person->two[l]) )
		nmiss++;
	    }

	}


      // Check run length?
      if (!run)
	{

	  bool accept = true;

	  if (par::homo_run_kb)
	    if ( locus[end]->bp - locus[start]->bp <
		 par::homo_run_length_kb * 1000 )
	      accept = false;

	  if (par::homo_run_snps)
	    if ( end - start +1 < par::homo_run_length_snps  )
	      accept = false;

	  if (accept)
	    {

	      HOM << setw(par::pp_maxfid) << person->fid << " "
		  << setw(par::pp_maxiid) << person->iid << " "
		  << setw(4) << locus[start]->chr << " "
		  << setw(par::pp_maxsnp) << locus[start]->name << " "
		  << setw(par::pp_maxsnp) << locus[end]->name << " "
		  << setw(12) <<locus[start]->bp << " "
		  << setw(12) <<locus[end]->bp << " "
		  << setw(10) << (double)(locus[end]->bp
					  - locus[start]->bp)/(double)1000 << " "
		  << setw(10) << end - start + 1 << " "
		  << setw(4) << nhet << " "
		  << setw(4) << nmiss << "\n";

	      Segment s;
	      s.p1 = s.p2 = person;
	      s.start = start;
	      s.finish = end;
	      segment.push_back(s);
	    }


	  //////////////////
	  // Clear counters

	  start = end = nmiss = 0;

	}

      ///////////////
      // Next locus

      l++;
    }

}


class HWindow {
public:

  int start, stop;
  bool leftHomozyg, leftMissing;
  bool rightHomozyg, rightMissing;
  bool finished;
  bool valid;

  Individual * person;
  Plink * P;

  int homCount, hetCount, misCount;

  // Constructor: must specify a person
  HWindow(Plink*,Individual*);

  // Set Window bounaries
  void set(int,int);

  // Full count
  void recount();

  // Shift update
  void shift();

};

HWindow::HWindow(Plink * plink, Individual * p)
{
  P = plink;
  person = p;
  finished = false;
  valid = true;
  start = stop = 0;
}

void HWindow::recount()
{
  // Assume individual-major mode

  // Reset counts
  homCount = hetCount = misCount = 0;

  for (int l = start; l <= stop; l++)
    {

      if ( person->one[l] )
	{
	  if ( person->two[l] )
	    homCount++;
	  else
	    misCount++;
	}
      else
	{
	  if ( person->two[l] )
	    hetCount++;
	  else
	    homCount++;
	}
    }

}

void HWindow::shift()
{

  // Find a new, valid (i.e. all on same autosomal chromosome)
  // window

  valid = false;

  bool moreThanOne = false;

  while ( ! valid )
    {
      set( ++start, ++stop );
      if ( ! valid )
	moreThanOne = true;
      if ( finished )
	return;
    }

  // Typically, we will just shift a single SNP, so we do not
  // need to recount eveything

  if ( ! moreThanOne )
    {
      // Update counts: remove left edge

      int trailing = start - 1;

      if ( person->one[trailing] == person->two[trailing] )
	--homCount;
      else if ( person->one[trailing] )
	--misCount;
      else
	--hetCount;

      // Add right edge; leading edge is now 'stop'

      if ( person->one[stop] == person->two[stop] )
	++homCount;
      else if ( person->one[stop] )
	++misCount;
      else
	++hetCount;
    }
  else
    {
      // A full recount
      recount();
    }

}

void HWindow::set(int s1, int s2)
{

  start = s1;
  stop = s2;

  // Can we set this window?

  if ( s1 < 0 || s2 >= P->nl_all || s1 > s2 ||
       P->locus[s1]->chr != P->locus[s2]->chr ||
       ( par::chr_sex[P->locus[s1]->chr] && person->sex ) ||
       par::chr_haploid[P->locus[s1]->chr] )
    {

      // Finished all SNPS?
      if ( s2 == P->nl_all )
	finished = true;

      valid = false;
      return;
    }

  // Set poisitions

  valid = true;

  return;
}


void Plink::findHomoWindow(Individual * person, ofstream & HOM)
{

  // Window properties

  // Only count a window if not too many missing genotypes
  // and if distance spanned if not too great. Then score
  // as 0 or 1 depending on whether we see too many hets.

  // Record homozygosity state

  vector<double> totalWindows(nl_all,0);
  vector<double> homozygWindows(nl_all,0);

  // Create an initial window and place just before
  // first SNP

  HWindow window(this,person);
  window.set(0,par::homo_windowSize-1);
  if ( ! window.valid )
    window.shift();
  window.recount();

  while ( 1 )
    {

      // End of genome

      if ( window.finished )
	break;

      // A valid window? Then record

      if ( window.valid )
	{

	  // Is this also homozygous enough ?

	  bool homozyg = false;

	  if ( window.misCount <= par::homo_windowAllowedMissing &&
	       window.hetCount <= par::homo_windowAllowedHet )
	    homozyg = true;

	  // Score for all SNPs

	  for ( int l = window.start ;
		l <= window.stop ;
		l++)
	    {
	      totalWindows[l]++;

	      if (homozyg)
		homozygWindows[l]++;
	    }
	}


      // Move window

      window.shift();

    }


  // Extract segments -- above threshold values

  for ( int l = 0 ; l < nl_all ; l++ )
    {
      if ( totalWindows[l] == 0 )
	homozygWindows[l] == 0;
      else
	homozygWindows[l] /= totalWindows[l];
    }

  if (par::verbose)
    HOM << "Segments for "
	<< person->fid << " "
	<< person->iid << "\n";

  int bp1 = 0;
  int l1 = 0;

  // Find segments

  int allSegs = 0;
  int incSegs = 0;

  bool inseg = false;

  for ( int l = 0 ; l < nl_all ; l++ )
    {

      if ( ( !inseg ) && homozygWindows[l] >= par::homo_threshold )
	{
	  inseg = true;
	  bp1 = locus[l]->bp;
	  l1 = l;
	}

      else if ( inseg )
	{

	  bool ending = homozygWindows[l] < par::homo_threshold;
	  bool bigGap = ( locus[l]->bp - locus[l-1]->bp ) > par::homo_run_gap;
	  bool newChr = locus[l]->chr != locus[l-1]->chr;
	  bool lastSNP = l == ( nl_all - 1) ;

	  if ( ending || bigGap || newChr || lastSNP )
	    {
	      inseg = false;

	      // Does this potential segment shape up?
	      // Length, number of SNPs, density of SNPs,
	      // largest gap

	      int l2 = l - 1;

	      // Check this is not the final, homozygous SNP

	      if ( lastSNP && ! ending )
		l2 = l;

	      // We might want to start a new segment on this SNP also

	      if ( ( newChr || bigGap ) && ! ending )
		--l;

	      double length = ( locus[l2]->bp - bp1 ) / 1000.0;
	      int snps = l2 - l1 + 1;
	      double density = length / (double)snps ;

	      if ( length >= par::homo_run_length_kb &&
		   snps >= par::homo_run_length_snps &&
		   length/(double)snps <= par::homo_run_density )
		{

		  incSegs++;

		  // Some sanity checks
		  double proHet = 0;
		  double proHom = 0;

		  for (int j= l1; j <= l2; j++)
		    {
		      bool s1 = person->one[j];
		      bool s2 = person->two[j];

		      if ( s1 == s2 )
			++proHom;
		      else if ( s2 )
			++proHet;
		    }

		  proHom /= (double)snps;
		  proHet /= (double)snps;


		  HOM << setw(par::pp_maxfid) << person->fid << " "
		      << setw(par::pp_maxiid) << person->iid << " "
		      << setw(8) << person->phenotype << " "
		      << setw(4) << locus[l1]->chr << " "
		      << setw(par::pp_maxsnp) << locus[l1]->name << " "
		      << setw(par::pp_maxsnp) << locus[l2]->name << " "
		      << setw(12) << locus[l1]->bp << " "
		      << setw(12) << locus[l2]->bp << " "
		      << setw(10) << length << " "
		      << setw(8) << snps << " "
		      << setw(8) << density << " "
		      << setw(8) << proHom << " "
		      << setw(8) << proHet << "\n";

		  Segment s;
		  s.p1 = s.p2 = person;
		  s.start = l1;
		  s.finish = l2;
		  segment.push_back(s);

		}
	      else if ( par::verbose )
		{

		  		  // Some sanity checks
		  double proHet = 0;
		  double proHom = 0;

		  for (int j= l1; j <=l2; j++)
		    {
		      bool s1 = person->one[j];
		      bool s2 = person->two[j];

		      if ( s1 == s2 )
			++proHom;
		      else if ( s2 )
			++proHet;
		    }

		  proHom /= (double)snps;
		  proHet /= (double)snps;


		  HOM <<"* " << setw(par::pp_maxfid) << person->fid << " "
		      << setw(par::pp_maxiid) << person->iid << " "
		      << setw(4) << person->phenotype << " "
		      << setw(4) << locus[l1]->chr << " "
		      << setw(par::pp_maxsnp) << locus[l1]->name << " "
		      << setw(par::pp_maxsnp) << locus[l2]->name << " "
		      << setw(12) << locus[l1]->bp << " "
		      << setw(12) << locus[l2]->bp << " "
		      << setw(10) << length << " "
		      << setw(8) << snps << " "
		      << setw(8) << density << " "
		      << setw(8) << proHom << " "
		      << setw(8) << proHet << "\n";


		}

	      allSegs++;

	    }
	}
    }

   if ( par::verbose )
     {
       for ( int l = 0 ; l < nl_all ; l++ )
 	{
 	  HOM << "SX " << locus[l]->chr << "\t"
	      << locus[l]->name << "\t"
	      << (double)locus[l]->bp / (1000.0*1000.0) << "\t"
	      << homozygWindows[l] << "\t";
 	  if ( homozygWindows[l] >= par::homo_threshold ) HOM << "1\t";
 	  else HOM << "0\t";

 	  if ( person->one[l] == person->two[l] ) HOM << "1\n";
 	  else if ( person->two[l] ) HOM << "0\n"; //het
 	  else HOM << "-1\n";

 	}
       HOM << "\n";
     }

}



bool segsOverlap(Segment * s1, Segment * s2)
{
  if ( s1->finish < s2->start ) return false;
  else if ( s2->finish < s1->start ) return false;
  else return true;
}




void Plink::homozygousSegmentPermutationTest(Perm & perm,
					     string f,
					     vector<int> & coverage_aff,
					     vector<int> & coverage_unaff )
{


  // Permutation test for excess of case homozygous segments in a
  // particular region (one-sided test)

  // Also applies for CNV data

  // Optionally allowed for this to operate on smoothed data (i.e.
  // average of event count over a KB window, forwards and backwards
  // from the given position)

  // Also, adds in summary statistics permutation

  double tot_aff=0;
  double tot_not=0;

  for (int i1=0; i1<n; i1++)
    if ( sample[i1]->aff ) tot_aff++;
    else tot_not++;

  printLOG(int2str((int)tot_aff)+" affected individuals out of "
	   +int2str(int(tot_not+tot_aff))+" in total\n");



  //////////////////////////////////////////
  // Test positons = MAP positions (nl_all)
  // Test positions = summed segment counts ( get from original counts )
  // Test position = aggregate statistics ( 7 tests)


  int nt = nl_all;

  if ( par::seg_test_region )
    nt = coverage_aff.size();
  else if ( par::cnv_indiv_perm )
    {
      nt = 7;
      if ( par::cnv_count_baseline )
	++nt;
    }


  // Option per-individual summary tests? (4 tests)

  // Cases - controls: total # segs
  //                   # people w/ 1+ seg
  //                   total kb length
  //                   mean segment length

  //                   gene-count
  //                   atleast-1-gene-count
  //                   gene-enrichment


  perm.setTests(nt);
  perm.setPermClusters(*this);
  perm.originalOrder();

  vector<double> original(nt);


  //////////////////////////////////////////////////
  // Genic/regional, or standard positional tests?

  if ( ! par::cnv_indiv_perm )
    {

//       // Get genome-wide means

//       double Ag = 0;
//       double Ug = 0;
//       for (int l=0; l<nt; l++)
// 	{
// 	  Ag += coverage_aff[l];
// 	  Ug += coverage_unaff[l];
// 	}
//       Ag /= (double)nt;
//       Ug /= (double)nt;

      for (int l=0; l<nt; l++)
	{

	  double A = coverage_aff[l];
	  double U = coverage_unaff[l];

	  //       // Correct for genome-wide average
	  //       A = A-Ag < 0 ? 0 : A - Ag;
	  //       U = U-Ug < 0 ? 0 : U - Ug;
	  //       double statistic = A / (double)tot_aff  -  U / (double)tot_not ;
	  //       original[l] = statistic < 0 ? 0 : statistic;

	  // Simple positional chi-sq
	  double A0 = tot_aff - A;
	  double U0 = tot_not - U;
	  double N = tot_aff + tot_not;
// 	  double rowHOM = A+U;
// 	  double rowNOM = A0+U0;

	  double E_A = ( tot_aff * (A+U) ) / N;
	  double E_A0 = ( tot_aff * (A0+U0) ) / N;
	  double E_U = ( tot_not * (A+U) ) / N;
	  double E_U0 = ( tot_not * (A0+U0) ) / N;

	  if ( E_A == 0 && E_U == 0 )
	    original[l] = 0;
	  else
	    {
	      if ( par::segment_test_1sided && A/tot_aff < U/tot_not )
		{
		  original[l] = 0;
		}
	      else
		{

		  original[l] =  ((A-E_A)*(A-E_A) ) / E_A +
		    ((A0-E_A0)*(A0-E_A0))/E_A0 +
		    ((U-E_U)*(U-E_U))/E_U +
		    ((U0-E_U0)*(U0-E_U0))/E_U0 ;

		  if ( ! realnum(original[l]) )
		    original[l] = 0;
		}
	    }

	}
    }


  ///////////////////////////
  // Individual summary tests

  if ( par::cnv_indiv_perm )
    {

      ofstream OUT;
      OUT.open( (par::output_file_name+".cnv.grp.summary").c_str(), ios::out);
      OUT.precision(4);

      printLOG("Writing group summary statistics to [ "
	       + par::output_file_name
	       + ".cnv.grp.summary ]\n");

      OUT << setw(12) << "TEST" << " "
	  << setw(8) << "GRP" << " "
	  << setw(12) << "AFF" << " "
	  << setw(12) << "UNAFF" << "\n";


      // Perform for all individuals
      CNVIndivReport a, u;

      summaryIndivSummaries( this, -1, segmentCount, segmentLength, a, u, original );


      OUT << setw(12) << "N" << " "
	  << setw(8) << "ALL" << " "
	  << setw(12) << a.segCount << " "
	  << setw(12) << u.segCount << "\n";

      OUT << setw(12) << "RATE" << " "
	  << setw(8) << "ALL" << " "
	  << setw(12) << a.t1 << " "
	  << setw(12) << u.t1 << "\n";

      OUT << setw(12) << "PROP" << " "
	  << setw(8) << "ALL" << " "
	  << setw(12) << a.t2 << " "
	  << setw(12) << u.t2 << "\n";

      OUT << setw(12) << "TOTKB" << " "
	  << setw(8) << "ALL" << " "
	  << setw(12) << a.t3 << " "
	  << setw(12) << u.t3 << "\n";

      OUT << setw(12) << "AVGKB" << " "
	  << setw(8) << "ALL" << " "
	  << setw(12) << a.t4 << " "
	  << setw(12) << u.t4 << "\n";

      if ( par::cnv_count )
	{
	  OUT << setw(12) << "GRATE" << " "
	      << setw(8) << "ALL" << " "
	      << setw(12) << a.t5 << " "
	      << setw(12) << u.t5 << "\n";

 	  OUT << setw(12) << "GPROP" << " "
 	      << setw(8) << "ALL" << " "
 	      << setw(12) << a.t6 << " "
 	      << setw(12) << u.t6 << "\n";


 	  OUT << setw(12) << "GRICH" << " "
 	      << setw(8) << "ALL" << " "
 	      << setw(12) << a.t7 << " "
 	      << setw(12) << u.t7 << "\n";


	  ////////////////////////////////////////////////////////////////////////////
	  // Second enrichment test: # of genes of interest / # of reference gene set

	  if ( par::cnv_count_baseline )
	    {
	      OUT << setw(12) << "GRICH2" << " "
		  << setw(8) << "ALL" << " "
		  << setw(12) << a.t8 << " "
		  << setw(12) << u.t8 << "\n";
	    }
	}

      // If we have sub-groups specified, then print things out separately
      // by these

      if ( nk>1 )
	{

	  vector<CNVIndivReport> ak(nk);
	  vector<CNVIndivReport> uk(nk);

	  for (int k=0; k<nk; k++)
	    {

	      vector_t dummy;

	      summaryIndivSummaries( this,
				     k,
				     segmentCount,
				     segmentLength,
				     ak[k],uk[k],
				     dummy );

	      // Within group comparisons often likely to
	      // have zero counts

	      if ( ak[k].segCount == 0 )
		{
		  ak[k].t1 = ak[k].t2 = ak[k].t3 =
		    ak[k].t4 = ak[k].t5 = ak[k].t6 =
		    ak[k].t7 = ak[k].t8 = -9;
		}
	      if ( uk[k].segCount == 0 )
		{
		  uk[k].t1 = uk[k].t2 = uk[k].t3 =
		    uk[k].t4 = uk[k].t5 = uk[k].t6 =
		    uk[k].t7 = uk[k].t8 = -9;
		}

	    }

	  // Output in order of test
	  for (int k=0; k<nk; k++)
	    OUT << setw(12) << "N" << " "
		<< setw(8) << kname[k] << " "
		<< setw(12) << ak[k].segCount << " "
		<< setw(12) << uk[k].segCount << "\n";

	  for (int k=0; k<nk; k++)
	    OUT << setw(12) << "RATE" << " "
		<< setw(8) << kname[k] << " "
		<< setw(12) << ak[k].t1 << " "
		<< setw(12) << uk[k].t1 << "\n";

	  for (int k=0; k<nk; k++)
	    OUT << setw(12) << "PROP" << " "
		<< setw(8) << kname[k] << " "
		<< setw(12) << ak[k].t2 << " "
		<< setw(12) << uk[k].t2 << "\n";

	  for (int k=0; k<nk; k++)
	    OUT << setw(12) << "KBTOT" << " "
		<< setw(8) << kname[k] << " "
		<< setw(12) << ak[k].t3 << " "
		<< setw(12) << uk[k].t3 << "\n";

	  for (int k=0; k<nk; k++)
	    OUT << setw(12) << "KBAVG" << " "
		<< setw(8) << kname[k] << " "
		<< setw(12) << ak[k].t4 << " "
		<< setw(12) << uk[k].t4 << "\n";

	  if ( par::cnv_count )
	    {
	      for (int k=0; k<nk; k++)
		OUT << setw(12) << "GRATE" << " "
		    << setw(8) << kname[k] << " "
		    << setw(12) << ak[k].t5 << " "
		    << setw(12) << uk[k].t5 << "\n";

	      for (int k=0; k<nk; k++)
		OUT << setw(12) << "GPROP" << " "
		    << setw(8) << kname[k] << " "
		    << setw(12) << ak[k].t6 << " "
		    << setw(12) << uk[k].t6 << "\n";

	      for (int k=0; k<nk; k++)
		{
		  OUT << setw(12) << "GRICH" << " "
		      << setw(8) << kname[k] << " "
		      << setw(12) << ak[k].t7 << " "
		      << setw(12) << uk[k].t7 << "\n";
		}

	      if ( par::cnv_count_baseline )
		for (int k=0; k<nk; k++)
		  {
		    OUT << setw(12) << "GRICH2" << " "
			<< setw(8) << kname[k] << " "
			<< setw(12) << ak[k].t8 << " "
			<< setw(12) << uk[k].t8 << "\n";
		  }
	    }
	}


      OUT.close();

    }


  //////////////////////
  // Begin permutations

  bool finished = false;
  while(!finished)
    {


      ///////////////////////////////
      // Permute

      vector<double> pr(nt);

      if ( par::cnv_pos_perm )
	{
	  // Offset randomisation method, i.e. shift
	  // gene list by a random constant
	  positionPermuteSegments();

	  // Recount for each individual
	  indivSegmentSummaryCalc(segmentCount, segmentLength, true, true);

// 	  for (int j=0; j<20; j++)
// 	    cout << j <<  " is " << segment[j].count << "\n";
// 	  cout << "----\n";
	}
      else
	{
	  // Label swapping for phenotype
	  perm.permuteInCluster();
	}


      //////////////////////////////
      // Retest permuted dataset

      if ( ! par::cnv_indiv_perm )
	{

	  vector<int> coverage_aff(nt,0);
	  vector<int> coverage_unaff(nt,0);


	  ///////////////////////////
	  // Re-calculate counts

	  if ( par::seg_test_region )
	    {
	      countCNVPerRegion(coverage_aff,
				coverage_unaff);

	    }
	  else
	    {

	      /////////////////////////
	      // Actual segments

	      vector<Segment>::iterator s = segment.begin();
	      while ( s != segment.end() )
		{

		  if ( s->p1->pperson->aff)
		    for (int l = s->start ; l <= s->finish; l++) coverage_aff[l]++;
		  else
		    for (int l = s->start ; l <= s->finish; l++) coverage_unaff[l]++;

		  s++;
		}

	      // Optionally, if we allow 'wings' to increase span
	      // of events (and so, each data point represents the
	      // number of events with X kb of that position)

	      if ( par::seg_test_window )
		{

		  vector<Segment>::iterator s = segment.begin();
		  while ( s != segment.end() )
		    {

		      // Shift left from start

		      int l = s->start;
		      Locus * loc1 = locus[s->start];
		      while ( 1 )
			{
			  --l;
			  if ( l < 0 ) break;
			  Locus * loc2 = locus[l];
			  if ( loc2->chr != loc1->chr ) break;
			  if ( loc1->bp - loc2->bp > par::seg_test_window_bp  ) break;
			  if ( s->p1->pperson->aff )
			    ++coverage_aff[l];
			  else
			    ++coverage_unaff[l];
			}

		      // Shift right from start
		      l = s->finish;
		      loc1 = locus[s->finish];
		      while ( 1 )
			{
			  ++l;
			  if ( l == nl_all ) break;
			  Locus * loc2 = locus[l];
			  if ( loc2->chr != loc1->chr ) break;
			  if ( loc2->bp - loc1->bp > par::seg_test_window_bp  ) break;
			  if ( s->p1->pperson->aff )
			    ++coverage_aff[l];
			  else
			    ++coverage_unaff[l];
			}

		      // Next segment
		      s++;
		    }
		}

	    }

	  //////////////////////////
	  // Get genome-wide average

// 	  double Ag = 0;
// 	  double Ug = 0;
// 	  for (int l=0; l<nt; l++)
// 	    {
// 	      Ag += coverage_aff[l];
// 	      Ug += coverage_unaff[l];
// 	    }

// 	  Ag /= (double)nt;
// 	  Ug /= (double)nt;

	  for (int l=0; l<nt; l++)
	    {

	      double A = coverage_aff[l];
	      double U = coverage_unaff[l];

	      // 	    // Correct for genome-wide average
	      // 	    A = (A-Ag < 0) ? 0 : A - Ag;
	      // 	    U = (U-Ug < 0) ? 0 : U - Ug;
	      // 	    double statistic = A / (double)tot_aff - (double)U / tot_not ;
	      // 	    pr[l] = statistic < 0 ? 0 : statistic;


	      // Simple positional chi-sq
	      double A0 = tot_aff - A;
	      double U0 = tot_not - U;
	      double N = tot_aff + tot_not;

	      double E_A = ( tot_aff * (A+U) ) / N;
	      double E_A0 = ( tot_aff * (A0+U0) ) / N;
	      double E_U = ( tot_not * (A+U) ) / N;
	      double E_U0 = ( tot_not * (A0+U0) ) / N;

	      if ( E_A == 0 && E_U == 0 )
		pr[l] = 0;
	      else
		{
		  if ( par::segment_test_1sided && A/tot_aff < U/tot_not )
		    {
		      pr[l] = 0;
		    }
		  else
		    {
		      pr[l] =  ((A-E_A)*(A-E_A) ) / E_A +
			((A0-E_A0)*(A0-E_A0))/E_A0 +
			((U-E_U)*(U-E_U))/E_U +
			((U0-E_U0)*(U0-E_U0))/E_U0 ;

		      if ( ! realnum(pr[l]) )
			pr[l] = 0;
		    }
		}
	    }

	}


      ///////////////////////////
      // Individual summary tests

      if ( par::cnv_indiv_perm )
	{
	  // Only perform for all individuals
	  CNVIndivReport a, b;
	  summaryIndivSummaries( this, -1, segmentCount, segmentLength, a,b, pr );
	}



      ////////////////////////////////
      // Standard permutation counting

      finished = perm.update(pr,original);

    }

  if (!par::silent)
    cout << "\n\n";


  ////////////////////////////
  // Display permuted p-values

  f += ".mperm";
  printLOG("Writing permuted results for segment test to [ "+f+" ]\n");


  ofstream PHOM;
  PHOM.open( f.c_str() , ios::out );

  PHOM << setw(4)  << "CHR" << " ";
  if ( par::seg_test_region )
    PHOM << setw(16) << "REGION" << " ";
  else
    PHOM << setw(par::pp_maxsnp) << "SNP" << " ";
  PHOM  << setw(12) << "EMP1" << " ";
  if ( !par::cnv_indiv_perm )
    PHOM << setw(12) << "EMP2" << " ";
  PHOM << "\n";

  vector<string> testname(8);
  testname[0] = "RATE";
  testname[1] = "PROP";
  testname[2] = "KBTOT";
  testname[3] = "KBAVG";
  testname[4] = "GRATE";
  testname[5] = "GPROP";
  testname[6] = "GRICH";
  testname[7] = "GRICH2";

  if ( par::seg_test_region )
    {
      set<Range>::iterator i1 = geneList.begin();
      int gCount = 0;
      while ( i1 != geneList.end() )
	{
	    PHOM << setw(4)  << i1->chr << " "
		 << setw(16) << i1->name << " ";
	    PHOM << setw(12) << perm.pvalue( gCount ) << " ";
	    PHOM << setw(12) << perm.max_pvalue( gCount ) << "\n";

	  ++i1;
	  ++gCount;
	}
    }
  else
    for (int l=0; l<nt; l++)
      {
	if ( par::cnv_indiv_perm )
	  {
	    PHOM << setw(4) << "S" << " "
		 << setw(par::pp_maxsnp) << testname[l] << " ";
	    PHOM << setw(12) << perm.pvalue(l) << "\n";
	  }
	else
	  {
	    PHOM << setw(4)  << locus[l]->chr << " "
		 << setw(par::pp_maxsnp) << locus[l]->name << " ";
	    PHOM << setw(12) << perm.pvalue(l) << " ";
	    PHOM << setw(12) << perm.max_pvalue(l) << "\n";
	  }

      }

  PHOM.close();

}
