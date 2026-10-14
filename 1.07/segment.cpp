

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
#include <cmath>
#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"
#include "fisher.h"
#include "stats.h"

void Plink::findSegments(int i1, int i2, vector_t & p, ofstream & SEG)
{

  Individual * p1 = sample[i1];
  Individual * p2 = sample[i2];
  
  bool inseg = false;
  
  Segment s;
  s.p1 = p1;
  s.p2 = p2;
  int npos = p.size();
		
  // Marker if minimal segment output file is selected
  bool already_seen_pair = false;
  
  for (int l=0; l<npos; l++)
    {
      
      ///////////////////////////
      // Start of a new segment?
      
      if ( (!inseg) && p[l] >= par::segment_threshold_start)
	{
	  inseg = true;

	  if (m1[l] != -1)
	    s.start = m1[l];
	  else
	    s.start = par::run_start;
	}
      
      
        ///////////////////////////
	// End of existing segment?
		    		    
	if ( inseg && 
	     ( p[l] < par::segment_threshold_finish || l == npos-1 ) )
	  {
	    inseg = false;

	    if (m2[l-1] != -1)
	      s.finish = m2[l-1];
	    else
	      s.finish = par::run_end;

	    // Do we like this segment?
	    if ( locus[s.finish]->bp - locus[s.start]->bp 
		 >= par::segment_length  
		 &&
		 s.finish-s.start+1 >= par::segment_snp  
		 )
	      {
		
		// Add segment to list
		segment.push_back(s);
		
		// Display?
		if (par::segment_output)
		  {
		    if ( par::segment_minimal )
		      {			
			if ( ! already_seen_pair )
			  {
			    already_seen_pair = true;
			    SEG << i1 << " " << i2 << " ";
			  }
			
			SEG << s.start << " "
			    << s.finish << " ";
		      }
		    else
		      {
			
			// Long format
			
			SEG << setw(par::pp_maxfid) << s.p1->fid << " "
			    << setw(par::pp_maxiid) << s.p1->iid << " "
			    << setw(par::pp_maxfid) << s.p2->fid << " "
			    << setw(par::pp_maxiid) << s.p2->iid << " ";
			
			if (par::bt)
			  {
			    if ( (!p1->aff) && (!p2->aff) ) 
			      SEG << setw(4) << "-1" << " ";
			    else if ( p1->aff && p2->aff )
			      SEG << setw(4) << "1" << " ";
			    else if ((!p1->aff) && p2->aff)
			      SEG << setw(4) << "0" << " ";
			    else if (p1->aff && !p2->aff)
			      SEG << setw(4) << "0" << " ";
			    else
			      SEG << setw(4) << "NA" << " ";
			  }
			else
			  SEG << setw(4) << "NA" << " ";
			
			Locus * start = locus[s.start];
			Locus * finish = locus[s.finish];
			
			SEG << setw(4) << par::run_chr << " "
			    << setw(10) << start->bp << " "
			    << setw(10) << finish->bp << " "
			    << setw(par::pp_maxsnp) << start->name << " "
			    << setw(par::pp_maxsnp) << finish->name << " "
			    << setw(6) << s.finish-s.start+1 << " "
			    << setw(10) << (double)(finish->bp - start->bp)/(double)1000 
			    << "\n";
			
			SEG.flush();
			
		      }
		  }
		
		
	      }
	    
	  }
    }

  // If we have found segments for this pair, then 
   // record a end-of-pair code if in minimal output mode

  if ( par::segment_minimal && already_seen_pair )
    SEG << "-1 -1\n";

}


void Plink::segmentIndividualTest(Perm & perm)
{
  
  printLOG("Writing individual-based segment tests to [ " 
	   + par::output_file_name + ".segtest1 ]\n");
  
  int total_cases = 0, total_controls = 0;

  for (int i=0; i<n; i++)
    if ( sample[i]->aff ) 
      total_cases++;
    else
      total_controls++;

  
  map<Individual*,int> ip;
  for (int i=0; i<n; i++)
    ip.insert(make_pair( sample[i],i ));
  

  if ( par::permute ) 
    {
      perm.setTests(nl_all);
      perm.setPermClusters(*this);
      perm.originalOrder();
    }
  
  
  
  /////////////////////////////////////
  // Generate original test statistics
  
  vector_t original = perm_segmentIndividualTest(perm,
						 true,
						 total_cases,
						 total_controls,
						 ip);
  
  if ( ! par::permute ) 
    return;
  

  //////////////////////////
  // Enter permutation stage
  
  bool finished = false;
  
  while(!finished)
    {
      
      perm.permuteInCluster();
      
      vector_t pr = perm_segmentIndividualTest(perm,
					       false,
					       total_cases,
					       total_controls,
					       ip);
      
      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);
      
    }
  
  if (!par::silent)
    cout << "\n\n";


  printLOG("Writing individual-based segment permutation results to [ " 
	   + par::output_file_name + ".segtest1.mperm ]\n");
  
  ofstream SEGS;
  
  SEGS.open( (par::output_file_name+".segtest1.mperm").c_str() , ios::out );

  SEGS << setw(4)  << "CHR" << " "
       << setw(par::pp_maxsnp) << "SNP" << " "
       << setw(12) << "EMP1" << " "      
       << setw(12) << "EMP2" << "\n";      
  
  
  for ( int l=0; l<nl_all; l++)
    {
      SEGS << setw(4) << locus[l]->chr << " "
	   << setw(par::pp_maxsnp) << locus[l]->name << " "
	   << setw(12) << perm.pvalue(l) << " "
	   << setw(12) << perm.max_pvalue(l) << "\n";      
    }

  SEGS.close();
 
}

vector_t Plink::perm_segmentIndividualTest(Perm & perm, 
					   bool display,
					   int total_cases,
					   int total_controls,
					   map<Individual*,int> & ip)
{
  
  vector_t results;

  ofstream SEGS;
  
  if ( display ) 
    {
      SEGS.open( (par::output_file_name+".segtest1").c_str() , ios::out );
      SEGS.precision(4);

      SEGS << setw(4)  << "CHR" << " "
	   << setw(par::pp_maxsnp) << "SNP" << " "
	   << setw(5) << "TEST" << " " 
	   << setw(8) << "AFF" << " " 
	   << setw(8) << "UNAFF" << " " 
	   << setw(8) << "PHAT" << " " 
	   << setw(8) << "P0" << " "
	   << setw(8) << "Z" << " "
	   << setw(8) << "P" << "\n";      
    }

  // Consider each SNP position
  // Count : number of cases with at least one segment
  //         number of controls with at least one segment
  //  versus : equivalent counts of those w/out a segment

  // Calculate genome-wide means for cases and controls

  // Tests:

  // 1) Number of cases with at least one segment versus 
  //    number of contorls with at least one segment

  // 2) Number of case segments versus numbr of control segments
  
  // Options:
  //    Make 1-sided (i.e. more case sharing expected)
  //    Adjust for genome-wide sharing level
  //    Use Fisher's exact versus use standard chi-square

  //    Use permutation or not (permute individuals)
  //    CHECK: okay to have floating point values for fisher's test?


  ///////////////////////////////////////////////////////
  // Calculate mean levels of sharing in cases and controls

  long int ncases = 0;
  long int ncontrols = 0;

  for (int l=0; l<nl_all; l++)
    {
      
      set<Individual*> people;
      
      vector<Segment>::iterator s = segment.begin();
      while ( s != segment.end() )
	{          
	  if ( s->start <= l && s->finish >= l )
	    {
	      
	      if ( people.find( s->p1 ) == people.end() )
		{
		  people.insert( s->p1 );
		  if ( s->p1->pperson->aff )
		    ncases++;
		  else
		    ncontrols++;
		}
	      
	      // Only consider second individual for IBD sharing test
	      if ( ! par::homo_run )
		{
		  if ( people.find( s->p2 ) == people.end() )

		    {
		      people.insert( s->p2 );
		      if ( s->p2->pperson->aff )
			ncases++;
		      else
			ncontrols++;
		    }
		} 
	    }
	  
	  s++;
	  
	} // next segment
    }
  
  double mean_case_sharing = (double)ncases / (double)(nl_all);
  double mean_control_sharing = (double)ncontrols / (double)(nl_all);
  
  
  
  /////////////////////////////////////////
  // Consider this position
  
  for (int l=0; l<nl_all; l++)
    {
      
      int ncases = 0, ncontrols = 0;
      set<Individual*> people;
      

      ///////////////////////////////////////////////////////
      // Consider all segments that might span this position
      
      vector<Segment>::iterator s = segment.begin();
      while ( s != segment.end() )
	{          
	  if ( s->start <= l && s->finish >= l )
	    {
	      
	      if ( people.find( s->p1 ) == people.end() )
		{
		  people.insert( s->p1 );
		  if ( s->p1->pperson->aff )
		    ncases++;
		  else
		    ncontrols++;
		}

	      // Only consider second individual for IBD sharing test
	      if ( ! par::homo_run )
		{
		  if ( people.find( s->p2 ) == people.end() )
		    {
		      people.insert( s->p2 );
		      if ( s->p2->pperson->aff )
			ncases++;
		      else
			ncontrols++;
		    }
		} 
	    }

	  s++;

	} // next segment
      

      //////////////////////////////////////////////////
      // Adjust table by the genome-wide means for case 
      // sharing and control sharing

      // Scale adjustment to be proportional to total amount of
      // sharing at this particular locus
      
//       double expected_cases = (ncases+ncontrols) * 
// 	( mean_case_sharing/(mean_case_sharing+mean_control_sharing) );
      
//       double expected_controls = (ncases+ncontrols) * 
// 	( mean_control_sharing/(mean_case_sharing+mean_control_sharing) );
      
//       double pvalue = 1;
      
//       if ( ncases > expected_cases )
// 	{
// 	  double chi1 = ncases - expected_cases;
// 	  chi1 *= chi1;
// 	  chi1 /= expected_cases;
	  
// 	  double chi2 = ncontrols - expected_controls;
// 	  chi2 *= chi2;
// 	  chi2 /= expected_controls;

// 	  double chisq = chi1 + chi2;
// 	  pvalue = chiprobP(chisq,1);
// 	}

      // Scale 
      

      //////////////////////////////////////////////////
      // Fisher's exact test, as counts might be small, 
      // or standard chi-sq

//       matrix_t t;
//       sizeMatrix(t,2,2);
//       t[0][0] = ncases_adj;
//       t[0][1] = ncontrols_adj;
//       t[1][0] = total_cases - ncases_adj;
//       t[1][1] = total_controls - ncontrols_adj;

//       cout << t[0][0] << "\t" 
// 	   << t[0][1] << "\t" 
// 	   << t[1][0] << "\t" 
// 	   << t[1][1] << "\n";

      //////////////////////////////////////////////////
      // Either chi-square of Fisher's exact for a p-value
      
//       double pvalue = 1;
      
//       // A one sided t-test:
      
//       if ( ncases_adj > ncontrols_adj )  
// 	   && (double)ncases/(double)total_cases 
// 	   > (double)ncontrols/(double)total_controls ) 
//	{	  
// 	  pvalue = par::segment_test_fisher 
// 	    ? 
// 	    fisher(t) : 
// 	    chiprobP(chi2x2(t),1);

// 	  pvalue = chiprobP(chi2x2(t),1);

// 	}
      

 
      //////////////////////////////////////////////////////////////
      // 1-sided test of proportions with the normal approximation
      
      // z = (p_hat - p0) / sqrt( ( p0*(1-p0) ) / N ) 
      // where N is total number of individuals
      // N*p should be > 5, 
      
      
      double p0 = mean_case_sharing/(mean_case_sharing+mean_control_sharing);
      double phat = (double)ncases/double(ncases+ncontrols);
      double z = (phat - p0) / sqrt( ( p0*(1-p0) ) / (ncases+ncontrols) ); 
      double pvalue = z < 0 ? 1 : normdist(-z);

      //      cout << "STAT: " << p0 << "\t" << phat << "\t" << z << "\n";


      ///////////////////
      // Store results?
      
      if ( ! realnum(p0) )
	{
	  pvalue = 1;
	}
      
      results.push_back(1-pvalue);

      
      ////////////////////
      // Display results?
      
      if ( display )
	{
	  if ( ! realnum(p0) )
	    SEGS << setw(4)  << locus[l]->chr << " " 
		 << setw(par::pp_maxsnp) << locus[l]->name << " "
		 << setw(5) << "ALL" << " " 
		 << setw(8) << ncases << " " 
		 << setw(8) << ncontrols << " " 
		 << setw(8) << "NA" << " "	    
		 << setw(8) << p0 << " "	    
		 << setw(8) << 0 << " " 
		 << setw(8) << 1 << "\n";      
	  else
	      SEGS << setw(4)  << locus[l]->chr << " " 
	       << setw(par::pp_maxsnp) << locus[l]->name << " "
	       << setw(5) << "ALL" << " " 
	       << setw(8) << ncases << " " 
	       << setw(8) << ncontrols << " " 
	       << setw(8) << phat << " "	    
	       << setw(8) << p0 << " "	    
	       << setw(8) << z << " " 
	       << setw(8) << pvalue << "\n";      
	}
      


      ////////////////////////////////////////////////////////////
      // Consider allele-specific segments spanning this position 
      // (groups of)

      // NEED TO ADJUST TEST STATISTIC BELOW, IF ABOVE FIX WORKS

      if ( false && par::segment_test_specific_segs ) 
	{
	  
	  groupSegmentsSpanning(l);
	  
	  map<int,int> groupCount;
	  
	  for ( int i = 0; i < n ; i++ )
	    {
	      
	      if ( indivSegmentGroup[i] == -1 )
		continue;
	      
	      map<int,int>::iterator gi = groupCount.find( indivSegmentGroup[i] );
	      
	      if ( gi == groupCount.end() ) 
		groupCount.insert(make_pair( indivSegmentGroup[i] , 1 ) );
	      else
		gi->second++;	  
	    }
	  
	  // Consider all allelically-matching groups with more than a 
	  // set number of segments
	  
	  map<int,int>::iterator gi = groupCount.begin();
	  
	  while ( gi != groupCount.end() )
	    {
	      if ( gi->second < 10 ) 
		{
		  gi++;
		  continue;
		}
	      
	      int group = gi->first;
	      
	      int ncases = 0, ncontrols = 0;
	      set<Individual*> people;
	      
	      // Consider all segments that might span this position
	      
	      vector<Segment>::iterator s = segment.begin();
	      while ( s != segment.end() )
		{          
		  if ( s->start <= l && s->finish >= l )
		    {
		      
		      if ( people.find( s->p1 ) == people.end() )
			{
			  if ( indivSegmentGroup[ ip.find( s->p1 )->second ] == group ) 
			    {
			      people.insert( s->p1 );		      
			      if ( s->p1->pperson->aff )
				ncases++;
			      else
				ncontrols++;
			    }
			}
		      
		      if ( ! par::homo_run )
			{
			  if ( people.find( s->p2 ) == people.end() )
			    {
			      if ( indivSegmentGroup[ ip.find( s->p2 )->second ] == group ) 
				{
				  people.insert( s->p2 );
				  if ( s->p2->pperson->aff )
				    ncases++;
				  else
				    ncontrols++;
				}
			    }
			} 
		    }
		  
		  s++;
		  
		} // next segment
	      
	      
	      // Fisher's exact test, as counts might be small
	      
	      table_t t;
	      sizeTable(t,2,2);
	      t[0][0] = ncases;
	      t[0][1] = ncontrols;
	      t[1][0] = total_cases - ncases;
	      t[1][1] = total_controls - ncontrols;

	      double pvalue = 1;

	      if ( ncases > 0 
		   && (double)ncases/(double)total_cases 
		   > (double)ncontrols/(double)total_controls ) 
		{	  
		  pvalue = par::segment_test_fisher ? fisher(t) : chiprobP(chi2x2(t),1);
		}

	      if ( display ) 
		{
		  SEGS << setw(4)  << locus[l]->chr << " "
		       << setw(par::pp_maxsnp) << locus[l]->name << " "
		       << setw(5) << group << " "
		       << setw(12) << ncases << " "
		       << setw(12) << ncontrols << " "
		       << setw(12) << pvalue << "\n";
		}

	      gi++;
	    }
	  
	} //end of specific segment section
    } // Next SNP  


  if (display)
    SEGS.close();


  return results;
  
}

void Plink::segmentPermutationTest(Perm & perm, bool ibd, string f,
				   vector<int> & coverage_conc_aff,
				   vector<int> & coverage_disc,
				   vector<int> & coverage_conc_unaff )
{
  

  // IBD or IBS segment permutation test for case/control data
  // Calculate total number of concordant/discordant pairs
  
  int tot_conc_aff=0;
  int tot_not=0;
  
  for (int i1=0; i1<n-1; i1++)
    for (int i2=i1+1; i2<n; i2++)
	if ( sample[i1]->aff && sample[i2]->aff ) tot_conc_aff++;
	else
	{
	    if ( par::segment_test_ignore_discordant ) 
	      {
		if ( (!sample[i1]->aff) && (!sample[i2]->aff) )
		  tot_not++;
	      }
	    else
	      tot_not++;
	}
  
  printLOG(int2str(tot_conc_aff)+" concordant affected pairs out of "
	   +int2str(tot_not+tot_conc_aff)+" in total\n");
  
  if ( par::segment_test_ignore_discordant ) 
      printLOG("Comparing case/case to control/control pairs\n");
  else
      printLOG("Comparing case/case to non-case/case pairs\n");


  int nt = nl_all;  // IBS mode
  if (ibd) { nt = nl; } 
  
  perm.setTests(nt);
  perm.setPermClusters(*this);
  perm.originalOrder();

  vector_t original(nt);
  vector_t origSC(nt);
  vector_t origSD(nt);

  // Get genome-wide means

  double SCg = 0;
  double SDg = 0;
  for (int l=0; l<nt; l++)
    {
      SCg += coverage_conc_aff[l];
      SDg += coverage_conc_unaff[l];
      if ( ! par::segment_test_ignore_discordant )
	  SDg += coverage_disc[l];
    }
  
  SCg /= (double)nt;
  SDg /= (double)nt;

  if ( par::segment_test_ignore_discordant )
    printLOG("Genome-wide average sharing (#segs, AA : UU) is " + dbl2str(SCg) + " : " + dbl2str(SDg) + "\n");
  else
    printLOG("Genome-wide average sharing (#segs, AA : AU+UU ) is " + dbl2str(SCg) + " : " + dbl2str(SDg) + "\n");

  for (int l=0; l<nt; l++)
    {
      
      double SC = coverage_conc_aff[l];
      double SD = coverage_conc_unaff[l];
      if ( ! par::segment_test_ignore_discordant )
	  SD += coverage_disc[l];

      // Correct for genome-wide average
//       SC = SC-SCg < 0 ? 0 : SC - SCg;
//       SD = SD-SDg < 0 ? 0 : SD - SDg;

      SC -= SCg;
      SD -= SDg;

      double statistic = SC / (double)tot_conc_aff  -  SD / (double)tot_not ;
//       original[l] = statistic < 0 ? 0 : statistic;
//       original[l] = (SC+SD) >0 ? statistic / (SC+SD) : statistic;
      original[l] = statistic;
      origSC[l] = SC;
      origSD[l] = SD;
    }
  
  
    //////////////////////
    // Begin permutations
    
    // Edit 'fringe' status if in IBD mode
    
     if (ibd)
       {
 	vector<Segment>::iterator s = segment.begin();
 	while ( s != segment.end() )
 	  {
 	    
	    if ( s->start == -1 ) s->start = par::run_start;
 	    if ( s->finish == -1 ) s->finish = par::run_end;
 	    
	    // Edit to make positions relative to start of this chromosome
 	    s->start -= par::run_start;
 	    s->finish -= par::run_start;
	    
	    s++;
 	  }
       }



    bool finished = false;
    while(!finished)
      {
	
	// Store permuted results
	
	vector<double> pr(nt);
	
	// Permute
	
	perm.permuteInCluster();
	
	// Retest
	
	vector<int> coverage_conc_aff(nt,0);
	vector<int> coverage_conc_unaff(nt,0);
	vector<int> coverage_disc(nt,0);
	
	vector<Segment>::iterator s = segment.begin();
	
	while ( s != segment.end() )
	  {    
	    
	    if (s->p1->pperson->aff == s->p2->pperson->aff)
	      {
		if ( s->p1->pperson->aff) 
		  for (int l = s->start ; l <= s->finish; l++) coverage_conc_aff[l]++;
		else
		  for (int l = s->start ; l <= s->finish; l++) coverage_conc_unaff[l]++;
	      }
	    else
	      for (int l = s->start ; l <= s->finish; l++) coverage_disc[l]++;
	    
	    s++;
	  }
	
	double SCg = 0;
	double SDg = 0;
	for (int l=0; l<nt; l++)
	{
	    SCg += coverage_conc_aff[l];
	    SDg += coverage_conc_unaff[l];
	    if ( ! par::segment_test_ignore_discordant )
		SDg += coverage_disc[l];
	}
	
	SCg /= (double)nt;
	SDg /= (double)nt;

	for (int l=0; l<nt; l++)
	  {
	    
	    double SC = coverage_conc_aff[l];
	    double SD = coverage_conc_unaff[l];
	    if ( ! par::segment_test_ignore_discordant )
		SD += coverage_disc[l];

	    // Correct for genome-wide average
//  	    SC = (SC-SCg < 0) ? 0 : SC - SCg;
//  	    SD = (SD-SDg < 0) ? 0 : SD - SDg;
	    
	    SC -= SCg;
	    SD -= SDg;
    
	    double statistic = SC / (double)tot_conc_aff - (double)SD / tot_not ;
// 	    pr[l] = statistic < 0 ? 0 : statistic;
// 	    pr[l] = (SC+SD) > 0 ? statistic / (SC+SD) : statistic;
	    pr[l] = statistic;
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
    printLOG("Writing segment test permutation results to [ "+f+" ]\n");
    
    ofstream SIBS;
    SIBS.open( f.c_str() , ios::out );
    SIBS.precision(4);

    SIBS << setw(4) << "CHR" << " "
	 << setw(par::pp_maxsnp) << "SNP" << " "
	 << setw(10) << "STAT" << " "
	 << setw(10) << "CONA" << " "
	 << setw(10) << "OTHER" << " "
	 << setw(10) << "EMP1" << " " 
	 << setw(10) << "EMP2" << "\n";

    for (int l=0; l<nt; l++)
      {	
	
	if (!ibd)
	  SIBS << setw(4)  << locus[l]->chr << " " 
	       << setw(par::pp_maxsnp) << locus[l]->name << " ";
	  
	  SIBS << setw(10) << original[l]  << " " 
	       << setw(10) << origSC[l] << " " 
	       << setw(10) << origSD[l] << " " 
	       << setw(10) << perm.pvalue(l) << " "
	       << setw(10) << perm.max_pvalue(l) << "\n";      
	
    }

  SIBS.close();

}



void Plink::testGenomeIBDByCovariate(Perm & perm)
{
  
  // Calculate case/case, control/control and case/control
  // Comparisons:

  //   1) Case/case versus all others
  //   2) Case/control versus all others
  //   3) Control/control versus all others

  //   4) Case/case versus case/control
  //   5) Control/Control versus case/control
  //   6) Case/case versus control/control

  // Test based on yes/no dichotomy for % pairs above certain PIHAT threshold 
  
  int tot_conc_aff=0;
  int tot_conc_unaff=0;
  int tot_disc=0;
  
  for (int i1=0; i1<n-1; i1++)
    for (int i2=i1+1; i2<n; i2++)
      {
	if ( sample[i1]->aff )
	  {
	    if ( sample[i2]->aff ) 
	      tot_conc_aff++;
	    else tot_disc++;
	  }
	else
	  {
	    if ( sample[i2]->aff ) 
	      tot_disc++;
	    else tot_conc_unaff++;
	  }
      }

  printLOG(int2str(tot_conc_aff)+" concordant affected pairs\n");
  printLOG(int2str(tot_conc_unaff)+" concordant unaffected pairs\n");
  printLOG(int2str(tot_disc)+" discordant pairs\n");
  printLOG(int2str(tot_conc_aff+tot_conc_unaff+tot_disc)+" total pairs\n");
  
  int nt = 6;  // six tests (above)

  perm.setTests(nt);
  perm.setPermClusters(*this);
  perm.originalOrder();


  ///////////////////////////
  // Original test statistics

  vector<double> original(nt);

  
  // Iterate over all pairs
  
  double prop_11 = 0;
  double prop_01 = 0;
  double prop_00 = 0;

  int c=0;
  for (int i1=0; i1<n-1; i1++)
    for (int i2=i1+1; i2<n; i2++)
      {
	Individual * p1 = sample[i1];
	Individual * p2 = sample[i2];
	
	int2 pair;
	pair.p1 = i1;
	pair.p2 = i2;

	// 	cout << c << " of " << saved_IBDg.size() << " check : \n";
	// 	cout << saved_IBDg[c].z1/2+saved_IBDg[c].z2 << "\n";
	// 	c++;
	//	cout <<"\n";


	if ( related.find(pair) != related.end() )
	  {
	    if ( sample[i1]->aff )
	      {
		if ( sample[i2]->aff ) 
		  prop_11++;
		else prop_01++;
	      }
	    else
	      {
		if ( sample[i2]->aff ) 
		  prop_01++;
		else prop_00++;
	      }
	  }
	
      }
  
  //   1) Case/case versus all others
  //   2) Case/control versus all others
  //   3) Control/control versus all others

  //   4) Case/case versus case/control
  //   5) Control/Control versus case/control
  //   6) Case/case versus control/control

  
  original[0] = (double)prop_11/(double)tot_conc_aff -   double(prop_01+prop_00)/(double)(tot_disc+tot_conc_unaff) ;
  original[1] = (double)prop_01/(double)tot_disc -       double(prop_11+prop_00)/(double)(tot_conc_aff+tot_conc_unaff) ;
  original[2] = (double)prop_00/(double)tot_conc_unaff - double(prop_01+prop_11)/(double)(tot_disc+tot_conc_aff) ;
  
  original[3] = (double)prop_11/(double)tot_conc_aff -   double(prop_01)/(double)(tot_disc) ;
  original[4] = (double)prop_00/(double)tot_conc_unaff - double(prop_01)/(double)(tot_disc) ;
  original[5] = (double)prop_11/(double)tot_conc_aff -   double(prop_00)/(double)(tot_conc_unaff) ;


  ////////////////////
  // Begin permutations
  
  bool finished = false;
  while(!finished)
    {
      
      // Store permuted results
      
      vector<double> pr(nt,0);
      
      // Permute
      
      perm.permuteInCluster();
      
      // Retest
      
      double prop_11 = 0;
      double prop_01 = 0;
      double prop_00 = 0;
      
      for (int i1=0; i1<n-1; i1++)
	for (int i2=i1+1; i2<n; i2++)
	  {
	    Individual * p1 = sample[i1];
	    Individual * p2 = sample[i2];
	    
	    int2 pair;
	    pair.p1 = i1;
	    pair.p2 = i2;
	    
	    if ( related.find(pair) != related.end() )
	      {
		if ( sample[i1]->pperson->aff )
		  {
		    if ( sample[i2]->pperson->aff ) 
		      prop_11++;
		    else prop_01++;
		  }
		else
		  {
		    if ( sample[i2]->pperson->aff ) 
		      prop_01++;
		    else prop_00++;
		  }
	      }
	    
	  }
      
      //   1) Case/case versus all others
      //   2) Case/control versus all others
      //   3) Control/control versus all others
      
      //   4) Case/case versus case/control
      //   5) Control/Control versus case/control
      //   6) Case/case versus control/control

      pr[0] = (double)prop_11/(double)tot_conc_aff -   double(prop_01+prop_00)/(double)(tot_disc+tot_conc_unaff) ;
      pr[1] = (double)prop_01/(double)tot_disc -       double(prop_11+prop_00)/(double)(tot_conc_aff+tot_conc_unaff) ;
      pr[2] = (double)prop_00/(double)tot_conc_unaff - double(prop_01+prop_11)/(double)(tot_disc+tot_conc_aff) ;
      
      pr[3] = (double)prop_11/(double)tot_conc_aff -   double(prop_01)/(double)(tot_disc) ;
      pr[4] = (double)prop_00/(double)tot_conc_unaff - double(prop_01)/(double)(tot_disc) ;
      pr[5] = (double)prop_11/(double)tot_conc_aff -   double(prop_00)/(double)(tot_conc_unaff) ;
      	
      
      ////////////////////////////////
      // Standard permutation counting
      
      finished = perm.update(pr,original);
	
      }
    
    if (!par::silent)
      cout << "\n\n";
  
  
    ////////////////////////////
    // Display permuted p-values
    
    string f = par::output_file_name + ".genome.mperm";
    printLOG("Writing permuted results for genome-wide IBD test to [ "+f+" ]\n");
    
    ofstream SIBS;
    SIBS.open( f.c_str() , ios::out );
    
    for (int l=0; l<nt; l++)
      {	
	
	SIBS << setw(4) << l << " "  
	     << setw(12) << original[l]  << " " 
	     << setw(12) << perm.pvalue(l) << " "
	     << setw(12) << perm.max_pvalue(l) << "\n";      
	
    }

  SIBS.close();

}

void Plink::validateSegments()
{
  // Individual major
  if (par::SNP_major) SNP2Ind();

  // Consider all segments

  string f = par::output_file_name + ".sibs";

  printLOG("Writing segmental IBS information to [ " + f + " ]\n");

  ofstream SIBS;
  SIBS.open(f.c_str(), ios::out);

  SIBS << setw(10) << "IBS0" << " "
       << setw(10) << "IBS1" << " "
       << setw(10) << "IBS2" << " "
       << setw(10) << "MISS" << "\n";
  
  vector<Segment>::iterator s = segment.begin();

  while ( s != segment.end() )
    {
      
      int c2 = 0;
      int c1 = 0;
      int c0 = 0;
      int cx = 0;

      for (int l = s->start; 
	   l <= s->finish; 
	   l++)
	{
	  
	  bool a1 = s->p1->one[l];
	  bool a2 = s->p1->two[l];

	  bool b1 = s->p2->one[l];
	  bool b2 = s->p2->two[l];

	  bool ibs0 = false;
	  bool ibs1 = false;
	  bool miss = false;

	  if ( a1 == a2 &&
	       b1 == b2 &&
	       a1 != b1  ) ibs0 = true;
	  else if ( a1 && !(a2) ) miss = true;
	  else if ( b1 && !(b2) ) miss = true;
	  else if ( a1 != b1 ||
		    a2 != b2 ) ibs1 = true;

	  if (ibs0) c0++;
	  else if (ibs1) c1++;
	  else if (miss) cx++;
	  else c2++;
	}

      SIBS << setw(10) << c0 << " " 
	   << setw(10) << c1 << " " 
	   << setw(10) << c2 << " " 
	   << setw(10) << cx << "\n"; 
    
      // Next segment 
      s++;

    }

  SIBS.close();
  

}



//////////////////////////////////////////////
// General segmental permutation test routine

void Plink::summaryIBSsegments(Perm & perm)
{

  vector<int> coverage_conc_aff(nl_all,0);
  vector<int> coverage_conc_unaff(nl_all,0);
  vector<int> coverage_disc(nl_all,0);
  
  vector<Segment>::iterator s = segment.begin();
  
  while ( s != segment.end() )
    {    
      
      if (s->p1->aff == s->p2->aff)
	{
	  if ( s->p1->aff) 
	    for (int l = s->start ; l <= s->finish; l++) 
	      coverage_conc_aff[l]++;
	  else
	    for (int l = s->start ; l <= s->finish; l++) 
	      coverage_conc_unaff[l]++;
	}
      else
	for (int l = s->start ; l <= s->finish; l++) 
	  coverage_disc[l]++;
      
      s++;
    }

  
  // Optionally, if we allow 'wings' to increase span 
  // of events (and so, each data point represents the 
  // number of events with X kb of that position)
  
  if ( ( par::homo_run || par::cnv_list ) && par::seg_test_window )
    {

      printLOG("Summarising segments within a window of " 
	       + int2str((int)(par::seg_test_window_bp/1000))
	       + " kb\n");

      
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
	      if ( s->p1->aff )
		++coverage_conc_aff[l];
	      else
		++coverage_conc_unaff[l];
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
	      if ( s->p1->aff )
		++coverage_conc_aff[l];
	      else
		++coverage_conc_unaff[l];
	    }
	  
	  // Next segment
	  s++;
	}	  
    }
  
  
  
  string f = par::output_file_name;

  if ( par::segment_output ) 
    f += ".segment.summary";
  else if ( par::homo_run )
    f += ".hom.summary";
  else if ( par::cnv_list )
    f += ".cnv.summary";
  else if ( par::ibs_2only )
    f += ".ibs2.summary";
  else
    f += ".ibs.summary";


  


  /////////////////////////////////////
  // Treat CNVs as homozygous segments
  // i.e. just a per-individual and not
  //      a per-pair phenomenon

  if ( par::cnv_list )
    par::homo_run = true;


  ofstream SIBS;
  SIBS.open( f.c_str() , ios::out );
  


  //////////////////////////////
  // And display

  if ( par::homo_run )
    SIBS << setw(4) << "CHR" << " " 
	 << setw(par::pp_maxsnp) << "SNP" << " " 
	 << setw(12) << "BP" << " "
	 << setw(8) << "AFF" << " " 
	 << setw(8) << "UNAFF" << "\n"; 
  else
    SIBS << setw(4) << "CHR" << " " 
	 << setw(par::pp_maxsnp) << "SNP" << " " 
	 << setw(12) << "BP" << " "
	 << setw(8) << "CONA" << " "
	 << setw(8) << "DISC" << " " 
	 << setw(8) << "CONU" << "\n";
  
  for (int l=0; l<nl_all; l++)
    {
      
      if ( par::homo_run )
	SIBS << setw(4) << locus[l]->chr << " " 
	     << setw(par::pp_maxsnp) << locus[l]->name << " " 
	     << setw(12) << locus[l]->bp << " "
	     << setw(8) << coverage_conc_aff[l] << " " 
	     << setw(8) << coverage_conc_unaff[l] << "\n"; 
      else
	SIBS << setw(4) << locus[l]->chr << " " 
	     << setw(par::pp_maxsnp) << locus[l]->name << " " 
	     << setw(12) << locus[l]->bp << " "
	     << setw(8) << coverage_conc_aff[l] << " "
	     << setw(8) << coverage_disc[l] << " " 
	     << setw(8) << coverage_conc_unaff[l] << "\n";       
    }
  
  SIBS.close();


  //////////////////////////////////
  // Permutation test, if C/C given?

  if (!par::permute) return;
  
  if (par::homo_run)
    homozygousSegmentPermutationTest(perm,f,coverage_conc_aff,coverage_conc_unaff);
  else
    segmentPermutationTest(perm,false,f,coverage_conc_aff,coverage_disc,coverage_conc_unaff);

  return;

}





void Plink::summaryIBDsegments(Perm & perm)
{

  // Number of map positions
  // (extra position for global IBD)
  int npos = m1.size() - 1;
  
  vector<int> coverage_conc_aff(nl,0);
  vector<int> coverage_conc_unaff(nl,0);
  vector<int> coverage_disc(nl,0);
  
  // pihat[pair][position]
  // check if, e.g. pihat>0.2 (i.e. atleast 80% chance of IBD1 sharing)
  // level set in options.cpp

      
  for (int l=0;l<npos;l++)
    {
      
      double p1, p2;
    
      if (m1[l]==-1) 
	{
	  p1 = locus[par::run_start]->pos - par::fringe;
	}
      else 
	{
	  p1 = locus[m1[l]]->pos;
	}
      
      if (m2[l]==-1) 
	{
	  p2 = locus[par::run_end]->pos + par::fringe;
	}
      else 
	{
	  p2 = locus[m2[l]]->pos;
	}
      
      if (m1[l]==-1 && m2[l]==-1)
	{
	  p1 = p2 = 0;
	}
      
      double d1 = p1 + pos[l] * (p2-p1);

      
      vector<Segment>::iterator s = segment.begin();
      
      while ( s != segment.end() )
	{    
	  
	  if ( locus[s->start]->chr != par::run_chr )
	    {
	      s++;
	      continue;
	    }
	  
	  
	  // In segment?
	  if (locus[s->start]->pos < d1 && locus[s->finish]->pos > d1  )
	    {
	      
	      if (s->p1->aff == s->p2->aff)
		{
		  if ( s->p1->aff) 
		    coverage_conc_aff[l]++;
		  else
		    coverage_conc_unaff[l]++;
		}
	      else
		coverage_disc[l]++;
	      
	    }
	  
	  s++;
	}
      
    }
  
  
  string f = par::output_file_name + ".segment.summary";

  ofstream SSEG;
  if (par::segment_output_started)
     SSEG.open( f.c_str() , ios::app );
  else
     SSEG.open( f.c_str() , ios::out );
  par::segment_output_started = true;

  //////////////////////////////
  // And display
  
  for (int l=0; l<npos; l++)
    {
      
      double p1, p2;
      string n1, n2;
      double fr=0;
      
      ////////////////////////////////
      // Single and multipoint output
      
      if (m1[l]==-1) 
	{
	  p1 = locus[par::run_start]->pos - par::fringe;
	  n1 = "fringe";
	}
      else 
	{
	  p1 = locus[m1[l]]->pos;
	  n1 = locus[m1[l]]->name;
	  fr = locus[m1[l]]->freq;
	}
      
      if (m2[l]==-1)
	{
	  p2 = locus[par::run_end]->pos + par::fringe;
	  n2 = "fringe";
	}
      else 
	{
	  p2 = locus[m2[l]]->pos;
	  n2 = locus[m2[l]]->name;
	}
      
      if (m1[l]==-1 && m2[l]==-1)
	{
	  p1 = p2 = 0;
	  n1 = "Genomewide";
	  n2 = "IBD";
	}

      double d1 = p1 + pos[l] * (p2-p1);

      SSEG << par::run_chr << " " 
	   << n1 << " " 
	   << n2 << " " 
	   << d1 << "  "
	   << fr << " " 
	   << coverage_conc_unaff[l] << " " 
	   << coverage_disc[l] << " " 
	   << coverage_conc_aff[l] << "\n"; 

    }

  SSEG.close();

  
  ////////////////////////
  // C/C permutation test?

  if (!par::permute) return;

  segmentPermutationTest(perm,true,f,coverage_conc_aff,coverage_disc,coverage_conc_unaff);

}



void Plink::summaryIBD()
{

  int npos = pihat[0].size();
  int npair = pihat.size();
  
  vector<int> coverage(npos,0);
  
  vector<int> coverage_conc_aff(npos,0);
  vector<int> coverage_conc_unaff(npos,0);
  vector<int> coverage_disc(npos,0);
  
  // pihat[pair][position]
  // check if, e.g. pihat>0.2 (i.e. atleast 80% chance of IBD1 sharing)
  // level set in options.cpp

  for (int i=0;i<pihat.size();i++)
    for (int j=0;j<npos;j++)
      {
	if (pihat[i][j] > par::IBD_threshold )
	  {
	    coverage[j]++;
	    
	    if (sample[pair1[i]]->phenotype 
		!= sample[pair2[i]]->phenotype) 
	      coverage_disc[j]++;
	    
	    if (sample[pair1[i]]->phenotype==1 
		&& sample[pair2[i]]->phenotype==1) 
	      coverage_conc_unaff[j]++;

	    if (sample[pair1[i]]->phenotype==2 
		&& sample[pair2[i]]->phenotype==2) 
	      coverage_conc_aff[j]++;
	  }
      }


  //////////////////////////////
  // And display
  
  for (int l=0; l<npos; l++)
    {
      
      double p1, p2;
      string n1, n2;
      double fr = 0;

      ////////////////////////////////
      // Single and multipoint output
      
      if (m1[l]==-1) 
	{
	  p1 = locus[par::run_start]->pos - par::fringe;
	  n1 = "fringe";
	}
      else 
	{
	  p1 = locus[m1[l]]->pos;
	  n1 = locus[m1[l]]->name;
	  fr = locus[m1[l]]->freq;
	}
      
      if (m2[l]==-1) 
	{
	  p2 = locus[par::run_end]->pos + par::fringe;
	  n2 = "fringe";
	}
      else 
	{
	  p2 = locus[m2[l]]->pos;
	  n2 = locus[m2[l]]->name;
	}
      
      if (m1[l]==-1 && m2[l]==-1)
	{
	  p1 = p2 = 0;
	  n1 = "Genomewide";
	  n2 = "IBD";
	}
      
      double d1 = p1 + pos[l] * (p2-p1);
      
      cout << "C " 
	   << par::run_chr << " " 
	   << n1 << " " 
	   << n2 << " " 
	   << d1 << "  "
	   << fr << " " 
	   << coverage[l] << " " 
	   << coverage_conc_unaff[l] << " " 
	   << coverage_disc[l] << " " 
	   << coverage_conc_aff[l] << " " 
	   << coverage[l]/(double)npair << "\n"; 
    }
}


void Plink::displayGMULTI(Individual * p1, Individual * p2, int l, ofstream & GMULTI)
{

  Locus * loc = locus[l];

  bool a1 = p1->one[l];
  bool a2 = p1->two[l];

  bool b1 = p2->one[l];
  bool b2 = p2->two[l];

  bool miss = false;
  
  GMULTI << setw(par::pp_maxfid) << p1->fid << " "
	 << setw(par::pp_maxiid) << p1->iid << " "
	 << setw(par::pp_maxfid) << p2->fid << " "
	 << setw(par::pp_maxiid) << p2->iid << " "
	 << setw(par::pp_maxsnp) << loc->name << " ";

  // Minor allele and frequency 

  GMULTI << setw(2) << loc->allele1 << " "
	 << setw(8) << loc->freq << " ";

  // Genotypes

  if ((!a1) && (!a2)) 
    GMULTI << setw(4) 
	   << loc->allele1+"/"+loc->allele1 
	   << " ";
  else if ((!a1) && a2) 
    GMULTI << setw(4) 
	   << loc->allele1+"/"+loc->allele2
	   << " ";
  else if (a1 && a2) 
    GMULTI << setw(4) 
	   << loc->allele2+"/"+loc->allele2
	   << " ";
  else 
    {
      GMULTI << setw(4) << "0/0" << " ";
      miss = true; 
    }
  
  if ((!b1) && (!b2)) 
    GMULTI << setw(4) 
	   << loc->allele1+"/"+loc->allele1;
  else if ((!b1) && b2) 
    GMULTI << setw(4) 
	   << loc->allele1+"/"+loc->allele2;
  else if (b1 && b2) 
    GMULTI << setw(4) 
	   << loc->allele2+"/"+loc->allele2;
  else 
    {
      GMULTI << setw(4) << "0/0";
      miss = true;
    }

  // IBS count

  if (miss)
    GMULTI << setw(3) << "NA" << " ";
  else
    {
      int ibs=0;
      if (a1==b1) ibs++;
      if (a2==b2) ibs++;
      GMULTI << setw(3) << ibs << " ";
    }

  GMULTI << "\n";
  
}


void Plink::indivSegmentSummaryCalc(map<indivPair,int> & segmentCount,
				    map<indivPair,double> & segmentLength,
				    bool countCases, bool countControls )
{

  segmentCount.clear();
  segmentLength.clear();
  segmentCount2.clear();

  if ( par::cnv_count_baseline )
    segmentCount2Baseline.clear();

  vector<Segment>::iterator s = segment.begin();
  
  while ( s != segment.end() )
    {    
      
      indivPair p;
      p.p1 = s->p1;
      p.p2 = s->p2;

      // For now, assume that case/control masking only applies
      // to homozygosity and CNV tests

      if ( ! countCases && s->p1->pperson->aff )
	continue;
      
      if ( ( !countControls )  && ! s->p1->pperson->aff )
	continue;


      // We have not yet seen this indiv/pair
      
      map<indivPair,int>::iterator ip = segmentCount.find(p);
      map<indivPair,double>::iterator il = segmentLength.find(p);
      map<indivPair,double>::iterator ic2 = segmentCount2.find(p);
      
      // KB length

      double l = (double)(locus[s->finish]->bp 
			  - locus[s->start]->bp)/(double)1000;
      
      if ( ip  == segmentCount.end() )
	{
	  segmentCount.insert( make_pair( p, 1 ) );
	  segmentLength.insert( make_pair( p, l ) );
	  	  
	  if ( par::cnv_weighted_gene_test )
	    {
	      segmentCount2.insert( make_pair( p, s->weightedCount ));
	      if ( par::cnv_count_baseline )
		segmentCount2Baseline.insert( make_pair( p, s->weightedBaseline ));
	    }
	  else
	    {
	      segmentCount2.insert( make_pair( p, s->count ));
	      if ( par::cnv_count_baseline )
		segmentCount2Baseline.insert( make_pair( p, s->baseline ));
	    }
	}
      else
	{
	  (ip->second)++;
	  (il->second) += l;
	  
	  if ( par::cnv_weighted_gene_test )
	    {
	      (ic2->second) += s->weightedCount;
	      if ( par::cnv_count_baseline )
		{
		  map<indivPair,double>::iterator ic2b = segmentCount2Baseline.find(p);
		  (ic2b->second) += s->weightedBaseline;
		}
	    }
	  else
	    {

	      (ic2->second) += s->count;

	      if ( par::cnv_count_baseline )
		{
		  map<indivPair,double>::iterator ic2b = segmentCount2Baseline.find(p);
		  (ic2b->second) += s->baseline;
		}
	    }

	}

      s++;
    }

}


void Plink::indivSegmentSummary()
{
  
  // Again, this function works for both homozygous and shared
  // segments (individuals and pairs)
  
  // Now these are part of Plink class
  //map<indivPair,int> segmentCount;
  //map<indivPair,double> segmentLength;
  
  // Generate basic summary file for all people
  indivSegmentSummaryCalc(segmentCount, segmentLength, true, true);
  

  ////////////
  // Output

  string f = par::output_file_name;
  
  if ( par::segment_output ) 
    f += ".segment.indiv";
  else if ( par::cnv_list )
    f += ".cnv.indiv";
  else if ( par::homo_run )
    f += ".hom.indiv";

  ofstream HOM;
  HOM.open( f.c_str() , ios::out );
  
  if ( par::homo_run )
    HOM << setw(par::pp_maxfid) << "FID" << " " 
	<< setw(par::pp_maxiid) << "IID" << " " 
	<< setw(4) << "PHE" << " "
	<< setw(8) << "NSEG" << " " 
	<< setw(8) << "KB" << " "
	<< setw(8) << "KBAVG" << "\n"; 
  else if ( par::cnv_list )
    {    
      HOM << setw(par::pp_maxfid) << "FID" << " " 
	  << setw(par::pp_maxiid) << "IID" << " " 
	  << setw(4) << "PHE" << " "
	  << setw(8) << "NSEG" << " "
	  << setw(8) << "KB" << " "
	  << setw(8) << "KBAVG" << " ";
      if ( par::cnv_count ) 
	HOM << setw(8) << "COUNT" << " ";
      HOM << "\n";
    }
  else
    HOM << setw(par::pp_maxfid) << "FID1" << " " 
	<< setw(par::pp_maxiid) << "IID2" << " " 
	<< setw(par::pp_maxfid) << "FID1" << " " 
	<< setw(par::pp_maxiid) << "IID2" << " " 
	<< setw(4) << "PHE" << " "
	<< setw(8) << "NSEG" << " " 
	<< setw(8) << "KB" << " " 
	<< setw(8) << "KBAVG" << "\n"; 

  
  // Output all individuals in homozygous segment mode

  if ( par::homo_run || par::cnv_list ) 
    {
      
      for ( int i = 0; i < n; i++)
	{
	
	  indivPair t;
	  t.p1 = t.p2 = sample[i];

	  map<indivPair,int>::iterator ic = segmentCount.find(t);
	  map<indivPair,double>::iterator il = segmentLength.find(t);
	  map<indivPair,double>::iterator ic2 = segmentCount2.find(t);
	  
	  indivPair p = ic->first;
	  
	  HOM << setw(par::pp_maxfid) << sample[i]->fid << " " 
	      << setw(par::pp_maxiid) << sample[i]->iid << " "
	      << setw(4) << sample[i]->phenotype << " ";
	  
	  if ( ic != segmentCount.end() ) 
	    {
	      HOM << setw(8) << ic->second << " " 
		  << setw(8) << il->second << " "
		  << setw(8) << il->second / (double)ic->second << " ";
	      if ( par::cnv_count )
		HOM << setw(8) << ic2->second << " ";
	      HOM << "\n";
	    }
	  else
	    {
	      HOM << setw(8) << 0 << " " 
		  << setw(8) << 0 << " ";
	      if ( par::cnv_count )
		HOM << setw(8) << 0 << " ";
	      HOM << setw(8) << 0 << "\n";
	    }

	}  // Next individual
    }
  else
    {
      // For now, just output obsevred pairs in segment mode
      // (i.e. typically too many pairs)

      map<indivPair,int>::iterator ic = segmentCount.begin();
      map<indivPair,double>::iterator il = segmentLength.begin();
      
      while ( ic != segmentCount.end() )
	{
	  indivPair p = ic->first;
	  
	  HOM << setw(par::pp_maxfid) << p.p1->fid << " " 
	      << setw(par::pp_maxiid) << p.p1->iid << " ";
	  
	  if ( p.p1 != p.p2 )
	    {
	      int pcode = 0;
	      if ( p.p1->aff ) 
		{
		  if ( p.p2->aff ) 
		    pcode = 1;
		  else 
		    pcode = 0;
		}
	      else
		{
		  if ( p.p2->aff ) 
		    pcode = 0;
		  else 
		    pcode = -1;
		}	    
	      
	      HOM << setw(par::pp_maxfid) << p.p2->fid << " " 
		  << setw(par::pp_maxiid) << p.p2->iid << " "
		  << setw(4) << pcode << " "; 
	    }
	  else
	    HOM << setw(4) << p.p1->phenotype << " ";
	  
	  HOM << setw(8) << ic->second << " " 
	      << setw(8) << il->second << " "
	      << setw(8) << il->second / (double)ic->second << "\n";
	  
	  ++ic;
	  ++il;
	}
    }
  
  HOM.close();
    
}

class SegmentSizeCmp {
public:  
  bool operator() (const Segment & s1, const Segment & s2) const 
  {
    int len1 = s1.finish - s1.start;
    int len2 = s2.finish - s2.start;
    
    if ( len1 < len2 ) 
      return true;
    if ( len1 > len2 ) 
      return false;

    if ( s1.p1 < s2.p1 ) 
      return true;
    
    if ( s1.p1 > s2.p1 ) 
      return false;

    return ( s1.start < s2.start );
  }
};




void Plink::displaySegmentsLong()
{
  
   printLOG("Writing long-format segment list to [ " 
 	   + par::output_file_name + ".cnv.seglist ]\n");

   ofstream SEG;
   string f = par::output_file_name + ".cnv.seglist";
   SEG.open( f.c_str(), ios::out );

   // Determine list of which chromosomes we need to report 
   // on

   set<int> chr;
   for (int l=0; l<nl_all; l++)
     chr.insert( locus[l]->chr );
   set<int>::iterator ichr = chr.begin();

   // Consider each chromosome, one at a time

   while ( ichr != chr.end() )
     {
    
       SEG << "\nChromosome " << *ichr << "\n\n";
    
       // Sort list of events in decreasing size
       // for this chromosome
    
       map<Segment,int,SegmentSizeCmp> smap;
       vector<Segment>::iterator s = segment.begin();      
       while ( s != segment.end() )
 	{
 	  if ( locus[ s->start ]->chr == *ichr )
 	    {
 	      smap.insert(make_pair(*s,0));
 	    }
 	  ++s;
 	}
    


       ///////////////////////////////////////////////////////
       // How many segments to consider on this chromosome?

       int nseg = smap.size();

       // For each segment, determine "height"
    
       map<Segment,int,SegmentSizeCmp>::iterator i1 = smap.begin();     

       while ( i1 != smap.end() )
 	{	  

 	  // Look at all smaller than this one
 	  map<Segment,int,SegmentSizeCmp>::iterator i2 = i1;

 	  while ( i2 != smap.end() )
 	    {
 	      if ( i1 == i2 )
 		{
 		  ++i2;
 		  continue;
 		}

 	      const Segment * s1 = &(i1->first);
 	      const Segment * s2 = &(i2->first);

//     	      cout << s1->start << " - " << s1->finish 
//     		   << " to " << s2->start << " - " << s2->finish << "\n";
//     	      cout << "overlap test\n";


 // 	      if ( ( s2->finish >= s1->start && s2->start <= s1->finish ) ||
 // 		   ( s1->finish >= s2->start && s1->start <= s2->finish ) )

 	      if ( s2->finish >= s1->start && s2->start <= s1->finish )     
 		{
 		  // Place this small segment one above the larger one
		  (i2->second) = (i1->second) + 1;
 		}
 	      ++i2;
 	    }
 	  ++i1;
 	}
    
       // Take smallest segment;
       // if any larger on overlaps

       // for (int i=1; i<nseg; i++)
       //  for (int j=i-1; j>=0; j--)
       // Use overlap <- function(x,y,k)
       //{ 
       // k$BP2[y] >= k$BP1[x] && k$BP1[y] <= k$BP2[x] 
       //}
       // if ( overlap(i,j,k) )
       // ++h[i];
  
 
       // Map into per-row logic

       // Find max height

       i1 = smap.begin();
       int maxh = 0;
       while ( i1 != smap.end() )
 	{
 	  if ( i1->second > maxh )
 	    maxh = i1->second;
 	  ++i1;
 	}

       ++maxh;

 //       cout << "Segs = \n";
 //       i1 = smap.begin();
 //       while ( i1 != smap.end() )
 // 	{
 // 	  const Segment * s = &(i1->first);
 // 	  cout << locus[ s->start ]->bp << " to " << locus[ s->finish ]->bp << ", h = " << i1->second << "\n";
 // 	  ++i1;
 // 	}
 //       cout << "-------\n";



       vector<bool> switches(maxh,false);
    
       // Display
       int chrmin = 0;
       int chrmax = 0;
       int tmp = 0;
       for (int l=0; l<nl_all; l++)
 	if ( locus[l]->chr == *ichr )
 	  {
 	    chrmin = chrmax = l;
 	    tmp = l;
 	    break;
 	  }
    
       for (int l=tmp; l<nl_all; l++)	
 	{
 	  if ( locus[l]->chr > *ichr )
 	    break;
 	  chrmax = l;
 	}

       int interval = 10 * 1000; // 10kb steps
       int fringe = 50 * 1000; // 50kb fringe

       bool done = false;
       int l = chrmin;
       int p = locus[ chrmin ]->bp - fringe;
       if ( p < 0 ) p = 0;

       // Track p (and move l along with it)
       bool firstPosition = true;

       while ( ! done )
 	{

 	  bool atLocus = false;
 	  if ( p == locus[l]->bp )
 	    atLocus = true;
	  
 	  vector<bool> starts(maxh,false);
 	  vector<bool> stops(maxh,false);	  
 	  vector<const Segment*> thisSeg(maxh);

 	  if ( atLocus )
 	    {
 	      i1 = smap.begin();
 	      while ( i1 != smap.end() )
 		{
 		  const Segment * s1 = &(i1->first);
 		  if ( s1->start == l ) 
 		    {
 		      starts[ i1->second ] = true;
 		      switches[ i1->second ] = true;
 		      thisSeg[ i1->second ] = s1;
 		    }
 		  else if ( s1->finish == l ) 
 		    {
 		      stops[ i1->second ] = true;
 		      switches[ i1->second ] = false;
 		      thisSeg[ i1->second ] = s1;
 		    }
 		  ++i1;
 		}
 	    }
	  

 	  //////////////
 	  // Output row
 	  bool display = atLocus;
	 
 	  if ( ! atLocus )
 	    {
 	      for (int j=0; j<maxh; j++)
 		if ( switches[j] )
 		  {
 		    display = true;
 		    break;
 		  }
 	    }
	  
 	  if ( display )
 	    {
 	      if ( atLocus )
 		SEG	<< setw(par::pp_maxsnp) << locus[l]->name << "     ";
 	      else
 		SEG << setw(par::pp_maxsnp) << "("+int2str(p)+")" << "     ";
	      
 	      // Symbols

 	      //  Duplication      +
 	      //  Deletion         -
 	      //  Case             A
               //  Control          U 

	      if ( par::cnv_write_freq )
		{
		  for (int j=0; j<maxh; j++)
		    {
		      if ( starts[j] ) 
			{
			  const Segment * seg = thisSeg[j];
			  if ( seg->freq > 9 ) 			    
			    SEG << seg->freq;	
			  else
			    SEG << seg->freq << " ";	
			}
		      else if ( stops[j] )
			{
			  const Segment * seg = thisSeg[j];
			  if ( seg->freq > 9 ) 			    
			    SEG << seg->freq;	
			  else
			    SEG << seg->freq << " ";	
			}
		      else if ( switches[j] )
			{
			  //const Segment * seg = thisSeg[j];
			  // if ( seg->freq > 9 ) 			    
// 			    SEG << seg->freq;	
// 			  else
			  // SEG << seg->freq << " ";	
			  SEG << "| ";
			}
		      else
			SEG << "  ";
		    }
		  SEG << "\n";
		}
	      else
		{
		  for (int j=0; j<maxh; j++)
		    {
		      if ( starts[j] ) 
			{
			  const Segment * seg = thisSeg[j];
			  if ( seg->type == 1 ) 
			    SEG << "+";
			  else
			    SEG << "-";
			}
		      else if ( stops[j] )
			{
			  const Segment * seg = thisSeg[j];
			  if ( seg->p1->aff ) 
			    SEG << "A";
			  else
			    SEG << "U";
			}
		      else if ( switches[j] )
			SEG << "|";
		      else
			SEG << " ";
		    }
		  SEG << "\n";
		}
	    }

 	  // Advance pointer
 	  p += interval;
	  
 	  // Are we done?	  
 	  if ( p > locus[chrmax]->bp + fringe )
 	    break;
	  
 	  // note -- we might have missed the first one...

 	  // Advance only as far as next SNP position
 	  if ( firstPosition && locus[l]->bp <= p )
 	    {
 	      p = locus[l]->bp;
 	      firstPosition = false;
 	    }
 	  else if ( l < chrmax && locus[l+1]->bp <= p )
 	    {
 	      ++l;
 	      firstPosition = false;
 	      p = locus[l]->bp;
 	    }
	  
 	}

       ++ichr;
     } // Next chromosome


   SEG.close();

}     




void Plink::displaySegmentsBED()
{


  //  BED format provides a flexible way to define the data lines that
  //  are displayed in an annotation track. BED lines have three
  //  required fields and nine additional optional fields. The number
  //  of fields per line must be consistent throughout any single set
  //  of data in an annotation track. The order of the optional fields
  //  is binding: lower-numbered fields must always be populated if
  //  higher-numbered fields are used.

  //The first three required BED fields are:

  // chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random)
  // or scaffold (e.g. scaffold10671).

  // chromStart - The starting position of the feature in the
  // chromosome or scaffold. The first base in a chromosome is
  // numbered 0.

  // chromEnd - The ending position of the feature in the chromosome
  // or scaffold. The chromEnd base is not included in the display of
  // the feature. For example, the first 100 bases of a chromosome are
  // defined as chromStart=0, chromEnd=100, and span the bases
  // numbered 0-99.

  //      The 9 additional optional BED fields are:

  // name - Defines the name of the BED line. This label is displayed
  // to the left of the BED line in the Genome Browser window when the
  // track is open to full display mode or directly to the left of the item in pack mode.

  //  score - A score between 0 and 1000. If the track line useScore
  //  attribute is set to 1 for this annotation data set, the score
  //  value will determine the level of gray in which this feature is
  //  displayed (higher numbers = darker gray).

  // strand - Defines the strand - either '+' or '-'.

  // thickStart - The starting position at which the feature is drawn
  // thickly (for example, the start codon in gene displays).

  // thickEnd - The ending position at which the feature is drawn
  // thickly (for example, the stop codon in gene displays).

  // itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the
  // track line itemRgb attribute is set to "On", this RBG value will
  // determine the display color of the data contained in this BED
  // line. NOTE: It is recommended that a simple color scheme (eight
  // colors or less) be used with this attribute to avoid overwhelming
  // the color resources of the Genome Browser and your Internet
  // browser.

  // blockCount - The number of blocks (exons) in the BED line.

  // blockSizes - A comma-separated list of the block sizes. The
  // number of items in this list should correspond to blockCount.

  // blockStarts - A comma-separated list of block starts. All of the
  // blockStart positions should be calculated relative to
  // chromStart. The number of items in this list should correspond to
  // blockCount.
	
  // Here's an example of an annotation track that uses a complete BED definition:

  //  track name=pairedReads description="Clone Paired Reads" useScore=1
  //  chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
  //  chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
  
  
  printLOG("Writing CNV information as BED track to [ " 
	   + par::output_file_name 
	   + ".cnv.bed ]\n");

  ofstream MOUT;

  MOUT.open( ( par::output_file_name + ".cnv.bed").c_str() , ios::out );
  
  MOUT << "track name=delCases description=\"Deletions, cases (PLINK CNV track)\" visibility=4 priority=1 itemRgb=\"On\"\n"; 
  vector<Segment>::iterator s = segment.begin();
  while ( s != segment.end() )
    {
      
      Individual * person = s->p1;
      
      if ( s->type == 1 && person->aff ) 
	{

	  int pos1 = locus[s->start]->bp;
	  int pos2 = locus[s->finish]->bp+1;
            
	  MOUT << "chr" << chromosomeName( locus[s->start]->chr ) << " " 
	       << pos1 << " " << pos2 << " "
	       << (person->fid + "_" + person->iid ) << " " 
	       << s->score << " " 
	       << ". "
	       << pos1 << " " << pos2 << " ";
	  // Colour
	  
	  if ( par::cnv_col == 0 )
	    MOUT << "255,0,0\n";
	  else if ( par::cnv_col = 1 ) 
	    MOUT << "0,0,255\n";
	  else if ( par::cnv_col = 2 ) 
	    MOUT << "0,255,0\n";
	  else if ( par::cnv_col = 3 ) 
	    MOUT << "255,0,0\n";
	}

      ++s;
    }


  MOUT << "track name=dupCases description=\"Duplications, cases (PLINK CNV track)\" visibility=4 priority=3 itemRgb=\"On\"\n"; 

  s = segment.begin();
  while ( s != segment.end() )
    {
      
      Individual * person = s->p1;      
      if ( s->type == 2 && person->aff ) 
	{

	  int pos1 = locus[s->start]->bp;
	  int pos2 = locus[s->finish]->bp+1;
	  
	  MOUT << "chr" << chromosomeName( locus[s->start]->chr ) << " " 
	       << pos1 << " " << pos2 << " "
	       << (person->fid + "_" + person->iid ) << " " 
	       << s->score << " " 
	       << ". " << pos1 << " " << pos2 << " ";
	  
	  // Colour
	  if ( par::cnv_col == 0 )
	    MOUT << "0,0,255\n";
	  else if ( par::cnv_col = 1 ) 
	    MOUT << "0,0,255\n";
	  else if ( par::cnv_col = 2 ) 
	    MOUT << "0,255,0\n";
	  else if ( par::cnv_col = 3 ) 
	    MOUT << "255,0,0\n";
	} 
      ++s;
    }
  
  MOUT << "track name=delControls description=\"Deletions, controls (PLINK CNV track)\" visibility=4 priority=2 itemRgb=\"On\"\n"; 
  s = segment.begin();
  while ( s != segment.end() )
    {
      
      Individual * person = s->p1;      
      if ( s->type == 1 && ! person->aff ) 
	{

	  int pos1 = locus[s->start]->bp;
	  int pos2 = locus[s->finish]->bp+1;
	  
	  MOUT << "chr" << chromosomeName( locus[s->start]->chr ) << " " 
	       << pos1 << " " << pos2 << " "
	       << (person->fid + "_" + person->iid ) << " " 
	       << s->score << " " 
	       << ". " << pos1 << " " << pos2 << " ";
	  
	  // Colour
	  if ( par::cnv_col == 0 )
	    MOUT << "128,0,0\n";
	  else if ( par::cnv_col = 1 ) 
	    MOUT << "0,0,255\n";
	  else if ( par::cnv_col = 2 ) 
	    MOUT << "0,255,0\n";
	  else if ( par::cnv_col = 3 ) 
	    MOUT << "255,0,0\n";
	}
      ++s;
    }
  
  MOUT << "track name=dupControls description=\"Duplications, controls (PLINK CNV track)\" visibility=4 priority=4 itemRgb=\"On\"\n"; 
  s = segment.begin();
  while ( s != segment.end() )
    {
      
      Individual * person = s->p1;      
      if ( s->type == 2  && ! person->aff ) 
	{

	  int pos1 = locus[s->start]->bp;
	  int pos2 = locus[s->finish]->bp+1;
	  
	  MOUT << "chr" << chromosomeName( locus[s->start]->chr ) << " " 
	       << pos1 << " " << pos2 << " "
	       << (person->fid + "_" + person->iid ) << " " 
	       << s->score << " " 
	       << ". " << pos1 << " " << pos2 << " ";
	  
	  // Colour
	  if ( par::cnv_col == 0 )
	    MOUT << "0,0,128\n";
	  else if ( par::cnv_col = 1 ) 
	    MOUT << "0,0,255\n";
	  else if ( par::cnv_col = 2 ) 
	    MOUT << "0,255,0\n";
	  else if ( par::cnv_col = 3 ) 
	    MOUT << "255,0,0\n";
	}
      ++s;
    }
  
  MOUT.close();
  
}
