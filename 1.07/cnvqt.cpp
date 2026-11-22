

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
#include <vector>
#include <map>
#include <iterator>
#include "plink.h"
#include "helper.h"
#include "options.h"
#include "perm.h"


extern Plink * PP;


void Plink::runTestCNVwithQT(Perm & perm)
{

  // Permutation test for mean difference in QT between people with
  // versus without a CNV.  By default two-sided, unless
  // par::segment_test_force_1sided = T
  
  // Optionally allowed for this to operate on smoothed data (i.e. 
  // average of event count over a KB window, forwards and backwards
  // from the given position)

  // Also performs genome-wide burden analyses for QTs -- is there 
  // an association between CNV size and QT, for example, etc.  These
  // are based on standard correlation
  
  int validN = 0;
  double grandMean = 0;
  
  for (int i=0; i<n; i++)
    {
      if ( !sample[i]->missing ) 
	{
	  grandMean += sample[i]->phenotype;
	  ++validN;
	}
    }
  
  grandMean /= (double)validN;
 
  printLOG("Total sample mean is " + dbl2str(grandMean) + ", based on " 
	   + int2str( validN ) + " individuals\n");


  //////////////////////////////////////////
  // Test positons = MAP positions (nl_all)
  // Test positions = summed segment counts ( get from original counts )
  // Test position = aggregate statistics ( 7 tests)

  int nt = nl_all;
  

  // IGNORE THIS FOR NOW...
  //   if ( par::seg_test_region ) 
  //     nt = coverage_aff.size();
  



  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Set up for individual burden tests?                            //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

//   if ( par::cnv_indiv_perm )
//     nt = 7;
  
  // Option per-individual summary tests? (4 tests)
  // Correlation between QT and these measures:

  //  total # segs 
  //  # people w/ 1+ seg
  //  total kb length
  //  mean segment length

  //  gene-count
  //  atleast-1-gene-count
  //  gene-enrichment



  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Initialise permutation procedures                              //
  //                                                                //
  ////////////////////////////////////////////////////////////////////
  
  perm.setTests(nt);
  perm.setPermClusters(*this);
  perm.originalOrder();
  
  vector_t original(nt);
  


  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Standard positional tests                                      //
  //                                                                //
  ////////////////////////////////////////////////////////////////////
  
  
  if ( par::cnv_indiv_perm )
    error("Not implemented --cnv-indiv-perm for QTs yet");
      

  // Test statistic is difference in QT bewteen people with 
  // versus without a CNV at this position

  vector_t count;
  vector_t m1;
  vector_t m0;

  original = testCNVwithQT(grandMean, validN, nt, count, m1, m0);
  

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Report to summary file                                         //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

  string f = par::output_file_name + ".cnv.qt.summary";
  printLOG("Writing CNV QT summary to [ "+f+" ]\n");
  
  ofstream FOUT;
  FOUT.open( f.c_str() , ios::out );
  FOUT.precision(4);

  FOUT << setw(4) << "CHR" << " "
       << setw(par::pp_maxsnp) << "SNP" << " "
       << setw(12) << "BP" << " "
       << setw(8) << "NCNV" << " "
       << setw(12) << "M1" << " "
       << setw(12) << "M0" << "\n";
    
  for (int l=0; l<nt; l++)
    {
      FOUT << setw(4) << locus[l]->chr << " "
	   << setw(par::pp_maxsnp) << locus[l]->name << " "
	   << setw(12) << locus[l]->bp << " "
	   << setw(8) << count[l] << " ";
      if ( count[l] > 0 ) 
	FOUT << setw(12) << m1[l] << " ";
      else
	FOUT << setw(12) << "NA" << " ";
      
      FOUT << setw(12) << m0[l] << "\n";
      
    }
  
  FOUT.close();

  


  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Run permutations                                               //
  //                                                                //
  ////////////////////////////////////////////////////////////////////


  bool finished = false;

  while(!finished)
    {      
      perm.permuteInCluster();
      vector_t pr = testCNVwithQT(grandMean, validN, nt, count, m1, m0);	  
      finished = perm.update(pr,original);      
    }
  
  if (!par::silent)
    cout << "\n\n";
  
  

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Display permuted results                                       //
  //                                                                //
  ////////////////////////////////////////////////////////////////////
  
  f += ".mperm";
  printLOG("Writing CNV QT permutation results to [ "+f+" ]\n");
  
  FOUT.open( f.c_str() , ios::out );
  FOUT.precision(4);

  FOUT << setw(4) << "CHR" << " "
       << setw(par::pp_maxsnp) << "SNP" << " "
       << setw(12) << "BP" << " "
       << setw(12) << "EMP1" << " "
       << setw(12) << "EMP2" << "\n";
    
  for (int l=0; l<nt; l++)
    {
      FOUT << setw(4) << locus[l]->chr << " "
	   << setw(par::pp_maxsnp) << locus[l]->name << " "
	   << setw(12) << locus[l]->bp << " "
	   << setw(12) << perm.pvalue( l ) << " " 
	   << setw(12) << perm.max_pvalue( l ) << "\n";
           
    }

  FOUT.close(); 

}



vector_t Plink::testCNVwithQT( double grandMean, int validN, int nt , 
			       vector_t & count,
			       vector_t & m1, 
			       vector_t & m0 )
{


  vector_t score(nt,0);

  m1.clear();
  m0.clear();
  count.clear();
  
  m1.resize(nt,0);
  m0.resize(nt, grandMean * validN) ;
  count.resize(nt,0);
  

  
  // Calculate QT mean for people with CNVs

  vector<Segment>::iterator s = segment.begin();  
  while ( s != segment.end() )
    {    
      for (int l = s->start ; l <= s->finish; l++) 
	{
	  ++count[ l ];
	  m1[ l ] += s->p1->pperson->phenotype;
	}
      ++s;
    }
  
  
  // Calculate QT mean for all other people, given grand mean
  
  for ( int l = 0 ; l < nl_all ; l++ ) 
    {     
      
      int k = validN - (int)count[l] ;
      
      m0[l] = k > 0 ? ( m0[l] - m1[l] ) / (double)k : 0 ;
      
      m1[l] = count[ l ] > 0 ? 
	m1[l] / count[ l ] : 
	0 ;
      
      score[ l ] = m1[l] - m0[l];
      
      if ( par::segment_test_force_1sided )
	{
	  if ( score[ l ] < 0 ) 
	    score[ l ] = 0;
	}
            
    }

  return score;

}
 

