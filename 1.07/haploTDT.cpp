

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

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"
#include "stats.h"


//////////////////////////////
// Unweighted tests

void HaploPhase::haplotypicTDT(map<int,int> & tests, int nt, bool display)
{
  
  // No implementation of TDT omnibus test yet
  
  if ( nt != 2 ) 
    {
      result = -9;
      pvalue = -9;
      odds = -9;
      return;
    }

  // This might be a haplotype-specific test (i.e. of a single
  // haplotype) or of a group of haplotypes versus the rest

  // When rescoring the T and U counts, we will have supplied a 'downcoding' 
  // map, which is the same as the 'tests' map, i.e. mapping each haplotype
  // onto a 0/1 space (i.e. nt==2)
 
  
  // Find test haplotype(s), if we are in display mode (i.e. here we 
  // know it is not a group of haplotypes, but a specific haplotype, and
  // we want the name of it, --hap-tdt;  for all other instances, we 
  // can not assume that we will be testing a specific haplotype; so 
  // we do not bother about the name (it might be a group).  We 
  // can assume a binary test though (nt==2), so always set hh to 0.
  
  int hh=0;
  double tr = 0;
  double un = 0;

  // This test is always of 1 haplotype/group versus all others --
  // i.e. the downcoding will have been performed previously, with the
  // transmissions being appropriately rescored before hand (and so
  // trans[]/untrans[] will already be downcoded.

  if ( display )
    {
      map<int,int>::iterator i1 = tests.begin();
      while ( i1 != tests.end() )
	{
	  if ( i1->second == 0 )
	    {	  
	      hh = i1->first;
	      break;
	    }
	  i1++;
	}
    }

  tr += trans[hh];
  un += untrans[hh];



  ///////////////////////////////////////////////////////
  // 'result' visible outside this class, 
  
  // Either use McNemar's chi-square (b-c)^2/(b+c) 
  // or normal approximation for test of transmission 
  // ratio equals 0.5 (with the empirical variance added
  // here)
  
  odds = tr/un;

//   if ( true || useEmpiricalVariance )
//     {

  result  = tr - un;
  result *= result;
  result /= tr + un;
  pvalue = chiprobP(result,1);
  case_freq = tr;
  control_freq = un;


//     }
//   else 
//     {
	  
//       // Calculate empirical variance of transmissions: the 
//       // relevant quantities will have been stored during the 
//       // transmission scoring routine

//       double transmissionProportion = tr / ( tr + un );
//       double transmissionCount = tr + un;
// //       double eHH = transmissionX2[hh] / ( transmissionTotal-1);
// //       double eH = transmissionX[hh] / ( transmissionTotal-1);
//       double eHH = transmissionX2[hh] / ( transmissionCount );
//       double eH = transmissionX[hh] / ( transmissionCount );
//       empiricalVariance = eHH - ( eH * eH ); 
//       empiricalVariance /= transmissionCount - 1 ;
//       double Z = ( transmissionProportion - 0.5 ) 
// 	/ ( sqrt( empiricalVariance * ( 1 / transmissionCount ) ) );
//       result = Z * Z;
//      pvalue = chiprobP( result , 1 );
//  }
  

  if ( display )
    HTEST << setw(10) << hname << " "
	  << setw(12) << haplotypeName(hh) << " "
	  << setw(10) << trans[hh] << " "
	  << setw(10) << untrans[hh] << " ";


  if ( display )
    {
      if ( realnum(result) )
	{
	  HTEST << setw(10) << result << " "
		<< setw(10) << pvalue << " ";
	}
      else
	{
	  HTEST << setw(10) << "NA" << " "
		<< setw(10) << "NA" << " ";
	}

      for (int snps=0; snps<ns-1; snps++)
	HTEST << P.locus[S[snps]]->name << "|";
      
      HTEST << P.locus[S[ns-1]]->name << "\n";
    }

  
  return;
}



void HaploPhase::haplotypicWeightedTDT()
{

  vector_t weights;
  
  for (int i=0; i<nh; i++)
    {
      map<string,double>::iterator 
	whap = new_pred_weighted_allele[current].find( haplotypeName(i) );
      
      if ( whap != new_pred_weighted_allele[current].end() )
	{
	  weights.push_back( whap->second );
	}
      else
	{
	  weights.push_back( 0 );
	}
    }

  
  double T = 0;
  double U = 0;
  
  for (int h=0; h<nh; h++)
    {
      T += trans[h] * weights[h];
      U += untrans[h] * weights[h];
    }
  
  double chisq = ( (T-U)*(T-U) )  / ( T+U ) ;
  
  HTEST << setw(10) << hname << " "
	<< setw(12) << new_map[current]->allele1 << " "
	<< setw(10) << T << " " 
	<< setw(10) << U << " ";
  
  if ( realnum(chisq) )
    {
      HTEST << setw(10) << chisq << " "
	    << setw(10) << chiprobP(chisq,1) << " ";
    }
  else
    {
      HTEST << setw(10) << "NA" << " "
	    << setw(10) << "NA" << " ";
    }
  
  for (int snps=0; snps<ns-1; snps++)
    HTEST << P.locus[S[snps]]->name << "|";
  
  HTEST << P.locus[S[ns-1]]->name << "\n";
  
  return;

}

