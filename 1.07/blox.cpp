

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
#include <string>
#include <cmath>
#include <cstdlib>

#include "options.h"
#include "helper.h"
#include "plink.h"
#include "phase.h"
#include "stats.h"

extern Plink * PP;



///////////////////////////////////////////////////////////////////////
//                                                                   //
// Haplotype block code, adapted from code courtesy of Jeff Barrett, //
// and HAPLOVIEW, following some very quick-and-dirty Java->C++...   //
//                                                                   //
/////////////////////////////////////////////////////////////////////// 



class LDPair
{

public:

  int s1;
  int s2;
  int dist;  

  LDPair(int s1_, int s2_, int dist_)
  { 
    s1 = s1_; s2 = s2_; dist = dist_; 
  }

  friend ostream & operator<<(ostream & out, LDPair & v)
  {
    out << "[" << v.s1 << " " << v.s2 << " " << v.dist << "]";    
    return out;
  }
  
};

struct Pair_cmp
{
  bool operator()(const LDPair & a, const LDPair & b) const
  {
    if ( a.dist < b.dist ) return true;
    if ( a.dist > b.dist ) return false;
    if ( a.s1 < b.s1 ) return true;
    if ( a.s1 > b.s1 ) return false;
    return ( a.s2 < b.s2 );
  }
};


class DPrime
{
public:
  double dp;
  double dpl;
  double dpu;
  double lod;
};

class PairwiseLinkage
{
public:
  PairwiseLinkage(int a_, int b_)
  {
    a=a_;
    b=b_;
    knownAA = knownAB = knownBA = knownBB = unknownDH = 0;
  }

  int a;
  int b;
  double dp, rsq;
  double dp_upper, dp_lower;
  double lod;
  void calculateCI();
  void calculateLD();
  int knownAA, knownAB, knownBA, knownBB, unknownDH;
  
};



map<Range,vector<int> > Plink::mkBlks(int null1, int null2 )
{
  

  // First SNP, vector of SNPs (inc. first)
  
  map< int, vector<int> > blocks;
  
  
  // Some constants
  
  const double cutHighCI = 0.98;
  const double cutLowCI = 0.70;
  const double cutLowCIVar [5] = {0,0,0.80,0.50,0.50};
  const double maxDist [5] = {0,0,20000,30000,1000000};
  const double recHighCI = 0.90;
  const double informFrac = 0.95;
  const double fourGameteCutoff = 0.01;
  const double mafThresh = 0.05;
  
  // Set to skip SNPs with low MAFs
  // Uses genome-wide reference number: need to allocate for all SNPs here
  
  vector<bool> skipMarker(nl_all,false);
  for (int x = 0; x < nl_all; x++)
    skipMarker[x] = locus[x]->freq < mafThresh;
  
  // Consider each chromosome one at a time; skip X for now
  
  int startChromosome = locus[ 0 ]->chr;
  int finalChromosome = locus[ nl_all - 1 ]->chr;
  
  for (int chr = startChromosome ; chr <= finalChromosome; chr++)
    {

      if ( scaffold.find(chr) == scaffold.end() )
	continue;

      int fromPosition = scaffold[chr].lstart;
      int toPosition = scaffold[chr].lstop;
      
      int nsnps = toPosition - fromPosition + 1;
      

      /////////////////////////////////////////////////////////////////////////
      // Make a list of marker pairs in "strong LD", sorted by distance apart 
      
      set<LDPair,Pair_cmp> strongPairs;
      map<int2,DPrime> dpStore;
      
      int numStrong = 0; 
      int numRec = 0; 
      int numInGroup = 0;
      
      // Each pair of markers
      
      for (int x = fromPosition; x < toPosition; x++)
	{

 	  if ( ! par::silent )
 	    {
	      std::cerr << "Chromosome " <<  locus[x]->chr
			<< ", position " << locus[x]->bp/1000000.0 
			<< "Mb                \r";
	    }
	  
	  for (int y = x+1; y <= toPosition; y++)
	    {

	      if ( locus[x]->chr != locus[y]->chr ) 
		continue;
	      
	      if ( ( locus[y]->bp - locus[x]->bp ) > par::disp_r_window_kb )
		{
		  continue;
		}

	      if ( locus[x]->freq == 0 || locus[y]->freq == 0 )
		continue;
	      
	      PairwiseLinkage thisPair(x,y);
	      thisPair.calculateLD();
	      thisPair.calculateCI();
	      
	      double lod = thisPair.lod;
	      double lowCI = thisPair.dp_lower;
	      double highCI = thisPair.dp_upper;
	      
	      int2 t(x,y);
	      DPrime d;
	      d.dp = thisPair.dp;
	      d.dpl = lowCI;
	      d.dpu = highCI;
	      d.lod = lod;
	      dpStore.insert( make_pair( t,d ) );
	      
	      // Is this pair in strong LD?
	      if (lod < -90) continue; //missing data
	      
	      if (highCI < cutHighCI || lowCI < cutLowCI) 
		continue; //must pass "strong LD" test
	
	      // Store this pair
	      LDPair p(x,y, abs( locus[x]->bp - locus[y]->bp ) );
	      
	      
	      strongPairs.insert( p );

	    }
	}

      
      // Now we have a list of SNPs in strong LD within this region
      // Now construct blocks based on this
      
      set<int> used;
      
      // #blocks:
      vector<vector<int> > blockArray;
      
      int cnt = 0;
      
      for ( set<LDPair>::reverse_iterator i = strongPairs.rbegin();
	    i != strongPairs.rend();
	    ++i )
	{

	  int numStrong = 0; 
	  int numRec = 0; 
	  int numInGroup = 0;
	  
	  vector<int> thisBlock;
	  
	  int first = i->s1; 
	  int last = i->s2;
	  long sep = i->dist;
	  

	  // See if this block overlaps with another:
	  
	  if ( used.find(first) != used.end() 
	       || used.find(last)  != used.end() ) 
	    {	      
	      continue;
	    }
	  
	  // Next, count the number of markers in the block.
	  // (nb. assume all SNPs belong)
	  
	  for (int x = first; x <=last ; x++)
	    {
	      if( !skipMarker[x] ) 
		numInGroup++;
	    }


	  // Skip it if it is too long in bases for it's size in markers
	  if (numInGroup < 4 && sep > maxDist[numInGroup]) 
	    {
	      continue;
	    }
	  
	  // Add first SNP
	  
	  thisBlock.push_back( first );

	  // Test block: requires 95% of informative markers to be "strong"

	  for (int y = first+1; y <= last; y++)
	    {
	      if (skipMarker[y]) 
		{
		  continue;
		}
	      
	      thisBlock.push_back(y);
	      
	      
	      //loop over columns in row y
	      
	      for (int x = first; x < y; x++)
		{
		  
		  if (skipMarker[x]) 
		    continue;
		  
		  double lod; 
		  double lowCI; 
		  double highCI;
		  
		  map<int2,DPrime>::iterator l = dpStore.find( int2(x,y) );
		  
		  if ( l == dpStore.end() ) 
		    {
		      // Recalculate
		      PairwiseLinkage thisPair(x,y);
		      thisPair.calculateLD();
		      thisPair.calculateCI();
		      
		      lod = thisPair.lod;
		      lowCI = thisPair.dp_lower;
		      highCI = thisPair.dp_upper;
		    }
		  else
		    {
		      // Get the right bits
		      
		      lod = l->second.lod;
		      lowCI = l->second.dpl;
		      highCI = l->second.dpu;
		    }
		  
		  
		  // Monomorphic marker error
		  if ( lod < -90)
		    continue;   
		  
		  
		  // Skip bad markers
		  if ( lod == 0 && lowCI == 0 && highCI == 0)
		    continue; 
		  
		  // For small blocks use different CI cutoffs
		  
		  if (numInGroup < 5)
		    {
		      if (lowCI > cutLowCIVar[numInGroup] && highCI >= cutHighCI) 
			numStrong++;
		    }
		  else
		    {
		      if (lowCI > cutLowCI &&  highCI >= cutHighCI) 
			numStrong++; //strong LD
		    }
		  
		  if (highCI < recHighCI) 
		    numRec++; //recombination
		  
		}
	    }
	  
	  
	  // Change the definition somewhat for small blocks
	  
	  if (numInGroup > 3)
	    {
	      if (numStrong + numRec < 6) 
		{
		  continue;
		}
	    }
	  else if (numInGroup > 2)
	    {
	      if (numStrong + numRec < 3) 
		{
		  continue;
		}
	    }
	  else
	    {
	      if (numStrong + numRec < 1) 
		{
		  continue;
		}
	    }
	  

	  // If this qualifies as a block, add to the block list, but in
	  // order by first marker number:
	  
	  if ( (double)numStrong/(double)(numStrong + numRec) > informFrac)
	    { 
	      blocks.insert( make_pair( first , thisBlock ));  
	      
	      // Track that these SNPs belong to a block
	      for (int u = first; u <= last; u++)
		used.insert(u);
	    }
	  
	  
	}
      
      
      // Next chromosome
    }


  if ( ! par::silent )
    cerr << "\n";

  map<int,vector<int> >::iterator j = blocks.begin();

  printLOG(int2str( blocks.size() ) 
	   + " blocks called, writing list to [ " 
	   + par::output_file_name + ".blocks ]\n");
  ofstream O1( (par::output_file_name+".blocks").c_str() , ios::out );
  
  printLOG("Writing extra block details to [ " + 
	   par::output_file_name + ".blocks.det ]\n");
  ofstream O2( (par::output_file_name+".blocks.det").c_str() , ios::out );

  O2 << setw(4) << "CHR" << " " 
     << setw(12) << "BP1" << " "
     << setw(12) << "BP2" << " "
     << setw(12) << "KB" << " "
     << setw(6) << "NSNPS" << " "
     << setw(4) << "SNPS" << "\n";

  while ( j != blocks.end() )
    {
      O1 << "*";
      vector<int> & b = j->second;
      for (int k=0; k<b.size(); k++)
	O1 << " " << PP->locus[b[k]]->name;
      O1 << "\n";
      
      O2 << setw(4) << PP->locus[b[0]]->chr << " " 
	 << setw(12) << PP->locus[b[0]]->bp << " "
	 << setw(12) << PP->locus[b[b.size()-1]]->bp << " "
	 << setw(12) << (double)(PP->locus[b[b.size()-1]]->bp - PP->locus[b[0]]->bp + 1)/1000.0 << " "
	 << setw(6) << b.size() << " ";
      for (int k=0; k<b.size(); k++)
	{
	  if ( k>0 )
	    O2 << "|" << PP->locus[b[k]]->name;
	  else
	    O2 << PP->locus[b[k]]->name;
	}
      O2 << "\n";

      ++j;

    }
  

  O1.close();
  O2.close();
  

  // List of blocks created here
  // (dummy; not used)

  map<Range,vector<int> > blocks0;
  return blocks0;
  
}


void PairwiseLinkage::calculateCI()
{
  
  // Get counts of observed, unambiguous haplotypes
  
  vector<vector<int> > t = two_locus_table(a,b);
  
  // Assume autosome
  
  knownAA = 2 * t[0][0] + t[0][1] + t[1][0];
  knownAB = 2 * t[0][2] + t[0][1] + t[1][2];
  knownBA = 2 * t[2][0] + t[1][0] + t[2][1];
  knownBB = 2 * t[2][2] + t[2][1] + t[1][2];
  unknownDH = t[1][1];

  int total_chroms = knownAA + knownAB + knownBA + knownBB + 2*unknownDH;

  
  // From Haploview code:
  // Likelihood surface
  
  vector_t lsurface(101);
  
  //   // Assumed
  //   // denom = of D'
  //   // 4 haplotype frequencies pA1, pA2, pB1, pB2
  
  const double LN10 = log(10.0);
  
  string sA1 = PP->locus[a]->allele1 + PP->locus[b]->allele1;
  string sA2 = PP->locus[a]->allele1 + PP->locus[b]->allele2;
  string sB1 = PP->locus[a]->allele2 + PP->locus[b]->allele1;
  string sB2 = PP->locus[a]->allele2 + PP->locus[b]->allele2;

  
  double pA1,pA2,pB1,pB2;
  for ( int i = 0 ; i < 4 ; i++ ) 
    {
      if ( PP->haplo->haplotypeName(i) == sA1 ) 
	pA1 = PP->haplo->f[i];
      else if ( PP->haplo->haplotypeName(i) == sA2 ) 
	pA2 = PP->haplo->f[i];
      else if ( PP->haplo->haplotypeName(i) == sB1 ) 
	pB1 = PP->haplo->f[i];
      else if ( PP->haplo->haplotypeName(i) == sB2 ) 
	pB2 = PP->haplo->f[i];
    }


  double pA = pA1 + pA2;
  double pB = 1 - pA;
  double p1 = pA1 + pB1;
  double p2 = 1 - p1;
  
  // Estimated haplotype counts 

  double D = pA1 - (pA*p1);
  
  if (D < 0) 
    {
      
      double tmp;
      
      /* flip matrix so we get the positive D' */
      /* flip AA with AB and BA with BB */
      
      tmp=pA1; pA1=pA2; pA2=tmp;
      tmp=pB2; pB2=pB1; pB1=tmp;
      
      /* flip frequency of second allele */
      
      tmp=p1; p1=p2; p2=tmp;
      
      /* flip known array for likelihood computation */

      int tmpi;
      tmpi=knownAA; knownAA=knownAB; knownAB=tmpi;
      tmpi=knownBB; knownBB=knownBA; knownBA=tmpi;
    }
  

  double dmax1 = pA * p2 ;
  double dmax2 = pB * p1 ;
  double denom = dmax1 < dmax2 ? dmax1 : dmax2;
  
  for (int i=0; i<=100; i++) 
    {
      double dpr = (double)i*0.01;
      
      double tmpAA = dpr*denom + pA*p1;
      double tmpAB = pA-tmpAA;
      double tmpBA = p1-tmpAA;
      double tmpBB = pB-tmpBA;
      
      if (i==100) 
       	{
	  /* one value will be 0 */
	  if (tmpAA < 1e-10) tmpAA=1e-10;
	  if (tmpAB < 1e-10) tmpAB=1e-10;
	  if (tmpBA < 1e-10) tmpBA=1e-10;
	  if (tmpBB < 1e-10) tmpBB=1e-10;
	}      

      lsurface[i] = ( knownAA * log( tmpAA ) +
		      knownAB * log( tmpAB ) + 
		      knownBA * log( tmpBA ) + 
		      knownBB * log( tmpBB ) + 
		      unknownDH * log( tmpAA*tmpBB + tmpAB*tmpBA)) / LN10;
      
    }
  
  double loglike1 = unknownDH * log( pA1*pB2 + pB1*pA2 );
  
  if ( pA1>0 ) loglike1 += knownAA * log( pA1 );
  if ( pA2>0 ) loglike1 += knownAB * log( pA2 );
  if ( pB1>0 ) loglike1 += knownBA * log( pB1 );
  if ( pB2>0 ) loglike1 += knownBB * log( pB2 );
  loglike1 /= LN10;
  
  double loglike0= (knownAA * log(pA*p1) + 
		    knownAB * log(pA*p2) + 
		    knownBA * log(pB*p1) + 
		    knownBB * log(pB*p2) + 
		    unknownDH * log(2*pA*pB*p1*p2))/LN10;
  
  lod = loglike1-loglike0;
  if ( lod < 0 ) lod = 0;
  
  
  double total_prob=0;
  double sum_prob=0;
  int high_i = 0;
  int low_i = 0;
  for (int i=0; i<=100; i++) 
    {
      lsurface[i] -= loglike1;
      lsurface[i] = pow(10.0,lsurface[i]);
      total_prob += lsurface[i];      
    }
  
  for (int i=0; i<=100; i++) 
    {
      sum_prob += lsurface[i];
      if (sum_prob > 0.05*total_prob &&
	  sum_prob-lsurface[i] < 0.05*total_prob) {
	low_i = i-1;
	break;
      }
    }

  sum_prob=0.0;
  for (int i=100; i>=0; i--) 
    {
      sum_prob += lsurface[i];
      if (sum_prob > 0.05*total_prob &&
	  sum_prob-lsurface[i] < 0.05*total_prob) {
	high_i = i+1;
	break;
      }
    }
  
  if (high_i > 100){ high_i = 100; }
  
  dp_lower = (double)low_i/100.0;
  dp_upper = (double)high_i/100.0;
  
  if ( par::verbose )
    {
      cout << PP->locus[ a ]->name << "   " 
	   << PP->locus[ b ]->name << " : ";
      cout << "Rsq= " << PP->haplo->rsq(a,b) << " : ";
      cout << "D' = " << dp << " CI = " << dp_lower 
	   << " to " << dp_upper << "; lod = " << lod << " " 
	   << loglike1 << " vs. " << loglike0 << "\n";
    }
}

void PairwiseLinkage::calculateLD()
{
  dp = PP->haplo->dprime( a, b );
}
