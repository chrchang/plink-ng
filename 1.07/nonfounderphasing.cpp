

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
#include <cassert>

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "genogroup.h"
#include "phase.h"
#include "haplowindow.h"

extern Plink * PP;

void displayFamTran(map<FamilyTransmissions,double> & pmap, int fi, HaploPhase * HP)
{

  cout << "FAMILY " << fi << " : " << PP->sample[fi]->fid << "\n";
  map<FamilyTransmissions,double>::iterator i = pmap.begin();
  cout << setw(12) << "PATERNAL" << " "
       << setw(12) << "MATERNAL" << " " 
       << " -> " 
       << setw(12) << "OFFSPRING" << " "
       << setw(8) << "PROB" << "\n";
  
  while ( i != pmap.end() )
    {
      const FamilyTransmissions * f = &(i->first);

      cout << setw(12) << (HP->haplotypeName( f->pt ) + "/" + HP->haplotypeName( f->pu ) ) << " " 
	   << setw(12) << (HP->haplotypeName( f->mt ) + "/" + HP->haplotypeName( f->mu ) ) << " " 
	   << " -> " 
	   << setw(12) << (HP->haplotypeName( f->pt ) + "/" + HP->haplotypeName( f->mt ) ) << " " 
	   << setw(8) << i->second << "\n";

      ++i;
    }
    
  cout << "\n";


}

void HaploPhase::validateNonfounder(int i, 
				    vector<bool> & s1, 
				    vector<bool> & s2)
{
  
  // Flipping allele-coding for homozygotes
  
  for (int s=0; s<ns; s++)
    {
      if (par::SNP_major)
	{
	  s1[s] = P.SNP[S[s]]->one[i];
	  s2[s] = P.SNP[S[s]]->two[i];
	}
      else
	{
	  s1[s] = P.sample[i]->one[S[s]];
	  s2[s] = P.sample[i]->two[S[s]];
	}
      
      if (s1[s] == s2[s])
	{
	  s1[s] = !s1[s];
	  s2[s] = !s2[s];
	}
      
    }
  

  //////////////////////////////////////////////////////////
  // Count amount of missing genotype data at this position
  
  int nm = 0;
  for (int s=0; s<ns; s++)
    if (s1[s] && !s2[s])
      nm++;
  

  // If any missing genotypes, this person counts 
  // as ambiguous
  
  if (nm>0)
    ambig[i] = true;
  
  // But if too much missing genotype data, then 
  // we should not even try to phase this individual
  // for this region; note -- females should always be
  // missing all genotypes for Y, so we don't need to 
  // worry about allowing for a special case here.

  
  
  if ( (double)nm/(double)ns >= par::hap_missing_geno )
    {
      include[i] = false;
    }


  ///////////////////////////////////////////////
  // 2 or more hets at any loci -> ambiguous
  // Haploid genotypes should never be heterozygous, 
  // so we are okay here w.r.t X chromosome
  
  int het=0;
  for (int s=0; s<ns; s++)
    if ( (!s1[s]) && s2[s])
      het++;
  if (het>1)
    ambig[i] = true;

  return;
}



bool HaploPhase::consistentNonfounderPhaseGivenGenotypes(vector<bool> & s1, 
							 vector<bool> & s2,
							 int h1, int h2)
{
  
  // This function works for autosomal, haploid and sex chromosomes

  // Template haplotypes
  
  vector<bool> & t1 = hap[h1];
  vector<bool> & t2 = hap[h2];
  
  for (int s=0; s<ns; s++)
    {
      
      // Ignore missing genotypes (observed; template will never be
      // missing)
      
      if ( s1[s] && !s2[s] ) 
	continue;
      
      // Template homozygous? 
      // (Haploid templates will always be homozygous)
      
      if ( t1[s] == t2[s] )
	{
	  if ( s1[s] != t1[s] ||
	       s2[s] != t2[s] )
	    return false;
	}
      else // heterozygous template
	{
	  if ( s1[s] == s2[s] )
	    return false;
	  
	}
    }

  // Looks like is does match
  return true;
}

bool HaploPhase::consistentNonfounderMalePhaseGivenXGenotypes(vector<bool> & s1, 
							      vector<bool> & s2, 
							      int h2)
{
  
  // This function works for haploid individuals (male X offspring);
  // only the mother transmitted the X
  
  // Template haplotypes
  
  vector<bool> & t1 = hap[h2];
  
  for (int s=0; s<ns; s++)
    {
      
      // Ignore missing genotypes (observed; template will never be
      // missing)
      
      if ( s1[s] && !s2[s] ) 
	continue;
      
      // Haploid templates will always be homozygous
      // Haploid genotype should always be homozygous too

      if ( s1[s] != t1[s] )
	return false;
      
    }

  // Looks like is does match
  return true;
}



bool HaploPhase::consistentNonfounderPhaseGivenParents(int i, 
						       int h1, int h2,      
						       int p1, int p2,      
						       int m1, int m2)
{
  
  // Given offspring haplotypes and parental haplotypes, identify
  // whether or not this offspring genotype is possible
  
  if ( X && P.sample[i]->sex )
    {

      // We should only have specified possible homozygous phases --
      // therefore, if male X chr, we only need to check that it is
      // consistent with at least one maternal X
      
      if ( h1 == m1 || h1 == m2 )
	return true;
      
    }
  
  return (h1 == p1 && h2 == m1 ) ||
    (h1 == p1 && h2 == m2 ) || 
    (h1 == p2 && h2 == m1 ) ||
    (h1 == p2 && h2 == m2 ) || 
    (h1 == m1 && h2 == p1 ) ||
    (h1 == m1 && h2 == p2 ) || 
    (h1 == m2 && h2 == p1 ) ||
    (h1 == m2 && h2 == p2 ); 

}


void HaploPhase::resolveWithKids(int i)
{

  // Consider the founders in each family, who 
  // have at least 1 child, and a genotyped spouse

  // We require a full family, with two parents
  //  if ( ! f->parents ) return;

  //   int pati = pat->

  //   Individual * pat =  f->pat;
  //   Individual * mat =  f->mat;

  //   A/a B/b A/a B/b -> A/A B/B
  //                      AB / AB 

  //   for (int i=0; i< P.family[f].size(); i++)
  //     cout << P.family[f]->fid << "\t"
  // 	 << P.family[f]->iid << "\t"
  // 	 << P.family[f]->pat->iid << "\t"
  // 	 << P.family[f]->mat->iid << "\n";

}


void HaploPhase::phaseAndScoreNonfounder(int i)
{

  //////////////////////////////////////////////
  // Always try to phase this offspring

  include[i] = true;

  //////////////////////////////////////////////
  // Link this individual up with their parents

  int father = P.sample[i]->ip;
  int mother = P.sample[i]->im;

  bool nofather = false;
  bool nomother = false;
  
  if (father==-1)
    nofather = true;
  else if (!include[father])
    nofather = true;

  if (mother==-1)
    nomother = true;
  else if (!include[mother])
    nomother = true;

  // For TDT purposes, we require both parents to be 'observed'
  // i.e. so we never we to consider the "AllPhases" list (so
  //      we now do not bother generating it, i.e. enumerateAllPhases()
  //      function call is commented out in the main loop above
  

  if ( nofather || nomother )
    {
      include[i] = false;
      return;
    }
  
  
  int pat_phases = hap1[father].size();
  int mat_phases = hap1[mother].size();

  // Too much ambiguity?
  
  if (pat_phases * mat_phases >= par::hap_max_nf_phases )
    {
      include[i] = false;
      return;
    }
  

  // Keep track of transmitted and non-transmitted 
  // haplotypes if performing a TDT-type analysis


  vector<vector<int> > trans1(0);
  vector<vector<int> > untrans1(0);

  /////////////////////////////////////////
  // Perform fill-in phasing for offspring

  // Step 1. Enumerate possible offspring phases
  //         or set to not include if too much missing, 
  //         and populate s1/s2 with genotype data for
  //         region

  vector<bool> s1(ns);
  vector<bool> s2(ns);


  validateNonfounder(i,s1,s2);


  //////////////////////////////////////////////
  // Do we want to attempt to reconstruct phase?

  if ( ! include[i] )
    {
      return;
    }
  

  ////////////////////////////////////////////////
  // Step 2. Joint distribution of parental phases

  double psum = 0;
  int pcnt=1;


  // Set offspring posterior probability list to nil

  pp[i].clear();
  
  
  // Consider all possible pairs of parental phases and implied
  // possible haplotypic transmissions

  // If no mother or father exists, we are using the standard
  // ph[] enumeration of all possible haplotypes: NOT SUPPORTED
  // CURRENTLY, BUT WE COULD IMPLEMENT AGAIN FOR SIBLINGS

  vector<int> & pathap1 = nofather ? ph_hap1 : hap1[father];
  vector<int> & pathap2 = nofather ? ph_hap2 : hap2[father];
  
  vector<int> & mathap1 = nomother ? ph_hap1 : hap1[mother];
  vector<int> & mathap2 = nomother ? ph_hap2 : hap2[mother];
  
  map<FamilyTransmissions,double> & pmap = phasemap[i];
  
  for (int z1=0; z1 < pat_phases ; z1++)
    for (int z2=0; z2 < mat_phases ; z2++)
      {

	int p1 = pathap1[z1];
	int p2 = pathap2[z1];
	
	int m1 = mathap1[z2];
	int m2 = mathap2[z2];
	
	// Legacy code: we no longer take this approach, but 
	// the code is left here to show how to call the function

	// 	// Obtain possible offspring phases, given offspring
	// 	// genotypes and parental haplotypes
	// 	enumerateNonfounderPhase(i,             // offspring individual
	// 				 s1, s2,        // offspring genotypes 
	// 				 p1, p2,        // paternal haplotypes
	// 				 m1, m2,        // maternal haplotypes
	// 				 phap1, phap2); // return possible offspring haplotypes
	

	// Given parental phases, there are four possible autosomal
	// offspring tranmissions (for autosomes). We should enumerate
	// these and see which are consistent with the observed
	// offspring genotypes
	
	// Autosome(PM)  Haploid*    X(->female)    X(->male)

	// 00            0          00             *0
	// 01            1          01             *1
	// 10
	// 11
	
	// * not implemented; i.e. for now these are skipped
	// (i.e. this function is never called) -- this will be added
	// in future versions; haploid genotypes are coded as
	// homozygous; but for the haploid case, we need special code
	// in place to indicate MT transmission, etc.
	
	// For the autosomal X case, when transmitting to males, we
	// need a special function, however, as we do not want to
	// consider at all the paternal (homozygous/haploid) X
	// genotype (i.e. Y was transmitted...)
	
	// For female offspring on the X, we do want to look at the
	// paternal X for concordance, but we should only look at one
	// copy (as the father should always be homozygous/haploid).
	  
	vector<FamilyTransmissions> offspring(4);
	vector<bool> possible(4,false);
	
	int npossible = 0;
	
	for (int tr_pat = 0; tr_pat < 2; tr_pat++)
	  for (int tr_mat = 0; tr_mat < 2; tr_mat++)
	    {
	      
	      // Handle special cases of non-autosomal chromosomes
	      
	      if ( X )
		{ 
		  // If boy, haploid and X must have come from mother
		  // (0 paternal possible transmissions) 
		  
		  // If girl, diploid, but father can only send one
		  // possible X (only 1 possible paternal
		  // transmission)
		  
		  if ( tr_pat == 1 )
		    continue;
		}	      
	      
	      // Offspring haplotypes
	      
	      // Store: c is 0..3 coding

	      //        first two are for maternal transmissions
	      //        so we can take a short cut and just consider
	      //        first two positions for X (where there is 
	      //        no variation in paternal transmission conditional 
	      //        on offspring sex) 

	      int c = tr_mat + tr_pat*2;
			      
	      // Paternal transmission
	      
	      if ( tr_pat == 0 ) 
		{
		  offspring[c].pt = p1;
		  offspring[c].pu = p2;
		}
	      else
		{
		  offspring[c].pt = p2;
		  offspring[c].pu = p1;
		}

	      // Maternal transmission

	      if ( tr_mat == 0 ) 
		{
		  offspring[c].mt = m1;
		  offspring[c].mu = m2;
		}
	      else
		{
		  offspring[c].mt = m2;
		  offspring[c].mu = m1;
		}
		

	      // Is this offspring phase compatible with the offspring
	      // genotypes?

	      if ( X && P.sample[i]->sex ) 
		{
		  
		  // Only consider maternal X transmission to male
		  if ( consistentNonfounderMalePhaseGivenXGenotypes(s1,s2,
								    offspring[c].mt) )
		    {
		      possible[c] = true;
		      npossible++;
		    }
		}	     
	      else if ( consistentNonfounderPhaseGivenGenotypes(s1,s2,
							   offspring[c].pt,
							   offspring[c].mt) )
		{
		  
		  // Add this to list of possible offspring phases,
		  // keeping track of frequency
		  
		  // If we were to revert to using absent parents,
		  // then ph_freq[] should have been populated by
		  // enumerateAllPhases()
		  
		  possible[c] = true;
		  npossible++;
		  
		}
	    
	    } // Next of 2/4 (max) possible parental transmissions		  
	

	// Need to scale these 0->4 possibilities to sum to correct
	// value
	
	// At least one possible phase?

	if ( npossible > 0 )
	  {
	    
	    double p = 1;
	    
	    if (ambig[father])
	      p *= pp[father][z1];
	    
	    if (ambig[mother])
	      p *= pp[mother][z2];
	    
	    // We explicitly consider both phases, so we remove this line
	    //		if (h1!=h2)
	    //		  p *= 2;
	    
	    p /= (double)npossible; 
	    
	    int numposs = 4;       // Autosomal
	    if ( X ) numposs = 2;  // X transmission

	    for (int j=0; j<numposs; j++)
	      {
		if ( possible[j] ) 
		  {
		    
		    map<FamilyTransmissions,double>::iterator 
		      ip = pmap.find( offspring[j] );
		    
		    if ( ip == pmap.end() )
		      pmap.insert( make_pair(offspring[j] , p) );
		    else
		      ip->second += p;
		    
		    // Keep track of total probability
		    
		    psum += p;
		    
		  }
		
	      }
	  }
	// Consider next possible parental phase
      }


  /////////////////////////////////////
  // Extract possible offspring phases and populate 
  // standard metrics
  
  pp[i].clear();
  hap1[i].clear();
  hap2[i].clear();


  ///////////////////////////////////////////////////////
  // Store the possible offspring phases (transmissions 
  // only), and keep track of the probabilities

  map<FamilyTransmissions,double>::iterator ip = pmap.begin();
  
  include[i] = ambig[i] = true;
  
  if ( pmap.size() == 0 )
    {
      include[i] = false;
    }
  else if ( pmap.size() == 1 )
    {
      ambig[i] = false;
      
      int h1 = ip->first.pt;
      int h2 = ip->first.mt;
      
      if ( h1 < h2 ) 
	{
	  hap1[i].push_back( h1 );
	  hap2[i].push_back( h2 );
	}
      else
	{
	  hap1[i].push_back( h2 );
	  hap2[i].push_back( h1 );
	}
    }
  else
    {
      // More than one possible phase for this offspring

      map<int2,int> mapBack;

      while ( ip != pmap.end() )
	{
	  int2 h;
	  h.p1 = ip->first.pt;
	  h.p2 = ip->first.mt;
      
	  if ( h.p2 < h.p1 ) 
	    {
	      int t = h.p1;
	      h.p1 = h.p2;
	      h.p2 = t;
	    }

	  // Have we already seen this pair of transmitted haplotypes?

	  map<int2,int>::iterator im = mapBack.find(h);

	  if ( im != mapBack.end() )
	    {
	      int k = im->second;
	      pp[i][k] += ip->second;
	    }
	  else
	    {
	      int t = pp[i].size();
	      mapBack.insert(make_pair(h,t));

	      pp[i].push_back( ip->second );

	      hap1[i].push_back( h.p1 );
	      hap2[i].push_back( h.p2 );	  
	      
	    }

	  // Next family transmission
	  ip++;
	}
    }

  
  ////////////////////////////
  // Normalise probabilities
  
  if (ambig[i])
    for (int z=0; z < pp[i].size(); z++)
      pp[i][z] /= psum;
  

  map<FamilyTransmissions,double>::iterator itp = pmap.begin();
  while ( itp != pmap.end() )
    {
      itp->second /= psum;
      ++itp;
    }


  
  ///////////////////////////////////////////////////////////
  // Score haplotype transmissions for this trio, and add to
  // tabulation of sample T and U counts

  if (par::test_hap_TDT || par::proxy_TDT)
    transmissionCount(i,pmap);

  return;
}

void HaploPhase::transmissionCount(int i, 
				   map<FamilyTransmissions,double> & pmap )
{
  

  // For debugging only:
  //  displayFamTran(pmap,i,this);


  map<FamilyTransmissions,double>::iterator ip = pmap.begin();

  int t = subhaplotypes ? nt : nh;
  
  
  ////////////////////////////////////
  // Consider each possible phase set
  
  while ( ip != pmap.end() )
    {
      
      vector<int> t1(t,0);
      vector<int> u1(t,0);
      
      FamilyTransmissions f = ip->first;
      double posterior = ip->second;

      // This function works fine for X chromosome as is.
      // i.e. fathers haploid/homozygous/uninformative;
      // son's/daughters genotype will always reflect X maternal
      // transmission
      
      int h1, h2, p1, p2, m1, m2;

      if ( subhaplotypes ) 
	{
	  
	  // Collapse from a 0..nh space to a 0..nt space, via
	  // downcoding<> We can assume the haplotype codes given here
	  // will always be valid (i.e. map between 0 and nh) and that the
	  // downcoding map will always have an appropriate key

	  //  AACCA  0   0  -AC--
	  //  ACCAC  1   1  -XX--
	  //  AACCC  2   0  -AC--
	  //  CCCCC  3   1  -XX--
	  
	  h1 = downcoding.find( f.pt )->second;
	  h2 = downcoding.find( f.mt )->second;

	  p1 = downcoding.find( f.pt )->second;
	  p2 = downcoding.find( f.pu )->second;

	  m1 = downcoding.find( f.mt )->second;
	  m2 = downcoding.find( f.mu )->second;
	  
	}
      else
	{
	  h1 = f.pt;
	  h2 = f.mt;

	  p1 = f.pt;
	  p2 = f.pu;

	  m1 = f.mt;
	  m2 = f.mu;
	}


      scoreTransmissions(h1,h2,p1,p2,m1,m2,t1,u1); 
            
      
      ///////////////////////////////////////////
      // Update sample totals for each haplotype
      // and also accumulate the empirical variance
      // of the transmissions
           
      for (int h=0; h<t; h++)
	{
	  trans[h] += t1[h] * ip->second;
	  untrans[h] += u1[h] * ip->second;
	}
      
      
      // Consider next family transmission set
      ++ip;
    }
    

}


void HaploPhase::scoreTransmissions(int h1, int h2,
				    int p1, int p2,
				    int m1, int m2,
				    vector<int> & t1, 
				    vector<int> & u1)
{

  // Return of vector of T and U (0,1,2) for each haplotype
  // for this particular trio
    
  // Father heterozygous? 
  if ( p1 != p2 )
    {
      // Mother homozygous?
      if ( m1 == m2 )
	{
	  // then select a kid allele that matches, 
	  // and score the other one for transmission
	  if ( h1 == m1 )
	    {
	      t1[h2]++;
	      if (p1==h2)
		u1[p2]++;
	      else
		u1[p1]++;
	    }
	  else
	    {
	      t1[h1]++;
	      if (p1==h1)
		u1[p2]++;
	      else
		u1[p1]++;
	      
	    }
	}
      else
	{
	  // Both parents are heterozygous, 
	  
	  // Transmitted alleles are unambiguous

	  t1[h1]++;
	  t1[h2]++;
	  
	  // Untransmitted alleles
	  // i.e. which two are left over 
	  // after accounting for the two
	  // transmitted alleles
	  
	  bool pat_accounted = false;
	  bool mat_accounted = false;

	  if (p1 != h1 && p1 != h2 )
	    {
	      u1[p1]++;
	      pat_accounted = true;
	    }
	  else if (p2 != h1 && p2 != h2 )
	    {
	      u1[p2]++;
	      pat_accounted = true;
	    }
	  
	  if (m1 != h1 && m1 != h2 )
	    {
	      u1[m1]++;
	      mat_accounted = true;
	    }
	  else if (m2 != h1 && m2 != h2 )
	    {
	      u1[m2]++;
	      mat_accounted = true;
	    }
	  
	  // This only happens with AB x AB -> AB
	  if ( ! (pat_accounted || mat_accounted ))
	    {
	      u1[h1]++;
	      u1[h2]++;
	    }
	  else if ( ( ( !pat_accounted ) &&   mat_accounted ) ||
		    (    pat_accounted   && (!mat_accounted) ) )
	    {

	      // If only 1 untransmitted allele accounted for, it must
	      // be the doubled allele that is untransmitted
	      
	      if (p1 == m1 || p1 == m2 )
		u1[p1]++;
	      else
		u1[p2]++;
	    }
	  
	  // AB  AB    BB -- 2 accounted for
	  //   AA
	  
	  // AB  AB    AB -- 0 accounted for
	  //   AB 
	  
	  // AB  AC    BA -- 1 accounted for
	  //   AC
	  
	  // AC  AB    BA -- 1 accounted for
	  //   AC
	  
	  // AB  CD    DB -- 2 accounted for 
	  //   AC
	  
	}
    }
  else if ( m1 != m2 )
    {
      // Mother heterozygous, father homozygous
      if (h1 == p1 )
	{
	  t1[h2]++;
	  if (m1==h2)
	    u1[m2]++;
	  else
	    u1[m1]++;
	}
      else
	{
	  t1[h1]++;
	  if (m1==h1)
	    u1[m2]++;
	  else
	    u1[m1]++;
	}
    }
    
  return;
}

