

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
#include <vector>
#include <map>
#include <assert.h>

#include "plink.h"
#include "options.h"
#include "phase.h"
#include "helper.h"
#include "genogroup.h"
#include "haplowindow.h"

extern ofstream LOG;

using namespace std;


string HaploPhase::haplotypeName(int h)
{
  string str;
  for (int s=0; s<ns; s++)
    {

      string a1 = P.locus[S[s]]->allele1;
      string a2 = P.locus[S[s]]->allele2;
      if ( a1 == "" ) a1 = "X";
      if ( a2 == "" ) a2 = "X";

      if (h == -1)
	str += "-"; // haploid gap
      else if (hap[h][s])
	str += a1;
      else
	str += a2;
    }
  return str;  
}



void HaploPhase::imputeAllHaplotypes()
{
  
  ///////////////////////////////////////////////
  // Impute missing SNPs -- create a new datafile

  // Imputation rules:

  // Missing predictor allele -> missing haplotype
  // P(H|G) < 0.8 (default) -> missing haplotype

  // Make space to new, imputed haplotype calls

  new_one.resize(P.n);
  new_two.resize(P.n);

  /////////////////////////////////
  // Phase all specified haplotypes

  phaseAllHaplotypes(true,*P.pperm);


  ///////////////////////////
  // Write new PED file

  P.printLOG("Imputing genotypes with P(H|G) threshold of " + dbl2str( par::hap_post_prob ) + "\n\n"); 

  string filename = par::output_file_name + ".impute.ped";

  P.printLOG("Writing imputed ped file to [ " + filename + " ] \n");

  ofstream PED(filename.c_str(), ios::out);
  PED.clear();

  for (int i=0; i<P.n; i++)
    {

      Individual * person = P.sample[i];
      PED << person->fid<< " "<< person->iid<< " "<< person->pat<< " "
	  << person->mat<< " "<< person->sexcode<< " ";

      if (par::bt)
	PED << (int)person->phenotype;
      else
	PED << person->phenotype;

      for (int l=0; l<new_one[i].size(); l++)
	{
	  if ( (!new_one[i][l]) && (!new_two[i][l]))
	    PED << par::recode_delimit<< actual_map[l]->allele1<< " "
		<< actual_map[l]->allele1;
	  else if ( (!new_one[i][l]) && new_two[i][l])
	    PED << par::recode_delimit<< actual_map[l]->allele1<< " "
		<< actual_map[l]->allele2;
	  else if (new_one[i][l] && new_two[i][l])
	    PED << par::recode_delimit<< actual_map[l]->allele2<< " "
		<< actual_map[l]->allele2;
	  else
	    PED << par::recode_delimit<< par::missing_genotype << " "
		<< par::missing_genotype;
	}
      PED << "\n";
    }

  PED.close();


  //////////////////
  // Write new map

  filename = par::output_file_name + ".impute.map";
  P.printLOG("Writing imputed map file to [ " + filename + " ] \n");

  ofstream MAP(filename.c_str(), ios::out);
  MAP.clear();

  for (int l=0; l < actual_map.size(); l++)
    {
      MAP << chromosomeName(actual_map[l]->chr) << "\t"
	  << actual_map[l]->name<< "\t"
	  << actual_map[l]->pos<< "\t"
	  << actual_map[l]->bp<< "\n";
    }

  MAP.close();

}

void HaploPhase::calculateHaplotypeFrequencies()
{

  string f = par::output_file_name + ".frq.hap";

  if (par::display_hap_freqs)
    {
      P.printLOG("Writing haplotype frequencies to [ " + f + " ]\n");
      HFRQ.open(f.c_str(), ios::out);
      HFRQ.precision(4);

      HFRQ << setw(10) << "LOCUS"<< " "<< setw(12) << "HAPLOTYPE"<< " "
	   << setw(10) << "F"<< "\n";
    }

  // Phase all SNPs (with frequency flag set, this routine
  // will write haplotype frequencies to HFRQ

  P.haplo->phaseAllHaplotypes(true,*P.pperm);

  if (par::display_hap_freqs)
    HFRQ.close();

  // So we do not re-write them
  par::display_hap_freqs = false;

}


void HaploPhase::imputeThisHaplotype(int l)
{

  ////////////////////////////////////////////////////
  // Impute for all individiuals, if common haplotype
  
  double w = 0;
  double c = 0;
  
  if (testHaplotypeFreq() >= par::min_hf && 
      testHaplotypeFreq() <= par::max_hf )
    {
      
      for (int i=0; i<P.n; i++)
	{
	  bool b1, b2;
	  
	  // Haplotype-weighting
	  
	  double t = imputeHaplotypes(i, b1, b2);
	  	  
	  if (t>=0)
	    {
	      w += t;
	      c++;
	    };
	  
	  // And set new imputed genotypes
	  new_one[i].push_back(b1);
	  new_two[i].push_back(b2);

	}
      
      /////////////////////////////////////////////
      // Add to map of actually imputed haplotypes
      
      Locus * loc = new Locus;
      loc->chr = new_map[l]->chr;
      loc->name = hname+"_"+haplotypeName(test_hap)+"_";
      loc->bp = new_map[l]->bp;
      loc->pos = new_map[l]->pos;
      loc->allele1 = new_map[l]->allele1;
      loc->allele2 = new_map[l]->allele2;
      actual_map.push_back( loc );
    }
}

/////////////////////////////////////////////
// Legacy function: now redundant

void HaploPhase::enumerateAllPhases()
{
  
  // Note: this function is no longer called
  
  // For individuals w/out parents: make a list of all possible
  // phases. Note: currently, we do not use this (i.e. we always
  // require father and mother to be 'observed' (i.e. genotyped and
  // adequately phased).
  
  // Also note: issue with representing heterozygote haplotypes twice
  // in list: previously we did not, but we now change this (seeing as
  // it is never used in any case...)
  
  // Also: now we build separate lists for diploid and haploid
  // chromosomes, not that we use either.
  
  // Diploid possible phases
  
  if ( !haploid )
    {
      for (int h1=0; h1<nh; h1++)
	for (int h2=0; h2<nh; h2++)
	  {
	    double freq = f[h1] * f[h2];
	    if (h1!=h2)
	      freq *= 2;
	    if (freq >= par::hap_min_phase_prob )
	      {
		ph_freq.push_back(freq );
		ph_hap1.push_back(h1 );
		ph_hap2.push_back(h2 );
	      }
	  }
    }
  
  // Haploid possible phases
  
  if (haploid || X )
	{
	  for (int h1=0; h1<nh; h1++)
	    {
	      double freq = f[h1];
	      if (freq >= par::hap_min_phase_prob )
		{
		  haploid_ph_freq.push_back(freq );
		  haploid_ph_hap1.push_back(h1 );
		}
	    }
	}
  
	// Original total number of phases 
  
  np = ph_hap1.size();
  haploid_np = haploid_ph_hap1.size();
  
}

// Legacy function: now redundant
/////////////////////////////////////////////


//////////////////////////////////////////////////
// Return a list of all possible haplotype names

vector<string> HaploPhase::returnHaplotypes(vector<int> & slist)
{

  vector<string> str;

  enumerateHaplotypes(slist );

  for (int h=0; h<hap.size(); h++)
    {
      string hstr;
      for (int s=0; s<ns; s++)
	if (!hap[h][s])
	  hstr += P.locus[S[s]]->allele1;
	else if (P.locus[S[s]]->allele2=="")
	  hstr += "0";
	else
	  hstr += P.locus[S[s]]->allele2;
      str.push_back(hstr);
    }
  return str;
}


//////////////////////////////////////////////////
// For multi-marker imputation and certain tasks, 
// we require a 'test' haplotype

void HaploPhase::setTestHaplotype(string t)
{

  // No specified test haplotype?
  if ( t == "" ) 
    {
      test_hap = -1;
      return;
    }
  
  // Create match template
  vector<bool> tmp(ns, false);
  for (int s=0; s<ns; s++)
    if (P.locus[S[s]]->allele1 == t.substr(s, 1) )
      tmp[s] = true;
  
  // Consider each haplotype 
  
  test_hap = -1;
  
  for (int h=0; h<hap.size(); h++)
    {
      bool match = true;
      
      for (int s=0; s<ns; s++)
	{
	  if (hap[h][s] != tmp[s])
	    {
	      match = false;
	      break;
	    }
	}
      
      if (match)
	{
	  test_hap = h;
	  break;
	}
    }
  
}

void HaploPhase::reportHaplotypeFrequencies()
{
  for (int h=0; h<nh; h++)
    {
      if (f[h] >= par::min_hf)
	{
	  HFRQ << setw(10) << hname << " "<< setw(12) << haplotypeName(h)
	       << " "<< setw(10) << f[h]<< "\n";
	}
    }
}

void HaploPhase::reportPhase()
{
  
  string fn = par::output_file_name+".phase-"+hname;
  ofstream PHASE(fn.c_str(), ios::out);

  P.printLOG("Writing phased haplotypes for "+ hname + " to [ "+ fn + " ]\n");
  
  PHASE << setw(par::pp_maxfid) << "FID" << " " 
	<< setw(par::pp_maxiid) << "IID"<< " "
	<< setw(4) << "PH"<< " "
	<< setw(10) << "HAP1"<< " "
	<< setw(10) << "HAP2"<< " "
	<< setw(12) << "POSTPROB"<< " "
    //  << setw(12) << "WEIGHT" << " "
	<< setw(6) << "BEST"<< " "<< "\n";
  
  PHASE.precision(4);
  
  for (int i = 0; i < P.n; i++)
    {
      
      if (include[i])
	{
	  
	  for (int z = 0; z < hap1[i].size(); z++)
	    {
	      
	      PHASE << setw(par::pp_maxfid) << P.sample[i]->fid<< " "
		    << setw(par::pp_maxiid) << P.sample[i]->iid<< " "
		    << setw(4) << z << " "
		    << setw(10) << haplotypeName(hap1[i][z]) << " ";
	      
	      if (haploid || (X && P.sample[i]->sex))
		PHASE << setw(10) << haplotypeName( -1 ) << " ";
	      else
		PHASE << setw(10) << haplotypeName(hap2[i][z]) << " ";
	      
	      if (ambig[i])
		{
		  PHASE << setw(12) << pp[i][z]<< " ";
		  int max_z = 0;
		  for (int z2=0; z2<hap1[i].size(); z2++)
		    max_z = pp[i][z2] > pp[i][max_z] ? z2 : max_z ;
		  
		  if (max_z == z)
		    PHASE << setw(6) << 1<< " "<< "  ";
		  else
		    PHASE << setw(6) << 0<< " "<< "  ";
		}
	      else
		PHASE << setw(12) << 1<< " "<< setw(6) << 1<< " "<< "  ";
	      
	      // Genotypes
	      //for (int s=0; s<ns; s++)
	      //PHASE << genotype(P, i, S[s]) << " ";
	      PHASE << "\n";
				
	    }
	}
      
      // Report also on excluded individuals
      // (Should be 0-size phase-set)
      else
	{
	  
	  PHASE << setw(par::pp_maxfid) << P.sample[i]->fid << " "
		<< setw(par::pp_maxiid) << P.sample[i]->iid << " "
		<< setw(4) << "NA" << " "
		<< setw(10) << "NA" << " "
		<< setw(10) << "NA" << " "
		<< setw(12) << "NA" << " "
		<< setw(6) << "NA"<< "   ";

	  // genotypes
	  // for (int s=0; s<ns; s++)
	  //   PHASE << genotype(P, i, S[s]) << " ";
	  
	  PHASE << "\n";
			
	}
    }
  
  PHASE.close();
  
}

void HaploPhase::reportPhaseWideFormat()
{
  string fn = par::output_file_name+".wphase-"+hname;
  ofstream PHASE(fn.c_str(), ios::out);
  
  P.printLOG("Writing wide-format phased haplotypes for "+ hname + " to [ "
	     + fn + " ]\n");
  
  PHASE << setw(par::pp_maxfid) << "FID"<< " "<< setw(par::pp_maxiid)
	<< "IID"<< " ";
  
  for (int h=0; h<nh; h++)
    if (f[h] >= par::min_hf)
      PHASE << setw(8) << "H_"+haplotypeName(h) << " ";
  
  PHASE << "\n";
  
  PHASE.precision(4);
  
  for (int i = 0; i < P.n; i++)
    {
      if (include[i])
	{
	  
	  PHASE << setw(par::pp_maxfid) << P.sample[i]->fid<< " "
		<< setw(par::pp_maxiid) << P.sample[i]->iid<< " ";
	  
	  vector_t hcnt(nh, 0);
	  
	  for (int z = 0; z < hap1[i].size(); z++)
	    {
	      
	      if (ambig[i])
		{
		  hcnt[hap1[i][z]] += pp[i][z];
		  if ( ! (haploid || (X && P.sample[i]->sex)))
		    hcnt[hap2[i][z]] += pp[i][z];
		}
	      else
		{
		  hcnt[hap1[i][z]] ++;
		  if ( ! (haploid || (X && P.sample[i]->sex)))
		    hcnt[hap2[i][z]] ++;
		}
	      
	    }
	  
	  for (int h=0; h<nh; h++)
	    if (f[h] >= par::min_hf)
	      PHASE << setw(8) << hcnt[h]<< " ";
	  PHASE << "\n";
	  
	}
      
      // Report also on excluded individuals
      // (Should be 0-size phase-set)
      
      else
	{
	  PHASE << setw(par::pp_maxfid) << P.sample[i]->fid<< " "
		<< setw(par::pp_maxiid) << P.sample[i]->iid<< " ";
	  for (int h=0; h<nh; h++)
	    if (f[h] >= par::min_hf)
	      PHASE << setw(8) << "NA"<< " ";
	  PHASE << "\n";
	}
    }
  
  PHASE.close();
  
}


map<int,int> HaploPhase::makeSubHaplotypeSet(boolvec_t & mask)
{
  map<int,int> t;
  map<boolvec_t,int> shap;
  int cnt=0;
  for (int h=0; h < nh; h++)
    {

      boolvec_t sh;
      
      for (int s = 0; s < ns ; s++)
	{
	  if (mask[s])
	    sh.push_back(hap[h][s]);
	}

      map<boolvec_t,int>::iterator si = shap.find(sh);
      if ( si == shap.end() )
	{
	  shap.insert(make_pair(sh,cnt));
	  t.insert(make_pair(h,cnt));
	  ++cnt;
	}
      else
	t.insert(make_pair(h, si->second));
    }
      
  return t;
  
}

map<int,int> HaploPhase::makeTestSet(boolvec_t & mask, boolvec_t & allele)
{
  
  map<int,int> tests;
  
  for (int h2=0; h2 < nh; h2++)
    {
      bool is_A = true;
      for (int s = 0; s < ns ; s++)
	{
	  if (mask[s] && hap[h2][s] != allele[s])
	    is_A = false;
	}
      
      if (is_A )
	tests.insert(make_pair(h2, 0));
      else
	tests.insert(make_pair(h2, 1));
    }
  
  return tests;
}

string HaploPhase::getSubHaplotypeName(boolvec_t & mask, boolvec_t & allele,
				       int blank)
{
  
  string str = "";
  
  for (int s=0; s < ns; s++)
    {
      if (s == blank )
	str += " ";
      else if (mask[s])
	{
	  if (allele[s])
	    str += P.locus[ S[s] ]->allele1;
	  else
	    str += P.locus[ S[s] ]->allele2;
	}
      else
	str += ".";
    }
  return str;
}

vector_t HaploPhase::imputeGenotype(int i, int l)
{
  
  // Probability of AA, AB and BB for position 'l' 
  // (of ns SNPs) for individual 'i'
  
  vector_t g(3);
  
  if (X || haploid )
    {
      g[0] = g[1] = g[2] = 0;
      return g;
      
      error("HaploPhase::imputeGenotypess() not yet set up for X \n");
    }
  
  // Not able to be imputed?
  
  if (!include[i])
    {
      g[0] = g[1] = g[2] = 0;
      return g;
    }
  
  // Unambiguous imputation?
  
  if (!ambig[i])
    {
      
      int h1 = hap1[i][0];
      int h2 = hap2[i][0];
      
      bool s1 = hap[h1][l];
      bool s2 = hap[h2][l];
      
      if (s1 != s2 )
	g[1] = 1;
      else if (s1 )
	g[0] = 1;
      else
	g[2] = 1;
      
      return g;
    }
  
  // Weighted, ambiguous imputation?
  
  for (int z=0; z<hap1[i].size(); z++)
    {
      
      // ?? include?? if (pp[i][max_z] >= par::hap_post_prob)
      
      int h1 = hap1[i][z];
      int h2 = hap2[i][z];
      
      bool s1 = hap[h1][l];
      bool s2 = hap[h2][l];
      
      if (s1 != s2 )
	g[1] += pp[i][z];
      else if (s1 )
	g[0] += pp[i][z];
      else
	g[2] += pp[i][z];
      
    } // next possible phase

  return g;
  
}

double HaploPhase::imputeHaplotypes(int i, bool & n1, bool & n2)
{

//   if ( X || haploid )
//     error("HaploPhase::imputeHaplotypes() not yet set up for X \n");
  
  bool actualX = X && P.sample[i]->sex;

  //////////////////////////////////////////////
  // Based on P(H|G) impute inferred haplotypes
  // for above-threshold individuals  
  
  // Not imputed

  double w = -1;

  if ( ! include[i] ) 
    {
      n1 = true;
      n2 = false;
      return w;
    }
  
  
  // First for individuals of unambiguous phase
  
  if (!ambig[i])
    {
      if (hap1[i][0] == test_hap)
	n1 = false;
      else
	n1 = true;
  
      if ( actualX || haploid ) 
	{
	  n2 = n1;
	}
      else
	{
	  if (hap2[i][0] == test_hap)
	    n2 = false;
	  else
	    n2 = true;
	}

      // Resolve potential het/missing coding confusion
      
      if (n1 && (!n2))
	{
	  n1 = false;
	  n2 = true;
	}
      
      // Unambiguous weighting (0, 1 or 2 copie of test_hap)
      if (!n1)
 	{
	  if ( actualX || haploid ) 
	    w=1;
	  else
	    {
	      if (!n2)
		w = 2;
	      else
		w = 1;
	    }
 	}
      else
 	w = 0;
    }
  else
    {
      
      // Second, for ambiguous individuals impute and assign weight
	  
      int max_z = 0;
      for (int z=0; z<hap1[i].size(); z++)
	max_z = pp[i][z] > pp[i][max_z] ? z : max_z ;
      
      // Set missing by default
      
      n1 = true;
      n2 = false;
      
      // Consider each phase z
      
      // Above threshold? 
      if (pp[i][max_z] >= par::hap_post_prob)
	{
	  
	  // Do we match 'test_hap' ( '1' allele ) 
	  // or not? ( '2' allele )
	  
	  if (hap1[i][max_z] == test_hap)
	    n1 = false;
	  else
	    n1 = true;

	  if ( actualX || haploid ) 
	    {
	      n2 = n1;
	    }
	  else
	    {
	      if (hap2[i][max_z] == test_hap)
		n2 = false;
	      else
		n2 = true;
	    }

	  // Resolve potential het/missing coding confusion
	  
	  if (n1 && (!n2))
	    {
	      n1 = false;
	      n2 = true;
	    }
	  
	  // Unambiguous weighting (1 or 2 copies of test_hap)
	  
	  // We are saying either 0, 1 or 2 copies
	  // Consider each haplotype
	  
	  // Number imputed / Actual number
	  for (int z=0; z<hap1[i].size(); z++)
	    {
	      if (hap1[i][z] == test_hap)
		w += pp[i][z];

	      if ( ! ( actualX || haploid )) 
		{
		  if (hap2[i][z] == test_hap)
		    w += pp[i][z];
		}
	    }
	  
	  w = pp[i][max_z]/ w;

	}
    }
  
  return w;

}



double HaploPhase::rsq_internal(int s1, int s2)
{
  // A convenience function for SNP x SNP r^2
  // i.e. here it does not matter which allele 
  // we consider, so just re-use mask

  if (s1 > ns || s2 > ns )
    error("Problem in rsq_internal(int,int)");

  boolvec_t m1(ns, false);
  boolvec_t m2(ns, false);

  m1[s1] = true;
  m2[s2] = true;

  return rsq_internal(m1, m1, m2, m2);
}

double HaploPhase::rsq_internal(boolvec_t & mask1, boolvec_t & alleles1,
				boolvec_t & mask2, boolvec_t & alleles2)
{
  
  // Assume f[] has been populated with sensible values
  // and hap[][] contains alleles
  
  if (mask1.size() != ns ||mask2.size() != ns ||alleles1.size() != ns
      ||alleles2.size() != ns )
    {
      cout << ns << " "
	   << mask1.size() << " "
	   << mask2.size() << " " 
	   << alleles1.size() << " " 
	   << alleles2.size() << "\n";

      error("Internal error in Phase::rsq");
    }
  //  ---X-X-  mask1
  //     0-0   alleles1
  
  //  ----X--  mask2
  //      1    alleles2

  // i.e. find r^2 between 00 haplotype made of SNPs 4 & 6 
  // from 7 SNP haplotype with allele 1 of SNP 5

  // Calculate frequency of first haplotype (fA)

  double fA = 0;
  double fB = 0;
  double fAB = 0, fAb = 0, faB = 0, fab = 0;

  for (int h = 0; h < nh; h++)
    {

      bool is_A = true;
      bool is_B = true;

      bool is_AB = true;
      bool is_Ab = true;
      bool is_aB = true;

      for (int s = 0; s < ns ; s++)
	{

	  if (mask1[s] && hap[h][s] != alleles1[s])
	    is_A = false;

	  if (mask2[s] && hap[h][s] != alleles2[s])
	    is_B = false;

	  if ( (mask1[s] && hap[h][s] != alleles1[s])
	       || (mask2[s] && hap[h][s] != alleles2[s]))
	    is_AB = false;
	  
	  if ( (mask1[s] && hap[h][s] != alleles1[s])
	       || (mask2[s] && hap[h][s] == alleles2[s]))
	    is_Ab = false;
	  
	  if ( (mask1[s] && hap[h][s] == alleles1[s])
	       || (mask2[s] && hap[h][s] != alleles2[s]))
	    is_aB = false;
	  
	}

      if (is_A )
	fA += f[h];

      if (is_B )
	fB += f[h];

      if (is_AB )
	fAB += f[h];
      else if (is_aB )
	faB += f[h];
      else if (is_Ab )
	fAb += f[h];
      else
	fab += f[h];

      // Next haplotype
    }

  double fa = 1 - fA;
  double fb = 1 - fB;

  // Calculate either r-sq or D'
  double D = fAB - fA * fB;

  if ( calculateDp ) 
    {
      double dmax1 = D > 0 ? fA * fb : fA * fB;
      double dmax2 = D > 0 ? fa * fB : fa * fb;
      double dmax = dmax1 < dmax2 ? dmax1 : dmax2;
      if ( dmax == 0 )
	return -1;
      return D / dmax;
    }
  else
    {      
      double denom = fA * fa * fB * fb;
      if (denom == 0)
	return -1;
      return (D*D) / denom;
    }
}

double HaploPhase::freq(boolvec_t & mask1, boolvec_t & alleles1)
{

  // Assume f[] has been populated with sensible values
  // and hap[][] contains alleles

  if (mask1.size() != ns ||alleles1.size() != ns )
    {
      cout << ns << " "
	   << mask1.size() << " "
	   << alleles1.size() << "\n"; 
	   
      error("Internal error in Phase::freq");
    }


  //  ---X-X-  mask1
  //     0-0   alleles1

  double fA = 0;

  for (int h = 0; h < nh; h++)
    {

      bool is_A = true;

      for (int s = 0; s < ns ; s++)
	{
	  if (mask1[s] && hap[h][s] != alleles1[s])
	    is_A = false;
	}

      if (is_A )
	fA += f[h];

      // Next haplotype
    }

  return fA;
}



double HaploPhase::rsq(int l1, int l2)
{

  reset();
  new_pred_locus.resize(1);
  new_map.resize(1);

  vector<int> twoSNPs(2);
  twoSNPs[0] = l1;
  twoSNPs[1] = l2;

  new_pred_locus[0] = twoSNPs;
  new_map[0] = P.locus[l1];

  bool old_silent = par::silent;
  par::silent = true;

  new_pred_allele = listPossibleHaplotypes(P, new_pred_locus[0]);

  phaseAllHaplotypes(true,*P.pperm);

  // hname = locus[l]->name;

  par::silent = old_silent;

  return rsq_internal(0, 1);

}

double HaploPhase::dprime(int l1, int l2)
{
  calculateDp = true;
  double dp = rsq(l1,l2);
  calculateDp = false;
  return fabs(dp);
}

void Plink::calcPairwiseLD()
{

  int l1 = getMarkerNumber(*this, par::ld_SNP1);
  int l2 = getMarkerNumber(*this, par::ld_SNP2);

  if (l1 == l2 )
    error("Cannot compute LD with self");

  if (l1 == -1)
    error("--ld {marker} {marker}: first marker not found");

  if (l2 == -1)
    error("--ld {marker} {marker}: second marker not found");

  printLOG("\nLD information for SNP pair [ "+ par::ld_SNP1 + " "
	   + par::ld_SNP2 + " ]\n\n");
  
  printLOG("   R-sq = " + dbl2str_fixed(haplo->rsq(l1, l2) , 3 ) + "     ");

  printLOG("D' = " + dbl2str_fixed(haplo->dprime(l1, l2) , 3 ) + "\n\n");
  
  printLOG("   Haplotype     Frequency    Expectation under LE\n");
  printLOG("   ---------     ---------    --------------------\n");

  for (int h=0; h < haplo->nh; h++)
    {

      printLOG("       " + haplo->haplotypeName(h) + "          " );

      printLOG(dbl2str_fixed( haplo->f[h]  ,3) + "            ");

      double e = 0;
      if ( haplo->haplotypeName(h) == locus[l1]->allele2 + locus[l2]->allele2 )
	e = (1 - locus[l1]->freq)*(1 - locus[l2]->freq);
      else if ( haplo->haplotypeName(h) == locus[l1]->allele1 + locus[l2]->allele2 )
	e = ( locus[l1]->freq)*(1 - locus[l2]->freq);
      else if ( haplo->haplotypeName(h) == locus[l1]->allele2 + locus[l2]->allele1 )
	e = (1 - locus[l1]->freq)*( locus[l2]->freq);
      else if ( haplo->haplotypeName(h) == locus[l1]->allele1 + locus[l2]->allele1 )
	e = ( locus[l1]->freq)*( locus[l2]->freq);

      printLOG(dbl2str_fixed( e ,3 ) + "\n");

    }
  printLOG("\n");

  int ch = 0;
  for (int h=0; h < haplo->nh; h++)
    if ( haplo->haplotypeName(h) == locus[l1]->allele2 + locus[l2]->allele2 )
      ch = h;

  // Is D positive or negative?

  string s;

  if ( haplo->f[ch] > (1 - locus[l1]->freq)*(1 - locus[l2]->freq) )
    s = locus[l1]->allele1 + locus[l2]->allele1 + "/" + locus[l1]->allele2 + locus[l2]->allele2;
  else
    s = locus[l1]->allele1 + locus[l2]->allele2 + "/" + locus[l1]->allele2 + locus[l2]->allele1;

  printLOG("   In phase alleles are " + s + "\n");

  
  return;
}



///////////////////////////////////////////////////////////////
//                                                           //
// For a particular pair of individuals, track the status    //
// of haplotype sharing across the chromosome/region; this   //
// is the driver function                                    //
//                                                           //
///////////////////////////////////////////////////////////////

void HaploPhase::trackSharedHaplotypes()
{

  // Find individual(s) to track

  p1 = -1;
  p2 = -1;
    
  for (int i=0; i<P.n; i++)
    {
      if ( P.sample[i]->fid == par::segment_haplotrack_fid1 && 
	   P.sample[i]->iid == par::segment_haplotrack_iid1 )
	{ 
	  p1 = i;
	}
      
      if ( P.sample[i]->fid == par::segment_haplotrack_fid2 && 
	   P.sample[i]->iid == par::segment_haplotrack_iid2 )
	{ 
	  p2 = i;
	}
      
      if ( p1 != -1 && p2 != -1 ) 
	break;
    }
  

  if ( p1 == -1 || p2 == -1 )
    {
      error("Problem finding individual(s) indicated in haplo-track option\n");
      return;
    }
  

  // Set whether looking at homozygosity of shared segments

  homozyg = p1 == p2;

  Individual * person1 = P.sample[p1];
  Individual * person2 = P.sample[p2];

  if ( homozyg ) 
    P.printLOG("\nReport for individual [ " + person1->fid + " " + person1->iid + " ]\n");
  else
    P.printLOG("\nReport for pair [ " + person1->fid + " " + person1->iid 
	       + ", "+ person2->fid + " " + person2->iid + " ]\n");    

  string f = par::output_file_name + ".shared";
  
  P.printLOG("Tracking shared haplotypes, writing output to [ " + f + " ]\n");
  HFRQ.open(f.c_str(), ios::out);
  HFRQ.precision(4);
  
  HFRQ << setw(4) << "CHR" << " "
       << setw(par::pp_maxsnp) << "SNP" << " "
       << "\n";

  trackedIBS.resize(P.nl_all);
  trackedN.resize(P.nl_all);

  ////////////////////
  // Do all the work

  P.haplo->phaseAllHaplotypes(true,*P.pperm);


  // Display

  for (int l=0; l<P.nl_all; l++)
    {
      HFRQ << P.locus[l]->name << "\t"
	   << P.locus[l]->bp << "\t"
	   << trackedIBS[l] << "\t"
	   << trackedN[l] << "\t"
	   << (double)trackedIBS[l]/(double)trackedN[l] << "\n";
    }

  if (par::display_hap_freqs)
    HFRQ.close();
 
}




///////////////////////////////////////////////////////////////
//                                                           //
// For a particular pair of individuals, track the status    //
// of haplotype sharing across the chromosome/region; this   //
// function does the actual work
//                                                           //
///////////////////////////////////////////////////////////////

void HaploPhase::trackThisSegment()
{
  
  // Are the chromosomes consistent with a shared segment at this position? 
  
  // No information? Then exit

  if ( ! ( include[p1] && include[p2] ) )
    return;

  if ( haploid 
       || (X && P.sample[p1]->sex) 
       || (X && P.sample[p2]->sex) )
    error("Cannot use haplo-track options on non-autosomal chromosomes yet");
  


  // Looking within an individual for homozygous segments?

  if ( homozyg ) 
    {

      double probHomozyg = 0;

      for (int z = 0; z < hap1[p1].size(); z++)
	{

	  // Is this region shared...

	  if ( hap1[p1][z] == hap2[p1][z] )
	    {

	      // ...and rare?

	      if ( f[ hap1[p1][z] ] < 0.02 ) 
		{
		  if ( ambig[p1] )
		    probHomozyg += pp[p1][z];
		  else
		    probHomozyg = 1;
		}
	    }
	}
    }
  
  else // ... or looking between individuals for shared segments?
    {

      double probShared = 0;
      int j=0;
      for (int z1 = 0; z1 < hap1[p1].size(); z1++)
	for (int z2 = 0; z2 < hap1[p2].size(); z2++)
        {
	  
	  // Figure IBS 0, 1 or 2
	  
	  int a1 = hap1[p1][z1];
	  int a2 = hap2[p1][z1];
	  
	  int b1 = hap1[p2][z2];
	  int b2 = hap2[p2][z2];
	  
	  double prob = 1;
	  if ( ambig[p1] )
	    prob *= pp[p1][z1];
	  if ( ambig[p2] )
	    prob *= pp[p2][z2];
	  
	  if ( a1 > a2 )
	    {
	      int tmp = a1;
	      a1 = a2;
	      a2 = tmp;
	    }
	  
	  if ( b1 > b2 )
	    {
	      int tmp = b1;
	      b1 = b2;
	      b2 = tmp;
	    }
	  
	  // Count up similar, rare haplotypes

	  int cnt =0 ;
	  if ( a1 == b1 && f[a1] < .2 ) 
	    {
	      probShared += 0.5 * prob;
	      cnt++;
	    }

	  if ( a2 == b2 && f[a2] < .2 ) 
	    {
	      probShared += 0.5 * prob;
	      cnt++;
	    }

	  
	  // Keep track
	  
	  for ( int s = 0; s < ns ; s++ ) 
	    {
	      trackedIBS[S[s]] += probShared;
	      trackedN[S[s]]++;
	    }
	  
	  // next pair of haplotypes
        } 

    }
}




///////////////////////////////////////////////////////////////
//                                                           //
// Return a set of haplotype number codes given the type of  //
// mask + allele template used in the proxy association      //
// procedures                                                //
//                                                           //
///////////////////////////////////////////////////////////////

set<int> HaploPhase::returnHaplotypeSet(boolvec_t & mask, 
				    boolvec_t & alleles)
{
  set<int> hs;
  for (int h = 0; h < nh; h++)
    {
      bool is_A = true;
      for (int s = 0; s < ns ; s++)
	{
	  if (mask[s] && hap[h][s] != alleles[s])
	    is_A = false;
	}
      if (is_A)
	hs.insert(h);
    }
  return hs;
}




void HaploPhase::calculateEmpiricalVariance(int h)
{
  set<int> hs;
  hs.insert(h);
  calculateEmpiricalVariance(hs);
}




///////////////////////////////////////////////////////////////
//                                                           //
// Post-phasing, for a group of haplotypes, return the       //
// empirical variance and ratio of this to asymptotic        //
// variance for all individuals                              //
//                                                           //
///////////////////////////////////////////////////////////////


void HaploPhase::calculateEmpiricalVariance(set<int> & hs)
{
  
  double frequency = 0;

  set<int>::iterator h = hs.begin();
  while ( h != hs.end() )
    {
      frequency += f[*h];
      ++h;
    }
 

  // Do we need to consider this haplotype/set of 
  // haplotypes?

  if( frequency < 0.0000001 )
    {
      ratio = 0;
      empiricalVariance = 0;
      return;
    }


  // Calculate theoretical variance of frequency given binomial
  
  double theoreticalVariance = frequency * ( 1 - frequency );
  
  double weightedVariance = 0;
  double dosageSSQ = 0;
  double haplotypeCount = 0;
  


  // Calculate empirical variance given imputed haplotype counts:

  int ncnt = 0; // for allele count
  int dosageCount = 0; // for dosage (individual) count

  for( int i = 0; i < P.n; i++ )
    {
      
      if ( include[i] )
	{
	  if ( (!X) || !P.sample[i]->sex )
	    {
	      if (!ambig[i])
		{
		  if( hs.find( hap1[i][0] ) != hs.end() )
		    haplotypeCount += 1;
		  if( hs.find( hap2[i][0] ) != hs.end() )
		    haplotypeCount += 1;
		  
		}
	      else
		for( int z = 0; z < pp[i].size(); z++ )
		  {
		    if( hs.find( hap1[i][z] ) != hs.end() )
		      haplotypeCount += pp[i][z];
		    if( hs.find( hap2[i][z] ) != hs.end() )
		      haplotypeCount += pp[i][z];
		  }		
	      
	      ncnt+=2;
	      dosageCount++;
	    }
	}
    }
  
  double mean =  haplotypeCount/(double)ncnt;
  double dmean = haplotypeCount/(double)dosageCount;
  
  // Ratio of variance of weighted versus variance of averages
  // (i.e. dosage -- this measures information loss, as the deflation
  // is only for the dosage))
  
  // Calculate variance:  S(x-mean)^2/(n-1)
  
  for( int i = 0; i < P.n; i++ )
    {
      if ( include[i] )
	{
	  if ( (!X) || !P.sample[i]->sex )
	    {
	      double dosage = 0;
	      
	      if (!ambig[i])
		{
		  
		  if( hs.find( hap1[i][0] ) != hs.end() )
		    {
		      weightedVariance += (1-mean) * (1-mean);
		      dosage++;
		    }
		  else
		    weightedVariance += mean*mean; // (0-mean)^2
		  
		  
		  if( hs.find( hap2[i][0] ) != hs.end() )
		    {
		      weightedVariance += (1-mean) * (1-mean);
		      dosage++;
		    }
		  else
		    weightedVariance += mean*mean;
		  
		  dosageSSQ += (dosage-dmean)*(dosage-dmean);
		}
	      else
		{
		  
		  // Variance based on weights
		  
		  for( int z = 0; z < pp[i].size(); z++ )
		    {
		      
		      if( hs.find( hap1[i][z] ) != hs.end() )
			{
			  weightedVariance += pp[i][z] * (1-mean) * (1-mean);
			  dosage += pp[i][z];
			}
		      
		      else
			weightedVariance += pp[i][z] * mean *mean;
		      
		      if( hs.find( hap2[i][z] ) != hs.end() )
			{
			  weightedVariance += pp[i][z] * (1-mean) * (1-mean);
			  dosage += pp[i][z];
			}
		      else
			weightedVariance += pp[i][z] * mean * mean;
		      
		    }
		  
		  dosageSSQ += (dosage-dmean)*(dosage-dmean);
		  
		}		
	    }
	}
    }
  

  // Use N, not N-1 denominator, as we are comparing to the expected
  // variance above (i.e. so ratio == 1 in case of complete
  // information)


  // Update variables in HaploPhase 

  weightedVariance /= (double)ncnt;
  empiricalVariance = dosageSSQ / ((double)dosageCount*2);
  ratio = theoreticalVariance > 0 ? empiricalVariance / theoreticalVariance : 0;
  
}




///////////////////////////////////////////////////////////////
//                                                           //
// Verbose display function for phasing routine              // 
//                                                           //
///////////////////////////////////////////////////////////////

void HaploPhase::verboseDisplayWindows(int i, bool use_ref )
{

  if ( ! include[i] ) 
    return;
  
  for (int w = startWindow; 
       w <= finishWindow ; w++)
    {
      
      int r = windows[w]->genoGroup[i]->reference;
      if ( ! use_ref ) 
	r = i;
      
      HaploWindow * thisWindow = windows[w];
      
      VPHASE << "WINDOW " << w << ": "
		    << windows[w]->start << " to " 
		    << windows[w]->stop 
		    << " ( " << windows[w]->ns << " SNPs )\n";

      for ( int s = 0 ; s < windows[w]->ns ; s++ ) 
	VPHASE << P.locus[ thisWindow->S[s]]->name << " ";
      VPHASE << "\n";

      // Display real genotypes

      VPHASE << setw(w) << " ";
      for (int s=0; s< thisWindow->ns; s++)
	{

	  bool s1 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->one[i] :
	    P.sample[i]->one[ thisWindow->S[s] ];
	  bool s2 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->two[i] :
	    P.sample[i]->two[ thisWindow->S[s] ];
	  
	  if ( s1 ) 
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele2 ;
	      else
		VPHASE << "-";
	    }
	  else
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele1 ;
	      else
		VPHASE << P.locus[ thisWindow->S[s] ]->allele1 ;	 	  
	    }
	}      
      
      VPHASE << " ";
      for (int s=0; s< thisWindow->ns; s++)
	{
	  bool s1 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->one[i] :
	    P.sample[i]->one[ thisWindow->S[s] ];
	  bool s2 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->two[i] :
	    P.sample[i]->two[ thisWindow->S[s] ];

	  
	  if ( s1 ) 
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele2 ;
	      else
		VPHASE << "-";
	    }
	  else 
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele2 ;
	      else
		VPHASE << P.locus[ thisWindow->S[s] ]->allele1 ;	 	  
	    }
	}

      VPHASE << "\n";
      
      VPHASE << setw(w) << " ";
      for (int s=0; s< thisWindow->ns; s++)
	{
	  bool s1 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->one[i] :
	    P.sample[i]->one[ thisWindow->S[s] ];
	  bool s2 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->two[i] :
	    P.sample[i]->two[ thisWindow->S[s] ];
	  
	  if ( s1 ) 
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele2 ;
	      else
		VPHASE << "-";
	    }
	  else
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele1 ;
	      else
		VPHASE << P.locus[ thisWindow->S[s] ]->allele1 ;	 	  
	    }
	}      
      VPHASE << " ";
      for (int s=0; s< thisWindow->ns; s++)
	{
	  bool s1 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->one[i] :
	    P.sample[i]->one[ thisWindow->S[s] ];
	  bool s2 = par::SNP_major ? 
	    P.SNP[ thisWindow->S[s] ]->two[i] :
	    P.sample[i]->two[ thisWindow->S[s] ];
	  
	  if ( s1 ) 
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele2 ;
	      else
		VPHASE << "-";
	    }
	  else 
	    {
	      if ( s2 ) 
		VPHASE << P.locus[ thisWindow->S[s] ]->allele2 ;
	      else
		VPHASE << P.locus[ thisWindow->S[s] ]->allele1 ;	 	  
	    }
	}

      VPHASE << "\n";


      for (int z = 0; z < windows[w]->hap1[r].size(); z++)
	{
	  
	  VPHASE << setw(w) << " "
			<< thisWindow->haplotypeName(thisWindow->hap1[r][z])
			<< "/"
			<< thisWindow->haplotypeName(thisWindow->hap2[r][z])
			<< " ";
	  VPHASE << "( " 
			<< thisWindow->f[ thisWindow->hap1[r][z] ] << " / " 
			<< thisWindow->f[ thisWindow->hap2[r][z] ] << " ) ";
	  
	  if ( thisWindow->hap1[r].size() == 1)
	    VPHASE << "[1]\n";
	  else
	    VPHASE << thisWindow->pp[r][z]<< "\n";
	  
	}
    }

}
