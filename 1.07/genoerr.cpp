

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2007 Shaun Purcell                  //
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

extern ofstream LOG;
using namespace std;


void HaploPhase::queryGenotype(int s)
{

  // Map within-region SNP coding back to within-genome SNP coding

  int l = S[s];
  

  // Right now, let's elect not to handle non-autosomal SNPs

  if ( par::chr_haploid[ P.locus[l]->chr ] ||
       par::chr_sex[ P.locus[l]->chr ] )
    error("--proxy-error not yet set up for non-autosomal SNPs");


  //////////////////////////////////////////////////////
  // Calculate sample genotype frequencies for test SNP
  
  int g0 = 0, g1 = 0, g2 = 0;

  for (int i = 0; i<P.n; i++)
    {
      
      bool s1 = par::SNP_major ? P.SNP[l]->one[i] : P.sample[i]->one[l];
      bool s2 = par::SNP_major ? P.SNP[l]->two[i] : P.sample[i]->two[l];
      
      if( s1 ) 
	{
	  if ( s2 ) 
	    ++g0;
	}
      else
	{
	  if ( s2 ) 
	    ++g1;
	  else
	    ++g2;
	}
    }

  vector_t geno_freq(3);  
  int total = g0 + g1 + g2;
  geno_freq[0] = (double)g0 / (double)total;
  geno_freq[1] = (double)g1 / (double)total;
  geno_freq[2] = (double)g2 / (double)total;
  


  //////////////////////////////////////////////////////
  // Check each genotype of each individual over region

  for (int i=0; i<P.n; i++)
    {
      if (include[i])
	{

	  bool s1 = par::SNP_major ? P.SNP[l]->one[i] : P.sample[i]->one[l];
	  bool s2 = par::SNP_major ? P.SNP[l]->two[i] : P.sample[i]->two[l];
	  
	  int g;

	  if ( s1 ) 
	    {
	      if ( s2 ) 
		g = 0;
	      else
		g = -9;
	    }
	  else
	    {
	      if ( s2 ) 
		g = 1;
	      else
		g = 2;
	    }
	  	

	  ///////////////////////////////////////////////////////////
	  // Determine most likely genotype given flanking haplotypes

	  queryThisGenotype(i, s, g, geno_freq);

	}
    }  
}


void HaploPhase::queryThisGenotype(int i, int s, int g, vector_t & geno_freq )
{
  
  // Do not do anything if the actual genotype is missing -- i.e. the
  // focus here is on correcting genotyping error rather than
  // imputation

  if ( g<0 ) 
    return;

  // For a given individual, for a given phased region and set of
  // sample haplotype frequencies, determine a) the possible set of
  // phased haplotypes consistent with region if the test SNP were in
  // fact missing, b) the relative probability that the observed
  // genotype is the true genotype given the new phases.
  
  ////////////////////////////////////////////////////////////////
  // These are the new possible phases, if the test genotype were
  // missing

  vector<int> new_hap1;
  vector<int> new_hap2;
  vector_t newpp;
  
  double psum = 0;
  int z2 = 0;

  // Consider 
  //    AAA-C-GGT / ACA-C-GGT

  //  Becomes up to three possible states:
  //    AAA-C-GGT / ACA-C-GGT
  //    AAA-C-GGT / ACA-T-GGT
  //    AAA-T-GGT / ACA-T-GGT

  for (int z=0; z<hap1[i].size(); z++)
    {
      
      vector<int> posshap1;
      vector<int> posshap2;
      
      vector<bool> h1 = hap[hap1[i][z]];
      vector<bool> h2 = hap[hap2[i][z]];
      vector<bool> h1_flip = hap[hap1[i][z]];
      vector<bool> h2_flip = hap[hap2[i][z]];
      
      // flip bit at SNP position s
      h1_flip[s] = !h1[s];
      h2_flip[s] = !h2[s];
      
      // Do we observed these possible haplotypes in the sample as a
      // whole? If not, no need to consider, as the posterior
      // probability will be 0.
      
      if( hapmapb.find(h1) != hapmapb.end() )
	posshap1.push_back(hapmapb[h1]);
      
      if( hapmapb.find(h2) != hapmapb.end() )
	posshap2.push_back(hapmapb[h2]);
      
      if( hapmapb.find(h1_flip) != hapmapb.end() )
	posshap1.push_back(hapmapb[h1_flip]);
      
      if( hapmapb.find(h2_flip) != hapmapb.end() )
	posshap2.push_back(hapmapb[h2_flip]);
      
      // get new probabilities for each possible new phasing
      for( int a = 0; a < posshap1.size(); a++ )
	for( int b = 0; b < posshap2.size(); b++ )
	  {
	    newpp.push_back(f[posshap1[a]] * f[posshap2[b]]);
	    new_hap1.push_back(posshap1[a]);
	    new_hap2.push_back(posshap2[b]);
	    
	    // We are already considering both explicitly	  
	    // 	  if (posshap1[a] != posshap2[b])
	    // 	    newpp[z2] *= 2;
	    
	    psum += newpp[z2];
	    z2++;	  
	  }
    }
  
  // adjust to sum to 1
  for (int z=0; z<z2; z++){
    newpp[z] /= psum;
    
  }


  // Display old and new here..

  int l = S[s];
  
  // Obtain the genotype of individual i at locus l (regional SNP s)
  // and compute sum of haplotype probabilities for each genotype
  
  ////////////////////////////////////
  // Consider each new possible phase
  
  double gh[3];
  gh[0] = gh[1] = gh[2] = 0;

  for (int z = 0; z < new_hap1.size(); z++)
    {

      // Implied genotype given this haplotype?
      // (remembering flipped allele coding)
      
      bool s1 = ! hap[ new_hap1[z]][s];
      bool s2 = ! hap[ new_hap2[z]][s];
      
      if ( s1 )
	{
	  if ( s2 ) 
	    gh[0] += newpp[z];
	  else
	    gh[1] += newpp[z];
	}
      else
	{
	  if ( s2 )
	    gh[1] += newpp[z];
	  else
	    gh[2] += newpp[z];
	}

    }
  
  
  //////////////////////////////////////////
  // Compare these to the observed genotype

  // Threshold is half the population frequency?
  // One-tenth of population frequency?

  double threshold = 0.1 * geno_freq[g];
  double impute_threshold = 0.9;

  if ( gh[g] < threshold ) 
    {
      HTEST << setw(par::pp_maxsnp) << P.locus[S[s]]->name << " "
	    << setw(par::pp_maxfid) << P.sample[i]->fid << " " 
	    << setw(par::pp_maxiid) << P.sample[i]->iid << " " 
	    << setw(6) << genotype(P, i, l) << " " 
	    << setw(8) << geno_freq[g] << " "
	    << setw(8) << gh[0] << " "
	    << setw(8) << gh[1] << " "
	    << setw(8) << gh[2] << " ";

      int ng = -9;
      for ( int j=0; j<=2; j++)
	if ( gh[j] > impute_threshold )
	  ng = j;
      
      if ( ng==0 )
	HTEST << setw(6) << P.locus[S[s]]->allele1+"/"
	  +P.locus[S[s]]->allele1 << "\n";
      else if ( ng==1 )
	HTEST << setw(6) << P.locus[S[s]]->allele1+"/"
	  +P.locus[S[s]]->allele2 << "\n";
      else if ( ng==2 )
	HTEST << setw(6) << P.locus[S[s]]->allele2+"/"
	  +P.locus[S[s]]->allele2 << "\n";
      else
	HTEST << setw(6) << par::missing_genotype+"/"
	  +par::missing_genotype << "\n";
      }


  ////////////////////////////////////////////////////////////
  // Only report dodgy looking genotypes even in verbose mode

  if ( par::proxy_full_report && gh[g] < threshold ) 
    {

      HTEST << "Individual " << P.sample[i]->fid << " " 
	    << P.sample[i]->iid << " ; locus " 
	    << P.locus[S[s]]->name << "\n";
      HTEST << "Observed genotype is " 
	    << genotype(P, i, S[s]) 
	    << " g = " << g << "\n";
      HTEST << "Prior genotypes probabilities = " 
	    << geno_freq[0] << " " 
	    << geno_freq[1] << " " 
	    << geno_freq[2] << endl;      
      HTEST << "Posterior genotypes probabilities = " 
	    << gh[0] << " " 
	    << gh[1] << " " 
	    << gh[2] << endl;      


      HTEST << "\nOriginal phases\n";
      for (int z = 0; z < hap1[i].size(); z++)
	{
	  
	  HTEST << setw(par::pp_maxfid) << P.sample[i]->fid<< " "
	       << setw(par::pp_maxiid) << P.sample[i]->iid<< " "
	       << setw(4) << z << " "
	       << setw(10) << haplotypeName(hap1[i][z]) << " "
	       << setw(10) << haplotypeName(hap2[i][z]) << "\t";      
	  
	  if (ambig[i])
	    {
	      HTEST << setw(12) << pp[i][z]<< " ";
	      int max_z = 0;
	      for (int z2=0; z2<hap1[i].size(); z2++)
		max_z = pp[i][z2] > pp[i][max_z] ? z2 : max_z ;
	      
	      if (max_z == z)
		HTEST << setw(6) << 1<< " "<< "  ";
	      else
		HTEST << setw(6) << 0<< " "<< "  ";
	    }
	  else
	    HTEST << setw(12) << 1<< " "<< setw(6) << 1<< " "<< "  ";
	  
	  // Genotypes
	  for (int j=0; j<ns; j++)
	    HTEST << genotype(P, i, S[j]) << " ";
	  HTEST << "\n";
	  HTEST << "\n";
	}
      
      HTEST << "New possible phases, if reference were missing\n";
      for (int z = 0; z < new_hap1.size(); z++)
	{
	  
	  HTEST << setw(par::pp_maxfid) << P.sample[i]->fid<< " "
	       << setw(par::pp_maxiid) << P.sample[i]->iid<< " "
	       << setw(4) << z << " "
	       << setw(10) << haplotypeName(new_hap1[z]) << " "
	       << setw(10) << haplotypeName(new_hap2[z]) << "\t";      
	  
	  HTEST << setw(12) << newpp[z]<< " ";
	  int max_z = 0;
	  for (int z2=0; z2<new_hap1.size(); z2++)
	    max_z = newpp[z2] > newpp[max_z] ? z2 : max_z ;
	  
	  if (max_z == z)
	    HTEST << setw(6) << 1<< " "<< "  ";
	  else
	    HTEST << setw(6) << 0<< " "<< "  ";
	  
	  // Genotypes
	  for (int j=0; j<ns; j++)
	    {
	      if ( j==s )
		HTEST << par::missing_genotype << "/" 
		      << par::missing_genotype << " ";
	      else
		HTEST << genotype(P, i, S[j]) << " ";
	    }
	  HTEST << "\n";
	  HTEST << "\n";
	}

      HTEST << "----------------------------------------------------------------------------\n";
      HTEST << endl;

      
    }
  
}


