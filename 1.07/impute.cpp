

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

extern ofstream LOG;
using namespace std;


class probabilisticGenotype{
public:
  probabilisticGenotype()
  {
    AA = AB = BA = BB = 0;
    phased = genotype = false;
    genotype = phased_genotype = -1;
    calculated = false;
  }

  bool calculated;
  double AA, AB, BA, BB;
  bool phased;
  bool genotyped;
  int genotype; // 0,1,2 = AA, AB, BB
  int phased_genotype; // 0,1,2,3 = AA, AB, BA, BB

};


void HaploPhase::updateForImputation()
{
  
  // Goal: given results now in HaploPhase, can we eliminate
  // any HaploWindow phases for any individuals
  
  // For each window involved, look at each phase of each
  // individual: was this supported by a phase in HaploPhase?
  
  // Also, start trying to order hap1 and hap2 across windows
  // to be consistent
  
  
  //////////////////////////////////////////////////////
  // Reconcile HaploPhase (waplotype) results back into
  // subhaplotypes of the HaploWindows
    
  int num_phase_0 = 0;
  int num_phase_1 = 0;
  
  for (int i=0; i<P.n; i++)
    {
      if ( ! include[i] ) 
	continue;
      if ( ! P.sample[i]->founder )
	continue;
      
      // If there were no available / possible haplotypes at this position, 
      // then just leave the windows as are.
      
      if ( hap1[i].size() == 0 )
	continue;
      
       
      for (int w=startWindow; w<=finishWindow; w++)
	{
	  
	  HaploWindow * currentWindow = windows[w];
		  
	  int wc = w - startWindow;

	  // Track original number of phases for this window
	  num_phase_0 +=  currentWindow->hap1[i].size();
	  
	  currentWindow->hap1[i].clear();
	  currentWindow->hap2[i].clear();
	  currentWindow->pp[i].clear();
	  
	  map<int2,double> added;
	  
	  for (int z=0; z<hap1[i].size(); z++)
	    {

	      if ( pp[i].size() == 0 ||
		   pp[i][z] > par::haplo_plem_meta_prune_phase )
		{			  
		  
		  // Add this to window, after checking to see
		  // if it already exists
		  
		  int2 subhaplotype;

		  subhaplotype.p1 = hapi[ hap1[i][z] ][wc];
		  subhaplotype.p2 = hapi[ hap2[i][z] ][wc];
		  
		  map<int2,double>::iterator ia = added.find( subhaplotype );
		  
		  if ( ia == added.end() )
		    {
		      int2 t;
		      t.p1 = hapi[ hap1[i][z] ][wc];
		      t.p2 = hapi[ hap2[i][z] ][wc];
		      
		      currentWindow->hap1[i].push_back( t.p1 );
		      currentWindow->hap2[i].push_back( t.p2 );
		      added.insert( make_pair( t, pp[i][z] ));
		    }
		  else
		    {
		      ia->second += pp[i][z];
		    }
		}
	    }
	  
	  
	  if ( currentWindow->hap1[i].size() > 1 ) 
	    {
	      currentWindow->pp[i].resize( currentWindow->hap1[i].size());
	      double psum = 0;
	      for (int z=0; z< currentWindow->hap1[i].size(); z++)
		{
		  int2 subhaplotype;
		  subhaplotype.p1 = currentWindow->hap1[i][z];
		  subhaplotype.p2 = currentWindow->hap2[i][z];
		  
		  map<int2,double>::iterator ia = added.find( subhaplotype );

		  if ( ia != added.end() ) // CAN REMOVE THIS CHECK
		    {
		      currentWindow->pp[i][z] = ia->second;
		      psum += ia->second;
		      
		    }		      
		}
	      
	      currentWindow->ambig[i] = true;
	      for (int z=0; z< currentWindow->pp[i].size(); z++)
		currentWindow->pp[i][z] /= psum;
	    }
	  else
	    {
	      currentWindow->pp[i].clear();
	      currentWindow->ambig[i] = false;
	    }

	  // Track updated number of phases for this window
	  num_phase_1 +=  currentWindow->hap1[i].size();
	  
	}
      
    }
  
  if (par::haplo_plem_verbose)
    {
      double reduction = (double) num_phase_1 / (double) num_phase_0 ;
      P.printLOG(dbl2str(reduction) 
		 + " pruning, from " 
		 + int2str(num_phase_0) 
		 + " to " + int2str( num_phase_1 ) 
		 + " phases     \n");
    }
}



void HaploPhase::mainImputation()
{
  

  //  P.printLOG("Entering final genotype imputation stage\n");
  

  //////////////////////////////////////////////////////////////////////////
  // Calculate information weights based on empirical variance for each SNP
  
  vector_t info(ns);
  vector_t infoCount(ns);

  for( int w = 0; w < nw; w++)
    {
      
      HaploWindow * currentWindow = windows[w];
      
      for (int s = 0; s < currentWindow->ns; s++)
	{
	  
	  int gs = currentWindow->start + s; 
	  
	  calculateEmpiricalVariance(gs);
	  
	  info[gs] += ratio;
	  infoCount[gs]++;	    	   
	}
    }
  

  // Normalise information score
  
  for (int s=0; s<ns; s++)
    info[s] /= infoCount[s];
  


  ///////////////////////////////////////////////////////////////////////
  // Consider each individual, imputing likely genotypes for SNPs with
  // high informativeness
  
  for (int i=0; i<P.n; i++)
    {
            
      //////////////////////////////////////////
      // Should we look at or skip this person?
      
      if ( ! include[i] ) 
	continue;
      if ( ! P.sample[i]->founder )
	continue;
      

      //////////////////////////////////////////////////////////
      // Store all imputed/phased genotypes for this individual
      
      vector<probabilisticGenotype> g(ns);
      

      ///////////////////////////////////
      // Consider each window
      
      for( int w = 0; w < nw; w++)
	{
	  

	  HaploWindow * currentWindow = windows[w];
	  

	  ////////////////////////////////
	  // Consider each possible phase
	  
	  for ( int z = 0 ; z < currentWindow->hap1[i].size(); z++)
	    {
	      

	      double posterior = currentWindow->ambig[i] ? 
		currentWindow->pp[i][z] : 1 ;
	      
	      // Consider each position
	      
	      for (int s = 0; s < currentWindow->ns; s++)
		{
		  
		  int gs = currentWindow->start + s; 
		  

		  // Do not attempt to impute low-confidence SNPs
		  
		  if ( info[ gs ] < par::proxy_info_threshold ) 
		    continue;
		  

		  // Otherwise calculate dosage
		  
		  if ( currentWindow->hap[currentWindow->hap1[i][z]][s] )
		    {
		      if ( currentWindow->hap[currentWindow->hap2[i][z]][s] )
			g[gs].AA += posterior;
		      else
			g[gs].AB += posterior;
		    }
		  else
		    {
		      if ( currentWindow->hap[currentWindow->hap2[i][z]][s] )
			g[gs].BA += posterior;
		      else
			g[gs].BB += posterior;
		      
		    }
		  
		  // Next SNP
		} 
	      
	    } // Next phase			  
	  
	} // Next window
      
      

      // Normalise dosage
      
      for (int s=0; s < g.size(); s++)
	{

	  double psum = g[s].AA + g[s].AB + g[s].BA + g[s].BB;
	  
	  if ( psum > 0 ) 
	    {
	      g[s].AA /= psum;
	      g[s].AB /= psum;
	      g[s].BA /= psum;
	      g[s].BB /= psum;		      	      
	    }
	  

	  // Impute into missing genotype data spaces; or give verbose
	  // output to a file
	  
	  if ( g[s].AA > par::proxy_impute_threshold )
	    {
	      g[s].genotype = 0;
	      g[s].phased_genotype = 0;
	      g[s].genotyped = g[s].phased = true;
	    }
	  else if ( g[s].BB > par::proxy_impute_threshold )
	    {
	      g[s].genotype = 2;
	      g[s].phased_genotype = 3;
	      g[s].genotyped = g[s].phased = true;
	    }
	  else if ( g[s].AB > par::proxy_impute_threshold )
	    {
	      g[s].genotype = 1;
	      g[s].phased_genotype = 1;
	      g[s].genotyped = g[s].phased = true;
	    }
	  else if ( g[s].BA > par::proxy_impute_threshold )
	    {
	      g[s].genotype = 1;
	      g[s].phased_genotype = 2;
	      g[s].genotyped = g[s].phased = true;
	    }		      
	  else if ( g[s].AB + g[s].BA > par::proxy_impute_threshold )
	    {
	      g[s].genotype = 1;
	      g[s].genotyped = true;
	      g[s].phased = false;
	    }
	  else
	    {
	      g[s].genotyped = false;
	      g[s].phased = false;
	    }
	  
	  
	  
	  ////////////////////////////////////
	  // Impute any missing genotype data
	  
	  bool s1 = par::SNP_major ? P.SNP[S[s]]->one[i] : P.sample[i]->one[S[s]] ;
	  bool s2 = par::SNP_major ? P.SNP[S[s]]->two[i] : P.sample[i]->two[S[s]] ;
	  
	  string original_genotype = genotype(P,i,S[s]);
		      

	  if ( s1 && ! s2 ) // Original data are missing
	    {
	      if ( g[s].genotyped ) 
		{
		  if ( par::SNP_major ) 
		    {
		      if ( g[s].genotype == 0 )
			{
			  P.SNP[S[s]]->one[i] = false;
			  P.SNP[S[s]]->two[i] = false;
			}
		      else if ( g[s].genotype == 1 )
			{
			  P.SNP[S[s]]->one[i] = false;
			  P.SNP[S[s]]->two[i] = true;
			}
		      else
			{
			  P.SNP[S[s]]->one[i] = true;
			  P.SNP[S[s]]->two[i] = true;
			}
		    }
		  else
		    {
		      if ( g[s].genotype == 0 )
			{
			  P.sample[i]->one[S[s]] = false;
			  P.sample[i]->two[S[s]] = false;
			}
		      else if ( g[s].genotype == 1 )
			{
			  P.sample[i]->one[S[s]] = false;
			  P.sample[i]->two[S[s]] = true;
			}
		      else
			{
			  P.sample[i]->one[S[s]] = true;
			  P.sample[i]->two[S[s]] = true;
			}
		    }
		}
	    }
		    


	  //////////////////////
	  // Verbose output mode
	  
	  if ( par::impute_verbose ) 
	    {

	      int l = S[s];			  		  
	      
	      HIMPUTE << P.sample[i]->fid << "\t" 
		      << P.sample[i]->iid << "\t"
		      << P.locus[l]->name << "\t"; 
	      

	      string g1 = P.locus[l]->allele1;
	      string g2 = P.locus[l]->allele2; 
	      
	      
	      // Assumption: par::proxy_impute_threshold must be at
	      // least 50%
	      
	      HIMPUTE << setw(8) << g[s].AA << " "
		      << setw(8) << g[s].AB << " "
		      << setw(8) << g[s].BA << " "
		      << setw(8) << g[s].BB << " "
		      << setw(8) << info[s] << " "
		      << setw(10) << g[s].AA + 0.5 * ( g[s].AB + g[s].BA ) << " ";
	      
	      if ( g[s].AA > par::proxy_impute_threshold )
		{
		  HIMPUTE << g1 << " " << g1 << "  "
			  << g1 << " " << g1 << "\t";
		}
	      else if ( g[s].BB > par::proxy_impute_threshold )
		{
		  HIMPUTE << g2 << " " << g2 << "  "
			  << g2 << " " << g2 << "\t";
		}
	      
	      else if ( g[s].AB > par::proxy_impute_threshold )
		{
		  HIMPUTE << g1 << " " << g2 << "  ";
		  
		  if ( g1 < g2 ) 
		    HIMPUTE << g1 << " " << g2 << "\t";
		  else
		    HIMPUTE << g2 << " " << g1 << "\t";
		}
	      
	      else if ( g[s].BA > par::proxy_impute_threshold )
		{
		  HIMPUTE << g2 << " " << g1 << "  ";
		  
		  if ( g1 < g2 ) 
		    HIMPUTE << g1 << " " << g2 << "\t";
		  else
		    HIMPUTE << g2 << " " << g1 << "\t";
		}
			  
	      else if ( g[s].AB + g[s].BA > par::proxy_impute_threshold )
		{
		  HIMPUTE << par::missing_genotype << " " 
			  << par::missing_genotype << "  ";
			      
		  if ( g1 < g2 ) 
		    HIMPUTE << g1 << " " << g2 << "\t";
		  else
		    HIMPUTE << g2 << " " << g1 << "\t";
		}
	      else
		{
		  HIMPUTE << par::missing_genotype << " " 
			  << par::missing_genotype << "  "
			  << par::missing_genotype << " " 
			  << par::missing_genotype << "\t";
		}
	      
	      HIMPUTE << original_genotype << "\n";	      	      
	      
	      // End of verbose output mode
	    }  
	  

	  // Next SNP
	} 
      

      // Next individual	  
    } 
	      
}

