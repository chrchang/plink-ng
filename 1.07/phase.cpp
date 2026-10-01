

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
#include <cassert>

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "genogroup.h"
#include "phase.h"
#include "haplowindow.h"

extern ofstream LOG;
using namespace std;


vector_t HaploPhase::phaseAllHaplotypes(bool display, Perm & perm )
{
  
  vector_t results;


  /////////////////////////////////////////////
  // Set individual to follow in verbose mode?

  if ( par::haplo_plem_follow )
    {
      par::haplo_plem_follow = false;
      for (int i=0; i<P.n; i++)
	{
	  if ( P.sample[i]->fid == par::haplo_plem_follow_fid &&
	       P.sample[i]->iid == par::haplo_plem_follow_iid )
	    {
	      P.printLOG("Following individual [ " 
			 + par::haplo_plem_follow_fid + " " 
			 + par::haplo_plem_follow_iid + " ] in EM\n");
	      par::haplo_plem_follow_ind = i;
	      par::haplo_plem_follow = true;
	    }
	}
    }


  //////////////////////////////
  // Begin phasing (and testing)
	
  int nms = new_map.size();


  for (int l=0; l<nms; l++)
    {

      if (display && ! par::silent)
	{
	  cerr << l+1<< " out of "<< nms 
	       << " haplotypes phased       \r";
	  cerr.flush();
	}

      // Set current haplotype
      current = l;

      // Get number of predictors
      int s = new_pred_locus[l].size();

      // X or haploid markers?
      X = par::chr_sex[new_map[l]->chr];
      haploid = par::chr_haploid[new_map[l]->chr];

      // Same window as previous, unless in weighted 
      // multimarker test mode, in which case do all

      bool same = true;
      if (l==0 || par::weighted_mm)
	same = false;
      else
	{
	  if (new_pred_locus[l].size() != new_pred_locus[l-1].size() )
	    {
	      same = false;
	    }
	  else
	    for (int j=0; j<s; j++)
	      {
		if (new_pred_locus[l][j] != new_pred_locus[l-1][j])
		  {
		    same = false;
		    break;
		  }
	      }
	}


  
      ///////////////////////////////////////////////////////////////
      // Calculate haplotype frequencies and posterior probabilities,
      // if different from previous round

      if (!same)
	{
	  
	  ///////////////////////////////////////////////////
	  // Initialise

	  reset();
	      
	  // Is this a wildcard haplotype?
	  string nm = new_map[l]->name;
	  if (nm.substr(nm.size()-1) == "_")
	    name(new_map[l]->name.substr(0, new_map[l]->name.find("_")));
	  else
	    name(new_map[l]->name);

	  S = new_pred_locus[l];
	  ns = S.size();


	  for (int i=0; i<P.n; i++)
	    {
	      Individual * person = P.sample[i];
	      if (person->founder)
		{
		  if (haploid || (X && person->sex))
		    validN++;
		  else
		    validN+=2;
		}

	      includeIndividuals(i);
	    }

 
	  ////////////////////////////////////////
	  // Verbose recording? (one window only)

	  if ( par::haplo_plem_verbose || par::haplo_plem_follow )
	    VPHASE.open("phased.verbose",ios::out);
	  

	  ///////////////////////////////////////////////////
	  // Define 'windows' within the 'region'

	  int pos = 0;

	  // Number of possible stub haplotypes
	  nsh = (int) pow(2.0, (double) par::haplo_plem_overlap );


	  // Check we have a valid overlap size
	  par::haplo_plem_overlap = par::haplo_plem_original_overlap;
	  if ( par::haplo_plem_window >= ns || 
	      par::haplo_plem_overlap >= par::haplo_plem_window )
	    par::haplo_plem_overlap = 0;


	  while (1)
	    {

	      HaploWindow * window = new HaploWindow( this , & P );


	      //////////////////////////////////
	      // Define SNP span for this window

	      window->start = pos;
	      window->stop = pos + par::haplo_plem_window - 1;
	      
	      if (window->stop >= ns )
		window->stop = ns-1;

	      /////////////////////////////////////////////////////////////////
	      // Enumerate haplotypes/genoGroups/phases and tally
	      // unambiguous

	      window->enumerateHaplotypes(new_pred_locus[l]);

	      window->enumerateGenogroups();

	      set<MultiLocusGenotype*>::iterator im 
		= window->genotypes.begin();

	      while (im != window->genotypes.end() )
		{
		  if (P.sample[ (*im)->reference ]->founder)
		    window->enumeratePhase( (*im)->reference );
		  ++im;
		}

	      window->tallyUnambiguousCounts();

	      ///////////////////////////////
	      // Add the window to the list

	      windows.push_back( window );
	      
	      ///////////////////////////
	      // Reached end of region?

	      if (window->stop == ns-1)
		break;


	      ///////////////////////////////////
	      // Otherwise advance by one window

	      pos += par::haplo_plem_window - par::haplo_plem_overlap;


  	      if (par::haplo_plem_verbose)
  		{
 		  cout << windows.size() << " windows added        \n";
 		  cout << window->stop << " ... " << ns-1 << "\n";
 		}
	      
	    }

  	  if ( par::haplo_plem_verbose )
  	    cout << "\n";

	  // Set number of windows;
	  nw = windows.size();

	  if (par::haplo_plem_verbose)
	    P.printLOG("Constructed "+int2str(nw)+" windows\n");


	  /////////////////////////////////////
	  // Use offspring to reduce ambiguity
	  
	  // To be added in; not yet implemented //
	  //  if (nonfounders)
	  //     for (int i=0; i<P.n; i++) 
	  //       if (P.sample[i]->founder && ambig[i])
	  // 	resolveWithChildren(i)
	  



	  ////////////////////////////////
	  // E-M phasing based on founders
	  
	  // Currently three versions of EM
	  
	  // performAlternEM()     drive chain-of-windows EM
	  // performEM()           individual-window EM
	  // performEM_original()  used for meta-EM	  
	  
	  performAlternEM();


	  //////////////////////////////////
	  // Free window storage

	  for (int w=0; w<nw; w++)
	    delete windows[w];
	  windows.clear();


	  //////////////////////////////////
	  // Prune unlikely regional phases

 	  for (int i=0; i<P.n; i++)
 	    if (P.sample[i]->founder)
 	      prunePhase(i);
	  

	  ////////////////////////////////////////
	  // Fill-in phasing for any non-founders

	  // These functions work at per-individual level, i.e. we do
	  // not use genogroups here


	  if ( nonfounders )
	    {
	      phasemap.clear();
	      phasemap.resize(P.n);

	      if (par::test_hap_TDT || par::proxy_TDT )
		{
		  trans.clear();
		  untrans.clear();
		  trans.resize(nh,0);
		  untrans.resize(nh,0);

		  transmissionX.clear();
		  transmissionX2.clear();
		  transmissionX.resize(nh,0);
		  transmissionX2.resize(nh,0);
		}
		  
	      if (haploid )
		error("Family-based haplotyping only for autosomes and X chromosome");
	      
	      // List all possible non-rare phases
	      //  enumerateAllPhases();

	      for (int i=0; i<P.n; i++)
		if (!P.sample[i]->founder)
		  {
	      
		    // Only phase affecteds, if performing TDT
		    if (par::test_hap_TDT || par::proxy_TDT )
		      if ( !P.sample[i]->aff)
			continue;

		    phaseAndScoreNonfounder(i);

		    prunePhase(i);

		  }
	    }
	  	  

	  if ( par::haplo_plem_verbose 
	       || par::haplo_plem_follow )
	    VPHASE.close();
	  
	  

	  ///////////////////////////////////////////////////////////////
	  //                                                           //
	  // Post-phasing actions                                      //
	  //                                                           //
	  ///////////////////////////////////////////////////////////////


	  ////////////////////////////////
	  // Haplotype frequencies

	  if (par::display_hap_freqs)
	    reportHaplotypeFrequencies();



	  ////////////////////////////////
	  // Haplotype association tests
	  
	  if (par::test_hap_CC ||
	      par::test_hap_GLM ||
	      par::test_hap_TDT ||
	      par::test_hap_QTL )
	    {

	      // Are we testing all haplotypes, or a specified
	      // pre-determined one?

	      if ( !par::phase_hap_all )
		setTestHaplotype(new_pred_allele[l]);
	      
	      // Perform the actual tests: ignore the results for now
	      // (these will be used when we incorporate permutation)
	      
	      vector_t t = performHaplotypeTests(display, perm );
	      for (int i=0; i<t.size(); i++)
		results.push_back( t[i] );
	      
	    }
	  


	  ////////////////////////////////
	  // Haplotype phase probabilities

	  if (par::display_phase_probs)
	    {
	      if (par::display_phase_probs_wide )
		reportPhaseWideFormat();
	      else
		reportPhase();
	    }



	  ////////////////////////////////////////////////////////
	  // Tracking haplotype sharing for a particular pair of
	  // individuals

	  if (par::segment_haplotrack)
	    {
	      trackThisSegment();
	    }
	  
	  // next unique set of SNPs to be phased 
	}
      

      
      // Impute both rare and non-rare
	
      if ( par::impute_tags ) 
	{
	  
	  // If we have specified a particular allele, impute that;
	  // otherwise, impute all above a certain frequency, and 
	  // add to map
	  
	  setTestHaplotype(new_pred_allele[l]);
	  
	  if ( test_hap > -1 ) 
	    imputeThisHaplotype(l);
	  else
	    {
	      for ( test_hap = 0; test_hap < nh; test_hap++ )
		if ( f[test_hap] >= par::min_hf 
		     && f[test_hap] <= par::max_hf )
		  {
		    imputeThisHaplotype(l);
		  }
	    }
	}



      ////////////////////////////////////
      // Next set of SNPs to be phased 
      
    }
  
  if ( par::haplo_plem_verbose )
    cout << "\n";


  return results;

}

void HaploPhase::enumerateHaplotypes(vector<int> & s)
{

  // Make list of haplotypes, and code as +/-
  
  S = s;
  ns = s.size();
  nh = (int)pow((double)2, ns);
  f.resize(nh);
  
  ph_hap1.clear();
  ph_hap2.clear();
  ph_freq.clear();

  // Optionally, for haplotype-based TDT
  if (par::test_hap_TDT || par::proxy_TDT )
    {
      trans.clear();
      untrans.clear();

      trans.resize(nh, 0);
      untrans.resize(nh, 0);
    }

  unsigned int h=0;

  while (h<nh)
    {
      vector<bool> m1;

      unsigned int p=1;
      for (int s=0; s<ns; s++)
	{
	  if (h & p )
	    m1.push_back(false);
	  else
	    m1.push_back(true);
	  p <<= 1;
	}

      // Add to list of haplotypes
      //hap.push_back(m1);

      // Add to HapMap
      //hapmap.insert(make_pair(m1, h));

      // Consider next haplotype
      h++;

    }

  // Product of allele frequencies
  for (int h=0; h<nh; h++)
    {
      f[h] = 1;
      for (int s=0; s<ns; s++)
	{
	  if ( hap[h][s] )
	    f[h] *= P.locus[S[s]]->freq;
	  else
	    f[h] *= ( 1 - P.locus[S[s]]->freq);
	}
    }
}

void HaploPhase::includeIndividuals(int i)
{

  // Do not look at non-reference individuals in some circumstances

  if ( reference_only && ! P.sample[i]->missing )
    {
      include[i] = false;

      if (P.sample[i]->founder)
	{
	  if (haploid || (X && P.sample[i]->sex))
	    validN--;
	  else
	    validN-=2;
	}
      
      return;
    }
  
  vector<bool> s1(ns);
  vector<bool> s2(ns);

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

  // But if too much missing genotype data, then 
  // we should not even try to phase this individual
  // for this region; note -- females should always be
  // missing all genotypes for Y, so we don't need to 
  // worry about allowing for a special case here.

  if ( (double)nm/(double)ns >= par::hap_missing_geno )
    {

      include[i] = false;

      if (P.sample[i]->founder)
	{
	  if (haploid || (X && P.sample[i]->sex))
	    validN--;
	  else
	    validN-=2;
	}

      return;
    }

}



void HaploPhase::performAlternEM()
{
  

  // Working variables for convergence for first 
  // pass meta-EM (chain of EMs)

  int iter = 0;
  bool converged;
  int num_converged = 0 ;

  // Start multiple EM runs, one per window, allowing each to update
  // the adjoining window's EM state

  startWindow = 0; 
  finishWindow = nw-1;

  bool verboseOutputInFirstRound = true;

  if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
    {
      VPHASE << "\n\nFOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
      VPHASE << " PRIOR TO START OF ANY EM\n";
      verboseDisplayWindows(par::haplo_plem_follow_ind);
    }
  
  
  do
    {

      if ( par::haplo_plem_verbose && verboseOutputInFirstRound )
	VPHASE << "OUTER LOOP ITERATION " << iter << "\n\n";
      
      // Iterate forward through windows
      int w = 0;

      while (w < nw )
	{

// 	  if ( ! par::silent ) 
// 	    {	      
// 	      cout << num_converged << " converged windows at iteration " 
// 		   << iter << " forwards chain, window " << w << "      \r";
// 	      cout.flush();
// 	    }

	  if ( par::haplo_plem_verbose && verboseOutputInFirstRound )
	    VPHASE << "\n\nFORWARD WINDOW " << w << "\n";
	  
	  HaploWindow * currentWindow = windows[w];
	  
	  if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
	    {
	      VPHASE << "\n\nFOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
	      VPHASE << "  FORWARD LOOP, WINDOW " << w << " PRIOR TO ADJUSTMENT \n";
	      verboseDisplayWindows(par::haplo_plem_follow_ind);
	    }

	  ///////////////////////////////////////////////////
	  // Only adjust freq if this window isn't converged

	  if ( !currentWindow->converged )
	    {
	      
	      if ( par::haplo_plem_verbose && verboseOutputInFirstRound )
		VPHASE << "WINDOW NOT YET CONVERGED\n";

	      
	      ///////////////////////////////////////////////////
	      // Only pass freq information once after convergence
	      

	      if ( false && w > 0 && !windows[w-1]->right_passed )
		{

		  if ( par::haplo_plem_verbose && verboseOutputInFirstRound )
		    VPHASE << "ADJUSTING FREQUENCIES BASED ON PREVIOUS WINDOW " << w-1 << "\n";

		  HaploWindow * previousWindow = windows[w-1];

		  // Populate stub frequencies for both windows

		  vector_t currentStub = currentWindow->leftStubFrequency();
		  vector_t previousStub = previousWindow->rightStubFrequency();

		  if ( par::haplo_plem_verbose )
		    {
		      for (int z=0; z<currentStub.size(); z++)
			{
			  VPHASE << "STUB " << z << "\t"
				 << currentStub[z] << " " 
				 << previousStub[z] << "\n";
			}		   
		    }
		  
		  vector_t of = currentWindow->f;
		  
		  // Adjust haplotype frequencies

		  for (int h = 0; h < currentWindow->nh; h++)
		    {
		      int stub = currentWindow->leftStub[h];
		      
		      if ( abs( previousStub[ stub ] - currentStub[ stub ] ) > 0.02 )
			{
			  double nf;
			  if (currentStub[ stub ] == 0)
			    nf = 0;
			  else
			    nf = currentWindow->f[h] * ( previousStub[ stub ] / currentStub[ stub ] );
			  
			  currentWindow->f[h] = ( currentWindow->f[h] + nf ) * 0.5; 
			}
		    }
		  

	      if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
		{
		  VPHASE << "FOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
		  VPHASE << "  FORWARD LOOP, WINDOW " << w << " POST ADJUSTMENT \n";
		  verboseDisplayWindows(par::haplo_plem_follow_ind);
		}

	      
	      if ( par::haplo_plem_verbose )
		{
		  
		  VPHASE << "OLD, ADJUSTED FREQS\n";
		  
		  for (int h = 0; h < currentWindow->nh; h++)
		    { 
		      if ( of[h] > 0.001 || currentWindow->f[h] > 0.001 ) 
			VPHASE << "ORIGINAL HAP " << currentWindow->haplotypeName(h) 
			       << "\t"
			       << currentWindow->leftStub[h] << " " 
			       << of[h] << "\t"
			       << currentWindow->f[h] << "\n";
		    }
		  VPHASE << "\n";
		}
	      
	      
	      // Once window has converged pass freq and keep track
	      if (previousWindow->converged)
		{
		  previousWindow->right_passed = true;
		  
		  if ( par::haplo_plem_verbose )
		    VPHASE << "\nPREVIOUS WINDOW CONVERGED AND PASSED ON: NOW FIXED\n";
		}
		}
	      
	      if ( par::haplo_plem_verbose )
		VPHASE << "\nENTERING INNER EM FOR THIS WINDOW\n";
	      
	      // Do some more EM iterations
	      currentWindow->performEM();

	      if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
		{
		  VPHASE << "FOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
		  VPHASE << "  FORWARD LOOP, WINDOW " << w << " POST EM PHASING \n";
		  verboseDisplayWindows(par::haplo_plem_follow_ind);
		}


	      if ( par::haplo_plem_verbose )
		VPHASE << "\nENTERING PRUNING FOR THIS WINDOW\n";
	      
	      // Get rid of unlikely phases 
	      currentWindow->pruneGenogroups();

	      if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
		{
		  VPHASE << "FOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
		  VPHASE << "  FORWARD LOOP, WINDOW " << w << " POST EM PRUNING \n";
		  verboseDisplayWindows(par::haplo_plem_follow_ind);
		}

	      
	      if ( par::haplo_plem_verbose )
		VPHASE << "\nDONE WITH WINDOW " << w << " IN FORWARDS LOOP\n";

	    }
	  

	  // Next window
	  ++w;
	}



      //////////////////////
      // Now move backwards

      w = windows.size() - 1;

      while (w >= 0)
	{

	  if ( par::haplo_plem_verbose ) 
 	    cout << num_converged << " converged windows at iteration " 
		 << iter << " backwards chain, window " << w << "      \r";
	  

	  if ( par::haplo_plem_verbose && verboseOutputInFirstRound )
	    VPHASE << "\n\nBACKWARD WINDOW " << w << "\n";

	  HaploWindow * currentWindow = windows[w];

	  if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
	    {
	      VPHASE << "\n\nFOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
	      VPHASE << "  BACKWARDS LOOP, WINDOW " << w << " PRIOR ADJUSTMENT \n";
	      verboseDisplayWindows(par::haplo_plem_follow_ind);
	    }


	  // Only adjust freq if window isn't converged
	  
	  if ( !currentWindow->converged)
	    {

	      if ( false && w+1< windows.size() && !windows[w+1]->left_passed)

		{
		  if ( par::haplo_plem_verbose && verboseOutputInFirstRound )
		    VPHASE << "ADJUSTING FREQUENCIES BASED ON PREVIOUS WINDOW " << w+1 << "\n";

		  HaploWindow * previousWindow = windows[w+1];

		  // Populate stub frequencies for both windows

		  vector_t currentStub = currentWindow->rightStubFrequency();
		  vector_t previousStub = previousWindow->leftStubFrequency();

		  vector_t of = currentWindow->f;

		  // Adjust haplotype frequencies

		  for (int h = 0; h < currentWindow->nh; h++)
		    {
		      int stub = currentWindow->rightStub[h];

		      if ( abs( previousStub[ stub ] - currentStub[ stub ] ) > 0.02 )
			{
			  double nf;
			  if (currentStub[ stub ] == 0)
			    nf = 0;
			  else
			    nf = currentWindow->f[h] * ( previousStub[ stub ] / currentStub[ stub ] );
			  currentWindow->f[h] = ( currentWindow->f[h] + nf ) * 0.5; 
			}
		    }

		  if ( par::haplo_plem_verbose )
		    {

		      VPHASE << "OLD, ADJUSTED FREQS\n";

		      for (int h = 0; h < currentWindow->nh; h++)
			{ 
			  if ( of[h] > 0.001 || currentWindow->f[h] > 0.001 ) 
			    VPHASE << "ORIGINAL HAP " << currentWindow->haplotypeName(h) 
				   << "\t"
				   << currentWindow->leftStub[h] << " " 
				   << of[h] << "\t"
				   << currentWindow->f[h] << "\n";
			}
		      VPHASE << "\n";
		    }

		  
		  // Once window has converged pass freq and keep
		  // track
		  
		  if (previousWindow->converged)
		    previousWindow->left_passed = true;

		}

	      if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
		{
		  VPHASE << "FOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
		  VPHASE << "  BACKWARDS LOOP, WINDOW " << w << " POST ADJUSTMENT \n";
		  verboseDisplayWindows(par::haplo_plem_follow_ind);
		}



	      // Do some more EM iterations
	      currentWindow->performEM();


	      if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
		{
		  VPHASE << "FOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
		  VPHASE << "  BACKWARDS LOOP, WINDOW " << w << " POST EM PHASING \n";
		  verboseDisplayWindows(par::haplo_plem_follow_ind);
		}

	      
	      // Get rid of unlikely phases 

	      currentWindow->pruneGenogroups();

	      if ( par::haplo_plem_follow && verboseOutputInFirstRound ) 
		{
		  VPHASE << "FOLLOWED INDIVIDUAL " << par::haplo_plem_follow_ind << "\n";
		  VPHASE << "  BACKWARDS LOOP, WINDOW " << w << " POST EM PRUNING \n";
		  verboseDisplayWindows(par::haplo_plem_follow_ind);
		}



	    }

	  // Next window
	  --w;
	}

      converged = true;
      num_converged = 0;
      for (int w = 0; w < nw; w++)
	if ( !windows[w]->converged)
	  converged = false;
	else
	  ++num_converged;
      
      // Next regional iteration
      iter++;

    } while ( !converged );

  
  if ( par::debug ) 
    cout << "Chain of EMs has convered, now putting together...                \n";
  
  if (par::haplo_plem_verbose)
    P.printLOG("Chain of E-Ms converged\n");


  ////////////////////////////////////////////////
  // Finished chain of window-EMs; now perform 
  // the meta-EM (to stitch windows together), 
  // if needed
  ////////////////////////////////////////////////


  //////////////////////////////////
  // Expand genoGroups

//   if ( ! par::silent ) 
//     cout << "Expanding genogroups...\n";

  for (int w = 0; w < nw; w++)
    {

      // A more intense pruning now at the end 

      windows[w]->pruneGenogroups( par::haplo_plem_meta_prune_phase );
      
      // before we expand out...
      windows[w]->expandGenogroups();
    }


  
  
  ////////////////////////////////////////////////////////////
  // If only a single window, just swap in the relevant variables from
  // the window
  
   if ( nw == 1 )
     {

       if ( par::meta_large_phase ) 
	 {
	   mainImputation();
	   return;
	 }

       HaploWindow * thisWindow = windows[0];
      
       f = thisWindow->f;
       hap = thisWindow->hap;
       hapmapb = thisWindow->hapmap;
       nh = thisWindow->nh;
       S = thisWindow->S;
       ns = thisWindow->ns;
       
       if ( par::test_hap_TDT )
 	{
 	  trans.clear();
 	  untrans.clear();
 	  trans.resize(nh,0);
 	  untrans.resize(nh,0);
 	}
      
       pp.resize(P.n);
       hap1.resize(P.n);
       hap2.resize(P.n);

       for (int i = 0; i < P.n; i++)
	 {	  
	   if ( include[i] ) 
	     {
	       pp[i] = thisWindow->pp[i];
	       hap1[i] = thisWindow->hap1[i];
	       hap2[i] = thisWindow->hap2[i];
	       
	       if ( hap1[i].size() == 1 )
		 ambig[i] = false;
	       else 
		 ambig[i] = true;
	       
	       if ( hap1[i].size() == 0 )
		 {
		   include[i] = false;
		 }
	       else 
		 include[i] = true;
	     }
 	}
      
       // Now we are done
       return;
      
     }
  
  
  
  ////////////////////////////////////////////////////////////
  // Combine haplotypes and define initial variables (f, pp)

  
  int waplotypeSize = par::haplo_plem_meta_window > nw ? 
    nw :
    par::haplo_plem_meta_window;
  
  if (par::haplo_plem_verbose)
    P.printLOG("Window size is " + int2str( waplotypeSize ) + " in meta-EM\n");

  startWindow = 0;
  finishWindow = startWindow + waplotypeSize - 1;
  
  if ( finishWindow >= nw ) 
    finishWindow = nw-1;

  actual_nw = finishWindow - startWindow + 1;


  if ( par::haplo_plem_follow ) 
    {
      VPHASE << "\n\n--------------\nENTERING NOW META-EM STAGE\n\n";
      if ( par::meta_large_phase ) 
	VPHASE << " ** IN IMPUTE MODE ** \n";
      else
	VPHASE << " ** IN HAPLOTYING MODE ** \n";
    }


  ////////////////////////////////////////////////////////
  // Double up phase possibilities for all windows beyond
  // the first
  
  for (int w=1; w<nw; w++)
    {
      for (int i=0; i<P.n; i++)
	{

	  if ( ! include[i] ) 
	    continue;
	  
	  if ( ! P.sample[i]->founder ) 
	    continue;
	  
	  HaploWindow * currentWindow = windows[w];
	  
	  if ( ! currentWindow->ambig[i] )
	    {
	      if ( currentWindow->hap1[i][0] != currentWindow->hap2[i][0] )
		{
		  currentWindow->hap1[i].push_back( currentWindow->hap2[i][0] );
		  currentWindow->hap2[i].push_back( currentWindow->hap1[i][0] );
		}	      
	      currentWindow->ambig[i] = true;
	      currentWindow->pp[i].resize(2,0.5);
	    }
	  else
	    {
	      int o = currentWindow->pp[i].size();
	      
	      for ( int z=0; z<o; z++)
		{
		  if ( currentWindow->hap1[i][z] != currentWindow->hap2[i][z] )
		    {
		      currentWindow->hap1[i].push_back( currentWindow->hap2[i][z] );
		      currentWindow->hap2[i].push_back( currentWindow->hap1[i][z] );
		      currentWindow->pp[i][z] /= 2;
		      currentWindow->pp[i].push_back( currentWindow->pp[i][z] );
		    }	      
		}
	    }
	}
    }

  

  //////////////////////////////////
  // Some temporary descriptions 

  long int cnt_phase = 0;
  long int cnt_hap  = 0;
  long int cnt_ambig = 0;
  for (int w = 0; w < nw; w++)
    {
      for (int i = 0; i < P.n; i++) 
	if ( include[i] && P.sample[i]->founder )
	  {
	    cnt_phase += windows[w]->hap1[i].size(); 
	    if ( windows[w]->pp[i].size() > 0 )
	      cnt_ambig++;
	  }
      cnt_hap += windows[w]->f.size();
    }

  // Average number of phases per person per window = 
  
  double avg_phase_depth = (double)cnt_phase / ( (double)P.n * (double)nw ) ;

  if (par::haplo_plem_verbose)
    {
      stringstream s2;
      s2 << "In " << nw << " windows, "
	 << cnt_phase << " total phases and " << cnt_hap << " haplotypes\n";
      s2 << "Of these, " << cnt_ambig << " individual/windows still need resolving\n";
      s2 << "Average phase depth = " << avg_phase_depth << "\n";  
      
      P.printLOG(s2.str());
    }


  ////////////////////////////////
  // Begin main loop for meta-EM

  bool completed = false;

  // For imputation mode, we might want to scan the windows multiple
  // times

  bool forwards_mode = true;
  int meta_iter = 0;
  int num_meta_iter = 1;

  while ( 1 ) 
    {
      

      /////////////////////////////////
      // Enumerate possible waplotypes

      for (int i = 0; i < P.n; i++)
	{
	  
	  // Always try to phase as much as possible when 
	  // in imputation mode
	  
	  if ( par::meta_large_phase ) 
	    include[i] = true;

	  if ( include[i] )
	    {

	      enumeratePhasedWindows(i);
	      
	      if ( hap1[i].size() > 1 )
		pp[i].resize(hap1[i].size());
	      else
		pp[i].resize(0);
	     	      
	      if ( par::haplo_plem_follow && par::haplo_plem_follow_ind == i )
		{
		  VPHASE << "AFTER ENUMERATED_PHASED DUMP, WINDOWS " << startWindow 
			 << " to " << finishWindow << "\n";
		  verboseDisplayWindows(i,false);
		  VPHASE << "\n\n";
		  
		}
	    }
	}
      

      // Generate starting values for haplotypes -- for now use a
      // uniform distribution
      
      if ( nh != hapmap.size() )
	error("weirdness...\n");
      
      if (par::haplo_plem_verbose)
	P.printLOG("Considering "+int2str(nh)
		   + " waplotypes in second-stage EM, "
		   + int2str(startWindow)
		   +" to "+int2str(finishWindow)+"\n");
      
      if ( par::haplo_plem_verbose )
	VPHASE << "ENTERING PLEM WITH " << nh << " WAPLOTYPES\n";
      
      // Set up haplotype space for HaploPhase
      
      f.resize( hapmap.size() );
      for (int h=0; h<nh; h++)
	f[h] = 1/(double)nh;
      

      // (Remove these in final version, until end?) Make overall
      // haplotype codes, etc

      S.clear();
      for(int w= startWindow; w <= finishWindow; w++)
        {
          HaploWindow * currentWindow = windows[w];
          int start = par::haplo_plem_overlap;
          if (w==startWindow) start = 0;
          for (int s = start; s < currentWindow->ns; s++)
            S.push_back( currentWindow->S[s] );
        }
      
      ns = S.size();

      if ( par::test_hap_TDT )
	{
	  trans.clear();
	  untrans.clear();
	  trans.resize(nh,0);
	  untrans.resize(nh,0);
	}
      
      
      
      if ( par::haplo_plem_follow ) 
	{

	  VPHASE << "inclusion status = " << include[ par::haplo_plem_follow_ind ] << "\n";

	  VPHASE << "HAPLOPHASE PRE_MERGING STATUS \n";
	  
	  for (int z=0; z<hap1[par::haplo_plem_follow_ind].size(); z++)
	    {
	      VPHASE << haplotypeName(hap1[par::haplo_plem_follow_ind][z]) << " / " 
		     << haplotypeName(hap2[par::haplo_plem_follow_ind][z]) << "\n";
	    }
	  
	  VPHASE << "\n\n";
	}
      

      //////////////////////////////////////////////////
      // Invoke basic EM algorithm on subset of windows 
      
//       if ( ! par::silent ) 
// 	{
// 	  cout << "Performing meta-EM on windows " << startWindow 
// 	       << " to " << finishWindow << " of " << nw << "         \r";
// 	  cout.flush();
// 	}
    
      performEM_original();
      

      if ( par::haplo_plem_verbose ) 
	{
	  
	  VPHASE << "Post meta-EM frequencies\n";
	  for (int h=0; h<nh; h++)
	    VPHASE << h << "\t" 
		   << f[h] << "\t"
		   << haplotypeName(h) << "\n";
	  
          VPHASE << "Post meta-EM frequencies of common haplotypes (1 MAF)\n";
          for (int h=0; h<nh; h++)
            if ( f[h] >= 0.01 ) 
	      VPHASE << h << "\t"
		     << f[h] << "\t"
		     << haplotypeName(h) << "\n";
	  
	  VPHASE << "Post meta-EM phase frequencies\n";
	  for (int i=0; i<P.n; i++)
	    {
	      for (int z=0; z<hap1[i].size(); z++)
		{
		  VPHASE << i << "\t"
			 << haplotypeName(hap1[i][z]) << " / " 
			 << haplotypeName(hap2[i][z]) << "\t";
		  if ( pp[i].size() == 0 ) 
		    VPHASE << "F" << ambig[i] << "\n";
		  else
		    VPHASE << pp[i][z] << " " << ambig[i] << "\n";
		}
	    }
	  VPHASE << "\n";
	}
      

      
      if ( par::haplo_plem_follow ) 
	{
	  VPHASE << "HAPLOPHASE MERGING RESULTS \n";
	  for (int z=0; z<hap1[par::haplo_plem_follow_ind].size(); z++)
	    {
	      if ( pp[par::haplo_plem_follow_ind].size() == 0 || 
		   pp[par::haplo_plem_follow_ind][z] > 0.01 )
		{
		  VPHASE << haplotypeName(hap1[par::haplo_plem_follow_ind][z]) << " / " 
			 << haplotypeName(hap2[par::haplo_plem_follow_ind][z]) << "\t";
		  VPHASE << "(freqs = " << f[hap1[par::haplo_plem_follow_ind][z]] << " / "
			 << f[hap2[par::haplo_plem_follow_ind][z]] << " )\t";		  
		  if ( pp[par::haplo_plem_follow_ind].size() == 0 ) 
		    VPHASE << "F" << ambig[par::haplo_plem_follow_ind] << "\n";
		  else
		    VPHASE << pp[par::haplo_plem_follow_ind][z] << " " 
			   << ambig[par::haplo_plem_follow_ind] << "\n";
		}
	    }
	  VPHASE << "\n\n";
	}
  


      ///////////////////////////////////////////////////////////////////
      // Given new results in phase, copy these back into the original
      // windows, applying appropriate pruning of phases
      

      if ( par::meta_large_phase ) 
	{
	  updateForImputation();
	}
      

      // In imputation mode, we mnight want to re-prune the windows
      // several times
      
      if ( par::meta_large_phase )
	{
	  
	  if ( finishWindow == nw-1 && forwards_mode ) 
	    {
	      if ( ++meta_iter == num_meta_iter )
		completed = true;
	      else
		forwards_mode = false;
	    }
	  
	  if ( startWindow == 0 && ! forwards_mode ) 
	    {
	      if ( ++meta_iter == num_meta_iter )
		completed = true;
	      else
		forwards_mode = true;
	    }
      
	}


      // In haplotyping mode, we only do a single run through
      // the windows

      if ( ! par::meta_large_phase ) 
	{
	  if ( finishWindow == nw-1 )
	    completed = true;
	}

      /////////////////////////////////////////////////////////////////
      // Have we finished? In this case, generate imputed solution, or
      // copy final solution into HaploPhase (if not already there?)
      
      if ( completed )
	{
	  
// 	  if ( !par::silent ) 
// 	    cout << "Finished meta-EM: creating the final solution...        \n";

	  ////////////////////////////////////////////////////////////
	  // Copy back SNPs for HaploPhase -- we shouldn't need to do
	  // this -- we only ever swapped things for debug purposes,
	  // therefore we can remove above and this change...
	  
	  
	  ////////////////////////////////////
	  // Imputation, or haplotyping mode?
	  
	  if ( par::meta_large_phase ) 
	    {

	      S.clear();
	      for( int w= 0; w < nw; w++)
		{
		  HaploWindow * currentWindow = windows[w];
		  int start = par::haplo_plem_overlap;
		  if (w==0) start = 0;
		  for (int s = start; s < currentWindow->ns; s++)
		    S.push_back( currentWindow->S[s] );
		}
	      
	      ns = S.size();
	      
	      mainImputation();
	    }
	  

	  /////////////////////////////////////////////////////////////
	  // Finish up haplotyping (copy frequencies, phases, etc back
	  // into HaploPhase)

	  if ( ! par::meta_large_phase )
	    {
	      
	      // The things we need / should check are all okay are:
	      
	      // f = thisWindow->f;
	      // hap = thisWindow->hap;
	      // nh = thisWindow->nh;
	      // S = thisWindow->S;
	      // ns = thisWindow->ns;
	      
	      // pp[i]
	      // hap1[i]
	      // hap2[i]
	      // ambig[i]
	      // include[i]
	      
	      ////////////////////////
	      // Some final things...
	      
	      if ( par::test_hap_TDT )
		{
		  trans.clear();
		  untrans.clear();
		  trans.resize(nh,0);
		  untrans.resize(nh,0);
		}
	      
	    }


	  //////////////////////////////////
	  // Some temporary descriptions 
	  
	  long int cnt_phase = 0;
	  long int cnt_hap  = 0;
	  long int cnt_ambig = 0;

	  for (int w = 0; w < nw; w++)
	    {
	      for (int i = 0; i < P.n; i++) 
		if ( include[i] && P.sample[i]->founder )
		  {
		    cnt_phase += windows[w]->hap1[i].size(); 
		    if ( windows[w]->pp[i].size() > 0 )
		      cnt_ambig++;
		  }
	      cnt_hap += windows[w]->f.size();
	    }
	  
	  if (par::haplo_plem_verbose)
	    {
	      stringstream s2;	  
	      s2 << "After meta-EM, in " << nw << " windows, "
		 << cnt_phase << " total phases and " << cnt_hap << " haplotypes\n";
	      s2 << "Of these, " << cnt_ambig << " individual/windows still need resolving\n";
	      P.printLOG(s2.str());
	    }
	  
	  ////////////////////
	  // And now finish

	  break;

	}

      
      //////////////////////////////////////////////////////////////////
      // Alternatively, in haplotyping mode, we want to keep track of
      // all non-rare haplotypes and accumulate

      //////////////////////////////////////////////////////////////
      // If haplotyping (estimating haplotype frequencies) as opposed
      // to imputation, then we use an alternate strategy: we copy
      // HaploPhase results back into last window in this set, and
      // shift window forwards (copying only common haplotypes
      // (default=?) and phases (default=10%) )
      
      if ( ! par::meta_large_phase ) 
	{
      
	  HaploWindow * thisWindow = windows[ finishWindow ];
	  
	  thisWindow->f.clear();
	  thisWindow->hap.clear();
	  thisWindow->nh = 0;

	  
	  /////////////////////////////////
	  // Save only non-rare haplotypes
	  
	  map<int,int> keptHaplotypes;
	  keptHaplotypes.clear();
	  int new_h = 0;
	  
	  for (int h=0; h<nh; h++)
	    {
	      if ( f[h] > par::haplo_plem_meta_prune_haplotype ) 
		{
		  thisWindow->f.push_back(f[h]);
		  thisWindow->hap.push_back(hap[h]);
		  thisWindow->nh++;
		  keptHaplotypes.insert( make_pair(h,new_h++) );
		}
	    }
	  
	  
	  /////////////////////////////////////////////////////////////////
	  // Update list of SNPs this HaploWindow now points to (handling
	  // any possible overlap)
	  
	  vector<int> tmp;
	  
	  for(int w = startWindow; w <= finishWindow; w++ )
	    {	  
	      HaploWindow * currentWindow = windows[w];
	      int start = par::haplo_plem_overlap;
	      if (w==startWindow) start = 0;
	      for (int s = start; s < currentWindow->ns; s++)
		tmp.push_back( currentWindow->S[s] );
	    }
	  
	  thisWindow->S = tmp;
	  
	  if ( hap[0].size() > 0 )
	    thisWindow->ns = hap[0].size();
	  else
	    thisWindow->ns = 0;
	  

	  ////////////////////////////////////////////////////
	  // For this new window, set the stub codes for each 
	  // new haplotype
	  
	  thisWindow->setStubCodes();
	  

	  /////////////////////////////////////////////////////////
	  // Copy per-person information (phases and probabilities)
	  
	  for (int i = 0; i < P.n; i++)
	    {
	      
	      thisWindow->pp[i].clear();
	      thisWindow->hap1[i].clear();
	      thisWindow->hap2[i].clear();
	      
	      if ( pp[i].size() == 0 ) 
		{
		  
		  thisWindow->pp[i] = pp[i];
		  if ( hap1[i].size() == 1 ) 
		    {
		      map<int,int>::iterator k1 = keptHaplotypes.find( hap1[i][0] );
		      map<int,int>::iterator k2 = keptHaplotypes.find( hap2[i][0] );
		      
		      if ( k1 != keptHaplotypes.end() 
			   &&
			   k2 != keptHaplotypes.end() )
			{
			  thisWindow->hap1[i].push_back( k1->second );
			  thisWindow->hap2[i].push_back( k2->second );      
			}
		    }
		  
		}
	      else
		{
		  for ( int z = 0 ; z < pp[i].size() ; z++ )
		    {
		      if (pp[i][z] > par::haplo_plem_meta_prune_phase )
			{
			  
			  map<int,int>::iterator k1 = keptHaplotypes.find( hap1[i][z] );
			  map<int,int>::iterator k2 = keptHaplotypes.find( hap2[i][z] );
			  
			  if ( k1 != keptHaplotypes.end() 
			       &&
			       k2 != keptHaplotypes.end() )
			    {
			      
			      thisWindow->pp[i].push_back(pp[i][z]);
			      thisWindow->hap1[i].push_back(k1->second);
			      thisWindow->hap2[i].push_back(k2->second);				      
			      
			    }
			}
		      
		    }
		  
		  // No need to rescale these, as we do not use posterior
		  // probs (or haplotype freqs) at this stage in any case
		  // yet -- but otherwise, we will want to adjust to sum
		  // to 1 if we edit out rare phases, or set ambig codes,
		  // etc
		  		  
		}	      
	    } // Next individual     
	  
	}



      if ( par::haplo_plem_follow ) 
	{
	  VPHASE << "POST WINDOW STITCHING... " 
		 << startWindow << " to " 
		 << finishWindow << "\n";
	  verboseDisplayWindows(par::haplo_plem_follow_ind,false);
	  VPHASE << "\n--------------------------------------\n";
	}
	
      
      /////////////////////////////////////////
      // Reset key variables in HaploPhase
      
      nh = 0;
      hap.clear();
      hapi.clear();
      f.clear();
      hapmap.clear();
      
      
      
      /////////////////////////////////////////
      // Advance to next set of windows       
      
      // Imputation?
      
      if ( par::meta_large_phase )
	{
	  if ( forwards_mode )
	    {  
	      ++startWindow;
	      ++finishWindow;
	    }
	  else
	    {
	      --startWindow;
	      --finishWindow;
	    }
	}

      
      // or haplotyping?

      else
	{
	  startWindow = finishWindow;
	  finishWindow = startWindow + waplotypeSize - 1;
	  if ( finishWindow >= nw ) 
	    finishWindow = nw-1;
	  actual_nw = finishWindow - startWindow + 1;
	}
      
    
    } // loop while not finished
  


  //////////////////////////////////
  // End of alternate EM algorithm

  return;
  
}




vector_t HaploPhase::performHaplotypeTests(bool display, Perm & perm)
{
  
  vector_t statistic;
  
  if ( par::test_hap_GLM )
    {                
      return P.glmHaplotypeTest(display,perm);
    }

  // Weighted multimarker test (test a single Haplotype)
  
  if ( par::weighted_mm )
    {

      if (par::test_hap_CC)
	haplotypicWeightedCC();
      else if (par::test_hap_TDT)
	haplotypicWeightedTDT();

      return statistic;
    }

  // Standard haplotype tests

  if ( !par::phase_hap_all )
    {

      map<int,int> tests;

      // Of the specific, prespecified haplotype? 
      // Test only the specific haplotype

      if (f[test_hap] >= par::min_hf)
	{
	  
	  for (int h2=0; h2 < nh; h2++)
	    {
	      if (f[h2] >= par::min_hf)
		{
		  if (test_hap==h2)
		    tests.insert(make_pair(h2, 0));
		  else
		    tests.insert(make_pair(h2, 1));
		}
	    }
	  
	  if (par::test_hap_CC)
	    haplotypicCC(tests, 2, true);
	  else if (par::test_hap_QTL)
	    haplotypicQTL(tests, 2, true);
	  else if (par::test_hap_TDT)
	    haplotypicTDT(tests, 2, true);
	
	}

    }
  else
    {
      // Perform an omnibus test and all haplotype-specific
      // tests, if the --hap-all flag has been set

      // Omnibus test

      // tests[0] = 0
      // tests[1] = 1
      // tests[2] = 2
      // ...
      // tests[h] = h

      map<int,int> tests;
      int nch=0;
      for (int h=0; h < nh; h++)
	if (f[h] >= par::min_hf)
	  {
	    tests.insert(make_pair(h, nch++));
	  }
      
      if (nch>2)
	{
	  if (par::test_hap_CC)
	    haplotypicCC(tests, nch, true);
	  else if (par::test_hap_QTL)
	    haplotypicQTL(tests, nch, true);
	  else if (par::test_hap_TDT)
	    haplotypicTDT(tests, nch, true);
	}
      
      
      // Haplotype-specific test

      // tests[0] = 0
      // tests[1] = 1,2,...,h

      // tests[0] = 1
      // tests[1] = 0,2,3,...,h

      // etc

      for (int h=0; h < nh; h++)
	{
	  
	  if (f[h] >= par::min_hf)
	    {
	      tests.clear();

	      for (int h2=0; h2 < nh; h2++)
		{
		  if (f[h2] >= par::min_hf)
		    {
		      if (h==h2)
			{
			  tests.insert(make_pair(h2, 0));
			}
		      else
			tests.insert(make_pair(h2, 1));
		    }
		}
	      
	      if (par::test_hap_CC)
		haplotypicCC(tests, 2, true);
	      else if (par::test_hap_QTL)
		haplotypicQTL(tests, 2, true);
	      else if (par::test_hap_TDT)
		haplotypicTDT(tests, 2, true);
	      
	    }
	}

    } // end of --hap-all routine


  // Dummy return vector (for now)
  return statistic;

}

void HaploPhase::prunePhase(int i)
{

  // Prune regional phases (HaploPhase)

  if ( (!include[i]) ||(!ambig[i]))
    return;


  double psum = 0;

  vector<double> new_pp(0);
  vector<int> new_h1(0);
  vector<int> new_h2(0);
  
  for (int z=0; z < hap1[i].size(); z++)
    {
      if (pp[i][z] >= par::hap_min_phase_prob)
	{
	  new_pp.push_back(pp[i][z]);
	  psum += pp[i][z];
	  new_h1.push_back(hap1[i][z]);
	  new_h2.push_back(hap2[i][z]);
	}
    }

  // Normalise?
  if (pp[i].size() > new_pp.size() )
    {
      for (int z=0; z < new_pp.size(); z++)
	new_pp[z] /= psum;
    }
  
  // Update
  pp[i] = new_pp;
  hap1[i] = new_h1;
  hap2[i] = new_h2;

}

set<int> HaploPhase::makeSetFromMap(map<int,int> & h)
{

  // Wrapper function for dosage that takes a set<int>
  // Extract indicated haplotypes (coded 0)

  set<int> tests;
  
  map<int,int>::iterator ih = h.begin();
  while ( ih != h.end() )
    {
      if ( ih->second == 0 ) 
	tests.insert( ih->first );
      ++ih;
    }

  return tests;
}


double HaploPhase::dosage(int i, set<int> & h)
{

  // Assume i and h are valid

  double d = 0;
  
  if ( ambig[i] )
    {

      vector_t & posterior = pp[i];
      vector<int> & h1 = hap1[i];
      vector<int> & h2 = hap2[i];

      for (int z = 0; z < h1.size(); z++)
	{
	  if ( h.find( h1[z] ) != h.end() ) 
	    d += posterior[z];
	  if ( ! (haploid || (X && P.sample[i]->sex)))
	    if ( h.find( h2[z] ) != h.end() ) 
	      d += posterior[z];
	}
    }
  else
    {
      if ( h.find( hap1[i][0] ) != h.end() ) 
	++d;
      if ( ! (haploid || (X && P.sample[i]->sex)))
	if ( h.find( hap2[i][0] ) != h.end() )  
	  ++d;
    }
    
  return d;

}
	  

