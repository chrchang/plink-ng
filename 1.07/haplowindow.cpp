
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
#include "stats.h"

#include "phase.h"
#include "genogroup.h"

#include "haplowindow.h"

class HaploPhase;


void verboseDisplayWindows2(HaploPhase * haplo, int i, bool use_ref = true )
{
  
  for (int w = haplo->startWindow; 
       w <= haplo->finishWindow ; w++)
    {
      
      int r = haplo->windows[w]->genoGroup[i]->reference;
      if ( ! use_ref ) 
	r = i;

      haplo->VPHASE << "WINDOW " << w << "\n";

      for (int z = 0; z < haplo->windows[w]->hap1[r].size(); z++)
	{
	  HaploWindow * thisWindow = haplo->windows[w];
	  
	  haplo->VPHASE << setw(w) << " "
			<< thisWindow->haplotypeName(thisWindow->hap1[r][z])
			<< "/"
			<< thisWindow->haplotypeName(thisWindow->hap2[r][z])
			<< " ";
	  haplo->VPHASE << "( " 
			<< thisWindow->f[ thisWindow->hap1[r][z] ] << " / " 
			<< thisWindow->f[ thisWindow->hap2[r][z] ] << " ) ";
	  
	  if ( thisWindow->hap1[r].size() == 1)
	    haplo->VPHASE << "[1]\n";
	  else
	    haplo->VPHASE << thisWindow->pp[r][z]<< "\n";
	  
	}
    }
}


string HaploWindow::haplotypeName(int h)
{
  string str;
  for (int s=0; s<ns; s++)
    {
      string a1 = haplo->P.locus[S[s]]->allele1;
      string a2 = haplo->P.locus[S[s]]->allele2;
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



HaploWindow::HaploWindow(HaploPhase * hp, Plink * plinkp)
{
  haplo = hp;
  P = plinkp;
  converged = false;
  left_passed = false;
  right_passed = false;
  genoGroup.resize(P->n, (MultiLocusGenotype*)0);
  
  ambig.resize(P->n, false);
  pp.resize(P->n);
  hap1.resize(P->n);
  hap2.resize(P->n);

}

HaploWindow::~HaploWindow()
{
  set<MultiLocusGenotype*>::iterator im = genotypes.begin();
  
  while (im != genotypes.end() )
    {
      delete *im;
      ++im;
    }
  
}

void HaploWindow::setStubCodes()
{

  /////////////////////////////////////////////////////////////
  // Set stub codes (if no overlap, all codes will be 0, okay)

  leftStub.clear();
  rightStub.clear();

  for (int h=0; h<nh; h++)
    {
      
      unsigned int leftCode = 0, rightCode = 0;
      
      unsigned int p=1;
      
      for (int s=0; s < par::haplo_plem_overlap; s++)
	{
	  if (hap[h][s])
	    leftCode += p;
	  p <<= 1;
	}
      
      p=1;
      for (int s=ns-par::haplo_plem_overlap; s < ns; s++)
	{
	  if (hap[h][s])
	    rightCode += p;
	  p <<= 1;
	}
  
      leftStub.push_back( leftCode );
      rightStub.push_back( rightCode );
    }

}

void HaploWindow::enumerateHaplotypes(intvec_t & Sall)
{

  // Make list of possible haplotypes for a window

  ns = stop - start + 1;

  S.clear();
  for (int s = start; s <= stop; s++)
    S.push_back(Sall[s]);
  

  // Number of possible haplotypes
  nh = (int)pow((double)2, ns);

  if ( nh < 0 ) 
    error("EM window size too large\n");

  f.resize(nh);
  zero.resize(nh, false);
  
  // Determine set of haplotypes

  unsigned int h=0;
  double fsum = 0;
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
      hap.push_back(m1);

      // Add to HapMap
      hapmap.insert(make_pair(m1, h));


      //////////////////////////////////
      // Product of allele frequencies
		
      f[h] = 1;
		
      for (int s=0; s<ns; s++)
	{
	  if ( hap[h][s] )
	    {
	      f[h] *= P->locus[S[s]]->freq;
	    }
	  else
	    {
	      f[h] *= ( 1 - P->locus[S[s]]->freq);
	    }
		    
	}
		
      fsum += f[h];
      // Consider next haplotype
      h++;

    }

  //////////////////
  // Set stub codes
  
  setStubCodes();

	
}

void HaploWindow::enumeratePhase(int i)
{

  vector<bool> s1(ns);
  vector<bool> s2(ns);
  
  hap1[i].clear();
  hap2[i].clear();

  // Flipping allele-coding for homozygotes

  for (int s = 0; s < ns; s++)
    {
      if (par::SNP_major)
	{
	  s1[s] = P->SNP[S[s]]->one[i];
	  s2[s] = P->SNP[S[s]]->two[i];
	}
      else
	{
	  s1[s] = P->sample[i]->one[S[s]];
	  s2[s] = P->sample[i]->two[S[s]];
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
  for (int s = 0; s < ns; s++)
    if (s1[s] && !s2[s])
      nm++;

  // If any missing genotypes, this person counts 
  // as ambiguous

  if (nm>0)
    ambig[i] = true;

  // We only worry about too much missing data at the HaploPhase
  // stage, not here


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

  //////////////////////////////////////
  // Construct list of consistent phases
  if (!ambig[i])
    {
      // Unambiguous means all no missing genotypes
      // and less than 2 hets

      // Match haplotype alleles: haploid individuals
      // will just be coded as homozygous here (but 
      // when considering phases, frequencies, etc
      // we will take care of this downstream)

      hap1[i].push_back(hapmap.find(s1)->second);
      hap2[i].push_back(hapmap.find(s2)->second);

    }
  else
    {

      // For individuals with ambiguity 

      // Which are the ambiguous sites
      // (missing or heterozygous)

      // We will not observe any hets for haploid 
      // individuals: but we do need to make sure 
      // that missing haploid genotypes are not 
      // allowed to be heterozygous

      vector<bool> het_site(ns, false);
      vector<bool> mis_site(ns, false);

      int firstHeterozygote = ns;

      int ambig_cnt=0;
      for (int s=0; s<ns; s++)
	{
	  // het
	  if ( (!s1[s]) && s2[s])
	    {
	      het_site[s] = true;

	      if (firstHeterozygote == ns )
		firstHeterozygote = s;

	      ambig_cnt++;

	    }
	  // missing
	  if (s1[s] && (!s2[s]))
	    {

	      mis_site[s] = true;

	      // haploid:  0- or 1-
	      // diplod :  00 or 01 or 11

	      if (haplo->haploid || (haplo->X && P->sample[i]->sex))
		ambig_cnt++;
	      else
		ambig_cnt+=2;

	    }
	}

      int ambig_nh = (int)pow((double)2, ambig_cnt);

      vector<bool> h1(ns);
      vector<bool> h2(ns);

      int original_firstHeterozygote = firstHeterozygote;

      int h=0;
      while (h<ambig_nh)
	{

	  vector<bool> m1;

	  unsigned int p=1;
	  for (int s=0; s<ambig_cnt; s++)
	    {
	      if (h & p )
		m1.push_back(false);
	      else
		m1.push_back(true);
	      p <<= 1;
	    }

	  // Splice m1-variant into h1, h2 to 
	  // reconstruct ambiguous sites

	  int ac=0;
	  bool skip = false;

	  // If missing homozygote imputed, the next 
	  // het can be fixed possibly; also, if not 
	  // already done, the first missing genotype
	  // can be fixed

	  int firstHeterozygote = original_firstHeterozygote ;

	  for (int s=0; s<ns; s++)
	    {

	      // Reconstruct heterozygous sites
	      // (except the first one)

	      if (het_site[s])
		{

		  if (m1[ac])
		    {
		      h1[s] = true;
		      h2[s] = false;
		    }
		  else
		    {
		      if (s != firstHeterozygote )
			{
			  h1[s] = false;
			  h2[s] = true;
			}
		      else
			{
			  h1[s] = false;
			  h2[s] = true;

			  skip = true;
			}
		    }

		  ac++;

		}

	      else if (mis_site[s])
		{
		  // Treat haploid and diploid differently for missing
		  // genotype data

		  if (haplo->haploid || (haplo->X && P->sample[i]->sex))
		    {
		      // select a hemizygote/homozygote
		      if (m1[ac])
			h1[s] = h2[s] = false;
		      else
			h1[s] = h2[s] = true;

		      ac++;
		    }
		  else
		    {
		      // Make het
		      if (m1[ac])
			{

			  if (m1[ac+1])
			    {
			      h1[s] = false;
			      h2[s] = true;
			    }
			  else
			    {
			      if (s < firstHeterozygote )
				{
				  skip = true;
				}
			      else
				{
				  h1[s] = true;
				  h2[s] = false;
				}
			    }

			  firstHeterozygote = s;

			}
		      else
			{

			  // otherwise, select a homozygote
			  if (m1[ac+1])
			    h1[s] = h2[s] = false;
			  else
			    h1[s] = h2[s] = true;
			}

		      ac+=2;

		    }
		}
	      else
		{
		  // Maintain unambigous site
		  // (which might be 1st het)
		  h1[s] = s1[s];
		  h2[s] = s2[s];
		}
	    }

	  if ( !skip )
	    {

	      // Add to (non-redundant) list?

	      int n1 = hapmap.find( h1 )->second;
	      int n2 = hapmap.find( h2 )->second;

	      hap1[i].push_back(n1 );
	      hap2[i].push_back(n2 );

	    }

	  // Consider next haplotype pair
	  h++;

	}
    }

  // Make space for posterior probabilities, and skip codes
  if (ambig[i])
    pp[i].resize(hap1[i].size());
  
  genoGroup[i]->skip.resize(hap1[i].size() , false);


}

void HaploWindow::pruneGenogroups(double t )
{
  set<MultiLocusGenotype*>::iterator im = genotypes.begin();
  while (im != genotypes.end() )
    {
      int i = (*im)->reference;
      if (P->sample[i]->founder)
	prunePhase(i,t);
      ++im;
    }
}

void HaploWindow::prunePhase(int i, double t)
{

  // pp[i][z]
  // hap1[i][z]
  // hap2[i][z]
  // skip[z] 

  if ( (!haplo->include[i]) ||(!ambig[i]))
    return;

  MultiLocusGenotype * mlg = genoGroup[i];

  double psum = 0;

  vector<double> new_pp(0);
  vector<int> new_h1(0);
  vector<int> new_h2(0);
  vector<bool> new_skip(0);

  for (int z=0; z < hap1[i].size(); z++)
    {

      if ( pp[i][z] >= t ) 
	{
	  new_pp.push_back(pp[i][z]);
	  psum += pp[i][z];
	  new_h1.push_back(hap1[i][z]);
	  new_h2.push_back(hap2[i][z]);
	  new_skip.push_back(mlg->skip[z]);
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
  mlg->skip = new_skip;
	

}

vector_t HaploWindow::leftStubFrequency()
{
  vector_t freq(haplo->nsh, 0);

  for (int h=0; h< nh; h++)
    freq[ leftStub[h] ] += f[h];
    
  return freq;
}

vector_t HaploWindow::rightStubFrequency()
{
  vector_t freq(haplo->nsh, 0);
	
  for (int h=0; h< nh; h++)
    freq[ rightStub[h] ] += f[h];
	
  return freq;
}

void HaploWindow::tallyUnambiguousCounts()
{

  uc.resize(nh, 0);
  
  set<MultiLocusGenotype*>::iterator im = genotypes.begin();

  while (im != genotypes.end() )
    {

      int i = (*im)->reference;

      if (!ambig[i])
	{
	  uc[hap1[i][0]] += (*im)->count;
	  if ( ! (haplo->haploid || (haplo->X && P->sample[i]->sex)))
	    uc[hap2[i][0]] += (*im)->count;
	}

      ++im;
    }

}

void HaploWindow::expandGenogroups()
{
  for (int i=0; i<P->n; i++)
    {
      
      if ( ! ( P->sample[i]->founder && haplo->include[i] ) )
	continue;
      
      int r = genoGroup[i]->reference;
      
      if (r != i )
	{
	  pp[i] = pp[r];
	  hap1[i] = hap1[r];
	  hap2[i] = hap2[r];
	}
    }
}

void HaploWindow::performEM()
{


  //////////////////  
  // Begin E-M

  if ( par::haplo_plem_verbose )
    haplo->VPHASE << "\nWINDOW spanning " << start << " to " << stop << "\n"
		  << "\nINNER EM LOOP FOR " << par::haplo_plem_iter << " ITERATIONS ";
    
  for (int j=0; j<=par::haplo_plem_iter; j++)
    {
      

      //////////////////////////
      // E-step for genoGroups
		
      set<MultiLocusGenotype*>::iterator im = genotypes.begin();

      while (im != genotypes.end() )
	{

	  int i = (*im)->reference;
	  
	  if (ambig[i])
	    {
	      double s=0;
	      
	      // Haploid phases...
	      if (haplo->haploid || (haplo->X && P->sample[i]->sex))
		{
		  for (int z=0; z<hap1[i].size(); z++)
		    {
		      pp[i][z] = f[hap1[i][z]];
		      s += pp[i][z];
		    }
		}
	      else // ... or diploid
		{
		  for (int z=0; z<hap1[i].size(); z++)
		    {

		      if ((*im)->skip[z])
			continue;

		      int h1 = hap1[i][z];
		      int h2 = hap2[i][z];
		      
		      if (zero[h1] || zero[h2])
			{
			  (*im)->skip[z] = true;
			  continue;
			}
		      
		      pp[i][z] = f[h1] * f[h2];
		      if (h1 != h2)
			pp[i][z] *= 2;
		      
		      s += pp[i][z];
		      
		    }
		}
	      

	      /////////////////////////////////////////
	      // Check for single phase with 0 probability

	      if ( s == 0 ) 
		{
		  if ( pp[i].size()==1 )
		    {
		      pp[i][0] = s = 1;
		      
		      if ( par::haplo_plem_verbose )
			haplo->VPHASE << "\n*** WARNING *** FIXED INDIVIDUAL "
			       << P->sample[i]->fid << " " 
			       << P->sample[i]->iid << " TO PP=1 FOR SINGLE IMPOSS PHASE\n";
		    }
		  else
		    {
		      if ( par::haplo_plem_verbose )
			{
			  haplo->VPHASE << "\n*** ERROR *** INDIVIDUAL "
					<< P->sample[i]->fid << " " 
					<< P->sample[i]->iid << " HAS >1 PHASE BUT PP SUMS TO 0\n";
			  
			  verboseDisplayWindows2(haplo,i,true);
			  
			  haplo->VPHASE.close();

			  error("See phased.verbose (--em-verbose) file");
			}
		    }
		}
	      

	      /////////////////////////////////////////
	      // Rescale haplotype phase probabilities

	      for (int z=0; z<hap1[i].size(); z++)
		{
		  if ( !(*im)->skip[z] )
		    {
		      pp[i][z] /= s;
		      
		      if ( par::haplo_plem_verbose )
			{
			  if ( (!realnum(pp[i][z])) ||  pp[i][z] < 0 || pp[i][z] > 1 )
			    haplo->VPHASE << "\n*** WARNING *** PROBLEM PP FOR INDIVIDUAL "
				   << P->sample[i]->fid << " " 
				   << P->sample[i]->iid << "\n";
			}
		      
		    }
		}
	    }
	  im++;
	}


      /////////////////////////////////////
      // M-step for pre-counted haplotypes

      // unambiguous counts

      for (int h=0; h<nh; h++)
	f[h] = uc[h];

      ////////////////////////////////////
      // M step for ambiguous genoGroups

      im = genotypes.begin();
      while (im != genotypes.end() )
	{
	  int i = (*im)->reference;
	  if (ambig[i])
	    {
	      if (haplo->haploid || (haplo->X && P->sample[i]->sex))
		{
		  for (int z=0; z<hap1[i].size(); z++)
		    {
		      f[hap1[i][z]] += pp[i][z] * (*im)->count;
		    }
		}
	      else
		{
		  for (int z=0; z<hap1[i].size(); z++)
		    {

		      if ((*im)->skip[z])
			{
			  continue;
			}

// 		      haplo->VPHASE << "considering " << haplotypeName( hap1[i][z] ) 
// 				    << " and " << haplotypeName( hap2[i][z] ) << " for " 
// 				    << P->sample[i]->fid << " " << P->sample[i]->iid << "\t"
// 				    << " times " << (*im)->count << " " << pp[i][z]  
// 				    << " and hap codes are " << hap1[i][z] << " " << hap2[i][z] << "\n" ;
		      
		      f[hap1[i][z]] += pp[i][z] * (*im)->count;
		      f[hap2[i][z]] += pp[i][z] * (*im)->count;
		      
		      
		      
		    }
		}
	    }
	  ++im;
	}


      // validN is the total number of *chromosomes*
      for (int h=0; h<nh; h++)
	f[h] /= (double)haplo->validN;
      

      //////////////////////////////////////////
      // Update likelihood (not every iteration)

      if ( j == par::haplo_plem_iter - 1 ||  
	   j % par::haplo_plem_likelihood_iter == 0)
	{
	  
	  
	  // Zero out unlikely haplotypes?

	  for (int h=0; h<nh; h++)
	    if ( !zero[h])
	      if (f[h] <= par::haplo_plem_zero_threshold )
		zero[h] = true;

	  if ( par::haplo_plem_nonzero_threshold )
	    {
	      double psum = 0;
	      for (int h=0; h<nh; h++)
		if ( ! zero[h] )
		  psum += f[h];
	      for (int h=0; h<nh; h++)
		if ( ! zero[h] )
		  f[h] /= psum;
	    }
			 
	  double lnl = 0;

	  // genoGroups

	  im = genotypes.begin();
	  while (im != genotypes.end() )
	    {
	      int i = (*im)->reference;

	      double lk = 0;

	      if (haplo->haploid || (haplo->X && P->sample[i]->sex))
		{
		  for (int z=0; z<hap1[i].size(); z++)
		    lk += f[hap1[i][z]];
		}
	      else
		for (int z=0; z<hap1[i].size(); z++)
		  {
		    if ((*im)->skip[z] || zero[hap1[i][z]]
			|| zero[hap2[i][z]])
		      continue;

		    lk += f[hap1[i][z]] * f[hap2[i][z]];
		    if (hap1[i][z] != hap2[i][z])
		      lk += f[hap1[i][z]] * f[hap2[i][z]];
		  }

	      if (lk > 0)
		lnl -= log(lk) * (*im)->count;

	      ++im;
	    }


	  if ( par::haplo_plem_verbose )
	    {
	      haplo->VPHASE << "INNER_LNL " << lnl << "\n";
	    }


	  if (j > 0 && sampleLogLikelihood - lnl < par::haplo_plem_window_tol )
	    {

	      if ( par::haplo_plem_verbose )
		haplo->VPHASE << "INNER_CONVERGED AT " << j << " ITERATIONS\n";

	      iter = j;
	      converged = true;
	      break;
	    }

	  sampleLogLikelihood = lnl;


	} // End of likelihood calculation

    } // Next EM iteration

  if ( par::haplo_plem_verbose )
    haplo->VPHASE << "INNER_EM HAS FINISHED/CONVERGED\n\n";

  if ( par::haplo_plem_verbose )
    {
      haplo->VPHASE << "INNER_FREQS ";
      for (int h=0; h<nh; h++)
	if ( f[h] > 0.001 ) haplo->VPHASE << h << " " << haplotypeName(h) << "\t" << f[h] << "\n";
      haplo->VPHASE << "\n--------------------\n";
    }
  

  // EM has converged/finished

}

void HaploWindow::reportPhase()
{
  string fn = par::output_file_name+".phase";
  ofstream PHASE(fn.c_str(), ios::out);

  //P.printLOG("Writing phased haplotypes for " + hname + " to [ "+ fn + " ]\n");

  cout << setw(12) << "FID"<< " "<< setw(12) << "IID"<< " "<< setw(4)<< "PH"
       << " "<< setw(10) << "HAP1"<< " "<< setw(10) << "HAP2"<< " "
       << setw(12) << "POSTPROB"<< " "<< setw(6) << "BEST"<< " "<< "\n";

  PHASE.precision(4);

  for (int i = 0; i < haplo->P.n; i++)
    {

      if (haplo->include[i])
	{

	  for (int z = 0; z < hap1[i].size(); z++)
	    {

	      cout << setw(12) << haplo->P.sample[i]->fid<< " "<< setw(12)
		   << haplo->P.sample[i]->iid<< " "<< setw(4) << z << " "
		   << setw(10) << haplotypeName(hap1[i][z]) << " ";

	      if (haplo->haploid || (haplo->X && haplo->P.sample[i]->sex))
		cout << setw(10) << haplotypeName( -1) << " ";
	      else
		cout << setw(10) << haplotypeName(hap2[i][z]) << " ";

	      if (ambig[i])
		{
		  cout << setw(12) << pp[i][z]<< " ";

		  int max_z = 0;
		  for (int z2=0; z2<hap1[i].size(); z2++)
		    max_z = pp[i][z2] > pp[i][max_z] ? z2 : max_z ;

		  // 		  int fac=1;
		  // 		  if ( hap1[i][max_z] != hap2[i][max_z] ) fac=2;
		  // 		  double w = ( pp[i][max_z] - fac*f[hap1[i][max_z]]*f[hap1[i][max_z]] ) 
		  // 		    / ( 1 - fac*f[hap1[i][max_z]]*f[hap1[i][max_z]] );

		  // Do not output weight for now
		  // if (max_z==z) PHASE << setw(12) << w << " " << setw(6) << 1 << " " << "\n";
		  // else PHASE << setw(12) << "."  << " " << setw(6) << 0  << " " << "\n";

		  if (max_z == z)
		    cout << setw(6) << 1<< " "<< "  ";
		  else
		    cout << setw(6) << 0<< " "<< "  ";
		}
	      else
		cout << setw(12) << 1<< " "<< setw(6) << 1<< " "<< "  ";

	      // Genotypes
	      //for (int s=0; s<ns; s++)
	      //	PHASE << genotype(P, i, S[s]) << " ";
	      cout << "\n";

	    }
	}

      // Report also on excluded individuals
      // (Should be 0-size phase-set)
      else
	{

	  cout << setw(12) << haplo->P.sample[i]->fid<< " "<< setw(12)
	       << haplo->P.sample[i]->iid<< " "<< setw(4) << "NA"<< " "
	       << setw(10) << "NA"<< " "<< setw(10) << "NA"<< " "
	       << setw(12) << "NA"<< " "<< setw(6) << "NA"<< "   ";

	  // genotypes
	  //for (int s=0; s<ns; s++)
	  //	PHASE << genotype(P, i, S[s]) << " ";
	  cout << "\n";

	}
    }

  PHASE.close();

}
