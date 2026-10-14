

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
#include <cmath>
#include <algorithm>

#include "perm.h"
#include "helper.h"
#include "stats.h"
#include "sets.h"


Perm::Perm(Plink & pref) : P(pref)
{
  
  if (par::adaptive_perm) 
    {
      adaptive = true;
      replicates = par::adaptive_max;
    }
  else
    {
      adaptive = false;
      replicates = par::replicates;
    }
  
  count = par::perm_count;
  
  min = par::adaptive_min;
  
  interval  = par::adaptive_interval;  
  performed = 0;

  dump_all = par::mperm_save_all;
  dump_best = par::mperm_save_best;

  if (dump_all)
    {
      PDUMP.open((par::output_file_name+".mperm.dump.all").c_str(),ios::out);
    }
  else if (dump_best)
    {
      PDUMP.open((par::output_file_name+".mperm.dump.best").c_str(),ios::out);
    }
}


void Perm::setTests(int x) 
{ 

  performed = 0;

  t = x;

  R.clear();
  N.clear();
  test.clear();
  snp_test.clear();
  maxR.clear();

  R.resize(t,0);
  
  if (adaptive)
    {
      N.resize(t,0);
      test.resize(t,true);
      snp_test.resize(t,true);
      
      if ( par::set_test )
	snp_test.resize(P.nl_all,true);
	
      // Given t tests, set the threshold to be 
      // p +/-  Phi^{-1} (1 - \gamma/2t ) sqrt( p(1-p)/N )
      zt = ltqnorm( 1 - par::adaptive_ci / ( 2 * t ) ) ; 
     
    }
  else
    {
      maxR.resize(t,0);
    }

  // For gene-dropping, set up some family-information
  if (par::perm_genedrop)
    preGeneDrop();

}

// Redundant
void Perm::setAdaptiveSetSNPs(int x)
{
//   snp_test.clear();
//   snp_test.resize(x,true);
}


void Perm::originalOrder()
{
  for (int i=0; i<P.n; i++)
    P.sample[i]->pperson = P.sample[i];
}


bool Perm::finished()
{
  if (performed>=replicates) return true;  
  else return false;
}


void Perm::permuteInCluster()
{

  // Store remapped IDs
  vector<vector< long int> > i(ns);
  
  // Permute phenotypes, within cluster
  for (int k=0; k<ns; k++)
    {
      vector<long int> p(s[k].size());
      permute(p);
      i[k]=p;
    }

  //////////////////////////
  // Post-permutation:
  // Iterate over clusters { s[][] }
  // i[][] holds the permuted codes
  // s[][] points to individuals (non-missing)
  
  // Genotype =           sample[s[j][k]];	    
  // Matching phenotype = sample[s[j][(int)i[j][k]]];	    

  // Create pheno[] with label-swapped codes
  for (int j=0; j<s.size(); j++)
    for (int k=0; k<s[j].size(); k++)
      P.sample[s[j][k]]->pperson = P.sample[s[j][(int)i[j][k]]];
  
}

void Perm::setPermClusters(Plink & P)
{
  
  // Permute within clusters only
  // (stored in sample[i]->sol) 
  
  // Get list of non-missing individuals, and number of solutions
  // These are always numbered 0,1,2,3,..
  
  // -1 indicates do not permute this individual
  // 0..ns indicate cluster numbers
  
  // Count the number of clusters: 'ns'
  
  ns=-1;
  for (int i=0; i<P.n; i++)
    if ((!P.sample[i]->missing) && P.sample[i]->sol > ns)
      ns=P.sample[i]->sol;
  ns++;
  
  // store set membership is 's' 
  
  s.resize(ns);
  for (int i=0; i<P.n; i++)
    if (!P.sample[i]->missing && P.sample[i]->sol>=0)
      s[P.sample[i]->sol].push_back(i);
  
  pheno.resize(P.n);

  if (par::permute && ! par::QTDT_test )
    P.printLOG("Set to permute within "+int2str(ns)+" cluster(s)\n");

}

void Perm::setOriginalRanking(vector_t & original)
{
  vector<Pair2> o;
  for (int i=0; i<original.size(); i++)
    {
      Pair2 p;
      p.p = original[i];
      p.l = i;
      o.push_back(p);
    }

  sort(o.begin(),o.end());

  order.clear();
  reorder.resize(original.size());
  for(int i=0; i<original.size(); i++)
    {
      order.push_back(o[i].l);
      reorder[o[i].l] = i;
    }
}



bool Perm::update(vector<double> & result, vector<double> & original)
{
  
  // Increment number of permutations performed
  performed++;
    
  // Finished all perms?
  bool done = false;
  
  //////////////////////////////
  // Update number of successes
  
  if (!adaptive)
    {
      for (int l=0; l<t; l++)
	if (result[l] >= original[l] || !realnum(original[l]) ) R[l]++;
    }
  else
    {
      for (int l=0; l<t; l++)
	if (test[l]) 
	  {
	    if (result[l] >= original[l]  || !realnum(original[l]) ) R[l]++;      
	    N[l]++;
	  }
    }
	

  // Stopping rules for adaptive permutation?
  int todo = 0;
  if (adaptive && 
      performed > min && 
      performed % interval == 0)
    {
      
      // Update interval
      interval = (int)(par::adaptive_interval + performed * par::adaptive_interval2);

      // Consider each test
      for (int l=0; l<t; l++)
	{
	  if (test[l])
	    {
	      
	      // Check for at least one success
	      if (R[l]>0)
		{
		  double pv = (double)(R[l]+1)/(double)(performed+1);
		  
		  double sd = sqrt( pv * (1-pv) / performed );
		  double lower = pv - zt * sd;
		  double upper = pv + zt * sd;
		  //double cv = sd/(performed*pv);
		  if (lower<0) lower = 0;
		  if (lower>1) upper = 0;

		  // Is lower bound greater than threshold, or 
		  // upper bound smaller than threshold?
		  if (lower > par::adaptive_alpha || upper < par::adaptive_alpha ) 
		    {
		      N[l] = performed;
		      test[l] = false;

		      if ( par::set_test )
			{
			  for (int j=0;j< P.pS->snpset[l].size();j++)
			    snp_test[ P.pS->snpset[l][j] ] = false;
			}
		      else 
			snp_test[l] = false;
		    }
		  
		  else
		    todo++;
		}
	      else todo++;
	    }
	}
      
      if (!par::silent)
	{
	  if ( par::set_test ) 
	  cout << "Adaptive permutation: "
	       << performed
	       << " of (max) " 
	       << replicates 
	       << " : " 
	       << todo 
	       << " sets left"
	       << "          \r";
	    else
	  cout << "Adaptive permutation: "
	       << performed
	       << " of (max) " 
	       << replicates 
	       << " : " 
	       << todo 
	       << " SNPs left"
	       << "          \r";
	}

      if (todo==0) done = true;
    }

  
  ///////////////////////////////////////////////////////////
  // For non-adaptive permutation, keep track of the maximum
  
  if (!adaptive)
    {
      
      if (dump_all)
	{
	  if (performed==1)
	    {
	      PDUMP << 0 << " ";
	      for (int l=0; l<t; l++)
		{
		  if( realnum( original[l] ) )
		    PDUMP << original[l] << " ";
		  else
		    PDUMP << "NA" << " ";
		}
	      PDUMP << "\n";
	    }
	  PDUMP << performed << " ";
	}
      
      
      if (dump_best)
	{
	  if (performed==1)
	    {
	      
	      PDUMP << 0 << " ";
	      
	      double mx=0;
	      for (int l=0; l<t; l++)
		if (original[l]>mx) 
		  if ( realnum(mx) ) 
		    mx=original[l];
	      PDUMP << mx << "\n";
	    }
	  PDUMP << performed << " ";
	}
      
      
      // Find maximum, or sort all results
      
      double mx=0;

      if ( par::mperm_rank )
	{
	  // Ranked permutation
	  // populate mx vector

	  // Set any NA to 0
	  for (int l=0; l<t; l++)
	    if ( ! realnum(result[l]) )
	      result[l] = 0;
	  
	  // Sort (ascending)
	  sort(result.begin(),result.end());

	  // Save max in any case
	  mx = result[result.size()-1];
	}
      else
	{
	  // Standard max(T)
	  // populate single mx variable
	  for (int l=0; l<t; l++)
	    if (result[l]>mx) 
	      if ( realnum(mx) ) 
		mx=result[l];
	}


      if (dump_best)
	PDUMP << mx << "\n";
      else if (dump_all)
	{
	  for (int l=0; l<t; l++)
	    {
	      if ( realnum( result[l] ) )
		PDUMP << result[l] << " ";
	      else
		PDUMP << "NA" << " ";
	    }
	  PDUMP << "\n";
	}
      

      // Max(T) permutation -- compare against best

      if ( ! par::mperm_rank )
	{
	  for (int l=0; l<t; l++)
	    if (mx >= original[l] || !realnum(original[l]) ) maxR[l]++;
	}
	  
      // Rank(T) permutation -- compare against similar rank
      
      else
	{
	  for (int l=0; l<t; l++)
	    {
	      //cout << "P " << order[l] << " " << original[order[l]] << " , pr = " << result[l] << " ";
	
	      if (result[l] >= original[order[l]] || !realnum(original[order[l]]) )
		maxR[order[l]]++;
		
	    }
	}



      if (!par::silent)
	{
	  cout << "maxT permutation: "
	       << performed
	       << " of " 
	       << par::replicates 
	       << "            \r";
	  cout.flush();
	}
    }  
  
  // Have we hit the maximum number of replicates?
  if (performed>=replicates) done = true;

  return done;
}

void Perm::nextSNP()
{
  // Reset Perm class for next SNP when in adaptive, SNP-by-SNP mode
  performed = 0;
  originalOrder();
}

bool Perm::updateSNP(double result, double original, int l)
{
  
  /////////////////////////////////////////
  // Single SNP adaptive permutation update
  // for QFAM -- do not allow set-based tests
  // here
  /////////////////////////////////////////

  // Increment number of permutations performed for this SNP
  performed++;
    
  // Finished all perms for this SNP?
  bool done = false;
  
  //////////////////////////////
  // Update number of successes

  if (test[l]) 
    {
      if (result >= original  || !realnum(original) ) R[l]++;      
      N[l]++;
    }
    

  // Stopping rules for adaptive permutation?

  if (adaptive && 
      performed > min && 
      performed % interval == 0)
    {
      
      // Update interval
      interval = (int)(par::adaptive_interval + performed * par::adaptive_interval2);

      // Consider this specific SNP
      if (test[l])
	{
	  
	  // Check for at least one success
	  if (R[l]>0)
	    {
	   
	      double pv = (double)(R[l]+1)/(double)(performed+1);
	      
	      double sd = sqrt( pv * (1-pv) / performed );
	      double lower = pv - zt * sd;
	      double upper = pv + zt * sd;
	      //double cv = sd/(performed*pv);
	      if (lower<0) lower = 0;
	      if (lower>1) upper = 0;
	      
	      // Is lower bound greater than threshold, or 
	      // upper bound smaller than threshold?
	      if (lower > par::adaptive_alpha || upper < par::adaptive_alpha ) 
		{
		  N[l] = performed;
		  test[l] = false;
		  done = true;
		}
	    }
	}
            
    }
  
  // Have we hit the maximum number of replicates?
  if (performed>=replicates) done = true;
  
  return done;
}


int Perm::rank(int l)
{
  if ( ! par::mperm_rank )
    return 0;
  else
    return t - reorder[l];
  
}

double Perm::pvalue(int l)
{
  if (count)
    return (double)R[l];
  if (adaptive)
    return (double)(R[l]+1) / (double)(N[l]+1);
  else
    return (double)(R[l]+1) / (double)(replicates+1);
}


double Perm::max_pvalue(int l)
{
  if (adaptive)
    return -1;
  else
    return (double)(maxR[l]+1) / (double)(replicates+1);
}

int Perm::reps_done(int l)
{
  if (adaptive)
    return N[l];
  else
    return replicates;
}

