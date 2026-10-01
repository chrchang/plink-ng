


//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
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

#include "plink.h"
#include "helper.h"
#include "options.h"


using namespace std;

void Plink::findIBSRuns(Individual * person1,
			Individual * person2,
			ofstream & IBS)
{
  int l=0; 
  int lastibs=0;
  int lastchr=-1;
  int last;
  int nmiss = 0;
  int nibs0 = 0;
  bool run = false;
  bool justfinished = false;
  int start = 0;
  int end = 0;
  int hetcnt = 0;
  
  vector<bool>::iterator a1 = person1->one.begin();
  vector<bool>::iterator a2 = person1->two.begin();
  vector<bool>::iterator b1 = person2->one.begin();
  vector<bool>::iterator b2 = person2->two.begin();

  while ( a1 != person1->one.end() )
    {
      
      // Skip haploid chromosomes, for now
      if ( ( par::chr_sex[locus[l]->chr] && ( person1->sex || person2->sex ) ) || 
	   par::chr_haploid[locus[l]->chr] ) 
	{
	  a1++;
	  a2++;
	  b1++;
	  b2++;
	  l++;
	  continue; 
	}
      
      bool miss = false;
      bool ibs0 = false;
      
      if ( *a1 == *a2 && 
	   *b1 == *b2 && 
	   *a1 != *b1  ) ibs0 = true;
      else if ( *a1 && !(*a2) ) miss = true;
      else if ( *b1 && !(*b2) ) miss = true;
      else if ( par::ibs_2only )
	{
	  // If we are only looking for IBS2, then make IBS1->0
	  if ( *a1 != *b1 ||
	       *a2 != *b2 )
	    ibs0 = true;
	}

      // Outside of a run?
      if (!run)
	{
	  // A new IBS run?
	  if ( ( ! miss) && ( ! ibs0 ) )
	    {
	      start = lastibs = l;
	      nmiss = nibs0 = 0;
	      run=true;
	    }
	}
      else // if already in a run, either end or increase length?
	{
	  
	  if ( ibs0 ) // ...found IBS0/error?
	    {
	      if (nibs0 == par::ibs_run_0)
		{
		  end = lastibs;
		  run = false;
		  justfinished = true;
		}
	      else
		{
		  nibs0++;
		}
	    }
	  else if ( miss ) // ...missing genotypes?
	    {
	      if (nmiss == par::ibs_run_missing)
		{
		  end = lastibs;
		  run = false;
		  justfinished = true;
		}
	      else
		{
		  nmiss++;
		}
	    }
	  
	  
	  if ( locus[l]->chr != locus[start]->chr ) // different chromosome?
	    {
	      end = l-1;
	      run = false;
	      justfinished = true;
	    }
	  else if ( l == (nl_all -1 ) ) // or end of all SNPs?
	    {
	      if ( ! ( miss || ibs0 ) ) 
		end = l;
	      else
		end = l-1;

	      run = false;
	      justfinished = true;
	    }
	  else if ( l>0 && 
		    locus[l]->chr == locus[l-1]->chr &&
		    ( locus[l]->bp - locus[l-1]->bp ) 
		    > par::ibs_inter_snp_distance )
	    // or too great a gap between SNPs?
	    {
	      end = lastibs;
	      run = false;
	      justfinished = true;
	    }
	  else // ...continue run, recording whether either locus is a het
	    {
	      if ( ! ( miss || ibs0 ) ) 
		{
		  lastibs=l;
		  if ( *a1 != *a2 ||
		       *b1 != *b2 ) 
		    hetcnt++;
		}
	    }
	}
      
      

      
      // Check run length?
      if (justfinished)
	{
	  if ( locus[end]->bp - locus[start]->bp >= 
	       par::ibs_run_length_kb * 1000  && 
	       end - start + 1 >= par::ibs_run_length_snps )
	    {

	      // Record
	      Segment s;
	      s.p1 = person1;
	      s.p2 = person2;
	      s.start = start;
	      s.finish = end;
	      segment.push_back(s);

	      // Display
	      IBS << setw(par::pp_maxfid) << person1->fid << " "
		  << setw(par::pp_maxiid) << person1->iid << " "
		  << setw(par::pp_maxfid) << person2->fid << " "
		  << setw(par::pp_maxiid) << person2->iid << " "
		  << setw(4) << locus[start]->chr << " "
		  << setw(par::pp_maxsnp) << locus[start]->name << " "
		  << setw(par::pp_maxsnp) << locus[end]->name << " "
		  << setw(12) << locus[start]->bp << " "
		  << setw(12) << locus[end]->bp << " "
		  << setw(10) 
		  << (double)(locus[end]->bp - locus[start]->bp) /(double)1000 
		  << " "
		  << setw(5) 
		  << (double)(locus[end]->pos - locus[start]->pos) 
		  << " "
		  << setw(5) << end - start +1 << " "
		  << setw(6) << nibs0 << " "
		  << setw(6) << (double)hetcnt/(double)(end-start+1) << " "
		  << setw(6) << nmiss << " ";
	      
	      if ( lastchr == locus[start]->chr ) 
		{
		  IBS << setw(6) 
		      << start - last - 1 << " ";
		  
		  IBS << setw(6) 
		      << (double)( locus[start]->bp - locus[last]->bp ) / 1000.0
		      << "\n";		  
		}
	      else
		{
		  IBS << setw(6) << "NA" << " ";
		  IBS << setw(6) << "NA" << "\n";
		  lastchr = locus[start]->chr;
		}
	      
	      last = end;

	    }
	
       	  
      
	  //////////////////
	  // Clear counters
	  
	  start = end = nmiss = hetcnt = 0;
	  justfinished=false;

	}
      
      
      ///////////////
	// Next locus
	
	a1++;
	a2++;
	b1++;
	b2++;
	l++;
    }
  
}


void Plink::preCalcPhenotypes()
{
  

  // For binary traits:
  //   SD = X + Y - 2XY
  //   CP = XY - KX - KY + K^2 
  //      = X + Y - (1/K)XY
  //   WT = X + Y - aXY
  //    a = (1-2K + K
  
  
  // Phenotype mean, and number of non-missing phenotypes

  int npheno=0;

  for (int i=0; i<n; i++)
    {
      Individual * person = sample[i];
      
      if (!person->missing)
	{
	  m_phenotype += person->phenotype;
	  npheno++;
	}
    } 
  

  // Test for no non-missing phenotypes: a problem if a subsequent 
  // test has been specified

  if (npheno==0 && 
      (par::assoc_test || 
       par::plink || 
       par::TDT_test || 
       par::ibs_sharing_test ||
       par::epistasis) ) 
    error("No nonmissing phenotypes / available individuals");
  

  // Option to fix the mean/variance or prevalence? Swap in here
  if (par::fix_prev)
    {
      m_phenotype = 1 + par::fixed_prev;
      v_phenotype = npheno * par::fixed_prev * (1-par::fixed_prev);
    }
  
  // Calculate mean of phenotype
  m_phenotype /= npheno;
  
  if (par::qt) 
    {
      stringstream s2;
      s2 << "Phenotype mean = " 
	 << m_phenotype << "\n";
      printLOG(s2.str());
    }
  printLOG("Final analysis contains " + 
	   int2str(npheno) + 
	   " non-missing individuals\n\n");
    
}



short Plink::calcPhenotypes(vector<double> & l, 
			    Individual *p1, 
			    Individual *p2)
{

  /////////////////////////////////////
  // Calculate pairwise phenotype score

  short skip = 0;
  
  if (p1->missing || p2->missing)
    skip = 1;
  else if (par::remove_unaffected_pairs 
	   && p1->phenotype == 1 
	   && p2->phenotype == 1) 
    l.push_back(-999);
  else
    {
      if (par::SD)
	l.push_back( (p1->phenotype - p2->phenotype) 
		     * (p1->phenotype - p2->phenotype ) );
      else if (par::CP) 
	l.push_back(-(p1->phenotype - m_phenotype) 
		    * (p2->phenotype - m_phenotype));
    }

  return skip;
}


void Plink::calcRegression(int chr)
{

  ///////////////////////////////////////////////
  // For a specific chromosome, regress the SD/CP
  // on all pi-hat values (multi- or single-point)

  // Permutes individuals & recalculates pairwise phenotypes 
  // Also, keep track of empirical p-values for the chromosome,
  // for subsequent minP correction
  
  // pihat[pair][position]

  // Number of positions
  int npos = pihat[0].size();

  // Number of replicates
  int R = par::replicates;
 

  vector<double> res;       // Partial correlations, for a chromosome
  vector<double> maxres;    // Largest correlation per replicate

  vector<double> pv(npos,0);  // Empirical p-values
  vector<int> pvalid(npos,0); // Number of valid partial correlations

  
  // Create phenotype list
  phenotype.resize(0);
  for (int i=0; i<pair1.size(); i++)
    calcPhenotypes(phenotype,sample[pair1[i]],sample[pair2[i]]);      
  
  // Precalculate mean and variance of pihat and phenotype
  preCalcRegression_PHENO(phenotype);
  preCalcRegression_PIHAT();  

  if(par::verbose)
    {  
      for (int i=0;i<v_pihat.size();i++)
	cout << i << " V pihat = " << v_pihat[i] << "\n";
      
      // Display entire dataset
      for (int i=0;i<pihat.size();i++)
	{
	  cout << "PIHAT " << i << "\t";
	  for (int j=0;j<pihat[i].size();j++)
	    cout << pihat[i][j] << " ";
	  cout << "\n";
	}
      
      cout << "\n";
      for (int i=0;i<phenotype.size();i++)
	cout << "PHENO " << i << "\t" << phenotype[i] << "\n";
      cout << "\n";


      for (int i=0;i<pihat_G.size();i++)
	cout << "PIHAT_G " << i << "\t" << pihat_G[i] << "\n";
      cout << "\n";

    }
 
 
  // Save original results for this chromosome
  res = doRegression(npos,phenotype);
  
  // Save original results for this chromosome
  original.push_back(res);
  

  // Get lists of permutable individuals
  // i.e. within homogeneous group, as specified
  //      by person->sol

  // Determine number of groups
  int ns=0;
  for (int i=0; i<in_anal.size(); i++)
    if (sample[in_anal[i]]->sol > ns)
      ns=sample[in_anal[i]]->sol;
  ns++;
  
  // Make 's' which records group membership
  // s[group][person]
  vector< vector<int> > s;
  s.resize(ns);
  for (int i=0; i<in_anal.size(); i++)
    s[sample[in_anal[i]]->sol].push_back(i);
  


  //////////////////////
  // Start permutations

  for (int p=1; p<=R; p++)
    {

      cout << "Regression permutation: "
	   << p
	   << " of " 
	   << R 
	   << "             \r";
      
      
      // Vector 'in_anal' contains a list of individuals who 
      // actually feature at least once in the main analysis
      
      // Store remapped IDs
      vector<vector< long int> > indx;
      
      // Permute phenotypes, within cluster
      for (int k=0; k<ns; k++)
        {
          vector<long int> p(s[k].size());
          permute(p);
          indx.push_back(p);
        }
      
      
      // Extract the new permuteds
      vector<int> pin_anal(in_anal.size());

      for (int j=0; j<s.size(); j++)
	for (int k=0; k<s[j].size(); k++)
	  for (int i=0; i<pin_anal.size(); i++)
	    pin_anal[s[j][k]] = in_anal[s[j][(int)indx[j][k]]];
      
      // Make lookup table (that includes missings) for analysis
      vector<int> pall(sample.size(),-1);
      for (int i=0; i<in_anal.size(); i++)
	pall[in_anal[i]] = pin_anal[i];
      
//        for (int i=0; i<pall.size();i++)
//   	cout << "PALL " << i << " " << pall[i] << "\n";
//        cout << "\n";


      // Recreate list of pairs
      
      vector<double> perm;
      perm.resize(0);
      
      for (int i=0; i<pair1.size(); i++)
	{

	  //  	  cout << "about to look up phenos for " 
	  //   	       << sample[pair1[i]]->fid
	  //  	       << "_"
	  //  	       << sample[pair1[i]]->iid
	  //  	       << " - "
	  //   	       << sample[pair2[i]]->fid 
	  //  	       << "_"
	  //  	       << sample[pair2[i]]->iid
	  //   	       << " becomes " 
	  //   	       << sample[pall[pair1[i]]]->fid
	  //  	       << "_"
	  //  	       << sample[pall[pair1[i]]]->iid << " - "
	  //   	       << sample[pall[pair2[i]]]->fid 
	  //  	       << "_"
	  //  	       << sample[pall[pair2[i]]]->iid
	  //   	       << "\n";
	  
	  calcPhenotypes(perm,sample[pall[pair1[i]]],sample[pall[pair2[i]]]);
	}      
      
      // Re-standardise pairs phenotype in the 
      // newly permuted list (as this will potentially 
      // have a different mean and variance compared to the 
      // original list, as each individual may now feature a 
      // different number of times
      
      preCalcRegression_PHENO(perm);


      ///////////////////
      // Perform analyses
      
      // Track the maximum observed test statistic
      
      double mx=0;

      // Get vector of test statistics
      
      vector<double> pres = doRegression(npos,perm);
	
      
      // Compare permuted test statistics against the originals
      double zero=0;

      for (int l=0; l<npos; l++)
	{
	  
	  // A valid partial correlation?
	  if (pres[l] == pres[l] && 
	      pres[l] != 1/zero  && 
	      pres[l] != -1/zero)
	    {
	      // Count this as a valid one
	      pvalid[l]++;
	      
	      // Does this exceed the original? (1-sided)
	      if (pres[l] >= res[l]) pv[l]++;
	      
	      // Is this the maximum?
	      if (pres[l]>mx) mx=pres[l];
	    }
	  
	}
      
     
      // Save maximum statistic
      maxres.push_back(mx);
     
      // Next permutation
    }

  //////////////////////////////////////////////
  // Finished permutations for this chromosome
  
  // Save max-values for genome-wide comparison
  
  for (int p=0; p<R; p++)
    if (maxres[p] >= maxr2[p]) maxr2[p]=maxres[p];
  
  cout << "\n";



  /////////////////
  // Output results
  
  ofstream PLO;
  string f = par::output_file_name + "-" + int2str(chr+1) + ".plink";
  PLO.open(f.c_str(),ios::out);
  printLOG("Writing main PLINK results to [ " + f + " ] \n");

  for (int l=0; l<npos; l++)
    {

      int maxpv=1;
      
      // Empircally-corrected p-value
      for (int p=0; p<R; p++)
	if (maxres[p] >= res[l]) maxpv++;

      double p1, p2;
      string n1, n2;
      double fr = 0, fr2 = 0;

      ////////////////////////////////
      // Single and multipoint output
      
      if (m1[l]==-1) 
	{
	  p1 = locus[par::run_start]->pos - par::fringe;
	  n1 = "fringe";
	}
      else 
	{
	  p1 = locus[m1[l]]->pos;
	  n1 = locus[m1[l]]->name;
	  fr = locus[m1[l]]->freq;
	  if (fr>0.5) fr2=1-fr;
	  else fr2=fr;
	}
      
      if (m2[l]==-1) 
	{
	  p2 = locus[par::run_end]->pos + par::fringe;
	  n2 = "fringe";
	}
      else 
	{
	  p2 = locus[m2[l]]->pos;
	  n2 = locus[m2[l]]->name;
	}
      
      if (m1[l]==-1 && m2[l]==-1)
	{
	  p1 = p2 = 0;
	  n1 = "Genomewide";
	  n2 = "IBD";
	}
      
      double d1 = p1 + pos[l] * (p2-p1);
      
      if (res[l] != res[l]) 
	PLO << "R " 
	    << par::run_chr << " " 
	    << n1 << " " 
	    << n2 << " " 
	    << d1 << "  "
	    << fr << " " << fr2 << " "
	    << "NaN NaN NaN NaN\n";
      else
	PLO << "R " 
	    << par::run_chr << " " 
	    << n1 << " " 
	    << n2 << " " 
	    << d1 << "  "
	    << fr << " " << fr2 << " "
	    << res[l] << " " 
	    << (double)(pv[l]+1)/(double)(pvalid[l]+1) << " "
	    << pvalid[l] << " "
	    << (double)maxpv/(double)(pvalid[l]+1) << "\n";
    }
  
  PLO.close();
}



void Plink::preCalcRegression_PHENO(vector<double> & pheno)
{

  /////////////////////////////////////
  // Precalculate means and variances
  
  // Number of pairs
  int N_pairs = pheno.size();
  
  // Mean
  double m_pair_phenotype=0;
  int ntmp=0;
  for (int i=0; i<N_pairs; i++)
    if (pheno[i]!=-999) // check code for exclusion of 1aff
      {      
	m_pair_phenotype += pheno[i];
	ntmp++;
      }
  m_pair_phenotype /= (double)ntmp;
  
 
  // Variance
  double v_pair_phenotype=0;
  for (int i=0; i<N_pairs; i++)
    if (pheno[i]!=-999)
      v_pair_phenotype += (pheno[i] - m_pair_phenotype) * (pheno[i] - m_pair_phenotype);
  v_pair_phenotype /= (double)ntmp-1;  
  
  // Standardise
  for (int i=0; i<N_pairs; i++)
    if (pheno[i]!=-999)
      pheno[i] = (pheno[i]-m_pair_phenotype)/sqrt(v_pair_phenotype);

//     for (int i=0; i<N_pairs; i++)
//       cout << pheno[i] << " ";
//     cout << "\n";
    
}



void Plink::preCalcRegression_PIHAT()
{

  /////////////////
  // Pi-hat

  // Number of pairs
  int N_pairs = pihat.size();
  
  // Number of positions
  int l = pihat[0].size();

  m_pihat.resize(0);
  v_pihat.resize(0);
  
  // For each locus
  for (int j=0; j<l; j++)
    {
      double m_p=0;
      double v_p=0;
      
      // Pihat mean
      for (int i=0; i<N_pairs; i++) 
	m_p += pihat[i][j]; 
      m_p /= (double)N_pairs;
      m_pihat.push_back(m_p); 
      
      // Pihat variance
      for (int i=0; i<N_pairs; i++)
	v_p += (pihat[i][j] - m_p) * ( pihat[i][j] - m_p);
      v_p /= (double)N_pairs-1;
      v_pihat.push_back(v_p);

      // And standardise
      for (int i=0; i<N_pairs; i++)
      	pihat[i][j] = (pihat[i][j]-m_p)/sqrt(v_p);

    }  
  
  ///////////////////////////
  // Global pi-hat
  double m_p=0;
  double v_p=0;

  // Pihat mean
  for (int i=0; i<N_pairs; i++) 
    m_p += pihat_G[i]; 
  m_p /= (double)N_pairs;
  
  // Pihat variance
  for (int i=0; i<N_pairs; i++)
    v_p += (pihat_G[i] - m_p) * ( pihat_G[i] - m_p);
  v_p /= (double)N_pairs-1;
 
 
  // And standardise
  for (int i=0; i<N_pairs; i++)
    pihat_G[i] = (pihat_G[i]-m_p)/sqrt(v_p);
  
 
}

vector<double> Plink::doRegression(int npos, vector<double> & ph)
{

  vector<double> rv;

  // Dependent variable = squared phenotype difference
  //                    = crossproduct

  // pihat[pair][position]

  // m_phenotype
  // v_phenotype

  // m_pihat[]
  // v_pihat[]

  // Partial corretion ( Phenotype ~ IBD local | IBD global )
  //  r_{PL|G} = ( r_{PL} - r_{PG}r_{LG} )  / sqrt( (1-r_PG^2)(1-r_LG^2)  )

  // r_PL : calculate as before
  // r_PG : calculate only once: make this the first item always
  // r_LG : need to also calculate this each time... 
  
  // Number of pairs
  int N_pairs = pihat.size();
  
  for (int l=0; l<npos; l++)
    {

      // Now we have mean-centered and standardized phenotype and pihat
      // Just use covariance as the statistic

      double r_XY = 0; // pheno ~ PI
      double r_XZ = 0; // PI ~ global
      double r_YZ = 0; // pheno ~ global

      int N_actual = 0;
      
      for (int i=0; i<N_pairs; i++)
	{

	  if (ph[i]>-998){
	    r_XY  += pihat[i][l] * ph[i];
	    r_XZ  += pihat[i][l] * pihat_G[i]; 
	    r_YZ  += ph[i]       * pihat_G[i];
	    N_actual++;
	  }
	}
      
      r_XY /= N_actual-1;
      r_XZ /= N_actual-1;
      r_YZ /= N_actual-1;
      
      if (par::FIXED) r_XZ = r_YZ = 0 ; 
      
      double partial = ( r_XY - r_XZ*r_YZ ) 
	/ sqrt( (1 - r_XZ*r_XZ ) * ( 1 - r_YZ*r_YZ ) ) ;

      
      //      cerr << r_XZ << "\n";
      //       cerr << "r_XY... " 
      // 	   << l << "\t"
      //       	   << r_XY << "\t" 
      //       	   << r_XZ << "\t"
      //       	   << r_YZ << "\t"
      //       	   << partial << "\n";


      // Check that we actually had locus-specific information? 
      // i.e. if the correlation between global and locus-specific
      //      IBD is too high, then we don't want to look at this 
      //      position (set to inf)
      double zero=0;
      if ( r_XZ > par::MAX_CORR_PIHAT_PIHAT_G ) 
	partial = 1/zero;
      else 
	{
	  
	  // Check for +/- inf or nan status -- for now, return 0 if so
	  if (partial != partial || partial == 1/zero || partial == -1/zero) 
	    {
	      
	      cerr << "WARNING: " 
		   << l << "\t" 
		   << partial << "\t"
		   << r_XY << "\t" 
		   << r_XZ << "\t" 
		   << r_YZ << "\n";
	      
	      partial = 1/zero;
	    }
	}

      //save (minus, as default is SD)
      // reverse sign at permutation counting stage
      // if a different test
      rv.push_back(-partial);

    }
  
  //  cout << "END";
  return rv;

  
}




void Plink::calcAssociationWithBootstrap()
{


  // TODO:
  //      1) how to handle missing phenotype and/or 
  //         genotype data in the bootstrap: i.e. ignore?
  //         adjust by factor of n/n* ?  impute? 


  // count number of non-missing individuals in sample
//   int nonmiss = 0;
//   for (int i=0; i<n; i++)
//     {
//       if (!sample[i]->missing) nonmiss++;
//       else if (sample[i]->phenotype!=1 || 
// 	       sample[i]->phenotype!=2 ) 
// 	error("Must be 1/2 coding for association test");
//     }  

//   vector<double> pv(nl_all);
//   vector<double> maxpv(nl_all);
//   vector<double> original;
  
//   int aff;
//   int unf;
//   vector<int> a1(nl_all);
//   vector<int> a2(nl_all);
//   vector<int> a0(nl_all);
  
//   vector<double> odds(nl_all);
 
//   vector<double> exp_afffreq1(nl_all);
//   vector<double> exp_afffreq2(nl_all);
//   vector<double> exp_unffreq1(nl_all);
//   vector<double> exp_unffreq2(nl_all);
  
//   // Original association results
  
//   vector<double> chisq = testAssoc(aff,unf,
// 				   a1,a2,a0,
// 				   odds,
// 				   exp_afffreq1, exp_afffreq2,
// 				   exp_unffreq1, exp_unffreq2,
// 				   perm);
  
  
//   for (int l=0; l< chisq.size(); l++)
//     {
      
//       // Z-score
//       double z = sqrt(chisq[l]);
      
//       // See direction of allele 1 as +ve
//       if ( odds[l] < 1)   // CHECK THIS...
// 	z *= -1;
      
      

// // Make "lower" labelled allele base
//       if (locus[l]->allele1 > locus[l]->allele2)
// 	z *= -1;
      
//       cout << z << " ";
      
//     }
//   cout << "\n";

  
//   // Bootstrap samples
//   //  C  00  AA  GC
//   //  U  AA  00  GG

//   // copy original phenotypes (1/2/0) 
//   vector<int> orig_pheno(sample.size());
//   for (int i=0; i<n; i++)
//     orig_pheno[i] = (int)sample[i]->phenotype;

//   // copy genotype data
//   vector< vector<bool> > orig_geno1(sample.size());
//   vector< vector<bool> > orig_geno2(sample.size());
//   for (int i=0; i<n; i++)
//     {
//       orig_geno1[i].resize(nl_all);
//       orig_geno2[i].resize(nl_all);
//       for (int l=0; l<nl_all; l++)
// 	{
// 	  orig_geno1[i][l] = sample[i]->one[l];
// 	  orig_geno2[i][l] = sample[i]->two[l];
// 	}  
//     }
  
 
//   for (int bs=0; bs<par::replicates; bs++)
//     {
 
//       // Create new bootstrap sample
//       // Ignore missingness, etc.
      
//       for (int i=0; i<n; i++)  
// 	{
// 	  int pick = CRandom::rand(n);
	  
// 	  sample[i]->phenotype = orig_pheno[pick];
	  
// 	  for (int l=0; l < nl_all; l++)
// 	    {
// 	      sample[i]->one[l] = orig_geno1[pick][l];
// 	      sample[i]->two[l] = orig_geno2[pick][l];
// 	    }
// 	}
      
      
//       // Perform association tests
      
//       chisq = testAssoc(aff,unf,
// 			a1,a2,a0,
// 			odds,
// 			exp_afffreq1, exp_afffreq2,
// 			exp_unffreq1, exp_unffreq2,
// 			perm);
      
      
//       // Output raw BS statistic (chisq)
      
//       for (int l=0; l< chisq.size(); l++)
// 	{

//       // Z-score
//       double z = sqrt(chisq[l]);
      
//       // See direction of allele 1 as +ve
//       if ( odds[l] < 1)   // CHECK THIS
// 	z *= -1;
      
//       // Make "lower" labelled allele base
//       if (locus[l]->allele1 > locus[l]->allele2)
// 	z *= -1;
      
//       cout << z << " ";
// 	}
//       cout << "\n";
      
//     }






//   //////////////
//   // Runs of IBS
  
//   if (par::ibs_run)
//     {
      
//       if (par::SNP_major) 
// 	P.SNP2Ind();

//       P.segment.resize(0);

//       ofstream IBS;
//       string f = par::output_file_name + ".ibs";
//       if (par::ibs_2only) f += "2";
//       IBS.open(f.c_str(),ios::out);
//       IBS.precision(4);

//       if (par::ibs_2only) 
// 	P.printLOG("Writing IBS(2)-run information to [ "+f+" ] \n");
//       else
// 	P.printLOG("Writing IBS(1+2)-run information to [ "+f+" ] \n");

//       P.printLOG("Run defined as "+int2str(par::ibs_run_length_kb) + " kb ");
//       P.printLOG("and "+int2str(par::ibs_run_length_snps) + " SNPs ");
//       P.printLOG("(maximum gap of "
// 		 +int2str( (int)((double)par::ibs_inter_snp_distance/1000))+" kb)\n");
      
//       if (par::ibs_2only) 
// 	P.printLOG("Allowing "+int2str(par::ibs_run_0)+" IBS 0/1 SNPs per run\n");
//       else 
// 	P.printLOG("Allowing "+int2str(par::ibs_run_0)+" IBS 0 SNPs per run\n");

//       P.printLOG("Allowing "
// 		 +int2str(par::ibs_run_missing)+" missing SNPs per run\n");

    
//       IBS << setw(par::pp_maxfid) << "FID1" << " "
// 	  << setw(par::pp_maxiid) << "IID1" << " "
// 	  << setw(par::pp_maxfid) << "FID2" << " "
// 	  << setw(par::pp_maxiid) << "IID2" << " "
// 	  << setw(4) << "CHR" << " "
// 	  << setw(par::pp_maxsnp) << "SNP1" << " "
// 	  << setw(par::pp_maxsnp) << "SNP2" << " "
// 	  << setw(12) << "POS1" << " "
// 	  << setw(12) << "POS2" << " "
// 	  << setw(10) << "KB" << " "
// 	  << setw(5) << "CM" << " "
// 	  << setw(5) << "NSNP" << " "
// 	  << setw(6) << "NIBS0" << " "
// 	  << setw(6) << "HET" << " "
// 	  << setw(6) << "NMISS" << " "
// 	  << setw(6) << "SNPGAP" << " "
// 	  << setw(6) << "KBGAP" << "\n";
            
//       int c=0;
//       for (int i1=0; i1<P.n-1; i1++)
// 	for (int i2=i1+1; i2<P.n; i2++)
// 	  {
// 	    if (!par::silent) 
// 	      cout << ++c << " of " << P.np << " pairs         \r";
// 	    P.findIBSRuns(P.sample[i1], P.sample[i2], IBS);
// 	  }
      
//       if (!par::silent) 
// 	cout << "\n\n";
      
//       IBS.close();

//       // Display summary (i.e. per locus)

//       P.printLOG("Found "+int2str(P.segment.size())+" segments\n");
      
//       f += ".summary";
//       P.printLOG("Writing segment summary to [ "+f+" ]\n\n");
//       P.summaryIBSsegments(perm);
      
//       shutdown();
//   }



//   ///////////////////////////////////
//   // Runs of missingness (deletions)
  
//   if (par::miss_run)
//     {

//       if (par::SNP_major) P.SNP2Ind();

//       ofstream RUN;
//       string f = par::output_file_name + ".rum";
//       RUN.open(f.c_str(),ios::out);
      
//       P.printLOG("Writing run-of-missings information to [ "+f+" ] \n");

//       string msg = "Run defined as " + int2str(par::miss_run_length);
//       if (par::miss_run_length_kb) msg += " kb\n";
//       else msg += " SNPs\n";
//       P.printLOG(msg);

//       stringstream s2;
//       s2 << "With at least " << par::miss_run_level << " missingness\n";
//       P.printLOG(s2.str());

//       for (int i1=0; i1<P.n; i1++)
// 	{
// 	  if (!par::silent)
// 	    cout << i1+1 << " of " << P.n << " individuals      \r";
// 	  P.findMissRuns(P.sample[i1],RUN);
// 	}

//       if (!par::silent)
// 	cout << "\n\n";
      
//       RUN.close();

//       shutdown();
//   }
  

}

vector<double> Plink::calcSinglePoint(vector<Z> & IBD, Z IBDg)
{
  vector<double> pihat;
  
  // Singlepoint individual loci
  for (int l=0; l<nl; l++)
    {	        
      pihat.push_back( (IBD[l].z1*0.5+IBD[l].z2) );      
      
      if (par::multi_output) 
	{
 	  cout << "S " 
	       << pairid << " " 
 	       << locus[m1[l]]->name << "\t"
	       << IBD[l].z0 << " "
	       << IBD[l].z1 << " "
	       << IBD[l].z2 << " => "
	       << IBD[l].z1*0.5+IBD[l].z2 << " "
	       << (IBDg.z1*0.5 + IBDg.z2) << "\n";	  
	}
    }
  
    // Final analysis is genome-wide IBD
  if (!par::done_global_pihat) pihat_G.push_back( (IBDg.z1*0.5 + IBDg.z2) );
  
  return pihat;
}

