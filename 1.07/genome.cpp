

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
#include <vector>
#include <map>
#include <cmath>
#include <iterator>
#include "plink.h"
#include "options.h"
#include "helper.h"
#include "stats.h"
#include "cfamily.h"
#include "zed.h"

void Plink::calcStratifiedAlleleFreqs()
{
  
  // w/ Modification to output counts by John Novembre

  // Assume SNP-major data

  if (!par::SNP_major) Ind2SNP();
  

  // This is called *after* any filters are applied 
  // (i.e. unlike the original --freq command)
  // This means that various things will have been fixed already, 
  // such as het haploid calls, and which allele is the minor one
  
  if (par::summ_nonfounders)
    printLOG("Writing stratified allele frequencies (all individuals) to [ " +
	     par::output_file_name + ".frq.strat ] \n");
  else
    printLOG("Writing stratified allele frequencies (founders-only) to [ " +
	     par::output_file_name + ".frq.strat ] \n");
  
  ofstream FRQ;
  string f = par::output_file_name + ".frq.strat";
  FRQ.open(f.c_str(), ifstream::out);
  FRQ.precision(4);

  FRQ << setw(4) << "CHR" << " "
      << setw(par::pp_maxsnp) << "SNP" << " "
      << setw(8) << "CLST" << " "
      << setw(4) << "A1" << " "
      << setw(4) << "A2" << " "
      << setw(8) << "MAF" << " "
      << setw(6) << "MAC" << " "
      << setw(8) << "NCHROBS"   // Modified by JN 5/24/07
      << "\n";	  
  


  //////////////////////////////////////////
  // Calculate allele frequencies, and FST
  
  vector<Locus*>::iterator loc = locus.begin();
  vector<CSNP*>::iterator s = SNP.begin();
  
  while ( loc != locus.end() ) 
    {
      
      // Track Fst

      vector_t het(nk);
      double tothet = 0;

      // Consider each cluster

      for (int k=0; k<nk; k++)
	{
	  
	  FRQ << setw(4) << (*loc)->chr << " "
	      << setw(par::pp_maxsnp) << (*loc)->name << " ";
	  
	  (*loc)->freq = 0;
	  
	  // count 1 per allele, for frequency
	  double nmfreq = 0; 
	  double totfreq = 0;
	  double count = 0;

	  // count 1 per genotype, for missingness
	  int geno_nm = 0; 
	  bool X = false;
	  bool haploid = false;
	  
	  // Determine type of SNP
	  if (par::chr_sex[(*loc)->chr]) X=true;
	  else if (par::chr_haploid[(*loc)->chr]) haploid=true;
	  
	  
	  ///////////////////////////////
	  // Iterate over each individual
	  
	  vector<bool>::iterator i1 = (*s)->one.begin();
	  vector<bool>::iterator i2 = (*s)->two.begin();
	  vector<Individual*>::iterator person = sample.begin();
	  
	  while ( person != sample.end() ) 
	    {
	      
	      bool s1 = *i1;
	      bool s2 = *i2;
	      
	      // For allele frequencies
	      // only consider founders?	
	      
	      if ( par::summ_nonfounders || (*person)->founder ) 
		{

		  if ( (*person)->sol == k ) 
		    {
		  
		      if ( haploid || ( X && (*person)->sex ) )
			{
			  
			  //////////////////
			  // Haploid counts
			  
			  // "1" allele count
			  
			  // Possible count of 1 allele
			  totfreq++;
			  
			  if ( (!s1) && (!s2) )   //  FF = hom(11)
			    {
			      (*loc)->freq++;
			      nmfreq++;
			    }	
			  else if ( s1 && s2 )   //  TT = hom(22)
			    {
			      nmfreq++;
			    }
			  
			}
		      else
			{
			  
			  //////////////////
			  // Autosomal count
			  
			  // "1" allele count
			  
			  // Possible count of 2 alleles
			  totfreq+=2;
			  
			  if (!s1)
			    { 
			      if (!s2)  //   00 = hom(11)
				{
				  (*loc)->freq+=2;
				  nmfreq+=2;
				}	
			      else                  //   01 = het(12)
				{
				  (*loc)->freq+=1;
				  nmfreq+=2;
				}
			    }
			  else if ( s2 ) // 11 = hom(22)
			    {
			      nmfreq+=2;
			    }
			}
		      
		    }
		  
		}

	      // Next individual
	      person++;
	      i1++;
	      i2++;
	      
	    }
	  
	  if (nmfreq>0){
	    count=(*loc)->freq;  // Added by JN 5/24/07	    
	    (*loc)->freq /= (double)nmfreq;
	  }

	  string a1 = (*loc)->allele1;
	  if (a1=="") a1="0";
	  FRQ << setw(8) << kname[k] << " "
	      << setw(4)  << a1  << " "
	      << setw(4)  << (*loc)->allele2  << " "
	      << setw(8) << (*loc)->freq << " "
              << setw(6) << count << " "
              << setw(8) << nmfreq << " " // Modified by JN 5/24/07
	      << "\n";
	  
	}
                
      
      // Next SNP
      loc++;
      s++;
    }

  FRQ.close();
  shutdown();
  
}




void Plink::findMissRuns(Individual * person, ofstream & RUN)
{

  int l=0; 
  int nmiss = 0;  // now means 'not missing'
  bool run = false;
  int start = 0;
  int end = 0;

  while ( l < nl_all )
    {
      
      // Outside of a run?
      if (!run)
	{
	  // A new run? (missing)
	  if (person->one[l] && (!person->two[l]))
	    {
	      start = l;
	      nmiss=1;
	      run=true;
	    }
	}
      else // if already in a run, either end or increase length?
	{
	  // found a non-missing? 00, 11, 01

	  if ( person->one[l] == person->two[l] || person->two[l])
	    {
	      nmiss++;

	      // Average non-missing rate now too high?, given we have at least 
	      // a certain number of SNPs in run
	      // (l-start) = number of SNPs in run currently
	      
	      if ((double)nmiss / (double)(l-start) >= par::miss_run_level && (l-start) >= par::miss_run_length)
		{
		  end = l-1;
		  run = false;
		}
	    }
	  else if ( locus[l]->chr != locus[start]->chr ) // different chromosome?
	    {
	      end = l-1;
	      run = false;
	    }
	  else if ( l == (nl_all -1) ) // or end of all SNPs?
	    {
	      end = l;
	      run = false;
	    }
	}


      // Check run length?
      if (!run)
	{
	  if (par::miss_run_length_kb)
	    {
	      if ( locus[end]->bp - locus[start]->bp >= par::miss_run_length * 1000 )
		RUN << setw(par::pp_maxfid) << person->fid << " "
		    << setw(par::pp_maxiid) << person->iid << " "
		    << setw(8) << person->phenotype << " "
		    << setw(4) << locus[start]->chr << " "
		    << setw(par::pp_maxsnp) << locus[start]->name << " "
		    << setw(par::pp_maxsnp) << locus[end]->name << " "
		    << setw(10) << (double)(locus[end]->bp - locus[start]->bp)/(double)1000 << " "
		    << setw(10) << end - start +1 << " "
		    << setw(10) << (double)nmiss/(double)(end - start + 1) << "\n";
	      
	    }
	  else 
	    {
	      if ( end - start +1 >= par::miss_run_length  )
		RUN << person->fid << "\t"
		    << person->iid << "\t"
		    << locus[start]->chr << "\t"
		    << locus[start]->name << "\t"
		    << locus[end]->name << "\t"
		    << (double)(locus[end]->bp - locus[start]->bp)/(double)1000 << "\t"
		    << end - start + 1 << "\t"
		    << (double)nmiss/(double)(end - start + 1) << "\n";
	    }
	  
	  
	  
	  //////////////////
	  // Clear counters
	  
	  start = end = nmiss = 0;
	  
	}

      ///////////////
      // Next locus
      
      l++;
      
    }

}


void Plink::sexCheck()
{

  // Get range of X chromosome markers

  if (par::SNP_major) 
    SNP2Ind();

  ofstream HET;
  string f = par::output_file_name + ".sexcheck";
  HET.open(f.c_str(),ios::out);
  HET.precision(4);
        
  printLOG("Writing X-chromosome sex check results to [ "+f+" ] \n");
      
  HET << setw(par::pp_maxfid) << "FID" << " "
      << setw(par::pp_maxiid) << "IID" << " "
      << setw(12) << "PEDSEX" << " " 
      << setw(12) << "SNPSEX" << " " 
      << setw(12) << "STATUS" << " " 
      << setw(12) << "F" << "\n";
            
  for (int i1=0; i1<n; i1++)
    calcInbreeding(sample[i1],0,nl_all-1,HET);
  
  HET.close();
  
  return;
}


void Plink::calcFst()
{

  ofstream FST;
  string f = par::output_file_name + ".fst";
  FST.open(f.c_str(),ios::out);
  FST.precision(4);
  
  printLOG("Writing Fst measures to [ "+f+" ] \n");
  
  FST << setw(4) << "CHR" << " "
      << setw(par::pp_maxsnp) << "SNP" << " "
      << setw(20) << "CLST" << " " 
      << setw(8) << "N" << " " 
      << setw(8) << "HET" << "\n"; 
            
  for (int l = 0 ; l < nl_all; l++)
    {

      bool skip = false;

      if ( par::chr_sex[locus[l]->chr] || 
	   par::chr_haploid[locus[l]->chr] ) skip = true;
      if ( locus[l]->nm <= 1 || locus[l]->freq < 1e-8 )
	skip = true;

      if ( skip ) 
	continue;

      
      double f = 0;
      
      vector_t hs(nk);
      double ht = 2 * locus[l]->freq * ( 1 - locus[l]->freq );

      // Calculate Fst here
      
//       double avg_hs = 0;
//       for (int k = 0; k < nk ; k++ )
// 	avg_hs += het[k];
//       avg_hs /= (double)nk;

//       double f = ( hs - ht ) / ht;


      for ( int k = 0 ; k < nk ; k++ )
	{
// 	  FST << setw(4) << locus[l]->chr  << " "
// 	      << setw(par::pp_maxsnp) << locus[l]->name << " "
// 	      << setw(20) << kname[k] << " " 
// 	      << setw(8) << kind[k] << " " 
// 	      << setw(12) << het[k] << "\n";
	}

      
      FST << setw(4) << locus[l]->chr << " "
	  << setw(par::pp_maxsnp) << locus[l]->name << " "
	  << setw(20) << "_FST_" << " " 
	  << setw(8) << locus[l]->nm << " " 
	  << setw(12) << f << "\n";
      
    }
  
  FST.close();
  
  return;
}




double Plink::calcInbreeding(Individual * p1, int m1, int m2, ofstream & HET)
{

  // P(Homo) = F + (1-F)P(Homo by chance)

  // P(Homo by chance) = p^2+q^2 for a biallelic locus.
  // For an individual with N genotyped loci, we
  //   1. count the total observed number of loci which are homozygous (O),
  //   2. calculate the total expected number of loci homozygous by chance (E)
  // Then, using the method of moments, we have
  //    O = NF + (1-F)E
  // Which rearranges to give
  //    F = (O-E)/(N-E)

  
  // Count of nonmissing loci
  double N=0;
  double O=0;
  double E=0;

  // Consider all loci
  for (int l=m1; l<=m2;l++)
    {

      ////////////////////////////////////////////////
      // Skip X and haploid chromosome markers, or not
      
      if ( par::check_sex ) // For sex-checks
	{
	  // only consider the X chromosome
	  if ( ! par::chr_sex[locus[l]->chr] ) continue; 
	}
      else // normal heterozygosity calculation
	{
	  // so skip haploid markers
	  if ( par::chr_sex[locus[l]->chr] || 
	       par::chr_haploid[locus[l]->chr] ) continue; 
	}
      
      // Skip monomorphic markers, uninformative markers
      if ( locus[l]->nm <= 1 || locus[l]->freq < 1e-8 )
	continue;

      //////////////////
      // Observed data
      
      // check not missing:
      if (!(p1->one[l] && (!p1->two[l])))
	{
	  // homozygous non-missing loci
	  if (p1->one[l] == p1->two[l]) O++;
	  // non-missing loci
	  N++;
	  

	  /////////////////////////
	  // Expected homozygousity
	  
	  // E = 2pq . 2N/(2N-1)
	  // (Using Nei's unbiased estimator)
	  
 	  E += 1
	    - ( 2 * locus[l]->freq * ( 1 - locus[l]->freq ) 
		* ( locus[l]->nm / ( locus[l]->nm - 1 ) ) );
	  

	}
           
    }

  double F = (O-E)/(N-E);

  if ( par::check_sex)
    {
      HET << setw(par::pp_maxfid) << p1->fid << " " 
	  << setw(par::pp_maxiid) << p1->iid << " "      
	  << setw(12) << p1->sexcode << " ";

      if ( F > par::sex_threshold_male ) 
	{
	  HET << setw(12) << 1 << " ";
	  if (p1->sexcode == "1")
	    HET << setw(12) << "OK" << " ";
	  else
	    {
	      HET << setw(12) << "PROBLEM" << " ";
	      if (par::impute_sex)
		{
		  p1->sexcode = "1";
		  p1->sex = true;
		}
	    }
	}
      else if ( F < par::sex_threshold_female )
	{
	  HET << setw(12) << 2 << " ";
	  if (p1->sexcode == "2")
	    HET << setw(12) << "OK" << " ";
	  else
	    {
	      HET << setw(12) << "PROBLEM" << " ";		  
	      if (par::impute_sex)
		{
		  p1->sexcode = "2";
		  p1->sex = false;
		}
	    }
	}
      else 
	{
	  HET << setw(12) << 0 << " "
	      << setw(12) << "PROBLEM" << " ";		  
	  
	  if (par::impute_sex)
	    {
	      p1->sexcode = "0";
	      p1->sex = false;
	      if (!par::ignore_missing_sex)
		p1->missing = true;	      
	    }
	  
	}


      HET << setw(12) << F << "\n";
    }
  else
    HET << setw(par::pp_maxfid) << p1->fid << " " 
	<< setw(par::pp_maxiid) << p1->iid << " "
	<< setw(12) << (int)O << " " 
	<< setw(12) << E << " "
	<< setw(12) << (int)N << " "
	<< setw(12) << F << "\n";
  
  return F;
  
}


Z Plink::calcGenomeIBS(Individual * p1, Individual * p2)
{

  // Vector of average genome-wide IBS
  Z IBSg;
  
  // Count of nonmissing loci
  int cnt=0;

  // Other metrics ( in Plink:: )
  pvIBS0 = 0;
  pvIBS2het = 0;
  int last_chr = -1;
  int last_bp = -1;
  

  // Consider all autosomal loci
  vector<bool>::iterator ia1 = p1->one.begin();
  vector<bool>::iterator ia2 = p1->two.begin();

  vector<bool>::iterator ib1 = p2->one.begin();
  vector<bool>::iterator ib2 = p2->two.begin();
  int l=0;
  
  while ( ia1 != p1->one.end() )
    {
      
      // Skip X and haploid chromosome markers
      if ( par::chr_sex[locus[l]->chr] || 
	   par::chr_haploid[locus[l]->chr] ) 
	{
	  l++;
	  ia1++;
	  ia2++;
	  ib1++;
	  ib2++;
	  continue; 
	}

      // Only count if both genotypes nonmissing
      bool a1 = *ia1;
      bool a2 = *ia2;
      if (a1 && !a2) 
	{
	  l++;
	  ia1++;
	  ia2++;
	  ib1++;
	  ib2++;
	  continue; 
	}

      bool b1 = *ib1;
      bool b2 = *ib2;
      if (b1 && !b2) 
	{
	  l++;
	  ia1++;
	  ia2++;
	  ib1++;
	  ib2++;
	  continue; 
	}

      
      // Calculate IBS from genotypes

      // 10 = missing
      // 00 = 11hom
      // 01 = 12het
      // 11 = 22hom
      
      if ( a1 == b1 && a2 == b2 ) IBSg.z2++;      // IBS 2
      else if ( a1 != b1 && a2 != b2 ) IBSg.z0++; // IBS 0
      else IBSg.z1++;                             // IBS 1
      
      cnt++;
      
      // Also calculate p-value binomial test
      if ( ! ( par::matrix || par::cluster || par::genome_output ) ) 
	{
	  ia1++;
	  ia2++;
	  ib1++;
	  ib2++;
	  l++;	  
	  continue;
	}
      
      if ( a1 != b1 && a2 != b2 ) // IBS 0 hom/hom
	{
	  // Can we count this?
	  if (locus[l]->chr != last_chr ||
	      locus[l]->bp > last_bp + par::ibstest_gap)
	    {
	      pvIBS0++;
	      last_chr = locus[l]->chr;
	      last_bp = locus[l]->bp;
	    }
	}
      else if ( a1 != a2 && b1 != b2 )  // IBS 2 het/het
	{
	  // Can we count this?
	  if (locus[l]->chr != last_chr ||
	      locus[l]->bp > last_bp + par::ibstest_gap)
	    {
	      pvIBS2het++;
	      last_chr = locus[l]->chr;
	      last_bp = locus[l]->bp;
	    }
	}
      
      // Next SNP
      ia1++;
      ia2++;
      ib1++;
      ib2++;
      l++;
      
    }
  
  

  
  if (cnt==0)
    { 
      string msg = "No nonmissing markers for individuals "
	+ p1->fid + " "
	+ p1->iid + " - " 
	+ p2->fid + " "
	+ p2->iid;
      error(msg);
    }


  // Standard genetic distance (proportion IBS 0 to 1) 
  if ( par::cluster_euclidean )
    dst = sqrt((IBSg.z1*0.5 + IBSg.z2*2)/(IBSg.z0+IBSg.z1+IBSg.z2*2));
  else 
    dst = (IBSg.z1*0.5 + IBSg.z2)/(IBSg.z0+IBSg.z1+IBSg.z2);


  // Also calculate p-value binomial test
  if ( par::cluster_missing || ! ( par::matrix || par::cluster || par::genome_output ) ) 
    return IBSg;

    
  // Calculate p-value for IBS test
  // IBS0 : IBS2(het) in 1:2 ratio 
  // n.b.  0.2222222 = 0.666*(1-0.666)
  double z = (pvIBS2het/(pvIBS0+pvIBS2het)-0.666666) 
    / (sqrt(0.2222222/(pvIBS0+pvIBS2het)));
  
  // Store p-value in Plink::pv
  pv = normdist(z);
  
  // Return counts
  return IBSg;
}


void Plink::preCalcGenomeIBD()
{
  
  // Take pairwise IBS and information on allele frequencies
  // to generate genome-wide IBD estimates
  
  // Expected IBS given IBD

  // Bias-corrected versison of estimator for expected proportion of
  // IBS SNP pairs given IBD status (cf. Nei bias-corrected
  // heterozygosity estimator).


  // All possible permutations of taking 4 alleles from 2N = 2N(2N-1)(2N-2)(2N-3)

  // x = count of allele 1
  // y = count of allele 2
  // p = x/2N, q=y/2N

  // Of these, 
  //            x(x-1)y(y-1) will be of order 11 22, 
  //            y(y-1)x(x-1) will be of order 22 11

  // Therefore the probability of IBS 0 given IBS 0 is 
  // E00 = 2 p*p*q*q * ( (x-1)/x * (y-1)/y * (2N/(2N-1)) * (2N/(2N-2)) * (2N/(2N-3)) )
  // and so on for E01, etc
  
  // E(IBD)(IBS)
  // IBS >= IBD
  
  E00=E10=E20=E01=E11=E21=E02=E12=E22=0;
  int cnt = 0;
  for (int l=0; l<nl_all; l++)
    {
      
      if ( par::chr_sex[locus[l]->chr] || 
	   par::chr_haploid[locus[l]->chr] ) continue;
      
      double p = locus[l]->freq;
      double q = 1 - p;

      double Na = locus[l]->nm; // = # alleles = 2N where N is number of individuals
      double x = p * Na;
      double y = q * Na;
      
      // Original, non bias-corrected versions
      //       E00 += 2*p*p*q*q; 
      //       E01 += 4*p*p*p*q+4*p*q*q*q;
      //       E02 += q*q*q*q + p*p*p*p + 4*p*p*q*q;
      //       E11 += 2*p*p*q + 2*p*q*q; 
      //       E12 += p*p*p + q*q*q + p*p*q + p*q*q;	
      

      double a00 = 2*p*p*q*q * ( (x-1)/x * (y-1)/y * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
      
      double a01 = 4*p*p*p*q * ( (x-1)/x * (x-2)/x *  (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) )
	+ 4*p*q*q*q * ( (y-1)/y * (y-2)/y *  (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );
      
      double a02 = q*q*q*q * ( (y-1)/y * (y-2)/y * (y-3)/y * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) )
	+ p*p*p*p  * ( (x-1)/x * (x-2)/x * (x-3)/x * (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) )
	+ 4*p*p*q*q * ( (x-1)/x * (y-1)/y *           (Na/(Na-1)) * (Na/(Na-2)) * (Na/(Na-3)) );

      double a11 = 2*p*p*q * ( (x-1)/x *  Na/(Na-1) * Na/(Na-2) )
	+ 2*p*q*q * ( (y-1)/y *  Na/(Na-1) * Na/(Na-2) );
      
      double a12 = p*p*p * ((x-1)/x * (x-2)/x *  Na/(Na-1) * Na/(Na-2)) 
	+ q*q*q * ( (y-1)/y * (y-2)/y *  Na/(Na-1) * Na/(Na-2))
	+ p*p*q * ( (x-1)/x * Na/(Na-1) * Na/(Na-2) )
	+ p*q*q * ((y-1)/y  * Na/(Na-1) * Na/(Na-2));

      if ( realnum(a00) &&
           realnum(a01) &&
           realnum(a02) &&
           realnum(a11) &&
           realnum(a12) )
        {
          E00 += a00;
          E01 += a01;
          E02 += a02;
          E11 += a11;
          E12 += a12;
          cnt++;
        }

      
    }
  
  E00 /= cnt; E10  = 0;   E20 = 0;
  E01 /= cnt; E11 /= cnt; E21 = 0;
  E02 /= cnt; E12 /= cnt; E22 = 1;
  
  if (par::verbose) { 
    cout << "P(IBS|IBD)  -- IBS row; IBD col\n";
    cout << E00 << "\t" << E10 << "\t" << E20 << "\n";
    cout << E01 << "\t" << E11 << "\t" << E21 << "\n";
    cout << E02 << "\t" << E12 << "\t" << E22 << "\n";
    cout << "\n";
  }
  
}



Z Plink::calcGenomeIBD(Individual * p1, Individual * p2, Z IBSg)
{
  
  Z z;
  
  double S = IBSg.z0 + IBSg.z1 + IBSg.z2;

  // E_IBS[row=IBS][col=IBD]
  // E(IBD)(IBS)
  
  double e00 = E00*S; 
  double e10 = E10*S; 
  double e20 = E20*S;
  
  double e01 = E01*S; 
  double e11 = E11*S; 
  double e21 = E21*S;
  
  double e02 = E02*S; 
  double e12 = E12*S; 
  double e22 = E22*S;
  
  z.z0 =  IBSg.z0 / e00;
  z.z1 = (IBSg.z1 - z.z0*e01) / e11;
  z.z2 = (IBSg.z2 - z.z0*e02 - z.z1*e12) / e22; 

  if (par::debug)
    cout << "DEBUG\t"
	 << z.z0 << " "
	 << z.z1 << " " 
	 << z.z2 << "\n";


  // Bound IBD estimates to sum to 1 
  // and fall within 0-1 range
  if (par::bound)
  {  
  if (z.z0>1) { z.z0=1; z.z1=z.z2=0; }
  if (z.z1>1) { z.z1=1; z.z0=z.z2=0; }
  if (z.z2>1) { z.z2=1; z.z0=z.z1=0; }
  
  if (z.z0<0) { double S=z.z1+z.z2; z.z1/=S; z.z2/=S; z.z0=0; }
  if (z.z1<0) { double S=z.z0+z.z2; z.z0/=S; z.z2/=S; z.z1=0; }
  if (z.z2<0) { double S=z.z0+z.z1; z.z0/=S; z.z1/=S; z.z2=0; }
  }
  
  // Possibly constrain IBD estimates to within possible triangle
  // i.e. 0.5 0.0 0.5 is invalid
  //
  // For purposes of sample checks, etc, we do not automatically do this (--genome)
  // For PLINK analysis we do (--plink)
  //
  // Constraint : z1^2 - 4 z0 z2 >= 0
  //            : x^2 - 2 pi x + z2  = 0
  //
  //              where pi = (z1 + 2 z2) / 2
  //
  // So the constaint can also be written as
  //
  //              pi^2 >=  z2

  
  double pihat =   z.z1/2 + z.z2 ;
  double impossible = false;
  
  if ( ( pihat * pihat ) < z.z2 ) 
    {
      impossible = true;
      
      // find new value for z1 (z1*) which satisfies the equation 
      //
      //     (z1* + 2 pi^2) / 2  = pi
      //
      // this gives
      //
      //      z1* = 2pi(1-pi)
      
      // the transformed IBD probabilities would be
      // 1 - 2pi(1-pi) - pi^2, 2pi(1-pi), pi^2
      
      if (par::nudge)
	{
	  z.z0 = ( 1 - pihat) * ( 1 - pihat);
	  z.z1 = 2 * pihat * (1-pihat);
	  z.z2 = pihat * pihat;
	}
    }


  if ( par::genome_output_minimal )
    {
      ZOUTFILE << dbl2str_fixed(dst,6) << " " 
	       << dbl2str(pv) << " " 
	       << dbl2str(z.z1/2 + z.z2) << "\n";     
    }
  else if (par::genome_output)
    {
      
      if ( (!par::pihat_filter) || 
	   ( pihat >= par::MIN_PIHAT && pihat <= par::MAX_PIHAT) )
	{
	  
	  ZOUTFILE << sw( p1->fid, par::pp_maxfid)
		   << sw( p1->iid, par::pp_maxiid)
		   << sw( p2->fid, par::pp_maxfid)
		   << sw( p2->iid, par::pp_maxiid);

	  string rt = relType(p1,p2);

	  ZOUTFILE << sw(rt,3);
	  
	  if ( rt == "UN" )
	    ZOUTFILE << sw("NA", 6);
	  else
	    ZOUTFILE << sw( genrel(p1,p2) , 6 );
	  

	  if (par::show_impossible_IBD || !impossible)
	    ZOUTFILE << sw(z.z0, 4,8) 
		     << sw(z.z1, 4,8) 
		     << sw(z.z2, 4,8);
	  else
	    ZOUTFILE << sw( -z.z0, 8)
		     << sw( -z.z1, 8)
		     << sw( -z.z2, 8);
	  
	  ZOUTFILE << sw( z.z1/2 + z.z2, 4, 8);
	  
	  if (par::bt)
	    {
	      if ( (!p1->aff) && (!p2->aff) ) 
		ZOUTFILE << sw("-1", 4);
	      else if ( p1->aff && p2->aff )
		ZOUTFILE << sw("1",4);
	      else if ((!p1->aff) && p2->aff)
		ZOUTFILE << sw("0",4);
	      else if (p1->aff && !p2->aff)
		ZOUTFILE << sw("0",4);
	      else
		ZOUTFILE << sw("NA",4);
	    }
	  else
	    ZOUTFILE << sw("NA",4);
	  
	  ZOUTFILE << sw(dst,6,10);
	  ZOUTFILE << sw(pv,4,8);
	  
	  double ov = (double)pvIBS2het / (double)pvIBS0;
	  if ( realnum(ov) )
	    ZOUTFILE << sw(ov,4,8);
	  else
	    ZOUTFILE << sw("NA",8);
	  
	  if ( par::genome_output_full )
	    ZOUTFILE << sw((int)IBSg.z0, 8)
		     << sw((int)IBSg.z1, 8) 
		     << sw((int)IBSg.z2, 8) 
		     << sw(pvIBS0, 4,8) 
		     << sw(pvIBS2het, 4, 8);
	  
	  ZOUTFILE << "\n";
	}
    } 
  return z;
}


void Plink::displayGenomeWideInfo()
{

  ///////////////////////////////////////
  // This is an individual-mode analysis

  if (par::SNP_major) SNP2Ind();

  string f = par::output_file_name + ".genome";
  if ( par::genome_output_minimal ) f += ".min";
  if ( par::compress_genome ) f += ".gz";

  if ( par::genome_output_minimal )
    printLOG("Writing minimal-format IBS information to [ " + f + " ] \n");      
  else
    printLOG("Writing whole genome IBS/IBD information to [ " + f + " ] \n");

  ZOUTFILE.open( f , par::compress_genome );

  stringstream s2;
  s2 << "Filtering output to include pairs with ( " 
     << par::MIN_PIHAT  << " <= PI-HAT <= " 
     << par::MAX_PIHAT << " )\n";
  printLOG(s2.str());
  
  if ( par::genome_output_minimal )
    {
      for (int i=0; i<n; i++)
	ZOUTFILE << sample[i]->fid << " " << sample[i]->iid << "\n";
      ZOUTFILE << "__END __END\n";
    }
  else
   {
     ZOUTFILE << sw("FID1",par::pp_maxfid)
	      << sw("IID1",par::pp_maxiid)
	      << sw("FID2",par::pp_maxfid)
 	      << sw("IID2",par::pp_maxiid)
 	      << sw("RT",3)
 	      << sw("EZ",6)
 	      << sw("Z0",8)
 	      << sw("Z1",8)
 	      << sw("Z2",8)
 	      << sw("PI_HAT",8)
 	      << sw("PHE",4)
 	      << sw("DST",10)
 	      << sw("PPC",8)
 	      << sw("RATIO",8);

     if ( par::genome_output_full )
       ZOUTFILE << sw("IBS0",8)
		<< sw("IBS1",8)
		<< sw("IBS2",8)
		<< sw("HOMHOM",8)
		<< sw("HETHET",8);
     
     ZOUTFILE << "\n"; 
   }

  
  int c=0;
  int c2=0;
  for (int i1=0; i1<n-1; i1++)
    for (int i2=i1+1; i2<n; i2++)
      {
	Individual * p1 = sample[i1];
	Individual * p2 = sample[i2];
	
	// Only update message every 100 iterations
	if ( (!par::silent ) && c==c2 || c==np)
	  {
	    cout << "IBD(g) calculation: " 
		 << c++ << " of " << np 
		 << "                  \r";
	    cout.flush();
	    c2+=100;
	  }
	else
	  ++c;

	
	/////////////////////////
	// Or perhaps we can skip?
	
	if ( par::genome_2sets )
	  {	    
	    if ( ! ( (gset1.find(p1) != gset1.end() && 
		      gset2.find(p2) != gset2.end() ) ||
		     (gset1.find(p2) != gset1.end() && 
		      gset2.find(p1) != gset2.end() ) ) )
	      continue;
	  }

	
	// Are we set to only consider people in the same family
	// (i.e. confirm known relatedness)?

	if ( par::genome_only_check_rels )
	  {
	    if ( relType(p1,p2)=="UN" )
	      continue;	    
	  }



	//////////////////////////////////////////
	// Perform main calculations for this pair
	
	Z IBSg = calcGenomeIBS(p1,p2);
	Z IBDg = calcGenomeIBD(p1,p2,IBSg);

	
	if ( par::genome_test ) 
	  {
	    if ( IBDg.z1/2 + IBDg.z2 >= par::genome_test_threshold )
	      {
		int2 pair;
		pair.p1 = i1;
		pair.p2 = i2;
		related.insert(pair);
	      }
	  }
      }	
  
  if (!par::silent)
    cout << "\n";
  
  ZOUTFILE.close();
  
}




void Plink::calcGenomeIBM(Individual * p1, Individual * p2)
{

  // Individual-major mode assumed
  
  int cnt = 0;
  
  // Consider all loci

  vector<bool>::iterator ia1 = p1->one.begin();
  vector<bool>::iterator ia2 = p1->two.begin();
  
  vector<bool>::iterator ib1 = p2->one.begin();
  vector<bool>::iterator ib2 = p2->two.begin();
  
  while ( ia1 != p1->one.end() )
    {          

      // Count discordantly missing SNPs
      
      if ( *ia1 && !*ia2 )
	{
	  if ( ! ( *ib1 && !*ib2 ) ) cnt++;
	}
      else if ( *ib1 && !*ib2 ) cnt++;
	
      // Next SNP      
      ia1++;
      ia2++;
      ib1++;
      ib2++;

    } 

  // IBM similiarity metric
  dst = 1.0 - ( (double)cnt / (double)nl_all );
}



void Plink::pruneLD()
{

  if (!par::SNP_major)
    Ind2SNP();

  printLOG("Performing LD-based pruning...\n");
  
  string f = par::output_file_name + ".prune.in";
  ofstream PIN(f.c_str(),ios::out);
  printLOG("Writing pruned-in SNPs to [ " + f + " ]\n");

  f = par::output_file_name + ".prune.out";
  ofstream POUT(f.c_str(),ios::out);
  printLOG("Writing pruned-out SNPs to [ " + f + " ]\n");

  int win_start = 0;
  int win_end = win_start + par::prune_ld_win;
  
  ////////////////////////
  // Scan each chromosome
  
  vector<int> chrs;
  if (par::run_chr==0) 
    {
      vector<int> r = getChromosomeRange(*this);

      printLOG("Scanning from chromosome "+
	       chromosomeName( r[0] ) +" to "+
	       chromosomeName( r[1] ) +"\n\n");

      for (int i=r[0];i<=r[1];i++) 
	{
	  if (seeChromosome(*this,i))
	    chrs.push_back(i);
	}
    } 
  else 
    chrs.push_back(par::run_chr);
  
  // Inclusion or no?
  vector<bool> include(nl_all,true);

  // Only consider founders (set flag)
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( (*person)->founder ) 
	(*person)->flag = true;
      else
	(*person)->flag = false;
      person++;
    }
  
  // Scan each chromosome
  for (int i=0;i<chrs.size();i++)
    {

      // Skip chromosome 0 
      if ( chrs[i] == 0 ) 
	{
	  printLOG("Skippng chromosome 0\n");
	  continue;
	}

      // Set chromosome
      par::run_chr = chrs[i];
      
      // Find scan range
      setMarkerRange();
      
      // go from par::run_start to par::run_end
      int s1 = par::run_start;
      int s2 = par::run_start + par::prune_ld_win - 1;

      // MW fix
      if (s2 > par::run_end ) 
      	s2 = par::run_end;

      while ( s2 <= par::run_end )
	{
	  
	  // calc VIF and set 
	  vector<int> nSNP(0);
	  for (int l=s1; l<=s2; l++)
	    if ( include[l] )
	      {
		nSNP.push_back( l );
	      }
	  
	  // Skip if we only have a single SNP left
	  if (nSNP.size() < 2)
	    {
	      if ( s2 == par::run_end ) 
		break;
	      
	      s1 += par::prune_ld_step;
	      s2 += par::prune_ld_step;
	      
	      if (s2 > par::run_end) 
		s2 = par::run_end;
	      
	      if ( s2-s1 < 1)
		break;
	      
	      continue;
	    }

	  vector<vector<double> > variance;
	  
	  if (!par::silent)
	    {
	      cout << "Pruning SNPs " << s1-par::run_start+1
		   << " to " << s2-par::run_start+1
		   << " of " << par::run_end - par::run_start+1 
		   << "         \r";
	      cout.flush();
	    }

	  // Calculate covariance matrices
	  variance = calcSetCovarianceMatrix(nSNP);
	  
	  // Calculate VIFs
	  vector<bool> cur = vif_prune(variance,par::prune_ld_vif,nSNP);

	  // Update main list
	  int k=0;
	  for (int l=s1; l<=s2; l++)
	    {
	      // Update main list, but do not get back
	      // already excluded SNPs
	      
	      if (include[l] && !cur[k++])
		include[l] = false;
	    }
	  
	  // Advance window
	  if ( s2 == par::run_end ) 
	    break;
	  
	  s1 += par::prune_ld_step;
	  s2 += par::prune_ld_step;
	  
 	  if (s2 > par::run_end) 
 	    s2 = par::run_end;

	  if ( s2-s1 < 1)
	    break;
	  
	} // next window
      
      if (!par::silent)
	cout << "\n";

      // Record what is in, what is out
      int cnt_in = 0, cnt_out = 0;
      for (int l=par::run_start; l<=par::run_end; l++)
	{
	  if (include[l]) 
	    {
	      PIN << locus[l]->name << "\n";
	      cnt_in++;
	    }
	  else
	    {
	      POUT << locus[l]->name << "\n";
	      cnt_out++;
	    }
	}

      printLOG("For chromosome "+int2str(par::run_chr)+
	       ", "+int2str(cnt_out)+
	       " SNPs pruned out, "+int2str(cnt_in)+" remaining\n");
      
    } // next chromosome


  PIN.close();
  POUT.close();

}
